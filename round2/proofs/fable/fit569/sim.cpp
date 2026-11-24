// Fast exact simulator for Isotropic Redirection (IR), Direct Attachment to Nucleus (DAN),
// and Prohibited Attachment to Nucleus (PAN) growing trees from the 3-path L_3.
//
// Compile: g++ -O3 -march=native -std=c++20 -pthread sim.cpp -o sim
//
// Usage:
//   ./sim validate [num_trees_small] [num_trees_n100]
//   ./sim run --model <IR|DAN|PAN> --mode <coupled|indep> --replicates <R> --trees <T>
//             --threads <1|2> --seed <S> --out <csv_path>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max();

// SplitMix64 for seeding
static inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// xoshiro256++ PRNG (Blackman & Vigna)
struct Rng {
    uint64_t s[4];
    uint64_t buf = 0;
    bool has_buf = false;

    explicit Rng(uint64_t seed) {
        uint64_t x = seed;
        for (int i = 0; i < 4; ++i) s[i] = splitmix64(x);
        has_buf = false;
    }

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    inline uint64_t next64() {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    inline uint32_t next32() {
        if (has_buf) {
            has_buf = false;
            return static_cast<uint32_t>(buf);
        }
        uint64_t v = next64();
        buf = v >> 32;
        has_buf = true;
        return static_cast<uint32_t>(v);
    }

    // Lemire's exact unbiased uniform integer in [0, range)
    inline uint32_t uniform(uint32_t range) {
        uint32_t x = next32();
        uint64_t m = uint64_t(x) * uint64_t(range);
        uint32_t l = uint32_t(m);
        if (__builtin_expect(l < range, 0)) {
            uint32_t t = (-range) % range;
            while (l < t) {
                x = next32();
                m = uint64_t(x) * uint64_t(range);
                l = uint32_t(m);
            }
        }
        return static_cast<uint32_t>(m >> 32);
    }
};

enum class Model { IR, DAN, PAN };

static Model parse_model(const std::string& s) {
    if (s == "IR") return Model::IR;
    if (s == "DAN") return Model::DAN;
    if (s == "PAN") return Model::PAN;
    std::cerr << "Unknown model: " << s << "\n";
    std::exit(1);
}

static std::string model_name(Model m) {
    switch (m) {
        case Model::IR: return "IR";
        case Model::DAN: return "DAN";
        case Model::PAN: return "PAN";
    }
    return "";
}

// 17-point quarter-decade grid from 10^2 to 10^6 (contains 9-point half-decade and 5-point decade grids)
static constexpr int NUM_GRID = 17;
static constexpr uint32_t GRID_N[NUM_GRID] = {
    100, 178, 316, 562,
    1000, 1778, 3162, 5623,
    10000, 17783, 31623, 56234,
    100000, 177828, 316228, 562341,
    1000000
};

// Reusable workspace for fast tree simulation up to max_n vertices
struct SimWorkspace {
    std::vector<uint32_t> leaf_parent; // [0..N1-1] -> core index
    std::vector<uint32_t> next_leaf;   // [0..N1-1] -> linked list next leaf of same core node
    std::vector<uint32_t> head_leaf;   // [0..C-1] -> head leaf slot of core node
    std::vector<uint32_t> leaf_deg;    // [0..C-1] -> number of leaf neighbors of core node
    std::vector<std::vector<uint32_t>> core_nbrs; // [0..C-1] -> core neighbors (IR only)

    explicit SimWorkspace(uint32_t max_n, bool need_core_nbrs)
        : leaf_parent(max_n), next_leaf(max_n), head_leaf(max_n), leaf_deg(max_n) {
        if (need_core_nbrs) {
            core_nbrs.resize(std::min<uint32_t>(max_n, 65536));
        }
    }

    inline void ensure_core_capacity(uint32_t c) {
        if (c >= core_nbrs.size()) {
            core_nbrs.resize(std::max<size_t>(c + 1, core_nbrs.size() * 2));
        }
    }

    // Grow one tree from L_3 up to target_n, recording C at each grid point <= target_n.
    // If out_c != nullptr, writes C at grid indices 0..max_g_idx.
    // Returns final C at target_n.
    template <Model M>
    uint32_t simulate_to(uint32_t target_n, Rng& rng, uint32_t* out_c, int num_grid_to_record) {
        // Initial state: L_3 (1 core node 0, 2 leaves in slots 0 and 1)
        uint32_t C = 1;
        uint32_t N1 = 2;
        leaf_parent[0] = 0;
        leaf_parent[1] = 0;
        next_leaf[0] = NONE;
        next_leaf[1] = 0;
        head_leaf[0] = 1;
        leaf_deg[0] = 2;
        if constexpr (M == Model::IR) {
            core_nbrs[0].clear();
        }

        int g_idx = 0;
        uint32_t next_g = (num_grid_to_record > 0) ? GRID_N[0] : (target_n + 1);

        uint32_t* __restrict lp = leaf_parent.data();
        uint32_t* __restrict nl = next_leaf.data();
        uint32_t* __restrict hl = head_leaf.data();
        uint32_t* __restrict ld = leaf_deg.data();

        for (uint32_t n = 3; n < target_n; ) {
            uint32_t r = rng.uniform(n);
            if (__builtin_expect(r < N1, 1)) {
                // Selected a leaf -> attach new leaf to its unique (rank-1) neighbor W
                uint32_t W = lp[r];
                uint32_t idx = N1++;
                lp[idx] = W;
                nl[idx] = hl[W];
                hl[W] = idx;
                ld[W]++;
                ++n;
            } else {
                // Selected core node c
                uint32_t c = r - N1;
                if constexpr (M == Model::IR) {
                    uint32_t ldeg = ld[c];
                    uint32_t cdeg = static_cast<uint32_t>(core_nbrs[c].size());
                    uint32_t k = rng.uniform(ldeg + cdeg);
                    if (k < ldeg) {
                        // Redirect to a leaf neighbor of c -> that leaf becomes a new core node c_new
                        uint32_t u = hl[c];
                        hl[c] = nl[u];
                        ld[c] = ldeg - 1;

                        uint32_t c_new = C++;
                        ensure_core_capacity(c_new);
                        core_nbrs[c_new].clear();
                        core_nbrs[c].push_back(c_new);
                        core_nbrs[c_new].push_back(c);

                        // Reuse leaf slot u for the new leaf attached to c_new
                        lp[u] = c_new;
                        hl[c_new] = u;
                        nl[u] = NONE;
                        ld[c_new] = 1;
                    } else {
                        // Redirect to a core neighbor W of c -> attach new leaf to W
                        uint32_t W = core_nbrs[c][k - ldeg];
                        uint32_t idx = N1++;
                        lp[idx] = W;
                        nl[idx] = hl[W];
                        hl[W] = idx;
                        ld[W]++;
                    }
                    ++n;
                } else if constexpr (M == Model::DAN) {
                    uint32_t ldeg = ld[c];
                    if (ldeg > 0) {
                        // Rank-1 node -> redirect to a leaf neighbor -> new core node c_new
                        uint32_t u = hl[c];
                        hl[c] = nl[u];
                        ld[c] = ldeg - 1;

                        uint32_t c_new = C++;
                        lp[u] = c_new;
                        hl[c_new] = u;
                        nl[u] = NONE;
                        ld[c_new] = 1;
                    } else {
                        // Nucleus node -> attach directly to c
                        uint32_t idx = N1++;
                        lp[idx] = c;
                        nl[idx] = NONE;
                        hl[c] = idx;
                        ld[c] = 1;
                    }
                    ++n;
                } else if constexpr (M == Model::PAN) {
                    uint32_t ldeg = ld[c];
                    if (ldeg > 0) {
                        // Rank-1 node -> redirect to a leaf neighbor -> new core node c_new
                        uint32_t u = hl[c];
                        hl[c] = nl[u];
                        ld[c] = ldeg - 1;

                        uint32_t c_new = C++;
                        lp[u] = c_new;
                        hl[c_new] = u;
                        nl[u] = NONE;
                        ld[c_new] = 1;
                        ++n;
                    } else {
                        // Nucleus node -> reject event (do not increment n)
                        continue;
                    }
                }
            }

            if (__builtin_expect(n == next_g, 0)) {
                out_c[g_idx++] = C;
                next_g = (g_idx < num_grid_to_record) ? GRID_N[g_idx] : (target_n + 1);
            }
        }
        return C;
    }
};

// Naive adjacency-list reference simulator for independent verification up to small/medium N
static uint32_t simulate_naive(Model m, uint32_t target_n, Rng& rng, std::vector<uint32_t>* out_trace = nullptr) {
    std::vector<std::vector<uint32_t>> adj(target_n);
    adj[0].push_back(1);
    adj[1].push_back(0);
    adj[1].push_back(2);
    adj[2].push_back(1);

    uint32_t n = 3;
    if (out_trace) {
        out_trace->assign(target_n + 1, 0);
        (*out_trace)[3] = 1;
    }

    while (n < target_n) {
        uint32_t v = rng.uniform(n);
        uint32_t w = NONE;
        if (m == Model::IR) {
            w = adj[v][rng.uniform(static_cast<uint32_t>(adj[v].size()))];
        } else {
            if (adj[v].size() == 1) {
                w = adj[v][0];
            } else {
                // Collect leaf neighbors
                std::vector<uint32_t> leaf_nbrs;
                for (uint32_t nb : adj[v]) {
                    if (adj[nb].size() == 1) leaf_nbrs.push_back(nb);
                }
                if (!leaf_nbrs.empty()) {
                    w = leaf_nbrs[rng.uniform(static_cast<uint32_t>(leaf_nbrs.size()))];
                } else {
                    if (m == Model::DAN) {
                        w = v;
                    } else {
                        // PAN: reject
                        continue;
                    }
                }
            }
        }
        adj[w].push_back(n);
        adj[n].push_back(w);
        ++n;
        if (out_trace) {
            uint32_t c = 0;
            for (uint32_t i = 0; i < n; ++i) {
                if (adj[i].size() >= 2) ++c;
            }
            (*out_trace)[n] = c;
        }
    }
    uint32_t c = 0;
    for (uint32_t i = 0; i < target_n; ++i) {
        if (adj[i].size() >= 2) ++c;
    }
    return c;
}

static void run_validation(uint64_t num_small, uint64_t num_n100) {
    // Exact E[C_N] for N=4..8 from validate_exact.py
    const double exact_IR[5]  = {4.0/3.0, 19.0/12.0, 647.0/360.0, 17189.0/8640.0, 4717423.0/2177280.0};
    const double exact_DAN[5] = {4.0/3.0, 5.0/3.0,   59.0/30.0,   9.0/4.0,        794.0/315.0};
    const double exact_PAN[5] = {4.0/3.0, 5.0/3.0,   119.0/60.0,  457.0/200.0,    972751.0/378000.0};

    auto check_model = [&](Model m, const double* exact) {
        std::cout << "=== Validating fast simulator: " << model_name(m) << " ===\n";
        SimWorkspace ws(128, m == Model::IR);
        Rng rng(123456789ULL + static_cast<uint64_t>(m) * 99991ULL);

        // 1. Compare against exact values for N = 4..8
        double sum_c[5] = {0}, sum_c2[5] = {0};
        uint64_t cnt_star[5] = {0}, cnt_line[5] = {0};
        for (uint64_t t = 0; t < num_small; ++t) {
            // Simulate step-by-step to N=8 using exact same logic
            // Just run to each N=8 by recording at N=4,5,6,7,8
            uint32_t C = 1, N1 = 2;
            ws.leaf_parent[0] = 0; ws.leaf_parent[1] = 0;
            ws.next_leaf[0] = NONE; ws.next_leaf[1] = 0;
            ws.head_leaf[0] = 1; ws.leaf_deg[0] = 2;
            if (m == Model::IR) ws.core_nbrs[0].clear();

            for (uint32_t n = 3; n < 8; ) {
                uint32_t r = rng.uniform(n);
                if (r < N1) {
                    uint32_t W = ws.leaf_parent[r];
                    uint32_t idx = N1++;
                    ws.leaf_parent[idx] = W;
                    ws.next_leaf[idx] = ws.head_leaf[W];
                    ws.head_leaf[W] = idx;
                    ws.leaf_deg[W]++;
                    ++n;
                } else {
                    uint32_t c = r - N1;
                    if (m == Model::IR) {
                        uint32_t ldeg = ws.leaf_deg[c];
                        uint32_t cdeg = static_cast<uint32_t>(ws.core_nbrs[c].size());
                        uint32_t k = rng.uniform(ldeg + cdeg);
                        if (k < ldeg) {
                            uint32_t u = ws.head_leaf[c];
                            ws.head_leaf[c] = ws.next_leaf[u];
                            ws.leaf_deg[c] = ldeg - 1;
                            uint32_t c_new = C++;
                            ws.core_nbrs[c_new].clear();
                            ws.core_nbrs[c].push_back(c_new);
                            ws.core_nbrs[c_new].push_back(c);
                            ws.leaf_parent[u] = c_new;
                            ws.head_leaf[c_new] = u;
                            ws.next_leaf[u] = NONE;
                            ws.leaf_deg[c_new] = 1;
                        } else {
                            uint32_t W = ws.core_nbrs[c][k - ldeg];
                            uint32_t idx = N1++;
                            ws.leaf_parent[idx] = W;
                            ws.next_leaf[idx] = ws.head_leaf[W];
                            ws.head_leaf[W] = idx;
                            ws.leaf_deg[W]++;
                        }
                        ++n;
                    } else if (m == Model::DAN) {
                        uint32_t ldeg = ws.leaf_deg[c];
                        if (ldeg > 0) {
                            uint32_t u = ws.head_leaf[c];
                            ws.head_leaf[c] = ws.next_leaf[u];
                            ws.leaf_deg[c] = ldeg - 1;
                            uint32_t c_new = C++;
                            ws.leaf_parent[u] = c_new;
                            ws.head_leaf[c_new] = u;
                            ws.next_leaf[u] = NONE;
                            ws.leaf_deg[c_new] = 1;
                        } else {
                            uint32_t idx = N1++;
                            ws.leaf_parent[idx] = c;
                            ws.next_leaf[idx] = NONE;
                            ws.head_leaf[c] = idx;
                            ws.leaf_deg[c] = 1;
                        }
                        ++n;
                    } else {
                        uint32_t ldeg = ws.leaf_deg[c];
                        if (ldeg > 0) {
                            uint32_t u = ws.head_leaf[c];
                            ws.head_leaf[c] = ws.next_leaf[u];
                            ws.leaf_deg[c] = ldeg - 1;
                            uint32_t c_new = C++;
                            ws.leaf_parent[u] = c_new;
                            ws.head_leaf[c_new] = u;
                            ws.next_leaf[u] = NONE;
                            ws.leaf_deg[c_new] = 1;
                            ++n;
                        } else {
                            continue;
                        }
                    }
                }
                int idx_n = static_cast<int>(n) - 4;
                sum_c[idx_n] += C;
                sum_c2[idx_n] += double(C) * double(C);
                if (C == 1) cnt_star[idx_n]++;
                if (C == n - 2) cnt_line[idx_n]++;
            }
        }
        for (int i = 0; i < 5; ++i) {
            uint32_t N = 4 + i;
            double mean = sum_c[i] / num_small;
            double var = sum_c2[i] / num_small - mean * mean;
            double se = std::sqrt(var / num_small);
            double z = (mean - exact[i]) / se;
            std::cout << std::fixed << std::setprecision(8)
                      << "  N=" << N << ": sim=" << mean << " exact=" << exact[i]
                      << " se=" << se << " z=" << std::setprecision(2) << z
                      << "  S_N=" << std::setprecision(6) << double(cnt_star[i])/num_small
                      << " (exact " << 2.0/(N-1) << ")"
                      << "  Q_N=" << double(cnt_line[i])/num_small << "\n";
        }

        // 2. Compare fast simulate_to() vs naive adjacency-list simulator at N=100
        double fast_s = 0, fast_s2 = 0, naive_s = 0, naive_s2 = 0;
        Rng rng_fast(987654321ULL + static_cast<uint64_t>(m) * 777ULL);
        Rng rng_naive(555555555ULL + static_cast<uint64_t>(m) * 777ULL);
        for (uint64_t t = 0; t < num_n100; ++t) {
            uint32_t cf = 0;
            if (m == Model::IR) cf = ws.simulate_to<Model::IR>(100, rng_fast, nullptr, 0);
            else if (m == Model::DAN) cf = ws.simulate_to<Model::DAN>(100, rng_fast, nullptr, 0);
            else cf = ws.simulate_to<Model::PAN>(100, rng_fast, nullptr, 0);
            fast_s += cf; fast_s2 += double(cf) * double(cf);

            uint32_t cn = simulate_naive(m, 100, rng_naive);
            naive_s += cn; naive_s2 += double(cn) * double(cn);
        }
        double mf = fast_s / num_n100, vf = fast_s2 / num_n100 - mf * mf;
        double mn = naive_s / num_n100, vn = naive_s2 / num_n100 - mn * mn;
        double se_diff = std::sqrt((vf + vn) / num_n100);
        double z_diff = (mf - mn) / se_diff;
        std::cout << std::fixed << std::setprecision(6)
                  << "  N=100 cross-check (T=" << num_n100 << "): "
                  << "fast=" << mf << " (sd=" << std::sqrt(vf) << ") vs "
                  << "naive=" << mn << " (sd=" << std::sqrt(vn) << "), z="
                  << std::setprecision(2) << z_diff << "\n\n";
    };

    check_model(Model::IR, exact_IR);
    check_model(Model::DAN, exact_DAN);
    check_model(Model::PAN, exact_PAN);
}

struct ReplicateResult {
    int rep_id = 0;
    double sum_c[NUM_GRID] = {0};
    double sum_c2[NUM_GRID] = {0};
    double sum_cross[NUM_GRID * NUM_GRID] = {0}; // only used in coupled mode
};

template <Model M>
static void run_replicate_coupled(int rep_id, int trees_per_rep, uint64_t base_seed,
                                  SimWorkspace& ws, ReplicateResult& res) {
    res.rep_id = rep_id;
    std::memset(res.sum_c, 0, sizeof(res.sum_c));
    std::memset(res.sum_c2, 0, sizeof(res.sum_c2));
    std::memset(res.sum_cross, 0, sizeof(res.sum_cross));

    // Exact injective seed mixer
    auto mix_seed = [](uint64_t base, uint64_t a, uint64_t b, uint64_t c) -> uint64_t {
        uint64_t x = base;
        x = splitmix64(x) ^ a;
        x = splitmix64(x) ^ b;
        x = splitmix64(x) ^ c;
        return splitmix64(x);
    };

    uint32_t c_vals[NUM_GRID];
    for (int t = 0; t < trees_per_rep; ++t) {
        uint64_t seed = mix_seed(base_seed, uint64_t(rep_id), 0ULL, uint64_t(t));
        Rng rng(seed);
        ws.simulate_to<M>(GRID_N[NUM_GRID - 1], rng, c_vals, NUM_GRID);
        for (int i = 0; i < NUM_GRID; ++i) {
            double ci = static_cast<double>(c_vals[i]);
            res.sum_c[i] += ci;
            res.sum_c2[i] += ci * ci;
            for (int j = i; j < NUM_GRID; ++j) {
                res.sum_cross[i * NUM_GRID + j] += ci * static_cast<double>(c_vals[j]);
            }
        }
    }
    for (int i = 0; i < NUM_GRID; ++i) {
        for (int j = 0; j < i; ++j) {
            res.sum_cross[i * NUM_GRID + j] = res.sum_cross[j * NUM_GRID + i];
        }
    }
}

template <Model M>
static void run_replicate_indep(int rep_id, int trees_per_rep, uint64_t base_seed,
                                SimWorkspace& ws, ReplicateResult& res) {
    res.rep_id = rep_id;
    std::memset(res.sum_c, 0, sizeof(res.sum_c));
    std::memset(res.sum_c2, 0, sizeof(res.sum_c2));
    std::memset(res.sum_cross, 0, sizeof(res.sum_cross));

    auto mix_seed = [](uint64_t base, uint64_t a, uint64_t b, uint64_t c) -> uint64_t {
        uint64_t x = base;
        x = splitmix64(x) ^ a;
        x = splitmix64(x) ^ b;
        x = splitmix64(x) ^ c;
        return splitmix64(x);
    };

    for (int g = 0; g < NUM_GRID; ++g) {
        uint32_t target_n = GRID_N[g];
        double sc = 0.0, sc2 = 0.0;
        for (int t = 0; t < trees_per_rep; ++t) {
            uint64_t seed = mix_seed(base_seed, uint64_t(rep_id), uint64_t(g + 1), uint64_t(t));
            Rng rng(seed);
            uint32_t c = ws.simulate_to<M>(target_n, rng, nullptr, 0);
            double cd = static_cast<double>(c);
            sc += cd;
            sc2 += cd * cd;
        }
        res.sum_c[g] = sc;
        res.sum_c2[g] = sc2;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: sim validate [num_small] [num_n100] | sim run ...\n";
        return 1;
    }
    std::string cmd = argv[1];
    if (cmd == "validate") {
        uint64_t num_small = (argc >= 3) ? std::stoull(argv[2]) : 2000000ULL;
        uint64_t num_n100  = (argc >= 4) ? std::stoull(argv[3]) : 200000ULL;
        run_validation(num_small, num_n100);
        return 0;
    }
    if (cmd != "run") {
        std::cerr << "Unknown command: " << cmd << "\n";
        return 1;
    }

    std::string model_str = "IR";
    std::string mode_str = "coupled";
    int num_reps = 40;
    int trees_per_rep = 2400;
    int num_threads = 2;
    uint64_t base_seed = 42;
    std::string out_prefix = "out";

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        auto next_val = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << "Missing value for " << arg << "\n"; std::exit(1); }
            return argv[++i];
        };
        if (arg == "--model") model_str = next_val();
        else if (arg == "--mode") mode_str = next_val();
        else if (arg == "--replicates") num_reps = std::stoi(next_val());
        else if (arg == "--trees") trees_per_rep = std::stoi(next_val());
        else if (arg == "--threads") num_threads = std::stoi(next_val());
        else if (arg == "--seed") base_seed = std::stoull(next_val());
        else if (arg == "--out") out_prefix = next_val();
        else { std::cerr << "Unknown arg: " << arg << "\n"; return 1; }
    }

    if (num_threads < 1) num_threads = 1;
    if (num_threads > 2) {
        std::cerr << "Capping threads to 2 per CPU budget.\n";
        num_threads = 2;
    }
    Model model = parse_model(model_str);
    bool is_coupled = (mode_str == "coupled");
    if (!is_coupled && mode_str != "indep") {
        std::cerr << "Invalid --mode (must be coupled or indep): " << mode_str << "\n";
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();
    std::vector<ReplicateResult> results(num_reps);
    std::atomic<int> next_rep{0};
    std::atomic<int> done_reps{0};

    auto worker = [&](int /*tid*/) {
        SimWorkspace ws(GRID_N[NUM_GRID - 1] + 16, model == Model::IR);
        while (true) {
            int r = next_rep.fetch_add(1);
            if (r >= num_reps) break;
            if (is_coupled) {
                if (model == Model::IR) run_replicate_coupled<Model::IR>(r, trees_per_rep, base_seed, ws, results[r]);
                else if (model == Model::DAN) run_replicate_coupled<Model::DAN>(r, trees_per_rep, base_seed, ws, results[r]);
                else run_replicate_coupled<Model::PAN>(r, trees_per_rep, base_seed, ws, results[r]);
            } else {
                if (model == Model::IR) run_replicate_indep<Model::IR>(r, trees_per_rep, base_seed, ws, results[r]);
                else if (model == Model::DAN) run_replicate_indep<Model::DAN>(r, trees_per_rep, base_seed, ws, results[r]);
                else run_replicate_indep<Model::PAN>(r, trees_per_rep, base_seed, ws, results[r]);
            }
            int d = done_reps.fetch_add(1) + 1;
            if (d % 10 == 0 || d == num_reps) {
                auto now = std::chrono::steady_clock::now();
                double sec = std::chrono::duration<double>(now - t0).count();
                std::cerr << "[" << model_str << "/" << mode_str << "] " << d << "/" << num_reps
                          << " replicates done (" << std::fixed << std::setprecision(1) << sec << " s)\n";
            }
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) threads.emplace_back(worker, i);
    for (auto& th : threads) th.join();

    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    // Write per-replicate means & sums CSV
    std::string rep_csv = out_prefix + "_reps.csv";
    {
        std::ofstream ofs(rep_csv);
        ofs << std::setprecision(15);
        ofs << "rep,trees";
        for (int g = 0; g < NUM_GRID; ++g) ofs << ",mean_" << GRID_N[g];
        for (int g = 0; g < NUM_GRID; ++g) ofs << ",sum2_" << GRID_N[g];
        ofs << "\n";
        for (int r = 0; r < num_reps; ++r) {
            ofs << r << "," << trees_per_rep;
            for (int g = 0; g < NUM_GRID; ++g) ofs << "," << (results[r].sum_c[g] / trees_per_rep);
            for (int g = 0; g < NUM_GRID; ++g) ofs << "," << results[r].sum_c2[g];
            ofs << "\n";
        }
    }

    // Write pooled covariance / summary CSV (exact pooled sums across all R * T trees)
    std::string cov_csv = out_prefix + "_pooled_cov.csv";
    {
        double total_trees = double(num_reps) * double(trees_per_rep);
        double tot_c[NUM_GRID] = {0}, tot_c2[NUM_GRID] = {0};
        double tot_cross[NUM_GRID * NUM_GRID] = {0};
        for (int r = 0; r < num_reps; ++r) {
            for (int i = 0; i < NUM_GRID; ++i) {
                tot_c[i] += results[r].sum_c[i];
                tot_c2[i] += results[r].sum_c2[i];
                for (int j = 0; j < NUM_GRID; ++j) {
                    tot_cross[i * NUM_GRID + j] += results[r].sum_cross[i * NUM_GRID + j];
                }
            }
        }
        std::ofstream ofs(cov_csv);
        ofs << std::setprecision(15);
        ofs << "i,N,total_trees,mean,sd,cv";
        for (int j = 0; j < NUM_GRID; ++j) ofs << ",cov_" << GRID_N[j];
        ofs << "\n";
        for (int i = 0; i < NUM_GRID; ++i) {
            double mi = tot_c[i] / total_trees;
            double vi = (tot_c2[i] - total_trees * mi * mi) / (total_trees - 1.0);
            double sdi = std::sqrt(std::max(0.0, vi));
            ofs << i << "," << GRID_N[i] << "," << static_cast<uint64_t>(total_trees)
                << "," << mi << "," << sdi << "," << (sdi / mi);
            for (int j = 0; j < NUM_GRID; ++j) {
                double cov_ij = 0.0;
                if (i == j) {
                    cov_ij = vi;
                } else if (is_coupled) {
                    double mj = tot_c[j] / total_trees;
                    cov_ij = (tot_cross[i * NUM_GRID + j] - total_trees * mi * mj) / (total_trees - 1.0);
                }
                ofs << "," << cov_ij;
            }
            ofs << "\n";
        }
    }

    std::cout << "Done " << model_str << " (" << mode_str << "): " << num_reps << " reps x "
              << trees_per_rep << " trees in " << std::fixed << std::setprecision(2) << elapsed << " s -> "
              << rep_csv << "\n";
    return 0;
}
