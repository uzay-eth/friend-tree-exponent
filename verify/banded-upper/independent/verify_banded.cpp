// Independent exact upper-bound verifier for the multi-band random friend tree
// certificate specified in attempts/strogatz-astra/banded/BOUND.md.
//
// Build (single CPU, low priority, UBSan):
//   g++ -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all verify_banded.cpp -o verify_banded
//
// Usage:
//   ./verify_banded <lambda_num> <lambda_den> <bands_file> <recipe_file>
//
// Architecture:
//   - Processes one central-degree band at a time from D = H-1 down to D = 2
//     (sliding window over D and D+1), keeping peak RAM low even at >200M states.
//   - Bags are stored in run-length encoded (RLE) CSR format in lexicographic order;
//     lookups use exact binary search directly on the RLE array (zero std::map).
//   - Every row finds the unique minimal positive integer z >= 1 satisfying the
//     exact 128-bit upward-rounded drift inequality against trial newborn rewards,
//     then independently asserts z >= 1 and exact 128-bit slack >= 0.
//   - Finally verifies actual newborn weights <= trial rewards for all c=1 labels
//     (which guarantees every finite inequality also holds with actual newborn
//     weights by monotonicity of L/D >= 0) and verifies the analytic D >= H tail
//     inequalities in exact 128-bit arithmetic against the actual pooled-H newborn.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

using i64 = int64_t;
using i128 = __int128_t;

static constexpr int WITNESS_CAP = 6;   // exact witness e in 2..5; pooled e >= 6
static constexpr int BRIDGE_DEG = 14;   // max sole-neighbour degree for witness
static constexpr int WITNESS_CORE = 4;  // witness active for k_w in 2..4

static int core_cap_for(int c) {
    if (c == 1) return 20;
    if (c == 2) return 7;
    return 5;
}

static i64 to_i64_checked(i128 x, const char* ctx) {
    if (x > numeric_limits<i64>::max() || x < numeric_limits<i64>::min()) {
        throw runtime_error(string("64-bit integer overflow in ") + ctx);
    }
    return static_cast<i64>(x);
}

static i128 ceil_div128(i128 a, i64 b) {
    if (b <= 0) throw runtime_error("nonpositive denominator in ceil_div128");
    i128 q = a / b;
    i128 r = a % b;
    if (r > 0) ++q;
    return q;
}

static i128 floor_div128(i128 a, i64 b) {
    if (b <= 0) throw runtime_error("nonpositive denominator in floor_div128");
    i128 q = a / b;
    i128 r = a % b;
    if (r < 0) --q;
    return q;
}

static string strip_comment(const string& line) {
    auto pos = line.find('#');
    return (pos == string::npos) ? line : line.substr(0, pos);
}

// -----------------------------------------------------------------------------
// Band geometry
// -----------------------------------------------------------------------------
struct BandSpec {
    int d_start = 0; // inclusive
    int d_end = 0;   // exclusive
    vector<int> h;   // h[c] for c >= 1 (h[0] = 2); missing entries default to 2

    int cap(int c) const {
        if (c <= 0) return 2;
        return (c < (int)h.size()) ? h[c] : 2;
    }
};

struct Geometry {
    vector<BandSpec> bands;
    int H = 0;
    int h1 = 0;
    vector<int> band_of_D;

    static Geometry load(const string& path) {
        ifstream in(path);
        if (!in) throw runtime_error("cannot open bands file: " + path);
        Geometry G;
        string raw;
        int line_no = 0;
        while (getline(in, raw)) {
            ++line_no;
            string s = strip_comment(raw);
            istringstream iss(s);
            int d_end = 0;
            if (!(iss >> d_end)) continue;
            BandSpec b;
            b.d_start = G.bands.empty() ? 2 : G.bands.back().d_end;
            b.d_end = d_end;
            if (b.d_end <= b.d_start) {
                throw runtime_error("band end must strictly exceed start at line " + to_string(line_no));
            }
            b.h.push_back(2); // c = 0
            int val = 0;
            while (iss >> val) {
                if (val < 2) throw runtime_error("cutoff < 2 at line " + to_string(line_no));
                b.h.push_back(val);
            }
            if (b.h.size() < 2) b.h.push_back(2);
            G.bands.push_back(move(b));
        }
        if (G.bands.empty()) throw runtime_error("empty bands file");
        G.H = G.bands.back().d_end;
        if (G.H < 4) throw runtime_error("H must be at least 4");
        G.h1 = G.bands[0].cap(1);

        int max_c_listed = 1;
        for (const auto& b : G.bands) max_c_listed = max(max_c_listed, (int)b.h.size());
        for (size_t bi = 0; bi < G.bands.size(); ++bi) {
            const auto& b = G.bands[bi];
            if (b.cap(1) != G.h1) {
                throw runtime_error("h_1 must be constant across all bands");
            }
            for (int c = 1; c <= max_c_listed; ++c) {
                int hc = b.cap(c);
                if (hc < 2 || hc > G.H) throw runtime_error("cutoff outside [2, H]");
                if (c >= 2 && hc > b.cap(c - 1)) throw runtime_error("cutoffs must be nonincreasing in c");
                if (bi > 0 && hc > G.bands[bi - 1].cap(c))
                    throw runtime_error("cutoffs must be nonincreasing across bands");
            }
        }

        G.band_of_D.assign(G.H, -1);
        for (size_t bi = 0; bi < G.bands.size(); ++bi) {
            for (int D = G.bands[bi].d_start; D < G.bands[bi].d_end; ++D) {
                G.band_of_D[D] = (int)bi;
            }
        }

        // Witness initialization constraint (BOUND.md):
        // Creator (D, c) with L = D - c >= 1 creating a witnessed newborn (D <= 14, D < h1, c in 1..3)
        // must have h_c(band(D)) >= 6 so degrees < 6 are exact and pooled labels certify e >= 6.
        int w_d_max = min({BRIDGE_DEG, G.h1 - 1, G.H - 1});
        for (int D = 2; D <= w_d_max; ++D) {
            int bi = G.band_of_D[D];
            for (int c = 1; c <= WITNESS_CORE - 1; ++c) {
                if (D - c >= 1 && G.bands[bi].cap(c) < WITNESS_CAP) {
                    throw runtime_error("witness cutoff violation: h_c(D) < 6 for creator D=" +
                                        to_string(D) + ", c=" + to_string(c));
                }
            }
        }
        return G;
    }

    int cap(int D, int c) const { return bands[band_of_D[D]].cap(c); }

    bool edge_observable(int D, int c, int d, int k, bool* used_boundary_cap = nullptr) const {
        if (D < 2 || D >= H || d < 2 || d >= H) return false;
        if (c < 1 || c > D || k < 1 || k > d) return false;
        int hc = cap(D, c);
        int hk = cap(d, k);
        if (d >= hc || D >= hk) return false;
        int Kc = core_cap_for(c);
        int Kk = core_cap_for(k);
        bool k_ok = (k < Kc) || (k == Kc && d == Kc);
        bool c_ok = (c < Kk) || (c == Kk && D == Kk);
        if (k_ok && c_ok && used_boundary_cap && ((k == Kc && d == Kc) || (c == Kk && D == Kk))) {
            *used_boundary_cap = true;
        }
        return k_ok && c_ok;
    }
};

// -----------------------------------------------------------------------------
// Layer (band, c): labels, RLE bags (sorted lexicographically), and CSR moves
// -----------------------------------------------------------------------------
struct Label {
    int d = 0;   // exact degree (< pool) or lower-bound mark (>= pool)
    int k = 0;   // capped core degree (< pool) or 0 (>= pool)
    int wd = 0;  // witness degree 2..5 (exact), 6 (pooled witness), 0 (none)
    int wk = 0;  // witness core degree 1..wd (exact), 0 (pooled or none)
};

struct Run {
    int type = 0;
    int count = 0;
};

struct Move {
    int target = -1;
    int count = 0;
    int type = -1;
    uint8_t flags = 0; // bit 0: grow, bit 1: witness
    bool is_grow() const { return (flags & 1) != 0; }
    bool is_witness() const { return (flags & 2) != 0; }
};

// Compare two RLE-encoded sorted multisets of the same total cardinality in
// lexicographic order of their expanded sorted arrays.
static int cmp_rle(const Run* a, int na, const Run* b, int nb) {
    int i = 0;
    int lim = min(na, nb);
    while (i < lim) {
        if (a[i].type != b[i].type) return (a[i].type < b[i].type) ? -1 : 1;
        if (a[i].count != b[i].count) {
            // Same type, different multiplicity: the one with larger multiplicity
            // stays at this smaller type while the other advances to a larger type.
            return (a[i].count > b[i].count) ? -1 : 1;
        }
        ++i;
    }
    return 0;
}

static void canonicalize_runs(vector<Run>& v, int pool) {
    if (v.empty()) return;
    sort(v.begin(), v.end(), [](const Run& x, const Run& y) { return x.type < y.type; });
    // Merge identical types
    size_t w = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        if (w > 0 && v[w - 1].type == v[i].type) {
            v[w - 1].count += v[i].count;
        } else {
            v[w++] = v[i];
        }
    }
    v.resize(w);

    // Keep at most ONE copy of the largest pooled mark > pool; reset all other copies > pool to pool.
    int reset_cnt = 0;
    int last_above = -1;
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i].type > pool) last_above = i;
    }
    if (last_above >= 0) {
        for (int i = 0; i < last_above; ++i) {
            if (v[i].type > pool) {
                reset_cnt += v[i].count;
                v[i].count = 0;
            }
        }
        if (v[last_above].count > 1) {
            reset_cnt += v[last_above].count - 1;
            v[last_above].count = 1;
        }
    }
    if (reset_cnt > 0) {
        // Remove zeroed runs and add reset_cnt to pool
        vector<Run> out;
        out.reserve(v.size() + 1);
        bool added_pool = false;
        for (const auto& r : v) {
            if (r.count == 0) continue;
            if (!added_pool && r.type > pool) {
                out.push_back({pool, reset_cnt});
                added_pool = true;
            }
            if (r.type == pool) {
                out.push_back({pool, r.count + reset_cnt});
                added_pool = true;
            } else {
                out.push_back(r);
            }
        }
        if (!added_pool) out.push_back({pool, reset_cnt});
        v.swap(out);
    }
}

struct Layer {
    int band = 0;
    int c = 0;
    int cap = 2;
    int k_limit = 5;
    int pool = 0;
    int H = 0;

    vector<Label> types;
    vector<vector<int>> ids;

    int n_rows = 0;
    vector<int> run_off; // length n_rows + 1
    vector<Run> runs;    // sorted lexicographically across rows!

    vector<int> convert_tgt; // row index in (band, c + 1)
    vector<int> band_up_tgt; // row index in (band + 1, c)
    vector<int> move_off;    // length n_rows + 1
    vector<Move> moves;

    vector<i64> arrow_max;   // A_max[t] for sender arrow charge

    void init_labels(int band_, int c_, int cap_, int H_) {
        band = band_;
        c = c_;
        cap = cap_;
        H = H_;
        k_limit = core_cap_for(c);

        types.clear();
        ids.assign(cap, {});
        for (int d = 2; d < cap; ++d) {
            int kmax = min(d, k_limit);
            ids[d].assign(kmax + 1, -1);
            for (int k = 1; k <= kmax; ++k) {
                ids[d][k] = (int)types.size();
                if (c == 1 && d <= BRIDGE_DEG && k >= 2 && k <= WITNESS_CORE) {
                    for (int e = 2; e < WITNESS_CAP; ++e) {
                        for (int j = 1; j <= e; ++j) types.push_back({d, k, e, j});
                    }
                    types.push_back({d, k, WITNESS_CAP, 0});
                } else {
                    types.push_back({d, k, 0, 0});
                }
            }
        }
        pool = (int)types.size();
        types.push_back({cap, 0, 0, 0});
        for (int m = 4; m < H; m *= 2) {
            if (m > cap) types.push_back({m, 0, 0, 0});
        }
        if (H > cap) types.push_back({H, 0, 0, 0});
    }

    void enumerate_bags() {
        n_rows = 0;
        run_off.clear();
        runs.clear();
        run_off.push_back(0);
        if (c == 0) {
            n_rows = 1;
            run_off.push_back(0);
            return;
        }
        if (pool == 0) {
            // Fast path when h_c == 2: all neighbours are pooled!
            // Bags in lexicographic order: [0,...,0] (mark h_c), then [0,...,0, m] for m = 1..types.size()-1.
            int nt = (int)types.size();
            for (int m = 0; m < nt; ++m) {
                if (m == 0) {
                    runs.push_back({0, c});
                } else {
                    if (c > 1) runs.push_back({0, c - 1});
                    runs.push_back({m, 1});
                }
                ++n_rows;
                run_off.push_back((int)runs.size());
            }
            return;
        }
        vector<Run> cur;
        cur.reserve(c);
        enumerate_rec(cur, 0, 0);
    }

    // Generates RLE bags in strictly increasing lexicographic order of their
    // expanded sorted multisets: at each position choose next type t >= start;
    // increasing t first (before decreasing multiplicity) maintains exact lex order!
    // Check: with cur_sz elements chosen, choosing t = start extends the previous run
    // (if cur non-empty and cur.back().type == t), which is lexicographically SMALLEST!
    void enumerate_rec(vector<Run>& cur, int start, int cur_sz) {
        if (cur_sz == c) {
            ++n_rows;
            for (const auto& r : cur) runs.push_back(r);
            run_off.push_back((int)runs.size());
            return;
        }
        int nt = (int)types.size();
        for (int t = start; t < nt; ++t) {
            if (t > pool && cur_sz + 1 < c) break;
            if (!cur.empty() && cur.back().type == t) {
                ++cur.back().count;
                enumerate_rec(cur, t, cur_sz + 1);
                --cur.back().count;
            } else {
                cur.push_back({t, 1});
                enumerate_rec(cur, t, cur_sz + 1);
                cur.pop_back();
            }
        }
    }

    int label_plain(int d, int k) const {
        if (d >= cap) {
            int t = pool;
            while (t + 1 < (int)types.size() && types[t + 1].d <= d) ++t;
            return t;
        }
        return ids[d].at(min(k, k_limit));
    }

    int label_bridge(int host_d, int host_k, int wd, int wk) const {
        int base = label_plain(host_d, host_k);
        if (c != 1 || host_d > BRIDGE_DEG || host_d >= cap ||
            host_k < 2 || host_k > WITNESS_CORE) {
            return base;
        }
        static constexpr int witness_extra = (WITNESS_CAP - 2) * (WITNESS_CAP + 1) / 2;
        if (wd >= WITNESS_CAP) return base + witness_extra;
        return base + (wd - 2) * (wd + 1) / 2 + (wk - 1);
    }

    int lookup_runs(vector<Run>& key) const {
        canonicalize_runs(key, pool);
        int lo = 0, hi = n_rows - 1;
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            const Run* rptr = runs.data() + run_off[mid];
            int rlen = run_off[mid + 1] - run_off[mid];
            int cmp = cmp_rle(rptr, rlen, key.data(), (int)key.size());
            if (cmp == 0) return mid;
            if (cmp < 0) lo = mid + 1;
            else hi = mid - 1;
        }
        throw runtime_error("internal error: canonical RLE bag not found in layer");
    }

    bool is_k_exact(int t) const {
        if (t >= pool) return false;
        const auto& lb = types[t];
        return (lb.k < k_limit) || (lb.k == k_limit && lb.d == k_limit);
    }
};

// -----------------------------------------------------------------------------
// Recipe data
// -----------------------------------------------------------------------------
using EdgeKey = array<int, 4>; // {D_rec, c_rec, d_send, k_send}
using ArrowKey = pair<int, int>; // {D, c}

struct Recipe {
    i64 scale = 0;
    int threshold = -1;
    vector<i64> rewards;
    map<EdgeKey, i64> edges;
    map<ArrowKey, i64> arrows;
    i64 a_star = 0;
    bool used_boundary_cap_edge = false;

    static Recipe load(const string& path, const Geometry& G, int n_c1_labels) {
        ifstream in(path);
        if (!in) throw runtime_error("cannot open recipe file: " + path);

        Recipe R;
        bool has_scale = false;
        bool has_threshold = false;
        vector<bool> seen_reward(n_c1_labels, false);
        R.rewards.assign(n_c1_labels, 0);

        string raw;
        int line_no = 0;
        while (getline(in, raw)) {
            ++line_no;
            string s = strip_comment(raw);
            istringstream iss(s);
            string kind;
            if (!(iss >> kind)) continue;

            if (kind == "scale") {
                i64 n = 0;
                if (!(iss >> n) || n <= 0)
                    throw runtime_error("invalid scale at line " + to_string(line_no));
                if (has_scale)
                    throw runtime_error("duplicate scale at line " + to_string(line_no));
                R.scale = n;
                has_scale = true;
            } else if (kind == "threshold") {
                int T = 0;
                if (!(iss >> T) || T < 2)
                    throw runtime_error("invalid threshold at line " + to_string(line_no));
                if (has_threshold)
                    throw runtime_error("duplicate threshold at line " + to_string(line_no));
                R.threshold = T;
                has_threshold = true;
            } else if (kind == "reward") {
                int idx = -1;
                i64 w = 0;
                if (!(iss >> idx >> w))
                    throw runtime_error("malformed reward at line " + to_string(line_no));
                if (idx < 0 || idx >= n_c1_labels)
                    throw runtime_error("reward index out of range at line " + to_string(line_no));
                if (seen_reward[idx])
                    throw runtime_error("duplicate reward index at line " + to_string(line_no));
                if (w < 1)
                    throw runtime_error("nonpositive reward weight at line " + to_string(line_no));
                seen_reward[idx] = true;
                R.rewards[idx] = w;
            } else if (kind == "edge") {
                int D = 0, c = 0, d = 0, k = 0;
                i64 v = 0;
                if (!(iss >> D >> c >> d >> k >> v))
                    throw runtime_error("malformed edge at line " + to_string(line_no));
                EdgeKey key = {D, c, d, k};
                if (R.edges.count(key))
                    throw runtime_error("duplicate edge key at line " + to_string(line_no));
                if (v != 0) {
                    if (!G.edge_observable(D, c, d, k, &R.used_boundary_cap_edge)) {
                        throw runtime_error("unobservable nonzero edge transfer at line " + to_string(line_no));
                    }
                    R.edges[key] = v;
                }
            } else if (kind == "arrow") {
                int D = 0, c = 0;
                i64 v = 0;
                if (!(iss >> D >> c >> v))
                    throw runtime_error("malformed arrow at line " + to_string(line_no));
                if (v < 0)
                    throw runtime_error("negative arrow value at line " + to_string(line_no));
                ArrowKey key = {D, c};
                if (R.arrows.count(key))
                    throw runtime_error("duplicate arrow key at line " + to_string(line_no));
                if (v > 0) {
                    if (D < 2 || D >= G.H || c < 1 || c > D) {
                        throw runtime_error("arrow entry outside finite core support at line " +
                                            to_string(line_no));
                    }
                    R.arrows[key] = v;
                    R.a_star = max(R.a_star, v);
                }
            } else {
                throw runtime_error("unknown recipe directive '" + kind + "' at line " + to_string(line_no));
            }
            string trailing;
            if (iss >> trailing) {
                throw runtime_error("unexpected trailing token at line " + to_string(line_no));
            }
        }

        if (!has_scale) throw runtime_error("missing 'scale' line in recipe");
        if (!has_threshold) {
            if (!R.arrows.empty()) throw runtime_error("missing 'threshold' line in recipe");
            R.threshold = G.H;
        }
        for (int i = 0; i < n_c1_labels; ++i) {
            if (!seen_reward[i]) {
                throw runtime_error("missing reward for c=1 index " + to_string(i));
            }
        }
        return R;
    }
};

// -----------------------------------------------------------------------------
// Sliding-window banded verifier
// -----------------------------------------------------------------------------
struct BandedVerifier {
    Geometry G;
    Recipe R;
    i64 lam_num = 0;
    i64 lam_den = 1;
    double lam = 0.0;

    Layer c1_ref; // reference c=1 layer (identical across all bands)

    struct Slot { i64 in = 0, out = 0; };
    vector<vector<vector<Slot>>> slots; // slots[D][c][t] (allocated only where nonzero)
    vector<vector<i64>> arrow_tab;      // arrow_tab[D][c]
    vector<i64> suffix_arrow_max;       // suffix_arrow_max[m] = max_{m<=D'<H, 1<=k'<=D'} arrow(D',k')

    i64 tail_A = 0;
    i64 tail_B = 0;

    void prepare_global_tables() {
        c1_ref.init_labels(0, 1, G.h1, G.H);
        c1_ref.enumerate_bags();

        // Arrow lookup tables
        arrow_tab.assign(G.H, {});
        for (int D = 0; D < G.H; ++D) arrow_tab[D].assign(D + 1, 0);
        for (const auto& [key, val] : R.arrows) arrow_tab[key.first][key.second] = val;

        suffix_arrow_max.assign(G.H + 1, 0);
        for (int D = G.H - 1; D >= 2; --D) {
            i64 mx = 0;
            for (int k = 1; k <= D; ++k) mx = max(mx, arrow_tab[D][k]);
            suffix_arrow_max[D] = max(suffix_arrow_max[D + 1], mx);
        }

        // Build per-(band, c) label tables transiently to construct and cross-verify slots[D][c]
        int nb = (int)G.bands.size();
        vector<vector<Layer>> label_only(nb);
        for (int b = 0; b < nb; ++b) {
            int max_c = G.bands[b].d_end - 1;
            label_only[b].resize(max_c + 1);
            for (int c = 1; c <= max_c; ++c) {
                label_only[b][c].init_labels(b, c, G.bands[b].cap(c), G.H);
            }
        }

        slots.resize(G.H);
        for (int D = 2; D < G.H; ++D) {
            int b = G.band_of_D[D];
            slots[D].resize(D + 1);
            for (int c = 1; c <= D; ++c) {
                const Layer& Lyr = label_only[b][c];
                bool any = false;
                for (int t = 0; t < Lyr.pool; ++t) {
                    if (!Lyr.is_k_exact(t)) continue;
                    const Label& lb = Lyr.types[t];
                    if (R.edges.count({D, c, lb.d, lb.k}) || R.edges.count({lb.d, lb.k, D, c})) {
                        any = true;
                        break;
                    }
                }
                if (!any) continue;
                slots[D][c].assign(Lyr.types.size(), Slot{0, 0});
                for (int t = 0; t < Lyr.pool; ++t) {
                    if (!Lyr.is_k_exact(t)) continue;
                    const Label& lb = Lyr.types[t];
                    auto it_in = R.edges.find({D, c, lb.d, lb.k});
                    if (it_in != R.edges.end()) slots[D][c][t].in = it_in->second;
                    auto it_out = R.edges.find({lb.d, lb.k, D, c});
                    if (it_out != R.edges.end()) slots[D][c][t].out = it_out->second;
                }
            }
        }

        // Mirror cross-check across EVERY recipe edge and EVERY pooled/ambiguous-saturated slot.
        for (const auto& [key, val] : R.edges) {
            int D = key[0], c = key[1], d = key[2], k = key[3];
            const Layer& rec_l = label_only[G.band_of_D[D]][c];
            const Layer& snd_l = label_only[G.band_of_D[d]][k];
            int rec_hits = 0, snd_hits = 0;
            for (int t = 0; t < rec_l.pool; ++t) {
                const Label& lb = rec_l.types[t];
                if (lb.d == d && lb.k == k) {
                    if (!rec_l.is_k_exact(t) || slots[D][c].empty() || slots[D][c][t].in != val)
                        throw runtime_error("edge mirror mismatch on recipient side");
                    ++rec_hits;
                }
            }
            for (int t = 0; t < snd_l.pool; ++t) {
                const Label& lb = snd_l.types[t];
                if (lb.d == D && lb.k == c) {
                    if (!snd_l.is_k_exact(t) || slots[d][k].empty() || slots[d][k][t].out != val)
                        throw runtime_error("edge mirror mismatch on sender side");
                    ++snd_hits;
                }
            }
            if (rec_hits == 0 || snd_hits == 0)
                throw runtime_error("edge transfer failed to match both endpoints");
        }
        for (int D = 2; D < G.H; ++D) {
            int b = G.band_of_D[D];
            for (int c = 1; c <= D; ++c) {
                if (slots[D][c].empty()) continue;
                const Layer& Lyr = label_only[b][c];
                for (int t = 0; t < (int)Lyr.types.size(); ++t) {
                    if (!Lyr.is_k_exact(t) && (slots[D][c][t].in != 0 || slots[D][c][t].out != 0))
                        throw runtime_error("nonzero transfer on pooled/ambiguous-saturated slot");
                }
            }
        }

        // Analytic tail coefficients:
        //   B = max(0, floor((s - a_star) / (1 + lambda)))
        //   A = max(1, ceil((B + a_star) / lambda), ceil(s / lambda) - B)
        i64 s_trial = R.rewards.back();
        i64 a_st = R.a_star;
        i128 B0 = (s_trial > a_st)
            ? floor_div128((i128)lam_den * (s_trial - a_st), lam_num + lam_den)
            : (i128)0;
        tail_B = to_i64_checked(B0, "tail_B");
        i128 A1 = ceil_div128((i128)lam_den * ((i128)tail_B + a_st), lam_num);
        i128 A2 = ceil_div128((i128)lam_den * s_trial, lam_num) - tail_B;
        tail_A = max<i64>(1, to_i64_checked(max(A1, A2), "tail_A"));

        if (tail_A < 1 || tail_B < 0 ||
            (i128)lam_den * ((i128)tail_B + a_st) > (i128)lam_num * tail_A ||
            (i128)lam_den * s_trial > (i128)lam_num * ((i128)tail_A + tail_B)) {
            throw runtime_error("internal error: analytic tail inequalities failed at construction");
        }
    }

    void build_band(int b, vector<Layer>& band_layers, const vector<Layer>* higher_band) const {
        int max_c = G.bands[b].d_end - 1;
        band_layers.clear();
        band_layers.resize(max_c + 1);
        for (int c = 0; c <= max_c; ++c) {
            Layer& Lyr = band_layers[c];
            Lyr.init_labels(b, c, G.bands[b].cap(c), G.H);
            Lyr.enumerate_bags();

            // Precompute sender arrow maxima A_max[t]
            Lyr.arrow_max.assign(Lyr.types.size(), 0);
            for (int t = 0; t < (int)Lyr.types.size(); ++t) {
                const Label& lb = Lyr.types[t];
                if (t >= Lyr.pool) {
                    Lyr.arrow_max[t] = (lb.d < G.H) ? suffix_arrow_max[lb.d] : 0;
                } else if (lb.k < Lyr.k_limit) {
                    Lyr.arrow_max[t] = arrow_tab[lb.d][lb.k];
                } else {
                    i64 mx = 0;
                    for (int k2 = Lyr.k_limit; k2 <= lb.d; ++k2) mx = max(mx, arrow_tab[lb.d][k2]);
                    Lyr.arrow_max[t] = mx;
                }
            }
        }

        // Wire conversions, cross-band projections, and intra-layer moves
        vector<Run> key;
        for (int c = 0; c <= max_c; ++c) {
            Layer& Lyr = band_layers[c];

            // Conversion target in (b, c + 1)
            if (c + 1 <= max_c) {
                const Layer& nxt = band_layers[c + 1];
                Lyr.convert_tgt.resize(Lyr.n_rows, -1);
                int nb_pendant = nxt.label_plain(2, 1);
                for (int r = 0; r < Lyr.n_rows; ++r) {
                    key.clear();
                    key.push_back({nb_pendant, 1});
                    for (int ri = Lyr.run_off[r]; ri < Lyr.run_off[r + 1]; ++ri) {
                        const Run& rn = Lyr.runs[ri];
                        const Label& lb = Lyr.types[rn.type];
                        key.push_back({nxt.label_plain(lb.d, lb.k), rn.count});
                    }
                    Lyr.convert_tgt[r] = nxt.lookup_runs(key);
                }
            }

            // Cross-band gain target in (b + 1, c)
            if (higher_band != nullptr) {
                const Layer& up = (*higher_band)[c];
                Lyr.band_up_tgt.resize(Lyr.n_rows, -1);
                if (c == 0) {
                    Lyr.band_up_tgt[0] = 0;
                } else if (c == 1) {
                    if (Lyr.n_rows != up.n_rows)
                        throw runtime_error("internal error: c=1 row mismatch across bands");
                    for (int r = 0; r < Lyr.n_rows; ++r) Lyr.band_up_tgt[r] = r;
                } else {
                    for (int r = 0; r < Lyr.n_rows; ++r) {
                        key.clear();
                        for (int ri = Lyr.run_off[r]; ri < Lyr.run_off[r + 1]; ++ri) {
                            const Run& rn = Lyr.runs[ri];
                            const Label& lb = Lyr.types[rn.type];
                            key.push_back({up.label_plain(lb.d, lb.k), rn.count});
                        }
                        Lyr.band_up_tgt[r] = up.lookup_runs(key);
                    }
                }
            }

            // Neighbour & witness moves (skip entirely when pool == 0!)
            Lyr.move_off.assign(Lyr.n_rows + 1, 0);
            if (Lyr.pool > 0) {
                for (int r = 0; r < Lyr.n_rows; ++r) {
                    Lyr.move_off[r] = (int)Lyr.moves.size();
                    int rb = Lyr.run_off[r];
                    int re = Lyr.run_off[r + 1];
                    for (int ri = rb; ri < re; ++ri) {
                        const Run& rn = Lyr.runs[ri];
                        int t = rn.type;
                        if (t >= Lyr.pool) break; // sorted ascending; remaining runs are pooled
                        const Label& lb = Lyr.types[t];

                        auto make_target = [&](int new_t) -> int {
                            key.clear();
                            for (int q = rb; q < re; ++q) {
                                if (q == ri) {
                                    if (Lyr.runs[q].count > 1) key.push_back({t, Lyr.runs[q].count - 1});
                                    key.push_back({new_t, 1});
                                } else {
                                    key.push_back(Lyr.runs[q]);
                                }
                            }
                            return Lyr.lookup_runs(key);
                        };

                        // 1. Neighbour gains a leaf
                        int t_grow = lb.wd ? Lyr.label_bridge(lb.d + 1, lb.k, lb.wd, lb.wk)
                                           : Lyr.label_plain(lb.d + 1, lb.k);
                        Lyr.moves.push_back({make_target(t_grow), rn.count, t, 1});

                        // 2. Neighbour converts a leaf
                        if (lb.k < lb.d && lb.k < Lyr.k_limit) {
                            int t_conv = (c == 1 && lb.d <= BRIDGE_DEG &&
                                          lb.k + 1 >= 2 && lb.k + 1 <= WITNESS_CORE)
                                ? Lyr.label_bridge(lb.d, lb.k + 1, 2, 1)
                                : Lyr.label_plain(lb.d, lb.k + 1);
                            Lyr.moves.push_back({make_target(t_conv), rn.count, t, 0});
                        }

                        // 3. Exact witness moves (c == 1, wk > 0)
                        if (lb.wk > 0) {
                            int t_wg = Lyr.label_bridge(lb.d, lb.k, lb.wd + 1, lb.wk);
                            Lyr.moves.push_back({t_wg, 1, t, 3}); // for c==1, row idx == type idx
                            if (lb.wk < lb.wd) {
                                int t_wc = Lyr.label_bridge(lb.d, lb.k, lb.wd, lb.wk + 1);
                                Lyr.moves.push_back({t_wc, 1, t, 2});
                            }
                        }
                    }
                }
                Lyr.move_off[Lyr.n_rows] = (int)Lyr.moves.size();
            }
        }
    }

    int newborn_c1_index(int D, int c, const Layer& Lyr, int row_idx) const {
        if (D <= BRIDGE_DEG && D < G.h1 && c >= 1 && c < WITNESS_CORE) {
            int first_t = Lyr.runs[Lyr.run_off[row_idx]].type;
            const Label& first = Lyr.types[first_t];
            return c1_ref.label_bridge(D, c + 1, first.d, first.k);
        }
        return c1_ref.label_plain(D, c + 1);
    }

    i64 tail_upper_at_H(int L_next) const {
        return to_i64_checked((i128)tail_A + ceil_div128((i128)tail_B * L_next, G.H), "tail_upper_at_H");
    }

    // Independent exact 128-bit upward-rounded row residual: floor(lam * z) - drift_up(z).
    i128 exact_row_slack(int D, int c, const Layer& Lyr, int row_idx, i64 z,
                         i64 up_z, i64 conv_sum, const i64* z_cur_c) const {
        int L = D - c;
        int rb = Lyr.run_off[row_idx];
        int re = Lyr.run_off[row_idx + 1];
        const auto& sl = slots[D][c];

        // 1. Amended central gain change U = (up_z - z) + sum_w g(d_w, k_w; D, c)
        i128 U = (i128)up_z - z;
        if (!sl.empty()) {
            for (int ri = rb; ri < re; ++ri) {
                const Run& rn = Lyr.runs[ri];
                if (sl[rn.type].out != 0) U += (i128)rn.count * sl[rn.type].out;
            }
        }
        i64 a_here = arrow_tab[D][c];

        // Leaf gain term L * U
        i128 drift = (i128)L * U;

        // Core-neighbour credit terms with incoming arrow subtractions
        for (int ri = rb; ri < re; ++ri) {
            const Run& rn = Lyr.runs[ri];
            int deg_or_mark = Lyr.types[rn.type].d;
            i128 adj = U - ((deg_or_mark >= R.threshold) ? (i128)a_here : (i128)0);
            if (rn.type < Lyr.pool) {
                drift += ceil_div128(adj * rn.count, deg_or_mark);
            } else if (adj > 0) {
                drift += ceil_div128(adj * rn.count, deg_or_mark);
            }
        }

        // Outgoing arrow charge at sender D >= T
        if (D >= R.threshold && c > 0) {
            i128 arrow_sum = 0;
            for (int ri = rb; ri < re; ++ri) {
                const Run& rn = Lyr.runs[ri];
                arrow_sum += (i128)rn.count * Lyr.arrow_max[rn.type];
            }
            if (arrow_sum > 0) drift += ceil_div128(arrow_sum, D);
        }

        // 2. Central conversion term (rate L / D)
        if (L > 0) {
            drift += ceil_div128(((i128)conv_sum - z) * L, D);
        }

        // 3. Neighbour and witness transitions
        int m_beg = Lyr.move_off[row_idx];
        int m_end = Lyr.move_off[row_idx + 1];
        for (int mi = m_beg; mi < m_end; ++mi) {
            const Move& m = Lyr.moves[mi];
            const Label& lb = Lyr.types[m.type];
            i128 v = (i128)z_cur_c[m.target] - z;
            if (m.is_grow() && !m.is_witness() && !sl.empty()) {
                v -= sl[m.type].in;
            }
            if (m.is_witness()) {
                if (m.is_grow()) {
                    i64 twice = 2 * (lb.wd - lb.wk) + (v > 0 ? (lb.wk - 1) : 0);
                    drift += ceil_div128(v * (twice * lb.d + 2), 2 * lb.d);
                } else {
                    drift += ceil_div128(v * (lb.wd - lb.wk), lb.wd);
                }
            } else if (lb.wd > 0 && m.is_grow()) {
                i64 num_r = (v > 0)
                    ? ((2 * lb.d - lb.k - 2) * D * lb.wd + 2 * lb.wd + 2 * D)
                    : (2 * (lb.d - lb.k) * D * lb.wd + 2 * lb.wd + (lb.wk ? 2 * D : 0));
                drift += ceil_div128(v * num_r, 2 * D * lb.wd);
            } else if (m.is_grow()) {
                i64 twice_D_rate = 2;
                if (v > 0) twice_D_rate += (i64)(2 * lb.d - lb.k - 1) * D;
                else if (lb.k < Lyr.k_limit) twice_D_rate += (i64)(2 * (lb.d - lb.k)) * D;
                drift += ceil_div128(v * ((i128)twice_D_rate * m.count), 2 * D);
            } else {
                drift += ceil_div128(v * ((i128)m.count * (lb.d - lb.k)), lb.d);
            }
        }

        return floor_div128((i128)lam_num * z, lam_den) - drift;
    }

    i64 float_guess(int D, int c, const Layer& Lyr, int row_idx,
                    i64 up_z, i64 conv_sum, const i64* z_cur_c) const {
        int L = D - c;
        int rb = Lyr.run_off[row_idx];
        int re = Lyr.run_off[row_idx + 1];
        const auto& sl = slots[D][c];

        double out_sum = 0.0;
        if (!sl.empty()) {
            for (int ri = rb; ri < re; ++ri) {
                const Run& rn = Lyr.runs[ri];
                if (sl[rn.type].out != 0) out_sum += (double)rn.count * (double)sl[rn.type].out;
            }
        }
        double up_eff = (double)up_z + out_sum;
        double a_here = (double)arrow_tab[D][c];

        double sender_arrow = 0.0;
        if (D >= R.threshold && c > 0) {
            i128 as = 0;
            for (int ri = rb; ri < re; ++ri) {
                const Run& rn = Lyr.runs[ri];
                as += (i128)rn.count * Lyr.arrow_max[rn.type];
            }
            sender_arrow = (double)as / (double)D;
        }

        double conv_rate = (double)L / (double)D;
        double conv_term = (L > 0) ? (conv_rate * (double)conv_sum) : 0.0;

        int m_beg = Lyr.move_off[row_idx];
        int m_end = Lyr.move_off[row_idx + 1];

        double z = max(1.0, (double)up_z);
        for (int it = 0; it < 5; ++it) {
            double U = up_eff - z;
            double num = (double)L * up_eff + conv_term + sender_arrow;
            double den = lam + (double)L + conv_rate;
            for (int ri = rb; ri < re; ++ri) {
                const Run& rn = Lyr.runs[ri];
                int d_or_m = Lyr.types[rn.type].d;
                double sub = (d_or_m >= R.threshold) ? a_here : 0.0;
                if (rn.type < Lyr.pool || (U - sub) > 0.0) {
                    double r_c = (double)rn.count / (double)d_or_m;
                    num += r_c * (up_eff - sub);
                    den += r_c;
                }
            }
            for (int mi = m_beg; mi < m_end; ++mi) {
                const Move& m = Lyr.moves[mi];
                const Label& lb = Lyr.types[m.type];
                double dst = (double)z_cur_c[m.target];
                if (m.is_grow() && !m.is_witness() && !sl.empty()) dst -= (double)sl[m.type].in;
                bool pos = (dst > z);
                double rate = 0.0;
                if (m.is_witness()) {
                    rate = m.is_grow()
                        ? (lb.wd - lb.wk + 1.0 / lb.d + (pos ? 0.5 * (lb.wk - 1) : 0.0))
                        : ((double)(lb.wd - lb.wk) / lb.wd);
                } else if (lb.wd > 0 && m.is_grow()) {
                    rate = (lb.d - lb.k) + 1.0 / D +
                           (pos ? (1.0 / lb.wd + 0.5 * (lb.k - 2)) : (lb.wk ? 1.0 / lb.wd : 0.0));
                } else if (m.is_grow()) {
                    double lo = (lb.k < Lyr.k_limit) ? (lb.d - lb.k) : 0.0;
                    double hi = lb.d - 0.5 * (lb.k + 1);
                    rate = m.count * (1.0 / D + (pos ? hi : lo));
                } else {
                    rate = m.count * ((double)(lb.d - lb.k) / lb.d);
                }
                num += rate * dst;
                den += rate;
            }
            double nz = max(1.0, num / den);
            if (fabs(nz - z) < 0.25) { z = nz; break; }
            z = nz;
        }
        if (!isfinite(z) || z > 1e17) return max<i64>(1, R.scale);
        return max<i64>(1, llround(z));
    }

    i64 solve_minimal_z(int D, int c, const Layer& Lyr, int row_idx,
                        i64 up_z, i64 conv_sum, const i64* z_cur_c) const {
        i64 z = float_guess(D, c, Lyr, row_idx, up_z, conv_sum, z_cur_c);
        i128 s = exact_row_slack(D, c, Lyr, row_idx, z, up_z, conv_sum, z_cur_c);
        if (s >= 0) {
            if (z == 1) return 1;
            if (exact_row_slack(D, c, Lyr, row_idx, z - 1, up_z, conv_sum, z_cur_c) < 0) return z;
            i64 step = 1, hi = z, lo = max<i64>(1, hi - step);
            while (lo > 1 && exact_row_slack(D, c, Lyr, row_idx, lo, up_z, conv_sum, z_cur_c) >= 0) {
                hi = lo;
                step = min<i64>(step * 2, hi - 1);
                lo = max<i64>(1, hi - step);
            }
            if (lo == 1 && exact_row_slack(D, c, Lyr, row_idx, 1, up_z, conv_sum, z_cur_c) >= 0) return 1;
            while (hi - lo > 1) {
                i64 mid = lo + (hi - lo) / 2;
                if (exact_row_slack(D, c, Lyr, row_idx, mid, up_z, conv_sum, z_cur_c) >= 0) hi = mid;
                else lo = mid;
            }
            return hi;
        } else {
            i64 lo = z;
            i64 step = max<i64>(1, to_i64_checked(ceil_div128((-s) * lam_den, lam_num), "step"));
            i64 hi = lo + step;
            while (exact_row_slack(D, c, Lyr, row_idx, hi, up_z, conv_sum, z_cur_c) < 0) {
                lo = hi;
                if (step > (numeric_limits<i64>::max() - hi) / 2)
                    throw runtime_error("weight overflow while searching for supersolution");
                step *= 2;
                hi = lo + step;
            }
            while (hi - lo > 1) {
                i64 mid = lo + (hi - lo) / 2;
                if (exact_row_slack(D, c, Lyr, row_idx, mid, up_z, conv_sum, z_cur_c) >= 0) hi = mid;
                else lo = mid;
            }
            return hi;
        }
    }

    int run() {
        auto t0 = chrono::steady_clock::now();
        prepare_global_tables();

        cout << "geometry: " << G.bands.size() << " bands, H=" << G.H
             << ", h1=" << G.h1 << ", c1_labels=" << c1_ref.n_rows << "\n";
        cout << "recipe: scale=" << R.scale << ", threshold=" << R.threshold
             << ", edges=" << R.edges.size() << ", arrows=" << R.arrows.size()
             << ", a_star=" << R.a_star << "\n";
        if (R.used_boundary_cap_edge) {
            cout << "note: recipe includes nonzero edge(s) at exact degree==core_cap boundary\n";
        }
        cout << "tail: A=" << tail_A << ", B=" << tail_B
             << " (for trial s=" << R.rewards.back() << ")\n" << flush;

        int nb = (int)G.bands.size();
        vector<Layer> cur_band, next_band;
        vector<vector<i64>> z_cur, z_nxt; // indexed by c = 0..D

        i64 min_z = numeric_limits<i64>::max();
        i64 max_z = 0;
        i128 min_slack = numeric_limits<i64>::max();
        int64_t total_states = 0;
        i64 root_z = 0, star_z = 0, p4_z = 0;
        vector<i64> actual_rewards;

        for (int b = nb - 1; b >= 0; --b) {
            auto tb = chrono::steady_clock::now();
            build_band(b, cur_band, (b + 1 < nb) ? &next_band : nullptr);
            int64_t band_states = 0;

            for (int D = G.bands[b].d_end - 1; D >= G.bands[b].d_start; --D) {
                bool at_band_top = (D == G.bands[b].d_end - 1);
                bool into_tail = (D + 1 == G.H);
                z_cur.assign(D + 1, {});

                for (int c = D; c >= 0; --c) {
                    int L = D - c;
                    const Layer& Lyr = cur_band[c];
                    z_cur[c].assign(Lyr.n_rows, 0);
                    i64* z_c = z_cur[c].data();
                    const i64* z_c1 = (L > 0) ? z_cur[c + 1].data() : nullptr;
                    const i64* z_up_arr = (!into_tail) ? z_nxt[c].data() : nullptr;
                    i64 tail_up = into_tail ? tail_upper_at_H(L + 1) : 0;

                    for (int r = Lyr.n_rows - 1; r >= 0; --r) {
                        i64 up_z = into_tail
                            ? tail_up
                            : z_up_arr[at_band_top ? Lyr.band_up_tgt[r] : r];
                        i64 conv_sum = 0;
                        if (L > 0) {
                            conv_sum = to_i64_checked(
                                (i128)z_c1[Lyr.convert_tgt[r]] +
                                R.rewards[newborn_c1_index(D, c, Lyr, r)],
                                "conv_sum");
                        }

                        i64 z = solve_minimal_z(D, c, Lyr, r, up_z, conv_sum, z_c);
                        // Independent check of exact 128-bit inequality on the finalized z
                        i128 slk = exact_row_slack(D, c, Lyr, r, z, up_z, conv_sum, z_c);
                        if (z < 1 || slk < 0) {
                            throw runtime_error("internal verification failure at D=" + to_string(D) +
                                                " c=" + to_string(c) + " r=" + to_string(r));
                        }
                        z_c[r] = z;
                        min_z = min(min_z, z);
                        max_z = max(max_z, z);
                        if (slk < min_slack) min_slack = slk;
                        ++band_states;
                        ++total_states;
                    }
                }

                if (D == 3) star_z = z_cur[0][0];
                if (D == 2) {
                    root_z = z_cur[0][0];
                    actual_rewards = z_cur[1];
                    p4_z = z_cur[1][c1_ref.label_plain(2, 1)];
                }

                z_nxt.swap(z_cur);
                z_cur.clear();
                // Once we have processed the single cross-band degree D = d_end[b] - 1, free next_band!
                if (at_band_top && !next_band.empty()) {
                    next_band.clear();
                    next_band.shrink_to_fit();
                }
            }

            double dt_b = chrono::duration<double>(chrono::steady_clock::now() - tb).count();
            cout << "band " << b << " [" << G.bands[b].d_start << ", " << G.bands[b].d_end
                 << "): " << band_states << " states verified in " << fixed << setprecision(2)
                 << dt_b << "s\n" << flush;

            next_band.swap(cur_band);
            cur_band.clear();
        }

        // Verify newborn consistency: actual_rewards[i] <= R.rewards[i] for all c=1 rows
        int n_c1 = c1_ref.n_rows;
        if ((int)actual_rewards.size() != n_c1) throw runtime_error("missing D=2, c=1 rewards");
        i64 min_nb_slack = numeric_limits<i64>::max();
        int worst_nb_idx = -1;
        for (int i = 0; i < n_c1; ++i) {
            i64 gap = R.rewards[i] - actual_rewards[i];
            if (gap < min_nb_slack) {
                min_nb_slack = gap;
                worst_nb_idx = i;
            }
        }
        cout << "newborn consistency: min(trial - actual) = " << min_nb_slack
             << " (worst c=1 index " << worst_nb_idx
             << ": actual=" << actual_rewards[worst_nb_idx]
             << ", trial=" << R.rewards[worst_nb_idx] << ")\n";
        if (min_nb_slack < 0) {
            cout << "FAIL: newborn consistency failed (actual weight exceeds trial reward)\n";
            return 1;
        }

        // Verify analytic tail inequalities against ACTUAL pooled-H newborn weight!
        i64 s_actual = actual_rewards.back();
        if (tail_A < 1 || tail_B < 0 ||
            (i128)lam_den * ((i128)tail_B + R.a_star) > (i128)lam_num * tail_A ||
            (i128)lam_den * s_actual > (i128)lam_num * ((i128)tail_A + tail_B)) {
            cout << "FAIL: analytic tail check failed against actual pooled-H newborn weight\n";
            return 1;
        }

        i64 global_min = min(min_z, tail_A);
        i64 global_max = max(max_z, to_i64_checked((i128)tail_A + tail_B, "global_max"));
        i128 drift3 = 2 * (i128)star_z + 2 * (i128)p4_z - 3 * (i128)root_z;
        double dt_tot = chrono::duration<double>(chrono::steady_clock::now() - t0).count();

        cout << "verified " << total_states << " finite states + analytic D>=" << G.H
             << " tail in " << fixed << setprecision(2) << dt_tot
             << "s; min finite slack = " << (long long)min_slack << "\n";
        cout << "weights (in recipe units): finite [" << min_z << ", " << max_z
             << "], tail [" << tail_A << ", " << (tail_A + tail_B)
             << "], global [" << global_min << ", " << global_max << "]\n";
        cout << "T3 exact scaled drift / root weight = " << (long long)drift3 << "/" << root_z
             << " (=" << setprecision(12) << (double)drift3 / (double)root_z << "); "
             << "30*drift - 17*root = " << (long long)(30 * drift3 - 17 * (i128)root_z) << "\n";
        cout << "PASS upper " << lam_num << "/" << lam_den
             << " prefactor " << root_z << "/" << global_min
             << " (=" << setprecision(12) << (double)root_z / (double)global_min << ")\n";
        return 0;
    }
};

int main(int argc, char** argv) {
    if (argc != 5) {
        cerr << "usage: " << argv[0] << " <lambda_num> <lambda_den> <bands_file> <recipe_file>\n";
        return 2;
    }
    try {
        BandedVerifier V;
        V.lam_num = stoll(argv[1]);
        V.lam_den = stoll(argv[2]);
        if (V.lam_num <= 0 || V.lam_den <= 0 || V.lam_num >= V.lam_den) {
            throw runtime_error("lambda must lie strictly in (0, 1)");
        }
        V.lam = (double)V.lam_num / (double)V.lam_den;
        V.G = Geometry::load(argv[3]);

        Layer tmp_c1;
        tmp_c1.init_labels(0, 1, V.G.h1, V.G.H);
        tmp_c1.enumerate_bags();
        V.R = Recipe::load(argv[4], V.G, tmp_c1.n_rows);

        return V.run();
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}
