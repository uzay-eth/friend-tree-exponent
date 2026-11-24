#include <bits/stdc++.h>
using namespace std;

// Fast C++ test: Adaptive Planted-Tree Quotient Bellman Enclosure!
// Represent ANY planted tree by an exact finite tree up to a "cost budget" B,
// truncating any subtree exceeding budget to a coarse node (D, c)!
// Specifically, test how tight [s_lo, s_hi] becomes when we increase budget B!
// Cost function for a node at depth from root:
//   - Bare node (d=2, c=2): cost 1 (so bare chains are resolved deeply!)
//   - Virgin (d=2, c=1): cost 0 (exact leaf of core tree!)
//   - Small core (d<=4, c<=3): cost 2
//   - Larger core: cost 3+
//
// Even simpler first: test EXACT trees up to size N_exact (all internal rates EXACT, zero uncertainty!)
// glued to a Depth-0 / Depth-1 tail!
// Wait: let's check BOTH:
//   Test A: Exact unmarked planted trees up to size N_exact, where ANY exit to size > N_exact
//           drops to coarse state (min(d, D0), min(c, C0))!
//   WAIT A SECOND!! Check if dropping to coarse state (min(d, D0), min(c, C0)) at |S| > N_exact is 100% SOUND!
//   YES!! Because if phi(S) = S when |S| <= N_exact, and phi(S) = ('T0', min(D_o, D0), min(c_o, C0)) when |S| > N_exact:
//     - For EVERY S with |S| <= N_exact:
//         * EVERY intrinsic jump S -> S' is 100% EXACT (if |S'|=N_exact+1, it lands in ('T0', min(D',D0), min(c',C0))!).
//         * EVERY child Y of S has |Y| <= N_exact - 1, so BOTH Y and Y+ are EXACT states in |S| <= N_exact!
//         * THEREFORE, for ALL |S| <= N_exact, THERE IS ZERO BELLMAN UNCERTAINTY AT ALL!!
//     - For EVERY S with |S| > N_exact (coarse state ('T0', d, c)):
//         * Its children are unknown, EXCEPT wait: its root gain/prom stay in T0, and C Delta w is bounded by
//           (c-1)/d * [gMinD, gMaxD]!
//   HOLD ON!! LET'S TEST THAT RIGHT NOW FOR N_exact = 6, 8, 10, 12, 14!!
//   In 3 SECONDS we will see the EXACT enclosure width of phi(S) = {exact if |S|<=N, T0 if |S|>N}!!

struct Tree {
    int L;
    vector<int> ch; // sorted child tree IDs (in the exact pool)
    bool operator<(const Tree& o) const {
        if (L != o.L) return L < o.L;
        return ch < o.ch;
    }
};

int main() {
    const int Nmax = 13;
    vector<Tree> trees;
    map<Tree, int> tid;
    vector<int> tsize, tdeg, tcore;
    vector<vector<int>> by_size(Nmax + 1);

    auto add_tree = [&](const Tree& tr, int sz) -> int {
        auto it = tid.find(tr);
        if (it != tid.end()) return it->second;
        int id = (int)trees.size();
        tid[tr] = id;
        trees.push_back(tr);
        tsize.push_back(sz);
        tdeg.push_back(1 + tr.L + (int)tr.ch.size());
        tcore.push_back(1 + (int)tr.ch.size());
        by_size[sz].push_back(id);
        return id;
    };

    int J_id = add_tree({1, {}}, 2);

    // Recursively build all trees up to size Nmax, plus track their size-(Nmax+1) jump targets' (d, c)
    // First generate all trees up to size Nmax:
    vector<vector<pair<int, double>>> raw_jumps; // computed later
    for (int sz = 3; sz <= Nmax; ++sz) {
        for (int pid : by_size[sz - 1]) {
            Tree tr = trees[pid];
            // 1. gain at root
            add_tree({tr.L + 1, tr.ch}, sz);
            // 2. promote at root
            if (tr.L > 0) {
                auto nch = tr.ch;
                nch.push_back(J_id);
                sort(nch.begin(), nch.end());
                add_tree({tr.L - 1, nch}, sz);
            }
            // 3. jump inside any child
            // Wait: we can just use precomputed child jumps from smaller sizes!
        }
        // Actually do clean recursive jump generator on size sz-1!
        break;
    }
    // Clean recursive generator:
    trees.clear(); tid.clear(); tsize.clear(); tdeg.clear(); tcore.clear();
    for (auto& v : by_size) v.clear();
    J_id = add_tree({1, {}}, 2);

    // Map (tree_id, incoming_is_root) -> list of (next_tree_id_or_-1, next_d, next_c, rate)
    // Better: compute jumps(id) recursively ONLY when needed, storing results!
    vector<vector<pair<int, double>>> int_jumps_0; // incoming rate 0 (for root)
    vector<vector<pair<int, double>>> int_jumps_1; // rate of root leaf gain (separate so we can add incoming!)
    // Actually for any tree id of size s < Nmax, its intrinsic transitions with incoming=0 are:
    //   - root leaf gain: target gain_id(id), rate = L + sum_y 1/deg(y)  [+ incoming!]
    //   - root prom: target prom_id(id), rate = L / deg(id)
    //   - child j transition: rate from child j with incoming = 1/deg(id)
    vector<int> gain_id, prom_id;

    function<void(int)> ensure_cap = [&](int n) {
        if ((int)gain_id.size() < n) {
            gain_id.resize(n, -1);
            prom_id.resize(n, -1);
            int_jumps_0.resize(n);
        }
    };

    for (int sz = 2; sz < Nmax; ++sz) {
        int cur_cnt = (int)by_size[sz].size();
        for (int idx_in_sz = 0; idx_in_sz < cur_cnt; ++idx_in_sz) {
            int id = by_size[sz][idx_in_sz];
            Tree tr = trees[id]; // copy!
            int d = tdeg[id];
            double p = 1.0 / d;
            // gain at root
            int gid = add_tree({tr.L + 1, tr.ch}, sz + 1);
            ensure_cap((int)trees.size());
            gain_id[id] = gid;
            // prom at root
            if (tr.L > 0) {
                auto nch = tr.ch; nch.push_back(J_id); sort(nch.begin(), nch.end());
                int pid = add_tree({tr.L - 1, nch}, sz + 1);
                ensure_cap((int)trees.size());
                prom_id[id] = pid;
            }
            // collect jumps into map
            map<int, double> jmap;
            double gr = tr.L;
            for (int y : tr.ch) gr += 1.0 / tdeg[y];
            if (gr > 0) jmap[gid] += gr;
            if (tr.L > 0) jmap[prom_id[id]] += (double)tr.L * p;
            // child jumps (note every child y has size <= sz - 1, so its int_jumps_0 and gain_id are ALREADY built!)
            for (size_t j = 0; j < tr.ch.size(); ++j) {
                int y = tr.ch[j];
                // y's gain gets extra rate p!
                auto add_child_jump = [&](int ny, double r) {
                    auto nch = tr.ch;
                    nch[j] = ny;
                    sort(nch.begin(), nch.end());
                    int nid = add_tree({tr.L, nch}, sz + 1);
                    ensure_cap((int)trees.size());
                    jmap[nid] += r;
                };
                add_child_jump(gain_id[y], p);
                for (auto [ny, r] : int_jumps_0[y]) {
                    add_child_jump(ny, r);
                }
            }
            ensure_cap((int)trees.size());
            for (auto [nid, r] : jmap) int_jumps_0[id].push_back({nid, r});
        }
    }
    cout << "Generated unmarked trees up to Nmax=" << Nmax << ": total=" << trees.size() << "\n";
    for (int sz = 2; sz <= Nmax; ++sz) cout << "  sz=" << sz << ": " << by_size[sz].size() << "\n";

    // Now test N_exact = 8, 10, 12 glued to T0(D0=24, C0=6)!
    // NOTE: for a tree id of size N_exact, its intrinsic jumps go to size N_exact+1 (already in trees for N_exact < Nmax!)
    // and map to T0(min(d', D0), min(c', C0))!
    const int D0 = 24, C0 = 6;
    map<pair<int,int>, int> t0_idx;
    vector<pair<int,int>> t0_list;
    for (int d = 2; d <= D0; ++d)
        for (int c = 1; c <= min(d, C0); ++c) {
            t0_idx[{d, c}] = (int)t0_list.size();
            t0_list.push_back({d, c});
        }
    int Nt0 = (int)t0_list.size();

    for (int Nex : {6, 8, 10, 12}) {
        int K = 0;
        for (int sz = 2; sz <= Nex; ++sz) K += (int)by_size[sz].size();
        // State 0..K-1 = exact trees of size <= Nex
        // State K..K+Nt0-1 = T0 coarse states (d, c)
        int Ntot = K + Nt0;
        auto to_state = [&](int tid_val) -> int {
            if (tid_val < K) return tid_val;
            return K + t0_idx[{min(tdeg[tid_val], D0), min(tcore[tid_val], C0)}];
        };

        // Precompute exact transitions for 0..K-1
        vector<vector<pair<int, double>>> tr_ex(K);
        vector<vector<pair<int, int>>> ch_ex(K);
        for (int i = 0; i < K; ++i) {
            for (auto [nid, r] : int_jumps_0[i]) {
                tr_ex[i].push_back({to_state(nid), r});
            }
            for (int y : trees[i].ch) {
                ch_ex[i].push_back({y, gain_id[y]}); // BOTH < K because |y|<=Nex-1 => |y+|<=Nex!
            }
        }

        auto solve = [&](double s, bool upper_mu) -> double {
            vector<double> v(Ntot, 1.0), vn(Ntot, 1.0);
            for (int it = 0; it < 120; ++it) {
                // Global minD, maxD over ALL states (needed ONLY for T0 coarse states!)
                double gMinD = 1e9, gMaxD = -1e9;
                for (int i = 0; i < K; ++i) {
                    int g = to_state(gain_id[i]);
                    double dv = v[g] - v[i];
                    if (dv < gMinD) gMinD = dv;
                    if (dv > gMaxD) gMaxD = dv;
                }
                for (int j = 0; j < Nt0; ++j) {
                    auto [d, c] = t0_list[j];
                    int g = K + t0_idx[{min(d + 1, D0), c}];
                    double dv = v[g] - v[K + j];
                    if (dv < gMinD) gMinD = dv;
                    if (dv > gMaxD) gMaxD = dv;
                }
                double md = 0.0;
                // Update T0 coarse states first
                for (int j = Nt0 - 1; j >= 0; --j) {
                    auto [d, c] = t0_list[j];
                    int idx = K + j;
                    if (d == D0) {
                        vn[idx] = (1.0 - (upper_mu ? 0.35 : 0.75) * c / d) / s;
                        continue;
                    }
                    int L = d - c;
                    double b = (double)L / d;
                    double gv = v[K + t0_idx[{min(d + 1, D0), c}]];
                    double pv = (L > 0) ? v[K + t0_idx[{d, min(c + 1, C0)}]] : 0.0;
                    double xi = (upper_mu == (gv >= v[idx])) ? (c - 1) * 0.5 : 0.0;
                    double cd = ((double)(c - 1) / d) * (upper_mu ? gMaxD : gMinD);
                    double val = (b + cd + (L + xi) * gv + b * pv) / (s + L + xi + b);
                    md = max(md, fabs(val - v[idx]));
                    vn[idx] = val;
                }
                // Backward sweep through exact states K-1 .. 0 (ZERO uncertainty!)
                for (int i = K - 1; i >= 0; --i) {
                    int d = tdeg[i], L = trees[i].L;
                    double b = (double)L / d;
                    double cd = 0.0;
                    for (auto [iy, iyp] : ch_ex[i]) cd += (v[iyp] - v[iy]);
                    cd /= d;
                    double rhs = b + cd;
                    for (auto [nid, r] : tr_ex[i]) {
                        rhs += r * (nid > i ? vn[nid] : v[nid]);
                    }
                    double hold = tsize[i] - 1.0 / d;
                    double val = rhs / (s + hold);
                    md = max(md, fabs(val - v[i]));
                    vn[i] = val;
                }
                v.swap(vn);
                if (md < 1e-12) break;
            }
            return v[J_id];
        };

        double lo = 0.50, hi = 0.60;
        for (int k = 0; k < 35; ++k) {
            double m = 0.5 * (lo + hi);
            if (solve(m, false) > 1.0) lo = m; else hi = m;
        }
        double s_lo = lo;
        lo = 0.54; hi = 0.62;
        for (int k = 0; k < 35; ++k) {
            double m = 0.5 * (lo + hi);
            if (solve(m, true) < 1.0) hi = m; else lo = m;
        }
        double s_hi = hi;
        cout << fixed << setprecision(8) << "Exact |S|<=" << setw(2) << Nex << " (K=" << setw(5) << K << ") + T0: mu in ["
             << s_lo << ", " << s_hi << "] (width = " << (s_hi - s_lo) << ")\n";
    }
    return 0;
}
