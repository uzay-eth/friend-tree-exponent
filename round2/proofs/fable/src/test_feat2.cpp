// test_feat2.cpp -- check that the analytic generator action (FastGen) equals the brute-force event replay (gen_bf)
// for every vertex of simulated trees.  Usage: test_feat2 SEED RUNS T [key=value ...]
#include "feat2.hpp"
#include <cstring>
#include <string>
static void tree_to_snap(Tree &tr, Snap &T) { T.L.swap(tr.L); T.c.swap(tr.c); T.parent.swap(tr.parent); T.first_child.swap(tr.first_child); T.next_sib.swap(tr.next_sib); }
int main(int argc, char **argv) {
    uint64_t seed = strtoull(argv[1], 0, 10); long runs = atol(argv[2]); double Tm = atof(argv[3]);
    F2Cfg G; int model = 0;
    for (int i = 4; i < argc; i++) { char *eq = strchr(argv[i], '='); if (!eq) continue; std::string k(argv[i], eq - argv[i]); if (k == "model") { model = atoi(eq + 1); continue; } if (!G.set(k, atoi(eq + 1))) { fprintf(stderr, "bad key\n"); return 1; } }
    Layout Y(G); int p = Y.total;
    Tree tr; Rng g; Snap T;
    std::vector<double> a1(p), b1(p), a2(p), b2(p); double maxd = 0, maxv = 0; long nv = 0, nnz = 0;
    std::vector<double> used(p, 0.0);
    for (long r = 0; r < runs; r++) {
        g.seed2(seed, (uint64_t)r); tr.init(); tr.model = model; tr.run_to(Tm, g);
        tree_to_snap(tr, T); T.build_cred();
        FastGen fg(Y, T, model);
        for (uint32_t u = 0; u < T.size(); u++) {
            std::fill(a1.begin(), a1.end(), 0); std::fill(b1.begin(), b1.end(), 0); std::fill(a2.begin(), a2.end(), 0); std::fill(b2.begin(), b2.end(), 0);
            VecSink sa1{a1.data()}, sb1{b1.data()}, sa2{a2.data()}, sb2{b2.data()};
            gen_bf(Y, T, u, sa1, sb1, model);
            fg.compute(u, sa2, sb2);
            for (int i = 0; i < p; i++) {
                double da = fabs(a1[i] - a2[i]), db = fabs(b1[i] - b2[i]);
                double sc = 1.0 + fabs(a1[i]);
                if (da / sc > maxd) maxd = da / sc; if (db > maxd) maxd = db;
                if (fabs(a1[i]) > maxv) maxv = fabs(a1[i]);
                if (a1[i] != 0 || b1[i] != 0) { nnz++; used[i] += 1; }
            }
            nv++;
        }
        tree_to_snap(tr, T);
    }
    int nused = 0; for (int i = 0; i < p; i++) if (used[i] > 0) nused++;
    printf("features p=%d (used %d)  vertices checked=%ld  nonzeros=%ld  max |a|=%.3g  max scaled difference=%.3e\n", p, nused, nv, nnz, maxv, maxd);
    return maxd < 1e-9 ? 0 : 2;
}
