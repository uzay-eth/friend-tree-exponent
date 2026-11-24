// ftvec.cpp -- training tool: simulate independent continuous-time friend trees and write, per run and checkpoint,
// the exact sums a = sum_u a_u and b = sum_u b_u of the feature generator action (feat2.hpp, analytic version).
// Usage: ftvec SEED RUNS OUT T1 [T2 ...] [key=value ...]
//   run r (run0 <= r < run0+RUNS, run0 given as run0=..., default 0) uses the stream Rng::seed2(SEED, r).
//   diffout=1 (two checkpoints): write only the difference record (t2 minus t1) per run.  Output OUT.meta.txt and OUT.rec.bin (float64; per run, per checkpoint:
//   [C, S, a(p), b(p)]),  p = number of features of the layout.
#include "feat2.hpp"
#include <cstring>
static void tree_to_snap(Tree &tr, Snap &T) { T.L.swap(tr.L); T.c.swap(tr.c); T.parent.swap(tr.parent); T.first_child.swap(tr.first_child); T.next_sib.swap(tr.next_sib); }
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s SEED RUNS OUT T1 [T2..] [key=value..]\n", argv[0]); return 1; }
    uint64_t seed = strtoull(argv[1], 0, 10); long runs = atol(argv[2]); std::string pre(argv[3]);
    std::vector<double> cps; F2Cfg G; int diffout = 0; long run0 = 0; int model = 0;
    for (int i = 4; i < argc; i++) {
        char *eq = strchr(argv[i], '=');
        if (!eq) { cps.push_back(atof(argv[i])); continue; }
        std::string k(argv[i], eq - argv[i]);
        if (k == "diffout") { diffout = atoi(eq + 1); continue; }
        if (k == "model") { model = atoi(eq + 1); continue; }
        if (k == "run0") { run0 = atol(eq + 1); continue; }
        if (!G.set(k, atoi(eq + 1))) { fprintf(stderr, "unknown key %s\n", argv[i]); return 1; }
    }
    if (diffout && cps.size() != 2) { fprintf(stderr, "diffout=1 needs exactly two checkpoints\n"); return 1; }
    Layout Y(G); size_t p = Y.total, ncp = cps.size(), rec = 2 + 2 * p;
    std::vector<double> buf(ncp * rec);
    FILE *fo = fopen((pre + ".rec.bin").c_str(), "wb"); if (!fo) { perror("out"); return 1; }
    Tree tr; Rng g; Snap T;
    for (long r = run0; r < run0 + runs; r++) {
        g.seed2(seed, (uint64_t)r); tr.init(); tr.model = model;
        std::fill(buf.begin(), buf.end(), 0.0);
        for (size_t k = 0; k < ncp; k++) {
            tr.run_to(cps[k], g);
            tree_to_snap(tr, T); T.build_cred();
            FastGen fg(Y, T, model);
            double *b = &buf[k * rec]; VecSink sa{b + 2}, sb{b + 2 + p};
            double S = 0; size_t C = T.size();
            for (uint32_t u = 0; u < C; u++) { S += fg.conv[u]; fg.compute(u, sa, sb); }
            b[0] = (double)C; b[1] = S;
            tree_to_snap(tr, T);
        }
        if (diffout) { for (size_t i = 0; i < rec; i++) buf[i] = buf[rec + i] - buf[i]; fwrite(buf.data(), 8, rec, fo); }
        else fwrite(buf.data(), 8, buf.size(), fo);
    }
    fclose(fo);
    FILE *fm = fopen((pre + ".meta.txt").c_str(), "w");
    fprintf(fm, "%zu %zu %ld %ld %d\n", p, diffout ? (size_t)1 : ncp, runs, runs, diffout);
    for (double t : cps) fprintf(fm, "%.6f ", t);
    fprintf(fm, "\n"); fclose(fm);
    return 0;
}
