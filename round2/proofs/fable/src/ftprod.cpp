// ftprod.cpp -- production tool: simulate independent continuous-time friend trees; at each checkpoint evaluate, for
// each given FIXED weight vector theta (feature layout of feat2.hpp + plain equation as last entry),
//     num = sum_u theta.a_u ,  den = sum_u theta.b_u ,
// where a_u, b_u are the exact generator action and value of the feature vector (analytic, feat2.hpp).
// For every fixed theta,  d/dt E[den] = E[num]  exactly; in the growth regime E[num] = mu E[den].
// Usage: ftprod SEED RUN0 RUNS OUT THETA1[,THETA2,...] T1 [T2 ...] [key=value ...] [diag=1]
//   run r (RUN0 <= r < RUN0+RUNS) uses the stream Rng::seed2(SEED, r).
//   OUT.txt: one line per run and checkpoint:  r t C S  num_1 den_1  num_2 den_2 ...
//   diag=1: additionally OUT.diag.txt with, for theta_1 at the last checkpoint, per own state (min(L,40), min(c,12)):
//           count, sum psi, sum psi^2  where psi_u = theta.(a_u - mu0 b_u), mu0 = 0.565723.
#include "feat2.hpp"
#include <cstring>
static void tree_to_snap(Tree &tr, Snap &T) { T.L.swap(tr.L); T.c.swap(tr.c); T.parent.swap(tr.parent); T.first_child.swap(tr.first_child); T.next_sib.swap(tr.next_sib); }
struct DotSink { const double *th; double acc; inline void add(int i, double x) { acc += th[i] * x; } };
struct MultiSink { int n; const double *th[8]; double acc[8]; inline void add(int i, double x) { for (int k = 0; k < n; k++) acc[k] += th[k][i] * x; } };
int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s SEED RUN0 RUNS OUT THETAS T1 [T2..] [key=value..]\n", argv[0]); return 1; }
    uint64_t seed = strtoull(argv[1], 0, 10); long run0 = atol(argv[2]), runs = atol(argv[3]); std::string pre(argv[4]);
    std::vector<std::string> thf; { std::string s(argv[5]); size_t pos = 0; while (true) { size_t q = s.find(',', pos); thf.push_back(s.substr(pos, q - pos)); if (q == std::string::npos) break; pos = q + 1; } }
    std::vector<double> cps; F2Cfg G; int diag = 0; int model = 0; int rfast = 1;
    for (int i = 6; i < argc; i++) {
        char *eq = strchr(argv[i], '=');
        if (!eq) { cps.push_back(atof(argv[i])); continue; }
        std::string k(argv[i], eq - argv[i]);
        if (k == "diag") { diag = atoi(eq + 1); continue; }
        if (k == "model") { model = atoi(eq + 1); continue; }
        if (k == "rfast") { rfast = atoi(eq + 1); continue; }
        if (!G.set(k, atoi(eq + 1))) { fprintf(stderr, "unknown key %s\n", argv[i]); return 1; }
    }
    Layout Y(G); size_t p = Y.total; int nth = (int)thf.size(); if (nth > 8) return 1;
    std::vector<std::vector<double>> th(nth, std::vector<double>(p + 1));
    for (int k = 0; k < nth; k++) { FILE *f = fopen(thf[k].c_str(), "rb"); if (!f || fread(th[k].data(), 8, p + 1, f) != p + 1) { fprintf(stderr, "bad theta file %s (layout p=%zu)\n", thf[k].c_str(), p); return 1; } fclose(f); }
    const double *thp[8]; for (int k = 0; k < nth; k++) thp[k] = th[k].data();
    bool useR = rfast && !diag && Y.NR > 0 && G.Rpow <= 3;
    RFast rfobj(Y, nth, thp);
    FILE *fo = fopen((pre + ".txt").c_str(), "w"); if (!fo) { perror("out"); return 1; }
    const int DL = 41, DC = 13; std::vector<double> dcnt(DL * DC, 0), dsum(DL * DC, 0), dsq(DL * DC, 0); const double mu0 = 0.565723;
    Tree tr; Rng g; Snap T;
    for (long r = run0; r < run0 + runs; r++) {
        g.seed2(seed, (uint64_t)r); tr.init(); tr.model = model;
        for (size_t k = 0; k < cps.size(); k++) {
            tr.run_to(cps[k], g);
            tree_to_snap(tr, T); T.build_cred();
            FastGen fg(Y, T, model); if (useR) { fg.rf = &rfobj; rfobj.reset(); }
            MultiSink sa, sb; sa.n = sb.n = nth; for (int i = 0; i < nth; i++) { sa.th[i] = sb.th[i] = th[i].data(); sa.acc[i] = sb.acc[i] = 0; }
            double S = 0; size_t C = T.size();
            bool dg = diag && k + 1 == cps.size();
            for (uint32_t u = 0; u < C; u++) {
                S += fg.conv[u];
                if (dg) {
                    double a0 = sa.acc[0], b0 = sb.acc[0];
                    fg.compute(u, sa, sb);
                    double psi = (sa.acc[0] - a0) + th[0][p] * fg.conv[u] - mu0 * ((sb.acc[0] - b0) + th[0][p]);
                    int li = (int)std::min<uint64_t>(T.L[u], DL - 1), ci = (int)std::min<uint32_t>(T.c[u], DC - 1);
                    dcnt[li * DC + ci] += 1; dsum[li * DC + ci] += psi; dsq[li * DC + ci] += psi * psi;
                } else fg.compute(u, sa, sb);
            }
            fprintf(fo, "%ld %.6g %zu %.9f", r, cps[k], C, S);
            if (useR) for (int i = 0; i < nth; i++) { sa.acc[i] += rfobj.A[i]; sb.acc[i] += rfobj.B[i]; }
            for (int i = 0; i < nth; i++) fprintf(fo, " %.12e %.12e", sa.acc[i] + th[i][p] * S, sb.acc[i] + th[i][p] * (double)C);
            fprintf(fo, "\n");
            tree_to_snap(tr, T);
        }
    }
    fclose(fo);
    if (diag) {
        FILE *fd = fopen((pre + ".diag.txt").c_str(), "w");
        for (int l = 0; l < DL; l++) for (int c = 0; c < DC; c++) if (dcnt[l * DC + c] > 0) fprintf(fd, "%d %d %.0f %.6e %.6e\n", l, c, dcnt[l * DC + c], dsum[l * DC + c], dsq[l * DC + c]);
        fclose(fd);
    }
    return 0;
}
