// residual_spec.cpp -- the next-generation kernel K_s of the renewal characterisation (see residual_sim.cpp for the
// objects and the simulation) evaluated on a grid of COMPLEX s from one set of simulated residual families.
// Purpose: look for zeros of det(I - K_s) other than s = mu in a strip Re s > sigma_0, i.e. for modes of the mean
// dynamics that would decay only slowly relative to the leading one.  Not a proof of anything.
// The simulation is the one of residual_sim.cpp (variant 0, or the capped 19-label variant with RES_VARIANT=4); only
// the bookkeeping differs: every reward is weighted with exp(-s t) for all grid points s = sigma_j + i omega_k.
// Output (binary doubles): header (T, batches, families, nr, NCOL, nsig, nom), sigma values, omega values, row marks,
// then per batch K, Klast, Kprev (each ns*nr*NCOL complex numbers) and Km (ns*NCOL complex numbers).
// Usage: residual_spec SEED BATCHES FAMILIES_PER_BATCH T OUT "sig1,sig2,..." "om1,om2,..."
#include "ftsim.hpp"
#include <cstring>
#include <string>

struct Fr {
    std::vector<uint64_t> L; std::vector<double> tl, mark; std::vector<int> par, c; std::vector<std::vector<int>> ch;
    int add(uint64_t l, double t, double m, int p, int cc) { L.push_back(l); tl.push_back(t); mark.push_back(m); par.push_back(p); c.push_back(cc); ch.emplace_back(); return (int)L.size() - 1; }
    size_t n() const { return L.size(); }
};
struct Particle { Fr f; double t0, a, tau; bool orig; };
static const int JX = 40;
static const int LB[5] = {41, 61, 101, 201, 501};
static inline int colof(uint64_t D) { if (D <= (uint64_t)JX) return (int)D - 2; for (int i = 4; i >= 0; i--) if (D >= (uint64_t)LB[i]) return JX - 1 + i; return 0; }
static const int NCOL = JX - 1 + 5;
static std::vector<double> parse(const char *s) { std::vector<double> v; const char *p = s; while (*p) { char *e; v.push_back(strtod(p, &e)); p = (*e == ',') ? e + 1 : e; if (e == p && *e) break; } return v; }

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: %s SEED BATCHES FAMILIES T OUT SIGMAS OMEGAS\n", argv[0]); return 1; }
    uint64_t seed = strtoull(argv[1], 0, 10); int nb = atoi(argv[2]); long nf = atol(argv[3]); double T = atof(argv[4]);
    std::vector<double> sig = parse(argv[6]), om = parse(argv[7]); int nsig = (int)sig.size(), nom = (int)om.size(), ns = nsig * nom;
    int variant = getenv("RES_VARIANT") ? atoi(getenv("RES_VARIANT")) : 0;
    std::vector<double> ra; for (int j = 2; j <= JX; j++) ra.push_back(1.0 / j);
    int ends[] = {41, 60, 61, 100, 101, 200, 201, 500, 501}; for (int e : ends) ra.push_back(1.0 / e); ra.push_back(0.0);
    int nr = (int)ra.size();
    FILE *fo = fopen(argv[5], "wb");
    double hdr[7] = {T, (double)nb, (double)nf, (double)nr, (double)NCOL, (double)nsig, (double)nom}; fwrite(hdr, 8, 7, fo);
    fwrite(sig.data(), 8, nsig, fo); fwrite(om.data(), 8, nom, fo); fwrite(ra.data(), 8, nr, fo);
    size_t blk = (size_t)ns * nr * NCOL * 2;
    std::vector<double> K(blk), Klast(blk), Kprev(blk), Km((size_t)ns * NCOL * 2), ew(nr), es(nsig), co(nom), si_(nom);
    Rng g; std::vector<Particle> stack; long long totpart = 0, totev = 0, totrew = 0;
    for (int b = 0; b < nb; b++) {
        std::fill(K.begin(), K.end(), 0.0); std::fill(Klast.begin(), Klast.end(), 0.0); std::fill(Kprev.begin(), Kprev.end(), 0.0); std::fill(Km.begin(), Km.end(), 0.0);
        for (long fam = 0; fam < nf; fam++) {
            g.seed2(seed, (uint64_t)b * 1000000007ULL + (uint64_t)fam);
            stack.clear(); { Particle p; p.f.add(1, 0.0, 0.0, -1, 1); p.t0 = 0; p.a = 0; p.tau = 0; p.orig = true; stack.push_back(std::move(p)); }
            while (!stack.empty()) {
                Particle P = std::move(stack.back()); stack.pop_back(); totpart++;
                Fr &f = P.f; double t = P.t0; double A = 0; for (int y : f.ch[0]) A += f.mark[y];
                for (;;) {
                    double kill = P.a;
                    double rate = (double)f.n() + kill + A;
                    t += g.expo() / rate; if (t >= T) break;
                    double u = g.unif() * rate; totev++;
                    if (u < kill) break;
                    u -= kill;
                    if (u < A) {
                        int y = -1; for (int cand : f.ch[0]) { if (u < f.mark[cand]) { y = cand; break; } u -= f.mark[cand]; }
                        if (y < 0) y = f.ch[0].back();
                        if (f.L[0]) f.L[0] = g.yule(f.L[0], t - f.tl[0]); f.tl[0] = t;
                        double p = 1.0 / (double)(f.L[0] + f.c[0]);
                        Particle Q; Q.t0 = t; Q.a = f.mark[y]; Q.orig = false; Q.tau = P.orig ? t : P.tau;
                        std::vector<int> map(f.n(), -1), order; order.push_back(y);
                        for (size_t i = 0; i < order.size(); i++) for (int w : f.ch[order[i]]) order.push_back(w);
                        for (int x : order) {
                            if (f.L[x]) f.L[x] = g.yule(f.L[x], t - f.tl[x]); f.tl[x] = t;
                            map[x] = Q.f.add(f.L[x], t, x == y ? Q.a : f.mark[x], x == y ? -1 : map[f.par[x]], f.c[x]);
                            if (x != y) Q.f.ch[map[f.par[x]]].push_back(map[x]);
                        }
                        if (g.unif() * f.mark[y] < p) Q.f.L[0] += 1;
                        stack.push_back(std::move(Q));
                        continue;
                    }
                    int x = (int)g.below(f.n());
                    if (f.L[x]) f.L[x] = g.yule(f.L[x], t - f.tl[x]); f.tl[x] = t;
                    uint64_t D = f.L[x] + f.c[x];
                    if (x == 0 && f.L[0] > 0) {                           // reward
                        double w = (double)f.L[0] / (double)D; int col = colof(variant == 4 && D > 20 ? 20 : D); double tw = P.orig ? t : P.tau; totrew++;
                        for (int r = 0; r < nr; r++) ew[r] = exp(-ra[r] * tw);
                        for (int j = 0; j < nsig; j++) es[j] = w * exp(-sig[j] * t);
                        for (int k = 0; k < nom; k++) { co[k] = cos(om[k] * t); si_[k] = -sin(om[k] * t); }
                        int last = t >= T - 1.0 ? 2 : (t >= T - 2.0 ? 1 : 0); double invD = 1.0 / (double)D;
                        for (int j = 0; j < nsig; j++) for (int k = 0; k < nom; k++) {
                            int s = j * nom + k; double zr = es[j] * co[k], zi = es[j] * si_[k];
                            double *row = &K[(((size_t)s * nr) * NCOL + col) * 2];
                            for (int r = 0; r < nr; r++) { row[(size_t)r * NCOL * 2] += zr * ew[r]; row[(size_t)r * NCOL * 2 + 1] += zi * ew[r]; }
                            if (last) { double *row2 = last == 2 ? &Klast[(((size_t)s * nr) * NCOL + col) * 2] : &Kprev[(((size_t)s * nr) * NCOL + col) * 2];
                                for (int r = 0; r < nr; r++) { row2[(size_t)r * NCOL * 2] += zr * ew[r]; row2[(size_t)r * NCOL * 2 + 1] += zi * ew[r]; } }
                            Km[((size_t)s * NCOL + col) * 2] += zr * invD; Km[((size_t)s * NCOL + col) * 2 + 1] += zi * invD;
                        }
                    }
                    uint64_t j = g.below(D);
                    if (j < f.L[x]) {
                        f.L[x]--; f.c[x]++;
                        int y = f.add(1, t, 1.0 / (double)(variant == 4 && D > 20 ? 20 : D), x, 1); f.ch[x].push_back(y);
                        if (x == 0) A += f.mark[y];
                    } else {
                        uint64_t k = j - f.L[x]; int w;
                        if (k == 0) { if (x == 0) continue; w = f.par[x]; }
                        else w = f.ch[x][k - 1];
                        if (f.L[w]) f.L[w] = g.yule(f.L[w], t - f.tl[w]); f.tl[w] = t; f.L[w]++;
                    }
                }
            }
        }
        double inv = 1.0 / (double)nf;
        for (auto *v : {&K, &Klast, &Kprev, &Km}) { for (double &x : *v) x *= inv; fwrite(v->data(), 8, v->size(), fo); }
        fflush(fo);
    }
    fprintf(stderr, "particles %lld events %lld rewards %lld\n", totpart, totev, totrew);
    fclose(fo); return 0;
}
