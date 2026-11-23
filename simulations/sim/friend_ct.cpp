// friend_ct.cpp -- exact continuous-time simulator of the random friend tree, core only.
//
// Model: every vertex fires at rate 1 and gives a new leaf to a uniformly chosen neighbour.
// The jump chain is the Krapivsky-Redner isotropic-redirection tree (r = 1) / random friend tree.
//
// Only core vertices (degree >= 2) are stored. For a core vertex u: L[u] = number of leaf neighbours,
// c[u] = number of core neighbours. Leaves of u fire at rate 1 each and always hit u, so between
// events that involve u the count L[u] is a Yule process (pure birth, rate L). It is therefore updated
// lazily: L(t) = L(s) + NegBin(L(s), exp(-(t-s))) whenever u is touched. Core vertices fire at total
// rate C; the firing vertex u picks a uniform neighbour among D = L + c:
//   leaf (prob L/D):  the leaf becomes a core vertex with (L,c) = (1,1); u -> (L-1, c+1)   ["conversion"]
//   core neighbour w: w gains a leaf.
// This is exact in distribution (up to the documented normal approximation of Poisson(lambda) for
// lambda > 1e9, which only concerns hubs with enormous degree).
//
// Start: T_3 (path on three vertices): one core vertex with L = 2, c = 0, at t = 0.
//
// Usage: friend_ct SEED RUNS TMAX DT [HIST]
//   prints, per run and per checkpoint t = k*DT:  run t n C S F Cv A Z
//   n = number of vertices, C = core size, S = sum_u L_u/D_u, and for the window (t-DT, t]:
//   F = core firings, Cv = conversions, A = sum over firings of L/D just before the firing,
//   Z = number of core vertices with L = 0 at the checkpoint.
//   With HIST = 1, the last checkpoint of every run also prints "H L c count sumcredit" lines
//   (L, c capped at 30): local-state histogram with the summed credit of the vertices in the cell.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <vector>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next_u64(void) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}
static uint64_t splitmix64(uint64_t &x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void seed_rng(uint64_t seed) { for (int i = 0; i < 4; i++) s[i] = splitmix64(seed); }
static inline double unif_pos(void) { return ((next_u64() >> 11) + 0.5) * 0x1.0p-53; }   // in (0,1)
// unbiased integer in [0, n)
static inline uint64_t rand_below(uint64_t n) {
    // Lemire's method
    __uint128_t m = (__uint128_t)next_u64() * n;
    uint64_t l = (uint64_t)m;
    if (l < n) {
        uint64_t t = (-n) % n;
        while (l < t) { m = (__uint128_t)next_u64() * n; l = (uint64_t)m; }
    }
    return (uint64_t)(m >> 64);
}
static inline double std_normal(void) {
    // Box-Muller (one value)
    double u1 = unif_pos(), u2 = unif_pos();
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586476925 * u2);
}
static inline double std_exp(void) { return -log(unif_pos()); }
// Gamma(shape a >= 1, scale 1), Marsaglia-Tsang (exact rejection)
static double gamma_ge1(double a) {
    double d = a - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do { x = std_normal(); v = 1.0 + c * x; } while (v <= 0.0);
        v = v * v * v;
        double u = unif_pos();
        if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v;
    }
}
// Poisson(lam)
static uint64_t poisson(double lam) {
    if (lam < 10.0) {
        // multiplication method
        double el = exp(-lam), p = 1.0; uint64_t k = 0;
        for (;;) { p *= unif_pos(); if (p <= el) return k; k++; }
    }
    if (lam > 1e9) {
        double x = lam + sqrt(lam) * std_normal();
        return x < 0 ? 0 : (uint64_t)llround(x);
    }
    // PTRS, Hormann 1993 (exact)
    double slam = sqrt(lam), loglam = log(lam);
    double b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
    double invalpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2.0);
    for (;;) {
        double U = unif_pos() - 0.5, V = unif_pos();
        double us = 0.5 - fabs(U);
        double kf = floor((2.0 * a / us + b) * U + lam + 0.43);
        if (us >= 0.07 && V <= vr) return (uint64_t)kf;
        if (kf < 0 || (us < 0.013 && V > us)) continue;
        if (log(V) + log(invalpha) - log(a / (us * us) + b) <= -lam + kf * loglam - lgamma(kf + 1.0))
            return (uint64_t)kf;
    }
}
// Yule process: population after time tau, starting from l >= 1 individuals.
static inline uint64_t yule(uint64_t l, double tau) {
    if (tau <= 0.0) return l;
    double q = -expm1(-tau);          // 1 - p, p = exp(-tau)
    if (q <= 0.0) return l;
    if (l <= 12) {
        double lq = log(q);           // < 0
        uint64_t add = 0;
        for (uint64_t i = 0; i < l; i++) add += (uint64_t)floor(log(unif_pos()) / lq);
        return l + add;
    }
    double p = exp(-tau);
    double lam = gamma_ge1((double)l) * (q / p);
    return l + poisson(lam);
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s SEED RUNS TMAX DT [HIST]\n", argv[0]); return 1; }
    uint64_t seed = strtoull(argv[1], 0, 10);
    long runs = atol(argv[2]);
    double tmax = atof(argv[3]), dt = atof(argv[4]);
    int hist = argc > 5 ? atoi(argv[5]) : 0;
    seed_rng(seed);
    std::vector<uint64_t> L;
    std::vector<double> tl;
    std::vector<uint32_t> c, parent, first_child, next_sib;
    const uint32_t NONE = 0xffffffffu;
    const int HC = 31;
    std::vector<double> hcnt(HC * HC), hcred(HC * HC);
    printf("# run t n C S F Cv A Z\n");
    for (long run = 0; run < runs; run++) {
        L.assign(1, 2); tl.assign(1, 0.0); c.assign(1, 0); parent.assign(1, NONE);
        first_child.assign(1, NONE); next_sib.assign(1, NONE);
        double t = 0.0, tnext = dt;
        uint64_t F = 0, Cv = 0; double A = 0.0;
        for (;;) {
            uint64_t C = L.size();
            double tev = t + std_exp() / (double)C;
            while (tev >= tnext) {
                // checkpoint at time tnext: bring every vertex up to date
                uint64_t n = C, Z = 0; double S = 0.0;
                for (uint64_t u = 0; u < C; u++) {
                    if (L[u]) { L[u] = yule(L[u], tnext - tl[u]); S += (double)L[u] / (double)(L[u] + c[u]); }
                    else Z++;
                    tl[u] = tnext;
                    n += L[u];
                }
                printf("%ld %.6g %llu %llu %.6f %llu %llu %.6f %llu\n", run, tnext, (unsigned long long)n,
                       (unsigned long long)C, S, (unsigned long long)F, (unsigned long long)Cv, A,
                       (unsigned long long)Z);
                F = 0; Cv = 0; A = 0.0;
                if (tnext + dt > tmax + 1e-9) {
                    if (hist) {
                        for (auto &x : hcnt) x = 0; for (auto &x : hcred) x = 0;
                        for (uint64_t u = 0; u < C; u++) {
                            double cred = 0.0;
                            if (parent[u] != NONE) cred += 1.0 / (double)(L[parent[u]] + c[parent[u]]);
                            for (uint32_t w = first_child[u]; w != NONE; w = next_sib[w])
                                cred += 1.0 / (double)(L[w] + c[w]);
                            int a = L[u] > 30 ? 30 : (int)L[u], b = c[u] > 30 ? 30 : (int)c[u];
                            hcnt[a * HC + b] += 1; hcred[a * HC + b] += cred;
                        }
                        for (int a = 0; a < HC; a++) for (int b = 0; b < HC; b++)
                            if (hcnt[a * HC + b] > 0)
                                printf("H %d %d %.0f %.6f\n", a, b, hcnt[a * HC + b], hcred[a * HC + b]);
                    }
                    goto done;
                }
                tnext += dt;
            }
            t = tev;
            uint32_t u = (uint32_t)rand_below(C);
            if (L[u]) L[u] = yule(L[u], t - tl[u]);
            tl[u] = t;
            uint64_t D = L[u] + c[u];
            F++;
            A += (double)L[u] / (double)D;
            uint64_t j = rand_below(D);
            if (j < L[u]) {
                // conversion
                Cv++;
                L[u]--; c[u]++;
                uint32_t w = (uint32_t)L.size();
                L.push_back(1); tl.push_back(t); c.push_back(1); parent.push_back(u);
                first_child.push_back(NONE); next_sib.push_back(first_child[u]);
                first_child[u] = w;
            } else {
                uint64_t k = j - L[u];      // index among the core neighbours, uniform in [0, c[u])
                uint32_t w;
                if (parent[u] != NONE) {
                    if (k == 0) w = parent[u];
                    else { w = first_child[u]; for (uint64_t i = 1; i < k; i++) w = next_sib[w]; }
                } else {
                    w = first_child[u]; for (uint64_t i = 0; i < k; i++) w = next_sib[w];
                }
                if (L[w]) L[w] = yule(L[w], t - tl[w]);
                tl[w] = t;
                L[w]++;
            }
        }
    done:;
        fflush(stdout);
    }
    return 0;
}
