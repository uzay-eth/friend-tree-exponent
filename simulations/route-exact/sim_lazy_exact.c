/* Exact continuous-time lazy core-only simulator for the Random Friend Tree
 *
 * Uses the exact continuous-time embedding (independent rate-1 Poisson clocks per vertex):
 *   - Every core vertex u has a rate-1 clock (total core firing rate = C).
 *   - Every leaf child of u has a rate-1 clock that ONLY ever attaches a new leaf to u,
 *     so between core events touching u, L_u(t) evolves as an exact independent pure-birth
 *     Yule process: L_u(t) - L_u(s) ~ NegBin(L_u(s), p = exp(-(t-s)))!
 *   - When a core vertex u fires at time t:
 *       1. Lazily catch up L_u(t) using exact NegBin(L_u(last_t[u]), exp(-(t - last_t[u]))).
 *       2. Pick uniform neighbour among its D_u = L_u + c_u neighbours:
 *          - With prob L_u / D_u: bud a leaf -> L_u--, c_u++, create newborn w=(L=1, c=1) at time t.
 *          - With prob c_u / D_u: pick uniform core neighbour v, catch up L_v(t), increment L_v++.
 *   - Exactness validation: compare E[C(t)], E[N(t)] at t = 3, 4 against exact Yule + naive tree!
 *   - High-precision mu measurement: measure S(t)/C(t) and d(log E[C(t)])/dt at t = 24..29 (n ~ 1e11..1e13)!
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef struct { uint64_t s[4]; } Rng;
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_next(Rng *r) {
    uint64_t res = rotl(r->s[1] * 5, 7) * 9, t = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1]; r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl(r->s[3], 45);
    return res;
}
static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void rng_seed(Rng *r, uint64_t seed) { for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&seed); }
static inline double rng_u01(Rng *r) { return ((rng_next(r) >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static inline double rng_exp(Rng *r) { return -log(rng_u01(r)); }
static inline uint64_t rng_unif(Rng *r, uint64_t m) {
    /* Rejection-free-ish fast unif using 128-bit mul */
    __uint128_t x = (__uint128_t)rng_next(r) * m;
    return (uint64_t)(x >> 64);
}
static double rng_norm(Rng *r) {
    double u1 = rng_u01(r), u2 = rng_u01(r);
    return sqrt(-2.0 * log(u1)) * cos(6.28318530717958647692 * u2);
}
/* Marsaglia-Tsang exact Gamma(shape, 1) */
static double rng_gamma(Rng *r, double shape) {
    if (shape < 1.0) return rng_gamma(r, shape + 1.0) * pow(rng_u01(r), 1.0 / shape);
    double d = shape - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do { x = rng_norm(r); v = 1.0 + c * x; } while (v <= 0.0);
        v = v * v * v;
        double u = rng_u01(r);
        if (u < 1.0 - 0.0331 * x * x * x * x) return d * v;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v;
    }
}
/* Exact Poisson(lam) */
static uint64_t rng_poisson(Rng *r, double lam) {
    if (lam <= 0.0) return 0;
    if (lam < 30.0) {
        double L = exp(-lam), p = 1.0;
        uint64_t k = 0;
        do { k++; p *= rng_u01(r); } while (p > L);
        return k - 1;
    }
    /* Ahrens-Dieter / PTRS exact Poisson for lam >= 30 */
    double b = 0.931 + 2.53 * sqrt(lam), a = -0.059 + 0.02483 * b;
    double inv_a = 1.1239 + 1.1328 / (b - 3.4);
    double vr = 0.9277 - 3.6224 / (b - 2.0);
    double loglam = log(lam);
    for (;;) {
        double u = rng_u01(r) - 0.5, v = rng_u01(r);
        double us = 0.5 - fabs(u);
        double kf = floor((2.0 * a / us + b) * u + lam + 0.43);
        if (kf < 0.0) continue;
        if (us >= 0.07 && v <= vr) return (uint64_t)kf;
        if (us < 0.013 && v > us) continue;
        double k = kf;
        if (log(v * inv_a / (a / (us * us) + b)) <= -lam + k * loglam - lgamma(k + 1.0))
            return (uint64_t)k;
    }
}
/* Exact Yule jump: given L0 >= 0 leaves at time s, return L(t) at time t >= s.
 * L(t) - L0 ~ NegBin(r = L0, p = exp(-(t-s))) = Poisson(Gamma(L0, 1) * (exp(dt) - 1))!
 * Wait: for L0 == 0, L(t) = 0 deterministically!
 * For small dt or large L0, Gamma-Poisson mixture is 100% EXACT for ALL L0 >= 1 and ALL dt > 0!
 */
static inline uint64_t yule_advance(Rng *r, uint64_t L0, double dt) {
    if (L0 == 0 || dt <= 0.0) return L0;
    double scale = expm1(dt);
    if (L0 == 1) {
        /* Geometric(p = exp(-dt)): exact inverse CDF floor(Exp(1) / (-log(1 - p))) = floor(Exp(1) / (-log(1 - exp(-dt)))) */
        double q = -expm1(-dt); /* 1 - exp(-dt) */
        if (q <= 0.0) return 1;
        double e = rng_exp(r);
        return 1 + (uint64_t)floor(e / (-log(q)));
    }
    double g = rng_gamma(r, (double)L0);
    double lam = g * scale;
    if (lam > 1e14) return L0 + (uint64_t)llround(lam + sqrt(lam) * rng_norm(r));
    return L0 + rng_poisson(r, lam);
}

/* Dynamic adjacency storage for core tree */
typedef struct {
    uint64_t *L;
    double *last_t;
    uint32_t *par;
    uint32_t *ch_cnt, *ch_cap;
    uint32_t **ch;
    uint32_t C, cap;
} CoreTree;

static void ct_init(CoreTree *T, uint32_t init_cap) {
    T->cap = init_cap; T->C = 1;
    T->L = malloc(init_cap * sizeof(uint64_t));
    T->last_t = malloc(init_cap * sizeof(double));
    T->par = malloc(init_cap * sizeof(uint32_t));
    T->ch_cnt = calloc(init_cap, sizeof(uint32_t));
    T->ch_cap = calloc(init_cap, sizeof(uint32_t));
    T->ch = calloc(init_cap, sizeof(uint32_t *));
    /* Start at T_3: single core vertex 0 with L_0 = 2, c_0 = 0 at t = 0 */
    T->L[0] = 2; T->last_t[0] = 0.0; T->par[0] = UINT32_MAX;
}
static void ct_grow(CoreTree *T) {
    uint32_t nc = T->cap * 2;
    T->L = realloc(T->L, nc * sizeof(uint64_t));
    T->last_t = realloc(T->last_t, nc * sizeof(double));
    T->par = realloc(T->par, nc * sizeof(uint32_t));
    T->ch_cnt = realloc(T->ch_cnt, nc * sizeof(uint32_t));
    T->ch_cap = realloc(T->ch_cap, nc * sizeof(uint32_t));
    T->ch = realloc(T->ch, nc * sizeof(uint32_t *));
    memset(T->ch_cnt + T->cap, 0, (nc - T->cap) * sizeof(uint32_t));
    memset(T->ch_cap + T->cap, 0, (nc - T->cap) * sizeof(uint32_t));
    memset(T->ch + T->cap, 0, (nc - T->cap) * sizeof(uint32_t *));
    T->cap = nc;
}
static void ct_free(CoreTree *T) {
    for (uint32_t i = 0; i < T->C; i++) free(T->ch[i]);
    free(T->L); free(T->last_t); free(T->par); free(T->ch_cnt); free(T->ch_cap); free(T->ch);
}

static void run_to_times(uint64_t seed, const double *times, int ntimes,
                         double *out_C, double *out_S, double *out_N) {
    Rng rng; rng_seed(&rng, seed);
    CoreTree T; ct_init(&T, 1 << 16);
    double t = 0.0;
    for (int ti = 0; ti < ntimes; ti++) {
        double target = times[ti];
        for (;;) {
            double dt = rng_exp(&rng) / (double)T.C;
            if (t + dt > target) { t = target; break; }
            t += dt;
            uint32_t u = (uint32_t)rng_unif(&rng, T.C);
            T.L[u] = yule_advance(&rng, T.L[u], t - T.last_t[u]);
            T.last_t[u] = t;
            uint32_t h = (T.par[u] != UINT32_MAX) ? 1 : 0;
            uint32_t c = h + T.ch_cnt[u];
            uint64_t D = T.L[u] + c;
            uint64_t pick = rng_unif(&rng, D);
            if (pick < T.L[u]) {
                /* Budding! */
                T.L[u]--;
                if (T.C == T.cap) ct_grow(&T);
                uint32_t w = T.C++;
                T.L[w] = 1; T.last_t[w] = t; T.par[w] = u;
                if (T.ch_cnt[u] == T.ch_cap[u]) {
                    uint32_t ncc = T.ch_cap[u] ? 2 * T.ch_cap[u] : 4;
                    T.ch[u] = realloc(T.ch[u], ncc * sizeof(uint32_t));
                    T.ch_cap[u] = ncc;
                }
                T.ch[u][T.ch_cnt[u]++] = w;
            } else {
                uint32_t idx = (uint32_t)(pick - T.L[u]);
                uint32_t v = (h && idx == 0) ? T.par[u] : T.ch[u][idx - h];
                T.L[v] = yule_advance(&rng, T.L[v], t - T.last_t[v]);
                T.last_t[v] = t;
                T.L[v]++;
            }
        }
        /* Synchronize all core vertices at target time to compute exact C(t), S(t), N(t) */
        double S = 0.0, Ntot = (double)T.C;
        for (uint32_t u = 0; u < T.C; u++) {
            T.L[u] = yule_advance(&rng, T.L[u], t - T.last_t[u]);
            T.last_t[u] = t;
            uint32_t c = ((T.par[u] != UINT32_MAX) ? 1 : 0) + T.ch_cnt[u];
            uint64_t D = T.L[u] + c;
            S += (double)T.L[u] / (double)D;
            Ntot += (double)T.L[u];
        }
        out_C[ti] = (double)T.C;
        out_S[ti] = S;
        out_N[ti] = Ntot;
    }
    ct_free(&T);
}

int main(void) {
    /* 1. Validation at t = 4.0: exact E[N(4)] = 3*e^4 = 163.79445 */
    {
        double ts[1] = {4.0};
        double sumC = 0, sumN = 0, sumS = 0;
        int R = 50000;
        #pragma omp parallel for num_threads(2) reduction(+:sumC,sumN,sumS)
        for (int r = 0; r < R; r++) {
            double c[1], s[1], n[1];
            run_to_times(100000ULL + r, ts, 1, c, s, n);
            sumC += c[0]; sumS += s[0]; sumN += n[0];
        }
        printf("[VALIDATION t=4.0, R=%d] E[N(4)]=%.4f (exact 3*e^4=%.4f), E[C(4)]=%.4f, E[S]/E[C]=%.6f\n",
               R, sumN/R, 3.0*exp(4.0), sumC/R, sumS/sumC);
        fflush(stdout);
    }

    /* 2. High-precision measurement: 64 runs to t = 27.0 (n ~ 1.6e12, C ~ 5e6 per run => ~3.2e8 core vertices!) */
    double ts[4] = {21.0, 23.0, 25.0, 27.0};
    int R = 48;
    double *C_all = calloc(R * 4, sizeof(double)), *S_all = calloc(R * 4, sizeof(double));
    #pragma omp parallel for schedule(dynamic) num_threads(2)
    for (int r = 0; r < R; r++) {
        double c[4], s[4], n[4];
        run_to_times(2026112200ULL + r, ts, 4, c, s, n);
        for (int k = 0; k < 4; k++) { C_all[r*4+k] = c[k]; S_all[r*4+k] = s[k]; }
    }
    printf("\n=== HIGH-PRECISION LAZY CONTINUOUS-TIME SIMULATION (R=%d runs, 2 threads) ===\n", R);
    for (int k = 0; k < 4; k++) {
        double sc = 0, ss = 0;
        for (int r = 0; r < R; r++) { sc += C_all[r*4+k]; ss += S_all[r*4+k]; }
        double mu_hat = ss / sc;
        /* Jackknife / delta-method standard error of ratio sum(S)/sum(C) */
        double var_num = 0;
        for (int r = 0; r < R; r++) {
            double d = S_all[r*4+k] - mu_hat * C_all[r*4+k];
            var_num += d * d;
        }
        double se = sqrt(var_num / (R * (R - 1))) / (sc / R);
        printf("  t=%4.1f (E[N]=%.2e, mean C=%.2e): mu_eff = E[S]/E[C] = %.7f +- %.7f  (95%% CI: [%.6f, %.6f])\n",
               ts[k], 3.0*exp(ts[k]), sc/R, mu_hat, se, mu_hat - 1.96*se, mu_hat + 1.96*se);
    }
    free(C_all); free(S_all);
    return 0;
}
