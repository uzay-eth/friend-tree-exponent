/*
 * Exact Lazy-Yule Core-Only Continuous-Time Friend-Tree Simulator (Validation)
 *
 * Exact continuous-time dynamics from T_3 (one core root with L=2 leaves):
 * - Every core vertex u has an independent rate-1 clock for all non-pure-Yule events:
 *     with prob L_u / D_u: promote a leaf (L_u--, attach new core child with L=1, k=0, h=1)
 *     with prob 1 / D_u for each core neighbor v: send a leaf to v (L_v++)
 * - Between events touching u, u's leaf cloud L_u evolves as an exact pure-birth Yule process:
 *     L_u(t) = L_u(tau) + NegBin(L_u(tau), exp(-(t - tau)))
 *   sampled in O(1) time via Gamma(L_u, exp(dt)-1) + Poisson(G)!
 *
 * Validation against exact small-t expectations (computed from exact_rationals_n21.txt):
 *   E_{T_3}[N(t)] = 3 e^t
 *   E_{T_3}[C(0.5)], E_{T_3}[C(1.0)], E_{T_3}[S(1.0)]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

typedef struct { uint64_t s[4]; } RNG;
static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_u64(RNG *r) {
    uint64_t res = rotl64(r->s[1] * 5, 7) * 9, t = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1]; r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl64(r->s[3], 45);
    return res;
}
static inline void rng_seed(RNG *r, uint64_t seed) {
    for (int i = 0; i < 4; i++) {
        uint64_t z = (seed += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        r->s[i] = z ^ (z >> 31);
    }
}
static inline double rng_unif(RNG *r) { return (rng_u64(r) >> 11) * 0x1.0p-53; }
static inline double rng_exp(RNG *r) { double u; do { u = rng_unif(r); } while (u <= 0.0); return -log(u); }
static inline uint32_t rng_int(RNG *r, uint32_t m) {
    uint64_t x = rng_u64(r) >> 32, prod = x * (uint64_t)m;
    uint32_t low = (uint32_t)prod;
    if (low < m) { uint32_t thr = (uint32_t)(-m) % m; while (low < thr) { x = rng_u64(r) >> 32; prod = x * (uint64_t)m; low = (uint32_t)prod; } }
    return (uint32_t)(prod >> 32);
}
static inline double rng_norm(RNG *r) {
    /* Polar Box-Muller */
    double u, v, s;
    do { u = 2.0 * rng_unif(r) - 1.0; v = 2.0 * rng_unif(r) - 1.0; s = u*u + v*v; } while (s >= 1.0 || s == 0.0);
    return u * sqrt(-2.0 * log(s) / s);
}

/* Exact Marsaglia-Tsang Gamma(shape = a >= 1, scale = 1) */
static inline double rng_gamma(RNG *r, double a) {
    if (a == 1.0) return rng_exp(r);
    double d = a - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do { x = rng_norm(r); v = 1.0 + c * x; } while (v <= 0.0);
        v = v * v * v;
        double u = rng_unif(r);
        if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v;
    }
}

/* Exact Poisson(lam) for any lam >= 0: Knuth for lam < 30, Hormann PTRS for lam >= 30 */
static inline uint64_t rng_poisson(RNG *r, double lam) {
    if (lam <= 0.0) return 0;
    if (lam < 30.0) {
        double L = exp(-lam), p = 1.0;
        uint64_t k = 0;
        do { k++; p *= rng_unif(r); } while (p > L);
        return k - 1;
    }
    /* Hormann (1993) PTRS exact transformed rejection */
    double b = 0.931 + 2.53 * sqrt(lam);
    double a = -0.059 + 0.02483 * b;
    double inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
    double v_r = 0.9277 - 3.6224 / (b - 2.0);
    double log_lam = log(lam);
    for (;;) {
        double U = rng_unif(r) - 0.5;
        double V = rng_unif(r);
        double us = 0.5 - fabs(U);
        double kf = floor((2.0 * a / us + b) * U + lam + 0.43);
        if (kf < 0.0) continue;
        if (us >= 0.07 && V <= v_r) return (uint64_t)kf;
        if (us < 0.013 && V > us) continue;
        double V2 = V * inv_alpha / (a / (us * us) + b);
        if (log(V2) <= kf * log_lam - lam - lgamma(kf + 1.0)) return (uint64_t)kf;
    }
}

/* Exact NegBin(L, exp(-dt)): number of NEW leaves born in a pure Yule cloud of size L over dt */
static inline uint64_t yule_add(RNG *r, uint64_t L, double dt) {
    if (L == 0 || dt <= 0.0) return 0;
    double scale = expm1(dt);
    double G = rng_gamma(r, (double)L) * scale;
    return rng_poisson(r, G);
}

typedef struct {
    uint64_t *L;
    double *tau;
    int32_t *par;
    int32_t *head, *next, *kcnt; /* linked list of core children! O(1) alloc, zero realloc per node! */
    int cap, C;
} CoreTree;

static void ct_init(CoreTree *T, int cap) {
    T->cap = cap;
    T->L = malloc(cap * sizeof(uint64_t));
    T->tau = malloc(cap * sizeof(double));
    T->par = malloc(cap * sizeof(int32_t));
    T->head = malloc(cap * sizeof(int32_t));
    T->next = malloc(cap * sizeof(int32_t));
    T->kcnt = malloc(cap * sizeof(int32_t));
}
static void ct_grow(CoreTree *T) {
    int nc = T->cap * 2;
    T->L = realloc(T->L, nc * sizeof(uint64_t));
    T->tau = realloc(T->tau, nc * sizeof(double));
    T->par = realloc(T->par, nc * sizeof(int32_t));
    T->head = realloc(T->head, nc * sizeof(int32_t));
    T->next = realloc(T->next, nc * sizeof(int32_t));
    T->kcnt = realloc(T->kcnt, nc * sizeof(int32_t));
    T->cap = nc;
}
static void ct_reset(CoreTree *T) {
    T->C = 1;
    T->L[0] = 2; T->tau[0] = 0.0; T->par[0] = -1; T->head[0] = -1; T->next[0] = -1; T->kcnt[0] = 0;
}
static inline void ct_sync(CoreTree *T, RNG *r, int u, double t) {
    double dt = t - T->tau[u];
    if (dt > 0.0 && T->L[u] > 0) T->L[u] += yule_add(r, T->L[u], dt);
    T->tau[u] = t;
}

static void ct_advance(CoreTree *T, RNG *r, double *t_cur, double t_end) {
    double t = *t_cur;
    for (;;) {
        double dt = rng_exp(r) / (double)T->C;
        if (t + dt >= t_end) {
            t = t_end;
            break;
        }
        t += dt;
        int u = (int)rng_int(r, (uint32_t)T->C);
        ct_sync(T, r, u, t);
        uint64_t Lu = T->L[u];
        uint32_t hu = (T->par[u] >= 0) ? 1 : 0;
        uint32_t cu = hu + (uint32_t)T->kcnt[u];
        uint64_t Du = Lu + cu;
        /* Pick uniform neighbor among Du: [0..Lu-1] = leaf (promote), [Lu..Du-1] = core neighbor */
        double pick = rng_unif(r) * (double)Du;
        if (pick < (double)Lu) {
            /* Promote a leaf! */
            T->L[u]--;
            if (T->C == T->cap) ct_grow(T);
            int v = T->C++;
            T->L[v] = 1; T->tau[v] = t; T->par[v] = u;
            T->head[v] = -1; T->next[v] = T->head[u]; T->head[u] = v;
            T->kcnt[v] = 0; T->kcnt[u]++;
        } else {
            uint32_t ci = rng_int(r, cu);
            int nbr;
            if (hu && ci == 0) nbr = T->par[u];
            else {
                uint32_t idx = hu ? (ci - 1) : ci;
                nbr = T->head[u];
                for (uint32_t s = 0; s < idx; s++) nbr = T->next[nbr];
            }
            ct_sync(T, r, nbr, t);
            T->L[nbr]++;
        }
    }
    for (int u = 0; u < T->C; u++) ct_sync(T, r, u, t_end);
    *t_cur = t_end;
}

int main(void) {
    /* Validation test 1: 2,000,000 runs at t = 0.5 and t = 1.0 */
    int NRUN = 2000000;
    double sumC05=0, sumN05=0, sumC10=0, sumN10=0, sumS10=0;
    #pragma omp parallel num_threads(4) reduction(+:sumC05,sumN05,sumC10,sumN10,sumS10)
    {
        int tid = omp_get_thread_num();
        RNG r; rng_seed(&r, 123456789ULL + 1000003ULL * tid);
        CoreTree T; ct_init(&T, 4096);
        #pragma omp for schedule(static)
        for (int run = 0; run < NRUN; run++) {
            ct_reset(&T);
            double t = 0.0;
            ct_advance(&T, &r, &t, 0.5);
            uint64_t N05 = T.C; for (int u = 0; u < T.C; u++) N05 += T.L[u];
            sumC05 += T.C; sumN05 += (double)N05;
            ct_advance(&T, &r, &t, 1.0);
            uint64_t N10 = T.C; double S10 = 0.0;
            for (int u = 0; u < T.C; u++) {
                N10 += T.L[u];
                double Du = (double)T.L[u] + (T.par[u] >= 0 ? 1 : 0) + T.kcnt[u];
                S10 += (double)T.L[u] / Du;
            }
            sumC10 += T.C; sumN10 += (double)N10; sumS10 += S10;
        }
        free(T.L); free(T.tau); free(T.par); free(T.head); free(T.next); free(T.kcnt);
    }
    printf("t=0.5: E[N]=%.6f (exact 3e^0.5=%.6f)  E[C]=%.6f\n", sumN05/NRUN, 3*exp(0.5), sumC05/NRUN);
    printf("t=1.0: E[N]=%.6f (exact 3e^1.0=%.6f)  E[C]=%.6f  E[S]=%.6f\n", sumN10/NRUN, 3*exp(1.0), sumC10/NRUN, sumS10/NRUN);

    /* Speed benchmark: 10 trees to t = 20 (N ~ 1.45 billion!) */
    double t0 = omp_get_wtime();
    RNG r; rng_seed(&r, 99999ULL);
    CoreTree T; ct_init(&T, 131072);
    double meanC20 = 0, meanN20 = 0;
    for (int i = 0; i < 10; i++) {
        ct_reset(&T); double t = 0.0;
        ct_advance(&T, &r, &t, 20.0);
        double N20 = T.C; for (int u = 0; u < T.C; u++) N20 += (double)T.L[u];
        meanC20 += T.C; meanN20 += N20;
    }
    double t1 = omp_get_wtime();
    printf("10 trees to t=20: mean C=%.1f, mean N=%.3e (exact 3e^20=%.3e), wall=%.3fs (%.1f ms/tree)\n",
           meanC20/10, meanN20/10, 3*exp(20.0), t1-t0, (t1-t0)*100.0);
    return 0;
}
