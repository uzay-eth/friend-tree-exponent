/* naive_chain.c -- independent direct discrete simulator of the random friend tree.
 *
 * Written cold for the simulation cross-check (Claude Opus 5.5); shares no code with sim/,
 * fixed-size/ or route-exact/.  Every vertex is stored explicitly.  Start: T_2, one edge.
 * Step: V uniform among the n vertices, W uniform among the neighbours of V, new leaf attached to W.
 * Neighbour list of x = (parent of x, if any) followed by its children, so a uniform neighbour is one
 * uniform index into that list.  RNG: SplitMix64 counter stream per run (run r uses seed ^ golden*r),
 * so results do not depend on the number of threads.
 *
 * Modes
 *   naive_chain fixed SEED RUNS NMAX [MODEL]
 *       mean and sd of C_n and S_n for every n = 2..NMAX (all runs reach NMAX).
 *   naive_chain sizes SEED RUNS n1 n2 ... nk [MODEL]
 *       mean and sd of C_n and S_n at the listed (increasing) sizes.
 *   naive_chain yule SEED RUNS T [MODEL]
 *       each run draws N = G1+G2+G3 with Gi iid geometric on {1,2,...}, success prob exp(-T)
 *       (exactly the law of the Yule size at time T from three vertices), grows the chain from
 *       T_2 to size N (T_3 is the deterministic path), and records C and S.  Prints mean/sd of
 *       N, C, S and writes one line "N C S" per run to stderr if the env var DUMP is set.
 *   naive_chain fixedN SEED RUNS T [MODEL]
 *       as yule but with the deterministic size N = round(3 exp(T)) (negative control).
 * MODEL: 0 = friend tree (default); 1 = attach to V itself (random recursive tree, negative control);
 *        2 = W chosen among neighbours of V proportionally to degree (negative control).
 * S = sum over core u of L_u / D_u (leaf neighbours over degree).
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline uint64_t sm64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline uint64_t below(uint64_t *st, uint64_t n) { /* unbiased, rejection */
    uint64_t lim = UINT64_MAX - UINT64_MAX % n;
    for (;;) { uint64_t r = sm64(st); if (r < lim) return r % n; }
}
static inline double unif01(uint64_t *st) { return ((sm64(st) >> 11) + 0.5) * 0x1.0p-53; }

typedef struct {
    int64_t cap, n;
    int64_t *parent, *deg;    /* deg = number of neighbours */
    int64_t *kid_off, *kid_cnt, *kid_cap; /* children stored in per-vertex blocks of an arena */
    int64_t *arena; int64_t arena_len, arena_cap;
    int64_t core;
} Tree;

static void tree_init(Tree *t, int64_t cap) {
    memset(t, 0, sizeof *t);
    t->cap = cap;
    t->parent = malloc(cap * sizeof(int64_t)); t->deg = malloc(cap * sizeof(int64_t));
    t->kid_off = malloc(cap * sizeof(int64_t)); t->kid_cnt = malloc(cap * sizeof(int64_t));
    t->kid_cap = malloc(cap * sizeof(int64_t));
    t->arena_cap = 4 * cap + 16; t->arena = malloc(t->arena_cap * sizeof(int64_t));
}
static void tree_grow_cap(Tree *t, int64_t need) {
    if (need <= t->cap) return;
    int64_t c = t->cap; while (c < need) c *= 2;
    t->parent = realloc(t->parent, c * sizeof(int64_t)); t->deg = realloc(t->deg, c * sizeof(int64_t));
    t->kid_off = realloc(t->kid_off, c * sizeof(int64_t)); t->kid_cnt = realloc(t->kid_cnt, c * sizeof(int64_t));
    t->kid_cap = realloc(t->kid_cap, c * sizeof(int64_t));
    t->cap = c;
}
static void add_kid(Tree *t, int64_t x, int64_t y) {
    if (t->kid_cnt[x] == t->kid_cap[x]) {           /* move block to a doubled block at arena end */
        int64_t nc = t->kid_cap[x] ? 2 * t->kid_cap[x] : 2;
        if (t->arena_len + nc > t->arena_cap) {
            while (t->arena_len + nc > t->arena_cap) t->arena_cap *= 2;
            t->arena = realloc(t->arena, t->arena_cap * sizeof(int64_t));
        }
        memcpy(t->arena + t->arena_len, t->arena + t->kid_off[x], t->kid_cnt[x] * sizeof(int64_t));
        t->kid_off[x] = t->arena_len; t->kid_cap[x] = nc; t->arena_len += nc;
    }
    t->arena[t->kid_off[x] + t->kid_cnt[x]++] = y;
}
static void tree_reset(Tree *t) {  /* T_2: vertex 0 -- vertex 1 */
    t->n = 2; t->arena_len = 0; t->core = 0;
    for (int i = 0; i < 2; i++) { t->kid_cnt[i] = 0; t->kid_cap[i] = 0; t->kid_off[i] = 0; }
    t->parent[0] = -1; t->deg[0] = 1;
    t->parent[1] = 0;  t->deg[1] = 1;
    add_kid(t, 0, 1);
}
static inline int64_t nbr(const Tree *t, int64_t x, int64_t j) {
    if (t->parent[x] >= 0) { if (j == 0) return t->parent[x]; j--; }
    return t->arena[t->kid_off[x] + j];
}
static void step(Tree *t, uint64_t *st, int model) {
    tree_grow_cap(t, t->n + 1);
    int64_t v = (int64_t)below(st, (uint64_t)t->n), w;
    if (model == 1) w = v;
    else if (model == 2) {           /* degree-biased neighbour */
        int64_t tot = 0;
        for (int64_t j = 0; j < t->deg[v]; j++) tot += t->deg[nbr(t, v, j)];
        int64_t r = (int64_t)below(st, (uint64_t)tot), j = 0;
        for (;; j++) { r -= t->deg[nbr(t, v, j)]; if (r < 0) break; }
        w = nbr(t, v, j);
    } else w = nbr(t, v, (int64_t)below(st, (uint64_t)t->deg[v]));
    int64_t y = t->n++;
    t->parent[y] = w; t->deg[y] = 1; t->kid_cnt[y] = 0; t->kid_cap[y] = 0; t->kid_off[y] = 0;
    add_kid(t, w, y);
    if (t->deg[w] == 1) t->core++;
    t->deg[w]++;
}
static double tree_S(const Tree *t, int64_t *L) {  /* sum over core of L_u/D_u */
    for (int64_t x = 0; x < t->n; x++) L[x] = 0;
    for (int64_t x = 0; x < t->n; x++) if (t->deg[x] == 1) L[nbr(t, x, 0)]++;
    double S = 0;
    for (int64_t x = 0; x < t->n; x++) if (t->deg[x] >= 2) S += (double)L[x] / (double)t->deg[x];
    return S;
}
static int64_t geom1(uint64_t *st, double p) { /* geometric on {1,2,...}, success prob p */
    if (p >= 1.0) return 1;
    return 1 + (int64_t)floor(log(unif01(st)) / log1p(-p));
}

typedef struct { double s, ss; } Acc;
static void acc(Acc *a, double x) { a->s += x; a->ss += x * x; }
static void pr(const char *name, Acc a, long R) {
    double m = a.s / R, var = a.ss / R - m * m; if (var < 0) var = 0;
    double sd = sqrt(var * R / (R - 1.0));
    printf(" %s_mean %.10g %s_sd %.10g %s_se %.6g", name, m, name, sd, name, sd / sqrt((double)R));
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "see header for usage\n"); return 2; }
    const char *mode = argv[1];
    uint64_t seed = strtoull(argv[2], 0, 10);
    long R = atol(argv[3]);
    int model = 0;
    int64_t *Lbuf = NULL; int64_t Lcap = 0;
    Tree t; tree_init(&t, 1 << 10);
    if (!strcmp(mode, "fixed")) {
        int nmax = atoi(argv[4]); if (argc > 5) model = atoi(argv[5]);
        Acc *aC = calloc(nmax + 1, sizeof(Acc)), *aS = calloc(nmax + 1, sizeof(Acc));
        Lbuf = malloc((nmax + 1) * sizeof(int64_t));
        for (long r = 0; r < R; r++) {
            uint64_t st = seed ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(r + 1));
            tree_reset(&t);
            for (;;) {
                acc(&aC[t.n], (double)t.core);
                acc(&aS[t.n], tree_S(&t, Lbuf));
                if (t.n == nmax) break;
                step(&t, &st, model);
            }
        }
        printf("# fixed seed=%llu runs=%ld nmax=%d model=%d\n", (unsigned long long)seed, R, nmax, model);
        for (int n = 2; n <= nmax; n++) { printf("n %d", n); pr("C", aC[n], R); pr("S", aS[n], R); printf("\n"); }
    } else if (!strcmp(mode, "sizes")) {
        int k = argc - 4; int64_t sz[64]; int nm = 0;
        for (int i = 0; i < k; i++) { sz[nm++] = atoll(argv[4 + i]); }
        /* optional trailing model flag written as m<digit> */
        if (nm && argv[argc - 1][0] == 'm') { model = argv[argc - 1][1] - '0'; nm--; }
        Acc aC[64] = {{0}}, aS[64] = {{0}};
        Lcap = sz[nm - 1] + 1; Lbuf = malloc(Lcap * sizeof(int64_t));
        for (long r = 0; r < R; r++) {
            uint64_t st = seed ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(r + 1));
            tree_reset(&t);
            for (int i = 0; i < nm; i++) {
                while (t.n < sz[i]) step(&t, &st, model);
                acc(&aC[i], (double)t.core); acc(&aS[i], tree_S(&t, Lbuf));
            }
        }
        printf("# sizes seed=%llu runs=%ld model=%d\n", (unsigned long long)seed, R, model);
        for (int i = 0; i < nm; i++) { printf("n %lld", (long long)sz[i]); pr("C", aC[i], R); pr("S", aS[i], R); printf("\n"); }
    } else if (!strcmp(mode, "yule") || !strcmp(mode, "fixedN")) {
        double T = atof(argv[4]); if (argc > 5) model = atoi(argv[5]);
        double p = exp(-T);
        int dump = getenv("DUMP") != NULL;
        Acc aN = {0}, aC = {0}, aS = {0};
        for (long r = 0; r < R; r++) {
            uint64_t st = seed ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(r + 1));
            int64_t N;
            if (mode[0] == 'y') N = geom1(&st, p) + geom1(&st, p) + geom1(&st, p);
            else N = llround(3.0 * exp(T));
            tree_reset(&t);
            while (t.n < N) step(&t, &st, model);
            if (Lcap < N + 1) { Lcap = 2 * N + 1; free(Lbuf); Lbuf = malloc(Lcap * sizeof(int64_t)); }
            double S = tree_S(&t, Lbuf);
            acc(&aN, (double)N); acc(&aC, (double)t.core); acc(&aS, S);
            if (dump) fprintf(stderr, "%lld %lld %.6f\n", (long long)N, (long long)t.core, S);
        }
        printf("# %s seed=%llu runs=%ld T=%g model=%d theory_meanN=%.6g\n", mode, (unsigned long long)seed, R, T, model, 3.0 * exp(T));
        printf("T %g", T); pr("N", aN, R); pr("C", aC, R); pr("S", aS, R); printf("\n");
    } else { fprintf(stderr, "unknown mode\n"); return 2; }
    return 0;
}
