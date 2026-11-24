/* Experimental Motif Solver v5: Pure 1-Knot 3-Scale 4-Vertex Potential (Single-Threaded, 1 CPU)
 *
 * ZERO multi-knot splines! ZERO left-ending corner checks! ZERO alpha_t!
 * Every level (L2, L3, L4) uses the EXACT SAME 1-knot tail:
 *   mu^(L)_v = (1 - t_L) delta_INF + t_L delta_{T1(kap_L)},   t_L = D0_L / D_v.
 *
 *   Level 1 (1V): f(L, k) exact for D < DF (=1400), exact analytic q-tail for D >= DF
 *   Level 2 (2V): g(p -> w) with D0_g=72, M_g=1, KT_g=24
 *   Level 3 (3V): hV, hC with D0_3=28, M_3=1, KT_3=12, DHU=24, KHU=12, DCU=22, KCU=10
 *   Level 4 (4V): hP4, hF4, hY4, hT4 (Trident) with D0_4=14, M_4=1, KT_4=8
 *
 * CLI: exp_v5 <mode:1=UP|2=LO> <lam_init> [D0_g=72] [D0_3=28] [DHU=24] [DCU=22] [D0_4=14]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define NKNOTS 15
#define MAX_D0 100

static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next64(void) {
    uint64_t r = rotl(xs[1] * 5, 7) * 9, t = xs[1] << 17;
    xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3];
    xs[2] ^= t; xs[3] = rotl(xs[3], 45);
    return r;
}
static uint64_t sm(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline uint32_t unif(uint32_t m) {
    uint64_t x = next64() >> 32, prod = x * (uint64_t)m;
    uint32_t low = (uint32_t)prod;
    if (low < m) {
        uint32_t thr = (uint32_t)(-m) % m;
        while (low < thr) { x = next64() >> 32; prod = x * (uint64_t)m; low = (uint32_t)prod; }
    }
    return (uint32_t)(prod >> 32);
}

typedef struct {
    int D0, KT;
    int INNER_N, INF_ID, TAIL_BASE, NG;
    int lut[MAX_D0][MAX_D0]; /* [k][L] -> inner_id */
    int *gL, *gK, *gD, *gKind; /* 0=inner, 1=INF, 2=T1 */
    double *invD;
    int *pars, npars, *p_to_pi;
    int *chs,  nchs,  *w_to_wi;
    size_t pair_n, tri_n;
    size_t *tri_off1, *tri_off2;
    int max_rank;
    int *p_by_r[128], p_rcnt[128];
    int *c_by_r[128], c_rcnt[128];
} Level;

static inline int lvl_inner(const Level *S, int L, int k) { return S->lut[k][L]; }
static inline int lvl_tail(const Level *S, int k) {
    int kap = (k > S->KT) ? (S->KT + 1) : k;
    return S->TAIL_BASE + kap;
}
static inline size_t lvl_pair(const Level *S, int w1, int w2) {
    int i1 = S->w_to_wi[w1], i2 = S->w_to_wi[w2];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    return (size_t)i1 * S->nchs - (size_t)i1 * (i1 - 1) / 2 + (i2 - i1);
}
static inline size_t lvl_tri(const Level *S, int w1, int w2, int w3) {
    int i1 = S->w_to_wi[w1], i2 = S->w_to_wi[w2], i3 = S->w_to_wi[w3];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    if (i2 > i3) { int t = i2; i2 = i3; i3 = t; }
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    return S->tri_off1[i1] + S->tri_off2[(size_t)i1 * S->nchs + i2] + (size_t)(i3 - i2);
}

static void init_level(Level *S, int D0, int KT, int build_tri) {
    memset(S, 0, sizeof(*S));
    S->D0 = D0; S->KT = KT;
    memset(S->lut, -1, sizeof(S->lut));
    int cnt = 0;
    for (int D = 2; D < D0; D++)
        for (int k = 0; k <= D - 1; k++) S->lut[k][D - 1 - k] = cnt++;
    S->INNER_N = cnt; S->INF_ID = cnt; S->TAIL_BASE = cnt + 1;
    S->NG = S->TAIL_BASE + (KT + 2);
    S->gL = calloc(S->NG, sizeof(int)); S->gK = calloc(S->NG, sizeof(int));
    S->gD = calloc(S->NG, sizeof(int)); S->gKind = calloc(S->NG, sizeof(int));
    S->invD = calloc(S->NG, sizeof(double));
    for (int D = 2; D < D0; D++)
        for (int k = 0; k <= D - 1; k++) {
            int L = D - 1 - k, id = S->lut[k][L];
            S->gL[id]=L; S->gK[id]=k; S->gD[id]=D; S->gKind[id]=0; S->invD[id]=1.0/D;
        }
    S->gKind[S->INF_ID]=1; S->gD[S->INF_ID]=0; S->gK[S->INF_ID]=-1; S->invD[S->INF_ID]=0.0;
    for (int kap = 0; kap <= KT + 1; kap++) {
        int id = S->TAIL_BASE + kap;
        S->gKind[id]=2; S->gK[id]=kap; S->gD[id]=D0; S->gL[id]=-1; S->invD[id]=1.0/D0;
    }
    S->pars = malloc(S->NG * sizeof(int)); S->chs = malloc(S->NG * sizeof(int));
    S->p_to_pi = malloc(S->NG * sizeof(int)); S->w_to_wi = malloc(S->NG * sizeof(int));
    for (int i = 0; i < S->NG; i++) { S->p_to_pi[i] = -1; S->w_to_wi[i] = -1; }
    for (int u = 0; u < S->NG; u++) {
        if (S->gKind[u] == 1 || S->gK[u] >= 1) { S->p_to_pi[u] = S->npars; S->pars[S->npars++] = u; }
        S->w_to_wi[u] = S->nchs; S->chs[S->nchs++] = u;
    }
    S->pair_n = (size_t)S->nchs * (S->nchs + 1) / 2;
    if (build_tri) {
        int n = S->nchs;
        S->tri_off1 = calloc(n + 1, sizeof(size_t));
        S->tri_off2 = calloc((size_t)n * n, sizeof(size_t));
        size_t tot = 0;
        for (int i1 = 0; i1 < n; i1++) {
            S->tri_off1[i1] = tot; size_t sub = 0;
            for (int i2 = i1; i2 < n; i2++) { S->tri_off2[(size_t)i1 * n + i2] = sub; sub += (size_t)(n - i2); }
            tot += sub;
        }
        S->tri_n = tot;
    }
    S->max_rank = D0 + 1; /* inner: D = 2..D0-1; T1: D0; INF: D0+1 */
    for (int r = 0; r <= S->max_rank; r++) { S->p_by_r[r] = malloc(S->NG * sizeof(int)); S->c_by_r[r] = malloc(S->NG * sizeof(int)); }
    for (int ip = 0; ip < S->npars; ip++) { int p = S->pars[ip], r = (S->gKind[p]==0)?S->gD[p]:((S->gKind[p]==2)?D0:(D0+1)); S->p_by_r[r][S->p_rcnt[r]++] = p; }
    for (int iw = 0; iw < S->nchs;  iw++) { int w = S->chs[iw],  r = (S->gKind[w]==0)?S->gD[w]:((S->gKind[w]==2)?D0:(D0+1)); S->c_by_r[r][S->c_rcnt[r]++] = w; }
    for (int r = 2; r < S->max_rank; r++) {
        for (int i=0;i<S->p_rcnt[r];i++) for (int j=i+1;j<S->p_rcnt[r];j++) if (S->gK[S->p_by_r[r][j]] > S->gK[S->p_by_r[r][i]]) { int t=S->p_by_r[r][i];S->p_by_r[r][i]=S->p_by_r[r][j];S->p_by_r[r][j]=t; }
        for (int i=0;i<S->c_rcnt[r];i++) for (int j=i+1;j<S->c_rcnt[r];j++) if (S->gK[S->c_by_r[r][j]] > S->gK[S->c_by_r[r][i]]) { int t=S->c_by_r[r][i];S->c_by_r[r][i]=S->c_by_r[r][j];S->c_by_r[r][j]=t; }
    }
}

static void free_level(Level *S) {
    free(S->gL); free(S->gK); free(S->gD); free(S->gKind); free(S->invD);
    free(S->pars); free(S->chs); free(S->p_to_pi); free(S->w_to_wi);
    free(S->tri_off1); free(S->tri_off2);
    for (int r = 0; r <= S->max_rank; r++) { free(S->p_by_r[r]); free(S->c_by_r[r]); }
}

static inline void get_mu(const Level *S, int L, int k, int ids[2], double wts[2], int *nc) {
    int D = L + k + 1;
    if (D < S->D0) { ids[0] = lvl_inner(S, L, k); wts[0] = 1.0; *nc = 1; return; }
    double t = (double)S->D0 / D;
    ids[0] = lvl_tail(S, k); wts[0] = t;
    ids[1] = S->INF_ID;      wts[1] = 1.0 - t;
    *nc = 2;
}

static inline void proj_state(const Level *F, const Level *C, int sf, int ids[2], double wts[2], int *nc) {
    if (F->gKind[sf] == 1) { ids[0] = C->INF_ID; wts[0] = 1.0; *nc = 1; return; }
    int D = F->gD[sf], k = F->gK[sf];
    if (F->gKind[sf] == 0 && D < C->D0) { ids[0] = lvl_inner(C, F->gL[sf], k); wts[0] = 1.0; *nc = 1; return; }
    double t = (double)C->D0 / D;
    ids[0] = lvl_tail(C, k); wts[0] = t;
    ids[1] = C->INF_ID;      wts[1] = 1.0 - t;
    *nc = 2;
}

typedef struct {
    int DF, NV, *f_off;
    Level L2, L3, L4;
    int DHU, KHU, DCU, KCU;
    int DP4_1, DP4_2, DP4_S, KP4_1, KP4_2;
    int DF4_1, DF4_2, DF4_S, KF4_1, KF4_2;
    int DY4, KY4, DT4, KT4;
    int *hu_list, *hu_to_act, NHU, NHU_ACT;
    int *cu_list, *cu_to_act, NCU, NCU_ACT;
    int *p4_map, *p4_u1, *p4_u2, NP4_ACT;
    int *f4_map, *f4_u1, *f4_u2, NF4_ACT;
    int *y4_map, *y4_u, NY4_ACT;
    int *t4_map, *t4_u, NT4_ACT;
} Layout;

static inline int vidx(const Layout *C, int L, int k) { return C->f_off[L + k + 1] + k; }
static inline int hu_idx(const Layout *C, int L, int k) {
    int D = L + k + 1; if (k < 1 || k > C->KHU || D >= C->DHU) return -1;
    return (k - 1) * C->DHU + L;
}
static inline int cu_idx(const Layout *C, int L, int k) {
    int D = L + k + 1; if (k < 2 || k > C->KCU || D >= C->DCU) return -1;
    return (k - 2) * C->DCU + L;
}
static inline int p4_act(const Layout *C, int L1, int k1, int L2, int k2) {
    int u1 = hu_idx(C, L1, k1), u2 = hu_idx(C, L2, k2);
    if (u1 < 0 || u2 < 0) return -1;
    return C->p4_map[C->hu_to_act[u1] * C->NHU_ACT + C->hu_to_act[u2]];
}
static inline int f4_act(const Layout *C, int L1, int k1, int L2, int k2) {
    int c1 = cu_idx(C, L1, k1), u2 = hu_idx(C, L2, k2);
    if (c1 < 0 || u2 < 0) return -1;
    return C->f4_map[C->cu_to_act[c1] * C->NHU_ACT + C->hu_to_act[u2]];
}
static inline int y4_act(const Layout *C, int L, int k) {
    int u = hu_idx(C, L, k); return (u >= 0) ? C->y4_map[C->hu_to_act[u]] : -1;
}
static inline int t4_act(const Layout *C, int L, int k) {
    int c = cu_idx(C, L, k); return (c >= 0) ? C->t4_map[C->cu_to_act[c]] : -1;
}
static inline size_t hv_off(const Layout *C, int p3, int hu_d, int w3) {
    return ((size_t)C->L3.p_to_pi[p3] * C->NHU_ACT + C->hu_to_act[hu_d]) * C->L3.nchs + C->L3.w_to_wi[w3];
}
static inline size_t hc_off(const Layout *C, int cu_d, int w1, int w2) {
    return (size_t)C->cu_to_act[cu_d] * C->L3.pair_n + lvl_pair(&C->L3, w1, w2);
}
static inline size_t p4_off(const Layout *C, int p4a, int a4, int b4) {
    return ((size_t)p4a * C->L4.npars + C->L4.p_to_pi[a4]) * C->L4.nchs + C->L4.w_to_wi[b4];
}
static inline size_t f4_off(const Layout *C, int f4a, int a4, int b4) {
    return ((size_t)f4a * C->L4.nchs + C->L4.w_to_wi[a4]) * C->L4.nchs + C->L4.w_to_wi[b4];
}
static inline size_t y4_off(const Layout *C, int y4a, int a4, int b1, int b2) {
    return ((size_t)y4a * C->L4.npars + C->L4.p_to_pi[a4]) * C->L4.pair_n + lvl_pair(&C->L4, b1, b2);
}
static inline size_t t4_off(const Layout *C, int t4a, int b1, int b2, int b3) {
    return (size_t)t4a * C->L4.tri_n + lvl_tri(&C->L4, b1, b2, b3);
}

static void init_layout(Layout *C, int DF, int D0_g, int D0_3, int DHU, int DCU, int D0_4) {
    memset(C, 0, sizeof(*C));
    C->DF = DF;
    C->f_off = calloc(DF + 2, sizeof(int));
    int nv = 0; for (int D = 2; D < DF; D++) { C->f_off[D] = nv; nv += D; }
    C->NV = nv;

    init_level(&C->L2, D0_g, 24, 0);
    init_level(&C->L3, D0_3, 12, 0);
    init_level(&C->L4, D0_4,  8, 1);

    C->DHU = DHU; C->KHU = 12; C->DCU = DCU; C->KCU = 10;
    C->NHU = C->KHU * C->DHU;
    C->hu_to_act = malloc(C->NHU * sizeof(int)); C->hu_list = malloc(C->NHU * sizeof(int));
    for (int i = 0; i < C->NHU; i++) C->hu_to_act[i] = -1;
    for (int k = 1; k <= C->KHU; k++) for (int L = 0; L < C->DHU; L++) {
        int ui = hu_idx(C, L, k);
        if (ui >= 0) { C->hu_to_act[ui] = C->NHU_ACT; C->hu_list[C->NHU_ACT++] = ui; }
    }
    C->NCU = (C->KCU - 1) * C->DCU;
    C->cu_to_act = malloc(C->NCU * sizeof(int)); C->cu_list = malloc(C->NCU * sizeof(int));
    for (int i = 0; i < C->NCU; i++) C->cu_to_act[i] = -1;
    for (int k = 2; k <= C->KCU; k++) for (int L = 0; L < C->DCU; L++) {
        int ci = cu_idx(C, L, k);
        if (ci >= 0) { C->cu_to_act[ci] = C->NCU_ACT; C->cu_list[C->NCU_ACT++] = ci; }
    }

    /* 4-vertex domains strictly inside D < D0_4! */
    C->DP4_1 = 11; C->DP4_2 = 11; C->DP4_S = 16; C->KP4_1 = 7; C->KP4_2 = 7;
    C->DF4_1 = 10; C->DF4_2 = 10; C->DF4_S = 15; C->KF4_1 = 6; C->KF4_2 = 6;
    C->DY4 = 10; C->KY4 = 6;
    C->DT4 =  8; C->KT4 = 5;
    if (C->DP4_1 > D0_4) C->DP4_1 = D0_4; if (C->DP4_2 > D0_4) C->DP4_2 = D0_4;
    if (C->DF4_1 > D0_4) C->DF4_1 = D0_4; if (C->DF4_2 > D0_4) C->DF4_2 = D0_4;
    if (C->DY4 > D0_4) C->DY4 = D0_4; if (C->DT4 > D0_4) C->DT4 = D0_4;

    C->p4_map = malloc((size_t)C->NHU_ACT * C->NHU_ACT * sizeof(int));
    for (int i = 0; i < C->NHU_ACT * C->NHU_ACT; i++) C->p4_map[i] = -1;
    C->p4_u1 = malloc((size_t)C->NHU_ACT * C->NHU_ACT * sizeof(int));
    C->p4_u2 = malloc((size_t)C->NHU_ACT * C->NHU_ACT * sizeof(int));
    for (int a1 = 0; a1 < C->NHU_ACT; a1++) {
        int ui1 = C->hu_list[a1], k1 = (ui1/DHU)+1, L1 = ui1%DHU, D1 = L1+k1+1;
        if (D1 >= C->DP4_1 || k1 > C->KP4_1) continue;
        for (int a2 = 0; a2 < C->NHU_ACT; a2++) {
            int ui2 = C->hu_list[a2], k2 = (ui2/DHU)+1, L2 = ui2%DHU, D2 = L2+k2+1;
            if (D2 >= C->DP4_2 || k2 > C->KP4_2 || D1 + D2 >= C->DP4_S) continue;
            int id = C->NP4_ACT++;
            C->p4_map[a1 * C->NHU_ACT + a2] = id;
            C->p4_u1[id] = lvl_inner(&C->L4, L1, k1); C->p4_u2[id] = lvl_inner(&C->L4, L2, k2);
        }
    }
    C->f4_map = malloc((size_t)C->NCU_ACT * C->NHU_ACT * sizeof(int));
    for (int i = 0; i < C->NCU_ACT * C->NHU_ACT; i++) C->f4_map[i] = -1;
    C->f4_u1 = malloc((size_t)C->NCU_ACT * C->NHU_ACT * sizeof(int));
    C->f4_u2 = malloc((size_t)C->NCU_ACT * C->NHU_ACT * sizeof(int));
    for (int c1 = 0; c1 < C->NCU_ACT; c1++) {
        int ci1 = C->cu_list[c1], k1 = (ci1/DCU)+2, L1 = ci1%DCU, D1 = L1+k1+1;
        if (D1 >= C->DF4_1 || k1 > C->KF4_1) continue;
        for (int a2 = 0; a2 < C->NHU_ACT; a2++) {
            int ui2 = C->hu_list[a2], k2 = (ui2/DHU)+1, L2 = ui2%DHU, D2 = L2+k2+1;
            if (D2 >= C->DF4_2 || k2 > C->KF4_2 || D1 + D2 >= C->DF4_S) continue;
            int id = C->NF4_ACT++;
            C->f4_map[c1 * C->NHU_ACT + a2] = id;
            C->f4_u1[id] = lvl_inner(&C->L4, L1, k1); C->f4_u2[id] = lvl_inner(&C->L4, L2, k2);
        }
    }
    C->y4_map = malloc(C->NHU_ACT * sizeof(int));
    for (int i = 0; i < C->NHU_ACT; i++) C->y4_map[i] = -1;
    C->y4_u = malloc(C->NHU_ACT * sizeof(int));
    for (int a1 = 0; a1 < C->NHU_ACT; a1++) {
        int ui1 = C->hu_list[a1], k1 = (ui1/DHU)+1, L1 = ui1%DHU, D1 = L1+k1+1;
        if (k1 < 2 || k1 > C->KY4 || D1 >= C->DY4 || cu_idx(C, L1, k1) < 0) continue;
        int id = C->NY4_ACT++; C->y4_map[a1] = id; C->y4_u[id] = lvl_inner(&C->L4, L1, k1);
    }
    C->t4_map = malloc(C->NCU_ACT * sizeof(int));
    for (int i = 0; i < C->NCU_ACT; i++) C->t4_map[i] = -1;
    C->t4_u = malloc(C->NCU_ACT * sizeof(int));
    for (int c1 = 0; c1 < C->NCU_ACT; c1++) {
        int ci1 = C->cu_list[c1], k1 = (ci1/DCU)+2, L1 = ci1%DCU, D1 = L1+k1+1;
        if (k1 < 3 || k1 > C->KT4 || D1 >= C->DT4) continue;
        int id = C->NT4_ACT++; C->t4_map[c1] = id; C->t4_u[id] = lvl_inner(&C->L4, L1, k1);
    }
}

static void free_layout(Layout *C) {
    free(C->f_off); free_level(&C->L2); free_level(&C->L3); free_level(&C->L4);
    free(C->hu_list); free(C->hu_to_act); free(C->cu_list); free(C->cu_to_act);
    free(C->p4_map); free(C->p4_u1); free(C->p4_u2);
    free(C->f4_map); free(C->f4_u1); free(C->f4_u2);
    free(C->y4_map); free(C->y4_u); free(C->t4_map); free(C->t4_u);
}

static inline double tail_gmin(const Level *S, int kap, double inv_nbr_lo) {
    if (kap > S->KT) return 0.0;
    double num = (double)(S->D0 - kap - 1) + inv_nbr_lo;
    return (num > 0.0) ? (num / (S->D0 + 1.0)) : 0.0;
}
static inline double tail_bmin(const Level *S, int kap) {
    if (kap > S->KT) return 0.0;
    int m = S->D0 - kap - 1;
    return (m > 0) ? ((double)m / S->D0) : 0.0;
}
static inline double nbr_inv_lo(const Level *S, int nbr) { return (S->gKind[nbr] == 0) ? (1.0 / S->gD[nbr]) : 0.0; }
static inline double ext_rate(int mode, double jump, double r_lo, double r_hi) {
    return (mode == 1) ? ((jump > 0.0) ? r_hi : r_lo) : ((jump < 0.0) ? r_hi : r_lo);
}

/* FY helpers */
static inline double eval_pl(const double *KNOTS, const double *vals, double x) {
    if (x <= KNOTS[0]) return vals[0];
    if (x >= KNOTS[NKNOTS - 1]) return vals[NKNOTS - 1];
    for (int i = 0; i < NKNOTS - 1; i++)
        if (x <= KNOTS[i + 1]) { double t = (x - KNOTS[i]) / (KNOTS[i + 1] - KNOTS[i]); return vals[i]*(1-t)+vals[i+1]*t; }
    return vals[NKNOTS - 1];
}
static inline double eval_astar(const double *KNOTS, const double *vals, int mode, double y) {
    double bv = (mode == 1) ? -1e300 : 1e300;
    for (int k = 0; k < NKNOTS; k++) { double c = KNOTS[k]*y - vals[k]; if ((mode==1&&c>bv)||(mode==2&&c<bv)) bv=c; }
    return bv;
}
static void project_convex(const double *KNOTS, double *a, double sgn) {
    int m = NKNOTS - 1, bs[NKNOTS], be[NKNOTS], nb = 0;
    double sl[NKNOTS], wt[NKNOTS];
    for (int i = 0; i < m; i++) {
        double dx = KNOTS[i+1] - KNOTS[i];
        sl[nb] = sgn*(a[i+1]-a[i])/dx; wt[nb] = dx; bs[nb] = i; be[nb] = i; nb++;
        while (nb >= 2 && sl[nb-2] > sl[nb-1]) {
            double nw = wt[nb-2] + wt[nb-1];
            sl[nb-2] = (sl[nb-2]*wt[nb-2] + sl[nb-1]*wt[nb-1]) / nw;
            wt[nb-2] = nw; be[nb-2] = be[nb-1]; nb--;
        }
    }
    double cur = a[0];
    for (int b = 0; b < nb; b++)
        for (int i = bs[b]; i <= be[b]; i++) { cur += sgn*sl[b]*(KNOTS[i+1]-KNOTS[i]); a[i+1] = cur; }
}
static void opt_knots(const double *KNOTS, double *a, double sgn,
                      const double *xs_v, const double *wxs, int nx,
                      const double *ys_v, const double *wys, int ny) {
    double ym = 0, wyt = 0, wxt = 0;
    for (int j = 0; j < ny; j++) { ym += wys[j]*ys_v[j]; wyt += wys[j]; }
    for (int i = 0; i < nx; i++) wxt += wxs[i];
    if (wyt <= 0 || wxt <= 0) { for (int k = 0; k < NKNOTS; k++) a[k] = 0; return; }
    ym /= wyt;
    for (int k = 0; k < NKNOTS; k++) a[k] = ym * KNOTS[k];
    double mu[NKNOTS] = {0};
    for (int i = 0; i < nx; i++) {
        double x = xs_v[i], w = wxs[i] / wxt;
        if (x <= KNOTS[0]) { mu[0] += w; continue; }
        if (x >= KNOTS[NKNOTS - 1]) { mu[NKNOTS - 1] += w; continue; }
        for (int k = 0; k < NKNOTS - 1; k++)
            if (x <= KNOTS[k + 1]) { double t = (x - KNOTS[k]) / (KNOTS[k+1]-KNOTS[k]); mu[k]+=w*(1-t); mu[k+1]+=w*t; break; }
    }
    double yvar = 0; for (int j = 0; j < ny; j++) { double d = ys_v[j] - ym; yvar += (wys[j]/wyt)*d*d; }
    double ysd = sqrt(yvar + 1e-12), best_a[NKNOTS], best_obj = 1e30;
    memcpy(best_a, a, sizeof(best_a));
    for (int step = 0; step < 60; step++) {
        double grad[NKNOTS], obj = 0;
        for (int k = 0; k < NKNOTS; k++) { grad[k] = sgn*mu[k]; obj += sgn*mu[k]*a[k]; }
        for (int j = 0; j < ny; j++) {
            double y = ys_v[j], w = wys[j]/wyt, bv = -1e30; int bk = 0;
            for (int k = 0; k < NKNOTS; k++) { double c = sgn*(KNOTS[k]*y - a[k]); if (c > bv) { bv=c; bk=k; } }
            obj += w * bv; grad[bk] -= sgn * w;
        }
        if (obj < best_obj) { best_obj = obj; memcpy(best_a, a, sizeof(best_a)); }
        double lr = 0.25 * ysd / sqrt(1.0 + 0.3 * step);
        for (int k = 0; k < NKNOTS; k++) a[k] -= lr * grad[k];
        double c0 = a[0]; for (int k = 0; k < NKNOTS; k++) a[k] -= c0;
        project_convex(KNOTS, a, sgn);
    }
    memcpy(a, best_a, sizeof(best_a));
}

/* Exact 1D affine envelope (with pure 1-knot F & C, min t = C->D0 / F->D0 >= 0.39 -> SUPER TIGHT!) */
static void build_env_1d(const Level *F, const Level *C, int is_par, const double *val_f, double *m_c, double *r_c) {
    for (int sc = 0; sc < C->NG; sc++) { m_c[sc] = 0.0; r_c[sc] = 0.0; }
    for (int sc = 0; sc < C->INNER_N; sc++) {
        if (is_par && C->gK[sc] < 1) continue;
        m_c[sc] = val_f[lvl_inner(F, C->gL[sc], C->gK[sc])]; r_c[sc] = 0.0;
    }
    m_c[C->INF_ID] = val_f[F->INF_ID]; r_c[C->INF_ID] = 0.0;
    double min_s[64], max_s[64]; int has_s[64] = {0};
    for (int kap = 0; kap <= C->KT + 1; kap++) { min_s[kap] = 1e300; max_s[kap] = -1e300; }
    int nf_cnt = is_par ? F->npars : F->nchs;
    const int *nf_arr = is_par ? F->pars : F->chs;
    for (int ii = 0; ii < nf_cnt; ii++) {
        int sf = nf_arr[ii];
        if (F->gKind[sf] == 1 || (F->gKind[sf] == 0 && F->gD[sf] < C->D0)) continue;
        int kap = (F->gK[sf] > C->KT) ? (C->KT + 1) : F->gK[sf];
        double t = (double)C->D0 / F->gD[sf];
        double y = m_c[C->INF_ID] + (val_f[sf] - m_c[C->INF_ID]) / t;
        if (!has_s[kap]) { min_s[kap] = max_s[kap] = y; has_s[kap] = 1; }
        else { if (y < min_s[kap]) min_s[kap] = y; if (y > max_s[kap]) max_s[kap] = y; }
    }
    for (int kap = 0; kap <= C->KT + 1; kap++) {
        int sc = lvl_tail(C, kap);
        if (has_s[kap]) { m_c[sc] = 0.5*(min_s[kap]+max_s[kap]); r_c[sc] = 0.5*(max_s[kap]-min_s[kap]) + 1e-14; }
    }
}

/* Exact 2D affine envelope F -> C */
static void build_env_2d(const Level *F, const Level *C, int a_is_par,
                         const double *val_f, double *m_c, double *r_c) {
    int na_c = a_is_par ? C->npars : C->nchs, nb_c = C->nchs;
    int na_f = a_is_par ? F->npars : F->nchs, nb_f = F->nchs;
    const int *a_arr_f = a_is_par ? F->pars : F->chs;
    const int *a_to_i_c = a_is_par ? C->p_to_pi : C->w_to_wi;
    memset(m_c, 0, (size_t)na_c * nb_c * sizeof(double));
    memset(r_c, 0, (size_t)na_c * nb_c * sizeof(double));
    int inf_a_c = a_to_i_c[C->INF_ID], inf_b_c = C->w_to_wi[C->INF_ID];
    int inf_a_f = a_is_par ? F->p_to_pi[F->INF_ID] : F->w_to_wi[F->INF_ID];
    int inf_b_f = F->w_to_wi[F->INF_ID];

    int rep_a[256], rep_b[256];
    for (int ia = 0; ia < na_c; ia++) {
        int sc = a_is_par ? C->pars[ia] : C->chs[ia];
        if (C->gKind[sc] == 0) rep_a[ia] = a_is_par ? F->p_to_pi[lvl_inner(F, C->gL[sc], C->gK[sc])] : F->w_to_wi[lvl_inner(F, C->gL[sc], C->gK[sc])];
        else if (C->gKind[sc] == 1) rep_a[ia] = inf_a_f;
        else { int k = C->gK[sc], L = C->D0 - 1 - k; if (L < 0) { k = C->D0 - 1; L = 0; } int sf = lvl_inner(F, L, k); rep_a[ia] = a_is_par ? F->p_to_pi[sf] : F->w_to_wi[sf]; }
    }
    for (int ib = 0; ib < nb_c; ib++) {
        int sc = C->chs[ib];
        if (C->gKind[sc] == 0) rep_b[ib] = F->w_to_wi[lvl_inner(F, C->gL[sc], C->gK[sc])];
        else if (C->gKind[sc] == 1) rep_b[ib] = inf_b_f;
        else { int k = C->gK[sc], L = C->D0 - 1 - k; if (L < 0) { k = C->D0 - 1; L = 0; } rep_b[ib] = F->w_to_wi[lvl_inner(F, L, k)]; }
    }
    for (int ia = 0; ia < na_c; ia++)
        for (int ib = 0; ib < nb_c; ib++)
            m_c[(size_t)ia * nb_c + ib] = val_f[(size_t)rep_a[ia] * nb_f + rep_b[ib]];

    for (int ia_f = 0; ia_f < na_f; ia_f++) {
        int saf = a_arr_f[ia_f], aid[2], na; double aw[2]; proj_state(F, C, saf, aid, aw, &na);
        for (int ib_f = 0; ib_f < nb_f; ib_f++) {
            int sbf = F->chs[ib_f], bid[2], nb; double bw[2]; proj_state(F, C, sbf, bid, bw, &nb);
            if ((na == 1 && nb == 1) || (na == 2 && nb == 2)) continue;
            double mv = 0.0;
            for (int x=0;x<na;x++) for (int y=0;y<nb;y++) mv += aw[x]*bw[y]*m_c[(size_t)a_to_i_c[aid[x]]*nb_c + C->w_to_wi[bid[y]]];
            double err = fabs(val_f[(size_t)ia_f * nb_f + ib_f] - mv);
            size_t idx = (size_t)a_to_i_c[aid[0]] * nb_c + C->w_to_wi[bid[0]];
            double need = err / ((na == 1) ? bw[0] : aw[0]) + 1e-14;
            if (need > r_c[idx]) r_c[idx] = need;
        }
    }
    for (int ia_f = 0; ia_f < na_f; ia_f++) {
        int saf = a_arr_f[ia_f], aid[2], na; double aw[2]; proj_state(F, C, saf, aid, aw, &na);
        if (na != 2) continue;
        for (int ib_f = 0; ib_f < nb_f; ib_f++) {
            int sbf = F->chs[ib_f], bid[2], nb; double bw[2]; proj_state(F, C, sbf, bid, bw, &nb);
            if (nb != 2) continue;
            double mv = 0.0;
            for (int x=0;x<2;x++) for (int y=0;y<2;y++) mv += aw[x]*bw[y]*m_c[(size_t)a_to_i_c[aid[x]]*nb_c + C->w_to_wi[bid[y]]];
            double err = fabs(val_f[(size_t)ia_f * nb_f + ib_f] - mv);
            double already = aw[0]*bw[1]*r_c[(size_t)a_to_i_c[aid[0]]*nb_c + inf_b_c]
                           + aw[1]*bw[0]*r_c[(size_t)inf_a_c*nb_c + C->w_to_wi[bid[0]]];
            if (err > already) {
                size_t idx = (size_t)a_to_i_c[aid[0]] * nb_c + C->w_to_wi[bid[0]];
                double need = (err - already) / (aw[0] * bw[0]) + 1e-14;
                if (need > r_c[idx]) r_c[idx] = need;
            }
        }
    }
}

typedef struct {
    int mode; double lam; Layout C; double KNOTS[NKNOTS];
    double *rpar, *rch, *xpar, *xch, *xcc; /* indexed by L2 inner id */
    double A_up, B_up, B1_lo, B2_lo, CD_lo, eta_lo;
    double *f, *g, *hV, *hC, *hP4, *hF4, *hY4, *hT4;
    double (*k_pc)[NKNOTS], (*k_cp)[NKNOTS], (*k_cc)[NKNOTS];
    double *emp_edge;
} Model;

static inline double eval_f(const Model *M, int L, int k) {
    int D = L + k + 1;
    if (D < M->C.DF) return M->f[vidx(&M->C, L, k)];
    double q = (double)L / D;
    return (M->mode == 1) ? (M->A_up + M->B_up * q) : (M->B1_lo*q - M->B2_lo*q*q - M->CD_lo*q/D);
}
static inline double eval_hv_at_L2(const Model *M, const double *hV_arr, int p2, int ui, int w3) {
    const Layout *C = &M->C; int ids[2], nc; double wts[2]; proj_state(&C->L2, &C->L3, p2, ids, wts, &nc);
    double s = 0.0; for (int i=0;i<nc;i++) s += wts[i]*hV_arr[hv_off(C, ids[i], ui, w3)]; return s;
}
static inline double eval_hc_at_L2(const Model *M, const double *hC_arr, int ci, int w2, int nb3) {
    const Layout *C = &M->C; int ids[2], nc; double wts[2]; proj_state(&C->L2, &C->L3, w2, ids, wts, &nc);
    double s = 0.0; for (int i=0;i<nc;i++) s += wts[i]*hC_arr[hc_off(C, ci, ids[i], nb3)]; return s;
}
static inline double eval_p4_at_L3(const Model *M, const double *hP4_arr, int p4a, int a3, int nb4) {
    const Layout *C = &M->C; int ids[2], nc; double wts[2]; proj_state(&C->L3, &C->L4, a3, ids, wts, &nc);
    double s = 0.0; for (int i=0;i<nc;i++) s += wts[i]*hP4_arr[p4_off(C, p4a, ids[i], nb4)]; return s;
}
static inline double eval_f4_at_L3(const Model *M, const double *hF4_arr, int f4a, int a3, int b4) {
    const Layout *C = &M->C; int ids[2], nc; double wts[2]; proj_state(&C->L3, &C->L4, a3, ids, wts, &nc);
    double s = 0.0; for (int i=0;i<nc;i++) s += wts[i]*hF4_arr[f4_off(C, f4a, ids[i], b4)]; return s;
}
static inline double eval_f4_b_at_L3(const Model *M, const double *hF4_arr, int f4a, int a4, int b3) {
    const Layout *C = &M->C; int ids[2], nc; double wts[2]; proj_state(&C->L3, &C->L4, b3, ids, wts, &nc);
    double s = 0.0; for (int i=0;i<nc;i++) s += wts[i]*hF4_arr[f4_off(C, f4a, a4, ids[i])]; return s;
}
static inline double eval_y4_at_L3(const Model *M, const double *hY4_arr, int y4a, int a3, int w3, int nb4) {
    const Layout *C = &M->C; int aid[2], wid[2], na, nw; double aw[2], ww[2];
    proj_state(&C->L3, &C->L4, a3, aid, aw, &na); proj_state(&C->L3, &C->L4, w3, wid, ww, &nw);
    double s = 0.0; for (int i=0;i<na;i++) for (int j=0;j<nw;j++) s += aw[i]*ww[j]*hY4_arr[y4_off(C, y4a, aid[i], wid[j], nb4)]; return s;
}
static inline double eval_t4_at_L3(const Model *M, const double *hT4_arr, int t4a, int w1_3, int w2_3, int nb4) {
    const Layout *C = &M->C; int a1[2], a2[2], n1, n2; double ww1[2], ww2[2];
    proj_state(&C->L3, &C->L4, w1_3, a1, ww1, &n1); proj_state(&C->L3, &C->L4, w2_3, a2, ww2, &n2);
    double s = 0.0; for (int i=0;i<n1;i++) for (int j=0;j<n2;j++) s += ww1[i]*ww2[j]*hT4_arr[t4_off(C, t4a, a1[i], a2[j], nb4)]; return s;
}

static void solve_all_motifs(const Model *M, double lam, const double *dLf_g, const double *A_g,
                             const double *hV_old, const double *hC_old,
                             const double *hP4_old, const double *hF4_old, const double *hY4_old, const double *hT4_old,
                             double *gn, double *hVn, double *hCn, double *hP4n, double *hF4n, double *hY4n, double *hT4n) {
    const Layout *C = &M->C; int mode = M->mode; double sgn = (mode == 1) ? 1.0 : -1.0;
    const Level *S2 = &C->L2, *S3 = &C->L3, *S4 = &C->L4;
    int N2 = S2->NG, nb3 = lvl_inner(S3, 1, 0), nb4 = lvl_inner(S4, 1, 0);

    /* ===== 1. SOLVE g (Level 2, PURE 1-knot tail, ZERO spline corners!) ===== */
    memset(gn, 0, (size_t)N2 * N2 * sizeof(double));
    for (int rsum = 2 * S2->max_rank; rsum >= 4; rsum--) {
        int rp_lo = rsum - S2->max_rank; if (rp_lo < 2) rp_lo = 2;
        int rp_hi = rsum - 2; if (rp_hi > S2->max_rank) rp_hi = S2->max_rank;
        for (int rp = rp_lo; rp <= rp_hi; rp++) {
            int rw = rsum - rp;
            for (int ip = 0; ip < S2->p_rcnt[rp]; ip++) {
                int p = S2->p_by_r[rp][ip], pk = S2->gKind[p], Lp = S2->gL[p], kp = S2->gK[p];
                int p_up  = (pk==0)?((S2->gD[p]+1<S2->D0)?lvl_inner(S2,Lp+1,kp):lvl_tail(S2,kp)):S2->INF_ID;
                int p_bud = (pk==0&&Lp>=1)?lvl_inner(S2,Lp-1,kp+1):((pk==2&&kp<=S2->KT)?lvl_tail(S2,kp+1):p);
                int p_in_hu = (pk==0&&hu_idx(C,Lp,kp)>=0), p_in_cu = (pk==0&&cu_idx(C,Lp,kp)>=0);
                int has_pc = (pk==0&&!p_in_hu), has_cc = (pk==0&&kp>=2&&!p_in_cu);
                int cu_pb = (pk==0&&Lp>=1)?cu_idx(C,Lp-1,kp+1):-1; double inv_p = S2->invD[p];

                for (int iw = 0; iw < S2->c_rcnt[rw]; iw++) {
                    int w = S2->c_by_r[rw][iw], wk = S2->gKind[w], Lw = S2->gL[w], kw = S2->gK[w];
                    if (pk == 1 && wk == 1) { gn[(size_t)p*N2+w] = 0.0; continue; }
                    int w_up  = (wk==0)?((S2->gD[w]+1<S2->D0)?lvl_inner(S2,Lw+1,kw):lvl_tail(S2,kw)):S2->INF_ID;
                    int w_bud = (wk==0&&Lw>=1)?lvl_inner(S2,Lw-1,kw+1):((wk==2&&kw<=S2->KT)?lvl_tail(S2,kw+1):w);
                    int w_in_hu = (wk==0&&hu_idx(C,Lw,kw)>=0), has_cp = (wk==0&&kw>=1&&!w_in_hu);
                    int hu_wb = (wk==0&&Lw>=1)?hu_idx(C,Lw-1,kw+1):-1; double inv_w = S2->invD[w];

                    double gp_plus = gn[(size_t)p_up*N2+w], gw_plus = gn[(size_t)p*N2+w_up];
                    double gp_bud = gn[(size_t)p_bud*N2+w], gw_bud = gn[(size_t)p*N2+w_bud];
                    double p_hcnb = (cu_pb>=0)?eval_hc_at_L2(M,hC_old,cu_pb,w,nb3):0.0;
                    double w_hvnb = (hu_wb>=0)?eval_hv_at_L2(M,hV_old,p,hu_wb,nb3):0.0;

                    double num0 = A_g[(size_t)p*N2+w], den0 = lam;
                    if (pk == 0) {
                        num0 += (inv_w - M->rch[p]) * dLf_g[p];
                        double ur = Lp + (wk==0?inv_w:0.0) + (p_in_hu?M->xpar[p]:0.0) + (p_in_cu?(kp-1)*M->xcc[p]:0.0);
                        num0 += ur * gp_plus; den0 += ur;
                        if (Lp>=1) { double bp=(double)Lp/S2->gD[p]; num0+=bp*(gp_bud+p_hcnb); den0+=bp; }
                    }
                    if (wk == 0) {
                        num0 += (inv_p - M->rpar[w]) * dLf_g[w];
                        double ur = Lw + (pk==0?inv_p:0.0) + (w_in_hu?kw*M->xch[w]:0.0);
                        num0 += ur * gw_plus; den0 += ur;
                        if (Lw>=1) { double bw=(double)Lw/S2->gD[w]; num0+=bw*(gw_bud+w_hvnb); den0+=bw; }
                    }
                    double pglo = (pk==2)?tail_gmin(S2,kp,nbr_inv_lo(S2,w)):0.0;
                    double wglo = (wk==2)?tail_gmin(S2,kw,nbr_inv_lo(S2,p)):0.0;
                    double d0_pw = (pk==2&&wk==0)?(gn[(size_t)S2->INF_ID*N2+w_up]-gn[(size_t)S2->INF_ID*N2+w]):0.0;
                    double d0_wp = (wk==2&&pk==0)?(gn[(size_t)p_up*N2+S2->INF_ID]-gn[(size_t)p*N2+S2->INF_ID]):0.0;
                    double x0 = 1.0 / S2->D0;

                    double v = num0 / den0;
                    for (int ns = 0; ns < 16; ns++) {
                        double num = num0, den = den0;
                        if (pk == 2) {
                            double gr = ext_rate(mode, gp_plus-v, pglo, 1.0); num+=gr*gp_plus; den+=gr;
                            if (kp<=S2->KT) { double br=ext_rate(mode,gp_bud-v,tail_bmin(S2,kp),1.0); num+=br*gp_bud; den+=br; }
                            if (wk==0) { double d1=gw_plus-v; if((mode==1&&d1>d0_pw)||(mode==2&&d1<d0_pw)){num+=x0*gw_plus;den+=x0;} else num+=x0*d0_pw; }
                        }
                        if (wk == 2) {
                            double gr = ext_rate(mode, gw_plus-v, wglo, 1.0); num+=gr*gw_plus; den+=gr;
                            if (kw<=S2->KT) { double br=ext_rate(mode,gw_bud-v,tail_bmin(S2,kw),1.0); num+=br*gw_bud; den+=br; }
                            if (pk==0) { double d1=gp_plus-v; if((mode==1&&d1>d0_wp)||(mode==2&&d1<d0_wp)){num+=x0*gp_plus;den+=x0;} else num+=x0*d0_wp; }
                        }
                        if (has_pc) { double yp=gp_plus-v, bv=(mode==1)?-1e30:1e30; int ik=0; for(int kk=0;kk<NKNOTS;kk++){double c=M->KNOTS[kk]*yp-M->k_pc[p][kk]; if((mode==1&&c>bv)||(mode==2&&c<bv)){bv=c;ik=kk;}} num+=M->KNOTS[ik]*gp_plus-M->k_pc[p][ik]; den+=M->KNOTS[ik]; }
                        if (has_cc) { double yp=gp_plus-v, bv=(mode==1)?-1e30:1e30; int ik=0; for(int kk=0;kk<NKNOTS;kk++){double c=M->KNOTS[kk]*yp-M->k_cc[p][kk]; if((mode==1&&c>bv)||(mode==2&&c<bv)){bv=c;ik=kk;}} num+=(kp-1)*(M->KNOTS[ik]*gp_plus-M->k_cc[p][ik]); den+=(kp-1)*M->KNOTS[ik]; }
                        if (has_cp) { double yw=gw_plus-v, bv=(mode==1)?-1e30:1e30; int ik=0; for(int kk=0;kk<NKNOTS;kk++){double c=M->KNOTS[kk]*yw-M->k_cp[w][kk]; if((mode==1&&c>bv)||(mode==2&&c<bv)){bv=c;ik=kk;}} num+=kw*(M->KNOTS[ik]*gw_plus-M->k_cp[w][ik]); den+=kw*M->KNOTS[ik]; }
                        double vn = num / den; if (fabs(vn-v)<1e-13) { v=vn; break; } v=vn;
                    }
                    gn[(size_t)p*N2+w] = v;
                }
            }
        }
    }

    /* ===== BUILD EXACT 1D ENVELOPES L2 -> L3 ===== */
    double *dguw_m = calloc((size_t)C->NHU_ACT * S3->NG, sizeof(double));
    double *dguw_r = calloc((size_t)C->NHU_ACT * S3->NG, sizeof(double));
    double *dgpu_m = calloc((size_t)C->NHU_ACT * S3->NG, sizeof(double));
    double *dgpu_r = calloc((size_t)C->NHU_ACT * S3->NG, sizeof(double));
    double *tmp2 = malloc(N2 * sizeof(double));
    for (int ua = 0; ua < C->NHU_ACT; ua++) {
        int ui = C->hu_list[ua], ku = (ui / C->DHU) + 1, Lu = ui % C->DHU;
        int gu2 = lvl_inner(S2, Lu, ku), gu2_up = lvl_inner(S2, Lu + 1, ku);
        for (int w2=0;w2<N2;w2++) tmp2[w2] = gn[(size_t)gu2_up*N2 + w2] - gn[(size_t)gu2*N2 + w2];
        build_env_1d(S2, S3, 0, tmp2, dguw_m + (size_t)ua*S3->NG, dguw_r + (size_t)ua*S3->NG);
        for (int p2=0;p2<N2;p2++) tmp2[p2] = gn[(size_t)p2*N2 + gu2_up] - gn[(size_t)p2*N2 + gu2];
        build_env_1d(S2, S3, 1, tmp2, dgpu_m + (size_t)ua*S3->NG, dgpu_r + (size_t)ua*S3->NG);
    }
    free(tmp2);

    /* ===== 2. SOLVE hV (Level 3) ===== */
    size_t hv_len = (size_t)S3->npars * C->NHU_ACT * S3->nchs;
    memset(hVn, 0, hv_len * sizeof(double));
    int ru_max = C->DHU - 1;
    for (int rsum = 2 * S3->max_rank + ru_max; rsum >= 6; rsum--)
        for (int ru = ru_max; ru >= 2; ru--)
            for (int rp = S3->max_rank; rp >= 2; rp--) {
                int rw = rsum - ru - rp;
                if (rw < 2 || rw > S3->max_rank || !S3->p_rcnt[rp] || !S3->c_rcnt[rw]) continue;
                for (int ku = C->KHU; ku >= 1; ku--) {
                    int Lu = ru - 1 - ku; if (Lu < 0 || Lu >= C->DHU) continue;
                    int ui = hu_idx(C, Lu, ku); if (ui < 0) continue;
                    int ua = C->hu_to_act[ui], gu2 = lvl_inner(S2, Lu, ku);
                    double iu = 1.0/ru, bu = (double)Lu/ru, xp_u = M->xpar[gu2], xw_u = M->xch[gu2];
                    int ui_up = hu_idx(C, Lu+1, ku), ui_bud = (Lu>=1)?hu_idx(C, Lu-1, ku+1):-1;
                    int u_in_y4 = (y4_act(C, Lu, ku) >= 0), ub_y4 = (Lu>=1)?y4_act(C, Lu-1, ku+1):-1;
                    double ru_keep = u_in_y4 ? ((ku-1)*M->rch[gu2]) : 0.0, ru_hi = u_in_y4 ? 0.0 : (0.5*(ku-1));
                    const double *mw = dguw_m + (size_t)ua*S3->NG, *rw_env = dguw_r + (size_t)ua*S3->NG;
                    const double *mp = dgpu_m + (size_t)ua*S3->NG, *rp_env = dgpu_r + (size_t)ua*S3->NG;

                    for (int ip = 0; ip < S3->p_rcnt[rp]; ip++) {
                        int p = S3->p_by_r[rp][ip], pk = S3->gKind[p], Lp = S3->gL[p], kp = S3->gK[p];
                        int p_up  = (pk==0)?((S3->gD[p]+1<S3->D0)?lvl_inner(S3,Lp+1,kp):lvl_tail(S3,kp)):S3->INF_ID;
                        int p_bud = (pk==0&&Lp>=1)?lvl_inner(S3,Lp-1,kp+1):((pk==2&&kp<=S3->KT)?lvl_tail(S3,kp+1):p);
                        double bp = (pk==0&&Lp>=1)?((double)Lp/S3->gD[p]):0.0, inv_p = S3->invD[p];
                        int p2 = (pk==0)?lvl_inner(S2,Lp,kp):-1;
                        int pu_p4 = (pk==0)?(p4_act(C,Lp,kp,Lu,ku)>=0):0, pu_f4 = (pk==0)?(f4_act(C,Lp,kp,Lu,ku)>=0):0;
                        int pbu_f4 = (pk==0&&Lp>=1)?f4_act(C,Lp-1,kp+1,Lu,ku):-1;
                        double rp_keep = (pk==0)?((pu_p4?M->rpar[p2]:0.0)+(pu_f4?(kp-1)*M->rch[p2]:0.0)):0.0;
                        double rp_hi   = (pk==0)?(0.5*((pu_p4?0:1)+(pu_f4?0:(kp-1)))):0.0;
                        double pglo = (pk==2)?tail_gmin(S3,kp,iu):0.0;
                        for (int iw = 0; iw < S3->c_rcnt[rw]; iw++) {
                            int w = S3->c_by_r[rw][iw], wk = S3->gKind[w], Lw = S3->gL[w], kw = S3->gK[w];
                            int w_up  = (wk==0)?((S3->gD[w]+1<S3->D0)?lvl_inner(S3,Lw+1,kw):lvl_tail(S3,kw)):S3->INF_ID;
                            int w_bud = (wk==0&&Lw>=1)?lvl_inner(S3,Lw-1,kw+1):((wk==2&&kw<=S3->KT)?lvl_tail(S3,kw+1):w);
                            double bw = (wk==0&&Lw>=1)?((double)Lw/S3->gD[w]):0.0, inv_w = S3->invD[w];
                            int w2 = (wk==0)?lvl_inner(S2,Lw,kw):-1;
                            int uw_p4 = (wk==0)?(p4_act(C,Lu,ku,Lw,kw)>=0):0, uwb_p4 = (wk==0&&Lw>=1)?p4_act(C,Lu,ku,Lw-1,kw+1):-1;
                            double rw_keep = (wk==0&&uw_p4)?(kw*M->rch[w2]):0.0, rw_hi = (wk==0&&!uw_p4)?(0.5*kw):0.0;
                            double wglo = (wk==2)?tail_gmin(S3,kw,iu):0.0;
                            size_t off = hv_off(C, p, ui, w);
                            double hp_up=hVn[hv_off(C,p_up,ui,w)], hu_up=(ui_up>=0)?hVn[hv_off(C,p,ui_up,w)]:0.0, hw_up=hVn[hv_off(C,p,ui,w_up)];
                            double hp_bud=(pk==0&&Lp>=1)||(pk==2&&kp<=S3->KT)?hVn[hv_off(C,p_bud,ui,w)]:0.0;
                            double hu_bud=(ui_bud>=0)?hVn[hv_off(C,p,ui_bud,w)]:0.0;
                            double hw_bud=(wk==0&&Lw>=1)||(wk==2&&kw<=S3->KT)?hVn[hv_off(C,p,ui,w_bud)]:0.0;
                            double nb_p=(pbu_f4>=0)?eval_f4_b_at_L3(M,hF4_old,pbu_f4,nb4,w):0.0;
                            double nb_u=(ub_y4>=0) ?eval_y4_at_L3(M,hY4_old,ub_y4,p,w,nb4):0.0;
                            double nb_w=(uwb_p4>=0)?eval_p4_at_L3(M,hP4_old,uwb_p4,p,nb4):0.0;
                            double src = (inv_p-xp_u)*mw[w]+sgn*fabs(inv_p-xp_u)*rw_env[w]
                                       + (inv_w-xw_u)*mp[p]+sgn*fabs(inv_w-xw_u)*rp_env[p];
                            double num0 = src, den0 = lam;
                            if(pk==0){double r=Lp+iu+rp_keep;num0+=r*hp_up;den0+=r;if(bp>0){num0+=bp*(hp_bud+nb_p);den0+=bp;}}
                            if(wk==0){double r=Lw+iu+rw_keep;num0+=r*hw_up;den0+=r;if(bw>0){num0+=bw*(hw_bud+nb_w);den0+=bw;}}
                            {double rub=Lu+(pk==0?inv_p:0.0)+(wk==0?inv_w:0.0)+ru_keep;num0+=rub*hu_up;den0+=rub;if(bu>0){num0+=bu*(hu_bud+nb_u);den0+=bu;}}
                            double d0_pu=(pk==2)?(((ui_up>=0)?hVn[hv_off(C,S3->INF_ID,ui_up,w)]:0.0)-hp_up):0.0;
                            double d0_wu=(wk==2)?(((ui_up>=0)?hVn[hv_off(C,p,ui_up,S3->INF_ID)]:0.0)-hw_up):0.0;
                            double v=num0/den0;
                            for (int ns=0; ns<12; ns++) {
                                double num=num0, den=den0;
                                if(rp_hi>0&&((mode==1&&hp_up>v)||(mode==2&&hp_up<v))){num+=rp_hi*hp_up;den+=rp_hi;}
                                if(ru_hi>0&&((mode==1&&hu_up>v)||(mode==2&&hu_up<v))){num+=ru_hi*hu_up;den+=ru_hi;}
                                if(rw_hi>0&&((mode==1&&hw_up>v)||(mode==2&&hw_up<v))){num+=rw_hi*hw_up;den+=rw_hi;}
                                if(pk==2){double gr=ext_rate(mode,hp_up-v,pglo,1.0);num+=gr*hp_up;den+=gr;
                                    if(kp<=S3->KT){double br=ext_rate(mode,hp_bud-v,tail_bmin(S3,kp),1.0);num+=br*hp_bud;den+=br;}
                                    double x0=1.0/S3->D0, d1=hu_up-v;
                                    if((mode==1&&d1>d0_pu)||(mode==2&&d1<d0_pu)){num+=x0*hu_up;den+=x0;} else num+=x0*d0_pu;}
                                if(wk==2){double gr=ext_rate(mode,hw_up-v,wglo,1.0);num+=gr*hw_up;den+=gr;
                                    if(kw<=S3->KT){double br=ext_rate(mode,hw_bud-v,tail_bmin(S3,kw),1.0);num+=br*hw_bud;den+=br;}
                                    double x0=1.0/S3->D0, d1=hu_up-v;
                                    if((mode==1&&d1>d0_wu)||(mode==2&&d1<d0_wu)){num+=x0*hu_up;den+=x0;} else num+=x0*d0_wu;}
                                double vn=num/den; if(fabs(vn-v)<1e-13){v=vn;break;} v=vn;
                            }
                            hVn[off] = v;
                        }
                    }
                }
            }

    /* ===== 3. SOLVE hC (Level 3) ===== */
    size_t hc_len = (size_t)C->NCU_ACT * S3->pair_n;
    memset(hCn, 0, hc_len * sizeof(double));
    int rcu_max = C->DCU - 1;
    for (int rsum = rcu_max + 2 * S3->max_rank; rsum >= 6; rsum--)
        for (int ru = rcu_max; ru >= 2; ru--)
            for (int r1 = S3->max_rank; r1 >= 2; r1--) {
                int r2 = rsum - ru - r1;
                if (r2 < r1 || r2 > S3->max_rank || !S3->c_rcnt[r1] || !S3->c_rcnt[r2]) continue;
                for (int ku = C->KCU; ku >= 2; ku--) {
                    int Lu = ru - 1 - ku; if (Lu < 0 || Lu >= C->DCU) continue;
                    int ci = cu_idx(C, Lu, ku); if (ci < 0) continue;
                    int ua = C->hu_to_act[hu_idx(C, Lu, ku)], gu2 = lvl_inner(S2, Lu, ku);
                    double iu = 1.0/ru, bu = (double)Lu/ru, xc_u = M->xcc[gu2];
                    int ci_up = cu_idx(C, Lu+1, ku), ci_bud = (Lu>=1)?cu_idx(C, Lu-1, ku+1):-1;
                    int u_in_y4 = (y4_act(C, Lu, ku) >= 0), u_in_t4 = (t4_act(C, Lu, ku) >= 0);
                    int ub_t4 = (Lu >= 1) ? t4_act(C, Lu-1, ku+1) : -1;
                    double ru_keep = (u_in_y4 ? M->rpar[gu2] : 0.0) + (u_in_t4 ? (ku - 2)*M->rch[gu2] : 0.0);
                    double ru_hi   = 0.5 * ((u_in_y4 ? 0 : 1) + (u_in_t4 ? 0 : (ku - 2)));
                    const double *mw = dguw_m + (size_t)ua*S3->NG, *rw_env = dguw_r + (size_t)ua*S3->NG;
                    for (int i1 = 0; i1 < S3->c_rcnt[r1]; i1++) {
                        int w1=S3->c_by_r[r1][i1], k1k=S3->gKind[w1], L1=S3->gL[w1], k1=S3->gK[w1];
                        int w1_up=(k1k==0)?((S3->gD[w1]+1<S3->D0)?lvl_inner(S3,L1+1,k1):lvl_tail(S3,k1)):S3->INF_ID;
                        int w1_bud=(k1k==0&&L1>=1)?lvl_inner(S3,L1-1,k1+1):((k1k==2&&k1<=S3->KT)?lvl_tail(S3,k1+1):w1);
                        double b1=(k1k==0&&L1>=1)?((double)L1/S3->gD[w1]):0.0, inv1=S3->invD[w1];
                        int w1_2=(k1k==0)?lvl_inner(S2,L1,k1):-1;
                        int uw1_f4=(k1k==0)?(f4_act(C,Lu,ku,L1,k1)>=0):0, uw1b_f4=(k1k==0&&L1>=1)?f4_act(C,Lu,ku,L1-1,k1+1):-1;
                        double r1_keep=(k1k==0&&uw1_f4)?(k1*M->rch[w1_2]):0.0, r1_hi=(k1k==0&&!uw1_f4)?(0.5*k1):0.0;
                        double g1lo=(k1k==2)?tail_gmin(S3,k1,iu):0.0;
                        for (int i2 = (r1==r2?i1:0); i2 < S3->c_rcnt[r2]; i2++) {
                            int w2=S3->c_by_r[r2][i2], k2k=S3->gKind[w2], L2=S3->gL[w2], k2=S3->gK[w2];
                            int w2_up=(k2k==0)?((S3->gD[w2]+1<S3->D0)?lvl_inner(S3,L2+1,k2):lvl_tail(S3,k2)):S3->INF_ID;
                            int w2_bud=(k2k==0&&L2>=1)?lvl_inner(S3,L2-1,k2+1):((k2k==2&&k2<=S3->KT)?lvl_tail(S3,k2+1):w2);
                            double b2=(k2k==0&&L2>=1)?((double)L2/S3->gD[w2]):0.0, inv2=S3->invD[w2];
                            int w2_2=(k2k==0)?lvl_inner(S2,L2,k2):-1;
                            int uw2_f4=(k2k==0)?(f4_act(C,Lu,ku,L2,k2)>=0):0, uw2b_f4=(k2k==0&&L2>=1)?f4_act(C,Lu,ku,L2-1,k2+1):-1;
                            double r2_keep=(k2k==0&&uw2_f4)?(k2*M->rch[w2_2]):0.0, r2_hi=(k2k==0&&!uw2_f4)?(0.5*k2):0.0;
                            double g2lo=(k2k==2)?tail_gmin(S3,k2,iu):0.0;
                            size_t off = hc_off(C, ci, w1, w2);
                            double cu_up=(ci_up>=0)?hCn[hc_off(C,ci_up,w1,w2)]:0.0, cu_bud=(ci_bud>=0)?hCn[hc_off(C,ci_bud,w1,w2)]:0.0;
                            double c1_up=hCn[hc_off(C,ci,w1_up,w2)], c2_up=hCn[hc_off(C,ci,w1,w2_up)];
                            double c1_bud=(k1k==0&&L1>=1)||(k1k==2&&k1<=S3->KT)?hCn[hc_off(C,ci,w1_bud,w2)]:0.0;
                            double c2_bud=(k2k==0&&L2>=1)||(k2k==2&&k2<=S3->KT)?hCn[hc_off(C,ci,w1,w2_bud)]:0.0;
                            double nb_u=(ub_t4>=0)?eval_t4_at_L3(M,hT4_old,ub_t4,w1,w2,nb4):0.0;
                            double nb_1=(uw1b_f4>=0)?eval_f4_at_L3(M,hF4_old,uw1b_f4,w2,nb4):0.0;
                            double nb_2=(uw2b_f4>=0)?eval_f4_at_L3(M,hF4_old,uw2b_f4,w1,nb4):0.0;
                            double src = (inv1-xc_u)*mw[w2]+sgn*fabs(inv1-xc_u)*rw_env[w2]
                                       + (inv2-xc_u)*mw[w1]+sgn*fabs(inv2-xc_u)*rw_env[w1];
                            double num0 = src, den0 = lam;
                            if(k1k==0){double r=L1+iu+r1_keep;num0+=r*c1_up;den0+=r;if(b1>0){num0+=b1*(c1_bud+nb_1);den0+=b1;}}
                            if(k2k==0){double r=L2+iu+r2_keep;num0+=r*c2_up;den0+=r;if(b2>0){num0+=b2*(c2_bud+nb_2);den0+=b2;}}
                            {double rub=Lu+(k1k==0?inv1:0.0)+(k2k==0?inv2:0.0)+ru_keep;num0+=rub*cu_up;den0+=rub;if(bu>0){num0+=bu*(cu_bud+nb_u);den0+=bu;}}
                            double d0_1=(k1k==2)?(((ci_up>=0)?hCn[hc_off(C,ci_up,S3->INF_ID,w2)]:0.0)-c1_up):0.0;
                            double d0_2=(k2k==2)?(((ci_up>=0)?hCn[hc_off(C,ci_up,w1,S3->INF_ID)]:0.0)-c2_up):0.0;
                            double v=num0/den0;
                            for (int ns=0; ns<12; ns++) {
                                double num=num0, den=den0;
                                if(ru_hi>0&&((mode==1&&cu_up>v)||(mode==2&&cu_up<v))){num+=ru_hi*cu_up;den+=ru_hi;}
                                if(r1_hi>0&&((mode==1&&c1_up>v)||(mode==2&&c1_up<v))){num+=r1_hi*c1_up;den+=r1_hi;}
                                if(r2_hi>0&&((mode==1&&c2_up>v)||(mode==2&&c2_up<v))){num+=r2_hi*c2_up;den+=r2_hi;}
                                if(k1k==2){double gr=ext_rate(mode,c1_up-v,g1lo,1.0);num+=gr*c1_up;den+=gr;
                                    if(k1<=S3->KT){double br=ext_rate(mode,c1_bud-v,tail_bmin(S3,k1),1.0);num+=br*c1_bud;den+=br;}
                                    double x0=1.0/S3->D0, d1=cu_up-v;
                                    if((mode==1&&d1>d0_1)||(mode==2&&d1<d0_1)){num+=x0*cu_up;den+=x0;} else num+=x0*d0_1;}
                                if(k2k==2){double gr=ext_rate(mode,c2_up-v,g2lo,1.0);num+=gr*c2_up;den+=gr;
                                    if(k2<=S3->KT){double br=ext_rate(mode,c2_bud-v,tail_bmin(S3,k2),1.0);num+=br*c2_bud;den+=br;}
                                    double x0=1.0/S3->D0, d1=cu_up-v;
                                    if((mode==1&&d1>d0_2)||(mode==2&&d1<d0_2)){num+=x0*cu_up;den+=x0;} else num+=x0*d0_2;}
                                double vn=num/den; if(fabs(vn-v)<1e-13){v=vn;break;} v=vn;
                            }
                            hCn[off] = v;
                        }
                    }
                }
            }
    free(dguw_m); free(dguw_r); free(dgpu_m); free(dgpu_r);

    /* ===== 4. SOLVE hP4 (Level 4) ===== */
    if (C->NP4_ACT > 0) {
        memset(hP4n, 0, (size_t)C->NP4_ACT * S4->npars * S4->nchs * sizeof(double));
        int *ord = malloc(C->NP4_ACT * sizeof(int));
        for (int i=0;i<C->NP4_ACT;i++) ord[i]=i;
        for (int i=0;i<C->NP4_ACT;i++) for (int j=i+1;j<C->NP4_ACT;j++) {
            int i1=ord[i],i2=ord[j],d1=S4->gD[C->p4_u1[i1]]+S4->gD[C->p4_u2[i1]],d2=S4->gD[C->p4_u1[i2]]+S4->gD[C->p4_u2[i2]];
            int k1=S4->gK[C->p4_u1[i1]]+S4->gK[C->p4_u2[i1]],k2=S4->gK[C->p4_u1[i2]]+S4->gK[C->p4_u2[i2]];
            if (d2>d1 || (d2==d1 && k2>k1)) { ord[i]=i2; ord[j]=i1; }
        }
        double *tmp3=malloc(S3->NG*sizeof(double)), *m1=malloc(S4->NG*sizeof(double)), *r1=malloc(S4->NG*sizeof(double));
        double *m2=malloc(S4->NG*sizeof(double)), *r2=malloc(S4->NG*sizeof(double));
        for (int oi = 0; oi < C->NP4_ACT; oi++) {
            int p4a=ord[oi], gu1_4=C->p4_u1[p4a], gu2_4=C->p4_u2[p4a];
            int L1=S4->gL[gu1_4],k1=S4->gK[gu1_4],D1=S4->gD[gu1_4], L2=S4->gL[gu2_4],k2=S4->gK[gu2_4],D2=S4->gD[gu2_4];
            int ui1=hu_idx(C,L1,k1), ui2=hu_idx(C,L2,k2);
            int gu1_3=lvl_inner(S3,L1,k1), gu1_3u=lvl_inner(S3,L1+1,k1), gu2_3=lvl_inner(S3,L2,k2), gu2_3u=lvl_inner(S3,L2+1,k2);
            double i1=1.0/D1, i2=1.0/D2, b1=(double)L1/D1, b2=(double)L2/D2, r1_hi=0.5*(k1-1), r2_hi=0.5*(k2-1);
            double z_up=M->rpar[lvl_inner(S2,L1,k1)], z_dn=M->rch[lvl_inner(S2,L2,k2)];
            int p4_1u=p4_act(C,L1+1,k1,L2,k2), p4_1b=(L1>=1)?p4_act(C,L1-1,k1+1,L2,k2):-1;
            int p4_2u=p4_act(C,L1,k1,L2+1,k2), p4_2b=(L2>=1)?p4_act(C,L1,k1,L2-1,k2+1):-1;
            for (int b3=0;b3<S3->NG;b3++) tmp3[b3]=hVn[hv_off(C,gu1_3u,ui2,b3)]-hVn[hv_off(C,gu1_3,ui2,b3)];
            build_env_1d(S3,S4,0,tmp3,m1,r1);
            memset(tmp3,0,S3->NG*sizeof(double));
            for (int ia=0;ia<S3->npars;ia++){int a3=S3->pars[ia]; tmp3[a3]=hVn[hv_off(C,a3,ui1,gu2_3u)]-hVn[hv_off(C,a3,ui1,gu2_3)];}
            build_env_1d(S3,S4,1,tmp3,m2,r2);
            for (int rsum=2*S4->max_rank;rsum>=4;rsum--) {
                int ra_lo=rsum-S4->max_rank; if(ra_lo<2)ra_lo=2; int ra_hi=rsum-2; if(ra_hi>S4->max_rank)ra_hi=S4->max_rank;
                for (int ra=ra_lo;ra<=ra_hi;ra++) {
                    int rb=rsum-ra;
                    for (int ia=0;ia<S4->p_rcnt[ra];ia++) {
                        int a=S4->p_by_r[ra][ia],ak=S4->gKind[a],La=S4->gL[a],ka=S4->gK[a];
                        int a_up=(ak==0)?((S4->gD[a]+1<S4->D0)?lvl_inner(S4,La+1,ka):lvl_tail(S4,ka)):S4->INF_ID;
                        int a_bud=(ak==0&&La>=1)?lvl_inner(S4,La-1,ka+1):((ak==2&&ka<=S4->KT)?lvl_tail(S4,ka+1):a);
                        double ba=(ak==0&&La>=1)?((double)La/S4->gD[a]):0.0, ra_hi_e=(ak==0)?(0.5*ka):0.0, inv_a=S4->invD[a];
                        double aglo=(ak==2)?tail_gmin(S4,ka,i1):0.0;
                        for (int ib=0;ib<S4->c_rcnt[rb];ib++) {
                            int b=S4->c_by_r[rb][ib],bk=S4->gKind[b],Lb=S4->gL[b],kb=S4->gK[b];
                            int b_up=(bk==0)?((S4->gD[b]+1<S4->D0)?lvl_inner(S4,Lb+1,kb):lvl_tail(S4,kb)):S4->INF_ID;
                            int b_bud=(bk==0&&Lb>=1)?lvl_inner(S4,Lb-1,kb+1):((bk==2&&kb<=S4->KT)?lvl_tail(S4,kb+1):b);
                            double bb=(bk==0&&Lb>=1)?((double)Lb/S4->gD[b]):0.0, rb_hi_e=(bk==0)?(0.5*kb):0.0, inv_b=S4->invD[b];
                            double bglo=(bk==2)?tail_gmin(S4,kb,i2):0.0;
                            size_t off=p4_off(C,p4a,a,b);
                            double va_up=hP4n[p4_off(C,p4a,a_up,b)], vb_up=hP4n[p4_off(C,p4a,a,b_up)];
                            double va_bud=(ak==0&&La>=1)||(ak==2&&ka<=S4->KT)?hP4n[p4_off(C,p4a,a_bud,b)]:0.0;
                            double vb_bud=(bk==0&&Lb>=1)||(bk==2&&kb<=S4->KT)?hP4n[p4_off(C,p4a,a,b_bud)]:0.0;
                            double v1_up=(p4_1u>=0)?hP4n[p4_off(C,p4_1u,a,b)]:0.0, v1_bud=(p4_1b>=0)?hP4n[p4_off(C,p4_1b,a,b)]:0.0;
                            double v2_up=(p4_2u>=0)?hP4n[p4_off(C,p4_2u,a,b)]:0.0, v2_bud=(p4_2b>=0)?hP4n[p4_off(C,p4_2b,a,b)]:0.0;
                            double num0=(inv_a-z_up)*m1[b]+sgn*fabs(inv_a-z_up)*r1[b]+(inv_b-z_dn)*m2[a]+sgn*fabs(inv_b-z_dn)*r2[a], den0=lam;
                            if(ak==0){double r=La+i1;num0+=r*va_up;den0+=r;if(ba>0){num0+=ba*va_bud;den0+=ba;}}
                            if(bk==0){double r=Lb+i2;num0+=r*vb_up;den0+=r;if(bb>0){num0+=bb*vb_bud;den0+=bb;}}
                            {double r1b=L1+(ak==0?inv_a:0.0)+i2;num0+=r1b*v1_up;den0+=r1b;if(b1>0){num0+=b1*v1_bud;den0+=b1;}}
                            {double r2b=L2+i1+(bk==0?inv_b:0.0);num0+=r2b*v2_up;den0+=r2b;if(b2>0){num0+=b2*v2_bud;den0+=b2;}}
                            double d0_a=(ak==2)?(((p4_1u>=0)?hP4n[p4_off(C,p4_1u,S4->INF_ID,b)]:0.0)-va_up):0.0;
                            double d0_b=(bk==2)?(((p4_2u>=0)?hP4n[p4_off(C,p4_2u,a,S4->INF_ID)]:0.0)-vb_up):0.0;
                            double v=num0/den0;
                            for (int ns=0;ns<12;ns++) {
                                double num=num0, den=den0;
                                if(ra_hi_e>0&&((mode==1&&va_up>v)||(mode==2&&va_up<v))){num+=ra_hi_e*va_up;den+=ra_hi_e;}
                                if(rb_hi_e>0&&((mode==1&&vb_up>v)||(mode==2&&vb_up<v))){num+=rb_hi_e*vb_up;den+=rb_hi_e;}
                                if(r1_hi>0  &&((mode==1&&v1_up>v)||(mode==2&&v1_up<v))){num+=r1_hi*v1_up;den+=r1_hi;}
                                if(r2_hi>0  &&((mode==1&&v2_up>v)||(mode==2&&v2_up<v))){num+=r2_hi*v2_up;den+=r2_hi;}
                                if(ak==2){double gr=ext_rate(mode,va_up-v,aglo,1.0);num+=gr*va_up;den+=gr;
                                    if(ka<=S4->KT){double br=ext_rate(mode,va_bud-v,tail_bmin(S4,ka),1.0);num+=br*va_bud;den+=br;}
                                    double x0=1.0/S4->D0, d1=v1_up-v;
                                    if((mode==1&&d1>d0_a)||(mode==2&&d1<d0_a)){num+=x0*v1_up;den+=x0;} else num+=x0*d0_a;}
                                if(bk==2){double gr=ext_rate(mode,vb_up-v,bglo,1.0);num+=gr*vb_up;den+=gr;
                                    if(kb<=S4->KT){double br=ext_rate(mode,vb_bud-v,tail_bmin(S4,kb),1.0);num+=br*vb_bud;den+=br;}
                                    double x0=1.0/S4->D0, d1=v2_up-v;
                                    if((mode==1&&d1>d0_b)||(mode==2&&d1<d0_b)){num+=x0*v2_up;den+=x0;} else num+=x0*d0_b;}
                                double vn=num/den; if(fabs(vn-v)<1e-13){v=vn;break;} v=vn;
                            }
                            hP4n[off]=v;
                        }
                    }
                }
            }
        }
        free(ord); free(tmp3); free(m1); free(r1); free(m2); free(r2);
    }

    /* ===== 5. SOLVE hF4 (Level 4) ===== */
    if (C->NF4_ACT > 0) {
        memset(hF4n, 0, (size_t)C->NF4_ACT * S4->nchs * S4->nchs * sizeof(double));
        int *ord = malloc(C->NF4_ACT * sizeof(int));
        for (int i=0;i<C->NF4_ACT;i++) ord[i]=i;
        for (int i=0;i<C->NF4_ACT;i++) for (int j=i+1;j<C->NF4_ACT;j++) {
            int i1=ord[i],i2=ord[j],d1=S4->gD[C->f4_u1[i1]]+S4->gD[C->f4_u2[i1]],d2=S4->gD[C->f4_u1[i2]]+S4->gD[C->f4_u2[i2]];
            int k1=S4->gK[C->f4_u1[i1]]+S4->gK[C->f4_u2[i1]],k2=S4->gK[C->f4_u1[i2]]+S4->gK[C->f4_u2[i2]];
            if (d2>d1 || (d2==d1 && k2>k1)) { ord[i]=i2; ord[j]=i1; }
        }
        double *tmp3=malloc(S3->NG*sizeof(double)), *m1=malloc(S4->NG*sizeof(double)), *r1=malloc(S4->NG*sizeof(double));
        double *m2=malloc(S4->NG*sizeof(double)), *r2=malloc(S4->NG*sizeof(double));
        for (int oi = 0; oi < C->NF4_ACT; oi++) {
            int f4a=ord[oi], gu1_4=C->f4_u1[f4a], gu2_4=C->f4_u2[f4a];
            int L1=S4->gL[gu1_4],k1=S4->gK[gu1_4],D1=S4->gD[gu1_4], L2=S4->gL[gu2_4],k2=S4->gK[gu2_4],D2=S4->gD[gu2_4];
            int ci1=cu_idx(C,L1,k1), ui2=hu_idx(C,L2,k2);
            int gu1_3=lvl_inner(S3,L1,k1), gu1_3u=lvl_inner(S3,L1+1,k1), gu2_3=lvl_inner(S3,L2,k2), gu2_3u=lvl_inner(S3,L2+1,k2);
            double i1=1.0/D1, i2=1.0/D2, b1=(double)L1/D1, b2=(double)L2/D2, r1_hi=0.5*(k1-1), r2_hi=0.5*(k2-1);
            double z_f1=M->rch[lvl_inner(S2,L1,k1)], z_f2=M->rch[lvl_inner(S2,L2,k2)];
            int f4_1u=f4_act(C,L1+1,k1,L2,k2), f4_1b=(L1>=1)?f4_act(C,L1-1,k1+1,L2,k2):-1;
            int f4_2u=f4_act(C,L1,k1,L2+1,k2), f4_2b=(L2>=1)?f4_act(C,L1,k1,L2-1,k2+1):-1;
            for (int b3=0;b3<S3->NG;b3++) tmp3[b3]=hVn[hv_off(C,gu1_3u,ui2,b3)]-hVn[hv_off(C,gu1_3,ui2,b3)];
            build_env_1d(S3,S4,0,tmp3,m1,r1);
            for (int a3=0;a3<S3->NG;a3++) tmp3[a3]=hCn[hc_off(C,ci1,a3,gu2_3u)]-hCn[hc_off(C,ci1,a3,gu2_3)];
            build_env_1d(S3,S4,0,tmp3,m2,r2);
            for (int rsum=2*S4->max_rank;rsum>=4;rsum--) {
                int ra_lo=rsum-S4->max_rank; if(ra_lo<2)ra_lo=2; int ra_hi=rsum-2; if(ra_hi>S4->max_rank)ra_hi=S4->max_rank;
                for (int ra=ra_lo;ra<=ra_hi;ra++) {
                    int rb=rsum-ra;
                    for (int ia=0;ia<S4->c_rcnt[ra];ia++) {
                        int a=S4->c_by_r[ra][ia],ak=S4->gKind[a],La=S4->gL[a],ka=S4->gK[a];
                        int a_up=(ak==0)?((S4->gD[a]+1<S4->D0)?lvl_inner(S4,La+1,ka):lvl_tail(S4,ka)):S4->INF_ID;
                        int a_bud=(ak==0&&La>=1)?lvl_inner(S4,La-1,ka+1):((ak==2&&ka<=S4->KT)?lvl_tail(S4,ka+1):a);
                        double ba=(ak==0&&La>=1)?((double)La/S4->gD[a]):0.0, ra_hi_e=(ak==0)?(0.5*ka):0.0, inv_a=S4->invD[a];
                        double aglo=(ak==2)?tail_gmin(S4,ka,i1):0.0;
                        for (int ib=0;ib<S4->c_rcnt[rb];ib++) {
                            int b=S4->c_by_r[rb][ib],bk=S4->gKind[b],Lb=S4->gL[b],kb=S4->gK[b];
                            int b_up=(bk==0)?((S4->gD[b]+1<S4->D0)?lvl_inner(S4,Lb+1,kb):lvl_tail(S4,kb)):S4->INF_ID;
                            int b_bud=(bk==0&&Lb>=1)?lvl_inner(S4,Lb-1,kb+1):((bk==2&&kb<=S4->KT)?lvl_tail(S4,kb+1):b);
                            double bb=(bk==0&&Lb>=1)?((double)Lb/S4->gD[b]):0.0, rb_hi_e=(bk==0)?(0.5*kb):0.0, inv_b=S4->invD[b];
                            double bglo=(bk==2)?tail_gmin(S4,kb,i2):0.0;
                            size_t off=f4_off(C,f4a,a,b);
                            double va_up=hF4n[f4_off(C,f4a,a_up,b)], vb_up=hF4n[f4_off(C,f4a,a,b_up)];
                            double va_bud=(ak==0&&La>=1)||(ak==2&&ka<=S4->KT)?hF4n[f4_off(C,f4a,a_bud,b)]:0.0;
                            double vb_bud=(bk==0&&Lb>=1)||(bk==2&&kb<=S4->KT)?hF4n[f4_off(C,f4a,a,b_bud)]:0.0;
                            double v1_up=(f4_1u>=0)?hF4n[f4_off(C,f4_1u,a,b)]:0.0, v1_bud=(f4_1b>=0)?hF4n[f4_off(C,f4_1b,a,b)]:0.0;
                            double v2_up=(f4_2u>=0)?hF4n[f4_off(C,f4_2u,a,b)]:0.0, v2_bud=(f4_2b>=0)?hF4n[f4_off(C,f4_2b,a,b)]:0.0;
                            double num0=(inv_a-z_f1)*m1[b]+sgn*fabs(inv_a-z_f1)*r1[b]+(inv_b-z_f2)*m2[a]+sgn*fabs(inv_b-z_f2)*r2[a], den0=lam;
                            if(ak==0){double r=La+i1;num0+=r*va_up;den0+=r;if(ba>0){num0+=ba*va_bud;den0+=ba;}}
                            if(bk==0){double r=Lb+i2;num0+=r*vb_up;den0+=r;if(bb>0){num0+=bb*vb_bud;den0+=bb;}}
                            {double r1b=L1+(ak==0?inv_a:0.0)+i2;num0+=r1b*v1_up;den0+=r1b;if(b1>0){num0+=b1*v1_bud;den0+=b1;}}
                            {double r2b=L2+i1+(bk==0?inv_b:0.0);num0+=r2b*v2_up;den0+=r2b;if(b2>0){num0+=b2*v2_bud;den0+=b2;}}
                            double d0_a=(ak==2)?(((f4_1u>=0)?hF4n[f4_off(C,f4_1u,S4->INF_ID,b)]:0.0)-va_up):0.0;
                            double d0_b=(bk==2)?(((f4_2u>=0)?hF4n[f4_off(C,f4_2u,a,S4->INF_ID)]:0.0)-vb_up):0.0;
                            double v=num0/den0;
                            for (int ns=0;ns<12;ns++) {
                                double num=num0, den=den0;
                                if(ra_hi_e>0&&((mode==1&&va_up>v)||(mode==2&&va_up<v))){num+=ra_hi_e*va_up;den+=ra_hi_e;}
                                if(rb_hi_e>0&&((mode==1&&vb_up>v)||(mode==2&&vb_up<v))){num+=rb_hi_e*vb_up;den+=rb_hi_e;}
                                if(r1_hi>0  &&((mode==1&&v1_up>v)||(mode==2&&v1_up<v))){num+=r1_hi*v1_up;den+=r1_hi;}
                                if(r2_hi>0  &&((mode==1&&v2_up>v)||(mode==2&&v2_up<v))){num+=r2_hi*v2_up;den+=r2_hi;}
                                if(ak==2){double gr=ext_rate(mode,va_up-v,aglo,1.0);num+=gr*va_up;den+=gr;
                                    if(ka<=S4->KT){double br=ext_rate(mode,va_bud-v,tail_bmin(S4,ka),1.0);num+=br*va_bud;den+=br;}
                                    double x0=1.0/S4->D0, d1=v1_up-v;
                                    if((mode==1&&d1>d0_a)||(mode==2&&d1<d0_a)){num+=x0*v1_up;den+=x0;} else num+=x0*d0_a;}
                                if(bk==2){double gr=ext_rate(mode,vb_up-v,bglo,1.0);num+=gr*vb_up;den+=gr;
                                    if(kb<=S4->KT){double br=ext_rate(mode,vb_bud-v,tail_bmin(S4,kb),1.0);num+=br*vb_bud;den+=br;}
                                    double x0=1.0/S4->D0, d1=v2_up-v;
                                    if((mode==1&&d1>d0_b)||(mode==2&&d1<d0_b)){num+=x0*v2_up;den+=x0;} else num+=x0*d0_b;}
                                double vn=num/den; if(fabs(vn-v)<1e-13){v=vn;break;} v=vn;
                            }
                            hF4n[off]=v;
                        }
                    }
                }
            }
        }
        free(ord); free(tmp3); free(m1); free(r1); free(m2); free(r2);
    }

    /* ===== 6. SOLVE hY4 (Level 4) ===== */
    if (C->NY4_ACT > 0) {
        memset(hY4n, 0, (size_t)C->NY4_ACT * S4->npars * S4->pair_n * sizeof(double));
        int ord[256]; for (int i=0;i<C->NY4_ACT;i++) ord[i]=i;
        for (int i=0;i<C->NY4_ACT;i++) for (int j=i+1;j<C->NY4_ACT;j++) {
            int i1=ord[i],i2=ord[j],d1=S4->gD[C->y4_u[i1]],d2=S4->gD[C->y4_u[i2]],k1=S4->gK[C->y4_u[i1]],k2=S4->gK[C->y4_u[i2]];
            if (d2>d1 || (d2==d1 && k2>k1)) { ord[i]=i2; ord[j]=i1; }
        }
        double *hv3=malloc((size_t)S3->npars*S3->nchs*sizeof(double)), *hc3=malloc((size_t)S3->nchs*S3->nchs*sizeof(double));
        double *mhv=malloc((size_t)S4->npars*S4->nchs*sizeof(double)), *rhv=malloc((size_t)S4->npars*S4->nchs*sizeof(double));
        double *mhc=malloc((size_t)S4->nchs*S4->nchs*sizeof(double)), *rhc=malloc((size_t)S4->nchs*S4->nchs*sizeof(double));
        for (int oi = 0; oi < C->NY4_ACT; oi++) {
            int y4a=ord[oi], gu4=C->y4_u[y4a], Lu=S4->gL[gu4], ku=S4->gK[gu4], Du=S4->gD[gu4];
            int ui=hu_idx(C,Lu,ku), ui_up=hu_idx(C,Lu+1,ku), ci=cu_idx(C,Lu,ku), ci_up=cu_idx(C,Lu+1,ku);
            int gu2=lvl_inner(S2,Lu,ku);
            double iu=1.0/Du, bu=(double)Lu/Du, ru_hi=0.5*(ku-2), z_y1=M->rch[gu2], z_y2=M->rpar[gu2];
            int y4_up=y4_act(C,Lu+1,ku), y4_bud=(Lu>=1)?y4_act(C,Lu-1,ku+1):-1;
            for (int ip=0;ip<S3->npars;ip++) { int a3=S3->pars[ip]; for(int iw=0;iw<S3->nchs;iw++){int b3=S3->chs[iw]; hv3[(size_t)ip*S3->nchs+iw]=((ui_up>=0)?hVn[hv_off(C,a3,ui_up,b3)]:0.0)-hVn[hv_off(C,a3,ui,b3)];} }
            for (int i1=0;i1<S3->nchs;i1++) { int b1=S3->chs[i1]; for(int i2=0;i2<S3->nchs;i2++){int b2=S3->chs[i2]; hc3[(size_t)i1*S3->nchs+i2]=((ci_up>=0)?hCn[hc_off(C,ci_up,b1,b2)]:0.0)-hCn[hc_off(C,ci,b1,b2)];} }
            build_env_2d(S3,S4,1,hv3,mhv,rhv); build_env_2d(S3,S4,0,hc3,mhc,rhc);
            for (int rsum=3*S4->max_rank;rsum>=6;rsum--)
                for (int ra=S4->max_rank;ra>=2;ra--) {
                    if (!S4->p_rcnt[ra]) continue;
                    for (int r1=S4->max_rank;r1>=2;r1--) {
                        int r2=rsum-ra-r1; if (r2<r1||r2>S4->max_rank||!S4->c_rcnt[r1]||!S4->c_rcnt[r2]) continue;
                        for (int ia=0;ia<S4->p_rcnt[ra];ia++) {
                            int a=S4->p_by_r[ra][ia],ak=S4->gKind[a],La=S4->gL[a],ka=S4->gK[a],pia=S4->p_to_pi[a];
                            int a_up=(ak==0)?((S4->gD[a]+1<S4->D0)?lvl_inner(S4,La+1,ka):lvl_tail(S4,ka)):S4->INF_ID;
                            int a_bud=(ak==0&&La>=1)?lvl_inner(S4,La-1,ka+1):((ak==2&&ka<=S4->KT)?lvl_tail(S4,ka+1):a);
                            double ba=(ak==0&&La>=1)?((double)La/S4->gD[a]):0.0, ra_hi_e=(ak==0)?(0.5*ka):0.0, inv_a=S4->invD[a], fac_a=inv_a-z_y2, afac_a=fabs(fac_a);
                            double aglo=(ak==2)?tail_gmin(S4,ka,iu):0.0;
                            const double *mhv_a=mhv+(size_t)pia*S4->nchs, *rhv_a=rhv+(size_t)pia*S4->nchs;
                            for (int i1=0;i1<S4->c_rcnt[r1];i1++) {
                                int b1=S4->c_by_r[r1][i1],k1k=S4->gKind[b1],L1=S4->gL[b1],k1=S4->gK[b1],wi1=S4->w_to_wi[b1];
                                int b1_up=(k1k==0)?((S4->gD[b1]+1<S4->D0)?lvl_inner(S4,L1+1,k1):lvl_tail(S4,k1)):S4->INF_ID;
                                int b1_bud=(k1k==0&&L1>=1)?lvl_inner(S4,L1-1,k1+1):((k1k==2&&k1<=S4->KT)?lvl_tail(S4,k1+1):b1);
                                double bb1=(k1k==0&&L1>=1)?((double)L1/S4->gD[b1]):0.0, r1_hi_e=(k1k==0)?(0.5*k1):0.0, inv1=S4->invD[b1], fac1=inv1-z_y1, afac1=fabs(fac1);
                                double g1lo=(k1k==2)?tail_gmin(S4,k1,iu):0.0, m1v=mhv_a[wi1], r1v=rhv_a[wi1];
                                for (int i2=(r1==r2?i1:0);i2<S4->c_rcnt[r2];i2++) {
                                    int b2=S4->c_by_r[r2][i2],k2k=S4->gKind[b2],L2=S4->gL[b2],k2=S4->gK[b2],wi2=S4->w_to_wi[b2];
                                    int b2_up=(k2k==0)?((S4->gD[b2]+1<S4->D0)?lvl_inner(S4,L2+1,k2):lvl_tail(S4,k2)):S4->INF_ID;
                                    int b2_bud=(k2k==0&&L2>=1)?lvl_inner(S4,L2-1,k2+1):((k2k==2&&k2<=S4->KT)?lvl_tail(S4,k2+1):b2);
                                    double bb2=(k2k==0&&L2>=1)?((double)L2/S4->gD[b2]):0.0, r2_hi_e=(k2k==0)?(0.5*k2):0.0, inv2=S4->invD[b2], fac2=inv2-z_y1, afac2=fabs(fac2);
                                    double g2lo=(k2k==2)?tail_gmin(S4,k2,iu):0.0;
                                    size_t off=y4_off(C,y4a,a,b1,b2);
                                    double va_up=hY4n[y4_off(C,y4a,a_up,b1,b2)], v1_up=hY4n[y4_off(C,y4a,a,b1_up,b2)], v2_up=hY4n[y4_off(C,y4a,a,b1,b2_up)];
                                    double va_bud=(ak==0&&La>=1)||(ak==2&&ka<=S4->KT)?hY4n[y4_off(C,y4a,a_bud,b1,b2)]:0.0;
                                    double v1_bud=(k1k==0&&L1>=1)||(k1k==2&&k1<=S4->KT)?hY4n[y4_off(C,y4a,a,b1_bud,b2)]:0.0;
                                    double v2_bud=(k2k==0&&L2>=1)||(k2k==2&&k2<=S4->KT)?hY4n[y4_off(C,y4a,a,b1,b2_bud)]:0.0;
                                    double vu_up=(y4_up>=0)?hY4n[y4_off(C,y4_up,a,b1,b2)]:0.0, vu_bud=(y4_bud>=0)?hY4n[y4_off(C,y4_bud,a,b1,b2)]:0.0;
                                    double num0 = fac2*m1v+sgn*afac2*r1v + fac1*mhv_a[wi2]+sgn*afac1*rhv_a[wi2]
                                                + fac_a*mhc[(size_t)wi1*S4->nchs+wi2]+sgn*afac_a*rhc[(size_t)wi1*S4->nchs+wi2], den0=lam;
                                    if(ak==0){double r=La+iu;num0+=r*va_up;den0+=r;if(ba>0){num0+=ba*va_bud;den0+=ba;}}
                                    if(k1k==0){double r=L1+iu;num0+=r*v1_up;den0+=r;if(bb1>0){num0+=bb1*v1_bud;den0+=bb1;}}
                                    if(k2k==0){double r=L2+iu;num0+=r*v2_up;den0+=r;if(bb2>0){num0+=bb2*v2_bud;den0+=bb2;}}
                                    {double rub=Lu+(ak==0?inv_a:0.0)+(k1k==0?inv1:0.0)+(k2k==0?inv2:0.0);num0+=rub*vu_up;den0+=rub;if(bu>0){num0+=bu*vu_bud;den0+=bu;}}
                                    double d0_a=(ak==2)?(((y4_up>=0)?hY4n[y4_off(C,y4_up,S4->INF_ID,b1,b2)]:0.0)-va_up):0.0;
                                    double d0_1=(k1k==2)?(((y4_up>=0)?hY4n[y4_off(C,y4_up,a,S4->INF_ID,b2)]:0.0)-v1_up):0.0;
                                    double d0_2=(k2k==2)?(((y4_up>=0)?hY4n[y4_off(C,y4_up,a,b1,S4->INF_ID)]:0.0)-v2_up):0.0;
                                    double v=num0/den0;
                                    for (int ns=0;ns<10;ns++) {
                                        double num=num0, den=den0;
                                        if(ru_hi>0  &&((mode==1&&vu_up>v)||(mode==2&&vu_up<v))){num+=ru_hi*vu_up;den+=ru_hi;}
                                        if(ra_hi_e>0&&((mode==1&&va_up>v)||(mode==2&&va_up<v))){num+=ra_hi_e*va_up;den+=ra_hi_e;}
                                        if(r1_hi_e>0&&((mode==1&&v1_up>v)||(mode==2&&v1_up<v))){num+=r1_hi_e*v1_up;den+=r1_hi_e;}
                                        if(r2_hi_e>0&&((mode==1&&v2_up>v)||(mode==2&&v2_up<v))){num+=r2_hi_e*v2_up;den+=r2_hi_e;}
                                        if(ak==2){double gr=ext_rate(mode,va_up-v,aglo,1.0);num+=gr*va_up;den+=gr;
                                            if(ka<=S4->KT){double br=ext_rate(mode,va_bud-v,tail_bmin(S4,ka),1.0);num+=br*va_bud;den+=br;}
                                            double x0=1.0/S4->D0, d1=vu_up-v;
                                            if((mode==1&&d1>d0_a)||(mode==2&&d1<d0_a)){num+=x0*vu_up;den+=x0;} else num+=x0*d0_a;}
                                        if(k1k==2){double gr=ext_rate(mode,v1_up-v,g1lo,1.0);num+=gr*v1_up;den+=gr;
                                            if(k1<=S4->KT){double br=ext_rate(mode,v1_bud-v,tail_bmin(S4,k1),1.0);num+=br*v1_bud;den+=br;}
                                            double x0=1.0/S4->D0, d1=vu_up-v;
                                            if((mode==1&&d1>d0_1)||(mode==2&&d1<d0_1)){num+=x0*vu_up;den+=x0;} else num+=x0*d0_1;}
                                        if(k2k==2){double gr=ext_rate(mode,v2_up-v,g2lo,1.0);num+=gr*v2_up;den+=gr;
                                            if(k2<=S4->KT){double br=ext_rate(mode,v2_bud-v,tail_bmin(S4,k2),1.0);num+=br*v2_bud;den+=br;}
                                            double x0=1.0/S4->D0, d1=vu_up-v;
                                            if((mode==1&&d1>d0_2)||(mode==2&&d1<d0_2)){num+=x0*vu_up;den+=x0;} else num+=x0*d0_2;}
                                        double vn=num/den; if(fabs(vn-v)<1e-13){v=vn;break;} v=vn;
                                    }
                                    hY4n[off]=v;
                                }
                            }
                        }
                    }
                }
        }
        free(hv3); free(hc3); free(mhv); free(rhv); free(mhc); free(rhc);
    }

    /* ===== 7. SOLVE hT4 (Trident, Level 4) ===== */
    if (C->NT4_ACT > 0) {
        memset(hT4n, 0, (size_t)C->NT4_ACT * S4->tri_n * sizeof(double));
        int ord[128]; for (int i=0;i<C->NT4_ACT;i++) ord[i]=i;
        for (int i=0;i<C->NT4_ACT;i++) for (int j=i+1;j<C->NT4_ACT;j++) {
            int i1=ord[i],i2=ord[j],d1=S4->gD[C->t4_u[i1]],d2=S4->gD[C->t4_u[i2]],k1=S4->gK[C->t4_u[i1]],k2=S4->gK[C->t4_u[i2]];
            if (d2>d1 || (d2==d1 && k2>k1)) { ord[i]=i2; ord[j]=i1; }
        }
        double *hc3=malloc((size_t)S3->nchs*S3->nchs*sizeof(double));
        double *mhc=malloc((size_t)S4->nchs*S4->nchs*sizeof(double)), *rhc=malloc((size_t)S4->nchs*S4->nchs*sizeof(double));
        for (int oi = 0; oi < C->NT4_ACT; oi++) {
            int t4a=ord[oi], gu4=C->t4_u[t4a], Lu=S4->gL[gu4], ku=S4->gK[gu4], Du=S4->gD[gu4];
            int ci=cu_idx(C,Lu,ku), ci_up=cu_idx(C,Lu+1,ku), gu2=lvl_inner(S2,Lu,ku);
            double iu=1.0/Du, bu=(double)Lu/Du, ru_hi=0.5*(1+ku-3), z_ch=M->rch[gu2];
            int t4_up=t4_act(C,Lu+1,ku), t4_bud=(Lu>=1)?t4_act(C,Lu-1,ku+1):-1;
            for (int i1=0;i1<S3->nchs;i1++) { int b1=S3->chs[i1]; for(int i2=0;i2<S3->nchs;i2++){int b2=S3->chs[i2]; hc3[(size_t)i1*S3->nchs+i2]=((ci_up>=0)?hCn[hc_off(C,ci_up,b1,b2)]:0.0)-hCn[hc_off(C,ci,b1,b2)];} }
            build_env_2d(S3,S4,0,hc3,mhc,rhc);
            for (int rsum=3*S4->max_rank;rsum>=6;rsum--)
                for (int r1=S4->max_rank;r1>=2;r1--) {
                    if (!S4->c_rcnt[r1]) continue;
                    for (int r2=r1;r2<=S4->max_rank;r2++) {
                        int r3=rsum-r1-r2; if (r3<r2||r3>S4->max_rank||!S4->c_rcnt[r2]||!S4->c_rcnt[r3]) continue;
                        for (int i1=0;i1<S4->c_rcnt[r1];i1++) {
                            int b1=S4->c_by_r[r1][i1],k1k=S4->gKind[b1],L1=S4->gL[b1],k1=S4->gK[b1],wi1=S4->w_to_wi[b1];
                            int b1_up=(k1k==0)?((S4->gD[b1]+1<S4->D0)?lvl_inner(S4,L1+1,k1):lvl_tail(S4,k1)):S4->INF_ID;
                            int b1_bud=(k1k==0&&L1>=1)?lvl_inner(S4,L1-1,k1+1):((k1k==2&&k1<=S4->KT)?lvl_tail(S4,k1+1):b1);
                            double bb1=(k1k==0&&L1>=1)?((double)L1/S4->gD[b1]):0.0, r1_hi_e=(k1k==0)?(0.5*k1):0.0, inv1=S4->invD[b1], fac1=inv1-z_ch, af1=fabs(fac1);
                            double g1lo=(k1k==2)?tail_gmin(S4,k1,iu):0.0;
                            for (int i2=(r1==r2?i1:0);i2<S4->c_rcnt[r2];i2++) {
                                int b2=S4->c_by_r[r2][i2],k2k=S4->gKind[b2],L2=S4->gL[b2],k2=S4->gK[b2],wi2=S4->w_to_wi[b2];
                                int b2_up=(k2k==0)?((S4->gD[b2]+1<S4->D0)?lvl_inner(S4,L2+1,k2):lvl_tail(S4,k2)):S4->INF_ID;
                                int b2_bud=(k2k==0&&L2>=1)?lvl_inner(S4,L2-1,k2+1):((k2k==2&&k2<=S4->KT)?lvl_tail(S4,k2+1):b2);
                                double bb2=(k2k==0&&L2>=1)?((double)L2/S4->gD[b2]):0.0, r2_hi_e=(k2k==0)?(0.5*k2):0.0, inv2=S4->invD[b2], fac2=inv2-z_ch, af2=fabs(fac2);
                                double g2lo=(k2k==2)?tail_gmin(S4,k2,iu):0.0;
                                double m12=mhc[(size_t)wi1*S4->nchs+wi2], r12=rhc[(size_t)wi1*S4->nchs+wi2];
                                for (int i3=(r2==r3?i2:0);i3<S4->c_rcnt[r3];i3++) {
                                    int b3=S4->c_by_r[r3][i3],k3k=S4->gKind[b3],L3=S4->gL[b3],k3=S4->gK[b3],wi3=S4->w_to_wi[b3];
                                    int b3_up=(k3k==0)?((S4->gD[b3]+1<S4->D0)?lvl_inner(S4,L3+1,k3):lvl_tail(S4,k3)):S4->INF_ID;
                                    int b3_bud=(k3k==0&&L3>=1)?lvl_inner(S4,L3-1,k3+1):((k3k==2&&k3<=S4->KT)?lvl_tail(S4,k3+1):b3);
                                    double bb3=(k3k==0&&L3>=1)?((double)L3/S4->gD[b3]):0.0, r3_hi_e=(k3k==0)?(0.5*k3):0.0, inv3=S4->invD[b3], fac3=inv3-z_ch, af3=fabs(fac3);
                                    double g3lo=(k3k==2)?tail_gmin(S4,k3,iu):0.0;
                                    size_t off=t4_off(C,t4a,b1,b2,b3);
                                    double v1_up=hT4n[t4_off(C,t4a,b1_up,b2,b3)], v2_up=hT4n[t4_off(C,t4a,b1,b2_up,b3)], v3_up=hT4n[t4_off(C,t4a,b1,b2,b3_up)];
                                    double v1_bud=(k1k==0&&L1>=1)||(k1k==2&&k1<=S4->KT)?hT4n[t4_off(C,t4a,b1_bud,b2,b3)]:0.0;
                                    double v2_bud=(k2k==0&&L2>=1)||(k2k==2&&k2<=S4->KT)?hT4n[t4_off(C,t4a,b1,b2_bud,b3)]:0.0;
                                    double v3_bud=(k3k==0&&L3>=1)||(k3k==2&&k3<=S4->KT)?hT4n[t4_off(C,t4a,b1,b2,b3_bud)]:0.0;
                                    double vu_up=(t4_up>=0)?hT4n[t4_off(C,t4_up,b1,b2,b3)]:0.0, vu_bud=(t4_bud>=0)?hT4n[t4_off(C,t4_bud,b1,b2,b3)]:0.0;
                                    double num0 = fac3*m12+sgn*af3*r12 + fac2*mhc[(size_t)wi1*S4->nchs+wi3]+sgn*af2*rhc[(size_t)wi1*S4->nchs+wi3]
                                                + fac1*mhc[(size_t)wi2*S4->nchs+wi3]+sgn*af1*rhc[(size_t)wi2*S4->nchs+wi3], den0=lam;
                                    if(k1k==0){double r=L1+iu;num0+=r*v1_up;den0+=r;if(bb1>0){num0+=bb1*v1_bud;den0+=bb1;}}
                                    if(k2k==0){double r=L2+iu;num0+=r*v2_up;den0+=r;if(bb2>0){num0+=bb2*v2_bud;den0+=bb2;}}
                                    if(k3k==0){double r=L3+iu;num0+=r*v3_up;den0+=r;if(bb3>0){num0+=bb3*v3_bud;den0+=bb3;}}
                                    {double rub=Lu+(k1k==0?inv1:0.0)+(k2k==0?inv2:0.0)+(k3k==0?inv3:0.0);num0+=rub*vu_up;den0+=rub;if(bu>0){num0+=bu*vu_bud;den0+=bu;}}
                                    double d0_1=(k1k==2)?(((t4_up>=0)?hT4n[t4_off(C,t4_up,S4->INF_ID,b2,b3)]:0.0)-v1_up):0.0;
                                    double d0_2=(k2k==2)?(((t4_up>=0)?hT4n[t4_off(C,t4_up,b1,S4->INF_ID,b3)]:0.0)-v2_up):0.0;
                                    double d0_3=(k3k==2)?(((t4_up>=0)?hT4n[t4_off(C,t4_up,b1,b2,S4->INF_ID)]:0.0)-v3_up):0.0;
                                    double v=num0/den0;
                                    for (int ns=0;ns<10;ns++) {
                                        double num=num0, den=den0;
                                        if(ru_hi>0  &&((mode==1&&vu_up>v)||(mode==2&&vu_up<v))){num+=ru_hi*vu_up;den+=ru_hi;}
                                        if(r1_hi_e>0&&((mode==1&&v1_up>v)||(mode==2&&v1_up<v))){num+=r1_hi_e*v1_up;den+=r1_hi_e;}
                                        if(r2_hi_e>0&&((mode==1&&v2_up>v)||(mode==2&&v2_up<v))){num+=r2_hi_e*v2_up;den+=r2_hi_e;}
                                        if(r3_hi_e>0&&((mode==1&&v3_up>v)||(mode==2&&v3_up<v))){num+=r3_hi_e*v3_up;den+=r3_hi_e;}
                                        if(k1k==2){double gr=ext_rate(mode,v1_up-v,g1lo,1.0);num+=gr*v1_up;den+=gr;
                                            if(k1<=S4->KT){double br=ext_rate(mode,v1_bud-v,tail_bmin(S4,k1),1.0);num+=br*v1_bud;den+=br;}
                                            double x0=1.0/S4->D0, d1=vu_up-v;
                                            if((mode==1&&d1>d0_1)||(mode==2&&d1<d0_1)){num+=x0*vu_up;den+=x0;} else num+=x0*d0_1;}
                                        if(k2k==2){double gr=ext_rate(mode,v2_up-v,g2lo,1.0);num+=gr*v2_up;den+=gr;
                                            if(k2<=S4->KT){double br=ext_rate(mode,v2_bud-v,tail_bmin(S4,k2),1.0);num+=br*v2_bud;den+=br;}
                                            double x0=1.0/S4->D0, d1=vu_up-v;
                                            if((mode==1&&d1>d0_2)||(mode==2&&d1<d0_2)){num+=x0*vu_up;den+=x0;} else num+=x0*d0_2;}
                                        if(k3k==2){double gr=ext_rate(mode,v3_up-v,g3lo,1.0);num+=gr*v3_up;den+=gr;
                                            if(k3<=S4->KT){double br=ext_rate(mode,v3_bud-v,tail_bmin(S4,k3),1.0);num+=br*v3_bud;den+=br;}
                                            double x0=1.0/S4->D0, d1=vu_up-v;
                                            if((mode==1&&d1>d0_3)||(mode==2&&d1<d0_3)){num+=x0*vu_up;den+=x0;} else num+=x0*d0_3;}
                                        double vn=num/den; if(fabs(vn-v)<1e-13){v=vn;break;} v=vn;
                                    }
                                    hT4n[off]=v;
                                }
                            }
                        }
                    }
                }
        }
        free(hc3); free(mhc); free(rhc);
    }
}

/* Exact analytic f tail for M_g = 1 (D0_g < DF): g(s_bud -> nb) = g_inf + (D0_g/D)*(g[T1,nb] - g_inf)! */
static void update_analytic_f_tail(Model *M, double lam) {
    const Layout *C = &M->C; const Level *S2 = &C->L2; int N2 = S2->NG, DF = C->DF, D0g = S2->D0;
    int nb2 = lvl_inner(S2, 1, 0), v_nb = vidx(C, 1, 0);
    double fnb = M->f[v_nb], g_inf = M->g[(size_t)S2->INF_ID*N2 + nb2], g_min = g_inf, g_max = g_inf;
    for (int kap = 1; kap <= S2->KT + 1; kap++) {
        double gv = M->g[(size_t)lvl_tail(S2, kap)*N2 + nb2];
        if (gv < g_min) g_min = gv;
        if (gv > g_max) g_max = gv;
    }
    if (M->mode == 1) {
        double B_est = (fnb + g_inf) / (1.0 + lam);
        double corr = fmax(0.0, (D0g * (g_max - g_inf) - B_est) / DF);
        double Cnb = fnb + g_inf + corr;
        M->B_up = Cnb / (1.0 + lam);
        M->A_up = Cnb / (lam * (1.0 + lam));
    } else {
        double Cinf = fnb + g_inf;
        double B1_0 = Cinf * (2.0 + lam) / (lam * (1.0 + lam));
        double B2_0 = Cinf / (lam * (1.0 + lam));
        double d_g = D0g * (g_min - g_inf), cF = (double)DF / (DF + 1.0);
        double lo = 0.0, hi = 0.5 * B2_0, best_eta = hi, best_cd = 0.0;
        for (int bi = 0; bi < 50; bi++) {
            double eta = (bi == 0) ? 0.0 : 0.5 * (lo + hi);
            double b1 = B1_0 - eta, b2 = B2_0 - eta;
            double cd = ((b1 - 2.0*b2) + b2/DF - d_g) / (cF + lam) + 1e-4;
            if (cd < 0.0) cd = 0.0;
            /* Check P0 Bernstein >= 0 AND P0 + P1/DF Bernstein >= 0 AND dLf >= 0! */
            double q1 = b1 + Cinf - lam*b1, q2 = -(b1+2.0*b2) + lam*b2, q3 = 2.0*b2;
            double invDF = 1.0 / DF;
            double p1 = q1 + invDF*(-b1-b2-cd+lam*cd-b1-fmax(0.0,b2-cd)/DF+d_g);
            double p2 = q2 + invDF*((b1+2.0*b2)+2.0*b2+cd*(cF+1.0)+2.0*b2);
            double p3 = q3 + invDF*(-3.0*b2);
            int ok0 = (q1 >= 0 && 2.0*q1+q2 >= 0 && q1+q2+q3 >= -1e-14);
            int ok1 = (p1 >= 0 && 2.0*p1+p2 >= 0 && p1+p2+p3 >= 0);
            int okd = (b1 - 2.0*b2 - (b2+cd)/(DF+1.0) > 0);
            int ok = ok0 && ok1 && okd;
            if (bi == 0 && ok) { best_eta = 0.0; best_cd = cd; break; }
            if (ok) { best_eta = eta; best_cd = cd; hi = eta; } else lo = eta;
        }
        M->B1_lo = B1_0 - best_eta; M->B2_lo = B2_0 - best_eta; M->CD_lo = best_cd; M->eta_lo = best_eta;
    }
}

/* Brute-force tree evaluator & verifier for ALL 8 motif families + 3 scales */
static double eval_phi_h1(const Model *M, int n, const int *par, const int *deg) {
    const Layout *C = &M->C; const Level *S2 = &C->L2, *S3 = &C->L3, *S4 = &C->L4;
    int *Lc=calloc(n,sizeof(int)), *Kc=calloc(n,sizeof(int)), *c_off=calloc(n+1,sizeof(int));
    for (int v=1;v<n;v++) if (par[v]>=1) { if (deg[v]==1) Lc[par[v]]++; else Kc[par[v]]++; }
    for (int v=1;v<n;v++) c_off[v+1]=c_off[v]+Kc[v];
    int *cch=malloc((c_off[n]+1)*sizeof(int)), *pos=calloc(n,sizeof(int));
    for (int v=1;v<n;v++) if (par[v]>=1 && deg[v]>=2) cch[c_off[par[v]]+pos[par[v]]++]=v;
    double s = 0.0;
    for (int u=1;u<n;u++) if (deg[u]>=2) s += eval_f(M, Lc[u], Kc[u]);
    for (int w=1;w<n;w++) if (deg[w]>=2) {
        int p = par[w], pid[2], wid[2], np, nw; double pw[2], ww[2];
        if (p == 0) { pid[0]=S2->INF_ID; pw[0]=1.0; np=1; } else get_mu(S2, Lc[p], Kc[p], pid, pw, &np);
        get_mu(S2, Lc[w], Kc[w], wid, ww, &nw);
        for (int a=0;a<np;a++) for (int b=0;b<nw;b++) s += pw[a]*ww[b]*M->g[(size_t)pid[a]*S2->NG+wid[b]];
        if (p >= 1) {
            int ui = hu_idx(C, Lc[p], Kc[p]);
            if (ui >= 0) {
                int gp = par[p];
                if (gp == 0) { pid[0]=S3->INF_ID; pw[0]=1.0; np=1; } else get_mu(S3, Lc[gp], Kc[gp], pid, pw, &np);
                get_mu(S3, Lc[w], Kc[w], wid, ww, &nw);
                for (int a=0;a<np;a++) for (int b=0;b<nw;b++) s += pw[a]*ww[b]*M->hV[hv_off(C,pid[a],ui,wid[b])];
            }
            int gp = par[p];
            if (gp >= 1) {
                int p4a = p4_act(C, Lc[gp], Kc[gp], Lc[p], Kc[p]);
                if (p4a >= 0) {
                    int ggp = par[gp];
                    if (ggp == 0) { pid[0]=S4->INF_ID; pw[0]=1.0; np=1; } else get_mu(S4, Lc[ggp], Kc[ggp], pid, pw, &np);
                    get_mu(S4, Lc[w], Kc[w], wid, ww, &nw);
                    for (int a=0;a<np;a++) for (int b=0;b<nw;b++) s += pw[a]*ww[b]*M->hP4[p4_off(C,p4a,pid[a],wid[b])];
                }
            }
        }
    }
    for (int u=1;u<n;u++) if (deg[u]>=2 && Kc[u]>=2) {
        int pu = par[u], ci = cu_idx(C, Lc[u], Kc[u]), y4a = y4_act(C, Lc[u], Kc[u]), t4a = t4_act(C, Lc[u], Kc[u]);
        for (int i=c_off[u];i<c_off[u+1];i++) {
            int w1 = cch[i];
            for (int j=i+1;j<c_off[u+1];j++) {
                int w2 = cch[j];
                if (ci >= 0) {
                    int a1[2],a2[2],n1,n2; double w1v[2],w2v[2];
                    get_mu(S3,Lc[w1],Kc[w1],a1,w1v,&n1); get_mu(S3,Lc[w2],Kc[w2],a2,w2v,&n2);
                    for (int x=0;x<n1;x++) for (int y=0;y<n2;y++) s += w1v[x]*w2v[y]*M->hC[hc_off(C,ci,a1[x],a2[y])];
                }
                if (y4a >= 0) {
                    int pa[2],a1[2],a2[2],np,n1,n2; double pw[2],w1v[2],w2v[2];
                    if (pu == 0) { pa[0]=S4->INF_ID; pw[0]=1.0; np=1; } else get_mu(S4,Lc[pu],Kc[pu],pa,pw,&np);
                    get_mu(S4,Lc[w1],Kc[w1],a1,w1v,&n1); get_mu(S4,Lc[w2],Kc[w2],a2,w2v,&n2);
                    for (int ip=0;ip<np;ip++) for (int x=0;x<n1;x++) for (int y=0;y<n2;y++)
                        s += pw[ip]*w1v[x]*w2v[y]*M->hY4[y4_off(C,y4a,pa[ip],a1[x],a2[y])];
                }
                int f4a1 = f4_act(C, Lc[u], Kc[u], Lc[w1], Kc[w1]);
                if (f4a1 >= 0) {
                    int aa[2], na; double aw[2]; get_mu(S4, Lc[w2], Kc[w2], aa, aw, &na);
                    for (int ic=c_off[w1]; ic<c_off[w1+1]; ic++) {
                        int b = cch[ic], bb[2], nb; double bw[2]; get_mu(S4, Lc[b], Kc[b], bb, bw, &nb);
                        for (int x=0;x<na;x++) for (int y=0;y<nb;y++) s += aw[x]*bw[y]*M->hF4[f4_off(C,f4a1,aa[x],bb[y])];
                    }
                }
                int f4a2 = f4_act(C, Lc[u], Kc[u], Lc[w2], Kc[w2]);
                if (f4a2 >= 0) {
                    int aa[2], na; double aw[2]; get_mu(S4, Lc[w1], Kc[w1], aa, aw, &na);
                    for (int ic=c_off[w2]; ic<c_off[w2+1]; ic++) {
                        int b = cch[ic], bb[2], nb; double bw[2]; get_mu(S4, Lc[b], Kc[b], bb, bw, &nb);
                        for (int x=0;x<na;x++) for (int y=0;y<nb;y++) s += aw[x]*bw[y]*M->hF4[f4_off(C,f4a2,aa[x],bb[y])];
                    }
                }
                if (t4a >= 0) {
                    int a1[2],a2[2],n1,n2; double w1v[2],w2v[2];
                    get_mu(S4,Lc[w1],Kc[w1],a1,w1v,&n1); get_mu(S4,Lc[w2],Kc[w2],a2,w2v,&n2);
                    for (int m=j+1; m<c_off[u+1]; m++) {
                        int w3 = cch[m], a3[2], n3; double w3v[2]; get_mu(S4,Lc[w3],Kc[w3],a3,w3v,&n3);
                        for (int x=0;x<n1;x++) for (int y=0;y<n2;y++) for (int z=0;z<n3;z++)
                            s += w1v[x]*w2v[y]*w3v[z]*M->hT4[t4_off(C,t4a,a1[x],a2[y],a3[z])];
                    }
                }
            }
        }
    }
    free(Lc); free(Kc); free(c_off); free(cch); free(pos);
    return s;
}

static int verify_h1_trees_bf(const Model *M) {
    int fails = 0, tested = 0; double worst = 1e30;
    double sgn = (M->mode == 1) ? 1.0 : -1.0;
    uint64_t seed = 9876543ULL; for (int i=0;i<4;i++) xs[i]=sm(&seed);
    for (int t = 0; t < 32; t++) {
        int n = 0; int *par = NULL, *deg = NULL;
        if (t < 22) {
            n = 4 + t * 5;
            par = calloc(n+1, sizeof(int)); deg = calloc(n+1, sizeof(int));
            int kids[140][140], kcnt[140] = {0};
            par[0]=-1; deg[0]=1000000; par[1]=0; deg[1]=2; par[2]=1; deg[2]=1;
            kids[1][0]=2; kcnt[1]=1;
            for (int cur=3; cur<n; cur++) {
                int V = 1 + (int)unif(cur - 1), W;
                if (deg[V] == 1) W = par[V];
                else { int r=(int)unif(deg[V]-(par[V]==0?1:0)); W=(par[V]!=0)?((r==0)?par[V]:kids[V][r-1]):kids[V][r]; }
                par[cur]=W; deg[cur]=1; deg[W]++; kids[W][kcnt[W]++]=cur;
            }
        } else if (t == 22) {
            int p_arr[] = {-1, 0, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 2, 14, 15};
            n = 18; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            for (int i=0;i<n;i++) par[i]=p_arr[i];
            deg[1]=1; for (int i=2;i<n;i++) { deg[i]++; deg[par[i]]++; }
        } else if (t == 23) {
            n = 70; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            int c1=idx++, c2=idx++, c3=idx++; par[c1]=1; par[c2]=1; par[c3]=1; deg[1]+=3; deg[c1]=deg[c2]=deg[c3]=1;
            par[idx]=1; deg[1]++; deg[idx++]=1;
            for(int i=0;i<16;i++){par[idx]=c1;deg[c1]++;deg[idx++]=1;}
            for(int i=0;i<30;i++){par[idx]=c2;deg[c2]++;deg[idx++]=1;}
            par[idx]=c3; deg[c3]++; deg[idx++]=1; n = idx;
        } else if (t == 24) {
            n = 180; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            int c1=idx++, c2=idx++, c3=idx++; par[c1]=1; par[c2]=1; par[c3]=1; deg[1]+=3; deg[c1]=deg[c2]=deg[c3]=1;
            for(int i=0;i<75;i++){par[idx]=1;deg[1]++;deg[idx++]=1;}
            for(int i=0;i<30;i++){par[idx]=c1;deg[c1]++;deg[idx++]=1;}
            for(int i=0;i<15;i++){par[idx]=c2;deg[c2]++;deg[idx++]=1;}
            par[idx]=c3; deg[c3]++; deg[idx++]=1; n = idx;
        } else if (t == 25) {
            int L = 250; n = 20 + L; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            for(int i=0;i<L;i++){par[idx]=1;deg[1]++;deg[idx++]=1;}
            int c1=idx++; par[c1]=1; deg[1]++; deg[c1]=1;
            int c2=idx++; par[c2]=c1; deg[c1]++; deg[c2]=1;
            int c3=idx++; par[c3]=c2; deg[c2]++; deg[c3]=1;
            par[idx]=c3; deg[c3]++; deg[idx++]=1; par[idx]=c1; deg[c1]++; deg[idx++]=1; n = idx;
        } else if (t == 26) {
            int L = M->C.DF + 30; n = 2 + L; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; for(int i=2;i<n;i++){par[i]=1;deg[1]++;deg[i]=1;}
        } else if (t == 27) {
            int nc = 30, nl = M->C.DF - 20; n = 2 + 2*nc + nl;
            par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            for(int i=0;i<nc;i++){int c=idx++; par[c]=1; deg[1]++; deg[c]=1; par[idx]=c; deg[c]++; deg[idx++]=1;}
            for(int i=0;i<nl;i++){par[idx]=1;deg[1]++;deg[idx++]=1;} n = idx;
        } else if (t == 28) {
            n = 170; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            int c1=idx++; par[c1]=1; deg[1]++; deg[c1]=1;
            for(int i=0;i<78;i++){par[idx]=1;deg[1]++;deg[idx++]=1;}
            for(int i=0;i<78;i++){par[idx]=c1;deg[c1]++;deg[idx++]=1;} n = idx;
        } else if (t == 29) {
            n = 160; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            int c1=idx++; par[c1]=1; deg[1]++; deg[c1]=1;
            for(int i=0;i<M->C.L2.D0-3;i++){par[idx]=1;deg[1]++;deg[idx++]=1;}
            for(int i=0;i<M->C.L2.D0-2;i++){par[idx]=c1;deg[c1]++;deg[idx++]=1;} n = idx;
        } else if (t == 30) {
            int nc = 28; n = 2 + 2*nc; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            for(int i=0;i<nc;i++){int c=idx++; par[c]=1; deg[1]++; deg[c]=1; par[idx]=c; deg[c]++; deg[idx++]=1;} n = idx;
        } else {
            n = 100; par=calloc(n+1,sizeof(int)); deg=calloc(n+1,sizeof(int));
            par[0]=-1; par[1]=0; deg[1]=1; int idx=2;
            for(int i=0;i<32;i++){par[idx]=1;deg[1]++;deg[idx++]=1;}
            int u1=idx++; par[u1]=1; deg[1]++; deg[u1]=1;
            int u2=idx++; par[u2]=u1; deg[u1]++; deg[u2]=1;
            int u3=idx++; par[u3]=u1; deg[u1]++; deg[u3]=1;
            int b1=idx++; par[b1]=u2; deg[u2]++; deg[b1]=1;
            for(int i=0;i<16;i++){par[idx]=b1;deg[b1]++;deg[idx++]=1;}
            for(int i=0;i<16;i++){par[idx]=u3;deg[u3]++;deg[idx++]=1;}
            par[idx]=u2; deg[u2]++; deg[idx++]=1; n = idx;
        }
        double phi0 = eval_phi_h1(M, n, par, deg);
        int *acnt=calloc(n,sizeof(int)), *aoff=calloc(n+1,sizeof(int));
        for (int v=1;v<n;v++) if (par[v]>=1) { acnt[v]++; acnt[par[v]]++; }
        for (int v=0;v<n;v++) aoff[v+1]=aoff[v]+acnt[v];
        int *adj=malloc((aoff[n]+1)*sizeof(int)), *pos=calloc(n,sizeof(int));
        for (int v=1;v<n;v++) if (par[v]>=1) { adj[aoff[v]+pos[v]++]=par[v]; adj[aoff[par[v]]+pos[par[v]]++]=v; }
        double drift = 0.0;
        for (int W = 1; W < n; W++) {
            double wt = 0.0;
            for (int j=aoff[W]; j<aoff[W+1]; j++) wt += 1.0 / deg[adj[j]];
            par[n] = W; deg[n] = 1; deg[W]++;
            double phi1 = eval_phi_h1(M, n+1, par, deg);
            deg[W]--;
            drift += wt * (phi1 - phi0);
        }
        double sl = sgn * (M->lam * phi0 - drift) / phi0;
        if (sl < worst) worst = sl;
        if (sl < -1e-7) {
            fails++;
            fprintf(stderr, "  H1 BF FAIL t=%d (n=%d): drift/phi=%.9f vs lam=%.9f (sl=%.3e)\n", t, n, drift/phi0, M->lam, sl);
        }
        tested++;
        free(par); free(deg); free(acnt); free(aoff); free(adj); free(pos);
    }
    printf("  H1 BF tree check (%d trees): fails=%d worst_rel_slack=%.3e\n", tested, fails, worst);
    return fails;
}

static void run_empirical(Model *M, int NMAX, int RUNS) {
    const Layout *C = &M->C; const Level *S2 = &C->L2; int N2 = S2->NG;
    double *cp=calloc(N2,sizeof(double)), *cc=calloc(N2,sizeof(double)), *rp=calloc(N2,sizeof(double)), *rc=calloc(N2,sizeof(double));
    uint64_t seed = 20261124ULL; for (int i=0;i<4;i++) xs[i]=sm(&seed);
    uint32_t *par=malloc(NMAX*sizeof(uint32_t)), *deg=malloc(NMAX*sizeof(uint32_t));
    uint32_t *kcnt=calloc(NMAX,sizeof(uint32_t)), *kcap=calloc(NMAX,sizeof(uint32_t));
    uint32_t **kids=calloc(NMAX,sizeof(uint32_t*)); int *st=malloc(NMAX*sizeof(int));
    for (int run=0; run<RUNS; run++) {
        for (int i=0;i<NMAX;i++) { kcnt[i]=0; if(kcap[i]>64){free(kids[i]);kids[i]=NULL;kcap[i]=0;} }
        par[0]=UINT32_MAX; deg[0]=2; par[1]=0; deg[1]=1; par[2]=0; deg[2]=1;
        if (kcap[0]<2){kids[0]=realloc(kids[0],2*sizeof(uint32_t));kcap[0]=2;}
        kids[0][0]=1; kids[0][1]=2; kcnt[0]=2;
        for (uint32_t n=3; n<(uint32_t)NMAX; n++) {
            uint32_t V=unif(n), d=deg[V], r=unif(d);
            uint32_t W=(par[V]!=UINT32_MAX)?((r==0)?par[V]:kids[V][r-1]):kids[V][r];
            par[n]=W; deg[n]=1; deg[W]++;
            if (kcnt[W]==kcap[W]){uint32_t nc=kcap[W]?2*kcap[W]:2;kids[W]=realloc(kids[W],nc*sizeof(uint32_t));kcap[W]=nc;}
            kids[W][kcnt[W]++]=n;
        }
        int *Lc=calloc(NMAX,sizeof(int)), *Kc=calloc(NMAX,sizeof(int));
        for (int u=0;u<NMAX;u++) if (deg[u]>=2) {
            for (uint32_t i=0;i<kcnt[u];i++){if (deg[kids[u][i]]==1) Lc[u]++; else Kc[u]++;}
            int D = (int)deg[u];
            if (D < S2->D0 && par[u]!=UINT32_MAX) st[u] = lvl_inner(S2, Lc[u], Kc[u]);
            else if (D >= S2->D0) st[u] = lvl_tail(S2, Kc[u]);
            else st[u] = -1;
        } else st[u] = -1;
        for (int w=1;w<NMAX;w++) if (deg[w]>=2 && st[w]>=0 && st[par[w]]>=0) {
            int p=par[w], sp=st[p], sw=st[w];
            M->emp_edge[(size_t)sp*N2+sw] += 1.0; rp[sw] += 1.0/deg[p]; cp[sw] += 1.0; rc[sp] += 1.0/deg[w]; cc[sp] += 1.0;
        }
        free(Lc); free(Kc);
    }
    for (int u=0; u<S2->INNER_N; u++) {
        int k=S2->gK[u], L=S2->gL[u];
        M->rpar[u] = (cp[u]>=5.0)?(rp[u]/cp[u]):(0.25/(1+0.15*L+0.35*k));
        M->rch[u]  = (cc[u]>=5.0)?(rc[u]/cc[u]):(k>0?(0.36/(1+0.04*L+0.08*k)):0.0);
        if (hu_idx(C,L,k)>=0) { M->xpar[u]=M->rpar[u]; M->xch[u]=M->rch[u]; }
        if (cu_idx(C,L,k)>=0) { M->xcc[u]=M->rch[u]; }
    }
    for (int i=0;i<NMAX;i++) free(kids[i]);
    free(par); free(deg); free(kcnt); free(kcap); free(kids); free(st); free(cp); free(cc); free(rp); free(rc);
}

static double eval_lambda(Model *M, double lam, int max_it, int opt_k_iters) {
    const Layout *C = &M->C; const Level *S2 = &C->L2, *S3 = &C->L3, *S4 = &C->L4;
    int N2 = S2->NG, DF = C->DF, mode = M->mode; double sgn = (mode == 1) ? 1.0 : -1.0;
    int v_nb = vidx(C, 1, 0), nb2 = lvl_inner(S2, 1, 0);
    size_t g_len  = (size_t)N2 * N2;
    size_t hv_len = (size_t)S3->npars * C->NHU_ACT * S3->nchs;
    size_t hc_len = (size_t)C->NCU_ACT * S3->pair_n;
    size_t p4_len = (size_t)C->NP4_ACT * S4->npars * S4->nchs;
    size_t f4_len = (size_t)C->NF4_ACT * S4->nchs * S4->nchs;
    size_t y4_len = (size_t)C->NY4_ACT * S4->npars * S4->pair_n;
    size_t t4_len = (size_t)C->NT4_ACT * S4->tri_n;

    double *fn=calloc(C->NV,sizeof(double)), *gn=calloc(g_len,sizeof(double));
    double *hVn=calloc(hv_len,sizeof(double)), *hCn=calloc(hc_len,sizeof(double));
    double *hP4n=calloc(p4_len+1,sizeof(double)), *hF4n=calloc(f4_len+1,sizeof(double));
    double *hY4n=calloc(y4_len+1,sizeof(double)), *hT4n=calloc(t4_len+1,sizeof(double));
    double *dLf_g=calloc(N2,sizeof(double)), *A_g=calloc(g_len,sizeof(double));

    double rho = 1.0;
    for (int it = 0; it < max_it; it++) {
        update_analytic_f_tail(M, lam);

        for (int D = DF - 1; D >= 2; D--) {
            for (int k = D - 1; k >= 0; k--) {
                int L = D - 1 - k, u = vidx(C, L, k);
                double f_up = (D + 1 < DF) ? fn[vidx(C, L+1, k)] : eval_f(M, L+1, k);
                double bL = (double)L / D, bv = 0.0;
                if (L >= 1) {
                    int a1[2], n1; double w1[2]; get_mu(S2, L-1, k+1, a1, w1, &n1);
                    double gnb = 0.0; for (int ii=0;ii<n1;ii++) gnb += w1[ii]*M->g[(size_t)a1[ii]*N2+nb2];
                    bv = fn[vidx(C, L-1, k+1)] + M->f[v_nb] + gnb;
                }
                if (D < S2->D0) {
                    int gu = lvl_inner(S2, L, k);
                    double up = L + M->rpar[gu] + k * M->rch[gu];
                    fn[u] = (up * f_up + (L >= 1 ? bL * bv : 0.0)) / (lam + up + (L >= 1 ? bL : 0.0));
                } else {
                    double r0 = L, r1 = L + 0.5 * (k + 1);
                    double v0 = (r0 * f_up + (L >= 1 ? bL * bv : 0.0)) / (lam + r0 + (L >= 1 ? bL : 0.0));
                    double v1 = (r1 * f_up + (L >= 1 ? bL * bv : 0.0)) / (lam + r1 + (L >= 1 ? bL : 0.0));
                    fn[u] = (mode == 1) ? fmax(v0, v1) : fmin(v0, v1);
                }
            }
        }

        memset(dLf_g, 0, N2 * sizeof(double));
        for (int u = 0; u < S2->INNER_N; u++) dLf_g[u] = fn[vidx(C, S2->gL[u]+1, S2->gK[u])] - fn[vidx(C, S2->gL[u], S2->gK[u])];

        int do_opt = (it >= 1 && it <= opt_k_iters);
        memset(A_g, 0, g_len * sizeof(double));
        double *xs_b=malloc(N2*sizeof(double)), *wxs_b=malloc(N2*sizeof(double));
        double *ys_b=malloc(N2*sizeof(double)), *wys_b=malloc(N2*sizeof(double));
        for (int u = 0; u < S2->INNER_N; u++) {
            int Lu=S2->gL[u], ku=S2->gK[u], up=(S2->gD[u]+1<S2->D0)?lvl_inner(S2,Lu+1,ku):lvl_tail(S2,ku);
            int in_hu = (hu_idx(C,Lu,ku)>=0), in_cu = (cu_idx(C,Lu,ku)>=0);
            if (!in_hu && ku >= 1) {
                if (do_opt) {
                    int nx=0, ny=0;
                    for (int ip=0;ip<S2->npars;ip++){int p=S2->pars[ip]; xs_b[nx]=S2->invD[p]; wxs_b[nx]=M->emp_edge[(size_t)p*N2+u]+1e-5/(S2->gKind[p]==0?(double)S2->gD[p]*S2->gD[p]:1e4); nx++;}
                    for (int iw=0;iw<S2->nchs;iw++){int w=S2->chs[iw]; ys_b[ny]=M->g[(size_t)up*N2+w]-M->g[(size_t)u*N2+w]; wys_b[ny]=M->emp_edge[(size_t)u*N2+w]+1e-5/(S2->gKind[w]==0?(double)S2->gD[w]*S2->gD[w]:1e4); ny++;}
                    opt_knots(M->KNOTS, M->k_pc[u], sgn, xs_b, wxs_b, nx, ys_b, wys_b, ny);
                }
                for (int ip=0;ip<S2->npars;ip++){int p=S2->pars[ip]; A_g[(size_t)p*N2+u] += ku * eval_pl(M->KNOTS, M->k_pc[u], S2->invD[p]);}
                if (do_opt) {
                    int nx=0, ny=0;
                    for (int iw=0;iw<S2->nchs;iw++){int w=S2->chs[iw]; xs_b[nx]=S2->invD[w]; wxs_b[nx]=M->emp_edge[(size_t)u*N2+w]+1e-5/(S2->gKind[w]==0?(double)S2->gD[w]*S2->gD[w]:1e4); nx++;}
                    for (int ip=0;ip<S2->npars;ip++){int p=S2->pars[ip]; ys_b[ny]=M->g[(size_t)p*N2+up]-M->g[(size_t)p*N2+u]; wys_b[ny]=M->emp_edge[(size_t)p*N2+u]+1e-5/(S2->gKind[p]==0?(double)S2->gD[p]*S2->gD[p]:1e4); ny++;}
                    opt_knots(M->KNOTS, M->k_cp[u], sgn, xs_b, wxs_b, nx, ys_b, wys_b, ny);
                }
                for (int iw=0;iw<S2->nchs;iw++){int w=S2->chs[iw]; A_g[(size_t)u*N2+w] += eval_pl(M->KNOTS, M->k_cp[u], S2->invD[w]);}
            }
            if (!in_cu && ku >= 2) {
                if (do_opt) {
                    int nx=0, ny=0;
                    for (int iw=0;iw<S2->nchs;iw++){int w=S2->chs[iw]; xs_b[nx]=S2->invD[w]; wxs_b[nx]=M->emp_edge[(size_t)u*N2+w]+1e-5/(S2->gKind[w]==0?(double)S2->gD[w]*S2->gD[w]:1e4); nx++; ys_b[ny]=M->g[(size_t)up*N2+w]-M->g[(size_t)u*N2+w]; wys_b[ny]=wxs_b[nx-1]; ny++;}
                    opt_knots(M->KNOTS, M->k_cc[u], sgn, xs_b, wxs_b, nx, ys_b, wys_b, ny);
                }
                for (int iw=0;iw<S2->nchs;iw++){int w=S2->chs[iw]; A_g[(size_t)u*N2+w] += (ku - 1) * eval_pl(M->KNOTS, M->k_cc[u], S2->invD[w]);}
            }
        }
        free(xs_b); free(wxs_b); free(ys_b); free(wys_b);

        solve_all_motifs(M, lam, dLf_g, A_g, M->hV, M->hC, M->hP4, M->hF4, M->hY4, M->hT4,
                         gn, hVn, hCn, hP4n, hF4n, hY4n, hT4n);

        rho = fn[v_nb] / M->f[v_nb];
        double inv_s = 1.0 / fn[v_nb], diff = 0.0;
        for (int i=0;i<C->NV;i++) { double v=fn[i]*inv_s; if(fabs(v-M->f[i])>diff) diff=fabs(v-M->f[i]); M->f[i]=v; }
        for (size_t i=0;i<g_len;i++) { double v=gn[i]*inv_s; if(fabs(v-M->g[i])>diff) diff=fabs(v-M->g[i]); M->g[i]=v; }
        for (size_t i=0;i<hv_len;i++) { double v=hVn[i]*inv_s; if(fabs(v-M->hV[i])>diff) diff=fabs(v-M->hV[i]); M->hV[i]=v; }
        for (size_t i=0;i<hc_len;i++) { double v=hCn[i]*inv_s; if(fabs(v-M->hC[i])>diff) diff=fabs(v-M->hC[i]); M->hC[i]=v; }
        for (size_t i=0;i<p4_len;i++) { double v=hP4n[i]*inv_s; if(fabs(v-M->hP4[i])>diff) diff=fabs(v-M->hP4[i]); M->hP4[i]=v; }
        for (size_t i=0;i<f4_len;i++) { double v=hF4n[i]*inv_s; if(fabs(v-M->hF4[i])>diff) diff=fabs(v-M->hF4[i]); M->hF4[i]=v; }
        for (size_t i=0;i<y4_len;i++) { double v=hY4n[i]*inv_s; if(fabs(v-M->hY4[i])>diff) diff=fabs(v-M->hY4[i]); M->hY4[i]=v; }
        for (size_t i=0;i<t4_len;i++) { double v=hT4n[i]*inv_s; if(fabs(v-M->hT4[i])>diff) diff=fabs(v-M->hT4[i]); M->hT4[i]=v; }
        M->A_up*=inv_s; M->B_up*=inv_s; M->B1_lo*=inv_s; M->B2_lo*=inv_s; M->CD_lo*=inv_s;
        for (int u=0;u<S2->INNER_N;u++) for (int k=0;k<NKNOTS;k++) { M->k_pc[u][k]*=inv_s; M->k_cp[u][k]*=inv_s; M->k_cc[u][k]*=inv_s; }
        if (it % 4 == 0 || it == max_it - 1 || diff < 1e-8) {
            printf("    it=%2d lam=%.7f rho=%.8f diff=%.2e\n", it, lam, rho, diff); fflush(stdout);
        }
        if (diff < 1e-8 && it >= opt_k_iters + 3) break;
    }
    update_analytic_f_tail(M, lam);
    M->lam = lam;
    free(fn); free(gn); free(hVn); free(hCn); free(hP4n); free(hF4n); free(hY4n); free(hT4n); free(dLf_g); free(A_g);
    return rho;
}

int main(int argc, char **argv) {
    int mode  = (argc >= 2) ? atoi(argv[1]) : 1;
    double lam = (argc >= 3) ? atof(argv[2]) : (mode == 1 ? 0.56579 : 0.56565);
    int D0_g  = (argc >= 4) ? atoi(argv[3]) : 72;
    int D0_3  = (argc >= 5) ? atoi(argv[4]) : 28;
    int DHU   = (argc >= 6) ? atoi(argv[5]) : 24;
    int DCU   = (argc >= 7) ? atoi(argv[6]) : 22;
    int D0_4  = (argc >= 8) ? atoi(argv[7]) : 14;
    int DF    = 1400;

    Model M; memset(&M, 0, sizeof(M)); M.mode = mode;
    init_layout(&M.C, DF, D0_g, D0_3, DHU, DCU, D0_4);
    Layout *C = &M.C; const Level *S2 = &C->L2, *S3 = &C->L3, *S4 = &C->L4;
    size_t g_len  = (size_t)S2->NG * S2->NG;
    size_t hv_len = (size_t)S3->npars * C->NHU_ACT * S3->nchs;
    size_t hc_len = (size_t)C->NCU_ACT * S3->pair_n;
    size_t p4_len = (size_t)C->NP4_ACT * S4->npars * S4->nchs;
    size_t f4_len = (size_t)C->NF4_ACT * S4->nchs * S4->nchs;
    size_t y4_len = (size_t)C->NY4_ACT * S4->npars * S4->pair_n;
    size_t t4_len = (size_t)C->NT4_ACT * S4->tri_n;
    size_t tot = C->NV + g_len + hv_len + hc_len + p4_len + f4_len + y4_len + t4_len;

    printf("=== EXP V5 (Pure 1-Knot 3-Scale + Full 4V): mode=%s D0=(%d,%d,%d) HU=%d CU=%d ===\n",
           mode==1?"UPPER":"LOWER", D0_g, D0_3, D0_4, DHU, DCU);
    printf("Levels: L2=%d L3=%d L4=%d | Act: HU=%d CU=%d P4=%d F4=%d Y4=%d T4=%d\n",
           S2->NG, S3->NG, S4->NG, C->NHU_ACT, C->NCU_ACT, C->NP4_ACT, C->NF4_ACT, C->NY4_ACT, C->NT4_ACT);
    printf("Tables: f=%d g=%zu hV=%zu hC=%zu hP4=%zu hF4=%zu hY4=%zu hT4=%zu (tot=%.2fM, RAM~%.0fMB)\n",
           C->NV, g_len, hv_len, hc_len, p4_len, f4_len, y4_len, t4_len, tot/1e6, 2.0*8.0*tot/(1024.0*1024.0));
    fflush(stdout);

    double kn[NKNOTS] = { 0.0, 1.0/S2->D0, 1.0/48.0, 1.0/32.0, 1.0/24.0, 1.0/18.0, 1.0/14.0, 1.0/10.0, 1.0/7.0, 1.0/5.0, 1.0/4.0, 1.0/3.0, 5.0/12.0, 11.0/24.0, 0.5 };
    for (int i=0;i<NKNOTS;i++) M.KNOTS[i] = kn[i];

    M.rpar=calloc(S2->NG,sizeof(double)); M.rch=calloc(S2->NG,sizeof(double));
    M.xpar=calloc(S2->NG,sizeof(double)); M.xch=calloc(S2->NG,sizeof(double)); M.xcc=calloc(S2->NG,sizeof(double));
    M.emp_edge=calloc(g_len,sizeof(double));
    run_empirical(&M, 500000, 30);

    M.f=calloc(C->NV,sizeof(double)); M.g=calloc(g_len,sizeof(double));
    M.hV=calloc(hv_len,sizeof(double)); M.hC=calloc(hc_len,sizeof(double));
    M.hP4=calloc(p4_len+1,sizeof(double)); M.hF4=calloc(f4_len+1,sizeof(double));
    M.hY4=calloc(y4_len+1,sizeof(double)); M.hT4=calloc(t4_len+1,sizeof(double));
    M.k_pc=calloc(S2->INNER_N,sizeof(*M.k_pc)); M.k_cp=calloc(S2->INNER_N,sizeof(*M.k_cp)); M.k_cc=calloc(S2->INNER_N,sizeof(*M.k_cc));
    for (int i=0;i<C->NV;i++) M.f[i] = 1.0;

    double rho = eval_lambda(&M, lam, 15, 7);
    double lam_prev = lam, rho_prev = rho;
    lam = lam + (rho - 1.0) * 0.50;
    for (int step = 0; step < 3; step++) {
        rho = eval_lambda(&M, lam, 8, 2);
        if (fabs(rho - 1.0) < 1e-6) break;
        double slope = (rho - rho_prev) / (lam - lam_prev);
        lam_prev = lam; rho_prev = rho;
        lam = (fabs(slope) > 0.1) ? (lam - (rho - 1.0) / slope) : (lam + (rho - 1.0) * 0.50);
    }
    double lam_feas = (mode == 1) ? (lam + 3e-6) : (lam - 3e-6);
    rho = eval_lambda(&M, lam_feas, 8, 0);
    printf("  FINAL [%s] lam*=%.7f | feasible lam=%.7f rho=%.8f | CD_lo=%.5f eta_lo=%.5f\n",
           mode==1?"UPPER":"LOWER", lam, lam_feas, rho, M.CD_lo, M.eta_lo);
    verify_h1_trees_bf(&M);

    free(M.rpar); free(M.rch); free(M.xpar); free(M.xch); free(M.xcc); free(M.emp_edge);
    free(M.f); free(M.g); free(M.hV); free(M.hC); free(M.hP4); free(M.hF4); free(M.hY4); free(M.hT4);
    free(M.k_pc); free(M.k_cp); free(M.k_cc);
    free_layout(&M.C);
    return 0;
}
