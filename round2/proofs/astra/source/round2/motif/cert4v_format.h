/* Frozen Deterministic Certificate Format for 3-Scale 4-Vertex Motif Brackets (h = 1 everywhere)
 *
 * Pure 1-knot tails at all three interaction levels (L2: g, L3: hV/hC, L4: hP4/hF4/hY4/hT4).
 * Zero floating-point arithmetic in this header; all values are exact little-endian integers.
 */
#ifndef CERT4V_FORMAT_H
#define CERT4V_FORMAT_H

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define CERT4V_MAGIC      0x4654325F34563031ULL /* "FT2_4V01" */
#define CERT4V_MAX_KNOTS  16
#define CERT4V_MAX_D0     128

typedef struct {
    uint64_t magic;
    int32_t  mode;          /* 1 = UPPER, 2 = LOWER */
    int32_t  P;             /* target exponent numerator (P > 0) */
    int32_t  Q;             /* target exponent denominator (Q > 0) */
    int32_t  reserved0;     /* must be 0 */

    int64_t  S;             /* weight scale (e.g. 10^10) */
    int64_t  S_r;           /* centering & alpha scale (e.g. 10^6) */

    /* Analytic f-tail parameters for D >= DF (scaled by S):
     *   UPPER: f(q)    = (f_p0 + f_p1 * q) / S,                  with f_p2 == 0
     *   LOWER: f(q, D) = (f_p0 * q - f_p1 * q^2 - f_p2 * q / D) / S
     */
    int64_t  f_p0;
    int64_t  f_p1;
    int64_t  f_p2;
    int64_t  alpha_lo_num;  /* LOWER core-tree shift alpha = alpha_lo_num / S_r (>= 0; 0 for UPPER) */

    /* Level 1 (f) cutoff */
    int32_t  DF;

    /* Level 2 (g): pure 1-knot tail */
    int32_t  D0_2;
    int32_t  KT_2;

    /* Level 3 (hV, hC): pure 1-knot tail */
    int32_t  D0_3;
    int32_t  KT_3;
    int32_t  DHU;
    int32_t  KHU;
    int32_t  DCU;
    int32_t  KCU;

    /* Level 4 (hP4, hF4, hY4, hT4): pure 1-knot tail */
    int32_t  D0_4;
    int32_t  KT_4;
    int32_t  DP4_1;
    int32_t  DP4_2;
    int32_t  DP4_S;
    int32_t  KP4_1;
    int32_t  KP4_2;
    int32_t  DF4_1;
    int32_t  DF4_2;
    int32_t  DF4_S;
    int32_t  KF4_1;
    int32_t  KF4_2;
    int32_t  DY4;
    int32_t  KY4;
    int32_t  DT4;
    int32_t  KT4;

    /* Piecewise-linear FY knot grid xi_j = knot_num[j] / Q_K on [0, 1/2] */
    int32_t  Q_K;
    int32_t  nknots;
    int32_t  knot_num[CERT4V_MAX_KNOTS];
    int32_t  reserved1;     /* must be 0 (for 8-byte struct alignment) */
} Cert4vHeader;

typedef struct {
    int D0, KT;
    int INNER_N, INF_ID, TAIL_BASE, NG;
    int lut[CERT4V_MAX_D0][CERT4V_MAX_D0]; /* [k][L] -> inner_id */
    int *gL, *gK, *gD, *gKind;             /* 0=inner, 1=INF, 2=T1 */
    int *pars, npars, *p_to_pi;
    int *chs,  nchs,  *w_to_wi;
    size_t pair_n, tri_n;
    size_t *tri_off1, *tri_off2;
} Cert4vLevel;

typedef struct {
    Cert4vHeader H;
    int NV;
    int *f_off;
    Cert4vLevel L2, L3, L4;
    int NHU, NHU_ACT, *hu_list, *hu_to_act;
    int NCU, NCU_ACT, *cu_list, *cu_to_act;
    int NP4_ACT, *p4_map, *p4_u1, *p4_u2;
    int NF4_ACT, *f4_map, *f4_u1, *f4_u2;
    int NY4_ACT, *y4_map, *y4_u;
    int NT4_ACT, *t4_map, *t4_u;
} Cert4vLayout;

static inline int c4_lvl_inner(const Cert4vLevel *S, int L, int k) { return S->lut[k][L]; }
static inline int c4_lvl_tail(const Cert4vLevel *S, int k) {
    int kap = (k > S->KT) ? (S->KT + 1) : k;
    return S->TAIL_BASE + kap;
}
static inline size_t c4_lvl_pair(const Cert4vLevel *S, int w1, int w2) {
    int i1 = S->w_to_wi[w1], i2 = S->w_to_wi[w2];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    return (size_t)i1 * (size_t)S->nchs - (size_t)i1 * (size_t)(i1 - 1) / 2 + (size_t)(i2 - i1);
}
static inline size_t c4_lvl_tri(const Cert4vLevel *S, int w1, int w2, int w3) {
    int i1 = S->w_to_wi[w1], i2 = S->w_to_wi[w2], i3 = S->w_to_wi[w3];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    if (i2 > i3) { int t = i2; i2 = i3; i3 = t; }
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    return S->tri_off1[i1] + S->tri_off2[(size_t)i1 * (size_t)S->nchs + (size_t)i2] + (size_t)(i3 - i2);
}

static inline int c4_vidx(const Cert4vLayout *C, int L, int k) { return C->f_off[L + k + 1] + k; }
static inline int c4_hu_idx(const Cert4vLayout *C, int L, int k) {
    int D = L + k + 1;
    if (L < 0 || k < 1 || k > C->H.KHU || D >= C->H.DHU) return -1;
    return (k - 1) * C->H.DHU + L;
}
static inline int c4_cu_idx(const Cert4vLayout *C, int L, int k) {
    int D = L + k + 1;
    if (L < 0 || k < 2 || k > C->H.KCU || D >= C->H.DCU) return -1;
    return (k - 2) * C->H.DCU + L;
}
static inline int c4_p4_act(const Cert4vLayout *C, int L1, int k1, int L2, int k2) {
    int u1 = c4_hu_idx(C, L1, k1), u2 = c4_hu_idx(C, L2, k2);
    if (u1 < 0 || u2 < 0) return -1;
    return C->p4_map[(size_t)C->hu_to_act[u1] * (size_t)C->NHU_ACT + (size_t)C->hu_to_act[u2]];
}
static inline int c4_f4_act(const Cert4vLayout *C, int L1, int k1, int L2, int k2) {
    int c1 = c4_cu_idx(C, L1, k1), u2 = c4_hu_idx(C, L2, k2);
    if (c1 < 0 || u2 < 0) return -1;
    return C->f4_map[(size_t)C->cu_to_act[c1] * (size_t)C->NHU_ACT + (size_t)C->hu_to_act[u2]];
}
static inline int c4_y4_act(const Cert4vLayout *C, int L, int k) {
    int u = c4_hu_idx(C, L, k);
    return (u >= 0) ? C->y4_map[C->hu_to_act[u]] : -1;
}
static inline int c4_t4_act(const Cert4vLayout *C, int L, int k) {
    int c = c4_cu_idx(C, L, k);
    return (c >= 0) ? C->t4_map[C->cu_to_act[c]] : -1;
}

static inline size_t c4_hv_off(const Cert4vLayout *C, int p3, int hu_d, int w3) {
    return ((size_t)C->L3.p_to_pi[p3] * (size_t)C->NHU_ACT + (size_t)C->hu_to_act[hu_d]) * (size_t)C->L3.nchs + (size_t)C->L3.w_to_wi[w3];
}
static inline size_t c4_hc_off(const Cert4vLayout *C, int cu_d, int w1, int w2) {
    return (size_t)C->cu_to_act[cu_d] * C->L3.pair_n + c4_lvl_pair(&C->L3, w1, w2);
}
static inline size_t c4_p4_off(const Cert4vLayout *C, int p4a, int a4, int b4) {
    return ((size_t)p4a * (size_t)C->L4.npars + (size_t)C->L4.p_to_pi[a4]) * (size_t)C->L4.nchs + (size_t)C->L4.w_to_wi[b4];
}
static inline size_t c4_f4_off(const Cert4vLayout *C, int f4a, int a4, int b4) {
    return ((size_t)f4a * (size_t)C->L4.nchs + (size_t)C->L4.w_to_wi[a4]) * (size_t)C->L4.nchs + (size_t)C->L4.w_to_wi[b4];
}
static inline size_t c4_y4_off(const Cert4vLayout *C, int y4a, int a4, int b1, int b2) {
    return ((size_t)y4a * (size_t)C->L4.npars + (size_t)C->L4.p_to_pi[a4]) * C->L4.pair_n + c4_lvl_pair(&C->L4, b1, b2);
}
static inline size_t c4_t4_off(const Cert4vLayout *C, int t4a, int b1, int b2, int b3) {
    return (size_t)t4a * C->L4.tri_n + c4_lvl_tri(&C->L4, b1, b2, b3);
}

/* Strict header structural validator (returns 0 on OK, nonzero error code otherwise) */
static int c4_validate_header(const Cert4vHeader *H) {
    if (H->magic != CERT4V_MAGIC) return 1;
    if (H->mode != 1 && H->mode != 2) return 2;
    if (H->P <= 0 || H->Q <= 0) return 3;
    if (H->reserved0 != 0 || H->reserved1 != 0) return 4;
    if (H->S <= 0 || H->S_r <= 0) return 5;
    if (H->mode == 1 && (H->f_p0 <= 0 || H->f_p1 <= 0 || H->f_p2 != 0 || H->alpha_lo_num != 0)) return 6;
    if (H->mode == 2 && (H->f_p0 <= 0 || H->f_p1 < 0 || H->f_p2 < 0 || H->alpha_lo_num <= 0)) return 7;
    if (H->D0_4 < 3 || H->D0_4 > H->D0_3 || H->D0_3 > H->D0_2 || H->D0_2 >= H->DF) return 8;
    if (H->D0_2 >= CERT4V_MAX_D0) return 9;
    if (H->KT_4 < 0 || H->KT_4 > H->KT_3 || H->KT_3 > H->KT_2) return 10;
    if (H->KT_2 + 1 > H->D0_2 || H->KT_3 + 1 > H->D0_3 || H->KT_4 + 1 > H->D0_4) return 11;
    if (H->DHU < 2 || H->DHU > H->D0_3 || H->KHU < 1) return 12;
    if (H->DCU < 2 || H->DCU > H->DHU || H->KCU < 2 || H->KCU > H->KHU) return 13;
    if (H->DP4_1 < 2 || H->DP4_1 > H->D0_4 || H->DP4_1 > H->DHU || H->KP4_1 < 1 || H->KP4_1 > H->KHU) return 14;
    if (H->DP4_2 < 2 || H->DP4_2 > H->D0_4 || H->DP4_2 > H->DHU || H->KP4_2 < 1 || H->KP4_2 > H->KHU) return 15;
    if (H->DF4_1 < 2 || H->DF4_1 > H->D0_4 || H->DF4_1 > H->DCU || H->KF4_1 < 2 || H->KF4_1 > H->KCU) return 16;
    if (H->DF4_2 < 2 || H->DF4_2 > H->D0_4 || H->DF4_2 > H->DHU || H->KF4_2 < 1 || H->KF4_2 > H->KHU) return 17;
    if (H->DY4 < 2 || H->DY4 > H->D0_4 || H->DY4 > H->DCU || H->KY4 < 2 || H->KY4 > H->KCU) return 18;
    if (H->DT4 < 2 || H->DT4 > H->DY4 || H->KT4 < 3 || H->KT4 > H->KY4) return 19;
    if (H->Q_K <= 0 || H->nknots < 2 || H->nknots > CERT4V_MAX_KNOTS) return 20;
    if (H->knot_num[0] != 0) return 21;
    if ((int64_t)H->knot_num[1] * (int64_t)H->D0_2 != (int64_t)H->Q_K) return 22;
    if ((int64_t)2 * (int64_t)H->knot_num[H->nknots - 1] != (int64_t)H->Q_K) return 23;
    for (int i = 0; i < H->nknots - 1; i++) if (H->knot_num[i] >= H->knot_num[i + 1]) return 24;
    for (int i = H->nknots; i < CERT4V_MAX_KNOTS; i++) if (H->knot_num[i] != 0) return 25;
    return 0;
}

static void c4_init_level(Cert4vLevel *S, int D0, int KT, int build_tri) {
    memset(S, 0, sizeof(*S));
    S->D0 = D0; S->KT = KT;
    memset(S->lut, -1, sizeof(S->lut));
    int cnt = 0;
    for (int D = 2; D < D0; D++)
        for (int k = 0; k <= D - 1; k++) S->lut[k][D - 1 - k] = cnt++;
    S->INNER_N = cnt; S->INF_ID = cnt; S->TAIL_BASE = cnt + 1;
    S->NG = S->TAIL_BASE + (KT + 2);
    S->gL = (int *)calloc((size_t)S->NG, sizeof(int));
    S->gK = (int *)calloc((size_t)S->NG, sizeof(int));
    S->gD = (int *)calloc((size_t)S->NG, sizeof(int));
    S->gKind = (int *)calloc((size_t)S->NG, sizeof(int));
    for (int D = 2; D < D0; D++)
        for (int k = 0; k <= D - 1; k++) {
            int L = D - 1 - k, id = S->lut[k][L];
            S->gL[id] = L; S->gK[id] = k; S->gD[id] = D; S->gKind[id] = 0;
        }
    S->gKind[S->INF_ID] = 1; S->gD[S->INF_ID] = 0; S->gK[S->INF_ID] = -1; S->gL[S->INF_ID] = -1;
    for (int kap = 0; kap <= KT + 1; kap++) {
        int id = S->TAIL_BASE + kap;
        S->gKind[id] = 2; S->gK[id] = kap; S->gD[id] = D0; S->gL[id] = -1;
    }
    S->pars = (int *)malloc((size_t)S->NG * sizeof(int));
    S->chs  = (int *)malloc((size_t)S->NG * sizeof(int));
    S->p_to_pi = (int *)malloc((size_t)S->NG * sizeof(int));
    S->w_to_wi = (int *)malloc((size_t)S->NG * sizeof(int));
    for (int i = 0; i < S->NG; i++) { S->p_to_pi[i] = -1; S->w_to_wi[i] = -1; }
    for (int u = 0; u < S->NG; u++) {
        if (S->gKind[u] == 1 || S->gK[u] >= 1) { S->p_to_pi[u] = S->npars; S->pars[S->npars++] = u; }
        S->w_to_wi[u] = S->nchs; S->chs[S->nchs++] = u;
    }
    S->pair_n = (size_t)S->nchs * (size_t)(S->nchs + 1) / 2;
    if (build_tri) {
        int n = S->nchs;
        S->tri_off1 = (size_t *)calloc((size_t)(n + 1), sizeof(size_t));
        S->tri_off2 = (size_t *)calloc((size_t)n * (size_t)n, sizeof(size_t));
        size_t tot = 0;
        for (int i1 = 0; i1 < n; i1++) {
            S->tri_off1[i1] = tot; size_t sub = 0;
            for (int i2 = i1; i2 < n; i2++) {
                S->tri_off2[(size_t)i1 * (size_t)n + (size_t)i2] = sub;
                sub += (size_t)(n - i2);
            }
            tot += sub;
        }
        S->tri_n = tot;
    }
}

static void c4_free_level(Cert4vLevel *S) {
    free(S->gL); free(S->gK); free(S->gD); free(S->gKind);
    free(S->pars); free(S->chs); free(S->p_to_pi); free(S->w_to_wi);
    free(S->tri_off1); free(S->tri_off2);
}

static void c4_init_layout(Cert4vLayout *C, const Cert4vHeader *H) {
    memset(C, 0, sizeof(*C));
    C->H = *H;
    C->f_off = (int *)calloc((size_t)(H->DF + 2), sizeof(int));
    int nv = 0; for (int D = 2; D < H->DF; D++) { C->f_off[D] = nv; nv += D; }
    C->NV = nv;

    c4_init_level(&C->L2, H->D0_2, H->KT_2, 0);
    c4_init_level(&C->L3, H->D0_3, H->KT_3, 0);
    c4_init_level(&C->L4, H->D0_4, H->KT_4, 1);

    C->NHU = H->KHU * H->DHU;
    C->hu_to_act = (int *)malloc((size_t)C->NHU * sizeof(int));
    C->hu_list   = (int *)malloc((size_t)C->NHU * sizeof(int));
    for (int i = 0; i < C->NHU; i++) C->hu_to_act[i] = -1;
    for (int k = 1; k <= H->KHU; k++) for (int L = 0; L < H->DHU; L++) {
        int ui = c4_hu_idx(C, L, k);
        if (ui >= 0) { C->hu_to_act[ui] = C->NHU_ACT; C->hu_list[C->NHU_ACT++] = ui; }
    }

    C->NCU = (H->KCU - 1) * H->DCU;
    C->cu_to_act = (int *)malloc((size_t)C->NCU * sizeof(int));
    C->cu_list   = (int *)malloc((size_t)C->NCU * sizeof(int));
    for (int i = 0; i < C->NCU; i++) C->cu_to_act[i] = -1;
    for (int k = 2; k <= H->KCU; k++) for (int L = 0; L < H->DCU; L++) {
        int ci = c4_cu_idx(C, L, k);
        if (ci >= 0) { C->cu_to_act[ci] = C->NCU_ACT; C->cu_list[C->NCU_ACT++] = ci; }
    }

    C->p4_map = (int *)malloc((size_t)C->NHU_ACT * (size_t)C->NHU_ACT * sizeof(int));
    for (size_t i = 0; i < (size_t)C->NHU_ACT * (size_t)C->NHU_ACT; i++) C->p4_map[i] = -1;
    C->p4_u1 = (int *)malloc((size_t)C->NHU_ACT * (size_t)C->NHU_ACT * sizeof(int));
    C->p4_u2 = (int *)malloc((size_t)C->NHU_ACT * (size_t)C->NHU_ACT * sizeof(int));
    for (int a1 = 0; a1 < C->NHU_ACT; a1++) {
        int ui1 = C->hu_list[a1], k1 = (ui1 / H->DHU) + 1, L1 = ui1 % H->DHU, D1 = L1 + k1 + 1;
        if (D1 >= H->DP4_1 || k1 > H->KP4_1) continue;
        for (int a2 = 0; a2 < C->NHU_ACT; a2++) {
            int ui2 = C->hu_list[a2], k2 = (ui2 / H->DHU) + 1, L2 = ui2 % H->DHU, D2 = L2 + k2 + 1;
            if (D2 >= H->DP4_2 || k2 > H->KP4_2 || D1 + D2 >= H->DP4_S) continue;
            int id = C->NP4_ACT++;
            C->p4_map[(size_t)a1 * (size_t)C->NHU_ACT + (size_t)a2] = id;
            C->p4_u1[id] = c4_lvl_inner(&C->L4, L1, k1);
            C->p4_u2[id] = c4_lvl_inner(&C->L4, L2, k2);
        }
    }

    C->f4_map = (int *)malloc((size_t)C->NCU_ACT * (size_t)C->NHU_ACT * sizeof(int));
    for (size_t i = 0; i < (size_t)C->NCU_ACT * (size_t)C->NHU_ACT; i++) C->f4_map[i] = -1;
    C->f4_u1 = (int *)malloc((size_t)C->NCU_ACT * (size_t)C->NHU_ACT * sizeof(int));
    C->f4_u2 = (int *)malloc((size_t)C->NCU_ACT * (size_t)C->NHU_ACT * sizeof(int));
    for (int c1 = 0; c1 < C->NCU_ACT; c1++) {
        int ci1 = C->cu_list[c1], k1 = (ci1 / H->DCU) + 2, L1 = ci1 % H->DCU, D1 = L1 + k1 + 1;
        if (D1 >= H->DF4_1 || k1 > H->KF4_1) continue;
        for (int a2 = 0; a2 < C->NHU_ACT; a2++) {
            int ui2 = C->hu_list[a2], k2 = (ui2 / H->DHU) + 1, L2 = ui2 % H->DHU, D2 = L2 + k2 + 1;
            if (D2 >= H->DF4_2 || k2 > H->KF4_2 || D1 + D2 >= H->DF4_S) continue;
            int id = C->NF4_ACT++;
            C->f4_map[(size_t)c1 * (size_t)C->NHU_ACT + (size_t)a2] = id;
            C->f4_u1[id] = c4_lvl_inner(&C->L4, L1, k1);
            C->f4_u2[id] = c4_lvl_inner(&C->L4, L2, k2);
        }
    }

    C->y4_map = (int *)malloc((size_t)C->NHU_ACT * sizeof(int));
    for (int i = 0; i < C->NHU_ACT; i++) C->y4_map[i] = -1;
    C->y4_u = (int *)malloc((size_t)C->NHU_ACT * sizeof(int));
    for (int a1 = 0; a1 < C->NHU_ACT; a1++) {
        int ui1 = C->hu_list[a1], k1 = (ui1 / H->DHU) + 1, L1 = ui1 % H->DHU, D1 = L1 + k1 + 1;
        if (k1 < 2 || k1 > H->KY4 || D1 >= H->DY4 || c4_cu_idx(C, L1, k1) < 0) continue;
        int id = C->NY4_ACT++;
        C->y4_map[a1] = id;
        C->y4_u[id] = c4_lvl_inner(&C->L4, L1, k1);
    }

    C->t4_map = (int *)malloc((size_t)C->NCU_ACT * sizeof(int));
    for (int i = 0; i < C->NCU_ACT; i++) C->t4_map[i] = -1;
    C->t4_u = (int *)malloc((size_t)C->NCU_ACT * sizeof(int));
    for (int c1 = 0; c1 < C->NCU_ACT; c1++) {
        int ci1 = C->cu_list[c1], k1 = (ci1 / H->DCU) + 2, L1 = ci1 % H->DCU, D1 = L1 + k1 + 1;
        if (k1 < 3 || k1 > H->KT4 || D1 >= H->DT4) continue;
        int id = C->NT4_ACT++;
        C->t4_map[c1] = id;
        C->t4_u[id] = c4_lvl_inner(&C->L4, L1, k1);
    }
}

static void c4_free_layout(Cert4vLayout *C) {
    free(C->f_off);
    c4_free_level(&C->L2); c4_free_level(&C->L3); c4_free_level(&C->L4);
    free(C->hu_list); free(C->hu_to_act); free(C->cu_list); free(C->cu_to_act);
    free(C->p4_map); free(C->p4_u1); free(C->p4_u2);
    free(C->f4_map); free(C->f4_u1); free(C->f4_u2);
    free(C->y4_map); free(C->y4_u); free(C->t4_map); free(C->t4_u);
}

/* Exact total int64_t element counts for payload verification (zero trailing bytes) */
typedef struct {
    size_t n_cent;      /* 5 * L2.INNER_N */
    size_t n_fy;        /* 3 * L2.INNER_N * nknots */
    size_t n_f;         /* NV */
    size_t n_g;         /* L2.NG * L2.NG */
    size_t n_hV;        /* L3.npars * NHU_ACT * L3.nchs */
    size_t n_hC;        /* NCU_ACT * L3.pair_n */
    size_t n_hP4;       /* NP4_ACT * L4.npars * L4.nchs */
    size_t n_hF4;       /* NF4_ACT * L4.nchs * L4.nchs */
    size_t n_hY4;       /* NY4_ACT * L4.npars * L4.pair_n */
    size_t n_hT4;       /* NT4_ACT * L4.tri_n */
    size_t n_env_23;    /* 4 * NHU_ACT * L3.NG (dguw_m, dguw_r, dgpu_m, dgpu_r) */
    size_t n_env_p4;    /* 4 * NP4_ACT * L4.NG (m1, r1, m2, r2) */
    size_t n_env_f4;    /* 4 * NF4_ACT * L4.NG (m1, r1, m2, r2) */
    size_t n_env_y4_hv; /* 2 * NY4_ACT * L4.npars * L4.nchs (mhv, rhv) */
    size_t n_env_y4_hc; /* 2 * NY4_ACT * L4.nchs * L4.nchs (mhc, rhc; shared with T4 since T4 subset Y4) */
    size_t total_i64;
    size_t total_bytes;
} Cert4vSizes;

static void c4_compute_sizes(const Cert4vLayout *C, Cert4vSizes *Z) {
    Z->n_cent      = (size_t)5 * (size_t)C->L2.INNER_N;
    Z->n_fy        = (size_t)3 * (size_t)C->L2.INNER_N * (size_t)C->H.nknots;
    Z->n_f         = (size_t)C->NV;
    Z->n_g         = (size_t)C->L2.NG * (size_t)C->L2.NG;
    Z->n_hV        = (size_t)C->L3.npars * (size_t)C->NHU_ACT * (size_t)C->L3.nchs;
    Z->n_hC        = (size_t)C->NCU_ACT * C->L3.pair_n;
    Z->n_hP4       = (size_t)C->NP4_ACT * (size_t)C->L4.npars * (size_t)C->L4.nchs;
    Z->n_hF4       = (size_t)C->NF4_ACT * (size_t)C->L4.nchs * (size_t)C->L4.nchs;
    Z->n_hY4       = (size_t)C->NY4_ACT * (size_t)C->L4.npars * C->L4.pair_n;
    Z->n_hT4       = (size_t)C->NT4_ACT * C->L4.tri_n;
    Z->n_env_23    = (size_t)4 * (size_t)C->NHU_ACT * (size_t)C->L3.NG;
    Z->n_env_p4    = (size_t)4 * (size_t)C->NP4_ACT * (size_t)C->L4.NG;
    Z->n_env_f4    = (size_t)4 * (size_t)C->NF4_ACT * (size_t)C->L4.NG;
    Z->n_env_y4_hv = (size_t)2 * (size_t)C->NY4_ACT * (size_t)C->L4.npars * (size_t)C->L4.nchs;
    Z->n_env_y4_hc = (size_t)2 * (size_t)C->NY4_ACT * (size_t)C->L4.nchs * (size_t)C->L4.nchs;
    Z->total_i64   = Z->n_cent + Z->n_fy + Z->n_f + Z->n_g + Z->n_hV + Z->n_hC
                   + Z->n_hP4 + Z->n_hF4 + Z->n_hY4 + Z->n_hT4
                   + Z->n_env_23 + Z->n_env_p4 + Z->n_env_f4 + Z->n_env_y4_hv + Z->n_env_y4_hc;
    Z->total_bytes = sizeof(Cert4vHeader) + Z->total_i64 * sizeof(int64_t);
}

#endif /* CERT4V_FORMAT_H */
