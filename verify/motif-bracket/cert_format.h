/* Shared binary certificate format & deterministic state indexing
 * for Analytic 1/D-Tail 3-Vertex + Cherry Potential Certificates.
 */
#ifndef CERT_FORMAT_H
#define CERT_FORMAT_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CERT_MAGIC 0x494E564433563031ULL /* "INVD3V01" */
#define CERT_NKNOTS 9
#define CERT_MAX_D0 64

typedef struct {
    uint64_t magic;
    int32_t mode;       /* 1 = UPPER (LHS <= (P/Q)*Phi), 2 = LOWER (LHS >= (P/Q)*Phi) */
    int32_t P;          /* target exponent numerator */
    int32_t Q;          /* target exponent denominator */
    int32_t reserved0;
    int64_t S;          /* common scale for f, g, hV, hC, knots (e.g. 10^10) */
    int64_t S_r;        /* common scale for rpar, rch, xpar, xch, xcc (e.g. 10^6) */
    int32_t LF;
    int32_t KF;
    int32_t D0;
    int32_t KT;
    int32_t DHU;
    int32_t KHU;
    int32_t DCU;
    int32_t KCU;
    int32_t Q_K;
    int32_t nknots;
    int32_t knot_num[CERT_NKNOTS];
    int32_t reserved1;
} CertHeader;

typedef struct {
    int LF, KF, NLF, NKF, NV, VSTAR0, VSTAR1;
    int D0, KT, INNER_N, INF_ID, T1_BASE, NG;
    int DHU, KHU, NHU, NHU_ACT;
    int DCU, KCU, NCU, NCU_ACT;
    int g_lut[2][CERT_MAX_D0][CERT_MAX_D0];
    int *gL, *gK, *gH, *gD, *gKind; /* 0=inner, 1=INF, 2=T1 */
    int *pars, npars, *p_to_pi;
    int *chs,  nchs,  *w_to_wi;
    int *hu_list, *hu_to_act; /* maps dense hu_idx -> [0..NHU_ACT-1] */
    int *cu_list, *cu_to_act; /* maps dense cu_idx -> [0..NCU_ACT-1] */
} CertLayout;

static inline int cl_vidx(const CertLayout *C, int L, int k, int h) {
    if (k > C->KF) return h ? C->VSTAR1 : C->VSTAR0;
    if (L > C->LF) L = C->LF;
    return (h ? 0 : (C->KF + 1) * C->NLF) + k * C->NLF + L;
}
static inline int cl_t1_idx(const CertLayout *C, int k, int h) {
    int kap = (k > C->KT) ? (C->KT + 1) : k;
    return C->T1_BASE + (h ? 0 : (C->KT + 2)) + kap;
}
static inline int cl_inner_idx(const CertLayout *C, int L, int k, int h) {
    return C->g_lut[h][k][L];
}
static inline int cl_hu_idx(const CertLayout *C, int L, int k, int h) {
    int D = L + k + h;
    if (h != 1 || k < 1 || k > C->KHU || D >= C->DHU) return -1;
    return (k - 1) * C->DHU + L;
}
static inline int cl_cu_idx(const CertLayout *C, int L, int k, int h) {
    int D = L + k + h;
    if (k < 2 || k > C->KCU || D >= C->DCU) return -1;
    return (h ? 0 : (C->KCU - 1) * C->DCU) + (k - 2) * C->DCU + L;
}

static inline size_t cl_hv_off(const CertLayout *C, int p, int hu_dense, int w) {
    int pi = C->p_to_pi[p], ua = C->hu_to_act[hu_dense], wi = C->w_to_wi[w];
    return ((size_t)pi * C->NHU_ACT + ua) * C->nchs + wi;
}
static inline size_t cl_hc_off(const CertLayout *C, int cu_dense, int w1, int w2) {
    int ca = C->cu_to_act[cu_dense], i1 = C->w_to_wi[w1], i2 = C->w_to_wi[w2];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    size_t pair_n = (size_t)C->nchs * (C->nchs + 1) / 2;
    size_t p_idx = (size_t)i1 * C->nchs - (size_t)i1 * (i1 - 1) / 2 + (i2 - i1);
    return (size_t)ca * pair_n + p_idx;
}

static void cl_init(CertLayout *C, const CertHeader *H) {
    memset(C, 0, sizeof(*C));
    C->LF = H->LF; C->KF = H->KF; C->NLF = H->LF + 1; C->NKF = H->KF + 1;
    C->VSTAR1 = 2 * C->NKF * C->NLF; C->VSTAR0 = C->VSTAR1 + 1; C->NV = C->VSTAR0 + 1;
    C->D0 = H->D0; C->KT = H->KT;
    memset(C->g_lut, -1, sizeof(C->g_lut));
    int cnt = 0;
    for (int h = 1; h >= 0; h--)
        for (int D = 2; D < C->D0; D++)
            for (int k = 0; k <= D - h; k++)
                C->g_lut[h][k][D - h - k] = cnt++;
    C->INNER_N = cnt;
    C->INF_ID = cnt;
    C->T1_BASE = cnt + 1;
    C->NG = C->T1_BASE + 2 * (C->KT + 2);

    C->gL = calloc(C->NG, sizeof(int));
    C->gK = calloc(C->NG, sizeof(int));
    C->gH = calloc(C->NG, sizeof(int));
    C->gD = calloc(C->NG, sizeof(int));
    C->gKind = calloc(C->NG, sizeof(int));
    for (int h = 1; h >= 0; h--)
        for (int D = 2; D < C->D0; D++)
            for (int k = 0; k <= D - h; k++) {
                int L = D - h - k, id = C->g_lut[h][k][L];
                C->gL[id]=L; C->gK[id]=k; C->gH[id]=h; C->gD[id]=D; C->gKind[id]=0;
            }
    C->gKind[C->INF_ID] = 1; C->gD[C->INF_ID] = 0; C->gK[C->INF_ID] = -1; C->gH[C->INF_ID] = -1; C->gL[C->INF_ID] = -1;
    for (int h = 1; h >= 0; h--)
        for (int kap = 0; kap <= C->KT + 1; kap++) {
            int id = cl_t1_idx(C, kap, h);
            C->gKind[id] = 2; C->gK[id] = kap; C->gH[id] = h; C->gD[id] = C->D0; C->gL[id] = -1;
        }

    C->pars = malloc(C->NG * sizeof(int));
    C->chs  = malloc(C->NG * sizeof(int));
    C->p_to_pi = malloc(C->NG * sizeof(int));
    C->w_to_wi = malloc(C->NG * sizeof(int));
    for (int i = 0; i < C->NG; i++) { C->p_to_pi[i] = -1; C->w_to_wi[i] = -1; }
    for (int u = 0; u < C->NG; u++) {
        int is_p = (C->gKind[u] == 1) || (C->gK[u] >= 1);
        int is_w = (C->gKind[u] == 1) || (C->gH[u] == 1);
        if (is_p) { C->p_to_pi[u] = C->npars; C->pars[C->npars++] = u; }
        if (is_w) { C->w_to_wi[u] = C->nchs;  C->chs[C->nchs++]   = u; }
    }

    C->DHU = H->DHU; C->KHU = H->KHU; C->NHU = C->KHU * C->DHU;
    C->hu_to_act = malloc(C->NHU * sizeof(int));
    C->hu_list   = malloc(C->NHU * sizeof(int));
    for (int i = 0; i < C->NHU; i++) C->hu_to_act[i] = -1;
    for (int k = 1; k <= C->KHU; k++)
        for (int L = 0; L < C->DHU; L++) {
            int ui = cl_hu_idx(C, L, k, 1);
            if (ui >= 0) { C->hu_to_act[ui] = C->NHU_ACT; C->hu_list[C->NHU_ACT++] = ui; }
        }

    C->DCU = H->DCU; C->KCU = H->KCU;
    C->NCU = 2 * (C->KCU - 1) * C->DCU;
    C->cu_to_act = malloc(C->NCU * sizeof(int));
    C->cu_list   = malloc(C->NCU * sizeof(int));
    for (int i = 0; i < C->NCU; i++) C->cu_to_act[i] = -1;
    for (int h = 1; h >= 0; h--)
        for (int k = 2; k <= C->KCU; k++)
            for (int L = 0; L < C->DCU; L++) {
                int ci = cl_cu_idx(C, L, k, h);
                if (ci >= 0) { C->cu_to_act[ci] = C->NCU_ACT; C->cu_list[C->NCU_ACT++] = ci; }
            }
}

static void cl_free(CertLayout *C) {
    free(C->gL); free(C->gK); free(C->gH); free(C->gD); free(C->gKind);
    free(C->pars); free(C->chs); free(C->p_to_pi); free(C->w_to_wi);
    free(C->hu_list); free(C->hu_to_act); free(C->cu_list); free(C->cu_to_act);
}

#endif
