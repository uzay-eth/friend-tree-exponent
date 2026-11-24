#ifndef CERT_B3F_FORMAT_H
#define CERT_B3F_FORMAT_H

/*
 * Binary certificate format for Unified Bernstein-3 + Full-Tail hV/hC + Level-4 (B3-FullH)
 * exact-arithmetic drift certificates for the random friend tree exponent mu.
 *
 * File layout (all little-endian fixed-width integers):
 *   1. CertB3FHeader (fixed struct)
 *   2. int64_t rpar[INNER_N]                (scaled by S_r)
 *   3. int64_t rch[INNER_N]                 (scaled by S_r)
 *   4. int64_t trpar[2*(KT+2)]              (scaled by S_r)
 *   5. int64_t trch[2*(KT+2)]               (scaled by S_r)
 *   6. int64_t c2_f[2*(KT+2)]               (scaled by S)
 *   7. int64_t c3_f[2*(KT+2)]               (scaled by S)
 *   8. int64_t Yp_ref[2*(KT+2) * nchs]      (scaled by S; indexed [tidx(kap,h)*nchs + wi])
 *   9. int64_t Yw_ref[npars * (KT+2)]       (scaled by S; indexed [pi*(KT+2) + kap], h=1)
 *  10. int64_t f[NV]                        (scaled by S)
 *  11. int64_t g[NG * NG]                   (scaled by S)
 *  12. int64_t hV[hv_len]                   (scaled by S)
 *  13. int64_t hC[hc_len]                   (scaled by S)
 *  14. int64_t qP4[p4_len]                  (scaled by S)
 *  15. int64_t qYt[yt_len]                  (scaled by S)
 *  16. int64_t qYb[yb_len]                  (scaled by S)
 *  17. int64_t qS3[s3_len]                  (scaled by S)
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define CERT_B3F_MAGIC    0x4E454F4E42334648ULL /* "NEONB3FH" */
#define CERT_B3F_MAX_D    64
/* Explicit generic-input range limits ensuring all pre-GMP native 64-bit arithmetic stays overflow-free */
#define CERT_B3F_MAX_D0   40
#define CERT_B3F_MAX_KT   64
#define CERT_B3F_MAX_LF   30000
#define CERT_B3F_MAX_KF   500
#define CERT_B3F_MAX_PQ   1000000000000LL       /* 1e12 */
#define CERT_B3F_MAX_S    100000000000000LL     /* 1e14 */
#define CERT_B3F_MAX_SR   1000000000000LL       /* 1e12 */
#define CERT_B3F_MAX_VAL  1000000000000000LL    /* 1e15 */
#define CERT_B3F_MAX_RATE 1000000000000LL       /* 1e12 */

static inline void *b3f_xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) { abort(); }
    return p;
}
static inline void *b3f_xcalloc(size_t n, size_t sz) {
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p) { abort(); }
    return p;
}
static inline void *b3f_xrealloc(void *ptr, size_t n) {
    void *p = realloc(ptr, n ? n : 1);
    if (!p) { abort(); }
    return p;
}
static inline int b3f_is_little_endian(void) {
    uint16_t e = 1;
    return *((const uint8_t *)&e) == 1;
}

typedef struct {
    uint64_t magic;
    int32_t  mode;          /* 1 = UPPER (L Phi <= (P/Q) Phi), 2 = LOWER (L Phi >= (P/Q) Phi) */
    int32_t  version;       /* = 1 */
    int64_t  P, Q;          /* exact rational exponent lambda = P / Q (P > 0, Q > 0) */
    int64_t  S;             /* exact value scale for f,g,hV,hC,qP4,qYt,qYb,qS3,c2_f,c3_f,Yp_ref,Yw_ref */
    int64_t  S_r;           /* exact rate scale for rpar,rch,trpar,trch */
    int32_t  LF, KF;        /* 1-body f truncation caps */
    int32_t  D0, KT;        /* Bernstein-3 cutoff D0 and tail core-child cap KT >= D0 - 1 */
    int32_t  DP4, KP4;      /* P4 inner cutoff */
    int32_t  DYT, KYT;      /* Ytop inner cutoff */
    int32_t  DYB, KYB;      /* Ybot inner cutoff */
    int32_t  DS3, KS3;      /* S3 inner cutoff */
} CertB3FHeader;
_Static_assert(sizeof(CertB3FHeader) == 96, "CertB3FHeader must be 96 bytes");
_Static_assert(sizeof(long) >= 8 && sizeof(size_t) >= 8, "Requires 64-bit LP64 ABI");

typedef struct {
    int LF, KF, NLF, NKF, NV, VSTAR0, VSTAR1;
    int D0, KT, NT, INNER_N, INF_ID, T1_BASE, T2_BASE, T3_BASE, NG;
    int DP4, KP4, DYT, KYT, DYB, KYB, DS3, KS3;
    int max_rank; /* inner: 2..D0-1; T3: D0; T2: D0+1; T1: D0+2; INF: D0+3 */

    int lut[2][CERT_B3F_MAX_D][CERT_B3F_MAX_D];
    int *L, *K, *H, *D, *kind; /* 0=inner, 1=INF, 2=T1, 3=T2, 4=T3 */
    int *pars, npars, *p_to_pi;
    int *chs,  nchs,  *w_to_wi;
    int *p_by_r[CERT_B3F_MAX_D + 6], p_rcnt[CERT_B3F_MAX_D + 6];
    int *c_by_r[CERT_B3F_MAX_D + 6], c_rcnt[CERT_B3F_MAX_D + 6];

    /* Full H_u (inner h=1,1<=k<=KT AND tail T1/T2/T3 h=1,1<=k<=KT)
     * and C_u (inner 2<=k<=KT AND tail T1/T2/T3 2<=k<=KT) */
    int *hu_of_g, *hu_g, NHU;
    int *cu_of_g, *cu_g, NCU;
    int *hu_by_r[CERT_B3F_MAX_D + 6], hu_rcnt[CERT_B3F_MAX_D + 6];
    int *cu_by_r[CERT_B3F_MAX_D + 6], cu_rcnt[CERT_B3F_MAX_D + 6];

    /* Level-4 inner-only domains */
    int NP4, p4_lut[CERT_B3F_MAX_D][CERT_B3F_MAX_D], *p4_L, *p4_K, *p4_hu, *p4_g;
    int NYT, yt_lut[CERT_B3F_MAX_D][CERT_B3F_MAX_D], *yt_L, *yt_K, *yt_hu, *yt_cu, *yt_g;
    int NYBU, ybu_lut[2][CERT_B3F_MAX_D][CERT_B3F_MAX_D], *ybu_L, *ybu_K, *ybu_H, *ybu_cu, *ybu_g;
    int NYBW, ybw_lut[CERT_B3F_MAX_D][CERT_B3F_MAX_D], *ybw_L, *ybw_K, *ybw_hu, *ybw_g;
    int NS3, s3_lut[2][CERT_B3F_MAX_D][CERT_B3F_MAX_D], *s3_L, *s3_K, *s3_H, *s3_cu, *s3_g;

    size_t yp_len, yw_len, g_len, hv_len, hc_len, p4_len, yt_len, yb_len, s3_len;
} B3FLayout;

static inline int b3f_vidx(const B3FLayout *C, int L, int k, int h) {
    if (k > C->KF) return h ? C->VSTAR1 : C->VSTAR0;
    if (L > C->LF) L = C->LF;
    return (h ? 0 : (C->KF + 1) * C->NLF) + k * C->NLF + L;
}
static inline int b3f_tidx(const B3FLayout *C, int kap, int h) {
    return (h ? 0 : (C->KT + 2)) + kap;
}
static inline int b3f_t1(const B3FLayout *C, int k, int h) {
    int kap = (k > C->KT) ? (C->KT + 1) : k;
    return C->T1_BASE + b3f_tidx(C, kap, h);
}
static inline int b3f_t2(const B3FLayout *C, int k, int h) {
    int kap = (k > C->KT) ? (C->KT + 1) : k;
    return C->T2_BASE + b3f_tidx(C, kap, h);
}
static inline int b3f_t3(const B3FLayout *C, int k, int h) {
    int kap = (k > C->KT) ? (C->KT + 1) : k;
    return C->T3_BASE + b3f_tidx(C, kap, h);
}
static inline int b3f_inner(const B3FLayout *C, int L, int k, int h) {
    return C->lut[h][k][L];
}
static inline int b3f_next_up(const B3FLayout *C, int u) {
    return (C->kind[u] == 0)
        ? ((C->D[u] + 1 < C->D0) ? b3f_inner(C, C->L[u] + 1, C->K[u], C->H[u])
                                 : b3f_t3(C, C->K[u], C->H[u]))
        : C->INF_ID;
}
static inline int b3f_next_bud(const B3FLayout *C, int u) {
    int uk = C->kind[u], k = C->K[u], h = C->H[u];
    if (uk == 0) return (C->L[u] >= 1) ? b3f_inner(C, C->L[u] - 1, k + 1, h) : u;
    if (uk == 2 && k <= C->KT) return b3f_t1(C, k + 1, h);
    if (uk == 3 && k <= C->KT) return b3f_t2(C, k + 1, h);
    if (uk == 4 && k <= C->KT) return b3f_t3(C, k + 1, h);
    return u;
}

static inline size_t b3f_yp_off(const B3FLayout *C, int tidx, int w) {
    return (size_t)tidx * C->nchs + C->w_to_wi[w];
}
static inline size_t b3f_yw_off(const B3FLayout *C, int p, int kap) {
    return (size_t)C->p_to_pi[p] * (C->KT + 2) + kap;
}
static inline size_t b3f_hv_off(const B3FLayout *C, int p, int ua, int w) {
    return ((size_t)C->p_to_pi[p] * C->NHU + ua) * C->nchs + C->w_to_wi[w];
}
static inline size_t b3f_hc_off(const B3FLayout *C, int ca, int w1, int w2) {
    int i1 = C->w_to_wi[w1], i2 = C->w_to_wi[w2];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    size_t pn = (size_t)C->nchs * (C->nchs + 1) / 2;
    return (size_t)ca * pn + (size_t)i1 * C->nchs - (size_t)i1 * (i1 - 1) / 2 + (size_t)(i2 - i1);
}
static inline size_t b3f_p4_off(const B3FLayout *C, int a, int b4, int c4, int d) {
    return (((size_t)C->p_to_pi[a] * C->NP4 + b4) * C->NP4 + c4) * C->nchs + C->w_to_wi[d];
}
static inline size_t b3f_yt_off(const B3FLayout *C, int p, int u4, int w1, int w2) {
    int i1 = C->w_to_wi[w1], i2 = C->w_to_wi[w2];
    if (i1 > i2) { int t = i1; i1 = i2; i2 = t; }
    size_t pn = (size_t)C->nchs * (C->nchs + 1) / 2;
    return ((size_t)C->p_to_pi[p] * C->NYT + u4) * pn + (size_t)i1 * C->nchs - (size_t)i1 * (i1 - 1) / 2 + (size_t)(i2 - i1);
}
static inline size_t b3f_yb_off(const B3FLayout *C, int u4, int w2_4, int w1, int c) {
    return (((size_t)u4 * C->NYBW + w2_4) * C->nchs + C->w_to_wi[w1]) * C->nchs + C->w_to_wi[c];
}
static inline size_t b3f_s3_off(const B3FLayout *C, int u4, int w1, int w2, int w3) {
    int a = C->w_to_wi[w1], b = C->w_to_wi[w2], c = C->w_to_wi[w3];
    if (a > b) { int t = a; a = b; b = t; }
    if (b > c) { int t = b; b = c; c = t; }
    if (a > b) { int t = a; a = b; b = t; }
    size_t n = (size_t)C->nchs, tn = n * (n + 1) * (n + 2) / 6;
    size_t na = n - (size_t)a, nb = n - (size_t)b;
    return (size_t)u4 * tn + (tn - na * (na + 1) * (na + 2) / 6)
         + (na * (na + 1) / 2 - nb * (nb + 1) / 2) + (size_t)(c - b);
}

static inline int b3f_validate_header(const CertB3FHeader *H) {
    if (!b3f_is_little_endian()) return -99;
    if (H->magic != CERT_B3F_MAGIC) return -1;
    if (H->mode != 1 && H->mode != 2) return -2;
    if (H->version != 1) return -3;
    if (H->P <= 0 || H->P > CERT_B3F_MAX_PQ || H->Q <= 0 || H->Q > CERT_B3F_MAX_PQ ||
        H->S <= 0 || H->S > CERT_B3F_MAX_S || H->S_r <= 0 || H->S_r > CERT_B3F_MAX_SR) return -4;
    if (H->D0 < 4 || H->D0 > CERT_B3F_MAX_D0) return -5;
    if (H->LF < H->D0 || H->LF > CERT_B3F_MAX_LF || H->KF < H->D0 || H->KF > CERT_B3F_MAX_KF) return -6;
    /* Require KT >= D0 - 1 so ALL inner states with k>=1,h=1 are in H_u and ALL with k>=2 are in C_u */
    if (H->KT < H->D0 - 1 || H->KT > CERT_B3F_MAX_KT || H->KT >= H->KF) return -7;
    if (H->DP4 < 3 || H->DP4 > H->D0 || H->KP4 < 1 || H->KP4 > H->KT) return -8;
    if (H->DYT < 3 || H->DYT > H->D0 || H->KYT < 2 || H->KYT > H->KT) return -9;
    if (H->DYB < 3 || H->DYB > H->D0 || H->KYB < 2 || H->KYB > H->KT) return -10;
    if (H->DS3 < 4 || H->DS3 > H->D0 || H->KS3 < 3 || H->KS3 > H->KT) return -11;
    return 0;
}

static inline void b3f_layout_init(B3FLayout *C, const CertB3FHeader *H) {
    memset(C, 0, sizeof(*C));
    C->LF = H->LF; C->KF = H->KF; C->NLF = H->LF + 1; C->NKF = H->KF + 1;
    C->VSTAR1 = 2 * C->NKF * C->NLF; C->VSTAR0 = C->VSTAR1 + 1; C->NV = C->VSTAR0 + 1;
    C->D0 = H->D0; C->KT = H->KT; C->NT = 2 * (H->KT + 2);
    C->DP4 = H->DP4; C->KP4 = H->KP4; C->DYT = H->DYT; C->KYT = H->KYT;
    C->DYB = H->DYB; C->KYB = H->KYB; C->DS3 = H->DS3; C->KS3 = H->KS3;
    C->max_rank = C->D0 + 3;

    memset(C->lut, -1, sizeof(C->lut));
    int cnt = 0;
    for (int h = 1; h >= 0; h--)
        for (int D = 2; D < C->D0; D++)
            for (int k = 0; k <= D - h; k++)
                C->lut[h][k][D - h - k] = cnt++;
    C->INNER_N = cnt;
    C->INF_ID  = cnt;
    C->T1_BASE = cnt + 1;
    C->T2_BASE = C->T1_BASE + C->NT;
    C->T3_BASE = C->T2_BASE + C->NT;
    C->NG      = C->T3_BASE + C->NT;

    C->L    = (int*)b3f_xcalloc(C->NG, sizeof(int));
    C->K    = (int*)b3f_xcalloc(C->NG, sizeof(int));
    C->H    = (int*)b3f_xcalloc(C->NG, sizeof(int));
    C->D    = (int*)b3f_xcalloc(C->NG, sizeof(int));
    C->kind = (int*)b3f_xcalloc(C->NG, sizeof(int));
    for (int h = 1; h >= 0; h--)
        for (int D = 2; D < C->D0; D++)
            for (int k = 0; k <= D - h; k++) {
                int L = D - h - k, id = C->lut[h][k][L];
                C->L[id] = L; C->K[id] = k; C->H[id] = h; C->D[id] = D; C->kind[id] = 0;
            }
    C->kind[C->INF_ID] = 1; C->D[C->INF_ID] = 0; C->K[C->INF_ID] = -1; C->H[C->INF_ID] = -1; C->L[C->INF_ID] = -1;
    for (int h = 1; h >= 0; h--)
        for (int kap = 0; kap <= C->KT + 1; kap++) {
            int i1 = b3f_t1(C, kap, h), i2 = b3f_t2(C, kap, h), i3 = b3f_t3(C, kap, h);
            C->kind[i1] = 2; C->K[i1] = kap; C->H[i1] = h; C->D[i1] = C->D0; C->L[i1] = -1;
            C->kind[i2] = 3; C->K[i2] = kap; C->H[i2] = h; C->D[i2] = C->D0; C->L[i2] = -1;
            C->kind[i3] = 4; C->K[i3] = kap; C->H[i3] = h; C->D[i3] = C->D0; C->L[i3] = -1;
        }

    C->pars    = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->chs     = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->p_to_pi = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->w_to_wi = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->hu_of_g = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->cu_of_g = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->hu_g    = (int*)b3f_xmalloc(C->NG * sizeof(int));
    C->cu_g    = (int*)b3f_xmalloc(C->NG * sizeof(int));
    for (int i = 0; i < C->NG; i++) { C->p_to_pi[i] = -1; C->w_to_wi[i] = -1; C->hu_of_g[i] = -1; C->cu_of_g[i] = -1; }
    for (int u = 0; u < C->NG; u++) {
        int uk = C->kind[u], k = C->K[u], h = C->H[u];
        if (uk == 1 || k >= 1) { C->p_to_pi[u] = C->npars; C->pars[C->npars++] = u; }
        if (uk == 1 || h == 1) { C->w_to_wi[u] = C->nchs;  C->chs[C->nchs++]   = u; }
        if ((uk == 0 || uk >= 2) && h == 1 && k >= 1 && k <= C->KT) { C->hu_of_g[u] = C->NHU; C->hu_g[C->NHU++] = u; }
        if ((uk == 0 || uk >= 2) && k >= 2 && k <= C->KT)           { C->cu_of_g[u] = C->NCU; C->cu_g[C->NCU++] = u; }
    }

    for (int r = 0; r <= C->max_rank; r++) {
        C->p_by_r[r]  = (int*)b3f_xmalloc(C->NG * sizeof(int));
        C->c_by_r[r]  = (int*)b3f_xmalloc(C->NG * sizeof(int));
        C->hu_by_r[r] = (int*)b3f_xmalloc(C->NG * sizeof(int));
        C->cu_by_r[r] = (int*)b3f_xmalloc(C->NG * sizeof(int));
    }
    for (int ip = 0; ip < C->npars; ip++) {
        int p = C->pars[ip], k = C->kind[p];
        int r = (k == 1) ? (C->D0 + 3) : ((k == 2) ? (C->D0 + 2) : ((k == 3) ? (C->D0 + 1) : ((k == 4) ? C->D0 : C->D[p])));
        C->p_by_r[r][C->p_rcnt[r]++] = p;
    }
    for (int iw = 0; iw < C->nchs; iw++) {
        int w = C->chs[iw], k = C->kind[w];
        int r = (k == 1) ? (C->D0 + 3) : ((k == 2) ? (C->D0 + 2) : ((k == 3) ? (C->D0 + 1) : ((k == 4) ? C->D0 : C->D[w])));
        C->c_by_r[r][C->c_rcnt[r]++] = w;
    }
    for (int ua = 0; ua < C->NHU; ua++) {
        int u = C->hu_g[ua], k = C->kind[u];
        int r = (k == 2) ? (C->D0 + 2) : ((k == 3) ? (C->D0 + 1) : ((k == 4) ? C->D0 : C->D[u]));
        C->hu_by_r[r][C->hu_rcnt[r]++] = ua;
    }
    for (int ca = 0; ca < C->NCU; ca++) {
        int u = C->cu_g[ca], k = C->kind[u];
        int r = (k == 2) ? (C->D0 + 2) : ((k == 3) ? (C->D0 + 1) : ((k == 4) ? C->D0 : C->D[u]));
        C->cu_by_r[r][C->cu_rcnt[r]++] = ca;
    }
    for (int r = 2; r <= C->D0 + 2; r++) {
        for (int i = 0; i < C->p_rcnt[r]; i++)
            for (int j = i + 1; j < C->p_rcnt[r]; j++)
                if (C->K[C->p_by_r[r][j]] > C->K[C->p_by_r[r][i]]) { int t = C->p_by_r[r][i]; C->p_by_r[r][i] = C->p_by_r[r][j]; C->p_by_r[r][j] = t; }
        for (int i = 0; i < C->c_rcnt[r]; i++)
            for (int j = i + 1; j < C->c_rcnt[r]; j++)
                if (C->K[C->c_by_r[r][j]] > C->K[C->c_by_r[r][i]]) { int t = C->c_by_r[r][i]; C->c_by_r[r][i] = C->c_by_r[r][j]; C->c_by_r[r][j] = t; }
        for (int i = 0; i < C->hu_rcnt[r]; i++)
            for (int j = i + 1; j < C->hu_rcnt[r]; j++)
                if (C->K[C->hu_g[C->hu_by_r[r][j]]] > C->K[C->hu_g[C->hu_by_r[r][i]]]) { int t = C->hu_by_r[r][i]; C->hu_by_r[r][i] = C->hu_by_r[r][j]; C->hu_by_r[r][j] = t; }
        for (int i = 0; i < C->cu_rcnt[r]; i++)
            for (int j = i + 1; j < C->cu_rcnt[r]; j++)
                if (C->K[C->cu_g[C->cu_by_r[r][j]]] > C->K[C->cu_g[C->cu_by_r[r][i]]]) { int t = C->cu_by_r[r][i]; C->cu_by_r[r][i] = C->cu_by_r[r][j]; C->cu_by_r[r][j] = t; }
    }

    memset(C->p4_lut, -1, sizeof(C->p4_lut));
    C->p4_L  = (int*)b3f_xmalloc((size_t)C->DP4 * C->KP4 * sizeof(int));
    C->p4_K  = (int*)b3f_xmalloc((size_t)C->DP4 * C->KP4 * sizeof(int));
    C->p4_hu = (int*)b3f_xmalloc((size_t)C->DP4 * C->KP4 * sizeof(int));
    C->p4_g  = (int*)b3f_xmalloc((size_t)C->DP4 * C->KP4 * sizeof(int));
    for (int k = 1; k <= C->KP4; k++)
        for (int L = 0; L + k + 1 < C->DP4; L++) {
            int id = C->NP4++, g = b3f_inner(C, L, k, 1);
            C->p4_lut[k][L] = id; C->p4_L[id] = L; C->p4_K[id] = k; C->p4_g[id] = g; C->p4_hu[id] = C->hu_of_g[g];
        }

    memset(C->yt_lut, -1, sizeof(C->yt_lut));
    C->yt_L  = (int*)b3f_xmalloc((size_t)C->DYT * C->KYT * sizeof(int));
    C->yt_K  = (int*)b3f_xmalloc((size_t)C->DYT * C->KYT * sizeof(int));
    C->yt_hu = (int*)b3f_xmalloc((size_t)C->DYT * C->KYT * sizeof(int));
    C->yt_cu = (int*)b3f_xmalloc((size_t)C->DYT * C->KYT * sizeof(int));
    C->yt_g  = (int*)b3f_xmalloc((size_t)C->DYT * C->KYT * sizeof(int));
    for (int k = 2; k <= C->KYT; k++)
        for (int L = 0; L + k + 1 < C->DYT; L++) {
            int id = C->NYT++, g = b3f_inner(C, L, k, 1);
            C->yt_lut[k][L] = id; C->yt_L[id] = L; C->yt_K[id] = k; C->yt_g[id] = g;
            C->yt_hu[id] = C->hu_of_g[g]; C->yt_cu[id] = C->cu_of_g[g];
        }

    memset(C->ybu_lut, -1, sizeof(C->ybu_lut));
    C->ybu_L  = (int*)b3f_xmalloc((size_t)2 * C->DYB * C->KYB * sizeof(int));
    C->ybu_K  = (int*)b3f_xmalloc((size_t)2 * C->DYB * C->KYB * sizeof(int));
    C->ybu_H  = (int*)b3f_xmalloc((size_t)2 * C->DYB * C->KYB * sizeof(int));
    C->ybu_cu = (int*)b3f_xmalloc((size_t)2 * C->DYB * C->KYB * sizeof(int));
    C->ybu_g  = (int*)b3f_xmalloc((size_t)2 * C->DYB * C->KYB * sizeof(int));
    for (int h = 1; h >= 0; h--)
        for (int k = 2; k <= C->KYB; k++)
            for (int L = 0; L + k + h < C->DYB; L++) {
                int id = C->NYBU++, g = b3f_inner(C, L, k, h);
                C->ybu_lut[h][k][L] = id; C->ybu_L[id] = L; C->ybu_K[id] = k; C->ybu_H[id] = h;
                C->ybu_g[id] = g; C->ybu_cu[id] = C->cu_of_g[g];
            }

    memset(C->ybw_lut, -1, sizeof(C->ybw_lut));
    C->ybw_L  = (int*)b3f_xmalloc((size_t)C->DYB * C->KYB * sizeof(int));
    C->ybw_K  = (int*)b3f_xmalloc((size_t)C->DYB * C->KYB * sizeof(int));
    C->ybw_hu = (int*)b3f_xmalloc((size_t)C->DYB * C->KYB * sizeof(int));
    C->ybw_g  = (int*)b3f_xmalloc((size_t)C->DYB * C->KYB * sizeof(int));
    for (int k = 1; k <= C->KYB; k++)
        for (int L = 0; L + k + 1 < C->DYB; L++) {
            int id = C->NYBW++, g = b3f_inner(C, L, k, 1);
            C->ybw_lut[k][L] = id; C->ybw_L[id] = L; C->ybw_K[id] = k; C->ybw_g[id] = g; C->ybw_hu[id] = C->hu_of_g[g];
        }

    memset(C->s3_lut, -1, sizeof(C->s3_lut));
    C->s3_L  = (int*)b3f_xmalloc((size_t)2 * C->DS3 * C->KS3 * sizeof(int));
    C->s3_K  = (int*)b3f_xmalloc((size_t)2 * C->DS3 * C->KS3 * sizeof(int));
    C->s3_H  = (int*)b3f_xmalloc((size_t)2 * C->DS3 * C->KS3 * sizeof(int));
    C->s3_cu = (int*)b3f_xmalloc((size_t)2 * C->DS3 * C->KS3 * sizeof(int));
    C->s3_g  = (int*)b3f_xmalloc((size_t)2 * C->DS3 * C->KS3 * sizeof(int));
    for (int h = 1; h >= 0; h--)
        for (int k = 3; k <= C->KS3; k++)
            for (int L = 0; L + k + h < C->DS3; L++) {
                int id = C->NS3++, g = b3f_inner(C, L, k, h);
                C->s3_lut[h][k][L] = id; C->s3_L[id] = L; C->s3_K[id] = k; C->s3_H[id] = h;
                C->s3_g[id] = g; C->s3_cu[id] = C->cu_of_g[g];
            }

    C->yp_len = (size_t)C->NT * C->nchs;
    C->yw_len = (size_t)C->npars * (C->KT + 2);
    C->g_len  = (size_t)C->NG * C->NG;
    C->hv_len = (size_t)C->npars * C->NHU * C->nchs;
    C->hc_len = (size_t)C->NCU * C->nchs * (C->nchs + 1) / 2;
    C->p4_len = (size_t)C->npars * C->NP4 * C->NP4 * C->nchs;
    C->yt_len = (size_t)C->npars * C->NYT * C->nchs * (C->nchs + 1) / 2;
    C->yb_len = (size_t)C->NYBU * C->NYBW * C->nchs * C->nchs;
    C->s3_len = (size_t)C->NS3 * C->nchs * (C->nchs + 1) * (C->nchs + 2) / 6;
}

static inline size_t b3f_expected_bytes(const B3FLayout *C) {
    size_t n_i64 = (size_t)2 * C->INNER_N + (size_t)4 * C->NT + C->yp_len + C->yw_len
                 + (size_t)C->NV + C->g_len + C->hv_len + C->hc_len
                 + C->p4_len + C->yt_len + C->yb_len + C->s3_len;
    return sizeof(CertB3FHeader) + n_i64 * sizeof(int64_t);
}

static inline void b3f_layout_free(B3FLayout *C) {
    free(C->L); free(C->K); free(C->H); free(C->D); free(C->kind);
    free(C->pars); free(C->chs); free(C->p_to_pi); free(C->w_to_wi);
    free(C->hu_of_g); free(C->cu_of_g); free(C->hu_g); free(C->cu_g);
    for (int r = 0; r <= C->max_rank; r++) {
        free(C->p_by_r[r]); free(C->c_by_r[r]);
        free(C->hu_by_r[r]); free(C->cu_by_r[r]);
    }
    free(C->p4_L); free(C->p4_K); free(C->p4_hu); free(C->p4_g);
    free(C->yt_L); free(C->yt_K); free(C->yt_hu); free(C->yt_cu); free(C->yt_g);
    free(C->ybu_L); free(C->ybu_K); free(C->ybu_H); free(C->ybu_cu); free(C->ybu_g);
    free(C->ybw_L); free(C->ybw_K); free(C->ybw_hu); free(C->ybw_g);
    free(C->s3_L); free(C->s3_K); free(C->s3_H); free(C->s3_cu); free(C->s3_g);
}

#endif
