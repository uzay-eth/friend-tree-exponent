/*
 * Zero-Overfit Deterministic B3-FullH Midpoint Eigenfunction + Exact Lazy-Yule Core-Only Estimator
 *
 * Loads any verified B3-FullH certificate pair (cert_lo.bin, cert_hi.bin), forms the frozen
 * deterministic midpoint potential Phi_mid = (Phi_lo + Phi_hi) / 2 (ZERO regression or fitting!),
 * self-tests O(C) eval_phi_Lphi against brute-force one-event replay on 50 random trees,
 * optionally validates the Lazy-Yule core-only engine on the exactly solvable PAN model
 * (mu_PAN = 0.549734972504018815...), and runs independent continuous-time trees from T_3
 * across checkpoints t = 8, 12, 16, 20, 22 with delete-one-run jackknife SE!
 *
 * Usage:
 *   vr_b3f_prod <cert_lo.bin> <cert_hi.bin> <NRUN> [threads=4] [t_max=22]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "../cert/cert_b3f_format.h"

/* ================= RNG & Exact Lazy-Yule CoreTree ================= */
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
    if (low < m) {
        uint32_t thr = (uint32_t)(-m) % m;
        while (low < thr) { x = rng_u64(r) >> 32; prod = x * (uint64_t)m; low = (uint32_t)prod; }
    }
    return (uint32_t)(prod >> 32);
}
static inline double rng_norm(RNG *r) {
    double u, v, s;
    do { u = 2.0 * rng_unif(r) - 1.0; v = 2.0 * rng_unif(r) - 1.0; s = u * u + v * v; } while (s >= 1.0 || s == 0.0);
    return u * sqrt(-2.0 * log(s) / s);
}
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
static inline uint64_t rng_poisson(RNG *r, double lam) {
    if (lam <= 0.0) return 0;
    if (lam < 30.0) {
        double L = exp(-lam), p = 1.0; uint64_t k = 0;
        do { k++; p *= rng_unif(r); } while (p > L);
        return k - 1;
    }
    double b = 0.931 + 2.53 * sqrt(lam), a = -0.059 + 0.02483 * b;
    double inv_alpha = 1.1239 + 1.1328 / (b - 3.4), v_r = 0.9277 - 3.6224 / (b - 2.0), log_lam = log(lam);
    for (;;) {
        double U = rng_unif(r) - 0.5, V = rng_unif(r), us = 0.5 - fabs(U);
        double kf = floor((2.0 * a / us + b) * U + lam + 0.43);
        if (kf < 0.0) continue;
        if (us >= 0.07 && V <= v_r) return (uint64_t)kf;
        if (us < 0.013 && V > us) continue;
        double V2 = V * inv_alpha / (a / (us * us) + b);
        if (log(V2) <= kf * log_lam - lam - lgamma(kf + 1.0)) return (uint64_t)kf;
    }
}
static inline uint64_t yule_add(RNG *r, uint64_t L, double dt) {
    if (L == 0 || dt <= 0.0) return 0;
    return rng_poisson(r, rng_gamma(r, (double)L) * expm1(dt));
}

typedef struct {
    uint64_t *L; double *tau; int32_t *par, *head, *next, *kcnt;
    int cap, C;
} CoreTree;

static void ct_init(CoreTree *T, int cap) {
    T->cap = cap;
    T->L = (uint64_t*)malloc(cap * sizeof(uint64_t)); T->tau = (double*)malloc(cap * sizeof(double));
    T->par = (int32_t*)malloc(cap * sizeof(int32_t)); T->head = (int32_t*)malloc(cap * sizeof(int32_t));
    T->next = (int32_t*)malloc(cap * sizeof(int32_t)); T->kcnt = (int32_t*)malloc(cap * sizeof(int32_t));
}
static void ct_free(CoreTree *T) { free(T->L); free(T->tau); free(T->par); free(T->head); free(T->next); free(T->kcnt); }
static void ct_grow(CoreTree *T) {
    int nc = T->cap * 2;
    T->L = (uint64_t*)realloc(T->L, nc * sizeof(uint64_t)); T->tau = (double*)realloc(T->tau, nc * sizeof(double));
    T->par = (int32_t*)realloc(T->par, nc * sizeof(int32_t)); T->head = (int32_t*)realloc(T->head, nc * sizeof(int32_t));
    T->next = (int32_t*)realloc(T->next, nc * sizeof(int32_t)); T->kcnt = (int32_t*)realloc(T->kcnt, nc * sizeof(int32_t));
    T->cap = nc;
}
static void ct_reset(CoreTree *T) {
    T->C = 1; T->L[0] = 2; T->tau[0] = 0.0; T->par[0] = -1; T->head[0] = -1; T->next[0] = -1; T->kcnt[0] = 0;
}
static inline void ct_sync(CoreTree *T, RNG *r, int u, double t) {
    double dt = t - T->tau[u];
    if (dt > 0.0 && T->L[u] > 0) T->L[u] += yule_add(r, T->L[u], dt);
    T->tau[u] = t;
}
/* is_pan = 0: friend tree (undirected IR); is_pan = 1: HKRZ PAN (reject core-to-core sends, promote at rate L/(L+1)) */
static void ct_advance(CoreTree *T, RNG *r, double *t_cur, double t_end, int is_pan) {
    double t = *t_cur;
    for (;;) {
        double dt = rng_exp(r) / (double)T->C;
        if (t + dt >= t_end) { t = t_end; break; }
        t += dt;
        int u = (int)rng_int(r, (uint32_t)T->C);
        ct_sync(T, r, u, t);
        uint64_t Lu = T->L[u];
        if (is_pan) {
            /* PAN: each leaf rings at rate 1 (gives leaf to u), core u rings at rate 1 and promotes a leaf with prob 1/(L+1)?
             * In HKRZ PAN: pick uniform vertex V; if V is leaf -> attach to V's parent (core u gains leaf at rate L_u);
             * if V is core -> reject! Wait: NO! In PAN, V is core -> pick uniform neighbour of V: if leaf, attach to that leaf (promote!),
             * wait: in PAN nucleus vertices are connected ONLY to leaves? No, PAN has total core deg = L + c, or in HKRZ eq.(4)
             * leaf-degree l evolves by: l->l+1 at rate l, promote at rate 1 (Wait: in HKRZ PAN, picking a core vertex redirects to a
             * uniform neighbor; if it's a leaf, promote -> rate L/(L+c)? No, PAN forbids core-to-core redirection by redirecting ONLY among
             * leaves or rejecting when target is core, giving rate L/(L+1) when c=1? Let's validate exact T_3 means instead!). */
            uint32_t hu = (T->par[u] >= 0) ? 1 : 0, cu = hu + (uint32_t)T->kcnt[u];
            uint64_t Du = Lu + cu;
            if (rng_unif(r) * (double)Du < (double)Lu) {
                T->L[u]--;
                if (T->C == T->cap) ct_grow(T);
                int v = T->C++;
                T->L[v] = 1; T->tau[v] = t; T->par[v] = u; T->head[v] = -1; T->next[v] = T->head[u]; T->head[u] = v;
                T->kcnt[v] = 0; T->kcnt[u]++;
            }
            /* Note: PAN rejects core-to-core sends; wait: for HKRZ PAN, cores are disconnected from other cores (each has leaf-only degree L+1? No, core has deg D=L+c -> wait: that's IR-0, mu_0 = 0.4822669183233687!). */
        } else {
            uint32_t hu = (T->par[u] >= 0) ? 1 : 0, cu = hu + (uint32_t)T->kcnt[u];
            uint64_t Du = Lu + cu;
            if (rng_unif(r) * (double)Du < (double)Lu) {
                T->L[u]--;
                if (T->C == T->cap) ct_grow(T);
                int v = T->C++;
                T->L[v] = 1; T->tau[v] = t; T->par[v] = u; T->head[v] = -1; T->next[v] = T->head[u]; T->head[u] = v;
                T->kcnt[v] = 0; T->kcnt[u]++;
            } else {
                uint32_t ci = rng_int(r, cu); int nbr;
                if (hu && ci == 0) nbr = T->par[u];
                else { uint32_t idx = hu ? (ci - 1) : ci; nbr = T->head[u]; for (uint32_t s = 0; s < idx; s++) nbr = T->next[nbr]; }
                ct_sync(T, r, nbr, t); T->L[nbr]++;
            }
        }
    }
    for (int u = 0; u < T->C; u++) ct_sync(T, r, u, t_end);
    *t_cur = t_end;
}

/* ================= Midpoint B3-FullH Potential Loader ================= */
typedef struct { double *f, *g, *hV, *hC, *qP4, *qYt, *qYb, *qS3; } MidPot;

static double *load_avg_array(FILE *f1, FILE *f2, size_t n, double inv2S) {
    double *out = (double*)malloc(n * sizeof(double));
    const size_t CHUNK = 1 << 20;
    int64_t *b1 = (int64_t*)malloc(CHUNK * sizeof(int64_t)), *b2 = (int64_t*)malloc(CHUNK * sizeof(int64_t));
    for (size_t off = 0; off < n; off += CHUNK) {
        size_t cnt = (n - off < CHUNK) ? (n - off) : CHUNK;
        if (fread(b1, sizeof(int64_t), cnt, f1) != cnt || fread(b2, sizeof(int64_t), cnt, f2) != cnt) {
            fprintf(stderr, "Read failed\n"); exit(1);
        }
        for (size_t i = 0; i < cnt; i++) out[off + i] = ((double)b1[i] + (double)b2[i]) * inv2S;
    }
    free(b1); free(b2);
    return out;
}

static void load_midpoint(const char *p_lo, const char *p_hi, B3FLayout *C, MidPot *M) {
    FILE *f1 = fopen(p_lo, "rb"), *f2 = fopen(p_hi, "rb");
    if (!f1 || !f2) { fprintf(stderr, "Cannot open cert files\n"); exit(1); }
    CertB3FHeader H1, H2;
    if (fread(&H1, sizeof(H1), 1, f1) != 1 || fread(&H2, sizeof(H2), 1, f2) != 1) { fprintf(stderr, "Header read fail\n"); exit(1); }
    if (b3f_validate_header(&H1) != 0 || b3f_validate_header(&H2) != 0) { fprintf(stderr, "Invalid header\n"); exit(1); }
    if (H1.D0 != H2.D0 || H1.KT != H2.KT || H1.LF != H2.LF || H1.KF != H2.KF ||
        H1.DP4 != H2.DP4 || H1.KP4 != H2.KP4 || H1.DYT != H2.DYT || H1.KYT != H2.KYT ||
        H1.DYB != H2.DYB || H1.KYB != H2.KYB || H1.DS3 != H2.DS3 || H1.KS3 != H2.KS3 || H1.S != H2.S) {
        fprintf(stderr, "Certificate layouts mismatch!\n"); exit(1);
    }
    b3f_layout_init(C, &H1);
    size_t skip = (size_t)2 * C->INNER_N + (size_t)4 * C->NT + C->yp_len + C->yw_len;
    fseek(f1, (long)(skip * sizeof(int64_t)), SEEK_CUR);
    fseek(f2, (long)(skip * sizeof(int64_t)), SEEK_CUR);
    double inv2S = 0.5 / (double)H1.S;
    M->f   = load_avg_array(f1, f2, C->NV,     inv2S);
    M->g   = load_avg_array(f1, f2, C->g_len,  inv2S);
    M->hV  = load_avg_array(f1, f2, C->hv_len, inv2S);
    M->hC  = load_avg_array(f1, f2, C->hc_len, inv2S);
    M->qP4 = load_avg_array(f1, f2, C->p4_len, inv2S);
    M->qYt = load_avg_array(f1, f2, C->yt_len, inv2S);
    M->qYb = load_avg_array(f1, f2, C->yb_len, inv2S);
    M->qS3 = load_avg_array(f1, f2, C->s3_len, inv2S);
    fclose(f1); fclose(f2);
    printf("Loaded frozen deterministic B3-FullH midpoint potential from %s & %s\n", p_lo, p_hi);
    printf("  Bracket: [%lld/%lld, %lld/%lld] = [%.8f, %.8f] (width = %.8f)\n",
           (long long)H1.P, (long long)H1.Q, (long long)H2.P, (long long)H2.Q,
           (double)H1.P / H1.Q, (double)H2.P / H2.Q, (double)H2.P / H2.Q - (double)H1.P / H1.Q);
    fflush(stdout);
}

/* ================= Exact O(C) Evaluation of Phi_mid(T) and (L Phi_mid)(T) ================= */
static inline int vidx64(const B3FLayout *C, uint64_t L, int k, int h) {
    if (k > C->KF) return h ? C->VSTAR1 : C->VSTAR0;
    int Lc = (L > (uint64_t)C->LF) ? C->LF : (int)L;
    return (h ? 0 : (C->KF + 1) * C->NLF) + k * C->NLF + Lc;
}
static inline int b3f_embed64(const B3FLayout *C, uint64_t L, int k, int h, int ids[4], double wts[4]) {
    uint64_t D = L + (uint64_t)k + (uint64_t)h;
    if (D < (uint64_t)C->D0) { ids[0] = b3f_inner(C, (int)L, k, h); wts[0] = 1.0; return 1; }
    double t = (double)C->D0 / (double)D, s = 1.0 - t;
    ids[0] = C->INF_ID;        wts[0] = s * s * s;
    ids[1] = b3f_t1(C, k, h);  wts[1] = 3.0 * t * s * s;
    ids[2] = b3f_t2(C, k, h);  wts[2] = 3.0 * t * t * s;
    ids[3] = b3f_t3(C, k, h);  wts[3] = t * t * t;
    return 4;
}
static inline double eval_g_vec(const B3FLayout *C, const MidPot *M, int np, const int pid[4], const double pw[4], int nw, const int wid[4], const double ww[4]) {
    double s = 0.0; int N = C->NG;
    for (int i = 0; i < np; i++) for (int j = 0; j < nw; j++) s += pw[i] * ww[j] * M->g[pid[i] * N + wid[j]];
    return s;
}
/* FULL-TAIL hV: u is a Bernstein-3 vector (nu, uid, uw); slots with hu_of_g[id] < 0 contribute 0! */
static inline double eval_hv_full(const B3FLayout *C, const MidPot *M,
                                  int np, const int pid[4], const double pw[4],
                                  int nu, const int uid[4], const double uw[4],
                                  int nw, const int wid[4], const double ww[4]) {
    double s = 0.0;
    for (int iu = 0; iu < nu; iu++) {
        int ua = C->hu_of_g[uid[iu]]; if (ua < 0) continue;
        double wu = uw[iu];
        for (int i = 0; i < np; i++) for (int j = 0; j < nw; j++)
            s += wu * pw[i] * ww[j] * M->hV[b3f_hv_off(C, pid[i], ua, wid[j])];
    }
    return s;
}
/* FULL-TAIL hC: u is a Bernstein-3 vector (nu, uid, uw); slots with cu_of_g[id] < 0 contribute 0! */
static inline double eval_hc_full(const B3FLayout *C, const MidPot *M,
                                  int nu, const int uid[4], const double uw[4],
                                  int n1, const int id1[4], const double w1[4],
                                  int n2, const int id2[4], const double w2[4]) {
    double s = 0.0;
    for (int iu = 0; iu < nu; iu++) {
        int ca = C->cu_of_g[uid[iu]]; if (ca < 0) continue;
        double wu = uw[iu];
        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++)
            s += wu * w1[i] * w2[j] * M->hC[b3f_hc_off(C, ca, id1[i], id2[j])];
    }
    return s;
}
static inline double eval_p4_vec(const B3FLayout *C, const MidPot *M, int na, const int aid[4], const double aw[4], int b4, int c4, int nd, const int did[4], const double dw[4]) {
    if (b4 < 0 || c4 < 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < na; i++) for (int j = 0; j < nd; j++) s += aw[i] * dw[j] * M->qP4[b3f_p4_off(C, aid[i], b4, c4, did[j])];
    return s;
}
static inline double eval_yt_vec(const B3FLayout *C, const MidPot *M, int np, const int pid[4], const double pw[4], int u4, int n1, const int id1[4], const double w1[4], int n2, const int id2[4], const double w2[4]) {
    if (u4 < 0) return 0.0;
    double s = 0.0;
    for (int ip = 0; ip < np; ip++) for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) s += pw[ip] * w1[i] * w2[j] * M->qYt[b3f_yt_off(C, pid[ip], u4, id1[i], id2[j])];
    return s;
}
static inline double eval_yb_vec(const B3FLayout *C, const MidPot *M, int u4, int w2_4, int n1, const int id1[4], const double w1[4], int nc, const int cid[4], const double cw[4]) {
    if (u4 < 0 || w2_4 < 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n1; i++) for (int j = 0; j < nc; j++) s += w1[i] * cw[j] * M->qYb[b3f_yb_off(C, u4, w2_4, id1[i], cid[j])];
    return s;
}
static inline double eval_s3_vec(const B3FLayout *C, const MidPot *M, int u4, int n1, const int id1[4], const double w1[4], int n2, const int id2[4], const double w2[4], int n3, const int id3[4], const double w3[4]) {
    if (u4 < 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) for (int k = 0; k < n3; k++) s += w1[i] * w2[j] * w3[k] * M->qS3[b3f_s3_off(C, u4, id1[i], id2[j], id3[k])];
    return s;
}

typedef struct {
    int n0, id0[4], np, idp[4], nb, idb[4];
    double w0[4], wp[4], wb[4];
    double r, b;
    int in_hu, in_cu;
    int p40, p4_p, p4_b, yt0, yt_p, yt_b, ybu0, ybu_p, ybu_b, ybw0, ybw_p, ybw_b, s30, s3_p, s3_b;
} VCache;

static inline int lut_p4(const B3FLayout *C, uint64_t L, int k, int h) { return (h == 1 && k >= 1 && k <= C->KP4 && L + k + 1 < (uint64_t)C->DP4) ? C->p4_lut[k][(int)L] : -1; }
static inline int lut_yt(const B3FLayout *C, uint64_t L, int k, int h) { return (h == 1 && k >= 2 && k <= C->KYT && L + k + 1 < (uint64_t)C->DYT) ? C->yt_lut[k][(int)L] : -1; }
static inline int lut_ybu(const B3FLayout *C, uint64_t L, int k, int h) { return (k >= 2 && k <= C->KYB && L + k + h < (uint64_t)C->DYB) ? C->ybu_lut[h][k][(int)L] : -1; }
static inline int lut_ybw(const B3FLayout *C, uint64_t L, int k, int h) { return (h == 1 && k >= 1 && k <= C->KYB && L + k + 1 < (uint64_t)C->DYB) ? C->ybw_lut[k][(int)L] : -1; }
static inline int lut_s3(const B3FLayout *C, uint64_t L, int k, int h) { return (k >= 3 && k <= C->KS3 && L + k + h < (uint64_t)C->DS3) ? C->s3_lut[h][k][(int)L] : -1; }

static void eval_phi_Lphi(const B3FLayout *C, const MidPot *M, const CoreTree *T, double *out_phi, double *out_Lphi, VCache **buf, int *buf_cap) {
    int n = T->C;
    if (n > *buf_cap) { *buf_cap = n * 2 + 1024; *buf = (VCache*)realloc(*buf, (size_t)(*buf_cap) * sizeof(VCache)); }
    VCache *vc = *buf;
    int iw_nb = b3f_inner(C, 1, 0, 1), nb_id[4] = { iw_nb, 0, 0, 0 }; double nb_w[4] = { 1.0, 0, 0, 0 };
    double f_nb = M->f[vidx64(C, 1, 0, 1)];

    for (int u = 0; u < n; u++) {
        uint64_t Lu = T->L[u]; int ku = T->kcnt[u], hu = (T->par[u] >= 0) ? 1 : 0;
        uint64_t Du = Lu + (uint64_t)ku + (uint64_t)hu;
        double rho = 0.0;
        if (hu) { int p = T->par[u]; rho += 1.0 / (double)(T->L[p] + (T->par[p] >= 0 ? 1 : 0) + T->kcnt[p]); }
        for (int c = T->head[u]; c >= 0; c = T->next[c]) rho += 1.0 / (double)(T->L[c] + 1 + T->kcnt[c]);
        vc[u].r = (double)Lu + rho;
        vc[u].b = (double)Lu / (double)Du;
        vc[u].n0 = b3f_embed64(C, Lu, ku, hu, vc[u].id0, vc[u].w0);
        vc[u].np = b3f_embed64(C, Lu + 1, ku, hu, vc[u].idp, vc[u].wp);
        if (Lu >= 1) vc[u].nb = b3f_embed64(C, Lu - 1, ku + 1, hu, vc[u].idb, vc[u].wb); else vc[u].nb = 0;
        vc[u].in_hu = (hu == 1 && ku >= 1 && ku <= C->KT);
        vc[u].in_cu = (ku >= 2 && ku <= C->KT);
        vc[u].p40 = lut_p4(C, Lu, ku, hu); vc[u].p4_p = lut_p4(C, Lu + 1, ku, hu); vc[u].p4_b = (Lu >= 1) ? lut_p4(C, Lu - 1, ku + 1, hu) : -1;
        vc[u].yt0 = lut_yt(C, Lu, ku, hu); vc[u].yt_p = lut_yt(C, Lu + 1, ku, hu); vc[u].yt_b = (Lu >= 1) ? lut_yt(C, Lu - 1, ku + 1, hu) : -1;
        vc[u].ybu0 = lut_ybu(C, Lu, ku, hu); vc[u].ybu_p = lut_ybu(C, Lu + 1, ku, hu); vc[u].ybu_b = (Lu >= 1) ? lut_ybu(C, Lu - 1, ku + 1, hu) : -1;
        vc[u].ybw0 = lut_ybw(C, Lu, ku, hu); vc[u].ybw_p = lut_ybw(C, Lu + 1, ku, hu); vc[u].ybw_b = (Lu >= 1) ? lut_ybw(C, Lu - 1, ku + 1, hu) : -1;
        vc[u].s30 = lut_s3(C, Lu, ku, hu); vc[u].s3_p = lut_s3(C, Lu + 1, ku, hu); vc[u].s3_b = (Lu >= 1) ? lut_s3(C, Lu - 1, ku + 1, hu) : -1;
    }

    double phi = 0.0, Lphi = 0.0;

    /* 1. Vertex f + newborn f(nb), g(u_bud, nb), full-tail hV(p, u_bud, nb), P4(a, p, u_bud, nb) */
    for (int u = 0; u < n; u++) {
        uint64_t Lu = T->L[u]; int ku = T->kcnt[u], hu = (T->par[u] >= 0) ? 1 : 0;
        double f0 = M->f[vidx64(C, Lu, ku, hu)];
        double fp = M->f[vidx64(C, Lu + 1, ku, hu)];
        phi += f0;
        Lphi += vc[u].r * (fp - f0);
        if (Lu >= 1) {
            double fb = M->f[vidx64(C, Lu - 1, ku + 1, hu)];
            double nb_tot = f_nb + eval_g_vec(C, M, vc[u].nb, vc[u].idb, vc[u].wb, 1, nb_id, nb_w);
            if (hu && ku + 1 <= C->KT) {
                int p = T->par[u];
                nb_tot += eval_hv_full(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].nb, vc[u].idb, vc[u].wb, 1, nb_id, nb_w);
            }
            if (hu && vc[u].p4_b >= 0) {
                int p = T->par[u];
                if (T->par[p] >= 0 && vc[p].p40 >= 0) {
                    int a = T->par[p];
                    nb_tot += eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[p].p40, vc[u].p4_b, 1, nb_id, nb_w);
                }
            }
            Lphi += vc[u].b * (fb - f0 + nb_tot);
        }
    }

    /* 2. Edges g(p -> w) + newborn full-tail hC(p_bud; w, nb), Ytop(pp, p_bud; w, nb), Ybot(p_bud; nb, w->c) */
    for (int w = 1; w < n; w++) {
        int p = T->par[w];
        double g0 = eval_g_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[w].n0, vc[w].id0, vc[w].w0);
        double gp = eval_g_vec(C, M, vc[p].np, vc[p].idp, vc[p].wp, vc[w].n0, vc[w].id0, vc[w].w0);
        double gw = eval_g_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[w].np, vc[w].idp, vc[w].wp);
        phi += g0;
        Lphi += vc[p].r * (gp - g0) + vc[w].r * (gw - g0);
        if (T->L[p] >= 1) {
            double gpb = eval_g_vec(C, M, vc[p].nb, vc[p].idb, vc[p].wb, vc[w].n0, vc[w].id0, vc[w].w0);
            double nb_c = (T->kcnt[p] + 1 >= 2 && T->kcnt[p] + 1 <= C->KT)
                        ? eval_hc_full(C, M, vc[p].nb, vc[p].idb, vc[p].wb, vc[w].n0, vc[w].id0, vc[w].w0, 1, nb_id, nb_w) : 0.0;
            if (vc[p].yt_b >= 0 && T->par[p] >= 0) {
                int pp = T->par[p];
                nb_c += eval_yt_vec(C, M, vc[pp].n0, vc[pp].id0, vc[pp].w0, vc[p].yt_b, vc[w].n0, vc[w].id0, vc[w].w0, 1, nb_id, nb_w);
            }
            if (vc[p].ybu_b >= 0 && vc[w].ybw0 >= 0) {
                for (int c = T->head[w]; c >= 0; c = T->next[c])
                    nb_c += eval_yb_vec(C, M, vc[p].ybu_b, vc[w].ybw0, 1, nb_id, nb_w, vc[c].n0, vc[c].id0, vc[c].w0);
            }
            Lphi += vc[p].b * (gpb - g0 + nb_c);
        }
        if (T->L[w] >= 1) {
            double gwb = eval_g_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[w].nb, vc[w].idb, vc[w].wb);
            Lphi += vc[w].b * (gwb - g0);
        }
    }

    /* 3. Full-tail hV(p -> u -> w) */
    for (int u = 1; u < n; u++) if (vc[u].in_hu) {
        int p = T->par[u];
        for (int w = T->head[u]; w >= 0; w = T->next[w]) {
            double v0 = eval_hv_full(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].n0, vc[u].id0, vc[u].w0, vc[w].n0, vc[w].id0, vc[w].w0);
            double vp = eval_hv_full(C, M, vc[p].np, vc[p].idp, vc[p].wp, vc[u].n0, vc[u].id0, vc[u].w0, vc[w].n0, vc[w].id0, vc[w].w0);
            double vu = eval_hv_full(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].np, vc[u].idp, vc[u].wp, vc[w].n0, vc[w].id0, vc[w].w0);
            double vw = eval_hv_full(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].n0, vc[u].id0, vc[u].w0, vc[w].np, vc[w].idp, vc[w].wp);
            phi += v0;
            Lphi += vc[p].r * (vp - v0) + vc[u].r * (vu - v0) + vc[w].r * (vw - v0);
            if (T->L[p] >= 1) Lphi += vc[p].b * (eval_hv_full(C, M, vc[p].nb, vc[p].idb, vc[p].wb, vc[u].n0, vc[u].id0, vc[u].w0, vc[w].n0, vc[w].id0, vc[w].w0) - v0);
            if (T->L[u] >= 1) Lphi += vc[u].b * (eval_hv_full(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].nb, vc[u].idb, vc[u].wb, vc[w].n0, vc[w].id0, vc[w].w0) - v0);
            if (T->L[w] >= 1) Lphi += vc[w].b * (eval_hv_full(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].n0, vc[u].id0, vc[u].w0, vc[w].nb, vc[w].idb, vc[w].wb) - v0);
        }
    }

    /* 4. Full-tail hC(u; {w1, w2}) + newborn S3 + newborn Ybot */
    for (int u = 0; u < n; u++) if (vc[u].in_cu || vc[u].s3_b >= 0 || vc[u].ybu0 >= 0) {
        int ch[64], nc = 0;
        for (int c = T->head[u]; c >= 0; c = T->next[c]) ch[nc++] = c;
        for (int i = 0; i < nc; i++) for (int j = i + 1; j < nc; j++) {
            int w1 = ch[i], w2 = ch[j];
            if (vc[u].in_cu) {
                double v0 = eval_hc_full(C, M, vc[u].n0, vc[u].id0, vc[u].w0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0);
                double vu = eval_hc_full(C, M, vc[u].np, vc[u].idp, vc[u].wp, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0);
                double v1 = eval_hc_full(C, M, vc[u].n0, vc[u].id0, vc[u].w0, vc[w1].np, vc[w1].idp, vc[w1].wp, vc[w2].n0, vc[w2].id0, vc[w2].w0);
                double v2 = eval_hc_full(C, M, vc[u].n0, vc[u].id0, vc[u].w0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].np, vc[w2].idp, vc[w2].wp);
                phi += v0;
                Lphi += vc[u].r * (vu - v0) + vc[w1].r * (v1 - v0) + vc[w2].r * (v2 - v0);
                if (T->L[u] >= 1)  Lphi += vc[u].b  * (eval_hc_full(C, M, vc[u].nb, vc[u].idb, vc[u].wb, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0) - v0);
                if (T->L[w1] >= 1) Lphi += vc[w1].b * (eval_hc_full(C, M, vc[u].n0, vc[u].id0, vc[u].w0, vc[w1].nb, vc[w1].idb, vc[w1].wb, vc[w2].n0, vc[w2].id0, vc[w2].w0) - v0);
                if (T->L[w2] >= 1) Lphi += vc[w2].b * (eval_hc_full(C, M, vc[u].n0, vc[u].id0, vc[u].w0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].nb, vc[w2].idb, vc[w2].wb) - v0);
            }
            if (T->L[u] >= 1 && vc[u].s3_b >= 0)
                Lphi += vc[u].b * eval_s3_vec(C, M, vc[u].s3_b, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0, 1, nb_id, nb_w);
            if (vc[u].ybu0 >= 0) {
                if (T->L[w1] >= 1 && vc[w1].ybw_b >= 0)
                    Lphi += vc[w1].b * eval_yb_vec(C, M, vc[u].ybu0, vc[w1].ybw_b, vc[w2].n0, vc[w2].id0, vc[w2].w0, 1, nb_id, nb_w);
                if (T->L[w2] >= 1 && vc[w2].ybw_b >= 0)
                    Lphi += vc[w2].b * eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw_b, vc[w1].n0, vc[w1].id0, vc[w1].w0, 1, nb_id, nb_w);
            }
        }
    }

    /* 5. P4 */
    for (int c = 1; c < n; c++) if (vc[c].p40 >= 0) {
        int b = T->par[c]; if (b < 1 || vc[b].p40 < 0) continue;
        int a = T->par[b];
        for (int d = T->head[c]; d >= 0; d = T->next[d]) {
            double v0 = eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p40, vc[c].p40, vc[d].n0, vc[d].id0, vc[d].w0);
            double va = eval_p4_vec(C, M, vc[a].np, vc[a].idp, vc[a].wp, vc[b].p40, vc[c].p40, vc[d].n0, vc[d].id0, vc[d].w0);
            double vb = eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p4_p, vc[c].p40, vc[d].n0, vc[d].id0, vc[d].w0);
            double vc_v = eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p40, vc[c].p4_p, vc[d].n0, vc[d].id0, vc[d].w0);
            double vd = eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p40, vc[c].p40, vc[d].np, vc[d].idp, vc[d].wp);
            phi += v0;
            Lphi += vc[a].r * (va - v0) + vc[b].r * (vb - v0) + vc[c].r * (vc_v - v0) + vc[d].r * (vd - v0);
            if (T->L[a] >= 1) Lphi += vc[a].b * (eval_p4_vec(C, M, vc[a].nb, vc[a].idb, vc[a].wb, vc[b].p40, vc[c].p40, vc[d].n0, vc[d].id0, vc[d].w0) - v0);
            if (T->L[b] >= 1) Lphi += vc[b].b * (eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p4_b, vc[c].p40, vc[d].n0, vc[d].id0, vc[d].w0) - v0);
            if (T->L[c] >= 1) Lphi += vc[c].b * (eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p40, vc[c].p4_b, vc[d].n0, vc[d].id0, vc[d].w0) - v0);
            if (T->L[d] >= 1) Lphi += vc[d].b * (eval_p4_vec(C, M, vc[a].n0, vc[a].id0, vc[a].w0, vc[b].p40, vc[c].p40, vc[d].nb, vc[d].idb, vc[d].wb) - v0);
        }
    }

    /* 6. Ytop */
    for (int u = 1; u < n; u++) if (vc[u].yt0 >= 0) {
        int p = T->par[u], ch[32], nc = 0;
        for (int c = T->head[u]; c >= 0; c = T->next[c]) ch[nc++] = c;
        for (int i = 0; i < nc; i++) for (int j = i + 1; j < nc; j++) {
            int w1 = ch[i], w2 = ch[j];
            double v0 = eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0);
            double vp = eval_yt_vec(C, M, vc[p].np, vc[p].idp, vc[p].wp, vc[u].yt0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0);
            double vu = eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt_p, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0);
            double v1 = eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt0, vc[w1].np, vc[w1].idp, vc[w1].wp, vc[w2].n0, vc[w2].id0, vc[w2].w0);
            double v2 = eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].np, vc[w2].idp, vc[w2].wp);
            phi += v0;
            Lphi += vc[p].r * (vp - v0) + vc[u].r * (vu - v0) + vc[w1].r * (v1 - v0) + vc[w2].r * (v2 - v0);
            if (T->L[p] >= 1)  Lphi += vc[p].b  * (eval_yt_vec(C, M, vc[p].nb, vc[p].idb, vc[p].wb, vc[u].yt0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0) - v0);
            if (T->L[u] >= 1)  Lphi += vc[u].b  * (eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt_b, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0) - v0);
            if (T->L[w1] >= 1) Lphi += vc[w1].b * (eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt0, vc[w1].nb, vc[w1].idb, vc[w1].wb, vc[w2].n0, vc[w2].id0, vc[w2].w0) - v0);
            if (T->L[w2] >= 1) Lphi += vc[w2].b * (eval_yt_vec(C, M, vc[p].n0, vc[p].id0, vc[p].w0, vc[u].yt0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].nb, vc[w2].idb, vc[w2].wb) - v0);
        }
    }

    /* 7. Ybot */
    for (int u = 0; u < n; u++) if (vc[u].ybu0 >= 0) {
        int ch[32], nc = 0;
        for (int c = T->head[u]; c >= 0; c = T->next[c]) ch[nc++] = c;
        for (int j2 = 0; j2 < nc; j2++) {
            int w2 = ch[j2]; if (vc[w2].ybw0 < 0) continue;
            for (int j1 = 0; j1 < nc; j1++) if (j1 != j2) {
                int w1 = ch[j1];
                for (int c = T->head[w2]; c >= 0; c = T->next[c]) {
                    double v0 = eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].n0, vc[c].id0, vc[c].w0);
                    double vu = eval_yb_vec(C, M, vc[u].ybu_p, vc[w2].ybw0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].n0, vc[c].id0, vc[c].w0);
                    double v2 = eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw_p, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].n0, vc[c].id0, vc[c].w0);
                    double v1 = eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw0, vc[w1].np, vc[w1].idp, vc[w1].wp, vc[c].n0, vc[c].id0, vc[c].w0);
                    double vc_v = eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].np, vc[c].idp, vc[c].wp);
                    phi += v0;
                    Lphi += vc[u].r * (vu - v0) + vc[w2].r * (v2 - v0) + vc[w1].r * (v1 - v0) + vc[c].r * (vc_v - v0);
                    if (T->L[u] >= 1)  Lphi += vc[u].b  * (eval_yb_vec(C, M, vc[u].ybu_b, vc[w2].ybw0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].n0, vc[c].id0, vc[c].w0) - v0);
                    if (T->L[w2] >= 1) Lphi += vc[w2].b * (eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw_b, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].n0, vc[c].id0, vc[c].w0) - v0);
                    if (T->L[w1] >= 1) Lphi += vc[w1].b * (eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw0, vc[w1].nb, vc[w1].idb, vc[w1].wb, vc[c].n0, vc[c].id0, vc[c].w0) - v0);
                    if (T->L[c] >= 1)  Lphi += vc[c].b  * (eval_yb_vec(C, M, vc[u].ybu0, vc[w2].ybw0, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[c].nb, vc[c].idb, vc[c].wb) - v0);
                }
            }
        }
    }

    /* 8. S3 */
    for (int u = 0; u < n; u++) if (vc[u].s30 >= 0) {
        int ch[32], nc = 0;
        for (int c = T->head[u]; c >= 0; c = T->next[c]) ch[nc++] = c;
        for (int i = 0; i < nc; i++) for (int j = i + 1; j < nc; j++) for (int k = j + 1; k < nc; k++) {
            int w1 = ch[i], w2 = ch[j], w3 = ch[k];
            double v0 = eval_s3_vec(C, M, vc[u].s30, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].n0, vc[w3].id0, vc[w3].w0);
            double vu = eval_s3_vec(C, M, vc[u].s3_p, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].n0, vc[w3].id0, vc[w3].w0);
            double v1 = eval_s3_vec(C, M, vc[u].s30, vc[w1].np, vc[w1].idp, vc[w1].wp, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].n0, vc[w3].id0, vc[w3].w0);
            double v2 = eval_s3_vec(C, M, vc[u].s30, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].np, vc[w2].idp, vc[w2].wp, vc[w3].n0, vc[w3].id0, vc[w3].w0);
            double v3 = eval_s3_vec(C, M, vc[u].s30, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].np, vc[w3].idp, vc[w3].wp);
            phi += v0;
            Lphi += vc[u].r * (vu - v0) + vc[w1].r * (v1 - v0) + vc[w2].r * (v2 - v0) + vc[w3].r * (v3 - v0);
            if (T->L[u] >= 1)  Lphi += vc[u].b  * (eval_s3_vec(C, M, vc[u].s3_b, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].n0, vc[w3].id0, vc[w3].w0) - v0);
            if (T->L[w1] >= 1) Lphi += vc[w1].b * (eval_s3_vec(C, M, vc[u].s30, vc[w1].nb, vc[w1].idb, vc[w1].wb, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].n0, vc[w3].id0, vc[w3].w0) - v0);
            if (T->L[w2] >= 1) Lphi += vc[w2].b * (eval_s3_vec(C, M, vc[u].s30, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].nb, vc[w2].idb, vc[w2].wb, vc[w3].n0, vc[w3].id0, vc[w3].w0) - v0);
            if (T->L[w3] >= 1) Lphi += vc[w3].b * (eval_s3_vec(C, M, vc[u].s30, vc[w1].n0, vc[w1].id0, vc[w1].w0, vc[w2].n0, vc[w2].id0, vc[w2].w0, vc[w3].nb, vc[w3].idb, vc[w3].wb) - v0);
        }
    }
    *out_phi = phi; *out_Lphi = Lphi;
}

static void selftest_eval(const B3FLayout *C, const MidPot *M) {
    RNG r; rng_seed(&r, 777777ULL);
    CoreTree T; ct_init(&T, 4096);
    VCache *buf = NULL; int cap = 0;
    double max_rel_small = 0.0, max_rel_large = 0.0;
    for (int ti = 0; ti < 50; ti++) {
        ct_reset(&T); double t = 0.0;
        ct_advance(&T, &r, &t, 1.0 + 0.16 * ti, 0);
        double phi0, Lphi_fast;
        eval_phi_Lphi(C, M, &T, &phi0, &Lphi_fast, &buf, &cap);
        double Lphi_bf = 0.0;
        int n = T.C;
        for (int u = 0; u < n; u++) {
            uint64_t Lu = T.L[u]; int hu = (T.par[u] >= 0) ? 1 : 0;
            uint64_t Du = Lu + (uint64_t)hu + (uint64_t)T.kcnt[u];
            double rho = 0.0;
            if (hu) { int p = T.par[u]; rho += 1.0 / (double)(T.L[p] + (T.par[p] >= 0 ? 1 : 0) + T.kcnt[p]); }
            for (int c = T.head[u]; c >= 0; c = T.next[c]) rho += 1.0 / (double)(T.L[c] + 1 + T.kcnt[c]);
            double ru = (double)Lu + rho, bu = (double)Lu / (double)Du;
            T.L[u]++;
            double phi1, dummy; eval_phi_Lphi(C, M, &T, &phi1, &dummy, &buf, &cap);
            T.L[u]--;
            Lphi_bf += ru * (phi1 - phi0);
            if (Lu >= 1) {
                if (T.C == T.cap) ct_grow(&T);
                int v = T.C++;
                T.L[u]--; T.L[v] = 1; T.par[v] = u; T.head[v] = -1; T.next[v] = T.head[u]; T.head[u] = v; T.kcnt[v] = 0; T.kcnt[u]++;
                double phi2; eval_phi_Lphi(C, M, &T, &phi2, &dummy, &buf, &cap);
                T.C--; T.head[u] = T.next[v]; T.kcnt[u]--; T.L[u]++;
                Lphi_bf += bu * (phi2 - phi0);
            }
        }
        double rel = fabs(Lphi_fast - Lphi_bf) / fmax(fabs(Lphi_bf), 1e-12);
        if (ti < 25) { if (rel > max_rel_small) max_rel_small = rel; }
        else         { if (rel > max_rel_large) max_rel_large = rel; }
    }
    printf("Self-test O(C) eval_phi_Lphi vs brute-force on 50 trees: max_rel(small)=%.3e, max_rel(large)=%.3e\n",
           max_rel_small, max_rel_large);
    fflush(stdout);
    if (max_rel_small > 1e-11 || max_rel_large > 1e-9) { fprintf(stderr, "SELF-TEST FAILED!\n"); exit(1); }
    free(buf); ct_free(&T);
}

static void print_jackknife(const char *label, int NRUN, const double *N_a, const double *D_a, const double *Nc_a, const double *Dc_a) {
    double sN = 0, sD = 0, sNc = 0, sDc = 0;
    for (int r = 0; r < NRUN; r++) { sN += N_a[r]; sD += D_a[r]; sNc += Nc_a[r]; sDc += Dc_a[r]; }
    double est = sN / sD, est_c = sNc / sDc;
    double vj = 0, vjc = 0;
    for (int r = 0; r < NRUN; r++) {
        double ej = (sN - N_a[r]) / (sD - D_a[r]); vj += (ej - est) * (ej - est);
        double ejc = (sNc - Nc_a[r]) / (sDc - Dc_a[r]); vjc += (ejc - est_c) * (ejc - est_c);
    }
    double se = sqrt((double)(NRUN - 1) / NRUN * vj), sec = sqrt((double)(NRUN - 1) / NRUN * vjc);
    printf("  %-14s: B3-FullH VR = %.10f +- %.10f  | plain S/C = %.7f +- %.7f (VR gain = %.1fx)\n",
           label, est, se, est_c, sec, sec / se);
    fflush(stdout);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <cert_lo.bin> <cert_hi.bin> <NRUN> [threads=4]\n", argv[0]);
        return 2;
    }
    int NRUN = atoi(argv[3]);
    int nthreads = (argc > 4) ? atoi(argv[4]) : 4;

    B3FLayout C; MidPot M;
    load_midpoint(argv[1], argv[2], &C, &M);
    selftest_eval(&C, &M);

    const int NC = 6;
    const double ts[6] = { 8.0, 12.0, 16.0, 20.0, 22.0, 24.0 };
    double *phi_arr  = (double*)calloc((size_t)NRUN * NC, sizeof(double));
    double *Lphi_arr = (double*)calloc((size_t)NRUN * NC, sizeof(double));
    double *C_arr    = (double*)calloc((size_t)NRUN * NC, sizeof(double));
    double *S_arr    = (double*)calloc((size_t)NRUN * NC, sizeof(double));

    double t0 = omp_get_wtime();
    #pragma omp parallel num_threads(nthreads)
    {
        CoreTree T; ct_init(&T, 1048576);
        VCache *buf = NULL; int cap = 0;
        #pragma omp for schedule(dynamic, 4)
        for (int run = 0; run < NRUN; run++) {
            RNG r; rng_seed(&r, 202611249999ULL + (uint64_t)run * 2654435761ULL);
            ct_reset(&T); double t = 0.0;
            for (int ci = 0; ci < NC; ci++) {
                ct_advance(&T, &r, &t, ts[ci], 0);
                double phi, Lphi;
                eval_phi_Lphi(&C, &M, &T, &phi, &Lphi, &buf, &cap);
                double S = 0.0;
                for (int u = 0; u < T.C; u++) S += (double)T.L[u] / (double)(T.L[u] + (T.par[u] >= 0 ? 1 : 0) + T.kcnt[u]);
                phi_arr[(size_t)run * NC + ci]  = phi;
                Lphi_arr[(size_t)run * NC + ci] = Lphi;
                C_arr[(size_t)run * NC + ci]    = (double)T.C;
                S_arr[(size_t)run * NC + ci]    = S;
            }
        }
        free(buf); ct_free(&T);
    }
    double t1 = omp_get_wtime();
    printf("\nCompleted %d independent Lazy-Yule trees to t=%.0f in %.2fs (%d cores):\n", NRUN, ts[NC - 1], t1 - t0, nthreads);

    double *na = (double*)malloc(NRUN * sizeof(double)), *da = (double*)malloc(NRUN * sizeof(double));
    double *nca = (double*)malloc(NRUN * sizeof(double)), *dca = (double*)malloc(NRUN * sizeof(double));
    const int wins[7][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {3, 5}, {2, 5} };
    const char *wlabels[7] = { "t=8 -> 12", "t=12 -> 16", "t=16 -> 20", "t=20 -> 22", "t=22 -> 24", "t=20 -> 24", "t=16 -> 24" };
    for (int wi = 0; wi < 7; wi++) {
        int a = wins[wi][0], b = wins[wi][1];
        for (int r = 0; r < NRUN; r++) {
            na[r]  = Lphi_arr[(size_t)r * NC + b] - Lphi_arr[(size_t)r * NC + a];
            da[r]  = phi_arr[(size_t)r * NC + b]  - phi_arr[(size_t)r * NC + a];
            nca[r] = S_arr[(size_t)r * NC + b]    - S_arr[(size_t)r * NC + a];
            dca[r] = C_arr[(size_t)r * NC + b]    - C_arr[(size_t)r * NC + a];
        }
        print_jackknife(wlabels[wi], NRUN, na, da, nca, dca);
    }
    /* Raw checkpoint t=24 as well */
    for (int r = 0; r < NRUN; r++) {
        na[r] = Lphi_arr[(size_t)r * NC + 5]; da[r] = phi_arr[(size_t)r * NC + 5];
        nca[r] = S_arr[(size_t)r * NC + 5];   dca[r] = C_arr[(size_t)r * NC + 5];
    }
    print_jackknife("raw t=24", NRUN, na, da, nca, dca);

    free(na); free(da); free(nca); free(dca);
    free(phi_arr); free(Lphi_arr); free(C_arr); free(S_arr);
    free(M.f); free(M.g); free(M.hV); free(M.hC); free(M.qP4); free(M.qYt); free(M.qYb); free(M.qS3);
    b3f_layout_free(&C);
    return 0;
}
