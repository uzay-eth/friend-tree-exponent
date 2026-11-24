/*
 * Deterministic Builder for B3-FullH + Level-4 Exact-Integer Drift Certificates
 *
 * Usage:
 *   cert_b3f_build --mode <1|2> --P <P> --Q <Q> --out <cert.bin> [options]
 *     --mode 1 (UPPER: L Phi <= (P/Q) Phi) | 2 (LOWER: L Phi >= (P/Q) Phi)
 *     --D0 16 --KT 15 --LF 12000 --KF 160
 *     --DP4 9 --KP4 6 --DYT 8 --KYT 5 --DYB 9 --KYB 6 --DS3 7 --KS3 5
 *     --nit 22 --threads 4 --S 100000000000 --Sr 1000000000
 *
 * Strategy:
 *   1. Deterministic xoshiro256** simulation generates empirical centering rates
 *      (rpar, rch, trpar, trch), quantized immediately to exact rationals / S_r
 *      so solver and verifier use IDENTICAL centering constants!
 *   2. Fixed-point iteration finds tight potential tables (f, g, hV, hC, P4, Yt, Yb, S3)
 *      and exact f-tail transfer coefficients (c2_f, c3_f) and T3 Y_ref tables.
 *   3. Final 3 sweeps run at the FIXED target lambda = P/Q with exact small per-level
 *      safety cushions (eps_L4 = 5e-10, eps_h = 2e-9, eps_g = 1e-8) so every quantized
 *      int64_t row satisfies the exact rational inequality with >= 50x margin over 1/S!
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include "cert_b3f_l5_format.h"
#include "l5_build_rows.h"

static inline double birth_value(const double *q, size_t i) { return i == SIZE_MAX ? 0 : q[i]; }

static int parse_i64_cli(const char *s, int64_t *out) {
    if (!s || !*s) return 0;
    errno = 0;
    char *end = NULL;
    long long v = strtoll(s, &end, 10);
    if (errno == ERANGE || !end || *end != '\0') return 0;
    *out = (int64_t)v;
    return 1;
}
static int parse_i32_cli(const char *s, int32_t *out) {
    int64_t v = 0;
    if (!parse_i64_cli(s, &v) || v < INT32_MIN || v > INT32_MAX) return 0;
    *out = (int32_t)v;
    return 1;
}
static int parse_dbl_cli(const char *s, double *out) {
    if (!s || !*s) return 0;
    errno = 0;
    char *end = NULL;
    double v = strtod(s, &end);
    if (errno == ERANGE || !end || *end != '\0' || !isfinite(v)) return 0;
    *out = v;
    return 1;
}

static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next64(void) {
    uint64_t r = rotl(xs[1] * 5, 7) * 9, t = xs[1] << 17;
    xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3];
    xs[2] ^= t; xs[3] = rotl(xs[3], 45);
    return r;
}
static uint64_t sm64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline uint32_t unif32(uint32_t m) {
    uint64_t x = next64() >> 32, prod = x * (uint64_t)m;
    uint32_t low = (uint32_t)prod;
    if (low < m) {
        uint32_t thr = (uint32_t)(-m) % m;
        while (low < thr) { x = next64() >> 32; prod = x * (uint64_t)m; low = (uint32_t)prod; }
    }
    return (uint32_t)(prod >> 32);
}

static inline int b3f_embed(const B3FLayout *C, int L, int k, int h, int ids[4], double wts[4]) {
    int D = L + k + h;
    if (D < C->D0) { ids[0] = b3f_inner(C, L, k, h); wts[0] = 1.0; return 1; }
    double t = (double)C->D0 / D, s = 1.0 - t;
    ids[0] = C->INF_ID;       wts[0] = s * s * s;
    ids[1] = b3f_t1(C, k, h); wts[1] = 3.0 * t * s * s;
    ids[2] = b3f_t2(C, k, h); wts[2] = 3.0 * t * t * s;
    ids[3] = b3f_t3(C, k, h); wts[3] = t * t * t;
    return 4;
}

static inline double invD_of(const B3FLayout *C, int u) {
    int k = C->kind[u];
    if (k == 0) return 1.0 / C->D[u];
    if (k == 1) return 0.0;
    if (k == 2) return 1.0 / (3.0 * C->D0);
    if (k == 3) return 2.0 / (3.0 * C->D0);
    return 1.0 / C->D0;
}

static inline void get_rin(const B3FLayout *C, int y, double *lo, double *hi) {
    if (C->kind[y] == 0) { *lo = *hi = 1.0 / C->D[y]; }
    else { *lo = 0.0; *hi = 1.0 / C->D0; }
}

static inline double eval_Y3_half(int D0, double F0, double F1, double F2, double F3) {
    double d1 = F0 - F1, d2 = F1 - F2, d3 = F2 - F3, dD0 = (double)D0, t = 0.5;
    double s = t / (dD0 + t);
    double da1 = -3.0 * (1.0 - t) * (dD0 + 1.0) / (dD0 * (dD0 + t));
    double da2 = +6.0 * (1.0 - t) / (dD0 + t);
    double da3 = 3.0 * s;
    double db0 = 3.0 * (1.0 - t) * s / (dD0 + t);
    double db1 = 3.0 * t * s / (dD0 + t);
    double dgm = s * t * t / ((dD0 + t) * (dD0 + t));
    return da1 * d1 + da2 * d2 + da3 * d3 + db0 * (d1 - d2) + db1 * (d2 - d3) + dgm * (d1 - 2.0 * d2 + d3);
}

static inline void b3_tail_const(int mode, int D0, int KT, int xk, int kx, int hx,
                                 double rho_min, double rho_max,
                                 double F0, double F1, double F1b, double F2b,
                                 int has_nbr, double d0_nbr, double d1_nbr,
                                 double *num0, double *den0) {
    if (xk == 2) {
        if (kx <= KT) { *num0 += F0 + F1b; *den0 += 2.0; }
        if (has_nbr) *num0 += d0_nbr / (3.0 * D0);
    } else if (xk == 3) {
        if (kx <= KT) {
            int c = kx + hx; double d1 = F0 - F1;
            double a_lo = -((double)(c + 1) - rho_min) / D0 * d1;
            double a_hi = -((double)(c + 1) - rho_max) / D0 * d1;
            *num0 += 2.0 * F1 + ((mode == 1) ? fmax(a_lo, a_hi) : fmin(a_lo, a_hi)) + F2b - ((double)c / D0) * (F1b - F1);
            *den0 += 3.0;
        }
        if (has_nbr) *num0 += (d1_nbr - d0_nbr / 3.0) / D0;
    }
}

static inline void b3_tail_var(int mode, int D0, int KT, int xk, int kx, int hx,
                               double r_in_lo, double r_in_hi, double rho_ext_max, double Rbar_ext, double Y_ref,
                               double F0, double F1, double F2, double F1b, double F2b, double F3b,
                               int n_nbr_in, double d0_nbr, double d1_nbr, double d2_nbr, double F3_up_sum,
                               double v, double *num, double *den, int *mb, int bit_shift) {
    if (xk == 2) {
        if (kx > KT && ((mode == 1 && F0 > v) || (mode == 2 && F0 < v))) {
            *num += F0; *den += 1.0; *mb |= (1 << bit_shift);
        }
    } else if (xk == 3) {
        if (kx > KT && ((mode == 1 && F1 > v) || (mode == 2 && F1 < v))) {
            *num += 2.0 * F1; *den += 2.0; *mb |= (1 << bit_shift);
        }
    } else if (xk == 4) {
        double d1 = F0 - F1, d2 = F1 - F2;
        if (kx <= KT) {
            int c = kx + hx;
            double dD0 = (double)D0, D02 = dD0 * dD0, D03 = D02 * dD0;
            double tmax = (c > D0) ? (dD0 / c) : 1.0, tm2 = tmax * tmax, tm3 = tm2 * tmax, tm4 = tm2 * tm2;
            double b0 = D03, b1 = 3 * D02, b2 = 3 * dD0, b3 = 1.0;
            double BB[5] = {
                b0,
                b0 + b1 * tmax * 0.25,
                b0 + b1 * tmax * 0.5 + b2 * tm2 / 6.0,
                b0 + b1 * tmax * 0.75 + b2 * tm2 * 0.5 + b3 * tm3 * 0.25,
                b0 + b1 * tmax + b2 * tm2 + b3 * tm3
            };
            double db = 3.0 * ((double)c / dD0) * ((F1b - F1) - (F2b - F2));
            double a0 = F3b + db, a1 = -((double)c / dD0) * F3b - db, e0 = 1.0, e1 = -(double)c / dD0;
            if (n_nbr_in > 0) {
                double nd0 = (d0_nbr - 3 * d1_nbr + 3 * d2_nbr) / dD0;
                a0 += nd0;
                a1 += (F3_up_sum / dD0) - nd0;
                e1 += (double)n_nbr_in / dD0;
            }
            double an[5] = { a0 * D03, (a1 * dD0 + 3 * a0) * D02, (3 * a1 * dD0 + 3 * a0) * dD0, 3 * a1 * dD0 + a0, a1 };
            double ad[5] = { e0 * D03, (e1 * dD0 + 3 * e0) * D02, (3 * e1 * dD0 + 3 * e0) * dD0, 3 * e1 * dD0 + e0, e1 };

            double rins[2] = { r_in_lo, r_in_hi };
            int nrin = (fabs(r_in_hi - r_in_lo) < 1e-15) ? 1 : 2;
            double rexts[2] = { 0.0, rho_ext_max };
            int nrext = (rho_ext_max > 1e-15) ? 2 : 1;

            double best_val = (mode == 1) ? -1e30 : 1e30, best_n = 0, best_d = 0;
            int best_id = 0;
            for (int ri = 0; ri < nrin; ri++) for (int re = 0; re < nrext; re++) {
                double rho = rins[ri] + rexts[re];
                double shift = (rexts[re] - Rbar_ext) * Y_ref;
                double x = (double)(c + 1) - rho;
                double n0 = an[0] + 3 * D03 * F2 + d1 * (3 * D02 * x + 3 * D02 + 3 * dD0 * x) + d2 * (-6 * D02 * x - 3 * D02);
                double n1 = an[1] + F2 * (6 * D02 - 3 * D02 * x) + d1 * (-3 * D02 * x - 3 * D02 + 6 * dD0 + 6 * x) + d2 * (6 * D02 * x + 6 * D02 - 9 * dD0 * x - 6 * dD0);
                double n2 = an[2] + F2 * (4 * dD0 - 3 * dD0 * x) + d1 * (-3 * dD0 * x - 5 * dD0 - 6 * x + 3 + 3 * x / dD0) + d2 * (6 * dD0 * x + 10 * dD0 - 3 * x - 3);
                double n3 = an[3] + F2 * (1 - x) + d1 * (-x - 2 - 3 * x / dD0) + d2 * (2 * x + 4);
                double n4 = an[4];
                double d0c = ad[0] + 3 * D03, d1c = ad[1] + 3 * D02 * (2 - x), d2c = ad[2] + dD0 * (4 - 3 * x), d3c = ad[3] + (1 - x), d4c = ad[4];
                double BN[5] = {
                    n0,
                    n0 + n1 * tmax * 0.25,
                    n0 + n1 * tmax * 0.5 + n2 * tm2 / 6.0,
                    n0 + n1 * tmax * 0.75 + n2 * tm2 * 0.5 + n3 * tm3 * 0.25,
                    n0 + n1 * tmax + n2 * tm2 + n3 * tm3 + n4 * tm4
                };
                double BD[5] = {
                    d0c,
                    d0c + d1c * tmax * 0.25,
                    d0c + d1c * tmax * 0.5 + d2c * tm2 / 6.0,
                    d0c + d1c * tmax * 0.75 + d2c * tm2 * 0.5 + d3c * tm3 * 0.25,
                    d0c + d1c * tmax + d2c * tm2 + d3c * tm3 + d4c * tm4
                };
                for (int j = 0; j < 5; j++) {
                    double rn = BN[j] / BB[j] - shift, rd = BD[j] / BB[j], val = rn - rd * v;
                    if ((mode == 1 && val > best_val) || (mode == 2 && val < best_val)) {
                        best_val = val; best_n = rn; best_d = rd; best_id = ((ri * nrext + re) * 5 + j) + 1;
                    }
                }
            }
            *num += best_n; *den += best_d; *mb ^= (best_id << bit_shift);
        } else {
            double dD0 = (double)D0, D02 = dD0 * dD0, BB2[3] = { D02, D02 + dD0, D02 + 2 * dD0 + 1.0 };
            double n0 = 3 * D02 * F2 + 3 * dD0 * d1 - 3 * dD0 * d2;
            double n1 = 3 * dD0 * F2 + d1 * (3 - 3 * dD0) + d2 * (6 * dD0 - 3);
            double n2 = F2 - 2 * d1 + 4 * d2;
            double d0c = 3 * D02, d1c = 3 * dD0, d2c = 1.0;
            double BN2[3] = { n0, n0 + 0.5 * n1, n0 + n1 + n2 };
            double BD2[3] = { d0c, d0c + 0.5 * d1c, d0c + d1c + d2c };
            double best_val = 0, best_n = 0, best_d = 0; int bid = 0;
            for (int j = 0; j < 3; j++) {
                double rn = BN2[j] / BB2[j], rd = BD2[j] / BB2[j], val = rn - rd * v;
                if ((mode == 1 && val > best_val) || (mode == 2 && val < best_val)) {
                    best_val = val; best_n = rn; best_d = rd; bid = j + 1;
                }
            }
            *num += best_n; *den += best_d; *mb ^= (bid << bit_shift);
            if (n_nbr_in > 0) {
                double c0 = (d0_nbr - 3 * d1_nbr + 3 * d2_nbr) / dD0, c1 = (F3_up_sum - n_nbr_in * v) / dD0;
                if ((mode == 1 && c1 > c0) || (mode == 2 && c1 < c0)) {
                    *num += F3_up_sum / dD0; *den += (double)n_nbr_in / dD0; *mb ^= (16 << bit_shift);
                } else {
                    *num += c0;
                }
            }
        }
    }
}

static void run_empirical(const B3FLayout *C, int NMAX, int RUNS, int64_t S_r,
                          int64_t *rpar_q, int64_t *rch_q, int64_t *trpar_q, int64_t *trch_q) {
    double *ep = (double*)b3f_xcalloc(C->INNER_N, sizeof(double)), *ec = (double*)b3f_xcalloc(C->INNER_N, sizeof(double));
    double *cp = (double*)b3f_xcalloc(C->INNER_N, sizeof(double)), *cc = (double*)b3f_xcalloc(C->INNER_N, sizeof(double));
    double *tep = (double*)b3f_xcalloc(C->NT, sizeof(double)), *tec = (double*)b3f_xcalloc(C->NT, sizeof(double));
    double *tcp = (double*)b3f_xcalloc(C->NT, sizeof(double)), *tcc = (double*)b3f_xcalloc(C->NT, sizeof(double));

    uint64_t seed = 20261124ULL;
    for (int i = 0; i < 4; i++) xs[i] = sm64(&seed);

    uint32_t *par = (uint32_t*)b3f_xmalloc(NMAX * sizeof(uint32_t));
    uint32_t *deg = (uint32_t*)b3f_xmalloc(NMAX * sizeof(uint32_t));
    uint32_t *kcnt = (uint32_t*)b3f_xcalloc(NMAX, sizeof(uint32_t));
    uint32_t *kcap = (uint32_t*)b3f_xcalloc(NMAX, sizeof(uint32_t));
    uint32_t **kids = (uint32_t**)b3f_xcalloc(NMAX, sizeof(uint32_t*));
    int *st = (int*)b3f_xmalloc(NMAX * sizeof(int)), *tk = (int*)b3f_xmalloc(NMAX * sizeof(int));

    for (int run = 0; run < RUNS; run++) {
        for (int i = 0; i < NMAX; i++) {
            kcnt[i] = 0;
            if (kcap[i] > 64) { free(kids[i]); kids[i] = NULL; kcap[i] = 0; }
        }
        par[0] = UINT32_MAX; deg[0] = 2; par[1] = 0; deg[1] = 1; par[2] = 0; deg[2] = 1;
        if (kcap[0] < 2) { kids[0] = (uint32_t*)b3f_xrealloc(kids[0], 2 * sizeof(uint32_t)); kcap[0] = 2; }
        kids[0][0] = 1; kids[0][1] = 2; kcnt[0] = 2;
        for (uint32_t n = 3; n < (uint32_t)NMAX; n++) {
            uint32_t V = unif32(n), d = deg[V], r = unif32(d);
            uint32_t W = (par[V] != UINT32_MAX) ? ((r == 0) ? par[V] : kids[V][r - 1]) : kids[V][r];
            par[n] = W; deg[n] = 1; deg[W]++;
            if (kcnt[W] == kcap[W]) {
                uint32_t nc = kcap[W] ? 2 * kcap[W] : 2;
                kids[W] = (uint32_t*)b3f_xrealloc(kids[W], nc * sizeof(uint32_t));
                kcap[W] = nc;
            }
            kids[W][kcnt[W]++] = n;
        }
        for (int u = 0; u < NMAX; u++) {
            if (deg[u] < 2) { st[u] = -1; tk[u] = -1; continue; }
            int L = 0;
            for (uint32_t i = 0; i < kcnt[u]; i++) if (deg[kids[u][i]] == 1) L++;
            int h = (par[u] != UINT32_MAX) ? 1 : 0, k = (int)deg[u] - L - h;
            if ((int)deg[u] < C->D0) { st[u] = b3f_inner(C, L, k, h); tk[u] = -1; }
            else { st[u] = -1; int kap = (k > C->KT) ? (C->KT + 1) : k; tk[u] = b3f_tidx(C, kap, h); }
        }
        for (int w = 1; w < NMAX; w++) {
            if (deg[w] < 2) continue;
            int p = (int)par[w];
            if (st[w] >= 0) { ep[st[w]] += 1.0 / deg[p]; cp[st[w]] += 1.0; }
            if (st[p] >= 0) { ec[st[p]] += 1.0 / deg[w]; cc[st[p]] += 1.0; }
            if (tk[w] >= 0) { tep[tk[w]] += 1.0 / deg[p]; tcp[tk[w]] += 1.0; }
            if (tk[p] >= 0) { tec[tk[p]] += 1.0 / deg[w]; tcc[tk[p]] += 1.0; }
        }
    }
    for (int u = 0; u < C->INNER_N; u++) {
        int k = C->K[u], L = C->L[u];
        double vp = (cp[u] >= 5.0) ? (ep[u] / cp[u]) : (C->H[u] ? (0.25 / (1.0 + 0.15 * L + 0.35 * k)) : 0.0);
        double vc = (cc[u] >= 5.0) ? (ec[u] / cc[u]) : (k > 0 ? (0.36 / (1.0 + 0.04 * L + 0.08 * k)) : 0.0);
        rpar_q[u] = (int64_t)llround(vp * (double)S_r);
        rch_q[u]  = (int64_t)llround(vc * (double)S_r);
    }
    for (int h = 1; h >= 0; h--) for (int kap = 0; kap <= C->KT + 1; kap++) {
        int idx = b3f_tidx(C, kap, h);
        double vp = (tcp[idx] >= 5.0) ? (tep[idx] / tcp[idx]) : (h ? 0.05 : 0.0);
        double vc = (tcc[idx] >= 5.0) ? (tec[idx] / tcc[idx]) : (kap > 0 ? 0.26 : 0.0);
        trpar_q[idx] = (int64_t)llround(vp * (double)S_r);
        trch_q[idx]  = (int64_t)llround(vc * (double)S_r);
    }
    for (int i = 0; i < NMAX; i++) free(kids[i]);
    free(par); free(deg); free(kcnt); free(kcap); free(kids); free(st); free(tk);
    free(ep); free(ec); free(cp); free(cc); free(tep); free(tec); free(tcp); free(tcc);
}

static int64_t quantize(double value, int64_t scale) {
    double v = value*(double)scale;
    if (!isfinite(v) || fabs(v) > CERT_B3F_MAX_VAL) { fputs("Quantization out of range\n",stderr); exit(1); }
    return llround(v);
}

static void write_i64_array(FILE *fp, const double *arr, size_t n, int64_t S) {
    const size_t CHUNK = 1 << 20;
    int64_t *buf = (int64_t*)b3f_xmalloc(CHUNK * sizeof(int64_t));
    for (size_t off = 0; off < n; off += CHUNK) {
        size_t cnt = (n - off < CHUNK) ? (n - off) : CHUNK;
        for (size_t i = 0; i < cnt; i++) buf[i] = quantize(arr[off+i],S);
        if (fwrite(buf, sizeof(int64_t), cnt, fp) != cnt) { fprintf(stderr, "fwrite failed\n"); exit(1); }
    }
    free(buf);
}

int main(int argc, char **argv) {
    CertB3FHeader H;
    memset(&H, 0, sizeof(H));
    H.magic   = CERT_B3F_MAGIC;
    H.version = 2;
    L5Header E = {6,2,3,0,0};
    H.mode    = 1;
    H.P       = 0;
    H.Q       = 1000000;
    H.S       = 100000000000LL; /* 1e11 */
    H.S_r     = 1000000000LL;   /* 1e9 */
    H.LF      = 12000; H.KF = 160;
    H.D0      = 16;    H.KT = 15;
    H.DP4     = 9;     H.KP4 = 6;
    H.DYT     = 8;     H.KYT = 5;
    H.DYB     = 9;     H.KYB = 6;
    H.DS3     = 7;     H.KS3 = 5;

    const char *out_path = NULL;
    int nit = 22, nthreads = 1, sizing = 0;
    double cushion = 2.0e-6;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--mode") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.mode)) { fprintf(stderr, "Invalid --mode\n"); return 2; } }
        else if (!strcmp(argv[i], "--P") && i + 1 < argc) { if (!parse_i64_cli(argv[++i], &H.P)) { fprintf(stderr, "Invalid --P\n"); return 2; } }
        else if (!strcmp(argv[i], "--Q") && i + 1 < argc) { if (!parse_i64_cli(argv[++i], &H.Q)) { fprintf(stderr, "Invalid --Q\n"); return 2; } }
        else if (!strcmp(argv[i], "--S") && i + 1 < argc) { if (!parse_i64_cli(argv[++i], &H.S)) { fprintf(stderr, "Invalid --S\n"); return 2; } }
        else if (!strcmp(argv[i], "--Sr") && i + 1 < argc) { if (!parse_i64_cli(argv[++i], &H.S_r)) { fprintf(stderr, "Invalid --Sr\n"); return 2; } }
        else if (!strcmp(argv[i], "--LF") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.LF)) { fprintf(stderr, "Invalid --LF\n"); return 2; } }
        else if (!strcmp(argv[i], "--KF") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.KF)) { fprintf(stderr, "Invalid --KF\n"); return 2; } }
        else if (!strcmp(argv[i], "--D0") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.D0)) { fprintf(stderr, "Invalid --D0\n"); return 2; } }
        else if (!strcmp(argv[i], "--KT") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.KT)) { fprintf(stderr, "Invalid --KT\n"); return 2; } }
        else if (!strcmp(argv[i], "--DP4") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.DP4)) { fprintf(stderr, "Invalid --DP4\n"); return 2; } }
        else if (!strcmp(argv[i], "--KP4") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.KP4)) { fprintf(stderr, "Invalid --KP4\n"); return 2; } }
        else if (!strcmp(argv[i], "--DYT") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.DYT)) { fprintf(stderr, "Invalid --DYT\n"); return 2; } }
        else if (!strcmp(argv[i], "--KYT") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.KYT)) { fprintf(stderr, "Invalid --KYT\n"); return 2; } }
        else if (!strcmp(argv[i], "--DYB") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.DYB)) { fprintf(stderr, "Invalid --DYB\n"); return 2; } }
        else if (!strcmp(argv[i], "--KYB") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.KYB)) { fprintf(stderr, "Invalid --KYB\n"); return 2; } }
        else if (!strcmp(argv[i], "--DS3") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.DS3)) { fprintf(stderr, "Invalid --DS3\n"); return 2; } }
        else if (!strcmp(argv[i], "--KS3") && i + 1 < argc) { if (!parse_i32_cli(argv[++i], &H.KS3)) { fprintf(stderr, "Invalid --KS3\n"); return 2; } }
        else if (!strcmp(argv[i], "--nit") && i + 1 < argc) { int32_t v; if (!parse_i32_cli(argv[++i], &v) || v < 1 || v > 1000) { fprintf(stderr, "Invalid --nit\n"); return 2; } nit = v; }
        else if (!strcmp(argv[i], "--cushion") && i + 1 < argc) { if (!parse_dbl_cli(argv[++i], &cushion) || cushion < 0.0 || cushion > 0.05) { fprintf(stderr, "Invalid --cushion\n"); return 2; } }
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) { int32_t v; if (!parse_i32_cli(argv[++i], &v) || v < 1 || v > 256) { fprintf(stderr, "Invalid --threads\n"); return 2; } nthreads = v; }
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--size")) sizing = 1;
        else if ((!strcmp(argv[i], "--Dc") || !strcmp(argv[i], "--buckets") || !strcmp(argv[i], "--l5mask")) && i + 1 < argc) {
            const char *key = argv[i]; int32_t v;
            if (!parse_i32_cli(argv[++i],&v) || v < 0) { fprintf(stderr,"Invalid %s\n",key); return 2; }
            if (!strcmp(key,"--Dc")) E.dc = (uint32_t)v;
            else if (!strcmp(key,"--buckets")) E.buckets = (uint32_t)v;
            else E.mask = (uint32_t)v;
        }
        else { fprintf(stderr, "Unknown arg: %s\n", argv[i]); return 2; }
    }
    if (!out_path && !sizing) { fprintf(stderr, "Missing --out <path>\n"); return 2; }
    int auto_pq = (H.P == 0);
    if (auto_pq) H.P = (H.mode == 1) ? 570000 : 560000;
    int hv = b3f_l5_validate_header(&H,&E);
    if (hv != 0) { fprintf(stderr, "Invalid header parameters: err=%d\n", hv); return 2; }

    B3FLayout C;
    b3f_layout_init(&C, &H);
    L5 *l5 = l5_init(&C,&E);
    E.rows = l5->rows;
    size_t bytes = l5_expected_bytes(l5);
    if (bytes == SIZE_MAX || bytes > LONG_MAX || !E.rows) { fputs("Invalid L5 payload size\n",stderr); return 2; }
    printf("L5 Dc=%u buckets=%u mask=%u rows=%zu file_bytes=%zu layout_bytes=%zu per_worker_cache_bytes=%zu\n",
           E.dc,E.buckets,E.mask,l5->rows,bytes,sizeof *l5,L5_CACHE*sizeof(L5FloatCache));
    if (sizing) {
        size_t lookup = sizeof C+sizeof *l5+sizeof(int)*((size_t)C.NG*
            (13+4*(C.max_rank+1)+2*l5->buckets) + 4*C.DP4*C.KP4 + 5*C.DYT*C.KYT +
            14*C.DYB*C.KYB + 10*C.DS3*C.KS3);
        size_t p4_stages = (size_t)C.NP4*C.NP4, yb_stages = (size_t)C.NYBU*C.NYBW;
        size_t stages = 2*sizeof(int)*(p4_stages > yb_stages ? p4_stages : yb_stages);
        size_t scratch = sizeof(double)*((size_t)C.NV+3*C.INNER_N+2*C.NT);
        size_t l5_cache = (size_t)nthreads*L5_CACHE*sizeof(L5FloatCache);
        printf("sizing lookup_bytes=%zu builder_sweep_bytes=%zu verifier_l5_bytes=%zu\n",lookup,
               bytes-120+lookup+scratch+(stages > l5_cache ? stages : l5_cache),
               bytes-120+lookup+sizeof(int64_t)*C.INNER_N+(size_t)nthreads*L5_CACHE*2*sizeof(int64_t));
        l5_free(l5); b3f_layout_free(&C); return 0;
    }
    l5_init_fibers(l5);
    double *q5 = b3f_xcalloc(l5->rows,sizeof *q5);
    printf("B3-FullH builder: mode=%d D0=%d KT=%d LF=%d KF=%d P4=(%d,%d) Yt=(%d,%d) Yb=(%d,%d) S3=(%d,%d)\n",
           H.mode, C.D0, C.KT, C.LF, C.KF, C.DP4, C.KP4, C.DYT, C.KYT, C.DYB, C.KYB, C.DS3, C.KS3);
    printf("  npars=%d nchs=%d NHU=%d NCU=%d | tot_entries=%.2fM | expected_bytes=%.2f MB\n",
           C.npars, C.nchs, C.NHU, C.NCU,
           (C.NV + C.g_len + C.hv_len + C.hc_len + C.p4_len + C.yt_len + C.yb_len + C.s3_len) / 1e6,
           b3f_expected_bytes(&C) / (1024.0 * 1024.0));
    fflush(stdout);

    int64_t *rpar_q  = (int64_t*)b3f_xcalloc(C.INNER_N, sizeof(int64_t));
    int64_t *rch_q   = (int64_t*)b3f_xcalloc(C.INNER_N, sizeof(int64_t));
    int64_t *trpar_q = (int64_t*)b3f_xcalloc(C.NT, sizeof(int64_t));
    int64_t *trch_q  = (int64_t*)b3f_xcalloc(C.NT, sizeof(int64_t));
    run_empirical(&C, 500000, 35, H.S_r, rpar_q, rch_q, trpar_q, trch_q);

    double *rpar  = (double*)b3f_xcalloc(C.INNER_N, sizeof(double));
    double *rch   = (double*)b3f_xcalloc(C.INNER_N, sizeof(double));
    double *trpar = (double*)b3f_xcalloc(C.NT, sizeof(double));
    double *trch  = (double*)b3f_xcalloc(C.NT, sizeof(double));
    for (int i = 0; i < C.INNER_N; i++) { rpar[i] = (double)rpar_q[i] / (double)H.S_r; rch[i] = (double)rch_q[i] / (double)H.S_r; }
    for (int i = 0; i < C.NT; i++)      { trpar[i] = (double)trpar_q[i] / (double)H.S_r; trch[i] = (double)trch_q[i] / (double)H.S_r; }

    int N = C.NG, D0 = C.D0, KT = C.KT, NT = C.NT, mode = H.mode;
    int v_nb = b3f_vidx(&C, 1, 0, 1), iw_nb = b3f_inner(&C, 1, 0, 1);

    double *f    = (double*)b3f_xcalloc(C.NV, sizeof(double)), *fn = (double*)b3f_xcalloc(C.NV, sizeof(double));
    double *g    = (double*)b3f_xcalloc(C.g_len, sizeof(double));
    double *hV   = (double*)b3f_xcalloc(C.hv_len, sizeof(double));
    double *hC   = (double*)b3f_xcalloc(C.hc_len, sizeof(double));
    double *qP4  = (double*)b3f_xcalloc(C.p4_len, sizeof(double));
    double *qYt  = (double*)b3f_xcalloc(C.yt_len, sizeof(double));
    double *qYb  = (double*)b3f_xcalloc(C.yb_len, sizeof(double));
    double *qS3  = (double*)b3f_xcalloc(C.s3_len, sizeof(double));
    double *dLf  = (double*)b3f_xcalloc(C.INNER_N, sizeof(double));
    double *c2_f = (double*)b3f_xcalloc(NT, sizeof(double)), *c3_f = (double*)b3f_xcalloc(NT, sizeof(double));
    double *Yp_ref = (double*)b3f_xcalloc(C.yp_len, sizeof(double)), *Yw_ref = (double*)b3f_xcalloc(C.yw_len, sizeof(double));

    double lam = 0.5657237;
    for (int i = 0; i < C.NV; i++) f[i] = 1.0;

    int total_it = nit + 5; /* last 5 sweeps: freeze f (except VSTAR), fix lam = P/Q, converge (g,hV,hC,L4) with strict cushions */
    double target_lam = (double)H.P / (double)H.Q;

    for (int it = 0; it < total_it; it++) {
        int final_phase = (it >= nit);
        if (it == nit) {
            if (auto_pq) {
                /* Pick clean rational with safety cushion outside converged float lam */
                double raw_p = (mode == 1) ? ceil((lam + cushion) * (double)H.Q)
                                           : floor((lam - cushion) * (double)H.Q);
                if (!isfinite(raw_p) || raw_p < 1.0 || raw_p > (double)CERT_B3F_MAX_PQ) {
                    fprintf(stderr, "Auto P/Q out of valid range: %.17g\n", raw_p); return 1;
                }
                H.P = (int64_t)raw_p;
                if (b3f_l5_validate_header(&H,&E) != 0) {
                    fprintf(stderr, "Finalized header invalid\n"); return 1;
                }
                target_lam = (double)H.P / (double)H.Q;
            }
            printf("  Converged float lam=%.8f -> fixing exact target P/Q = %lld/%lld = %.10f\n",
                   lam, (long long)H.P, (long long)H.Q, target_lam);
            fflush(stdout);
            lam = target_lam;
        }

        double eps_g  = final_phase ? ((mode == 1) ? +1.5e-8 : -1.5e-8) : 0.0;
        double eps_h  = final_phase ? ((mode == 1) ? +3.0e-9 : -3.0e-9) : 0.0;
        double eps_L4 = final_phase ? ((mode == 1) ? +6.0e-10 : -6.0e-10) : 0.0;
        /* L4 consumes the previous L5 sweep. Refresh L5 against the newly
         * solved L4 below; fixed-lambda sweeps make the lag smaller than slack. */

        /* Freeze c2_f, c3_f, Yp_ref, Yw_ref before final_phase so all last sweeps use exact quantized constants */
        if (it >= 1 && it <= nit - 4) {
            for (int h = 1; h >= 0; h--) for (int kap = 0; kap <= KT; kap++) {
                int idx = b3f_tidx(&C, kap, h), D1 = (D0 > kap + h + 1) ? D0 : (kap + h + 1), D2 = 2 * D1;
                double d1 = f[b3f_vidx(&C, D1 - kap - h + 1, kap, h)] - f[b3f_vidx(&C, D1 - kap - h, kap, h)];
                double d2 = f[b3f_vidx(&C, D2 - kap - h + 1, kap, h)] - f[b3f_vidx(&C, D2 - kap - h, kap, h)];
                double y1 = (double)D1 * D1 * d1, y2 = (double)D2 * D2 * d2;
                c2_f[idx] = (D2 * y2 - D1 * y1) / (double)(D2 - D1);
                c3_f[idx] = (double)D1 * D2 * (y1 - y2) / (double)(D2 - D1);
            }
        }
        if (it <= nit - 1) {
            for (int h = 1; h >= 0; h--) for (int kap = 1; kap <= KT; kap++) {
                int idx = b3f_tidx(&C, kap, h), i1 = b3f_t1(&C, kap, h), i2 = b3f_t2(&C, kap, h), i3 = b3f_t3(&C, kap, h);
                for (int iw = 0; iw < C.nchs; iw++) {
                    int w = C.chs[iw];
                    Yp_ref[b3f_yp_off(&C, idx, w)] = eval_Y3_half(D0, g[C.INF_ID * N + w], g[i1 * N + w], g[i2 * N + w], g[i3 * N + w]);
                }
            }
            for (int kap = 1; kap <= KT; kap++) {
                int i1 = b3f_t1(&C, kap, 1), i2 = b3f_t2(&C, kap, 1), i3 = b3f_t3(&C, kap, 1);
                for (int ip = 0; ip < C.npars; ip++) {
                    int p = C.pars[ip];
                    Yw_ref[b3f_yw_off(&C, p, kap)] = eval_Y3_half(D0, g[p * N + C.INF_ID], g[p * N + i1], g[p * N + i2], g[p * N + i3]);
                }
            }
        }
        /* Snap f, c2_f, c3_f, Yp_ref, Yw_ref to exact 1/S multiples on entering final_phase */
        if (it == nit) {
            double invS = 1.0 / (double)H.S;
            for (int i = 0; i < C.NV; i++) f[i] = (double)quantize(f[i],H.S) * invS;
            for (int i = 0; i < NT; i++) { c2_f[i] = (double)quantize(c2_f[i],H.S) * invS; c3_f[i] = (double)quantize(c3_f[i],H.S) * invS; }
            for (size_t i = 0; i < C.yp_len; i++) Yp_ref[i] = (double)quantize(Yp_ref[i],H.S) * invS;
            for (size_t i = 0; i < C.yw_len; i++) Yw_ref[i] = (double)quantize(Yw_ref[i],H.S) * invS;
            memcpy(fn, f, C.NV * sizeof(double));
        }

        /* 1. Solve f (skipped during final_phase so f stays fixed with strict lam_target-lam_raw slack!) */
        if (!final_phase) {
        if (mode == 2) {
            fn[C.VSTAR0] = fn[C.VSTAR1] = 0.0;
        } else {
            double g0 = g[C.INF_ID * N + iw_nb];
            for (int hb = 0; hb <= 1; hb++) {
                double g1 = g[b3f_t1(&C, KT + 1, hb) * N + iw_nb];
                double g2 = g[b3f_t2(&C, KT + 1, hb) * N + iw_nb];
                double g3 = g[b3f_t3(&C, KT + 1, hb) * N + iw_nb];
                double tm = (double)D0 / (C.KF + 1 + hb), sm = 1.0 - tm;
                double q0 = g0, q1 = sm * g0 + tm * g1, q2 = sm * sm * g0 + 2 * tm * sm * g1 + tm * tm * g2;
                double q3 = sm * sm * sm * g0 + 3 * tm * sm * sm * g1 + 3 * tm * tm * sm * g2 + tm * tm * tm * g3;
                fn[hb ? C.VSTAR1 : C.VSTAR0] = fmax(0.0, f[v_nb] + fmax(fmax(q0, q1), fmax(q2, q3))) / lam;
            }
        }
        for (int hb = 1; hb >= 0; hb--) for (int k = C.KF; k >= 0; k--) for (int L = C.LF; L >= 0; L--) {
            int u = b3f_vidx(&C, L, k, hb), D = L + k + hb;
            if (D < 2) { fn[u] = 0.0; continue; }
            double bL = (double)L / D;
            if (L == C.LF) {
                int ub = b3f_vidx(&C, L - 1, k + 1, hb), ubt = b3f_vidx(&C, C.LF, k + 1, hb);
                int pid[4]; double pw[4];
                int np = b3f_embed(&C, L - 1, k + 1, hb, pid, pw);
                double gnb = 0.0; for (int i = 0; i < np; i++) gnb += pw[i] * g[pid[i] * N + iw_nb];
                double snb = f[v_nb] + gnb;
                double g0 = g[C.INF_ID * N + iw_nb], g1 = g[b3f_t1(&C, k + 1, hb) * N + iw_nb];
                double g2 = g[b3f_t2(&C, k + 1, hb) * N + iw_nb], g3 = g[b3f_t3(&C, k + 1, hb) * N + iw_nb];
                double tm = (double)D0 / (C.LF + 1 + k + hb), sm = 1.0 - tm;
                double q0 = g0, q1 = sm * g0 + tm * g1, q2 = sm * sm * g0 + 2 * tm * sm * g1 + tm * tm * g2;
                double q3 = sm * sm * sm * g0 + 3 * tm * sm * sm * g1 + 3 * tm * tm * sm * g2 + tm * tm * tm * g3;
                double snbt = f[v_nb] + ((mode == 1) ? fmax(fmax(q0, q1), fmax(q2, q3)) : fmin(fmin(q0, q1), fmin(q2, q3)));
                double blo = (double)(C.LF + 1) / (C.LF + 1 + k + hb), xfer_adj = 0.0;
                if (k <= KT) {
                    int idx = b3f_tidx(&C, k, hb);
                    double x0 = 1.0 / D, p2 = c2_f[idx] * x0 * x0 / 3.0, p3 = c2_f[idx] * x0 * x0 + c3_f[idx] * x0 * x0 * x0;
                    double pmin = fmin(0.0, fmin(p2, p3)), pmax = fmax(0.0, fmax(p2, p3));
                    double rr_lo = hb * trpar[idx] + k * trch[idx] - 0.5 * (hb + k), rr_hi = hb * trpar[idx] + k * trch[idx];
                    double z1 = rr_lo * pmin, z2 = rr_lo * pmax, z3 = rr_hi * pmin, z4 = rr_hi * pmax;
                    xfer_adj = (mode == 1) ? fmax(fmax(z1, z2), fmax(z3, z4)) : fmin(fmin(z1, z2), fmin(z3, z4));
                }
                double v1 = (bL * (fn[ub] + snb) + xfer_adj) / (lam + bL);
                double v2 = (blo * (fn[ubt] + snbt) + xfer_adj) / (lam + blo);
                double v3 = (fn[ubt] + snbt + xfer_adj) / (lam + 1.0);
                fn[u] = (mode == 1) ? fmax(v1, fmax(v2, v3)) : fmin(v1, fmin(v2, v3));
                continue;
            }
            double f_up = fn[b3f_vidx(&C, L + 1, k, hb)];
            int ub = (L >= 1) ? b3f_vidx(&C, L - 1, k + 1, hb) : u;
            double bv = 0.0;
            if (L >= 1) {
                int pid[4]; double pw[4];
                int np = b3f_embed(&C, L - 1, k + 1, hb, pid, pw);
                double gnb = 0.0; for (int i = 0; i < np; i++) gnb += pw[i] * g[pid[i] * N + iw_nb];
                bv = fn[ub] + f[v_nb] + gnb;
            }
            if (D < D0) {
                int gu = b3f_inner(&C, L, k, hb);
                double up = L + hb * rpar[gu] + k * rch[gu];
                fn[u] = (up * f_up + (L >= 1 ? bL * bv : 0.0)) / (lam + up + (L >= 1 ? bL : 0.0));
            } else if (k <= KT) {
                int idx = b3f_tidx(&C, k, hb);
                double dD = (double)D, psi = c2_f[idx] / (dD * dD) + c3_f[idx] / (dD * dD * dD);
                double base_r = hb * trpar[idx] + k * trch[idx];
                double nc = L * f_up + base_r * psi + (L >= 1 ? bL * bv : 0.0);
                double dc = lam + L + (L >= 1 ? bL : 0.0), rh = 0.5 * (hb + k);
                double v0 = nc / dc, v1 = (nc + rh * (f_up - psi)) / (dc + rh);
                fn[u] = (mode == 1) ? fmax(v0, v1) : fmin(v0, v1);
            } else {
                double r0 = (double)L, r1 = L + 0.5 * (hb + k);
                double v0 = (r0 * f_up + (L >= 1 ? bL * bv : 0.0)) / (lam + r0 + (L >= 1 ? bL : 0.0));
                double v1 = (r1 * f_up + (L >= 1 ? bL * bv : 0.0)) / (lam + r1 + (L >= 1 ? bL : 0.0));
                fn[u] = (mode == 1) ? fmax(v0, v1) : fmin(v0, v1);
            }
        }
        } else {
            /* Update only VSTAR0/1 for mode=1 in case g shifted slightly, keeping all finite f[u] frozen */
            if (mode == 1) {
                double g0 = g[C.INF_ID * N + iw_nb], invS = 1.0 / (double)H.S;
                for (int hb = 0; hb <= 1; hb++) {
                    double g1 = g[b3f_t1(&C, KT + 1, hb) * N + iw_nb];
                    double g2 = g[b3f_t2(&C, KT + 1, hb) * N + iw_nb];
                    double g3 = g[b3f_t3(&C, KT + 1, hb) * N + iw_nb];
                    double tm = (double)D0 / (C.KF + 1 + hb), sm = 1.0 - tm;
                    double q0 = g0, q1 = sm * g0 + tm * g1, q2 = sm * sm * g0 + 2 * tm * sm * g1 + tm * tm * g2;
                    double q3 = sm * sm * sm * g0 + 3 * tm * sm * sm * g1 + 3 * tm * tm * sm * g2 + tm * tm * tm * g3;
                    double val = (fmax(0.0, f[v_nb] + fmax(fmax(q0, q1), fmax(q2, q3))) + 1e-7) / lam;
                    fn[hb ? C.VSTAR1 : C.VSTAR0] = (double)quantize(val,H.S) * invS;
                }
            }
        }
        for (int u = 0; u < C.INNER_N; u++)
            dLf[u] = fn[b3f_vidx(&C, C.L[u] + 1, C.K[u], C.H[u])] - fn[b3f_vidx(&C, C.L[u], C.K[u], C.H[u])];

        /* 2. Solve g */
        memset(g, 0, C.g_len * sizeof(double));
        int mr = C.max_rank;
        for (int rsum = 2 * mr; rsum >= 4; rsum--) {
            int rp_lo = rsum - mr; if (rp_lo < 2) rp_lo = 2;
            int rp_hi = rsum - 2;  if (rp_hi > mr) rp_hi = mr;
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int rp = rp_lo; rp <= rp_hi; rp++) {
                int rw = rsum - rp;
                for (int ip = 0; ip < C.p_rcnt[rp]; ip++) {
                    int p = C.p_by_r[rp][ip], pk = C.kind[p], Lp = C.L[p], kp = C.K[p], hp = C.H[p];
                    int p_up = b3f_next_up(&C, p), p_bud = b3f_next_bud(&C, p);
                    int p_t1 = (pk >= 3) ? b3f_t1(&C, kp, hp) : p, p_t2 = (pk >= 4) ? b3f_t2(&C, kp, hp) : p;
                    int p_t1b = (pk >= 3 && kp <= KT) ? b3f_t1(&C, kp + 1, hp) : p_t1;
                    int p_t2b = (pk >= 4 && kp <= KT) ? b3f_t2(&C, kp + 1, hp) : p_t2;
                    int p_in_hu = (C.hu_of_g[p] >= 0), p_in_cu = (C.cu_of_g[p] >= 0);
                    int pidx_t = (pk >= 2 && kp <= KT) ? b3f_tidx(&C, kp, hp) : 0;
                    double inv_p = invD_of(&C, p);
                    int cu_p1b = (pk >= 2 && kp + 1 <= KT && kp + 1 >= 2) ? C.cu_of_g[b3f_t1(&C, kp + 1, hp)] : -1;
                    int cu_p2b = (pk >= 3 && kp + 1 <= KT && kp + 1 >= 2) ? C.cu_of_g[b3f_t2(&C, kp + 1, hp)] : -1;
                    int cu_p3b = (pk == 4 && kp + 1 <= KT && kp + 1 >= 2) ? C.cu_of_g[b3f_t3(&C, kp + 1, hp)] : ((pk == 0 && Lp >= 1) ? C.cu_of_g[p_bud] : -1);

                    for (int iw = 0; iw < C.c_rcnt[rw]; iw++) {
                        int w = C.c_by_r[rw][iw], wk = C.kind[w], Lw = C.L[w], kw = C.K[w];
                        if (pk == 1 && wk == 1) { g[p * N + w] = 0.0; continue; }
                        int w_up = b3f_next_up(&C, w), w_bud = b3f_next_bud(&C, w);
                        int w_t1 = (wk >= 3) ? b3f_t1(&C, kw, 1) : w, w_t2 = (wk >= 4) ? b3f_t2(&C, kw, 1) : w;
                        int w_t1b = (wk >= 3 && kw <= KT) ? b3f_t1(&C, kw + 1, 1) : w_t1;
                        int w_t2b = (wk >= 4 && kw <= KT) ? b3f_t2(&C, kw + 1, 1) : w_t2;
                        int w_in_hu = (C.hu_of_g[w] >= 0);
                        int widx_t = (wk >= 2 && kw <= KT) ? b3f_tidx(&C, kw, 1) : 0;
                        double inv_w = invD_of(&C, w);
                        int hu_w1b = (wk >= 2 && kw + 1 <= KT) ? C.hu_of_g[b3f_t1(&C, kw + 1, 1)] : -1;
                        int hu_w2b = (wk >= 3 && kw + 1 <= KT) ? C.hu_of_g[b3f_t2(&C, kw + 1, 1)] : -1;
                        int hu_w3b = (wk == 4 && kw + 1 <= KT) ? C.hu_of_g[b3f_t3(&C, kw + 1, 1)] : ((wk == 0 && Lw >= 1) ? C.hu_of_g[w_bud] : -1);

                        double gp_plus = g[p_up * N + w], gw_plus = g[p * N + w_up];
                        double gp_bud = g[p_bud * N + w] + ((pk == 0 && cu_p3b >= 0) ? hC[b3f_hc_off(&C, cu_p3b, w, iw_nb)] : 0.0);
                        double gw_bud = g[p * N + w_bud] + ((wk == 0 && hu_w3b >= 0) ? hV[b3f_hv_off(&C, p, hu_w3b, iw_nb)] : 0.0);

                        double pF0 = (pk >= 2) ? g[C.INF_ID * N + w] : 0, pF1 = (pk >= 3) ? g[p_t1 * N + w] : 0, pF2 = (pk >= 4) ? g[p_t2 * N + w] : 0;
                        double pF1b = (pk >= 2) ? (g[(pk == 2 ? p_bud : p_t1b) * N + w] + (cu_p1b >= 0 ? hC[b3f_hc_off(&C, cu_p1b, w, iw_nb)] : 0.0)) : 0;
                        double pF2b = (pk >= 3) ? (g[(pk == 3 ? p_bud : p_t2b) * N + w] + (cu_p2b >= 0 ? hC[b3f_hc_off(&C, cu_p2b, w, iw_nb)] : 0.0)) : 0;
                        double pF3b = (pk == 4) ? (g[p_bud * N + w] + (cu_p3b >= 0 ? hC[b3f_hc_off(&C, cu_p3b, w, iw_nb)] : 0.0)) : 0;
                        double pd0 = (pk >= 2 && wk == 0) ? (g[C.INF_ID * N + w_up] - pF0) : 0;
                        double pd1 = (pk >= 3 && wk == 0) ? (g[p_t1 * N + w_up] - pF1) : 0;
                        double pd2 = (pk >= 4 && wk == 0) ? (g[p_t2 * N + w_up] - pF2) : 0;

                        double wF0 = (wk >= 2) ? g[p * N + C.INF_ID] : 0, wF1 = (wk >= 3) ? g[p * N + w_t1] : 0, wF2 = (wk >= 4) ? g[p * N + w_t2] : 0;
                        double wF1b = (wk >= 2) ? (g[p * N + (wk == 2 ? w_bud : w_t1b)] + (hu_w1b >= 0 ? hV[b3f_hv_off(&C, p, hu_w1b, iw_nb)] : 0.0)) : 0;
                        double wF2b = (wk >= 3) ? (g[p * N + (wk == 3 ? w_bud : w_t2b)] + (hu_w2b >= 0 ? hV[b3f_hv_off(&C, p, hu_w2b, iw_nb)] : 0.0)) : 0;
                        double wF3b = (wk == 4) ? (g[p * N + w_bud] + (hu_w3b >= 0 ? hV[b3f_hv_off(&C, p, hu_w3b, iw_nb)] : 0.0)) : 0;
                        double wd0 = (wk >= 2 && pk == 0) ? (g[p_up * N + C.INF_ID] - wF0) : 0;
                        double wd1 = (wk >= 3 && pk == 0) ? (g[p_up * N + w_t1] - wF1) : 0;
                        double wd2 = (wk >= 4 && pk == 0) ? (g[p_up * N + w_t2] - wF2) : 0;

                        double num0 = eps_g, den0 = lam;
                        if (pk == 0) num0 += (inv_w - rch[p]) * dLf[p];
                        else if (pk == 3 && kp <= KT) num0 += (inv_w - trch[pidx_t]) * (c2_f[pidx_t] / (3.0 * D0 * D0));
                        else if (pk == 4 && kp <= KT) num0 += (inv_w - trch[pidx_t]) * (c2_f[pidx_t] / ((double)D0 * D0) + c3_f[pidx_t] / ((double)D0 * D0 * D0));
                        if (wk == 0) num0 += (inv_p - rpar[w]) * dLf[w];
                        else if (wk == 3 && kw <= KT) num0 += (inv_p - trpar[widx_t]) * (c2_f[widx_t] / (3.0 * D0 * D0));
                        else if (wk == 4 && kw <= KT) num0 += (inv_p - trpar[widx_t]) * (c2_f[widx_t] / ((double)D0 * D0) + c3_f[widx_t] / ((double)D0 * D0 * D0));

                        double rp_rem = 0, rw_rem = 0;
                        if (pk == 0) {
                            double ur = Lp + (wk == 0 ? inv_w : 0) + (p_in_hu ? rpar[p] : 0) + (p_in_cu ? (kp - 1) * rch[p] : 0);
                            num0 += ur * gp_plus; den0 += ur;
                            if (Lp >= 1) { double b = (double)Lp / C.D[p]; num0 += b * gp_bud; den0 += b; }
                            if (hp && !p_in_hu) rp_rem += 0.5;
                            if (kp >= 2 && !p_in_cu) rp_rem += 0.5 * (kp - 1);
                        }
                        if (wk == 0) {
                            double ur = Lw + (pk == 0 ? inv_p : 0) + (w_in_hu ? kw * rch[w] : 0);
                            num0 += ur * gw_plus; den0 += ur;
                            if (Lw >= 1) { double b = (double)Lw / C.D[w]; num0 += b * gw_bud; den0 += b; }
                            if (kw >= 1 && !w_in_hu) rw_rem += 0.5 * kw;
                        }

                        double pr_in_lo = 0, pr_in_hi = 0, wr_in_lo = 0, wr_in_hi = 0;
                        if (pk >= 2) get_rin(&C, w, &pr_in_lo, &pr_in_hi);
                        if (wk >= 2) get_rin(&C, p, &wr_in_lo, &wr_in_hi);
                        double p_Rbar = (pk >= 2 && kp <= KT) ? (hp * trpar[pidx_t] + (kp - 1) * trch[pidx_t]) : 0.0;
                        double w_Rbar = (wk >= 2 && kw <= KT) ? (kw * trch[widx_t]) : 0.0;
                        double p_Yref = (pk == 4 && kp <= KT) ? Yp_ref[b3f_yp_off(&C, pidx_t, w)] : 0.0;
                        double w_Yref = (wk == 4 && kw <= KT && kw >= 1) ? Yw_ref[b3f_yw_off(&C, p, kw)] : 0.0;

                        b3_tail_const(mode, D0, KT, pk, kp, hp, pr_in_lo + p_Rbar, pr_in_hi + p_Rbar, pF0, pF1, pF1b, pF2b, wk == 0, pd0, pd1, &num0, &den0);
                        b3_tail_const(mode, D0, KT, wk, kw, 1,  wr_in_lo + w_Rbar, wr_in_hi + w_Rbar, wF0, wF1, wF1b, wF2b, pk == 0, wd0, wd1, &num0, &den0);

                        double v = num0 / den0;
                        if (rp_rem > 0 || rw_rem > 0 || (pk == 2 && kp > KT) || (pk == 3 && kp > KT) || pk == 4 || (wk == 2 && kw > KT) || (wk == 3 && kw > KT) || wk == 4) {
                            int mp = -1;
                            for (int ns = 0; ns < 16; ns++) {
                                double num = num0, den = den0; int m = 0;
                                if (rp_rem > 0 && ((mode == 1 && gp_plus > v) || (mode == 2 && gp_plus < v))) { num += rp_rem * gp_plus; den += rp_rem; m |= 1; }
                                if (rw_rem > 0 && ((mode == 1 && gw_plus > v) || (mode == 2 && gw_plus < v))) { num += rw_rem * gw_plus; den += rw_rem; m |= 2; }
                                b3_tail_var(mode, D0, KT, pk, kp, hp, pr_in_lo, pr_in_hi, 0.5 * (kp + hp - 1), p_Rbar, p_Yref, pF0, pF1, pF2, pF1b, pF2b, pF3b, wk == 0, pd0, pd1, pd2, gw_plus, v, &num, &den, &m, 2);
                                b3_tail_var(mode, D0, KT, wk, kw, 1,  wr_in_lo, wr_in_hi, 0.5 * kw,            w_Rbar, w_Yref, wF0, wF1, wF2, wF1b, wF2b, wF3b, pk == 0, wd0, wd1, wd2, gp_plus, v, &num, &den, &m, 10);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                        }
                         if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                        g[p * N + w] = v;
                    }
                }
            }
        }

        /* 3. Solve FULL hV */
        memset(hV, 0, C.hv_len * sizeof(double));
        int ru_max = D0 + 2;
        for (int rsum = 2 * mr + ru_max; rsum >= 6; rsum--) {
            int tasks[16384][3], nt = 0;
            for (int ru = 2; ru <= ru_max; ru++) if (C.hu_rcnt[ru] > 0)
                for (int rp = 2; rp <= mr; rp++) {
                    int rw = rsum - ru - rp;
                    if (rw >= 2 && rw <= mr && C.p_rcnt[rp] > 0 && C.c_rcnt[rw] > 0) {
                        tasks[nt][0] = rp; tasks[nt][1] = ru; tasks[nt][2] = rw; nt++;
                    }
                }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ti = 0; ti < nt; ti++) {
                int rp = tasks[ti][0], ru = tasks[ti][1], rw = tasks[ti][2];
                for (int iua = 0; iua < C.hu_rcnt[ru]; iua++) {
                    int ua = C.hu_by_r[ru][iua], gu = C.hu_g[ua], uk = C.kind[gu], Lu = C.L[gu], ku = C.K[gu];
                    int gup = (uk == 0) ? b3f_next_up(&C, gu) : C.INF_ID;
                    int ua_up = (uk == 0) ? C.hu_of_g[gup] : -1;
                    int gu_bud = b3f_next_bud(&C, gu);
                    int ua_bud = (uk == 0 && Lu == 0) ? -1 : ((ku < KT) ? C.hu_of_g[gu_bud] : -1);
                    int u_t1 = (uk >= 3) ? C.hu_of_g[b3f_t1(&C, ku, 1)] : ua;
                    int u_t2 = (uk >= 4) ? C.hu_of_g[b3f_t2(&C, ku, 1)] : ua;
                    int u_t1b = (uk >= 3 && ku < KT) ? C.hu_of_g[b3f_t1(&C, ku + 1, 1)] : -1;
                    int u_t2b = (uk >= 4 && ku < KT) ? C.hu_of_g[b3f_t2(&C, ku + 1, 1)] : -1;
                    int uidx_t = (uk >= 2) ? b3f_tidx(&C, ku, 1) : 0;
                    double iu = (uk == 0) ? (1.0 / ru) : 0.0, bu = (uk == 0) ? ((double)Lu / ru) : 0.0;
                    double urpar = (uk == 0) ? rpar[gu] : trpar[uidx_t];
                    double urch  = (uk == 0) ? rch[gu]  : trch[uidx_t];

                    int u_p4  = (uk == 0 && ku <= C.KP4 && ru < C.DP4) ? C.p4_lut[ku][Lu] : -1;
                    int u_yt  = (uk == 0 && ku >= 2 && ku <= C.KYT && ru < C.DYT) ? C.yt_lut[ku][Lu] : -1;
                    int ub_yt = (uk == 0 && Lu >= 1 && ku + 1 >= 2 && ku + 1 <= C.KYT && ru < C.DYT) ? C.yt_lut[ku + 1][Lu - 1] : -1;
                    int u_ybw = (uk == 0 && ku <= C.KYB && ru < C.DYB) ? C.ybw_lut[ku][Lu] : -1;

                    for (int ip = 0; ip < C.p_rcnt[rp]; ip++) {
                        int p = C.p_by_r[rp][ip], pk = C.kind[p], Lp = C.L[p], kp = C.K[p], hp = C.H[p];
                        int p_up = b3f_next_up(&C, p), p_bud = b3f_next_bud(&C, p);
                        int p_t1 = (pk >= 3) ? b3f_t1(&C, kp, hp) : p, p_t2 = (pk >= 4) ? b3f_t2(&C, kp, hp) : p;
                        int p_t1b = (pk >= 3 && kp <= KT) ? b3f_t1(&C, kp + 1, hp) : p_t1;
                        int p_t2b = (pk >= 4 && kp <= KT) ? b3f_t2(&C, kp + 1, hp) : p_t2;
                        double inv_p = invD_of(&C, p), dgpu = 0.0;
                        if (uk == 0) dgpu = g[p * N + gup] - g[p * N + gu];
                        else if (uk == 3) dgpu = (g[p * N + C.INF_ID] - g[p * N + b3f_t1(&C, ku, 1)]) / D0;
                        else if (uk == 4) dgpu = Yw_ref[b3f_yw_off(&C, p, ku)];

                        int p_p4   = (pk == 0 && hp == 1 && kp >= 1 && kp <= C.KP4 && C.D[p] < C.DP4) ? C.p4_lut[kp][Lp] : -1;
                        int p_ybu  = (pk == 0 && kp >= 2 && kp <= C.KYB && C.D[p] < C.DYB) ? C.ybu_lut[hp][kp][Lp] : -1;
                        int pb_ybu = (pk == 0 && Lp >= 1 && kp + 1 >= 2 && kp + 1 <= C.KYB && C.D[p] < C.DYB) ? C.ybu_lut[hp][kp + 1][Lp - 1] : -1;
                        int x_pp = (p_p4 >= 0 && u_p4 >= 0), x_pc = (p_ybu >= 0 && u_ybw >= 0);

                        for (int iw = 0; iw < C.c_rcnt[rw]; iw++) {
                            int w = C.c_by_r[rw][iw], wk = C.kind[w], Lw = C.L[w], kw = C.K[w];
                            int w_up = b3f_next_up(&C, w), w_bud = b3f_next_bud(&C, w);
                            int w_t1 = (wk >= 3) ? b3f_t1(&C, kw, 1) : w, w_t2 = (wk >= 4) ? b3f_t2(&C, kw, 1) : w;
                            int w_t1b = (wk >= 3 && kw <= KT) ? b3f_t1(&C, kw + 1, 1) : w_t1;
                            int w_t2b = (wk >= 4 && kw <= KT) ? b3f_t2(&C, kw + 1, 1) : w_t2;
                            double inv_w = invD_of(&C, w), dguw = 0.0;
                            if (uk == 0) dguw = g[gup * N + w] - g[gu * N + w];
                            else if (uk == 3) dguw = (g[C.INF_ID * N + w] - g[b3f_t1(&C, ku, 1) * N + w]) / D0;
                            else if (uk == 4) dguw = Yp_ref[b3f_yp_off(&C, uidx_t, w)];

                            int w_p4  = (wk == 0 && kw >= 1 && kw <= C.KP4 && C.D[w] < C.DP4) ? C.p4_lut[kw][Lw] : -1;
                            int wb_p4 = (wk == 0 && Lw >= 1 && kw + 1 <= C.KP4 && C.D[w] < C.DP4) ? C.p4_lut[kw + 1][Lw - 1] : -1;
                            int x_wc = (u_p4 >= 0 && w_p4 >= 0), x_uz = (u_yt >= 0);

                            double hp_up = hV[b3f_hv_off(&C, p_up, ua, w)];
                            double hu_up = (ua_up >= 0) ? hV[b3f_hv_off(&C, p, ua_up, w)] : 0.0;
                            double hw_up = hV[b3f_hv_off(&C, p, ua, w_up)];
                            double hp_bud = hV[b3f_hv_off(&C, p_bud, ua, w)];
                            double hu_bud = (ua_bud >= 0) ? hV[b3f_hv_off(&C, p, ua_bud, w)] : 0.0;
                            double hw_bud = hV[b3f_hv_off(&C, p, ua, w_bud)];
                            double nb_p = (pb_ybu >= 0 && u_ybw >= 0) ? qYb[b3f_yb_off(&C, pb_ybu, u_ybw, iw_nb, w)] : 0.0;
                            double nb_u = (ub_yt >= 0) ? qYt[b3f_yt_off(&C, p, ub_yt, w, iw_nb)] : 0.0;
                            double nb_w = (u_p4 >= 0 && wb_p4 >= 0) ? qP4[b3f_p4_off(&C, p, u_p4, wb_p4, iw_nb)] : 0.0;

                            double num0 = (inv_p - urpar) * dguw + (inv_w - urch) * dgpu + eps_h, den0 = lam;
                            double rp_hi = 0, ru_hi = 0, rw_hi = 0;
                            if (pk == 0) {
                                double r = Lp + (uk == 0 ? iu : 0) + (x_pp ? rpar[p] : 0) + (x_pc ? (kp - 1) * rch[p] : 0);
                                num0 += r * hp_up; den0 += r;
                                if (Lp >= 1) { double b = (double)Lp / C.D[p]; num0 += b * (hp_bud + nb_p); den0 += b; }
                                if (hp && !x_pp) rp_hi += 0.5;
                                if (kp >= 2 && !x_pc) rp_hi += 0.5 * (kp - 1);
                            }
                            if (wk == 0) {
                                double r = Lw + (uk == 0 ? iu : 0) + (x_wc ? kw * rch[w] : 0);
                                num0 += r * hw_up; den0 += r;
                                if (Lw >= 1) { double b = (double)Lw / C.D[w]; num0 += b * (hw_bud + nb_w); den0 += b; }
                                if (kw >= 1 && !x_wc) rw_hi += 0.5 * kw;
                            }
                            if (uk == 0) {
                                double r = Lu + (pk == 0 ? inv_p : 0) + (wk == 0 ? inv_w : 0) + (x_uz ? (ku - 1) * rch[gu] : 0);
                                num0 += r * hu_up; den0 += r;
                                if (bu > 0) { num0 += bu * (hu_bud + nb_u); den0 += bu; }
                                if (ku >= 2 && !x_uz) ru_hi += 0.5 * (ku - 1);
                            }

                            double pr_in_lo = 0, pr_in_hi = 0; if (pk >= 2) get_rin(&C, gu, &pr_in_lo, &pr_in_hi);
                            double pF0 = (pk >= 2) ? hV[b3f_hv_off(&C, C.INF_ID, ua, w)] : 0;
                            double pF1 = (pk >= 3) ? hV[b3f_hv_off(&C, p_t1, ua, w)] : 0;
                            double pF2 = (pk >= 4) ? hV[b3f_hv_off(&C, p_t2, ua, w)] : 0;
                            double pF1b = (pk >= 3) ? hV[b3f_hv_off(&C, p_t1b, ua, w)] : hp_bud;
                            double pF2b = (pk >= 4) ? hV[b3f_hv_off(&C, p_t2b, ua, w)] : hp_bud;
                            double pd0 = (pk >= 2 && uk == 0) ? (((ua_up >= 0) ? hV[b3f_hv_off(&C, C.INF_ID, ua_up, w)] : 0) - pF0) : 0;
                            double pd1 = (pk >= 3 && uk == 0) ? (((ua_up >= 0) ? hV[b3f_hv_off(&C, p_t1, ua_up, w)] : 0) - pF1) : 0;
                            double pd2 = (pk >= 4 && uk == 0) ? (((ua_up >= 0) ? hV[b3f_hv_off(&C, p_t2, ua_up, w)] : 0) - pF2) : 0;

                            double wr_in_lo = 0, wr_in_hi = 0; if (wk >= 2) get_rin(&C, gu, &wr_in_lo, &wr_in_hi);
                            double wF0 = (wk >= 2) ? hV[b3f_hv_off(&C, p, ua, C.INF_ID)] : 0;
                            double wF1 = (wk >= 3) ? hV[b3f_hv_off(&C, p, ua, w_t1)] : 0;
                            double wF2 = (wk >= 4) ? hV[b3f_hv_off(&C, p, ua, w_t2)] : 0;
                            double wF1b = (wk >= 3) ? hV[b3f_hv_off(&C, p, ua, w_t1b)] : hw_bud;
                            double wF2b = (wk >= 4) ? hV[b3f_hv_off(&C, p, ua, w_t2b)] : hw_bud;
                            double wd0 = (wk >= 2 && uk == 0) ? (((ua_up >= 0) ? hV[b3f_hv_off(&C, p, ua_up, C.INF_ID)] : 0) - wF0) : 0;
                            double wd1 = (wk >= 3 && uk == 0) ? (((ua_up >= 0) ? hV[b3f_hv_off(&C, p, ua_up, w_t1)] : 0) - wF1) : 0;
                            double wd2 = (wk >= 4 && uk == 0) ? (((ua_up >= 0) ? hV[b3f_hv_off(&C, p, ua_up, w_t2)] : 0) - wF2) : 0;

                            double ur_in_lo = 0, ur_in_hi = 0;
                            if (uk >= 2) {
                                double lp, hp_, lw, hw_; get_rin(&C, p, &lp, &hp_); get_rin(&C, w, &lw, &hw_);
                                ur_in_lo = lp + lw; ur_in_hi = hp_ + hw_;
                            }
                            double uF0 = 0.0, uF1 = (uk >= 3) ? hV[b3f_hv_off(&C, p, u_t1, w)] : 0.0, uF2 = (uk >= 4) ? hV[b3f_hv_off(&C, p, u_t2, w)] : 0.0;
                            double uF1b = (uk == 2) ? hu_bud : ((uk >= 3 && u_t1b >= 0) ? hV[b3f_hv_off(&C, p, u_t1b, w)] : 0.0);
                            double uF2b = (uk == 3) ? hu_bud : ((uk >= 4 && u_t2b >= 0) ? hV[b3f_hv_off(&C, p, u_t2b, w)] : 0.0);
                            int u_nin = (uk >= 2) ? ((pk == 0 ? 1 : 0) + (wk == 0 ? 1 : 0)) : 0;
                            double ud0 = 0.0;
                            double ud1 = (uk >= 3 && u_nin > 0) ? ((pk == 0 ? (hV[b3f_hv_off(&C, p_up, u_t1, w)] - uF1) : 0.0) + (wk == 0 ? (hV[b3f_hv_off(&C, p, u_t1, w_up)] - uF1) : 0.0)) : 0.0;
                            double ud2 = (uk >= 4 && u_nin > 0) ? ((pk == 0 ? (hV[b3f_hv_off(&C, p_up, u_t2, w)] - uF2) : 0.0) + (wk == 0 ? (hV[b3f_hv_off(&C, p, u_t2, w_up)] - uF2) : 0.0)) : 0.0;
                            double uF3_up_sum = (uk == 4 && u_nin > 0) ? ((pk == 0 ? hp_up : 0.0) + (wk == 0 ? hw_up : 0.0)) : 0.0;

                            b3_tail_const(mode, D0, KT, pk, kp, hp, pr_in_lo, pr_in_hi + 0.5 * (kp + hp - 1), pF0, pF1, (pk == 2) ? hp_bud : pF1b, hp_bud, uk == 0, pd0, pd1, &num0, &den0);
                            b3_tail_const(mode, D0, KT, wk, kw, 1,  wr_in_lo, wr_in_hi + 0.5 * kw,            wF0, wF1, (wk == 2) ? hw_bud : wF1b, hw_bud, uk == 0, wd0, wd1, &num0, &den0);
                            b3_tail_const(mode, D0, KT, uk, ku, 1,  ur_in_lo, ur_in_hi + 0.5 * (ku - 1),      uF0, uF1, uF1b, uF2b, u_nin > 0, ud0, ud1, &num0, &den0);

                            double v = num0 / den0; int mp = -1;
                            for (int ns = 0; ns < 12; ns++) {
                                double num = num0, den = den0; int m = 0;
                                if (rp_hi > 0 && ((mode == 1 && hp_up > v) || (mode == 2 && hp_up < v))) { num += rp_hi * hp_up; den += rp_hi; m |= 1; }
                                if (ru_hi > 0 && ((mode == 1 && hu_up > v) || (mode == 2 && hu_up < v))) { num += ru_hi * hu_up; den += ru_hi; m |= 2; }
                                if (rw_hi > 0 && ((mode == 1 && hw_up > v) || (mode == 2 && hw_up < v))) { num += rw_hi * hw_up; den += rw_hi; m |= 4; }
                                b3_tail_var(mode, D0, KT, pk, kp, hp, pr_in_lo, pr_in_hi, 0.5 * (kp + hp - 1), 0, 0, pF0, pF1, pF2, pF1b, pF2b, hp_bud, uk == 0, pd0, pd1, pd2, hu_up, v, &num, &den, &m, 4);
                                b3_tail_var(mode, D0, KT, wk, kw, 1,  wr_in_lo, wr_in_hi, 0.5 * kw,            0, 0, wF0, wF1, wF2, wF1b, wF2b, hw_bud, uk == 0, wd0, wd1, wd2, hu_up, v, &num, &den, &m, 10);
                                b3_tail_var(mode, D0, KT, uk, ku, 1,  ur_in_lo, ur_in_hi, 0.5 * (ku - 1),      0, 0, uF0, uF1, uF2, uF1b, uF2b, hu_bud, u_nin, ud0, ud1, ud2, uF3_up_sum, v, &num, &den, &m, 16);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                             if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                            hV[b3f_hv_off(&C, p, ua, w)] = v;
                        }
                    }
                }
            }
        }

        /* 4. Solve FULL hC */
        memset(hC, 0, C.hc_len * sizeof(double));
        for (int rsum = ru_max + 2 * mr; rsum >= 6; rsum--) {
            int tasks[16384][3], nt = 0;
            for (int ru = 2; ru <= ru_max; ru++) if (C.cu_rcnt[ru] > 0)
                for (int r1 = 2; r1 <= mr; r1++) {
                    int r2 = rsum - ru - r1;
                    if (r2 >= r1 && r2 <= mr && C.c_rcnt[r1] > 0 && C.c_rcnt[r2] > 0) {
                        tasks[nt][0] = ru; tasks[nt][1] = r1; tasks[nt][2] = r2; nt++;
                    }
                }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ti = 0; ti < nt; ti++) {
                int ru = tasks[ti][0], r1 = tasks[ti][1], r2 = tasks[ti][2];
                for (int ica = 0; ica < C.cu_rcnt[ru]; ica++) {
                    int ca = C.cu_by_r[ru][ica], gu = C.cu_g[ca], uk = C.kind[gu], Lu = C.L[gu], ku = C.K[gu], hu = C.H[gu];
                    int gup = (uk == 0) ? b3f_next_up(&C, gu) : C.INF_ID;
                    int ca_up = (uk == 0) ? C.cu_of_g[gup] : -1;
                    int gu_bud = b3f_next_bud(&C, gu);
                    int ca_bud = (uk == 0 && Lu == 0) ? -1 : ((ku < KT) ? C.cu_of_g[gu_bud] : -1);
                    int u_t1 = (uk >= 3) ? C.cu_of_g[b3f_t1(&C, ku, hu)] : ca;
                    int u_t2 = (uk >= 4) ? C.cu_of_g[b3f_t2(&C, ku, hu)] : ca;
                    int u_t1b = (uk >= 3 && ku < KT) ? C.cu_of_g[b3f_t1(&C, ku + 1, hu)] : -1;
                    int u_t2b = (uk >= 4 && ku < KT) ? C.cu_of_g[b3f_t2(&C, ku + 1, hu)] : -1;
                    int uidx_t = (uk >= 2) ? b3f_tidx(&C, ku, hu) : 0;
                    double iu = (uk == 0) ? (1.0 / ru) : 0.0, bu = (uk == 0) ? ((double)Lu / ru) : 0.0;
                    double urch = (uk == 0) ? rch[gu] : trch[uidx_t];

                    int u_yt  = (uk == 0 && hu == 1 && ku <= C.KYT && ru < C.DYT) ? C.yt_lut[ku][Lu] : -1;
                    int u_ybu = (uk == 0 && ku <= C.KYB && ru < C.DYB) ? C.ybu_lut[hu][ku][Lu] : -1;
                    int u_s3  = (uk == 0 && ku >= 3 && ku <= C.KS3 && ru < C.DS3) ? C.s3_lut[hu][ku][Lu] : -1;
                    int ub_s3 = (uk == 0 && Lu >= 1 && ku + 1 >= 3 && ku + 1 <= C.KS3 && ru < C.DS3) ? C.s3_lut[hu][ku + 1][Lu - 1] : -1;
                    int x_up = (u_yt >= 0), x_uz = (u_s3 >= 0);

                    for (int i1 = 0; i1 < C.c_rcnt[r1]; i1++) {
                        int w1 = C.c_by_r[r1][i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                        int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1);
                        int w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                        int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1;
                        int w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                        double inv1 = invD_of(&C, w1), dg1 = 0.0;
                        if (uk == 0) dg1 = g[gup * N + w1] - g[gu * N + w1];
                        else if (uk == 3) dg1 = (g[C.INF_ID * N + w1] - g[b3f_t1(&C, ku, hu) * N + w1]) / D0;
                        else if (uk == 4) dg1 = Yp_ref[b3f_yp_off(&C, uidx_t, w1)];

                        int w1_ybw  = (k1k == 0 && k1 >= 1 && k1 <= C.KYB && C.D[w1] < C.DYB) ? C.ybw_lut[k1][L1] : -1;
                        int w1b_ybw = (k1k == 0 && L1 >= 1 && k1 + 1 <= C.KYB && C.D[w1] < C.DYB) ? C.ybw_lut[k1 + 1][L1 - 1] : -1;
                        int x_w1 = (u_ybu >= 0 && w1_ybw >= 0);

                        for (int i2 = (r1 == r2 ? i1 : 0); i2 < C.c_rcnt[r2]; i2++) {
                            int w2 = C.c_by_r[r2][i2], k2k = C.kind[w2], L2 = C.L[w2], k2 = C.K[w2];
                            int w2_up = b3f_next_up(&C, w2), w2_bud = b3f_next_bud(&C, w2);
                            int w2_t1 = (k2k >= 3) ? b3f_t1(&C, k2, 1) : w2, w2_t2 = (k2k >= 4) ? b3f_t2(&C, k2, 1) : w2;
                            int w2_t1b = (k2k >= 3 && k2 <= KT) ? b3f_t1(&C, k2 + 1, 1) : w2_t1;
                            int w2_t2b = (k2k >= 4 && k2 <= KT) ? b3f_t2(&C, k2 + 1, 1) : w2_t2;
                            double inv2 = invD_of(&C, w2), dg2 = 0.0;
                            if (uk == 0) dg2 = g[gup * N + w2] - g[gu * N + w2];
                            else if (uk == 3) dg2 = (g[C.INF_ID * N + w2] - g[b3f_t1(&C, ku, hu) * N + w2]) / D0;
                            else if (uk == 4) dg2 = Yp_ref[b3f_yp_off(&C, uidx_t, w2)];

                            int w2_ybw  = (k2k == 0 && k2 >= 1 && k2 <= C.KYB && C.D[w2] < C.DYB) ? C.ybw_lut[k2][L2] : -1;
                            int w2b_ybw = (k2k == 0 && L2 >= 1 && k2 + 1 <= C.KYB && C.D[w2] < C.DYB) ? C.ybw_lut[k2 + 1][L2 - 1] : -1;
                            int x_w2 = (u_ybu >= 0 && w2_ybw >= 0);

                            double cu_up = (ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, w1, w2)] : 0.0;
                            double cu_bud = (ca_bud >= 0) ? hC[b3f_hc_off(&C, ca_bud, w1, w2)] : 0.0;
                            double c1_up = hC[b3f_hc_off(&C, ca, w1_up, w2)], c2_up = hC[b3f_hc_off(&C, ca, w1, w2_up)];
                            double c1_bud = hC[b3f_hc_off(&C, ca, w1_bud, w2)], c2_bud = hC[b3f_hc_off(&C, ca, w1, w2_bud)];
                            double nb_u = (ub_s3 >= 0) ? qS3[b3f_s3_off(&C, ub_s3, w1, w2, iw_nb)] : 0.0;
                            double nb_1 = (u_ybu >= 0 && w1b_ybw >= 0) ? qYb[b3f_yb_off(&C, u_ybu, w1b_ybw, w2, iw_nb)] : 0.0;
                            double nb_2 = (u_ybu >= 0 && w2b_ybw >= 0) ? qYb[b3f_yb_off(&C, u_ybu, w2b_ybw, w1, iw_nb)] : 0.0;

                            double num0 = (inv1 - urch) * dg2 + (inv2 - urch) * dg1 + eps_h, den0 = lam;
                            double ru_hi = 0, e1_hi = 0, e2_hi = 0;
                            if (k1k == 0) {
                                double r = L1 + (uk == 0 ? iu : 0) + (x_w1 ? k1 * rch[w1] : 0);
                                num0 += r * c1_up; den0 += r;
                                if (L1 >= 1) { double b = (double)L1 / C.D[w1]; num0 += b * (c1_bud + nb_1); den0 += b; }
                                if (k1 >= 1 && !x_w1) e1_hi += 0.5 * k1;
                            }
                            if (k2k == 0) {
                                double r = L2 + (uk == 0 ? iu : 0) + (x_w2 ? k2 * rch[w2] : 0);
                                num0 += r * c2_up; den0 += r;
                                if (L2 >= 1) { double b = (double)L2 / C.D[w2]; num0 += b * (c2_bud + nb_2); den0 += b; }
                                if (k2 >= 1 && !x_w2) e2_hi += 0.5 * k2;
                            }
                            if (uk == 0) {
                                double r = Lu + (k1k == 0 ? inv1 : 0) + (k2k == 0 ? inv2 : 0) + (x_up ? rpar[gu] : 0) + (x_uz ? (ku - 2) * rch[gu] : 0);
                                num0 += r * cu_up; den0 += r;
                                if (bu > 0) { num0 += bu * (cu_bud + nb_u); den0 += bu; }
                                if (hu && !x_up) ru_hi += 0.5;
                                if (ku >= 3 && !x_uz) ru_hi += 0.5 * (ku - 2);
                            }

                            double r1_lo = 0, r1_hi = 0, r2_lo = 0, r2_hi = 0;
                            if (k1k >= 2) get_rin(&C, gu, &r1_lo, &r1_hi);
                            if (k2k >= 2) get_rin(&C, gu, &r2_lo, &r2_hi);
                            double f1_0 = (k1k >= 2) ? hC[b3f_hc_off(&C, ca, C.INF_ID, w2)] : 0;
                            double f1_1 = (k1k >= 3) ? hC[b3f_hc_off(&C, ca, w1_t1, w2)] : 0;
                            double f1_2 = (k1k >= 4) ? hC[b3f_hc_off(&C, ca, w1_t2, w2)] : 0;
                            double f1_1b = (k1k >= 3) ? hC[b3f_hc_off(&C, ca, w1_t1b, w2)] : c1_bud;
                            double f1_2b = (k1k >= 4) ? hC[b3f_hc_off(&C, ca, w1_t2b, w2)] : c1_bud;
                            double d1_0 = (k1k >= 2 && uk == 0) ? (((ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, C.INF_ID, w2)] : 0) - f1_0) : 0;
                            double d1_1 = (k1k >= 3 && uk == 0) ? (((ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, w1_t1, w2)] : 0) - f1_1) : 0;
                            double d1_2 = (k1k >= 4 && uk == 0) ? (((ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, w1_t2, w2)] : 0) - f1_2) : 0;

                            double f2_0 = (k2k >= 2) ? hC[b3f_hc_off(&C, ca, w1, C.INF_ID)] : 0;
                            double f2_1 = (k2k >= 3) ? hC[b3f_hc_off(&C, ca, w1, w2_t1)] : 0;
                            double f2_2 = (k2k >= 4) ? hC[b3f_hc_off(&C, ca, w1, w2_t2)] : 0;
                            double f2_1b = (k2k >= 3) ? hC[b3f_hc_off(&C, ca, w1, w2_t1b)] : c2_bud;
                            double f2_2b = (k2k >= 4) ? hC[b3f_hc_off(&C, ca, w1, w2_t2b)] : c2_bud;
                            double d2_0 = (k2k >= 2 && uk == 0) ? (((ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, w1, C.INF_ID)] : 0) - f2_0) : 0;
                            double d2_1 = (k2k >= 3 && uk == 0) ? (((ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, w1, w2_t1)] : 0) - f2_1) : 0;
                            double d2_2 = (k2k >= 4 && uk == 0) ? (((ca_up >= 0) ? hC[b3f_hc_off(&C, ca_up, w1, w2_t2)] : 0) - f2_2) : 0;

                            double ur_in_lo = 0, ur_in_hi = 0;
                            if (uk >= 2) {
                                double l1, h1_, l2, h2_; get_rin(&C, w1, &l1, &h1_); get_rin(&C, w2, &l2, &h2_);
                                ur_in_lo = l1 + l2; ur_in_hi = h1_ + h2_;
                            }
                            double uF0 = 0.0, uF1 = (uk >= 3) ? hC[b3f_hc_off(&C, u_t1, w1, w2)] : 0.0, uF2 = (uk >= 4) ? hC[b3f_hc_off(&C, u_t2, w1, w2)] : 0.0;
                            double uF1b = (uk == 2) ? cu_bud : ((uk >= 3 && u_t1b >= 0) ? hC[b3f_hc_off(&C, u_t1b, w1, w2)] : 0.0);
                            double uF2b = (uk == 3) ? cu_bud : ((uk >= 4 && u_t2b >= 0) ? hC[b3f_hc_off(&C, u_t2b, w1, w2)] : 0.0);
                            int u_nin = (uk >= 2) ? ((k1k == 0 ? 1 : 0) + (k2k == 0 ? 1 : 0)) : 0;
                            double ud0 = 0.0;
                            double ud1 = (uk >= 3 && u_nin > 0) ? ((k1k == 0 ? (hC[b3f_hc_off(&C, u_t1, w1_up, w2)] - uF1) : 0.0) + (k2k == 0 ? (hC[b3f_hc_off(&C, u_t1, w1, w2_up)] - uF1) : 0.0)) : 0.0;
                            double ud2 = (uk >= 4 && u_nin > 0) ? ((k1k == 0 ? (hC[b3f_hc_off(&C, u_t2, w1_up, w2)] - uF2) : 0.0) + (k2k == 0 ? (hC[b3f_hc_off(&C, u_t2, w1, w2_up)] - uF2) : 0.0)) : 0.0;
                            double uF3_up_sum = (uk == 4 && u_nin > 0) ? ((k1k == 0 ? c1_up : 0.0) + (k2k == 0 ? c2_up : 0.0)) : 0.0;

                            b3_tail_const(mode, D0, KT, k1k, k1, 1, r1_lo, r1_hi + 0.5 * k1, f1_0, f1_1, (k1k == 2) ? c1_bud : f1_1b, c1_bud, uk == 0, d1_0, d1_1, &num0, &den0);
                            b3_tail_const(mode, D0, KT, k2k, k2, 1, r2_lo, r2_hi + 0.5 * k2, f2_0, f2_1, (k2k == 2) ? c2_bud : f2_1b, c2_bud, uk == 0, d2_0, d2_1, &num0, &den0);
                            b3_tail_const(mode, D0, KT, uk, ku, hu, ur_in_lo, ur_in_hi + 0.5 * (hu + ku - 2), uF0, uF1, uF1b, uF2b, u_nin > 0, ud0, ud1, &num0, &den0);

                            double v = num0 / den0; int mp = -1;
                            for (int ns = 0; ns < 12; ns++) {
                                double num = num0, den = den0; int m = 0;
                                if (ru_hi > 0 && ((mode == 1 && cu_up > v) || (mode == 2 && cu_up < v))) { num += ru_hi * cu_up; den += ru_hi; m |= 1; }
                                if (e1_hi > 0 && ((mode == 1 && c1_up > v) || (mode == 2 && c1_up < v))) { num += e1_hi * c1_up; den += e1_hi; m |= 2; }
                                if (e2_hi > 0 && ((mode == 1 && c2_up > v) || (mode == 2 && c2_up < v))) { num += e2_hi * c2_up; den += e2_hi; m |= 4; }
                                b3_tail_var(mode, D0, KT, k1k, k1, 1, r1_lo, r1_hi, 0.5 * k1, 0, 0, f1_0, f1_1, f1_2, f1_1b, f1_2b, c1_bud, uk == 0, d1_0, d1_1, d1_2, cu_up, v, &num, &den, &m, 4);
                                b3_tail_var(mode, D0, KT, k2k, k2, 1, r2_lo, r2_hi, 0.5 * k2, 0, 0, f2_0, f2_1, f2_2, f2_1b, f2_2b, c2_bud, uk == 0, d2_0, d2_1, d2_2, cu_up, v, &num, &den, &m, 10);
                                b3_tail_var(mode, D0, KT, uk, ku, hu, ur_in_lo, ur_in_hi, 0.5 * (hu + ku - 2), 0, 0, uF0, uF1, uF2, uF1b, uF2b, cu_bud, u_nin, ud0, ud1, ud2, uF3_up_sum, v, &num, &den, &m, 16);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                             if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                            hC[b3f_hc_off(&C, ca, w1, w2)] = v;
                        }
                    }
                }
            }
        }

        /* 5a. Solve P4 */
        memset(qP4, 0, C.p4_len * sizeof(double));
        int (*p4_pairs)[2] = (int(*)[2])b3f_xmalloc((size_t)C.NP4 * C.NP4 * sizeof(*p4_pairs));
        for (int rsum = 2 * mr + 2 * (C.DP4 - 1); rsum >= 8; rsum--) for (int ksum = 2 * C.KP4; ksum >= 2; ksum--) {
            int npairs = 0;
            for (int b4 = 0; b4 < C.NP4; b4++) for (int c4 = 0; c4 < C.NP4; c4++)
                if (C.p4_K[b4] + C.p4_K[c4] == ksum) { p4_pairs[npairs][0] = b4; p4_pairs[npairs][1] = c4; npairs++; }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int pi = 0; pi < npairs; pi++) {
                int b4 = p4_pairs[pi][0], c4 = p4_pairs[pi][1];
                int Lb = C.p4_L[b4], kb = C.p4_K[b4], Db = Lb + kb + 1, b_hu = C.p4_hu[b4], b_g = C.p4_g[b4], b_g_up = b3f_next_up(&C, b_g);
                int b4_up = (Db + 1 < C.DP4) ? C.p4_lut[kb][Lb + 1] : -1, b4_bud = (Lb >= 1 && kb + 1 <= C.KP4) ? C.p4_lut[kb + 1][Lb - 1] : -1;
                double ib = 1.0 / Db, bb = (double)Lb / Db, rb_ext = 0.5 * (kb - 1);
                int Lc = C.p4_L[c4], kc = C.p4_K[c4], Dc = Lc + kc + 1, c_hu = C.p4_hu[c4], c_g = C.p4_g[c4], c_g_up = b3f_next_up(&C, c_g);
                int c4_up = (Dc + 1 < C.DP4) ? C.p4_lut[kc][Lc + 1] : -1, c4_bud = (Lc >= 1 && kc + 1 <= C.KP4) ? C.p4_lut[kc + 1][Lc - 1] : -1;
                double ic = 1.0 / Dc, bc = (double)Lc / Dc, rc_ext = 0.5 * (kc - 1);
                int rem = rsum - Db - Dc;
                for (int ra = 2; ra <= mr; ra++) {
                    int rd = rem - ra; if (rd < 2 || rd > mr) continue;
                    for (int ia = 0; ia < C.p_rcnt[ra]; ia++) {
                        int a = C.p_by_r[ra][ia], ak = C.kind[a], La = C.L[a], ka = C.K[a], ha = C.H[a];
                        int a_up = b3f_next_up(&C, a), a_bud = b3f_next_bud(&C, a), a_t1 = (ak >= 3) ? b3f_t1(&C, ka, ha) : a, a_t2 = (ak >= 4) ? b3f_t2(&C, ka, ha) : a;
                        int a_t1b = (ak >= 3 && ka <= KT) ? b3f_t1(&C, ka + 1, ha) : a_t1, a_t2b = (ak >= 4 && ka <= KT) ? b3f_t2(&C, ka + 1, ha) : a_t2;
                        double ra_ext = (ak == 0) ? 0.5 * (ha + ka - 1) : 0.0, inv_a = invD_of(&C, a);
                        double dh_abc = hV[b3f_hv_off(&C, a, b_hu, c_g_up)] - hV[b3f_hv_off(&C, a, b_hu, c_g)];
                        for (int id = 0; id < C.c_rcnt[rd]; id++) {
                            int d = C.c_by_r[rd][id], dk = C.kind[d], Ld = C.L[d], kd = C.K[d];
                            int d_up = b3f_next_up(&C, d), d_bud = b3f_next_bud(&C, d), d_t1 = (dk >= 3) ? b3f_t1(&C, kd, 1) : d, d_t2 = (dk >= 4) ? b3f_t2(&C, kd, 1) : d;
                            int d_t1b = (dk >= 3 && kd <= KT) ? b3f_t1(&C, kd + 1, 1) : d_t1, d_t2b = (dk >= 4 && kd <= KT) ? b3f_t2(&C, kd + 1, 1) : d_t2;
                            double rd_ext = (dk == 0) ? 0.5 * kd : 0.0, inv_d = invD_of(&C, d);
                            double dh_bcd = hV[b3f_hv_off(&C, b_g_up, c_hu, d)] - hV[b3f_hv_off(&C, b_g, c_hu, d)];
                            double qa_up = qP4[b3f_p4_off(&C, a_up, b4, c4, d)], qb_up = (b4_up >= 0) ? qP4[b3f_p4_off(&C, a, b4_up, c4, d)] : 0;
                            double qc_up = (c4_up >= 0) ? qP4[b3f_p4_off(&C, a, b4, c4_up, d)] : 0, qd_up = qP4[b3f_p4_off(&C, a, b4, c4, d_up)];
                            double qa_b = qP4[b3f_p4_off(&C, a_bud, b4, c4, d)], qb_b = (b4_bud >= 0) ? qP4[b3f_p4_off(&C, a, b4_bud, c4, d)] : 0;
                            double qc_b = (c4_bud >= 0) ? qP4[b3f_p4_off(&C, a, b4, c4_bud, d)] : 0, qd_b = qP4[b3f_p4_off(&C, a, b4, c4, d_bud)];
                            size_t bi[4]; int ed[4];
                            l5_couple(l5,(L5Tree){4,{-1,0,1,2},{a,b_g,c_g,d}},iw_nb,bi,ed);
                            qa_b += birth_value(q5,bi[0]); qb_b += birth_value(q5,bi[1]);
                            qc_b += birth_value(q5,bi[2]); qd_b += birth_value(q5,bi[3]);
                            double num0 = (inv_a - rpar[b_g]) * dh_bcd + (inv_d - rch[c_g]) * dh_abc + eps_L4, den0 = lam;
                            if (ak == 0) { double r = La + ib; num0 += r * qa_up; den0 += r; if (La >= 1) { double b = (double)La / C.D[a]; num0 += b * qa_b; den0 += b; } }
                            if (dk == 0) { double r = Ld + ic; num0 += r * qd_up; den0 += r; if (Ld >= 1) { double b = (double)Ld / C.D[d]; num0 += b * qd_b; den0 += b; } }
                            { double r = Lb + (ak == 0 ? inv_a : 0) + ic; num0 += r * qb_up; den0 += r; if (bb > 0) { num0 += bb * qb_b; den0 += bb; } }
                            { double r = Lc + ib + (dk == 0 ? inv_d : 0); num0 += r * qc_up; den0 += r; if (bc > 0) { num0 += bc * qc_b; den0 += bc; } }
                            double aF0 = (ak >= 2) ? qP4[b3f_p4_off(&C, C.INF_ID, b4, c4, d)] : 0, aF1 = (ak >= 3) ? qP4[b3f_p4_off(&C, a_t1, b4, c4, d)] : 0, aF2 = (ak >= 4) ? qP4[b3f_p4_off(&C, a_t2, b4, c4, d)] : 0;
                            double aF1b = (ak >= 3) ? qP4[b3f_p4_off(&C, a_t1b, b4, c4, d)] : qa_b, aF2b = (ak >= 4) ? qP4[b3f_p4_off(&C, a_t2b, b4, c4, d)] : qa_b;
                            double ad0 = (ak >= 2) ? (((b4_up >= 0) ? qP4[b3f_p4_off(&C, C.INF_ID, b4_up, c4, d)] : 0) - aF0) : 0;
                            double ad1 = (ak >= 3) ? (((b4_up >= 0) ? qP4[b3f_p4_off(&C, a_t1, b4_up, c4, d)] : 0) - aF1) : 0;
                            double ad2 = (ak >= 4) ? (((b4_up >= 0) ? qP4[b3f_p4_off(&C, a_t2, b4_up, c4, d)] : 0) - aF2) : 0;
                            double dF0 = (dk >= 2) ? qP4[b3f_p4_off(&C, a, b4, c4, C.INF_ID)] : 0, dF1 = (dk >= 3) ? qP4[b3f_p4_off(&C, a, b4, c4, d_t1)] : 0, dF2 = (dk >= 4) ? qP4[b3f_p4_off(&C, a, b4, c4, d_t2)] : 0;
                            double dF1b = (dk >= 3) ? qP4[b3f_p4_off(&C, a, b4, c4, d_t1b)] : qd_b, dF2b = (dk >= 4) ? qP4[b3f_p4_off(&C, a, b4, c4, d_t2b)] : qd_b;
                            double dd0 = (dk >= 2) ? (((c4_up >= 0) ? qP4[b3f_p4_off(&C, a, b4, c4_up, C.INF_ID)] : 0) - dF0) : 0;
                            double dd1 = (dk >= 3) ? (((c4_up >= 0) ? qP4[b3f_p4_off(&C, a, b4, c4_up, d_t1)] : 0) - dF1) : 0;
                            double dd2 = (dk >= 4) ? (((c4_up >= 0) ? qP4[b3f_p4_off(&C, a, b4, c4_up, d_t2)] : 0) - dF2) : 0;
                            b3_tail_const(mode, D0, KT, ak, ka, ha, ib, ib + 0.5 * (ka + ha - 1), aF0, aF1, (ak == 2) ? qa_b : aF1b, qa_b, 1, ad0, ad1, &num0, &den0);
                            b3_tail_const(mode, D0, KT, dk, kd, 1,  ic, ic + 0.5 * kd,            dF0, dF1, (dk == 2) ? qd_b : dF1b, qd_b, 1, dd0, dd1, &num0, &den0);
                            double v = num0 / den0; int mp = -1;
                            for (int ns = 0; ns < 12; ns++) {
                                double num = num0, den = den0; int m = 0;
                                L5_CORNER(ra_ext,qa_up,0); L5_CORNER(rb_ext,qb_up,1);
                                L5_CORNER(rc_ext,qc_up,2); L5_CORNER(rd_ext,qd_up,3);
                                b3_tail_var(mode, D0, KT, ak, ka, ha, ib, ib, 0.5 * (ka + ha - 1), 0, 0, aF0, aF1, aF2, aF1b, aF2b, qa_b, 1, ad0, ad1, ad2, qb_up, v, &num, &den, &m, 4);
                                b3_tail_var(mode, D0, KT, dk, kd, 1,  ic, ic, 0.5 * kd,            0, 0, dF0, dF1, dF2, dF1b, dF2b, qd_b, 1, dd0, dd1, dd2, qc_up, v, &num, &den, &m, 10);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                             if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                            qP4[b3f_p4_off(&C, a, b4, c4, d)] = v;
                        }
                    }
                }
            }
        }

        free(p4_pairs);

        /* 5b. Solve Ytop */
        memset(qYt, 0, C.yt_len * sizeof(double));
        for (int rsum = 3 * mr + C.DYT - 1; rsum >= 8; rsum--) for (int ku = C.KYT; ku >= 2; ku--) for (int u4 = 0; u4 < C.NYT; u4++) {
            if (C.yt_K[u4] != ku) continue;
            int Lu = C.yt_L[u4], Du = Lu + ku + 1, u_hu = C.yt_hu[u4], u_cu = C.yt_cu[u4], u_g = C.yt_g[u4], u_g_up = b3f_next_up(&C, u_g);
            int u_hu_up = C.hu_of_g[u_g_up], u_cu_up = C.cu_of_g[u_g_up];
            int u4_up = (Du + 1 < C.DYT) ? C.yt_lut[ku][Lu + 1] : -1, u4_bud = (Lu >= 1 && ku + 1 <= C.KYT) ? C.yt_lut[ku + 1][Lu - 1] : -1;
            double iu = 1.0 / Du, bu = (double)Lu / Du, ru_ext = 0.5 * (ku - 2);
            int rem = rsum - Du, tasks[4096][3], nt = 0;
            for (int rp = 2; rp <= mr; rp++) for (int r1 = 2; r1 <= mr; r1++) {
                int r2 = rem - rp - r1;
                if (r2 >= r1 && r2 <= mr && C.p_rcnt[rp] > 0 && C.c_rcnt[r1] > 0 && C.c_rcnt[r2] > 0) {
                    tasks[nt][0] = rp; tasks[nt][1] = r1; tasks[nt][2] = r2; nt++;
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ti = 0; ti < nt; ti++) {
                int rp = tasks[ti][0], r1 = tasks[ti][1], r2 = tasks[ti][2];
                for (int ip = 0; ip < C.p_rcnt[rp]; ip++) {
                    int p = C.p_by_r[rp][ip], pk = C.kind[p], Lp = C.L[p], kp = C.K[p], hp = C.H[p];
                    int p_up = b3f_next_up(&C, p), p_bud = b3f_next_bud(&C, p), p_t1 = (pk >= 3) ? b3f_t1(&C, kp, hp) : p, p_t2 = (pk >= 4) ? b3f_t2(&C, kp, hp) : p;
                    int p_t1b = (pk >= 3 && kp <= KT) ? b3f_t1(&C, kp + 1, hp) : p_t1, p_t2b = (pk >= 4 && kp <= KT) ? b3f_t2(&C, kp + 1, hp) : p_t2;
                    double rp_ext = (pk == 0) ? 0.5 * (hp + kp - 1) : 0.0, inv_p = invD_of(&C, p);
                    for (int i1 = 0; i1 < C.c_rcnt[r1]; i1++) {
                        int w1 = C.c_by_r[r1][i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                        int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1), w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                        int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                        double r1_ext = (k1k == 0) ? 0.5 * k1 : 0.0, inv1 = invD_of(&C, w1);
                        double dhv1 = ((u_hu_up >= 0) ? hV[b3f_hv_off(&C, p, u_hu_up, w1)] : 0) - hV[b3f_hv_off(&C, p, u_hu, w1)];
                        for (int i2 = (r1 == r2 ? i1 : 0); i2 < C.c_rcnt[r2]; i2++) {
                            int w2 = C.c_by_r[r2][i2], k2k = C.kind[w2], L2 = C.L[w2], k2 = C.K[w2];
                            int w2_up = b3f_next_up(&C, w2), w2_bud = b3f_next_bud(&C, w2), w2_t1 = (k2k >= 3) ? b3f_t1(&C, k2, 1) : w2, w2_t2 = (k2k >= 4) ? b3f_t2(&C, k2, 1) : w2;
                            int w2_t1b = (k2k >= 3 && k2 <= KT) ? b3f_t1(&C, k2 + 1, 1) : w2_t1, w2_t2b = (k2k >= 4 && k2 <= KT) ? b3f_t2(&C, k2 + 1, 1) : w2_t2;
                            double r2_ext = (k2k == 0) ? 0.5 * k2 : 0.0, inv2 = invD_of(&C, w2);
                            double dhv2 = ((u_hu_up >= 0) ? hV[b3f_hv_off(&C, p, u_hu_up, w2)] : 0) - hV[b3f_hv_off(&C, p, u_hu, w2)];
                            double dhc  = ((u_cu_up >= 0) ? hC[b3f_hc_off(&C, u_cu_up, w1, w2)] : 0) - hC[b3f_hc_off(&C, u_cu, w1, w2)];
                            double qp_up = qYt[b3f_yt_off(&C, p_up, u4, w1, w2)], qu_up = (u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, w1, w2)] : 0;
                            double q1_up = qYt[b3f_yt_off(&C, p, u4, w1_up, w2)], q2_up = qYt[b3f_yt_off(&C, p, u4, w1, w2_up)];
                            double qp_b = qYt[b3f_yt_off(&C, p_bud, u4, w1, w2)], qu_b = (u4_bud >= 0) ? qYt[b3f_yt_off(&C, p, u4_bud, w1, w2)] : 0;
                            double q1_b = qYt[b3f_yt_off(&C, p, u4, w1_bud, w2)], q2_b = qYt[b3f_yt_off(&C, p, u4, w1, w2_bud)];
                            size_t bi[4]; int ed[4];
                            l5_couple(l5,(L5Tree){4,{-1,0,1,1},{p,u_g,w1,w2}},iw_nb,bi,ed);
                            qp_b += birth_value(q5,bi[0]); qu_b += birth_value(q5,bi[1]);
                            q1_b += birth_value(q5,bi[2]); q2_b += birth_value(q5,bi[3]);
                            double num0 = (inv_p - rpar[u_g]) * dhc + (inv2 - rch[u_g]) * dhv1 + (inv1 - rch[u_g]) * dhv2 + eps_L4, den0 = lam;
                            if (pk == 0)  { double r = Lp + iu; num0 += r * qp_up; den0 += r; if (Lp >= 1) { double b = (double)Lp / C.D[p]; num0 += b * qp_b; den0 += b; } }
                            if (k1k == 0) { double r = L1 + iu; num0 += r * q1_up; den0 += r; if (L1 >= 1) { double b = (double)L1 / C.D[w1]; num0 += b * q1_b; den0 += b; } }
                            if (k2k == 0) { double r = L2 + iu; num0 += r * q2_up; den0 += r; if (L2 >= 1) { double b = (double)L2 / C.D[w2]; num0 += b * q2_b; den0 += b; } }
                            { double r = Lu + (pk == 0 ? inv_p : 0) + (k1k == 0 ? inv1 : 0) + (k2k == 0 ? inv2 : 0); num0 += r * qu_up; den0 += r; if (bu > 0) { num0 += bu * qu_b; den0 += bu; } }
                            double pF0 = (pk >= 2) ? qYt[b3f_yt_off(&C, C.INF_ID, u4, w1, w2)] : 0, pF1 = (pk >= 3) ? qYt[b3f_yt_off(&C, p_t1, u4, w1, w2)] : 0, pF2 = (pk >= 4) ? qYt[b3f_yt_off(&C, p_t2, u4, w1, w2)] : 0;
                            double pF1b = (pk >= 3) ? qYt[b3f_yt_off(&C, p_t1b, u4, w1, w2)] : qp_b, pF2b = (pk >= 4) ? qYt[b3f_yt_off(&C, p_t2b, u4, w1, w2)] : qp_b;
                            double pd0 = (pk >= 2) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, C.INF_ID, u4_up, w1, w2)] : 0) - pF0) : 0;
                            double pd1 = (pk >= 3) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p_t1, u4_up, w1, w2)] : 0) - pF1) : 0;
                            double pd2 = (pk >= 4) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p_t2, u4_up, w1, w2)] : 0) - pF2) : 0;
                            double f10 = (k1k >= 2) ? qYt[b3f_yt_off(&C, p, u4, C.INF_ID, w2)] : 0, f11 = (k1k >= 3) ? qYt[b3f_yt_off(&C, p, u4, w1_t1, w2)] : 0, f12 = (k1k >= 4) ? qYt[b3f_yt_off(&C, p, u4, w1_t2, w2)] : 0;
                            double f11b = (k1k >= 3) ? qYt[b3f_yt_off(&C, p, u4, w1_t1b, w2)] : q1_b, f12b = (k1k >= 4) ? qYt[b3f_yt_off(&C, p, u4, w1_t2b, w2)] : q1_b;
                            double d10 = (k1k >= 2) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, C.INF_ID, w2)] : 0) - f10) : 0;
                            double d11 = (k1k >= 3) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, w1_t1, w2)] : 0) - f11) : 0;
                            double d12 = (k1k >= 4) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, w1_t2, w2)] : 0) - f12) : 0;
                            double f20 = (k2k >= 2) ? qYt[b3f_yt_off(&C, p, u4, w1, C.INF_ID)] : 0, f21 = (k2k >= 3) ? qYt[b3f_yt_off(&C, p, u4, w1, w2_t1)] : 0, f22 = (k2k >= 4) ? qYt[b3f_yt_off(&C, p, u4, w1, w2_t2)] : 0;
                            double f21b = (k2k >= 3) ? qYt[b3f_yt_off(&C, p, u4, w1, w2_t1b)] : q2_b, f22b = (k2k >= 4) ? qYt[b3f_yt_off(&C, p, u4, w1, w2_t2b)] : q2_b;
                            double d20 = (k2k >= 2) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, w1, C.INF_ID)] : 0) - f20) : 0;
                            double d21 = (k2k >= 3) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, w1, w2_t1)] : 0) - f21) : 0;
                            double d22 = (k2k >= 4) ? (((u4_up >= 0) ? qYt[b3f_yt_off(&C, p, u4_up, w1, w2_t2)] : 0) - f22) : 0;
                            b3_tail_const(mode, D0, KT, pk, kp, hp, iu, iu + 0.5 * (kp + hp - 1), pF0, pF1, (pk == 2) ? qp_b : pF1b, qp_b, 1, pd0, pd1, &num0, &den0);
                            b3_tail_const(mode, D0, KT, k1k, k1, 1, iu, iu + 0.5 * k1,            f10, f11, (k1k == 2) ? q1_b : f11b, q1_b, 1, d10, d11, &num0, &den0);
                            b3_tail_const(mode, D0, KT, k2k, k2, 1, iu, iu + 0.5 * k2,            f20, f21, (k2k == 2) ? q2_b : f21b, q2_b, 1, d20, d21, &num0, &den0);
                            double v = num0 / den0; int mp = -1;
                            for (int ns = 0; ns < 12; ns++) {
                                double num = num0, den = den0; int m = 0;
                                L5_CORNER(rp_ext,qp_up,0); L5_CORNER(ru_ext,qu_up,1);
                                L5_CORNER(r1_ext,q1_up,2); L5_CORNER(r2_ext,q2_up,3);
                                b3_tail_var(mode, D0, KT, pk, kp, hp, iu, iu, 0.5 * (kp + hp - 1), 0, 0, pF0, pF1, pF2, pF1b, pF2b, qp_b, 1, pd0, pd1, pd2, qu_up, v, &num, &den, &m, 4);
                                b3_tail_var(mode, D0, KT, k1k, k1, 1, iu, iu, 0.5 * k1,            0, 0, f10, f11, f12, f11b, f12b, q1_b, 1, d10, d11, d12, qu_up, v, &num, &den, &m, 10);
                                b3_tail_var(mode, D0, KT, k2k, k2, 1, iu, iu, 0.5 * k2,            0, 0, f20, f21, f22, f21b, f22b, q2_b, 1, d20, d21, d22, qu_up, v, &num, &den, &m, 16);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                             if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                            qYt[b3f_yt_off(&C, p, u4, w1, w2)] = v;
                        }
                    }
                }
            }
        }

        /* 5c. Solve Ybot */
        memset(qYb, 0, C.yb_len * sizeof(double));
        int (*yb_pairs)[2] = (int(*)[2])b3f_xmalloc((size_t)C.NYBU * C.NYBW * sizeof(*yb_pairs));
        for (int rsum = 2 * mr + 2 * (C.DYB - 1); rsum >= 8; rsum--) for (int ksum = 2 * C.KYB; ksum >= 3; ksum--) {
            int npairs = 0;
            for (int u4 = 0; u4 < C.NYBU; u4++) for (int w2_4 = 0; w2_4 < C.NYBW; w2_4++)
                if (C.ybu_K[u4] + C.ybw_K[w2_4] == ksum) { yb_pairs[npairs][0] = u4; yb_pairs[npairs][1] = w2_4; npairs++; }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int pi = 0; pi < npairs; pi++) {
                int u4 = yb_pairs[pi][0], w2_4 = yb_pairs[pi][1];
                int Lu = C.ybu_L[u4], ku = C.ybu_K[u4], hu = C.ybu_H[u4], Du = Lu + ku + hu, u_cu = C.ybu_cu[u4], u_g = C.ybu_g[u4], u_g_up = b3f_next_up(&C, u_g);
                int u4_up = (Du + 1 < C.DYB) ? C.ybu_lut[hu][ku][Lu + 1] : -1, u4_bud = (Lu >= 1 && ku + 1 <= C.KYB) ? C.ybu_lut[hu][ku + 1][Lu - 1] : -1;
                double iu = 1.0 / Du, bu = (double)Lu / Du, ru_ext = 0.5 * (hu + ku - 2);
                int L2 = C.ybw_L[w2_4], k2 = C.ybw_K[w2_4], D2 = L2 + k2 + 1, w2_hu = C.ybw_hu[w2_4], w2_g = C.ybw_g[w2_4], w2_g_up = b3f_next_up(&C, w2_g);
                int w2_up = (D2 + 1 < C.DYB) ? C.ybw_lut[k2][L2 + 1] : -1, w2_bud = (L2 >= 1 && k2 + 1 <= C.KYB) ? C.ybw_lut[k2 + 1][L2 - 1] : -1;
                double i2 = 1.0 / D2, b2 = (double)L2 / D2, r2_ext = 0.5 * (k2 - 1);
                int rem = rsum - Du - D2;
                for (int r1 = 2; r1 <= mr; r1++) {
                    int rc = rem - r1; if (rc < 2 || rc > mr) continue;
                    for (int i1 = 0; i1 < C.c_rcnt[r1]; i1++) {
                        int w1 = C.c_by_r[r1][i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                        int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1), w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                        int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                        double r1_ext = (k1k == 0) ? 0.5 * k1 : 0.0, inv1 = invD_of(&C, w1);
                        double dhc = hC[b3f_hc_off(&C, u_cu, w1, w2_g_up)] - hC[b3f_hc_off(&C, u_cu, w1, w2_g)];
                        for (int ic = 0; ic < C.c_rcnt[rc]; ic++) {
                            int c = C.c_by_r[rc][ic], ck = C.kind[c], Lc = C.L[c], kc = C.K[c];
                            int c_up = b3f_next_up(&C, c), c_bud = b3f_next_bud(&C, c), c_t1 = (ck >= 3) ? b3f_t1(&C, kc, 1) : c, c_t2 = (ck >= 4) ? b3f_t2(&C, kc, 1) : c;
                            int c_t1b = (ck >= 3 && kc <= KT) ? b3f_t1(&C, kc + 1, 1) : c_t1, c_t2b = (ck >= 4 && kc <= KT) ? b3f_t2(&C, kc + 1, 1) : c_t2;
                            double rc_ext = (ck == 0) ? 0.5 * kc : 0.0, inv_c = invD_of(&C, c);
                            double dhv = hV[b3f_hv_off(&C, u_g_up, w2_hu, c)] - hV[b3f_hv_off(&C, u_g, w2_hu, c)];
                            double qu_up = (u4_up >= 0) ? qYb[b3f_yb_off(&C, u4_up, w2_4, w1, c)] : 0, q2_up = (w2_up >= 0) ? qYb[b3f_yb_off(&C, u4, w2_up, w1, c)] : 0;
                            double q1_up = qYb[b3f_yb_off(&C, u4, w2_4, w1_up, c)], qc_up = qYb[b3f_yb_off(&C, u4, w2_4, w1, c_up)];
                            double qu_b = (u4_bud >= 0) ? qYb[b3f_yb_off(&C, u4_bud, w2_4, w1, c)] : 0, q2_b = (w2_bud >= 0) ? qYb[b3f_yb_off(&C, u4, w2_bud, w1, c)] : 0;
                            double q1_b = qYb[b3f_yb_off(&C, u4, w2_4, w1_bud, c)], qc_b = qYb[b3f_yb_off(&C, u4, w2_4, w1, c_bud)];
                            size_t bi[4]; int ed[4];
                            l5_couple(l5,(L5Tree){4,{-1,0,0,2},{u_g,w1,w2_g,c}},iw_nb,bi,ed);
                            qu_b += birth_value(q5,bi[0]); q1_b += birth_value(q5,bi[1]);
                            q2_b += birth_value(q5,bi[2]); qc_b += birth_value(q5,bi[3]);
                            double num0 = (inv1 - rch[u_g]) * dhv + (inv_c - rch[w2_g]) * dhc + eps_L4, den0 = lam;
                            if (k1k == 0) { double r = L1 + iu; num0 += r * q1_up; den0 += r; if (L1 >= 1) { double b = (double)L1 / C.D[w1]; num0 += b * q1_b; den0 += b; } }
                            if (ck == 0)  { double r = Lc + i2; num0 += r * qc_up; den0 += r; if (Lc >= 1) { double b = (double)Lc / C.D[c]; num0 += b * qc_b; den0 += b; } }
                            { double r = Lu + (k1k == 0 ? inv1 : 0) + i2; num0 += r * qu_up; den0 += r; if (bu > 0) { num0 += bu * qu_b; den0 += bu; } }
                            { double r = L2 + iu + (ck == 0 ? inv_c : 0); num0 += r * q2_up; den0 += r; if (b2 > 0) { num0 += b2 * q2_b; den0 += b2; } }
                            double f10 = (k1k >= 2) ? qYb[b3f_yb_off(&C, u4, w2_4, C.INF_ID, c)] : 0, f11 = (k1k >= 3) ? qYb[b3f_yb_off(&C, u4, w2_4, w1_t1, c)] : 0, f12 = (k1k >= 4) ? qYb[b3f_yb_off(&C, u4, w2_4, w1_t2, c)] : 0;
                            double f11b = (k1k >= 3) ? qYb[b3f_yb_off(&C, u4, w2_4, w1_t1b, c)] : q1_b, f12b = (k1k >= 4) ? qYb[b3f_yb_off(&C, u4, w2_4, w1_t2b, c)] : q1_b;
                            double d10 = (k1k >= 2) ? (((u4_up >= 0) ? qYb[b3f_yb_off(&C, u4_up, w2_4, C.INF_ID, c)] : 0) - f10) : 0;
                            double d11 = (k1k >= 3) ? (((u4_up >= 0) ? qYb[b3f_yb_off(&C, u4_up, w2_4, w1_t1, c)] : 0) - f11) : 0;
                            double d12 = (k1k >= 4) ? (((u4_up >= 0) ? qYb[b3f_yb_off(&C, u4_up, w2_4, w1_t2, c)] : 0) - f12) : 0;
                            double fc0 = (ck >= 2) ? qYb[b3f_yb_off(&C, u4, w2_4, w1, C.INF_ID)] : 0, fc1 = (ck >= 3) ? qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t1)] : 0, fc2 = (ck >= 4) ? qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t2)] : 0;
                            double fc1b = (ck >= 3) ? qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t1b)] : qc_b, fc2b = (ck >= 4) ? qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t2b)] : qc_b;
                            double dc0 = (ck >= 2) ? (((w2_up >= 0) ? qYb[b3f_yb_off(&C, u4, w2_up, w1, C.INF_ID)] : 0) - fc0) : 0;
                            double dc1 = (ck >= 3) ? (((w2_up >= 0) ? qYb[b3f_yb_off(&C, u4, w2_up, w1, c_t1)] : 0) - fc1) : 0;
                            double dc2 = (ck >= 4) ? (((w2_up >= 0) ? qYb[b3f_yb_off(&C, u4, w2_up, w1, c_t2)] : 0) - fc2) : 0;
                            b3_tail_const(mode, D0, KT, k1k, k1, 1, iu, iu + 0.5 * k1, f10, f11, (k1k == 2) ? q1_b : f11b, q1_b, 1, d10, d11, &num0, &den0);
                            b3_tail_const(mode, D0, KT, ck,  kc, 1, i2, i2 + 0.5 * kc, fc0, fc1, (ck == 2) ? qc_b : fc1b, qc_b, 1, dc0, dc1, &num0, &den0);
                            double v = num0 / den0; int mp = -1;
                            for (int ns = 0; ns < 12; ns++) {
                                double num = num0, den = den0; int m = 0;
                                L5_CORNER(ru_ext,qu_up,0); L5_CORNER(r1_ext,q1_up,1);
                                L5_CORNER(r2_ext,q2_up,2); L5_CORNER(rc_ext,qc_up,3);
                                b3_tail_var(mode, D0, KT, k1k, k1, 1, iu, iu, 0.5 * k1, 0, 0, f10, f11, f12, f11b, f12b, q1_b, 1, d10, d11, d12, qu_up, v, &num, &den, &m, 4);
                                b3_tail_var(mode, D0, KT, ck,  kc, 1, i2, i2, 0.5 * kc, 0, 0, fc0, fc1, fc2, fc1b, fc2b, qc_b, 1, dc0, dc1, dc2, q2_up, v, &num, &den, &m, 10);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                             if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                            qYb[b3f_yb_off(&C, u4, w2_4, w1, c)] = v;
                        }
                    }
                }
            }
        }

        free(yb_pairs);

        /* 5d. Solve S3 */
        memset(qS3, 0, C.s3_len * sizeof(double));
        for (int rsum = 3 * mr + C.DS3 - 1; rsum >= 8; rsum--) for (int ku = C.KS3; ku >= 3; ku--) for (int u4 = 0; u4 < C.NS3; u4++) {
            if (C.s3_K[u4] != ku) continue;
            int Lu = C.s3_L[u4], hu = C.s3_H[u4], Du = Lu + ku + hu, u_cu = C.s3_cu[u4], u_g = C.s3_g[u4], u_cu_up = C.cu_of_g[b3f_next_up(&C, u_g)];
            int u4_up = (Du + 1 < C.DS3) ? C.s3_lut[hu][ku][Lu + 1] : -1, u4_bud = (Lu >= 1 && ku + 1 <= C.KS3) ? C.s3_lut[hu][ku + 1][Lu - 1] : -1;
            double iu = 1.0 / Du, bu = (double)Lu / Du, ru_ext = 0.5 * (hu + ku - 3);
            int rem = rsum - Du, tasks[4096][3], nt = 0;
            for (int r1 = 2; r1 <= mr; r1++) for (int r2 = r1; r2 <= mr; r2++) {
                int r3 = rem - r1 - r2;
                if (r3 >= r2 && r3 <= mr && C.c_rcnt[r1] > 0 && C.c_rcnt[r2] > 0 && C.c_rcnt[r3] > 0) {
                    tasks[nt][0] = r1; tasks[nt][1] = r2; tasks[nt][2] = r3; nt++;
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ti = 0; ti < nt; ti++) {
                int r1 = tasks[ti][0], r2 = tasks[ti][1], r3 = tasks[ti][2];
                for (int i1 = 0; i1 < C.c_rcnt[r1]; i1++) {
                    int w1 = C.c_by_r[r1][i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                    int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1), w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                    int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                    double r1_ext = (k1k == 0) ? 0.5 * k1 : 0.0, inv1 = invD_of(&C, w1);
                    for (int i2 = (r1 == r2 ? i1 : 0); i2 < C.c_rcnt[r2]; i2++) {
                        int w2 = C.c_by_r[r2][i2], k2k = C.kind[w2], L2 = C.L[w2], k2 = C.K[w2];
                        int w2_up = b3f_next_up(&C, w2), w2_bud = b3f_next_bud(&C, w2), w2_t1 = (k2k >= 3) ? b3f_t1(&C, k2, 1) : w2, w2_t2 = (k2k >= 4) ? b3f_t2(&C, k2, 1) : w2;
                        int w2_t1b = (k2k >= 3 && k2 <= KT) ? b3f_t1(&C, k2 + 1, 1) : w2_t1, w2_t2b = (k2k >= 4 && k2 <= KT) ? b3f_t2(&C, k2 + 1, 1) : w2_t2;
                        double r2_ext = (k2k == 0) ? 0.5 * k2 : 0.0, inv2 = invD_of(&C, w2);
                        double dh12 = ((u_cu_up >= 0) ? hC[b3f_hc_off(&C, u_cu_up, w1, w2)] : 0) - hC[b3f_hc_off(&C, u_cu, w1, w2)];
                        for (int i3 = (r2 == r3 ? i2 : 0); i3 < C.c_rcnt[r3]; i3++) {
                            int w3 = C.c_by_r[r3][i3], k3k = C.kind[w3], L3 = C.L[w3], k3 = C.K[w3];
                            int w3_up = b3f_next_up(&C, w3), w3_bud = b3f_next_bud(&C, w3), w3_t1 = (k3k >= 3) ? b3f_t1(&C, k3, 1) : w3, w3_t2 = (k3k >= 4) ? b3f_t2(&C, k3, 1) : w3;
                            int w3_t1b = (k3k >= 3 && k3 <= KT) ? b3f_t1(&C, k3 + 1, 1) : w3_t1, w3_t2b = (k3k >= 4 && k3 <= KT) ? b3f_t2(&C, k3 + 1, 1) : w3_t2;
                            double r3_ext = (k3k == 0) ? 0.5 * k3 : 0.0, inv3 = invD_of(&C, w3);
                            double dh13 = ((u_cu_up >= 0) ? hC[b3f_hc_off(&C, u_cu_up, w1, w3)] : 0) - hC[b3f_hc_off(&C, u_cu, w1, w3)];
                            double dh23 = ((u_cu_up >= 0) ? hC[b3f_hc_off(&C, u_cu_up, w2, w3)] : 0) - hC[b3f_hc_off(&C, u_cu, w2, w3)];
                            double qu_up = (u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, w2, w3)] : 0, qu_b = (u4_bud >= 0) ? qS3[b3f_s3_off(&C, u4_bud, w1, w2, w3)] : 0;
                            double q1_up = qS3[b3f_s3_off(&C, u4, w1_up, w2, w3)], q2_up = qS3[b3f_s3_off(&C, u4, w1, w2_up, w3)], q3_up = qS3[b3f_s3_off(&C, u4, w1, w2, w3_up)];
                            double q1_b = qS3[b3f_s3_off(&C, u4, w1_bud, w2, w3)], q2_b = qS3[b3f_s3_off(&C, u4, w1, w2_bud, w3)], q3_b = qS3[b3f_s3_off(&C, u4, w1, w2, w3_bud)];
                            size_t bi[4]; int ed[4];
                            l5_couple(l5,(L5Tree){4,{-1,0,0,0},{u_g,w1,w2,w3}},iw_nb,bi,ed);
                            qu_b += birth_value(q5,bi[0]); q1_b += birth_value(q5,bi[1]);
                            q2_b += birth_value(q5,bi[2]); q3_b += birth_value(q5,bi[3]);
                            double num0 = (inv3 - rch[u_g]) * dh12 + (inv2 - rch[u_g]) * dh13 + (inv1 - rch[u_g]) * dh23 + eps_L4, den0 = lam;
                            if (k1k == 0) { double r = L1 + iu; num0 += r * q1_up; den0 += r; if (L1 >= 1) { double b = (double)L1 / C.D[w1]; num0 += b * q1_b; den0 += b; } }
                            if (k2k == 0) { double r = L2 + iu; num0 += r * q2_up; den0 += r; if (L2 >= 1) { double b = (double)L2 / C.D[w2]; num0 += b * q2_b; den0 += b; } }
                            if (k3k == 0) { double r = L3 + iu; num0 += r * q3_up; den0 += r; if (L3 >= 1) { double b = (double)L3 / C.D[w3]; num0 += b * q3_b; den0 += b; } }
                            { double r = Lu + (k1k == 0 ? inv1 : 0) + (k2k == 0 ? inv2 : 0) + (k3k == 0 ? inv3 : 0); num0 += r * qu_up; den0 += r; if (bu > 0) { num0 += bu * qu_b; den0 += bu; } }
                            double f10 = (k1k >= 2) ? qS3[b3f_s3_off(&C, u4, C.INF_ID, w2, w3)] : 0, f11 = (k1k >= 3) ? qS3[b3f_s3_off(&C, u4, w1_t1, w2, w3)] : 0, f12 = (k1k >= 4) ? qS3[b3f_s3_off(&C, u4, w1_t2, w2, w3)] : 0;
                            double f11b = (k1k >= 3) ? qS3[b3f_s3_off(&C, u4, w1_t1b, w2, w3)] : q1_b, f12b = (k1k >= 4) ? qS3[b3f_s3_off(&C, u4, w1_t2b, w2, w3)] : q1_b;
                            double d10 = (k1k >= 2) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, C.INF_ID, w2, w3)] : 0) - f10) : 0;
                            double d11 = (k1k >= 3) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1_t1, w2, w3)] : 0) - f11) : 0;
                            double d12 = (k1k >= 4) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1_t2, w2, w3)] : 0) - f12) : 0;
                            double f20 = (k2k >= 2) ? qS3[b3f_s3_off(&C, u4, w1, C.INF_ID, w3)] : 0, f21 = (k2k >= 3) ? qS3[b3f_s3_off(&C, u4, w1, w2_t1, w3)] : 0, f22 = (k2k >= 4) ? qS3[b3f_s3_off(&C, u4, w1, w2_t2, w3)] : 0;
                            double f21b = (k2k >= 3) ? qS3[b3f_s3_off(&C, u4, w1, w2_t1b, w3)] : q2_b, f22b = (k2k >= 4) ? qS3[b3f_s3_off(&C, u4, w1, w2_t2b, w3)] : q2_b;
                            double d20 = (k2k >= 2) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, C.INF_ID, w3)] : 0) - f20) : 0;
                            double d21 = (k2k >= 3) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, w2_t1, w3)] : 0) - f21) : 0;
                            double d22 = (k2k >= 4) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, w2_t2, w3)] : 0) - f22) : 0;
                            double f30 = (k3k >= 2) ? qS3[b3f_s3_off(&C, u4, w1, w2, C.INF_ID)] : 0, f31 = (k3k >= 3) ? qS3[b3f_s3_off(&C, u4, w1, w2, w3_t1)] : 0, f32 = (k3k >= 4) ? qS3[b3f_s3_off(&C, u4, w1, w2, w3_t2)] : 0;
                            double f31b = (k3k >= 3) ? qS3[b3f_s3_off(&C, u4, w1, w2, w3_t1b)] : q3_b, f32b = (k3k >= 4) ? qS3[b3f_s3_off(&C, u4, w1, w2, w3_t2b)] : q3_b;
                            double d30 = (k3k >= 2) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, w2, C.INF_ID)] : 0) - f30) : 0;
                            double d31 = (k3k >= 3) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, w2, w3_t1)] : 0) - f31) : 0;
                            double d32 = (k3k >= 4) ? (((u4_up >= 0) ? qS3[b3f_s3_off(&C, u4_up, w1, w2, w3_t2)] : 0) - f32) : 0;
                            b3_tail_const(mode, D0, KT, k1k, k1, 1, iu, iu + 0.5 * k1, f10, f11, (k1k == 2) ? q1_b : f11b, q1_b, 1, d10, d11, &num0, &den0);
                            b3_tail_const(mode, D0, KT, k2k, k2, 1, iu, iu + 0.5 * k2, f20, f21, (k2k == 2) ? q2_b : f21b, q2_b, 1, d20, d21, &num0, &den0);
                            b3_tail_const(mode, D0, KT, k3k, k3, 1, iu, iu + 0.5 * k3, f30, f31, (k3k == 2) ? q3_b : f31b, q3_b, 1, d30, d31, &num0, &den0);
                            double v = num0 / den0; int mp = -1;
                            for (int ns = 0; ns < 12; ns++) {
                                double num = num0, den = den0; int m = 0;
                                L5_CORNER(ru_ext,qu_up,0); L5_CORNER(r1_ext,q1_up,1);
                                L5_CORNER(r2_ext,q2_up,2); L5_CORNER(r3_ext,q3_up,3);
                                b3_tail_var(mode, D0, KT, k1k, k1, 1, iu, iu, 0.5 * k1, 0, 0, f10, f11, f12, f11b, f12b, q1_b, 1, d10, d11, d12, qu_up, v, &num, &den, &m, 4);
                                b3_tail_var(mode, D0, KT, k2k, k2, 1, iu, iu, 0.5 * k2, 0, 0, f20, f21, f22, f21b, f22b, q2_b, 1, d20, d21, d22, qu_up, v, &num, &den, &m, 10);
                                b3_tail_var(mode, D0, KT, k3k, k3, 1, iu, iu, 0.5 * k3, 0, 0, f30, f31, f32, f31b, f32b, q3_b, 1, d30, d31, d32, qu_up, v, &num, &den, &m, 16);
                                double vn = num / den; if (m == mp) { v = vn; break; } mp = m; v = vn;
                            }
                             if (final_phase) v = (double)quantize(v,H.S) / (double)H.S;
                            qS3[b3f_s3_off(&C, u4, w1, w2, w3)] = v;
                        }
                    }
                }
            }
        }

        double start_l5 = omp_get_wtime();
        l5_sweep((L5Solve){l5,{qP4,qYt,qYb,qS3},q5,lam,eps_L4,mode,NULL},nthreads);
        printf("  sweep=%d lambda=%.12g L5_wall=%.3fs\n",it+1,lam,omp_get_wtime()-start_l5);
        fflush(stdout);
        double rho = fn[v_nb] / f[v_nb];
        if (!final_phase) {
            if (it >= 2) lam += 0.55 * (rho - 1.0) * lam;
            double inv_s = 1.0 / fn[v_nb];
            #pragma omp parallel for num_threads(nthreads)
            for (int i = 0; i < C.NV; i++) f[i] = fn[i] * inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.g_len; i++) g[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.hv_len; i++) hV[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.hc_len; i++) hC[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.p4_len; i++) qP4[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.yt_len; i++) qYt[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.yb_len; i++) qYb[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < C.s3_len; i++) qS3[i] *= inv_s;
            #pragma omp parallel for num_threads(nthreads)
            for (size_t i = 0; i < l5->rows; ++i) q5[i] *= inv_s;
        } else {
            /* Do NOT renormalize in final_phase: check monotone contraction/expansion toward strict feasibility */
            printf("  final_it=%d lam=%.10f f[v_nb]=%.10f->%.10f (ratio=%.10f)\n", it - nit, lam, f[v_nb], fn[v_nb], rho);
            fflush(stdout);
            memcpy(f, fn, C.NV * sizeof(double));
        }
    }

    FILE *fp = fopen(out_path, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", out_path); return 1; }
    if (fwrite(&H, sizeof(H), 1, fp) != 1 || !l5_header_io(fp,&E,1) ||
        fwrite(rpar_q,  sizeof(int64_t), C.INNER_N, fp) != (size_t)C.INNER_N ||
        fwrite(rch_q,   sizeof(int64_t), C.INNER_N, fp) != (size_t)C.INNER_N ||
        fwrite(trpar_q, sizeof(int64_t), C.NT, fp) != (size_t)C.NT ||
        fwrite(trch_q,  sizeof(int64_t), C.NT, fp) != (size_t)C.NT) {
        fprintf(stderr, "Header/rate fwrite failed\n"); return 1;
    }
    write_i64_array(fp, c2_f,   C.NT,     H.S);
    write_i64_array(fp, c3_f,   C.NT,     H.S);
    write_i64_array(fp, Yp_ref, C.yp_len, H.S);
    write_i64_array(fp, Yw_ref, C.yw_len, H.S);
    write_i64_array(fp, f,      C.NV,     H.S);
    write_i64_array(fp, g,      C.g_len,  H.S);
    write_i64_array(fp, hV,     C.hv_len, H.S);
    write_i64_array(fp, hC,     C.hc_len, H.S);
    write_i64_array(fp, qP4,    C.p4_len, H.S);
    write_i64_array(fp, qYt,    C.yt_len, H.S);
    write_i64_array(fp, qYb,    C.yb_len, H.S);
    write_i64_array(fp, qS3,    C.s3_len, H.S);
    write_i64_array(fp, q5, l5->rows, H.S);
    if (fclose(fp) != 0) { fprintf(stderr, "fclose failed\n"); return 1; }
    printf("Wrote %s (%zu bytes), exact P/Q = %lld/%lld (%.10f)\n",
           out_path, bytes, (long long)H.P, (long long)H.Q, (double)H.P / (double)H.Q);

    free(rpar_q); free(rch_q); free(trpar_q); free(trch_q);
    free(rpar); free(rch); free(trpar); free(trch);
    free(f); free(fn); free(g); free(hV); free(hC);
    free(qP4); free(qYt); free(qYb); free(qS3);
    free(dLf); free(c2_f); free(c3_f); free(Yp_ref); free(Yw_ref);
    free(q5); l5_free(l5);
    b3f_layout_free(&C);
    return 0;
}
