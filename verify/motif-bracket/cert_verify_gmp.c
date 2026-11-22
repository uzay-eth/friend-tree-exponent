/* 100% Exact Arbitrary-Precision GMP Verifier for Analytic 1/D-Tail Certificates
 *
 * ZERO floating-point arithmetic is used in ANY verification inequality!
 * Every single inequality across all ~17 million states (f, g, h_V, h_C),
 * global pointwise comparison (Phi(T) >= c0*C(T) for UPPER, <= C0*C(T) for LOWER),
 * and brute-force finite-tree exact-rational drift checks (small trees + extreme tail trees)
 * is checked in exact GMP mpz_t / mpq_t arithmetic.
 *
 * Also supports --override-pq P Q and --perturb <kind> for built-in negative controls.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <gmp.h>
#include <omp.h>
#include "cert_format.h"

_Static_assert(sizeof(long) == 8, "Requires 64-bit signed long for GMP *_si");

/* Exact rational evaluation of piecewise-linear function A(x) at x = x_num / x_den:
 * Returns exact mpq_t out = A(x_num/x_den), where knot values are k_q[j] / S
 * and knot abscissae are knot_num[j] / Q_K.
 */
static void eval_pl_q(mpq_t out, const CertHeader *H, const int64_t *k_q, long x_num, long x_den) {
    /* Compare x_num/x_den against knot_num[j]/Q_K via x_num*Q_K vs knot_num[j]*x_den */
    long xQ = x_num * H->Q_K;
    int m = H->nknots;
    if (xQ <= (long)H->knot_num[0] * x_den) {
        mpq_set_si(out, (long)k_q[0], (unsigned long)H->S);
        mpq_canonicalize(out);
        return;
    }
    if (xQ >= (long)H->knot_num[m - 1] * x_den) {
        mpq_set_si(out, (long)k_q[m - 1], (unsigned long)H->S);
        mpq_canonicalize(out);
        return;
    }
    for (int i = 0; i < m - 1; i++) {
        long k0d = (long)H->knot_num[i] * x_den;
        long k1d = (long)H->knot_num[i + 1] * x_den;
        if (xQ <= k1d) {
            /* t = (xQ - k0d) / (k1d - k0d), 1 - t = (k1d - xQ) / (k1d - k0d) */
            mpz_t num, den, tmp;
            mpz_inits(num, den, tmp, NULL);
            mpz_set_si(num, (long)k_q[i]);
            mpz_mul_si(num, num, k1d - xQ);
            mpz_set_si(tmp, (long)k_q[i + 1]);
            mpz_mul_si(tmp, tmp, xQ - k0d);
            mpz_add(num, num, tmp);
            mpz_set_si(den, (long)H->S);
            mpz_mul_si(den, den, k1d - k0d);
            mpz_set(mpq_numref(out), num);
            mpz_set(mpq_denref(out), den);
            mpq_canonicalize(out);
            mpz_clears(num, den, tmp, NULL);
            return;
        }
    }
}

/* Exact rational eval_g(p -> w) for arbitrary (Lp, kp, hp, Lw, kw) */
static void eval_g_q(mpq_t out, const CertHeader *H, const CertLayout *C, const int64_t *g_q,
                     int Lp, int kp, int hp, int Lw, int kw) {
    int Dp = Lp + kp + hp, Dw = Lw + kw + 1, D0 = C->D0, N = C->NG;
    if (Dp >= D0 && Dw >= D0) { mpq_set_ui(out, 0, 1); return; }
    if (Dp < D0 && Dw < D0) {
        int ip = cl_inner_idx(C, Lp, kp, hp), iw = cl_inner_idx(C, Lw, kw, 1);
        mpq_set_si(out, (long)g_q[ip * N + iw], (unsigned long)H->S);
        mpq_canonicalize(out);
        return;
    }
    mpz_t num, den, t; mpz_inits(num, den, t, NULL);
    if (Dp >= D0) {
        int iw = cl_inner_idx(C, Lw, kw, 1), t1 = cl_t1_idx(C, kp, hp);
        /* ((Dp - D0)*g[INF, iw] + D0*g[t1, iw]) / (S * Dp) */
        mpz_set_si(num, (long)g_q[C->INF_ID * N + iw]); mpz_mul_si(num, num, Dp - D0);
        mpz_set_si(t,   (long)g_q[t1 * N + iw]);        mpz_mul_si(t, t, D0);
        mpz_add(num, num, t);
        mpz_set_si(den, (long)H->S); mpz_mul_si(den, den, Dp);
    } else {
        int ip = cl_inner_idx(C, Lp, kp, hp), t1 = cl_t1_idx(C, kw, 1);
        mpz_set_si(num, (long)g_q[ip * N + C->INF_ID]); mpz_mul_si(num, num, Dw - D0);
        mpz_set_si(t,   (long)g_q[ip * N + t1]);        mpz_mul_si(t, t, D0);
        mpz_add(num, num, t);
        mpz_set_si(den, (long)H->S); mpz_mul_si(den, den, Dw);
    }
    mpz_set(mpq_numref(out), num); mpz_set(mpq_denref(out), den);
    mpq_canonicalize(out);
    mpz_clears(num, den, t, NULL);
}

static void eval_hV_q(mpq_t out, const CertHeader *H, const CertLayout *C, const int64_t *hV_q,
                      int Lp, int kp, int hp, int Lu, int ku, int hu, int Lw, int kw) {
    int ui = cl_hu_idx(C, Lu, ku, hu);
    if (ui < 0) { mpq_set_ui(out, 0, 1); return; }
    int Dp = Lp + kp + hp, Dw = Lw + kw + 1, D0 = C->D0;
    int pid[2], np; long pw[2], pden;
    if (Dp < D0) { pid[0] = cl_inner_idx(C,Lp,kp,hp); pw[0] = 1; pden = 1; np = 1; }
    else { pid[0] = C->INF_ID; pw[0] = Dp - D0; pid[1] = cl_t1_idx(C,kp,hp); pw[1] = D0; pden = Dp; np = 2; }
    int wid[2], nw; long ww[2], wden;
    if (Dw < D0) { wid[0] = cl_inner_idx(C,Lw,kw,1); ww[0] = 1; wden = 1; nw = 1; }
    else { wid[0] = C->INF_ID; ww[0] = Dw - D0; wid[1] = cl_t1_idx(C,kw,1); ww[1] = D0; wden = Dw; nw = 2; }
    mpz_t num, den, t; mpz_inits(num, den, t, NULL);
    mpz_set_ui(num, 0);
    for (int a = 0; a < np; a++)
        for (int b = 0; b < nw; b++) {
            mpz_set_si(t, (long)hV_q[cl_hv_off(C, pid[a], ui, wid[b])]);
            mpz_mul_si(t, t, pw[a] * ww[b]);
            mpz_add(num, num, t);
        }
    mpz_set_si(den, (long)H->S); mpz_mul_si(den, den, pden * wden);
    mpz_set(mpq_numref(out), num); mpz_set(mpq_denref(out), den);
    mpq_canonicalize(out);
    mpz_clears(num, den, t, NULL);
}

static void eval_hC_q(mpq_t out, const CertHeader *H, const CertLayout *C, const int64_t *hC_q,
                      int Lu, int ku, int hu, int L1, int k1, int L2, int k2) {
    int ci = cl_cu_idx(C, Lu, ku, hu);
    if (ci < 0) { mpq_set_ui(out, 0, 1); return; }
    int D1 = L1 + k1 + 1, D2 = L2 + k2 + 1, D0 = C->D0;
    int aid[2], na; long aw[2], aden;
    if (D1 < D0) { aid[0] = cl_inner_idx(C,L1,k1,1); aw[0] = 1; aden = 1; na = 1; }
    else { aid[0] = C->INF_ID; aw[0] = D1 - D0; aid[1] = cl_t1_idx(C,k1,1); aw[1] = D0; aden = D1; na = 2; }
    int bid[2], nb; long bw[2], bden;
    if (D2 < D0) { bid[0] = cl_inner_idx(C,L2,k2,1); bw[0] = 1; bden = 1; nb = 1; }
    else { bid[0] = C->INF_ID; bw[0] = D2 - D0; bid[1] = cl_t1_idx(C,k2,1); bw[1] = D0; bden = D2; nb = 2; }
    mpz_t num, den, t; mpz_inits(num, den, t, NULL);
    mpz_set_ui(num, 0);
    for (int a = 0; a < na; a++)
        for (int b = 0; b < nb; b++) {
            mpz_set_si(t, (long)hC_q[cl_hc_off(C, ci, aid[a], bid[b])]);
            mpz_mul_si(t, t, aw[a] * bw[b]);
            mpz_add(num, num, t);
        }
    mpz_set_si(den, (long)H->S); mpz_mul_si(den, den, aden * bden);
    mpz_set(mpq_numref(out), num); mpz_set(mpq_denref(out), den);
    mpq_canonicalize(out);
    mpz_clears(num, den, t, NULL);
}

/* Exact rational evaluation of Phi(T) for a rooted tree (n, par, deg) */
static void eval_phi_tree_q(mpq_t out, const CertHeader *H, const CertLayout *C,
                            const int64_t *f_q, const int64_t *g_q, const int64_t *hV_q, const int64_t *hC_q,
                            int n, const int *par, const int *deg, int *core_out) {
    int *Lc = calloc(n, sizeof(int)), *Kc = calloc(n, sizeof(int)), *Hc = calloc(n, sizeof(int));
    int *chn = calloc(n, sizeof(int)), **chl = calloc(n, sizeof(int*));
    for (int v = 0; v < n; v++) if (deg[v] >= 2) { Hc[v] = (par[v] != -1) ? 1 : 0; chl[v] = malloc((deg[v]+1)*sizeof(int)); }
    for (int v = 1; v < n; v++) {
        if (deg[v] == 1) Lc[par[v]]++;
        else { Kc[par[v]]++; chl[par[v]][chn[par[v]]++] = v; }
    }
    mpq_t sum, tmp; mpq_inits(sum, tmp, NULL);
    mpq_set_ui(sum, 0, 1);
    int core = 0;
    for (int v = 0; v < n; v++) if (deg[v] >= 2) {
        core++;
        mpq_set_si(tmp, (long)f_q[cl_vidx(C, Lc[v], Kc[v], Hc[v])], (unsigned long)H->S);
        mpq_canonicalize(tmp);
        mpq_add(sum, sum, tmp);
    }
    for (int w = 1; w < n; w++) if (deg[w] >= 2) {
        int p = par[w];
        eval_g_q(tmp, H, C, g_q, Lc[p], Kc[p], Hc[p], Lc[w], Kc[w]);
        mpq_add(sum, sum, tmp);
        if (par[p] != -1) {
            int gp = par[p];
            eval_hV_q(tmp, H, C, hV_q, Lc[gp], Kc[gp], Hc[gp], Lc[p], Kc[p], 1, Lc[w], Kc[w]);
            mpq_add(sum, sum, tmp);
        }
    }
    for (int u = 0; u < n; u++) if (deg[u] >= 2 && chn[u] >= 2)
        for (int i = 0; i < chn[u]; i++)
            for (int j = i + 1; j < chn[u]; j++) {
                eval_hC_q(tmp, H, C, hC_q, Lc[u], Kc[u], Hc[u], Lc[chl[u][i]], Kc[chl[u][i]], Lc[chl[u][j]], Kc[chl[u][j]]);
                mpq_add(sum, sum, tmp);
            }
    mpq_set(out, sum);
    if (core_out) *core_out = core;
    mpq_clears(sum, tmp, NULL);
    for (int v = 0; v < n; v++) free(chl[v]);
    free(Lc); free(Kc); free(Hc); free(chn); free(chl);
}

/* Brute-force exact rational drift check on a single tree T:
 * Computes exact D(T) = sum_{V=0..n-1} (1/deg[V]) * sum_{W~V} (Phi(T + leaf@W) - Phi(T))
 * and checks: UPPER: D(T) <= (P/Q)*Phi(T) and Phi(T) >= c0*C(T)
 *             LOWER: D(T) >= (P/Q)*Phi(T) and Phi(T) <= C0*C(T)
 */
static int check_tree_exact(const CertHeader *H, const CertLayout *C,
                            const int64_t *f_q, const int64_t *g_q, const int64_t *hV_q, const int64_t *hC_q,
                            const mpq_t bound_c, int n, const int *par, const int *deg, const char *label) {
    mpq_t phi0, phi2, drift, diff, term, rhs, ratio, c_cmp;
    mpq_inits(phi0, phi2, drift, diff, term, rhs, ratio, c_cmp, NULL);
    int core0 = 0;
    eval_phi_tree_q(phi0, H, C, f_q, g_q, hV_q, hC_q, n, par, deg, &core0);

    /* Build adjacency */
    int *cnt = calloc(n, sizeof(int)), **adj = calloc(n, sizeof(int*));
    for (int v = 0; v < n; v++) adj[v] = malloc((deg[v] + 1) * sizeof(int));
    for (int v = 1; v < n; v++) { adj[v][cnt[v]++] = par[v]; adj[par[v]][cnt[par[v]]++] = v; }

    /* Group by attachment target W to speed up: weight(W) = sum_{V ~ W} 1/deg[V] */
    mpq_set_ui(drift, 0, 1);
    int *p2 = malloc((n + 1) * sizeof(int)), *d2 = malloc((n + 1) * sizeof(int));
    memcpy(p2, par, n * sizeof(int)); memcpy(d2, deg, n * sizeof(int));
    for (int W = 0; W < n; W++) {
        mpq_t wtW, invd; mpq_inits(wtW, invd, NULL);
        mpq_set_ui(wtW, 0, 1);
        for (int i = 0; i < cnt[W]; i++) {
            mpq_set_ui(invd, 1, (unsigned long)deg[adj[W][i]]);
            mpq_add(wtW, wtW, invd);
        }
        p2[n] = W; d2[n] = 1; d2[W]++;
        eval_phi_tree_q(phi2, H, C, f_q, g_q, hV_q, hC_q, n + 1, p2, d2, NULL);
        d2[W]--;
        mpq_sub(diff, phi2, phi0);
        mpq_mul(term, wtW, diff);
        mpq_add(drift, drift, term);
        mpq_clears(wtW, invd, NULL);
    }
    /* Compare drift vs (P/Q)*phi0 */
    mpq_set_si(ratio, H->P, (unsigned long)H->Q); mpq_canonicalize(ratio);
    mpq_mul(rhs, ratio, phi0);
    int ok_drift = (H->mode == 1) ? (mpq_cmp(drift, rhs) <= 0) : (mpq_cmp(drift, rhs) >= 0);

    /* Compare phi0 vs bound_c * core0 */
    mpq_set_si(c_cmp, core0, 1); mpq_mul(c_cmp, c_cmp, bound_c);
    int ok_bound = (H->mode == 1) ? (mpq_cmp(phi0, c_cmp) >= 0) : (mpq_cmp(phi0, c_cmp) <= 0);

    if (!ok_drift || !ok_bound) {
        fprintf(stderr, "FAIL exact tree check [%s] (n=%d, C=%d): ok_drift=%d ok_bound=%d drift=%.10f rhs=%.10f phi=%.10f bound=%.10f\n",
                label, n, core0, ok_drift, ok_bound, mpq_get_d(drift), mpq_get_d(rhs), mpq_get_d(phi0), mpq_get_d(c_cmp));
    }
    for (int v = 0; v < n; v++) free(adj[v]);
    free(adj); free(cnt); free(p2); free(d2);
    mpq_clears(phi0, phi2, drift, diff, term, rhs, ratio, c_cmp, NULL);
    return ok_drift && ok_bound;
}

/* Helper: exact GMP check of sign of LHS - (P/Q)*v.
 * In all inequalities, we clear denominators cleanly into an exact mpz_t LHS and exact mpz_t RHS!
 * For UPPER: pass iff LHS <= RHS (diff = LHS - RHS <= 0).
 * For LOWER: pass iff LHS >= RHS (diff = LHS - RHS >= 0).
 */
static inline int check_sign(int mode, const mpz_t lhs, const mpz_t rhs) {
    int cmp = mpz_cmp(lhs, rhs);
    return (mode == 1) ? (cmp <= 0) : (cmp >= 0);
}

/* Helper: exact rational extreme choice for rate r in {r_lo, r_hi} multiplying (target - v).
 * Returns exact mpz_t num contribution where common denominator is DEN (so rate * DEN is exact integer!).
 */
static inline long choose_rate(int mode, int64_t target, int64_t v, long r_lo_scaled, long r_hi_scaled) {
    if (mode == 1) return (target > v) ? r_hi_scaled : r_lo_scaled;
    else           return (target < v) ? r_hi_scaled : r_lo_scaled;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <cert.bin> [--override-pq P Q] [--perturb f|g|hV|hC] [--skip-trees]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    int ov_p = -1, ov_q = -1, skip_trees = 0;
    const char *perturb = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--override-pq") && i + 2 < argc) { ov_p = atoi(argv[++i]); ov_q = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--perturb") && i + 1 < argc) { perturb = argv[++i]; }
        else if (!strcmp(argv[i], "--skip-trees")) { skip_trees = 1; }
    }

    FILE *fp = fopen(path, "rb");
    if (!fp) { perror("fopen"); return 1; }
    CertHeader H;
    if (fread(&H, sizeof(H), 1, fp) != 1) { fprintf(stderr, "Short header\n"); return 1; }
    if (H.magic != CERT_MAGIC) { fprintf(stderr, "Bad magic\n"); return 1; }
    if (ov_p > 0 && ov_q > 0) { H.P = ov_p; H.Q = ov_q; }

    /* Structural header checks */
    if (H.mode != 1 && H.mode != 2) { fprintf(stderr, "Bad mode\n"); return 1; }
    if (H.P <= 0 || H.Q <= 0 || H.S <= 0 || H.S_r <= 0 || H.Q_K <= 0) { fprintf(stderr, "Bad scales\n"); return 1; }
    if (H.DHU > H.D0 || H.DCU > H.D0 || H.KT + 1 > H.KF || H.D0 > H.LF) { fprintf(stderr, "Bad domain nesting\n"); return 1; }
    if (H.nknots != CERT_NKNOTS || H.knot_num[0] != 0 || (long)H.knot_num[1] * H.D0 != H.Q_K || H.knot_num[H.nknots - 1] * 2 != H.Q_K) {
        fprintf(stderr, "Bad knots\n"); return 1;
    }
    for (int i = 0; i < H.nknots - 1; i++) if (H.knot_num[i] >= H.knot_num[i + 1]) { fprintf(stderr, "Non-increasing knots\n"); return 1; }

    CertLayout C;
    cl_init(&C, &H);
    int N = C.NG, D0 = C.D0, KT = C.KT;
    size_t hv_len = (size_t)C.npars * C.NHU_ACT * C.nchs;
    size_t hc_len = (size_t)C.NCU_ACT * C.nchs * (C.nchs + 1) / 2;

    int64_t *rpar_q = malloc(C.INNER_N * sizeof(int64_t));
    int64_t *rch_q  = malloc(C.INNER_N * sizeof(int64_t));
    int64_t *xpar_q = malloc(C.INNER_N * sizeof(int64_t));
    int64_t *xch_q  = malloc(C.INNER_N * sizeof(int64_t));
    int64_t *xcc_q  = malloc(C.INNER_N * sizeof(int64_t));
    int64_t (*k_pc_q)[CERT_NKNOTS] = malloc((size_t)C.INNER_N * sizeof(*k_pc_q));
    int64_t (*k_cp_q)[CERT_NKNOTS] = malloc((size_t)C.INNER_N * sizeof(*k_cp_q));
    int64_t (*k_cc_q)[CERT_NKNOTS] = malloc((size_t)C.INNER_N * sizeof(*k_cc_q));
    int64_t *f_q  = malloc(C.NV * sizeof(int64_t));
    int64_t *g_q  = malloc((size_t)N * N * sizeof(int64_t));
    int64_t *hV_q = malloc(hv_len * sizeof(int64_t));
    int64_t *hC_q = malloc(hc_len * sizeof(int64_t));

    size_t ok = 1;
    ok &= (fread(rpar_q, sizeof(int64_t), C.INNER_N, fp) == (size_t)C.INNER_N);
    ok &= (fread(rch_q,  sizeof(int64_t), C.INNER_N, fp) == (size_t)C.INNER_N);
    ok &= (fread(xpar_q, sizeof(int64_t), C.INNER_N, fp) == (size_t)C.INNER_N);
    ok &= (fread(xch_q,  sizeof(int64_t), C.INNER_N, fp) == (size_t)C.INNER_N);
    ok &= (fread(xcc_q,  sizeof(int64_t), C.INNER_N, fp) == (size_t)C.INNER_N);
    ok &= (fread(k_pc_q, sizeof(int64_t), (size_t)C.INNER_N * CERT_NKNOTS, fp) == (size_t)C.INNER_N * CERT_NKNOTS);
    ok &= (fread(k_cp_q, sizeof(int64_t), (size_t)C.INNER_N * CERT_NKNOTS, fp) == (size_t)C.INNER_N * CERT_NKNOTS);
    ok &= (fread(k_cc_q, sizeof(int64_t), (size_t)C.INNER_N * CERT_NKNOTS, fp) == (size_t)C.INNER_N * CERT_NKNOTS);
    ok &= (fread(f_q,    sizeof(int64_t), C.NV, fp) == (size_t)C.NV);
    ok &= (fread(g_q,    sizeof(int64_t), (size_t)N * N, fp) == (size_t)N * N);
    ok &= (fread(hV_q,   sizeof(int64_t), hv_len, fp) == hv_len);
    ok &= (fread(hC_q,   sizeof(int64_t), hc_len, fp) == hc_len);
    fclose(fp);
    if (!ok) { fprintf(stderr, "Truncated cert file\n"); return 1; }

    /* Optional perturbation for negative control */
    if (perturb) {
        int64_t delta = (H.mode == 1) ? -(H.S / 1000) : (H.S / 1000);
        if (!strcmp(perturb, "f"))  f_q[cl_vidx(&C, 1, 0, 1)] += delta;
        if (!strcmp(perturb, "g"))  g_q[cl_inner_idx(&C,1,1,0)*N + cl_inner_idx(&C,1,0,1)] += delta;
        if (!strcmp(perturb, "hV")) hV_q[0] += delta;
        if (!strcmp(perturb, "hC")) hC_q[0] += delta;
    }

    printf("=== EXACT GMP VERIFIER ===\n");
    printf("File: %s | Mode: %s | Target: %d/%d = %.8f | S=%lld S_r=%lld\n",
           path, H.mode==1?"UPPER":"LOWER", H.P, H.Q, (double)H.P/H.Q, (long long)H.S, (long long)H.S_r);
    printf("Domains: f(%d,%d) D0=%d(KT=%d,NG=%d) hV(%d,%d;act=%d) hC(%d,%d;act=%d)\n",
           C.LF, C.KF, D0, KT, N, C.DHU, C.KHU, C.NHU_ACT, C.DCU, C.KCU, C.NCU_ACT);
    fflush(stdout);

    /* Check g(tail, tail) == 0 */
    for (int p = C.INNER_N; p < N; p++)
        for (int w = C.INNER_N; w < N; w++)
            if (g_q[p * N + w] != 0) { fprintf(stderr, "FAIL: g(tail,tail) != 0 at (%d,%d)\n", p, w); return 1; }

    /* Precompute exact rational A_in[p, w] and A_out[p, w] as mpq_t?
     * Even cleaner and faster: note every x = gInvD[v] is 0, 1/D (2<=D<D0), or 1/D0.
     * Precompute exact mpq_t A_tot[p*N + w] for the ~160k active (p, w) pairs!
     */
    mpq_t *A_tot = malloc((size_t)N * N * sizeof(mpq_t));
    for (size_t i = 0; i < (size_t)N * N; i++) { mpq_init(A_tot[i]); mpq_set_ui(A_tot[i], 0, 1); }
    {
        mpq_t tv, mv; mpq_inits(tv, mv, NULL);
        for (int u = 0; u < C.INNER_N; u++) {
            int Lu = C.gL[u], ku = C.gK[u], hu = C.gH[u];
            int in_hu = (cl_hu_idx(&C, Lu, ku, hu) >= 0);
            int in_cu = (cl_cu_idx(&C, Lu, ku, hu) >= 0);
            if (!in_hu && hu == 1 && ku >= 1) {
                for (int ip = 0; ip < C.npars; ip++) {
                    int p = C.pars[ip];
                    long xn = (C.gKind[p] == 1) ? 0 : 1, xd = (C.gKind[p] == 1) ? 1 : C.gD[p];
                    eval_pl_q(tv, &H, k_pc_q[u], xn, xd);
                    mpq_set_si(mv, ku, 1); mpq_mul(tv, tv, mv);
                    mpq_add(A_tot[p * N + u], A_tot[p * N + u], tv);
                }
                for (int iw = 0; iw < C.nchs; iw++) {
                    int w = C.chs[iw];
                    long xn = (C.gKind[w] == 1) ? 0 : 1, xd = (C.gKind[w] == 1) ? 1 : C.gD[w];
                    eval_pl_q(tv, &H, k_cp_q[u], xn, xd);
                    mpq_add(A_tot[u * N + w], A_tot[u * N + w], tv);
                }
            }
            if (!in_cu && ku >= 2) {
                for (int iw = 0; iw < C.nchs; iw++) {
                    int w = C.chs[iw];
                    long xn = (C.gKind[w] == 1) ? 0 : 1, xd = (C.gKind[w] == 1) ? 1 : C.gD[w];
                    eval_pl_q(tv, &H, k_cc_q[u], xn, xd);
                    mpq_set_si(mv, ku - 1, 1); mpq_mul(tv, tv, mv);
                    mpq_add(A_tot[u * N + w], A_tot[u * N + w], tv);
                }
            }
        }
        mpq_clears(tv, mv, NULL);
    }

    int v_nb = cl_vidx(&C, 1, 0, 1), iw_nb = cl_inner_idx(&C, 1, 0, 1);

    /* ===== CHECK 1: ALL VERTEX INEQUALITIES f (EXACT mpq_t) ===== */
    long long f_checks = 0, f_fails = 0;
    {
        mpq_t lhs, rhs, tmp, gnb, gnb_t0, gnb_t1, bq, rq, lamq, vq, fu, fup, fub, fnb;
        mpq_inits(lhs, rhs, tmp, gnb, gnb_t0, gnb_t1, bq, rq, lamq, vq, fu, fup, fub, fnb, NULL);
        mpq_set_si(lamq, H.P, (unsigned long)H.Q); mpq_canonicalize(lamq);
        mpq_set_si(fnb, (long)f_q[v_nb], (unsigned long)H.S); mpq_canonicalize(fnb);

        /* VSTAR0 and VSTAR1 (k > KF): D >= KF + 1 + h, b in [0, 1], t in [0, D0/(KF+1+h)] */
        for (int hb = 0; hb <= 1; hb++) {
            int vs = hb ? C.VSTAR1 : C.VSTAR0;
            mpq_set_si(fu, (long)f_q[vs], (unsigned long)H.S); mpq_canonicalize(fu);
            mpq_mul(rhs, lamq, fu);
            /* g(u' -> nb) at t=0 (INF) and t=D0/(KF+1+h) */
            mpq_set_si(gnb_t0, (long)g_q[C.INF_ID * N + iw_nb], (unsigned long)H.S); mpq_canonicalize(gnb_t0);
            eval_g_q(gnb_t1, &H, &C, g_q, 0, C.KF + 1, hb, 1, 0); /* D = KF+1+h, kap = KT+1 */
            mpq_t *g_ends[2] = { &gnb_t0, &gnb_t1 };
            for (int bi = 0; bi <= 1; bi++) {
                for (int gi = 0; gi < 2; gi++) {
                    if (bi == 0) mpq_set_ui(lhs, 0, 1);
                    else { mpq_add(lhs, fnb, *g_ends[gi]); }
                    f_checks++;
                    int ok_c = (H.mode == 1) ? (mpq_cmp(lhs, rhs) <= 0) : (mpq_cmp(lhs, rhs) >= 0);
                    if (!ok_c) {
                        f_fails++;
                        if (f_fails <= 5) fprintf(stderr, "FAIL f[VSTAR%d] bi=%d gi=%d lhs=%.10f rhs=%.10f\n", hb, bi, gi, mpq_get_d(lhs), mpq_get_d(rhs));
                    }
                }
            }
        }

        /* All (L, k, h) in [0..LF] x [0..KF] x {0, 1} with D >= 2 */
        for (int hb = 0; hb <= 1; hb++) {
            for (int k = 0; k <= C.KF; k++) {
                for (int L = 0; L <= C.LF; L++) {
                    int D = L + k + hb;
                    if (D < 2) continue;
                    int u = cl_vidx(&C, L, k, hb);
                    mpq_set_si(fu, (long)f_q[u], (unsigned long)H.S); mpq_canonicalize(fu);
                    mpq_mul(rhs, lamq, fu);

                    if (L == C.LF) {
                        /* (a) Exact L == LF */
                        int ub = cl_vidx(&C, L - 1, k + 1, hb);
                        mpq_set_si(fub, (long)f_q[ub], (unsigned long)H.S); mpq_canonicalize(fub);
                        eval_g_q(gnb, &H, &C, g_q, L - 1, k + 1, hb, 1, 0);
                        mpq_sub(tmp, fub, fu); mpq_add(tmp, tmp, fnb); mpq_add(tmp, tmp, gnb);
                        mpq_set_si(bq, L, (unsigned long)D); mpq_canonicalize(bq);
                        mpq_mul(lhs, bq, tmp);
                        f_checks++;
                        int ok_c = (H.mode == 1) ? (mpq_cmp(lhs, rhs) <= 0) : (mpq_cmp(lhs, rhs) >= 0);
                        if (!ok_c) { f_fails++; if (f_fails<=5) fprintf(stderr, "FAIL f(LF=%d,k=%d,h=%d) exact\n", L, k, hb); }

                        /* (b) Tail L >= LF + 1: b in [(LF+1)/(LF+1+k+h), 1], t in [0, D0/(LF+1+k+h)] */
                        int ubt = cl_vidx(&C, C.LF, k + 1, hb);
                        mpq_set_si(fub, (long)f_q[ubt], (unsigned long)H.S); mpq_canonicalize(fub);
                        mpq_set_si(gnb_t0, (long)g_q[C.INF_ID * N + iw_nb], (unsigned long)H.S); mpq_canonicalize(gnb_t0);
                        eval_g_q(gnb_t1, &H, &C, g_q, C.LF, k + 1, hb, 1, 0); /* D = LF+1+k+h */
                        mpq_t *g_ends[2] = { &gnb_t0, &gnb_t1 };
                        for (int bi = 0; bi <= 1; bi++) {
                            if (bi == 0) { mpq_set_si(bq, C.LF + 1, (unsigned long)(C.LF + 1 + k + hb)); mpq_canonicalize(bq); }
                            else mpq_set_ui(bq, 1, 1);
                            for (int gi = 0; gi < 2; gi++) {
                                mpq_sub(tmp, fub, fu); mpq_add(tmp, tmp, fnb); mpq_add(tmp, tmp, *g_ends[gi]);
                                mpq_mul(lhs, bq, tmp);
                                f_checks++;
                                int ok2 = (H.mode == 1) ? (mpq_cmp(lhs, rhs) <= 0) : (mpq_cmp(lhs, rhs) >= 0);
                                if (!ok2) { f_fails++; if (f_fails<=5) fprintf(stderr, "FAIL f(tail LF,k=%d,h=%d) bi=%d gi=%d\n", k, hb, bi, gi); }
                            }
                        }
                    } else {
                        int u_up = cl_vidx(&C, L + 1, k, hb);
                        mpq_set_si(fup, (long)f_q[u_up], (unsigned long)H.S); mpq_canonicalize(fup);
                        mpq_t bud_term; mpq_init(bud_term);
                        if (L >= 1) {
                            int ub = cl_vidx(&C, L - 1, k + 1, hb);
                            mpq_set_si(fub, (long)f_q[ub], (unsigned long)H.S); mpq_canonicalize(fub);
                            eval_g_q(gnb, &H, &C, g_q, L - 1, k + 1, hb, 1, 0);
                            mpq_sub(tmp, fub, fu); mpq_add(tmp, tmp, fnb); mpq_add(tmp, tmp, gnb);
                            mpq_set_si(bq, L, (unsigned long)D); mpq_canonicalize(bq);
                            mpq_mul(bud_term, bq, tmp);
                        } else mpq_set_ui(bud_term, 0, 1);

                        if (D < D0) {
                            int gu = cl_inner_idx(&C, L, k, hb);
                            /* rate = L + hb*rpar[gu]/S_r + k*rch[gu]/S_r */
                            long rnum = (long)L * H.S_r + (long)hb * rpar_q[gu] + (long)k * rch_q[gu];
                            mpq_set_si(rq, rnum, (unsigned long)H.S_r); mpq_canonicalize(rq);
                            mpq_sub(tmp, fup, fu); mpq_mul(lhs, rq, tmp); mpq_add(lhs, lhs, bud_term);
                            f_checks++;
                            int ok_c = (H.mode == 1) ? (mpq_cmp(lhs, rhs) <= 0) : (mpq_cmp(lhs, rhs) >= 0);
                            if (!ok_c) { f_fails++; if (f_fails<=5) fprintf(stderr, "FAIL f_in(L=%d,k=%d,h=%d) lhs=%.10f rhs=%.10f\n", L, k, hb, mpq_get_d(lhs), mpq_get_d(rhs)); }
                        } else {
                            /* rate in {L, L + (hb+k)/2} */
                            for (int ri = 0; ri <= 1; ri++) {
                                if (ri == 0) mpq_set_si(rq, L, 1);
                                else { mpq_set_si(rq, 2 * L + hb + k, 2); mpq_canonicalize(rq); }
                                mpq_sub(tmp, fup, fu); mpq_mul(lhs, rq, tmp); mpq_add(lhs, lhs, bud_term);
                                f_checks++;
                                int ok_c = (H.mode == 1) ? (mpq_cmp(lhs, rhs) <= 0) : (mpq_cmp(lhs, rhs) >= 0);
                                if (!ok_c) { f_fails++; if (f_fails<=5) fprintf(stderr, "FAIL f_out(L=%d,k=%d,h=%d) ri=%d lhs=%.10f rhs=%.10f\n", L, k, hb, ri, mpq_get_d(lhs), mpq_get_d(rhs)); }
                            }
                        }
                        mpq_clear(bud_term);
                    }
                }
            }
        }
        mpq_clears(lhs, rhs, tmp, gnb, gnb_t0, gnb_t1, bq, rq, lamq, vq, fu, fup, fub, fnb, NULL);
    }
    printf("  [1/5] f checks:  %10lld  fails: %lld  (%s)\n", f_checks, f_fails, f_fails==0?"PASS":"FAIL");
    fflush(stdout);

    /* Precompute exact integer dLf_g_q[u] = f_q(L+1,k,h) - f_q(L,k,h) (scale S) */
    int64_t *dLf_q = calloc(N, sizeof(int64_t));
    for (int u = 0; u < C.INNER_N; u++)
        dLf_q[u] = f_q[cl_vidx(&C, C.gL[u] + 1, C.gK[u], C.gH[u])] - f_q[cl_vidx(&C, C.gL[u], C.gK[u], C.gH[u])];

    /* ===== CHECK 2: ALL DIRECTED EDGE INEQUALITIES g(p -> w) (EXACT mpq_t) ===== */
    long long g_checks = 0, g_fails = 0;
    #pragma omp parallel num_threads(3) reduction(+:g_checks,g_fails)
    {
        mpq_t lhs, rhs, tmp, q1, q2, lamq, gval, d0, d1, best, cand;
        mpq_inits(lhs, rhs, tmp, q1, q2, lamq, gval, d0, d1, best, cand, NULL);
        mpq_set_si(lamq, H.P, (unsigned long)H.Q); mpq_canonicalize(lamq);

        #pragma omp for schedule(dynamic)
        for (int ip = 0; ip < C.npars; ip++) {
            int p = C.pars[ip], pk = C.gKind[p], Lp = C.gL[p], kp = C.gK[p], hp = C.gH[p];
            int p_up  = (pk == 0) ? ((C.gD[p]+1 < D0) ? cl_inner_idx(&C,Lp+1,kp,hp) : cl_t1_idx(&C,kp,hp)) : C.INF_ID;
            int p_bud = (pk == 0 && Lp>=1) ? cl_inner_idx(&C,Lp-1,kp+1,hp) : ((pk == 2 && kp <= KT) ? cl_t1_idx(&C,kp+1,hp) : p);
            int p_in_hu = (pk == 0 && cl_hu_idx(&C,Lp,kp,hp) >= 0);
            int p_in_cu = (pk == 0 && cl_cu_idx(&C,Lp,kp,hp) >= 0);
            int has_pc = (pk == 0 && hp == 1 && !p_in_hu);
            int has_cc = (pk == 0 && kp >= 2 && !p_in_cu);
            int cu_pb = (pk == 0 && Lp >= 1) ? cl_cu_idx(&C, Lp-1, kp+1, hp) : -1;

            for (int iw = 0; iw < C.nchs; iw++) {
                int w = C.chs[iw], wk = C.gKind[w], Lw = C.gL[w], kw = C.gK[w];
                if (pk != 0 && wk != 0) continue; /* (Tail, Tail) == 0 already verified */
                int w_up  = (wk == 0) ? ((C.gD[w]+1 < D0) ? cl_inner_idx(&C,Lw+1,kw,1) : cl_t1_idx(&C,kw,1)) : C.INF_ID;
                int w_bud = (wk == 0 && Lw>=1) ? cl_inner_idx(&C,Lw-1,kw+1,1) : ((wk == 2 && kw <= KT) ? cl_t1_idx(&C,kw+1,1) : w);
                int w_in_hu = (wk == 0 && cl_hu_idx(&C,Lw,kw,1) >= 0);
                int has_cp = (wk == 0 && kw >= 1 && !w_in_hu);
                int hu_wb = (wk == 0 && Lw >= 1) ? cl_hu_idx(&C, Lw-1, kw+1, 1) : -1;

                int64_t gv = g_q[p * N + w];
                int64_t gp_plus = g_q[p_up * N + w], gw_plus = g_q[p * N + w_up];
                int64_t gp_bud  = g_q[p_bud * N + w], gw_bud = g_q[p * N + w_bud];
                int64_t dgp = gp_plus - gv, dgw = gw_plus - gv;

                mpq_set(lhs, A_tot[p * N + w]);

                /* p contributions */
                if (pk == 0) {
                    /* (inv_w - rch[p]) * dLf_g[p] */
                    if (wk == 1) mpq_set_si(q1, -rch_q[p], (unsigned long)H.S_r);
                    else { mpq_set_si(q1, (long)H.S_r - (long)C.gD[w] * rch_q[p], (unsigned long)((long)C.gD[w] * H.S_r)); }
                    mpq_canonicalize(q1);
                    mpq_set_si(q2, (long)dLf_q[p], (unsigned long)H.S); mpq_canonicalize(q2);
                    mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);

                    /* ur_p * dgp: ur_p = Lp + (wk==0 ? 1/Dw : 0) + (p_in_hu ? xpar[p] : 0) + (p_in_cu ? (kp-1)*xcc[p] : 0) */
                    long rnum = (long)Lp * H.S_r + (p_in_hu ? xpar_q[p] : 0) + (p_in_cu ? (long)(kp - 1) * xcc_q[p] : 0);
                    mpq_set_si(q1, rnum, (unsigned long)H.S_r); mpq_canonicalize(q1);
                    if (wk == 0) { mpq_set_ui(q2, 1, (unsigned long)C.gD[w]); mpq_add(q1, q1, q2); }
                    mpq_set_si(q2, (long)dgp, (unsigned long)H.S); mpq_canonicalize(q2);
                    mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);

                    /* p budding */
                    if (Lp >= 1) {
                        int64_t hc_nb = (cu_pb >= 0) ? hC_q[cl_hc_off(&C, cu_pb, w, iw_nb)] : 0;
                        mpq_set_si(q1, Lp, (unsigned long)C.gD[p]); mpq_canonicalize(q1);
                        mpq_set_si(q2, (long)(gp_bud - gv + hc_nb), (unsigned long)H.S); mpq_canonicalize(q2);
                        mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);
                    }
                } else if (pk == 2) {
                    /* Tail T1(kp, hp) leaf gain: gamma_p in {gamma_min, 1} * (g[INF, w] - g[T1, w]) */
                    int m = (kp <= KT && D0 - kp - hp > 0) ? (D0 - kp - hp) : 0;
                    if ((H.mode == 1 && dgp > 0) || (H.mode == 2 && dgp < 0)) mpq_set_ui(q1, 1, 1);
                    else { mpq_set_si(q1, m, (unsigned long)(D0 + 1)); mpq_canonicalize(q1); }
                    mpq_set_si(q2, (long)dgp, (unsigned long)H.S); mpq_canonicalize(q2);
                    mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);

                    /* Tail T1(kp, hp) budding (if kp <= KT): b_p in {b_min, 1} * (g[T1(kp+1), w] - g[T1, w]) */
                    if (kp <= KT) {
                        int64_t db = gp_bud - gv;
                        if ((H.mode == 1 && db > 0) || (H.mode == 2 && db < 0)) mpq_set_ui(q1, 1, 1);
                        else { mpq_set_si(q1, m, (unsigned long)D0); mpq_canonicalize(q1); }
                        mpq_set_si(q2, (long)db, (unsigned long)H.S); mpq_canonicalize(q2);
                        mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);
                    }
                    /* Quadratic term (1/D0) * ext(g[INF,w+]-g[INF,w], g[T1,w+]-g[T1,w]) */
                    int64_t d0_w = g_q[C.INF_ID * N + w_up] - g_q[C.INF_ID * N + w];
                    int64_t dext = (H.mode == 1) ? (dgw > d0_w ? dgw : d0_w) : (dgw < d0_w ? dgw : d0_w);
                    mpq_set_si(q1, (long)dext, (unsigned long)((long)D0 * H.S)); mpq_canonicalize(q1);
                    mpq_add(lhs, lhs, q1);
                }

                /* w contributions */
                if (wk == 0) {
                    if (pk == 1) mpq_set_si(q1, -rpar_q[w], (unsigned long)H.S_r);
                    else { mpq_set_si(q1, (long)H.S_r - (long)C.gD[p] * rpar_q[w], (unsigned long)((long)C.gD[p] * H.S_r)); }
                    mpq_canonicalize(q1);
                    mpq_set_si(q2, (long)dLf_q[w], (unsigned long)H.S); mpq_canonicalize(q2);
                    mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);

                    long rnum = (long)Lw * H.S_r + (w_in_hu ? (long)kw * xch_q[w] : 0);
                    mpq_set_si(q1, rnum, (unsigned long)H.S_r); mpq_canonicalize(q1);
                    if (pk == 0) { mpq_set_ui(q2, 1, (unsigned long)C.gD[p]); mpq_add(q1, q1, q2); }
                    mpq_set_si(q2, (long)dgw, (unsigned long)H.S); mpq_canonicalize(q2);
                    mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);

                    if (Lw >= 1) {
                        int64_t hv_nb = (hu_wb >= 0) ? hV_q[cl_hv_off(&C, p, hu_wb, iw_nb)] : 0;
                        mpq_set_si(q1, Lw, (unsigned long)C.gD[w]); mpq_canonicalize(q1);
                        mpq_set_si(q2, (long)(gw_bud - gv + hv_nb), (unsigned long)H.S); mpq_canonicalize(q2);
                        mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);
                    }
                } else if (wk == 2) {
                    int m = (kw <= KT && D0 - kw - 1 > 0) ? (D0 - kw - 1) : 0;
                    if ((H.mode == 1 && dgw > 0) || (H.mode == 2 && dgw < 0)) mpq_set_ui(q1, 1, 1);
                    else { mpq_set_si(q1, m, (unsigned long)(D0 + 1)); mpq_canonicalize(q1); }
                    mpq_set_si(q2, (long)dgw, (unsigned long)H.S); mpq_canonicalize(q2);
                    mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);

                    if (kw <= KT) {
                        int64_t db = gw_bud - gv;
                        if ((H.mode == 1 && db > 0) || (H.mode == 2 && db < 0)) mpq_set_ui(q1, 1, 1);
                        else { mpq_set_si(q1, m, (unsigned long)D0); mpq_canonicalize(q1); }
                        mpq_set_si(q2, (long)db, (unsigned long)H.S); mpq_canonicalize(q2);
                        mpq_mul(tmp, q1, q2); mpq_add(lhs, lhs, tmp);
                    }
                    int64_t d0_p = g_q[p_up * N + C.INF_ID] - g_q[p * N + C.INF_ID];
                    int64_t dext = (H.mode == 1) ? (dgp > d0_p ? dgp : d0_p) : (dgp < d0_p ? dgp : d0_p);
                    mpq_set_si(q1, (long)dext, (unsigned long)((long)D0 * H.S)); mpq_canonicalize(q1);
                    mpq_add(lhs, lhs, q1);
                }

                /* Fenchel-Young conjugate A^*(d) = ext_k ((knot_num[k]/Q_K)*d - k_q[k]) */
                if (has_pc) {
                    for (int kk = 0; kk < CERT_NKNOTS; kk++) {
                        /* ((knot_num[kk]*dgp) - Q_K*k_pc_q[p][kk]) / (Q_K * S) */
                        mpz_t zn; mpz_init(zn);
                        mpz_set_si(zn, (long)dgp); mpz_mul_si(zn, zn, H.knot_num[kk]);
                        mpz_t zt; mpz_init(zt);
                        mpz_set_si(zt, (long)k_pc_q[p][kk]); mpz_mul_si(zt, zt, H.Q_K);
                        mpz_sub(zn, zn, zt);
                        mpz_set(mpq_numref(cand), zn);
                        mpz_set_si(mpq_denref(cand), (long)H.Q_K * H.S);
                        mpq_canonicalize(cand);
                        mpz_clears(zn, zt, NULL);
                        if (kk == 0 || (H.mode == 1 ? mpq_cmp(cand, best) > 0 : mpq_cmp(cand, best) < 0)) mpq_set(best, cand);
                    }
                    mpq_add(lhs, lhs, best);
                }
                if (has_cc) {
                    for (int kk = 0; kk < CERT_NKNOTS; kk++) {
                        mpz_t zn, zt; mpz_inits(zn, zt, NULL);
                        mpz_set_si(zn, (long)dgp); mpz_mul_si(zn, zn, H.knot_num[kk]);
                        mpz_set_si(zt, (long)k_cc_q[p][kk]); mpz_mul_si(zt, zt, H.Q_K);
                        mpz_sub(zn, zn, zt);
                        mpz_set(mpq_numref(cand), zn);
                        mpz_set_si(mpq_denref(cand), (long)H.Q_K * H.S);
                        mpq_canonicalize(cand);
                        mpz_clears(zn, zt, NULL);
                        if (kk == 0 || (H.mode == 1 ? mpq_cmp(cand, best) > 0 : mpq_cmp(cand, best) < 0)) mpq_set(best, cand);
                    }
                    mpq_set_si(q1, kp - 1, 1); mpq_mul(best, best, q1);
                    mpq_add(lhs, lhs, best);
                }
                if (has_cp) {
                    for (int kk = 0; kk < CERT_NKNOTS; kk++) {
                        mpz_t zn, zt; mpz_inits(zn, zt, NULL);
                        mpz_set_si(zn, (long)dgw); mpz_mul_si(zn, zn, H.knot_num[kk]);
                        mpz_set_si(zt, (long)k_cp_q[w][kk]); mpz_mul_si(zt, zt, H.Q_K);
                        mpz_sub(zn, zn, zt);
                        mpz_set(mpq_numref(cand), zn);
                        mpz_set_si(mpq_denref(cand), (long)H.Q_K * H.S);
                        mpq_canonicalize(cand);
                        mpz_clears(zn, zt, NULL);
                        if (kk == 0 || (H.mode == 1 ? mpq_cmp(cand, best) > 0 : mpq_cmp(cand, best) < 0)) mpq_set(best, cand);
                    }
                    mpq_set_si(q1, kw, 1); mpq_mul(best, best, q1);
                    mpq_add(lhs, lhs, best);
                }

                mpq_set_si(gval, (long)gv, (unsigned long)H.S); mpq_canonicalize(gval);
                mpq_mul(rhs, lamq, gval);
                g_checks++;
                int ok_c = (H.mode == 1) ? (mpq_cmp(lhs, rhs) <= 0) : (mpq_cmp(lhs, rhs) >= 0);
                if (!ok_c) {
                    g_fails++;
                    #pragma omp critical
                    if (g_fails <= 5) fprintf(stderr, "FAIL g(p=%d[k=%d,D=%d] -> w=%d[k=%d,D=%d]) lhs=%.10f rhs=%.10f\n",
                                              p, pk, C.gD[p], w, wk, C.gD[w], mpq_get_d(lhs), mpq_get_d(rhs));
                }
            }
        }
        mpq_clears(lhs, rhs, tmp, q1, q2, lamq, gval, d0, d1, best, cand, NULL);
    }
    printf("  [2/5] g checks:  %10lld  fails: %lld  (%s)\n", g_checks, g_fails, g_fails==0?"PASS":"FAIL");
    fflush(stdout);

    /* ===== CHECK 3: ALL VERTICAL 2-PATH INEQUALITIES h_V (EXACT mpz_t) ===== */
    /* Common exact integer denominator per (p, u, w):
     *   DEN_r = 2 * Dp_eff * Du * Dw_eff * (D0 + 1) * S_r
     *   Then multiply entire inequality by Q * DEN_r * S -> EVERY term is an exact integer (mpz_t)!
     *   (Wait: each hV_q, g_q already has denominator S, so multiplying (LHS - (P/Q)*hV) by S * Q * DEN_r
     *    cancels S cleanly and gives exact integer LHS_Z vs RHS_Z = P * DEN_r * hV_q!)
     */
    long long hv_checks = 0, hv_fails = 0;
    #pragma omp parallel num_threads(3) reduction(+:hv_checks,hv_fails)
    {
        mpz_t lhs_z, rhs_z, tmp_z;
        mpz_inits(lhs_z, rhs_z, tmp_z, NULL);

        #pragma omp for schedule(dynamic)
        for (int ua = 0; ua < C.NHU_ACT; ua++) {
            int ui = C.hu_list[ua];
            int ku = (ui / C.DHU) + 1, Lu = ui % C.DHU, Du = Lu + ku + 1;
            int gu = cl_inner_idx(&C, Lu, ku, 1);
            int gup = (Du + 1 < D0) ? cl_inner_idx(&C, Lu + 1, ku, 1) : cl_t1_idx(&C, ku, 1);
            int ui_up  = cl_hu_idx(&C, Lu + 1, ku, 1);
            int ui_bud = (Lu >= 1) ? cl_hu_idx(&C, Lu - 1, ku + 1, 1) : -1;
            int64_t xp_q = xpar_q[gu], xw_q = xch_q[gu];

            for (int ip = 0; ip < C.npars; ip++) {
                int p = C.pars[ip], pk = C.gKind[p], Lp = C.gL[p], kp = C.gK[p], hp = C.gH[p];
                long Dp_eff = (pk == 1) ? 1 : C.gD[p]; /* 1 for INF, D_p for inner, D0 for T1 */
                int p_up  = (pk == 0) ? ((C.gD[p]+1 < D0) ? cl_inner_idx(&C,Lp+1,kp,hp) : cl_t1_idx(&C,kp,hp)) : C.INF_ID;
                int p_bud = (pk == 0 && Lp>=1) ? cl_inner_idx(&C,Lp-1,kp+1,hp) : ((pk == 2 && kp <= KT) ? cl_t1_idx(&C,kp+1,hp) : p);
                int64_t dgpu = g_q[p * N + gup] - g_q[p * N + gu];
                int mp_t1 = (pk == 2 && kp <= KT && D0 - kp - hp > 0) ? (D0 - kp - hp) : 0;

                for (int iw = 0; iw < C.nchs; iw++) {
                    int w = C.chs[iw], wk = C.gKind[w], Lw = C.gL[w], kw = C.gK[w];
                    long Dw_eff = (wk == 1) ? 1 : C.gD[w];
                    int w_up  = (wk == 0) ? ((C.gD[w]+1 < D0) ? cl_inner_idx(&C,Lw+1,kw,1) : cl_t1_idx(&C,kw,1)) : C.INF_ID;
                    int w_bud = (wk == 0 && Lw>=1) ? cl_inner_idx(&C,Lw-1,kw+1,1) : ((wk == 2 && kw <= KT) ? cl_t1_idx(&C,kw+1,1) : w);
                    int64_t dguw = g_q[gup * N + w] - g_q[gu * N + w];
                    int mw_t1 = (wk == 2 && kw <= KT && D0 - kw - 1 > 0) ? (D0 - kw - 1) : 0;

                    int64_t v = hV_q[cl_hv_off(&C, p, ui, w)];
                    int64_t hp_up  = hV_q[cl_hv_off(&C, p_up, ui, w)];
                    int64_t hu_up  = (ui_up >= 0) ? hV_q[cl_hv_off(&C, p, ui_up, w)] : 0;
                    int64_t hw_up  = hV_q[cl_hv_off(&C, p, ui, w_up)];
                    int64_t hp_bud = (pk == 0 && Lp >= 1) || (pk == 2 && kp <= KT) ? hV_q[cl_hv_off(&C, p_bud, ui, w)] : 0;
                    int64_t hu_bud = (ui_bud >= 0) ? hV_q[cl_hv_off(&C, p, ui_bud, w)] : 0;
                    int64_t hw_bud = (wk == 0 && Lw >= 1) || (wk == 2 && kw <= KT) ? hV_q[cl_hv_off(&C, p, ui, w_bud)] : 0;

                    /* Exact integer common denominator D_all = 2 * Dp_eff * Du * Dw_eff * (D0+1) * S_r
                     * Check: 2*24*16*24*25*1e6 = 4.608e11 < 9.22e18 (fits in signed 64-bit long!)
                     */
                    long D_base = 2L * Dp_eff * Du * Dw_eff * (D0 + 1);
                    long D_all  = D_base * H.S_r;

                    /* Accumulate exact integer sum = (LHS * S) * D_all into mpz_t lhs_z! */
                    mpz_set_ui(lhs_z, 0);

                    /* Source 1: (inv_p - xp_q/S_r) * dguw */
                    long c_s1 = ((pk == 1) ? 0 : (D_all / Dp_eff)) - D_base * xp_q;
                    mpz_set_si(tmp_z, (long)dguw); mpz_mul_si(tmp_z, tmp_z, c_s1); mpz_add(lhs_z, lhs_z, tmp_z);

                    /* Source 2: (inv_w - xw_q/S_r) * dgpu */
                    long c_s2 = ((wk == 1) ? 0 : (D_all / Dw_eff)) - D_base * xw_q;
                    mpz_set_si(tmp_z, (long)dgpu); mpz_mul_si(tmp_z, tmp_z, c_s2); mpz_add(lhs_z, lhs_z, tmp_z);

                    /* p leaf gain & budding */
                    if (pk == 0) {
                        /* rate = Lp + 1/Du + ext(0, (hp+kp-1)/2) */
                        long r_lo = Lp * D_all + D_all / Du;
                        long r_hi = r_lo + (hp + kp - 1) * (D_all / 2);
                        long r = choose_rate(H.mode, hp_up, v, r_lo, r_hi);
                        mpz_set_si(tmp_z, (long)(hp_up - v)); mpz_mul_si(tmp_z, tmp_z, r); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (Lp >= 1) {
                            long rb = Lp * (D_all / Dp_eff);
                            mpz_set_si(tmp_z, (long)(hp_bud - v)); mpz_mul_si(tmp_z, tmp_z, rb); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                    } else if (pk == 2) {
                        long g_lo = mp_t1 * (D_all / (D0 + 1)), g_hi = D_all;
                        long gr = choose_rate(H.mode, hp_up, v, g_lo, g_hi);
                        mpz_set_si(tmp_z, (long)(hp_up - v)); mpz_mul_si(tmp_z, tmp_z, gr); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (kp <= KT) {
                            long b_lo = mp_t1 * (D_all / D0), b_hi = D_all;
                            long br = choose_rate(H.mode, hp_bud, v, b_lo, b_hi);
                            mpz_set_si(tmp_z, (long)(hp_bud - v)); mpz_mul_si(tmp_z, tmp_z, br); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                        int64_t d0_pu = ((ui_up>=0)?hV_q[cl_hv_off(&C,C.INF_ID,ui_up,w)]:0) - hV_q[cl_hv_off(&C,C.INF_ID,ui,w)];
                        int64_t d1_pu = hu_up - v;
                        int64_t dext = (H.mode == 1) ? (d1_pu > d0_pu ? d1_pu : d0_pu) : (d1_pu < d0_pu ? d1_pu : d0_pu);
                        mpz_set_si(tmp_z, (long)dext); mpz_mul_si(tmp_z, tmp_z, D_all / D0); mpz_add(lhs_z, lhs_z, tmp_z);
                    }

                    /* w leaf gain & budding */
                    if (wk == 0) {
                        long r_lo = Lw * D_all + D_all / Du;
                        long r_hi = r_lo + kw * (D_all / 2);
                        long r = choose_rate(H.mode, hw_up, v, r_lo, r_hi);
                        mpz_set_si(tmp_z, (long)(hw_up - v)); mpz_mul_si(tmp_z, tmp_z, r); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (Lw >= 1) {
                            long rb = Lw * (D_all / Dw_eff);
                            mpz_set_si(tmp_z, (long)(hw_bud - v)); mpz_mul_si(tmp_z, tmp_z, rb); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                    } else if (wk == 2) {
                        long g_lo = mw_t1 * (D_all / (D0 + 1)), g_hi = D_all;
                        long gr = choose_rate(H.mode, hw_up, v, g_lo, g_hi);
                        mpz_set_si(tmp_z, (long)(hw_up - v)); mpz_mul_si(tmp_z, tmp_z, gr); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (kw <= KT) {
                            long b_lo = mw_t1 * (D_all / D0), b_hi = D_all;
                            long br = choose_rate(H.mode, hw_bud, v, b_lo, b_hi);
                            mpz_set_si(tmp_z, (long)(hw_bud - v)); mpz_mul_si(tmp_z, tmp_z, br); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                        int64_t d0_wu = ((ui_up>=0)?hV_q[cl_hv_off(&C,p,ui_up,C.INF_ID)]:0) - hV_q[cl_hv_off(&C,p,ui,C.INF_ID)];
                        int64_t d1_wu = hu_up - v;
                        int64_t dext = (H.mode == 1) ? (d1_wu > d0_wu ? d1_wu : d0_wu) : (d1_wu < d0_wu ? d1_wu : d0_wu);
                        mpz_set_si(tmp_z, (long)dext); mpz_mul_si(tmp_z, tmp_z, D_all / D0); mpz_add(lhs_z, lhs_z, tmp_z);
                    }

                    /* u leaf gain & budding */
                    {
                        long ru_lo = Lu * D_all + (pk == 0 ? (D_all / Dp_eff) : 0) + (wk == 0 ? (D_all / Dw_eff) : 0);
                        long ru_hi = ru_lo + (ku - 1) * (D_all / 2);
                        long ru = choose_rate(H.mode, hu_up, v, ru_lo, ru_hi);
                        mpz_set_si(tmp_z, (long)(hu_up - v)); mpz_mul_si(tmp_z, tmp_z, ru); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (Lu >= 1) {
                            long rb = Lu * (D_all / Du);
                            mpz_set_si(tmp_z, (long)(hu_bud - v)); mpz_mul_si(tmp_z, tmp_z, rb); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                    }

                    /* Compare Q * lhs_z vs P * D_all * v */
                    mpz_mul_si(lhs_z, lhs_z, H.Q);
                    mpz_set_si(rhs_z, (long)v); mpz_mul_si(rhs_z, rhs_z, D_all); mpz_mul_si(rhs_z, rhs_z, H.P);
                    hv_checks++;
                    if (!check_sign(H.mode, lhs_z, rhs_z)) {
                        hv_fails++;
                        #pragma omp critical
                        if (hv_fails <= 5) fprintf(stderr, "FAIL hV(p=%d, u=(L=%d,k=%d), w=%d)\n", p, Lu, ku, w);
                    }
                }
            }
        }
        mpz_clears(lhs_z, rhs_z, tmp_z, NULL);
    }
    printf("  [3/5] hV checks: %10lld  fails: %lld  (%s)\n", hv_checks, hv_fails, hv_fails==0?"PASS":"FAIL");
    fflush(stdout);

    /* ===== CHECK 4: ALL CHERRY INEQUALITIES h_C (EXACT mpz_t) ===== */
    long long hc_checks = 0, hc_fails = 0;
    #pragma omp parallel num_threads(3) reduction(+:hc_checks,hc_fails)
    {
        mpz_t lhs_z, rhs_z, tmp_z;
        mpz_inits(lhs_z, rhs_z, tmp_z, NULL);

        #pragma omp for schedule(dynamic)
        for (int ca = 0; ca < C.NCU_ACT; ca++) {
            int ci = C.cu_list[ca];
            int block = (C.KCU - 1) * C.DCU;
            int hu = (ci < block) ? 1 : 0;
            int rem = (hu == 1) ? ci : (ci - block);
            int ku = (rem / C.DCU) + 2, Lu = rem % C.DCU, Du = Lu + ku + hu;
            int gu = cl_inner_idx(&C, Lu, ku, hu);
            int gup = (Du + 1 < D0) ? cl_inner_idx(&C, Lu + 1, ku, hu) : cl_t1_idx(&C, ku, hu);
            int ci_up  = cl_cu_idx(&C, Lu + 1, ku, hu);
            int ci_bud = (Lu >= 1) ? cl_cu_idx(&C, Lu - 1, ku + 1, hu) : -1;
            int64_t xc_q = xcc_q[gu];

            for (int i1 = 0; i1 < C.nchs; i1++) {
                int w1 = C.chs[i1], k1_kind = C.gKind[w1], L1 = C.gL[w1], k1 = C.gK[w1];
                long D1_eff = (k1_kind == 1) ? 1 : C.gD[w1];
                int w1_up  = (k1_kind == 0) ? ((C.gD[w1]+1<D0)?cl_inner_idx(&C,L1+1,k1,1):cl_t1_idx(&C,k1,1)) : C.INF_ID;
                int w1_bud = (k1_kind == 0 && L1>=1) ? cl_inner_idx(&C,L1-1,k1+1,1) : ((k1_kind==2&&k1<=KT)?cl_t1_idx(&C,k1+1,1):w1);
                int64_t dg1 = g_q[gup * N + w1] - g_q[gu * N + w1];
                int m1_t1 = (k1_kind == 2 && k1 <= KT && D0 - k1 - 1 > 0) ? (D0 - k1 - 1) : 0;

                for (int i2 = i1; i2 < C.nchs; i2++) {
                    int w2 = C.chs[i2], k2_kind = C.gKind[w2], L2 = C.gL[w2], k2 = C.gK[w2];
                    long D2_eff = (k2_kind == 1) ? 1 : C.gD[w2];
                    int w2_up  = (k2_kind == 0) ? ((C.gD[w2]+1<D0)?cl_inner_idx(&C,L2+1,k2,1):cl_t1_idx(&C,k2,1)) : C.INF_ID;
                    int w2_bud = (k2_kind == 0 && L2>=1) ? cl_inner_idx(&C,L2-1,k2+1,1) : ((k2_kind==2&&k2<=KT)?cl_t1_idx(&C,k2+1,1):w2);
                    int64_t dg2 = g_q[gup * N + w2] - g_q[gu * N + w2];
                    int m2_t1 = (k2_kind == 2 && k2 <= KT && D0 - k2 - 1 > 0) ? (D0 - k2 - 1) : 0;

                    int64_t v = hC_q[cl_hc_off(&C, ci, w1, w2)];
                    int64_t cu_up  = (ci_up >= 0)  ? hC_q[cl_hc_off(&C, ci_up, w1, w2)] : 0;
                    int64_t cu_bud = (ci_bud >= 0) ? hC_q[cl_hc_off(&C, ci_bud, w1, w2)] : 0;
                    int64_t c1_up  = hC_q[cl_hc_off(&C, ci, w1_up, w2)];
                    int64_t c2_up  = hC_q[cl_hc_off(&C, ci, w1, w2_up)];
                    int64_t c1_bud = (k1_kind==0&&L1>=1)||(k1_kind==2&&k1<=KT) ? hC_q[cl_hc_off(&C, ci, w1_bud, w2)] : 0;
                    int64_t c2_bud = (k2_kind==0&&L2>=1)||(k2_kind==2&&k2<=KT) ? hC_q[cl_hc_off(&C, ci, w1, w2_bud)] : 0;

                    long D_base = 2L * D1_eff * Du * D2_eff * (D0 + 1);
                    long D_all  = D_base * H.S_r;
                    mpz_set_ui(lhs_z, 0);

                    /* Sources: (inv1 - xc)*dg2 + (inv2 - xc)*dg1 */
                    long c_s1 = ((k1_kind == 1) ? 0 : (D_all / D1_eff)) - D_base * xc_q;
                    mpz_set_si(tmp_z, (long)dg2); mpz_mul_si(tmp_z, tmp_z, c_s1); mpz_add(lhs_z, lhs_z, tmp_z);
                    long c_s2 = ((k2_kind == 1) ? 0 : (D_all / D2_eff)) - D_base * xc_q;
                    mpz_set_si(tmp_z, (long)dg1); mpz_mul_si(tmp_z, tmp_z, c_s2); mpz_add(lhs_z, lhs_z, tmp_z);

                    /* w1 leaf gain & budding */
                    if (k1_kind == 0) {
                        long r_lo = L1 * D_all + D_all / Du, r_hi = r_lo + k1 * (D_all / 2);
                        long r = choose_rate(H.mode, c1_up, v, r_lo, r_hi);
                        mpz_set_si(tmp_z, (long)(c1_up - v)); mpz_mul_si(tmp_z, tmp_z, r); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (L1 >= 1) {
                            long rb = L1 * (D_all / D1_eff);
                            mpz_set_si(tmp_z, (long)(c1_bud - v)); mpz_mul_si(tmp_z, tmp_z, rb); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                    } else if (k1_kind == 2) {
                        long gr = choose_rate(H.mode, c1_up, v, m1_t1 * (D_all / (D0 + 1)), D_all);
                        mpz_set_si(tmp_z, (long)(c1_up - v)); mpz_mul_si(tmp_z, tmp_z, gr); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (k1 <= KT) {
                            long br = choose_rate(H.mode, c1_bud, v, m1_t1 * (D_all / D0), D_all);
                            mpz_set_si(tmp_z, (long)(c1_bud - v)); mpz_mul_si(tmp_z, tmp_z, br); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                        int64_t d0 = ((ci_up>=0)?hC_q[cl_hc_off(&C,ci_up,C.INF_ID,w2)]:0) - hC_q[cl_hc_off(&C,ci,C.INF_ID,w2)];
                        int64_t d1 = cu_up - v;
                        int64_t dext = (H.mode == 1) ? (d1 > d0 ? d1 : d0) : (d1 < d0 ? d1 : d0);
                        mpz_set_si(tmp_z, (long)dext); mpz_mul_si(tmp_z, tmp_z, D_all / D0); mpz_add(lhs_z, lhs_z, tmp_z);
                    }

                    /* w2 leaf gain & budding */
                    if (k2_kind == 0) {
                        long r_lo = L2 * D_all + D_all / Du, r_hi = r_lo + k2 * (D_all / 2);
                        long r = choose_rate(H.mode, c2_up, v, r_lo, r_hi);
                        mpz_set_si(tmp_z, (long)(c2_up - v)); mpz_mul_si(tmp_z, tmp_z, r); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (L2 >= 1) {
                            long rb = L2 * (D_all / D2_eff);
                            mpz_set_si(tmp_z, (long)(c2_bud - v)); mpz_mul_si(tmp_z, tmp_z, rb); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                    } else if (k2_kind == 2) {
                        long gr = choose_rate(H.mode, c2_up, v, m2_t1 * (D_all / (D0 + 1)), D_all);
                        mpz_set_si(tmp_z, (long)(c2_up - v)); mpz_mul_si(tmp_z, tmp_z, gr); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (k2 <= KT) {
                            long br = choose_rate(H.mode, c2_bud, v, m2_t1 * (D_all / D0), D_all);
                            mpz_set_si(tmp_z, (long)(c2_bud - v)); mpz_mul_si(tmp_z, tmp_z, br); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                        int64_t d0 = ((ci_up>=0)?hC_q[cl_hc_off(&C,ci_up,w1,C.INF_ID)]:0) - hC_q[cl_hc_off(&C,ci,w1,C.INF_ID)];
                        int64_t d1 = cu_up - v;
                        int64_t dext = (H.mode == 1) ? (d1 > d0 ? d1 : d0) : (d1 < d0 ? d1 : d0);
                        mpz_set_si(tmp_z, (long)dext); mpz_mul_si(tmp_z, tmp_z, D_all / D0); mpz_add(lhs_z, lhs_z, tmp_z);
                    }

                    /* u leaf gain & budding */
                    {
                        long ru_lo = Lu * D_all + (k1_kind == 0 ? (D_all / D1_eff) : 0) + (k2_kind == 0 ? (D_all / D2_eff) : 0);
                        long ru_hi = ru_lo + (hu + ku - 2) * (D_all / 2);
                        long ru = choose_rate(H.mode, cu_up, v, ru_lo, ru_hi);
                        mpz_set_si(tmp_z, (long)(cu_up - v)); mpz_mul_si(tmp_z, tmp_z, ru); mpz_add(lhs_z, lhs_z, tmp_z);
                        if (Lu >= 1) {
                            long rb = Lu * (D_all / Du);
                            mpz_set_si(tmp_z, (long)(cu_bud - v)); mpz_mul_si(tmp_z, tmp_z, rb); mpz_add(lhs_z, lhs_z, tmp_z);
                        }
                    }

                    mpz_mul_si(lhs_z, lhs_z, H.Q);
                    mpz_set_si(rhs_z, (long)v); mpz_mul_si(rhs_z, rhs_z, D_all); mpz_mul_si(rhs_z, rhs_z, H.P);
                    hc_checks++;
                    if (!check_sign(H.mode, lhs_z, rhs_z)) {
                        hc_fails++;
                        #pragma omp critical
                        if (hc_fails <= 5) fprintf(stderr, "FAIL hC(u=(L=%d,k=%d,h=%d), w1=%d, w2=%d)\n", Lu, ku, hu, w1, w2);
                    }
                }
            }
        }
        mpz_clears(lhs_z, rhs_z, tmp_z, NULL);
    }
    printf("  [4/5] hC checks: %10lld  fails: %lld  (%s)\n", hc_checks, hc_fails, hc_fails==0?"PASS":"FAIL");
    fflush(stdout);

    /* ===== CHECK 5: GLOBAL POINTWISE VERTEX-GROUPED COMPARISON (EXACT int64_t / mpq_t) ===== */
    int64_t *g_ext  = calloc(N, sizeof(int64_t));
    int64_t *hv_ext = calloc(N, sizeof(int64_t));
    int64_t *hc_ext = calloc(C.NCU, sizeof(int64_t));
    for (int iw = 0; iw < C.nchs; iw++) {
        int w = C.chs[iw];
        int64_t ev = 0;
        for (int ip = 0; ip < C.npars; ip++) {
            int64_t v = g_q[C.pars[ip] * N + w];
            if ((H.mode == 1 && v < ev) || (H.mode == 2 && v > ev)) ev = v;
        }
        g_ext[w] = ev;
        int64_t he = 0;
        for (int ip = 0; ip < C.npars; ip++)
            for (int ua = 0; ua < C.NHU_ACT; ua++) {
                int64_t v = hV_q[cl_hv_off(&C, C.pars[ip], C.hu_list[ua], w)];
                if ((H.mode == 1 && v < he) || (H.mode == 2 && v > he)) he = v;
            }
        hv_ext[w] = he;
    }
    for (int ca = 0; ca < C.NCU_ACT; ca++) {
        int ci = C.cu_list[ca];
        int64_t ce = 0;
        for (int i1 = 0; i1 < C.nchs; i1++)
            for (int i2 = i1; i2 < C.nchs; i2++) {
                int64_t v = hC_q[cl_hc_off(&C, ci, C.chs[i1], C.chs[i2])];
                if ((H.mode == 1 && v < ce) || (H.mode == 2 && v > ce)) ce = v;
            }
        hc_ext[ci] = ce;
    }
    int64_t worst_W = (H.mode == 1) ? INT64_MAX : INT64_MIN;
    for (int h = 0; h <= 1; h++) for (int k = 0; k <= C.KF; k++) for (int L = 0; L <= C.LF; L++) {
        int D = L + k + h; if (D < 2) continue;
        int64_t W = f_q[cl_vidx(&C, L, k, h)];
        if (h == 1) {
            if (D < D0) { int w = cl_inner_idx(&C, L, k, 1); W += g_ext[w] + hv_ext[w]; }
            else {
                int w1 = C.INF_ID, w2 = cl_t1_idx(&C, k, 1);
                W += (H.mode == 1) ? (g_ext[w1] < g_ext[w2] ? g_ext[w1] : g_ext[w2]) : (g_ext[w1] > g_ext[w2] ? g_ext[w1] : g_ext[w2]);
                W += (H.mode == 1) ? (hv_ext[w1] < hv_ext[w2] ? hv_ext[w1] : hv_ext[w2]) : (hv_ext[w1] > hv_ext[w2] ? hv_ext[w1] : hv_ext[w2]);
            }
        }
        int ci = cl_cu_idx(&C, L, k, h);
        if (ci >= 0) W += ((int64_t)k * (k - 1) / 2) * hc_ext[ci];
        if ((H.mode == 1 && W < worst_W) || (H.mode == 2 && W > worst_W)) worst_W = W;
    }
    /* VSTAR0 and VSTAR1 */
    {
        int64_t W0 = f_q[C.VSTAR0];
        if ((H.mode == 1 && W0 < worst_W) || (H.mode == 2 && W0 > worst_W)) worst_W = W0;
        int w1 = C.INF_ID, w2 = cl_t1_idx(&C, KT + 1, 1);
        int64_t W1 = f_q[C.VSTAR1]
            + ((H.mode == 1) ? (g_ext[w1] < g_ext[w2] ? g_ext[w1] : g_ext[w2]) : (g_ext[w1] > g_ext[w2] ? g_ext[w1] : g_ext[w2]))
            + ((H.mode == 1) ? (hv_ext[w1] < hv_ext[w2] ? hv_ext[w1] : hv_ext[w2]) : (hv_ext[w1] > hv_ext[w2] ? hv_ext[w1] : hv_ext[w2]));
        if ((H.mode == 1 && W1 < worst_W) || (H.mode == 2 && W1 > worst_W)) worst_W = W1;
    }
    int64_t phi_T3_q = f_q[cl_vidx(&C, 2, 0, 0)];
    int ok_global = (worst_W > 0) && (phi_T3_q > 0);
    mpq_t bound_c, prefactor; mpq_inits(bound_c, prefactor, NULL);
    mpq_set_si(bound_c, (long)worst_W, (unsigned long)H.S); mpq_canonicalize(bound_c);
    mpq_set_si(prefactor, (long)phi_T3_q, (unsigned long)worst_W); mpq_canonicalize(prefactor);
    gmp_printf("  [5/5] Global comparison: Phi(T3)=%lld/%lld, %s W(u)=%lld/%lld, exact prefactor A = %Qd (%.8f)  (%s)\n",
               (long long)phi_T3_q, (long long)H.S, H.mode==1?"min":"max", (long long)worst_W, (long long)H.S,
               prefactor, mpq_get_d(prefactor), ok_global ? "PASS" : "FAIL");

    /* ===== BONUS CHECK 6: EXACT RATIONAL BRUTE-FORCE TREE DRIFT CHECKS ===== */
    int tree_fails = 0, tree_cnt = 0;
    if (!skip_trees && f_fails == 0 && g_fails == 0 && hv_fails == 0 && hc_fails == 0 && ok_global) {
        /* 1. T_3 (star K_{1,2}) */
        { int par[3]={-1,0,0}, deg[3]={2,1,1}; tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,3,par,deg,"T3")) tree_fails++; }
        /* 2. T_4 star K_{1,3} and path P_4 */
        { int par[4]={-1,0,0,0}, deg[4]={3,1,1,1}; tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,4,par,deg,"K1,3")) tree_fails++; }
        { int par[4]={-1,0,0,1}, deg[4]={2,2,1,1}; tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,4,par,deg,"P4")) tree_fails++; }
        /* 3. Extreme adversarial trees hitting D=D0-1, D=D0, D=40, L=0 chains, k>KT, L=LF+2, k=KF+2! */
        /* 3a. Hub star K_{1,30} with two core children of degrees D0-1 and D0+5 */
        {
            int n = 1 + 30 + 25 + 30;
            int *par = calloc(n, sizeof(int)), *deg = calloc(n, sizeof(int));
            par[0] = -1; int idx = 1;
            int c1 = idx++; par[c1] = 0; deg[0]++; deg[c1]++;
            int c2 = idx++; par[c2] = 0; deg[0]++; deg[c2]++;
            for (int i = 0; i < 28; i++) { par[idx] = 0; deg[0]++; deg[idx++] = 1; }
            for (int i = 0; i < D0 - 2; i++) { par[idx] = c1; deg[c1]++; deg[idx++] = 1; } /* D(c1)=D0-1 */
            for (int i = 0; i < D0 + 4; i++) { par[idx] = c2; deg[c2]++; deg[idx++] = 1; } /* D(c2)=D0+5 */
            tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,idx,par,deg,"HubTail_D0")) tree_fails++;
            free(par); free(deg);
        }
        /* 3b. Dead core chain of L=0 vertices ending in leaves + cherry */
        {
            int par[10] = {-1, 0, 0, 1, 3, 4, 2, 6, 4, 6};
            int deg[10] = {0};
            for (int i = 1; i < 10; i++) { deg[i]++; deg[par[i]]++; }
            tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,10,par,deg,"DeadChain")) tree_fails++;
        }
        /* 3c. VSTAR boundary: root with k = KF + 2 core children (each deg 2) + 2 leaves! */
        {
            int n = 1 + 2 * (C.KF + 2) + 2;
            int *par = calloc(n, sizeof(int)), *deg = calloc(n, sizeof(int));
            par[0] = -1; int idx = 1;
            for (int i = 0; i < C.KF + 2; i++) {
                int ch = idx++; par[ch] = 0; deg[0]++; deg[ch]++;
                int lf = idx++; par[lf] = ch; deg[ch]++; deg[lf]++;
            }
            for (int i = 0; i < 2; i++) { par[idx] = 0; deg[0]++; deg[idx++] = 1; }
            tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,idx,par,deg,"VSTAR_k_gt_KF")) tree_fails++;
            free(par); free(deg);
        }
        /* 3d. LF boundary: root with L = LF + 2 leaves! */
        {
            int n = 1 + C.LF + 2;
            int *par = calloc(n, sizeof(int)), *deg = calloc(n, sizeof(int));
            par[0] = -1;
            for (int i = 1; i < n; i++) { par[i] = 0; deg[0]++; deg[i] = 1; }
            tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,n,par,deg,"LeafTail_L_gt_LF")) tree_fails++;
            free(par); free(deg);
        }
        /* 4. 30 random friend trees of sizes n = 5 .. 45 */
        uint64_t rng = 31415926535ULL;
        for (int t = 0; t < 30; t++) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            int n = 5 + (int)((rng >> 33) % 41);
            int par[64], deg[64], kids[64][64], kcnt[64] = {0};
            par[0] = -1; deg[0] = 2; par[1] = 0; deg[1] = 1; par[2] = 0; deg[2] = 1;
            kids[0][0] = 1; kids[0][1] = 2; kcnt[0] = 2;
            for (int cur = 3; cur < n; cur++) {
                rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
                int V = (int)((rng >> 33) % cur);
                rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
                int r = (int)((rng >> 33) % deg[V]);
                int W = (par[V] != -1) ? ((r == 0) ? par[V] : kids[V][r - 1]) : kids[V][r];
                par[cur] = W; deg[cur] = 1; deg[W]++; kids[W][kcnt[W]++] = cur;
            }
            char lbl[32]; snprintf(lbl, sizeof(lbl), "Rand_%d_n%d", t, n);
            tree_cnt++; if (!check_tree_exact(&H,&C,f_q,g_q,hV_q,hC_q,bound_c,n,par,deg,lbl)) tree_fails++;
        }
        printf("  [6/6] Exact brute-force trees: %d tested, fails: %d (%s)\n", tree_cnt, tree_fails, tree_fails==0?"PASS":"FAIL");
    }

    long long total_checks = f_checks + g_checks + hv_checks + hc_checks;
    long long total_fails  = f_fails + g_fails + hv_fails + hc_fails + (!ok_global) + tree_fails;
    printf("=== SUMMARY: %lld exact GMP inequality checks, %lld failures => %s ===\n",
           total_checks, total_fails, total_fails == 0 ? "VERIFIED OK" : "REJECTED");

    mpq_clears(bound_c, prefactor, NULL);
    for (size_t i = 0; i < (size_t)N * N; i++) mpq_clear(A_tot[i]);
    free(A_tot); free(dLf_q); free(g_ext); free(hv_ext); free(hc_ext);
    free(rpar_q); free(rch_q); free(xpar_q); free(xch_q); free(xcc_q);
    free(k_pc_q); free(k_cp_q); free(k_cc_q); free(f_q); free(g_q); free(hV_q); free(hC_q);
    cl_free(&C);
    return (total_fails == 0) ? 0 : 2;
}
