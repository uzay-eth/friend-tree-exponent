/*
 * Standalone Pure-GMP Exact-Arithmetic Verifier for B3-FullH + Level-4 Drift Certificates
 *
 * ZERO floating point is used anywhere in the mathematical verification!
 * Every single inequality across all 8 motif tables (f, g, hV, hC, P4, Ytop, Ybot, S3),
 * the comparison bounds c0 C(T) <= Phi(T) <= C0 C(T), Phi(T_3) > 0, and a built-in
 * exact-rational brute-force drift verification on concrete trees are checked in
 * arbitrary-precision GMP integer/rational arithmetic (mpz_t / mpq_t).
 *
 * Usage:
 *   cert_b3f_verify_gmp <cert.bin> [--threads 4] [--override-pq P Q] [--tree-check N]
 *
 * Exit code:
 *   0  iff ALL checks pass (header, file size, g(INF,INF)==0, all 8 drift tables,
 *      c0 > 0, C0 > 0, Phi(T_3) > 0, and exact-rational concrete tree checks).
 *   1  if ANY check fails.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>
#include <gmp.h>
#include <omp.h>
#include "cert_b3f_l5_format.h"

/* Exact fraction with strictly positive denominator den > 0 (no gcd in hot loops). */
typedef struct { mpz_t num, den; } ZFrac;

typedef struct {
    ZFrac acc, best, cur, tmp;
    mpz_t z1, z2, z3, z4; /* private scratch exclusively for zf_* primitives */
    mpz_t t1, t2, t3, t4, M, A0, A1, E0, E1;
    mpz_t P[5], BP[5], BB[5], sh_n, sh_d;
} GmpWork;

static void gw_init(GmpWork *W) {
    mpz_init(W->acc.num);  mpz_init(W->acc.den);
    mpz_init(W->best.num); mpz_init(W->best.den);
    mpz_init(W->cur.num);  mpz_init(W->cur.den);
    mpz_init(W->tmp.num);  mpz_init(W->tmp.den);
    mpz_init(W->z1); mpz_init(W->z2); mpz_init(W->z3); mpz_init(W->z4);
    mpz_init(W->t1); mpz_init(W->t2); mpz_init(W->t3); mpz_init(W->t4);
    mpz_init(W->M); mpz_init(W->A0); mpz_init(W->A1); mpz_init(W->E0); mpz_init(W->E1);
    for (int i = 0; i < 5; i++) { mpz_init(W->P[i]); mpz_init(W->BP[i]); mpz_init(W->BB[i]); }
    mpz_init(W->sh_n); mpz_init(W->sh_d);
}
static void gw_clear(GmpWork *W) {
    mpz_clear(W->acc.num);  mpz_clear(W->acc.den);
    mpz_clear(W->best.num); mpz_clear(W->best.den);
    mpz_clear(W->cur.num);  mpz_clear(W->cur.den);
    mpz_clear(W->tmp.num);  mpz_clear(W->tmp.den);
    mpz_clear(W->z1); mpz_clear(W->z2); mpz_clear(W->z3); mpz_clear(W->z4);
    mpz_clear(W->t1); mpz_clear(W->t2); mpz_clear(W->t3); mpz_clear(W->t4);
    mpz_clear(W->M); mpz_clear(W->A0); mpz_clear(W->A1); mpz_clear(W->E0); mpz_clear(W->E1);
    for (int i = 0; i < 5; i++) { mpz_clear(W->P[i]); mpz_clear(W->BP[i]); mpz_clear(W->BB[i]); }
    mpz_clear(W->sh_n); mpz_clear(W->sh_d);
}

static inline void zf_set_si(ZFrac *F, int64_t n, uint64_t d) {
    mpz_set_si(F->num, n); mpz_set_ui(F->den, d);
}
static inline void zf_add_z(ZFrac *F, const mpz_t n2, const mpz_t d2, GmpWork *W) {
    if (mpz_sgn(n2) == 0) return;
    if (mpz_cmp(F->den, d2) == 0) {
        mpz_add(F->num, F->num, n2);
    } else {
        mpz_mul(W->z1, F->num, d2);
        mpz_mul(W->z2, n2, F->den);
        mpz_add(F->num, W->z1, W->z2);
        mpz_mul(F->den, F->den, d2);
    }
}
static inline void zf_add_si(ZFrac *F, int64_t n2, uint64_t d2, GmpWork *W) {
    if (n2 == 0) return;
    if (mpz_cmp_ui(F->den, d2) == 0) {
        if (n2 > 0) mpz_add_ui(F->num, F->num, (uint64_t)n2);
        else        mpz_sub_ui(F->num, F->num, (uint64_t)(-n2));
    } else {
        mpz_mul_ui(F->num, F->num, d2);
        mpz_set_si(W->z2, n2);
        mpz_addmul(F->num, W->z2, F->den);
        mpz_mul_ui(F->den, F->den, d2);
    }
}
static inline void zf_add_mul_si(ZFrac *F, int64_t a, int64_t b, uint64_t d2, GmpWork *W) {
    if (a == 0 || b == 0) return;
    mpz_set_si(W->z3, a);
    mpz_mul_si(W->z3, W->z3, b);
    mpz_set_ui(W->z4, d2);
    zf_add_z(F, W->z3, W->z4, W);
}
static inline int zf_cmp(const ZFrac *A, const ZFrac *B, GmpWork *W) {
    mpz_mul(W->z1, A->num, B->den);
    mpz_mul(W->z2, B->num, A->den);
    return mpz_cmp(W->z1, W->z2);
}
static inline int zf_fail(const ZFrac *A, int mode) {
    int s = mpz_sgn(A->num); /* den > 0 always */
    return (mode == 1) ? (s > 0) : (s < 0);
}
static inline int64_t ext_i64(int mode, int64_t a, int64_t b) {
    return (mode == 1) ? ((a >= b) ? a : b) : ((a <= b) ? a : b);
}

/* Exact invD(u) = num / den (den > 0) */
static inline void invD_rat(const B3FLayout *C, int u, int64_t *num, int64_t *den) {
    int k = C->kind[u];
    if (k == 0)      { *num = 1; *den = C->D[u]; }
    else if (k == 1) { *num = 0; *den = 1; }
    else if (k == 2) { *num = 1; *den = 3 * C->D0; }
    else if (k == 3) { *num = 2; *den = 3 * C->D0; }
    else             { *num = 1; *den = C->D0; }
}

/* Exact (invD(u) - rate_q / S_r) * diff_val added into F */
static inline void add_centered_source(ZFrac *F, const B3FLayout *C, int u, int64_t rate_q, int64_t S_r, int64_t diff_val, GmpWork *W) {
    if (diff_val == 0) return;
    int64_t inum, iden;
    invD_rat(C, u, &inum, &iden);
    /* (inum/iden - rate_q/S_r) * diff_val = (inum*S_r - rate_q*iden) * diff_val / (iden * S_r) */
    mpz_set_si(W->t3, inum);
    mpz_mul_si(W->t3, W->t3, S_r);
    mpz_set_si(W->t4, rate_q);
    mpz_mul_si(W->t4, W->t4, iden);
    mpz_sub(W->t3, W->t3, W->t4);
    mpz_mul_si(W->t3, W->t3, diff_val);
    mpz_set_ui(W->t4, (uint64_t)iden);
    mpz_mul_ui(W->t4, W->t4, (uint64_t)S_r);
    zf_add_z(F, W->t3, W->t4, W);
}

/* Exact r_in interval [lo_n/lo_d, hi_n/hi_d] for 1 neighbor y, or sum of 2 neighbors y1, y2 */
static inline void get_rin1(const B3FLayout *C, int y, int64_t *ln, int64_t *ld, int64_t *hn, int64_t *hd) {
    if (C->kind[y] == 0) { *ln = 1; *ld = C->D[y]; *hn = 1; *hd = C->D[y]; }
    else                 { *ln = 0; *ld = 1;       *hn = 1; *hd = C->D0; }
}
static inline void get_rin2(const B3FLayout *C, int y1, int y2, int64_t *ln, int64_t *ld, int64_t *hn, int64_t *hd) {
    int64_t l1n, l1d, h1n, h1d, l2n, l2d, h2n, h2d;
    get_rin1(C, y1, &l1n, &l1d, &h1n, &h1d);
    get_rin1(C, y2, &l2n, &l2d, &h2n, &h2d);
    *ln = l1n * l2d + l2n * l1d; *ld = l1d * l2d;
    *hn = h1n * h2d + h2n * h1d; *hd = h1d * h2d;
}

/*
 * Exact pure-GMP Bernstein-3 tail contribution for ONE tail vertex x (kind in {2,3,4} = T1, T2, T3).
 * Adds worst-case rational drift contribution (max for mode=1, min for mode=2) directly into acc!
 *   xfer_ext = 1 on g (for kx <= KT): T2 uses rho_in + Rbar_q/S_r; T3 subtracts (r_ext - Rbar_q/S_r)*Y_ref.
 *   xfer_ext = 0 on hV/hC/L4: T2 uses rho in [r_in_lo, r_in_hi + ext_half/2]; T3 has shift = 0.
 */
static void add_tail_gmp(ZFrac *acc, GmpWork *W, int mode, int D0, int KT, int64_t S_r,
                         int xk, int kx, int hx,
                         int64_t rlo_n, int64_t rlo_d, int64_t rhi_n, int64_t rhi_d,
                         int ext_half, int xfer_ext, int64_t Rbar_q, int64_t Y_ref,
                         int64_t F0, int64_t F1, int64_t F2, int64_t F3,
                         int64_t F1b, int64_t F2b, int64_t F3b,
                         int n_nbr_in, int64_t d0_nbr, int64_t d1_nbr, int64_t d2_nbr, int64_t F3_up_sum) {
    if (xk < 2) return;

    if (xk == 2) {
        /* T1: current value v = F1 */
        if (kx <= KT) {
            /* F0 + F1b - 2*F1 */
            mpz_set_si(W->t3, F0);
            mpz_set_si(W->t4, F1b); mpz_add(W->t3, W->t3, W->t4);
            mpz_set_si(W->t4, F1);  mpz_submul_ui(W->t3, W->t4, 2);
            mpz_set_ui(W->t4, 1);
            zf_add_z(acc, W->t3, W->t4, W);
        } else {
            /* overflow k > KT: ext(0, F0 - F1) */
            mpz_set_si(W->t3, F0);
            mpz_set_si(W->t4, F1);
            mpz_sub(W->t3, W->t3, W->t4);
            int s = mpz_sgn(W->t3);
            if ((mode == 1 && s > 0) || (mode == 2 && s < 0)) {
                mpz_set_ui(W->t4, 1);
                zf_add_z(acc, W->t3, W->t4, W);
            }
        }
        if (n_nbr_in > 0) zf_add_si(acc, d0_nbr, (uint64_t)(3 * D0), W);
        return;
    }

    if (xk == 3) {
        /* T2: current value v = F2 */
        if (kx <= KT) {
            int c = kx + hx;
            /* Fixed part: 2*F1 + F2b - 3*F2 - (c/D0)*(F1b - F1) */
            mpz_set_si(W->t3, F1);  mpz_mul_ui(W->t3, W->t3, 2);
            mpz_set_si(W->t4, F2b); mpz_add(W->t3, W->t3, W->t4);
            mpz_set_si(W->t4, F2);  mpz_submul_ui(W->t3, W->t4, 3);
            mpz_mul_ui(W->t3, W->t3, (uint64_t)D0);
            mpz_set_si(W->t4, F1b); mpz_set_si(W->t1, F1); mpz_sub(W->t4, W->t4, W->t1);
            mpz_submul_ui(W->t3, W->t4, (uint64_t)c);
            mpz_set_ui(W->t4, (uint64_t)D0);
            zf_add_z(acc, W->t3, W->t4, W);

            /* rho term: (-(c+1) + rho) * (F0 - F1) / D0 */
            mpz_set_si(W->t1, F0); mpz_set_si(W->t2, F1); mpz_sub(W->t1, W->t1, W->t2); /* d1 */
            if (mpz_sgn(W->t1) != 0) {
                /* Want max (mode=1) or min (mode=2) of rho * d1 -> pick hi_rho iff (mode==1 ? d1>0 : d1<0) */
                int pick_hi = (mode == 1) ? (mpz_sgn(W->t1) > 0) : (mpz_sgn(W->t1) < 0);
                int64_t rn = pick_hi ? rhi_n : rlo_n;
                int64_t rd = pick_hi ? rhi_d : rlo_d;
                if (xfer_ext) {
                    /* rho = rn/rd + Rbar_q/S_r; coeff = (-(c+1)*rd*S_r + rn*S_r + Rbar_q*rd) / (D0 * rd * S_r) */
                    mpz_set_si(W->t3, -(int64_t)(c + 1) * rd + rn);
                    mpz_mul_si(W->t3, W->t3, S_r);
                    mpz_set_si(W->t4, Rbar_q);
                    mpz_addmul_ui(W->t3, W->t4, (uint64_t)rd);
                    mpz_mul(W->t3, W->t3, W->t1);
                    mpz_set_ui(W->t4, (uint64_t)D0);
                    mpz_mul_ui(W->t4, W->t4, (uint64_t)rd);
                    mpz_mul_ui(W->t4, W->t4, (uint64_t)S_r);
                    zf_add_z(acc, W->t3, W->t4, W);
                } else {
                    /* rho = rn/rd + (pick_hi ? ext_half/2 : 0) */
                    int eh = pick_hi ? ext_half : 0;
                    /* (-(c+1)*2*rd + 2*rn + eh*rd) * d1 / (2 * D0 * rd) */
                    mpz_set_si(W->t3, (-(int64_t)(c + 1) * 2 + (int64_t)eh) * rd + 2 * rn);
                    mpz_mul(W->t3, W->t3, W->t1);
                    mpz_set_ui(W->t4, (uint64_t)(2 * D0) * (uint64_t)rd);
                    zf_add_z(acc, W->t3, W->t4, W);
                }
            }
        } else {
            /* overflow k > KT: ext(0, 2*(F1 - F2)) */
            mpz_set_si(W->t3, F1); mpz_set_si(W->t4, F2); mpz_sub(W->t3, W->t3, W->t4);
            mpz_mul_ui(W->t3, W->t3, 2);
            int s = mpz_sgn(W->t3);
            if ((mode == 1 && s > 0) || (mode == 2 && s < 0)) {
                mpz_set_ui(W->t4, 1);
                zf_add_z(acc, W->t3, W->t4, W);
            }
        }
        if (n_nbr_in > 0) {
            /* (3*d1_nbr - d0_nbr) / (3*D0) */
            mpz_set_si(W->t3, d1_nbr); mpz_mul_ui(W->t3, W->t3, 3);
            mpz_set_si(W->t4, d0_nbr); mpz_sub(W->t3, W->t3, W->t4);
            mpz_set_ui(W->t4, (uint64_t)(3 * D0));
            zf_add_z(acc, W->t3, W->t4, W);
        }
        return;
    }

    /* xk == 4 (T3): current value v = F3 */
    mpz_t d1z, d2z;
    mpz_init(d1z); mpz_init(d2z);
    mpz_set_si(d1z, F0); mpz_set_si(W->t1, F1); mpz_sub(d1z, d1z, W->t1);
    mpz_set_si(d2z, F1); mpz_set_si(W->t1, F2); mpz_sub(d2z, d2z, W->t1);

    if (kx <= KT) {
        int c = kx + hx;
        uint64_t dD0 = (uint64_t)D0, D02 = dD0 * dD0, D03 = D02 * dD0;
        uint64_t tn = (c > D0) ? (uint64_t)D0 : 1ULL;
        uint64_t td = (c > D0) ? (uint64_t)c  : 1ULL;
        uint64_t td2 = td * td, td3 = td2 * td, td4 = td2 * td2;
        uint64_t tn2 = tn * tn, tn3 = tn2 * tn, tn4 = tn2 * tn2;

        /* Precompute exact positive integer BB[0..4] = 12 * td^4 * B_j */
        uint64_t b0 = D03, b1 = 3 * D02, b2 = 3 * dD0, b3 = 1;
        uint64_t c0 = 12 * td4, c1 = 3 * td3 * tn, c2 = 2 * td2 * tn2, c3 = 3 * td * tn3, c4 = 12 * tn4;
        mpz_set_ui(W->BB[0], c0); mpz_mul_ui(W->BB[0], W->BB[0], b0);
        mpz_set_ui(W->BB[1], c0); mpz_mul_ui(W->BB[1], W->BB[1], b0);
        mpz_set_ui(W->t1, c1);    mpz_mul_ui(W->t1, W->t1, b1); mpz_add(W->BB[1], W->BB[1], W->t1);
        mpz_set_ui(W->BB[2], c0); mpz_mul_ui(W->BB[2], W->BB[2], b0);
        mpz_set_ui(W->t1, 2*c1);  mpz_mul_ui(W->t1, W->t1, b1); mpz_add(W->BB[2], W->BB[2], W->t1);
        mpz_set_ui(W->t1, c2);    mpz_mul_ui(W->t1, W->t1, b2); mpz_add(W->BB[2], W->BB[2], W->t1);
        mpz_set_ui(W->BB[3], c0); mpz_mul_ui(W->BB[3], W->BB[3], b0);
        mpz_set_ui(W->t1, 3*c1);  mpz_mul_ui(W->t1, W->t1, b1); mpz_add(W->BB[3], W->BB[3], W->t1);
        mpz_set_ui(W->t1, 3*c2);  mpz_mul_ui(W->t1, W->t1, b2); mpz_add(W->BB[3], W->BB[3], W->t1);
        mpz_set_ui(W->t1, c3);    mpz_mul_ui(W->t1, W->t1, b3); mpz_add(W->BB[3], W->BB[3], W->t1);
        mpz_set_ui(W->BB[4], c0); mpz_mul_ui(W->BB[4], W->BB[4], b0);
        mpz_set_ui(W->t1, 4*c1);  mpz_mul_ui(W->t1, W->t1, b1); mpz_add(W->BB[4], W->BB[4], W->t1);
        mpz_set_ui(W->t1, 6*c2);  mpz_mul_ui(W->t1, W->t1, b2); mpz_add(W->BB[4], W->BB[4], W->t1);
        mpz_set_ui(W->t1, 4*c3);  mpz_mul_ui(W->t1, W->t1, b3); mpz_add(W->BB[4], W->BB[4], W->t1);

        int64_t rin_ns[2] = { rlo_n, rhi_n }, rin_ds[2] = { rlo_d, rhi_d };
        int nrin = (rlo_n * rhi_d == rhi_n * rlo_d) ? 1 : 2;
        int rext_2s[2] = { 0, ext_half };
        int nrext = (ext_half > 0) ? 2 : 1;

        /* Precompute (D0 * a0), (D0 * a1), (D0 * e0), (D0 * e1) BEFORE multiplying by x_den! */
        /* D0*a0 = D0*F3b + 3*c*((F1b-F1)-(F2b-F2)) + (d0 - 3*d1 + 3*d2) [if n_nbr_in>0] */
        mpz_t da0, da1;
        mpz_init(da0); mpz_init(da1);
        /* db_scaled = 3*c*((F1b - F1) - (F2b - F2)) */
        mpz_set_si(W->t1, F1b); mpz_set_si(W->t2, F1); mpz_sub(W->t1, W->t1, W->t2);
        mpz_set_si(W->t2, F2b); mpz_set_si(W->t3, F2); mpz_sub(W->t2, W->t2, W->t3);
        mpz_sub(W->t1, W->t1, W->t2);
        mpz_mul_ui(W->t1, W->t1, (uint64_t)(3 * c)); /* 3*c*(...) */
        /* nd0_scaled = d0 - 3*d1 + 3*d2 (if n_nbr_in > 0) */
        mpz_set_ui(W->t2, 0);
        if (n_nbr_in > 0) {
            mpz_set_si(W->t2, d0_nbr);
            mpz_set_si(W->t3, d1_nbr); mpz_submul_ui(W->t2, W->t3, 3);
            mpz_set_si(W->t3, d2_nbr); mpz_addmul_ui(W->t2, W->t3, 3);
        }
        /* da0 = D0*F3b + W->t1 + W->t2 */
        mpz_set_si(da0, F3b); mpz_mul_ui(da0, da0, dD0);
        mpz_add(da0, da0, W->t1); mpz_add(da0, da0, W->t2);
        /* da1 = -c*F3b - W->t1 + (n_nbr_in>0 ? F3_up_sum : 0) - W->t2 */
        mpz_set_si(da1, F3b); mpz_mul_si(da1, da1, -(int64_t)c);
        mpz_sub(da1, da1, W->t1); mpz_sub(da1, da1, W->t2);
        if (n_nbr_in > 0) { mpz_set_si(W->t3, F3_up_sum); mpz_add(da1, da1, W->t3); }
        /* Subtract (D0 * e0)*F3 and (D0 * e1)*F3 directly into da0, da1!
         * Since P(t) = N(t) - D(t)*F3 and an[k] - ad[k]*F3 is linear in (a0 - e0*F3, a1 - e1*F3)! */
        /* D0*(a0 - e0*F3) = da0 - D0*F3; D0*(a1 - e1*F3) = da1 - (-c + n_nbr_in)*F3 */
        mpz_set_si(W->t3, F3); mpz_submul_ui(da0, W->t3, dD0);
        mpz_set_si(W->t3, F3); mpz_set_si(W->t4, (int64_t)(n_nbr_in - c)); mpz_submul(da1, W->t3, W->t4);

        /* F2_m_F3 = F2 - F3 */
        mpz_t f23; mpz_init(f23);
        mpz_set_si(f23, F2); mpz_set_si(W->t3, F3); mpz_sub(f23, f23, W->t3);

        int have_best = 0;
        for (int ri = 0; ri < nrin; ri++) for (int re = 0; re < nrext; re++) {
            /* x = (c + 1) - rin_ns[ri]/rin_ds[ri] - rext_2s[re]/2 = x_num / x_den */
            int64_t x_den = 2 * rin_ds[ri];
            int64_t x_num = (int64_t)(c + 1) * x_den - 2 * rin_ns[ri] - (int64_t)rext_2s[re] * rin_ds[ri];
            uint64_t M_u = dD0 * (uint64_t)x_den; /* positive integer multiplier */

            /* A0 = x_den * (da0 - D0*F3), A1 = x_den * (da1 - (n_in-c)*F3) */
            mpz_mul_ui(W->A0, da0, (uint64_t)x_den);
            mpz_mul_ui(W->A1, da1, (uint64_t)x_den);

            /* P[k] = M * (n_k - d_k * F3) for k = 0..4!
             * Base from (A0, A1):
             *   P[0] = D0^3 * A0
             *   P[1] = D0^2 * (D0 * A1 + 3 * A0)
             *   P[2] = D0   * (3 * D0 * A1 + 3 * A0)
             *   P[3] = 3 * D0 * A1 + A0
             *   P[4] = A1
             */
            mpz_mul_ui(W->P[0], W->A0, D03);
            mpz_mul_ui(W->P[1], W->A1, dD0); mpz_addmul_ui(W->P[1], W->A0, 3); mpz_mul_ui(W->P[1], W->P[1], D02);
            mpz_mul_ui(W->P[2], W->A1, 3 * dD0); mpz_addmul_ui(W->P[2], W->A0, 3); mpz_mul_ui(W->P[2], W->P[2], dD0);
            mpz_mul_ui(W->P[3], W->A1, 3 * dD0); mpz_add(W->P[3], W->P[3], W->A0);
            mpz_set(W->P[4], W->A1);

            /* Add F2 - F3 terms multiplied by M = D0 * x_den (so M*1 = D0*x_den, M*x = D0*x_num):
             *   k=0: +3*D0^3 * M * f23
             *   k=1: +3*D0^2 * (2*M - M*x) * f23
             *   k=2: +D0 * (4*M - 3*M*x) * f23
             *   k=3: +(M - M*x) * f23
             */
            int64_t M_s = (int64_t)M_u, Mx = (int64_t)dD0 * x_num;
            mpz_set_si(W->t1, 3 * (int64_t)D03); mpz_mul_si(W->t1, W->t1, M_s); mpz_addmul(W->P[0], f23, W->t1);
            mpz_set_si(W->t1, 2 * M_s - Mx);     mpz_mul_ui(W->t1, W->t1, 3 * D02); mpz_addmul(W->P[1], f23, W->t1);
            mpz_set_si(W->t1, 4 * M_s - 3 * Mx); mpz_mul_ui(W->t1, W->t1, dD0);     mpz_addmul(W->P[2], f23, W->t1);
            mpz_set_si(W->t1, M_s - Mx);                                             mpz_addmul(W->P[3], f23, W->t1);

            /* Add d1z terms multiplied by M = D0 * x_den (note M * x / D0 = x_num):
             *   k=0: (3*D0^2*Mx + 3*D0^2*M + 3*D0*Mx)
             *   k=1: (-3*D0^2*Mx - 3*D0^2*M + 6*D0*M + 6*Mx)
             *   k=2: (-3*D0*Mx - 5*D0*M - 6*Mx + 3*M + 3*x_num)
             *   k=3: (-Mx - 2*M - 3*x_num)
             */
            mpz_set_si(W->t1, 3 * (int64_t)D02 * (Mx + M_s) + 3 * (int64_t)dD0 * Mx); mpz_addmul(W->P[0], d1z, W->t1);
            mpz_set_si(W->t1, -3 * (int64_t)D02 * (Mx + M_s) + 6 * (int64_t)dD0 * M_s + 6 * Mx); mpz_addmul(W->P[1], d1z, W->t1);
            mpz_set_si(W->t1, -3 * (int64_t)dD0 * Mx - 5 * (int64_t)dD0 * M_s - 6 * Mx + 3 * M_s + 3 * x_num); mpz_addmul(W->P[2], d1z, W->t1);
            mpz_set_si(W->t1, -Mx - 2 * M_s - 3 * x_num); mpz_addmul(W->P[3], d1z, W->t1);

            /* Add d2z terms multiplied by M:
             *   k=0: (-6*D0^2*Mx - 3*D0^2*M)
             *   k=1: (6*D0^2*Mx + 6*D0^2*M - 9*D0*Mx - 6*D0*M)
             *   k=2: (6*D0*Mx + 10*D0*M - 3*Mx - 3*M)
             *   k=3: (2*Mx + 4*M)
             */
            mpz_set_si(W->t1, -6 * (int64_t)D02 * Mx - 3 * (int64_t)D02 * M_s); mpz_addmul(W->P[0], d2z, W->t1);
            mpz_set_si(W->t1, 6 * (int64_t)D02 * (Mx + M_s) - 9 * (int64_t)dD0 * Mx - 6 * (int64_t)dD0 * M_s); mpz_addmul(W->P[1], d2z, W->t1);
            mpz_set_si(W->t1, 6 * (int64_t)dD0 * Mx + 10 * (int64_t)dD0 * M_s - 3 * Mx - 3 * M_s); mpz_addmul(W->P[2], d2z, W->t1);
            mpz_set_si(W->t1, 2 * Mx + 4 * M_s); mpz_addmul(W->P[3], d2z, W->t1);

            /* Convert P[0..4] to Bernstein-4 BP[0..4] = 12 * td^4 * B_j(P) */
            mpz_mul_ui(W->BP[0], W->P[0], c0);
            mpz_mul_ui(W->BP[1], W->P[0], c0); mpz_addmul_ui(W->BP[1], W->P[1], c1);
            mpz_mul_ui(W->BP[2], W->P[0], c0); mpz_addmul_ui(W->BP[2], W->P[1], 2 * c1); mpz_addmul_ui(W->BP[2], W->P[2], c2);
            mpz_mul_ui(W->BP[3], W->P[0], c0); mpz_addmul_ui(W->BP[3], W->P[1], 3 * c1); mpz_addmul_ui(W->BP[3], W->P[2], 3 * c2); mpz_addmul_ui(W->BP[3], W->P[3], c3);
            mpz_mul_ui(W->BP[4], W->P[0], c0); mpz_addmul_ui(W->BP[4], W->P[1], 4 * c1); mpz_addmul_ui(W->BP[4], W->P[2], 6 * c2); mpz_addmul_ui(W->BP[4], W->P[3], 4 * c3); mpz_addmul_ui(W->BP[4], W->P[4], c4);

            /* Optional shift = -(rext_2s[re]/2 - Rbar_q/S_r) * Y_ref = ((2*Rbar_q - rext_2s[re]*S_r) * Y_ref) / (2*S_r) */
            int has_shift = (xfer_ext && Y_ref != 0);
            if (has_shift) {
                mpz_set_si(W->sh_n, Rbar_q); mpz_mul_ui(W->sh_n, W->sh_n, 2);
                mpz_set_si(W->t1, S_r);      mpz_submul_ui(W->sh_n, W->t1, (uint64_t)rext_2s[re]);
                mpz_mul_si(W->sh_n, W->sh_n, Y_ref);
                mpz_set_ui(W->sh_d, (uint64_t)(2 * S_r));
            }

            for (int j = 0; j < 5; j++) {
                mpz_set(W->cur.num, W->BP[j]);
                mpz_mul_ui(W->cur.den, W->BB[j], M_u);
                if (has_shift) zf_add_z(&W->cur, W->sh_n, W->sh_d, W);
                if (!have_best) {
                    mpz_set(W->best.num, W->cur.num); mpz_set(W->best.den, W->cur.den); have_best = 1;
                } else {
                    int cmp = zf_cmp(&W->cur, &W->best, W);
                    if ((mode == 1 && cmp > 0) || (mode == 2 && cmp < 0)) {
                        mpz_set(W->best.num, W->cur.num); mpz_set(W->best.den, W->cur.den);
                    }
                }
            }
        }
        zf_add_z(acc, W->best.num, W->best.den, W);
        mpz_clear(da0); mpz_clear(da1); mpz_clear(f23);
    } else {
        /* overflow k > KT: degree-2 Bernstein hull on [0, 1] vs 0 */
        uint64_t dD0 = (uint64_t)D0, D02 = dD0 * dD0;
        mpz_t f23; mpz_init(f23);
        mpz_set_si(f23, F2); mpz_set_si(W->t1, F3); mpz_sub(f23, f23, W->t1);

        /* P0 = 3*D0^2*f23 + 3*D0*(d1 - d2) */
        mpz_mul_ui(W->P[0], f23, 3 * D02);
        mpz_sub(W->t1, d1z, d2z); mpz_addmul_ui(W->P[0], W->t1, 3 * dD0);
        /* P1 = 3*D0*f23 + (3 - 3*D0)*d1 + (6*D0 - 3)*d2 */
        mpz_mul_ui(W->P[1], f23, 3 * dD0);
        mpz_set_si(W->t1, 3 - 3 * (int64_t)dD0); mpz_addmul(W->P[1], d1z, W->t1);
        mpz_set_si(W->t1, 6 * (int64_t)dD0 - 3); mpz_addmul(W->P[1], d2z, W->t1);
        /* P2 = f23 - 2*d1 + 4*d2 */
        mpz_set(W->P[2], f23);
        mpz_submul_ui(W->P[2], d1z, 2);
        mpz_addmul_ui(W->P[2], d2z, 4);

        /* 2 * BP2[0..2] and 2 * BB2[0..2] */
        mpz_mul_ui(W->BP[0], W->P[0], 2);
        mpz_mul_ui(W->BP[1], W->P[0], 2); mpz_add(W->BP[1], W->BP[1], W->P[1]);
        mpz_mul_ui(W->BP[2], W->P[0], 2); mpz_addmul_ui(W->BP[2], W->P[1], 2); mpz_addmul_ui(W->BP[2], W->P[2], 2);
        uint64_t bb2[3] = { 2 * D02, 2 * D02 + 2 * dD0, 2 * D02 + 4 * dD0 + 2 };

        zf_set_si(&W->best, 0, 1); /* alpha = 0 corner */
        for (int j = 0; j < 3; j++) {
            mpz_set(W->cur.num, W->BP[j]); mpz_set_ui(W->cur.den, bb2[j]);
            int cmp = zf_cmp(&W->cur, &W->best, W);
            if ((mode == 1 && cmp > 0) || (mode == 2 && cmp < 0)) {
                mpz_set(W->best.num, W->cur.num); mpz_set(W->best.den, W->cur.den);
            }
        }
        zf_add_z(acc, W->best.num, W->best.den, W);
        if (n_nbr_in > 0) {
            /* ext(d0 - 3*d1 + 3*d2, F3_up_sum - n_nbr_in*F3) / D0 */
            mpz_set_si(W->t1, d0_nbr);
            mpz_set_si(W->t3, d1_nbr); mpz_submul_ui(W->t1, W->t3, 3);
            mpz_set_si(W->t3, d2_nbr); mpz_addmul_ui(W->t1, W->t3, 3);
            mpz_set_si(W->t2, F3_up_sum);
            mpz_set_si(W->t3, F3);     mpz_submul_ui(W->t2, W->t3, (uint64_t)n_nbr_in);
            int cmp = mpz_cmp(W->t2, W->t1);
            mpz_set_ui(W->t4, dD0);
            if ((mode == 1 && cmp > 0) || (mode == 2 && cmp < 0)) zf_add_z(acc, W->t2, W->t4, W);
            else                                                   zf_add_z(acc, W->t1, W->t4, W);
        }
        mpz_clear(f23);
    }
    mpz_clear(d1z); mpz_clear(d2z);
}

/* Exact rational evaluation of Phi(T) in mpq_t for built-in concrete tree verification */
typedef struct {
    const B3FLayout *C;
    const CertB3FHeader *H;
    const L5 *X;
    const int64_t *q5;
    const int64_t *rpar, *rch, *trpar, *trch, *c2_f, *c3_f, *Yp_ref, *Yw_ref;
    const int64_t *f, *g, *hV, *hC, *qP4, *qYt, *qYb, *qS3;
} CertData;
#include "l5_verify_rows.h"

/* Grow by core edges, in increasing physical vertex order. Parents precede
 * children, so each connected subset has exactly one growth path. */
static void tree_l5(const CertData *d, int n, const int *par, const int *deg,
                    const int *L, const int *k, const int *h, int vertices[5],
                    int len, long long *sum, mpz_t total) {
    if (len == 5) {
        L5Tree t = {.n=5,.p={-1}};
        for (int i = 1; i < 5; ++i) {
            int j = 0;
            while (vertices[j] != par[vertices[i]]) ++j;
            t.p[i] = j;
        }
        for (int i = 0; i < 5; ++i) {
            int u = vertices[i];
            if (deg[u] < d->X->dc) t.g[i] = b3f_inner(d->C,L[u],k[u],h[u]);
            else {
                if (l5_degree(&t,i) != 1) return;
                t.g[i] = l5_big(i != 0,d->X->buckets == 2 && deg[u] >= d->C->D0);
            }
        }
        int64_t v = l5_qint(d->X,d->q5,&t);
        /* Keep pending additions within native range; spill exactly to GMP. */
        if (*sum > INT64_MAX/2 || *sum < INT64_MIN/2) {
            if (*sum >= 0) mpz_add_ui(total,total,(unsigned long)*sum);
            else mpz_sub_ui(total,total,(unsigned long)-*sum);
            *sum = 0;
        }
        *sum += v;
        return;
    }
    for (int u = vertices[len-1]+1; u < n; ++u) if (deg[u] >= 2) {
        for (int j = 0; j < len; ++j) if (par[u] == vertices[j]) {
            vertices[len] = u;
            tree_l5(d,n,par,deg,L,k,h,vertices,len+1,sum,total);
            break;
        }
    }
}

static void b3f_embed_q(const B3FLayout *C, int L, int k, int h, int ids[4], mpq_t wts[4], int *nw) {
    int D = L + k + h;
    if (D < C->D0) {
        ids[0] = b3f_inner(C, L, k, h);
        mpq_set_ui(wts[0], 1, 1);
        *nw = 1;
        return;
    }
    mpq_t t, s, t2, s2;
    mpq_init(t); mpq_init(s); mpq_init(t2); mpq_init(s2);
    mpq_set_ui(t, (uint64_t)C->D0, (uint64_t)D); mpq_canonicalize(t);
    mpq_set_ui(s, (uint64_t)(D - C->D0), (uint64_t)D); mpq_canonicalize(s);
    mpq_mul(t2, t, t); mpq_mul(s2, s, s);
    ids[0] = C->INF_ID;       mpq_mul(wts[0], s2, s);
    ids[1] = b3f_t1(C, k, h); mpq_mul(wts[1], t, s2);  mpz_mul_ui(mpq_numref(wts[1]), mpq_numref(wts[1]), 3); mpq_canonicalize(wts[1]);
    ids[2] = b3f_t2(C, k, h); mpq_mul(wts[2], t2, s);  mpz_mul_ui(mpq_numref(wts[2]), mpq_numref(wts[2]), 3); mpq_canonicalize(wts[2]);
    ids[3] = b3f_t3(C, k, h); mpq_mul(wts[3], t2, t);
    *nw = 4;
    mpq_clear(t); mpq_clear(s); mpq_clear(t2); mpq_clear(s2);
}

static void add_scaled_q(mpq_t tot, const mpq_t w, int64_t val, mpq_t tmp) {
    if (val == 0) return;
    mpq_set_si(tmp, val, 1);
    mpq_mul(tmp, tmp, w);
    mpq_add(tot, tot, tmp);
}

/* Evaluates S * Phi(T) in exact mpq_t */
static void eval_Sphi_q(const CertData *D, int n, const int *par, const int *deg, mpq_t out) {
    const B3FLayout *C = D->C;
    int *Lv = (int*)b3f_xcalloc(n, sizeof(int)), *kv = (int*)b3f_xcalloc(n, sizeof(int)), *hv = (int*)b3f_xcalloc(n, sizeof(int)), *isC = (int*)b3f_xcalloc(n, sizeof(int));
    int *ck_cnt = (int*)b3f_xcalloc(n, sizeof(int));
    int **ck = (int**)b3f_xcalloc(n, sizeof(int*));
    for (int u = 0; u < n; u++) if (deg[u] >= 2) { isC[u] = 1; hv[u] = (par[u] >= 0) ? 1 : 0; ck[u] = (int*)b3f_xmalloc((size_t)deg[u] * sizeof(int)); }
    for (int v = 1; v < n; v++) {
        int p = par[v];
        if (deg[v] == 1) Lv[p]++;
        else ck[p][ck_cnt[p]++] = v;
    }
    for (int u = 0; u < n; u++) kv[u] = ck_cnt[u];

    mpq_set_ui(out, 0, 1);
    mpq_t tmp, w12, w123, w_a[4], w_b[4], w_c[4], w_d[4];
    mpq_init(tmp); mpq_init(w12); mpq_init(w123);
    for (int i = 0; i < 4; i++) { mpq_init(w_a[i]); mpq_init(w_b[i]); mpq_init(w_c[i]); mpq_init(w_d[i]); }
    int id_a[4], id_b[4], id_c[4], id_d[4], na, nb, nc, nd;
    int N = C->NG;

    /* f */
    for (int u = 0; u < n; u++) if (isC[u]) {
        mpq_set_si(tmp, D->f[b3f_vidx(C, Lv[u], kv[u], hv[u])], 1);
        mpq_add(out, out, tmp);
    }
    /* g */
    for (int w = 1; w < n; w++) if (isC[w]) {
        int p = par[w];
        b3f_embed_q(C, Lv[p], kv[p], hv[p], id_a, w_a, &na);
        b3f_embed_q(C, Lv[w], kv[w], hv[w], id_b, w_b, &nb);
        for (int i = 0; i < na; i++) for (int j = 0; j < nb; j++) {
            mpq_mul(w12, w_a[i], w_b[j]);
            add_scaled_q(out, w12, D->g[id_a[i] * N + id_b[j]], tmp);
        }
    }
    /* hV */
    for (int u = 1; u < n; u++) if (isC[u] && kv[u] >= 1 && kv[u] <= C->KT) {
        int p = par[u];
        b3f_embed_q(C, Lv[u], kv[u], 1, id_b, w_b, &nb);
        b3f_embed_q(C, Lv[p], kv[p], hv[p], id_a, w_a, &na);
        for (int j = 0; j < ck_cnt[u]; j++) {
            int w = ck[u][j];
            b3f_embed_q(C, Lv[w], kv[w], 1, id_c, w_c, &nc);
            for (int iu = 0; iu < nb; iu++) {
                int ua = C->hu_of_g[id_b[iu]]; if (ua < 0) continue;
                for (int ia = 0; ia < na; ia++) for (int ic = 0; ic < nc; ic++) {
                    mpq_mul(w12, w_b[iu], w_a[ia]); mpq_mul(w123, w12, w_c[ic]);
                    add_scaled_q(out, w123, D->hV[b3f_hv_off(C, id_a[ia], ua, id_c[ic])], tmp);
                }
            }
        }
    }
    /* hC */
    for (int u = 0; u < n; u++) if (isC[u] && kv[u] >= 2 && kv[u] <= C->KT) {
        b3f_embed_q(C, Lv[u], kv[u], hv[u], id_a, w_a, &na);
        for (int j1 = 0; j1 < ck_cnt[u]; j1++) for (int j2 = j1 + 1; j2 < ck_cnt[u]; j2++) {
            int w1 = ck[u][j1], w2 = ck[u][j2];
            b3f_embed_q(C, Lv[w1], kv[w1], 1, id_b, w_b, &nb);
            b3f_embed_q(C, Lv[w2], kv[w2], 1, id_c, w_c, &nc);
            for (int iu = 0; iu < na; iu++) {
                int ca = C->cu_of_g[id_a[iu]]; if (ca < 0) continue;
                for (int ib = 0; ib < nb; ib++) for (int ic = 0; ic < nc; ic++) {
                    mpq_mul(w12, w_a[iu], w_b[ib]); mpq_mul(w123, w12, w_c[ic]);
                    add_scaled_q(out, w123, D->hC[b3f_hc_off(C, ca, id_b[ib], id_c[ic])], tmp);
                }
            }
        }
    }
    /* P4 */
    for (int c = 1; c < n; c++) if (isC[c] && kv[c] >= 1 && kv[c] <= C->KP4 && deg[c] < C->DP4) {
        int b = par[c];
        if (hv[b] == 1 && kv[b] >= 1 && kv[b] <= C->KP4 && deg[b] < C->DP4) {
            int b4 = C->p4_lut[kv[b]][Lv[b]], c4 = C->p4_lut[kv[c]][Lv[c]], a = par[b];
            b3f_embed_q(C, Lv[a], kv[a], hv[a], id_a, w_a, &na);
            for (int j = 0; j < ck_cnt[c]; j++) {
                int d = ck[c][j];
                b3f_embed_q(C, Lv[d], kv[d], 1, id_d, w_d, &nd);
                for (int ia = 0; ia < na; ia++) for (int id = 0; id < nd; id++) {
                    mpq_mul(w12, w_a[ia], w_d[id]);
                    add_scaled_q(out, w12, D->qP4[b3f_p4_off(C, id_a[ia], b4, c4, id_d[id])], tmp);
                }
            }
        }
    }
    /* Ytop */
    for (int u = 1; u < n; u++) if (isC[u] && kv[u] >= 2 && kv[u] <= C->KYT && deg[u] < C->DYT) {
        int u4 = C->yt_lut[kv[u]][Lv[u]], p = par[u];
        b3f_embed_q(C, Lv[p], kv[p], hv[p], id_a, w_a, &na);
        for (int j1 = 0; j1 < ck_cnt[u]; j1++) for (int j2 = j1 + 1; j2 < ck_cnt[u]; j2++) {
            int w1 = ck[u][j1], w2 = ck[u][j2];
            b3f_embed_q(C, Lv[w1], kv[w1], 1, id_b, w_b, &nb);
            b3f_embed_q(C, Lv[w2], kv[w2], 1, id_c, w_c, &nc);
            for (int ia = 0; ia < na; ia++) for (int ib = 0; ib < nb; ib++) for (int ic = 0; ic < nc; ic++) {
                mpq_mul(w12, w_a[ia], w_b[ib]); mpq_mul(w123, w12, w_c[ic]);
                add_scaled_q(out, w123, D->qYt[b3f_yt_off(C, id_a[ia], u4, id_b[ib], id_c[ic])], tmp);
            }
        }
    }
    /* Ybot */
    for (int u = 0; u < n; u++) if (isC[u] && kv[u] >= 2 && kv[u] <= C->KYB && deg[u] < C->DYB) {
        int u4 = C->ybu_lut[hv[u]][kv[u]][Lv[u]];
        for (int j2 = 0; j2 < ck_cnt[u]; j2++) {
            int w2 = ck[u][j2];
            if (!(kv[w2] >= 1 && kv[w2] <= C->KYB && deg[w2] < C->DYB)) continue;
            int w2_4 = C->ybw_lut[kv[w2]][Lv[w2]];
            for (int j1 = 0; j1 < ck_cnt[u]; j1++) if (j1 != j2) {
                int w1 = ck[u][j1];
                b3f_embed_q(C, Lv[w1], kv[w1], 1, id_b, w_b, &nb);
                for (int jc = 0; jc < ck_cnt[w2]; jc++) {
                    int c = ck[w2][jc];
                    b3f_embed_q(C, Lv[c], kv[c], 1, id_c, w_c, &nc);
                    for (int ib = 0; ib < nb; ib++) for (int ic = 0; ic < nc; ic++) {
                        mpq_mul(w12, w_b[ib], w_c[ic]);
                        add_scaled_q(out, w12, D->qYb[b3f_yb_off(C, u4, w2_4, id_b[ib], id_c[ic])], tmp);
                    }
                }
            }
        }
    }
    /* S3 */
    for (int u = 0; u < n; u++) if (isC[u] && kv[u] >= 3 && kv[u] <= C->KS3 && deg[u] < C->DS3) {
        int u4 = C->s3_lut[hv[u]][kv[u]][Lv[u]];
        for (int j1 = 0; j1 < ck_cnt[u]; j1++) for (int j2 = j1 + 1; j2 < ck_cnt[u]; j2++) for (int j3 = j2 + 1; j3 < ck_cnt[u]; j3++) {
            int w1 = ck[u][j1], w2 = ck[u][j2], w3 = ck[u][j3];
            b3f_embed_q(C, Lv[w1], kv[w1], 1, id_b, w_b, &nb);
            b3f_embed_q(C, Lv[w2], kv[w2], 1, id_c, w_c, &nc);
            b3f_embed_q(C, Lv[w3], kv[w3], 1, id_d, w_d, &nd);
            for (int ib = 0; ib < nb; ib++) for (int ic = 0; ic < nc; ic++) for (int id = 0; id < nd; id++) {
                mpq_mul(w12, w_b[ib], w_c[ic]); mpq_mul(w123, w12, w_d[id]);
                add_scaled_q(out, w123, D->qS3[b3f_s3_off(C, u4, id_b[ib], id_c[ic], id_d[id])], tmp);
            }
        }
    }

    long long sum = 0;
    mpq_set_ui(tmp,0,1);
    for (int u = 0; u < n; ++u) if (isC[u]) {
        int vertices[5] = {u};
        tree_l5(D,n,par,deg,Lv,kv,hv,vertices,1,&sum,mpq_numref(tmp));
    }
    if (sum >= 0) mpz_add_ui(mpq_numref(tmp),mpq_numref(tmp),(unsigned long)sum);
    else mpz_sub_ui(mpq_numref(tmp),mpq_numref(tmp),(unsigned long)-sum);
    mpq_add(out,out,tmp);
    mpq_clear(tmp); mpq_clear(w12); mpq_clear(w123);
    for (int i = 0; i < 4; i++) { mpq_clear(w_a[i]); mpq_clear(w_b[i]); mpq_clear(w_c[i]); mpq_clear(w_d[i]); }
    for (int u = 0; u < n; u++) free(ck[u]);
    free(Lv); free(kv); free(hv); free(isC); free(ck_cnt); free(ck);
}

static int check_concrete_tree_q(const CertData *D, int n, const int *par_in,
                                 const mpq_t lower, const mpq_t upper) {
    int *par = (int*)b3f_xmalloc((size_t)(n + 1) * sizeof(int));
    int *deg = (int*)b3f_xcalloc((size_t)(n + 1), sizeof(int));
    memcpy(par, par_in, (size_t)n * sizeof(int));
    for (int v = 1; v < n; v++) { deg[v]++; deg[par[v]]++; }
    mpq_t phi0, phi1, lphi, rate, diff, rhs, tmp;
    mpq_init(phi0); mpq_init(phi1); mpq_init(lphi); mpq_init(rate); mpq_init(diff); mpq_init(rhs); mpq_init(tmp);
    eval_Sphi_q(D, n, par, deg, phi0);
    int core = 0;
    for (int w = 0; w < n; ++w) core += deg[w] >= 2;
    mpq_set_ui(rate,core,1);
    mpq_mul(tmp,lower,rate);
    int bounds_ok = mpq_cmp(tmp,phi0) <= 0;
    mpq_mul(tmp,upper,rate);
    bounds_ok &= mpq_cmp(tmp,phi0) >= 0;
    mpq_set_ui(lphi, 0, 1);
    for (int w = 0; w < n; w++) {
        mpq_set_ui(rate, 0, 1);
        if (par[w] >= 0) { mpq_set_ui(tmp, 1, (uint64_t)deg[par[w]]); mpq_add(rate, rate, tmp); }
        for (int v = 1; v < n; v++) if (par[v] == w) { mpq_set_ui(tmp, 1, (uint64_t)deg[v]); mpq_add(rate, rate, tmp); }
        par[n] = w; deg[n] = 1; deg[w]++;
        eval_Sphi_q(D, n + 1, par, deg, phi1);
        deg[w]--;
        mpq_sub(diff, phi1, phi0);
        mpq_mul(tmp, rate, diff);
        mpq_add(lphi, lphi, tmp);
    }
    mpq_set_si(rhs,D->H->P,(uint64_t)D->H->Q); mpq_canonicalize(rhs);
    mpq_mul(rhs, rhs, phi0);
    int cmp = mpq_cmp(lphi, rhs);
    int ok = (D->H->mode == 1) ? (cmp <= 0) : (cmp >= 0);
    mpq_clear(phi0); mpq_clear(phi1); mpq_clear(lphi); mpq_clear(rate); mpq_clear(diff); mpq_clear(rhs); mpq_clear(tmp);
    free(par); free(deg);
    return ok && bounds_ok;
}

static int64_t *read_i64_array(FILE *fp, size_t n, int64_t max_abs, const char *name) {
    int64_t *arr = l5_array(n,sizeof *arr);
    if (!arr) { fprintf(stderr, "FAIL: OOM allocating %zu int64s for %s\n", n, name); exit(1); }
    if (n > 0 && fread(arr, sizeof(int64_t), n, fp) != n) {
        fprintf(stderr, "FAIL: unexpected EOF reading %zu int64s for %s\n", n, name);
        free(arr); exit(1);
    }
    for (size_t i = 0; i < n; i++) {
        if (arr[i] < -max_abs || arr[i] > max_abs) {
            fprintf(stderr, "FAIL: %s[%zu] = %lld exceeds verified native bound %lld\n",
                    name, i, (long long)arr[i], (long long)max_abs);
            free(arr); exit(1);
        }
    }
    return arr;
}

/* Exact rational comparison bounds c0 C(T) <= Phi(T) <= C0 C(T) via exact vertex charging */
static int verify_comparison_bounds(const CertData *D, mpq_t lower, mpq_t upper) {
    const B3FLayout *C = D->C;
    const CertB3FHeader *H = D->H;
    int N = C->NG, D0 = C->D0, KT = C->KT;

    /* 1. For each child state w in chs: min/max of g[p,w], hV[p,ua,w], P4[a,b4,c4,w], and (k_u4-1)*Yb[u4,w2_4,w1,w] */
    int64_t *g_min = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t)), *g_max = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t));
    int64_t *hv_min = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t)), *hv_max = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t));
    int64_t *p4_min = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t)), *p4_max = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t));
    int64_t *yb_min = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t)), *yb_max = (int64_t*)b3f_xmalloc(C->nchs * sizeof(int64_t));

    for (int iw = 0; iw < C->nchs; iw++) {
        int w = C->chs[iw];
        int64_t mn = D->g[C->pars[0] * N + w], mx = mn;
        for (int ip = 1; ip < C->npars; ip++) {
            int64_t v = D->g[C->pars[ip] * N + w];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        g_min[iw] = mn; g_max[iw] = mx;
        hv_min[iw] = 0; hv_max[iw] = 0;
        for (int ip = 0; ip < C->npars; ip++) for (int ua = 0; ua < C->NHU; ua++) {
            int64_t v = D->hV[b3f_hv_off(C, C->pars[ip], ua, w)];
            if (v < hv_min[iw]) hv_min[iw] = v;
            if (v > hv_max[iw]) hv_max[iw] = v;
        }
        p4_min[iw] = 0; p4_max[iw] = 0;
        for (int ip = 0; ip < C->npars; ip++) for (int b4 = 0; b4 < C->NP4; b4++) for (int c4 = 0; c4 < C->NP4; c4++) {
            int64_t v = D->qP4[b3f_p4_off(C, C->pars[ip], b4, c4, w)];
            if (v < p4_min[iw]) p4_min[iw] = v;
            if (v > p4_max[iw]) p4_max[iw] = v;
        }
        yb_min[iw] = 0; yb_max[iw] = 0;
        for (int u4 = 0; u4 < C->NYBU; u4++) {
            int64_t mult = (int64_t)(C->ybu_K[u4] - 1);
            for (int w2_4 = 0; w2_4 < C->NYBW; w2_4++) for (int i1 = 0; i1 < C->nchs; i1++) {
                int64_t v = mult * D->qYb[b3f_yb_off(C, u4, w2_4, C->chs[i1], w)];
                if (v < yb_min[iw]) yb_min[iw] = v;
                if (v > yb_max[iw]) yb_max[iw] = v;
            }
        }
    }

    /* 2. Apex bounds: hC[ca] per pair, Yt[u4] per pair, S3[u4] per triple */
    int64_t *hc_min = (int64_t*)b3f_xmalloc(C->NCU * sizeof(int64_t)), *hc_max = (int64_t*)b3f_xmalloc(C->NCU * sizeof(int64_t));
    for (int ca = 0; ca < C->NCU; ca++) {
        int64_t mn = D->hC[b3f_hc_off(C, ca, C->chs[0], C->chs[0])], mx = mn;
        for (int i1 = 0; i1 < C->nchs; i1++) for (int i2 = i1; i2 < C->nchs; i2++) {
            int64_t v = D->hC[b3f_hc_off(C, ca, C->chs[i1], C->chs[i2])];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        hc_min[ca] = mn; hc_max[ca] = mx;
    }
    int64_t *yt_min = (int64_t*)b3f_xmalloc(C->NYT * sizeof(int64_t)), *yt_max = (int64_t*)b3f_xmalloc(C->NYT * sizeof(int64_t));
    for (int u4 = 0; u4 < C->NYT; u4++) {
        int64_t mn = D->qYt[b3f_yt_off(C, C->pars[0], u4, C->chs[0], C->chs[0])], mx = mn;
        for (int ip = 0; ip < C->npars; ip++) for (int i1 = 0; i1 < C->nchs; i1++) for (int i2 = i1; i2 < C->nchs; i2++) {
            int64_t v = D->qYt[b3f_yt_off(C, C->pars[ip], u4, C->chs[i1], C->chs[i2])];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        yt_min[u4] = mn; yt_max[u4] = mx;
    }
    int64_t *s3_min = (int64_t*)b3f_xmalloc(C->NS3 * sizeof(int64_t)), *s3_max = (int64_t*)b3f_xmalloc(C->NS3 * sizeof(int64_t));
    for (int u4 = 0; u4 < C->NS3; u4++) {
        int64_t mn = D->qS3[b3f_s3_off(C, u4, C->chs[0], C->chs[0], C->chs[0])], mx = mn;
        for (int i1 = 0; i1 < C->nchs; i1++) for (int i2 = i1; i2 < C->nchs; i2++) for (int i3 = i2; i3 < C->nchs; i3++) {
            int64_t v = D->qS3[b3f_s3_off(C, u4, C->chs[i1], C->chs[i2], C->chs[i3])];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        s3_min[u4] = mn; s3_max[u4] = mx;
    }

    /* 3. Tail convex-hull child bounds for each kap in 0..KT+1 (h=1) */
    int64_t *tg_min = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t)), *tg_max = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t));
    int64_t *thv_min = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t)), *thv_max = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t));
    int64_t *tp4_min = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t)), *tp4_max = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t));
    int64_t *tyb_min = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t)), *tyb_max = (int64_t*)b3f_xmalloc((KT + 2) * sizeof(int64_t));
    for (int kap = 0; kap <= KT + 1; kap++) {
        int ids[4] = { C->INF_ID, b3f_t1(C, kap, 1), b3f_t2(C, kap, 1), b3f_t3(C, kap, 1) };
        int i0 = C->w_to_wi[ids[0]];
        tg_min[kap] = g_min[i0]; tg_max[kap] = g_max[i0];
        thv_min[kap] = hv_min[i0]; thv_max[kap] = hv_max[i0];
        tp4_min[kap] = p4_min[i0]; tp4_max[kap] = p4_max[i0];
        tyb_min[kap] = yb_min[i0]; tyb_max[kap] = yb_max[i0];
        for (int j = 1; j < 4; j++) {
            int iw = C->w_to_wi[ids[j]];
            if (g_min[iw] < tg_min[kap]) tg_min[kap] = g_min[iw];
            if (g_max[iw] > tg_max[kap]) tg_max[kap] = g_max[iw];
            if (hv_min[iw] < thv_min[kap]) thv_min[kap] = hv_min[iw];
            if (hv_max[iw] > thv_max[kap]) thv_max[kap] = hv_max[iw];
            if (p4_min[iw] < tp4_min[kap]) tp4_min[kap] = p4_min[iw];
            if (p4_max[iw] > tp4_max[kap]) tp4_max[kap] = p4_max[iw];
            if (yb_min[iw] < tyb_min[kap]) tyb_min[kap] = yb_min[iw];
            if (yb_max[iw] > tyb_max[kap]) tyb_max[kap] = yb_max[iw];
        }
    }

    /* Tail apex hC convex-hull bounds (0 included because INF has hC = 0) */
    int64_t *thc_min = (int64_t*)b3f_xcalloc(C->NT, sizeof(int64_t)), *thc_max = (int64_t*)b3f_xcalloc(C->NT, sizeof(int64_t));
    for (int h = 0; h <= 1; h++) for (int kap = 2; kap <= KT; kap++) {
        int idx = b3f_tidx(C, kap, h);
        int cas[3] = { C->cu_of_g[b3f_t1(C, kap, h)], C->cu_of_g[b3f_t2(C, kap, h)], C->cu_of_g[b3f_t3(C, kap, h)] };
        int64_t mn = 0, mx = 0;
        for (int j = 0; j < 3; j++) {
            if (hc_min[cas[j]] < mn) mn = hc_min[cas[j]];
            if (hc_max[cas[j]] > mx) mx = hc_max[cas[j]];
        }
        thc_min[idx] = mn; thc_max[idx] = mx;
    }

    /* 4. Sweep EVERY core vertex state (inner, finite tail, L=LF tail, VSTAR0, VSTAR1) using exact mpz_t!
     * Exact zero-sum edge re-charging for the lower bound c0:
     *   sum_{u in K(T)} (1/4) * (k_u - h_u) == (|E(K)| - |E(K)|) / 4 == 0 on EVERY rooted tree T!
     * We accumulate exact 4*S*c0 (with +(k - h)*S on finite states and +(KF + 1 - h)*S on VSTAR_h,
     * valid since k >= KF + 1 on VSTAR_h and +1/4 > 0), while keeping beta = 0 for exact S*C0. */
    mpz_t c0z4, C0z, loz, hiz, loz4, t;
    mpz_init(c0z4); mpz_init(C0z); mpz_init(loz); mpz_init(hiz); mpz_init(loz4); mpz_init(t);
    int first = 1;
    int64_t (*l5min)[L5_SHAPES] = b3f_xcalloc(C->INNER_N,sizeof *l5min);
    int64_t (*l5max)[L5_SHAPES] = b3f_xcalloc(C->INNER_N,sizeof *l5max);
    for (int i = 0; i < D->X->count; ++i) if (D->X->s[i].n == 5)
        for (int h = 0; h < 2; ++h) for (size_t r = 0; r < D->X->s[i].len[h]; ++r) {
            L5Tree row = {0};
            l5_decode(D->X,&row,i,h,-1,r);
            int g = row.g[l5_charge(&row)];
            assert(g >= 0 && C->D[g] < D->X->dc);
            int64_t v = D->q5[D->X->s[i].off[h]+r];
            if (v < l5min[g][i]) l5min[g][i] = v;
            if (v > l5max[g][i]) l5max[g][i] = v;
        }

    /* Finite (h, k <= KF, L <= LF) with D = L + k + h >= 2 */
    for (int h = 0; h <= 1; h++) for (int k = 0; k <= C->KF; k++) for (int L = 0; L <= C->LF; L++) {
        int Dv = L + k + h; if (Dv < 2) continue;
        int64_t fv = D->f[b3f_vidx(C, L, k, h)];
        mpz_set_si(loz, fv); mpz_set_si(hiz, fv);

        /* Child charges (g, hV, P4, Yb) when h == 1 */
        if (h == 1) {
            if (Dv < D0) {
                int iw = C->w_to_wi[b3f_inner(C, L, k, 1)];
                mpz_set_si(t, g_min[iw]);  mpz_add(loz, loz, t); mpz_set_si(t, g_max[iw]);  mpz_add(hiz, hiz, t);
                mpz_set_si(t, hv_min[iw]); mpz_add(loz, loz, t); mpz_set_si(t, hv_max[iw]); mpz_add(hiz, hiz, t);
                mpz_set_si(t, p4_min[iw]); mpz_add(loz, loz, t); mpz_set_si(t, p4_max[iw]); mpz_add(hiz, hiz, t);
                mpz_set_si(t, yb_min[iw]); mpz_add(loz, loz, t); mpz_set_si(t, yb_max[iw]); mpz_add(hiz, hiz, t);
            } else {
                int kap = (k > KT) ? (KT + 1) : k;
                mpz_set_si(t, tg_min[kap]);  mpz_add(loz, loz, t); mpz_set_si(t, tg_max[kap]);  mpz_add(hiz, hiz, t);
                mpz_set_si(t, thv_min[kap]); mpz_add(loz, loz, t); mpz_set_si(t, thv_max[kap]); mpz_add(hiz, hiz, t);
                mpz_set_si(t, tp4_min[kap]); mpz_add(loz, loz, t); mpz_set_si(t, tp4_max[kap]); mpz_add(hiz, hiz, t);
                mpz_set_si(t, tyb_min[kap]); mpz_add(loz, loz, t); mpz_set_si(t, tyb_max[kap]); mpz_add(hiz, hiz, t);
            }
        }
        /* Apex charges (hC, Yt, S3) */
        if (k >= 2 && k <= KT) {
            uint64_t c2 = (uint64_t)k * (uint64_t)(k - 1) / 2;
            if (Dv < D0) {
                int ca = C->cu_of_g[b3f_inner(C, L, k, h)];
                mpz_set_si(t, hc_min[ca]); mpz_addmul_ui(loz, t, c2);
                mpz_set_si(t, hc_max[ca]); mpz_addmul_ui(hiz, t, c2);
            } else {
                int idx = b3f_tidx(C, k, h);
                mpz_set_si(t, thc_min[idx]); mpz_addmul_ui(loz, t, c2);
                mpz_set_si(t, thc_max[idx]); mpz_addmul_ui(hiz, t, c2);
            }
        }
        if (h == 1 && k >= 2 && k <= C->KYT && Dv < C->DYT) {
            uint64_t c2 = (uint64_t)k * (uint64_t)(k - 1) / 2;
            int u4 = C->yt_lut[k][L];
            mpz_set_si(t, yt_min[u4]); mpz_addmul_ui(loz, t, c2);
            mpz_set_si(t, yt_max[u4]); mpz_addmul_ui(hiz, t, c2);
        }
        if (k >= 3 && k <= C->KS3 && Dv < C->DS3) {
            uint64_t c3 = (uint64_t)k * (uint64_t)(k - 1) * (uint64_t)(k - 2) / 6;
            int u4 = C->s3_lut[h][k][L];
            mpz_set_si(t, s3_min[u4]); mpz_addmul_ui(loz, t, c3);
            mpz_set_si(t, s3_max[u4]); mpz_addmul_ui(hiz, t, c3);
        }
        if (Dv < D->X->dc) {
            int g = b3f_inner(C,L,k,h);
            unsigned long sum_cap = 0;
            for (int i = 0; i < D->X->count; ++i) if (D->X->s[i].n == 5) {
                const L5Shape *s = &D->X->s[i];
                /* When the motif's apex is an endpoint, it must be the
                 * charge vertex's unique parent and so costs no choices. */
                unsigned long b = s->a == 1 ? (h ? l5_shape_charge_cap(D->X,s->child[0],k) : 0)
                                           : l5_shape_charge_cap(D->X,i,k);
                sum_cap += b;
                mpz_set_si(t,l5min[g][i]); mpz_addmul_ui(loz,t,b);
                mpz_set_si(t,l5max[g][i]); mpz_addmul_ui(hiz,t,b);
            }
            assert(sum_cap <= l5_charge_cap(Dv,D->X->dc));
        }
        /* Apply exact zero-sum +(1/4)*(k - h) edge re-charging to lower bound loz4 = 4*loz + (k - h)*S */
        mpz_mul_ui(loz4, loz, 4);
        mpz_set_si(t, (int64_t)(k - h));
        mpz_mul_si(t, t, H->S);
        mpz_add(loz4, loz4, t);

        if (first) { mpz_set(c0z4, loz4); mpz_set(C0z, hiz); first = 0; }
        else {
            if (mpz_cmp(loz4, c0z4) < 0) mpz_set(c0z4, loz4);
            if (mpz_cmp(hiz, C0z) > 0) mpz_set(C0z, hiz);
        }
    }
    /* Overflow k > KF (VSTAR0, VSTAR1): always D >= KF+1 >= D0, k > KT, so hC=Yt=S3=0; h=1 has tail kap=KT+1 charges.
     * Since k >= KF + 1 and +1/4 > 0, (1/4)*(k - h) >= (1/4)*(KF + 1 - h) for all k > KF. */
    for (int h = 0; h <= 1; h++) {
        int64_t fv = D->f[h ? C->VSTAR1 : C->VSTAR0];
        mpz_set_si(loz, fv); mpz_set_si(hiz, fv);
        if (h == 1) {
            int kap = KT + 1;
            mpz_set_si(t, tg_min[kap]);  mpz_add(loz, loz, t); mpz_set_si(t, tg_max[kap]);  mpz_add(hiz, hiz, t);
            mpz_set_si(t, thv_min[kap]); mpz_add(loz, loz, t); mpz_set_si(t, thv_max[kap]); mpz_add(hiz, hiz, t);
            mpz_set_si(t, tp4_min[kap]); mpz_add(loz, loz, t); mpz_set_si(t, tp4_max[kap]); mpz_add(hiz, hiz, t);
            mpz_set_si(t, tyb_min[kap]); mpz_add(loz, loz, t); mpz_set_si(t, tyb_max[kap]); mpz_add(hiz, hiz, t);
        }
        mpz_mul_ui(loz4, loz, 4);
        mpz_set_si(t, (int64_t)(C->KF + 1 - h));
        mpz_mul_si(t, t, H->S);
        mpz_add(loz4, loz4, t);
        if (mpz_cmp(loz4, c0z4) < 0) mpz_set(c0z4, loz4);
        if (mpz_cmp(hiz, C0z) > 0) mpz_set(C0z, hiz);
    }

    int64_t phi_T3 = D->f[b3f_vidx(C, 2, 0, 0)];
    gmp_printf("  [COMPARISON] 4*S*c0 = %Zd, S*C0 = %Zd, S*Phi(T_3) = %lld (S = %lld)\n",
               c0z4, C0z, (long long)phi_T3, (long long)H->S);
    mpz_set(mpq_numref(lower),c0z4); mpz_set_ui(mpq_denref(lower),4); mpq_canonicalize(lower);
    mpz_set(mpq_numref(upper),C0z); mpz_set_ui(mpq_denref(upper),1);
    int ok = 1;
    if (phi_T3 <= 0) { fprintf(stderr, "  FAIL: Phi(T_3) <= 0!\n"); ok = 0; }
    if (mpz_sgn(C0z) <= 0) { fprintf(stderr, "  FAIL: C0 <= 0!\n"); ok = 0; }
    if (H->mode == 1 && mpz_sgn(c0z4) <= 0) { fprintf(stderr, "  FAIL: UPPER cert requires c0 > 0!\n"); ok = 0; }
    if (H->mode == 2 && mpz_sgn(c0z4) <= 0) {
        printf("  [NOTE] LOWER cert has c0 <= 0 at some state; valid because lower bound uses ONLY Phi(T) <= C0*C(T) with C0 > 0 and Phi(T_3) > 0.\n");
    }

    free(g_min); free(g_max); free(hv_min); free(hv_max); free(p4_min); free(p4_max); free(yb_min); free(yb_max);
    free(hc_min); free(hc_max); free(yt_min); free(yt_max); free(s3_min); free(s3_max);
    free(tg_min); free(tg_max); free(thv_min); free(thv_max); free(tp4_min); free(tp4_max); free(tyb_min); free(tyb_max);
    free(thc_min); free(thc_max);
    free(l5min); free(l5max);
    mpz_clear(c0z4); mpz_clear(C0z); mpz_clear(loz); mpz_clear(hiz); mpz_clear(loz4); mpz_clear(t);
    return ok;
}

static int parse_i64_arg(const char *s, int64_t *out) {
    if (!s || !*s) return 0;
    errno = 0;
    char *end = NULL;
    long long v = strtoll(s, &end, 10);
    if (errno == ERANGE || !end || *end != '\0') return 0;
    *out = (int64_t)v;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <cert.bin> [--threads N] [--override-pq P Q] [--tree-check N]\n", argv[0]);
        return 2;
    }
    const char *path = argv[1];
    int nthreads = 1, override_pq = 0, ntrees_check = 30;
    int64_t ov_P = 0, ov_Q = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--threads") && i + 1 < argc) {
            int64_t tmp;
            if (!parse_i64_arg(argv[++i], &tmp) || tmp < 1 || tmp > 64) { fprintf(stderr, "Invalid --threads\n"); return 2; }
            nthreads = (int)tmp;
        } else if (!strcmp(argv[i], "--tree-check") && i + 1 < argc) {
            int64_t tmp;
            if (!parse_i64_arg(argv[++i], &tmp) || tmp < 0 || tmp > 10000) { fprintf(stderr, "Invalid --tree-check\n"); return 2; }
            ntrees_check = (int)tmp;
        } else if (!strcmp(argv[i], "--override-pq") && i + 2 < argc) {
            override_pq = 1;
            if (!parse_i64_arg(argv[++i], &ov_P) || !parse_i64_arg(argv[++i], &ov_Q)) { fprintf(stderr, "Invalid --override-pq\n"); return 2; }
        } else { fprintf(stderr, "Unknown arg: %s\n", argv[i]); return 2; }
    }

    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    if (fseek(fp, 0, SEEK_END) != 0) { fprintf(stderr, "FAIL: fseek end failed\n"); fclose(fp); return 1; }
    long fsz = ftell(fp);
    if (fsz < 0 || fseek(fp, 0, SEEK_SET) != 0) { fprintf(stderr, "FAIL: ftell/fseek failed\n"); fclose(fp); return 1; }
    CertB3FHeader H;
    L5Header E;
    if (fread(&H, sizeof(H), 1, fp) != 1 || !l5_header_io(fp,&E,0)) { fprintf(stderr, "Header read failed\n"); fclose(fp); return 1; }
    int hv = b3f_l5_validate_header(&H,&E);
    if (hv != 0) { fprintf(stderr, "FAIL: invalid header (err=%d)\n", hv); fclose(fp); return 1; }
    if (override_pq) {
        H.P = ov_P; H.Q = ov_Q;
        hv = b3f_l5_validate_header(&H,&E);
        if (hv != 0) { fprintf(stderr, "FAIL: invalid override header (err=%d)\n", hv); fclose(fp); return 1; }
    }

    B3FLayout C;
    b3f_layout_init(&C, &H);
    L5 *x = l5_init(&C,&E);
    size_t exp_bytes = l5_expected_bytes(x);
    if (!x->rows || E.rows != x->rows || exp_bytes == SIZE_MAX || (size_t)fsz != exp_bytes) {
        fprintf(stderr, "FAIL: file size/rows (%ld/%" PRIu64 ") != expected (%zu/%zu)\n", fsz,E.rows,exp_bytes,x->rows);
        fclose(fp); l5_free(x); b3f_layout_free(&C); return 1;
    }
    l5_init_fibers(x);

    CertData D;
    D.C = &C; D.H = &H; D.X = x;
    D.rpar   = read_i64_array(fp, C.INNER_N, CERT_B3F_MAX_RATE, "rpar");
    D.rch    = read_i64_array(fp, C.INNER_N, CERT_B3F_MAX_RATE, "rch");
    D.trpar  = read_i64_array(fp, C.NT,      CERT_B3F_MAX_RATE, "trpar");
    D.trch   = read_i64_array(fp, C.NT,      CERT_B3F_MAX_RATE, "trch");
    D.c2_f   = read_i64_array(fp, C.NT,      CERT_B3F_MAX_VAL,  "c2_f");
    D.c3_f   = read_i64_array(fp, C.NT,      CERT_B3F_MAX_VAL,  "c3_f");
    D.Yp_ref = read_i64_array(fp, C.yp_len,  CERT_B3F_MAX_VAL,  "Yp_ref");
    D.Yw_ref = read_i64_array(fp, C.yw_len,  CERT_B3F_MAX_VAL,  "Yw_ref");
    D.f      = read_i64_array(fp, C.NV,      CERT_B3F_MAX_VAL,  "f");
    D.g      = read_i64_array(fp, C.g_len,   CERT_B3F_MAX_VAL,  "g");
    D.hV     = read_i64_array(fp, C.hv_len,  CERT_B3F_MAX_VAL,  "hV");
    D.hC     = read_i64_array(fp, C.hc_len,  CERT_B3F_MAX_VAL,  "hC");
    D.qP4    = read_i64_array(fp, C.p4_len,  CERT_B3F_MAX_VAL,  "qP4");
    D.qYt    = read_i64_array(fp, C.yt_len,  CERT_B3F_MAX_VAL,  "qYt");
    D.qYb    = read_i64_array(fp, C.yb_len,  CERT_B3F_MAX_VAL,  "qYb");
    D.qS3    = read_i64_array(fp, C.s3_len,  CERT_B3F_MAX_VAL,  "qS3");
    D.q5 = read_i64_array(fp,x->rows,CERT_B3F_MAX_VAL,"q5");
    if (fgetc(fp) != EOF || ferror(fp) || fclose(fp) != 0) { fputs("FAIL: certificate read/close\n",stderr); return 1; }

    int mode = H.mode, D0 = C.D0, KT = C.KT, N = C.NG;
    int64_t P = H.P, Q = H.Q, S_r = H.S_r;
    int v_nb = b3f_vidx(&C, 1, 0, 1), iw_nb = b3f_inner(&C, 1, 0, 1);

    printf("Verifying %s | mode=%s | lambda = %lld/%lld | D0=%d KT=%d LF=%d KF=%d | threads=%d\n",
           path, (mode == 1) ? "UPPER" : "LOWER", (long long)P, (long long)Q, D0, KT, C.LF, C.KF, nthreads);
    fflush(stdout);

    long long fails_f = 0, fails_g = 0, fails_hv = 0, fails_hc = 0, fails_p4 = 0, fails_yt = 0, fails_yb = 0, fails_s3 = 0;

    /* 0. Check g(INF, INF) == 0 */
    if (D.g[C.INF_ID * N + C.INF_ID] != 0) {
        fprintf(stderr, "FAIL: g[INF, INF] = %lld != 0\n", (long long)D.g[C.INF_ID * N + C.INF_ID]);
        fails_g++;
    }

    /* 1. Verify f (exact GMP) */
    {
        GmpWork W; gw_init(&W);
        /* VSTAR0, VSTAR1 */
        for (int hb = 0; hb <= 1; hb++) {
            int64_t fv = D.f[hb ? C.VSTAR1 : C.VSTAR0];
            /* b = 0 corner: -lambda * fv */
            zf_set_si(&W.acc, 0, 1); zf_add_mul_si(&W.acc, -P, fv, (uint64_t)Q, &W);
            if (zf_fail(&W.acc, mode)) fails_f++;
            /* b = 1 corners: f[v_nb] + q_j - lambda * fv for j = 0..3 on [0, D0/(KF+1+hb)] */
            int64_t g0 = D.g[C.INF_ID * N + iw_nb];
            int64_t g1 = D.g[b3f_t1(&C, KT + 1, hb) * N + iw_nb];
            int64_t g2 = D.g[b3f_t2(&C, KT + 1, hb) * N + iw_nb];
            int64_t g3 = D.g[b3f_t3(&C, KT + 1, hb) * N + iw_nb];
            uint64_t den_m = (uint64_t)(C.KF + 1 + hb), tm_n = (uint64_t)D0, sm_n = den_m - tm_n;
            uint64_t d2 = den_m * den_m, d3 = d2 * den_m;
            for (int j = 0; j < 4; j++) {
                zf_set_si(&W.acc, D.f[v_nb], 1);
                zf_add_mul_si(&W.acc, -P, fv, (uint64_t)Q, &W);
                if (j == 0) zf_add_si(&W.acc, g0, 1, &W);
                else if (j == 1) { zf_add_mul_si(&W.acc, (int64_t)sm_n, g0, den_m, &W); zf_add_mul_si(&W.acc, (int64_t)tm_n, g1, den_m, &W); }
                else if (j == 2) {
                    zf_add_mul_si(&W.acc, (int64_t)(sm_n * sm_n), g0, d2, &W);
                    zf_add_mul_si(&W.acc, (int64_t)(2 * tm_n * sm_n), g1, d2, &W);
                    zf_add_mul_si(&W.acc, (int64_t)(tm_n * tm_n), g2, d2, &W);
                } else {
                    zf_add_mul_si(&W.acc, (int64_t)(sm_n * sm_n * sm_n), g0, d3, &W);
                    zf_add_mul_si(&W.acc, (int64_t)(3 * tm_n * sm_n * sm_n), g1, d3, &W);
                    zf_add_mul_si(&W.acc, (int64_t)(3 * tm_n * tm_n * sm_n), g2, d3, &W);
                    zf_add_mul_si(&W.acc, (int64_t)(tm_n * tm_n * tm_n), g3, d3, &W);
                }
                if (zf_fail(&W.acc, mode)) fails_f++;
            }
        }

        /* Finite (hb, k <= KF, L <= LF) */
        #pragma omp parallel num_threads(nthreads) reduction(+:fails_f)
        {
            GmpWork Wt; gw_init(&Wt);
            #pragma omp for schedule(dynamic)
            for (int kk = 0; kk < 2 * (C.KF + 1); kk++) {
                int hb = kk / (C.KF + 1), k = kk % (C.KF + 1);
                for (int L = 0; L <= C.LF; L++) {
                    int Dv = L + k + hb; if (Dv < 2) continue;
                    int u = b3f_vidx(&C, L, k, hb);
                    int64_t fv = D.f[u];
                    if (L == C.LF) {
                        /* Exact xfer_adj extremal bound over p in {0, c2/(3 D^2), c2/D^2 + c3/D^3} and r in {Rbar - (k+h)/2, Rbar} */
                        ZFrac xadj; mpz_init(xadj.num); mpz_init(xadj.den); zf_set_si(&xadj, 0, 1);
                        if (k <= KT) {
                            int idx = b3f_tidx(&C, k, hb);
                            int64_t c2 = D.c2_f[idx], c3 = D.c3_f[idx];
                            int64_t Rbar_q = hb * D.trpar[idx] + k * D.trch[idx];
                            /* 2 values of rr = (2*Rbar_q - m*S_r) / (2*S_r) for m in {0, k+hb} */
                            /* 3 values of p(x) on [0, 1/Dv]: 0, c2/(3*Dv^2), (c2*Dv + c3)/Dv^3 */
                            uint64_t D2 = (uint64_t)Dv * (uint64_t)Dv, D3 = D2 * (uint64_t)Dv;
                            int have_x = 0;
                            for (int ri = 0; ri < 2; ri++) {
                                int m = (ri == 0) ? 0 : (k + hb);
                                mpz_set_si(Wt.t1, Rbar_q); mpz_mul_ui(Wt.t1, Wt.t1, 2);
                                mpz_set_si(Wt.t2, S_r);    mpz_submul_ui(Wt.t1, Wt.t2, (uint64_t)m);
                                for (int pi = 0; pi < 3; pi++) {
                                    if (pi == 0) { zf_set_si(&Wt.cur, 0, 1); }
                                    else if (pi == 1) {
                                        mpz_mul_si(Wt.cur.num, Wt.t1, c2);
                                        mpz_set_ui(Wt.cur.den, (uint64_t)(2 * S_r)); mpz_mul_ui(Wt.cur.den, Wt.cur.den, 3 * D2);
                                    } else {
                                        mpz_set_si(Wt.t2, c2); mpz_mul_ui(Wt.t2, Wt.t2, (uint64_t)Dv);
                                        mpz_set_si(Wt.t3, c3); mpz_add(Wt.t2, Wt.t2, Wt.t3);
                                        mpz_mul(Wt.cur.num, Wt.t1, Wt.t2);
                                        mpz_set_ui(Wt.cur.den, (uint64_t)(2 * S_r)); mpz_mul_ui(Wt.cur.den, Wt.cur.den, D3);
                                    }
                                    if (!have_x) { mpz_set(xadj.num, Wt.cur.num); mpz_set(xadj.den, Wt.cur.den); have_x = 1; }
                                    else {
                                        int cmp = zf_cmp(&Wt.cur, &xadj, &Wt);
                                        if ((mode == 1 && cmp > 0) || (mode == 2 && cmp < 0)) { mpz_set(xadj.num, Wt.cur.num); mpz_set(xadj.den, Wt.cur.den); }
                                    }
                                }
                            }
                        }
                        /* Case 1: L == LF exact: bL*(f[ub] + f[v_nb] + g_embed - fv) + xadj - lam*fv */
                        int ub = b3f_vidx(&C, L - 1, k + 1, hb);
                        zf_set_si(&Wt.acc, 0, 1);
                        zf_add_z(&Wt.acc, xadj.num, xadj.den, &Wt);
                        zf_add_mul_si(&Wt.acc, -P, fv, (uint64_t)Q, &Wt);
                        /* bL * (f[ub] + f[v_nb] - fv) */
                        mpz_set_si(Wt.t3, D.f[ub]); mpz_set_si(Wt.t4, D.f[v_nb]); mpz_add(Wt.t3, Wt.t3, Wt.t4);
                        mpz_set_si(Wt.t4, fv);      mpz_sub(Wt.t3, Wt.t3, Wt.t4);
                        mpz_mul_ui(Wt.t3, Wt.t3, (uint64_t)L);
                        mpz_set_ui(Wt.t4, (uint64_t)Dv);
                        zf_add_z(&Wt.acc, Wt.t3, Wt.t4, &Wt);
                        /* bL * g(embed(L-1, k+1, hb), v_nb) -- always Dv >= LF >= D0 so tail! */
                        {
                            uint64_t s_n = (uint64_t)(Dv - D0), t_n = (uint64_t)D0;
                            uint64_t D3 = (uint64_t)Dv * (uint64_t)Dv * (uint64_t)Dv;
                            uint64_t D4 = D3 * (uint64_t)Dv;
                            int64_t g0 = D.g[C.INF_ID * N + iw_nb], g1 = D.g[b3f_t1(&C, k + 1, hb) * N + iw_nb];
                            int64_t g2 = D.g[b3f_t2(&C, k + 1, hb) * N + iw_nb], g3 = D.g[b3f_t3(&C, k + 1, hb) * N + iw_nb];
                            mpz_set_si(Wt.t3, g0); mpz_mul_ui(Wt.t3, Wt.t3, s_n * s_n * s_n);
                            mpz_set_si(Wt.t4, g1); mpz_addmul_ui(Wt.t3, Wt.t4, 3 * t_n * s_n * s_n);
                            mpz_set_si(Wt.t4, g2); mpz_addmul_ui(Wt.t3, Wt.t4, 3 * t_n * t_n * s_n);
                            mpz_set_si(Wt.t4, g3); mpz_addmul_ui(Wt.t3, Wt.t4, t_n * t_n * t_n);
                            mpz_mul_ui(Wt.t3, Wt.t3, (uint64_t)L);
                            mpz_set_ui(Wt.t4, D4);
                            zf_add_z(&Wt.acc, Wt.t3, Wt.t4, &Wt);
                        }
                        if (zf_fail(&Wt.acc, mode)) fails_f++;

                        /* Case 2: L >= LF + 1: b in [blo, 1], g in de Casteljau q0..q3 on [0, D0/(LF+1+k+hb)] */
                        int ubt = b3f_vidx(&C, C.LF, k + 1, hb);
                        int64_t g0 = D.g[C.INF_ID * N + iw_nb], g1 = D.g[b3f_t1(&C, k + 1, hb) * N + iw_nb];
                        int64_t g2 = D.g[b3f_t2(&C, k + 1, hb) * N + iw_nb], g3 = D.g[b3f_t3(&C, k + 1, hb) * N + iw_nb];
                        uint64_t dm = (uint64_t)(C.LF + 1 + k + hb), tn = (uint64_t)D0, sn = dm - tn;
                        uint64_t dm2 = dm * dm, dm3 = dm2 * dm;
                        uint64_t bn[2] = { (uint64_t)(C.LF + 1), 1 }, bd[2] = { dm, 1 };
                        for (int bi = 0; bi < 2; bi++) for (int j = 0; j < 4; j++) {
                            zf_set_si(&Wt.acc, 0, 1);
                            zf_add_z(&Wt.acc, xadj.num, xadj.den, &Wt);
                            zf_add_mul_si(&Wt.acc, -P, fv, (uint64_t)Q, &Wt);
                            /* b * (f[ubt] + f[v_nb] - fv + q_j) */
                            mpz_set_si(Wt.cur.num, D.f[ubt]); mpz_set_si(Wt.t1, D.f[v_nb]); mpz_add(Wt.cur.num, Wt.cur.num, Wt.t1);
                            mpz_set_si(Wt.t1, fv); mpz_sub(Wt.cur.num, Wt.cur.num, Wt.t1);
                            mpz_set_ui(Wt.cur.den, 1);
                            if (j == 0) zf_add_si(&Wt.cur, g0, 1, &Wt);
                            else if (j == 1) { zf_add_mul_si(&Wt.cur, (int64_t)sn, g0, dm, &Wt); zf_add_mul_si(&Wt.cur, (int64_t)tn, g1, dm, &Wt); }
                            else if (j == 2) {
                                zf_add_mul_si(&Wt.cur, (int64_t)(sn * sn), g0, dm2, &Wt);
                                zf_add_mul_si(&Wt.cur, (int64_t)(2 * tn * sn), g1, dm2, &Wt);
                                zf_add_mul_si(&Wt.cur, (int64_t)(tn * tn), g2, dm2, &Wt);
                            } else {
                                zf_add_mul_si(&Wt.cur, (int64_t)(sn * sn * sn), g0, dm3, &Wt);
                                zf_add_mul_si(&Wt.cur, (int64_t)(3 * tn * sn * sn), g1, dm3, &Wt);
                                zf_add_mul_si(&Wt.cur, (int64_t)(3 * tn * tn * sn), g2, dm3, &Wt);
                                zf_add_mul_si(&Wt.cur, (int64_t)(tn * tn * tn), g3, dm3, &Wt);
                            }
                            mpz_mul_ui(Wt.cur.num, Wt.cur.num, bn[bi]);
                            mpz_mul_ui(Wt.cur.den, Wt.cur.den, bd[bi]);
                            zf_add_z(&Wt.acc, Wt.cur.num, Wt.cur.den, &Wt);
                            if (zf_fail(&Wt.acc, mode)) fails_f++;
                        }
                        mpz_clear(xadj.num); mpz_clear(xadj.den);
                        continue;
                    }

                    /* L < LF */
                    int64_t f_up = D.f[b3f_vidx(&C, L + 1, k, hb)], df_up = f_up - fv;
                    /* Common base: L*df_up + bL*(f[ub] + f[v_nb] + g_embed - fv) - (P/Q)*fv */
                    zf_set_si(&Wt.acc, 0, 1);
                    zf_add_mul_si(&Wt.acc, (int64_t)L, df_up, 1, &Wt);
                    zf_add_mul_si(&Wt.acc, -P, fv, (uint64_t)Q, &Wt);
                    if (L >= 1) {
                        int ub = b3f_vidx(&C, L - 1, k + 1, hb);
                        mpz_set_si(Wt.t3, D.f[ub]); mpz_set_si(Wt.t4, D.f[v_nb]); mpz_add(Wt.t3, Wt.t3, Wt.t4);
                        mpz_set_si(Wt.t4, fv);      mpz_sub(Wt.t3, Wt.t3, Wt.t4);
                        mpz_mul_ui(Wt.t3, Wt.t3, (uint64_t)L);
                        mpz_set_ui(Wt.t4, (uint64_t)Dv);
                        zf_add_z(&Wt.acc, Wt.t3, Wt.t4, &Wt);
                        if (Dv < D0) {
                            int64_t gnb = D.g[b3f_inner(&C, L - 1, k + 1, hb) * N + iw_nb];
                            zf_add_mul_si(&Wt.acc, (int64_t)L, gnb, (uint64_t)Dv, &Wt);
                        } else {
                            uint64_t s_n = (uint64_t)(Dv - D0), t_n = (uint64_t)D0;
                            uint64_t D4 = (uint64_t)Dv * (uint64_t)Dv * (uint64_t)Dv * (uint64_t)Dv;
                            int64_t g0 = D.g[C.INF_ID * N + iw_nb], g1 = D.g[b3f_t1(&C, k + 1, hb) * N + iw_nb];
                            int64_t g2 = D.g[b3f_t2(&C, k + 1, hb) * N + iw_nb], g3 = D.g[b3f_t3(&C, k + 1, hb) * N + iw_nb];
                            mpz_set_si(Wt.t3, g0); mpz_mul_ui(Wt.t3, Wt.t3, s_n * s_n * s_n);
                            mpz_set_si(Wt.t4, g1); mpz_addmul_ui(Wt.t3, Wt.t4, 3 * t_n * s_n * s_n);
                            mpz_set_si(Wt.t4, g2); mpz_addmul_ui(Wt.t3, Wt.t4, 3 * t_n * t_n * s_n);
                            mpz_set_si(Wt.t4, g3); mpz_addmul_ui(Wt.t3, Wt.t4, t_n * t_n * t_n);
                            mpz_mul_ui(Wt.t3, Wt.t3, (uint64_t)L);
                            mpz_set_ui(Wt.t4, D4);
                            zf_add_z(&Wt.acc, Wt.t3, Wt.t4, &Wt);
                        }
                    }
                    if (Dv < D0) {
                        int gu = b3f_inner(&C, L, k, hb);
                        int64_t base_q = hb * D.rpar[gu] + k * D.rch[gu];
                        zf_add_mul_si(&Wt.acc, base_q, df_up, (uint64_t)S_r, &Wt);
                        if (zf_fail(&Wt.acc, mode)) fails_f++;
                    } else if (k <= KT) {
                        int idx = b3f_tidx(&C, k, hb);
                        int64_t base_q = hb * D.trpar[idx] + k * D.trch[idx];
                        uint64_t D3 = (uint64_t)Dv * (uint64_t)Dv * (uint64_t)Dv;
                        /* psi = (c2*Dv + c3) / D^3 */
                        mpz_set_si(Wt.t1, D.c2_f[idx]); mpz_mul_ui(Wt.t1, Wt.t1, (uint64_t)Dv);
                        mpz_set_si(Wt.t2, D.c3_f[idx]); mpz_add(Wt.t1, Wt.t1, Wt.t2); /* psi_n */
                        /* Add base_r * psi = (base_q * psi_n) / (S_r * D^3) */
                        mpz_mul_si(Wt.t3, Wt.t1, base_q);
                        mpz_set_ui(Wt.t4, (uint64_t)S_r); mpz_mul_ui(Wt.t4, Wt.t4, D3);
                        zf_add_z(&Wt.acc, Wt.t3, Wt.t4, &Wt);
                        /* Check rho = 0 corner */
                        if (zf_fail(&Wt.acc, mode)) fails_f++;
                        /* Check rho = (hb + k) / 2 corner: adds ((hb+k)/2) * (df_up - psi) */
                        if (hb + k > 0) {
                            mpz_set_si(Wt.t3, df_up); mpz_mul_ui(Wt.t3, Wt.t3, D3);
                            mpz_sub(Wt.t3, Wt.t3, Wt.t1); /* D^3*df_up - psi_n */
                            mpz_mul_ui(Wt.t3, Wt.t3, (uint64_t)(hb + k));
                            mpz_set_ui(Wt.t4, 2 * D3);
                            zf_add_z(&Wt.acc, Wt.t3, Wt.t4, &Wt);
                            if (zf_fail(&Wt.acc, mode)) fails_f++;
                        }
                    } else {
                        /* k > KT: rho in [0, (hb+k)/2] */
                        if (zf_fail(&Wt.acc, mode)) fails_f++;
                        if (hb + k > 0) {
                            zf_add_mul_si(&Wt.acc, (int64_t)(hb + k), df_up, 2, &Wt);
                            if (zf_fail(&Wt.acc, mode)) fails_f++;
                        }
                    }
                }
            }
            gw_clear(&Wt);
        }
        gw_clear(&W);
    }
    printf("  [f ] checked %d states: fails = %lld\n", C.NV, fails_f); fflush(stdout);

    /* Precompute exact dLf[u] = f[L+1,k,h] - f[L,k,h] for inner states */
    int64_t *dLf = (int64_t*)b3f_xcalloc(C.INNER_N, sizeof(int64_t));
    for (int u = 0; u < C.INNER_N; u++)
        dLf[u] = D.f[b3f_vidx(&C, C.L[u] + 1, C.K[u], C.H[u])] - D.f[b3f_vidx(&C, C.L[u], C.K[u], C.H[u])];

    /* 2. Verify g */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_g)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int ip = 0; ip < C.npars; ip++) {
            int p = C.pars[ip], pk = C.kind[p], Lp = C.L[p], kp = C.K[p], hp = C.H[p];
            int p_up = b3f_next_up(&C, p), p_bud = b3f_next_bud(&C, p);
            int p_t1 = (pk >= 3) ? b3f_t1(&C, kp, hp) : p, p_t2 = (pk >= 4) ? b3f_t2(&C, kp, hp) : p;
            int p_t1b = (pk >= 3 && kp <= KT) ? b3f_t1(&C, kp + 1, hp) : p_t1;
            int p_t2b = (pk >= 4 && kp <= KT) ? b3f_t2(&C, kp + 1, hp) : p_t2;
            int p_in_hu = (C.hu_of_g[p] >= 0), p_in_cu = (C.cu_of_g[p] >= 0);
            int pidx_t = (pk >= 2 && kp <= KT) ? b3f_tidx(&C, kp, hp) : 0;
            int cu_p1b = (pk >= 2 && kp + 1 <= KT && kp + 1 >= 2) ? C.cu_of_g[b3f_t1(&C, kp + 1, hp)] : -1;
            int cu_p2b = (pk >= 3 && kp + 1 <= KT && kp + 1 >= 2) ? C.cu_of_g[b3f_t2(&C, kp + 1, hp)] : -1;
            int cu_p3b = (pk == 4 && kp + 1 <= KT && kp + 1 >= 2) ? C.cu_of_g[b3f_t3(&C, kp + 1, hp)] : ((pk == 0 && Lp >= 1) ? C.cu_of_g[p_bud] : -1);

            for (int iw = 0; iw < C.nchs; iw++) {
                int w = C.chs[iw], wk = C.kind[w], Lw = C.L[w], kw = C.K[w];
                if (pk == 1 && wk == 1) continue;
                int w_up = b3f_next_up(&C, w), w_bud = b3f_next_bud(&C, w);
                int w_t1 = (wk >= 3) ? b3f_t1(&C, kw, 1) : w, w_t2 = (wk >= 4) ? b3f_t2(&C, kw, 1) : w;
                int w_t1b = (wk >= 3 && kw <= KT) ? b3f_t1(&C, kw + 1, 1) : w_t1;
                int w_t2b = (wk >= 4 && kw <= KT) ? b3f_t2(&C, kw + 1, 1) : w_t2;
                int w_in_hu = (C.hu_of_g[w] >= 0);
                int widx_t = (wk >= 2 && kw <= KT) ? b3f_tidx(&C, kw, 1) : 0;
                int hu_w1b = (wk >= 2 && kw + 1 <= KT) ? C.hu_of_g[b3f_t1(&C, kw + 1, 1)] : -1;
                int hu_w2b = (wk >= 3 && kw + 1 <= KT) ? C.hu_of_g[b3f_t2(&C, kw + 1, 1)] : -1;
                int hu_w3b = (wk == 4 && kw + 1 <= KT) ? C.hu_of_g[b3f_t3(&C, kw + 1, 1)] : ((wk == 0 && Lw >= 1) ? C.hu_of_g[w_bud] : -1);

                int64_t v = D.g[p * N + w], gp_plus = D.g[p_up * N + w], gw_plus = D.g[p * N + w_up];
                int64_t gp_bud = D.g[p_bud * N + w] + ((pk == 0 && cu_p3b >= 0) ? D.hC[b3f_hc_off(&C, cu_p3b, w, iw_nb)] : 0);
                int64_t gw_bud = D.g[p * N + w_bud] + ((wk == 0 && hu_w3b >= 0) ? D.hV[b3f_hv_off(&C, p, hu_w3b, iw_nb)] : 0);

                zf_set_si(&W.acc, 0, 1);
                zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);

                /* Source terms from f */
                if (pk == 0) add_centered_source(&W.acc, &C, w, D.rch[p], S_r, dLf[p], &W);
                else if (pk == 3 && kp <= KT) {
                    /* (inv_w - trch[pidx_t]) * c2_f[pidx_t] / (3 * D0^2) */
                    int64_t inum, iden; invD_rat(&C, w, &inum, &iden);
                    mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                    mpz_set_si(W.t4, D.trch[pidx_t]); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                    mpz_mul_si(W.t3, W.t3, D.c2_f[pidx_t]);
                    mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r); mpz_mul_ui(W.t4, W.t4, (uint64_t)(3 * D0 * D0));
                    zf_add_z(&W.acc, W.t3, W.t4, &W);
                } else if (pk == 4 && kp <= KT) {
                    int64_t inum, iden; invD_rat(&C, w, &inum, &iden);
                    mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                    mpz_set_si(W.t4, D.trch[pidx_t]); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                    mpz_set_si(W.t1, D.c2_f[pidx_t]); mpz_mul_ui(W.t1, W.t1, (uint64_t)D0);
                    mpz_set_si(W.t2, D.c3_f[pidx_t]); mpz_add(W.t1, W.t1, W.t2);
                    mpz_mul(W.t3, W.t3, W.t1);
                    mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r); mpz_mul_ui(W.t4, W.t4, (uint64_t)(D0 * D0 * D0));
                    zf_add_z(&W.acc, W.t3, W.t4, &W);
                }

                if (wk == 0) add_centered_source(&W.acc, &C, p, D.rpar[w], S_r, dLf[w], &W);
                else if (wk == 3 && kw <= KT) {
                    int64_t inum, iden; invD_rat(&C, p, &inum, &iden);
                    mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                    mpz_set_si(W.t4, D.trpar[widx_t]); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                    mpz_mul_si(W.t3, W.t3, D.c2_f[widx_t]);
                    mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r); mpz_mul_ui(W.t4, W.t4, (uint64_t)(3 * D0 * D0));
                    zf_add_z(&W.acc, W.t3, W.t4, &W);
                } else if (wk == 4 && kw <= KT) {
                    int64_t inum, iden; invD_rat(&C, p, &inum, &iden);
                    mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                    mpz_set_si(W.t4, D.trpar[widx_t]); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                    mpz_set_si(W.t1, D.c2_f[widx_t]); mpz_mul_ui(W.t1, W.t1, (uint64_t)D0);
                    mpz_set_si(W.t2, D.c3_f[widx_t]); mpz_add(W.t1, W.t1, W.t2);
                    mpz_mul(W.t3, W.t3, W.t1);
                    mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r); mpz_mul_ui(W.t4, W.t4, (uint64_t)(D0 * D0 * D0));
                    zf_add_z(&W.acc, W.t3, W.t4, &W);
                }

                /* Inner transitions */
                if (pk == 0) {
                    int64_t dp = gp_plus - v;
                    zf_add_mul_si(&W.acc, (int64_t)Lp, dp, 1, &W);
                    if (wk == 0) zf_add_si(&W.acc, dp, (uint64_t)C.D[w], &W);
                    int64_t rq = (p_in_hu ? D.rpar[p] : 0) + (p_in_cu ? (int64_t)(kp - 1) * D.rch[p] : 0);
                    if (rq != 0) zf_add_mul_si(&W.acc, rq, dp, (uint64_t)S_r, &W);
                    if (Lp >= 1) zf_add_mul_si(&W.acc, (int64_t)Lp, gp_bud - v, (uint64_t)C.D[p], &W);
                    int rem_h = (hp && !p_in_hu ? 1 : 0) + (kp >= 2 && !p_in_cu ? (kp - 1) : 0);
                    if (rem_h > 0 && ((mode == 1 && dp > 0) || (mode == 2 && dp < 0))) zf_add_mul_si(&W.acc, (int64_t)rem_h, dp, 2, &W);
                }
                if (wk == 0) {
                    int64_t dw = gw_plus - v;
                    zf_add_mul_si(&W.acc, (int64_t)Lw, dw, 1, &W);
                    if (pk == 0) zf_add_si(&W.acc, dw, (uint64_t)C.D[p], &W);
                    int64_t rq = w_in_hu ? (int64_t)kw * D.rch[w] : 0;
                    if (rq != 0) zf_add_mul_si(&W.acc, rq, dw, (uint64_t)S_r, &W);
                    if (Lw >= 1) zf_add_mul_si(&W.acc, (int64_t)Lw, gw_bud - v, (uint64_t)C.D[w], &W);
                    int rem_h = (kw >= 1 && !w_in_hu) ? kw : 0;
                    if (rem_h > 0 && ((mode == 1 && dw > 0) || (mode == 2 && dw < 0))) zf_add_mul_si(&W.acc, (int64_t)rem_h, dw, 2, &W);
                }

                /* Tail p */
                if (pk >= 2) {
                    int64_t ln, ld, hn, hd; get_rin1(&C, w, &ln, &ld, &hn, &hd);
                    int64_t pF0 = D.g[C.INF_ID * N + w], pF1 = (pk >= 3) ? D.g[p_t1 * N + w] : v, pF2 = (pk >= 4) ? D.g[p_t2 * N + w] : v;
                    int64_t pF1b = D.g[(pk == 2 ? p_bud : p_t1b) * N + w] + (cu_p1b >= 0 ? D.hC[b3f_hc_off(&C, cu_p1b, w, iw_nb)] : 0);
                    int64_t pF2b = (pk >= 3) ? (D.g[(pk == 3 ? p_bud : p_t2b) * N + w] + (cu_p2b >= 0 ? D.hC[b3f_hc_off(&C, cu_p2b, w, iw_nb)] : 0)) : 0;
                    int64_t pF3b = (pk == 4) ? (D.g[p_bud * N + w] + (cu_p3b >= 0 ? D.hC[b3f_hc_off(&C, cu_p3b, w, iw_nb)] : 0)) : 0;
                    int64_t pd0 = (wk == 0) ? (D.g[C.INF_ID * N + w_up] - pF0) : 0;
                    int64_t pd1 = (pk >= 3 && wk == 0) ? (D.g[p_t1 * N + w_up] - pF1) : 0;
                    int64_t pd2 = (pk >= 4 && wk == 0) ? (D.g[p_t2 * N + w_up] - pF2) : 0;
                    int64_t p_Rbar = (kp <= KT) ? (hp * D.trpar[pidx_t] + (int64_t)(kp - 1) * D.trch[pidx_t]) : 0;
                    int64_t p_Yref = (pk == 4 && kp <= KT) ? D.Yp_ref[b3f_yp_off(&C, pidx_t, w)] : 0;
                    add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, pk, kp, hp, ln, ld, hn, hd,
                                 kp + hp - 1, (kp <= KT), p_Rbar, p_Yref,
                                 pF0, pF1, pF2, v, pF1b, pF2b, pF3b, wk == 0, pd0, pd1, pd2, gw_plus);
                }
                /* Tail w */
                if (wk >= 2) {
                    int64_t ln, ld, hn, hd; get_rin1(&C, p, &ln, &ld, &hn, &hd);
                    int64_t wF0 = D.g[p * N + C.INF_ID], wF1 = (wk >= 3) ? D.g[p * N + w_t1] : v, wF2 = (wk >= 4) ? D.g[p * N + w_t2] : v;
                    int64_t wF1b = D.g[p * N + (wk == 2 ? w_bud : w_t1b)] + (hu_w1b >= 0 ? D.hV[b3f_hv_off(&C, p, hu_w1b, iw_nb)] : 0);
                    int64_t wF2b = (wk >= 3) ? (D.g[p * N + (wk == 3 ? w_bud : w_t2b)] + (hu_w2b >= 0 ? D.hV[b3f_hv_off(&C, p, hu_w2b, iw_nb)] : 0)) : 0;
                    int64_t wF3b = (wk == 4) ? (D.g[p * N + w_bud] + (hu_w3b >= 0 ? D.hV[b3f_hv_off(&C, p, hu_w3b, iw_nb)] : 0)) : 0;
                    int64_t wd0 = (pk == 0) ? (D.g[p_up * N + C.INF_ID] - wF0) : 0;
                    int64_t wd1 = (wk >= 3 && pk == 0) ? (D.g[p_up * N + w_t1] - wF1) : 0;
                    int64_t wd2 = (wk >= 4 && pk == 0) ? (D.g[p_up * N + w_t2] - wF2) : 0;
                    int64_t w_Rbar = (kw <= KT) ? ((int64_t)kw * D.trch[widx_t]) : 0;
                    int64_t w_Yref = (wk == 4 && kw <= KT && kw >= 1) ? D.Yw_ref[b3f_yw_off(&C, p, kw)] : 0;
                    add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, wk, kw, 1, ln, ld, hn, hd,
                                 kw, (kw <= KT), w_Rbar, w_Yref,
                                 wF0, wF1, wF2, v, wF1b, wF2b, wF3b, pk == 0, wd0, wd1, wd2, gp_plus);
                }

                if (zf_fail(&W.acc, mode)) fails_g++;
            }
        }
        gw_clear(&W);
    }
    printf("  [g ] checked %zu states: fails = %lld\n", (size_t)C.npars * C.nchs, fails_g); fflush(stdout);

    /* 3. Verify FULL hV */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_hv)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int ua = 0; ua < C.NHU; ua++) {
            int gu = C.hu_g[ua], uk = C.kind[gu], Lu = C.L[gu], ku = C.K[gu], Du = C.D[gu];
            int gup = (uk == 0) ? b3f_next_up(&C, gu) : C.INF_ID;
            int ua_up = (uk == 0) ? C.hu_of_g[gup] : -1;
            int gu_bud = b3f_next_bud(&C, gu);
            int ua_bud = (uk == 0 && Lu == 0) ? -1 : ((ku < KT) ? C.hu_of_g[gu_bud] : -1);
            int u_t1 = (uk >= 3) ? C.hu_of_g[b3f_t1(&C, ku, 1)] : ua;
            int u_t2 = (uk >= 4) ? C.hu_of_g[b3f_t2(&C, ku, 1)] : ua;
            int u_t1b = (uk >= 3 && ku < KT) ? C.hu_of_g[b3f_t1(&C, ku + 1, 1)] : -1;
            int u_t2b = (uk >= 4 && ku < KT) ? C.hu_of_g[b3f_t2(&C, ku + 1, 1)] : -1;
            int uidx_t = (uk >= 2) ? b3f_tidx(&C, ku, 1) : 0;
            int64_t urpar_q = (uk == 0) ? D.rpar[gu] : D.trpar[uidx_t];
            int64_t urch_q  = (uk == 0) ? D.rch[gu]  : D.trch[uidx_t];

            int u_p4  = (uk == 0 && ku <= C.KP4 && Du < C.DP4) ? C.p4_lut[ku][Lu] : -1;
            int u_yt  = (uk == 0 && ku >= 2 && ku <= C.KYT && Du < C.DYT) ? C.yt_lut[ku][Lu] : -1;
            int ub_yt = (uk == 0 && Lu >= 1 && ku + 1 >= 2 && ku + 1 <= C.KYT && Du < C.DYT) ? C.yt_lut[ku + 1][Lu - 1] : -1;
            int u_ybw = (uk == 0 && ku <= C.KYB && Du < C.DYB) ? C.ybw_lut[ku][Lu] : -1;

            for (int ip = 0; ip < C.npars; ip++) {
                int p = C.pars[ip], pk = C.kind[p], Lp = C.L[p], kp = C.K[p], hp = C.H[p];
                int p_up = b3f_next_up(&C, p), p_bud = b3f_next_bud(&C, p);
                int p_t1 = (pk >= 3) ? b3f_t1(&C, kp, hp) : p, p_t2 = (pk >= 4) ? b3f_t2(&C, kp, hp) : p;
                int p_t1b = (pk >= 3 && kp <= KT) ? b3f_t1(&C, kp + 1, hp) : p_t1;
                int p_t2b = (pk >= 4 && kp <= KT) ? b3f_t2(&C, kp + 1, hp) : p_t2;
                int p_p4   = (pk == 0 && hp == 1 && kp >= 1 && kp <= C.KP4 && C.D[p] < C.DP4) ? C.p4_lut[kp][Lp] : -1;
                int p_ybu  = (pk == 0 && kp >= 2 && kp <= C.KYB && C.D[p] < C.DYB) ? C.ybu_lut[hp][kp][Lp] : -1;
                int pb_ybu = (pk == 0 && Lp >= 1 && kp + 1 >= 2 && kp + 1 <= C.KYB && C.D[p] < C.DYB) ? C.ybu_lut[hp][kp + 1][Lp - 1] : -1;
                int x_pp = (p_p4 >= 0 && u_p4 >= 0), x_pc = (p_ybu >= 0 && u_ybw >= 0);

                int64_t dgpu_num = 0; int dgpu_divD0 = 0;
                if (uk == 0) dgpu_num = D.g[p * N + gup] - D.g[p * N + gu];
                else if (uk == 3) { dgpu_num = D.g[p * N + C.INF_ID] - D.g[p * N + b3f_t1(&C, ku, 1)]; dgpu_divD0 = 1; }
                else if (uk == 4) dgpu_num = D.Yw_ref[b3f_yw_off(&C, p, ku)];

                for (int iw = 0; iw < C.nchs; iw++) {
                    int w = C.chs[iw], wk = C.kind[w], Lw = C.L[w], kw = C.K[w];
                    int w_up = b3f_next_up(&C, w), w_bud = b3f_next_bud(&C, w);
                    int w_t1 = (wk >= 3) ? b3f_t1(&C, kw, 1) : w, w_t2 = (wk >= 4) ? b3f_t2(&C, kw, 1) : w;
                    int w_t1b = (wk >= 3 && kw <= KT) ? b3f_t1(&C, kw + 1, 1) : w_t1;
                    int w_t2b = (wk >= 4 && kw <= KT) ? b3f_t2(&C, kw + 1, 1) : w_t2;
                    int w_p4  = (wk == 0 && kw >= 1 && kw <= C.KP4 && C.D[w] < C.DP4) ? C.p4_lut[kw][Lw] : -1;
                    int wb_p4 = (wk == 0 && Lw >= 1 && kw + 1 <= C.KP4 && C.D[w] < C.DP4) ? C.p4_lut[kw + 1][Lw - 1] : -1;
                    int x_wc = (u_p4 >= 0 && w_p4 >= 0), x_uz = (u_yt >= 0);

                    int64_t dguw_num = 0; int dguw_divD0 = 0;
                    if (uk == 0) dguw_num = D.g[gup * N + w] - D.g[gu * N + w];
                    else if (uk == 3) { dguw_num = D.g[C.INF_ID * N + w] - D.g[b3f_t1(&C, ku, 1) * N + w]; dguw_divD0 = 1; }
                    else if (uk == 4) dguw_num = D.Yp_ref[b3f_yp_off(&C, uidx_t, w)];

                    int64_t v = D.hV[b3f_hv_off(&C, p, ua, w)];
                    int64_t hp_up = D.hV[b3f_hv_off(&C, p_up, ua, w)], hu_up = (ua_up >= 0) ? D.hV[b3f_hv_off(&C, p, ua_up, w)] : 0, hw_up = D.hV[b3f_hv_off(&C, p, ua, w_up)];
                    int64_t hp_bud = D.hV[b3f_hv_off(&C, p_bud, ua, w)], hu_bud = (ua_bud >= 0) ? D.hV[b3f_hv_off(&C, p, ua_bud, w)] : 0, hw_bud = D.hV[b3f_hv_off(&C, p, ua, w_bud)];
                    int64_t nb_p = (pb_ybu >= 0 && u_ybw >= 0) ? D.qYb[b3f_yb_off(&C, pb_ybu, u_ybw, iw_nb, w)] : 0;
                    int64_t nb_u = (ub_yt >= 0) ? D.qYt[b3f_yt_off(&C, p, ub_yt, w, iw_nb)] : 0;
                    int64_t nb_w = (u_p4 >= 0 && wb_p4 >= 0) ? D.qP4[b3f_p4_off(&C, p, u_p4, wb_p4, iw_nb)] : 0;

                    zf_set_si(&W.acc, 0, 1);
                    zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);

                    /* Sources: (inv_p - urpar)*dguw + (inv_w - urch)*dgpu */
                    if (dguw_num != 0) {
                        int64_t inum, iden; invD_rat(&C, p, &inum, &iden);
                        mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                        mpz_set_si(W.t4, urpar_q); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                        mpz_mul_si(W.t3, W.t3, dguw_num);
                        mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r);
                        if (dguw_divD0) mpz_mul_ui(W.t4, W.t4, (uint64_t)D0);
                        zf_add_z(&W.acc, W.t3, W.t4, &W);
                    }
                    if (dgpu_num != 0) {
                        int64_t inum, iden; invD_rat(&C, w, &inum, &iden);
                        mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                        mpz_set_si(W.t4, urch_q); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                        mpz_mul_si(W.t3, W.t3, dgpu_num);
                        mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r);
                        if (dgpu_divD0) mpz_mul_ui(W.t4, W.t4, (uint64_t)D0);
                        zf_add_z(&W.acc, W.t3, W.t4, &W);
                    }

                    if (pk == 0) {
                        int64_t dp = hp_up - v;
                        zf_add_mul_si(&W.acc, (int64_t)Lp, dp, 1, &W);
                        if (uk == 0) zf_add_si(&W.acc, dp, (uint64_t)Du, &W);
                        int64_t rq = (x_pp ? D.rpar[p] : 0) + (x_pc ? (int64_t)(kp - 1) * D.rch[p] : 0);
                        if (rq != 0) zf_add_mul_si(&W.acc, rq, dp, (uint64_t)S_r, &W);
                        if (Lp >= 1) zf_add_mul_si(&W.acc, (int64_t)Lp, hp_bud + nb_p - v, (uint64_t)C.D[p], &W);
                        int rh = (hp && !x_pp ? 1 : 0) + (kp >= 2 && !x_pc ? (kp - 1) : 0);
                        if (rh > 0 && ((mode == 1 && dp > 0) || (mode == 2 && dp < 0))) zf_add_mul_si(&W.acc, (int64_t)rh, dp, 2, &W);
                    }
                    if (wk == 0) {
                        int64_t dw = hw_up - v;
                        zf_add_mul_si(&W.acc, (int64_t)Lw, dw, 1, &W);
                        if (uk == 0) zf_add_si(&W.acc, dw, (uint64_t)Du, &W);
                        int64_t rq = x_wc ? (int64_t)kw * D.rch[w] : 0;
                        if (rq != 0) zf_add_mul_si(&W.acc, rq, dw, (uint64_t)S_r, &W);
                        if (Lw >= 1) zf_add_mul_si(&W.acc, (int64_t)Lw, hw_bud + nb_w - v, (uint64_t)C.D[w], &W);
                        int rh = (kw >= 1 && !x_wc) ? kw : 0;
                        if (rh > 0 && ((mode == 1 && dw > 0) || (mode == 2 && dw < 0))) zf_add_mul_si(&W.acc, (int64_t)rh, dw, 2, &W);
                    }
                    if (uk == 0) {
                        int64_t du = hu_up - v;
                        zf_add_mul_si(&W.acc, (int64_t)Lu, du, 1, &W);
                        if (pk == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[p], &W);
                        if (wk == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w], &W);
                        int64_t rq = x_uz ? (int64_t)(ku - 1) * D.rch[gu] : 0;
                        if (rq != 0) zf_add_mul_si(&W.acc, rq, du, (uint64_t)S_r, &W);
                        if (Lu >= 1) zf_add_mul_si(&W.acc, (int64_t)Lu, hu_bud + nb_u - v, (uint64_t)Du, &W);
                        int rh = (ku >= 2 && !x_uz) ? (ku - 1) : 0;
                        if (rh > 0 && ((mode == 1 && du > 0) || (mode == 2 && du < 0))) zf_add_mul_si(&W.acc, (int64_t)rh, du, 2, &W);
                    }

                    if (pk >= 2) {
                        int64_t ln, ld, hn, hd; get_rin1(&C, gu, &ln, &ld, &hn, &hd);
                        int64_t pF0 = D.hV[b3f_hv_off(&C, C.INF_ID, ua, w)], pF1 = (pk >= 3) ? D.hV[b3f_hv_off(&C, p_t1, ua, w)] : v, pF2 = (pk >= 4) ? D.hV[b3f_hv_off(&C, p_t2, ua, w)] : v;
                        int64_t pF1b = (pk >= 3) ? D.hV[b3f_hv_off(&C, p_t1b, ua, w)] : hp_bud, pF2b = (pk >= 4) ? D.hV[b3f_hv_off(&C, p_t2b, ua, w)] : hp_bud;
                        int64_t pd0 = (uk == 0) ? (((ua_up >= 0) ? D.hV[b3f_hv_off(&C, C.INF_ID, ua_up, w)] : 0) - pF0) : 0;
                        int64_t pd1 = (pk >= 3 && uk == 0) ? (((ua_up >= 0) ? D.hV[b3f_hv_off(&C, p_t1, ua_up, w)] : 0) - pF1) : 0;
                        int64_t pd2 = (pk >= 4 && uk == 0) ? (((ua_up >= 0) ? D.hV[b3f_hv_off(&C, p_t2, ua_up, w)] : 0) - pF2) : 0;
                        add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, pk, kp, hp, ln, ld, hn, hd, kp + hp - 1, 0, 0, 0,
                                     pF0, pF1, pF2, v, (pk == 2) ? hp_bud : pF1b, (pk == 3) ? hp_bud : pF2b, hp_bud, uk == 0, pd0, pd1, pd2, hu_up);
                    }
                    if (wk >= 2) {
                        int64_t ln, ld, hn, hd; get_rin1(&C, gu, &ln, &ld, &hn, &hd);
                        int64_t wF0 = D.hV[b3f_hv_off(&C, p, ua, C.INF_ID)], wF1 = (wk >= 3) ? D.hV[b3f_hv_off(&C, p, ua, w_t1)] : v, wF2 = (wk >= 4) ? D.hV[b3f_hv_off(&C, p, ua, w_t2)] : v;
                        int64_t wF1b = (wk >= 3) ? D.hV[b3f_hv_off(&C, p, ua, w_t1b)] : hw_bud, wF2b = (wk >= 4) ? D.hV[b3f_hv_off(&C, p, ua, w_t2b)] : hw_bud;
                        int64_t wd0 = (uk == 0) ? (((ua_up >= 0) ? D.hV[b3f_hv_off(&C, p, ua_up, C.INF_ID)] : 0) - wF0) : 0;
                        int64_t wd1 = (wk >= 3 && uk == 0) ? (((ua_up >= 0) ? D.hV[b3f_hv_off(&C, p, ua_up, w_t1)] : 0) - wF1) : 0;
                        int64_t wd2 = (wk >= 4 && uk == 0) ? (((ua_up >= 0) ? D.hV[b3f_hv_off(&C, p, ua_up, w_t2)] : 0) - wF2) : 0;
                        add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, wk, kw, 1, ln, ld, hn, hd, kw, 0, 0, 0,
                                     wF0, wF1, wF2, v, (wk == 2) ? hw_bud : wF1b, (wk == 3) ? hw_bud : wF2b, hw_bud, uk == 0, wd0, wd1, wd2, hu_up);
                    }
                    if (uk >= 2) {
                        int64_t ln, ld, hn, hd; get_rin2(&C, p, w, &ln, &ld, &hn, &hd);
                        int64_t uF1 = (uk >= 3) ? D.hV[b3f_hv_off(&C, p, u_t1, w)] : v, uF2 = (uk >= 4) ? D.hV[b3f_hv_off(&C, p, u_t2, w)] : v;
                        int64_t uF1b = (uk == 2) ? hu_bud : ((uk >= 3 && u_t1b >= 0) ? D.hV[b3f_hv_off(&C, p, u_t1b, w)] : 0);
                        int64_t uF2b = (uk == 3) ? hu_bud : ((uk >= 4 && u_t2b >= 0) ? D.hV[b3f_hv_off(&C, p, u_t2b, w)] : 0);
                        int u_nin = (pk == 0 ? 1 : 0) + (wk == 0 ? 1 : 0);
                        int64_t ud1 = (uk >= 3 && u_nin > 0) ? ((pk == 0 ? (D.hV[b3f_hv_off(&C, p_up, u_t1, w)] - uF1) : 0) + (wk == 0 ? (D.hV[b3f_hv_off(&C, p, u_t1, w_up)] - uF1) : 0)) : 0;
                        int64_t ud2 = (uk >= 4 && u_nin > 0) ? ((pk == 0 ? (D.hV[b3f_hv_off(&C, p_up, u_t2, w)] - uF2) : 0) + (wk == 0 ? (D.hV[b3f_hv_off(&C, p, u_t2, w_up)] - uF2) : 0)) : 0;
                        int64_t uF3_up = (uk == 4 && u_nin > 0) ? ((pk == 0 ? hp_up : 0) + (wk == 0 ? hw_up : 0)) : 0;
                        add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, uk, ku, 1, ln, ld, hn, hd, ku - 1, 0, 0, 0,
                                     0, uF1, uF2, v, uF1b, uF2b, hu_bud, u_nin, 0, ud1, ud2, uF3_up);
                    }

                    if (zf_fail(&W.acc, mode)) fails_hv++;
                }
            }
        }
        gw_clear(&W);
    }
    printf("  [hV] checked %zu states: fails = %lld\n", C.hv_len, fails_hv); fflush(stdout);

    /* 4. Verify FULL hC */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_hc)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int ca = 0; ca < C.NCU; ca++) {
            int gu = C.cu_g[ca], uk = C.kind[gu], Lu = C.L[gu], ku = C.K[gu], hu = C.H[gu], Du = C.D[gu];
            int gup = (uk == 0) ? b3f_next_up(&C, gu) : C.INF_ID;
            int ca_up = (uk == 0) ? C.cu_of_g[gup] : -1;
            int gu_bud = b3f_next_bud(&C, gu);
            int ca_bud = (uk == 0 && Lu == 0) ? -1 : ((ku < KT) ? C.cu_of_g[gu_bud] : -1);
            int u_t1 = (uk >= 3) ? C.cu_of_g[b3f_t1(&C, ku, hu)] : ca;
            int u_t2 = (uk >= 4) ? C.cu_of_g[b3f_t2(&C, ku, hu)] : ca;
            int u_t1b = (uk >= 3 && ku < KT) ? C.cu_of_g[b3f_t1(&C, ku + 1, hu)] : -1;
            int u_t2b = (uk >= 4 && ku < KT) ? C.cu_of_g[b3f_t2(&C, ku + 1, hu)] : -1;
            int uidx_t = (uk >= 2) ? b3f_tidx(&C, ku, hu) : 0;
            int64_t urch_q = (uk == 0) ? D.rch[gu] : D.trch[uidx_t];

            int u_yt  = (uk == 0 && hu == 1 && ku <= C.KYT && Du < C.DYT) ? C.yt_lut[ku][Lu] : -1;
            int u_ybu = (uk == 0 && ku <= C.KYB && Du < C.DYB) ? C.ybu_lut[hu][ku][Lu] : -1;
            int u_s3  = (uk == 0 && ku >= 3 && ku <= C.KS3 && Du < C.DS3) ? C.s3_lut[hu][ku][Lu] : -1;
            int ub_s3 = (uk == 0 && Lu >= 1 && ku + 1 >= 3 && ku + 1 <= C.KS3 && Du < C.DS3) ? C.s3_lut[hu][ku + 1][Lu - 1] : -1;
            int x_up = (u_yt >= 0), x_uz = (u_s3 >= 0);

            for (int i1 = 0; i1 < C.nchs; i1++) {
                int w1 = C.chs[i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1);
                int w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                int w1_ybw = (k1k == 0 && k1 >= 1 && k1 <= C.KYB && C.D[w1] < C.DYB) ? C.ybw_lut[k1][L1] : -1;
                int w1b_ybw = (k1k == 0 && L1 >= 1 && k1 + 1 <= C.KYB && C.D[w1] < C.DYB) ? C.ybw_lut[k1 + 1][L1 - 1] : -1;
                int x_w1 = (u_ybu >= 0 && w1_ybw >= 0);
                int64_t dg1_num = 0; int dg1_divD0 = 0;
                if (uk == 0) dg1_num = D.g[gup * N + w1] - D.g[gu * N + w1];
                else if (uk == 3) { dg1_num = D.g[C.INF_ID * N + w1] - D.g[b3f_t1(&C, ku, hu) * N + w1]; dg1_divD0 = 1; }
                else if (uk == 4) dg1_num = D.Yp_ref[b3f_yp_off(&C, uidx_t, w1)];

                for (int i2 = i1; i2 < C.nchs; i2++) {
                    int w2 = C.chs[i2], k2k = C.kind[w2], L2 = C.L[w2], k2 = C.K[w2];
                    int w2_up = b3f_next_up(&C, w2), w2_bud = b3f_next_bud(&C, w2);
                    int w2_t1 = (k2k >= 3) ? b3f_t1(&C, k2, 1) : w2, w2_t2 = (k2k >= 4) ? b3f_t2(&C, k2, 1) : w2;
                    int w2_t1b = (k2k >= 3 && k2 <= KT) ? b3f_t1(&C, k2 + 1, 1) : w2_t1, w2_t2b = (k2k >= 4 && k2 <= KT) ? b3f_t2(&C, k2 + 1, 1) : w2_t2;
                    int w2_ybw = (k2k == 0 && k2 >= 1 && k2 <= C.KYB && C.D[w2] < C.DYB) ? C.ybw_lut[k2][L2] : -1;
                    int w2b_ybw = (k2k == 0 && L2 >= 1 && k2 + 1 <= C.KYB && C.D[w2] < C.DYB) ? C.ybw_lut[k2 + 1][L2 - 1] : -1;
                    int x_w2 = (u_ybu >= 0 && w2_ybw >= 0);
                    int64_t dg2_num = 0; int dg2_divD0 = 0;
                    if (uk == 0) dg2_num = D.g[gup * N + w2] - D.g[gu * N + w2];
                    else if (uk == 3) { dg2_num = D.g[C.INF_ID * N + w2] - D.g[b3f_t1(&C, ku, hu) * N + w2]; dg2_divD0 = 1; }
                    else if (uk == 4) dg2_num = D.Yp_ref[b3f_yp_off(&C, uidx_t, w2)];

                    int64_t v = D.hC[b3f_hc_off(&C, ca, w1, w2)];
                    int64_t cu_up = (ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, w1, w2)] : 0, cu_bud = (ca_bud >= 0) ? D.hC[b3f_hc_off(&C, ca_bud, w1, w2)] : 0;
                    int64_t c1_up = D.hC[b3f_hc_off(&C, ca, w1_up, w2)], c2_up = D.hC[b3f_hc_off(&C, ca, w1, w2_up)];
                    int64_t c1_bud = D.hC[b3f_hc_off(&C, ca, w1_bud, w2)], c2_bud = D.hC[b3f_hc_off(&C, ca, w1, w2_bud)];
                    int64_t nb_u = (ub_s3 >= 0) ? D.qS3[b3f_s3_off(&C, ub_s3, w1, w2, iw_nb)] : 0;
                    int64_t nb_1 = (u_ybu >= 0 && w1b_ybw >= 0) ? D.qYb[b3f_yb_off(&C, u_ybu, w1b_ybw, w2, iw_nb)] : 0;
                    int64_t nb_2 = (u_ybu >= 0 && w2b_ybw >= 0) ? D.qYb[b3f_yb_off(&C, u_ybu, w2b_ybw, w1, iw_nb)] : 0;

                    zf_set_si(&W.acc, 0, 1);
                    zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);

                    if (dg2_num != 0) {
                        int64_t inum, iden; invD_rat(&C, w1, &inum, &iden);
                        mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                        mpz_set_si(W.t4, urch_q); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                        mpz_mul_si(W.t3, W.t3, dg2_num);
                        mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r);
                        if (dg2_divD0) mpz_mul_ui(W.t4, W.t4, (uint64_t)D0);
                        zf_add_z(&W.acc, W.t3, W.t4, &W);
                    }
                    if (dg1_num != 0) {
                        int64_t inum, iden; invD_rat(&C, w2, &inum, &iden);
                        mpz_set_si(W.t3, inum); mpz_mul_si(W.t3, W.t3, S_r);
                        mpz_set_si(W.t4, urch_q); mpz_submul_ui(W.t3, W.t4, (uint64_t)iden);
                        mpz_mul_si(W.t3, W.t3, dg1_num);
                        mpz_set_ui(W.t4, (uint64_t)iden); mpz_mul_ui(W.t4, W.t4, (uint64_t)S_r);
                        if (dg1_divD0) mpz_mul_ui(W.t4, W.t4, (uint64_t)D0);
                        zf_add_z(&W.acc, W.t3, W.t4, &W);
                    }

                    if (k1k == 0) {
                        int64_t d1 = c1_up - v;
                        zf_add_mul_si(&W.acc, (int64_t)L1, d1, 1, &W);
                        if (uk == 0) zf_add_si(&W.acc, d1, (uint64_t)Du, &W);
                        if (x_w1) zf_add_mul_si(&W.acc, (int64_t)k1 * D.rch[w1], d1, (uint64_t)S_r, &W);
                        if (L1 >= 1) zf_add_mul_si(&W.acc, (int64_t)L1, c1_bud + nb_1 - v, (uint64_t)C.D[w1], &W);
                        if (k1 >= 1 && !x_w1 && ((mode == 1 && d1 > 0) || (mode == 2 && d1 < 0))) zf_add_mul_si(&W.acc, (int64_t)k1, d1, 2, &W);
                    }
                    if (k2k == 0) {
                        int64_t d2 = c2_up - v;
                        zf_add_mul_si(&W.acc, (int64_t)L2, d2, 1, &W);
                        if (uk == 0) zf_add_si(&W.acc, d2, (uint64_t)Du, &W);
                        if (x_w2) zf_add_mul_si(&W.acc, (int64_t)k2 * D.rch[w2], d2, (uint64_t)S_r, &W);
                        if (L2 >= 1) zf_add_mul_si(&W.acc, (int64_t)L2, c2_bud + nb_2 - v, (uint64_t)C.D[w2], &W);
                        if (k2 >= 1 && !x_w2 && ((mode == 1 && d2 > 0) || (mode == 2 && d2 < 0))) zf_add_mul_si(&W.acc, (int64_t)k2, d2, 2, &W);
                    }
                    if (uk == 0) {
                        int64_t du = cu_up - v;
                        zf_add_mul_si(&W.acc, (int64_t)Lu, du, 1, &W);
                        if (k1k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w1], &W);
                        if (k2k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w2], &W);
                        int64_t rq = (x_up ? D.rpar[gu] : 0) + (x_uz ? (int64_t)(ku - 2) * D.rch[gu] : 0);
                        if (rq != 0) zf_add_mul_si(&W.acc, rq, du, (uint64_t)S_r, &W);
                        if (Lu >= 1) zf_add_mul_si(&W.acc, (int64_t)Lu, cu_bud + nb_u - v, (uint64_t)Du, &W);
                        int rh = (hu && !x_up ? 1 : 0) + (ku >= 3 && !x_uz ? (ku - 2) : 0);
                        if (rh > 0 && ((mode == 1 && du > 0) || (mode == 2 && du < 0))) zf_add_mul_si(&W.acc, (int64_t)rh, du, 2, &W);
                    }

                    if (k1k >= 2) {
                        int64_t ln, ld, hn, hd; get_rin1(&C, gu, &ln, &ld, &hn, &hd);
                        int64_t f10 = D.hC[b3f_hc_off(&C, ca, C.INF_ID, w2)], f11 = (k1k >= 3) ? D.hC[b3f_hc_off(&C, ca, w1_t1, w2)] : v, f12 = (k1k >= 4) ? D.hC[b3f_hc_off(&C, ca, w1_t2, w2)] : v;
                        int64_t f11b = (k1k >= 3) ? D.hC[b3f_hc_off(&C, ca, w1_t1b, w2)] : c1_bud, f12b = (k1k >= 4) ? D.hC[b3f_hc_off(&C, ca, w1_t2b, w2)] : c1_bud;
                        int64_t d10 = (uk == 0) ? (((ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, C.INF_ID, w2)] : 0) - f10) : 0;
                        int64_t d11 = (k1k >= 3 && uk == 0) ? (((ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, w1_t1, w2)] : 0) - f11) : 0;
                        int64_t d12 = (k1k >= 4 && uk == 0) ? (((ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, w1_t2, w2)] : 0) - f12) : 0;
                        add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k1k, k1, 1, ln, ld, hn, hd, k1, 0, 0, 0,
                                     f10, f11, f12, v, (k1k == 2) ? c1_bud : f11b, (k1k == 3) ? c1_bud : f12b, c1_bud, uk == 0, d10, d11, d12, cu_up);
                    }
                    if (k2k >= 2) {
                        int64_t ln, ld, hn, hd; get_rin1(&C, gu, &ln, &ld, &hn, &hd);
                        int64_t f20 = D.hC[b3f_hc_off(&C, ca, w1, C.INF_ID)], f21 = (k2k >= 3) ? D.hC[b3f_hc_off(&C, ca, w1, w2_t1)] : v, f22 = (k2k >= 4) ? D.hC[b3f_hc_off(&C, ca, w1, w2_t2)] : v;
                        int64_t f21b = (k2k >= 3) ? D.hC[b3f_hc_off(&C, ca, w1, w2_t1b)] : c2_bud, f22b = (k2k >= 4) ? D.hC[b3f_hc_off(&C, ca, w1, w2_t2b)] : c2_bud;
                        int64_t d20 = (uk == 0) ? (((ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, w1, C.INF_ID)] : 0) - f20) : 0;
                        int64_t d21 = (k2k >= 3 && uk == 0) ? (((ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, w1, w2_t1)] : 0) - f21) : 0;
                        int64_t d22 = (k2k >= 4 && uk == 0) ? (((ca_up >= 0) ? D.hC[b3f_hc_off(&C, ca_up, w1, w2_t2)] : 0) - f22) : 0;
                        add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k2k, k2, 1, ln, ld, hn, hd, k2, 0, 0, 0,
                                     f20, f21, f22, v, (k2k == 2) ? c2_bud : f21b, (k2k == 3) ? c2_bud : f22b, c2_bud, uk == 0, d20, d21, d22, cu_up);
                    }
                    if (uk >= 2) {
                        int64_t ln, ld, hn, hd; get_rin2(&C, w1, w2, &ln, &ld, &hn, &hd);
                        int64_t uF1 = (uk >= 3) ? D.hC[b3f_hc_off(&C, u_t1, w1, w2)] : v, uF2 = (uk >= 4) ? D.hC[b3f_hc_off(&C, u_t2, w1, w2)] : v;
                        int64_t uF1b = (uk == 2) ? cu_bud : ((uk >= 3 && u_t1b >= 0) ? D.hC[b3f_hc_off(&C, u_t1b, w1, w2)] : 0);
                        int64_t uF2b = (uk == 3) ? cu_bud : ((uk >= 4 && u_t2b >= 0) ? D.hC[b3f_hc_off(&C, u_t2b, w1, w2)] : 0);
                        int u_nin = (k1k == 0 ? 1 : 0) + (k2k == 0 ? 1 : 0);
                        int64_t ud1 = (uk >= 3 && u_nin > 0) ? ((k1k == 0 ? (D.hC[b3f_hc_off(&C, u_t1, w1_up, w2)] - uF1) : 0) + (k2k == 0 ? (D.hC[b3f_hc_off(&C, u_t1, w1, w2_up)] - uF1) : 0)) : 0;
                        int64_t ud2 = (uk >= 4 && u_nin > 0) ? ((k1k == 0 ? (D.hC[b3f_hc_off(&C, u_t2, w1_up, w2)] - uF2) : 0) + (k2k == 0 ? (D.hC[b3f_hc_off(&C, u_t2, w1, w2_up)] - uF2) : 0)) : 0;
                        int64_t uF3_up = (uk == 4 && u_nin > 0) ? ((k1k == 0 ? c1_up : 0) + (k2k == 0 ? c2_up : 0)) : 0;
                        add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, uk, ku, hu, ln, ld, hn, hd, hu + ku - 2, 0, 0, 0,
                                     0, uF1, uF2, v, uF1b, uF2b, cu_bud, u_nin, 0, ud1, ud2, uF3_up);
                    }

                    if (zf_fail(&W.acc, mode)) fails_hc++;
                }
            }
        }
        gw_clear(&W);
    }
    printf("  [hC] checked %zu states: fails = %lld\n", C.hc_len, fails_hc); fflush(stdout);

    /* 5a. Verify P4 */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_p4)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int b4 = 0; b4 < C.NP4; b4++) {
            int Lb = C.p4_L[b4], kb = C.p4_K[b4], Db = Lb + kb + 1, b_hu = C.p4_hu[b4], b_g = C.p4_g[b4], b_g_up = b3f_next_up(&C, b_g);
            int b4_up = (Db + 1 < C.DP4) ? C.p4_lut[kb][Lb + 1] : -1, b4_bud = (Lb >= 1 && kb + 1 <= C.KP4) ? C.p4_lut[kb + 1][Lb - 1] : -1;
            for (int c4 = 0; c4 < C.NP4; c4++) {
                int Lc = C.p4_L[c4], kc = C.p4_K[c4], Dc = Lc + kc + 1, c_hu = C.p4_hu[c4], c_g = C.p4_g[c4], c_g_up = b3f_next_up(&C, c_g);
                int c4_up = (Dc + 1 < C.DP4) ? C.p4_lut[kc][Lc + 1] : -1, c4_bud = (Lc >= 1 && kc + 1 <= C.KP4) ? C.p4_lut[kc + 1][Lc - 1] : -1;
                for (int ia = 0; ia < C.npars; ia++) {
                    int a = C.pars[ia], ak = C.kind[a], La = C.L[a], ka = C.K[a], ha = C.H[a];
                    int a_up = b3f_next_up(&C, a), a_bud = b3f_next_bud(&C, a), a_t1 = (ak >= 3) ? b3f_t1(&C, ka, ha) : a, a_t2 = (ak >= 4) ? b3f_t2(&C, ka, ha) : a;
                    int a_t1b = (ak >= 3 && ka <= KT) ? b3f_t1(&C, ka + 1, ha) : a_t1, a_t2b = (ak >= 4 && ka <= KT) ? b3f_t2(&C, ka + 1, ha) : a_t2;
                    int64_t dh_abc = D.hV[b3f_hv_off(&C, a, b_hu, c_g_up)] - D.hV[b3f_hv_off(&C, a, b_hu, c_g)];
                    for (int id = 0; id < C.nchs; id++) {
                        int d = C.chs[id], dk = C.kind[d], Ld = C.L[d], kd = C.K[d];
                        int d_up = b3f_next_up(&C, d), d_bud = b3f_next_bud(&C, d), d_t1 = (dk >= 3) ? b3f_t1(&C, kd, 1) : d, d_t2 = (dk >= 4) ? b3f_t2(&C, kd, 1) : d;
                        int d_t1b = (dk >= 3 && kd <= KT) ? b3f_t1(&C, kd + 1, 1) : d_t1, d_t2b = (dk >= 4 && kd <= KT) ? b3f_t2(&C, kd + 1, 1) : d_t2;
                        int64_t dh_bcd = D.hV[b3f_hv_off(&C, b_g_up, c_hu, d)] - D.hV[b3f_hv_off(&C, b_g, c_hu, d)];
                        int64_t v = D.qP4[b3f_p4_off(&C, a, b4, c4, d)];
                        int64_t qa_up = D.qP4[b3f_p4_off(&C, a_up, b4, c4, d)], qb_up = (b4_up >= 0) ? D.qP4[b3f_p4_off(&C, a, b4_up, c4, d)] : 0;
                        int64_t qc_up = (c4_up >= 0) ? D.qP4[b3f_p4_off(&C, a, b4, c4_up, d)] : 0, qd_up = D.qP4[b3f_p4_off(&C, a, b4, c4, d_up)];
                        int64_t qa_b = D.qP4[b3f_p4_off(&C, a_bud, b4, c4, d)], qb_b = (b4_bud >= 0) ? D.qP4[b3f_p4_off(&C, a, b4_bud, c4, d)] : 0;
                        int64_t qc_b = (c4_bud >= 0) ? D.qP4[b3f_p4_off(&C, a, b4, c4_bud, d)] : 0, qd_b = D.qP4[b3f_p4_off(&C, a, b4, c4, d_bud)];
                        size_t bi[4]; int ed[4];
                        l5_couple(x,(L5Tree){4,{-1,0,1,2},{a,b_g,c_g,d}},iw_nb,bi,ed);
                        qa_b += newborn_int(&D,bi[0]); qb_b += newborn_int(&D,bi[1]);
                        qc_b += newborn_int(&D,bi[2]); qd_b += newborn_int(&D,bi[3]);

                        zf_set_si(&W.acc, 0, 1);
                        zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);
                        add_centered_source(&W.acc, &C, a, D.rpar[b_g], S_r, dh_bcd, &W);
                        add_centered_source(&W.acc, &C, d, D.rch[c_g],  S_r, dh_abc, &W);

                        if (ak == 0) {
                            int64_t da = qa_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)La, da, 1, &W); zf_add_si(&W.acc, da, (uint64_t)Db, &W);
                            if (La >= 1) zf_add_mul_si(&W.acc, (int64_t)La, qa_b - v, (uint64_t)C.D[a], &W);
                            if (ha + ka - 1 > 0 && ((mode == 1 && da > 0) || (mode == 2 && da < 0))) zf_add_mul_si(&W.acc, (int64_t)(ha + ka - 1), da, ed[0], &W);
                        }
                        if (dk == 0) {
                            int64_t dd = qd_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Ld, dd, 1, &W); zf_add_si(&W.acc, dd, (uint64_t)Dc, &W);
                            if (Ld >= 1) zf_add_mul_si(&W.acc, (int64_t)Ld, qd_b - v, (uint64_t)C.D[d], &W);
                            if (kd > 0 && ((mode == 1 && dd > 0) || (mode == 2 && dd < 0))) zf_add_mul_si(&W.acc, (int64_t)kd, dd, ed[3], &W);
                        }
                        {
                            int64_t db = qb_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lb, db, 1, &W); zf_add_si(&W.acc, db, (uint64_t)Dc, &W);
                            if (ak == 0) zf_add_si(&W.acc, db, (uint64_t)C.D[a], &W);
                            if (Lb >= 1) zf_add_mul_si(&W.acc, (int64_t)Lb, qb_b - v, (uint64_t)Db, &W);
                            if (kb - 1 > 0 && ((mode == 1 && db > 0) || (mode == 2 && db < 0))) zf_add_mul_si(&W.acc, (int64_t)(kb - 1), db, ed[1], &W);
                        }
                        {
                            int64_t dc = qc_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lc, dc, 1, &W); zf_add_si(&W.acc, dc, (uint64_t)Db, &W);
                            if (dk == 0) zf_add_si(&W.acc, dc, (uint64_t)C.D[d], &W);
                            if (Lc >= 1) zf_add_mul_si(&W.acc, (int64_t)Lc, qc_b - v, (uint64_t)Dc, &W);
                            if (kc - 1 > 0 && ((mode == 1 && dc > 0) || (mode == 2 && dc < 0))) zf_add_mul_si(&W.acc, (int64_t)(kc - 1), dc, ed[2], &W);
                        }
                        if (ak >= 2) {
                            int64_t aF0 = D.qP4[b3f_p4_off(&C, C.INF_ID, b4, c4, d)], aF1 = (ak >= 3) ? D.qP4[b3f_p4_off(&C, a_t1, b4, c4, d)] : v, aF2 = (ak >= 4) ? D.qP4[b3f_p4_off(&C, a_t2, b4, c4, d)] : v;
                            int64_t aF1b = (ak >= 3) ? D.qP4[b3f_p4_off(&C, a_t1b, b4, c4, d)] : qa_b, aF2b = (ak >= 4) ? D.qP4[b3f_p4_off(&C, a_t2b, b4, c4, d)] : qa_b;
                            int64_t ad0 = ((b4_up >= 0) ? D.qP4[b3f_p4_off(&C, C.INF_ID, b4_up, c4, d)] : 0) - aF0;
                            int64_t ad1 = (ak >= 3) ? (((b4_up >= 0) ? D.qP4[b3f_p4_off(&C, a_t1, b4_up, c4, d)] : 0) - aF1) : 0;
                            int64_t ad2 = (ak >= 4) ? (((b4_up >= 0) ? D.qP4[b3f_p4_off(&C, a_t2, b4_up, c4, d)] : 0) - aF2) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, ak, ka, ha, 1, Db, 1, Db, ka + ha - 1, 0, 0, 0,
                                         aF0, aF1, aF2, v, (ak == 2) ? qa_b : aF1b, (ak == 3) ? qa_b : aF2b, qa_b, 1, ad0, ad1, ad2, qb_up);
                        }
                        if (dk >= 2) {
                            int64_t dF0 = D.qP4[b3f_p4_off(&C, a, b4, c4, C.INF_ID)], dF1 = (dk >= 3) ? D.qP4[b3f_p4_off(&C, a, b4, c4, d_t1)] : v, dF2 = (dk >= 4) ? D.qP4[b3f_p4_off(&C, a, b4, c4, d_t2)] : v;
                            int64_t dF1b = (dk >= 3) ? D.qP4[b3f_p4_off(&C, a, b4, c4, d_t1b)] : qd_b, dF2b = (dk >= 4) ? D.qP4[b3f_p4_off(&C, a, b4, c4, d_t2b)] : qd_b;
                            int64_t dd0 = ((c4_up >= 0) ? D.qP4[b3f_p4_off(&C, a, b4, c4_up, C.INF_ID)] : 0) - dF0;
                            int64_t dd1 = (dk >= 3) ? (((c4_up >= 0) ? D.qP4[b3f_p4_off(&C, a, b4, c4_up, d_t1)] : 0) - dF1) : 0;
                            int64_t dd2 = (dk >= 4) ? (((c4_up >= 0) ? D.qP4[b3f_p4_off(&C, a, b4, c4_up, d_t2)] : 0) - dF2) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, dk, kd, 1, 1, Dc, 1, Dc, kd, 0, 0, 0,
                                         dF0, dF1, dF2, v, (dk == 2) ? qd_b : dF1b, (dk == 3) ? qd_b : dF2b, qd_b, 1, dd0, dd1, dd2, qc_up);
                        }
                        if (zf_fail(&W.acc, mode)) fails_p4++;
                    }
                }
            }
        }
        gw_clear(&W);
    }
    printf("  [P4] checked %zu states: fails = %lld\n", C.p4_len, fails_p4); fflush(stdout);

    /* 5b. Verify Ytop */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_yt)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int u4 = 0; u4 < C.NYT; u4++) {
            int Lu = C.yt_L[u4], ku = C.yt_K[u4], Du = Lu + ku + 1, u_hu = C.yt_hu[u4], u_cu = C.yt_cu[u4], u_g = C.yt_g[u4], u_g_up = b3f_next_up(&C, u_g);
            int u_hu_up = C.hu_of_g[u_g_up], u_cu_up = C.cu_of_g[u_g_up];
            int u4_up = (Du + 1 < C.DYT) ? C.yt_lut[ku][Lu + 1] : -1, u4_bud = (Lu >= 1 && ku + 1 <= C.KYT) ? C.yt_lut[ku + 1][Lu - 1] : -1;
            for (int ip = 0; ip < C.npars; ip++) {
                int p = C.pars[ip], pk = C.kind[p], Lp = C.L[p], kp = C.K[p], hp = C.H[p];
                int p_up = b3f_next_up(&C, p), p_bud = b3f_next_bud(&C, p), p_t1 = (pk >= 3) ? b3f_t1(&C, kp, hp) : p, p_t2 = (pk >= 4) ? b3f_t2(&C, kp, hp) : p;
                int p_t1b = (pk >= 3 && kp <= KT) ? b3f_t1(&C, kp + 1, hp) : p_t1, p_t2b = (pk >= 4 && kp <= KT) ? b3f_t2(&C, kp + 1, hp) : p_t2;
                for (int i1 = 0; i1 < C.nchs; i1++) {
                    int w1 = C.chs[i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                    int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1), w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                    int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                    int64_t dhv1 = ((u_hu_up >= 0) ? D.hV[b3f_hv_off(&C, p, u_hu_up, w1)] : 0) - D.hV[b3f_hv_off(&C, p, u_hu, w1)];
                    for (int i2 = i1; i2 < C.nchs; i2++) {
                        int w2 = C.chs[i2], k2k = C.kind[w2], L2 = C.L[w2], k2 = C.K[w2];
                        int w2_up = b3f_next_up(&C, w2), w2_bud = b3f_next_bud(&C, w2), w2_t1 = (k2k >= 3) ? b3f_t1(&C, k2, 1) : w2, w2_t2 = (k2k >= 4) ? b3f_t2(&C, k2, 1) : w2;
                        int w2_t1b = (k2k >= 3 && k2 <= KT) ? b3f_t1(&C, k2 + 1, 1) : w2_t1, w2_t2b = (k2k >= 4 && k2 <= KT) ? b3f_t2(&C, k2 + 1, 1) : w2_t2;
                        int64_t dhv2 = ((u_hu_up >= 0) ? D.hV[b3f_hv_off(&C, p, u_hu_up, w2)] : 0) - D.hV[b3f_hv_off(&C, p, u_hu, w2)];
                        int64_t dhc  = ((u_cu_up >= 0) ? D.hC[b3f_hc_off(&C, u_cu_up, w1, w2)] : 0) - D.hC[b3f_hc_off(&C, u_cu, w1, w2)];
                        int64_t v = D.qYt[b3f_yt_off(&C, p, u4, w1, w2)];
                        int64_t qp_up = D.qYt[b3f_yt_off(&C, p_up, u4, w1, w2)], qu_up = (u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, w1, w2)] : 0;
                        int64_t q1_up = D.qYt[b3f_yt_off(&C, p, u4, w1_up, w2)], q2_up = D.qYt[b3f_yt_off(&C, p, u4, w1, w2_up)];
                        int64_t qp_b = D.qYt[b3f_yt_off(&C, p_bud, u4, w1, w2)], qu_b = (u4_bud >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_bud, w1, w2)] : 0;
                        int64_t q1_b = D.qYt[b3f_yt_off(&C, p, u4, w1_bud, w2)], q2_b = D.qYt[b3f_yt_off(&C, p, u4, w1, w2_bud)];
                        size_t bi[4]; int ed[4];
                        l5_couple(x,(L5Tree){4,{-1,0,1,1},{p,u_g,w1,w2}},iw_nb,bi,ed);
                        qp_b += newborn_int(&D,bi[0]); qu_b += newborn_int(&D,bi[1]);
                        q1_b += newborn_int(&D,bi[2]); q2_b += newborn_int(&D,bi[3]);

                        zf_set_si(&W.acc, 0, 1);
                        zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);
                        add_centered_source(&W.acc, &C, p,  D.rpar[u_g], S_r, dhc,  &W);
                        add_centered_source(&W.acc, &C, w2, D.rch[u_g],  S_r, dhv1, &W);
                        add_centered_source(&W.acc, &C, w1, D.rch[u_g],  S_r, dhv2, &W);

                        if (pk == 0) {
                            int64_t dp = qp_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lp, dp, 1, &W); zf_add_si(&W.acc, dp, (uint64_t)Du, &W);
                            if (Lp >= 1) zf_add_mul_si(&W.acc, (int64_t)Lp, qp_b - v, (uint64_t)C.D[p], &W);
                            if (hp + kp - 1 > 0 && ((mode == 1 && dp > 0) || (mode == 2 && dp < 0))) zf_add_mul_si(&W.acc, (int64_t)(hp + kp - 1), dp, ed[0], &W);
                        }
                        if (k1k == 0) {
                            int64_t d1 = q1_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L1, d1, 1, &W); zf_add_si(&W.acc, d1, (uint64_t)Du, &W);
                            if (L1 >= 1) zf_add_mul_si(&W.acc, (int64_t)L1, q1_b - v, (uint64_t)C.D[w1], &W);
                            if (k1 > 0 && ((mode == 1 && d1 > 0) || (mode == 2 && d1 < 0))) zf_add_mul_si(&W.acc, (int64_t)k1, d1, ed[2], &W);
                        }
                        if (k2k == 0) {
                            int64_t d2 = q2_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L2, d2, 1, &W); zf_add_si(&W.acc, d2, (uint64_t)Du, &W);
                            if (L2 >= 1) zf_add_mul_si(&W.acc, (int64_t)L2, q2_b - v, (uint64_t)C.D[w2], &W);
                            if (k2 > 0 && ((mode == 1 && d2 > 0) || (mode == 2 && d2 < 0))) zf_add_mul_si(&W.acc, (int64_t)k2, d2, ed[3], &W);
                        }
                        {
                            int64_t du = qu_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lu, du, 1, &W);
                            if (pk == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[p], &W);
                            if (k1k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w1], &W);
                            if (k2k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w2], &W);
                            if (Lu >= 1) zf_add_mul_si(&W.acc, (int64_t)Lu, qu_b - v, (uint64_t)Du, &W);
                            if (ku - 2 > 0 && ((mode == 1 && du > 0) || (mode == 2 && du < 0))) zf_add_mul_si(&W.acc, (int64_t)(ku - 2), du, ed[1], &W);
                        }
                        if (pk >= 2) {
                            int64_t pF0 = D.qYt[b3f_yt_off(&C, C.INF_ID, u4, w1, w2)], pF1 = (pk >= 3) ? D.qYt[b3f_yt_off(&C, p_t1, u4, w1, w2)] : v, pF2 = (pk >= 4) ? D.qYt[b3f_yt_off(&C, p_t2, u4, w1, w2)] : v;
                            int64_t pF1b = (pk >= 3) ? D.qYt[b3f_yt_off(&C, p_t1b, u4, w1, w2)] : qp_b, pF2b = (pk >= 4) ? D.qYt[b3f_yt_off(&C, p_t2b, u4, w1, w2)] : qp_b;
                            int64_t pd0 = ((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, C.INF_ID, u4_up, w1, w2)] : 0) - pF0;
                            int64_t pd1 = (pk >= 3) ? (((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p_t1, u4_up, w1, w2)] : 0) - pF1) : 0;
                            int64_t pd2 = (pk >= 4) ? (((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p_t2, u4_up, w1, w2)] : 0) - pF2) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, pk, kp, hp, 1, Du, 1, Du, kp + hp - 1, 0, 0, 0,
                                         pF0, pF1, pF2, v, (pk == 2) ? qp_b : pF1b, (pk == 3) ? qp_b : pF2b, qp_b, 1, pd0, pd1, pd2, qu_up);
                        }
                        if (k1k >= 2) {
                            int64_t f10 = D.qYt[b3f_yt_off(&C, p, u4, C.INF_ID, w2)], f11 = (k1k >= 3) ? D.qYt[b3f_yt_off(&C, p, u4, w1_t1, w2)] : v, f12 = (k1k >= 4) ? D.qYt[b3f_yt_off(&C, p, u4, w1_t2, w2)] : v;
                            int64_t f11b = (k1k >= 3) ? D.qYt[b3f_yt_off(&C, p, u4, w1_t1b, w2)] : q1_b, f12b = (k1k >= 4) ? D.qYt[b3f_yt_off(&C, p, u4, w1_t2b, w2)] : q1_b;
                            int64_t d10 = ((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, C.INF_ID, w2)] : 0) - f10;
                            int64_t d11 = (k1k >= 3) ? (((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, w1_t1, w2)] : 0) - f11) : 0;
                            int64_t d12 = (k1k >= 4) ? (((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, w1_t2, w2)] : 0) - f12) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k1k, k1, 1, 1, Du, 1, Du, k1, 0, 0, 0,
                                         f10, f11, f12, v, (k1k == 2) ? q1_b : f11b, (k1k == 3) ? q1_b : f12b, q1_b, 1, d10, d11, d12, qu_up);
                        }
                        if (k2k >= 2) {
                            int64_t f20 = D.qYt[b3f_yt_off(&C, p, u4, w1, C.INF_ID)], f21 = (k2k >= 3) ? D.qYt[b3f_yt_off(&C, p, u4, w1, w2_t1)] : v, f22 = (k2k >= 4) ? D.qYt[b3f_yt_off(&C, p, u4, w1, w2_t2)] : v;
                            int64_t f21b = (k2k >= 3) ? D.qYt[b3f_yt_off(&C, p, u4, w1, w2_t1b)] : q2_b, f22b = (k2k >= 4) ? D.qYt[b3f_yt_off(&C, p, u4, w1, w2_t2b)] : q2_b;
                            int64_t d20 = ((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, w1, C.INF_ID)] : 0) - f20;
                            int64_t d21 = (k2k >= 3) ? (((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, w1, w2_t1)] : 0) - f21) : 0;
                            int64_t d22 = (k2k >= 4) ? (((u4_up >= 0) ? D.qYt[b3f_yt_off(&C, p, u4_up, w1, w2_t2)] : 0) - f22) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k2k, k2, 1, 1, Du, 1, Du, k2, 0, 0, 0,
                                         f20, f21, f22, v, (k2k == 2) ? q2_b : f21b, (k2k == 3) ? q2_b : f22b, q2_b, 1, d20, d21, d22, qu_up);
                        }
                        if (zf_fail(&W.acc, mode)) fails_yt++;
                    }
                }
            }
        }
        gw_clear(&W);
    }
    printf("  [Yt] checked %zu states: fails = %lld\n", C.yt_len, fails_yt); fflush(stdout);

    /* 5c. Verify Ybot */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_yb)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int u4 = 0; u4 < C.NYBU; u4++) {
            int Lu = C.ybu_L[u4], ku = C.ybu_K[u4], hu = C.ybu_H[u4], Du = Lu + ku + hu, u_cu = C.ybu_cu[u4], u_g = C.ybu_g[u4], u_g_up = b3f_next_up(&C, u_g);
            int u4_up = (Du + 1 < C.DYB) ? C.ybu_lut[hu][ku][Lu + 1] : -1, u4_bud = (Lu >= 1 && ku + 1 <= C.KYB) ? C.ybu_lut[hu][ku + 1][Lu - 1] : -1;
            for (int w2_4 = 0; w2_4 < C.NYBW; w2_4++) {
                int L2 = C.ybw_L[w2_4], k2 = C.ybw_K[w2_4], D2 = L2 + k2 + 1, w2_hu = C.ybw_hu[w2_4], w2_g = C.ybw_g[w2_4], w2_g_up = b3f_next_up(&C, w2_g);
                int w2_up = (D2 + 1 < C.DYB) ? C.ybw_lut[k2][L2 + 1] : -1, w2_bud = (L2 >= 1 && k2 + 1 <= C.KYB) ? C.ybw_lut[k2 + 1][L2 - 1] : -1;
                for (int i1 = 0; i1 < C.nchs; i1++) {
                    int w1 = C.chs[i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                    int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1), w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                    int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                    int64_t dhc = D.hC[b3f_hc_off(&C, u_cu, w1, w2_g_up)] - D.hC[b3f_hc_off(&C, u_cu, w1, w2_g)];
                    for (int ic = 0; ic < C.nchs; ic++) {
                        int c = C.chs[ic], ck = C.kind[c], Lc = C.L[c], kc = C.K[c];
                        int c_up = b3f_next_up(&C, c), c_bud = b3f_next_bud(&C, c), c_t1 = (ck >= 3) ? b3f_t1(&C, kc, 1) : c, c_t2 = (ck >= 4) ? b3f_t2(&C, kc, 1) : c;
                        int c_t1b = (ck >= 3 && kc <= KT) ? b3f_t1(&C, kc + 1, 1) : c_t1, c_t2b = (ck >= 4 && kc <= KT) ? b3f_t2(&C, kc + 1, 1) : c_t2;
                        int64_t dhv = D.hV[b3f_hv_off(&C, u_g_up, w2_hu, c)] - D.hV[b3f_hv_off(&C, u_g, w2_hu, c)];
                        int64_t v = D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c)];
                        int64_t qu_up = (u4_up >= 0) ? D.qYb[b3f_yb_off(&C, u4_up, w2_4, w1, c)] : 0, q2_up = (w2_up >= 0) ? D.qYb[b3f_yb_off(&C, u4, w2_up, w1, c)] : 0;
                        int64_t q1_up = D.qYb[b3f_yb_off(&C, u4, w2_4, w1_up, c)], qc_up = D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c_up)];
                        int64_t qu_b = (u4_bud >= 0) ? D.qYb[b3f_yb_off(&C, u4_bud, w2_4, w1, c)] : 0, q2_b = (w2_bud >= 0) ? D.qYb[b3f_yb_off(&C, u4, w2_bud, w1, c)] : 0;
                        int64_t q1_b = D.qYb[b3f_yb_off(&C, u4, w2_4, w1_bud, c)], qc_b = D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c_bud)];
                        size_t bi[4]; int ed[4];
                        l5_couple(x,(L5Tree){4,{-1,0,0,2},{u_g,w1,w2_g,c}},iw_nb,bi,ed);
                        qu_b += newborn_int(&D,bi[0]); q1_b += newborn_int(&D,bi[1]);
                        q2_b += newborn_int(&D,bi[2]); qc_b += newborn_int(&D,bi[3]);

                        zf_set_si(&W.acc, 0, 1);
                        zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);
                        add_centered_source(&W.acc, &C, w1, D.rch[u_g],  S_r, dhv, &W);
                        add_centered_source(&W.acc, &C, c,  D.rch[w2_g], S_r, dhc, &W);

                        if (k1k == 0) {
                            int64_t d1 = q1_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L1, d1, 1, &W); zf_add_si(&W.acc, d1, (uint64_t)Du, &W);
                            if (L1 >= 1) zf_add_mul_si(&W.acc, (int64_t)L1, q1_b - v, (uint64_t)C.D[w1], &W);
                            if (k1 > 0 && ((mode == 1 && d1 > 0) || (mode == 2 && d1 < 0))) zf_add_mul_si(&W.acc, (int64_t)k1, d1, ed[1], &W);
                        }
                        if (ck == 0) {
                            int64_t dc = qc_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lc, dc, 1, &W); zf_add_si(&W.acc, dc, (uint64_t)D2, &W);
                            if (Lc >= 1) zf_add_mul_si(&W.acc, (int64_t)Lc, qc_b - v, (uint64_t)C.D[c], &W);
                            if (kc > 0 && ((mode == 1 && dc > 0) || (mode == 2 && dc < 0))) zf_add_mul_si(&W.acc, (int64_t)kc, dc, ed[3], &W);
                        }
                        {
                            int64_t du = qu_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lu, du, 1, &W); zf_add_si(&W.acc, du, (uint64_t)D2, &W);
                            if (k1k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w1], &W);
                            if (Lu >= 1) zf_add_mul_si(&W.acc, (int64_t)Lu, qu_b - v, (uint64_t)Du, &W);
                            if (hu + ku - 2 > 0 && ((mode == 1 && du > 0) || (mode == 2 && du < 0))) zf_add_mul_si(&W.acc, (int64_t)(hu + ku - 2), du, ed[0], &W);
                        }
                        {
                            int64_t d2 = q2_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L2, d2, 1, &W); zf_add_si(&W.acc, d2, (uint64_t)Du, &W);
                            if (ck == 0) zf_add_si(&W.acc, d2, (uint64_t)C.D[c], &W);
                            if (L2 >= 1) zf_add_mul_si(&W.acc, (int64_t)L2, q2_b - v, (uint64_t)D2, &W);
                            if (k2 - 1 > 0 && ((mode == 1 && d2 > 0) || (mode == 2 && d2 < 0))) zf_add_mul_si(&W.acc, (int64_t)(k2 - 1), d2, ed[2], &W);
                        }
                        if (k1k >= 2) {
                            int64_t f10 = D.qYb[b3f_yb_off(&C, u4, w2_4, C.INF_ID, c)], f11 = (k1k >= 3) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1_t1, c)] : v, f12 = (k1k >= 4) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1_t2, c)] : v;
                            int64_t f11b = (k1k >= 3) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1_t1b, c)] : q1_b, f12b = (k1k >= 4) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1_t2b, c)] : q1_b;
                            int64_t d10 = ((u4_up >= 0) ? D.qYb[b3f_yb_off(&C, u4_up, w2_4, C.INF_ID, c)] : 0) - f10;
                            int64_t d11 = (k1k >= 3) ? (((u4_up >= 0) ? D.qYb[b3f_yb_off(&C, u4_up, w2_4, w1_t1, c)] : 0) - f11) : 0;
                            int64_t d12 = (k1k >= 4) ? (((u4_up >= 0) ? D.qYb[b3f_yb_off(&C, u4_up, w2_4, w1_t2, c)] : 0) - f12) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k1k, k1, 1, 1, Du, 1, Du, k1, 0, 0, 0,
                                         f10, f11, f12, v, (k1k == 2) ? q1_b : f11b, (k1k == 3) ? q1_b : f12b, q1_b, 1, d10, d11, d12, qu_up);
                        }
                        if (ck >= 2) {
                            int64_t fc0 = D.qYb[b3f_yb_off(&C, u4, w2_4, w1, C.INF_ID)], fc1 = (ck >= 3) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t1)] : v, fc2 = (ck >= 4) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t2)] : v;
                            int64_t fc1b = (ck >= 3) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t1b)] : qc_b, fc2b = (ck >= 4) ? D.qYb[b3f_yb_off(&C, u4, w2_4, w1, c_t2b)] : qc_b;
                            int64_t dc0 = ((w2_up >= 0) ? D.qYb[b3f_yb_off(&C, u4, w2_up, w1, C.INF_ID)] : 0) - fc0;
                            int64_t dc1 = (ck >= 3) ? (((w2_up >= 0) ? D.qYb[b3f_yb_off(&C, u4, w2_up, w1, c_t1)] : 0) - fc1) : 0;
                            int64_t dc2 = (ck >= 4) ? (((w2_up >= 0) ? D.qYb[b3f_yb_off(&C, u4, w2_up, w1, c_t2)] : 0) - fc2) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, ck, kc, 1, 1, D2, 1, D2, kc, 0, 0, 0,
                                         fc0, fc1, fc2, v, (ck == 2) ? qc_b : fc1b, (ck == 3) ? qc_b : fc2b, qc_b, 1, dc0, dc1, dc2, q2_up);
                        }
                        if (zf_fail(&W.acc, mode)) fails_yb++;
                    }
                }
            }
        }
        gw_clear(&W);
    }
    printf("  [Yb] checked %zu states: fails = %lld\n", C.yb_len, fails_yb); fflush(stdout);

    /* 5d. Verify S3 */
    #pragma omp parallel num_threads(nthreads) reduction(+:fails_s3)
    {
        GmpWork W; gw_init(&W);
        #pragma omp for schedule(dynamic)
        for (int u4 = 0; u4 < C.NS3; u4++) {
            int Lu = C.s3_L[u4], ku = C.s3_K[u4], hu = C.s3_H[u4], Du = Lu + ku + hu, u_cu = C.s3_cu[u4], u_g = C.s3_g[u4], u_cu_up = C.cu_of_g[b3f_next_up(&C, u_g)];
            int u4_up = (Du + 1 < C.DS3) ? C.s3_lut[hu][ku][Lu + 1] : -1, u4_bud = (Lu >= 1 && ku + 1 <= C.KS3) ? C.s3_lut[hu][ku + 1][Lu - 1] : -1;
            for (int i1 = 0; i1 < C.nchs; i1++) {
                int w1 = C.chs[i1], k1k = C.kind[w1], L1 = C.L[w1], k1 = C.K[w1];
                int w1_up = b3f_next_up(&C, w1), w1_bud = b3f_next_bud(&C, w1), w1_t1 = (k1k >= 3) ? b3f_t1(&C, k1, 1) : w1, w1_t2 = (k1k >= 4) ? b3f_t2(&C, k1, 1) : w1;
                int w1_t1b = (k1k >= 3 && k1 <= KT) ? b3f_t1(&C, k1 + 1, 1) : w1_t1, w1_t2b = (k1k >= 4 && k1 <= KT) ? b3f_t2(&C, k1 + 1, 1) : w1_t2;
                for (int i2 = i1; i2 < C.nchs; i2++) {
                    int w2 = C.chs[i2], k2k = C.kind[w2], L2 = C.L[w2], k2 = C.K[w2];
                    int w2_up = b3f_next_up(&C, w2), w2_bud = b3f_next_bud(&C, w2), w2_t1 = (k2k >= 3) ? b3f_t1(&C, k2, 1) : w2, w2_t2 = (k2k >= 4) ? b3f_t2(&C, k2, 1) : w2;
                    int w2_t1b = (k2k >= 3 && k2 <= KT) ? b3f_t1(&C, k2 + 1, 1) : w2_t1, w2_t2b = (k2k >= 4 && k2 <= KT) ? b3f_t2(&C, k2 + 1, 1) : w2_t2;
                    int64_t dh12 = ((u_cu_up >= 0) ? D.hC[b3f_hc_off(&C, u_cu_up, w1, w2)] : 0) - D.hC[b3f_hc_off(&C, u_cu, w1, w2)];
                    for (int i3 = i2; i3 < C.nchs; i3++) {
                        int w3 = C.chs[i3], k3k = C.kind[w3], L3 = C.L[w3], k3 = C.K[w3];
                        int w3_up = b3f_next_up(&C, w3), w3_bud = b3f_next_bud(&C, w3), w3_t1 = (k3k >= 3) ? b3f_t1(&C, k3, 1) : w3, w3_t2 = (k3k >= 4) ? b3f_t2(&C, k3, 1) : w3;
                        int w3_t1b = (k3k >= 3 && k3 <= KT) ? b3f_t1(&C, k3 + 1, 1) : w3_t1, w3_t2b = (k3k >= 4 && k3 <= KT) ? b3f_t2(&C, k3 + 1, 1) : w3_t2;
                        int64_t dh13 = ((u_cu_up >= 0) ? D.hC[b3f_hc_off(&C, u_cu_up, w1, w3)] : 0) - D.hC[b3f_hc_off(&C, u_cu, w1, w3)];
                        int64_t dh23 = ((u_cu_up >= 0) ? D.hC[b3f_hc_off(&C, u_cu_up, w2, w3)] : 0) - D.hC[b3f_hc_off(&C, u_cu, w2, w3)];
                        int64_t v = D.qS3[b3f_s3_off(&C, u4, w1, w2, w3)];
                        int64_t qu_up = (u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, w2, w3)] : 0, qu_b = (u4_bud >= 0) ? D.qS3[b3f_s3_off(&C, u4_bud, w1, w2, w3)] : 0;
                        int64_t q1_up = D.qS3[b3f_s3_off(&C, u4, w1_up, w2, w3)], q2_up = D.qS3[b3f_s3_off(&C, u4, w1, w2_up, w3)], q3_up = D.qS3[b3f_s3_off(&C, u4, w1, w2, w3_up)];
                        int64_t q1_b = D.qS3[b3f_s3_off(&C, u4, w1_bud, w2, w3)], q2_b = D.qS3[b3f_s3_off(&C, u4, w1, w2_bud, w3)], q3_b = D.qS3[b3f_s3_off(&C, u4, w1, w2, w3_bud)];
                        size_t bi[4]; int ed[4];
                        l5_couple(x,(L5Tree){4,{-1,0,0,0},{u_g,w1,w2,w3}},iw_nb,bi,ed);
                        qu_b += newborn_int(&D,bi[0]); q1_b += newborn_int(&D,bi[1]);
                        q2_b += newborn_int(&D,bi[2]); q3_b += newborn_int(&D,bi[3]);

                        zf_set_si(&W.acc, 0, 1);
                        zf_add_mul_si(&W.acc, -P, v, (uint64_t)Q, &W);
                        add_centered_source(&W.acc, &C, w3, D.rch[u_g], S_r, dh12, &W);
                        add_centered_source(&W.acc, &C, w2, D.rch[u_g], S_r, dh13, &W);
                        add_centered_source(&W.acc, &C, w1, D.rch[u_g], S_r, dh23, &W);

                        if (k1k == 0) {
                            int64_t d1 = q1_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L1, d1, 1, &W); zf_add_si(&W.acc, d1, (uint64_t)Du, &W);
                            if (L1 >= 1) zf_add_mul_si(&W.acc, (int64_t)L1, q1_b - v, (uint64_t)C.D[w1], &W);
                            if (k1 > 0 && ((mode == 1 && d1 > 0) || (mode == 2 && d1 < 0))) zf_add_mul_si(&W.acc, (int64_t)k1, d1, ed[1], &W);
                        }
                        if (k2k == 0) {
                            int64_t d2 = q2_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L2, d2, 1, &W); zf_add_si(&W.acc, d2, (uint64_t)Du, &W);
                            if (L2 >= 1) zf_add_mul_si(&W.acc, (int64_t)L2, q2_b - v, (uint64_t)C.D[w2], &W);
                            if (k2 > 0 && ((mode == 1 && d2 > 0) || (mode == 2 && d2 < 0))) zf_add_mul_si(&W.acc, (int64_t)k2, d2, ed[2], &W);
                        }
                        if (k3k == 0) {
                            int64_t d3 = q3_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)L3, d3, 1, &W); zf_add_si(&W.acc, d3, (uint64_t)Du, &W);
                            if (L3 >= 1) zf_add_mul_si(&W.acc, (int64_t)L3, q3_b - v, (uint64_t)C.D[w3], &W);
                            if (k3 > 0 && ((mode == 1 && d3 > 0) || (mode == 2 && d3 < 0))) zf_add_mul_si(&W.acc, (int64_t)k3, d3, ed[3], &W);
                        }
                        {
                            int64_t du = qu_up - v;
                            zf_add_mul_si(&W.acc, (int64_t)Lu, du, 1, &W);
                            if (k1k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w1], &W);
                            if (k2k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w2], &W);
                            if (k3k == 0) zf_add_si(&W.acc, du, (uint64_t)C.D[w3], &W);
                            if (Lu >= 1) zf_add_mul_si(&W.acc, (int64_t)Lu, qu_b - v, (uint64_t)Du, &W);
                            if (hu + ku - 3 > 0 && ((mode == 1 && du > 0) || (mode == 2 && du < 0))) zf_add_mul_si(&W.acc, (int64_t)(hu + ku - 3), du, ed[0], &W);
                        }
                        if (k1k >= 2) {
                            int64_t f10 = D.qS3[b3f_s3_off(&C, u4, C.INF_ID, w2, w3)], f11 = (k1k >= 3) ? D.qS3[b3f_s3_off(&C, u4, w1_t1, w2, w3)] : v, f12 = (k1k >= 4) ? D.qS3[b3f_s3_off(&C, u4, w1_t2, w2, w3)] : v;
                            int64_t f11b = (k1k >= 3) ? D.qS3[b3f_s3_off(&C, u4, w1_t1b, w2, w3)] : q1_b, f12b = (k1k >= 4) ? D.qS3[b3f_s3_off(&C, u4, w1_t2b, w2, w3)] : q1_b;
                            int64_t d10 = ((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, C.INF_ID, w2, w3)] : 0) - f10;
                            int64_t d11 = (k1k >= 3) ? (((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1_t1, w2, w3)] : 0) - f11) : 0;
                            int64_t d12 = (k1k >= 4) ? (((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1_t2, w2, w3)] : 0) - f12) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k1k, k1, 1, 1, Du, 1, Du, k1, 0, 0, 0,
                                         f10, f11, f12, v, (k1k == 2) ? q1_b : f11b, (k1k == 3) ? q1_b : f12b, q1_b, 1, d10, d11, d12, qu_up);
                        }
                        if (k2k >= 2) {
                            int64_t f20 = D.qS3[b3f_s3_off(&C, u4, w1, C.INF_ID, w3)], f21 = (k2k >= 3) ? D.qS3[b3f_s3_off(&C, u4, w1, w2_t1, w3)] : v, f22 = (k2k >= 4) ? D.qS3[b3f_s3_off(&C, u4, w1, w2_t2, w3)] : v;
                            int64_t f21b = (k2k >= 3) ? D.qS3[b3f_s3_off(&C, u4, w1, w2_t1b, w3)] : q2_b, f22b = (k2k >= 4) ? D.qS3[b3f_s3_off(&C, u4, w1, w2_t2b, w3)] : q2_b;
                            int64_t d20 = ((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, C.INF_ID, w3)] : 0) - f20;
                            int64_t d21 = (k2k >= 3) ? (((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, w2_t1, w3)] : 0) - f21) : 0;
                            int64_t d22 = (k2k >= 4) ? (((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, w2_t2, w3)] : 0) - f22) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k2k, k2, 1, 1, Du, 1, Du, k2, 0, 0, 0,
                                         f20, f21, f22, v, (k2k == 2) ? q2_b : f21b, (k2k == 3) ? q2_b : f22b, q2_b, 1, d20, d21, d22, qu_up);
                        }
                        if (k3k >= 2) {
                            int64_t f30 = D.qS3[b3f_s3_off(&C, u4, w1, w2, C.INF_ID)], f31 = (k3k >= 3) ? D.qS3[b3f_s3_off(&C, u4, w1, w2, w3_t1)] : v, f32 = (k3k >= 4) ? D.qS3[b3f_s3_off(&C, u4, w1, w2, w3_t2)] : v;
                            int64_t f31b = (k3k >= 3) ? D.qS3[b3f_s3_off(&C, u4, w1, w2, w3_t1b)] : q3_b, f32b = (k3k >= 4) ? D.qS3[b3f_s3_off(&C, u4, w1, w2, w3_t2b)] : q3_b;
                            int64_t d30 = ((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, w2, C.INF_ID)] : 0) - f30;
                            int64_t d31 = (k3k >= 3) ? (((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, w2, w3_t1)] : 0) - f31) : 0;
                            int64_t d32 = (k3k >= 4) ? (((u4_up >= 0) ? D.qS3[b3f_s3_off(&C, u4_up, w1, w2, w3_t2)] : 0) - f32) : 0;
                            add_tail_gmp(&W.acc, &W, mode, D0, KT, S_r, k3k, k3, 1, 1, Du, 1, Du, k3, 0, 0, 0,
                                         f30, f31, f32, v, (k3k == 2) ? q3_b : f31b, (k3k == 3) ? q3_b : f32b, q3_b, 1, d30, d31, d32, qu_up);
                        }
                        if (zf_fail(&W.acc, mode)) fails_s3++;
                    }
                }
            }
        }
        gw_clear(&W);
    }
    printf("  [S3] checked %zu states: fails = %lld\n", C.s3_len, fails_s3); fflush(stdout);

    long long fails_l5 = verify_l5(&D,nthreads);
    printf("  [L5] checked %zu states: fails = %lld\n",x->rows,fails_l5); fflush(stdout);
    /* 6. Exact comparison bounds */
    mpq_t lower, upper; mpq_inits(lower,upper,NULL);
    int comp_ok = verify_comparison_bounds(&D,lower,upper);

    /* 7. Built-in exact-rational GMP concrete tree drift verification (boundary + random trees!) */
    int tree_fails = 0;
    if (ntrees_check > 0) {
        /* Hand-crafted boundary trees hitting T_3, star at D0-1, D0, D0+1, broom/caterpillar hitting KT+1, deep path */
        {
            int par_T3[3] = { -1, 0, 0 };
            if (!check_concrete_tree_q(&D,3,par_T3,lower,upper)) tree_fails++;
            /* Core path, apex at D0-1, last endpoint at D0; births
             * exercise the split crossing and exact core promoters. */
            int ns = 2*D0+5, *ps = b3f_xmalloc(ns*sizeof *ps);
            ps[0] = -1;
            for (int i = 1; i < 6; ++i) ps[i] = i-1;
            for (int i = 6; i < D0+4; ++i) ps[i] = 0;
            for (int i = D0+4; i < ns; ++i) ps[i] = i >= ns-2 ? 2 : 5;
            if (!check_concrete_tree_q(&D,ns,ps,lower,upper)) tree_fails++;
            /* High core-degree hub with k = KT + 2 core children! */
            int k_big = KT + 2, nh = 1 + 2 * k_big;
            int *ph = (int*)b3f_xmalloc(nh * sizeof(int)); ph[0] = -1;
            for (int i = 0; i < k_big; i++) { ph[1 + 2 * i] = 0; ph[2 + 2 * i] = 1 + 2 * i; }
            if (!check_concrete_tree_q(&D,nh,ph,lower,upper)) tree_fails++;
            /* Deep path hitting P4, Ytop, Ybot, S3 + tail middle hV/hC */
            int np = D0 + 14; int *pp = (int*)b3f_xmalloc(np * sizeof(int)); pp[0] = -1; pp[1] = 0; pp[2] = 0; pp[3] = 1; pp[4] = 1; pp[5] = 3; pp[6] = 3; pp[7] = 5;
            for (int i = 8; i < np; i++) pp[i] = (i % 3 == 0) ? 1 : ((i % 3 == 1) ? 3 : 0);
            if (!check_concrete_tree_q(&D,np,pp,lower,upper)) tree_fails++;
            free(ps); free(ph); free(pp);
        }
        /* Deterministic PRNG random friend trees up to n = 60 */
        uint64_t rng = 0x123456789ABCDEF0ULL;
        for (int ti = 0; ti < ntrees_check; ti++) {
            int n = 4 + (ti * 3) % 55;
            int *pr = (int*)b3f_xmalloc(n * sizeof(int));
            int *dg = (int*)b3f_xcalloc(n, sizeof(int));
            pr[0] = -1; pr[1] = 0; pr[2] = 0; dg[0] = 2; dg[1] = 1; dg[2] = 1;
            for (int j = 3; j < n; j++) {
                rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                int V = (int)(rng % (uint64_t)j);
                rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                int r = (int)(rng % (uint64_t)dg[V]), Wv = -1;
                if (pr[V] >= 0) { if (r == 0) Wv = pr[V]; else r--; }
                if (Wv < 0) { for (int u = 1; u < j; u++) if (pr[u] == V) { if (r == 0) { Wv = u; break; } r--; } }
                pr[j] = Wv; dg[j] = 1; dg[Wv]++;
            }
            if (!check_concrete_tree_q(&D,n,pr,lower,upper)) tree_fails++;
            free(pr); free(dg);
        }
        printf("  [TREES] exact-rational GMP drift check on %d concrete/boundary trees: fails = %d\n", ntrees_check + 4, tree_fails);
        fflush(stdout);
    }

    long long total_fails = fails_f + fails_g + fails_hv + fails_hc + fails_p4 + fails_yt + fails_yb + fails_s3 + fails_l5 + (!comp_ok) + tree_fails;
    printf("=== VERDICT: %s (total_fails = %lld, exact P/Q = %lld/%lld) ===\n",
           (total_fails == 0) ? "PASS" : "FAIL", total_fails, (long long)P, (long long)Q);

    free(dLf);
    free((void*)D.rpar); free((void*)D.rch); free((void*)D.trpar); free((void*)D.trch);
    free((void*)D.c2_f); free((void*)D.c3_f); free((void*)D.Yp_ref); free((void*)D.Yw_ref);
    free((void*)D.f); free((void*)D.g); free((void*)D.hV); free((void*)D.hC);
    free((void*)D.qP4); free((void*)D.qYt); free((void*)D.qYb); free((void*)D.qS3);
    free((void*)D.q5); l5_free(x); mpq_clears(lower,upper,NULL);
    b3f_layout_free(&C);
    return (total_fails == 0) ? 0 : 1;
}
