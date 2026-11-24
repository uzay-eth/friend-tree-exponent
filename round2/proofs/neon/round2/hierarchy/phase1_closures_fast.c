/*
 * Fast, exact-normalization Phase 1 closures:
 * 1. IR-0 (rho=0)
 * 2. MF-sym (rho = alpha * c, exact E[k]=1 normalization!)
 * 3. MF-dir (rho = alpha_par + k * alpha_ch, exact E[k]=1 normalization!)
 * 4. Cavity-2body (directed core edge p -> w) with fast 1-body bisection (Dmax1=400 + 3-term EM tail)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    double M;
    double mass;
    double E_invD_ch;
    double E_invD_par; /* uses exact sum_k = mass! */
    double E_invD_sym; /* uses exact sum_(k+1) = 2*mass! */
    double last_bD;
} OneBodyRes;

static OneBodyRes solve_1body(double mu, int Dmax, double ap, double ac, int mode_sym) {
    double *prev = calloc(Dmax + 2, sizeof(double));
    double *curr = calloc(Dmax + 2, sizeof(double));
    double sum_b = 0.0, sum_q = 0.0;
    double sum_invD = 0.0, sum_k_invD = 0.0, sum_c_invD = 0.0;
    double last_bD = 0.0, last_invD = 0.0, last_kinvD = 0.0;

    for (int D = 2; D <= Dmax; D++) {
        double bD = 0.0, qD = 0.0, kiD = 0.0;
        for (int k = 0; k <= D - 1; k++) {
            int L = D - 1 - k;
            double rho = mode_sym ? (ap * (k + 1)) : (ap + ac * k);
            double rhs = (D == 2 && k == 0) ? 1.0 : 0.0;
            if (L >= 1) {
                double rho0 = mode_sym ? (ap * (k + 1)) : (ap + ac * k);
                rhs += ((L - 1) + rho0) * prev[k];
            }
            if (k >= 1) rhs += ((double)(L + 1) / D) * curr[k - 1];
            double q = rhs / (mu + L + rho + (double)L / D);
            curr[k] = q;
            bD += ((double)L / D) * q;
            qD += q;
            kiD += (k * q) / D;
        }
        sum_b += bD; sum_q += qD;
        sum_invD += qD / D;
        sum_k_invD += kiD;
        sum_c_invD += kiD + qD / D;
        last_bD = bD; last_invD = qD / D; last_kinvD = kiD;
        double *t = prev; prev = curr; curr = t;
    }
    free(prev); free(curr);

    double tail_b = last_bD * (Dmax / mu - 0.5 + (1.0 + mu) / (12.0 * Dmax));
    double tail_q = tail_b;
    double tail_invD = last_invD * (Dmax / (1.0 + mu) - 0.5);
    /* k*q/D ~ A * log(D) * D^{-(2+mu)}, tail approx last_kinvD * Dmax / (1 + mu) */
    double tail_kinvD = last_kinvD * (Dmax / (1.0 + mu) - 0.5);

    double tot_q = sum_q + tail_q;
    OneBodyRes r;
    r.M = sum_b + tail_b;
    r.mass = tot_q;
    r.E_invD_ch  = (sum_invD + tail_invD) / tot_q;
    /* Exact tree identity: sum k q = sum q = tot_q, so sum (k+1) q = 2 * tot_q! */
    r.E_invD_par = (sum_k_invD + tail_kinvD) / tot_q;
    r.E_invD_sym = (sum_c_invD + tail_invD + tail_kinvD) / (2.0 * tot_q);
    r.last_bD = last_bD;
    return r;
}

static double find_mu_1body(int Dmax, double ap, double ac, int mode_sym, OneBodyRes *out) {
    double lo = 0.2, hi = 0.9;
    for (int it = 0; it < 70; it++) {
        double mid = 0.5 * (lo + hi);
        OneBodyRes r = solve_1body(mid, Dmax, ap, ac, mode_sym);
        if (r.M > 1.0) lo = mid; else hi = mid;
    }
    double mu = 0.5 * (lo + hi);
    if (out) *out = solve_1body(mu, Dmax, ap, ac, mode_sym);
    return mu;
}

/* Fast 2-body cavity closure */
static void solve_cavity_2body(int Dcut, int Dmax1) {
    int *off = calloc(Dcut + 2, sizeof(int));
    int Nst = 0;
    for (int D = 2; D <= Dcut; D++) { off[D] = Nst; Nst += D; }

    double *rpar = calloc(Nst, sizeof(double));
    double *rch  = calloc(Nst, sizeof(double));
    for (int i = 0; i < Nst; i++) { rpar[i] = 0.14; rch[i] = 0.29; }

    double *E = calloc((size_t)Nst * Nst, sizeof(double));
    double *q1 = calloc(Nst, sizeof(double));
    double *q_prev = calloc(Dmax1 + 2, sizeof(double));
    double *q_curr = calloc(Dmax1 + 2, sizeof(double));
    double *E_inf_w = calloc(Nst, sizeof(double));
    double *E_p_inf = calloc(Nst, sizeof(double));

    double mu = 0.5657, rpar_tail = 0.04, rch_tail = 0.26;

    for (int outer = 0; outer < 30; outer++) {
        double lo = 0.54, hi = 0.59;
        for (int bi = 0; bi < 45; bi++) {
            double mid = 0.5 * (lo + hi);
            double sum_b = 0.0, last_bD = 0.0;
            for (int D = 2; D <= Dmax1; D++) {
                double bD = 0.0;
                for (int k = 0; k <= D - 1; k++) {
                    int L = D - 1 - k;
                    double rp = (D <= Dcut) ? rpar[off[D] + k] : rpar_tail;
                    double rc = (D <= Dcut) ? rch[off[D] + k]  : rch_tail;
                    double rho = rp + k * rc;
                    double rhs = (D == 2 && k == 0) ? 1.0 : 0.0;
                    if (L >= 1) {
                        double rp0 = (D - 1 <= Dcut) ? rpar[off[D - 1] + k] : rpar_tail;
                        double rc0 = (D - 1 <= Dcut) ? rch[off[D - 1] + k]  : rch_tail;
                        rhs += ((L - 1) + rp0 + k * rc0) * q_prev[k];
                    }
                    if (k >= 1) rhs += ((double)(L + 1) / D) * q_curr[k - 1];
                    q_curr[k] = rhs / (mid + L + rho + (double)L / D);
                    bD += ((double)L / D) * q_curr[k];
                }
                sum_b += bD; last_bD = bD;
                double *t = q_prev; q_prev = q_curr; q_curr = t;
            }
            double M = sum_b + last_bD * (Dmax1 / mid - 0.5 + (1.0 + mid) / (12.0 * Dmax1));
            if (M > 1.0) lo = mid; else hi = mid;
        }
        mu = 0.5 * (lo + hi);

        /* Compute q1 and tail 1/D averages */
        double tail_p_num = 0.0, tail_p_den = 0.0; /* k-weighted for parents D > Dcut */
        double tail_w_num = 0.0, tail_w_den = 0.0; /* unweighted for children D > Dcut */
        double sum_q_cut = 0.0, sum_b_cut = 0.0, last_bD = 0.0;
        for (int D = 2; D <= Dmax1; D++) {
            double bD = 0.0;
            for (int k = 0; k <= D - 1; k++) {
                int L = D - 1 - k;
                double rp = (D <= Dcut) ? rpar[off[D] + k] : rpar_tail;
                double rc = (D <= Dcut) ? rch[off[D] + k]  : rch_tail;
                double rho = rp + k * rc;
                double rhs = (D == 2 && k == 0) ? 1.0 : 0.0;
                if (L >= 1) {
                    double rp0 = (D - 1 <= Dcut) ? rpar[off[D - 1] + k] : rpar_tail;
                    double rc0 = (D - 1 <= Dcut) ? rch[off[D - 1] + k]  : rch_tail;
                    rhs += ((L - 1) + rp0 + k * rc0) * q_prev[k];
                }
                if (k >= 1) rhs += ((double)(L + 1) / D) * q_curr[k - 1];
                q_curr[k] = rhs / (mu + L + rho + (double)L / D);
                bD += ((double)L / D) * q_curr[k];
                if (D <= Dcut) { q1[off[D] + k] = q_curr[k]; sum_q_cut += q_curr[k]; }
                else {
                    tail_w_num += q_curr[k] / D; tail_w_den += q_curr[k];
                    tail_p_num += k * q_curr[k] / D; tail_p_den += k * q_curr[k];
                }
            }
            if (D <= Dcut) sum_b_cut += bD;
            last_bD = bD;
            double *t = q_prev; q_prev = q_curr; q_curr = t;
        }
        /* Include D > Dmax1 in tail_w_den (tot_q = 1/mu) and tail_p_den (exact sum k q = 1/mu)! */
        double tot_q = 1.0 / mu;
        double sum_k_cut = 0.0;
        for (int D = 2; D <= Dcut; D++)
            for (int k = 1; k <= D - 1; k++)
                sum_k_cut += k * q1[off[D] + k];
        double tail_w_den_exact = tot_q - sum_q_cut;
        double tail_p_den_exact = tot_q - sum_k_cut;
        double invD_w_tail = tail_w_num / tail_w_den_exact;
        double invD_p_tail = tail_p_num / tail_p_den_exact;
        double birth_from_tail = 1.0 - sum_b_cut;
        if (birth_from_tail < 0) birth_from_tail = 0;

        /* 2-body forward sweep */
        memset(E, 0, (size_t)Nst * Nst * sizeof(double));
        for (int Dp = 2; Dp <= Dcut; Dp++) {
            double invDp = 1.0 / Dp;
            for (int kp = 1; kp <= Dp - 1; kp++) {
                int Lp = Dp - 1 - kp;
                int sp = off[Dp] + kp;
                double base_p = Lp + rpar[sp] + (kp - 1) * rch[sp];
                double bp = (double)Lp * invDp;
                double birth_p = ((double)(Lp + 1) * invDp) * q1[off[Dp] + (kp - 1)];
                int sp_Lm = (Lp >= 1 && Dp >= 3) ? (off[Dp - 1] + kp) : -1;
                double base_p0 = (sp_Lm >= 0) ? ((Lp - 1) + rpar[sp_Lm] + (kp - 1) * rch[sp_Lm]) : 0.0;
                int sp_km = (kp >= 2) ? (off[Dp] + (kp - 1)) : -1;
                double b_pkm = (sp_km >= 0) ? ((double)(Lp + 1) * invDp) : 0.0;

                double *Ep = E + (size_t)sp * Nst;
                const double *EpLm = (sp_Lm >= 0) ? (E + (size_t)sp_Lm * Nst) : NULL;
                const double *Epkm = (sp_km >= 0) ? (E + (size_t)sp_km * Nst) : NULL;

                for (int Dw = 2; Dw <= Dcut; Dw++) {
                    double invDw = 1.0 / Dw;
                    for (int kw = 0; kw <= Dw - 1; kw++) {
                        int Lw = Dw - 1 - kw;
                        int sw = off[Dw] + kw;
                        double rhs = (Dw == 2 && kw == 0) ? birth_p : 0.0;
                        if (EpLm) rhs += (base_p0 + invDw) * EpLm[sw];
                        if (Epkm) rhs += b_pkm * Epkm[sw];
                        if (Lw >= 1 && Dw >= 3) {
                            int sw0 = off[Dw - 1] + kw;
                            rhs += ((Lw - 1) + invDp + kw * rch[sw0]) * Ep[sw0];
                        }
                        if (kw >= 1) {
                            int sw0 = off[Dw] + (kw - 1);
                            rhs += ((double)(Lw + 1) * invDw) * Ep[sw0];
                        }
                        double denom = mu + (base_p + invDw) + bp + (Lw + invDp + kw * rch[sw]) + (double)Lw * invDw;
                        Ep[sw] = rhs / denom;
                    }
                }
            }
        }

        /* E_inf_w and E_p_inf */
        memset(E_inf_w, 0, Nst * sizeof(double));
        for (int Dw = 2; Dw <= Dcut; Dw++) {
            double invDw = 1.0 / Dw;
            for (int kw = 0; kw <= Dw - 1; kw++) {
                int Lw = Dw - 1 - kw;
                int sw = off[Dw] + kw;
                double rhs = (Dw == 2 && kw == 0) ? birth_from_tail : 0.0;
                for (int kp = 1; kp <= Dcut - 1; kp++) {
                    int Lp = Dcut - 1 - kp;
                    int sp = off[Dcut] + kp;
                    rhs += (Lp + invDw + rpar[sp] + (kp - 1) * rch[sp]) * E[(size_t)sp * Nst + sw];
                }
                if (Lw >= 1 && Dw >= 3) {
                    int sw0 = off[Dw - 1] + kw;
                    rhs += ((Lw - 1) + invD_p_tail + kw * rch[sw0]) * E_inf_w[sw0];
                }
                if (kw >= 1) rhs += ((double)(Lw + 1) * invDw) * E_inf_w[off[Dw] + (kw - 1)];
                E_inf_w[sw] = rhs / (mu + Lw + invD_p_tail + kw * rch[sw] + (double)Lw * invDw);
            }
        }
        memset(E_p_inf, 0, Nst * sizeof(double));
        for (int Dp = 2; Dp <= Dcut; Dp++) {
            double invDp = 1.0 / Dp;
            for (int kp = 1; kp <= Dp - 1; kp++) {
                int Lp = Dp - 1 - kp;
                int sp = off[Dp] + kp;
                double rhs = 0.0;
                for (int kw = 0; kw <= Dcut - 1; kw++) {
                    int Lw = Dcut - 1 - kw;
                    int sw = off[Dcut] + kw;
                    rhs += (Lw + invDp + kw * rch[sw]) * E[(size_t)sp * Nst + sw];
                }
                if (Lp >= 1 && Dp >= 3) {
                    int sp0 = off[Dp - 1] + kp;
                    rhs += ((Lp - 1) + invD_w_tail + rpar[sp0] + (kp - 1) * rch[sp0]) * E_p_inf[sp0];
                }
                if (kp >= 2) rhs += ((double)(Lp + 1) * invDp) * E_p_inf[off[Dp] + (kp - 1)];
                E_p_inf[sp] = rhs / (mu + Lp + invD_w_tail + rpar[sp] + (kp - 1) * rch[sp] + (double)Lp * invDp);
            }
        }

        double maxdiff = 0.0;
        for (int D = 2; D <= Dcut; D++) {
            for (int k = 0; k <= D - 1; k++) {
                int s = off[D] + k;
                double num_p = invD_p_tail * E_inf_w[s], den_p = E_inf_w[s];
                for (int Dp = 2; Dp <= Dcut; Dp++) {
                    double invDp = 1.0 / Dp;
                    for (int kp = 1; kp <= Dp - 1; kp++) {
                        double e = E[(size_t)(off[Dp] + kp) * Nst + s];
                        num_p += invDp * e; den_p += e;
                    }
                }
                double nrp = num_p / den_p;
                if (fabs(nrp - rpar[s]) > maxdiff) maxdiff = fabs(nrp - rpar[s]);
                rpar[s] = 0.7 * nrp + 0.3 * rpar[s];

                if (k >= 1) {
                    double num_c = invD_w_tail * E_p_inf[s], den_c = E_p_inf[s];
                    for (int Dw = 2; Dw <= Dcut; Dw++) {
                        double invDw = 1.0 / Dw;
                        for (int kw = 0; kw <= Dw - 1; kw++) {
                            double e = E[(size_t)s * Nst + (off[Dw] + kw)];
                            num_c += invDw * e; den_c += e;
                        }
                    }
                    double nrc = num_c / den_c;
                    if (fabs(nrc - rch[s]) > maxdiff) maxdiff = fabs(nrc - rch[s]);
                    rch[s] = 0.7 * nrc + 0.3 * rch[s];
                } else rch[s] = 0.0;
            }
        }
        double tp_num = 0, tp_den = 0, tc_num = 0, tc_den = 0;
        for (int k = 0; k <= Dcut - 1; k++) {
            int s = off[Dcut] + k;
            tp_num += rpar[s] * q1[s]; tp_den += q1[s];
            if (k >= 1) { tc_num += rch[s] * k * q1[s]; tc_den += k * q1[s]; }
        }
        rpar_tail = tp_num / tp_den;
        rch_tail  = tc_num / tc_den;

        if (outer == 29 || maxdiff < 1e-11) {
            printf("Cavity-2body (Dcut=%2d): it=%2d  mu = %.9f  maxdiff=%.1e  rpar(1,0)=%.6f rpar(0,1)=%.6f rch(0,1)=%.6f q(1,0)=%.6f q(0,1)=%.6f\n",
                   Dcut, outer, mu, maxdiff,
                   rpar[off[2]+0], rpar[off[2]+1], rch[off[2]+1],
                   q1[off[2]+0] * mu, q1[off[2]+1] * mu);
            break;
        }
    }
    free(off); free(rpar); free(rch); free(E); free(q1); free(q_prev); free(q_curr); free(E_inf_w); free(E_p_inf);
}

int main(void) {
    OneBodyRes r;
    printf("=== 1. MF-sym (exact E[k]=1 normalization) ===\n");
    for (int Dmax = 1000; Dmax <= 4000; Dmax *= 2) {
        double alpha = 0.2, mu = 0.56;
        for (int it = 0; it < 30; it++) {
            mu = find_mu_1body(Dmax, alpha, 0.0, 1, &r);
            alpha = r.E_invD_sym;
        }
        printf("Dmax=%4d: mu_MF_sym = %.10f  alpha = %.10f\n", Dmax, mu, alpha);
    }
    printf("\n=== 2. MF-dir (exact E[k]=1 normalization) ===\n");
    for (int Dmax = 1000; Dmax <= 4000; Dmax *= 2) {
        double ap = 0.14, ac = 0.29, mu = 0.56;
        for (int it = 0; it < 30; it++) {
            mu = find_mu_1body(Dmax, ap, ac, 0, &r);
            ap = r.E_invD_par; ac = r.E_invD_ch;
        }
        printf("Dmax=%4d: mu_MF_dir = %.10f  alpha_par = %.10f  alpha_ch = %.10f\n", Dmax, mu, ap, ac);
    }
    printf("\n=== 3. Cavity-2body (directed core edge) ===\n");
    solve_cavity_2body(20, 400);
    solve_cavity_2body(35, 400);
    solve_cavity_2body(50, 400);
    solve_cavity_2body(65, 500);
    return 0;
}
