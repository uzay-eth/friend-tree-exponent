#include "marginal_format.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <string>

using namespace std;

/*
 * Untrusted proposal builder & exact-integer envelope exporter for FT2_M5V1.
 * Usage: build_marginal <tiny|med|high> <obs_id 1..8> <lower|upper>
 *                       <mu_lo_P> <mu_lo_Q> <mu_hi_P> <mu_hi_Q> <out.bin>
 *
 * Solves the signed unoriented 5-vertex coboundary DP at lam = 0.5*(mu_lo + mu_hi)
 * with gauge f_nb = 0 and drives f(J) -> 0 via secant updates on theta, then exports
 * quantized int64_t tables + exact 128-bit cross-level envelopes.
 */

static void proj_1d_fp(const Q5ProjMap& PM, const double* val_fine, double* mid, double* rad) {
    double v_inf = val_fine[PM.F_INF];
    mid[PM.C_INF] = v_inf; rad[PM.C_INF] = 0.0;
    for (int e = 0; e < PM.NC; ++e)
        if (PM.exact_fine[e] >= 0) { mid[e] = val_fine[PM.exact_fine[e]]; rad[e] = 0.0; }
    for (int e : PM.non_exact) {
        double zmin = 1e30, zmax = -1e30;
        for (auto fe : PM.buckets[e]) {
            double t = (double)fe.num / fe.den;
            double z = (val_fine[fe.fw] - (1.0 - t) * v_inf) / t;
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        mid[e] = 0.5 * (zmin + zmax);
        rad[e] = 0.5 * (zmax - zmin);
    }
}

static void proj_2d_gen_fp(const Q5ProjMap& PM, const vector<double>& Mf, vector<double>& mid, vector<double>& rad) {
    int NC = PM.NC, NF = PM.NF, CINF = PM.C_INF, FINF = PM.F_INF;
    mid.assign((size_t)NC * NC, 0.0); rad.assign((size_t)NC * NC, 0.0);
    double v_ii = Mf[(size_t)FINF * NF + FINF];
    mid[(size_t)CINF * NC + CINF] = v_ii;
    vector<double> r_inf(NF), c_inf(NF), mr(NC), rr(NC), mc(NC), rc(NC);
    for (int w = 0; w < NF; ++w) { r_inf[w] = Mf[(size_t)w * NF + FINF]; c_inf[w] = Mf[(size_t)FINF * NF + w]; }
    proj_1d_fp(PM, r_inf.data(), mr.data(), rr.data());
    proj_1d_fp(PM, c_inf.data(), mc.data(), rc.data());
    for (int e = 0; e < NC; ++e) {
        mid[(size_t)e * NC + CINF] = mr[e]; rad[(size_t)e * NC + CINF] = rr[e];
        mid[(size_t)CINF * NC + e] = mc[e]; rad[(size_t)CINF * NC + e] = rc[e];
    }
    vector<double> m1(NC), r1(NC), col(NF);
    for (int e1 = 0; e1 < NC; ++e1) {
        if (PM.exact_fine[e1] < 0) continue;
        proj_1d_fp(PM, &Mf[(size_t)PM.exact_fine[e1] * NF], m1.data(), r1.data());
        for (int e2 = 0; e2 < NC; ++e2) { mid[(size_t)e1 * NC + e2] = m1[e2]; rad[(size_t)e1 * NC + e2] = r1[e2]; }
    }
    for (int e2 = 0; e2 < NC; ++e2) {
        if (PM.exact_fine[e2] < 0) continue;
        int fw2 = PM.exact_fine[e2];
        for (int fw1 = 0; fw1 < NF; ++fw1) col[fw1] = Mf[(size_t)fw1 * NF + fw2];
        proj_1d_fp(PM, col.data(), m1.data(), r1.data());
        for (int e1 : PM.non_exact) { mid[(size_t)e1 * NC + e2] = m1[e1]; rad[(size_t)e1 * NC + e2] = r1[e1]; }
    }
    for (int e1 : PM.non_exact) for (int e2 : PM.non_exact) {
        double zmin = 1e30, zmax = -1e30;
        for (auto f1 : PM.buckets[e1]) {
            double t1 = (double)f1.num / f1.den;
            for (auto f2 : PM.buckets[e2]) {
                double t2 = (double)f2.num / f2.den;
                double res = Mf[(size_t)f1.fw * NF + f2.fw] - (1.0 - t1)*(1.0 - t2)*v_ii
                           - t1*(1.0 - t2)*mr[e1] - (1.0 - t1)*t2*mc[e2];
                double z = res / (t1 * t2);
                if (z < zmin) zmin = z;
                if (z > zmax) zmax = z;
            }
        }
        double m = 0.5 * (zmin + zmax), r = 0.0;
        for (auto f1 : PM.buckets[e1]) {
            double t1 = (double)f1.num / f1.den;
            for (auto f2 : PM.buckets[e2]) {
                double t2 = (double)f2.num / f2.den;
                double res = Mf[(size_t)f1.fw * NF + f2.fw] - (1.0 - t1)*(1.0 - t2)*v_ii
                           - t1*(1.0 - t2)*mr[e1] - (1.0 - t1)*t2*mc[e2];
                double slack = t1*(1.0 - t2)*rr[e1] + (1.0 - t1)*t2*rc[e2];
                double need = (fabs(res - t1 * t2 * m) - slack) / (t1 * t2);
                if (need > r) r = need;
            }
        }
        mid[(size_t)e1 * NC + e2] = m; rad[(size_t)e1 * NC + e2] = r;
    }
}

static void proj_3d_sym_fp(const Q5ProjMap& PM, const vector<int>& tF, const vector<int>& tC, int NCt,
                           const vector<double>& Tf, vector<double>& mid, vector<double>& rad) {
    int NF = PM.NF, NC = PM.NC, CINF = PM.C_INF, FINF = PM.F_INF;
    mid.assign(NCt, 0.0); rad.assign(NCt, 0.0);
    auto getF = [&](int a, int b, int c) { return Tf[tF[(a * NF + b) * NF + c]]; };
    auto setC = [&](int a, int b, int c, double m, double r) { int id = tC[(a * NC + b) * NC + c]; mid[id] = m; rad[id] = r; };
    auto gCm = [&](int a, int b, int c) { return mid[tC[(a * NC + b) * NC + c]]; };
    auto gCr = [&](int a, int b, int c) { return rad[tC[(a * NC + b) * NC + c]]; };

    vector<double> M2((size_t)NF * NF), m2, r2;
    for (int u = 0; u < NF; ++u) for (int v = 0; v < NF; ++v) M2[(size_t)u * NF + v] = getF(FINF, u, v);
    proj_2d_gen_fp(PM, M2, m2, r2);
    for (int u = 0; u < NC; ++u) for (int v = u; v < NC; ++v) setC(CINF, u, v, m2[(size_t)u * NC + v], r2[(size_t)u * NC + v]);
    for (int e1 = 0; e1 < NC; ++e1) {
        if (PM.exact_fine[e1] < 0) continue;
        int fw1 = PM.exact_fine[e1];
        for (int u = 0; u < NF; ++u) for (int v = 0; v < NF; ++v) M2[(size_t)u * NF + v] = getF(fw1, u, v);
        proj_2d_gen_fp(PM, M2, m2, r2);
        for (int u = e1; u < NC; ++u) for (int v = u; v < NC; ++v) setC(e1, u, v, m2[(size_t)u * NC + v], r2[(size_t)u * NC + v]);
    }
    int nn = (int)PM.non_exact.size();
    for (int i1 = 0; i1 < nn; ++i1) {
        int e1 = PM.non_exact[i1];
        for (int i2 = i1; i2 < nn; ++i2) {
            int e2 = PM.non_exact[i2];
            for (int i3 = i2; i3 < nn; ++i3) {
                int e3 = PM.non_exact[i3];
                double cm[2][2][2], cr[2][2][2];
                for (int b1=0;b1<2;++b1) for (int b2=0;b2<2;++b2) for (int b3=0;b3<2;++b3) {
                    if (b1 && b2 && b3) { cm[1][1][1] = cr[1][1][1] = 0.0; continue; }
                    cm[b1][b2][b3] = gCm(b1?e1:CINF, b2?e2:CINF, b3?e3:CINF);
                    cr[b1][b2][b3] = gCr(b1?e1:CINF, b2?e2:CINF, b3?e3:CINF);
                }
                double zmin = 1e30, zmax = -1e30;
                for (auto f1 : PM.buckets[e1]) { double t1=(double)f1.num/f1.den, w1[2]={1-t1,t1};
                for (auto f2 : PM.buckets[e2]) { double t2=(double)f2.num/f2.den, w2[2]={1-t2,t2};
                for (auto f3 : PM.buckets[e3]) { double t3=(double)f3.num/f3.den, w3[2]={1-t3,t3};
                    double pred = 0.0;
                    for (int b1=0;b1<2;++b1) for (int b2=0;b2<2;++b2) for (int b3=0;b3<2;++b3)
                        if (!(b1&&b2&&b3)) pred += w1[b1]*w2[b2]*w3[b3]*cm[b1][b2][b3];
                    double z = (getF(f1.fw, f2.fw, f3.fw) - pred) / (t1*t2*t3);
                    if (z < zmin) zmin = z;
                    if (z > zmax) zmax = z;
                }}}
                double m = 0.5*(zmin+zmax), r = 0.0;
                for (auto f1 : PM.buckets[e1]) { double t1=(double)f1.num/f1.den, w1[2]={1-t1,t1};
                for (auto f2 : PM.buckets[e2]) { double t2=(double)f2.num/f2.den, w2[2]={1-t2,t2};
                for (auto f3 : PM.buckets[e3]) { double t3=(double)f3.num/f3.den, w3[2]={1-t3,t3};
                    double pred = 0.0, slack = 0.0;
                    for (int b1=0;b1<2;++b1) for (int b2=0;b2<2;++b2) for (int b3=0;b3<2;++b3)
                        if (!(b1&&b2&&b3)) { double w=w1[b1]*w2[b2]*w3[b3]; pred += w*cm[b1][b2][b3]; slack += w*cr[b1][b2][b3]; }
                    double need = (fabs(getF(f1.fw, f2.fw, f3.fw) - pred - t1*t2*t3*m) - slack) / (t1*t2*t3);
                    if (need > r) r = need;
                }}}
                setC(e1, e2, e3, m, r);
            }
        }
    }
}

static void proj_1d_exact(const Q5ProjMap& PM, const int64_t* vf, int64_t* mid, int64_t* rad) {
    int64_t v_inf = vf[PM.F_INF];
    mid[PM.C_INF] = v_inf; rad[PM.C_INF] = 0;
    for (int e = 0; e < PM.NC; ++e)
        if (PM.exact_fine[e] >= 0) { mid[e] = vf[PM.exact_fine[e]]; rad[e] = 0; }
    for (int e : PM.non_exact) {
        q5_i128 zmin_lo = ((q5_i128)1 << 100), zmax_hi = -((q5_i128)1 << 100);
        for (auto fe : PM.buckets[e]) {
            q5_i128 num = (q5_i128)fe.den * vf[fe.fw] - (q5_i128)(fe.den - fe.num) * v_inf;
            q5_i128 zlo = q5_floor_div(num, fe.num), zhi = q5_ceil_div(num, fe.num);
            if (zlo < zmin_lo) zmin_lo = zlo;
            if (zhi > zmax_hi) zmax_hi = zhi;
        }
        int64_t m = q5_to_i64(q5_floor_div(zmin_lo + zmax_hi, 2), "proj_1d_m");
        q5_i128 r = 0;
        for (auto fe : PM.buckets[e]) {
            q5_i128 diff = q5_abs128((q5_i128)fe.den * vf[fe.fw] - (q5_i128)(fe.den - fe.num) * v_inf - (q5_i128)fe.num * m);
            q5_i128 need = q5_ceil_div(diff, fe.num);
            if (need > r) r = need;
        }
        mid[e] = m; rad[e] = q5_to_i64(r, "proj_1d_r");
    }
}

static void proj_2d_exact(const Q5ProjMap& PM, const vector<int64_t>& Mf, vector<int64_t>& mid, vector<int64_t>& rad) {
    int NC = PM.NC, NF = PM.NF, CINF = PM.C_INF, FINF = PM.F_INF;
    mid.assign((size_t)NC * NC, 0); rad.assign((size_t)NC * NC, 0);
    int64_t v_ii = Mf[(size_t)FINF * NF + FINF];
    mid[(size_t)CINF * NC + CINF] = v_ii;
    vector<int64_t> r_inf(NF), c_inf(NF), mr(NC), rr(NC), mc(NC), rc(NC);
    for (int w = 0; w < NF; ++w) { r_inf[w] = Mf[(size_t)w * NF + FINF]; c_inf[w] = Mf[(size_t)FINF * NF + w]; }
    proj_1d_exact(PM, r_inf.data(), mr.data(), rr.data());
    proj_1d_exact(PM, c_inf.data(), mc.data(), rc.data());
    for (int e = 0; e < NC; ++e) {
        mid[(size_t)e * NC + CINF] = mr[e]; rad[(size_t)e * NC + CINF] = rr[e];
        mid[(size_t)CINF * NC + e] = mc[e]; rad[(size_t)CINF * NC + e] = rc[e];
    }
    vector<int64_t> m1(NC), r1(NC), col(NF);
    for (int e1 = 0; e1 < NC; ++e1) {
        if (PM.exact_fine[e1] < 0) continue;
        proj_1d_exact(PM, &Mf[(size_t)PM.exact_fine[e1] * NF], m1.data(), r1.data());
        for (int e2 = 0; e2 < NC; ++e2) { mid[(size_t)e1 * NC + e2] = m1[e2]; rad[(size_t)e1 * NC + e2] = r1[e2]; }
    }
    for (int e2 = 0; e2 < NC; ++e2) {
        if (PM.exact_fine[e2] < 0) continue;
        int fw2 = PM.exact_fine[e2];
        for (int fw1 = 0; fw1 < NF; ++fw1) col[fw1] = Mf[(size_t)fw1 * NF + fw2];
        proj_1d_exact(PM, col.data(), m1.data(), r1.data());
        for (int e1 : PM.non_exact) { mid[(size_t)e1 * NC + e2] = m1[e1]; rad[(size_t)e1 * NC + e2] = r1[e1]; }
    }
    for (int e1 : PM.non_exact) for (int e2 : PM.non_exact) {
        q5_i128 zmin_lo = ((q5_i128)1 << 100), zmax_hi = -((q5_i128)1 << 100);
        for (auto f1 : PM.buckets[e1]) {
            q5_i128 w10 = f1.den - f1.num, w11 = f1.num;
            for (auto f2 : PM.buckets[e2]) {
                q5_i128 w20 = f2.den - f2.num, w21 = f2.num;
                q5_i128 res = (q5_i128)f1.den * f2.den * Mf[(size_t)f1.fw * NF + f2.fw]
                             - w10 * w20 * v_ii - w11 * w20 * mr[e1] - w10 * w21 * mc[e2];
                q5_i128 d12 = w11 * w21;
                q5_i128 zlo = q5_floor_div(res, d12), zhi = q5_ceil_div(res, d12);
                if (zlo < zmin_lo) zmin_lo = zlo;
                if (zhi > zmax_hi) zmax_hi = zhi;
            }
        }
        int64_t m = q5_to_i64(q5_floor_div(zmin_lo + zmax_hi, 2), "proj_2d_m");
        q5_i128 r = 0;
        for (auto f1 : PM.buckets[e1]) {
            q5_i128 w10 = f1.den - f1.num, w11 = f1.num;
            for (auto f2 : PM.buckets[e2]) {
                q5_i128 w20 = f2.den - f2.num, w21 = f2.num;
                q5_i128 res = (q5_i128)f1.den * f2.den * Mf[(size_t)f1.fw * NF + f2.fw]
                             - w10 * w20 * v_ii - w11 * w20 * mr[e1] - w10 * w21 * mc[e2];
                q5_i128 slack = w11 * w20 * rr[e1] + w10 * w21 * rc[e2];
                q5_i128 d12 = w11 * w21;
                q5_i128 diff = q5_abs128(res - d12 * m);
                if (diff > slack) {
                    q5_i128 need = q5_ceil_div(diff - slack, d12);
                    if (need > r) r = need;
                }
            }
        }
        mid[(size_t)e1 * NC + e2] = m; rad[(size_t)e1 * NC + e2] = q5_to_i64(r, "proj_2d_r");
    }
}

static void proj_3d_sym_exact(const Q5ProjMap& PM, const vector<int>& tF, const vector<int>& tC, int NCt,
                              const vector<int64_t>& Tf, int64_t* mid, int64_t* rad) {
    int NF = PM.NF, NC = PM.NC, CINF = PM.C_INF, FINF = PM.F_INF;
    memset(mid, 0, (size_t)NCt * sizeof(int64_t));
    memset(rad, 0, (size_t)NCt * sizeof(int64_t));
    auto getF = [&](int a, int b, int c) { return Tf[tF[((size_t)a * NF + b) * NF + c]]; };
    auto setC = [&](int a, int b, int c, int64_t m, int64_t r) { int id = tC[((size_t)a * NC + b) * NC + c]; mid[id] = m; rad[id] = r; };
    auto gCm = [&](int a, int b, int c) { return mid[tC[((size_t)a * NC + b) * NC + c]]; };
    auto gCr = [&](int a, int b, int c) { return rad[tC[((size_t)a * NC + b) * NC + c]]; };

    vector<int64_t> M2((size_t)NF * NF), m2, r2;
    for (int u = 0; u < NF; ++u) for (int v = 0; v < NF; ++v) M2[(size_t)u * NF + v] = getF(FINF, u, v);
    proj_2d_exact(PM, M2, m2, r2);
    for (int u = 0; u < NC; ++u) for (int v = u; v < NC; ++v) setC(CINF, u, v, m2[(size_t)u * NC + v], r2[(size_t)u * NC + v]);
    for (int e1 = 0; e1 < NC; ++e1) {
        if (PM.exact_fine[e1] < 0) continue;
        int fw1 = PM.exact_fine[e1];
        for (int u = 0; u < NF; ++u) for (int v = 0; v < NF; ++v) M2[(size_t)u * NF + v] = getF(fw1, u, v);
        proj_2d_exact(PM, M2, m2, r2);
        for (int u = e1; u < NC; ++u) for (int v = u; v < NC; ++v) setC(e1, u, v, m2[(size_t)u * NC + v], r2[(size_t)u * NC + v]);
    }
    int nn = (int)PM.non_exact.size();
    for (int i1 = 0; i1 < nn; ++i1) {
        int e1 = PM.non_exact[i1];
        for (int i2 = i1; i2 < nn; ++i2) {
            int e2 = PM.non_exact[i2];
            for (int i3 = i2; i3 < nn; ++i3) {
                int e3 = PM.non_exact[i3];
                int64_t cm[2][2][2], cr[2][2][2];
                for (int b1=0;b1<2;++b1) for (int b2=0;b2<2;++b2) for (int b3=0;b3<2;++b3) {
                    if (b1 && b2 && b3) { cm[1][1][1] = cr[1][1][1] = 0; continue; }
                    cm[b1][b2][b3] = gCm(b1?e1:CINF, b2?e2:CINF, b3?e3:CINF);
                    cr[b1][b2][b3] = gCr(b1?e1:CINF, b2?e2:CINF, b3?e3:CINF);
                }
                q5_i128 zmin_lo = ((q5_i128)1 << 100), zmax_hi = -((q5_i128)1 << 100);
                for (auto f1 : PM.buckets[e1]) { q5_i128 w1[2]={f1.den-f1.num, f1.num};
                for (auto f2 : PM.buckets[e2]) { q5_i128 w2[2]={f2.den-f2.num, f2.num};
                for (auto f3 : PM.buckets[e3]) { q5_i128 w3[2]={f3.den-f3.num, f3.num};
                    q5_i128 pred = 0;
                    for (int b1=0;b1<2;++b1) for (int b2=0;b2<2;++b2) for (int b3=0;b3<2;++b3)
                        if (!(b1&&b2&&b3)) pred += w1[b1]*w2[b2]*w3[b3]*cm[b1][b2][b3];
                    q5_i128 res = (q5_i128)f1.den*f2.den*f3.den*getF(f1.fw, f2.fw, f3.fw) - pred;
                    q5_i128 d123 = w1[1]*w2[1]*w3[1];
                    q5_i128 zlo = q5_floor_div(res, d123), zhi = q5_ceil_div(res, d123);
                    if (zlo < zmin_lo) zmin_lo = zlo;
                    if (zhi > zmax_hi) zmax_hi = zhi;
                }}}
                int64_t m = q5_to_i64(q5_floor_div(zmin_lo + zmax_hi, 2), "proj_3d_m");
                q5_i128 r = 0;
                for (auto f1 : PM.buckets[e1]) { q5_i128 w1[2]={f1.den-f1.num, f1.num};
                for (auto f2 : PM.buckets[e2]) { q5_i128 w2[2]={f2.den-f2.num, f2.num};
                for (auto f3 : PM.buckets[e3]) { q5_i128 w3[2]={f3.den-f3.num, f3.num};
                    q5_i128 pred = 0, slack = 0;
                    for (int b1=0;b1<2;++b1) for (int b2=0;b2<2;++b2) for (int b3=0;b3<2;++b3)
                        if (!(b1&&b2&&b3)) { q5_i128 w=w1[b1]*w2[b2]*w3[b3]; pred += w*cm[b1][b2][b3]; slack += w*cr[b1][b2][b3]; }
                    q5_i128 res = (q5_i128)f1.den*f2.den*f3.den*getF(f1.fw, f2.fw, f3.fw) - pred;
                    q5_i128 d123 = w1[1]*w2[1]*w3[1];
                    q5_i128 diff = q5_abs128(res - d123 * m);
                    if (diff > slack) {
                        q5_i128 need = q5_ceil_div(diff - slack, d123);
                        if (need > r) r = need;
                    }
                }}}
                setC(e1, e2, e3, m, q5_to_i64(r, "proj_3d_r"));
            }
        }
    }
}

struct MarginalSolver {
    const Q5Layout* C = nullptr;
    int nk = 0;
    int64_t obs_id = 1, obs_sign = 1;
    double KNOTS[Q5_MAX_KNOTS] = {0};
    double lam_lo = 0.5657, lam_hi = 0.56575, lam = 0.565725;
    double A_tail = 0.0, B_tail = 0.0;
    vector<double> r1, knots, f, g, h3, p4, s3, p5, b5, s4;

    void init(const Q5Layout* c_ptr, int64_t oid, int64_t osgn, double llo, double lhi) {
        C = c_ptr; nk = (int)C->h.nknots;
        obs_id = oid; obs_sign = osgn;
        lam_lo = llo; lam_hi = lhi; lam = 0.5 * (llo + lhi);
        for (int i = 0; i < nk; ++i) KNOTS[i] = (double)C->h.knot_num[i] / (double)C->h.Q_K;
        r1.assign(C->n_r1, 0.28);
        knots.assign(C->n_knots, 0.0);
        f.assign(C->n_f, 0.0);
        g.assign(C->n_g, 0.0);
        h3.assign(C->n_h3, 0.0);
        p4.assign(C->n_p4, 0.0);
        s3.assign(C->n_s3, 0.0);
        p5.assign(C->n_p5, 0.0);
        b5.assign(C->n_b5, 0.0);
        s4.assign(C->n_s4, 0.0);
    }

    static inline double invD_end(const Q5EndSpace& E, int u) { return (E.kind[u] == 2) ? 0.0 : (1.0 / E.D[u]); }
    static inline double b_lo_end(const Q5EndSpace& E, int u) {
        if (E.kind[u] == 0) return (E.L[u] >= 1) ? (double)E.L[u] / E.D[u] : 0.0;
        if (E.kind[u] == 1 && E.C[u] <= E.Kcut) return max(0.0, (double)(E.dcut - E.C[u]) / E.dcut);
        return 0.0;
    }
    static inline double b_hi_end(const Q5EndSpace& E, int u) {
        if (E.kind[u] == 0) return (E.L[u] >= 1) ? (double)E.L[u] / E.D[u] : 0.0;
        if (E.kind[u] == 1 && E.C[u] <= E.Kcut) return 1.0;
        return 0.0;
    }
    static inline double g_lo_end(const Q5EndSpace& E, int u, double inv_nbr) {
        if (E.kind[u] == 1 && E.C[u] <= E.Kcut) return max(0.0, ((double)(E.dcut - E.C[u]) + inv_nbr) / (E.dcut + 1.0));
        return 0.0;
    }
    static inline double r_ovf_end(const Q5EndSpace& E, int u) {
        return (E.kind[u] == 3) ? ((double)E.D[u] - 0.5 * (E.Kcut + 2)) : 0.0;
    }

    inline double eval_pl(const double* a, double x) const {
        if (x <= KNOTS[0]) return a[0];
        if (x >= KNOTS[nk - 1]) return a[nk - 1];
        for (int i = 0; i < nk - 1; ++i)
            if (x <= KNOTS[i + 1]) { double t = (x - KNOTS[i]) / (KNOTS[i + 1] - KNOTS[i]); return a[i] * (1.0 - t) + a[i + 1] * t; }
        return a[nk - 1];
    }

    void project_convex(double* a) const {
        int m = nk - 1, bs[Q5_MAX_KNOTS], be[Q5_MAX_KNOTS], nb = 0;
        double sl[Q5_MAX_KNOTS], wt[Q5_MAX_KNOTS];
        for (int i = 0; i < m; ++i) {
            double dx = KNOTS[i + 1] - KNOTS[i];
            sl[nb] = (a[i + 1] - a[i]) / dx;
            wt[nb] = dx; bs[nb] = i; be[nb] = i; nb++;
            while (nb >= 2 && sl[nb - 2] > sl[nb - 1]) {
                double nw = wt[nb - 2] + wt[nb - 1];
                sl[nb - 2] = (sl[nb - 2] * wt[nb - 2] + sl[nb - 1] * wt[nb - 1]) / nw;
                wt[nb - 2] = nw; be[nb - 2] = be[nb - 1]; nb--;
            }
        }
        double cur = a[0];
        for (int b = 0; b < nb; ++b)
            for (int i = bs[b]; i <= be[b]; ++i) { cur += sl[b] * (KNOTS[i + 1] - KNOTS[i]); a[i + 1] = cur; }
    }

    void opt_knots(double* a, const vector<double>& xs, const vector<double>& wxs,
                   const vector<double>& ys, const vector<double>& wys, int nsteps = 40) const {
        double ym = 0, wyt = 0, wxt = 0;
        for (size_t j = 0; j < ys.size(); ++j) { ym += wys[j] * ys[j]; wyt += wys[j]; }
        for (size_t i = 0; i < xs.size(); ++i) wxt += wxs[i];
        if (wyt <= 0 || wxt <= 0) { for (int k = 0; k < nk; ++k) a[k] = 0; return; }
        ym /= wyt;
        for (int k = 0; k < nk; ++k) a[k] = ym * KNOTS[k];
        double mu[Q5_MAX_KNOTS] = {0};
        for (size_t i = 0; i < xs.size(); ++i) {
            double x = xs[i], w = wxs[i] / wxt;
            if (x <= KNOTS[0]) { mu[0] += w; continue; }
            if (x >= KNOTS[nk - 1]) { mu[nk - 1] += w; continue; }
            for (int k = 0; k < nk - 1; ++k)
                if (x <= KNOTS[k + 1]) { double t = (x - KNOTS[k]) / (KNOTS[k + 1] - KNOTS[k]); mu[k] += w*(1-t); mu[k+1] += w*t; break; }
        }
        double yvar = 0;
        for (size_t j = 0; j < ys.size(); ++j) { double d = ys[j] - ym; yvar += (wys[j] / wyt) * d * d; }
        double ysd = sqrt(yvar + 1e-12), best_a[Q5_MAX_KNOTS], best_obj = 1e30;
        memcpy(best_a, a, nk * sizeof(double));
        for (int step = 0; step < nsteps; ++step) {
            double grad[Q5_MAX_KNOTS], obj = 0;
            for (int k = 0; k < nk; ++k) { grad[k] = mu[k]; obj += mu[k] * a[k]; }
            for (size_t j = 0; j < ys.size(); ++j) {
                double y = ys[j], w = wys[j] / wyt, bv = -1e30; int bk = 0;
                for (int k = 0; k < nk; ++k) {
                    double c = KNOTS[k] * y - a[k];
                    if (c > bv) { bv = c; bk = k; }
                }
                obj += w * bv; grad[bk] -= w;
            }
            if (obj < best_obj) { best_obj = obj; memcpy(best_a, a, nk * sizeof(double)); }
            double lr = 0.25 * ysd / sqrt(1.0 + 0.5 * step);
            for (int k = 0; k < nk; ++k) a[k] -= lr * grad[k];
            double c0 = a[0]; for (int k = 0; k < nk; ++k) a[k] -= c0;
            project_convex(a);
        }
        memcpy(a, best_a, nk * sizeof(double));
    }

    inline double eval_f_ana(int D, int c) const {
        double q = (double)(D - c) / D;
        return A_tail + B_tail * q;
    }

    // Signed linear analytic tail w(q) = A_tail + B_tail * q valid at BOTH lam_lo and lam_hi!
    void update_analytic_f_tail(double theta, double fnb) {
        int DF = (int)C->h.DF, d2 = (int)C->h.d2, J_g = C->G2.idx[2][1];
        double g_inf = g[C->g_off(C->G2.INF, J_g)], g_max = g_inf;
        for (int kap = 1; kap <= (int)C->h.K2 + 1; ++kap) {
            double gv = g[C->g_off(C->G2.t_id(kap), J_g)];
            if (gv > g_max) g_max = gv;
        }
        double Cnb = fnb + g_inf;
        double d_g = d2 * (g_max - g_inf);
        double b = (Cnb >= 0.0) ? (Cnb / (1.0 + lam)) : (Cnb / lam);
        for (int it = 0; it < 6; ++it) {
            double C_eff = Cnb + max(0.0, d_g - b) / (double)DF;
            b = (C_eff >= 0.0) ? (C_eff / (1.0 + lam)) : (C_eff / lam);
        }
        double a = (b >= 0.0) ? ((b - (double)obs_sign * theta) / lam) : (-(double)obs_sign * theta / lam);
        // Cushion so cap_f_ana <= 0 at BOTH endpoints lam_lo and lam_hi!
        double dl = 0.5 * (lam_hi - lam_lo);
        double max_abs_w = max(fabs(a), fabs(a + b)) + 0.1;
        a += (dl * max_abs_w * 2.5 + 1e-8) / lam_lo;
        A_tail = a;
        B_tail = b;
    }

    double eval_theta(double theta, const vector<double>& emp_edge, int max_it, int opt_k_iters, int knot_steps = 40) {
        const Q5Layout& L = *C;
        const auto& h = L.h;
        double lm = lam;
        int nt = (L.n_p4 < 500000) ? 1 : 8;
        int NG = L.G2.N, N3 = L.G3.N, N4 = L.G4.N, N5 = L.G5.N;
        int J_g = L.G2.idx[2][1], J_e3 = L.G3.idx[2][1], J_e4 = L.G4.idx[2][1], J_e5 = L.G5.idx[2][1];

        int maxR4 = L.G4.maxR, maxR5 = L.G5.maxR;
        vector<vector<vector<int>>> tri4_rc(3 * maxR4 + 1);
        for (int rt = 0; rt <= 3 * maxR4; ++rt) tri4_rc[rt].assign(3 * ((int)h.K4e + 2) + 1, {});
        for (int tid = 0; tid < L.N_G4_tri; ++tid) {
            auto [a, b, c] = L.tri_G4[tid];
            int ra = L.G4.is_inner(a)?L.G4.D[a]:((L.G4.kind[a]==1)?(int)h.d4e:maxR4);
            int rb = L.G4.is_inner(b)?L.G4.D[b]:((L.G4.kind[b]==1)?(int)h.d4e:maxR4);
            int rc = L.G4.is_inner(c)?L.G4.D[c]:((L.G4.kind[c]==1)?(int)h.d4e:maxR4);
            tri4_rc[ra+rb+rc][L.G4.C[a]+L.G4.C[b]+L.G4.C[c]].push_back(tid);
        }
        vector<vector<vector<int>>> quad5_rc(4 * maxR5 + 1);
        for (int rq = 0; rq <= 4 * maxR5; ++rq) quad5_rc[rq].assign(4 * ((int)h.K5e + 2) + 1, {});
        for (int qid = 0; qid < L.N_G5_quad; ++qid) {
            auto [a, b, c, d] = L.quad_G5[qid];
            int ra = L.G5.is_inner(a)?L.G5.D[a]:((L.G5.kind[a]==1)?(int)h.d5e:maxR5);
            int rb = L.G5.is_inner(b)?L.G5.D[b]:((L.G5.kind[b]==1)?(int)h.d5e:maxR5);
            int rc = L.G5.is_inner(c)?L.G5.D[c]:((L.G5.kind[c]==1)?(int)h.d5e:maxR5);
            int rd = L.G5.is_inner(d)?L.G5.D[d]:((L.G5.kind[d]==1)?(int)h.d5e:maxR5);
            quad5_rc[ra+rb+rc+rd][L.G5.C[a]+L.G5.C[b]+L.G5.C[c]+L.G5.C[d]].push_back(qid);
        }

        double fJ_res = 0.0;
        for (int it = 0; it < max_it; ++it) {
            update_analytic_f_tail(theta, 0.0);
            vector<double>& fn = f;
            fill(fn.begin(), fn.end(), 0.0);
            for (int D = (int)h.DF - 1; D >= 2; --D) {
                for (int c = D; c >= 0; --c) {
                    int Lv = D - c; double b = (double)Lv / D;
                    double f_up = (D + 1 < h.DF) ? fn[L.f_id(D + 1, c)] : eval_f_ana(D + 1, c);
                    double f_bd = (Lv > 0) ? fn[L.f_id(D, c + 1)] : 0.0;
                    double g_nb = 0.0;
                    if (Lv > 0) {
                        if (D < h.d2) {
                            int id = (c + 1 <= h.K2) ? L.G2.idx[D][c + 1] : L.G2.ovf_idx[D];
                            g_nb = g[L.g_off(id, J_g)];
                        } else {
                            double t = (double)h.d2 / D;
                            g_nb = t * g[L.g_off(L.G2.t_id(c + 1), J_g)] + (1.0 - t) * g[L.g_off(L.G2.INF, J_g)];
                        }
                    }
                    // Gauge f_nb = 0 (at root f[J] == 0, exact match!)
                    double nb = (Lv > 0) ? g_nb : 0.0;
                    double src_f = (double)obs_sign * ((double)m5_eval_obs(obs_id, D, c) - theta);
                    if (D < h.d2 && c >= 1 && c <= h.K2) {
                        double ur = Lv + c * r1[L.G2.idx[D][c]];
                        fn[L.f_id(D, c)] = (ur * f_up + b * (f_bd + nb) + src_f) / (lm + ur + b);
                    } else {
                        double r0 = (double)Lv, r1r = Lv + 0.5 * c;
                        double v0 = (r0 * f_up + b * (f_bd + nb) + src_f) / (lm + r0 + b);
                        double v1 = (r1r * f_up + b * (f_bd + nb) + src_f) / (lm + r1r + b);
                        fn[L.f_id(D, c)] = max(v0, v1);
                    }
                }
            }
            fJ_res = fn[L.f_id(2, 1)];

            vector<double> dLf(L.G2.Nin_ex);
            for (int u = 0; u < L.G2.Nin_ex; ++u) dLf[u] = fn[L.f_id(L.G2.D[u] + 1, L.G2.C[u])] - fn[L.f_id(L.G2.D[u], L.G2.C[u])];

            if (it >= 1 && it <= opt_k_iters) {
                #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1)
                for (int u = 0; u < L.G2.Nin_ex; ++u) {
                    bool in_m3 = (L.G2.D[u] < h.d3m && L.G2.C[u] >= 2 && L.G2.C[u] <= h.K3m);
                    double* ku = &knots[(size_t)u * nk];
                    if (L.G2.C[u] < 2 || in_m3) { fill(ku, ku + nk, 0.0); continue; }
                    vector<double> xs, wxs, ys, wys;
                    int u_up = L.G2.up(u);
                    for (int v = 0; v < NG; ++v) {
                        double w = emp_edge[(size_t)u * NG + v];
                        if (w <= 0) continue;
                        xs.push_back(invD_end(L.G2, v)); wxs.push_back(w);
                        ys.push_back(g[L.g_off(u_up, v)] - g[L.g_off(u, v)]); wys.push_back(w);
                    }
                    opt_knots(ku, xs, wxs, ys, wys, knot_steps);
                }
            }
            vector<double> A_tot(L.n_g, 0.0);
            #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1)
            for (int u = 0; u < NG; ++u) {
                bool fy_u = (L.G2.kind[u] == 0 && L.G2.C[u] >= 2 && !(L.G2.D[u] < h.d3m && L.G2.C[u] <= h.K3m));
                for (int v = u; v < NG; ++v) {
                    bool fy_v = (L.G2.kind[v] == 0 && L.G2.C[v] >= 2 && !(L.G2.D[v] < h.d3m && L.G2.C[v] <= h.K3m));
                    double val = 0.0;
                    if (fy_u) val += (L.G2.C[u] - 1) * eval_pl(&knots[(size_t)u * nk], invD_end(L.G2, v));
                    if (fy_v) val += (L.G2.C[v] - 1) * eval_pl(&knots[(size_t)v * nk], invD_end(L.G2, u));
                    A_tot[L.g_off(u, v)] = val;
                }
            }

            // === STEP 2: g ===
            auto eval_h3_nb = [&](int u_bud, int v) -> double {
                if (L.G2.kind[u_bud] != 0) return 0.0;
                int D = L.G2.D[u_bud], c = L.G2.C[u_bud];
                if (D >= h.d3m || c < 2 || c > h.K3m) return 0.0;
                int m = L.M3.idx[D][c], ce, num, den;
                L.G3.proj_int(L.G2, v, ce, num, den);
                double w = (double)num / den;
                return w * h3[L.h3_off(m, ce, J_e3)] + (1.0 - w) * h3[L.h3_off(m, L.G3.INF, J_e3)];
            };
            vector<double>& gn = g;
            fill(gn.begin(), gn.end(), 0.0);
            int maxR2 = L.G2.maxR;
            for (int rsum = 2 * maxR2; rsum >= 4; --rsum) {
                int ru_lo = max(2, rsum - maxR2), ru_hi = min(maxR2, rsum / 2);
                #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1 && (ru_hi - ru_lo) >= 4)
                for (int ru = ru_lo; ru <= ru_hi; ++ru) {
                    int rv = rsum - ru;
                    if (ru == maxR2 && rv == maxR2) continue;
                    for (int iu = (int)L.G2.by_r[ru].size() - 1; iu >= 0; --iu) {
                        int u = L.G2.by_r[ru][iu];
                        int uk = L.G2.kind[u], Du = L.G2.D[u], cu = L.G2.C[u], Lu = L.G2.L[u];
                        bool u_in = L.G2.is_inner(u);
                        int u_up = L.G2.up(u), u_bud = L.G2.bud(u);
                        double bu_lo = b_lo_end(L.G2, u), bu_hi = b_hi_end(L.G2, u), ru_ovf = r_ovf_end(L.G2, u), inv_u = invD_end(L.G2, u);
                        bool u_m3 = (uk == 0 && Du < h.d3m && cu >= 2 && cu <= h.K3m);
                        bool u_fy = (uk == 0 && cu >= 2 && !u_m3);

                        int iv_lo = (ru == rv) ? iu : 0;
                        for (int iv = (int)L.G2.by_r[rv].size() - 1; iv >= iv_lo; --iv) {
                            int v = L.G2.by_r[rv][iv];
                            int vk = L.G2.kind[v], Dv = L.G2.D[v], cv = L.G2.C[v], Lv = L.G2.L[v];
                            bool v_in = L.G2.is_inner(v);
                            int v_up = L.G2.up(v), v_bud = L.G2.bud(v);
                            double bv_lo = b_lo_end(L.G2, v), bv_hi = b_hi_end(L.G2, v), rv_ovf = r_ovf_end(L.G2, v), inv_v = invD_end(L.G2, v);
                            double gu_lo = g_lo_end(L.G2, u, v_in ? inv_v : 0.0);
                            double gv_lo = g_lo_end(L.G2, v, u_in ? inv_u : 0.0);
                            bool v_m3 = (vk == 0 && Dv < h.d3m && cv >= 2 && cv <= h.K3m);
                            bool v_fy = (vk == 0 && cv >= 2 && !v_m3);

                            double gu_plus = gn[L.g_off(u_up, v)], gv_plus = gn[L.g_off(u, v_up)];
                            double gu_bud_v = gn[L.g_off(u_bud, v)], gv_bud_v = gn[L.g_off(u, v_bud)];
                            double u_nb = (uk == 0 && Lu >= 1) ? eval_h3_nb(u_bud, v) : 0.0;
                            double v_nb = (vk == 0 && Lv >= 1) ? eval_h3_nb(v_bud, u) : 0.0;

                            double num0 = A_tot[L.g_off(u, v)], den0 = lm;
                            if (uk == 0) {
                                num0 += (inv_v - r1[u]) * dLf[u];
                                double ur = Lu + (v_in ? inv_v : 0.0) + (u_m3 ? (cu - 1) * r1[u] : 0.0);
                                num0 += ur * gu_plus; den0 += ur;
                                if (bu_lo > 0) { num0 += bu_lo * (gu_bud_v + u_nb); den0 += bu_lo; }
                            } else if (uk == 3 && v_in) { num0 += inv_v * gu_plus; den0 += inv_v; }
                            if (vk == 0) {
                                num0 += (inv_u - r1[v]) * dLf[v];
                                double vr = Lv + (u_in ? inv_u : 0.0) + (v_m3 ? (cv - 1) * r1[v] : 0.0);
                                num0 += vr * gv_plus; den0 += vr;
                                if (bv_lo > 0) { num0 += bv_lo * (gv_bud_v + v_nb); den0 += bv_lo; }
                            } else if (vk == 3 && u_in) { num0 += inv_u * gv_plus; den0 += inv_u; }
                            double d0_v = (uk == 1 && v_in) ? (gn[L.g_off(L.G2.INF, v_up)] - gn[L.g_off(L.G2.INF, v)]) : 0.0;
                            double d0_u = (vk == 1 && u_in) ? (gn[L.g_off(u_up, L.G2.INF)] - gn[L.g_off(u, L.G2.INF)]) : 0.0;

                            double val = num0 / den0;
                            if (u_fy || v_fy || uk == 1 || vk == 1 || uk == 3 || vk == 3) {
                                int i1p = -1, i2p = -1, mp = -1;
                                for (int ns = 0; ns < 16; ++ns) {
                                    double num = num0, den = den0; int i1 = 0, i2 = 0, mask = 0;
                                    if (uk == 3 && gu_plus > val) { num += ru_ovf * gu_plus; den += ru_ovf; mask |= 64; }
                                    if (vk == 3 && gv_plus > val) { num += rv_ovf * gv_plus; den += rv_ovf; mask |= 128; }
                                    if (uk == 1) {
                                        double gr = (gu_plus > val) ? (mask |= 1, 1.0) : gu_lo;
                                        num += gr * gu_plus; den += gr;
                                        if (bu_hi > 0) { double br = (gu_bud_v > val) ? (mask |= 2, bu_hi) : bu_lo; num += br * gu_bud_v; den += br; }
                                        if (v_in) { double d1 = gv_plus - val; if (d1 > d0_v) { mask |= 4; num += inv_u * gv_plus; den += inv_u; } else num += inv_u * d0_v; }
                                    }
                                    if (vk == 1) {
                                        double gr = (gv_plus > val) ? (mask |= 8, 1.0) : gv_lo;
                                        num += gr * gv_plus; den += gr;
                                        if (bv_hi > 0) { double br = (gv_bud_v > val) ? (mask |= 16, bv_hi) : bv_lo; num += br * gv_bud_v; den += br; }
                                        if (u_in) { double d1 = gu_plus - val; if (d1 > d0_u) { mask |= 32; num += inv_v * gu_plus; den += inv_v; } else num += inv_v * d0_u; }
                                    }
                                    if (u_fy) {
                                        double yp = gu_plus - val, bv = -1e30;
                                        const double* ku = &knots[(size_t)u * nk];
                                        for (int kk = 0; kk < nk; ++kk) { double cc = KNOTS[kk] * yp - ku[kk]; if (cc > bv) { bv = cc; i1 = kk; } }
                                        num += (cu - 1) * (KNOTS[i1] * gu_plus - ku[i1]); den += (cu - 1) * KNOTS[i1];
                                    }
                                    if (v_fy) {
                                        double yv = gv_plus - val, bv = -1e30;
                                        const double* kv = &knots[(size_t)v * nk];
                                        for (int kk = 0; kk < nk; ++kk) { double cc = KNOTS[kk] * yv - kv[kk]; if (cc > bv) { bv = cc; i2 = kk; } }
                                        num += (cv - 1) * (KNOTS[i2] * gv_plus - kv[i2]); den += (cv - 1) * KNOTS[i2];
                                    }
                                    double vn = num / den;
                                    if (i1 == i1p && i2 == i2p && mask == mp) { val = vn; break; }
                                    i1p = i1; i2p = i2; mp = mask; val = vn;
                                }
                            }
                            gn[L.g_off(u, v)] = val;
                        }
                    }
                }
            }

            // === STEP 3: H3 ===
            vector<double> mid23(L.n_e23), rad23(L.n_e23);
            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
            for (int m = 0; m < L.M3.N; ++m) {
                int gu = L.G2.idx[L.M3.D[m]][L.M3.C[m]], gup = L.G2.up(gu);
                vector<double> dg(NG);
                for (int w = 0; w < NG; ++w) dg[w] = gn[L.g_off(gup, w)] - gn[L.g_off(gu, w)];
                proj_1d_fp(L.PM23, dg.data(), &mid23[(size_t)m * N3], &rad23[(size_t)m * N3]);
            }
            auto eval_p4_nb = [&](int a_bud_e3, int m_m3, int b_e3) -> double {
                if (L.G3.kind[a_bud_e3] != 0) return 0.0;
                int Da = L.G3.D[a_bud_e3], ca = L.G3.C[a_bud_e3], Dm = L.M3.D[m_m3], cm = L.M3.C[m_m3];
                if (Da >= h.d4i || ca < 2 || ca > h.K4i || Dm >= h.d4i || cm < 2 || cm > h.K4i || Da + Dm > h.SP4) return 0.0;
                int ma = L.M4v.idx[Da][ca], mm = L.M4v.idx[Dm][cm], ce, num, den;
                L.G4.proj_int(L.G3, b_e3, ce, num, den);
                double w = (double)num / den;
                return w * p4[L.p4_off(ma, mm, J_e4, ce)] + (1.0 - w) * p4[L.p4_off(ma, mm, J_e4, L.G4.INF)];
            };
            auto eval_s3_nb = [&](int Dm, int cm_next, int a_e3, int b_e3) -> double {
                if (Dm >= h.dS3 || cm_next < 3 || cm_next > h.KS3) return 0.0;
                int mo = L.MS3.idx[Dm][cm_next], ca, na, da, cb, nb, db;
                L.G4.proj_int(L.G3, a_e3, ca, na, da); L.G4.proj_int(L.G3, b_e3, cb, nb, db);
                int ai[2]={ca, L.G4.INF}, bi[2]={cb, L.G4.INF};
                double aw[2]={(double)na/da, 1.0-(double)na/da}, bw[2]={(double)nb/db, 1.0-(double)nb/db};
                double res = 0.0;
                for (int p=0;p<2;++p) if (aw[p]>0) for (int q=0;q<2;++q) if (bw[q]>0)
                    res += aw[p] * bw[q] * s3[L.s3_off(mo, ai[p], bi[q], J_e4)];
                return res;
            };
            vector<double>& h3n = h3;
            fill(h3n.begin(), h3n.end(), 0.0);
            int maxR3 = L.G3.maxR, rm_max = (int)h.d3m - 1;
            for (int rsum = rm_max + 2 * maxR3; rsum >= 6; --rsum) {
                vector<array<int, 3>> tasks;
                for (int rm = 2; rm <= rm_max; ++rm)
                    for (int ra = 2; ra <= maxR3; ++ra) {
                        int rb = rsum - rm - ra;
                        if (rb >= ra && rb <= maxR3 && !L.G3.by_r[ra].empty() && !L.G3.by_r[rb].empty())
                            tasks.push_back({rm, ra, rb});
                    }
                #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1 && tasks.size() >= 8)
                for (size_t ti = 0; ti < tasks.size(); ++ti) {
                    auto [rm, ra, rb] = tasks[ti];
                    for (int cm = min(rm, (int)h.K3m); cm >= 2; --cm) {
                        int m = L.M3.idx[rm][cm], Lm = rm - cm;
                        double im = 1.0 / rm, bm = (double)Lm / rm, rm_c = r1[L.G2.idx[rm][cm]];
                        int m_up = L.M3.up(m), m_bud = L.M3.bud(m);
                        int m_m4 = (rm < h.d4i && cm <= h.K4i) ? L.M4v.idx[rm][cm] : -1;
                        bool m_in_s3 = (rm < h.dS3 && cm >= 3 && cm <= h.KS3);
                        double rm_hi = (!m_in_s3) ? 0.5 * (cm - 2) : 0.0;

                        for (int ia = (int)L.G3.by_r[ra].size() - 1; ia >= 0; --ia) {
                            int a = L.G3.by_r[ra][ia];
                            int ak = L.G3.kind[a], Da = L.G3.D[a], ca = L.G3.C[a], La = L.G3.L[a];
                            bool a_in = L.G3.is_inner(a);
                            int a_up = L.G3.up(a), a_bud = L.G3.bud(a);
                            double ba_lo = b_lo_end(L.G3, a), ba_hi = b_hi_end(L.G3, a), ga_lo = g_lo_end(L.G3, a, im), inv_a = invD_end(L.G3, a);
                            double mid_a = mid23[(size_t)m * N3 + a], rad_a = rad23[(size_t)m * N3 + a];
                            bool a_to_p4 = (m_m4 >= 0 && ak == 0 && Da < h.d4i && ca >= 2 && ca <= h.K4i && Da + rm <= h.SP4);
                            double ra_c = (ak == 0) ? r1[L.G2.idx[Da][ca]] : 0.0;
                            double ra_hi = (ak == 0 && !a_to_p4) ? 0.5 * (ca - 1) : r_ovf_end(L.G3, a);

                            int ib_lo = (ra == rb) ? ia : 0;
                            for (int ib = (int)L.G3.by_r[rb].size() - 1; ib >= ib_lo; --ib) {
                                int b = L.G3.by_r[rb][ib];
                                int bk = L.G3.kind[b], Db = L.G3.D[b], cb = L.G3.C[b], Lb = L.G3.L[b];
                                bool b_in = L.G3.is_inner(b);
                                int b_up = L.G3.up(b), b_bud = L.G3.bud(b);
                                double bb_lo = b_lo_end(L.G3, b), bb_hi = b_hi_end(L.G3, b), gb_lo = g_lo_end(L.G3, b, im), inv_b = invD_end(L.G3, b);
                                double mid_b = mid23[(size_t)m * N3 + b], rad_b = rad23[(size_t)m * N3 + b];
                                bool b_to_p4 = (m_m4 >= 0 && bk == 0 && Db < h.d4i && cb >= 2 && cb <= h.K4i && Db + rm <= h.SP4);
                                double rb_c = (bk == 0) ? r1[L.G2.idx[Db][cb]] : 0.0;
                                double rb_hi = (bk == 0 && !b_to_p4) ? 0.5 * (cb - 1) : r_ovf_end(L.G3, b);

                                double src = (inv_a - rm_c) * mid_b + fabs(inv_a - rm_c) * rad_b
                                           + (inv_b - rm_c) * mid_a + fabs(inv_b - rm_c) * rad_a;

                                double ha_up = h3n[L.h3_off(m, a_up, b)], hb_up = h3n[L.h3_off(m, a, b_up)];
                                double hm_up = (m_up >= 0) ? h3n[L.h3_off(m_up, a, b)] : 0.0;
                                double ha_bud = (ba_hi > 0) ? h3n[L.h3_off(m, a_bud, b)] : 0.0;
                                double hb_bud = (bb_hi > 0) ? h3n[L.h3_off(m, a, b_bud)] : 0.0;
                                double hm_bud = (m_bud >= 0) ? h3n[L.h3_off(m_bud, a, b)] : 0.0;

                                double a_nb = (ak == 0 && La >= 1) ? eval_p4_nb(a_bud, m, b) : 0.0;
                                double b_nb = (bk == 0 && Lb >= 1) ? eval_p4_nb(b_bud, m, a) : 0.0;
                                double m_nb = (Lm >= 1) ? eval_s3_nb(rm, cm + 1, a, b) : 0.0;

                                double num0 = src, den0 = lm;
                                if (ak == 0) { double r = La + im + (a_to_p4 ? (ca - 1) * ra_c : 0.0); num0 += r * ha_up; den0 += r; if (ba_lo > 0) { num0 += ba_lo * (ha_bud + a_nb); den0 += ba_lo; } }
                                else if (ak == 3) { num0 += im * ha_up; den0 += im; }
                                if (bk == 0) { double r = Lb + im + (b_to_p4 ? (cb - 1) * rb_c : 0.0); num0 += r * hb_up; den0 += r; if (bb_lo > 0) { num0 += bb_lo * (hb_bud + b_nb); den0 += bb_lo; } }
                                else if (bk == 3) { num0 += im * hb_up; den0 += im; }
                                { double r = Lm + (a_in ? inv_a : 0.0) + (b_in ? inv_b : 0.0) + (m_in_s3 ? (cm - 2) * rm_c : 0.0); num0 += r * hm_up; den0 += r; if (bm > 0) { num0 += bm * (hm_bud + m_nb); den0 += bm; } }

                                double d0_am = (ak == 1) ? (((m_up >= 0) ? h3n[L.h3_off(m_up, L.G3.INF, b)] : 0.0) - h3n[L.h3_off(m, L.G3.INF, b)]) : 0.0;
                                double d0_bm = (bk == 1) ? (((m_up >= 0) ? h3n[L.h3_off(m_up, a, L.G3.INF)] : 0.0) - h3n[L.h3_off(m, a, L.G3.INF)]) : 0.0;

                                double val = num0 / den0; int mp = -1;
                                for (int ns = 0; ns < 12; ++ns) {
                                    double num = num0, den = den0; int mask = 0;
                                    if (ra_hi > 0 && ha_up > val) { num += ra_hi * ha_up; den += ra_hi; mask |= 1; }
                                    if (rb_hi > 0 && hb_up > val) { num += rb_hi * hb_up; den += rb_hi; mask |= 2; }
                                    if (rm_hi > 0 && hm_up > val) { num += rm_hi * hm_up; den += rm_hi; mask |= 4; }
                                    if (ak == 1) {
                                        double gr = (ha_up > val) ? (mask |= 8, 1.0) : ga_lo;
                                        num += gr * ha_up; den += gr;
                                        if (ba_hi > 0) { double br = (ha_bud > val) ? (mask |= 16, ba_hi) : ba_lo; num += br * ha_bud; den += br; }
                                        double d1 = hm_up - val;
                                        if (d1 > d0_am) { mask |= 32; num += inv_a * hm_up; den += inv_a; } else num += inv_a * d0_am;
                                    }
                                    if (bk == 1) {
                                        double gr = (hb_up > val) ? (mask |= 64, 1.0) : gb_lo;
                                        num += gr * hb_up; den += gr;
                                        if (bb_hi > 0) { double br = (hb_bud > val) ? (mask |= 128, bb_hi) : bb_lo; num += br * hb_bud; den += br; }
                                        double d1 = hm_up - val;
                                        if (d1 > d0_bm) { mask |= 256; num += inv_b * hm_up; den += inv_b; } else num += inv_b * d0_bm;
                                    }
                                    double vn = num / den; if (mask == mp) { val = vn; break; } mp = mask; val = vn;
                                }
                                h3n[L.h3_off(m, a, b)] = val;
                            }
                        }
                    }
                }
            }

            // === STEP 4A: P4 ===
            vector<double> mid34_p(L.n_e34p), rad34_p(L.n_e34p);
            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
            for (int pid = 0; pid < L.NP4_ord; ++pid) {
                auto [b, cc] = L.p4_ord[pid];
                int b_e3 = L.G3.idx[L.M4v.D[b]][L.M4v.C[b]], b_e3up = L.G3.up(b_e3);
                int c_m3 = L.M3.idx[L.M4v.D[cc]][L.M4v.C[cc]];
                vector<double> dh(N3);
                for (int d3 = 0; d3 < N3; ++d3) dh[d3] = h3n[L.h3_off(c_m3, b_e3up, d3)] - h3n[L.h3_off(c_m3, b_e3, d3)];
                proj_1d_fp(L.PM34, dh.data(), &mid34_p[(size_t)pid * N4], &rad34_p[(size_t)pid * N4]);
            }
            auto eval_p5_nb_from_p4 = [&](int a_bud_e4, int b_m4, int c_m4, int d_e4) -> double {
                if (L.G4.kind[a_bud_e4] != 0) return 0.0;
                int Da = L.G4.D[a_bud_e4], ca = L.G4.C[a_bud_e4];
                int Db = L.M4v.D[b_m4], cb = L.M4v.C[b_m4], Dc = L.M4v.D[c_m4], cc = L.M4v.C[c_m4];
                if (Da >= h.d5i || ca < 2 || ca > h.K5i || Db >= h.d5i || cb > h.K5i || Dc >= h.d5i || cc > h.K5i) return 0.0;
                if (Da + Db > h.SP4 || Db + Dc > h.SP4 || Da + Db + Dc > h.SP5) return 0.0;
                int ia = L.I5_2.idx[Da][ca], ib = L.I5_2.idx[Db][cb], ic = L.I5_2.idx[Dc][cc], ce, num, den;
                L.G5.proj_int(L.G4, d_e4, ce, num, den);
                double w = (double)num / den;
                return w * p5[L.p5_off(ia, ib, ic, J_e5, ce)] + (1.0 - w) * p5[L.p5_off(ia, ib, ic, J_e5, L.G5.INF)];
            };
            auto eval_b5_nb_from_p4 = [&](int Do, int co_next, int a_e4, int c_m4, int d_e4) -> double {
                if (Do >= h.d5o || co_next < 3 || co_next > h.K5o) return 0.0;
                int Dc = L.M4v.D[c_m4], cc = L.M4v.C[c_m4];
                if (Dc >= h.d5c || cc > h.K5c || Do + Dc > h.SB5) return 0.0;
                int bid = L.b5_int_map[(size_t)L.I5_3.idx[Do][co_next] * L.I5_2c.N + L.I5_2c.idx[Dc][cc]];
                int ca, na, da, cd, nd, dd;
                L.G5.proj_int(L.G4, a_e4, ca, na, da); L.G5.proj_int(L.G4, d_e4, cd, nd, dd);
                int ai[2]={ca, L.G5.INF}, di[2]={cd, L.G5.INF};
                double aw[2]={(double)na/da, 1.0-(double)na/da}, dw[2]={(double)nd/dd, 1.0-(double)nd/dd};
                double res = 0.0;
                for (int p=0;p<2;++p) if (aw[p]>0) for (int q=0;q<2;++q) if (dw[q]>0)
                    res += aw[p] * dw[q] * b5[L.b5_off(bid, ai[p], J_e5, di[q])];
                return res;
            };
            vector<double>& p4n = p4;
            fill(p4n.begin(), p4n.end(), 0.0);
            int ri4 = (int)h.d4i - 1;
            for (int rsum = (int)h.SP4 + 2 * maxR4; rsum >= 8; --rsum) {
                vector<array<int, 4>> tasks;
                for (int rb = 2; rb <= ri4; ++rb)
                    for (int rc = rb; rc <= min(ri4, (int)h.SP4 - rb); ++rc)
                        for (int ra = 2; ra <= maxR4; ++ra) {
                            int rd = rsum - rb - rc - ra;
                            if (rd >= 2 && rd <= maxR4 && !L.G4.by_r[ra].empty() && !L.G4.by_r[rd].empty()) {
                                if (rb == rc && ra > rd) continue;
                                tasks.push_back({ra, rb, rc, rd});
                            }
                        }
                #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1 && tasks.size() >= 8)
                for (size_t ti = 0; ti < tasks.size(); ++ti) {
                    auto [ra, rb, rc, rd] = tasks[ti];
                    for (int cb = min(rb, (int)h.K4i); cb >= 2; --cb) {
                        int b = L.M4v.idx[rb][cb], Lb = rb - cb;
                        double ib_inv = 1.0 / rb, bb = (double)Lb / rb, r_b = r1[L.G2.idx[rb][cb]];
                        int b_up = L.M4v.up(b), b_bud = L.M4v.bud(b);
                        int cc_lo = (rb == rc && ra == rd) ? cb : 2;
                        for (int cc = min(rc, (int)h.K4i); cc >= cc_lo; --cc) {
                            int c = L.M4v.idx[rc][cc], Lc = rc - cc;
                            double ic_inv = 1.0 / rc, bc = (double)Lc / rc, r_c = r1[L.G2.idx[rc][cc]];
                            int c_up = L.M4v.up(c), c_bud = L.M4v.bud(c);

                            int pid = L.p4_ord_map[(size_t)b * L.M4v.N + c], pid_rev = L.p4_ord_map[(size_t)c * L.M4v.N + b];
                            bool has_bup = (b_up >= 0 && rb + 1 + rc <= h.SP4);
                            bool has_cup = (c_up >= 0 && rb + rc + 1 <= h.SP4);
                            bool cov_b = (rb < h.d5o && cb >= 3 && cb <= h.K5o && rc < h.d5c && cc <= h.K5c && rb + rc <= h.SB5);
                            bool cov_c = (rc < h.d5o && cc >= 3 && cc <= h.K5o && rb < h.d5c && cb <= h.K5c && rb + rc <= h.SB5);
                            bool bc_in_p5 = (rb < h.d5i && cb <= h.K5i && rc < h.d5i && cc <= h.K5i);
                            double rb_hi = (!cov_b) ? 0.5 * (cb - 2) : 0.0, rc_hi = (!cov_c) ? 0.5 * (cc - 2) : 0.0;

                            for (int ia = (int)L.G4.by_r[ra].size() - 1; ia >= 0; --ia) {
                                int a = L.G4.by_r[ra][ia];
                                int ak = L.G4.kind[a], Da = L.G4.D[a], ca = L.G4.C[a], La = L.G4.L[a];
                                bool a_in = L.G4.is_inner(a);
                                int a_up = L.G4.up(a), a_bud = L.G4.bud(a);
                                double ba_lo = b_lo_end(L.G4, a), ba_hi = b_hi_end(L.G4, a), ga_lo = g_lo_end(L.G4, a, ib_inv), inv_a = invD_end(L.G4, a);
                                bool cov_a = (bc_in_p5 && ak == 0 && Da < h.d5i && ca >= 2 && ca <= h.K5i && Da + rb <= h.SP4 && Da + rb + rc <= h.SP5);
                                double r_a = (ak == 0) ? r1[L.G2.idx[Da][ca]] : 0.0;
                                double ra_hi = (ak == 0 && !cov_a) ? 0.5 * (ca - 1) : r_ovf_end(L.G4, a);

                                int id_lo = (rb == rc && cb == cc && ra == rd) ? ia : 0;
                                for (int id = (int)L.G4.by_r[rd].size() - 1; id >= id_lo; --id) {
                                    int d = L.G4.by_r[rd][id];
                                    int dk = L.G4.kind[d], Dd = L.G4.D[d], cd = L.G4.C[d], Ld = L.G4.L[d];
                                    bool d_in = L.G4.is_inner(d);
                                    int d_up = L.G4.up(d), d_bud = L.G4.bud(d);
                                    double bd_lo = b_lo_end(L.G4, d), bd_hi = b_hi_end(L.G4, d), gd_lo = g_lo_end(L.G4, d, ic_inv), inv_d = invD_end(L.G4, d);
                                    bool cov_d = (bc_in_p5 && dk == 0 && Dd < h.d5i && cd >= 2 && cd <= h.K5i && Dd + rc <= h.SP4 && rb + rc + Dd <= h.SP5);
                                    double r_d = (dk == 0) ? r1[L.G2.idx[Dd][cd]] : 0.0;
                                    double rd_hi = (dk == 0 && !cov_d) ? 0.5 * (cd - 1) : r_ovf_end(L.G4, d);

                                    size_t o_bc = (size_t)pid * N4 + d, o_cb = (size_t)pid_rev * N4 + a;
                                    double src = (inv_a - r_b) * mid34_p[o_bc] + fabs(inv_a - r_b) * rad34_p[o_bc]
                                               + (inv_d - r_c) * mid34_p[o_cb] + fabs(inv_d - r_c) * rad34_p[o_cb];

                                    double pa_up = p4n[L.p4_off(b, c, a_up, d)], pd_up = p4n[L.p4_off(b, c, a, d_up)];
                                    double pb_up = has_bup ? p4n[L.p4_off(b_up, c, a, d)] : 0.0;
                                    double pc_up = has_cup ? p4n[L.p4_off(b, c_up, a, d)] : 0.0;
                                    double pa_bud = (ba_hi > 0) ? p4n[L.p4_off(b, c, a_bud, d)] : 0.0;
                                    double pd_bud = (bd_hi > 0) ? p4n[L.p4_off(b, c, a, d_bud)] : 0.0;
                                    double pb_bud = (b_bud >= 0) ? p4n[L.p4_off(b_bud, c, a, d)] : 0.0;
                                    double pc_bud = (c_bud >= 0) ? p4n[L.p4_off(b, c_bud, a, d)] : 0.0;

                                    double a_nb = (ak == 0 && La >= 1) ? eval_p5_nb_from_p4(a_bud, b, c, d) : 0.0;
                                    double d_nb = (dk == 0 && Ld >= 1) ? eval_p5_nb_from_p4(d_bud, c, b, a) : 0.0;
                                    double b_nb = (Lb >= 1) ? eval_b5_nb_from_p4(rb, cb + 1, a, c, d) : 0.0;
                                    double c_nb = (Lc >= 1) ? eval_b5_nb_from_p4(rc, cc + 1, d, b, a) : 0.0;

                                    double num0 = src, den0 = lm;
                                    if (ak == 0) { double r = La + ib_inv + (cov_a ? (ca - 1) * r_a : 0.0); num0 += r * pa_up; den0 += r; if (ba_lo > 0) { num0 += ba_lo * (pa_bud + a_nb); den0 += ba_lo; } }
                                    else if (ak == 3) { num0 += ib_inv * pa_up; den0 += ib_inv; }
                                    if (dk == 0) { double r = Ld + ic_inv + (cov_d ? (cd - 1) * r_d : 0.0); num0 += r * pd_up; den0 += r; if (bd_lo > 0) { num0 += bd_lo * (pd_bud + d_nb); den0 += bd_lo; } }
                                    else if (dk == 3) { num0 += ic_inv * pd_up; den0 += ic_inv; }
                                    { double r = Lb + (a_in ? inv_a : 0.0) + ic_inv + (cov_b ? (cb - 2) * r_b : 0.0); num0 += r * pb_up; den0 += r; if (bb > 0) { num0 += bb * (pb_bud + b_nb); den0 += bb; } }
                                    { double r = Lc + ib_inv + (d_in ? inv_d : 0.0) + (cov_c ? (cc - 2) * r_c : 0.0); num0 += r * pc_up; den0 += r; if (bc > 0) { num0 += bc * (pc_bud + c_nb); den0 += bc; } }

                                    double d0_ab = (ak == 1) ? ((has_bup ? p4n[L.p4_off(b_up, c, L.G4.INF, d)] : 0.0) - p4n[L.p4_off(b, c, L.G4.INF, d)]) : 0.0;
                                    double d0_dc = (dk == 1) ? ((has_cup ? p4n[L.p4_off(b, c_up, a, L.G4.INF)] : 0.0) - p4n[L.p4_off(b, c, a, L.G4.INF)]) : 0.0;

                                    double val = num0 / den0; int mp = -1;
                                    for (int ns = 0; ns < 12; ++ns) {
                                        double num = num0, den = den0; int mask = 0;
                                        if (ra_hi > 0 && pa_up > val) { num += ra_hi * pa_up; den += ra_hi; mask |= 1; }
                                        if (rd_hi > 0 && pd_up > val) { num += rd_hi * pd_up; den += rd_hi; mask |= 2; }
                                        if (rb_hi > 0 && pb_up > val) { num += rb_hi * pb_up; den += rb_hi; mask |= 4; }
                                        if (rc_hi > 0 && pc_up > val) { num += rc_hi * pc_up; den += rc_hi; mask |= 8; }
                                        if (ak == 1) {
                                            double gr = (pa_up > val) ? (mask |= 16, 1.0) : ga_lo;
                                            num += gr * pa_up; den += gr;
                                            if (ba_hi > 0) { double br = (pa_bud > val) ? (mask |= 32, ba_hi) : ba_lo; num += br * pa_bud; den += br; }
                                            double d1 = pb_up - val;
                                            if (d1 > d0_ab) { mask |= 64; num += inv_a * pb_up; den += inv_a; } else num += inv_a * d0_ab;
                                        }
                                        if (dk == 1) {
                                            double gr = (pd_up > val) ? (mask |= 128, 1.0) : gd_lo;
                                            num += gr * pd_up; den += gr;
                                            if (bd_hi > 0) { double br = (pd_bud > val) ? (mask |= 256, bd_hi) : bd_lo; num += br * pd_bud; den += br; }
                                            double d1 = pc_up - val;
                                            if (d1 > d0_dc) { mask |= 512; num += inv_d * pc_up; den += inv_d; } else num += inv_d * d0_dc;
                                        }
                                        double vn = num / den; if (mask == mp) { val = vn; break; } mp = mask; val = vn;
                                    }
                                    p4n[L.p4_off(b, c, a, d)] = val;
                                }
                            }
                        }
                    }
                }
            }

            // === STEP 4B: S3 ===
            int NS3m = L.MS3.N;
            vector<vector<double>> mid34_s(NS3m), rad34_s(NS3m);
            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
            for (int o = 0; o < NS3m; ++o) {
                int om3 = L.M3.idx[L.MS3.D[o]][L.MS3.C[o]], om3_up = L.M3.up(om3);
                vector<double> Mf((size_t)N3 * N3);
                for (int u = 0; u < N3; ++u) for (int v = 0; v < N3; ++v)
                    Mf[(size_t)u * N3 + v] = ((om3_up >= 0) ? h3n[L.h3_off(om3_up, u, v)] : 0.0) - h3n[L.h3_off(om3, u, v)];
                proj_2d_gen_fp(L.PM34, Mf, mid34_s[o], rad34_s[o]);
            }
            auto eval_b5_nb_from_s3 = [&](int o_ms3, int c_bud_e4, int a_e4, int b_e4) -> double {
                if (L.G4.kind[c_bud_e4] != 0) return 0.0;
                int Do = L.MS3.D[o_ms3], co = L.MS3.C[o_ms3], Dc = L.G4.D[c_bud_e4], cc = L.G4.C[c_bud_e4];
                if (Do >= h.d5o || co > h.K5o || Dc >= h.d5c || cc < 2 || cc > h.K5c || Do + Dc > h.SB5) return 0.0;
                int bid = L.b5_int_map[(size_t)L.I5_3.idx[Do][co] * L.I5_2c.N + L.I5_2c.idx[Dc][cc]];
                int ca, na, da, cb, nb, db;
                L.G5.proj_int(L.G4, a_e4, ca, na, da); L.G5.proj_int(L.G4, b_e4, cb, nb, db);
                int ai[2]={ca, L.G5.INF}, bi[2]={cb, L.G5.INF};
                double aw[2]={(double)na/da, 1.0-(double)na/da}, bw[2]={(double)nb/db, 1.0-(double)nb/db};
                double res = 0.0;
                for (int p=0;p<2;++p) if (aw[p]>0) for (int q=0;q<2;++q) if (bw[q]>0)
                    res += aw[p] * bw[q] * b5[L.b5_off(bid, ai[p], bi[q], J_e5)];
                return res;
            };
            auto eval_s4_nb_from_s3 = [&](int Do, int co_next, int a_e4, int b_e4, int c_e4) -> double {
                if (Do >= h.dS4 || co_next < 4 || co_next > h.KS4) return 0.0;
                int o4 = L.MS4.idx[Do][co_next], ca, na, da, cb, nb, db, cc, nc, dc;
                L.G5.proj_int(L.G4, a_e4, ca, na, da); L.G5.proj_int(L.G4, b_e4, cb, nb, db); L.G5.proj_int(L.G4, c_e4, cc, nc, dc);
                int ai[2]={ca, L.G5.INF}, bi[2]={cb, L.G5.INF}, ci[2]={cc, L.G5.INF};
                double aw[2]={(double)na/da, 1.0-(double)na/da}, bw[2]={(double)nb/db, 1.0-(double)nb/db}, cw[2]={(double)nc/dc, 1.0-(double)nc/dc};
                double res = 0.0;
                for (int p=0;p<2;++p) if (aw[p]>0) for (int q=0;q<2;++q) if (bw[q]>0) for (int r=0;r<2;++r) if (cw[r]>0)
                    res += aw[p] * bw[q] * cw[r] * s4[L.s4_off(o4, ai[p], bi[q], ci[r], J_e5)];
                return res;
            };
            vector<double>& s3n = s3;
            fill(s3n.begin(), s3n.end(), 0.0);
            int rs3_max = (int)h.dS3 - 1;
            for (int rsum = rs3_max + 3 * maxR4; rsum >= 9; --rsum) {
                for (int ro = min(rs3_max, rsum - 6); ro >= 3; --ro) {
                    int rt = rsum - ro;
                    if (rt < 6 || rt > 3 * maxR4) continue;
                    for (int co = min(ro, (int)h.KS3); co >= 3; --co) {
                        int o = L.MS3.idx[ro][co], Lo = ro - co;
                        double io = 1.0 / ro, bo = (double)Lo / ro, r_o = r1[L.G2.idx[ro][co]];
                        int o_up = L.MS3.up(o), o_bud = L.MS3.bud(o);
                        bool cov_o = (ro < h.dS4 && co >= 4 && co <= h.KS4);
                        double ro_hi = (!cov_o) ? 0.5 * (co - 3) : 0.0;

                        for (int csum = (int)tri4_rc[rt].size() - 1; csum >= 0; --csum) {
                            const auto& bucket = tri4_rc[rt][csum];
                            if (bucket.empty()) continue;
                            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1 && bucket.size() > 1500)
                            for (size_t idx_b = 0; idx_b < bucket.size(); ++idx_b) {
                                int tid = bucket[idx_b];
                                auto [a, b, c] = L.tri_G4[tid];
                                int e[3] = {a, b, c};
                                double src = 0.0;
                                for (int k = 0; k < 3; ++k) {
                                    int d1 = e[(k+1)%3], d2 = e[(k+2)%3], don = e[k];
                                    double idon = invD_end(L.G4, don);
                                    src += (idon - r_o) * mid34_s[o][(size_t)d1 * N4 + d2]
                                         + fabs(idon - r_o) * rad34_s[o][(size_t)d1 * N4 + d2];
                                }
                                double ho_up = (o_up >= 0) ? s3n[(size_t)o_up * L.N_G4_tri + tid] : 0.0;
                                double ho_bud = (o_bud >= 0) ? s3n[(size_t)o_bud * L.N_G4_tri + tid] : 0.0;
                                double o_nb = (Lo >= 1) ? eval_s4_nb_from_s3(ro, co + 1, a, b, c) : 0.0;

                                double num0 = src, den0 = lm, rm_b = Lo + (cov_o ? (co - 3) * r_o : 0.0);
                                double he_up[3], he_bud[3], r_ext[3], b_lo[3], b_hi[3], g_lo[3], d0_eo[3];
                                int ek[3];
                                for (int k = 0; k < 3; ++k) {
                                    int u = e[k], u1 = e[(k+1)%3], u2 = e[(k+2)%3];
                                    ek[k] = L.G4.kind[u];
                                    int u_up = L.G4.up(u), u_bud = L.G4.bud(u);
                                    b_lo[k] = b_lo_end(L.G4, u); b_hi[k] = b_hi_end(L.G4, u); g_lo[k] = g_lo_end(L.G4, u, io);
                                    he_up[k] = s3n[L.s3_off(o, u_up, u1, u2)];
                                    he_bud[k] = (b_hi[k] > 0) ? s3n[L.s3_off(o, u_bud, u1, u2)] : 0.0;
                                    double u_nb = (ek[k] == 0 && L.G4.L[u] >= 1) ? eval_b5_nb_from_s3(o, u_bud, u1, u2) : 0.0;
                                    bool cov_u = (ro < h.d5o && co <= h.K5o && ek[k] == 0 && L.G4.D[u] < h.d5c && L.G4.C[u] >= 2 && L.G4.C[u] <= h.K5c && ro + L.G4.D[u] <= h.SB5);
                                    r_ext[k] = (ek[k] == 0 && !cov_u) ? 0.5 * (L.G4.C[u] - 1) : r_ovf_end(L.G4, u);
                                    if (ek[k] == 0) {
                                        double r = L.G4.L[u] + io + (cov_u ? (L.G4.C[u] - 1) * r1[L.G2.idx[L.G4.D[u]][L.G4.C[u]]] : 0.0);
                                        num0 += r * he_up[k]; den0 += r;
                                        if (b_lo[k] > 0) { num0 += b_lo[k] * (he_bud[k] + u_nb); den0 += b_lo[k]; }
                                        rm_b += invD_end(L.G4, u);
                                    } else if (ek[k] == 3) {
                                        num0 += io * he_up[k]; den0 += io; rm_b += invD_end(L.G4, u);
                                    }
                                    d0_eo[k] = (ek[k] == 1) ? (((o_up >= 0) ? s3n[L.s3_off(o_up, L.G4.INF, u1, u2)] : 0.0) - s3n[L.s3_off(o, L.G4.INF, u1, u2)]) : 0.0;
                                }
                                num0 += rm_b * ho_up; den0 += rm_b;
                                if (bo > 0) { num0 += bo * (ho_bud + o_nb); den0 += bo; }

                                double val = num0 / den0; int mp = -1;
                                for (int ns = 0; ns < 12; ++ns) {
                                    double num = num0, den = den0; int mask = 0;
                                    if (ro_hi > 0 && ho_up > val) { num += ro_hi * ho_up; den += ro_hi; mask |= 1; }
                                    for (int k = 0; k < 3; ++k) {
                                        int bit = 1 << (1 + 4 * k);
                                        if (r_ext[k] > 0 && he_up[k] > val) { num += r_ext[k] * he_up[k]; den += r_ext[k]; mask |= bit; }
                                        if (ek[k] == 1) {
                                            double gr = (he_up[k] > val) ? (mask |= (bit << 1), 1.0) : g_lo[k];
                                            num += gr * he_up[k]; den += gr;
                                            if (b_hi[k] > 0) { double br = (he_bud[k] > val) ? (mask |= (bit << 2), b_hi[k]) : b_lo[k]; num += br * he_bud[k]; den += br; }
                                            double d1 = ho_up - val, iv = invD_end(L.G4, e[k]);
                                            if (d1 > d0_eo[k]) { mask |= (bit << 3); num += iv * ho_up; den += iv; } else num += iv * d0_eo[k];
                                        }
                                    }
                                    double vn = num / den; if (mask == mp) { val = vn; break; } mp = mask; val = vn;
                                }
                                s3n[(size_t)o * L.N_G4_tri + tid] = val;
                            }
                        }
                    }
                }
            }

            // === STEP 5A: P5 ===
            vector<double> mid45_p(L.n_e45p), rad45_p(L.n_e45p);
            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
            for (int pid5 = 0; pid5 < L.NP5_ord; ++pid5) {
                auto [b, c, d] = L.p5_ord[pid5];
                int b_e4 = L.G4.idx[L.I5_2.D[b]][L.I5_2.C[b]], b_e4up = L.G4.up(b_e4);
                int c_m4 = L.M4v.idx[L.I5_2.D[c]][L.I5_2.C[c]], d_m4 = L.M4v.idx[L.I5_2.D[d]][L.I5_2.C[d]];
                vector<double> dp(N4);
                for (int e4 = 0; e4 < N4; ++e4) dp[e4] = p4n[L.p4_off(c_m4, d_m4, b_e4up, e4)] - p4n[L.p4_off(c_m4, d_m4, b_e4, e4)];
                proj_1d_fp(L.PM45, dp.data(), &mid45_p[(size_t)pid5 * N5], &rad45_p[(size_t)pid5 * N5]);
            }
            vector<double>& p5n = p5;
            fill(p5n.begin(), p5n.end(), 0.0);
            int ri5 = (int)h.d5i - 1, N52 = L.I5_2.N;
            for (int rsum = (int)h.SP5 + 2 * maxR5; rsum >= 10; --rsum) {
                vector<array<int, 5>> tasks;
                for (int rb = 2; rb <= ri5; ++rb)
                    for (int rc = 2; rc <= min(ri5, (int)h.SP4 - rb); ++rc)
                        for (int rd = rb; rd <= min({ri5, (int)h.SP4 - rc, (int)h.SP5 - rb - rc}); ++rd)
                            for (int ra = 2; ra <= maxR5; ++ra) {
                                int re = rsum - rb - rc - rd - ra;
                                if (re >= 2 && re <= maxR5 && !L.G5.by_r[ra].empty() && !L.G5.by_r[re].empty()) {
                                    if (rb == rd && ra > re) continue;
                                    tasks.push_back({ra, rb, rc, rd, re});
                                }
                            }
                #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1 && tasks.size() >= 8)
                for (size_t ti = 0; ti < tasks.size(); ++ti) {
                    auto [ra, rb, rc, rd, re] = tasks[ti];
                    for (int cb = min(rb, (int)h.K5i); cb >= 2; --cb) {
                        int b = L.I5_2.idx[rb][cb], b_up = L.I5_2.up(b), b_bud = L.I5_2.bud(b), Lb = rb - cb;
                        double ib_i = 1.0 / rb, bb = (double)Lb / rb, r_b = r1[L.G2.idx[rb][cb]];
                        for (int cc = min(rc, (int)h.K5i); cc >= 2; --cc) {
                            int c = L.I5_2.idx[rc][cc], c_up = L.I5_2.up(c), c_bud = L.I5_2.bud(c), Lc = rc - cc;
                            double ic_i = 1.0 / rc, bc = (double)Lc / rc;
                            int cd_lo = (rb == rd && ra == re) ? cb : 2;
                            for (int cd = min(rd, (int)h.K5i); cd >= cd_lo; --cd) {
                                int d = L.I5_2.idx[rd][cd], d_up = L.I5_2.up(d), d_bud = L.I5_2.bud(d), Ld = rd - cd;
                                double id_i = 1.0 / rd, bd = (double)Ld / rd, r_d = r1[L.G2.idx[rd][cd]];

                                int pid = L.p5_ord_map[((size_t)b * N52 + c) * N52 + d];
                                int pid_rev = L.p5_ord_map[((size_t)d * N52 + c) * N52 + b];
                                bool has_bup = (b_up >= 0 && rb + 1 + rc <= h.SP4 && rb + 1 + rc + rd <= h.SP5);
                                bool has_cup = (c_up >= 0 && rb + rc + 1 <= h.SP4 && rc + 1 + rd <= h.SP4 && rb + rc + 1 + rd <= h.SP5);
                                bool has_dup = (d_up >= 0 && rc + rd + 1 <= h.SP4 && rb + rc + rd + 1 <= h.SP5);

                                for (int ia = (int)L.G5.by_r[ra].size() - 1; ia >= 0; --ia) {
                                    int a = L.G5.by_r[ra][ia];
                                    int ak = L.G5.kind[a], ca = L.G5.C[a], La = L.G5.L[a];
                                    bool a_in = L.G5.is_inner(a);
                                    int a_up = L.G5.up(a), a_bud = L.G5.bud(a);
                                    double ba_lo = b_lo_end(L.G5, a), ba_hi = b_hi_end(L.G5, a), ga_lo = g_lo_end(L.G5, a, ib_i), inv_a = invD_end(L.G5, a);

                                    int ie_lo = (rb == rd && cb == cd && ra == re) ? ia : 0;
                                    for (int ie = (int)L.G5.by_r[re].size() - 1; ie >= ie_lo; --ie) {
                                        int e = L.G5.by_r[re][ie];
                                        int ek = L.G5.kind[e], ce = L.G5.C[e], Le = L.G5.L[e];
                                        bool e_in = L.G5.is_inner(e);
                                        int e_up = L.G5.up(e), e_bud = L.G5.bud(e);
                                        double be_lo = b_lo_end(L.G5, e), be_hi = b_hi_end(L.G5, e), ge_lo = g_lo_end(L.G5, e, id_i), inv_e = invD_end(L.G5, e);

                                        size_t o_f = (size_t)pid * N5 + e, o_r = (size_t)pid_rev * N5 + a;
                                        double src = (inv_a - r_b) * mid45_p[o_f] + fabs(inv_a - r_b) * rad45_p[o_f]
                                                   + (inv_e - r_d) * mid45_p[o_r] + fabs(inv_e - r_d) * rad45_p[o_r];

                                        double up[5] = {
                                            p5n[L.p5_off(b, c, d, a_up, e)],
                                            has_bup ? p5n[L.p5_off(b_up, c, d, a, e)] : 0.0,
                                            has_cup ? p5n[L.p5_off(b, c_up, d, a, e)] : 0.0,
                                            has_dup ? p5n[L.p5_off(b, c, d_up, a, e)] : 0.0,
                                            p5n[L.p5_off(b, c, d, a, e_up)]
                                        };
                                        double bd_a = (ba_hi > 0) ? p5n[L.p5_off(b, c, d, a_bud, e)] : 0.0;
                                        double bd_e = (be_hi > 0) ? p5n[L.p5_off(b, c, d, a, e_bud)] : 0.0;
                                        double bd_b = (b_bud >= 0) ? p5n[L.p5_off(b_bud, c, d, a, e)] : 0.0;
                                        double bd_c = (c_bud >= 0) ? p5n[L.p5_off(b, c_bud, d, a, e)] : 0.0;
                                        double bd_d = (d_bud >= 0) ? p5n[L.p5_off(b, c, d_bud, a, e)] : 0.0;

                                        double num0 = src, den0 = lm;
                                        if (ak == 0) { double r = La + ib_i; num0 += r * up[0]; den0 += r; if (ba_lo > 0) { num0 += ba_lo * bd_a; den0 += ba_lo; } }
                                        else if (ak == 3) { num0 += ib_i * up[0]; den0 += ib_i; }
                                        { double r = Lb + (a_in ? inv_a : 0.0) + ic_i; num0 += r * up[1]; den0 += r; if (bb > 0) { num0 += bb * bd_b; den0 += bb; } }
                                        { double r = Lc + ib_i + id_i; num0 += r * up[2]; den0 += r; if (bc > 0) { num0 += bc * bd_c; den0 += bc; } }
                                        { double r = Ld + ic_i + (e_in ? inv_e : 0.0); num0 += r * up[3]; den0 += r; if (bd > 0) { num0 += bd * bd_d; den0 += bd; } }
                                        if (ek == 0) { double r = Le + id_i; num0 += r * up[4]; den0 += r; if (be_lo > 0) { num0 += be_lo * bd_e; den0 += be_lo; } }
                                        else if (ek == 3) { num0 += id_i * up[4]; den0 += id_i; }

                                        double r_ext[5] = {
                                            (ak == 0) ? 0.5 * (ca - 1) : r_ovf_end(L.G5, a),
                                            0.5 * (cb - 2), 0.5 * (cc - 2), 0.5 * (cd - 2),
                                            (ek == 0) ? 0.5 * (ce - 1) : r_ovf_end(L.G5, e)
                                        };
                                        double d0_ab = (ak == 1) ? ((has_bup ? p5n[L.p5_off(b_up, c, d, L.G5.INF, e)] : 0.0) - p5n[L.p5_off(b, c, d, L.G5.INF, e)]) : 0.0;
                                        double d0_ed = (ek == 1) ? ((has_dup ? p5n[L.p5_off(b, c, d_up, a, L.G5.INF)] : 0.0) - p5n[L.p5_off(b, c, d, a, L.G5.INF)]) : 0.0;

                                        double val = num0 / den0; int mp = -1;
                                        for (int ns = 0; ns < 12; ++ns) {
                                            double num = num0, den = den0; int mask = 0;
                                            for (int k = 0; k < 5; ++k)
                                                if (r_ext[k] > 0 && up[k] > val) { num += r_ext[k] * up[k]; den += r_ext[k]; mask |= (1 << k); }
                                            if (ak == 1) {
                                                double gr = (up[0] > val) ? (mask |= 32, 1.0) : ga_lo;
                                                num += gr * up[0]; den += gr;
                                                if (ba_hi > 0) { double br = (bd_a > val) ? (mask |= 64, ba_hi) : ba_lo; num += br * bd_a; den += br; }
                                                double d1 = up[1] - val;
                                                if (d1 > d0_ab) { mask |= 128; num += inv_a * up[1]; den += inv_a; } else num += inv_a * d0_ab;
                                            }
                                            if (ek == 1) {
                                                double gr = (up[4] > val) ? (mask |= 256, 1.0) : ge_lo;
                                                num += gr * up[4]; den += gr;
                                                if (be_hi > 0) { double br = (bd_e > val) ? (mask |= 512, be_hi) : be_lo; num += br * bd_e; den += br; }
                                                double d1 = up[3] - val;
                                                if (d1 > d0_ed) { mask |= 1024; num += inv_e * up[3]; den += inv_e; } else num += inv_e * d0_ed;
                                            }
                                            double vn = num / den; if (mask == mp) { val = vn; break; } mp = mask; val = vn;
                                        }
                                        p5n[L.p5_off(b, c, d, a, e)] = val;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // === STEP 5B: B5 ===
            int NB5 = L.NB5_int;
            vector<vector<double>> mid_bp(NB5), rad_bp(NB5), mid_bs(NB5), rad_bs(NB5);
            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
            for (int bid = 0; bid < NB5; ++bid) {
                auto [o, c] = L.b5_int[bid];
                int Do = L.I5_3.D[o], co = L.I5_3.C[o], Dc = L.I5_2c.D[c], cc = L.I5_2c.C[c];
                int o_m4 = L.M4v.idx[Do][co], o_m4up = L.M4v.up(o_m4), c_m4 = L.M4v.idx[Dc][cc];
                bool has_up = (o_m4up >= 0 && Do + 1 + Dc <= h.SP4);
                vector<double> Mf((size_t)N4 * N4);
                for (int u = 0; u < N4; ++u) for (int v = 0; v < N4; ++v)
                    Mf[(size_t)u * N4 + v] = (has_up ? p4n[L.p4_off(o_m4up, c_m4, u, v)] : 0.0) - p4n[L.p4_off(o_m4, c_m4, u, v)];
                proj_2d_gen_fp(L.PM45, Mf, mid_bp[bid], rad_bp[bid]);

                int o_ms3 = L.MS3.idx[Do][co], c_e4 = L.G4.idx[Dc][cc], c_e4up = L.G4.up(c_e4);
                for (int u = 0; u < N4; ++u) for (int v = 0; v < N4; ++v)
                    Mf[(size_t)u * N4 + v] = s3n[L.s3_off(o_ms3, u, v, c_e4up)] - s3n[L.s3_off(o_ms3, u, v, c_e4)];
                proj_2d_gen_fp(L.PM45, Mf, mid_bs[bid], rad_bs[bid]);
            }
            vector<double>& b5n = b5;
            fill(b5n.begin(), b5n.end(), 0.0);
            int ro5 = (int)h.d5o - 1, rc5 = (int)h.d5c - 1;
            for (int rsum = (int)h.SB5 + 3 * maxR5; rsum >= 11; --rsum) {
                vector<array<int, 5>> tasks;
                for (int ro = 3; ro <= ro5; ++ro)
                    for (int rc = 2; rc <= min(rc5, (int)h.SB5 - ro); ++rc)
                        for (int ra = 2; ra <= maxR5; ++ra)
                            for (int rb = ra; rb <= maxR5; ++rb) {
                                int rd = rsum - ro - rc - ra - rb;
                                if (rd >= 2 && rd <= maxR5 && !L.G5.by_r[ra].empty() && !L.G5.by_r[rb].empty() && !L.G5.by_r[rd].empty())
                                    tasks.push_back({ro, rc, ra, rb, rd});
                            }
                #pragma omp parallel for schedule(dynamic) num_threads(nt) if(nt > 1 && tasks.size() >= 8)
                for (size_t ti = 0; ti < tasks.size(); ++ti) {
                    auto [ro, rc, ra, rb, rd] = tasks[ti];
                    for (int co = min(ro, (int)h.K5o); co >= 3; --co) {
                        int o = L.I5_3.idx[ro][co], o_up = L.I5_3.up(o), o_bud = L.I5_3.bud(o), Lo = ro - co;
                        double io = 1.0 / ro, bo = (double)Lo / ro, r_o = r1[L.G2.idx[ro][co]];
                        for (int cc = min(rc, (int)h.K5c); cc >= 2; --cc) {
                            int c = L.I5_2c.idx[rc][cc], c_up = L.I5_2c.up(c), c_bud = L.I5_2c.bud(c), Lc = rc - cc;
                            double ic = 1.0 / rc, bc = (double)Lc / rc, r_c = r1[L.G2.idx[rc][cc]];

                            int bid = L.b5_int_map[(size_t)o * L.I5_2c.N + c];
                            int b_oup = (o_up >= 0) ? L.b5_int_map[(size_t)o_up * L.I5_2c.N + c] : -1;
                            int b_cup = (c_up >= 0) ? L.b5_int_map[(size_t)o * L.I5_2c.N + c_up] : -1;
                            int b_obd = (o_bud >= 0) ? L.b5_int_map[(size_t)o_bud * L.I5_2c.N + c] : -1;
                            int b_cbd = (c_bud >= 0) ? L.b5_int_map[(size_t)o * L.I5_2c.N + c_bud] : -1;

                            for (int ia = (int)L.G5.by_r[ra].size() - 1; ia >= 0; --ia) {
                                int a = L.G5.by_r[ra][ia];
                                int ak = L.G5.kind[a], ca = L.G5.C[a], La = L.G5.L[a];
                                bool a_in = L.G5.is_inner(a);
                                int a_up = L.G5.up(a), a_bud = L.G5.bud(a);
                                double ba_lo = b_lo_end(L.G5, a), ba_hi = b_hi_end(L.G5, a), ga_lo = g_lo_end(L.G5, a, io), inv_a = invD_end(L.G5, a);

                                int ib_lo = (ra == rb) ? ia : 0;
                                for (int ib = (int)L.G5.by_r[rb].size() - 1; ib >= ib_lo; --ib) {
                                    int b = L.G5.by_r[rb][ib];
                                    int bk = L.G5.kind[b], cb = L.G5.C[b], Lb = L.G5.L[b];
                                    bool b_in = L.G5.is_inner(b);
                                    int b_up = L.G5.up(b), b_bud = L.G5.bud(b);
                                    double bb_lo = b_lo_end(L.G5, b), bb_hi = b_hi_end(L.G5, b), gb_lo = g_lo_end(L.G5, b, io), inv_b = invD_end(L.G5, b);

                                    double m_bs = mid_bs[bid][(size_t)a * N5 + b], r_bs = rad_bs[bid][(size_t)a * N5 + b];
                                    for (int id = (int)L.G5.by_r[rd].size() - 1; id >= 0; --id) {
                                        int d = L.G5.by_r[rd][id];
                                        int dk = L.G5.kind[d], cd = L.G5.C[d], Ld = L.G5.L[d];
                                        bool d_in = L.G5.is_inner(d);
                                        int d_up = L.G5.up(d), d_bud = L.G5.bud(d);
                                        double bd_lo = b_lo_end(L.G5, d), bd_hi = b_hi_end(L.G5, d), gd_lo = g_lo_end(L.G5, d, ic), inv_d = invD_end(L.G5, d);

                                        double src = (inv_b - r_o) * mid_bp[bid][(size_t)a * N5 + d] + fabs(inv_b - r_o) * rad_bp[bid][(size_t)a * N5 + d]
                                                   + (inv_a - r_o) * mid_bp[bid][(size_t)b * N5 + d] + fabs(inv_a - r_o) * rad_bp[bid][(size_t)b * N5 + d]
                                                   + (inv_d - r_c) * m_bs + fabs(inv_d - r_c) * r_bs;

                                        double up[5] = {
                                            (b_oup >= 0) ? b5n[L.b5_off(b_oup, a, b, d)] : 0.0,
                                            (b_cup >= 0) ? b5n[L.b5_off(b_cup, a, b, d)] : 0.0,
                                            b5n[L.b5_off(bid, a_up, b, d)],
                                            b5n[L.b5_off(bid, a, b_up, d)],
                                            b5n[L.b5_off(bid, a, b, d_up)]
                                        };
                                        double bd_o = (b_obd >= 0) ? b5n[L.b5_off(b_obd, a, b, d)] : 0.0;
                                        double bd_c = (b_cbd >= 0) ? b5n[L.b5_off(b_cbd, a, b, d)] : 0.0;
                                        double bd_a = (ba_hi > 0) ? b5n[L.b5_off(bid, a_bud, b, d)] : 0.0;
                                        double bd_b = (bb_hi > 0) ? b5n[L.b5_off(bid, a, b_bud, d)] : 0.0;
                                        double bd_d = (bd_hi > 0) ? b5n[L.b5_off(bid, a, b, d_bud)] : 0.0;

                                        double num0 = src, den0 = lm;
                                        { double r = Lo + (a_in ? inv_a : 0.0) + (b_in ? inv_b : 0.0) + ic; num0 += r * up[0]; den0 += r; if (bo > 0) { num0 += bo * bd_o; den0 += bo; } }
                                        { double r = Lc + io + (d_in ? inv_d : 0.0); num0 += r * up[1]; den0 += r; if (bc > 0) { num0 += bc * bd_c; den0 += bc; } }
                                        if (ak == 0) { double r = La + io; num0 += r * up[2]; den0 += r; if (ba_lo > 0) { num0 += ba_lo * bd_a; den0 += ba_lo; } }
                                        else if (ak == 3) { num0 += io * up[2]; den0 += io; }
                                        if (bk == 0) { double r = Lb + io; num0 += r * up[3]; den0 += r; if (bb_lo > 0) { num0 += bb_lo * bd_b; den0 += bb_lo; } }
                                        else if (bk == 3) { num0 += io * up[3]; den0 += io; }
                                        if (dk == 0) { double r = Ld + ic; num0 += r * up[4]; den0 += r; if (bd_lo > 0) { num0 += bd_lo * bd_d; den0 += bd_lo; } }
                                        else if (dk == 3) { num0 += ic * up[4]; den0 += ic; }

                                        double r_ext[5] = {
                                            0.5 * (co - 3), 0.5 * (cc - 2),
                                            (ak == 0) ? 0.5 * (ca - 1) : r_ovf_end(L.G5, a),
                                            (bk == 0) ? 0.5 * (cb - 1) : r_ovf_end(L.G5, b),
                                            (dk == 0) ? 0.5 * (cd - 1) : r_ovf_end(L.G5, d)
                                        };
                                        double d0_ao = (ak == 1) ? (((b_oup >= 0) ? b5n[L.b5_off(b_oup, L.G5.INF, b, d)] : 0.0) - b5n[L.b5_off(bid, L.G5.INF, b, d)]) : 0.0;
                                        double d0_bo = (bk == 1) ? (((b_oup >= 0) ? b5n[L.b5_off(b_oup, a, L.G5.INF, d)] : 0.0) - b5n[L.b5_off(bid, a, L.G5.INF, d)]) : 0.0;
                                        double d0_dc = (dk == 1) ? (((b_cup >= 0) ? b5n[L.b5_off(b_cup, a, b, L.G5.INF)] : 0.0) - b5n[L.b5_off(bid, a, b, L.G5.INF)]) : 0.0;

                                        double val = num0 / den0; int mp = -1;
                                        for (int ns = 0; ns < 12; ++ns) {
                                            double num = num0, den = den0; int mask = 0;
                                            for (int k = 0; k < 5; ++k)
                                                if (r_ext[k] > 0 && up[k] > val) { num += r_ext[k] * up[k]; den += r_ext[k]; mask |= (1 << k); }
                                            if (ak == 1) {
                                                double gr = (up[2] > val) ? (mask |= 32, 1.0) : ga_lo;
                                                num += gr * up[2]; den += gr;
                                                if (ba_hi > 0) { double br = (bd_a > val) ? (mask |= 64, ba_hi) : ba_lo; num += br * bd_a; den += br; }
                                                double d1 = up[0] - val;
                                                if (d1 > d0_ao) { mask |= 128; num += inv_a * up[0]; den += inv_a; } else num += inv_a * d0_ao;
                                            }
                                            if (bk == 1) {
                                                double gr = (up[3] > val) ? (mask |= 256, 1.0) : gb_lo;
                                                num += gr * up[3]; den += gr;
                                                if (bb_hi > 0) { double br = (bd_b > val) ? (mask |= 512, bb_hi) : bb_lo; num += br * bd_b; den += br; }
                                                double d1 = up[0] - val;
                                                if (d1 > d0_bo) { mask |= 1024; num += inv_b * up[0]; den += inv_b; } else num += inv_b * d0_bo;
                                            }
                                            if (dk == 1) {
                                                double gr = (up[4] > val) ? (mask |= 2048, 1.0) : gd_lo;
                                                num += gr * up[4]; den += gr;
                                                if (bd_hi > 0) { double br = (bd_d > val) ? (mask |= 4096, bd_hi) : bd_lo; num += br * bd_d; den += br; }
                                                double d1 = up[1] - val;
                                                if (d1 > d0_dc) { mask |= 8192; num += inv_d * up[1]; den += inv_d; } else num += inv_d * d0_dc;
                                            }
                                            double vn = num / den; if (mask == mp) { val = vn; break; } mp = mask; val = vn;
                                        }
                                        b5n[L.b5_off(bid, a, b, d)] = val;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // === STEP 5C: S4 ===
            int NS4m = L.MS4.N;
            vector<vector<double>> mid_s4(NS4m), rad_s4(NS4m);
            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
            for (int o = 0; o < NS4m; ++o) {
                int oms3 = L.MS3.idx[L.MS4.D[o]][L.MS4.C[o]], oms3_up = L.MS3.up(oms3);
                vector<double> Tf(L.N_G4_tri);
                for (int tid = 0; tid < L.N_G4_tri; ++tid)
                    Tf[tid] = ((oms3_up >= 0) ? s3n[(size_t)oms3_up * L.N_G4_tri + tid] : 0.0) - s3n[(size_t)oms3 * L.N_G4_tri + tid];
                proj_3d_sym_fp(L.PM45, L.tri_G4_map, L.tri_G5_map, L.N_G5_tri, Tf, mid_s4[o], rad_s4[o]);
            }
            vector<double>& s4n = s4;
            fill(s4n.begin(), s4n.end(), 0.0);
            int rs4_max = (int)h.dS4 - 1;
            for (int rsum = rs4_max + 4 * maxR5; rsum >= 12; --rsum) {
                for (int ro = min(rs4_max, rsum - 8); ro >= 4; --ro) {
                    int rq = rsum - ro;
                    if (rq < 8 || rq > 4 * maxR5) continue;
                    for (int co = min(ro, (int)h.KS4); co >= 4; --co) {
                        int o = L.MS4.idx[ro][co], Lo = ro - co;
                        double io = 1.0 / ro, bo = (double)Lo / ro, r_o = r1[L.G2.idx[ro][co]];
                        int o_up = L.MS4.up(o), o_bud = L.MS4.bud(o);
                        double ro_hi = 0.5 * (co - 4);
                        for (int csum = (int)quad5_rc[rq].size() - 1; csum >= 0; --csum) {
                            const auto& bucket = quad5_rc[rq][csum];
                            if (bucket.empty()) continue;
                            #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1 && bucket.size() > 1500)
                            for (size_t qi = 0; qi < bucket.size(); ++qi) {
                                int qid = bucket[qi];
                                auto [a, b, c, d] = L.quad_G5[qid];
                                int e[4] = {a, b, c, d};
                                double src = 0.0;
                                for (int k = 0; k < 4; ++k) {
                                    int u1 = e[(k+1)%4], u2 = e[(k+2)%4], u3 = e[(k+3)%4];
                                    int tid5 = L.tri_G5_map[(u1 * N5 + u2) * N5 + u3];
                                    double idon = invD_end(L.G5, e[k]);
                                    src += (idon - r_o) * mid_s4[o][tid5] + fabs(idon - r_o) * rad_s4[o][tid5];
                                }
                                double ho_up = (o_up >= 0) ? s4n[(size_t)o_up * L.N_G5_quad + qid] : 0.0;
                                double ho_bud = (o_bud >= 0) ? s4n[(size_t)o_bud * L.N_G5_quad + qid] : 0.0;

                                double num0 = src, den0 = lm, rm_b = Lo;
                                double he_up[4], he_bud[4], r_ext[4], b_lo[4], b_hi[4], g_lo[4], d0_eo[4];
                                int ek[4];
                                for (int k = 0; k < 4; ++k) {
                                    int u = e[k], u1 = e[(k+1)%4], u2 = e[(k+2)%4], u3 = e[(k+3)%4];
                                    ek[k] = L.G5.kind[u];
                                    int u_up = L.G5.up(u), u_bud = L.G5.bud(u);
                                    b_lo[k] = b_lo_end(L.G5, u); b_hi[k] = b_hi_end(L.G5, u); g_lo[k] = g_lo_end(L.G5, u, io);
                                    he_up[k] = s4n[L.s4_off(o, u_up, u1, u2, u3)];
                                    he_bud[k] = (b_hi[k] > 0) ? s4n[L.s4_off(o, u_bud, u1, u2, u3)] : 0.0;
                                    r_ext[k] = (ek[k] == 0) ? 0.5 * (L.G5.C[u] - 1) : r_ovf_end(L.G5, u);
                                    if (ek[k] == 0) {
                                        double r = L.G5.L[u] + io; num0 += r * he_up[k]; den0 += r;
                                        if (b_lo[k] > 0) { num0 += b_lo[k] * he_bud[k]; den0 += b_lo[k]; }
                                        rm_b += invD_end(L.G5, u);
                                    } else if (ek[k] == 3) { num0 += io * he_up[k]; den0 += io; rm_b += invD_end(L.G5, u); }
                                    d0_eo[k] = (ek[k] == 1) ? (((o_up >= 0) ? s4n[L.s4_off(o_up, L.G5.INF, u1, u2, u3)] : 0.0) - s4n[L.s4_off(o, L.G5.INF, u1, u2, u3)]) : 0.0;
                                }
                                num0 += rm_b * ho_up; den0 += rm_b;
                                if (bo > 0) { num0 += bo * ho_bud; den0 += bo; }

                                double val = num0 / den0; int mp = -1;
                                for (int ns = 0; ns < 12; ++ns) {
                                    double num = num0, den = den0; int mask = 0;
                                    if (ro_hi > 0 && ho_up > val) { num += ro_hi * ho_up; den += ro_hi; mask |= 1; }
                                    for (int k = 0; k < 4; ++k) {
                                        int bit = 1 << (1 + 4 * k);
                                        if (r_ext[k] > 0 && he_up[k] > val) { num += r_ext[k] * he_up[k]; den += r_ext[k]; mask |= bit; }
                                        if (ek[k] == 1) {
                                            double gr = (he_up[k] > val) ? (mask |= (bit << 1), 1.0) : g_lo[k];
                                            num += gr * he_up[k]; den += gr;
                                            if (b_hi[k] > 0) { double br = (he_bud[k] > val) ? (mask |= (bit << 2), b_hi[k]) : b_lo[k]; num += br * he_bud[k]; den += br; }
                                            double d1 = ho_up - val, iv = invD_end(L.G5, e[k]);
                                            if (d1 > d0_eo[k]) { mask |= (bit << 3); num += iv * ho_up; den += iv; } else num += iv * d0_eo[k];
                                        }
                                    }
                                    double vn = num / den; if (mask == mp) { val = vn; break; } mp = mask; val = vn;
                                }
                                s4n[(size_t)o * L.N_G5_quad + qid] = val;
                            }
                        }
                    }
                }
            }
        }
        update_analytic_f_tail(theta, f[L.f_id(2, 1)]);
        return fJ_res;
    }
};

static double run_emp_marginal(const Q5Layout& L, int64_t obs_id, int runs, int nmax,
                               vector<double>& emp_r1, vector<double>& emp_edge) {
    int NG = L.G2.N;
    emp_r1.assign(L.G2.Nin_ex, 0.28);
    emp_edge.assign((size_t)L.G2.Nin_ex * NG, 0.0);
    vector<double> cnt(L.G2.Nin_ex, 0.0), rsum(L.G2.Nin_ex, 0.0);
    double tot_core = 0.0, tot_obs = 0.0;
    mt19937_64 rng(20261124ULL);
    vector<int> deg(nmax); vector<vector<int>> adj(nmax);
    for (int run = 0; run < runs; ++run) {
        for (int i = 0; i < nmax; ++i) adj[i].clear();
        adj[0] = {1, 2}; adj[1] = {0}; adj[2] = {0}; deg[0] = 2; deg[1] = 1; deg[2] = 1;
        for (int n = 3; n < nmax; ++n) {
            int v = (int)(rng() % n), w = adj[v][rng() % adj[v].size()];
            adj[w].push_back(n); adj[n].push_back(w); deg[w]++; deg[n] = 1;
        }
        vector<int> st(nmax, -1);
        for (int u = 0; u < nmax; ++u) {
            if (deg[u] < 2) continue;
            int c = 0; for (int v : adj[u]) if (deg[v] >= 2) c++;
            tot_core += 1.0;
            if (deg[u] < L.h.DF) tot_obs += (double)m5_eval_obs(obs_id, deg[u], c);
            if (c == 0) continue;
            if (deg[u] < L.h.d2) st[u] = (c <= L.h.K2) ? L.G2.idx[deg[u]][c] : L.G2.ovf_idx[deg[u]];
            else st[u] = L.G2.t_id(c);
        }
        for (int u = 0; u < nmax; ++u) {
            if (st[u] < 0 || st[u] >= L.G2.Nin_ex) continue;
            for (int v : adj[u]) {
                if (st[v] < 0) continue;
                emp_edge[(size_t)st[u] * NG + st[v]] += 1.0;
                rsum[st[u]] += 1.0 / deg[v];
                cnt[st[u]] += 1.0;
            }
        }
    }
    for (int u = 0; u < L.G2.Nin_ex; ++u)
        emp_r1[u] = (cnt[u] >= 5.0) ? (rsum[u] / cnt[u]) : (0.30 / (1.0 + 0.05 * L.G2.D[u]));
    return (tot_core > 0.0) ? (tot_obs / tot_core) : 0.1;
}

static int64_t q_val(double x, int64_t sc) {
    assert(isfinite(x) && fabs(x) < 16.0);
    return llround(x * (double)sc);
}

int main(int argc, char** argv) {
    if (argc < 9) {
        cerr << "Usage: build_marginal <tiny|med|high> <obs_id 1..8> <lower|upper> <mu_lo_P> <mu_lo_Q> <mu_hi_P> <mu_hi_Q> <out.bin>\n";
        return 1;
    }
    string preset = argv[1];
    int64_t obs_id = stoll(argv[2]);
    string side = argv[3];
    int64_t obs_sign = (side == "upper") ? 1 : ((side == "lower") ? -1 : 0);
    if (obs_id < 1 || obs_id > M5_MAX_OBS || !obs_sign) { cerr << "invalid obs_id or side\n"; return 1; }
    int64_t mu_lo_P = stoll(argv[4]), mu_lo_Q = stoll(argv[5]);
    int64_t mu_hi_P = stoll(argv[6]), mu_hi_Q = stoll(argv[7]);
    string out_path = argv[8];

    MarginalHeader mh;
    memset(&mh, 0, sizeof(mh));
    mh.m_magic = M5_MAGIC;
    mh.obs_id = obs_id;
    mh.obs_sign = obs_sign;
    mh.theta_Q = 100000000;
    mh.theta_P = 10000000;
    mh.f_p3 = 0;
    mh.mu_lo_P = mu_lo_P; mh.mu_lo_Q = mu_lo_Q;
    mh.mu_hi_P = mu_hi_P; mh.mu_hi_Q = mu_hi_Q;

    QuintetHeader& hdr = mh.base;
    hdr.magic = Q5_MAGIC;
    hdr.mode = 1;
    hdr.Q = 100000000;
    hdr.P = llround(0.5 * ((double)mu_lo_P / mu_lo_Q + (double)mu_hi_P / mu_hi_Q) * (double)hdr.Q);
    hdr.S = Q5_WEIGHT_SCALE;
    hdr.S_r = Q5_CENTER_SCALE;
    hdr.alpha_lo_num = 0;
    hdr.Q_K = 10080;

    int n_outer = 9, emp_runs = 24, emp_n = 120000;
    if (preset == "tiny") {
        hdr.DF=120; hdr.d2=16; hdr.K2=6; hdr.d3m=10; hdr.K3m=5; hdr.d3e=10; hdr.K3e=4;
        hdr.d4i=8; hdr.K4i=4; hdr.SP4=12; hdr.dS3=8; hdr.KS3=4; hdr.d4e=8; hdr.K4e=4;
        hdr.d5i=6; hdr.K5i=4; hdr.SP5=11; hdr.d5o=6; hdr.K5o=4; hdr.d5c=6; hdr.K5c=4; hdr.SB5=10; hdr.dS4=6; hdr.KS4=4; hdr.d5e=6; hdr.K5e=3;
        int ks[] = {0, 630, 840, 1120, 1440, 1680, 2016, 2520, 2880, 3360, 3780, 4200, 5040};
        hdr.nknots = (int)(sizeof(ks)/sizeof(ks[0]));
        for (int i = 0; i < hdr.nknots; ++i) hdr.knot_num[i] = ks[i];
        n_outer = 9; emp_runs = 20; emp_n = 4000;
    } else if (preset == "med") {
        hdr.DF=1000; hdr.d2=72; hdr.K2=16; hdr.d3m=24; hdr.K3m=10; hdr.d3e=28; hdr.K3e=10;
        hdr.d4i=12; hdr.K4i=6; hdr.SP4=16; hdr.dS3=11; hdr.KS3=6; hdr.d4e=13; hdr.K4e=6;
        hdr.d5i=8; hdr.K5i=5; hdr.SP5=14; hdr.d5o=8; hdr.K5o=5; hdr.d5c=7; hdr.K5c=5; hdr.SB5=12; hdr.dS4=8; hdr.KS4=5; hdr.d5e=9; hdr.K5e=4;
        int ks[] = {0, 140, 280, 420, 560, 840, 1120, 1440, 1680, 2016, 2520, 2880, 3360, 3780, 4200, 5040};
        hdr.nknots = (int)(sizeof(ks)/sizeof(ks[0]));
        for (int i = 0; i < hdr.nknots; ++i) hdr.knot_num[i] = ks[i];
        n_outer = 9; emp_runs = 24; emp_n = 120000;
    } else if (preset == "high") {
        hdr.DF=1600; hdr.d2=120; hdr.K2=20; hdr.d3m=32; hdr.K3m=12; hdr.d3e=36; hdr.K3e=12;
        hdr.d4i=14; hdr.K4i=7; hdr.SP4=18; hdr.dS3=12; hdr.KS3=7; hdr.d4e=15; hdr.K4e=7;
        hdr.d5i=9; hdr.K5i=5; hdr.SP5=16; hdr.d5o=8; hdr.K5o=5; hdr.d5c=8; hdr.K5c=5; hdr.SB5=13; hdr.dS4=8; hdr.KS4=5; hdr.d5e=10; hdr.K5e=4;
        int ks[] = {0, 84, 168, 252, 360, 504, 672, 840, 1120, 1440, 1680, 2016, 2520, 2880, 3360, 3780, 4200, 5040};
        hdr.nknots = (int)(sizeof(ks)/sizeof(ks[0]));
        for (int i = 0; i < hdr.nknots; ++i) hdr.knot_num[i] = ks[i];
        n_outer = 9; emp_runs = 28; emp_n = 150000;
    } else {
        cerr << "Unknown preset: " << preset << "\n"; return 1;
    }

    int vh = m5_validate_header(mh);
    if (vh != 0) { cerr << "Invalid header code=" << vh << "\n"; return 1; }
    Q5Layout L; L.init(hdr);
    double mu_lo = (double)mu_lo_P / (double)mu_lo_Q;
    double mu_hi = (double)mu_hi_P / (double)mu_hi_Q;
    cout << "Building FT2_M5V1 [" << preset << " / " << m5_obs_name(obs_id) << " / " << side << "]: total_i64="
         << L.total_i64 << " (" << fixed << setprecision(1) << (L.total_bytes + 128) / 1e6 << " MB)\n" << flush;

    double t0 = omp_get_wtime();
    vector<double> emp_r1, emp_edge;
    double th_emp = run_emp_marginal(L, obs_id, emp_runs, emp_n, emp_r1, emp_edge);

    MarginalSolver S; S.init(&L, obs_id, obs_sign, mu_lo, mu_hi);
    S.r1 = emp_r1;
    auto clamp01 = [](double x) { return max(1e-7, min(1.0 - 1e-7, x)); };
    double th0 = clamp01(th_emp + (double)obs_sign * 0.002);
    double f0 = S.eval_theta(th0, emp_edge, n_outer, 5, 40);
    double th1 = clamp01(th0 + (double)obs_sign * f0 * 0.5);
    double f1 = S.eval_theta(th1, emp_edge, 6, 1, 40);
    double th2 = (fabs(f1 - f0) > 1e-12) ? clamp01(th1 - f1 * (th1 - th0) / (f1 - f0)) : th1;
    double f2 = S.eval_theta(th2, emp_edge, 6, 0, 40);
    double th3 = (fabs(f2 - f1) > 1e-12) ? clamp01(th2 - f2 * (th2 - th1) / (f2 - f1)) : th2;
    double f3 = S.eval_theta(th3, emp_edge, 5, 0, 40);
    double th_root = (fabs(f3 - f2) > 1e-12) ? clamp01(th3 - f3 * (th3 - th2) / (f3 - f2)) : th3;
    double f_fin = S.eval_theta(th_root, emp_edge, 6, 0, 40);
    cout << "  Solver converged in " << fixed << setprecision(2) << omp_get_wtime() - t0
         << "s: emp=" << setprecision(6) << th_emp << " theta=" << setprecision(10) << th_root
         << " (f[J]=" << scientific << setprecision(3) << f_fin << ")\n" << flush;

    mh.theta_P = (obs_sign == 1) ? (int64_t)ceil(th_root * (double)mh.theta_Q)
                                 : (int64_t)floor(th_root * (double)mh.theta_Q);
    if (mh.theta_P < 0) mh.theta_P = 0;
    if (mh.theta_P > mh.theta_Q) mh.theta_P = mh.theta_Q;
    hdr.f_p0 = q_val(S.A_tail, hdr.S);
    hdr.f_p1 = q_val(S.B_tail, hdr.S);
    hdr.f_p2 = 0;
    mh.f_p3 = 0;
    assert(m5_validate_header(mh) == 0);

    auto quant_move = [](vector<double>& src, int64_t sc) -> vector<int64_t> {
        vector<int64_t> dst(src.size());
        for (size_t i = 0; i < src.size(); ++i) dst[i] = q_val(src[i], sc);
        vector<double>().swap(src);
        return dst;
    };
    vector<int64_t> q_r1    = quant_move(S.r1, hdr.S_r);
    vector<int64_t> q_knots = quant_move(S.knots, hdr.S);
    vector<int64_t> q_f     = quant_move(S.f, hdr.S);
    vector<int64_t> q_g     = quant_move(S.g, hdr.S);
    q_g[L.g_off(L.G2.INF, L.G2.INF)] = 0;
    vector<int64_t> q_h3    = quant_move(S.h3, hdr.S);
    vector<int64_t> q_p4    = quant_move(S.p4, hdr.S);
    vector<int64_t> q_s3    = quant_move(S.s3, hdr.S);
    vector<int64_t> q_p5    = quant_move(S.p5, hdr.S);
    vector<int64_t> q_b5    = quant_move(S.b5, hdr.S);
    vector<int64_t> q_s4    = quant_move(S.s4, hdr.S);

    int nt = (L.n_p4 < 500000) ? 1 : 8;
    const int64_t* g_i  = q_g.data();
    const int64_t* h3_i = q_h3.data();
    const int64_t* p4_i = q_p4.data();
    const int64_t* s3_i = q_s3.data();
    int NG = L.G2.N, N3 = L.G3.N, N4 = L.G4.N, N5 = L.G5.N;
    size_t g4_p2 = (size_t)N4 * (N4 + 1) / 2, g5_p2 = (size_t)N5 * (N5 + 1) / 2;

    vector<int64_t> qe23_m(L.n_e23, 0),   qe23_r(L.n_e23, 0);
    vector<int64_t> qe34p_m(L.n_e34p, 0), qe34p_r(L.n_e34p, 0);
    vector<int64_t> qe34s_m(L.n_e34s, 0), qe34s_r(L.n_e34s, 0);
    vector<int64_t> qe45p_m(L.n_e45p, 0), qe45p_r(L.n_e45p, 0);
    vector<int64_t> qe45bp_m(L.n_e45bp, 0), qe45bp_r(L.n_e45bp, 0);
    vector<int64_t> qe45bs_m(L.n_e45bs, 0), qe45bs_r(L.n_e45bs, 0);
    vector<int64_t> qe45s4_m(L.n_e45s4, 0), qe45s4_r(L.n_e45s4, 0);

    #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
    for (int m = 0; m < L.M3.N; ++m) {
        int gu = L.G2.idx[L.M3.D[m]][L.M3.C[m]], gup = L.G2.up(gu);
        vector<int64_t> dg(NG);
        for (int w = 0; w < NG; ++w) dg[w] = g_i[L.g_off(gup, w)] - g_i[L.g_off(gu, w)];
        proj_1d_exact(L.PM23, dg.data(), &qe23_m[(size_t)m * N3], &qe23_r[(size_t)m * N3]);
    }
    #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
    for (int pid = 0; pid < L.NP4_ord; ++pid) {
        auto [b, cc] = L.p4_ord[pid];
        int b_e3 = L.G3.idx[L.M4v.D[b]][L.M4v.C[b]], b_e3up = L.G3.up(b_e3);
        int c_m3 = L.M3.idx[L.M4v.D[cc]][L.M4v.C[cc]];
        vector<int64_t> dh(N3);
        for (int d3 = 0; d3 < N3; ++d3) dh[d3] = h3_i[L.h3_off(c_m3, b_e3up, d3)] - h3_i[L.h3_off(c_m3, b_e3, d3)];
        proj_1d_exact(L.PM34, dh.data(), &qe34p_m[(size_t)pid * N4], &qe34p_r[(size_t)pid * N4]);
    }
    #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
    for (int o = 0; o < L.MS3.N; ++o) {
        int om3 = L.M3.idx[L.MS3.D[o]][L.MS3.C[o]], om3_up = L.M3.up(om3);
        vector<int64_t> Mf((size_t)N3 * N3), m2, r2;
        for (int u = 0; u < N3; ++u) for (int v = 0; v < N3; ++v)
            Mf[(size_t)u * N3 + v] = ((om3_up >= 0) ? h3_i[L.h3_off(om3_up, u, v)] : 0) - h3_i[L.h3_off(om3, u, v)];
        proj_2d_exact(L.PM34, Mf, m2, r2);
        for (int u = 0; u < N4; ++u) for (int v = u; v < N4; ++v) {
            size_t so = (size_t)o * g4_p2 + Q5Layout::sym_off(u, v);
            qe34s_m[so] = m2[(size_t)u * N4 + v]; qe34s_r[so] = r2[(size_t)u * N4 + v];
        }
    }
    #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
    for (int pid5 = 0; pid5 < L.NP5_ord; ++pid5) {
        auto [b, c, d] = L.p5_ord[pid5];
        int b_e4 = L.G4.idx[L.I5_2.D[b]][L.I5_2.C[b]], b_e4up = L.G4.up(b_e4);
        int c_m4 = L.M4v.idx[L.I5_2.D[c]][L.I5_2.C[c]], d_m4 = L.M4v.idx[L.I5_2.D[d]][L.I5_2.C[d]];
        vector<int64_t> dp(N4);
        for (int e4 = 0; e4 < N4; ++e4) dp[e4] = p4_i[L.p4_off(c_m4, d_m4, b_e4up, e4)] - p4_i[L.p4_off(c_m4, d_m4, b_e4, e4)];
        proj_1d_exact(L.PM45, dp.data(), &qe45p_m[(size_t)pid5 * N5], &qe45p_r[(size_t)pid5 * N5]);
    }
    #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
    for (int bid = 0; bid < L.NB5_int; ++bid) {
        auto [o, c] = L.b5_int[bid];
        int Do = L.I5_3.D[o], co = L.I5_3.C[o], Dc = L.I5_2c.D[c], cc = L.I5_2c.C[c];
        int o_m4 = L.M4v.idx[Do][co], o_m4up = L.M4v.up(o_m4), c_m4 = L.M4v.idx[Dc][cc];
        bool has_up = (o_m4up >= 0 && Do + 1 + Dc <= hdr.SP4);
        vector<int64_t> Mf((size_t)N4 * N4), m2, r2;
        for (int u = 0; u < N4; ++u) for (int v = 0; v < N4; ++v)
            Mf[(size_t)u * N4 + v] = (has_up ? p4_i[L.p4_off(o_m4up, c_m4, u, v)] : 0) - p4_i[L.p4_off(o_m4, c_m4, u, v)];
        proj_2d_exact(L.PM45, Mf, m2, r2);
        memcpy(&qe45bp_m[(size_t)bid * N5 * N5], m2.data(), (size_t)N5 * N5 * sizeof(int64_t));
        memcpy(&qe45bp_r[(size_t)bid * N5 * N5], r2.data(), (size_t)N5 * N5 * sizeof(int64_t));

        int o_ms3 = L.MS3.idx[Do][co], c_e4 = L.G4.idx[Dc][cc], c_e4up = L.G4.up(c_e4);
        for (int u = 0; u < N4; ++u) for (int v = 0; v < N4; ++v)
            Mf[(size_t)u * N4 + v] = s3_i[L.s3_off(o_ms3, u, v, c_e4up)] - s3_i[L.s3_off(o_ms3, u, v, c_e4)];
        proj_2d_exact(L.PM45, Mf, m2, r2);
        for (int u = 0; u < N5; ++u) for (int v = u; v < N5; ++v) {
            size_t so = (size_t)bid * g5_p2 + Q5Layout::sym_off(u, v);
            qe45bs_m[so] = m2[(size_t)u * N5 + v]; qe45bs_r[so] = r2[(size_t)u * N5 + v];
        }
    }
    #pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
    for (int o = 0; o < L.MS4.N; ++o) {
        int oms3 = L.MS3.idx[L.MS4.D[o]][L.MS4.C[o]], oms3_up = L.MS3.up(oms3);
        vector<int64_t> Tf(L.N_G4_tri);
        for (int tid = 0; tid < L.N_G4_tri; ++tid)
            Tf[tid] = ((oms3_up >= 0) ? s3_i[(size_t)oms3_up * L.N_G4_tri + tid] : 0) - s3_i[(size_t)oms3 * L.N_G4_tri + tid];
        proj_3d_sym_exact(L.PM45, L.tri_G4_map, L.tri_G5_map, L.N_G5_tri, Tf,
                          &qe45s4_m[(size_t)o * L.N_G5_tri], &qe45s4_r[(size_t)o * L.N_G5_tri]);
    }

    FILE* fp = fopen(out_path.c_str(), "wb");
    if (!fp) { perror("fopen"); return 1; }
    if (fwrite(&mh, sizeof(mh), 1, fp) != 1) { perror("fwrite header"); fclose(fp); return 1; }
    auto write_vec = [&](const vector<int64_t>& v, size_t expected) {
        assert(v.size() == expected);
        if (!v.empty() && fwrite(v.data(), sizeof(int64_t), v.size(), fp) != v.size()) {
            perror("fwrite payload"); fclose(fp); std::exit(1);
        }
    };
    write_vec(q_r1, L.n_r1);       write_vec(q_knots, L.n_knots);
    write_vec(q_f, L.n_f);         write_vec(q_g, L.n_g);
    write_vec(q_h3, L.n_h3);       write_vec(q_p4, L.n_p4);
    write_vec(q_s3, L.n_s3);       write_vec(q_p5, L.n_p5);
    write_vec(q_b5, L.n_b5);       write_vec(q_s4, L.n_s4);
    write_vec(qe23_m, L.n_e23);     write_vec(qe23_r, L.n_e23);
    write_vec(qe34p_m, L.n_e34p);   write_vec(qe34p_r, L.n_e34p);
    write_vec(qe34s_m, L.n_e34s);   write_vec(qe34s_r, L.n_e34s);
    write_vec(qe45p_m, L.n_e45p);   write_vec(qe45p_r, L.n_e45p);
    write_vec(qe45bp_m, L.n_e45bp); write_vec(qe45bp_r, L.n_e45bp);
    write_vec(qe45bs_m, L.n_e45bs); write_vec(qe45bs_r, L.n_e45bs);
    write_vec(qe45s4_m, L.n_e45s4); write_vec(qe45s4_r, L.n_e45s4);
    fclose(fp);
    cout << "  Wrote " << out_path << " (" << (L.total_bytes + 128) << " bytes), theta="
         << mh.theta_P << "/" << mh.theta_Q << ", tail=(" << hdr.f_p0 << ", " << hdr.f_p1 << ")\n" << flush;
    return 0;
}
