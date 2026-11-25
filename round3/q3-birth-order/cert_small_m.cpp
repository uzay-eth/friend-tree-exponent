/*
 * cert_small_m.cpp — 100% Exact Outward-Rounded Integer (__int128) Certificate
 * Proves pre-placement conditional stochastic dominance:
 *   P(Z_{n+1} > x | T_n) > P(Z_{n+2} > x | T_n)
 * for EVERY finite tree T_n of size n in {2, 3, ..., 197} (i.e. b = n+1 in {3, ..., 198},
 * m = n+2 in {4, ..., 199}) and EVERY threshold x in [0, 1)!
 *
 * Together with cert_all_m.cpp (all m >= 200, all x in [0, 1)) and the elementary
 * founder seam b = 2 (Z_1 =^d Z_2 >_st Z_3), this completes a 100% self-contained,
 * zero-floating-point-trust proof of BOTH:
 *   (1) Universal pre-placement dominance P(Z_{n+1} > x | T_n) > P(Z_{n+2} > x | T_n)
 *       for ALL n >= 2, ALL finite trees T_n, and ALL x in [0, 1);
 *   (2) Unconditional strict first-order stochastic dominance Z_1 =^d Z_2 >_st Z_3 >_st ...
 *       (Addario-Berry et al. 2024/2026 §7 Question 3, second half).
 */

#include <bits/stdc++.h>
using namespace std;

using i128 = __int128_t;

static constexpr int64_t S   = 1000000000000LL;     // 10^12 (grid probability scale)
static constexpr int64_t S2  = 1000000000000000LL;  // 10^15 (internal transcendental scale)
static constexpr int64_t LSC = 1000000000LL;        // 10^9  (normalized ledger scale)
static constexpr int H_DEN   = 128;                 // h = 1/128 for m0 = 200 solver

static inline i128 div_fl(i128 a, i128 b) {
    assert(b > 0);
    if (a >= 0) return a / b;
    return -((-a + b - 1) / b);
}
static inline i128 div_cl(i128 a, i128 b) {
    assert(b > 0);
    if (a >= 0) return (a + b - 1) / b;
    return -((-a) / b);
}

struct Iv { i128 lo, hi; };

static inline Iv iv_from_rat(i128 p, i128 q) {
    assert(q > 0);
    return {div_fl(p * S2, q), div_cl(p * S2, q)};
}
static inline Iv iv_add(Iv a, Iv b) { return {a.lo + b.lo, a.hi + b.hi}; }
static inline Iv iv_sub(Iv a, Iv b) { return {a.lo - b.hi, a.hi - b.lo}; }
static inline Iv iv_mul_nn(Iv a, Iv b) {
    assert(a.lo >= 0 && b.lo >= 0);
    return {div_fl(a.lo * b.lo, S2), div_cl(a.hi * b.hi, S2)};
}
static inline Iv iv_div_pos(Iv a, Iv b) {
    assert(a.lo >= 0 && b.lo > 0);
    return {div_fl(a.lo * S2, b.hi), div_cl(a.hi * S2, b.lo)};
}

// Exact outward enclosure of exp(x) for 0 <= x <= 16 (overflow-free in i128)
static Iv iv_exp_pos(Iv x) {
    assert(x.lo >= 0 && x.hi <= (i128)16 * S2);
    int k = 0;
    Iv u = x;
    while (u.hi > S2 / 8) {
        u = {div_fl(u.lo, 2), div_cl(u.hi, 2)};
        k++;
    }
    Iv term = {S2, S2}, sum = {S2, S2};
    for (int j = 1; j <= 15; j++) {
        term = {div_fl(term.lo * u.lo, (i128)S2 * j), div_cl(term.hi * u.hi, (i128)S2 * j)};
        sum = iv_add(sum, term);
    }
    sum.hi += div_cl(term.hi * u.hi * 2, (i128)S2 * 16);
    for (int i = 0; i < k; i++) sum = iv_mul_nn(sum, sum);
    return sum;
}

// Exact outward enclosure of exp(-x) for ANY x >= 0 via chunking by 4
static Iv iv_exp_neg(Iv x) {
    assert(x.lo >= 0);
    if (x.lo >= (i128)45 * S2) return {0, 1};
    static Iv em4 = {0, 0};
    if (em4.lo == 0) {
        Iv e4 = iv_exp_pos({(i128)4 * S2, (i128)4 * S2});
        em4 = {div_fl((i128)S2 * S2, e4.hi), div_cl((i128)S2 * S2, e4.lo)};
    }
    auto sc = [&](i128 val, bool is_hi) -> i128 {
        if (val >= (i128)45 * S2) return is_hi ? 1 : 0;
        int chunks = (int)(val / ((i128)4 * S2));
        i128 rem = val - (i128)chunks * 4 * S2;
        Iv er = iv_exp_pos({rem, rem});
        Iv res = {div_fl((i128)S2 * S2, er.hi), div_cl((i128)S2 * S2, er.lo)};
        for (int c = 0; c < chunks; c++) res = iv_mul_nn(res, em4);
        return is_hi ? res.hi : res.lo;
    };
    return {sc(x.hi, false), sc(x.lo, true)};
}

static Iv iv_log_ge1_raw(Iv y) {
    assert(y.lo >= S2);
    auto log_unit = [](Iv z) -> Iv {
        Iv num = {max((i128)0, z.lo - S2), z.hi - S2};
        Iv den = {z.lo + S2, z.hi + S2};
        Iv w = iv_div_pos(num, den);
        Iv w2 = iv_mul_nn(w, w);
        Iv poww = w, s = {0, 0};
        for (int j = 0; j < 20; j++) {
            int d = 2 * j + 1;
            s.lo += div_fl(poww.lo, d);
            s.hi += div_cl(poww.hi, d);
            poww = iv_mul_nn(poww, w2);
        }
        s.hi += div_cl(poww.hi * 9, 41 * 8);
        return {2 * s.lo, 2 * s.hi};
    };
    static Iv ln2 = {0, 0};
    if (ln2.lo == 0) ln2 = log_unit({2 * S2, 2 * S2});
    auto sc = [&](i128 val, bool is_hi) -> i128 {
        int k = 0;
        Iv z = {val, val};
        while (z.hi >= 2 * S2) {
            z = {div_fl(z.lo, 2), div_cl(z.hi, 2)};
            k++;
        }
        Iv res = iv_add({(i128)k * ln2.lo, (i128)k * ln2.hi}, log_unit(z));
        return is_hi ? res.hi : res.lo;
    };
    return {sc(y.lo, false), sc(y.hi, true)};
}

// Exact outward enclosure of y^{p/q} for y > 0 and 0 < p < q
static Iv iv_pow_rat(Iv y, int p, int q) {
    assert(y.lo > 0 && p > 0 && p < q);
    auto sc = [&](i128 val, bool is_hi) -> i128 {
        if (val == S2) return S2;
        if (val > S2) {
            Iv lg = iv_log_ge1_raw({val, val});
            Iv prod = {div_fl(lg.lo * p, q), div_cl(lg.hi * p, q)};
            Iv res = iv_exp_pos(prod);
            return is_hi ? res.hi : res.lo;
        } else {
            Iv inv = {div_fl((i128)S2 * S2, val), div_cl((i128)S2 * S2, val)};
            Iv lg = iv_log_ge1_raw(inv);
            Iv prod = {div_fl(lg.lo * p, q), div_cl(lg.hi * p, q)};
            Iv res = iv_exp_neg(prod);
            return is_hi ? res.hi : res.lo;
        }
    };
    return {sc(y.lo, false), sc(y.hi, true)};
}

static map<pair<int,int>, Iv> gamma_cache;
static Iv iv_gamma_rat(int p, int q) {
    assert(p > 0 && p < q);
    int g = std::gcd(p, q);
    p /= g; q /= g;
    auto key = make_pair(p, q);
    auto it = gamma_cache.find(key);
    if (it != gamma_cache.end()) return it->second;

    const int N = 60;
    Iv x = iv_from_rat((i128)N * q + p, q);
    Iv lnx = iv_log_ge1_raw(x);
    Iv base = iv_sub(iv_mul_nn(iv_sub(x, {S2 / 2, S2 / 2}), lnx), x);
    Iv ln2pi = iv_log_ge1_raw({6283185307179586LL, 6283185307179587LL});
    Iv half_ln2pi = {div_fl(ln2pi.lo, 2), div_cl(ln2pi.hi, 2)};
    i128 j_hi = div_cl((i128)S2 * S2, 12 * x.lo);
    i128 x3_lo = div_fl(div_fl(x.lo * x.lo, S2) * x.lo, S2);
    i128 j_lo = div_fl((i128)S2 * S2, 12 * x.hi) - div_cl((i128)S2 * S2, 360 * x3_lo);
    Iv lnG = iv_add(iv_add(base, half_ln2pi), {j_lo, j_hi});
    Iv sum_ln = {0, 0};
    for (int k = 1; k < N; k++) {
        sum_ln = iv_add(sum_ln, iv_log_ge1_raw(iv_from_rat((i128)k * q + p, q)));
    }
    Iv ln_qp = iv_log_ge1_raw(iv_from_rat(q, p));
    Iv ln_res = iv_add(iv_sub(lnG, sum_ln), ln_qp);
    assert(ln_res.lo >= 0);
    Iv res = iv_exp_pos(ln_res);
    gamma_cache[key] = res;
    return res;
}

static i128 exp_gamma_upper_hi_frac(int p, int q, i128 v_lo) {
    assert(p > 0 && p < q && v_lo > 0);
    Iv gam = iv_gamma_rat(p, q);
    Iv v_iv = {v_lo, v_lo};
    Iv va = iv_pow_rat(v_iv, p, q);
    if (v_lo <= 2 * S2) {
        i128 s_lo = 0;
        Iv pow_vn = {S2, S2};
        for (int n = 0; n <= 23; n++) {
            i128 den = (i128)p + (i128)n * q;
            if (n % 2 == 0) s_lo += div_fl(pow_vn.lo * q, den);
            else s_lo -= div_cl(pow_vn.hi * q, den);
            if (n < 23) {
                pow_vn = {
                    div_fl(pow_vn.lo * v_lo, (i128)S2 * (n + 1)),
                    div_cl(pow_vn.hi * v_lo, (i128)S2 * (n + 1))
                };
            }
        }
        if (s_lo < 0) s_lo = 0;
        i128 gam_inc_lo = div_fl(va.lo * s_lo, S2);
        i128 Q_hi = S2 - div_fl(gam_inc_lo * S2, gam.hi);
        if (Q_hi < 0) Q_hi = 0;
        Iv ev = iv_exp_pos(v_iv);
        return div_cl(ev.hi * Q_hi, S2);
    } else {
        const int K = 45;
        i128 D_odd_lo = v_lo;
        for (int n = K; n >= 1; n--) {
            i128 D_even_hi = S2 + div_cl((i128)n * S2 * S2, D_odd_lo);
            i128 a_even_num = (i128)n * q - p;
            D_odd_lo = v_lo + div_fl(a_even_num * S2 * S2, (i128)q * D_even_hi);
        }
        i128 ev_Gam_hi = div_cl(va.hi * S2, D_odd_lo);
        return div_cl(ev_Gam_hi * S2, gam.lo);
    }
}

static i128 exp_gamma_upper_hi_general(int int_part, int p, int q, i128 v_lo, i128 v_hi) {
    if (int_part == 0) return exp_gamma_upper_hi_frac(p, q, v_lo);
    if (p == 0) {
        if (int_part == 1) return S2;
        if (int_part == 2) return S2 + v_hi;
        assert(false);
    }
    i128 base = exp_gamma_upper_hi_frac(p, q, v_hi);
    Iv gam_a = iv_gamma_rat(p, q);
    i128 gam1_lo = div_fl(gam_a.lo * p, q);
    Iv va = iv_pow_rat({v_hi, v_hi}, p, q);
    i128 term1 = div_cl(va.hi * S2, gam1_lo);
    if (int_part == 1) return base + term1;
    i128 gam2_lo = div_fl(gam1_lo * (p + q), q);
    i128 v1a_hi = div_cl(va.hi * v_hi, S2);
    i128 term2 = div_cl(v1a_hi * S2, gam2_lo);
    return base + term1 + term2;
}

static vector<Iv> EXP_GRID;
static void init_exp_grid(int max_i) {
    EXP_GRID.resize(max_i + 1);
    for (int i = 0; i <= max_i; i++) EXP_GRID[i] = iv_exp_pos(iv_from_rat(i, H_DEN));
}

// --- Exact fixed-integer-N primitives (for N >= 4, ANY y = y_num/y_den in [0, 2]) ---

// Exact outward LOWER bound (scale S) on H_{N, j}(x) where x = 1 - exp(-y_hi), y_hi = yn/yd in [0, 2]!
// Uses the exact identity B_{N, r}(x) = e^{-(N-1)y} * sum_{i=0}^{r-1} C(N-1, i) (e^y - 1)^i.
static int64_t H_exact_N_lower(int N, int j, int64_t yn, int64_t yd) {
    if (j <= 0) return 0;
    if (j >= N) return S; // impossible in tree of size N, but B_{N, N} = 1
    static i128 inv_e_lo = 0;
    if (inv_e_lo == 0) inv_e_lo = iv_exp_neg({S2, S2}).lo;
    if (yn == 0) {
        i128 sum = S2, fk = S2;
        for (int k = 1; k < j; k++) { fk = div_fl(fk, k); sum += fk; }
        i128 ans = div_fl(inv_e_lo * sum, S2);
        return (int64_t)min((i128)S, div_fl(ans * S, S2));
    }
    // L = (N - 1) * y_hi
    Iv y_iv = iv_from_rat(yn, yd);
    i128 L_hi = div_cl((i128)(N - 1) * yn * S2, yd);
    if (L_hi >= (i128)45 * S2) return 0;
    // z = e^{y_hi} - 1 (need LOWER bound z_lo!)
    Iv ey = iv_exp_pos(y_iv);
    i128 z_lo = max((i128)0, ey.lo - S2);

    int exp_chunks = (int)(L_hi / ((i128)4 * S2));
    i128 exp_rem = L_hi - (i128)exp_chunks * 4 * S2;
    static i128 e4_hi = 0;
    if (e4_hi == 0) e4_hi = iv_exp_pos({(i128)4 * S2, (i128)4 * S2}).hi;
    i128 erem_hi = iv_exp_pos({exp_rem, exp_rem}).hi;

    int shift = 0;
    i128 cur = S2, sum_val = S2;
    i128 B_val[36];
    int B_sh[36];
    B_val[1] = sum_val; B_sh[1] = 0;
    for (int i = 1; i < j; i++) {
        // multiply cur by (N - i) * z_lo / (i * S2)
        while (cur > ((i128)1 << 60)) {
            cur >>= 1;
            sum_val >>= 1;
            shift++;
        }
        i128 fac_lo = div_fl((i128)(N - i) * z_lo, i);
        cur = div_fl(cur * fac_lo, S2);
        sum_val += cur;
        while (sum_val > ((i128)1 << 60)) {
            cur >>= 1;
            sum_val >>= 1;
            shift++;
        }
        B_val[i + 1] = sum_val;
        B_sh[i + 1] = shift;
    }
    i128 B_lo[36];
    for (int r = 1; r <= j; r++) {
        i128 v = B_val[r];
        int sh = B_sh[r];
        v = div_fl(v * S2, erem_hi);
        for (int c = 0; c < exp_chunks; c++) {
            int take = min(sh, 6);
            v <<= take; sh -= take;
            v = div_fl(v * S2, e4_hi);
        }
        while (sh > 0 && v < S2) { v <<= 1; sh--; }
        if (sh > 0 || v > S2) v = S2;
        B_lo[r] = v;
    }
    i128 ans = 0, fk = S2;
    for (int k = 0; k < j; k++) {
        if (k > 0) fk = div_fl(fk, k);
        ans += div_fl(fk * B_lo[j - k], S2);
    }
    ans = div_fl(inv_e_lo * ans, S2);
    return (int64_t)min((i128)S, div_fl(ans * S, S2));
}

// Precomputed log(k) for k = 1..2050 for exact small-N Beta prefactor!
static vector<Iv> LOG_INT;
static void init_log_int(int max_k) {
    LOG_INT.resize(max_k + 1);
    LOG_INT[1] = {0, 0};
    for (int k = 2; k <= max_k; k++) {
        LOG_INT[k] = iv_log_ge1_raw({(i128)k * S2, (i128)k * S2});
    }
}

// Precomputed quadrature nodes v_0 < v_1 < ... < v_K for exact small-N Beta integral
static vector<Iv> Q_V, Q_EMV;
static void init_quad_nodes() {
    // v_0 = 1/64, geometric up to 1, then step 1/4 up to 8, then step 1/2 up to 28
    vector<pair<int,int>> rats;
    rats.push_back({1, 64});
    rats.push_back({1, 32});
    rats.push_back({1, 16});
    rats.push_back({1, 8});
    rats.push_back({1, 4});
    rats.push_back({1, 2});
    rats.push_back({3, 4});
    rats.push_back({1, 1});
    for (int k = 5; k <= 32; k++) rats.push_back({k, 4});    // 1.25 .. 8.0
    for (int k = 17; k <= 56; k++) rats.push_back({k, 2});   // 8.5 .. 28.0
    Q_V.resize(rats.size());
    Q_EMV.resize(rats.size());
    for (size_t i = 0; i < rats.size(); i++) {
        Q_V[i] = iv_from_rat(rats[i].first, rats[i].second);
        Q_EMV[i] = iv_exp_neg(Q_V[i]);
    }
}

// Exact prefactor UPPER bound on Gamma(N) / (Gamma(r) * Gamma(N - r) * B^r) where B = N - r,
// for integer N >= 4 and r = int_part + p/q (0 < p < q, N > r)!
// Wait: we actually return UPPER bound on Gamma(N) / (Gamma(N - r) * B^r) (without 1/Gamma(r)),
// OR with 1/Gamma(r) included! Let's return UPPER bound on C_{N,r} = Gamma(N) / (Gamma(r) * Gamma(N - r) * B^r)!
static map<tuple<int,int,int,int>, i128> beta_pref_cache;
static i128 beta_pref_hi(int N, int int_part, int p, int q) {
    assert(p > 0 && p < q);
    int g = std::gcd(p, q);
    p /= g; q /= g;
    auto key = make_tuple(N, int_part, p, q);
    auto it = beta_pref_cache.find(key);
    if (it != beta_pref_cache.end()) return it->second;

    // B = N - int_part - p/q = ((N - int_part)*q - p) / q = M_int + (q - p)/q where M_int = N - int_part - 1 >= 1!
    int M_int = N - int_part - 1;
    assert(M_int >= 1);
    // log Gamma(N) = sum_{k=1}^{N-1} log(k)
    // log Gamma(B) = log Gamma((q-p)/q) + sum_{j=0}^{M_int - 1} log(j + (q-p)/q)
    // Wait: combine sum_{j=1}^{M_int-1} [log(j) - log(j + (q-p)/q)] to avoid cancellation? No, direct sum in Iv has error ~ N*10^{-14} << 10^{-10}!
    // log Gamma(r) = log Gamma(p/q) + sum_{j=0}^{int_part-1} log(j + p/q)
    Iv ln_num = {0, 0};
    for (int k = 2; k < N; k++) ln_num = iv_add(ln_num, LOG_INT[k]);

    Iv gam_a = iv_gamma_rat(p, q);
    Iv gam_1ma = iv_gamma_rat(q - p, q);
    Iv ln_den = iv_add(iv_log_ge1_raw(gam_a), iv_log_ge1_raw(gam_1ma));
    // Add log(j + p/q) for j = 0..int_part-1 (note j=0 is p/q < 1 -> subtract log(q/p)!)
    if (int_part >= 1) {
        ln_den = iv_sub(ln_den, iv_log_ge1_raw(iv_from_rat(q, p)));
    }
    if (int_part >= 2) {
        ln_den = iv_add(ln_den, iv_log_ge1_raw(iv_from_rat(q + p, q)));
    }
    // Add log(j + (q-p)/q) for j = 0..M_int-1 (j=0 is (q-p)/q < 1 -> subtract log(q/(q-p))!)
    ln_den = iv_sub(ln_den, iv_log_ge1_raw(iv_from_rat(q, q - p)));
    for (int j = 1; j < M_int; j++) {
        ln_den = iv_add(ln_den, iv_log_ge1_raw(iv_from_rat((i128)j * q + (q - p), q)));
    }
    // Add r * log(B) where B = ((N - int_part)*q - p) / q > 1
    Iv B_iv = iv_from_rat((i128)(N - int_part) * q - p, q);
    Iv lnB = iv_log_ge1_raw(B_iv);
    i128 r_num = (i128)int_part * q + p;
    Iv r_lnB = {div_fl(lnB.lo * r_num, q), div_cl(lnB.hi * r_num, q)};
    ln_den = iv_add(ln_den, r_lnB);

    // ln_C = ln_num - ln_den. Could be negative (since 1/Gamma(r) <= 1 and Gamma(B+r)/(Gamma(B)B^r) ~ 1)!
    Iv ln_C = iv_sub(ln_num, ln_den);
    i128 res_hi;
    if (ln_C.hi >= 0) {
        res_hi = iv_exp_pos({max((i128)0, ln_C.lo), ln_C.hi}).hi;
    } else {
        res_hi = iv_exp_neg({-ln_C.hi, -ln_C.lo}).hi;
    }
    beta_pref_cache[key] = res_hi;
    return res_hi;
}

// Exact outward UPPER bound (scale S) on B_{N, r}(x) for fixed integer N >= 4,
// r = int_part + p/q (int_part in {0,1,2}, 0 <= p < q, N > r), and y_lo = yn/yd in (0, 2]!
static int64_t B_exact_N_upper(int N, int int_part, int p, int q, int64_t yn, int64_t yd) {
    if (yn <= 0) return S;
    assert(yn <= 2 * yd); // y <= 2
    // Exact L = (N - 1) * y_lo
    i128 L_lo_fl = div_fl((i128)(N - 1) * yn * S2, yd);
    if (L_lo_fl >= (i128)45 * S2) return 1; // <= exp(-45) * 2 * (1 + 150 + 150^2/2) < 10^{-15} < 1/S
    i128 exp_q_hi = iv_exp_neg({L_lo_fl, L_lo_fl}).hi;

    // Special case r == 1 (exact!)
    if (int_part == 1 && p == 0) {
        return (int64_t)min((i128)S, div_cl(exp_q_hi * S, S2));
    }
    assert(p > 0 && p < q);
    int g = std::gcd(p, q);
    p /= g; q /= g;

    // Exact B = N - r = ((N - int_part)*q - p) / q
    i128 B_num = (i128)(N - int_part) * q - p;
    assert(B_num > q); // B > 1
    Iv B_iv = iv_from_rat(B_num, q);
    Iv y_iv = iv_from_rat(yn, yd);
    Iv ey = iv_exp_pos(y_iv);
    Iv z_iv = {max((i128)1, ey.lo - S2), ey.hi - S2}; // e^{y_lo} - 1
    Iv w_iv = iv_mul_nn(B_iv, z_iv);

    // Fast path for N >= 60 (B >= 57): Stieltjes CF + Gautschi/exact prefactor!
    if (N >= 60) {
        if (int_part == 0) {
            int Vcut = 24;
            i128 alpha_lo = S2 - div_cl((i128)Vcut * S2 * S2, 2 * B_iv.lo);
            i128 w_eff_lo = div_fl(w_iv.lo * S2, alpha_lo);
            i128 inv_alpha_hi = div_cl((i128)S2 * S2, alpha_lo);
            i128 alpha_pow_hi = iv_pow_rat({inv_alpha_hi, inv_alpha_hi}, q - p, q).hi;
            i128 main_part = div_cl(alpha_pow_hi * exp_gamma_upper_hi_frac(p, q, w_eff_lo), S2);

            i128 base_tail_lo = w_iv.lo + alpha_lo * Vcut;
            i128 inv_bt_hi = div_cl((i128)S2 * S2, base_tail_lo);
            i128 bt_pow_hi = iv_pow_rat({inv_bt_hi, inv_bt_hi}, q - p, q).hi;
            i128 evcut_hi = iv_exp_neg({(i128)Vcut * S2, (i128)Vcut * S2}).hi;
            Iv gam = iv_gamma_rat(p, q);
            i128 tail_part = div_cl(div_cl(bt_pow_hi * evcut_hi, S2) * S2, gam.lo);

            i128 ans = div_cl(exp_q_hi * (main_part + tail_part), S2);
            return (int64_t)min((i128)S, div_cl(ans * S, S2));
        } else {
            i128 pref_hi;
            if (int_part == 1) {
                pref_hi = S2 + div_cl((i128)p * S2 * S2, (i128)q * B_iv.lo);
            } else {
                i128 f1 = S2 + div_cl((i128)(p + q) * S2 * S2, (i128)q * B_iv.lo);
                i128 f2 = S2 + div_cl((i128)p * S2 * S2, (i128)q * B_iv.lo);
                pref_hi = div_cl(f1 * f2, S2);
            }
            i128 val = div_cl(pref_hi * exp_gamma_upper_hi_general(int_part, p, q, w_iv.lo, w_iv.hi), S2);
            i128 ans = div_cl(exp_q_hi * val, S2);
            return (int64_t)min((i128)S, div_cl(ans * S, S2));
        }
    }

    // Exact small-N (4 <= N < 60) monotone quadrature with exact Gamma prefactor!
    // Integral I = int_0^inf (w + B*(1 - exp(-v/B)))^{r-1} e^{-v} dv
    i128 C_hi = beta_pref_hi(N, int_part, p, q); // upper bound on Gamma(N) / (Gamma(r)*Gamma(N-r)*B^r)
    int K = (int)Q_V.size();
    i128 I_hi = 0;

    if (int_part == 0) {
        // r = p/q in (0, 1): g(v) = (w + B*(1 - exp(-v/B)))^{r-1} is DECREASING in v!
        // 1. [0, v_0] with v_0 = 1/64: alpha_0 = 1 - v_0 / (2*B_lo)
        i128 v0_hi = Q_V[0].hi;
        i128 alpha0_lo = S2 - div_cl(v0_hi * S2, 2 * B_iv.lo);
        // Integral <= ((w_lo + alpha0_lo * v0_hi)^r - w_lo^r) / (r * alpha0_lo)
        // Wait: is ((w + A)^r - w^r) decreasing in w for 0 < r < 1? YES! d/dw = r*((w+A)^{r-1} - w^{r-1}) < 0!
        // So evaluate at w = w_iv.lo (lower bound on w -> UPPER bound on ((w+A)^r - w^r)!)
        i128 A_hi = div_cl(alpha0_lo * v0_hi, S2); // Wait: need UPPER bound on int_0^{v_0} (w_lo + alpha0_lo * v)^{r-1} dv
        // Exact antiderivative of (w_lo + alpha0_lo * v)^{r-1} on [0, v0_hi] is ((w_lo + alpha0_lo * v0_hi)^r - w_lo^r) / (r * alpha0_lo)!
        i128 p_hi = iv_pow_rat({w_iv.lo + A_hi, w_iv.lo + A_hi}, p, q).hi;
        i128 p_lo = iv_pow_rat({w_iv.lo, w_iv.lo}, p, q).lo;
        i128 diff = max((i128)0, p_hi - p_lo);
        // Divide by (p/q) * alpha0_lo
        I_hi += div_cl(diff * q * S2, (i128)p * alpha0_lo);

        // 2. [v_i, v_{i+1}] for i = 0..K-2: g(v_i) * (e^{-v_i} - e^{-v_{i+1}})
        for (int i = 0; i < K - 1; i++) {
            // Lower bound on B*(1 - exp(-v_i/B)):
            // exp(-v_i/B)_hi using v_i.lo / B.hi (lower bound on v_i/B -> upper bound on exp(-v_i/B))!
            i128 arg_lo = div_fl(Q_V[i].lo * S2, B_iv.hi);
            i128 em_hi = iv_exp_neg({arg_lo, arg_lo}).hi;
            i128 one_m_lo = max((i128)0, S2 - em_hi);
            i128 base_lo = w_iv.lo + div_fl(B_iv.lo * one_m_lo, S2);
            i128 inv_b_hi = div_cl((i128)S2 * S2, base_lo);
            i128 g_hi = iv_pow_rat({inv_b_hi, inv_b_hi}, q - p, q).hi;
            i128 de_hi = max((i128)0, Q_EMV[i].hi - Q_EMV[i + 1].lo);
            I_hi += div_cl(g_hi * de_hi, S2);
        }
        // 3. Tail [v_{K-1}, inf): g(v_{K-1}) * e^{-v_{K-1}}
        {
            i128 arg_lo = div_fl(Q_V[K - 1].lo * S2, B_iv.hi);
            i128 em_hi = iv_exp_neg({arg_lo, arg_lo}).hi;
            i128 base_lo = w_iv.lo + div_fl(B_iv.lo * max((i128)0, S2 - em_hi), S2);
            i128 inv_b_hi = div_cl((i128)S2 * S2, base_lo);
            i128 g_hi = iv_pow_rat({inv_b_hi, inv_b_hi}, q - p, q).hi;
            I_hi += div_cl(g_hi * Q_EMV[K - 1].hi, S2);
        }
    } else {
        // r = int_part + p/q > 1: g(v) = (w + B*(1 - exp(-v/B)))^{r-1} is INCREASING in v!
        // Helper for upper bound on g(v) at v_hi:
        auto g_hi_at = [&](i128 v_hi) -> i128 {
            // Upper bound on B*(1 - exp(-v/B)): use v_hi / B_lo -> lower bound on exp(-v/B) -> upper bound on 1 - exp!
            i128 arg_hi = div_cl(v_hi * S2, B_iv.lo);
            i128 em_lo = iv_exp_neg({arg_hi, arg_hi}).lo;
            i128 base_hi = w_iv.hi + div_cl(B_iv.hi * (S2 - em_lo), S2);
            // base_hi^{r-1} = base_hi^{int_part - 1} * base_hi^{p/q}
            i128 frac_hi = iv_pow_rat({base_hi, base_hi}, p, q).hi;
            if (int_part == 1) return frac_hi;
            return div_cl(base_hi * frac_hi, S2);
        };
        // [0, v_0]: g(v_0) * (1 - e^{-v_0})
        I_hi += div_cl(g_hi_at(Q_V[0].hi) * (S2 - Q_EMV[0].lo), S2);
        // [v_i, v_{i+1}] for i = 0..K-2: g(v_{i+1}) * (e^{-v_i} - e^{-v_{i+1}})
        for (int i = 0; i < K - 1; i++) {
            i128 de_hi = max((i128)0, Q_EMV[i].hi - Q_EMV[i + 1].lo);
            I_hi += div_cl(g_hi_at(Q_V[i + 1].hi) * de_hi, S2);
        }
        // Tail [v_{K-1}, inf): (w + B)^{r-1} * e^{-v_{K-1}}
        {
            i128 base_hi = w_iv.hi + B_iv.hi;
            i128 frac_hi = iv_pow_rat({base_hi, base_hi}, p, q).hi;
            i128 g_max = (int_part == 1) ? frac_hi : div_cl(base_hi * frac_hi, S2);
            I_hi += div_cl(g_max * Q_EMV[K - 1].hi, S2);
        }
    }

    i128 val = div_cl(C_hi * I_hi, S2);
    i128 ans = div_cl(exp_q_hi * val, S2);
    return (int64_t)min((i128)S, div_cl(ans * S, S2));
}

// --- m0 = 200 piecewise-linear initializer (identical to cert_all_m.cpp, returning full state arrays at N=200!) ---

static int64_t H_raw_lower_200(int j, int t_idx, int64_t tau_num, int64_t tau_den, int m0) {
    if (j <= 0) return 0;
    static i128 inv_e_lo = 0;
    if (inv_e_lo == 0) inv_e_lo = iv_exp_neg({S2, S2}).lo;
    if (tau_num == 0) {
        i128 sum = S2, fk = S2;
        for (int k = 1; k < j; k++) { fk = div_fl(fk, k); sum += fk; }
        i128 ans = div_fl(inv_e_lo * sum, S2);
        return (int64_t)min((i128)S, div_fl(ans * S, S2));
    }
    Iv et = EXP_GRID[t_idx];
    i128 a_plus_lo = div_fl((i128)m0 * et.lo - S2, m0 - 1);
    i128 a_plus_hi = div_cl((i128)m0 * et.hi - S2, m0 - 1);
    i128 a_plus_shift_hi = a_plus_hi + div_cl(S2, m0 - 1);
    i128 L_exp_hi = div_cl(a_plus_shift_hi * tau_num, tau_den);
    if (L_exp_hi >= (i128)42 * S2) return 0;
    int exp_chunks = (int)(L_exp_hi / ((i128)4 * S2));
    i128 exp_rem = L_exp_hi - (i128)exp_chunks * 4 * S2;
    static i128 e4_hi = 0;
    if (e4_hi == 0) e4_hi = iv_exp_pos({(i128)4 * S2, (i128)4 * S2}).hi;
    i128 erem_hi = iv_exp_pos({exp_rem, exp_rem}).hi;

    int shift = 0;
    i128 cur = S2, sum_val = S2;
    i128 B_val[32]; int B_sh[32];
    B_val[1] = sum_val; B_sh[1] = 0;
    for (int i = 1; i < j; i++) {
        i128 f_num_lo = a_plus_lo - div_cl((i128)(i - 1) * S2, m0 - 1);
        if (f_num_lo <= 0) cur = 0;
        else {
            i128 fac_lo = div_fl(f_num_lo * tau_num, (i128)tau_den * i);
            while (cur > ((i128)1 << 62) && fac_lo > 0) { cur >>= 1; sum_val >>= 1; shift++; }
            cur = div_fl(cur * fac_lo, S2);
            sum_val += cur;
            while (sum_val > ((i128)1 << 62)) { cur >>= 1; sum_val >>= 1; shift++; }
        }
        B_val[i + 1] = sum_val; B_sh[i + 1] = shift;
    }
    i128 B_lo[32];
    for (int r = 1; r <= j; r++) {
        i128 v = B_val[r]; int sh = B_sh[r];
        v = div_fl(v * S2, erem_hi);
        for (int c = 0; c < exp_chunks; c++) {
            int take = min(sh, 6); v <<= take; sh -= take;
            v = div_fl(v * S2, e4_hi);
        }
        while (sh > 0 && v < S2) { v <<= 1; sh--; }
        if (sh > 0 || v > S2) v = S2;
        B_lo[r] = v;
    }
    i128 ans = 0, fk = S2;
    for (int k = 0; k < j; k++) {
        if (k > 0) fk = div_fl(fk, k);
        ans += div_fl(fk * B_lo[j - k], S2);
    }
    ans = div_fl(inv_e_lo * ans, S2);
    return (int64_t)min((i128)S, div_fl(ans * S, S2));
}

static int64_t B_raw_upper_200(int int_part, int p, int q, int t_idx, int64_t tau_num, int64_t tau_den, int m0) {
    if (tau_num <= 0) return S;
    assert(tau_num * 2 <= tau_den * (m0 - 1));
    Iv et = EXP_GRID[t_idx];
    i128 L_lo_fl = div_fl(et.lo * tau_num, tau_den);
    i128 L_lo_cl = div_cl(et.hi * tau_num, tau_den);
    if (L_lo_fl >= (i128)42 * S2) return 1;
    i128 exp_q_hi = iv_exp_neg({L_lo_fl, L_lo_fl}).hi;
    i128 r_S2 = (i128)int_part * S2 + div_cl((i128)p * S2, q);
    i128 B0_lo = (i128)(m0 - 1) * et.lo + S2 - r_S2;
    assert(B0_lo > 2 * S2);
    if (int_part == 0) {
        int Vcut = (int)min((i128)24, B0_lo / (2 * S2));
        i128 alpha_lo = S2 - div_cl((i128)Vcut * S2 * S2, 2 * B0_lo);
        i128 w_eff_lo = div_fl(L_lo_fl * S2, alpha_lo);
        i128 inv_alpha_hi = div_cl((i128)S2 * S2, alpha_lo);
        i128 alpha_pow_hi = iv_pow_rat({inv_alpha_hi, inv_alpha_hi}, q - p, q).hi;
        i128 main_part = div_cl(alpha_pow_hi * exp_gamma_upper_hi_frac(p, q, w_eff_lo), S2);
        i128 base_tail_lo = L_lo_fl + alpha_lo * Vcut;
        i128 inv_bt_hi = div_cl((i128)S2 * S2, base_tail_lo);
        i128 bt_pow_hi = iv_pow_rat({inv_bt_hi, inv_bt_hi}, q - p, q).hi;
        i128 evcut_hi = iv_exp_neg({(i128)Vcut * S2, (i128)Vcut * S2}).hi;
        Iv gam = iv_gamma_rat(p, q);
        i128 tail_part = div_cl(div_cl(bt_pow_hi * evcut_hi, S2) * S2, gam.lo);
        i128 ans = div_cl(exp_q_hi * (main_part + tail_part), S2);
        return (int64_t)min((i128)S, div_cl(ans * S, S2));
    } else {
        Iv u0 = iv_from_rat(tau_num, tau_den * (m0 - 1));
        i128 ex_ratio_hi;
        if (u0.hi <= S2 / 2) {
            i128 u = u0.hi;
            i128 u2 = div_cl(u * u, S2);
            i128 u3 = div_cl(u2 * u, S2);
            ex_ratio_hi = S2 + div_cl(u, 2) + div_cl(u2, 6) + div_cl(u3, 12);
        } else {
            Iv eu0 = iv_exp_pos(u0);
            ex_ratio_hi = div_cl((eu0.hi - S2) * S2, u0.lo);
        }
        i128 w_hi = div_cl(L_lo_cl * ex_ratio_hi, S2);
        i128 pref_hi;
        if (int_part == 1) pref_hi = S2 + div_cl((i128)p * S2 * S2, (i128)q * B0_lo);
        else {
            i128 f1 = S2 + div_cl((i128)(p + q) * S2 * S2, (i128)q * B0_lo);
            i128 f2 = S2 + div_cl((i128)p * S2 * S2, (i128)q * B0_lo);
            pref_hi = div_cl(f1 * f2, S2);
        }
        i128 val = div_cl(pref_hi * exp_gamma_upper_hi_general(int_part, p, q, L_lo_fl, w_hi), S2);
        i128 ans = div_cl(exp_q_hi * val, S2);
        return (int64_t)min((i128)S, div_cl(ans * S, S2));
    }
}

static vector<vector<int64_t>> solve_static_hub_upper_int(int Kmax) {
    vector<vector<int64_t>> U(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    for (int j = 0; j < Kmax; j++) U[Kmax + 1][j] = (int64_t)div_cl(S, j + 1);
    U[Kmax + 1][Kmax] = min(S, (int64_t)div_cl(2 * S, Kmax + 1));
    for (int k = Kmax; k >= 2; k--) {
        for (int j = 0; j < k; j++) {
            int64_t lp = (j > 0) ? U[k][j - 1] : 0;
            int64_t gp = U[k + 1][j + 1];
            int64_t best = 0;
            for (int g2 : {2 * (j + 1), k + j + 1}) {
                i128 num = (i128)2 * S + (i128)2 * j * lp + (i128)k * g2 * gp;
                i128 den = 2 * (1 + j) + k * g2;
                int64_t cand = (int64_t)div_cl(num, den);
                if (cand > best) best = cand;
            }
            int64_t cap = (j <= k - 2) ? (int64_t)div_cl(S, j + 1) : min(S, (int64_t)div_cl(2 * S, j + 1));
            U[k][j] = min(best, cap);
        }
    }
    return U;
}

struct FullState {
    vector<vector<int64_t>> V, R, U;
    vector<int64_t> H, P;
};

static FullState init_at_200(int64_t ylo_n, int64_t ylo_d, int64_t yhi_n, int64_t yhi_d,
                             int M_steps, int Dmax, int Rmax, int Kmax, int Hcut,
                             const vector<vector<int64_t>>& U_hub) {
    const int m0 = 200;
    int64_t tlo_n = ylo_n * (m0 - 1), tlo_d = ylo_d;
    int64_t thi_n = yhi_n * (m0 - 1), thi_d = yhi_d;

    vector<vector<int64_t>> Vn(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
    vector<vector<int64_t>> Rn(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
    vector<int64_t> Hn(Hcut + 2, 0);
    vector<vector<int64_t>> gVn(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
    vector<vector<int64_t>> gRn(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
    vector<int64_t> gHn(Hcut + 2, 0);
    vector<vector<bool>> stopVn(Dmax + 2, vector<bool>(Dmax + 2, true));
    vector<bool> stopHn(Hcut + 2, true);

    for (int d = 2; d <= Dmax; d++)
        for (int j = 1; j < d; j++)
            Vn[d][j] = H_raw_lower_200(j, M_steps + 1, thi_n, thi_d, m0);
    for (int p = 2; p <= Hcut + 1; p++) Hn[p] = Vn[2][1];

    for (int i = M_steps - 1; i >= 0; i--) {
        int t_next = i + 1;
        i128 dt_cor_hi = div_cl((i128)S2 * S2, (i128)m0 * EXP_GRID[t_next].lo - S2);
        vector<vector<int64_t>> Vc(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
        vector<vector<int64_t>> Rc(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
        vector<int64_t> Hc(Hcut + 2, 0);
        vector<vector<int64_t>> gVc(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
        vector<vector<int64_t>> gRc(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
        vector<int64_t> gHc(Hcut + 2, 0);
        vector<vector<bool>> stopVc(Dmax + 2, vector<bool>(Dmax + 2, false));
        vector<bool> stopHc(Hcut + 2, false);

        vector<int64_t> H_floor(Dmax, 0);
        for (int j = 1; j < Dmax; j++) H_floor[j] = H_raw_lower_200(j, t_next, thi_n, thi_d, m0);
        for (int j = 1; j < Dmax; j++) {
            assert(H_floor[j] >= Vn[Dmax][j]);
            Vc[Dmax][j] = H_floor[j];
            stopVc[Dmax][j] = true;
            gVc[Dmax][j] = (int64_t)((i128)(Vc[Dmax][j] - Vn[Dmax][j]) * H_DEN);
        }
        for (int d = Dmax - 1; d >= 1; d--) {
            for (int k = Rmax; k >= 2; k--) {
                i128 kd = (i128)k * d, r1 = d, r2 = (i128)k * (d * (k - 1) + 1), rtot = r1 + r2;
                i128 num_R = r1 * Vn[d + 1][1] + r2 * Rn[d][k + 1] - rtot * Rn[d][k];
                i128 g_R = (num_R > 0) ? div_fl(num_R, kd) : 0;
                i128 num_L = (r1 * Vc[d + 1][1] + r2 * Rc[d][k + 1] - rtot * Rn[d][k]) * H_DEN;
                i128 g_L = (num_L > 0) ? div_fl(num_L, kd * H_DEN + rtot) : 0;
                int64_t delta = (int64_t)div_fl(max((i128)0, min(g_R, g_L)), H_DEN);
                Rc[d][k] = Rn[d][k] + delta;
                gRc[d][k] = (int64_t)((i128)delta * H_DEN);
            }
            if (d >= 2) {
                for (int j = 1; j < d; j++) {
                    int64_t lt_n = (j > 1) ? Vn[d][j - 1] : Rn[d][2];
                    int64_t lt_c = (j > 1) ? Vc[d][j - 1] : Rc[d][2];
                    int64_t glt_n = (j > 1) ? gVn[d][j - 1] : gRn[d][2];
                    int64_t gt_n = Vn[d + 1][j + 1], gt_c = Vc[d + 1][j + 1], ggt_n = gVn[d + 1][j + 1];
                    i128 rg = (i128)d * j, rl = j, rtot = rg + rl;
                    bool can_cont = true;
                    if (stopVn[d][j]) {
                        i128 bk = -(i128)d * gVn[d][j] + rg * (gt_n - Vn[d][j]) + rl * (lt_n - Vn[d][j]);
                        i128 sd = rg * (ggt_n - gVn[d][j]) + rl * (glt_n - gVn[d][j]);
                        i128 pen = (sd > 0) ? div_cl(sd * dt_cor_hi, S2) : 0;
                        if (bk - pen < 0) can_cont = false;
                    }
                    int64_t v_ode = -1;
                    if (can_cont) {
                        i128 num_R = rg * gt_n + rl * lt_n - rtot * Vn[d][j];
                        i128 num_L = (rg * gt_c + rl * lt_c - rtot * Vn[d][j]) * H_DEN;
                        if (num_R >= 0 && num_L >= 0) {
                            i128 gamma = min(div_fl(num_R, d), div_fl(num_L, (i128)d * H_DEN + rtot));
                            v_ode = Vn[d][j] + (int64_t)div_fl(gamma, H_DEN);
                        } else can_cont = false;
                    }
                    if (!can_cont || H_floor[j] > v_ode) {
                        assert(H_floor[j] >= Vn[d][j]);
                        Vc[d][j] = H_floor[j]; stopVc[d][j] = true;
                    } else {
                        Vc[d][j] = v_ode; stopVc[d][j] = false;
                    }
                    gVc[d][j] = (int64_t)((i128)(Vc[d][j] - Vn[d][j]) * H_DEN);
                }
            }
        }
        Hc[Hcut + 1] = Vc[2][1]; gHc[Hcut + 1] = gVc[2][1]; stopHc[Hcut + 1] = stopVc[2][1];
        for (int p = Hcut; p >= 2; p--) {
            i128 rg = 2 * (p + 1), rl = p;
            bool can_cont = true;
            i128 best_g = (i128)1 << 120;
            for (i128 rh : {(i128)0, (i128)p * (2 * p - 1)}) {
                i128 rtot = rg + rl + rh;
                if (stopHn[p]) {
                    i128 bk = -(i128)2 * p * gHn[p] + rg * (Vn[3][2] - Hn[p]) + rl * (Rn[2][2] - Hn[p]) + rh * (Hn[p + 1] - Hn[p]);
                    i128 sd = rg * (gVn[3][2] - gHn[p]) + rl * (gRn[2][2] - gHn[p]) + rh * (gHn[p + 1] - gHn[p]);
                    i128 pen = (sd > 0) ? div_cl(sd * dt_cor_hi, S2) : 0;
                    if (bk - pen < 0) can_cont = false;
                }
                i128 num_R = rg * Vn[3][2] + rl * Rn[2][2] + rh * Hn[p + 1] - rtot * Hn[p];
                i128 num_L = (rg * Vc[3][2] + rl * Rc[2][2] + rh * Hc[p + 1] - rtot * Hn[p]) * H_DEN;
                if (num_R < 0 || num_L < 0) can_cont = false;
                else {
                    i128 gc = min(div_fl(num_R, 2 * p), div_fl(num_L, (i128)2 * p * H_DEN + rtot));
                    if (gc < best_g) best_g = gc;
                }
            }
            int64_t h_ode = can_cont ? (Hn[p] + (int64_t)div_fl(best_g, H_DEN)) : -1;
            if (!can_cont || H_floor[1] > h_ode) {
                assert(H_floor[1] >= Hn[p]);
                Hc[p] = H_floor[1]; stopHc[p] = true;
            } else {
                Hc[p] = h_ode; stopHc[p] = false;
            }
            gHc[p] = (int64_t)((i128)(Hc[p] - Hn[p]) * H_DEN);
        }
        for (int d = 2; d < Dmax; d++)
            for (int j = 1; j < d; j++) assert(Vc[d + 1][j + 1] >= Vc[d][j]);
        for (int d = 1; d < Dmax; d++) {
            for (int k = 2; k < Rmax; k++) assert(Rc[d][k] >= Rc[d][k + 1]);
            assert(Vc[d + 1][1] >= Rc[d][2]);
        }
        for (int p = 2; p <= Hcut; p++) assert(Hc[p] >= Hc[p + 1]);
        Vn.swap(Vc); Rn.swap(Rc); Hn.swap(Hc);
        gVn.swap(gVc); gRn.swap(gRc); gHn.swap(gHc);
        stopVn.swap(stopVc); stopHn.swap(stopHc);
    }

    FullState st;
    st.V = move(Vn); st.R = move(Rn); st.H = move(Hn);
    if (tlo_n == 0) {
        st.U = U_hub;
        st.P.assign(Kmax + 2, S);
        return st;
    }

    vector<int64_t> P_stop(Kmax + 2, 0);
    vector<vector<int64_t>> U_stop(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    for (int k = 2; k <= Kmax + 1; k++) {
        P_stop[k] = B_raw_upper_200(1, 1, k, M_steps, tlo_n, tlo_d, m0);
        int64_t cap = B_raw_upper_200(0, 1, k, M_steps, tlo_n, tlo_d, m0);
        for (int j = 0; j < k; j++) U_stop[k][j] = min(cap, U_hub[k][j]);
    }
    for (int k = Kmax; k >= 2; k--) P_stop[k] = max(P_stop[k], P_stop[k + 1]);

    i128 dt_term_hi = div_cl((i128)S2 * S2, (i128)(m0 - 1) * EXP_GRID[M_steps - 1].lo);
    vector<int64_t> Pn(Kmax + 2, 0);
    vector<vector<int64_t>> Un(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    Pn[Kmax + 1] = max(P_stop[Kmax + 1], B_raw_upper_200(1, 1, Kmax + 1, M_steps - 1, tlo_n, tlo_d, m0));
    {
        int64_t cap_m1 = B_raw_upper_200(0, 1, Kmax + 1, M_steps - 1, tlo_n, tlo_d, m0);
        for (int j = 0; j <= Kmax; j++) Un[Kmax + 1][j] = max(U_stop[Kmax + 1][j], min(cap_m1, U_hub[Kmax + 1][j]));
    }
    for (int k = Kmax; k >= 2; k--) {
        int64_t bg_M = B_raw_upper_200(2, 1, k, M_steps, tlo_n, tlo_d, m0);
        int64_t bl_M = (k == 2) ? B_raw_upper_200(1, 0, 1, M_steps, tlo_n, tlo_d, m0)
                                : B_raw_upper_200(0, k + 2, 2 * k, M_steps, tlo_n, tlo_d, m0);
        i128 md = 0;
        for (int gw2 : {1, 2 * k - 1}) {
            i128 num = (i128)2 * (k + 1) * (bg_M - P_stop[k]) + (i128)k * (bl_M - P_stop[k]) + (i128)k * gw2 * (P_stop[k + 1] - P_stop[k]);
            if (num > 0) md = max(md, div_cl(num, 2 * k));
        }
        Pn[k] = max(P_stop[k] + (int64_t)div_cl(md * dt_term_hi, S2), Pn[k + 1]);
    }
    for (int k = Kmax; k >= 2; k--) {
        for (int j = 0; j < k; j++) {
            int64_t lp = (j > 0) ? U_stop[k][j - 1] : 0;
            i128 md = 0;
            for (int g2 : {2 * (j + 1), k + j + 1}) {
                i128 num = (i128)2 * (P_stop[k] - U_stop[k][j]) + (i128)2 * j * (lp - U_stop[k][j]) + (i128)k * g2 * (U_stop[k + 1][j + 1] - U_stop[k][j]);
                if (num > 0) md = max(md, div_cl(num, 2 * k));
            }
            Un[k][j] = U_stop[k][j] + (int64_t)div_cl(md * dt_term_hi, S2);
        }
    }
    for (int i = M_steps - 1; i >= 0; i--) {
        vector<vector<int64_t>> Uc(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
        vector<int64_t> Pc(Kmax + 2, 0);
        int t_bnd = max(0, i - 1);
        Pc[Kmax + 1] = max(Pn[Kmax + 1], B_raw_upper_200(1, 1, Kmax + 1, t_bnd, tlo_n, tlo_d, m0));
        int64_t cap_top = B_raw_upper_200(0, 1, Kmax + 1, t_bnd, tlo_n, tlo_d, m0);
        for (int j = 0; j <= Kmax; j++) Uc[Kmax + 1][j] = max(Un[Kmax + 1][j], min(cap_top, U_hub[Kmax + 1][j]));
        for (int k = Kmax; k >= 2; k--) {
            int64_t bg_k = B_raw_upper_200(2, 1, k, i, tlo_n, tlo_d, m0);
            int64_t bl_k = (k == 2) ? B_raw_upper_200(1, 0, 1, i, tlo_n, tlo_d, m0)
                                    : B_raw_upper_200(0, k + 2, 2 * k, i, tlo_n, tlo_d, m0);
            i128 beta = 0;
            for (int gw2 : {1, 2 * k - 1}) {
                i128 rg = 2 * (k + 1), rl = k, rh = (i128)k * gw2, rtot = rg + rl + rh;
                i128 nR = rg * bg_k + rl * bl_k + rh * Pn[k + 1] - rtot * Pn[k];
                if (nR > 0) beta = max(beta, div_cl(nR, 2 * k));
                i128 nL = (rg * bg_k + rl * bl_k + rh * Pc[k + 1] - rtot * Pn[k]) * H_DEN;
                if (nL > 0) beta = max(beta, div_cl(nL, (i128)2 * k * H_DEN + rtot));
            }
            Pc[k] = max(Pn[k] + (int64_t)div_cl(beta, H_DEN), Pc[k + 1]);
        }
        for (int k = Kmax; k >= 2; k--) {
            for (int j = 0; j < k; j++) {
                int64_t ln = (j > 0) ? Un[k][j - 1] : 0, lc = (j > 0) ? Uc[k][j - 1] : 0;
                i128 beta = 0;
                for (int g2 : {2 * (j + 1), k + j + 1}) {
                    i128 rtot = 2 * (1 + j) + (i128)k * g2;
                    i128 nR = (i128)2 * Pn[k] + (i128)2 * j * ln + (i128)k * g2 * Un[k + 1][j + 1] - rtot * Un[k][j];
                    if (nR > 0) beta = max(beta, div_cl(nR, 2 * k));
                    i128 nL = ((i128)2 * Pc[k] + (i128)2 * j * lc + (i128)k * g2 * Uc[k + 1][j + 1] - rtot * Un[k][j]) * H_DEN;
                    if (nL > 0) beta = max(beta, div_cl(nL, (i128)2 * k * H_DEN + rtot));
                }
                Uc[k][j] = Un[k][j] + (int64_t)div_cl(beta, H_DEN);
            }
        }
        Un.swap(Uc); Pn.swap(Pc);
    }
    for (int k = 2; k <= Kmax + 1; k++) {
        Pn[k] = min(Pn[k], B_raw_upper_200(1, 1, k, 0, tlo_n, tlo_d, m0));
        int64_t cap0 = B_raw_upper_200(0, 1, k, 0, tlo_n, tlo_d, m0);
        for (int j = 0; j < k; j++) Un[k][j] = min({Un[k][j], cap0, U_hub[k][j]});
    }
    st.U = move(Un); st.P = move(Pn);
    return st;
}

// --- Exact integer-N 1-step backward induction (N+1 -> N) ---
static void step_backward_N(FullState& st, int N, int64_t ylo_n, int64_t ylo_d, int64_t yhi_n, int64_t yhi_d,
                            int Dmax, int Rmax, int Kmax, int Hcut, const vector<vector<int64_t>>& U_hub) {
    // 1. Lower bound step (exact at size N, y_hi)
    vector<vector<int64_t>> Vc(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
    vector<vector<int64_t>> Rc(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
    vector<int64_t> Hc(Hcut + 2, 0);

    vector<int64_t> H_fl(Dmax + 1, 0);
    for (int j = 1; j < min(Dmax, N); j++) H_fl[j] = H_exact_N_lower(N, j, yhi_n, yhi_d);
    if (Dmax <= N - 1) {
        for (int j = 1; j < Dmax; j++) Vc[Dmax][j] = H_fl[j];
    }

    int d_lim = min(Dmax - 1, N - 1);
    for (int d = d_lim; d >= 1; d--) {
        int k_lim = min(Rmax, N - d);
        for (int k = k_lim; k >= 2; k--) {
            i128 kd = (i128)k * d;
            i128 r1 = d, r2 = (i128)k * (d * (k - 1) + 1);
            i128 rstay = (i128)N * kd - r1 - r2;
            assert(rstay >= 0);
            i128 num = rstay * st.R[d][k] + r1 * st.V[d + 1][1] + r2 * st.R[d][k + 1];
            Rc[d][k] = (int64_t)div_fl(num, (i128)N * kd);
        }
        if (d >= 2) {
            for (int j = 1; j < d; j++) {
                int64_t lt = (j > 1) ? st.V[d][j - 1] : st.R[d][2];
                int64_t gt = st.V[d + 1][j + 1];
                i128 rg = (i128)d * j, rl = j;
                i128 rstay = (i128)N * d - rg - rl;
                assert(rstay >= 0);
                int64_t v_step = (int64_t)div_fl(rstay * st.V[d][j] + rg * gt + rl * lt, (i128)N * d);
                Vc[d][j] = max(v_step, H_fl[j]);
            }
        }
    }
    int p_lim = min(Hcut, N - 2);
    for (int p = Hcut + 1; p > p_lim; p--) Hc[p] = Vc[2][1];
    for (int p = p_lim; p >= 2; p--) {
        i128 rg = 2 * (p + 1), rl = p;
        int64_t best = S;
        for (i128 rh : {(i128)0, (i128)p * (2 * p - 1)}) {
            i128 rstay = (i128)2 * p * N - rg - rl - rh;
            assert(rstay >= 0);
            int64_t cand = (int64_t)div_fl(rstay * st.H[p] + rg * st.V[3][2] + rl * st.R[2][2] + rh * st.H[p + 1], (i128)2 * p * N);
            if (cand < best) best = cand;
        }
        Hc[p] = max(best, H_fl[1]);
    }

    // Verify exact pathwise monotonicity on ALL reachable states at size N!
    for (int d = 2; d <= min(Dmax - 1, N - 2); d++)
        for (int j = 1; j < d; j++) assert(Vc[d + 1][j + 1] >= Vc[d][j]);
    for (int d = 1; d <= min(Dmax - 1, N - 2); d++) {
        int kl = min(Rmax, N - d);
        for (int k = 2; k < kl; k++) assert(Rc[d][k] >= Rc[d][k + 1]);
        if (kl >= 2) assert(Vc[d + 1][1] >= Rc[d][2]);
    }
    for (int p = 2; p < p_lim; p++) assert(Hc[p] >= Hc[p + 1]);
    if (p_lim >= 2 && p_lim == Hcut) assert(Hc[Hcut] >= Hc[Hcut + 1]);

    st.V.swap(Vc); st.R.swap(Rc); st.H.swap(Hc);

    // 2. Upper bound step (exact at size N, y_lo)
    if (ylo_n == 0) {
        // Static U_hub is already exact at y_lo = 0!
        return;
    }
    vector<vector<int64_t>> Uc(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    vector<int64_t> Pc(Kmax + 2, 0);

    if (Kmax + 1 <= N - 2) {
        Pc[Kmax + 1] = B_exact_N_upper(N, 1, 1, Kmax + 1, ylo_n, ylo_d);
    }
    if (Kmax + 1 <= N - 1) {
        int64_t cap = B_exact_N_upper(N, 0, 1, Kmax + 1, ylo_n, ylo_d);
        for (int j = 0; j <= Kmax; j++) Uc[Kmax + 1][j] = min(cap, U_hub[Kmax + 1][j]);
    }

    int pk_lim = min(Kmax, N - 2);
    for (int k = pk_lim; k >= 2; k--) {
        int64_t bg = B_exact_N_upper(N + 1, 2, 1, k, ylo_n, ylo_d);
        int64_t bl = (k == 2) ? B_exact_N_upper(N + 1, 1, 0, 1, ylo_n, ylo_d)
                              : B_exact_N_upper(N + 1, 0, k + 2, 2 * k, ylo_n, ylo_d);
        int64_t best = 0;
        for (int gw2 : {1, 2 * k - 1}) {
            i128 rg = 2 * (k + 1), rl = k, rh = (i128)k * gw2;
            i128 rstay = (i128)2 * k * N - rg - rl - rh;
            assert(rstay >= 0);
            int64_t cand = (int64_t)div_cl(rstay * st.P[k] + rg * bg + rl * bl + rh * st.P[k + 1], (i128)2 * k * N);
            if (cand > best) best = cand;
        }
        int64_t cap = B_exact_N_upper(N, 1, 1, k, ylo_n, ylo_d);
        Pc[k] = min(best, cap);
    }

    int uk_lim = min(Kmax, N - 1);
    for (int k = uk_lim; k >= 2; k--) {
        int64_t cap = B_exact_N_upper(N, 0, 1, k, ylo_n, ylo_d);
        for (int j = 0; j < k; j++) {
            int64_t lp = (j > 0) ? st.U[k][j - 1] : 0;
            int64_t best = 0;
            for (int g2 : {2 * (j + 1), k + j + 1}) {
                i128 rp = 2, rl = 2 * j, rg = (i128)k * g2;
                i128 rstay = (i128)2 * k * N - rp - rl - rg;
                assert(rstay >= 0);
                // Note: if k == N - 1, st.U[k+1][j+1] is at N+1 where k+1 == N == (N+1)-1 (reachable!)
                int64_t cand = (int64_t)div_cl(rstay * st.U[k][j] + rp * st.P[k] + rl * lp + rg * st.U[k + 1][j + 1], (i128)2 * k * N);
                if (cand > best) best = cand;
            }
            Uc[k][j] = min({best, cap, U_hub[k][j]});
        }
    }
    st.U.swap(Uc); st.P.swap(Pc);
}

// --- Exact All-Tree Signed Ledger for fixed m in {4, ..., 199} (max degree in T_{m-2} is m - 3!) ---
struct Tariffs {
    int64_t a, b_base, fee_lo, fee_hi, fee_t, xfer;
    int D_hi;
};

static bool verify_ledger_for_m(const FullState& st, int m, int64_t ylo_n, int64_t ylo_d,
                                int Rmax, int Kmax, int Hcut, string& fail_reason) {
    int max_deg = m - 3; // max possible degree of any vertex in T_n (n = m - 2)!
    int64_t M_lo = st.V[2][1];
    if (M_lo <= 0) { fail_reason = "M_lo<=0"; return false; }

    // Special exact check for m = 4 (T_2 = single edge) AND Founder Seam b = 2 (Z_2 >_st Z_3 via S_4 & P_4)!
    if (m == 4) {
        if (st.H[2] > st.U[2][0] && st.V[3][2] > st.U[3][2]) return true;
        fail_reason = (st.H[2] <= st.U[2][0]) ? "m4_seed" : "founder_S4";
        return false;
    }

    auto H = [&](int p) -> int64_t {
        if (p <= min(Hcut, m - 2)) return (int64_t)div_fl((i128)max(st.H[p], M_lo) * LSC, M_lo);
        return LSC;
    };
    auto L = [&](int k) -> int64_t {
        if (k <= min(Rmax, m - 1)) return (int64_t)div_fl((i128)st.R[1][k] * LSC, M_lo);
        return 0;
    };

    int64_t u_tail_beta = (ylo_n > 0 && Kmax + 2 < m) ? (int64_t)div_cl((i128)B_exact_N_upper(m, 0, 1, Kmax + 2, ylo_n, ylo_d) * LSC, M_lo)
                                                      : (int64_t)div_cl((i128)S * LSC, M_lo);
    int64_t u_ktop_0 = (Kmax + 1 <= m - 1) ? (int64_t)div_cl((i128)st.U[Kmax + 1][0] * LSC, M_lo) : u_tail_beta;
    int64_t u_tail0_cap = min(u_tail_beta, u_ktop_0);

    int k_avail = min(Kmax + 1, m - 1);
    vector<vector<int64_t>> Utab(k_avail + 1, vector<int64_t>(k_avail + 1, 0));
    for (int k = 2; k <= k_avail; k++)
        for (int j = 0; j < k; j++)
            Utab[k][j] = (int64_t)div_cl((i128)st.U[k][j] * LSC, M_lo);

    auto U = [&](int k, int j) -> int64_t {
        if (k <= k_avail) return Utab[k][min(j, k - 1)];
        int64_t c_hub = (int64_t)div_cl((i128)S * LSC, (i128)(j + 1) * M_lo);
        return min(c_hub, u_tail0_cap);
    };

    int64_t u20 = U(2, 0);
    if (H(2) <= u20) { fail_reason = "seed"; return false; }

    auto site_credit_lo = [&](int d, int c, i128 r_num, i128 r_den, i128& s_num_out, i128& s_den_out) -> i128 {
        int l = d - c;
        i128 sn = (i128)l * r_den + r_num;
        i128 sd = r_den * d * (d + 1);
        s_num_out = sn; s_den_out = sd;
        int64_t hd = H(d + 1) - u20;
        int64_t h2 = H(2) - u20;
        i128 C = div_fl((i128)d * sn * hd, sd) + div_fl((i128)l * h2, 2 * d);
        if (l > 0) {
            i128 diff = sd - 2 * d * sn;
            if (diff > 0) {
                int64_t slope = max((int64_t)0, u20 - L(d + 1));
                C -= div_cl((i128)l * diff * slope, 2 * d * sd);
            } else if (diff < 0) {
                int64_t slope = max((int64_t)0, U(d + 1, l - 1) - L(2));
                C -= div_cl((i128)l * (-diff) * slope, 2 * d * sd);
            }
        }
        return C;
    };

    // Special exact check for m = 5 (T_3 = path P_3 = star_2 ONLY!):
    if (m == 5) {
        i128 sn, sd;
        if (site_credit_lo(2, 0, 0, 1, sn, sd) > 0) return true;
        fail_reason = "m5_star2";
        return false;
    }

    // One-line universal inequality using true global U3_max = max_{k>=3, 0<=j<k} U(k, j) and u20 >= U3_max!
    int64_t U3_max = (max_deg + 1 > k_avail) ? u_tail0_cap : 0;
    for (int k = 3; k <= k_avail; k++)
        for (int j = 0; j < k; j++)
            U3_max = max(U3_max, Utab[k][j]);
    if (u20 >= U3_max && (i128)3 * LSC > (i128)4 * u20 + (i128)2 * U3_max) return true;

    vector<Tariffs> tlist = {
        { 9200000, 22500000, 4800000, 8200000, 27500000, 4600000, 10},
        { 8800000, 22500000, 4800000, 8100000, 27500000, 4800000, 10},
        {10000000, 22000000, 5200000, 8200000, 27500000, 4800000, 10},
        {11000000, 22000000, 5500000, 8400000, 27500000, 5000000, 10},
        {12500000, 21500000, 6000000, 8600000, 27500000, 5500000, 10},
        {15000000, 20000000, 7500000, 9000000, 27500000, 6000000, 10},
        {20000000, 20000000,10000000,12000000, 35000000, 8000000, 10},
        {30000000, 25000000,15000000,15000000, 45000000,10000000, 10},
    };

    string last_fail;
    for (size_t ti = 0; ti < tlist.size(); ti++) {
        const auto& T = tlist[ti];
        bool ok = true;

        // 1. Stars d = 2..max_deg
        for (int d = 2; d <= max_deg; d++) {
            i128 sn, sd;
            if (site_credit_lo(d, 0, 0, 1, sn, sd) <= 0) {
                ok = false; last_fail = "star_" + to_string(d); break;
            }
        }
        if (!ok) continue;

        // 2. Finite sites d = 2..max_deg, c = 1..d (excluding bare deg-2 d=c=2)
        // Note: in a tree of size n = m - 2, a vertex of degree d with c non-leaf neighbors
        // has at least d + 1 + c <= n vertices, so c <= n - 1 - d = max_deg - d!
        for (int d = 2; d <= max_deg && ok; d++) {
            int c_max = min(d, max_deg - d); // exact! each non-leaf neighbor adds >= 1 vertex beyond the d+1 closed neighborhood!
            for (int c = 1; c <= c_max && ok; c++) {
                if (d == 2 && c == 2) continue;
                int l = d - c;
                int64_t u_adv = U(d + 1, l);
                if (d >= T.D_hi) {
                    i128 s_max_num = 2 * l + c;
                    i128 s_max_den = (i128)2 * d * (d + 1);
                    if (div_cl(s_max_num * u_adv, s_max_den) <= T.b_base) {
                        int verts[3][2] = {{0, 0}, {c, 0}, {0, c}};
                        for (int vi = 0; vi < 3 && ok; vi++) {
                            int k = verts[vi][0], pp = verts[vi][1];
                            i128 r_nums[2] = {2 * k + pp, 2 * c};
                            for (int ri = 0; ri < 2 && ok; ri++) {
                                i128 sn, sd;
                                i128 C = site_credit_lo(d, c, r_nums[ri], 4, sn, sd);
                                i128 R = C + (i128)T.a * (c - 2) - T.b_base - (i128)k * T.fee_hi - (i128)pp * T.xfer;
                                if (R < 0) { ok = false; last_fail = "hi_d" + to_string(d) + "_c" + to_string(c); }
                            }
                        }
                    } else {
                        for (int k = 0; k <= c && ok; k++) {
                            int pps[2] = {0, c - k};
                            for (int pi = 0; pi < 2 && ok; pi++) {
                                int pp = pps[pi];
                                i128 r_nums[2] = {2 * k + pp, 2 * c};
                                for (int ri = 0; ri < 2 && ok; ri++) {
                                    i128 sn, sd;
                                    i128 C = site_credit_lo(d, c, r_nums[ri], 4, sn, sd);
                                    i128 R = C + (i128)T.a * (c - 2) - T.b_base - (i128)k * T.fee_hi - (i128)pp * T.xfer;
                                    i128 exc = (i128)(c - k) * max((i128)0, div_cl(sn * u_adv, sd) - T.b_base);
                                    if (R < exc) { ok = false; last_fail = "hi_d" + to_string(d) + "_c" + to_string(c); }
                                }
                            }
                        }
                    }
                } else {
                    int64_t F_ch = (d == 2) ? T.fee_t : T.fee_lo;
                    bool coarse_ok = true;
                    for (int k = 0; k <= c && coarse_ok; k++) {
                        i128 r_nums[2] = {k, c};
                        for (int ri = 0; ri < 2 && coarse_ok; ri++) {
                            i128 sn, sd;
                            i128 C = site_credit_lo(d, c, r_nums[ri], 2, sn, sd);
                            i128 R = C + (i128)T.a * (c - 2) - T.b_base - (i128)k * F_ch;
                            i128 exc = (i128)(c - k) * max((i128)0, div_cl(sn * u_adv, sd) - T.b_base);
                            if (R < exc) coarse_ok = false;
                        }
                    }
                    if (coarse_ok) continue;

                    const i128 LCM = 2520;
                    vector<int> es;
                    int e_hi_bound = min(T.D_hi - 1, max_deg);
                    bool has_hi_nbr = (max_deg >= T.D_hi);
                    function<void(int, int, int)> dfs = [&](int k, int rem, int start_idx) {
                        if (!ok) return;
                        if (rem == 0) {
                            int n_hi = 0;
                            i128 r_known_num = (i128)k * (LCM / 2);
                            for (int e : es) {
                                if (e == 0) n_hi++;
                                else r_known_num += LCM / e;
                            }
                            i128 r_unk_opts[2] = {0, (i128)n_hi * (LCM / T.D_hi)};
                            for (int ui = 0; ui < 2 && ok; ui++) {
                                i128 sn, sd;
                                i128 C = site_credit_lo(d, c, r_known_num + r_unk_opts[ui], LCM, sn, sd);
                                i128 R = C + (i128)T.a * (c - 2) - T.b_base - (i128)k * F_ch;
                                if (d <= 4) R += (i128)n_hi * T.xfer;
                                i128 exc = 0;
                                for (int e : es) {
                                    if (e == 0) {
                                        exc += max((i128)0, div_cl(sn * u_adv, sd) - T.b_base);
                                    } else {
                                        i128 diff_n = sn * e * (e + 1) - LCM * (d + 1);
                                        if (diff_n > 0) {
                                            int64_t slope = max((int64_t)0, u_adv - L(e + 1));
                                            i128 j_uv = div_cl(diff_n * slope, sd * e * (e + 1));
                                            exc += max((i128)0, j_uv - T.b_base);
                                        }
                                    }
                                }
                                if (R < exc) { ok = false; last_fail = "lo_d" + to_string(d) + "_c" + to_string(c); }
                            }
                            return;
                        }
                        for (int idx = start_idx; idx <= e_hi_bound; idx++) {
                            if (idx == 0 && !has_hi_nbr) continue;
                            if (idx == 1) continue;
                            es.push_back(idx);
                            dfs(k, rem - 1, idx);
                            es.pop_back();
                            if (!ok) return;
                        }
                    };
                    for (int k = 0; k <= c && ok; k++) dfs(k, c - k, 0);
                }
            }
        }
        if (!ok) continue;

        // 3. Bare degree-2 chains of lengths 1, 2, 3 (only if a bare deg-2 vertex can even exist: needs n >= 5, i.e. m >= 7!)
        if (m >= 7) {
            const i128 X = 600000000LL;
            struct End {
                i128 q_lo, q_hi, s_lo, s_hi;
                int64_t w_ch;
                i128 own_or_we, fee_XL;
                bool is_cap;
            };
            vector<End> ends;
            int64_t u30 = U(3, 0), L3 = L(3);
            for (int d = 2; d <= min(T.D_hi - 1, max_deg); d++) {
                int64_t F_ch = (d == 2) ? T.fee_t : T.fee_lo;
                i128 ql = div_fl(X, 6 * d), qh = div_cl(X, 6 * d);
                int64_t wch = max((int64_t)0, u30 - L(d + 1));
                int c_max = min(d, max_deg - d);
                for (int c = 1; c <= c_max; c++) {
                    if (d == 2 && c == 2) continue;
                    int l = d - c;
                    int64_t we = max((int64_t)0, U(d + 1, l) - L3);
                    for (int r2 : {1, c}) {
                        i128 sl = div_fl((i128)(2 * l + r2) * X, 2 * d * (d + 1));
                        i128 sh = div_cl((i128)(2 * l + r2) * X, 2 * d * (d + 1));
                        ends.push_back({ql, qh, sl, sh, wch, we, (i128)F_ch * X, false});
                    }
                }
            }
            for (int d = T.D_hi; d <= min(Kmax, max_deg); d++) {
                i128 ql = div_fl(X, 6 * d), qh = div_cl(X, 6 * d);
                i128 sl = div_fl(X, 2 * d * (d + 1));
                int64_t wch = max((int64_t)0, u30 - L(d + 1));
                i128 own_cap = 0;
                for (int l = 0; l < d; l++) {
                    int64_t we = max((int64_t)0, U(d + 1, l) - L3);
                    own_cap = max(own_cap, div_cl((i128)(d + l) * X * we, 2 * d * (d + 1)));
                }
                ends.push_back({ql, qh, sl, 0, wch, own_cap, (i128)T.fee_hi * X, true});
            }
            if (max_deg > Kmax) {
                i128 oc1 = div_cl((i128)u_tail0_cap * X, 2 * (Kmax + 2));
                i128 oc2a = div_cl((i128)(Kmax + 2) * S * LSC * X, (i128)4 * (Kmax + 1) * (Kmax + 2) * M_lo);
                i128 oc2b = div_cl((i128)u_tail0_cap * X, Kmax + 2);
                i128 oc_tail = max(oc1, min(oc2a, oc2b));
                ends.push_back({0, 0, 0, 0, u30, oc_tail, (i128)T.fee_hi * X, true});
                i128 qtl = div_fl(X, 6 * (Kmax + 1)), qth = div_cl(X, 6 * (Kmax + 1));
                ends.push_back({qtl, qth, 0, 0, u30, oc_tail, (i128)T.fee_hi * X, true});
            }
            int64_t v_int = max((int64_t)0, u30 - L3);
            int64_t cred_site = 2 * (H(3) - u20);
            i128 base_XL = (i128)T.b_base * X;
            auto port_debt = [&](i128 sa_lo, i128 sa_hi, const End& e) -> i128 {
                i128 d1 = (sa_hi > e.s_lo) ? (sa_hi - e.s_lo) * e.w_ch : 0;
                i128 d2 = e.is_cap ? e.own_or_we : ((e.s_hi > sa_lo) ? (e.s_hi - sa_lo) * e.own_or_we : 0);
                return max({(i128)0, d1, d2}) - e.fee_XL;
            };
            int ne = (int)ends.size();
            i128 X12 = X / 12, X6 = X / 6;
            for (int len = 1; len <= 3 && ok; len++) {
                if (m - 2 < len + 4) break; // a chain of length `len` between two non-bare non-leaf endpoints needs >= len + 4 vertices!
                for (int ix = 0; ix < ne && ok; ix++) {
                    const auto& ex = ends[ix];
                    for (int iy = ix; iy < ne && ok; iy++) {
                        const auto& ey = ends[iy];
                        i128 debt;
                        if (len == 1) {
                            i128 s1l = ex.q_lo + ey.q_lo, s1h = ex.q_hi + ey.q_hi;
                            debt = port_debt(s1l, s1h, ex) + port_debt(s1l, s1h, ey) - s1l * cred_site;
                        } else if (len == 2) {
                            i128 s1l = X12 + ex.q_lo, s1h = X12 + ex.q_hi;
                            i128 s2l = X12 + ey.q_lo, s2h = X12 + ey.q_hi;
                            debt = port_debt(s1l, s1h, ex) + port_debt(s2l, s2h, ey)
                                 + max(s2h - s1l, s1h - s2l) * v_int - (s1l + s2l) * cred_site;
                        } else {
                            i128 s1l = X12 + ex.q_lo, s1h = X12 + ex.q_hi;
                            i128 s2l = X6, s2h = X6;
                            i128 s3l = X12 + ey.q_lo, s3h = X12 + ey.q_hi;
                            debt = port_debt(s1l, s1h, ex) + port_debt(s3l, s3h, ey)
                                 + (max(s2h - s1l, s1h - s2l) + max(s3h - s2l, s2h - s3l)) * v_int
                                 - (s1l + s2l + s3l) * cred_site;
                        }
                        if (debt >= base_XL) { ok = false; last_fail = "chain_len" + to_string(len); }
                    }
                }
            }
        }
        if (ok) return true;
    }
    fail_reason = last_fail;
    return false;
}

int main() {
    init_exp_grid(1600);
    init_log_int(2100);
    init_quad_nodes();

    // Exact rational boxes on y = log(1/(1-x)) in [0, 2], denominator YQ = 19900000!
    // (199 * 100000 so tau_200 = 199*y has exact rational denominator 100000!)
    const int64_t YQ = 19900000LL;
    const int64_t ONE = YQ; // y = 1.0

    // Regime 1: y in [0, 6/199] (tau_200 in [0, 6]), uses init_at_200 + backward sweep 199..4!
    // Regime 2: y in [6/199, 2.0] (tau_200 in [6, 398]), uses exact direct integer-N sweep from N_top down to 4!
    vector<pair<int64_t, int64_t>> boxes;
    // Regime 1: tau_200 endpoints (in units of 1/100000 = 1/YQ * 199)
    vector<int64_t> t200_pts = {
        0, 25, 100, 300, 800, 2000, 5000, 12000, 25000, 45000,
        75000, 110000, 150000, 200000, 260000, 330000, 410000, 500000, 600000
    };
    for (size_t i = 0; i + 1 < t200_pts.size(); i++) {
        boxes.push_back({t200_pts[i], t200_pts[i + 1]});
    }
    // Regime 2: y from 6/199 (= 600000/YQ) up to 2.0 (= 2*YQ/YQ)
    vector<int64_t> y_pts = {
        600000, 750000, 950000, 1200000, 1500000, 1900000, 2400000, 3000000,
        3800000, 4800000, 6000000, 7500000, 9500000, 12000000, 15000000,
        19000000, 24000000, 30000000, 2 * ONE
    };
    for (size_t i = 0; i + 1 < y_pts.size(); i++) {
        boxes.push_back({y_pts[i], y_pts[i + 1]});
    }
    reverse(boxes.begin(), boxes.end());

    // Precompute static U_hub for max Kmax = 32
    auto U_hub32 = solve_static_hub_upper_int(32);

    int n_pass_boxes = 0, n_fail_boxes = 0;
    long long total_m_checks = 0;
    int64_t covered_up_to = 0;

    while (!boxes.empty()) {
        auto [a, b] = boxes.back();
        boxes.pop_back();

        // Which m in {4, ..., 199} still need certification on [a/YQ, b/YQ]?
        // Need (m - 1) * (a / YQ) < 6  <=>  (m - 1) * a < 6 * YQ!
        int m_max = 199;
        if (a > 0) {
            i128 lim = ((i128)6 * YQ - 1) / a + 1;
            if (lim < m_max) m_max = (int)lim;
        }
        if (m_max < 4) {
            assert(a == covered_up_to);
            covered_up_to = b;
            n_pass_boxes++;
            continue;
        }

        double ylo = (double)a / YQ, yhi = (double)b / YQ;
        double t200_hi = 199.0 * yhi;
        int Dmax = 22, Rmax = 24, Kmax = 26, Hcut = 40;
        FullState st;
        int N_start;

        if (b <= 600000) {
            // Regime 1: initialize at N = 200 via exact piecewise-linear solver!
            int M_steps;
            if (a == 0) { M_steps = 11 * H_DEN; Dmax = 24; Rmax = 32; Kmax = 32; }
            else if (t200_hi <= 0.01) { M_steps = 9 * H_DEN; Dmax = 24; Rmax = 30; Kmax = 30; }
            else if (t200_hi <= 0.15) { M_steps = 6 * H_DEN; Dmax = 22; Rmax = 26; Kmax = 28; }
            else { M_steps = 4 * H_DEN; Dmax = 22; Rmax = 24; Kmax = 26; }
            auto U_hub = solve_static_hub_upper_int(Kmax);
            st = init_at_200(a, YQ, b, YQ, M_steps, Dmax, Rmax, Kmax, Hcut, U_hub);
            N_start = 200;
        } else {
            // Regime 2: ylo >= 6/199 ~ 0.03015; start directly at N_start where (N_start - 1)*ylo >= 22!
            // Choose N_start = max(m_max + 4, min(900, ceil(24 / ylo) + 2))
            int n_cut = (int)div_cl((i128)24 * YQ, a) + 2;
            N_start = max(m_max + 4, min(900, n_cut));
            Dmax = min(22, N_start - 1);
            Rmax = min(24, N_start - 1);
            Kmax = min(26, N_start - 2);
            Hcut = min(40, N_start - 2);
            auto U_hub = solve_static_hub_upper_int(Kmax);
            st.V.assign(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
            st.R.assign(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
            st.H.assign(Hcut + 2, 0);
            st.U.assign(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
            st.P.assign(Kmax + 2, 0);
            for (int d = 2; d <= Dmax; d++)
                for (int j = 1; j < d; j++)
                    st.V[d][j] = H_exact_N_lower(N_start, j, b, YQ);
            for (int p = 2; p <= Hcut + 1; p++) st.H[p] = st.V[2][1];
            for (int k = 2; k <= Kmax + 1; k++) {
                st.P[k] = B_exact_N_upper(N_start, 1, 1, k, a, YQ);
                int64_t cap = B_exact_N_upper(N_start, 0, 1, k, a, YQ);
                for (int j = 0; j < k; j++) st.U[k][j] = min(cap, U_hub[k][j]);
            }
        }

        auto U_hub = solve_static_hub_upper_int(Kmax);
        bool all_m_ok = true;
        int failed_m = -1;
        string fail_reason;

        // Sweep N backward from N_start - 1 down to 4!
        for (int N = N_start - 1; N >= 4; N--) {
            step_backward_N(st, N, a, YQ, b, YQ, Dmax, Rmax, Kmax, Hcut, U_hub);
            if (N <= m_max) {
                total_m_checks++;
                string fr;
                if (!verify_ledger_for_m(st, N, a, YQ, Rmax, Kmax, Hcut, fr)) {
                    all_m_ok = false;
                    failed_m = N;
                    fail_reason = fr;
                    break;
                }
            }
        }

        if (all_m_ok) {
            assert(a == covered_up_to);
            covered_up_to = b;
            n_pass_boxes++;
            cout << "  PASS y=[" << fixed << setprecision(6) << ylo << ", " << yhi
                 << "] (m=4.." << m_max << ", N_start=" << N_start << ")\n" << flush;
        } else {
            if (b - a <= 4) {
                n_fail_boxes++;
                cout << "  HARD FAIL y=[" << fixed << setprecision(6) << ylo << ", " << yhi
                     << "] at m=" << failed_m << " reason=" << fail_reason
                     << " M=" << (double)st.V[2][1] / S << " U20=" << (double)st.U[2][0] / S
                     << " H2=" << (double)st.H[2] / S << "\n" << flush;
            } else {
                int64_t mid = (a + b) / 2;
                boxes.push_back({mid, b});
                boxes.push_back({a, mid});
            }
        }
    }

    assert(n_fail_boxes > 0 || covered_up_to == 2 * ONE);
    cout << "\n=== SMALL-m CERTIFICATE SUMMARY (ALL m in {4..199} + Founder Seam b=2, ALL x in [0, 1)) ===\n";
    cout << "Boxes PASSED: " << n_pass_boxes << ", FAILED: " << n_fail_boxes
         << ", total (box, m) checks: " << total_m_checks
         << ", exact contiguous y-coverage: [0, " << covered_up_to << "/" << YQ << "]\n";
    return (n_fail_boxes == 0) ? 0 : 1;
}
