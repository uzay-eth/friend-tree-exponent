/*
 * cert_all_m.cpp — 100% Exact Outward-Rounded Integer (__int128) Certificate
 * Proves pre-placement conditional stochastic dominance:
 *   P(Z_{n+1} > x | T_n) > P(Z_{n+2} > x | T_n)
 * for EVERY finite tree T_n of size n >= m0 - 2 (m0 = 200, i.e. b = n+1 >= 199)
 * and EVERY threshold x in [0, 1), covering tau = (m-1)*log(1/(1-x)) in [0, 6]
 * via exact discrete-time piecewise-linear sub/superharmonics + normalized all-tree
 * signed ledger, and tau >= 6 via exact rational Poisson-Beta inequality.
 *
 * ZERO floating-point trust in any inequality! All arithmetic uses __int128
 * outward-rounded fixed-point intervals at internal scale S2 = 10^15 and grid scale S = 10^12.
 */

#include <bits/stdc++.h>
using namespace std;

using i128 = __int128_t;
using u128 = __uint128_t;

static constexpr int64_t S  = 1000000000000LL;       // 10^12 (grid probability scale)
static constexpr int64_t S2 = 1000000000000000LL;    // 10^15 (internal transcendental scale)

// Exact floor/ceil integer division for non-negative numerator, positive denominator
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

// Interval [lo, hi] representing real numbers in [lo/S2, hi/S2], always non-negative unless noted
struct Iv {
    i128 lo, hi;
};

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

// Exact outward enclosure of exp(x) for 0 <= x <= 16 (guaranteed overflow-free in i128!)
// Uses range reduction x = u * 2^k with u <= 1/8, 16-term Taylor series + exact remainder:
//   For 0 <= u <= 1/8: sum_{j=0}^{15} u^j/j! <= exp(u) <= sum_{j=0}^{15} u^j/j! + 2 * u^{16}/16!
static Iv iv_exp_pos(Iv x) {
    assert(x.lo >= 0 && x.hi <= (i128)16 * S2);
    int k = 0;
    Iv u = x;
    while (u.hi > S2 / 8) {
        u = {div_fl(u.lo, 2), div_cl(u.hi, 2)};
        k++;
    }
    Iv term = {S2, S2};
    Iv sum = {S2, S2};
    for (int j = 1; j <= 15; j++) {
        term = {div_fl(term.lo * u.lo, (i128)S2 * j), div_cl(term.hi * u.hi, (i128)S2 * j)};
        sum = iv_add(sum, term);
    }
    i128 rem = div_cl(term.hi * u.hi * 2, (i128)S2 * 16);
    sum.hi += rem;
    for (int i = 0; i < k; i++) {
        sum = iv_mul_nn(sum, sum);
    }
    return sum;
}

// Exact outward enclosure of exp(-x) for ANY x >= 0 via chunking x = q*4 + r (r in [0, 4]),
// completely overflow-free for all x >= 0!
static Iv iv_exp_neg(Iv x) {
    assert(x.lo >= 0);
    if (x.lo >= (i128)42 * S2) {
        return {0, 1};
    }
    static Iv em4 = {0, 0};
    if (em4.lo == 0) {
        Iv e4 = iv_exp_pos({(i128)4 * S2, (i128)4 * S2});
        em4 = {div_fl((i128)S2 * S2, e4.hi), div_cl((i128)S2 * S2, e4.lo)};
    }
    auto exp_neg_scalar = [&](i128 val, bool is_hi) -> i128 {
        if (val >= (i128)42 * S2) return is_hi ? 1 : 0;
        int chunks = (int)(val / ((i128)4 * S2));
        i128 rem = val - (i128)chunks * 4 * S2;
        Iv er = iv_exp_pos({rem, rem});
        Iv res = {div_fl((i128)S2 * S2, er.hi), div_cl((i128)S2 * S2, er.lo)};
        for (int c = 0; c < chunks; c++) {
            res = iv_mul_nn(res, em4);
        }
        return is_hi ? res.hi : res.lo;
    };
    // exp(-x) is decreasing in x: lo comes from x.hi, hi comes from x.lo!
    return {exp_neg_scalar(x.hi, false), exp_neg_scalar(x.lo, true)};
}

// Exact outward enclosure of log(y) for y >= 1 (y as Iv >= S2)
// Uses range reduction y = 2^k * z with z in [1, 2), then w = (z-1)/(z+1) in [0, 1/3]:
//   log(z) = 2 * sum_{j=0}^\infty w^{2j+1} / (2j+1)
// For w <= 1/3, tail after j=18 (power 39) <= 2 * w^{39}/39 * (9/8) < 10^{-18}!
static Iv iv_log_ge1_raw(Iv y) {
    assert(y.lo >= S2);
    // Compute log(2) once via w = 1/3
    auto log_unit = [](Iv z) -> Iv {
        // z in [1, 2], w = (z-1)/(z+1) in [0, 1/3]
        Iv num = {max((i128)0, z.lo - S2), z.hi - S2};
        Iv den = {z.lo + S2, z.hi + S2};
        Iv w = iv_div_pos(num, den);
        Iv w2 = iv_mul_nn(w, w);
        Iv poww = w;
        Iv s = {0, 0};
        for (int j = 0; j < 20; j++) {
            int d = 2 * j + 1;
            s.lo += div_fl(poww.lo, d);
            s.hi += div_cl(poww.hi, d);
            poww = iv_mul_nn(poww, w2);
        }
        // Tail <= poww.hi / 41 * 1/(1 - 1/9) = 9/8 * poww.hi / 41
        s.hi += div_cl(poww.hi * 9, 41 * 8);
        return {2 * s.lo, 2 * s.hi};
    };
    static Iv ln2 = {0, 0};
    if (ln2.lo == 0) {
        ln2 = log_unit({2 * S2, 2 * S2});
    }
    // Separate lo and hi range reductions for tight enclosure!
    auto log_scalar = [&](i128 val, bool is_hi) -> i128 {
        int k = 0;
        Iv z = {val, val};
        while (z.hi >= 2 * S2) {
            z = {div_fl(z.lo, 2), div_cl(z.hi, 2)};
            k++;
        }
        Iv res = iv_add({(i128)k * ln2.lo, (i128)k * ln2.hi}, log_unit(z));
        return is_hi ? res.hi : res.lo;
    };
    return {log_scalar(y.lo, false), log_scalar(y.hi, true)};
}

// Exact outward enclosure of y^r for y > 0 (Iv) and rational r = p/q in (0, 1)
static Iv iv_pow_rat(Iv y, int p, int q) {
    assert(y.lo > 0 && p > 0 && p < q);
    auto pow_scalar = [&](i128 val, bool is_hi) -> i128 {
        if (val == S2) return S2;
        if (val > S2) {
            Iv lg = iv_log_ge1_raw({val, val});
            Iv prod = {div_fl(lg.lo * p, q), div_cl(lg.hi * p, q)};
            Iv res = iv_exp_pos(prod);
            return is_hi ? res.hi : res.lo;
        } else {
            // val < S2 => inv = 1/val > 1, y^r = exp(-r * log(1/y))
            Iv inv = {div_fl((i128)S2 * S2, val), div_cl((i128)S2 * S2, val)};
            Iv lg = iv_log_ge1_raw(inv);
            Iv prod = {div_fl(lg.lo * p, q), div_cl(lg.hi * p, q)};
            Iv res = iv_exp_neg(prod);
            return is_hi ? res.hi : res.lo;
        }
    };
    return {pow_scalar(y.lo, false), pow_scalar(y.hi, true)};
}

// Exact outward enclosure of Gamma(1 + p/q) for 0 < p < q, then Gamma(p/q) = (q/p) * Gamma(1 + p/q)
// Method: split integral int_0^inf t^r e^{-t} dt into [0, 1] + [1, 24] + [24, inf).
//   On [0, 1]: exact alternating Taylor series of e^{-t} up to m=17 (even=upper, odd=lower)!
//   On [1, 24]: Riemann rectangle bounds with step 1/32 (736 subintervals):
//     Since t^r is increasing for r > 0, on [a, b]:
//       a^r * (e^{-a} - e^{-b}) <= int_a^b t^r e^{-t} dt <= b^r * (e^{-a} - e^{-b})!
//     (Even better: use width 1/64 -> error < 10^{-3}? Wait: 1/64 on t^r gives relative error ~ r/(64*t) ~ 0.005!
//      No, let's get 10^{-11} accuracy! Use Taylor/quadrature on [0, 12] via Gamma(1+r) = Gamma(N+1+r) / prod_{j=1}^N (j+r)
//      with Stirling's EXACT asymptotic enclosure for Gamma(N+1+r) at N = 20!!
//      WAIT A SECOND! Check Stirling for x = N + 1 + r >= 21:
//      Exact Binet/Stirling bounds (Robbins / Nemes / Binet's first formula, 100% elementary!):
//        For ALL real x > 0:
//          0 < log Gamma(x) - (x - 1/2)log(x) + x - (1/2)log(2*pi) - 1/(12x) + 1/(360 x^3) < 1/(1260 x^5)!
//        At x >= 21: 1/(1260 * 21^5) < 2 * 10^{-10}!!
//        Even simpler: at x >= 10, Binet's exact bound:
//          1/(12x) - 1/(360 x^3) < log Gamma(x) - ((x - 1/2)log x - x + (1/2)log(2pi)) < 1/(12x) !!
//        Wait: at x = 60: 1/(360 * 60^3) = 1.28 * 10^{-11}!!
//        AND we ONLY need a 1-line check of pi (e.g. 3141592653589793 / 10^15 < pi < 3141592653589794 / 10^15)
//        and Binet's elementary integral:
//          J(x) = int_0^inf (1/(e^t - 1) - 1/t + 1/2) * (e^{-xt}/t) dt,
//        where for all t > 0, t/12 - t^3/720 < 1/(e^t - 1) - 1/t + 1/2 < t/12 (elementary Bernoulli inequality!),
//        so 1/(12x) - 1/(360 x^3) < J(x) < 1/(12x) FOR ALL x > 0!!
//      OH! THAT IS SO CLEAN, SO ELEMENTARY, AND GIVES 11-DIGIT ACCURACY IN 50 OPERATIONS!!)
static map<pair<int,int>, Iv> gamma_cache;
static Iv iv_gamma_rat(int p, int q) {
    assert(p > 0 && p < q);
    int g = std::gcd(p, q);
    p /= g; q /= g;
    auto key = make_pair(p, q);
    auto it = gamma_cache.find(key);
    if (it != gamma_cache.end()) return it->second;

    // Shift to x = N + p/q with N = 60, so x in (60, 61)
    const int N = 60;
    // x = (N*q + p) / q
    i128 x_num = (i128)N * q + p;
    Iv x = iv_from_rat(x_num, q);
    Iv lnx = iv_log_ge1_raw(x);
    // (x - 1/2)*ln(x) - x
    Iv xmhalf = iv_sub(x, {S2 / 2, S2 / 2});
    Iv term1 = iv_mul_nn(xmhalf, lnx);
    Iv base = iv_sub(term1, x);
    // (1/2) * log(2*pi): 2*pi in [6283185307179586, 6283185307179587] / 10^15
    Iv two_pi = {6283185307179586LL, 6283185307179587LL};
    Iv ln2pi = iv_log_ge1_raw(two_pi);
    Iv half_ln2pi = {div_fl(ln2pi.lo, 2), div_cl(ln2pi.hi, 2)};
    // J(x) in [1/(12x) - 1/(360 x^3), 1/(12x)]
    i128 j_hi = div_cl((i128)S2 * S2, 12 * x.lo);
    i128 x3_lo = div_fl(div_fl(x.lo * x.lo, S2) * x.lo, S2);
    i128 j_lo = div_fl((i128)S2 * S2, 12 * x.hi) - div_cl((i128)S2 * S2, 360 * x3_lo);
    Iv lnG = iv_add(iv_add(base, half_ln2pi), {j_lo, j_hi});
    // Now divide Gamma(N + p/q) by prod_{k=0}^{N-1} (k + p/q) IN LOG SPACE to avoid any overflow!
    // Sum log(k + p/q) for k = 1..N-1 (all >= 1!), plus log(p/q) = -log(q/p)!
    Iv sum_ln = {0, 0};
    for (int k = 1; k < N; k++) {
        Iv v = iv_from_rat((i128)k * q + p, q);
        sum_ln = iv_add(sum_ln, iv_log_ge1_raw(v));
    }
    Iv ln_qp = iv_log_ge1_raw(iv_from_rat(q, p));
    // log Gamma(p/q) = lnG - sum_ln + ln(q/p)
    Iv ln_res = iv_add(iv_sub(lnG, sum_ln), ln_qp);
    // Note Gamma(p/q) >= Gamma(1) = 1 for all p/q in (0, 1) (since Gamma(1)=1 and Gamma'(x)<0 on (0,1])!
    // Actually Gamma(r) > 1 for all r in (0, 1), so ln_res.lo >= 0!
    assert(ln_res.lo >= 0);
    Iv res = iv_exp_pos(ln_res);
    gamma_cache[key] = res;
    return res;
}

// Exact outward UPPER bound on e^v * Q(a, v) = e^v * P(Gamma(a, 1) > v) for rational a = p/q in (0, 1) and v_iv > 0!
// Since e^v * Q(a, v) = (1/Gamma(a)) int_0^inf (v + t)^{a-1} e^{-t} dt is STRICTLY DECREASING in v for a in (0, 1),
// we evaluate at v = v_iv.lo (lower bound on v -> UPPER bound on e^v Q(a, v))!
// Two regimes:
//   1. v <= 2: Q(a, v) = 1 - gamma(a, v)/Gamma(a), where gamma(a, v) >= v^a * sum_{n=0}^{2K-1} (-1)^n v^n / (n!(a+n)) (odd truncation!)
//      Multiply upper bound on Q(a, v) by exp(v)_hi!
//   2. v > 2: Stieltjes continued fraction C_{2K+1} (with K = 45) for e^v * Gamma(a, v), divided by Gamma(a)_lo!
static i128 exp_gamma_upper_hi_frac(int p, int q, i128 v_lo) {
    assert(p > 0 && p < q && v_lo > 0);
    Iv gam = iv_gamma_rat(p, q);
    Iv v_iv = {v_lo, v_lo};
    Iv va = iv_pow_rat(v_iv, p, q);

    if (v_lo <= 2 * S2) {
        // Alternating Taylor series for gamma(a, v) / v^a = sum_{n=0}^\infty (-1)^n v^n / (n! * (p/q + n))
        // Stop at odd index 2K - 1 = 23 (since v <= 2, 2^23 / 23! < 3e-16 << 1/S2!)
        i128 s_lo = 0;
        // term_n = v^n / n! (track both lo and hi!)
        Iv pow_vn = {S2, S2}; // n = 0
        for (int n = 0; n <= 23; n++) {
            // Add (-1)^n * pow_vn * q / (p + n*q)
            i128 den = (i128)p + (i128)n * q;
            if (n % 2 == 0) {
                // positive term -> round DOWN
                s_lo += div_fl(pow_vn.lo * q, den);
            } else {
                // negative term -> round UP magnitude (subtract ceil)
                s_lo -= div_cl(pow_vn.hi * q, den);
            }
            if (n < 23) {
                pow_vn = {
                    div_fl(pow_vn.lo * v_lo, (i128)S2 * (n + 1)),
                    div_cl(pow_vn.hi * v_lo, (i128)S2 * (n + 1))
                };
            }
        }
        if (s_lo < 0) s_lo = 0;
        i128 gam_inc_lo = div_fl(va.lo * s_lo, S2);
        // Q(a, v) <= 1 - gam_inc_lo / gam.hi
        i128 Q_hi = S2 - div_fl(gam_inc_lo * S2, gam.hi);
        if (Q_hi < 0) Q_hi = 0;
        Iv ev = iv_exp_pos(v_iv);
        return div_cl(ev.hi * Q_hi, S2);
    } else {
        // Stieltjes continued fraction C_{2K+1} with K = 45 (odd convergent -> STRICT UPPER BOUND on e^v * Gamma(a, v)!)
        // a_{2n} = n - p/q, b_{2n} = 1; a_{2n+1} = n, b_{2n+1} = v.
        // Start at D_{2K+1}_lo = v_lo (odd -> lower bound!)
        const int K = 45;
        i128 D_odd_lo = v_lo;
        for (int n = K; n >= 1; n--) {
            // j = 2n (even -> UPPER bound): D_{2n}_hi = 1 + n / D_{2n+1}_lo
            i128 D_even_hi = S2 + div_cl((i128)n * S2 * S2, D_odd_lo);
            // j = 2n - 1 (odd -> LOWER bound): D_{2n-1}_lo = v + (n - p/q) / D_{2n}_hi
            // (n - p/q) = (n*q - p) / q
            i128 a_even_num = (i128)n * q - p;
            D_odd_lo = v_lo + div_fl(a_even_num * S2 * S2, (i128)q * D_even_hi);
        }
        // e^v * Gamma(a, v) <= va.hi / D_1_lo
        i128 ev_Gam_hi = div_cl(va.hi * S2, D_odd_lo);
        // Divide by gam.lo
        return div_cl(ev_Gam_hi * S2, gam.lo);
    }
}

// Exact outward UPPER bound on e^v * Q(r, v) for r = int_part + p/q, where int_part in {0, 1, 2} and 0 <= p < q!
// Note:
//   - For r < 1 (int_part == 0, p > 0): e^v * Q(r, v) is DECREASING in v -> use v_lo!
//   - For r == 1 (int_part == 1, p == 0): e^v * Q(1, v) == 1 identically!
//   - For r > 1: e^v * Q(r, v) is INCREASING in v -> use v_hi!
//     And by exact Gamma recurrence:
//       e^v Q(1 + a, v) = e^v Q(a, v) + v^a / Gamma(1 + a)
//       e^v Q(2 + a, v) = e^v Q(a, v) + v^a / Gamma(1 + a) + v^{1+a} / Gamma(2 + a)
//     (Wait: for e^v Q(a, v) inside that sum at v = v_hi, evaluate exp_gamma_upper_hi_frac(p, q, v_hi)!)
static i128 exp_gamma_upper_hi_general(int int_part, int p, int q, i128 v_lo, i128 v_hi) {
    if (int_part == 0) {
        assert(p > 0 && p < q);
        return exp_gamma_upper_hi_frac(p, q, v_lo);
    }
    if (p == 0) {
        if (int_part == 1) return S2; // exact 1
        if (int_part == 2) return S2 + v_hi; // exact 1 + v
        assert(false);
    }
    // int_part is 1 or 2, 0 < p < q: evaluate at v = v_hi!
    // Check: is e^v Q(1+a, v) strictly increasing in v?
    // d/dv [int_0^inf (v+t)^a e^{-t} dt] = a int_0^inf (v+t)^{a-1} e^{-t} dt > 0! YES!
    i128 base = exp_gamma_upper_hi_frac(p, q, v_hi);
    Iv gam_a = iv_gamma_rat(p, q);
    // Gamma(1 + a) = a * Gamma(a) = (p/q) * gam_a
    i128 gam1_lo = div_fl(gam_a.lo * p, q);
    Iv va = iv_pow_rat({v_hi, v_hi}, p, q);
    i128 term1 = div_cl(va.hi * S2, gam1_lo);
    if (int_part == 1) return base + term1;
    // int_part == 2: add v^{1+a} / Gamma(2+a) = (v / (1+a)) * (v^a / Gamma(1+a))
    // Gamma(2+a) = (1 + p/q) * Gamma(1+a) = (p+q)/q * gam1
    i128 gam2_lo = div_fl(gam1_lo * (p + q), q);
    i128 v1a_hi = div_cl(va.hi * v_hi, S2);
    i128 term2 = div_cl(v1a_hi * S2, gam2_lo);
    return base + term1 + term2;
}

// --- Exact outward-rounded grid primitives (scale S = 10^12) ---

static constexpr int H_DEN = 128; // h = 1/128

// Precomputed exp(i / 128) for i = 0..1600
static vector<Iv> EXP_GRID;
static void init_exp_grid(int max_i) {
    EXP_GRID.resize(max_i + 1);
    for (int i = 0; i <= max_i; i++) {
        EXP_GRID[i] = iv_exp_pos(iv_from_rat(i, H_DEN));
    }
}

// Exact outward LOWER bound (scale S) on RAW H_{N, j}(x) for all m >= m0, tau <= tau_hi (rational p_tau/q_tau),
// and log(N/m) <= t_idx / H_DEN!
static int64_t H_raw_lower_int(int j, int t_idx, int64_t tau_num, int64_t tau_den, int m0) {
    if (j <= 0) return 0;
    Iv inv_e = iv_exp_neg({S2, S2});
    if (tau_num == 0) {
        // e^{-1} * sum_{k=0}^{j-1} 1/k!
        i128 sum = S2;
        i128 fk = S2;
        for (int k = 1; k < j; k++) {
            fk = div_fl(fk, k);
            sum += fk;
        }
        i128 ans = div_fl(inv_e.lo * sum, S2);
        return (int64_t)min((i128)S, div_fl(ans * S, S2));
    }
    // a_plus = (m0 * exp(t_hi) - 1) / (m0 - 1)
    Iv et = EXP_GRID[t_idx];
    // For a_plus: we need UPPER bound in exponent L_exp = tau_hi * (a_plus + 1/(m0-1)),
    // and LOWER bound in polynomial factor f_num = a_plus - (i-1)/(m0-1)!
    // WAIT! Be careful! Both use the SAME fixed real number a_plus = (m0 * exp(t_hi) - 1) / (m0 - 1):
    //   a_plus_lo <= a_plus <= a_plus_hi
    i128 a_plus_lo = div_fl((i128)m0 * et.lo - S2, m0 - 1);
    i128 a_plus_hi = div_cl((i128)m0 * et.hi - S2, m0 - 1);
    // L_exp_hi = tau_hi * (a_plus_hi + 1/(m0 - 1))
    i128 a_plus_shift_hi = a_plus_hi + div_cl(S2, m0 - 1);
    i128 L_exp_hi = div_cl(a_plus_shift_hi * tau_num, tau_den);
    if (L_exp_hi >= (i128)42 * S2) return 0;
    // Split L_exp_hi = chunks * 4 + rem so exp(rem) and exp(4) never exceed 55 * S2 (overflow-free!)
    int exp_chunks = (int)(L_exp_hi / ((i128)4 * S2));
    i128 exp_rem = L_exp_hi - (i128)exp_chunks * 4 * S2;
    static i128 e4_hi = 0;
    if (e4_hi == 0) e4_hi = iv_exp_pos({(i128)4 * S2, (i128)4 * S2}).hi;
    i128 erem_hi = iv_exp_pos({exp_rem, exp_rem}).hi;
    // Clean overflow-free exact __int128 loop with dynamic power-of-2 scaling!
    // Represent values as val * 2^shift / S2.
    int shift = 0;
    i128 cur = S2;
    i128 B_val[32];
    int B_sh[32];
    i128 sum_val = S2;
    B_val[1] = sum_val; B_sh[1] = 0;
    for (int i = 1; i < j; i++) {
        i128 f_num_lo = a_plus_lo - div_cl((i128)(i - 1) * S2, m0 - 1);
        if (f_num_lo <= 0) {
            cur = 0;
        } else {
            i128 fac_lo = div_fl(f_num_lo * tau_num, (i128)tau_den * i);
            while (cur > ((i128)1 << 62) && fac_lo > 0) {
                cur >>= 1;
                sum_val >>= 1;
                shift++;
            }
            cur = div_fl(cur * fac_lo, S2);
            sum_val += cur;
            while (sum_val > ((i128)1 << 62)) {
                cur >>= 1;
                sum_val >>= 1;
                shift++;
            }
        }
        B_val[i + 1] = sum_val;
        B_sh[i + 1] = shift;
    }
    // Divide each B_val[r] * 2^{B_sh[r]} by erem_hi * (e4_hi)^exp_chunks, consuming shift as we go!
    i128 B_lo[32];
    for (int r = 1; r <= j; r++) {
        i128 v = B_val[r];
        int sh = B_sh[r];
        // First divide by erem_hi: v = floor(v * S2 / erem_hi)
        v = div_fl(v * S2, erem_hi);
        for (int c = 0; c < exp_chunks; c++) {
            // Pull up to 6 bits of shift before dividing by e4_hi (~54.6 * S2 < 64 * S2) to preserve precision!
            int take = min(sh, 6);
            v <<= take;
            sh -= take;
            v = div_fl(v * S2, e4_hi);
        }
        while (sh > 0 && v < S2) {
            v <<= 1;
            sh--;
        }
        if (sh > 0 || v > S2) v = S2;
        B_lo[r] = v;
    }
    i128 ans = 0;
    i128 fk = S2;
    for (int k = 0; k < j; k++) {
        if (k > 0) fk = div_fl(fk, k);
        ans += div_fl(fk * B_lo[j - k], S2);
    }
    ans = div_fl(inv_e.lo * ans, S2);
    return (int64_t)min((i128)S, div_fl(ans * S, S2));
}

// Exact outward UPPER bound (scale S) on RAW B_{N, r}(x) for all m >= m0, tau >= tau_lo (p/q),
// and log((N-1)/(m-1)) >= t_idx / H_DEN!
// Here r = int_part + p/q (with int_part in {0, 1, 2}, 0 <= p < q).
static int64_t B_raw_upper_int(int int_part, int p, int q, int t_idx, int64_t tau_num, int64_t tau_den, int m0) {
    if (tau_num <= 0) return S;
    Iv et = EXP_GRID[t_idx];
    // L_lo = tau_lo * exp(t_lo): lower bound on exponent in q^{N-1} = exp(-tau * a)
    assert(tau_num * 2 <= tau_den * (m0 - 1)); // ensures u0 <= 1/2 so w_hi <= 1.3 * L_lo and e^{-42} tail < 1/S
    i128 L_lo_fl = div_fl(et.lo * tau_num, tau_den);
    i128 L_lo_cl = div_cl(et.hi * tau_num, tau_den);
    if (L_lo_fl >= (i128)42 * S2) return 1; // <= exp(-42) * 1.1 * (1 + 55 + 55^2/2) < 10^{-14} < 1/S
    i128 exp_q_hi = iv_exp_neg({L_lo_fl, L_lo_fl}).hi;

    // B0 = (m0 - 1) * a_lo + 1 - r >= (m0 - 1)*et.lo + S2 - r*S2
    i128 r_S2 = (i128)int_part * S2 + div_cl((i128)p * S2, q);
    i128 B0_lo = (i128)(m0 - 1) * et.lo + S2 - r_S2;
    assert(B0_lo > 2 * S2);

    if (int_part == 0) {
        // 0 < r < 1: Vcut = min(24, B0/2) (take integer Vcut = min(24, floor(B0/(2*S2))))
        int Vcut = (int)min((i128)24, B0_lo / (2 * S2));
        // alpha_lo = 1 - Vcut / (2 * B0_lo)
        i128 alpha_lo = S2 - div_cl((i128)Vcut * S2 * S2, 2 * B0_lo);
        // w_eff_lo = L_lo_fl / alpha_hi? Wait: we want LOWER bound on w_eff = w_lo / alpha,
        // where alpha <= 1, and actually (w_lo + alpha*v)^{r-1} = alpha^{r-1} (w_lo/alpha + v)^{r-1},
        // which is DECREASING in alpha for fixed w_lo! So replacing alpha by alpha_lo <= alpha is an UPPER bound!
        // Specifically: for alpha_lo <= alpha, (w_lo + alpha*v)^{r-1} <= (w_lo + alpha_lo*v)^{r-1} = alpha_lo^{r-1}*(w_lo/alpha_lo + v)^{r-1}!
        // So we use w_eff_lo = floor(L_lo_fl * S2 / alpha_lo) >= L_lo_fl!
        i128 w_eff_lo = div_fl(L_lo_fl * S2, alpha_lo);
        // alpha_lo^{r-1} = (1 / alpha_lo)^{1 - r} = (1 / alpha_lo)^{(q-p)/q}
        i128 inv_alpha_hi = div_cl((i128)S2 * S2, alpha_lo);
        i128 alpha_pow_hi = iv_pow_rat({inv_alpha_hi, inv_alpha_hi}, q - p, q).hi;
        i128 main_part = div_cl(alpha_pow_hi * exp_gamma_upper_hi_frac(p, q, w_eff_lo), S2);

        // Tail part: (w_lo + B0*(1 - exp(-Vcut/B0)))^{r-1} * exp(-Vcut) / Gamma(r)
        // Lower bound on B0*(1 - exp(-Vcut/B0)) >= B0_lo * (1 - exp(-Vcut/B0_lo)_hi) >= alpha_lo * Vcut!
        // Wait: check B0*(1 - e^{-V/B0}) >= V*(1 - V/(2 B0)) = alpha_lo * V! YES! (Since 1 - e^{-u} >= u - u^2/2 for u >= 0!)
        i128 base_tail_lo = L_lo_fl + div_fl(alpha_lo * Vcut, 1);
        // base_tail_lo^{r-1} = (1/base_tail_lo)^{1-r}
        i128 inv_bt_hi = div_cl((i128)S2 * S2, base_tail_lo);
        i128 bt_pow_hi = iv_pow_rat({inv_bt_hi, inv_bt_hi}, q - p, q).hi;
        i128 evcut_hi = iv_exp_neg({(i128)Vcut * S2, (i128)Vcut * S2}).hi;
        Iv gam = iv_gamma_rat(p, q);
        i128 tail_part = div_cl(div_cl(bt_pow_hi * evcut_hi, S2) * S2, gam.lo);

        i128 ans = div_cl(exp_q_hi * (main_part + tail_part), S2);
        return (int64_t)min((i128)S, div_cl(ans * S, S2));
    } else {
        // r >= 1: need w_hi = L_lo * (e^{u0} - 1) / u0 where u0 = tau_lo / (m0 - 1)
        Iv u0 = iv_from_rat(tau_num, tau_den * (m0 - 1));
        i128 ex_ratio_hi;
        if (u0.hi == 0) {
            ex_ratio_hi = S2;
        } else {
            Iv eu0 = iv_exp_pos(u0);
            // (e^u - 1)/u is increasing in u, bound above using u0.hi:
            // Wait: if u0.lo == 0 (impossible since tau_num > 0 and tau_den*(m0-1) < 10^12), use Taylor:
            // (e^u - 1)/u <= 1 + u/2 + u^2/6 + u^3/24 * 2 (for u <= 1/2)!
            if (u0.hi <= S2 / 2) {
                i128 u = u0.hi;
                i128 u2 = div_cl(u * u, S2);
                i128 u3 = div_cl(u2 * u, S2);
                ex_ratio_hi = S2 + div_cl(u, 2) + div_cl(u2, 6) + div_cl(u3, 12);
            } else {
                ex_ratio_hi = div_cl((eu0.hi - S2) * S2, u0.lo);
            }
        }
        i128 w_hi = div_cl(L_lo_cl * ex_ratio_hi, S2);
        i128 pref_hi;
        if (int_part == 1) {
            // 1 + (p/q) / B0
            pref_hi = S2 + div_cl((i128)p * S2 * S2, (i128)q * B0_lo);
        } else {
            // (1 + (1 + p/q)/B0) * (1 + (p/q)/B0)
            i128 f1 = S2 + div_cl((i128)(p + q) * S2 * S2, (i128)q * B0_lo);
            i128 f2 = S2 + div_cl((i128)p * S2 * S2, (i128)q * B0_lo);
            pref_hi = div_cl(f1 * f2, S2);
        }
        i128 val = div_cl(pref_hi * exp_gamma_upper_hi_general(int_part, p, q, L_lo_fl, w_hi), S2);
        i128 ans = div_cl(exp_q_hi * val, S2);
        return (int64_t)min((i128)S, div_cl(ans * S, S2));
    }
}

// Static time-independent hubhood upper bound U_hub[k][j] (exact outward-rounded integer at scale S!)
// Uses analytic supermartingale cap 1/(j+1) when k >= j+2 (permanent non-leaf neighbor) and min(1, 2/(j+1)) when j == k-1.
static vector<vector<int64_t>> solve_static_hub_upper_int(int Kmax) {
    vector<vector<int64_t>> U(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    for (int j = 0; j < Kmax; j++) {
        U[Kmax + 1][j] = (int64_t)div_cl(S, j + 1);
    }
    U[Kmax + 1][Kmax] = min(S, (int64_t)div_cl(2 * S, Kmax + 1));
    for (int k = Kmax; k >= 2; k--) {
        for (int j = 0; j < k; j++) {
            int64_t loss_pay = (j > 0) ? U[k][j - 1] : 0;
            int64_t gain_pay = U[k + 1][j + 1];
            int64_t best = 0;
            for (int g2 : {2 * (j + 1), k + j + 1}) {
                i128 num = (i128)2 * S + (i128)2 * j * loss_pay + (i128)k * g2 * gain_pay;
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

struct BoxResult {
    int64_t M_lo;
    vector<int64_t> H_host; // index 2..Hcut+1
    vector<int64_t> L_lo;   // index 2..Rmax
    vector<vector<int64_t>> U_hi; // [k][j], k=2..Kmax+1
    int64_t U_tail_cap_0;   // B_raw_upper(1/(Kmax+2), 0)
};

// Exact outward-rounded discrete-time piecewise-linear solver for ALL m >= m0 on [tau_lo, tau_hi]!
// Uses the exact Affine Knot-Crossing Lemma in DAG order for BOTH lower submartingale and upper supermartingale.
static BoxResult solve_box_int(int64_t tlo_n, int64_t tlo_d, int64_t thi_n, int64_t thi_d,
                               int m0, int M_steps, int Dmax, int Rmax, int Kmax, int Hcut) {
    assert(m0 - 1 >= H_DEN); // guarantees Delta_t <= 1/(m0-1) <= h = 1/128!

    // --- 1. LOWER SUBMARTINGALE at tau = tau_hi ---
    // Continuous piecewise-linear on [0, T + h], constant on [T, T + h] (gamma_M == 0).
    // Terminal payoff at t >= T is H_raw_lower(j, M_steps + 1), valid for all t' in [T, T + h].
    vector<vector<int64_t>> Vn(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
    vector<vector<int64_t>> Rn(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
    vector<int64_t> Hn(Hcut + 2, 0);
    vector<vector<int64_t>> gVn(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
    vector<vector<int64_t>> gRn(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
    vector<int64_t> gHn(Hcut + 2, 0);
    vector<vector<bool>> stopVn(Dmax + 2, vector<bool>(Dmax + 2, true));
    vector<bool> stopHn(Hcut + 2, true);

    for (int d = 2; d <= Dmax; d++) {
        for (int j = 1; j < d; j++) {
            Vn[d][j] = H_raw_lower_int(j, M_steps + 1, thi_n, thi_d, m0);
        }
    }
    for (int p = 2; p <= Hcut + 1; p++) Hn[p] = Vn[2][1];

    for (int i = M_steps - 1; i >= 0; i--) {
        int t_next = i + 1;
        // Exact upper bound on Delta_t when crossing knot t_{i+1}:
        // t' > t_{i+1} => N + 1 > m0 * exp(t_{i+1}) => N > m0 * exp(t_{i+1}) - 1 => dt < 1 / (m0*exp(t_{i+1}) - 1).
        i128 denom_cor = (i128)m0 * EXP_GRID[t_next].lo - S2;
        assert(denom_cor > 0);
        i128 dt_cor_hi = div_cl((i128)S2 * S2, denom_cor);

        vector<vector<int64_t>> Vc(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
        vector<vector<int64_t>> Rc(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
        vector<int64_t> Hc(Hcut + 2, 0);
        vector<vector<int64_t>> gVc(Dmax + 2, vector<int64_t>(Dmax + 2, 0));
        vector<vector<int64_t>> gRc(Dmax + 2, vector<int64_t>(Rmax + 2, 0));
        vector<int64_t> gHc(Hcut + 2, 0);
        vector<vector<bool>> stopVc(Dmax + 2, vector<bool>(Dmax + 2, false));
        vector<bool> stopHc(Hcut + 2, false);

        vector<int64_t> H_floor(Dmax, 0);
        for (int j = 1; j < Dmax; j++) {
            H_floor[j] = H_raw_lower_int(j, t_next, thi_n, thi_d, m0);
        }
        // Boundary d = Dmax: always stopped! Continuous linear from H_floor(t_{i+1}) down to Vn[Dmax][j].
        for (int j = 1; j < Dmax; j++) {
            assert(H_floor[j] >= Vn[Dmax][j]);
            Vc[Dmax][j] = H_floor[j];
            stopVc[Dmax][j] = true;
            gVc[Dmax][j] = (int64_t)((i128)(Vc[Dmax][j] - Vn[Dmax][j]) * H_DEN);
        }

        // Exact DAG topological order: d descending from Dmax-1 down to 1!
        for (int d = Dmax - 1; d >= 1; d--) {
            // 1a. Return states R[d][k] for k = Rmax down to 2 (targets V[d+1][1] and R[d][k+1] ALREADY computed!)
            for (int k = Rmax; k >= 2; k--) {
                i128 kd = (i128)k * d;
                i128 r1 = d;
                i128 r2 = (i128)k * (d * (k - 1) + 1);
                i128 rtot = r1 + r2;
                // (L-R) at t_{i+1}: gamma <= (r1*Vn[d+1][1] + r2*Rn[d][k+1] - rtot*Rn[d][k]) / kd
                i128 num_R = r1 * Vn[d + 1][1] + r2 * Rn[d][k + 1] - rtot * Rn[d][k];
                i128 g_R = (num_R > 0) ? div_fl(num_R, kd) : 0;
                // (L-L) at t_i: gamma <= (r1*Vc[d+1][1] + r2*Rc[d][k+1] - rtot*Rn[d][k]) * H_DEN / (kd*H_DEN + rtot)
                i128 num_L = (r1 * Vc[d + 1][1] + r2 * Rc[d][k + 1] - rtot * Rn[d][k]) * H_DEN;
                i128 den_L = kd * H_DEN + rtot;
                i128 g_L = (num_L > 0) ? div_fl(num_L, den_L) : 0;
                // Note: at i == M_steps - 1, Rn == 0, gRn == 0, gVn == 0, Vn >= 0 => (L-KNOT) holds trivially!
                // For i < M_steps - 1, R is never stopped on slab i+1!
                i128 gamma = min(g_R, g_L);
                if (gamma < 0) gamma = 0;
                int64_t delta = (int64_t)div_fl(gamma, H_DEN);
                Rc[d][k] = Rn[d][k] + delta;
                gRc[d][k] = (int64_t)((i128)delta * H_DEN);
            }

            // 1b. Positive stock V[d][j] for j = 1 up to d-1 (targets V[d+1][j+1] and V[d][j-1]/R[d][2] ALREADY computed!)
            if (d >= 2) {
                for (int j = 1; j < d; j++) {
                    int64_t lt_n = (j > 1) ? Vn[d][j - 1] : Rn[d][2];
                    int64_t lt_c = (j > 1) ? Vc[d][j - 1] : Rc[d][2];
                    int64_t glt_n = (j > 1) ? gVn[d][j - 1] : gRn[d][2];
                    int64_t gt_n = Vn[d + 1][j + 1];
                    int64_t gt_c = Vc[d + 1][j + 1];
                    int64_t ggt_n = gVn[d + 1][j + 1];

                    i128 rg = (i128)d * j;
                    i128 rl = j;
                    i128 rtot = rg + rl;

                    // Check (L-KNOT) if (i+1, d, j) was stopped!
                    bool can_continue = true;
                    if (stopVn[d][j]) {
                        // Need -d*gVn[d][j] + rg*(gt_n - Vn) + rl*(lt_n - Vn) - delta2*(rg*(ggt_n - gVn) + rl*(glt_n - gVn)) >= 0
                        // for delta2 in {0, dt_cor_hi}!
                        i128 base_knot = -(i128)d * gVn[d][j] + rg * (gt_n - Vn[d][j]) + rl * (lt_n - Vn[d][j]);
                        i128 slope_diff = rg * (ggt_n - gVn[d][j]) + rl * (glt_n - gVn[d][j]);
                        i128 pen = (slope_diff > 0) ? div_cl(slope_diff * dt_cor_hi, S2) : 0;
                        if (base_knot - pen < 0) can_continue = false;
                    }

                    int64_t v_ode = -1;
                    if (can_continue) {
                        i128 num_R = rg * gt_n + rl * lt_n - rtot * Vn[d][j];
                        i128 g_R = (num_R > 0) ? div_fl(num_R, d) : 0;
                        i128 num_L = (rg * gt_c + rl * lt_c - rtot * Vn[d][j]) * H_DEN;
                        i128 den_L = (i128)d * H_DEN + rtot;
                        i128 g_L = (num_L > 0) ? div_fl(num_L, den_L) : 0;
                        // Wait: if num_R < 0, continuing violates (L-R) even at gamma = 0!
                        if (num_R >= 0 && num_L >= 0) {
                            i128 gamma = min(g_R, g_L);
                            int64_t delta = (int64_t)div_fl(gamma, H_DEN);
                            v_ode = Vn[d][j] + delta;
                        } else {
                            can_continue = false;
                        }
                    }

                    if (!can_continue || H_floor[j] > v_ode) {
                        assert(H_floor[j] >= Vn[d][j]);
                        Vc[d][j] = H_floor[j];
                        stopVc[d][j] = true;
                        gVc[d][j] = (int64_t)((i128)(Vc[d][j] - Vn[d][j]) * H_DEN);
                    } else {
                        Vc[d][j] = v_ode;
                        stopVc[d][j] = false;
                        gVc[d][j] = (int64_t)((i128)(Vc[d][j] - Vn[d][j]) * H_DEN);
                    }
                }
            }
        }

        // 1c. Serial host bonus H[p] for p = Hcut down to 2 (targets V[3][2], R[2][2], H[p+1] ALREADY computed!)
        Hc[Hcut + 1] = Vc[2][1];
        gHc[Hcut + 1] = gVc[2][1];
        stopHc[Hcut + 1] = stopVc[2][1];
        for (int p = Hcut; p >= 2; p--) {
            i128 r_gain = 2 * (p + 1);
            i128 r_loss = p;
            bool can_continue = true;
            i128 best_gamma = (i128)1 << 120;
            for (i128 r_host : {(i128)0, (i128)p * (2 * p - 1)}) {
                i128 rtot = r_gain + r_loss + r_host;
                if (stopHn[p]) {
                    i128 base_knot = -(i128)2 * p * gHn[p]
                                     + r_gain * (Vn[3][2] - Hn[p])
                                     + r_loss * (Rn[2][2] - Hn[p])
                                     + r_host * (Hn[p + 1] - Hn[p]);
                    i128 slope_diff = r_gain * (gVn[3][2] - gHn[p])
                                      + r_loss * (gRn[2][2] - gHn[p])
                                      + r_host * (gHn[p + 1] - gHn[p]);
                    i128 pen = (slope_diff > 0) ? div_cl(slope_diff * dt_cor_hi, S2) : 0;
                    if (base_knot - pen < 0) can_continue = false;
                }
                i128 num_R = r_gain * Vn[3][2] + r_loss * Rn[2][2] + r_host * Hn[p + 1] - rtot * Hn[p];
                i128 num_L = (r_gain * Vc[3][2] + r_loss * Rc[2][2] + r_host * Hc[p + 1] - rtot * Hn[p]) * H_DEN;
                if (num_R < 0 || num_L < 0) {
                    can_continue = false;
                } else {
                    i128 g_R = div_fl(num_R, 2 * p);
                    i128 den_L = (i128)2 * p * H_DEN + rtot;
                    i128 g_L = div_fl(num_L, den_L);
                    i128 gc = min(g_R, g_L);
                    if (gc < best_gamma) best_gamma = gc;
                }
            }
            int64_t h_ode = -1;
            if (can_continue) {
                int64_t delta = (int64_t)div_fl(best_gamma, H_DEN);
                h_ode = Hn[p] + delta;
            }
            if (!can_continue || H_floor[1] > h_ode) {
                assert(H_floor[1] >= Hn[p]);
                Hc[p] = H_floor[1];
                stopHc[p] = true;
                gHc[p] = (int64_t)((i128)(Hc[p] - Hn[p]) * H_DEN);
            } else {
                Hc[p] = h_ode;
                stopHc[p] = false;
                gHc[p] = (int64_t)((i128)(Hc[p] - Hn[p]) * H_DEN);
            }
        }

        // Strict knot monotonicity assertions on slab i (implies all real t in [t_i, t_{i+1}] by linearity!):
        for (int d = 2; d < Dmax; d++) {
            for (int j = 1; j < d; j++) {
                assert(Vc[d + 1][j + 1] >= Vc[d][j]);
            }
        }
        for (int d = 1; d < Dmax; d++) {
            for (int k = 2; k < Rmax; k++) assert(Rc[d][k] >= Rc[d][k + 1]);
            assert(Vc[d + 1][1] >= Rc[d][2]);
        }
        for (int p = 2; p <= Hcut; p++) assert(Hc[p] >= Hc[p + 1]);

        Vn.swap(Vc); Rn.swap(Rc); Hn.swap(Hc);
        gVn.swap(gVc); gRn.swap(gRc); gHn.swap(gHc);
        stopVn.swap(stopVc); stopHn.swap(stopHc);
    }

    BoxResult res;
    res.M_lo = max(Vn[2][1], H_raw_lower_int(1, 0, thi_n, thi_d, m0));
    res.H_host.assign(Hcut + 2, res.M_lo);
    for (int p = 2; p <= Hcut + 1; p++) res.H_host[p] = max(Hn[p], res.M_lo);
    res.L_lo.assign(Rmax + 2, 0);
    for (int k = 2; k <= Rmax; k++) res.L_lo[k] = Rn[1][k];

    // --- 2. UPPER SUPERMARTINGALE at tau = tau_lo ---
    auto U_hub = solve_static_hub_upper_int(Kmax);
    if (tlo_n == 0) {
        res.U_hi = U_hub;
        res.U_tail_cap_0 = S;
        return res;
    }

    // Terminal stopping caps valid for all t >= t_M = M_steps / H_DEN:
    vector<int64_t> P_stop(Kmax + 2, 0);
    vector<vector<int64_t>> U_stop(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    for (int k = 2; k <= Kmax + 1; k++) {
        P_stop[k] = B_raw_upper_int(1, 1, k, M_steps, tlo_n, tlo_d, m0);
        int64_t cap = B_raw_upper_int(0, 1, k, M_steps, tlo_n, tlo_d, m0);
        for (int j = 0; j < k; j++) U_stop[k][j] = min(cap, U_hub[k][j]);
    }
    for (int k = Kmax; k >= 2; k--) P_stop[k] = max(P_stop[k], P_stop[k + 1]);

    // Exact 1-step terminal knot lift at t_M^- so stepping from t < t_M to t' >= t_M is a valid supermartingale!
    // For t in [t_{M-1}, t_M), 1/N < 1/(N-1) <= 1 / ((m0 - 1) * exp(t_{M-1})).
    i128 dt_term_hi = div_cl((i128)S2 * S2, (i128)(m0 - 1) * EXP_GRID[M_steps - 1].lo);
    vector<int64_t> Pn(Kmax + 2, 0);
    vector<vector<int64_t>> Un(Kmax + 2, vector<int64_t>(Kmax + 2, 0));

    // Boundary k = Kmax + 1 at t_M^- uses M_steps - 1 (valid for all t in [t_{M-1}, t_M]!)
    Pn[Kmax + 1] = max(P_stop[Kmax + 1], B_raw_upper_int(1, 1, Kmax + 1, M_steps - 1, tlo_n, tlo_d, m0));
    {
        int64_t cap_top_m1 = B_raw_upper_int(0, 1, Kmax + 1, M_steps - 1, tlo_n, tlo_d, m0);
        for (int j = 0; j <= Kmax; j++) {
            Un[Kmax + 1][j] = max(U_stop[Kmax + 1][j], min(cap_top_m1, U_hub[Kmax + 1][j]));
        }
    }
    for (int k = Kmax; k >= 2; k--) {
        int64_t bg_M = B_raw_upper_int(2, 1, k, M_steps, tlo_n, tlo_d, m0);
        int64_t bl_M = (k == 2) ? B_raw_upper_int(1, 0, 1, M_steps, tlo_n, tlo_d, m0)
                                : B_raw_upper_int(0, k + 2, 2 * k, M_steps, tlo_n, tlo_d, m0);
        i128 max_drift = 0;
        for (int gw2 : {1, 2 * k - 1}) {
            // Multiply by 2*k: gain=2*(k+1), leaf=k, host=k*gw2
            i128 num = (i128)2 * (k + 1) * (bg_M - P_stop[k])
                     + (i128)k * (bl_M - P_stop[k])
                     + (i128)k * gw2 * (P_stop[k + 1] - P_stop[k]);
            if (num > 0) {
                i128 d = div_cl(num, 2 * k);
                if (d > max_drift) max_drift = d;
            }
        }
        Pn[k] = P_stop[k] + (int64_t)div_cl(max_drift * dt_term_hi, S2);
        if (Pn[k] < Pn[k + 1]) Pn[k] = Pn[k + 1];
    }
    for (int k = Kmax; k >= 2; k--) {
        for (int j = 0; j < k; j++) {
            int64_t lp = (j > 0) ? U_stop[k][j - 1] : 0;
            i128 max_drift = 0;
            for (int g2 : {2 * (j + 1), k + j + 1}) {
                i128 num = (i128)2 * (P_stop[k] - U_stop[k][j])
                         + (i128)2 * j * (lp - U_stop[k][j])
                         + (i128)k * g2 * (U_stop[k + 1][j + 1] - U_stop[k][j]);
                if (num > 0) {
                    i128 d = div_cl(num, 2 * k);
                    if (d > max_drift) max_drift = d;
                }
            }
            Un[k][j] = U_stop[k][j] + (int64_t)div_cl(max_drift * dt_term_hi, S2);
        }
    }

    // Backward induction over slabs i = M_steps - 1 down to 0
    for (int i = M_steps - 1; i >= 0; i--) {
        vector<vector<int64_t>> Uc(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
        vector<int64_t> Pc(Kmax + 2, 0);

        // Stopped boundary k = Kmax + 1: evaluate at max(0, i-1) so value at t_i >= B(t_{i-1}) valid on [t_{i-1}, t_i]!
        int t_bnd = max(0, i - 1);
        Pc[Kmax + 1] = max(Pn[Kmax + 1], B_raw_upper_int(1, 1, Kmax + 1, t_bnd, tlo_n, tlo_d, m0));
        int64_t cap_top = B_raw_upper_int(0, 1, Kmax + 1, t_bnd, tlo_n, tlo_d, m0);
        for (int j = 0; j <= Kmax; j++) {
            Uc[Kmax + 1][j] = max(Un[Kmax + 1][j], min(cap_top, U_hub[Kmax + 1][j]));
        }

        // Post-promotion P[k] for k = Kmax down to 2 (targets bg_i[k], bl_i[k], P[k+1] ALREADY computed!)
        for (int k = Kmax; k >= 2; k--) {
            int64_t bg_k = B_raw_upper_int(2, 1, k, i, tlo_n, tlo_d, m0);
            int64_t bl_k = (k == 2) ? B_raw_upper_int(1, 0, 1, i, tlo_n, tlo_d, m0)
                                    : B_raw_upper_int(0, k + 2, 2 * k, i, tlo_n, tlo_d, m0);
            i128 beta = 0;
            for (int gw2 : {1, 2 * k - 1}) {
                // Multiply rates by 2*k: rg = 2*(k+1), rl = k, rhost = k*gw2, rtot = 3*k + 2 + k*gw2
                i128 rg = 2 * (k + 1);
                i128 rl = k;
                i128 rh = (i128)k * gw2;
                i128 rtot = rg + rl + rh;
                // (U-R) at t_{i+1}:
                i128 num_R = rg * bg_k + rl * bl_k + rh * Pn[k + 1] - rtot * Pn[k];
                if (num_R > 0) {
                    i128 bR = div_cl(num_R, 2 * k);
                    if (bR > beta) beta = bR;
                }
                // (U-L) at t_i:
                i128 num_L = (rg * bg_k + rl * bl_k + rh * Pc[k + 1] - rtot * Pn[k]) * H_DEN;
                if (num_L > 0) {
                    i128 den_L = (i128)2 * k * H_DEN + rtot;
                    i128 bL = div_cl(num_L, den_L);
                    if (bL > beta) beta = bL;
                }
            }
            int64_t delta = (int64_t)div_cl(beta, H_DEN);
            int64_t pc = Pn[k] + delta;
            if (pc < Pc[k + 1]) pc = Pc[k + 1];
            Pc[k] = pc;
        }

        // Pending U[k][j] for k = Kmax down to 2, j = 0 up to k-1 (targets P[k], U[k][j-1], U[k+1][j+1] ALREADY computed!)
        for (int k = Kmax; k >= 2; k--) {
            for (int j = 0; j < k; j++) {
                int64_t prom_n = Pn[k], prom_c = Pc[k];
                int64_t loss_n = (j > 0) ? Un[k][j - 1] : 0;
                int64_t loss_c = (j > 0) ? Uc[k][j - 1] : 0;
                int64_t gain_n = Un[k + 1][j + 1];
                int64_t gain_c = Uc[k + 1][j + 1];
                i128 beta = 0;
                for (int g2 : {2 * (j + 1), k + j + 1}) {
                    i128 rtot = 2 * (1 + j) + (i128)k * g2;
                    // (U-R) at t_{i+1}:
                    i128 num_R = (i128)2 * prom_n + (i128)2 * j * loss_n + (i128)k * g2 * gain_n - rtot * Un[k][j];
                    if (num_R > 0) {
                        i128 bR = div_cl(num_R, 2 * k);
                        if (bR > beta) beta = bR;
                    }
                    // (U-L) at t_i:
                    i128 num_L = ((i128)2 * prom_c + (i128)2 * j * loss_c + (i128)k * g2 * gain_c - rtot * Un[k][j]) * H_DEN;
                    if (num_L > 0) {
                        i128 den_L = (i128)2 * k * H_DEN + rtot;
                        i128 bL = div_cl(num_L, den_L);
                        if (bL > beta) beta = bL;
                    }
                }
                int64_t delta = (int64_t)div_cl(beta, H_DEN);
                Uc[k][j] = Un[k][j] + delta;
            }
        }

        for (int k = 2; k <= Kmax; k++) assert(Pc[k] >= Pc[k + 1]);
        Un.swap(Uc); Pn.swap(Pc);
    }

    res.U_hi.assign(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    for (int k = 2; k <= Kmax + 1; k++) {
        int64_t cap0 = B_raw_upper_int(0, 1, k, 0, tlo_n, tlo_d, m0);
        for (int j = 0; j < k; j++) res.U_hi[k][j] = min({Un[k][j], cap0, U_hub[k][j]});
    }
    res.U_tail_cap_0 = B_raw_upper_int(0, 1, Kmax + 2, 0, tlo_n, tlo_d, m0);
    return res;
}

// --- Exact Integer All-Tree Signed Ledger Verifier ---
// All dimensionless values normalized by M_lo in integer scale L_SC = 10^9!
// Outward rounding:
//   H_tilde[p] = floor(H_host[p] * L_SC / M_lo) (lower!)
//   L_tilde[k] = floor(L_lo[k] * L_SC / M_lo)   (lower!)
//   U_tilde[k][j] = ceil(U_hi[k][j] * L_SC / M_lo) (upper!)
static constexpr int64_t LSC = 1000000000LL; // 10^9

struct Tariffs {
    int64_t a, b_base, fee_lo, fee_hi, fee_t, xfer; // in units of 1/LSC
    int D_hi;
};

static bool verify_all_tree_ledger(const BoxResult& B,
                                   int Kmax, int Rmax, int Hcut, string& info) {
    int64_t M_lo = B.M_lo;
    if (M_lo <= 0) { info = "M_lo<=0"; return false; }

    auto H = [&](int p) -> int64_t {
        if (p <= Hcut + 1) return (int64_t)div_fl((i128)B.H_host[p] * LSC, M_lo);
        return LSC;
    };
    auto L = [&](int k) -> int64_t {
        if (k <= Rmax) return (int64_t)div_fl((i128)B.L_lo[k] * LSC, M_lo);
        return 0;
    };
    int64_t u_tail_beta = (int64_t)div_cl((i128)B.U_tail_cap_0 * LSC, M_lo);
    int64_t u_ktop_0 = (int64_t)div_cl((i128)B.U_hi[Kmax + 1][0] * LSC, M_lo);
    int64_t u_tail0_cap = min(u_tail_beta, u_ktop_0);

    // Precompute U_tab[k][j] for k <= Kmax+1
    vector<vector<int64_t>> Utab(Kmax + 2, vector<int64_t>(Kmax + 2, 0));
    for (int k = 2; k <= Kmax + 1; k++) {
        for (int j = 0; j < k; j++) {
            Utab[k][j] = (int64_t)div_cl((i128)B.U_hi[k][j] * LSC, M_lo);
        }
    }
    auto U = [&](int k, int j) -> int64_t {
        if (k <= Kmax + 1) return Utab[k][min(j, k - 1)];
        // For k > Kmax+1: min(1/((j+1)*M_lo), u_tail0_cap)
        // 1/((j+1)*M_lo) in units of LSC is ceil(S * LSC / ((j+1)*M_lo))
        int64_t c_hub = (int64_t)div_cl((i128)S * LSC, (i128)(j + 1) * M_lo);
        return min(c_hub, u_tail0_cap);
    };

    int64_t u20 = U(2, 0);
    int64_t g_base = LSC - u20;
    if (g_base <= 0) { info = "g_base<=0"; return false; }

    // Fast check: one-line universal inequality 3*LSC > 4*u20 + 2*U3_max (with u20 >= U3_max >= max_{k>=3} U(k,j))!
    int64_t U3_max = u_tail0_cap;
    for (int k = 3; k <= Kmax + 1; k++)
        for (int j = 0; j < k; j++)
            U3_max = max(U3_max, Utab[k][j]);
    if (u20 >= U3_max && (i128)3 * LSC > (i128)4 * u20 + (i128)2 * U3_max) {
        info = "ONE-LINE EXACT PASS!";
        return true;
    }

    vector<Tariffs> tlist = {
        { 9200000, 22500000, 4800000, 8200000, 27500000, 4600000, 10},
        { 8800000, 22500000, 4800000, 8100000, 27500000, 4800000, 10},
        {10000000, 22000000, 5200000, 8200000, 27500000, 4800000, 10},
        {11000000, 22000000, 5500000, 8400000, 27500000, 5000000, 10},
        {12500000, 21500000, 6000000, 8600000, 27500000, 5500000, 10},
        {15000000, 20000000, 7500000, 9000000, 27500000, 6000000, 10},
    };

    // Exact rational site credit lower bound C_lo(d, c, r_num, r_den) in units of 1/LSC!
    // s = (l*r_den + r_num) / (r_den * d * (d+1)), t = 1 / (2*d)
    // C = d*s*(H(d+1) - u20) + l*t*(H(2) - u20) - l*(t-s)_+*(u20 - L(d+1))_+ - l*(s-t)_+*(U(d+1,l-1) - L(2))_+
    auto site_credit_lo = [&](int d, int c, i128 r_num, i128 r_den, i128& s_num_out, i128& s_den_out) -> i128 {
        int l = d - c;
        i128 sn = (i128)l * r_den + r_num;
        i128 sd = r_den * d * (d + 1);
        s_num_out = sn; s_den_out = sd;
        int64_t hd = H(d + 1) - u20;
        int64_t h2 = H(2) - u20;
        i128 C = div_fl((i128)d * sn * hd, sd) + div_fl((i128)l * h2, 2 * d);
        if (l > 0) {
            // Compare t = 1/(2d) vs s = sn/sd: t - s = (sd - 2*d*sn) / (2*d*sd)
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

    // Precompute exact D_tail = 300
    const int D_tail = 300;
    string last_fail;

    for (size_t ti = 0; ti < tlist.size(); ti++) {
        const auto& T = tlist[ti];
        if (H(2) <= u20) { last_fail = "seed"; continue; }
        bool ok = true;

        // 1. Stars d = 2..D_tail
        for (int d = 2; d <= D_tail; d++) {
            i128 sn, sd;
            if (site_credit_lo(d, 0, 0, 1, sn, sd) <= 0) {
                ok = false; last_fail = "star_" + to_string(d); break;
            }
        }
        if (!ok) continue;
        // Analytic star tail d > D_tail: (14/10)*g_base - (1/2)*U(D_tail+2, D_tail) > 0
        if ((i128)14 * g_base <= (i128)5 * U(D_tail + 2, D_tail)) {
            last_fail = "star_tail"; continue;
        }

        // 2. Finite sites d = 2..D_tail, c = 1..d (excluding bare deg-2 d=c=2)
        for (int d = 2; d <= D_tail && ok; d++) {
            for (int c = 1; c <= d && ok; c++) {
                if (d == 2 && c == 2) continue;
                int l = d - c;
                int64_t u_adv = U(d + 1, l);
                if (d >= T.D_hi) {
                    // Check if s_max * u_adv <= T.b_base: if so, exc == 0 identically, so R is jointly concave
                    // on the (k, p_pay, r) polytope and 3 simplex vertices suffice! Otherwise loop ALL k = 0..c
                    // (for fixed k, R - exc IS jointly concave in (p_pay, r), so 2x2 endpoints are exact!).
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
                                if (R < 0) {
                                    ok = false;
                                    last_fail = "hi_site_d" + to_string(d) + "_c" + to_string(c);
                                }
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
                                    i128 s_u_hi = div_cl(sn * u_adv, sd);
                                    i128 exc = (i128)(c - k) * max((i128)0, s_u_hi - T.b_base);
                                    if (R < exc) {
                                        ok = false;
                                        last_fail = "hi_site_d" + to_string(d) + "_c" + to_string(c);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    int64_t F_ch = (d == 2) ? T.fee_t : T.fee_lo;
                    // Coarse test first
                    bool coarse_ok = true;
                    for (int k = 0; k <= c && coarse_ok; k++) {
                        i128 r_nums[2] = {k, c};
                        for (int ri = 0; ri < 2 && coarse_ok; ri++) {
                            i128 sn, sd;
                            i128 C = site_credit_lo(d, c, r_nums[ri], 2, sn, sd);
                            i128 R = C + (i128)T.a * (c - 2) - T.b_base - (i128)k * F_ch;
                            i128 s_u_hi = div_cl(sn * u_adv, sd);
                            i128 exc = (i128)(c - k) * max((i128)0, s_u_hi - T.b_base);
                            if (R < exc) coarse_ok = false;
                        }
                    }
                    if (coarse_ok) continue;

                    // Detailed multiset enumeration for direct neighbors in {0 (>=D_hi), 2..D_hi-1}
                    // Common denominator for exact receive r: lcm(2, 3, ..., D_hi)
                    const i128 LCM = 2520; // divisible by 2,3,4,5,6,7,8,9,10!
                    assert(T.D_hi <= 10);
                    vector<int> es;
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
                                        i128 su = div_cl(sn * u_adv, sd);
                                        exc += max((i128)0, su - T.b_base);
                                    } else {
                                        // s - 1/(d*e*(e+1)) = (sn*e*(e+1) - LCM*(d+1)) / (sd*e*(e+1))
                                        // Wait: sd = LCM * d * (d+1), so 1/(d*e*(e+1)) = LCM*(d+1) / (sd*e*(e+1))!
                                        i128 diff_n = sn * e * (e + 1) - LCM * (d + 1);
                                        if (diff_n > 0) {
                                            int64_t slope = max((int64_t)0, u_adv - L(e + 1));
                                            i128 j_uv = div_cl(diff_n * slope, sd * e * (e + 1));
                                            exc += max((i128)0, j_uv - T.b_base);
                                        }
                                    }
                                }
                                if (R < exc) {
                                    ok = false;
                                    last_fail = "lo_site_d" + to_string(d) + "_c" + to_string(c);
                                }
                            }
                            return;
                        }
                        // Choices: 0 (meaning >= D_hi), 2, 3, ..., D_hi - 1
                        for (int idx = start_idx; idx < T.D_hi; idx++) {
                            if (idx == 1) continue;
                            es.push_back(idx);
                            dfs(k, rem - 1, idx);
                            es.pop_back();
                            if (!ok) return;
                        }
                    };
                    for (int k = 0; k <= c && ok; k++) {
                        dfs(k, c - k, 0);
                    }
                }
            }
        }
        if (!ok) continue;

        // 3. Analytic high-degree site tail d > D_tail (exact integer inequalities!)
        //   (i) Outgoing exposure s*U(d+1,l) <= U(D_tail+2, 0)/(D_tail+2) <= b_base
        int64_t u_tail_dt = U(D_tail + 2, 0);
        if (div_cl(u_tail_dt, D_tail + 2) > T.b_base) { last_fail = "hi_tail_exc"; continue; }
        if (T.a < max(T.fee_hi, T.xfer)) { last_fail = "hi_tail_a"; continue; }
        //   (ii) Core-side leaf debit T_core_leaf = min(ceil(S*LSC/(2*(D_tail+2)*M_lo)), ceil(u_tail_dt/2))
        i128 T_cl_1 = div_cl((i128)S * LSC, (i128)2 * (D_tail + 2) * M_lo);
        i128 T_cl_2 = div_cl(u_tail_dt, 2);
        i128 T_cl = min(T_cl_1, T_cl_2);
        //   (iii) Case A (l <= d/2, c >= d/2 >= D_tail/2):
        //     z_lo = floor(D_tail*g_base/(D_tail+1)) + floor((g_base - u20)/2)
        i128 z_lo = div_fl((i128)D_tail * g_base, D_tail + 1) + div_fl((i128)(g_base - u20), 2);
        i128 stock_A = min((i128)0, div_fl(z_lo, 2));
        i128 net_A = stock_A - T_cl + (i128)(T.a - max(T.fee_hi, T.xfer)) * (D_tail / 2) - 2 * T.a - T.b_base;
        //   (iv) Case B (l > d/2):
        //     serial >= floor(D_tail*g_base/(2*(D_tail+1))) + floor(g_base/4)
        //     leaf-side debit <= ceil(u20 / (2*(D_tail+2)))
        i128 ser_B = div_fl((i128)D_tail * g_base, 2 * (D_tail + 1)) + div_fl(g_base, 4);
        i128 net_B = ser_B - div_cl(u20, 2 * (D_tail + 2)) - T_cl - 2 * T.a - T.b_base;
        if (net_A <= 0 || net_B <= 0) { last_fail = "hi_tail_AB"; continue; }

        // 4. Bare degree-2 chains of lengths 1, 2, 3 (exact integer arithmetic!)
        // Scale s and q by X = 6 * 10^8 (divisible by 12, 6, ...), values by LSC = 10^9!
        // Note: all endpoints use EXACT intervals [q_lo, q_hi] and [s_lo, s_hi] in units of 1/X,
        // with high tail split into its TWO exact convex corners q = 0 and q = 1/(6*(Kmax+1))!
        const i128 X = 600000000LL; // 6 * 10^8
        struct End {
            i128 q_lo, q_hi;     // 1/(6d) in units of 1/X
            i128 s_lo, s_hi;     // s_e in units of 1/X
            int64_t w_ch;        // (U(3,0) - L(d+1))_+ in units of 1/LSC
            i128 own_or_we;      // if is_cap: own_cap in units of 1/(X*LSC); else (U(d+1,l)-L(3))_+ in 1/LSC
            i128 fee_XL;         // fee * X in units of 1/(X*LSC)
            bool is_cap;
        };
        vector<End> ends;
        int64_t u30 = U(3, 0);
        int64_t L3 = L(3);
        for (int d = 2; d < T.D_hi; d++) {
            int64_t F_ch = (d == 2) ? T.fee_t : T.fee_lo;
            i128 ql = div_fl(X, 6 * d), qh = div_cl(X, 6 * d);
            int64_t wch = max((int64_t)0, u30 - L(d + 1));
            for (int c = 1; c <= d; c++) {
                if (d == 2 && c == 2) continue;
                int l = d - c;
                int64_t we = max((int64_t)0, U(d + 1, l) - L3);
                for (int r2 : {1, c}) {
                    // s_e = (2*l + r2) / (2*d*(d+1))
                    i128 sl = div_fl((i128)(2 * l + r2) * X, 2 * d * (d + 1));
                    i128 sh = div_cl((i128)(2 * l + r2) * X, 2 * d * (d + 1));
                    ends.push_back({ql, qh, sl, sh, wch, we, (i128)F_ch * X, false});
                }
            }
        }
        for (int d = T.D_hi; d <= Kmax; d++) {
            i128 ql = div_fl(X, 6 * d), qh = div_cl(X, 6 * d);
            i128 sl = div_fl(X, 2 * d * (d + 1));
            int64_t wch = max((int64_t)0, u30 - L(d + 1));
            i128 own_cap = 0;
            for (int l = 0; l < d; l++) {
                int64_t we = max((int64_t)0, U(d + 1, l) - L3);
                // ((d+l)/(2*d*(d+1))) * X * we
                i128 cand = div_cl((i128)(d + l) * X * we, 2 * d * (d + 1));
                if (cand > own_cap) own_cap = cand;
            }
            ends.push_back({ql, qh, sl, 0, wch, own_cap, (i128)T.fee_hi * X, true});
        }
        // High tail d > Kmax: exact stock-coupled own_cap_tail in units of 1/(X*LSC)
        i128 oc1 = div_cl((i128)u_tail0_cap * X, 2 * (Kmax + 2));
        i128 oc2a = div_cl((i128)(Kmax + 2) * S * LSC * X, (i128)4 * (Kmax + 1) * (Kmax + 2) * M_lo);
        i128 oc2b = div_cl((i128)u_tail0_cap * X, Kmax + 2);
        i128 oc_tail = max(oc1, min(oc2a, oc2b));
        // Two exact convex corners q = 0 and q = 1/(6*(Kmax+1))!
        ends.push_back({0, 0, 0, 0, u30, oc_tail, (i128)T.fee_hi * X, true});
        i128 qtl = div_fl(X, 6 * (Kmax + 1)), qth = div_cl(X, 6 * (Kmax + 1));
        ends.push_back({qtl, qth, 0, 0, u30, oc_tail, (i128)T.fee_hi * X, true});

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
                        i128 mid_diff = max(s2h - s1l, s1h - s2l);
                        debt = port_debt(s1l, s1h, ex) + port_debt(s2l, s2h, ey)
                             + mid_diff * v_int - (s1l + s2l) * cred_site;
                    } else {
                        i128 s1l = X12 + ex.q_lo, s1h = X12 + ex.q_hi;
                        i128 s2l = X6, s2h = X6;
                        i128 s3l = X12 + ey.q_lo, s3h = X12 + ey.q_hi;
                        i128 d12 = max(s2h - s1l, s1h - s2l);
                        i128 d23 = max(s3h - s2l, s2h - s3l);
                        debt = port_debt(s1l, s1h, ex) + port_debt(s3l, s3h, ey)
                             + (d12 + d23) * v_int - (s1l + s2l + s3l) * cred_site;
                    }
                    if (debt >= base_XL) {
                        ok = false;
                        last_fail = "chain_len" + to_string(len);
                    }
                }
            }
        }
        if (ok) {
            info = "EXACT PASS (tariff#" + to_string(ti) + ")";
            return true;
        }
    }
    info = "FAIL:" + last_fail;
    return false;
}

// Self-contained exact rational verification of the high-tau tail tau >= 6 (for ALL b >= 3)!
// Proves: for all b >= 3 and 0 < x < 1 with b*log(1/(1-x)) >= 6 (hence b*x/(1-x) >= 6):
//   e^{-1} - (4/3)*B_{m,1/2}(x)/q^b - (2/3)*B_{m,1/3}(x)/q^b > 19/10000 > 0!
static void verify_high_tau_tail_exact() {
    // 1. e^{-1} > 367879/1000000
    Iv inv_e = iv_exp_neg({S2, S2});
    i128 M_lo = div_fl(inv_e.lo * 1000000, S2); // 367879
    assert(M_lo >= 367879);
    // 2. For r in {1/2, 1/3}, s = 1 - r, B = m - r >= 7/2 (since m = b+1 >= 4), v = q/(Bx) <= 1/6:
    //    B_{m,r}(x)/q^b <= (v^s / Gamma(r)) * (1 - d_s * v), where d_s = (7/9)*s*(1 - (1+s)/6)
    //    Since v^s * (1 - d_s*v) increases on [0, 1/6] (deriv = v^{s-1}*(s - (s+1)*d_s*v) >= 0 since (s+1)*d_s/6 < s!),
    //    evaluate at v = 1/6!
    // Check r = 1/2 (s = 1/2): d_{1/2} = (7/9)*(1/2)*(3/4) = 7/24. (1 - d_{1/2}/6) = 137/144.
    //   Gamma(1/2) = sqrt(pi): check via iv_gamma_rat(1, 2)!
    Iv gam12 = iv_gamma_rat(1, 2);
    Iv gam13 = iv_gamma_rat(1, 3);
    Iv v16_12 = iv_pow_rat(iv_from_rat(1, 6), 1, 2); // 6^{-1/2}
    Iv v16_23 = iv_pow_rat(iv_from_rat(1, 6), 2, 3); // 6^{-2/3}
    // Upper bound on (4/3) * B_{m,1/2}/q^b:
    //   (4/3) * (v16_12.hi / gam12.lo) * (137/144)
    i128 u12_hi = div_cl(div_cl(v16_12.hi * S2, gam12.lo) * 4 * 137, 3 * 144);
    // For r = 1/3 (s = 2/3): d_{2/3} = (7/9)*(2/3)*(1 - 5/18) = (14/27)*(13/18) = 91/243.
    //   1 - d_{2/3}/6 = 1 - 91/1458 = 1367/1458.
    i128 u13_hi = div_cl(div_cl(v16_23.hi * S2, gam13.lo) * 2 * 1367, 3 * 1458);
    i128 margin = inv_e.lo - u12_hi - u13_hi;
    assert(margin > (i128)19 * S2 / 10000);
    cout << "[HIGH-TAU EXACT] PASS tau >= 6 (all b >= 3, all trees): inv_e="
         << fixed << setprecision(7) << (double)inv_e.lo / S2
         << " - debit_hi=" << (double)(u12_hi + u13_hi) / S2
         << " = margin +" << (double)margin / S2 << " > 0.0019\n";
}

int main() {
    init_exp_grid(1600);

    // Sanity cross-checks of exact transcendental primitives
    {
        Iv g12 = iv_gamma_rat(1, 2);
        Iv g13 = iv_gamma_rat(1, 3);
        double sqpi = sqrt(M_PI);
        assert((double)g12.lo / S2 <= sqpi && sqpi <= (double)g12.hi / S2);
        double g13_ref = tgamma(1.0 / 3.0);
        assert((double)g13.lo / S2 <= g13_ref && g13_ref <= (double)g13.hi / S2);
        cout << "[PRIMITIVES] Gamma(1/2) in [" << fixed << setprecision(12)
             << (double)g12.lo / S2 << ", " << (double)g12.hi / S2 << "]\n";
        cout << "[PRIMITIVES] Gamma(1/3) in ["
             << (double)g13.lo / S2 << ", " << (double)g13.hi / S2 << "]\n";
    }

    verify_high_tau_tail_exact();

    const int m0 = 200;
    // Adaptive exact rational box cover of [0, 6] (using rational endpoints p/Q with Q = 100000!)
    const int64_t Q = 100000;
    vector<pair<int64_t, int64_t>> stack_boxes = {
        {0, 25}, {25, 100}, {100, 300}, {300, 800}, {800, 2000},
        {2000, 5000}, {5000, 12000}, {12000, 25000}, {25000, 45000},
        {45000, 75000}, {75000, 110000}, {110000, 150000}, {150000, 200000},
        {200000, 260000}, {260000, 330000}, {330000, 410000}, {410000, 500000},
        {500000, 600000}
    };
    reverse(stack_boxes.begin(), stack_boxes.end());

    int n_pass = 0, n_fail = 0;
    int64_t covered_up_to = 0;
    while (!stack_boxes.empty()) {
        auto [a, b] = stack_boxes.back();
        stack_boxes.pop_back();
        double thi = (double)b / Q;
        // Choose horizon T and cutoffs based on box regime
        int M_steps, Dmax, Rmax, Kmax, Hcut = 40;
        if (a == 0) {
            M_steps = 11 * H_DEN; Dmax = 24; Rmax = 32; Kmax = 32;
        } else if (thi <= 0.01) {
            M_steps = 9 * H_DEN; Dmax = 24; Rmax = 30; Kmax = 30;
        } else if (thi <= 0.15) {
            M_steps = 6 * H_DEN; Dmax = 22; Rmax = 26; Kmax = 28;
        } else {
            M_steps = 4 * H_DEN; Dmax = 22; Rmax = 24; Kmax = 26;
        }
        BoxResult br = solve_box_int(a, Q, b, Q, m0, M_steps, Dmax, Rmax, Kmax, Hcut);
        string info;
        bool ok = verify_all_tree_ledger(br, Kmax, Rmax, Hcut, info);
        if (ok) {
            assert(a == covered_up_to);
            covered_up_to = b;
            n_pass++;
            cout << "  PASS [" << fixed << setprecision(5) << (double)a / Q << ", " << (double)b / Q
                 << "] | M=" << setprecision(5) << (double)br.M_lo / S
                 << " U20/M=" << setprecision(4) << (double)br.U_hi[2][0] / br.M_lo
                 << " | " << info << "\n" << flush;
        } else {
            if (b - a <= 2) {
                n_fail++;
                cout << "  HARD FAIL [" << fixed << setprecision(5) << (double)a / Q << ", " << (double)b / Q
                     << "] | M=" << setprecision(5) << (double)br.M_lo / S
                     << " U20/M=" << setprecision(4) << (double)br.U_hi[2][0] / br.M_lo
                     << " | " << info << "\n" << flush;
            } else {
                int64_t mid = (a + b) / 2;
                stack_boxes.push_back({mid, b});
                stack_boxes.push_back({a, mid});
            }
        }
    }
    assert(n_fail > 0 || covered_up_to == 6 * Q);
    cout << "\n=== CERTIFICATE SUMMARY (ALL m >= " << m0 << ", ALL x in [0, 1)) ===\n";
    cout << "Boxes PASSED: " << n_pass << ", FAILED: " << n_fail
         << ", exact contiguous tau-coverage: [0, " << covered_up_to << "/" << Q << "]\n";
    return (n_fail == 0) ? 0 : 1;
}
