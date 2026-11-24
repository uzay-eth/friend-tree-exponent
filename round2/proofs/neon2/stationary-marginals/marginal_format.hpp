#pragma once
#include "../quintet-cert/quintet_format.hpp"

/*
 * FT2_M5V1: Exact zero-floating-point unoriented 5-vertex stationary-marginal
 * coboundary certificate format.
 *
 * Wraps the exact 512-byte QuintetHeader (with base.mode == 1, alpha_lo_num == 0,
 * and f_p0/f_p1/f_p2 signed in [-Q5_MAX_W, Q5_MAX_W]) plus a 128-byte extension
 * header specifying the bounded vertex observable psi(D, c) in {0, 1}, sign
 * (+1 upper bound, -1 lower bound), exact rational nominal target theta = theta_P / theta_Q,
 * 4th analytic tail coefficient f_p3, and certified mu bracket [mu_lo_P/mu_lo_Q, mu_hi_P/mu_hi_Q].
 *
 * Signed analytic singleton tail for D >= DF (q = L/D, x = 1/D):
 *   w(q, x) = S * f(D, c) = f_p0 + f_p1 * q + f_p2 * q^2 + f_p3 * q * x.
 */

static constexpr uint64_t M5_MAGIC = 0x4654325F4D355631ULL; // "FT2_M5V1"
static constexpr int64_t  M5_MAX_OBS = 8;

struct MarginalHeader {
    uint64_t m_magic;       // M5_MAGIC
    int64_t  obs_id;        // 1..8
    int64_t  obs_sign;      // +1 (upper bound on pi(psi)) or -1 (lower bound on pi(psi))
    int64_t  theta_P;       // nominal target theta = theta_P / theta_Q
    int64_t  theta_Q;       // 1 <= theta_Q <= Q5_MAX_Q
    int64_t  f_p3;          // coefficient of q*x in w(q,x) = S*f(D,c)
    int64_t  mu_lo_P;       // certified lower endpoint of mu
    int64_t  mu_lo_Q;
    int64_t  mu_hi_P;       // certified upper endpoint of mu
    int64_t  mu_hi_Q;
    int64_t  m_res[6];      // must be 0
    QuintetHeader base;     // exact 512-byte unoriented 5-vertex table layout (base.mode == 1)
};
static_assert(sizeof(MarginalHeader) == 640, "MarginalHeader must be exactly 640 bytes");

inline int m5_eval_obs(int64_t obs_id, int D, int c) {
    // Only evaluated for core vertices D >= 2, 0 <= c <= D.
    // All observables vanish identically on D >= DF (for obs_id == 8, bare cores with D >= DF
    // are bounded separately by the exact analytic tail (40 / mu_lo) * 2^(-DF)).
    switch (obs_id) {
        case 1: return (D == 2) ? 1 : 0;                  // p_2 = pi(D_o = 2)
        case 2: return (D == 3) ? 1 : 0;                  // p_3 = pi(D_o = 3)
        case 3: return (D == 4) ? 1 : 0;                  // p_4 = pi(D_o = 4)
        case 4: return (D == 5) ? 1 : 0;                  // p_5 = pi(D_o = 5)
        case 5: return (D == 6) ? 1 : 0;                  // p_6 = pi(D_o = 6)
        case 6: return (D == 2 && c == 1) ? 1 : 0;        // pi(J) = pi(D_o = 2, c_o = 1)
        case 7: return (D == 2 && c == 2) ? 1 : 0;        // pi(2, 2) = bare degree-2 fraction
        case 8: return (D == c) ? 1 : 0;                  // bare cores L_o = 0 (restricted to D < DF)
        default: return 0;
    }
}

inline const char* m5_obs_name(int64_t obs_id) {
    switch (obs_id) {
        case 1: return "p2_D2";
        case 2: return "p3_D3";
        case 3: return "p4_D4";
        case 4: return "p5_D5";
        case 5: return "p6_D6";
        case 6: return "pi_J_2_1";
        case 7: return "pi_bare_2_2";
        case 8: return "pi_bare_L0_fin";
        default: return "unknown";
    }
}

inline int m5_validate_header(const MarginalHeader& mh) {
    if (mh.m_magic != M5_MAGIC) return 101;
    if (mh.obs_id < 1 || mh.obs_id > M5_MAX_OBS) return 102;
    if (mh.obs_sign != 1 && mh.obs_sign != -1) return 103;
    if (mh.theta_Q <= 0 || mh.theta_Q > Q5_MAX_Q || mh.theta_P < 0 || mh.theta_P > mh.theta_Q) return 104;
    if (mh.f_p3 < -Q5_MAX_W || mh.f_p3 > Q5_MAX_W) return 105;
    if (mh.mu_lo_P <= 0 || mh.mu_lo_Q <= 0 || mh.mu_lo_P >= mh.mu_lo_Q || mh.mu_lo_Q > Q5_MAX_Q) return 106;
    if (mh.mu_hi_P <= 0 || mh.mu_hi_Q <= 0 || mh.mu_hi_P >= mh.mu_hi_Q || mh.mu_hi_Q > Q5_MAX_Q) return 107;
    if ((q5_i128)mh.mu_lo_P * mh.mu_hi_Q > (q5_i128)mh.mu_hi_P * mh.mu_lo_Q) return 108;
    for (int i = 0; i < 6; ++i) if (mh.m_res[i] != 0) return 109;
    // Validate base QuintetHeader (note: base.f_p2 is allowed nonzero for signed quadratic tail,
    // so validate base with f_p2 temporarily zeroed, then check base.f_p2 range directly).
    if (mh.base.mode != 1 || mh.base.alpha_lo_num != 0) return 110;
    if (mh.base.f_p2 < -Q5_MAX_W || mh.base.f_p2 > Q5_MAX_W) return 111;
    QuintetHeader tmp = mh.base;
    tmp.f_p2 = 0;
    return q5_validate_header(tmp);
}

inline int m5_compute_expected_bytes(const MarginalHeader& mh, uint64_t& out_bytes) {
    int vh = m5_validate_header(mh);
    if (vh != 0) return vh;
    QuintetHeader tmp = mh.base;
    tmp.f_p2 = 0;
    uint64_t base_bytes = 0;
    int vb = q5_compute_expected_bytes(tmp, base_bytes);
    if (vb != 0) return vb;
    out_bytes = base_bytes + (sizeof(MarginalHeader) - sizeof(QuintetHeader));
    return 0;
}
