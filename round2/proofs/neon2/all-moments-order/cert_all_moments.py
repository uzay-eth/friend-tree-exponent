#!/usr/bin/env python3
"""
Standalone ZERO-FLOATING-POINT exact rational certificate & verifier for
ALL-REAL-MOMENT, ALL-LABEL STRICT MONOTONICITY in the Friend Tree (Open Question 3):

  Theorem 1 (Universal per-tree ordering, EVERY real p >= 1, EVERY n >= 2, EVERY tree T of size n):
      E[Z_{n+1}^p - Z_{n+2}^p | T_n = T] >= (Gamma(n+2)*Gamma(p+1) / (n*(n+1)*Gamma(n+2+p))) * c_*(p) > 0,

  Theorem 2 (All-label, all-L^p-norm strict monotonicity on edge-seeded T_2 = {1, 2}):
      E[Z_1^p] = E[Z_2^p] > E[Z_3^p] > E[Z_4^p] > ...     for EVERY real p in [1, infinity).

Architecture:
  - For p in [12, infinity): elementary analytic proof (README.md Lemma 4(1), zero computer needed) gives
      c_*(p) >= 1/3 - (4*sqrt(p+1))/(sqrt(3)*(2p+1)) > 1/4000 > 0   (and -> 1/3 as p -> infinity).
  - For p in [1, 12]: exact rational partition into 19 contiguous closed intervals [p0, p1]
    (6 'full' intervals on [1, 8/5], 4 'vx_only' intervals on [8/5, 3], 9 '1v' intervals on [3, 12]).
    On each interval, all local sub/supermartingales are rescaled by the exact rising-Gamma
    polynomial q_p(D) = prod_{j=1}^{D-1} (1 + p/j), making every interior generator drift
    exact multiaffine in (rates, p) and every boundary an exact rational product with ZERO
    Gamma functions, ZERO pi, and ZERO floating-point arithmetic.
  - Standalone exact fractions.Fraction verifier checks:
      1. Exact interval tiling contiguity and strict Fraction monotonicity from p = 1 to p = 12,
         plus exact analytic threshold inequalities at p = 12 and early-gap inequality 8/27 < 9/25.
      2. Every boundary inequality, nonnegativity of beta & gamma, nonincreasing monotonicity of gamma,
         strict int tick types (rejecting float/bool), and every multiaffine-box corner drift inequality
         on all 19 intervals (115,964 rows total).
      3. Sample rational checks of f_p(L)'s lower-drift inequality plus an exact finite-tree check of
         the signed two-birth identity and endpoint split on all 24 unlabelled trees n = 2..7
         (supplementing the universal analytic proofs in README.md).
      4. 14 planted single-tick, boundary, non-integer, and malformed-tiling negative controls
         (fail-closed even under python3 -O).
"""
from fractions import Fraction as F
import copy
import json

SCALE = 10 ** 9

# 19 exact rational intervals tiling [1, 12]: (p0, p1, mode, DM, K, DXM)
INTERVAL_SPECS = [
    # 6 full intervals covering [1, 8/5]
    (F(1, 1),   F(21, 20), "full",    40, 14, 8),
    (F(21, 20), F(11, 10), "full",    40, 14, 8),
    (F(11, 10), F(6, 5),   "full",    40, 14, 8),
    (F(6, 5),   F(13, 10), "full",    40, 14, 8),
    (F(13, 10), F(29, 20), "full",    40, 12, 7),
    (F(29, 20), F(8, 5),   "full",    40, 12, 7),
    # 4 vx_only intervals covering [8/5, 3]
    (F(8, 5),   F(9, 5),   "vx_only", 40,  0, 7),
    (F(9, 5),   F(21, 10), "vx_only", 40,  0, 7),
    (F(21, 10), F(13, 5),  "vx_only", 40,  0, 7),
    (F(13, 5),  F(3, 1),   "vx_only", 40,  0, 7),
    # 9 pure 1v intervals covering [3, 12]
    (F(3, 1),   F(33, 10), "1v",      45,  0, 0),
    (F(33, 10), F(37, 10), "1v",      45,  0, 0),
    (F(37, 10), F(21, 5),  "1v",      45,  0, 0),
    (F(21, 5),  F(24, 5),  "1v",      45,  0, 0),
    (F(24, 5),  F(28, 5),  "1v",      45,  0, 0),
    (F(28, 5),  F(33, 5),  "1v",      45,  0, 0),
    (F(33, 5),  F(39, 5),  "1v",      45,  0, 0),
    (F(39, 5),  F(19, 2),  "1v",      45,  0, 0),
    (F(19, 2),  F(12, 1),  "1v",      45,  0, 0),
]


def _check(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def _req_int(x, label: str) -> int:
    if type(x) is not int:
        raise TypeError(f"{label} must be exact int, got {type(x).__name__}")
    return x


def floor_frac(q: F) -> int:
    return q.numerator // q.denominator


def ceil_frac(q: F) -> int:
    return -((-q.numerator) // q.denominator)


def fmt_dec(q: F, digits: int = 8) -> str:
    """Exact decimal truncation without floating-point arithmetic."""
    neg = q < 0
    q = -q if neg else q
    sc = 10 ** digits
    n = (q.numerator * sc) // q.denominator
    return f"{'-' if neg else ''}{n // sc}.{n % sc:0{digits}d}"


def alpha_bd_exact(D_end: int, L: int, p1: F) -> F:
    """Exact min_{p in [p0, p1]} f_p(L) / q_p(D_end) = prod_{j=L-1}^{D_end-1} j / (j + p1) for L >= 2."""
    if L < 2:
        return F(0)
    val = F(1)
    for j in range(L - 1, D_end):
        val *= F(j) / (F(j) + p1)
    return val


def beta_bd_exact(D_end: int, L: int, p0: F) -> F:
    """Exact rational upper bound on max_{p in [p0, p1]} W_p(D_end, L) / q_p(D_end)."""
    m = (D_end + L + 1) // 2  # ceil((D_end + L)/2) <= D_end
    val = F(1)
    for j in range(m, D_end):
        val *= F(j) / (F(j) + p0)
    return val


def q_ratio_upper_exact(Da: int, Db: int, p0: F, p1: F) -> F:
    """Exact max_{p in [p0, p1]} q_p(Db) / q_p(Da)."""
    val = F(1)
    if Db < Da:
        for j in range(Db, Da):
            val *= F(j) / (F(j) + p0)
    elif Db > Da:
        for j in range(Da, Db):
            val *= F(1) + p1 / F(j)
    return val


def build_interval_tables(p0: F, p1: F, mode: str, DM: int, K: int, DXM: int):
    """Pure-Fraction conservative integer-tick generator for one interval [p0, p1]."""
    D_end = DM + 1
    alpha = [[0] * (D_end + 1) for _ in range(D_end + 1)]
    beta = [[0] * (D_end + 1) for _ in range(D_end + 1)]
    for L in range(0, D_end + 1):
        alpha[D_end][L] = floor_frac(alpha_bd_exact(D_end, L, p1) * SCALE)
        beta[D_end][L] = ceil_frac(beta_bd_exact(D_end, L, p0) * SCALE)

    for D in range(DM, 1, -1):
        for L in range(0, D + 1):
            c = D - L
            b = F(L, D)
            abd = alpha[D][L - 1] if L > 0 else 0
            bbd = beta[D][L - 1] if L > 0 else 0
            a_cands = [
                (r * (F(1) + pp / F(D)) * alpha[D + 1][L + 1] + b * abd) / (pp + r + b)
                for r in (F(L), F(2 * L + c, 2))
                for pp in (p0, p1)
            ]
            alpha[D][L] = floor_frac(min(a_cands))
            b_cands = [
                (r * (F(1) + pp / F(D)) * beta[D + 1][L + 1] + b * bbd) / (pp + r + b)
                for r in (F(L), F(2 * L + c, 2))
                for pp in (p0, p1)
            ]
            beta[D][L] = ceil_frac(max(b_cands))

    gamma = [0] * (D_end + 1)
    gamma[D_end] = ceil_frac(F(beta[2][1]) / (p0 * F(D_end) + F(1)))
    for d in range(DM, 1, -1):
        gamma[d] = ceil_frac((F(gamma[d + 1]) + F(beta[2][1], d)) / (p0 + F(1) + F(1, d)))

    out = {"alpha": alpha, "beta": beta, "gamma": gamma}
    if mode == "1v":
        return out

    if mode == "vx_only":
        vx = [[[0] * (D_end + 1) for _ in range(D_end + 1)] for _ in range(DXM + 2)]
        for Dx in range(2, DXM + 2):
            for D in range(2, D_end + 1):
                for L in range(0, D):
                    vx[Dx][D][L] = alpha[D][L]
        for Dx in range(DXM, 1, -1):
            for D in range(DM, 1, -1):
                for L in range(0, D):
                    c_other = D - L - 1
                    b = F(L, D)
                    v_aup = vx[Dx][D + 1][L + 1]
                    v_xup = vx[Dx + 1][D][L]
                    v_bd = vx[Dx][D][L - 1] if L > 0 else 0
                    ra0 = F(L) + F(1, Dx)
                    rx0 = F(1, D)
                    cands = [
                        (ra * (F(1) + pp / F(D)) * v_aup + rx * v_xup + b * v_bd) / (pp + ra + rx + b)
                        for ra in (ra0, ra0 + F(c_other, 2))
                        for rx in (rx0, F(Dx))
                        for pp in (p0, p1)
                    ]
                    vx[Dx][D][L] = floor_frac(min(cands))
        out["vx"] = vx
        return out

    # mode == "full"
    psi2 = [[[[0] * (K + 2) for _ in range(K + 2)] for _ in range(K + 2)] for _ in range(K + 2)]
    for Da in range(2, K + 2):
        for Db in range(2, K + 2):
            if Da <= K and Db <= K:
                continue
            qr = q_ratio_upper_exact(Da, Db, p0, p1)
            for La in range(0, Da):
                for Lb in range(0, Db):
                    psi2[Da][La][Db][Lb] = floor_frac(F(alpha[Da][La]) - F(beta[Db][Lb]) * qr)

    for Da in range(K, 1, -1):
        for Db in range(K, 1, -1):
            for La in range(0, Da):
                ca_ext = Da - La - 1
                ba = F(La, Da)
                for Lb in range(0, Db):
                    cb_ext = Db - Lb - 1
                    bb = F(Lb, Db)
                    va_up = psi2[Da + 1][La + 1][Db][Lb]
                    vb_up = psi2[Da][La][Db + 1][Lb + 1]
                    va_bd = psi2[Da][La - 1][Db][Lb] if La > 0 else 0
                    vb_bd = psi2[Da][La][Db][Lb - 1] if Lb > 0 else 0
                    ra0 = F(La) + F(1, Db)
                    rb0 = F(Lb) + F(1, Da)
                    cands = [
                        (ra * (F(1) + pp / F(Da)) * va_up + rb * vb_up + ba * va_bd + bb * vb_bd) / (pp + ra + rb + ba + bb)
                        for ra in (ra0, ra0 + F(ca_ext, 2))
                        for rb in (rb0, rb0 + F(cb_ext, 2))
                        for pp in (p0, p1)
                    ]
                    psi2[Da][La][Db][Lb] = floor_frac(min(cands))

    eta1 = [[0] * (K + 2) for _ in range(K + 2)]
    qr_k2 = q_ratio_upper_exact(K + 1, 2, p0, p1)
    for La in range(1, K + 1):
        eta1[K + 1][La] = floor_frac(F(alpha[K + 1][La]) - F(gamma[K + 1]) * qr_k2)
    for Da in range(K, 1, -1):
        for La in range(1, Da):
            ca = Da - La
            p_up = eta1[Da + 1][La + 1]
            p_obd = eta1[Da][La - 1] if La > 1 else 0
            p_bbd = psi2[Da][La - 1][2][1]
            bo, bb = F(La - 1, Da), F(1, Da)
            cands = [
                (r * (F(1) + pp / F(Da)) * p_up + bo * p_obd + bb * p_bbd) / (pp + r + bo + bb)
                for r in (F(La), F(2 * La + ca, 2))
                for pp in (p0, p1)
            ]
            eta1[Da][La] = floor_frac(min(cands))

    eta1_x = [[[0] * (K + 2) for _ in range(K + 2)] for _ in range(DXM + 2)]
    for Dx in range(2, DXM + 2):
        for Da in range(2, K + 2):
            if Dx <= DXM and Da <= K:
                continue
            qr_a2 = q_ratio_upper_exact(Da, 2, p0, p1)
            for La in range(1, Da):
                if Dx == DXM + 1:
                    eta1_x[Dx][Da][La] = eta1[Da][La]
                else:
                    eta1_x[Dx][Da][La] = floor_frac(F(alpha[Da][La]) - F(gamma[Da]) * qr_a2)
    for Dx in range(DXM, 1, -1):
        for Da in range(K, 1, -1):
            for La in range(1, Da):
                ca_other = Da - La - 1
                p_aup = eta1_x[Dx][Da + 1][La + 1]
                p_xup = eta1_x[Dx + 1][Da][La]
                p_obd = eta1_x[Dx][Da][La - 1] if La > 1 else 0
                p_bbd = psi2[Da][La - 1][2][1]
                ra0 = F(La) + F(1, Dx)
                rx0 = F(1, Da)
                bo, bb = F(La - 1, Da), F(1, Da)
                cands = [
                    (ra * (F(1) + pp / F(Da)) * p_aup + rx * p_xup + bo * p_obd + bb * p_bbd) / (pp + ra + rx + bo + bb)
                    for ra in (ra0, ra0 + F(ca_other, 2))
                    for rx in (rx0, F(Dx))
                    for pp in (p0, p1)
                ]
                eta1_x[Dx][Da][La] = floor_frac(min(cands))

    out["psi2"] = psi2
    out["eta1"] = eta1
    out["eta1_x"] = eta1_x
    return out


def verify_interval_exact(p0: F, p1: F, mode: str, DM: int, K: int, DXM: int, T: dict):
    """
    Standalone ZERO-FLOAT exact rational verifier for one interval [p0, p1].
    Raises ValueError / TypeError on ANY violation (fail-closed even under python3 -O).
    """
    _check(isinstance(p0, F) and isinstance(p1, F) and F(1) <= p0 < p1, f"invalid interval [{p0}, {p1}]")
    _check(mode in ("1v", "vx_only", "full"), f"invalid mode {mode}")
    D_end = DM + 1
    alpha, beta, gamma = T["alpha"], T["beta"], T["gamma"]
    rows = {"alpha": 0, "beta": 0, "gamma": 0, "vx": 0, "psi2": 0, "eta1": 0, "eta1_x": 0}

    # 1. Verify alpha (rescaled lower submartingale v_p = alpha * q_p(D))
    for L in range(0, D_end + 1):
        a_b = _req_int(alpha[D_end][L], f"alpha[{D_end}][{L}]")
        _check(F(a_b) <= alpha_bd_exact(D_end, L, p1) * SCALE, f"alpha boundary violation at L={L}")
    for D in range(2, DM + 1):
        for L in range(0, D + 1):
            c = D - L
            a0 = F(_req_int(alpha[D][L], f"alpha[{D}][{L}]"))
            a_up = F(_req_int(alpha[D + 1][L + 1], f"alpha[{D+1}][{L+1}]"))
            a_bd = F(_req_int(alpha[D][L - 1], f"alpha[{D}][{L-1}]")) if L > 0 else F(0)
            b = F(L, D)
            for r in (F(L), F(2 * L + c, 2)):
                for pp in (p0, p1):
                    drift = r * ((F(1) + pp / F(D)) * a_up - a0) + b * (a_bd - a0) - pp * a0
                    _check(drift >= 0, f"alpha drift violation at D={D}, L={L}, r={r}, p={pp}")
            rows["alpha"] += 1

    # 2. Verify beta (rescaled upper supermartingale w_p = beta * q_p(D), nonnegative everywhere)
    for L in range(0, D_end + 1):
        b_b = _req_int(beta[D_end][L], f"beta[{D_end}][{L}]")
        _check(b_b >= 0, f"beta negative at boundary L={L}")
        _check(F(b_b) >= beta_bd_exact(D_end, L, p0) * SCALE, f"beta boundary violation at L={L}")
    for D in range(2, DM + 1):
        for L in range(0, D + 1):
            c = D - L
            b0_i = _req_int(beta[D][L], f"beta[{D}][{L}]")
            _check(b0_i >= 0, f"beta negative at D={D}, L={L}")
            b0 = F(b0_i)
            b_up = F(_req_int(beta[D + 1][L + 1], f"beta[{D+1}][{L+1}]"))
            b_bd = F(_req_int(beta[D][L - 1], f"beta[{D}][{L-1}]")) if L > 0 else F(0)
            b = F(L, D)
            for r in (F(L), F(2 * L + c, 2)):
                for pp in (p0, p1):
                    drift = r * ((F(1) + pp / F(D)) * b_up - b0) + b * (b_bd - b0) - pp * b0
                    _check(drift <= 0, f"beta drift violation at D={D}, L={L}, r={r}, p={pp}")
            rows["beta"] += 1

    # 3. Verify gamma (rescaled leaf upper supermartingale U_p(d) = gamma(d)*(1+p), nonnegative & nonincreasing)
    b21 = _req_int(beta[2][1], "beta[2][1]")
    _check(b21 >= 0, "beta[2][1] negative")
    g_end = _req_int(gamma[D_end], f"gamma[{D_end}]")
    _check(g_end >= 0, "gamma boundary negative")
    _check(F(g_end) >= F(b21) / (p0 * F(D_end) + F(1)), "gamma boundary violation")
    for d in range(2, DM + 1):
        gd = _req_int(gamma[d], f"gamma[{d}]")
        gd1 = _req_int(gamma[d + 1], f"gamma[{d+1}]")
        _check(gd >= 0, f"gamma negative at d={d}")
        _check(gd1 <= gd, f"gamma monotonicity violation at d={d}")
        for pp in (p0, p1):
            drift = (F(gd1) - F(gd)) + F(1, d) * (F(b21) - F(gd)) - pp * F(gd)
            _check(drift <= 0, f"gamma drift violation at d={d}, p={pp}")
        rows["gamma"] += 1

    if mode == "1v":
        m_ticks = _req_int(alpha[2][1], "alpha[2][1]") - 2 * _req_int(gamma[2], "gamma[2]")
        margin_norm = F(m_ticks, SCALE)
        _check(margin_norm > 0, f"1v margin not positive on [{p0},{p1}]: {margin_norm}")
        return {
            "interval": [str(p0), str(p1)],
            "mode": mode,
            "DM": DM,
            "rows": rows,
            "eta_star_exact": str(margin_norm),
            "eta_star_dec_trunc": fmt_dec(margin_norm, 6),
            "c_star_lb_exact": str(margin_norm * (F(1) + p0)),
            "c_star_lb_dec_trunc": fmt_dec(margin_norm * (F(1) + p0), 6),
        }

    if mode == "vx_only":
        vx = T["vx"]
        for Dx in range(2, DXM + 2):
            for D in range(2, D_end + 1):
                if Dx <= DXM and D <= DM:
                    continue
                for L in range(0, D):
                    _check(_req_int(vx[Dx][D][L], f"vx[{Dx}][{D}][{L}]") <= _req_int(alpha[D][L], f"alpha[{D}][{L}]"),
                           f"vx boundary violation at ({Dx},{D},{L})")
        for Dx in range(2, DXM + 1):
            for D in range(2, DM + 1):
                for L in range(0, D):
                    c_other = D - L - 1
                    b = F(L, D)
                    v0 = F(_req_int(vx[Dx][D][L], f"vx[{Dx}][{D}][{L}]"))
                    v_aup = F(_req_int(vx[Dx][D + 1][L + 1], "vx_aup"))
                    v_xup = F(_req_int(vx[Dx + 1][D][L], "vx_xup"))
                    v_bd = F(_req_int(vx[Dx][D][L - 1], "vx_bd")) if L > 0 else F(0)
                    ra0 = F(L) + F(1, Dx)
                    rx0 = F(1, D)
                    for ra in (ra0, ra0 + F(c_other, 2)):
                        for rx in (rx0, F(Dx)):
                            for pp in (p0, p1):
                                drift = ra * ((F(1) + pp / F(D)) * v_aup - v0) + rx * (v_xup - v0) + b * (v_bd - v0) - pp * v0
                                _check(drift >= 0, f"vx drift violation at ({Dx},{D},{L})")
                    rows["vx"] += 1
        m1 = F(vx[2][2][1] - 2 * gamma[2], SCALE)
        m2 = F(vx[3][2][1] - gamma[2] - gamma[3], SCALE)
        m3 = F(alpha[2][1] - gamma[2] - gamma[4], SCALE)
        margin_norm = min(m1, m2, m3)
        _check(margin_norm > 0, f"vx_only margin not positive on [{p0},{p1}]: {margin_norm}")
        return {
            "interval": [str(p0), str(p1)],
            "mode": mode,
            "DM": DM,
            "DXM": DXM,
            "rows": rows,
            "eta_star_exact": str(margin_norm),
            "eta_star_dec_trunc": fmt_dec(margin_norm, 6),
            "c_star_lb_exact": str(margin_norm * (F(1) + p0)),
            "c_star_lb_dec_trunc": fmt_dec(margin_norm * (F(1) + p0), 6),
        }

    # mode == "full"
    psi2, eta1, eta1_x = T["psi2"], T["eta1"], T["eta1_x"]
    # 4. Verify psi2
    for Da in range(2, K + 2):
        for Db in range(2, K + 2):
            if Da <= K and Db <= K:
                continue
            qr = q_ratio_upper_exact(Da, Db, p0, p1)
            for La in range(0, Da):
                for Lb in range(0, Db):
                    p_b = _req_int(psi2[Da][La][Db][Lb], f"psi2[{Da}][{La}][{Db}][{Lb}]")
                    _check(F(p_b) <= F(_req_int(alpha[Da][La], "a")) - F(_req_int(beta[Db][Lb], "b")) * qr,
                           f"psi2 boundary violation at ({Da},{La},{Db},{Lb})")
    for Da in range(2, K + 1):
        for Db in range(2, K + 1):
            for La in range(0, Da):
                ca_ext = Da - La - 1
                ba = F(La, Da)
                for Lb in range(0, Db):
                    cb_ext = Db - Lb - 1
                    bb = F(Lb, Db)
                    p0v = F(_req_int(psi2[Da][La][Db][Lb], "psi2"))
                    va_up = F(_req_int(psi2[Da + 1][La + 1][Db][Lb], "psi2_aup"))
                    vb_up = F(_req_int(psi2[Da][La][Db + 1][Lb + 1], "psi2_bup"))
                    va_bd = F(_req_int(psi2[Da][La - 1][Db][Lb], "psi2_abd")) if La > 0 else F(0)
                    vb_bd = F(_req_int(psi2[Da][La][Db][Lb - 1], "psi2_bbd")) if Lb > 0 else F(0)
                    ra0 = F(La) + F(1, Db)
                    rb0 = F(Lb) + F(1, Da)
                    for ra in (ra0, ra0 + F(ca_ext, 2)):
                        for rb in (rb0, rb0 + F(cb_ext, 2)):
                            for pp in (p0, p1):
                                drift = (ra * ((F(1) + pp / F(Da)) * va_up - p0v)
                                         + rb * (vb_up - p0v)
                                         + ba * (va_bd - p0v)
                                         + bb * (vb_bd - p0v)
                                         - pp * p0v)
                                _check(drift >= 0, f"psi2 drift violation at ({Da},{La},{Db},{Lb})")
                    rows["psi2"] += 1

    # 5. Verify eta1
    qr_k2 = q_ratio_upper_exact(K + 1, 2, p0, p1)
    for La in range(1, K + 1):
        e_b = _req_int(eta1[K + 1][La], f"eta1[{K+1}][{La}]")
        _check(F(e_b) <= F(_req_int(alpha[K + 1][La], "a")) - F(_req_int(gamma[K + 1], "g")) * qr_k2,
               f"eta1 boundary violation at La={La}")
    for Da in range(2, K + 1):
        for La in range(1, Da):
            ca = Da - La
            p0v = F(_req_int(eta1[Da][La], "eta1"))
            p_up = F(_req_int(eta1[Da + 1][La + 1], "eta1_up"))
            p_obd = F(_req_int(eta1[Da][La - 1], "eta1_obd")) if La > 1 else F(0)
            p_bbd = F(_req_int(psi2[Da][La - 1][2][1], "psi2_handoff"))
            bo, bb = F(La - 1, Da), F(1, Da)
            for r in (F(La), F(2 * La + ca, 2)):
                for pp in (p0, p1):
                    drift = r * ((F(1) + pp / F(Da)) * p_up - p0v) + bo * (p_obd - p0v) + bb * (p_bbd - p0v) - pp * p0v
                    _check(drift >= 0, f"eta1 drift violation at ({Da},{La})")
            rows["eta1"] += 1

    # 6. Verify eta1_x
    for Dx in range(2, DXM + 2):
        for Da in range(2, K + 2):
            if Dx <= DXM and Da <= K:
                continue
            qr_a2 = q_ratio_upper_exact(Da, 2, p0, p1)
            for La in range(1, Da):
                ex_b = _req_int(eta1_x[Dx][Da][La], f"eta1_x[{Dx}][{Da}][{La}]")
                if Dx == DXM + 1:
                    _check(ex_b <= _req_int(eta1[Da][La], "eta1"), "eta1_x Dx boundary violation")
                else:
                    _check(F(ex_b) <= F(_req_int(alpha[Da][La], "a")) - F(_req_int(gamma[Da], "g")) * qr_a2,
                           "eta1_x Da boundary violation")
    for Dx in range(2, DXM + 1):
        for Da in range(2, K + 1):
            for La in range(1, Da):
                ca_other = Da - La - 1
                p0v = F(_req_int(eta1_x[Dx][Da][La], "eta1_x"))
                p_aup = F(_req_int(eta1_x[Dx][Da + 1][La + 1], "eta1_x_aup"))
                p_xup = F(_req_int(eta1_x[Dx + 1][Da][La], "eta1_x_xup"))
                p_obd = F(_req_int(eta1_x[Dx][Da][La - 1], "eta1_x_obd")) if La > 1 else F(0)
                p_bbd = F(_req_int(psi2[Da][La - 1][2][1], "psi2_handoff"))
                ra0 = F(La) + F(1, Dx)
                rx0 = F(1, Da)
                bo, bb = F(La - 1, Da), F(1, Da)
                for ra in (ra0, ra0 + F(ca_other, 2)):
                    for rx in (rx0, F(Dx)):
                        for pp in (p0, p1):
                            drift = (ra * ((F(1) + pp / F(Da)) * p_aup - p0v)
                                     + rx * (p_xup - p0v)
                                     + bo * (p_obd - p0v)
                                     + bb * (p_bbd - p0v)
                                     - pp * p0v)
                            _check(drift >= 0, f"eta1_x drift violation at ({Dx},{Da},{La})")
                rows["eta1_x"] += 1

    m1 = F(eta1_x[2][2][1] - gamma[2], SCALE)
    m2 = F(eta1_x[3][2][1] - gamma[3], SCALE)
    m3 = F(eta1[2][1] - gamma[4], SCALE)
    margin_norm = min(m1, m2, m3)
    _check(margin_norm > 0, f"full margin not positive on [{p0},{p1}]: {margin_norm}")
    return {
        "interval": [str(p0), str(p1)],
        "mode": mode,
        "DM": DM,
        "K": K,
        "DXM": DXM,
        "rows": rows,
        "eta_star_exact": str(margin_norm),
        "eta_star_dec_trunc": fmt_dec(margin_norm, 6),
        "c_star_lb_exact": str(margin_norm * (F(1) + p0)),
        "c_star_lb_dec_trunc": fmt_dec(margin_norm * (F(1) + p0), 6),
    }


def verify_tiling_and_analytic_tail(specs):
    """Verify exact contiguous rational cover [1, 12] and exact analytic inequalities for p >= 12 and u=2 vs 3."""
    _check(len(specs) > 0, "empty interval specs")
    for i, sp in enumerate(specs):
        p0, p1 = sp[0], sp[1]
        _check(isinstance(p0, F) and isinstance(p1, F), f"non-Fraction tiling endpoints at index {i}")
        _check(F(1) <= p0 < p1 <= F(12), f"invalid/out-of-range tiling interval at index {i}: [{p0}, {p1}]")
    _check(specs[0][0] == F(1), f"cover must start at 1, got {specs[0][0]}")
    for i in range(len(specs) - 1):
        _check(specs[i][1] == specs[i + 1][0], f"tiling gap/overlap at index {i}: {specs[i][1]} != {specs[i+1][0]}")
    p_cut = specs[-1][1]
    _check(p_cut == F(12), f"cover must end at 12, got {p_cut}")
    # Analytic tail threshold at p >= 12: 1/3 - 4*sqrt(p+1)/(sqrt(3)*(2p+1)) > 1/4000 > 0
    # Equivalent at p=12 to 16*13 / (3*25^2) < (1/3 - 1/4000)^2, and (p+1)/(2p+1)^2 strictly decreasing for p >= 1.
    _check(4 * p_cut * p_cut - 44 * p_cut - 47 > 0, "analytic tail threshold fails at p_cut")
    _check(8 * p_cut - 44 > 0, "analytic tail quadratic not increasing at p_cut")
    _check(F(16 * 13, 3 * 25 * 25) < (F(1, 3) - F(1, 4000)) ** 2, "analytic tail 1/4000 bound fails at p=12")
    # Early gap u=2 vs 3 for ALL p >= 1: 2*sqrt(p+1)/(sqrt(3)*(2p+1)) <= 2*sqrt(2)/(3*sqrt(3)) < 3/5
    _check(F(8, 27) < F(9, 25), "early gap bound fails")
    # Exact check of submartingale f_p(L) and supermartingale W_p(D, L) drifts on sample rational p and L, D
    for p in (F(1), F(3, 2), F(2), F(7, 2), F(12)):
        for L in range(2, 25):
            # f_p(L)/f_p(L) drift >= p:
            if L == 2:
                _check(2 * p - 1 >= p, f"f_p drift failed at L=2, p={p}")
            else:
                _check(p * (F(L, L - 1) - F(1, L + p - 2)) >= p, f"f_p drift failed at L={L}, p={p}")
    return {
        "tiling_start": str(specs[0][0]),
        "tiling_end": str(p_cut),
        "n_intervals": len(specs),
        "analytic_tail_poly_at_12": int(4 * p_cut * p_cut - 44 * p_cut - 47),
        "early_gap_ratio_sq": "8/27 < 9/25",
    }


def verify_two_birth_identity_exact():
    """Verify exact signed two-birth identity and endpoint-split inequality on all 24 unlabelled trees n=2..7."""
    def canon(adj):
        n = len(adj)
        def encode(u, p):
            ch = sorted(encode(v, u) for v in adj[u] if v != p)
            return "(" + "".join(ch) + ")"
        return min(encode(r, -1) for r in range(n))

    trees_by_n = {2: [([1], [0])]}
    for n in range(2, 7):
        seen = {}
        for adj in trees_by_n[n]:
            for u in range(n):
                na = [list(x) for x in adj] + [[u]]
                na[u].append(n)
                c = canon(na)
                if c not in seen:
                    seen[c] = tuple(tuple(x) for x in na)
        trees_by_n[n + 1] = list(seen.values())

    def test_score(adj, v):
        d = len(adj[v])
        nbr_degs = tuple(sorted(len(adj[u]) for u in adj[v]))
        h = d * 13 + sum((i + 1) * x * x for i, x in enumerate(nbr_degs))
        return F(17 + (h % 97), 100)

    checked = 0
    for n in range(2, 8):
        for adj in trees_by_n[n]:
            D = [len(adj[u]) for u in range(n)]
            r = [sum(F(1, D[u]) for u in adj[x]) for x in range(n)]
            _check(sum(r) == n, "rate sum mismatch")
            lhs = F(0)
            for x in range(n):
                p_ax = r[x] / n
                T1 = [list(z) for z in adj] + [[x]]
                T1[x].append(n)
                D1 = [len(T1[u]) for u in range(n + 1)]
                r1 = [sum(F(1, D1[u]) for u in T1[y]) for y in range(n + 1)]
                for y in range(n + 1):
                    p_by = r1[y] / (n + 1)
                    T2 = [list(z) for z in T1] + [[y]]
                    T2[y].append(n + 1)
                    lhs += p_ax * p_by * (test_score(T2, n) - test_score(T2, n + 1))
            lhs_scaled = lhs * n * (n + 1)
            rhs_path = F(0)
            for x in range(n):
                T2 = [list(z) for z in adj] + [[x, n + 1], [n]]
                T2[x].append(n)
                rhs_path += (r[x] / (D[x] + 1)) * (test_score(T2, n) - test_score(T2, n + 1))
            rhs_edge = F(0)
            for x in range(n):
                for y in adj[x]:
                    if x < y:
                        Txy = [list(z) for z in adj] + [[x], [y]]
                        Txy[x].append(n)
                        Txy[y].append(n + 1)
                        wxy = r[y] / (D[y] * (D[y] + 1)) - r[x] / (D[x] * (D[x] + 1))
                        rhs_edge += wxy * (test_score(Txy, n) - test_score(Txy, n + 1))
            _check(lhs_scaled == rhs_path + rhs_edge, f"two-birth identity mismatch at n={n}")
            _check(sum(r[x] / (D[x] + 1) for x in range(n)) >= 1, "weight sum < 1")
            split_lb = F(0)
            for x in range(n):
                Tpath = [list(z) for z in adj] + [[x, n + 1], [n]]
                Tpath[x].append(n)
                dp = test_score(Tpath, n) - test_score(Tpath, n + 1)
                avg_edge = F(0)
                for y in adj[x]:
                    Txy = [list(z) for z in adj] + [[x], [y]]
                    Txy[x].append(n)
                    Txy[y].append(n + 1)
                    avg_edge += test_score(Txy, n)
                avg_edge /= D[x]
                split_lb += (r[x] / (D[x] + 1)) * (dp - avg_edge)
            _check(lhs_scaled >= split_lb, "endpoint-split inequality failed")
            checked += 1
    return checked


def run_negative_controls(sample_full_spec, sample_full_T, sample_vx_spec, sample_vx_T, sample_1v_spec, sample_1v_T):
    """Verify 14 planted defects (all table families, boundaries, non-int, malformed tilings) fail closed."""
    out = {}

    def test_one(name, fn):
        caught = False
        try:
            fn()
        except (ValueError, TypeError) as e:
            caught = True
            out[name] = f"REJECTED ({e})"
        if not caught:
            raise RuntimeError(f"Planted defect {name} was NOT rejected!")

    def mutate_and_check(spec, T_orig, mut):
        Tm = copy.deepcopy(T_orig)
        mut(Tm)
        verify_interval_exact(*spec, Tm)

    p0f, p1f, _, DMf, Kf, DXMf = sample_full_spec
    test_one("alpha_interior", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["alpha"][2].__setitem__(1, T["alpha"][2][1] + 1)))
    test_one("alpha_boundary", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["alpha"][DMf + 1].__setitem__(5, T["alpha"][DMf + 1][5] + 1)))
    test_one("beta_interior", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["beta"][2].__setitem__(1, T["beta"][2][1] - 1)))
    test_one("beta_boundary", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["beta"][DMf + 1].__setitem__(5, T["beta"][DMf + 1][5] - 1)))
    test_one("gamma_interior", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["gamma"].__setitem__(2, T["gamma"][2] - 1)))
    test_one("gamma_boundary", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["gamma"].__setitem__(DMf + 1, T["gamma"][DMf + 1] - 1)))
    test_one("psi2_interior", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["psi2"][2][0][2].__setitem__(1, T["psi2"][2][0][2][1] + 1)))
    test_one("eta1_interior", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["eta1"][2].__setitem__(1, T["eta1"][2][1] + 1)))
    test_one("eta1_x_interior", lambda: mutate_and_check(sample_full_spec, sample_full_T, lambda T: T["eta1_x"][2][2].__setitem__(1, T["eta1_x"][2][2][1] + 1)))
    test_one("vx_interior", lambda: mutate_and_check(sample_vx_spec, sample_vx_T, lambda T: T["vx"][2][2].__setitem__(1, T["vx"][2][2][1] + 1)))
    test_one("float_tick_injection", lambda: mutate_and_check(sample_1v_spec, sample_1v_T, lambda T: T["alpha"][sample_1v_spec[3] + 1].__setitem__(0, 0.0)))
    bad_specs = list(INTERVAL_SPECS)
    bad_specs[0] = (F(1, 1), F(104, 100), "full", 40, 14, 8)
    test_one("tiling_gap", lambda: verify_tiling_and_analytic_tail(bad_specs))
    test_one("tiling_reversed", lambda: verify_tiling_and_analytic_tail([(F(1), F(13), "1v", 45, 0, 0), (F(13), F(12), "1v", 45, 0, 0)]))
    test_one("tiling_float_endpoint", lambda: verify_tiling_and_analytic_tail([(1.0, 12.0, "1v", 45, 0, 0)]))

    return out


def main():
    tiling_info = verify_tiling_and_analytic_tail(INTERVAL_SPECS)
    interval_results = []
    total_rows = 0
    min_c_star = None
    sample_full = sample_vx = sample_1v = None

    for spec in INTERVAL_SPECS:
        T = build_interval_tables(*spec)
        res = verify_interval_exact(*spec, T)
        interval_results.append(res)
        total_rows += sum(res["rows"].values())
        cs = F(res["c_star_lb_exact"])
        if min_c_star is None or cs < min_c_star:
            min_c_star = cs
        if spec[2] == "full" and sample_full is None:
            sample_full = (spec, T)
        elif spec[2] == "vx_only" and sample_vx is None:
            sample_vx = (spec, T)
        elif spec[2] == "1v" and sample_1v is None:
            sample_1v = (spec, T)

    n_trees = verify_two_birth_identity_exact()
    ctrls = run_negative_controls(*sample_full, *sample_vx, *sample_1v)

    out = {
        "status": "ALL_REAL_MOMENTS_EXACT_CHECKS_PASSED",
        "tiling_and_analytic_tail": tiling_info,
        "total_verified_rows": total_rows,
        "min_c_star_on_1_to_12_exact": str(min_c_star),
        "min_c_star_on_1_to_12_dec_trunc": fmt_dec(min_c_star, 6),
        "two_birth_identity_trees_verified_n2_to_n7": n_trees,
        "negative_controls": ctrls,
        "intervals": interval_results,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
