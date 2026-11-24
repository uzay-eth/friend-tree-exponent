#!/usr/bin/env python3
"""
Standalone ZERO-FLOATING-POINT exact rational/integer certificate & verifier for
FRACTIONAL-MOMENT, ALL-LABEL STRICT MONOTONICITY in the Friend Tree
(Open Question 3 for ALL real p in [1/2, 1], uniting with Claim 9 for ALL p in [1/2, infinity)):

  Theorem 1 (Universal per-tree ordering, EVERY real p in [1/2, 1], EVERY n >= 2, EVERY tree T):
      E[Z_{n+1}^p - Z_{n+2}^p | T_n = T] >= (Gamma(n+2)*Gamma(p+1) / (n*(n+1)*Gamma(n+2+p))) * c_0(p) > 0.

  Theorem 2 (All-label strict monotonicity on edge-seeded T_2 = {1, 2}, EVERY real p in [1/2, 1]):
      E[Z_1^p] = E[Z_2^p] > E[Z_3^p] > E[Z_4^p] > ...     for EVERY real p in [1/2, 1],
      with early gap E[Z_2^p] - E[Z_3^p] > 2 / (5 * (p + 2)).

Architecture:
  - Partition [1/2, 1] into 14 contiguous closed rational intervals [p0, p1].
  - On each interval, rescale all core potentials by the exact rising-Gamma polynomial
    q_p(D) = prod_{j=1}^{D-1} (1 + p/j), with:
      1. Analytic core lower tail a_p^tail(L) = (L/(L+2)) * q_p(L) (exact submartingale for ALL D >= 12,
         0 <= L <= D, p in [0, 1]), handing off to alpha[D][L] on 2 <= D <= DM;
      2. Analytic core upper tail W_p(D, L) = Gamma((D+L)/2 + p) / (Gamma((D+L)/2)*Gamma(p+1)),
         handing off to beta[D][L] on 2 <= D <= DM;
      3. Exact predictable-hazard leaf tail beta[2][1]/(p0*(DM+1)+1), handing off to 2D leaf
         supersolution Q[D][L] (2 <= D <= DM, 1 <= L <= D) and 1D tail leaf supersolution gamma[d];
      4. Coupled 2-vertex post-promotion submartingale psi2[Da][La][Db][Lb] (2 <= Da, Db <= K);
      5. Pre-promotion core-leaf difference eta1[Da][La] (2 <= Da <= K);
      6. Parent-degree-tracked layer etax[Dx][Da][La] (3 <= Dx <= DXM) for core x (D_T(x) >= 2), plus
         bare-parent layer etax2_bare[Da][La] (Dx = 2, Lx = 0) for leaf x (D_T(x) = 1 in any tree of size n >= 3;
         at n = 2, T_2 = {1, 2} has identically zero edge-pairing term so eta1[2][1] > 0 directly!).
  - Exact arbitrary-precision integer cross-multiplication (ZERO floats, ZERO Gamma approximations,
    fail-closed even under python3 -O), plus concrete physical tree checks on all 24 unlabelled trees
    n = 2..7 and 14 planted negative controls.
"""
from fractions import Fraction as F
import copy
import json
import math

SCALE = 10 ** 9

# 14 contiguous exact rational intervals tiling [1/2, 1]: (p0, p1, DM, K, DXM)
INTERVAL_SPECS = [
    (F(1, 2),    F(101, 200), 100, 19, 8),
    (F(101, 200), F(51, 100), 100, 19, 8),
    (F(51, 100), F(103, 200), 100, 19, 8),
    (F(103, 200), F(13, 25),  100, 18, 8),
    (F(13, 25),  F(53, 100),  100, 18, 8),
    (F(53, 100), F(27, 50),   100, 18, 8),
    (F(27, 50),  F(11, 20),   100, 18, 8),
    (F(11, 20),  F(29, 50),    90, 18, 8),
    (F(29, 50),  F(31, 50),    90, 18, 8),
    (F(31, 50),  F(67, 100),   80, 16, 8),
    (F(67, 100), F(73, 100),   80, 16, 8),
    (F(73, 100), F(4, 5),      80, 16, 8),
    (F(4, 5),    F(9, 10),     80, 16, 8),
    (F(9, 10),   F(1, 1),      80, 16, 8),
]


def _check(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def _req_int(x, label: str) -> int:
    if type(x) is not int:
        raise TypeError(f"{label} must be exact int, got {type(x).__name__}")
    return x


def floor_div(num: int, den: int) -> int:
    _check(den > 0, "floor_div requires positive denominator")
    return num // den


def ceil_div(num: int, den: int) -> int:
    _check(den > 0, "ceil_div requires positive denominator")
    return -((-num) // den)


def floor_frac(q: F) -> int:
    return q.numerator // q.denominator


def ceil_frac(q: F) -> int:
    return -((-q.numerator) // q.denominator)


def fmt_dec(q: F, digits: int = 8) -> str:
    """Exact outward-truncated decimal formatting without floats."""
    neg = q < 0
    q = -q if neg else q
    sc = 10 ** digits
    n = (q.numerator * sc) // q.denominator
    return f"{'-' if neg else ''}{n // sc}.{n % sc:0{digits}d}"


def alpha_bd_exact(D_end: int, L: int, p1: F) -> F:
    """Exact min_{p in [p0, p1]} a_p^tail(L) / q_p(D_end) = (L/(L+2)) * prod_{j=L}^{D_end-1} j/(j+p1)."""
    if L < 1:
        return F(0)
    val = F(L, L + 2)
    for j in range(L, D_end):
        val *= F(j) / (F(j) + p1)
    return val


def beta_bd_exact(D_end: int, L: int, p0: F) -> F:
    """Exact rational upper bound on max_{p in [p0, p1]} W_p(D_end, L) / q_p(D_end)."""
    m = (D_end + L + 1) // 2
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


def _p_common(p0: F, p1: F):
    qp = math.lcm(p0.denominator, p1.denominator)
    pn0 = p0.numerator * (qp // p0.denominator)
    pn1 = p1.numerator * (qp // p1.denominator)
    return pn0, pn1, qp


def build_interval_tables(p0: F, p1: F, DM: int, K: int, DXM: int):
    """Exact pure-integer conservative generator for one interval [p0, p1] subset [1/2, 1]."""
    D_end = DM + 1
    pn0, pn1, qp = _p_common(p0, p1)
    pns = (pn0, pn1)

    alpha = [[0] * (D_end + 1) for _ in range(D_end + 1)]
    beta = [[0] * (D_end + 1) for _ in range(D_end + 1)]
    for L in range(0, D_end + 1):
        alpha[D_end][L] = floor_frac(alpha_bd_exact(D_end, L, p1) * SCALE)
        beta[D_end][L] = ceil_frac(beta_bd_exact(D_end, L, p0) * SCALE)

    # 1. alpha and beta (multiply drift by 2 * D * qp)
    for D in range(DM, 1, -1):
        for L in range(0, D + 1):
            rns = (2 * L, D + L)  # 2*r in {2L, D+L}
            abd = alpha[D][L - 1] if L > 0 else 0
            bbd = beta[D][L - 1] if L > 0 else 0
            a_up = alpha[D + 1][L + 1]
            b_up = beta[D + 1][L + 1]
            a_min = None
            b_max = None
            for rn in rns:
                for pn in pns:
                    den = 2 * D * pn + rn * D * qp + 2 * L * qp
                    a_num = rn * (D * qp + pn) * a_up + 2 * L * qp * abd
                    b_num = rn * (D * qp + pn) * b_up + 2 * L * qp * bbd
                    av = floor_div(a_num, den)
                    bv = ceil_div(b_num, den)
                    if a_min is None or av < a_min:
                        a_min = av
                    if b_max is None or bv > b_max:
                        b_max = bv
            alpha[D][L] = a_min
            beta[D][L] = b_max

    # 2. Leaf upper bounds: 2D Q[D][L] and 1D gamma[d] (in units of q_p(2), worst-case p = p0)
    b21 = beta[2][1]
    q_bd = ceil_frac(F(b21) / (p0 * F(D_end) + F(1)))
    Q = [[q_bd] * (D_end + 1) for _ in range(D_end + 1)]
    for D in range(DM, 1, -1):
        Q[D][0] = 0
        for L in range(1, D + 1):
            rns = (2 * L, D + L)
            q_up = Q[D + 1][L + 1]
            q_o = Q[D][L - 1] if L > 1 else 0
            q_mx = None
            for rn in rns:
                # Multiply drift by 2 * D * qp:
                den = 2 * D * pn0 + rn * D * qp + 2 * L * qp
                num = rn * D * qp * q_up + 2 * (L - 1) * qp * q_o + 2 * qp * b21
                qv = ceil_div(num, den)
                if q_mx is None or qv > q_mx:
                    q_mx = qv
            Q[D][L] = q_mx

    gamma = [0] * (D_end + 1)
    gamma[D_end] = q_bd
    for d in range(DM, 1, -1):
        gamma[d] = ceil_frac((F(gamma[d + 1]) + F(b21, d)) / (p0 + F(1) + F(1, d)))

    # 3. Coupled post-promotion 2-vertex psi2[Da][La][Db][Lb] (units of q_p(Da))
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
                ca_e = Da - La - 1
                ra_n0 = 2 * La * Db + 2
                ra_n1 = ra_n0 + ca_e * Db
                ra_ns = (ra_n0,) if ca_e == 0 else (ra_n0, ra_n1)
                for Lb in range(0, Db):
                    cb_e = Db - Lb - 1
                    rb_n0 = 2 * Lb * Da + 2
                    rb_n1 = rb_n0 + cb_e * Da
                    rb_ns = (rb_n0,) if cb_e == 0 else (rb_n0, rb_n1)
                    va_u = psi2[Da + 1][La + 1][Db][Lb]
                    vb_u = psi2[Da][La][Db + 1][Lb + 1]
                    va_b = psi2[Da][La - 1][Db][Lb] if La > 0 else 0
                    vb_b = psi2[Da][La][Db][Lb - 1] if Lb > 0 else 0
                    c_bd = 2 * La * Db * qp * va_b + 2 * Lb * Da * qp * vb_b
                    d_bd = 2 * La * Db * qp + 2 * Lb * Da * qp
                    best = None
                    for ra_n in ra_ns:
                        t_ra_d = ra_n * Da * qp
                        for rb_n in rb_ns:
                            t_rb_d = rb_n * Db * qp
                            t_rb_u = t_rb_d * vb_u
                            for pn in pns:
                                den = 2 * Da * Db * pn + t_ra_d + t_rb_d + d_bd
                                num = ra_n * (Da * qp + pn) * va_u + t_rb_u + c_bd
                                val = floor_div(num, den)
                                if best is None or val < best:
                                    best = val
                    psi2[Da][La][Db][Lb] = best

    # 4. Pre-promotion eta1[Da][La] (units of q_p(Da))
    eta1 = [[0] * (K + 2) for _ in range(K + 2)]
    qr_k2 = q_ratio_upper_exact(K + 1, 2, p0, p1)
    for La in range(1, K + 1):
        eta1[K + 1][La] = floor_frac(F(alpha[K + 1][La]) - F(Q[K + 1][La]) * qr_k2)
    for Da in range(K, 1, -1):
        for La in range(1, Da):
            rns = (2 * La, Da + La)
            p_up = eta1[Da + 1][La + 1]
            p_ob = eta1[Da][La - 1] if La > 1 else 0
            p_bb = psi2[Da][La - 1][2][1]
            best = None
            for rn in rns:
                for pn in pns:
                    den = 2 * Da * pn + rn * Da * qp + 2 * La * qp
                    num = rn * (Da * qp + pn) * p_up + 2 * (La - 1) * qp * p_ob + 2 * qp * p_bb
                    val = floor_div(num, den)
                    if best is None or val < best:
                        best = val
            eta1[Da][La] = best

    # 5. Parent-tracked etax[Dx][Da][La] for core x (3 <= Dx <= DXM, handoff at DXM+1 to eta1)
    etax = [[[eta1[Da][La] for La in range(K + 2)] for Da in range(K + 2)] for _ in range(DXM + 2)]
    for Dx in range(DXM, 2, -1):
        for Da in range(K, 1, -1):
            for La in range(1, Da):
                ca_o = Da - La - 1
                ra_n0 = 2 * La * Dx + 2
                ra_n1 = ra_n0 + ca_o * Dx
                ra_ns = (ra_n0,) if ca_o == 0 else (ra_n0, ra_n1)
                rx_ns = (1, Dx * Da)  # rx * Da in {1, Dx * Da}
                v_au = etax[Dx][Da + 1][La + 1]
                v_xu = etax[Dx + 1][Da][La]
                v_ob = etax[Dx][Da][La - 1] if La > 1 else 0
                v_bb = psi2[Da][La - 1][2][1]
                # Multiply drift by 2 * Da * Dx * qp:
                c_b = 2 * (La - 1) * Dx * qp * v_ob + 2 * Dx * qp * v_bb
                d_b = 2 * La * Dx * qp
                best = None
                for ra_n in ra_ns:
                    t_ra_d = ra_n * Da * qp
                    for rx_n in rx_ns:
                        t_rx_d = 2 * rx_n * Dx * qp
                        t_rx_u = t_rx_d * v_xu
                        for pn in pns:
                            den = 2 * Da * Dx * pn + t_ra_d + t_rx_d + d_b
                            num = ra_n * (Da * qp + pn) * v_au + t_rx_u + c_b
                            val = floor_div(num, den)
                            if best is None or val < best:
                                best = val
                etax[Dx][Da][La] = best

    # 6. Bare-parent layer etax2_bare[Da][La] for leaf x in any tree T of size n >= 3:
    #    In T_{x-a-b}, x starts at (Dx=2, Lx=0) with core neighbors a (deg Da) and y (deg Dy >= 2),
    #    so ra in [La + 1/2, La + 1/2 + (Da-La-1)/2] and rx = 1/Da + 1/Dy in [1/Da, 1/Da + 1/2];
    #    on x's first leaf gain, Dx jumps 2 -> 3, handing off to etax[3][Da][La]!
    etax2_bare = [[eta1[Da][La] for La in range(K + 2)] for Da in range(K + 2)]
    for Da in range(K, 1, -1):
        for La in range(1, Da):
            ca_o = Da - La - 1
            ra_ns = (2 * La + 1, 2 * La + 1 + ca_o)  # 2 * ra
            rx_ns = (2, 2 + Da)                      # 2 * Da * rx in {2, 2 + Da}
            v_au = etax2_bare[Da + 1][La + 1]
            v_xu = etax[3][Da][La]
            v_ob = etax2_bare[Da][La - 1] if La > 1 else 0
            v_bb = psi2[Da][La - 1][2][1]
            # Multiply drift by 2 * Da * qp:
            c_b = 2 * (La - 1) * qp * v_ob + 2 * qp * v_bb
            d_b = 2 * La * qp
            best = None
            for ra_n in ra_ns:
                t_ra_d = ra_n * Da * qp
                for rx_n in rx_ns:
                    t_rx_u = rx_n * qp * v_xu
                    for pn in pns:
                        den = 2 * Da * pn + t_ra_d + rx_n * qp + d_b
                        num = ra_n * (Da * qp + pn) * v_au + t_rx_u + c_b
                        val = floor_div(num, den)
                        if best is None or val < best:
                            best = val
            etax2_bare[Da][La] = best

    return {
        "alpha": alpha,
        "beta": beta,
        "Q": Q,
        "gamma": gamma,
        "psi2": psi2,
        "eta1": eta1,
        "etax": etax,
        "etax2_bare": etax2_bare,
    }


def verify_interval_exact(p0: F, p1: F, DM: int, K: int, DXM: int, tbl: dict):
    """Standalone exact zero-float verifier for one interval [p0, p1] subset [1/2, 1]."""
    _check(isinstance(p0, F) and isinstance(p1, F), "p0, p1 must be exact Fractions")
    _check(F(1, 2) <= p0 < p1 <= F(1, 1), f"interval [{p0}, {p1}] outside [1/2, 1]")
    _req_int(DM, "DM"); _req_int(K, "K"); _req_int(DXM, "DXM")
    _check(DM >= 12 and 2 <= K < DM and 3 <= DXM < DM, "invalid DM/K/DXM parameters")

    D_end = DM + 1
    pn0, pn1, qp = _p_common(p0, p1)
    pns = (pn0, pn1)
    alpha = tbl["alpha"]; beta = tbl["beta"]; Q = tbl["Q"]; gamma = tbl["gamma"]
    psi2 = tbl["psi2"]; eta1 = tbl["eta1"]; etax = tbl["etax"]; etax2_bare = tbl["etax2_bare"]
    rows = 0

    # 1. Boundaries & interior of alpha and beta
    for L in range(0, D_end + 1):
        a_b = _req_int(alpha[D_end][L], f"alpha[{D_end}][{L}]")
        b_b = _req_int(beta[D_end][L], f"beta[{D_end}][{L}]")
        _check(F(a_b) <= alpha_bd_exact(D_end, L, p1) * SCALE, f"alpha boundary fail L={L}")
        _check(F(b_b) >= beta_bd_exact(D_end, L, p0) * SCALE, f"beta boundary fail L={L}")
        _check(b_b >= 0, f"negative beta boundary L={L}")

    for D in range(2, DM + 1):
        for L in range(0, D + 1):
            a0 = _req_int(alpha[D][L], f"alpha[{D}][{L}]")
            b0 = _req_int(beta[D][L], f"beta[{D}][{L}]")
            _check(b0 >= 0, f"negative beta[{D}][{L}]")
            a_up = _req_int(alpha[D + 1][L + 1], "a_up")
            b_up = _req_int(beta[D + 1][L + 1], "b_up")
            a_bd = _req_int(alpha[D][L - 1], "a_bd") if L > 0 else 0
            b_bd = _req_int(beta[D][L - 1], "b_bd") if L > 0 else 0
            for rn in (2 * L, D + L):
                for pn in pns:
                    den = 2 * D * pn + rn * D * qp + 2 * L * qp
                    _check(rn * (D * qp + pn) * a_up + 2 * L * qp * a_bd >= den * a0,
                           f"alpha drift fail D={D}, L={L}")
                    _check(rn * (D * qp + pn) * b_up + 2 * L * qp * b_bd <= den * b0,
                           f"beta drift fail D={D}, L={L}")
            rows += 2

    # 2. Leaf upper bounds Q[D][L] and gamma[d]
    b21 = _req_int(beta[2][1], "beta[2][1]")
    q_bd_req = F(b21) / (p0 * F(D_end) + F(1))
    for L in range(1, D_end + 1):
        ql = _req_int(Q[D_end][L], f"Q[{D_end}][{L}]")
        _check(ql >= 0 and F(ql) >= q_bd_req, f"Q boundary fail L={L}")
    for D in range(2, DM + 1):
        for L in range(1, D + 1):
            q0 = _req_int(Q[D][L], f"Q[{D}][{L}]")
            _check(q0 >= 0, f"negative Q[{D}][{L}]")
            q_up = _req_int(Q[D + 1][L + 1], "q_up")
            q_o = _req_int(Q[D][L - 1], "q_o") if L > 1 else 0
            for rn in (2 * L, D + L):
                den = 2 * D * pn0 + rn * D * qp + 2 * L * qp
                num = rn * D * qp * q_up + 2 * (L - 1) * qp * q_o + 2 * qp * b21
                _check(num <= den * q0, f"Q drift fail D={D}, L={L}")
            rows += 1

    g_end = _req_int(gamma[D_end], f"gamma[{D_end}]")
    _check(g_end >= 0 and F(g_end) >= q_bd_req, "gamma boundary fail")
    for d in range(2, DM + 1):
        gd = _req_int(gamma[d], f"gamma[{d}]")
        gd1 = _req_int(gamma[d + 1], f"gamma[{d+1}]")
        _check(gd >= 0, f"negative gamma[{d}]")
        _check(gd >= gd1, f"gamma not nonincreasing at d={d}")
        _check((F(gd1) - F(gd)) + F(b21 - gd, d) - p0 * F(gd) <= 0, f"gamma drift fail d={d}")
        rows += 1

    # 3. psi2 boundaries & interior
    for Da in range(2, K + 2):
        for Db in range(2, K + 2):
            if Da <= K and Db <= K:
                continue
            qr = q_ratio_upper_exact(Da, Db, p0, p1)
            for La in range(0, Da):
                for Lb in range(0, Db):
                    pv = _req_int(psi2[Da][La][Db][Lb], "psi2_bd")
                    _check(F(pv) <= F(alpha[Da][La]) - F(beta[Db][Lb]) * qr,
                           f"psi2 boundary fail ({Da},{La},{Db},{Lb})")
                    rows += 1

    for Da in range(2, K + 1):
        for Db in range(2, K + 1):
            for La in range(0, Da):
                ca_e = Da - La - 1
                ra_ns = (2 * La * Db + 2,) if ca_e == 0 else (2 * La * Db + 2, 2 * La * Db + 2 + ca_e * Db)
                for Lb in range(0, Db):
                    cb_e = Db - Lb - 1
                    rb_ns = (2 * Lb * Da + 2,) if cb_e == 0 else (2 * Lb * Da + 2, 2 * Lb * Da + 2 + cb_e * Da)
                    v0 = _req_int(psi2[Da][La][Db][Lb], "psi2_v0")
                    va_u = _req_int(psi2[Da + 1][La + 1][Db][Lb], "psi2_au")
                    vb_u = _req_int(psi2[Da][La][Db + 1][Lb + 1], "psi2_bu")
                    va_b = _req_int(psi2[Da][La - 1][Db][Lb], "psi2_ab") if La > 0 else 0
                    vb_b = _req_int(psi2[Da][La][Db][Lb - 1], "psi2_bb") if Lb > 0 else 0
                    c_bd = 2 * La * Db * qp * va_b + 2 * Lb * Da * qp * vb_b
                    d_bd = 2 * La * Db * qp + 2 * Lb * Da * qp
                    for ra_n in ra_ns:
                        t_ra_d = ra_n * Da * qp
                        for rb_n in rb_ns:
                            t_rb_d = rb_n * Db * qp
                            t_rb_u = t_rb_d * vb_u
                            for pn in pns:
                                den = 2 * Da * Db * pn + t_ra_d + t_rb_d + d_bd
                                num = ra_n * (Da * qp + pn) * va_u + t_rb_u + c_bd
                                _check(num >= den * v0, f"psi2 drift fail ({Da},{La},{Db},{Lb})")
                    rows += 1

    # 4. eta1 boundaries & interior
    qr_k2 = q_ratio_upper_exact(K + 1, 2, p0, p1)
    for La in range(1, K + 1):
        ev = _req_int(eta1[K + 1][La], f"eta1[{K+1}][{La}]")
        _check(F(ev) <= F(alpha[K + 1][La]) - F(Q[K + 1][La]) * qr_k2, f"eta1 boundary fail La={La}")
        rows += 1
    for Da in range(2, K + 1):
        for La in range(1, Da):
            v0 = _req_int(eta1[Da][La], f"eta1[{Da}][{La}]")
            p_up = _req_int(eta1[Da + 1][La + 1], "eta1_up")
            p_ob = _req_int(eta1[Da][La - 1], "eta1_ob") if La > 1 else 0
            p_bb = _req_int(psi2[Da][La - 1][2][1], "eta1_bb")
            for rn in (2 * La, Da + La):
                for pn in pns:
                    den = 2 * Da * pn + rn * Da * qp + 2 * La * qp
                    num = rn * (Da * qp + pn) * p_up + 2 * (La - 1) * qp * p_ob + 2 * qp * p_bb
                    _check(num >= den * v0, f"eta1 drift fail Da={Da}, La={La}")
            rows += 1

    # 5. etax boundaries & interior (3 <= Dx <= DXM)
    for Dx in range(3, DXM + 2):
        for Da in range(2, K + 2):
            for La in range(1, Da):
                if Dx <= DXM and Da <= K:
                    continue
                ev = _req_int(etax[Dx][Da][La], "etax_bd")
                _check(ev <= _req_int(eta1[Da][La], "eta1_ref"), f"etax boundary fail ({Dx},{Da},{La})")
                rows += 1
    for Dx in range(3, DXM + 1):
        for Da in range(2, K + 1):
            for La in range(1, Da):
                v0 = _req_int(etax[Dx][Da][La], "etax_v0")
                v_au = _req_int(etax[Dx][Da + 1][La + 1], "etax_au")
                v_xu = _req_int(etax[Dx + 1][Da][La], "etax_xu")
                v_ob = _req_int(etax[Dx][Da][La - 1], "etax_ob") if La > 1 else 0
                v_bb = _req_int(psi2[Da][La - 1][2][1], "etax_bb")
                ca_o = Da - La - 1
                ra_ns = (2 * La * Dx + 2,) if ca_o == 0 else (2 * La * Dx + 2, 2 * La * Dx + 2 + ca_o * Dx)
                c_b = 2 * (La - 1) * Dx * qp * v_ob + 2 * Dx * qp * v_bb
                d_b = 2 * La * Dx * qp
                for ra_n in ra_ns:
                    t_ra_d = ra_n * Da * qp
                    for rx_n in (1, Dx * Da):
                        t_rx_d = 2 * rx_n * Dx * qp
                        t_rx_u = t_rx_d * v_xu
                        for pn in pns:
                            den = 2 * Da * Dx * pn + t_ra_d + t_rx_d + d_b
                            num = ra_n * (Da * qp + pn) * v_au + t_rx_u + c_b
                            _check(num >= den * v0, f"etax drift fail ({Dx},{Da},{La})")
                rows += 1

    # 6. etax2_bare boundary & interior (Dx=2, Lx=0 for leaf x in size n >= 3 trees)
    for La in range(1, K + 1):
        ev = _req_int(etax2_bare[K + 1][La], f"etax2_bare[{K+1}][{La}]")
        _check(ev <= _req_int(eta1[K + 1][La], "eta1_k1"), f"etax2_bare boundary fail La={La}")
        rows += 1
    for Da in range(2, K + 1):
        for La in range(1, Da):
            v0 = _req_int(etax2_bare[Da][La], f"etax2_bare[{Da}][{La}]")
            v_au = _req_int(etax2_bare[Da + 1][La + 1], "bare_au")
            v_xu = _req_int(etax[3][Da][La], "bare_xu")
            v_ob = _req_int(etax2_bare[Da][La - 1], "bare_ob") if La > 1 else 0
            v_bb = _req_int(psi2[Da][La - 1][2][1], "bare_bb")
            ca_o = Da - La - 1
            c_b = 2 * (La - 1) * qp * v_ob + 2 * qp * v_bb
            d_b = 2 * La * qp
            for ra_n in (2 * La + 1, 2 * La + 1 + ca_o):
                t_ra_d = ra_n * Da * qp
                for rx_n in (2, 2 + Da):
                    t_rx_u = rx_n * qp * v_xu
                    for pn in pns:
                        den = 2 * Da * pn + t_ra_d + rx_n * qp + d_b
                        num = ra_n * (Da * qp + pn) * v_au + t_rx_u + c_b
                        _check(num >= den * v0, f"etax2_bare drift fail Da={Da}, La={La}")
            rows += 1

    # 7. Universal per-vertex margin over ALL D_x >= 1 (n >= 3) AND T_2 (n = 2) AND T_3 early gap!
    #    - n = 2 (T_2 = {1,2}): edge term is identically 0, margin is eta1[2][1] * q_p(2) / SCALE!
    #    - n >= 3, leaf x (D_T(x) = 1 -> D_new = 2, bare in T_{x-a-b}): margin is (etax2_bare[2][1] - max_L Q[2][L]) * q_p(2) / SCALE
    #    - n >= 3, core x (2 <= D_T(x) <= DXM-1 -> 3 <= D_new <= DXM): (etax[D_new][2][1] - max_L Q[D_new][L]) * q_p(2) / SCALE
    #    - n >= 3, core x (D_T(x) >= DXM -> D_new >= DXM+1): (eta1[2][1] - gamma[DXM+1]) * q_p(2) / SCALE
    qmax = [0] * (DXM + 1)
    for d in range(2, DXM + 1):
        qmax[d] = max(_req_int(Q[d][L], f"Q[{d}][{L}]") for L in range(1, d + 1))

    m_t2 = _req_int(eta1[2][1], "eta1[2][1]")
    m_leaf = _req_int(etax2_bare[2][1], "etax2_bare[2][1]") - qmax[2]
    m_core_fin = min(_req_int(etax[d][2][1], f"etax[{d}][2][1]") - qmax[d] for d in range(3, DXM + 1))
    m_core_tail = _req_int(eta1[2][1], "eta1[2][1]") - _req_int(gamma[DXM + 1], f"gamma[{DXM+1}]")
    min_scaled_ticks = min(m_t2, m_leaf, m_core_fin, m_core_tail)
    _check(min_scaled_ticks > 0, f"nonpositive margin on [{p0}, {p1}]: ticks={min_scaled_ticks}")

    # Multiply by min_{p in [p0, p1]} q_p(2) = 1 + p0!
    c_interval = (F(1) + p0) * F(min_scaled_ticks, SCALE)
    # T_3 early gap: E[Z_2^p] - E[Z_3^p] >= (alpha[2][2] - Q[2][2]) / (SCALE * (p + 2))
    t3_ticks = _req_int(alpha[2][2], "alpha[2][2]") - _req_int(Q[2][2], "Q[2][2]")
    _check(t3_ticks > 0, f"nonpositive T_3 early gap on [{p0}, {p1}]")
    t3_coeff = F(t3_ticks, SCALE)

    return {
        "p0": str(p0),
        "p1": str(p1),
        "DM": DM,
        "K": K,
        "DXM": DXM,
        "rows_verified": rows,
        "c_interval_exact": str(c_interval),
        "c_interval_dec": fmt_dec(c_interval, 8),
        "t3_coeff_dec": fmt_dec(t3_coeff, 6),
        "_c_frac": c_interval,
        "_t3_frac": t3_coeff,
    }


def verify_tiling(specs):
    """Verify exact contiguous tiling of [1/2, 1] with strict Fraction endpoints."""
    _check(isinstance(specs, (list, tuple)) and len(specs) >= 1, "specs must be nonempty sequence")
    for idx, item in enumerate(specs):
        _check(isinstance(item, (list, tuple)) and len(item) == 5, f"malformed spec {idx}")
        p0, p1, DM, K, DXM = item
        _check(isinstance(p0, F) and isinstance(p1, F), f"spec {idx} endpoints must be exact Fractions")
        _check(p0 < p1, f"spec {idx} must satisfy p0 < p1")
        _req_int(DM, "DM"); _req_int(K, "K"); _req_int(DXM, "DXM")
    _check(specs[0][0] == F(1, 2), f"tiling must start at 1/2, got {specs[0][0]}")
    _check(specs[-1][1] == F(1, 1), f"tiling must end at 1, got {specs[-1][1]}")
    for idx in range(len(specs) - 1):
        _check(specs[idx][1] == specs[idx + 1][0], f"tiling discontinuity at {idx}")


def verify_concrete_trees_and_tail():
    """
    Independently verify:
      1. Exact rational tail drift of a_p^tail(L) = (L/(L+2)) * q_p(L) across D in {12, 20, 80, 100},
         all 0 <= L <= D, and p in {1/2, 3/4, 1}.
      2. Exact signed two-birth pairing identity (Lemma 3), zero edge-term at n=2 (T_2), bare-leaf
         neighbor property D_T(y) >= 2 for every leaf x in every tree of size n = 3..7, and
         sum_x r_T(x)/(D_x+1) >= 1 on ALL 24 unlabelled trees n = 2..7.
    """
    tail_checks = 0
    for p in (F(1, 2), F(3, 4), F(1, 1)):
        for D in (12, 20, 80, 100):
            for L in range(1, D + 1):
                gain = F(2) * (L + p) / ((L + 2) * (L + 3))
                loss = F(L, D) * (F(L, L + 2) - (F(L - 1, L + 1) * F(L - 1, L - 1 + p) if L > 1 else F(0)))
                _check(gain - loss > 0, f"tail drift failed p={p}, D={D}, L={L}")
                tail_checks += 1

    # Enumerate all 24 unlabelled trees n = 2..7
    def canon(adj):
        n = len(adj)
        def enc(u, par):
            return "(" + "".join(sorted(enc(v, u) for v in adj[u] if v != par)) + ")"
        return min(enc(r, -1) for r in range(n))

    trees = {2: [([1], [0])]}
    for n in range(2, 7):
        seen = {}
        for adj in trees[n]:
            for u in range(n):
                na = [list(x) for x in adj] + [[u]]
                na[u].append(n)
                c = canon(na)
                if c not in seen:
                    seen[c] = tuple(tuple(x) for x in na)
        trees[n + 1] = list(seen.values())

    checked = 0
    for n in range(2, 8):
        for adj in trees[n]:
            deg = [len(adj[u]) for u in range(n)]
            rT = [sum(F(1, deg[v]) for v in adj[u]) for u in range(n)]
            _check(sum(rT[x] / (deg[x] + 1) for x in range(n)) >= 1, f"weight sum < 1 at n={n}")
            if n == 2:
                _check(rT[1] / (deg[1] * (deg[1] + 1)) - rT[0] / (deg[0] * (deg[0] + 1)) == 0, "T_2 edge term != 0")
            else:
                for x in range(n):
                    if deg[x] == 1:
                        y = adj[x][0]
                        _check(deg[y] >= 2, f"leaf neighbor has deg < 2 at n={n}")
            # Exact two-birth coefficient check for unordered edges
            for x in range(n):
                for y in adj[x]:
                    if x < y:
                        p_xy = rT[x] * (rT[y] - F(1, deg[x] * (deg[x] + 1)))
                        p_yx = rT[y] * (rT[x] - F(1, deg[y] * (deg[y] + 1)))
                        _check(p_xy - p_yx == rT[y] / (deg[y] * (deg[y] + 1)) - rT[x] / (deg[x] * (deg[x] + 1)),
                               "pairing identity mismatch")
            checked += 1
    return {"tail_sample_checks": tail_checks, "unlabelled_trees_n2_to_n7": checked}


def run_negative_controls(spec0, tbl0):
    """Verify 14 planted defects across every layer, type check, and tiling fail closed."""
    p0, p1, DM, K, DXM = spec0
    out = {}

    def expect_fail(name, fn):
        caught = False
        try:
            fn()
        except (ValueError, TypeError) as e:
            caught = True
            out[name] = f"REJECTED ({e})"
        if not caught:
            raise RuntimeError(f"Planted defect {name} was NOT rejected!")

    def mut(fn):
        tm = copy.deepcopy(tbl0)
        fn(tm)
        verify_interval_exact(p0, p1, DM, K, DXM, tm)

    expect_fail("alpha_boundary_tick", lambda: mut(lambda t: t["alpha"][DM + 1].__setitem__(5, t["alpha"][DM + 1][5] + 1)))
    expect_fail("beta_boundary_tick", lambda: mut(lambda t: t["beta"][DM + 1].__setitem__(5, t["beta"][DM + 1][5] - 1)))
    expect_fail("alpha_interior_tick", lambda: mut(lambda t: t["alpha"][2].__setitem__(1, t["alpha"][2][1] + 1)))
    expect_fail("beta_interior_tick", lambda: mut(lambda t: t["beta"][2].__setitem__(1, t["beta"][2][1] - 1)))
    expect_fail("Q_boundary_tick", lambda: mut(lambda t: t["Q"][DM + 1].__setitem__(1, t["Q"][DM + 1][1] - 1)))
    expect_fail("Q_interior_tick", lambda: mut(lambda t: t["Q"][2].__setitem__(1, t["Q"][2][1] - 1)))
    expect_fail("gamma_interior_tick", lambda: mut(lambda t: t["gamma"].__setitem__(5, t["gamma"][5] - 1)))
    expect_fail("psi2_interior_tick", lambda: mut(lambda t: t["psi2"][2][1][2].__setitem__(1, t["psi2"][2][1][2][1] + 1)))
    expect_fail("eta1_interior_tick", lambda: mut(lambda t: t["eta1"][2].__setitem__(1, t["eta1"][2][1] + 1)))
    expect_fail("etax_interior_tick", lambda: mut(lambda t: t["etax"][3][2].__setitem__(1, t["etax"][3][2][1] + 1)))
    expect_fail("etax2_bare_interior_tick", lambda: mut(lambda t: t["etax2_bare"][2].__setitem__(1, t["etax2_bare"][2][1] + 1)))
    expect_fail("float_tick_rejected", lambda: mut(lambda t: t["alpha"][2].__setitem__(1, float(t["alpha"][2][1]))))
    expect_fail("bool_tick_rejected", lambda: mut(lambda t: t["gamma"].__setitem__(2, True)))
    bad_tiling = list(INTERVAL_SPECS)
    bad_tiling[1] = (F(506, 1000), F(51, 100), 100, 19, 8)
    expect_fail("tiling_gap_rejected", lambda: verify_tiling(bad_tiling))
    return out


def main():
    verify_tiling(INTERVAL_SPECS)
    results = []
    total_rows = 0
    min_c = None
    min_t3 = None
    first_tbl = None

    for idx, (p0, p1, DM, K, DXM) in enumerate(INTERVAL_SPECS):
        tbl = build_interval_tables(p0, p1, DM, K, DXM)
        res = verify_interval_exact(p0, p1, DM, K, DXM, tbl)
        if idx == 0:
            first_tbl = tbl
        total_rows += res["rows_verified"]
        c_f = res.pop("_c_frac")
        t3_f = res.pop("_t3_frac")
        if min_c is None or c_f < min_c:
            min_c = c_f
        if min_t3 is None or t3_f < min_t3:
            min_t3 = t3_f
        results.append(res)

    _check(min_c > F(1, 2500), f"min_c={min_c} <= 1/2500")
    _check(min_t3 > F(2, 5), f"min_t3={min_t3} <= 2/5")
    tree_res = verify_concrete_trees_and_tail()
    ctrls = run_negative_controls(INTERVAL_SPECS[0], first_tbl)

    out = {
        "status": "FRACTIONAL_MOMENTS_HALF_TO_ONE_EXACT_CHECKS_PASSED",
        "scale": SCALE,
        "intervals_count": len(results),
        "total_rows_verified": total_rows,
        "min_c_on_half_to_one_exact": str(min_c),
        "min_c_on_half_to_one_dec": fmt_dec(min_c, 8),
        "certifies_min_c_gt_1_over_2500": min_c > F(1, 2500),
        "min_t3_coeff_exact": str(min_t3),
        "min_t3_coeff_dec": fmt_dec(min_t3, 6),
        "certifies_t3_gt_2_over_5": min_t3 > F(2, 5),
        "concrete_checks": tree_res,
        "intervals": results,
        "negative_controls": ctrls,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
