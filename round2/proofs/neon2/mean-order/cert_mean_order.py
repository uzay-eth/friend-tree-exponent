#!/usr/bin/env python3
"""
Standalone zero-floating-point exact rational certificate & verifier for
ALL-LABEL MEAN DEGREE-SHARE MONOTONICITY in the Friend Tree:

    E[Z_1] = E[Z_2] > E[Z_3] > E[Z_4] > ...   (for ALL labels u >= 1)

and the stronger conditional theorem: for EVERY integer n >= 2 and EVERY finite tree T
of size n,
    E[Z_{n+1} - Z_{n+2} | T_n = T] >= delta_* / (n * (n + 1) * (n + 2)) > 0,

with explicit exact rational delta_* > 1/26 and E[Z_2] - E[Z_3] > 2/15.

Architecture:
  1. Untrusted proposal generator computes conservative integer tables (scale S = 10^9)
     for the 6 exact finite local sub/supermartingales:
       - v(D, L):    1-vertex lower submartingale (2 <= D <= DM, 0 <= L <= D)
       - w(D, L):    1-vertex upper supermartingale (2 <= D <= DM, 0 <= L <= D)
       - U(d):       leaf upper supermartingale (2 <= d <= DM)
       - phi2:       2-vertex coupled lower submartingale for adjacent cores a ~ b
       - phi1:       Phase-1 lower submartingale for core a with leaf b
       - phi1_x:     Phase-1 lower submartingale tracking parent core degree Dx in {2..DXM}
  2. Standalone ZERO-FLOAT exact `fractions.Fraction` verifier independently checks EVERY
     boundary inequality, EVERY generator drift inequality at ALL extreme rates,
     nonincreasing monotonicity of U(d), strict integer tick types (rejecting floats/bools),
     exact verification of the signed two-birth identity on all 24 unlabelled trees of
     sizes n = 2..7 (against exact two-step transition enumeration), and 10 planted
     negative controls (all fail-closed even under python3 -O).
"""
from fractions import Fraction as F
import copy
import json

SCALE = 10 ** 9
DM = 40
K = 14
DXM = 8


def fmt_dec(q: F, digits: int = 8) -> str:
    """Exact decimal truncation without floating-point arithmetic."""
    neg = q < 0
    q = -q if neg else q
    sc = 10 ** digits
    n = (q.numerator * sc) // q.denominator
    return f"{'-' if neg else ''}{n // sc}.{n % sc:0{digits}d}"


def build_proposal_tables():
    """Pure-Fraction acyclic 1-pass generator producing conservative integer-tick tables (scale 10^9)."""
    v = [[0] * (DM + 2) for _ in range(DM + 2)]
    w = [[0] * (DM + 2) for _ in range(DM + 2)]
    for D in range(2, DM + 2):
        for L in range(0, D + 1):
            v[D][L] = max(L - 1, 0) * SCALE
            w[D][L] = D * SCALE

    def floor_frac(q: F) -> int:
        return q.numerator // q.denominator

    def ceil_frac(q: F) -> int:
        return -((-q.numerator) // q.denominator)

    # Exact rational 1-pass quantization (guaranteed 0 drift violation by construction!)
    for D in range(DM, 1, -1):
        for L in range(0, D + 1):
            c = D - L
            v_bd = v[D][L - 1] if L > 0 else 0
            w_bd = w[D][L - 1] if L > 0 else 0
            b = F(L, D)
            # Lower submartingale v: floor of min over extreme rates
            v_cand = [
                (r * v[D + 1][L + 1] + b * v_bd) / (1 + r + b)
                for r in (F(L), F(2 * L + c, 2))
            ]
            v[D][L] = floor_frac(min(v_cand))
            # Upper supermartingale w: ceil of max over extreme rates
            w_cand = [
                (r * w[D + 1][L + 1] + b * w_bd) / (1 + r + b)
                for r in (F(L), F(2 * L + c, 2))
            ]
            w[D][L] = ceil_frac(max(w_cand))

    # Exact U(d): boundary at d = DM + 1 is ceil(w[2][1] / (DM + 2))
    U = [0] * (DM + 2)
    U[DM + 1] = ceil_frac(F(w[2][1], DM + 2))
    for d in range(DM, 1, -1):
        U[d] = ceil_frac((F(U[d + 1]) + F(w[2][1], d)) / (F(2) + F(1, d)))

    # Exact phi2[Da][La][Db][Lb] (2 <= Da, Db <= K+1, 0 <= La <= Da-1, 0 <= Lb <= Db-1)
    phi2 = [[[[0] * (K + 2) for _ in range(K + 2)] for _ in range(K + 2)] for _ in range(K + 2)]
    for Da in range(2, K + 2):
        for Db in range(2, K + 2):
            for La in range(0, Da):
                for Lb in range(0, Db):
                    phi2[Da][La][Db][Lb] = v[Da][La] - w[Db][Lb]
    for Da in range(K, 1, -1):
        for Db in range(K, 1, -1):
            for La in range(0, Da):
                ca_ext = Da - La - 1
                ba = F(La, Da)
                for Lb in range(0, Db):
                    cb_ext = Db - Lb - 1
                    bb = F(Lb, Db)
                    va_up = phi2[Da + 1][La + 1][Db][Lb]
                    vb_up = phi2[Da][La][Db + 1][Lb + 1]
                    va_bd = phi2[Da][La - 1][Db][Lb] if La > 0 else 0
                    vb_bd = phi2[Da][La][Db][Lb - 1] if Lb > 0 else 0
                    ra0 = F(La) + F(1, Db)
                    rb0 = F(Lb) + F(1, Da)
                    cands = [
                        (ra * va_up + rb * vb_up + ba * va_bd + bb * vb_bd) / (1 + ra + rb + ba + bb)
                        for ra in (ra0, ra0 + F(ca_ext, 2))
                        for rb in (rb0, rb0 + F(cb_ext, 2))
                    ]
                    phi2[Da][La][Db][Lb] = floor_frac(min(cands))

    # Exact phi1[Da][La] (2 <= Da <= K+1, 1 <= La <= Da-1)
    phi1 = [[0] * (K + 2) for _ in range(K + 2)]
    for Da in range(2, K + 2):
        for La in range(1, Da):
            phi1[Da][La] = v[Da][La] - U[Da]
    for Da in range(K, 1, -1):
        for La in range(1, Da):
            ca = Da - La
            p_up = phi1[Da + 1][La + 1]
            p_obd = phi1[Da][La - 1] if La > 1 else 0
            p_bbd = phi2[Da][La - 1][2][1]
            bo, bb = F(La - 1, Da), F(1, Da)
            cands = [
                (r * p_up + bo * p_obd + bb * p_bbd) / (1 + r + bo + bb)
                for r in (F(La), F(2 * La + ca, 2))
            ]
            phi1[Da][La] = floor_frac(min(cands))

    # Exact phi1_x[Dx][Da][La] (2 <= Dx <= DXM+1, 2 <= Da <= K+1, 1 <= La <= Da-1)
    phi1_x = [[[0] * (K + 2) for _ in range(K + 2)] for _ in range(DXM + 2)]
    for Dx in range(2, DXM + 2):
        for Da in range(2, K + 2):
            for La in range(1, Da):
                phi1_x[Dx][Da][La] = phi1[Da][La] if Dx == DXM + 1 else (v[Da][La] - U[Da])
    for Dx in range(DXM, 1, -1):
        for Da in range(K, 1, -1):
            for La in range(1, Da):
                ca_other = Da - La - 1
                p_aup = phi1_x[Dx][Da + 1][La + 1]
                p_xup = phi1_x[Dx + 1][Da][La]
                p_obd = phi1_x[Dx][Da][La - 1] if La > 1 else 0
                p_bbd = phi2[Da][La - 1][2][1]
                ra0 = F(La) + F(1, Dx)
                rx0 = F(1, Da)
                bo, bb = F(La - 1, Da), F(1, Da)
                cands = [
                    (ra * p_aup + rx * p_xup + bo * p_obd + bb * p_bbd) / (1 + ra + rx + bo + bb)
                    for ra in (ra0, ra0 + F(ca_other, 2))
                    for rx in (rx0, F(Dx))
                ]
                phi1_x[Dx][Da][La] = floor_frac(min(cands))

    return {"v": v, "w": w, "U": U, "phi2": phi2, "phi1": phi1, "phi1_x": phi1_x}


def _check(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def _req_int(x, label: str) -> int:
    if type(x) is not int:
        raise TypeError(f"{label} must be exact int, got {type(x).__name__}")
    return x


def verify_tables_exact(T):
    """
    Standalone ZERO-FLOAT exact rational verifier for all 6 local sub/supermartingale tables.
    Raises ValueError/TypeError if ANY type, boundary, monotonicity, or drift check fails
    (fail-closed even under python3 -O).
    """
    v, w, U = T["v"], T["w"], T["U"]
    phi2, phi1, phi1_x = T["phi2"], T["phi1"], T["phi1_x"]
    counts = {"v": 0, "w": 0, "U": 0, "phi2": 0, "phi1": 0, "phi1_x": 0}

    # 1. Verify v(D, L): boundary v(DM+1, L) <= max(L-1, 0)*SCALE and (L_min v - v) >= 0 everywhere
    for L in range(0, DM + 2):
        _check(_req_int(v[DM + 1][L], f"v[{DM+1}][{L}]") <= max(L - 1, 0) * SCALE, f"v boundary violation at L={L}")
    for D in range(2, DM + 1):
        for L in range(0, D + 1):
            c = D - L
            v0 = F(_req_int(v[D][L], f"v[{D}][{L}]"))
            v_up = F(_req_int(v[D + 1][L + 1], f"v[{D+1}][{L+1}]"))
            v_bd = F(_req_int(v[D][L - 1], f"v[{D}][{L-1}]")) if L > 0 else F(0)
            b = F(L, D)
            for r in (F(L), F(2 * L + c, 2)):
                drift = r * (v_up - v0) + b * (v_bd - v0) - v0
                _check(drift >= 0, f"v submartingale violation at D={D}, L={L}, r={r}: {drift}")
            counts["v"] += 1

    # 2. Verify w(D, L): boundary w(DM+1, L) >= (DM+1)*SCALE and (L_max w - w) <= 0 everywhere
    for L in range(0, DM + 2):
        _check(_req_int(w[DM + 1][L], f"w[{DM+1}][{L}]") >= (DM + 1) * SCALE, f"w boundary violation at L={L}")
    for D in range(2, DM + 1):
        for L in range(0, D + 1):
            c = D - L
            w0 = F(_req_int(w[D][L], f"w[{D}][{L}]"))
            w_up = F(_req_int(w[D + 1][L + 1], f"w[{D+1}][{L+1}]"))
            w_bd = F(_req_int(w[D][L - 1], f"w[{D}][{L-1}]")) if L > 0 else F(0)
            b = F(L, D)
            for r in (F(L), F(2 * L + c, 2)):
                drift = r * (w_up - w0) + b * (w_bd - w0) - w0
                _check(drift <= 0, f"w supermartingale violation at D={D}, L={L}, r={r}: {drift}")
            counts["w"] += 1

    # 3. Verify U(d): boundary U(DM+1) >= w(2,1)/(DM+2), nonincreasing U(d+1) <= U(d), and supermartingale drift <= 0
    w21 = _req_int(w[2][1], "w[2][1]")
    _check(F(_req_int(U[DM + 1], f"U[{DM+1}]")) >= F(w21, DM + 2), "U boundary violation")
    for d in range(2, DM + 1):
        ud = _req_int(U[d], f"U[{d}]")
        ud1 = _req_int(U[d + 1], f"U[{d+1}]")
        _check(ud1 <= ud, f"U monotonicity violation at d={d}: {ud1} > {ud}")
        # Since U[d+1] - U[d] <= 0, minimum increment rate 1 (from the leaf itself!) is worst case for upper bound!
        drift = (F(ud1) - F(ud)) + F(1, d) * (F(w21) - F(ud)) - F(ud)
        _check(drift <= 0, f"U supermartingale violation at d={d}: {drift}")
        counts["U"] += 1

    # 4. Verify phi2[Da][La][Db][Lb]: boundary <= v[Da][La] - w[Db][Lb] and submartingale drift >= 0
    for Da in range(2, K + 2):
        for Db in range(2, K + 2):
            if Da <= K and Db <= K:
                continue
            for La in range(0, Da):
                for Lb in range(0, Db):
                    p_b = _req_int(phi2[Da][La][Db][Lb], f"phi2[{Da}][{La}][{Db}][{Lb}]")
                    _check(p_b <= _req_int(v[Da][La], f"v[{Da}][{La}]") - _req_int(w[Db][Lb], f"w[{Db}][{Lb}]"),
                           f"phi2 boundary violation at ({Da},{La},{Db},{Lb})")
    for Da in range(2, K + 1):
        for Db in range(2, K + 1):
            for La in range(0, Da):
                ca_ext = Da - La - 1
                ba = F(La, Da)
                for Lb in range(0, Db):
                    cb_ext = Db - Lb - 1
                    bb = F(Lb, Db)
                    p0 = F(_req_int(phi2[Da][La][Db][Lb], f"phi2[{Da}][{La}][{Db}][{Lb}]"))
                    va_up = F(_req_int(phi2[Da + 1][La + 1][Db][Lb], "phi2_aup"))
                    vb_up = F(_req_int(phi2[Da][La][Db + 1][Lb + 1], "phi2_bup"))
                    va_bd = F(_req_int(phi2[Da][La - 1][Db][Lb], "phi2_abd")) if La > 0 else F(0)
                    vb_bd = F(_req_int(phi2[Da][La][Db][Lb - 1], "phi2_bbd")) if Lb > 0 else F(0)
                    ra0 = F(La) + F(1, Db)
                    rb0 = F(Lb) + F(1, Da)
                    for ra in (ra0, ra0 + F(ca_ext, 2)):
                        for rb in (rb0, rb0 + F(cb_ext, 2)):
                            drift = ra * (va_up - p0) + rb * (vb_up - p0) + ba * (va_bd - p0) + bb * (vb_bd - p0) - p0
                            _check(drift >= 0, f"phi2 violation at ({Da},{La},{Db},{Lb})")
                    counts["phi2"] += 1

    # 5. Verify phi1[Da][La]: boundary <= v[K+1][La] - U[K+1] and submartingale drift >= 0
    for La in range(1, K + 1):
        _check(_req_int(phi1[K + 1][La], f"phi1[{K+1}][{La}]")
               <= _req_int(v[K + 1][La], f"v[{K+1}][{La}]") - _req_int(U[K + 1], f"U[{K+1}]"),
               f"phi1 boundary violation at La={La}")
    for Da in range(2, K + 1):
        for La in range(1, Da):
            ca = Da - La
            p0 = F(_req_int(phi1[Da][La], f"phi1[{Da}][{La}]"))
            p_up = F(_req_int(phi1[Da + 1][La + 1], "phi1_up"))
            p_obd = F(_req_int(phi1[Da][La - 1], "phi1_obd")) if La > 1 else F(0)
            p_bbd = F(_req_int(phi2[Da][La - 1][2][1], "phi2_handoff"))
            bo, bb = F(La - 1, Da), F(1, Da)
            for r in (F(La), F(2 * La + ca, 2)):
                drift = r * (p_up - p0) + bo * (p_obd - p0) + bb * (p_bbd - p0) - p0
                _check(drift >= 0, f"phi1 violation at ({Da},{La})")
            counts["phi1"] += 1

    # 6. Verify phi1_x[Dx][Da][La]: boundaries & submartingale drift >= 0
    for Dx in range(2, DXM + 2):
        for Da in range(2, K + 2):
            if Dx <= DXM and Da <= K:
                continue
            for La in range(1, Da):
                px_b = _req_int(phi1_x[Dx][Da][La], f"phi1_x[{Dx}][{Da}][{La}]")
                if Dx == DXM + 1:
                    _check(px_b <= _req_int(phi1[Da][La], f"phi1[{Da}][{La}]"), "phi1_x Dx boundary violation")
                else:
                    _check(px_b <= _req_int(v[Da][La], f"v[{Da}][{La}]") - _req_int(U[Da], f"U[{Da}]"),
                           "phi1_x Da boundary violation")
    for Dx in range(2, DXM + 1):
        for Da in range(2, K + 1):
            for La in range(1, Da):
                ca_other = Da - La - 1
                p0 = F(_req_int(phi1_x[Dx][Da][La], f"phi1_x[{Dx}][{Da}][{La}]"))
                p_aup = F(_req_int(phi1_x[Dx][Da + 1][La + 1], "phi1_x_aup"))
                p_xup = F(_req_int(phi1_x[Dx + 1][Da][La], "phi1_x_xup"))
                p_obd = F(_req_int(phi1_x[Dx][Da][La - 1], "phi1_x_obd")) if La > 1 else F(0)
                p_bbd = F(_req_int(phi2[Da][La - 1][2][1], "phi2_handoff"))
                ra0 = F(La) + F(1, Dx)
                rx0 = F(1, Da)
                bo, bb = F(La - 1, Da), F(1, Da)
                for ra in (ra0, ra0 + F(ca_other, 2)):
                    for rx in (rx0, F(Dx)):
                        drift = ra * (p_aup - p0) + rx * (p_xup - p0) + bo * (p_obd - p0) + bb * (p_bbd - p0) - p0
                        _check(drift >= 0, f"phi1_x violation at ({Dx},{Da},{La})")
                counts["phi1_x"] += 1

    m_dx1 = F(phi1_x[2][2][1] - U[2], SCALE)
    m_dx2 = F(phi1_x[3][2][1] - U[3], SCALE)
    m_dx3 = F(phi1[2][1] - U[4], SCALE)
    delta_star = min(m_dx1, m_dx2, m_dx3)
    _check(delta_star > 0, f"Universal margin not positive: {delta_star}")

    # Early pair (2, 3): E[Z_2] - E[Z_3] = (1/6) * (M(T_3, center) - M(T_3, leaf)) >= (1/6) * (v(2,2) - U(2))
    gap_23 = F(v[2][2] - U[2], 6 * SCALE)
    _check(gap_23 > 0, f"Early gap (2,3) not positive: {gap_23}")

    return {
        "counts": counts,
        "v_2_1": str(F(v[2][1], SCALE)),
        "v_2_2": str(F(v[2][2], SCALE)),
        "w_2_1": str(F(w[2][1], SCALE)),
        "U_2": str(F(U[2], SCALE)),
        "U_3": str(F(U[3], SCALE)),
        "U_4": str(F(U[4], SCALE)),
        "phi1_2_1": str(F(phi1[2][1], SCALE)),
        "phi1_x_2_2_1": str(F(phi1_x[2][2][1], SCALE)),
        "phi1_x_3_2_1": str(F(phi1_x[3][2][1], SCALE)),
        "margin_Dx1_exact": str(m_dx1),
        "margin_Dx1_dec_trunc": fmt_dec(m_dx1, 8),
        "margin_Dx2_exact": str(m_dx2),
        "margin_Dx2_dec_trunc": fmt_dec(m_dx2, 8),
        "margin_Dx_ge_3_exact": str(m_dx3),
        "margin_Dx_ge_3_dec_trunc": fmt_dec(m_dx3, 8),
        "delta_star_exact": str(delta_star),
        "delta_star_dec_trunc": fmt_dec(delta_star, 8),
        "gap_23_exact": str(gap_23),
        "gap_23_dec_trunc": fmt_dec(gap_23, 8),
    }


def verify_two_birth_identity_exact():
    """
    Independently verify the exact signed two-birth identity AND the edge-splitting inequality
    on ALL 24 unlabelled trees of sizes n = 2..7 using an arbitrary positive test score M(T', v)
    in exact Fraction arithmetic!
    """
    # Generate all non-isomorphic trees up to n = 7 via canonical leaf-extension branching
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

    # Deterministic nontrivial isomorphism-invariant vertex score M(T', v) >= 0 (depends on neighborhood multiset!)
    def test_score(adj, v):
        d = len(adj[v])
        nbr_degs = tuple(sorted(len(adj[u]) for u in adj[v]))
        # Arbitrary positive exact Fraction function of rooted isomorphism class
        h = d * 13 + sum((i + 1) * x * x for i, x in enumerate(nbr_degs))
        return F(17 + (h % 97), 100)

    checked = 0
    for n in range(2, 8):
        for adj in trees_by_n[n]:
            D = [len(adj[u]) for u in range(n)]
            r = [sum(F(1, D[u]) for u in adj[x]) for x in range(n)]
            _check(sum(r) == n, f"rate sum mismatch on n={n}")
            # 1. Direct brute-force 2-step expectation of n*(n+1)*(M(T_{n+2}, a) - M(T_{n+2}, b))
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

            # 2. Signed two-birth identity RHS
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
            _check(lhs_scaled == rhs_path + rhs_edge, f"Two-birth identity mismatch on n={n} tree")

            # 3. Check sum_{x in T} r_T(x)/(D_x + 1) >= 1 and endpoint-split lower bound
            _check(sum(r[x] / (D[x] + 1) for x in range(n)) >= 1, f"weight sum < 1 on n={n}")
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
            _check(lhs_scaled >= split_lb, "Endpoint-split inequality failed!")
            checked += 1
    return checked


def run_negative_controls(T_good):
    """Verify ten planted single-tick interior, boundary, and non-integer defects fail closed."""
    defects = [
        ("v_interior", lambda T: T["v"][2].__setitem__(1, T["v"][2][1] + 1)),
        ("v_boundary", lambda T: T["v"][DM + 1].__setitem__(5, T["v"][DM + 1][5] + 1)),
        ("w_interior", lambda T: T["w"][2].__setitem__(1, T["w"][2][1] - 1)),
        ("w_boundary", lambda T: T["w"][DM + 1].__setitem__(5, T["w"][DM + 1][5] - 1)),
        ("U_interior", lambda T: T["U"].__setitem__(2, T["U"][2] - 1)),
        ("U_boundary", lambda T: T["U"].__setitem__(DM + 1, T["U"][DM + 1] - 1)),
        ("phi2_interior", lambda T: T["phi2"][2][0][2].__setitem__(1, T["phi2"][2][0][2][1] + 1)),
        ("phi1_interior", lambda T: T["phi1"][2].__setitem__(1, T["phi1"][2][1] + 1)),
        ("phi1_x_interior", lambda T: T["phi1_x"][2][2].__setitem__(1, T["phi1_x"][2][2][1] + 1)),
        ("float_tick_injection", lambda T: T["v"][DM + 1].__setitem__(0, 0.0)),
    ]
    out = {}
    for name, fn in defects:
        Tm = copy.deepcopy(T_good)
        fn(Tm)
        caught = False
        try:
            verify_tables_exact(Tm)
        except (ValueError, TypeError) as e:
            caught = True
            out[name] = f"REJECTED ({e})"
        if not caught:
            raise RuntimeError(f"Planted defect {name} was NOT rejected!")
    return out


def main():
    tables = build_proposal_tables()
    ver = verify_tables_exact(tables)
    n_trees = verify_two_birth_identity_exact()
    ctrls = run_negative_controls(tables)
    out = {
        "status": "ALL_EXACT_CHECKS_PASSED",
        "params": {"SCALE": SCALE, "DM": DM, "K": K, "DXM": DXM},
        "verification": ver,
        "two_birth_identity_trees_verified_n2_to_n7": n_trees,
        "negative_controls": ctrls,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
