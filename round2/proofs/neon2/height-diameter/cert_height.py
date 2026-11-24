#!/usr/bin/env python3
"""
Standalone ZERO-FLOATING-POINT exact rational certificate & verifier for
LOGARITHMIC HEIGHT & DIAMETER BOUNDS IN THE RANDOM FRIEND TREE:

  Theorem 1 (Certified Almost-Sure Height & Diameter Enclosures — Zero Astra Dependencies):
      Let (T_n)_{n >= 2} be the edge-seeded random friend tree, H(T_n) = max_{v in T_n} dist(1, v)
      its height, and diam(T_n) = max_{u,v in T_n} dist(u, v) its diameter. Then almost surely:

          5000/5067 (0.986777)  <=  liminf H(T_n)/log n    <=  limsup H(T_n)/log n    <=  1.462248,
         10000/5067 (1.973554)  <=  liminf diam(T_n)/log n <=  limsup diam(T_n)/log n <=  2.924496,

      with explicit non-asymptotic finite-n exponential upper tail for every n >= 3 and every c > c_hi:
          P(H_core(T_n) >= c * log n) <= (Phi_z(T_3) / f_min) * n^{-(c * log z - lambda)}.

Architecture:
  1. Upper bound (pure discrete-time depth-tilted 1-vertex supermartingale):
     For z = 31/10, lam = 165439/100000, DM = 100, scale S = 10^12, integer-tick table
     f[D][L] (2 <= D <= DM, 0 <= L <= D, extended by f = S for D >= DM+1) satisfies:
       - f[D][L] >= 1 (positive lower bound f_min > 0),
       - z * f[2][1] <= lam * S (exact far-tail inequality for ALL D >= DM+1 simultaneously),
       - max_{r in {L, L+c/2}} [ r*(f[D+1][L+1] - f[D][L]) + (L/D)*(f[D][L-1] - f[D][L] + z*f[2][1]) ] <= lam*f[D][L]
         at all 5,148 interior states.
     Exact rational positive-term series log z > 2 * sum_{k=0}^{25} ((z-1)/(z+1))^{2k+1}/(2k+1) certifies c_hi > lam/log z.
  2. Lower bound (pathwise-coupled reflected (L <= 50, j <= 50) CMJ branching process + Biggins 1977):
     Integer-tick sub/super-solutions u_lo <= phi(s) <= u_hi and v_lo <= -phi'(s) <= v_hi at
     s = 0 (proving phi(0) = E[N] > 12 > 1, supercriticality) and at two bracket points
     theta_1 = 865/1000, theta_2 = 867/1000 with beta_hi = 10134/10000 certify:
       - h'(theta_1) <= 0 <= h'(theta_2) for convex h(theta) = log phi(theta) + beta_hi * theta,
       - h(theta) > 0 everywhere on [0, infinity) via exact alternating Taylor upper bound on e^{-x},
       hence beta_* = sup_{theta>0} (-log phi(theta))/theta <= beta_hi and speed 1/beta_* >= 5000/5067.
  3. Concrete physical generator check on all 23 in-domain unlabelled trees n = 3..7 (confirming
     L(Phi_z)(T) <= lam * Phi_z(T) directly by enumerating every physical directed arrow), plus
     11 planted single-tick, far-tail, boundary, non-int, and convexity-violation negative controls (fail-closed
     even under python3 -O).
"""
from fractions import Fraction as F
import copy
import json

SCALE = 10 ** 12

# Upper-bound parameters
Z_UP = F(31, 10)
LAM_UP = F(165439, 100000)
DM_UP = 100

# Lower-bound parameters
LMAX = 50
JMAX = 50
THETA_1 = F(865, 1000)
THETA_2 = F(867, 1000)
BETA_HI = F(10134, 10000)


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


def fmt_dec(q: F, digits: int = 6, round_up: bool = False) -> str:
    """Exact decimal outward formatting without floating-point arithmetic."""
    _check(q >= 0, "fmt_dec expects nonnegative Fraction")
    sc = 10 ** digits
    scaled = q * sc
    n = ceil_frac(scaled) if round_up else floor_frac(scaled)
    return f"{n // sc}.{n % sc:0{digits}d}"


def log_rational_lower(z: F, terms: int = 26) -> F:
    """Exact strict rational lower bound on log(z) for z > 1 via 2 * sum_{k=0}^{terms-1} w^{2k+1}/(2k+1)."""
    _check(isinstance(z, F) and z > 1, f"z must be Fraction > 1, got {z}")
    w = (z - 1) / (z + 1)
    w2 = w * w
    pw = w
    total = F(0)
    for k in range(terms):
        total += pw / (2 * k + 1)
        pw *= w2
    return 2 * total


def exp_neg_rational_upper(x: F, half_deg: int = 16) -> F:
    """
    Exact strict rational upper bound on exp(-x) for x > 0 via alternating Taylor sum
    sum_{k=0}^{2*half_deg} (-x)^k / k!, valid whenever 2*half_deg + 2 > x.
    """
    _check(isinstance(x, F) and x > 0, f"x must be positive Fraction, got {x}")
    _check(2 * half_deg + 2 > x, f"half_deg={half_deg} too small for x={x}")
    total = F(1)
    term = F(1)
    for k in range(1, 2 * half_deg + 1):
        term = -term * x / k
        total += term
    return total


def build_upper_table():
    """Pure-Fraction backward DP building integer-tick supermartingale table f[D][L] (scale SCALE)."""
    f21_cap = floor_frac((LAM_UP / Z_UP) * SCALE)
    f = [[SCALE] * (DM_UP + 2) for _ in range(DM_UP + 2)]
    for D in range(DM_UP, 1, -1):
        for L in range(0, D + 1):
            c = D - L
            b = F(L, D)
            fbd = f[D][L - 1] if L > 0 else 0
            cands = [
                (r * f[D + 1][L + 1] + b * (F(fbd) + Z_UP * f21_cap)) / (LAM_UP + r + b)
                for r in (F(L), F(2 * L + c, 2))
            ]
            f[D][L] = ceil_frac(max(cands))
    return f


def verify_upper_exact(f, z: F = Z_UP, lam: F = LAM_UP, dm: int = DM_UP):
    """Standalone zero-float verifier for the depth-tilted upper supermartingale table."""
    _check(isinstance(z, F) and z > 1, "z must be Fraction > 1")
    _check(isinstance(lam, F) and lam > 0, "lam must be Fraction > 0")
    _req_int(dm, "dm")
    _check(dm >= 2, "dm must be >= 2")
    _check(isinstance(f, (list, tuple)) and len(f) == dm + 2, f"f must have length {dm+2}")
    for D in range(2, dm + 2):
        _check(isinstance(f[D], (list, tuple)) and len(f[D]) >= D + 1, f"f[{D}] row too short")
    for L in range(0, dm + 2):
        _check(_req_int(f[dm + 1][L], f"f[{dm+1}][{L}]") == SCALE, f"boundary f[{dm+1}][{L}] != SCALE")
    f21 = _req_int(f[2][1], "f[2][1]")
    _check(f21 > 0, "f[2][1] must be positive")
    # Far-tail condition for ALL D >= dm + 1: z * f(2,1) <= lam * 1
    _check(z * f21 <= lam * SCALE, f"far-tail condition failed: z*f[2][1]={z*f21} > lam*S={lam*SCALE}")

    f_min = SCALE
    rows = 0
    for D in range(2, dm + 1):
        for L in range(0, D + 1):
            f0_i = _req_int(f[D][L], f"f[{D}][{L}]")
            _check(f0_i > 0, f"nonpositive f[{D}][{L}]={f0_i}")
            if f0_i < f_min:
                f_min = f0_i
            c = D - L
            b = F(L, D)
            f0 = F(f0_i)
            f_up = F(_req_int(f[D + 1][L + 1], f"f[{D+1}][{L+1}]"))
            f_bd = F(_req_int(f[D][L - 1], f"f[{D}][{L-1}]")) if L > 0 else F(0)
            for r in (F(L), F(2 * L + c, 2)):
                drift = r * (f_up - f0) + b * (f_bd - f0 + z * f21) - lam * f0
                _check(drift <= 0, f"upper supermartingale violation at D={D}, L={L}, r={r}: {drift}")
            rows += 1

    log_z_lb = log_rational_lower(z, terms=26)
    c_hi = lam / log_z_lb
    f_T3 = _req_int(f[2][2], "f[2][2]")
    prefactor = F(f_T3, f_min)
    return {
        "z": str(z),
        "lam": str(lam),
        "DM": dm,
        "rows_verified": rows,
        "f_2_1_ticks": f21,
        "f_2_2_ticks": f_T3,
        "f_min_ticks": f_min,
        "prefactor_Phi_T3_over_f_min_exact": str(prefactor),
        "prefactor_dec_up": fmt_dec(prefactor, 6, round_up=True),
        "log_z_lb_dec_dn": fmt_dec(log_z_lb, 12, round_up=False),
        "height_upper_exact_dec_up": fmt_dec(c_hi, 6, round_up=True),
        "diam_upper_exact_dec_up": fmt_dec(2 * c_hi, 6, round_up=True),
    }


def build_cmj_bounds(s: F, lmax: int = LMAX, jmax: int = JMAX, want_deriv: bool = True):
    """
    Pure-Fraction acyclic 1-pass generator for integer-tick sub/super-solutions
    u_lo <= phi(s) <= u_hi and (if want_deriv) v_lo <= -phi'(s) <= v_hi on (1..lmax, 1..jmax).
    """
    u_lo = [[0] * (jmax + 2) for _ in range(lmax + 1)]
    u_hi = [[0] * (jmax + 2) for _ in range(lmax + 1)]
    v_lo = [[0] * (jmax + 2) for _ in range(lmax + 1)] if want_deriv else None
    v_hi = [[0] * (jmax + 2) for _ in range(lmax + 1)] if want_deriv else None

    for j in range(jmax, 0, -1):
        # Reflected top boundary L = lmax (no gain jump!)
        rhs_lo = F(lmax * SCALE, lmax + j) + u_lo[lmax - 1][j + 1]
        rhs_hi = F(lmax * SCALE, lmax + j) + u_hi[lmax - 1][j + 1]
        u_lo[lmax][j] = floor_frac(rhs_lo / (s + 1))
        u_hi[lmax][j] = ceil_frac(rhs_hi / (s + 1))
        if want_deriv:
            v_lo[lmax][j] = floor_frac((F(u_lo[lmax][j]) + v_lo[lmax - 1][j + 1]) / (s + 1))
            v_hi[lmax][j] = ceil_frac((F(u_hi[lmax][j]) + v_hi[lmax - 1][j + 1]) / (s + 1))
        for L in range(lmax - 1, 0, -1):
            denom = s + L + 1
            u_lo[L][j] = floor_frac((L * u_lo[L + 1][j] + F(L * SCALE, L + j) + u_lo[L - 1][j + 1]) / denom)
            u_hi[L][j] = ceil_frac((L * u_hi[L + 1][j] + F(L * SCALE, L + j) + u_hi[L - 1][j + 1]) / denom)
            if want_deriv:
                v_lo[L][j] = floor_frac((u_lo[L][j] + L * v_lo[L + 1][j] + v_lo[L - 1][j + 1]) / denom)
                v_hi[L][j] = ceil_frac((u_hi[L][j] + L * v_hi[L + 1][j] + v_hi[L - 1][j + 1]) / denom)
    return {"u_lo": u_lo, "u_hi": u_hi, "v_lo": v_lo, "v_hi": v_hi}


def verify_cmj_table_exact(s: F, tbl: dict, lmax: int = LMAX, jmax: int = JMAX, check_deriv: bool = True):
    """
    Standalone zero-float verifier for sub/super-solution tables of phi(s) and -phi'(s).
    Verifies every boundary zero and every linear inequality on all lmax*jmax states.
    """
    _check(isinstance(s, F) and s >= 0, f"s must be nonnegative Fraction, got {s}")
    _req_int(lmax, "lmax")
    _req_int(jmax, "jmax")
    _check(lmax >= 2 and jmax >= 1, "require lmax >= 2 and jmax >= 1")
    u_lo, u_hi = tbl["u_lo"], tbl["u_hi"]
    v_lo, v_hi = tbl["v_lo"], tbl["v_hi"]
    for name, arr in (("u_lo", u_lo), ("u_hi", u_hi)) + ((("v_lo", v_lo), ("v_hi", v_hi)) if check_deriv else ()):
        _check(isinstance(arr, (list, tuple)) and len(arr) == lmax + 1, f"{name} must have {lmax+1} rows")
        for L in range(0, lmax + 1):
            _check(isinstance(arr[L], (list, tuple)) and len(arr[L]) >= jmax + 2, f"{name}[{L}] too short")

    # Boundary zero checks: L = 0 and j = jmax + 1
    for j in range(1, jmax + 2):
        _check(_req_int(u_lo[0][j], f"u_lo[0][{j}]") == 0, "u_lo[0][j] != 0")
        _check(_req_int(u_hi[0][j], f"u_hi[0][{j}]") == 0, "u_hi[0][j] != 0")
        if check_deriv:
            _check(_req_int(v_lo[0][j], f"v_lo[0][{j}]") == 0, "v_lo[0][j] != 0")
            _check(_req_int(v_hi[0][j], f"v_hi[0][{j}]") == 0, "v_hi[0][j] != 0")
    for L in range(0, lmax + 1):
        _check(_req_int(u_lo[L][jmax + 1], f"u_lo[{L}][end]") == 0, "u_lo[L][jmax+1] != 0")
        _check(_req_int(u_hi[L][jmax + 1], f"u_hi[{L}][end]") == 0, "u_hi[L][jmax+1] != 0")
        if check_deriv:
            _check(_req_int(v_lo[L][jmax + 1], f"v_lo[{L}][end]") == 0, "v_lo[L][jmax+1] != 0")
            _check(_req_int(v_hi[L][jmax + 1], f"v_hi[{L}][end]") == 0, "v_hi[L][jmax+1] != 0")

    rows = 0
    for j in range(1, jmax + 1):
        for L in range(1, lmax + 1):
            ul = _req_int(u_lo[L][j], f"u_lo[{L}][{j}]")
            uh = _req_int(u_hi[L][j], f"u_hi[{L}][{j}]")
            _check(0 <= ul <= uh, f"u ordering violation at ({L},{j})")
            ul_dn = _req_int(u_lo[L - 1][j + 1], "ul_dn")
            uh_dn = _req_int(u_hi[L - 1][j + 1], "uh_dn")
            emit = F(L * SCALE, L + j)
            if L == lmax:
                _check((s + 1) * ul <= emit + ul_dn, f"u_lo violation at top ({L},{j})")
                _check((s + 1) * uh >= emit + uh_dn, f"u_hi violation at top ({L},{j})")
            else:
                ul_up = _req_int(u_lo[L + 1][j], "ul_up")
                uh_up = _req_int(u_hi[L + 1][j], "uh_up")
                _check((s + L + 1) * ul <= L * ul_up + emit + ul_dn, f"u_lo violation at ({L},{j})")
                _check((s + L + 1) * uh >= L * uh_up + emit + uh_dn, f"u_hi violation at ({L},{j})")
            if check_deriv:
                vl = _req_int(v_lo[L][j], f"v_lo[{L}][{j}]")
                vh = _req_int(v_hi[L][j], f"v_hi[{L}][{j}]")
                _check(0 <= vl <= vh, f"v ordering violation at ({L},{j})")
                vl_dn = _req_int(v_lo[L - 1][j + 1], "vl_dn")
                vh_dn = _req_int(v_hi[L - 1][j + 1], "vh_dn")
                if L == lmax:
                    _check((s + 1) * vl <= ul + vl_dn, f"v_lo violation at top ({L},{j})")
                    _check((s + 1) * vh >= uh + vh_dn, f"v_hi violation at top ({L},{j})")
                else:
                    vl_up = _req_int(v_lo[L + 1][j], "vl_up")
                    vh_up = _req_int(v_hi[L + 1][j], "vh_up")
                    _check((s + L + 1) * vl <= ul + L * vl_up + vl_dn, f"v_lo violation at ({L},{j})")
                    _check((s + L + 1) * vh >= uh + L * vh_up + vh_dn, f"v_hi violation at ({L},{j})")
            rows += 1

    res = {"rows": rows, "u_lo_11": F(u_lo[1][1], SCALE), "u_hi_11": F(u_hi[1][1], SCALE)}
    if check_deriv:
        res["v_lo_11"] = F(v_lo[1][1], SCALE)
        res["v_hi_11"] = F(v_hi[1][1], SCALE)
    return res


def verify_lower_exact(tbl0, tbl1, tbl2, th1: F = THETA_1, th2: F = THETA_2, beta_hi: F = BETA_HI):
    """
    Standalone zero-float verifier proving:
      1. Supercriticality: phi(0) = E[N] >= u_lo_0 > 1.
      2. Convex global positivity of h(theta) = log phi(theta) + beta_hi * theta on [0, infinity):
         - h'(th1) <= 0  (checked via v_lo(th1)/u_hi(th1) >= beta_hi),
         - h'(th2) >= 0  (checked via v_hi(th2)/u_lo(th2) <= beta_hi),
         - h(th2) > 0    (checked via u_lo(th2) > exp_neg_rational_upper(beta_hi * th2)),
         - min_{theta in [th1, th2]} h(theta) >= h(th1) + h'(th1)*(th2 - th1) > 0
           (checked via u_lo(th1) > exp_neg_rational_upper(beta_hi * th1 - slope_hi * (th2 - th1))
            where slope_hi = v_hi(th1)/u_lo(th1) - beta_hi >= -h'(th1) >= 0).
    """
    _check(isinstance(th1, F) and isinstance(th2, F) and isinstance(beta_hi, F), "params must be Fractions")
    _check(0 < th1 < th2 and beta_hi > 0, "require 0 < th1 < th2 and beta_hi > 0")

    r0 = verify_cmj_table_exact(F(0), tbl0, check_deriv=False)
    _check(r0["u_lo_11"] > 1, f"CMJ tree not supercritical: E[N] >= {r0['u_lo_11']} <= 1")

    r1 = verify_cmj_table_exact(th1, tbl1, check_deriv=True)
    r2 = verify_cmj_table_exact(th2, tbl2, check_deriv=True)

    # 1. h'(th1) <= 0  <=>  -phi'(th1)/phi(th1) >= beta_hi  (implied by v_lo(th1)/u_hi(th1) >= beta_hi)
    _check(r1["v_lo_11"] >= beta_hi * r1["u_hi_11"], "h'(th1) <= 0 check failed")
    # 2. h'(th2) >= 0  <=>  -phi'(th2)/phi(th2) <= beta_hi  (implied by v_hi(th2)/u_lo(th2) <= beta_hi)
    _check(r2["v_hi_11"] <= beta_hi * r2["u_lo_11"], "h'(th2) >= 0 check failed")

    # 3. h(th2) > 0  <=>  phi(th2) > exp(-beta_hi * th2)
    e2_hi = exp_neg_rational_upper(beta_hi * th2, half_deg=16)
    _check(r2["u_lo_11"] > e2_hi, f"h(th2) > 0 check failed: {r2['u_lo_11']} <= {e2_hi}")

    # 4. Uniform positivity on [th1, th2] via tangent line at th1:
    #    For theta in [th1, th2], h(theta) >= h(th1) + h'(th1)*(theta - th1) >= h(th1) - slope_hi*(th2 - th1)
    #    where slope_hi = v_hi(th1)/u_lo(th1) - beta_hi >= -h'(th1) >= 0.
    #    So h(theta) > 0 on [th1, th2] (and on [0, th1] since h'(th1) <= 0!) whenever
    #    phi(th1) > exp(-(beta_hi * th1 - slope_hi * (th2 - th1))).
    slope_hi = r1["v_hi_11"] / r1["u_lo_11"] - beta_hi
    _check(slope_hi >= 0, "unexpected negative slope_hi")
    x1_eff = beta_hi * th1 - slope_hi * (th2 - th1)
    _check(x1_eff > 0, f"effective exponent nonpositive: {x1_eff}")
    e1_hi = exp_neg_rational_upper(x1_eff, half_deg=16)
    _check(r1["u_lo_11"] > e1_hi, f"h([th1,th2]) > 0 check failed: {r1['u_lo_11']} <= {e1_hi}")

    c_lo = F(1) / beta_hi
    return {
        "LMAX": LMAX,
        "JMAX": JMAX,
        "theta_1": str(th1),
        "theta_2": str(th2),
        "beta_hi": str(beta_hi),
        "rows_verified_total": r0["rows"] + r1["rows"] + r2["rows"],
        "phi_0_lb_dec_dn": fmt_dec(r0["u_lo_11"], 6, round_up=False),
        "phi_th1_dec": [fmt_dec(r1["u_lo_11"], 8, False), fmt_dec(r1["u_hi_11"], 8, True)],
        "phi_th2_dec": [fmt_dec(r2["u_lo_11"], 8, False), fmt_dec(r2["u_hi_11"], 8, True)],
        "height_lower_exact": str(c_lo),
        "height_lower_dec_dn": fmt_dec(c_lo, 6, round_up=False),
        "diam_lower_exact": str(2 * c_lo),
        "diam_lower_dec_dn": fmt_dec(2 * c_lo, 6, round_up=False),
    }


def verify_physical_trees_exact(f, z: F = Z_UP, lam: F = LAM_UP, dm: int = DM_UP):
    """
    Independently check the exact physical continuous-time generator inequality
        (L Phi_z)(T) <= lam * Phi_z(T)
    on ALL 23 unlabelled trees of sizes n = 3..7 (every rooted choice of core root o in K(T))
    by directly enumerating every physical directed arrow (caller -> neighbor) at rate 1/D_caller.
    """
    def canon(adj):
        n = len(adj)
        def enc(u, p):
            return "(" + "".join(sorted(enc(v, u) for v in adj[u] if v != p)) + ")"
        return min(enc(r, -1) for r in range(n))

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

    def eval_phi(adj, root):
        n = len(adj)
        deg = [len(adj[u]) for u in range(n)]
        # BFS distances from root
        dist = [-1] * n
        dist[root] = 0
        q = [root]
        for u in q:
            for v in adj[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    q.append(v)
        val = F(0)
        for u in range(n):
            if deg[u] >= 2:
                L = sum(1 for v in adj[u] if deg[v] == 1)
                fw = F(f[deg[u]][L]) if deg[u] <= dm else F(SCALE)
                val += (z ** dist[u]) * fw
        return val

    trees_checked = 0
    rooted_checked = 0
    for n in range(3, 8):
        for adj in trees_by_n[n]:
            deg = [len(adj[u]) for u in range(n)]
            cores = [u for u in range(n) if deg[u] >= 2]
            for root in cores:
                phi0 = eval_phi(adj, root)
                drift = F(0)
                for caller in range(n):
                    for target in adj[caller]:
                        na = [list(x) for x in adj] + [[target]]
                        na[target].append(n)
                        phi1 = eval_phi(na, root)
                        drift += F(1, deg[caller]) * (phi1 - phi0)
                _check(drift <= lam * phi0, f"physical tree drift violation at n={n}, root={root}")
                rooted_checked += 1
            trees_checked += 1
    return {"unlabelled_trees_n3_to_n7": trees_checked, "rooted_core_configs_checked": rooted_checked}


def run_negative_controls(f_good, tbl0_good, tbl1_good, tbl2_good):
    """Verify 11 planted defects across upper, lower, convexity, and type checks fail closed."""
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

    # Upper controls
    def mut_f(fn):
        fm = copy.deepcopy(f_good)
        fn(fm)
        verify_upper_exact(fm)

    expect_fail("upper_interior_tick", lambda: mut_f(lambda fm: fm[2].__setitem__(0, fm[2][0] - 1)))
    expect_fail("upper_far_tail_tick", lambda: mut_f(lambda fm: fm[2].__setitem__(1, floor_frac((LAM_UP / Z_UP) * SCALE) + 1)))
    expect_fail("upper_boundary_tick", lambda: mut_f(lambda fm: fm[DM_UP + 1].__setitem__(5, SCALE - 1)))
    expect_fail("upper_too_tight_lam", lambda: verify_upper_exact(f_good, lam=F(165, 100)))
    expect_fail("upper_float_tick", lambda: mut_f(lambda fm: fm[2].__setitem__(0, float(fm[2][0]))))

    # Lower controls
    def mut_t1(fn):
        t1m = copy.deepcopy(tbl1_good)
        fn(t1m)
        verify_lower_exact(tbl0_good, t1m, tbl2_good)

    expect_fail("lower_u_lo_inflated", lambda: mut_t1(lambda t: t["u_lo"][1].__setitem__(1, t["u_lo"][1][1] + 1)))
    expect_fail("lower_v_hi_deflated", lambda: mut_t1(lambda t: t["v_hi"][1].__setitem__(1, t["v_hi"][1][1] - 1)))
    expect_fail("lower_boundary_nonzero", lambda: mut_t1(lambda t: t["u_lo"][0].__setitem__(1, 1)))
    expect_fail("lower_bool_tick", lambda: mut_t1(lambda t: t["u_lo"][0].__setitem__(1, False)))
    expect_fail("lower_too_small_beta_hi", lambda: verify_lower_exact(tbl0_good, tbl1_good, tbl2_good, beta_hi=F(1, 1)))
    expect_fail("lower_non_supercritical", lambda: verify_lower_exact(
        {"u_lo": [[0] * (JMAX + 2) for _ in range(LMAX + 1)], "u_hi": tbl0_good["u_hi"], "v_lo": None, "v_hi": None},
        tbl1_good, tbl2_good))
    return out


def main():
    f_up = build_upper_table()
    up_res = verify_upper_exact(f_up)
    tree_res = verify_physical_trees_exact(f_up)

    tbl0 = build_cmj_bounds(F(0), want_deriv=False)
    tbl1 = build_cmj_bounds(THETA_1, want_deriv=True)
    tbl2 = build_cmj_bounds(THETA_2, want_deriv=True)
    lo_res = verify_lower_exact(tbl0, tbl1, tbl2)

    ctrls = run_negative_controls(f_up, tbl0, tbl1, tbl2)
    out = {
        "status": "HEIGHT_DIAMETER_EXACT_CHECKS_PASSED",
        "scale": SCALE,
        "upper_bound": up_res,
        "physical_tree_checks": tree_res,
        "lower_bound": lo_res,
        "negative_controls": ctrls,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
