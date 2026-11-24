#!/usr/bin/env python3
"""
Standalone exact rational tensor-Bernstein verifier for the FT2_U5V1 analytic vertex tail (D >= DF).
Zero floating-point arithmetic: reads the binary snapshot directly via struct, checks the exact
polynomial identity against concrete integer states (D, L), and bounds the residual over ALL
q in [0, 1], x = 1/D in [0, 1/DF] via rational Bernstein panels.

NOTE: Standalone execution verifies ONLY the analytic singleton tail (D >= DF) and its
comparability extrema, NOT the full FT2_U5V1 certificate. Use cert_check.py for full verification.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import json
import struct

Q5_MAGIC = 0x4654325F55355631
Q5_MAX_W = 16 * (1 << 40)


def read_header_bytes(raw):
    if len(raw) < 512:
        raise ValueError("short header")
    vals = struct.unpack("<Q63q", raw[:512])
    names = (
        "magic mode P Q S S_r alpha_lo_num f_p0 f_p1 f_p2 "
        "DF d2 K2 d3m K3m d3e K3e d4i K4i SP4 dS3 KS3 d4e K4e "
        "d5i K5i SP5 d5o K5o d5c K5c SB5 dS4 KS4 d5e K5e "
        "Q_K nknots"
    ).split()
    h = dict(zip(names, vals[:38]))
    h["knot_num"] = list(vals[38:62])
    h["reserved0"] = vals[62]
    h["reserved1"] = vals[63]
    if h["magic"] != Q5_MAGIC or h["mode"] not in (1, 2):
        raise ValueError("invalid FT2_U5V1 header")
    if not (0 < h["P"] < h["Q"] <= (1 << 40)):
        raise ValueError("invalid P/Q in header")
    if not (20 <= h["DF"] <= 4096 and 6 <= h["d2"] < h["DF"] and 2 <= h["K2"] < h["d2"]):
        raise ValueError("invalid DF/d2/K2 in header")
    if not (3 <= h["nknots"] <= 24):
        raise ValueError("invalid nknots in header")
    for k in ("f_p0", "f_p1", "f_p2"):
        if not (-Q5_MAX_W <= h[k] <= Q5_MAX_W):
            raise ValueError(f"analytic parameter {k} out of range")
    return h


def read_header(fp):
    return read_header_bytes(fp.read(512))


def read_tail_constants_bytes(raw):
    h = read_header_bytes(raw)
    d2, K2, DF, nknots = h["d2"], h["K2"], h["DF"], h["nknots"]
    # Reconstruct exact G2 layout independently of C++
    nin_ex = sum(min(d, K2) for d in range(2, d2))
    n_ovf = max(0, d2 - 1 - K2)
    tail_base = nin_ex + n_ovf
    g2_inf = tail_base + (K2 + 1)
    # f(D, c) index: (D-2)*(D+3)//2 + c; J = (2, 1) -> idx 1
    nf = (DF - 2) * (DF + 3) // 2
    off_f = nin_ex + nin_ex * nknots
    off_g = off_f + nf

    def item(ix):
        pos = 512 + 8 * ix
        if pos + 8 > len(raw):
            raise ValueError("truncated snapshot")
        return struct.unpack_from("<q", raw, pos)[0]

    def sym_off(a, b):
        if a > b:
            a, b = b, a
        return b * (b + 1) // 2 + a

    # J in G2 is (D=2, c=1) -> exact index 0!
    j_g2 = 0
    f_J = item(off_f + 1)  # f(2, 1)
    g_inf_J = item(off_g + sym_off(g2_inf, j_g2))
    # After promotion of a tail vertex (D >= DF, c -> c+1 >= 1), c+1 clamps to 1..K2+1
    tail_vals = [item(off_g + sym_off(tail_base + c - 1, j_g2)) for c in range(1, K2 + 2)]
    ext_val = max(tail_vals) if h["mode"] == 1 else min(tail_vals)
    d_g = d2 * (ext_val - g_inf_J)
    return h, f_J + g_inf_J, d_g


def read_tail_constants(path):
    return read_tail_constants_bytes(Path(path).read_bytes())


class Poly(dict):
    def __add__(self, other):
        other = poly(other)
        return Poly({k: self.get(k, 0) + other.get(k, 0) for k in self.keys() | other.keys()})
    __radd__ = __add__
    def __neg__(self):
        return Poly({k: -v for k, v in self.items()})
    def __sub__(self, other):
        return self + (-poly(other))
    def __rsub__(self, other):
        return poly(other) + (-self)
    def __mul__(self, other):
        ans = Poly()
        for (i, j), v in self.items():
            for (k, l), w in poly(other).items():
                key = (i + k, j + l)
                ans[key] = ans.get(key, 0) + v * w
        return ans
    __rmul__ = __mul__
    def __pow__(self, n):
        res = poly(1)
        for _ in range(n):
            res = res * self
        return res
    def at(self, q, x):
        ans = 0
        for (i, j), v in self.items():
            ans = ans + v * (q ** i) * (x ** j)
        return ans


def poly(v):
    return v if isinstance(v, Poly) else Poly({(0, 0): v})


Q_VAR = Poly({(1, 0): 1})
X_VAR = Poly({(0, 1): 1})


def bernstein_2d(p_, q0, q1, xmax):
    p = {k: v for k, v in p_.items() if v} or {(0, 0): 0}
    nq = max(i for i, _ in p)
    nx = max(j for _, j in p)
    power = {}
    for (i, j), v in p.items():
        for a in range(i + 1):
            key = (a, j)
            power[key] = power.get(key, 0) + v * comb(i, a) * (q0 ** (i - a)) * ((q1 - q0) ** a) * (xmax ** j)
    return [
        sum(
            c * F(comb(i, a), comb(nq, a)) * F(comb(j, b), comb(nx, b))
            for (a, b), c in power.items()
            if a <= i and b <= j
        )
        for i in range(nq + 1)
        for j in range(nx + 1)
    ]


def verify_bytes(raw, target=None, panels=16):
    if not isinstance(panels, int) or panels < 1:
        raise ValueError(f"panels must be a positive integer, got {panels!r}")
    h, C_nb, d_g = read_tail_constants_bytes(raw)
    lam = target if target is not None else F(h["P"], h["Q"])
    upper = (h["mode"] == 1)
    sgn = 1 if upper else -1
    a, b, c, e = ((h["f_p0"], h["f_p1"], 0, 0) if upper else (0, h["f_p0"], -h["f_p1"], -h["f_p2"]))
    # Unsigned analytic weight in integer ticks: w(q, x) = a + b*q + c*q^2 + e*q*x where q = L/D, x = 1/D.
    w = a + b * Q_VAR + c * Q_VAR * Q_VAR + e * Q_VAR * X_VAR
    den = (1 + X_VAR) * (1 + X_VAR)
    # Exact identities:
    #   gain = (1+x)^2 * (w((q+x)/(1+x), x/(1+x)) - w(q, x)) / x
    #   bud  = w(q-x, x) - w(q, x)
    gain = b * (1 - Q_VAR) * (1 + X_VAR) + c * (1 - Q_VAR) * (2 * Q_VAR + (1 + Q_VAR) * X_VAR) + e * ((1 - 2 * Q_VAR) * X_VAR - Q_VAR * X_VAR * X_VAR)
    bud = -b * X_VAR + c * (-2 * Q_VAR * X_VAR + X_VAR * X_VAR) - e * X_VAR * X_VAR

    cap = None
    xmax = F(1, h["DF"])
    # For unoriented (D, c), L = D - c = q*D, c = (1-q)*D, gain rate r in [L, L + c/2] = [q/x, (1+q)/(2x)]!
    for extreme, rho in enumerate((Q_VAR, (1 + Q_VAR) * F(1, 2))):
        p = sgn * (rho * gain + (Q_VAR * (bud + C_nb + d_g * X_VAR) - w * lam) * den)
        # Independent exact oracle check at concrete integer states (D, L)
        for D in (h["DF"], h["DF"] + 1, 3 * h["DF"] + 17):
            x = F(1, D)
            for L in (0, 1, D // 3, D - 1, D):
                q = L * x
                r = F(L) if extreme == 0 else F(D + L, 2)
                direct = r * (w.at(F(L + 1, D + 1), F(1, D + 1)) - w.at(q, x))
                direct += q * (w.at(q - x, x) - w.at(q, x) + C_nb + d_g * x) - lam * w.at(q, x)
                if p.at(q, x) != sgn * direct * (1 + x) ** 2:
                    raise ArithmeticError(f"gain/bud polynomial oracle mismatch at D={D}, L={L}")
        for j in range(panels):
            top = max(bernstein_2d(p, F(j, panels), F(j + 1, panels), xmax))
            bound = top if top >= 0 else top / ((1 + xmax) ** 2)
            cap = bound if cap is None else max(cap, bound)

    # Also verify exact analytic comparison bounds on D >= DF (q in [0, 1], x in [0, 1/DF]):
    # Upper: w(q, x) = a + b*q (a >= 0, b >= 0) -> min = a, max = a + b
    # Lower: w(q, x) = q*(f_p0 - f_p1*q - f_p2*x): check monotonicity / nonnegativity!
    w_min = min(bernstein_2d(w, F(0), F(1), xmax))
    w_max = max(bernstein_2d(w, F(0), F(1), xmax))
    # For c = 0 (single-core star, q = 1): min_{x in [0, xmax]} w(1, x)
    w_q1_min = min(w.at(F(1), F(0)), w.at(F(1), xmax))
    # For c = 1 (core leaf, q = 1 - x): min_{x in [0, xmax]} w(1-x, x)
    w_c1_poly = w.at(1 - X_VAR, X_VAR)
    w_c1_min = min(bernstein_2d(w_c1_poly, F(0), F(1), xmax))
    # For c = 2 (q = 1 - 2*x): min_{x in [0, xmax]} w(1-2x, x)
    w_c2_poly = w.at(1 - 2 * X_VAR, X_VAR)
    w_c2_min = min(bernstein_2d(w_c2_poly, F(0), F(1), xmax))

    cap_int = -((-cap.numerator) // cap.denominator)
    return {
        "family": "f_analytic",
        "mode": h["mode"],
        "DF": h["DF"],
        "panels": panels,
        "P": lam.numerator,
        "Q": lam.denominator,
        "cap_ticks": cap_int,
        "cap_rational": str(cap),
        "w_min_ticks": int(w_min.numerator // w_min.denominator),
        "w_max_ticks": -((-w_max.numerator) // w_max.denominator),
        "w_c0_min_ticks": int(w_q1_min.numerator // w_q1_min.denominator),
        "w_c1_min_ticks": int(w_c1_min.numerator // w_c1_min.denominator),
        "w_c2_min_ticks": int(w_c2_min.numerator // w_c2_min.denominator),
    }


def verify(path, target=None, panels=16):
    return verify_bytes(Path(path).read_bytes(), target=target, panels=panels)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("snapshot")
    ap.add_argument("--target", type=F, default=None)
    ap.add_argument("--panels", type=int, default=16)
    args = ap.parse_args()
    if args.panels < 1:
        ap.error("--panels must be >= 1")
    print(json.dumps(verify(args.snapshot, args.target, args.panels), sort_keys=True))
