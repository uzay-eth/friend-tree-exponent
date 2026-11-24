#!/usr/bin/env python3
"""
Standalone exact rational tensor-Bernstein verifier for FT2_M5V1 signed coboundary
analytic singleton tails (D >= DF), evaluated at BOTH rational endpoints [mu_lo, mu_hi].
Zero floating-point arithmetic.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import json
import struct
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "quintet-cert"))
import verify_tail_bernstein as vtb

M5_MAGIC = 0x4654325F4D355631
Q5_MAX_W = 16 * (1 << 40)
Q5_MAX_Q = 1 << 40


def read_marginal_header_bytes(raw):
    if len(raw) < 640:
        raise ValueError("short FT2_M5V1 header")
    mvals = struct.unpack("<Q15q", raw[:128])
    mnames = (
        "m_magic obs_id obs_sign theta_P theta_Q f_p3 "
        "mu_lo_P mu_lo_Q mu_hi_P mu_hi_Q "
        "m_res0 m_res1 m_res2 m_res3 m_res4 m_res5"
    ).split()
    mh = dict(zip(mnames, mvals))
    if mh["m_magic"] != M5_MAGIC:
        raise ValueError("invalid FT2_M5V1 magic")
    if not (1 <= mh["obs_id"] <= 8) or mh["obs_sign"] not in (-1, 1):
        raise ValueError("invalid obs_id or obs_sign")
    if not (0 < mh["theta_Q"] <= Q5_MAX_Q and 0 <= mh["theta_P"] <= mh["theta_Q"]):
        raise ValueError("invalid theta_P/theta_Q")
    if not (-Q5_MAX_W <= mh["f_p3"] <= Q5_MAX_W):
        raise ValueError("f_p3 out of range")
    if not (0 < mh["mu_lo_P"] < mh["mu_lo_Q"] <= Q5_MAX_Q):
        raise ValueError("invalid mu_lo")
    if not (0 < mh["mu_hi_P"] < mh["mu_hi_Q"] <= Q5_MAX_Q):
        raise ValueError("invalid mu_hi")
    if F(mh["mu_lo_P"], mh["mu_lo_Q"]) > F(mh["mu_hi_P"], mh["mu_hi_Q"]):
        raise ValueError("mu_lo > mu_hi")
    for i in range(6):
        if mh[f"m_res{i}"] != 0:
            raise ValueError("nonzero reserved field in MarginalHeader")
    base = vtb.read_header_bytes(raw[128:640])
    if base["mode"] != 1:
        raise ValueError("base.mode must be 1 for FT2_M5V1")
    mh["base"] = base
    return mh


def read_marginal_tail_constants_bytes(raw):
    mh = read_marginal_header_bytes(raw)
    h = mh["base"]
    d2, K2, DF, nknots = h["d2"], h["K2"], h["DF"], h["nknots"]
    nin_ex = sum(min(d, K2) for d in range(2, d2))
    n_ovf = max(0, d2 - 1 - K2)
    tail_base = nin_ex + n_ovf
    g2_inf = tail_base + (K2 + 1)
    nf = (DF - 2) * (DF + 3) // 2
    off_f = nin_ex + nin_ex * nknots
    off_g = off_f + nf

    def item(ix):
        pos = 640 + 8 * ix
        if pos + 8 > len(raw):
            raise ValueError("truncated FT2_M5V1 snapshot")
        return struct.unpack_from("<q", raw, pos)[0]

    def sym_off(a, b):
        if a > b:
            a, b = b, a
        return b * (b + 1) // 2 + a

    j_g2 = 0
    f_J = item(off_f + 1)
    g_inf_J = item(off_g + sym_off(g2_inf, j_g2))
    tail_vals = [item(off_g + sym_off(tail_base + c - 1, j_g2)) for c in range(1, K2 + 2)]
    # Always mode == 1 (UPPER bound on L Phi - lam Phi + sgn*(psi - theta))
    ext_val = max(tail_vals)
    d_g = d2 * (ext_val - g_inf_J)
    return mh, f_J + g_inf_J, d_g


def verify_bytes(raw, mu_lo=None, mu_hi=None, theta_override=None, panels=16):
    if not isinstance(panels, int) or panels < 1:
        raise ValueError("panels must be positive int")
    mh, C_nb, d_g = read_marginal_tail_constants_bytes(raw)
    h = mh["base"]
    S = h["S"]
    lam_lo = mu_lo if mu_lo is not None else F(mh["mu_lo_P"], mh["mu_lo_Q"])
    lam_hi = mu_hi if mu_hi is not None else F(mh["mu_hi_P"], mh["mu_hi_Q"])
    theta = theta_override if theta_override is not None else F(mh["theta_P"], mh["theta_Q"])
    obs_sign = mh["obs_sign"]
    # Unsigned ticks polynomial: w(q, x) = a + b*q + c*q^2 + e*q*x
    a, b, c, e = h["f_p0"], h["f_p1"], h["f_p2"], mh["f_p3"]
    Q_VAR, X_VAR = vtb.Q_VAR, vtb.X_VAR
    w = a + b * Q_VAR + c * Q_VAR * Q_VAR + e * Q_VAR * X_VAR
    den = (1 + X_VAR) * (1 + X_VAR)
    gain = (
        b * (1 - Q_VAR) * (1 + X_VAR)
        + c * (1 - Q_VAR) * (2 * Q_VAR + (1 + Q_VAR) * X_VAR)
        + e * ((1 - 2 * Q_VAR) * X_VAR - Q_VAR * X_VAR * X_VAR)
    )
    bud = -b * X_VAR + c * (-2 * Q_VAR * X_VAR + X_VAR * X_VAR) - e * X_VAR * X_VAR
    src_const = -obs_sign * theta * S  # psi(D, c) == 0 on D >= DF
    xmax = F(1, h["DF"])

    caps = {}
    for label, lam in (("lo", lam_lo), ("hi", lam_hi)):
        cap = None
        for extreme, rho in enumerate((Q_VAR, (1 + Q_VAR) * F(1, 2))):
            p = rho * gain + (Q_VAR * (bud + C_nb + d_g * X_VAR) - w * lam + src_const) * den
            # Independent exact oracle check at concrete integer states (D, L)
            for D in (h["DF"], h["DF"] + 1, 3 * h["DF"] + 17):
                x = F(1, D)
                for L in (0, 1, D // 3, D - 1, D):
                    q = L * x
                    r = F(L) if extreme == 0 else F(D + L, 2)
                    direct = r * (w.at(F(L + 1, D + 1), F(1, D + 1)) - w.at(q, x))
                    direct += q * (w.at(q - x, x) - w.at(q, x) + C_nb + d_g * x) - lam * w.at(q, x) + src_const
                    if p.at(q, x) != direct * (1 + x) ** 2:
                        raise ArithmeticError(f"tail oracle mismatch at D={D}, L={L}")
            for j in range(panels):
                top = max(vtb.bernstein_2d(p, F(j, panels), F(j + 1, panels), xmax))
                bound = top if top >= 0 else top / ((1 + xmax) ** 2)
                cap = bound if cap is None else max(cap, bound)
        caps[label] = -((-cap.numerator) // cap.denominator)

    w_min = min(vtb.bernstein_2d(w, F(0), F(1), xmax))
    w_max = max(vtb.bernstein_2d(w, F(0), F(1), xmax))
    return {
        "family": "f_analytic_marginal",
        "obs_id": mh["obs_id"],
        "obs_sign": obs_sign,
        "DF": h["DF"],
        "panels": panels,
        "cap_ticks_lo": caps["lo"],
        "cap_ticks_hi": caps["hi"],
        "w_min_ticks": w_min.numerator // w_min.denominator,
        "w_max_ticks": -((-w_max.numerator) // w_max.denominator),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Tail-only (D >= DF) exact rational 2D tensor-Bernstein verifier for FT2_M5V1; does not validate finite tables or full header bounds."
    )
    ap.add_argument("bin_path")
    ap.add_argument("--panels", type=int, default=16)
    args = ap.parse_args()
    print(json.dumps(verify_bytes(Path(args.bin_path).read_bytes(), panels=args.panels), indent=2))
