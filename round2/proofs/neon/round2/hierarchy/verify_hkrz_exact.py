#!/usr/bin/env python3
"""Verify exact HKRZ Fig. 6 extracted points and reproduce their 0.569, 0.773, 0.545 fits!"""
import math
import numpy as np

# Exact coordinates from Stream 0 of arXiv:2604.01540 CvsN.pdf:
y10 = 52.707678
dy_dec = 32.12199166666667  # (149.073653 - 52.707678) / 3

def decode_series(y0, dys):
    ys = [y0]
    for d in dys:
        ys.append(ys[-1] + d)
    vals = [10.0**(1.0 + (y - y10) / dy_dec) for y in ys]
    slopes = [d / dy_dec for d in dys]
    return np.array(vals), np.array(slopes)

N5 = np.array([1e2, 1e3, 1e4, 1e5, 1e6])

pan_vals, pan_sl = decode_series(62.8074835276, [18.5890982579, 18.1390397472, 17.779898585, 17.4722471226])
dan_vals, dan_sl = decode_series(64.3735205558, [24.9380511219, 24.8067352083, 24.6625314941, 24.8462998472])
ir_vals,  ir_sl  = decode_series(52.2076523769, [18.3894157741, 18.2331034833, 18.006979511, 18.3370259431])

def nls_power(N, y):
    logN = np.log(N)
    lo, hi = 0.3, 0.95
    gr = (math.sqrt(5) - 1) / 2
    def loss(b):
        x = np.exp(b * logN)
        a = np.dot(x, y) / np.dot(x, x)
        return np.sum((y - a * x)**2)
    c = hi - gr * (hi - lo); d = lo + gr * (hi - lo)
    fc, fd = loss(c), loss(d)
    for _ in range(90):
        if fc < fd: hi, d, fd = d, c, fc; c = hi - gr * (hi - lo); fc = loss(c)
        else:       lo, c, fc = c, d, fd; d = lo + gr * (hi - lo); fd = loss(d)
    b = 0.5 * (lo + hi)
    x = np.exp(b * logN)
    a = np.dot(x, y) / np.dot(x, x)
    return a, b

def log_ls(N, y):
    x = np.log(N); ly = np.log(y); xm = x.mean()
    return np.dot(x - xm, ly) / np.dot(x - xm, x - xm)

for name, vals, sl, true_mu, hkrz_mu in [
    ("DAN", dan_vals, dan_sl, 0.7711923, 0.773),
    ("PAN", pan_vals, pan_sl, 0.5497350, 0.545),
    ("IR ", ir_vals,  ir_sl,  0.5657230, 0.569),
]:
    a, b_nls = nls_power(N5, vals)
    b_log = log_ls(N5, vals)
    print(f"=== {name} (true={true_mu:.6f}, HKRZ caption={hkrz_mu}) ===")
    print("  Extracted <C(N)>:", [round(float(v), 3) for v in vals])
    print("  Decade slopes:   ", [round(float(s), 6) for s in sl])
    print(f"  Unweighted NLS fit:  b = {b_nls:.6f}  (a = {a:.5f})  <-- HKRZ reports {hkrz_mu}!")
    print(f"  Log-log OLS fit:     b = {b_log:.6f}")

# Also check the exact slopes of the 3 straight lines drawn in CvsN.pdf!
print("\n=== Slopes of the 3 reference lines drawn in CvsN.pdf ===")
print("  PAN red line  :", (81.761297 - 64.102714) / dy_dec)
print("  DAN blue line :", (88.892760 - 64.120537) / dy_dec)
print("  IR  green line:", (70.350950 - 52.169903) / dy_dec)
