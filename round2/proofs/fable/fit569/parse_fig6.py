#!/usr/bin/env python3
"""Extract the exact data points and fitted curves from Fig. 6 (CvsN.pdf) of arXiv:2604.01540
and check which fitting method reproduces the reported exponents (0.569, 0.773, 0.545).
"""
import numpy as np
from scipy.optimize import curve_fit

# From cvsn_stream.txt (exact PDF vector coordinates from CvsN.pdf):
# X-axis ticks: 10^2 .. 10^6
x_ticks = {2: 63.827267, 3: 114.630298, 4: 165.433330, 5: 216.236362, 6: 267.039394}
x_per_dec = (x_ticks[6] - x_ticks[2]) / 4.0

# Y-axis ticks: 10^1 .. 10^4
y_ticks = {1: 52.707678, 2: 84.829670, 3: 116.951662, 4: 149.073653}
y_per_dec = (y_ticks[4] - y_ticks[1]) / 3.0
y0 = y_ticks[1] - 1.0 * y_per_dec  # y for <C> = 1

def y_to_C(y):
    return 10.0 ** ((y - y0) / y_per_dec)

# Cumulative marker translations (cm) for each series in Fig. 6:
series_deltas = {
    "PAN (red /M0)": [62.8074835276, 18.5890982579, 18.1390397472, 17.779898585, 17.4722471226],
    "DAN (blue /M1)": [64.3735205558, 24.9380511219, 24.8067352083, 24.6625314941, 24.8462998472],
    "IR (green /M2)": [52.2076523769, 18.3894157741, 18.2331034833, 18.006979511, 18.3370259431],
}

# Plotted fit lines (y at N=10^2 and N=10^6):
fit_lines = {
    "PAN (red /M0)": (64.102714, 134.737046),
    "DAN (blue /M1)": (64.120537, 163.209429),
    "IR (green /M2)": (52.169903, 124.894092),
}

Ns = np.array([1e2, 1e3, 1e4, 1e5, 1e6])

print(f"x_per_dec = {x_per_dec:.10f}, y_per_dec = {y_per_dec:.10f}")
print()
for name, deltas in series_deltas.items():
    ys = np.cumsum(deltas)
    Cs = y_to_C(ys)
    local_slopes = np.diff(ys) / y_per_dec

    y_line_start, y_line_end = fit_lines[name]
    line_slope = (y_line_end - y_line_start) / (4.0 * y_per_dec)
    line_a = y_to_C(y_line_start) / (100.0 ** line_slope)

    # Fit (i): unweighted NLS in linear space: C = a * N^b
    popt_nls, _ = curve_fit(lambda x, a, b: a * (x ** b), Ns, Cs, p0=[1.0, 0.6], maxfev=20000)
    # Fit (ii): unweighted OLS in log-log space
    b_log, log_a = np.polyfit(np.log(Ns), np.log(Cs), 1)
    # Fit (iii): NLS with additive constant: C = a * N^b + c
    popt_nlsc, _ = curve_fit(lambda x, a, b, c: a * (x ** b) + c, Ns, Cs, p0=[popt_nls[0], popt_nls[1], 0.0], maxfev=20000)
    # Fit (iv): NLS with relative weighting (sigma = C) -> approx log-log
    popt_rel, _ = curve_fit(lambda x, a, b: a * (x ** b), Ns, Cs, sigma=Cs, p0=[1.0, 0.6], maxfev=20000)

    print(f"=== {name} ===")
    print("  Extracted <C> at N=[1e2, 1e3, 1e4, 1e5, 1e6]:")
    print("   ", [round(float(c), 4) for c in Cs])
    print("  Adjacent decade local slopes (1e2-1e3, 1e3-1e4, 1e4-1e5, 1e5-1e6):")
    print("   ", [round(float(s), 6) for s in local_slopes])
    print(f"  Plotted straight line in Fig. 6:  slope = {line_slope:.6f}, a = {line_a:.6f}")
    print(f"  Fit (i) linear NLS a*N^b:         b = {popt_nls[1]:.6f}, a = {popt_nls[0]:.6f}")
    print(f"  Fit (ii) log-log OLS:             b = {b_log:.6f}, a = {np.exp(log_a):.6f}")
    print(f"  Fit (iii) linear NLS a*N^b + c:   b = {popt_nlsc[1]:.6f}, a = {popt_nlsc[0]:.6f}, c = {popt_nlsc[2]:.4f}")
    print(f"  Fit (iv) rel-weighted NLS:        b = {popt_rel[1]:.6f}, a = {popt_rel[0]:.6f}")
    print()
