#!/usr/bin/env python3
"""Pool naive_chain 'fixed' outputs and compare with exact rationals (cross-check, Claude Opus 5.5).

Exact file columns used (identified cold from n=3,4,5 by hand): n, #trees, E[C_n], E[S_n], E[S_n]/E[C_n], Var(C_n).
For each n reports z-scores of the pooled sample mean of C, the pooled mean of S, and the pooled sample
variance of C (its SE estimated by the delta method from the 4th central moment is not available here,
so the variance check uses the normal-theory SE sqrt(2/(R-1)) * Var, flagged as approximate).
Usage: compare_fixed.py EXACT_FILE out1 [out2 ...]
"""
import sys, math
from fractions import Fraction

exact = {}
for line in open(sys.argv[1]):
    f = line.split()
    if len(f) >= 6 and f[0].isdigit():
        exact[int(f[0])] = tuple(float(Fraction(x)) for x in f[2:6])

pool = {}
for path in sys.argv[2:]:
    R = None
    for line in open(path):
        if line.startswith('#'):
            R = int(line.split('runs=')[1].split()[0]); continue
        f = line.split(); n = int(f[1]); d = dict(zip(f[2::2], map(float, f[3::2])))
        p = pool.setdefault(n, [0, 0.0, 0.0, 0.0, 0.0])
        for k, (m, sd) in enumerate(((d['C_mean'], d['C_sd']), (d['S_mean'], d['S_sd']))):
            p[1 + 2 * k] += R * m
            p[2 + 2 * k] += (R - 1) * sd * sd + R * m * m
        p[0] += R

print(f"{'n':>3} {'runs':>10} {'meanC':>11} {'exactC':>11} {'zC':>7} {'meanS':>11} {'exactS':>11} {'zS':>7} {'varC':>10} {'exactVar':>10} {'zVar~':>7}")
worst = 0.0; chi = 0.0; k = 0
for n in sorted(pool):
    if n not in exact: continue
    R, sc, ssc, ss, sss = pool[n]
    mc = sc / R; vc = (ssc - R * mc * mc) / (R - 1)
    ms = ss / R; vs = (sss - R * ms * ms) / (R - 1)
    eC, eS, _, eV = exact[n]
    zc = (mc - eC) / math.sqrt(vc / R) if vc > 0 else 0.0
    zs = (ms - eS) / math.sqrt(vs / R) if vs > 0 else 0.0
    zv = (vc - eV) / (eV * math.sqrt(2.0 / (R - 1))) if eV > 0 else 0.0
    if vc > 0:
        worst = max(worst, abs(zc), abs(zs)); chi += zc * zc; k += 1
    print(f"{n:>3} {R:>10} {mc:>11.6f} {eC:>11.6f} {zc:>7.2f} {ms:>11.6f} {eS:>11.6f} {zs:>7.2f} {vc:>10.5f} {eV:>10.5f} {zv:>7.2f}")
print(f"max |z| over mean C and mean S (n with nonzero variance): {worst:.2f}; sum of zC^2 = {chi:.2f} over {k} sizes (correlated across n)")
