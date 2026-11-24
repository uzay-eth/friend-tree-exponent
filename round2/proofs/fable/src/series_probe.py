#!/usr/bin/env python3
"""Probe: how smooth is the exact sequence c_n = E[C_n] (n = 3..21)?  Local exponents and Neville tables."""
import sys
from fractions import Fraction
from mpmath import mp, mpf, log
mp.dps = 60
path = sys.argv[1]
c = {}
for line in open(path):
    f = line.split()
    if f and f[0].isdigit():
        c[int(f[0])] = Fraction(f[2])
ns = sorted(c)
# local exponent from one-step ratio: c_{n+1}/c_n = 1 + mu_n/n  (exact drift identity: mu_n = E[S_n]/E[C_n])
def F2m(q): return mpf(q.numerator)/mpf(q.denominator)
mu1 = {n: F2m((c[n+1]/c[n] - 1) * n) for n in ns[:-1]}
print("n   mu_n = n(c_{n+1}/c_n - 1)")
for n in ns[:-1]:
    print(n, mp.nstr(mu1[n], 20))
# Neville (polynomial in 1/n) extrapolation tables
def neville(seq, ns, shift=0.0):
    # seq[i] at x_i = 1/(n_i + shift); extrapolate to x=0
    T = [list(seq)]
    xs = [1/(mpf(n) + shift) for n in ns]
    k = 1
    while len(T[-1]) > 1:
        prev = T[-1]
        cur = []
        for i in range(len(prev) - 1):
            x0, xk = xs[i], xs[i + k]
            cur.append((x0 * prev[i + 1] - xk * prev[i]) / (x0 - xk))
        T.append(cur); k += 1
    return T
nn = ns[:-1]
for start in (3, 6, 9):
    sub = [n for n in nn if n >= start]
    T = neville([mu1[n] for n in sub], sub)
    print(f"\nNeville in 1/n from n >= {start}: last entries of each column")
    for k, col in enumerate(T):
        print(k, " ".join(mp.nstr(v, 12) for v in col[-3:]))
