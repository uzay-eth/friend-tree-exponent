#!/usr/bin/env python3
"""Fit exact c_n = E[C_n] (n<=21) to A*G(n,mu) + B + (further terms), G(n,a)=Gamma(n+a)/Gamma(n), mu fixed."""
import sys
from fractions import Fraction
from mpmath import mp, mpf, gamma, matrix, lu_solve
mp.dps = 40
c = {}
for line in open(sys.argv[1]):
    f = line.split()
    if f and f[0].isdigit():
        q = Fraction(f[2]); c[int(f[0])] = mpf(q.numerator)/mpf(q.denominator)
mu = mpf(sys.argv[2]) if len(sys.argv) > 2 else mpf('0.565723')
G = lambda n, a: gamma(n + a)/gamma(n)
def fit(basis, ns):
    # exact solve using the last len(basis) points of ns
    pts = ns[-len(basis):]
    Am = matrix([[b(n) for b in basis] for n in pts]); rhs = matrix([c[n] for n in pts])
    return lu_solve(Am, rhs)
models = {
 "A G(mu) + B": [lambda n: G(n, mu), lambda n: mpf(1)],
 "A G(mu) + B + D G(mu-1)": [lambda n: G(n, mu), lambda n: mpf(1), lambda n: G(n, mu-1)],
 "A G(mu) + B + D G(mu-1) + E/n": [lambda n: G(n, mu), lambda n: mpf(1), lambda n: G(n, mu-1), lambda n: 1/mpf(n)],
 "A G(mu) + D G(mu-1) (no constant)": [lambda n: G(n, mu), lambda n: G(n, mu-1)],
 "A G(mu) + D G(mu-1) + F G(mu-2) (no constant)": [lambda n: G(n, mu), lambda n: G(n, mu-1), lambda n: G(n, mu-2)],
}
for name, basis in models.items():
    print(name)
    for last in (15, 17, 19, 21):
        ns = list(range(3, last + 1))
        sol = fit(basis, ns)
        print("  last n =", last, " coeffs:", [mp.nstr(x, 9) for x in sol])
