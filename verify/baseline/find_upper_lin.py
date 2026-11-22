import sys, certs6
from fractions import Fraction as Fr
Lmax, cmax, K, Pmax = map(int, sys.argv[1:5]); lo, hi = float(sys.argv[5]), float(sys.argv[6])
M = certs6.Model('upper', Lmax, cmax, K, Pmax)
for _ in range(int(sys.argv[7]) if len(sys.argv) > 7 else 10):
    mid = (lo + hi) / 2
    out, rho = certs6.upper_lin_table(M, mid)
    lo, hi = (lo, mid) if out is not None else (mid, hi)
    print(f"  lam={mid:.5f} rho={rho:.5f}", flush=True)
print('upper', Lmax, cmax, K, Pmax, 'threshold in', (lo, hi), flush=True)
if len(sys.argv) > 8:
    lam = Fr(sys.argv[8])
    lam_build = float(sys.argv[9]) if len(sys.argv) > 9 else float(lam)   # optional: build at another lambda
    (f, G), rho = certs6.upper_lin_table(M, lam_build)
    snb, sup, n = certs6.verify(M, f, lam, G)
    x0 = (2, 0, tuple([0] * K), 0, 3)
    print(f"UPPER(K={K},Pmax={Pmax}) certificate OK: lambda={lam} (={float(lam):.4f}) Lmax={Lmax} cmax={cmax} "
          f"sup f={float(sup):.4f} f(T3 centre)={f[x0]:.4f} G={G:.4f}; {n} checks", flush=True)
    from mpmath import mp, gamma, mpf
    mp.dps = 30
    kap = Fr(f[x0])
    const = mpf(kap.numerator) / kap.denominator * 2 / gamma(3 + mpf(lam.numerator) / lam.denominator)
    print(f"  => E[C_n] <= {float(kap):.8f} * prod_(m=3..n-1)(1+lambda/m) = {float(const):.8f} * Gamma(n+lambda)/Gamma(n)  (n>=3)")
