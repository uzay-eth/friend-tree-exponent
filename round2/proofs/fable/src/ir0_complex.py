#!/usr/bin/env python3
"""Complex roots of the IR-0 characteristic equation H_s(2,1) = 1 (floating point, exploratory; not a proof).

IR-0 = the friend tree in which a ring that picks a core neighbour is discarded.  One core vertex is then a Markov
chain on (d, c): (d,c) -> (d+1,c) at rate L = d - c and (d,c) -> (d,c+1) at rate b = L/d (birth of a newborn in state
(2,1)).  H_s(d,c) = E int e^{-st} b dt solves (s d/L + d + 1) H(d,c) = 1 + d H(d+1,c) + H(d,c+1), H(d,d) = 0; for d
beyond the cut the asymptotic form 1/s - [c/(s+1) + 1/(s+1)^2]/d is used.  The exponent mu_0 is the real root of
H_s(2,1) = 1; this script scans -0.9 < Re s < 0.6, 0 <= Im s <= 40 for further zeros and refines them by Newton steps.
Usage: ir0_complex.py   (about 4 minutes on one core)
"""
import sys
import numpy as np

def H21(s, Dmax=6000, cmax=70):
    # asymptotic boundary for d >= Dmax: H = 1/s - [c/(s+1) + 1/(s+1)^2]/d
    def asym(d, c): return 1/s - (c/(s+1) + 1/(s+1)**2)/d
    Hnext = np.array([asym(d, cmax + 1) for d in range(1, Dmax + 2)], dtype=complex); Hnext = np.concatenate([[0j], Hnext])   # H(., c+1), indexed by d
    for c in range(cmax, 0, -1):
        H = np.zeros(Dmax + 2, dtype=complex); H[Dmax + 1] = asym(Dmax + 1, c)
        for d in range(Dmax, c, -1):
            L = d - c
            H[d] = (1 + d * H[d + 1] + (Hnext[d] if d > c + 1 else 0.0)) / (s * d / L + d + 1)   # H(d,c+1)=0 when d=c+1
        Hnext = H
    return Hnext[2]
if __name__ == '__main__':
    print("check real root: H(0.4822669183233687) - 1 =", H21(0.4822669183233687 + 0j) - 1)
    print("boundary sensitivity:", H21(0.3 + 2j, 3000, 50), H21(0.3 + 2j, 6000, 70))

def H21v(S, Dmax=4000, cmax=60):
    S = np.asarray(S, dtype=complex)
    def asym(d, c): return 1/S - (c/(S+1) + 1/(S+1)**2)/d
    Hnext = [None] + [asym(d, cmax + 1) for d in range(1, Dmax + 2)]
    for c in range(cmax, 0, -1):
        H = [None] * (Dmax + 2); H[Dmax + 1] = asym(Dmax + 1, c)
        for d in range(Dmax, c, -1):
            L = d - c
            H[d] = (1 + d * H[d + 1] + (Hnext[d] if d > c + 1 else 0.0)) / (S * d / L + d + 1)
        Hnext = H
    return Hnext[2]
def scan():
    re = np.arange(-0.90, 0.6001, 0.025); im = np.arange(0.0, 40.001, 0.05)
    RE, IM = np.meshgrid(re, im, indexing='ij'); S = RE + 1j * IM
    F = H21v(S.ravel(), 2500, 45).reshape(S.shape) - 1
    # candidate zeros: cells where Re F and Im F both change sign
    cand = []
    for i in range(len(re) - 1):
        for j in range(len(im) - 1):
            q = F[i:i+2, j:j+2]
            if q.real.min() < 0 < q.real.max() and q.imag.min() <= 0 <= q.imag.max(): cand.append(S[i, j] + 0.0125 + 0.025j)
    roots = []
    for s0 in cand:
        s = s0
        for it in range(40):
            h = 1e-5; f = H21(s) - 1; df = (H21(s + h) - H21(s - h)) / (2 * h)
            step = f / df; s = s - step
            if abs(step) < 1e-11: break
        if abs(H21(s) - 1) < 1e-7 and not any(abs(s - r) < 1e-6 for r in roots): roots.append(s)
    roots.sort(key=lambda z: -z.real)
    print("zeros of H_s(2,1) - 1 found in -0.9 < Re s < 0.6, 0 <= Im s <= 40 (and complex conjugates):")
    for r in roots: print(f"   s = {r.real:+.8f} {r.imag:+.8f} i    |H-1| = {abs(H21(r) - 1):.1e}")
    print("min |F| on the lines Re s = 0.05..0.55, Im s >= 0.5:", {round(float(re[i]), 3): float(np.abs(F[i, 10:]).min()) for i in range(len(re)) if 0.04 < re[i] < 0.56 and i % 4 == 0})
if __name__ == '__main__': scan()
