#!/usr/bin/env python3
"""Lead's independent recomputation for Astra claim 28 (four-reward route, c488e00).
Inputs (read from a clean archive of the Astra commit, path given as argv[1]): the authenticated Gold2 receipt (endpoints,
m_max, err_per_core, nominal) and the local upper table u(D, c) of host/conditioned24.json. Everything else is recomputed
here with my own code: the stopped discounted lifetime V(D, c) of B = m_max - u(D, c) on 2 <= D <= cut, the potential
ceiling M = ceil(m_max - lo V(2,1)), the lower endpoint lambda + r/M, the margin at 99/175, and the Farey statement.
The generator inequality is checked at EVERY rational gain rate on a grid of 61 points in [L, L + c/2] (not only at the
two ends), as an extra test of the affine-in-gain argument.
"""
import json, sys, hashlib
from fractions import Fraction as Q
from math import floor, ceil
from pathlib import Path
A = Path(sys.argv[1]).resolve()
rec = json.loads((A / 'round2/attacks/neon2-quintet/gold2/production.json').read_bytes())
cells = json.loads((A / 'round2/occupation-quintet-gold2/host/conditioned24.json').read_bytes())
print('receipt sha256', hashlib.sha256((A / 'round2/attacks/neon2-quintet/gold2/production.json').read_bytes()).hexdigest()[:16])
lo, hi = Q(rec['lower']['certified_exact']), Q(rec['upper']['certified_exact'])
Mstar, err, lam = rec['lower']['m_max'], rec['lower']['err_per_core'], Q(rec['lower']['nominal'])
u = {(d, c): int(cells['results']['lower']['cells'][f'{d},{c-1}'][1]) for d in range(2, 25) for c in range(1, d + 1)}
print('lo', lo, float(lo), ' hi', hi, float(hi), ' m_max', Mstar, ' err', err, ' nominal', lam)
floors = [207990, 395261, 52341, 251635]          # the four reward floors as stated (not recomputed here)
r = sum(floors) - err
print('sum of floors', sum(floors), ' net r', r)
for cut in (6, 18, 24):
    V = {}
    for D in range(cut, 1, -1):
        for c in range(D, 0, -1):
            L = D - c; B = Mstar - u[D, c]; assert B >= 0
            p = Q(L, D); Vg = V.get((D + 1, c), 0); Vp = V.get((D, c + 1), 0)
            best = None
            for g in (Q(L), Q(L) + Q(c, 2)):
                val = (B + g * Vg + p * Vp) / (hi + g + p)
                best = val if best is None or val < best else best
            V[D, c] = floor(best)
            for k in range(61):                       # generator inequality on a grid of gain rates
                g = Q(L) + Q(c, 2) * Q(k, 60)
                assert B + g * (Vg - V[D, c]) + p * (Vp - V[D, c]) - hi * V[D, c] >= 0, (cut, D, c, g)
    M = ceil(Mstar - lo * V[2, 1]); endpoint = lam + Q(r, M); margin = r - (Q(99, 175) - lam) * M
    print(f'cut {cut}: states {len(V)}, V(2,1) = {V[2,1]}, M = {M}, lower endpoint = {float(endpoint):.12f}, '
          f'margin at 99/175 = {float(margin):.1f} ticks/core; floors needed: r > {float((Q(99,175)-lam)*M):.1f}')
    for i, f in enumerate(floors):
        print(f'   without reward {i+1} ({f}): margin {float(margin - f):.1f}')
a, b = Q(99, 175), Q(142, 251)
print('Farey determinant', b.numerator * a.denominator - b.denominator * a.numerator, '; 142/251 - hi =', float(b - hi))
