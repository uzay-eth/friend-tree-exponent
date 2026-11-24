"""Small independent rational checks of the endpoint stencil and broad interpolation."""
from fractions import Fraction as F
from math import ceil, floor
from pathlib import Path
import re

from physical_trigger_bounds import D, UNIT, table


def compare_rows():
    lam = F(1133, 2000) + F(1, 3)
    toll = [12345 * (d + 2) for d in range(D + 1)]
    tail = lambda c: F(UNIT, 5 + c)
    checks = 0
    for parent in (0, 7):
        for upper in (False, True):
            directed = ceil if upper else floor
            extrema = max if upper else min
            direct = {(D, c): directed(tail(c)) for c in range(1, D + 1)}
            actual = table(lam, parent, toll, None, tail, upper)
            for d in range(D - 1, 1, -1):
                for c in range(d, 0, -1):
                    b = F(d - c, d)
                    gain = d - c
                    cap = F(c - 1, 2) + (F(1, parent) if parent else 0)
                    continuation = direct.get((d, c + 1), 0) + toll[d]
                    value = extrema((b * continuation + (gain + extra) * direct[d + 1, c])
                                    / (lam + b + gain + extra) for extra in (0, cap))
                    direct[d, c] = directed(value)
                    assert actual[d][c] == direct[d, c]
                    checks += 1
    print('independent Fraction endpoint comparisons', checks)


def broad_interpolation():
    logs = Path(__file__).resolve().parents[2] / 'outputs'
    critical = (logs / 'finite-physical-trigger-nineteen.log').read_text()
    root = (logs / 'finite-physical-trigger-final-bounds.log').read_text()
    upper = lambda pattern, data: int(re.search(pattern, data).group(1))
    anchor = F(upper(r'anchor at .462, k 4 <= (\d+)', root), 10**6)
    crt = F(upper(r'Perron-norm contraction <= (\d+)', critical), UNIT)
    assert anchor <= F('3.489') and crt <= F('.970148')
    assert upper(r'blocks 2=(\d+)', critical) >= 126364920098
    for s, target in (('.56564', '.980421'), ('.565', '.988137'),
                      ('.5649', '.989348'), ('.5642', '.997867'), ('.5641', '.999090')):
        p = (F('.5665') - F(s)) / (F('.5665') - F('.462'))
        numerator, denominator = p.as_integer_ratio()
        assert (F('3.489')**numerator * F('.970148')**(denominator - numerator)
                < F(target)**denominator)
        print('using only [.565,.5665] at', s, 'p <=', p, 'complement norm <', target)
    low = (logs / 'finite-physical-trigger-sector-six.log').read_text()
    _, hi = map(int, re.search(r'scalar in units 1e-12 (\d+) (\d+)', low).groups())
    assert UNIT - hi == 746245301832
    print('H=6, K=3 critical missing return >=', UNIT - hi, '/ 1e12')


def size_bound():
    # Replay the two endpoint inequalities including the child-copy gain.
    damp, q, delta = F('.18'), F('.268'), F('.11764')
    tail = F(ceil(10**9 / ((1 + damp*D)*q)), 10**9)
    h = {(D, c): tail for c in range(1, D+1)}
    for d in range(D-1, 1, -1):
        for c in range(d, 0, -1):
            L, b = d-c, F(d-c, d)
            vals = [((L+e)*h[d+1, c] + b*(h.get((d,c+1), 0) + 1/(1+damp*d)))
                    / (q + L + e + b) for e in (0, F(c+1, 2)+damp)]
            h[d, c] = F(ceil(10**9*max(F(1,4), *vals)), 10**9)
    h[2, 1] = F(1)
    for d in range(2, D):
        for c in range(1, d+1):
            for e in (0, F(c+1, 2)+damp):
                assert ((d-c+e)*(h[d+1,c]-h[d,c])
                        + F(d-c,d)*(h.get((d,c+1),0) + 1/(1+damp*d) - h[d,c])
                        <= q*h[d,c])
    hmin = min(h.values())
    assert hmin == F(81522947,125000000) and q*tail >= 1/(1+damp*D)
    assert delta*q*max(h.values())/(F(1,20)+damp) <= hmin
    assert damp + q + delta == F('.56564')

    log = (Path(__file__).resolve().parents[2] / 'outputs' /
           'finite-physical-trigger-nineteen-fine.log').read_text()
    m, theta = F('.12510370'), F('.976434')
    assert F(int(re.search(r'blocks 2=(\d+)', log).group(1)), UNIT) >= m
    assert F(int(re.search(r'Perron-norm contraction <= (\d+)', log).group(1)), UNIT) <= theta
    birth_tail = 1 + F('.006')/(1-theta)
    assert birth_tail < F('1.255')
    assert birth_tail*theta**900 < F('5.986e-10')
    assert birth_tail*theta**1024 < F('3.111e-11')
    Csize = birth_tail/(m*hmin*(F(1,20)+damp)*(1-theta))
    assert Csize < 2837
    print('189 degree pairs x 2 endpoints and high tail verified, size exponent', delta,
          'return size-exit coefficient < 2837')
    print('untriggered generation tail at .56564+: m=900 <5.986e-10, m=1024 <3.111e-11')


if __name__ == '__main__':
    compare_rows()
    broad_interpolation()
    size_bound()
