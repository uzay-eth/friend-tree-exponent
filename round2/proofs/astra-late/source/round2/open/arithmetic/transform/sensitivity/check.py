#!/usr/bin/env python3
"""Sharp interval arithmetic, independent of the old three-credit lower bound."""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import gcd, isqrt, lcm
from pathlib import Path
import runpy


GOLD2_LO = Q(35083855339718511, 62017003637500000)
GOLD2_HI = Q(127993535256830043, 226242120050000000)
ROOT_TEST = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'check.py'))
# Reuse only exact polynomial helpers, never the frozen writer's interval.
excludes = ROOT_TEST['excludes']
value = ROOT_TEST['value']
irreducible_prime = ROOT_TEST['irreducible_prime']
WITNESSES = {
    2: (35, -63, 2),
    3: (-15, 15, 13, 13),
    4: (-5, 6, 3, 3, 1),
    5: (-3, 2, 3, 3, 3, 1),
}


def possible_polynomials(degree, height, lo, hi):
    den = lcm(lo.denominator, hi.denominator)
    a, b = (int(x*den) for x in (lo, hi))
    ai = [a**i * den**(degree-i) for i in range(1, degree+1)]
    bi = [b**i * den**(degree-i) for i in range(1, degree+1)]
    scale = den**degree
    for lead in range(1, height+1):
        for suffix in product(range(-height, height+1), repeat=degree-1):
            coeff = suffix + (lead,)
            top = sum(c*(y if c >= 0 else x) for c,x,y in zip(coeff, ai, bi))
            bot = sum(c*(x if c >= 0 else y) for c,x,y in zip(coeff, ai, bi))
            for j in range(max(-height, -((-bot)//scale)), min(height, top//scale)+1):
                yield (-j,) + coeff


def sharp_height(degree, ceiling, lo, hi):
    candidates = []
    for f in possible_polynomials(degree, ceiling, lo, hi):
        if gcd(*f) == 1 and not excludes(f, lo, hi):
            candidates.append((max(map(abs, f)), f))
    if not candidates:
        return None
    height, f = min(candidates)
    # This certifies sharpness if one among the lowest unresolved height is
    # a true root witness. Do not silently discard inconclusive inputs.
    for h, witness in sorted(candidates):
        if h != height:
            break
        if value(witness, lo)*value(witness, hi) < 0:
            try:
                return height, witness, irreducible_prime(witness)
            except AssertionError:
                continue
    raise AssertionError(('sharpness unresolved', degree, height, f))


def first_rational(lo, hi, cap, square):
    for q in range(1, cap+1):
        first = -((-lo.numerator*q)//lo.denominator)
        last = hi.numerator*q//hi.denominator
        for p in range(first, last+1):
            if gcd(p,q) != 1:
                continue
            delta = p*p+4*q*q
            if not square or isqrt(delta)**2 == delta:
                return Q(p,q)
    raise AssertionError(('no witness by denominator cap', cap, square))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', nargs=2, metavar=('LO', 'HI'), type=Q,
                        help='exact externally validated endpoints; provenance not checked here')
    parser.add_argument('--cap', type=int, default=94251)
    args = parser.parse_args()
    lo, hi = args.interval or (GOLD2_LO, GOLD2_HI)
    assert 0 < lo <= hi < 1
    if not args.interval:
        assert lo < Q(99,175) < hi
    print(f'Exact interval [{lo}, {hi}]')
    for degree, witness in WITNESSES.items():
        height = max(map(abs, witness))
        answer = sharp_height(degree, height, lo, hi)
        if answer:
            h, f, p = answer
            print(f'PASS degree {degree}: sharp height {h}, witness {f} mod {p}')
        else:
            print(f'PASS degree {degree}: height > {height}; sharpness unchecked')
    for square in (False, True):
        r = first_rational(lo, hi, args.cap, square)
        if not args.interval:
            assert r == (Q(53320,94251) if square else Q(99,175))
        print(f'PASS {"square-discriminant" if square else "unrestricted"} sharp denominator {r.denominator}: {r}')


if __name__ == '__main__':
    assert __debug__
    main()
