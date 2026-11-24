#!/usr/bin/env python3
"""Exact interval exclusions and a growing-denominator countermodel."""
from fractions import Fraction as Q
from itertools import product
from math import gcd, isqrt

LO = Q(126905728962161948969119589499374909399,
       224328301684105796116575000000000000000)
HI = Q(127993535256830043, 226242120050000000)


def heights():
    # Test polynomials with positive leading coefficient. Integer terms are
    # bounded by monotonicity of x**k; this encloses even nonmonotone polynomials.
    for degree, ceiling, witness in (
        (2, 63, (35, -63, 2)),
        (3, 15, (-15, 15, 13, 13)),
        (4, 6, (-5, 6, 3, 3, 1)),
        (5, 3, (-3, 2, 3, 3, 3, 1)),
    ):
        h = ceiling - 1
        # Common integer scale for both endpoint power lists.
        den = LO.denominator * HI.denominator
        xs = [int(x * den) for x in (LO, HI)]
        assert all(Q(i, den) == x for i, x in zip(xs, (LO, HI)))
        l, u = [[x**k * den**(degree-k) for k in range(1, degree+1)] for x in xs]
        scale = den**degree
        count = 0
        for lead in range(1, h+1):
            for inner in product(range(-h, h+1), repeat=degree-1):
                coeffs = inner + (lead,)
                upper = sum(c*(b if c >= 0 else a) for c, a, b in zip(coeffs, l, u))
                lower = sum(c*(a if c >= 0 else b) for c, a, b in zip(coeffs, l, u))
                # Need some -c0 in [-h,h] with -c0*scale in [lower,upper].
                for opposite in range(max(-h, -((-lower)//scale)), min(h, upper//scale)+1):
                    f = (-opposite,) + coeffs
                    if monotone(f, LO, HI):
                        assert value(f, LO)*value(f, HI) > 0, f
                    else:
                        assert excludes(f, LO, HI), f
                count += 1
        vals = [sum(c*x**i for i,c in enumerate(witness)) for x in (LO, HI)]
        assert vals[0]*vals[1] < 0 and max(map(abs, witness)) == ceiling
        print(f'PASS degree {degree}: all heights < {ceiling} excluded ({count} tests); root witness')
        p = irreducible_prime(witness)
        print(f'  irreducible witness verified modulo {p}')


def value(f, x):
    return sum(c*x**i for i,c in enumerate(f))


def bounds(f, lo, hi):
    lower = sum(c*(lo if c >= 0 else hi)**i for i,c in enumerate(f))
    upper = sum(c*(hi if c >= 0 else lo)**i for i,c in enumerate(f))
    return lower, upper


def monotone(f, lo, hi):
    a, b = bounds([i*c for i,c in enumerate(f) if i], lo, hi)
    return a >= 0 or b <= 0


def excludes(f, lo, hi, depth=12):
    a, b = bounds(f, lo, hi)
    if a > 0 or b < 0:
        return True
    if monotone(f, lo, hi):
        return value(f, lo)*value(f, hi) > 0
    m = (lo+hi)/2
    return depth and excludes(f, lo, m, depth-1) and excludes(f, m, hi, depth-1)


def irreducible_prime(f):
    degree = len(f)-1
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        if f[-1] % p == 0:
            continue
        def reducible():
            for k in range(1, degree//2+1):
                for g in product(range(p), repeat=k):
                    r = list(f)
                    for j in range(degree, k-1, -1):
                        a = r[j] % p
                        for i,b in enumerate(g):
                            r[j-k+i] -= a*b
                    if all(c % p == 0 for c in r[:k]):
                        return True
            return False
        if not reducible():
            return p
    raise AssertionError(f)


def discriminant():
    cap = 94251
    out = []
    for v in range(1, isqrt(cap)+1):
        for u in range(v+1, min(2*v, cap//v)+1):
            if gcd(u,v) != 1:
                continue
            s = Q(u*u-v*v,u*v)
            assert s.denominator == u*v
            if LO <= s <= HI:
                out.append((u,v,s))
    assert out == [(353, 267, Q(53320, 94251))]
    print('PASS sharp square-discriminant denominator 94251')


def countermodel():
    r = Q(241,426)
    assert LO < r < HI
    prefix = Q(1,2)/(r+2)
    c = 1-prefix
    value = prefix
    for n in range(1, 257):
        w = Q(1,n*n) - Q(1,(n+1)**2)
        value += c*w*(r+n+2)/(r+n+2)
        assert 1-value == c/(n+1)**2
        assert (1-value).denominator == c.denominator*(n+1)**2//gcd(c.numerator,(n+1)**2)
    print('PASS 256 rational-root truncations; exact reduced denominator grows quadratically')


def physical_tail():
    # Independent telescoping control for the star and reset lower bounds.
    for s in (LO, HI, Q(1)):
        prod = Q(1)
        for d in range(2, 257):
            alpha = d+Q(1,2)-Q(1,d)
            prod *= (d-1)/(s+alpha)
            m = d+1
            assert prod >= Q(6,(m-1)*m*(m+1))
            if m >= 20:
                reset = (1-Q(1,m))/(s+m+Q(1,2)-Q(1,m))
                end = Q(1,2)/(s+Q(31,20))
                assert prod*reset*end >= Q(1,m*m*(m+1)*(m+2))
    for n in range(20,257):
        assert sum(Q(1,m*m*(m+1)*(m+2)) for m in range(n+1,2*n+1)) >= Q(1,16*(n+1)**3)
    print('PASS physical star/reset tail lower envelopes (3 discounts; 237 blocks)')


if __name__ == '__main__':
    assert __debug__
    heights()
    discriminant()
    countermodel()
    physical_tail()
