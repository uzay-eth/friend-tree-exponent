#!/usr/bin/env python3
"""Exact interval exclusions and fixed-value pole-sign countermodel controls."""
from fractions import Fraction as Q
from math import comb, factorial, gcd, isqrt, lcm

LO = Q(191976090122981726641370633087, 339373453939523193250000000000)
HI = Q(24172191945911554486185483421, 42724375658564528593750000000)


def ceil(x):
    return -((-x.numerator)//x.denominator)


def rational():
    for q in range(1, 10001):
        p = ceil(q*LO)
        if p <= q*HI:
            return Q(p, q)
    raise AssertionError('search budget')


def square(limit):
    hits = []
    for v in range(1, isqrt(limit)+1):
        # Binary search for first numerator with u/v - v/u >= LO.
        start, end = v+1, min(2*v, limit//v)+1
        lo, hi = start, end
        while lo < hi:
            u = (lo+hi)//2
            if Q(u*u-v*v, u*v) < LO:
                lo = u+1
            else:
                hi = u
        for u in range(lo, end):
            z = Q(u*u-v*v, u*v)
            if z > HI:
                break
            if gcd(u, v) == 1:
                assert z.denominator == u*v
                hits.append(z)
    return min(hits, key=lambda z: z.denominator)


def quadratic(limit):
    hits = []
    for a in range(1, limit+1):
        for b in range(-limit, limit+1):
            points = [LO, HI]
            if LO < Q(-b, 2*a) < HI:
                points.append(Q(-b, 2*a))
            vals = [-a*x*x-b*x for x in points]
            for c in range(max(-limit, ceil(min(vals))), min(limit, max(vals)//1)+1):
                disc = b*b-4*a*c
                if gcd(gcd(a, b), c) == 1 and disc > 0 and isqrt(disc)**2 != disc:
                    hits.append((max(a, abs(b), abs(c)), a, b, c))
    return min(hits)


def prime(p):
    return p >= 2 and all(p % j for j in range(2, isqrt(p)+1))


def spike(p, order, n):
    return (Q(factorial(p-2)*(p-1), 2*p) * (-factorial(p))**(n-order)
            * comb(n-1, order-1)) if n >= order else Q(0)


def digits_at_rational(r):
    # Enclose the infinite spike value by a dyadic tail. Verified digits
    # give a finite crossing bracket for the exact greedy construction.
    P = 31
    base = 1+r
    center = Q(1, 2)/(r+2)
    low = center + sum(Q(factorial(p-2)*(p-1), 2*p)/(r+factorial(p))**p
                       for p in range(5, P+1) if prime(p))
    high = low + Q(1, 2**P)
    a, b = 1-high, 1-low
    acc, digits = Q(0), []
    for n in range(1, 16):
        assert 0 <= a <= b < 1
        a, b = a*base, b*base
        assert a//1 == b//1
        digit = a//1
        assert digit in (0, 1)
        digits.append(digit)
        a, b = a-digit, b-digit
        acc += digit/base**n
    assert low+acc < 1 < high+acc+base**(-len(digits))
    return digits


def main():
    assert Q('0.565677980686114') <= LO < HI <= Q('0.565770513279953')
    r, sq, q2 = rational(), square(500000), quadratic(100)
    print(f'least rational: {r}; least quadratic (H,a,b,c): {q2}; least square case: {sq}')
    assert r == Q(99, 175) and q2 == (7, 3, -7, 3) and sq == Q(32361, 57200)
    assert Q(2603, 4602) < LO < sq < HI
    assert 3*LO*LO-7*LO+3 > 0 > 3*HI*HI-7*HI+3
    assert Q(56, 99) < LO <= r <= HI < Q(43, 76) and 43*99 - 56*76 == 1
    assert lcm(*range(1, 21))*(LO-Q('0.5641')) > 2
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        for shift, denom in ((0, 2), (-1, 1)):
            n = p+shift
            # Integer-coefficient digit correction, with an integer pole.
            regular = Q((-2)**(n-1), 2) + sum((-1)**(n-j)*comb(n-1, j-1)
                                             for j in range(1, n+1, 2))
            a = regular + (2//denom)*sum(spike(q, q+shift, n) for q in range(5, p+1) if prime(q))
            target = Q(factorial(p-2)*(p-1), denom*p)
            assert (a-target).denominator % p and a.denominator % p == 0
            assert (p*a).denominator % p
    print('PASS: 18 exact first-prime coefficient congruences with positive rational chains')
    print('rational crossing prefix digits:', ''.join(map(str, digits_at_rational(r))))

    # At positive rational trial discount a reciprocal degree p in a
    # holding denominator turns into a numerator zero modulo p.
    p = 23
    holding = 1/(r+Q(1, 2)+p-Q(1, p))
    assert holding.numerator % p == 0 and holding.denominator % p != 0
    assert (Q(p-1, p)*holding).denominator % p != 0
    print('PASS: exact p=23 holding-factor regularization at the rational trial')


if __name__ == '__main__':
    main()
