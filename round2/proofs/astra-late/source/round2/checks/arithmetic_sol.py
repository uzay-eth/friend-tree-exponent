#!/usr/bin/env python3
"""Exact cold controls for the original trigger arithmetic audit."""
from collections import Counter
from fractions import Fraction as Q
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'finite'))
import residual_entries as physical

physical.MARK_CAP = 20
ALPHA = Q(31, 20)


def exclusions():
    archive = (ROOT / 'outputs/occupation-affine-verified.txt').read_text()
    lo, hi = [Q(re.search(rf'^{side} fraction (\d+/\d+) ', archive, re.M)[1])
              for side in ('lower', 'upper')]
    a, b = Q(13, 23), Q(30, 53)
    assert a < lo <= Q(43, 76) <= hi < b
    assert b.numerator*a.denominator - a.numerator*b.denominator == 1
    # Different from the author's u,v search: test the discriminant directly.
    hits = []
    for q in range(1, 57201):
        p0 = -(-lo.numerator*q//lo.denominator)
        p1 = hi.numerator*q//hi.denominator
        for p in range(p0, p1+1):
            if gcd(p, q) == 1:
                d = p*p + 4*q*q
                if isqrt(d)**2 == d:
                    hits.append(Q(p, q))
    assert hits == [Q(32361, 57200)]
    assert Q(275, 208) - Q(208, 275) == hits[0]
    assert Q(2603, 4602) < lo and 2603**2 + 4*4602**2 == 9565**2

    # Any root in the interval lies on a monotone piece, split at the vertex.
    quadratics = []
    for a in range(1, 8):
        for b in range(-7, 8):
            for c in range(-7, 8):
                d = b*b - 4*a*c
                if gcd(a, gcd(b, c)) != 1 or d <= 0 or isqrt(d)**2 == d:
                    continue
                points = sorted({lo, hi, min(hi, max(lo, Q(-b, 2*a)))})
                vals = [a*x*x+b*x+c for x in points]
                if any(x*y <= 0 for x, y in zip(vals, vals[1:])):
                    quadratics.append((a, b, c))
    assert quadratics == [(3, -7, 3)]
    print('PASS: exact archived endpoints; Farey denominator 76, quadratic height 7, direct discriminant denominator 57200')
    return hi


def row(t):
    n, degree, _ = physical.statistics(t)
    p = Q(1, degree)
    rates = Counter(dict(physical.intrinsic_moves(t)))
    assert sum(rates.values()) == n-p
    for ch in t[2]:
        rates[ch] += physical.mark(ch[0])-p
        rates[physical.plus(ch)] += p
    if degree > 2:
        rates[physical.newborn(min(20, degree))] += t[1]*p
    assert all(r >= 0 for r in rates.values())
    return Q(1, t[0])+n-p, +rates


def path_checks(hi):
    for label in (2, 20):
        j = physical.newborn(label)
        stars = [(label, d-1, ()) for d in range(2, 21)]
        assert stars[0] == j
        assert [row(t)[1][u] for t, u in zip(stars, stars[1:])] == list(range(1, 19))
        assert row(stars[-1])[1][physical.newborn(20)] == Q(19, 20)
        assert [row(t)[0] for t in stars] == [Q(1, label)+d-Q(1, d) for d in range(2, 21)]
        assert all(row(t)[0] > ALPHA for t in stars[1:])
    # A concrete leading coefficient for each repeated-cycle witness.
    def cycle_coeff(i):
        return Q(19, 20)*prod(Q(d-1, 1)/(Q(1, i)+d-Q(1, d)-ALPHA)
                              for d in range(2 if i == 2 else 3, 21))
    assert cycle_coeff(2)*cycle_coeff(20) > 0
    for m in (0, 18, 19, 37, 38, 85, 150):
        n = max(20, m+2)
        exact = Q(1, 2)/(hi+ALPHA) * Q(n-1, n)/(hi+Q(1, 2)+n-Q(1, n))
        exact *= prod(Q(d-1)/(hi+Q(1, 2)+d-Q(1, d)) for d in range(2, n))
        assert n-1 > m and exact > Q(1, n*n*(n+1)*(n+2))
    print('PASS: capped-generator reset edges, positive witness Laurent coefficients, 7 exact omitted star paths')


def times_linear(poly, a):
    return [a*x+y for x, y in zip(poly+[0], [0]+poly)]


def divide_linear(poly, a):
    quot = [Q(0)]*(len(poly)-1)
    carry = poly[-1]
    for j in range(len(poly)-2, -1, -1):
        quot[j] = carry
        carry = poly[j]-a*carry
    return quot, carry


def reduced(terms):
    common = Counter()
    for _, factors in terms:
        common |= factors
    den = [Q(1)]
    for a in common.elements():
        den = times_linear(den, a)
    num = den[:]
    for coeff, factors in terms:
        term = [coeff]
        for a in (common-factors).elements():
            term = times_linear(term, a)
        for j, value in enumerate(term):
            num[j] -= value
    kept = common.copy()
    for a in common.elements():
        quot, rem = divide_linear(num, a)
        if rem == 0:
            num = quot
            kept[a] -= 1
    c = lcm(*(x.denominator for x in num))
    prim = [int(c*x) for x in num]
    assert num[-1] == 1 and gcd(*prim) == 1 and prim[-1] == c >= 1
    return prim, +kept


def reduction_checks():
    # Exact finite prefixes of the actual generator. These small depths do
    # not reach J20; the arbitrary-depth order is a proof, not a test.
    frontier = [(physical.newborn(2), Q(1), Counter())]
    terms = []
    for depth in range(5):
        nxt = []
        for t, weight, visited in frontier:
            a, edges = row(t)
            factors = visited + Counter([a])
            if t == physical.newborn(t[0]):
                terms.append((weight/2, factors))
            nxt.extend((u, weight*r, factors) for u, r in edges.items())
        frontier = nxt
        poly, poles = reduced(terms)
        smallest = min(min(fs) for _, fs in terms)
        assert min(poles) == smallest
        assert poles[smallest] == max(fs[smallest] for _, fs in terms)
        for s in (Q(43, 76), Q(3, 5)):
            tm = sum(w/prod((s+a) for a in fs.elements()) for w, fs in terms)
            assert sum(x*s**j for j, x in enumerate(poly)) == poly[-1]*prod(
                s+a for a in poles.elements())*(1-tm)
        if depth == 0:
            assert poly == [3, 2]
    # Positive histories can cancel a nonminimal pole: 1/[(s+a)(s+b)]
    # + 1/[(b-a)(s+b)] = 1/[(b-a)(s+a)].
    _, poles = reduced([(Q(1), Counter([2, 3])), (Q(1), Counter([3]))])
    assert poles == Counter({2: 1})
    # Positive toy histories with arbitrarily large evaluated cancellation,
    # even at a fixed positive rational input and with primitive polynomials.
    for m in (2, 20, 10**50):
        poly, poles = reduced([(Q(1, m), Counter([2])), (Q(9, 2)-Q(3, m), Counter([2, 2]))])
        assert poly[-1] == m and poles == Counter({2: 2})
        assert sum(poly) == 9*m//2 and Q(sum(poly), 9*m) == Q(1, 2)
    print('PASS: 5 physical prefix reductions; nonminimal-pole and scalar-gcd adversarial examples')


if __name__ == '__main__':
    path_checks(exclusions())
    reduction_checks()
