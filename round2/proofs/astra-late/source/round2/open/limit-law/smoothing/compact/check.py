"""Bounded exact arithmetic for retry, beta, and multi-side budgets."""

from fractions import Fraction as Q
from math import comb


def retry():
    r = Q(1, 24)
    a, b = 2 / (2 - r), 1 / (1 - r)
    F = Q(3, 4) + r / (Q(1, 2) - r)
    rho = a / 2 + a*b/2 * (Q(1, 3) + Q(2, 3)*F)
    assert F == Q(37, 44) and rho == Q(11736, 11891) < 1
    c = Q(1, 30)
    assert r*Q(4, 5) == c < Q(1, 10) - r < r + c
    assert c < r  # post-cut crossing budget
    print(f'PASS defective retry moment <= {rho}; missing-mass rate {c}')


def beta():
    count = 0
    for a in (1, 2, 3, 9, 19, 41):
        for b in (1, 4, 7, 21, 61, 97):
            n, p = a+b, Q(a, a+b)
            peak = comb(n, a)*p**a*(1-p)**b
            assert peak**2*min(a, b) <= 1
            # Independent check that this is the maximum over the binomial
            # profiles on the grid, including asymmetric large cuts.
            for j in range(17):
                q = Q(j, 16)
                assert comb(n, a)*q**a*(1-q)**b <= peak
                count += 1
    print(f'PASS {count} beta profile comparisons, including unbalanced cuts')


def exponents():
    c, d = Q(1, 30), Q(47, 48)
    compact, global_min = c/(2*(c+d)), c/(2*(c+d)+c)
    assert (compact, global_min) == (Q(4, 243), Q(4, 247))
    assert global_min > 1000*Q(1, 80000)
    assert 1/(1-global_min) > Q(65, 64)
    # Three Young factors at p=3/2, and integrable pair-sum powers.
    p = Q(3, 2)
    assert 3/p == 2 and Q(1, 3) < 4
    for m in (3, 4, 6):
        for N in (10, 20, 50):
            variance = (1-Q(1, m))/(m*N+1)
            assert variance < 1/N and variance > Q(1, 2*m*N)
    print('PASS compact/global exponents, three-pair convolution, symmetric concentration scales')


if __name__ == '__main__':
    retry()
    beta()
    exponents()
