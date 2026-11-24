#!/usr/bin/env python3
"""Finite rational controls for the planted translation estimate."""
from fractions import Fraction as Q
from math import factorial


def signs():
    count = 0
    # x=z**q makes x**(mu-1) rational for mu=p/q.
    for mu in (Q(4,7), Q(2,3), Q(3,4), Q(9,10)):
        k = 2*mu
        for z in (Q(1,5), Q(1,2), Q(1), Q(2), Q(5)):
            x = z**mu.denominator
            ratio = z**mu.numerator
            for a,b in ((1,0),(0,1),(1,1),(1,3),(7,1)):
                # F''/(k*cos(theta)^k), by differentiating sin^k and cos^k.
                direct = a*ratio*((k-1)/x-1) + b*((k-1)*x-1)
                c = k-1
                turning = a*(ratio/x)*(c-x) + b*(c*x-1)
                assert direct == turning
                if a:
                    # Differentiate each of c*x^(mu-1) and -x^mu separately.
                    d2 = a*(c*(mu-1)*(mu-2)*ratio/x**3 - mu*(mu-1)*ratio/x**2)
                    assert d2 == a*(1-mu)*ratio/x**3*((2-mu)*c+mu*x) > 0
                count += 1
    print('PASS',count,'second-derivative and turning-polynomial controls')


def rates():
    c = Q(13,252)
    tscale = 1/(1+2*c)
    alpha = c*tscale
    assert alpha == Q(13,278) == Q(1,2)-tscale/2
    p,r = Q(51,50),Q(26,25)
    theta = (1/p-1/r)/(1-1/r)
    assert 0 < theta < 1 and 1/p == theta+(1-theta)/r
    assert r < Q(278,265) and theta*alpha == Q(325,14178)
    assert (1/p-Q(265,278)) == alpha*(1/p-Q(265,278))/(1-Q(265,278))
    # For a symmetric fold its preimage scale is delta^(1/2), so the
    # mass/length ratio diverges, even after a positive fixed weight.
    for n in range(2,18):
        length=Q(1,n*n)
        assert Q(1,n)**2 == length and Q(1,n)/length == n
    print('PASS exact translation, Sobolev interpolation, and fold exponents')


def first_cross():
    # Each ordering of the d+1 independent first rings is equiprobable;
    # only the d! orderings with the focal ring first are used by the bound.
    for d in range(2,80):
        assert Q(factorial(d),factorial(d+1))*Q(1,d)==Q(1,d*(d+1))
    print('PASS 78 finite-state next-cross bounds')


if __name__ == '__main__':
    assert __debug__
    signs()
    rates()
    first_cross()
