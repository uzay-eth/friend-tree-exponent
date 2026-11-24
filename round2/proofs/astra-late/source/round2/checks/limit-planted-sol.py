"""Rational checks of angular slopes, weighted translation and cut balances."""
from fractions import Fraction as Q


def angular():
    rows = 0
    # atan(t) with rational t is sampled via x=tan(theta)^2; choose x so x^mu is rational.
    for mu in (Q(3, 5), Q(2, 3), Q(3, 4), Q(9, 10)):
        c = 2*mu-1
        for z in (Q(1, 10), Q(1, 2), Q(1), Q(2), Q(10)):
            x, power = z**mu.denominator, z**mu.numerator
            for a, b in ((0, 1), (1, 0), (1, 1), (1, 20), (20, 1)):
                f = a*power+b  # divided by cos(theta)^(2mu)
                slope = a*power/x-b  # additionally divide F' by 2mu*tan(theta)
                direct = a*power*(c/x-1)+b*(c*x-1)
                assert direct == -2*(1-mu)*f+c*(1-x)*slope
                if not a or not b:
                    assert (slope > 0) == bool(a)
                if a == b:
                    assert (slope > 0) == (x < 1)
                rows += 1
    return rows


def translate(xs, h):
    return {x+h: w for x, w in xs.items()}


def weighted():
    k = {0: Q(1, 3), 1: Q(1, 3), 2: Q(1, 3)}
    l = Q(2, 5)
    retained = {x: l*w for x, w in k.items()}
    translated = translate(retained, 1)
    assert translated == {x: l*w for x, w in translate(k, 1).items()}
    gated = {x: w for x, w in k.items() if x == 0}
    assert translate(gated, 1) == {1: Q(1, 3)}
    assert not {x: w for x, w in translate(k, 1).items() if x == 0}


def exponents():
    c = Q(13, 252)
    alpha = c/(1+2*c)
    assert alpha == Q(13, 278)
    rp = Q(278, 265)
    for p in (Q(1), Q(51, 50), Q(103, 100)):
        interior = (p+rp)/2
        theta = (1/p-1/interior)/(1-1/interior)
        assert (0 < theta <= 1) and (theta*alpha <= 1/p-Q(265, 278))
    for mu in (Q(3, 5), Q(7, 8)):
        # Equate exp(-ct) and n^a exp(t/2) sqrt(delta).
        a = (1+mu)/2
        size_cost = c*a/(c+Q(1, 2))
        assert size_cost == c*(1+mu)/(1+2*c)
        assert (1-c/(c+Q(1, 2)))*a == a/(1+2*c)
    for r in (Q(21, 20), Q(3, 2), Q(2)):
        assert 2/r-1 == (2-r)/r


if __name__ == '__main__':
    print(f'PASS {angular()} angular curvature and monotonicity controls including zero coefficients')
    weighted()
    print('PASS scalar-weight translation and adverse output-gate fixture')
    exponents()
    print('PASS Besov, strict Sobolev, and gated-side refinement balances')
