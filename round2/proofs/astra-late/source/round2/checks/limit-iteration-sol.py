"""Independent rational controls for angular and six-side smoothing."""
from fractions import Fraction as Q
from math import comb, factorial


def retry():
    r = Q(13, 200)
    weights = [Q(j+1, 2**(j+2)) for j in range(1, 9)]
    remain = Q(3, 4) - sum(weights)
    f = Q(3, 4) + r*(sum(w/(Q(j, 2)-r) for j, w in enumerate(weights, 1))
                       + remain/(Q(9, 2)-r))
    rho = (1+(1+2*f)/(3*(1-r)))/(2-r)
    assert f == Q(5715644171206736442597875, 7046497926887658264039168)
    assert rho == Q(11951004048194066431207745075, 11951906448537490802583997242) < 1
    for k in range(101):
        w = Q(k, 100)
        q = w/(4-2*w)
        assert q+q*q <= sum(v*w**j for j, v in enumerate(weights, 1))+remain*w**9


def angular():
    count = 0
    for n in (4, 5, 12, 32, 128):
        for a in range(2, n-1):
            b = n-a
            p = Q(2*a-1, 2*(n-1))
            # Squared angular mode; no fractional powers and no balanced-size assumption.
            peak2 = (2*(n-1)*comb(n-2, a-1))**2 * p**(2*a-1) * (1-p)**(2*b-1)
            assert peak2 <= 4*n
            count += 1
    c, alpha, r = Q(13, 252), Q(100, 126), Q(13, 200)
    assert c == r*alpha == (1-alpha)/2-c
    for mu in (Q(1, 2), Q(57, 100), Q(9999, 10000), Q(1)):
        xi = c/(1+2*c+mu*c)
        assert xi == 13/(278+13*mu)
        diff = 1/(1-xi) - Q(291, 278)
        assert diff == 169*(1-mu)/(278*(265+13*mu))
        assert (diff > 0) == (mu < 1)
    assert c/(1+2*c) == Q(13, 278)
    return count


def six_side():
    rows = 0
    for q in [Q(1, 3), Q(1), Q(2), Q(7, 2), Q(10)]:
        for l in range(5, 80):
            if l < 4+2*q:
                continue
            assert -q*l/(l-2)+q/(l-3-q) <= -q
            rows += 1
    # Integer-order analogue of inverse Gamma extraction and pair-sum factors:
    # E[P_i^-1], i=1..3, jointly under Dirichlet with integer totals.
    product_checks = 0
    for sizes in [(4, 4, 4), (4, 20, 300), (32, 45, 17)]:
        n = sum(sizes)
        den = factorial(n-4)
        moment = Q(factorial(n-1), den)
        for k in sizes:
            moment *= Q(factorial(k-2), factorial(k-1))
        assert moment == Q((n-1)*(n-2)*(n-3), (sizes[0]-1)*(sizes[1]-1)*(sizes[2]-1))
        product_checks += 1
    return rows, product_checks


if __name__ == '__main__':
    retry()
    print('PASS 101 defective-tail envelopes and exact retry moment')
    print(f'PASS {angular()} angular mode bounds and strict L^p endpoint arithmetic')
    a, b = six_side()
    print(f'PASS {a} stopped-stock drift checks and {b} inverse Dirichlet product controls')
