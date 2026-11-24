"""Independent arithmetic controls for deterministic-time compact smoothing."""
from fractions import Fraction as Q
from math import comb, factorial


def retry(r, failure):
    return (1 + (Q(1, 3) + Q(2, 3) * failure) / (1 - r)) / (2 - r)


def extinction():
    for k in range(101):
        w = Q(k, 100)
        q = w / (4 - 2*w)
        minimum = (Q(1, 2) + q)**2 - Q(1, 4)
        assert minimum == q + q*q
        assert minimum <= w/4 + 3*w*w/16 + 5*w**3/16 <= 3*w/4 <= w
    r = Q(1, 24)
    f = Q(3, 4) + r / (Q(1, 2) - r)
    assert f == Q(37, 44) and retry(r, f) == Q(11736, 11891) < 1
    r = Q(1, 16)
    f = Q(3, 4) + r * sum(q / (rate - r) for q, rate in (
        (Q(1, 4), Q(1, 2)), (Q(3, 16), Q(1)), (Q(5, 16), Q(3, 2))))
    assert f == Q(1307, 1610) and retry(r, f) == Q(373424, 374325) < 1


def beta_modes():
    checks = 0
    for n in (4, 8, 16, 40, 80, 128):
        for a in range(2, n-1):
            b = n - a
            mode = Q(a-1, n-2)
            density = Q(factorial(n-1), factorial(a-1)*factorial(b-1)) * mode**(a-1) * (1-mode)**(b-1)
            binomial = comb(n-2, a-1) * mode**(a-1) * (1-mode)**(b-1)
            assert density == (n-1)*binomial
            assert density**2 * (min(a,b)-1) <= n*n
            checks += 1
    return checks


def balances():
    for r, gamma, c, d, denom in (
        (Q(1, 24), Q(1, 24), Q(1, 30), Q(47, 48), Q(243, 4)),
        (Q(1, 16), Q(1, 20), Q(1, 20), Q(39, 40), Q(41)),
    ):
        a = Q(4, 5)
        assert c == min(a*r, (1-a)/2-gamma, gamma) and d == 1-gamma/2
        assert Q(2)*(c+d)/c == denom
        for mu in (Q(1, 4), Q(113, 200), Q(1133, 2000), Q(999, 1000)):
            xi = 1/(denom+mu)
            assert Q(1, 2) - d*xi/c - mu*xi/2 == xi
            assert (Q(247, 243) if denom > 50 else Q(42, 41)) < 1/(1-xi)
    assert Q(65, 64) < Q(247, 243)


if __name__ == '__main__':
    extinction()
    print('PASS 101 defective minimum-tail controls and both exact retry budgets')
    print(f'PASS {beta_modes()} integer Beta mode and stock-aware bounds')
    balances()
    print('PASS original and improved compact/global balance and strong p endpoints')
