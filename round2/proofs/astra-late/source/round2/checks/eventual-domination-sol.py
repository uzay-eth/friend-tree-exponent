"""Exact algebra and endpoint mutations for the Claim 24 synthesis."""

from fractions import Fraction as Q
from math import comb, factorial


def rising(a, k):
    result = Q(1)
    for i in range(k):
        result *= a + i
    return result


def beta_over_pi(a, b):
    """B(a,b)/pi for positive half-odd-integer arguments by gamma values."""
    i, j = int(a - Q(1, 2)), int(b - Q(1, 2))
    assert a == i + Q(1, 2) and b == j + Q(1, 2) and min(i, j) >= 0
    gamma_i = Q(comb(2 * i, i) * factorial(i), 4**i)
    gamma_j = Q(comb(2 * j, j) * factorial(j), 4**j)
    return gamma_i * gamma_j / factorial(i + j)


def placement_checks():
    cases = 0
    for v in (3, 4, 5, 9, 16, 31):
        a, b = Q(1, 2), v - Q(1, 2)
        # Dividing normalizations leaves 2v; the kernel adds one theta.
        normalizer = beta_over_pi(a, b) / beta_over_pi(a + 1, b)
        assert normalizer == 2 * v and normalizer != 2 * (v - 1)
        for k in range(7):
            tilted = beta_over_pi(a + 1 + k, b) / beta_over_pi(a + 1, b)
            source = beta_over_pi(a + 1 + k, b) / beta_over_pi(a, b)
            # Use the actual n(n+1) and n placement denominators.
            n = v - 1
            assert tilted / (2 * n * (n + 1)) == source / n
            assert tilted / (2 * n**2) != source / n
            assert tilted == 2 * v * rising(a, k + 1) / rising(v, k + 1)
            cases += 1
    return cases


def power_checks():
    cases = 0
    for v in (3, 8, 31, 96):
        for alpha in (2, 3, 7, 23, 129):
            a, b = alpha - Q(1, 2), v + Q(1, 2)
            base = beta_over_pi(a, b)
            gap = beta_over_pi(a + 1, b)
            later = beta_over_pi(a, b + 1)
            assert base == gap + later
            assert gap / base == (alpha - Q(1, 2)) / (alpha + v)
            assert later / base == (v + Q(1, 2)) / (alpha + v)
            assert 0 < gap < base and 0 < later < base
            for red in (a, a + 1):
                inverse = beta_over_pi(red - 1, b) / (v * beta_over_pi(red, b))
                assert inverse == Q(1, v) + b / (v * (red - 1))
            cases += 1
    return cases


def log_correction_bounds(v, y):
    """Enclose log(exp(y)*(1-y/v)**(v-1/2)) by a rational series."""
    x = Q(y, v)
    terms = sum((x**k / k for k in range(1, 5)), Q(0))
    tail = x**5 / (5 * (1 - x))
    power = v - Q(1, 2)
    return y - power * (terms + tail), y - power * terms


def endpoint_checks():
    cases = 0
    for y in (8, 16, 32, 64, 128):
        # vx→∞ and x→0 alone do not make the exponential correction one.
        lo, hi = log_correction_bounds(y**2, y)
        assert Q(-3, 5) < lo <= hi < Q(-2, 5)
        lo, hi = log_correction_bounds(y**3, y)
        assert -Q(1, y) < lo <= hi < 0
        cases += 1
    for k in (4, 8, 16, 32, 64):
        # Near one, moving the v-1 into the exponent is not a small change.
        epsilon = Q(1, 2**k)
        v = k + 3
        exponent_error_squared = epsilon ** (2 * (v - 1) - 1) / epsilon ** (2 * v - 1)
        assert exponent_error_squared == 2 ** (2 * k)

        # Model relative one-label errors tending to zero: their difference
        # can reverse the signed gap in the small-threshold transition.
        x, error = Q(1, k**2), Q(1, k)
        assert (1 - x) * (1 + error) > 1
        # Away from zero, the same error cannot reverse the comparison.
        assert (1 - Q(1, 2)) * (1 + error) < 1
        cases += 1
    return cases


if __name__ == "__main__":
    print(f"PASS {placement_checks()} exact two-size placement/tilt checks")
    print(f"PASS {power_checks()} exact two-power-weight beta ratio and inverse-mean checks")
    print(f"PASS {endpoint_checks()} rational overlap and endpoint controls")
