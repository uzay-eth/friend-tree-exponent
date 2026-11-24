"""Independent rational/exponent controls for the all-threshold cutoff."""

from fractions import Fraction as Q


def wide():
    h = 512
    r = Q(7, 8)
    # Sum the polynomial envelope by derivatives of the geometric series.
    second = r * (1 + r) / (1 - r)**3
    first = r / (1 - r)**2
    assert 6 + second + 10 * first + 25 / (1 - r) < 2048
    assert Q(1, 8) + Q(9 + 2048, h*h) < Q(1, 7)
    assert Q(384, 191) * Q(22, 7)**2 / 6 < 4
    assert Q(6, 7 * 3**6) > Q(1, 1024)
    assert Q(1, 4096) - Q(24, 2**17) == Q(1, 16384)
    assert Q(113, 200) - Q(1, 2**21) - Q(2, 128) > Q(1, 2)
    assert Q(9, 28) + Q(1, 3) > Q(1, 2)  # gamma(4/3) lower
    assert Q(2**13*h*h, 2**39) + Q((2**7 + 2**14)*h*h, 2**78) < 1


def finite_read():
    # Worst-case shifts in units of delta. Lower read threshold is y+5δ,
    # upper is y−5δ, both allowed a further signed approximation of 2δ.
    low, high = 5, -5
    for gamma in (-2, 2):
        for predictor in (-1, 1):
            for rounding in (-2, 2):
                assert low + rounding - predictor >= gamma
                assert high + rounding - predictor <= gamma
                assert -8 <= high + rounding - predictor
                assert low + rounding - predictor <= 8
    assert 4 - 2 - 1 < 2  # the adversarial rounding of the original wording


def compact():
    a = 2**119
    b = a + 100
    # Upper/lower exponent certificates: quantities at level 2 never expanded.
    log_t, log_s = 3*b + 6, 3*b + 13
    log_factorial_index_lower = 5*b - 1
    assert log_factorial_index_lower > log_s + 3
    assert (42*b + 32 + b).bit_length() < log_factorial_index_lower
    assert (256*b).bit_length() < b - 4
    # log2 of stationary multiplier is <= 6BR+5B+1.
    log_log_multiplier = 5*b + (6*b + 2).bit_length()
    assert log_log_multiplier + 1 < 6*b - 3 < 7*b - 3
    assert (b + 12).bit_length() < 6*b - 4
    # Cost at most (1+t)/C up to constant, plus low stock and density bands.
    assert log_t + 2 - 8*b < -b
    assert 32*b + 3 - 60*b < -b
    assert 8*b + 3 - (16*b - 1) < -b
    assert 7*b + 1 < 12*b
    # The main margin uses disjoint budgets: stationary loss and clock errors.
    assert Q(1, 2**88) + Q(2, 16) < 1
    # Centered spectral coefficient < 2^976; clock factor <= 2^64 J.
    # Floors and +1 in k=R+1 absorbed in strict power comparisons.
    log2log2n = 16*a
    assert log2log2n == 2**123
    for rhs_upper_exponent in (358, 348 + 5*b + 10, 12*b + 12, log_t + 4):
        assert rhs_upper_exponent < log2log2n
    # A simple sufficient margin in stationary construction, in units 2^-d*.
    assert Q(113, 200) * (1-Q(1133, 2000)) / 256 > Q(1, 2048)
    assert 2**118 + 11 < a


if __name__ == '__main__':
    wide()
    print('PASS wide likelihood, beta separation and stationary flux budgets')
    finite_read()
    print('PASS 8 extreme readout shifts; original unqualified rounding fails')
    compact()
    print('PASS compact margins, signed-transfer exponents, log2(log2(N)) = 2^123')
