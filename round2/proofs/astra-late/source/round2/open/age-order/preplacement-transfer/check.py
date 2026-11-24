"""Exact scalar barriers and sufficient-condition arithmetic."""

from fractions import Fraction as Q
from functools import cache


@cache
def p(d, j, stop):
    """Upper bound on first send to a designated leaf."""
    if d == stop:
        return Q(1, j + 1)
    return (1 + d * (j + 1) * p(d + 1, j + 1, stop)
            + (j * p(d, j - 1, stop) if j else 0)) / (1 + d * (j + 1) + j)


@cache
def q(d, j, stop):
    """Lower bound on hubhood from j tagged leaves."""
    if not j:
        return Q(0)
    if d == stop:
        return 1 - Q(1, d)**j
    return (d * q(d + 1, j + 1, stop) + q(d, j - 1, stop)) / (d + 1)


def main():
    stop = 112
    rows = 0
    for d in range(2, stop):
        for j in range(d - 1):
            assert p(d, j, stop) >= p(d + 1, j + 1, stop)
            if j:
                assert p(d, j, stop) <= p(d, j - 1, stop)
            rows += 1
    lo = q(2, 1, 16)
    assert lo > Q(613, 1000)
    print(f"scalar drift signs checked at {rows:,} (degree, stock) rows")
    serial = p(2, 0, stop)
    assert serial < Q(611, 1000) and lo - serial > Q(29, 10000)
    m = 1024
    certificate = (Q(29, 10000) - Q(2, m)) * (1 - Q(3, m + 1)) - Q(611, 1000 * (m + 1))
    assert certificate > 0
    print(f"middle hub lower: {float(lo):.9f}; serial-tip upper: {float(serial):.9f}")
    print(f"uniform serial gap lower: {float(lo - serial):.9f}")
    print(f"1,024-leaf criterion: numerator per core > {float(certificate):.9f}")


if __name__ == "__main__":
    main()
