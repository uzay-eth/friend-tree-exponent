"""Exact final arithmetic for the cutoff in README.md, §3."""

from argparse import ArgumentParser
from fractions import Fraction as Q


L = Q(113, 200)
U = Q(1133, 2000)
BETA = (1 - U) * L / 2


def envelopes(c, eta, h_lo, h_hi, radius, degree_cap):
    if c <= 0 or eta <= 0 or not 0 < h_lo <= h_hi:
        raise ValueError("need C, eta > 0 and 0 < h_lo <= h_hi")
    if radius < 0 or degree_cap < 2:
        raise ValueError("need radius >= 0 and degree cap >= 2")
    j = 96 * (radius + 1) * (degree_cap + 1) ** (radius + 1) + 6
    return (
        (2 * c / h_lo, eta, Q(1)),
        (c * 4**radius * 6 * 32 * (h_hi + c) / h_lo, eta / 20, BETA / 8),
        (16 * j / h_lo, L - Q(21, 40), BETA / 8),
    )


def certifies(power, rows):
    return power >= 3 and all(
        coefficient < budget * 2 ** int(power * decay)
        for coefficient, decay, budget in rows
    )


def cutoff_power(rows):
    """Choose s with N0 = 2**s, using floor(s * decay) in each row."""
    result = 3
    for coefficient, decay, budget in rows:
        log_bound = 0
        while coefficient >= budget * 2**log_bound:
            log_bound += 1
        result = max(result, -(-log_bound // decay))
    return result


def self_check():
    assert BETA == Q(97971, 800000)
    assert L - Q(21, 40) == Q(1, 25)
    rows = envelopes(Q(2), Q(1), Q(1), Q(2), 1, 2)
    s = cutoff_power(rows)
    assert certifies(s, rows)
    assert not certifies(s - 1, rows)
    assert not certifies(150, [rows[-1]])
    assert cutoff_power([(Q(1), Q(2, 7), Q(1, 2))]) == 7
    print(f"rational arithmetic checks passed; fixture N0 = 2^{s}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    for name in ("c", "eta", "h-lo", "h-hi"):
        parser.add_argument(f"--{name}", type=Q)
    parser.add_argument("--radius", type=int)
    parser.add_argument("--degree-cap", type=int)
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        values = (args.c, args.eta, args.h_lo, args.h_hi, args.radius, args.degree_cap)
        if None in values:
            parser.error("all six certified inputs are required")
        rows = envelopes(*values)
        s = cutoff_power(rows)
        assert certifies(s, rows)
        print(f"N0 = 2^{s}; p_(n+1) > p_(n+2) for n >= N0")
