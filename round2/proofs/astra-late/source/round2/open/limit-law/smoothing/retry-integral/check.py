#!/usr/bin/env python3
"""Exact arithmetic controls for the retry integral (no simulation)."""
from fractions import Fraction as Q


def clock():
    r = Q(747, 10000)
    weights = [Q(j + 1, 2 ** (j + 2)) for j in range(1, 9)]
    tail_at_zero = Q(1, 2) + Q(1, 2) ** 2
    remaining = tail_at_zero - sum(weights)
    assert remaining == Q(11, 1024)
    upper = tail_at_zero + r * (
        sum(w / (Q(j, 2) - r) for j, w in enumerate(weights, 1))
        + remaining / (Q(9, 2) - r)
    )
    lower = sum(Q(j + 1, 2 ** (j + 2)) * j / (j - 2 * r) for j in range(1, 81))
    assert lower < upper < 1
    rho = (1 + upper * Q(2, 3) / (Q(2, 3) - r)) / (2 - r)
    assert rho < Q(99991, 100000)
    print(f"PASS exact defective-clock envelope: rho <= {rho} < 99991/100000")

    sigma = r / (1 + 2 * r)
    assert sigma == Q(1, 2) - Q(1, 2) / (1 + 2 * r) == Q(747, 11494)
    assert Q(13, 278) < Q(13, 226) < sigma
    assert sigma - Q(13, 278) == Q(14561, 798833)
    assert 1 / (1 - sigma) == Q(11494, 10747)
    print("PASS full-law cutoff 747/11494; increase 14561/798833")


def reciprocal():
    for k in range(1, 17):
        q = Q(1, k)
        v = 1 - q
        # Harmonic block: first k terms, then dyadic blocks, then geometric tail.
        # Each block (2^i k, 2^(i+1) k] contributes at most q v^(2^i k).
        cutoff = 8 * k
        partial = sum(q * v ** (j - 1) / (j + 1) for j in range(1, cutoff + 1))
        tail = v ** cutoff / (cutoff + 2)
        harmonic = sum(Q(1, j) for j in range(1, k + 1))
        envelope = q * (harmonic + sum(v ** (2 ** i * k) for i in range(3))) + tail
        assert partial + tail <= envelope <= q * (harmonic + 4)
    print("PASS 16 rational geometric-reciprocal envelopes")


def fold():
    # At x=0, density 1/(2 sqrt(m)) under m=y^2 is uniform in y.
    # On (4^(-k-1), 4^(-k)] the kernel 1/sqrt(m) is at least 2^k.
    for k in range(32):
        y_hi, y_lo = Q(1, 2 ** k), Q(1, 2 ** (k + 1))
        assert (y_hi - y_lo) / y_hi == Q(1, 2)
    print("PASS 32 critical fold shells: lower contribution 1/2 per shell")


if __name__ == "__main__":
    assert __debug__
    clock()
    reciprocal()
    fold()
