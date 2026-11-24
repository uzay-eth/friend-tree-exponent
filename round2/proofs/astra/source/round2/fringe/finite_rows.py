#!/usr/bin/env python3
"""Rational enclosures for the 19-type renewal matrix at finite fringe size."""

import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import cache
from math import lcm

from check import RESIDUAL_BOUND, degree, frac_mark, gain, intrinsic_transitions, newborn, node, weight

LABELS = tuple(range(2, 21))
SCALE = 10**12


@cache
def cap_marks(state):
    leaves, mark, children = state
    return node(leaves, min(mark, 20), (cap_marks(child) for child in children))


@cache
def size(state):
    return 1 + state[0] + sum(size(y) for y in state[2])


def step(state):
    inv = F(1, degree(state, 1))
    row = Counter()
    for y, rate in intrinsic_transitions(state).items():
        row[cap_marks(y)] += rate
    for y in state[2]:
        row[gain(y)] += inv
        if frac_mark(y[1]) > inv:
            row[y] += frac_mark(y[1]) - inv
    return row


def ceil_fraction(x):
    return -(-x.numerator // x.denominator)


def prepare(cutoff, discount):
    initial = [newborn(i) for i in LABELS]
    queue = deque(initial)
    index = {s: i for i, s in enumerate(initial)}
    rows = []
    for_expected_starts = list(index)
    assert for_expected_starts == initial

    while queue:
        x = queue.popleft()
        assert index[x] == len(rows)
        den = discount + frac_mark(x[1]) + size(x) - F(1, degree(x, 1))
        direct = F(x[0], degree(x, 1)) / den
        successors = []
        boundary = F(0)
        for y, rate in step(x).items():
            coefficient = rate / den
            if size(y) > cutoff:
                assert size(y) == cutoff + 1
                boundary += coefficient * weight(y) / (discount - RESIDUAL_BOUND)
            else:
                if y not in index:
                    index[y] = len(index)
                    queue.append(y)
                successors.append((index[y], coefficient))

        # The resulting Bellman bounds use exclusively integers. One common
        # denominator per source also covers the explicit exit value.
        denominator = lcm(direct.denominator, boundary.denominator, *(c.denominator for _, c in successors))
        arcs = [(j, int(c * denominator)) for j, c in successors]
        source = int(SCALE * direct * denominator)
        exit_reward = int(SCALE * boundary * denominator)
        rows.append((arcs, min(degree(x, 1), 20) - 2, source, exit_reward, denominator))

    initial_upper = [ceil_fraction(SCALE * weight(x) / (discount - RESIDUAL_BOUND)) for x in index]
    return rows, initial_upper


def kernel(rows, bounds, iterations):
    # Near a small cutoff only the first labels have nonzero direct source.
    # The upper continuation is common to every output. Keeping it as one
    # extra column saves the all-zero lower columns.
    visible = 1 + max(label for _, label, _, _, _ in rows)
    width = visible + 1
    low = [[0] * width for _ in rows]
    high = [[bound] * width for bound in bounds]
    for _ in range(iterations):
        new_low, new_high = [], []
        for pos, (arcs, label, reward, boundary, den) in enumerate(rows):
            lower = [0] * width
            upper = [boundary] * width
            lower[label] = reward
            upper[label] += reward
            for j, c in arcs:
                for k in range(width):
                    lower[k] += c * low[j][k]
                    upper[k] += c * high[j][k]
            new_low.append([max(low[pos][k], lower[k] // den) for k in range(width)])
            new_high.append([min(high[pos][k], -(-upper[k] // den)) for k in range(width)])
        low, high = new_low, new_high
    lo_matrix = [row[:visible] + [0] * (len(LABELS) - visible) for row in low[:len(LABELS)]]
    hi_matrix = [row[:visible] + [row[-1]] * (len(LABELS) - visible) for row in high[:len(LABELS)]]
    assert all(0 <= a <= b for lo, hi in zip(lo_matrix, hi_matrix) for a, b in zip(lo, hi))
    return lo_matrix, hi_matrix


def perron_enclosure(matrix):
    # Floating arithmetic proposes a test column. The returned interval is
    # obtained by reapplying the integer matrix to a positive integer column.
    g = [1.0] * len(matrix)
    for _ in range(64):
        v = [sum(x * y for x, y in zip(row, g)) for row in matrix]
        high = max(v)
        assert high > 0
        g = [x / high for x in v]
    trial = [max(1, round(x * SCALE)) for x in g]
    y = [sum(x * w for x, w in zip(row, trial)) for row in matrix]
    ratios = [F(a, SCALE * b) for a, b in zip(y, trial)]
    return min(ratios), max(ratios)


def endpoints(lo, hi):
    unit = 10**8
    return f"[{(lo * unit).__floor__()}/{unit}, {ceil_fraction(hi * unit)}/{unit}]"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoffs", nargs="+", type=int, default=[2, 3, 4, 5])
    parser.add_argument("--discount", type=F, default=F(113, 200))
    parser.add_argument("--steps", type=int, default=96)
    args = parser.parse_args()
    assert args.discount > RESIDUAL_BOUND and args.steps > 0
    print(f"discount = {args.discount}, integer unit = 1/{SCALE}, Bellman iterations = {args.steps}")
    print("Intervals below enclose r(K(s)) at the stated discount.")
    for n in args.cutoffs:
        if n < 2:
            raise ValueError("cutoff must be at least two")
        rows, bound = prepare(n, args.discount)
        low, high = kernel(rows, bound, args.steps)
        lo = perron_enclosure(low)[0]
        hi = perron_enclosure(high)[1]
        print(f"N={n}: {len(rows)} marked fringes, {endpoints(lo, hi)}")
        # At N=2 the lower matrix has just one nonzero column. The only
        # occupied state at a lower successful birth has remained at J_i.
        if n == 2:
            target = F(1, 2) / (args.discount + 2)
            assert lo <= target and low[0][0] <= SCALE * target < low[0][0] + 1


if __name__ == "__main__":
    main()
