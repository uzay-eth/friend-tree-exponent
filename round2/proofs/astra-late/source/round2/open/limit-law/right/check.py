#!/usr/bin/env python3
"""Exact finite controls for the embedded gluing comparison."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product


def urn(k, m):
    steps = k * (m - 2)
    table = defaultdict(Q)
    for draws in product(range(k), repeat=steps):
        counts = [2] * k
        w = Q(1)
        for i, j in enumerate(draws):
            w *= Q(counts[j], 5*k+i)
            counts[j] += 1
        table[tuple(counts)] += w
    no_scaffold = Q(1)
    for i in range(steps):
        no_scaffold *= Q(2*k+i, 5*k+i)
    assert sum(table.values()) == no_scaffold
    balanced = table[(m,) * k]
    assert balanced == max(table.values())
    for counts, w in table.items():
        ratio = Q(1)
        for n in counts:
            ratio *= Q(n - 1, m - 1)
        assert w == balanced * ratio
    return len(table)


def main():
    rows = sum(urn(k, m) for k in (2, 3) for m in (3, 4, 5))
    print(f'PASS: {rows} exact urn compositions; scaffold avoidance and balance')

    # Two edge blocks joined to a scaffold path of size six. In one block
    # the root sends to the original leaf, then the original leaf sends to
    # the root, then the root sends again at free degree two.
    sizes = (10, 11, 12)
    d_free = (1, 2, 2)
    d_physical = (2, 2, 3)
    physical = reference = Q(1)
    for n, d, degree in zip(sizes, d_free, d_physical):
        physical *= Q(1, n * degree)
        reference *= Q(1, n * d)
    assert physical / reference == Q(1, 3)

    channels = 0
    for d in range(1, 12):
        for c in range(1, 8):
            deficit = sum(Q(1, d) - Q(1, d+c) for _ in range(d))
            extra = sum(Q(1, d+c) for _ in range(c))
            assert deficit == extra == Q(c, d+c)
            channels += 1
    print(f'PASS: exact directed-arrow likelihood and {channels} cut-source generator rows')

    # From T3 the first birth is a promotion with probability 1/3, and
    # the waiting time is Exp(3), independent of that winner. At t=log 2
    # its survival is 1/8. Conditioning on departure at size 4 by 2t
    # truncates the clock even though the same size alone does not.
    survival = Q(1, 8)
    single_step_departure = Q(1, 3)
    early = single_step_departure * (1 - survival)
    by_cutoff = single_step_departure * (1 - survival**2)
    assert early / by_cutoff == Q(8, 9) != Q(7, 8)
    print('COUNTEREXAMPLE: stopping by physical time changes the clock probability from 7/8 to 8/9')


if __name__ == '__main__':
    main()
