#!/usr/bin/env python3
"""Exact counts of dense, symmetry-packed, pure-inner rooted Level-5 tables."""
from collections import Counter
from functools import cache
from itertools import combinations_with_replacement
from math import comb


@cache
def shapes(n):
    if n == 1:
        return ((),)
    candidates = [(s, x) for s in range(1, n) for x in shapes(s)]
    return tuple(tuple(x for _, x in xs) for k in range(1, n)
                 for xs in combinations_with_replacement(candidates, k)
                 if sum(s for s, _ in xs) == n - 1)


def label(shape):
    return '(' + ''.join(map(label, shape)) + ')'


@cache
def rows(shape, d5, h):
    # A vertex must have at least as many core children as the motif includes.
    own = sum(max(0, d - h - len(shape) + 1) for d in range(2, d5))
    for child, mult in Counter(shape).items():
        own *= comb(rows(child, d5, 1) + mult - 1, mult)
    return own


for d5 in (6, 7, 8):
    print('L4 candidate all-inner rows before L4 cap filtering:',
          [(label(s), rows(s, d5, 0)+rows(s, d5, 1)) for s in shapes(4)])
    total = 0
    for shape in shapes(5):
        root, nonroot = (rows(shape, d5, h) for h in (0, 1))
        total += root + nonroot
        print(d5, label(shape), root, nonroot, root + nonroot)
    print(f'TOTAL D5={d5}: entries={total} bytes={8*total} double_buffer_bytes={16*total}')
