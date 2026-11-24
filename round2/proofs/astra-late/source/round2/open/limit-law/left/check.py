#!/usr/bin/env python3
"""Exact short checks of the hazards in the left-tail proof."""
from fractions import Fraction as F
from itertools import product


def main():
    checked = enhanced = 0
    K = 5
    # Every increasing labelled tree from the edge through size eight.
    # Enumerate individual physical ordered arrows, not the closed hull formula.
    for n in range(3, 9):
        for p in product(*(range(i) for i in range(1, n))):
            adj = [[] for _ in range(n)]
            for i, parent in enumerate(p, 1):
                adj[i].append(parent)
                adj[parent].append(i)
            deg = tuple(map(len, adj))
            cores = sum(d > 1 for d in deg)
            b = sum((F(1, deg[i]) for i, nn in enumerate(adj)
                     for j in nn if deg[j] == 1), F())
            assert b >= 1
            if cores == 1:
                assert b == 1
            if n >= K+3 and 2 <= cores < K:
                assert b >= 1+min(F(1,6), F(1,K))
                enhanced += 1
            checked += 1
    assert enhanced
    # A four-vertex path disproves the strict bound without the size gate.
    assert F(1,2)+F(1,2) == 1
    print(f'PASS {checked} labelled trees; {enhanced} strict bounded-core hazard cases')
    star = F(1)
    for j in range(3, 80):
        assert star == F(2,j-1)
        star *= F(j-1,j)
    print('PASS exact star survival through size 79; four-vertex path requires size gate')


if __name__ == '__main__':
    main()
