#!/usr/bin/env python3
"""Exact independent arrow and marked-Poisson controls; no sample paths."""
from collections import Counter, deque
from fractions import Fraction as F
from math import comb, factorial


def A(d):
    return F(d * factorial(d + 1), 12 * 4**(d - 2))


def V(d, l):
    return A(d) * F(2, d + 1)**l


def labelled_trees(planted, cap):
    first = (0, 1) if planted else (1, 1)
    queue = deque([first])  # parent of node i+2; physical founder is 1
    while queue:
        tree = queue.popleft()
        yield tree
        if len(tree) < cap:
            queue.extend(tree + (y,) for y in range(1, len(tree) + 2))


def physical(planted, cap=6, defect=None):
    rows = mixed = 0
    for tree in labelled_trees(planted, cap):
        edges = {i: [] for i in range(len(tree) + 2)}
        for y, x in enumerate(tree, 2):
            edges[y].append(x)
            edges[x].append(y)
        if planted:
            edges[1].append(0)
        degree = {x: len(nb) for x, nb in edges.items()}
        leaves = {x: sum(y != 0 and degree[y] == 1 for y in nb)
                  for x, nb in edges.items()}
        assert degree[1] >= 2
        for v in range(1, len(tree) + 2):
            d, l = degree[v], leaves[v]
            if d < 2:
                continue
            directed = Counter()
            for x in range(1, len(tree) + 2):
                for y in edges[x]:
                    if y == 0:
                        continue
                    # Appending a leaf at y changes v's degree iff y=v;
                    # otherwise v's leaf count only notices an adjacent leaf.
                    dd = d + (y == v)
                    ll = l + (y == v) - (y in edges[v] and degree[y] == 1)
                    if (dd, ll) != (d, l):
                        directed[dd, ll] += F(1, degree[x])
            incoming = l + sum(F(1, degree[y]) for y in edges[v]
                               if y != 0 and degree[y] >= 2)
            p = F(l, d + (defect == 'promotion-denominator'))
            predicted = Counter()
            if incoming:
                predicted[d + 1, l + 1] += incoming
            if p:
                predicted[d + (defect == 'promotion-degree'), l - 1] += p
            assert directed == predicted, (tree, v, directed, predicted)
            assert sum(rate * (V(*nxt) - V(d, l))
                       for nxt, rate in directed.items()) <= 0
            if l:
                assert V(d + 1, l + 1) / V(d, l) == (
                    F(d + 1, 2*d) * F(d + 1, d + 2)**l)
                assert V(d, l - 1) / V(d, l) == F(d + 1, 2)
            rows += 1
            mixed += bool(incoming > l and p)
    return rows, mixed


def pgf_control(defect=None):
    # Divide a rate-one clock on [0,log 4] into equal time halves.
    # Marks 0,1,2 have probabilities 1/4,1/2,1/4. At thresholds 1/2
    # and 1 the selected time lengths are exact multiples of log 2.
    # Only intensity integrated over time is used; endpoint atoms have zero time.
    # Categories: >=1 from mark 2, [1/2,1) from mark 1 early or mark 2 late,
    # then the rest. Conditional on n points the counts must be multinomial.
    weights = ((0, F(1,4)), (1, F(1,2)), (2, F(1,4)))
    counts = [F(0)]*3
    for discount in (F(3,4), F(3,8)):  # interior representative of each half
        for mark, mass in weights:
            x = mark*discount
            counts[0 if x >= 1 else 1 if x >= F(1,2) else 2] += mass/2
    cells = tuple(counts)
    if defect == 'condition-on-mark':
        cells = (F(1, 6), F(1, 2), F(1, 3))
    if defect == 'overlap-categories':
        cells = (F(1, 8), F(1, 2), F(1, 2))
    n, x, y = 5, F(2, 3), F(3, 5)
    actual = sum(comb(n, i)*comb(n-i, j)*cells[0]**i*cells[1]**j
                 * cells[2]**(n-i-j)*x**i*y**j
                 for i in range(n+1) for j in range(n-i+1))
    assert actual == (1 + (x-1)/8 + 3*(y-1)/8)**n
    # Coefficients of the exponential generating series for the unnormalised
    # Poisson pgf factor, separately check the independent exponential factors.
    for n in range(9):
        factored = sum((F(1, 8)*x)**i * (F(3, 8)*y)**j * F(1, 2)**(n-i-j)
                       / (factorial(i)*factorial(j)*factorial(n-i-j))
                       for i in range(n+1) for j in range(n-i+1))
        assert factored == ((x+3*y)/8+F(1, 2))**n / factorial(n)


def reject(name, test):
    try:
        test()
    except AssertionError:
        print(f'rejected {name}')
    else:
        raise RuntimeError(f'undetected planted fault: {name}')


if __name__ == '__main__':
    for planted in (False, True):
        rows, mixed = physical(planted)
        print(f'planted={planted}: {rows} physical core rows; {mixed} mixed rows')
    for d in range(2, 65):
        assert A(d+1) / A(d) == F((d+1)*(d+2), 4*d) > 1
        star = F(1)
        for j in range(2, d):
            star *= F(j, j+1)
        assert star == F(2, d)
        contest = F(1)
        for l in range(1, d):
            contest *= F(l, d) / (l + F(l, d))
        assert contest == F(1, d+1)**(d-1)
    pgf_control()
    print('factorial ratios and silent-star contests through degree 64; Poisson allocation verified')
    for fault in ('promotion-denominator', 'promotion-degree'):
        reject(fault, lambda: physical(True, defect=fault))
    for fault in ('condition-on-mark', 'overlap-categories'):
        reject(fault, lambda: pgf_control(fault))
