#!/usr/bin/env python3
"""Independent exact controls for the finite mechanisms in the density proof."""
from fractions import Fraction as F
from itertools import product


def add(adj, v):
    adj = [set(row) for row in adj] + [{v}]
    adj[v].add(len(adj) - 1)
    return adj


def side(adj, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for w in adj[stack.pop()] - seen - {blocked}:
            seen.add(w)
            stack.append(w)
    return seen


def check_cut(adj, x, y):
    left, right = side(adj, x, y), side(adj, y, x)
    assert left.isdisjoint(right) and left | right == set(range(len(adj)))
    cross = F(1, len(adj[x])) + F(1, len(adj[y]))
    total = F(0)
    for vertices, root in ((left, x), (right, y)):
        assert len(vertices) >= 2
        rate = F(0)
        for v in vertices:
            degree = len(adj[v] & vertices) + (v == root)
            assert degree == len(adj[v])
            for w in adj[v] & vertices:
                rate += F(1, degree)
                if v == root and len(adj[w]) == 1:
                    assert F(1, 2) * F(2, degree) == F(1, degree)
        assert rate == len(vertices) - F(1, len(adj[root]))
        total += rate
    assert total + cross == len(adj)
    return total, cross


def main():
    trees = cuts = 0
    for n in range(3, 8):
        for parents in product(*(range(i) for i in range(1, n))):
            adj = [set()]
            for p in parents:
                adj = add(adj, p)
            x = next(v for v in range(n) if len(adj[v]) >= 2
                     and sum(len(adj[w]) >= 2 for w in adj[v]) <= 1)
            grown = add(add(adj, x), x)
            leaves = [w for w in grown[x] if len(grown[w]) == 1]
            assert F(len(leaves), len(grown[x])) >= F(3, 4)
            for y in leaves:
                joined = add(grown, y)
                assert all(any(len(joined[w]) == 1 for w in joined[v])
                           for v in (x, y))
                check_cut(joined, x, y)
                cuts += 1
            trees += 1
    print(f'PASS: {trees} labelled seeds, {cuts} promoted-leaf cuts and exact arrow rates')

    # P4 split at the middle: inert internal rates total 3. Two possible
    # first internal arrows have equal prior weight. Sending from the left
    # leaf grows its root degree to 3; sending from the root grows its leaf.
    # The future probability of no cross over one fixed silent-internal
    # interval is exp(-5h/6) versus exp(-h), a path-dependent tilt.
    path = [{1}, {0, 2}, {1, 3}, {2}]
    internal, cross = check_cut(path, 1, 2)
    hazards = [F(1, len(add(path, w)[1])) + F(1, 2) for w in (1, 0)]
    assert hazards == [F(5, 6), F(1)]
    assert internal == 3 and cross / (internal + cross) == F(1, 4)
    print('COUNTEREXAMPLE: equal-prior path prefixes have cross hazards 5/6 and 1')

    # At mu=2/3 set u=r^(1/3). F(r)^3 is rational at rational u,
    # and a turning point lies in (0,1) exactly when B*C < A*D.
    for a, b, c, d in ((2, 1, 1, 1), (1, 2, 1, 1), (1, 1, 1, 1)):
        def cube(u):
            return (b*u*u+d)**3 / (a*u**3+c)**2

        grid = [F(k, 16) for k in range(17)]
        critical = F(b*c, a*d)
        for lo, hi in zip(grid, grid[1:]):
            if hi <= critical:
                assert cube(hi) > cube(lo)
            elif lo >= critical:
                assert cube(hi) < cube(lo)
    print('PASS: exact two-branch, monotone, and endpoint-critical smoothing fixtures')


if __name__ == '__main__':
    main()
