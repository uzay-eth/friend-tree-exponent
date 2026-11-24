#!/usr/bin/env python3
"""Exact signed coefficients of future laws on rooted unlabelled trees."""
from collections import defaultdict
from fractions import Fraction as F
from math import lcm


def rooted(adj, x, parent=-1):
    return '(' + ''.join(sorted(rooted(adj, y, x) for y in adj[x] if y != parent)) + ')'


def labelled(N):
    """Increasing-parent histories, with common exact denominator at each size."""
    adj = [[1], [0]] + [[] for _ in range(N - 2)]
    scales = {n: lcm(*range(1, n)) for n in range(2, N)}

    def visit(n, w, den):
        yield n, w, den, adj
        if n == N:
            return
        s = scales[n]
        # Tally directed arrows by target, independently of the rooted DP.
        rates = [0] * n
        for x in range(n):
            for y in adj[x]:
                rates[y] += s // len(adj[x])
        for x, a in enumerate(rates):
            adj[x].append(n)
            adj[n].append(x)
            yield from visit(n + 1, w * a, den * s * n)
            adj[n].pop()
            adj[x].pop()

    yield from visit(2, 1, 1)


def decode(key):
    adj, stack = [], []
    for c in key:
        if c == '(':
            i = len(adj)
            adj.append([])
            if stack:
                p = stack[-1]
                adj[i].append(p)
                adj[p].append(i)
            stack.append(i)
        else:
            stack.pop()
    return adj


def steps(key, s):
    adj = decode(key)
    out = defaultdict(int)
    for x, nbrs in enumerate(adj):
        a = sum(s // len(adj[y]) for y in nbrs)
        nbrs.append(len(adj))
        adj.append([x])
        out[rooted(adj, 0)] += a
        adj.pop()
        nbrs.pop()
    assert sum(out.values()) == len(adj) * s
    return out


def advance(weights, n):
    s = lcm(*range(1, n))
    out = defaultdict(int)
    for key, w in weights.items():
        for k, a in steps(key, s).items():
            out[k] += w * a
    return {k: w for k, w in out.items() if w}, n * s


def expected():
    # Hand-derived coefficients; see REPORT.md.
    return {
        3: {'(()())': F(1, 2), '((()))': -F(1, 2)},
        4: {'((())())': F(1, 6), '(((())))': -F(1, 6)},
        5: dict(zip(('((((()))))', '(((()())))', '(((())()))', '(((()))())', '((()())())'),
                    (F(-1,24), F(-1,24), F(-1,72), F(1,24), F(1,18)))),
        6: dict(zip(('(((((())))))', '((((()()))))', '((((())())))', '((((())))())', '(((()()())))',
                     '(((()()))())', '(((())()()))', '(((())())())', '((()()())())'),
                    (F(-1,120), F(-1,60), F(-1,60), F(1,120), F(-1,60), F(1,60), F(-1,120), F(1,60), F(1,40)))),
    }


def check(N=9):
    nums = defaultdict(lambda: defaultdict(int))
    denoms = {}
    for n, w, den, adj in labelled(N):
        denoms[n] = den
        for v in range(2, min(n, 6) + 1):
            p, q = rooted(adj, v - 2), rooted(adj, v - 1)
            if p != q:
                nums[n, v][p] += w
                nums[n, v][q] -= w
    actual = {}
    for v in range(2, 7):
        w = {k: x for k, x in nums[v, v].items() if x}
        d = denoms[v]
        actual[v] = {k: F(x, d) for k, x in w.items()}
        assert sum(w.values()) == 0
        print(f'{v-1} vs {v}: {dict(sorted(actual[v].items()))}')
        for n in range(v, N):
            w, factor = advance(w, n)
            d *= factor
            assert d == denoms[n + 1]
            assert w == {k: x for k, x in nums[n + 1, v].items() if x}
    for v, law in expected().items():
        assert actual[v] == law
        assert actual[v] != {k: -x for k, x in law.items()}
    assert not actual[2]
    print(f'Labelled vs rooted propagation agrees through n={N}; reversed-sign controls rejected.')
    return actual


if __name__ == '__main__':
    check()
