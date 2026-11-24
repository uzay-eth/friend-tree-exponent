#!/usr/bin/env python3
"""Exact finite checks. No third-party dependencies; tuples are unordered trees."""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import lcm


def encode(par):
    children = [[] for _ in par]
    for v in range(len(par) - 1, 0, -1):
        children[par[v]].append(tuple(sorted(children[v])))
    return tuple(sorted(children[0]))


def decode(tree):
    par = []

    def visit(branch, parent):
        v = len(par)
        par.append(parent)
        for child in branch:
            visit(child, v)

    visit(tree, -1)
    return tuple(par)


def levels(max_size):
    out = {1: {()}}
    for n in range(2, max_size + 1):
        out[n] = {encode(decode(t) + (v,)) for t in out[n - 1] for v in range(n - 1)}
    return out


def adjacency(par):
    adj = [[] for _ in par]
    for v in range(1, len(par)):
        adj[v].append(par[v])
        adj[par[v]].append(v)
    return adj


def weight_check(trees):
    # Label i means mark 1/i. Use a common integer denominator for z and
    # for 1/(a_root+d), and use h in hundredths. No float decisions.
    labels = range(2, 21)
    d = Q(3, 20)
    z = {i: Q(1, i) / (Q(1, i) + d) for i in labels}
    r = {i: 1 / (Q(1, i) + d) for i in labels}
    zden = lcm(*(x.denominator for x in z.values()))
    rden = lcm(*(x.denominator for x in r.values()))
    zi = {i: int(x * zden) for i, x in z.items()}
    ri = tuple(int(x * rden) for x in r.values())
    # Since root label contributes only the normalization, collapse all 19
    # root choices to the exact worst case; count each covered assignment.
    min_root = min(ri)
    checked = 0
    for n in range(2, 7):
        for t in trees[n]:
            par = decode(t)
            adj = adjacency(par)
            D = [len(a) + (v == 0) for v, a in enumerate(adj)]
            cores = [v for v in range(n) if D[v] >= 2]
            core_set = set(cores)
            children = [v for v in cores if par[v] == 0]
            # For each N, max_{|f|<=1}|phi| is attained by f=1.
            phi = max((sum(D[v] <= N for v in children) for N in range(2, D[0])), default=0)
            ch = [min(4, sum(w in core_set for w in adj[v]) + (v == 0)) for v in cores]
            h = tuple((200, 100, 85, 78, 70)[c] for c in ch)
            parents = [cores.index(par[v]) for v in cores[1:]]
            den = 100 * rden * zden ** (len(cores) - 1)
            for marks in product(labels, repeat=len(cores) - 1):
                weights = [zden ** (len(cores) - 1)]
                for parent, mark in zip(parents, marks):
                    weights.append(weights[parent] * zi[mark] // zden)
                num = min_root * sum(a * b for a, b in zip(weights, h))
                if 7 * den * phi > 26 * num:
                    raise AssertionError((t, marks, phi, Q(num, den)))
                checked += len(ri)
    # Reachable capped star-of-core-children, plus two deliberately false
    # variants: removing the cap and replacing 26/7 by 1.
    # Root label 2; leaves can be gained first, then each child promoted.
    def star_weight(k, birth_parent_degree, capped):
        a = Q(1, min(birth_parent_degree, 20) if capped else birth_parent_degree)
        return (Q(7, 10) + k * a / (a + d)) / (Q(1, 2) + d)

    k, deg = 99, 100
    assert k <= Q(26, 7) * star_weight(k, deg, True)
    assert k > Q(26, 7) * star_weight(k, deg, False)
    assert k > star_weight(k, deg, True)
    print(f'W: {checked:,} marked cases (including extremal root reduction); 2 planted negative controls passed')


def intrinsic_rates(tree):
    par = decode(tree)
    adj = adjacency(par)
    rates = defaultdict(Q)
    for v, neighbors in enumerate(adj):
        d = len(neighbors) + (v == 0)
        for w in neighbors:
            rates[encode(par + (w,))] += Q(1, d)
    return rates, Q(1, len(adj[0]) + 1)


def holding_check(trees):
    # Marks do not enter intrinsic rates, so lump them out exactly. For L^0
    # killed at size 9, solve the row system x(-L^0)=e_J by forward
    # substitution (the matrix is triangular in vertex count). x is the J
    # row of the resolvent (-L^0)^(-1), i.e. expected occupation times.
    inflow = defaultdict(Q, {((),): Q(1)})
    print('Holding: v, reachable states, exact mass (rounded in output), upper bound')
    for v in range(2, 9):
        assert sum(inflow[t] for t in trees[v]) == 1, ('entrance mass', v)
        occupation = Q(0)
        reachable = 0
        for t in trees[v]:
            rates, null = intrinsic_rates(t)
            rate = sum(rates.values())
            assert rate == v - null >= Q(2 * v - 1, 2)
            x = inflow[t] / rate
            assert x * rate == inflow[t]  # Exact diagonal residual of the solve.
            reachable += bool(x)
            occupation += x
            for nxt, q in rates.items():
                inflow[nxt] += x * q
        assert occupation <= Q(2, 2 * v - 1)
        assert (occupation == Q(2, 2 * v - 1)) == (v == 2)
        print(f'  {v}: {reachable:3d}  {float(occupation):.10f} <= {Q(2, 2*v-1)}')


def components(par, N):
    # Preserve vertex identities to check the stronger partition invariance.
    D = list(map(len, adjacency(par)))
    groups = defaultdict(list)
    top = {}
    for v, p in enumerate(par):
        if D[v] <= N:
            top[v] = top.get(p, v)
            groups[top[v]].append(v)
    return {root: tuple(vs) for root, vs in groups.items() if D[root] >= 2}


def attachment_check(trees):
    checked = 0
    for n in range(3, 9):
        for t in trees[n]:
            par = decode(t)
            if par.count(0) < 2:
                continue  # Physical founder is a core.
            for x, adj in enumerate(adjacency(par)):
                for N in range(2, len(adj)):
                    assert components(par, N) == components(par + (x,), N), (par, x, N)
                    checked += 1
    # x=1 crosses the cutoff, and the surviving clean founder component
    # loses its dirty descendant and the descendant's original leaf.
    before = (-1, 0, 0, 1)
    after = before + (1,)
    assert components(before, 2) == {0: (0, 1, 2, 3)}
    assert components(after, 2) == {0: (0, 2)}
    # Isolated *original* cores can survive without a real leaf or child.
    p = (-1, 0, 1, 0, 0, 2, 2, 0)
    assert components(p, 2) == {1: (1,)}
    print(f'Hub insertion: {checked:,} invariant cases through size 8; crossing and singleton controls passed')


if __name__ == '__main__':
    trees = levels(8)
    assert [len(trees[n]) for n in range(1, 9)] == [1, 1, 2, 4, 9, 20, 48, 115]
    weight_check(trees)
    holding_check(trees)
    attachment_check(trees)
