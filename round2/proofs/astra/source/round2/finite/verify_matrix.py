"""Independently check the published matrices, rational residuals and polynomials."""

import json
import sys
from collections import defaultdict
from fractions import Fraction as F


def unpack(code):
    # Store leaf and inactive-slot counts, with all full degrees reconstructed.
    graph, leaf, slot, degree = {}, {}, {}, {}

    def read(item, parent):
        total, count, children = item
        i = len(graph)
        graph[i] = set()
        leaf[i], degree[i] = count, total
        if parent is not None:
            graph[parent].add(i)
            graph[i].add(parent)
        for child in children:
            read(child, i)
    read(code, None)
    for i in graph:
        slot[i] = degree[i] - leaf[i] - len(graph[i])
        assert slot[i] >= 0
    return graph, leaf, slot


def connected(graph):
    todo = set(graph)
    parts = []
    while todo:
        piece, stack = set(), [next(iter(todo))]
        while stack:
            i = stack.pop()
            if i in todo:
                todo.remove(i)
                piece.add(i)
                stack.extend(graph[i])
        parts.append(piece)
    return parts


def signature(graph, leaf, slot, vertices):
    # Strip peripheral layers to the one- or two-vertex center. This does not
    # share the builder's trial-root minimization or its (D,L) labels.
    remaining = set(vertices)
    while len(remaining) > 2:
        tips = {u for u in remaining if len(graph[u] & remaining) <= 1}
        assert tips and tips != remaining
        remaining -= tips

    def half(i, previous):
        return (leaf[i], slot[i], tuple(sorted(half(j, i) for j in graph[i] & vertices
                                                if j != previous)))
    center = sorted(remaining)
    if len(center) == 1:
        return ('single', half(center[0], None))
    u, v = center
    return ('edge', tuple(sorted((half(u, v), half(v, u)))))


def rebuild_rows(case):
    types = case['types']
    cutoff, cap = case['M'], case['B']
    hub = types.index('hub')
    index = {}
    for i, code in enumerate(types):
        if i == hub:
            continue
        graph, leaf, slot = unpack(code)
        assert 0 < len(graph) <= cap
        assert all(2 <= leaf[u]+slot[u]+len(ns) < cutoff for u, ns in graph.items())
        assert len(connected(graph)) == 1
        index[signature(graph, leaf, slot, set(graph))] = i
    assert len(index) + 1 == len(types)

    rows = []
    for i, code in enumerate(types):
        row = defaultdict(F)
        if i == hub:
            row[index[('single', (1, 1, ()))]] = F(1)
            rows.append(dict(row))
            continue
        graph, leaf, slot = unpack(code)
        degree = {u: leaf[u]+slot[u]+len(ns) for u, ns in graph.items()}
        choices = []
        for u, neighbors in graph.items():
            # Enumerate individual sources and neighbor choices. Ghost targets
            # of a low core are null rings, and do not enter the generator.
            for _ in range(leaf[u]):
                choices.append(('gain', u, F(1)))  # ring of a leaf
                choices.append(('bud', u, F(1, degree[u])))  # chosen by its parent
            for neighbor in neighbors:
                choices.append(('gain', neighbor, F(1, degree[u])))
        for event, u, rate in choices:
            row[i] -= rate
            g = {x: ns.copy() for x, ns in graph.items()}
            stock, ports = leaf.copy(), slot.copy()
            if event == 'gain':
                stock[u] += 1
                if degree[u]+1 == cutoff:
                    row[hub] += rate
                    for neighbor in g.pop(u):
                        g[neighbor].remove(u)
                        ports[neighbor] += 1
            else:
                child = max(g)+1
                stock[u] -= 1
                stock[child], ports[child] = 1, 0
                g[child] = {u}
                g[u].add(child)
            if len(g) <= cap:
                for part in connected(g):
                    row[index[signature(g, stock, ports, part)]] += rate
                continue

            # Independently test each undirected edge by explicitly deleting it,
            # rebuilding both slots and computing the resulting components.
            options = []
            for u1 in g:
                for v1 in g[u1]:
                    if u1 >= v1:
                        continue
                    cut = {x: ns.copy() for x, ns in g.items()}
                    cut[u1].remove(v1)
                    cut[v1].remove(u1)
                    pieces = connected(cut)
                    assert len(pieces) == 2
                    if min(map(len, pieces)) * (cutoff-1) < cap:
                        continue
                    ghost = ports.copy()
                    ghost[u1] += 1
                    ghost[v1] += 1
                    options.append([index[signature(cut, stock, ghost, part)] for part in pieces])
            assert options
            for outcome in options:
                for target in outcome:
                    row[target] += rate/len(options)
        rows.append({j: x for j, x in row.items() if x})
    return rows, index[('single', (2, 0, ()))]


def reachable(rows, start):
    visited, todo = set(), [start]
    while todo:
        i = todo.pop()
        if i not in visited:
            visited.add(i)
            todo.extend(j for j, q in rows[i].items() if j != i and q > 0)
    return visited


def det_bareiss(matrix):
    n = len(matrix)
    previous, sign = 1, 1
    for k in range(n-1):
        pivotrow = next((i for i in range(k, n) if matrix[i][k]), None)
        if pivotrow is None:
            return 0
        if pivotrow != k:
            matrix[k], matrix[pivotrow] = matrix[pivotrow], matrix[k]
            sign *= -1
        pivot = matrix[k][k]
        for i in range(k+1, n):
            factor = matrix[i][k]
            for j in range(k+1, n):
                value = pivot*matrix[i][j] - factor*matrix[k][j]
                quotient, remainder = divmod(value, previous)
                assert not remainder
                matrix[i][j] = quotient
            matrix[i][k] = 0
        previous = pivot
    return sign*matrix[-1][-1]


def value(coefficients, x):
    v = 0
    for c in coefficients:
        v = v*x + c
    return v


def verify(case):
    rows, founder = rebuild_rows(case)
    encoded = [{j: F(a, b) for j, a, b in raw} for raw in case['rows']]
    assert encoded == rows, 'single-ring reconstruction differs'
    assert all(i == j or q > 0 for i, row in enumerate(rows) for j, q in row.items())
    cores = [1 if x == 'hub' else len(unpack(x)[0]) for x in case['types']]
    minimum = (case['B'] + case['M'] - 2)//(case['M']-1)
    potential = [0 if x == 'hub' else max(0, size-minimum)
                 for x, size in zip(case['types'], cores, strict=True)]
    for i, x in enumerate(case['types']):
        drift = sum(q*cores[j] for j, q in rows[i].items())
        cuts = drift if x != 'hub' and cores[i] == case['B'] else 0
        assert sum(q*potential[j] for j, q in rows[i].items()) + minimum*cuts <= drift
        if x == 'hub':
            assert drift == 1
        else:
            graph, leaf, slots = unpack(x)
            assert drift == sum(F(leaf[u], leaf[u]+slots[u]+len(ns)) for u, ns in graph.items())
    accessible = [reachable(rows, i) for i in range(len(rows))]
    assert len(accessible[founder]) == len(rows)
    equivalence = {frozenset(j for j in a if i in accessible[j])
                   for i, a in enumerate(accessible)}
    listed = [frozenset(g['indices']) for g in case['classes']]
    assert len(listed) == len(set(listed)) and set(listed) == equivalence

    estimates = []
    for component in case['classes']:
        members = component['indices']
        weights = dict(zip(members, map(int, component['weights']), strict=True))
        assert min(weights.values()) > 0
        ratios = [sum(q*weights[j] for j, q in rows[i].items() if j in weights)/F(weights[i])
                  for i in members]
        assert max(ratios)-min(ratios) < F(1, 10**16)
        estimates.append((min(ratios), max(ratios)))
    d = case['bracket_denominator']
    lo, hi = F(case['lower_numerator'], d), F(case['upper_numerator'], d)
    best = case['winner']
    assert lo <= estimates[best][0] <= estimates[best][1] <= hi
    assert all(i == best or upper < lo for i, (_, upper) in enumerate(estimates))

    members = case['classes'][best]['indices']
    scale = case['char_scale']
    coeff = list(map(int, case['char_coefficients']))
    n = len(members)
    assert coeff[0] == 1 and len(coeff) == n+1 and scale > 0
    a = [[rows[i].get(j, F(0))*scale for j in members] for i in members]
    assert all(x.denominator == 1 for row in a for x in row)
    a = [[int(x) for x in row] for row in a]
    for x in range(n+1):
        mat = [[(x if i == j else 0)-a[i][j] for j in range(n)] for i in range(n)]
        assert det_bareiss(mat) == value(coeff, x)
    left, right = scale*lo, scale*hi
    assert value(coeff, left) <= 0 <= value(coeff, right)
    # Interval Horner for the derivative; strict positivity isolates the root.
    lower = upper = F(0)
    for i, c in enumerate(coeff[:-1]):
        cand = (lower*left, lower*right, upper*left, upper*right)
        lower, upper = min(cand)+(n-i)*c, max(cand)+(n-i)*c
    assert lower > 0
    print(f'M={case["M"]} B={case["B"]}: {len(rows)} rows and seam drifts checked; '
          f'{len(listed)} classes verified; spread < 1e-16; '
          f'{n+1} exact determinants; one root in printed interval')
    return len(rows), n+1


if __name__ == '__main__':
    with open(sys.argv[1], encoding='utf-8') as source:
        data = json.load(source)
    totals = [verify(case) for case in data]
    print(f'PASS: {len(data)} matrices, {sum(x for x, _ in totals)} reachable rows, '
          f'{sum(y for _, y in totals)} characteristic-polynomial evaluations')
