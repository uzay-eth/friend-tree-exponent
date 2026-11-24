#!/usr/bin/env python3
"""Coefficientwise cancellation tests (a free symbol per decorated motif)."""

from collections import Counter
from fractions import Fraction as Q
from functools import cache
from itertools import combinations, product


# State: (adjacency, number of ordinary leaves at each core, founder or None).
def seed(edges, leaves, root=None):
    adj = [[] for _ in leaves]
    for v, w in edges:
        adj[v].append(w)
        adj[w].append(v)
    return tuple(map(tuple, adj)), tuple(leaves), root


def gain(t, v):
    adj, leaves, root = t
    return adj, leaves[:v] + (leaves[v] + 1,) + leaves[v + 1 :], root


def promote(t, v):
    adj, leaves, root = t
    assert leaves[v]
    return (adj[:v] + (adj[v] + (len(adj),),) + adj[v + 1 :] + ((v,),),
            leaves[:v] + (leaves[v] - 1,) + leaves[v + 1 :] + (1,), root)


def degree(t, v):
    return len(t[0][v]) + t[1][v]


def distances(adj, v):
    ds = {v: 0}
    queue = [v]
    for w in queue:
        for z in adj[w]:
            if z not in ds:
                ds[z] = ds[w] + 1
                queue.append(z)
    return ds


def update(out, src, q):
    for k, val in src.items():
        out[k] += q * val
        if not out[k]:
            del out[k]


def diff(fn, t, *edits):
    out = Counter()
    for flags in product((0, 1), repeat=len(edits)):
        state = t
        for flag, (edit, v) in zip(flags, edits):
            if flag:
                state = edit(state, v)
        update(out, fn(state), (-1) ** (len(edits) - sum(flags)))
    return out


def local(d):
    @cache
    def score(t):
        adj, _, root = t
        parents = distances(adj, root) if root is not None else None
        metric = [distances(adj, v) for v in range(len(adj))]
        result = Counter()
        for n in range(1, len(adj) + 1):
            for vertices in combinations(range(len(adj)), n):
                if any(metric[v][w] > d for v, w in combinations(vertices, 2)):
                    continue
                if sum(w in vertices for v in vertices for w in adj[v]) != 2 * (n - 1):
                    continue

                def canon(v, prev):
                    branches = [(0 if parents is None else (1 if parents[w] > parents[v] else -1),
                                 canon(w, v)) for w in adj[v] if w != prev and w in vertices]
                    return degree(t, v), len(adj[v]), v == root, tuple(sorted(branches))

                result[min(canon(v, None) for v in vertices)] += 1
        return result
    return score


def terms(t, fn, blank=None, profile='physical'):
    adj, leaves, _ = t
    out = {}
    for v, nbrs in enumerate(adj):
        out['leaf', v] = (leaves[v], diff(fn, t, (gain, v)))
        for w in nbrs:
            out['send', v, w] = (Q(1, degree(t, w)), diff(fn, t, (gain, v)))
        if profile == 'constant':
            def add_core(state, w):
                return promote(gain(state, w), w)
            out['birth', v] = (Q(1), diff(fn, t, (add_core, v)))
        elif leaves[v]:
            out['birth', v] = (Q(leaves[v], degree(t, v)), diff(fn, t, (promote, v)))
    total = Counter()
    for key, (rate, delta) in out.items():
        if key != blank:
            update(total, delta, rate)
    return total


def term_diffs(t, fn, edits, blank=None):
    ledger = {}
    states = list(product((0, 1), repeat=len(edits)))
    for flags in states:
        state = t
        for flag, (edit, v) in zip(flags, edits):
            if flag:
                state = edit(state, v)
        n = len(state[0])
        keys = [('leaf', v) for v in range(n)] + [('birth', v) for v in range(n)]
        keys += [('send', v, w) for v, nbrs in enumerate(state[0]) for w in nbrs]
        sign = (-1) ** (len(edits) - sum(flags))
        for key in keys:
            update(ledger.setdefault(key, Counter()), single(state, fn, key, blank), sign)
    return {k: v for k, v in ledger.items() if v}


def single(t, fn, key, blank):
    out = Counter()
    if key == blank:
        return out
    kind, v, *source = key
    if kind == 'birth' and not t[1][v]:
        return out
    edit = promote if kind == 'birth' else gain
    factor = (Q(t[1][v], degree(t, v)) if kind == 'birth' else
              Q(1, degree(t, source[0])) if kind == 'send' else t[1][v])
    update(out, diff(fn, t, (edit, v)), factor)
    return out


def scaled(fn, q):
    out = Counter()
    update(out, fn, q)
    return out


def dp(t, v):
    d = degree(t, v)
    return -Q(1, d * (d + 1))


def main():
    cases = 0
    for d in (1, 2, 3):
        y, u, v = 0, 1, d + 1
        for root, endpoint_tail, spare in product((None, 1), (False, True), (0, 2)):
            z = v + 1
            edges = [(i, i + 1) for i in range(z)]
            leaves = [1] * (z + 1)
            if endpoint_tail:
                edges.append((z, z + 1))
                leaves.append(1)
            leaves[y] = leaves[z] = spare if endpoint_tail else spare + 1
            if leaves[y] == 0:
                edges.append((y, len(leaves)))
                leaves.append(1)
            leaves[u] = 2
            t = seed(edges, leaves, root)
            phi = local(d)
            edits = ((gain, z), (promote, u))
            assert not diff(phi, t, *edits)
            expected = {('send', v, z): scaled(diff(phi, t, (gain, v), (promote, u)), dp(t, z))}
            assert term_diffs(t, phi, edits) == {k: x for k, x in expected.items() if x}
            edits = ((gain, y), (gain, v))
            assert not diff(phi, t, *edits)
            expected = {
                ('send', u, y): scaled(diff(phi, t, (gain, u), (gain, v)), dp(t, y)),
                ('send', v - 1, v): scaled(diff(phi, t, (gain, v - 1), (gain, y)), dp(t, v)),
            }
            assert term_diffs(t, phi, edits) == {k: x for k, x in expected.items() if x}
            cases += 1
    print(f'D=1,2,3: {cases} decorated cases, all first-round AA and AB coefficients match term by term')

    for d in (1, 2, 3):
        # A symmetric LCA when length is even, central founder-fixed fork
        # around the whole even corridor otherwise. Only the even corridor
        # has an oriented endpoint-flipping automorphism.
        edges = [(i, i + 1) for i in range(2 * d)] + [(d, 2 * d + 1), (2 * d + 1, 2 * d + 2)]
        t = seed(edges, [1] * (2 * d + 3), 2 * d + 2)
        phi = local(d)
        assert phi(gain(t, 0)) == phi(gain(t, 2 * d))
        phi = local(2 * d)
        assert diff(phi, t, (gain, 0), (gain, 2 * d))  # reflection alone is not antisymmetry
        t = t[0], t[1], d  # founder may also be the LCA
        assert phi(gain(t, 0)) == phi(gain(t, 2 * d))
    print('Founder-oriented V symmetry: 3 explicit reflections; reflection alone does not force zero')

    t = seed([(0, 1), (1, 2), (2, 3)], [1, 2, 1, 1], 1)
    n = lambda x: Counter({'n': len(x[0]) + sum(x[1])})
    c = lambda x: Counter({'c': len(x[0])})
    one = lambda x: Counter({'one': 1})
    assert terms(t, n) == n(t) and not terms(t, one)
    assert terms(t, c, profile='constant') == c(t)
    bare = seed([(0, 1), (1, 2)], [1, 0, 1], 1)
    assert terms(bare, c, profile='constant') == c(bare)
    assert terms(bare, c) != c(bare)
    phi = local(1)
    edits = ((gain, 0), (gain, 2))
    assert term_diffs(t, phi, edits) != term_diffs(t, phi, edits, blank=('send', 1, 0))
    print('Controls: size at s=1, constant at s=0, fake constant promotion at s=1, dropped send detected')

    # Degree-free decorated range 2: P(1) holds but BB need not vanish.
    def topology(x):
        adj = x[0]
        return Counter({'topology': len(adj) + sum(min(len(adj[i]), 3) * min(len(adj[j]), 3)
                                                  for nbrs in adj if len(nbrs) == 2 for i, j in [nbrs])})
    ledger = term_diffs(t, topology, ((gain, 2), (promote, 0)))
    assert set(ledger) == {('birth', 2)} and ledger['birth', 2]['topology'] > 0
    print(f'Lower-round gap witness: surviving (Δ_z b_z) B_z B_u = {ledger["birth", 2]["topology"]}')

    # Compare local motifs differing only in a core decoration off the tested path.
    # The path endpoints' mixed increment sees that decoration, despite P(2).
    def fork(x):
        adj = x[0]
        value = 0
        for nbrs in adj:
            for arms in combinations(nbrs, 3):
                for i, j in combinations(arms, 2):
                    w = next(v for v in arms if v not in (i, j))
                    value += degree(x, i) * degree(x, j) * len(adj[w])
        return Counter({'fork': value})

    f1 = seed([(1, 0), (1, 2), (1, 3)], [1, 1, 1, 2], 1)
    f2 = seed([(1, 0), (1, 2), (1, 3), (3, 4)], [1, 1, 1, 1, 1], 1)
    edits = ((gain, 0), (gain, 2))
    increments = tuple(diff(fork, x, *edits)['fork'] for x in (f1, f2))
    assert increments == (1, 2)
    print(f'Diameter-2 fork: same decorated path, mixed increments {increments}')

    @cache
    def first(x):
        return terms(x, fork)

    @cache
    def second(x):
        return terms(x, first)

    t = seed([(0, 1), (1, 2), (2, 3), (1, 4)], [1] * 5, 1)
    edits = ((gain, 0), (gain, 3))
    assert not diff(fork, t, *edits)
    print('D=2 fixed distance 3: sample (AA Φ, AA LΦ, AA L²Φ) =',
          tuple(str(diff(fn, t, *edits)['fork']) for fn in (fork, first, second)))


if __name__ == '__main__':
    main()
