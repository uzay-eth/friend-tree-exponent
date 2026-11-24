"""Independent finite witnesses for the low-forest cut and young Palm law."""

from collections import Counter, defaultdict
from fractions import Fraction as Q


def tree(edges):
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def cut(adj, n):
    live = {v for v in adj if len(adj[v]) <= n}
    unseen = set(live)
    components = []
    while unseen:
        v = unseen.pop()
        seen, todo = {v}, [v]
        while todo:
            for w in adj[todo.pop()] & unseen:
                unseen.remove(w)
                seen.add(w)
                todo.append(w)
        components.append(frozenset(seen))
    return components


def attached(adj, target):
    out = {v: set(nb) for v, nb in adj.items()}
    new = max(out) + 1
    out[target].add(new)
    out[new] = {target}
    return out


def rates(adj):
    return {v: sum((Q(1, len(adj[w])) for w in nb), Q(0))
            for v, nb in adj.items()}


def degrees(adj, vertices):
    return tuple(sorted(len(adj[v] & vertices) for v in vertices))


def check_cut():
    # High cores 0 and 1, isolated core 2 between them, a two-core path
    # via 3 and 4, and a one-core path via 5. All other vertices are leaves.
    adj = tree([(0, 2), (2, 1), (0, 3), (3, 4), (4, 6), (1, 5),
                (5, 7), (0, 8), (1, 9)])
    components = cut(adj, 2)
    assert set(components) == {frozenset(x) for x in [(2,), (3, 4, 6),
                                                     (5, 7), (8,), (9,)]}
    cores = {v for v in adj if len(adj[v]) >= 2}
    direct, weighted, sampled = Counter(), Counter(), Counter()
    for comp in components:
        k = len(comp & cores)
        if k:
            direct[len(comp)] += Q(1, len(cores))
            for _ in comp & cores:
                weighted[len(comp)] += Q(1, len(cores) * k)
                sampled[len(comp)] += Q(1, len(cores))
    assert direct == weighted == {1: Q(1, 6), 2: Q(1, 6), 3: Q(1, 6)}
    assert sampled == {1: Q(1, 6), 2: Q(1, 6), 3: Q(1, 3)}
    assert sum(len(c) > 1 for c in components) == 2

    # Old anchor 0 is already high at N=3. Its two new core branches
    # cannot connect through it, but one has a path turning below it.
    old = tree([(0, 1), (0, 2), (0, 3), (0, 4)])
    now = tree([(0, w) for w in old[0]] + [(0, 5), (0, 6), (5, 7),
               (5, 8), (7, 9), (8, 10), (6, 11)])
    from_seven = next(c for c in cut(now, 3) if 7 in c)
    assert from_seven == {5, 7, 8, 9, 10}
    assert all(len(now[v]) >= 2 for v in (7, 5, 8))
    assert {w for v in from_seven for w in now[v] - from_seven} == {0}
    return len(components)


def check_palm(d):
    # This is a full physical tree, including all of the parent's leaves.
    # Enumerate by attachment target, allowing both exterior and interior
    # births. With two total births the size-four cut cannot use an exterior one.
    adj = tree([(0, 1), (1, 2)] + [(0, i) for i in range(3, d + 2)])
    totals = defaultdict(Q)
    for x, a in rates(adj).items():
        middle = attached(adj, x)
        for y, b in rates(middle).items():
            final = attached(middle, y)
            comp = next(c for c in cut(final, 4) if 1 in c)
            if len(comp) == 4:
                totals[degrees(final, comp)] += a * b
    u = 1 + Q(1, d)
    star, path = (1, 1, 1, 3), (1, 1, 2, 2)
    assert totals == {star: u * (u + 1) + Q(3, 4), path: 7 * u / 6}
    ratio = totals[star] / sum(totals.values())
    assert ratio > Q(2, 3)
    return ratio


if __name__ == '__main__':
    print(f'cut: {check_cut()} components; inverse-core weighting and high-anchor turn pass')
    three = tree([(0, 1), (1, 2)])
    ordinary = rates(three)[1] / sum(rates(three).values())
    assert ordinary == Q(2, 3)
    results = {d: check_palm(d) for d in (5, 8, 20, 64)}
    assert results[20] == Q(1161, 1651)
    print(f'Palm: {results}; ordinary size-four star={ordinary}')
