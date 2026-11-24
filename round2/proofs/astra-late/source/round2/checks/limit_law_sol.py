#!/usr/bin/env python3
"""Independent exact controls for the free coupling, cells and reachable histories."""
from collections import Counter, defaultdict
from fractions import Fraction as F


def adjacency(parents):
    g = {v: [] for v in parents}
    for v, p in parents.items():
        if p is not None:
            g[v].append(p)
            g[p].append(v)
    return g


def physical(parents):
    g = adjacency(parents)
    return {(v, w): F(1, len(ns)) for v, ns in g.items() for w in ns}


def coupled(trees, swap=False):
    a, b = map(physical, trees)
    for edge in sorted(a.keys() | b.keys()):
        x, y = a.get(edge, F(0)), b.get(edge, F(0))
        both = min(x, y)
        if both:
            yield (0, 1), edge, both
        for side, rate in enumerate((x - both, y - both)):
            if rate:
                yield (1 - side if swap and both else side,), edge, rate


def marginal(events, side):
    out = Counter()
    for sides, edge, rate in events:
        if side in sides:
            out[edge] += rate
    return out


def advance(trees, owners, event):
    sides, (_, target), _ = event
    newborn = 1 + max(set(trees[0]) | set(trees[1]))
    infect = len(sides) == 1 and target in trees[1 - sides[0]] and owners[sides[0]][target] is None
    new_trees, new_owners = tuple(map(dict.copy, trees)), tuple(map(dict.copy, owners))
    for side in sides:
        own = owners[side][target]
        new_trees[side][newborn] = target
        new_owners[side][newborn] = own if own is not None else (
            newborn if infect or degrees_differ(trees, target) else None
        )
    return new_trees, new_owners


def degrees_differ(trees, x):
    return x in trees[0] and x in trees[1] and (
        len(adjacency(trees[0])[x]) != len(adjacency(trees[1])[x])
    )


def verify_partition(trees, owners):
    common = trees[0].keys() & trees[1].keys()
    for x in common:
        assert (owners[0][x] is None) == (owners[1][x] is None)
        assert trees[0][x] == trees[1][x]
    infected = {x for x in common if degrees_differ(trees, x)}
    assert 0 in infected
    # Every path toward the seed preserves infection and parity-favoured side.
    for x in infected:
        p = trees[0][x]
        if p is not None:
            assert p in infected
    graphs = tuple(map(adjacency, trees))
    for side, (tree, own, g) in enumerate(zip(trees, owners, graphs)):
        for x in common:
            path = x
            depth = 0
            while tree[path] is not None:
                path = tree[path]
                depth += 1
            if any(v not in common for v in g[x]):
                assert side == depth % 2
        accounts = defaultdict(set)
        for x, account in own.items():
            if account is not None:
                accounts[account].add(x)
            else:
                assert x in common
        for root, members in accounts.items():
            assert root in members
            cut = {(x, y) for x in members for y in g[x] if y not in members}
            assert cut == {(root, tree[root])}
            assert tree[root] in infected and own[tree[root]] is None
    count = [sum(len(g[x]) >= 2 for x in g) for g in graphs]
    tracked = [sum(len(g[x]) >= 2 for x in g if own[x] is not None or x in infected)
               for g, own in zip(graphs, owners)]
    assert count[0] - count[1] == tracked[0] - tracked[1]


def pair_checks():
    small = {0: None, 1: 0, 2: 1}
    trees = (small | {3: 0}, small)
    owners = ({0: None, 1: None, 2: None, 3: 3}, {v: None for v in small})
    original = trees, owners
    wrong = list(coupled(trees, swap=True))
    assert (marginal(wrong, 0)[(0, 1)], marginal(wrong, 1)[(0, 1)]) == (1, F(1, 2))
    assert (physical(trees[0])[(0, 1)], physical(trees[1])[(0, 1)]) == (F(1, 2), 1)
    layer = [original]
    total = 0
    for depth in range(4):
        nxt = []
        for ts, os in layer:
            events = list(coupled(ts))
            for side in (0, 1):
                assert marginal(events, side) == physical(ts[side])
            verify_partition(ts, os)
            total += 1
            if depth < 3:
                nxt.extend(advance(ts, os, ev) for ev in events)
        layer = nxt
    print(f'PASS: {total} paired prefixes; swapped 0->1 rates (1,1/2) != (1/2,1)')

    def step(sides, edge):
        nonlocal trees, owners
        ev, = (ev for ev in coupled(trees) if ev[:2] == (sides, edge))
        trees, owners = advance(trees, owners, ev)
        verify_partition(trees, owners)
        return ev[2]

    assert step((0, 1), (2, 1)) == 1  # 4: newer than time zero, older than infection of 1
    assert step((1,), (0, 1)) == F(1, 2)  # 5: incoming private account at anchor 1
    assert step((0,), (1, 4)) == F(1, 12)  # 6: 4 is a retained singleton, its child is owned
    assert owners[0][4] is None and owners[0][6] == 6
    assert step((0, 1), (2, 1)) == 1  # 7: paired accounts at already-infected anchor 1
    assert step((0,), (1, 7)) == F(1, 20)  # 8: first promotion in A's account 7
    assert owners[0][7] == owners[0][8] == owners[1][7] == 7
    assert degrees_differ(trees, 7) and len(adjacency(trees[0])[7]) == 2
    assert step((0, 1), (2, 1)) == 1  # 9: another paired account after 7's promotion
    assert step((0,), (1, 9)) == F(1, 30)  # 10: independent accounting entry, not a retained child
    assert owners[0][9] == owners[0][10] == owners[1][9] == 9
    assert degrees_differ(trees, 9)
    assert step((0,), (7, 8)) == F(1, 2)  # 11: promotion *inside* account 7
    assert owners[0][11] == 7 and len(adjacency(trees[0])[8]) == 2
    # Omitting the recurring class opened by paired births at an infected
    # anchor drops both terms (7 and 9); one of them already has two cores.
    assert sum(owners[0][v] in (7, 9) and len(adj) >= 2
               for v, adj in adjacency(trees[0]).items()) == 3
    print('PASS: pre-infection birth 4 retained at 1/12; recurrent paired accounts 7 and 9 promoted at 1/20 and 1/30, missing cost 3')


def branch(tree, root):
    out = {root}
    for x in tree:  # parents precede children in every labelled fixture
        if tree[x] in out:
            out.add(x)
    return out


def cell_checks():
    # Anchor 0, old leaves 1 and 2; branch 3 has just been promoted by 0.
    seed = {0: None, 1: 0, 2: 0, 3: 0, 4: 3}

    def generator(f, live=True):
        def apply(tree):
            got = F(0)
            for edge, rate in physical(tree).items():
                if not live and edge == (0, 3):
                    continue
                got += rate * (f(tree | {max(tree) + 1: edge[1]}) - f(tree))
            return got
        return apply

    def cores(tree):
        g = adjacency(tree)
        return sum(len(g[v]) >= 2 for v in branch(tree, 3))

    live, cut = generator(cores), generator(cores, False)
    assert live(seed) == cut(seed) == F(1, 2)
    assert generator(live)(seed) - generator(cut, False)(seed) == F(1, 18)
    assert physical(seed)[0, 3] == F(1, 3)
    assert physical(seed | {5: 0})[0, 3] == F(1, 4)
    layer = [(seed, {0: 0, 3: 3})]
    total = 0
    for depth in range(4):
        nxt = []
        for tree, labels in layer:
            g = adjacency(tree)
            assert {v for v, label in labels.items() if label == 3} == {
                v for v in branch(tree, 3) if len(g[v]) >= 2
            }
            total += 1
            if depth == 3:
                continue
            for source, target in physical(tree):
                new_labels = labels.copy()
                if len(g[target]) == 1:
                    new_labels[target] = labels[source]
                nxt.append((tree | {max(tree) + 1: target}, new_labels))
        layer = nxt
    print(f'PASS: live cell labels on {total} arrow prefixes; muted comparison has equal first derivative 1/2, second-derivative gap 1/18; incoming rate 1/3 -> 1/4')


def reachability_checks():
    layer = [({0: None, 1: 0, 2: 0}, F(1))]
    total = 0
    for n in range(3, 9):
        next_layer = []
        for tree, weight in layer:
            g = adjacency(tree)
            c = Counter(map(len, g.values()))
            assert max(c) <= 3 and c[3] == c[1] - 2 and 2*(n-c[1]) >= n-2
            probs = Counter()
            for (_, target), rate in physical(tree).items():
                probs[target] += rate / n
            prefer = 1 if 4*c[1] >= n else 2
            allowed = [v for v in g if len(g[v]) == prefer]
            kept = sum(probs[v] for v in allowed)
            assert kept >= F(1, 12)
            next_layer.extend((tree | {n: v}, weight*probs[v]) for v in allowed)
            total += 1
        next_mass = sum(w for _, w in next_layer)
        assert next_mass >= sum(w for _, w in layer) / 12
        layer = next_layer
    print(f'PASS: {total} adaptively selected degree-cap prefixes; size-9 retained mass {next_mass}')

    n0 = 16
    tree = {v: (0 if v else None) for v in range(n0)}
    # All leaf-driven star choices and fixed arm steps are actual positive-probability transitions.
    star_prob = F(1)
    for n in range(3, n0):
        star_prob *= sum(r/n for (_, y), r in physical({v: (0 if v else None)
                                                         for v in range(n)}).items() if y == 0)
    assert star_prob == F(2, n0 - 1)
    endpoint = 1
    for n in range(n0, n0 + 24):
        assert sum(len(x) >= 2 for x in adjacency(tree).values()) == 1 + n - n0
        source, = adjacency(tree)[endpoint]
        assert physical(tree)[source, endpoint]/n == F(1, n*((n0-1) if n == n0 else 2))
        tree[n], endpoint = endpoint, n

    # Deterministic, exact mean-preserving control: replace a reachable near-x
    # conditional law by endpoints while keeping its conditional mean. Neither
    # mean nor a fine grid of means forces positive mass in the target interval.
    a, b, target, radius = F(0), F(2), F(1), F(1, 8)
    for i in range(29, 36):
        y = F(i, 32)
        p_right = (y-a)/(b-a)
        assert (1-p_right)*a+p_right*b == y
        assert abs(y-target) < radius and min(abs(a-target), abs(b-target)) > radius
        var = (1-p_right)*(a-y)**2+p_right*(b-y)**2
        assert var == (y-a)*(b-y) >= F(1015, 1024)
    # Conversely, deterministic support with tiny variance at the only visited
    # states can omit a neighbourhood if the interpolating prefix is unreachable.
    for y in (a, b):
        assert abs(y-target) > radius and (y-y)**2 == 0
    print('PASS: star probability 2/15 and 24 arm steps; exact mean-preserving two-point gap control')


if __name__ == '__main__':
    pair_checks()
    cell_checks()
    reachability_checks()
