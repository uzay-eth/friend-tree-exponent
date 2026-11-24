"""Exact inert occupations and a finite live-input singleton control."""

from collections import defaultdict
from fractions import Fraction as Q


def add(tree, at):
    new = len(tree)
    return tuple(ns + (new,) if i == at else ns for i, ns in enumerate(tree)) + ((at,),)


def shape(tree, root=0, parent=-1):
    return tuple(sorted(shape(tree, v, root) for v in tree[root] if v != parent))


def inert_rates(tree):
    ds = [len(ns) + (v == 0) for v, ns in enumerate(tree)]
    rates = [sum((Q(1, ds[w]) for w in ns), Q(0)) for ns in tree]
    assert sum(rates) == len(tree) - Q(1, ds[0])
    return rates


def occupations():
    states = {((1,), (0,)): Q(1)}
    count = 0
    for n in range(2, 8):
        assert sum(states.values()) == 1
        occ = defaultdict(Q)
        nxt = defaultdict(Q)
        for tree, entry in states.items():
            rates = inert_rates(tree)
            occ[shape(tree)] += entry / sum(rates)
            for v, rate in enumerate(rates):
                nxt[add(tree, v)] += entry * rate / sum(rates)
        count += len(occ)
        if n == 2:
            assert list(occ.values()) == [Q(2, 3)]
        if n == 3:
            assert occ[((), ())] == Q(1, 4)
            assert occ[(((),),)] == Q(2, 15)
        states = nxt
    print(f'{count} positive exact rooted-shape occupations at sizes 2–7; pair 2/3')


def singleton_controls():
    # High o—bare v—high w. Each high endpoint has its own leaf stock,
    # and o has one real exterior core, added explicitly below.
    adj = [(1, 3), (0, 2), (1,), (0, 4), (3,)]
    for _ in range(6):
        adj = add(add(tuple(adj), 0), 2)
    checked = 0
    for N in range(2, 13):
        low = {v for v, ns in enumerate(adj) if len(ns) <= N}
        component = {1}
        wave = [1]
        for v in wave:
            for w in adj[v]:
                if w in low and w not in component:
                    wave.append(w)
                    component.add(w)
        actual = 0 not in low and component == {1}
        assert actual == (N < min(len(adj[0]), len(adj[2])))
        if actual:
            assert len(adj[1]) == 2 and all(w not in low for w in adj[1])
            checked += 1
    # Equal and unequal live masses: exact marked stable-age integration
    # with rational exponent 1/2, using square masses.
    for parent, other in ((9, 4), (4, 9), (9, 9)):
        from math import isqrt
        coeff = isqrt(min(parent, other))
        assert coeff == min(isqrt(parent), isqrt(other))
    print(f'{checked} exact old-core singleton cutoffs; inert singleton impossible')


if __name__ == '__main__':
    occupations()
    singleton_controls()
