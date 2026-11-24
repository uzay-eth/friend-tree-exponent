"""Exact finite graph checks for the parent census and the cylinder jump bound."""

from collections import Counter
from fractions import Fraction as Q


def graph(par):
    adj = [set() for _ in range(len(par) + 1)]
    for v, p in enumerate(par, 1):
        adj[p].add(v)
        adj[v].add(p)
    return adj


def census(par, cutoff, height, cap):
    adj = graph(par)
    low = {v for v, ns in enumerate(adj) if len(ns) <= cutoff}
    cores = {v for v, ns in enumerate(adj) if len(ns) >= 2}

    def piece(v, previous):
        kids = [piece(w, v) for w in adj[v] & low if w != previous]
        return (v in cores, tuple(sorted(kids)))

    def downward(v):
        return (v in cores, tuple(sorted(downward(w) for w in adj[v] & low if w > v)))

    def volume(mark):
        return int(mark[0]) + sum(volume(sub) for sub in mark[1])

    def core_height(mark):
        return 1 + max((core_height(sub) for sub in mark[1] if sub[0]), default=0)

    direct = Counter()
    survivors = Counter()
    omitted = 0
    unseen = set(low)
    while unseen:
        r = min(unseen)
        walk = [r]
        unseen.remove(r)
        for v in walk:
            for w in adj[v] & unseen:
                unseen.remove(w)
                walk.append(w)
        count = len(cores.intersection(walk))
        if count == 0:
            continue
        shape = piece(r, -1)
        assert count == volume(shape)
        if not r or core_height(shape) > height or len(adj[par[r - 1]]) > cap:
            omitted += count
            continue
        mark = (Q(len(adj[par[r - 1]]), cutoff), shape)
        direct[mark] += 1
        survivors[mark] += count

    parent = Counter()
    weighted = Counter()
    mutants = 0
    for v in cores - low:
        if len(adj[v]) > cap:
            continue
        for w in adj[v] & low & cores:
            if w <= v:
                mutants += 1
                continue
            shape = downward(w)
            if core_height(shape) > height:
                continue
            mark = (Q(len(adj[v]), cutoff), shape)
            parent[mark] += 1
            weighted[mark] += volume(shape)
    assert parent == direct and weighted == survivors
    assert sum(weighted.values()) + omitted == len(low & cores)
    return parent, weighted, omitted, mutants


def main():
    trees = [(0, 0)]  # T3 rooted at its centre
    checked = founder_or_tail = wrong_orientation = changes = 0
    for n in range(3, 8):
        for par in trees:
            for cutoff in (2, 3, 4):
                for height, cap in ((1, 4), (2, 8), (3, 16)):
                    counts, weights, omitted, bad = census(par, cutoff, height, cap)
                    checked += 1
                    founder_or_tail += omitted > 0
                    wrong_orientation += bad
                    if n < 7:
                        for p in range(n):
                            nxt, nxt_w, *_ = census(par + (p,), cutoff, height, cap)
                            for before, after, bound in (
                                (counts, nxt, cap),
                                (weights, nxt_w, cap * (cutoff + 1) ** (height + 2)),
                            ):
                                norm = sum(abs(before[z] - after[z]) for z in before.keys() | after.keys())
                                assert norm <= 6 * (height + 4) * bound
                                changes += norm > 0
        trees = [par + (p,) for par in trees for p in range(n)]
    assert min(founder_or_tail, wrong_orientation, changes) > 0
    print(f'{checked} exact parent/survivor censuses; {founder_or_tail} nonempty residuals; '
          f'{wrong_orientation} reversed-edge controls; {changes} nonzero one-birth changes')


if __name__ == '__main__':
    main()
