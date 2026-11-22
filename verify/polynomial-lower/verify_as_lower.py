"""Enumerate one-step moves in exact arithmetic on sampled and targeted trees."""

import argparse
from fractions import Fraction as Q
import os
from pathlib import Path
import random
import sys


class Tree:
    def __init__(self, parents=(-1, 0, 0), labels=(3, 3, 3)):
        self.parents, self.labels = parents, labels
        self.adj = [[] for _ in parents]
        for v, p in enumerate(parents[1:], 1):
            self.adj[p].append(v)
            self.adj[v].append(p)
        self.degree = tuple(map(len, self.adj))
        self.core = {v for v, d in enumerate(self.degree) if d >= 2}
        self.leaves = [sum(w not in self.core for w in row) for row in self.adj]
        self.counts = [d - l for d, l in zip(self.degree, self.leaves)]

    def attach(self, receiver):
        labels = list(self.labels) + [3]
        if receiver not in self.core:
            p = self.parents[receiver]
            after = (self.leaves[p] - 1, self.counts[p] + 1)
            if after not in ((1, 1), (0, 2)):
                labels[receiver] = min(self.degree[p], 20)
        return Tree(self.parents + (receiver,), tuple(labels))

    def branch(self, founder):
        inside = {founder}
        for v in range(founder + 1, len(self.parents)):
            if self.parents[v] in inside:
                inside.add(v)
        return inside

    def classes(self, u, pendant_cap, degree_cap):
        pendants = [0] * pendant_cap
        tracked = [0] * (degree_cap - 3)
        class_b = 0
        for v in self.adj[u]:
            if v not in self.core:
                continue
            l, c, d = self.leaves[v], self.counts[v], self.degree[v]
            if c == 1 and 1 <= l <= pendant_cap:
                pendants[l - 1] += 1
            elif (l, c) == (0, 2):
                class_b += 1
            elif d < degree_cap:
                tracked[d - 3] += 1
        return tuple(pendants), class_b, tuple(tracked)

    def states(self, u):
        base = (self.leaves[u], self.counts[u])
        tracked = base + self.classes(u, 2, 6)
        pendants, class_b, _ = self.classes(u, 1, 3)
        static = base + (pendants, class_b, self.labels[u])
        p = self.parents[u]
        parent = None if p < 0 else "top" if self.degree[p] >= 7 else self.peek(p)
        children = sorted(self.peek(v) for v in self.adj[u]
                          if v != p and 2 <= self.degree[v] < 5)
        rooted = base + (parent, tuple(children))
        return tracked, static, rooted

    def peek(self, v):
        return self.leaves[v], self.counts[v]


def probe(state, cutoff):
    # Positive, deliberately nonmonotone rational weights for testing structural bounds.
    l, c, *rest = state
    checksum = sum((i + 1) * ord(ch) for i, ch in enumerate(repr((min(l, cutoff), *rest))))
    return Q(0) if c > 10 else Q(1 + checksum % 17, 17)


def report_tables():
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "baseline"))
    import certs6
    import certs9

    lower_model = certs9.Model(60, 10, 2, 6)
    lower, lower_root = certs9.lower_table(lower_model, 0.5433)
    upper_model = certs6.Model("upper", 60, 10, 1, 20)
    (upper, outer), upper_root = certs6.upper_lin_table(upper_model, 0.586)
    assert lower_root > 1 and upper_root < 1
    return lower, upper, Q(outer)


def values(tree, tables):
    tracked, rooted, report_lower, report_upper = {}, {}, {}, {}
    for u in tree.core:
        lo, up, root = tree.states(u)
        # L_max = D_cap - 1 exercises the edge case in the sharp jump bound.
        tracked[u], rooted[u] = probe(lo, 5), probe(root, 4)
        if tables:
            lower, upper, outer = tables
            l, c, *rest = lo
            report_lower[u] = Q(lower[min(l, 60), c, *rest]) if c <= 10 else Q(0)
            l, c, *rest = up
            report_upper[u] = Q(upper[min(l, 60), c, *rest]) if c <= 10 else outer
    return (tracked, rooted, report_lower, report_upper) if tables else (tracked, rooted)


def verify_rate_intervals(tree, u, rho, rates):
    degree = tree.degree[u]
    for pendant_cap, degree_cap in ((2, 6), (1, 3)):
        pendants, class_b, tracked = tree.classes(u, pendant_cap, degree_cap)
        known = sum(Q(k, l + 1) for l, k in enumerate(pendants, 1)) + Q(class_b, 2)
        known += sum(Q(k, d) for d, k in enumerate(tracked, 3))
        unknown = tree.counts[u] - sum(pendants) - class_b - sum(tracked)
        upper = known + Q(unknown, degree_cap)
        label = tree.labels[u]
        if pendant_cap == 1 and label >= 4:
            assert unknown >= 1 and tree.degree[tree.parents[u]] >= label
            upper += Q(1, label) - Q(1, 3)
        assert known <= rho[u] <= upper
    for v in tree.adj[u]:
        if v in tree.core:
            if tree.peek(v) == (0, 2):
                assert Q(1, degree) <= rates[v] <= Q(1, degree) + Q(1, 2)
            assert 0 <= rates[v] <= tree.degree[v]


def verify_tree(tree, founders, tables):
    n = len(tree.parents)
    regions = [set(range(n))] + [tree.branch(v) for v in founders]
    rates = [Q(0)] * n
    for v, row in enumerate(tree.adj):
        for w in row:
            rates[w] += Q(1, tree.degree[v])
    assert sum(rates) == n
    rho = dict.fromkeys(tree.core, Q(0))
    for u in tree.core:
        rho[u] = sum(Q(1, tree.degree[v]) for v in tree.adj[u] if v in tree.core)
        assert rates[u] == tree.leaves[u] + rho[u]
        assert 0 <= rho[u] <= Q(tree.counts[u], 2)
        if tree.leaves[u]:
            assert rates[u] >= 1 and rates[u] >= Q(tree.leaves[u], tree.degree[u])
            assert Q(tree.leaves[u] + 1, tree.degree[u] + 1) >= Q(tree.leaves[u], tree.degree[u])
    assert sum(rho.values()) == sum(Q(tree.counts[u], tree.degree[u]) for u in tree.core)
    assert sum(rho.values()) <= len(tree.core)
    for w in regions[0] - tree.core:
        assert rates[w] == Q(1, tree.degree[tree.parents[w]])
    for u in tree.core:
        verify_rate_intervals(tree, u, rho, rates)
    thresholds = (5, 7, 60, 60)
    entry_limits = (6, 8, 6, 3)
    jump_bounds = (6, 8, 6, 18)
    base = values(tree, tables)
    if tables:
        assert all(0 <= x <= 1 for x in base[2].values())
        assert all(1 <= x <= 6 for x in base[3].values())
    local = [dict.fromkeys(tree.core, Q(0)) for _ in base]
    drift = [[Q(0)] * len(base) for _ in regions]
    square = [[Q(0)] * len(base) for _ in regions]
    for w in range(n):
        after = tree.attach(w)
        for inside in regions[1:]:
            if w not in inside:
                assert all(tree.states(u) == after.states(u) for u in inside & tree.core)
        future = values(after, tables)
        creator = tree.parents[w] if w not in tree.core else None
        for i, (old, new) in enumerate(zip(base, future)):
            changes = {u: x - old.get(u, 0) for u, x in new.items() if x != old.get(u, 0)}
            assert len(changes) <= entry_limits[i]
            if w in tree.core and tree.leaves[w] >= thresholds[i]:
                assert not changes, (i, n, w, changes)
            assert abs(sum(changes.values())) <= jump_bounds[i]
            for u, change in changes.items():
                local[i][creator if u == w and creator is not None else u] += rates[w] * change
            for j, inside in enumerate(regions):
                inside_changes = {u: change for u, change in changes.items() if u in inside}
                if w not in inside:
                    assert not inside_changes, (i, n, w, inside_changes)
                jump = sum(inside_changes.values())
                assert abs(jump) <= jump_bounds[i]
                drift[j][i] += rates[w] * jump
                square[j][i] += rates[w] * jump**2

    if tables:
        for u in tree.core:
            assert local[2][u] >= Q(543, 1000) * base[2][u]
            assert local[3][u] <= Q(293, 500) * base[3][u]
    for j, inside in enumerate(regions):
        core = tree.core & inside
        promotions = sum(rates[w] for w in inside - core)
        assert promotions == sum(Q(tree.leaves[u], tree.degree[u]) for u in core) <= len(core)
        credits = sum(Q(1, tree.degree[v]) for u in core for v in tree.adj[u] if v in core)
        if j:
            p = tree.parents[founders[j - 1]]
            assert tree.degree[p] >= 60
            credits += Q(1, tree.degree[p])
            assert min(rates[w] for w in inside) >= Q(1, len(inside))
            assert sum(rates[w] for w in inside) <= len(inside) + 1
        assert sum(rho[u] for u in core) == credits <= len(core) + bool(j)
        for i, old in enumerate(base):
            # Sum the exact local terms, with newborn weight assigned to its creator.
            assert drift[j][i] == sum(local[i][u] for u in core)
            threshold = thresholds[i]
            active_gains = sum(rates[u] for u in core if tree.leaves[u] < threshold)
            assert active_gains <= (threshold + 1 + bool(j)) * len(core)
            assert square[j][i] <= jump_bounds[i]**2 * (promotions + active_gains)
            assert square[j][i] - drift[j][i]**2 / n >= 0
            potential = sum(old[u] for u in core)
            if tables and i == 2:
                assert drift[j][i] >= Q(543, 1000) * potential
            if tables and i == 3:
                assert drift[j][i] <= Q(293, 500) * potential


def random_tree(rng, size, tree=None):
    tree = tree or Tree()
    while len(tree.parents) < size:
        v = rng.randrange(len(tree.parents))
        tree = tree.attach(rng.choice(tree.adj[v]))
    return tree


def verify_graft(pattern, star):
    birth = len(star.parents)
    graft = star.attach(1).attach(1)
    mapped = [1, birth, birth + 1]
    for w in pattern.parents[3:]:
        mapped.append(len(graft.parents))
        graft = graft.attach(mapped[w])
    for u in pattern.core:
        if pattern.parents[u] > 0:
            assert graft.states(mapped[u]) == pattern.states(u)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-tables", action="store_true", help="also check the tables from REPORT.md")
    args = parser.parse_args()
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    tables = report_tables() if args.report_tables else None
    rng = random.Random(543586)
    star = Tree()
    while len(star.parents) < 66:
        star = star.attach(0)
    assert star.degree[0] >= 60
    small = Tree().attach(1)
    assert small.states(1) != small.attach(0).states(1)  # an unsaturated boundary can change
    founders = tuple(range(1, 5))
    siblings = star
    for v in founders:
        siblings = siblings.attach(v)
    verify_tree(star, (), tables)
    verify_tree(siblings, founders, tables)
    for _ in range(80):
        tree = random_tree(rng, rng.randrange(3, 71))
        verify_tree(tree, (), tables)
        verify_graft(tree, star)
    for _ in range(20):
        tree = random_tree(rng, rng.randrange(71, 101), siblings)
        verify_tree(tree, founders, tables)
    print("PASS: exact receiver rates, local drift sums, jumps, and second moments on 102 trees")
    print("PASS: exterior invariance and replay rate bounds on 84 saturated branches; 80 grafts agree")
    if tables:
        print("PASS: exact local and restricted report-table drift at delta=543/1000 and lambda=293/500")


if __name__ == "__main__":
    main()
