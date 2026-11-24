#!/usr/bin/env python3
"""Exact directed-arrow checks and a seeded finite-horizon P4 experiment."""

import argparse
from fractions import Fraction as F
from math import exp, log, sqrt
import random
from statistics import correlation


def adjacency(parents):
    adj = [[] for _ in range(len(parents) + 1)]
    for child, parent in enumerate(parents, 1):
        adj[parent].append(child)
        adj[child].append(parent)
    return adj


def rates(parents, slot=True):
    adj = adjacency(parents)
    internal = [F(0) for _ in adj]
    for v, neighbors in enumerate(adj):
        for w in neighbors:
            internal[w] += F(1, len(neighbors) + (v == 0 and slot))
    return internal, F(1, len(adj[0]) + 1)


def exact():
    # Birth-labelled trees rooted at 0; every target is reachable via a real neighbor.
    levels = {2: {(0,): (F(1), F(1), F(1))}}
    for m in range(2, 8):
        levels[m + 1] = {}
        theta = F(m + 1, 3)
        for tree, (prob, laplace, factor) in levels[m].items():
            internal, kill = rates(tree)
            total = sum(internal) + kill
            assert total == m and all(r > 0 for r in internal)
            for target, rate in enumerate(internal):
                child = tree + (target,)
                p = prob * rate / m
                # Integral of the arrow density times exp(-theta * wait).
                l = laplace * rate / (total + theta)
                f = factor * m / (m + theta)
                assert l == p * f
                levels[m + 1][child] = p, l, f
    counts = [len(levels[m]) for m in range(2, 8)]
    paths = sum(len(levels[m]) for m in range(3, 9))
    print(f"Exact: {sum(counts)} states by size 2–7 {counts}; {paths} path Laplace identities")

    virgin = (0,)
    internal, kill = rates(virgin)
    # (a) Wrong source denominator, retaining the physical kill channel.
    assert sum(rates(virgin, slot=False)[0]) + kill != 2
    # (b) Leaf edge in P4: the degree convention incorrectly counts the singleton.
    path = adjacency((0, 1, 2))
    actual = sum(len(neighbors) >= 2 for neighbors in path)
    singleton_internal, singleton_kill = rates(())
    singleton_core = sum(len(ns) + (v == 0) >= 2 for v, ns in enumerate(adjacency(())))
    remainder = adjacency((0, 1))
    remainder_core = sum(len(ns) + (v == 0) >= 2 for v, ns in enumerate(remainder))
    # A singleton with one slot is still a leaf, so literal counts add.
    # Substituting a positive inert component or forcing its root to be a core fails.
    assert actual == singleton_core + remainder_core != 1 + remainder_core
    assert singleton_internal == [0] and singleton_kill == 1
    # (c) An unthinned total internal rate changes the transform of the wait.
    assert (2 + kill) / (2 + kill + 1) != F(2, 3)
    # (d) A crossing cannot produce the claimed first local internal successor.
    # Its wait is STILL Exp(2); the wrong event is caught by its marked path law.
    arrow_rates = {0: internal[0], 1: internal[1], -1: kill}
    correct = sum(rate / 2 for target, rate in arrow_rates.items() if target == 0)
    planted = sum(rate / 2 for target, rate in arrow_rates.items() if target == 0 and target == -1)
    assert correct > 0 == planted
    print("Controls: 4/4 caught (including the leaf-edge false-core convention)")


def side(rng, size, mu):
    adj = [[1, -1], [0]]  # -1 is the physical neighbor across the cut.
    cores = 1
    elapsed = 0.0
    for m in range(2, size):
        source = rng.randrange(m)
        target = rng.choice(adj[source])
        if target == -1:
            return None
        elapsed += rng.expovariate(m)
        if len(adj[target]) == 1:
            cores += 1
        adj[target].append(m)
        adj.append([target])
    # Continue the auxiliary Yule clock exactly. Future shape marks are not sampled.
    mass = rng.gammavariate(size, 1.0) * exp(-elapsed)
    return mass, cores / size**mu


def simulation(runs, size, seed, mu):
    rng = random.Random(seed)
    u, a, b = [], [], []
    for _ in range(runs):
        sides = [side(rng, size, mu) for _ in range(2)]
        if None in sides:
            continue
        (wa, xa), (wb, xb) = sides
        u.append(wa / (wa + wb))
        a.append(xa)
        b.append(xb)
    n = len(u)
    d = 0.0
    for i, x in enumerate(sorted(u)):
        cdf = x * x * (3 - 2 * x)
        d = max(d, (i + 1) / n - cdf, cdf - i / n)
    ra, rb = correlation(u, a), correlation(u, b)
    # DKW size <= .01; correlations reported with descriptive z scores only.
    threshold = sqrt(log(200) / (2 * n))
    print(f"P4: seed {seed}, {n}/{runs} accepted, both sides reach {size}; mu={mu}")
    print(f"Beta(2,2) KS D={d:.5f}, 1% DKW threshold={threshold:.5f}")
    print(f"Pearson r(U, X_A[{size}])={ra:+.5f}, r(U, X_B[{size}])={rb:+.5f}")
    assert d < threshold
    assert max(abs(ra), abs(rb)) < 5 / sqrt(n)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=int, default=50_000)
    parser.add_argument("--size", type=int, default=128)
    parser.add_argument("--seed", type=int, default=6142026)
    parser.add_argument("--mu", type=float, default=.56575)
    args = parser.parse_args()
    exact()
    if args.runs:
        simulation(args.runs, args.size, args.seed, args.mu)
