#!/usr/bin/env python3
"""Exact physical-tree check of L4->L5 source and promotion ownership."""
from fractions import Fraction as Q
from itertools import combinations
from random import Random


def motifs(adj, size):
    return [vs for vs in combinations(range(len(adj)), size)
            if sum(y in adj[x] for x, y in combinations(vs, 2)) == size - 1]


def edges(parents):
    adj = [set() for _ in range(len(parents)+1)]
    for y, x in enumerate(parents, 1):
        adj[x].add(y)
        adj[y].add(x)
    return adj


def q4(m, degree):
    # Arbitrary signed, nonconstant function; degree jumps remain visible.
    return Q(sum((u+3)*degree[u]**2 - (u+1)*degree[u]**3 for u in m), 739)


rng = Random(155)
trees = gains = transfers = births = cutoff_exits = 0
for size in range(5, 14):
    for rep in range(25):
        parents = [rng.randrange(v) for v in range(1, size)]
        adj = edges(parents)
        deg = [max(2, len(adj[u])) + rng.randrange(3) for u in range(size)]
        for d5 in (6, 7, 8):
            eligible = lambda m: all(deg[u] < d5 for u in m)
            qgain = lambda m, u: q4(m, deg[:u]+[deg[u]+1]+deg[u+1:]) - q4(m, deg)
            ms = motifs(adj, 4)
            ns = motifs(adj, 5)
            direct = kept = moved = Q(0)
            for m in ms:
                for u in m:
                    for z in adj[u].difference(m):
                        a = Q(1, deg[z])
                        b = min(a, Q(1, d5)) if eligible(m) else a
                        delta = qgain(m, u)
                        direct += a*delta
                        kept += b*delta
                        gains += 1
            for n in ns:
                if not eligible(n):
                    continue
                for z in n:
                    near = adj[z].intersection(n)
                    if len(near) == 1:
                        (u,) = near
                        m = tuple(v for v in n if v != z)
                        moved += (Q(1, deg[z])-Q(1, d5))*qgain(m, u)
                        cutoff_exits += deg[u]+1 == d5
                        transfers += 1
            assert direct == kept + moved
            for u in range(size):
                if not eligible([u]) or deg[u] == len(adj[u]):
                    continue
                post = edges(parents+[u])
                direct_births = {n for n in motifs(post, 5) if size in n and eligible(n[:-1])}
                transferred_births = {m+(size,) for m in ms if u in m and eligible(m)}
                assert direct_births == transferred_births
                births += len(direct_births)
            trees += 1

assert transfers and cutoff_exits and births
print(f'Transfer checks PASS: {trees} tree/cutoff pairs, {gains} external credits, '
      f'{transfers} exact transfers, {cutoff_exits} cutoff-crossing gains, {births} newborn motifs.')
