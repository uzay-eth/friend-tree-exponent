#!/usr/bin/env python3
"""Lead's own exact test: does the 'alternating matching' order certify FULL stochastic domination of Z_u over Z_(u+1)?

Order on rooted trees of equal size (A, a) >= (B, b):
   GE0(alpha, beta)  (even level): the children of beta inject into the children of alpha, beta' -> phi(beta'), with
                     LE1(phi(beta'), beta') for each;
   LE1(alpha, beta)  (odd level):  the children of alpha inject into the children of beta, alpha' -> psi(alpha'), with
                     GE0(alpha', psi(alpha')) for each.
Coupling lemma (proved in the note next to this file): if (A, a) >= (B, b) then the two friend-tree chains can be
coupled so that the relation holds at all later sizes; in particular deg(a) >= deg(b) for ever, so Z_a >= Z_b.

Certificate for labels u, u+1 at size m: a transport plan between q_u = law of (T_m, u) and q_(u+1) = law of
(T_m, u+1) supported on comparable pairs.  It exists iff the max-flow of the bipartite network with capacities q_u(s),
q_(u+1)(t) and arcs s -> t for s >= t equals 1.  All capacities are exact integers (common denominator).
Usage: domination_transport.py [M]     (M = tree size, default 8; all histories are enumerated)
"""
import sys
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd
M = int(sys.argv[1]) if len(sys.argv) > 1 else 8

def canon(adj, v, parent=0):
    return tuple(sorted(canon(adj, w, v) for w in adj[v - 1] if w != parent))

def has_injection(small, big, ok):
    """is there an injection of the list small into the list big with ok(s, b) on every pair (bipartite matching)"""
    if len(small) > len(big): return False
    match = [-1] * len(big)
    def augment(i, seen):
        for j in range(len(big)):
            if j not in seen and ok(small[i], big[j]):
                seen.add(j)
                if match[j] < 0 or augment(match[j], seen): match[j] = i; return True
        return False
    return all(augment(i, set()) for i in range(len(small)))
@lru_cache(None)
def GE0(a, b):          # a on the A side, b on the B side, even level: children(b) inject into children(a) with LE1
    return has_injection(list(b), list(a), lambda bb, aa: LE1(aa, bb))
@lru_cache(None)
def LE1(a, b):          # odd level: children(a) inject into children(b) with GE0
    return has_injection(list(a), list(b), lambda aa, bb: GE0(aa, bb))

# all histories with exact probabilities
level = [(Q(1), ((2,), (1,)))]
for n in range(2, M):
    nxt = []
    for p, adj in level:
        deg = [len(x) for x in adj]
        for w in range(1, n + 1):
            pw = sum(Q(1, deg[u - 1]) for u in adj[w - 1]) / n
            new = [list(x) for x in adj]; new[w - 1].append(n + 1); new.append([w])
            nxt.append((p * pw, tuple(tuple(x) for x in new)))
    level = nxt
assert sum(p for p, _ in level) == 1
q = [dict() for _ in range(M + 1)]
for p, adj in level:
    for u in range(1, M + 1):
        s = canon(adj, u); q[u][s] = q[u].get(s, Q(0)) + p
shapes = sorted(set().union(*[set(d) for d in q[1:]]))
print(f"size {M}: {len(level)} histories, {len(shapes)} rooted shapes")
den = 1
for u in range(1, M + 1):
    for v in q[u].values(): den = den * v.denominator // gcd(den, v.denominator)

def maxflow(sup, dem):
    """exact max-flow source -> supply shapes -> demand shapes -> sink; integer capacities; Dinic"""
    S = list(sup); T = list(dem); ns, nt = len(S), len(T); N = ns + nt + 2; src, snk = N - 2, N - 1
    head = [[] for _ in range(N)]; to = []; cap = []
    def add(a, b, c): head[a].append(len(to)); to.append(b); cap.append(c); head[b].append(len(to)); to.append(a); cap.append(0)
    INF = den * 4
    for i, s in enumerate(S): add(src, i, int(sup[s] * den))
    for j, t in enumerate(T): add(ns + j, snk, int(dem[t] * den))
    arcs = 0
    for i, s in enumerate(S):
        for j, t in enumerate(T):
            if GE0(s, t): add(i, ns + j, INF); arcs += 1
    flow = 0
    while True:
        lvl = [-1] * N; lvl[src] = 0; qu = [src]
        for x in qu:
            for e in head[x]:
                if cap[e] > 0 and lvl[to[e]] < 0: lvl[to[e]] = lvl[x] + 1; qu.append(to[e])
        if lvl[snk] < 0: break
        it = [0] * N
        def dfs(x, f):
            if x == snk: return f
            while it[x] < len(head[x]):
                e = head[x][it[x]]; y = to[e]
                if cap[e] > 0 and lvl[y] == lvl[x] + 1:
                    g = dfs(y, min(f, cap[e]))
                    if g > 0: cap[e] -= g; cap[e ^ 1] += g; return g
                it[x] += 1
            return 0
        while True:
            g = dfs(src, INF)
            if g == 0: break
            flow += g
    return Q(flow, den), arcs

sys.setrecursionlimit(100000)
for u in range(1, M):
    f, arcs = maxflow(q[u], q[u + 1])
    same = q[u] == q[u + 1]
    print(f"labels ({u},{u+1}): supports {len(q[u])} and {len(q[u+1])} shapes, {arcs} comparable pairs, max-flow = {f} = {float(f):.9f}"
          + ("   FULL STOCHASTIC DOMINATION CERTIFIED" if f == 1 else f"   not certified at this size (missing {float(1 - f):.3e})") + ("   (equal laws)" if same else ""))
