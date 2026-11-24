"""Finite controls for the proportional-cutoff proof; standard library only."""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial


def add(parent, at):
    return parent + (at,)


def adjacency(parent):
    adj = [[] for _ in parent]
    for v, p in enumerate(parent[1:], 1):
        adj[v].append(p)
        adj[p].append(v)
    return adj


def census(parent, high):
    adj = adjacency(parent)
    out = []
    seen = set(high)
    for v in range(len(adj)):
        if v in seen:
            continue
        group = [v]
        seen.add(v)
        for w in group:
            for u in adj[w]:
                if u not in seen:
                    seen.add(u)
                    group.append(u)
        if any(len(adj[w]) >= 2 for w in group):
            out.append((min(group), frozenset(group)))
    return out


def cut_controls():
    checked = 0
    for n in range(3, 8):
        for prefix in product(*(range(i) for i in range(1, n))):
            parent = (-1,) + prefix
            adj = adjacency(parent)
            by_degree = [frozenset(i for i in range(n) if len(adj[i]) > cap)
                         for cap in range(1, n)]
            for high in set(by_degree):
                ancestor = set()
                for h in high:
                    while h >= 0:
                        ancestor.add(h)
                        h = parent[h]
                descendants = [{v} for v in range(n)]
                for v in range(n-1, 0, -1):
                    descendants[parent[v]].update(descendants[v])
                groups = census(parent, high)
                singles = {v for v, g in groups if len(g) == 1}
                assert singles == {v for v in range(n) if len(adj[v]) >= 2 and
                                   v not in high and set(adj[v]) <= high}
                for v, g in groups:
                    if v != 0:
                        assert parent[v] in high
                    if v != 0 and v not in ancestor:
                        assert g == descendants[v]
                        assert len(g) >= 2
                    if len(g) == 1:
                        assert any(v in (set(adj[h]) & set(adj[j]))
                                   for h,j in combinations(high, 2))
                checked += 1
    # Original label 1 may be a bounded core with no orientation parent.
    parent = (-1,0,0,1,2,1,2)
    d = list(map(len,adjacency(parent)))
    high = {v for v,x in enumerate(d) if x>2}
    assert high == {1,2} and census(parent,high) == [(0,frozenset({0}))]
    assert census(add(parent,0),high) == [(0,frozenset({0,7}))]
    print(checked,'exact original-degree masks; orientation-root singleton verified')


def poisson_control():
    # Conditional on six proposal points at one hub, a deterministic-window
    # mark takes exactly one retained type or the residual category. The
    # factorial cancellation below is what produces independent Poisson cells.
    p = (Q(2,9),Q(1,9),Q(4,9),Q(2,9))
    h = Q(3)
    checked = 0
    for n in range(7):
        dist = Counter()
        for marks in product(range(len(p)), repeat=n):
            w = Q(1)
            for m in marks:
                w *= p[m]
            dist[tuple(marks.count(i) for i in range(len(p)))] += w
        assert sum(dist.values()) == 1
        for ks, prob in dist.items():
            rhs = Q(1)
            for pi,k in zip(p,ks):
                rhs *= (h*pi)**k/factorial(k)
            assert h**n/factorial(n)*prob == rhs
            checked += 1
    # Integrate the envelope for finite-size occupancy by the gamma integrals
    # ∫exp(-λu)(λu)^j/j! du = 1/λ. Include immigration and occupancy jumps.
    for s in range(2,8):
        area = Q(2*(s-1),3)
        assert sum((Q(2,3) for _ in range(s-1)),Q(0)) == area
        for k in (1,2,5):
            eps = Q(1,100)
            assert 2*eps*k*(1+s*area) == Q(k*(3+2*s*(s-1)),150)
    print(checked,'exact marked Poisson factorial cancellations; 18 clock-shift budgets')


def split_control():
    # Keep one original cut slot at the root of each independent side.
    side = adjacency((-1,0,1,1))
    degree = [len(ns)+(v==0) for v,ns in enumerate(side)]
    for q in range(len(side),len(side)+19):
        rates = [sum((Q(1,degree[w]) for w in ns),Q(0)) for ns in side]
        users = [Q(0)]*len(side)
        for w,ends in enumerate(side):
            for v in ends:
                users[v] += q*Q(1,q)*Q(1,degree[w])
        assert rates == users and sum(users) == len(side)-Q(1,degree[0])
        assert q-sum(users) == q-len(side)+Q(1,degree[0])
    # In the inert degree/size ratio, a fixed embedded trajectory has a single
    # crossing at any possible positive level, as the gamma split varies.
    for a,b,d in product((Q(1,4), Q(1,2), Q(3,4)), repeat=3):
        if a>b:
            continue
        for z in (Q(1,10),Q(1,4),Q(1,2)):
            if a<=z*b:
                continue
            w = z*d/(a-z*b+z*d)
            ratio = lambda y: a*y/(b*y+d*(1-y))
            assert ratio(w/2)< z == ratio(w) < ratio((1+w)/2)
    print('19 exact inert Yule slot generators and strict Dirichlet ratio crossings')


if __name__ == '__main__':
    cut_controls()
    poisson_control()
    split_control()
