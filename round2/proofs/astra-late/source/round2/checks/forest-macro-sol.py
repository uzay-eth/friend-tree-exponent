#!/usr/bin/env python3
"""Independent exact finite controls for the macroscopic forest audit."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial


def graph(parent):
    adj = [set() for _ in parent]
    for j,i in enumerate(parent[1:],1):
        adj[j].add(i)
        adj[i].add(j)
    return adj


def blocks(parent, high):
    """Union by substituting component labels along retained edges."""
    labels = list(range(len(parent)))
    for j,i in enumerate(parent[1:],1):
        if i not in high and j not in high:
            labels[j] = labels[i]
    out = defaultdict(set)
    for i,root in enumerate(labels):
        if i not in high:
            out[root].add(i)
    return dict(out)


def chain(parent,j):
    while j>=0:
        yield j
        j = parent[j]


def cut():
    masks = singleton = edges = 0
    for n in range(2,7):
        for suffix in product(*(range(i) for i in range(1,n))):
            t = (-1,)+suffix
            adj = graph(t)
            core = {i for i,a in enumerate(adj) if len(a)>=2}
            dpaths = [set(chain(t,i)) for i in range(n)]
            # Future hub sets, after all their labels are present as cores,
            # need not be a degree threshold of this finite fixture.
            for m in range(len(core)+1):
                for hs in combinations(core,m):
                    high = set(hs)
                    cutpaths = set().union(*(dpaths[i] for i in high))
                    census = blocks(t,high)
                    for root, comp in census.items():
                        assert root==min(comp) and (root==0 or t[root] in high)
                        if root not in cutpaths:
                            assert comp=={j for j,path in enumerate(dpaths) if root in path}
                        if len(comp)==1 and root in core:
                            assert adj[root]<=high and len(adj[root])>=2
                            singleton+=1
                        for i in comp:
                            for j in adj[i]&comp:
                                assert adj[i]-comp<=high and adj[j]-comp<=high
                                edges+=1
                    for dst in range(n):
                        post = blocks(t+(dst,),high)
                        assert all(v<=post[k] for k,v in census.items())
                    masks+=1
    # 2 can be the hub and have 1 as a leaf even after the tree has grown.
    t = (-1,0,1,1,1,1,1)
    assert graph(t)[0]=={1} and len(graph(t)[1])==6
    t += (0,)  # Hub 2 sends to its inherited orientation-parent leaf.
    assert blocks(t,{1})[0]=={0,7} and t[0]==-1
    # The orientation-root pair can look like a virgin, but has no high parent.
    # After 1 promotes a child (the first action to 7), there can be a root singleton.
    t += (7,7,7,7,7)
    assert blocks(t,{1,7})[0]=={0} and len(graph(t)[0])==2
    # A non-root singleton with an original degree of two and a missing parent.
    t = (-1,0,1,2,0,0,3,3)
    assert blocks(t,{1,3})[2]=={2} and t[2]==1
    print('PASS',masks,'fixed-cut geometries;',singleton,'core singleton incidences;',edges,'retained directed edges; reversed initial promotion')


def arrows():
    # The old hubs 1 and 2 are adjacent; 1 has an old promoted child (2 in
    # zero-based labels), 2 has a leaf to promote. The new leaf arriving from
    # a return from 6 may become the next promoted child of hub 2.
    t = (-1,0,0,2,0,1)
    hb = {0,1}
    t += (5,)  # 2 -> 6 in physical labels: promote its existing leaf
    t += (1,)  # 6 -> 2: genuine return, newborn leaf is outside its branch
    t += (7,)  # 2 -> 8: promote the returning leaf with a fresh mark
    bounds = checks = 0
    for t,roots in ((t[:7],(5,)),(t[:8],(5,)),(t,(5,7)),
                    (t+(6,6,8),(5,7)),(t+(5,8,9),(5,7))):
        adj = graph(t)
        ds = tuple(len(a) for a in adj)
        owners = {}
        for root in roots:
            assert t[root] in hb and ds[root]>=2
            for j in range(len(t)):
                if root in chain(t,j):
                    assert j not in owners and j not in hb
                    owners[j] = root
        for root in roots:
            branch = {j for j,r in owners.items() if r==root}
            assert set().union(*(adj[j]-branch for j in branch))=={t[root]}
            true = Counter()
            for sender in branch:
                for receiver in adj[sender]&branch:
                    true[receiver]+=Q(1,ds[sender])
            internal = {j:sum(Q(1,ds[i]) for i in adj[j]&branch) for j in branch}
            assert true==internal and sum(true.values())==len(branch)-Q(1,ds[root])>=Q(3,2)
            checks+=len(branch)
            for dummy in range(5):
                q = len(branch)+dummy
                assert q-sum(true.values())==dummy+Q(1,ds[root])
                assert q*Q(1,q)/ds[root]==Q(1,ds[root])
                bounds+=1
        failure = sum(Q(sum(ds[i]>=2 for i in adj[h]),ds[h]) for h in hb)
        incoming = sum(Q(1,ds[t[r]]) for r in roots)
        assert failure>=incoming
    print('PASS',checks,'recipient-rate identities;',bounds,'proposal-slot checks; donor return and distinct later promotion')


def selection():
    # A finite marking control with heterogeneous sites and mark categories.
    for p in ((Q(1,3),Q(1,6),Q(1,2)), (Q(1,5),Q(3,5),Q(1,5))):
        for n in range(6):
            totals = Counter()
            for marks in product(range(3),repeat=n):
                w = Q(1)
                for m in marks:
                    w *= p[m]
                totals[tuple(marks.count(j) for j in range(3))] += w
            for k,w in totals.items():
                split = Q(factorial(n))
                for a,b in zip(p,k):
                    split *= a**b/factorial(b)
                assert w==split
    # A future-selection toy gate sharing the same marks: two independent
    # Bernoulli indicators become negatively correlated given 'at least one'.
    cols = list(product((0,1),repeat=2))
    keep = [v for v in cols if any(v)]
    assert Q(sum(a*b for a,b in keep),len(keep))==Q(1,3)!=Q(2,3)**2
    # At a deterministic shifted time, stopping when a mark changes would
    # bias its law. Control by any relevant jump on the entire bracket.
    for s in range(2,14):
        area = Q(2,3)*(s-1)
        for k in (1,4,19):
            eps = Q(1,256)
            assert 2*eps*(k+s*k*area)==Q(k*(3+2*s*(s-1)),384)
    print('PASS two independent mark splits, future-selection negative control, 36 clock bracket budgets')


def threshold():
    for v in (3,4,11,29,100):
        # Beta(1,v-1) raw moment and exact sum of the maximal certificates.
        beta2 = 1-Q(2*(v-1),v)+Q(v-1,v+1)
        assert beta2==Q(2,v*(v+1))
        n = 2*v
        assert sum(Q(2,j*(j+1)) for j in range(v,n+1))==Q(2,v)-Q(2,n+1)
    for a,b,d in product((Q(1,7),Q(2,5),Q(3,4)),repeat=3):
        # Endpoints may have different side size ratios / offsets.
        for x in (Q(1,13),Q(1,7),Q(1,2),Q(1)):
            if a<=x*b:
                continue
            w = x*d/(a-x*b+x*d)
            f = lambda z: a*z/(b*z+d*(1-z))
            assert 0<w<1 and f(w/2)<x==f(w)<f((w+1)/2)
    assert 1/(Q(2)-Q(1,2))==Q(2,3)
    print('PASS Beta maximal moment/telescoping, conditional gamma-ratio crossing, pair intensity 2/3')


if __name__ == '__main__':
    cut()
    arrows()
    selection()
    threshold()
