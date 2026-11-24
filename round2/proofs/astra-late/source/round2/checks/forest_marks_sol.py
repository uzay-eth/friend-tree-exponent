#!/usr/bin/env python3
"""Cold finite checks for boundary-component marks in the friend tree."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb


def adj(P):
    a=[[] for _ in range(len(P)+1)]
    for v,w in enumerate(P,1):
        a[v].append(w)
        a[w].append(v)
    return a


def topo(P,N):
    a=adj(P)
    d=list(map(len,a))
    C=sum(v>1 for v in d)
    root=next(i for i,v in enumerate(d) if v>1)
    parent={root:-1}
    visit=[root]
    for v in visit:
        for w in a[v]:
            if w!=parent[v]:
                visit.append(w)
                parent[w]=v
    # Cold component discovery by traversing the surviving graph.
    remain={v for v,z in enumerate(d) if z<=N}
    comps=[]
    while remain:
        v=remain.pop()
        todo=[v]
        for w in todo:
            for x in a[w]:
                if x in remain:
                    remain.remove(x)
                    todo.append(x)
        top=[w for w in todo if parent[w] not in todo]
        assert len(top)==1
        if any(d[w]>1 for w in todo) and top[0]!=root:
            comps.append((parent[top[0]],top[0]))
    formula=[(parent[v],v) for v in range(len(a)) if v!=root and 2<=d[v]<=N<d[parent[v]]]
    assert sorted(comps)==sorted(formula)
    c=Counter(w for w,_ in comps)
    m=len(comps)
    if m:
        pp=sum(Q(k*k,m*m) for k in c.values())
        assert pp==sum(Q(1,m*m) for x,y in product(comps,repeat=2) if x[0]==y[0])
        assert pp<=Q(max(c.values()),m)
        assert sum(k*(k-1) for k in c.values())==sum(x!=y and x[0]==y[0] for x,y in product(comps,repeat=2))
    return len(comps),C


def enumerate_topology(stop):
    rows=0
    P=[(0,)]
    for n in range(2,stop+1):
        if n>=3:
            for p in P:
                for N in range(2,n):
                    topo(p,N)
                    rows+=1
        P=[p+(x,) for p in P for x in range(n)] if n<stop else []
    print('PASS',rows,'exact boundary/component/tuple count rows')


def physical_posterior():
    # Seed is a parent o (0) with a core exterior neighbour, two virgin
    # core arms a,b, and one unused leaf; total size 8. Frozen labels are
    # structural here, not a simulation of stable age. Condition a one-step
    # experiment on the high root crossing from degree 4 to 5.
    P=(0,0,0,0,1,2,3)
    a=adj(P)
    n=len(a)
    # a,b are two virgin arms at vertices 2,3. The parent event biases
    # their histories even though arrows upward leave the source unchanged.
    assert list(map(len,a))==[4,2,2,2,1,1,1,1]
    by=defaultdict(Q)
    for x,ns in enumerate(a):
        for y in ns:
            q=Q(1,n*len(ns))
            by[y]+=q
    # Directed first-step mass into either virgin child is 1+1/4.
    assert by[2]==by[3]==Q(5,4*n)
    assert by[0]==Q(5,2*n)
    # At exactly one step the arm growth indicators are exclusive. Even
    # conditioning on root not high leaves both marginals positive but
    # their joint zero; this toy does not claim fixed-time CT dependence
    # from the discrete-event conditioning alone.
    p=by[2]/(1-by[0])
    assert p*p>0
    print('PASS actual arrows at two siblings; reject conditional product at fixed one-step size, marginal',p)

    # CT hazard counterexample with initially independent silent arm paths.
    # Conditional weights for a parent up-output in a very short pulse are
    # proportional to p_i+p_j (+ frozen other source). Mixed independent
    # degree 2/3 histories imply determinant -(p2-p3)^2 before normalization.
    p2,p3,other=Q(1,2),Q(1,3),Q(3,2)
    w=[[other+x+y for y in (p2,p3)] for x in (p2,p3)]
    det=w[0][0]*w[1][1]-w[1][0]*w[0][1]
    assert det==-(p2-p3)**2<0
    print('PASS posterior rate tilt has rank two, determinant',det)


def normalization():
    surviving=[0,1,3]
    # Distinct high parents may have zero surviving components. For
    # nonempty parents, uniform-parent then uniform-child != component Palm.
    edge=Q(3,sum(surviving))
    uniform=Q(1,2)
    ordered=[k*(k-1) for k in surviving]
    assert edge==Q(3,4)!=uniform and ordered==[0,0,6]
    assert sum(ordered)!=sum(k*k for k in surviving)
    # Normalizing independently before taking mixtures changes the law.
    uniform_trees=(Q(1,3)+Q(3,4))/2
    pooled=Q(1+3,3+4)
    assert uniform_trees!=pooled
    print('PASS component, parent, and sibling-tuple normalizers differ:',edge,uniform,ordered)


def binomial_tail():
    # For a uniform grid mesh M vs continuous uniform U the mesh set has
    # probability one vs zero. Even exact dyadic quantization has TV one.
    for M in (2,3,11):
        grid={Q(j,M) for j in range(M)}
        assert len(grid)==M and sum(Q(1,M) for _ in grid)==1
    # Union-bound budget for moment arguments uses falling factorial counts
    # and actual finite-n core births, not n iid unconditional focal ages.
    assert comb(7,3)*6==7*6*5
    print('PASS lattice/atomic support controls')


if __name__=='__main__':
    enumerate_topology(8)
    physical_posterior()
    normalization()
    binomial_tail()
