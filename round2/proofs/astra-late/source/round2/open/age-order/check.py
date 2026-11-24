#!/usr/bin/env python3
"""Exact laws for deterministic age labels; independent directed-arrow replay."""
from collections import defaultdict, Counter
from fractions import Fraction as Q
from math import comb, lcm


def adjacency(parents):
    a = [[] for _ in range(len(parents)+1)]
    for i,p in enumerate(parents,1):
        a[i].append(p)
        a[p].append(i)
    return tuple(tuple(ns) for ns in a)


def star(a,v):
    return sum(Q(1,len(a[x])) for x in a[v])


def transitions(p):
    a = adjacency(p)
    return tuple(star(a,v)/len(a) for v in range(len(a)))


def independent(p):
    # Does not call adjacency or star: rebuild by edges and enumerate arrows.
    arcs = tuple((v,u) for i,j in enumerate(p,1) for v,u in [(i,j),(j,i)])
    deg = Counter(v for v,_ in arcs)
    law = [Q(0)]*(len(p)+1)
    for v,u in arcs:
        law[u]+=Q(1,len(law)*deg[v])
    return tuple(law)


def beta_tail(n, d, x):
    return sum(comb(n-1,j)*x**j*(1-x)**(n-1-j) for j in range(d))


def promotion_cap(k):
    if k<4:
        return Q(1)
    return min(Q(1),(k+Q(1,2))/(Q(k,2)-1)**2+Q(2,k))


def enumerate_exact(N):
    weights = {(0,):1}
    den=1
    for n in range(2,N+1):
        masses=[[0]*n for _ in range(n)]
        for ps,w in weights.items():
            ds=Counter(ps)
            for i in range(n):
                masses[i][ds[i]+(i>0)]+=w
        assert all(sum(a)==den for a in masses)
        assert masses[0]==masses[1]
        tails=[[sum(a[k:]) for k in range(n)] for a in masses]
        crosses=[(u+1,u+2,k,Q(tails[u][k],den),Q(tails[u+1][k],den))
                 for u in range(n-1) for k in range(2,n) if tails[u][k]<tails[u+1][k]]
        print('n',n,'trees',len(weights),'adjacent reversals',len(crosses),flush=True)
        if crosses:
            print('  first',crosses[0],flush=True)
        if n==N:
            for i in range(n):
                print('  label',i+1,'P(D>=2)',Q(tails[i][2],den),
                      'P(D>=3)',Q(tails[i][3],den),'mean',Q(sum(k*x for k,x in enumerate(masses[i])),den))
            return
        m=lcm(*range(1,n))
        new={}
        for ps,w in weights.items():
            a=adjacency(ps)
            scaled=[sum(m//len(a[x]) for x in ns) for ns in a]
            assert sum(scaled)==n*m
            if n<=7:
                assert tuple(Q(x,n*m) for x in scaled)==transitions(ps)==independent(ps)
            for v,r in enumerate(scaled):
                new[ps+(v,)]=w*r
        weights=new
        den*=n*m


def swap(p,x,y):
    a=adjacency(p)
    n=len(a)
    first=transitions(p)
    secondx=transitions(p+(x,))
    secondy=transitions(p+(y,))
    f,g=first[x]*secondx[y],first[y]*secondy[x]
    if x!=y and y in a[x]:
        dx,dy=len(a[x]),len(a[y])
        want=(star(a,y)*Q(1,dy*(dy+1))-star(a,x)*Q(1,dx*(dx+1)))/(n*(n+1))
    else:
        want=0
    assert f-g==want
    return f,g


def fixtures():
    # Four-path middle edge: degree tied, attachment sums differ
    # on an asymmetrically extended branch. Also look for a pair where
    # the newer leaf initially has the higher gain hazard.
    frontier=[(0,)]
    best=None
    checks=0
    for n in range(2,8):
        for p in frontier:
            for x in range(n):
                for y in range(n):
                    swap(p,x,y)
                    checks+=1
            if n>=4:
                a=adjacency(p)
                u,v=n-2,n-1
                if len(a[u])==len(a[v])==1 and star(a,u)<star(a,v) and best is None:
                    best=(p,star(a,u),star(a,v))
        frontier=[p+(x,) for p in frontier for x in range(n)]
    assert best
    print('PASS',checks,'two-order likelihood identities; later leaf wins local hazard',best)
    for p in ((0,0,0),(0,0,0,1,2,3,1,2,3)):
        a=adjacency(p)
        x,y=0,1
        f,g=swap(p,x,y)
        print('adjacent fixture degrees',(len(a[x]),len(a[y])),
              'rates',(star(a,x),star(a,y)),'path order probabilities',f,g)
    # Seed edge special event: if x was born first, the next birth to x
    # is impossible under label swap; merged unlabelled two-path is real.
    assert transitions((0,0))[2]==Q(1,6)
    print('PASS younger-to-older newborn attachment and impossible reverse order')

    # Conditional limiting counterexample, not an unconditional age reversal:
    # Adjacent old u=n-1 sits at the large hub 1; final v sits at vertex 2.
    k=200
    ps=(0,)+(0,)*(k+1)+(1,)
    a=adjacency(ps)
    u,v=len(a)-2,len(a)-1
    assert len(a[1])==2 and a[u]==(0,) and a[v]==(1,)
    tags=sum(len(a[x])==1 and x!=u for x in a[0])
    exact_upper=promotion_cap(k)
    assert tags==k and exact_upper<Q(16,tags)<Q(1,8)
    print('conditional old/new hub bounds at size',len(a),Q(16,tags),'<',Q(1,10),
          '; exact upper',exact_upper,'versus newer lower 1/8; stock',tags)
    for n in range(3,12):
        assert beta_tail(n,1,Q(1,3))==Q(2,3)**(n-1)
        assert beta_tail(n,n-1,Q(1,3))==1-Q(1,3)**(n-1)
    print('PASS integer beta tail identities')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--size',type=int,default=8)
    args=p.parse_args()
    assert 4<=args.size<=10
    fixtures()
    enumerate_exact(args.size)
