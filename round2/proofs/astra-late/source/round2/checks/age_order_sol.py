#!/usr/bin/env python3
"""Exact physical hazard and endpoint arithmetic controls for the Q3 partial order."""

from collections import defaultdict
from fractions import Fraction as Q
from math import factorial


def directed(edges):
    graph = defaultdict(list)
    for a,b in edges:
        graph[a].append(b)
        graph[b].append(a)
    rate = defaultdict(Q)
    for sender,targets in graph.items():
        for target in targets:
            rate[sender,target] += Q(1,len(targets))
    return graph,rate


def conditioning():
    k = 200
    n = k+4
    u,v = n-1,n
    edges = [(1,j) for j in range(2,n)]+[(2,n)]
    graph,rate = directed(edges)
    assert sum(rate.values())==n
    assert [rate[1,u],rate[2,v],rate[v,2],rate[1,2]]==[Q(1,202),Q(1,2),Q(1),Q(1,202)]
    assert sum(len(graph[j])==1 and j!=u for j in graph[1])==k
    # The specified labelled shape has one possible attachment target per step.
    mass = Q(1)
    for j in range(2,n):
        _,rates = directed(edges[:j-1])
        target = edges[j-1][0]
        mass *= sum(r for (a,b),r in rates.items() if b==target)/j
    assert mass==Q(1,203*202**2)
    old = (k+Q(1,2))/(Q(k,2)-1)**2+Q(2,k)
    young = (Q(1,2)/(Q(3,2)+rate[1,2]))/2
    assert old==Q(29851,980100)<Q(1,8)<young==Q(101,608)
    assert mass*(young-old)>0 and mass*old+1-mass > mass*young
    print('PASS n=204: event',mass,'hazards 1/202,1/2; old upper',old,'new lower',young)


def ceil(x):
    return -(-x.numerator//x.denominator)


def endpoints():
    checks = 0
    for u in (1,2,3,4,20,200):
        m,pi = (3,Q(1,2)) if u<=2 else (u+1,Q(1,u*(u-1)**2))
        B=m-1
        A=m+17+8*(B+1)
        c=pi*factorial(m-1)/Q(2*A**B)
        L=c/2**B
        K=max(16,ceil(128/L))
        assert Q(32,K)<=L/4
        # For v-1 >= T**6, the low-threshold power-bound term is <= L/4.
        T=ceil(8*(K+5)/L)
        assert Q(2*(K+5),T)<=L/4
        # v >= 2+B+ceil(log2(1/c)); integer arithmetic can avoid the log.
        p=0
        q=c
        while q<1:
            q*=2
            p+=1
        for eps in (Q(1,2),Q(1,3),Q(1,10**6),Q(1,10**30)):
            n=max(m,16,ceil(8*(B+1)/eps))
            assert n<=A/eps
            assert Q(4,n)+2*(B+1)/(eps*n)<=Q(1,2)
            # (1-x)^p <= (1/2)^p for every threshold in the upper band.
            assert Q(1,2**(p+1))<c
            checks+=1
    for k in range(8,500):
        b=(k+Q(1,2))/(Q(k,2)-1)**2+Q(2,k)
        assert b<Q(16,k)
    print('PASS',checks,'exact endpoint-budget controls, including 1-x=10^-30; stock bound 16/k')


def core_bound():
    # Enumerate all chronological trees through size seven; group their
    # directed arrows by the target's degree to compute the exact drift.
    rows=0
    frontier=[[(1,2)]]
    for n in range(2,8):
        next_frontier=[]
        for edges in frontier:
            graph,rates=directed(edges)
            deg={v:len(neighbours) for v,neighbours in graph.items()}
            C=sum(d>=2 for d in deg.values())
            H=sum(d>=3 for d in deg.values())
            if n>=4:
                a=sum(r for (_,b),r in rates.items() if deg[b]==1)
                b=sum(r for (_,b),r in rates.items() if deg[b]==2)
                assert a<=Q(C+H,2) and b<=Q(3*(C-H),2)
                assert a+Q(2,9)*b<=Q(5,6)*(C+Q(2,9)*H)
                rows+=1
            if n<7:
                next_frontier.extend(edges+[(v,n+1)] for v in graph)
        frontier=next_frontier
    print('PASS core-count drift on',rows,'chronological trees at sizes 4–7')


if __name__=='__main__':
    conditioning()
    endpoints()
    core_bound()
