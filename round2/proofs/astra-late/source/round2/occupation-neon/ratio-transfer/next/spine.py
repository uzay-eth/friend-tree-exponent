#!/usr/bin/env python3
"""Signed Poisson upper bound for one realizable D18 P4 extreme."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from replay import Pair, Cert, J, tree

T=Q('1.2415')
LO=Q('.565674722414')
HI=Q('.565773666188')
GRID=10**21
FOUR=((0,1,1),)*4
BASE=Q('.565772451624')  # first body + <=12 higher + <=4-core suffix
DEN=Q('1.518769126087491')


def up(x):
    return Q(-(-x*GRID//1),GRID)


def coefficients(pair):
    aa,ff=zip(*(e.row('P4',FOUR) for e in pair.e))
    assert all(tuple(x[:3] for x in f)==((0,True,1),(3,False,1)) for f in ff)
    scale=pair.c[0].S
    a,p,w=((x-T*y)/scale for x,y in zip((aa[0],*(i[3] for i in ff[0])),
                                         (aa[1],*(i[3] for i in ff[1]))))
    assert (a,p,w)==(Q(287589315471350942439,5*10**22),Q('-.005553926199875'),Q('-.02262119047951'))
    for d,e in product((2,7),(2,3,8)):
        crown=tree(0,tree(0,tree(0,tree(0,tree(e-1)))))
        rr=[both for f,_,both in pair.charged(crown) if f=='P4']
        assert len(rr)==1
        (x,dx),(y,dy)=rr[0]
        assert x+dx/d-T*(y+dy/d)==a+p/d+w/e
        pair.full(tree(d-1,crown),True)
    return a,p,w


def super_poisson(coeff, cap):
    a,p,w=coeff
    # Beyond the tracked degree gains and promotions leave U constant; the
    # face reward <=a and the four-vertex exit intensity is always >=3.
    tail=up(a/(LO+3))
    assert tail*(LO+3)>=a and p<0 and w<0
    U={}
    visits=0
    for d in range(cap,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1
            b=Q(l,d)
            G=U[d+1,k] if d<cap else tail
            B=U[d,k+1] if l else Q(0)
            rows=[]
            for s,r,g in product((LO,HI),(Q(0),Q(1,2)),(l+Q(1,2),l+Q(k+1,2))):
                den=s+3+Q(1,d)+r+g+b
                rhs=a+p*r+w/d+g*G+b*B
                rows.append((den,rhs))
            U[d,k]=up(max(rhs/den for den,rhs in rows))
            assert all(den*U[d,k]>=rhs for den,rhs in rows)
            visits+=len(rows)
    return U,tail,visits


def birth_path(hi):
    # Promote the successive tip five times, monitoring all current cores.
    # Then permit six gains at the tip below the bare arm. Stop on any
    # gain at the protected five or on any promotion of the tip's leaves.
    return (prod(Q(1,2)/(hi+n) for n in range(2,7))
            *prod(l/(hi+l+Q(13,2)) for l in range(1,7)))


def actual_first_events(hi):
    # Independently reconstruct monitored receiver hazards from literal
    # adjacency; outside parent has degree two for the upper denominator.
    lo=Q(1)
    for m,l,which in [(i,1,'promotion') for i in range(1,6)] + [(6,i,'gain') for i in range(1,7)]:
        adj = [set() for _ in range(m+l)]
        for u,v in [(i,i+1) for i in range(m-1)]+[(m-1,m+i) for i in range(l)]:
            adj[u].add(v); adj[v].add(u)
        deg=[len(nb)+(i==0) for i,nb in enumerate(adj)]
        hazard=[sum(Q(1,deg[j]) for j in nb)+Q(i==0,2) for i,nb in enumerate(adj)]
        assert sum(hazard)==m+l-Q(1,deg[0])+Q(1,2)
        wanted = hazard[m] if which=='promotion' else hazard[m-1]
        lo *= wanted/(hi+sum(hazard))
    return lo


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('lower',type=Path)
    ap.add_argument('upper',type=Path)
    args=ap.parse_args()
    pair=Pair(Cert(args.lower),Cert(args.upper))
    c=coefficients(pair)
    for cap in (4,8,12):
        U,tail,count=super_poisson(c,cap)
        assert U[2,0]<0 and U[2,1]<0
        print('cap',cap,'endpoint inequalities',count,'virgin and bare arm U',U[2,0],U[2,1],
              'upper tail',tail)
    # Start a bare P4 episode at crown size 13 after the selected path.
    # Supported charges are disjoint from both previous higher ledgers.
    q=birth_path(HI)
    physical_q=actual_first_events(HI)
    assert physical_q>q
    numer=-LO*physical_q*U[2,1]/T
    gain=numer/DEN
    assert gain>Q('1e-13')
    edge=BASE-gain
    print('selected discounted birth prefix, coarse and graph lower',float(q),float(physical_q))
    print('disjoint stationary upper-witness slack, lower',float(numer))
    print('strict endpoint gain vs arms ledger',float(gain))
    print('new endpoint outward at 18 places',Q(-(-edge*10**18//1),10**18))


if __name__=='__main__':
    main()
