#!/usr/bin/env python3
"""Backward all-age branching registration of negative physical path rows."""
import argparse
from fractions import Fraction as Q
from itertools import combinations, product
from screen import HI, LO, PATH, RECEIPT, Rows, change

GRID=2**20
BRANCHES={
    'fork': ((0,0),(1,0,2)),
    # In a bent P4 distinguish whether the extended arm is older or younger
    # and whether the grandchild arrived before the other child of the apex.
    'bent-older': ((0,0,1),(3,1,0,2)),
    'bent-younger': ((0,0,2),(1,0,2,3)),
    'bent-early': ((0,1,0),(2,1,0,3)),
}


def solve_row(a, b, boxes, den, num, successors, upper=False):
    """Directional floor, preserving the common reward/gain exterior hazard."""
    nz=[i for i,(x,y) in enumerate(boxes) if x!=y]
    corners=list(product(*(boxes[i] for i in nz)))
    # Piecewise affine numerator on the two halfboxes separated by a+b*rho=0.
    # Add vertices on the cut hyperplane (each has at most one free coordinate).
    for k,i in enumerate(nz):
        if not b[i]:
            continue
        for v in product(*(boxes[j] for j in nz if j!=i)):
            x=-(a+sum(b[j]*t for j,t in zip((j for j in nz if j!=i),v)))/b[i]
            if boxes[i][0]<x<boxes[i][1]:
                corners.append(v[:k]+(x,)+v[k:])
    inequalities=[(num+sum(t*successors[i] for i,t in zip(nz,v))+
                   GRID*max(0,-a-sum(b[i]*t for i,t in zip(nz,v))),den+sum(v)) for v in corners]
    if upper:
        value=max(-(-x//z) for x,z in inequalities)
        assert value>=0 and all(value*z>=x for x,z in inequalities)
    else:
        value=min(x//z for x,z in inequalities)
        assert value>=0 and all(value*z<=x for x,z in inequalities)
    return value,len(corners)


def price(r,caps,projected=False,extra=0,shape='chain',upper=False):
    # Every existing prefix continues through an allowed promotion. The same
    # event may spawn one new longer prefix, with birth order fixed by shape.
    nxt={}
    n=len(caps)
    parents,order=BRANCHES.get(shape,(tuple(range(n-1)),tuple(range(n))))
    for k in range(n,0,-1):
        deg=[sum(x==i for x in parents[:k-1])+(i>0) for i in range(k)]
        minimum=tuple(x+(i==0) for i,x in enumerate(deg))
        table={}
        corners=0
        positive=0
        best=Q(0)
        for D in product(*(range(cap,max(2,c)-1,-1) for cap,c in zip(caps[:k],minimum))):
            for c in product(*(range(min(d,m+extra),m-1,-1) for d,m in zip(D,minimum))):
                p=tuple((D[i],c[i]) for i in order) if k==n else ()
                if k==n and not r.s.support(PATH[n],p[1:-1]):
                    continue
                boxes=tuple((Q(0),Q(c[i]-deg[i],2)) for i in range(k))
                if k==n:
                    A,B,E,_=r.row(p)
                    a=A+(E if projected else 0)
                    b=tuple(B[order.index(i)] for i in range(n))
                    best=max(best,-a-sum(min(x*y for y in box) for x,box in zip(b,boxes)))
                    positive+=a+sum(max(x*y for y in box) for x,box in zip(b,boxes))<0
                else:
                    a,b=0,(0,)*k
                nb=Q(0)
                if k<n:
                    i=parents[k-1]
                    nb=Q(D[i]-c[i],D[i])*nxt.get((D+(2,),change(c,i,c[i]+1)+(1,)),0)
                y=[table.get((change(D,i,d+1),c),0) for i,d in enumerate(D)]
                base_d,base_n=(LO if upper else HI),nb
                for i,d in enumerate(D):
                    ij=[j for j in range(1,k) if parents[j-1]==i]
                    if i:
                        ij.append(parents[i-1])
                    intern=sum(Q(1,D[j]) for j in ij)
                    L=d-c[i]
                    hazard=L+intern
                    birth=Q(L,d)
                    base_d+=hazard+birth
                    base_n+=hazard*y[i]+birth*table.get((D,change(c,i,c[i]+1)),0)
                v, checks=solve_row(a,b,boxes,base_d,base_n,y,upper)
                table[D,c]=v
                corners+=checks
        print(r.side,'projected' if projected else 'physical','prefix',k,'states',len(table),
              'corners',corners,'robust positive',positive,'max ticks',float(best),
              'minimum-label D=2 value ticks',float(Q(table.get(((2,)*k,minimum),0),GRID)),flush=True)
        nxt=table
    return (HI if upper else LO)*Q(nxt[(2,),(1,)],GRID)


def check_registration():
    # Chronological labels: parent 0 is the founder. Generate by extending
    # prefixes at the prescribed parent and compare to physical vertex sets.
    true_parent=(None,0,1,1,2,3,2,4,1,8,8,9,4)
    templates={**BRANCHES, 'chain3':((0,1),(0,1,2)), 'chain4':((0,1,2),(0,1,2,3))}
    got={key:[set() for _ in range(len(order))] for key,(pa,order) in templates.items()}
    counts={key:0 for key in templates}
    for j in range(1,len(true_parent)):
        for key,(pa,_) in templates.items():
            rows=got[key]
            for k in range(len(pa),0,-1):
                for prefix in rows[k-1]:
                    if prefix[pa[k-1]]==true_parent[j]:
                        rows[k].add((*prefix,j))
            rows[0].add((j,))
    for n in (3,4):
        seen=set()
        for vs in combinations(range(1,len(true_parent)),n):
            idx={x:i for i,x in enumerate(vs)}
            if any(true_parent[x] not in idx for x in vs[1:]):
                continue
            parent=tuple(idx[true_parent[x]] for x in vs[1:])
            find=[key for key,(pa,_) in templates.items() if pa==parent]
            if not find:
                assert n==4 and parent in ((0,0,0),(0,1,1))  # S3 centered at first or second birth
                continue
            key,=find
            assert vs in got[key][-1] and vs not in seen
            seen.add(vs)
            counts[key]+=1
        assert len(seen)==sum(len(got[key][-1]) for key,(pa,_) in templates.items() if len(pa)==n-1)
    assert all(counts.values())
    print('birth registration oracle',counts,'all nonfounder embeddings unique')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--caps',type=int,nargs='+')
    p.add_argument('--check-registration',action='store_true')
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--projected',action='store_true')
    p.add_argument('--upper-control',action='store_true',help='ceiling for the same registered class')
    p.add_argument('--extra',type=int,default=0)
    p.add_argument('--shape',choices=('chain',*BRANCHES),default='chain')
    a=p.parse_args()
    if a.check_registration:
        check_registration()
        return
    if a.caps is None:
        p.error('--caps is required when pricing')
    assert len(a.caps) in (3,4) and min(a.caps)>=2
    assert a.extra>=0 and (a.shape=='chain' or len(a.caps)==len(BRANCHES[a.shape][1]))
    for side in a.side or ('lower','upper'):
        r=Rows(side)
        v=price(r,a.caps,a.projected,a.extra,a.shape,a.upper_control)
        receipt=RECEIPT[side]
        print(side,a.shape,a.caps,'extra',a.extra,'ceiling' if a.upper_control else 'floor',
              'reward ticks/core',float(v),'fraction',v,
              'potential units/core',float(v/r.s.h['S']),
              'ticks divided by old maximum comparison',float(v/receipt['m_max']),flush=True)


if __name__=='__main__':
    main()
