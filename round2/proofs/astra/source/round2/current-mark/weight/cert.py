#!/usr/bin/env python3
"""Exact edge-coefficient certificate for a current-mark capped gauge."""
from fractions import Fraction as Q
from math import ceil
import argparse

RHO=Q(141,250)
DAMP=Q(8,25)
GRID=10**10


def grids(cap):
    # Endpoints cover the true, uncapped reciprocals in the final chamber.
    for j in range(2,cap):
        yield j,Q(1,j)
    yield cap,Q(1,cap)
    yield cap,Q(0)


def basic(cap, d):
    z=lambda D: 1/(1+d*min(D,cap))
    k=lambda D: (z(D)-z(D+1))/d
    return z,k


def propose(cap,d,rho, tight=False):
    z,k=basic(cap,float(d)); a=.5+float(d)
    Z=[z(D) for D in range(cap+2)]
    js=[(j,float(p)) for j,p in grids(cap)]
    if not tight:
        Z[cap]=max([z(cap)+a/ cap * k(j)/z(j) for j,p in js]); Z[cap+1]=Z[cap]
    for _ in range(100):
        old=Z[:]
        if tight:
            Z[cap]=max(z(cap),*(z(cap)+(a*k(j)/(a*old[j]+(old[j]-old[j+1])/cap))/cap for j,p in js))
            Z[cap+1]=Z[cap]
        for D in range(cap-1,1,-1):
            Z[D]=max(z(D),Z[D+1],max(
                (a*(z(D)+r*k(D)+(1/D)*k(j)/old[j])+r*Z[D+1])
                /(a+r+(1/D)*(old[j]-old[j+1])/old[j])
                for j,p in js for r in (p,D-2+p)))
        if max(abs(x-y) for x,y in zip(Z,old))<1e-13:
            break
    # Uniform relative inflation makes all edge defects strictly positive:
    # each left product is quadratic while the uninflated RHS is at most linear.
    Z=[Q(ceil(x*1.000005*GRID),GRID) for x in Z]
    Z[cap+1]=Z[cap]
    return Z


def prove(cap,d,rho,Z,plant=None):
    z,k=basic(cap,d)
    A=Q(1,2); g=rho-d
    need=lambda p,msg: None if p else (_ for _ in ()).throw(ValueError(msg))
    need(all(1>=Z[D]>=max(Z[D+1],z(D)) for D in range(2,cap+1)), 'edge factors')
    m=None; how=0
    for D in range(2,cap+1):
        for pd in ((Q(1,D),) if D<cap else (Q(0),Q(1,cap))):
            for j,pj in grids(cap):
                for r in ((pj,D-2+pj) if D<cap else (0,)):
                    for a in (Q(1,cap),A):
                        u=a+d
                        # Negative rate terms are source root gain and its
                        # true send into the first child, respectively.
                        ideal=z(D) if plant!='cap20-rates' else 1/(1+d*min(D,20))
                        margin=(u*(Z[D]-ideal-r*k(D))*Z[j]
                                +r*(Z[D]-Z[D+1])*Z[j]
                                +pd*(Z[j]-Z[j+1])*Z[D]-u*pd*k(j))
                        if plant=='wrong-send':
                            margin -= u*(Q(1,2)-pd)*k(j)
                        need(margin>=0,f'edge corner {(D,pd,j,pj,r,a)}: {margin}')
                        m=margin if m is None else min(m,margin)
                        how+=1
        for r in ((Q(0),Q(D)-Q(3,2)) if D<cap else (Q(0),)):
            for a in (Q(1,cap),A):
                margin=(a+d)*(Z[D]-z(D)-r*k(D))+r*(Z[D]-Z[D+1])
                need(margin>=0,f'first-edge corner {(D,r,a)}: {margin}')
                m=min(m,margin)
                how+=1
    # Tail h may be smaller than Z/g on a standalone root, but the uniform
    # internal profile below keeps the same setting on all descendants.
    H={(cap,c):Z[cap]/g for c in range(1,cap+1)}
    def h(D,c): return H[min(cap,D),min(cap,c)]
    for D in range(cap-1,1,-1):
        for c in range(D,0,-1):
            L=D-c; b=Q(L,D)
            B=b*(h(D,c+1)+Z[D]) if L else 0
            p=max(Q(1,50),*((r*h(D+1,c)+B)/(g+r+b) for r in (L,L+Q(c,2)+A+d)))
            H[D,c]=Q(ceil(p*GRID),GRID)
    virgin=h(2,1)
    need(virgin<=1,f'birth normalization needs {float(virgin):.10f}>1')
    H[2,1]=Q(1)
    for D in range(2,cap):
        for c in range(1,D+1):
            L=D-c;b=Q(L,D)
            B=b*(h(D,c+1)-h(D,c)+Z[D]) if L else 0
            for r in (L,L+Q(c,2)+A+d):
                need(r*(h(D+1,c)-h(D,c))+B<=g*h(D,c),f'local {(D,c,r)}')
    return H,how,m,virgin


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cap',type=int,default=36)
    p.add_argument('--plant',choices=('wrong-send','cap20-rates'))
    p.add_argument('--loose-tail',action='store_true')
    p.add_argument('--rho',type=Q,default=RHO)
    p.add_argument('--d',type=Q,default=DAMP)
    args=p.parse_args()
    Z=propose(args.cap,args.d,args.rho,not args.loose_tail)
    H,n,m,v=prove(args.cap,args.d,args.rho,Z,args.plant)
    print(f'PASS cap={args.cap} rho={args.rho} d={args.d} edge corners={n} degree rows={len(H)}')
    print(f'min edge slack={m} virgin required={v} h_min={min(H.values())} h_max={max(H.values())}')


if __name__=='__main__':
    main()
