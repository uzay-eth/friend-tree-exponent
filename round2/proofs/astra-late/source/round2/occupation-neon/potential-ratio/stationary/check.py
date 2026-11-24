#!/usr/bin/env python3
"""Exact stationary upper separator for sign*(Phi_hi-q*Phi_lo)."""
import argparse
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from check import Cert, PINS, pin, reductions, local, differences

LOW,HIGH=Q('.565673'),Q('.565775')
GRID=10**15


class Envelope:
    def __init__(self,lo,hi,q,sign):
        self.c,self.scale=lo,q.denominator*lo.S
        self.q,self.sign=q,sign
        r=reductions(lo,hi,q.numerator,q.denominator)
        if sign<0:
            r=(-r[0], *[(-hi,-lo) for lo,hi in r[1:]])
        self.r=r
        self.delta=differences(q.numerator,q.denominator)
        self.gmat=sign*self.delta(lo.g,hi.g)
        self.base=local(lo,r,1)

    def parts(self,d,k):
        c=self.c
        f,g,hv,p4,yb,hc,yt,s3=self.r
        l=d-k-1
        s=(l,k,1)
        out=[Q(int(f[c.vi(*s)]))]+[Q(0)]*7
        for x,w in c.embed(s):
            j=c.child[x]
            out[1]+=w*int(g[1][j])
            for i,tab in ((2,hv),(4,p4),(6,yb)):
                out[i]+=w*max(0,int(tab[1][j]))
            if x in c.hc and 2<=k<=c.K:
                out[3]+=w*comb(k,2)*int(hc[1][c.hc[x]])
        for i,ind,x,m in ((5,c.yt,yt,2),(7,c.s3,s3,3)):
            if s in ind:
                out[i]+=comb(k,m)*int(x[1][ind[s]])
        assert sum(out)<=self.base[min(k,c.KF+1),1]
        return tuple(v/self.scale for v in out)

    def tail(self,t):
        assert t<=0  # all larger arities in the overflow bucket
        return max(Q(v,self.scale)+t*(k-1) for (k,h),v in self.base.items() if h==1)

    def edge_gap(self,p,x):
        """Raw shared g cap minus the exact physical paired edge."""
        c=self.c
        e=c.embed(x)
        v=sum(w*int(self.r[1][1][c.child[a]]) for a,w in e)
        assert p[1]>=1
        w=sum(pa*xb*int(self.gmat[c.state[a],c.state[b]])
              for a,pa in c.embed(p) for b,xb in e)
        assert v>=w
        return (v-w)/self.scale


def edge_credit(e,lo,hi,cap,full):
    """Stopped physical newborn solve for disjoint known allocation gaps."""
    from functools import cache
    pth=HERE.parents[1]/'check'
    sys.path.insert(0,str(pth))
    from crowns import J,gain,moves,n
    sys.path.insert(0,str(HERE.parents[1]/'ratio-transfer'))
    from replay import Pair, embedding, shape, tree
    from extract import Cert as PhysicalCert
    pair=Pair(PhysicalCert(lo),PhysicalCert(hi)) if full else None

    @cache
    def physical(z):
        return tree(z[0], *(physical(y) for y in z[1]))

    @cache
    def parts(x):
        return e.parts(sum(x),x[1])

    def owned_gap(z):
        if not full:
            l,k=z[0],len(z[1])
            return sum(e.edge_gap((l,k,1),(y[0],len(y[1]),1)) for y in z[1])
        p,x=shape(physical(z),1)
        # A: branching groups owned by this tracked root. B: bottom groups
        # with their ancestor exactly at this root. This assigns each gap to
        # just one ancestor; no visited fringe counts it a second time.
        dep=[0]
        for parent in p[1:]:
            dep.append(dep[parent]+1)
        depth={1:1,2:2,4:3,5:1,6:2}
        gap={(0,i):parts(x[0])[i] for i in (3,7)}
        for j,xj in enumerate(x):
            for i,d in depth.items():
                if dep[j]==d:
                    gap[j,i]=parts(xj)[i]
        kinds=('f','g','hV','hC','P4','Yt','Yb','S3')
        for kind,ids in embedding(p):
            i=kinds.index(kind)
            if i==0:
                continue
            j=ids[1] if i==5 else ids[-1] if i in depth else ids[0]
            if (j,i) in gap:
                a,b=pair.weight(kind,tuple(x[k] for k in ids))
                gap[j,i]-=e.sign*(b-e.q*a)
        assert all(v>=0 for v in gap.values()),(z,gap)
        return sum(gap.values())

    @cache
    def visit(z):
        if n(z)>cap:
            return Q(0)
        reward=owned_gap(z)
        jumps=moves(z)
        rate=sum(r for _,r in jumps)
        value=reward+sum(r*visit(y) for y,r in jumps)
        rec=visit(gain(z))
        v=min((value+s*rec)/(HIGH+rate+s) for s in (Q(0),Q(1,2)))
        floor=Q(v*GRID//1,GRID)
        assert all((HIGH+rate+s)*floor<=value+s*rec for s in (Q(0),Q(1,2)))
        return floor

    return LOW*visit(J),visit.cache_info().currsize


def bellman(p,cap,upper):
    """Floats only propose a recharge; final arithmetic uses outward rationals."""
    cast=Q if upper else float
    lo,hi=cast(LOW),cast(HIGH)
    stop=cast(p[-1][0])
    nxt=[stop]*(cap+1)
    for d in range(cap,1,-1):
        cur=[cast(0)]*d
        for k in range(d-1,-1,-1):
            l=d-k-1
            b=cast(l)/d
            val=cast(p[d][k])+(b*cur[k+1] if l else 0)
            x=max((val+g*nxt[k])/(s+b+g) for s in (lo,hi) for g in (cast(l),cast(l)+cast(k+1)/2))
            cur[k]=Q(-(-x*GRID//1),GRID) if upper else x
            if upper:
                assert all((s+b+g)*cur[k]>=val+g*nxt[k]
                           for s in (lo,hi) for g in (Q(l),Q(l)+Q(k+1,2)))
        nxt=cur
    return cur[0]


def run(lo,hi,paths,q,sign,cap,edge_cap,full):
    e=Envelope(lo,hi,q,sign)
    bare=[None,None]+[[sum(e.parts(d,k)) for k in range(d)] for d in range(2,cap+1)]

    def prices(t):
        tail=e.tail(t)
        return [None,None]+[[v+t*(k-1) for k,v in enumerate(bare[d])] for d in range(2,cap+1)]+[[max(tail/LOW,tail/HIGH)]]

    t=min((Q(-i,1000) for i in range(251)),key=lambda t:bellman(prices(t),cap,False))
    v=bellman(prices(t),cap,True)
    credit,n=edge_credit(e,*paths,edge_cap,full) if edge_cap else (Q(0),0)
    bound=max(LOW*v,HIGH*v)-credit
    print('ratio',q,'sign',sign,'cap',cap,'recharge',t,'all-degree cap',float(e.tail(t)),
          'lifetime upper',float(v),'motif gap credit/states',float(credit),n,
          'per-core outward upper',float(bound),'PASS' if bound<0 else 'unresolved',flush=True)
    if bound<0:
        assert bound<=Q(-1,10**7)
    return bound


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lo',type=Path)
    p.add_argument('hi',type=Path)
    p.add_argument('trials',nargs='+',type=Q)
    p.add_argument('--cap',type=int,default=40)
    p.add_argument('--edge-cap',type=int,default=9)
    p.add_argument('--full-credit',action='store_true')
    p.add_argument('--preauthenticated',action='store_true')
    a=p.parse_args()
    paths=(a.lo,a.hi)
    def stamps():
        return tuple((s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns) for s in (x.stat() for x in paths))
    stamp=stamps()
    if not a.preauthenticated:
        assert tuple(map(pin,paths))==PINS
        print('PASS both SHA256 pins',flush=True)
    lo,hi=map(Cert,paths)
    assert (lo.mode,hi.mode)==(2,1) and (lo.lam,hi.lam)==(LOW,HIGH) and lo.S==10**11
    for key in ('S','D','K','LF','KF','state','par','child','hv','hc','p4','yt','ybu','ybw','s3'):
        assert getattr(lo,key)==getattr(hi,key)
    for q in a.trials:
        run(lo,hi,paths,q,1 if q>=1 else -1,a.cap,a.edge_cap,a.full_credit)
    assert stamps()==stamp
