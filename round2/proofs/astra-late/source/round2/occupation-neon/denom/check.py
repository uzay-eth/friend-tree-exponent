#!/usr/bin/env python3
"""All-degree upper charge and exact stable-age deficit from Neon D18."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb
from pathlib import Path
import sys

import numpy as np

ROUND2=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROUND2/'attacks'/'neon1'),str(ROUND2/'finite')]
from audit import Cert
import residual_entries as residual

LOW=Q('.565673570322');HIGH=Q('.565775');GRID=10**18
C0={1:195820675182,2:184124847315}


def down(x):return Q(x*GRID//1,GRID)


def show(x,places=12,upper=True):
    unit=10**places;t=-(-x*unit//1) if upper else x*unit//1
    return f'{t//unit}.{t%unit:0{places}d}'


class Charge:
    def __init__(self,c):
        self.c=c
        # Group by owner. Edge, P3, P4 and branch paths are charged to the
        # bottom vertex. Hc, Ytail and S3 are charged to their branching core.
        lim=max(max(abs(int(x.min())),abs(int(x.max()))) for x in
                (c.f,c.g,c.hvtab,c.hctab,c.p4tab,c.yttab,c.ybtab,c.s3tab))
        assert lim*(4+max(s[1] for s in c.ybu)+2*comb(c.K,2)+comb(c.K,3))<2**62
        g=c.g[np.ix_([c.state[x] for x in c.par],[c.state[x] for x in c.child])]
        def bottom(a):return a.reshape(-1,a.shape[-1]).max(axis=0)
        bend=np.maximum.reduce([(s[1]-1)*bottom(c.ybtab[i]) for s,i in c.ybu.items()])
        self.edge=bottom(g)+np.maximum(0,bottom(c.hvtab))+np.maximum(0,bottom(c.p4tab))+np.maximum(0,bend)
        self.hc=c.hctab.max(axis=1)
        self.yt=c.yttab.max(axis=(0,2));self.s3=c.s3tab.max(axis=1)
        self.by_k={}
        for h,k in product((0,1),range(c.KF+2)):
            start=max(0,2-k-h)
            arr=np.array(c.f[c.vi(0,k,h):c.vi(c.LF,k,h)+1]) if k<=c.KF else np.array([c.f[c.vi(0,k,h)]])
            tail=max(0,c.D-k-h)
            if h:
                arr[tail:]+=max(self.edge[c.child[t]] for t in c.tail(0,k,h))
                for l in range(start,tail):arr[l]+=self.edge[c.child[l,k,h]]
            if 2<=k<=c.K:
                arr[tail:]+=comb(k,2)*max(self.hc[c.hc[t]] if t in c.hc else 0
                                             for t in c.tail(0,k,h))
                for l in range(start,tail):arr[l]+=comb(k,2)*self.hc[c.hc[l,k,h]]
            for index,values,m in ((c.yt,self.yt,2),(c.s3,self.s3,3)):
                for (l,j,bh),i in index.items():
                    if (j,bh)==(k,h):arr[l]+=comb(k,m)*values[i]
            self.by_k[h,k]=int(arr[start:].max())

    @cache
    def exact(self,l,k,h):
        c=self.c;s=(l,k,h);out=Q(int(c.f[c.vi(*s)]))
        for t,w in c.embed(s):
            if h:out+=w*int(self.edge[c.child[t]])
            if t in c.hc:out+=comb(k,2)*w*int(self.hc[c.hc[t]])
        if s in c.yt:out+=comb(k,2)*int(self.yt[c.yt[s]])
        if s in c.s3:out+=comb(k,3)*int(self.s3[c.s3[s]])
        assert out<=self.by_k[h,min(k,c.KF+1)]
        return out

    def cap(self,t):
        assert t<=0
        return max(v+t*(k-h) for (h,k),v in self.by_k.items())


def fixtures(c,charge):
    from audit import build, physical_view
    cases=[((-1,),[2]),((-1,0,1,2),[2,0,0,1]),
           ((-1,0,0,0,1,1,2,4),[2,0,0,1,0,1,1,18]),
           ((-1,0,1,1),[3,c.LF+2,1,1]),
           ((-1,0)+(1,)*(c.KF+1),[1,0]+[1]*(c.KF+1))]
    modes=set()
    for parents,leaves in cases:
        ch,states,_=physical_view(build(parents,leaves))
        parts=c.potential(ch,states,parts=True)
        assert sum(k-h for _,k,h in states)==0
        assert sum(parts)<=sum(charge.exact(*s) for s in states)
        modes.update(i for i,x in enumerate(parts) if x)
    # Add a Y-tail below an actual parent.
    ch,states,_=physical_view(build((-1,0,1,1),[2,0,1,1]))
    parts=c.potential(ch,states,parts=True)
    assert sum(parts)<=sum(charge.exact(*s) for s in states)
    modes.update(i for i,x in enumerate(parts) if x)
    assert modes==set(range(8)),modes
    return len(cases)+1


def physical_deficit(reward,cap):
    # Root census progression, arbitrary predictable core-neighbour gain.
    v=defaultdict(Q)
    for d in range(cap,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1;b=Q(l,d)
            a=reward[d,k]+b*v[d,k+1]
            v[d,k]=down(min((a+g*v[d+1,k])/(HIGH+g+b)
                            for g in (l,l+Q(k+1,2))))
    return LOW*v[2,0]


def renewal_deficit(reward,cap,steps):
    # Critical left Perron mixture, with unknown nonnegative weights: keep
    # the worst of the nineteen genuine label rows. Normalize occupation
    # with nu(b)=mu, since each core promotion adds one core.
    residual.MARK_CAP=20
    vals=[]
    for i in range(2,21):
        mass={residual.newborn(i):Q(1)}
        value=Q(0)
        for _ in range(steps+1):
            nxt=defaultdict(Q)
            for tree,m in mass.items():
                den,moves,*_=residual.system(tree,HIGH)
                w=down(m/den)
                _,l,ch=tree;k=len(ch);d=l+k+1
                value+=w*reward[d,k]
                for y,r in moves.items():
                    if residual.statistics(y)[0]<=cap:nxt[y]+=down(w*r)
            mass=nxt
        vals.append(value)
    # Denominator ell*G1 <= 1/LOW since ell*K*1=1 and nu(b)=mu.
    return LOW*min(vals),vals


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='+',type=Path)
    p.add_argument('--degree',type=int,default=64)
    p.add_argument('--size',type=int,default=8)
    p.add_argument('--steps',type=int,default=10)
    a=p.parse_args();assert 2<=a.size<=a.degree
    for path in a.paths:
        c=Cert(path);charge=Charge(c);S=c.S
        print(path.name, 'lambda',c.lam,'all eight motifs in',fixtures(c,charge),'physical charge fixtures',flush=True)
        zero=charge.cap(Q(0)); print('C0 from grouped audit',show(zero/S),flush=True)
        assert zero<=C0[c.mode]
        best=None
        for j in range(0,51):
            t=-Q(j*S,1000)
            cap=charge.cap(t)
            rewards={(d,k):(cap-charge.exact(d-k-1,k,1)-t*(k-1))/S
                     for d in range(2,a.degree+1) for k in range(d)}
            assert min(rewards.values())>=0
            reserve=physical_deficit(rewards,a.degree)
            v=cap/S-reserve
            if best is None or v<best[0]:best=(v,t,cap,rewards,reserve)
        v,t,cap,rewards,r=best
        rr,labels=renewal_deficit(rewards,a.size,a.steps)
        result=cap/S-max(rr,r)
        rounded=Q(-(-result*10**12//1),10**12)
        print('recharge',t/S,'all-degree cap',show(cap/S),'physical census credit',show(r,upper=False),
              '19-label prefix credit',show(rr,upper=False),'best label',labels.index(max(labels))+2,
              'worst label',labels.index(min(labels))+2,flush=True)
        print('piPhi <=',show(rounded),'rational',rounded,'C0 reduction',show(Q(C0[c.mode],S)-rounded),
              'relative reduction percent',show(100*(1-rounded/Q(C0[c.mode],S)),8,False),flush=True)
        assert rounded<Q(C0[c.mode],S)
        # Hold the earlier certified numerator fixed. No new row-slack
        # occupation is inferred from the denominator computation.
        old_edge={1:Q('.565773985336'),2:Q('.565674338794')}[c.mode]
        gap=abs(old_edge-c.lam)*Q(C0[c.mode],S)/rounded
        signed_edge=c.lam+(-gap if c.mode==1 else gap)
        print('same-numerator endpoint',show(signed_edge,upper=c.mode==1),flush=True)


if __name__=='__main__':main()
