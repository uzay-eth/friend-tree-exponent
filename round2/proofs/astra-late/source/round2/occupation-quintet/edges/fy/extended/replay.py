#!/usr/bin/env python3
"""Rolling directed degree/age FY bounds, disjoint from the pinned masks."""
import argparse
from fractions import Fraction as Q
import json
from math import lcm
from pathlib import Path
import sys
from time import process_time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from replay import Edges, Snapshot, PRIVATE, PINS, LO, HI, J

GRID = 10**5  # fixed-point ticks in the continuation and reward


def solve(reward, terms, degrees, discount, upper):
    """Integer Bellman sweep and exact adverse-corner check.

    Each term is (known rate numerator, variable width numerator, value).
    All rates share the given denominator. The discount is exact.
    """
    d, q = degrees, discount.denominator
    a = reward*d*q
    b = discount.numerator*d
    flex = []
    for r,w,v in terms:
        a += q*r*v
        b += q*r
        if w:
            flex.append((v,w*q))
    aa,bb = a,b
    for v,w in sorted(flex,reverse=upper):
        if (v*b <= a if upper else v*b >= a):
            break
        a += w*v
        b += w
    v = -(-a//b) if upper else a//b
    gap = aa-bb*v+sum(w*(max if upper else min)(0,z-v) for z,w in flex)
    assert (gap<=0 if upper else gap>=0) and v>=0
    return v


def states(d,c):
    return [(n,k) for n in range(2,d+1) for k in range(1,min(n,c)+1)]


def axes(ss):
    idx = {s:i for i,s in enumerate(ss)}
    return idx,[(idx.get((d+1,c),-1),idx.get((d,c+1),-1)) for d,c in ss]


class Pass:
    def __init__(self,e,a,upper):
        self.e,self.a,self.upper=e,a,upper
        self.disc=LO if upper else HI
        self.scale=e.snap.scale
        self.rs,self.ts=states(a.receiver,a.rstock),states(a.donor,a.donor)
        self.ri,self.rm=axes(self.rs); self.ti,self.tm=axes(self.ts)
        # Preserve the full receiver degree/arity in the transition. Bin only
        # the reciprocal in the reward, with no early arrival at D=6.
        stops=sorted({5,a.receiver,*[v for v in a.bins if 5<v<a.receiver]})
        self.bins=dict((d,(x+1,y)) for x,y in zip(stops,stops[1:]) for d in range(x+1,y+1))
        self.nt=len(self.ts); self.n=len(self.rs)*self.nt
        self.zero=[0]*(self.n+1)
        self.r0=self.ri[J]; self.t0=self.ti[J]
        self.path_count=0
        self.bands=0
        self.band_range=0

    def round(self,r):
        n=r.numerator*GRID; d=r.denominator
        return -(-n//d) if self.upper else n//d

    def rewards(self,h):
        """Exact discrete reciprocal-bin extrema for a physical donor degree."""
        d,c=h
        e=self.e
        if c<2 or h in e.mids:
            return {}
        if h not in e.snap.inner:
            return {}
        knots=e.knots
        ka=[e.fy(h,t) for t in knots]
        ta={i:e.fy(h,Q(1,i)) for i in range(2,self.a.donor+1)}
        by={}
        for r in self.rs:
            k,rc=r
            if k<=5 or max(d,k)<=18:
                continue
            z=e.g((d+1,c),r)-e.g(h,r)
            m=max(e.sig*(t*z-ak) for t,ak in zip(knots,ka,strict=True))
            for t in ta:
                v=m-e.sig*(Q(1,t)*z-ta[t])
                assert v>=0
                key=(self.bins[k],rc,t)
                # Also key on whether the old edge-degree exclusion is met.
                # No corner value from inside the old mask enters this bin.
                if key in by:
                    lo,hi=by[key]
                    by[key]=min(lo,v),max(hi,v)
                else:
                    by[key]=v,v
        self.bands+=len(by)
        self.band_range+=sum(hi-lo for lo,hi in by.values())
        return {k:self.round(v[self.upper]) for k,v in by.items()}

    def step(self,ss,adj,gs,bs,reward=0):
        # Minimum common scale for all reciprocals and promotion denominators.
        d=lcm(2,*(s[0] for s in ss))
        terms=[]
        for i,(n,c) in enumerate(ss):
            l=n-c
            r=l*d+sum(d//ss[j][0] for j in adj[i])
            w=(c-len(adj[i]))*(d//2)
            assert w>=0
            terms.append((r,w,gs[i]))
            if l:
                terms.append((l*(d//n),0,bs[i]))
        return solve(reward,terms,d,self.disc,self.upper)

    def layer(self,h,pg,pb):
        d,c=h
        if c<2:
            return self.zero
        reward=self.rewards(h)
        out=[0]*(self.n+1)
        nt=self.nt
        for i in range(len(self.rs)-1,-1,-1):
            r=self.rs[i]; rd,rc=r; rg,rb=self.rm[i]
            bin_=self.bins.get(rd)
            for j in range(nt-1,-1,-1):
                t=self.ts[j]; tg,tb=self.tm[j]
                k=i*nt+j
                w=reward.get((bin_,rc,t[0]),0) if max(d,rd)>18 else 0
                out[k]=self.step((h,r,t),((1,2),(0,),(0,)),
                                 (pg[k],out[rg*nt+j] if rg>=0 else 0,out[i*nt+tg] if tg>=0 else 0),
                                 (pb[k],out[rb*nt+j] if rb>=0 else 0,out[i*nt+tb] if tb>=0 else 0),w)
        self.path_count+=self.n
        return out

    def edge(self,h,ss,moves,eg,eb,pb,old_is_receiver):
        out=[0]*(len(ss)+1)
        for i in range(len(ss)-1,-1,-1):
            g,b=moves[i]
            new=pb[i*self.nt+self.t0] if old_is_receiver else pb[self.r0*self.nt+i]
            out[i]=self.step((h,ss[i]),((1,),(0,)),(eg[i],out[g]),(eb[i]+new,out[b]))
        return out

    def singleton(self,ss,moves,birth):
        out=[0]*(len(ss)+1)
        for i in range(len(ss)-1,-1,-1):
            g,b=moves[i]
            out[i]=self.step((ss[i],),((),),(out[g],),(out[b]+birth[i],))
        return out[0]

    def run(self):
        a=self.a
        zero_r=[0]*(len(self.rs)+1); zero_t=[0]*(len(self.ts)+1)
        pg={}; rg={}; tg={}; sg={}
        for d in range(a.host,1,-1):
            pn={}; rn={}; tn={}; sn={}
            for c in range(min(d,a.stock),0,-1):
                h=(d,c)
                pb=pn.get(c+1,self.zero)
                er=rn.get(c+1,zero_r); et=tn.get(c+1,zero_t)
                pn[c]=self.layer(h,pg.get(c,self.zero),pb)
                rn[c]=self.edge(h,self.rs,self.rm,rg.get(c,zero_r),er,pb,True)
                tn[c]=self.edge(h,self.ts,self.tm,tg.get(c,zero_t),et,pb,False)
                sn[c]=self.step((h,),((),),(sg.get(c,0),),
                                (sn.get(c+1,0)+er[self.r0]+et[self.t0],))
            pg,rg,tg,sg=pn,rn,tn,sn
        # Mother–child paths: the first stage follows the future mother from
        # her own birth; her promotion creates the host and its mother edge.
        def maternal(ss,edge_at_j):
            idx,moves=axes(ss)
            return self.singleton(ss,moves,[edge_at_j[idx.get((d,c+1),-1)] for d,c in ss])
        src=[maternal(self.rs,rg[1]),maternal(self.ts,tg[1]),sg[1]]
        mu=HI if self.upper else LO
        return [mu*Q(x,GRID*self.scale) for x in src]


def selfcheck():
    for u in (False,True):
        for terms in ([(2,3,5),(2,3,1)],[(2,0,3),(0,7,9),(0,3,1)]):
            from itertools import product
            ratios=[Q(7+sum((r+x)*v for (r,_,v),x in zip(terms,ws)),
                      HI+sum(r+x for (r,_,_),x in zip(terms,ws)))
                    for ws in product(*[(0,w) for _,w,_ in terms])]
            best=(max if u else min)(ratios)
            assert solve(7,terms,1,HI,u)==(-(-best.numerator//best.denominator) if u else best.numerator//best.denominator)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=100)
    p.add_argument('--receiver',type=int,default=14)
    p.add_argument('--donor',type=int,default=3)
    p.add_argument('--stock',type=int,default=8)
    p.add_argument('--rstock',type=int,default=5)
    p.add_argument('--bins',type=int,nargs='+',default=[6,8,11,15,18,23,29,36,48,64])
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--ceiling-only',action='store_true')
    p.add_argument('--hash',action='store_true')
    a=p.parse_args()
    assert __debug__ and 38<=a.host<180 and 6<=a.receiver<=179 and 2<=a.donor<=18
    assert 2<=a.stock<=26 and 1<=a.rstock<=26
    selfcheck()
    receipt=json.loads((PRIVATE/'gold.json').read_bytes())['lower']
    assert receipt['ok'] and receipt['verify']['sha256']==PINS['lower'][2]
    with Snapshot(PRIVATE/'quintet_lower.bin','lower',verify=a.hash) as s:
        e=Edges(s)
        for upper in ((True,) if a.ceiling_only else (False,True) if a.ceiling else (False,)):
            start=process_time(); price=Pass(e,a,upper)
            print('START','upper bound' if upper else 'lower bound',vars(a),flush=True)
            vs=price.run()
            for name,v in zip(('mother pays','child pays with mother donor','child-child, both directed orders'),vs,strict=True):
                print(name,v,float(v))
            print('TOTAL',sum(vs),float(sum(vs)),'path states',price.path_count,'reciprocal bands',price.bands,
                  'mean conditional bin width (ticks)',float(price.band_range/max(1,price.bands)),
                  'CPU seconds',round(process_time()-start,3),flush=True)
        s.check_unchanged()
    print('PASS exact integer Bellman corner checks and pinned copy unchanged')


if __name__=='__main__':main()
