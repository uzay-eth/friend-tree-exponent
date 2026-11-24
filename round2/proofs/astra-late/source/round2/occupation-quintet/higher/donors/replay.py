#!/usr/bin/env python3
"""Exact H3/P4 donor incidence fixtures and all-age payment/box-gap prices."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from price import BRANCHES, GRID, price, solve_row
from screen import HI, LO, Rows, PATH, layout, gain, change


class Donors(Rows):
    @lru_cache(None)
    def parts(self,p):
        # Reuse only the narrow physical row evaluator. Recover the individual
        # receipt terms separately for raw fixtures and inequality checks.
        a,b,e,moved=Rows.row(self,p)
        boxes=self.boxes(p)
        ext=sum(max(z*t for t in box) for z,box in zip(b,boxes))
        return a,b,e,ext,moved

    def details(self,p):
        n=len(p)
        fam='23' if n==3 else '34p'
        co=self.s.G[n]
        _,j,_,_=self.s.envspec[fam]
        fine=self.s.G[j]
        combined=Q(0)
        for pp in (p,p[::-1]):
            w,x,*rem=pp
            z=self.path((gain(x),*rem))-self.path(tuple(pp[1:]))
            r=self.s.r1(x)
            donor=co.dist(w)
            kept=co.dist(rem[-1])
            # Sender retains r1*z per outside incidence. The true receive term
            # closes that exact RHS; the envelope excess does not alter sender.
            assert r*z+(Q(1,w[0])-r)*z==Q(1,w[0])*z
            if fam=='23':
                assert x in self.s.M['h3']
            print('  incidence',fam,'receiver',x,'donor',w,'retained',rem[-1],
                  'r1_int',r*self.s.h['Sr'],'fine corners/delta',
                  [(fine.labels[i],q,layout.source_delta(self.s,fam,pp[1:-1],(fine.labels[i],)))
                   for i,q in fine.dist(rem[-1]).items()])
            for (u,wu),(v,wv) in product(donor.items(),kept.items()):
                m,rad=self.s.envelope(fam,pp[1:-1],(v,))
                alpha=co.inv(u)-r
                print('    coarse',co.labels[u],co.labels[v],'weight',wu*wv,
                      'm/r',m,rad,'signed charge',self.sig*alpha*m+abs(alpha)*rad)
            ex=self.sig*(Q(1,w[0])-r)*z
            bound=self.project(fam,pp[1:-1],(w,),(rem[-1],),z,x)
            assert bound>=ex
            combined+=bound-ex
            print('    true',ex,'bound',bound,'gap',bound-ex)
        return combined

    def fixtures(self):
        fixtures=[(((48,2),(2,2),(2,1)),(Q(1,4),0,0)),
                  (((49,2),(2,2),(2,1)),(Q(1,4),0,0)),
                  (((2,2),(16,2),(48,1)),(Q(1,7),0,0)),
                  (((2,2),(12,9),(2,1)),(Q(1,5),Q(7,4),0)),
                  (((18,2),(2,2),(2,2),(2,1)),(Q(1,4),0,0,0)),
                  (((2,2),(8,8),(2,2),(18,1)),(Q(1,4),Q(3,2),0,0))]
        for p,rho in fixtures:
            a,b,e,ext,m=self.parts(p)
            out=a+sum(z*x for z,x in zip(b,rho))
            bar=a+e+ext
            assert e>=0 and bar>=out
            assert bar<=({'h3':(225,408),'p4':(42,41)}[PATH[len(p)]][self.sig==1])
            print(self.side,'row',p,'exterior',rho,'centered incidences',m,'raw slopes',b,
                  'physical',out,'envelope',e,'exterior gap',ext-sum(z*x for z,x in zip(b,rho)),
                  'box',bar,'negative direct',max(0,-out),'increment beyond old negative box',
                  max(0,-out)-max(0,-bar))
            assert self.details(p)==e
            direct, negative_box, gap=max(0,-out),max(0,-bar),bar-out
            inc=direct-negative_box
            union=negative_box+gap
            assert 0<=inc<=gap and union==max(direct,gap)
            assert out+union==max(0,bar)


class Reward(Donors):
    def __init__(self,side,mode,old,virgin):
        super().__init__(side)
        self.mode,self.old,self.virgin=mode,old,virgin

    @lru_cache(None)
    def row(self,p):
        # Reward masks are data-independent. In particular 'virgin' below
        # means low tip state, not true age; no timer or discounted clock test.
        if sorted((p[0][0],p[-1][0]))[1]<self.old or min(p[0][0],p[-1][0])>self.virgin:
            return Q(0),(Q(0),)*len(p),Q(0),0
        a,b,e,ext,moved=self.parts(p)
        if self.mode=='gap':
            # Unbooked diagnostic. Production virtual upper row minus the
            # physical one; not a credit above a block cap or host reserve.
            a=-e-ext
        elif self.mode=='box':
            a,b=a+e+ext,(0,)*len(p)
        elif self.mode=='increment':
            # Increment after a hypothetical booking of max(0,-box).
            a+=max(0,-(a+e+ext))
        elif self.mode=='union':
            # gap + negative-box credit = max(direct, gap). Both terms have
            # disjoint algebra and pay against the same nonnegative cap.
            a-=max(0,a+e+ext)
        else:
            assert self.mode=='physical'
        return a,b,Q(0),moved


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fixtures',action='store_true')
    ap.add_argument('--caps',type=int,nargs='+')
    ap.add_argument('--shape',choices=('chain',*BRANCHES),default='chain')
    ap.add_argument('--extra',type=int,default=0)
    ap.add_argument('--old',type=int,default=2)
    ap.add_argument('--virgin',type=int,default=10**6)
    ap.add_argument('--mode',choices=('physical','increment','box','gap','union'),default='physical')
    ap.add_argument('--ceiling',action='store_true')
    ap.add_argument('--mid-arity',type=int,help='restrict the H3 middle to this terminal arity')
    ap.add_argument('--mid-min',type=int,help='include all terminal middle arities from this value to --mid-arity')
    a=ap.parse_args()
    if not a.fixtures:
        assert a.caps and len(a.caps) in (3,4) and min(a.caps)>=2
        assert a.shape=='chain' or len(a.caps)==len(BRANCHES[a.shape][1])
    for side in ('lower','upper'):
        if a.fixtures:
            Donors(side).fixtures()
            continue
        r=Reward(side,a.mode,a.old,a.virgin)
        if a.mid_arity:
            assert a.shape=='chain' and len(a.caps)==3 and not a.extra
            v=fixed_middle(r,a.caps,a.mid_arity,a.ceiling,a.mid_min)
        else:
            v=price(r,a.caps,extra=a.extra,shape=a.shape,upper=a.ceiling)
        print(side,a.mode,a.shape,a.caps,a.extra,'old/virgin',a.old,a.virgin,'ceiling' if a.ceiling else 'floor',
              'ticks/core',float(v),'exact ticks/core',v,'units/core',float(v/r.s.h['S']),
              'of target',float(v*200000/r.s.h['S']),flush=True)


def fixed_middle(r,caps,M,upper,minimum=None):
    """H3 chain, only middle may have side children; count each future last tip.

    At a promotion of a prefix tip in stage 2, the vertex both continues in
    stage 2 and registers one possible H3. Existing terminal rows may gain
    additional siblings. Other roles use strict caps, so the stage-1 seed
    can only use its first child and each terminal last tip is a core leaf.
    """
    assert 2<=M<=caps[1] and M<=r.s.h['K3m'] and caps[1]<r.s.h['d3m']
    minimum=M if minimum is None else minimum
    assert 2<=minimum<=M
    discount=LO if upper else HI
    nxt={}
    for k in (3,2,1):
        table={}
        cnt=pos=0
        for D in product(*(range(cap,1,-1) for cap in caps[:k])):
            middle=range(min(M,D[1]),k-2,-1) if k>1 else (1,)
            for cm in middle:
                c=(1,) if k==1 else (2,cm) if k==2 else (2,cm,1)
                p=tuple(zip(D,c))
                boxes=tuple((0,Q(cj-int(i>0)-int(i<k-1),2)) for i,cj in enumerate(c))
                a,b,_,_=r.row(p) if k==3 and cm>=minimum else (0,(0,)*k,0,0)
                pos+=a+sum(max(x*y for y in z) for x,z in zip(b,boxes))<0
                den,num=discount,Q(0)
                gains=[table.get((change(D,i,d+1),cm),0) for i,d in enumerate(D)]
                for i,(d,cj) in enumerate(p):
                    birth=Q(d-cj,d)
                    haz=d-cj+sum(Q(1,D[j]) for j in (i-1,i+1) if 0<=j<k)
                    den+=haz+birth
                    num+=haz*gains[i]
                    if i==1:
                        num+=birth*table.get((D,cm+1),0)
                    if i==k-1 and k<3:
                        num+=birth*nxt.get((D+(2,),1 if k==1 else cm+1),0)
                val,checks=solve_row(a,b,boxes,den,num,gains,upper)
                table[D,cm]=val
                cnt+=checks
        print(r.side,'middle arity band',(minimum,M),'prefix',k,'states',len(table),'checks',cnt,
              'robust positive',pos,flush=True)
        nxt=table
    return (HI if upper else LO)*Q(nxt[(2,),1],GRID)


if __name__=='__main__':
    main()
