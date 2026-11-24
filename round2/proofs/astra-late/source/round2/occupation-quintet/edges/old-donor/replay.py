#!/usr/bin/env python3
"""Old-donor FY: exact occupations, discrete donor bands only in the reward."""
import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from pathlib import Path
from time import process_time


spec = spec_from_file_location('fy_extended', Path(__file__).resolve().parents[1]/'fy/extended/replay.py')
ext = module_from_spec(spec)
spec.loader.exec_module(ext)


class OldDonor(ext.Pass):
    def __init__(self,e,a,upper):
        super().__init__(e,a,upper)
        self.ts=ext.states(a.donor,a.tstock)
        self.ti,self.tm=ext.axes(self.ts)
        self.nt=len(self.ts)
        self.n=len(self.rs)*self.nt
        self.zero=[0]*(self.n+1)
        self.t0=self.ti[ext.J]
        ends=sorted({5,a.donor,*[d for d in a.bins if 5<d<a.donor]})
        self.donor_bands=tuple(range(x+1,y+1) for x,y in zip(ends,ends[1:]))
        self.band_id={d:i for i,ds in enumerate(self.donor_bands) for d in ds}
        self.donor_id=[self.band_id.get(d,-1) for d,_ in self.ts]
        self.reward_states=0
        self.positive=0
        self.max_width=Q(0)

    def rewards(self,h):
        d,c=h
        e=self.e
        if d<=18 or c<2 or h in e.mids or h not in e.snap.inner:
            return {}
        knots=e.knots
        ak=tuple(e.fy(h,t) for t in knots)
        reciprocals={t:(Q(1,t),e.fy(h,Q(1,t))) for t in self.band_id}
        by={}
        for i,r in enumerate(self.rs):
            z=e.g((d+1,c),r)-e.g(h,r)
            top=max(e.sig*(t*z-a) for t,a in zip(knots,ak,strict=True))
            for b,ds in enumerate(self.donor_bands):
                vals=[top-e.sig*(u*z-a) for u,a in (reciprocals[t] for t in ds)]
                lo,hi=min(vals),max(vals)
                assert lo>=0
                self.bands+=1
                self.reward_states+=len(ds)
                self.positive+=sum(v>0 for v in vals)
                self.band_range+=hi-lo
                self.max_width=max(self.max_width,hi-lo)
                by[i,b]=self.round(hi if self.upper else lo)
        return by

    def layer(self,h,pg,pb):
        if h[1]<2:
            return self.zero
        reward=self.rewards(h)
        out=[0]*(self.n+1)
        nt=self.nt
        for i in range(len(self.rs)-1,-1,-1):
            r=self.rs[i]
            rg,rb=self.rm[i]
            for j in range(nt-1,-1,-1):
                t=self.ts[j]
                tg,tb=self.tm[j]
                k=i*nt+j
                out[k]=self.step((h,r,t),((1,2),(0,),(0,)),
                    (pg[k],out[rg*nt+j] if rg>=0 else 0,out[i*nt+tg] if tg>=0 else 0),
                    (pb[k],out[rb*nt+j] if rb>=0 else 0,out[i*nt+tb] if tb>=0 else 0),
                    reward.get((i,self.donor_id[j]),0))
        self.path_count+=self.n
        return out


def fixture(e):
    """Independent knot-piece evaluation and exact raw payload controls."""
    s=e.snap
    h=(38,2)
    for r in ((2,1),(4,1)):
        gain=e.g((39,2),r)
        old=e.g(h,r)
        z=gain-old
        raw=[e.fy(h,u) for u in e.knots]
        vals=[]
        for d in (5,6,8,18,32,90):
            x=Q(1,d)
            a,b=next((i,i+1) for i in range(len(raw)-1) if e.knots[i]<=x<=e.knots[i+1])
            p=(x-e.knots[a])/(e.knots[b]-e.knots[a])
            ax=(1-p)*raw[a]+p*raw[b]
            assert ax==e.fy(h,x)
            score=max(e.sig*((t-x)*z-v+ax) for t,v in zip(e.knots,raw))
            assert score>=0
            # Only the degree gate determines whether this exact score is paid.
            vals.append((d,score if d>5 else 0))
        print(s.side,h,r,'g pre/post',old,gain,'FY knot ticks',raw,'eligible reward',vals)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=90)
    p.add_argument('--stock',type=int,default=5)
    p.add_argument('--receiver',type=int,choices=(5,),default=5)
    p.add_argument('--rstock',type=int,default=5)
    p.add_argument('--donor',type=int,default=24)
    p.add_argument('--tstock',type=int,default=5)
    p.add_argument('--bins',type=int,nargs='*',default=[6,8,11,15,18,24,32,48,64,90,120,150])
    p.add_argument('--exact',action='store_true')
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--ceiling-only',action='store_true')
    p.add_argument('--fixtures',action='store_true')
    p.add_argument('--hash',action='store_true')
    a=p.parse_args()
    assert __debug__ and 38<=a.host<=179 and 6<=a.donor and 2<=a.stock<=26
    assert 1<=a.rstock<=5 and 1<=a.tstock<=a.donor
    if a.exact:
        a.bins=list(range(6,a.donor))
    ext.selfcheck()
    receipt=json.loads((ext.PRIVATE/'gold.json').read_bytes())
    for side in a.side or ('lower','upper'):
        assert receipt[side]['ok'] and receipt[side]['verify']['sha256']==ext.PINS[side][2]
        with ext.Snapshot(ext.PRIVATE/f'quintet_{side}.bin',side,verify=a.hash) as s:
            e=ext.Edges(s)
            if a.fixtures:
                fixture(e)
                continue
            for upper in ((True,) if a.ceiling_only else (False,True) if a.ceiling else (False,)):
                t=process_time(); engine=OldDonor(e,a,upper)
                print(side,'ceiling' if upper else 'floor',vars(a),flush=True)
                cs=engine.run()
                for name,v in zip(('mother pays','mother donates','two children, both orders'),cs,strict=True):
                    print(name,v,float(v))
                print('TOTAL',sum(cs),float(sum(cs)),'states',engine.path_count,'band rows',engine.bands,
                      'positive/exact degree rewards',engine.positive,engine.reward_states,
                      'mean/max conditional bin width (ticks)',float(engine.band_range/max(1,engine.bands)),
                      float(engine.max_width),'CPU seconds',round(process_time()-t,3),flush=True)
            s.check_unchanged()
    print('PASS exact corner checks and pinned copy unchanged')


if __name__=='__main__':
    main()
