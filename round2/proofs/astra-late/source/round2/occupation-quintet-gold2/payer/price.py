#!/usr/bin/env python3
"""Original GOLD2 lower FY with paying-neighbour degree strictly above five."""
import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
from time import process_time

ROOT=Path(__file__).resolve().parents[1]
spec=spec_from_file_location('gold2_edges',ROOT/'edges/fy.py')
gold=module_from_spec(spec)
spec.loader.exec_module(gold)
ext=gold.ext


class Price(ext.Pass):
    def __init__(self,e,a,ceiling):
        super().__init__(e,a,ceiling)
        self.disc=gold.LO if ceiling else gold.HI
        self.ts=ext.states(a.donor,a.tstock)
        self.ti,self.tm=ext.axes(self.ts)
        self.nt=len(self.ts); self.n=len(self.rs)*self.nt
        self.zero=[0]*(self.n+1); self.t0=self.ti[ext.J]

    def rewards(self,h):
        d,c=h; e=self.e
        if c<2 or h in e.mids or h not in e.snap.inner:
            return {}
        knots=e.knots
        ak=[e.fy(h,x) for x in knots]
        donors=range(2,self.a.donor+1) if self.a.mask=='small' else range(6,self.a.donor+1)
        ta={t:e.fy(h,Q(1,t)) for t in donors}
        # The original discrete reciprocal, not an extremum over a payer bin.
        out={}
        for r in self.rs:
            k,rc=r
            if k<=5 or max(d,k)<=18:
                continue
            z=e.g((d+1,c),r)-e.g(h,r)
            best=max(e.sig*(x*z-akv) for x,akv in zip(knots,ak,strict=True))
            for t,av in ta.items():
                gap=best-e.sig*(Q(1,t)*z-av)
                assert gap>=0
                out[((k,k),rc,t)]=self.round(gap)
        self.bands+=len(out)
        return out

    def fixtures(self):
        count=0
        for h in ((15,2),(19,2),(38,2),(50,5),(179,3)):
            rows=self.rewards(h)
            if not rows:
                continue
            e=self.e
            r=(8,2) if h[0]>18 else (19,2)
            if r not in self.ri:
                continue
            best=lambda x:e.sig*(x*(e.g((h[0]+1,h[1]),r)-e.g(h,r))-e.fy(h,x))
            for d in range(2,self.a.donor+1):
                key=((r[0],r[0]),r[1],d)
                active=(d<=5 if self.a.mask=='small' else d>5)
                if active:
                    assert rows[key]==self.round(max(map(best,e.knots))-best(Q(1,d)))
                    count+=1
                else:
                    assert key not in rows
            assert all(k[0][0]>5 and max(h[0],k[0][0])>18 for k in rows)
        assert count
        print('fixtures',count,'direct exact-knot equalities and exclusions',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=100)
    p.add_argument('--stock',type=int,default=8)
    p.add_argument('--receiver',type=int,default=24)
    p.add_argument('--rstock',type=int,default=6)
    p.add_argument('--donor',type=int,default=5)
    p.add_argument('--tstock',type=int,default=5)
    p.add_argument('--mask',choices=('small','both-high'),default='small')
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--fixtures',action='store_true')
    a=p.parse_args()
    assert __debug__ and 38<=a.host<=179 and 8<=a.receiver<=179
    assert 2<=a.stock<=27 and 2<=a.rstock<=27 and 2<=a.tstock<=a.donor
    assert 2<=a.donor<=5 if a.mask=='small' else 6<=a.donor<=179
    a.bins=list(range(6,a.receiver+1))
    source=gold.PRIVATE/'source/round2/neon2/quintet-cert/certificate_gold2.json'
    receipt=json.loads(source.read_bytes())
    assert receipt['lower_sha256']==receipt['lower']['verify']['sha256']==gold.DIGEST
    assert list(map(Q,receipt['certified_interval_dec']))==[gold.LO,gold.HI]
    ext.LO,ext.HI=gold.LO,gold.HI
    ext.selfcheck()
    with gold.Snapshot(verify=True) as s:
        assert a.host<s.d2
        e=ext.Edges(s)
        assert len(e.knots)==20 and e.knots[0]==0 and e.knots[-1]==Q(1,2)
        print('SHA-authenticated GOLD2 lower',s.lam,vars(a),flush=True)
        b=Price(e,a,a.ceiling)
        b.fixtures()
        b.bands=0
        if not a.fixtures:
            t=process_time()
            v=b.run()
            for name,x in zip(('mother pays','mother donates','both child orders'),v,strict=True):
                print(name,x,float(x))
            print('TOTAL',sum(v),float(sum(v)),'states',b.path_count,'exact reward labels',b.bands,
                  'CPU',round(process_time()-t,3),flush=True)
        s.check_unchanged()
    print('PASS pinned lower bytes and exact Bellman inequalities',flush=True)


if __name__=='__main__':
    main()
