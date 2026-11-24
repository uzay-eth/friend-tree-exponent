#!/usr/bin/env python3
"""Gold2 upper: exact large-payer / small-donor FY floor and matched ceiling."""
import argparse
from fractions import Fraction as Q
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
from time import process_time


def load(name,path):
    spec=spec_from_file_location(name,path)
    mod=module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


HERE=Path(__file__).resolve().parent
payer=load('original_payer',HERE.parent/'price.py')
upper=load('gold2_upper_edges',HERE.parents[1]/'edges/upper/replay.py')


def physical_checks(e, b):
    n=0
    for h,p in (((80,2),(9,1)),((26,11),(12,3)),((251,3),(251,2)),
                ((26,11),(252,2)),((26,11),(1998,2))):
        assert h in e.snap.inner and h not in e.mids
        z=e.g((h[0]+1,h[1]),p)-e.g(h,p)
        words=[e.snap.item(len(e.snap.inner)+e.snap.inner[h]*20+i) for i in range(20)]
        assert words==[e.fy(h,k) for k in e.knots]
        for d in range(2,6):
            t=Q(1,d)
            j=next(i for i in range(19) if e.knots[i]<=t<=e.knots[i+1])
            l,r=e.knots[j:j+2]
            a=((r-t)*words[j]+(t-l)*words[j+1])/(r-l)
            gap=max((k-t)*z-w+a for k,w in zip(e.knots,words,strict=True))
            assert gap>=0 and a==e.fy(h,t)
            row=e.row(h,p,[t]*(h[1]-1),[Q(1,5)]*(p[1]-1))
            assert row+(h[1]-1)*gap<=e.cap
            if p in b.ri and d<=b.a.donor:
                assert b.rewards(h)[((p[0],p[0]),p[1],d)]==b.round(gap)
            n+=1
    print('PASS',n,'raw-knot interpolation, cap and physical incidence controls',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=135)
    p.add_argument('--stock',type=int,default=10)
    p.add_argument('--receiver',type=int,default=30)
    p.add_argument('--rstock',type=int,default=6)
    p.add_argument('--donor',type=int,default=5)
    p.add_argument('--mode',choices=('floor','ceiling','both','checks'),default='both')
    a=p.parse_args()
    assert __debug__ and 38<=a.host<=251 and 8<=a.receiver<=251
    assert 2<=a.stock<=22 and 2<=a.rstock<=22 and 2<=a.donor<=5
    a.mask='small'; a.tstock=a.donor; a.bins=list(range(6,a.receiver+1))
    receipt=(upper.base.PRIVATE/'source/round2/neon2/quintet-cert/certificate_gold2.json').read_bytes()
    assert hashlib.sha256(receipt).hexdigest()==upper.RECEIPT_SHA
    q=json.loads(receipt)
    assert list(map(Q,q['certified_interval_dec']))==[payer.gold.LO,payer.gold.HI]
    row=q['upper']; assert row['ok'] and row['verify']['sha256']==q['upper_sha256']==upper.SHA
    payer.ext.LO,payer.ext.HI=payer.gold.LO,payer.gold.HI
    payer.ext.selfcheck()
    with upper.Upper(True) as s:
        e=payer.ext.Edges(s); e.cap=row['cap_g']
        assert e.sig==1 and e.cap==463 and s.lam==Q(row['nominal'])
        assert len(e.knots)==20 and e.knots[0]==0 and e.knots[-1]==Q(1,2)
        print('SHA-authenticated GOLD2 upper',s.lam,'cap_g',e.cap,vars(a),flush=True)
        for ceiling in ((False,True) if a.mode=='both' else (a.mode=='ceiling',)):
            b=payer.Price(e,a,ceiling)
            b.fixtures(); physical_checks(e,b); b.bands=0
            if a.mode=='checks':
                break
            t=process_time(); v=b.run(); total=sum(v)
            print('COMPONENTS',*map(str,v),flush=True)
            print('TOTAL','ceiling' if ceiling else 'floor',total,float(total),
                  'ticks',total*s.scale, 'path states',b.path_count,'reward labels',b.bands,
                  'CPU',round(process_time()-t,3),flush=True)
        s.check_unchanged()
    print('PASS upper bytes unchanged and exact integer Bellman checks',flush=True)


if __name__=='__main__':
    main()
