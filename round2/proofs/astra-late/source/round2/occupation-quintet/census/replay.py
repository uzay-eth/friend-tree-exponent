#!/usr/bin/env python3
"""Export exact original singleton census rows and test the source conventions."""
import argparse
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess

from census import PINS, SOURCE, Snapshot

BASE = 'round2/neon2/quintet-cert/'
BLOBS = {'quintet_format.hpp':'647a892cfbcbea11efd8c2bc5dc9898f8d564e06',
         'verify_quintet.cpp':'6945724f2981a8d4a698e1d8330c0023323cab84',
         'verify_tail_bernstein.py':'461f48541fe968083ddf0bad9031ca63bf88d256'}


def source_pins(path):
    if subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip()!=SOURCE:
        raise ValueError('wrong detached source commit')
    for name,pin in BLOBS.items():
        with (path/BASE/name).open('rb') as f:
            raw=f.read()
        assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==pin


def controls(s):
    cases=((2,0),(2,1),(2,2),(s.K2+1,s.K2+1),(s.d2-1,s.K2+1),
           (s.d2-1,2),(s.d2,2),(s.DF-1,1),(s.DF-1,s.DF-1),(s.DF,1),(100000,4))
    wrong=[]
    for d,c in cases:
        a,b=s.row_ticks(d,c); f=s.f(d,c); delta=s.f(d+1,c)-f; l=d-c
        new=0 if not l else Q(l,d)*(s.f(d,c+1)+s.f(2,1)+s.newborn_edge(d,c+1)-f)
        for rho in (0,Q(c,2)):
            center=s.center(d,c)
            physical=(l+rho)*delta+new-s.lam*f
            outgoing=(rho-c*center)*delta if center is not None else 0
            assert s.sign*(physical-outgoing)==a+rho*b
            if center is None and c and delta and rho:
                # Applying an unconditional frozen rate c/4 to a genuinely
                # uncentered row loses physical slope even at this endpoint.
                planted=s.sign*((l+Q(c,4))*delta+new-s.lam*f)
                assert planted!=a+rho*b
                wrong.append((d,c,planted-(a+rho*b)))
    assert wrong
    print(s.side,'exact physical cases',len(cases),'wrong-centering at',wrong[0][0],wrong[0][1],
          'planted minus physical ticks',str(wrong[0][2]),flush=True)


def export(s,degree,kmax,receipt=None):
    stocks={c:s.stock_floor(c) for c in range(kmax+1)}
    rows=[]
    for d in range(2,degree+1):
        for c in range(d+1):
            a,b=s.row(d,c)
            rows.append((d,c,a,b))
    doc={'source':SOURCE,'source_blobs':BLOBS,'sha256':PINS[s.side][2],
            'orientation':'inward +lower / -upper; not adjusted to the certified endpoint',
            'header':s.header,'newborn_f_ticks':s.newborn_f,'newborn_inf_edge_ticks':s.ginf,
            'newborn_tail_edge_ticks':s.tail_edges,'degree_cap':degree,
            'stock_columns':['finite_min','argmin_D','analytic_tail_floor','union_floor'],
            'stock_bands':stocks,'row_columns':['D','core_degree','intercept','total_receive_slope'],'rows':rows}
    if receipt is not None:
        from host import HostCaps
        cap=HostCaps(s,receipt); h=s.header
        # Independently enumerate ordered internal partners satisfying the
        # degree-sum masks; verify the closed-form multiplicity maxima.
        inner=lambda a,k,m: [(d,c) for d in range(m,a) for c in range(m,min(d,k)+1)]
        p4=inner(h['d4i'],h['K4i'],2); p5=inner(h['d5i'],h['K5i'],2); b5=inner(h['d5c'],h['K5c'],2)
        from math import comb
        for d in range(2,max(h['d4i'],h['d5i'],h['d5o'])):
            for c in range(2,d+1):
                _,b,_,v,w,_=cap.multiplicities(d,c)
                assert b==(comb(c,2)*max((j-1 for i,j in p4 if d+i<=h['SP4']),default=0) if (d,c) in p4 else 0)
                assert v==(comb(c,2)*max(((j-1)*(l-1) for (i,j),(k,l) in product(p5,repeat=2)
                                         if d+i+k<=h['SP5']),default=0) if (d,c) in p5 else 0)
                assert w==(c*comb(c-1,2)*max((j-1 for i,j in b5 if d+i<=h['SB5']),default=0)
                           if d<h['d5o'] and 3<=c<=h['K5o'] else 0)
        doc['host_columns']=['D','core_degree','positive_higher_cap_ticks','recharged_comparison_low_ticks','unrecharged_comparison_high_ticks']
        doc['host'] = [(d,c,cap.higher_ticks(d,c),*cap.comparison_ticks(d,c))
                       for d,c,*_ in rows]
        assert all(0<=v<=cap.high and (lo<=hi or cap.alpha>0 and c>2) for d,c,v,lo,hi in doc['host'])
        print(s.side,'checked host degree masks; host cap range',min(r[2] for r in doc['host']),max(r[2] for r in doc['host']),flush=True)
    return doc


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('payload',type=Path); p.add_argument('side',choices=PINS)
    p.add_argument('--source',type=Path,required=True); p.add_argument('--degree',type=int,default=96)
    p.add_argument('--kmax',type=int,default=12); p.add_argument('--out',type=Path)
    p.add_argument('--compare',type=Path)
    p.add_argument('--receipt',type=Path,help='optional original golden receipt, enables sparse host/weight export')
    a=p.parse_args()
    if not __debug__:p.error('assertions required')
    if not (2<=a.degree<=512 and 0<=a.kmax<=32):p.error('unsupported export cut')
    source_pins(a.source)
    with Snapshot(a.payload,a.side) as s:
        print('SHA256',s.side,PINS[s.side][2],flush=True)
        controls(s); data=export(s,a.degree,a.kmax,a.receipt); s.check_unchanged()
        if a.out:
            with a.out.open('x') as o:json.dump(data,o,default=str,sort_keys=True,separators=(',',':'));o.write('\n')
        if a.compare:
            assert json.loads(a.compare.read_bytes())==json.loads(json.dumps(data,default=str))
            print('exact export reproduced',a.compare,flush=True)
        print('rows',len(data['rows']),'stock bands',[(k, float(v[0]),v[1],float(v[2]),float(v[3])) for k,v in data['stock_bands'].items()],flush=True)


if __name__=='__main__':main()
