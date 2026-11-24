#!/usr/bin/env python3
"""Cold old-donor checks, followed by one serial matched owner replay."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

R = Path(__file__).resolve().parents[1]
E = R/'occupation-quintet/edges'
S = 2**40
PINS = {
    'old-donor/README.md': '3c9bc4b6dd59fa59adc7842bf7b3a170253041b1cd63929f2277109a9991859a',
    'old-donor/replay.py': '93fad1188fd39638a83e7c904b466015e96f2dd2072ee318ef6273c1266a7ce2',
    'fy/extended/replay.py': '4000aa83a0395eeb2324ef9fb97d2a64e0795c9fa08ab58f4774a0cf9d2d3aee',
    'fy/extended/cut36-floor.log': '56b31c406cfe1d794a85f6dea700211eaedbf297e5e1cd00470896e0bd28a824',
    'replay.py': 'b0c4f4a9af8238b7980d61d83882a81c37e5519b6c6a5534813037e69a853223',
    'replay.log': '4d2cc123b2e9826db5cc8315e86d283b36b1e1231ed6f4d3adff0bfacf555c0d',
    'fy/stock8.log': 'deb154ef476a27a643b7c4ff5f1d12cca07f74a474cd989d193c84032632a236',
    '../census/census.py': 'b6c133acc92ef36d651a92f4684852fc5ed241e8c0e88d3ac0fdad5554c2550b',
}


def import_file(name,path):
    spec = spec_from_file_location(name,path)
    mod = module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def support():
    counts = Counter()
    for n in range(2,8):
        for parents in product(*(range(i) for i in range(1,n))):
            adj = [set() for _ in range(n)]; paths=Counter()
            for i,h in enumerate(parents,1):
                for p in adj[h]:
                    paths[h,p,i]+=1; paths[h,i,p]+=1
                adj[h].add(i); adj[i].add(h)
                assert paths==Counter((v,p,d) for v,nb in enumerate(adj) for p,d in permutations(nb,2))
                for j in range(7):
                    dd = [max(len(nb),(2,5,6,18,19,24,37)[(q+j)%7]) for q,nb in enumerate(adj)]
                    for v,p,d in paths:
                        h,r,t=(dd[x] for x in (v,p,d))
                        gates = (max(h,r)<=18, max(h,r)>18 and r<=5 and t<=5,
                                 max(h,r)>18 and r>5 and t<=5, h>18 and r<=5 and t>5)
                        assert sum(gates)<=1
                        for k,on in enumerate(gates): counts[k]+=on
    assert len(counts)==4 and all(counts.values())
    print('live-mask ordered cases D18 / small-small / old-payer / old-donor:',*[counts[x] for x in range(4)],flush=True)


def controls(mod):
    cfg=SimpleNamespace(host=100,stock=7,receiver=5,rstock=5,donor=36,tstock=6,bins=[6,8,11,15,18,24,32,48,64,90,120,150])
    for side,sign,cap in (('lower',-1,401),('upper',1,946)):
        with mod.ext.Snapshot(mod.ext.PRIVATE/f'quintet_{side}.bin',side,verify=False) as s:
            e=mod.ext.Edges(s); assert e.sig==sign
            objs=[mod.OldDonor(e,cfg,u) for u in (False,True)]
            num=0; physical=0
            for h in ((15,15),(16,16),(18,17),(19,17),(37,3),(38,2),(100,7)):
                got=[p.rewards(h) for p in objs]; expected={}
                if h[0]>18 and h[1]>=2 and h in s.inner and h not in e.mids:
                    k=e.knots; start=len(s.inner)+len(k)*s.inner[h]; raw=[s.item(start+i) for i in range(len(k))]
                    def A(x):
                        i=next(i for i in range(len(k)-1) if k[i]<=x<=k[i+1]); l,r=k[i:i+2]
                        return ((r-x)*raw[i]+(x-l)*raw[i+1])/(r-l)
                    for i,r in enumerate(objs[0].rs):
                        d=s.g(s.inner[h[0]+1,h[1]],s.inner[r])-s.g(s.inner[h],s.inner[r])
                        def gap(t):
                            x=Q(1,t)
                            return max(sign*((u-x)*d-a+A(x)) for u,a in zip(k,raw))
                        for j,ds in enumerate(objs[0].donor_bands):
                            vs=[gap(t)*mod.ext.GRID for t in ds]
                            assert min(vs)>=0
                            expected[i,j]=(min(vs)//1,-(-max(vs)//1))
                        for t in (5,6,18,36,90):
                            # Physical tree fixture: other incidences at each endpoint
                            # can be completed by distinct rooted branches.
                            donors=[Q(1,t)]*(h[1]-1)
                            row=e.row(h,r,donors,[Q(1,4)]*(r[1]-1))
                            assert row+(h[1]-1)*(gap(t) if t>5 else 0)<=cap
                            physical+=1
                        if h==(38,2) and r==(2,1):
                            assert tuple(gap(t) for t in (6,90))==((Q(2816730),Q(286280377,15)) if sign<0 else (Q(156300613,20),Q(776383,60)))
                assert got==[{k:v[j] for k,v in expected.items()} for j in (0,1)]
                num+=len(expected)
            assert all(objs[0].donor_id[j]==-1 for j,(d,c) in enumerate(objs[0].ts) if d<=5)
            s.check_unchanged()
            print(side,'independent signed bands / raw row payout cases',num,physical,flush=True)


def run(mod,wide=False):
    params = ['--host','170','--stock','12','--donor','70','--tstock','7','--side','lower'] if wide else ['--host','100','--stock','7','--donor','36','--tstock','6','--ceiling']
    cmd=[sys.executable,'-Bu',str(E/'old-donor/replay.py'),*params]
    text=subprocess.check_output(cmd,text=True)
    assert text.splitlines()[-1]=='PASS exact corner checks and pinned copy unchanged'
    result={}
    for side,mode,body in re.findall(r'^(lower|upper) (floor|ceiling) \{.*?\}\n(.*?)(?=^(?:lower|upper) (?:floor|ceiling)|^PASS)',text,re.M|re.S):
        rows=body.splitlines(); assert len(rows)==4
        parts=[Q(row.split()[-2]) for row in rows[:3]]
        total=Q(rows[3].split()[1]); assert total==sum(parts)
        n=int(re.search(r' states (\d+) ',rows[3])[1]); assert n==(11819808 if wide else 1621200)
        disc=mod.ext.LO if mode=='floor' else mod.ext.HI
        assert all((v*S*mod.ext.GRID/disc).denominator==1 for v in parts)
        if wide: assert total==Q('80878351158765429162857/109951162777600000000000000000')
        result[side,mode]=total
        print(side,mode,'wide' if wide else 'matched',*parts,'TOTAL',total,'states',n,flush=True)
    assert len(result)==(1 if wide else 4)
    return result


def main():
    args=argparse.ArgumentParser();args.add_argument('--wide',action='store_true'); opt=args.parse_args()
    for p,h in PINS.items(): assert sha256((E/p).read_bytes()).hexdigest()==h,p
    support()
    mod=import_file('old_donor',E/'old-donor/replay.py')
    controls(mod)
    result=run(mod)
    # Distinguish the two original signs from the Bellman floor / ceiling directions.
    for side in ('lower','upper'):
        assert 0<result[side,'floor']<result[side,'ceiling']
    best=run(mod,wide=True) if opt.wide else result
    g=[Q(x)/S for x in re.findall(r'mask degree 18 .*?credit ticks/core (\S+)',(E/'replay.log').read_text())]
    fy=[Q(x)/S for x in re.findall(r'^TOTAL (\S+)',(E/'fy/stock8.log').read_text(),re.M)]
    payer=Q(re.search(r'^TOTAL (\S+)',(E/'fy/extended/cut36-floor.log').read_text(),re.M)[1])
    for i,side in enumerate(('lower','upper')):
        v=(best if side=='lower' else result)[side,'floor']
        print(side,'combined normalized edge floor',g[i]+fy[i]+v+(payer if side=='lower' else 0),flush=True)
    print('PASS old-donor original payout, disjointness and normalization',flush=True)


if __name__=='__main__': main()
