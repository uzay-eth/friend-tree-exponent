#!/usr/bin/env python3
"""Disjoint upper Gold2 H3 mother-chain floor and matched ceiling."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path
from time import process_time
import json
import os

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def load(name,path):
    spec=spec_from_file_location(name,path)
    mod=module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


base=load('frontier_base',HERE.parent/'replay.py')


class Frontier(base.H3):
    def layer(self,h,pg,pb):
        if h[1]<2: return self.zero
        adj=((1,2),(0,),(0,)); nt=self.nt
        out=[0]*(self.n+1)
        for i in range(len(self.rs)-1,-1,-1):
            r=self.rs[i]; rg,rb=self.rm[i]
            for j in range(nt-1,-1,-1):
                t=self.ts[j]; tg,tb=self.tm[j]; k=i*nt+j
                gs=(pg[k],out[rg*nt+j] if rg>=0 else 0,out[i*nt+tg] if tg>=0 else 0)
                bs=(pb[k],out[rb*nt+j] if rb>=0 else 0,out[i*nt+tb] if tb>=0 else 0)
                p=h,r,t
                if not self.mask(*p):
                    out[k]=self.step(p,adj,gs,bs); continue
                outward,slopes=self.e.triple(*p)
                a=-outward
                widths=tuple(-s*Q(c-len(ns),2) for s,(_,c),ns in zip(slopes,p,adj))
                lo=a+sum(min(0,w) for w in widths); hi=a+sum(max(0,w) for w in widths)
                self.paid+=1
                if hi<=0:
                    out[k]=self.step(p,adj,gs,bs); continue
                d=lcm(2,*(x[0] for x in p)); terms=[]
                for (n,c),ns,y,z,w in zip(p,adj,gs,bs,widths):
                    L=n-c; width=(c-len(ns))*(d//2)
                    shift=self.round(w/width) if width else 0
                    terms.append((L*d+sum(d//p[x][0] for x in ns),width,y+shift,y))
                    if L: terms.append((L*(d//n),0,z,z))
                v=base.affine(a,terms,d,self.disc,self.upper)
                if lo<0: v=max(v,self.step(p,adj,gs,bs))
                assert v>=0
                out[k]=v
        self.path_count+=self.n
        return out

    def mask(self,h,r,t):
        a=self.a
        if h not in self.e.mids or not a.old<=r[0]<=a.last or r[1]<a.mother_min_stock:
            return False
        # These inequalities apply at payment, not at a tracked vertex's birth.
        if r[0]<=79 and r[1]<=8 and h[0]<=5 and h[1]<=4 and t[0]<=18 and t[1]<=3:
            return False
        return h[0]>=a.middle_min and t[0]>=a.child_min


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--host',type=int,default=5,help='middle degree cap')
    ap.add_argument('--stock',type=int,default=4)
    ap.add_argument('--receiver',type=int,default=125,help='mother degree cap')
    ap.add_argument('--rstock',type=int,default=10)
    ap.add_argument('--donor',type=int,default=8,help='child degree cap')
    ap.add_argument('--tstock',type=int,default=3)
    ap.add_argument('--old',type=int,default=80)
    ap.add_argument('--last',type=int)
    ap.add_argument('--mother-min-stock',type=int,default=2)
    ap.add_argument('--middle-min',type=int,default=2)
    ap.add_argument('--child-min',type=int,default=2)
    ap.add_argument('--ceiling',action='store_true')
    ap.add_argument('--fixtures',action='store_true')
    a=ap.parse_args(); a.reserve=False; a.last=a.last or a.receiver
    assert 19<=a.old<=a.last<=a.receiver and 2<=a.mother_min_stock<=a.rstock
    path=base.gold.PRIVATE
    assert (path/'hashes.log').read_text().splitlines()==['quintet_lower.bin: OK','quintet_upper.bin: OK']
    data=(path/'source/round2/neon2/quintet-cert/certificate_gold2.json').read_bytes()
    assert sha256(data).hexdigest()=='5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
    docs=json.loads(data); lo,hi=(Q(docs[s]['certified_exact']) for s in ('lower','upper'))
    base.ext.LO,base.ext.HI=lo,hi
    wire=load('gold2_frontier_wire',ROOT/'attacks/neon2-quintet/gold2/physical/check.py')
    t=process_time()
    p=wire.Gold(path/'quintet_upper.bin')
    try:
        with base.Snapshot('upper') as s:
            ident=lambda st:(st.st_dev,st.st_ino,st.st_size,st.st_ctime_ns,st.st_mtime_ns)
            assert ident(p.stat)==ident(s.stat) and s.lam==Q(docs['upper']['nominal'])
            e=base.Rows(s,p,docs['upper'])
            v=Frontier(e,a,a.ceiling,lo,hi)
            if a.fixtures: base.fixtures(e,wire.old)
            print('original upper Gold2',vars(a),flush=True)
            total=sum(v.run()); ticks=total*s.scale
            print('normalized',total,float(total),'ticks',ticks,float(ticks),
                  'states',v.path_count,'paid',v.paid,'CPU seconds',process_time()-t,flush=True)
            s.check_unchanged(); assert ident(os.fstat(p.fd))==ident(p.stat)
    finally:
        os.close(p.fd)


if __name__=='__main__': main()
