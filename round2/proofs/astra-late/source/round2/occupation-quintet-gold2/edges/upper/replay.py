#!/usr/bin/env python3
"""Gold2 upper witness: disjoint small-payer FY and degree-18 edge rebates."""
import argparse
from fractions import Fraction as Q
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
import json
import mmap
import os
from pathlib import Path
import struct
from math import lcm
from time import process_time

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('gold2_edges', HERE.parent/'fy.py')
base = module_from_spec(spec)
spec.loader.exec_module(base)
ext = base.ext
SHA = 'f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876'
RECEIPT_SHA = '5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'


class Upper(base.Snapshot):
    def __init__(self, verify):
        self.path = base.PRIVATE/'quintet_upper.bin'
        self.side = 'upper'
        with self.path.open('rb') as fp:
            self.stat = os.fstat(fp.fileno())
            assert self.stat.st_size == 6344902168
            raw = struct.unpack('<Q63q', fp.read(512))
            h = self.header = dict(zip(ext.Snapshot.__init__.__globals__['FIELDS'], raw))
            assert tuple(h[k] for k in ('magic','mode','P','Q','S','nknots','d2','K2','d3m','K3m')) == (
                0x4654325F55355631,1,56573686,100000000,1<<40,20,252,22,80,10)
            if verify:
                fp.seek(0)
                assert hashlib.file_digest(fp, 'sha256').hexdigest() == SHA
            self.data = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
        self.lam, self.scale, self.sign = Q(h['P'], h['Q']), h['S'], -1
        self.DF, self.d2, self.K2 = (h[k] for k in ('DF','d2','K2'))
        self.inner = ext.Edges.exacts(self.d2, self.K2, 1)
        self.tail_base = len(self.inner)+self.d2-1-self.K2
        self.inf = self.tail_base+self.K2+1
        self.off_f = len(self.inner)*(1+h['nknots'])
        self.off_g = self.off_f+self.fid(self.DF,0)
        self.newborn_f = self.f(2,1)
        self.ginf = self.g(self.inf,0)
        self.tail_edges = tuple(self.g(i,0) for i in range(self.tail_base,self.inf))
        self.check_unchanged()


def audit(e):
    assert e.sig == 1 and len(e.knots) == 20
    assert e.knots[0] == 0 and e.knots[-1] == Q(1,2)
    for h in ((80,2),(26,11),(251,3)):
        assert h in e.snap.inner and h not in e.mids
        for payer in ((2,1),(4,2),(5,1)):
            old = e.g(h,payer)
            z = e.g((h[0]+1,h[1]),payer)-old
            values = [e.fy(h,x) for x in e.knots]
            for d in (2,5,6,12,36,251):
                t = Q(1,d)
                i = next(i for i in range(19) if e.knots[i]<=t<=e.knots[i+1])
                x,y = e.knots[i:i+2]
                a = ((y-t)*values[i]+(t-x)*values[i+1])/(y-x)
                gap = max((x-t)*z-k+a for x,k in zip(e.knots,values,strict=True))
                assert gap >= 0 and a == e.fy(h,t)
                donors = [t]*(h[1]-1)
                row = e.row(h,payer,donors,[Q(1,5)]*(payer[1]-1))
                assert row+(h[1]-1)*gap <= e.cap
    print('PASS 54 upper physical charged-edge plus FY checks')


class Ceiling(base.Price):
    """Same physical mask, adverse high rates and upper reward envelopes."""
    def __init__(self,e,a):
        super().__init__(e,a)
        self.upper, self.disc = True, base.LO
        self.tail_birth = None

    def rewards(self,h):
        d,c=h; e=self.e
        if d<=18 or c<2 or h in e.mids or h not in e.snap.inner:return {}
        out={}
        for i,r in enumerate(self.rs):
            z=e.g((d+1,c),r)-e.g(h,r)
            def f(x):return e.sig*(x*z-e.fy(h,x))
            top=max(map(f,e.knots))
            for j,band in enumerate(self.donor_bands):
                out[i,j]=self.round(top-min(f(Q(1,d)) for d in band))
            t=Q(1,self.tail[0])
            out[i,len(self.donor_bands)]=self.round(top-min(map(f,(Q(0),t,*(x for x in e.knots if 0<x<t)))))
        return out

    def step(self,ss,adj,gs,bs,reward=0):
        if self.tail not in ss:
            return ext.Pass.step(self,ss,adj,gs,bs,reward)
        d=lcm(2,*(x[0] for x in ss)); terms=[]
        for i,s in enumerate(ss):
            n,c=s
            if s==self.tail:
                assert not bs[i]
                continue
            exact=[ss[j] for j in adj[i] if ss[j]!=self.tail]
            l=n-c
            w=(c-len(adj[i]))*(d//2)+(len(adj[i])-len(exact))*(d//self.tail[0])
            terms.append((l*d+sum(d//x for x,_ in exact),w,gs[i]))
            if l:terms.append((l*(d//n),0,bs[i]))
        return ext.solve(reward,terms,d,self.disc,True)

    def edge(self,h,ss,moves,eg,eb,pb,old_is_receiver):
        v=super().edge(h,ss,moves,eg,eb,pb,old_is_receiver)
        if h==ext.J and not old_is_receiver:self.tail_birth=v[-2]
        return v

    def singleton(self,ss,moves,birth):
        if ss[-1]!=self.tail:return super().singleton(ss,moves,birth)
        assert self.tail_birth is not None
        out=[0]*(len(ss)+1)
        # Tail promotion leaves the donor band unchanged, but can create a
        # host with the old mother edge. Its actual hazard is at most one.
        out[-2]=ext.solve(self.tail_birth,(),1,self.disc,True)
        for i in range(len(ss)-2,-1,-1):
            g,b=self.tm[i]
            out[i]=self.step((ss[i],),((),),(out[g],),(out[b]+birth[i],))
        return out[0]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--donor',type=int,nargs='+',default=[8,12])
    p.add_argument('--stock',type=int,default=17)
    p.add_argument('--fy-only',action='store_true')
    p.add_argument('--edge-only',action='store_true')
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--hash',action='store_true')
    a = p.parse_args()
    assert __debug__ and 2<=a.stock<=22 and all(7<=d<=32 for d in a.donor)
    receipt = (base.PRIVATE/'source/round2/neon2/quintet-cert/certificate_gold2.json').read_bytes()
    assert hashlib.sha256(receipt).hexdigest() == RECEIPT_SHA
    gold = json.loads(receipt)
    assert list(map(Q,gold['certified_interval_dec'])) == [base.LO,base.HI]
    row = gold['upper']
    assert row['ok'] and row['verify']['sha256'] == gold['upper_sha256'] == SHA and row['cap_g'] == 463
    assert (base.PRIVATE/'hashes.log').read_text().splitlines().count('quintet_upper.bin: OK') == 1
    ext.LO,ext.HI = base.LO,base.HI
    ext.Edges.price.__globals__.update(LO=base.LO,HI=base.HI)
    ext.selfcheck()
    with Upper(a.hash) as s:
        e = ext.Edges(s)
        e.cap = row['cap_g']
        print('Gold2 upper',s.lam,'cap_g',e.cap,'knots',len(e.knots),flush=True)
        audit(e)
        edge = None if a.fy_only else e.price(18)/s.scale
        if edge is not None:
            print('EDGE18',edge,float(edge),flush=True)
        a.host,a.receiver,a.rstock,a.exclude,a.mask,a.absorb = 251,5,5,18,'both',True
        a.bins = [6,8,11,15,18,24,32]
        for donor in (() if a.edge_only else a.donor):
            a.donor,a.tstock = donor,5
            t = process_time()
            b = (Ceiling if a.ceiling else base.Price)(e,a)
            v = b.run()
            print('FY',donor,'stock',a.stock,'mother pays / mother donates / two child orders',*map(str,v))
            total = sum(v)
            print('TOTAL FY', 'ceiling' if a.ceiling else 'floor',total,float(total),'states',b.path_count,'CPU',round(process_time()-t,3),flush=True)
            if edge is not None and not a.ceiling:
                print('TOTAL WITH EDGE18',total+edge,float(total+edge),flush=True)
        s.check_unchanged()
    print('PASS upper physical floors and unchanged authenticated copy')


if __name__ == '__main__':
    main()
