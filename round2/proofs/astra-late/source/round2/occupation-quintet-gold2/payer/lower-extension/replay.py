#!/usr/bin/env python3
"""Original Gold2 lower FY wider payer: exact rates and 20-knot rewards."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from itertools import product
from math import lcm
from pathlib import Path
from time import process_time

spec=spec_from_file_location('payer_base',Path(__file__).resolve().parents[1]/'price.py')
base=module_from_spec(spec); spec.loader.exec_module(base)
ext=base.ext


class Wide(base.Price):
    def __init__(self,e,a):
        super().__init__(e,a,False)
        self.disc=ext.HI
        self.tail=(a.donor,0) if a.absorb else None
        if self.tail:
            self.ts=ext.states(a.donor-1,a.tstock)+[self.tail]
            self.ti,self.tm=ext.axes(self.ts)
            for i,(d,c) in enumerate(self.ts[:-1]):
                if d+1==a.donor:
                    self.tm[i]=(len(self.ts)-1,self.tm[i][1])
            self.nt=len(self.ts); self.t0=self.ti[ext.J]
            self.n=len(self.rs)*self.nt; self.zero=[0]*(self.n+1)
        self.unit=lcm(self.disc.denominator,*range(2,max(a.host,a.receiver,a.donor)+1))
        self.inv=[0]+[self.unit//d for d in range(1,max(a.host,a.receiver,a.donor)+1)]
        self.discount=self.disc.numerator*(self.unit//self.disc.denominator)

    def rewards(self,h):
        d,c=h; e=self.e
        if c<2 or h in e.mids or h not in e.snap.inner:
            return {}
        ts=range(2,self.a.donor+1) if self.a.mask=='small' else range(6,self.a.donor+1)
        knot=[(x,e.fy(h,x)) for x in e.knots]
        donor=[(t,Q(1,t),e.fy(h,Q(1,t))) for t in ts if not self.tail or t!=self.tail[0]]
        if self.tail:
            end=Q(1,self.tail[0]); xs=[Q(0),end,*[x for x in e.knots if 0<x<end]]
            donor += [(self.tail[0],x,e.fy(h,x)) for x in xs]
        scale=lcm(*(z.denominator for x,y in knot for z in (x,y)),
                  *(z.denominator for _,x,y in donor for z in (x,y)))
        def i(q): return q.numerator*(scale//q.denominator)
        knot=[(i(x),i(y)) for x,y in knot]
        donor=[(t,i(x),i(y)) for t,x,y in donor]
        out={}
        for r in self.rs:
            k,rc=r
            if k<=5 or max(d,k)<=18:
                continue
            z=e.g((d+1,c),r)-e.g(h,r)
            assert z.denominator==1
            z=z.numerator
            best=max(y-x*z for x,y in knot)
            for t,x,y in donor:
                gap=best+x*z-y
                assert gap>=0
                key=((k,k),rc,t)
                v=gap*ext.GRID//scale
                out[key]=min(v,out.get(key,v))
        self.bands+=len(out)
        return out

    def step(self, ss, adj, gs, bs, reward=0):
        num,den=reward*self.unit,self.discount
        flex=[]
        for i,(d,c) in enumerate(ss):
            if ss[i]==self.tail:
                continue
            lo=(d-c)*self.unit+sum(self.inv[ss[j][0]] for j in adj[i] if ss[j]!=self.tail); birth=(d-c)*self.inv[d]
            width=(c-len(adj[i]))*self.inv[2]+sum(self.inv[ss[j][0]] for j in adj[i] if ss[j]==self.tail)
            assert width>=0
            den+=lo+birth; num+=lo*gs[i]+birth*bs[i]
            if width:
                flex.append((gs[i],width))
        a,b=num,den
        for v,w in sorted(flex):
            if v*b>=a:
                break
            a+=v*w; b+=w
        z=a//b
        assert z>=0 and num-den*z+sum(w*min(0,v-z) for v,w in flex)>=0
        return z

    def check(self):
        n=0
        for d,c,t in product((3,9,19),(1,2,3),(2,5)):
            ss=((d,max(2,c)),(8,2),(t,1)); adj=((1,2),(0,),(0,))
            gs=(13,2,61); bs=(3,7,11); reward=12
            vals=[]
            widths=[Q(k-len(row),2) for (_,k),row in zip(ss,adj)]
            for cor in product(*((0,w) for w in widths)):
                g=[D-k+sum(Q(1,ss[j][0]) for j in adj[i])+cor[i] for i,(D,k) in enumerate(ss)]
                b=[Q(D-k,D) for D,k in ss]
                vals.append((reward+sum(x*v+y*w for x,y,v,w in zip(g,b,gs,bs)))/(self.disc+sum(g)+sum(b)))
            assert self.step(ss,adj,gs,bs,reward)==min(vals)//1
            n+=1
        if not self.tail:
            self.fixtures()
        else:
            h=(80,2); r=(19,2); e=self.e
            z=e.g((81,2),r)-e.g(h,r)
            f=lambda x:e.fy(h,x)-x*z
            xmax=Q(1,self.tail[0]); xs=[0,xmax,*[x for x in e.knots if 0<x<xmax]]
            reward=self.rewards(h)
            assert reward[((19,19),2,self.tail[0])]==self.round(max(map(f,e.knots))-max(map(f,xs)))
            assert all(k[0][0]>5 and k[2]>=6 for k in reward)
            assert all(self.tm[i][0]==len(self.ts)-1 for i,(d,c) in enumerate(self.ts) if d==self.tail[0]-1)
            print('PASS absorbing donor reciprocal seam, raw spline reward and mask',flush=True)
        print('PASS',n,'independent exact rate-corner checks',flush=True)

    def layer(self,h,pg,pb):
        if h[1]==1 and h[0]%25==0:
            print('host layers through',h[0],'path states',self.path_count,flush=True)
        return super().layer(h,pg,pb)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=251)
    p.add_argument('--stock',type=int,default=16)
    p.add_argument('--receiver',type=int,default=52)
    p.add_argument('--rstock',type=int,default=9)
    p.add_argument('--donor',type=int,default=5)
    p.add_argument('--tstock',type=int,default=5)
    p.add_argument('--mask',choices=('small','both-high'),default='small')
    p.add_argument('--absorb',action='store_true')
    a=p.parse_args()
    assert __debug__ and 38<=a.host<=251 and 19<=a.receiver<=251 and 2<=a.stock<=22 and 2<=a.rstock<=22
    assert 2<=a.tstock<=a.donor and (2<=a.donor<=5 if a.mask=='small' else 6<=a.donor<=251)
    assert not a.absorb or (a.mask=='both-high' and a.tstock>=a.donor-1)
    a.bins=list(range(6,a.receiver+1))
    src=base.gold.PRIVATE/'source/round2/neon2/quintet-cert/certificate_gold2.json'
    raw=src.read_bytes()
    assert sha256(raw).hexdigest()=='5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
    q=json.loads(raw)
    assert q['lower']['ok'] and q['lower']['verify']['sha256']==q['lower_sha256']==base.gold.DIGEST
    ext.LO,ext.HI=(Q(q[s]['certified_exact']) for s in ('lower','upper'))
    with base.gold.Snapshot(verify=True) as s:
        e=ext.Edges(s); assert e.sig==-1 and s.lam==Q(q['lower']['nominal'])
        assert a.host<s.d2 and len(e.knots)==20 and e.knots[0]==0 and e.knots[-1]==Q(1,2)
        b=Wide(e,a)
        print('GOLD2 lower SHA passed',vars(a),flush=True)
        b.check(); b.bands=0
        t=process_time(); vs=b.run(); v=sum(vs)
        s.check_unchanged()
        print('COMPONENTS',*map(str,vs),flush=True)
        print('TOTAL',v,float(v),'ticks',v*s.scale,float(v*s.scale),
              'states',b.path_count,'reward labels',b.bands,'CPU',process_time()-t,flush=True)
    print('PASS exact lower replay and unchanged authenticated payload',flush=True)


if __name__=='__main__': main()
