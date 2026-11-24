#!/usr/bin/env python3
"""FY both-high occupation with invariant absorbing neighbour degree bands."""
import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path
import sys
from time import process_time

ROOT=Path(__file__).resolve().parents[1]
spec=spec_from_file_location('fy_extended',ROOT/'extended/replay.py')
base=module_from_spec(spec);spec.loader.exec_module(base)
Pass,Edges,Snapshot,PRIVATE,LO,HI,J,GRID,states,axes,solve=(
    getattr(base,k) for k in 'Pass Edges Snapshot PRIVATE LO HI J GRID states axes solve'.split())


class Invariant(Pass):
    def __init__(self,e,a,upper):
        self.e,self.a,self.upper=e,a,upper
        self.disc=LO if upper else HI
        self.scale=e.snap.scale
        self.tail=(a.cut,0)
        ss=states(a.cut-1,a.cut-1)+[self.tail]
        self.rs=self.ts=ss
        self.ri,self.rm=axes(ss); self.ti,self.tm=self.ri,self.rm
        for i,(d,c) in enumerate(ss[:-1]):
            if d+1==a.cut:
                self.rm[i]=(len(ss)-1,self.rm[i][1])
        self.nt=len(ss);self.n=self.nt**2;self.zero=[0]*(self.n+1)
        self.r0=self.t0=self.ri[J]
        self.path_count=0; self.cells=0
        s=e.snap
        self.tail_ids=tuple(i for (d,c),i in s.inner.items() if d>=a.cut)+tuple(range(len(s.inner),s.inf+1))
        self.tail19=tuple(i for (d,c),i in s.inner.items() if d>=19)+tuple(range(len(s.inner),s.inf+1))

    def reward(self,h):
        d,c=h;e=self.e;s=e.snap
        if h in e.mids or c<2:return {}
        u=s.inner[h];v=s.inner.get((d+1,c),s.tail_base+c-1)
        xx=e.knots; aa=[e.fy(h,x) for x in xx]; sign=e.sig
        tail=self.tail_ids if d>18 else self.tail19
        zlo,zhi=0,0
        if d>18 or self.upper:
            zz=[s.g(v,i)-s.g(u,i) for i in tail]
            zlo,zhi=min(zz),max(zz)
        # The signed donor gap is convex in the paying slope. For a donor
        # reciprocal band, bound at every knot within the band and its ends.
        # For its minimum allow every slope between discrete extrema. The
        # minimum of the upper envelope minus a line lies at a hull breakpoint.
        hull=[]
        for x,a in sorted(zip(xx,aa),reverse=sign<0):
            m,b=sign*x,-sign*a
            while hull:
                pm,pb,pz=hull[-1]; z=(pb-b)/(m-pm)
                if pz is None or z>pz:break
                hull.pop()
            hull.append((m,b,(hull[-1][1]-b)/(m-hull[-1][0]) if hull else None))
        def f(z,ts):
            best=max(m*z+b for m,b,_ in hull)
            actual=[sign*(t*z-e.fy(h,t)) for t in ts]
            return best-(min if self.upper else max)(actual)
        out={}
        for i,r in enumerate(self.rs):
            if r[0]<6:continue
            if r==self.tail:
                if d<=18 and not self.upper:continue
                zs=[zlo,zhi]+([] if self.upper else [z for _,_,z in hull if z is not None and zlo<z<zhi])
            else:
                if max(d,r[0])<=18:continue
                zs=[e.g((d+1,c),r)-e.g(h,r)]
            for j in range(self.nt):
                t=self.ts[j]
                if t[0]<6:continue
                ts=(Q(0),Q(1,t[0]),*(x for x in xx if 0<x<Q(1,t[0]))) if t==self.tail else (Q(1,t[0]),)
                val=(max if self.upper else min)(f(z,ts) for z in zs)
                assert val>=0
                out[i*self.nt+j]=self.round(val)
        self.cells+=len(out)
        return out

    def step(self,ss,adj,gs,bs,reward=0):
        d=lcm(2,*(s[0] for s in ss)); terms=[]
        for i,s in enumerate(ss):
            n,c=s
            if s==self.tail:
                # Tail stock and degree no longer affect the state; only a
                # promotion that creates a newly tracked path can contribute.
                assert not bs[i]
                continue
            l=n-c
            exact=[ss[j] for j in adj[i] if ss[j]!=self.tail]
            r=l*d+sum(d//x for x,_ in exact)
            w=(c-len(adj[i]))*(d//2)+(len(adj[i])-len(exact))*(d//self.a.cut)
            assert w>=0
            terms.append((r,w,gs[i]))
            if l:terms.append((l*(d//n),0,bs[i]))
        return solve(reward,terms,d,self.disc,self.upper)

    def layer(self,h,pg,pb):
        if h[1]<2:return self.zero
        r=self.reward(h); out=[0]*(self.n+1); n=self.nt
        for i in range(n-1,-1,-1):
            g,b=self.rm[i]
            for j in range(n-1,-1,-1):
                x,y=self.tm[j];k=i*n+j
                out[k]=self.step((h,self.rs[i],self.ts[j]),((1,2),(0,),(0,)),
                                 (pg[k],out[g*n+j] if g>=0 else 0,out[i*n+x] if x>=0 else 0),
                                 (pb[k],out[b*n+j] if b>=0 else 0,out[i*n+y] if y>=0 else 0),r.get(k,0))
        self.path_count+=self.n
        return out

    def run(self):
        a=self.a
        zr=[0]*(self.nt+1); pg={};rg={};tg={};sg={}
        for d in range(a.host,1,-1):
            pn={};rn={};tn={};sn={}
            for c in range(min(d,a.stock),0,-1):
                h=(d,c);pb=pn.get(c+1,self.zero);er=rn.get(c+1,zr);et=tn.get(c+1,zr)
                pn[c]=self.layer(h,pg.get(c,self.zero),pb)
                rn[c]=self.edge(h,self.rs,self.rm,rg.get(c,zr),er,pb,True)
                tn[c]=self.edge(h,self.ts,self.tm,tg.get(c,zr),et,pb,False)
                sn[c]=self.step((h,),((),),(sg.get(c,0),),(sn.get(c+1,0)+er[self.r0]+et[self.t0],))
            pg,rg,tg,sg=pn,rn,tn,sn
        def mother(v):
            out=[0]*(self.nt+1)
            # At tail, use b∈[0,1] on the positive edge-birth source; its
            # self-loop cancels, leaving a constant super/subsolution.
            out[-2]=-(-(v[-2]*self.disc.denominator)//self.disc.numerator) if self.upper else 0
            for i in range(self.nt-2,-1,-1):
                g,b=self.rm[i]
                out[i]=self.step((self.rs[i],),((),),(out[g],),(out[b]+v[b],))
            return out[0]
        return [(HI if self.upper else LO)*Q(v,GRID*self.scale) for v in (mother(rg[1]),mother(tg[1]),sg[1])]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cut',type=int,default=8)
    p.add_argument('--host',type=int,default=179)
    p.add_argument('--stock',type=int,default=8)
    p.add_argument('--side',choices=('lower','upper','both'),default='both')
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--hash',action='store_true')
    a=p.parse_args(); assert __debug__ and 7<=a.cut<=19 and 38<=a.host<=179 and 2<=a.stock<=26
    for side in ('lower','upper') if a.side=='both' else (a.side,):
        with Snapshot(PRIVATE/f'quintet_{side}.bin',side,verify=a.hash) as s:
            start=process_time(); q=Invariant(Edges(s),a,a.ceiling);v=q.run()
            print(side,'ceiling' if a.ceiling else 'floor',vars(a),flush=True)
            for k,x in zip(('mother pays','child pays with mother donor','both child orders'),v,strict=True):
                print(k,x,float(x))
            print('TOTAL',sum(v),float(sum(v)),'states',q.path_count,'cells',q.cells,'CPU',round(process_time()-start,3),flush=True)
            s.check_unchanged()
    print('PASS')


if __name__=='__main__':main()
