#!/usr/bin/env python3
"""Gold2 post-FY g slack and original signed H3 excess on old/young masks."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import lcm
from pathlib import Path
from time import process_time
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('gold2_fy',HERE.parent/'edges/fy.py')
gold = module_from_spec(spec)
spec.loader.exec_module(gold)
ext = gold.ext
GRID, J = ext.GRID, ext.J


class Snapshot(gold.Snapshot):
    def __init__(self,side):
        if side=='lower':
            super().__init__(False)
            return
        import mmap
        import os
        import struct
        self.path=gold.PRIVATE/'quintet_upper.bin'; self.side=side
        with self.path.open('rb') as f:
            self.stat=os.fstat(f.fileno()); assert self.stat.st_size==6344902168
            fields=self.fields=struct.unpack('<Q63q',f.read(512))
            h=self.header=dict(zip(ext.Snapshot.__init__.__globals__['FIELDS'],fields))
            assert (h['magic'],h['mode'],h['S'],h['nknots'])==(0x4654325F55355631,1,1<<40,20)
            self.data=mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ)
        self.lam=Q(h['P'],h['Q']); self.scale=h['S']; self.sign=-1
        self.DF,self.d2,self.K2=(h[x] for x in ('DF','d2','K2'))
        self.inner=ext.Edges.exacts(self.d2,self.K2,1)
        self.tail_base=len(self.inner)+self.d2-1-self.K2
        self.inf=self.tail_base+self.K2+1
        self.off_f=len(self.inner)*(1+h['nknots'])
        self.off_g=self.off_f+self.fid(self.DF,0)
        self.newborn_f=self.f(2,1); self.ginf=self.g(self.inf,0)
        self.tail_edges=tuple(self.g(i,0) for i in range(self.tail_base,self.inf))
        self.check_unchanged()


class Rows(ext.Edges):
    def __init__(self,s,p,r):
        super().__init__(s)
        self.p=p; self.cap=r['cap_g']
        self.hcap=r['verify']['caps']['h3']

    @lru_cache(maxsize=24000)
    def end_relaxed(self,x,y):
        a,b,fy=self.end(x,y)
        n=x[1]-1
        a=self.sig*a
        if fy:
            a+=n*max(self.sig*(t*b-self.fy(x,t)) for t in self.knots)
        return a,0 if fy else self.sig*b

    def edge(self,x,y):
        a,b=self.end_relaxed(x,y); z,w=self.end_relaxed(y,x)
        return self.cap+self.sig*self.snap.lam*self.g(x,y)-a-z,(-b,-w)

    @lru_cache(maxsize=32000)
    def triple(self,h,r,t):
        """Original outward H3 affine row, center first. Endpoints retained physically."""
        w=self.h3(r,h,t)
        a=-self.snap.lam*w
        slopes=[]
        p=(h,r,t)
        adj=((1,2),(0,),(0,))
        for i,(d,c) in enumerate(p):
            j=self.h3(*reorder(replace(p,i,(d+1,c))))-w
            n=c-len(adj[i]); L=d-c
            a+=(L+sum(Q(1,p[k][0]) for k in adj[i]))*j
            # Coverage for an actual external donor is independent of its label.
            ext_supported=self.extension(i,p,J) is not None if n else False
            if ext_supported:
                a+=n*self.snap.center(d,c)*j
            slopes.append(self.sig*j if n and not ext_supported else 0)
            if L:
                pp=replace(p,i,(d,c+1))
                a+=Q(L,d)*(self.h3(*reorder(pp))-w+(self.extension(i,pp,J) or 0))
        center=self.snap.center(*h)
        for donor,keep in ((r,t),(t,r)):
            j=self.g((h[0]+1,h[1]),keep)-self.g(h,keep)
            a+=(Q(1,donor[0])-center)*j
        return self.sig*a,tuple(slopes)

    def extension(self,i,p,new):
        h,r,t=p
        if i==0:
            return self.p.row('s3',(0,),(1,2,3),dict(enumerate((*p,new))))
        far=t if i==1 else r
        arm=p[i]
        return self.p.row('p4',(0,1),(2,3),dict(enumerate((h,arm,far,new))))


def replace(t,i,x): return t[:i]+(x,)+t[i+1:]
def reorder(p): return p[1],p[0],p[2]


class H3(ext.Pass):
    def __init__(self,e,a,upper,lo,hi):
        self.e,self.a,self.upper=e,a,upper
        self.disc=lo if upper else hi
        self.scale=e.snap.scale
        self.rs,self.ts=ext.states(a.receiver,a.rstock),ext.states(a.donor,a.tstock)
        self.ri,self.rm=ext.axes(self.rs); self.ti,self.tm=ext.axes(self.ts)
        self.nt=len(self.ts); self.n=len(self.rs)*self.nt
        self.zero=[0]*(self.n+1); self.r0=self.ri[J]; self.t0=self.ti[J]
        self.path_count=self.paid=0
        self.reserve=a.reserve

    def run(self):
        # The mother of the middle is the designated old receiver. Each
        # promotion seeds one continuation independently of its siblings.
        a=self.a
        pnext={}; rnext={}; Z=[0]*(len(self.rs)+1)
        for d in range(a.host,1,-1):
            P={}; R={}
            for c in range(min(d,a.stock),0,-1):
                h=(d,c); pb=P.get(c+1,self.zero)
                P[c]=self.layer(h,pnext.get(c,self.zero),pb)
                R[c]=self.edge(h,self.rs,self.rm,rnext.get(c,Z),R.get(c+1,Z),pb,True)
            pnext,rnext=P,R
        plain=self.singleton(self.rs,self.rm,[rnext[1][self.ri.get((d,c+1),-1)] for d,c in self.rs])
        return [(ext.HI if self.upper else ext.LO)*Q(plain,GRID*self.scale)]

    def layer(self,h,pg,pb):
        if h[1]<2: return self.zero
        adj=((1,2),(0,),(0,)); adjlen=(2,1,1)
        out=[0]*(self.n+1); nt=self.nt
        for i in range(len(self.rs)-1,-1,-1):
            r=self.rs[i]; rg,rb=self.rm[i]
            for j in range(nt-1,-1,-1):
                t=self.ts[j]; tg,tb=self.tm[j]; k=i*nt+j
                gs=(pg[k],out[rg*nt+j] if rg>=0 else 0,out[i*nt+tg] if tg>=0 else 0)
                bs=(pb[k],out[rb*nt+j] if rb>=0 else 0,out[i*nt+tb] if tb>=0 else 0)
                p=(h,r,t)
                if h not in self.e.mids or r[0]<self.a.old or t[0]>self.a.young:
                    out[k]=self.step(p,adj,gs,bs); continue
                # Pay the true signed negative part, or unused nonnegative
                # production family-cap budget if explicitly requested.
                outward,slopes=self.e.triple(*p)
                a=(self.e.hcap if self.reserve else 0)-outward
                widths=tuple(-s*Q(c-n,2) for s,(_,c),n in zip(slopes,p,adjlen))
                lo,hi=a+sum(min(0,w) for w in widths),a+sum(max(0,w) for w in widths)
                if hi<=0:
                    assert not self.reserve or (hi==0 and lo==0)
                    out[k]=self.step(p,adj,gs,bs); continue
                if self.reserve: assert lo>=0
                d=lcm(2,*(x[0] for x in p)); terms=[]
                for z,((dd,c),ns,v,b,w) in enumerate(zip(p,adj,gs,bs,widths)):
                    L=dd-c; width=(c-len(ns))*(d//2)
                    shift=self.round(w/width) if width else 0
                    terms.append((L*d+sum(d//p[y][0] for y in ns),width,v+shift,v))
                    if L: terms.append((L*(d//dd),0,b,b))
                val=affine(a,terms,d,self.disc,self.upper)
                if lo<0 and not self.reserve:
                    zero=self.step(p,adj,gs,bs)
                    # max(0, affine) on payoff. Sup over ports commutes
                    # with max, while max of infs remains a valid floor.
                    # Under a floor the affine subreward may be negative.
                    val=max(zero,val)
                assert val>=0
                out[k]=val; self.paid+=1
        self.path_count+=self.n
        return out


def affine(reward,terms,d,discount,upper):
    """Piecewise rational adverse corner with reward and gain on same port."""
    q=discount.denominator
    x=reward*GRID*d*q
    a=(-(-x.numerator//x.denominator) if upper else x.numerator//x.denominator)
    b=discount.numerator*d
    flex=[]
    for rate,w,shifted,raw in terms:
        a+=q*rate*raw; b+=q*rate
        if w: flex.append((shifted,w*q))
    aa,bb=a,b
    for v,w in sorted(flex,reverse=upper):
        if (v*b<=a if upper else v*b>=a): break
        a+=w*v; b+=w
    v=-(-a//b) if upper else a//b
    gap=aa-bb*v+sum(w*(max if upper else min)(0,z-v) for z,w in flex)
    assert gap<=0 if upper else gap>=0
    return v


def finite(e,limit,wing,arm,upper,lo,hi):
    """All arities in the closed square plus a rectangular medium/young wing."""
    disc=lo if upper else hi
    marks=ext.states(wing,wing)[::-1]; low=ext.states(arm,arm)[::-1]
    small=ext.states(limit,limit)[::-1]
    V={}; count=0; peak=Q(0)
    def val(x,y): return V.get((max(x,y),min(x,y)),0)
    for x in marks:
        if x[0]>limit and not low: continue
        for y in low if x[0]>limit else small:
            if y>x: continue
            a,b=e.edge(x,y)
            p=(x,y); d=lcm(2,x[0],y[0]); terms=[]
            for i,(n,c) in enumerate(p):
                w=(c-1)*(d//2); j=b[i]*Q(c-1,2)
                k=j*GRID/w if w else Q(0)
                sh=-(-k.numerator//k.denominator) if upper else k.numerator//k.denominator
                g=val((n+1,c),p[1-i]); bb=val((n,c+1),p[1-i]); L=n-c
                terms.extend(((L*d+d//p[1-i][0],w,g+sh,g),(L*(d//n),0,bb,bb)))
            lower=a+sum(min(0,z*Q(c-1,2)) for z,(_,c) in zip(b,p))
            assert lower>=0,(p,lower)
            peak=max(peak,a+sum(max(0,z*Q(c-1,2)) for z,(_,c) in zip(b,p)))
            V[x,y]=affine(a,terms,d,disc,upper); count+=1
    W={}
    for d,c in marks:
        L=d-c; bb=W.get((d,c+1),0)+val((d,c+1),J)
        g=W.get((d+1,c),0)
        W[d,c]=ext.solve(0,[(L*2*d,c*d,g),(2*L,0,bb)],2*d,disc,upper)
    v=(hi if upper else lo)*Q(W[J],GRID)
    print('EDGE',count,'peak ticks',float(peak),'ticks/core',float(v),'exact',v,'normalized',float(v/e.snap.scale),flush=True)
    return v


def fixtures(e,physical):
    counts=Counter()
    for D in ((38,3,3,8,2,4),(65,5,4,8,2,11),(79,5,5,8,3,9),
              (19,4,18,3,2,13),(11,4,5,8,2,7),(81,4,5,8,2,7),
              (252,4,5,8,2,7),(95,4,5,8,2,7)):
        edges=((0,1),(1,2),(0,3),(1,4),(2,5))
        adj={i:set() for i in range(len(D))}
        for x,y in edges: adj[x].add(y); adj[y].add(x)
        p={i:(d,len(adj[i])) for i,d in enumerate(D)}
        pre=physical.subsets(e.p,adj,p)
        n=len(D)
        plus=[physical.subsets(e.p,*physical.grow(adj,p,x,None)) for x in adj]
        birth=[physical.subsets(e.p,*physical.grow(adj,p,x,n))
               if D[x]>len(adj[x]) else {} for x in adj]
        def extract(rows,k,ids):
            return rows.get((k,tuple(sorted(ids))),(0,))[0]
        for (k,ids),(v,inner) in pre.items():
            if k not in ('g','h3'): continue
            S=set(ids); val=-e.snap.lam*v; gaps=Q(0)
            for x in S:
                d,c=p[x]; L=d-c
                j=extract(plus[x],k,S)-v
                port=adj[x]-S
                val+=(L+sum(Q(1,D[y]) for y in adj[x]))*j
                if L:
                    next_size=('h3',) if k=='g' else ('p4','s3')
                    val+=Q(L,d)*(extract(birth[x],k,S)-v+
                                      sum(extract(birth[x],t,S|{n}) for t in next_size))
                center=e.snap.center(d,c)
                if k=='g':
                    y,=S-{x}; iy=Q(1,D[y]);
                    if center is not None:
                        val+=(iy-center)*(extract(plus[x],'f',{x})-extract(pre,'f',{x}))
                    if center is not None and p[x] not in e.mids and c>=2:
                        for w in port:
                            z=Q(1,D[w]); val+=e.fy(p[x],iy)-e.fy(p[x],z)
                            gap=max(e.sig*((knot-z)*j-e.fy(p[x],knot)+e.fy(p[x],z))
                                    for knot in e.knots)
                            assert gap>=0
                            gaps+=gap; counts['FY gaps']+=1
                for y in port:
                    ext=('h3',) if k=='g' else ('p4','s3')
                    if any((t,tuple(sorted(S|{y}))) in pre for t in ext):
                        val-=(Q(1,D[y])-center)*j
                        counts[k+' outgoing']+=1
            if k=='h3':
                h,=inner; r,t=sorted(S-{h}); z=e.snap.center(*p[h])
                for w,y in ((r,t),(t,r)):
                    pair={h,y}
                    val+=(Q(1,D[w])-z)*(extract(plus[h],'g',pair)-extract(pre,'g',pair))
                a,b=e.triple(*(p[x] for x in (h,r,t)))
                assert a+sum(bb*sum(Q(1,D[y]) for y in adj[x]-S) for bb,x in zip(b,(h,r,t)))==e.sig*val
                assert e.sig*val<=e.hcap
            else:
                x,y=ids; a,b=e.edge(p[x],p[y]);
                assert a+sum(bb*sum(Q(1,D[w]) for w in adj[u]-S) for bb,u in zip(b,ids))==e.cap-e.sig*val-gaps
                assert e.cap-e.sig*val-gaps>=0
            counts[k]+=1
    assert all(counts.values()) and counts['FY gaps'] and counts['h3 outgoing']
    # Extra births preserve every older parent-child registration. Compare
    # seeded paths with physical sets on an intentionally branching tree.
    parents=(None,0,1,1,2,2,1,3,3,5,2,6,7,1)
    pair=set(); triples=set()
    for x,p in enumerate(parents[1:],1):
        if p:
            pair.add((p,x))
        triples.update((r,h,x) for r,h in pair if h==p)
    expected={(parents[p],p,x) for x,p in enumerate(parents[1:],1) if p and parents[p]}
    assert triples==expected and len(triples)==len({frozenset(t) for t in triples})
    print('literal physical checks',dict(counts),'distinct registered paths',len(triples),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--side',choices=('lower','upper'),default='lower')
    ap.add_argument('--kind',choices=('g','h3','threshold'),required=True)
    ap.add_argument('--ceiling',action='store_true')
    ap.add_argument('--limit',type=int,default=12)
    ap.add_argument('--wing',type=int,default=32)
    ap.add_argument('--arm',type=int,default=4)
    ap.add_argument('--host',type=int,default=42)
    ap.add_argument('--stock',type=int,default=6)
    ap.add_argument('--receiver',type=int,default=55)
    ap.add_argument('--rstock',type=int,default=4)
    ap.add_argument('--donor',type=int,default=7)
    ap.add_argument('--tstock',type=int,default=4)
    ap.add_argument('--old',type=int,default=19)
    ap.add_argument('--young',type=int,default=7)
    ap.add_argument('--reserve',action='store_true')
    ap.add_argument('--fixtures',action='store_true')
    a=ap.parse_args(); start=process_time()
    if a.kind=='threshold':
        assert not a.ceiling and not a.reserve
        a.host,a.stock,a.receiver,a.rstock=5,4,79,8
        a.donor,a.tstock,a.old,a.young=18,3,19,18
        a.fixtures=True
    base=gold.PRIVATE
    assert (base/'hashes.log').read_text().splitlines()==['quintet_lower.bin: OK','quintet_upper.bin: OK']
    assert subprocess.check_output(['git','-C',str(base/'source'),'rev-parse','HEAD'],text=True).strip()=='b435b4c824e2ba5835f58b85e67a27de7d6a0d49'
    raw=(base/'source/round2/neon2/quintet-cert/certificate_gold2.json').read_bytes()
    assert sha256(raw).hexdigest()=='5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
    docs=json.loads(raw); lo,hi=(Q(docs[s]['certified_exact']) for s in ('lower','upper'))
    ext.LO,ext.HI=lo,hi
    wire=HERE.parents[1]/'attacks/neon2-quintet/gold2/physical/check.py'
    sp=spec_from_file_location('gold2_physical',wire); mod=module_from_spec(sp); sp.loader.exec_module(mod)
    p=mod.Gold(base/f'quintet_{a.side}.bin')
    try:
        with Snapshot(a.side) as s:
            ident=lambda st:(st.st_dev,st.st_ino,st.st_size,st.st_ctime_ns,st.st_mtime_ns)
            assert ident(p.stat)==ident(s.stat)
            e=Rows(s,p,docs[a.side]); assert s.lam==Q(docs[a.side]['nominal'])
            print('original Gold2',vars(a),flush=True)
            if a.fixtures: fixtures(e,mod.old)
            if a.kind=='g':
                assert 2<=a.arm<=a.limit<=a.wing<s.d2
                finite(e,a.limit,a.wing,a.arm,a.ceiling,lo,hi)
            else:
                assert 2<=a.young<a.old<=a.receiver and a.young==a.donor
                vals=H3(e,a,a.ceiling,lo,hi); total=sum(vals.run())
                if a.kind=='threshold' and a.side=='lower':
                    assert total==Q(1084423030952180046412413909,1704710415481443450880000000000000)
                print('H3 normalized',total,float(total),'ticks/core',float(total*s.scale),
                      'states',vals.path_count,'paid rows',vals.paid,flush=True)
            s.check_unchanged()
            assert ident(os.fstat(p.fd))==ident(p.stat)
    finally:
        os.close(p.fd)
    print('CPU seconds',process_time()-start,flush=True)


if __name__=='__main__': main()
