#!/usr/bin/env python3
"""Original Gold2 upper S4, actual mother and three ordered hub children."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, chain, product
from math import lcm
from pathlib import Path
from time import process_time
import json
import os

HERE=Path(__file__).resolve().parent


def load(name,path):
    sp=spec_from_file_location(name,path); m=module_from_spec(sp); sp.loader.exec_module(m)
    return m


base=load('s4_base',HERE.parent/'post-fy/replay.py')
J=base.J; GRID=base.GRID


def change(p,i,x): return p[:i]+(x,)+p[i+1:]


def port_shift(slope):
    # w/d is the physical rate width. The Bellman numerator multiplies
    # BOTH reward and transition rates by d, so their ratio has no d.
    x=Q(slope)*GRID
    return x.numerator//x.denominator


def affine_fixture():
    # A fractional discount and nontrivial common denominator expose both
    # slope signs. The exact physical corner ratios do not use base.affine.
    for d in (6,30,420):
        a=Q(301,3); disc=Q(13,23); slopes=(Q(-129,7),Q(41,11))
        g=(6*GRID,2*GRID); rate=(Q(1,3),Q(5,6)); width=(Q(1,2),Q(1))
        terms=[(int(r*d),int(w*d),v+port_shift(s),v) for r,w,v,s in zip(rate,width,g,slopes)]
        exact=min((GRID*(a+sum(z*s for z,s in zip(x,slopes)))+sum((r+z)*v for r,z,v in zip(rate,x,g))) /
                  (disc+sum(rate)+sum(x)) for x in product(*((0,w) for w in width)))
        result=base.affine(a,terms,d,disc,False)
        assert exact-2<result<=exact,(d,result,exact)
        wrong=[(r,w,v+port_shift(s/d),v) for (r,w,_,v),s in zip(terms,slopes)]
        assert base.affine(a,wrong,d,disc,False)>exact
    print('affine scaling: 12 exact physical corners, omitted-factor control rejected',flush=True)


class Row:
    def __init__(self,p,s,receipt):
        self.p,self.s=p,s
        self.cap=receipt['verify']['caps']['s4']; assert self.cap==50

    def val(self,p):
        return self.value('s4',p[0],tuple(sorted(p[1:])))

    @lru_cache(maxsize=120000)
    def value(self,key,hub,tips):
        i=self.p.ins[key.upper()].get(hub)
        if i is None: return 0
        return self.p.value(key,(i,),tips)

    def signed(self,p):
        """Upper signed row, affine in each exterior reciprocal sum."""
        v=self.val(p); a=-self.s.lam*v; b=[]
        for i,(D,c) in enumerate(p):
            adj=range(1,5) if i==0 else (0,)
            adj=tuple(adj); stock=D-c
            dz=self.val(change(p,i,(D+1,c)))-v
            a+=(stock+sum(Q(1,p[j][0]) for j in adj))*dz
            if stock:
                a+=Q(stock,D)*(self.val(change(p,i,(D,c+1)))-v)
            b.append(dz if c>len(adj) else 0)
        D,c=p[0]; center=self.s.center(D,c)
        for i in range(1,5):
            # Physical deletion once per label, even where the remaining
            # endpoint tuples are numerically identical.
            q=p[:i]+p[i+1:]
            tips=tuple(sorted(q[1:]))
            a+=(Q(1,p[i][0])-center)*(self.value('s3',(D+1,c),tips)-self.value('s3',(D,c),tips))
        return a,tuple(b)

    def global_cap_check(self,hub):
        """Fixture check against production's largest local cap; no rebooking."""
        assert hub in self.p.ins['S4']
        return self.cap


def fixture(r,physical):
    affine_fixture()
    checked=repeated=exterior=0
    for D in ((8,2,2,2,2,2),(10,14,14,14,14,14),(9,37,2,2,2,5),(10,180,2,14,13,3),
              (9,37,2,2,3,5,13,40,3)):
        adj={0:set(range(1,6)),**{i:{0} for i in range(1,6)}}
        if len(D)>6:
            for i,j in ((1,6),(1,7),(4,8)):
                adj[i].add(j); adj[j]={i}
        p={i:(d,len(adj[i])) for i,d in enumerate(D)}
        newborn=len(D)
        pre=physical.subsets(r.p,adj,p)
        up=[physical.subsets(r.p,*physical.grow(adj,p,i,None)) for i in adj]
        nb=[physical.subsets(r.p,*physical.grow(adj,p,i,newborn)) for i in adj]
        def val(rows,k,S): return rows.get((k,tuple(sorted(S))),(0,))[0]
        for arms in combinations(range(1,6),4):
            S={0,*arms}; old=val(pre,'s4',S); z=-r.s.lam*old
            for i in S:
                d,c=p[i]; L=d-c
                z+=(L+sum(Q(1,D[j]) for j in adj[i]))*(val(up[i],'s4',S)-old)
                z+=Q(L,d)*(val(nb[i],'s4',S)-old)
            terms=[]
            for i in arms:
                T=S-{i}; j=val(up[0],'s3',T)-val(pre,'s3',T)
                terms.append((Q(1,D[i])-r.s.center(*p[0]))*j)
            order=(0,*arms); a,b=r.signed(tuple(p[i] for i in order))
            assert a+sum(v*sum(Q(1,D[j]) for j in adj[i]-S) for i,v in zip(order,b))==z+sum(terms)<=r.global_cap_check(p[0])
            exterior+=any(v and adj[i]-S for i,v in zip(order[1:],b[1:]))
            repeated+=len(set(terms))==1 and bool(terms[0]); checked+=1
    # A promotion grows all earlier registrations and adds each new subset
    # of that age exactly once; equal labels do not determine identity.
    subsets=[{()}]+[set() for _ in range(3)]
    for child in range(7):
        for k in (3,2,1): subsets[k].update(x+(child,) for x in subsets[k-1])
    assert subsets[3]==set(combinations(range(7),3)) and repeated and exterior
    print('fixtures S4 rows',checked,'four equal nonzero arrows',repeated,'nonzero exterior tip rows',exterior,
          'unique crowns',len(subsets[3]),flush=True)


def states(d,c):
    return tuple((n,k) for n in range(d,1,-1) for k in range(min(n,c),0,-1))


def price(r,a,lo,hi):
    """Sorted child states quotient storage only, not physical promotions."""
    child=states(a.child,a.child_stock)
    young=tuple(x for x in child if x[0]<=a.young)
    older=tuple(x for x in child if x[0]>a.young)
    hub=states(a.hub,a.stock); mother=states(a.mother,a.mother_stock)
    childsets=[((),)]+[tuple(chain(((x,*cs) for x in older for cs in combinations_with_replacement(young,k-1)),
                                combinations_with_replacement(young,k))) for k in range(1,4)]
    # Child descriptors are sorted descending. Each child change strictly
    # increases the sorted tuple lexicographically and is already available.
    old_layers={}; root={}; counts=[0]*4; paid=0
    for m in mother:
        D,c=m; cur={}
        for h in hub:
            H,C=h
            for k in (3,2,1,0):
                if C<k+1: continue
                for cs in childsets[k]:
                    p=(m,h,*cs); terms=[]
                    scale=lcm(2,*(x[0] for x in p))
                    adj=((1,),tuple(i for i in range(k+2) if i!=1))+((1,),)*k
                    A=Q(0); slopes=(0,)*(k+2)
                    if k==3 and a.old<=D<=a.through and h in r.p.ins['S4'] and (not a.bare or H==C):
                        q=(h,m,*cs)
                        A,slopes=r.signed(q); A=-A
                        slopes=tuple(-slopes[i] for i in (1,0,2,3,4)); paid+=1
                    for i,((n,z),neigh) in enumerate(zip(p,adj)):
                        leaf=n-z; w=(z-len(neigh))*(scale//2)
                        def cont(x):
                            if i==0: return old_layers.get(x,{}).get((h,cs),0)
                            if i==1: return cur.get((x,cs),0)
                            kids=tuple(sorted(change(cs,i-2,x),reverse=True))
                            return cur.get((h,kids),0)
                        g=cont((n+1,z)); b=cont((n,z+1))
                        if i==1 and k<3:
                            b+=cur.get(((n,z+1),tuple(sorted((*cs,J),reverse=True))),0)
                        rate=leaf*scale+sum(scale//p[j][0] for j in neigh)
                        terms.append((rate,w,g))
                        if leaf:
                            terms.append((leaf*(scale//n),0,b))
                    out=base.ext.solve(0,terms,scale,hi,False)
                    if k==3 and (A or any(slopes)):
                        offsets=[port_shift(s) for s in slopes]
                        shifted=[]
                        for rate,w,g in terms:
                            # Each vertex contributes one gain followed by
                            # its nonzero promotion (if any).
                            shifted.append((rate,w,g,g))
                        j=0
                        for (n,z),shift in zip(p,offsets):
                            rate,w,_,g=shifted[j]; shifted[j]=(rate,w,g+shift,g)
                            j+=1+(n>z)
                        assert j==len(shifted)
                        out=max(out,base.affine(A,shifted,scale,hi,False))
                    cur[h,cs]=out; counts[k]+=1
        # Keep only successor degree slice and the current degree's arities.
        old_layers[m]=cur
        root[m]=base.ext.solve(0,[((D-c)*2*D,c*D,root.get((D+1,c),0)),
            (2*(D-c),0,root.get((D,c+1),0)+old_layers.get((D,c+1),{}).get((J,()),0))],2*D,hi,False)
        expired=(D+1,c)
        old_layers.pop(expired,None)
    result=lo*Q(root[J],GRID)
    print('prefix states',counts,'paid mask states',paid,'ticks/core',result,float(result),
          'normalized',result/r.s.scale,float(result/r.s.scale),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mother',type=int,default=48)
    ap.add_argument('--mother-stock',type=int,default=7)
    ap.add_argument('--hub',type=int,default=7)
    ap.add_argument('--stock',type=int,default=6)
    ap.add_argument('--child',type=int,default=3)
    ap.add_argument('--child-stock',type=int,default=2)
    ap.add_argument('--young',type=int,help='at most one included child may exceed this degree')
    ap.add_argument('--old',type=int,default=19)
    ap.add_argument('--through',type=int,help='last paid mother degree (defaults to continuation cap)')
    ap.add_argument('--bare',action='store_true',help='pay only when the hub has no ordinary leaves')
    ap.add_argument('--fixtures',action='store_true')
    ap.add_argument('--fixtures-only',action='store_true')
    a=ap.parse_args(); start=process_time(); a.through=a.through or a.mother; a.young=a.young or a.child
    assert 2<=a.old<=a.through<=a.mother and 4<=a.hub<=10 and 4<=a.stock<=6 and 2<=a.young<=a.child
    path=base.gold.PRIVATE
    assert (path/'hashes.log').read_text().splitlines()==['quintet_lower.bin: OK','quintet_upper.bin: OK']
    raw=(path/'source/round2/neon2/quintet-cert/certificate_gold2.json').read_bytes()
    assert sha256(raw).hexdigest()=='5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
    docs=json.loads(raw); lo,hi=(Q(docs[s]['certified_exact']) for s in ('lower','upper'))
    wire=load('s4_wire',HERE.parents[1]/'attacks/neon2-quintet/gold2/physical/check.py')
    p=wire.Gold(path/'quintet_upper.bin')
    try:
        with base.Snapshot('upper') as s:
            tag=lambda x:(x.st_dev,x.st_ino,x.st_size,x.st_mtime_ns,x.st_ctime_ns)
            assert tag(s.stat)==tag(p.stat) and s.lam==Q(docs['upper']['nominal'])
            r=Row(p,s,docs['upper']); print('original upper S4',vars(a),flush=True)
            if a.fixtures or a.fixtures_only: fixture(r,wire.old)
            if not a.fixtures_only: price(r,a,lo,hi)
            s.check_unchanged(); assert tag(os.fstat(p.fd))==tag(p.stat)
    finally: os.close(p.fd)
    print('CPU seconds',process_time()-start,flush=True)


if __name__=='__main__': main()
