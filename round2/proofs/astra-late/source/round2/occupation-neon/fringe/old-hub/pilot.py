#!/usr/bin/env python3
"""Conditional first-arm occupation; integrates all other arms as one port."""
import argparse
from collections import defaultdict
from functools import cache
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from charges import Cert, Charges, NB, Q, embedding, shape
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'sol-kernel-verifier'))
import verify as v

LO,HI=Q('.565676634491'),Q('.565771874464')
Z=Q(0); GRID=10**21; J=(2,1,()); NAMES=('hV','P4','Yt')
DEN=(Q('1.356833021311'),Q('1.367169995210'))


def down(x):return Q(x*GRID//1,GRID)


def arms(cap):
    # One dummy mark throughout. v.row's physical intrinsic law never
    # reads marks; its copy/kill fields are NOT used in this physical DP.
    def plain(x):return (2,x[1],tuple(sorted(plain(y) for y in x[2])))
    seen={J}; todo=[J]; transitions={}
    while todo:
        x=todo.pop()
        q=defaultdict(Q)
        for y,r in v.row(x)[3].items():q[plain(y)]+=r
        assert sum(q.values())==v.size(x)-Q(1,v.degree(x))
        transitions[x]=q
        for y in q:
            if y not in seen and v.size(y)<=cap:seen.add(y);todo.append(y)
    return sorted(seen,key=lambda t:(-v.size(t),t)),transitions


class Toll:
    def __init__(self,c,xs):
        self.c=c; self.stencils={}
        for x in xs:
            p,s=shape(x,1); ts=[]
            for f,ids in embedding((-1,0,*(i+1 for i in p[1:]))):
                if f not in NAMES or ids[0]!=0:continue
                xx=tuple(s[i-1] for i in ids[1:])
                if NB not in xx:ts.append((f,ids,xx))
            self.stencils[x]=(p,s,ts)

    @cache
    def row(self,d,k,x):
        p,s,ts=self.stencils[x]
        reward=[[Z,Z] for _ in NAMES]
        for f,ids,xx in ts:
            if max(map(sum,xx))>d:continue
            st=((d-k-1,k,1),*xx)
            if not self.c.supported(f,st):continue
            a,free=self.c.terms(f,st); slope=None
            for j,par,num,coef in free:
                if j==0:
                    assert num==(1 if par else k-1)
                    # Tail apex has no outgoing transfer in these families.
                    # Parent and siblings each multiply the same derivative.
                    if slope is not None:assert slope==coef
                    slope=coef
                else:
                    ii=ids[j]-1
                    ds=[p[ii]] if par else [i for i,r in enumerate(p) if r==ii and i+1 not in ids]
                    assert len(ds)==num and min(ds,default=0)>=0
                    a+=coef*sum(Q(1,sum(s[i])) for i in ds)
            assert slope is not None
            q=reward[NAMES.index(f)]; q[0]+=a;q[1]+=slope
        return tuple(tuple(r) for r in reward)


def host(hi,cap,kcap,pay):
    # Immediate toll is conditional on the birth of one arm. All actual
    # births get one such toll, including two births in the same d-bin.
    nxt={}
    for d in range(cap,1,-1):
        cur={}
        for k in range(min(kcap,d-1),-1,-1):
            L=d-k-1; b=Q(L,d)
            p=b*(cur.get(k+1,Z)+pay.get((d,k+1),Z))
            cur[k]=down(min((p+g*nxt.get(k,Z))/(hi+b+g) for g in (L,L+Q(k+1,2))))
        nxt=cur
    return LO*nxt[0]


def solve(c,cap,kcap,acap,details=False):
    xs,arcs=arms(acap); t=Toll(c,xs)
    # An arm unable to support any future motif before the cut has zero
    # continuation. Work backwards once on the arm DAG to prove that.
    seen=set()
    for x in xs:
        ax=v.plus(x)
        if t.stencils[x][2] or ax in seen or any(y in seen for y in arcs[x]):seen.add(x)
    xs=[x for x in xs if x in seen]
    assert J in xs
    z=(Z,)*4; nxt={}; pay=[{} for _ in z]; count=0
    for d in range(cap,1,-1):
        cur={}
        for k in range(min(kcap,d-1),0,-1):
            L=d-k-1; p=Q(1,d); b=L*p
            for x in xs:
                rates=arcs[x]; ad=v.degree(x); base=HI+b+sum(rates.values())+p+L+Q(1,ad)
                a=cur.get((k,v.plus(x)),z)
                w=cur.get((k+1,x),z); grow=nxt.get((k,x),z)
                r=t.row(d,k,x) if d>=65 and k>=5 else ((Z,Z),)*3
                rr=(tuple(map(sum,zip(*r))),*r) if details else (tuple(map(sum,zip(*r))),)
                closed=[b*w[j]+p*a[j] for j in range(len(rr))]
                for y,rate in rates.items():
                    ar=cur.get((k,y),z)
                    for j in range(len(rr)):closed[j]+=rate*ar[j]
                hi=Q(k,2)
                def bell(pair,q,gg):
                    u0,u1=pair; pts=[Z,hi]
                    if u1 and 0<-u0/u1<hi:pts.append(-u0/u1)
                    return down(min((max(Z,u0+u*u1)+q+u*gg)/(base+u) for u in pts))
                cur[k,x]=tuple(bell(pair,q,gg) for pair,q,gg in zip(rr,closed,grow))
                count+=1
            for table,f in zip(pay,cur[k,J]):table[d,k]=f
        nxt=cur
    return [host(HI,cap,kcap,p) for p in pay if p],len(xs),count


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blobs',nargs=2,type=Path)
    p.add_argument('--degree',type=int,default=100);p.add_argument('--arity',type=int,default=10)
    p.add_argument('--arm',type=int,default=6);p.add_argument('--detail',action='store_true')
    a=p.parse_args()
    if a.degree<65 or a.arity<5 or a.arm<6:p.error('requires degree >=65, arity >=5 and arm size >=6')
    endpoints=[]
    for i,blob in enumerate(a.blobs):
        c=Charges(Cert(blob))
        assert c.sign==(1,-1)[i] and c.c.lam==(Q('.565673'),Q('.565775'))[i] and c.c.D==18
        val,n,r=solve(c,a.degree,a.arity,a.arm,a.detail)
        reward=max(val[0],sum(val[1:]))
        assert reward>0
        gain=reward/DEN[i]
        print(c.sign,'tracked arm geometries',n,'Bellman rows',r,'lower per core',
              *[v.bracket(x,x,23) for x in val], 'inward move',v.bracket(gain,gain,23),flush=True)
        endpoints.append((LO,HI)[i]+c.sign*gain)
        Toll.row.cache_clear()
    print('outward',v.bracket(*endpoints,16))


if __name__=='__main__':main()
