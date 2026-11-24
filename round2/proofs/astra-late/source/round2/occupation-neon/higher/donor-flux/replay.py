#!/usr/bin/env python3
"""Exact directed donor transposition and finite-snapshot flux diagnostics."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
import hashlib
from pathlib import Path
from random import Random
from importlib.util import module_from_spec, spec_from_file_location
_spec=spec_from_file_location('older_local',Path(__file__).resolve().parents[1]/'older'/'replay.py')
_older=module_from_spec(_spec); _spec.loader.exec_module(_older)
Cert,DEN,Local,NB,PAYLOAD_SHA,SPEC,bud,label,pick,sample,up=(getattr(_older,s) for s in (
    'Cert','DEN','Local','NB','PAYLOAD_SHA','SPEC','bud','label','pick','sample','up'))

TARGET=('hV','P4','Yt')


def joint_box(total,slots):
    """Sharp affine endpoints at fixed total; slots=(signed price,lo,hi)."""
    assert sum(l for _,l,h in slots)<=total<=sum(h for _,l,h in slots)
    def fill(reverse):
        left=total-sum(l for _,l,h in slots); value=sum(c*l for c,l,h in slots)
        for c,l,h in sorted(slots,reverse=reverse):
            used=min(left,h-l); value+=c*used; left-=used
        assert left==0
        return value
    return fill(False),fill(True)


class Flux:
    def __init__(self,E):self.E=E

    @cache
    def slopes(self,f,x):
        E=self.E; t=SPEC[f]; ds=E.terms(f,x)[1]
        out=[]
        for i,(_,k,h) in enumerate(x):
            for parent,n in ((True,h-int(i!=0)),(False,k-t.count(i))):
                if n:
                    shift=E.transfer(f,x,i) if E.moved(f,x,i,parent) else 0
                    out.append((i,parent,n,E.sign*(ds[i]-shift)/E.c.S))
        return tuple(out)

    def snapshot(self,p,x,raw,ch,add):
        by_receiver=Q(0); donor_load=defaultdict(Q); joint=defaultdict(Q)
        sz=[1+s[0] for s in x]; nc=[1]*len(x)
        for j in range(len(x)-1,0,-1):sz[p[j]]+=sz[j];nc[p[j]]+=nc[j]
        for f,ids,row in raw:
            a=ids[0]; D=sum(x[a]); k=x[a][1]
            if f not in TARGET or not a or nc[a]<5 or sz[a]<=12:continue
            if any(x[j]==NB for j in ids):continue
            region='13-17' if 13<=D<18 and k>=3 else '18-64:k3-4' if 18<=D<=64 and 3<=k<=4 else (
                '18-64:k5+' if 18<=D<=64 and k>=5 else '65+:k5+' if D>=65 and k>=5 else 'other')
            y=pick(x,ids); base=self.E.charge(f,y,[(0,0)]*len(ids)); predicted=base
            lo=base; hi=base; total=Q(0); slots=[]; residual=[]; low_degree_value=Q(0); low_degree_total=Q(0)
            for i,parent,n,d in self.slopes(f,y):
                j=ids[i]; donors=([p[j]] if parent else [v for v in ch[j] if v not in ids])
                assert len(donors)==n
                r=sum((Q(1,sum(x[v])) for v in donors),Q(0)); predicted+=d*r
                lo+=min(0,d*n/2); hi+=max(0,d*n/2)
                total+=r; slots.append((d,Q(0),Q(n,2)))
                small=[v for v in donors if not parent and sum(x[v])<=4]
                retained=sum((Q(1,sum(x[v])) for v in small),Q(0))
                low_degree_total+=retained; low_degree_value+=d*retained
                residual.append((d,Q(0),Q(n-len(small),2 if parent else 5)))
                if parent:continue
                # Upper/lower slots keep the actual single signed coefficient.
                # No donor is used twice for a fixed role/embedding.
                key=(region,f,i,1 if d>=0 else -1)
                by_receiver+=d*r
                for v in donors:
                    wt=Q(1,sum(x[v])); donor_load[v]+=d
                    add(key,y,n,d,wt,sum(x[v]),x[v]==NB)
            assert predicted==row and 0<=lo<=row<=hi,(f,y,predicted,row)
            a,b=joint_box(total,slots)
            c,d=joint_box(total-low_degree_total,residual)
            c+=low_degree_value;d+=low_degree_value
            assert lo<=base+a<=base+c<=row<=base+d<=base+b<=hi
            joint[f,'independent']+=hi-lo
            joint[f,'joint total']+=b-a
            joint[f,'joint and cheap edges']+=d-c
        by_donor=sum((v/sum(x[j]) for j,v in donor_load.items()),Q(0))
        assert by_receiver==by_donor
        return joint


def mask_on(x,receiver,cut):
    return int(sum(x[receiver])>cut)


def cross_control(F):
    """Adjoint gain identity including newly born embedding and transverse seams."""
    E=F.E; f='P4'; t=SPEC[f]; i=1; cut=4
    x=((8,1,1),(1,2,1),(0,1,1),(1,0,1)); assert sum(x[i])==cut
    # A receiver-only threshold gives the exact crossing identity. Refine
    # by a compact transverse aperture, which introduces signed seam terms.
    aperture=lambda y:int(sum(y[2])==2 and y[0][1]==1)
    phi=lambda y:mask_on(y,i,cut)*aperture(y)
    for y in (x,label(x,i,up(x[i])),label(x,2,up(x[2]))):
        rho=(Q(1,7),Q(1,5),Q(0),Q(1,2)); g=[]; b=[]
        for j,(l,k,h) in enumerate(y):
            g.append(l+sum(Q(1,sum(y[w])) for w in range(len(y)) if t[w]==j or t[j]==w)+rho[j])
            b.append(Q(l,l+k+h))
        gain=sum(g[j]*(phi(label(y,j,up(z)))-phi(y)) for j,z in enumerate(y))
        prom=sum(b[j]*(phi(label(y,j,bud(z)))-phi(y)) for j,z in enumerate(y) if z[0])
        seam=sum(g[j]*(phi(label(y,j,up(z)))-phi(y)) for j,z in enumerate(y) if j!=i)
        crossing=aperture(y)*int(sum(y[i])==cut)
        inside=g[i]-rho[i]
        assert gain+prom-inside*crossing-seam-prom==rho[i]*crossing
    print('exact projected crossing/side-seam controls 3',flush=True)
    # Total incoming at a receiver can come from its parent or extra child.
    # Shift stocks in the omitted vertices; both actual total receives agree,
    # while one-edge downward/child charge remains different.
    slots=[(Q(0),Q(0),Q(1,2)),(Q(1,1000),Q(0),Q(1,2))]
    a,b=joint_box(Q(3,4),slots)
    assert (a,b)==(Q(1,4000),Q(1,2000))
    assert a==Q(1,1000)*Q(1,4) and b==Q(1,1000)*Q(1,2)
    print('joint total alone: physical D_parent,D_child=(2,4) vs (4,2), width',str(b-a),flush=True)


def coarsened_mark_control():
    # Receiver's two incoming sources have the same total in each mode;
    # a shared donor edge has a single joint endpoint degree. Neither
    # marginal residence nor crossing counts encode its pairing.
    lower=upper=Q(0)
    for weight,c in ((Q(1,100),Q(3,1000)),(Q(1,200),-Q(3,1000))):
        a,b=joint_box(Q(3,4),[(c,Q(0),Q(1,2)),(Q(0),Q(0),Q(1,2))])
        lower+=weight*a; upper+=weight*b
    assert lower==0 and upper==Q(9,800000)
    print('fabricated signed paired bins, exact total pressure, residual interval',str(lower),str(upper),flush=True)


def run(E,side,n,reps,seed):
    F=Flux(E); cross_control(F); coarsened_mark_control(); rng=Random(seed); num=0
    D=defaultdict(lambda:[Q(0),Q(0),Q(0),Q(0),0,0]); desc=defaultdict(lambda:[Q(0),Q(0),Q(0),0]);
    reserve=defaultdict(lambda:[Q(0),Q(0),Q(0),0,0]); widths=defaultdict(Q)
    def take(key,y,m,d,w,degree,virgin):
        # positive and negative signs stay distinct until after integration
        q=D[key]; q[0]+=w;q[1]+=abs(d)*w;q[2]+=abs(d)*w if degree<=4 else 0
        q[3]+=abs(d)*w if virgin else 0;q[4]+=1;q[5]+=degree<=4
        z=desc[key,y];z[0]+=abs(d)*w;z[1]+=w;z[2]=abs(d);z[3]=m
        # Conservative complement keeps the independent external corners.
        # For a negative slope only its shortfall from half can be booked.
        e=abs(d)*(w if d>=0 else Q(1,2)-w)
        k=key[:2]; p=reserve[k]; p[0]+=e
        if d>0 and d>=Q(1,1000) and degree<=4:p[1]+=e;p[2]+=w;p[3]+=1
        if d>0 and d>=Q(1,1000):p[4]+=1
    for _ in range(reps):
        p,x=sample(n,rng); raw,ch=E.tree(p,x);num+=len(x)
        for k,v in F.snapshot(p,x,raw,ch,take).items():widths[k]+=v
    print('sample transpose/affine equalities passed trees',reps,'n',n,'seed',seed,'cores',num,flush=True)
    for key,v in sorted(D.items(),key=lambda kv:kv[1][1],reverse=True)[:18]:
        m,price,cheap,virgin,cnt,cnt4=v
        print('sample flux',key,'per core',float(m/num),'absolute price',float(price/num),
              'D<=4 share',float(cheap/price) if price else 0,'virgin share',float(virgin/price) if price else 0,
              'D<=4 uniform edge share',cnt4/cnt,flush=True)
    target=DEN[side]/10**5
    for (key,x),(price,w,c,m) in sorted(desc.items(),key=lambda kv:kv[1][0],reverse=True)[:9]:
        if c==0:continue
        print('fixed',key,x,'coeff',float(c),'open slots',m,'sample price/core',float(price/num),
              'rate mass for +1e-6 move',float(target/(10*c)),'for +1e-5 move',float(target/c),
              'max price/occurrence',float(c*m/2),flush=True)
    for f in TARGET:
        v=[sum((v[1] for key,v in D.items() if key[1]==f and key[-1]==sign),Q(0))/num for sign in (1,-1)]
        print('sample all-region signed child pressure',f,*(float(x) for x in (*v,v[0]-v[1])),flush=True)
        print('oracle sample width/core',f,*(float(widths[f,k]/num) for k in
              ('independent','joint total','joint and cheap edges')),flush=True)
    for k,v in sorted(reserve.items(),key=lambda kv:kv[1][0],reverse=True):
        if v[0]:print('sample nonnegative donor excess',k,'all',float(v[0]/num),'cheap strong',float(v[1]/num),
                      'cheap strong reciprocal mass',float(v[2]/num),'cheap/strong edge hits',v[3:5],flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blob',type=Path)
    p.add_argument('--n',type=int,default=16000);p.add_argument('--reps',type=int,default=4)
    p.add_argument('--seed',type=int,default=122419);p.add_argument('--verify-hash',action='store_true');a=p.parse_args()
    assert __debug__ and a.reps>=1 and a.n>=3
    c=Cert(a.blob);side='upper' if c.mode==1 else 'lower';assert c.D==18 and c.S==10**11
    if a.verify_hash:
        with a.blob.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        assert digest==PAYLOAD_SHA[side];print('pinned',side,digest,flush=True)
    run(Local(c),side,a.n,a.reps,a.seed)


if __name__=='__main__':main()
