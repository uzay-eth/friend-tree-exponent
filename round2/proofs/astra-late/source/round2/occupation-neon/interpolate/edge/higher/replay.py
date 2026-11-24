#!/usr/bin/env python3
"""Authenticated staged hV/P4/Yt correction above the frozen edge replay."""
import argparse
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
FROZEN='e90651922acb6b3a427872cf3801f4b9097d2b83c7abe0ab75add6e60c78350f'
assert sha256((HERE.parent/'replay.py').read_bytes()).hexdigest()==FROZEN
spec=spec_from_file_location('frozen_edge',HERE.parent/'replay.py')
E=module_from_spec(spec);spec.loader.exec_module(E)
from charges import SPEC

Q,Z,J=E.Q,E.ZERO,E.J
FAMILIES=('hV','P4','Yt')
EDGE_POINT={2:Q('.565677980686114'),1:Q('.565770513279953')}


def reps(x,i,t):return x[:i]+(t,)+x[i+1:]


def gate(x):
    # Connected planted fragment: each included child replaces one
    # unit of degree at its parent. Outside children add at least one
    # ordinary vertex and one core apiece to the full fringe.
    k=sum(t[1] for t in x)
    return sum(t[0] for t in x)-len(x)+1>12 and (k>=4 or max(t[0] for t in x)>=17)


def normal(f,x):
    if f=='Yt' and x[2]>x[3]:return x[:2]+(x[3],x[2])
    return x


def row(s,a,g,b,box,coef,ga,ba,round_up=False):
    # Affine reward; GMP nonnegative box checked in the extractor.
    a+=sum(t*z for t,z in zip(b,ba))+sum(t*z for t,z in zip(g,ga))
    den=s+sum(b)+sum(g)
    # For ratio R, a corner minimises iff u=max when coef+ga<R.
    # Sorting breakpoints visits every possible minimising corner.
    corners=sorted((c+z,u) for c,z,u in zip(coef,ga,box) if u)
    best=a/den;raw,base=a,den
    for v,u in corners:
        if v>=best:break
        a+=u*v;den+=u;best=a/den
    w=-E.D.down(-best) if round_up else E.D.down(best)
    assert w>=0
    if round_up:assert w*den>=a
    else:
        # Independent full-box minimum of the cross-multiplied residual.
        assert raw-base*w+sum(min(Z,u*(v-w)) for v,u in corners)>=0
    return w


def layer(c,s,f,dom,children,*,upper=False,replace=False):
    p=SPEC[f];adj=tuple(tuple(j for j in range(len(p)) if p[j]==i or p[i]==j) for i in range(len(p)))
    values={};paid=0;flow=0;by_k={}
    for x in product(*dom):
        if x!=normal(f,x):continue
        b=[Q(d-k-1,d) for d,k in x]
        box=tuple(Q(k-p.count(i)+int(i==0),2) for i,(_,k) in enumerate(x))
        g=tuple(d-k-1+sum(Q(1,x[j][0]) for j in adj[i]) for i,(d,k) in enumerate(x))
        gains=[E.moves(t)[0] for t in x];buds=[E.moves(t)[1] for t in x]
        ga=[values.get(normal(f,reps(x,i,y)),Z) for i,y in enumerate(gains)]
        ba=[values.get(normal(f,reps(x,i,y)),Z) for i,y in enumerate(buds)]
        if f=='hV':
            for i,f2 in ((1,'Yt'),(2,'P4')):
                z=children.get(f2,{}).get(normal(f2,reps(x,i,buds[i])+(J,)),Z)
                ba[i]+=z;flow+=bool(b[i] and z)
        a=Z;coef=(Z,)*len(x)
        if replace or gate(x):
            st=tuple(map(E.desc,x))
            if c.supported(f,st):
                a,free=c.terms(f,st);coeff={};count={}
                for i,_,n,v in free:
                    if i in coeff:assert coeff[i]==v
                    coeff[i]=v;count[i]=count.get(i,0)+n
                assert tuple(Q(count.get(i,0),2) for i in range(len(x)))==box
                coef=tuple(coeff.get(i,Z) for i in range(len(x)))
                worst=a+sum(min(Z,v*u) for v,u in zip(coef,box));assert worst>=0,(f,x,worst)
                if a or any(coef):
                    paid+=1
                    j=min(5,x[0][1]);by_k[j]=by_k.get(j,0)+1
        values[x]=row(s,a,g,b,box,coef,ga,ba,upper)
    print(f,'states',len(values),'paid',paid,'paid by apex k (5=5+)',sorted(by_k.items()),
          'live post-promotion sources',flow,'sup',float(max(values.values(),default=Z)),flush=True)
    return values


def edge(e,s,tri,*,upper=False):
    v={};births=0
    for p,w in product(e.p,e.w):
        pg,pp=E.moves(p);wg,wp=E.moves(w);pb,wb,pts=e.rows[p,w]
        birth=tri.get((p,wp,J),Z);births+=bool(wb and birth)
        a=pb*v.get((pp,w),Z)+wb*(v.get((p,wp),Z)+birth)
        vals=tuple((r+a+g*v.get((pg,w),Z)+h*v.get((p,wg),Z),s+pb+wb+g+h) for g,h,r in pts)
        best=min(x/y for x,y in vals)
        z=-E.D.down(-best) if upper else E.D.down(min(Q('.004'),best))
        assert (any if upper else all)(z*y>=x if upper else z*y<=x for x,y in vals)
        v[p,w]=z
    print('edge mapped hV births',births,flush=True)
    return v


def domains(cfg):
    d,k,m,mk,t,tk,d4,k4=cfg
    def zz(d,k,n=0):return E.states(d,n,k)
    return {'hV':(zz(d,k,1),zz(m,mk,1),zz(t,tk)),
            'P4':(zz(d4,k4,1),zz(min(m,5),mk,1),zz(min(t,5),tk,1),zz(min(t,5),tk)),
            'Yt':(zz(d4,k4,1),zz(min(m,5),mk,2),zz(min(t,5),tk),zz(min(t,5),tk))}


def report(c,cs,s,stage,v,H,*,upper=False,replace=False):
    r=cs.solve(s,v,upper=upper);m,_,z,tail=r
    h=max(c.sign*(q-z[J]) for q in (*z.values(),*tail.values()))
    if c.sign>0:h+=max(v.values())+H
    reserve=({1:Q('2.799995635483e-7'),2:Q('5.520306761551e-7')}[c.c.mode]
             if replace else E.RETAINED[c.c.mode])
    paid=m+reserve
    print('stage',stage,'s',s,'local margin',float(m),'composed numerator',float(paid),flush=True)
    if not upper and paid>=0:
        end=s+c.sign*paid/(E.DEN[c.c.mode]+h)
        print('stage',stage,'outward',E.D.show(end,15,upper=c.sign<0),
              'gain from frozen edge',float(c.sign*(end-EDGE_POINT[c.c.mode])),flush=True)
    return r


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('blob',type=Path);ap.add_argument('--s',type=Q)
    ap.add_argument('--config',type=int,nargs=8,default=(24,6,6,3,6,3,16,5))
    stage=ap.add_mutually_exclusive_group()
    stage.add_argument('--only-hv',action='store_true');stage.add_argument('--full-only',action='store_true')
    ap.add_argument('--robust-upper',action='store_true')
    ap.add_argument('--replace',action='store_true',help='admit every small row; retain only frozen high-parent block')
    a=ap.parse_args();assert not sys.flags.optimize
    for path,h in E.SOURCE_HASHES.items():assert sha256((HERE.parents[3]/path).read_bytes()).hexdigest()==h,path
    c=E.Charges(E.D.Cert(a.blob));s=a.s or EDGE_POINT[c.c.mode]+c.sign*Q('5e-6')
    assert c.sign*(s-c.c.lam)>=0 and s>0
    with a.blob.open('rb') as blob:
        h=sha256()
        for data in iter(lambda:blob.read(2**23),b''):h.update(data)
    assert h.hexdigest()==E.PAYLOAD_HASHES[c.c.mode]
    print('payload',h.hexdigest(),'config',a.config,flush=True)
    cs=E.Census(c.c);e=E.Edge(c,*( (40,8) if c.sign>0 else (32,7)),mask='disjoint')
    d,k,m,mk,t,tk,d4,k4=a.config
    assert 2<=d4<=d<=max(p[0] for p in e.p) and k4<=k<=max(p[1] for p in e.p)
    assert m<=max(w[0] for w in e.w) and mk<=max(w[1] for w in e.w)
    report(c,cs,s,'no higher',e.solve(s),Z)
    dom=domains(a.config)
    for high in (False,) if a.only_hv else (True,) if a.full_only else (False,True):
        v={f:layer(c,s,f,dom[f],{},replace=a.replace) for f in FAMILIES[1:]} if high else {}
        v['hV']=layer(c,s,'hV',dom['hV'],v,replace=a.replace)
        ev=edge(e,s,v['hV'])
        # Norm via terminal ownership for paths, and bounded middle arity
        # for forks; sum of the three individual family sup bounds.
        H=sum((a.config[3] if f=='Yt' else 1)*max(t.values()) for f,t in v.items())
        stage='hV+P4+Yt' if high else 'hV'
        r=report(c,cs,s,stage,ev,H,replace=a.replace)
        assert r[1]+H<Q('.1')
        print('absolute norm/core',float(r[1]+H),flush=True)
        from verify_tree import verify
        verify(E,normal,c,cs,s,ev,v,r)
    if a.robust_upper:
        v={f:layer(c,s,f,dom[f],{},upper=True,replace=a.replace) for f in FAMILIES[1:]} if not a.only_hv else {}
        v['hV']=layer(c,s,'hV',dom['hV'],v,upper=True,replace=a.replace)
        r=report(c,cs,s,'same finite robust class upper',edge(e,s,v['hV'],upper=True),Z,upper=True,replace=a.replace)
        cap=(Q('5.53e-7') if c.sign>0 else Q('2.80e-7')) if a.replace else (Q('2.42e-6') if c.sign>0 else Q('2.25e-6'))
        print('robust class shortfall even granting frozen reserve ceiling',float(-r[0]-cap),flush=True)


if __name__=='__main__':main()
