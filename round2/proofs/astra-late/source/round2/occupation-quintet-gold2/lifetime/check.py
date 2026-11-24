#!/usr/bin/env python3
"""Exact Gold2-only original singleton lifetime and all-stock control."""
import argparse
import ast
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from json import loads
from math import ceil, floor
from math import lcm
from itertools import product
from mmap import mmap, ACCESS_READ
from pathlib import Path
from struct import unpack_from
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
R2=HERE.parents[1]


def load(name,path):
    spec=spec_from_file_location(name,path); m=module_from_spec(spec);spec.loader.exec_module(m);return m


tail=load('gold2_tail',R2/'attacks/neon2-quintet/gold2/tail/check.py')
Poly=load('polynomial',R2/'occupation-quintet/census/census.py').Poly
GRID=10**12


def neighbour_algorithm(lo,hi,gold):
    # Reuse only the payload-agnostic function definitions; do not execute
    # the Gold1 module's top-level setup or open its receipt.
    path=R2/'occupation-quintet/lifetime/conditioned/replay.py'
    funcs={'rnd','coefficients','extremum','conditional','birth_mix'}
    ns=dict(Q=Q,lcm=lcm,product=product,G=GRID,base=SimpleNamespace(LO=lo,HI=hi,GOLD=gold))
    tree=ast.parse(path.read_bytes())
    exec(compile(ast.Module([n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in funcs],[]),str(path),'exec'),ns)
    return SimpleNamespace(**{f:ns[f] for f in funcs})


def row(t,D,c,cap):
    L=D-c
    def w(d,k): return t.f(d,k) if d<t.df else t.w(Q(d-k,d),Q(1,d))
    old=w(D,c); jump=w(D+1,c)-old
    if L:
        if D<t.d:
            lab=(D,c+1) if c<t.k else (D,'over')
            nb=t.f(2,1)+t.g(lab,(2,1))
        else:
            nb=t.const+Q(t.d,D)*(t.tj[min(c+1,t.k+1)]-t.gi)
        bud=Q(L,D)*(w(D,c+1)-old+nb)
    else: bud=0
    if D<t.d and c<=t.k:
        idx=t.ids[D,c]; center=Q(unpack_from('<q',t.raw,512+8*idx)[0],2**30)
        return cap-t.s*((L+c*center)*jump+bud-t.lam*old), Q(0)
    return cap-t.s*(L*jump+bud-t.lam*old), -t.s*jump


def suffix(t,D,C,disc,cap,up):
    X,one=Poly([0,1]),Poly([1]); plus=one+X; xmax=Q(1,D)
    A,B,M,E=t.A,t.B,t.C,t.E
    h=[Q(0)]*(C+2); const=h.copy()
    for c in reversed(range(1,C+1)):
        q=one-c*X; w=A+B*q+M*q*q+E*q*X
        gain=B*(one-q)*plus+M*(one-q)*(2*q+(one+q)*X)+E*((one-2*q)*X-q*X*X)
        bud=-B*X+M*(X*X-2*q*X)-E*X*X
        nb=t.const+t.d*(t.tj[min(c+1,t.k+1)]-t.gi)*X
        vals=[]
        for ep,r in enumerate((q,q+Q(c,2)*X)):
            P=(cap-t.s*(q*(bud+nb)-t.lam*w))*plus*plus-t.s*r*gain
            assert P[0]>=0
            for d in (D,D+1,10**8):
                a,b=row(t,d,c,cap)
                assert P(Q(1,d))==(a+ep*Q(c,2)*b)*(1+Q(1,d))**2
            if up: vals.append((P[0],max(0,*Poly(P[1:]).bernstein(xmax))))
            else: vals.append(max(0,min((Poly(P[1:])+D*P[0]).bernstein(xmax)))/(1+xmax)**2)
        rates=(1-c*xmax,1)
        if up:
            k0=max(y[0] for y in vals); k1=max(y[1] for y in vals); decay=Q(D-c,D+1)
            const[c]=Q(ceil(GRID*max((k0+b*const[c+1])/(disc+b) for b in rates)),GRID)
            h[c]=Q(ceil(GRID*max((k1+b*h[c+1])/(disc+decay+b) for b in rates)),GRID)
            assert all((disc+b)*const[c]>=k0+b*const[c+1] and (disc+decay+b)*h[c]>=k1+b*h[c+1] for b in rates)
        else:
            k=min(vals)
            h[c]=Q(floor(GRID*min((k+b*h[c+1])/(disc+1+b) for b in rates)),GRID)
            assert all((disc+1+b)*h[c]<=k+b*h[c+1] for b in rates)
    rnd=ceil if up else floor
    return [rnd(GRID*(a+b/D)) for a,b in zip(const,h)]


def price(t,H,C,disc,cap,up,saved=None):
    nxt=suffix(t,H+1,C,disc,cap,up) if H+1>=t.df else [0]*(C+2)
    assert H+1>=t.df or not up
    q,p=disc.denominator,disc.numerator; rnd=ceil if up else floor; n=0
    for D in reversed(range(2,H+1)):
        cur=[0]*(C+2); w=2*D*q
        for c in reversed(range(1,min(D,C)+1)):
            L=D-c; a,b=row(t,D,c,cap)
            assert min(a,a+Q(c,2)*b)>=0
            num=(rnd(a*GRID)+L*nxt[c])*w+2*L*q*cur[c+1]
            den=2*D*p+2*L*q+L*w
            dx=c*D*q; num2=num+dx*(rnd(b*GRID)+nxt[c]); den2=den+dx
            v=max(-(-num//den),-(-num2//den2)) if up else min(num//den,num2//den2)
            assert v>=0 and (v*den>=num and v*den2>=num2 if up else v*den<=num and v*den2<=num2)
            if D in (2,t.d-1,t.d,t.df-1,t.df,H) and c<=4:
                z=[((a+b*r)*GRID+(L+r)*nxt[c]+Q(L,D)*cur[c+1])/(disc+L+r+Q(L,D)) for r in (0,Q(c,2))]
                assert (max(z)<=v if up else min(z)>=v)
            cur[c]=v;n+=1
        if saved is not None and D<=saved[0]: saved[1][D]=cur
        nxt=cur
    return Q(nxt[1],GRID),n


def fmt(x,up,places=15):
    scale=10**places; k=(ceil if up else floor)(x*scale)
    return f'{k//scale}.{k%scale:0{places}}'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--degree',type=int,default=7999);p.add_argument('--arity',type=int,default=30)
    p.add_argument('--side',choices=('lower','upper'),action='append');p.add_argument('--upper-control',action='store_true')
    p.add_argument('--conditioned',type=int,metavar='D',help='condition on real mother and latest child through degree D')
    a=p.parse_args(); assert a.degree>=2 and 1<=a.arity<=min(a.degree,100)
    assert (tail.BASE/'hashes.log').read_text().splitlines()==['quintet_lower.bin: OK','quintet_upper.bin: OK']
    gold=loads(tail.RECEIPT.read_text()); lo,hi=[Q(gold[x]['certified_exact']) for x in ('lower','upper')]
    assert (gold['lower']['cap_f'],gold['upper']['cap_f'])==(55,2359)
    for side in a.side or ('lower','upper'):
        rec=gold[side]; up=a.upper_control
        with (tail.BASE/f'quintet_{side}.bin').open('rb') as fp,mmap(fp.fileno(),0,access=ACCESS_READ) as raw:
            assert len(raw)==tail.SIZE
            t=tail.Snapshot(raw,rec); disc=lo if up else hi
            saved=(a.conditioned+1,{}) if a.conditioned else None
            v,n=price(t,a.degree,a.arity,disc,rec['cap_f'],up,saved)
            if a.conditioned:
                assert 2<=a.conditioned<a.degree and a.arity>=18
                m=neighbour_algorithm(lo,hi,gold)
                def ticks(d,c):
                    x,y=row(t,d,c,rec['cap_f']); return x-rec['cap_f'],y
                roots,states,checks=m.conditional(SimpleNamespace(side=side,row_ticks=ticks),saved[1],a.conditioned,18,2,3,up)
                z=m.birth_mix(roots,up)
                assert z<=v if up else z>=v
                print(side,'conditioned',a.conditioned,'states/corners',states,checks,'birth roots',roots,'mixed',z,flush=True)
                v=z
            reward=(hi if up else lo)*v
            if up: reward+=130000*2**40/(1+lo)**a.arity
            net=reward-rec['err_per_core']; den=rec['m_max'] if net>=0 else rec['m_min']
            endpoint=t.lam-t.s*net/den
            inward=t.s*(Q(rec['certified_exact'])-endpoint)
            print(side,'upper control' if up else 'lower price','n',n,'V ticks',v,'reward ticks',fmt(reward,up,9),
                  'den',den,'endpoint',fmt(endpoint,(t.s>0)!=up),'inward',fmt(inward,up),flush=True)


if __name__=='__main__': main()
