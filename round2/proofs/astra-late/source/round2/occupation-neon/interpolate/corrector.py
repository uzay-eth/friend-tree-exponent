#!/usr/bin/env python3
"""Compile sparse census Bellman correction into a pointwise drift witness."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from importlib.util import spec_from_file_location, module_from_spec
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'compose'))
from sources import authenticated, bands
spec=spec_from_file_location('denominator',HERE.parent/'denom'/'check.py')
D=module_from_spec(spec);spec.loader.exec_module(D)
_BANDS=bands()
FLOOR={1:_BANDS['upper'],2:_BANDS['lower']}


def operands(c,n):
    name='upper' if c.mode==1 else 'lower'
    manifest=json.loads(authenticated('rows/charged-f.json'))
    assert manifest['payload_sha256']=={
        'lower':'9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
        'upper':'1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41'}
    dump=manifest[name]
    assert Q(dump['lambda'])==c.lam and dump['S']==c.S
    rows={}
    for z in dump['rows']:
        l,k,h=z['state']; d=l+k+h
        if d>n:continue
        # Reconstruct exact affine charged row and recheck both pinned corners.
        f=z['f']; gain=z['f_gain']-f
        a=l*gain-c.lam*f
        if l:
            g=sum(i*Q(w) for i,w in z['g_bud_newborn'])
            a+=Q(l,d)*(z['f_bud']+z['f_newborn']+g-f)
        p,r=map(Q,z['center']); center=p+k*r
        if d<c.D:a+=center*gain;b=Q(0)
        else:
            psi=Q(z['c2']*d+z['c3'],d**3);a+=center*psi;b=gain-psi
        sign=-1 if c.mode==1 else 1
        a*=Q(sign,c.S);b*=Q(sign,c.S)
        assert list(map(Q,z['signed_corners']))==([a] if d<c.D else [a,a+b*Q(k+1,2)])
        assert min(a,a+b*Q(k+1,2))>=0
        rows[d,k]=a,b
    assert len(rows)==sum(range(2,n+1))
    return rows


def solve(c,charge,rows,n,s,flat=True):
    shift=abs(s-c.lam);t=Q(-9 if c.mode==1 else -2,1000)*c.S
    M=charge.cap(t)/c.S
    # Fix the correction's birth value to zero so no birth term remains.
    floors=tuple(map(Q,FLOOR[c.mode]));tail={}
    w=Q(0)
    for k in range(c.KF+1,-1,-1):
        deficit=M-(charge.by_k[1,k]+t*(k-1))/c.S
        assert deficit>=0
        a=((floors[k] if k<len(floors) else 0)+shift*deficit) if flat else 0
        corners=(max(0,1-Q(k+1,n+1)),Q(1))
        tail[k]=D.down(a/s if k==c.KF+1 else min((a+b*w)/(s+b) for b in corners))
        assert tail[k]>=0
        if k<c.KF+1:assert all((s+b)*tail[k]<=a+b*w for b in corners)
        else:assert s*tail[k]<=a
        w=tail[k]
    v=defaultdict(Q)
    for d in range(n,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1;b=Q(l,d);a,r=rows[d,k]
            deficit=M-(charge.exact(l,k,1)+t*(k-1))/c.S
            assert deficit>=0
            a+=shift*deficit
            u=v[d+1,k] if d<n else tail.get(k,Q(0));p=v[d,k+1]
            v[d,k]=D.down(min((a+rho*r+(l+rho)*u+b*p)/(s+l+rho+b)
                            for rho in (0,Q(k+1,2))))
            assert all((s+l+rho+b)*v[d,k] <= a+rho*r+(l+rho)*u+b*p
                       for rho in (0,Q(k+1,2)))
    L=s*v[2,0]-shift*M
    norm=max([*v.values(),*tail.values()])+v[2,0]+shift*M/s
    assert norm<Q(1,100)
    # Choose V=z-z(J), so V(J)=0 exactly and every nonfounder pays L.
    # The founder contributes no correction; its cap-only loss shift*M
    # is paid by adding the signed global constant -shift*M/s.
    return L,v[2,0],M,norm,tail


def all8_K(c):
    """Conservative |(generator-s)Phi|/C, including unbounded stock/arity."""
    # Every existing motif count <= c_n C. Every promotion creates at most
    # c_n copies containing its newborn; the same ownership proof applies.
    # The newborn is a tip. Edge, paths, pairs, Y-tail and triple creation
    # per promoter are bounded by their count multipliers. In a bend it is
    # either the bottom tip (<=5 sibling choices) or the side tip (<=5
    # branches times <=6 children each). Allow sixteen times the count
    # multiplier uniformly. Roles with unbounded degree are end slots.
    from math import comb
    mult=(1,1,1,comb(c.K,2),1,comb(7,2),5,comb(6,3))
    tables=(c.f,c.g,c.hvtab,c.hctab,c.p4tab,c.yttab,c.ybtab,c.s3tab)
    arity=(1,2,3,3,4,4,4,4)
    total=Q(0)
    for i,(a,n,E) in enumerate(zip(mult,arity,tables)):
        m=Q(max(abs(int(E.min())),abs(int(E.max()))),c.S)
        # gain <= D, and at finite stock/census transitions D<=LF+KF+2
        # covers every f change and every switch of support/inner controls.
        # Pure Bernstein D>=18: |d/dd basis weights|_1 <= 6*18/D².
        g=2*(c.LF+c.KF+2) if i==0 else 2*(c.D+1)
        total+=a*m*(1+n*(g+2)+16)
    return total


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='+',type=Path)
    p.add_argument('--degree',type=int,default=64)
    p.add_argument('--s',type=Q)
    p.add_argument('--mix',action='store_true')
    p.add_argument('--mode',type=int,choices=(1,2))
    a=p.parse_args();assert 5<=a.degree<=64 and sys.flags.optimize==0
    certs=[D.Cert(path) for path in a.paths]
    for c in certs:
        if a.mode and c.mode!=a.mode:continue
        charge=D.Charge(c);rows=operands(c,a.degree)
        outward={1:1,2:-1}[c.mode]
        if a.s is not None:
            targets=[a.s]
        else:
            # Coarse decimal grid, then bisect the last successful shift.
            low=0;hi=10000000  # ticks at 1e-12, huge relative to the sparse correction
            while hi-low>1:
                mid=(low+hi)//2
                s=c.lam-outward*Q(mid,10**12)
                if solve(c,charge,rows,a.degree,s)[0]>=0:low=mid
                else:hi=mid
            targets=[c.lam-outward*Q(low,10**12)]
        for s in targets:
            if a.mix and a.s is None:s+=outward*Q('1e-10')
            assert 0<=outward*(c.lam-s) and 0<s<1
            floor,v0,m,norm,tail=solve(c,charge,rows,a.degree,s)
            killed=solve(c,charge,rows,a.degree,s,False)[0]
            print('mode',c.mode,'s',D.show(s),'exact per-core drift floor',floor,
                  'killed-tail floor',killed,'z(J)',D.show(v0,16),
                  'cap',D.show(m),'correction abs/core <=',D.show(norm,16),flush=True)
            if floor>=0:
                print('CERTIFIED physical pointwise sign', 'upper' if c.mode==1 else 'lower',flush=True)
                if a.mix:
                    other,= [x for x in certs if x.mode!=c.mode]
                    K=all8_K(other);eps=floor/(2*K)
                    assert eps>0 and norm+eps*K/s<Q('.01')
                    # Founder payment leaves floor*(C-1). Pay the extra
                    # worst-case founder residual with a signed constant.
                    print('all8 scale K upper',D.show(K,6),'strict opposite admixture',eps,
                          'signed constant added',-eps*K/s,flush=True)


if __name__=='__main__':main()
