#!/usr/bin/env python3
"""Authenticate D18 input, stream exact P4 directional charge operands."""
import argparse
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
S=HERE.parents[1]/'replay.py'
assert sha256(S.read_bytes()).hexdigest()=='e90651922acb6b3a427872cf3801f4b9097d2b83c7abe0ab75add6e60c78350f'
spec=spec_from_file_location('frozen_transport_edge',S)
E=module_from_spec(spec);spec.loader.exec_module(E)
Q=E.Q
V=10**18
LO,HI=Q('.565673'),Q('.565775')


def states(d,k):return [(i,j) for i in range(d,1,-1) for j in range(min(i-1,k),-1,-1)]


def clear_frozen_p4(d,k,x):
    """Count the whole fringe from below using the four tagged path vertices."""
    path=((d,k),*x)
    off=sum(j for _,j in path)-3
    size=4+sum(i-j-1 for i,j in path)+off
    cores=4+off
    return size>12 and (size>20 or cores>4 or max(i for i,_ in path)>16)


def coefficient(c,d,k,x,family):
    z=tuple(map(E.desc,((d,k),*x)))
    i=3
    if family=='Yb':
        z=(z[0],z[3],z[1],z[2]);i=1
    if not c.supported(family,z):return Q(0)
    a,free=c.terms(family,z)
    assert a+sum(min(0,n*v/2) for _,_,n,v in free)>=0
    t=tuple(v for j,p,n,v in free if j==i and not p)
    assert len(t)==bool(x[-1][1])
    return t[0] if t else Q(0)


def outward(x):return (x*V//1,-(-x*V//1))


def stream(proc,args,tables,weights):
    da,ka,db,kb,dc,kc,dt,kt,dd,kd=args.config
    paidmax=args.paidmax or da
    print(da,ka,db,kb,dc,kc,dt,kt,dd,kd,paidmax,int(args.family=='Yb'),file=proc.stdin)
    by_k={};peak=Q(0)
    for d in range(paidmax,17 if args.family=='P4' else 1,-1):
        slab=[]
        for k,x in product(range(ka+1),product(states(db,kb),states(dc,kc),states(dt,kt))):
            if k<2 or k>=d or x[0][1]<1 or not x[2][1] or (x[1][1]<1 if args.family=='P4' else x[2][0]<18):a=b=Q(0)
            elif args.family=='P4' and not clear_frozen_p4(d,k,x):a=b=Q(0)
            else:
                if args.family=='P4':
                    assert d>=18 and d+sum(i for i,_ in x)-3>=21 and k+sum(j for _,j in x)+1>=6
                ss=[coefficient(c,d,k,x,args.family) if w else Q(0) for c,w in zip(tables,weights)]
                a=sum(w*max(Q(0),-t)/2 for w,t in zip(weights,ss));b=sum(w*t for w,t in zip(weights,ss))
                if a or b:
                    by_k[k]=by_k.get(k,0)+1
                    peak=max(peak,abs(a),abs(a+b/2))
            slab.append(' '.join(map(str,(*outward(a),*outward(b)))))
        proc.stdin.write('\n'.join(slab)+'\n')
    proc.stdin.close()
    return by_k,peak


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lower',type=Path);p.add_argument('upper',type=Path);p.add_argument('--native',type=Path,required=True)
    p.add_argument('--config',nargs=10,type=int,default=(48,5,5,2,5,2,7,2,8,2));p.add_argument('--paidmax',type=int)
    p.add_argument('--objective',choices=('lower','upper','both','contrast'),default='contrast');p.add_argument('--theta',type=Q,default=Q('1.2415'))
    p.add_argument('--family',choices=('P4','Yb'),default='P4')
    a=p.parse_args();assert not sys.flags.optimize and a.theta>=0
    a.native=a.native.resolve()
    da,ka,*remainder=a.config
    assert (18 if a.family=='P4' else 3)<=(a.paidmax or da)<=da<=2048 and 2<=ka<da and all(2<=d<=64 and 1<=k<d for d,k in zip(remainder[::2],remainder[1::2]))
    for path,h in E.SOURCE_HASHES.items():assert sha256((HERE.parents[4]/path).read_bytes()).hexdigest()==h,path
    c=[]
    for p,mode in zip((a.lower,a.upper),(2,1)):
        h=sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(2**23),b''):h.update(chunk)
        assert h.hexdigest()==E.PAYLOAD_HASHES[mode]
        c.append(E.Charges(E.D.Cert(p)))
        assert (c[-1].c.mode,c[-1].c.D,c[-1].c.S)==(mode,18,10**11)
        print('payload',h.hexdigest(),flush=True)
    for objective in ('lower','upper') if a.objective=='both' else (a.objective,):
        w={'lower':(1,0),'upper':(0,1),'contrast':(1,-a.theta)}[objective]
        with subprocess.Popen([str(a.native.resolve())],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True) as proc:
            rows,peak=stream(proc,a,c,w)
            raw=proc.stdout.readline().split();assert proc.wait()==0 and raw
        l,u=map(int,raw[:2]);bounds=(min(LO*l,HI*l)/V,max(LO*u,HI*u)/V)
        print('family',a.family,'objective',objective,'config',a.config,'nonzero toll rows by apex arity',sorted(rows.items()),
              'max one-donor endpoint toll',float(peak),flush=True)
        print('native',*raw,'physical normalized interval',*[float(q) for q in bounds],flush=True)
        print('exact normalized interval',*bounds,flush=True)
        for t in c:
            t.terms.cache_clear();t.weight.cache_clear()


if __name__=='__main__':main()
