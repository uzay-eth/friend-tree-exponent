#!/usr/bin/env python3
"""Small exact labelled-clock replay of the integer producer, no C++ row code."""
from fractions import Fraction as F
from functools import cache
from itertools import product
import subprocess
import sys

Q=10**10
s=F('0.565658')


def up(x): return -(-x*Q//1)


@cache
def coupons(d,m):
    step=64
    disc=F(up(1/(1+F(int(s*Q/step),Q))),Q)
    x=F(up(F(d+1,step*d)),Q)
    t=F(1); g=t
    for j in range(1,13):
        t=F(up(t*x/j),Q); g+=t
    g+=2*F(up(t*x/13),Q)
    age=F(1); toll=F(up(F(2*m,d)),Q); v=F(0)
    while toll<2:
        toll=F(up(toll*g),Q)
        v+=age*min(2,toll)/step
        age=F(up(age*disc),Q)
    return F(up(v)+up(2*age/s),Q)


def check(exe,N,plant=None):
    @cache
    def one(d,j):
        if d>N: return F(up(min(1/s,(j+1/s)/(1+s*d))),Q)
        b=F(d-1-j,d)
        v=[]
        for ext in (0,F(j,2)):
            r=d-1-j+ext
            v.append((F(j,d)+r*one(d+1,j)+(b*one(d,j+1) if b else 0))/(s+r+b))
        return F(up(max(v)),Q)

    @cache
    def value(u,v,c,e,L,K):
        Y=(u,v,c,e,L,K)
        B=1+u+v+c+e; A=B+L+K
        if A>N: return min(F(up(2/s),Q),one(B,v+c+e)+one(A,c+e+K),coupons(B,min(B,v+e+L+K)))
        def append(i):
            z=list(Y); z[i]+=1; return tuple(z)
        fixed=[(u,append(0)),(L,append(4))]
        private=['l']*L+['k']*K
        paired=['u']*u+['v']*v+['c']*c+['e']*e+['slot']
        for name in paired:
            z=list(Y)
            if name=='u': z[0]-=1;z[2]+=1
            if name=='v': z[1]-=1;z[3]+=1
            fixed.append((F(1,A),tuple(z)))
        for aa,bb in product(private,paired):
            z=list(Y)
            if aa=='l': z[4]-=1;z[5]+=1
            if bb=='u': z[0]-=1;z[1]+=1
            if bb=='c': z[2]-=1;z[3]+=1
            if plant=='split-root' and aa=='l' and bb in ('u','c'):
                y=list(z);y[4]+=1;y[5]-=1
                w=list(Y);w[4]-=1;w[5]+=1
                fixed.extend([(F(1,A*B),tuple(y)),(F(1,A*B),tuple(w))]);continue
            fixed.append((F(1,A*B),tuple(z)))
        C,P=append(0),append(4)
        alternatives=[([(0,C),(v,P)],[(F(v,2),C),(F(v,2),P)]),([],[(F(c,2),C)]),
                      ([],[(F(e,2),C)],[(F(e,2),P)]),([],[(F(K,2),P)] if plant!='drop-private-core' else [])]
        r=F(c+e+v,B)-F(c-e-K,A)
        res=[]
        for blocks in product(*alternatives):
            edges=[(t,z) for t,z in fixed+sum(blocks,[]) if t and z!=Y]
            res.append((r+sum(t*value(*z) for t,z in edges))/(s+sum(t for t,z in edges)))
        return F(up(max(res)),Q)
    output=subprocess.run([exe,str(N),'all'],check=True,text=True,capture_output=True).stdout
    for b,u,score in (map(int,l.split()) for l in output.splitlines()):
        ref=min(value(u,0,b-u-1,0,1,0),coupons(b,1))
        if score!=up(ref):
            raise ValueError(f'd={b} L={u}: expected {score}, independent {up(ref)}')
    print(f'PASS cap {N}: labelled physical clock replay and exact rational Bellman rows')


if __name__=='__main__':
    try:check(sys.argv[1],int(sys.argv[2]),sys.argv[3] if len(sys.argv)>3 else None)
    except ValueError as e:sys.exit(f'REJECT: {e}')
