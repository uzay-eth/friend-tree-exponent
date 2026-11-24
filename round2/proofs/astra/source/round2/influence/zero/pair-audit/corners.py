#!/usr/bin/env python3
"""Exact independent Bellman corner replay for the six-count relaxation."""

from fractions import Fraction as Q
from functools import cache
from itertools import product
import argparse
import subprocess

from audit import jump

S=Q(282829,500000)
TICKS=10**10


def up(t):
    return (t*TICKS).__ceil__()/Q(TICKS)


@cache
def coupon(d,m):
    if m>=d: return up(2/S)
    x=up((1+Q(1,d))/64)
    t=g=Q(1)
    for i in range(1,13):
        t=up(t*x/i); g+=t
    g+=2*up(t*x/13)
    disc=up(1/(1+(S*TICKS/64).__floor__()/Q(TICKS)))
    age=Q(1); rate=up(Q(2*m,d)); answer=0
    while rate<2:
        rate=up(rate*g)
        answer+=age*min(2,rate)/64
        age=up(age*disc)
    return up(answer)+up(age*2/S)


def rate(s,p,k,f,lc=Q(0)):
    u,v,c,e,l,z=s
    q=Q(1,1+u+v+c+e)
    a=Q(1,1+sum(s))
    R={}
    def pay(w,*d):
        t=jump(s,d)
        if w and t!=s: R[t]=R.get(t,0)+w
    pay(u+lc+c*p+(Q(e,2) if k==2 else 0),1,0,0,0,0,0)
    pay(l+v-lc+(Q(e,2) if k==1 else 0)+f*z,0,0,0,0,1,0)
    pay(a*u,-1,0,1,0,0,0)
    pay(a*v,0,-1,0,1,0,0)
    pay(a*q*u*l,-1,1,0,0,-1,1)
    pay(a*q*c*l,0,0,-1,1,-1,1)
    pay(a*q*(1+v+e)*l,0,0,0,0,-1,1)
    pay(a*q*u*z,-1,1,0,0,0,0)
    pay(a*q*c*z,0,0,-1,1,0,0)
    return R


def table(cap,side_bound=True):
    @cache
    def single(d,j):
        if d>=cap+1:
            if side_bound in ('discounted-arity','production'):
                return min(up(1/S),up((j+1/S)/(1+S*d)))
            return min(up(1/S),up((j/S+1/S**2)/d))
        promote=Q(d-1-j,d); leaves=d-1-j
        return max(up((Q(j,d)+(promote*single(d,j+1) if promote else 0)
                       +(leaves+h)*single(d+1,j))/(S+leaves+h+promote)) for h in (0,Q(j,2)))

    @cache
    def pair(s):
        u,v,c,e,l,k=s
        d=1+sum(s); b=1+u+v+c+e
        if d>=cap+1:
            v1=min(up(2/S),single(d,c+e+k)+single(b,v+c+e)) if side_bound else up(2/S)
            return min(v1,coupon(b,v+e+l+k)) if side_bound=='production' else v1
        reward=c*(Q(1,b)-Q(1,d))+e*(Q(1,b)+Q(1,d))+Q(v,b)+Q(k,d)
        ans=0
        for p,q,f,lc in product((Q(0),Q(1,2)),range(3),(Q(0),Q(1,2)),(Q(0),Q(v,2))):
            R=rate(s,p,q,f,lc)
            ans=max(ans,up((reward+sum(w*pair(t) for t,w in R.items()))/(S+sum(R.values()))))
        return ans
    return pair


def test(cap,exe=None):
    f,g,h=table(cap),table(cap,False),table(cap,'discounted-arity')
    prod=table(cap,'production')
    seeds={}
    for d in range(2,cap):
        for leaves in range(d):
            a=(leaves,0,d-1-leaves,0,1,0)
            seeds[d,leaves]=int(min(prod(a),coupon(d,1))*TICKS)
    if exe:
        proc=subprocess.run([exe,str(cap),'check'],capture_output=True,text=True,check=True)
        theirs={(d,l):w for d,l,w in (map(int,line.split()) for line in proc.stdout.splitlines())}
        assert theirs==seeds,(theirs,seeds)
        print('C++ cross-replay:',len(seeds),'seed bounds exactly equal independent corners')
    for a in [(0,0,0,1,1,0),(1,0,0,1,1,0),(0,1,0,1,1,0)]:
        print(a,'controlled terminal',float(f(a)),'arity terminal',float(h(a)),'constant terminal',float(g(a)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cap',type=int,default=6)
    p.add_argument('--exe',type=str)
    a=p.parse_args()
    test(a.cap,a.exe)
