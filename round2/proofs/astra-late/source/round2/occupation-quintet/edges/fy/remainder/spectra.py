#!/usr/bin/env python3
"""Pointwise FY remainder envelopes on original physical endpoint domains."""
from fractions import Fraction as Q
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from replay import Edges, Snapshot, PRIVATE, LO, HI


class Spectrum:
    def __init__(self,e):
        self.e=e
        s=e.snap
        self.high=tuple(i for (d,c),i in s.inner.items() if d>=6)
        self.low=tuple(i for (d,c),i in s.inner.items() if d<=5)
        self.all_high=self.high+tuple(range(len(s.inner),s.inf+1))
        self.high19=tuple(i for (d,c),i in s.inner.items() if d>=19)+tuple(range(len(s.inner),s.inf+1))

    def slopes(self,h):
        d,c=h;s=self.e.snap
        a=s.inner[h];b=s.inner.get((d+1,c),s.tail_base+c-1)
        return tuple(s.g(b,i)-s.g(a,i) for i in range(s.inf+1))

    def width(self,h,kind,zz):
        d,c=h;e=self.e
        if kind=='hh':
            ids=self.all_high if d>18 else self.high19
            left,right=Q(0),Q(1,6)
        elif kind=='all':
            ids=self.low+self.all_high if d>18 else self.high19
            left,right=Q(0),Q(1,2)
        else:raise ValueError(kind)
        lo,hi=min(zz[i] for i in ids),max(zz[i] for i in ids)
        ts=(left,right,*(t for t in e.knots if left<t<right))
        aks=tuple(e.fy(h,t) for t in e.knots)
        ats=tuple(e.fy(h,t) for t in ts)
        def payoff(z):
            return max(e.sig*(t*z-a) for t,a in zip(e.knots,aks))-min(e.sig*(t*z-a) for t,a in zip(ts,ats))
        return max(payoff(lo),payoff(hi))

    def upper(self,r,max_stock=26):
        """Physical singleton Bellman ceiling for host-only bounded reward."""
        # High exterior rate is never assumed jointly attainable across hosts.
        # Monotone degree and stock mean stopped boundaries have zero value.
        v={};grid=10**5
        for d in range(179,1,-1):
            for c in range(min(d,max_stock),0,-1):
                b=Q(d-c,d);g=v.get((d+1,c),0);p=v.get((d,c+1),0)
                a=r.get((d,c),0)*grid+b*p
                vs=[(a+t*g)/(LO+b+t) for t in (d-c,d-c+Q(c,2))]
                m=max(vs);u=-(-m.numerator//m.denominator)
                assert all(u>=x for x in vs)
                v[d,c]=u
        return HI*Q(v[2,1],grid*self.e.snap.scale)


def main():
    for side in ('lower','upper'):
        with Snapshot(PRIVATE/f'quintet_{side}.bin',side,verify=False) as s:
            e=Edges(s);p=Spectrum(e)
            rs=[{},{}]
            for d,c in s.inner:
                if c<2 or (d,c) in e.mids:continue
                q=c*(c-1)
                zz=p.slopes((d,c))
                rs[0][d,c]=q*p.width((d,c),'hh',zz)
                if c>=9:rs[1][d,c]=q*p.width((d,c),'all',zz)
            for tag,r in zip(('both-high all arities','all pairs arity >=9'),rs):
                print(side,tag,'rows',len(r),'worst',sorted(((float(x/s.scale),k) for k,x in r.items()),reverse=True)[:4],flush=True)
                for lo,hi in ((2,7),(8,12),(9,12),(10,26),(11,26),(13,15),(16,18),(19,26),(13,26),(2,26)):
                    v=p.upper({k:x for k,x in r.items() if lo<=k[1]<=hi})
                    print('arity',lo,hi,'ceiling',v,float(v),flush=True)
            s.check_unchanged()
    print('PASS')


if __name__=='__main__':main()
