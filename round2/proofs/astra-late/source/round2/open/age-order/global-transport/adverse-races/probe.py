"""Exploratory (not outward) limiting continuous-time Bellman solver."""
from math import exp, lgamma, log
from pathlib import Path
from runpy import run_path
from fractions import Fraction as Q
from types import FunctionType
from argparse import ArgumentParser

HERE=Path(__file__).resolve().parent
F=run_path(str(HERE.parent/'annealed/extension/third-physical/scaled/check.py'))


def gamma_tail(a,v):
    # Incomplete gamma continued fraction, used only at v>=1.8.
    b=v+1-a
    c=1e100
    d=f=1/b
    for i in range(1,150):
        an=-i*(i-a)
        b+=2
        d=1/(an*d+b)
        c=b+an/c
        fac=c*d
        f*=fac
        if abs(fac-1)<1e-14:
            break
    return exp(a*log(v)-v-lgamma(a))*f


def lim(y,D=32,N=4096):
    T=2
    row=[[0.]*k for k in range(D+2)]
    for i in range(N,-1,-1):
        v=y*exp(i*T/N)
        for k in range(2,D+2):
            locked=gamma_tail(1/k,v)
            if k==D+1 or i==N:
                row[k]=[locked]*k
                continue
            s=row[k]
            sup=gamma_tail(1+1/k,v)/k
            for j in range(k-1,-1,-1):
                g=row[k+1][j+1]-s[j]
                s[j]=min(locked,s[j]+T/N*(sup-s[j]/k+j/k*((s[j-1] if j else 0)-s[j])+max((j+1)*g,(k+j+1)/2*g)))*1
            # Need ascending degree: gains read the next row. Loss also reads next row.
        if i==0:
            return {(k,j):float(exp(y)*row[k][j]) for k in range(2,D+2) for j in range(k)}


def signed_finite(M,U,b=33):
    A=run_path(str(HERE.parent/'annealed/extension/third-physical/check.py'))
    return FunctionType(A['ledger'].__code__,dict(A,tips=lambda _:{p:Q(u) for p,u in U.items()}))(b,Q(M),*map(Q,('.004','.008','.006','.004')))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('y', nargs='*', type=float, default=[2.,2.5,3.,3.5,4.])
    p.add_argument('--steps', type=int,default=4096)
    for y in p.parse_args().y:
        U=lim(y,N=p.parse_args().steps)
        print(y, {key:round(U[key],8) for key in [(2,0),(3,0),(3,1),(4,0),(4,2)]},flush=True)
        for M in ('.36787','.462'):
            print(M,signed_finite(M,U),flush=True)
