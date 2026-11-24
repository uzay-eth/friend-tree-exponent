#!/usr/bin/env python3
"""Exact signed child-occupation rows for the killed zero-mark renewal."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'current-mark'/'eval'))
from gauge import J, deg, plus, size, intrinsic, shapes


def occupations(N,s):
    sh=shapes(N)
    rows={}
    for S in sh:
        p=Q(1,deg(S)); div=s+size(S)-p
        burden=defaultdict(Q)
        for Y in S[1]: burden[Y]+=p
        b=S[0]*p
        for Z,r in intrinsic(S).items():
            if Z not in rows: continue
            f,cs=rows[Z]; b+=r*f
            for Y,t in cs.items(): burden[Y]+=r*t
        rows[S]=(b/div,{Y:t/div for Y,t in burden.items()})
    domain=sorted({Y for _,row in rows.values() for Y in row})
    atom={Z:i for i,Z in enumerate(domain)}
    out=[]; forcing=[]
    for S in domain:
        b,A=rows[plus(S)]
        a,B=rows[S]
        out.append({atom[Y]:A.get(Y,0)-B.get(Y,0) for Y in A.keys()|B.keys()
                    if A.get(Y,0)!=B.get(Y,0)})
        forcing.append(b-a)
    return domain, rows, out, forcing


def enclosure(rows,f,weights):
    n=len(rows)
    y=[0.]*n
    ar=[[(j,float(t)) for j,t in r.items()] for r in rows]
    bf=list(map(float,f))
    for _ in range(80):
        y=[a+sum(t*y[j] for j,t in r) for a,r in zip(bf,ar)]
    scale=10**18
    x=[Q(round(t*scale),scale) for t in y]
    if all(w==1 for w in weights):
        drift=[1-sum(map(abs,r.values())) for r in rows]
    else:
        drift=[weights[i]-sum(abs(t)*weights[j] for j,t in r.items()) for i,r in enumerate(rows)]
    if not drift or min(drift)<=0:
        raise ValueError('weighted absolute norm not contractive')
    err=max(abs(x[i]-f[i]-sum(t*x[j] for j,t in r.items()))/drift[i] for i,r in enumerate(rows))
    return x,err


def upper(x):
    k=-(-x*10**12//1)
    return f'{k//10**12}.{k%10**12:012d}'


def lower(x):
    k=x*10**12//1
    return f'{k//10**12}.{k%10**12:012d}'


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--cutoff',type=int,default=8)
    a.add_argument('--s',default='0.565728')
    a.add_argument('--unweighted',action='store_true')
    a.add_argument('--split',action='store_true',help='separate h at a virgin child; exact Sherman-Morrison')
    args=a.parse_args()
    s=Q(args.s)
    if args.cutoff<3 or s<=0: a.error('need cutoff >= 3 and s > 0')
    domain,V,T,f=occupations(args.cutoff,s)
    w=[Q(1)]*len(T)
    norm=max(sum(map(abs,row.values())) for row in T)
    print('N',args.cutoff,'intrinsic shapes',len(V),'difference rows',len(T),'TV max',upper(norm))
    # For a shaping proposal use (I-|T|)^-1 1 if the iteration is positive.
    for _ in range(0 if args.unweighted else 12):
        w=[1+sum(abs(t)*w[j] for j,t in r.items()) for r in T]
        w=[Q(-(-v*10**12//1),10**12) for v in w]
    c=max(sum(abs(t)*w[j] for j,t in r.items())/w[i] for i,r in enumerate(T))
    print('shaped TV',upper(c),'worst index',max(range(len(T)),key=lambda i:sum(map(abs,T[i].values()))))
    width=max(max(sum(t for t in r.values() if t>0),-sum(t for t in r.values() if t<0)) for r in T)
    print('max oscillation coefficient on 0<=h<=1',upper(width))
    x,e=enclosure(T,f,w)
    reward,row=V[J]; idx={S:i for i,S in enumerate(domain)}
    est=reward+sum(t*x[idx[S]] for S,t in row.items())
    rad=e*sum(t*w[idx[S]] for S,t in row.items())
    print('killed U0(J)',lower(est-rad),upper(est+rad),'H residual correction',upper(e))
    if args.split:
        pivot=idx[J]
        sums=[sum(r.values()) for r in T]
        A=[dict(row) for row in T]
        for r,t in zip(A,sums): r[pivot]=r.get(pivot,0)-t
        c0=max(sum(map(abs,r.values())) for r in A)
        xa,ea=enclosure(A,f,[Q(1)]*len(T))
        z,ez=enclosure(A,sums,[Q(1)]*len(T))
        denlo=1-z[pivot]-ez; denhi=1-z[pivot]+ez
        if denlo<=0: raise ValueError('uncertified rank-one denominator')
        ends=[yy/dd for yy in (xa[pivot]-ea,xa[pivot]+ea) for dd in (denlo,denhi)]
        a0,b0=min(ends),max(ends)
        low=high=reward
        for S,t in row.items():
            i=idx[S]
            v=[xx*yy for xx in (z[i]-ez,z[i]+ez) for yy in (a0,b0)]
            low+=t*(xa[i]-ea+min(v)); high+=t*(xa[i]+ea+max(v))
        if high<est-rad or low>est+rad: raise ValueError('rank-one action disjoint')
        print('rank-one remainder TV',upper(c0),'denominator',lower(denlo),upper(denhi),
              'killed U0(J)',lower(low),upper(high))


if __name__=='__main__':
    main()
