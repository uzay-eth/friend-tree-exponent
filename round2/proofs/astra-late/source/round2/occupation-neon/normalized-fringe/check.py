#!/usr/bin/env python3
"""Exact physical counterexample to an eight-digit bulk-normalized separator."""
import argparse
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ratio-transfer'))
from replay import Pair, Cert, tree

A,B=Q('.565673'),Q('.565775')
EPS=Q('1e-8')


def chain(stock,n):
    z=tree(stock)
    for _ in range(n-1):z=tree(stock,z)
    return z


def bulk(pair,stock):
    vals=[]
    for n in (8,9,10):
        d,p=pair.full(chain(stock,n),True)
        vals.append(d+p)
    r=tuple(y-x for x,y in zip(vals[0],vals[1]))
    assert r==tuple(y-x for x,y in zip(vals[1],vals[2]))
    assert min(r)>0
    return r


def sign(row,s,t):
    d,e,p,q=row
    # Match bulk on these same physical inputs; potential and opposite
    # slack are combined before division by a positive denominator.
    base=p+t*q
    return (d-t*e-(s-A)*p+t*(B-s)*q)/base


def root(x,t):
    return A+t*(B-A)*x[3]/(x[2]+t*x[3])+(x[0]-t*x[1])/(x[2]+t*x[3])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('payloads',nargs=2,type=Path)
    a=p.parse_args()
    c=Pair(*(Cert(f) for f in a.payloads))
    rows=[bulk(c,L) for L in (14001,10**6)]
    for L,(d,e,p,q) in zip((14001,10**6),rows):
        lr,ur=A+d/p,B-e/q
        print('leaves/core',L,'forcing/Phi-',float(d/p),'forcing/Phi+',float(e/q),
              'induced roots',float(lr),float(ur),flush=True)
    # For t>=0 each normalized residual is a convex combination of the
    # t=0 and t=infinity normalized residuals. Strict interval separation
    # proves the obstruction for every such t, including minimax pairing.
    first,second=rows
    lo=max(A+first[0]/first[2],B-first[1]/first[3])
    hi=min(A+second[0]/second[2],B-second[1]/second[3])
    assert hi-lo > 100*EPS
    print('all positive multipliers: separated induced roots by at least',float(hi-lo))
    for t in (Q(0),Q(9,10),Q(1),Q(6,5),Q(5)):
        v,w=(root(x,t) for x in rows)
        assert sign(first,(v+w)/2,t)<-EPS and sign(second,(v+w)/2,t)>EPS
        # A common constant centering cannot eliminate this separation.
        lower=abs(w-v)/2
        print('theta',t,'roots',float(v),float(w),'uniform centered radius >=',float(lower))


if __name__=='__main__':main()
