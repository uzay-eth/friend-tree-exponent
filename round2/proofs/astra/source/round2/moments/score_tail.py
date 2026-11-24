"""Affine upper bounds for unsigned tail-cell scores at the actual birth intensity."""

from array import array
from bisect import bisect_left
from fractions import Fraction as F
from math import ceil
import sys


def build_lines(path, header, end, degree_bands, child_bands, cells, tail):
    side,_,_,unit,df,p0,p1,p2,kt=header
    if sys.byteorder!='little' or p1<=0 or p2<0:
        raise ValueError('unsupported nominal potential')
    slopes={}
    for j in range(16):
        # Include the tangent to the analytic singleton at the bin centre.
        derivative=p1 if side==1 else p0-p1*(2*j+1)//16
        slopes[j]=tuple(sorted({unit//4,unit//2,unit,2*unit,3*unit,max(0,derivative)}-{0}))
    weights=array('q')
    with path.open('rb') as st:
        weights.fromfile(st,df*(df-1)//2-1)
        if st.read(1): raise ValueError('unexpected comparison profile length')
    rows={}
    for d in range(end+1,df):
        db=bisect_left(degree_bands,d)
        start=d*(d-1)//2-1
        for k in range(d):
            leaves=d-k-1
            j=16*leaves//d
            key=db,bisect_left(child_bands,k),j
            price=weights[start+k]
            hh=rows.get(key)
            if hh is None:
                hh=rows[key]=[price-s*leaves//d for s in slopes[j]]
            else:
                for i,s in enumerate(slopes[j]):
                    # ceil(price - slope * leaves/d), entirely in integer ticks.
                    a=price-s*leaves//d
                    if a>hh[i]: hh[i]=a
    for key,(dl,dh,kl,kh,*_) in cells.items():
        if dh: continue
        j=key[-1]
        q0,q1=F(j,16),F(j+1,16)
        if kh is not None: q0=max(q0,1-F(kh+1,dl))
        bottom=max(intercept+max(0,F(slope,df)) for k,intercept,slope in tail
                   if min(kl,kt+1)<=k<=min(kh if kh is not None else kt+1,kt+1))
        entries=[]
        for s in slopes[j]:
            if side==1:
                val=p0+max((p1-s)*q0,(p1-s)*q1)
            else:
                q=min(q1,max(q0,F(p0-s,2*p1)))
                val=(p0-s)*q-p1*q*q
            entries.append(ceil(val+bottom)+4)
        rows[key]=entries
    if rows.keys()!=cells.keys(): raise ValueError('incomplete unsigned score cells')
    return {key:tuple(zip(slopes[key[-1]],v)) for key,v in rows.items()}
