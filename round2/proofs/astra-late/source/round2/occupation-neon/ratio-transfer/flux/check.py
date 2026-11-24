#!/usr/bin/env python3
"""Exact incident-normalized observation countermodel for the D18 S3 hinge."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'fringe'))
from charges import Charges,Cert

THETA=Q('1.2415')
HOST=(2,3,1)
TIPS=((1,0,1),(2,0,1))


def fixture(bits):
    # Actual rooted tree with one degree-three founder, degree-six host,
    # and three named sterile tips. The row uses true reciprocal degrees;
    # its only exterior donor is the founder.
    p=(-1,0,1,1,1)
    x=((2,1,0),HOST,*(TIPS[b] for b in bits))
    return p,x


def measure(bits):
    p,x=fixture(bits)
    one=defaultdict(Q);sterile=defaultdict(Q);pair=defaultdict(Q);incident=defaultdict(Q);star2=defaultdict(Q)
    for state in x:one[state]+=Q(1,len(x))
    for state in x[2:]:sterile[state]+=Q(1,3)
    for j,i in enumerate(p[1:],1):
        pair[x[i],x[j]]+=Q(1,len(x))
        for a,b,down in ((i,j,True),(j,i,False)):
            incident[x[a],x[b],down]+=Q(1,len(x)*sum(x[a]))
    for i,j in product(range(2,5),repeat=2):
        if i!=j:star2[x[i],x[j]]+=Q(1,6)
    return (one,sterile,pair,incident,star2)


def charge(E,bits):
    p,x=fixture(bits)
    out=[]
    for e in E:
        a,free=e.terms('S3',x[1:])
        assert all(i==0 and par and n==1 for i,par,n,_ in free)
        a+=sum(c/Q(sum(x[0])) for _,_,_,c in free)
        assert a>=0
        out.append(a)
    return out[0]-THETA*out[1]


def combine(items):
    out=[defaultdict(Q) for _ in range(5)]
    for bits in items:
        for dest,src in zip(out,measure(bits),strict=True):
            for k,v in src.items():dest[k]+=v/len(items)
    norm=sum(out[3].values())
    assert norm==sum(m*Q(s[1]+s[2],sum(s)) for s,m in out[0].items())
    assert sum(out[0].values())==sum(out[1].values())==sum(out[4].values())==1
    out[3]={key:v/norm for key,v in out[3].items()}
    assert sum(out[3].values())==1
    return out,norm


def main():
    E=[Charges(Cert(Path(s))) for s in sys.argv[1:]]
    assert len(E)==2 and [e.c.mode for e in E]==[2,1]
    F=[charge(E,(1,)*j+(0,)*(3-j)) for j in range(4)]
    even,odd=([b for b in product((0,1),repeat=3) if sum(b)%2==r] for r in (0,1))
    obs,norm=combine(even)
    assert (obs,norm)==combine(odd)
    print('identical uniform, sterile, oriented edge, normalized incident, and two-arm laws; raw incident mass',norm)
    for j,w in enumerate(F):print('exact S3 contrast at',j,'degree-three arms',w,float(w))
    means=[]
    for name,p in [('even',even),('odd',odd)]:
        vals=[F[sum(bits)] for bits in p]
        m=sum(vals)/len(vals)
        h=sum(max(Q(0),v) for v in vals)/len(vals)
        second=sum(v*v for v in vals)/len(vals)
        means.append(m)
        # Even the exact optimal scalar quadratic majorant using its true
        # moments gives (m+sqrt(second))/2. Check the rational form instead.
        assert 2*h-m>=0 and (2*h-m)**2<=second
        print(name,'signed mean',float(m),'positive hinge mean',float(h),'raw second moment',float(second))
    assert means[0]<Q('-1.91e-6')<0<Q('7.57e-6')<means[1]
    # The fixture has five cores and a single eligible S3 apex. In a
    # shared background with mixture mass eta, multiply this gap by eta.
    print('per-core S3 mean ambiguity at unit fixture mass',float((means[1]-means[0])/5))
    print('opposite signs with identical retained observations: PASS')
    # Separate temporal countermodel with identical complete one-step joint
    # X laws and both stationary marginals uniform. P is either a shift on
    # 0011 or a complete resample. Right-limit variance of the time integral
    # cannot follow from a one-step pairing, even if it matches exactly.
    ring=(Q(-1),Q(-1),Q(1),Q(1))
    cyclic=defaultdict(Q)
    for i,x in enumerate(ring):cyclic[x,ring[(i+1)%4]]+=Q(1,4)
    assert cyclic=={(x,y):Q(1,4) for x,y in product((-1,1),repeat=2)}
    for horizon in (4,8,16):
        # Frozen hidden sign is used only at alternating phases. Consecutive
        # observed pairs remain independent fair signs, but the persistent
        # component produces an O(horizon^2) accumulated variance.
        n=horizon//2
        independent=Q(horizon)
        periodic=sum(sum(ring[(i+j)%4] for j in range(horizon))**2 for i in range(4))/4
        persistent=Q(n*n+n)
        assert periodic==0 and persistent>independent
        print('same one-step law, block',horizon,'variances periodic, iid, persistent',periodic,independent,persistent)


if __name__=='__main__':main()
