#!/usr/bin/env python3
"""Coarse rational all-arm bounds without scanning any motif table."""
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check import LO, HI, FAMILIES, Cert, common_age, upward  # noqa: E402


def polygon(m):
    """Coefficients of E[(j+N)^m] as polynomial in j, N geometric."""
    stirling = [1]
    moments = [Q(1)]
    for k in range(1, m+1):
        stirling = [0] + [i*(stirling[i] if i < k else 0)+stirling[i-1] for i in range(1,k+1)]
        moments.append(sum(v*factorial(i)/LO**i for i,v in enumerate(stirling)))
    return [comb(m,i)*moments[m-i] for i in range(m+1)]


def arity(K, r):
    return sum(x*(K+2)**i for i,x in enumerate(polygon(r))) / (1+LO)**(K+1)


def arm_size(B):
    # Size+1 starts at three; dominate by a rate-n Yule chain and observe
    # at the actual Exp(mu) birth age, monotonically bounded by Exp(LO).
    p = Q(1)
    for n in range(3,B+2):
        p *= n/(n+LO)
    return p


def norms(c):
    # Authenticated outward control maxima (embedding-flow/full-tail).
    g,hv,hc,p4,yt,yb = (Q(a,c.S) for a in
                       (26340094623,10079031623,13860049440,2377223860,4095258789,2871334141))
    center = max(Q(abs(int(a)),c.Sr) for t in (c.rpar,c.rch,c.trpar,c.trch) for a in t)
    gref = max(Q(abs(int(a)),c.S) for t in (c.yp,c.yw) for a in t.flat)
    gmove = max(2*g,6*g/c.D+gref)
    # At maximal arity j all receive rates are at most (j+1)/2.
    # Per tag: leaf gain <=34*w; receive increment <=(j+1)*w;
    # ordinary promotion <=2*w. At most one extension per indicated role
    # (the sum of all possible replacements is bounded by len(ext)*max).
    for family,r,w,extension,route,parents in (
            ('hV',3,hv,3*max(p4,yt,yb),2*hv,2*gmove),
            ('P4',4,p4,0,0,4*hv),
            ('Yt',4,yt,0,0,4*hv+2*hc)):
        a = HI*w+r*37*w+r*extension+(Q(1,2)+center)*parents
        b = r*(w+(Q(1,2)+center)*route)
        yield family,r,a,b


def main(blob):
    c = Cert(Path(blob))
    assert (c.D,c.S)==(18,10**11)
    D = 101
    # On max arity <=K, an occurrence rooted at the offending first arm has
    # one apex and at most K^(r-2) descendants; at the apex K^(r-1).
    P = common_age(LO,D,0)[0]
    for K,B in ((12,9),(30,1000),(60,10_000)):
        Psize = arm_size(B)
        for family,r,a,b in norms(c):
            high_arity = r*((a+b)*arity(K,r))
            large_arm = (a+b*K)*min(Psize*K**(r-2),P*K**(r-1))
            print(c.mode,K,B,family,'omitted arity',upward(high_arity),
                  'omitted arm',upward(large_arm),flush=True)


if __name__ == '__main__':
    for blob in sys.argv[1:]:
        main(blob)
