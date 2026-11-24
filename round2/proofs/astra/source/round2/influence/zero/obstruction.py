#!/usr/bin/env python3
"""Rational tests of the proposed constant split and a size-cut artifact."""
from fractions import Fraction as Q


def reciprocal(d,lambda_):
    # A degree-d source has gain rate <= its degree. 1/(D-1) has
    # drift >= minus itself; Jensen gives the required 1/D lower bound.
    return 1/(d*(lambda_+1))


def artificial(N,s):
    # Root bare, with one virgin marked-by-shape branch and arbitrary
    # padding in that branch. At n=N plus(S) is already outside the cut.
    # More generally just pay the exact first-interval contribution.
    return Q(1,2)/(s+N-Q(1,2))


def high_bound(d,s):
    """Rounded rational power bound; min(2, 5 exp(t)/d) integrated."""
    if d<=Q(5,2): return 2/s
    # Replace s by a modest denominator lower bound in the power only.
    a,b=(Q(int(1000*s),1000)).as_integer_ratio()
    v=Q(5,2*d)**a
    scale=10**12
    lo,hi=0,scale
    while lo+1<hi:
        i=(lo+hi)//2
        if Q(i,scale)**b>=v: hi=i
        else: lo=i
    return 2*Q(hi,scale)/(s*(1-s))-Q(5,d)/(1-s)


if __name__=='__main__':
    s=Q('0.565728')
    # p c=1-b-p, where c=number of actual children; the parent slot matters.
    for L in range(16):
        for c in range(16):
            d=1+L+c
            p=Q(1,d)
            assert c*p==1-L*p-p
            if d>1: assert c*p!=1-L*p  # plant omission of the inert slot
    print('PASS 256 exact local rank-one split rows; dropped parent slot rejected')
    # Polynomial-in-time first-discrepancy estimate must be integrable. Even
    # taking the given 3 exp(t)-t-2 at face value does not give a TV integral.
    assert s<1
    print('infection majorant has nonintegrable e^{(1-s)t} term at the target')
    print('reciprocal lower integrated d=2:',reciprocal(2,s))
    print('size-cut example first child occupation',artificial(8,s))
    # Before an unmatched root selection, the exact total root mismatch rate
    # is (2D-1)(1/D-1/(D+delta)); shared mass allocation is verified independently.
    for D in range(2,61):
        for delta in range(1,61):
            ca=Q(1,D+delta); cb=Q(1,D)
            hazard=(D-1)*(cb-ca)+delta*ca
            assert hazard==(2*D-1)*(cb-ca)<=2*Q(delta,D+delta)
            assert (D-1)*(cb-ca)<=Q(delta,D+delta)
            if delta: assert hazard>(D-1)*(cb-ca) # reject omitting private-leaf promotions
    for d in (32,64,128,256,1024):
        v=high_bound(d,s)
        tick=-(-v*10**9//1)
        print('unlimited TV bound D >=',d, f'{tick//10**9}.{tick%10**9:09d}')
    assert high_bound(128,s)<1
