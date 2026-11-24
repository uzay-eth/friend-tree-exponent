#!/usr/bin/env python3
"""Exact compact-mark interpolation constants and reassignment controls."""
from fractions import Fraction as Q
from math import factorial

D=Q(73,200)
RHO=Q(53,100)
GAMMA=RHO-D
S=Q('0.56572365')


def interp_bounds(s, M, m, hmin):
    a=Q(1,M)
    # The root toll is <=1. The other forcing terms are bounded by
    # their total weight times ||U/W||; the proof is in UNCAPPED.md.
    k=s-GAMMA
    # For the accepted table, P_i >= z_i and P_i <=1.
    # Compare each copy/dormant source to its own depth-1 or depth-2 V arm.
    emission=1/hmin+(2+Q(1,2)/D)/D
    F0=(Q(1,2)+D)/(hmin*(s-RHO))
    force=1/hmin+emission*F0
    return a/k, force/k*(a/k)**m


def linear_rem(a,b,x,n,k):
    f=lambda u: Q(factorial(n),1)/(k+u)**(n+1)
    v=((b-x)*f(a)+(x-a)*f(b))/(b-a)-f(x)
    lim=Q(factorial(n+2),1)*(x-a)*(b-x)/(2*(k+a)**(n+3))
    assert 0<=v<=lim


def main():
    for j in range(2,80):
        a,b=Q(1,j+1),Q(1,j)
        for n in range(4):
            linear_rem(a,b,(a+b)/2,n,S-GAMMA)
    print('checked 312 exact Laplace interpolation inequalities')
    h=Q(86271979,10**9)
    for M,m in ((20,8),(20,12),(20,16),(36,8),(36,12)):
        scale,err=interp_bounds(S,M,m,h)
        print('M',M,'terms',m,'ratio',float(scale),'analytic row bound',float(err))


if __name__=='__main__': main()
