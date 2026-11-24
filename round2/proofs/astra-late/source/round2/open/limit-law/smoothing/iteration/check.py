#!/usr/bin/env python3
"""Exact finite controls for retry rates, angular curvature and stock gates."""
from fractions import Fraction as Q
from math import comb, factorial


def retry():
    # Multiply the two series q and q^2 without using the closed coefficient.
    for j in range(1,65):
        assert Q(1,2**(j+1))+sum(Q(1,2**(i+1))*Q(1,2**(j-i+1)) for i in range(1,j)) == Q(j+1,2**(j+2))
    r=Q(13,200)
    a=[Q(j+1,2**(j+2)) for j in range(1,9)]
    F=Q(3,4)+r*(sum(v/(Q(j,2)-r) for j,v in enumerate(a,1))+(Q(3,4)-sum(a))/(Q(9,2)-r))
    rho=(1+(1+2*F)/(3*(1-r)))/(2-r)
    assert F==Q(5715644171206736442597875,7046497926887658264039168)
    assert rho==Q(11951004048194066431207745075,11951906448537490802583997242)<1
    alpha=1/(1+4*r); c=r*alpha
    assert (1-alpha)/2-c==c==Q(13,252)
    assert c/(1+2*c)==Q(13,278)
    for mu in (Q(1,2),Q(3,5),Q(99,100),Q(1)):
        x=13/(278+13*mu)
        assert x==Q(1,2)-x/(2*c)-mu*x/2
        assert (1/(1-x)>=Q(291,278)) and (mu==1 or 1/(1-x)>Q(291,278))
    for m in range(2,13):
        assert sum(Q(comb(m,i),2**m) for i in range(m+1))==1
    print('PASS 64 defective-tail coefficients; exact retry bound',rho,'and all-set powers')


def kernel():
    count=0
    for k in (Q(1),Q(113,100),Q(6,5),Q(3,2),Q(2)):
        for p in (Q(1,500),Q(1,10),Q(1,2),Q(9,10),Q(499,500)):
            # Verified separately on sin^k and cos^k after dividing
            # by that monomial. c=cot(theta), t=tan(theta), c*t=1.
            for x2,y2 in ((1-p,p),(p,1-p)):
                g2=k*(k-1)*x2/y2-k
                transport=k*(k-1)*(x2/y2-1)
                assert g2-transport==k*(k-2)
                count+=1
    maxima=0
    for a in range(2,33):
        for b in (2,3,4,7,16,32,65,128,256):
            n=a+b;p=Q(2*a-1,2*(n-1))
            # Exact squared angular density at its mode, divided by N.
            peak2=4*Q(factorial(n-1),factorial(a-1)*factorial(b-1))**2*p**(2*a-1)*(1-p)**(2*b-1)
            assert peak2<n
            maxima+=1
    print('PASS',count,'angular curvature controls and',maxima,'integer Beta peak controls')


def stock():
    n=0
    for q in (Q(1,3),Q(1),Q(2),Q(7,2),Q(10)):
        for l in range(5,102):
            if l<4+2*q: continue
            plus=-q/(l-2)
            minus=q/(l-3-q)
            assert l*plus+minus<=-q and 1+plus>0
            n+=1
    print('PASS',n,'inverse-stock generator controls')


if __name__=='__main__':
    assert __debug__
    retry(); kernel(); stock()
