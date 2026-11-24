#!/usr/bin/env python3
"""Independent exact checks for the original Q3 fixed-threshold comparison."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial


def rates(par):
    adj=[[] for _ in par]
    for y,x in enumerate(par[1:],1):adj[x].append(y);adj[y].append(x)
    d=list(map(len,adj));r=[Q(0)]*len(par)
    for y,nn in enumerate(adj):
        for x in nn:r[x]+=Q(1,d[y])
    assert sum(r)==len(par)
    L=[sum(d[y]==1 for y in nn) if d[x]>1 else 0 for x,nn in enumerate(adj)]
    return adj,d,r,L


def entrance():
    rows=0
    for n in range(3,8):
        for suf in product(*(range(j) for j in range(3,n))):
            p=(-1,0,0)+suf;adj,d,r,L=rates(p)
            actual=Q(0);hostile=0
            for x in range(n):
                post=p+(x,);_,dn,rn,Ln=rates(post)
                assert rn[n]==Q(1,d[x]+1)
                promoted=post+(n,);_,de,_,Le=rates(promoted)
                assert de[n]==2
                enough=Le[x]>=1
                if not enough:hostile+=1
                else:actual+=r[x]*rn[n]/(n*(n+1))
                # In a live postpromotion state the old core contributes <=1/2.
                assert rates(promoted)[2][n]<=Q(3,2)
                assert (enough and d[x]>1)==bool(L[x])
            birth=sum(Q(l,d0) for l,d0 in zip(L,d))
            assert actual>=Q(2,3)*birth/(n*(n+1))
            assert sum(r[x]/(d[x]+1) for x in range(n))<=Q(3,2)*sum(z>1 for z in d)
            assert actual==sum(r[x]/(d[x]+1) for x in range(n) if L[x])/(n*(n+1))
            assert actual>=Q(1,3*n*(n+1))
            assert hostile
            rows+=1
    # A birth onto an old leaf, and a birth onto a bare middle of a path,
    # both leave the exterior with no tagged leaf after immediate promotion.
    for p,x in (((-1,0,0),1),((-1,0,0,1,2),0)):
        n=len(p);assert rates(p+(x,n))[3][x]==0
        restored=p+(x,x,n)
        assert rates(restored)[3][x]==1
    print('PASS',rows,'exact labelled entrance/promotion rows; leaf-target and bare-parent hostile plants')


def bayes_and_phantoms():
    rows=0
    for a,b in ((1,4),(3,7),(9,28)):
        for j,k,u in product(range(1,8),range(1,12),(Q(1,7),Q(9,5))):
            if k<j:continue
            nb=comb(k-1,j-1)*Q(a,b)**j*(1-Q(a,b))**(k-j)
            lr=Q(b**k,a**j)*u**(k-j)*Q(factorial(j-1),factorial(k-1))
            assert nb*lr==((b-a)*u)**(k-j)/factorial(k-j)
            rows+=1
    for n in range(1,30):
        for h in range(n+1):
            assert Q(h,n)*Q(h+1,n+1)+Q(n-h,n)*Q(h,n+1)==Q(h,n)
    print('PASS',rows,'exact conditional Yule transitions and 464 phantom-fraction rows')


def denominator():
    x=Q(2,5);a=x/(1-x);D=Q(7,3)
    for S,U in product(map(Q,(0,3,19)),map(Q,(1,7,40))):
        # Suffices without independence and without substituting the
        # proposal size for the actual numerator or side denominator.
        if U>a*S+(1+a)*D:
            W=U-D;Smark=W+D
            assert W>x*(Smark+S)
    for x in (Q(1,10**6),Q(1,2),1-Q(1,10**6)):
        for n in (3,17,50):
            # Laplace normalizer of a gamma built from n independent
            # exponentials; alpha S+ after tilt scales as x Gamma(n,1).
            a=x/(1-x)
            assert (1+a)**(-n)==(1-x)**n and a/(1+a)==x
            assert 1/(1+a*(1-x))>=Q(1,2)
        # A bounded oscillation can preserve both envelopes and destroy sign.
        q=1-x
        assert all(q*q**v < q**(v+1)/q for v in (4,6,8))
    print('PASS actual-ratio, auxiliary tilt, endpoint scales and adjacent-oscillation controls')


def near_one():
    widths=[]
    for v in (3,4,8,16,32,64):
        A=Q((v+1)**2,v)
        # Bound Γ(v+2)/(Γ(3/2)Γ(v+3/2)) using π>3.
        u=Q(2*4**(v+1)*factorial(v+1)**2,3*factorial(2*v+2))
        lower=Q(factorial(v),v*(v-1)*(v+1))/A**v
        upper=u/(2*(v+1))
        w=min(Q(1,2),(lower/(2*upper))**2)
        assert w*upper**2 <= (1-w)**2*lower**2
        widths.append((v,float(w)))
        for q in (Q(1,5),Q(1,10**6),w):
            n=-(-(A/q)//1)
            assert (n-1)*q<A and 1-Q(v+1,n)/q>=Q(1,v+1)
    print('PASS explicit elementary near-one widths',widths)


if __name__=='__main__':
    entrance()
    bayes_and_phantoms()
    denominator()
    near_one()
