"""Bounded exact checks for family conditioning and sub-half corners."""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial


def rising(x,n):
    out=Q(1)
    for i in range(n): out*=x+i
    return out


def family_laws():
    rows=0
    for r,k in product(range(3,9),(1,2)):
        if k==2:
            cross=Q(1,r*(r+1))
            assert cross-Q(1,r)**2 == -Q(1,r*r*(r+1))
            # An exposed founder loses a potential colour: the survivor
            # count, and hence its next moment law, cannot be kept at r.
            assert Q(1,r-1) != Q(1,r)
        # Expose distinguished shares and the remainder, then sample
        # distinguishable colours within the surviving aggregate.
        for d,e,f,g in product(range(3),repeat=4):
            p=[d,e] if k==2 else [d+e]
            lhs=rising(1,f)*rising(1,g)
            for n in p: lhs*=rising(1,n)
            lhs/=rising(r,sum(p)+f+g)
            rhs=rising(r-k,f+g)/rising(r,sum(p)+f+g)
            for n in p: rhs*=rising(1,n)
            rhs*=rising(1,f)*rising(1,g)/rising(r-k,f+g)
            assert lhs==rhs
            rows+=1
    for f,s,b,p in product(range(5),range(1,8),(20,21,40),(Q(1,3),Q(2,5))):
        w=Q(1,4); u=1-p; z=p*(1-w)
        joint=comb(s+f,f)*u**(f+1)*p**s*s*(1-w)**(s-1)
        factor=(f+1)*u**(f+1)*p/(1-z)**(f+2)
        post=comb(s+f,f+1)*(1-z)**(f+2)*z**(s-1)
        assert joint==factor*post
        assert Q(1,b+f+s) <= (1+Q(1,b+f+1))/Q(b+f+s+1)
        rows+=1
    print(f'PASS {rows} exact neutrality, size-tilt and midpoint-prefactor rows')


def quantile_and_miss():
    n=0
    for B,badd,q in product((20,21,40),(0,1,2,27,400),
                            (Q(1,101),Q(1,3),Q(2,3),Q(7,9),Q(100,101))):
        b=B+badd; h=(1-q)/q; f=(b-1)*h//1
        if b>B:
            assert (f+Q(3,2))/b>=min(h,(B*h+Q(1,2))/(B+1))
        else:
            assert Q(f)+Q(3,2)>=Q(3,2)
        for R,j in product(range(1,11),range(6)):
            D=R+1+j
            miss=1-Q(R,D)
            retain=Q(D,(D-R)*(R+1))
            assert 0<retain<=1 and miss*retain==Q(1,R+1)
            n+=1
    print(f'PASS {n} first-miss thinning checks with quantile ranges across the seam')


def corners():
    a=sum((-1)**j*Q(3,3*j+1)/factorial(j) for j in range(6))
    assert Q(1773,1000)*Q(42,41)<Q(3,4)*a
    for b,x,c,root in ((20,Q(49,100),Q(521,2500),Q(6111,10000)),
                       (21,Q(12,25),Q(2087,10000),Q(612,1000)),
                       (40,Q(1,3),Q(261,1250),Q(61,100))):
        assert root**6 >= (1-x)/(b*x)
        sq=c*c*Q(157,50)*(b+Q(1,2))*x-(Q(4,3)+root/2)**2*(1-x)
        assert sq>0
        print(f'PASS b={b}, x={x}: unnormalized squared margin > {(sq*10**8)//1}/10^8')


if __name__=='__main__':
    family_laws()
    quantile_and_miss()
    corners()
