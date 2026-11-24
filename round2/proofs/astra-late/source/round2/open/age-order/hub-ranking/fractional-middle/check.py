"""Outward rational fractional W-moment bounds and a signed ledger."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import factorial

S = 10**13


def down(x):
    return Q(x*S//1,S)


def up(x):
    return -down(-x)


def middle(lo, hi, D=40, K=96):
    assert 0 <= lo <= hi <= 1
    P = [Q(1)]
    for j in range(1,D):
        P.append(P[-1]*(1+lo/j))
    e_inv = sum(Q((-1)**j,factorial(j)) for j in range(18))
    after = [0]+[int(S*down(e_inv*sum(P[j-1-k]/factorial(k) for k in range(j))))
                 for j in range(1,D+1)]
    q,z = hi.numerator,hi.denominator
    for d in range(D-1,1,-1):
        rescue = 0
        for k in range(K-1,1,-1):
            g = k*(d*(k-1)+1)
            rescue = z*(d*S+g*rescue)//(z*(d+g)+d*k*q)
        row = [rescue*after[1]//S]
        for j in range(1,d+1):
            row.append(j*z*(d*after[j+1]+row[-1])//(d*q+j*(d+1)*z))
        after = row
    return Q(after[1],S)


@cache
def gamma_upper(a,lo,hi,terms=48):
    """Euler product plus directed tail; a is in [1/2,3/2]."""
    if a <= 1:
        # Each factor is decreasing in p for a<=1. The tail is <=1.
        p = lo
    else:
        p = hi
    v = Q(1)
    for k in range(terms):
        v = up(v*(1+p*(a-1)/((k+1)*(k+a+p))))
    if a > 1:
        x = p*(a-1)/terms
        v = up(v/(1-x))
    return v


def tips(lo,hi,K=64):
    def cap(k):
        a = 1+Q(1,k)
        large = gamma_upper(a,lo,hi)
        small = (gamma_upper(a-Q(1,2),lo,hi)+gamma_upper(1+Q(1,k+1),lo,hi))/2
        return (1-1/(hi+a+1))*large+small/(lo+a+1)
    U = {K:up(cap(K)/(1+K*lo))}
    for k in range(K-1,1,-1):
        U[k] = max(U[k+1],up((cap(k)/k+U[k+1])/(1+lo+Q(1,k))))
    return U


def side_bounds(lo,hi,M):
    U = tips(lo,hi)
    L = {32:Q(0)}
    for k in range(31,1,-1):
        L[k] = down((M/k+k*L[k+1])/(hi+k+Q(1,k)))
    return (lambda k: U.get(k,1/(lo*k))), (lambda k: L.get(k,Q(0)))


def ledger(lo,hi,M):
    U,L = side_bounds(lo,hi,M)
    g = M-U(2)
    # A skeleton baseline, core-degree transfer, regular and terminal chain bills.
    A,B,F,FT = map(Q,('.0085','.008','.004','.0131'))
    D = 32
    u = 1/(lo*(D+1))
    z = Q(D,D+1)*g+(g-U(2))/2
    if min(B-u/(D+1),L(2)-u,z-(2*A+B),
           (A-F)*D/2-(2*A+B),g,Q(D,D+1)*(g-u)) <= 0:
        return False,'tail or seed'

    def profile(d,c,r):
        l = d-c
        w = (l+r)/(d+1)
        s,t = w/d,Q(1,2*d)
        credit = (w+l*t)*g-l*(max(0,t-s)*max(0,U(2)-L(d+1))
                             +max(0,s-t)*max(0,U(d+1)-L(2)))
        return s,credit

    for d in range(2,D):
        if profile(d,0,Q(0))[1]<=0:
            return False,('star',d)
        for c in range(1,d+1):
            if d==c==2:
                continue
            bill = FT if d==2 else F
            for k in ((0,c) if U(d+1)/(d+1)<=B else range(c+1)):
                for r in (Q(k,2),Q(c,2)):
                    s,C = profile(d,c,r)
                    if C+A*(c-2)-B-k*bill-(c-k)*max(0,s*U(d+1)-B)<0:
                        return False,('retained',d,c,k,r)

    # Scale the chain endpoints outwards to bounded integers, then compare
    # the same signed port bill on lengths 1, 2, 3.
    X = 6*10**9
    Y = 10**9
    end = []
    def endpoint(recip,lower,s,own,bill):
        end.append((recip*X/6, max(0,U(3)-lower),s,own,bill))
    for d in range(2,10):
        for c in range(1,d+1):
            if d==c==2:
                continue
            for r in (Q(1,2),Q(c,2)):
                s,_ = profile(d,c,r)
                endpoint(Q(1,d),L(d+1),s,U(d+1),FT if d==2 else F)
    for d in range(10,31):
        endpoint(Q(1,d),L(d+1),None,U(d+1)/(d+1),F)
    for r in (Q(0),Q(1,31)):
        endpoint(r,0,None,1/(lo*32**2),F)
    ef = []
    def ceil(x): return -((-x)//1)
    for r,f,s,own,bill in end:
        ef.append((r//1,ceil(r),ceil(Y*f),None if s is None else (X*s//1,ceil(X*s)),
                   ceil(Y*own),ceil(Y*max(0,own-L(3))),bill*X*Y//1))

    @cache
    def port(i,j,e):
        _,_,f,t,own,rev,bill = e
        return (max(j*f,X*own) if t is None else max(0,(j-t[0])*f,(t[1]-i)*rev))-bill

    slope = ceil(Y*max(0,U(3)-L(3)))
    credit = 2*g*Y//1
    for length in (1,2,3):
        for x,y in product(ef,repeat=2):
            if length==1:
                Ls = ((x[0]+y[0],x[1]+y[1]),)
            else:
                Ls = ((X//12+x[0],X//12+x[1]),)+((X//6,X//6),)*(length-2)+((X//12+y[0],X//12+y[1]),)
            total = port(*Ls[0],x)+port(*Ls[-1],y)
            total += slope*sum(max(j-k,l-i) for (i,j),(k,l) in zip(Ls,Ls[1:]))
            total -= credit*sum(i for i,_ in Ls)
            if total>B*X*Y:
                return False,('chain',length)
    return True,None


if __name__ == '__main__':
    for j in range(256):
        lo,hi = Q(j,256),Q(j+1,256)
        m = middle(lo,hi)
        assert m>Q(55,100)
        if j>=32:
            assert m>tips(lo,hi)[2],(lo,hi)
    print('PASS 256 boxes: physical serial middle > .55 on 0<p<=1',flush=True)
    print('PASS serial-only margin on 224 boxes covering [1/8,1]',flush=True)
    for j in range(473,512):
        lo,hi = Q(j,512),Q(j+1,512)
        ok,where = ledger(lo,hi,middle(lo,hi))
        assert ok,(lo,hi,where)
    print('PASS signed whole-tree margin on 39 boxes covering [473/512,1]')
