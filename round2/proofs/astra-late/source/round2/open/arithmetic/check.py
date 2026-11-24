#!/usr/bin/env python3
"""Small exact controls for bracket exclusions and trigger clearers."""
from fractions import Fraction as Q
from math import gcd, isqrt, prod
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'finite'))
import residual_entries as physical

physical.MARK_CAP=20

LO=Q(565662790716854568223,1000008252270910000000)
HI=Q(1698046074520368376961,3001152517940340000000)
LONG=Q(2603,4602)
ALPHA=Q(31,20)


def first_rational():
    for q in range(1,500):
        p=(LO.numerator*q+LO.denominator-1)//LO.denominator
        if p<=HI.numerator*q//HI.denominator:
            return Q(p,q)
    raise AssertionError('search cap')


def first_square(limit):
    # s=t−1/t, t=u/v>1 reduced; q=uv with no further cancellation.
    best=None
    for v in range(1,isqrt(limit)+1):
        for u in range(v+1,limit//v+1):
            if gcd(u,v)!=1:
                continue
            p,q=u*u-v*v,u*v
            if LO*q<=p<=HI*q and (best is None or q<best.denominator):
                best=Q(p,q)
    return best


def first_quadratic_height(limit):
    # The interval is positive, so solve for the integer constant term.
    hits=[]
    for a in range(1,limit+1):
        for b in range(-limit,limit+1):
            p=lambda x: a*x*x+b*x
            vals=[p(LO),p(HI)]
            v=Q(-b,2*a)
            if LO<v<HI:
                vals.append(p(v))
            lo,hi=-max(vals),-min(vals)
            for c in range(max(-limit,-(-lo.numerator//lo.denominator)),
                           min(limit,hi.numerator//hi.denominator)+1):
                D=b*b-4*a*c
                if gcd(gcd(a,b),c)==1 and D>0 and isqrt(D)**2!=D:
                    hits.append((max(a,abs(b),abs(c)),a,b,c))
    return min(hits) if hits else None


def gain(d,s):
    return Q(d-1)/(s+Q(1,2)+d-Q(1,d))


def omission(N,s):
    # Exact all-gain path to a degree-N star, reset to J20, two
    # additional intrinsic changes and one plain copy to the trigger target.
    out=Q(1)
    for d in range(2,N):
        out*=gain(d,s)
    out*=Q(N-1,N)/(s+Q(1,2)+N-Q(1,N))
    out*=Q(1,2)/(s+Q(31,20))
    # The extra gain after a promotion is the promoted child's up-send.
    out*=Q(1,2)/(s+Q(51,20))
    out*=(Q(1,2)-Q(1,3))/(s+Q(1,20)+4-Q(1,3))
    return out*Q(1,2)/(s+2)


def loop_poles():
    # A 19-transition reset cycle from J20: 18 gains then direct birth.
    assert min(Q(n)-Q(1,n)+Q(1,20) for n in range(3,50))>ALPHA
    # A fast exit J20 -> root degree 2 with core child -> gain ->
    # plain copy to J2. The mark at the child is 1/2, not 1/20.
    assert Q(1,2)-Q(1,3)==Q(1,6)
    assert Q(2)-Q(1,2)+Q(1,20)==ALPHA
    assert Q(2)-Q(1,2)+Q(1,2)>ALPHA
    for N in (22,44,85):
        assert prod(gain(d,HI) for d in range(2,N))>=Q(6,(N-1)*N*(N+1))
        assert omission(N,HI)>Q(1,10000*N**4)
    print('exact trigger loop/exit rates and omitted-path lower controls passed')


def rate(t,u):
    n,D,_=physical.statistics(t)
    p=Q(1,D)
    b=Q(t[1],D)
    value=dict(physical.intrinsic_moves(t)).get(u,Q(0))
    value+=sum((physical.mark(y[0])-p)*(u==y)+p*(u==physical.plus(y)) for y in t[2])
    if D>2 and u==physical.newborn(min(D,20)):
        value+=b
    return value,physical.mark(t[0])+n-p


def real_paths():
    j2,j20=map(physical.newborn,(2,20))
    after=(20,0,(j2,))
    grown=(20,1,(j2,))
    tail=(j20,after,grown,j2)
    assert [rate(t,u)[0] for t,u in zip(tail,tail[1:])]==[Q(1,2),Q(1,2),Q(1,6)]
    assert [rate(t,u)[1] for t,u in zip(tail,tail[1:])]==[ALPHA,Q(51,20),Q(223,60)]
    for m in (20,22,44,85):
        stars=[(2,d-1,()) for d in range(2,m+1)]
        path=stars+list(tail)
        exact=prod(r/(HI+a) for t,u in zip(path,path[1:]) for r,a in (rate(t,u),))*Q(1,2)/(HI+2)
        assert exact==omission(m,HI)
        short=stars+[j20]
        exact=prod(r/(HI+a) for t,u in zip(short,short[1:]) for r,a in (rate(t,u),))*Q(1,2)/(HI+ALPHA)
        assert exact>Q(1,m*m*(m+1)*(m+2))
    cycle=[(20,d-1,()) for d in range(2,21)]+[j20]
    assert len(cycle)-1==19
    assert all(rate(t,u)[0]>0 for t,u in zip(cycle,cycle[1:]))
    print('actual capped generator paths: 19-step cycle, exit, 4 omitted histories')


def positive_counterexample():
    r=Q(43,76)
    # F(r)=1 by a geometric telescope; the 19-label first return is
    # F/(19-18F). This is a model counterexample, not the physical K.
    for n in (1,4,20):
        f=sum(Q(r+j,r+j)*Q(1,2**j) for j in range(1,n+1))
        assert f+Q(1,2**n)==1
        assert f/(19-18*f)<1
    print('positive transcendental-characteristic control: exact rational crossing',r)


if __name__=='__main__':
    r=first_rational()
    sq=first_square(57200)
    print('least rational denominator:',r)
    print('least rational square-discriminant denominator:',sq)
    assert LONG<LO and isqrt(LONG.numerator**2+4*LONG.denominator**2)==9565
    assert LONG.numerator**2+4*LONG.denominator**2==9565**2
    assert r==Q(43,76) and sq==Q(32361,57200)
    q2=first_quadratic_height(30)
    print('least primitive quadratic height (height,a,b,c):',q2)
    assert q2==(7,3,-7,3)
    loop_poles()
    real_paths()
    positive_counterexample()
