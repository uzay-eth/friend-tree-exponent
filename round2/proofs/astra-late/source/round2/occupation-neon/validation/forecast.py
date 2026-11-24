"""Exact conditional forecast and Yule exit-cost replay; standard library only."""
from fractions import Fraction as Q
from itertools import product
from math import prod

LO,HI=Q('.565674722414'),Q('.565773666188')
A,B,C=Q(287589315471350942439,5*10**22),Q('-.005553926199875'),Q('-.02262119047951')
GRID=10**23


def rounded(v,upper):
    return Q(-(-v*GRID//1) if upper else v*GRID//1,GRID)


def tail(cap,upper,affine):
    # Constant is the limiting drift at x=0 with the adverse parent input.
    # V=u-v*x. Its gain increment is v/(d*(d+1)); promotion has zero delta.
    p,s=(Q(0),LO) if upper else (Q(1,2),HI)
    u=(A+B*p)/(s+3+p)
    assert u>0
    if not affine:
        if upper:return u,Q(0)
        x=Q(1,cap+1)
        return (A+B*p+C*x)/(s+3+p+x),Q(0)
    # Upper use g<=d, lower g>=0 (all arities). Sign follows from
    # low/high extrema in p,s for the same positive candidate V.
    x=Q(1,cap+1)
    v=(u-C)/(s+3+p)
    if upper:
        # desired surplus / x = u-C - v*(s+3+p+x+1/(1+x)).
        # x+1/(1+x) increases, so choose slope at the cap, not at zero.
        v=(u-C)/(s+3+p+x+1/(1+x))
        assert u-C-v*(s+3+p+x+1/(1+x))==0
    else:
        assert u-C-v*(s+3+p)==0
    assert v>0 and u-v*x>0 and B<0
    return u,v


def solve(cap,upper,affine):
    u,v=tail(cap,upper,affine)
    vals={}
    checks=0
    for d in range(cap,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1;b=Q(l,d)
            G=vals[d+1,k] if d<cap else u-v/(cap+1)
            P=vals[d,k+1] if l else 0
            rows=[]
            for s,p,g in product((LO,HI),(Q(0),Q(1,2)),(l+Q(1,2),l+Q(k+1,2))):
                den=s+3+Q(1,d)+p+g+b
                rhs=A+B*p+C/d+g*G+b*P
                rows.append((den,rhs))
            choose=max if upper else min
            out=rounded(choose(rhs/den for den,rhs in rows),upper)
            assert all((den*out>=rhs if upper else den*out<=rhs) for den,rhs in rows)
            vals[d,k]=out
            checks+=len(rows)
    return vals,checks


def exit_mass(cap,s,delta):
    return prod((n+delta)/(n+delta+s) for n in range(2,cap+1))


def cost():
    for cap in (12,20,64,256,1024):
        a,b=exit_mass(cap,HI,Q(-1,2)),exit_mass(cap,LO,Q(1,2))
        assert 0<a<b<1
        print(f'crown cap {cap}: discounted exit mass [{float(a):.9g},{float(b):.9g}]')
    # Take constant +/- E terminal continuation (width 2E) and ask for
    # total width <= T after multiplying by the stable-age factor.
    E=Q('0.001'); T=Q('2.931916e-6'); cap=10000
    # log(1+y)<=y and sum_{n=2}^K 1/(n-1/2)<=2/3+log((K-1/2)/(3/2)).
    # Use HI<3/5, exp(-2/5)>3/5 by exp(-x)>=1-x.
    # At K=10000 use (K-.5)/1.5 < 6^5, hence mass > 3/(5*6^3).
    assert HI<Q(3,5) and (cap-Q(1,2))/Q(3,2)<6**5
    floor=Q(3,5*6**3)
    assert 2*E*LO*floor>T
    # Directed fixed-point multiplication finds a much stronger necessary cap.
    # Each step is integer rounding down, so the printed obstruction is rigorous.
    hi_num,hi_den=HI.numerator,HI.denominator
    tick=10**20; high=tick; low=tick
    stop=T/(2*E*LO)
    for n in range(2,100001):
        a=(2*n-1)*hi_den; b=a+2*hi_num
        low=low*a//b
        high=(high*a+b-1)//b
        if Q(high,tick)<=stop:
            assert Q(old_low,tick)>stop
            break
        old_low=low
    else:raise AssertionError('missing threshold')
    print(f'two-sided ±.001 terminal-only error: Yule lower obstruction persists through cap {n-1}, '
          f'clears at cap {n} (necessary, not sufficient)')


def main():
    for cap in (8,16,32):
        l,n=solve(cap,False,True);u,_=solve(cap,True,True)
        l0,_=solve(cap,False,False);u0,_=solve(cap,True,False)
        assert all(l[x]<=l0[x]<=u[x]<=u0[x] for x in l)
        for d,k in ((2,0),(cap,3)):
            print(f'cap {cap}, entry D={d},k={k}, {2*n} exact corners: '
                  f'best [{float(max(l[d,k],l0[d,k])):+.10g},{float(u[d,k]):+.10g}]; '
                  f'constant [{float(l0[d,k]):+.10g},{float(u0[d,k]):+.10g}]')
    cost()


if __name__=='__main__':main()
