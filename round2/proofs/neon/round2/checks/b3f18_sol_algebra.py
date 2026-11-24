"""Exact referee identities. No imports from the certificate implementation."""
from fractions import Fraction as Q
from itertools import product
from math import comb, prod
from random import Random


def B(a, t):
    n = len(a) - 1
    return sum(comb(n,j)*t**j*(1-t)**(n-j)*v for j,v in enumerate(a))


def poly(a, t):
    return sum(v*t**i for i,v in enumerate(a))


def controls(a,n,tmax):
    return [sum(Q(comb(j,i),comb(n,i))*a[i]*tmax**i
                for i in range(min(j+1,len(a)))) for j in range(n+1)]


def restrict(a,m):
    out = [a[0]]
    while len(a)>1:
        a = [(1-m)*x+m*y for x,y in zip(a,a[1:])]
        out.append(a[0])
    return out


def c_hull4(a,m):
    n,d=m.numerator,m.denominator
    c0,c1,c2,c3,c4=12*d**4,3*d**3*n,2*d*d*n*n,3*d*n**3,12*n**4
    A,B,C,D,E=a
    return [v/c0 for v in (c0*A,c0*A+c1*B,c0*A+2*c1*B+c2*C,
            c0*A+3*c1*B+3*c2*C+c3*D,c0*A+4*c1*B+6*c2*C+4*c3*D+c4*E)]


def tail(d,c,F,G,credit,ri,re,rb,Y,transfer):
    f0,f1,f2,f3 = F
    g1,g2,g3 = G[1:]
    e0,e1,e2,e3 = credit
    x = c+1-ri-re
    a,b,p = f0-f1,f1-f2,f2-f3
    u = 3*c*((g1-f1)-(g2-f2))/d
    v = (e0-3*e1+3*e2)/d
    A = g3-f3+u+v
    H = -Q(c,d)*(g3-f3)-u+e3/d-v
    # Transcribe the three additive coefficient blocks in C, cancelling M.
    cp = [3*d**3, 3*d*d*(2-x), d*(4-3*x), 1-x, 0]
    ca = [3*d*d*(x+1)+3*d*x, -3*d*d*(x+1)+6*d+6*x,
          -3*d*x-5*d-6*x+3+3*x/d, -x-2-3*x/d, 0]
    cb = [-6*d*d*x-3*d*d, 6*d*d*(x+1)-9*d*x-6*d,
          6*d*x+10*d-3*x-3, 2*x+4, 0]
    P = [s+p*i+a*j+b*k for s,i,j,k in zip(
        [d**3*A,d*d*(d*H+3*A),d*(3*d*H+3*A),3*d*H+A,H],cp,ca,cb)]
    one = f0+g1-2*f1+e0/(3*d)
    two = 2*f1+g2-3*f2-Q(c,d)*(g1-f1)-(c+1-ri-(rb if transfer else re))*a/d+(3*e1-e0)/(3*d)
    return one,two,P, (rb-re)*Y if transfer else Q(0)


def overflow(d,F):
    f0,f1,f2,f3 = F
    a,b,p = f0-f1,f1-f2,f2-f3
    return a,2*b,[3*d*d*p+3*d*(a-b),3*d*p+(3-3*d)*a+(6*d-3)*b,p-2*a+4*b]


def main():
    rng = Random(669523)
    vec = lambda n: [Q(rng.randrange(-31,32),rng.randrange(1,32)) for _ in range(n)]
    count=0
    for d,c,n,transfer in product((4,10,18,23), (*range(19),24), range(3), (False, True)):
        if c < n: continue
        m=min(Q(1),Q(d,c)) if c else Q(1)
        den=[d**3,3*d*d,3*d,1]
        # Unit directions cover every F, bud, and credit variable; Y_ref tested separately.
        for axis in range(12):
            a=[Q(i==axis) for i in range(12)]
            F=a[:4]; G=[F[0]]+a[4:7]; credit=a[7:11] if n else [Q(0)]*4
            Y=a[11]; rb=Q(c-n,7)
            for ri,re in product((Q(0),Q(n,2)),(Q(0),Q(c-n,2))):
                one,two,P,shift=tail(d,c,F,G,credit,ri,re,rb,Y,transfer)
                bp=controls(P,4,m); bd=controls(den,4,m)
                assert c_hull4(list(map(Q,P)),m)==bp and c_hull4(list(map(Q,den))+[Q(0)],m)==bd
                assert min(bd)>0
                # Eight distinct nonzero t's: after clearing denominators the
                # residual is a polynomial of degree <= 7, hence identically zero.
                for i in range(1,9):
                    t=m*i/8; s=1-t; D=d/t; z=d/(D+1)
                    lhs=(D-c+ri+re)*(B(F,z)-B(F,t))+(1-c/D)*(B(G,t)-B(F,t))+B(credit,t)/D
                    if transfer: lhs-=(re-rb)*(3*t*t*s*(F[0]-F[1])/d+t**3*Y)
                    assert lhs==3*t*s*s*one+3*t*t*s*two+t**3*(poly(P,t)/poly(den,t)+shift)
                    assert B(bp,t/m)==poly(P,t) and B(bd,t/m)==poly(den,t)
                    ratios=[u/v+shift for u,v in zip(bp,bd)]
                    assert min(ratios)<=poly(P,t)/poly(den,t)+shift<=max(ratios)
                    count+=1
        for _ in range(3):
            F=vec(4); credit=vec(4) if n else [Q(0)]*4
            aa,bb,P=overflow(d,F); den2=[d*d,2*d,1]
            bp,bd=controls(P,2,Q(1)),controls(den2,2,Q(1))
            assert [2*v for v in bd]==[2*d*d,2*d*d+2*d,2*d*d+4*d+2]
            assert [2*v for v in bp]==[2*P[0],2*P[0]+P[1],2*sum(P)]
            ratios=[Q(0)]+[u/v for u,v in zip(bp,bd)]
            x0,x1,x2,x3=credit; A=x0/(3*d); H=(3*x1-x0)/(3*d)
            lohi=[(x0-3*x1+3*x2)/d,x3/d]
            for i in range(1,9):
                t=m*i/8; s=1-t; D=d/t; z=d/(D+1)
                rho=Q(i,8)*c/2; rate=D-c+rho
                if c: assert 0<=rate<=D-c/2<D+1
                r=poly(P,t)/poly(den2,t); ell=(1-t)*lohi[0]+t*lohi[1]
                assert (D+1)*(B(F,z)-B(F,t))==3*t*s*s*aa+3*t*t*s*bb+t**3*r
                assert B(credit,t)/D==3*t*s*s*A+3*t*t*s*H+t**3*ell
                alpha=rate/(D+1)
                assert min(ratios)<=alpha*r<=max(ratios) and min(lohi)<=ell<=max(lohi)
    print(f'tail coefficients and hulls: {count} exact nonoverflow basis identities; overflow and inverse-degree credit passed')

    for d,D in product((4,10,18),(18,19,51,14018,14181,10**18)):
        t=Q(d,D); c2,c3=vec(2)
        assert B([0,Q(1,3*d),Q(2,3*d),Q(1,d)],t)==Q(1,D)
        assert B([0,0,c2/(3*d*d),c2/d**2+c3/d**3],t)==c2/D**2+c3/D**3
    for h,k,LF in product(range(2),(0,1,17,18,180,181,182,300),(4000,14000)):
        d=18; c=k+h; kf=180
        F=vec(4); c2,c3=vec(2); R,a,lamf=vec(3)
        # a = one-body bud jump plus newborn f; take separate values at cap and above.
        for L in (LF,LF+1,LF+5,10**9,0,1) if k>kf else (LF,LF+1,LF+5,10**9):
            D=L+c; t=Q(d,D); s=D-d; b=Q(L,D)
            if k<=kf:
                m=Q(d,LF+1+c); blo=Q(LF+1,LF+1+c); x=Q(1,LF+c)
                cb=controls([0,0,c2,c3],3,x)
                assert cb==[0,0,c2*x*x/3,c2*x*x+c3*x**3]
                xs=[(R-rho)*v for rho,v in product((0,Q(c,2)),cb)] if k<=17 else [Q(0)]
                ex=(R-Q(c,3))*(c2/D**2+c3/D**3) if k<=17 else Q(0)
                assert min(xs)<=ex<=max(xs)
                if L==LF:
                    nums=[s**3,3*d*s*s,3*d*d*s,d**3]
                    assert b*B(F,t)==Q(L,D**4)*sum(f*g for f,g in zip(F,nums))
                    continue
            else:
                m=Q(d,kf+1+h); blo=Q(0); xs=[Q(0)]; ex=Q(0)
            cp=restrict(F,m)
            assert cp==[B(F[:j+1],m) for j in range(4)]
            assert B(cp,t/m)==B(F,t)
            vals=[b0*(a+q)+xx-lamf for b0,q,xx in product((blo,1),cp,xs)]
            assert min(vals)<=b*(a+B(F,t))+ex-lamf<=max(vals)
    lo=Q(565673,10**6); hi=Q(22631,40000)
    assert hi-lo==Q(102,10**6)
    lowpref=Q(146042017510,184124847315); highpref=Q(584179227080,115462837236)
    assert lowpref>Q('0.793168') and highpref<Q('5.059457')
    for lam in (lo,hi):
        vals=[prod(1+lam/j for j in range(3,n)) for n in range(3,90)]
        assert vals[0]==1 and all(v*(1+lam/j)==w for j,(v,w) in enumerate(zip(vals,vals[1:]),3))
    print('invD, f tail transfers, cap/VSTAR Bernstein controls, endpoint width and finite-size constants passed')


if __name__=='__main__':
    main()
