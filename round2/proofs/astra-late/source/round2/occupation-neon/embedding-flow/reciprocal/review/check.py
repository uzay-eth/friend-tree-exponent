"""Independent symbolic and rational mutation controls; no LP or payloads."""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial, prod

import sympy as S


def bern(poly, axes, n):
    """Exact tensor coefficients: x^j = sum_i (i)_j/(n)_j B_i^n."""
    terms = S.Poly(poly,*axes).terms()
    assert all(max(j)<=n for j,c in terms if c)
    return {i:sum(c*prod(S.Rational(comb(k,d),comb(n,d)) if k>=d else 0
                          for k,d in zip(i,j)) for j,c in terms)
            for i in product(range(n+1),repeat=len(axes))}


def contract(coefs, ys, n):
    return sum(c*prod(comb(n,i)*y**i*(1-y)**(n-i) for i,y in zip(idx,ys))
               for idx,c in coefs.items())


def closure():
    y,z=S.symbols('y z')
    H=7
    a,b=1+y/H,1+z/H
    count=0
    for m in range(1,6):
        for j,l in product(range(m+1),repeat=2):
            P=y**j*z**l
            W=a**m*b**m
            difference=S.expand(y**j*(a**(m-j)-a**m)*z**l*b**m)
            assert S.expand(difference).subs(y,0)==0
            gain=S.cancel((H/y-3+z/H)*difference)
            promo=S.expand((1-3*y/H)*W*P)
            for g in (gain,promo,difference):
                assert max(S.Poly(g,y,z).degree_list())<=2*m+1 or g==0
            assert S.cancel(W*(P.subs(y,y/a)-P)-difference)==0
            count+=1
        # One simultaneous two-high gain read checked through the positive coordinates.
        delta=S.expand(y**m*(1-a**m)*z**m*b**m)
        G=S.cancel((H/y-3+z/H)*delta)
        ys=(Q(7,12),Q(7,19))
        assert contract(bern(G,(y,z),2*m+2),ys,2*m+2)==G.subs(dict(zip((y,z),ys)))
        assert S.Poly(G,y,z).degree(z)==2*m+1
    # Separate gain terms have a real pole for the constant test: only the difference closes.
    assert S.denom(S.cancel(7*a/y)).has(y)
    # Maximal stored test degree is not a legal projection/generator test degree.
    assert S.Poly(a**3*y**8,y).degree()==11>8
    return count


def physical():
    H,m=7,3
    # Independent direct physical event list for a fixed tracked edge.
    def value(d,k,e,j):
        x=Q(H,d) if d>=H else Q(2,d)
        y=Q(H,e) if e>=H else Q(1,e)
        return (k+1)*x**2*(j+1)*y**3

    def rhs(d,k,e,j,far_p,far_x,mutation=None):
        ld,le=d-k-1,e-j-1
        base=value(d,k,e,j)
        events=((value(d+1,k,e,j),ld+Q(1,e)+sum(far_p)),
                (value(d,k,e+1,j),le+Q(1,d)+sum(far_x)),
                (value(d,k+1,e,j),Q(ld,d)),
                (value(d,k,e,j+1),Q(le,e)))
        direct=sum((next-base)*r for next,r in events)
        Y,Z=S.symbols('Y Z')
        axes=[]; vals=[]; W=1
        for D,axis in ((d,Y),(e,Z)):
            if D>=H:
                axes.append(axis); vals.append(Q(H,D)); W*=(1+axis/H)**m
        if d<H:
            Y=S.Rational(2,d)
        if e<H:
            Z=S.Rational(1,e)
        x1=(Y/(1+Y/H) if d>=H else S.Rational(1 if d==H-1 else 2,d+1)*(H if d==H-1 else 1))
        y1=(Z/(1+Z/H) if e>=H else S.Rational(H if e==H-1 else 1,e+1))
        if mutation=='seam':
            x1=S.Rational(2,d+1)
        if mutation=='shift':
            y1=Z/(1+Z/(H+1))
        p0=(k+1)*Y**2*(j+1)*Z**3
        Dp=(k+1)*x1**2*(j+1)*Z**3-p0
        Dx=(k+1)*Y**2*(j+1)*y1**3-p0
        rp=Y/H if d>=H else Q(1,d)
        rx=Z/H if e>=H else Q(1,e)
        mass=S.cancel(W*((1/rp-k-1+rx)*Dp+(1/rx-j-1+rp)*Dx
                        + (1-(k+1)*rp)*p0/(k+1)+(1-(j+1)*rx)*p0/(j+1)))
        f1=S.cancel(W*Dp); f2=S.cancel(W*Dx)
        sub=dict(zip(axes,vals))
        denom=S.sympify(W).subs(sub)
        read=lambda P:contract(bern(P,axes,2*m+2),vals,2*m+2)/denom
        if mutation=='port-weight':
            return direct,read(mass)+(read(f1)*sum(far_p)+read(f2)*sum(far_x))*denom
        return direct,read(mass)+read(f1)*sum(far_p)+read(f2)*sum(far_x)

    fp,fx=[Q(1,17),Q(1,8)],[Q(1,23)]
    for d,e in ((6,7),(7,6),(7,7),(11,19),(23,9)):
        assert len(fp)==2 and len(fx)==1
        a,b=rhs(d,2,e,1,fp,fx)
        assert a==b
    for d,e,mut in ((6,7,'seam'),(11,19,'shift'),(11,19,'port-weight')):
        try:
            a,b=rhs(d,2,e,1,fp,fx,mut)
        except S.PolynomialError:
            assert mut=='shift'  # wrong shift can prevent cancellation entirely
        else:
            assert a!=b,mut
    # Newborn edge source uses root storage even when eventual edge has two high slots.
    source=S.Rational(1,2)*(1+S.Rational(1,2))**m
    assert source/(1+S.Rational(1,2))**m==S.Rational(1,2)
    assert source*(1+S.Rational(1,2))**m/(1+S.Rational(1,2))**m!=S.Rational(1,2)


def projection_and_repair():
    y,z=S.symbols('y z')
    H,m,n=7,3,8
    W=(1+y/H)**m; new=(1+z/H)**m
    B=bern(W*y**3,(y,),n)
    E=bern(W*new*y**3,(y,z),n)
    R=bern(W*new*y**3*z/H,(y,z),n)
    wrong=bern(W*y**3*z/H,(y,z),n)
    k=Q(7,8)
    pool=[(Q(1,3),Q(7,9)),(Q(1,6),Q(7,14))]
    w=W.subs(y,k)
    base=contract(B,(k,),n)/w
    extension=lambda coef:sum(p*contract(coef,(k,i),n)/(w*new.subs(z,i)) for p,i in pool)
    assert base==k**3 and extension(E)==base/2
    exact_flux=extension(R)+Q(1,30)*base/2
    assert extension(R)<=exact_flux<=extension(R)+(base-extension(E))/H
    assert extension(wrong)<extension(R)
    assert not exact_flux<=extension(R)  # dropping positive omitted flux fails
    # General rational group repair on a nonoptimal dual, with a split group.
    x=(Q(1,4),Q(3,4),Q(1,2)); c=(Q(1),Q(-2),Q(3))
    ydual,zdual=Q(1,5),Q(-1,7)
    # A=[1,1,0], b=1; G=[0,0,1], h=2/3.
    r=(c[0]-ydual,c[1]-ydual,c[2]-zdual)
    bound=ydual+Q(2,3)*zdual-max(0,-r[0],-r[1])-Q(2,3)*max(0,-r[2])
    assert 0>=bound and sum(a*b for a,b in zip(x,c))>=bound
    # Occupancy budget: E (k)_j at an Exp(s) age for dominating Poisson(k).
    s=Q(565673,10**6)
    assert sum(comb(2,j)*factorial(j)/s**j for j in range(3))>1


if __name__=='__main__':
    print('PASS symbolic high-shift cancellation and degree counts:',closure(),'monomial rows')
    physical()
    print('PASS five exact physical edge states; rejected wrong seam, shifted gain and flux weight')
    projection_and_repair()
    print('PASS weighted projection, omitted reciprocal envelope, newborn source weight and rational dual repair')
