"""Independent exact root counts, finite-field tests and denominator scan."""

from fractions import Fraction as Q
from itertools import product
from math import comb, gcd, isqrt, lcm, prod

L=Q(126905728962161948969119589499374909399,224328301684105796116575000000000000000)
U=Q(127993535256830043,226242120050000000)
CASES=[((35,-63,2),3),((-15,15,13,13),11),((-5,6,3,3,1),17),((-3,2,3,3,3,1),2)]


def counts():
    # Descartes after x=(L+U*t)/(1+t), which maps t>0 to the open interval.
    den=lcm(L.denominator,U.denominator)
    a,b=int(L*den),int(U*den)
    for f,_ in CASES:
        d=len(f)-1
        h=max(map(abs,f))
        base=[]
        for i in range(d+1):
            base.append([den**(d-i)*sum(comb(i,k)*a**(i-k)*b**k*comb(d-i,j-k)
                                     for k in range(max(0,j-d+i),min(i,j)+1))
                         for j in range(d+1)])
        def coeff(c):
            return [sum(x*base[i][j] for i,x in enumerate(c)) for j in range(d+1)]
        def sign_changes(v):
            signs=[x>0 for x in v if x]
            return sum(x!=y for x,y in zip(signs,signs[1:]))
        assert sign_changes(coeff(f))==1
        for H in (h-1,h):
            hits=[]
            for leading in range(1,H+1):
                for inner in product(range(-H,H+1),repeat=d-1):
                    c=(0,)+inner+(leading,)
                    v=coeff(c)
                    # Constant shifts all transformed coefficients by a common
                    # amount after division by the corresponding binomial.
                    q=[Q(x,base[0][j]) for j,x in enumerate(v)]
                    mn,mx=min(q),max(q)
                    start=max(-H,(-mx).__ceil__())
                    end=min(H,(-mn).__floor__())
                    for c0 in range(start,end+1):
                        assert v[0]+c0*base[0][0] and v[-1]+c0*base[0][-1]
                        changes=sign_changes([x+c0*y for x,y in zip(v,base[0])])
                        assert changes<=1
                        if changes:
                            hits.append((c0,)+c[1:])
            if H<h:
                assert not hits
            else:
                assert f in hits
                print(f'PASS degree {d}: zero roots below height {h}, {len(hits)} polynomial(s) at height {h}: {hits}')


def irreducible(f,p):
    # A d-degree polynomial is irreducible iff gcd(X^(p^i)-X,f)=1
    # for i=1,...,floor(d/2), since a reducible f has a small factor.
    def trim(a):
        a=[x%p for x in a]
        while a and not a[-1]: a.pop()
        return a

    def rem(a,b):
        a=trim(a)
        while len(a)>=len(b):
            q=a[-1]*pow(b[-1],-1,p)%p
            k=len(a)-len(b)
            for j,x in enumerate(b): a[j+k]=(a[j+k]-q*x)%p
            a=trim(a)
        return a

    assert gcd(*f)==1 and f[-1]%p
    f=trim(f)
    for i in range(1,(len(f)-1)//2+1):
        a,b=f,rem([0,-1]+[0]*(p**i-2)+[1],f)
        while b: a,b=b,rem(a,b)
        if len(a)!=1: return False
    return True


def discriminants():
    assert Q(99,175)<L<Q(241,426)<U<Q(142,251)
    ordinary=[]
    square=[]
    for q in range(1,94252):
        for p in range((L*q).__ceil__(),(U*q).__floor__()+1):
            if gcd(p,q)!=1: continue
            if q<=426: ordinary.append(Q(p,q))
            z=p*p+4*q*q
            if isqrt(z)**2==z: square.append(Q(p,q))
    assert ordinary==[Q(241,426)] and square==[Q(53320,94251)]
    print('PASS direct reduced-denominator/isqrt scan through 94251')


def tails():
    # The integral lower bound for the block n+1 through 2n is uniform
    # since (n+1)/(n+3)>=21/23 and (n+1)/(2n+3)<1/2.
    assert (Q(21,23)**3-Q(1,2)**3)/3>Q(1,16)
    assert Q(10,51)>Q(1,6) and L>Q('0.56564')
    for m in range(21,161):
        star=prod(Q(d-1,d+2) for d in range(2,m))
        assert star==Q(6,(m-1)*m*(m+1))
        assert star*Q(m-1,m*(m+2))*Q(1,6)==Q(1,m*m*(m+1)*(m+2))
    r=Q(241,426)
    c=1-1/(2*(r+2))
    for n in (1,2,17,100,1001,94251):
        e=c/(n+1)**2
        g=gcd(c.numerator,(n+1)**2)
        assert (e.numerator,e.denominator)==(c.numerator//g,c.denominator*(n+1)**2//g)
    print('PASS uniform physical block inequality and exact toy denominator reductions')


if __name__=='__main__':
    counts()
    for f,p in CASES:
        assert irreducible(f,p)
    print('PASS all four irreducibility witnesses by Frobenius gcd')
    assert irreducible((-5,6,6,-4,4),11)
    print('PASS second quartic attainer irreducible modulo 11')
    discriminants()
    tails()
