"""Marked cut identity and one shared positive all-copy comparison."""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
from argparse import ArgumentParser
from functools import cache
import sys

R2=Path(__file__).resolve().parents[7]
sys.path.insert(0,str(R2/'occupation-neon'/'ratio-transfer'))
from replay import Pair, Cert, shape, tree
sys.path.insert(0,str(R2/'finite'))
import residual_entries as phy
phy.MARK_CAP=20
S=Q('.56572365')


def expo(x,n=24):
    assert 0<=x<n+2
    return sum(x**j/factorial(j) for j in range(n+1))+x**(n+1)/factorial(n+1)/(1-x/(n+2))


def source_budget():
    # Use the audited post-EDGE prior, separately from the pointwise test.
    L,U=Q('.565678056632067'),Q('.565770469400619')
    sys.path.insert(0,str(R2/'finite'/'gauge'))
    from physical_trigger_bounds import blocks,UNIT
    assert sum(Q('1.402')**j/factorial(j) for j in range(17))>Q('4.062')
    def anchor(t):
        r=expo(Q('1.402')*(U-t)/(U-Q('.4481')))
        m=Q(blocks(U,r,16,tuple(range(2,21)))[2],UNIT)
        theta=r-m/(2*(t+Q(1,20)+Q(3,2)))
        assert theta<1
        return 1+(r-1)/(1-theta)
    f=anchor(S)
    rows=[]
    for t in map(Q,('.5646','.565')):
        a=anchor(t)
        rows.append((1-(S-L)*(a-1)/(L-t),t))
    l,t=max(rows)
    print('prior-assisted source interval',float(l),float(f),'lower anchor',t,flush=True)
    return l,f


def degree_weight():
    # Reconstruct the published decreasing .19 degree weight. The source
    # comparison below also verifies the finite and saturated inequalities.
    d,q=Q('.19'),Q('.272');D=20;unit=10**9
    z=lambda n:1/(1+d*n)
    from math import ceil
    h={(D,c):Q(ceil(unit*z(D)/q),unit) for c in range(1,D+1)}
    for n in range(D-1,1,-1):
        for c in range(n,0,-1):
            L=n-c;b=Q(L,n)
            v=(L*h[n+1,c]+b*(h.get((n,c+1),0)+z(n)))/(q+L+b)
            h[n,c]=Q(ceil(unit*max(Q(1,4),h[n+1,c],v)),unit)
    h[2,1]=Q(1)
    for n in range(2,D):
        for c in range(1,n+1):
            assert h[n+1,c]<=h[n,c]
            assert (n-c)*(h[n+1,c]-h[n,c])+Q(n-c,n)*(h.get((n,c+1),0)+z(n)-h[n,c])<=q*h[n,c]
    assert q*h[D,1]>=z(D) and min(h.values())==Q('.765931373')
    return d,q,h


def excursion():
    d,q,h=degree_weight();rho=d+q
    # W applied to the whole exit tree; at the virgin impose the terminal
    # reward (rather than charging W there). For all states the fluid part
    # obeys (s-R)W >= h + (s-rho)*W. This supplies a strictly smaller scale.
    m=min(h.values());gamma=S-rho
    c=1/(m*(1+gamma/(Q(1,2)+d)))
    assert all(c/(Q(1,i)+d)>=1 for i in range(2,21))
    @cache
    def W(t):
        i,L,ch=t;n=1+L+len(ch);k=1+len(ch)
        return (h[min(n,20),min(k,20)]+sum(Q(1,y[0])*W(y) for y in ch))/(Q(1,i)+d)
    old=(2,2,((3,0,((2,1,()),)),(4,3,())))
    # All copy modes act on the same W and the same old subtree; gain gives
    # W(Y+)<=W(Y), so the a*W comparison preserves that pairing.
    assert W(old)>=W(phy.plus(old))
    assert c<1/m
    print('all-copy short-action envelope scale',float(c),'old 1/hmin scale',float(1/m),
          'heterogeneous example envelope',float(c*W(old)),flush=True)
    cap=Q(280);Wmax=1/(Q(1,20)+d)
    omitted=c*Wmax/sum((gamma*cap)**j/factorial(j) for j in range(160))
    print('one-virgin-row omitted at 280 for max terminal 1',float(omitted),flush=True)


def cut(pair):
    @cache
    def P(t,h):return pair.local(*shape(t,h))[2][1]
    @cache
    def A(t):return P(t,1)-sum(P(y,1) for y in t[2])
    J=phy.newborn(2);v=A(J);assert v>0

    def rooted(t):
        i,L,ch=t;d=1+L+len(ch);p=Q(1,d);a=Q(1,i);b=Q(L,d)
        intrinsic=sum(r*(A(u)-A(t)) for u,r in phy.intrinsic_moves(t))
        copy=sum((Q(1,y[0])-p)*A(y)+p*A(phy.plus(y)) for y in ch)
        direct=b*A(phy.newborn(min(d,20))) if d>2 else 0
        drift=intrinsic-a*A(t)+copy+direct-S*A(t)
        # Intrinsic drift of P at the root minus that at detached children,
        # with actual sending pressure into each previous child restored.
        D=lambda x:sum(r*(P(u,1)-P(x,1)) for u,r in phy.intrinsic_moves(x))
        inward=p*sum(P(phy.plus(y),1)-P(y,1) for y in ch)
        birth=b*P(phy.newborn(min(d,20)),1)
        assert intrinsic==D(t)-sum(D(y) for y in ch)-inward-birth
        return A(t)/v,drift/v,(intrinsic/v,copy/v,direct/v,-a*A(t)/v)

    def old(i,L,*ch):return (i,L,tuple(sorted(ch)))
    for i in range(2,21):
        val,r,terms=rooted(phy.newborn(i))
        if i in (2,20):print('virgin mark',i,'normalized weight',float(val),'residual',float(r),flush=True)
    xs=[old(2,L,old(2,4,old(3,0,old(2,1)),old(4,3)),old(3,2,old(2,1))) for L in (1,4,12,24,70)]
    xs += [old(i,L) for i in (2,20) for L in (2,3,4,5,12,24,70)]
    out=[(rooted(t)[1],t,rooted(t)) for t in xs]
    for key in (min,max):
        r,t,(w,r,terms)=key(out)
        print('extremum',t,'weight',float(w),'residual exact',r,'float',float(r),'terms',tuple(map(float,terms)),flush=True)
    print('positive test rows',sum(r>0 for r,_,_ in out),'of',len(out),'positive with grandchildren',sum(rooted(t)[1]>0 for t in xs[:5]),flush=True)
    assert max(out)[0]>0
    # A free-tree bound on the closed potential also has no root-mark term.
    # Check the *same* geometry and physical inward/copy arrows with two
    # admissible exterior marks to isolate that missing transport term.
    t=xs[0];u=(20,*t[1:]);w,d,_=rooted(t);w2,d2,_=rooted(u)
    assert w==w2 and d2-d==Q(9,20)*w
    # This tests the direct transport whose terminal profile is the virgin
    # trace of A; it is not a claim about arbitrary terminal corrections.
    val,dr,_=rooted(t)
    assert dr>0 and val>0
    print('same geometry mark-transport residual increment',float(d2-d),flush=True)


if __name__=='__main__':
    p=ArgumentParser();p.add_argument('lo',type=Path);p.add_argument('hi',type=Path);a=p.parse_args()
    cut(Pair(Cert(a.lo),Cert(a.hi)))
    excursion()
    source_budget()
