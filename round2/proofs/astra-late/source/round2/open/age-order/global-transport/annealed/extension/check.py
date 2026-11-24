"""Exact margins for paired edge debits and the old-jump occupation charge."""
from fractions import Fraction as Q
from math import factorial


def attach(a,u):
    return a[:u]+(a[u]+(len(a),),)+a[u+1:]+((u,),)


def profile(a):
    d=tuple(map(len,a))
    r=[sum(Q(1,d[v]) for v in row) for row in a]
    w=[q/(k+1) for q,k in zip(r,d)]
    s=[q/k for q,k in zip(w,d)]
    h=[sum(d[v]>1 for v in row) for row in a]
    return d,r,w,s,h


def observables(a):
    d,r,w,s,h=profile(a)
    A=sum(w)
    F=sum(w[u] for u,k in enumerate(d) if k==1)*2
    edge_hi=Q(0); edge_lo=Q(0)
    for u,row in enumerate(a):
        for v in row:
            if s[u]<=s[v]: continue
            if d[u]==1: edge_hi+=s[u]-s[v]
            else: edge_lo+=s[u]-s[v]
    assert edge_hi<=A/3 and edge_lo+edge_hi<=A
    assert edge_hi<=sum((1-Q(2,d[a[u][0]]+1))*w[u] for u,k in enumerate(d) if k==1)
    assert A <= 3*sum(k>=2 for k in d)/Q(2)
    return A,F,edge_hi,edge_lo


def occupation(a):
    d,r,w,s,h=profile(a)
    n=len(a); A,F,*_=observables(a)
    drift=Q(0)
    for u,prob in enumerate(r):
        newF=observables(attach(a,u))[1]
        L=d[u]-h[u]
        direct=(Q(L+1,d[u]+1)-Q(L,d[u])) if d[u]>1 else Q(1,2)-Q(1,d[a[u][0]])
        assert newF-F==direct
        drift+=prob*direct
    assert drift>=0
    return A,F,drift


def coef(b):
    t=Q(5,6)
    l=sum(((-1)**(j+1)*Q(5,11)**j/j for j in range(1,12)),Q())
    e=sum(((-t)**j/factorial(j) for j in range(12)),Q())
    return Q(5*(b-1),16*b+5)*(1-l/t)*e*(1-t*t/(2*(b-t)))


def budget(Y,rho,kappa,b=None):
    # kappa is the coefficient for the smaller-beta/large-beta comparison
    # after factoring out Y^-1/6. Round the sixth root outward rationally.
    den=10**6
    root,hi=den,4*den
    while hi-root>1:
        mid=(root+hi)//2
        if mid**6<=den**6*Y: root=mid
        else: hi=mid
    theta=kappa*den/root
    assert theta<1
    cost=1+rho+(1-rho)*theta
    assert coef(Y if b is None else b)**2*Q(157,50)*Y > cost**2
    return cost


def fixed_threshold(x,bmin):
    # Keep both the core-tip density factor and the high-tip sqrt(1-x).
    assert Q(3,2)**6<=bmin*x
    theta=Q(4,5)*Q(2,3)
    cost=Q(4,3)+Q(2,3)*theta
    assert coef(bmin)**2*Q(157,50)*bmin*x > cost**2*(1-x)


if __name__=='__main__':
    seeds=[((1,),(0,2),(1,)),tuple(tuple(v for v in (u-1,u+1) if 0<=v<12) for u in range(12))]
    star=seeds[0]
    for _ in range(15): star=attach(star,1)
    seeds += [star,attach(attach(star,0),2)]
    for t in seeds:
        occupation(t)
    a=((1,2,3),(0,4),(0,5),(0,6),(1,),(2,),(3,))
    b=((1,2,3),(0,4),(0,),(0,),(1,5),(4,6),(5,))
    assert sorted(map(len,a))==sorted(map(len,b))
    assert all(any(v<u for v in t[u]) for t in (a,b) for u in range(1,len(t)))
    Aa,_,Xa,_=observables(a); Ab,_,Xb,_=observables(b)
    assert (Xa/Aa,Xb/Ab)==(Q(2,59),0)
    # Integral lower for Gamma(1/3): odd alternating Taylor polynomial
    # for exp(-t) on [0,1]; use only this interval.
    gamma_third=sum(((-1)**j*Q(3,3*j+1)/factorial(j) for j in range(6)),Q())
    assert gamma_third > Q(12,5)
    # sqrt(pi)<9/5; Jensen ratio distortion <= (b+1)/(b+1/2).
    assert Q(22,7) < Q(9,5)**2
    assert Q(9,5)/gamma_third*Q(282,281)<Q(3,4)
    assert Q(9,5)/gamma_third*Q(242,241)<Q(3,4)
    assert Q(9,5)/gamma_third*Q(42,41)<Q(4,5)
    budget(141,Q(1,3),Q(3,4))
    budget(140,Q(19,60),Q(3,4))
    budget(135,Q(1,4),Q(3,4))
    budget(120,Q(1,8),Q(3,4))
    fixed_threshold(Q(10,11),20)
    assert coef(141)>Q(73,1000)
    print('PASS true old-jump stock flux identity and paired edge caps on 4 physical fixtures')
    print('PASS degree-marginal obstruction: equal degree sequences, X/A = 2/59 versus 0')
    print('PASS paired cutoff 141; annealed cutoff tests (140,19/60), (135,1/4), (120,1/8)')
    print('PASS fixed threshold x>=10/11 for every n+1>=20')
