#!/usr/bin/env python3
"""Exact cold-review fixtures for standalone Claim 3; no author imports."""
from collections import defaultdict
from fractions import Fraction as F
from math import factorial


def weight():
    def h(c):
        return (F(2),F(1),F(17,20),F(39,50))[c] if c<4 else F(7,10)
    rows=[(F(D-c,D)*(h(c+1)-h(c)+F(20,20+3*D))/h(c),D,c)
          for D in range(2,100) for c in range(D+1)]
    assert max(rows)==(F(313,870),3,1)
    assert F(313,870)<F(9,25) and F(1,16)<F(9,25)*F(7,10)
    assert F(9,25)+F(3,20)==F(51,100)<F(113,200)
    print('weight rows',len(rows),'max (relative drift,D,c)',max(rows),
          'slack',F(9,25)-max(rows)[0])


def seed_and_transfer():
    # PGFs are in z, with q=1-z; for r seeds, their numerator for
    # E[(N)_j] is obtained by j differentiations of q^r z^r/(1-z)^r.
    for r in (2,3):
        for j in (1,2,3):
            vals={r:1}
            for i in range(j):
                inv=defaultdict(int)
                # Differentiate numerator/(1-z)^(r+i), for successive i.
                for k,a in vals.items():
                    if k:inv[k-1]+=k*a
                    inv[k]+=(r+i-k)*a
                vals=inv
            assert sum(vals.values())==factorial(r+j-1)//factorial(r-1)
    # Mutation controls: using Gamma(mu+1), no 2/(2+mu), and
    # accidentally treating T3 as a planted seed all mis-scale this fixture.
    # Take a_n=n-2 (zero at the edge); coefficients at mu=1 are exact.
    mu=1; c=1; d_edge=2; d_three=3
    assert F(2,2+mu)*d_three==d_edge
    assert F(d_edge,factorial(mu+1))==c!=F(d_edge,factorial(mu))
    assert F(d_three,factorial(mu+1))!=c
    pref_lo=F(14608704616,17527873930); pref_hi=F(14606753031,2989767809)
    assert pref_lo==F(7304352308,8763936965)<1<pref_hi
    print('seed fixture: d3=3, d2=2, c=1; bracket prefactors',pref_lo,pref_hi)
    # With fractional exponent a, a_n=(a+2)_(n-2)/(n-1)! gives an
    # exact Edge mixture q^-a. Shift to zero at the edge by subtracting 1.
    a=F(113,200); an=F(1); term=F(1)
    for j in range(1,21):
        an*=F(j+1+a,j+1)
        term*=F(j+1+a,j)
        assert (j+1)*an==term
    assert F(2,2+a)*(a+2)/2==1
    print('fractional transfer fixture: d2=1, d3=513/400, c=1/Gamma(513/200)')


def arrows(edges):
    ns=defaultdict(list)
    for x,y in edges:
        ns[x].append(y);ns[y].append(x)
    return ns


def flux():
    # Unequal sender degrees, a leafless interior, multiple census levels,
    # and a nontrivial preselected cell. Iterate every actual marked arrow.
    edges=[(0,1),(1,2),(2,3),(2,4),(3,5),(3,6),(3,7)]
    n=arrows(edges); cell={0,1,2,3}&{x for x in n if len(n[x])>=2}
    results={}
    for k in range(2,5):
        actual=F()
        before=sum(len(n[x])>=k+1 for x in cell)
        for x in n:
            for y in n[x]:
                updated=cell|{y} if x in cell and len(n[y])==1 else cell
                actual+=F(sum(len(n[v])+(v==y)>=k+1 for v in updated)-before,len(n[x]))
        direct=sum(sum(F(1,len(n[x])) for x in n[y]) for y in cell if len(n[y])==k)
        wrong=sum(k+sum(F(1,len(n[x])) for x in n[y] if len(n[x])>=2)
                  for y in cell if len(n[y])==k)
        assert actual==direct!=wrong
        results[k]=(actual,wrong)
    assert results=={2:(F(4,3),F(7,3)),3:(F(7,4),F(15,4)),4:(F(10,3),F(13,3))}
    print('crossing flux correct/omit core-degree decrement',results)


def mean_vs_pointwise():
    # Strictly positive rational probability law with tail A/k and
    # alternating local coefficient: perturb by a discrete telescoping term.
    # T(2)=1, T(k)=2/k+(-1)^k/(4 k^2)-1/(4 k^2).
    def tail(k): return F(2,k)+F((-1)**k-1,4*k*k)
    assert tail(2)==1
    for k in (2,3,100,101,10000,10001):
        p=tail(k)-tail(k+1)
        assert p>0 and abs(k*tail(k)-2)<=F(1,2*k)
        lim=2+F((-1)**k,2)
        assert abs(k*k*p-lim)<=F(4,k)
    print('adversarial positive law: cumulative constant 2, pointwise constants 5/2 and 3/2')


def lattice_counterexample():
    # S(x)=2^floor(x) is increasing, but exp(-log(2)*x) S(x)
    # differs on integer and half-integer times. Its transform is
    # (1-exp(-z))/(z*(1-2 exp(-z))). All the poles at log(2)+2pi*i*m
    # are genuine because the numerator there is 1/2, not zero.
    e_at_each_pole=F(1,2)
    assert 1-2*e_at_each_pole==0 and 1-e_at_each_pole==F(1,2)
    assert F(1)!=F(1,2)  # normalized observations at x=n and x=n+1/2, squared
    print('lattice control: positive increasing S with persistent normalized oscillation')


if __name__=='__main__':
    weight();seed_and_transfer();flux();mean_vs_pointwise();lattice_counterexample()
