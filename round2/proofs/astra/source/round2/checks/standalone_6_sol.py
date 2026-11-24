#!/usr/bin/env python3
"""Independent exact endpoint controls for the Claim 6 cold audit."""

from fractions import Fraction as Q


def H(c):
    return (Q(2), Q(1), Q(17,20), Q(39,50), Q(7,10))[min(c,4)]


def crown(deg, core):
    """Founder with core pendant children, each bearing a single leaf."""
    stock = deg-core
    d = Q(3,20)
    a = Q(1,2)
    child_w = 1/(a+d)
    child_factor = a/(a+d)
    founder_z = H(core)+core*child_factor
    direct = Q(stock,deg)/(Q(1,min(deg,20))+d)
    immigration = direct+core*a*child_w
    assert immigration == direct+founder_z-H(core)
    assert immigration <= 1/d+2*(core+1)
    assert core*child_w <= Q(26,3)*(core+1)
    return founder_z, immigration, core*child_w


def variation(deg, c, *, capped):
    """Exact worst root gain for opposite successive label scores."""
    l = deg-c
    rate = l+Q(c-1,2)  # all core children virgin; exterior slot silent
    b, bp = Q(l,deg),Q(l+1,deg+1)
    label_changes = not capped or deg<20
    return rate*(b+bp if label_changes else bp-b)


def reconstruct_replay():
    # Fixed one-core planted stars at arbitrarily high leaf stock, and both
    # sparse and dense founder extremes; no allocation proportional to D.
    largest = Q(0)
    count = 0
    for deg in (*range(2,65),10**6,10**12,10**60):
        for c in sorted({1,min(2,deg),max(1,deg//2),deg-1,deg}):
            v = variation(deg,c,capped=True)
            assert 0<=v<=40
            if deg>=20:
                assert v== (deg-c+Q(c-1,2))*Q(c,deg*(deg+1)) < 1
            largest=max(largest,v)
            crown(deg,c)
            count+=1
    assert variation(10**60,1,capped=False)>10**60
    assert variation(10**60,1,capped=True)<Q(1,10**59)
    # Same core skeleton with 10^60 more leaves: old-fringe contribution
    # and Z really are unchanged, even though initial event rates diverge.
    a=crown(40,7)
    b=crown(10**60+40,7)
    assert a[0]==b[0] and a[2]==b[2] and a[1]<b[1]
    print(f'PASS: {count} analytic star/crown endpoint rows through D=10^60; max capped gain {largest}; uncapped plant rejected')


def finite_rate_replay():
    # Exact exponents for a conservative grid of admissible spectral margins.
    for mu in (Q(113,200), Q(1133,2000)):
        for eta in (Q(1,10**9),Q(1,100), (mu-Q(51,100))/3):
            delta=2*mu-1
            a=delta/8
            beta=min(delta-(4-2*mu)*a,2*eta*a)
            gamma=min(delta-3*a,2*eta*a,beta*(1-a))
            r=gamma/4
            assert a<1 and beta>0 and gamma>0 and 0<r<Q(1,2)
            assert gamma-3*r>0
            # Mutant: a too-fine mesh uses more than the full second-moment
            # exponent once its cardinality and error threshold both count.
            assert gamma-3*(gamma/2)<0
    print('PASS: 6 exact finite-rate budget rows (including eta=10^-9); oversampled mesh plant rejected')


if __name__=='__main__':
    reconstruct_replay()
    finite_rate_replay()
