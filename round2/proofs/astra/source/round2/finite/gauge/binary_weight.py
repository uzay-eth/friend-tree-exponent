"""An exact two-mark residual supersolution; all entries are rational."""
from fractions import Fraction as F
from math import ceil

A, BETA, DAMP = F(1,2), F(1,16), F(4,25)
UNIT = 10**9

def mark(degree):
    return A if degree < 16 else BETA

def solve(cutoff,q,verbose=True):
    z=lambda D:mark(D)/(mark(D)+DAMP)
    tail=F(ceil(UNIT*z(cutoff)/q),UNIT)
    h={(cutoff,c):tail for c in range(1,cutoff+1)}
    for D in range(cutoff-1,1,-1):
        for c in range(D,0,-1):
            L=D-c
            b=F(L,D)
            bnext=h.get((D,c+1),F(0))+z(D)
            vnext=h[D+1,c]
            candidates=[((L+e)*vnext+b*bnext)/(q+L+e+b)
                        for e in (F(0),F(c+1,2)+DAMP)]
            h[D,c]=F(ceil(UNIT*max(candidates)),UNIT)
    original=h[2,1]
    assert 0 < original <= 1
    h[2,1]=F(1)
    for D in range(2,cutoff):
        for c in range(1,D+1):
            L=D-c
            b=F(L,D)
            for e in (F(0),F(c+1,2)+DAMP):
                drift=(L+e)*(h[D+1,c]-h[D,c])+b*(h.get((D,c+1),0)-h[D,c]+z(D))
                assert drift <= q*h[D,c]
    assert z(cutoff)<=q*tail
    h0=max(F(2),max(h[D,1]+z(D) for D in range(2,cutoff+1)))
    hmin,hmax=min(h.values()),max(h.values())
    assert hmin>0
    ratio=max(F(1),max(h[D+1,c]/h[D,c] for D in range(2,cutoff) for c in range(1,D+1)))
    delta=F(113,200)-q-DAMP
    assert delta*q*hmax/(BETA+DAMP)<=hmin
    if verbose:
        print(f'Dmax={cutoff} checked_rows={sum(2*D for D in range(2,cutoff))} '
              f'h21_before={original} min={hmin} max={hmax} hcentre={h0} '
              f'one_leaf={ratio} d={DAMP} q={q} rho={q+DAMP} delta={delta}')
    return h,ratio

if __name__=='__main__':
    solve(16,F(2,5))
    solve(64,F(49,125))
