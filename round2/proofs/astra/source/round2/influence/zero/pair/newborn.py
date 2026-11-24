#!/usr/bin/env python3
"""Certified true newborn interval with a finite defect and the global TV cap."""
from fractions import Fraction as F
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from solver import occupations, enclosure, lower, upper


def down(x):
    return '-'+upper(-x) if x<0 else lower(x)


def main(N,S):
    # Production cap 48 certifies q uniformly for S >= .565658.
    q=F(8600511260,10**10)
    if S<F('0.565658'): raise ValueError('discount below proved threshold')
    dom,V,T,f=occupations(N,S)
    cand,_=enclosure(T,f,[F(1)]*len(T))
    # Zero extension of the killed H candidate. No sign assumption on H.
    norm=max(map(abs,cand))
    from solver import size, intrinsic, plus, deg, J
    idx={z:i for i,z in enumerate(dom)}
    # Bound |V_killed g - V g| by first-exit mass times sup|g|/s.
    # For g=b+C h, b+Cc<=1 and |h|<=nu give |g|<=max(1,nu).
    esc={}
    for z in V:
        div=S+size(z)-F(1,deg(z))
        esc[z]=sum(r*esc.get(y,1) for y,r in intrinsic(z).items())/div
    defect=max((abs(x-f[i]-sum(t*cand[j] for j,t in T[i].items()))+
                (esc[z]+esc[plus(z)])*max(1,norm)/S for i,(x,z) in enumerate(zip(cand,dom))))
    # For S outside the derivative mask, use ||ΔV b|| <= q + ||ΔV p||,
    # or directly 1/s (positive b); plus q||candidate||. The direct bound is
    # sufficient for the pilot and is shape-independent.
    outside=1/S+q*norm
    err=max(outside,defect)/(1-q)
    reward,occ=V[J]
    approx=reward+sum(t*cand[idx[y]] for y,t in occ.items())
    # Full C occupancy >= killed occupancy, <= 1/s. Use inflight mass to
    # pay future rewards on the union of exact and extended scores.
    radius=(sum(occ.values())+esc[J]/S)*err + esc[J]*max(1,norm)/S
    print('finite defect upper',upper(defect),'global defect upper',upper(max(outside,defect)))
    print('newborn candidate',lower(approx),upper(approx),'certified radius upper',upper(radius))
    print('true signed newborn', down(approx-radius),upper(approx+radius))


if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 8,F(sys.argv[2]) if len(sys.argv)>2 else F('0.565728'))
