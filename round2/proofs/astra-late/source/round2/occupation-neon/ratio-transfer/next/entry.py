#!/usr/bin/env python3
"""Admit all live-arm histories until the disjoint bare-P4 continuation opens."""
from functools import cache
from pathlib import Path
from collections import defaultdict
import sys

from spine import Pair, Cert, Q, T, LO, HI, BASE, DEN, GRID, coefficients, super_poisson

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]/'check'))
from crowns import n, degree, moves


def floor(x):
    return Q(x*GRID//1,GRID)


def birth_price():
    # Four consecutive tip promotions open the face with a virgin arm.
    # A leaf gain before the fourth promotion would prevent a bare tip.
    out=Q(1)
    for size in range(2,6):
        out*=Q(1,2)/(HI+size)
    return out


def continuation(U, threshold):
    terminal = {}
    @cache
    def V(z):
        d,k=degree(z),len(z[1])
        if n(z)>=threshold:
            out=max(Q(0),-U.get((d,k),0)/T)
            terminal[z]=out
            return out
        jump=moves(z,Q(1,2))
        rate=sum(r for _,r in jump)
        assert rate==n(z)-Q(1,d)+Q(1,2)
        # Killing all four bare vertices costs 3+1/d+p. The physical
        # contribution of the same bottom arm occurs in both hazards.
        den=HI+rate+3+Q(1,d)+Q(1,2)
        assert den==HI+n(z)+4
        v=floor(sum(r*V(y) for y,r in jump)/den)
        assert v>=0
        return v
    value=V((1,()))
    # Independently check the backwards inequality by a forward exact
    # discounted flux to the same geometric boundary; only roundoff may
    # separate these answers. No optimizer or per-side charges here.
    mass={(1,()):Q(1)}
    for size in range(2,threshold):
        nxt=defaultdict(Q)
        for z,x in mass.items():
            assert n(z)==size
            for y,a in moves(z,Q(1,2)):
                nxt[y]+=x*a/(HI+size+4)
        mass=nxt
    direct=sum(x*terminal[z] for z,x in mass.items())
    assert 0<=direct-value<=(threshold-2)/Q(GRID)
    return value,terminal,V.cache_info().currsize


def main():
    c=Pair(*(Cert(Path(a)) for a in sys.argv[1:]))
    U,_,_=super_poisson(coefficients(c),12)
    # Four selected tip-promotion events create four bare spine cores
    # followed by one virgin. Adjacency confirms rate 1/2; exterior <=1/2.
    born=birth_price()
    for threshold in (9,10):
        arm,ends,num=continuation(U,threshold)
        active=[z for z,v in ends.items() if v]
        assert active and max(n(z) for z in active)==threshold
        whole=LO*born*arm
        gain=whole/DEN
        print('arm threshold',threshold,'contexts',num,'positive exits',len(active),'/',len(ends),
              'largest active descendant degree',max(max_degree(z) for z in active),flush=True)
        print('  arm value',float(arm),'pi delta+ lower',float(whole),
              'gain vs arms baseline',float(gain),'outward endpoint ticks at 1e-18',
              -(-(BASE-gain)*10**18//1),flush=True)
        if threshold==9:
            assert gain>Q('5.7011e-11')


def max_degree(z):
    return max([degree(z),*(max_degree(y) for y in z[1])])


if __name__=='__main__':main()
