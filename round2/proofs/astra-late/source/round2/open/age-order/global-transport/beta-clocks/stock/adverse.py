"""Uniform finite-size upper prices for promotion before a sharp share threshold."""

from fractions import Fraction as Q
from functools import cache
from pathlib import Path
from runpy import run_path

POS = run_path(str(Path(__file__).with_name('check.py')))
BASE = POS['BASE']
ceil = POS['ceil']
S = 10**12


def high_constant(y, minimum=21, panels=128):
    """A/(k-1) bounds every locked price k>=32, including the large-degree tail."""
    assert 0<y<2 and minimum==21 and Q(9,8)**32 > y+minimum-1
    return Q(9,8)*sum(1/(y+max(0,w)) for w in BASE['wait_cells'](panels,minimum-1))/panels


def terminal_grid(left, right, minimum, stop, rate, times, mesh):
    """Event-time expectations of the pending and the promoted locked upper."""
    assert 0<left<=right<minimum and stop<minimum
    unit=10**100
    ratio=Q(mesh+1,mesh)
    xi=Q(1,2*mesh+1)
    loglo=2*(xi+xi**3/3+xi**5/5)
    loghi=loglo+2*xi**7/(7*(1-xi**2))
    step_lo=unit//ratio**rate
    step_hi=ceil(unit/ratio**rate)
    surv_hi=[unit]*times
    pl=ph=elo=ehi=unit
    rows=[]
    cut=left if left >= Q(1,2) else left*(1-Q(1,2*minimum))
    pcap=BASE['price_caps'](cut,minimum=minimum)
    mesh_width=10**12
    at_zero=[pcap(k,cut) for k in range(2,stop+1)]
    gammas=[BASE['gamma_lower'](k)/k for k in range(2,stop+1)]
    minimum_q=1-right/minimum
    prom_decay=left if left>=Q(3,2) else left*(1-Q(3,2*minimum))
    def exp_hi(z):
        # A genuine rational upper on exp(-z), via a lower on exp(z).
        squares=0
        while z>Q(1,4):
            z/=2
            squares+=1
        den=1+z+z*z/2+z**3/6+z**4/24
        u=ceil(mesh_width/den)
        for _ in range(squares):
            u=ceil(Q(u*u,mesh_width))
        return Q(u,mesh_width)
    def prices(v,w):
        a=exp_hi(prom_decay*v)
        b=exp_hi(left*(1-Q(1,minimum))*v)
        bound=[ceil(S*b*item) for item in at_zero]
        # Jensen for the residual Beta(1,m-r) expectation, with r=1+1/k.
        for k,gamma in zip(range(2,stop+1),gammas):
            base=(1+right*(1+w))/minimum_q
            root=Q(1+BASE['root_ticks'](base,k,mesh_width),mesh_width)
            bound.append(ceil(S*a*root/gamma))
        return bound
    i=0
    while Q(pl,unit)-1<40/left:
        v=Q(pl,unit)-1
        pl=pl*(mesh+1)//mesh
        ph=ceil(Q(ph*(mesh+1),mesh))
        w=Q(ph,unit)-1
        i+=1
        elo=elo*step_lo//unit
        ehi=ceil(Q(ehi*step_hi,unit))
        cdf=POS['poisson_cdf'](times,elo,ehi,i*rate*loglo,i*rate*loghi,unit)
        mass=[max(0,prev-l) for prev,(l,h) in zip(surv_hi,cdf)]
        surv_hi=[h for _,h in cdf]
        rows.append((mass,prices(v,w)))
    # The true beta tails decrease in the random threshold. Apply the last
    # cell's bound at its left endpoint to the entire unfinished time tail.
    rows.append((surv_hi,rows[-1][1]))
    levels=[[0]*(2*(stop-1)) for _ in range(times)]
    for mass,bound in rows:
        for t,weight in enumerate(mass):
            if weight:
                for k,b in enumerate(bound):
                    levels[t][k]+=weight*b
    return [[ceil(Q(v,unit)) for v in level] for level in levels], at_zero


def prices(left,right,minimum=21,stop=14,rate=16,times=65,mesh=64):
    """Upper on each virtual state at the initial common completion size."""
    assert stop+1<=rate<=minimum
    levels,at_zero=terminal_grid(left,right,minimum,stop,rate,times,mesh)
    nextval={(k,j):levels[-1][k-2] for k in range(2,stop+1) for j in range(k)}
    for t,level in reversed(list(enumerate(levels))):
        val={}
        prev=levels[t-1] if t else [ceil(S*x) for x in at_zero]
        for k in range(2,stop):
            up=level[stop-1+k-2]
            for j in range(k):
                gain=(j+1)*k
                hold=rate*k-gain-j-1
                assert hold>=0
                total=gain*nextval[k+1,j+1]+j*nextval[k,max(0,j-1)]+up+hold*nextval[k,j]
                val[k,j]=min(prev[k-2],ceil(Q(total,rate*k)))
        for j in range(stop):
            val[stop,j]=prev[stop-2]
        nextval=val
    return {key:Q(v,S) for key,v in nextval.items()}
