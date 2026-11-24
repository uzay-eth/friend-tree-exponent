"""Exact slack for the frozen edge ledger and compact retry time budgets."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parents[1] / 'hostile'
E = run_path(str(root / 'edge-dual/check.py'))
M = E['M']
H = M['D']['H']
G,U,P,A,B,F,D2 = (E[k] for k in ('G','U','P','A','B','F','D2'))
delta = Q(1,200000)


def check():
    high = min(E['row'](d,c,*P(d,(e,)*c))-c*cost
               for d in range(10,113) for c in range(2,d+1)
               for cost,e in ((0,0),(0,2),(F,4),(F,2),(D2,2)))
    bracket=M['g']*Q(113,114)+G(2)/2-M['p']/2
    tail=min(bracket/2+(A-F)*Q(113,2)-2*A-B,
             M['g']*Q(113,228)-2*F-B, M['g']*Q(113,228)-2*D2-B)
    terminal=Q(9,11)*M['g']+Q(9,20)*G(2)-Q(1,22)-A-D2
    assert min(high,tail,terminal)>delta
    print('PASS every high endpoint retains slack > 1/200000')

    for d in range(2,10):
        worst=Q(1)
        bill=Q('.021') if d==2 else Q('.0037')
        for c in range(1,d+1):
            if d==c==2: continue
            for ds in combinations_with_replacement((0,*range(2,10)),c):
                for r in (0,10):
                    s,C,u=P(d,tuple(e or r for e in ds))
                    cost=sum(max(0,E['upper_out'](d,e,s,u)-B,bill if e==2 else 0) for e in ds)
                    receipts=F*ds.count(0) if (d,c) in ((3,2),(3,3),(4,4)) else 0
                    worst=min(worst,C+A*(c-2)-B-cost+receipts)
        assert worst>delta
        print(f'PASS low endpoint degree {d}: slack >= {float(worst):.12f}')
    assert G(3)/3>delta/4
    # Worst suppressed-chain margins in the frozen exact replays.
    assert B+2*D2-Q(70627717243,2400000000000)>delta
    assert B-Q('.016905')>delta and B-Q('.010694')>delta
    for d in (2,3):
        lo,_=H['score'](H['spider'](d,1))
        assert lo>Q(1,100)
    assert M['g']*Q(4,5)>Q(1,100)

    # S tail: E exp(S/2)=2^m; S>32m costs (2e^-16)^m.
    # e>2 suffices; here m>=3.
    assert 3*(16-1)>40
    # Time bound: 700 stages, mgf at 1/64 at most 2 per stage.
    assert 2**16//64-700>40
    # Extra age 128: e^(-(19/20)*128) < 2^-121. At most 20 families.
    assert Q(19,20)*128>121
    assert 2*20*Q(2**32,2**121)+Q(8,2**32)+Q(2,2**40)<Q(1,2**27)
    assert 131400>2*(2**16+128)+5
    assert Q(6,2**27)<delta
    assert G(2)>Q(1,2**27)
    assert Q(4*2**32,2**64)+Q(8,2**32)+Q(1,2**40)<Q(1,2**27)
    assert Q(3,4)-1+Q(108,305)==Q(127,1220)>Q(1,2**27)
    print('PASS star, time, true-amplitude and total signed-loss budgets')


if __name__=='__main__':
    check()
