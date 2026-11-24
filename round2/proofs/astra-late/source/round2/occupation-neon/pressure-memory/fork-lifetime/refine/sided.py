#!/usr/bin/env python3
"""Reprice frozen, checked full-mask forecasts with the sided physical law."""
from fractions import Fraction as Q
from math import ceil, floor
from pathlib import Path
from replay import B, GROUPS, GRID, load

# Exact discounted ticks from `--caps 0 3 2 --mother 2 --tip32 --forecast`.
FORECAST = {
    '9/10': ((-29005083669,27156926218),(-32501793863,43654974690),
             (-19019618764,42981770416),(-8658913568,42422893768),(-3253555020,42122543191)),
    '1': ((-33164288162,26188101184),(-37545852063,43003642040),
          (-23231591562,42324809936),(-12139928146,41762153471),(-6424315612,41460093733)),
}


def main():
    s=load('birth_law', B.NEON/'pressure-memory'/'birth-law'/'replay.py')
    a=load('p4_kernel', B.NEON.parent/'checks'/'stationary-p4-sol.py')
    lp=s.build(24,True,a.class_kernel())
    print('sided physical LP',len(lp.rows),'rows',len(lp.cols),'columns',flush=True)
    for theta,table in FORECAST.items():
        def entry(d, side):
            i=max(i for i,g in enumerate(GROUPS) if g<=d)
            return Q(table[i][side],GRID)
        x,y=s.price(lp,lambda d,k:entry(d,0),lambda d,k:entry(d,1))
        l,h=map(Q, ('-.000001532406105','.000002265035818') if theta=='9/10' else
                         ('-.000001790570406','.000002224397548'))
        x,y=max(x,l),min(y,h)
        assert x<=y
        f=B.decimal15
        print(theta, f(floor(x*10**15)),f(ceil(y*10**15)),
              'width<=',f(ceil((y-x)*10**15)),flush=True)


if __name__=='__main__':
    main()
