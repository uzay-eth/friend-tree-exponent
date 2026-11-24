"""Source-specific losses at the degree-two direct terminal."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path


P = run_path(str(Path(__file__).with_name('check.py')))
F = P['FROZEN']
G,U,L = F['G'],F['U'],F['E']['L']
A,B = F['A'],F['B']


def main():
    atom,I = P['stock_cost'](16)
    serial,parallel = P['source_times'](16)
    root = Q(1,80_000_000)
    base = atom+8/P['SHARP']['exp_lo'](Q(24))
    loss2,loss3,edge3 = (base+4*I*root*t for t in (serial(2),serial(3),parallel(3)))
    for reciprocal in (Q(0),Q(1,31)):
        w = (1+reciprocal)/3
        s,t = w/2,Q(1,4)
        assert s<t
        old = w*G(3)+t*G(2)-(t-s)*max(0,U(2,0)-L(3))-A-max(B,s*U(3,1))
        cost = w*loss3+t*loss2+(t-s)*edge3
        assert cost<old
        print(f'PASS reciprocal {reciprocal}: frozen slack {float(old):.9g}, '
              f'source-specific charge {float(cost):.9g}')
    print(f'PASS degree-two direct terminal against degree >=31 for y <= {root*root}')


if __name__ == '__main__':
    main()
