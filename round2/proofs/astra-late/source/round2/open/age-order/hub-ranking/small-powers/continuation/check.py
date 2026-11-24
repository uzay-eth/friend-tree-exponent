"""Exact discounted stock boundary and physical middle continuation."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
from runpy import run_path


FM = run_path(str(Path(__file__).resolve().parents[2] / 'fractional-middle/check.py'))
S = FM['S']


def stock_boundary(lo, hi, d, cutoff):
    """Value of reaching cutoff before zero, then cashing the harmonic moment."""
    assert 0 <= lo <= hi <= 1 and 2 <= d <= cutoff
    pol = [Q(1)]
    for j in range(1, cutoff):
        pol.append(pol[-1] * (1 + lo/j))
    e_inv = sum(Q((-1)**j, factorial(j)) for j in range(18))
    terminal = e_inv * sum(pol[cutoff-1-k]/factorial(k) for k in range(cutoff))
    ratios = [Q(0)]
    for j in range(1, cutoff):
        death = min(Q(j,d), 1)
        ratios.append(j / (j+death+hi-death*ratios[-1]))
        assert 0 < ratios[-1] <= 1
    row = [Q(0)]*cutoff + [terminal]
    for j in range(cutoff-1, 0, -1):
        row[j] = ratios[j]*row[j+1]
    for j in range(1, cutoff):
        death = min(Q(j,d), 1)
        assert j*row[j+1]+death*row[j-1] == (j+death+hi)*row[j]
    return [int(S*x) for x in row[:d+1]]


def middle(lo, hi, d, cutoff=None, sender=96):
    cutoff = d if cutoff is None else cutoff
    after = stock_boundary(lo, hi, d, cutoff)
    q,z = hi.numerator,hi.denominator
    for k in range(d-1, 1, -1):
        ret = 0
        for h in range(sender-1, 1, -1):
            g = h*(k*(h-1)+1)
            ret = z*(k*S+g*ret)//(z*(k+g)+k*h*q)
        row = [ret*after[1]//S]
        for j in range(1,k+1):
            row.append(j*z*(k*after[j+1]+row[-1])//(k*q+j*(k+1)*z))
        after = row
    return Q(after[1],S)


def tips(lo, hi, stop=32):
    """Discount and sibling stock jointly at the physical cutoff."""
    gm, up = FM['gamma_upper'], FM['up']
    rows = {}
    for k in range(stop, 1, -1):
        a = 1+Q(1,k)
        unshifted = min(a, gm(a,lo,hi))
        cap = min(unshifted, up((1-1/(hi+a+1))*gm(a,lo,hi)
                  +(gm(a-Q(1,2),lo,hi)+gm(1+Q(1,k+1),lo,hi))/(2*(lo+a+1))))
        for j in range(k-1):
            if k == stop:
                rows[k,j] = up(unshifted/(j+1+lo*k))
                continue
            gain = k*(j+1)
            child = rows[k+1,j+1]
            rows[k,j] = max(child,up((cap+gain*child+j*rows.get((k,j-1),0))
                                      /(1+gain+j+lo*k)))
        assert all(rows[k,j] <= rows[k,j-1] for j in range(1,k-1))
    return rows


def main():
    # Independent endpoint identity at zero: gambler's ruin with death/birth 1/d.
    d=8
    row=stock_boundary(Q(0),Q(0),d,d)
    assert abs(row[1]-row[d]*(1-Q(1,d))/(1-Q(1,d)**d)) < 1
    worst = Q(1)
    for j in range(128):
        lo,hi=Q(j,512),Q(j+1,512)
        d = 8 if j < 24 else 12 if j < 54 else 16 if j < 92 else 20
        old,new=FM['middle'](lo,hi),middle(lo,hi,d)
        worst = min(worst,new-old)
        if j in (0, 31, 63, 127):
            print(f'PASS [{lo},{hi}]: serial floor {float(old):.9f} -> {float(new):.9f} (D={d})')
    assert worst > Q(1,4000)
    print(f'PASS 128 boxes covering (0,1/4]: smallest improvement {worst}')

    for j in range(16):
        for x in range(16):
            k = j+x+2
            for p in (Q(0),Q(1,512),Q(1,8),Q(1)):
                z = j+p*k
                direct = k*(j+1)*(1+p)*z-j*(z+1)*(z+2+p)
                expanded = j*(j+1)*x+p*((j+2)**2*x+(j+1)*x*x+3*j+4)+p*p*k*(x+2)
                assert direct == expanded >= 0
    baseline = run_path(str(Path(__file__).resolve().parents[2] / 'mean-ledger/fractional/check.py'))
    for j in (0, 7, 15, 31, 63, 127):
        lo,hi=Q(j,512),Q(j+1,512)
        before=baseline['Bounds'](lo,hi,32)
        after=tips(lo,hi)
        assert all(v<=before.U(k,l) for (k,l),v in after.items())
        assert after[2,0]<before.U(2,0)
        print(f'PASS [{lo},{hi}]: tip cap {float(before.U(2,0)):.9f} -> {float(after[2,0]):.9f}')


if __name__ == '__main__':
    main()
