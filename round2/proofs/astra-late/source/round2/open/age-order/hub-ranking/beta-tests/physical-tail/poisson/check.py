"""Apply the Poisson-stock serial floor to the pinned physical ledgers."""
from fractions import Fraction as Q
from argparse import ArgumentParser
from functools import cache
from hashlib import sha1
from math import factorial
from pathlib import Path
from runpy import run_path

EXT = Path(__file__).resolve().parents[6] / 'open/age-order/global-transport/annealed/extension/third-physical'
M = Q(36787, 100000)


def pin(name, expected):
    data = (EXT / name).read_bytes()
    assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest() == expected, name


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=20)
    parser.add_argument('--stop', type=int, default=35)
    parser.add_argument('--finite-only', action='store_true')
    args = parser.parse_args()
    assert 20 <= args.start <= args.stop <= 35
    for name, expected in (
        ('check.py', 'e39f8bfdf66f19e34a59247fa3649ccaf884358b'),
        ('band/check.py', 'e2aa75a98bf62981a91e33c3b3664e2f55d248c6'),
        ('scaled/check.py', '85e9ce571c1b6a842e6d3aa8faa71dbda147b2df'),
    ):
        pin(name, expected)
    assert sum(Q((-1)**j, factorial(j)) for j in range(10)) > M
    Y = Q(9,2)
    if not args.finite_only:
        symbolic(Y)

    band = run_path(str(EXT / 'band/check.py'))
    cuts = tuple(map(Q, ('1/8','77/600','2/15','1/7','1/6','1/4')))
    for b in range(args.start, args.stop):
        seam = Y/(b+1)
        for left, right in zip(cuts, cuts[1:]):
            if left >= seam:
                break
            right = min(right, seam)
            ok, kind = band['ledger'](b, M, left, right)
            assert ok, (b, left, right, kind)
        print(f'PASS b={b} on [1/8,{seam}] to symbolic seam', flush=True)
    assert 36*cuts[0] == Y


def symbolic(Y):
    scaled = run_path(str(EXT / 'scaled/check.py'))
    chains, plain = scaled['chains'], scaled['cap']
    assert chains(5, M) == (True, 735)
    print('PASS locked ledger at (b+1)x >= 5 (hence quarter for b>=20)', flush=True)

    S = scaled['S']
    # Keep the inherited high-degree cut with proxy Y0=4, and use the
    # Taylor cap at the true cutoff Y for smaller degrees.
    assert 2*3**32 > 4*2**33 and Q(33,32)*Q(3,2) < 2

    @cache
    def tighter(k, y0):
        if k >= 32:
            return plain(k, y0)
        w = Y**(k-1)
        lo, hi = S, 5*S
        assert hi**k*w.denominator > S**k*w.numerator
        while hi-lo > 1:
            z = (lo+hi)//2
            if z**k*w.denominator <= S**k*w.numerator:
                lo = z
            else:
                hi = z
        s = 1-Q(1,k)
        c = Q(20,21)-(1+s)/Y
        assert 0 < c and (1+s)*c < Y
        return (1-s*c/Y)/(Q(lo,S)*scaled['gamma_lower'](k))

    chains.__globals__['cap'] = tighter
    u = tighter(2, 4)
    g = M-u
    assert Q(32,33)*(g-Q(2,128))+(g-u)/2 >= 0
    assert chains(4, M) == (True, 735)
    print('PASS Taylor-tightened symbolic (b+1)x >= 9/2 for b>=20', flush=True)

if __name__ == '__main__':
    main()
