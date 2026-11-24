"""Exact uniform normalized seed/star numerator on 0 < p <= 1.

The seed has zero paired debit and total serial weight one. For a star,
use the pinned old-leaf-adjusted profile at (d, c, r)=(d, 0, 0).
Above its finite degree cutoff, retain just the pinned analytic centre-side
credit minus the upper on its outgoing debit. The first exponent box uses
zero only as an endpoint of an enclosure for strictly positive exponents.
"""
from fractions import Fraction as Q
from hashlib import sha1
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
FLOOR = Q(9, 5000)


def authenticate(path, expected):
    data = path.read_bytes()
    blob = sha1(f'blob {len(data)}\0'.encode() + data).hexdigest()
    assert blob == expected, (path, blob)


def main():
    lower_path = HERE.parents[1] / 'fractional-middle/check.py'
    authenticate(lower_path, 'd8111c15376ad27f2fe0a583cabf8192f925ec88')
    authenticate(HERE / 'check.py', '76bbdfc337be47cf8a0a91138ee5d191264017ab')
    fractional = run_path(str(HERE / 'check.py'))
    Bounds, lower = fractional['Bounds'], fractional['FM']
    print('PASS both frozen checker blobs', flush=True)
    overall = {}

    for start, stop in ((0, 36), (36, 128), (128, 473), (473, 512)):
        minima = {}

        def book(kind, value, j, degree):
            assert value >= FLOOR, (kind, j, degree, value)
            row = value, j, degree
            minima[kind] = min(minima.get(kind, row), row)
            overall[kind] = min(overall.get(kind, row), row)

        for j in range(start, stop):
            lo, hi = Q(j, 512), Q(j + 1, 512)
            if j < 473:
                tip_cut, degree_cut = (48, 96) if j < 36 else (32, 64)
                box = Bounds(lo, hi, tip_cut)
                if j < 128:
                    box.boost()
                book('seed', box.g, j, 1)
                for d in range(2, degree_cut):
                    book('finite star', box.profile(d, 0, 0)[1], j, d)
                tail = Q(degree_cut, degree_cut+1)*box.g - (1+Q(1, degree_cut+1))/(degree_cut+1)
                book('high star line', tail, j, degree_cut)
                Bounds.profile.cache_clear()
            else:
                # Exact specialization of fractional-middle's local profile:
                # s=1/(d+1)>t=1/(2d), so there is no old-leaf-side debit.
                degree_cut = 32
                M = lower['middle'](lo, hi)
                U, L = lower['side_bounds'](lo, hi, M)
                g = M-U(2)
                book('seed', g, j, 1)
                for d in range(2, degree_cut):
                    w, s, t = Q(d,d+1), Q(1,d+1), Q(1,2*d)
                    C = (w+d*t)*g-d*(s-t)*max(0,U(d+1)-L(2))
                    book('finite star', C, j, d)
                u = 1/(lo*(degree_cut+1))
                assert g > u
                book('high star line', Q(degree_cut, degree_cut+1)*(g-u), j, degree_cut)
            lower['gamma_upper'].cache_clear()

        print(f'PASS cells {start}..{stop-1}:', flush=True)
        for kind, (value, j, degree) in minima.items():
            print(f'  {kind}: {value} at j={j}, d={degree}', flush=True)

    print('Global minima (j denotes [j/512,(j+1)/512]):')
    for kind, (value, j, degree) in overall.items():
        print(f'  {kind}: {value} at j={j}, d={degree}')
    print('PASS all 512 boxes: normalized seed and star numerator >= 9/5000')


if __name__ == '__main__':
    main()
