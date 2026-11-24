#!/usr/bin/env python3
"""Correct the pinned 835c upper H3 slope rounding in memory and replay its mask."""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

GOLD = Path(__file__).resolve().parents[2]
SOURCE = GOLD/'higher/cohort.py'
PIN = 'b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c'
OLD = b'vg-(-GRID*affine[3][i]//1)'
NEW = b'vg+(-GRID*affine[3][i]//1)'


def main():
    assert __debug__
    raw = SOURCE.read_bytes()
    assert sha256(raw).hexdigest() == PIN and raw.count(OLD) == 1
    data = raw.replace(OLD, NEW)
    # For outward row q=A+Bρ, the affine branch of [-q]+ uses
    # -B per unit of unknown exterior receive, rounded toward -infinity.
    for b in (Q(-2, 3), Q(), Q(2, 3)):
        assert (-2**20*b)//1 <= -2**20*b < (-2**20*b)//1+1
    m = module_from_spec(spec_from_file_location('gold2_corrected_cohort', SOURCE))
    exec(compile(data, str(SOURCE), 'exec'), m.__dict__)
    a = SimpleNamespace(side='upper', min_mother=80, mother=480, mother_stock=12,
                        middle=2, middle_stock=2, young=12, young_stock=3,
                        sibling=2, upstream=2, outside_threshold=True, exclude_strict=False)
    with m.host.open_host('upper') as h:
        print('corrected cohort source', PIN, 'executed', sha256(data).hexdigest(), flush=True)
        print('upper original Gold2', h.docs['upper_sha256'], a, flush=True)
        r = m.Rows(h)
        m.registration()
        m.physical_checks(r)
        solver = m.Cohort(r, a)
        assert solver.mask((79, 8), (2, 2), (12, 3)) is None
        assert solver.mask((80, 8), (2, 2), (12, 3)) is not None
        floor = solver.sweep(('hmy',), 'physical')
        assert floor > 0
    print('corrected upper mother-high ticks/core', floor*2**40)
    print('PASS corrected original Gold2 upper H3')


if __name__ == '__main__':
    main()
