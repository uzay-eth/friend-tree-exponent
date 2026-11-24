#!/usr/bin/env python3
"""Extract round-tripped binary64 endpoints and exact widths from finished logs."""
from decimal import Decimal as D, localcontext
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
REF = D('0.00012210130022050069698025254183448851108551025390625')
CLIP = D('0.00010392855440299886282673469395376741886138916015625')
BASE = D('0.00013360681999630142')


def endpoint(case, mode):
    text = (HERE / f'{case}m{mode}.log').read_text()
    lam, = re.findall(r'^FLOAT.* lam=(\S+)$', text, re.M)
    wall, user, system, rss = re.findall(r'^TIME wall=(\S+) user=(\S+) sys=(\S+) peak_KiB=(\d+) exit=0$', text, re.M)[0]
    sweeps = re.findall(r'^SWEEP.* it=(\d+) lam=(\S+) rho=', text, re.M)
    assert [int(i) for i,_ in sweeps] == list(range(1,19))
    shift = abs(D.from_float(float(sweeps[-1][1]))-D.from_float(float(sweeps[-2][1])))
    return lam, wall, user, system, rss, shift


with localcontext() as ctx:
    ctx.prec = 80
    for case in sys.argv[1:]:
        lo, hi = (endpoint(case, m) for m in (2, 1))
        width = D.from_float(float(hi[0]))-D.from_float(float(lo[0]))
        improvement = REF-width
        print(f'{case}: lo={lo[0]} hi={hi[0]} width={width} wall={lo[1]}/{hi[1]}s '
              f'user={lo[2]}/{hi[2]}s sys={lo[3]}/{hi[3]}s RSS={lo[4]}/{hi[4]}KiB '
              f'last_shift={max(lo[5],hi[5]):.3E}')
        print(f'  vs prior D6 single: {improvement:+}, {improvement/REF:.6%} width reduction; '
              f'{improvement/(REF-CLIP):.6%} of prior clip gap; vs no L5: {(BASE-width)/BASE:.6%} width reduction')
