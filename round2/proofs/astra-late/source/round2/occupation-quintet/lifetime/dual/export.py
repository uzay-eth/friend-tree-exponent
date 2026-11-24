#!/usr/bin/env python3
"""Export the exact low slices and reduced rewards of the frozen lifetime solve."""
import argparse
import hashlib
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from json import dumps
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMMIT = 'b04cc50c1e391e0420eb7b006ac781d484c68543'
RECEIPT_HASH = 'd26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb'
EXPECTED = {'lower': 525907178733179898, 'upper': 991668031155825523}
H, C, CUT = 12000, 30, 64


def load(name, path):
    spec = spec_from_file_location(name, path)
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def sweep(s, base):
    # Reuse the identical tail boundary; independently implement the scaled
    # quotient of the unchanged producer. Only the low slices are retained.
    nxt = base.tail_suffix(s, H+1, C)
    cap, hi, grid = base.GOLD[s.side]['cap_f'], base.HI, base.GRID
    values = {}
    hd = hi.denominator
    for d in range(H, 1, -1):
        cur = [0]*(C+2)
        w = 2*d*hd
        for c in range(min(d,C), 0, -1):
            a,b = s.row_ticks(d,c)
            a += cap
            L = d-c
            assert min(a,a+b*Q(c,2)) >= 0
            n = ((a*grid//1)+L*nxt[c])*w + 2*L*hd*cur[c+1]
            k = 2*d*hi.numerator+2*L*hd+L*w
            dx = c*d*hd
            m = n + dx*((b*grid//1)+nxt[c])
            v = min(n//k,m//(k+dx))
            assert 0 <= v and v*k <= n and v*(k+dx) <= m
            cur[c] = v
            if d <= CUT+1:
                values[d,c] = v
        nxt = cur
    assert values[2,1] == EXPECTED[s.side] and len(values) == 1514
    return values


def export(s, base):
    grid,hi = base.GRID,base.HI
    cap = base.GOLD[s.side]['cap_f']
    values = sweep(s,base)
    rows, checks = [], []
    for (d,c), v in sorted(values.items()):
        if d > CUT:
            continue
        L = d-c
        a,b = s.row_ticks(d,c)
        W = Q(v,grid)
        dv = Q(values[d+1,c],grid)-W
        # For c=d the bud rate vanishes; for c=C the continuation is zero.
        dw = Q(values.get((d,c+1),0),grid)-W
        intercept = cap+a+L*dv+Q(L,d)*dw-hi*W
        slope = b+dv
        for rho in (Q(0),Q(c,2)):
            exact = (cap+a+b*rho)+(L+rho)*dv+Q(L,d)*dw-hi*W
            assert exact == intercept+slope*rho >= 0
            checks.append(exact)
        rows.append([d,c,str(intercept),str(slope)])
    assert len(rows) == 1484
    print(s.side, 'V_J', Q(values[2,1],grid), 'values',len(values),'rows',len(rows),
          'corner checks',len(checks),'min',min(checks),'max',max(checks),flush=True)
    return dict(schema='quintet-f-dual-v1',side=s.side,source_commit=COMMIT,
                payload_sha256=base.GOLD[s.side]['verify']['sha256'],receipt_sha256=RECEIPT_HASH,
                scale=s.scale,cap_f_ticks=cap,**{'lambda':str(s.lam)},LO=str(base.LO),HI=str(base.HI),
                V_J_ticks=str(Q(values[2,1],grid)),V_grid=grid,degree=H,arity=C,cut=CUT,
                values=[[d,c,v] for (d,c),v in sorted(values.items())],rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replay',action='store_true',help='require byte-identical existing exports instead of writing')
    a=p.parse_args()
    base=load('quintet_lifetime_frozen',HERE.parent/'check.py')
    receipt=base.BASE/'source/round2/neon2/quintet-cert/certificate_gold1.json'
    assert hashlib.sha256(receipt.read_bytes()).hexdigest() == RECEIPT_HASH
    for side in ('lower','upper'):
        with base.Snapshot(base.BASE/f'quintet_{side}.bin',side) as s:
            ans=export(s,base)
            s.check_unchanged()
        # One row per line, without fragmenting an affine tuple over lines.
        lines = dumps({k:v for k,v in ans.items() if k not in ('values','rows')},indent=2).rstrip()[:-2]
        for key in ('values','rows'):
            lines += ',\n  '+dumps(key)+': [\n'+',\n'.join('    '+dumps(row) for row in ans[key])+'\n  ]'
        raw = (lines+'\n}\n').encode()
        dest=HERE/f'{side}.json'
        if a.replay:
            assert dest.read_bytes()==raw
        else:
            dest.write_bytes(raw)
        print(dest.name,hashlib.sha256(raw).hexdigest(),flush=True)


if __name__=='__main__':
    main()
