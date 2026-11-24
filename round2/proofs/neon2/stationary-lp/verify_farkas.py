"""Exact box-Farkas checker. The verification kernel uses only the stdlib.

For Ax <= r and 0 <= x <= U, y >= 0 gives
    y A x >= sum_j min(0,(y A)_j) U_j.
A strict excess over y r is a contradiction. This is ordinary Farkas with
exact multipliers max(0,-(y A)_j) on the redundant box inequalities.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def fingerprint(rows, caps):
    h = hashlib.sha256()
    h.update((' '.join(map(str,caps))+'\n').encode())
    for name, expr, rhs in rows:
        h.update((f'{name!r}:{rhs}:'+','.join(f'{j}={x}' for j,x in sorted(expr.items()))+'\n').encode())
    return h.hexdigest()


def check(rows, caps, dual):
    if not dual or any(i<0 or y<0 for i,y in dual.items()):
        raise ValueError('invalid nonnegative multipliers')
    used, z, rhs = set(), {}, F(0)
    for i, (_, row, r) in enumerate(rows):
        if i not in dual:
            continue
        used.add(i)
        y = dual[i]
        rhs += y*r
        for j,x in row.items():
            if not 0<=j<len(caps) or isinstance(x,float):
                raise ValueError('invalid exact matrix')
            z[j] = z.get(j,0)+y*x
    if used != dual.keys():
        raise ValueError('unknown dual row')
    lhs = sum(x*caps[j] for j,x in z.items() if x<0)
    return rhs-lhs, rhs, lhs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('path')
    ap.add_argument('--exact',action='store_true',help='print the exact rational contradiction gap')
    ap.add_argument('--plant',choices=['promotion','donor','inside'])
    ap.add_argument('--replay-inside',action='store_true',help='test the multipliers at the interior point without trusting the saved hash')
    args = ap.parse_args()
    cert = json.loads(Path(args.path).read_text())
    a,b = map(F,cert['slice'])
    plant = args.plant
    if plant=='inside' or args.replay_inside:
        a=b=F('.56572366')
    from lp_level1 import Model
    from lp_level2 import MotifModel
    model = (MotifModel if cert['D3'] else Model)(cert['N'],cert['D3'],plant)
    if not args.replay_inside and fingerprint(model.exact_rows(a,b),model.cap)!=cert['fingerprint']:
        print('REJECT: exact physical system or slice differs from certificate')
        return 1
    dual = {i:F(y) for i,y in cert['dual']}
    if len(dual)!=len(cert['dual']):
        raise ValueError('duplicate dual row')
    margin, rhs, lhs = check(model.exact_rows(a,b),model.cap,dual)
    if args.replay_inside and margin<0:
        raise RuntimeError('control failed: exact contradiction at interior point')
    print(f'{"CERTIFIED" if margin<0 else "REJECT"}: N={model.N}, D3={model.D3}, '
          f'slice=[{a},{b}], multipliers={len(dual)}, exact gap~{float(margin):.12g}')
    if args.exact:
        print(margin)
    return int(margin>=0)


if __name__=='__main__':
    raise SystemExit(main())
