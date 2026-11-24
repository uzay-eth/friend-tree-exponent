"""Native exact donor_lift replay using the authenticated mixed D18 tables."""

import argparse
from contextlib import redirect_stdout
import ctypes as C
from fractions import Fraction as Q
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import re
import sys
from time import process_time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from donor_lift import GRID, HIGH, LOW, Reward, eligible, interval, lifetime, original


PARTS = ('minimum', 'mother', '0', '1', '2', '3')
I, U = C.c_int64, C.c_uint64
LOAD = C.CFUNCTYPE(C.c_int, C.POINTER(C.c_int), C.POINTER(C.c_int), C.POINTER(I))
INDEX = C.POINTER(C.c_int)


class Tables(C.Structure):
    _fields_ = [(n, C.POINTER(I)) for n in ('p4', 'hv', 'rp', 'rc')] + [
        (n, INDEX) for n in ('parent', 'child', 'middle', 'triple', 'state')] + [
        (n, I) for n in ('nc', 'n4', 'n3', 'scale', 'rscale', 'lam', 'lscale', 'sign')]


def native_tables(tables):
    keep = []

    def encode(index):
        flat = (C.c_int * (18*18))(*([-1] * (18*18)))
        for label, j in index.items():
            if len(label) == 3 and label[2] == 1 and isinstance(label[0], int):
                l,k,h = label
                flat[(l+k+h)*18+k] = j
        keep.append(flat)
        return flat

    def wire(t, expected_sign):
        c = t.c
        if c.D != 18 or t.sign != expected_sign or c.lam != (LOW if t.sign == 1 else HIGH):
            raise ValueError('unexpected original D18 header')
        data = (c.p4tab,c.hvtab,c.rpar,c.rch)
        if sys.byteorder != 'little' or any(a.dtype != np.dtype('<i8') or not a.flags.c_contiguous
                                           for a in data):
            raise ValueError('original arrays must be contiguous little-endian signed i64')
        return Tables(*(a.ctypes.data_as(C.POINTER(I)) for a in data),
                      *(encode(index) for index in (c.par,c.child,c.p4,c.hv,c.state)),
                      len(c.child),len(c.p4),len(c.hv),c.S,c.Sr,int(c.lam*10**6),10**6,t.sign)

    return (Tables*2)(*(wire(t,s) for t,s in zip(tables,(1,-1)))), keep


def endpoints(entries):
    return [(max(Q(0), min(LOW*lo, HIGH*lo)/GRID), max(LOW*hi, HIGH*hi)/GRID)
            for j in (0, 2) for lo, hi in [(min(x[j] for x in entries), max(x[j+1] for x in entries))]]


def library(path):
    dll = C.CDLL(str(path.resolve(strict=True)))
    dll.prepare.argtypes = (INDEX, INDEX, C.c_int, C.c_int, C.POINTER(Tables), C.c_int, LOAD)
    dll.prepare.restype = C.c_void_p
    dll.price.argtypes = (C.c_void_p, C.c_int, U, C.POINTER(I), C.POINTER(U), C.POINTER(I))
    dll.price.restype = C.c_int
    dll.release.argtypes = (C.c_void_p,)
    dll.release.restype = None
    return dll


def wire_reward(tables, exclude_previous):
    """Price only reached terminal censuses; preserve the original mixed tables."""
    errors = []
    counts = [0, 0]

    @LOAD
    def load(d, k, dst):
        try:
            counts[0] += 1
            path = tuple((d[i], k[i]) for i in range(4))
            if not eligible(path, exclude_previous):
                return -1
            x = tuple((di-ki-1, ki, 1) for di, ki in path)
            supported = tuple(t.supported('P4', x) for t in tables)
            if supported[0] != supported[1]:
                raise ValueError('mismatched D18 P4 membership')
            if not supported[0]:
                return -1
            counts[1] += 1
            minima = []
            slopes = []
            for table in tables:
                base, free = table.row('P4', x)
                minimum = base / table.c.S
                slots = [Q(0)] * 5
                for i, parent, count, weight in free:
                    if (parent and i != 0) or count != (1 if parent else k[i]-int(i<3)):
                        raise ValueError('unexpected D18 exterior port')
                    slope = slots[0 if parent else i+1] = Q(weight, table.c.S)
                    minimum += count*min(Q(0), slope)/2
                if minimum < 0:
                    raise ValueError('negative physical pointwise minimum')
                minima.extend(interval(minimum))
                slopes.append(slots)
            out = minima
            for left, right in zip(*slopes):
                out.extend((*interval(max(0, -left)/2), *interval(max(0, -right)/2),
                            *interval(left), *interval(right)))
            if any(abs(z) >= 1 << 60 for z in out):
                raise OverflowError('reward interval outside native range')
            for j, z in enumerate(out):
                dst[j] = z
            return 1
        except Exception as exc:
            errors.append(exc)
            return 0

    return load, errors, counts


def parity(scores, part, a, native, stats):
    with redirect_stdout(receipt := StringIO()):
        expected = lifetime(scores, part, a.degrees, a.arities, a.side_cap, a.mother_cap, a.budget)
    match = re.fullmatch(r'component \S+ censuses (\d+) live states (\d+) residuals (\d+)\s*',
                         receipt.getvalue())
    if not match or tuple(map(int, match.groups())) != (stats[0], stats[1], stats[4]) or native != expected:
        raise ValueError(f'Python parity failed for {part}: {receipt.getvalue()} {expected} vs {native}')
    print('python parity', part, 'exact fractions + counts', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lower', type=Path, required=True)
    ap.add_argument('--upper', type=Path, required=True)
    ap.add_argument('--library', type=Path, required=True)
    ap.add_argument('--degrees', type=int, nargs=4, default=[6, 5, 4, 6])
    ap.add_argument('--arities', type=int, nargs=4, default=[2, 3, 2, 1])
    ap.add_argument('--side-cap', type=int, default=3)
    ap.add_argument('--mother-cap', type=int, default=3)
    ap.add_argument('--parts', nargs='+', choices=PARTS, default=PARTS)
    ap.add_argument('--budget', type=int, default=3000000)
    ap.add_argument('--include-old-band', action='store_true')
    ap.add_argument('--python-parity', action='store_true')
    ap.add_argument('--row-parity', action='store_true', help='compare every reached original charge with Extract')
    ap.add_argument('--entries', type=Path, help='write exact conditional birth-mother entries as JSON')
    args = ap.parse_args()
    if len(set(args.parts)) != len(args.parts) or not 0 < args.budget < 1 << 63:
        ap.error('invalid parts or budget')
    if any(not (0 <= k < d <= 16 and d >= 2 and k >= (i < 3))
           for i, (d, k) in enumerate(zip(args.degrees, args.arities))):
        ap.error('invalid path caps')
    if any(not (m == 0 or 2 <= m <= 23) for m in (args.side_cap, args.mother_cap)):
        ap.error('invalid memory caps')
    if GRID != 10**18 or (LOW, HIGH) != (Q('.565673'), Q('.565775')):
        raise ValueError('native constants no longer match donor_lift')
    t0 = process_time()
    tables = tuple(original(p, side) for side, p in [('lower', args.lower), ('upper', args.upper)])
    print('authenticated both D18 originals', 'degrees', args.degrees, 'arities', args.arities,
          'side', args.side_cap, 'mother', args.mother_cap, 'old included', args.include_old_band,
          'library_sha256', sha256(args.library.read_bytes()).hexdigest(),
          f'auth_cpu={process_time()-t0:.3f}', flush=True)
    load, errors, counts = (wire_reward(tables, not args.include_old_band)
                            if args.python_parity or args.row_parity else (LOAD(), [], None))
    scores = Reward(tables, not args.include_old_band) if args.python_parity else None
    dll = library(args.library)
    dims = C.c_int * 4
    tab, indices = native_tables(tables)  # keep index arrays alive while native pointers are borrowed
    handle = dll.prepare(dims(*args.degrees), dims(*args.arities), args.side_cap, args.mother_cap,
                         tab, not args.include_old_band, load)
    if not handle:
        raise RuntimeError('native prepare failed')
    totals = [[Q(0), Q(0)] for _ in range(2)]
    conditional = {}
    try:
        for part in args.parts:
            t0 = process_time()
            entries, stats, span = (I*(4*max(1,args.mother_cap)))(), (U*9)(), (I*2)()
            if dll.price(handle, PARTS.index(part)-2, args.budget, entries, stats, span):
                raise RuntimeError(f'native price failed for {part}') from (errors[-1] if errors else None)
            rows = list(zip(*[iter(entries)]*4))
            boxes = endpoints(rows)
            conditional[part] = dict(zip(('overflow', *map(str, range(2, args.mother_cap+1))), rows))
            npath, ns, ne, nv, nr, slots, checksum, requested, supported = stats
            print(f'component {part} censuses {npath} live states {ns} edges {ne} variable {nv} '
                  f'residuals {nr} slots {slots} range {tuple(span)} checksum {checksum:016x} '
                  f'cpu={process_time()-t0:.3f}', flush=True)
            for j, (lo, hi) in enumerate(boxes):
                totals[j][0] += lo
                totals[j][1] += hi
                print(('lower','upper')[j], lo, hi, f'[{float(lo):.12g}, {float(hi):.12g}]', flush=True)
            if scores is not None:
                parity(scores, part, args, boxes, stats)
    finally:
        dll.release(handle)
    if counts is not None and counts != [requested,supported]:
        raise ValueError('native/Python row counts differ')
    print('unique terminal charge requests', requested, 'eligible and supported', supported)
    for name,(lo,hi) in zip(('lower','upper'),totals):
        print('sum of listed parts', name, lo, hi, f'[{float(lo):.12g}, {float(hi):.12g}]')
    if args.entries:
        payload = dict(degrees=args.degrees, arities=args.arities, side_cap=args.side_cap,
                       mother_cap=args.mother_cap, include_old_band=args.include_old_band,
                       grid=GRID, discount=list(map(str,(LOW,HIGH))),
                       channels=['lower_lo','lower_hi','upper_lo','upper_hi'],
                       library_sha256=sha256(args.library.read_bytes()).hexdigest(),
                       conditional=conditional)
        with args.entries.open('x') as f:
            json.dump(payload, f, indent=2)
            f.write('\n')


if __name__ == '__main__':
    main()
