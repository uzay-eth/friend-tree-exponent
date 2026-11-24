#!/usr/bin/env python3
"""Authenticate the hV replay ledger and compose it with the exact e911a1c ledger."""

import argparse
import ast
from fractions import Fraction as Q
from hashlib import file_digest
from pathlib import Path
import re
import runpy
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
OCC = Path('round2/occupation-neon')
PRIOR = OCC / 'pressure-memory/combined'
SOURCE = OCC / 'pressure-memory/hv-lifetime/old'
BATCHES = {
    'apex20': ((20, 5, 5), (6, 2, 2), 4, 2, (('apex', 170010),)),
    'apex48': ((48, 5, 5), (6, 2, 2), 4, 21, (('apex', 463058),)),
    'parent20': ((20, 5, 5), (6, 2, 2), 2, 2, (('parent', 18470),)),
    'parent48': ((48, 5, 5), (6, 2, 2), 2, 21, (('parent', 49438),)),
    'sides': ((11, 10, 10), (3, 2, 2), 2, 2, (('middle', 76382), ('tip', 85103))),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pins(manifest):
    files = set()
    for line in manifest.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        expected, name = line.split()
        path = Path(name)
        require(not path.is_absolute() and '..' not in path.parts and path not in files,
                f'invalid manifest path: {name}')
        with (ROOT / path).open('rb') as f:
            require(file_digest(f, 'sha256').hexdigest() == expected, f'hash mismatch: {name}')
        files.add(path)
    require(files, f'empty manifest: {manifest}')
    return files


def one(pattern, text):
    match = re.findall(pattern, text, re.M)
    require(len(match) == 1, f'expected exactly one receipt row: {pattern}')
    return match[0]


def receipt(name, text, parent_weight):
    caps, arity, cap, low, roles = BATCHES[name]
    require(one(r'^wall=[\d.]+ rss=\d+ exit=(\d+)$', text) == '0', f'{name}: replay failed')
    require(len(re.findall(r'outward 1e-16 ticks', text)) == len(roles) + (roles[0][0] == 'parent'),
            f'{name}: unexpected receipt rows')
    out = []
    for role, states in roles:
        prefix = (f'{caps} {list(arity)} old min 8 atom {role} degree cap {cap} disjoint False '
                  f'ead947b complement True min apex {low} states {states} outward 1e-16 ticks ')
        ranges = ast.literal_eval(one(r'^' + re.escape(prefix) + r'(.+)$', text))
        require(len(ranges) == 3 and all(type(a) is int and type(b) is int and a <= b for a, b in ranges),
                f'{name}: invalid exact tick intervals')
        floor = tuple(pair[0] for pair in ranges[:2])
        if role == 'parent':
            got = tuple(map(int, one(r'^parent D_birth>=3 unconditional floors, outward 1e-16 ticks (\d+) (\d+)$', text)))
            # The Palm weight is applied to the lower conditional ticks, with directed rounding.
            require(got == tuple(parent_weight * value // 1 for value in floor), f'{name}: wrong Palm adjustment')
            floor = got
        require(all(x > 0 for x in floor), f'{name}: nonpositive floor')
        out.append((role, low, caps[0], floor, tuple(ranges)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--payloads', nargs=2, type=lambda x: Path(x).resolve())
    ap.add_argument('--compare-replays', type=Path, metavar='DIR')
    args = ap.parse_args()
    files = pins(HERE / 'inputs.sha256')
    require(PRIOR/'release.sha256' in files, 'missing pinned e911a1c release manifest')
    pins(ROOT / PRIOR / 'release.sha256')
    output = subprocess.check_output(
        [sys.executable, '-EBu', str(ROOT/PRIOR/'verify.py'),
         *(['--payloads', *map(str, args.payloads)] if args.payloads else [])], cwd=ROOT, text=True)
    r = r'(\d+/\d+)'
    old = tuple(Q(one(rf'^exact {side} {r}$', output)) for side in ('lower', 'upper'))
    prior = runpy.run_path(str(ROOT / PRIOR / 'verify.py'))
    base = prior['base']((ROOT / prior['MEM'] / 'clean-replay.log').read_text())
    den = tuple(row[3] for row in base)
    require(den == tuple(Q(one(rf'^{side} EDGE net {r} old joint floor {r} denominator {r}$', output)[2])
                         for side in ('lower', 'upper')), 'e911a1c denominators disagree')
    # The hV source imports this lower discount from higher/extract.py.
    text = (ROOT / OCC / 'higher/extract.py').read_text()
    lo = Q(one(r"^LO,HI=Q\('([.\d]+)'\),Q\('[.\d]+'\)$", text))
    weight = 1 - Q(1, 2)/(lo + Q(3, 2))
    rows = []
    for name in BATCHES:
        path = (HERE / f'{name}.log').relative_to(ROOT)
        require(path in files, f'unpinned replay receipt: {path}')
        paid = receipt(name, (ROOT/path).read_text(), weight)
        if args.compare_replays:
            require(paid == receipt(name, (args.compare_replays/path.name).read_text(), weight),
                    f'{name}: fresh replay differs')
        rows.extend(paid)
    require([(role, low, high) for role, low, high, _, _ in rows] == [
        ('apex', 2, 20), ('apex', 21, 48), ('parent', 2, 20), ('parent', 21, 48),
        ('middle', 2, 11), ('tip', 2, 11)], 'hV role/band partition differs')
    ticks = tuple(sum(row[3][i] for row in rows) for i in (0, 1))
    original = SOURCE/'EAD947B.md'
    require(original in files, 'missing pinned producer ledger')
    original_ticks = tuple(map(int, one(r'^\| \*\*Total\*\* \| \| \*\*(\d+)\*\* \| \*\*(\d+)\*\* \|',
                                        (ROOT/original).read_text())))
    require(ticks == original_ticks, 'hV floors differ from published source ledger')
    ends = []
    for i, (side, sgn) in enumerate((('lower', 1), ('upper', -1))):
        trial, edge_net, memory, d, _ = base[i]
        booked = tuple(Q(one(rf'^{side} {name} {r}$', output)) for name in
                       ('root P4', 'complementary P4 lift', 'corrected S3', 'unused Yb'))
        net = edge_net + memory + sum(booked)
        require(trial + sgn*net/d == old[i], f'{side}: e911a1c support ledger disagrees')
        hv = Q(ticks[i], 10**16)
        end = trial + sgn*(net+hv)/d
        require(end == old[i] + sgn*hv/d, f'{side}: exact residual arithmetic disagrees')
        print(side, 'e911a1c net numerator', net, 'denominator', d, 'hV floor', hv)
        ends.append(end)
    require(old[0] < ends[0] < ends[1] < old[1], 'hV does not improve both exact endpoints')
    print(f'PASS {len(files)} hV pins and e911a1c manifest; all six hV role/band rows agree')
    if args.payloads:
        print('PASS original D18 payload SHA-256s')
    print('exact lower', ends[0])
    print('exact upper', ends[1])
    print(f'outward [{prior["outward"](ends[0], False)}, {prior["outward"](ends[1], True)}]')


if __name__ == '__main__':
    main()
