#!/usr/bin/env python3
"""Verify the exact low-apex P4 replacement on the f7cc2e1 ledger."""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import runpy
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
M = Path('round2/occupation-neon/pressure-memory')
H = M / 'hv-lifetime/release'
C = M / 'combined'
R = M / 'path-fast'
PARTS = ('minimum', 'mother', '0', '1', '2', '3')
BATCHES = {
    'parity': PARTS, 'overflow': PARTS, 'medium': PARTS,
    'head': PARTS[:4], 'tail': PARTS[4:], 'mother12': ('mother',),
}
ENTRIES = {'medium': 'medium', 'head': 'wide-head', 'tail': 'wide-tail', 'mother12': 'wide-mother12'}
ROWS = r'^component (\S+) censuses (\d+) live states (\d+) edges (\d+) variable (\d+) residuals (\d+) slots (\d+) range (\([^\n]+\)) checksum ([0-9a-f]{16}) cpu=[\d.]+$'
RAT = r'([+-]?\d+(?:/\d+)?)'


def require(ok, what):
    if not ok:
        raise ValueError(what)


def one(pattern, text):
    result = re.findall(pattern, text, re.M)
    require(len(result) == 1, f'expected one row: {pattern}')
    return result[0]


def finished(text):
    require('terminated by signal' not in text and 'Command exited with non-zero status' not in text, 'interrupted receipt')
    require(one(r'^wall=[\d.]+ cpu=[\d.]+\+[\d.]+ rss=\d+ exit=(\d+)$', text) == '0', 'failed receipt')
    return text


def receipt(tag, text, digest, timed=True):
    parts = BATCHES[tag]
    wide = tag in ('head', 'tail', 'mother12')
    dims = [16, 9, 9, 16] if wide else [9, 7, 6, 10] if tag == 'medium' else [6, 5, 4, 6]
    arities = [4, 4, 4, 3] if wide else [3, 3, 3, 2] if tag == 'medium' else [2, 3, 2, 1]
    cap = 12 if tag == 'mother12' else 0 if tag == 'overflow' else 3
    header = (f'authenticated both D18 originals degrees {dims} arities {arities} side '
              f'{0 if tag == "overflow" else 3} mother {cap} old included {tag == "overflow"} library_sha256 {digest}')
    one(r'^' + re.escape(header) + r' auth_cpu=[\d.]+$', finished(text) if timed else text)
    matches = list(re.finditer(ROWS, text, re.M))
    require(tuple(w[1] for w in matches) == parts, f'{tag}: missing/repeated parts')
    results = []
    for i, (part, m) in enumerate(zip(parts, matches)):
        block = text[m.end():matches[i+1].start() if i+1 < len(matches) else len(text)]
        census, states, edges, var, checks, slots = map(int, m.groups()[1:7])
        require(0 < census <= states == slots and 0 <= var <= edges and checks == 4*states,
                f'{tag}/{part}: inconsistent state/drift counts')
        ends = tuple(tuple(map(Q, one(rf'^{side} {RAT} {RAT} \[', block))) for side in ('lower', 'upper'))
        require(all(0 <= lo <= hi for lo, hi in ends), f'{tag}/{part}: invalid endpoint order')
        if tag in ('parity', 'overflow'):
            one(rf'^python parity {part} exact fractions \+ counts$', block)
        results.append(((census, states, checks, *ends), (edges, var, slots, m[8], m[9])))
    requested, supported = map(int, one(r'^unique terminal charge requests (\d+) eligible and supported (\d+)$', text))
    require(0 < supported <= requested, f'{tag}: implausible terminal coverage')
    for i, side in enumerate(('lower', 'upper')):
        require(tuple(map(Q, one(rf'^sum of listed parts {side} {RAT} {RAT} \[', text))) ==
                tuple(sum(result[0][3+i][j] for result in results) for j in (0, 1)), f'{tag}: total mismatch')
    return results, (requested, supported)


def pool(text):
    one(r'^degrees \[16, 9, 9, 16\] arities \[4, 4, 4, 3\] mother 12 LP rows 462 cols 3820$', finished(text))
    box = tuple(tuple(map(Q, one(rf'^{side} shared-class uniform \[[^\n]+\] exact {RAT} {RAT}$', text)))
                for side in ('lower', 'upper'))
    require(all(0 < lo <= hi for lo, hi in box), 'pooled LP endpoints crossed')
    return box


def confirm_entries(data, tag, observed, digest):
    wide = tag != 'medium'
    cap = 12 if tag == 'mother12' else 3
    require((data['degrees'], data['arities'], data['mother_cap'], data['side_cap'], data['include_old_band']) ==
            ([16, 9, 9, 16] if wide else [9, 7, 6, 10], [4, 4, 4, 3] if wide else [3, 3, 3, 2], cap, 3, False),
            f'{tag}: wrong conditional support')
    lo, hi = Q('.565673'), Q('.565775')
    grid = 10**18
    require(data['library_sha256'] == digest and data['grid'] == grid and tuple(map(Q, data['discount'])) == (lo, hi) and
            data['channels'] == ['lower_lo', 'lower_hi', 'upper_lo', 'upper_hi'], f'{tag}: wrong conditional provenance')
    rows = data['conditional']
    require(tuple(rows) == BATCHES[tag], f'{tag}: incorrect conditional parts')
    for part, (summary, _) in zip(rows, observed):
        require(set(rows[part]) == {'overflow', *map(str, range(2, cap+1))}, f'{tag}/{part}: incorrect mother classes')
        values = rows[part].values()
        require(all(len(v) == 4 and all(type(x) is int and 0 <= x < 2**60 for x in v) and
                    v[0] <= v[1] and v[2] <= v[3] for v in values), f'{tag}/{part}: invalid conditional ticks')
        bounds = tuple((lo*min(v[j] for v in values)/grid, hi*max(v[j+1] for v in values)/grid) for j in (0, 2))
        require(bounds == summary[3:], f'{tag}/{part}: conditional ticks do not reproduce printed intervals')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--payloads', nargs=2, type=lambda p: Path(p).resolve())
    ap.add_argument('--compare-replays', type=Path)
    args = ap.parse_args()
    hv = runpy.run_path(str(ROOT / H / 'verify.py'))
    prior = runpy.run_path(str(ROOT / C / 'verify.py'))
    pins = hv['pins'](HERE / 'inputs.sha256')
    required = {H/'release.sha256', C/'release.sha256',
                *[R/f'{name}-entries.json' for name in ENTRIES.values()],
                *[(HERE/name).relative_to(ROOT) for name in ('full.log', 'pool.log', 'cold.log',
                  *[f'{tag}.log' for tag in BATCHES],
                  *[f'cold-{tag}.log' for tag in ('head', 'tail', 'mother12', 'medium-ubsan', 'wire')])]}
    require(required <= pins, 'baseline or receipt pins missing')
    hv['pins'](ROOT / H / 'release.sha256')
    hv['pins'](ROOT / C / 'release.sha256')
    baseline = subprocess.check_output(
        [sys.executable, '-EBu', str(ROOT / H / 'verify.py'),
         *(['--payloads', *map(str, args.payloads)] if args.payloads else [])], cwd=ROOT, text=True)
    digest = one(r'^([0-9a-f]{64})  .*path-fast\.so$', (HERE / 'full.log').read_text())
    full = (HERE / 'full.log').read_text()
    one(r'^PASS four byte-identical conditional receipt files$', full)
    one(r'^PASS pooled fresh receipts$', full)
    batches = {tag: receipt(tag, (HERE / f'{tag}.log').read_text(), digest) for tag in BATCHES}
    for tag in ('head', 'tail', 'mother12'):
        require(receipt(tag, (HERE/f'cold-{tag}.log').read_text(), digest, timed=False) == batches[tag],
                f'{tag}: isolated and cold native states disagree')
    ubsan = (HERE/'cold-medium-ubsan.log').read_text()
    ubsan_digest = one(r'^authenticated both D18 originals .* library_sha256 ([0-9a-f]{64}) .*$', ubsan)
    require(receipt('medium', ubsan, ubsan_digest, timed=False) == batches['medium'], 'cold UBSan medium states disagree')
    one(r'^PASS independent wide C\+\+/Fraction wire and UBSan exit/rounding fixtures: .+$', (HERE/'cold-wire.log').read_text())
    for tag, name in ENTRIES.items():
        confirm_entries(json.loads((ROOT / R / f'{name}-entries.json').read_text()), tag, batches[tag][0], digest)
    if args.compare_replays:
        for tag, saved in batches.items():
            require(receipt(tag, (args.compare_replays / f'{tag}.log').read_text(), digest) == saved,
                    f'{tag}: fresh counters, checksums or intervals differ')
        for name in ENTRIES.values():
            path = f'{name}-entries.json'
            require((args.compare_replays/path).read_bytes() == (ROOT/R/path).read_bytes(), f'fresh {path} differs')
    old, medium = prior['lift'](prior['finished']((ROOT/C/'donor-lift.log').read_text()))
    require([r[0] for r in batches['medium'][0]] == medium, 'medium historical endpoints or counts differ')
    pooled = pool((HERE / 'pool.log').read_text())
    if args.compare_replays:
        require(pool((args.compare_replays/'pool.log').read_text()) == pooled, 'fresh exact pool differs')
    cold = finished((HERE / 'cold.log').read_text())
    cold_price = one(r'^exact shared mother-only cap 12: .+$', cold)
    require(tuple(Q(int(n), int(d)) for n, d in re.findall(r'Fraction\((\d+), (\d+)\)', cold_price)) ==
            tuple(x for pair in pooled for x in pair), 'cold pool differs')
    cold_ends = tuple(map(Q, one(rf'^exact f7cc2e1 replacement: {RAT} {RAT} outward: [\d.]+ [\d.]+$', cold)))
    base_rows = prior['base']((ROOT / prior['MEM'] / 'clean-replay.log').read_text())
    ends, frozen = [], []
    for i, (side, sign) in enumerate((('lower', 1), ('upper', -1))):
        v = Q(one(rf'^exact {side} {RAT}$', baseline))
        prev_net, den, h = map(Q, one(rf'^{side} e911a1c net numerator {RAT} denominator {RAT} hV floor {RAT}$', baseline))
        trial = base_rows[i][0]
        require(den == base_rows[i][3] and v == trial + sign*(prev_net+h)/den, f'{side}: bad baseline ledger')
        replacement = pooled[i][0]
        require(replacement > old[i], f'{side}: replacement does not improve booked medium floor')
        net = prev_net + h - old[i] + replacement
        end = trial + sign*net/den
        require(end == v + sign*(replacement-old[i])/den == cold_ends[i], f'{side}: exact composition disagrees')
        print(side, 'old complementary P4', old[i], 'new pooled P4', replacement, 'net numerator', net, 'denominator', den)
        frozen.append(v)
        ends.append(end)
    require(frozen[0] < ends[0] < ends[1] < frozen[1], 'replacement endpoints do not strictly improve')
    print(f'PASS {len(pins)} low-apex input pins, six replay batches, exact pooled LP and independent composition')
    if args.payloads:
        print('PASS original D18 payload SHA-256s')
    print('exact lower', ends[0])
    print('exact upper', ends[1])
    print(f'outward [{prior["outward"](ends[0], False)}, {prior["outward"](ends[1], True)}]')


if __name__ == '__main__':
    main()
