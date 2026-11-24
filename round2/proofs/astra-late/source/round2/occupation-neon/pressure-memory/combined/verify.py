#!/usr/bin/env python3
"""Verify the pinned multi-family receipt ledger using exact arithmetic."""

import argparse
from fractions import Fraction as Q
from hashlib import file_digest
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OCC = Path('round2/occupation-neon')
EDGE = OCC / 'interpolate/edge'
COMPOSE = OCC / 'interpolate/compose'
MEM = COMPOSE / 'residual/edge-projection/combined'
SIGNS = ('lower', 'upper')
PAYLOADS = (
    '9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
    '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41',
)
PARTS = ('minimum', 'mother', '0', '1', '2', '3')
RATIONAL = r'(-?\d+(?:/\d+)?)'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open('rb') as f:
        return file_digest(f, 'sha256').hexdigest()


def manifest(path):
    pins = set()
    for line in (ROOT / path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        h, rel = line.split()
        p = Path(rel)
        require(not p.is_absolute() and '..' not in p.parts and p not in pins,
                f'unsafe or duplicate manifest entry: {rel}')
        require(digest(ROOT / p) == h, f'digest differs: {rel}')
        pins.add(p)
    require(pins, f'empty manifest: {path}')
    return pins


def unique(pattern, text):
    found = re.findall(pattern, text, re.M)
    require(len(found) == 1, f'expected one match, got {len(found)}: {pattern}')
    return found[0]


def base(text):
    text = text.split('=== producer-lower.log ===', 1)[1]
    rows = []
    for side, mode, sgn, sha in zip(SIGNS, (2, 1), (1, -1), PAYLOADS):
        trial, margin, reserve, den, edge = map(Q, unique(
            rf'^{side} authenticated {sha} trial {RATIONAL} margin {RATIONAL} '
            rf'reserved {RATIONAL} denominator {RATIONAL} exact EDGE {RATIONAL}$', text))
        old = Q(unique(rf'^{side} memory joint lower {RATIONAL} dual repair ', text))
        cold_old = Q(unique(rf'^{mode} memory exact joint floor {RATIONAL} repair numerator ', text))
        producer = Q(unique(rf'^{side} CANDIDATE exact {RATIONAL} outward ', text))
        cold_den, cold_end = map(Q, unique(
            rf'^{mode} den {RATIONAL} exact endpoint {RATIONAL} outward ', text))
        require(margin < 0 < old == cold_old and reserve > 0 and den == cold_den > 0,
                f'{side}: invalid original margin, reserve, price or denominator')
        net = margin + reserve
        require(net > 0 and trial + sgn * net / den == edge, f'{side}: EDGE composition')
        require(trial + sgn * (net + old) / den == producer == cold_end,
                f'{side}: memory composition')
        rows.append((trial, net, old, den, producer))
    return rows


def p4(text):
    unique(r'^original payloads authenticated; mask \[12, 8, 6, 16\] 3 donor/mother 2 4 theta 1$', text)
    unique(r'^physical D18 rows 14850 unpaid rows 12346 live states 196460\b', text)
    return [Q(unique(rf'^unpaid-{s} {RATIONAL} {RATIONAL} \[', text)[0]) for s in SIGNS]


def lift(text):
    unique(r'^authenticated both D18 originals degrees \[9, 7, 6, 10\] arities \[3, 3, 3, 2\] '
           r'side 3 mother 3 old included False$', text)
    headers = list(re.finditer(r'^component (\S+) censuses (\d+) live states (\d+) residuals (\d+)$', text, re.M))
    require(tuple(m[1] for m in headers) == PARTS, 'need all six complementary lift components exactly once')
    totals = [[Q(0), Q(0)] for _ in SIGNS]
    components = []
    for i, header in enumerate(headers):
        block = text[header.end():headers[i+1].start() if i+1 < len(headers) else len(text)]
        census, states, checks = map(int, header.groups()[1:])
        require(0 < census <= states <= 3000000 and checks == 4*states, 'invalid lift census or residual count')
        component = [census, states, checks]
        for side, total in zip(SIGNS, totals):
            lo, hi = map(Q, unique(rf'^{side} {RATIONAL} {RATIONAL} \[', block))
            require(0 <= lo <= hi, f'invalid {side} component price')
            component.append((lo, hi))
            total[0] += lo
            total[1] += hi
        components.append(tuple(component))
    for side, total in zip(SIGNS, totals):
        require(list(map(Q, unique(rf'^sum of listed parts {side} {RATIONAL} {RATIONAL} \[', text))) == total,
                f'{side}: lift totals do not match the six exact rows')
    return [total[0] for total in totals], components


def s3(text):
    unique(r'^PASS ordinary-tree support: 9/10/12 rejected, 13/15 accepted;', text)
    blocks = re.split(r'^(lower|upper) true host labels 10 row labels including tip age 40$', text, flags=re.M)
    require(tuple(blocks[1::2]) == SIGNS, 'expected corrected S3 lower then upper')
    floors = []
    for section in blocks[2::2]:
        unique(r'^unpaid k>=4 and true size>12 ' + RATIONAL + r' ', section)
        floors.append(Q(unique(r'^broader tip cap 5, 3360 true row classes ' + RATIONAL + r' ', section)))
    return floors


def s3_cold(text, side):
    require(re.findall(r'^(lower|upper) authenticated$', text, re.M) == [side],
            f'{side}: wrong independent S3 sign')
    return Q(unique(r'^stationary broad ' + RATIONAL + r' large-only ', finished(text)))


def yb(text):
    unique(r'^independent Bellman states 126813 row labels .*7056.*5865\}', text)
    return [Q(unique(rf'^disjoint {side} {RATIONAL} {RATIONAL} \[', text)[0]) for side in SIGNS]


def finished(text):
    require(unique(r'^wall=.* exit=(\d+)$', text) == '0', 'replay did not exit successfully')
    return text


def outward(value, upper):
    ticks = value * 10**15
    n = -((-ticks) // 1) if upper else ticks // 1
    return f'{n // 10**15}.{n % 10**15:015d}'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--payloads', nargs=2, type=Path, metavar=('LOWER', 'UPPER'))
    ap.add_argument('--compare-replays', type=Path, metavar='DIR', help='compare fresh exact replay receipts')
    ap.add_argument('--inputs-only', action='store_true', help='check manifests and original ledger without composing')
    args = ap.parse_args()
    pins = manifest(HERE.relative_to(ROOT) / 'inputs.sha256')
    inherited = set()
    for path in (MEM/'release-inputs.sha256', MEM/'inputs.sha256', EDGE/'inputs.sha256', COMPOSE/'inputs.txt'):
        require(path in pins, f'unpinned inherited manifest: {path}')
        inherited.update(manifest(path))

    def read(path):
        require(path in pins or path in inherited, f'unpinned receipt: {path}')
        return (ROOT / path).read_text()

    old = base(read(MEM / 'clean-replay.log'))
    if args.payloads:
        for p, sha in zip(args.payloads, PAYLOADS):
            require(p.stat().st_size == 3495091968 and digest(p) == sha, f'wrong original witness: {p}')
        print('PASS both original D18 payload digests')
    print(f'PASS {len(inherited)} inherited and {len(pins)} package input pins; exact EDGE and ead947b ledgers')
    if args.inputs_only:
        return

    package = HERE.relative_to(ROOT)
    sources = (
        ('root P4', p4, package/'p4.log', 'p4.log'),
        ('complementary P4 lift', lift, package/'donor-lift.log', 'donor-lift.log'),
        ('corrected S3', s3, OCC/'s3-lifetime/tagged-repair.log', 's3.log'),
        ('unused Yb', yb, package/'yb.log', 'yb.log'),
    )
    new = {}
    for name, parse, pin, fresh in sources:
        pinned = parse(finished(read(pin)))
        if args.compare_replays:
            observed = parse(finished((args.compare_replays/fresh).read_text()))
            require(observed == pinned, f'{name}: exact fresh receipt differs')
        new[name] = pinned[0] if parse is lift else pinned
    require([s3_cold(read(package/f's3-cold-{s}.log'), s) for s in SIGNS] == new['corrected S3'],
            'broad S3 floors disagree with independent replay')
    require(s3(finished(read(OCC/'s3-lifetime/disjoint.log'))) == new['corrected S3'],
            'tagged repair changed the historical broad S3 floor')
    unique(r'^PASS wide-row census: 14850 supported, 12346 unpaid, 118800 exact nonnegative sign corners$',
           finished(read(package/'p4-cold.log')))
    require(lift(read(OCC/'pressure-memory/p4-donor-medium.log')) == lift(read(package/'donor-lift.log')),
            'original and isolated medium lift receipts disagree')
    repaired = finished(read(package/'s3.log'))
    unique(r'^PASS 13-vertex upper negative regression:', repaired)
    require(s3(repaired) == new['corrected S3'],
            'isolated S3 producer floors disagree')
    if args.compare_replays:
        print('PASS all four fresh replay floor pairs match the pinned receipts exactly')
    ends = []
    for i, (side, sgn) in enumerate(zip(SIGNS, (1, -1))):
        trial, net, memory, den, prior = old[i]
        require(all(len(pair) == 2 and pair[i] > 0 for pair in new.values()), f'{side}: missing or nonpositive floor')
        delta = sum(pair[i] for pair in new.values())
        end = trial + sgn * (net + memory + delta) / den
        require(end == prior + sgn * delta / den, f'{side}: inconsistent numerator')
        print(side, 'EDGE net', net, 'old joint floor', memory, 'denominator', den)
        for name, pair in new.items():
            print(side, name, pair[i])
        ends.append(end)
    require(old[0][-1] < ends[0] < ends[1] < old[1][-1], 'enclosure does not improve both sides')
    print('exact lower', ends[0])
    print('exact upper', ends[1])
    print(f'outward [{outward(ends[0], False)}, {outward(ends[1], True)}]')


if __name__ == '__main__':
    main()
