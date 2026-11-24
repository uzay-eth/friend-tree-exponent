#!/usr/bin/env python3
"""Replay the two D18 occupation bounds from frozen sources and external payloads."""

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import file_digest
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import re
import sys
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OCC = 'round2/occupation-neon'
BASE = 'round2/attacks/neon1'
ARCHIVE = f'{BASE}/outputs/neon1-production-rehearsal.tar.gz'
MODULES = {
    'crowns': f'{OCC}/check/crowns.py',
    'renewal': f'{OCC}/renewal/check.py',
    'young': f'{OCC}/renewal/denominator.py',
    'higher': f'{OCC}/higher/extract.py',
    'denom': f'{OCC}/denom/check.py',
    'extract': f'{OCC}/rows/extract.py',
    'stock': f'{OCC}/rows/check.py',
}
PAYLOADS = {
    'lower': (2, '565673/1000000', 113465364604, 184124847315,
              '9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec'),
    'upper': (1, '565775/1000000', 115462837236, 195820675182,
              '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41'),
}
SIZE = 3_495_091_968
SCALE = 10**11
PRIOR = Q('.565673570322'), Q('.565775')
HIGHER = frozenset(('hV', 'hC', 'P4', 'Yt', 'Yb', 'S3'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def check_sources():
    expected = {*MODULES.values(), ARCHIVE, f'{BASE}/audit.py',
                'round2/finite/residual_entries.py'}
    pins = {}
    for line in (HERE / 'inputs.txt').read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        sha, name = line.split()
        require(name in expected and name not in pins and re.fullmatch('[a-f0-9]{64}', sha),
                f'invalid source manifest entry: {name}')
        pins[name] = sha
    require(pins.keys() == expected, 'incomplete source manifest')
    for name, sha in pins.items():
        require(digest(ROOT / name) == sha, f'source hash mismatch: {name}')


def check_receipts():
    """Read only the two archived GMP receipts and their production session."""
    with tarfile.open(ROOT / ARCHIVE, 'r:gz') as archive:
        def receipt(name):
            with archive.extractfile(f'outputs/production-final/{name}.log') as stream:
                return stream.read().decode()

        session = receipt('session')
        require('cec3855cc52ab4265ce3e36b8a7a38f1ae257f77' in session and
                'Both pinned D0=18 production verifications and independent oracles passed.' in session
                and session.rstrip().endswith('run exit=0'), 'unexpected base production session')
        for side, (_, trial, c0_lower, c0_upper, _) in PAYLOADS.items():
            log = receipt('lo' if side == 'lower' else 'hi')
            rows = re.findall(r'^\s*\[(\w+)\s*\] checked ([1-9]\d*) states: fails = 0$', log, re.M)
            require(len(rows) == 8 and {name for name, _ in rows} == {'f', 'g', *HIGHER},
                    f'{side}: incomplete GMP row receipt')
            require(f'exact P/Q = {trial})' in log and 'PASS (total_fails = 0,' in log and
                    f'4*S*c0 = {c0_lower}, S*C0 = {c0_upper},' in log and
                    'on 34 concrete/boundary trees: fails = 0' in log and
                    log.rstrip().endswith('Exit status: 0'), f'{side}: invalid GMP receipt')
    print('Base receipts: both D18 signs, eight zero-failure sections and 34 tree checks each', flush=True)


def load(name, relative):
    spec = spec_from_file_location(f'release_{name}', ROOT / relative)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fingerprint(path):
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


@dataclass(frozen=True)
class Reward:
    label: str
    amount: Q
    families: frozenset[str]
    min_size: int = 2
    max_size: int | None = None


def combine(trial, direction, ceiling, rewards):
    require(direction in (-1, 1) and ceiling > 0, 'invalid sign or denominator')
    for i, term in enumerate(rewards):
        require(term.amount > 0 and bool(term.families) and term.min_size >= 2 and
                (term.max_size is None or term.max_size >= term.min_size), 'invalid reward')
        for other in rewards[:i]:
            apart = ((term.max_size is not None and term.max_size < other.min_size) or
                     (other.max_size is not None and other.max_size < term.min_size))
            require(apart or not term.families & other.families,
                    f'overlapping rewards: {other.label}, {term.label}')
    return trial + direction * sum(term.amount for term in rewards) / ceiling


def outward(edge, direction):
    ticks = direction * (direction * edge * 10**12 // 1)
    return f'{ticks // 10**12}.{ticks % 10**12:012d}'


def first_body(m, cert, side):
    # Re-extract the degree-64 operands and the entire supported stock slabs.
    slab = m['extract'].portable(cert)
    stock = m['stock'].table(slab)
    band = tuple(m['stock'].envelope(stock[k]) for k in range(5))
    r = m['renewal']
    require(band == tuple(map(Q, r.BANDS[side])), f'{side}: stock floors disagree')
    pay = r.row_rewards({side: m['extract'].buffer(cert)}, side)
    direct = r.tracked(pay, band)
    a, b = band[:2]
    require(a > 2*b > 0, f'{side}: invalid affine tail floor')
    floor = lambda k: a - (a-b)*k
    excess = {state: tuple(x-floor(state[1]) for x in pair) for state, pair in pay.items()}
    tail = tuple(x-floor(k) for k, x in enumerate(band))
    require(all(min(pair) >= 0 for pair in excess.values()) and min(tail) >= 0,
            f'{side}: negative affine residual')
    affine = b + r.tracked(excess, tail)
    bins = r.CEILING_BINS
    ceilings = tuple(sum(m['young'].charges(cert, (d-k-1, k, 1))) for d, k in bins)
    require(ceilings == r.YOUNG_CEILINGS[side], f'{side}: seven-bin extraction disagrees')
    gaps = {key: (r.C0[side] - Q(n, SCALE),)*2 for key, n in zip(bins, ceilings, strict=True)}
    seven = r.C0[side] - r.tracked(gaps, limit=5)
    print(f'{side}: extracted {len(pay)} first-body rows, 5 all-stock floors, 7 young ceilings', flush=True)
    return max(direct, affine), seven, frozenset(pay)


def improved_denominator(m, cert):
    d = m['denom']
    charge = d.Charge(cert)
    count = d.fixtures(cert, charge)
    c0 = d.C0[cert.mode]
    require(charge.cap(Q(0)) <= c0, 'grouped universal ceiling exceeds base C0')
    best = None
    for j in range(51):
        recharge = -Q(j*SCALE, 1000)
        cap = charge.cap(recharge)
        reward = {(degree, k): (cap - charge.exact(degree-k-1, k, 1) - recharge*(k-1))/SCALE
                  for degree in range(2, 65) for k in range(degree)}
        require(min(reward.values()) >= 0, 'negative potential-deficit reward')
        physical = d.physical_deficit(reward, 64)
        bound = cap/SCALE - physical
        if best is None or bound < best[0]:
            best = bound, cap, reward, recharge
    bound, cap, reward, recharge = best
    prefix, labels = d.renewal_deficit(reward, 8, 10)
    require(len(labels) == 19, 'incomplete residual occupation check')
    d.Charge.exact.cache_clear()
    return min(bound, cap/SCALE - prefix), count, recharge/SCALE


def run_side(m, side, path):
    mode, trial_text, _, c0, sha = PAYLOADS[side]
    trial = Q(trial_text)
    direction = 1 if side == 'lower' else -1
    # The sparse check also authenticates the entire original mapped payload.
    crowns = m['crowns']
    require(crowns.PINS[mode] == sha, f'{side}: payload pin disagreement')
    price = crowns.Price(path)
    try:
        require((price.lam, price.sign, price.M) == (trial, -direction, Q(c0, SCALE)),
                f'{side}: incorrect payload role or baseline')
        baseline = crowns.root_value(price, 64, 2)
        bonus = sum(crowns.rich_birth_mass(degree) * (crowns.root_value(price, 64, degree)-baseline)
                    for degree in range(3, 9))
        crown = crowns.LOW*(baseline + bonus)
        require(baseline > 0 and bonus >= 0 and crown > 0, f'{side}: invalid crown comparator')
        require(outward(trial + direction*crown/price.M, direction) ==
                {'lower': '0.565674338794', 'upper': '0.565773985336'}[side],
                f'{side}: crown comparator did not reproduce')
    finally:
        price.data.close()
        price.file.close()

    cert = m['denom'].Cert(path)
    require((cert.mode, cert.lam, cert.D, cert.K, cert.S) == (mode, trial, 18, 17, SCALE),
            f'{side}: unexpected numeric layout')
    first, seven, support = first_body(m, cert, side)
    high = m['higher']
    require(set(high.SPEC) - {'f', 'g'} == HIGHER, 'higher-family support changed')
    high.physical.MARK_CAP = 2
    extractor = high.Extract(cert)
    require(extractor.sign == direction, f'{side}: higher-row sign mismatch')
    high.validate(extractor, 8)
    higher = high.birth(extractor, 12)
    ceiling, fixtures, recharge = improved_denominator(m, cert)
    require(0 < ceiling < seven < Q(c0, SCALE), f'{side}: denominator did not improve')
    rewards = (Reward('first-body', first, frozenset({'f'})),
               Reward('size-at-most-12 higher', higher, HIGHER, max_size=12))
    previous = combine(trial, direction, seven, rewards)
    edge = combine(trial, direction, ceiling, rewards)
    require(direction * (edge-previous) > 0, f'{side}: combined endpoint did not improve')
    print(f'{side}: first={first}; higher={higher}; seven-bin ceiling={seven}', flush=True)
    print(f'{side}: improved ceiling={ceiling}; recharge={recharge}; {fixtures} charge fixtures', flush=True)
    print(f'{side}: exact combined endpoint={edge}; outward={outward(edge, direction)}', flush=True)
    return edge, support, tuple((r.families, r.min_size, r.max_size) for r in rewards)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lower', required=True, type=Path)
    parser.add_argument('--upper', required=True, type=Path)
    args = parser.parse_args()
    require(__debug__, 'run with assertions enabled')
    require(sys.version_info >= (3, 11) and sys.maxsize > 2**32, 'need 64-bit Python 3.11+')
    sys.dont_write_bytecode = True
    import numpy  # noqa: F401 — used by the pinned readers and reducers

    check_sources()
    check_receipts()
    paths = {side: getattr(args, side).resolve(strict=True) for side in PAYLOADS}
    before = {side: fingerprint(path) for side, path in paths.items()}
    require(all(path.is_file() and before[side][2] == SIZE for side, path in paths.items()),
            'both payloads must be regular original D18 files of 3495091968 bytes')
    require(before['lower'][:2] != before['upper'][:2], 'payloads must be different files')
    m = {name: load(name, path) for name, path in MODULES.items()}
    require(all(prior == PRIOR for prior in ((m['renewal'].LO, m['renewal'].HI),
                (m['crowns'].LOW, m['crowns'].HIGH), (m['higher'].LO, m['higher'].HI),
                (m['denom'].LOW, m['denom'].HIGH))), 'the occupation priors disagree')
    results = {side: run_side(m, side, path) for side, path in paths.items()}
    lower, upper = results['lower'], results['upper']
    require(lower[1:] == upper[1:], 'priced support differs between signs')
    require(PRIOR[0] < lower[0] < upper[0] < PRIOR[1], 'combined interval does not sharpen the prior')
    require(before == {side: fingerprint(path) for side, path in paths.items()},
            'payload file changed during the replay')
    check_sources()
    print(f'PASS: combined outward interval [{outward(lower[0], 1)}, {outward(upper[0], -1)}]')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, AssertionError, ImportError) as error:
        raise SystemExit(f'verify: {error}') from error
