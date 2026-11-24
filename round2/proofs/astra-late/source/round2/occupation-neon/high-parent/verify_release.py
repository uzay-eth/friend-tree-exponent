#!/usr/bin/env python3
"""Replay the reviewed size-20 and high-parent extension of Claim 21."""

import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import re
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OCC = 'round2/occupation-neon'
FILES = {
    'base': f'{OCC}/release/verify.py',
    'manifest': f'{OCC}/release/inputs.txt',
    'edges': f'{OCC}/rows/check_edges.py',
    'suffix': f'{OCC}/renewal/next/replay.py',
    'general': f'{OCC}/fringe/charges.py',
    'high_parent': f'{OCC}/high-parent/replay.py',
    'audit': 'round2/checks/neon_high_parent_sol.py',
}
EXPECTED = {
    'lower': (Q(58961819710386489, 10**23), Q(55203067615510485475692014463, 10**35),
              Q('76752873706217324196974475692014463/135683302131100000000000000000000000')),
    'upper': (Q(46635843154908967, 10**23), Q(139999781774155931356847859909, 5*10**35),
              Q('386753165449918989299233643152140091/683584997605000000000000000000000000')),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(key):
    spec = spec_from_file_location(f'extension_{key}', ROOT / FILES[key])
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def check_sources():
    from hashlib import file_digest

    pins = {}
    for line in (HERE / 'inputs.txt').read_text().splitlines():
        sha, name = line.split()
        require(name in FILES.values() and name not in pins and re.fullmatch('[a-f0-9]{64}', sha),
                f'invalid extension manifest entry: {name}')
        pins[name] = sha
    require(pins.keys() == set(FILES.values()), 'incomplete extension manifest')
    for name, sha in pins.items():
        with (ROOT / name).open('rb') as stream:
            require(file_digest(stream, 'sha256').hexdigest() == sha, f'source hash mismatch: {name}')


def run_side(base, modules, edges, suffix, hp, side, path, masses, ell):
    mode, trial, _, _, sha = base.PAYLOADS[side]
    require(base.digest(path) == sha, f'{side}: payload hash mismatch')
    print(f'{side}: original D18 SHA-256 verified: {sha}', flush=True)
    cert = modules['denom'].Cert(path)
    trial = Q(trial)
    require((cert.mode, cert.lam, cert.D, cert.K, cert.S) == (mode, trial, 18, 17, 10**11),
            f'{side}: unexpected witness layout')
    direction = 1 if side == 'lower' else -1

    # Reconstruct the first release's inputs, then take the conservative
    # floors and upward decimal ceiling used in the reviewed extension.
    first, _, _ = base.first_body(modules, cert, side)
    higher = modules['higher']
    higher.physical.MARK_CAP = suffix.P.MARK_CAP = 2
    higher_pay = higher.birth(higher.Extract(cert), 12)
    ceiling, fixtures, recharge = base.improved_denominator(modules, cert)
    require(first >= suffix.FIRST[side] > 0 and higher_pay >= suffix.OLD[side] > 0,
            f'{side}: inherited reward floors do not follow')
    den = suffix.DEN[side]
    require(den == Q(-(-ceiling*10**12 // 1), 10**12), f'{side}: wrong denominator rounding')
    print(f'{side}: inherited floors verified; ceiling={ceiling}, used={den}, '
          f'recharge={recharge}, {fixtures} charge fixtures', flush=True)

    edge_rows = edges.rows(modules['extract'].edges_buffer(cert))
    require(len(edge_rows) == 43776 and min(edge_rows.values()) == suffix.EDGE[side] > 0,
            f'{side}: inner-edge floor does not reproduce')
    del edge_rows
    # This joint solve pays higher rows only at sizes 13–20 and the actual
    # inner→inner edge floor throughout the retained cutoff, even below 13.
    size20 = suffix.floor(suffix.evaluate(suffix.H.Extract(cert), side, 20, masses, ell))

    # Mixed g has tail parent and a virgin child. hC has the same tail
    # parent and at least four children, hence at least five cores there.
    physical = hp.Physical(cert)
    rows = {(d, k): hp.charges(physical, d, k) for d in range(18, 192) for k in range(1, 11)}
    edge = lambda d, k, j: rows[d, k][0]*j if j else 0
    pair = lambda d, k, j: rows[d, k][1]*comb(j, 2) if j >= 2 else 0
    joint = hp.solve(18, 191, 10, lambda d, k, j: edge(d, k, j) + pair(d, k, j))
    del rows
    physical.value.cache_clear()

    paid = suffix.FIRST[side] + suffix.OLD[side] + size20
    endpoint = trial + direction*(paid + joint)/den
    require((size20, joint, endpoint) == EXPECTED[side], f'{side}: reviewed ledger did not reproduce')
    print(f'{side}: size20={size20}; high-parent joint={joint}; paid={paid + joint}', flush=True)
    print(f'{side}: exact endpoint={endpoint}; outward={base.outward(endpoint, direction)}', flush=True)
    return endpoint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lower', type=Path, required=True)
    parser.add_argument('--upper', type=Path, required=True)
    args = parser.parse_args()
    require(__debug__ and sys.version_info >= (3, 11) and sys.maxsize > 2**32,
            'need 64-bit Python 3.11+ with assertions enabled')
    sys.dont_write_bytecode = True
    check_sources()
    base = load('base')
    base.check_sources()
    base.check_receipts()
    paths = {side: getattr(args, side).resolve(strict=True) for side in EXPECTED}
    before = {side: base.fingerprint(path) for side, path in paths.items()}
    require(all(path.is_file() and before[side][2] == base.SIZE for side, path in paths.items()),
            'both payloads must be original D18 regular files')
    require(before['lower'][:2] != before['upper'][:2], 'payloads must be distinct')

    modules = {name: base.load(name, path) for name, path in base.MODULES.items()}
    edges, suffix, hp = map(load, ('edges', 'suffix', 'high_parent'))
    require((suffix.LO, suffix.HI) == (hp.LO, hp.HI) == base.PRIOR, 'physical prior mismatch')
    require(set(suffix.H.SPEC) - {'f', 'g'} == base.HIGHER, 'higher row family mismatch')
    masses, ell = suffix.newborn_kernel(False), suffix.newborn_kernel(True)
    lo, hi = (run_side(base, modules, edges, suffix, hp, s, p, masses, ell) for s, p in paths.items())
    require(before == {side: base.fingerprint(path) for side, path in paths.items()}, 'payload changed')
    check_sources()
    base.check_sources()
    print(f'PASS: size-20 + high-parent outward interval [{base.outward(lo, 1)}, {base.outward(hi, -1)}]')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, AssertionError, ImportError) as error:
        raise SystemExit(f'verify_release: {error}') from error
