#!/usr/bin/env python3
"""Replay the disjoint two-sign Gold2 stationary frontier."""
import argparse
from fractions import Fraction as Q
from hashlib import file_digest, sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from mmap import ACCESS_READ, mmap
from os import fstat
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
GOLD = HERE.parents[1]
R2 = GOLD.parent
TARGET = Q(241, 426)
S = 2**40
UNIT = 10**6
SINGLETON = {
    'lower': Q(49262475234292952786766775195773, 148301530437500000000000000),
    'upper': Q(117851807982474570278951451121, 272874816005000000000000),
}
PINS = {
    'host/reader.py': '6be2e8aa5835ad0f1087179af5d84664d1176ed3ab071ca64f58af8d628fc910',
    'host/conditioned24.json': '3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790',
    'lifetime/pressure/frontier/replay.py': '1add7809bf6f509432b3931f2470e089eb8faf55ef0ca97a8828f18314ffa7df',
    'lifetime/pressure/frontier/marginals/lower-floor.json': '95b05a81de95b5ef559414d006a23acc82f4790126db69d8509fe488277f979a',
    'lifetime/pressure/frontier/marginals/upper-floor.json': '8b1b54a47a621e0f15006dfa5ec94aa26db53ce022e3502f793938ba271aae08',
    'lifetime/pressure/replay.py': 'a2d3a053ba785abad1a338bc929eae8c70cbdde367d4edce08dd93fc632a962c',
    'lifetime/check.py': 'ddae6be5f70649fcd778f872215d5d0b9826aecea061eefc198463369f756eca',
    '../attacks/neon2-quintet/gold2/tail/check.py': '76da7a850c85cc08d36d2900626e86ce4b05927823f313647ba5c0fee9a6f279',
    '../attacks/neon2-quintet/tail/check.py': 'aa2f7f59b693a2a1a83bc7f6df2929afcff134ddfc5b3e1c4f1dd0450534026d',
    '../occupation-quintet/census/census.py': 'b6c133acc92ef36d651a92f4684852fc5ed241e8c0e88d3ac0fdad5554c2550b',
    'edges/fy.py': '8263cce07bab4459cb6841bfa99b273cd2e3eff3da5860251baa119481028608',
    'edges/upper/replay.py': 'd33dd951fe030543e5bc2a9e6e23077e4fd64a41848fd4c7aaf8596a0077259e',
    'payer/price.py': '1a897871fa8713ba48690f5c6343ce516a83328c974fc5af3dbcd3749a0c822f',
    'payer/upper/replay.py': '59f1df6af82afa7280fb469d88523533e3abbb0658f7377d7ff27a212bcc5e20',
    'higher/cohort.py': 'b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c',
    'post-fy/upper-frontier/frontier.py': '1e2071a42f1b39c52c647f6a3b0f0205b4f12309d107aa8874e4b22afd7a9592',
    'post-fy/replay.py': 'd8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23',
    'edges/edge18.log': '1061d512219da5a38f264c07b4082f1db391b213a885f4cd955ea691946e2948',
    'edges/absorb12.log': 'eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521',
    'payer/README.md': '0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f',
    'edges/upper/ceiling8.log': 'c98252f27826882c085ece54e13e6294afb4dba2fcbf5977ca2416eb49e47bb3',
    'edges/upper/absorb15.log': '26c3491057f745ccbb7b86e87ea1b1a0c45311ea6f4d3afb498680aa78fe582c',
    'payer/upper/README.md': 'd2662a426eabb34f3d7046aea6945dc87a7f98479e2a9bd50ba5d0903d48b697',
    'host/frontier/replay-upper-h3.py': '49a4deeb6f8c3cf215ea7d093e2a2ab4c34d2a5af241bda30e0a3e3cf23ad0fe',
    'host/frontier/upper-h3.log': 'a195178330e6bc3addc6311315ab9c4cbf52d32fcec9b05e103f6575d3de402b',
    'host/frontier/replay-upper-increments.py': '480e62f6d24e83faa61f5a9b8d0f25693e6f29d07c773bdbf75a16bc9c937d2a',
    'host/frontier/upper-increments.log': '8842de55ece9fcf4b97681e57a550783b6bfd1ea426b252d9fbd227bf88117b8',
    '../checks/gold2-lower-h3-alt-sol.py': '609c89a5087786cf0a08db8f78b29771ca337ea0e9cf1a0a777695f5dcdfbbac',
    '../checks/gold2-lower-h3-alt-sol.md': '23e9bd2bda07ab5f31c5e95651459c55b0c2368a45ad30861e35a7b6400265c0',
    'diagnostic/upper-correction/threshold-upper.log': 'cad93785d631f3f8e95971f2a17398ab5929127a4e9a6e1448816ef55fb71415',
    'diagnostic/upper-correction/replay-threshold.py': '125b7b106273d8ef4796bf0183a98b081ccb6206652d57cc11cda93ef81c4c56',
    'true-drift/check.py': '1b95ebdbe76249b7afa07368202b71a2cf8fe2c2d3bf43cf7697cfc128210704',
}
ROOT_DIR = 'lifetime/pressure/frontier/marginals'
LP = '7dd40ac41e5d95ae41843da51de436fd45b781f1'
BASE_SHA = '9f9402eb89f711e5c12729a3fcfe50e31ee2e3ad9ee3e862a9dcfd55dbf3b03b'
SIBLING_PINS = {
    'diagnostic/lower-sibling/compact.log': 'c3d7d979fe82c4ac626cbdbef5f461411e966abd4e53ab35529a9cd5bb75c02a',
    'diagnostic/lower-sibling/price.py': '589e5d8d14f91621d0f5ad6bfecfdc00e4018a05bb3100416d0de8b22a4587f1',
}


def load(name, path, data):
    mod = module_from_spec(spec_from_file_location(name, path))
    exec(compile(data, str(path), 'exec'), mod.__dict__)
    return mod


def ledgers(data, profile, original):
    def ticks(path, pattern, scale=1):
        return scale*Q(re.search(pattern, data[path].decode(), re.M)[1])

    result = {
        'lower': {
            'edge18': ticks('edges/edge18.log', r'^GOLD2 edge18 normalized (\d+/\d+) ', S),
            'fy_small': ticks('edges/absorb12.log', r'^TOTAL (\d+/\d+) ', S),
            'fy_wide': ticks('payer/README.md', r'Exact normalized value: `(\d+/\d+)`', S),
        },
        'upper': {
            'edge18': ticks('edges/upper/ceiling8.log', r'^EDGE18 (\d+/\d+) ', S),
            'fy_small': ticks('edges/upper/absorb15.log', r'^TOTAL FY floor (\d+/\d+) ', S),
            'fy_wide': ticks('payer/upper/README.md', r'^\| `251/12` \| `35/6` \| floor \| `(\d+/\d+)`'),
        },
    }
    assert sum(result['lower'].values()) == Q(profile['lower']['edge_floor'])
    assert sum(result['upper'].values()) == Q(2288760217567599688149934647, 5*10**21)
    seam = load('double_tail_transform', GOLD/'true-drift/check.py', data['true-drift/check.py'])
    lower, upper = (Q(original[s]['certified_exact']) for s in ('lower', 'upper'))
    assert Q(seam.LO, seam.X) <= lower <= upper <= Q(seam.HI, seam.X)
    a, b = seam.host(252, False), seam.host(252, True)
    flux = max(Q(0), a-(10*b+seam.geom(10)/(1-seam.geom(1)))/252)
    result['upper']['edge_double_tail'] = original['upper']['cap_g']*flux*a/Q(seam.HI, seam.X)
    assert Q('0.3606267611') < result['upper']['edge_double_tail'] < Q('0.3606267612')
    audit = '../checks/gold2-lower-h3-alt-sol.md'
    result['lower']['h3_threshold'] = ticks(audit, r'^\| Mother `19\.\.79.*? \| `(\d+/\d+)`')
    result['lower']['h3_seam'] = ticks(audit, r'^\| Mother `80\.\.100.*? \| `(\d+/\d+)`')
    assert result['lower']['h3_threshold']+result['lower']['h3_seam'] == Q(919076249761939897477242027,65029541606195200000000)
    result['upper']['h3_threshold'] = ticks('diagnostic/upper-correction/threshold-upper.log',
                                           r'^corrected exact ticks/core (\d+/\d+)$')
    cohort = 'host/frontier/upper-h3.log'
    result['upper']['h3_mother_high'] = ticks(cohort, r'^corrected upper mother-high ticks/core (\d+/\d+)$')
    assert result['upper']['h3_mother_high'] == ticks(cohort, r'^order floor physical hmy (\d+/\d+) ', S)
    inc = data['host/frontier/upper-increments.log'].decode()
    assert 'PASS corrected original Gold2 upper H3 increments' in inc
    for name in ('middle3', 'stock'):
        normalized, floor = re.search(rf'^corrected upper increment {name}\n(?:[^\n]+\n)*?'
                                       r'normalized (\d+/\d+) \S+ ticks (\d+/\d+) ', inc, re.M).groups()
        assert S*Q(normalized) == Q(floor) > 0
        result['upper'][f'h3_mother_high_{name}'] = Q(floor)
    return result


def roots_for(side, rec, hi, frontier, reader):
    life = frontier.life
    path = reader.BASE/f'quintet_{side}.bin'

    def identity(st):
        return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns

    with path.open('rb') as fp:
        key = identity(fstat(fp.fileno()))
        if key[2] != life.tail.SIZE or file_digest(fp, 'sha256').hexdigest() != rec['verify']['sha256']:
            raise ValueError(f'wrong original {side} payload')
        with mmap(fp.fileno(), 0, access=ACCESS_READ) as raw:
            t = life.tail.Snapshot(raw, rec)
            assert t.s == (-1 if side == 'lower' else 1)
            W = {}
            life.price(t, 7999, 30, hi, rec['cap_f'], False, (701, W))
            roots, states, corners = frontier.pressure.live(t, rec['cap_f'], W, 700, 12, 5, 3, 2, hi, False)
        if identity(fstat(fp.fileno())) != key or identity(path.stat()) != key:
            raise ValueError(f'original {side} payload changed during replay')
    assert (states, corners) == (2800083, 17262)
    return roots


def frozen_roots(data, side, row, rec, original_sha):
    doc = json.loads(data[f'{ROOT_DIR}/{side}-floor.json'])
    assert (doc['side'], doc['bound'], doc['gold2_receipt_sha256'], doc['gold2_payload_sha256']) == (
        side, 'floor', original_sha, rec['verify']['sha256'])
    assert tuple(map(Q, doc['prior'])) == tuple(map(Q, row['prior']))
    assert Q(doc['discount']) == Q(row['prior'][1]) and Q(doc['overflow_ticks']) == 0
    assert (doc['states'], doc['corners']) == (2800083, 17262)
    assert doc['options'] == dict(mother=5, live_cut=700, live_arity=12, tip=3, donors=2, tag_limit=12)
    for name in ('lifetime/pressure/frontier/replay.py', 'lifetime/pressure/replay.py', 'lifetime/check.py'):
        assert doc['sources'][f'{GOLD.name}/{name}'] == PINS[name]
    roots = {tuple(map(int, k.split(','))): Q(v) for k, v in doc['roots'].items()}
    assert set(roots) == {(0, 0)} | {(D, c) for D in range(2, 6) for c in range(2, D+1)}
    assert all(v >= 0 for v in roots.values())
    assert Q(doc['results']['control']['singleton_ticks']) == SINGLETON[side]
    return roots


def solve(side, row, rec, ledger, frontier, reader, lp_blob, roots):
    d = 1 if side == 'lower' else -1
    lo, hi = map(Q, row['prior'])
    lam, E, H, M = Q(rec['nominal']), rec['err_per_core'], rec['higher_per_core'], rec['m_max']
    D = Q(row['denominator'][1])
    saved = Q(row['saved_plus_reduced'])
    cells = {tuple(map(int, k.split(','))): tuple(map(int, v)) for k, v in row['cells'].items()}
    assert len(cells) == 307 and 0 < D <= M and 0 <= saved <= H
    # The lower envelope includes the alpha recharge, the upper does not.
    assert all(rec['m_min'] <= l and u <= M and 0 <= h <= H for l, u, h in cells.values())
    B = sum(ledger.values())-E

    def build(low, high):
        model = reader.pinned('gold2_frontier_lp', R2/'occupation-neon/pressure-memory/birth-law/replay.py', lp_blob)
        model.LO, model.HI, model.SAT = low, high, 7
        lp = model.build(24, sided=True)
        frontier.tagged_births(lp, 12, low, high)
        assert (len(lp.rows), len(lp.cols)) == (623, 5564)
        return lp

    def mark(x):
        degree, k, b, *_ = x
        assert not b or degree > 5 or (degree, k+2) in roots
        return b*roots.get((degree, k+2), roots[0, 0])

    def price(lp, low, high, fn):
        z = UNIT*lp.lower(lambda kind, x: fn(x)/UNIT if kind == 'occ' else 0)
        return min(low*z, high*z)

    base = build(lo, hi)
    F = price(base, lo, hi, mark)
    assert F == SINGLETON[side] and F > Q(row['f_floor']) and B+F > 0

    def trial(s, restricted=False):
        assert lo <= s <= hi
        low, high = ((lo, s) if d == 1 else (s, hi)) if restricted else (lo, hi)
        a = d*(s-lam)
        assert a >= 0

        def objective(x):
            _, u, h = cells[x[:2]]
            assert u <= M and 0 <= h <= H
            return mark(x)+H-h-a*u

        lp = build(low, high) if restricted else base
        return B+price(lp, low, high, objective)

    scalar = lam+d*(B+F)/M
    separate = lam+d*(B+F+saved)/D
    margin = trial(separate)
    uniform = separate+d*max(0, margin)/D
    # Proposed points are grid rationals; only positive exact tests certify them.
    fixed = Q(d*((d*uniform*10**12).__ceil__()+1), 10**12)
    fixed_margin = trial(fixed, True)
    chosen = fixed if fixed_margin > 0 else uniform
    target_uniform_margin = trial(TARGET)
    target_margin = trial(TARGET, True)
    secondary = Q(99, 175)
    lower_secondary_margin = trial(secondary, True) if side == 'lower' else None
    print(side, 'uniform margin ticks/core', float(margin), 'fixed margin', float(fixed_margin),
          'target uniform/restricted margins', float(target_uniform_margin), float(target_margin),
          'outward', frontier.pressure.dec(chosen, d == -1, 15), flush=True)
    return dict(ledger={k: str(v) for k, v in ledger.items()}, singleton=str(F), replaced_singleton=row['f_floor'],
                singleton_roots={f'{degree},{c}': str(v) for (degree, c), v in roots.items()},
                original_prior=row['prior'], nominal=str(lam), payload_sha256=rec['verify']['sha256'],
                higher_cap=H, h3_family_cap=rec['verify']['caps']['h3'], original_upper_denominator=M,
                error=E, net=str(B+F), d24_upper_denominator=str(D), standalone_host_savings=str(saved),
                scalar_endpoint=str(scalar), separate_endpoint=str(separate), uniform_margin=str(margin),
                scalar_target_margin=str(B+F-d*(TARGET-lam)*M),
                uniform_endpoint=str(uniform), fixed_point=str(fixed), fixed_margin=str(fixed_margin),
                final_endpoint=str(chosen), outward=frontier.pressure.dec(chosen, d == -1, 15),
                target_uniform_margin=str(target_uniform_margin), target_margin=str(target_margin),
                lower_99_175_margin=None if lower_secondary_margin is None else str(lower_secondary_margin),
                inward_target_gap=str(d*(TARGET-chosen)))


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--replay-singleton', action='store_true')
    parser.add_argument('--live-sibling', action='store_true', help='separate lower replacement receipt')
    args = parser.parse_args()
    extra = SIBLING_PINS if args.live_sibling else {}
    pins = PINS | extra
    data = {name: (GOLD/name).read_bytes() for name in pins}
    for name, value in data.items():
        if sha256(value).hexdigest() != pins[name]:
            raise ValueError(f'changed frozen input: {name}')
    reader = load('reader', GOLD/'host/reader.py', data['host/reader.py'])
    frontier = load('marked_singleton', GOLD/'lifetime/pressure/frontier/replay.py', data['lifetime/pressure/frontier/replay.py'])
    original = reader.receipt()
    host = json.loads(data['host/conditioned24.json'])
    assert host['original_receipt_sha256'] == reader.RECEIPT
    profile = host['results']
    prior = tuple(Q(original[s]['certified_exact']) for s in ('lower', 'upper'))
    ledger = ledgers(data, profile, original)
    results = {}
    if args.live_sibling:
        frozen_base = (HERE/'receipt.json').read_bytes()
        assert sha256(frozen_base).hexdigest() == BASE_SHA
        results['upper'] = json.loads(frozen_base)['results']['upper']
        sibling = data['diagnostic/lower-sibling/compact.log'].decode()
        assert sibling.startswith(f"original GOLD2 lower {original['lower_sha256']} ")
        amounts = {name: Q(re.search(rf'^exact {name} ticks/core (\d+/\d+) ', sibling, re.M)[1])
                   for name in ('sibling', 'upstream')}
        assert amounts == {
            'sibling': Q(16363762149897743918379249111, 65029541606195200000000),
            'upstream': Q(1742793847212361294561574577, 65029541606195200000000),
        }
        del ledger['lower']['h3_threshold']
        ledger['lower'].update({f'h3_live_{name}': ticks for name, ticks in amounts.items()})
    for side, row in profile.items():
        if side in results:
            continue
        rec = original[side]
        assert (row['limit'], row['saturation'], row['reduced_f']) == (24, 7, False)
        assert tuple(map(Q, row['prior'])) == prior
        assert (row['payload_sha256'], row['nominal'], row['cap_error']) == (
            rec['verify']['sha256'], rec['nominal'], rec['err_per_core'])
        roots = frozen_roots(data, side, row, rec, reader.RECEIPT)
        if args.replay_singleton:
            assert roots == roots_for(side, rec, prior[1], frontier, reader)
        results[side] = solve(side, row, rec, ledger[side], frontier, reader, LP, roots)
    L, U = (Q(results[side]['final_endpoint']) for side in ('lower', 'upper'))
    assert prior[0] < L < TARGET < U < prior[1]
    doc = dict(target=str(TARGET), target_inside=L<TARGET<U, results=results, frozen_inputs=pins,
               original_receipt_sha256=reader.RECEIPT, original_source=reader.SOURCE,
               host_sha256=PINS['host/conditioned24.json'], lp_git_blob=LP)
    if args.live_sibling:
        assert Q(results['lower']['lower_99_175_margin']) > 0 and Q(99, 175) < L
        doc['base_receipt_sha256'] = BASE_SHA
    path = HERE/('receipt-live-sibling.json' if args.live_sibling else 'receipt.json')
    raw = (json.dumps(doc, indent=2, sort_keys=True)+'\n').encode()
    if args.write:
        path.write_bytes(raw)
    else:
        assert path.read_bytes() == raw
    print('outward', *(results[side]['outward'] for side in ('lower', 'upper')), 'receipt SHA-256', sha256(raw).hexdigest())
    print('PASS exact Gold2 disjoint frontier; 241/426 strictly inside the certified interval')


if __name__ == '__main__':
    main()
