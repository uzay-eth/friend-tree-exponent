"""Positive occupations and signed endpoint arithmetic for the pinned D18 rows."""
import json
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'finite'))
import residual_entries as physical

physical.MARK_CAP = 20
LO, HI = Q('0.565673570322'), Q('0.565775')
BASE = {'lower': Q('0.565673'), 'upper': HI}
C0 = {'lower': Q('1.84124847315'), 'upper': Q('1.95820675182')}
CEILING_BINS = ((2, 0), (3, 0), (4, 0), (5, 0), (2, 1), (3, 1), (4, 1))
YOUNG_CEILINGS = {
    'lower': (107297548088, 133110464918, 147478867063, 155384464355,
              59239266005, 100778864189, 125439375502),
    'upper': (107685711916, 134965659268, 149021514791, 156536518169,
              59817213941, 101372600114, 126953590473),
}
BANDS = {
    'lower': ('256918177771749/100000000000000000000', '4272535757583/4000000000000000000',
              '6713770119873/5000000000000000000', '932799090941/625000000000000000',
              '78795139786299/50000000000000000000'),
    'upper': ('20697884432099/10000000000000000000', '861071198001/1000000000000000000',
              '5412445357533/5000000000000000000', '6016226911043/5000000000000000000',
              '31764287097597/25000000000000000000'),
}
UNIT = 10**21


def down(x):
    return Q(x * UNIT // 1, UNIT)


def row_rewards(dump, name):
    item = dump[name]
    sign = {'lower': 1, 'upper': -1}[name]
    assert (Q(item['lambda']), item['sign'], item['S'], item['D0']) == (BASE[name], sign, 10**11, 18)
    out = {}
    for row in item['rows']:
        l, k, h = row['state']
        d = l + k + h
        assert h == 1 and (d, k) not in out
        v = row['f']
        gain = row['f_gain'] - v
        z = l * gain - BASE[name] * v
        if l:
            g = sum(y * Q(x) for y, x in row['g_bud_newborn'])
            z += Q(l, d) * (row['f_bud'] + row['f_newborn'] + g - v)
        p, c = map(Q, row['center'])
        center = p + k*c
        if d < 18:
            corners = [z + center * gain]
        else:
            t = Q(row['c2'] * d + row['c3'], d**3)
            corners = [z + center*t + rho*(gain-t) for rho in (0, Q(k+1, 2))]
        signed = [sign * x / item['S'] for x in corners]
        assert signed == list(map(Q, row['signed_corners'])) and min(signed) >= 0
        out[d, k] = (signed[0], signed[-1])
    assert set(out) == {(d, k) for d in range(2, 65) for k in range(d)}
    return out


def tracked(reward, band=(), limit=64):
    """Physical stable-age root; worst case over adaptive exterior input."""
    v = defaultdict(Q)
    tail = defaultdict(Q)
    for k in reversed(range(len(band))):
        r = band[k]
        tail[k] = down(min((r+b*tail[k+1])/(HI+b)
                           for b in (max(0, 1-Q(k+1, limit+1)), 1)))
        v[limit+1, k] = tail[k]
    for d in range(limit, 1, -1):
        for k in range(d-1, -1, -1):
            l = d-k-1
            b = Q(l, d)
            r = reward.get((d, k), (0, 0))
            v[d, k] = down(min((a + b*v[d, k+1] + g*v[d+1, k])/(HI+g+b)
                               for a, g in zip(r, (l, l+Q(k+1, 2)))))
    return LO*v[2, 0]


def fringe_occupation(cap):
    """Finite chronological prefix of residual paths staying at size <= cap."""
    # Repeated residual copy resets can return to smaller sizes. Sum all paths
    # through a chronological cutoff rather than sweeping by size.
    mass = {physical.newborn(2): Q(1)}
    spent = defaultdict(Q)
    for _ in range(cap + 3):
        nxt = defaultdict(Q)
        for t, m in mass.items():
            den, moves, *_ = physical.system(t, HI)
            w = down(m / den)
            spent[t] += w
            for child, rate in moves.items():
                if rate and physical.statistics(child)[0] <= cap:
                    nxt[child] += down(w*rate)
        mass = nxt
    return spent


def show(x):
    ticks = x*10**15 // 1
    return f'{ticks // 10**15}.{ticks % 10**15:015d}'


def main():
    dump = json.loads((HERE.parent/'rows'/'charged-f.json').read_text())
    assert dump['payload_sha256'] == {
        'lower': '9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
        'upper': '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41'}
    bins = defaultdict(Q)
    tip = defaultdict(Q)
    cherries = defaultdict(Q)
    trees = fringe_occupation(8)
    # Critical normalization A(b)=mu*A(1)=1 for ell(1)=1. The J_2
    # occupation is the least of the 19 rows for mark-blind scores.
    normalized = LO
    for t, x in trees.items():
        _, l, children = t
        d, k = l+len(children)+1, len(children)
        bins[d, k] += normalized*x
        tips = sum(ch[1:] == (1, ()) for ch in children)
        tip[d, k] += normalized*x*tips
        cherries[d, k] += normalized*x*tips*(tips-1)/2
    print(f'residual size 8, paths 0..10, {len(trees)} terminal states')
    for key in [(2, 0), (3, 0), (4, 0), (2, 1), (3, 2), (4, 2)]:
        print(key, 'root', show(bins[key]), 'virgin-child edge', show(tip[key]),
              'two-virgin cherry', show(cherries[key]))
    assert bins[2, 0] >= normalized*down(1/(HI+2))
    assert 0 < cherries[3, 2] <= bins[3, 2]
    endpoints = {}
    for name in ('lower', 'upper'):
        pay = row_rewards(dump, name)
        first = {(l+1, 0): pay[l+1, 0] for l in (1, 2, 3)}
        nominal = tracked(first)
        res = sum(min(pay[b])*m for b, m in bins.items())
        # Combine disjoint root bins by taking the best certified occupation
        # per bin or, independently, solve for the entire nonnegative reward.
        hybrid = sum(min(r)*max(bins[b], tracked({b: (1, 1)}, limit=b[0]))
                     for b, r in pay.items() if b[0] <= 8)
        root = tracked({b: (min(r),)*2 for b, r in pay.items()})
        joint_rho = tracked(pay)
        band = tuple(map(Q, BANDS[name]))
        reserve = tracked(pay, band)
        assert reserve >= joint_rho >= root
        print(name, 'first three', show(nominal), 'renewal', show(res),
              'hybrid D<=8', show(hybrid), 'tracked all bins', show(root),
              'joint rho', show(joint_rho), 'with all-time tail', show(reserve))
        a, b = band[:2]
        def floor(k):
            return a - (a-b)*k
        excess = {state: tuple(x-floor(state[1]) for x in r) for state, r in pay.items()}
        extail = tuple(x-floor(k) for k, x in enumerate(band))
        assert all(x >= 0 for r in excess.values() for x in r) and min(extail) >= 0
        affine = b + tracked(excess, extail)
        print(name, 'affine floor with discounted excess', show(affine))
        reserve = max(reserve, affine)
        gaps = {key: (C0[name]-Q(v, 10**11),)*2
                for key, v in zip(CEILING_BINS, YOUNG_CEILINGS[name])}
        assert all(x[0] > 0 for x in gaps.values())
        denom = C0[name]-tracked(gaps, limit=5)
        assert 0 < denom < C0[name]
        print(name, 'universal C0', show(C0[name]), 'occupation ceiling (rounded up)',
              show(-(-denom*10**15 // 1)/Q(10**15)),
              'gain using C0', show(reserve/C0[name]))
        move = reserve/denom
        edge = BASE[name] + (move if name == 'lower' else -move)
        ticks = edge*10**12
        rounded = Q(ticks // 1 if name == 'lower' else -(-ticks // 1), 10**12)
        endpoints[name] = rounded
        print(name, 'gain at Neon endpoint at least', show(move),
              'outward endpoint', show(rounded), 'reserve', show(reserve))
    joint = max(LO, endpoints['lower']), min(HI, endpoints['upper'])
    assert LO < joint[0] < joint[1] < HI
    print('joint:', *(show(x) for x in joint))


if __name__ == '__main__':
    main()
