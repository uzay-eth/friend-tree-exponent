"""Independent cold checker for the excursion refinement (deliverable excursions-hy4).

Reads only logs. No author module is imported; all arithmetic below is written here.
"""
from fractions import Fraction as F
from math import ceil, floor
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
GAUGE = HERE.parent / 'finite' / 'gauge'
UNIT = 10**14
ROU = re.compile(r'--s (\S+) --test (\S+)')
ROW = re.compile(r'^(true row 1/2|true row 1/16|true positive dormant row) \[([^,]+),([^\]]+)\]')


def parse(path):
    out, active = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('$ '):
            m = ROU.search(line)
            active = (F(m[1]), m[2]) if m else None
        elif active:
            m = ROW.fullmatch(line)
            if m:
                out[active + (m[1],)] = (F(m[2]), F(m[3]))
            n = re.search(r'retained multi-core states (\d+) copy contraction <= (\d+) / (\d+)', line)
            if n:
                out[active + ('states',)] = (int(n[1]), F(int(n[2]), int(n[3])))
    return out


def ratio(lam, k, good=True):
    return F(k * k) / (k * k + (2 + lam) * k + lam) if good else F(k * k) / (k * k + (1 + lam) * k + lam)


def star_weights(lam, j, count, good=True):
    b = F(j - 1, j)
    w = b / (lam + j - 1 + b)
    yield w
    for k in range(j, j + count):
        w *= ratio(lam, k, good)
        yield w


def gamma_identity(lam, k, good=True):
    r = ratio(lam, k, good)
    g = lam + F(2) / (lam + 2)
    return (k + g) - (k + 1 + g) * r == (1 + lam) - lam * lam * r / ((lam + 2) * k * k)


def moments(lam, K, END=96):
    mass = F(1)
    pt = pm = F(0)
    for k in range(K, END):
        pt += mass
        pm += mass / k
        mass *= ratio(lam, k)
    g = lam + F(2) / (lam + 2)
    tt = mass * (END + g)
    tlo = tt / (1 + lam)
    thi = tt / (1 + lam - lam * lam / ((lam + 2) * END * END))
    mlo = mass / (2 + lam)
    mhi = mass / (END * (1 - ratio(lam, END)))
    return (pt + tlo, pt + thi, pm + mlo, pm + mhi)


def direct_sum(lam, K, BIG=2_000_000):
    lf = float(lam)
    mass = 1.0
    T = M = 0.0
    for k in range(K, BIG):
        T += mass
        M += mass / k
        mass *= k * k / (k * k + (2 + lf) * k + lf)
    return T, M


def check_excursion(fresh_path, recorded_path):
    fresh, old = parse(fresh_path), parse(recorded_path)
    assert set(fresh) == set(old), (set(fresh) ^ set(old))
    for key in sorted(fresh, key=str):
        assert fresh[key] == old[key], (key, fresh[key], old[key])
    print('fresh rows match recorded log exactly:', len(fresh), 'parsed quantities')

    for (s, test, label), (lo, hi) in sorted(fresh.items(), key=lambda kv: str(kv[0])):
        if label in ('true row 1/2', 'true row 1/16', 'true positive dormant row'):
            assert F(0) <= lo <= hi
    for key in [k for k in fresh if k[2] == 'states']:
        n, theta = fresh[key]
        assert theta < 1
        assert n == (954 if key[0] == F('.5665') or key[0] in (F('.54001'), F('.565615'), F('.565809')) else 156)
    print('all six cases: nonnegative ordered rows, contraction < 1, 954/156 state counts OK')

    for lam in (F('.5665') + F(1, 2), F('.5665') + F(1, 16), F(1) + F(1, 2)):
        for j in (2, 3, 16):
            ws = list(star_weights(lam, j, 12))
            assert abs(ws[0] - (F(j - 1, j) / (lam + j - 1 + F(j - 1, j)))) == 0
            if j == 2:
                assert ws[0] == F(1) / (3 + 2 * lam)
            assert all(a >= b for a, b in zip(ws, ws[1:]))
        for k in (32, 33, 95):
            assert gamma_identity(lam, k)
        iv = moments(lam, 32)
        direct = direct_sum(lam, 32)
        assert float(iv[0]) - 1e-8 <= direct[0] <= float(iv[1]) + 1e-8
        assert float(iv[2]) - 1e-8 <= direct[1] <= float(iv[3]) + 1e-8
        ak = [k * (1 - ratio(lam, k)) for k in (32, 40, 64, 95)]
        assert ak == sorted(ak) and all(a < 2 + lam for a in ak)
    print('star weights, gamma telescoping identity, moment enclosures, monotone a_k: PASS')

    assert not gamma_identity(F(1) + F(1, 2), 32, good=False)
    assert list(star_weights(F(1) + F(1, 2), 2, 1, good=False))[0] == F(1) / (3 + 2 * (F(1) + F(1, 2)))
    iv = moments(F(1) + F(1, 2), 32)
    direct = direct_sum(F(1) + F(1, 2), 32)
    swapped = (iv[1], iv[0], iv[3], iv[2])
    assert not (swapped[0] <= direct[0] <= swapped[1])
    print('negative controls (wrong ratio / swapped moment endpoints) rejected')

    kA = fresh[F('.5665'), 'A', 'true row 1/2'][0]
    kB = fresh[F('.5665'), 'B', 'true row 1/2'][0]
    dA = fresh[F('.5665'), 'A', 'true positive dormant row'][0]
    anchor = fresh[F('.54001'), 'B', 'true positive dormant row'][1]
    kA, kB, dA = (F(floor(v * 10**7), 10**7) for v in (kA, kB, dA))
    assert anchor < F('4.002')
    assert (kA, kB, dA) == (F('0.5225363'), F('0.0706297'), F('0.2995629'))
    critical = 1 - (1 + dA) * kB / (1 - kA)
    assert critical == F(38567596224187, 47746370000000) and critical < F('.807760')
    print('four fresh inputs -> critical excess', critical, '< .807760')

    u, a, s, cap = F('.5665'), F('.54001'), F('.5630'), F('.9980')
    p = (u - s) / (u - a)
    assert p == F(350, 2649)
    lhs = F('4.002')**p.numerator * F('.807760')**(p.denominator - p.numerator)
    assert lhs < cap**p.denominator
    w = F('.974') / (F(1, 16) + F('.14'))
    hmin = F('.501218471')
    factor = 1 + w / ((1 - cap) * hmin)
    assert ceil(factor) == 4800
    print('Holder at .5630 gives d_beta <', cap, '; multiplier ceil =', ceil(factor))

    assert 1 - kB / (1 - kA) > F('.807760')
    bad = 1 + w / ((1 - F('.9981')) * hmin)
    assert ceil(bad) > 4800
    badlhs = F('4.1')**p.numerator * F('.807760')**(p.denominator - p.numerator)
    assert not badlhs < cap**p.denominator
    print('negative controls (drop (1+d_A), cap .9981, anchor 4.1) all fail as required')


if __name__ == '__main__':
    fresh = HERE / 'excursions_hy4_replay.log'
    recorded = GAUGE / 'excursion-evidence.log'
    check_excursion(sys.argv[1] if len(sys.argv) > 1 else fresh,
                    sys.argv[2] if len(sys.argv) > 2 else recorded)
    print('EXCURSIONS COLD CHECK: PASS')
