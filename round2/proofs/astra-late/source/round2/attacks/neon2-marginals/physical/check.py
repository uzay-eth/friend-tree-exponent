#!/usr/bin/env python3
"""Independent exact arithmetic and physical arrows; no production payload scan."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import subprocess


def exact_table(root):
    archive = root / 'source/round2/neon2/stationary-marginals'
    data = json.loads((archive / 'certificates_all.json').read_text())
    auth = json.loads((root / 'authenticated.json').read_text())
    mu0, mu1 = map(F, data['mu_bracket_exact'])
    assert (mu0, mu1) == (F(565713486333, 10**12), F(565736986679, 10**12))
    iv = {}
    for i in range(1, 9):
        obs = data['observables'][str(i)]
        iv[i] = tuple(map(F, obs['certified_interval_exact']))
        for sign, side in [(-1, 'lower'), (1, 'upper')]:
            x = obs[side]
            v, t = x['verify'], x['tail']
            assert v['sha256'] == obs[side + '_sha256'] == auth[f'obs{i}_{side}.bin']['sha256']
            assert (v['obs_id'], v['obs_sign'], v['S'], v['DF']) == (i, sign, 2**40, 1600)
            err = []
            for endpoint in ['lo', 'hi']:
                row = v['caps_' + endpoint]
                cap = max(row['f_in'], row['f_tab'], t['cap_ticks_' + endpoint])
                err.append(max(0, cap + max(0, row['g']) + row['higher_per_core']))
            eps = F(max(err), 2**40)
            tail = F(40, 2**1600) / mu0 if i == 8 and sign == 1 else F(0)
            fin = F(v['theta_P'], v['theta_Q']) + sign*eps
            assert (F(x['eps_exact']), F(x['certified_fin_exact']), F(x['bare_tail_add_exact'])) == (eps, fin, tail)
            assert fin + tail == F(x['certified_exact']) == iv[i][sign == 1]
    print('PASS 16 authenticated receipt joins and exact endpoint combinations (no payload scan)')
    cor = data['corollaries']
    j0, j1 = iv[6]
    b0, b1 = iv[7]
    z0, z1 = iv[8]
    computed = {
        'pi_D_ge_7': (1-sum(iv[i][1] for i in range(1, 6)), 1-sum(iv[i][0] for i in range(1, 6))),
        'pi_bare_D_ge_3': (z0-b1, z1-b0),
        'E_invD_parent_given_J': (mu0*(1/j1-1)-F(3, 2), mu1*(1/j0-1)-F(3, 2)),
        'E_sum_invD_nbrs_given_2_2': (j0/(2*b1)-mu1, j1/(2*b0)-mu0),
    }
    for name, (lo, hi) in computed.items():
        assert (lo, hi) == tuple(map(F, cor[name]['exact']))
        left, right = map(F, cor[name]['dec'])
        assert 0 <= lo-left < F(1, 10**8) and 0 <= right-hi < F(1, 10**8)
        print(name, *cor[name]['dec'])
    assert j0+b0 <= iv[1][0] <= iv[1][1] <= j1+b1
    for k in range(1, 5):
        gap = iv[k][0] - iv[k+1][1]
        shown = F(cor['strict_monotonicity_p2_to_p6']['pairs'][f'p{k+1}_gt_p{k+2}']['certified_gap_dec'])
        assert 0 < shown <= gap < shown+F(1, 10**8)
    assert F('0.9015') < b0/z1 <= b1/z0 <= F('0.9038')
    assert F(40, 2**1600) / mu0 < F(1, 10**479)
    print('PASS degree-2 decomposition, four monotonicity gaps, bare fraction and analytic tail')
    # Consequence of 0 <= rho(J) <= 1/2, using only the published interval for J.
    assert F(3, 2)*j0/(1-j0) < mu0 < mu1 < 2*j1/(1-j1)
    print('J and reciprocal support alone allow mu in', float(F(3, 2)*j0/(1-j0)), float(2*j1/(1-j1)))


def marks(graph):
    return [(len(a), sum(len(graph[v]) >= 2 for v in a)) for a in graph]


def values(graph):
    m = marks(graph)
    counts = Counter(m)
    return (sum(D >= 2 for D, c in m), counts[2, 1], counts[2, 2],
            sum(D == c and D >= 2 for D, c in m))


def physical_rows():
    trees = arrows = 0
    for n in range(3, 8):
        for parents in product(*(range(v) for v in range(1, n))):
            g = [[] for _ in range(n)]
            for v, u in enumerate(parents, 1):
                g[u].append(v)
                g[v].append(u)
            before = values(g)
            drift = [F(0)] * 4
            cut = [F(0)] * (n+1)
            before_stock = Counter(D-c for D, c in marks(g) if D >= 2)
            for u, adj in enumerate(g):
                for v in adj:
                    after = [a[:] for a in g] + [[v]]
                    after[v].append(n)
                    for i, (z, x) in enumerate(zip(values(after), before)):
                        drift[i] += F(z-x, len(adj))
                    after_stock = Counter(D-c for D, c in marks(after) if D >= 2)
                    for j in range(n+1):
                        cut[j] += F(sum(after_stock[k]-before_stock[k] for k in range(j+1)), len(adj))
                    arrows += 1
            m = marks(g)
            rows = [(D, c, sum((F(1, len(g[v])) for v in g[u] if len(g[v]) >= 2), F(0)))
                    for u, (D, c) in enumerate(m) if D >= 2]
            birth = sum(F(D-c, D) for D, c, rho in rows)
            zj = sum(rho for D, c, rho in rows if (D, c) == (2, 1))
            zb = sum(rho for D, c, rho in rows if (D, c) == (2, 2))
            correction = int((2, 0) in m)
            assert drift[:3] == [birth, birth-F(3, 2)*before[1]-zj+correction, F(1, 2)*before[1]-zb]
            for j in range(n+1):
                expected = birth*(j >= 1) + sum(F(D-c, D)*(D-c == j+1)-(D-c+rho)*(D-c == j)
                                               for D, c, rho in rows)
                assert cut[j] == expected
            assert drift[3] == cut[0]
            trees += 1
    print(f'PASS {trees} increasing trees, {arrows} directed physical arrows, degree-2 and leaf-stock rows')


def plant(root, tiny):
    cmd = [str(root / 'verify_marginal'), str(tiny), '--threads', '1']
    base = json.loads(subprocess.check_output(cmd))
    hurt = json.loads(subprocess.check_output(cmd + ['--perturb', 'f']))
    assert hurt['perturbed'] and hurt['unperturbed_sha256'] == base['sha256'] != hurt['sha256']
    assert base['rows'] == hurt['rows'] and base['envelope_checks'] == hurt['envelope_checks']
    for k in ['lo', 'hi']:
        before, after = (v['caps_' + k] for v in (base, hurt))
        assert after['f_in'] - before['f_in'] > 10**9
        assert after['g'] > before['g']
        print(f'PASS J singleton plant {k}: f_in {before["f_in"]} -> {after["f_in"]}; g {before["g"]} -> {after["g"]}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('archive', type=Path)
    ap.add_argument('--tiny', type=Path)
    args = ap.parse_args()
    exact_table(args.archive)
    physical_rows()
    if args.tiny is not None:
        plant(args.archive, args.tiny)
