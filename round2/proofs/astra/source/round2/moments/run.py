#!/usr/bin/env python3
"""Build, solve or replay exact stationary occupation duals from the repository root."""

import argparse
import gzip
import json
from fractions import Fraction as F
from pathlib import Path

from dual import Solver, verify
from full_channel import WholeChannel
from marginal_lift import HighDegrees
from occupation import Occupation, TIP, combine
from proofcheck import IntegerProof
from rewards import Rewards
from resources import resources


def objectives(model, reward, shared, side):
    unit = {j: 1 for j in model.v.values()}
    if reward is None:
        tip = {model.v[TIP]: 1}
        return {'tip-hub-min': (model.tip_hub(), tip, 1),
                'tip-hub-max': (model.tip_hub(), tip, -1),
                'tip-parent-input-min': (model.tip_input(), tip, 1),
                'tip-parent-input-max': (model.tip_input(), tip, -1)}
    part = reward.parts['tip-donor']
    if model.channel:
        total = combine(*[(1, row) for row in reward.parts.values()])
        full_lower = combine((1, part), (1, reward.parts['high-middle-donor']))
        return {'entire-donor-surplus-min': (full_lower, unit, 1),
                'entire-donor-shift-max': (reward.full_upper, reward.den_lower, -1),
                'denominator-min': (reward.den_lower, unit, 1),
                'total-shift-min': (total, reward.den, 1),
                'potential-per-birth-max': (reward.den, model.nu, -1)}
    rows = {'donor-surplus-min': (part, unit, 1)}
    if shared:
        total = combine(*[(1, row) for row in reward.parts.values()])
        rows.update({'donor-surplus-max': (reward.tip_upper, unit, -1),
                     'denominator-min': (reward.den_lower, unit, 1),
                     'donor-shift-max': (reward.tip_upper, reward.den_lower, -1),
                     'total-shift-min': (total, reward.den, 1)})
        if model.lift:
            rows['potential-per-birth-max'] = reward.den, model.nu, -1
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['solve', 'check'])
    ap.add_argument('certificate', type=Path)
    ap.add_argument('--cutoff', type=int, default=4)
    ap.add_argument('--free-donor', action='store_true')
    ap.add_argument('--no-chain-balance', action='store_true')
    ap.add_argument('--vertex-cutoff', type=int, default=0)
    ap.add_argument('--full-channel', action='store_true')
    ap.add_argument('--sides', nargs='+', choices=('geometry', 'lower', 'upper'),
                    default=['geometry', 'lower', 'upper'])
    args = ap.parse_args()
    resources()
    solving = args.action == 'solve'
    if solving:
        sheet = dict(version=2, cutoff=args.cutoff, shared=not args.free_donor, chain=not args.no_chain_balance,
                     vertex_cutoff=args.vertex_cutoff, full_channel=args.full_channel,
                     sides=args.sides, hashes={}, duals={})
    else:
        with gzip.open(args.certificate, 'rt') as stream:
            sheet = json.load(stream)
        if sheet['version'] != 2:
            raise ValueError('unknown dual format')
    model = Occupation(sheet['cutoff'], sheet['shared'], sheet['chain'])
    if sheet.get('vertex_cutoff'):
        HighDegrees(model, sheet['vertex_cutoff'])
    if sheet.get('full_channel'):
        WholeChannel(model)
    print(f"cutoff={model.cutoff} vertex_cutoff={sheet.get('vertex_cutoff', 0)} "
          f"shared={model.shared} chain={model.chain} columns={len(model.names)} "
          f"equalities={len(model.equal)} inequalities={len(model.atleast)}", flush=True)
    solver = Solver(model) if solving else None
    integer_proof = IntegerProof(model) if not solving else None
    for side in sheet['sides']:
        reward = None if side == 'geometry' else Rewards(model, side)
        if reward:
            if solving:
                sheet['hashes'][side] = reward.digest
            if sheet['hashes'][side] != reward.digest:
                raise ValueError('certificate hash mismatch')
        rows = objectives(model, reward, model.shared, side)
        proven = {}
        for name, (numerator, denominator, direction) in rows.items():
            key = side + '/' + name
            if solving:
                sheet['duals'][key] = solver.ratio(numerator, denominator, direction)
            w = sheet['duals'][key]
            q = F(w['bound'])
            if direction != w['direction']:
                raise ValueError('wrong direction')
            obj = combine((direction, numerator), (-direction*q, denominator))
            margin = verify(model, obj, w)
            if margin < 0 or str(margin) != w['margin']:
                raise ArithmeticError(('negative or mismatched certificate margin', key, margin))
            if integer_proof and integer_proof.evaluate(obj, w) != margin:
                raise ArithmeticError(('integer verifier disagrees', key))
            proven[name] = q
            print(f'{key} {float(q):.9f} = {q} margin {float(margin):.3g}', flush=True)
        if 'total-shift-min' in proven:
            if proven['denominator-min'] <= 0:
                raise ArithmeticError('positive normalization required for shift ceiling')
            direction = 1 if side == 'lower' else -1
            print(f'{side} physical endpoint '
                  f'{reward.c.lam + direction*proven["total-shift-min"]}', flush=True)
    if solving:
        with args.certificate.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as out:
                out.write((json.dumps(sheet, separators=(',', ':'))+'\n').encode())


if __name__ == '__main__':
    main()
