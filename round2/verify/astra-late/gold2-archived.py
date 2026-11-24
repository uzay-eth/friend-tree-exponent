#!/usr/bin/env python3
"""Replay the frozen Gold2 rational comparison from authenticated archival tables."""
from fractions import Fraction as Q
import json
from math import ceil
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parents[2] / 'proofs/astra-late/source/round2'
g = root / 'occupation-quintet-gold2'
writer = run_path(str(g / 'rationals/repaired/check.py'))
cold = run_path(str(root / 'checks/gold2-rational-repair-sol.py'))
require, pinned = cold['require'], writer['read_pinned']
receipt = json.loads(pinned(root / 'attacks/neon2-quintet/gold2/production.json', writer['RECEIPT_PIN']))
conditioned = json.loads(pinned(g / 'host/conditioned24.json', writer['CELLS_PIN']))
require(conditioned['original_receipt_sha256'] == writer['RECEIPT_PIN'], 'original Gold2 receipt')
row = conditioned['results']['lower']
upper = {(d, c): int(row['cells'][f'{d},{c-1}'][1])
         for d in range(2, 25) for c in range(1, d+1)}
small = pinned(g / 'edges/absorb12.log', 'eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521')
wide = pinned(root.parents[1] / 'authenticated/round2/occupation-quintet-gold2/payer/README.md.txt',
              '0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f')
sibling = pinned(g / 'diagnostic/lower-sibling/compact.log', 'c3d7d979fe82c4ac626cbdbef5f461411e966abd4e53ab35529a9cd5bb75c02a')
match = writer['match_fraction']
scale = writer['SCALE']
floors = tuple(map(int, (Q(row['f_floor']), scale * match(small, r'^TOTAL (\d+/\d+) '),
                        scale * match(wide, r'Exact normalized value: `(\d+/\d+)`'),
                        match(sibling, r'^exact sibling ticks/core (\d+/\d+) '))))
require(floors == cold['FLOORS'], 'four archived reward floors')
lo, hi = (Q(receipt[side]['certified_exact']) for side in ('lower', 'upper'))
lower = receipt['lower']
ceiling, error, nominal = (lower['m_max'], lower['err_per_core'], Q(lower['nominal']))
require((ceiling, error, nominal) == (1808603193693, 19217, Q(7071421, 12500000)), 'Gold2 lower data')
net = sum(floors) - error
require(net == 888010, 'net reward')
count, arrows = cold['physical']()
print(f'PASS {count} literal core drift fixtures, {arrows} arrows')
ends = []
for cut, needed in ((6, 4700), (18, 25700), (24, 27680)):
    v, rows = cold['bellman'](upper, cut, ceiling, hi)
    require(v == writer['occupation_saving'](cut, ceiling, upper, hi), 'matching Bellman values')
    bound = ceil(ceiling - lo*v)
    margin = net - (cold['TARGET'] - nominal)*bound
    lower = nominal + Q(net, bound)
    require(0 < bound < 1459000000000 and margin > needed, 'repaired margin')
    require(cold['TARGET'] < lower <= hi < Q(142, 251) and 142*175-99*251 == 1, 'Farey interval')
    require(cold['fractions'](lower, hi, 425) == [] and
            cold['fractions'](lower, hi, 426) == [Q(241, 426)], 'denominator scan')
    ends.append(lower)
    print(f'PASS D<={cut}: {rows} Bellman corners; M<={bound}; margin={margin}')
require(ends == sorted(ends) and ends[-1] == Q(10043824822832044427, 17754235186087500000),
        'nested cutoffs and pinned repair endpoint')
print(f'PASS archived-input rational exclusion through 425: [{ends[-1]}, {hi}]')

frontier = g / 'host/frontier'
live = json.loads(pinned(frontier / 'receipt-live-sibling.json',
                         '7a6ddda08b83574f6b745f11b3e5445ba80a60f47d7ef75c625b093a56cad20f'))
base = json.loads(pinned(frontier / 'receipt.json',
                         '9f9402eb89f711e5c12729a3fcfe50e31ee2e3ad9ee3e862a9dcfd55dbf3b03b'))
require(live['original_receipt_sha256'] == writer['RECEIPT_PIN'] and
        live['host_sha256'] == writer['CELLS_PIN'], 'frozen frontier dependencies')
require(live['results']['upper'] == base['results']['upper'], 'unchanged corrected upper result')
lf, uf = (Q(live['results'][side]['final_endpoint']) for side in ('lower', 'upper'))
require(lo < Q(99, 175) < Q('0.565714467904') == lf < Q(241, 426) < uf ==
        Q('0.565735974375') < hi < Q(142, 251), 'archived frontier inclusion and decimal endpoints')
require(all(Q(live['results'][side]['fixed_margin']) > 0 for side in ('lower', 'upper')),
        'recorded strict margins')
print('PASS archived frontier receipt hashes and decimal endpoints: 0.565714467904, 0.565735974375')
