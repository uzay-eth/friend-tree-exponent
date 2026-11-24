"""Paired delta-method error in the two observed ratios of pooled totals."""
from math import sqrt
from pathlib import Path
from statistics import mean, stdev

pilot = Path(__file__).resolve().parents[2] / 'outputs' / 'depth-pilot.log'
by_time = {16: {}, 22: {}}
with pilot.open() as data:
    for line in data:
        row = line.split()
        if row[:1] == ['D'] and row[2] in ('16', '22') and row[3] == '1.000':
            by_time[int(row[2])][int(row[1])] = (int(row[4]), int(row[6]))
seeds = sorted(by_time[16])
assert len(seeds) == 350 and set(seeds) == set(by_time[22])
ratios, influence = [], []
for t in (16, 22):
    counts, depths = zip(*(by_time[t][i] for i in seeds))
    ratio = sum(depths)/sum(counts)
    ratios.append(ratio)
    influence.append([(x-ratio*n)/mean(counts) for n,x in zip(counts,depths)])
value = (ratios[1]-ratios[0])/6
error = stdev((b-a)/6 for a,b in zip(*influence))/sqrt(len(seeds))
print('paired seeds',len(seeds),'growth in pooled mean depth per time',f'{value:.12f}')
print('paired delta SE',f'{error:.12f}','standard errors above one half',f'{(value-.5)/error:.4f}')
