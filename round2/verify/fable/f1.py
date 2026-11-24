"""Replay F1's published finite-window statistics from selected chunk moments."""
from collections import defaultdict
from csv import DictReader
from decimal import Decimal, ROUND_CEILING
from math import exp, fsum, hypot, sqrt
from pathlib import Path
import sys

data = Path(sys.argv[1])
samples = defaultdict(list)
with data.open() as stream:
    for row in DictReader(stream, delimiter='\t'):
        rows = samples[row['sample']]
        assert int(row['chunk']) == len(rows)
        fields = [tuple(map(float, row[k].split(','))) for k in ('sums', 'cross_products')]
        rows.append((int(row['runs']), float(row['mu_ref']), *fields))


class Cohort:
    def __init__(self, name):
        rows = samples[name]
        self.ns, refs, self.sums, cross = zip(*rows)
        self.n = sum(self.ns)
        self.k = len(rows)
        self.ref, = set(refs)
        self.s = tuple(map(fsum, zip(*self.sums)))
        self.q = tuple(map(fsum, zip(*cross)))
        assert len(self.q) == len(self.s)**2

    def se(self, v):
        d = len(v)
        ss = fsum(v[i]*v[j]*self.q[i*d+j] for i in range(d) for j in range(d))
        mean = fsum(x*y for x, y in zip(v, self.s))
        return sqrt((ss-mean**2/self.n)*self.n/(self.n-1))

    def ratio(self, numerator, denominator, *, offset=True):
        e, d = (fsum(self.s[i] for i in ix) for ix in (numerator, denominator))
        assert d > 0
        ratio = e/d
        v = tuple(((i in numerator)-ratio*(i in denominator))/d for i in range(len(self.s)))
        leave = [(e-fsum(row[i] for i in numerator))/(d-fsum(row[i] for i in denominator)) for row in self.sums]
        mean = fsum(leave)/self.k
        return ratio + (self.ref if offset else 0), self.se(v), sqrt((self.k-1)/self.k*fsum((x-mean)**2 for x in leave)), v


m, a, b, l, c = map(Cohort, ('main','r1216','r1620','late','mainA'))
p = m.ratio((0,2), (1,3))
old = a.ratio((0,), (1,))
middle = b.ratio((0,), (1,))
late = l.ratio((4,), (5,))
w = m.ratio((1,), (1,3), offset=False)
for label, sample, r, expected in (
    ('primary', m, p, ('0.5657236606','6.15e-09','5.95e-09',19000)),
    ('old rung', a, old, ('0.5657237030','1.84e-08','1.27e-08',260000)),
    ('middle rung', b, middle, ('0.5657236424','1.88e-08','1.86e-08',25500)),
    ('late', l, late, ('0.5657236688','5.98e-09','6.26e-09',2840)),
):
    out = f'{r[0]:.10f}', f'{r[1]:.2e}', f'{r[2]:.2e}', sample.n
    assert out == expected, (out, expected)
    print(label, *out, '(ratio, run SE, leave-chunk SE, runs)')

w_up = w[0] + 1.96*w[1]
assert f'{w_up:.9f}' == '0.094247339'
print(f'w={w[0]:.9f}; upper normal limit={w_up:.9f}')


def bound(coefficient, target):
    k = Decimal.from_float(coefficient).quantize(Decimal('0.000001'), rounding=ROUND_CEILING)
    x = Decimal.from_float(float(k)*(abs(old[0]-target[0]) + 2.2414*hypot(old[1],target[1])))
    return f'{x.quantize(Decimal(1).scaleb(x.adjusted()-2), rounding=ROUND_CEILING):.2e}'


for shift, expected in [(0, ('1.73e-9','8.83e-11')), (0.30, ('1.42e-8','3.37e-9'))]:
    r = exp(-4*(0.565-shift))
    q = r*(w_up+(1-w_up)*r)
    out = bound(q/(1-q),p), bound(r**3/(1-r**3),late)
    assert out == expected, (out, expected)
    print(f'conditional bias (shift {shift:g}), primary/late:', *out)

diff = late[0]-p[0], hypot(late[1], p[1])
assert tuple(f'{x:.3e}' for x in diff) == ('8.222e-09','8.580e-09')
for sample, first, second, expected in [(c, (0,1), (2,3), '-0.16'), (l, (0,1), (2,3), '+1.18')]:
    x = sample.ratio((first[0],), (first[1],))
    y = sample.ratio((second[0],), (second[1],))
    z = (x[0]-y[0])/sample.se(tuple(v-u for v,u in zip(x[3],y[3])))
    assert f'{z:+.2f}' == expected, z
    print(f'ANPQ hold-out control: {z:+.2f} SE')
