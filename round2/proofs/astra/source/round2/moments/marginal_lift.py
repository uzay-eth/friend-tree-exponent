"""Refine high vertex degrees and both incoming reciprocal flows on the same edge law."""

from collections import defaultdict
from fractions import Fraction as F

from occupation import ONE, ZERO, add, axpy, combine


class HighDegrees:
    def __init__(self, model, end):
        if not model.cutoff < end:
            raise ValueError('vertex cutoff must exceed the joint cutoff')
        model.lift, self.end = self, end
        K = model.cutoff
        self.states = [(d, c) for d in range(K+1, end+1) for c in range(1, d+1)]
        self.v, self.p, self.ch = {}, {}, {}
        self.residual = {name: {} for name in ('mass', 'slots', 'inv', 'slotinv',
                                              'parent', 'child', 'parentinv', 'childinv')}
        rem = self.residual
        v, e, r = model.v, model.e, model.inv
        for h in [ZERO, ONE]:
            add(rem['mass'], v[h], 1)
            add(rem['inv'], v[h], r[h])
            for a in model.parents:
                add(rem['parent'], e[a, h], r[a])
                add(rem['parentinv'], e[a, h], r[a]*r[h])
            for b in model.labels:
                j = e[h, b]
                add(rem['slots'], j, 1)
                add(rem['slotinv'], j, r[h])
                add(rem['child'], j, r[b])
                add(rem['childinv'], j, r[b]*r[h])
        flow = defaultdict(dict)
        for c in range(1, K+1):
            axpy(flow[K+1, c], model.vertex_gain[K, c])
        for s in self.states:
            d, c = s
            # Every conditional child input is <= total edge reciprocal <= 1/2
            # after normalization; the per-vertex cap can be larger.
            j = self.v[s] = model.variable(('high-vertex', s))
            p = self.p[s] = model.variable(('high-parent-input', s), F(1, 2))
            kids = {}
            model.ge(('high-parent-cap', s), {j: F(1, 2), p: -1})
            if c > 1:
                h = self.ch[s] = model.variable(('high-child-input', s), F(1, 2))
                kids[h] = 1
                model.ge(('high-child-cap', s), {j: F(c-1, 2), h: -1})
            for name, coeff in (('mass', 1), ('slots', c-1), ('inv', F(1, d)),
                                ('slotinv', F(c-1, d))):
                add(rem[name], j, -coeff)
            add(rem['parent'], p, -1)
            add(rem['parentinv'], p, -F(1, d))
            axpy(rem['child'], kids, -1)
            axpy(rem['childinv'], kids, -F(1, d))
            gain = combine((1, {j: d-c, p: 1}), (1, kids))
            axpy(flow[s], gain, -1)
            axpy(flow[d+1, c] if d < end else flow['rest'], gain)
            if c < d:
                add(flow[s], j, -F(d-c, d))
                add(flow[d, c+1], j, F(d-c, d))
        for name, expr in rem.items():
            model.ge(('high-remainder', name), expr)
        for a, b, invcap in (
                ('inv', 'mass', end+1), ('slotinv', 'slots', end+1),
                ('parent', 'mass', 2), ('child', 'slots', 2),
                ('parentinv', 'parent', end+1), ('childinv', 'child', end+1),
                ('parentinv', 'inv', 2), ('childinv', 'slotinv', 2)):
            model.ge(('remainder-cap', a, b), combine((1, rem[b]), (-invcap, rem[a])))
        model.ge('remainder-birth', combine((1, rem['mass']), (-1, rem['inv']), (-1, rem['slotinv'])))
        for a, b, slots, moment in [('parent', 'parentinv', 'mass', 'inv'),
                                    ('child', 'childinv', 'slots', 'slotinv')]:
            model.ge(('remainder-complement', a),
                     combine((1, rem[slots]), (-(end+1), rem[moment]), (-2, rem[a]), (2*(end+1), rem[b])))
        for s in [*self.states, 'rest']:
            key = ('lift', s)
            occ = rem['mass'] if s == 'rest' else {self.v[s]: 1}
            model.mass[key], model.flow[key] = occ, flow[s]
            model.ge(('balance-low', key), combine((1, flow[s]), (-model.low, occ)))
            model.ge(('balance-high', key), combine((-1, flow[s]), (model.high, occ)))
