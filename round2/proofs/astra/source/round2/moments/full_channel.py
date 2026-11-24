"""Transport all remaining active middle hubs onto the same bare-bottom path marginal."""

import lzma
import struct
from fractions import Fraction as F
from pathlib import Path

from occupation import BARE, ONE, ZERO, add, axpy, combine


class WholeChannel:
    def __init__(self, m):
        root = Path(__file__).resolve().parents[2]
        support = set()
        for side in ('upper', 'lower'):
            with lzma.open(root / 'verify/motif-bracket/reference' / (side+'.bin.xz')) as f:
                header = struct.unpack('<Qiiiiqqiiiiiiiiii9ii', f.read(120))
            support.add(tuple(header[11:13]))
        if len(support) != 1:
            raise ValueError('certificate path supports differ')
        degree, children = support.pop()
        if m.lift is None or m.lift.end < degree-1:
            raise ValueError('the vertex lift must cover every active path middle')
        self.states = [s for s in m.lift.states if s[0] < degree and 1 <= s[1]-1 <= children]
        self.rates = sorted(set(m.inv.values()))
        self.z = {}
        m.channel = self
        for s in self.states:
            for a in m.parents:
                for r in self.rates:
                    self.z[a, s, r] = m.variable(('high-middle-path', a, s, r))

        # Remove the resolved middle hubs from the upper end of the observed
        # high -> bare -> donor triples. A remaining parent either has inactive
        # degree or too many core children to support the archived vertical term.
        min_inactive = max(m.cutoff+1, min(degree, children+2))
        for r in self.rates:
            mass, reciprocal = {}, {}
            for h in (ZERO, ONE):
                for w in m.labels:
                    if m.inv[w] == r:
                        j = m.t[h, BARE, w]
                        add(mass, j, 1)
                        add(reciprocal, j, m.inv[h])
            for s in self.states:
                for a in m.parents:
                    j = self.z[a, s, r]
                    add(mass, j, -1)
                    add(reciprocal, j, -F(1, s[0]))
            m.ge(('inactive-mass', r), mass)
            m.ge(('inactive-reciprocal', r), reciprocal)
            m.ge(('inactive-degree', r), combine((1, mass), (-min_inactive, reciprocal)))
        for s in self.states:
            k = s[1]-1
            slots = {m.lift.v[s]: k}
            incoming = {m.lift.p[s]: k}
            bare_mass = {}
            for a in m.parents:
                for r in self.rates:
                    j = self.z[a, s, r]
                    add(slots, j, -1)
                    add(incoming, j, -m.inv[a])
                    add(bare_mass, j, 1)
            m.ge(('unused-hub-slots', s), slots)
            m.ge(('unused-hub-parent-input', s), incoming)
            m.ge(('unused-hub-parent-cap', s), combine((F(1, 2), slots), (-1, incoming)))
            m.ge(('bare-child-input', s), combine((1, {m.lift.ch[s]: 1}), (-F(1, 2), bare_mass)))
        for a in m.parents:
            for h in (ZERO, ONE):
                row = {m.e[a, h]: 1}
                for s in self.states:
                    dq = F(m.cutoff+1, s[0])
                    factor = (dq if h == ONE else 1-dq)/(s[1]-1)
                    for r in self.rates:
                        add(row, self.z[a, s, r], -factor)
                m.ge(('hub-parent-edge', a, h), row)

    def rewards(self, reward):
        m, c = reward.model, reward.c
        from rewards import donor, gain, state
        sb = state(BARE)
        upper, lower = {}, {}
        for s in self.states:
            sm = state(s)
            unresolved = reward.tail_delta(sm)
            for a in m.parents:
                if a[0]:
                    z = reward.sign*(c.HV(state(a), sm, gain(sb))-c.HV(state(a), sm, sb))
                    lo, hi = z, z
                else:
                    lo, hi = unresolved
                for r in self.rates:
                    j = self.z[a, s, r]
                    small = donor(lo, r) if lo > 0 else donor(hi, r) if hi < 0 else 0
                    add(lower, j, small)
                    add(upper, j, max(donor(lo, r), donor(hi, r)))
        reward.parts['high-middle-donor'] = lower
        reward.full_upper = dict(reward.tip_upper)
        axpy(reward.full_upper, upper)
