"""Physical certificate surplus and a pointwise upper ceiling for its denominator."""

import hashlib
import lzma
import os
import shutil
import sys
from fractions import Fraction as F
from functools import cache
from pathlib import Path

from ceiling import ScoreRange
from occupation import BARE, ONE, TIP, ZERO, add

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'verify/tree-drift'))
from motif_identity import Cert


def certificate(side):
    path = ROOT / 'verify/motif-bracket/reference' / (side + '.bin.xz')
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    with os.fdopen(os.memfd_create('moments-motif'), 'w+b') as raw, lzma.open(path) as src:
        shutil.copyfileobj(src, raw)
        raw.flush()
        cert = Cert('/proc/self/fd/' + str(raw.fileno()))
    if (side, cert.mode) not in {('lower', 2), ('upper', 1)}:
        raise ValueError('certificate direction mismatch')
    return cert, digest


def state(s):
    d, c = s
    return d-c, c-1, 1


def gain(s):
    return s[0]+1, s[1], s[2]


def bud(s):
    return s[0]-1, s[1]+1, s[2]


def donor(a, r):
    return max(a, 0) * (F(1, 2) - r) + max(-a, 0) * r


class Rewards:
    def __init__(self, model, side):
        self.c, self.digest = certificate(side)
        self.sign = 1 if side == 'upper' else -1
        self.model = model
        c = self.c
        if not 3 <= model.cutoff <= min(c.KCU+1, c.DCU-2, c.DHU-2):
            raise ValueError('reward projection requires all finite branching states active')
        for name in ('F', 'G', 'HV', 'HC'):
            setattr(c, name, cache(getattr(c, name)))
        self.parts = {name: {} for name in ('vertex', 'edge', 'path', 'cherry', 'tip-donor')}
        self.den = {}
        self.den_lower = {}
        self.tip_upper = {}
        self.fill()

    def cent(self, table, s):
        c = self.c
        return F(table[c.inner[s]], c.Sr)

    def vertex(self, s):
        c = self.c
        L, k, _ = s
        kept = self.cent(c.rpar, s) + k*self.cent(c.rch, s)
        drift = (L+kept)*(c.F(gain(s))-c.F(s))
        if L:
            drift += F(L, sum(s))*(c.F(bud(s))-c.F(s)+c.F(state(TIP))+c.G(bud(s), state(TIP)))
        return self.sign*(c.lam*c.F(s)-drift)

    def edge(self, a, b):
        c = self.c
        ia, ib = F(1, sum(a)), F(1, sum(b))
        h = c.G(a, b)
        dp, dw = c.G(gain(a), b)-h, c.G(a, gain(b))-h
        drift = (ib-self.cent(c.rch, a))*(c.F(gain(a))-c.F(a))
        drift += (ia-self.cent(c.rpar, b))*(c.F(gain(b))-c.F(b))
        assert a in c.hu and (a[1] < 2 or a in c.cu) and (b[1] == 0 or b in c.hu)
        drift += (a[0]+ib+self.cent(c.xpar, a)+(a[1]-1)*self.cent(c.xcc, a))*dp
        drift += (b[0]+ia+b[1]*self.cent(c.xch, b))*dw
        if a[0]:
            drift += F(a[0], sum(a))*(c.G(bud(a), b)-h+c.HC(bud(a), b, state(TIP)))
        if b[0]:
            drift += F(b[0], sum(b))*(c.G(a, bud(b))-h+c.HV(a, bud(b), state(TIP)))
        return self.sign*(c.lam*h-drift)

    def triple(self, family, s):
        c = self.c
        a, m, b = s
        ip, im, iw = (F(1, sum(x)) for x in s)
        if family == 'path':
            fn = c.HV
            transfer = (ip-self.cent(c.xpar, m))*(c.G(gain(m), b)-c.G(m, b))
            transfer += (iw-self.cent(c.xch, m))*(c.G(a, gain(m))-c.G(a, m))
            inside, count = (im, ip+iw, im), (1, 2, 1)
        else:
            # Here the first argument is the apex; m and b are its children.
            fn = c.HC
            bias = self.cent(c.xcc, a)
            transfer = (im-bias)*(c.G(gain(a), b)-c.G(a, b))
            transfer += (iw-bias)*(c.G(gain(a), m)-c.G(a, m))
            inside, count = (im+iw, ip, ip), (2, 1, 1)
        h = fn(*s)
        drift = transfer
        for j, v in enumerate(s):
            delta = fn(*(s[:j]+(gain(v),)+s[j+1:]))-h
            outside = (v[1]+1-count[j]) * F(1, 2) if self.sign*delta > 0 else 0
            drift += (v[0]+inside[j]+outside)*delta
            if v[0]:
                drift += F(v[0], sum(v))*(fn(*(s[:j]+(bud(v),)+s[j+1:]))-h)
        return self.sign*(c.lam*h-drift)

    def tail_delta(self, m):
        """All possible signed bottom increments for an unresolved top parent."""
        c = self.c
        b = state(BARE)
        def delta(p):
            return self.sign*(c.HV(p, m, gain(b))-c.HV(p, m, b))
        vals = [delta(state((d, k))) for d in range(self.model.cutoff+1, c.D0+1)
                for k in range(2, d+1)]
        # At infinity the top interpolation has the same entry for all child counts.
        base = (c.pi[c.INF]*c.nHU+c.hu[m])*c.nW
        vals.append(F(self.sign*(c.hV[base+c.wi[c.inner[gain(b)]]]-
                                 c.hV[base+c.wi[c.inner[b]]]), c.S))
        return min(vals), max(vals)

    def fill(self):
        model, c = self.model, self.c
        self.den_lower, self.den = ScoreRange(model, c).bounds()
        for s, j in model.v.items():
            if s[0]:
                add(self.parts['vertex'], j, self.vertex(state(s)))
        for (p, w), j in model.e.items():
            if p[0] and w[0]:
                add(self.parts['edge'], j, self.edge(state(p), state(w)))
        for (a, m, b), j in model.t.items():
            if a[0] and b[0]:
                row = tuple(map(state, (a, m, b)))
                add(self.parts['path'], j, self.triple('path', row))
        for (m, a, b), j in model.ch.items():
            if a[0] and b[0]:
                row = tuple(map(state, (m, a, b)))
                add(self.parts['cherry'], j, self.triple('cherry', row))
        intervals = {}
        for m in model.inner_parents:
            sm, sb = state(m), state(BARE)
            tail = self.tail_delta(sm)
            for a in model.parents:
                if a[0]:
                    h = self.sign*(c.HV(state(a), sm, gain(sb))-c.HV(state(a), sm, sb))
                    intervals[a, m] = (h, h)
                else:
                    intervals[a, m] = tail
        for (a, m, _, w), j in model.q.items():
            lo, hi = intervals[a, m]
            # The fourth site is the unique child of a bare degree-two bottom.
            credit = donor(lo, model.inv[w]) if lo > 0 else donor(hi, model.inv[w]) if hi < 0 else 0
            add(self.parts['tip-donor'], j, credit)
            add(self.tip_upper, j, max(donor(lo, model.inv[w]), donor(hi, model.inv[w])))
        for name, row in self.parts.items():
            if any(q < 0 for q in row.values()):
                raise AssertionError(('negative lower surplus', name, min(row.values())))
        if model.lift:
            self.lift_vertices()
        if model.channel:
            model.channel.rewards(self)

    def lift_vertices(self):
        m, c = self.model, self.c
        lift = m.lift
        lo, hi = ScoreRange(m, c).high_vertex(lift.end)
        lo, hi = F(lo, c.S), F(hi, c.S)
        for tail in (ZERO, ONE):
            self.den_lower[m.v[tail]], self.den[m.v[tail]] = lo, hi
        row = self.parts['vertex']
        for s, j in lift.v.items():
            sv = state(s)
            f = c.F(sv)
            add(self.den_lower, j, f-lo)
            add(self.den, j, f-hi)
            if c.isinner(sv):
                add(row, j, self.vertex(sv))
                continue
            d, core = s
            L = d-core
            df = c.F(gain(sv))-f
            drift = L*df
            if L:
                drift += F(L, d)*(c.F(bud(sv))-f+c.F(state(TIP))+c.G(bud(sv), state(TIP)))
            add(row, j, self.sign*(c.lam*f-drift))
            add(row, lift.p[s], -self.sign*df)
            if s in lift.ch:
                add(row, lift.ch[s], -self.sign*df)
