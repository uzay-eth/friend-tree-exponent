"""Pointwise range of each archived motif term on projected endpoint classes."""

from fractions import Fraction as F
from functools import cache

from occupation import TAIL, coarse


class ScoreRange:
    def __init__(self, model, cert):
        self.c, self.model = cert, model
        self.vertices = [*model.finite, TAIL]
        self.parents = [*model.inner_parents, TAIL]
        self.p_ids, self.w_ids = {}, {}
        for s in model.finite:
            d, c = s
            i = cert.inner[d-c, c-1, 1]
            self.w_ids[s] = [i]
            if c >= 2:
                self.p_ids[s] = [i]
        ps, ws = {cert.INF}, {cert.INF}
        for d in range(model.cutoff+1, cert.D0+1):
            for c in range(1, d+1):
                ids = [i for i, _ in cert.mu((d-c, c-1, 1))]
                ws.update(ids)
                if c >= 2:
                    ps.update(ids)
        self.p_ids[TAIL], self.w_ids[TAIL] = sorted(ps), sorted(ws)
        self.edge = cache(self.edge)
        self.path = cache(self.path)
        self.cherry = cache(self.cherry)

    @staticmethod
    def extrema(rows):
        lo, hi = None, None
        for x in rows:
            if lo is None or x < lo:
                lo = x
            if hi is None or x > hi:
                hi = x
        return lo, hi

    def high_vertex(self, after=None):
        c = self.c
        cutoff = self.model.cutoff if after is None else after
        if cutoff >= min(c.LF, c.KF):
            raise ValueError('vertex lift beyond the range of this finite state scan')
        vals = (c.f[k*(c.LF+1)+l] for k in range(c.KF+1) for l in range(c.LF+1)
                if l+k+1 > cutoff)
        lo, hi = self.extrema(vals)
        overflow = c.f[2*(c.KF+1)*(c.LF+1)]
        return min(lo, overflow), max(hi, overflow)

    def edge(self, a, b):
        c = self.c
        return self.extrema(0 if i >= c.NI and j >= c.NI else c.g[i*c.NG+j]
                            for i in self.p_ids[a] for j in self.w_ids[b])

    def path(self, mid, a, b):
        c = self.c
        starts = [(c.pi[i]*c.nHU+c.hu[mid])*c.nW for i in self.p_ids[a]]
        ends = [c.wi[j] for j in self.w_ids[b]]
        return self.extrema(c.hV[start+j] for start in starts for j in ends)

    def cherry(self, apex, a, b):
        c = self.c
        base = c.cu[apex]*c.pairn
        norm = c.nW
        ix = [c.wi[i] for i in self.w_ids[a]]
        jx = [c.wi[j] for j in self.w_ids[b]]
        def rows():
            for x in ix:
                row = base+x*norm-x*(x-1)//2-x
                for y in jx:
                    yield c.hC[row+y if x <= y else base+y*norm-y*(y-1)//2+x-y]
        return self.extrema(rows())

    def hidden_middles(self, w):
        # A bottom has at most one vertical triple. Unordered cherries at a
        # hidden apex split equally over their two outgoing child incidences.
        c = self.c
        plo = phi = clo = chi = 0
        for mid in c.hu:
            if sum(mid) > self.model.cutoff:
                for a in self.parents:
                    lo, hi = self.path(mid, a, w)
                    plo, phi = min(plo, lo), max(phi, hi)
        for top in c.cu:
            if top[2] == 1 and sum(top) > self.model.cutoff:
                k = top[1]
                for b in self.vertices:
                    lo, hi = self.cherry(top, w, b)
                    clo, chi = min(clo, (k-1)*lo), max(chi, (k-1)*hi)
        return plo+F(clo, 2), phi+F(chi, 2)

    def bounds(self):
        m, c = self.model, self.c
        high = self.high_vertex()
        hidden = {w: self.hidden_middles(w) for w in self.vertices}
        bounds = {}
        for s, j in m.v.items():
            if s[0]:
                val = c.f[(s[1]-1)*(c.LF+1)+s[0]-s[1]]
                bounds[j] = val, val
            else:
                bounds[j] = high
        for (a, w), j in m.e.items():
            lo, hi = self.edge(coarse(a), coarse(w))
            if not a[0]:
                extra = hidden[coarse(w)]
                lo, hi = lo+extra[0], hi+extra[1]
            bounds[j] = lo, hi
        for (a, mid, w), j in m.t.items():
            d, k = mid
            bounds[j] = self.path((d-k, k-1, 1), coarse(a), coarse(w))
        for (mid, a, b), j in m.ch.items():
            d, k = mid
            bounds[j] = self.cherry((d-k, k-1, 1), coarse(a), coarse(b))
        return ({j: F(v[0], c.S) for j, v in bounds.items()},
                {j: F(v[1], c.S) for j, v in bounds.items()})
