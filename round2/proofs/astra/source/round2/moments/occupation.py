"""Rational projected stationary occupation of physical core vertices and their overlaps."""

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement


LOW, HIGH = F(113, 200), F(1133, 2000)
TIP, BARE = (2, 1), (2, 2)
TAIL, ZERO, ONE = (0, -1), (0, 0), (0, 1)


def coarse(s):
    return s if s[0] else TAIL


def add(out, key, value):
    if value:
        out[key] = out.get(key, 0) + value
        if not out[key]:
            del out[key]


def axpy(out, other, scale=1):
    for key, value in other.items():
        add(out, key, scale * value)


def combine(*terms):
    result = {}
    for scale, expr in terms:
        axpy(result, expr, scale)
    return result


class Problem:
    def __init__(self):
        self.names = []
        self.cap = []
        self.equal = []
        self.atleast = []

    def variable(self, name, cap=1):
        col = len(self.names)
        self.names.append(name)
        self.cap.append(cap)
        return col

    def eq(self, name, expr, rhs=0):
        self.equal.append((name, expr, rhs))

    def ge(self, name, expr, rhs=0):
        self.atleast.append((name, expr, rhs))


class Occupation(Problem):
    """Tail colours interpolate reciprocal degree. Drift tests see their union."""

    def __init__(self, cutoff, shared=True, chain=True, low=LOW, high=HIGH):
        super().__init__()
        if cutoff < 3 or not 0 < low <= high < 1:
            raise ValueError("invalid cutoff or growth interval")
        self.cutoff, self.shared, self.chain = cutoff, shared, chain
        self.low, self.high = low, high
        self.finite = [(d, c) for d in range(2, cutoff + 1) for c in range(1, d + 1)]
        self.labels = [*self.finite, ZERO, ONE]
        self.inner_parents = [s for s in self.finite if s[1] >= 2]
        self.parents = [*self.inner_parents, ZERO, ONE]
        self.inv = {s: F(1, s[0]) for s in self.finite} | {ZERO: 0, ONE: F(1, cutoff+1)}
        self.v = {s: self.variable(('V', s)) for s in self.labels}
        self.e = {(a, b): self.variable(('E', a, b)) for a in self.parents for b in self.labels}
        self.t = {(a, s, b): self.variable(('T', a, s, b))
                  for s in self.inner_parents for a in self.parents for b in self.labels}
        self.ch = {(s, a, b): self.variable(('C', s, a, b), cutoff)
                   for s in self.inner_parents if s[1] >= 3
                   for a, b in combinations_with_replacement(self.labels, 2)}
        self.q = {(a, m, BARE, w): self.variable(('Q', a, m, BARE, w))
                  for m in self.inner_parents for a in self.parents for w in self.labels}
        self.mass = {}
        self.flow = {}
        self.extra_flux = {}
        self.lift = None
        self.channel = None
        self._marginals()
        self._balance()

    def _marginals(self):
        v, e, t = self.v, self.e, self.t
        self.eq('unit', {col: 1 for col in v.values()}, 1)
        for b in self.labels:
            row = {e[a, b]: 1 for a in self.parents}
            row[v[b]] = -1
            self.eq(('child', b), row)
        for s in self.inner_parents:
            k = s[1] - 1
            row = {e[s, b]: 1 for b in self.labels}
            row[v[s]] = -k
            self.eq(('parent', s), row)
            for a in self.parents:
                row = {t[a, s, b]: 1 for b in self.labels}
                row[e[a, s]] = -k
                self.eq(('path-left', a, s), row)
            for b in self.labels:
                row = {t[a, s, b]: 1 for a in self.parents}
                row[e[s, b]] = -1
                self.eq(('path-right', s, b), row)
        siblings = defaultdict(dict)
        for (s, a, b), j in self.ch.items():
            add(siblings[s, a], j, 1)
            add(siblings[s, b], j, 1)
        for (s, a), row in siblings.items():
            row[e[s, a]] = -(s[1] - 2)
            self.eq(('sibling', s, a), row)
        for m in self.inner_parents:
            for a in self.parents:
                row = {self.q[a, m, BARE, w]: 1 for w in self.labels}
                row[t[a, m, BARE]] = -1
                self.eq(('four-left', a, m), row)
            if self.shared:
                for w in self.labels:
                    row = {self.q[a, m, BARE, w]: 1 for a in self.parents}
                    row[t[m, BARE, w]] = -1
                    self.eq(('four-right', m, w), row)

    def gains_to(self, s):
        d, c = s
        return (d+1, c) if d < self.cutoff else TAIL

    def _balance(self):
        v, e, r = self.v, self.e, self.inv
        mass, flow = defaultdict(dict), defaultdict(dict)
        parent_gain, child_gain = defaultdict(dict), defaultdict(dict)
        vg = {s: {v[s]: s[0]-s[1]} for s in self.finite}
        self.birth = {s: {v[s]: F(s[0]-s[1], s[0])} for s in self.finite}
        # c/D = (one incoming slot + k outgoing slots)/D for a noninitial core.
        tail_birth = {v[ZERO]: 1, v[ONE]: 1-r[ONE]}
        for w in self.labels:
            add(tail_birth, e[ONE, w], -r[ONE])
        self.birth[TAIL] = tail_birth
        self.ge('tail-birth-nonnegative', tail_birth)
        self.nu = combine(*[(1, b) for b in self.birth.values()])
        self.ge('growth-low', self.nu, self.low)
        self.ge('growth-high', combine((-1, self.nu)), -self.high)

        for s, j in v.items():
            mass['V', coarse(s)][j] = 1
        for (a, b), j in e.items():
            p, w = coarse(a), coarse(b)
            mass['E', p, w][j] = 1
            if a[0]:
                add(vg[a], j, r[b])
                add(parent_gain[p, w], j, a[0]-a[1]+r[b])
            if b[0]:
                add(vg[b], j, r[a])
                add(child_gain[p, w], j, b[0]-b[1]+r[a])
        for (a, s, b), j in self.t.items():
            add(parent_gain[s, coarse(b)], j, r[a])
            add(child_gain[coarse(a), s], j, r[b])
        for (s, a, b), j in self.ch.items():
            add(parent_gain[s, coarse(a)], j, r[b])
            add(parent_gain[s, coarse(b)], j, r[a])

        def move(old, new, flux):
            axpy(flow[old], flux, -1)
            axpy(flow[new], flux)

        for s in self.finite:
            d, c = s
            move(('V', s), ('V', self.gains_to(s)), vg[s])
            if c < d:
                move(('V', s), ('V', (d, c+1)), self.birth[s])
                axpy(flow['E', (d, c+1), TIP], self.birth[s])
        axpy(flow['V', TIP], self.nu)
        axpy(flow['E', TAIL, TIP], tail_birth)
        for key, occ in list(mass.items()):
            if key[0] != 'E':
                continue
            _, p, w = key
            if p != TAIL:
                d, c = p
                move(key, ('E', self.gains_to(p), w), parent_gain[p, w])
                if c < d:
                    move(key, ('E', (d, c+1), w), combine((F(d-c, d), occ)))
            if w != TAIL:
                d, c = w
                move(key, ('E', p, self.gains_to(w)), child_gain[p, w])
                if c < d:
                    move(key, ('E', p, (d, c+1)), combine((F(d-c, d), occ)))
        self.vertex_gain = vg
        self.parent_gain = parent_gain
        if self.chain:
            self._chain_balances(mass, flow, parent_gain, child_gain)
        self.mass, self.flow = dict(mass), dict(flow)
        for key, occ in mass.items():
            self.ge(('balance-low', key), combine((1, flow[key]), (-self.low, occ)))
            self.ge(('balance-high', key), combine((-1, flow[key]), (self.high, occ)))

    def _chain_balances(self, mass, flow, parent_gain, child_gain):
        # H_n is an oriented path from a top through n bare cores to an endpoint.
        # Both tests use the union of the two tail colours at either end.
        halves = F(1, 2)
        labels = [*self.finite, TAIL]
        parents = [*self.inner_parents, TAIL]
        ends = {}
        for n in (1, 2):
            for a in self.parents:
                for w in self.labels:
                    key = ('chain', coarse(a), n, coarse(w))
                    col = self.t[a, BARE, w] if n == 1 else self.q[a, BARE, BARE, w]
                    add(mass[key], col, 1)
                    # Each existing bare disappears upon a gain. Its two
                    # internal donations sum to r_top + (n-1) + r_end in total.
                    add(flow[key], col, -(self.inv[a]+n-1+self.inv[w]))
                    if w == TIP:
                        old = self.e[a, TIP] if n == 1 else self.t[a, BARE, TIP]
                        add(flow[key], old, halves)
            for a in parents:
                for w in labels:
                    key = ('chain', a, n, w)
                    occ = mass[key]
                    for site, s in (('top', a), ('end', w)):
                        if s == TAIL:
                            continue
                        d, c = s
                        extra = {}
                        if site == 'end' and n == 1 and w == BARE and self.shared:
                            # The known fourth site is the unique outside child.
                            for acol in ((ZERO, ONE) if a == TAIL else (a,)):
                                for z in self.labels:
                                    add(extra, self.q[acol, BARE, BARE, z], self.inv[z])
                        elif c > 1:
                            bound = (c-1)*halves
                            j = self.variable(('outside', site, a, n, w), bound)
                            self.extra_flux[site, a, n, w] = j
                            extra = {j: 1}
                            self.ge(('outside-cap', site, a, n, w), combine((bound, occ), (-1, extra)))
                        ends[site, a, n, w] = extra
                        gain = combine((d-c+halves, occ), (1, extra))
                        def move(to, rate):
                            new = ('chain', to, n, w) if site == 'top' else ('chain', a, n, to)
                            axpy(flow[key], rate, -1)
                            axpy(flow[new], rate)
                        move(self.gains_to(s), gain)
                        if c < d:
                            move((d, c+1), combine((F(d-c, d), occ)))
        for a in self.inner_parents:
            d, c = a
            edge = self.e[a, BARE]
            prev = combine((1, parent_gain[a, BARE]), (-(d-c+halves), {edge: 1}))
            for n in (1, 2):
                row = combine((-1, prev), *[(1, ends['top', a, n, w]) for w in labels])
                self.eq(('top-flux-overlap', a, n), row)
                prev = ends['top', a, 1, BARE]
        for w in self.finite:
            d, c = w
            edge = self.e[BARE, w]
            prev = combine((1, child_gain[BARE, w]), (-(d-c+halves), {edge: 1}))
            for n in (1, 2):
                row = combine((-1, prev), *[(1, ends['end', a, n, w]) for a in parents])
                self.eq(('end-flux-overlap', w, n), row)
                prev = ends['end', BARE, 1, w]

        # Birth of the fourth site from a pre-existing tip. For a general
        # middle hub additional entries are favourable. Bound every exit of the
        # already present configuration, aggregating the persistent top tail.
        for m in self.inner_parents:
            if m == BARE:
                continue
            dm, cm = m
            for p in [*self.inner_parents, TAIL]:
                labels = (ZERO, ONE) if p == TAIL else (p,)
                row = {}
                for a in labels:
                    born = self.t[a, m, TIP]
                    curr = self.q[a, m, BARE, TIP]
                    add(row, born, -halves)
                    # Known internal donors are a, m, the bare and the tip.
                    exit_rate = (dm-cm+self.inv[a]+halves+(cm-2)*halves
                                 + F(dm-cm, dm) + F(1, dm)+halves + 2)
                    if p != TAIL:
                        da, ca = p
                        exit_rate += da-ca+F(1, dm)+(ca-1)*halves+F(da-ca, da)
                    add(row, curr, self.high+exit_rate)
                self.ge(('four-entry', p, m), row)
            # The same birth balance one level earlier forces tip occupancy
            # by the postpromotion middle type, using its actual incoming edge.
            before = (dm, cm-1)
            for p in [*self.inner_parents, TAIL]:
                labels = (ZERO, ONE) if p == TAIL else (p,)
                row = {}
                for a in labels:
                    add(row, self.e[a, before], -F(dm-cm+1, dm))
                    exit_rate = dm-cm+self.inv[a]+halves+(cm-2)*halves+F(dm-cm, dm)
                    exit_rate += 1 + F(1, dm) + halves
                    if p != TAIL:
                        da, ca = p
                        exit_rate += da-ca+F(1, dm)+(ca-1)*halves+F(da-ca, da)
                    add(row, self.t[a, m, TIP], self.high+exit_rate)
                self.ge(('three-entry', p, m), row)

    def tip_hub(self):
        return {j: 1 for (a, w), j in self.e.items() if w == TIP and a != BARE}

    def tip_input(self):
        return {j: self.inv[a] for (a, w), j in self.e.items() if w == TIP}
