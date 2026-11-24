"""Signed bounds for the analytic core row, with physical gain-rate corners."""

from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import factorial
from pathlib import Path
from runpy import run_path
from subprocess import check_output

raw = run_path(str(Path(__file__).resolve().parents[1] / 'check.py'))
PIN, TRIAL, WORK, Payload, ledger, physical_core, physical_tree, geom_weight, up, subsets = (
    raw[k] for k in ('PIN', 'TRIAL', 'WORK', 'Payload', 'ledger', 'physical_core',
                     'physical_tree', 'geom_weight', 'up', 'subsets'))

LO = Q(2825181729593748861, 4994114995300000000)
HI = Q(439124781075521909, 776183857175000000)
S0 = Q(565702, 10**6)
M, K, C = 2000, 106, 55
GRID = 10**17


def fmt(q):
    return f'{float(q):+.12e}'


class Row:
    def __init__(self, ps):
        self.ps = ps
        p, u = ps
        S = p.h['S']
        assert p.h['mode'] == 2 and u.h['mode'] == 1
        assert all(x.h['DF'] == M and x.h['S'] == S for x in ps)
        self.a = Q(p.h['f0']-p.h['f1']+u.h['f0']+u.h['f1'], S)
        self.k = Q(p.h['f0']-2*p.h['f1']+u.h['f1'], S)
        self.h = Q(p.h['f2'], S)
        self.j = Q(p.h['f1'], S)
        self.gcap = max(q.ends[2].k+1 for q in ps)
        # At a known high D, the affine interpolation has two unknowns; two
        # sparse reads (knot and twice knot) give the exact INF and slope.
        def newborn(D, c):
            marks = {0: (D, c), 1: (2, 1)}
            return sum(Q(q.value('f', ((2, 1),), ()) + q.row('g', (), (0, 1), marks), S) for q in ps)
        d = p.ends[2].d
        assert u.ends[2].d == d
        self.n = 2*newborn(2*d, 1)-newborn(d, 1)
        self.z = tuple(d*(newborn(d, min(c+1, self.gcap))-self.n) for c in range(K+1))
        assert self.k > 0 and self.h > 0 and self.j > 0
        self.scale = S
        self.ticks = tuple(int(S*v) for v in (self.a, self.k, self.h, self.j, self.n))
        self.z_ticks = tuple(int(S*v) for v in self.z)
        assert all(Q(n, S) == v for n, v in zip(self.ticks+self.z_ticks, (self.a, self.k, self.h, self.j, self.n)+self.z))

    def coeff(self, c):
        return -self.k*c-self.h, -self.j*c*c+self.h*c

    def value(self, D, c):
        a1, a2 = self.coeff(c)
        return self.a + a1/D + a2/D**2

    def linear(self, c):
        return self.z[c] - self.k + (1+TRIAL)*self.h + c*((1+TRIAL)*self.k-self.n)

    def base(self):
        return self.n-TRIAL*self.a

    def terms(self, D, c):
        """Return non-gain reward, exact Δf and polynomial Δf from lane C."""
        a1, a2 = self.coeff(c)
        x = Q(1, D)
        f = self.a + a1*x + a2*x*x
        df = -a1*x*x/(1+x)-a2*x**3*(2+x)/(1+x)**2
        poly = -a1*x*x*(1-x+x*x)-a2*x**3*(2+x)*(1-2*x+3*x*x)
        promo = self.n+(self.z[c]-self.k)*x+(-self.j*(2*c+1)+self.h)*x*x
        return (1-c*x)*promo-TRIAL*f, df, poly

    def quadratic_envelope(self):
        """|raw - (base + linear(c)/D)| <= (A+B c+H c²)/D² for c<=K."""
        z = max(abs(v-self.k) for v in self.z)
        # Bound the gain error after extracting -a1/D, then promotion and discount.
        A = self.h + self.j + self.h
        B = self.k + self.h + 3*self.h + 2*self.j + z + TRIAL*self.h
        H = self.k + (3+TRIAL)*self.j
        # b*promotion coefficient is multiplied by b<=1, so no cubic is needed.
        return A, B, H

    def poly_reward_ticks(self, D, c):
        """2*rate corner input; reward numerator is in GRID units."""
        a, k, h, j, n = self.ticks
        u, v = -k*c-h, -j*c*c+h*c
        # D^6 clears the polynomial derivative: Δ(1/D)^2 term
        # contains (2+1/D)(1−2/D+3/D²)/D³.
        d2, d3 = D*D, D**3
        d6 = d3*d3
        shift = -u*(d2-D+1)*d2-v*(2*D+1)*(d2-2*D+3)
        pr = n*d2+(self.z_ticks[c]-k)*D-j*(2*c+1)+h
        disc = a*d2+u*D+v
        sden = TRIAL.denominator
        base = 2*((D-c)*pr*d3*sden-TRIAL.numerator*disc*d2*d2)
        den = 2*self.scale*d6*sden
        # promotion term needs D^3, discount D^2, derivative D^6.
        return tuple((GRID*(base+gr*shift*sden), den) for gr in (2*(D-c), 2*D-c))


@cache
def common_age(D, order=2):
    """Upper E[(N)_i; Y>=D] for independent proposals at one shared age."""
    p = Q(1)
    h = [up(S0**(-i)) for i in range(order+1)]
    for d in range(2, D):
        rate = d-Q(1, 2)
        p = up(p*rate/(rate+S0))
        for i in range(1, order+1):
            h[i] = up(h[i]+h[i-1]/(rate+S0))
    return tuple(p*factorial(i)*v for i, v in enumerate(h))


def exterior(row, D):
    """Expected range for the signed unbounded-degree suffix on c<=C."""
    w0, w1, w2 = common_age(D)
    # c=1+N envelope; moments come from that same age and degree hit.
    a, b, h = row.quadratic_envelope()
    e = (a*w0+b*(w0+w1)+h*(w0+3*w1+w2))/D**2
    lo = min(Q(0), min(row.linear(c) for c in range(1, C+1)))
    # Positive linear coefficients occur only in the rare higher arities.
    hits = [c for c in range(1, C+1) if row.linear(c) > 0]
    positive = max((row.linear(c) for c in hits), default=Q(0))
    pos_mass = min(w0, geom_weight(min(hits)-1, 0)) if hits else Q(0)
    return (min(Q(0), row.base()+lo/D)*w0-e,
            max(Q(0), row.base())*w0+positive*pos_mass/D+e)


def overflow(row):
    """Common target excludes c>K only; this pays the computation's c>C gap."""
    a, b, h = row.quadratic_envelope()
    l = max(abs(row.linear(c)) for c in range(C+1, K+1))
    return ((abs(row.base())+l/M+a/M**2)*geom_weight(C, 0)
            +b/M**2*geom_weight(C, 1)+h/M**2*geom_weight(C, 2))


def solve(row, stop):
    """Integer-grid Bellman certificates on the finite part, same target both signs."""
    # Use S0 and a decimal upper rational, with product checked separately.
    low, high = S0, Q(565749, 10**6)
    assert low <= LO <= HI <= high
    den = 10**6
    sl, sh = int(low*den), int(high*den)
    nxt = [[0]*(C+2), [0]*(C+2)]
    for D in range(stop-1, 1, -1):
        cur = [[0]*(C+2), [0]*(C+2)]
        for c in range(min(C, D), 0, -1):
            if D >= M:
                # Rate-weighted reward is affine in the very same gain used
                # as the transition. Represent rewards on exact integer scale.
                rn = row.poly_reward_ticks(D, c)
            else:
                rn = [(0, 1), (0, 1)]
            promo = 2*(D-c)
            for side in (0, 1):
                bounds = []
                for gr, (r, rd) in zip((2*(D-c), 2*D-c), rn):
                    numerator = (gr*D*nxt[side][c] + promo*cur[side][c+1])*rd + 2*D*r
                    n = numerator*den
                    d = (gr*D+promo)*den
                    for s in (sl, sh):
                        divisor = (d+2*D*s)*rd
                        bounds.append((n, divisor))
                v = (min(n//d for n, d in bounds) if side == 0 else
                     max(-((-n)//d) for n, d in bounds))
                assert all(v*d <= n if side == 0 else v*d >= n for n, d in bounds)
                cur[side][c] = v
        nxt = cur
    raw = Q(nxt[0][1], GRID), Q(nxt[1][1], GRID)
    return min(LO*raw[0], HI*raw[0]), max(LO*raw[1], HI*raw[1])


def split(row, graph, poly):
    core, marks = physical_core(graph)
    out = defaultdict(Q)
    Dset = {v for v, (D, c) in marks.items() if D >= M and c <= K}
    for v in Dset:
        D, c = marks[v]
        b = Q(D-c, D)
        born = sum(Q(p.value('f', ((2, 1),), ()), p.h['S']) for p in row.ps)
        base, gain, taylor = row.terms(D, c)
        receive = sum(Q(1, marks[u][0]) for u in core[v] if marks[u][1] <= K)
        out['paid'] += base+(D-c+receive)*(taylor if poly else gain)
        # The exactly observed complement groups (not estimated by the bound).
        out['f_birth'] += born*b
        out['g_birth'] += b*(row.n-born+row.z[c]/D)
        out['f_discount'] -= TRIAL*row.value(D, c)
        out['f_promo'] += b*(row.value(D, c+1)-row.value(D, c))
        delta = taylor if poly else gain
        out['f_leaf'] += (D-c)*delta
        out['f_external'] += receive*delta
    assert out['paid'] == sum(v for k, v in out.items() if k != 'paid')
    return out


def direct_complement(row, graph):
    """Directed ordinary births with just the selected high-action changes removed."""
    old_core, old_marks = physical_core(graph)
    high = {u for u, (D, _) in old_marks.items() if D >= M}

    def weights():
        pairs = defaultdict(Q)
        for p in row.ps:
            for k, (v, _) in subsets(p, *physical_core(graph)).items():
                pairs[k] += Q(v, p.h['S'])
        return pairs

    pre = weights()
    result = -TRIAL*sum(v for (k, roles), v in pre.items() if k != 'f' or roles[0] not in high)
    rates = defaultdict(Q)
    for u, ns in enumerate(graph):
        for v in ns:
            rates[v] += Q(1, len(ns))
    for v, rate in rates.items():
        if v in old_core:
            excluded = {('f', (v,))} if v in high else set()
        else:
            u, = graph[v]
            excluded = {('f', (u,)), ('f', (v,)), ('g', tuple(sorted((u, v))))} if u in high else set()
        graph[v].add(len(graph))
        graph.append({v})
        post = weights()
        graph.pop()
        graph[v].remove(len(graph))
        result += rate*sum(post[k]-pre[k] for k in (pre.keys() | post.keys()) - excluded)
    return result


def fixtures(row):
    # Finite checks against raw ledger: split high and low promoters/receivers
    # without relying on per-endpoint charged singleton allocations.
    for D in (1999, 2000, 2001, 4096, 65536):
        for c in (1, 2, 12, C):
            base, _, df = row.terms(D, c)
            for (n, d), gr in zip(row.poly_reward_ticks(D, c), (D-c, D-Q(c,2))):
                assert Q(n, d) == GRID*(base+gr*df)
        for graph in (physical_tree([], [D]), physical_tree([(0, 1), (1, 2), (1, 3)], [D, 4, 2, 2])):
            for poly in (False, True):
                parts = split(row, graph, poly)
                raw = defaultdict(Q)
                for p in row.ps:
                    for (k, ch), v in ledger(p, graph, TRIAL, kcut=K,
                                              jcut=M-1 if poly else None,
                                              order=2 if poly else None).items():
                        raw[k+'_'+ch] += v
                # On one-core stars the chosen piece is the entire ledger.
                if len(physical_core(graph)[0]) == 1 and D >= M:
                    assert parts['paid'] == sum(raw.values())
                    assert parts['g_birth'] == raw['g_birth']
                if D < M:
                    assert not parts['paid']
                if not poly and D <= 2001:
                    assert direct_complement(row, graph) == sum(raw.values())-parts['paid']
            print(f'fixture D={D} cores={len(physical_core(graph)[0])}: paid={fmt(parts["paid"])} signed complement={fmt(sum(raw.values())-parts["paid"])}')
    # An overflow donor adjacent to a high receiver is excluded even though
    # the high receiver itself survives. The only remaining inflow is leaves.
    graph = physical_tree([(1, i) for i in range(K+2) if i != 1], [M, K+1]+[2]*K)
    assert physical_core(graph)[1][1] == (K+1, K+1)
    assert split(row, graph, False)['f_external'] == 0
    assert split(row, graph, True)['f_external'] == 0


def main():
    assert check_output(['git', '-C', str(WORK/'source'), 'rev-parse', 'HEAD'], text=True).strip() == PIN
    row = Row([Payload(WORK/f'quintet_{side}.bin') for side in ('lower', 'upper')])
    fixtures(row)
    stop = 10000
    a, b = solve(row, stop)
    lo, hi = exterior(row, stop)
    e = overflow(row)
    rom = common_age(M)
    z = (abs(row.coeff(K)[0])*rom[0]/M**4
         + 12*(row.j*K*K+row.h*K)*rom[0]/M**5)
    # The finite Bellman comparison uses complete receives; the retained
    # definition omits a send from c>K. Re-root at the actual donor. Its
    # total outgoing rate to cores is c/D<=1, so no arity multiplier remains.
    delta = abs(row.coeff(K)[0])/M**2 + 3*(row.j*K*K+row.h*K)/M**3
    donors = delta*geom_weight(K, 0)
    assert z < Q('1.947e-13') and donors < Q('8.277e-26')
    print(f'signed finite [{fmt(a)}, {fmt(b)}]')
    print(f'signed degree suffix [{fmt(lo)}, {fmt(hi)}]; arity envelope {fmt(e)}; raw-to-retained {fmt(z)}')
    print(f'overflow donor correction {fmt(donors)}')
    # The finite block is polynomial; the degree/arity suffixes above were
    # bounded for the exact block. The finite residual joins either target;
    # a second uniform error is safe for the retained variant.
    lower, upper = a+lo-e-2*z-donors, b+hi+e+2*z+donors
    print(f'common raw / retained [{fmt(lower)}, {fmt(upper)}], width {fmt(upper-lower)}')
    assert Q('-5.829193e-7') <= lower <= upper <= Q('-1.116081e-7')
    assert upper-lower < Q('4.714e-7')
    print('sparse payload words:', *(len(p.reads) for p in row.ps))


if __name__ == '__main__':
    main()
