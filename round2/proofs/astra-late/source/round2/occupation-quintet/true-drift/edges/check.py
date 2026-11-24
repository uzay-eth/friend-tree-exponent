"""Rational edge balance with explicit physical seam and sender ownership."""

from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
from subprocess import check_output

m = run_path(str(Path(__file__).resolve().parents[1] / 'moments/check.py'))
M, K, S0, LO, HI, s = (m[k] for k in ('M', 'K', 'S0', 'LO', 'HI', 'TRIAL'))
S = 2**40
UP = Q('0.565749')
GRID = 10**18
CAP = Q(350000000000, S)  # sum of original verified far-edge absolute caps, rounded up


def fmt(x):
    return f'{float(x):+.10e}'


def report(name, lo, hi):
    assert lo <= hi
    print(f'{name}: [{fmt(lo)}, {fmt(hi)}]  width {fmt(hi-lo)}')


class Pair:
    def __init__(self):
        self.ps = [m['Payload'](m['WORK'] / f'quintet_{side}.bin') for side in ('lower', 'upper')]
        self.e = self.ps[0].ends[2]
        other = self.ps[1].ends[2]
        assert (self.e.d, self.e.k, self.e.labels) == (other.d, other.k, other.labels)
        self.row = m['Row'](self.ps)

    def word(self, i, j):
        a, b = sorted((i, j))
        return sum(p.word(p.offset['g']+b*(b+1)//2+a) for p in self.ps)

    def f(self, x):
        return sum(Q(p.value('f', (x,), ()), S) for p in self.ps)

    def g(self, x, y):
        return sum(Q(p.value('g', (), (x, y)), S) for p in self.ps)

    def verify_faces(self):
        e = self.e
        t = [e.ids[e.d, c] for c in range(1, e.k+2)] + [e.inf]
        assert all(self.word(i, j) == 0 for i in t for j in t)
        # Full INF column only. This gives an exact interval when its partner
        # interpolates. Other tail rows need not be scanned.
        assert all(-CAP <= Q(self.word(e.inf, j), S) <= 0 for j in range(e.n))
        # A sparse set of occupied low-leaf classes, uniformly over ALL high c.
        limit = -Q(1, 100)
        for d in (2, 3, 4, 5):
            j = e.ids[d, 1]
            for i in t:
                w = Q((M-e.d)*self.word(e.inf, j)+e.d*self.word(i, j), M*S)
                assert max(w, Q(self.word(e.inf, j), S)) <= limit
        return -limit

    def directional(self):
        """Exact extrema of high-side coefficients and the reverse donor jump."""
        e = self.e
        slots = [e.ids[e.d, c] for c in range(1, e.k+2)]
        gain, promo, reverse, l_gain, l_promo = ([0, 0] for _ in range(5))
        def extend(r, x):
            r[0], r[1] = min(r[0], x), max(r[1], x)
        for j, (D,c) in enumerate(e.labels):
            inf = self.word(e.inf, j)
            for i, slot in enumerate(slots):
                edge = self.word(slot, j)
                extend(gain, e.d*(inf-edge))
                extend(promo, e.d*(self.word(slots[min(i+1,len(slots)-1)],j)-edge))
                if D < e.d:
                    j2 = e.ids[D+1,c]
                    delta_inf = self.word(e.inf,j2)-inf
                    delta_knot = self.word(slot,j2)-edge
                    for delta in (delta_inf, Q((M-e.d)*delta_inf+e.d*delta_knot,M)):
                        extend(reverse, delta)
                        extend(l_gain, delta*Q(2*D-c,2))
                    if c<D and c<=e.k:
                        j3 = e.ids[D,c+1]
                        db0 = self.word(e.inf,j3)-inf
                        db1 = self.word(slot,j3)-edge
                        for delta in (db0, Q((M-e.d)*db0+e.d*db1,M)):
                            extend(l_promo, delta*Q(D-c,D))
        return tuple(tuple(Q(x,S) for x in r) for r in (gain,promo,reverse,l_gain,l_promo))


def newborn_bounds(pair, top=34):
    """Bound E_pi[1_high b] using a monotone two-coordinate birth calculation.

    Values below M are conditional hitting/continuation rewards; the high
    suffix is a signed-safe constant bound using occupation and common age.
    """
    # At a hit with c<=top, future promotions of the same core may take c
    # above top. Pay a constant 1 on the comparison rather than absorbing here:
    # bound pi(1_high) and separately deduct E[1_high c/D].
    nxt = [[GRID]*(top+2), [GRID]*(top+2)]
    for d in range(M-1, 1, -1):
        cur = [[0]*(top+2), [GRID]*(top+2)]
        for c in range(min(d, top), 0, -1):
            b = 2*(d-c)
            for side in (0, 1):
                quotients = []
                for g in (2*(d-c), 2*d-c):
                    n = (g*d*nxt[side][c]+b*cur[side][c+1])*10**6
                    for discount in (565702, 565749):
                        den = (g*d+b)*10**6+2*d*discount
                        quotients.append((n, den))
                v = (min(n//z for n,z in quotients) if side == 0 else
                     max(-(-n//z) for n,z in quotients))
                assert all(v*z <= n if side == 0 else v*z >= n for n,z in quotients)
                cur[side][c] = v
        nxt = cur
    p0, p1, p2 = m['common_age'](M)
    lo, hi = Q(nxt[0][1], GRID), min(Q(nxt[1][1], GRID), p0)
    c1, c2 = p0+p1, p0+3*p1+p2
    assert 0 < lo <= hi
    return lo, hi, max(Q(0), lo-c1/M), hi, c1, c2


def seam_floor(low_degree_weight, mass, target):
    """Flux at a pre-seam core directly supporting a young actual child.

    From D=M-1 with c<=target, demand a promotion before its next gain;
    then demand that gain before the new child's first change. This is
    a subset of the seam, not of high promotions.
    """
    # All paths go from first entry into M-1. For the old root the desired
    # promotion hazard >= 1-target/(M-1); root gain <= M-1. New child
    # changes at rate <=2.5; root's desired gain >= M-2-target.
    b = 1-Q(target, M-1)
    ratio = b/(UP+M)
    g = M-2-target
    hit = ratio*g/(UP+g+Q(7,2))
    return -low_degree_weight*LO*mass*hit


def scalar(pair):
    wlow = pair.verify_faces()
    p_lo, p_hi, b_lo, b_hi, c1, c2 = newborn_bounds(pair)
    report('physical high-core mass', p_lo, p_hi)
    report('physical high-promotion intensity', b_lo, b_hi)

    age = m['geom_weight']
    # Birth edge N-F(J). Enclose each coefficient at D=M and infinity.
    born = pair.f((2,1))
    n = pair.row.n-born
    clo = min(n, *(n+z/M for z in pair.row.z))
    chi = max(n, *(n+z/M for z in pair.row.z))
    assert chi < 0
    birth = (clo*b_hi- CAP*age(K,0), chi*max(Q(0), b_lo-age(K,0)))
    report('high newborn g (already in H)', *birth)

    # Sum of weights of all old edges with at least one high endpoint.
    # Re-root at a high endpoint, once or twice (absolute bound only).
    wg = CAP*c1
    # Each admitted high promotion supplies an actual virgin edge. It remains
    # selected until the child changes or its high parent crosses arity K.
    # Union bounding those hazards by 5/2+1 is valid even when correlated.
    # The edge has weight <= chi throughout virgin residence (current parent
    # class differs from the entry class, but the same uniform bound works).
    # for births the actual parent's post-c is always at least two.
    # At a birth from c=K the successor is not selected. Remove that flux.
    w_floor = (-chi)*max(Q(0), b_lo-age(K-1,0))/(HI+Q(7,2))
    assert s > LO and s < HI
    stationary = (-(HI-s)*wg, (s-LO)*wg)  # g<=0
    report('selected high-edge potential mass', -wg, -w_floor)
    report('(mu-s) times selected g mass', *stationary)

    # At crossing, only neighbors below d_g can contribute (double tails
    # exactly zero). Without assuming this number from a marginal, use c.
    # Conditional proposal moments at D>=M-1, and crossing identity with
    # c-weight: on the first crossing c is dominated by the Poisson at
    # age A on {D>=M}; the normalized stationary flux is mu times that
    # discounted first-crossing expectation (pre C equals post C).
    target = 34
    mass = max(Q(0), p_lo-age(target-1,0))
    seam = (-CAP*HI*c1, seam_floor(wlow, mass, target))
    report('old-edge seam inflow (signed)', *seam)

    # Restore arity-excluded exterior core arrows into a selected endpoint.
    # At an overflow sender u, c/D<=1 and any retained receiver has <=K
    # incident edges, each |delta|<=CAP because high-face values are in [-CAP,0].
    # Loss of edge weight itself at an endpoint promotion crossing K is
    # a *different* exit; bounded by CAP*K*pi(c>=K).
    # The matched birth is paid at pre-c=K in H, even though the newborn
    # edge is excluded from the capped potential because its parent is now
    # at c=K+1. Both this flux and the old losses are covered here.
    overflow = CAP*(K*age(K,0)+(K+1)*age(K-1,0))
    old = (stationary[0]-birth[1]-seam[1]-overflow,
           stationary[1]-birth[0]-seam[0]+overflow)
    report('retained old high-edge subtotal', *old)

    # H and B are not priced independently when we add the edge sector.
    # H-B retains f offspring, discount, own gain and old promotion on
    # exactly the same high vertices. Use a near-constant holding reward
    # and their certified physical mass, with quadratic error common-age.
    r = pair.row
    a,b,h = r.quadratic_envelope()
    b += max(Q(0), r.k-max(abs(z-r.k) for z in r.z))
    e = (a*p_hi+b*c1+h*c2)/M**2
    # F(J)-s a leading constant; linear coefficient without any newborn g.
    l0 = -r.k+(1+s)*r.h
    l1 = (1+s)*r.k-born
    h0 = born-s*r.a
    low = h0*p_lo + min(Q(0), l0)*p_hi/M + min(Q(0), l1)*c1/M-e
    high = h0*p_hi + max(Q(0), l0)*p_hi/M + max(Q(0), l1)*c1/M+e
    # Arity-overflow norm at this mask (not same as old-edge overflow).
    loss = (abs(h0)+abs(l0)/M+a/M**2)*age(K,0) + (abs(l1)/M+b/M**2)*age(K,1)+h/M**2*age(K,2)
    gain_edge = (r.k*K+r.h)/M**2+3*(r.j*K*K+r.h*K)/M**3
    loss += gain_edge*age(K,0)
    # Intersect direct H-B with the independently bounded full H minus B.
    hb = (max(low-loss, Q('-5.829193e-7')-birth[1]),
          min(high+loss, Q('-1.116081e-7')-birth[0]))
    report('high H minus already-counted g birth', *hb)
    broad = (hb[0]+stationary[0]-seam[1]-overflow,
             hb[1]+stationary[1]-seam[0]+overflow)
    report('H plus old high edges (matched scalar)', *broad)

    ga, pr, rev, low_ga, low_pr = pair.directional()
    # Sum over actual incident edges (double-high edges have zero weight).
    # G<=D and promotion b<=1, hence one factor 1/D on these directions.
    report('old edge: gains at high endpoint, all senders', ga[0]*c1/M, ga[1]*c1/M)
    report('old edge: promotions at high endpoint', pr[0]*c1/M, pr[1]*c1/M)
    # A high->low arrow has one reciprocal of the actual high sender; a
    # low->high arrow pays the derivative at the high receiver, and its
    # total low input is bounded by c/2. Reroot the latter at the receiver.
    report('old edge: internal high->low send', rev[0]*c1/M, rev[1]*c1/M)
    report('old edge: low core sends into high receiver (internal + exterior)',
           ga[0]*c2/(2*M*M), ga[1]*c2/(2*M*M))
    report('old edge: gains at low endpoint, all senders', low_ga[0]*c1, low_ga[1]*c1)
    report('old edge: promotions at low endpoint', low_pr[0]*c1, low_pr[1]*c1)

    # Sender-to-low receiver, plus high donor appended to old low-low g.
    # Max raw jump budget on that receiver: singleton <=64 and low-low g
    # <=64 per incidence. Re-root at the high donor (all clock rates 1/D).
    # Exact virgin sender correction: there are no low-low old edges here.
    df = pair.f((3,1))-pair.f((2,1))
    assert df>0
    virgin = df*c1/M
    report('high sender into virgin low singleton (signed slice)', 0, virgin)
    send = (64+64*K)*c1/M
    report('remaining high sender into low f and old low-low g', -send, send)
    matched = broad[0], broad[1]+virgin
    report('broader matched scalar including virgin sends', *matched)

    # Independent total paired drift bounds from original all-tree comparison.
    # Remaining sector = full raw - exact capped broad scalar; omitted arity
    # rows and all higher-family terms are intentionally still inside it.
    weight_min = Q(99882299906+155236771435,S)
    weight_max = Q(1808571640916+1808494729937,S)
    total = ((LO-s)*weight_max, (HI-s)*weight_max)
    assert weight_min > 0
    report('full raw scalar (gold bracket only)', *total)
    low_sector = total[0]-matched[1], total[1]-matched[0]
    report('remaining signed sector incl. sender and routing', *low_sector)
    assert Q('9.1908e-4') <= matched[0] <= matched[1] <= Q('3.00631e-2')
    assert Q('-3.01349e-2') <= low_sector[0] <= low_sector[1] <= Q('-8.3879e-4')

    # Numerical comparison without treating dependent signed intervals as
    # independent estimates of the total. The final inequality is pointwise
    # between the physical expectations at their actual common mu.
    ratio = max(-total[0], total[1])/matched[0]
    assert ratio < Q('0.088')
    print('maximum relative imbalance |full| / matched:', fmt(ratio))


def fixture(pair):
    for D, E, F in ((1999, 3, 2),(2000, 3, 2),(2001, 3, 2),(2000, 2000, 2),(2000, 2, 2000)):
        graph = m['physical_tree']([(0,1),(1,2)],[D,E,F])
        core, marks = m['physical_core'](graph)
        assert all(c<=K for _,c in marks.values())

        def selected(graph):
            co, ma = m['physical_core'](graph)
            return sum(pair.g(ma[u], ma[v]) for u in co for v in co[u]
                       if u<v and max(ma[u][0],ma[v][0])>=M)

        base = selected(graph)
        old, birth, seam, send, direct, full_diff = (Q(0),)*6
        # Independent full f+g potential, used to test the low-sector partition.
        def fg(graph):
            co, ma = m['physical_core'](graph)
            return sum(pair.f(x) for x in ma.values()) + sum(pair.g(ma[u], ma[v])
                    for u in co for v in co[u] if u<v)
        whole = fg(graph)
        virgin = sum(Q(1,marks[u][0])*(pair.f((3,1))-pair.f((2,1)))
                     for v,x in marks.items() if x==(2,1)
                     for u in core[v] if marks[u][0]>=M)
        rates = defaultdict(Q)
        leaves = {}
        for ns in graph:
            for v in ns:
                if v not in core:
                    parent, = graph[v]
                    v = leaves.setdefault(parent, v)
                rates[v] += Q(1,len(ns))
        for v, rate in rates.items():
            Dv = len(graph[v])
            if v in core:
                for u in core[v]:
                    x,y = marks[v],marks[u]
                    if max(x[0],y[0])>=M:
                        old += rate*(pair.g((x[0]+1,x[1]),y)-pair.g(x,y))
                    elif x[0]==M-1:
                        seam += rate*pair.g((M,x[1]),y)
                for u in core[v]:
                    if marks[u][0]>=M and Dv<M:
                        x=marks[v]; xx=(x[0]+1,x[1]); inv=Q(1,marks[u][0])
                        send += inv*(pair.f(xx)-pair.f(x))
                        for w in core[v]:
                            if marks[w][0]<M:
                                send += inv*(pair.g(xx,marks[w])-pair.g(x,marks[w]))
            else:
                u,=graph[v]
                x=marks[u]; xx=(x[0],x[1]+1)
                for w in core[u]:
                    if max(x[0],marks[w][0])>=M:
                        old += rate*(pair.g(xx,marks[w])-pair.g(x,marks[w]))
                if x[0]>=M:
                    birth += rate*pair.g(xx,(2,1))
            graph[v].add(len(graph));graph.append({v})
            direct += rate*(selected(graph)-base)
            full_diff += rate*(fg(graph)-whole)
            graph.pop();graph[v].remove(len(graph))
        assert direct==old+birth+seam
        if D == M-1:
            assert seam<0 and not old and not birth
        if D >= M:
            assert old-s*base>0 and birth<0
        if (D,E,F)==(M,2,M):
            assert not virgin and send>0  # two distinct high senders into a bare low core
        H = m['split'](pair.row, graph, False)['paid']
        broad = H+old-s*base+virgin
        residual = full_diff-s*whole-broad
        assert residual == (full_diff-direct)-s*(whole-base)-(H-birth-seam)-virgin
        print('fixture', (D,E,F), 'old',fmt(old-s*base),'birth',fmt(birth),'seam',fmt(seam),'exterior send',fmt(send))


if __name__ == '__main__':
    assert check_output(['git','-C',str(m['WORK']/'source'),'rev-parse','HEAD'],text=True).strip()==m['PIN']
    p=Pair()
    fixture(p)
    scalar(p)
    print('sparse payload words:', *(len(x.reads) for x in p.ps))
