#!/usr/bin/env python3
"""Search floating LP duals, then certify them with rational reduced costs."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'renewal'))
from check import BASE, BANDS, LO, HI, row_rewards, tracked

UNIT = 10**6
PRICE_LIMIT = 64
LIMIT = 36
SAT = 12
START = (2, 0)


class LP:
    def __init__(self):
        self.rows = {}
        self.columns = []
        self.rhs = {}

    def col(self, kind, values, meta):
        for row in values:
            self.rows.setdefault(row, len(self.rows))
        self.columns.append((kind, {r: Q(v) for r, v in values.items() if v}, meta))

    def bound(self, objective):
        from scipy.optimize import linprog
        from scipy.sparse import coo_matrix
        costs = [Q(objective(kind, meta)) for kind, _, meta in self.columns]
        ii, jj, vv = [], [], []
        indexed = []
        for j, (kind, col, _) in enumerate(self.columns):
            terms = [(self.rows[x], v) for x, v in col.items()]
            indexed.append((kind, terms))
            for i, v in terms:
                ii.append(i)
                jj.append(j)
                vv.append(float(v))
        m = len(self.rows)
        A = coo_matrix((vv, (ii, jj)), shape=(m, len(costs))).tocsc()
        rhs = [Q(0)]*m
        for row, v in self.rhs.items():
            rhs[self.rows[row]] = Q(v)
        opt = linprog(list(map(float, costs)), A_eq=A, b_eq=list(map(float, rhs)),
                      method='highs-ipm', options={'presolve': True, 'threads': 1})
        assert opt.success, opt.message
        grid = 10**11
        y = [Q(round(v*grid), grid) for v in opt.eqlin.marginals]
        # Repair any floating feasibility error using *exact* normalization
        # rows: sum occupation = 1, sum edge mass = 1, excess stubs <= 1.
        max_defect = defaultdict(Q)
        for c, (kind, terms) in zip(costs, indexed):
            error = sum((v*y[i] for i, v in terms), Q(0))-c
            max_defect[kind] = max(max_defect[kind], error)
        assert all(sum((v*y[i] for i,v in terms),Q(0))-max_defect[kind]<=c
                   for c,(kind,terms) in zip(costs,indexed))
        margin = sum(max_defect.values(), Q(0))
        result = sum((v*w for v, w in zip(rhs, y)), Q(0))-margin
        return result, margin


def band(d, k, detailed):
    # Same endpoint cells on both sides; retain leafy tips, intermediate,
    # and many-arm cores. LIMIT+1 represents the unbounded degree tail.
    ends = (*[d for d in (2, 3, 4, 5, 6, 8, 11, 14, 17, 21, 27, 36, 48, 64) if d<LIMIT],LIMIT,LIMIT+1)
    lo = 2
    for hi in ends:
        if d <= hi:
            return lo, hi, (min(k, 1) if not detailed else 0 if k == 0 else 1 if k < 4 else 4)
        lo = hi+1
    raise AssertionError(d)


def recips(cell):
    lo, hi, _ = cell
    return sorted({Q(1, lo), Q(0) if lo > LIMIT else Q(1, hi)})


def model(flow, detailed=True, lo=LO, hi=HI):
    lp = LP()
    cells, senders = set(), set()
    states = [(d, k) for d in range(2, LIMIT+1) for k in range(d)]
    states += [(LIMIT+1,k) for k in range(SAT+1)]
    lp.rhs = {'mass': 1, 'k': 1}
    if flow:
        lp.rhs['edges'] = 1
    for d, k in states:
        x = (d, k)
        cell = band(d, k, detailed)
        cells.add(cell)
        if k:
            senders.add(cell)
        tail = d > LIMIT
        if not tail:
            physical = [(Q(d-k-1, d), Q(1,d))]
        elif k < SAT:
            physical = [(1-(k+1)*p, p) for p in (Q(0), Q(1,d))]
        else:
            # Exact b+p+kp=1; relaxation omits a lower coupling to the
            # excess unweighted stubs of this aggregated row. Its upper
            # bound b+6p<=1 still follows from k_actual >= 5.
            physical = [(Q(0),Q(0)), (Q(1),Q(0)),
                        (Q(0),Q(1,d)), (1-Q(SAT+1,d),Q(1,d))]
        for s, (b,p), r in product((lo,hi), physical, (Q(0),Q(k+1,2)) if x!=(LIMIT+1,SAT) else (Q(0),)):
            col = defaultdict(Q, {'mass': 1, 'k': k, 'b': b-s})
            col['rho'] += r+s-1
            col['age', x] += s
            col['age', START] -= s
            if not tail:
                g = d-k-1+r
                col['age', x] += g+b
                col['age', (d+1,min(k,SAT) if d==LIMIT else k)] -= g
                if b:
                    col['age', (d,k+1)] -= b
            elif k<SAT:
                col['age', x] += b
                col['age', (d,k+1)] -= b
            if flow:
                col['child_n',cell] += 1
                col['child_p',cell] += p
                if k:
                    col['parent_n',cell] += k
                    col['parent_p',cell] += 1-b-p
                col['recv',cell] += r
            lp.col('occ', col, (x,b,r))
    # Above the census cutoff and k>=5, keep the full unknown additional
    # number of children and full incoming rate, not just those of k=5.
    high = band(LIMIT+1, SAT, detailed)
    col = {'k':1}
    if flow:
        col['parent_n',high] = 1
    lp.col('stub', col, None)
    col = {'rho': 1}
    if flow:
        col['recv',high] = 1
    lp.col('tailin', col, None)
    if flow:
        for a, b in product(sorted(senders), sorted(cells)):
            for p, q in product(recips(a), recips(b)):
                col = defaultdict(Q, {'edges':1, ('parent_n',a):-1,
                                      ('parent_p',a):-p, ('child_n',b):-1,
                                      ('child_p',b):-q})
                col['recv',a] -= q
                col['recv',b] -= p
                lp.col('edge', col, (a,b,p,q))
    return lp


def tail_price(pay, tail):
    # Price each remaining degree using the exact D<=64 rows and the
    # unbounded bands thereafter.  Receiver dependence beyond the census
    # seam is relaxed, not discarded.
    return tuple(min([tail[k] if k < len(tail) else Q(0)] +
                     [min(pay[d,k]) for d in range(LIMIT+1, PRICE_LIMIT+1) if k<d])
                 for k in range(SAT+1))


def pay_objective(pay, tail):
    def obj(kind, meta):
        if kind != 'occ':
            return 0
        (d,k), _, r = meta
        if d > LIMIT:
            return UNIT*tail[k] if k < len(tail) else 0
        u,v = pay[d,k]
        return UNIT*(u+(v-u)*r/Q(k+1,2))
    return obj


def outward(x, up=False, digits=12):
    n = 10**digits
    v = -(-x*n//1) if up else x*n//1
    sign = '-' if v<0 else ''
    a,b = divmod(abs(v),n)
    return f'{sign}{a}.{b:0{digits}d}'


def controls():
    # Parent with degree 8 and six core children, three each of degree 2
    # and 4. Reverse traffic has the same six incidences, not a new sample.
    p = Q(1,8)
    joint = {2: 3, 4: 3}
    up = sum(n*Q(1,d) for d,n in joint.items())
    down = sum(joint.values())*p
    assert (up,down) == (Q(9,4),Q(3,4))
    assert 6*Q(1,2)>up  # taking all upward rates at their pointwise cap
    # These directed incidence margins both have six entries, and match
    # total parent reciprocal traffic. But their child-degree columns are
    # inconsistent and cannot represent one table of physical edges.
    fake = {2:6,4:0}
    assert sum(fake.values())==sum(joint.values())
    assert sum(fake.values())*p==down and fake!=joint
    # A genuine common table also carries exact receiver mass-transport.
    # Tree of seven cores: top + six children, add d-1 leaves to
    # each child, two leaves at the free top to make its degree eight.
    degrees = [8,2,2,2,4,4,4]
    directed = [(0,i) for i in range(1,7)] + [(i,0) for i in range(1,7)]
    rho = [sum((Q(1,degrees[i]) for i,j in directed if j==v),Q(0)) for v in range(7)]
    assert sum(rho) == sum(Q(sum(i==v for i,_ in directed),degrees[v]) for v in range(7))
    for threshold in (2,4,8):
        assert sum(r for v,r in enumerate(rho) if degrees[v]>=threshold)==sum(
            Q(1,degrees[i]) for i,j in directed if degrees[j]>=threshold)
    assert sum(rho[v] for v in range(7) if degrees[v]>=8)!=Q(6,8)
    # Paired opposing traffic for a high-degree, many-arm receiver.
    # D=32,k=15 with an incoming parent of degree 16; children of degree
    # 2/8 occur 5/10 times. Opposite traffic can neither fill 15 virgin
    # edges nor use 15 degree-eight edges independently.
    assert Q(5,2)+Q(10,8)+Q(1,16)==Q(61,16)
    assert Q(15,2)+Q(1,16)>Q(61,16)
    # Stronger: both permutations separately satisfy *all* endpoint mass
    # and weighted-stub margins. Reusing the degree-two receiving/donating
    # class on different parent edges makes their separate maxima invalid.
    a, b = [(16,2),(32,8)], [(16,8),(32,2)]
    def ports(edges):
        return (sum(Q(1,y) for x,y in edges if x==32),
                sum(Q(1,x) for x,y in edges if y==2))
    assert ports(a)==(Q(1,8),Q(1,16)) and ports(b)==(Q(1,2),Q(1,32))
    assert max(sum(ports(e)) for e in (a,b))==Q(17,32)<Q(9,16)
    print('PASS: shared endpoint counts, reversed receiver weights, and independent-direction mutant 9/16 > 17/32')


def main():
    global LIMIT
    import warnings
    warnings.filterwarnings('ignore', message='Unrecognized options detected:.*')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse',action='store_true')
    parser.add_argument('--diagnostic',action='store_true',help='price receive flux into moderate/high multi-arm exit classes')
    parser.add_argument('--limit',type=int,default=LIMIT,choices=(18,24,28,36,48,64))
    args = parser.parse_args()
    LIMIT = args.limit
    controls()
    dump = json.loads((HERE.parent/'rows'/'charged-f.json').read_text())
    assert dump['payload_sha256']=={
        'lower':'9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
        'upper':'1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41'}
    models = [model(yes,not args.coarse) for yes in (False,True)]
    for m in models:
        print('LP',len(m.rows),'rows',len(m.columns),'columns',flush=True)
    denom = {'lower':Q('1.356833021311'), 'upper':Q('1.367169995210')}
    # Conservative floors obtained from size-12 published directed
    # endpoints and the same denominators used in that replay.
    higher = {
        'lower': (Q('.565675597658')-Q('.565674722414'))*Q('1.458909469275164'),
        'upper': (Q('.565773666188')-Q('.565772677869'))*Q('1.518769126087491'),
    }
    for side in ('lower','upper'):
        pay = row_rewards(dump,side)
        tail = tuple(map(Q,BANDS[side]))
        r0 = UNIT*tracked(pay,tail)
        results = [m.bound(pay_objective(pay,tail_price(pay,tail))) for m in models]
        print(side,'old census',outward(r0), 'normalized census, joint edge',
              *(outward(v) for v,_ in results),'exact repair',*(outward(e,True) for _,e in results),flush=True)
        price = max(r0,*(v for v,_ in results))/UNIT
        sign = 1 if side=='lower' else -1
        print(' first-body endpoint',outward(BASE[side]+sign*price/denom[side],side=='upper'),flush=True)
        # The floor is inferred from directed endpoints, not from the
        # approximate displays of the higher-order rewards.
        reward = higher[side]//Q('1e-11')*Q('1e-11')
        print(' with disjoint higher-order floor',outward(BASE[side]+sign*(price+reward)/denom[side],side=='upper'),flush=True)
    if args.diagnostic:
        for d0,k0 in ((18,4),(LIMIT,4),(LIMIT,8)):
            def reward(kind, meta):
                if kind == 'occ':
                    (d,k),_,r = meta
                    if d0<=d<=LIMIT and k>=k0:
                        return -r
                return 0
            results = [m.bound(reward) for m in models]
            print(f'pi(rho; {d0}<=D<={LIMIT},k>={k0}) upper',*(outward(-x,True) for x,_ in results),flush=True)


if __name__=='__main__':
    main()
