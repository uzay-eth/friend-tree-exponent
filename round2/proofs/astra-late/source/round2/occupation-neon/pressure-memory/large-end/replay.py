#!/usr/bin/env python3
"""Exact D18 endpoint probes and common-age stationary tail budgets."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
import sys

OCC = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(OCC / 'all-contrast'), str(OCC / 'higher' / 'older')]
from scan import counts, pairs, receives, up, label
from extract import HI, SPEC, embedding
from tail import bernstein
from replay import Local


S = Q(565673, 10**6)
GRID = 10**24
MAX_DEGREE = 262144
ROLES = (('P4', 0), ('P4', 3), ('Yb', 1), ('Yb', 3))
# Compatible occurrences when re-rooted at the end: factor times k**power.
MULT = {('P4', 0): (49, 1), ('P4', 3): (1, 0),
        ('Yb', 1): (42, 0), ('Yb', 3): (6, 0)}


def ceil(n, d):
    return -(-n // d)


def factorial_bounds(cutoffs, order):
    """Upper bounds on E[(N)_m 1{Y>J}] for m<=order, positive rounding."""
    p = GRID
    h = [GRID] + [ceil(GRID * S.denominator**m, S.numerator**m)
                  for m in range(1, order + 1)]
    out = {}
    for d in range(2, max(cutoffs) + 1):
        a = (2*d-1)*S.denominator
        b = a + 2*S.numerator
        p = ceil(p*a, b)
        # h_m bounds E[(T')^m]/m!: add Exp((d-1/2)+S).
        for m in range(1, order + 1):
            h[m] += ceil(h[m-1]*2*S.denominator, b)
        if d in cutoffs:
            out[d] = tuple(Q(ceil(p*v, GRID)*factorial(m), GRID)
                           for m, v in enumerate(h))
    return out


def power_moment(f, power):
    return (f[0], f[1], f[1]+f[2])[power]


def poisson_tail(k, power):
    """E[N**power 1{N>k}] at the same Exp(S) clock, geometric mixture."""
    r = 1/(1+S)
    f = (Q(1), 1/S, 2/S**2)
    return r**(k+1)*sum(comb(power,j)*(k+1)**(power-j)*power_moment(f,j)
                       for j in range(power+1))


def remainder(role, j, k, f, coefficient=None, *, exterior=False):
    n, power = MULT[role]
    # Rate 1/D, or external endpoint-gain slope times arity, respectively.
    if not exterior:
        return n*(Q(1) if coefficient is None else coefficient)*(
            power_moment(f,power)/(j+1) + poisson_tail(k,power)/max(18,k+2))
    power += 1
    d = max(18,k+2)
    return n*coefficient*(power_moment(f,power)/((j+1)*(j+2))
                          + poisson_tail(k,power)/(d*(d+1)))


def controls(a, b, role):
    """Global Lipschitz factor 18 sup |W'(t)| for the signed row contrast."""
    import numpy as np

    family, slot = role
    c = a.c
    name = {'P4': 'p4tab', 'Yb': 'ybtab'}[family]
    axis = {('P4',0): 0, ('P4',3): 3, ('Yb',1): 2, ('Yb',3): 3}[role]
    where = c.par if role == ('P4',0) else c.child
    u, v = getattr(c, name), getattr(b.c, name)
    assert u.shape == v.shape and c.S == b.c.S and c.mode == 2 and b.c.mode == 1
    # Integer difference of differences is safe given an independently
    # checked per-slice bound on each original entry before arithmetic.
    best = 0
    for k in range(int(role == ('P4',0)), c.K+2):
        indices = [where[('inf',)]] + [where[(f'T{i}',k,1)] for i in (1,2,3)]
        prev = None
        for i in indices:
            x, y = np.take(u,i,axis=axis), np.take(v,i,axis=axis)
            for z in (x,y):
                assert int(z.min()) > -2**60 and int(z.max()) < 2**60
            # Signed row contrast is +drift(lo) - (-drift(hi)).
            value = x+y
            if prev is not None:
                diff = value-prev
                best = max(best, int(diff.max()), -int(diff.min()))
            prev = value
    return Q(54*best,c.S)


def transfer_max(a, b, role):
    """Maximum paired connected-deletion coefficient over compatible controls."""
    c = a.c
    if role[0] == 'P4':
        domain = c.p4
        if role[1] == 0:
            def line(x,y): return c.par[x],c.hv[y],slice(None)
            left, right, name = domain, domain, 'hvtab'
        else:
            def line(x,y): return slice(None),c.hv[x],c.child[y]
            left, right, name = domain, domain, 'hvtab'
    else:
        left, right = [x for x in c.ybu if x[2]], c.ybw
        if role[1] == 1:
            def line(x,y): return c.par[x],c.hv[y],slice(None)
            name = 'hvtab'
        else:
            def line(x,y): return c.hc[x],[c.tp(z,y) for z in c.child]
            name = 'hctab'
    u,v = getattr(c,name),getattr(b.c,name)
    best = 0
    for x in left:
        for y in right:
            x2,y2 = (up(x),y) if role[1] in (0,1) else (x,up(y))
            i,j = line(x,y),line(x2,y2)
            z0,z1,z2,z3 = u[i],u[j],v[i],v[j]
            assert all(int(w.min()) > -2**60 and int(w.max()) < 2**60 for w in (z0,z1,z2,z3))
            diff = (z1-z0)+(z3-z2)
            best = max(best,int(diff.max()),-int(diff.min()))
    return Q(best,c.S)


def probe(a, b, role):
    f,i = role
    x = [(0,1,1),(1,1,1),(1,1,1),(1,0,1)] if f == 'P4' else [
        (1,2,1),(1,0,1),(1,1,1),(1,0,1)]
    if role == ('P4',0):
        # Stage hV/P4/Yt pays only when all non-apex P4 slots have D <= 5.
        x[1] = (4,1,1)
    # Add one off-motif child at the high end; its presence forces a fifth core.
    k = x[i][1]+1
    g, ids, donor, recv, parent = next(e for e in receives(f) if e[2] == i)
    assert donor == i and i not in ids
    locals_ = Local(a.c), Local(b.c)

    def eval_at(d):
        y = label(tuple(x),i,(d-k-1,k,1))
        n = counts(f,y)
        ext = tuple((Q(n[2*j],8),Q(n[2*j+1],3)) for j in range(4))
        flat = sum(ext,())
        result = []
        for r, local in zip((a,b),locals_):
            z, c = r.coefficients(f,y)
            rvalue = z+sum(p*q for p,q in zip(c,flat))
            assert rvalue == local.charge(f,y,ext)
            suffix = tuple(y[j] for j in ids)
            assert r.routed(g,suffix,recv,parent)
            d3 = r.value(g,label(suffix,recv,up(suffix[recv])))-r.value(g,suffix)
            # This older transferred gain is independent of the sender census.
            result.append((rvalue,c[2*i+1],d3*r.sign/r.c.S))
        return tuple(p-q for p,q in zip(*result))

    def den(t): return (18+t)**3

    ds = tuple(range(18,26))
    ts = [Q(18,d) for d in ds]
    values = [eval_at(d) for d in ds]
    polynomials = [bernstein([(t,den(t)*v) for t,v in zip(ts,col)])
                   for col in zip(*values)]
    for d in (27,43,113):
        t = Q(18,d)
        v = tuple(sum(comb(7,j)*t**j*(1-t)**(7-j)*p[j] for j in range(8))
                  /den(t) for p in polynomials)
        assert v == eval_at(d)
    assert all(p[2] == values[0][2] for p in values)
    assert polynomials[1][0] == polynomials[1][1] == 0
    y18,y19,y20 = values[:3]
    second = y20[0]-2*y19[0]+y18[0]
    print(f'{f}:{i} D18 paired row {float(y18[0]):+.10e}; '
          f'end child coefficient {float(y18[1]):+.10e}; '
          f'transfer coefficient {float(y18[2]):+.10e}; '
          f'row forward curvature {float(second):+.10e}',flush=True)
    final = label(tuple(x),i,(17-k,k,1))
    assert sum(final[i]) == 18 and sum(v[1] for v in final) == 4
    if f == 'P4':
        assert final[0][1] < 5 and max(map(sum,final[1:])) > 5
    floor(a,b,role,final)


def floor(a,b,role,x):
    """Lower physical occupation of one explicitly prescribed planted tree."""
    f, end = role
    par = SPEC[f]+(end,)
    goal = tuple(sum(v) for v in x)+(3,)
    d, k = [2], [0]

    def rates():
        gains = [v-u-1 for v,u in zip(d,k)]
        for i, j in enumerate(par[1:len(d)],1):
            gains[i] += Q(1,d[j])
            gains[j] += Q(1,d[i])
        buds = [Q(v-u-1,v) for v,u in zip(d,k)]
        # Worst root-parent rate appears in holding, not in the numerator.
        return gains, buds, HI+Q(1,2)+sum(gains)+sum(buds)

    p = S
    for i, wanted in enumerate(goal):
        if i:
            j = par[i]
            # Every parent is processed before its designated child.
            _, buds, hold = rates()
            assert buds[j]>0
            p *= buds[j]/hold
            k[j] += 1
            d.append(2)
            k.append(0)
        while d[i]<wanted:
            gains, _, hold = rates()
            p *= gains[i]/hold
            d[i] += 1
    assert tuple((v-u-1,u,1) for v,u in zip(d,k)) == x+((2,0,1),)
    p /= rates()[2]
    # Shared real parent port: both corners have the same endpoint donor 1/3.
    lows = []
    for r in (a,b):
        z, slopes = r.coefficients(f,x)
        z += slopes[2*end+1]/3
        lo = min(z,z+slopes[0]/2)
        assert lo>=0
        lows.append(lo)
    print(f'  one scripted cylinder: mass/core>={float(p):.10e}, '
          f'full lo/hi floors/core={float(p*lows[0]):.10e}/{float(p*lows[1]):.10e}',flush=True)


def optimize(roles, budgets, coefficients, target, *, exterior=False):
    """Least rectangular label count over the integer grid, proof by enumeration."""
    best = None
    for k in range(0,97):
        def cost(j):
            return sum(remainder(r,j,k,budgets[j],coefficients.get(r),exterior=exterior)
                       for r in roles)
        lo, hi = 17, max(budgets)
        if cost(hi) > target:
            continue
        while hi-lo>1:
            mid = (hi+lo)//2
            if cost(mid) <= target:
                hi = mid
            else:
                lo = mid
        candidate = ((hi-17)*(k+1),hi,k,cost(hi))
        if best is None or candidate < best:
            best = candidate
    return best


def multiplicity_controls():
    total = 0
    for size in range(4,8):
        for prefix in product(*(range(i) for i in range(1,size))):
            par = (-1,)+prefix
            k = Counter(prefix)
            hits = Counter((f,i,ids[i]) for f,ids in embedding(par)
                           for fam,i in ROLES if fam==f)
            for (f,i,y), n in hits.items():
                mul,power = MULT[f,i]
                assert n <= mul*k[y]**power
                off = [w for w in range(1,size) if par[w] == y]
                assert n*(len(off)-SPEC[f].count(i)) <= mul*k[y]**(power+1)
            total += 1
    print('multiplicity fixture:',total,'chronological core trees',flush=True)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument('--lower',type=Path,default=Path('b3f18_lo.bin'))
    p.add_argument('--upper',type=Path,default=Path('b3f18_hi.bin'))
    args = p.parse_args()
    a,b = pairs(args.lower,args.upper)
    multiplicity_controls()
    assert a.c.D == b.c.D == 18
    for r in a,b:
        for domain in r.c.p4,r.c.ybu,r.c.ybw:
            assert max(map(sum,domain)) == 9 and max(x[1] for x in domain) == 7
    coefs, transfers = {}, {}
    for role in ROLES:
        probe(a,b,role)
        coefs[role] = controls(a,b,role)
        transfers[role] = transfer_max(a,b,role)
        print('  global exterior gain curvature factor',str(coefs[role]),
              'internal transfer coefficient ceiling',str(transfers[role]),flush=True)
    budgets = factorial_bounds(range(18,MAX_DEGREE+1),2)
    # Independent exact small-cap reference: laplace tilt plus sums of
    # independent exponentials. The second moment is mean^2 + variance.
    for j in (18,24,43):
        p0, mean, variance = Q(1), 1/S, 1/S**2
        for d in range(2,j+1):
            rate = Q(2*d-1,2)
            p0 *= rate/(rate+S)
            mean += 1/(rate+S)
            variance += 1/(rate+S)**2
        for exact, rounded in zip((p0,p0*mean,p0*(mean**2+variance)),budgets[j]):
            assert 0<=rounded-exact<Q(1,10**19)
    for m in range(3):
        assert all(budgets[j][m] >= budgets[j+1][m] for j in range(18,MAX_DEGREE))
    for name, coefficients, targets, exterior in (
            ('raw internal end-send mass', {}, ('1e-2','1e-3','1e-4'), False),
            ('paired internal transfer',transfers,('1e-5','1e-6','1e-7'), False),
            ('paired exterior end-child price',coefs,('1e-6','1e-7','1e-8'), True)):
        for target in targets:
            result = optimize(ROLES,budgets,coefficients,Q(target),exterior=exterior)
            if result is None:
                print(name,target,'not met on grid')
            else:
                cells,j,k,err = result
                print(name,target,f'J={j} K={k} rectangle={cells} remainder<={float(err):.10e}')


if __name__ == '__main__':
    main()
