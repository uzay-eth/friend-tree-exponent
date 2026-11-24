#!/usr/bin/env python3
"""Rational controls for normalized motif-lineage truncation."""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'higher'))
from extract import Cert, NB

LO, HI = Q('.565673570322'), Q('.565775')
THETA, GMIN = Q('.976434'), Q('.12510370')
EPS = Q('1e-8')


def genealogy():
    """Check every uniquely newest-node assignment on labelled increasing trees."""
    from itertools import combinations, product
    from extract import embedding

    counts = dict.fromkeys(('f','g','hV','hC','P4','Yt','Yb','S3'), 0)
    trees = 0
    for n in range(1, 8):
        for ps in product(*(range(i) for i in range(1,n))):
            p = (-1,)+ps
            rebuilt = set()
            for newest in range(n):
                # Independent connected subsets; leaf deletion is verified
                # below without using the motif extractor to generate them.
                for m in range(1,min(newest+1,4)+1):
                    for old in combinations(range(newest),m-1):
                        ids = old+(newest,)
                        if sum(p[i] not in ids for i in ids) != 1:
                            continue
                        if m > 1:
                            assert p[newest] in old
                            assert sum(p[i] not in old for i in old) == 1
                        rebuilt.add(frozenset(ids))
            em = embedding(p)
            assert rebuilt == {frozenset(ids) for _,ids in em}
            assert len(rebuilt) == len(em)
            for name,_ in em:
                counts[name] += 1
            trees += 1
    assert all(counts.values())
    print('genealogy', trees, 'trees;', sum(counts.values()), 'unique motif births')


def normalized():
    # omega_m(1) includes the original virgin residence. At any s in the
    # prior, bound the distance from the normalized m-generation law.
    prefactor = (HI+2)*Q(50,7)/(GMIN*(1-THETA))
    for osc in (1, Q('1e-4')):
        m = 1
        while prefactor*osc*THETA**m > EPS:
            m += 1
        assert prefactor*osc*THETA**(m-1) > EPS
        print('generation', m, 'oscillation', float(osc), 'normalized error <=', float(prefactor*osc*THETA**m))

    # With the degree-dependent residual W and C=number of cores, the sum
    # of the weights of intrinsic promotion successors is <=(C+.401) W.
    rho,delta = Q('.448'), Q('.117')
    sys.path.insert(0, str(HERE.parents[1]/'finite'))
    from residual_weight import weight_table, z
    d,q = Q('.18'), Q('.268')
    h = weight_table(d,q,False)
    hmin,hmax = min(h.values()),max(h.values())
    gamma = max([z(20,d)/h[20,1]] + [Q(D-c,D)*(h[D,c+1]-h[D,c]+z(D,d))/h[D,c]
                for D in range(2,20) for c in range(1,D)])
    assert gamma == Q(14283051730,35661533271) < Q('.401')
    assert gamma*delta*hmax/Q('.23') <= hmin
    assert rho+delta <= LO
    # First-exit factor already needs this N at unit multipler. This tests
    # the exact power inequality with no rounded log arithmetic.
    assert (10**68)**117 < (10**8)**1000 < (10**69)**117
    print('core-count envelope exponent', float(delta), 'unit-prefactor 1e-8 requires N > 1e68')



def pressure_fixture():
    # Founder 0, bare parent 1, virgin tracked vertex 2. Enumerate every
    # physical directed arrow, then read the whole post-event pressure.
    def evaluate(leaves):
        adj = [[1]+list(range(4,4+leaves)), [0,2], [1,3], [2]]+[[0] for _ in range(leaves)]

        def read(g):
            return sum(Q(1,len(g[v])) for v in g[2] if len(g[v]) >= 2)

        before = read(adj)
        drift = Q(0)
        for ns in adj:
            for receiver in ns:
                after = [a[:] for a in adj]+[[receiver]]
                after[receiver].append(len(adj))
                drift += (read(after)-before)/len(ns)
        return before,drift

    a,b = (evaluate(n) for n in (1,3))
    assert a == (Q(1,2),Q(1,12)) and b == (Q(1,2),Q(1,8))
    print('one-step pressure defect', b[1]-a[1])


def virgin_tail(path, expected_mode):
    """Check the true first-body row at k=0, L>=LF+1, to all degrees."""
    c = Cert(path)
    assert (c.mode,c.D,c.S) == (expected_mode,18,10**11)
    sign = 1 if c.mode == 2 else -1
    f = int(c.f[c.vi(c.LF,0,1)])
    fb = int(c.f[c.vi(c.LF,1,1)])
    fn = int(c.f[c.vi(*NB)])
    new = [int(c.g[c.state[x],c.state[NB]]) for x in c.tail(0,1,1)]
    asym = Q(sign,c.S)*(fb+fn+new[0]-(1+c.lam)*f)
    # |B3(t)-B3(0)| <= 3t max |new_i-new_0|, t=18/d.
    err1 = abs(fb+fn+new[0]-f)+54*max(abs(x-new[0]) for x in new)
    center = Q(int(c.trpar[0]),c.Sr)
    d0 = c.LF+2
    err2 = max(abs(center),abs(center-Q(1,2)))*(abs(int(c.f2[0]))+Q(abs(int(c.f3[0])),d0))
    def outer(d):
        return Q(err1,d*c.S) + err2/(d*d*c.S)

    # Also output the exact feasible signed row range at the first
    # all-stock point. The large-D envelope above uses no actual samples.
    d = d0
    stock = Q(d-1,d)*(fb+fn+sum(w*int(c.g[c.state[x],c.state[NB]])
                                for x,w in c.embed((d-2,1,1)))-f)-c.lam*f
    psi = Q(int(c.f2[0]),d*d)+Q(int(c.f3[0]),d**3)
    vals = [Q(sign,c.S)*(stock+(center-r)*psi) for r in (0,Q(1,2))]
    assert all(abs(v-asym)<=outer(d) for v in vals)
    a,b,g,h = new
    poly = [fb+fn+a-f, 54*(b-a), 3*18**2*(a-2*b+g), 18**3*(-a+3*b-3*g+h)]
    # Expand (1-1/d) times the cubic in 1/d. Add the pressure-dependent
    # transfer only at orders two and three.
    p = [poly[0]-c.lam*f]+[y-x for x,y in zip(poly,poly[1:])]+[-poly[-1]]
    for d in (d0, 2*d0, 10**6):
        for r in (Q(0),Q(1,3),Q(1,2)):
            expected = sum(y*Q(1,d)**j for j,y in enumerate(p))
            expected += (center-r)*(Q(int(c.f2[0]),d**2)+Q(int(c.f3[0]),d**3))
            direct = Q(d-1,d)*(fb+fn-f+sum(w*int(c.g[c.state[x],c.state[NB]])
                    for x,w in c.embed((d-2,1,1))))-c.lam*f
            direct += (center-r)*(Q(int(c.f2[0]),d**2)+Q(int(c.f3[0]),d**3))
            assert expected == direct
    print('mode',c.mode,'k=0 limiting signed f',float(asym),
          'range at',d0,*(float(v) for v in sorted(vals)),
          'rigorous error at D>=1e6',float(outer(10**6)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('payload', type=Path, nargs='*')
    a = p.parse_args()
    genealogy()
    normalized()
    pressure_fixture()
    if a.payload:
        assert len(a.payload) == 2
        for path,mode in zip(a.payload,(2,1)):
            virgin_tail(path,mode)
