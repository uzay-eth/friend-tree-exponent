"""Disjoint D18 suffix: marked newborn classes and physical crown occupation."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('higher', HERE.parents[1]/'higher'/'extract.py')
H = module_from_spec(spec)
spec.loader.exec_module(H)
P = H.physical
LO, HI = H.LO, H.HI
GRID = 10**23
LABELS = (2, 3, 5, 10, 20)
ENDS = LABELS[1:] + (21,)
FIRST = {'lower': Q('0.000002512846393'), 'upper': Q('0.000002025752646')}
OLD = {'lower': Q('0.000001276902849680'), 'upper': Q('0.000001501029458674')}
DEN = {'lower': Q('1.356833021311'), 'upper': Q('1.367169995210')}
EDGE = {'lower': Q('25227936004759/1700000000000000000000'),
        'upper': Q('5046841434259/340000000000000000000')}


def floor(x):
    return Q(x*GRID // 1, GRID)


def newborn_kernel(residual):
    """Grouped lower direct-birth kernel, with separate beta and ell laws."""
    mat = []
    for i in LABELS:
        row = []
        for a, z in zip(LABELS, ENDS):
            V = defaultdict(Q)
            for d in range(26, 1, -1):
                for k in range(d-1, -1, -1):
                    l = d-k-1
                    b = Q(l, d)
                    r = b * (a <= min(d,20) < z)
                    kill = Q(1,i) if residual else 0
                    end = l+Q(k,2)+(0 if residual else Q(1,i))
                    val = min((r + b*V[d,k+1] + g*V[d+1,k])/(HI+kill+b+g)
                              for g in (l,end))
                    V[d,k] = floor(val)
            row.append(V[2,0])
        mat.append(row)
    m = [min(row[j] for row in mat) for j in range(len(LABELS))]
    for _ in range(80):
        left = 1-sum(m)
        assert left >= 0
        new = [max(m[j], sum(p*row[j] for p,row in zip(m,mat)) + left*min(row[j] for row in mat))
               for j in range(len(m))]
        m = list(map(floor,new))
    return m, mat


def stationary_dual(values, mat):
    """Exact optimal value for row-wise kernel lower bounds / free exits."""
    n = len(values)
    # Complete the missing mass to class zero. Solve
    # alpha = f_i + (A v)_i - v_i with v_0=0.
    rows = [[Q(1), *(Q(i==j)-r[j] for j in range(1,n)), f]
            for i,(r,f) in enumerate(zip(mat,values))]
    for i in range(n):
        k = next(k for k in range(i,n) if rows[k][i])
        rows[i],rows[k] = rows[k],rows[i]
        rows[i] = [x/rows[i][i] for x in rows[i]]
        for k in range(n):
            if k!=i:
                rows[k] = [a-rows[k][i]*b for a,b in zip(rows[k],rows[i])]
    alpha = rows[0][-1]
    v = [Q(0)] + [r[-1] for r in rows[1:]]
    if min(v)<0:
        return Q(0)
    assert all(f+sum(x*y for x,y in zip(r,v))-w == alpha
               for r,w,f in zip(mat,v,values))
    return alpha


@cache
def ncore(t):
    return 1+sum(map(ncore,t[2]))


def evaluate(E, name, cap, masses, ell, detail=False):
    @cache
    def geometry(t):
        n,d,_ = P.statistics(t)
        p = Q(1,d)
        return n, p, P.intrinsic_moves(t)

    @cache
    def score(t):
        n,p,_ = geometry(t)
        if n<=12 or ncore(t)<3:
            bonus, slope = Q(0), Q(0)
        else:
            rows, slope = E.charges(t)
            bonus = sum(rows.values())
            assert bonus+min(0,slope/2) >= 0
        g = EDGE[name] * sum(1/p < 18 and P.statistics(ch)[1]<18 for ch in t[2])
        return bonus, slope, g

    def solve(i, residual=False, what='both'):
        @cache
        def V(t, mark, copies):
            n,p,moves = geometry(t)
            if n>cap or ncore(t)>4:
                return Q(0)
            rmax = Q(1,mark)
            reward,slope,g = score(t)
            if what in ('higher', 'three', 'four'): g = 0
            if ((what == 'three' and ncore(t)!=3) or (what == 'four' and ncore(t)!=4)):
                reward,slope = Q(0), Q(0)
            if what == 'edge': reward,slope = 0,0
            reward += g
            z = sum(q*V(child,mark,copies) for child,q in moves)
            den = HI + n-p
            if residual:
                # Plus-copy retains the true 1/D; after cutting to a child
                # use the largest possible kill/score bound for its mark.
                # Deliberately omit plain copies and all further direct births.
                if copies and any(t[2]):
                    z += p*sum(V(P.plus(y),2,copies-1) for y in t[2])
                return floor((reward+min(0,slope*rmax)+z)/(den+rmax))
            adv = V(P.plus(t),mark,copies) if n < cap else 0
            v = floor(min((reward+z+r*(slope+adv))/(den+r) for r in (0,rmax)))
            assert (den*v <= reward+z and (den+rmax)*v <= reward+z+rmax*(slope+adv))
            return v

        ans = V(P.newborn(2),i,2 if residual else 0)
        V.cache_clear()
        return ans

    def payoff(values, kernel):
        weights, mat = kernel
        ref = values[0]
        assert all(a>=ref for a in values)
        return LO*max(ref+sum(w*(v-ref) for w,v in zip(weights,values)), stationary_dual(values,mat))

    pieces = {}
    for what in (('higher','three','four','edge','both') if detail else ('both',)):
        vals = [solve(i,what=what) for i in LABELS]
        print(name, what, 'physical base', float(LO*vals[0]), 'stratified', float(payoff(vals,masses)), flush=True)
        pieces[what] = payoff(vals,masses)
    physical = max(pieces['both'], pieces.get('higher',0)+pieces.get('edge',0),
                   sum(pieces.get(x,0) for x in ('three','four','edge')))
    renewal = Q(0)
    if detail:
        rv = [solve(i,True) for i in LABELS]
        renewal = payoff(rv,ell)
        print(name, 'marked residual subset, <=2 plus-copy resets', float(renewal), 'geometries', geometry.cache_info().currsize, flush=True)
    return max(physical,renewal)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('blob', type=Path, nargs='+')
    p.add_argument('--cap', type=int, default=18, choices=range(13,21))
    p.add_argument('--detail', action='store_true')
    a = p.parse_args()
    masses = newborn_kernel(False)
    ell = newborn_kernel(True)
    print('beta mass floors', LABELS, *(float(x) for x in masses[0]), flush=True)
    print('ell mass floors', LABELS, *(float(x) for x in ell[0]), flush=True)
    P.MARK_CAP=2  # mark-blind geometry; true birth weights computed above
    paid = {}
    for blob in a.blob:
        c = H.Cert(blob)
        name = ('upper' if c.mode==1 else 'lower')
        assert (c.D,c.K,c.S,c.lam) == (18,17,10**11,Q('.565775' if c.mode==1 else '.565673'))
        next_pay = floor(evaluate(H.Extract(c),name,a.cap,masses,ell,a.detail))
        base = FIRST[name]+OLD[name]
        total = base+next_pay
        sign = 1 if name=='lower' else -1
        base_endpoint = c.lam+sign*base/DEN[name]
        print(name,'base',str(base),'base endpoint',float(base_endpoint),
              'new',str(floor(next_pay)),'total',str(floor(total)),flush=True)
        result = sign*((sign*(c.lam+sign*total/DEN[name])*10**12)//1)
        print(name,'outward',Q(result,10**12),'decimal',f'{result//10**12}.{result%10**12:012d}',
              'gain against cumulative base',float(next_pay/DEN[name]),flush=True)
        paid[name] = (Q(result,10**12),total,c.lam)
    if len(paid)==2:
        lo,hi = paid['lower'][0],paid['upper'][0]
        if a.cap>=18:
            assert lo >= Q('0.565676153882') and hi <= Q('0.565772151582')
        mid=(lo+hi)/2
        for name,(edge,total,lam) in paid.items():
            width = (hi-lam) if name=='lower' else (lam-lo)
            print(name,'omitted slack budget full interval',float(width*DEN[name]-total),
                  'at midpoint',float(abs(mid-lam)*DEN[name]-total),flush=True)
        sums = [sum(abs(s-lam)*DEN[name]-total for name,(_,total,lam) in paid.items())
                for s in (lo,hi)]
        print('common omitted sum cap',float(max(sums)),flush=True)


if __name__=='__main__':
    main()
