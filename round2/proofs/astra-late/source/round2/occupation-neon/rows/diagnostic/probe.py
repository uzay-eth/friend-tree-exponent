"""Sparse exact eight-body residuals on seeded physical friend-tree trajectories."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import prod
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / 'attacks/neon1'))
from audit import Cert, arrows, physical_view  # noqa: E402

NAMES = ('f', 'g', 'hV', 'hC', 'p4', 'yT', 'yB', 's3')


def motifs(c, a, children, st):
    yield 0, (a,)
    sa = st[a]
    ch = children[a]
    for b in ch:
        yield 1, (a,b)
        sb = st[b]
        for d in children[b]:
            if sb[1] <= c.K:
                yield 2, (a,b,d)
            if sb in c.p4 and st[d] in c.p4:
                for e in children[d]:
                    yield 4, (a,b,d,e)
            if sa in c.ybu and sb in c.ybw:
                for w in ch:
                    if w != b:
                        yield 6, (a,b,w,d)
        if sb in c.yt:
            for d,e in combinations(children[b],2):
                yield 5, (a,b,d,e)
    if 2 <= sa[1] <= c.K:
        for b,d in combinations(ch,2):
            yield 3, (a,b,d)
    if sa in c.s3:
        for b,d,e in combinations(ch,3):
            yield 7, (a,b,d,e)


def evaluator(c):
    @lru_cache(maxsize=30000)
    def value(i, st):
        if i == 0:
            return Q(int(c.f[c.vi(*st[0])]))

        def cell(x):
            a,b,*r = x
            if i == 1: return c.g[c.state[a],c.state[b]]
            d,*e = r
            if i == 2: return c.hvtab[c.par[a],c.hv[b],c.child[d]] if b in c.hv else 0
            if i == 3: return c.hctab[c.hc[a],c.tp(b,d)] if a in c.hc else 0
            e, = e
            if i == 4: return c.p4tab[c.par[a],c.p4[b],c.p4[d],c.child[e]]
            if i == 5: return c.yttab[c.par[a],c.yt[b],c.tp(d,e)]
            if i == 6: return c.ybtab[c.ybu[a],c.ybw[b],c.child[d],c.child[e]]
            return c.s3tab[c.s3[a],c.tt(b,d,e)]

        return sum((int(cell(x))*prod(w) for row in product(*(c.embed(s) for s in st))
                    for x,w in [tuple(zip(*row))]), Q())
    return value


def sparse(c, value, par):
    children, st, core = physical_view(par)
    C = len(core)
    up = [-1]*C
    incidence = [[] for _ in core]
    old = {}
    for a,ch in enumerate(children):
        for b in ch: up[b] = a
        for key in motifs(c,a,children,st):
            i, ids = key
            old[key] = value(i,tuple(st[j] for j in ids))
            for j in ids: incidence[j].append(key)
    base = [Q() for _ in NAMES]
    by_anchor = [[Q() for _ in NAMES] for _ in core]
    drift = base.copy()
    for (i,ids),v in old.items():
        base[i] += v
        by_anchor[ids[0]][i] -= c.lam*v
    for j,(L,k,h) in enumerate(st):
        rho = sum((Q(1,sum(st[u])) for u in children[j]), Q())
        if h: rho += Q(1,sum(st[up[j]]))
        anchors = [j]
        while up[anchors[-1]] >= 0 and len(anchors)<4: anchors.append(up[anchors[-1]])
        for birth, rate in ((False,L+rho),(True,Q(L,L+k+h))):
            if not rate: continue
            if birth:
                st[j] = (L-1,k+1,h)
                st.append((1,0,1)); children.append([]); children[j].append(C)
            else:
                st[j] = (L+1,k,h)
            new = {key:value(key[0],tuple(st[u] for u in key[1]))
                   for a in anchors for key in motifs(c,a,children,st)
                   if j in key[1] or (birth and C in key[1])}
            if birth: new[0,(C,)] = value(0,((1,0,1),))
            for key in incidence[j]: new[key] = new.get(key,0)-old[key]
            for (i,ids),delta in new.items():
                gain = rate*delta
                drift[i] += gain
                by_anchor[j if ids[0]==C else ids[0]][i] += gain
            st[j] = (L,k,h)
            if birth: st.pop(); children.pop(); children[j].pop()
    sign = 1 if c.mode == 2 else -1
    for i in range(8):
        assert sum(row[i] for row in by_anchor) == drift[i]-c.lam*base[i]
    return base, drift, st, [[sign*v/c.S for v in row] for row in by_anchor]


def replay(c, value, par):
    phi = c.potential(*physical_view(par)[:2],parts=True)
    drift = [Q() for _ in NAMES]
    for j,r in arrows(par).items():
        after = c.potential(*physical_view(par+(j,))[:2],parts=True)
        drift = [d+r*(a-b) for d,a,b in zip(drift,after,phi)]
    assert sparse(c,value,par)[:2] == (phi,drift)


def snapshots(seed, sizes):
    rng = random.Random(seed)
    par = [-1,0,0]
    g = [[1,2],[0],[0]]
    for n in sizes:
        while len(par)<n:
            donor = rng.randrange(len(par))
            target = rng.choice(g[donor])
            g[target].append(len(par)); g.append([target]); par.append(target)
        yield tuple(par)


def report(c, trees, sizes):
    value = evaluator(c)
    # Fixed reachable witnesses: mixed tail, deep path, branch and fork motifs.
    core = (-1,0,0,0,1,1,2,4,7)
    for fixture in [(-1,0,0), core+tuple(range(9))+(0,)*18+(4,)*18,
                    next(snapshots(199,(32,)))]:
        replay(c,value,fixture)
    for n in sizes:
        phi = Q(); core_n = 0
        totals = [Q() for _ in NAMES]
        classes = defaultdict(lambda: [0,Q()])
        samples = []
        for par in trees[n]:
            base,_,st,items = sparse(c,value,par)
            p = sum(base)/c.S
            phi += p; core_n += len(st)
            local = sum(map(sum,items)); assert local >= 0
            samples.append(float(local/p))
            for (L,k,h),row in zip(st,items):
                for i,v in enumerate(row): totals[i] += v
                for key in (('state',L,k,h),('band',min(L,2),min(k,3),h),
                            ('degree','inner' if L+k+h<c.D else 'tail',h)):
                    classes[key][0] += 1; classes[key][1] += sum(row)
        print(f'n={n} trajectories={len(samples)} mean_C={core_n/len(samples):.3f} '
              f'R/phi={float(sum(totals)/phi):.10g} R/C={float(sum(totals)/core_n):.10g} '
              f'per_tree_R/phi_range=[{min(samples):.7g},{max(samples):.7g}]',flush=True)
        print('  raw eight R/phi:', ' '.join(f'{kind}={float(v/phi):+.8g}' for kind,v in zip(NAMES,totals)))
        for group in ('state','band','degree'):
            for key,(count,total) in sorted((item for item in classes.items() if item[0][0]==group),
                                            key=lambda x:abs(x[1][1]),reverse=True)[:10]:
                print(f'  {key} hits={count} mean={float(total/count):+.8g} R/phi={float(total/phi):+.8g}')
        sys.stdout.flush()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--runs', type=int, default=12)
    p.add_argument('--seed', type=int, default=20261124)
    p.add_argument('--sizes', type=int, nargs='+', default=[256,2048,16384])
    a = p.parse_args()
    assert a.runs>0 and a.sizes == sorted(set(a.sizes)) and min(a.sizes)>=3
    trees = {n:[] for n in a.sizes}
    for i in range(a.runs):
        for n,par in zip(a.sizes,snapshots(a.seed+i,a.sizes)): trees[n].append(par)
    print('SIMULATION DIAGNOSTIC, NOT A BOUND; exact per-tree residual. seed0=',a.seed,flush=True)
    for sign in ('lo','hi'):
        cert = Cert(Path(f'b3f18_{sign}.bin'))
        print(sign,'lambda',cert.lam,flush=True)
        report(cert,trees,a.sizes)


if __name__ == '__main__': main()
