#!/usr/bin/env python3
"""Paired production slacks and directed stable-newborn occupation bounds."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
import gzip
import hashlib
import json
from math import lcm
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'rows'))
import check_edges
sys.path.insert(0, str(HERE.parent / 'renewal'))
import check as root

LO, HI = root.LO, root.HI
UNIT = 10**21
ZERO = Q(0)
TIP = (1, 0, 1)
F_CLASSES = {
    'f-young': lambda d, k: d <= 5 and k <= 1,
    'f-inner': lambda d, k: d <= 17,
    'f-tail-small-k': lambda d, k: 18 <= d <= 64 and k <= 4,
    'f-finite': lambda d, k: d <= 64,
}
G_CLASSES = {
    'g-virgin': lambda p, w: w == TIP,
    'g-young': lambda p, w: sum(p) <= 7 and sum(w) <= 7,
    'g-inner': lambda p, w: True,
}


def text(x, digits=8):
    return f'{float(x):.{digits}g}'  # descriptive only; bound output uses brackets()


def rounded(x, upper=False):
    y = x * UNIT
    return -(-y // 1) if upper else y // 1


def div(n, d, upper):
    return -(-n // d) if upper else n // d


def brackets(lo, hi, places=15):
    def fixed(x, up):
        k = div((x * 10**places).numerator, (x * 10**places).denominator, up)
        return f'{"-" if k < 0 else ""}{abs(k)//10**places}.{abs(k)%10**places:0{places}d}'
    return f'[{fixed(lo, False)}, {fixed(hi, True)}]'


def affine(a, u):
    return a[0] + (a[1] - a[0]) * u


def hinge(a, b, t, independent=False):
    if independent:
        return tuple((u, max(ZERO, max(a)-t*min(b))) for u in (ZERO, Q(1)))
    x, y = (p-t*q for p, q in zip(a, b))
    points = [(ZERO, max(ZERO, x)), (Q(1), max(ZERO, y))]
    if x*y < 0:
        points.append((x/(x-y), ZERO))
    return tuple(points)


def solve_root(pay, upper, limit, robust=False):
    """Nonnegative reward specified at shared-rho endpoints and hinge kinks."""
    discounts = (LO, HI) if robust else (LO if upper else HI,)
    op = max if upper else min
    v = defaultdict(int)
    for d in range(limit, 1, -1):
        for k in range(d-1, -1, -1):
            l = d-k-1
            b = Q(l, d)
            points = pay.get((d, k), ((ZERO, ZERO), (Q(1), ZERO)))
            # Include all supplied kinks in both directions (needed for
            # the lower of a hinge and the upper of concave overlap).
            x = op((rounded(r, upper)+b*v[d, k+1]+g*v[d+1, k])/(discount+b+g)
                   for u, r in points for g in (l+u*Q(k+1, 2),) for discount in discounts)
            v[d, k] = div(x.numerator, x.denominator, upper)
    return (1 if robust else HI if upper else LO)*Q(v[2, 0], UNIT)


def root_range(pay, limit):
    return tuple(solve_root(pay, u, limit) for u in (False, True))


def edge_grid(pmax=17, wmax=17, virgin=False):
    # Dense topological order, both degrees and both child counts decrease.
    # Parent h=1 and k>=1, child h=1. Founder births have zero limiting mass.
    def states(cap):
        return tuple((d, k) for d in range(cap, 1, -1) for k in range(d-1, -1, -1))
    ws = ((2,0),) if virgin else states(wmax)
    keys = tuple((d, k, e, j) for d, k in states(pmax) if k for e, j in ws)
    indices = {key: i for i, key in enumerate(keys)}
    dead = len(keys)
    compiled = []
    for (d, k, e, j) in keys:
        l, n = d-k-1, e-j-1
        next_indices = tuple(indices.get(key, dead) for key in (
            (d+1,k,e,j), (d,k,e+1,j), (d,k+1,e,j), (d,k,e,j+1)))
        den = lcm(2, d, e)
        a = ((l*e+1)*(den//e), (n*d+1)*(den//d), l*(den//d), n*(den//e))
        compiled.append((next_indices, den, a, k*(den//2), j*(den//2)))
    return keys, indices, compiled, pmax


def edge_solve(pay, grid, upper, robust=False):
    """Post-birth integral on a tracked live parent-child pair, then Palm flux."""
    keys, indices, compiled, pmax = grid
    discounts = (LO, HI) if robust else (LO if upper else HI,)
    op = max if upper else min
    v = [0]*(1+len(keys))
    for i, ((ip, iw, jp, jw), den, (a,b,c,d), pext, wext) in enumerate(compiled):
        force = rounded(pay.get(keys[i], ZERO), upper)*den
        vals = (v[ip], v[iw], v[jp], v[jw])
        base = force+a*vals[0]+b*vals[1]+c*vals[2]+d*vals[3]
        v[i] = op(div((base+x*vals[0]+y*vals[1])*s.denominator,
                      (a+b+c+d+x+y)*s.denominator+den*s.numerator, upper)
                  for x in (0, pext) for y in (0, wext) for s in discounts)
    birth = {}
    for d in range(2, pmax+1):
        for k in range(d-1):
            x = v[indices[d,k+1,2,0]] * Q(d-k-1, d*UNIT)
            birth[d,k] = ((ZERO,x), (Q(1),x))
    return solve_root(birth, upper, pmax, robust)


def edge_range(pay, grid, occupation):
    lo, hi = (edge_solve(pay, grid, u) for u in (False, True))
    if occupation:
        lo = max(lo, sum(occupation.get(key, ZERO)*v for key, v in pay.items()))
    return lo, hi


def fuse_rows(f, g):
    return {key: (f[key], g[key]) for key in f}


def samples(fp, gp):
    for key in ((2,0), (3,0), (4,0), (2,1), (5,1), (17,0), (18,0),
                (18,4), (22,2), (25,0), (32,4), (64,0)):
        a,b = fp[key]
        print('f case', key, 'slacks ppm', *(text(x*10**6) for x in (*a,*b)),
              'ratios', *(text(x/y) if y else 'inf' for x,y in zip(a,b)))
    for key in (((0,1,1),TIP), ((1,1,1),TIP), ((1,2,1),TIP), ((1,1,1),(2,0,1))):
        a,b = gp[key]
        print('g case', key, 'slacks ppb', text(a*10**9), text(b*10**9), 'ratio', text(a/b))


def describe(name, pairs, edge=False):
    values = [(x, y, key, u) for key, (a,b) in pairs.items() for u,(x,y) in
              enumerate(zip((a,) if edge else a, (b,) if edge else b))]
    assert all(x>=0 and y>=0 for x,y,_,_ in values)
    def brief(x):
        a,b,k,u = x
        return text(a/b) if b else 'inf', k, u, text((a+b)*10**6)
    extrema = sorted((x for x in values if x[0] or x[1]),
                     key=lambda x:(not x[1], x[0]/x[1] if x[1] else ZERO))
    print(name, 'count',len(pairs),'joint zero corners',sum(not(a or b) for a,b,_,_ in values),
          'ratio min/max (ratio,key,corner,sum ppm)',
          brief(extrema[0]),brief(extrema[-1]),'largest sum',brief(max(values,key=lambda x:x[0]+x[1])))
    xbar = sum(float(x) for x,_,_,_ in values)/len(values)
    ybar = sum(float(y) for _,y,_,_ in values)/len(values)
    dx, dy = zip(*((float(x)-xbar,float(y)-ybar) for x,y,_,_ in values))
    corr = sum(x*y for x,y in zip(dx,dy))/(sum(x*x for x in dx)*sum(y*y for y in dy))**.5
    print('  unweighted grid Pearson',text(corr))


def print_bounds(name, pairs, t, occupation, edge=False, independent=False, ratio=None):
    if t == 0:
        a = {k:x if edge else ((ZERO,x[0]),(Q(1),x[1])) for k,(x,_) in pairs.items()}
        print(name, 't 0 H-',brackets(*occupation(a)), 'H+ [0,0], uncancelled share 1')
        return
    tr = (lambda a,b: ((ZERO,a),(Q(1),b)))
    if edge:
        plus = {k:max(ZERO,a-t*b) for k,(a,b) in pairs.items()}
        minus = {k:max(ZERO,t*b-a) for k,(a,b) in pairs.items()}
        overlap = {k:min(a,t*b) for k,(a,b) in pairs.items()}
        total = {k:a+t*b for k,(a,b) in pairs.items()}
    else:
        plus = {k:hinge(a,b,t) for k,(a,b) in pairs.items()}
        minus = {k:tuple((u,t*r) for u,r in hinge(b,a,1/t)) for k,(a,b) in pairs.items()}
        total = {k:tr(a[0]+t*b[0],a[1]+t*b[1]) for k,(a,b) in pairs.items()}
        overlap = {}
        for k,(a,b) in pairs.items():
            overlap[k] = tuple((u,min(affine(a,u),t*affine(b,u))) for u,_ in plus[k])
        if independent:
            sep = {k:hinge(a,b,t,True) for k,(a,b) in pairs.items()}
            sepm = {k:tuple((u,t*r) for u,r in hinge(b,a,1/t,True)) for k,(a,b) in pairs.items()}
            print(name, 'independent upper hinges', brackets(ZERO,occupation(sep)[1]),
                  brackets(ZERO,occupation(sepm)[1]))
    pp,mm,ss,oo = (occupation(r) for r in (plus,minus,total,overlap))
    # 2 min(A,B)/(A+B) = cancellation. Combining safe intervals for the
    # two nonnegative sides avoids a signed discount-monotonicity argument.
    hl,hu = pp[0]+mm[0], pp[1]+mm[1]
    assert ss[0] > 0 and oo[0] > 0
    loss_lo = max(ZERO, hl/(hl+2*oo[1]), hl/ss[1], 1-2*oo[1]/ss[0])
    loss_hi = min(Q(1), hu/(hu+2*oo[0]), hu/ss[0], 1-2*oo[0]/ss[1])
    assert loss_lo <= loss_hi
    if ratio:
        for upper in (False, True):
            a, b = loss_lo, loss_hi
            for _ in range(8):
                test = (a+b)/2
                # Round the scalar candidate exactly before building the score.
                test = Q(test * 10**10 // 1, 10**10)
                if not a < test < b:
                    break
                if edge:
                    score = {k:abs(x-t*y)-test*(x+t*y) for k,(x,y) in pairs.items()}
                else:
                    score = {k:tuple((u,abs(affine(x,u)-t*affine(y,u)) -
                                      test*(affine(x,u)+t*affine(y,u))) for u,_ in plus[k])
                             for k,(x,y) in pairs.items()}
                result = ratio(score, upper)
                secures_upper = result <= 0 if upper else result < 0
                if secures_upper:
                    b = test
                else:
                    a = test
            if upper:
                loss_hi = b
            else:
                loss_lo = a
    print(name, 't',str(t), 'H-',brackets(*pp),'H+',brackets(*mm),
          'sum',brackets(*ss),'overlap',brackets(*oo),
          'uncancelled share',brackets(loss_lo,loss_hi,8))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--t', type=Q, help='common nonnegative multiplier (default: 1.2415 for f, 1 for g)')
    parser.add_argument('--prefix', action='store_true', help='also use physical positive fringe prefix for g lower')
    parser.add_argument('--tight-ratio', action='store_true', help='signed Bellman bisection of the occupation shares')
    parser.add_argument('--only', nargs='+', choices=tuple(F_CLASSES)+tuple(G_CLASSES))
    args = parser.parse_args()
    assert args.t is None or args.t >= 0
    classes = set(args.only or (*F_CLASSES, *G_CLASSES))
    fpath = HERE.parent/'rows'/'charged-f.json'
    gpath = HERE.parent/'rows'/'pinned-fg-slabs.json.gz'
    assert hashlib.sha256(gpath.read_bytes()).hexdigest() == '63aab8b10d07b43abd80e5852e60b72f5513bc23618d83ce6b0fd3ba3d5ab321'
    dump = json.loads(fpath.read_text())
    with gzip.open(gpath,'rt') as stream:
        pack = json.load(stream)
    assert dump['payload_sha256'] == pack['payload_sha256']
    fp = fuse_rows(*(root.row_rewards(dump,n) for n in ('lower','upper')))
    gp = fuse_rows(*(check_edges.rows(pack[n]['g']) for n in ('lower','upper')))
    samples(fp,gp)
    a,b = fp[22,4]
    t = Q('1.2415')
    same = min(x+t*y for x,y in zip(a,b))
    split = min(a)+t*min(b)
    assert same > 4*split
    # A lower of either hinge must check the interior kink, not just the endpoints.
    assert min(r for _,r in hinge((Q(-1),Q(1)), (Q(0),Q(0)), t)) == 0
    print('D22 k4 shared vs separate weighted sum floors (ppm)',text(same*10**6),text(split*10**6))
    for name, accept in F_CLASSES.items():
        if name in classes:
            pairs = {k:v for k,v in fp.items() if accept(*k)}
            describe(name,pairs)
            limit=max(d for d,k in pairs)
            ratio = (lambda r,u:solve_root(r,u,limit,True)) if args.tight_ratio else None
            print_bounds(name,pairs,Q('1.2415') if args.t is None else args.t,
                         lambda r:root_range(r,limit), independent=limit>=18,ratio=ratio)
    if not classes.intersection(G_CLASSES):
        return
    fringe = defaultdict(Q)
    if args.prefix:
        for tree, x in root.fringe_occupation(8).items():
            _, l, children = tree
            for _, u, grand in children:
                fringe[l+len(children)+1,len(children),u+len(grand)+1,len(grand)] += LO*x
    for name, accept in G_CLASSES.items():
        if name in classes:
            pairs = {key:v for key,v in gp.items() if key[0][2] and accept(*key)}
            describe(name,pairs,True)
            d = {(sum(p),p[1],sum(w),w[1]):v for (p,w),v in pairs.items()}
            grid = edge_grid(max(k[0] for k in d), max(k[2] for k in d), name == 'g-virgin')
            ratio = (lambda r,u:edge_solve(r,grid,u,True)) if args.tight_ratio else None
            print_bounds(name,d,Q(1) if args.t is None else args.t,
                         lambda r:edge_range(r,grid,fringe),True,ratio=ratio)


if __name__ == '__main__':
    main()
