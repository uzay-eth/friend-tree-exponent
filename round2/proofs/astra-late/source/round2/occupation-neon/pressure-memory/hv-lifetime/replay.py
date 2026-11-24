#!/usr/bin/env python3
"""All-age stationary signed hV rows, with branching physical path registration."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cache
import hashlib
from math import floor, ceil, lcm
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'higher' / 'older'))
from replay import Cert, Local, LO, HI, PAYLOAD_SHA, sample

GRID = 10**19
VIRGIN = (2, 0)
TAIL = (0, 0)
ZERO = (0,) * 6
RATE = lcm(10**12, *range(2, 8))
NO_CREDIT = (Q(0),)*3


def paid(path, disjoint):
    if len(path) != 3:
        return False
    if not disjoint:
        return True
    n = sum(k for _, k in path) - 2
    # Every side child was born with a two-vertex fringe. The full fringe
    # never loses vertices, including on promotions of those leaves.
    return 3 + n >= 5 and sum(d - k for d, k in path) + 2*n > 12


@cache
def row(E, path):
    x = tuple((d-k-1, k, 1) for d, k in path)
    intercept = E.charge('hV', x, [(0, 0)]*3)
    raw = E.terms('hV', x)[1]
    slopes = []
    for i, parent, n in ((0, True, 1), *((j, False, x[j][1]-(j<2)) for j in range(3))):
        c = Q(0)
        if n:
            shift = E.transfer('hV', x, i) if E.moved('hV', x, i, parent) else 0
            c = E.sign * (raw[i]-shift) / E.c.S
        slopes.append(c)
    low, high, _ = E.extrema('hV', x)
    assert low == intercept + sum(min(0, c*n/2) for c, n in zip(slopes, (1, *(k-(i<2) for i, (_,k) in enumerate(path)))))
    return intercept, *slopes


@cache
def rate(q):
    n,d=q.as_integer_ratio()
    assert RATE%d==0
    return n*(RATE//d)


@cache
def fixed_reward(rewards):
    den=lcm(*(q.denominator for q in rewards))
    scale=GRID*RATE
    return den, [r.numerator*(den//r.denominator)*scale for r in rewards for _ in (0,1)]


@cache
def scaled_credit(credit, den):
    return tuple(c.numerator*(den//c.denominator)*GRID for c in credit for _ in (0,1))


def integer_row(rewards, edges, lo, hi):
    """Max/min affine-ratio rows simultaneously using exact integer arithmetic."""
    d0, r0 = fixed_reward(rewards)
    den = lcm(d0, *(q.denominator for _,_,_,c in edges if c is not NO_CREDIT for q in c))
    base = rate(lo)
    swing = rate(hi)-base
    nums = [v*(den//d0) for v in r0]
    opts = [[(0,swing)] if swing else [] for _ in nums]
    for a,b,v,c in edges:
        A=rate(a); B=rate(b)-A
        if not (A or B):continue
        base += A
        credit = ZERO if c is NO_CREDIT else scaled_credit(c,den)
        for j in range(6):
            f = v[j]*den + credit[j]
            nums[j] += A*f
            if B:
                opts[j].append((f,B))
    assert base > 0
    out = []
    for j,(num,choices) in enumerate(zip(nums, opts)):
        tot = base
        for f,b in sorted(choices, reverse=bool(j%2)):
            if (f*tot <= num if j%2 else f*tot >= num):
                break
            num += f*b
            tot += b
        divisor = tot*den
        out.append(-((-num)//divisor) if j%2 else num//divisor)
    return tuple(out)


def check_arithmetic():
    from itertools import product
    rewards = (Q(1,17),-Q(2,13),Q(1,17)+Q(2,13))
    edges = [(Q(1,3),Q(2,3),(-30,31,40,44,-80,-75), (Q(1,11),Q(2,19),-Q(3,7))),
             (Q(1,7),Q(1,2),(-50,-48,800,803,-850,-820),NO_CREDIT)]
    oracle = []
    for j in range(6):
        values = [(rewards[j//2]*GRID+sum(r*(v[j]+c[j//2]*GRID) for r,(_,_,v,c) in zip(rs,edges))) /
                  (mu+sum(rs)) for mu,*rs in product((LO,HI),*((a,b) for a,b,_,_ in edges))]
        oracle.append(ceil(max(values)) if j%2 else floor(min(values)))
    assert tuple(oracle)==integer_row(rewards,edges,LO,HI)


def enclosure(E, degree, arity, donor, mother_cap, disjoint, lo, hi):
    new = VIRGIN if donor else TAIL
    tail = Q(1, max(2,donor+1))
    zero = NO_CREDIT

    @cache
    def future(path, sides, mother):
        n = len(path)
        if any(d > degree or k > arity for d,k in path) or (n<3 and path[-1][1] == arity):
            return ZERO
        if paid(path, disjoint):
            a,b = (row(e,path) for e in E)
            scores = tuple((x,y,x-y) for x,y in zip(a,b))
        else:
            scores = (zero,)*5
        edges = []
        for i,(d,k) in enumerate(path):
            l = d-k-1
            gain = future(path[:i]+((d+1,k),)+path[i+1:],sides,mother)
            intrinsic = Q(l) + sum((Q(1,path[j][0]) for j in (i-1,i+1) if 0<=j<n),Q(0))
            edges.append((intrinsic,intrinsic,gain,zero))
            if not i:
                p = Q(1,mother) if mother else Q(0)
                edges.append((p,p if mother else Q(1,max(2,mother_cap+1)),gain,scores[1]))
            registry = sides[i]
            exact = sum((Q(1,dd) for dd,kk in registry if dd),Q(0))
            edges.append((exact,exact+registry.count(TAIL)*tail,gain,scores[i+2]))
            if l:
                post = path[:i]+((d,k+1),)+path[i+1:]
                reg = sides[:i]+(tuple(sorted(registry+(new,))),)+sides[i+1:]
                value = future(post,reg,mother)
                if i == n-1 and n<3:
                    extension = future(post+(VIRGIN,),sides+((),),mother)
                    value = tuple(x+y for x,y in zip(value,extension))
                edges.append((Q(l,d),Q(l,d),value,zero))
            for (dd,kk), mult in Counter(registry).items():
                if not dd:
                    continue
                leaves = dd-kk-1
                following = (dd+1,kk) if dd<donor else TAIL
                def change(new_state):
                    q=list(registry); q.remove((dd,kk)); q.append(new_state)
                    return future(path,sides[:i]+(tuple(sorted(q)),)+sides[i+1:],mother)
                gain_rate = leaves+Q(1,d)
                edges.append((mult*gain_rate,mult*(gain_rate+Q(kk,2)),change(following),zero))
                if leaves:
                    birth=Q(mult*leaves,dd)
                    edges.append((birth,birth,change((dd,kk+1)),zero))
        if mother:
            p=Q(1,path[0][0])
            later=mother+1 if mother<mother_cap else 0
            edges.append((p,p+mother-1,future(path,sides,later),zero))
        return integer_row(scores[0], edges, lo, hi)

    start = {d: future((VIRGIN,),((),),d) for d in (0,*range(2,mother_cap+1))}
    bounds = []
    for j in (0,2,4):
        v=min(s[j] for s in start.values()); w=max(s[j+1] for s in start.values())
        bounds.append((min(lo*v,hi*v)/GRID,max(lo*w,hi*w)/GRID))
    n=future.cache_info().currsize
    future.cache_clear()
    return bounds,n


def finite_check(E, degree, arity, disjoint, size, seed):
    from random import Random
    p,x=sample(size,Random(seed)); totals=[Q(0)]*2
    kids=[[] for _ in x]
    for j,i in enumerate(p[1:],1):kids[i].append(j)
    count=0
    for a in range(1,len(x)):
        for b in kids[a]:
            for c in kids[b]:
                ids=(a,b,c); path=tuple((sum(x[i]),x[i][1]) for i in ids)
                if any(d>degree or k>arity for d,k in path) or not paid(path,disjoint):continue
                if disjoint:
                    fringe=[a]
                    for v in fringe:fringe.extend(kids[v])
                    assert len(fringe)>=5 and sum(1+x[v][0] for v in fringe)>12
                outside=[Q(1,sum(x[p[a]]))]+[sum((Q(1,sum(x[j])) for j in kids[i] if j not in ids),Q(0)) for i in ids]
                for j,e in enumerate(E):
                    v,*coeff=row(e,path)
                    v+=sum(a*b for a,b in zip(coeff,outside))
                    assert v==e.charge('hV',tuple(x[i] for i in ids),[(outside[0],outside[1]),(0,outside[2]),(0,outside[3])])>=0
                    totals[j]+=v
                count+=1
    return len(x),count,totals


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lower',type=Path); p.add_argument('upper',type=Path)
    p.add_argument('--degree',type=int,choices=range(3,8),default=7)
    p.add_argument('--arity',type=int,choices=(2,3),default=2)
    p.add_argument('--donors',type=int,nargs='+',choices=(0,2,3,4),default=[0,2])
    p.add_argument('--mother',type=int,choices=(0,2,3,4,5,6),default=0)
    p.add_argument('--disjoint',action='store_true');p.add_argument('--sample',type=int,default=16000)
    p.add_argument('--reps',type=int,default=4)
    p.add_argument('--mu-cells',type=int,choices=(1,2,4,8),default=1)
    a=p.parse_args()
    check_arithmetic()
    for side, f in [('lower',a.lower),('upper',a.upper)]:
        with f.open('rb') as stream: digest=hashlib.file_digest(stream,'sha256').hexdigest()
        if digest!=PAYLOAD_SHA[side]: p.error(f'{side} blob hash mismatch')
        print(side,'pinned',digest,flush=True)
    E=tuple(Local(Cert(f)) for f in (a.lower,a.upper))
    if a.sample:
        samples=[finite_check(E,a.degree,a.arity,a.disjoint,a.sample,122419+j) for j in range(a.reps)]
        cores=sum(n for n,_,_ in samples); count=sum(n for _,n,_ in samples)
        totals=[sum(v[j] for _,_,v in samples) for j in range(2)]
        print('finite sample only',cores,'cores',count,'paid hV',*(float(x/cores) for x in (*totals,totals[0]-totals[1])),flush=True)
    for d in a.donors:
        cells=[enclosure(E,a.degree,a.arity,d,a.mother,a.disjoint,LO+(HI-LO)*Q(j,a.mu_cells),
                         LO+(HI-LO)*Q(j+1,a.mu_cells)) for j in range(a.mu_cells)]
        values=[(min(v[j][0] for v,_ in cells),max(v[j][1] for v,_ in cells)) for j in range(3)]
        print('D,k,donor,mother,disjoint',a.degree,a.arity,d,a.mother,a.disjoint,'mu cells',a.mu_cells,
              'states',sum(n for _,n in cells),flush=True)
        for name,(lo,hi) in zip(('lower','upper','lower-upper'),values):
            print(name,'outward 1e-16 ticks',floor(lo*10**16),ceil(hi*10**16),flush=True)


if __name__=='__main__':main()
