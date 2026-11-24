#!/usr/bin/env python3
"""Independent Gold2 lower H3 corner, mask, and literal seam checks."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import lcm, floor
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
R2 = HERE.parent
BASE = run_path(str(HERE/'gold2-h3-sol.py'))
WIRE, WRITER = BASE['wire'], BASE['writer']
P = R2/'occupation-quintet-gold2/higher/cohort.py'
GRID = 2**20
FULL = Q(5988193558027003984698925131,8937592143119350239749734400000000)
SEAM = Q(5820589581257439366906303,388590962744319575641292800000000)
OLD = Q(1084423030952180046412413909,1704710415481443450880000000000000)
SAFE_SEAM = Q(229762209993424478073686403,17875184286238700479499468800000000)
SAFE_FULL = Q(5481957648397036499283,14300147428990960383599575040000000)


def slope_controls():
    # One flexible gain, unit physical width, d=2 units per unit rate,
    # discount 1, successor 0. True uncut reward is 2-rho.
    # The genuine adverse corner is rho=1 and pays value 1/2.
    assert WRITER.affine(Q(2),[(0,2,-WRITER.GRID,0)],2,Q(1),False)==WRITER.GRID//2
    wrong_scale = WRITER.affine(Q(2),[(0,2,-WRITER.GRID//2,0)],2,Q(1),False)
    assert wrong_scale==3*WRITER.GRID//4 and 1-2*Q(wrong_scale,WRITER.GRID)<0
    b = Q(1,3)
    assert -(-GRID*b//1)>0>floor(-GRID*b)
    print('planted controls: wrong sign and unscaled slope both overpay',flush=True)


def states(D,C,k):
    return [(d,c) for d in range(D,1,-1) for c in range(min(d,C),k-1,-1)]


def solve(e, lo, hi, low, high, young):
    """Rolling tables; enumerate the adverse corners rather than sort slopes."""
    ds = (high,5,young)
    den = lcm(hi.denominator,*range(2,max(ds)+1))
    inv = [0,0]+[den//d for d in range(2,max(ds)+1)]
    discount = hi.numerator*(den//hi.denominator)
    mother = states(high,8,1)
    child = states(young,3,1)
    mid = states(5,4,1)
    counts = Counter()

    def corners(p, adj, gain, bud, a=0, slopes=()):
        gs, ws, bs = [], [], []
        for (d,c), ns in zip(p,adj):
            gs.append((d-c)*den+sum(inv[p[j][0]] for j in ns))
            ws.append((c-len(ns))*inv[2])
            bs.append((d-c)*inv[d])
        n = sum(g*v+b*w for g,v,b,w in zip(gs,gain,bs,bud))
        z = sum(gs)+sum(bs)+discount
        floor_zero = min((n+sum(w*v for w,v,k in zip(ws,gain,bits) if k))//
                         (z+sum(w for w,k in zip(ws,bits) if k))
                         for bits in product(*[(0,1) if w else (0,) for w in ws]))
        if not slopes:
            return floor_zero
        counts['paying candidates'] += 1
        n += floor(-GRID*den*a)
        gains = [v+floor(-GRID*b) for v,b in zip(gain,slopes)]
        floor_signed = min((n+sum(w*v for w,v,k in zip(ws,gains,bits) if k))//
                           (z+sum(w for w,k in zip(ws,bits) if k))
                           for bits in product(*[(0,1) if w else (0,) for w in ws]))
        return max(floor_zero,floor_signed)

    def up(x): return x[0]+1,x[1]
    def born(x): return x[0],x[1]+1

    paths = {}
    edges = {}
    solo = {}
    J = (2,1)
    for r in mother:
        triples = {}
        for h in mid:
            if h[1]<2: continue
            for t in child:
                g = (paths.get((up(r),h,t),0),triples.get((up(h),t),0),triples.get((h,up(t)),0))
                b = (paths.get((born(r),h,t),0),triples.get((born(h),t),0),triples.get((h,born(t)),0))
                a, s = 0, ()
                if low<=r[0]<=high:
                    a, ss = e.triple(h,r,t)
                    ss = (ss[1],ss[0],ss[2])
                    widths = (r[1]-1,h[1]-2,t[1]-1)
                    if a+sum(min(0,v)*Q(n,2) for v,n in zip(ss,widths))<0:
                        s = ss
                triples[h,t] = corners((r,h,t),((1,),(0,2),(1,)),g,b,a,s)
        paths.update({(r,*k):v for k,v in triples.items()})
        for h in mid:
            edges[r,h] = corners((r,h),((1,),(0,)),(edges.get((up(r),h),0),edges.get((r,up(h)),0)),
                                 (edges.get((born(r),h),0),edges.get((r,born(h)),0)+triples.get((born(h),J),0)))
        solo[r] = corners((r,),((),),(solo.get(up(r),0),),(solo.get(born(r),0)+edges.get((born(r),J),0),))
    print('independent terminal/edge/single counts:',len(paths),len(edges),len(solo),dict(counts),flush=True)
    return lo*Q(solo[J],GRID*e.snap.scale)


def literal_seam(e, p, alt):
    # Reuse the earlier independent ordinary-arrow oracle, widening only its
    # test mask and respecting the smaller Young cap in the added seam slice.
    original = BASE['literal']
    ns = original.__globals__
    old = ns['continuation']
    def cont(path, marks):
        if marks[path[0]][0]>79:
            if marks[path[0]][0]>100 or (len(path)==3 and marks[path[-1]][0]>5): return 0
            marks = dict(marks)
            marks[path[0]] = (79, marks[path[0]][1])
        return old(path,marks)
    ns['continuation'] = cont
    counts = Counter()
    try:
        for d in (79,80,87,88,89,100,101):
            original(p,e,(None,0,1,2,2,1,3),(7,d,5,4,2,5,2),counts)
            for t in ((2,1),(4,3),(5,2)):
                a,b=e.triple((5,3),(d,3),t)
                aa,bb=alt.affine(((d,3),(5,3),t))
                assert (a,(b[1],b[0],b[2]))==(aa,bb)
    finally:
        ns['continuation'] = old
    assert counts['rejected donor-rate mutations'] and counts['nonzero newborn']
    print('literal mixed-tail and boundary controls:',dict(counts),flush=True)


def main():
    slope_controls()
    assert sha256(P.read_bytes()).hexdigest()=='b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c'
    w=run_path(str(P))
    with w['host'].open_host('lower') as h, WRITER.Snapshot('lower') as snap:
        e=WRITER.Rows(snap,h.p,h.golden)
        assert (e.sig,e.hcap)==(-1,408)
        assert h.golden['verify']['sha256']==WIRE.DIGEST['lower']
        literal_seam(e,h.p,w['Rows'](h))
        amounts = []
        for low,high,young,want in ((80,100,5,SEAM),(19,79,18,FULL)):
            got=solve(e,*h.prior,low,high,young)
            assert 0<got<want,(got,want)
            assert got==(SAFE_SEAM if low==80 else SAFE_FULL)
            amounts.append(got)
            print('exact normalized',low,high,got,'ticks/core',got*2**40,float(got*2**40),flush=True)
        snap.check_unchanged()
    print('corrected disjoint floors combined ticks/core:',sum(amounts)*2**40,float(sum(amounts)*2**40))


if __name__=='__main__':
    main()
