"""Exact degree/stock discharge for the signed hubhood numerator."""
from fractions import Fraction as Q
from collections import defaultdict
from functools import cache
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

H = run_path(str(Path(__file__).resolve().parents[1] / 'check.py'))
PREC = 10**12


@cache
def G(d):
    return Q((PREC*H['serial_lower'](d))//1,PREC)


@cache
def U(d,j):
    return Q(-((-PREC*H['upper'](d,j))//1),PREC)


@cache
def L(d):
    return Q((PREC*H['rescue'](2,1)*H['return_lower'](1,d))//1,PREC) if d < 32 else Q(0)


def positive(x):
    return max(Q(0),x)


@cache
def profile(d, ds):
    """Lower charge at an old core, after paying only its old leaves."""
    c = len(ds)
    l = d-c
    w = (l+sum((Q(1,e) for e in ds if e),Q(0)))/(d+1)
    s = w/d
    out = w*G(d+1)+Q(l,2*d)*G(2)
    if l:
        sl = Q(1,2*d)
        out -= l*(positive(sl-s)*positive(U(2,0)-L(d+1))
                  +positive(s-sl)*positive(U(d+1,l-1)-L(2)))
    return s, out, U(d+1,l)


TICKS = 10**8


@cache
def edge_slopes(d,l,e,j):
    return (positive(U(d+1,l)-L(e+1)),positive(U(e+1,j)-L(d+1)))


def run(cap, envelope=False, stocked=False):
    assert not stocked or cap<10
    # Symbol 10 stands for all d>=10 cores with exactly one core slot.
    degrees = tuple(range(2,cap+1)) + ((10,) if stocked else ())
    families = defaultdict(list)
    profiles = []
    for d in degrees:
        for c in range(1, (1 if stocked and d==10 else d)+1):
            for ds in combinations_with_replacement(degrees,c):
                if 10 in ds and d==10:
                    continue
                a = (d,ds)
                profiles.append(a)
                for e in set(ds):
                    families[d,e].append(a)
    def variants(a):
        d,ds = a
        return (ds,tuple(0 if e==10 else e for e in ds)) if stocked and 10 in ds else (ds,)

    def lower_credit(a):
        d,ds = a
        if stocked and d==10:
            return Q(9,11)*Q(16694235,10**8)+Q(9,20)*G(2)
        return min(profile(d,row)[1] for row in variants(a))

    credit = {a: (TICKS*lower_credit(a))//1 for a in profiles}
    endpoints = {a: ( (0,TICKS//11) if stocked and a[0]==10 else
                      tuple((TICKS*profile(a[0],ds)[0])//1 for ds in variants(a)),
                      a[0]-len(a[1])) for a in profiles}

    @cache
    def slopes(d,l,e,j):
        if stocked and (d==10 or e==10):
            def u(d,l): return Q(1,10) if d==10 else U(d+1,l)
            def lwr(d): return Q(0) if d==10 else L(d+1)
            v = positive(u(d,l)-lwr(e)), positive(u(e,j)-lwr(d))
        else:
            v = edge_slopes(d,l,e,j)
        return tuple(-((-TICKS*x)//1) for x in v)

    def raw(d,l,x,b):
        ys,j = endpoints[b]
        f,g = slopes(d,l,b[0],j)
        # Each rounded difference can underestimate by at most one tick.
        return max((x-min(ys)+1)*f,(max(ys)-x+1)*g,0)

    def give(a,b):
        xs,l = endpoints[a]
        return (max(raw(a[0],l,x,b) for x in xs) + TICKS-1)//TICKS

    vals = {(a,e):credit[e,(a[0],)]-give(a,(e,(a[0],)))
            for a in profiles for e in set(a[1])}
    totals = {a:credit[a]+sum(vals[a,e] for e in a[1]) for a in profiles}
    worst = (0,None,None)
    if not envelope:
        for (a,e), v in vals.items():
            d = a[0]
            for b in families[e,d]:
                residual = totals[b]-vals[b,d]-give(a,b)-v
                if residual < worst[0]:
                    worst = residual,a,b
    else:
        ranges = {}
        for a,e in vals:
            xs,l = endpoints[a]
            assert len(xs)==1
            x = xs[0]
            key = a[0],e,l
            lo,hi = ranges.get(key,(x,x))
            ranges[key] = min(lo,x),max(hi,x)
        for (d,e,l),(lo,hi) in ranges.items():
            t = (e,(d,))
            for b in families[e,d]:
                if b == t:
                    continue
                # raw(b)-raw(t) is piecewise affine in the parent s tick.
                # Its slope changes only within one tick of either endpoint.
                points = {lo,hi}
                for a in (t,b):
                    (mid,),_ = endpoints[a]
                    points.update(x for x in (mid-1,mid,mid+1) if lo <= x <= hi)
                delta = max(raw(d,l,x,b)-raw(d,l,x,t) for x in points)
                residual = totals[b]-vals[b,d]-credit[t]-(delta+TICKS-1)//TICKS
                if residual < worst[0]:
                    worst = residual,(d,e,l),b
    root = min((t,a) for a,t in totals.items())
    print(f'cap={cap} stocked={stocked} profiles={len(profiles)} messages={len(vals)} '
          f'worst={worst} root={root}',flush=True)
    return worst[0] >= 0 and root[0] >= (22124577 if stocked else 27404273)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('cap',nargs='*',type=int,default=[7])
    p.add_argument('--pairs',action='store_true')
    p.add_argument('--stocked',action='store_true')
    args = p.parse_args()
    # A star with d>=4 has signed core-leaf debit <=0, since the
    # adverse centre-parent first-send bound is <=1/d <= L(2).
    assert L(2)>Q(1,4)
    assert H['rescue'](2,1)-H['upper'](2,0)>Q(16694235,10**8)
    assert 2*(Q(9,11)*Q(16694235,10**8)+Q(9,20)*G(2))-Q(1,110)>Q(22124577,TICKS)
    assert all(H['score'](H['spider'](d,1))[0]>0 for d in (2,3))
    print('PASS all stars (unbounded degree)',flush=True)
    for cap in args.cap:
        assert run(cap,not (args.pairs or args.stocked),args.stocked)
