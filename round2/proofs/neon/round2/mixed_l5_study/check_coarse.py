#!/usr/bin/env python3
"""Exact source, gain, bud and global drift checks on concrete rooted trees."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from itertools import product
from random import Random
from zlib import crc32

D5 = DC = 6
D0, KT = 12, 11
INF = (-3, 0, 0)
BIG = (-1, 0, 0)
ZERO = Q(0)


def adjacencies(par):
    adj = [set() for _ in par]
    for i, p in enumerate(par[1:], 1):
        adj[i].add(p)
        adj[p].add(i)
    return adj


def motifs(adj):
    s = {(i,) for i in range(len(adj))}
    for n in range(2, 6):
        s = {tuple(sorted(m + (v,))) for m in s for u in m for v in adj[u] if v not in m}
        if n >= 4:
            yield s


def structure(m, par):
    return tuple(m.index(par[u]) if par[u] in m else -1 for u in m)


def degrees(p):
    return tuple(p.count(i) + bool(i) for i in range(len(p)))


@cache
def factor(s, d, apex):
    return Q(crc32(repr((s,d,apex)).encode()) % 37 - 18, 17)


@cache
def row(p, st):
    if len(p) == 4:
        return prod(factor(s,d,i==0) for i,(s,d) in enumerate(zip(st,degrees(p))))
    def tree(i):
        return st[i], tuple(sorted(tree(u) for u in range(i+1, len(p)) if p[u] == i))
    return Q(crc32(repr(tree(0)).encode()) % 199 - 99, 37)


def embed(s):
    l, k, h = s
    d = sum(s)
    if d < D0:
        return [(s, Q(1))]
    t = Q(D0, d)
    f = 1-t
    return [(INF, f**3), ((-2, min(k, KT+1), h), 3*t*f*f),
            ((-3, min(k, KT+1), h), 3*t*t*f), ((-4, min(k, KT+1), h), t**3)]


def active4(p, st):
    deg = degrees(p)
    if p == (-1, 0, 1, 2):
        cut = 10, 8
    elif p.count(1) == 2:
        cut = 9, 7
    elif p.count(0) == 2:
        cut = 10, 8
    else:
        cut = 8, 6
    return all(sum(s) < cut[0] and s[1] <= cut[1] for j,s in enumerate(st) if deg[j] > 1)


@cache
def q4(p, st):
    if not active4(p, st):
        return ZERO
    return prod(sum(w*factor(v,d,i==0) for v,w in embed(s))
                for i,(s,d) in enumerate(zip(st,degrees(p))))


def prod(seq):
    x = Q(1)
    for t in seq:
        x *= t
    return x


def project(p, st):
    a = tuple(s if sum(s) < DC else BIG for s in st)
    if any(s == BIG and d != 1 for s,d in zip(a,degrees(p))):
        return None
    return a


def q5(p, st):
    a = project(p, st)
    return row(p, a) if a is not None else ZERO


def event(st, u, dl, dk):
    l,k,h = st[u]
    return st[:u]+((l+dl,k+dk,h),)+st[u+1:]


@cache
def fiber(is_apex):
    # All physical inner BIG controls plus all appropriate tail controls. For
    # tail controls include infeasible corners too (a conservative superset).
    states = {INF}
    for h in ((0,1) if is_apex else (1,)):
        for k in range(int(is_apex), KT+2):
            states.update((l,k,h) for d in range(DC,D0) if (l:=d-k-h)>=0)
            states.update((i,k,h) for i in (-2,-3,-4))
    return tuple(states)


@cache
def extrema(p, a, u):
    if BIG not in a:
        x = q4(p, event(a,u,1,0))-q4(p,a)
        return x,x
    # Eligibility of the two preimage rows is decided by interiors only;
    # the enumerated labels replace endpoints only.
    up = event(a,u,1,0)
    before_ok, after_ok = active4(p,a), active4(p,up)
    ds = degrees(p)
    c = int(after_ok)*factor(up[u],ds[u],u==0)-int(before_ok)*factor(a[u],ds[u],u==0)
    terms = []
    for i,s in enumerate(a):
        if i==u:
            continue
        if s!=BIG:
            c *= factor(s,ds[i],i==0)
        else:
            terms.append(factor_bounds(ds[i],i==0))
    vals = [c*prod(v) for v in product(*terms)]
    return min(vals), max(vals)


@cache
def factor_bounds(d, apex):
    vals = [factor(s,d,apex) for s in fiber(apex)]
    return min(vals),max(vals)


def supported4(p, st, u, mask):
    return (bool(mask & (1 if degrees(p)[u]>1 else 2))
            and active4(p, st) and sum(st[u]) < DC and project(p, st) is not None)


def check(par, leaves, mask, counts):
    adj = adjacencies(par)
    st = tuple((leaves[i], len(a)-bool(i), int(bool(i))) for i,a in enumerate(adj))
    assert all(sum(s)>=2 and s[0]>=0 for s in st)
    ms, ns = motifs(adj)
    sources, reversed_sources = {}, {}
    direct = rem = lo = hi = ZERO
    buds = ZERO
    born = {}
    old_by_bud = {}
    for m in ms:
        p = structure(m,par)
        s = tuple(st[u] for u in m)
        v = q4(p,s)
        for i,u in enumerate(m):
            delta = q4(p,event(s,i,1,0))-v
            for z in adj[u].difference(m):
                r = Q(1,sum(st[z]))
                retained = min(r,Q(1,D5)) if supported4(p,s,i,mask) else r
                x = r-retained
                direct += r*delta
                rem += retained*delta
                counts['credits'] += 1
                if x:
                    sources[m,u,z] = x*delta
                    counts['mixed_credits'] += project(p,s).count(BIG)>0
                    counts['gain_crosses_dc'] += sum(st[u])+1==DC
            if st[u][0]:
                a = event(s,i,-1,1)+((1,0,1),)
                pa = p+(i,)
                if project(pa,a) is not None:
                    assert active4(p,s) and sum(st[u]) < DC and st[u][1] <= KT
                    n = m+(len(adj),)
                    assert (u,n) not in old_by_bud
                    old_by_bud[u,n] = Q(st[u][0],sum(st[u]))*q5(pa,a)
    l5actual = l5lo = l5hi = ZERO
    for n in ns:
        p = structure(n,par)
        s = tuple(st[u] for u in n)
        a = project(p,s)
        if a is None:
            # An unsupported large interior can never return to support.
            assert all(project(p,event(s,i,1,0)) is None and
                       project(p,event(s,i,-1,1)) is None for i in range(5))
            continue
        sublo = subhi = subactual = ZERO
        hits = 0
        for z in n:
            near = adj[z].intersection(n)
            if len(near) != 1 or sum(st[z]) >= D5:
                continue
            u, = near
            m = tuple(v for v in n if v != z)
            pp = structure(m,par)
            ss = tuple(st[v] for v in m)
            j = m.index(u)
            if not supported4(pp,ss,j,mask):
                continue
            x = Q(1,sum(st[z])) - Q(1,D5)
            exact = x*(q4(pp,event(ss,j,1,0))-q4(pp,ss))
            key = m,u,z
            assert key not in reversed_sources
            reversed_sources[key] = exact
            hits += 1
            f0,f1 = extrema(pp,project(pp,ss),j)
            assert x*f0 <= exact <= x*f1
            counts['strict_source'] += x*f0 < exact < x*f1
            sublo += x*f0
            subhi += x*f1
            subactual += exact
        assert sublo <= subactual <= subhi
        counts['multiple_deletions'] += hits>=2
        lo += sublo
        hi += subhi
        v = q5(p,s)
        for i,u in enumerate(n):
            g = q5(p,event(s,i,1,0))-v
            b = q5(p,event(s,i,-1,1))-v if st[u][0] else ZERO
            if a[i] == BIG:
                assert not g and not b
                counts['big_stays_big'] += 1
                continue
            r = Q(st[u][0]) + sum(Q(1,sum(st[w])) for w in adj[u] if w in n and a[n.index(w)]!=BIG)
            slack = Q(sum(w in n and a[n.index(w)]==BIG for w in adj[u]),DC) + Q(len(adj[u].difference(n)),2)
            rr = Q(st[u][0]) + sum(Q(1,sum(st[w])) for w in adj[u])
            br = Q(st[u][0],sum(st[u]))*b
            low, high = r*g+min(ZERO,slack*g)+br, r*g+max(ZERO,slack*g)+br
            exact = rr*g+br
            assert low<=exact<=high
            l5actual += exact
            l5lo += low
            l5hi += high
            counts['l5_rates'] += 1
    assert sources == reversed_sources and direct == rem + sum(sources.values())
    # Independent event enumeration: enumerate just the newborn's connected
    # motifs (expand from the newborn), without using the L4 preimage list.
    nb = len(adj)
    for u,s in enumerate(st):
        l,k,h = s
        if not l:
            continue
        nbr = adj + [{u}]
        nbr[u] = nbr[u] | {nb}
        frontier = {(nb,)}
        for _ in range(4):
            frontier = {tuple(sorted(n+(w,))) for n in frontier for v in n for w in nbr[v] if w not in n}
        post = event(st,u,-1,1)+((1,0,1),)
        for n in frontier:
            p = structure(n,par+[u])
            x = tuple(post[i] for i in n)
            if project(p,x) is not None:
                born[u,n] = Q(l,sum(s))*q5(p,x)
                counts['birth_with_big'] += BIG in project(p,x)
            elif sum(s)>=DC:
                counts['big_birth_excluded'] += 1
                counts['overflow_birth_excluded'] += k>KT
    assert born == old_by_bud
    buds = sum(born.values())
    assert rem+lo+l5lo+buds <= direct+l5actual+buds <= rem+hi+l5hi+buds
    counts.update(trees=1, transfers=len(sources), births=len(born))


if __name__ == '__main__':
    rng = Random(5162026)
    counts = Counter()
    for rep in range(320):
        n = rng.randrange(6,13)
        p = [-1]+[rng.randrange(i) for i in range(1,n)]
        # Explicit overflow endpoints, both apex and child; place a path
        # next to the overflow hub so it can be an endpoint of a motif.
        if rep%80==0:
            p = [-1,0,1,2,3,4] + [1]*12
        elif rep%80==1:
            p = [-1,0,1,2,3,4] + [0]*12
        a = adjacencies(p)
        leaves = [max(2-len(q),0)+rng.choice((0,0,0,1,1,2,3,5,9,20)) for q in a]
        if rep%80<2:
            leaves[1 if rep%80==0 else 0] = 3
        check(p,leaves,1+rep%3,counts)
    for k in ('strict_source','multiple_deletions','mixed_credits','birth_with_big','overflow_birth_excluded','gain_crosses_dc'):
        assert counts[k], k
    print('Exact coarse tree inequalities PASS:', ', '.join(f'{k}={v}' for k,v in sorted(counts.items())))
