#!/usr/bin/env python3
"""Cold Fraction oracle; does not import builder, verifier or shipped Python."""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
import json
from math import comb, prod
import mmap
from pathlib import Path
import random
import struct
import subprocess


FIELDS = ('magic mode P Q S Sr alpha p0 p1 p2 DF d2 K2 d3m K3m d3e K3e '
          'd4i K4i SP4 dS3 KS3 d4e K4e d5i K5i SP5 d5o K5o d5c K5c SB5 '
          'dS4 KS4 d5e K5e QK nk').split()
EDGES = dict(f=(), g=((0, 1),), h3=((0, 1), (1, 2)),
             p4=((0, 1), (1, 2), (2, 3)), s3=((0, 1), (0, 2), (0, 3)),
             p5=((0, 1), (1, 2), (2, 3), (3, 4)),
             b5=((0, 1), (0, 2), (0, 3), (3, 4)),
             s4=((0, 1), (0, 2), (0, 3), (0, 4)))
INS = dict(f=(), g=(), h3=(1,), p4=(1, 2), s3=(0,), p5=(1, 2, 3),
           b5=(0, 3), s4=(0,))
ENV = dict(e23=('g', 'h3', 1), e34p=('h3', 'p4', 1), e34s=('h3', 's3', 2),
           e45p=('p4', 'p5', 1), e45bp=('p4', 'b5', 2),
           e45bs=('s3', 'b5', 2), e45s4=('s3', 's4', 3))
INF = (2, 0, 0)
J = (0, 2, 1)


def adj(edges, n):
    a = [set() for _ in range(n)]
    for u, v in edges:
        a[u].add(v)
        a[v].add(u)
    return a


ADJ = {k: adj(e, len(e)+1) for k, e in EDGES.items()}


def shape(a):
    """Canonical roles by graph topology, not by stored corner labels."""
    tips = [u for u in range(len(a)) if len(a[u]) == 1]
    if len(a) == 1:
        return 'f', (0,)
    if len(tips) == 2:
        path = [tips[0]]
        while path[-1] != tips[1]:
            path.append(next(iter(a[path[-1]] - set(path[-2:]))))
        return {2: 'g', 3: 'h3', 4: 'p4', 5: 'p5'}[len(a)], tuple(path)
    hub = next(u for u in range(len(a)) if len(a[u]) >= 3)
    if len(tips) + 1 == len(a):
        return {3: 's3', 4: 's4'}[len(tips)], (hub, *tips)
    arm = next(u for u in a[hub] if len(a[u]) == 2)
    return 'b5', (hub, *sorted(a[hub] - {arm}), arm, next(iter(a[arm] - {hub})))


def pair(a, b):
    a, b = sorted((a, b))
    return comb(b+1, 2) + a


def inv(s):
    return F(0) if s == INF else F(1, s[1])


class Space:
    def __init__(self, d, k):
        self.d, self.k = d, k
        self.ex = [(0, D, c) for D in range(2, d) for c in range(1, min(D, k)+1)]
        self.states = (self.ex + [(3, D, k+1) for D in range(k+1, d)]
                       + [(1, d, c) for c in range(1, k+2)] + [INF])
        self.ids = {s: i for i, s in enumerate(self.states)}
        self.n = len(self.states)
        self.pairs = comb(self.n+1, 2)

    def project(self, s):
        kind, D, c = s
        if s == INF:
            return [(INF, F(1))]
        c = min(c, self.k+1)
        if D < self.d:
            return [((0 if c <= self.k else 3, D, c), F(1))]
        t = F(self.d, D)
        return [((1, self.d, c), t)] + ([(INF, 1-t)] if t < 1 else [])

    def up(self, s):
        if s[0] in (1, 2):
            return INF
        return self.project((0, s[1]+1, s[2]))[0][0]

    def bud(self, s):
        k, D, c = s
        if k in (2, 3) or (k == 0 and c == D) or c > self.k:
            return s
        return self.project((0, D, c+1))[0][0]


class Snapshot:
    def __init__(self, path):
        with Path(path).open('rb') as fp:
            self.raw = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
        vals = struct.unpack_from('<Q63q', self.raw)
        h = self.h = dict(zip(FIELDS, vals))
        assert h['magic'] == 0x4654325F55355631 and not any(vals[62:])
        self.xis = [F(n, h['QK']) for n in vals[38:38+h['nk']]]
        self.lam = F(h['P'], h['Q'])
        self.sign = 1 if h['mode'] == 1 else -1
        self.weight = lru_cache(maxsize=32_768)(self.weight)
        g2, g3, g4, g5 = [Space(h['d'+j], h['K'+j]) for j in ('2', '3e', '4e', '5e')]
        self.space = dict(g=g2, h3=g3, p4=g4, s3=g4, p5=g5, b5=g5, s4=g5)

        def inside(j, low):
            return [(0,D,c) for D in range(low, h['d'+j]) for c in range(low, min(D,h['K'+j])+1)]

        m3, m4, o3, m5, ob, cb, o4 = [inside(j,c) for j,c in
            [('3m',2), ('4i',2), ('S3',3), ('5i',2), ('5o',3), ('5c',2), ('S4',4)]]
        self.inners = dict(g=[()], h3=[(m,) for m in m3],
            p4=[(b,c) for b,c in product(m4,repeat=2) if b[1]+c[1] <= h['SP4']],
            s3=[(o,) for o in o3],
            p5=[(b,c,d) for b,c,d in product(m5,repeat=3)
                if max(b[1]+c[1], c[1]+d[1]) <= h['SP4'] and b[1]+c[1]+d[1] <= h['SP5']],
            b5=[(o,c) for o,c in product(ob,cb) if o[1]+c[1] <= h['SB5']],
            s4=[(o,) for o in o4])
        self.ids = {k:{s:i for i,s in enumerate(v)} for k,v in self.inners.items()}
        self.off = {}
        n = 0

        def alloc(k, size):
            nonlocal n
            self.off[k] = n
            n += size

        alloc('r',len(g2.ex))
        alloc('fy',len(g2.ex)*len(self.xis))
        alloc('f',comb(h['DF']+1,2)-3)
        self.blocks = {}
        for k, G in self.space.items():
            for I in self.inners[k]:
                if k in ('p4','p5') and I[0] > I[-1]:
                    continue
                a = len(EDGES[k])+1-len(I)
                size = (G.n**2 if k in ('p4','p5') and I[0] != I[-1] else
                        G.pairs*G.n if k == 'b5' else comb(G.n+a-1,a))
                self.blocks[k,I] = n
                n += size
        self.envsize = {}
        for k, (_, high, a) in ENV.items():
            G = self.space[high]
            size = G.n**a if k == 'e45bp' else comb(G.n+a-1,a)
            self.envsize[k] = size
            for col in ('m','r'):
                alloc(k+col, size*len(self.inners[high]))
        assert 512 + 8*n == len(self.raw), (512 + 8*n,len(self.raw))

    def i64(self, i):
        return struct.unpack_from('<q',self.raw,512+8*i)[0]

    def center(self, s):
        return F(self.i64(self.off['r']+self.space['g'].ids[s]),self.h['Sr'])

    def fy(self, s):
        return [self.i64(self.off['fy']+self.space['g'].ids[s]*len(self.xis)+i)
                for i in range(len(self.xis))]

    def f(self,D,c):
        h = self.h
        if D < h['DF']:
            return self.i64(self.off['f']+comb(D+1,2)-3+c)
        q = F(D-c,D)
        return h['p0']+h['p1']*q if self.sign == 1 else q*(h['p0']-h['p1']*q-F(h['p2'],D))

    @staticmethod
    def sorted_rank(indices, n):
        # Lex rank of combinations with replacement, via stars and bars.
        ids = sorted(indices)
        m = len(ids)
        prev = -1
        rank = 0
        for i, a in enumerate(ids):
            j = a+i
            rank += comb(n+m-prev-2,m-i)-comb(n+m-j-1,m-i)
            prev = j
        return rank

    def value(self,k,states):
        if k == 'f':
            return self.f(*states[0][1:])
        G = self.space[k]
        I = tuple(states[i] for i in INS[k])
        O = [G.ids[s] for i,s in enumerate(states) if i not in INS[k]]
        if k in ('p4','p5') and (I[0] > I[-1]):
            I, O = I[::-1], O[::-1]
        base = self.blocks.get((k,I))
        if base is None:
            return 0
        ix = (O[0]*G.n+O[1] if k in ('p4','p5') and I[0] != I[-1] else
              pair(*O[:2])*G.n+O[2] if k == 'b5' else pair(*O) if len(O) == 2 else
              self.sorted_rank(O,G.n))
        return self.i64(base+ix)

    def weight(self,k,s):
        if k == 'f':
            return self.value(k,s)
        return sum(prod(w for _,w in row)*self.value(k,[x for x,_ in row])
                   for row in product(*[[(x,F(1))] if i in INS[k] else self.space[k].project(x)
                                        for i,x in enumerate(s)]))

    def env(self,k,I,s,col):
        G = self.space[ENV[k][1]]
        O = [G.ids[x] for x in s]
        ix = (O[0]*G.n+O[1] if k == 'e45bp' else pair(*O) if len(O) == 2 else
              O[0] if len(O) == 1 else self.sorted_rank(O,G.n))
        return self.i64(self.off[k+col]+self.ids[ENV[k][1]][I]*self.envsize[k]+ix)

    def source(self,k,states,don):
        """Delete donor and locate its neighbor in the smaller motif."""
        a = ADJ[k]
        host = next(iter(a[don]))
        keep = [i for i in range(len(a)) if i != don]
        low, order = shape([{keep.index(j) for j in a[i] if j != don} for i in keep])
        order = tuple(keep[i] for i in order)
        i = order.index(host)
        if low == 'f':
            x = states[host]
            return (inv(states[don])-self.center(x))*(self.f(x[1]+1,x[2])-self.f(x[1],x[2]))
        fam = next(e for e,(l,h,_) in ENV.items() if (l,h) == (low,k))
        ordered = states
        if ((k in ('p4','p5') and don != 0)):
            ordered = states[::-1]
        I = tuple(ordered[i] for i in INS[k])
        O = [states[j] for j in order if j != host and order.index(j) not in INS[low]]
        if fam == 'e45bp':
            O = [states[3-don],states[4]]  # other short tip, then arm tip
        z = inv(states[don])-self.center(states[host])
        return self.sign*z*self.env(fam,I,O,'m')+abs(z)*self.env(fam,I,O,'r')


def replace(xs,i,x):
    return (*xs[:i],x,*xs[i+1:])


def extend(k,i):
    a = [x.copy() for x in ADJ[k]] + [{i}]
    a[i].add(len(a)-1)
    return shape(a)


def supported(snap,k,s):
    return k == 'f' or tuple(s[i] for i in INS[k]) in snap.ids[k]


def finite_row(snap,k,states):
    """Compute the signed cap by topology, projecting births and upper-bounding rates."""
    G = snap.space.get(k)
    sig = snap.sign
    a = ADJ[k]
    val = snap.value(k,states)
    ans = -sig*snap.lam*val

    def up(i):
        s = states[i]
        return (0,s[1]+1,s[2]) if k == 'f' or i in INS[k] else G.up(s)

    def changed(i,x):
        return snap.value(k,replace(states,i,x))

    if k == 'g':
        ans += sig*sum(snap.source(k,states,d) for d in range(2)
                       if states[1-d][0] == 0)
    elif k != 'f':
        ans += sum(snap.source(k,states,d) for d in range(len(a)) if len(a[d]) == 1)

    for i,s in enumerate(states):
        typ,D,c = s
        if typ == 2:
            continue
        delta = sig*(changed(i,up(i))-val)
        inner = typ in (0,3)
        # Kind-1 to inner donor credit is charged instead at the tail coordinate.
        linear = sum(inv(states[j]) for j in a[i] if states[j][0] in (0,3))
        if k == 'f':
            linear = 0
        if typ == 1:
            nb = next(iter(a[i]))
            gamma = max(0,(D-c+linear)/F(D+1)) if c <= G.k else 0
            ans += delta*(1 if delta > 0 else gamma)
            delta_b = sig*(changed(i,G.bud(s))-val)
            ans += delta_b*(1 if delta_b > 0 else max(0,F(D-c,D)))
            if states[nb][0] in (0,3):
                no_i = replace(states,i,INF)
                d0 = sig*(snap.value(k,replace(no_i,nb,up(nb)))-snap.value(k,no_i))
                d1 = sig*(changed(nb,up(nb))-val)
                ans += max(d0,d1)/F(D)
        elif inner:
            ans += linear*delta
            if typ == 3:
                ans += max(0,delta)*F(2*D-G.k-2,2)
                continue
            L, ext = D-c, c-len(a[i])
            ans += L*delta
            if ext < 0:
                raise AssertionError((k,states))
            nxt = extend(k,i) if len(a) < 5 else None
            if nxt:
                high, order = nxt
                post = (*states,J)
                covers = (k != 'f' or s in snap.space['g'].ids) and supported(snap,high,[post[j] for j in order])
            else:
                covers = False
            fy = k == 'g' and not covers and c >= 2
            if fy:
                xi = inv(states[1-i])
                xs, ys = snap.xis, snap.fy(s)
                j = next(j for j in range(len(xs)-1) if xs[j] <= xi <= xs[j+1])
                A = ys[j]+(ys[j+1]-ys[j])*(xi-xs[j])/(xs[j+1]-xs[j])
                ans += ext*(sig*A+max(x*delta-sig*y for x,y in zip(xs,ys)))
            elif covers:
                ans += ext*snap.center(s)*delta
            else:
                ans += F(ext,2)*max(0,delta)
            if L:
                dst = (0,D,c+1) if k == 'f' or i in INS[k] else G.bud(s)
                newborn = 0
                if nxt:
                    post = (*replace(states,i,(0,D,c+1)),J)
                    born = tuple(post[j] for j in order)
                    newborn = snap.weight(high,born) if supported(snap,high,born) else 0
                if k == 'f':
                    newborn += snap.f(2,1)
                ans += sig*F(L,D)*(changed(i,dst)-val+newborn)
    return ans


def samples(items, rng):
    return list(dict.fromkeys([items[0], items[-1], *rng.choices(items,k=min(7,len(items)))]))


def corners(G, rng):
    corners = [J, (0, G.d-1, min(G.d-1,G.k)), (1,G.d,1), (1,G.d,G.k+1), INF]
    if (3,G.d-1,G.k+1) in G.ids:
        corners.append((3,G.d-1,G.k+1))
    return [*corners,*rng.choices(G.states,k=10)]


def fill(k, I, O):
    inside, outside = iter(I), iter(O)
    return tuple(next(inside) if i in INS[k] else next(outside) for i in range(len(ADJ[k])))


def sampled_rows(s, rec, rng):
    caps = rec['verify']['caps']
    hits = Counter()
    h = s.h
    for D in sorted({2,3,*(h[j]+d for j in ('d2','d3e','d4e','d5e','DF') for d in (-1,0,1))}):
        if D >= h['DF']:
            continue
        for c in {0,1,2,min(D,h['K2']),min(D,h['K2']+1),D}:
            cap = caps['f_in' if D < h['d2'] and 1 <= c <= h['K2'] else 'f_tab']
            r = finite_row(s,'f',((0,D,c),))
            assert r <= cap, ('f',D,c,r,cap)
            hits['f'] += 1
    for k, G in s.space.items():
        E = corners(G,rng)
        arity = len(ADJ[k])-len(INS[k])
        rows = [(a,)*arity for a in E] + [(E[i], E[-1-i], *([J]*(arity-2))) for i in range(len(E))]
        rows.extend(tuple(rng.choices(E,k=arity)) for _ in range(29))
        for I in samples(s.inners[k],rng):
            for O in rows:
                r = finite_row(s,k,fill(k,I,O))
                cap = 0 if k == 'g' and all(x==INF for x in O) else caps[k]
                assert r <= cap, (k,I,O,r,cap)
                hits[k] += 1
    return dict(hits)


def sampled_envelopes(s, rng):
    counts = Counter()
    for e,(low,high,a) in ENV.items():
        FINE, COARSE = s.space[low], s.space[high]
        # Pick a canonical deletion that yields the source in forward order.
        keep, host = dict(e23=((0,1),1), e34p=((1,2,3),0),
            e34s=((2,0,3),1), e45p=((1,2,3,4),0), e45bp=((1,0,3,4),1),
            e45bs=((0,1,2,3),3), e45s4=((0,2,3,4),0))[e]
        C = corners(FINE,rng)
        rows = [(x,)*a for x in C] + [tuple(rng.choices(C,k=a)) for _ in range(39)]
        for I in samples(s.inners[high],rng):
            high_states = fill(high,I,[INF]*(len(ADJ[high])-len(I)))
            for O in rows:
                it = iter(O)
                low_states = tuple(high_states[j] if j in INS[high] else next(it) for j in keep)
                up = low_states[host]
                up = (0,up[1]+1,up[2]) if host in INS[low] else FINE.up(up)
                d = s.value(low,replace(low_states,host,up)) - s.value(low,low_states)
                m = r = F(0)
                for row in product(*(COARSE.project(x) for x in O)):
                    w, state = prod(p for _,p in row), [x for x,_ in row]
                    radius = s.env(e,I,state,'r')
                    assert radius >= 0
                    m += w*s.env(e,I,state,'m')
                    r += w*radius
                assert abs(d-m) <= r, (e,I,O,d,m,r)
                counts[e] += 1
    return dict(counts)


def tail(s, rec):
    h, hits = s.h, 0
    for D in (h['DF'], h['DF']+1, 10**6):
        for c in sorted({0,1,2,3,h['K2'],h['K2']+1,D//2,D-1,D}):
            w = s.f(D,c)
            for R in (D-c,F(2*D-c,2)):
                bud = 0 if D == c else F(D-c,D)*(s.f(D,c+1)-w+s.f(2,1)+s.weight('g',((0,D,c+1),J)))
                r = s.sign*(R*(s.f(D+1,c)-w)+bud-s.lam*w)
                assert r <= rec['cap_f_ana'], (D,c,R,r,rec['cap_f_ana'])
                hits += 1
            assert rec['tail']['w_min_ticks'] <= w <= rec['tail']['w_max_ticks']
            if c < 3:
                assert rec['tail'][f'w_c{c}_min_ticks'] <= w
    return hits


def connected_sets(a):
    groups = {(i,) for i in range(len(a))}
    for _ in range(5):
        yield from sorted(groups)
        groups = {tuple(sorted((*g,v))) for g in groups for u in g for v in a[u] if v not in g}


@lru_cache(maxsize=48)
def topology(edges):
    a = [set(nb) for nb in edges]
    result = []
    for g in connected_sets(a):
        k, order = shape([{g.index(v) for v in a[u] if v in g} for u in g])
        result.append((k,tuple((g[i],len(a[g[i]])) for i in order)))
    return result


def potential(s,a,D,cases):
    values = Counter()
    for k, pos in topology(tuple(tuple(sorted(nb)) for nb in a)):
        st = tuple((0,D[u],c) for u,c in pos)
        if supported(s,k,st):
            cases[k] += 1
            values[k,st] += 1
    return sum((n*s.weight(k,st) for (k,st),n in values.items()), F(0))


def adversaries(s, rec):
    h = s.h
    tests = []
    chain = [(0,1),(1,2),(2,3),(3,4)]
    for j in ('2','3e','4e','5e','F'):
        d = h['DF' if j == 'F' else 'd'+j]
        for p,q in product((d-1,d,d+1),repeat=2):
            tests.append(([(0,1)],[p,q]))
            tests.append((chain,[p,2,2,3,q]))
    # Adjacent hubs, both exact overflow or both degree-tail; keep core small by adding only leaves needed for c.
    for j in ('2','3e','4e','5e'):
        d,k = h['d'+j],h['K'+j]
        for c in (k,k+1):
            star = [(0,i) for i in range(1,c+1)]
            for D in sorted({c,max(c,d-1),max(c,d),h['DF']-1,h['DF']}):
                tests.append((star,[D]+[2+(i % 2)*(d-2) for i in range(c)]))
    for X,Y in ((h['d2'],h['d2']+1),(h['DF'],10**6),(h['DF']-1,h['DF'])):
        tests.append(([(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7)],
                      [X,Y,3,2,3,2,h['d4e'],h['d5e']]))
    # Mixed kind 0/3/1/INF outer arms (INF is reached by finite interpolation).
    k,d = h['K5e'],h['d5e']
    mixed = [(0,1),(0,2),(0,3),(0,4), *((1,i) for i in range(5,5+k))]
    tests.append((mixed,[4,k+1,2,d,h['DF']]+[2]*k))
    tests.extend((EDGES[k],[max(2,len(a)) for a in ADJ[k]]) for k in ADJ)
    cases = Counter()
    worst = None
    for edges,D in tests:
        a = adj(edges,len(D))
        assert all(D[u] >= max(2,len(a[u])) for u in range(len(a)))
        phi = potential(s,a,D,cases)
        drift = F(0)
        for u,x in enumerate(D):
            L = x-len(a[u])
            gain = L+sum(F(1,D[v]) for v in a[u])
            drift += gain*(potential(s,a,replace(D,u,x+1),cases)-phi)
            if L:
                aa = [y.copy() for y in a]+[{u}]
                aa[u].add(len(a))
                drift += F(L,x)*(potential(s,aa,(*D,2),cases)-phi)
        residual = s.sign*(drift-s.lam*phi)/len(a)
        assert rec['m_min']*len(a) <= phi <= rec['m_max']*len(a), (edges,D,phi)
        assert residual <= rec['err_per_core'], (edges,D,residual,rec['err_per_core'])
        worst = residual if worst is None else max(worst,residual)
    assert set(cases) == set(EDGES)
    return len(tests), str(worst)


def projection(s):
    hits = 0
    G = [s.space[k] for k in ('g','h3','p4','p5')]
    for fine,coarse in zip(G,G[1:]):
        for D in (2, fine.d-1, fine.d, fine.d+1, coarse.d-1, coarse.d, 10**6):
            for c in {1,min(D,fine.k),min(D,fine.k+1),min(D,coarse.k),min(D,coarse.k+1),D}:
                x = (0,D,c)
                f, direct = fine.project(x), dict(coarse.project(x))
                comp = Counter()
                for y,p in f:
                    for z,q in coarse.project(y):
                        comp[z] += p*q
                assert dict(comp) == direct and sum(p*inv(y) for y,p in f) == F(1,D)
                hits += 1
    return hits


def negative_controls(s, path, rec, verifier):
    cmd = [str(verifier), str(path), '--threads', '4']
    for mutation, why in [('ginf', 'g[INF,INF]'), ('env23', 'env23'), ('env45s4', 'env45s4')]:
        result = subprocess.run([*cmd, '--perturb', mutation], capture_output=True, text=True)
        assert result.returncode == 2 and why in result.stderr, (mutation, result)
        print(mutation, 'rejected:', result.stderr.strip(), flush=True)
    P, Q = s.h['P'], s.h['Q']
    shifted = s.lam-s.sign*F(1,10_000)
    assert Q % 10_000 == 0
    result = subprocess.run([*cmd,'--override-pq',str(P-s.sign*(Q//10_000)),str(Q)],
                            capture_output=True, text=True, check=True)
    v = json.loads(result.stdout)
    assert F(v['P'], v['Q']) == shifted
    # Ignoring the analytic cap gives a lower bound on the aggregate penalty.
    finite_err = max(0,max(v['caps']['f_in'],v['caps']['f_tab'])+max(0,v['caps']['g'])+v['higher_per_core'])
    assert finite_err > rec['err_per_core'] + 10**7
    best = shifted+s.sign*F(finite_err,rec['m_min'])
    assert s.sign*(best-F(rec['certified_exact'])) > 0
    print('tighter P/Q finite-only penalty and best possible endpoint:', finite_err, best, flush=True)


def binding_probe(lower, upper, verifier):
    """Schedule a pathname substitution around the real subprocess in the shipped wrapper."""
    import contextlib
    import hashlib
    import io
    import sys
    from tempfile import TemporaryDirectory
    from unittest.mock import patch
    import cert_check as wrapper  # Deliberately test the wrapper; independent oracles above do not import it.

    good = lower.read_bytes()
    snap = Snapshot(lower)
    damaged = bytearray(good)
    struct.pack_into('<q', damaged, 512 + 8*snap.off['f'], snap.f(2,0)+10**9)
    snap.raw.close()
    damaged = bytes(damaged)
    run_real = subprocess.run
    with TemporaryDirectory() as directory:
        p = Path(directory)/'probe.bin'
        p.write_bytes(damaged)
        direct = wrapper.combine_side(p,verifier,threads=1,expected_mode=2)
        print('direct damaged lower:',direct['certified_exact'],direct['err_per_core'])

        def interleave(cmd,*args,**kwargs):
            if str(p) not in cmd:
                return run_real(cmd,*args,**kwargs)
            p.write_bytes(good)
            try:
                return run_real(cmd,*args,**kwargs)
            finally:
                p.write_bytes(damaged)

        cli = ['cert_check.py','--verify-bin',str(verifier),'--lower',str(p),'--upper',str(upper),
               '--threads','1','--skip-trees']
        out, err = io.StringIO(), io.StringIO()
        with (patch.object(wrapper.subprocess,'run',interleave),patch.object(sys,'argv',cli),
              contextlib.redirect_stdout(out),contextlib.redirect_stderr(err)):
            try:
                wrapper.main()
            except SystemExit as stop:
                if stop.code == 2:
                    print('pathname substitution rejected',file=sys.__stdout__)
                    return 0
                raise
        receipt = json.loads(out.getvalue())
        snap = Snapshot(p)
        r = finite_row(snap,'f',((0,2,0),))
        snap.raw.close()
        claimed = receipt['lower']['err_per_core']
        assert receipt['lower_sha256'] == hashlib.sha256(damaged).hexdigest()
        print('reported lower:',receipt['lower']['certified_exact'], 'claims SHA256 of damaged input')
        print('damaged C=1 residual:',r,'reported cap:',claimed)
        return 1 if r > claimed else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('receipt',type=Path)
    ap.add_argument('lower',type=Path)
    ap.add_argument('upper',type=Path)
    ap.add_argument('--negative-controls',type=Path,metavar='VERIFIER')
    ap.add_argument('--binding-probe',type=Path,metavar='VERIFIER')
    args = ap.parse_args()
    if args.binding_probe:
        raise SystemExit(binding_probe(args.lower,args.upper,args.binding_probe))
    receipt = json.loads(args.receipt.read_text())
    for side in ('lower','upper'):
        snap = Snapshot(getattr(args,side))
        rec = receipt[side]
        assert snap.h['mode'] == rec['mode'] == (2 if side == 'lower' else 1)
        exact = snap.lam+snap.sign*F(rec['err_per_core'],rec['m_min'])
        assert exact == F(rec['certified_exact'])
        rng = random.Random(514)
        print(side, 'projection',projection(snap),'rows',sampled_rows(snap,rec,rng),
              'envelopes',sampled_envelopes(snap,rng),'tail',tail(snap,rec),flush=True)
        print(side, 'adversaries (count, worst signed residual/core):',adversaries(snap,rec),flush=True)
        if args.negative_controls:
            negative_controls(snap,getattr(args,side),rec,args.negative_controls)
        snap.raw.close()


if __name__ == '__main__':
    main()
