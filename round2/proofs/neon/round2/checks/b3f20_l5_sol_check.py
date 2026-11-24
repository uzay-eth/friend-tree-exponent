#!/usr/bin/env python3
"""Cold v2 audit: independent layout, charge reduction and literal tree generator.

Run without -O. Requires Python 3.10+, LP64 little endian, and gcc for the
one-pass integer extrema kernel. No production source/checker is imported.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement, product
from math import comb, prod
import mmap
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
CERT = HERE.parent / 'cert_l5'
HDR = struct.Struct('<Qii4q12i4IQ')
PROFILE = (16000, 200, 20, 19, 12, 9, 11, 8, 12, 9, 10, 8, 7, 2, 3)
POT = (6432404, 273529, 27746004, 13880592, 354856788, 540925385,
       331822656, 131094480, 14215993)
EXPECTED = {
    2: (5656929, 76840591160, 184481683938, 146042051850, 184481683938),
    1: (5657570, 88296060680, 200835990548, 146044250013, 200701592409),
}


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def chunks(xs):
    return sorted(Counter(xs).items())


class Layout:
    def __init__(self, path):
        self.path = Path(path)
        with self.path.open('rb') as fp:
            (magic, self.mode, version, self.p, self.q, self.scale, sr,
             self.lf, self.kf, self.d0, self.kt, dp, kp, dt, kt,
             db, kb, ds, ks, self.dc, self.buckets, self.mask, reserved, rows
             ) = HDR.unpack(fp.read(HDR.size))
            require(magic == 0x4e454f4e42334648 and version == 2 and not reserved and
                    self.mode in (1, 2) and self.mask in (1, 2, 3), 'v2 header')
            require(all(0 < x <= 10**12 for x in (self.p, self.q, sr)) and
                    0 < self.scale <= 10**14, 'scale or exponent')
            require(4 <= self.dc <= self.d0 <= 40 and self.d0-1 <= self.kt <= 64 and
                    max(self.kt+1, self.d0) <= self.kf <= 500 and
                    self.d0 <= self.lf <= 30000 and self.buckets in (1, 2) and
                    (self.buckets == 1 or self.dc < self.d0), 'base cutoffs')
            require(all(self.dc <= d <= self.d0 and self.dc-1 <= k <= self.kt
                        for d, k in [(dp, kp), (dt, kt), (db, kb), (ds, ks)]),
                    'L4/L5 compatibility')
            self.profile = (self.lf, self.kf, self.d0, self.kt, dp, kp, dt, kt,
                            db, kb, ds, ks, self.dc, self.buckets, self.mask)
            self.sr = sr
            self.states = [(h, k, d-h-k) for h in (1, 0) for d in range(2, self.d0)
                           for k in range(d-h+1)]
            self.exact = {s: i for i, s in enumerate(self.states)}
            self.inf = len(self.states)
            self.nt = 2*(self.kt+2)
            controls = self.states + [None] + [(h, k, -1) for _ in range(3)
                       for h in (1, 0) for k in range(self.kt+2)]
            self.ng = len(controls)
            self.par = {g:i for i,g in enumerate(g for g,s in enumerate(controls)
                                                  if s is None or s[1] >= 1)}
            self.child = {g:i for i,g in enumerate(g for g,s in enumerate(controls)
                                                    if s is None or s[0] == 1)}
            self.hu = {g:i for i,g in enumerate(g for g,s in enumerate(controls)
                          if s is not None and s[0] == 1 and 1 <= s[1] <= self.kt)}
            self.cu = {g:i for i,g in enumerate(g for g,s in enumerate(controls)
                          if s is not None and 2 <= s[1] <= self.kt)}
            self.lookup = {}
            for tag, d, kmax, hs, kmin in [('p', dp, kp, (1,), 1),
                    ('t', dt, kt, (1,), 2), ('b', db, kb, (1,0), 2),
                    ('w', db, kb, (1,), 1), ('s', ds, ks, (1,0), 3)]:
                self.lookup[tag] = {s:i for i,s in enumerate((h,k,l) for h in hs
                       for k in range(kmin,kmax+1) for l in range(d-k-h))}

            # Build shape catalog by unordered forests of each total size.
            self.shapes, sizes = [], []
            self.domain, self.rank, self.length, self.offset = {}, {}, {}, {}
            nr = 0
            for size in range(1, 6):
                motifs = sorted(cs for a in range(size)
                                for cs in combinations_with_replacement(range(len(sizes)), a)
                                if sum(sizes[c] for c in cs) == size-1)
                for cs in motifs:
                    i = len(sizes)
                    self.shapes.append(cs)
                    sizes.append(size)
                    for h in (0, 1):
                        a = len(cs)
                        dom = [j for j,s in enumerate(self.states)
                               if s[0] == h and s[1] >= a and sum(s) < self.dc]
                        if a+h == 1:
                            dom += [('bucket', h, b) for b in range(self.buckets)]
                        self.domain[i,h] = dom
                        self.rank[i,h] = {g:r for r,g in enumerate(dom)}
                        count = len(dom)*prod(comb(self.length[c,1]+m-1,m)
                                              for c,m in chunks(cs))
                        self.length[i,h] = count
                        if size == 5:
                            self.offset[i,h] = nr
                            nr += count
            require((len(sizes), sizes.count(5), nr) == (17, 9, rows) and nr > 0,
                    'rooted shapes and row count')
            self.shape = {cs:i for i,cs in enumerate(self.shapes)}
            self.n = n = len(self.child)
            self.pairs = comb(n+1, 2)
            self.triples = comb(n+2, 3)
            np, hp, cp = map(len, (self.par, self.hu, self.cu))
            lp, lt, lb, lw, ls = (len(self.lookup[tag]) for tag in 'ptbws')
            self.stride = (self.kf+1)*(self.lf+1)
            self.aux = 2*self.inf+4*self.nt+self.nt*n+np*(self.kt+2)
            self.counts = (2*self.stride+2, self.ng**2, np*hp*n, cp*self.pairs,
                           np*lp*lp*n, np*lt*self.pairs, lb*lw*n*n, ls*self.triples,nr)
            self.checked = sum(self.counts)-self.ng**2+np*n
            self.bytes = HDR.size+8*(self.aux+sum(self.counts))
            require(os.fstat(fp.fileno()).st_size == self.bytes, 'certificate length')
            require(sys.byteorder == 'little', 'little-endian memoryview required')
            self.mapping = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
        raw = memoryview(self.mapping)[HDR.size:].cast('q')
        self.tables, self.starts = {}, {}
        i = self.aux
        for tag,size in zip(('f','g','hv','hc','p','t','b','s','5'), self.counts):
            self.starts[tag] = i
            self.tables[tag] = raw[i:i+size]
            i += size
        require(self.tables['g'][self.inf*(self.ng+1)] == 0, 'g(INF,INF)')
        if self.profile == PROFILE:
            require((self.q,self.scale,sr) == (10**7,10**11,10**9) and
                    self.p == EXPECTED[self.mode][0], 'production scales and exponent')
            require(self.counts == POT and (self.aux,self.bytes,self.checked) ==
                    (21687,11370156264,1421095995), 'production counts')
        print(f'{self.path}: mode={self.mode} P/Q={self.p}/{self.q} '
              f'counts={self.counts} aux={self.aux} bytes={self.bytes} '
              f'logged_states={self.checked}', flush=True)

    def tail(self, h, k):
        first = self.inf+1+(1-h)*(self.kt+2)+min(k,self.kt+1)
        return [self.inf] + [first+j*self.nt for j in range(3)]

    def find_f(self, s):
        h,k,l = s
        return (2*self.stride+1-h if k > self.kf else
                (1-h)*self.stride+k*(self.lf+1)+min(l,self.lf))

    def embed(self, s):
        d = sum(s)
        if d < self.d0:
            return ((self.exact[s], F(1)),)
        t = F(self.d0,d)
        return tuple((g,w) for j,g in enumerate(self.tail(*s[:2]))
                     if (w := comb(3,j)*t**j*(1-t)**(3-j)))

    def unordered(self, gs):
        order = sorted(self.child[g] for g in gs)
        n, m = self.n, len(order)
        r, prev = 0, 0
        for j,x in enumerate(order):
            r += comb(n-prev+m-j-1,m-j)-comb(n-x+m-j-1,m-j)
            prev = x
        return r

    def connected_value(self, adj, state, verts, controls=None):
        children = {u: sorted(v for v in adj[u] & verts if v > u) for u in verts}
        root = min(verts)
        cs = children[root]
        size = len(verts)
        tbl = self.tables
        if size == 1:
            return tbl['f'][self.find_f(state[root])], 'f'

        def weighted(vertices, ix):
            opts = ((controls[v],) if controls is not None and v in controls else
                    self.embed(state[v]) for v in vertices)
            return sum((prod(pair[1] for pair in ctrl)*ix(*(pair[0] for pair in ctrl))
                         for ctrl in product(*opts)), F(0))

        if size == 2:
            return weighted((root,cs[0]), lambda a,b: tbl['g'][a*self.ng+b]), 'g'
        if size == 3 and len(cs) == 1:
            u = cs[0]
            def row(a,b,c):
                return (tbl['hv'][(self.par[a]*len(self.hu)+self.hu[b])*self.n+self.child[c]]
                        if b in self.hu else 0)
            return weighted((root,u,children[u][0]),row), 'hv'
        if size == 3:
            return weighted((root,*cs),lambda a,b,c: tbl['hc'][self.cu[a]*self.pairs+
                     self.unordered((b,c))] if a in self.cu else 0), 'hc'
        if size == 4:
            lut = self.lookup
            if len(cs) == 3:
                i = lut['s'].get(state[root])
                return (0 if i is None else weighted(cs,lambda *g: tbl['s'][i*self.triples+
                                                      self.unordered(g)])), 's'
            if len(cs) == 2:
                u = next(c for c in cs if children[c])
                w = next(c for c in cs if c != u)
                i,j = lut['b'].get(state[root]), lut['w'].get(state[u])
                return (0 if None in (i,j) else weighted((w,children[u][0]),lambda a,b:
                        tbl['b'][((i*len(lut['w'])+j)*self.n+self.child[a])*self.n+self.child[b]])), 'b'
            u = cs[0]
            if len(children[u]) == 2:
                i = lut['t'].get(state[u])
                return (0 if i is None else weighted((root,*children[u]),lambda a,b,c:
                        tbl['t'][(self.par[a]*len(lut['t'])+i)*self.pairs+self.unordered((b,c))])), 't'
            mid = children[u][0]
            i,j = lut['p'].get(state[u]), lut['p'].get(state[mid])
            return (0 if None in (i,j) else weighted((root,children[mid][0]),lambda a,b:
                    tbl['p'][((self.par[a]*len(lut['p'])+i)*len(lut['p'])+j)*self.n+self.child[b]])), 'p'

        if any(sum(state[u]) >= self.dc and len(adj[u] & verts) != 1 for u in verts):
            return 0, None

        def rank(u):
            arr = sorted(rank(v)[:2] for v in children[u])
            ids = tuple(i for i,_ in arr)
            i = self.shape[ids]
            h = state[u][0]
            if sum(state[u]) >= self.dc:
                h = int(u != root)
                b = int(self.buckets == 2 and sum(state[u]) >= self.d0)
                token = ('bucket', h, b)
            else:
                token = self.exact[state[u]]
            r = self.rank[i,h][token]
            start = 0
            for c,m in chunks(ids):
                r *= comb(self.length[c,1]+m-1,m)
                r += sum(comb(arr[start+j][1]+j,j+1) for j in range(m))
                start += m
            return i,r,h
        i,r,h = rank(root)
        return tbl['5'][self.offset[i,h]+r], i

    def close(self):
        self.tables.clear()
        self.mapping.close()


# A standalone, shape-agnostic scanner: sequential native-endian int64 input,
# each consecutive group of stride*width maps into the same width output bins.
# No ranks, cutoffs, arithmetic on coefficients, or certificate logic in C.
SCAN = r'''
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
int main(int argc, char **argv) {
    if (argc != 2) return 2;
    int fd = open(argv[1], O_RDONLY);
    struct stat st;
    if (fd < 0 || fstat(fd, &st)) return 2;
    const int64_t *data = mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (data == MAP_FAILED) return 2;
    size_t first, repeat, bins, width, skip;
    while (scanf("%zu %zu %zu %zu %zu", &first,&repeat,&bins,&width,&skip) == 5) {
        if (!repeat || !bins || !width || first+(repeat-1)*skip+bins*width >
            (size_t)st.st_size/8) return 2;
        int64_t *lo = malloc(bins*sizeof(*lo)), *hi = malloc(bins*sizeof(*hi));
        if (!lo || !hi) return 2;
        for (size_t b=0;b<bins;b++) lo[b]=INT64_MAX, hi[b]=INT64_MIN;
        for (size_t r=0;r<repeat;r++) {
            const int64_t *p=data+first+r*skip;
            for (size_t b=0;b<bins;b++) {
                int64_t l=lo[b], h=hi[b];
                for (size_t j=0;j<width;j++,p++) {
                    if (*p<l) l=*p;
                    if (*p>h) h=*p;
                }
                lo[b]=l; hi[b]=h;
            }
        }
        for (size_t b=0;b<bins;b++) printf("%ld %ld\n", (long)lo[b],(long)hi[b]);
        free(lo); free(hi);
    }
    return ferror(stdout) || !feof(stdin);
}
'''


def comparison(c, scanner):
    queries = []

    def request(tag, start, repeat, bins, width, skip=None):
        i = len(queries)
        queries.append((HDR.size//8+c.starts[tag]+start,repeat,bins,width,
                        bins*width if skip is None else skip))
        return i

    n, ntp = c.n, len(c.par)
    sizes = {tag:len(index) for tag,index in c.lookup.items()}
    hvc = request('hv',0,ntp*len(c.hu),n,1)
    pc = request('p',0,ntp*sizes['p']**2,n,1)
    ybc = [request('b',i*sizes['w']*n*n,sizes['w']*n,n,1)
           for i in range(sizes['b'])]
    hc = request('hc',0,1,len(c.cu),c.pairs)
    yt = request('t',0,ntp,sizes['t'],c.pairs)
    s3 = request('s',0,1,sizes['s'],c.triples)
    l5 = defaultdict(list)
    for (i,h),off in c.offset.items():
        dom = c.domain[i,h]
        if not c.length[i,h]:
            continue
        span = c.length[i,h]//len(dom)
        cs = c.shapes[i]
        if len(cs) > 1:
            for r,g in enumerate(dom):
                l5[g,i].append(request('5',off+r*span,1,1,span))
        else:
            dom2 = c.domain[cs[0],1]
            width = span//len(dom2)
            for r,g in enumerate(dom2):
                l5[g,i].append(request('5',off+r*width,len(dom),1,width,span))

    output = subprocess.run([str(scanner),str(c.path)], check=True, text=True,
                            input=''.join(' '.join(map(str,r))+'\n' for r in queries),
                            stdout=subprocess.PIPE).stdout.splitlines()
    it = iter(output)
    stats = [[tuple(map(int,next(it).split())) for _ in range(q[2])] for q in queries]
    require(next(it,None) is None, 'extrema kernel result length')

    def merge(ps, zero=False):
        if zero:
            ps = [(0,0),*ps]
        return min(a for a,b in ps), max(b for a,b in ps)

    def plus(ps):
        return tuple(map(sum,zip(*ps)))

    def scale(p,k):
        return p[0]*k,p[1]*k

    gm = c.tables['g']
    child_pieces = {}
    for g,j in c.child.items():
        gv = [gm[p*c.ng+g] for p in c.par]
        yv = merge([scale(stats[q][j],s[1]-1) for q,s in zip(ybc,c.lookup['b'])],True)
        child_pieces[g] = [(min(gv),max(gv)), merge([stats[hvc][j]],True),
                           merge([stats[pc][j]],True), yv]
    single = {g:plus(ps) for g,ps in child_pieces.items()}
    sep, joint, t_hc = {}, {}, {}
    for k in range(c.kt+2):
        gs = c.tail(1,k)
        sep[k] = plus([merge(list(ps)) for ps in zip(*(child_pieces[g] for g in gs))])
        joint[k] = merge([single[g] for g in gs])
        if 2 <= k <= c.kt:
            for h in (0,1):
                t_hc[h,k] = merge([stats[hc][c.cu[g]] for g in c.tail(h,k)[1:]],True)

    env = {}
    shapes = tuple(c.shapes)
    for s,g in c.exact.items():
        if sum(s) >= c.dc:
            continue
        h,k,_ = s
        terms,caps = [],[]
        for i in range(8,17):
            cs = shapes[i]
            count = (shape_cap(shapes,i,k,c.dc-2) if len(cs)>1 else
                     h*shape_cap(shapes,cs[0],k,c.dc-2))
            caps.append(count)
            terms.append(scale(merge([stats[q][0] for q in l5.get((g,i),[])],True),count))
        require(4*sum(caps) <= sum(s)*comb(sum(s)-1+4*(c.dc-2),3), 'L5 outer cap')
        env[g] = plus(terms)

    lo4, hi, jlo4, jhi = [], [], [], []
    for h in (0,1):
        for k in range(c.kf+2):
            j = min(c.kt+1,k)
            if k > c.kf:
                v = c.tables['f'][c.find_f((h,k,0))]
                a = b = v
            else:
                for l in range(max(0,2-k-h),c.d0-k-h):
                    s = h,k,l
                    g = c.exact[s]
                    v = c.tables['f'][c.find_f(s)]
                    ps = [(v,v),env.get(g,(0,0))]
                    if h:
                        ps.append(single[g])
                    if g in c.cu:
                        ps.append(scale(stats[hc][c.cu[g]],comb(k,2)))
                    for tag,q,arity in [('t',yt,2),('s',s3,3)]:
                        if s in c.lookup[tag]:
                            ps.append(scale(stats[q][c.lookup[tag][s]],comb(k,arity)))
                    a,b = plus(ps)
                    lo4.append(4*a+(k-h)*c.scale)
                    hi.append(b)
                    jlo4.append(lo4[-1]); jhi.append(b)
                start = c.find_f((h,k,max(0,c.d0-k-h)))
                block = c.tables['f'][start:c.find_f((h,k,c.lf))+1]
                a,b = min(block),max(block)
            x,y = scale(t_hc.get((h,k),(0,0)),comb(k,2))
            l,u = sep[j] if h else (0,0)
            jl,ju = joint[j] if h else (0,0)
            lo4.append(4*(a+x+l)+(k-h)*c.scale); hi.append(b+y+u)
            jlo4.append(4*(a+x+jl)+(k-h)*c.scale); jhi.append(b+y+ju)
    phi = c.tables['f'][c.find_f((0,0,2))]
    result = min(lo4),max(hi),phi,max(jhi)
    print(f'4Sc0={result[0]} SC0={result[1]} SPhi(T3)={phi} '
          f'joint_4Sc0={min(jlo4)} joint_SC0={result[3]}', flush=True)
    if c.profile == PROFILE:
        require(result == EXPECTED[c.mode][1:], 'production comparison integers')
    return result


def adjacency(par):
    adj = [set() for _ in par]
    for v,p in enumerate(par[1:],1):
        adj[p].add(v)
        adj[v].add(p)
    return adj


@cache
def shape_cap(shapes, i, k, forwards):
    result = 1
    for j,m in chunks(shapes[i]):
        if m > k:
            return 0
        result *= comb(k,m)*shape_cap(shapes,j,forwards,forwards)**m
        k -= m
    return result


def potential(c, adj):
    core = frozenset(u for u,nb in enumerate(adj) if len(nb) >= 2)
    state = {u: (int(u != 0), len(adj[u] & core)-int(u != 0), len(adj[u]-core))
             for u in core}
    require(0 in core, 'physical root must be core')
    level = {frozenset((u,)) for u in core}
    value = F(0)
    terms = {}
    charged = Counter()
    for n in range(1,6):
        for subset in level:
            val,tag = c.connected_value(adj,state,subset)
            value += val
            if n == 5 and tag is not None:
                terms[subset] = val,tag
                apex = min(subset)
                ch = adj[apex] & subset
                anchor = apex if len(ch) > 1 else next(iter(ch))
                charged[anchor,tag] += 1
        if n < 5:
            level = {s | {u} for s in level for v in s for u in (adj[v] & core)-s}
            if n == 4:
                level = {s for s in level if all(len(adj[u]) < c.dc or
                          len(adj[u] & s) == 1 for u in s)}
    shapes = tuple(c.shapes)
    for (anchor,i),count in charged.items():
        h,k,_ = state[anchor]
        cs = shapes[i]
        cap = (h*shape_cap(shapes,cs[0],k,c.dc-2) if len(cs) == 1 else
               shape_cap(shapes,i,k,c.dc-2))
        require(sum(state[anchor]) < c.dc and count <= cap, 'observed per-shape charge cap')
    return value,terms,state


def check_tree(c, par, bounds, crossing_cap=None):
    adj = adjacency(par)
    phi,old,states = potential(c,adj)
    require(bounds[0]*len(states) <= 4*phi <= 4*bounds[1]*len(states), 'comparison on tree')
    drift = F(0)
    known_leaves = set()
    coverage = Counter()
    for u,nb in enumerate(adj):
        if len(nb) == 1:
            parent, = nb
            if parent in known_leaves:
                continue
            known_leaves.add(parent)
            mult = sum(len(adj[v]) == 1 for v in adj[parent])
        else:
            mult = 1
        rate = sum(F(1,len(adj[v])) for v in nb)
        after = par+(u,)
        newadj = adjacency(after)
        val,new,_ = potential(c,newadj)
        drift += mult*rate*(val-phi)
        if len(nb) == 1:
            # Enumerate pre-event four-cores via undirected expansion, then
            # extend only through the unique promoted leaf.
            expected = set()
            if len(adj[parent]) < c.dc:
                expected = {frozenset((parent,))}
                for _ in range(3):
                    expected = {s | {x} for s in expected for w in s
                                for x in (adj[w] & states.keys())-s}
                expected = {s | {u} for s in expected if all(len(adj[x]) < c.dc or
                              len(adj[x] & s) == 1 for x in s)}
            require(new.keys()-old.keys() == expected, 'newborn ownership')
            coverage['birth'] += len(expected)
            coverage['large_promoter'] += len(adj[parent]) >= c.dc
        elif c.buckets == 2 and len(nb) == c.d0-1:
            cap = c.d0-1 if crossing_cap is None else crossing_cap
            for s,(v,i) in old.items():
                if u not in s:
                    continue
                delta = new[s][0]-v
                require(min(0,cap*delta) <= rate*delta <= max(0,cap*delta),
                        'BIG_INNER crossing corner')
                coverage['nonzero_cross'] += bool(delta)
                coverage['apex_cross' if u == min(s) else 'child_cross'] += bool(delta)
                coverage['nonroot_apex_cross'] += bool(delta) and u == min(s) and u != 0
                coverage['shape_'+str(i)] += bool(delta)
    residual = drift-F(c.p,c.q)*phi
    require(residual <= 0 if c.mode == 1 else residual >= 0,
            f'concrete drift failure: {residual} on parents {par}')
    return coverage


def case_list(c):
    def decorate(par, ds=None):
        ds = ds or {}
        adj = adjacency(par)
        return par+tuple(u for u,nb in enumerate(adj) for _ in
                         range(ds.get(u,max(2,len(nb)))-len(nb)))

    yield 'T3', (-1,0,0)
    for d in (c.d0-1,c.d0,c.d0+1,109):
        yield f'star-{d}',(-1,)+(0,)*d
    for label, par, endpoints in [
        ('long-path',(-1,0,1,2,3,4,5,6),(0,7)),
        ('fork',(-1,0,0,1,1,3,4,2),(5,6)),
        ('sub-star',(-1,0,1,1,1,1,2,6),(5,7)),
        ('nonroot-apex',(-1,0,1,2,3,4,5,6),(1,7)),
    ]:
        for d in (c.dc-1,c.dc,c.d0-1,c.d0,c.d0+1,109):
            yield f'{label}-{d}',decorate(par,{u:d for u in endpoints})
    for root in (True,False):
        h = int(not root)
        p = (-1,) if root else (-1,0)
        p += (h,)*(c.kt+2)
        p += (h+1,c.kt+h+2,c.kt+h+3)  # One deep branch through the hub.
        yield f'overflow-k{c.kt+2}-h{h}',decorate(p)
    # Core root with all physical neighbors in the core at D=Dc-1.
    yield 'exact-apex-L0',decorate((-1,)+(0,)*(c.dc-1))
    # Hit each of the nine serialized shapes even if its particular crossing
    # jump happens to be zero. Every physical leaf is attached after the core.
    for i in range(8,17):
        p = []

        def emit(i,parent):
            u = len(p)
            p.append(parent)
            for child in c.shapes[i]:
                emit(child,u)

        emit(i,-1)
        endpoints = {u:c.d0-1 for u,nb in enumerate(adjacency(p)) if len(nb) == 1}
        yield f'rooted-shape-{i}',decorate(tuple(p),endpoints)
    for d in (c.lf-1,c.lf,c.lf+1):
        yield f'large-L-{d}',decorate((-1,0,1,2,3),{0:d+1})
    rng = random.Random(713893)
    for n in (8,12,16,20,24,29,35,42):
        par = [-1,0,0]
        for j in range(3,n):
            par.append(rng.choice(sorted(adjacency(par)[rng.randrange(j)])))
        yield f'friend-{n}',tuple(par)


def trees(c,bounds):
    totals = Counter()
    cases = list(case_list(c))
    require(len(set(t for _,t in cases)) == len(cases), 'duplicate trees')
    for name,t in cases:
        try:
            totals.update(check_tree(c,t,bounds))
        except ValueError as err:
            raise ValueError(f'{name}: {err}') from err
    require(all(totals[k] for k in ('birth','nonzero_cross','large_promoter',
                'apex_cross','child_cross','nonroot_apex_cross')), 'missing boundary coverage')
    print(f'{len(cases)}/{len(cases)} literal Fraction tree drifts passed: '
          f'{dict(sorted(totals.items()))}',flush=True)


def transfers(c):
    """Test the physical zero-sum edge identity and exact extrema of single-spectator fibers."""
    samples = Counter()
    for name,par in case_list(c):
        if not name.startswith(('long-path-', 'fork-', 'sub-star-', 'nonroot-apex-')):
            continue
        adj = adjacency(par)
        _,motifs,state = potential(c,adj)
        core = frozenset(state)
        # Independent 4-sets obtained by all edge-frontier expansions.
        level = {frozenset((u,)) for u in core}
        for _ in range(3):
            level = {s | {x} for s in level for z in s for x in (adj[z]&core)-s}
        moved = {}
        for m in level:
            degrees = {u:len(adj[u]&m) for u in m}
            if any(degrees[u] > 1 and sum(state[u]) >= c.dc for u in m):
                continue
            for u in m:
                r = degrees[u]
                if sum(state[u]) >= c.dc or not c.mask & (1 if r > 1 else 2):
                    continue
                h,k,l = state[u]
                next_state = {**state,u:(h,k,l+1)}
                jump = (c.connected_value(adj,next_state,m)[0]-
                        c.connected_value(adj,state,m)[0])
                ext = (adj[u]&core)-m
                rates = {z:F(1,sum(state[z])) for z in ext}
                clamped = sum(min(r,F(1,c.dc)) for r in rates.values())
                residuals = {z:r-F(1,c.dc) for z,r in rates.items() if r > F(1,c.dc)}
                require(clamped+sum(residuals.values()) == sum(rates.values()) and
                        0 <= clamped <= F(len(ext),c.dc), 'external split')
                for z,r in residuals.items():
                    moved[m,z] = r*jump
                spectators = [w for w in m if sum(state[w]) >= c.dc]
                if residuals and len(spectators) == 1:
                    w = spectators[0]
                    role = 'apex' if w == min(m) else 'child'
                    tail = sum(state[w]) >= c.d0
                    table_domain = c.par if role == 'apex' else c.child
                    fiber = [g for g in table_domain if (g >= c.inf and (tail or c.buckets == 1)) or
                             (g < c.inf and (not tail or c.buckets == 1) and sum(c.states[g]) >= c.dc)]
                    def gain(g):
                        override = {w:(g,F(1))}
                        return (c.connected_value(adj,next_state,m,override)[0]-
                                c.connected_value(adj,state,m,override)[0])
                    values = {g:gain(g) for g in fiber}
                    require(jump == sum(weight*values[g] for g,weight in c.embed(state[w])) and
                            min(values.values()) <= jump <= max(values.values()),'paired fiber hull')
                    samples[role+'_'+('tail' if tail else 'inner')] += 1

        received = {}
        for m5 in motifs:
            for z in m5:
                edge = adj[z] & m5
                if len(edge) != 1 or sum(state[z]) >= c.dc:
                    continue
                u, = edge
                m = m5-{z}
                if not c.mask & (1 if len(adj[u]&m)>1 else 2):
                    continue
                h,k,l = state[u]
                after = {**state,u:(h,k,l+1)}
                jump = c.connected_value(adj,after,m)[0]-c.connected_value(adj,state,m)[0]
                received[m,z] = (F(1,sum(state[z]))-F(1,c.dc))*jump
        require(moved == received, 'edge-wise L4/L5 transfer bijection and cancellation')
        samples['transfers'] += len(moved)
    if c.buckets == 2:
        require(all(samples[k] for k in ('apex_inner','apex_tail','child_inner','child_tail')),
                'fiber sample coverage')
    print(f'exact split and paired fiber check: {dict(sorted(samples.items()))}',flush=True)


def finite_sizes():
    # Propagate probability over labeled graphs with root fixed to the center
    # of T3. Aggregate attachment rates as a separate direct vertex-neighbor
    # loop. No tree is sampled or rejected.
    dist = {(-1,0,0):F(1)}
    lo,hi = F(8113447325,10248982441),F(146044250013,22074015170)
    a,b = F(5656929,10**7),F(565757,10**6)
    require(b-a == F(641,10**7), 'width')
    require(lo == F(EXPECTED[2][3],EXPECTED[2][2]) and
            hi == F(4*EXPECTED[1][3],EXPECTED[1][1]), 'product prefactors')
    for n in range(3,10):
        expectation = F(0)
        nxt = defaultdict(F)
        for par,mass in dist.items():
            adj = adjacency(par)
            expectation += mass*sum(len(nb)>=2 for nb in adj)
            if n < 9:
                for nb in adj:
                    for u in nb:
                        nxt[par+(u,)] += mass/(n*len(nb))
        require(sum(dist.values()) == 1 and lo <= expectation <= hi,'finite-size bound')
        print(f'n={n}: histories={len(dist)} E[C_n]={expectation} exact products bracket',flush=True)
        dist = nxt
        lo *= 1+a/n
        hi *= 1+b/n


def clean_run(args, expect=None):
    print('+',' '.join(map(str,args)),flush=True)
    out = subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    print(out.stdout,end='',flush=True)
    print(f'exit={out.returncode}',flush=True)
    require('runtime error:' not in out.stdout, 'UBSan failure')
    if expect is None:
        require(out.returncode == 0,'unexpected process failure')
    else:
        require(out.returncode == 1 and expect in out.stdout,'expected clean negative control failure')
    return out.stdout


def controls(scratch):
    flags = ['gcc','-std=c11','-fopenmp','-Wall','-Wextra','-Werror']
    builder = scratch / 'sol_b3f20_l5_build'
    san = scratch / 'sol_b3f20_l5_ubsan'
    planted = scratch/'wrong_bucket'
    src = CERT/'cert_b3f_l5_verify_gmp.c'
    sanitizer = ['-O2','-g','-fsanitize=undefined','-fno-sanitize-recover=all']
    clean_run([*flags,'-O3','-ffp-contract=off',CERT/'cert_b3f_l5_build.c','-lm','-o',builder])
    clean_run([*flags,*sanitizer,src,'-lgmp','-o',san])
    common = '--D0 10 --KT 9 --LF 4000 --KF 80 --DP4 8 --KP4 6 --DYT 7 --KYT 5 --DYB 8 --KYB 6 --DS3 6 --KS3 5 --Dc 6 --buckets 2 --l5mask 3 --Q 10000000 --cushion 5e-7 --nit 22 --threads 1'.split()
    try:
        for side in (2,1):
            path = scratch/f'small-{side}.bin'
            clean_run([builder,'--mode',str(side),*common,'--out',path])
            c = Layout(path)
            clean_run([san,path,'--threads','1','--tree-check','30'])
            scanner = scratch/'extrema'
            bounds = comparison(c,scanner)
            trees(c,bounds)
            clean_run([san,path,'--threads','1','--tree-check','0','--override-pq',
                       str(c.p+(20 if side == 2 else -20)),str(c.q)],'VERDICT: FAIL')
            if side == 1:
                bad = scratch/'bad.bin'
                shutil.copyfile(path,bad)
                with bad.open('r+b') as out:
                    out.seek(92)
                    out.write(struct.pack('<i',4))
                clean_run([san,bad,'--threads','1','--tree-check','0'],'invalid header (err=-21)')
                with bad.open('r+b') as out:
                    out.seek(92)
                    out.write(struct.pack('<i',5))
                    out.seek(-8,os.SEEK_END)
                    out.write(struct.pack('<q',10**15+1))
                clean_run([san,bad,'--threads','1','--tree-check','0'],'exceeds verified native bound')
                with bad.open('r+b') as out:
                    out.seek(-8,os.SEEK_END)
                    out.write(struct.pack('<q',-10**15))
                clean_run([san,bad,'--threads','1','--tree-check','0'],'VERDICT: FAIL')
                clean_run([*flags,*sanitizer,'-DL5_PLANT_BUCKET',src,'-lgmp','-o',planted])
                clean_run([planted,path,'--threads','1','--tree-check','0'],'VERDICT: FAIL')
                try:
                    for _,case in case_list(c):
                        check_tree(c,case,bounds,crossing_cap=0)
                except ValueError as err:
                    require(str(err) == 'BIG_INNER crossing corner','wrong crossing failure')
                    print('Zero crossing cap rejected by literal tree audit:',err,flush=True)
                else:
                    raise ValueError('zero crossing cap survived')
            c.close()
    finally:
        for p in (builder,san):
            p.unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('certificates',type=Path,nargs='*')
    ap.add_argument('--controls',action='store_true')
    ap.add_argument('--sizes-only',action='store_true')
    args = ap.parse_args()
    with tempfile.TemporaryDirectory(dir='.',prefix='sol_b3f20_l5_audit_') as tmp:
        scratch = Path(tmp)
        if (args.certificates and not args.sizes_only) or args.controls:
            subprocess.run(['gcc','-x','c','-','-std=c11','-O3','-Wall','-Wextra','-Werror',
                            '-o',str(scratch/'extrema')],input=SCAN,text=True,check=True)
        for path in args.certificates:
            c = Layout(path)
            if not args.sizes_only:
                trees(c,comparison(c,scratch/'extrema'))
                transfers(c)
            c.close()
        finite_sizes()
        if args.controls:
            controls(scratch)


if __name__ == '__main__':
    main()
