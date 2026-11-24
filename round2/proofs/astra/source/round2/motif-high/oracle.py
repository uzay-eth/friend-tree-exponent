"""Independent physical evaluator for FT2_4V01, with integer weights in snapshot units."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
import mmap
import struct
import sys

if not __debug__:
    raise RuntimeError("the independent physical tests require Python assertions")

FAMILIES = ("E", "V", "C", "P", "F", "Y", "T")
ARITY = dict(E=2, V=3, C=3, P=4, F=4, Y=4, T=4)
FIELDS = """DF D2 K2 D3 K3 DH KH DC KC D4 K4 DP1 DP2 DPS KP1 KP2
            DF1 DF2 DFS KF1 KF2 DY KY DT KT QK nk""".split()
HEADER = struct.Struct("<Q4i6q44i")
INF = None
NB = (1, 0)


def degree(s):
    l, k = s
    return l + k + 1


def inv(s):
    return Q(0) if s is INF else Q(1, degree(s))


def gain(s):
    return s[0] + 1, s[1]


def bud(s):
    assert s[0] > 0
    return s[0] - 1, s[1] + 1


class Level:
    def __init__(self, d0, kt):
        self.d0, self.kt = d0, kt
        self.internal = [(d-1-k, k) for d in range(2, d0) for k in range(d)]
        self.id = {s: i for i, s in enumerate(self.internal)}
        self.inf = len(self.internal)
        self.n = self.inf + kt + 3
        parent = [i for i in range(self.n)
                  if i == self.inf or (self.internal[i][1] if i < self.inf else i-self.inf-1) >= 1]
        self.parent = {v: i for i, v in enumerate(parent)}
        self.pairs = self.n*(self.n+1)//2
        self.tris = self.pairs*(self.n+2)//3

    @lru_cache(maxsize=8192)
    def mu(self, s):
        if s is INF:
            return ((self.inf, Q(1)),)
        if s in self.id:
            return ((self.id[s], Q(1)),)
        t = Q(self.d0, degree(s))
        v = self.inf + 1 + min(s[1], self.kt+1)
        return ((v, t), (self.inf, 1-t)) if t < 1 else ((v, t),)

    def pair(self, i, j):
        a, b = sorted((i, j))
        return self.pairs - (self.n-a)*(self.n-a+1)//2 + b-a

    def tri(self, i, j, k):
        a, b, c = sorted((i, j, k))
        left = self.n-a
        first = self.tris-left*(left+1)*(left+2)//6
        second = left*(left+1)//2 - (self.n-b)*(self.n-b+1)//2
        return first + second + c-b


class Certificate:
    def __init__(self, path):
        assert sys.byteorder == "little"
        with open(path, "rb") as f:
            self.mapping = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        h = HEADER.unpack_from(self.mapping)
        assert h[0] == 0x4654325F34563031 and h[1] in (1, 2)
        self.sign = 1 if h[1] == 1 else -1
        self.lam, self.S, self.Sr = Q(h[2], h[3]), h[5], h[6]
        self.par = dict(zip(FIELDS, h[11:38]))
        self.levels = {j: Level(self.par['D'+str(j)],self.par['K'+str(j)]) for j in (2,3,4)}
        h = self.par
        g2, g3, g4 = self.levels.values()
        def middle(dmax, kmin, kmax):
            return [(l,k) for k in range(kmin,kmax+1) for l in range(dmax) if l+k+1 < dmax]
        hv, hc = middle(h['DH'],1,h['KH']), middle(h['DC'],2,h['KC'])
        def pairs(centers, suffix):
            return [(u,v) for u in centers for v in hv
                    if degree(u) < h['D'+suffix+'1'] and degree(v) < h['D'+suffix+'2']
                    and degree(u)+degree(v) < h['D'+suffix+'S']
                    and u[1] <= h['K'+suffix+'1'] and v[1] <= h['K'+suffix+'2']]
        domains = dict(V=hv, C=hc, P=pairs(hv,'P'), F=pairs(hc,'F'),
                       Y=[u for u in hv if u in set(hc) and u[1] <= h['KY'] and degree(u) < h['DY']],
                       T=[u for u in hc if 3 <= u[1] <= h['KT'] and degree(u) < h['DT']])
        self.act = {f: {s:a for a,s in enumerate(dom)} for f,dom in domains.items()}
        self.data = memoryview(self.mapping)[HEADER.size:].cast('q')
        self.sections = {}
        cursor = 0
        def take(name, size):
            nonlocal cursor
            self.sections[name] = (cursor, size)
            a = self.data[cursor:cursor+size]
            cursor += size
            return a
        self.cent = [take(n, g2.inf) for n in ('rpar','rch','xpar','xch','xcc')]
        take('fy', 3*g2.inf*h['nk'])
        take('f', h['DF']*(h['DF']-1)//2-1)
        lengths = [g2.n**2, len(g3.parent)*len(hv)*g3.n, len(hc)*g3.pairs,
                   len(domains['P'])*len(g4.parent)*g4.n, len(domains['F'])*g4.n**2,
                   len(domains['Y'])*len(g4.parent)*g4.pairs, len(domains['T'])*g4.tris]
        self.tables = {f: take(f,n) for f,n in zip(FAMILIES,lengths)}
        for a,b in [('e23',len(hv)*g3.n), ('ep',len(domains['P'])*g4.n),
                    ('ef',len(domains['F'])*g4.n), ('eyv',len(domains['Y'])*len(g4.parent)*g4.n),
                    ('eyc',len(domains['Y'])*g4.n*g4.n)]:
            for j in range(2 if a.startswith('ey') else 4): take(a+str(j),b)
        assert cursor == len(self.data)

    @staticmethod
    def level(f):
        return ARITY[f]

    def active(self, f, s):
        if f == 'E': return 0
        if f in ('V','Y'): key = s[1]
        elif f in ('C','T'): key = s[0]
        elif f == 'P': key = (s[1],s[2])
        else: key = (s[0],s[2])
        return self.act[f].get(key, -1)

    def offset(self, f, a, x):
        g = self.levels[self.level(f)]
        if f == 'E': return x[0]*g.n+x[1]
        if f == 'V': return (g.parent[x[0]]*len(self.act[f])+a)*g.n+x[2]
        if f == 'C': return a*g.pairs+g.pair(x[1],x[2])
        if f == 'P': return (a*len(g.parent)+g.parent[x[0]])*g.n+x[3]
        if f == 'F': return (a*g.n+x[1])*g.n+x[3]
        if f == 'Y': return (a*len(g.parent)+g.parent[x[0]])*g.pairs+g.pair(x[2],x[3])
        return a*g.tris+g.tri(*x[1:])

    @lru_cache(maxsize=65536)
    def corners(self, f, s):
        g = self.levels[self.level(f)]
        res = []
        for items in product(*(g.mu(q) for q in s)):
            weights = Q(1)
            for _,w in items: weights *= w
            res.append((tuple(i for i,_ in items),weights))
        return tuple(res)

    @lru_cache(maxsize=131072)
    def weight(self, f, s):
        a = self.active(f,s)
        if a < 0: return Q(0)
        return sum((w*self.tables[f][self.offset(f,a,x)] for x,w in self.corners(f,s)),Q(0))

    def center(self, idx, s):
        return Q(self.cent[idx][self.levels[2].id[s]],self.Sr)


class Tree:
    # 0 is the inert apex. The other entries are physical h=1 vertices.
    def __init__(self, parents, leaves):
        assert parents[0] == -1 and len(parents) == len(leaves)
        self.p, self.l = tuple(parents), tuple(leaves)
        ch = [[] for _ in parents]
        for i in range(1,len(parents)):
            assert 0 <= parents[i] < i
            ch[parents[i]].append(i)
        assert len(ch[0]) == 1
        self.ch = ch
        self.s = (INF,)+tuple((leaves[i],len(ch[i])) for i in range(1,len(parents)))
        assert all(s[0] >= 0 and degree(s) >= 2 for s in self.s[1:])
        self.r = [Q(0)]+[s[0]+inv(self.s[self.p[i]])+sum(inv(self.s[j]) for j in ch[i])
                           for i,s in enumerate(self.s) if i]
        self.b = (Q(0),)+tuple(Q(s[0],degree(s)) for s in self.s[1:])

    def after(self, x, promotion):
        l = list(self.l)
        l[x] += -1 if promotion else 1
        return Tree(self.p+(x,),l+[1]) if promotion else Tree(self.p,l)

    def embeddings(self):
        for u,children in enumerate(self.ch):
            for a in children:
                yield 'E',(u,a)
                for b in self.ch[a]:
                    yield 'V',(u,a,b)
                    for c in self.ch[b]: yield 'P',(u,a,b,c)
                    for side in children:
                        if side != a: yield 'F',(u,side,a,b)
            for v,w in combinations(children,2):
                yield 'C',(u,v,w)
                if u: yield 'Y',(self.p[u],u,v,w)
            for a,b,c in combinations(children,3):
                yield 'T',(u,a,b,c)

    def active(self, cert):
        result = {}
        for f,verts in self.embeddings():
            states = tuple(self.s[v] for v in verts)
            if cert.active(f,states) < 0: continue
            support = frozenset(verts)
            assert support not in result
            result[support] = (f,verts,states,cert.weight(f,states))
        return result

    def potential(self, cert):
        return sum((w for f,v,s,w in self.active(cert).values()),Q(0))


def replace(s, i, val):
    return s[:i]+(val,)+s[i+1:]


def newborn_extensions(f, post_states, pos):
    a = post_states
    if f == 'E':
        return ('C',(a[0],a[1],NB)) if pos == 0 else ('V',(a[0],a[1],NB))
    if f == 'V':
        if pos == 0: return 'F',(a[0],NB,a[1],a[2])
        if pos == 1: return 'Y',(a[0],a[1],a[2],NB)
        return 'P',a+(NB,)
    if f == 'C':
        if pos == 0: return 'T',a+(NB,)
        return 'F',(a[0],a[3-pos],a[pos],NB)
    raise ValueError(f)


# Cut specifications are manual and independent of the outside-neighbour walk.
# donor, contact, lower family and its ordered retained coordinates, centering index.
CUTS = {
 'V': [(0,1,'E',(1,2),2),(2,1,'E',(0,1),3)],
 'C': [(1,0,'E',(0,2),4),(2,0,'E',(0,1),4)],
 'P': [(0,1,'V',(1,2,3),0),(3,2,'V',(0,1,2),1)],
 'F': [(1,0,'V',(0,2,3),1),(3,2,'C',(0,1,2),1)],
 'Y': [(0,1,'C',(1,2,3),0),(2,1,'V',(0,1,3),1),(3,1,'V',(0,1,2),1)],
 'T': [(1,0,'C',(0,2,3),1),(2,0,'C',(0,1,3),1),(3,0,'C',(0,1,2),1)]
}
