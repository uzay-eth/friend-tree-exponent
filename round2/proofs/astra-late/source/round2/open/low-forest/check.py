"""Independent labelled-arrow fixtures and rational tail constants."""

from collections import Counter
from fractions import Fraction as F
from itertools import product


def tree(parents):
    out = [set() for _ in range(len(parents)+1)]
    for v,u in enumerate(parents,1):
        out[v].add(u)
        out[u].add(v)
    return out


def moves(g, inside):
    for u,nu in enumerate(g):
        for v in nu & inside:
            out = [set(s) for s in g]+[{v}]
            out[v].add(len(g))
            yield out, F(1,len(nu))


def walks(g, inside, k):
    if not k:
        yield tuple(sorted(len(g[v] & inside) for v in inside)), F(1)
        return
    for out,r in moves(g,inside):
        for key,w in walks(out,inside|{len(g)},k-1):
            yield key,r*w


def short_shapes():
    star, path = (1,1,1,3), (1,1,2,2)
    for d in (4,6,20,100):
        # 0 is the exterior parent, arm is 1--2. Other children of 0
        # are not part of the tested fringe. Exterior updates do not enter
        # the leading coefficient for exactly two inside births.
        g = tree([0,1]+[0]*(d-1))
        assert len(g[0]) == d
        out = Counter()
        for key,r in walks(g,{1,2},2):
            out[key] += r
        a = 1+F(1,d)
        assert out[star] == a*(a+1)+F(3,4)
        assert out[path] == F(7,6)*a
        p = out[star]/sum(out.values())
        assert p > F(1,2)
        print(f'planted D_parent={d}: star={p}, path={1-p}')
    out = Counter()
    for key,r in walks(tree([0]),{0,1},2):
        out[key] += r
    assert out[star] == 4 and out[path] == 2
    assert F(33,47)>F(2,3)
    print('ordinary four-vertex RFT: star=2/3; silent-parent small-age star=33/47')


def split_fixtures():
    num = nontriv = 0
    for n in range(3,8):
        for parents in product(*(range(v) for v in range(1,n))):
            g = tree(parents)
            d = list(map(len,g))
            for N in range(2,6):
                low = {v for v in range(n) if d[v]<=N}
                raw = components(g,low)
                cores = {v for v in low if d[v]>=2}
                pieces = components(g,cores)
                restored = {frozenset(frozenset(p).union(*(set(w for w in g[v] if d[w]==1) for v in p))) for p in pieces}
                assert restored == {frozenset(p) for p in raw if p & cores}
                assert sum(len(p & cores) for p in raw) == len(cores)
                for p in restored:
                    assert all(d[w]>N for v in p for w in g[v]-p)
                for p in pieces:
                    # Reconstruct the highest surviving core from a genuine
                    # core orientation (first core is the root, not label 0).
                    root = next(v for v in range(n) if d[v]>=2)
                    parent = {root: None}
                    order = [root]
                    for u in order:
                        for v in g[u]:
                            if d[v]>=2 and v not in parent:
                                parent[v] = u
                                order.append(v)
                    tops = [v for v in p if parent[v] not in p]
                    assert len(tops)==1
                    top = tops[0]
                    assert parent[top] is None or d[parent[top]]>N
                    nontriv += 1
                num += 1
    print(f'forest fixtures: {num} threshold/tree pairs, {nontriv} surviving pieces, exact boundary/top PASS')


def components(g,vertices):
    left = set(vertices)
    out = []
    while left:
        v = left.pop()
        q = [v]
        for u in q:
            unseen = g[u] & left
            left -= unseen
            q.extend(unseen)
        out.append(set(q))
    return out


def tail_constants():
    # k=3r=16Nt, t>=4N^2. The analytic bound is <=5 exp(-.08k)
    # after the stated old exponent comparison; e^(-x)<=2^(-floor(x)).
    for N in (2,5,20,100):
        k = 3*((64*N**3+2)//3)
        t,j = F(k,16*N),k//3
        assert t >= 4*N*N
        # log N <= N/2 for N>=2; t>=4N^2 suffices for the union absorption.
        assert j >= 16*N*(F(87,200)+F(32,25)*N)
        assert F(87,200)/(16*N)+F(1,3)*(-F(1,4)-F(1,32)) <= -F(2,25)
        print(f'N={N} radius={k} bound <= 5 * 2^(-{(2*k)//25})')


if __name__ == '__main__':
    short_shapes()
    split_fixtures()
    tail_constants()
