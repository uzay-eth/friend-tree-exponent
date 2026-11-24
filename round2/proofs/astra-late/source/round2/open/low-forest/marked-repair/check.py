"""Exact controls of the alphabet witness, two censuses and the radius split."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product


def shape(bits):
    if len(bits) == 1:
        return bits[0]
    h = len(bits)//2
    return tuple(sorted((shape(bits[:h]), shape(bits[h:]))))


def graph(par):
    g = [set() for _ in range(len(par)+1)]
    for v,u in enumerate(par,1):
        assert 0 <= u < v
        g[u].add(v)
        g[v].add(u)
    return g


def binary_fixture(colors, N):
    # Vertex zero is the deleted parent; its first child is the skeleton top.
    par = [0]
    for _ in range(N):
        par.append(0)
    layer = [1]
    while len(layer) < len(colors):
        nxt = []
        for u in layer:
            nxt.extend((len(par)+1,len(par)+2))
            par.extend((u,u))
        layer = nxt
    for u,k in zip(layer,colors,strict=True):
        par.extend([u]*k)
    # Starting from 0--1, every later target has an existing adjacent sender.
    g = graph(par)
    assert len(g[0]) == N+1 and all(1 <= len(s) <= 5 for s in g[1:])
    return tuple(par)


def compare(par,N,J,B):
    g = graph(par)
    deg = list(map(len,g))
    core = {v for v,d in enumerate(deg) if d >= 2}
    low = {v for v,d in enumerate(deg) if d <= N}

    def marked(v,p):
        return (v in core,tuple(sorted(marked(w,v) for w in g[v]&low if w != p)))

    direct, weighted = Counter(),Counter()
    seen = set()
    count = dropped = dropped_weight = 0
    for v in sorted(low):
        if v in seen:
            continue
        comp = [v]
        seen.add(v)
        for u in comp:
            for w in g[u]&low-seen:
                comp.append(w)
                seen.add(w)
        K = len(set(comp)&core)
        if not K:
            continue
        count += 1
        if v == 0 or len(comp)>J or deg[par[v-1]]>B:
            dropped += 1
            dropped_weight += K
        else:
            mark = Q(deg[par[v-1]],N),marked(v,-1),K
            direct[mark] += 1
            weighted[mark] += K

    parent, survivor = Counter(),Counter()
    total_tops = 0
    for u in core-low:
        for v in g[u]&low&core:
            if v == 0 or par[v-1] != u:
                continue
            total_tops += 1
            if deg[u]>B:
                continue
            todo = [v]
            for w in todo:
                todo.extend(x for x in g[w]&low if x>w)
                if len(todo)>J:
                    break
            if len(todo)>J:
                continue
            K = len(set(todo)&core)
            mark = Q(deg[u],N),marked(v,u),K
            parent[mark] += 1
            survivor[mark] += K
    assert parent == direct and survivor == weighted
    assert total_tops + int(0 in core&low) == count
    assert sum(parent.values()) + dropped == count
    assert sum(survivor.values()) + dropped_weight == len(core&low)
    return dropped_weight>dropped


def main():
    expected = 4
    for h in range(4):
        L = 2**h
        assert expected >= 2**(L+1)
        if h<=2:
            exact = {shape(seq) for seq in product(range(1,5),repeat=L)}
            actual = set()
            for seq in product(range(1,5),repeat=L):
                g = graph(binary_fixture(seq,5))
                def view(v,p):
                    return tuple(sorted(view(w,v) for w in g[v] if w!=p))
                actual.add(view(1,0))
            assert len(exact)==len(actual)==expected
        expected = expected*(expected+1)//2
    print('PASS binary alphabet: exact orbit counts 4, 10, 55; symbolic lower bound and reachable realizations')

    trees = [(0,0)]
    count = strict = 0
    for n in range(3,8):
        for par in trees:
            for N,J,B in product((2,3,5),(1,3,8),(3,8)):
                strict += compare(par,N,J,B)
                count += 1
        trees = [par+(u,) for par in trees for u in range(n)]
    assert strict>0
    print(f'PASS {count} exact size-capped censuses; {strict} strict survivor-versus-component losses')

    for N in (2,3,5,19,64,256):
        for H in (1,100,10000):
            t = 8*N*N+H
            r = 32*N
            assert t>=4*N*N and t>=8*N*(N-1)
            # Integer upper bound log N <= N-1 in the recent budget.
            ell = 8*(t+r*(N-1)+H)
            assert Q(3*t,ell)<=Q(3,8)
            assert t-r*Q(t,16*N)<=-H
            assert t+Q(r*(N-1),2)-Q(ell,4)<=-H
    print('PASS asymmetric path split: 18 integer exponent controls')


if __name__ == '__main__':
    main()
