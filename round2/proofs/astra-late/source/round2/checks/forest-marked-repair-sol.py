#!/usr/bin/env python3
"""Independent controls for the marked-forest truncation and sign alphabet."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product


def binary_control():
    types = {((),)*j for j in range(1,5)}
    for h in range(4):
        leaves = 2**h
        assert len(types) >= 2**(leaves+1)
        if h < 3:
            nxt = {tuple(sorted((a,b))) for a,b in product(types,repeat=2)}
            assert len(nxt) == len(types)*(len(types)+1)//2
            types = nxt
    return len(types)


def graph(par):
    g = [set() for _ in par]
    for w,v in enumerate(par[1:],1):
        g[w].add(v)
        g[v].add(w)
    return g


def direct(par,N,B,J,j,K):
    g=graph(par)
    low={i for i,a in enumerate(g) if len(a)<=N}
    cores={i for i,a in enumerate(g) if len(a)>1}
    def shape(v,u,V,point):
        return (v in cores,v==point,tuple(sorted(shape(w,v,V,point) for w in g[v]&V if w!=u)))
    left=set(low)
    count,seen=Counter(),Counter()
    missed=[0,0]
    while left:
        top=min(left)
        q=[top]; d={top:0}; left.remove(top)
        for v in q:
            for w in g[v]&left:
                left.remove(w); q.append(w); d[w]=d[v]+1
        V=set(q); C=V&cores
        if not C:continue
        p=par[top]
        if not top or len(g[p])>B or len(g[p]&cores)>K or len(V)>J or max(d[v] for v in C)>=j:
            missed[0]+=1; missed[1]+=len(C); continue
        assert len(g[p])>N and top in C
        mark=Q(len(g[p]),N)
        count[mark,shape(top,-1,V,-1)]+=1
        for point in C:seen[mark,shape(top,-1,V,point)]+=1
    return count,seen,tuple(missed)


def at_parent(par,N,B,J,j,K):
    g=graph(par)
    kids=defaultdict(list)
    for v,u in enumerate(par[1:],1):kids[u].append(v)
    degree=list(map(len,g)); core=[d>=2 for d in degree]
    count,seen=Counter(),Counter()
    raw_count=0
    for v,d in enumerate(degree):
        if d<=N:continue
        tops=[w for w in kids[v] if core[w] and degree[w]<=N]
        raw_count+=len(tops)
        if d>B or sum(core[w] for w in g[v])>K:continue
        for w in tops:
            total=0; accepted=[]
            def go(x,h):
                nonlocal total
                total+=1
                if total>J or (core[x] and h>=j):raise OverflowError
                if core[x]:accepted.append(x)
                return x,tuple(go(y,h+1) for y in kids[x] if degree[y]<=N)
            try:t=go(w,0)
            except OverflowError:continue
            def rooted(t,point):
                x,a=t
                return core[x],x==point,tuple(sorted(rooted(y,point) for y in a))
            mark=Q(d,N)
            count[mark,rooted(t,-1)]+=1
            for point in accepted:seen[mark,rooted(t,point)]+=1
    root_component=int(core[0] and degree[0]<=N)
    return count,seen,(raw_count+root_component-sum(count.values()),
                       sum(1 for d in degree if 2<=d<=N)-sum(seen.values()))


def path_split_control():
    """Enumerate every path in small old/new trees with an arbitrary old prefix."""
    cases = [(-1, 0, 0)]
    counts = Counter()
    for n in range(3, 8):
        for par in cases:
            g = graph(par)
            for m in range(1, n):
                anchor = list(range(m))
                depth = [0] * m
                for u in par[m:]:
                    anchor.append(anchor[u])
                    depth.append(depth[u] + 1)
                for start in range(n):
                    paths = [(start,)]
                    for p in paths:
                        paths.extend(p + (v,) for v in g[p[-1]] if len(p) == 1 or v != p[-2])
                        if start == p[-1]:
                            continue
                        old = [i for i, v in enumerate(p) if v < m]
                        if old:
                            a, b = old[0], old[-1]
                            assert old == list(range(a, b+1)) and p[a] == anchor[start]
                            assert depth[start] == a and depth[p[-1]] == len(p)-1-b
                            counts['old'] += 1
                        else:
                            assert all(anchor[v] == anchor[start] for v in p)
                            i = min(range(len(p)), key=lambda j: depth[p[j]])
                            assert i <= depth[start] and len(p)-1-i <= depth[p[-1]]
                            counts['recent'] += 1
        cases = [p+(v,) for p in cases for v in range(n)]
    return counts['old'], counts['recent']


def main():
    print('binary decorated skeleton types at height 3:',binary_control())
    cases=[(-1,0,0)]
    matches=omissions=changes=0
    for n in range(3,8):
        for p in cases:
            for N,(B,J,j,K) in product((2,3,4),((4,3,1,2),(8,4,3,3),(16,20,20,16))):
                x=direct(p,N,B,J,j,K)
                assert x==at_parent(p,N,B,J,j,K),(p,N,B,J,j,K)
                matches+=1; omissions+=bool(sum(x[2]))
                if n<7:
                    for r in range(n):
                        y=at_parent(p+(r,),N,B,J,j,K)
                        for a,b in zip(x[:2],y[:2]):
                            difference={z:b[z]-a[z] for z in a.keys()|b.keys()}
                            norm=sum(map(abs,difference.values()))
                            assert norm==sum(v*(1 if v>=0 else -1) for v in difference.values())
                            assert norm<=6*(min(J,j)+5)*min(B,K)*J
                            changes+=norm>0
        cases=[p+(v,) for p in cases for v in range(n)]
    print(f'{matches} exact component/pointed-survivor censuses, {omissions} nonempty remainders, {changes} nonzero signed changes')
    old, recent = path_split_control()
    assert old and recent
    print(f'{old} paths through the old tree, {recent} wholly recent paths')


if __name__=='__main__':main()
