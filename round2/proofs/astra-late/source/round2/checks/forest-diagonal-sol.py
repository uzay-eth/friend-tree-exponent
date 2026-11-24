"""Independent exact checks of the diagonal forest transfer."""

from collections import Counter
from fractions import Fraction as Q
from itertools import product


def graph(par):
    g = [set() for _ in range(len(par)+1)]
    for i,p in enumerate(par,1):
        g[i].add(p)
        g[p].add(i)
    return g


def shape(g, low, v, p=-1, distinguished=-1):
    return v == distinguished, tuple(sorted(shape(g,low,u,v,distinguished) for u in g[v]&low if u!=p))


def compare(par,N,j,B):
    g = graph(par)
    d = list(map(len,g))
    low = {i for i,k in enumerate(d) if k<=N}
    core = {i for i,k in enumerate(d) if k>=2}
    remaining = set(low)
    direct, seen = Counter(), Counter()
    dropped = founder = orientation = 0
    while remaining:
        r = min(remaining)
        distance = {r:0}
        queue = [r]
        remaining.remove(r)
        for v in queue:
            for w in g[v]&remaining:
                distance[w] = distance[v]+1
                queue.append(w)
                remaining.remove(w)
        cc = set(queue)&core
        if not cc:
            continue
        if r == 0 or d[par[r-1]]>B or max(distance[v] for v in cc)>=j:
            dropped += len(cc)
            founder += len(cc) if r==0 else 0
            continue
        ratio = Q(d[par[r-1]],N)
        direct[ratio,shape(g,low,r),len(cc)] += 1
        for v in cc:
            seen[ratio,shape(g,low,r,distinguished=v),len(cc)] += 1
    counts, survivors = Counter(), Counter()
    for p in core-low:
        if d[p]>B:
            continue
        for r in g[p]&core&low:
            if r==0 or par[r-1]!=p:
                orientation += 1
                continue
            todo = [(r,0)]
            members = []
            rejected = False
            while todo:
                v,k = todo.pop()
                if k==j:
                    rejected = True
                    break
                members.append(v)
                todo.extend((u,k+1) for u in g[v]&core&low if u and par[u-1]==v)
            if rejected:
                continue
            ratio = Q(d[p],N)
            counts[ratio,shape(g,low,r),len(members)] += 1
            for v in members:
                survivors[ratio,shape(g,low,r,distinguished=v),len(members)] += 1
    assert counts==direct and survivors==seen
    assert sum(survivors.values())+dropped==len(low&core)
    return counts,survivors,founder,orientation


def norm(a,b):
    hi = sum(max(0,a[x]-b[x]) for x in a.keys()|b.keys())
    lo = sum(max(0,b[x]-a[x]) for x in a.keys()|b.keys())
    assert hi+lo == sum((1 if a[x]>=b[x] else -1)*(a[x]-b[x]) for x in a.keys()|b.keys())
    return hi+lo


def hidden_kind_fixture():
    # Same parent degree 3 and observed rooted edge at N=2; the distal
    # vertex is originally a leaf in A, a core in B with a high child.
    a = (0,0,0,3)
    b = (0,0,0,3,4,5,5)
    x, *_ = compare(a,2,3,8)
    y, *_ = compare(b,2,3,8)
    (xkey,), (ykey,) = tuple(x),tuple(y)
    assert xkey[:2]==ykey[:2] and (xkey[2],ykey[2])==(1,2)
    # Mass norm (not the probability-convention half) sees both atoms.
    assert norm(x,y)==2


def arithmetic():
    assert Q(99,100)+2*Q(11,20)-2 == Q(9,100)
    assert Q(99,100)+4*Q(1,100) == Q(103,100)
    assert 2*Q(11,20)-Q(103,100)==Q(7,100)
    assert Q(11,20)+2*Q(1,100)<Q(103,100)
    # Worst radius at N=Nmax; write the log term as a nonnegative
    # rational stand-in for log log n, and retain exact roundings.
    def ceil(x):
        return -(-x//1)
    nchecks=0
    for N,L in product(range(2,65),(Q(0),Q(1,13),Q(1),Q(79,10),Q(25))):
        j=3*ceil((132*N**3+100*L)/3)
        j0=3*(j//6)
        assert j0>=64*N**3
        assert j0>=66*N**3+50*L-3
        assert Q(2,25)*j0>=4*L
        nchecks+=1
    return nchecks


if __name__=='__main__':
    hidden_kind_fixture()
    states=changes=founder=wrong=0
    for n in range(3,8):
        for rest in product(*(range(i) for i in range(3,n))):
            p=(0,0)+rest
            for N,j,B in product((2,3,4),(1,2,3),(3,5,8)):
                a,b,f,o=compare(p,N,j,B)
                founder+=f>1
                wrong+=o
                states+=1
                if n<7:
                    for i in range(n):
                        x,y,*_=compare(p+(i,),N,j,B)
                        na,nb=norm(a,x),norm(b,y)
                        assert na<=6*(j+4)*B and nb<=6*(j+4)*B*(N+1)**(j+3)
                        changes+=bool(na or nb)
    assert founder and wrong and changes
    print(f'{states} exact censuses; {changes} nonzero one-birth changes; {founder} multi-core founder residuals; {wrong} ancestor exclusions')
    print(f'{arithmetic()} radius controls; clock exponents 9/100, 103/100, 7/100; same observed edge with original K=1 and K=2')
