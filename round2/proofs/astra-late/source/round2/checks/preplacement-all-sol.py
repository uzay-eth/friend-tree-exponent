"""Independent physical contraction and exact adverse-host controls."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parents[1] / 'open/age-order/global-transport/hostile'
D = run_path(str(ROOT / 'discharge/check.py'))
G, U, L = [D[k] for k in ('G', 'U', 'L')]
A, B, F = Q('.00634'), Q('.017'), Q('.00356')


def physical(links, sizes):
    tree = [[] for _ in sizes]
    for i, j in links:
        tree[i].append(j)
        tree[j].append(i)
    for i, d in enumerate(sizes):
        assert len(tree[i]) <= d
        for _ in range(d-len(tree[i])):
            tree[i].append(len(tree))
            tree.append([i])
    return tuple(map(tuple, tree))


def append(t, u):
    return tuple(x+(len(t),) if k == u else x for k, x in enumerate(t))+((u,),)


def check_tree(t):
    d = list(map(len, t))
    core = {u for u,k in enumerate(d) if k>1}
    nbr = {u: set(t[u]) & core for u in core}
    w = [sum(Q(1,d[j]) for j in row)/(d[u]+1) for u,row in enumerate(t)]
    s = [w[u]/d[u] for u in range(len(t))]

    def debit(u,v):
        if s[v]>s[u]:
            u,v = v,u
        final=append(append(t,u),v)
        j=sum(len(final[a])==1 for a in final[u])-1
        return (s[u]-s[v])*max(0,U(d[u]+1,j)-L(d[v]+1))

    credit = {u:w[u]*G(d[u]+1)+sum(w[j]*G(2)-debit(u,j) for j in t[u] if j not in core)
              for u in core}
    keep = {u for u in core if d[u]!=2 or len(nbr[u])!=2}
    paths = []
    for u in keep:
        for v in nbr[u]:
            route = [u,v]
            while route[-1] not in keep:
                (j,) = nbr[route[-1]]-{route[-2]}
                route.append(j)
            if u < route[-1]:
                paths.append(route)
    assert len(paths)==len(keep)-1
    assert sorted(x for route in paths for x in route[1:-1])==sorted(core-keep)
    assert sum(len(nbr[u])-2 for u in keep)==-2
    reduced = {u: set() for u in keep}
    paid={u: Q(0) for u in keep}
    moved=0
    simultaneous=set()
    violations=[]
    chain_count=0

    def bill(u):
        return Q('.021') if d[u]==2 else Q('.0037') if d[u]<10 else Q('.00625')

    def excess(u,v):
        e=d[v]
        return max(0, max(0,s[u]-Q(1,d[u]*e*(e+1))) *
                   max(0,U(d[u]+1,d[u]-len(nbr[u]))-L(e+1))-B)

    cover=Q(0)
    for route in paths:
        u,v = route[0],route[-1]
        reduced[u].add(v)
        reduced[v].add(u)
        net=sum(debit(i,j) for i,j in zip(route,route[1:]))-sum(credit[j] for j in route[1:-1])
        budget=B
        if len(route)==2:
            for x,y in ((u,v),(v,u)):
                charge=excess(x,y)
                budget+=charge
                paid[x]+=charge
                if d[x]>=10 and d[y] in (3,4):
                    paid[x]+=F
                    paid[y]-=F
                    simultaneous.add((x,'direct'))
                    moved+=1
        else:
            chain_count+=1
            for x in (u,v):
                budget+=bill(x)
                paid[x]+=bill(x)
                simultaneous.add((x,'chain'))
            if net > budget-min(bill(u),bill(v)):
                violations.append(('omit-one-end',tuple(route)))
        assert net<=budget
        cover+=net
    assert all(len(reduced[u])==len(nbr[u]) for u in keep)
    seen=set()
    todo=[next(iter(keep))]
    while todo:
        u=todo.pop()
        if u not in seen:
            seen.add(u)
            todo.extend(reduced[u]-seen)
    assert seen==keep
    for u in keep:
        assert credit[u]+A*(len(nbr[u])-2)-B-paid[u]>=0
    assert sum(credit[u] for u in keep)-cover>=B+2*A
    return len(paths),chain_count,moved,sum((x,'direct') in simultaneous for x,k in simultaneous if k=='chain'),len(violations)


def host_tail():
    g=Q('.16694235')
    assert D['H']['rescue'](2,1)-D['H']['upper'](2,0)>g
    assert U(2,0)<Q('.611') and G(2)>Q(2,113)
    m=g*Q(113,114)+G(2)/2-Q('.611')/2
    assert m<0 and m/2+(A-F)*Q(109,2)-2*F-B>0
    assert g*Q(113,228)-2*F-B>0
    assert g*Q(113,228)-2*Q('.00625')-B>0
    control=0
    for d in range(10,42):
        for l in range(d):
            assert Q(d+l,2*d*(d+1))*U(d+1,l)<Q('.012')
            control+=1
    assert Q(1,86)<Q('.012')
    return control


if __name__=='__main__':
    early=D['H']['EARLY']
    assert early['adjacent']((0,),4)[0]>0
    assert (early['silent_bound'](1,4)-Q(1,4))/2==Q(127,2440)>0
    assert all(D['H']['score'](D['H']['spider'](d,1))[0]>0 for d in (2,3,4,9,36))
    fixtures=[physical([(0,1)],[2,2]),physical([(0,1)],[2,113])]
    for k in (1,2,3,4,12):
        fixtures.append(physical([(j,j+1) for j in range(k+1)],[2]*(k+2)))
        fixtures.append(physical([(j,j+1) for j in range(k+1)],[10]+[2]*k+[42]))
    # High donor with competing direct and chain ports; low recipient with both types.
    fixtures.append(physical([(0,1),(0,2),(2,3),(1,4),(4,5),(0,6),(3,7),(3,8),(5,9)],
                             [10,3,2,113,2,4,3,2,10,2]))
    rows=[check_tree(t) for t in fixtures]
    # Sharp high/high chain: both ends' bills are essential to this bound.
    # Use the permissible infinite-degree endpoint of the high-tail interval.
    u=U(3,0); k=G(3); ell=L(3)
    x,y=Q(1,10),Q(0)
    s,t=(Q(1,2)+x)/6,(Q(1,2)+y)/6
    raw=max(Q('.012'),s*(u-L(11)))+max(Q('.012'),t*u)+(u-ell)*abs(s-t)-2*k*(s+t)
    assert B+Q('.00625')<raw<B+2*Q('.00625')
    assert sum(r[3] for r in rows)>0
    print(f'PASS {len(rows)} physical mixed/path fixtures, {sum(r[0] for r in rows)} reduced edges, {sum(r[1] for r in rows)} suppressed chains, {sum(r[2] for r in rows)} direct transfers')
    print('PASS seed, five stars, simultaneous payments and planted missing-chain-end error')
    print(f'PASS {host_tail()} actual-stock high chain exposures and analytic high-degree corners')
