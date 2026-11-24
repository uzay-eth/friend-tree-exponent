"""Independent exact shape growth, adjacent signed laws and Hall tests."""
from collections import defaultdict, deque
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import lcm
import argparse


def rooted(t, v, p=-1):
    return tuple(sorted(rooted(t, w, v) for w in t[v] if w != p))


def decode(s):
    g = []
    def dfs(s, p):
        v = len(g)
        g.append([] if p < 0 else [p])
        if p >= 0:
            g[p].append(v)
        for t in s:
            dfs(t, v)
    dfs(s, -1)
    return g


def unrooted(t):
    return min(rooted(t, v) for v in range(len(t)))


def rates(t):
    return [sum((F(1, len(t[w])) for w in v), F(0)) for v in t]


def add(t, v):
    g = [x[:] for x in t]
    g[v].append(len(g))
    g.append([v])
    return g


@cache
def step(s):
    t = decode(s)
    out = defaultdict(F)
    for v, r in enumerate(rates(t)):
        out[unrooted(add(t, v))] += r / len(t)
    assert sum(out.values()) == 1
    return dict(out)


def advance(law):
    new = defaultdict(F)
    for s, p in law.items():
        for t, q in step(s).items():
            new[t] += p*q
    return new


@cache
def pieces(s):
    t = decode(s)
    n = len(t)
    r = rates(t)
    d = list(map(len, t))
    rho = [r[v]/(d[v]*(d[v]+1)) for v in range(n)]
    pos, neg, universal = defaultdict(F), defaultdict(F), F(0)
    for x in range(n):
        g = add(t, x)
        u = r[x]/(d[x]+1)
        universal += u
        neg[rooted(add(g, n), n+1)] += u
        pay = F(0)
        for y in t[x]:
            delta = rho[y]-rho[x]
            if delta <= 0:
                continue
            h = add(g, y)
            a, b = rooted(h, n), rooted(h, n+1)
            if d[x] == 1:
                # Explicitly cancel i via the budget at y; applicable before mixing seeds.
                assert a == rooted(add(add(t, y), n), n+1)
            pos[a] += delta
            neg[b] += delta
        for y in t[x]:
            if d[y] == 1:
                pay += max(F(0), rho[x] - rho[y])
        assert pay <= u
    return pos, neg, universal


def signed(law):
    a = defaultdict(F)
    u = F(0)
    for s, p in law.items():
        tp, tn, tu = pieces(s)
        u += p*tu
        for t, q in tp.items():
            a[t] += p*q
        for t, q in tn.items():
            a[t] -= p*q
    assert sum(a.values()) == -u
    return ({t:q for t,q in a.items() if q>0}, {t:-q for t,q in a.items() if q<0}, u)


def rawclass(law):
    p,q=defaultdict(F),defaultdict(F)
    u=F(0)
    for s,prob in law.items():
        a,b,v=pieces(s)
        u += prob*v
        for x,w in a.items():
            k=len(x[0]); p[k] += prob*w
            if k==1: # Positive x was formerly a leaf, so cancel the same K shape.
                assert x[0].count(()) == 0
                q[0] -= prob*w
        for y,w in b.items():
            q[y[0].count(())] += prob*w
    p[1]=F(0)
    gaps=[(sum(w for j,w in q.items() if j>=k)-sum(w for j,w in p.items() if j>=k), k)
          for k,v in p.items() if v>0]
    return min(gaps) if gaps else (F(0),0),u


def direct(s):
    t = decode(s)
    a = defaultdict(F)
    for x, r in enumerate(rates(t)):
        g = add(t, x)
        for y, w in enumerate(rates(g)):
            h = add(g, y)
            if y != len(t):
                a[rooted(h, len(t))] += r*w
            a[rooted(h, len(t)+1)] -= r*w
    p, q, _ = signed({s:F(1)})
    assert {s:w for s,w in a.items() if w>0} == p
    assert {s:-w for s,w in a.items() if w<0} == q


def inject(left, right, rel):
    if len(left) > len(right):
        return False
    adj = [[j for j, b in enumerate(right) if rel(a, b)] for a in left]
    used = [-1]*len(right)
    def dfs(i, seen):
        for j in adj[i]:
            if j in seen:
                continue
            seen.add(j)
            if used[j] < 0 or dfs(used[j], seen):
                used[j] = i
                return True
        return False
    return all(dfs(i, set()) for i in sorted(range(len(left)), key=lambda i:len(adj[i])))


@cache
def ge(a,b):
    return a == b or inject(b, a, ge)


def stars_only(a,b):
    return a == b or (not any(b) and a.count(()) >= len(b))


@cache
def chops(s, parity, depth=0):
    out = set()
    for i,x in enumerate(s):
        if i and x == s[i-1]:
            continue
        rest = s[:i]+s[i+1:]
        if depth and depth%2==parity:
            out.add(rest)
        for y in chops(x, parity, depth+1):
            out.add(tuple(sorted(rest+(y,))))
    return frozenset(out)


def edge(a,b,f):
    x,y=a[0],b[0]
    if len(x)>len(y):
        return False
    leaves = y.count(()) >= len(x)
    if f=='leaves':
        return leaves
    if f=='oneswap':
        return bool(chops(a,0)&chops(b,1))
    if f=='both':
        return leaves or bool(chops(a,0)&chops(b,1))
    if f=='swapleaf':
        return leaves or bool(chops(a,0)&chops(b,1)) or inject(x,y,stars_only)
    if f=='stars':
        return inject(x,y,stars_only)
    if f=='full':
        return inject(x,y,ge)
    raise ValueError(f)


def flow(p,q,nb,force=-1):
    """Integer Dinic. If force>=0, force that supply into the source side of the cut."""
    m,n=len(p),len(q)
    S,T=m+n,m+n+1
    g=[[] for _ in range(T+1)]
    inf=sum(p)+sum(q)+1
    def add(v,w,x):
        arc=[w,x,None]
        rev=[v,0,arc]
        arc[2]=rev
        g[v].append(arc)
        g[w].append(rev)
        return arc
    for i,x in enumerate(p):
        add(S,i,x if i!=force else inf)
    arcs=[]
    for i,js in enumerate(nb):
        arcs.append([(j,add(i,m+j,inf)) for j in js])
    for j,x in enumerate(q):
        add(m+j,T,x)
    total=0
    while True:
        level=[-1]*(T+1)
        level[S]=0
        todo=[S]
        for v in todo:
            for w,c,r in g[v]:
                if c and level[w]<0:
                    level[w]=level[v]+1
                    todo.append(w)
        if level[T]<0:
            return total, [i for i in range(m) if level[i]>=0], [
                (i,j,r[2][1]) for i,a in enumerate(arcs) for j,r in a if r[2][1]]
        it=[0]*(T+1)
        def send(v,lim):
            if v==T:
                return lim
            while it[v]<len(g[v]):
                arc=g[v][it[v]]
                w,c,r=arc
                if c and level[w]==level[v]+1:
                    x=send(w,min(lim,c))
                    if x:
                        arc[1]-=x
                        r[1]+=x
                        return x
                it[v]+=1
            return 0
        while (x:=send(S,inf)):
            total+=x


def hall(A,B,f,margin=False):
    a,b=list(A),list(B)
    den=lcm(*(x.denominator for x in (*A.values(),*B.values())))
    p,q=([int(x*den) for x in D.values()] for D in (A,B))
    if not p:
        return 0,None,[],[]
    if f=='leaves' and margin:
        ks=[len(x[0]) for x in a]
        ls=[y[0].count(()) for y in b]
        slack,k=min((sum(v for j,v in zip(ls,q) if j>=k)-sum(v for j,v in zip(ks,p) if j>=k), k)
                    for k in set(ks))
        return F(max(0,-slack),den), F(slack,den), [x for x,j in zip(a,ks) if j>=k], []
    neigh=[[j for j,y in enumerate(b) if edge(x,y,f)] for x in a]
    got,cut,paths=flow(p,q,neigh)
    slack=got-sum(p)
    if not slack and margin:
        base=sum(p)
        if f in ('both','stars','swapleaf','full'):
            # For a set with min k, all negative shapes with l>=k are neighbors.
            # Thus if every threshold <= k has slack >= the current best, skip this force.
            ks=[len(x[0]) for x in a]
            ls=[y[0].count(()) for y in b]
            bound={}
            m=sum(q)
            for k in sorted(set(ks)):
                m=min(m,sum(v for j,v in zip(ls,q) if j>=k)-sum(v for j,v in zip(ks,p) if j>=k))
                bound[k]=m
            visit=sorted(range(len(p)),key=lambda i:(bound[ks[i]],-ks[i]))
        else:
            ks=[0]*len(p);bound={0:0};visit=range(len(p))
        slack=den+sum(q)
        for i in visit:
            if slack<=bound[ks[i]]:
                break
            value,sub,_=flow(p,q,neigh,i)
            if value-base<slack:
                slack,cut=value-base,sub
    return F(sum(p)-got,den), F(slack,den), [a[i] for i in cut], [(a[i],b[j],F(w,den)) for i,j,w in paths]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('mode', choices=['exact','conditional','raw'])
    parser.add_argument('n', type=int)
    parser.add_argument('--margin',action='store_true')
    parser.add_argument('--details', type=int, default=0)
    args=parser.parse_args()
    law={((),):F(1)}
    for n in range(2,args.n+1):
        if args.mode=='exact':
            if n<=8:
                for seed in law:
                    direct(seed)
            a,b,u=signed(law)
            print('birth',n+2,'census',len(law),len(a),len(b),'supply/universal',float(sum(a.values())/u),flush=True)
            for f in ['leaves','oneswap','both','stars','swapleaf','full']:
                miss,slack,cut,paths=hall(a,b,f,args.margin)
                print(f,'miss',miss,'slack',slack,'slack/U',None if slack is None else float(slack/u),flush=True)
                if args.details==n+2 and (miss or f=='full'):
                    for x in cut:
                        print('cut',x,a[x],flush=True)
                    for x,y,w in paths:
                        print('path',x,y,w,flush=True)
        elif args.mode=='raw':
            (gap,k),u=rawclass(law)
            print('raw',n,n+2,len(law),k,gap,float(gap/u),flush=True)
        else:
            for k in range(2,5):
                bad=0
                for seed in law:
                    cond={seed:F(1)}
                    for _ in range(k-2):
                        cond=advance(cond)
                    a,b,u=signed(cond)
                    miss,slack,cut,_=hall(a,b,'full')
                    if miss:
                        bad+=1
                        if bad==1:
                            print('first',n,k,seed,'miss',miss,'non-normalized','cut',cut,flush=True)
                print('seed size',n,'lookback',k,'census',len(law),'bad',bad,flush=True)
        if n<args.n:
            law=advance(law)


if __name__=='__main__':
    main()
