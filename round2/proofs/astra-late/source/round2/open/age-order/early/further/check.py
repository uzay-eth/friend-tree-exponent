"""Exact fixed-point certificate with unlabelled tagged trees."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import lcm
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

spec=spec_from_file_location('early_bounds',Path(__file__).resolve().parents[1]/'check.py')
base=module_from_spec(spec)
spec.loader.exec_module(base)
silent_bound, core_bound, out=base.silent_bound,base.core_bound,base.out

UNIT = 10**12


def down(x):
    return (UNIT*x)//1


def up(x):
    return -down(-x)


def encode(adj, a, b=None):
    """Root at first tag, with the second tag marked regardless of ancestry."""
    def visit(u,p):
        return (int(u==b), tuple(sorted(visit(v,u) for v in adj[u] if v!=p)))
    return visit(a,-1)


def decode(t):
    adj=[]; b=None
    def visit(x,p):
        nonlocal b
        i=len(adj); adj.append([])
        if p>=0:
            adj[i].append(p); adj[p].append(i)
        if x[0]: b=i
        for y in x[1]: visit(y,i)
    visit(t,-1)
    return adj,b


def grow(adj,u):
    nxt=[x.copy() for x in adj]
    nxt[u].append(len(adj)); nxt.append([u])
    return nxt


@cache
def absolute(d,l,ns):
    if d==1:
        D,L=ns[0]
        return (down(core_bound(2,1,12)/(1+D*D)), up(1-silent_bound(L-1,max(32,L-1),1,D-L+1)))
    lo=core_bound(d,l,max(d,12))
    hi=Q(1) if l else min(Q(1),sum(1-silent_bound(L,max(32,L),0,D-L) for D,L in ns))
    return down(lo),up(hi)


def endpoint(adj,i):
    def dl(j): return len(adj[j]),sum(len(adj[k])==1 for k in adj[j])
    d,l=dl(i)
    return absolute(d,l,tuple(sorted(dl(j) for j in adj[i])))


def rates(adj):
    return tuple(sum(Q(1,len(adj[v])) for v in ns) for ns in adj)


@cache
def bound(t,h):
    adj,b=decode(t)
    swapped=encode(adj,b,0)
    if t==swapped: return 0,0
    if t>swapped:
        l,u=bound(swapped,h)
        return -u,-l
    al,au=endpoint(adj,0); bl,bu=endpoint(adj,b)
    lo,hi=al-bu,au-bl
    if h:
        den=lcm(*(len(ns) for ns in adj))
        l=u=0
        for i,ns in enumerate(adj):
            w=sum(den//len(adj[v]) for v in ns)
            x,y=bound(encode(grow(adj,i),0,b),h-1)
            l+=w*x;u+=w*y
        den*=len(adj)
        lo=max(lo,l//den);hi=min(hi,-((-u)//den))
    assert lo<=hi
    return lo,hi


def placements(adj,i,j,h):
    n=len(adj)
    return tuple(Q(x,UNIT) for x in bound(encode(grow(grow(adj,i),j),n,n+1),h))


def scale(t,w):
    return (w*t[0],w*t[1]) if w>=0 else (w*t[1],w*t[0])


def add(terms):
    lo=hi=Q(0)
    for x,y in terms: lo+=x;hi+=y
    return lo,hi


def compare(adj,h):
    n=len(adj); r=rates(adj)
    w=[r[i]/(len(ns)+1) for i,ns in enumerate(adj)]
    delta=[w[i]/len(ns) for i,ns in enumerate(adj)]
    it=[placements(adj,i,n,h) for i in range(n)]
    ed={(i,j):placements(adj,i,j,h) for i,ns in enumerate(adj) for j in ns if i<j}
    whole=add([scale(t,delta[j]-delta[i]) for (i,j),t in ed.items()]+[scale(t,x) for t,x in zip(it,w)])
    rows=[]
    for i,ns in enumerate(adj):
        rows.append(add([scale(it[i],w[i])]+[scale(ed[min(i,j),max(i,j)],-delta[i] if i<j else delta[i]) for j in ns]))
    owner=[i if len(ns)>1 else ns[0] for i,ns in enumerate(adj)] if n>2 else [0,0]
    grouped=defaultdict(list)
    for i,t in enumerate(it): grouped[owner[i]].append(scale(t,w[i]))
    for (i,j),t in ed.items():
        if owner[i]==owner[j]:
            grouped[owner[i]].append(scale(t,delta[j]-delta[i]))
        else:
            grouped[owner[i]].append(scale(t,-delta[i]))
            grouped[owner[j]].append(scale(t,delta[j]))
    return scale(whole,Q(1,n*(n+1))),rows,{u:add(ts) for u,ts in grouped.items()}


def classes(n):
    """Generate all free tree shapes by one-vertex extensions."""
    lev={((0, ((0, ()),))):Q(1)}
    for _ in range(2,n):
        nxt=defaultdict(Q)
        for t,m in lev.items():
            adj,_=decode(t)
            for i,w in enumerate(rates(adj)):
                y=grow(adj,i)
                nxt[min(encode(y,j) for j in range(len(y)))]+=m*w/len(adj)
        lev=nxt
    assert sum(lev.values())==1
    return lev


def fmt(t): return '['+out(t[0])+', '+out(t[1],True)+']'


def early(n,h):
    tot=[]
    for t,m in classes(n).items():
        a,_=decode(t)
        tot.append(scale(compare(a,h)[0],m))
    return add(tot)


def control():
    # Cross-check rate totals, all shape weights, and signed regrouping on
    # independently enumerated increasing-parent trees.
    from itertools import product
    for n in range(2,7):
        agg=defaultdict(Q)
        for par in product(*(range(i) for i in range(1,n))):
            adj=[[] for _ in range(n)]
            for i,u in enumerate(par,1): adj[i].append(u);adj[u].append(i)
            a=[[]]; mass=Q(1)
            for u in par:
                if len(a)>1: mass*=rates(a)[u]/len(a)
                a=grow(a,u)
            agg[min(encode(adj,i) for i in range(n))]+=mass
        assert agg==classes(n)
    # A first send followed by the finite core lower bound cannot exceed
    # the direct leaf upper bound over these tested degree/stock pairs.
    for D in range(2,65):
        for L in range(1,D+1):
            lo,hi=absolute(1,0,((D,L),))
            assert 0<lo<=hi<=UNIT
    print('shape laws through size 6 match labelled enumeration; leaf bounds ordered')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--horizon',type=int,default=5)
    p.add_argument('--n',nargs='+',type=int,default=[4,5,6]); p.add_argument('--scan',action='store_true');p.add_argument('--control',action='store_true')
    p.add_argument('--fixture',action='store_true')
    args=p.parse_args()
    if args.control: control()
    if args.fixture:
        a=[[1,2],[0],[0,3,5],[2,4],[3],[2,6],[5]]
        x=compare(a,6)[0]; assert x[0]>0
        print('size-seven tripod with arms of length two:',fmt(x))
        bound.cache_clear()
    for n in args.n:
        if args.scan:
            low=[]; vertex=[]; group=[]; counts=[0,0,0]
            for t in classes(n):
                a,_=decode(t); gap,rows,groups=compare(a,args.horizon)
                counts[0]+=gap[0]>0
                counts[1]+=min(x[0] for x in rows)>0
                counts[2]+=min(x[0] for x in groups.values())>0
                low.append((gap,t));vertex.extend((x,t,u) for u,x in enumerate(rows))
                group.extend((x,t,u) for u,x in groups.items())
            print('scan',n,len(low),'proved positive whole/all-rows/all-groups',counts,
                  'negative whole/vertex/group',[sum(x[0][1]<0 for x in z) for z in (low,vertex,group)])
            for k,z in zip(('whole','vertex','group'),(low,vertex,group)):
                x,t,*u=min(z,key=lambda x:x[0][1]); print(k,fmt(x),'adj',decode(t)[0],u)
            for gap,t in low:
                if gap[0]<=0: print('unresolved whole',fmt(gap),decode(t)[0])
        else:
            v=early(n,args.horizon)
            if args.horizon==5 and 4<=n<=9: assert v[0]>0
            print(f'p{n+1}-p{n+2}, n={n}, h={args.horizon}: {fmt(v)}')
        print('cache',bound.cache_info(),flush=True)
        bound.cache_clear()
