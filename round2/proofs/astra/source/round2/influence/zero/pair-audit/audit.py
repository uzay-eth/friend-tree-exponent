#!/usr/bin/env python3
"""Rebuild synchronized source clocks on concrete labelled planted trees."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product

SILENT=-1
ROOT=0


def neighbors(T):
    A=defaultdict(set)
    for x,y in T:
        A[x].add(y)
        if y!=SILENT:
            A[y].add(x)
    return A


def arity(T):
    return neighbors(T)


def star(kids):
    """A and B edges; i-th child receives literal lists of pendant labels."""
    outs=[{(ROOT,SILENT)} for _ in range(2)]
    for i,pair in enumerate(kids,1):
        for T,ls in zip(outs,pair):
            if ls is None:
                continue
            T.add((ROOT,i))
            for j in ls:
                T.add((i,j))
    return tuple(map(frozenset,outs))


def source_rows(T,U,flaw=None):
    A,B=map(arity,(T,U))
    for x in sorted(set(A)|set(B)):
        Na,Nb=A.get(x,()),B.get(x,())
        if not Na or not Nb:
            for y in Na: yield x,y,None,Q(1,len(Na))
            for y in Nb: yield x,None,y,Q(1,len(Nb))
            continue
        pa,pb=Q(1,len(Na)),Q(1,len(Nb))
        a,b={y:pa for y in Na},{y:pb for y in Nb}
        total=0
        for y in sorted(set(a)&set(b)):
            m=min(a[y],b[y]); a[y]-=m; b[y]-=m; total+=m
            yield x,y,y,m
        if total!=1:
            for y,z in product(a,b):
                w=a[y]*b[z]/(1-total)
                if w:
                    if flaw=='split-clock' and x==ROOT and y!=SILENT and z!=SILENT:
                        yield x,y,None,w
                        yield x,None,z,w
                    else:
                        yield x,y,z,w


def grow(T,U,y,z):
    j=1+max(set(arity(T))|set(arity(U)))
    if y==z and y not in (None,SILENT):
        return T|{(y,j)},U|{(z,j)}
    return (T if y in (None,SILENT) else T|{(y,j)},
            U if z in (None,SILENT) else U|{(z,j+1)})


def sub(T,x):
    # All edges are recorded in ancestry order.
    got={x}; R=set()
    for y,z in sorted(T):
        if y in got and z!=SILENT:
            got.add(z); R.add((y,z))
    return frozenset(R)


def canon(T,x):
    return tuple(sorted(canon(T,z) for y,z in T if y==x and z!=SILENT))


def six(T,U):
    A,B=map(arity,(T,U))
    depth={ROOT:0}
    for x,y in sorted(T&U):
        if y!=SILENT: depth[y]=depth[x]+1
    for x,r in depth.items():
        innocent,dominant=(A[x],B[x]) if r%2 else (B[x],A[x])
        assert innocent<=dominant, (x,r,innocent,dominant)
    da,db=A[0]-{SILENT},B[0]-{SILENT}
    assert db<=da and len(da)>len(db)
    out=[0]*6
    for x in da:
        if x not in db:
            g=4 if len(A[x])==1 else 5
        elif len(A[x])==len(B[x])==1:
            g=0
        elif len(A[x])==1:
            g=1
        else:
            assert len(B[x])>1 and len(A[x])<=len(B[x])
            g=2 if sub(T,x)==sub(U,x) else 3
        out[g]+=1
    return tuple(out)  # u,v,c,e,l,k


def jump(a,inc):
    return tuple(t+d for t,d in zip(a,inc))


def assumed(T,U,flaw=None):
    # Direct decomposition of concrete outgoing clocks, using local rates.
    a=six(T,U); u,v,c,e,l,k=a
    A,B=map(arity,(T,U)); pa,pb=Q(1,len(A[0])),Q(1,len(B[0]))
    out=defaultdict(Q)
    def pay(w,*d):
        b=jump(a,d)
        if w and b!=a: out[b]+=w
    pay(u,1,0,0,0,0,0)
    pay(l,0,0,0,0,1,0)
    for x in A[0]-{SILENT}:
        if x not in B[0]:
            if len(A[x])>1 and flaw!='drop-private':
                pay(Q(1,len(A[x])),0,0,0,0,1,0)
        elif len(A[x])==1<len(B[x]):
            q=Q(1,len(B[x]))
            pay(q,1,0,0,0,0,0); pay(1-q,0,0,0,0,1,0)
        elif len(A[x])>1:
            s,t=Q(1,len(A[x])),Q(1,len(B[x]))
            pay(t,1,0,0,0,0,0); pay(s-t,0,0,0,0,1,0)
    pay(pa*u,-1,0,1,0,0,0)
    pay(pa*v,0,-1,0,1,0,0)
    for number,dk in ((l,(0,0,0,0,-1,1)),(k,(0,)*6)):
        for count,dc in [(u,(-1,1,0,0,0,0)),(c,(0,0,-1,1,0,0)),(v+e+1,(0,)*6)]:
            pay(pa*pb*number*count,*(x+y for x,y in zip(dk,dc)))
    return dict(out)


def compare(T,U,flaw=None):
    a=six(T,U); A,B=map(arity,(T,U))
    direct=defaultdict(Q); ma=defaultdict(Q); mb=defaultdict(Q)
    for x,y,z,w in source_rows(T,U,flaw):
        if y is not None: ma[x,y]+=w
        if z is not None: mb[x,z]+=w
        a1=six(*grow(T,U,y,z))
        if a1!=a: direct[a1]+=w
    assert ma=={(x,y):Q(1,len(ns)) for x,ns in A.items() for y in ns}
    assert mb=={(x,y):Q(1,len(ns)) for x,ns in B.items() for y in ns}
    assert direct==assumed(T,U,flaw),(a,direct,assumed(T,U,flaw))
    u,v,c,e,l,k=a; p,q=Q(1,len(A[0])),Q(1,len(B[0]))
    mu=defaultdict(Q)
    for tree,ns,fac in ((U,B,q),(T,A,-p)):
        for x in ns[0]-{SILENT}:
            if len(ns[x])>1: mu[canon(tree,x)]+=fac
    exact=sum(map(abs,mu.values()))
    cap=c*(q-p)+e*(p+q)+v*q+k*p
    assert exact<=cap
    return exact,cap,len(direct)


def fixtures():
    # One clean and one dirty common core; LC, LL, private L and private core.
    H=star([([],[]),([],[41]),([42],[42]),([43],[43,44]),([],None),([45],None)])
    X=star([([41],[41,42]),([],None)]) # one dirty, low degrees
    Y=star([([],[]),([],None)])
    return H,X,Y


def heldout():
    """Two uniformization steps including every descendant source clock."""
    def reward(T,U):
        out=defaultdict(Q)
        for tree,sign in ((T,-1),(U,1)):
            N=arity(tree)
            for x in N[ROOT]-{SILENT}:
                if len(N[x])>1: out[canon(tree,x)]+=Q(sign,len(N[ROOT]))
        return out

    # Start from a genuine unmatched root send into a shared deep branch.
    T,U=star([([41],[41]),([],None)])
    T,U=T|{(41,42)},U|{(41,42)}
    T,U=grow(T,U,2,1)
    a=six(T,U)
    assert a==(0,0,0,1,0,1)
    lam=Q(len(arity(T))+len(arity(U))+4)

    def push(batch):
        out=defaultdict(Q)
        for (P,H),mass in batch.items():
            compare(P,H)
            total=0
            for x,y,z,w in source_rows(P,H):
                out[grow(P,H,y,z)]+=mass*w/lam
                total+=w
            assert total<=lam
            out[P,H]+=mass*(1-total/lam)
        assert sum(out.values())==1
        return out

    batch={(T,U):Q(1)}
    for age in range(3):
        mu=defaultdict(Q); envelope=0
        for (P,H),w in batch.items():
            exact,cap,_=compare(P,H)
            envelope+=w*cap
            for x,z in reward(P,H).items(): mu[x]+=w*z
        tv=sum(map(abs,mu.values()))
        assert tv<=envelope
        print('held-out jump age',age,'states',len(batch),'signed TV',tv,'six envelope',envelope)
        if age<2: batch=push(batch)


if __name__=='__main__':
    for T,U in fixtures(): print(six(T,U),'true TV, allowed TV, active successors:',*compare(T,U))
    for defect in ('split-clock','drop-private'):
        try: compare(*fixtures()[0],defect)
        except AssertionError: print('rejected',defect)
        else: raise AssertionError('accepted '+defect)
    heldout()
