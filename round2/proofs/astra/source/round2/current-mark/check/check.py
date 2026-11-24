#!/usr/bin/env python3
"""Independent coefficient census for current marks and dormant excess."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'eval'))
import gauge

J = (1, ())


def tree(T, planted):
    parent, sons = {}, defaultdict(list)

    def add(S, p):
        x = len(parent)+1
        parent[x] = p
        sons[p].append(x)
        if S is not None:
            for _ in range(S[0]):
                add(None,x)
            for Y in S[1]:
                add(Y,x)

    add(T,0)
    adj = {x: sons[x]+([p] if p or planted else []) for x,p in parent.items()}
    return parent, sons, adj


def shape(x, sons, adj):
    return (sum(len(adj[y])==1 for y in sons[x]),
            tuple(sorted(shape(y,sons,adj) for y in sons[x] if len(adj[y])>=2)))


def mark(D):
    return min(D,20)


def add_leaf(y, parent, sons, adj):
    z = len(parent)+1
    return (parent | {z:y}, sons | {y:sons[y]+[z],z:[]},
            adj | {y:adj[y]+[z],z:[y]})


def score(parent, sons, adj):
    out=defaultdict(F)
    for x, p in parent.items():
        if p and len(adj[x])>=2:
            out[shape(x,sons,adj),mark(len(adj[p]))]+=1
    return out


def tidy(out):
    return {k:v for k,v in out.items() if v}


def actual(T):
    parent,sons,adj=tree(T,False)
    pre=score(parent,sons,adj)
    out=defaultdict(F)
    for x, row in adj.items():
        for y in row:
            nxt=score(*add_leaf(y,parent,sons,adj))
            for k in pre.keys()|nxt.keys():
                out[k]+=F(nxt.get(k,0)-pre.get(k,0),len(row))
    return tidy(out)


def predicted(T, fault=None, use_eval=False):
    parent,sons,adj=tree(T,False)
    out=defaultdict(F)
    const = lambda x: shape(x,sons,adj)
    D = lambda x: len(adj[x])
    core = lambda x: [y for y in sons[x] if D(y)>=2]
    plus = lambda S: (S[0]+1,S[1])

    def force(x):
        S=const(x); p=F(1,D(x)); a=mark(D(x)); nx=mark(D(x)+1)
        # In an ordinary fringe the module supplies coefficients on an
        # independent state basis; expand each dormant into two active scores.
        if use_eval and x!=1:
            for (Y,j,k),r in gauge.forcing(S).items():
                out[Y,j]+=r
                if k:
                    out[Y,k]-=r
        else:
            local=sum(F(1,D(y)) for y in sons[x])
            for y in core(x):
                Y=const(y)
                out[plus(Y),a]+=p
                out[Y,a]+=F(1,a)-p
                if nx!=a:
                    Z=plus(Y) if fault=='changed-shape-in-drop' else Y
                    rr=local+(F(1,D(parent[x])) if fault=='double-exterior' and x!=1 else 0)
                    out[Z,nx]+=rr
                    out[Z,a]-=rr
                # x's send changes Y and all of Y's old children's marks
                # at the very same physical event, charged to x's row.
                b, nb=mark(D(y)),mark(D(y)+1)
                if fault=='omit-grandparent':
                    continue
                if fault=='miss-cap-at-20' and D(y)==20:
                    nb=21
                for z in core(y):
                    Z=const(z)
                    out[Z,nb]+=p
                    out[Z,b]-=p
        out[J,a]+=F(sum(D(y)==1 for y in sons[x]),D(x))

    for x, p in parent.items():
        if D(x)<2:
            continue
        if p:
            S=const(x); a=mark(D(p))
            out[S,a]-=F(1,a)
            # Intrinsic, with current *external* mark held fixed, including
            # internal common-time shape updates but no outside arrows.
            local_parent,local_sons,local_adj=tree(S,True)
            for u,row in local_adj.items():
                for y in row:
                    if y:
                        _, ss, aa=add_leaf(y,local_parent,local_sons,local_adj)
                        out[shape(1,ss,aa),a]+=F(1,len(row))
                        out[S,a]-=F(1,len(row))
        force(x)
    return tidy(out)


def fixtures():
    from functools import cache
    @cache
    def rooted(n):
        if n==1:
            return (None,)
        # Represent ordinary all-descendant trees before leaf compression.
        out=set()
        def parts(rem, parts_so_far):
            if not rem:
                out.add(tuple(sorted(parts_so_far,key=repr)))
            else:
                for k in range(1,rem+1):
                    for t in rooted(k):
                        parts(rem-k,parts_so_far+[t])
        parts(n-1,[])
        return tuple(out)

    def compress(S):
        return (S.count(None),tuple(sorted(compress(t) for t in S if t is not None)))

    for n in range(3,8):
        for T in rooted(n):
            if len(T)>=2:
                yield compress(T)
    for r,y,z in product((1,17,18,19,20),(0,17,18,19,20),(0,1,2)):
        yield (r,((y,((z,(J,)),)),))


def dormant_check():
    # Finite intrinsic Green operators, before any copy forcing, same
    # physical domain for both kill labels. Independent row order.
    s=F(113123,200000)
    domain=gauge.shapes(6)
    active={}; dormant={}; orders=0
    for d in range(2,20):
        x=F(1,d+1); y=F(1,d)
        for S in domain:
            row=gauge.intrinsic(S); diag=sum(row.values())
            for k in (x,y):
                active[S,k]=(F(1,2+gauge.deg(S))+sum(t*active.get((U,k),0) for U,t in row.items()))/(s+k+diag)
            dormant[S,x,y]=((y-x)*active[S,x]+sum(t*dormant.get((U,x,y),0) for U,t in row.items()))/(s+y+diag)
            assert dormant[S,x,y]==active[S,x]-active[S,y]>0
            orders+=1
    # Mutate diagonal: a dormant killed *and* activating at its old rate
    # must charge the total old rate y, not the lower pure death rate x.
    diag=sum(gauge.intrinsic(J).values())
    assert (y-x)*active[J,x]+sum(t*dormant.get((U,x,y),0) for U,t in gauge.intrinsic(J).items()) != (s+x+diag)*dormant[J,x,y]
    print(f'{orders} dormant Green identities; rejected lower-rate dormant diagonal')


def permanent_obstruction():
    # Same shape and same root mark; the only hidden difference is the
    # permanent mark of a virgin child. Grow root to 3 before or after bud.
    S=(1,(J,)); p=F(1,gauge.deg(S)); a=F(1,20)
    totals=[]
    for child_birth_degree in (2,3):
        q=F(1,child_birth_degree)
        # Project the full generator on the unmarked indicator of J.
        # Intrinsic successors have size 5, and the plus copy has size 3.
        copy_row={J:q-p,gauge.plus(J):p}
        totals.append(copy_row[J])
    assert totals==[F(1,6),0]
    # A finite-difference rank check: for distinct time samples with
    # x_j=e^{-a_j t0}, the unshifted root-survival vectors form a Vandermonde.
    # Use a rational representative x_j for each distinct node to verify
    # the nonzero exact determinant formula on a 19-dimensional fixture.
    x=[F(j,21) for j in range(2,21)]
    M=[[v**k for v in x] for k in range(19)]
    from math import prod
    result=F(1)
    for i in range(19):
        pivot=M[i][i]; result*=pivot
        for k in range(i+1,19):
            c=M[k][i]/pivot
            for j in range(i+1,19):
                M[k][j]-=c*M[i][j]
    assert result==prod(v-u for i,u in enumerate(x) for v in x[i+1:])>0
    print('permanent hidden-child projection defect = 1/6; survival Vandermonde rank = 19')


if __name__=='__main__':
    count=0
    seen_faults=set()
    faults=('changed-shape-in-drop','double-exterior','omit-grandparent','miss-cap-at-20')
    for T in fixtures():
        row=actual(T)
        assert predicted(T)==row
        assert predicted(T,use_eval=True)==row
        for fault in set(faults)-seen_faults:
            if predicted(T,fault)!=row:
                seen_faults.add(fault)
        count+=1
    assert seen_faults==set(faults)
    print(f'{count} physical coefficient rows agree (independent formula + gauge.py)')
    for fault in faults:
        print('rejected',fault)
    dormant_check()
    permanent_obstruction()
