"""Whole-tree, finite-matrix and mark-support checks for the physical trigger."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import comb
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import residual_entries as physical
from dormant_check import solve

physical.MARK_CAP=20
S=F(3,5)
A=F(1,2)


def virgin(t):
    return t[1:]==(1,())


@cache
def rates(t):
    n,d,_=physical.statistics(t)
    p=F(1,d)
    moves=defaultdict(F,physical.intrinsic_moves(t))
    for child in t[2]:
        moves[child]+=physical.mark(child[0])-p
        moves[physical.plus(child)]+=p
    b=F(t[1],d)
    assert (b if d==2 else 0)==F(virgin(t),2)
    return S+physical.mark(t[0])+n-p,dict(moves),min(d,20),b


@cache
def walk(t,k):
    den,moves,j,b=rates(t)
    out=F(virgin(t),2)
    if k:
        out+=sum(rate*walk(u,k-1) for u,rate in moves.items() if rate)
        if j!=2 and b:
            out+=b*walk(physical.newborn(j),k-1)
    return out/den


def finite(cap):
    labels=range(2,cap+1)
    states={physical.newborn(i) for i in labels}
    todo=list(states)
    while todo:
        t=todo.pop()
        for u,p in rates(t)[1].items():
            if p and physical.statistics(u)[0]<=cap and u not in states:
                states.add(u)
                todo.append(u)
    states=sorted(states)
    idx={t:i for i,t in enumerate(states)}
    m,n=len(labels),len(states)
    M=[]
    B=[]
    killed=[]
    scalar=[]
    for t in states:
        den,moves,j,b=rates(t)
        row=[F(0)]*n
        row[idx[t]]=den
        for u,p in moves.items():
            if u in idx:
                row[idx[u]]-=p
        M.append(row)
        B.append([b*int(j==i) for i in labels])
        red=row.copy()
        if j!=2 and b:
            red[idx[physical.newborn(j)]]-=b
        killed.append(red)
        scalar.append([F(virgin(t),2)])
    U=solve(M,B)
    K=[U[idx[physical.newborn(i)]] for i in labels]
    C=[[F(int(i==j))-(v if j else F(0)) for j,v in enumerate(k)] for i,k in enumerate(K)]
    first=[z[0] for z in solve(C,[[row[0]] for row in K])]
    actual=[z[0] for z in solve(killed,scalar)]
    fresh=[i for i,t in enumerate(states) if virgin(t)]
    old=[i for i,t in enumerate(states) if not virgin(t)]
    assert all(not killed[i][j] for i in fresh for j in fresh if i!=j)
    schur=[[killed[i][j]-sum(killed[i][v]*killed[v][j]/killed[v][v] for v in fresh)
            for j in old] for i in old]
    toll=[[scalar[i][0]-sum(killed[i][v]*scalar[v][0]/killed[v][v] for v in fresh)] for i in old]
    assert [x[0] for x in solve(schur,toll)]==[actual[i] for i in old]
    for i,t in enumerate(states):
        rhs=U[i][0]
        rhs+=sum(x*y for x,y in zip(U[i][1:],first[1:]))
        assert rhs==actual[i]
    for i in labels:
        assert actual[idx[physical.newborn(i)]]==first[i-2]
    return n,first[0]


def marks_check():
    ncases=0
    for root in (2,20):
        for depth in (1,2,3):
            def two(x,y):
                ch=tuple(sorted(((x,1,()),(y,2,()))))
                return walk((root,2,ch),depth)
            assert two(3,4)-two(3,5)-two(5,4)+two(5,5)==0
            assert two(3,4)!=two(5,4)
            ncases+=1
    def chain(x,y):
        return walk((2,4,((x,3,((y,1,()),)),)),2)
    defect=chain(3,3)-chain(3,5)-chain(5,3)+chain(5,5)
    assert defect
    print('incomparable mixed differences exactly zero',ncases,
          'nonzero comparable-pair difference',defect)


def multiset(types,remaining,limit):
    out={(0,0):1}
    for size,num in enumerate(types,1):
        new=defaultdict(int)
        for (n,m),val in out.items():
            for j in range(min((remaining-n)//size,limit-m)+1):
                new[n+j*size,m+j]+=val*comb(num+j-1,j)
        out=new
    return out


def counts(H,kmax):
    bare,chain,full=[],[],[]
    for k in range(1,kmax+1):
        U=C=P=0
        plain=multiset(bare,k-1,H-1)
        marked_one=multiset(full,k-1,H-1)
        eligible={(0,0):1}
        for d in range(2,H+1):
            new=defaultdict(int)
            for (n,m),w in eligible.items():
                for (extra,j),ways in marked_one.items():
                    if n+extra<k and m+j<d:
                        new[n+extra,m+j]+=w*ways
            eligible=new
            P+=sum(val for (n,m),val in eligible.items() if n==k-1)
            for m in range(d):
                U+=plain.get((k-1,m),0)
                if m:
                    C+=sum(((d-1)*chain[j-1]-bare[j-1])*plain.get((k-1-j,m-1),0)
                           for j in range(1,k))
        bare.append(U)
        chain.append(U+C)
        full.append(P)
    return list(zip(bare,chain,full))


if __name__=='__main__':
    for cap in (3,4,5):
        n,ret=finite(cap)
        print('physical size',cap,'marked states',n,'all state scalar/resolvent identities',n,
              'virgin Schur elimination checked, killed return',ret)
    marks_check()
    for k,triple in enumerate(counts(6,6),1):
        print('D<=6 cores',k,'unmarked, chain basis, reachable nineteen-mark counts',*triple)
