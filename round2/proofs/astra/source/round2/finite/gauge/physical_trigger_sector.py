"""Exact-integer powers for the killed virgin return using ancestral mark paths."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import ceil, floor, lcm
from time import monotonic

from physical_trigger_check import counts, physical, virgin


@cache
def geometry(t):
    d=1+t[1]+len(t[2])
    parts=[geometry(y) for y in t[2]]
    return max([d]+[x for x,k in parts]),1+sum(k for x,k in parts)


@cache
def late(t):
    tag,L,ch=t
    j=min(20,1+L+len(ch))
    return tag,L,tuple(sorted((j,*late(y)[1:]) for y in ch))


@cache
def retain_path(t,path):
    if not path:
        return late(t)
    tag,L,ch=t
    j=min(20,1+L+len(ch))
    return tag,L,tuple(sorted(retain_path(y,path[1:]) if i==path[0] else (j,*late(y)[1:])
                              for i,y in enumerate(ch)))


def vertices(t,prefix=()):
    j=min(20,1+t[1]+len(t[2]))
    for idx,y in enumerate(t[2]):
        path=prefix+(idx,)
        if y[0]!=j:
            yield path
        yield from vertices(y,path)


@cache
def expansion(t):
    weights=defaultdict(int,{late(t):1})
    for path in vertices(t):
        weights[retain_path(t,path)]+=1
        weights[retain_path(t,path[:-1])]-=1
    return tuple((u,v) for u,v in weights.items() if v)


def row(t,H,K,rate,unit,basis):
    n,d,_=physical.statistics(t)
    root=physical.mark(t[0])
    p=F(1,d)
    moves=defaultdict(F,physical.intrinsic_moves(t))
    moves[t]+=rate-root-n+p
    assert moves[t]>=0
    if d!=2 and t[1]:
        moves[physical.newborn(min(d,20))]+=F(t[1],d)
    for y in t[2]:
        moves[y]+=physical.mark(y[0])-p
        moves[physical.plus(y)]+=p
    combined=defaultdict(int)
    for y,weight in moves.items():
        D,k=geometry(y)
        if D>H or k>K or not weight:
            continue
        scaled=weight*unit
        assert scaled.denominator==1 and scaled>=0
        for state,coeff in (expansion(y) if basis else ((y,1),)):
            combined[state]+=scaled.numerator*coeff
    return tuple((u,v) for u,v in combined.items() if v)


def solve(H,K,s,m,basis):
    stamp=monotonic()
    rate=K*(H-1)+2
    unit=lcm(s.denominator,*range(2,H+1))
    den=int((rate+s)*unit)
    todo=[physical.newborn(i) for i in range(2,H+1)]
    rows={}
    seen=set(todo)
    while todo:
        t=todo.pop()
        outgoing=row(t,H,K,rate,unit,basis)
        rows[t]=outgoing
        for nxt,coefficient in outgoing:
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    states=list(rows)
    index={t:i for i,t in enumerate(states)}
    transitions=[[(index[u],p) for u,p in rows[t]] for t in states]
    del rows
    expected=(H-1)*sum(triple[1 if basis else 2] for triple in counts(H,K))
    assert len(states)==expected,(len(states),expected)
    x=[0]*len(states)
    z=[1]*len(states)
    source=[int(virgin(t))*unit//2 for t in states]
    power=1
    for _ in range(m):
        xn=[]
        zn=[]
        for edges,b in zip(transitions,source):
            xi=b*power
            zi=0
            for j,p in edges:
                xi+=p*x[j]
                zi+=p*z[j]
            xn.append(xi)
            zn.append(zi)
        x,z=xn,zn
        power*=den
    assert min(x)>=0 and min(z)>=0
    leaves=max(1,K-1) if basis else 1
    theta=F(leaves*max(z),power)
    J=index[physical.newborn(2)]
    low=F(x[J],power)
    print('representation','ancestral paths' if basis else 'all labels','states',len(states),'edges',
          sum(map(len,transitions)),'negative coefficients',sum(p<0 for edges in transitions for j,p in edges),
          'seconds',round(monotonic()-stamp,3))
    if theta<1:
        delta=F(z[J],power)*F(leaves*max(x),power)/(1-theta)
        print('m-step norm <=',ceil(theta*10**12),'/ 1e12, killed scalar in units 1e-12',
              floor(low*10**12),ceil((low+delta)*10**12))
    else:
        print('m-step norm bound >=1, killed scalar lower in units 1e-12',floor(low*10**12))
    return {t:(x[i],z[i]) for i,t in enumerate(states)}


if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('--degree',type=int,default=4)
    p.add_argument('--cores',type=int,default=3)
    p.add_argument('--steps',type=int,default=96)
    p.add_argument('--mode',choices=('both','basis','full'),default='both')
    args=p.parse_args()
    assert 2<=args.degree<=19 and 1<=args.cores
    s=F('.56564')
    print('s',s,'degree cap',args.degree,'core cap',args.cores,'uniformized steps',args.steps)
    if args.mode in ('both','basis'):
        compressed=solve(args.degree,args.cores,s,args.steps,True)
    if args.mode in ('both','full'):
        full=solve(args.degree,args.cores,s,args.steps,False)
    if args.mode=='both':
        for t,values in full.items():
            assert values==tuple(sum(p*compressed[y][j] for y,p in expansion(t)) for j in (0,1))
        print('all positive path sums and remaining unit-mass numerators agree exactly',2*len(full))
