#!/usr/bin/env python3
"""Local physical count corrections, exact replay and a blind-cell control."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import cache
from itertools import combinations, product
from math import comb, factorial, prod
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from replay import CAP, FIXED, FLOOR, NB, S, Cert, Rows, add, audit, catalog, cert, children, independent, raw, stamp

# Frozen from the previous sum-one global contrast fit. No refitting here.
W = tuple(map(Q, ('.00190648','.28049101','.03306576','.68453675')))
GRID = 10**14
LIMIT = 25


@cache
def occurrences(p):
    """Rooted connected subsets, once each, of at most five physical cores."""
    ch=children(p)
    def down(i):
        v=[(i,)]
        for j in ch[i]:
            a=down(j)
            v += [x+y for x in v for y in a if len(x)+len(y)<=5]
        return v
    return tuple(ids for i in range(len(p)) for ids in down(i))


def canon(p,x,ids):
    keep=set(ids)
    def encode(i):
        return x[i], tuple(sorted(encode(j) for j in ids if p[j]==i and j in keep))
    return encode(ids[0])


def count(p,x,allow=None):
    out=Counter()
    for ids in occurrences(p):
        if all(sum(x[i])<=LIMIT for i in ids):
            cell=canon(p,x,ids)
            if allow is None or cell in allow:
                out[cell]+=1
    return out


def score(p,ls,allow=None):
    ch=children(p)
    x=tuple((l,len(k),int(i>0)) for i,(l,k) in enumerate(zip(ls,ch)))
    ds=tuple(map(sum,x))
    old=count(p,x,allow)
    flow=defaultdict(Q)
    total=Q(0)
    def step(rate, new):
        nonlocal total
        total+=rate
        for c,n in new.items():flow[c]+=rate*n
    for i,(l,k,h) in enumerate(x):
        ns=ch[i]+([p[i]] if h else [])
        step(l+sum(Q(1,ds[j]) for j in ns),count(p,add(x,i,(l+1,k,h)),allow))
        if l:
            step(Q(l,ds[i]),count(p+(i,),add(x,i,(l-1,k+1,h))+(NB,),allow))
    for c,n in old.items():flow[c]-=(total+S)*n
    return {c:v/len(p) for c,v in flow.items() if v},old


def atom_count(tree, allow):
    ch,x,ids=audit.physical_view(tree)
    p=(-1,*(ids.index(tree[v]) for v in ids[1:]))
    total=Counter()
    # Independent subset enumeration and connectivity test (one top role).
    for k in range(1,min(5,len(p))+1):
        for a in combinations(range(len(p)),k):
            if any(p[j] not in a for j in a[1:]) or any(sum(x[j])>LIMIT for j in a):
                continue
            q=canon(p,x,a)
            if q in allow:total[q]+=1
    return total


def marked_flow(p,ls,z):
    """Pre-event marked occurrences; actual outside arrows remain separate."""
    ch=children(p)
    x=tuple((l,len(k),int(i>0)) for i,(l,k) in enumerate(zip(ls,ch)))
    get=lambda p,y,ids: z.get(canon(p,y,ids),Q(0)) if all(sum(y[i])<=LIMIT for i in ids) else Q(0)
    amounts=defaultdict(Q)
    for ids in occurrences(p):
        old=get(p,x,ids)
        amounts['discount']-=S*old
        for i in ids:
            l,k,h=x[i]
            delta=get(p,add(x,i,(l+1,k,h)),ids)-old
            amounts['leaf gain']+=l*delta
            for donor in ch[i]+([p[i]] if h else []):
                amounts['internal gain' if donor in ids else 'external gain']+=delta/sum(x[donor])
            if l:
                rate=Q(l,sum(x[i]))
                post=add(x,i,(l-1,k+1,h))
                amounts['retained promotion']+=rate*(get(p,post,ids)-old)
                if len(ids)<5:
                    amounts['extension birth']+=rate*get(p+(i,),post+(NB,),ids+(len(p),))
    newborn=z.get((NB,()),Q(0))
    amounts['singleton birth']+=sum(Q(l,sum(s)) for s in x for l in s[:1])*newborn
    return amounts


def verify(p,ls,allow,z):
    tree=audit.build(p,ls)
    def val(t):return sum(z[c]*n for c,n in atom_count(t,allow).items())
    old=val(tree)
    result=sum(r*(val(tree+(j,))-old) for j,r in audit.arrows(tree).items())-S*old
    g,occ=score(p,ls,allow)
    assert (result/len(p),old)==(sum(z[c]*a for c,a in g.items()),sum(z[c]*n for c,n in occ.items()))
    assert result==sum(marked_flow(p,ls,z).values())


def multiplicity(c):
    d=comb(c[0][1],len(c[1]))*factorial(len(c[1]))
    for child,n in Counter(c[1]).items():
        # Divide the multinomial after multiplying all factorials.
        d*=multiplicity(child)**n
    return d//prod(factorial(n) for n in Counter(c[1]).values())


def fit(rows, y, cells, arity, cap, live=False, budget=None):
    zlist=[c for c in cells if arity(c)<=cap and (not live or c[0][2])]
    a=np.array([[r.get(c,0) for c in zlist] for r in rows],dtype=float)
    b=np.array(y,dtype=float)*10**6
    m=len(zlist)
    a=np.column_stack((a,-a,np.ones(len(rows))))
    ub=np.column_stack((np.vstack((a,-a)),-np.ones(2*len(rows))))
    rhs=np.r_[-b,b]
    box=[(0,None)]*(2*m)+[(None,None),(0,None)]
    cost=[multiplicity(c) if c[0][2] else 0 for c in zlist]*2+[0,0]
    t=linprog([0]*(2*m+1)+[1],A_ub=np.vstack((ub,cost)) if budget else ub,
              b_ub=np.r_[rhs,float(budget)*10**6] if budget else rhs,bounds=box,method='highs')
    assert t.success,t.message
    # Among near-minimax fits, choose a sparse min-l1 correction. Both
    # minimizations are proposals; all published values are replayed in Q.
    u=linprog([1]*(2*m)+[0],A_ub=np.vstack((a,-a,cost[:-1])) if budget else np.vstack((a,-a)),
              b_ub=np.r_[rhs+t.fun+1e-6,float(budget)*10**6] if budget else rhs+t.fun+1e-6,
              bounds=[(0,None)]*(2*m)+[(None,None)],method='highs')
    assert u.success,u.message
    grid=lambda x:Q(round(float(x)*GRID/10**6),GRID)
    return {c:grid(v) for c,v in zip(zlist,u.x[:m]-u.x[m:2*m]) if grid(v)},grid(u.x[-1])


def size(c):return 1+sum(map(size,c[1]))


def join(w,row):return sum(a*b for a,b in zip(w,row))


def test_family():
    """Not in discovery catalog: two adjacent branching centres, asymmetric tips."""
    p=(-1,0,0,0,1,1,4,4)
    for a,b,c in product((0,2,5),(0,1,4),(1,3,8)):
        yield 'broom',p,(a,b,c,c,b,c,c,c+1)
    for n in (4,5,7,8):
        for label,p in (('blind chain',(-1,*range(n-1))),('blind star',(-1,*(0,)*(n-1))),
                        ('blind shifted star',(-1,0,*(1,)*(n-2)))):
            for a,b in ((7,9),(8,12),(10,7),(12,13),(11,11),(7,13)):
                yield label,p,tuple((a,b)[i%2] for i in range(n))
    for n in (3,5,8):
        for label,p in (('overflow chain',(-1,*range(n-1))),('overflow star',(-1,*(0,)*(n-1)))):
            for a,b in ((27,38),(51,77),(64,101),(26,161),(65,65),(99,99)):
                yield label,p,tuple((a,b)[i%2] for i in range(n))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--blobs',nargs=4,type=Path,default=list(map(Path,('b3f18_lo.bin','b3f18_hi.bin','b3f20_lo.bin','b3f20_hi.bin'))))
    args=ap.parse_args()
    pins=list(map(stamp,args.blobs))
    cs=[*(audit.Cert(p) for p in args.blobs[:2]),*(Cert(p,s,verify=False) for p,s in zip(args.blobs[2:],('lower','upper')))]
    assert [(c.D,c.mode,c.S) for c in cs]==[(18,2,10**11),(18,1,10**11),(20,2,10**11),(20,1,10**11)]
    es=list(map(Rows,cs))
    cat=list(dict.fromkeys(catalog()))
    train=[cat[i] for i in sorted(set().union(*FIXED.values()))]
    test=list(dict.fromkeys(test_family()))
    assert not ({(p,ls) for _,p,ls in cat} & {(p,ls) for _,p,ls in test})
    def four(seq):return [tuple(raw(c,e,p,ls)[1]/len(p) for c,e in zip(cs,es)) for _,p,ls in seq]
    y,yt=four(train),four(test)
    gx=[score(p,ls)[0] for _,p,ls in train]
    cells=sorted(set().union(*gx))
    print('trial',S,'train',len(train),'adversarial',len(test),'cell columns',len(cells),flush=True)
    allow=set(cells)
    gy=[score(p,ls,allow)[0] for _,p,ls in test]
    assert all(cells)
    mass=join(W,FLOOR)
    diameter=lambda v:(max(v)-min(v))/2
    u=[join(W,r) for r in y]; v=[join(W,r) for r in yt]
    pairs=[]
    for i,j in combinations(range(len(train)),2):
        slope=max((abs(gx[i].get(c,0)-gx[j].get(c,0))/multiplicity(c)
                   for c in cells if c[0][2]), default=Q(0))
        pairs.append(((abs(u[i]-u[j])-slope/50)/2,i,j,slope))
    lower,i,j,slope=max(pairs)
    print('exact budget-1/50 pair barrier',float(lower),'pair',i,j,'max contrast per budget',float(slope),flush=True)
    for cap,live,budget in ((1,False,None),(3,False,None),(5,False,None),(5,True,None),(5,True,Q('.02'))):
        z,a=fit(gx,[join(W,r) for r in y],cells,size,cap,live,budget)
        f=lambda t,g:t+sum(z.get(c,0)*x for c,x in g.items())
        trainval=[f(t,g) for t,g in zip(u,gx)]
        testval=[f(t,g) for t,g in zip(v,gy)]
        # min comparison of corrected denominator using a cap on occurrences
        # per assigned top core, computed only for cells with nonzero value.
        sup=sum(abs(val)*multiplicity(c) for c,val in z.items() if c[0][2])
        print('order',cap,'nonfounder-only',live,'budget',budget,'used cells',len(z),'by size',dict(Counter(map(size,z))),
              'train raw/corrected',float(diameter(u)),float(diameter(trainval)),
              'holdout raw/corrected',float(diameter(v)),float(diameter(testval)),
              'frozen-offset holdout',float(max(abs(q+a) for q in testval)),flush=True)
        print(' proposal', 'offset',str(a),'coefficient l1',float(sum(map(abs,z.values()))),
              'absolute stationary denominator budget',float(sup),'raw denominator floor',float(mass),flush=True)
        if budget:
            assert sup<budget+Q(1,10**9)<mass
            for case in train+test[:4]+test[-3:]:verify(case[1],case[2],set(z),z)
            assert sum(marked_flow(train[0][1],train[0][2],z).values())==sum(z.get(c,0)*v for c,v in gx[0].items())*len(train[0][1])
            print('worst heldout',test[max(range(len(test)),key=lambda i:abs(testval[i]+a))],flush=True)
    # No retained census is visible on these overflow trees before a jump.
    # A newborn singleton is visible after promotions; retain its true rate
    # as an extra column and annihilate it as well as the per-core offset.
    newborn=(NB,())
    blind=[r+(g.get(newborn,Q(0)),) for r,g in zip(yt,gy) if not set(g)-{newborn}]
    assert blind
    blind_cases=[(p,ls,r) for (_,p,ls),g,r in zip(test,gy,yt) if not set(g)-{newborn}]
    costs=[]
    support=set()
    for target in product(*zip(FLOOR,CAP)):
        norm,ids,_=cert(blind,target+(Q(0),Q(0)),Q('1e-8'),True)
        support.update(ids)
        costs.append(norm)
    for i in support:
        p,ls,r=blind_cases[i]
        assert all(independent(c,p,ls)[1]/len(p)==v for c,v in zip(cs,r))
    n=(10**12/max(costs))//1
    assert n>10**4 and Q(n,10**12)*max(costs)<=1
    print('blind physical trees',len(blind),'exact corner count',len(costs),'obstruction 1e-12 units',n,
          'fixed-mixture blind halfwidth before newborn correction',float(diameter([join(W,r) for r in blind])),flush=True)
    # A fixed-mixture barrier after arbitrary newborn and constant corrections.
    wnorm,ids,_=cert(blind,(Q(1),Q(1),Q(1),Q(1),Q(0),Q(0)),Q('1e-8'),True)
    assert sum(W)==1 and wnorm
    for i in set(ids)-support:
        p,ls,r=blind_cases[i]
        assert all(independent(c,p,ls)[1]/len(p)==v for c,v in zip(cs,r))
    print('blind support union',len(support|set(ids)),'fixed-gauge barrier 1e-12 units',(10**12/wnorm)//1,flush=True)
    assert pins==list(map(stamp,args.blobs))


if __name__=='__main__':main()
