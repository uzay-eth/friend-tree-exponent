"""Exact shared-send and rectangular Bellman bounds on one-high three-core exits."""
from argparse import ArgumentParser
from ast import literal_eval
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import ceil,floor
from pathlib import Path

from sector import BETA, UNIT as MASS, geometry, hy, info, mark, newborn, physical, plus

UNIT=10**12
S=F('.565809')
CAP=0


def nodes(tree,path=()):
    yield path,tree
    for i,y in enumerate(tree[2]):
        yield from nodes(y,path+(i,))


def replace(tree,path,value):
    if not path:
        return value
    tag,l,ch=tree
    i=path[0]
    updated=replace(ch[i],path[1:],value)
    return tag,l,tuple(sorted((*ch[:i],updated,*ch[i+1:])))


@cache
def leaf_bounds(t):
    lo,hi=hy.ex.boundary(t,BETA)
    return floor(UNIT*lo),ceil(UNIT*hi)


def extremum(numer,den,choices,upper):
    # A linear-fractional box extremum occurs at a sorted prefix of the
    # switch ratios. Each choice is scaled to endpoints zero and one.
    ratios=[]
    for a,b in choices:
        if b:
            ratios.append((a/b,a,b))
        elif (a>0 if upper else a<0):
            numer+=a
    for ratio,a,b in sorted(ratios,reverse=upper):
        if (a*den>numer*b if upper else a*den<numer*b):
            numer+=a
            den+=b
        else:
            break
    return numer/den


def row_terms(tree,val):
    all_nodes=dict(nodes(tree))
    gain={v:val(replace(tree,v,plus(t))) for v,t in all_nodes.items()}
    d,n=info(tree)
    den=S+mark(tree)
    N=[F(UNIT*tree[1]*hy.rows.g(BETA,d),d)]*2
    terms=[]
    for path,t in all_nodes.items():
        D,_=info(t)
        bud=F(t[1],D)
        den+=t[1]+bud
        if t[1]:
            fresh=newborn(2 if D<16 else 16)
            promoted=t[0],t[1]-1,tuple(sorted((*t[2],fresh)))
            after=val(replace(tree,path,promoted))
            for k in (0,1):
                N[k]+=t[1]*gain[path][k]+bud*after[k]
        if path:
            parent=path[:-1]
            for sender,receiver in ((path,parent),(parent,path)):
                rate=F(1,info(all_nodes[sender])[0])
                terms.append((sender,rate*gain[receiver][0],rate*gain[receiver][1],rate))
    for child in tree[2]:
        unshifted,shifted=val(child),val(plus(child))
        for k in (0,1):
            N[k]+=mark(child)*unshifted[k]
        terms.append(((),F(shifted[0]-unshifted[0],d),F(shifted[1]-unshifted[1],d),F(0)))
    assert den+sum(x[3] for x in terms)==S+mark(tree)+n-F(1,d)
    return N,den,terms


@cache
def backup(tree,depth,mode):
    old=leaf_bounds(tree)
    if not depth or geometry(tree)[1:]!=(3,1) or CAP and geometry(tree)[0]>CAP:
        return old
    def val(t):
        return backup(t,depth-1 if depth>0 else -1,mode)
    N,den,terms=row_terms(tree,val)
    if mode=='tied':
        grouped=defaultdict(lambda:[F(0),F(0),F(0)])
        for sender,a,b,c in terms:
            row=grouped[sender]
            row[0]+=a
            row[1]+=b
            row[2]+=c
        terms=[(v,*row) for v,row in grouped.items()]
    lower,upper=old
    for k in (0,1):
        choices=[(term[k+1],term[3]) for term in terms]
        if mode=='actual':
            x=(N[k]+sum(a for a,b in choices))/(den+sum(b for a,b in choices))
        else:
            x=extremum(N[k],den,choices,k==1)
        if k:
            upper=min(upper,ceil(x))
        else:
            lower=max(lower,floor(x))
    assert lower<=upper
    return lower,upper


if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('--depth',type=int,default=1)
    p.add_argument('--part',type=int,default=0)
    p.add_argument('--parts',type=int,default=1)
    p.add_argument('--modes',default='box,tied,actual')
    p.add_argument('--cap',type=int,default=0)
    p.add_argument('--position',choices=('all','root','below'),default='all')
    p.add_argument('--tag',type=int,choices=(0,2,16),default=0)
    p.add_argument('--mass',type=Path,default=Path('round2/outputs/finite-low-entrances.tsv'))
    args=p.parse_args()
    CAP=args.cap
    depth=-1 if CAP else args.depth
    assert 0<=CAP<=31 and (CAP or 1<=depth<=15)
    hy.rows.S=hy.ex.S=hy.tc.S=S
    hy.rows.direct_pair=hy.tc.improved_direct_pair
    total=defaultdict(lambda:[0,0])
    marginal=defaultdict(lambda:[0,0])
    count=0
    for i,line in enumerate(args.mass.read_text().splitlines()[1:]):
        if i%args.parts!=args.part:
            continue
        l,h,state=line.split('\t')
        tree=literal_eval(state)
        mass=int(l)
        pos='root' if info(tree)[0]>6 else 'below'
        if args.position not in ('all',pos) or args.tag not in (0,tree[0]):
            continue
        old=leaf_bounds(tree)
        last=old
        for mode in args.modes.split(','):
            lo,hi=backup(tree,depth,mode)
            assert last[0]<=lo<=hi<=last[1]
            acc=total[pos,mode]
            # The original boundary uses exact fractions before its display
            # rounding. Discard one unit before claiming a gain against it.
            acc[0]+=mass*max(0,lo-old[0]-1)
            acc[1]+=mass*max(0,old[1]-hi-1)
            step=marginal[pos,mode]
            step[0]+=mass*max(0,lo-last[0]-1)
            step[1]+=mass*max(0,last[1]-hi-1)
            last=lo,hi
        count+=1
    print('depth',depth,'cap',CAP,'part',args.part,'of',args.parts,'states',count)
    for (pos,mode),(l,u) in sorted(total.items()):
        print('gain',pos,mode,'lower',l//MASS,'upper',u//MASS,'units 1e-12')
        dl,du=marginal[pos,mode]
        print('extra',pos,mode,'lower',dl//MASS,'upper',du//MASS,'units 1e-12')
