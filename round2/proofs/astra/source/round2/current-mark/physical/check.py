#!/usr/bin/env python3
"""Exact transfer checks, including a current-vs-permanent kernel witness."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import chain


def q(d):
    return F(1,min(d,20))


def tree(l,kids=()):
    return l,tuple(sorted(kids))


J=tree(1)


def gain(t):
    return tree(t[0]+1,t[1])


def d(t,root=1):
    return root+t[0]+len(t[1])


@cache
def physical(t,slot=1):
    """Direct per-arrow enumeration on chronological paths, no mark metadata."""
    nodes={}
    def visit(t,k):
        nodes[k]=t
        for i,y in enumerate(t[1]):
            visit(y,k+(i,))
    visit(t,())

    def change(k,fn,x=t):
        if not k:
            return fn(x)
        l,ys=x
        i=k[0]
        return tree(l,(*ys[:i],change(k[1:],fn,ys[i]),*ys[i+1:]))

    r=defaultdict(F)
    for k,x in nodes.items():
        D=d(x,slot if not k else 1)
        if x[0]:
            r[change(k,gain)]+=x[0]
            r[change(k,lambda u:tree(u[0]-1,(*u[1],J)))]+=F(x[0],D)
        if k:
            r[change(k[:-1],gain)]+=F(1,D)
        for i in range(len(x[1])):
            r[change(k+(i,),gain)]+=F(1,D)
    return dict(r)


def active(a,t):
    return a,a,t


def entries(t,slot=1,emit=True,scaled=False,plant=None):
    """Independent births attached to an active or physical centre."""
    D=d(t,slot)
    a=q(D)
    a1=q(D+1)
    p=F(1,D)
    r0=t[0]+sum(F(1,d(y)) for y in t[1])
    out=defaultdict(F)
    if emit and t[0]:
        out[active(a,J)]+=F(t[0],D)

    def dormant(lo,hi,y,r):
        if hi>lo and r:
            out[lo,hi,y]+=r*((hi-lo) if scaled else 1)

    for y in t[1]:
        out[active(a,gain(y))]+=p
        if a>p:
            out[active(a,y)]+=a-p
        if plant!='omit-root-drop':
            dormant(a1,a,y,r0)
        if plant!='omit-grandchild':
            for z in y[1]:
                D_y=d(y)
                dormant(q(D_y+1),q(D_y),z,p if plant!='reverse-rate' else F(1,D_y))
    return out


def transfer_tests():
    def test(t,centre,plant=None):
        def fn(a,x):
            return F((x[0]+2)**2,1+len(x[1]))+a*(1+x[0]+sum(y[0] for y in x[1]))

        def wrap(a,b,x):
            return fn(a,x) if a==b else fn(a,x)-fn(b,x)

        def total(t,a,slot):
            return (fn(a,t) if slot else 0)+sum(total(y,q(d(t,slot)),1) for y in t[1])

        def rhs(t,a,slot):
            result=sum(r*wrap(*z) for z,r in entries(t,slot,plant=plant).items())
            if slot:
                result+=sum(r*(fn(a,y)-fn(a,t)) for y,r in physical(t).items())-a*fn(a,t)
            return result+sum(rhs(y,q(d(t,slot)),1) for y in t[1])

        slot=int(not centre)
        a=F(0) if centre else F(1,7)
        lhs=sum(r*(total(y,a,slot)-total(t,a,slot)) for y,r in physical(t,slot).items())
        expected=rhs(t,a,slot)+(a*fn(a,t) if slot else 0)
        return lhs==expected

    levels={J}
    corpus=set(levels)
    for _ in range(5):
        levels={y for t in levels for y in physical(t)}
        corpus|=levels
    # All centre starts must have degree at least two. Include exact cap seam
    # and two equal old subtrees, each with a grandchild.
    branch=tree(1,(J,))
    corpus |= {tree(k,(branch,branch)) for k in (2,16,17,18,19,20,23)}
    errors=Counter()
    cases=0
    for t in corpus:
        for centre in (False,True):
            if centre and d(t,0)<2:
                continue
            assert test(t,centre)
            cases+=1
            for p in ('omit-root-drop','omit-grandchild','reverse-rate'):
                errors[p]+=not test(t,centre,p)
    assert min(errors.values())>0
    print('physical centre/plant transfers',cases,'mutations rejected',dict(errors))


def residual(x,f,scaled=False,plant=None):
    a,b,t=x
    L=sum(r*(f((a,b,y))-f(x)) for y,r in physical(t).items())
    if a!=b:
        if scaled:
            return L-b*f(x)+f(active(a,t))
        return L-(a if plant=='wrong-total-departure' else b)*f(x)+(b-a)*f(active(a,t))
    return L-a*f(x)+sum(r*f(y) for y,r in entries(t,emit=False,scaled=scaled).items())


def dormant_tests():
    t=tree(2,(tree(1,(J,)),tree(2)))
    a,b=F(1,4),F(1,3)
    x=(a,b,t)

    @cache
    def moment(k,x,scaled=False,plant=None):
        a,b,t=x
        return residual(x,lambda y:moment(k-1,y,scaled,plant),scaled,plant) if k else F(a==b)*F(t[0],d(t))

    rejected=0
    for k in range(5):
        delta=moment(k,active(a,t))-moment(k,active(b,t))
        assert delta==moment(k,x)==(b-a)*moment(k,x,True)
        assert moment(k,active(a,t))==moment(k,active(a,t),True)
        rejected+=delta!=moment(k,x,False,'wrong-total-departure')
    assert rejected
    print('dormant differences/scaled derivatives',5,'wrong total departure rejected',rejected)


@cache
def permanent(t,parent=0):
    """Labelled trees: root tag plus unmarked local leaves and tagged children."""
    return parent,t[0],tuple(permanent(y,parent) for y in t[1])


def pd(x):
    return 1+x[1]+len(x[2])


def pj(a):
    return a,1,()


def pp(x):
    return x[0],x[1]+1,x[2]


@cache
def pmoves(x,input=F(0)):
    a,l,k=x
    p=F(1,pd(x))
    out=defaultdict(F)
    rate=l+input+sum(F(1,pd(y)) for y in k)
    if rate:
        out[pp(x)]+=rate
    if l:
        out[a,l-1,tuple(sorted((*k,pj(q(pd(x))))) )]+=l*p
    for j,y in enumerate(k):
        for yy,r in pmoves(y,p).items():
            out[a,l,tuple(sorted((*k[:j],yy,*k[j+1:])))]+=r
    return dict(out)


def derivatives():
    @cache
    def old(n,x):
        a,l,k=x
        p=F(1,pd(x))
        if not n:
            return l*p
        v=old(n-1,x)
        return (sum(r*(old(n-1,y)-v) for y,r in pmoves(x).items())-a*v
                +sum((y[0]-p)*old(n-1,y)+p*old(n-1,pp(y)) for y in k))

    @cache
    def new(n,x):
        a,b,t=x
        if not n:
            return F(t[0],d(t)) if a==b else F(0)
        return residual(x,lambda y:new(n-1,y))

    aa=F(1,2)
    ds=[]
    for j in range(7):
        ds.append(new(j,active(aa,J))-old(j,pj(aa)))
    assert any(ds)
    print('new minus original newborn unit-reward derivatives, orders 0..6:',','.join(map(str,ds)))


if __name__=='__main__':
    transfer_tests()
    dormant_tests()
    derivatives()
