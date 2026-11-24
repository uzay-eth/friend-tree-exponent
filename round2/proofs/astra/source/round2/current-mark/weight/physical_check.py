#!/usr/bin/env python3
"""Cold labelled event enumeration for the current-mark weight."""
from fractions import Fraction as F
from functools import cache

from cert import propose, prove, RHO, DAMP, basic

CAP=36
J=(1,())
Z=propose(CAP,DAMP,RHO,True)
H,*_=prove(CAP,DAMP,RHO,Z)
z,k=basic(CAP,DAMP)


def tree(l,cs): return l,tuple(sorted(cs))
def deg(t): return 1+t[0]+len(t[1])
def plus(t): return t[0]+1,t[1]
def q(D): return F(1,min(D,CAP))


@cache
def w(t):
    return H[min(CAP,deg(t)),min(CAP,1+len(t[1]))]+Z[min(CAP,deg(t))]*sum(map(w,t[1]))


def intrinsic(t):
    nodes={}
    def walk(x,p=()):
        nodes[p]=x
        for i,y in enumerate(x[1]): walk(y,p+(i,))
    walk(t)
    def change(x,p,bud=False):
        l,cs=x
        if not p: return tree(l-1,cs+(J,)) if bud else plus(x)
        i=p[0]
        return tree(l,cs[:i]+(change(cs[i],p[1:],bud),)+cs[i+1:])
    for p,x in nodes.items():
        l,cs=x
        if l:
            yield F(l),change(t,p)
            yield F(l,deg(x)),change(t,p,True)
        for dst in ([p[:-1]] if p else [])+[p+(i,) for i in range(len(cs))]:
            yield F(1,deg(x)),change(t,dst)


def gains(t): return t[0]+sum(F(1,deg(y)) for y in t[1])


def active(t,a):
    u=a+DAMP
    p=F(1,deg(t)); hi=q(deg(t));lo=q(deg(t)+1)
    def dormant(y,l,h): return (1/(l+DAMP)-1/(h+DAMP))*w(y)
    total=sum(r*(w(v)-w(t)) for r,v in intrinsic(t))/u-a*w(t)/u
    for y in t[1]:
        total+=((hi-p)*w(y)+p*w(plus(y)))/(hi+DAMP)
        total+=gains(t)*dormant(y,lo,hi)
        a1=q(deg(y)+1);a0=q(deg(y))
        for x in y[1]: total+=p*dormant(x,a1,a0)
    return total-RHO*w(t)/u


def main():
    states={J}
    for _ in range(5): states.update(v for t in tuple(states) for _,v in intrinsic(t))
    states.update([tree(0,(J,)*j) for j in (2,9,18,34,35,36,60)])
    states.update([tree(1,(tree(0,(J,)*j),)) for j in (2,18,34,35,36,42)])
    states.add((0,((0,((0,(J,)),)),)))
    for t in states:
        for a in (F(1,CAP),F(1,2)):
            assert active(t,a)<=0,(t,a,active(t,a))
        l,u=F(1,CAP),F(1,CAP-1)
        b=1/(l+DAMP)-1/(u+DAMP)
        d=sum(r*(w(v)-w(t)) for r,v in intrinsic(t))*b-u*b*w(t)+(u-l)*w(t)/(l+DAMP)
        assert d<=RHO*b*w(t)
        bad=d+(u-l)*b*w(t)
        # The activation consumes the old dormant; this is its exact drift
        # shift if the disappearance term is omitted.
        assert bad-d==(u-l)*b*w(t)>0
    print('cold active endpoints and dormant activation:',len(states),'shapes')
    for fault in ('wrong-send','cap20-rates'):
        try:
            prove(CAP,DAMP,RHO,Z,fault)
        except ValueError as e:
            assert 'edge corner' in str(e)
            print('rejected',fault,str(e).split(':')[0])
        else: raise AssertionError('planted defect accepted')


if __name__=='__main__': main()
