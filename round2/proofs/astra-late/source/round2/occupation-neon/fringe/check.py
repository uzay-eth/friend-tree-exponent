#!/usr/bin/env python3
"""Independent physical drift and all-eight allocation regression."""
import argparse
from functools import cache
from pathlib import Path
from charges import Cert, Charges, NB, SPEC, Q, shape, label, up, bud


def tree(l,*kids):return (0,l,tuple(sorted(kids)))


def reach(cap):
    @cache
    def nxt(z):
        _,l,k=z; out={tree(l+1,*k)}
        if l:out.add(tree(l-1,*k,tree(1)))
        for i,y in enumerate(k):
            out.update(tree(l,*k[:i],yy,*k[i+1:]) for yy in nxt(y))
        return out
    layer={tree(2)}
    for _ in range(3,cap+1):
        yield from sorted(layer)
        layer={y for x in layer for y in nxt(x)}


def drift(c,z):
    p,x=shape(z,0); k=[[] for _ in p]
    for i,j in enumerate(p[1:],1):k[j].append(i)
    f=c.potential
    old=f(k,x); a=-c.lam*old
    for i,y in enumerate(x):
        L=y[0]
        a+=(L+sum(Q(1,sum(x[j])) for j in k[i]+([p[i]] if i else [])))*(f(k,label(x,i,up(y)))-old)
        if L:
            kk=[list(ch) for ch in k]+[[]];kk[i].append(len(x))
            a+=Q(L,sum(y))*(f(kk,label(x,i,bud(y))+(NB,))-old)
    return a/c.S


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blobs',nargs=2,type=Path);a=p.parse_args()
    J=tree(1); three=tree(1,J,J,J); four=tree(1,tree(1,J),three,J)
    mixed=[tree(17,four,J),tree(1,tree(18,four,J),tree(16,J)),
           tree(17,tree(18,three,J),tree(19,J)),tree(14001,tree(14000,J)),
           tree(20,*([J]*19)),tree(1,tree(14000,*([J]*18)),four)]
    controls=list(reach(8))+mixed

    class Omit(Charges):
        def transferred(self,f,x,i,upward):
            r,q=super().transferred(f,x,i,upward)
            if f=='g' and sum(x[i])>=18:return Q(0),Q(0)
            return r,q

        def terms(self,f,x):
            # One-sided mutation: drop the transferred amount from the
            # source edge but leave the hV/hC receive credit intact.
            if f in ('hV','hC'):
                return pristine.terms(f,x)
            return super().terms(f,x)

    for b in a.blobs:
        c=Cert(b);pristine=Charges(c); mutant=Omit(c)
        seen=set(); bad=0
        for z in controls:
            rs=pristine.rows(z,0,False)
            assert not any(s for _,s in rs.values())
            signed=pristine.sign*drift(c,z)
            assert sum(r for r,_ in rs.values())==signed
            seen.update(f for f,(r,_) in rs.items() if r)
        for z in mixed:
            try:
                rr=mutant.rows(z,0,False)
                bad+=sum(r for r,_ in rr.values())!=pristine.sign*drift(c,z)
            except AssertionError:bad+=1
        assert seen==set(SPEC) and bad>0
        print('sign',pristine.sign,'whole-tree exact physical equalities',len(controls),'all 8 nonzero','tail-transfer omissions caught',bad)


if __name__=='__main__':main()
