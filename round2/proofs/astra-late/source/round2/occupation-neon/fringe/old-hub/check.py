#!/usr/bin/env python3
"""Compare the first-arm toll and physical hazards on full high-degree trees."""
import sys
from pathlib import Path
from collections import defaultdict

from pilot import Cert, Charges, NAMES, NB, Q, Toll, arms, shape, v

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'higher'/'older'))
from replay import Local


def tree(l,*kids):return (2,l,tuple(sorted(kids)))


def fixtures(blob):
    c=Charges(Cert(blob)); independent=Local(c.c)
    tip=tree(1); old=tree(2); tall=tree(5); bare=tree(0,tip)
    arm_fixtures=(tree(0,bare),tree(0,tree(0,old)),tree(1,bare,old),
                  tree(0,old,old),tree(0,tree(0,tree(0,old)),old))
    other_fixtures=((tip,tip,tip,tip),(tip,old,bare,tall),(tall,tall,old,bare,tip))
    arcstates,arcs=arms(9)
    assert set(arm_fixtures)<=set(arcstates)
    def plain(x):return tree(x[1],*(plain(y) for y in x[2]))
    check=Toll(c,arm_fixtures); seen=set(); n=0; changed=0
    for d in (65,79,100):
        for outsiders in other_fixtures:
            k=len(outsiders)+1
            for x in arm_fixtures:
                for parent_degree in (2,10,25):
                    host=(2,d-k-1,(x,*outsiders))
                    whole=(2,parent_degree-2,(host,))
                    p,s=shape(whole,1)
                    ch=[[i for i,j in enumerate(p) if j==a] for a in range(len(p))]
                    expected=defaultdict(Q)
                    for f,ids in independent.occurrences(ch,s):
                        if f not in NAMES or ids[:2]!=(1,2):continue
                        st=tuple(s[i] for i in ids)
                        if NB in st or max(map(sum,st))>d:continue
                        ext=[(Q(1,sum(s[p[i]])) if p[i]>=0 and p[i] not in ids else Q(0),
                              sum(Q(1,sum(s[j])) for j in ch[i] if j not in ids)) for i in ids]
                        q=independent.charge(f,st,ext)
                        assert q>=0
                        expected[f]+=q
                        if q:seen.add(f)
                    u=Q(1,parent_degree)+sum(Q(1,v.degree(y)) for y in outsiders)
                    for f,(a,b) in zip(NAMES,check.row(d,k,x)):
                        assert a+b*u==expected[f],(f,d,k,x,u,a,b,expected[f])
                        changed+=b!=0
                    # Sum exact full-tree arrows by their observable effect.
                    # The verifier sorts successor siblings; use a unique
                    # arm label solely to recover the occurrence being tracked.
                    tagged=(3,x[1],x[2]); full=(2,parent_degree-2,((2,d-k-1,tuple(sorted((tagged,*outsiders)))),))
                    projected=defaultdict(Q)
                    for y,r in v.row(full)[3].items():
                        h=y[2][0]; child,=(a for a in h[2] if a[0]==3)
                        state=(v.degree(h),len(h[2]),plain(child))
                        if state!=(d,k,x):projected[state]+=r
                    want=defaultdict(Q,{(d+1,k,x):d-k-1+u+Q(1,v.degree(x)),
                                       (d,k+1,x):Q(d-k-1,d)})
                    want[d,k,v.plus(x)]+=Q(1,d)
                    for y,r in arcs[x].items():want[d,k,y]+=r
                    # A host promotion can produce permanent label 20;
                    # newborns inside the tracked arm also need de-marking.
                    # Tagged root is unique since its host never has D=3.
                    assert projected==want,(d,k,x,projected,want)
                    n+=1
    assert seen==set(NAMES) and changed
    print(c.sign,n,'full-crown charge and hazard fixtures; all three families nonzero')


if __name__=='__main__':
    for b in sys.argv[1:]:fixtures(Path(b))
