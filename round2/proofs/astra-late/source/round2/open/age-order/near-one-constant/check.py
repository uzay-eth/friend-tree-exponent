#!/usr/bin/env python3
"""Exact finite controls for first-defect law and mixture normalization."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product


def neighbors(tree):
    a=[set() for _ in tree]
    for j,i in enumerate(tree[1:],1):
        a[i].add(j)
        a[j].add(i)
    return a


def partition(a,v):
    blocks={}
    for root in a[v]:
        seen={v,root}; todo=[root]
        for i in todo:
            fresh=a[i]-seen
            todo.extend(fresh)
            seen.update(fresh)
        if len(todo)>1:
            blocks[root]=set(todo)
    return blocks


def census(a,v):
    d=len(a[v]); n=len(a)
    blocks=partition(a,v)
    N=n-d
    A=sum((1-Q(1,len(a[y])) for y in blocks),Q(0))
    assert N==1+sum(len(B)-1 for B in blocks.values())
    assert 1+sum(len(B)-Q(1,len(a[y])) for y,B in blocks.items())==N+A
    return N,A,blocks


def physical():
    rows=arrows=0
    for n in range(2,7):
        for parents in product(*(range(i) for i in range(1,n))):
            seed=(-1,)+parents
            for v in range(n):
                original=neighbors(seed)
                N0,A0,old=census(original,v)
                beta=N0+A0
                for extra in (0,1,4):
                    tree=seed+(v,)*extra
                    a=neighbors(tree); N,A,blocks=census(a,v)
                    assert (N,A,blocks)==(N0,A0,old)
                    got=Counter(); weights=Counter()
                    for x,adj in enumerate(a):
                        for z in adj:
                            weight=Q(1,len(adj))
                            if z==v:
                                got['tag gain']+=weight
                                continue
                            post=neighbors(tree+(z,))
                            NN,AA,new=census(post,v)
                            assert NN==N+1 and AA>=A
                            assert NN+A0==beta+1
                            if x==v:
                                case='into core' if z in blocks else 'into leaf'
                                if case=='into leaf':
                                    assert new[z]=={z,len(tree)} and len(new)==len(blocks)+1
                                else:
                                    assert new[z]==blocks[z]|{len(tree)}
                            else:
                                root,=tuple(y for y,B in blocks.items() if x in B and z in B)
                                assert new[root]==blocks[root]|{len(tree)}
                                case='old internal'
                            weights[case]+=weight/beta
                            arrows+=1
                    d=len(a[v]); c=len(blocks)
                    assert got['tag gain']==len(tree)-beta
                    assert weights==Counter({'old internal':(beta-1)/beta,
                                             'into leaf':Q(d-c,d)/beta,'into core':Q(c,d)/beta})
                    assert sum(weights.values())==1
                    rows+=1
    print(f'PASS {rows} perfect prefixes and {arrows} directed outside arrows: deficit, Palm weights and postdefect shift')


def sums():
    runs=0
    for n0,alpha in ((2,Q(1)),(3,Q(1,2)),(4,Q(1,3)),(5,Q(5,2)),(6,Q(11,6))):
        beta=n0-alpha
        s=Q(1); total=Q(0)
        for m in range(n0,100):
            w=s*beta/m
            s*=1-beta/m
            total+=w
            assert total+s==1
            assert m-beta>0 and (beta+1)+(m-beta)==m+1
            if beta==1:
                assert w==Q(n0-1,m*(m-1))
            runs+=1
    # Integrate the exact limiting size kernel against a two-type amplitude.
    # Integer beta so all weighted powers and survival tests are rational.
    for b in (1,2,5):
        atoms=((Q(2,5),Q(1,2)),(Q(3,5),Q(3)))
        inv=sum(p*y**(-b) for p,y in atoms)
        assert sum(p*y**(-b)/inv for p,y in atoms)==1
        for r in (Q(1),Q(3,2),Q(4)):
            # Integral of b*u^(-b-1) above ry, divided by total mass.
            assert sum(p*(r*y)**(-b) for p,y in atoms)/inv==r**(-b)
    print(f'PASS {runs} exact hazard/telescoping steps and inverse-tilted Pareto normalization')


if __name__=='__main__':
    physical()
    sums()
