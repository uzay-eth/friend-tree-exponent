#!/usr/bin/env python3
"""Independent exact fixtures for the multicut likelihood and quantitative mesh."""
from fractions import Fraction as Q
from math import comb


def rates(adj, cut):
    good, cross = {}, {}
    for i, nbrs in enumerate(adj):
        for j in nbrs:
            table = cross if frozenset((i,j)) in cut else good
            table[i,j]=Q(1,len(nbrs))
    return good,cross


def cuts():
    # Three anchors each with a private real leaf; two slots at the middle.
    adj = [[1,3],[0,2,4],[1,5],[0],[1],[2]]
    cut={frozenset((0,1)),frozenset((1,2))}
    good,cross=rates(adj,cut)
    assert cross=={(0,1):Q(1,2),(1,0):Q(1,3),(1,2):Q(1,3),(2,1):Q(1,2)}
    assert sum(good.values())+sum(cross.values())==6
    roots=(0,1,2)
    owner=[0,1,2,0,1,2]
    w=Q(1)
    for i in (1,0,1,2,2,1,0,1):
        good,cross=rates(adj,cut)
        n=owner.count(i)
        k=sum(r for (src,_),r in cross.items() if owner[src]==i)
        lam=sum(r for (src,_),r in good.items() if owner[src]==i)
        assert lam+k==n and lam>0
        # Wait+cross-killing exponent is -n for every side, even the
        # uncompleted waits; completed jump multipliers are lambda/n.
        assert -(lam+k)==-n
        victim=roots[i]
        leaf=next(j for j in adj[victim] if len(adj[j])==1)
        physical=good[leaf,victim]
        ref=n*physical/lam
        assert physical/ref==1-k/n
        w*=physical/ref
        adj[victim].append(len(adj))
        adj.append([victim])
        owner.append(i)
    assert w==Q(28,135)


def birth_death():
    for z in (Q(1),Q(1,2),Q(1,64),Q(1,1024)):
        p0=(1-z)/(2-z)
        p=z/(2-z)
        active=1-p0
        # Future eventual extinction conditional on current stock k is 2^-k.
        assert active*p/(1+p)==Q(1,2)-p0==z/(2*(2-z))
        for k in (2,3,8,100):
            assert active*(1-(1-p)**(k-1))<=k*z


def beta_mesh():
    for n in range(2,12):
        for m in range(2,12):
            for x in (Q(0),Q(1,5),Q(1,2),Q(1)):
                a=comb(n+m-2,n-1)*x**(n-1)*(1-x)**(m-1)
                assert 0<=a<=1
                assert sum(comb(n+m-2,j)*x**j*(1-x)**(n+m-2-j)
                           for j in range(n+m-1))==1
    # Main mesh term: 2/16 -(1-1/4)/2 is -1/4 in log(1/|A|).
    assert Q(2,16)-(1-Q(1,4))/2==-Q(1,4)


if __name__=='__main__':
    cuts()
    birth_death()
    beta_mesh()
    print('PASS: two-slot middle, eight completed jump factors, unfinished waits')
    print('PASS: defective extinction formula, geometric stock bounds, 400 Beta caps, mesh exponents')
