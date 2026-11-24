#!/usr/bin/env python3
"""Exact local checks of Claim 5(d,e); no Fable transport code used."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import prod


def graph(t):
    adj=[[] for _ in t]
    for j,i in enumerate(t[1:],1):
        adj[i].append(j)
        adj[j].append(i)
    return adj


def rates(a):
    return tuple(sum((Q(1,len(a[w])) for w in v),Q(0)) for v in a)


def physical():
    rows=steps=leaves=0
    conditional=set()
    for n in range(2,8):
        for p in product(*(range(i) for i in range(1,n))):
            t=(-1,)+p
            a=graph(t); r=rates(a)
            assert sum(r)==n
            delta=[len(v)-q for v,q in zip(a,r)]
            if n>2:
                assert r[-1]==Q(1,len(a[t[-1]]))<=Q(1,2)
                conditional.add((n,n-r[-1]))
            for dst in range(n):
                b=graph(t+(dst,)); s=rates(b)
                for v in range(n):
                    assert len(b[v])-s[v]>=delta[v]
                    if v==dst:
                        assert s[v]==r[v]+1
                    steps+=1
            for v in range(n):
                l=sum(len(a[w])==1 for w in a[v])
                if len(a[v])>=2:
                    assert r[v]-Q(l,len(a[v]))>=l-1
                    rows+=1
                if len(a[v])==1:
                    u,=a[v]
                    j=sum(w!=v and len(a[w])==1 for w in a[u])
                    if j>=2:
                        # Force the tagged leaf to be promoted at reward one.
                        drift=(j+1)*(Q(1,j+1)-Q(1,j))
                        drift+=Q(j,len(a[u]))*(Q(1,j-1)-Q(1,j))
                        drift+=(1-Q(1,j))/len(a[u])
                        assert drift<=0
                        leaves+=1
    assert (4,Q(7,2)) in conditional and (4,Q(11,3)) in conditional
    print('PASS',steps,'physical deficit transitions;',rows,'leaf-stock drift rows;',leaves,'tagged-leaf barrier rows')


def beta_urn():
    def rising(x,k):
        return prod(x+j for j in range(k))
    count=0
    for m,a in ((2,Q(1)),(3,Q(1,2)),(4,Q(1,3)),(5,Q(7,3)),(7,Q(3,4))):
        b=m-a
        for k in range(9):
            weights=Counter()
            for seq in product((0,1),repeat=k):
                red=a; p=Q(1)
                for j,x in enumerate(seq):
                    p*= (red if x else m+j-red)/(m+j)
                    red+=x
                weights[sum(seq)]+=p
                assert p==rising(a,sum(seq))*rising(b,k-sum(seq))/rising(m,k)
                count+=1
            assert sum(weights.values())==1
            assert weights[k]==rising(a,k)/rising(m,k)
    print('PASS',count,'fractional Pólya paths; all-red path matches actual all-attachment product')


def residual():
    # A third vertex at the A root makes A a 2-star; B is a path rooted
    # at a leaf. The A leaf matched at odd depth and B's parent stay fixed.
    tA=(-1,0,0); tB=(-1,0,1)
    pa=pb=0
    for k in range(24):
        A,B=graph(tA),graph(tB)
        r,s=rates(A),rates(B)
        common0=s[0]; common1=r[1]
        assert r[0]-common0>=1 and s[1]-common1>=Q(1,2)
        rho=len(A)-(common0+common1)
        assert rho>0 and (r[0]-common0)*(s[1]-common1)/(rho*len(A))>0
        pa=sum(len(A[w])==1 for w in A[0])
        pb=sum(w!=0 and len(B[w])==1 for w in B[1])
        assert pa>=k+2 and pb>=k+1
        assert len(B[0])==1
        tA+=(0,);tB+=(1,)
    print('PASS 24 positive residual forcing steps with a fixed matching and protected B leaf')


if __name__=='__main__':
    physical()
    beta_urn()
    residual()
