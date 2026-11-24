"""Rational physical bounds for early deterministic-label hub differences."""

from fractions import Fraction as Q
from functools import cache
import argparse


@cache
def silent_bound(j, cap, immigration=1, base=1):
    """Lower probability that a specified neighbour never gets the arrow.

    j tracks OTHER leaves at its parent, with the tag still initiating.
    Compare with gain j+immigration, loss j/(j+1), and killing 1/(j+1).
    """
    ratios = []
    for i in range(cap):
        den=max(2,i+base)
        gain, loss, kill = i+immigration, Q(i,den), Q(1,den)
        ratios.append(gain/(gain+loss+kill-(loss*ratios[-1] if i else 0)))
    v = max(Q(0), 1-Q(2,cap+1))
    for i in range(cap-1, j-1, -1):
        v *= ratios[i]
    return max(v, 1-Q(2,j+1))


@cache
def core_bound(d, l, cap):
    """Track an autonomous lower stock until extinction or degree cap."""
    if not l:
        return Q(0)
    if d >= cap:
        return 1-Q(1,d)**l
    return (d*core_bound(d+1,l+1,cap)+core_bound(d,l-1,cap))/(d+1)


def check_barriers():
    def f(j):
        return min(Q(1),Q(2,j+1))
    for j in range(1,65):
        # No immigration and minimal base are worst for this decreasing
        # probability barrier; a gain of the tag sets the score to one.
        drift=j*(f(j+1)-f(j))+Q(j,j+1)*(f(j-1)-f(j))+(1-f(j))/(j+1)
        assert drift==(Q(-1,3) if j==1 else -Q(j-2,(j+1)*(j+2)))<=0
    for a in (0,1):
        for base in range(1,7):
            values=[silent_bound(j,32,a,base) for j in range(33)]
            assert values==sorted(values)


@cache
def graph(par):
    a = [[] for _ in range(len(par)+1)]
    for v, u in enumerate(par,1):
        a[u].append(v)
        a[v].append(u)
    d = tuple(map(len,a))
    r = tuple(sum((Q(1,d[w]) for w in ns),Q(0)) for ns in a)
    return a,d,r


@cache
def abs_hub(par,u):
    a,d,_ = graph(par)
    if d[u]==1:
        (p,) = a[u]
        j = sum(d[w]==1 for w in a[p])-1
        # Once promoted the tag starts at degree two with one leaf.
        base=d[p]-j
        survival = silent_bound(j,max(j,32),1,base)
        return Q(0),1-survival
    l = sum(d[w]==1 for w in a[u])
    lower = core_bound(d[u],l,max(d[u],12))
    # Becoming a hub requires a NEXT gain; each candidate send is along
    # one of the existing edges. A leaf neighbour can always be charged 1.
    if l:
        return lower,Q(1)
    upper=Q(0)
    for p in a[u]:
        j=sum(d[w]==1 for w in a[p]); base=d[p]-j
        upper+=1-silent_bound(j,32,0,base)
    return lower,min(Q(1),upper)


def with_weight(lo_hi,c):
    lo,hi=lo_hi
    return (c*lo,c*hi) if c>=0 else (c*hi,c*lo)


@cache
def difference(par,u,v,h):
    """Pull back one signed difference through h exact physical jumps."""
    a,d,r=graph(par)
    if sorted(u if x==v else v if x==u else x for x in a[u])==sorted(a[v]) and (
        # Transposition preserves the whole graph iff the two incidence
        # lists agree after transposition; other incidences then agree too.
        len(a[u])==len(a[v])
    ):
        return Q(0),Q(0)
    lu,uu=abs_hub(par,u)
    lv,uv=abs_hub(par,v)
    lo,hi=lu-uv,uu-lv
    if h:
        pull_l=pull_u=Q(0)
        for x,rate in enumerate(r):
            a1,b1=difference(par+(x,),u,v,h-1)
            pull_l+=rate*a1/len(r)
            pull_u+=rate*b1/len(r)
        lo,hi=max(lo,pull_l),min(hi,pull_u)
    assert lo<=hi
    return lo,hi


def adjacent(par,h):
    a,d,r=graph(par)
    n=len(r)
    lo=hi=Q(0)
    w=[r[i]/(d[i]*(d[i]+1)) for i in range(n)]
    for i in range(n):
        for j,c in [(n,r[i]/(d[i]+1))]+[(j,w[j]-w[i]) for j in a[i] if i<j]:
            l,u=with_weight(difference(par+(i,j),n,n+1,h),c)
            lo+=l;hi+=u
    return lo/(n*(n+1)),hi/(n*(n+1))


def out(x, upper=False, places=12):
    unit=10**places
    y=-((-x*unit)//1) if upper else (x*unit)//1
    sign='-' if y<0 else ''
    return f'{sign}{abs(y)//unit}.{abs(y)%unit:0{places}}'


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--horizon',type=int,default=4)
    h=p.parse_args().horizon
    check_barriers()
    # A stand-alone small rational certificate, no finite-size fitting.
    freeze=silent_bound(1,4)
    gap=(freeze-Q(1,4))/2
    assert gap==Q(127,2440)>0
    assert core_bound(2,2,4)==Q(437,512)
    print(f'p1=p2; p2-p3 >= {gap} > {out(gap,places=9)}')
    l,u=difference((0,0),0,1,h)
    print(f'p2-p3, h={h}: [{out(l/2)}, {out(u/2,True)}]')
    # The unlabelled T3 shape is deterministic; p3/p4 is the adjacent
    # two-birth calculation from T2, p4/p5 is from T3.
    for par in ((0,), (0,0)):
        lo,hi=adjacent(par,h)
        n=len(par)+1
        print(f'p{n+1}-p{n+2}, h={h}: [{out(lo)}, {out(hi,True)}]')
        if h>=4:
            assert lo>0
