#!/usr/bin/env python3
"""Independent exact fixtures for the finite-block right-tail argument."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import comb, factorial


def degrees(edges):
    d = defaultdict(int)
    for a, b in edges:
        d[a] += 1
        d[b] += 1
    return d


def setup(scaffold, anchors):
    k = len(anchors)
    n0 = 3*k
    old = degrees(scaffold)
    assert len(old) == n0 and max(old.values()) <= 3
    donors = {a: min(v if u == a else u for u,v in scaffold if a in (u,v))
              for a in anchors}
    mass = Q(1)
    edges = list(scaffold)
    # For a used set, label order is irrelevant to rates. Check all subsets
    # including adjacent anchors that are each other's designated source.
    for stage in (0, 1):
        base = list(edges)
        for j in range(k):
            for used in combinations(range(k), j):
                ss = base + [(anchors[i], n0+i) if not stage else (n0+i,n0+k+i)
                             for i in used]
                d = degrees(ss)
                q = sum(Q(1, len(d)*d[anchors[i] if stage else donors[anchors[i]]])
                        for i in range(k) if i not in used)
                assert q >= Q(k-j,(20 if stage else 16)*k)
            mass *= Q(k-j,(20 if stage else 16)*k)
        edges += [(a,n0+i) if not stage else (n0+i,n0+k+i)
                  for i,a in enumerate(anchors)]
    assert mass == Q(factorial(k)**2,(320*k*k)**k)
    return edges


def cut(edges,colors):
    d=degrees(edges);local=degrees([e for e in edges if colors[e[0]]==colors[e[1]]])
    tv=Q(0);cross=0
    for a,b in edges:
        if colors[a]==colors[b]:
            tv+=Q(1,local[a])-Q(1,d[a])+Q(1,local[b])-Q(1,d[b])
        else:
            cross+=1;tv+=Q(1,d[a])+Q(1,d[b])
    assert tv==2*sum(Q(d[x]-local[x],d[x]) for x in d)<=4*cross
    assert sum(x>=2 for x in d.values())-sum(x>=2 for x in local.values())<=2*cross
    return d,local


def block(n):
    """Rare comb: original root stays a leaf and supplies each fresh tooth."""
    e=[(0,1)];prob=Q(1);visits=0
    for j in range(2,n):
        source,target=(0,1) if j%2==0 else (1,j-1)
        d=degrees(e)
        assert any({u,v}=={source,target} for u,v in e)
        prob*=Q(1,j*d[source]);visits+=source==0;e.append((target,j))
    assert visits==(n-1)//2 and degrees(e)[0]==1
    assert sum(x>=2 for x in degrees(e).values())>=visits
    return e,prob,visits


def schedules(left):
    if not any(left):
        yield ()
    for i,v in enumerate(left):
        if v:
            for s in schedules(left[:i]+(v-1,)+left[i+1:]):yield (i,)+s


def labelled_glue():
    sc=[(i,i+1) for i in range(8)]
    edges=setup(sc,(0,1,2));k=3;n0=9;m=5
    colors={i: -1 if i<n0 else (i-n0)%k for i in range(5*k)}
    cut(edges,colors)
    # Multiple cut edges at one scaffold vertex also obey the TV bound.
    multi=[(u if colors[u]>=0 else 0,v) if colors[u]!=colors[v] else (u,v)
           for u,v in edges]
    cut(multi,colors)
    refs=actual=Q(0);seen=set()
    for s in schedules((m-2,)*k):
        e=edges[:];cs=colors.copy();w=v=Q(1);used=[0]*k
        label=[[n0+i,n0+k+i] for i in range(k)];order=[]
        for j in s:
            d,f=cut(e,cs);size=len(d);loc=used[j]+2
            i,t=(0,1) if loc%2==0 else (1,loc-1)
            a,b=label[j][i],label[j][t]
            w*=Q(1,size*f[a]);v*=Q(1,size*d[a]);order.append((a,b))
            assert {a,b} in ({x,y} for x,y in e)
            used[j]+=1;label[j].append(size);cs[size]=j;e.append((b,size))
        assert tuple(order) not in seen
        seen.add(tuple(order));refs+=w;actual+=v
    avoidance=Q(1)
    for j in range(k*(m-2)):avoidance*=Q(2*k+j,5*k+j)
    balance=Q((m-1)**k,comb(k*m-1,2*k-1))
    _,h,r=block(m)
    assert refs==avoidance*balance*h**k and actual==refs/2**(k*r)
    print('PASS',len(seen),'distinct interleaved physical histories, adjacent anchors, and shared-anchor cut budget')


def biased_trim():
    n=80
    _,h,r=block(n)
    # No conditional Bernoulli argument works on this event: it forces r
    # root visits and a linear core count. Compare the exact unconditional
    # Bernoulli tail with the single history's mass.
    row=[1]
    for j in range(2,n):
        nxt=[0]*(len(row)+1)
        for i,w in enumerate(row):nxt[i]+=w*(j-1);nxt[i+1]+=w
        row=nxt
    den=factorial(n-1)
    assert sum(row)==den and row[0]==den//(n-1)
    assert sum(i*v for i,v in enumerate(row))*Q(1,den)==sum(Q(1,j) for j in range(2,n))
    assert Q(sum(row[r:]),den)>h and r==39
    # A moderate rare event can still force a substantial bias: root on
    # all the first 15 source draws, while afterward retaining every draw.
    forced=Q(1,factorial(16))
    assert forced>h and forced<Q(1,10**12)
    print('PASS biased-source stress: linear core comb forces 39 root visits; unconditional trim cannot be used at arbitrary event cost')


def forecasts():
    # Exact optimizer/scale controls at a=1/3, mu=2/3: use t=z**3.
    a,mu,p=Q(1,3),Q(2,3),3
    assert 2+p*(2*mu-1)==p
    B=Q(5,4);z=Q(1,3)/B
    assert all(z*(1-B*z)**2>=u*(1-B*u)**2 for u in (Q(i,1000)/B for i in range(1001)))
    m=900;g=Q(3,2);b=Q(1,4)
    assert ((m+3)*((b*m-g)/(m+3))**p)/(m*b**p)==(1-g/(m*b))**p*(1+Q(3,m))**(1-p)
    print('PASS repetition overhead, tail scaling and Hoeffding optimizer algebra')


if __name__=='__main__':
    labelled_glue()
    biased_trim()
    forecasts()
