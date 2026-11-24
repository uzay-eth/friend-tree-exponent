#!/usr/bin/env python3
"""Physical labelled clock, all-degree coupon, and unlimited sign controls."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import comb

from current_gauge_sol import J, degree, moves


S0 = Q(565658,10**6)
KAPPA = Q(8600511260,10**10)


def labels(edges):
    adj = defaultdict(set, {0:{-1}})  # inert selection
    for i,j in edges:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def append(state,i,j):
    A,B=state
    w=1+max(max(edges,default=(0,0))[0] for edges in state)
    if i>=0:
        A=A+((w,i),)
        w+=i!=j
    if j>=0:
        B=B+((w,j),)
    return A,B


def clocks(state):
    A,B = map(labels,state)
    for v in sorted(A.keys()|B.keys()):
        x,y = A.get(v,set()),B.get(v,set())
        pa={i:Q(1,len(x)) for i in x}
        pb={i:Q(1,len(y)) for i in y}
        for i in sorted(x&y):
            r=min(pa[i],pb[i]); pa[i]-=r; pb[i]-=r
            yield v,i,i,r
        if not y:
            yield from ((v,i,-1,r) for i,r in pa.items() if r)
        elif not x:
            yield from ((v,-1,j,r) for j,r in pb.items() if r)
        else:
            rem=sum(pa.values())
            assert rem==sum(pb.values())
            if rem:
                yield from ((v,i,j,r*t/rem) for i,r in pa.items() for j,t in pb.items() if r*t)


def census(state):
    A,B = map(labels,state)
    def descendants_match(v,parent):
        return (A[v]==B[v] and all(descendants_match(w,v) for w in A[v] if w!=parent))
    out=[0]*6  # u,v,c,e,l,k
    for v in A[0]-{-1}:
        if v not in B[0]:
            ix=4+(len(A[v])>1)
        elif len(A[v])==1:
            ix=int(len(B[v])>1)
        else:
            assert len(B[v])>1
            ix=2+int(not descendants_match(v,0))
        out[ix]+=1
    return tuple(out)


def check_state(state):
    a,b = map(labels,state)
    before=census(state)
    u,v,c,e,l,k=before
    A,B=len(a[0]),len(b[0]); m=v+e+l+k
    row=defaultdict(Q); ma,mb=defaultdict(Q),defaultdict(Q)
    for src,i,j,r in clocks(state):
        if i>=0 or src==0 and i in a[0]: ma[src,i]+=r
        if j>=0 or src==0 and j in b[0]: mb[src,j]+=r
        out=append(state,i,j)
        d=census(out)
        for z in set(labels(out[0])) & set(labels(out[1])):
            adj0,adj1 = map(labels,out)
            par=dict(out[0]); age=0; zz=z
            while zz: zz=par[zz]; age+=1
            assert (adj1[z]<=adj0[z]) if age%2==0 else (adj0[z]<=adj1[z])
        if d!=before:
            row[src==0,d]+=r
    for adj,marg in ((a,ma),(b,mb)):
        for src, targets in adj.items():
            for dst in targets:
                if dst>=0:
                    assert marg[src,dst]==Q(1,len(targets))
    dm=sum(r*(sum(z[i] for i in (1,3,4,5))-m) for (_,z),r in row.items())
    assert dm <= m*(1+Q(1,B))
    # Compare non-root transitions at the actual rates, not fitted endpoints.
    common=private=Q(0)
    for i in a[0]-{-1}:
        a_i=Q(1,len(a[i])); b_i=Q(1,len(b[i])) if i in b[0] else Q(0)
        common+=b_i; private+=a_i-b_i
    targets={tuple(u+1 if i==0 else x for i,x in enumerate(before)):common,
             tuple(l+1 if i==4 else x for i,x in enumerate(before)):private}
    assert {z:r for (root,z),r in row.items() if not root}=={z:r for z,r in targets.items() if r}
    return row


def pair_fixture():
    # Two clean child branches (one deep), one common leaf, one private A leaf.
    B=((1,0),(2,0),(3,0),(4,2),(5,4),(6,3),(7,3))
    seed=(B+((8,0),),B)
    first=append(seed,8,2)
    assert (0,8,2,Q(1,20)) in set(clocks(seed))
    assert census(seed)==(1,0,2,0,1,0)
    assert census(first)==(1,0,1,1,0,1)
    assert (8,0,-1,Q(1,2)) in set(clocks(first))
    middle=append(first,0,-1)
    a,b=map(labels,middle)
    new=max(a)
    assert (0,new,1,Q(1,24)) in set(clocks(middle))
    last=append(middle,new,1)
    assert census(last)==(0,1,1,1,0,2)
    for state in (seed,first,middle,last):
        row=check_state(state)
        a,b=map(labels,state)
        before=census(state)
        arcs=defaultdict(Q)
        for i in sorted(a[0]-b[0]):
            for j in sorted(b[0]):
                z=census(append(state,i,j))
                if z!=before: arcs[z]+=Q(1,len(a[0])*len(b[0]))
        # Subtract the synchronous choices (LL and LC promotion only).
        yes={z:r for (root,z),r in row.items() if root}
        for j in sorted((a[0]&b[0])-{-1}):
            z=census(append(state,j,j))
            if z!=before: yes[z]-=Q(1,len(a[0]))
        assert dict(arcs)=={z:r for z,r in yes.items() if r}
    row=check_state(seed)
    joint=census(first)
    assert row[True,joint]==Q(1,10)  # either clean B target
    assert joint not in {census(append(seed,8,-1)),census(append(seed,-1,2))}
    assert (False,(2,0,1,1,0,1)) in check_state(first)  # paired returns
    assert check_state(first)[False,(1,0,1,1,1,1)]>0  # private core plus dirty return
    print('actual successive residual pair censuses:',*(census(x) for x in (seed,first,middle,last)))
    print('rejected split-clock substitution and dropping the private-core return')


def ceil(a,b):
    return -(-a//b)


def coupons_and_integer_limits():
    # New envelope: finer grid, exact Taylor polynomial (no intermediate
    # tick roundings), rational Padé-free discount bounds.
    d=48; N=256; step=Q(1,N); x=(1+Q(1,d))*step
    term=g=Q(1)
    for j in range(1,16):
        term*=x/j; g+=term
    g+=2*term*x/16
    tick=10**15
    gx=ceil(g*tick,1); disc=ceil(tick/(1+S0*step),1)
    age=tick; v=ceil(2*tick,d); ans=0; loops=0
    while v<2*tick:
        v=ceil(v*gx,tick)
        ans+=age*min(2*tick,v)
        age=ceil(age*disc,tick)
        loops+=1
    bound=Q(ans,tick*tick*N)+2*Q(age,tick)/S0
    assert bound<Q('.846')<Q('.8484503829')<KAPPA
    # True large-d integral is decreasing: both initial magnitude and
    # exponential rate decrease.  Checking d=48 covers all d>=48.
    assert d==48 and loops<1024
    max_cap=150; dm=2*max_cap**2
    tvticks=ceil(2*10**10,S0)+1
    # Clock total is <=2 DA and rate scale <=2 DA^2. At the production cap
    # the coupon runs fewer than 512 steps (even d=151 and m=1).
    assert (1+Q(1,64))**512>151
    rates=2*max_cap*dm
    rhs=rates*tvticks+2*dm*10**10
    assert rhs*10**6<2**127 and comb(299,6)<2**64 and tvticks<2**64
    assert 512*2*(10**10)**2<2**127
    print('independent all-degree coupon via d=48:',float(bound),'substeps',loops)


def global_sign(s):
    # Global signed derivative envelope, then true (not killed) low-degree
    # values, paying both continuation signs at the first physical exit.
    Dcap=7
    shape_size=lambda S:1+S[0]+sum(shape_size(Y) for Y in S[1])
    H=1/(s*(1-KAPPA))
    for d in (2,3,Dcap,48,500,10**8):
        for j in (0,(d-1)//2,d-1):
            F=(j+1/s)/(1+s*d)
            assert s*F-Q(d-1-j,d)/(1+s*d)-Q(j,d)==1/(d*(1+s*d))
            # Negative control: charging the inert slot as a leaf loses
            # exactly the advertised strict margin.
            assert s*F-Q(d-j,d)/(1+s*d)-Q(j,d)==0
    @cache
    def visit(S,plant):
        d=degree(S); j=len(S[1]); p=Q(1,d)
        if shape_size(S)>=Dcap:
            if plant: return Q(0),Q(0)
            # All-degree barriers: fitted cutoff at d<=6 is irrelevant.
            C=min(1/s,(j+1/s)/(1+s*d))
            b0=max(Q(0),1/s-(j+1+1/s)/(1+s*d))
            return b0-H*C,1/s+H*C
        row=moves(S); den=s+sum(r for r,_ in row)
        charge=S[0]*p
        lo,hi=charge-p*j*H,charge+p*j*H
        for r,Y in row:
            a,b=visit(Y,plant);lo+=r*a;hi+=r*b
        return lo/den,hi/den
    lo,hi=visit(J,False)
    assert lo<=hi
    if s==2:
        assert hi<1
    else:
        assert lo<1<hi
        a,b=visit(J,True)
        assert lo<a and b<hi  # the killed boundary nontrivially narrows the interval
    print('actual global low-degree fixture s=',s,'interval',float(lo),float(hi))


def critical_constants():
    u=Q('.565798'); lo=Q(565662790716854568223,1000008252270910000000)
    hi=Q(1698046074520368376961,3001152517940340000000)
    assert S0 < lo <= hi < u
    assert Q('.5656581227') < lo and hi < Q('.5657979941')
    assert S0<Q('.56572345')<Q('.56572385')<u
    wide=Q('.565615')
    assert Q('.5649')<wide<S0<lo<=hi<u<Q('.565809')
    assert KAPPA+2*(1/wide-1/S0)<1  # norm continuity down to original rung-3 target
    assert Q(50,7)**4<Q('1.006')**1395
    P,M,m,delta = Q(200,21),Q(50,7),Q('.117'),Q('.021')
    v=1+Q('.006')/delta
    A=P*(1+M/(m*delta)); B=P*v*M/(m*delta)
    factor=2*u*(u+2)**2*(A+B)
    assert 237358<factor<471522  # check the displayed premises, not the printed conclusion
    assert u*(u+2)*P*M/m<845
    assert 1/(1-KAPPA)<Q('7.146')
    print('critical constants: prior and slope pass; displayed full-window factor is',float(factor))

    # Pole-clearing sensitivity control with a rational positive return
    # F=1/(s+1/2), positive resident A and nonzero immigration residue B.
    for s,expected in ((Q(49,100),1),(Q(1,2),0),(Q(51,100),-1)):
        F=1/(s+Q(1,2)); A=1/s; B=Q(2)
        D=(1-F)*A+B
        U=1-(1-F)/(s*D)
        assert s*D*(1-U)==1-F
        assert (U>1)-(U<1)==expected
        if s==Q(1,2):
            assert D==B>0 and (1-F)*A==0  # loss of B destroys the sign bridge


if __name__ == '__main__':
    pair_fixture()
    coupons_and_integer_limits()
    global_sign(Q('0.565728'))
    global_sign(Q(2))
    critical_constants()
