#!/usr/bin/env python3
"""Compressed stationary H3 joint mother-port recovery."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import lcm
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'higher/full-h3'))
from cohort import GRID, HI, LO, J, Reader, optimum


def young_parent_birth_ceiling(K):
    """All births at parental D<K, including the diagnostic root's actual rate box."""
    v={}
    for d in range(K-1,1,-1):
        for c in range(d,0,-1):
            L=d-c
            vb=v.get((d,c+1),0)
            vg=v.get((d+1,c),0)
            b=Q(L,d)
            v[d,c]=optimum(GRID*b,[(Q(L),Q(L)+Q(c,2),vg),(b,b,vb)],LO,True)
    # v has subgrid units of births, not witness ticks. A founder is omitted
    # by uniform-core stationarity and is not counted in the Palm birth rate.
    return HI*Q(v.get(J,0),GRID)


def run(r,a):
    C=a.arity; Y=a.young; H=a.high; K=a.upstream
    discount=LO if a.ceiling else HI
    den=lcm(discount.denominator,*(range(2,max(H,Y,K)+2)))
    q=discount.numerator*(den//discount.denominator)
    inv=[0]+[den//j for j in range(1,max(H,Y,K)+2)]
    # Count immature external children (D<QD). Their departures across the
    # threshold have rate bounded by (QD-1) per child; QD=3 gives exact D=2.
    QD=3
    threshold=inv[QD]

    @lru_cache(maxsize=25000)
    def ports(h,m,y):
        d,c=h; e,f=m; u,v=y
        if d<a.low or u>=d or (c,f,v)==(2,2,1) and d<=65 and u<=65 and e<=5:
            return (0,)*(1 if a.upstream_only else c-1)
        p,b=r.h3_physical((h,m,y))
        if not b[0]:return (0,)*(1 if a.upstream_only else c-1)
        # pivot = row with adverse entire mother's exterior port. Upstream,
        # middle and young remain physical in the identity; use the adverse
        # maximum of the latter two to minorize the difference of clips.
        other=min if a.ceiling else max
        B=p+max(0,b[0])*Q(c-1,2)+sum(other(0,z)*Q(n,2) for z,n in zip(b[1:],(f-2,v-1)))
        b=b[0]
        B=max(0,B)
        if a.upstream_only:
            z=max(0,b)* (Q(1,2)-Q(1,K))
            return (max(0,(GRID*(z-B))//1),)
        if a.ceiling:
            z0=b*Q(c-1,2) if b>0 else -b*(Q(c-2,QD)+Q(1,K))
            dz=-b*Q(1,2) if b>0 else -b*(Q(1,2)-Q(1,QD))
        else:
            z0=b*(Q(c-1,2)-Q(c-2,QD)-Q(1,K)) if b>0 else Q(0)
            dz=b*(Q(1,QD)-Q(1,2)) if b>0 else -b*Q(1,2)
        return tuple(max(0, -(-GRID*(z0+k*dz-B)//1) if a.ceiling else (GRID*(z0+k*dz-B))//1)
                     for k in range(c-1))

    def register(hi,mu,yu):
        out={}
        final=yu is not None
        has_mid=mu is not None
        # d and c increase. k decreases on a donor's gain, increases on a
        # mother birth (and on neither internal registered child's birth).
        for d in range(H,1,-1):
            for c in range(min(d,C),1 if has_mid else 0,-1):
                for k in range(1 if a.upstream_only else c-int(has_mid)):
                    out[d,c,k]=row(hi,out,(d,c,k),mu,yu,final,has_mid)
        return out

    def row(nxt,cur,h,m,y,final,has_mid):
        d,c,k=h
        n=c-1-int(has_mid)
        # Mother h has n exterior children. Each immature child has D=2;
        # each mature child contributes at most 1/3; the parent at most 1/K.
        l=k*inv[2]; u=l+(n-k)*(inv[2] if a.upstream_only else threshold)+inv[K]
        ext=(l,u-l)
        mass,total=q,0
        flex=[]
        def add(lo,width,v):
            nonlocal mass,total
            mass+=lo;total+=lo*v
            if width:flex.append((v,width))
        if final:
            total=ports((d,c),m,y)[k]*den
        # A visible immature sibling crosses 2→3 with at least 1/D_h.
        # At D=2,c=1 it has one leaf and only h, and at D=2,c=2 it
        # has no leaf and one other neighbor contributing ≤1/2.
        if k:add(k*inv[d],k*den,cur[d,c,k-1])
        gain=cur.get((d+1,c,k),0)
        add((d-c)*den+(inv[m[0]] if has_mid else 0)+ext[0],ext[1],gain)
        born=cur.get((d,c+1,k+int(not a.upstream_only)),0)
        if not has_mid:born+=nxt.get((J,None),{}).get((d,c+1,k),0)
        add((d-c)*inv[d],0,born)
        if has_mid:
            e,f=m
            g=nxt.get(((e+1,f),y),{}).get(h,0) # global already-built next degree
            b=nxt.get(((e,f+1),y),{}).get(h,0)
            if not final:b+=terminal.get((e,f+1),{}).get(h,0)
            add((e-f)*den+inv[d]+(inv[y[0]] if final else 0),(f-1-int(final))*inv[2],g)
            add((e-f)*inv[e],0,b)
        if final:
            e,f=y
            add((e-f)*den+inv[m[0]],(f-1)*inv[2],nxt.get((m,(e+1,f)),{}).get(h,0))
            add((e-f)*inv[e],0,nxt.get((m,(e,f+1)),{}).get(h,0))
        x,denom=total,mass
        for v,w in sorted(flex,reverse=a.ceiling):
            if (v*denom<=x if a.ceiling else v*denom>=x):break
            x+=v*w;denom+=w
        z=-(-x//denom) if a.ceiling else x//denom
        assert isinstance(z,int)
        lhs=total-mass*z+sum(w*(max if a.ceiling else min)(0,v-z) for v,w in flex)
        assert z>=0 and (lhs<=0 if a.ceiling else lhs>=0),(h,m,y,z,total,mass,flex)
        return z

    # Sparse tensor sweep: retain only two adjacent young degree planes,
    # plus the newborn slice that feeds the 2-vertex prefixes.
    tail={}
    mid=tuple((d,c) for d in range(a.middle,1,-1) for c in range(min(d,a.mid_arity),1,-1))
    count=0
    for dy in range(Y,1,-1):
        curr={}
        for y in ((dy,c) for c in range(min(dy,a.young_arity),0,-1)):
            for m in mid:
                # Read current plane + previous degree plane without copying.
                curr[m,y]=register(View(curr,tail),m,y)
                count+=len(curr[m,y])
        tail=curr
    terminal={m:v for (m,y),v in tail.items() if y==J}
    prefix={}
    for dm in range(a.middle,1,-1):
        for cm in range(min(dm,a.mid_arity),0,-1):
            m=(dm,cm)
            prefix[m,None]=register(prefix,m,None)
    birth=Q(register(prefix,None,None)[2,1,0],GRID*r.s.h['S'])
    flow=HI if a.ceiling else max(0,LO-young_parent_birth_ceiling(K))
    print(r.side,'mask',a,'states',count,'birth flux',flow,float(flow),'Palm lifetime',birth,float(birth),
          'ceiling' if a.ceiling else 'floor',flow*birth,float(flow*birth),flush=True)


class View:
    def __init__(self,here,next):self.here,self.next=here,next
    def get(self,k,default):return self.here.get(k,self.next.get(k,default))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=('lower','upper'),required=True)
    p.add_argument('--low',type=int,default=66)
    p.add_argument('--high',type=int,default=100)
    p.add_argument('--arity',type=int,default=8)
    p.add_argument('--middle',type=int,default=5)
    p.add_argument('--mid-arity',type=int,default=4)
    p.add_argument('--young',type=int,default=18)
    p.add_argument('--young-arity',type=int,default=3)
    p.add_argument('--upstream',type=int,default=3)
    p.add_argument('--ceiling',action='store_true',help='run with upstream floor D=2 to cover all births')
    p.add_argument('--upstream-only',action='store_true',help='price second clip, keeping all siblings adverse')
    a=p.parse_args()
    assert 19<=a.low<=a.high<=1998 and 2<=a.upstream<=a.high
    assert 2<=a.middle<=5 and 2<=a.mid_arity<=4 and a.arity>=2 and a.young>=2 and a.young_arity>=1
    assert not a.ceiling or a.upstream==2 and not a.upstream_only
    r=Reader(a.side);run(r,a);r.s.close()
