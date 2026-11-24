#!/usr/bin/env python3
"""Stationary bounds on the additional-child-only H3 recovery."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, permutations
from math import lcm
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'higher/full-h3'))
from cohort import GRID, HI, LO, J, Reader

PARENT={'m':'h','y':'m','s':'h'}
ORDERS=tuple(('h',)+p for p in permutations('mys') if p.index('m')<p.index('y'))


def registration_check():
    from itertools import combinations
    parents=(None,0,1,1,2,1,3,2,1,5,8,5,6,3,2,9)
    found=[]
    for order in ORDERS:
        prefixes=[set() for _ in order]
        for birth in range(1,len(parents)):
            for n in range(3,0,-1):
                idx=order.index(PARENT[order[n]])
                prefixes[n].update(p+(birth,) for p in prefixes[n-1] if p[idx]==parents[birth])
            prefixes[0].add((birth,))
        found.extend(tuple(p[order.index(i)] for i in 'hmys') for p in prefixes[-1])
    expected=[]
    for group in combinations(range(1,len(parents)),4):
        for h,m,y,s in permutations(group):
            if parents[m]==h and parents[y]==m and parents[s]==h:
                expected.append((h,m,y,s))
    assert set(found)==set(expected) and len(found)==len(set(found)) and expected
    print('registration oracle',len(expected),'unique distinguished incidences')


def units(r,cut,min_degree,ceiling,orders):
    # (cut+1,0) is a permanent reciprocal interval, not a discarded donor.
    over=(cut+1,0)
    D={'h':65,'m':5,'y':18,'s':cut}
    C={'h':7,'m':4,'y':3,'s':cut}
    axes={j:tuple((d,c) for d in range(2,D[j]+1) for c in range(1,min(d,C[j])+1))
          +( (over,) if j=='s' else ()) for j in D}
    discount=LO if ceiling else HI
    den=lcm(discount.denominator,*range(2,66))
    def exact(q):
        assert den%q.denominator==0
        return q.numerator*(den//q.denominator)
    discount=exact(discount)
    inv=[0]+[den//j for j in range(1,66)]

    @lru_cache(maxsize=30000)
    def affine(h,m,y):
        dh,ch=h; dm,cm=m; dy,cy=y
        if dh<min_degree or (ch,cm,cy)==(2,2,1):
            return None
        a,b=r.h3_physical((h,m,y))
        if not b[0]:return None
        n=ch-2
        assert n>0
        other=(1,cm-2,cy-1)
        direction=min if ceiling else max
        B=a+n*max(0,b[0]/2)+sum(direction(0,z*k/2) for z,k in zip(b,other))
        return B,b[0],n

    @lru_cache(maxsize=100000)
    def reward(h,m,y,s):
        x=affine(h,m,y)
        if x is None:return 0
        B,b,n=x
        lo,hi=(Q(0),Q(1,cut+1)) if s==over else (Q(1,s[0]),)*2
        z=max(0,b/2)-b*((lo if b>0 else hi) if ceiling else (hi if b>0 else lo))
        assert z>=0
        w=(max(0,n*z-B)-max(0,-B))/n if ceiling else max(0,z-B)-max(0,-B)
        assert w>=0
        return -(-GRID*w//1) if ceiling else GRID*w//1

    result=[]
    terminal=None
    for order in orders:
        nxt={}
        for n in range(4,0,-1):
            roles=order[:n]
            adj=tuple(tuple(k for k,j in enumerate(roles) if PARENT.get(j)==i or PARENT.get(i)==j) for i in roles)
            nxt_parent=roles.index(PARENT[order[n]]) if n<4 else -1
            pay=tuple(roles.index(j) for j in 'hmys') if n==4 else ()
            if n==4 and terminal is not None:
                nxt=terminal
                continue
            eff_axes=[tuple(p for p in reversed(axes[j]) if p==over and j=='s' or p[1]>=len(adj[i])+(i==0)) for i,j in enumerate(roles)]
            cur={};positive=0
            for p in product(*eff_axes):
                val=reward(*(p[i] for i in pay)) if pay else 0
                positive+=val>0
                mass,total=discount,val*den
                flex=[]
                for i,(d,c) in enumerate(p):
                    if roles[i]=='s' and (d,c)==over:
                        # All physical gains and promotions leave this abstract donor unchanged.
                        continue
                    g=over if roles[i]=='s' and d==cut else (d+1,c)
                    visible=[p[j] for j in adj[i]]
                    lo=(d-c)*den+sum(inv[e] for e,k in visible if k)
                    width=(c-len(visible))*inv[2]+sum(inv[e] for e,k in visible if not k)
                    vg=cur.get(p[:i]+(g,)+p[i+1:],0)
                    pb=p[:i]+((d,c+1),)+p[i+1:]
                    birth=(d-c)*inv[d]
                    vb=cur.get(pb,0)
                    if i==nxt_parent:
                        child=pb+(J,)
                        if n==3:child=tuple(child[order.index(j)] for j in 'hmys')
                        vb+=nxt.get(child,0)
                    mass+=lo+birth
                    total+=lo*vg+birth*vb
                    if width:flex.append((vg,width))
                x,y=total,mass
                for v,w in sorted(flex,reverse=ceiling):
                    if (v*y<=x if ceiling else v*y>=x):break
                    x+=v*w;y+=w
                z=-(-x//y) if ceiling else x//y
                lhs=total-mass*z+sum(w*(max if ceiling else min)(0,v-z) for v,w in flex)
                assert z>=0 and (lhs<=0 if ceiling else lhs>=0)
                cur[tuple(p[i] for i in pay) if pay else p]=z
            print(r.side,'ceiling' if ceiling else 'floor',''.join(order),n,len(cur),positive,flush=True)
            nxt=cur
            if n==4:terminal=cur
        v=(HI if ceiling else LO)*Q(nxt[(J,)],GRID*r.s.h['S'])
        result.append(v)
        print('order',r.side,''.join(order),v,float(v),flush=True)
    return sum(result)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=('lower','upper'),required=True)
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--cut',type=int,default=2)
    p.add_argument('--min-degree',type=int,choices=(19,38),default=19)
    p.add_argument('--order',choices=tuple(map(''.join,ORDERS)))
    a=p.parse_args()
    assert 2<=a.cut<=8
    registration_check()
    r=Reader(a.side)
    v=units(r,a.cut,a.min_degree,a.ceiling,tuple(o for o in ORDERS if not a.order or ''.join(o)==a.order))
    print('total',a.side,'ceiling' if a.ceiling else 'floor',a.cut,a.min_degree,v,float(v),flush=True)
    r.s.close()


if __name__=='__main__':main()
