#!/usr/bin/env python3
"""Price exact FY gaps on ordered, birth-owned physical two-edge paths."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from replay import Edges, Snapshot, PRIVATE, PINS, CAP, LO, HI, J

GRID=10**6


class Price:
    def __init__(self,e,host,tip,stock):
        self.e=e
        self.host,self.tip,self.stock=host,tip,stock

    @lru_cache(None)
    def extremes(self,h,r):
        e=self.e
        d,c=h
        z=e.g((d+1,c),r)-e.g(h,r)
        vals=[e.sig*(t*z-e.fy(h,t)) for t in e.knots]
        return z,min(vals),max(vals)

    def gap(self,h,r,t,disjoint=True):
        e=self.e
        if h not in e.snap.inner or h in e.mids or h[1]<2:
            return Q(0),Q(0)
        if disjoint and max(h[0],r[0])<=18:
            return Q(0),Q(0)
        z,lo,hi=self.extremes(h,r)
        value=e.sig*(Q(1,t[0])*z-e.fy(h,Q(1,t[0])))
        assert lo<=value<=hi
        return hi-value,hi-lo

    def bellman(self,states,adj,reward,step,new=None):
        # Each tracked vertex has known incoming reciprocals from the path;
        # all other neighbors contribute an interval [0, 1/2] each.
        num=Q(reward)*GRID; den=HI; adjustable=[]
        for i,(d,c) in enumerate(states):
            l=d-c
            r=l+sum(Q(1,states[j][0]) for j in adj[i]); w=Q(c-len(adj[i]),2)
            assert w>=0
            g=step(i,(d+1,c))
            if l:
                b=Q(l,d)
                num+=b*(step(i,(d,c+1))+(0 if new is None else new(i)))
                den+=b
            num+=r*g; den+=r
            if w: adjustable.append((g,w))
        base_num,base_den=num,den
        best=num/den
        # In a linear-fractional minimum, add any high rate whose continuation
        # lies below the current ratio, in increasing continuation order.
        for g,w in sorted(adjustable):
            if g>=best: break
            num+=w*g; den+=w; best=num/den
        v=best.numerator//best.denominator
        assert v>=0
        # Check the exact adverse affine corner at the rounded test value.
        assert base_num-base_den*v+sum(w*min(0,g-v) for g,w in adjustable)>=0
        return v

    def ok(self,s,host=False):
        d,c=s
        return d<=(self.host if host else self.tip) and c<=self.stock

    def path(self,h,a,b,mode):
        if not mode and a>b:
            a,b=b,a
        return self.path_cached(h,a,b,mode)

    @lru_cache(None)
    def path_cached(self,h,a,b,mode):
        if not(self.ok(h,True) and self.ok(a) and self.ok(b)): return 0
        # Mode 0: both ordered directions. Mode 1: first neighbor pays.
        # Mode 2: second neighbor pays. The later-born neighbor is always b.
        reward=sum(self.gap(h,r,t)[0] for i,(r,t) in enumerate(((a,b),(b,a)),1)
                   if mode in (0,i))
        # Nonzero support cannot be reached if c<2, but all physical paths
        # already have two incident neighbors.
        return self.bellman((h,a,b),((1,2),(0,),(0,)),reward,
                            lambda i,x:self.path(*replace((h,a,b),i,x),mode))

    @lru_cache(None)
    def edge(self,h,a,mode):
        if not(self.ok(h,True) and self.ok(a)): return 0
        return self.bellman((h,a),((1,),(0,)),0,
                            lambda i,x:self.edge(*replace((h,a),i,x),mode),
                            lambda i:self.path((h[0],h[1]+1),a,J,mode) if i==0 else 0)

    def root(self,kind,mode):
        @lru_cache(None)
        def v(s):
            children=kind=='children'
            if not self.ok(s,children):return 0
            d,c=s
            return self.bellman((s,),((),),0,lambda i,x:v(x),lambda i:
                self.edge((d,c+1),J,mode) if children else self.edge(J,(d,c+1),mode))
        return LO*Q(v(J),GRID)

    def audit(self):
        # Paying-edge endpoint degrees, not the other donor's degree, decide
        # overlap with the previous mask. The FY arity seam yields both cases.
        overlap=0; outside=0
        for h,r,t in product(((16,16),(18,17),(19,17),(38,3),(42,3)),
                             ((2,1),(8,2),(18,2),(19,2),(180,1),(237,1)),
                             ((2,1),(4,1),(9,2),(180,1),(219,1))):
            gap,width=self.gap(h,r,t,False)
            if max(h[0],r[0])<=18:
                assert self.gap(h,r,t)==(0,0)
                overlap+=gap>0
            else:
                assert self.gap(h,r,t)==(gap,width)
                outside+=gap>0
        assert overlap and outside
        print('PASS positive mask overlap exists; excluded by paying edge:',overlap,'inside,',outside,'outside fixtures')
        paid=0
        for h,r,t in product(((16,16),(18,17),(19,17),(38,2),(45,3)),
                             ((2,1),(8,2),(18,2),(19,2)),((2,1),(8,1),(219,1))):
            dh=[Q(1,t[0])]*(h[1]-1); dr=[Q(1,9)]*(r[1]-1)
            actual=self.e.row(h,r,dh,dr)
            gap=self.gap(h,r,t,False)[0]*(h[1]-1)
            assert 0<=gap<=self.e.cap-actual
            paid+=1
        print('PASS exact charged row plus FY increment under cap:',paid,'physical stars')
        for h,r in (((16,16),(2,1)),((38,2),(2,1)),((38,2),(4,1)),((45,3),(5,2)),((38,2),(237,1))):
            vals=[self.gap(h,r,(d,1),False) for d in (2,3,4,6,9,18,48,180,237)]
            print('conditional host/receiver',h,r,'donor D=2,3,4,6,9,18,48,180,237',
                  'gaps ticks',[round(float(g),3) for g,_ in vals], 'box width ticks',float(vals[0][1]))


def replace(t,i,x): return (*t[:i],x,*t[i+1:])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=48)
    p.add_argument('--tip',type=int,default=4)
    p.add_argument('--stock',type=int,default=5)
    p.add_argument('--side',choices=('lower','upper','both'),default='both')
    p.add_argument('--children-only',action='store_true',help='replay only the two-child contribution')
    p.add_argument('--hash',action='store_true')
    a=p.parse_args()
    assert __debug__ and 38<=a.host<180 and 2<=a.tip<=a.host and 2<=a.stock<=26
    receipt=json.loads((PRIVATE/'gold.json').read_bytes())
    for side in CAP if a.side=='both' else (a.side,):
        assert receipt[side]['verify']['sha256']==PINS[side][2]
        with Snapshot(PRIVATE/f'quintet_{side}.bin',side,verify=a.hash) as s:
            print(side,'caps host/tip/stock',a.host,a.tip,a.stock,flush=True)
            q=Price(Edges(s),a.host,a.tip,a.stock); q.audit()
            total=Q(0)
            for name,kind,mode in (('mother pays','mother',1),('child pays, mother donates','mother',2),
                                   ('two children, both ordered directions','children',0)):
                if a.children_only and kind!='children':continue
                v=q.root(kind,mode); total+=v
                print(name,'ticks/core',v,float(v),'normalized',float(v/s.scale),flush=True)
            print('TOTAL',total,float(total/s.scale),'path classes',q.path_cached.cache_info(),flush=True)
            q.path_cached.cache_clear(); q.edge.cache_clear(); q.extremes.cache_clear()
            s.check_unchanged()
    print('PASS')


if __name__=='__main__':main()
