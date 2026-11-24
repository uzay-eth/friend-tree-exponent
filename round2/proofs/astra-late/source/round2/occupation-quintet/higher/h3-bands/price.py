#!/usr/bin/env python3
"""Adverse high-port H3 prices with elastic reward bands and exact degree clocks."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'full-h3'))
from cohort import GRID, HI, LO, J, I, Reader, states
from math import lcm


def bands(first, last, step):
    """Geometric degree panels after the exact G2 seam."""
    while first <= last:
        end = min(last, first + max(1, first // step) - 1)
        yield first, end
        first = end + 1


def cubic_bound(p, lo, hi):
    """Range enclosure of cubic p(x)/(1+x) on a nonnegative interval."""
    a,b,c,d=p
    def value(x):
        return a+x*(b+x*(c+x*d))
    def slope(x):
        return b+x*(2*c+3*x*d)
    w=(hi-lo)/3
    v,u=value(lo),value(hi)
    z=(v/(1+lo),(v+w*slope(lo))/(1+lo+w),
       (u-w*slope(hi))/(1+hi-w),u/(1+hi))
    return min(z),max(z)


class Replay:
    def __init__(self, r, a, upper):
        self.r, self.a, self.upper = r, a, upper
        self.parents = (0, 0) if a.orientation.startswith('crown') else (0, 1)
        self.order = (1, 0, 2) if self.parents == (0, 0) else (0, 1, 2)
        self.high = self.order[0 if a.orientation in ('mother', 'crown-old') else 2]
        self.mid = self.order[1]
        self.other = next(i for i in range(3) if i not in (self.high, self.mid))
        self.caps = tuple(a.last if i == self.high else a.mid if i == self.mid else a.other for i in range(3))
        stock = (a.stock if i == self.high else a.mid_stock if i == self.mid else a.other_stock for i in range(3))
        self.axes = [tuple(reversed(states(d, c, i == self.other))) for i, (d, c) in enumerate(zip(self.caps, stock))]
        self.intervals = {d: (lo, hi) for lo, hi in bands(max(180, a.first), a.last, a.step) for d in range(lo, hi+1)}
        self.discount = LO if upper else HI
        self.unit = lcm(self.discount.denominator, *range(2, max(self.caps+(a.other+1,))+1))
        self.inv = [0]+[self.unit//d for d in range(1,max(self.caps+(a.other+1,))+1)]
        self.disc = self.discount.numerator*(self.unit//self.discount.denominator)
        self.round = (lambda x: -(-GRID*x//1)) if upper else (lambda x: GRID*x//1)
        self.reward = lru_cache(maxsize=32000)(self.reward)
        self.gain = lru_cache(maxsize=50000)(self.gain)
        self.polynomial = lru_cache(maxsize=10000)(self.polynomial)

    def polynomial(self,low,middle,c,e,f,g):
        # The row times (1+x) is cubic: all affine payload jumps and the
        # finite-neighbour gain terms contribute degree <=3, and the tail
        # gain's reciprocal denominator is exactly (1+x).
        t=180*3
        v=[]
        for j in range(4):
            x=Q(j,t)
            row=self.r.tail_row(low,middle,c,I(x),tuple(I(z) for z in (e,f,g)))
            assert row.a==row.b
            v.append((1+x)*row.a)
        u=v[1]-v[0]; w=v[2]-2*v[1]+v[0]; z=v[3]-3*v[2]+3*v[1]-v[0]
        return v[0], t*(u-w/2+z/3), t*t*(w-z)/2,t**3*z/6

    def reward(self, low, middle, c, left, right):
        """Bound the same clipped reward on each panel, on both witnesses.

        For a fixed degree, freezing the high rate means maximizing the row in
        that coordinate. Other exterior rates are live in the target; use
        their adverse/favorable extrema for the floor/ceiling respectively.
        """
        lim = (Q(low[1]-1,2), Q(middle[1]-2,2), Q(c-1,2))
        if left == right and left < 180:
            v, b = self.r.h3_physical((low, middle, (left,c)))
            fn = min if self.upper else max
            row = v + max(0,b[2]*lim[2]) + sum(fn(0,z*w) for z,w in zip(b[:2],lim[:2]))
        else:
            x = Q(1,right), Q(1,left)
            ext = [(Q(0),) if not w else (Q(0), w) for w in lim]
            reduced = []
            for e,f in product(*ext[:2]):
                parts = [cubic_bound(self.polynomial(low,middle,c,e,f,g),*x) for g in ext[2]]
                reduced.append(max(z[not self.upper] for z in parts))
            row = (min if self.upper else max)(reduced)
        return self.round(max(0,-row))

    def gain(self, mark, adjacent):
        d,c = mark
        lo = (d-c)*self.unit + sum(self.inv[e] for e in adjacent)
        hi = lo + (c-len(adjacent))*self.inv[2] + adjacent.count(0)*self.inv[self.a.other+1]
        return lo,hi

    def value(self, reward, terms):
        a,b=reward*self.unit,self.disc
        flex=[]
        for lo,hi,v in terms:
            a+=lo*v; b+=lo
            if hi>lo:
                flex.append((v,hi-lo))
        x,y=a,b
        for v,w in sorted(flex,reverse=self.upper):
            if (v*b<=a if self.upper else v*b>=a):
                break
            a+=w*v; b+=w
        z=-(-a//b) if self.upper else a//b
        err=x-y*z+sum(w*(max if self.upper else min)(0,v-z) for v,w in flex)
        assert z>=0 and (err<=0 if self.upper else err>=0)
        return z

    def earning(self,p):
        h,c=p[self.high]; m=p[self.mid]; low=p[self.other]
        if not self.a.first<=h<=self.a.last or not low[0]:
            return None
        if self.a.orientation == 'mother' and h<=65 and m[0]<=5 and low[0]<=65 and (c,m[1],low[1])==(2,2,1):
            return None
        if self.a.orientation == 'reverse' and h<=65 and m[0]<=5 and (low[1],m[1],c)==(2,2,1):
            return None
        return low,m,c,*self.intervals.get(h,(h,h))

    def solve(self):
        nxt={}
        for k in (3,2,1):
            adj=tuple(tuple(j for j in range(k) if (j>0 and self.parents[j-1]==i) or (i>0 and self.parents[i-1]==j)) for i in range(k))
            axes=tuple(tuple(s for s in self.axes[i] if s[1]>=len(adj[i])+(i==0)) for i in range(k))
            # Use only two successive exact high-degree layers. The next
            # shorter prefix needs exclusively entries whose last mark is J.
            layers=({},{}); saved={}
            order=([self.high] if self.high<k else [])+[i for i in range(k) if i!=self.high]
            prev_d=None
            count=pay=pos=0
            for q in product(*(axes[i] for i in order)):
                p=tuple(q[order.index(i)] for i in range(k))
                d0=p[self.high][0] if self.high<k else None
                if d0!=prev_d:
                    layers=({},layers[0]); prev_d=d0
                cur,next_degree=layers
                key=self.earning(p) if k==3 else None
                reward=self.reward(*key) if key else 0
                pay+=key is not None
                pos+=reward>0
                # A retained nonzero potential is looked up only after a
                # physical monotone gain or promotion (or the newborn event).
                terms=[]
                for i,(d,c) in enumerate(p):
                    if d:
                        target=(0,c) if i==self.other and d==self.caps[i] else (d+1,c)
                        g=p[:i]+(target,)+p[i+1:]
                        lo,hi=self.gain((d,c),tuple(sorted(p[j][0] for j in adj[i])))
                        terms.append((lo,hi,(next_degree if i==self.high else cur).get(g,0)))
                        lo=hi=(d-c)*self.inv[d]
                    else:
                        lo,hi=max(0,self.a.other+1-c)*self.inv[self.a.other+1],self.unit
                    b=p[:i]+((d,c+1),)+p[i+1:]
                    value=cur.get(b,0)
                    if k<3 and i==self.parents[k-1]:
                        value+=nxt.get(b+(J,),0)
                    terms.append((lo,hi,value))
                value=self.value(reward,terms)
                if value:
                    cur[p]=value
                    if p[-1]==J:
                        saved[p]=value
                count+=1
            print(self.r.side,'ceiling' if self.upper else 'floor','prefix',k,'states',count,'pay/positive',pay,pos,'birth slice',len(saved),flush=True)
            nxt=saved
        return (HI if self.upper else LO)*Q(nxt.get((J,),0),GRID*self.r.s.h['S'])


def check(r):
    for m,l,c in product(((2,2),(4,3)),((2,1),(8,3)),(3,7,15)):
        for lo,hi in ((180,190),(330,499),(1800,1998)):
            fake=argparse.Namespace(first=lo,last=hi,step=16,orientation='mother',mid=5,other=8,stock=15,mid_stock=4,other_stock=3)
            floor=Replay(r,fake,False).reward(l,m,c,lo,hi)
            ceil=Replay(r,fake,True).reward(l,m,c,lo,hi)
            for d in (lo,(lo+hi)//2,hi):
                a,b=r.h3_physical((l,m,(d,c)))
                a+=max(0,b[2]*Q(c-1,2))
                for x,y in product((0,Q(l[1]-1,2)),(0,Q(m[1]-2,2))):
                    assert floor<=GRID*max(0,-a-b[0]*x-b[1]*y)<=ceil
    print(r.side,'432 exact panel/port inequalities',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--orientation',choices=('mother','reverse','crown-old','crown-young'),default='mother')
    p.add_argument('--first',type=int,default=66)
    p.add_argument('--last',type=int,default=1998)
    p.add_argument('--mid',type=int,default=5)
    p.add_argument('--mid-stock',type=int,default=4)
    p.add_argument('--other',type=int,default=8)
    p.add_argument('--other-stock',type=int,default=3)
    p.add_argument('--stock',type=int,default=10)
    p.add_argument('--step',type=int,default=16)
    p.add_argument('--ceiling',action='store_true')
    p.add_argument('--check',action='store_true')
    a=p.parse_args()
    assert 19<=a.first<=a.last<=1998 and 2<=a.mid<=37 and 2<=a.mid_stock<=4
    assert 2<=a.other<=18 and 1<=a.other_stock<=a.other and 2<=a.stock<=27 and a.step>=2
    # The older 38..65 H3 baseline is frozen. The old crown-control masks
    # did not enter the combined paid ledger. All young passes require low<19.
    assert a.orientation!='mother' or a.last<=37 or a.first>=66
    for side in a.side or ('lower','upper'):
        r=Reader(side)
        assert all(r.s.support('h3',((d,c),)) for d,c in states(a.mid,a.mid_stock,False) if c>=2)
        if a.check:
            check(r)
        else:
            value=Replay(r,a,a.ceiling).solve()
            print('total',side,a.orientation,'ceiling' if a.ceiling else 'floor',value,float(value),flush=True)
        r.s.close()
        r.path.cache_clear(); r.affine.cache_clear(); r._h3_physical.cache_clear()


if __name__=='__main__':
    main()
