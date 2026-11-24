"""Price unordered FY pairs with an older donor through a two-vertex test."""

from fractions import Fraction as Q
from pathlib import Path
import argparse
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from replay import Edges, J, LO, HI, PRIVATE, Snapshot


class Exchange:
    def __init__(self, edge, tip, tip_stock):
        self.edge = edge
        self.snap = s = edge.snap
        self.tip, self.tip_stock = tip, tip_stock
        assert 6 <= tip < s.d2 and 1 <= tip_stock <= s.K2
        self.finite = [(d,c) for d in range(2,s.d2) for c in range(1,min(d,s.K2+1)+1)]
        self.labels = np.array([next(iter(edge.weights((d,c),s.inner,s.d2,s.K2)))[0]
                                for d,c in self.finite], dtype=np.int64)
        self.tags = [u for u in self.finite if self.marked(u)]
        self.position = {u:i for i,u in enumerate(self.finite)}
        self.degrees = np.array([d for d,c in self.finite], dtype=np.int64)
        self.high = self.degrees >= 6
        self.both = np.array([self.marked(u) for u in self.finite])
        self.raw = np.ndarray(((s.inf+1)*(s.inf+2)//2,), dtype="<i8", buffer=s.data,
                              offset=512+8*s.off_g)
        self.q = s.header['Q_K']
        self.knots = tuple(int(k*self.q) for k in edge.knots)
        assert edge.knots[0] == 0 and edge.knots[1] == Q(1,s.d2)

    def marked(self,u):
        d,c = u
        return 6 <= d <= self.tip and c <= self.tip_stock

    def rows(self,h):
        """Lower integer ticks for each marked neighbor, for every other neighbor."""
        e, s = self.edge, self.snap
        d,c = h
        if h not in s.inner or h in e.mids or c < 2 or d <= 18:
            return {}
        if (d+1,c) not in s.inner:
            raise ValueError('the host gain must have an exact table row')
        def payload_row(u):
            n = np.arange(s.inf+1,dtype=np.int64)
            j = s.inner[u]
            hi = np.maximum(j,n)
            return self.raw[hi*(hi+1)//2+np.minimum(j,n)]
        old, new = payload_row(h), payload_row((d+1,c))
        # Arithmetic in the 64-bit vector code is justified for the input words.
        assert max(abs(int(z)) for a in (old,new) for z in (a.min(),a.max())) <= 2**44
        slopes = e.sig*(new-old)
        z = slopes[self.labels]
        knot_a = [int(e.sig*e.fy(h,k)) for k in e.knots]
        assert max(map(abs,knot_a)) <= 2**44 and self.q <= 2**14
        m = np.maximum.reduce([(k*z-self.q*a)//self.q for k,a in zip(self.knots,knot_a)])
        by_degree = {d0: e.sig*e.fy(h,Q(1,d0)) for d0 in range(2,s.d2+1)}
        at = np.array([int(by_degree[d0]//1) for d0,_ in self.finite],dtype=np.int64)
        special = slopes[s.tail_base:s.inf]
        zinfty = int(slopes[s.inf])
        # In the tail, each slope is affine between its seam and infinity
        # words. Minimize each supporting line on that segment first. The
        # denominator bound also covers rounds at negative slopes.
        tail_forward = {}
        for dt in range(6,self.tip+1):
            candidates = []
            for k,ak in zip(self.knots,knot_a):
                coef = k*dt-self.q
                small = np.minimum(coef*special,coef*zinfty)
                candidates.append(small//(self.q*dt)-ak)
            tail_forward[dt] = max(0, int(np.maximum.reduce(candidates).min()) + by_degree[dt]//1)
        outputs = {}
        for tag in self.tags:
            dt,_ = tag
            j = self.position[tag]
            # j donates to i for any unmarked i. If i is also a high donor,
            # include its reverse direction. Two marked ends share this sum.
            active = np.maximum(0, m+at[j]+(-z)//dt)
            reverse = np.maximum(0, m[j]+at+(-z[j])//self.degrees)
            counts = np.where(self.both, (active+reverse)//2,
                              active+np.where(self.high,reverse,0))
            finite = int(counts.min())
            if finite == 0:
                outputs[tag] = 0
                continue
            # The reverse term is affine on [0,1/d2], the first knot segment.
            reverse_tail = max(0, int(m[j])+min(knot_a[0],knot_a[1]+int(-z[j])//s.d2))
            outputs[tag] = min(finite, tail_forward[dt]+reverse_tail)
        return outputs


GRID = 10**8


def integer_bellman(reward, rate_den, fixed, variable):
    """Exact upper-discount subsolution, grid times ticks per tagged state."""
    total = sum(rate for rate,_ in fixed)
    gain = reward * GRID * rate_den + sum(rate*value for rate,value in fixed)
    numerator = gain * HI.denominator
    denominator = total * HI.denominator + rate_den * HI.numerator
    # Add unresolved rate only if its destination is below the current value.
    # Sorting the destinations finds the extremizing corner of a fraction.
    a,b = numerator,denominator
    for width,value in sorted(variable,key=lambda item:item[1]):
        if value*b >= a:
            break
        a += width*HI.denominator*value
        b += width*HI.denominator
    estimate = a//b
    drift = numerator-denominator*estimate + HI.denominator*sum(
        width*min(0,value-estimate) for width,value in variable)
    assert estimate >= 0 and drift >= 0
    return estimate


def censuses(degree,stock):
    return [(d,c) for d in range(degree,1,-1) for c in range(min(d,stock),0,-1)]


def direct_check(ex,h,prices):
    """Compare to raw knots on mixed, repeated, overflow and limiting geometries."""
    e = ex.edge
    sample = ((2,1),(5,5),(6,1),(8,2),(15,3),(27,27),(179,48),
              (180,4),(181,30),(10**12+39,2))
    def gap(p,t):
        slope=e.sig*(e.g((h[0]+1,h[1]),p)-e.g(h,p))
        dual=max(q*slope-e.sig*e.fy(h,q) for q in e.knots)
        return dual-slope*t+e.sig*e.fy(h,t)
    tests=0
    for u,paid in list(prices.items())[::max(1,len(prices)//7)]:
        for v in sample+(u,):
            rhs = gap(v,Q(1,u[0]))
            if v[0]>=6:
                rhs += gap(u,Q(1,v[0]))
            assert paid+prices.get(v,0) <= rhs
            tests+=1
    return tests


def lifetime(ex,dmax,stock):
    hs,ys = censuses(dmax,stock),censuses(ex.tip,ex.tip_stock)
    values = {}
    extrema = (0,0)
    positive=states=tests=0
    for h in hs:
        d,c=h
        reward=ex.rows(h)
        if reward and h in ((38,2),(41,3),(dmax,stock)):
            tests+=direct_check(ex,h,reward)
        for y in ys:
            e,k=y
            gap = reward.get(y,0)
            if gap:
                positive+=1
                extrema=max(extrema,(gap,h[0]))
            # Rate denominator = 2*d*e. Retain both internal sends. The
            # remaining c-1 and k-1 neighbors contribute at most a half each.
            n=2*d*e
            gh=values.get(((d+1,c),y),0)
            gy=values.get((h,(e+1,k)),0)
            fixed=((n*(d-c)+2*d,gh),(n*(e-k)+2*e,gy),
                   (2*e*(d-c),values.get(((d,c+1),y),0)),
                   (2*d*(e-k),values.get((h,(e,k+1)),0)))
            variable=((d*e*(c-1),gh),(d*e*(k-1),gy))
            values[h,y]=integer_bellman((c-1)*gap,n,fixed,variable)
            states+=1
    # Each physical core birth creates one edge. Before its birth use the
    # promoter continuation; at birth sum both choices of endpoint as host.
    singles={}
    for d,c in censuses(max(dmax,ex.tip),max(stock,ex.tip_stock)):
        gain=singles.get((d+1,c),0)
        mother=(d,c+1)
        birth=values.get((mother,J),0)+values.get((J,mother),0)
        fixed=((2*d*(d-c),gain),(2*(d-c),singles.get(mother,0)+birth))
        singles[d,c]=integer_bellman(0,2*d,fixed,((c*d,gain),))
    exact=LO*Q(singles[J],GRID)
    print('edge states',states,'positive rewards',positive,'peak ticks and D',extrema,
          'physical pair comparisons',tests,'exact ticks/core',exact,'per core',
          float(exact/ex.snap.scale),flush=True)
    return exact


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',type=int,default=64)
    p.add_argument('--host-stock',type=int,default=5)
    p.add_argument('--tip',type=int,default=14)
    p.add_argument('--tip-stock',type=int,default=4)
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--hash',action='store_true')
    p.add_argument('--pilot',action='store_true')
    args=p.parse_args()
    for side in args.side or ('lower','upper'):
        with Snapshot(PRIVATE/f'quintet_{side}.bin',side,verify=args.hash) as s:
            assert 38 <= args.host < s.d2-1 and 2 <= args.host_stock <= s.K2
            ex=Exchange(Edges(s),args.tip,args.tip_stock)
            if args.pilot:
                for h in ((38,2),(38,3),(40,3),(48,3),(64,5),(100,3),(170,5)):
                    x=ex.rows(h)
                    print(side,h,'positive',sum(v>0 for v in x.values()),'of',len(x),
                          'max',max(x.values(),default=0),'ticks',flush=True)
            else:
                print(side,'host and neighbor caps',args.host,args.host_stock,args.tip,args.tip_stock,flush=True)
                lifetime(ex,args.host,args.host_stock)
            s.check_unchanged()
    print('PASS',flush=True)


if __name__=='__main__':
    main()
