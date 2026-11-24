#!/usr/bin/env python3
"""Original quintet g residual after removing every physical FY gap."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from replay import Edges, J, PRIVATE, Snapshot

GRID = 10**8
GOLD = 'd26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb'


def solve(corners):
    assert all(isinstance(n,(int,Q)) and isinstance(d,Q) and n >= 0 and d > 0 for n,d in corners), corners
    v = Q(min((n*GRID)//d for n,d in corners), GRID)
    assert all(v*d <= n for n,d in corners)
    return v


class PostFY:
    def __init__(self, edge, hi):
        self.e, self.s, self.hi = edge, edge.snap, hi
        self.d, self.k, self.sig = self.s.d2, self.s.K2, edge.sig
        self.cap = edge.cap
        assert self.d == 180 and self.k == 26 and edge.knots[:2] == (0,Q(1,self.d))
        assert edge.h['d3e'] <= self.d
        self.tails = range(1,self.k+2)
        # Check all stored double-tail corners, not just the INF,INF word.
        assert all(self.s.g(i,j) == 0 for i,j in product(range(self.s.tail_base,self.s.inf+1),repeat=2))

    @cache
    def knots(self, x):
        return tuple((z, self.sig*self.e.fy(x,z)) for z in self.e.knots)

    def maximum(self, x, jump):
        return max(z*self.sig*jump-a for z,a in self.knots(x))

    def relaxed(self, x, y):
        """Return exact signed worst FY row at physical x,y; leave other donors live."""
        base = -self.sig*self.s.lam*self.e.g(x,y)
        slopes = []
        for p,q in ((x,y),(y,x)):
            a,b,fy = self.e.end(p,q)
            n = p[1]-1
            base += self.sig*a + (n*self.maximum(p,b) if fy else 0)
            slopes.append(0 if fy else self.sig*n*b)
        return base, slopes

    def finite(self, limit, arm):
        """Pay max(D)>18 in a downward-closed exact-degree pair mask."""
        assert 19 <= limit < self.d and 0 <= arm <= limit
        # Optional wider rectangular wing, one endpoint through arm, all arities there.
        wing = self.d-1 if arm else limit
        states = [(d,c) for d in range(wing,1,-1) for c in range(d,0,-1)]
        low = [(d,c) for d in range(arm,1,-1) for c in range(d,0,-1)]
        inside = [(d,c) for d in range(limit,1,-1) for c in range(d,0,-1)]
        value = {}

        def get(x,y):
            return value.get((x,y) if x>=y else (y,x),0)

        n = 0
        for i,x in enumerate(states):
            d,c = x
            ys = low if d>limit else inside[i-(len(states)-len(inside)):]
            for y in ys:
                f,k = y
                bx,by = Q(d-c,d),Q(f-k,f)
                fixed = bx*get((d,c+1),y)+by*get(x,(f,k+1))
                vx,vy = get((d+1,c),y),get(x,(f+1,k))
                rx,ry = d-c+Q(1,f),f-k+Q(1,d)
                if d > 18:
                    worst,(sx,sy) = self.relaxed(x,y)
                    reward = self.cap-worst
                    assert reward >= (max(0,sx)+max(0,sy))/2, (x,y,reward,sx,sy)
                else:
                    reward = sx = sy = Q(0)
                a = [(rx+(c-1)*u, sx*u) for u in ((0,Q(1,2)) if c>1 else (0,))]
                b = [(ry+(k-1)*v, sy*v) for v in ((0,Q(1,2)) if k>1 else (0,))]
                value[x,y] = solve([(reward+fixed+r*vx+s*vy-dx-dy, self.hi+bx+by+r+s)
                                      for (r,dx),(s,dy) in product(a,b)])
                n += 1
        entrants = {(d,c): get((d,c+1), J) for d,c in states}
        return entrants,n

    def tail_row(self, c, y):
        """Mixed verifier knot with live tail gamma,beta and finite donor mean."""
        e,s = self.e,self.s
        dy,cy = y
        x = self.d,c
        assert 2 <= dy < self.d and 1 <= cy <= dy
        g1 = e.g(x,y)
        j = e.weights(y,s.inner,self.d,self.k)[0][0]
        j_up = e.weights((dy+1,cy),s.inner,self.d,self.k)[0][0]
        # At the seam dy+1=d2, the INF,tail word is zero as checked above.
        g0 = s.g(s.inf,j)
        jump0 = s.g(s.inf,j_up)-g0
        base,slope,fy = e.end(y,x)
        jump1 = e.g((dy+1,cy),x)-g1
        # end() included the adjacent tail's physical send at t=1. In the
        # verifier lift this send is carried by the tail column, not INF.
        base = self.sig*base-Q(self.sig*jump1,self.d)
        base += Q(max(self.sig*jump0,self.sig*jump1),self.d)
        if fy:
            base += (cy-1)*self.maximum(y,jump1)
            slope = 0
        constant = base-self.sig*s.lam*g1
        dg = self.sig*(g0-g1)
        db = self.sig*(e.g((self.d,c+1),y)-g1) if c<=self.k else Q(0)
        return constant,dg,db,self.sig*(cy-1)*slope

    def tail_inf(self, y):
        """INF corner: finite endpoint only, excluding the tail's adjacent send."""
        d,c = y
        s,e = self.s,self.e
        def g(p):
            i = e.weights(p,s.inner,self.d,self.k)[0][0]
            return s.g(i,s.inf)
        old = g(y)
        j = g((d+1,c))-old
        center = s.center(d,c)
        base = (d-c)*j-s.lam*old
        if d>c:
            post = (d,c+1)
            m = e.mids.get(post)
            h3 = 0 if m is None else s.item(e.off_h3 + m*e.hstride + (e.n3-1)*e.n3//2)
            base += Q(d-c,d)*(g(post)-old+h3)
        if center is not None:
            base -= center*(s.f(d+1,c)-s.f(d,c))
        if y in e.mids:
            base += (c-1)*center*j
            return self.sig*base,0
        if center is not None and c>=2:
            return self.sig*base+(c-1)*(self.sig*e.fy(y,0)+self.maximum(y,j)),0
        return self.sig*base,self.sig*(c-1)*j

    def tail(self, tip):
        """t U(c,y) on one unbounded G2 tail joined to a finite low endpoint."""
        assert 2 <= tip < self.d
        marks = [(d,c) for d in range(tip,1,-1) for c in range(d,0,-1)]
        u = {}
        for y in marks:
            f,k = y
            upper = Q(f-k,f)
            r0 = self.tail_inf(y)
            for c in reversed(self.tails):
                z = c<=self.k
                gammas = (Q(self.d-c+Q(1,f),self.d+1),Q(1)) if z else (0,1)
                betas = (Q(self.d-c,self.d),Q(1)) if z else (0,)
                base,dg,db,slope = self.tail_row(c,y)
                uy = u.get((c,(f+1,k)),0)
                by = upper*u.get((c,(f,k+1)),0)
                corners = []
                for gamma,beta,zext,t in product(gammas,betas,(0,Q(1,2)),(0,1)):
                    assert self.cap >= r0[0]+zext*r0[1]
                    reward = self.cap-base-gamma*dg-beta*db-zext*slope
                    assert reward >= 0, (c,y,gamma,beta,zext,reward)
                    r = f-k+(k-1)*zext+Q(t,self.d)
                    count = reward + by + beta*u.get((c+1,y),0)+r*uy
                    corners.append((count,self.hi+gamma+beta+upper+r))
                u[c,y] = solve(corners)
        return u

    def bridge(self, u, tip, stock):
        """Carry the mixed-tail continuation back across 179 -> 180."""
        assert 1 <= stock <= self.k
        states = [(d,c) for d in range(self.d-1,1,-1) for c in range(min(d,stock),0,-1)]
        tips = [(d,c) for d in range(tip,1,-1) for c in range(d,0,-1)]
        # Both-low states need both full arities, since either end can grow
        # first; once a side grows past tip it uses the stock cap.
        order = [x for x in states if x[0]>tip]+tips
        value = {}

        def get(x,y):
            if y>x:
                x,y = y,x
            if x[0]==self.d:
                return u.get((min(x[1],self.k+1),y),0)
            return value.get((x,y),0)

        for i,x in enumerate(order):
            d,c = x
            for y in tips if d>tip else order[i:]:
                f,k = y
                bx,by = Q(d-c,d),Q(f-k,f)
                fixed = bx*get((d,c+1),y)+by*get(x,(f,k+1))
                vx,vy = get((d+1,c),y),get(x,(f+1,k))
                rx,ry = d-c+Q(1,f),f-k+Q(1,d)
                value[x,y] = solve([(fixed+r*vx+s*vy,self.hi+bx+by+r+s)
                                     for r,s in product((rx,rx+Q(c-1,2)),(ry,ry+Q(k-1,2)))])
        return {(d,c):get((d,c+1),J) for d,c in order},len(value)

    def fixtures(self):
        checked = Counter()
        for x,y in product(((19,1),(27,27),(38,2),(179,2),(179,27),(180,1),(180,27),(181,2),
                            (370,2),(370,27),(370,180),(20000,299)), repeat=2):
            g, slopes = self.relaxed(x,y)
            assert self.relaxed(y,x) == (g,slopes[::-1])
            for z,w in product((0,Q(1,180),Q(1,2)),repeat=2):
                raw = self.e.row(x,y,[z]*(x[1]-1),[w]*(y[1]-1))
                gaps = 0
                for p,q,r in ((x,y,z),(y,x,w)):
                    _,jump,fy = self.e.end(p,q)
                    if fy:
                        gaps += (p[1]-1)*(self.maximum(p,jump)-self.sig*(r*jump-self.e.fy(p,r)))
                actual = g+z*slopes[0]+w*slopes[1]
                assert gaps >= 0 and actual == raw+gaps <= self.cap
                if min(x[0],y[0])>=self.d:
                    assert actual == 0
                    checked['double'] += 1
                elif x[0]>=self.d>y[0]:
                    d,c = x
                    c0 = min(c,self.k+1)
                    t = Q(self.d,d)
                    gamma = (d-c+Q(1,y[0])+(c-1)*z)/(d+1)
                    beta = Q(d-c,d) if c<=self.k else 0
                    base,dg,db,sy = self.tail_row(c0,y)
                    inf,iy = self.tail_inf(y)
                    enveloped = t*(base+gamma*dg+beta*db+w*sy)+(1-t)*(inf+w*iy)
                    assert actual <= enveloped <= self.cap
                    checked['mixed'] += 1
                else:
                    checked['other'] += 1
        return dict(checked)


def promoter(hi, lo, exits, newborn_edge, terminal=None):
    """Physical birth input, with an optional unbounded tail subsolution."""
    value = {}
    for d,c in exits:
        b = Q(d-c,d)
        birth = newborn_edge(d,c)
        gain = value.get((d+1,c),0 if terminal is None else terminal(d+1,c))
        bud = value.get((d,c+1),0)
        value[d,c] = solve([(b*(birth+bud)+r*gain,hi+b+r) for r in (d-c,d-c+Q(c,2))])
    return lo*value[J]


def tail_promoter(post, lo, u, entries):
    """At physical D>=d2, use (d2/D)*W(c); normalize births only once."""
    D,k,hi = post.d, post.k, post.hi
    w = {k+1:Q(0)}
    # The high-arity continuation is safely left unpaid. The pair subsolution
    # still covers the entire overflow class; this only clips birth credit.
    for c in range(k,0,-1):
        b_range = (Q(D-c,D),Q(1))
        newborn = u[c+1,J]
        w[c] = solve([(b*(w[c+1]+newborn),hi+b+gamma) for b,gamma in product(b_range,(0,1))])
    states = [(d,c) for d in range(D-1,1,-1) for c in range(min(d,k),0,-1)]
    return promoter(hi,lo,states,lambda d,c:entries.get((d,c),0),lambda d,c:w[c] if d==D else 0)


def double_tail_floor(post,lo):
    """Both endpoints in G2 tails: constant payoff, with a killed entry route."""
    D,hi = post.d,post.hi
    # One end already lies in the tail. Keep the second end at c=1 until
    # it crosses D, killing on its first promotion. The first cannot leave.
    v = Q(post.cap,hi)
    for d in range(D-1,1,-1):
        b = Q(d-1,d)
        old = v
        v = Q(d-1)*old/(hi+b+d-1)
        assert 0 < v < old  # increasing the actual gain rate only helps
    # Register edges on tail-promoter births only. Use c=1 until promotion
    # for this auxiliary singleton; the first promotion kills the test.
    w = Q(D-1,D)*v/(hi+Q(D-1,D))
    assert 0 < w < v
    for d in range(D-1,1,-1):
        old = w
        w = Q(d-1)*old/(hi+Q(d-1,d)+d-1)
        assert 0 < w < old
    return lo*w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--finite',type=int,default=24)
    p.add_argument('--wing',type=int,default=0)
    p.add_argument('--tip',type=int,default=8)
    p.add_argument('--bridge-stock',type=int,default=8)
    p.add_argument('--side',choices=('lower','upper','both'),default='both')
    args = p.parse_args()
    assert __debug__
    raw = (PRIVATE/'gold.json').read_bytes()
    assert sha256(raw).hexdigest() == GOLD
    gold = json.loads(raw)
    lo,hi = (Q(gold[s]['certified_exact']) for s in ('lower','upper'))
    for side in ('lower','upper') if args.side=='both' else (args.side,):
        with Snapshot(PRIVATE/f'quintet_{side}.bin',side) as snapshot:
            post = PostFY(Edges(snapshot),hi)
            assert post.cap == gold[side]['cap_g']
            print(side,'fixtures',post.fixtures(),flush=True)
            entrants, n = post.finite(args.finite,args.wing)
            states = sorted((k for k,v in entrants.items() if k[0] <= (post.d-1 if args.wing else args.finite)),reverse=True)
            finite = promoter(hi,lo,states,lambda d,c:entrants[d,c])
            print(side,'finite',args.finite,'wing',args.wing,'edge states',n,'ticks/core',finite,float(finite),flush=True)
            u = post.tail(args.tip)
            entries,n = post.bridge(u,args.tip,args.bridge_stock) if args.bridge_stock else ({},0)
            tail = tail_promoter(post,lo,u,entries)
            print(side,'mixed unbounded, tip',args.tip,'edge states',len(u),'bridge states',n,
                  'ticks/core',tail,float(tail),flush=True)
            double = double_tail_floor(post,lo)
            print(side,'double-tail normalized floor',float(double/snapshot.scale),flush=True)
            # Round the long exact product outward before printing a compact receipt.
            total = finite+tail+Q(double*10**12//1,10**12)
            print(side,'total ticks/core',total,'normalized',float(total/snapshot.scale),
                  'fraction of 5e-6',float(total/(snapshot.scale*Q(5,10**6))),flush=True)
            snapshot.check_unchanged()


if __name__ == '__main__':
    main()
