#!/usr/bin/env python3
"""Paired U5 pendant-fringe renewal, with exact checked Bellman bounds."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from json import loads
from math import comb, prod
from os import close
from pathlib import Path
from subprocess import check_output
import sys

R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / 'attacks/neon2-quintet/physical'))
from check import PIN, WORK, Payload, physical_core, subsets

LO, HI = Q('0.565702177913'), Q('0.565748407438')
TRIAL = Q('0.565724')
BOX = (Q(0), Q(1, 2))
GRID = 2**64
FAMILIES = ('f','g','h3','p4','s3','p5','b5','s4')
SHAPES = {
    'p5': ((0,1), (1,2), (2,3), (3,4)),
    'b5': ((0,1), (1,2), (2,3), (2,4)),
    's4': ((0,1), (1,2), (1,3), (1,4)),
}


def adjacent(edges, n):
    a = {i: set() for i in range(n)}
    for u,v in edges:
        a[u].add(v)
        a[v].add(u)
    return a


def bump(ds, i):
    return ds[:i] + (ds[i]+1,) + ds[i+1:]


def mul(a,b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out


def fit(xs,ys):
    out = [Q(0)]*len(xs)
    for x,y in zip(xs,ys):
        fac = [y]
        for u in xs:
            if u!=x:
                fac = mul(fac,[-u/(x-u),1/(x-u)])
        out = [a+b for a,b in zip(out,fac)]
    return out


def peval(poly,x):
    out = Q(0)
    for a in reversed(poly):
        out = out*x+a
    return out


def bern(poly,a,b):
    n = len(poly)-1
    shifted = [(b-a)**i*sum(poly[j]*comb(j,i)*a**(j-i) for j in range(i,n+1))
               for i in range(n+1)]
    return [sum(shifted[i]*Q(comb(j,i),comb(n,i)) for i in range(j+1)) for j in range(n+1)]


class Fringe:
    def __init__(self, edges, witnesses):
        self.a = adjacent(edges, 5)
        self.c = tuple(len(self.a[i])+(i==0) for i in range(5))
        self.low = tuple(max(2,c) for c in self.c)
        self.ps = witnesses
        self.norm = [Q(p.row('f',(0,),(),{0:(2,1)})) for p in witnesses]
        assert min(self.norm) > 0

    @lru_cache(None)
    def weights(self, j, ds, promoter=-1):
        a = self.a
        m = dict(enumerate(zip(ds,self.c)))
        p = self.ps[j]
        if promoter >= 0:
            a = dict(a)
            a[5] = {promoter}
            a[promoter] = a[promoter] | {5}
            m[promoter] = ds[promoter],self.c[promoter]+1
            m[5] = 2,1
        w = defaultdict(Q)
        for (kind,_),(value,_) in subsets(p,a,m).items():
            w[kind] += value/self.norm[j]
        return w

    def hazards(self, ds, x):
        g = tuple(d-c+sum((Q(1,ds[v]) for v in self.a[i]),Q())+(x if i==0 else 0)
                  for i,(d,c) in enumerate(zip(ds,self.c)))
        b = tuple(Q(d-c,d) for d,c in zip(ds,self.c))
        return g,b

    @lru_cache(None)
    def raw(self, j, ds, x):
        g,b = self.hazards(ds,x)
        w = self.weights(j,ds)
        # Sum true jumps of the within-fringe potential, including within-fringe
        # newborns. No row centering, cap slack, or clipping of successors.
        nxt = [self.weights(j,bump(ds,i)) for i in range(5)]
        new = [self.weights(j,ds,i) if b[i] else w for i in range(5)]
        ans = {}
        for k in FAMILIES:
            ans[k] = (-TRIAL*w[k] + sum(g[i]*(nxt[i][k]-w[k])+b[i]*(new[i][k]-w[k])
                                       for i in range(5)))
        return ans

    @lru_cache(None)
    def values(self, ds, x):
        a,b = (self.raw(j,ds,x) for j in range(2))
        A,B = sum(a.values()),sum(b.values())
        return A,B,A+B,sum(a[k]+b[k] for k in FAMILIES[-3:]),Q(1),abs(A+B),abs(A)+abs(B)

    def solve(self, d0, d1, extra, reward_start):
        """Discount killed on any promotion, cap-crossing or gain past d1.

        For abs(sum) add its exact zero to the x corners. There is a single
        common x in every signed reward, continuation hazard and denominator.
        """
        values = {}
        size = 0
        zero = (Q(0),)*7
        for ds in product(range(d1,d0-1,-1),*(range(d+extra,d-1,-1) for d in self.low[1:])):
            v = [[],[]]
            nex = [values.get(bump(ds,i),(zero,zero)) for i in range(5)]
            corners = list(BOX)
            if ds[0]>=reward_start:
                for k in range(3):
                    a,b = (self.values(ds,x)[k] for x in BOX)
                    if a*b<0:
                        corners.append(BOX[1]*a/(a-b))
            for x,disc in product(corners,(LO,HI)):
                g,b = self.hazards(ds,x)
                den = disc+sum(g)+sum(b)
                row = self.values(ds,x) if ds[0]>=reward_start else zero
                for side in range(2):
                    v[side].append(tuple((r+sum(g[i]*nex[i][side][k] for i in range(5)))/den
                                         for k,r in enumerate(row)))
                size += 1
            lower = tuple(Q(min(row[k] for row in v[0])*GRID // 1,GRID) for k in range(7))
            upper = tuple(Q(-(-max(row[k] for row in v[1])*GRID // 1),GRID) for k in range(7))
            assert all(all(x<=y for x,y in zip(lower,r)) for r in v[0])
            assert all(all(x>=y for x,y in zip(upper,r)) for r in v[1])
            values[ds] = lower,upper
        return values,size

    def old_edge(self, ds, x):
        """Raw existing 0--1 edge action, excluding any newborn motif."""
        g,b = self.hazards(ds,x)
        def w(p,d,c):
            return Q(p.row('g',(),(0,1),dict(enumerate(zip(d[:2],c[:2])))),p.h['S'])
        ans = []
        for p in self.ps:
            old = w(p,ds,self.c)
            ans.append(-TRIAL*old+sum(g[i]*(w(p,bump(ds,i),self.c)-old)
                                      +(b[i]*(w(p,ds,bump(self.c,i))-old) if b[i] else 0) for i in (0,1)))
        return sum(ans)

    def tail(self, stock):
        def coefficients(fn,x):
            def exact(d):
                return fn((d,*stock),x)*(1+Q(1,d))**2
            ds = (2000,2001,2003,2011,2077,2221)
            poly = fit([Q(1,d) for d in ds],[exact(d) for d in ds])
            for d in (2012,2027,3000,8000,31001):
                assert peval(poly,Q(1,d))==exact(d)
            return poly

        for label,fn in (('all',lambda ds,x:self.values(ds,x)[2]),
                         ('quintet',lambda ds,x:self.values(ds,x)[3]),('edge01',self.old_edge)):
            a,b = [coefficients(fn,x) for x in BOX]
            limit = a[0]
            assert b[0]==limit
            # r(D,x)-r(infty) = y*(linearized correction numerator)/(1+y)^2.
            hulls = [bern([c[1]-2*limit,c[2]-limit,*c[3:]],Q(0),Q(1,2000)) for c in (a,b)]
            slope = [z/den for h in hulls for z in (min(h),max(h)) for den in (1,(1+Q(1,2000))**2)]
            finite = [bern(c,Q(1,2011),Q(1,2004)) for c in (a,b)]
            far = [bern(c,Q(0),Q(1,100000)) for c in (a,b)]
            if label=='all':
                assert max(map(max,finite))<0<min(map(min,far))
            print('tail',label,'limit',fmt(limit),'linear slope hull',fmt(min(slope)),fmt(max(slope)),
                  'near and far numerator sign hull',*[tuple(map(fmt,(min(map(min,h)),max(map(max,h))))) for h in (finite,far)])

    def entrance(self):
        """One prescribed history from the newborn, uniformly in parent receipt."""
        a = {0:set()}
        ds = [2]
        factors = []
        def step(i,promotion):
            cs = [len(a[u])+(u==0) for u in a]
            gains = [ds[u]-cs[u]+sum((Q(1,ds[v]) for v in a[u]),Q()) for u in a]
            buds = [Q(ds[u]-cs[u],ds[u]) for u in a]
            rate = buds[i] if promotion else gains[i]
            assert rate>0
            factors.append(rate/(HI+sum(gains)+sum(buds)+BOX[1]))
            if promotion:
                j=len(a)
                a[j]={i}; a[i].add(j); ds.append(2)
            else:
                ds[i]+=1

        for j in range(1,5):
            i, = self.a[j].intersection(a)
            if ds[i]==len(a[i])+(i==0):
                step(i,False)
            step(i,True)
        assert tuple(ds)==self.low and a==self.a
        while ds[0]<1999:
            step(0,False)
        return prod(factors)

    def arrow_control(self, ds, parent_d):
        # Literal ordinary tree. The parent and its own ordinary leaves are
        # excluded from the reward, but its physical arrows are still counted.
        tree = [set(self.a[i]) for i in self.a] + [{0}]
        tree[0].add(5)
        for i,d in enumerate((*ds,parent_d)):
            for _ in range(d-len(tree[i])):
                tree[i].add(len(tree)); tree.append({i})
        rates = [Q(0)]*len(tree)
        for row in tree:
            for v in row:
                rates[v] += Q(1,len(row))
        assert sum(rates)==len(tree)
        def weights(p,t,keep):
            core,m = physical_core(t)
            core = {i:core[i]&keep for i in keep}
            w = defaultdict(Q)
            for (kind,_),(value,_) in subsets(p,core,m).items():
                w[kind] += Q(value,p.h['S'])
            return w

        keep = set(range(5))
        for j,p in enumerate(self.ps):
            w = weights(p,tree,keep)
            act = {k:-TRIAL*w[k] for k in FAMILIES}
            for i in keep:
                # All ordinary leaves of i have the same action/hazard.
                leaves = [v for v in tree[i] if len(tree[v])==1]
                groups = [(i,1)] + ([(leaves[0],len(leaves))] if leaves else [])
                for v,mult in groups:
                    t = tree.copy()+[{v}]
                    t[v] = t[v] | {len(tree)}
                    nxt = weights(p,t,keep | ({v} if v!=i else set()))
                    for k in FAMILIES:
                        act[k] += mult*rates[v]*(nxt[k]-w[k])
            assert act==self.raw(j,ds,Q(1,parent_d))


def fmt(v):
    return f'{float(v):+.12e}'


def report_bound(bounds, k):
    lo,hi = bounds[0][k],bounds[1][k]
    return f'[{fmt(lo)}, {fmt(hi)}] exact=[{lo}, {hi}]'


def main():
    assert check_output(['git','-C',WORK/'source','rev-parse','HEAD'],text=True).strip()==PIN
    paths = [WORK/f'quintet_{name}.bin' for name in ('lower','upper')]
    stats = [x.stat() for x in paths]
    ps = [Payload(x) for x in paths]
    assert [(p.h['mode'],p.h['P'],p.h['Q'],p.h['DF']) for p in ps]==[(2,56570237,10**8,2000),(1,56574828,10**8,2000)]
    gold = (WORK/'gold.json').read_bytes()
    from hashlib import sha256
    assert sha256(gold).hexdigest()=='d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb'
    assert loads(gold)
    for shape,edges in SHAPES.items():
        f = Fringe(edges,ps)
        for d,parent_d in ((1999,2),(2000,2000),(2004,7)):
            f.arrow_control((d,*f.low[1:]),parent_d)
        print(shape,'48 literal family/side/seam matches')
        print(shape,'stock',f.c,'bare',f.low,'seed words',f.norm)
        paired,checks = f.solve(1999,2003,1,2000)
        omitted,more = f.solve(1999,2011,1,2004)
        print('corner/discount checks',checks,more)
        for stock in (f.low[1:], tuple(d+1 for d in f.low[1:])):
            start = (1999,*stock)
            print('entry',start, 'raw A B sum quintet at x=0 and 1/2')
            for d in (1999,2000,2003,2004,2011):
                print(d, *[tuple(map(fmt,f.values((d,*stock),x)[:4])) for x in BOX])
            for k,name in enumerate(('A','B','sum','quintets','time','abs-sum','sum-abs')):
                print(name, report_bound(paired[start],k))
            print('omitted sum',report_bound(omitted[start],2))
            print('omitted abs-sum',report_bound(omitted[start],5))
            w = paired[start]
            assert omitted[start][0][5] > w[1][5]
            print('paired vs separate width',fmt(w[1][2]-w[0][2]),fmt(w[1][0]+w[1][1]-w[0][0]-w[0][1]))
            f.tail(stock)
        entry=f.entrance()
        absfloor=omitted[(1999,*f.low[1:])][0][5]
        assert absfloor>0
        print('discounted entry floor',fmt(entry),'omitted absolute per-core floor via prescribed path',fmt(LO*entry*absfloor))
        f.values.cache_clear(); f.raw.cache_clear(); f.weights.cache_clear()
    for p,path,s in zip(ps,paths,stats):
        now=path.stat()
        assert (now.st_dev,now.st_ino,now.st_size,now.st_mtime_ns,now.st_ctime_ns)==(s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)
        print(path.name,'sparse words',len(p.reads))
        close(p.fd)


if __name__=='__main__':
    main()
