#!/usr/bin/env python3
"""Physical H3 occupation with a permanent, arity-tracked overflow arm."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import lcm
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from screen import Rows, LO, HI, J, bud, gain
from price import check_registration, solve_row

GRID=2**20


class I:
    """Closed rational interval; only used for the persistent reward arm."""
    def __init__(self,a,b=None):
        self.a,self.b=Q(a),Q(a if b is None else b)
        assert self.a<=self.b

    def __add__(self,o):
        o=iv(o)
        return I(self.a+o.a,self.b+o.b)

    __radd__=__add__

    def __neg__(self):
        return I(-self.b,-self.a)

    def __sub__(self,o):
        return self+-iv(o)

    def __rsub__(self,o):
        return iv(o)+-self

    def __mul__(self,o):
        o=iv(o)
        x=[a*b for a,b in product((self.a,self.b),(o.a,o.b))]
        return I(min(x),max(x))

    __rmul__=__mul__


def iv(x):
    return x if isinstance(x,I) else I(x)


class Reader(Rows):
    @lru_cache(maxsize=131072)
    def path(self,p):
        return super().path.__wrapped__(self,p)

    def h3_physical(self,p):
        if p[0]>p[2]:
            a,b=self._h3_physical(p[::-1])
            return a,b[::-1]
        return self._h3_physical(p)

    @lru_cache(maxsize=131072)
    def _h3_physical(self,p):
        """Use the same physical charges as Rows.row, without box payload reads."""
        s=self.s
        m=p[1]
        h=self.path(p)
        a=-self.lam*h
        B=[]
        for i,(d,c) in enumerate(p):
            adj=p[max(0,i-1):i]+p[i+1:i+2]
            ext=c-len(adj)
            jump=self.path(p[:i]+(gain(p[i]),)+p[i+1:])-h
            covered=ext and s.support(*self.extension(p,i,J)[:2])
            rate=d-c+sum(Q(1,D) for D,C in adj)+(ext*s.r1(p[i]) if covered else 0)
            a+=rate*jump
            B.append(0 if covered or not ext else self.sig*jump)
            if d>c:
                post=p[:i]+(bud(p[i]),)+p[i+1:]
                a+=Q(d-c,d)*(self.path(post)-h+s.value(*self.extension(post,i,J)))
        r=s.r1(m)
        for (d,c),tip in ((p[0],p[2]),(p[2],p[0])):
            a+=(Q(1,d)-r)*(self.path((gain(m),tip))-self.path((m,tip)))
        return self.sig*a,tuple(B)
    @lru_cache(maxsize=32768)
    def affine(self,key,inner,ends,axis,c):
        """a+b x where the unknown arm is at or past this level's cutoff."""
        j={'g':2,'h3':3,'p4':4,'s3':4}[key]
        g=self.s.G[j]
        e=list(ends)
        e[axis]=None
        a=self.s.value(key,inner,tuple(e))
        e[axis]=(g.d,c)
        b=g.d*(self.s.value(key,inner,tuple(e))-a)
        return a,b

    def tail_row(self,tracked,mid,c,x,rho):
        """Interval row when the other endpoint is beyond G2's cutoff.

        evaluate a gain jump as a difference of affine coefficients, so no
        D*width singularity enters the infinite stock limit.
        """
        s=self.s
        fit=lambda m,t,cc=c:self.affine('h3',(m,),(t,None),1,cc)
        def diff(f,g):
            return f[0]-g[0]+(f[1]-g[1])*x
        f=fit(mid,tracked)
        D,C=tracked; M,K=mid
        # Degree-180+ unknown tip cannot send an H3->P4 transfer or own a
        # newborn P4; its promotion can still change its H3 clamp class.
        cov_t=C>1 and s.support('p4',(tracked,mid))
        cov_m=K>2 and s.support('s3',(mid,))
        v=-self.lam*(f[0]+f[1]*x)
        v+=(D-C+Q(1,M)+((C-1)*s.r1(tracked) if cov_t else rho[0]))*diff(fit(mid,gain(tracked)),f)
        v+=(M-K+Q(1,D)+x+((K-2)*s.r1(mid) if cov_m else rho[1]))*diff(fit(gain(mid),tracked),f)
        # gamma = (1−c*x+x*(1/M+rho))/(1+x); the reciprocal is monotone.
        gamma=(1-c*x+x*(Q(1,M)+rho[2]))*I(1/(1+x.b),1/(1+x.a))
        v-=f[1]*x*gamma
        if D>C:
            nb=self.affine('p4',(bud(tracked),mid),(J,None),1,c)
            v+=Q(D-C,D)*(diff(fit(mid,bud(tracked)),f)+nb[0]+nb[1]*x)
        if M>K:
            nb=self.affine('s3',(bud(mid),),(tracked,None,J),1,c)
            v+=Q(M-K,M)*(diff(fit(bud(mid),tracked),f)+nb[0]+nb[1]*x)
        v+=(1-c*x)*diff(fit(mid,tracked,c+1),f)
        z=self.path((gain(mid),tracked))-self.path((mid,tracked))
        zz=diff(self.affine('g',(),(gain(mid),None),1,c),self.affine('g',(),(mid,None),1,c))
        center=s.r1(mid)
        v+=(x-center)*z+(Q(1,D)-center)*zz
        return self.sig*v

    def envelope_bound(self,dt,ct,m,k,do,co):
        """Exact separated robust upper row with explicit labeled endpoints."""
        a,b=self.h3_physical(((dt,ct),(m,k),(do,co)))
        return a+sum(max(0,x)*Q(c-d,2) for x,c,d in zip(b,(ct,k,co),(1,2,1)))


def optimum(reward,terms,discount,upper=False):
    # Reuse the production-derived FY adverse-corner ordering, but clear all
    # denominators once so per-state exact checks do not need Fraction sums.
    scale=lcm(*(z.denominator for lo,hi,_ in terms for z in (lo,hi)))
    q=discount.denominator
    a,b=reward*scale*q,discount.numerator*scale
    flex=[]
    for lo,hi,v in terms:
        ll=lo.numerator*(scale//lo.denominator)*q
        hh=hi.numerator*(scale//hi.denominator)*q
        a+=ll*v
        b+=ll
        if hi>lo:
            flex.append((v,hh-ll))
    x,y=a,b
    for v,w in sorted(flex,reverse=upper):
        if (v*b<=a if upper else v*b>=a):
            break
        a+=w*v
        b+=w
    z=-(-a//b) if upper else a//b
    assert z>=0
    g=x-y*z+sum(w*(max if upper else min)(0,k-z) for k,w in flex)
    assert g<=0 if upper else g>=0
    return z


def coupled(a,slopes,widths,den,num,successors,upper):
    """Fast affine pieces; delegate a genuine sign crossing to solve_row."""
    low=a+sum(min(0,b*w) for b,w in zip(slopes,widths))
    high=a+sum(max(0,b*w) for b,w in zip(slopes,widths))
    if low<0<high:
        return solve_row(a,slopes,[(0,w) for w in widths],den,num,successors,upper)[0]
    if high<=0:
        num-=GRID*a
        successors=tuple(v-GRID*b if w else v for v,b,w in zip(successors,slopes,widths))
    flexible=sorted(((v,w) for v,w in zip(successors,widths) if w),reverse=upper)
    a,b=num,den
    for v,w in flexible:
        if (v*b<=a if upper else v*b>=a):
            break
        a+=w*v
        b+=w
    z=-(-a//b) if upper else a//b
    g=num-den*z+sum(w*(max if upper else min)(0,v-z) for v,w in flexible)
    assert z>=0 and (g<=0 if upper else g>=0)
    return z


def states(D,C,overflow):
    return [(d,c) for d in range(2,D+1) for c in range(1,min(d,C)+1)]+(
        [(0,c) for c in range(1,C+1)] if overflow else [])


class Cohort:
    def __init__(self,r,a,shape):
        self.r,self.a,self.shape=r,a,shape
        self.crown=shape=='crown'
        self.parents=(0,0) if self.crown else (0,1)
        # chronological role index: old tip, middle, young tip for a chain;
        # middle, early child, late child for a crown.
        self.order=(1,0,2) if self.crown else (0,1,2)
        if a.upstream:
            self.parents=(0,1,2)
            self.order=(1,2,3)
        self.other= self.order[0 if a.tracked=='young' else 2]
        self.tracked= self.order[2 if a.tracked=='young' else 0]
        self.mid=self.order[1]
        caps=([a.upstream] if a.upstream else [])+[a.mid,a.mid,a.mid]
        caps[self.other],caps[self.tracked]=a.other,a.tip
        self.caps=caps
        self.persistent={self.other}|({0} if a.upstream_overflow else set())
        self.axes=[states(d,a.upstream_stock if a.upstream and i==0 else a.center_stock if i==self.mid else a.tracked_stock if i==self.tracked else a.tip_stock,i in self.persistent)
                   for i,d in enumerate(caps)]
        self.masks=0
        self.pos=0
        # Bound caches by sampled labels, not by the degree-cut volume.
        self.exact_cached=lru_cache(maxsize=16000)(r.envelope_bound)

    def gap(self,p):
        # Disjoint from the wider previous strict-chain payment (65,5,65).
        if self.crown or tuple(c for d,c in p)!=(2,2,1):
            return True
        return not (p[1][0]<=5 and all(0<d<=65 or (d==0 and self.a.other<65) for d,c in p[::2]))

    def pays(self,p):
        D,C=p[self.tracked]
        U=p[self.other][0]
        return (self.a.min_tip<=D<=1998 and C>=self.a.min_tracked_stock and 2<=U<=1998 and
                (not self.a.upstream_overflow or p[0][0]==0) and
                (not self.a.split or self.a.tracked=='old' or U<19) and
                self.gap(tuple(p[i] for i in self.order)))

    @lru_cache(maxsize=131072)
    def reward(self,p):
        D,C=p[self.tracked]; U,K=p[self.other]; M,B=p[self.mid]
        q=self.exact_cached(D,C,M,B,U,K)
        self.masks+=1
        self.pos+=q<0
        return max(0,-q)

    def exact_reward(self,p):
        out,slopes=self.r.h3_physical(tuple(p[i] for i in self.order))
        if self.a.upstream and p[0][0]:
            out+=slopes[0]*Q(1,p[0][0])
        if self.a.sibling_box:
            out+=max(0,slopes[0]*Q(p[self.tracked][1]-2,2))
        bmap=dict(zip(self.order,slopes))
        self.masks+=1
        return out,tuple(bmap.get(i,0) for i in range(len(p)))

    def accept(self,p,adj):
        return all(c>=len(adj[i])+(i==0) for i,(d,c) in enumerate(p))

    def run(self):
        nxt={}
        discount=LO if self.a.ceiling else HI
        n=len(self.order)+bool(self.a.upstream)
        for k in range(n,0,-1):
            adj=[[] for _ in range(k)]
            for j,i in enumerate(self.parents[:k-1],1):
                adj[i].append(j)
                adj[j].append(i)
            cur={}
            for p in product(*(reversed(ax) for ax in self.axes[:k])):
                if not self.accept(p,adj):
                    continue
                active=k==n and self.pays(p)
                exact=self.exact_reward(p) if active and self.a.coupled else None
                reward=self.reward(p) if active and exact is None else 0
                terms=[]
                gs=[]
                for i,(d,c) in enumerate(p):
                    if d:
                        g=(d+1,c)
                        if i in self.persistent and d==self.caps[i]:
                            g=(0,c)
                        lo=hi=Q(d-c)
                        for j in adj[i]:
                            dj=p[j][0]
                            if dj:
                                lo+=Q(1,dj); hi+=Q(1,dj)
                            else:
                                hi+=Q(1,self.caps[j]+1)
                        hi+=Q(c-len(adj[i]),2)
                        v=cur.get(p[:i]+(g,)+p[i+1:],0)
                        terms.append((lo,hi,v))
                        gs.append((lo,hi,v))
                        b0=b1=Q(d-c,d)
                    else:
                        # A degree gain stays within the permanent interval.
                        # Its value jump is identically zero, with no rate cap.
                        b0,b1=Q(max(0,self.caps[i]+1-c),self.caps[i]+1),Q(1)
                        gs.append((Q(0),Q(0),0))
                    pb=p[:i]+((d,c+1),)+p[i+1:]
                    b=cur.get(pb,0)
                    if k<n and i==self.parents[k-1]:
                        b+=nxt.get(pb+(J,),0)
                    terms.append((b0,b1,b))
                # Ignore rate variables with zero width naturally. Gains
                # and births strictly advance the represented state.
                if exact is None:
                    cur[p]=optimum((Q(reward)*GRID)//1,terms,discount,self.a.ceiling)
                else:
                    # solve_row uses 2**20 fractional ticks. Match that grid
                    # throughout all stages, including overflow and prefixes.
                    assert GRID==2**20
                    widths=[hi-lo for lo,hi,_ in gs]
                    a,b=exact
                    # Exact gain widths are already in gs. At a paying state,
                    # only the explicitly overflowed upstream birth adds a
                    # flexible variable, with zero slope in the H3 reward.
                    if self.a.upstream_overflow:
                        d,c=p[0]
                        assert d==0 and terms[0][0]==max(Q(0),1-Q(c,self.caps[0]+1))
                        lo,hi,v=terms[0]
                        widths.append(hi-lo)
                        b=(*b,0)
                        vs=[v0 for _,_,v0 in gs]+[v]
                    else:
                        assert sum(hi>lo for lo,hi,_ in terms)==sum(w>0 for w in widths)
                        vs=[v for _,_,v in gs]
                    if self.a.sibling_box:
                        # Freeze only the sibling contribution to the reward
                        # at its adverse corner; keep the actual sibling rate
                        # as an independent gain variable with the same jump.
                        i=self.tracked
                        w=Q(p[i][1]-2,2)
                        lo,hi,_=gs[i]
                        assert 0<=w<=hi-lo and widths[i]==hi-lo
                        widths[i]-=w
                        widths.append(w)
                        b=(*b,0)
                        vs.append(vs[i])
                    self.pos+=a+sum(max(0,z*w) for z,w in zip(b,widths))<0
                    cur[p]=coupled(a,b,widths,discount+sum(lo for lo,_,_ in terms),
                                   sum(lo*v for lo,_,v in terms),vs,self.a.ceiling)
            print(self.r.side,self.shape,'prefix',k,'states',len(cur),flush=True)
            nxt=cur
        return (HI if self.a.ceiling else LO)*Q(nxt[(J,)],GRID*self.r.s.h['S'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mid',type=int,default=6)
    p.add_argument('--tip',type=int,default=32)
    p.add_argument('--min-tip',type=int,default=19)
    p.add_argument('--other',type=int,default=8)
    p.add_argument('--tip-stock',type=int,default=3)
    p.add_argument('--tracked-stock',type=int)
    p.add_argument('--min-tracked-stock',type=int,default=1)
    p.add_argument('--center-stock',type=int,default=4)
    p.add_argument('--shape',choices=('chain','crown'),action='append')
    p.add_argument('--tracked',choices=('old','young'),default='old')
    p.add_argument('--split',action='store_true',help='young-target rows require older endpoint degree below 19')
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--check-tail',action='store_true')
    p.add_argument('--coupled',action='store_true')
    p.add_argument('--ceiling',action='store_true',help='upper bound the specified exact-tip payment, with zero reward after overflow')
    p.add_argument('--upstream',type=int,default=0,help='track the high mother’s nonfounder parent up to this degree')
    p.add_argument('--upstream-stock',type=int,default=2)
    p.add_argument('--upstream-overflow',action='store_true',help='pay only after that parent crosses its degree cutoff')
    p.add_argument('--sibling-box',action='store_true',help='price the baseline with sibling rates adverse in the reward only')
    a=p.parse_args()
    a.tracked_stock=a.tip_stock if a.tracked_stock is None else a.tracked_stock
    assert 2<=a.other<=1998 and 19<=a.min_tip<=a.tip<=1998 and 2<=a.mid<38
    assert 1<=a.tip_stock<=15 and 1<=a.tracked_stock<=26 and 2<=a.center_stock<=4
    assert 1<=a.min_tracked_stock<=a.tracked_stock
    assert not a.ceiling or a.coupled
    assert not a.upstream or (a.shape==['chain'] and a.tracked=='old' and a.coupled and
                             2<=a.upstream<=1998 and 2<=a.upstream_stock<=a.upstream)
    assert not a.upstream_overflow or a.upstream
    assert not a.sibling_box or (a.shape==['chain'] and a.tracked=='old' and a.coupled)
    assert not a.split or a.tracked=='old' or a.min_tip==19
    if a.check_tail:
        check_registration()
        solver_check()
    for side in a.side or ('lower','upper'):
        r=Reader(side)
        assert all(r.s.support('h3',((d,c),)) for d,c in states(a.mid,a.center_stock,False) if c>=2)
        if a.check_tail:
            tail_check(r)
            r.s.close()
            continue
        for shape in a.shape or ('chain','crown'):
            c=Cohort(r,a,shape)
            if shape=='crown' and min(a.mid,a.center_stock)<3:
                continue
            v=c.run()
            print(side,shape,a.tracked,'ceiling' if a.ceiling else 'floor',v,float(v),
                  'reward rows/negative throughout box',c.masks,c.pos,flush=True)
            c.reward.cache_clear()
        r.row.cache_clear()
        r._h3_physical.cache_clear()
        r.path.cache_clear()
        r.affine.cache_clear()
        r.s.close()


def tail_check(r):
    count=0
    for t,m,c,D in product(((19,2),(48,3),(180,3)),((3,2),(6,4)),(1,3,17),(180,181,257,1998,99999)):
        p=(t,m,(D,c))
        a,b,_,_=r.row(p)
        assert (a,b)==r.h3_physical(p)
        rho=tuple(Q(C-k,4) for (_,C),k in zip(p,(1,2,1)))
        exact=a+sum(x*y for x,y in zip(b,rho))
        u=r.tail_row(t,m,c,I(Q(1,D)),tuple(I(z) for z in rho))
        assert u.a==u.b==exact,(p,exact,u.a,u.b)
        U=r.tail_row(t,m,c,I(0,Q(1,180)),tuple(I(0,2*z) for z in rho))
        assert U.a<=exact<=U.b
        count+=1
    print(r.side,'tail physical algebra',count,'point equalities and interval bounds',flush=True)


def solver_check():
    count=0
    for a,b,w,upper in product((-2,0,2),product((-3,0,3),repeat=3),
                                product((Q(0),Q(1,2)),repeat=3),(False,True)):
        den,num,vs=Q(7,4),Q(3,4),(9,3,15)
        assert coupled(a,b,w,den,num,vs,upper)==solve_row(
            a,b,[(0,x) for x in w],den,num,vs,upper)[0]
        count+=1
    print('coupled solver',count,'exact comparisons with sign-hyperplane oracle',flush=True)


if __name__=='__main__':
    main()
