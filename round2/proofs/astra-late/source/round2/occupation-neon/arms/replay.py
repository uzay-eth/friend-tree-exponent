#!/usr/bin/env python3
"""Positive occupation of raw higher rows on small physical fringes."""
from argparse import ArgumentParser
from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb,prod
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'higher'))
from extract import Extract,Cert,NB,LO,HI,F_EDGE,F_DEN,SPEC,embedding,label,up,bud,physical,Q


class Extended(Extract):
    @lru_cache(maxsize=60000)
    def value(self,f,x):
        if all(sum(y)<self.c.D for y in x):return self.weight(f,x)
        return sum(self.weight(f,tuple(v for v,w in a))*prod(w for v,w in a)
                   for a in product(*(self.c.embed(y) for y in x)))

    def live(self,f,x):
        if f=='hV':return x[1][2]==1 and 1<=x[1][1]<=self.c.K
        if f=='hC':return 2<=x[0][1]<=self.c.K
        return self.supported(f,x)

    def center(self,y,par):
        c=self.c
        if sum(y)<c.D:
            return Q(int((c.rpar if par else c.rch)[c.state[y]]),c.Sr)
        k,h=y[1:];i=(1-h)*(c.K+2)+min(k,c.K+1)
        return Q(int((c.trpar if par else c.trch)[i]),c.Sr)

    def moved(self,f,x,i,par):
        c=self.c;y=x[i]
        if f=='g' and sum(y)>=c.D:
            k,h=y[1:]
            return k<=c.K and (i==0 and (k>=1 if par else k>=2) or i==1 and k>=1)
        return self.route(f,x,i,par)

    def transfer(self,f,x,i,par):
        """Transferred derivative; tail edge receivers use T2/T3 controls."""
        c=self.c;y=x[i];d=sum(y)
        if f!='g' or d<c.D:return self.value(f,label(x,i,up(y)))-self.value(f,x)
        t=Q(c.D,d);k,h=y[1:];a=(1-h)*(c.K+2)+k
        z=x[1-i]
        amp=Q(0)
        for v,w in c.embed(z):
            f0=int(c.g[c.state[('inf',)],c.state[v]] if i==0 else c.g[c.state[v],c.state[('inf',)]])
            t1=('T1',k,h)
            f1=int(c.g[c.state[t1],c.state[v]] if i==0 else c.g[c.state[v],c.state[t1]])
            ref=int(c.yp[a,c.child[v]] if i==0 else c.yw[c.par[v],k])
            amp+=w*(3*t*t*(1-t)*Q(f0-f1,c.D)+t**3*ref)
        return amp

    @lru_cache(maxsize=50000)
    def extended(self,p,x):
        A=B=Q(0);c=self.c
        for f,ids in embedding(p):
            if len(ids)<3 or ids[0]!=0:continue
            xx=tuple(x[i] for i in ids)
            if not self.live(f,xx):continue
            topo=SPEC[f];n=len(topo);own=self.value(f,xx);a=-c.lam*own;b=Q(0)
            for i,z in enumerate(xx):
                L,k,h=z;j=ids[i];D=sum(z)
                neigh=[v for v in range(n) if topo[v]==i]+([topo[i]] if i else [])
                delta=self.value(f,label(xx,i,up(z)))-own
                a+=(L+sum(Q(1,sum(xx[v])) for v in neigh))*delta
                if L:
                    born=label(xx,i,bud(z))+(NB,)
                    new=sum(self.value(g,tuple(born[q] for q in jj)) for g,jj in embedding(topo+(i,)) if len(jj)==n+1)
                    a+=Q(L,D)*(self.value(f,born[:-1])-own+new)
                for par,donors in ((True,[p[j]] if i==0 else []),
                                   (False,[v for v in range(len(x)) if p[v]==j and v not in ids])):
                    for donor in donors:
                        if self.moved(f,xx,i,par):a+=self.center(z,par)*delta
                        elif donor<0:b+=delta
                        else:a+=Q(1,sum(x[donor]))*delta
            for donor in range(n):
                for g,keep in embedding(topo):
                    if len(keep)!=n-1 or donor in keep:continue
                    j,=[u for u in keep if topo[u]==donor or topo[donor]==u]
                    small=tuple(xx[u] for u in keep); i=keep.index(j); par=topo[j]==donor
                    if not self.live(g,small):continue
                    assert self.moved(g,small,i,par)
                    a+=(Q(1,sum(xx[donor]))-self.center(xx[j],par))*self.transfer(g,small,i,par)
            A+=a*self.sign/c.S;B+=b*self.sign/c.S
        return A,B

    @lru_cache(maxsize=2000)
    def face(self,p,x,i):
        """All degrees >=18 in coordinate i; other stocks fixed until exit."""
        d0=self.c.D;k,h=x[i][1:]; m=7
        # On t=18/D, the exact source times (18+t)^3 is polynomial of
        # degree <=7: one cubic embedded slot; D+1 after a gain; affine
        # promotion/reciprocal; the physical linear gain cancels at t=0.
        # Interpolate in t with exact arithmetic, then use Bernstein.
        ts=[Q(d0,d0+j) for j in range(m+1)]
        rewards=[]
        for j,t in enumerate(ts):
            a,b=self.extended(p,label(x,i,(d0+j-k-h,k,h)))
            rewards.append(((d0+t)**3*a,(d0+t)**3*(a+b/2)))
        def coeff(j):
            poly=[Q(1)]
            for q,t in enumerate(ts):
                if q==j:continue
                b=[Q(0)]*(len(poly)+1)
                for l,y in enumerate(poly):b[l]-=t*y/(ts[j]-t);b[l+1]+=y/(ts[j]-t)
                poly=b
            return poly
        cs=[coeff(j) for j in range(m+1)]
        polys=[[sum(a[l]*r[end] for a,r in zip(cs,rewards)) for l in range(m+1)] for end in (0,1)]
        for d in (27,43):
            t=Q(d0,d);a,b=self.extended(p,label(x,i,(d-k-h,k,h)))
            assert all(sum(y*t**l for l,y in enumerate(poly))==(d0+t)**3*r for poly,r in zip(polys,(a,a+b/2)))
        for power in range(4):
            if any(any(y for y in poly[:power]) for poly in polys):break
            n=m-power
            low=min(sum(poly[power+l]*Q(comb(j,l),comb(n,l)) for l in range(j+1))
                    for poly in polys for j in range(n+1))/(d0+1)**3
            if low<=0:continue
            # Kill on every other gain or any promotion. Tail gains can
            # continue forever. Its discounted drift costs at most power.
            rate=Q(1,2)+len(x)+sum(y[0] for j,y in enumerate(x) if j!=i)
        rate+=sum(Q(1,sum(x[v])) if v!=i else Q(1,d0)
                      for j in range(len(x)) if j!=i for v in range(len(x)) if p[j]==v or p[v]==j)
            return low/(HI+rate+power),power
        return Q(0),0


def solve(E,cap,degree,tail,unbounded):
    from extract import shape
    U=10**23; out=lambda x:Q(x*U//1,U)
    payments=defaultdict(int)
    @lru_cache(maxsize=None)
    def V(z):
        p,x=shape(z,1); n=len(x);mx=max(map(sum,x)); N=sum(y[0]+1 for y in x)
        if n>4:return Q(0)
        if N>cap or mx>degree:
            if unbounded and n>=3 and mx>=18 and sum(sum(y)>=18 for y in x)==1:
                i=next(j for j,y in enumerate(x) if sum(y)>=18)
                if max(sum(y) for j,y in enumerate(x) if j!=i)<=5:
                    proxy=label(x,i,(18-x[i][1]-x[i][2],*x[i][1:]))
                    coef,power=E.face(p,proxy,i)
                    if coef:
                        payments['unbounded_face']+=1
                        return out(coef*Q(18,mx)**power)
            return Q(0)
        a=b=Q(0)
        if n>=3 and N>12:
            if mx+1<E.c.D:
                t,b=E.charges(z);a=sum(t.values())
                assert (a,b)==E.extended(p,x)
            elif tail:a,b=E.extended(p,x)
            assert min(a,a+b/2)>=0,(z,a,b)
            if a or b:payments['tail' if mx>=E.c.D else 'inner']+=1
        rates=physical.intrinsic_moves(z);den=HI+sum(r for _,r in rates)
        assert den-HI==N-Q(1,sum(x[0]))
        term=a+sum(r*V(y) for y,r in rates); plus=V(physical.plus(z))+b
        v=out(min(term/den,(term+plus/2)/(den+Q(1,2))))
        assert 0<=v and den*v<=term and (den+Q(1,2))*v<=term+plus/2
        return v
    v=LO*V(physical.newborn(2));print('stable-age addition outward',float(v),'exact',v,'contexts',V.cache_info().currsize,'paid',dict(payments),flush=True)
    return v


def main():
    p=ArgumentParser(description=__doc__);p.add_argument('blob',type=Path);p.add_argument('--cap',type=int,default=17)
    p.add_argument('--degree',type=int,default=19);p.add_argument('--tail',action='store_true');p.add_argument('--unbounded',action='store_true');a=p.parse_args()
    c=Cert(a.blob);side='lower' if c.mode==2 else 'upper';physical.MARK_CAP=2
    assert (c.D,c.K,c.S,c.lam)==(18,17,10**11,Q('.565673' if side=='lower' else '.565775'))
    v=solve(Extended(c),a.cap,a.degree,a.tail,a.unbounded)
    base=Q('.565675597658' if side=='lower' else '.565772677869')
    print(side,'additional endpoint move',float(v/F_DEN[side]),'exact updated endpoint',base+(1 if side=='lower' else -1)*v/F_DEN[side],flush=True)
    new=base+(1 if side=='lower' else -1)*v/F_DEN[side]; q=new*10**12
    t=q//1 if side=='lower' else -(-q//1)
    print('outward endpoint',f'{t//10**12}.{t%10**12:012d}')


if __name__=='__main__':main()
