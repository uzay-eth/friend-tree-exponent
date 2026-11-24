#!/usr/bin/env python3
"""Lower GOLD2, real extra-child H3 recovery with all other ports adverse."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path

GOLD=Path(__file__).resolve().parents[2]
SOURCE=GOLD/'higher/cohort.py'
PIN='b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c'
raw=SOURCE.read_bytes(); assert sha256(raw).hexdigest()==PIN
base=module_from_spec(spec_from_file_location('gold2_sibling_cohort',SOURCE))
exec(compile(raw,str(SOURCE),'exec'),base.__dict__)
GRID=base.GRID


class Cohort(base.Cohort):
    def __init__(self,r,a):
        super().__init__(r,a)
        self.exact_reward=lru_cache(maxsize=100000)(self.exact_reward)

    def exact_reward(self,h,m,y,s,u,kind):
        row=self.mask(h,m,y)
        if row is None:return Q(0)
        B,b,*_=row
        def gain(role,q):
            if q is None:return Q(0)
            x=(Q(1,self.D[role]+1) if b>=0 else Q(0)) if q==base.OVER else Q(1,q[0])
            return max(0,b/2)-b*x
        clip=lambda z:max(0,z-B)
        zu=gain('u',u); zs=gain('s',s)
        if kind=='box':return clip(0)
        if kind=='upstream':return clip(zu)-clip(0)
        assert kind=='sibling' and h[1]>2
        # Jensen: summing n copies with n*z_i in the argument can be an
        # upper, not lower, bound. Use one z_i and subtract the same base.
        return clip(zs)-clip(0)

    def reward(self,tags,kind):
        q=self.exact_reward(*(tags[r] for r in 'hmy'),tags.get('s'),tags.get('u'),kind)
        assert q>=0
        return GRID*q//1


def literal(rows,c):
    """Build actual donor stars; compare physical cap accounting and clock units."""
    phys=rows.h.physical
    for dh,ds,du in product((19,38,65,79),(2,3),(2,3)):
        # u--h--m--y, with two real siblings of m attached at h.
        edges=((6,0),(0,1),(1,2),(2,3),(1,4),(1,5))
        adj={i:set() for i in range(7)}
        for i,j in edges:adj[i].add(j);adj[j].add(i)
        marks={i:(d,len(adj[i])) for i,d in enumerate((du,dh,2,2,ds,3,2))}
        h,m,y,s,u=(marks[i] for i in (1,2,3,4,0))
        if not c.mask(h,m,y):continue
        old=phys.subsets(rows.p,adj,marks); path=(1,2,3)
        def w(data,k,ids):return data.get((k,tuple(sorted(ids))),(Q(0),))[0]
        f=w(old,'h3',path); P=-rows.h.lam*f
        for i in path:
            d,arity=marks[i]; L=d-arity
            ext=adj[i]-set(path)
            plus=phys.subsets(rows.p,*phys.grow(adj,marks,i,None))
            jump=w(plus,'h3',path)-f
            P+=(L+sum(Q(1,marks[j][0]) for j in adj[i]))*jump
            for j in ext:
                if any((k,tuple(sorted((*path,j)))) in old for k in ('p4','s3')):
                    P-=(Q(1,marks[j][0])-rows.h.center(d,arity))*jump
            if L:
                born=phys.subsets(rows.p,*phys.grow(adj,marks,i,7))
                P+=Q(L,d)*(w(born,'h3',path)-f+sum(w(born,k,(*path,7)) for k in ('p4','s3')))
            if i==2:
                for donor,keep in ((1,3),(3,1)):
                    P+=(Q(1,marks[donor][0])-rows.h.center(*m))*(w(plus,'g',(i,keep))-w(old,'g',(i,keep)))
        P*=rows.sig
        B,b,A,slopes=c.mask(h,m,y)
        rho=sum(Q(1,marks[i][0]) for i in (0,4,5))
        assert P==A+slopes[0]*rho
        total=c.exact_reward(h,m,y,None,None,'box')+c.exact_reward(h,m,y,None,u,'upstream')
        for i in (4,5):total+=c.exact_reward(h,m,y,marks[i],None,'sibling')
        assert 0<=total<=max(0,-P)
        den=c.den; fixed=dh-h[1]+Q(1,m[0])+Q(1,ds)
        residual=Q(1,du)+Q(1,3)
        assert (den*fixed).denominator==(den*residual).denominator==1
        assert den*(fixed+residual)==den*(dh-h[1]+sum(Q(1,marks[j][0]) for j in adj[1]))
        # For the constant robust subreward, scaling is exactly GRID*den;
        # no receive slope is shifted onto a gain continuation.
        reward=c.exact_reward(h,m,y,s,None,'sibling')
        if reward>0 and ds<=c.D['s'] and du<=c.D['u']:
            scaled=(GRID*reward//1)*den
            assert Q(scaled,GRID*den)<=reward and Q(scaled,GRID*den*rows.h.scale)<=reward/rows.h.scale
            print('literal original marks h,m,y,s1,s2,u',h,m,y,s,marks[5],u,'signed row ticks',P,
                  'box ticks',B,'s1 reward ticks',reward,'GRID',GRID,'potential scale',rows.h.scale,
                  'rate denominator',den,'fixed gain',fixed,'exterior gain',residual,
                  'scaled reward numerator',scaled,flush=True)
            return
    raise AssertionError('no nonzero literal sibling fixture')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key,v in [('mother',79),('min-mother',19),('mother-stock',8),('middle',5),('middle-stock',4),
                  ('young',9),('young-stock',3),('sibling',3),('upstream',3)]:
        p.add_argument('--'+key,type=int,default=v)
    p.add_argument('--kind',choices=('sibling','upstream','box','all'),default='sibling')
    a=p.parse_args()
    assert __debug__ and 19<=a.min_mother<=a.mother<=79
    a.side='lower'; a.outside_threshold=a.exclude_strict=False
    with base.host.open_host('lower') as h:
        assert h.docs['lower_sha256']=='41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413'
        print('original GOLD2 lower',h.docs['lower_sha256'],a,flush=True)
        r=base.Rows(h); base.registration(); base.physical_checks(r)
        c=Cohort(r,a)
        literal(r,c)
        result=Q(0)
        for k,orders in [('sibling',base.REG),('upstream',('uhmy',)),('box',('hmy',))]:
            if a.kind in (k,'all'):
                v=c.sweep(orders,k); result+=v
                print('exact',k,'ticks/core',v*h.scale,float(v*h.scale),flush=True)
        print('TOTAL exact ticks/core',result*h.scale,float(result*h.scale),'units/core',result,flush=True)


if __name__=='__main__':main()
