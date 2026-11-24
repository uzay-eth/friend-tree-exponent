#!/usr/bin/env python3
"""Stationary physical core marginals: sparse LP proposal, exact dual repair."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import argparse
import gzip
import json
import math
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'renewal'),str(HERE.parent/'rows')]
from check import row_rewards
from check_edges import rows as edge_rewards

A,B=Q('.565673'),Q('.565775')
UNIT=10**12
C0lo,C0hi=Q('1.84124847315'),Q('1.95820675182')


class Model:
    def __init__(self,s,N,age=False):
        self.s,self.N=s,N
        types=[(d,k) for d in range(2,N+1) for k in range(d)]+[(0,0)]
        self.types=types; self.tail=types[-1]; self.virgin=types[0]
        self.bounds=[]; self.names=[]; self.eq=[]; self.le=[]
        self.x={u:self.var(('vertex',u),1) for u in types}
        self.y={}
        self.a={}
        for u,v in product(types,repeat=2):
            if u[0] and not u[1]:continue
            self.y[u,v]=self.var(('edge',u,v),1)
        # Two independent overflow donor moments per orientation when needed.
        # Each is the reciprocal of the appropriate *same* endpoint, integrated
        # against the edge cell, never a free replacement for a finite degree.
        for (u,v),i in self.y.items():
            for j,w in enumerate((u,v)):
                if w==self.tail:
                    z=self.var(('reciprocal',u,v,j),Q(1,N+1))
                    self.a[u,v,j]=z
                    self.add(self.le,{z:1,i:-Q(1,N+1)},0)
        self.add(self.eq,{i:1 for i in self.x.values()},1)
        for u in types:
            self.add(self.eq,{**{i:1 for (p,v),i in self.y.items() if v==u},self.x[u]:-1},0)
            if u[0] and u[1]:
                self.add(self.eq,{**{i:1 for (p,v),i in self.y.items() if p==u},self.x[u]:-u[1]},0)
        # Birth of cores = s. In the overflow donor bin only 0 <= b <= 1
        # is used; the extra constraints below separately drain degree shells.
        self.birth={u:self.bud(u) for u in types if u[0]}
        birth={self.x[u]:r for u,r in self.birth.items() if r}
        self.add(self.le,birth,s)
        self.add(self.le,{**{i:-v for i,v in birth.items()},self.x[self.tail]:-1},-s)
        for u in types[:-1]:
            # Projection of the full stationary generator onto the vertex.
            c=defaultdict(Q)
            for v in types[:-1]:
                if v==u:c[self.x[v]]-=s
                for nxt,r in ((self.gain(v),self.leaves(v)),(self.prom(v),self.bud(v))):
                    c[self.x[v]]+=r*(int(nxt==u)-int(v==u))
            for (p,w),i in self.y.items():
                for j,r in enumerate((p,w)):
                    if r[0]:
                        coef=int(self.gain(r)==u)-int(r==u)
                        if coef:self.donor(c,p,w,1-j,coef)
            self.add(self.eq,c,-s if u==self.virgin else 0)
        self.shell_flows()
        if age:self.age_cut()

    def var(self,n,u):
        j=len(self.names); self.names.append(n); self.bounds.append(Q(u)); return j

    @staticmethod
    def add(rows,c,b):rows.append(({j:Q(v) for j,v in c.items() if v},Q(b)))

    @staticmethod
    def leaves(u):return u[0]-u[1]-1

    def bud(self,u):return Q(self.leaves(u),u[0])

    def gain(self,u):return (u[0]+1,u[1]) if u[0]<self.N else self.tail

    def prom(self,u):return (u[0],u[1]+1)

    def donor(self,c,p,w,j,coef):
        u=(p,w)[j]
        if u[0]:c[self.y[p,w]]+=Q(coef,u[0])
        else:c[self.a[p,w,j]]+=coef

    def shell_flows(self):
        # Prefix stationary edge balances for each exact pair of low-endpoint
        # types. Boundary donor from a third core is exposed as a *flux*;
        # constrain it by projectivity back to the vertex gain marginal.
        extra=defaultdict(list)
        flux={}
        for (u,v),i in self.y.items():
            for j,(r,d) in enumerate(((u,v),(v,u))):
                if not r[0]:continue
                slots=r[1]+1-1  # one actual included edge, irrespective of orientation
                if slots:
                    z=self.var(('outside gain',u,v,j),Q(slots,2));flux[u,v,j]=z
                    self.add(self.le,{z:1,i:-Q(slots,2)},0)
                    extra[j,r].append(z)
        # Sum over children of a fixed parent: k times total core input,
        # minus the included child inputs. Sum over the single parent of a
        # fixed child: total child input. These constraints share the same y.
        for j,r in product((0,1),self.types[:-1]):
            c=defaultdict(Q,{z:1 for z in extra[j,r]})
            for u,v in self.y:
                if u==r and j==0:self.donor(c,u,v,1,1-r[1])
                if v==r and j==0:self.donor(c,u,v,0,-r[1])
                if u==r and j==1:self.donor(c,u,v,1,-1)
            self.add(self.eq,c,0)
        eqs={key:defaultdict(Q,{i:-self.s}) for key,i in self.y.items()
             if key[0][0] and key[1][0]}
        for (u,v),i in self.y.items():
            key=u,v
            for j,r in enumerate(key):
                if not r[0]:continue
                for z,rate in ((self.gain(r),self.leaves(r)),(self.prom(r),self.bud(r))):
                    n=(z,v) if j==0 else (u,z)
                    if rate and key in eqs:eqs[key][i]-=rate
                    if rate and n in eqs:eqs[n][i]+=rate
                n=(self.gain(r),v) if j==0 else (u,self.gain(r))
                for k,sgn in ((key,-1),(n,1)):
                    if k in eqs:
                        self.donor(eqs[k],u,v,1-j,sgn)
                        if (u,v,j) in flux:eqs[k][flux[u,v,j]]+=sgn
        for u,r in self.birth.items():
            if r:eqs[self.prom(u),self.virgin][self.x[u]]+=r
        for c in eqs.values():self.add(self.eq,c,0)

    def age_cut(self):
        # The exact virgin balance already implies x_virgin >= s/(s+2).
        # A strictly positive uniform extra constraint: to reach D>N within
        # age t, the root must make N-1 gains. Dominate by a rate-D Yule.
        # Discounted successive holding times give the following upper.
        ub=math.prod(Q(d)/(d+self.s) for d in range(2,self.N+1))
        self.add(self.le,{self.x[self.tail]:1},ub)

    def objective(self,fp,gp,t,side):
        c=defaultdict(Q)
        for u,i in self.x.items():
            if not u[0]:continue
            plus,minus=fp[u]
            # Pay the finite row using its shared affine rho. Columns with
            # overflow neighbours carry their actual reciprocal moments.
            z0=plus[0]-t*minus[0]; z1=plus[1]-t*minus[1]
            c[i]+=side*z0
            coef=side*(z1-z0)/Q(u[1]+1,2)
            if coef:
                for p,w in self.y:
                    if u==p:self.donor(c,p,w,1,coef)
                    if u==w:self.donor(c,p,w,0,coef)
        for (u,v),i in self.y.items():
            if 0<u[0]<18 and 0<v[0]<18:
                p=(self.leaves(u),u[1],1); w=(self.leaves(v),v[1],1)
                a,b=gp[p,w];c[i]+=side*(a-t*b)
        return dict(c)


def optimize(m,c):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix
    def pack(rows):
        rr=[];cc=[];vv=[]
        for k,(r,_) in enumerate(rows):
            for i,x in r.items():rr.append(k);cc.append(i);vv.append(float(x))
        return csr_matrix((vv,(rr,cc)),shape=(len(rows),len(m.names))),[float(b) for _,b in rows]
    eq,b=pack(m.eq); le,d=pack(m.le)
    # Scale ppm contrasts before asking the floating LP for a candidate dual.
    sc=10**6
    f=np.array([float(sc*c.get(i,0)) for i in range(len(m.names))])
    r=linprog(f,A_ub=le,b_ub=d,A_eq=eq,b_eq=b,
              bounds=[(0,float(x)) for x in m.bounds],method='highs')
    if not r.success:return r,None
    # Weak duality holds for these rounded multipliers and every exact
    # feasible primal vector, regardless of LP quality. Pay reduced costs
    # that point outside the dual cone against the finite variable caps.
    rem=defaultdict(Q,c); center=Q(0); norm=Q(0)
    for rows,vals in ((m.eq,r.eqlin.marginals),(m.le,r.ineqlin.marginals)):
        for (row,rhs),y in zip(rows,vals):
            q=Q(round(float(y)*UNIT),UNIT*sc)
            if rows is m.le:q=min(q,Q(0))
            center+=q*rhs;norm+=abs(q)
            for i,x in row.items():rem[i]-=q*x
    repair=sum(min(0,x)*m.bounds[i] for i,x in rem.items())
    amps=[]
    for rows,vals in ((m.eq,r.eqlin.marginals),(m.le,r.ineqlin.marginals)):
        amps.extend(abs(Q(round(float(y)*UNIT),UNIT*sc)) for y in vals)
    return r,(center+repair,-repair,norm,max(amps,default=Q(0)))


def fixture(m):
    """Verify every generator-row coefficient against direct physical jumps.

    The first three old vertices form a directed test cycle solely to give
    everyone a parent. No cycle law is used in the LP or in the target process.
    """
    parents=(2,0,1,0,0,3,4,4,4,4,4,4,4,4,4)
    leaves=(2,1,0,0,1,2,1,2,1,3,1,1,1,1,1)
    def counts(pa,ll):
        k=[pa.count(i) for i in range(len(pa))]
        ds=[1+l+c for l,c in zip(ll,k)]
        us=[(d,c) if d<=m.N else m.tail for d,c in zip(ds,k)]
        v=defaultdict(Q)
        for u in us:v['vertex',u]+=1
        for i,j in enumerate(pa):
            u,w=us[j],us[i]
            v['edge',u,w]+=1
            for side,(r,q) in enumerate(((j,i),(i,j))):
                if us[r]==m.tail:v['reciprocal',u,w,side]+=Q(1,ds[r])
                else:
                    other=[b for a,b in enumerate(pa) if a==r and b!=q]
                    other += [a for a,b in enumerate(pa) if b==r and a!=q]
                    v['outside gain',u,w,side]+=sum(Q(1,ds[a]) for a in other)
        return ds,us,v
    ds,us,v=counts(parents,leaves)
    delta=defaultdict(Q)
    def jump(pa,ll,r):
        if r:
            _,_,w=counts(pa,ll)
            for key in w.keys()|v.keys():delta[key]+=r*(w[key]-v[key])
    for i,d in enumerate(ds):
        ll=list(leaves);ll[i]+=1
        receive=leaves[i]+Q(1,ds[parents[i]])+sum(Q(1,ds[j]) for j,b in enumerate(parents) if b==i)
        jump(parents,tuple(ll),receive)
        ll[i]-=2
        jump((*parents,i),(*ll,1),Q(leaves[i],d))
    ntypes=len(m.types)
    nproj=1+ntypes+sum(bool(d and k) for d,k in m.types)
    start_edge=nproj+3*(ntypes-1)
    keys=[('vertex',u) for u in m.types[:-1]]
    birth=sum(Q(l,d) for l,d in zip(leaves,ds))
    for (r,_),key in zip(m.eq[nproj:nproj+ntypes-1],keys,strict=True):
        calc=sum(x*v[m.names[j]] for j,x in r.items())
        assert calc==delta[key]-m.s*v[key]-(birth if key[1]==m.virgin else 0)
    for r,_ in m.eq[nproj+ntypes-1:start_edge]:
        assert sum(x*v[m.names[j]] for j,x in r.items())==0
    keys=[('edge',u,w) for u,w in m.y if u[0] and w[0]]
    for (r,_),key in zip(m.eq[start_edge:],keys,strict=True):
        assert sum(x*v[m.names[j]] for j,x in r.items())==delta[key]-m.s*v[key]


def run():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--caps',nargs='+',type=int,default=[6,8])
    p.add_argument('--s',type=Q,default=Q('.56572365'))
    p.add_argument('--t',type=Q,default=Q('1.24'))
    p.add_argument('--full-budget',action='store_true',help='for remainders also optimize individual paid floors')
    args=p.parse_args()
    pmin,pmax=Q('.999'),Q('1.12')
    odds=lambda s:(s-A)/(B-s)
    cross=lambda r:(A+B*r)/(1+r)
    print('full ratio if trial is critical (existing potential-ratio bounds):',
          float(odds(args.s)/pmax),float(odds(args.s)/pmin),flush=True)
    print('mu if full ratio equals t:',float(cross(args.t*pmin)),float(cross(args.t*pmax)),flush=True)
    dump=json.loads((HERE.parent/'rows'/'charged-f.json').read_text())
    with gzip.open(HERE.parent/'rows'/'pinned-fg-slabs.json.gz','rt') as f:pack=json.load(f)
    assert dump['payload_sha256']==pack['payload_sha256']
    f0,f1=(row_rewards(dump,n) for n in ('lower','upper'))
    g0,g1=(edge_rewards(pack[n]['g']) for n in ('lower','upper'))
    fp={k:(f0[k],f1[k]) for k in f0};gp={k:(g0[k],g1[k]) for k in g0}
    for n in args.caps:
        assert 2<=n<=64
        m=Model(args.s,n,True)
        fixture(m)
        print('D',n,'vars / eq / le',len(m.names),len(m.eq),len(m.le),flush=True)
        objectives=[('tail mass',{m.x[m.tail]:Q(1)}),('negative tail mass',{m.x[m.tail]:Q(-1)}),
                    ('contrast lower',m.objective(fp,gp,args.t,1)),
                    ('contrast negative lower',m.objective(fp,gp,args.t,-1))]
        if args.full_budget:
            zero=((Q(0),Q(0)),(Q(0),Q(0)))
            zf={k:(x,zero[1]) for k,(x,y) in fp.items()}
            zg={k:(x,Q(0)) for k,(x,y) in gp.items()}
            wf={k:(zero[0],y) for k,(x,y) in fp.items()}
            wg={k:(Q(0),y) for k,(x,y) in gp.items()}
            objectives += [('paid minus',m.objective(zf,zg,0,1)),('paid plus',m.objective(wf,wg,1,-1))]
        results={}
        for name,c in objectives:
            r,out=optimize(m,c)
            if out:
                v,repair,norm,maxmult=out
                results[name]=v
                print(name,'exact lower',math.floor(v*10**16)/10**16,'repair <=',
                      math.ceil(repair*10**18)/10**18,'dual L1',float(norm),
                      'max multiplier',float(maxmult),flush=True)
            else:print(name,'floating solver only:',r.message,flush=True)
        if args.full_budget and len(results)==len(objectives):
            # Conditional on this trial being the physical exponent, both
            # remainders (including higher families and all overflow) are >=0.
            # These whole-witness upper bounds use no discarded motif estimate.
            upper_minus=(args.s-A)*C0lo-results['paid minus']
            upper_plus=(B-args.s)*C0hi-results['paid plus']
            L=results['contrast lower']-args.t*upper_plus
            U=-results['contrast negative lower']+upper_minus
            print('unpaid minus/plus ceilings',math.ceil(upper_minus*10**16)/10**16,
                  math.ceil(upper_plus*10**16)/10**16,'all contrast outer',
                  math.floor(L*10**16)/10**16,math.ceil(U*10**16)/10**16,flush=True)


if __name__=='__main__':run()
