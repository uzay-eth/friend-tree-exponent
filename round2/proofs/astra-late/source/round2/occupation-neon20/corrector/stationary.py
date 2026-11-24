"""D20 stationary census/edge incidence relaxation with exact dual repair."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product

from extract import LO, HI


class Stationary:
    def __init__(self, limit):
        import numpy as np
        from scipy.sparse import coo_matrix
        self.names={};self.cols=[];self.meta=[]
        self.rhs={'mass':1,'arity':1,'edges':1}
        cap=8
        end=limit+1
        ends=sorted({*range(2,9),11,15,19,24,30,40,limit,end}&set(range(2,end+1)))
        bands={}
        for d in range(2,end+1):
            i=next(i for i,h in enumerate(ends) if h>=d)
            bands[d]=(2 if i==0 else ends[i-1]+1,ends[i])
        def cell(d,k):return (*bands[d],min(k,3))
        def add(kind,m,terms):
            self.cols.append(tuple((self.names.setdefault(key,len(self.names)),Q(a)) for key,a in terms.items() if a))
            self.meta.append((kind,m))
        nodes=set(); parents=set()
        for d in range(2,end+1):
            for k in range(d if d<end else cap+1):
                x=(d,k);c=cell(d,k);nodes.add(c)
                if k:parents.add(c)
                sat=x==(end,cap)
                if d<end:
                    bp=((Q(d-k-1,d),Q(1,d)),)
                elif not sat:
                    bp=((Q(1),Q(0)),(1-Q(k+1,d),Q(1,d)))
                else:
                    bp=((0,0),(1,0),(0,Q(1,d)),(1-Q(cap+1,d),Q(1,d)))
                for mu,(b,p),u,v in product((LO,HI),bp,(0,Q(1,2)),(0,Q(k,2)) if k and not sat else (0,)):
                    z=defaultdict(Q,{'mass':1,'arity':k,'birth':b-mu,'receive':u+v+mu-1,
                                    ('child-n',c):1,('child-p',c):p,('parent-n',c):k,
                                    ('parent-p',c):1-b-p,('up',c):v,('down',c):u})
                    z['age',x]+=mu+b;z['age',(2,0)]-=mu
                    if b:z['age',(d,k+1) if not sat else x]-=b
                    if d<end:
                        g=d-k-1+u+v
                        z['age',x]+=g;z['age',(d+1,k if d+1<end else min(k,cap))]-=g
                    add('occ',(x,b,u+v),z)
        c=cell(end,cap)
        add('stubs',None,{'arity':1,('parent-n',c):1})
        add('tailup',None,{'receive':1,('up',c):1})
        def reciprocal(c):
            return {Q(1,c[0]),Q(1,c[1]) if c[1]<end else Q(0)}
        for a,b in product(sorted(parents),sorted(nodes)):
            for p,q in product(sorted(reciprocal(a)),sorted(reciprocal(b))):
                add('edge',None,{'edges':1,('parent-n',a):-1,('parent-p',a):-p,
                                 ('child-n',b):-1,('child-p',b):-q,('up',a):-q,('down',b):-p})
        rr,cc,vv=[],[],[]
        for j,z in enumerate(self.cols):
            for i,v in z:rr.append(i);cc.append(j);vv.append(float(v))
        self.A=coo_matrix((vv,(rr,cc)),shape=(len(self.names),len(self.cols))).tocsc()
        self.b=np.zeros(len(self.names))
        for r,v in self.rhs.items():self.b[self.names[r]]=v
        self.limit=limit
        self.cap=cap

    def bound(self,cs,E,delta):
        from scipy.optimize import linprog
        scale=10**7
        tail={}
        for k in range(self.cap+1):
            ops=[(z+min(0,w*Q(k1+1,2)),cs.deficit[d,k1])
                 for (d,k1),(z,w) in cs.pay.items() if d>self.limit and min(k1,self.cap)==k]
            ops.extend((cs.bands[j] if j<len(cs.bands) else Q(0),cs.tail[j])
                       for j in range(k,cs.K+2) if min(j,self.cap)==k)
            tail[k]=min(a+delta*v for a,v in ops)
        def obj(kind,m):
            if kind!='occ':return Q(0)
            (d,k),b,r=m
            if d<=self.limit:
                z,w=cs.pay[d,k]; a=z+w*r+b*E.get(((d,k+1),(2,0)),0)
                v=cs.deficit[d,k]
            else:
                return scale*tail[k]
            return scale*(a+delta*v)
        costs=[obj(*m) for m in self.meta]
        res=linprog(list(map(float,costs)),A_eq=self.A,b_eq=self.b,method='highs-ipm',
                    options={'ipm_optimality_tolerance':1e-11})
        assert res.success,res.message
        grid=10**12;y=[Q(round(f*grid),grid) for f in res.eqlin.marginals]
        defects=defaultdict(Q)
        paired=[sum(y[i]*a for i,a in t) for t in self.cols]
        for (kind,_),v,c in zip(self.meta,paired,costs):
            defects[kind]=max(defects[kind],v-c)
        assert all(v<=c+defects[kind] for (kind,_),v,c in zip(self.meta,paired,costs))
        repair=sum(defects.values(),Q(0))/scale
        return sum(Q(v)*y[self.names[r]] for r,v in self.rhs.items())/scale-repair,repair
