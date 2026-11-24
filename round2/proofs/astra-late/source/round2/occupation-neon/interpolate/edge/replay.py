#!/usr/bin/env python3
"""Exact two-ended edge Schur correction for the pinned Neon witnesses."""
import argparse
from collections import defaultdict
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import product
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent),str(HERE.parents[1]/'fringe')]
from corrector import D, operands
from charges import Charges, NB, Q, shape, up, bud, label

ZERO=Q(0)
J=(2,0)
RELEASE={1:Q('.565771634792'),2:Q('.565676930574')}
RETAINED={1:Q('1.501029458674e-6')+Q('4.663584315490e-7')+Q('2.799995635483e-7'),
          2:Q('1.276902849680e-6')+Q('5.896181971038e-7')+Q('5.520306761551e-7')}
DEN={1:Q('1.367169995210'),2:Q('1.356833021311')}
INNER={1:Q('5046841434259/340000000000000000000'),2:Q('25227936004759/1700000000000000000000')}
SOURCE_HASHES={
    'occupation-neon/fringe/charges.py':'44ed53fc7a0fc9ad25d941cd4160fe43aa9d00da334adb4e7de91a9a1a7c97e8',
    'occupation-neon/higher/extract.py':'9f40a50a3623756b277243d192f672e76c45b4b9784715ed4d0d8fb27d36131c',
    'attacks/neon1/audit.py':'d548a9aafa899ecab2b00a7f6ad0a4e6f828600088780f0a15f08927f17f0bfe',
    'occupation-neon/interpolate/corrector.py':'127a091d44996bd519963dc2720b37f76f93a9c9db4b6cd2f1b42d0d4c172f78',
    'occupation-neon/denom/check.py':'58f81c0c759ab64599ce81f05e0771cc82a82bdfbd1e98215e21ca584f0389c7',
    'occupation-neon/high-parent/replay.py':'497667552ff79080a460e5e99b0a4ab03f9e0e5120af341a6ec55973eb4e3d7e',
    'occupation-neon/interpolate/compose/lower.log':'a1a918a152ba0940ecd33d55e8042421d6eb76180310295179b38fae3af03c5d',
    'occupation-neon/interpolate/compose/upper.log':'054d18e6941a253735795bbe0561eee7f61022cc027bc29aacda695433b444c5',
}
PAYLOAD_HASHES={1:'1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41',
                2:'9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec'}


def pinned_rows(n,mode):
    path=HERE.parents[1]/'renewal'/'check.py'
    assert sha256(path.read_bytes()).hexdigest()=='423cfdddf34559c705b55b7f50f8c4a57f655915226d23cf7170028a0ff81388'
    spec=spec_from_file_location('pinned_neon_renewal',path)
    ref=module_from_spec(spec);spec.loader.exec_module(ref)
    data=(HERE.parents[1]/'rows'/'charged-f.json').read_bytes()
    assert sha256(data).hexdigest()=='415a712878e03b71c431a36bf4cf8eb11a047c2a47cafb5d9298fe23bc1561da'
    name={1:'upper',2:'lower'}[mode]
    frozen=ref.row_rewards(json.loads(data),name)
    return tuple(map(Q,ref.BANDS[name])),{x:r for x,r in frozen.items() if x[0]<=n}


def states(n,k0=0,kmax=None):
    return [(d,k) for d in range(n,1,-1) for k in range(min(d-1,d if kmax is None else kmax),k0-1,-1)]


def desc(x):
    d,k=x
    return d-k-1,k,1


def moves(x):
    d,k=x
    return (d+1,k),(d,k+1)


def ports(c,p,w):
    """Affine g in the two exterior reciprocal sums; neither is centered in A."""
    a,free=c.terms('g',(desc(p),desc(w)))
    coeff=[None,None]
    for i,_,_,b in free:
        if coeff[i] is not None:assert coeff[i]==b
        coeff[i]=b
    return a,tuple(b or ZERO for b in coeff)


class Edge:
    def __init__(self,c,n,arity,mask='full'):
        self.c=c;self.p=states(n,1,arity);self.w=states(n,0,arity)
        self.rows={}
        if mask=='disjoint':
            spec=spec_from_file_location('pinned_high_parent',HERE.parents[1]/'high-parent'/'replay.py')
            hp=module_from_spec(spec);spec.loader.exec_module(hp)
            physical=hp.Physical(c.c)
            booked={p:hp.charges(physical,*p)[0] for p in self.p if 18<=p[0]<=191 and p[1]<=10}
            physical.value.cache_clear()
        for p,w in product(self.p,self.w):
            a,coef=ports(c,p,w)
            pg=desc(p)[0]+Q(1,w[0]);wg=desc(w)[0]+Q(1,p[0]);pb=Q(desc(p)[0],p[0]);wb=Q(desc(w)[0],w[0])
            # Cancellation holds affine in u,v. Nonnegative true g also
            # permits using the hinge when the relaxed box is larger.
            box=(Q(p[1],2),Q(w[1],2))
            if mask=='disjoint':
                if p[0]<18 and w[0]<18:
                    a-=INNER[c.c.mode]
                elif p in booked and w==J:
                    assert a+min(ZERO,coef[0]*box[0])==booked[p]
                    a-=booked[p]
            pts=list(product((ZERO,box[0]),(ZERO,box[1])))
            if mask=='disjoint' and ((p[0]<18 and w[0]<18) or (18<=p[0]<=191 and p[1]<=10 and w==J)):
                assert min(a+coef[0]*u+coef[1]*v for u,v in pts)>=0
            if min(a+coef[0]*u+coef[1]*v for u,v in pts)<0:
                for axis in range(2):
                    b=coef[axis];j=1-axis
                    if not b:continue
                    for z in (ZERO,box[j]):
                        t=-(a+coef[j]*z)/b
                        if 0<t<box[axis]:pts.append((t,z) if axis==0 else (z,t))
            self.rows[p,w]=(pb,wb,tuple((pg+u,wg+v,max(ZERO,a+coef[0]*u+coef[1]*v)) for u,v in pts))

    def solve(self,s,mode='joint',limit=Q('0.004'),*,upper=False):
        values={}
        for p in self.p:
            pg,pp=moves(p)
            for w in self.w:
                wg,wp=moves(w);pb,wb,pts=self.rows[p,w]
                a=values.get((pg,w),ZERO);b=values.get((p,wg),ZERO)
                z=pb*values.get((pp,w),ZERO)+wb*values.get((p,wp),ZERO)
                if mode=='flat':
                    r=min(t for _,_,t in pts)
                    pts=tuple((g,h,r) for g,h,_ in pts)
                vals=tuple((r+g*a+h*b+z,s+pb+wb+g+h) for g,h,r in pts)
                best=min(x/y for x,y in vals)
                v=-D.down(-best) if upper else D.down(min(limit,best))
                assert v>=0
                if upper:assert any(v*y>=x for x,y in vals)
                else:assert all(v*y<=x for x,y in vals)
                values[p,w]=v
        return values


class Census:
    def __init__(self,c,n=64):
        self.c=c;self.n=n;self.rows=operands(c,n);self.charge=D.Charge(c)
        self.floors,pinned=pinned_rows(n,c.mode)
        assert {x:(a,a+b*Q(x[1]+1,2)) for x,(a,b) in self.rows.items()}==pinned
        self.t=Q(-9 if c.mode==1 else -2,1000)*c.S
        self.M=self.charge.cap(self.t)/c.S
        self.deficit={(d,k):self.M-(self.charge.exact(d-k-1,k,1)+self.t*(k-1))/c.S
                      for d,k in self.rows}
        assert min(self.deficit.values())>=0

    def solve(self,s,edges,*,floors=None,legacy=False,upper=False):
        c,n=self.c,self.n;delta=abs(s-c.lam);m=self.M
        floors=self.floors if floors is None else floors
        tail={};v={}
        rnd=(lambda x:-D.down(-x)) if upper else D.down
        for k in range(c.KF+1,-1,-1):
            deficit=m-(self.charge.by_k[1,k]+self.t*(k-1))/c.S
            assert deficit>=0
            a=(floors[k] if k<len(floors) else ZERO)+delta*deficit
            nxt=tail.get(k+1,ZERO)
            rates=(max(0,1-Q(k+1,n+1)),Q(1)) if k<=c.KF else (ZERO,)
            tail[k]=rnd(min((a+b*nxt)/(s+b) for b in rates))
            assert (any if upper else all)((s+b)*tail[k]>=a+b*nxt if upper else (s+b)*tail[k]<=a+b*nxt for b in rates)
        for (d,k) in states(n):
            l=d-k-1;b=Q(l,d);a,r=self.rows[d,k];a+=delta*self.deficit[d,k]
            u=v[d+1,k] if d<n else (tail.get(k,ZERO) if legacy else tail[min(k,c.KF+1)])
            source=edges.get(((d,k+1),J),ZERO)
            p=v.get((d,k+1),ZERO)+source
            vals=tuple((a+rho*r+(l+rho)*u+b*p,s+l+rho+b) for rho in (ZERO,Q(k+1,2)))
            x=rnd(min(z/h for z,h in vals))
            assert (any if upper else all)(h*x>=z if upper else h*x<=z for z,h in vals)
            v[d,k]=x
        norm=max(abs(x-v[J]) for x in (*v.values(),*tail.values()))+max(edges.values(),default=ZERO)+delta*m/s
        assert norm<Q('.01')
        return s*v[J]-delta*m,norm,v,tail


def fixtures(c,cs,s,V,result):
    """Full potential differences versus local incidence propagation."""
    from check import reach, tree, drift
    j=tree(1);bare=tree(0,j);old=tree(4)
    controls=list(reach(8))+[
        tree(2,tree(2,bare,old),bare),tree(1,tree(16,bare),tree(17,old)),
        tree(62,tree(63,bare,old),tree(65,j)),tree(1,tree(14001,old),old),
        tree(2,tree(1,*([j]*(c.c.KF+2))))]
    margin,_,z,tail=result;delta=abs(s-c.c.lam)
    def val(p,x):
        def census(y):return sum(y),y[1]
        q=[census(y) for y in x]
        ans=-delta*cs.M/s
        for i in range(1,len(x)):
            d,k=q[i];ans+=z.get((d,k),tail[min(k,c.c.KF+1)])-z[J]
            if p[i]>0:ans+=V.get((q[p[i]],q[i]),ZERO)
        return ans
    births=0;minslack=None
    for tr in controls:
        p,x=shape(tr,0);ch=[[j for j,a in enumerate(p) if a==i] for i in range(len(p))]
        base=val(p,x);corr=-s*base;local=-s*base
        for i,y in enumerate(x):
            l,k,h=y;d=sum(y);b=Q(l,d)
            g=l+sum(Q(1,sum(x[a])) for a in ch[i]+([p[i]] if i else []))
            corr+=g*(val(p,label(x,i,up(y)))-base)
            if b:corr+=b*(val(p+(i,),label(x,i,bud(y))+(NB,))-base)
            get=lambda y:(sum(y),y[1])
            if i:
                def zv(t):
                    dd,kk=get(t);return z.get((dd,kk),tail[min(kk,c.c.KF+1)])
                local+=g*(zv(up(y))-zv(y))
                if b:
                    local+=b*(zv(bud(y))-zv(y)+V.get((get(bud(y)),J),ZERO))
                    births+=V.get((get(bud(y)),J),ZERO)>0
            for a in ch[i]+([p[i]] if i else []):
                parent,w=(i,a) if p[a]==i else (a,i)
                if parent==0:continue
                def ev(t):return V.get((get(t),get(x[a])) if i==parent else (get(x[a]),get(t)),ZERO)
                local+=g*(ev(up(y))-ev(y))
                if b:local+=b*(ev(bud(y))-ev(y))
        assert corr==local,(p,x,corr,local)
        physical=c.sign*drift(c.c,tr)-delta*c.c.potential(ch,x)/c.c.S+corr
        slack=physical-margin*(len(x)-1)
        assert slack>=0,(p,x,float(c.sign*drift(c.c,tr)),float(corr),float(physical),float(slack))
        minslack=slack if minslack is None else min(minslack,slack)
    assert births>0
    print(c.sign,len(controls),'full-tree residual checks; nonzero edge births',births,'minimum residual surplus',float(minslack),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('blob',nargs='+',type=Path)
    ap.add_argument('--degree',type=int,default=18);ap.add_argument('--arity',type=int,default=8)
    ap.add_argument('--s',type=Q);ap.add_argument('--ticks',type=int,default=1000)
    ap.add_argument('--flat',action='store_true')
    ap.add_argument('--challenge',action='store_true',help='test release and 1e-6 targets with overlap-safe and optimistic accounting')
    ap.add_argument('--mask',choices=('full','disjoint'),default='full')
    a=ap.parse_args();assert 4<=a.degree<=64 and 1<=a.arity<a.degree and not sys.flags.optimize
    for name,h in SOURCE_HASHES.items():assert sha256((HERE.parents[2]/name).read_bytes()).hexdigest()==h,name
    for path in a.blob:
        c=Charges(D.Cert(path));assert c.c.D==18 and c.c.lam==(Q('.565673') if c.sign==1 else Q('.565775'))
        with path.open('rb') as blob:
            h=sha256()
            for chunk in iter(lambda:blob.read(2**23),b''):h.update(chunk)
        assert h.hexdigest()==PAYLOAD_HASHES[c.c.mode]
        print(c.sign,'authenticated payload',h.hexdigest(),flush=True)
        cs=Census(c.c);e=Edge(c,a.degree,a.arity,a.mask);mode='flat' if a.flat else 'joint'
        if a.challenge:
            side=c.c.mode
            for s in ((a.s,) if a.s else (RELEASE[side],RELEASE[side]+c.sign*Q('1e-6'))):
                v=e.solve(s,mode);r=cs.solve(s,v);f=cs.solve(s,{})[0]
                high=e.solve(s,'joint',upper=True);ub=cs.solve(s,high,upper=True)[0]
                assert r[0]<=ub
                print(c.sign,a.mask,len(v),'trial',s,'baseline',float(f),'margin',float(r[0]),
                      'robust class upper',float(ub),'norm',float(r[1]),'edge norm',float(max(v.values())),
                      'pointwise shortfall',float(max(ZERO,-ub)),flush=True)
                fixtures(c,cs,s,v,r)
                if a.mask=='disjoint':
                    z,tail=r[2:];H=max(c.sign*(w-z[J]) for w in (*z.values(),*tail.values()))
                    if c.sign==1:H+=max(v.values())
                    paid=r[0]+RETAINED[side]
                    if paid>=0:
                        endpoint=s+c.sign*paid/(DEN[side]+H)
                        print('disjoint composition outward',D.show(endpoint,15,upper=c.sign<0),
                              'gain against release',float(c.sign*(endpoint-RELEASE[side])),flush=True)
                # Upper ceiling on everything the frozen ledger books,
                # used only to challenge this local proof even if overlap
                # were waived. It is not a bound on actual physical slack.
                if s==RELEASE[side]+c.sign*Q('1e-6'):
                    entitlement_ceiling=Q('2.75e-6') if side==1 else Q('2.42e-6')
                    print('shortfall after granting entire frozen numerator ceiling',float(-ub-entitlement_ceiling),flush=True)
            continue
        # At s0 greater than any trial, nonnegative edge values remain
        # subsolutions for all candidate discounts (s-A)V <= (s0-A)V.
        edge=e.solve(Q('.565775'),mode)
        def search(table):
            low,high=0,100_000_000//a.ticks
            assert high>1
            while high-low>1:
                mid=(high+low)//2;s=c.c.lam+c.sign*Q(mid*a.ticks,10**12)
                if cs.solve(s,table)[0]>0:low=mid
                else:high=mid
            return c.c.lam+c.sign*Q(low*a.ticks,10**12)
        if a.s is None:
            s=search(edge);base=search({})
            print(c.sign,'grid baseline',D.show(base),'grid corrected',D.show(s),'inward',float(c.sign*(s-base)),flush=True)
        else:s=a.s
        assert 0<s<=Q('.565775') and c.sign*(s-c.c.lam)>=0
        result=cs.solve(s,edge);baseline=cs.solve(s,{})[0]
        margin,norm,*_=result
        print(c.sign,mode,'edge rows',len(edge),'s',s,'margin',float(margin),'baseline',float(baseline),
              'norm',float(norm),'edge norm',float(max(edge.values())),flush=True)
        if margin>0:fixtures(c,cs,s,edge,result)


if __name__=='__main__':main()
