#!/usr/bin/env python3
"""Project only the unbooked f and g slack of the frozen two-ended corrector."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
OCC=HERE.parents[3]
ROOT=OCC.parents[1]
PINS={
 'round2/occupation-neon/interpolate/edge/replay.py':'e90651922acb6b3a427872cf3801f4b9097d2b83c7abe0ab75add6e60c78350f',
 'pressure-memory/marginal/replay.py':'b35c7e4faacd697fd85d5957f33fa27ab0790ff7c111a4cc578328940d684a70',
}
SPEC={2:('.565677930574',40,8),1:('.565770634792',32,7)}
Z=Q(0)


def module(path,name):
    spec=spec_from_file_location(name,ROOT/path);m=module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m)
    return m


def make_model(mem,limit,exact,memory):
    base_cell=mem.cells_for
    def cells_for(n):
        base=base_cell(n)
        def cell(d,k):return (d,d,k) if d<=exact else base(d,k)
        return cell
    mem.cells_for=cells_for
    try:orig=mem.model(limit,memory,2)
    finally:mem.cells_for=base_cell
    lp=mem.F.LP();lp.rhs=orig.rhs;cell=cells_for(limit)
    for kind,terms,x in orig.columns:
        if kind=='occ':
            (d,k),_,rho,_=x
            if d<=exact:
                terms=dict(terms)
                # Exact same receive rate as the f cost and census balance.
                terms['port-in',cell(d,k)]=rho
                if k:terms['port-out',cell(d,k)]=k*rho
        if kind!='edge':lp.col(kind,terms,x);continue
        a,b,p,q=x
        ad=a[1]<=exact;bd=b[1]<=exact
        U=(Z,Q(a[2],2)) if ad else (Z,)
        V=(Z,Q(b[2],2)) if bd else (Z,)
        for u,v in product(sorted(set(U)),sorted(set(V))):
            c=dict(terms)
            if ad:c['port-out',a]=-(u+q)
            if bd:
                c['port-in',b]=-(p+v)
                # For exact child census w=(L+p+v)/D. Substituting for
                # the old independently lifted child-w column keeps the
                # two receive consumers tied. Its memory source rescales
                # linearly by replacing w in the same q hull corner.
                new=q*(b[0]-b[2]-1+p+v)
                assert 0<=new<=1
                c['child_w',b]=-new
                # At exact q both Jensen bounds coincide, including when
                # the old independent w corner was zero.
                if memory and a[1]<=limit and a[2]>0:c['B-age',a]=-new*q/(1+q)
            lp.col(kind,c,(a,b,p,q,u,v))
    return lp


def freward(E,c,cs,s,V,z,tail,cut):
    delta=abs(s-c.c.lam);out={}
    for (d,k),(a,r) in cs.rows.items():
        l=d-k-1;b=Q(l,d);v=z[d,k];up=z[d+1,k] if d<64 else tail[k]
        bud=z.get((d,k+1),Z)+V.get(((d,k+1),E.J),Z)
        A=a+delta*cs.deficit[d,k]+l*(up-v)+b*(bud-v)-s*v
        B=r+up-v
        assert 0<=min(A,A+B*Q(k+1,2))<=(s+l+b+Q(k+1,2))/10**18
        if d<=cut:out[d,k]=A,B
    # In the original tail D>64 no edge can be newly retained (edge cap
    # <=40); same monotone support and tail ledger as the frozen witness.
    return out


def greward(E,c,edge,s,V,exact):
    get=lambda x,y:V.get((x,y),Z)
    out={}
    # exact <=17 restricts this price to the both-inner reserve; high-parent
    # reserve remains entirely untouched and paid by its existing ledger.
    assert exact<18
    for p in edge.p:
        if p[0]>exact:continue
        for w in edge.w:
            if w[0]>exact:continue
            a,(b,d)=E.ports(c,p,w);a-=E.INNER[c.c.mode]
            pp,pb=E.moves(p);wp,wb=E.moves(w)
            v=get(p,w);dv=get(pp,w)-v;dw=get(p,wp)-v
            # No newborn term here: the old edge survives the promotions.
            a+=(E.desc(p)[0]+Q(1,w[0]))*dv+(E.desc(w)[0]+Q(1,p[0]))*dw
            a+=Q(E.desc(p)[0],p[0])*(get(pb,w)-v)+Q(E.desc(w)[0],w[0])*(get(p,wb)-v)-s*v
            b+=dv;d+=dw
            assert min(a+b*u+d*v0 for u,v0 in product((Z,Q(p[1],2)),(Z,Q(w[1],2))))>=0
            out[p,w]=a,b,d
    return out


def objective(f,g):
    def pay(kind,x):
        if kind=='occ' and x[0] in f:
            a,b=f[x[0]];return 10**6*(a+b*x[2])
        if kind=='edge':
            p,w,_,_,u,v=x
            coeff=g.get(((p[0],p[2]),(w[0],w[2]))) if p[0]==p[1] and w[0]==w[1] else None
            if coeff is not None:
                a,b,d=coeff;return 10**6*(a+b*u+d*v)
        return Z
    return pay


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('payloads',nargs='+',type=Path)
    p.add_argument('--limit',type=int,default=18)
    p.add_argument('--exact',type=int,default=5)
    a=p.parse_args();assert __debug__ and 2<=a.exact<=8 and 18<=a.limit<=32
    for path,sha in PINS.items():assert sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    E,mem=(module(path,name) for path,name in zip(PINS,('frozen_edge','memory_projection'),strict=True))
    R=module(str((HERE.parent/'replay.py').relative_to(ROOT)),'frozen_residual')
    for name,sha in E.SOURCE_HASHES.items():assert sha256((ROOT/'round2'/name).read_bytes()).hexdigest()==sha,name
    off=make_model(mem,a.limit,a.exact,False);on=make_model(mem,a.limit,a.exact,True)
    assert off.rhs==on.rhs and len(off.columns)==len(on.columns)
    assert all((i,x,a)==(j,{r:v for r,v in y.items() if v and not (isinstance(r,tuple) and r[0]=='B-age')},b)
               for (i,x,a),(j,y,b) in zip(off.columns,on.columns,strict=True))
    print('finite shared-port models',[(len(m.rows),len(m.columns)) for m in (off,on)],flush=True)
    for path in a.payloads:
        c=E.Charges(E.D.Cert(path));mode=c.c.mode;s,n,k=SPEC[mode];s=Q(s)
        assert E.PAYLOAD_HASHES[mode]==digest(path)
        edge=E.Edge(c,n,k,'disjoint');cs=E.Census(c.c)
        V=edge.solve(s);m,_,z,tail=cs.solve(s,V)
        print('frozen mode',mode,'pointwise margin',m,'edge rows',len(V),flush=True)
        f=freward(E,c,cs,s,V,z,tail,a.limit);g=greward(E,c,edge,s,V,a.exact)
        for lab,lp in (('flow',off),('memory',on)):
            paid=[]
            for fam,tab1,tab2 in (('f',f,{}),('g',{},g),('f+g',f,g)):
                obj=objective(tab1,tab2)
                assert all(obj(j,x)>=0 for j,_,x in lp.columns)
                val,err=R.price(mem,lp,obj,(mode,lab,fam));paid.append(max(Z,val/10**6))
                print(mode,lab,fam,'lower',val/10**6,'display',mem.F.outward(val/10**6,False,16),'repair',err/10**6,flush=True)
            print(mode,lab,'joint excess over separate',paid[2]-sum(paid[:2]),flush=True)


def digest(path):
    from hashlib import file_digest
    with path.open('rb') as f:return file_digest(f,'sha256').hexdigest()


if __name__=='__main__':main()
