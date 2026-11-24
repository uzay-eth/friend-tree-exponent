#!/usr/bin/env python3
"""Exact post-EDGE residual, physical column lifts and dual repair audit."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import file_digest, sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
O=ROOT/'round2/occupation-neon'; R=O/'interpolate/compose/residual'
EDGE=O/'interpolate/edge'; MEM=ROOT/'pressure-memory/marginal/replay.py'
PROJ=R/'edge-projection/replay.py'; LIMIT=18; EXACT=4; SAT=12; SCALE=10**6
Z=Q(0); J=(2,0)
PINS={PROJ:'65210a3f2d5fa4ae5631b2f57b98e1e6fed44b24b2d73b080288c4c1706748e9',
      MEM:'b35c7e4faacd697fd85d5957f33fa27ab0790ff7c111a4cc578328940d684a70',
      O/'flow/check.py':'28365603eafd33eb07901c308cadcc5ac3e0f6cfc39d883a02dbedbda18f4e35'}
BASE={2:Q('191976090122981726641370633087/339373453939523193250000000000'),
      1:Q('24172191945911554486185483421/42724375658564528593750000000')}


def load(name,path):
    spec=spec_from_file_location(name,path);m=module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def digest(path):
    with path.open('rb') as f:return file_digest(f,'sha256').hexdigest()


def mixes(x,a,b):
    assert a<=x<=b,(x,a,b)
    if a==b:return [(a,Q(1))]
    t=Q(x-a,b-a)
    return [(v,w) for v,w in ((a,1-t),(b,t)) if w]


def aggregate(col,rows,scale):
    for r,x in rows.items():col[r]+=scale*x


def clean(rows):return {r:x for r,x in rows.items() if x}


def physical_check(mem,lp):
    """Rebuild an actual tree's local contributions and find convex column lifts."""
    index=defaultdict(list)
    for kind,col,m in lp.columns:
        key=m[0] if kind=='occ' else m[:2] if kind=='edge' else None
        index[kind,key].append(col)
    # Multiple asymmetric exact bins; both unsaturated and saturated tails;
    # huge ordinary degree and core-degree>=SAT with nonzero incoming traffic.
    p=[-1,0,1,2,2,3,1,6,0,8]+[1]*14+[8]*14
    k=[p.count(i) for i in range(len(p))]
    d=[2,10**10,4,3,2,4,18,23,19,2]+[2,3]*14
    d[0]=3
    def cell(d,k):return (d,d,k) if d<=EXACT else base(min(d,LIMIT+1),k)
    base=mem.cells_for(LIMIT); C=[cell(a,b) for a,b in zip(d,k)]
    assert all(D>=max(2,K+(i>0)) for i,(D,K) in enumerate(zip(d,k)))
    q=[Q(1,D) for D in d];kids=[[v for v,j in enumerate(p) if j==i] for i in range(len(p))]
    up=[sum((q[j] for j in vs),Z) for vs in kids];rho=[up[i]+(q[p[i]] if i else 0) for i in range(len(p))]
    L=[D-K-(i>0) for i,(D,K) in enumerate(zip(d,k))];W=[(l+r)*a for l,r,a in zip(L,rho,q)]
    mu=(mem.F.LO+mem.F.HI)/2
    for i in range(1,len(p)):
        D,K,c,pp=d[i],k[i],C[i],q[p[i]]; b=L[i]*q[i];w=W[i];r=rho[i]
        x=(D,K) if D<=LIMIT else (LIMIT+1,min(K,SAT)); sat=x==(LIMIT+1,SAT)
        B=Q(K,2)-up[i]
        want=defaultdict(Q,{'mass':1,'k':K,'birth':b-mu,'edges-rho':r+mu-1,
                            ('child_n',c):1,('child_q',c):q[i],('parent_n',c):K,
                            ('parent_q',c):K*q[i],('down',c):pp,('up',c):up[i],('child_w',c):w})
        want['age',x]+=mu;want['age',J]-=mu
        if D<=LIMIT:
            g=L[i]+r;want['age',x]+=g+b
            want['age',(D+1,min(K,SAT) if D==LIMIT else K)]-=g
            if b:want['age',(D,K+1)]-=b
            if K:
                want['B-age',c]+=(mu+g+b)*B
                if D<LIMIT:want['B-age',cell(D+1,K)]-=g*B
                if b:want['B-age',cell(D,K+1)]-=b*B
            if D<=EXACT:want['port-in',c]+=r;want['port-out',c]+=K*r
            A=Q(K,2);lo,hi=(0,A/2) if B<=A/2 else (A/2,A)
            if hi==lo:bc=[((Z,Z),Q(1))]
            else:
                t=Q(B-lo,hi-lo);bc=[((lo,lo**2),(1-t)**2),(((lo+hi)/2,lo*hi),2*t*(1-t)),((hi,hi**2),t*t)]
            phys=[(b,q[i],Q(1))]
        else:
            bc=[((Z,Z),Q(1))]
            if not sat:
                want['age',x]+=b;want['age',(x[0],K+1)]-=b
                phys=[(1-(K+1)*v,v,t) for v,t in mixes(q[i],0,Q(1,LIMIT+1))]
            else:
                cap=Q(1,LIMIT+1);t=b/(1-(SAT+1)*q[i]);physical=[]
                for v,a in mixes(q[i],0,cap):physical.extend([(Z,v,a*(1-t)),(1-(SAT+1)*v,v,a*t)])
                phys=physical
        got=defaultdict(Q)
        _,j0=x
        # Build the index order from the explicitly derived coordinates.
        phys_order=([(b,q[i])] if D<=LIMIT else [(1-(j0+1)*v,v) for v in (Z,Q(1,LIMIT+1))]
                    if not sat else [(Z,Z),(Q(1),Z),(Z,Q(1,LIMIT+1)),(1-Q(SAT+1,LIMIT+1),Q(1,LIMIT+1))])
        cap=Q(j0,2) if not sat else Z
        B_order=mem.controls(Z,cap,2,D<=LIMIT)
        cs=index['occ',x]
        for m,a in mixes(mu,mem.F.LO,mem.F.HI):
            for pp0,h in mixes(pp,0,Q(1,2)):
                for b0,p0,t in phys:
                    for B0,j in bc:
                        for w0,f in ([(w,Q(1))] if D<=LIMIT else [(v,a) for v,a in mixes(w,b,(1+b)/2)]):
                            if not a*h*t*j*f:continue
                            row=(((int(m==mem.F.HI)*2+int(pp0!=0))*len(phys_order)+phys_order.index((b0,p0)))*len(B_order)+B_order.index(B0))
                            ix=row if D<=LIMIT else 2*row+int(w0!=b)
                            aggregate(got,cs[ix],a*h*t*j*f)
        if sat:
            # Physical saturated pp still belongs to the down marginal;
            # only upward receive is put in the extra-budget column.
            aggregate(got,index['tailup',None][0],up[i]);aggregate(got,index['stub',None][0],K-SAT)
            # Reconstruct pp from the actual corresponding column difference.
            pp_row=len(cs)//4
            diff=defaultdict(Q);aggregate(diff,cs[pp_row],1);aggregate(diff,cs[0],-1)
            assert clean(diff)=={'edges-rho':Q(1,2),('down',c):Q(1,2)}
        assert clean(got)==clean(want),(i,clean(got),clean(want))

    for i in range(2,len(p)):
        j=p[i]
        if not j:continue
        a,c=C[j],C[i];u=rho[j]-q[i];v=rho[i]-q[j]
        want=defaultdict(Q,{'edges':1,('parent_n',a):-1,('parent_q',a):-q[j],
                             ('child_n',c):-1,('child_q',c):-q[i],('down',c):-q[j],
                             ('up',a):-q[i],('child_w',c):-W[i]})
        if d[j]<=EXACT:want['port-out',a]=-rho[j]
        if d[i]<=EXACT:want['port-in',c]=-rho[i]
        if d[j]<=LIMIT:want['B-age',a]=-W[i]*q[i]/(1+q[i])
        cs=index['edge',(a,c)];got=defaultdict(Q)
        ap=sorted(set(mem.rrange(a,LIMIT)));cq=sorted(set(mem.rrange(c,LIMIT)))
        kmin=c[2];kmax=0 if not kmin else 3 if kmin==1 else None
        ws=sorted({max(Z,1-(kmax+1)*cq[-1]) if kmax is not None else Z,1-Q(kmin+1,2)*cq[0]})
        us=sorted({Z,Q(k[j],2)}) if d[j]<=EXACT else [Z]
        vs=sorted({Z,Q(k[i],2)}) if d[i]<=EXACT else [Z]
        for x,y in product(mixes(q[j],ap[0],ap[-1]),mixes(q[i],cq[0],cq[-1])):
            pp,h=x;rr,t=y
            ss=sorted(set(mem.source(rr,cq[0],cq[-1]))) if d[j]<=LIMIT else [Z]
            # Chord and tangent are affine globally in q. Use a single
            # interpolation weight at the physical q before corner lifting.
            bounds=mem.source(q[i],cq[0],cq[-1]) if d[j]<=LIMIT else (Z,Z)
            vm=q[i]/(1+q[i]) if d[j]<=LIMIT else Z
            theta=Q(vm-bounds[0],bounds[1]-bounds[0]) if bounds[0]!=bounds[1] else Z
            for w,f in ([(ws[0],Q(1))] if d[i]<=EXACT else mixes(W[i],ws[0],ws[-1])):
                for si,e in ([(0,Q(1))] if len(ss)==1 else [(0,1-theta),(1,theta)]):
                    for uu,g in (mixes(u,us[0],us[-1]) if d[j]<=EXACT else [(Z,Q(1))]):
                        for vv,b in (mixes(v,vs[0],vs[-1]) if d[i]<=EXACT else [(Z,Q(1))]):
                            offset=(sum(len(set(mem.source(q0,cq[0],cq[-1]))) if d[j]<=LIMIT else 1 for q0 in cq)*ap.index(pp)
                                    +sum(len(set(mem.source(q0,cq[0],cq[-1]))) if d[j]<=LIMIT else 1 for q0 in cq[:cq.index(rr)]))*len(ws)
                            ix=(((offset+ws.index(w)*len(ss)+si)*len(us)+us.index(uu))*len(vs)+vs.index(vv))
                            aggregate(got,cs[ix],h*t*f*e*g*b)
        assert clean(got)==clean(want),(j,i,clean(got),clean(want))
    print('physical column lifts:',len(p)-1,'vertex and',sum(p[i]>0 for i in range(1,len(p))),'edge controls',flush=True)


def residual(E,c,s,V,cs,z,tail,oracle):
    def ez(D,k):return z[D,k] if D<=64 else tail[min(k,c.c.KF+1)]
    def edge(x,y):return V.get((x,y),Z)
    def slopes(fun,n):
        a=fun(*([Z]*n));b=[fun(*(Q(i==j) for i in range(n)))-a for j in range(n)]
        return a,*b
    eps=c.sign;delta=eps*(s-c.c.lam);f={};g={};missing=double=unpaid=0
    for x in z:
        D,k=x;L=D-k-1;bud=(D,k+1);b=Q(L,D)
        a,r=cs.rows[x];new=b*edge(bud,J);v=z[x]
        def rf(rho):
            return (a+r*rho+delta*cs.deficit[x]+(L+rho)*(ez(D+1,k)-v)
                    +(b*(z[bud]-v) if b else 0)+new-s*v)
        checks=[rf(y) for y in (Z,Q(k+1,2))];assert min(checks)>=0
        if D<=LIMIT:
            f[x]=slopes(rf,1);missing+=min(checks)<new
    for k,v in tail.items():
        a=(cs.floors[k] if k<len(cs.floors) else Z)+delta*(cs.M-(cs.charge.by_k[1,k]+cs.t*(k-1))/c.c.S)
        step=tail.get(k+1,v)-v
        rates=(Z,) if k>c.c.KF else (max(Z,1-Q(k+1,65)),Q(1))
        assert all(a+b*step>=s*v for b in rates)
    for x,y in V:
        D,k=x;H,j=y
        if D>EXACT or H>EXACT:continue
        L=D-k-1;M=H-j-1;a,ss,_,_=oracle.row(c.c,(Q(D),k,1),(Q(H),j,1));v=V[x,y]
        assert ss[0]==ss[1]
        def rg(u,w):
            gain=(L+Q(1,H)+u)*(edge((D+1,k),y)-v)+(M+Q(1,D)+w)*(edge(x,(H+1,j))-v)
            old=Q(L,D)*(edge((D,k+1),y)-v)+Q(M,H)*(edge(x,(H,j+1))-v)
            return eps*(a+ss[0]*u+(ss[2]*w if j else 0))-E.INNER[c.c.mode]+gain+old-s*v
        assert all(rg(u,w)>=0 for u,w in product((Z,Q(k,2)),(Z,Q(j,2))))
        g[x,y]=slopes(rg,2);unpaid+=bool(E.INNER[c.c.mode]);double+=bool(L*edge((D,k+1),J)+M*edge((H,j+1),J))
    assert len(f)==170 and len(g)==54 and missing and double and unpaid==len(g)
    print('residual paid supports:',len(f),len(g),'birth omission fails',missing,
          'old-edge double-birth controls',double,'same-edge reserve controls',unpaid,flush=True)
    return f,g


def dual(lp,objective):
    import scipy.optimize as so
    from scipy.sparse import coo_matrix
    ii,jj,vv=[],[],[];rhs=[Z]*len(lp.rows);cost=[]
    for r,v in lp.rhs.items():rhs[lp.rows[r]]=Q(v)
    for j,(kind,col,meta) in enumerate(lp.columns):
        cost.append(objective(kind,meta));assert cost[-1]>=0
        for r,x in col.items():ii.append(lp.rows[r]);jj.append(j);vv.append(float(x))
    A=coo_matrix((vv,(ii,jj)),shape=(len(rhs),len(cost))).tocsc()
    best=None
    for method,tol,S in (('highs-ipm',1e-9,1),('highs-ipm',1e-7,10),('highs-ds',1e-7,1)):
        opt=so.linprog([float(c*S) for c in cost],A_eq=A,b_eq=list(map(float,rhs)),method=method,
                      options={'primal_feasibility_tolerance':tol,'dual_feasibility_tolerance':tol})
        if not opt.success:continue
        grid=10**11;y={r:Q(round(opt.eqlin.marginals[j]*grid),grid*S) for r,j in lp.rows.items()}
        errors=[sum((x*y[r] for r,x in col.items()),Z)-c for (_,col,_),c in zip(lp.columns,cost,strict=True)]
        cap=defaultdict(Q)
        for (kind,_,_),e in zip(lp.columns,errors,strict=True):cap[kind]=max(cap[kind],e)
        assert set(cap)=={'occ','edge','stub','tailup'}
        assert all(e<=cap[kind] for (kind,_,_),e in zip(lp.columns,errors,strict=True))
        result=(sum((y[r]*x for r,x in lp.rhs.items()),Z)-sum(cap.values()))/SCALE
        if best is None or result>best[0]:best=(result,cap)
        if result>0:break
    assert best is not None
    return best


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('payloads',type=Path,nargs='*');a=ap.parse_args()
    assert __debug__
    for p,h in PINS.items():assert digest(p)==h,p
    for line in (EDGE/'inputs.sha256').read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            h,p=line.split();assert digest(ROOT/p)==h,p
    P=load('projection',PROJ);E=load('frozen_edge',EDGE/'replay.py');mem=load('physical_memory',MEM)
    oracle=load('independent_physical_g',O/'all-contrast/check/row_oracle.py')
    assert (mem.F.LO,mem.F.HI)==(Q('.565673570322'),Q('.565775'))
    off=P.make_model(mem,LIMIT,EXACT,False);on=P.make_model(mem,LIMIT,EXACT,True)
    assert off.rhs==on.rhs=={'mass':1,'edges':1,'k':1}
    counts=Counter(j for j,_,_ in on.columns);assert len(on.columns)==9674 and len(on.rows)-len(off.rows)==20
    for (i,c,x),(j,d,y) in zip(off.columns,on.columns,strict=True):
        assert (i,c,x)==(j,{r:v for r,v in d.items() if not (isinstance(r,tuple) and r[0]=='B-age')},y)
    print('matched models',len(off.rows),len(on.rows),'columns',dict(counts),flush=True)
    physical_check(mem,on)
    for path in a.payloads:
        c=E.Charges(E.D.Cert(path));md=c.c.mode;assert digest(path)==E.PAYLOAD_HASHES[md]
        ss,d,k=P.SPEC[md];s=Q(ss);e=E.Edge(c,d,k,'disjoint');cs=E.Census(c.c);V=e.solve(s);m,_,z,tail=cs.solve(s,V)
        f,g=residual(E,c,s,V,cs,z,tail,oracle)
        assert f==P.freward(E,c,cs,s,V,z,tail,LIMIT)
        assert g==P.greward(E,c,e,s,V,EXACT)
        def obj(kind,x):
            if kind=='occ' and x[0] in f:
                a,b=f[x[0]];return SCALE*(a+b*x[2])
            if kind=='edge':
                p,w,_,_,u,v=x
                if p[1]<=EXACT and w[1]<=EXACT:
                    a,b,d=g[(p[0],p[2]),(w[0],w[2])];return SCALE*(a+b*u+d*v)
            return Z
        ref=P.objective(f,g)
        assert all(obj(j,x)==ref(j,x) for j,_,x in on.columns)
        obtained=[]
        for label,lp in (('flow',off),('memory',on)):
            val,cap=dual(lp,obj);obtained.append(val)
            print(md,label,'exact joint floor',val,'repair numerator',sum(cap.values())/SCALE,flush=True)
        lone,cap=dual(on,lambda j,x:obj(j,x) if j=='occ' else Z)
        assert lone>0
        print(md,'memory omit same-edge objective',lone,'repair numerator',sum(cap.values())/SCALE,flush=True)
        paid=max(obtained);net=m+E.RETAINED[md];eps=c.sign
        H=max(eps*(v-z[J]) for v in (*z.values(),*tail.values()))+(max(V.values()) if eps==1 else Z)
        den=E.DEN[md]+H;old=s+eps*net/den
        assert old==BASE[md] and paid>0 and net>0
        end=s+eps*(net+paid)/den
        assert eps*(end-old)==paid/den
        print(md,'den',den,'exact endpoint',end,'outward',E.D.show(end,15,upper=eps<0),flush=True)


if __name__=='__main__':main()
