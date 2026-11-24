#!/usr/bin/env python3
"""Authenticate the frozen corrector and price only its unused f residual."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
OCC=HERE.parents[2]
ROOT=OCC.parents[1]


def load(name,path):
    spec=spec_from_file_location(name,path);mod=module_from_spec(spec)
    sys.modules[name]=mod;spec.loader.exec_module(mod)
    return mod


def residual(c,charge,rows,floors,s,recharge,z,tail,M,cut):
    delta=abs(s-c.lam);assert 2<=cut<=64
    assert (s-c.lam)*(1 if c.mode==2 else -1)==delta
    pay={}
    for (d,k),(a,r) in rows.items():
        l=d-k-1;b=Q(l,d);v=z[d,k];y=z[d+1,k] if d<64 else tail[k]
        w=z[d,k+1] if b else 0
        deficit=M-charge.exact(l,k,1)/c.S-recharge*(k-1)
        base=a+delta*deficit+l*(y-v)+b*(w-v)-s*v
        coef=r+y-v
        worst=min(base,base+Q(k+1,2)*coef)
        assert worst>=0 and worst<=(s+l+Q(k+1,2)+b)/10**18
        if d<=cut:pay[d,k]=base,coef
    for k,v in tail.items():
        a=(floors[k] if k<len(floors) else Q(0))+delta*(M-Q(charge.by_k[1,k],c.S)-recharge*(k-1))
        dv=tail[k+1]-v if k<c.KF+1 else 0
        assert min(a+b*dv-s*v for b in (max(0,1-Q(k+1,65)),Q(1)))>=0
    assert len(pay)==sum(range(2,cut+1))
    return pay


def cost(coeff):
    def obj(kind,meta):
        if kind!='occ' or meta[0] not in coeff:return 0
        a,b=coeff[meta[0]];return 10**6*(a+b*meta[2])
    return obj


def price(mem,lp,reward,key):
    import scipy.optimize
    func=scipy.optimize.linprog
    for method,tol,scale in (('highs-ipm',1e-9,1),('highs-ipm',1e-7,10),('highs-ds',1e-7,1)):
        def alt(*args,**opts):
            opts['method']=method
            opts['options']={'primal_feasibility_tolerance':tol,'dual_feasibility_tolerance':tol}
            return func(*args,**opts)
        scipy.optimize.linprog=alt
        try:
            val,err=mem.optimize(lp,lambda k,m:scale*reward(k,m),key)
            return val/scale,err/scale
        except AssertionError as err:
            if 'HiGHS Status' not in str(err):raise
            print(key,'proposal failed',method,tol,flush=True)
        finally:scipy.optimize.linprog=func
    # The objectives used here are nonnegative on every column. This is
    # itself an exact rational zero dual, useful when a proposal fails.
    assert all(reward(k,m)>=0 for k,_,m in lp.columns)
    return Q(0),Q(0)


def physical_floor(floors,cut,lo,hi):
    z={}
    for d in range(cut,1,-1):
        for k in reversed(range(d)):
            l=d-k-1;b=Q(l,d);a,r=floors[d,k]
            def expr(rho):
                g=l+rho
                return (a+r*rho+g*z.get((d+1,k),0)+b*z.get((d,k+1),0))/(hi+g+b)
            z[d,k]=Q(min(expr(x) for x in (0,Q(k+1,2)))*10**21//1,10**21)
            assert z[d,k]>=0
    return lo*z[2,0]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lower',type=Path,required=True);p.add_argument('--upper',type=Path,required=True)
    p.add_argument('--side',choices=['lower','upper']);p.add_argument('--cut',type=int,default=24,choices=[18,24,32,48])
    p.add_argument('--fast',action='store_true',help='diagnostic only: use pinned previously independently replayed non-f prices and denominator')
    a=p.parse_args();assert __debug__ and sys.maxsize>2**32
    compose=load('frozen_compose',HERE.parent/'replay.py')
    ext=compose.load('compose_extension',OCC/'high-parent'/'verify_release.py')
    base=ext.load('base');compose.check_sources(base);ext.check_sources();base.check_sources();base.check_receipts()
    orig=compose.load('compose_corrector',HERE.parents[1]/'corrector.py')
    source=ROOT/'pressure-memory'/'marginal'/'replay.py'
    assert sha256(source.read_bytes()).hexdigest()=='b35c7e4faacd697fd85d5957f33fa27ab0790ff7c111a4cc578328940d684a70'
    mem=load('frozen_memory',source);mem.corners(a.cut,2)
    off=mem.model(a.cut,False,2);on=mem.model(a.cut,True,2)
    assert off.rhs==on.rhs and len(off.columns)==len(on.columns)
    assert all((x,c,m)==(y,{r:v for r,v in d.items() if not (isinstance(r,tuple) and r[0]=='B-age')},n)
               for (x,c,m),(y,d,n) in zip(off.columns,on.columns,strict=True))
    print('matched flow/memory:',len(off.columns),'columns,',len(on.rows)-len(off.rows),'added memory rows',flush=True)
    modules={n:base.load(n,x) for n,x in base.MODULES.items()}
    edges,suffix,hp=map(ext.load,('edges','suffix','high_parent'))
    assert (suffix.LO,suffix.HI)==(hp.LO,hp.HI)==base.PRIOR==(mem.F.LO,mem.F.HI)
    masses,ell=suffix.newborn_kernel(False),suffix.newborn_kernel(True)
    paths={side:getattr(a,side).resolve(strict=True) for side in ('lower','upper')}
    before={side:base.fingerprint(path) for side,path in paths.items()}
    assert all(x[2]==base.SIZE for x in before.values()) and before['lower'][:2]!=before['upper'][:2]
    results={}
    for side,path in paths.items():
        if a.side and a.side!=side:continue
        mode,lam,c0,_,sha=base.PAYLOADS[side];sign=1 if side=='lower' else -1
        c=modules['denom'].Cert(path)
        assert (c.mode,c.lam,c.D,c.K,c.S)==(mode,Q(lam),18,17,base.SCALE)
        if not a.fast:
            old=compose.side_run(side,path,ext,base,modules,edges,suffix,hp,masses,ell,orig)
        else:
            assert base.digest(path)==sha
            old=compose.ENDPOINT[side]
        charge=modules['denom'].Charge(c);stock=modules['stock'];extr=modules['extract']
        slab=stock.table(extr.portable(c));floors=tuple(stock.envelope(slab[k]) for k in range(5))
        assert floors==tuple(map(Q,orig.FLOOR[mode]))
        rows=orig.operands(c,64)
        extracted=extr.buffer(c)['rows'];assert len(extracted)==len(rows)
        for raw in extracted:
            l,k,h=raw['state'];d=l+k+h;x,y=rows[d,k]
            assert h==1 and (x,x+y*Q(k+1,2))==tuple(map(Q,raw['signed_corners']))*(2 if d<18 else 1)
        recharge=Q('-.002' if side=='lower' else '-.009');s=compose.TRIAL[side]
        m,z,tail,M=compose.census(c,charge,rows,floors,s,recharge)
        assert m>=0
        _,occ=compose.upper_correction(z,tail,sign,base.PRIOR,64);den=suffix.DEN[side]+occ
        assert den>0
        res=residual(c,charge,rows,floors,s,recharge,z,tail,M,a.cut)
        direct=physical_floor(res,a.cut,*base.PRIOR)
        vals=[]
        for lp in (off,on):
            val,defect=price(mem,lp,cost(res),(side,'residual'))
            vals.append(val/10**6)
            print(side,'unused residual lower',val/10**6,'repair',defect/10**6,flush=True)
        # Same objective mask and off-flow constraints, in old witness units.
        first,_=price(mem,off,cost({x:rows[x] for x in res}),(side,'original'))
        paid=sign*(old-s)*den-m
        assert paid>0
        prev_flow=c.lam+sign*(first/10**6+paid)/suffix.DEN[side]
        reserve=max(direct,*vals)
        edge=old+sign*reserve/den
        assert sign*(edge-old)>0
        print(side,'original f (same mask, flow)',first/10**6,'previous flow comparator',base.outward(prev_flow,sign),
              'direct residual',direct,'seven frozen',paid,'den',den,'old',base.outward(old,sign),
              'NEW exact',edge,'outward',base.outward(edge,sign),flush=True)
        results[side]=base.outward(edge,sign)
    assert before=={side:base.fingerprint(path) for side,path in paths.items()}
    compose.check_sources(base);ext.check_sources();base.check_sources()
    print('PASS residual composition' if not a.fast else 'PASS fast frozen-ledger comparison',results,flush=True)


if __name__=='__main__':main()
