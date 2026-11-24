#!/usr/bin/env python3
"""Authenticated composition of the first-body corrector and disjoint g/higher rewards."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
OCC=HERE.parents[1]
WITHDRAWN=Q('.565770220167')
REPAIRED_DRIFT=Q('3418020050304145090365704218566209189345801573143/219577377380801149075491840000000000000000000000000000')
TRIAL={'lower':Q('.565675148075'),'upper':Q('.565773278615')}
ENDPOINT={
    'lower':Q('15350577873033193312515081376688304908311/27136651758904633839660987238280000000000'),
    'upper':Q('15470132261475507150976552475495052519651/27343421462238657786289217467400000000000'),
}
INPUTS=(HERE.parent/'corrector.py',HERE/'sources.py',HERE/'lower.log',HERE/'upper.log',
        OCC/'high-parent'/'verify_release.py',OCC/'high-parent'/'inputs.txt',
        OCC.parents[0]/'checks'/'interpolate_sol.py')


def load(name,path):
    spec=spec_from_file_location(name,path);mod=module_from_spec(spec)
    sys.modules[name]=mod;spec.loader.exec_module(mod)
    return mod


def check_sources(base):
    expected={str(p.relative_to(base.ROOT)) for p in INPUTS}
    seen=set()
    for line in (HERE/'inputs.txt').read_text().splitlines():
        sha,name=line.split()
        assert name in expected-seen and re.fullmatch('[a-f0-9]{64}',sha),name
        assert base.digest(base.ROOT/name)==sha,name
        seen.add(name)
    assert seen==expected


def census(c,chg,rows,floors,s,t,n=64,tail_kind='full',deficit=True,reward=True):
    sign=1 if c.mode==2 else -1
    delta=sign*(s-c.lam);M=chg.cap(t*c.S)/c.S
    assert delta>=0
    z=defaultdict(Q);tail=defaultdict(Q)
    for k in range(c.KF+1,-1,-1):
        a=Q(0)
        if tail_kind in ('stock','full') and k<len(floors) and reward:a+=floors[k]
        if tail_kind in ('deficit','full') and deficit:
            a+=delta*(M-Q(chg.by_k[1,k],c.S)-t*(k-1))
        assert a>=0
        nxt=k+1 if k<c.KF+1 else k
        bmin=max(0,1-Q(k+1,n+1)); bmax=Q(1)
        tail[k]=down(a/s if k==nxt else min((a+b*tail[nxt])/(s+b) for b in (bmin,bmax)))
        assert all((s+b)*tail[k] <= a+b*tail[nxt] for b in (bmin,bmax))
    for d in range(n,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1;b=Q(l,d)
            a,r=rows[d,k] if reward else (0,0)
            deficit_x=M-chg.exact(l,k,1)/c.S-t*(k-1)
            assert deficit_x>=0
            if deficit:a+=delta*deficit_x
            y=z[d+1,k] if d<n else tail[k]
            u=z[d,k+1]
            rhs=[(a+r*rho+(l+rho)*y+b*u)/(s+l+rho+b) for rho in (0,Q(k+1,2))]
            z[d,k]=down(min(rhs))
            assert 0<=z[d,k]<=min(rhs)
    margin=s*z[2,0]-delta*M
    return margin,z,tail,M


def down(x):return Q(x*10**18//1,10**18)


def forward(c,charge,rows,s,recharge,z,tail,M):
    """Forward discounted occupation under the selected finite census policy."""
    unit=10**32
    delta=abs(s-c.lam)
    lo,hi=defaultdict(int),defaultdict(int)
    lo[2,0]=hi[2,0]=unit
    low=high=Q(0)
    for d in range(2,65):
        for k in range(d):
            l=d-k-1;b=Q(l,d)
            a,r=rows[d,k]
            a+=delta*(M-charge.exact(l,k,1)/c.S-recharge*(k-1))
            y=z[d+1,k] if d<64 else tail[k]
            u=z[d,k+1]
            rho=min((0,Q(k+1,2)),key=lambda p:(a+r*p+(l+p)*y+b*u)/(s+l+p+b))
            g=l+rho; den=s+g+b
            payout=(a+r*rho+(g*y if d==64 else 0))/den
            assert payout>=0
            low+=lo[d,k]*payout/unit;high+=hi[d,k]*payout/unit
            for rate,dest in ((b,(d,k+1)),(g,(d+1,k)) if d<64 else (0,None)):
                if rate:
                    lo[dest]+=lo[d,k]*rate//den
                    hi[dest]+=-(-hi[d,k]*rate//den)
    # The backward rounding loses at most one 10^-18 tick per visited
    # census state. The positive flow bounds bracket the unrounded score.
    assert z[2,0]<=high and low<=high<=z[2,0]+Q(len(rows)+1,10**18)
    assert high-low<Q(1,10**26)
    print('PASS finite forward discounted census',len(rows),'states; value',fmt(low),flush=True)


def upper_correction(z,tail,sign,prior,n):
    v0=z[2,0]; vals=[sign*(v-v0) for v in (*z.values(),*tail.values())]
    H=max(vals); assert H>=0
    def gap(d,k):return H-sign*((z[d,k] if d<=n else tail[min(k,maxk)])-v0)
    # Physical stable-age payment of nonnegative gaps, with degree-cut
    # tail at zero. This sharpens one-sided sup without integrating signs.
    maxk=max(tail);lo,hi=prior
    w=defaultdict(Q)
    for d in range(n,1,-1):
        for k in range(d-1,-1,-1):
            l=d-k-1;b=Q(l,d);a=gap(d,k)+b*w[d,k+1]
            assert gap(d,k)>=0
            w[d,k]=down(min((a+g*w[d+1,k])/(hi+g+b) for g in (l,l+Q(k+1,2))))
    return H,H-lo*w[2,0]


def fmt(v):return f'{float(v):.14g}'


def counterexample(c,charge,rows,floors,recharge,repaired):
    if c.mode!=1:return
    parents=(-1,0,1,2,1,0,5);leaves=(62,63,0,1,4,65,1)
    x=tuple((l,parents.count(i),int(i!=0)) for i,l in enumerate(leaves))
    children=[[j for j in range(len(x)) if parents[j]==i] for i in range(len(x))]
    from audit import build, arrows, physical_view
    physical=build(parents,leaves)
    assert tuple(physical_view(physical)[1])==x
    for planted,s in ((True,WITHDRAWN),(False,repaired)):
        f=list(floors)
        if planted:f[1]*=1000
        m,z,tail,M=census(c,charge,rows,f,s,recharge)
        def phi(p):
            ch,st,_=physical_view(p)
            V=sum((z[sum(y),y[1]] if sum(y)<=64 else tail[min(y[1],c.KF+1)])-z[2,0] for y in st[1:])
            return c.potential(ch,st)/c.S-V+(c.lam-s)*M/s
        old=phi(physical)
        signed=s*old-sum(r*(phi(physical+(i,))-old) for i,r in arrows(physical).items())
        exact=Q('-902875978288208312083656815181293852064686619898589/1097886886904005745377459200000000000000000000000000000')
        if planted:assert signed==exact and m>0
        else:assert signed==REPAIRED_DRIFT and signed>=m*(len(x)-1)>0
        print('whole-tree regression', 'planted x1000' if planted else 'authenticated floors',
              'target',s,'local margin',fmt(m),'physical upper drift',signed,flush=True)


def best_trial(c,charge,rows,floors,recharge):
    sign=1 if c.mode==2 else -1
    low,high=0,10**7
    while high-low>1:
        j=(low+high)//2;s=c.lam+sign*Q(j,10**12)
        if census(c,charge,rows,floors,s,recharge)[0]>=0:low=j
        else:high=j
    return c.lam+sign*Q(low,10**12)


def side_run(side,path,ext,base,modules,edges,suffix,hp,masses,ell,orig):
    sign=1 if side=='lower' else -1
    mode,lam,c0,_,sha=base.PAYLOADS[side]
    assert base.digest(path)==sha
    print(side,'payload SHA-256',sha,flush=True)
    c=modules['denom'].Cert(path)
    assert (c.mode,c.lam,c.D,c.K,c.S)==(mode,Q(lam),18,17,base.SCALE)
    # Independently re-extract f only to verify the *pointwise* correction's
    # source and all-stock floors; never book an f occupation numerator.
    stock=modules['stock'];extr=modules['extract']
    table=stock.table(extr.portable(c))
    floors=tuple(stock.envelope(table[k]) for k in range(5))
    assert floors==tuple(map(Q,orig.FLOOR[mode]))
    raw=extr.buffer(c)
    rows=orig.operands(c,64)
    for row in raw['rows']:
        l,k,h=row['state'];d=l+k+h;a,r=rows[d,k]
        expected=tuple(map(Q,row['signed_corners']))
        assert h==1 and (a,a+r*Q(k+1,2))==(expected if d>=18 else expected*2)
    print(side,'authenticated pointwise f rows',len(rows),'and 5 stock floors',flush=True)
    ceiling,fixtures,recharge=base.improved_denominator(modules,c)
    den=suffix.DEN[side]
    assert den==Q(-(-ceiling*10**12//1),10**12)
    charge=modules['denom'].Charge(c)
    s=best_trial(c,charge,rows,floors,recharge)
    assert s==TRIAL[side]
    counterexample(c,charge,rows,floors,recharge,s)
    print(side,'repaired pointwise trial',base.outward(s,sign),flush=True)
    m,z,tail,M=census(c,charge,rows,floors,s,recharge)
    forward(c,charge,rows,s,recharge,z,tail,M)
    ref=orig.solve(c,charge,rows,64,s)
    assert (m,z[2,0],M,dict(tail))==(ref[0],ref[1],ref[2],ref[4]) and m>=0
    sup,occ=upper_correction(z,tail,sign,base.PRIOR,64)
    assert den+occ>0 and occ<=sup
    assert Q(c0,4*c.S)>max(abs(v-z[2,0]) for v in (*z.values(),*tail.values()))+abs(s-c.lam)*M/s
    for kind,deficit,reward in (('zero',False,True),('zero',True,True),('stock',True,True),
                                ('deficit',True,True),('full',True,False),('full',True,True)):
        a,zz,tt,_=census(c,charge,rows,floors,s,recharge,tail_kind=kind,deficit=deficit,reward=reward)
        assert (a==m>0) if kind=='full' and reward else a<0
        print(side,'ablation',kind,'deficit',deficit,'f',reward,'margin',fmt(a),
              'zJ',fmt(zz[2,0]),'tail k=0,1,2',*(fmt(tt[k]) for k in range(3)),flush=True)
    print(side,'new denominator',fmt(den+occ),'old',den,'correction sup/occupation',fmt(sup),fmt(occ),flush=True)

    higher=modules['higher'];higher.physical.MARK_CAP=suffix.P.MARK_CAP=2
    high=higher.birth(higher.Extract(c),12)
    assert high>=suffix.OLD[side]>0
    print(side,'higher ≤12',fmt(high),flush=True)
    rr=edges.rows(extr.edges_buffer(c));assert len(rr)==43776 and min(rr.values())==suffix.EDGE[side]>0
    del rr
    size20=suffix.floor(suffix.evaluate(suffix.H.Extract(c),side,20,masses,ell))
    assert size20==ext.EXPECTED[side][0]
    print(side,'joint size 13–20 and actual inner-edge floor',size20,flush=True)
    phy=hp.Physical(c)
    mapped={(d,k):hp.charges(phy,d,k) for d in range(18,192) for k in range(1,11)}
    mixed=hp.solve(18,191,10,lambda d,k,j:mapped[d,k][0]*j+mapped[d,k][1]*comb(j,2) if j else 0)
    assert mixed==ext.EXPECTED[side][1]
    phy.value.cache_clear()
    assert set(suffix.H.SPEC)-{'f','g'}==base.HIGHER
    for cap,bound in (('uniform',den+sup),('occupation',den+occ)):
        paid=high+size20+mixed
        full=s+sign*(m+paid)/bound
        print(side,cap,'non-f',fmt(paid),'high-parent',mixed,'exact endpoint',full,
              'outward',base.outward(full,sign),flush=True)
    endpoint=s+sign*(m+high+size20+mixed)/(den+occ)
    assert endpoint==ENDPOINT[side]
    return endpoint


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lower',type=Path,required=True);p.add_argument('--upper',type=Path,required=True)
    p.add_argument('--side',choices=['lower','upper'])
    a=p.parse_args();assert __debug__ and sys.version_info>=(3,11) and sys.maxsize>2**32
    sys.dont_write_bytecode=True
    ext=load('compose_extension',OCC/'high-parent'/'verify_release.py')
    base=ext.load('base');check_sources(base);ext.check_sources();base.check_sources();base.check_receipts()
    orig=load('compose_corrector',HERE.parent/'corrector.py')
    modules={n:base.load(n,x) for n,x in base.MODULES.items()}
    edges,suffix,hp=map(ext.load,('edges','suffix','high_parent'))
    assert (suffix.LO,suffix.HI)==(hp.LO,hp.HI)==base.PRIOR
    masses,ell=suffix.newborn_kernel(False),suffix.newborn_kernel(True)
    paths={side:getattr(a,side).resolve(strict=True) for side in ('lower','upper')}
    before={side:base.fingerprint(path) for side,path in paths.items()}
    assert all(path.is_file() and before[side][2]==base.SIZE for side,path in paths.items())
    assert before['lower'][:2]!=before['upper'][:2]
    if a.side is None:
        subprocess.run([sys.executable,'-Bu',str(INPUTS[-1]),str(paths['lower']),str(paths['upper'])],check=True)
    out={side:side_run(side,path,ext,base,modules,edges,suffix,hp,masses,ell,orig)
         for side,path in paths.items() if a.side is None or a.side==side}
    assert before=={side:base.fingerprint(path) for side,path in paths.items()}
    base.check_sources();ext.check_sources();check_sources(base)
    print('PASS composition', {side:base.outward(v,1 if side=='lower' else -1) for side,v in out.items()},flush=True)


if __name__=='__main__':main()
