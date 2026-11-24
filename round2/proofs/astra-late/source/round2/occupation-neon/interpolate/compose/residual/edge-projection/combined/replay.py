#!/usr/bin/env python3
"""Compose only unspent post-EDGE slack with the frozen 02d0bd8 ledger."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[6]
PROJ=HERE.parent/'replay.py'
BASE={2:Q('191976090122981726641370633087/339373453939523193250000000000'),
      1:Q('24172191945911554486185483421/42724375658564528593750000000')}
MIN={2:Q('1.03096e-7'),1:Q('5.99909e-8')}


def pin_manifest(path):
    seen=set()
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        h,rel=line.split();assert rel not in seen and len(h)==64
        seen.add(rel)
        assert sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
    assert seen


def load(path,name):
    spec=spec_from_file_location(name,path);m=module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lower',type=Path,required=True);ap.add_argument('--upper',type=Path,required=True)
    ap.add_argument('--limit',type=int,choices=(18,24,32),default=18)
    ap.add_argument('--exact',type=int,choices=range(4,8),default=4)
    ap.add_argument('--side',choices=('lower','upper'))
    opt=ap.parse_args();assert __debug__
    pin_manifest(HERE/'inputs.sha256')
    edgehome=ROOT/'round2/occupation-neon/interpolate/edge'
    pin_manifest(edgehome/'inputs.sha256')
    pin_manifest(edgehome.parent/'compose'/'inputs.txt')
    P=load(PROJ,'post_edge_projection')
    assert P.ROOT==ROOT
    for path,h in P.PINS.items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    E,mem=(P.module(path,name) for path,name in zip(P.PINS,('combined_edge','combined_memory'),strict=True))
    R=load(HERE.parents[1]/'replay.py','combined_residual')
    mem.corners(opt.limit,2)
    off=P.make_model(mem,opt.limit,opt.exact,False);on=P.make_model(mem,opt.limit,opt.exact,True)
    assert off.rhs==on.rhs and len(off.columns)==len(on.columns)
    assert all((i,x,a)==(j,{r:v for r,v in y.items() if not (isinstance(r,tuple) and r[0]=='B-age')},b)
               for (i,x,a),(j,y,b) in zip(off.columns,on.columns,strict=True))
    print('matched shared-donor models',len(on.columns),'columns,',len(on.rows)-len(off.rows),'memory rows',flush=True)
    out={}
    for side,path in (('upper',opt.upper),('lower',opt.lower)):
        if opt.side and side!=opt.side:continue
        before=path.stat();c=E.Charges(E.D.Cert(path));mode=1 if side=='upper' else 2
        assert mode==c.c.mode and E.PAYLOAD_HASHES[mode]==P.digest(path)
        s,n,k=P.SPEC[mode];s=Q(s);eps=c.sign
        assert eps==(1 if side=='lower' else -1)
        edge=E.Edge(c,n,k,'disjoint');cs=E.Census(c.c)
        V=edge.solve(s);m,_,z,tail=cs.solve(s,V)
        H=max(eps*(w-z[E.J]) for w in (*z.values(),*tail.values()))
        if eps==1:H+=max(V.values())
        den=E.DEN[mode]+H;reserved=E.RETAINED[mode];net=m+reserved
        assert m<0<net and den>0
        old=s+eps*net/den;assert old==BASE[mode]
        print(side,'authenticated',E.PAYLOAD_HASHES[mode],'trial',s,'margin',m,'reserved',reserved,
              'denominator',den,'exact EDGE',old,flush=True)
        f=P.freward(E,c,cs,s,V,z,tail,opt.limit);g=P.greward(E,c,edge,s,V,opt.exact)
        assert len(f)==sum(range(2,opt.limit+1))
        assert len(g)==sum(range(1,opt.exact))*sum(range(2,opt.exact+1))
        assert all(p[0]<=opt.exact<18 and w[0]<=opt.exact for p,w in g)
        assert all((p,w) in edge.rows for p,w in g)
        obj=P.objective(f,g);fo=P.objective(f,{});go=P.objective({},g)
        for lp in (off,on):
            for kind,_,meta in lp.columns:
                a,b=fo(kind,meta),go(kind,meta)
                assert min(a,b)>=0 and obj(kind,meta)==a+b
                if kind=='edge' and b:
                    p,w,*_=meta
                    assert p[0]==p[1]<=opt.exact and w[0]==w[1]<=opt.exact
                if kind=='occ' and a:assert meta[0][0]<=opt.limit
                if kind not in ('occ','edge'):assert not (a or b)
        print(side,'mask rows',len(f),len(g),'higher families and mixed g untouched',flush=True)
        vals=[]
        for label,lp in (('flow',off),('memory',on)):
            val,repair=R.price(mem,lp,obj,(side,label,'joint'))
            vals.append(val/10**6)
            print(side,label,'joint lower',vals[-1],'dual repair',repair/10**6,flush=True)
        if opt.limit==18 and opt.exact==4:assert vals[-1]>=MIN[mode]
        paid=max(vals);assert paid>0
        end=s+eps*(net+paid)/den;gain=eps*(end-old)
        assert gain==paid/den
        render=lambda x:E.D.show(x,15,upper=eps<0)
        print(side,'CANDIDATE exact',end,'outward',render(end),'gain',E.D.show(gain,18,upper=False),
              'reward needed for 1e-6',E.D.show(den/10**6,16),
              'required multiple',E.D.show((den/10**6)/paid,8),flush=True)
        after=path.stat();assert (before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)==(
            after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)
        out[side]=render(end)
    pin_manifest(HERE/'inputs.sha256');pin_manifest(edgehome/'inputs.sha256')
    print('PASS post-EDGE composition candidate (pending cold review)',out,flush=True)


if __name__=='__main__':main()
