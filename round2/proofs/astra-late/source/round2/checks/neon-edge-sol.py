#!/usr/bin/env python3
"""Independent EDGE source/row/physical-law checks; --full repeats exact composition."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
O=ROOT/'round2/occupation-neon'; HERE=O/'interpolate/edge'
J=(2,0); NB=(1,0,1)
PIN=('ff7803d6d39078d14fe47d464be77c677cdcad381e06c4b1caa054cd0ef5cda6',
     'e90651922acb6b3a427872cf3801f4b9097d2b83c7abe0ab75add6e60c78350f',
     '7a5979e12e1c15b20995d65c228140aba277b750424423062f820ee7564ccb5a')
BAD=Q('-902875978288208312083656815181293852064686619898589/1097886886904005745377459200000000000000000000000000000')
SEVEN=((-1,0,1,2,1,0,5),(62,63,0,1,4,65,1))


def load(name,path):
    spec=spec_from_file_location(name,path); m=module_from_spec(spec);spec.loader.exec_module(m);return m


def topology(p):return [[j for j,q in enumerate(p) if q==i] for i in range(len(p))]


def fixture(c,gen,embedding,oracle,s,z,tail,values,delta,M,reserve):
    eps=gen.sign
    def census(a):return sum(a),a[1]
    def score(p,x):
        tally=-delta*M/s
        for i,a in enumerate(x[1:],1):
            d,k=census(a);tally+=z.get((d,k),tail[min(k,c.KF+1)])-z[J]
            if p[i]>0:tally+=values.get((census(x[p[i]]),(d,k)),0)
        return tally

    cases=[SEVEN,((-1,), (2,)),((-1,), (10**9,)),
           ((-1,0,1,1,2,0,5),(2,5,4,1,1,5,1)),
           ((-1,0,1,1,2,0,5),(13,14,15,1,2,37,1)),
           ((-1,0,1,1,2,0,5),(64,1,1,1,2,14001,1)),
           ((-1,0,1,2,1),(1,2,1,1,2))]
    born=pre=founder=0
    for p,l in cases:
        kids=topology(p);x=tuple((l[i],len(ys),int(i>0)) for i,ys in enumerate(kids))
        old=score(p,x); corr=-s*old; base=c.potential(kids,x)/c.S; orig=-c.lam*base
        charged=defaultdict(Q); stockp=[]
        for f,ids in embedding(p):
            a,free=gen.terms(f,tuple(x[i] for i in ids))
            for j,h,n,b in free:
                i=ids[j]; ds=[p[i]] if h else [k for k in kids[i] if k not in ids]
                assert len(ds)==n and min(ds)>=0
                a+=b*sum(Q(1,sum(x[k])) for k in ds)
            assert a>=0
            charged[f]+=a
            if f=='g':
                u,w=ids; xx,yy=x[u],x[w]; a0,b,_,_=oracle.row(c,(Q(sum(xx)),*xx[1:]),(Q(sum(yy)),*yy[1:]))
                actual=(Q(1,sum(x[p[u]])) if u else 0,
                        sum(Q(1,sum(x[k])) for k in kids[u] if k!=w),
                        sum(Q(1,sum(x[k])) for k in kids[w]))
                assert eps*(a0+sum(t*r for t,r in zip(b,actual)))==a
                if u:
                    keep=reserve(census(xx),census(yy))
                    assert 0<=keep<=a
                    stockp.append(keep)
        for i,(L,k,h) in enumerate(x):
            d=L+k+h
            # Rebuild rates from individual physical arrows; external sums
            # are not new clocks, nor independent copies of included sends.
            rr=sum(Q(1,sum(x[w])) for w in kids[i]+([p[i]] if i else []))
            assert 0<=rr<=Q(k+h,2)
            xx=x[:i]+((L+1,k,h),)+x[i+1:]
            corr+=(L+rr)*(score(p,xx)-old)
            orig+=(L+rr)*(c.potential(kids,xx)/c.S-base)
            if not L:continue
            out=x[:i]+((L-1,k+1,h),)+x[i+1:]+(NB,); pp=p+(i,)
            after=census(out[i]); before=census(x[i])
            assert (after,before)==((d,k+1),(d,k)) and census(out[-1])==J
            new=values.get((after,J),0) if i else 0
            wrong=values.get((before,J),0) if i else 0
            born+=bool(new); pre+=new!=wrong;founder+=not i and bool(values.get((after,J),0))
            corr+=Q(L,d)*(score(pp,out)-old)
            orig+=Q(L,d)*(c.potential(topology(pp),out)/c.S-base)
        assert eps*orig==sum(charged.values())
        signed=eps*orig-delta*base+corr
        margin=s*z[J]-delta*M
        if s==Q('.565770220167'):
            if (p,l)==SEVEN:assert signed==BAD and margin>0
        else:
            unspent=sum(stockp)+sum(v for f,v in charged.items() if f not in ('f','g'))
            assert signed>=margin*(len(x)-1)+unspent
    if values:assert born and pre and founder
    return born,pre,founder


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('blob',type=Path,nargs='+');ap.add_argument('--full',action='store_true');a=ap.parse_args()
    for name,h in zip(('README.md','replay.py','regression.py'),PIN):assert sha256((HERE/name).read_bytes()).hexdigest()==h
    R=load('edge_under_audit',HERE/'replay.py')
    for path,h in R.SOURCE_HASHES.items():assert sha256((ROOT/'round2'/path).read_bytes()).hexdigest()==h,path
    oracle_path=O/'all-contrast/check/row_oracle.py'
    assert sha256(oracle_path.read_bytes()).hexdigest()=='9e79d6307332aac8b2e8d5a74a6cc03350dee3f4d7b52d2d82d25ff2869f2461'
    I=load('cold_census',ROOT/'round2/checks/interpolate_sol.py');oracle=load('cold_edge',oracle_path)
    G=load('cold_general',O/'fringe/charges.py'); X=load('cold_slabs',O/'rows/extract.py'); B=load('cold_stock',O/'rows/check.py')
    E=load('cold_inner',O/'rows/check_edges.py');HP=load('cold_high_parent',O/'high-parent/replay.py')
    for blob in a.blob:
        c=R.D.Cert(blob); gen=G.Charges(c); eps=gen.sign; mode=c.mode
        s=Q('.565677930574' if eps==1 else '.565770634792');delta=eps*(s-c.lam); assert delta>0
        cs=R.Census(c)
        floors=tuple(B.envelope(b) for b in B.table(X.portable(c)).values());assert floors==cs.floors
        for (d,k),(v,w) in cs.rows.items():
            base,terms=gen.terms('f',((d-k-1,k,1),));assert (base,{t[3] for t in terms})==(v,{w})
        inner=E.rows(X.edges_buffer(c));assert len(inner)==43776 and min(inner.values())==R.INNER[mode]
        cap,kmax=((40,8) if eps==1 else (32,7)) if a.full else (19,3)
        hp=HP.Physical(c); reserved={}
        # All high-parent virgin reserve rows within the production edge cut.
        dmax,kprod=(40,8) if eps==1 else (32,7)
        for d,k in product(range(18,dmax+1),range(1,kprod+1)):
            u=(Q(d),k,1);w=(Q(2),0,1);v,b,_,_=oracle.row(c,u,w)
            signed=eps*v+min(0,eps*b[0]*Q(k,2));assert b[0]==b[1]
            assert signed==HP.charges(hp,d,k)[0]
            reserved[d,k]=signed
        def reserve(p,w):
            return R.INNER[mode] if max(p[0],w[0])<18 else (reserved.get(p,0) if w==J else 0)
        # Independently contract sampled mixed/both-tail ports, and plant
        # the historical omitted f→g source; it changes the intercept.
        mixed=badrows=0
        for d,k,e,j in product((2,5,17,18,19,37,64),(1,3),(2,17,18,37),(0,2)):
            if not (k<d and j<e):continue
            v,b,_,trace=oracle.row(c,(Q(d),k,1),(Q(e),j,1))
            base,slopes=R.ports(gen,(d,k),(e,j))
            assert base==eps*v and slopes==(eps*b[0],eps*b[2] if j else 0)
            mixed+=1
            badrows+=eps*(v-sum(t[0] for t in trace))<0
        assert badrows
        edge=R.Edge(gen,cap,kmax,'disjoint');V=edge.solve(s)
        m,norm,z,tail=cs.solve(s,V)
        if eps<0:
            wrong=list(floors);wrong[1]*=1000; bads=Q('.565770220167')
            _,_,bz,bt=cs.solve(bads,{},floors=wrong,legacy=True)
            fixture(c,gen,G.embedding,oracle,bads,bz,bt,{},c.lam-bads,cs.M,lambda p,w:0)
            assert cs.solve(bads,{})[0]<0
        counts=fixture(c,gen,G.embedding,oracle,s,z,tail,V,delta,cs.M,reserve)
        # Frozen exact source fractions, not the old endpoints or optimistic
        # upper approximations. The same physical witness supplies all terms.
        side='lower' if eps==1 else 'upper'; high,size,hp0,den=I.LEDGER[side]
        assert den==R.DEN[mode]
        floors_paid={1:('1.501029458674e-6','4.663584315490e-7','2.799995635483e-7'),
                     2:('1.276902849680e-6','5.896181971038e-7','5.520306761551e-7')}
        terms=tuple(map(Q,floors_paid[mode]));assert all(v<=cap for v,cap in zip(terms,(high,size,hp0)))
        assert sum(terms)==R.RETAINED[mode]
        H=max(eps*(v-z[J]) for v in (*z.values(),*tail.values()))+(max(V.values()) if eps==1 else 0)
        paid=m+R.RETAINED[mode]
        print(side,'f',len(cs.rows),'stock floors',len(floors),'inner g',len(inner),'mixed g',mixed,'omitted-source negatives',badrows,
              'mixed-reserve',len(reserved),'edge states',len(V),'birth/pre-state/founder controls',counts,flush=True)
        if a.full:
            assert m<0<paid and H>=0 and norm<Q('.01')
            endpoint=s+eps*Q(paid,den+H)
            ticks=eps*((eps*endpoint*10**15)//1)
            assert ticks==(565677980686114 if eps==1 else 565770513279953)
            print(side,'exact margin',m,'paid',paid,'denominator',den+H,'exact endpoint',endpoint,'outward',f'0.{ticks:015d}',flush=True)


if __name__=='__main__':main()
