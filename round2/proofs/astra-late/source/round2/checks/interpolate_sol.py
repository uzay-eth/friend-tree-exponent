#!/usr/bin/env python3
"""Cold replay of Neon D18 pointwise correction and retained-charge incidence."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from hashlib import file_digest
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import argparse
import sys

BASE = Path(__file__).resolve().parents[1] / 'occupation-neon'
J = (1,0,1)
UNIT = 10**18
PRIOR = Q('.565673570322'), Q('.565775')
LEDGER = {
    'lower': (Q('798064281050388986950514709/625000000000000000000000000000000'),
              Q(58961819710386489,10**23), Q(55203067615510485475692014463,10**35), Q('1.356833021311')),
    'upper': (Q('1876286823343708393724246511/1250000000000000000000000000000000'), Q(46635843154908967,10**23),
              Q(139999781774155931356847859909,5*10**35), Q('1.367169995210')),
}
TICKS = {'lower': (565676930572,565676930574), 'upper': (565771634797,565771634792)}


def load(name,path):
    spec=spec_from_file_location(name,path)
    mod=module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


class Envelope:
    """Build the bound by independent raw motif incidences, before recharge."""
    def __init__(self,c):
        self.c=c; self.bottom={}; self.crown={}
        # The bottom's parent/grandparent slots are live; maximize over
        # their full table axes, not over the bottom's own ancestry flag.
        g=c.g[[c.state[x] for x in c.par],:]
        for x,j in c.child.items():
            self.bottom[x]=int(g[:,c.state[x]].max())+sum(max(0,v) for v in (
                int(c.hvtab[...,j].max()),int(c.p4tab[...,j].max()),
                max((u[1]-1)*int(c.ybtab[i,...,j].max()) for u,i in c.ybu.items()),
            ))
        for x,i in c.hc.items(): self.crown[x]=comb(x[1],2)*int(c.hctab[i].max())
        self.local=defaultdict(int)
        for x,i in c.yt.items(): self.local[x]+=comb(x[1],2)*int(c.yttab[:,i,:].max())
        for x,i in c.s3.items(): self.local[x]+=comb(x[1],3)*int(c.s3tab[i].max())
        self.slab={}
        for h in (0,1):
            for k in range(c.KF+2):
                lo=max(0,2-k-h)
                start=max(lo,c.D-k-h)
                lo_stock=range(lo,min(start,c.LF+1))
                finite=[self.raw((l,k,h)) for l in lo_stock]
                controls=c.tail(0,k,h)
                fmax=int(c.f[c.vi(start,k,h):c.vi(c.LF,k,h)+1].max())
                tail=fmax+(max(self.bottom[u] for u in controls) if h else 0)
                if 2<=k<=c.K: tail+=max(self.crown.get(u,0) for u in controls)
                self.slab[h,k]=max([tail,*finite])

    @cache
    def raw(self,x):
        l,k,h=x; c=self.c
        return (int(c.f[c.vi(*x)])+self.local[x]+sum(
            w*((self.bottom[u] if h else 0)+self.crown.get(u,0)) for u,w in c.embed(x)))


def solve(c,E,general,floor,s,t,flat=True,legacy_exit=False):
    eps=1 if c.mode==2 else -1
    delta=eps*(s-c.lam); assert delta>=0
    M=max(Q(v,c.S)+t*(k-h) for (h,k),v in E.slab.items())
    assert t<=0 and c.K<=c.KF
    cutoff=c.KF+1
    tail={}
    tests=0
    for k in reversed(range(cutoff+1)):
        r=delta*(M-Q(E.slab[1,k],c.S)-t*(k-1))+(floor[k] if k<len(floor) else 0)
        assert r>=0
        rates=(max(0,1-Q(k+1,65)),Q(1)) if k<cutoff else (Q(0),)
        u=tail.get(k+1,0)
        v=Q(min((r+b*u)*UNIT/(s+b) for b in rates)//1,UNIT) if flat else Q(0)
        assert all((s+b)*v<=r+b*u for b in rates)
        tail[k]=v; tests+=len(rates)
    z={}
    for d in reversed(range(2,65)):
        for k in reversed(range(d)):
            l=d-k-1; x=(l,k,1); b=Q(l,d)
            r,free=general.terms('f',(x,))
            slope={w for _,_,_,w in free}; assert len(slope)==1
            r+=delta*(M-E.raw(x)/c.S-t*(k-1))
            prev=(tail.get(k,Q(0)) if legacy_exit else tail[min(k,cutoff)]) if d==64 else z[d+1,k]
            bud=z[d,k+1] if b else 0
            rates=(Q(0),Q(k+1,2))
            vals=[(r+x*next(iter(slope))+(l+x)*prev+b*bud,s+l+x+b) for x in rates]
            z[d,k]=Q(min(u*UNIT/v for u,v in vals)//1,UNIT)
            assert all(v*z[d,k]<=u for u,v in vals); tests+=2
    margin=s*z[2,0]-delta*M
    return z,tail,margin,delta,M,tests


def check_trees(c,E,gen,embedding,s,t,sol,opposed=False):
    z,tail,m,delta,M,_=sol
    eps=1 if c.mode==2 else -1
    j=z[2,0]
    def v(x):
        l,k,h=x; d=sum(x)
        if not h:return 0
        return (z[d,k] if d<=64 else tail[min(k,c.KF+1)])-j
    cases=[((-1,),[2]),((-1,),[10**12]),((-1,0,1,2),[2,0,0,1]),
           ((-1,0,0,0,1,1,2,4),[2,0,0,1,0,1,1,18]),
           ((-1,0,1,1),[2,0,1,1]),((-1,0,1,1),[47,14000,1,1]),
           ((-1,0,1,2),[14001]*4),((-1,0,1,2),[10**6]*4),
           ((-1,0,1,2,1,0,5),[62,63,0,1,4,65,1])]
    for d in (17,18,63,64,65,66):
        cases.append(((-1,0,1,1),[2,d-3,1,2]))
    for k in (c.K-1,c.K,c.K+1,c.KF-1,c.KF,c.KF+1,c.KF+5):
        for l in (0,1,10**9):
            for h in (0,1):
                a=(l,k,h)
                assert E.raw(a)/c.S+t*(k-h)<=M
                if k>=c.KF+1:
                    assert E.raw(a)<=E.slab[h,c.KF+1]
            if k<=c.K+1 and l==1:
                cases.append(((-1,0)+(1,)*k,[2,l]+[1]*k))
    if opposed:
        cases=[((-1,0,1,2,1,0,5),[62,63,0,1,4,65,1])]
    seen=set(); count=0
    for p,l in cases:
        kids=[[u for u,w in enumerate(p) if w==i] for i in range(len(p))]
        x=tuple((l[i],len(ch),int(i>0)) for i,ch in enumerate(kids))
        old=c.potential(kids,x)/c.S
        if len(x)==1 and l==[2]:
            assert old>LEDGER['lower' if eps==1 else 'upper'][-1], 'average ceiling is not pointwise'
        assert old<=sum(E.raw(a)/c.S for a in x)
        assert sum(k-h for _,k,h in x)==0
        for a in x:
            _,k,h=a
            assert E.raw(a)/c.S+t*(k-h)<=M
        assert v(J)==0
        def potential(kk,xx,sign=eps,constant=True):
            return c.potential(kk,xx)/c.S+sign*(sum(map(v,xx))-(delta*M/s if constant else 0))
        def arrow(fun,discount):
            old=fun(kids,x); out=-discount*old
            for i,(l,k,h) in enumerate(x):
                rho=sum(Q(1,sum(x[w])) for w in kids[i]+([p[i]] if i else []))
                out+=(l+rho)*(fun(kids,x[:i]+((l+1,k,h),)+x[i+1:])-old)
                if l:
                    kk=[ch[:] for ch in kids]+[[]]; kk[i].append(len(x))
                    xx=x[:i]+((l-1,k+1,h),)+x[i+1:]+(J,)
                    out+=Q(l,sum(x[i]))*(fun(kk,xx)-old)
            return out
        charge=defaultdict(Q)
        for f,ids in embedding(p):
            selected=tuple(x[i] for i in ids)
            a,terms=gen.terms(f,selected)
            for j0,parent,n,w in terms:
                j1=ids[j0]; donors=([p[j1]] if parent else [u for u in kids[j1] if u not in ids])
                assert len(donors)==n and min(donors)>=0
                a+=w*sum(Q(1,sum(x[u])) for u in donors)
            assert a>=0
            charge[f]+=a
        seen.update(f for f,a in charge.items() if a)
        assert sum(charge.values())==eps*arrow(lambda kk,xx:c.potential(kk,xx)/c.S,c.lam)
        kept=sum(a for f,a in charge.items() if f!='f')
        new=eps*arrow(potential,s)
        if len(x)==1 and l==[2]:
            reversed_constant=2*eps*arrow(lambda kk,xx:potential(kk,xx,constant=False),s)-new
            assert reversed_constant<0  # a founder-only tree needs the constant with the displayed sign
        if opposed:
            assert new==Q('-902875978288208312083656815181293852064686619898589/1097886886904005745377459200000000000000000000000000000')
            assert new<0
            print('FAULTY whole-tree upper residual',new,'parents',p,'leaves',l,flush=True)
            return
        else:
            assert new>=m*(len(x)-1)+kept
        assert new-eps*arrow(lambda kk,xx:potential(kk,xx,constant=False),s)==delta*M
        count+=1
    assert not opposed, 'negative physical control not found'
    assert seen=={'f','g','hV','hC','P4','Yt','Yb','S3'}
    print('PASS',count,'exact whole-tree generator and all-eight incidence cases',flush=True)


def composition(side,s,eps,sol,higher):
    z,tail,m,_,_,_=sol
    j=z[2,0]; H=max(eps*(v-j) for v in (*z.values(),*tail.values()))
    lo,hi=PRIOR; nxt={}
    for d in range(64,1,-1):
        cur={}
        for k in reversed(range(d)):
            l=d-k-1; b=Q(l,d)
            gap=H-eps*(z[d,k]-j)
            assert gap>=0
            r=gap+b*cur.get(k+1,0)
            corners=[(r+g*nxt.get(k,0),hi+b+g) for g in (l,l+Q(k+1,2))]
            cur[k]=Q(min(a*UNIT/v for a,v in corners)//1,UNIT)
            assert all(a>=v*cur[k] for a,v in corners)
        nxt=cur
    shift=H-lo*nxt[0]
    old_higher,suffix,hp,den=LEDGER[side]
    assert higher>=old_higher
    for i,(label,D) in enumerate((('uniform',den+H),('occupation',den+shift))):
        endpoint=s+eps*(m+higher+suffix+hp)/D
        ticks=eps*((eps*endpoint*10**12)//1)
        assert ticks==TICKS[side][i]
        print(side,label,'signed denominator correction',D-den,'endpoint',endpoint,
              f'outward 0.{ticks:012d}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('lower',type=Path);p.add_argument('upper',type=Path)
    a=p.parse_args();sys.path.insert(0,str(BASE/'fringe'))
    G=load('cold_charges',BASE/'fringe/charges.py')
    C=load('original_corrector',BASE/'interpolate/corrector.py')
    R=load('slab_extract',BASE/'rows/extract.py')
    B=load('stock_floors',BASE/'rows/check.py')
    certs=[]
    for name,path,s,t in (('lower',a.lower,Q('.565675148075'),-Q(1,500)),
                          ('upper',a.upper,Q('.565770220167'),-Q(9,1000))):
        with path.open('rb') as stream: assert file_digest(stream,'sha256').hexdigest()==R.HASH[name]
        c=G.Cert(path);certs.append(c)
        assert (c.S,c.Sr,c.LF,c.D,c.K)==(10**11,10**9,14000,18,17)
        bands=B.table(R.portable(c)); floor=tuple(B.envelope(bands[k]) for k in range(5))
        expected=tuple(map(Q,('256918177771749/100000000000000000000','4272535757583/4000000000000000000','6713770119873/5000000000000000000','932799090941/625000000000000000','78795139786299/50000000000000000000') if name=='lower' else ('20697884432099/10000000000000000000','861071198001/1000000000000000000','5412445357533/5000000000000000000','6016226911043/5000000000000000000','31764287097597/25000000000000000000')))
        assert floor==expected
        print(name,'authenticated all-stock floors',*(str(x) for x in floor),flush=True)
        E=Envelope(c); gen=G.Charges(c)
        author=C.D.Charge(c); assert author.by_k==E.slab
        rows=C.operands(c,64)
        for d,k in rows:
            x=(d-k-1,k,1); r,ext=gen.terms('f',(x,))
            assert (r,{w for _,_,_,w in ext})==(rows[d,k][0],{rows[d,k][1]})
            assert E.raw(x)==author.exact(*x)
        q=solve(c,E,gen,floor,s,t,legacy_exit=True)
        if name=='upper':
            assert q[2]<0
            faulty=(floor[0],1000*floor[1],*floor[2:])
            bad=solve(c,E,gen,faulty,s,t,legacy_exit=True)
            assert bad[2]==Q('1167842799204166437/1000000000000000000000000000000')
            print('WITHDRAWN',s,'authenticated margin',q[2],'faulty margin',bad[2],flush=True)
            check_trees(c,E,gen,G.embedding,s,t,bad,opposed=True)
            s=Q('.565773278615')
            q=solve(c,E,gen,floor,s,t,legacy_exit=True)
        failed=solve(c,E,gen,floor,s+Q(1 if name=='lower' else -1,10**12),t,legacy_exit=True)
        assert failed[2]<0
        killed=solve(c,E,gen,floor,s,t,False,legacy_exit=True)
        z,tail,m,delta,M,n=q
        assert m>0>killed[2]
        assert solve(c,E,gen,floor,s,t)[2]>=m
        eps=1 if name=='lower' else -1
        vals=[*z.values(),*tail.values()]
        sup=max(eps*(x-z[2,0]) for x in vals)
        assert max(abs(x-z[2,0]) for x in vals)+delta*M/s<Q('.01')
        print(name,'s',s,'rows',len(rows),'Bellman corners',n,'cap',M,'margin',m,'killed',killed[2],
              'sup signed correction',sup,flush=True)
        check_trees(c,E,gen,G.embedding,s,t,q)
        mod=sys.modules[G.Extract.__module__]; mod.physical.MARK_CAP=2
        comp=mod.birth(mod.Extract(c),12)
        composition(name,s,eps,q,comp)
    # An opposed example: both raw witnesses have the same sign on each
    # of two trees at a common interior target, with opposite signs across trees.
    for l,sign in ((14001,-1),(10**6,1)):
        states=((l,1,0),(l,1,1),(l,0,1)); kids=[[1],[2],[]]
        residuals=[]
        for c in certs:
            def phi(x, ch): return c.potential(ch,x)/c.S
            old=phi(states,kids); residual=-Q('.56572365')*old
            for i,(stock,k,h) in enumerate(states):
                neighbors=kids[i]+([i-1] if i else [])
                d=stock+k+h
                gain=stock+sum(Q(1,sum(states[j])) for j in neighbors)
                residual+=gain*(phi(states[:i]+((stock+1,k,h),)+states[i+1:],kids)-old)
                after=[row[:] for row in kids]+[[]]; after[i].append(3)
                born=states[:i]+((stock-1,k+1,h),)+states[i+1:]+(J,)
                residual+=Q(stock,d)*(phi(born,after)-old)
            assert old>0 and sign*residual>0
            residuals.append(float(residual/old))
        print('opposed frozen chain',l,'normalized residuals',*residuals,flush=True)


if __name__=='__main__': main()
