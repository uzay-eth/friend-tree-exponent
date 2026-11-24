#!/usr/bin/env python3
"""Independent slice/physical-arrow controls and exact Quintet reduced LP replay."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from pathlib import Path
from unittest.mock import patch

R=Path(__file__).resolve().parents[1]
D=R/'occupation-quintet/lifetime/dual'
PINS={
    'occupation-quintet/lifetime/dual/export.py':'17ebb88215bba02d2de549797c05df67a1adbd326775c92acb764b60d10b42c1',
    'occupation-quintet/lifetime/dual/sided.py':'4328c139ec0cb162b756b60181ce9ee625ae20a431fa5250786b7bff7959bbd9',
    'occupation-quintet/lifetime/dual/lower.json':'d698df2bdd20fd6d5a802f84d18124a1fcd6f491bb4dde8ae5038ca5d27f6169',
    'occupation-quintet/lifetime/dual/upper.json':'a23fe15ac518c35f9e7dcc92bc1fa497253f16e7677a9d94e572815682253fa8',
    'occupation-quintet/lifetime/check.py':'f2331f6a43799b6dd5eebebd655e9984f441ef72201dc75cc5624b933c472b7a',
    'occupation-quintet/census/census.py':'b6c133acc92ef36d651a92f4684852fc5ed241e8c0e88d3ac0fdad5554c2550b',
    'checks/quintet-lifetime-sol.py':'41f258fe768c300778a1657065dcd9ca4cf7caf6fac5d68c3c3b2fa2221b7d7b',
    'checks/stationary-p4-sol.py':'7432857b1b9fc05b7850cb9d2a6f22fc8200d22dc77f134bbd75eb69da2f0653',
    'occupation-neon/pressure-memory/birth-law/replay.py':'92a99ef12d8567798c4ed4dc3beb6535d3f18ef3ad5028aff0ec84e3684a5ba8',
}
SOURCE='b04cc50c1e391e0420eb7b006ac781d484c68543'
RECEIPT='d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb'
PRICE={
    'lower':Q(1276783006042809110549818991696409322453507199781,13644757069600662706896800000000000000000000),
    'upper':Q(212108873684137139278951921799899102453819851509,1240432460872787518808800000000000000000000),
}


def load(path):
    s=spec_from_file_location(Path(path).stem,R/path)
    m=module_from_spec(s); s.loader.exec_module(m); return m


def physical(s,A,W,rows,cap,hi):
    n=Counter()
    # Founder 0, focal nonfounder 1, and true side children with unequal
    # degrees; add all remaining stocks as literal ordinary vertices.
    for d,c in [(2,1),(2,2),(3,1),(3,2),(4,3),(26,26),(27,26),(27,27),
                (30,30),(31,30),(31,31),(63,29),(64,29),(64,30),(64,31),(65,30),
                (179,27),(180,27),(1999,30),(2000,30),(12000,30),(12001,30)]:
        k=c-1; L=d-c
        p=(-1,0)+(1,)*k+(0,)*3+tuple(i for i in range(2,k+2) for _ in range(1+(i%3)))+(1,)*L
        nb=[[] for _ in p]
        for j,i in enumerate(p[1:],1):nb[i].append(j);nb[j].append(i)
        deg=list(map(len,nb))
        assert deg[1]==d and sum(deg[j]>1 for j in nb[1])==c
        u=sum((Q(1,deg[j]) for j in nb[1] if j!=0 and deg[j]>1),Q())
        v=Q(1,deg[0]); receive=u+v
        F=cap-s.sign*s.lam*s.f(d,c)
        drift=Q(0); gain=bud=Q(0)
        values=(W.get((d,c)),W.get((d+1,c)),W.get((d,c+1),0))
        old=s.f(d,c); df=s.f(d+1,c)-old
        for a in nb[1]:
            rate=Q(1,deg[a]); gain+=rate
            F+=s.sign*rate*df
            if values[0] is not None and values[1] is not None:drift+=rate*(values[1]-values[0])
        for a in nb[1]:
            if deg[a]!=1:continue
            rate=Q(1,d); bud+=rate
            F+=s.sign*rate*(s.f(d,c+1)-old+s.j+s.nb(d,c+1))
            if values[0] is not None:drift+=rate*(values[2]-values[0])
        if d<180 and c<=26:
            F-=s.sign*(receive-Q(c*s.tick(s.slot(d,c)),A.SR))*df
        intercept,slope=s.row(d,c)
        assert F==cap+intercept+slope*receive and F>=0
        assert gain==L+receive and bud==Q(L,d)
        if (d,c) in rows:
            a,b=rows[d,c]; want=F+drift-hi*W[d,c]
            assert want==a+b*receive>=0
            n['masked arrows']+=1
            n['drop child receive rejected']+=a+b*v!=want
            n['drop mother receive rejected']+=a+b*u!=want
            if d==64:
                assert values[1]>0 and want-gain*values[1]!=want
                n['real degree-65 successor']+=1
        else:
            n['unpaid seam/overflow']+=1
            if (d,c)==(31,31):assert u>0 and v>0 and values[0] is None
        n['physical fixtures']+=1
    assert all(n[z] for z in ('drop child receive rejected','drop mother receive rejected','real degree-65 successor'))
    return dict(n)


def rows_check(side,A,lo,hi,g):
    z=json.loads((D/f'{side}.json').read_bytes())
    assert (z['schema'],z['side'],z['source_commit'],z['receipt_sha256'])==('quintet-f-dual-v1',side,SOURCE,RECEIPT)
    assert (Q(z['LO']),Q(z['HI']),z['cap_f_ticks'],z['payload_sha256'])==(lo,hi,g['cap_f'],g['verify']['sha256'])
    assert (z['scale'],z['V_grid'],z['degree'],z['arity'],z['cut'])==(A.S,A.G,12000,30,64)
    cap=z['cap_f_ticks']; s=A.Raw(side)
    assert s.lam==Q(z['lambda']) and s.sign==(1 if side=='lower' else -1)
    W={(d,c):Q(w,A.G) for d,c,w in z['values']}
    rows={(d,c):(Q(a),Q(b)) for d,c,a,b in z['rows']}
    cells=lambda H:{(d,c) for d in range(2,H+1) for c in range(1,min(d,30)+1)}
    assert len(z['values'])==len(W)==1514 and set(W)==cells(65)
    assert len(z['rows'])==len(rows)==1484 and set(rows)==cells(64)
    assert W[2,1]==Q(z['V_J_ticks'])==Q((525907178733179898,991668031155825523)[s.side],A.G)
    mini=None
    for (d,c),(ra,rb) in rows.items():
        old=W[d,c]; a,b=s.row(d,c); L=d-c
        bud=W.get((d,c+1),0); gain=W[d+1,c]
        for r in (Q(0),Q(c,2),Q(c,7)):
            expect=cap+a+b*r+(L+r)*(gain-old)+Q(L,d)*(bud-old)-hi*old
            assert ra+rb*r==expect>=0
            mini=expect if mini is None else min(mini,expect)
        # Direct Fraction quotient on the directed grid, independent of export's integer factorization.
        nums=Q(L,d)*bud
        pay=lambda r:(Q((cap+a)*A.G//1,A.G)+r*Q(b*A.G//1,A.G)+nums+(L+r)*gain)/(hi+Q(L,d)+L+r)
        assert old==Q(min(pay(0),pay(Q(c,2)))*A.G//1,A.G)
    counts=physical(s,A,W,rows,cap,hi)
    s.finish()
    print(side,'1514 values / 1484 exact rows / 4452 residual checks; minimum',mini,'physical',counts,flush=True)
    return rows,W[2,1]


def palm_check(a,lo,hi):
    assert lo<=hi<=a.HI
    K=a.class_kernel()
    assert len(K)==12 and all(len(r)==12 and all(x>=0 for x in r) and sum(r)<=1 for r in K)
    # Reconstruct by the discounted one-jump inequality at both endpoints;
    # last mother class uses p<=1/13, last reward class includes 13..34.
    for i,j in product(range(12),repeat=2):
        vals={}
        for d in reversed(range(2,35)):
            for k in reversed(range(d)):
                l=d-1-k; b=Q(l,d); nxt=vals.get((d+1,k),0)
                pay=b*((min(d,13)==j+2)+vals.get((d,k+1),0))
                def f(r):return (pay+r*nxt)/(a.HI+b+r)
                g=l+Q(k,2)+Q(1,i+2)
                value=Q(min(f(l),f(g))*a.G//1,a.G)
                assert 0<=value<=f(Q(l+g,2))
                vals[d,k]=value
        assert K[i][j]==vals[2,0]
    print('144 nonnegative substochastic Palm kernel entries; original kernel discount >= quintet HI',flush=True)
    return K


def price(mod,K,side,lo,hi,rows,V):
    from scipy.optimize import linprog
    mod.LO,mod.HI=lo,hi
    lp=mod.build(24,True,K); sols=[]
    assert (len(lp.rows),len(lp.cols))==(613,5554)
    def cost(kind,m):
        if kind!='occ' or m[0]>24:return Q(0)
        d,k,_,_,u,v=m; a,b=rows[d,k+1]
        return a+b*(u+v)
    for kind,col,m in lp.cols:
        z=dict(col)
        def get(key):return z.get(lp.rows[key],Q(0))
        if kind=='occ':
            d,k,b,p,u,v=m
            if 2<d<=24:
                target=(d+1,k if d<24 else min(7,k))
                assert get(('age',target))==-(d-k-1+u+v)
            if (d,k)==(25,7):
                assert u==0 and cost(kind,m)==0 and get(('up',(25,3)))==0
        if kind=='tail-up':
            assert get('rho')==get(('up',(25,3)))==1 and cost(kind,m)==0
    def capture(*a,**kw):
        s=linprog(*a,**kw); sols.append(s); return s
    with patch('scipy.optimize.linprog',capture):
        got=lp.lower(lambda k,m:cost(k,m)/10**7)*10**7
    q=[Q(round(x*10**12),10**5) for x in sols[-1].eqlin.marginals]
    defects=defaultdict(Q); n=Counter()
    for k,col,m in lp.cols:
        a=cost(k,m); assert a>=0
        defects[k]=max(defects[k],sum(q[i]*v for i,v in col)-a)
        n[k]+=1
    assert got==q[lp.rows['birth']]-sum(defects.values())/lo
    assert got*lo==PRICE[side]
    print(side,'LP columns',dict(n),'exact grouped repair (occupation ticks)',sum(defects.values())/lo,
          'stationary reduced ticks',mod.dec(lo*got),'full singleton ticks',mod.dec(lo*(V+got)),flush=True)


if __name__=='__main__':
    assert __debug__
    for f,h in PINS.items():assert sha256((R/f).read_bytes()).hexdigest()==h,f
    A=load('checks/quintet-lifetime-sol.py')
    raw=(A.P/'source/round2/neon2/quintet-cert/certificate_gold1.json').read_bytes()
    assert sha256(raw).hexdigest()==RECEIPT
    gold=json.loads(raw); lo,hi=(Q(gold[s]['certified_exact']) for s in PRICE)
    B=load('occupation-neon/pressure-memory/birth-law/replay.py')
    K=palm_check(load('checks/stationary-p4-sol.py'),lo,hi)
    for side in PRICE:
        rows,V=rows_check(side,A,lo,hi,gold[side])
        price(B,K,side,lo,hi,rows,V)
