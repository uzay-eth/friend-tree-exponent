"""Independent physical shifts, zero-box barrier, and all-power seam checks."""

from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from runpy import run_path
import sys


def barrier():
    rows = 0
    for k in (2,3,4,32,48,64,96,257):
        for j in range(k-1):
            f = Q(1,j+1)
            drift = (j+1)*(Q(1,j+2)-f)+Q(1,k)*(1-f)
            if j:
                drift += Q(j,k)*(Q(1,j)-f)
            assert drift == -Q(1,j+2)+(Q(1,k) if j else 0) <= 0
            for d in (k-1,k,k+1):
                if j<d:
                    smax=Q(d+j,2*d*(d+1))
                    assert smax/(j+1)<=Q(1,2*(d+1))
            rows += 1
    print(f'PASS {rows} exact first-send and stock-exposure rows, including lo=0')


def fixtures():
    rows=0
    for n in range(2,8):
        for parents in product(*(range(i) for i in range(1,n))):
            g=[set() for _ in range(n)]
            for i,v in enumerate(parents,1):
                g[i].add(v); g[v].add(i)
            d=list(map(len,g))
            r=[sum(Q(1,d[v]) for v in vs) for vs in g]
            s=[r[u]/(d[u]*(d[u]+1)) for u in range(n)]
            assert sum(r)==n
            for u,v in ((u,v) for u in range(n) for v in g[u] if v>u):
                h=[vs.copy() for vs in g]+[{u},{v}]
                h[u].add(n); h[v].add(n+1)
                way=r[u]*(r[v]-Q(1,d[u]*(d[u]+1)))
                other=r[v]*(r[u]-Q(1,d[v]*(d[v]+1)))
                assert way-other==s[v]-s[u]
                for x,y,tip in ((u,v,n),(v,u,n+1)):
                    stock=sum(len(g[z])==1 for z in g[x])-(d[y]==1)
                    actual=sum(z!=tip and len(h[z])==1 for z in h[x])
                    assert len(h[x])==d[x]+1 and actual==stock<=len(h[x])-2
                    if d[x]>=2 and d[y]>=2:
                        assert s[y]>=Q(1,d[x]*d[y]*(d[y]+1))
                    rows+=1
            cores=[u for u in range(n) if d[u]>1]
            retained=[u for u in cores if not (d[u]==2 and all(v in cores for v in g[u]))]
            if len(cores)>=2:
                c=[sum(v in cores for v in g[u]) for u in retained]
                assert len(retained)>=2 and sum(x-2 for x in c)==-2
                for a,b in ((Q('.015'),Q('.015')),(Q('.011'),Q('.011'))):
                    assert b*len(c)-(len(c)-1)*b-a*sum(x-2 for x in c)==b+2*a
    print(f'PASS {rows} oriented two-placement stock shifts and skeleton census fixtures')


def seams():
    parts=[(Q(0),Q(473,512)),(Q(473,512),Q(1)),(Q(1),Q(3)),(Q(3),Q(7))]
    assert all(a[1]==b[0] for a,b in zip(parts,parts[1:]))
    for n,p in product(range(2,17),range(1,10)):
        from math import factorial as f
        B=lambda n:Q(f(n-1)*f(p),f(n+p))
        assert (B(n)-B(n+1))/(p+1)==Q(f(p)*f(n-1),f(n+p+1))
    assert min(Q(8,11)/(1-Q(108,305)),Q(48,121)/Q(387,1000))>1
    assert Q(1892360505093,10**15)>Q(9,5000)
    assert Q(22890923,1058158080)>Q(21,1000)
    print('PASS closed exponent seams, 135 beta telescopes and fractional founder endpoints')


def stars():
    root=Path(__file__).resolve().parents[1]/'open/age-order/hub-ranking/mean-ledger/fractional'
    source=run_path(str(root/'check.py'))
    inner=source['FM']
    minima=[Q(1)]*3
    for j in range(512):
        lo,hi=Q(j,512),Q(j+1,512)
        if j<473:
            D=96 if j<36 else 64
            box=source['Bounds'](lo,hi,48 if j<36 else 32)
            if j<128:
                box.boost()
            g=box.m-box.U(2,0)
            costs=[]
            for d in range(2,D):
                centre,leaf=Q(d,d+1),Q(1,2)
                outgoing=centre-leaf
                value=(centre+leaf)*g+centre*box.bonuses.get(d+1,0)+leaf*box.bonuses.get(2,0)
                value-=outgoing*max(0,box.U(d+1,d-1)-box.L(2))
                assert value==box.profile(d,0,0)[1]
                costs.append(value)
            tail=Q(D,D+1)*g-Q(D+2,(D+1)**2)
            box.profile.cache_clear()
        else:
            D=32
            M=inner['middle'](lo,hi)
            U,L=inner['side_bounds'](lo,hi,M)
            g=M-U(2)
            costs=[]
            for d in range(2,D):
                centre,leaf=Q(d,d+1),Q(1,2)
                costs.append((centre+leaf)*g-(centre-leaf)*max(0,U(d+1)-L(2)))
            tail=Q(D,D+1)*(g-1/(lo*(D+1)))
        minima=[min(x,y) for x,y in zip(minima,(g,min(costs),tail))]
        inner['gamma_upper'].cache_clear()
    assert minima==[Q(406542520567,2500000000000),Q(2085176012479,12000000000000),Q(265041348327659,2012097656250000)]
    assert min(minima)>Q(9,5000)
    print(f'PASS 512 independent star formulas and seed/tail minima: {min(minima)}')


if __name__=='__main__':
    barrier()
    fixtures()
    seams()
    if '--stars' in sys.argv:
        stars()
