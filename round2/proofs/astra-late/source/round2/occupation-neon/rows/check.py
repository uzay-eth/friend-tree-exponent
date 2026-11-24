"""Standard-library replay of all L for the supported k bands; no giant files."""
import gzip
import json
import hashlib
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path

HERE=Path(__file__).resolve().parent


def f(k,L,v):
    return v['f'][k][L]


def embed(vals,x):
    return sum(y*comb(3,i)*x**i*(1-x)**(3-i) for i,y in enumerate(vals))


def table(v, mutate=False):
    cap=v['LF']; edge=v['D0']; scale=v['S']; sign=v['sign']; lam=Q(v['lam'])
    nb=f(0,1,v); out={}
    for k in range(5):
        floor=max(0,1-k); best=None; tail=None
        rp,rc,c2,c3=v['tr'][k]; center=Q(rp+k*rc,v['Sr']); c=k+1
        for L in range(floor,cap):
            d=L+c; prev=f(k,L,v); gain=f(k,L+1,v)-prev
            r0=L*gain-lam*prev
            if L:
                newborn=v['e'][k][L-1] if d<edge else embed(v['t'][k],Q(edge,d))
                r0+=Q(L,d)*(f(k+1,L-1,v)+nb+newborn-prev)
            if d<edge:
                p,q=v['r'][k][L-floor]; rows=[r0+Q(p+k*q,v['Sr'])*gain]
            else:
                correction=Q(c2*d+c3,d**3)
                if mutate: correction=-correction
                rows=[r0+center*correction+rho*(gain-correction) for rho in (0,Q(c,2))]
            signed=min(sign*x/scale for x in rows)
            if best is None or signed<best[0]: best=(signed,L)
            if d>=edge and (tail is None or signed<tail[0]): tail=(signed,L)
        d=cap+c; a=Q(c2,d*d); b=Q(c3,d**3); prev=f(k,cap,v)
        external=[r*x for r,x in product((center,center-Q(c,2)),(0,a/3,a+b))]
        atcap=[-lam*prev+Q(cap,d)*(f(k+1,cap-1,v)+nb+embed(v['t'][k],Q(edge,d))-prev)+e for e in external]
        # Restricted Bernstein coefficients: j independent draws with success
        # probability t; hence the j-th coefficient is a binomial-j average.
        t=Q(edge,d+1)
        restricted=[sum(comb(j,i)*t**i*(1-t)**(j-i)*v['t'][k][i] for i in range(j+1)) for j in range(4)]
        over=[-lam*prev+prob*(f(k+1,cap,v)+nb+x-prev)+e for prob,x,e in
              product((Q(cap+1,d+1),1),restricted,external)]
        out[k]=(best,min(sign*x/scale for x in atcap),min(sign*x/scale for x in over),tail)
    return out


def envelope(a):
    return min(a[0][0],a[1],a[2])


if __name__=='__main__':
    from check_edges import rows
    with gzip.open(HERE/'pinned-fg-slabs.json.gz','rt') as f0: d=json.load(f0)
    for name in ('lower','upper'):
        raw=d[name]['f']
        assert (raw['D0'],raw['LF'],raw['S'],raw['Sr'])==(18,14000,10**11,10**9)
        assert len(raw['f'])==6 and all(len(x)==14001 for x in raw['f'])
        assert [len(x) for x in raw['r']]==[16,16,15,14,13]
        assert [len(x) for x in raw['e']]==[16,15,14,13,12]
        assert len(raw['t'])==len(raw['tr'])==5 and all(len(x)==4 for x in raw['t']+raw['tr'])
        for kind in ('f','g'):
            packed=json.dumps(d[name][kind],sort_keys=True,separators=(',',':')).encode()
            print(name,kind,'canonical SHA-256',hashlib.sha256(packed).hexdigest(),flush=True)
        r=table(d[name]['f']); bad=table(d[name]['f'],True)
        for k,a in r.items():
            assert envelope(a)>0
            print(name,k,a[0][1],str(envelope(a)),f'{float(envelope(a)):.13g}',flush=True)
        a,b=(envelope(r[k]) for k in (0,1))
        assert a>2*b and envelope(bad[0])<0
        print(f'{name}: affine all-tree floor {b}*(C-1); planted reversed f-tail transfer gives k=0 minimum {envelope(bad[0])}',flush=True)
        g=min(rows(d[name]['g']).values())
        assert a-b>=g and 5*(a-b)-a>=6*g
        for k,(_,cap1,cap2,tail) in r.items():
            surplus=min(cap1,cap2,tail[0])-(a-(a-b)*k)
            assert surplus>g*(k+1)
            print(name,'tail coverage',k,'min surplus',surplus)
        C0=Q({'lower':184124847315,'upper':195820675182}[name],10**11)
        z=Q(d[name]['f']['lam'])+d[name]['f']['sign']*(b+g)/C0
        print(f'{name}: all-tree combined floor {b+g}*(C-1); endpoint={z} ({float(z):.16g}); incremental g shift={float(g/C0):.10g}',flush=True)
