#!/usr/bin/env python3
"""Joint physical mother reciprocal during the live D18 Yt (3,2) fork."""
import argparse
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
import importlib.util
from itertools import product
from math import ceil, floor, lcm
from pathlib import Path

HERE=Path(__file__).resolve().parent
FORK=HERE.parents[1]/'fork-lifetime'
PIN='db0fba60406930918139d59731e6cac545bc877955508aa8f819967def980999'
GRID=10**16
LO,HI=Q('.565673'),Q('.565775')
NAMES=('full','mother','mass','mass*p')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def raw():
    path=FORK/'replay.py'
    assert sha256(path.read_bytes()).hexdigest()==PIN
    m=load('fork_reference',path)
    assert sha256((m.NEON/'all-contrast'/'scan.py').read_bytes()).hexdigest()==m.SCAN_SHA
    rows=m.pairs(Path('b3f18_lo.bin'),Path('b3f18_hi.bin'))
    domain={}
    for z in m.labels(rows[0].c,'Yt',5):
        if z[0][2]!=1 or any(k>3 for _,k,_ in z) or (sum(z[1]),z[1][1])!=(3,2):
            continue
        a,b=(r.coefficients('Yt',z) for r in rows)
        base=a[0]-b[0]
        c=tuple(x-y for x,y in zip(a[1],b[1]))
        assert c[2:4]==(0,0)
        domain[z]=base,c
    assert domain
    return domain


def bellman(initial, edges, ds, du, nr):
    den0=ds+sum(lo for lo,_,_,_ in edges)
    result=[]
    for j,f in enumerate(initial):
        out=[]
        for side in (0,1):
            num,den=f[side]*nr,den0
            optional=[(0,du)]
            for lo,hi,vals,credit in edges:
                z=vals[j][side]+credit[j][side]
                num+=lo*z
                if hi>lo:optional.append((z,hi-lo))
            for z,r in sorted(optional,reverse=bool(side)):
                if (z*den<=num if side else z*den>=num):break
                num+=r*z
                den+=r
            v=-(-num//den) if side else num//den
            rho=f[side]*nr-(ds+du*int((v<0)==bool(side)))*v
            for lo,hi,vals,credit in edges:
                z=vals[j][side]+credit[j][side]-v
                rho+=(hi if (z>0)==bool(side) else lo)*z
            assert rho<=0 if side else rho>=0
            out.append(v)
        assert out[0]<=out[1]
        result.append(tuple(out))
    return tuple(result)


def occupations(domain, mother, arity=False, channels=NAMES, tilt=Q(0)):
    """Mom controls: coarse degree or joint (degree, number of other cores)."""
    nr=lcm(10**6,tilt.denominator,*range(2,max(5,mother+1)+1))
    rate=[0]+[nr//d for d in range(1,max(5,mother+1)+1)]
    ds,du=int((LO+tilt)*nr),int((HI-LO)*nr)
    zero=((0,0),)*len(channels)
    stop=((5,3),(3,2),(5,3),(5,3))
    # 0->1 stem, 1->2 older tip, 1->3 younger tip.
    parents=(0,0,1,1)
    hit=set()

    def ticks(x): return floor(x*GRID),ceil(x*GRID)

    @cache
    def toll(path):
        z=tuple((d-k-1,k,1) for d,k in path)
        ordering=lambda z:(z[1],z[0],z[2])
        if ordering(z[-2])>ordering(z[-1]):
            canon=z[:-2]+(z[-1],z[-2]); swap=True
        else:
            canon=z;swap=False
        a,c=domain[canon]
        hit.add(canon)
        if swap:c=c[:-4]+c[-2:]+c[-4:-2]
        initial=tuple(ticks(a if j=='full' else 1 if j=='mass' else 0) for j in channels)
        coeff=tuple(tuple(ticks(c[i] if j=='full' or j=='mother' and i==0
                                else 1 if j=='mass*p' and i==0 else 0)
                          for j in channels) for i in range(8))
        return initial,coeff

    @cache
    def V(path,m,a):
        n=len(path)
        if any(d>D or k>K for (d,k),(D,K) in zip(path,stop)):
            return zero
        if n<4 and (path[parents[n]][1]==stop[parents[n]][1]
                    or n>=2 and path[1][1]>n-2):
            return zero
        on=n==4 and path[1]==(3,2)
        initial,coeff=toll(path) if on else (zero,(zero,)*8)
        edges=[]
        def edge(l,h,z,c=zero):
            assert l<=h
            if h:edges.append((l,h,z,c))
        for i,(d,k) in enumerate(path):
            stock=d-k-1
            kids=sum(parents[j]==i for j in range(i+1,n))
            assert stock>=0 and k>=kids
            fixed=stock*nr+sum(rate[path[j][0]] for j in range(i+1,n) if parents[j]==i)
            if i:fixed+=rate[path[parents[i]][0]]
            gain=V(path[:i]+((d+1,k),)+path[i+1:],m,a)
            edge(fixed,fixed,gain)
            if not i:edge(rate[m],rate[m] if m else rate[mother+1],gain,coeff[0])
            # Any unregistered fork child would irreversibly violate (3,2)
            # before completion. Other roles use the matched opaque ports.
            if i==1:assert k==kids
            edge(0,(k-kids)*rate[2],gain,coeff[2*i+1])
            if stock:
                z=path[:i]+((d,k+1),)+path[i+1:]
                u=V(z,m,a)
                if n<4 and i==parents[n]:
                    w=V(z+((2,0),),m,a)
                    u=tuple((l+p,h+q) for (l,h),(p,q) in zip(u,w))
                edge(stock*rate[d],stock*rate[d],u)
        if m:
            nxt=m+1 if m<mother else 0
            own=rate[path[0][0]]
            if arity:
                stock=m-a-2
                assert stock>=0
                # Actual mother: ancestor + tagged apex + a other cores.
                edge(own+stock*nr,own+stock*nr+(a+1)*rate[2],V(path,nxt,a if nxt else -1))
                edge(stock*rate[m],stock*rate[m],V(path,m,a+1) if stock else zero)
            else:
                edge(own,own+m*nr-3*rate[2],V(path,nxt,-1))
        return bellman(initial,edges,ds,du,nr)

    keys=[(0,-1)]+[(m,k) for m in range(2,mother+1)
                   for k in (range(m-1) if arity else (-1,))]
    entries={x:V(((2,0),),*x) for x in keys}
    assert hit==domain.keys(),(len(hit),len(domain))
    print('Bellman',mother,'joint mother stock' if arity else 'degree only',
          'rows',len(domain),'states',V.cache_info().currsize,flush=True)
    V.cache_clear()
    return entries


def instant(ent,channels=NAMES):
    return {j:(min(s*Q(v[i][0],GRID) for v in ent.values() for s in (LO,HI)),
               max(s*Q(v[i][1],GRID) for v in ent.values() for s in (LO,HI)))
            for i,j in enumerate(channels)}


def pooled(bl,lp,ent,mother,joint,channels=NAMES):
    def row(d,k):
        if joint and d<=mother and k==d-1:
            # Bare mothers never give birth, so this column has zero b.
            return ent[0,-1]
        return ent.get((d,k if joint else -1),ent[0,-1])
    # Well-scaled small integer magnitudes for HiGHS, then exact verification.
    return {j: price_channel(lp,row,i,j) for i,j in enumerate(channels)}


def price_channel(lp,row,i,name):
    unit=10**6 if name in ('full','mother') else 10**3
    l=lp.lower(lambda k,m:unit*m[2]*Q(row(*m[:2])[i][0],GRID) if k=='occ' else 0)/unit
    u=-lp.lower(lambda k,m:-unit*m[2]*Q(row(*m[:2])[i][1],GRID) if k=='occ' else 0)/unit
    ends=[s*x for s,x in product((LO,HI),(l,u))]
    return min(ends),max(ends)


def past_discount_tail(lam,T):
    """Upper e^-lam*T by the reciprocal of a positive exponential sum."""
    x=lam*T
    a=p=Q(1)
    for k in range(1,40):
        p*=x/k
        a+=p
    return 1/a


def age_test(bl,domain,box,m,joint):
    # Discount tilting pays only the selected live crowns, with no bound by
    # the age of an exterior tree or by an unconditioned mean parent law.
    a=Q(1,4)
    ent=occupations(domain,m,joint,('mass*p',),-a)
    upper=instant(ent,('mass*p',))['mass*p'][1]
    for t in (4,8,12,16):
        tail=min(box[1],upper*past_discount_tail(a,t))
        young=(max(Q(0),box[0]-tail),box[1])
        print('apex age >',t,'moment ceiling',bl.dec(tail,True),
              'apex age <=',t,'moment',bl.fmt(young),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mother',type=int,default=6)
    p.add_argument('--limit',type=int,default=18)
    p.add_argument('--channels',nargs='+',choices=NAMES,default=NAMES)
    p.add_argument('--no-pool',action='store_true')
    p.add_argument('--age',action='store_true')
    p.add_argument('--birth-bias',action='store_true')
    p.add_argument('--variant',choices=('both','degree','stock'),default='both')
    p.add_argument('--skip-baseline',action='store_true')
    a=p.parse_args()
    assert 2<=a.mother<=12 and a.limit in (18,24)
    bl=load('birth_law',HERE.parent/'replay.py')
    ref=load('p4_audit',bl.R2/'checks'/'stationary-p4-sol.py')
    lp=bl.build(a.limit,True,ref.class_kernel())
    d=raw()
    channels=a.channels
    jobs=[('baseline',2,False)] if not a.skip_baseline else []
    if a.variant!='stock':jobs.append(('mother-degree',a.mother,False))
    if a.variant!='degree':jobs.append(('mother-stock',a.mother,True))
    for label,m,joint in jobs:
        entries=occupations(d,m,joint,channels)
        uniform=instant(entries,channels)
        got=uniform.copy() if a.no_pool else pooled(bl,lp,entries,m,joint,channels)
        for k in channels:
            z=got[k]=max(got[k][0],uniform[k][0]),min(got[k][1],uniform[k][1])
            assert z[0]<=z[1]
            print(label,k,'uniform',bl.fmt(uniform[k]),'bound',bl.fmt(z),
                  'width',bl.dec(z[1]-z[0],True),flush=True)
        if 'mass' in got and 'mass*p' in got and got['mass'][0]>0:
            a0,b0=got['mass*p']; c,e=got['mass']
            print(label,'incidence-weighted actual reciprocal',bl.fmt((max(0,a0)/e,min(Q(1,2),b0/c))),flush=True)
            if a.birth_bias and joint and not a.no_pool:
                i,j=channels.index('mass'),channels.index('mass*p')
                def bias(d,k):
                    z=entries.get((d,k),entries[0,-1])
                    if (d,k) not in entries:return ((0,ceil(Q(z[i][1],m+1))),)
                    return ((max(0,floor(Q(z[i][0],d)-z[j][1])),
                             ceil(Q(z[i][1],d)-z[j][0])),)
                loss=price_channel(lp,bias,0,'loss')
                print(label,'joint discounted loss versus same mothers at birth',
                      bl.fmt((max(0,loss[0]),loss[1])),flush=True)
    if a.age and 'mass*p' in got:age_test(bl,d,got['mass*p'],m,joint)


if __name__=='__main__':main()
