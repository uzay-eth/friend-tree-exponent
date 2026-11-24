#!/usr/bin/env python3
"""Same-mask Yt refinement: asymmetric side caps and physical birth mixing."""
import argparse
import json
from fractions import Fraction as Q
from functools import cache
from importlib.util import module_from_spec, spec_from_file_location
from math import ceil, floor, lcm
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load('fork_base', HERE.parent/'replay.py')
LOW, HIGH, GRID, THETA = B.LOW, B.HIGH, B.GRID, B.THETA
CAP, ARITY = 5, 3
PARENT = (-1,0,1,1)
VIRGIN = (2,0,())
TAIL = (0,0)
GROUPS = (2,3,5,10,20)


def birth_floors():
    # Recompute the physical (not residual/Perron) direct birth kernel.
    source = load('physical_birth_kernel', B.NEON/'renewal'/'next'/'replay.py')
    assert (source.LABELS, source.HI)==(GROUPS,HIGH)
    masses, _ = source.newborn_kernel(False)
    assert all(a>=Q(b) for a,b in zip(masses,('.20490922','.16627050','.13980026','.08556734','.02650492')))
    print('physical birth floors', *(B.decimal15(floor(a*10**15)) for a in masses), flush=True)
    return masses


def prepare(rows, diagnostic):
    domain = {}
    zero = (0,)*9
    names = [str(t) for t in THETA]
    if diagnostic:
        names += ['fork32/'+s for s in ('total','intercept','apex-parent','apex-children','tips-children')]
    for x in B.labels(rows[0].c,'Yt',CAP):
        if x[0][2]!=1 or any(k>ARITY for _,k,_ in x):
            continue
        a,b=(row.coefficients('Yt',x) for row in rows)
        q=[(a[0]-t*b[0],*(u-t*v for u,v in zip(a[1],b[1]))) for t in THETA]
        if diagnostic:
            u=q[1]
            q += [u, *(tuple(v if i in inds else 0 for i,v in enumerate(u))
                       for inds in ({0},{1},{2},{6,8}))] if x[1]==(0,2,1) else [zero]*5
        domain[x]=tuple(q)
    assert len(domain)==4095
    return names,domain


def solve(names,domain,caps,mom_cap,masses,budget,barrier=False,only32=False,forecast=False,sided=False,tip32=False):
    columns=len(names)+int(barrier)
    zr=(0,)*(2*columns)
    max_degree=max(GROUPS[-1],CAP,mom_cap+1,max(caps)+1)
    rr=lcm(10**6,*range(2,max_degree+1))
    inv=[0]+[rr//j for j in range(1,max_degree+1)]
    ds,du=int(LOW*rr),int((HIGH-LOW)*rr)
    # Include a two-sided Bellman column for the signed parent source under
    # a discrete antiderivative in apex degree. Optional, diagnostics only.
    if barrier:
        assert names[-1]=='fork32/tips-children'
    if only32:
        names=names[2:]
        columns-=2
        zr=(0,)*(2*columns)
    states=0
    visited=set()

    def ticks(q):
        return floor(q*GRID),ceil(q*GRID)

    def canon(path):
        return path if len(path)<4 or path[2]<=path[3] else path[:2]+(path[3],path[2])

    @cache
    def score(x):
        if x in domain:
            q=domain[x]; visited.add(x)
        else:
            y=x[:2]+(x[3],x[2]); visited.add(y)
            q=tuple(u[:-4]+u[-2:]+u[-4:-2] for u in domain[y])
        if barrier:
            q=q+(q[4],)
        if only32:
            q=q[2:]
        return tuple(tuple(w for u in q for w in ticks(u[i])) for i in range(9))

    @cache
    def phi(path):
        # Killed exact antiderivative: Phi(D+1)-Phi(D) cancels the signed
        # pressure slope. The post-cap potential is identically zero.
        if not barrier or len(path)!=4 or path[1][:2]!=(3,2):
            return Q(0)
        x=tuple((d-k-1,k,1) for d,k,_ in path)
        d,k,h=x[0]
        return sum((domain[( (j-k-1,k,1),)+x[1:]][4][1] if ((j-k-1,k,1),)+x[1:] in domain
                    else domain[((j-k-1,k,1),x[1],x[3],x[2])][4][1]
                for j in range(sum(x[0]),CAP+1)),Q(0))

    def future(path,mom):
        if any(d>CAP or k>ARITY for d,k,_ in path) or (only32 and len(path)>1 and path[1][0]>3):
            return zr
        if only32 and 1<len(path)<4 and path[1][1]>len(path)-2:
            return zr  # an extra fork child leaves too few core slots to complete
        if len(path)<4 and path[PARENT[len(path)]][1]==ARITY:
            return zr
        if tip32 and len(path)>2 and path[1][0]>3:
            path=path[:2]+tuple((d,k,(TAIL,)*len(s)) for d,k,s in path[2:])
        return rec(canon(path),mom)

    @cache
    def rec(path,mom):
        nonlocal states
        states+=1
        if states>budget:
            raise RuntimeError('retained state budget exceeded')
        n=len(path)
        credit=score(tuple((d-k-1,k,1) for d,k,_ in path)) if n==4 else (zr,)*9
        edges=[]
        pot=phi(path) if barrier else 0

        def put_edge(lo,hi,y,mark=mom,src=zr):
            if not hi:
                return
            val=future(y,mark)
            if barrier:
                nxt=phi(y) if all(d<=CAP and k<=ARITY for d,k,_ in y) else 0
                if src is credit[1]:
                    # Exact cancellation before directed rounding. Discount
                    # and all the other arrivals retain their drift terms.
                    src=src[:-2]+(0,0)
                elif nxt!=pot:
                    a,b=ticks(nxt-pot)
                    src=src[:-2]+(src[-2]+a,src[-1]+b)
            edges.append((lo,hi,tuple(a+b for a,b in zip(val,src)) if src!=zr else val))

        for i,(d,k,sides) in enumerate(path):
            cap=caps[i]
            if tip32 and i>=2 and path[1][0]>3:
                cap=0
            l=d-k-1
            ch=tuple(j for j in range(i+1,n) if PARENT[j]==i)
            assert k==len(ch)+len(sides)

            def put(d1,k1,s1):
                return path[:i]+((d1,k1,s1),)+path[i+1:]

            y=put(d+1,k,sides)
            fixed=l*rr+sum(inv[path[j][0]] for j in ch)+(inv[path[PARENT[i]][0]] if i else 0)
            put_edge(fixed,fixed,y)
            if not i:
                d0,k0=mom
                a=inv[d0] if d0>0 else 0
                b=a if d0>0 else inv[-d0]
                put_edge(a,b,y,src=credit[1])
            known=sum(inv[s[0]] for s in sides)
            put_edge(known,known+sides.count(TAIL)*inv[cap+1 if cap else 2],y,src=credit[2+2*i])
            if l:
                plain=put(d,k+1,tuple(sorted(sides+(((2,0) if cap else TAIL),))))
                if n<4 and i==PARENT[n]:
                    extend=put(d,k+1,sides)+(VIRGIN,)
                    v,w=future(plain,mom),future(extend,mom)
                    rew=ticks(phi(extend)) if barrier and all(dd<=CAP and kk<=ARITY for dd,kk,_ in extend) else (0,0)
                    value=tuple(a+b for a,b in zip(v,w))
                    value=value[:-2]+(value[-2]+rew[0],value[-1]+rew[1]) if barrier else value
                    edges.append((l*inv[d],l*inv[d],value))
                else:
                    put_edge(l*inv[d],l*inv[d],plain)
            for j,(dd,kk) in enumerate(sides):
                if not dd:
                    continue
                stock=dd-kk-1

                def change(s):
                    return put(d,k,tuple(sorted(sides[:j]+(s,)+sides[j+1:])))

                a=stock*rr+inv[d]
                put_edge(a,a+kk*inv[2],change((dd+1,kk) if dd<cap else TAIL))
                put_edge(stock*inv[dd],stock*inv[dd],change((dd,kk+1)))
        d,k=mom
        if d>0:
            l=d-k-1
            a=l*rr+inv[path[0][0]]
            b=a+k*inv[2]  # remaining k-1 children and the mother's parent
            put_edge(a,b,path,(d+1,k) if d<mom_cap else (-(mom_cap+1),0))
            put_edge(l*inv[d],l*inv[d],path,(d,k+1))

        base=ds+sum(a for a,b,v in edges)
        toll=credit[0]
        discount_credit=(0,0)
        if barrier:
            a,b=ticks(-LOW*pot)
            toll=toll[:-2]+(toll[-2]+a,toll[-1]+b)
            discount_credit=ticks(-pot)
        ans=[]
        for j,f in enumerate(toll):
            upper=j%2
            num,den=f*rr,base
            discount_y=discount_credit[upper] if barrier and j>=len(zr)-2 else 0
            terms=[(discount_y,du)]
            for a,b,v in edges:
                q=v[j]; num+=a*q
                if b>a:
                    terms.append((q,b-a))
            for q,r in sorted(terms,reverse=bool(upper)):
                if (q*den<=num if upper else q*den>=num):
                    break
                num+=q*r; den+=r
            x=(num+den-1)//den if upper else num//den
            dx=discount_y-x
            res=f*rr-ds*x+(du*dx if (dx>0)==bool(upper) else 0)
            res+=sum((b if ((v[j]-x>0)==bool(upper)) else a)*(v[j]-x) for a,b,v in edges)
            assert res<=0 if upper else res>=0
            ans.append(x)
        return tuple(ans)

    entries={}
    for d in range(2,mom_cap+1):
        entries[d]=tuple(future((VIRGIN,),(d,k)) for k in range(1,d))
    # The tail respects the actual permanent birth-bin lower degree, rather
    # than passing every unmonitored birth through the same artificial state.
    starts=[]
    for i,g in enumerate(GROUPS):
        stop=GROUPS[i+1] if i+1<len(GROUPS) else 21
        vals=[v for d in range(g,min(stop,mom_cap+1)) for v in entries[d]]
        if stop>mom_cap+1:
            vals.append(future((VIRGIN,),(-max(g,mom_cap+1,2),0)))
        starts.append(vals)
    expected={x for x in domain if not only32 or x[1]==(0,2,1)}
    assert visited==expected, (len(visited),len(expected))

    def extrema(values,j):
        return (min(v[2*j] for v in values),max(v[2*j+1] for v in values))

    residual=1-sum(masses)
    if sided:
        sfile=B.NEON/'pressure-memory'/'birth-law'/'replay.py'
        S=load('birth_law',sfile)
        A=load('p4_birth_kernel',B.NEON.parent/'checks'/'stationary-p4-sol.py')
        LP=S.build(24,True,A.class_kernel())
    print('caps',caps,'mother',mom_cap,'rows',len(visited),'states',states,
          'residuals',states*len(zr),flush=True)
    for j,name in enumerate(names+(['fork32/apex-parent-drift'] if barrier else [])):
        bins=[extrema(v,j) for v in starts]
        lo=min(v[0] for v in bins); hi=max(v[1] for v in bins)
        mix=(sum(w*v[0] for w,v in zip(masses,bins))+residual*lo,
             sum(w*v[1] for w,v in zip(masses,bins))+residual*hi)
        def norm(b):
            vals=[s*v/GRID for s in (LOW,HIGH) for v in b]
            return min(vals),max(vals)
        l,u=norm(mix); a,b=norm((lo,hi))
        fmt=B.decimal15
        print(f'{name}: [{fmt(floor(l*10**15))}, {fmt(ceil(u*10**15))}] '
              f'width≤{fmt(ceil((u-l)*10**15))} unweighted≤{fmt(ceil((b-a)*10**15))}',flush=True)
        if j in (1,2,4):
            print('  birth-bin width at high discount:',*(fmt(ceil(HIGH*(b-a)*10**15/GRID)) for a,b in bins),flush=True)
        if forecast:
            print('forecast',json.dumps(dict(name=name,grid=GRID,groups=GROUPS,bins=bins,entries=entries,
                                             caps=caps,mother=mom_cap)),flush=True)
        if sided and (name in ('1','9/10','fork32/total','fork32/apex-parent','fork32/apex-parent-drift')):
            def cv(d,k,side):
                if d<=mom_cap:
                    v=entries[d][min(max(1,k+1),d-1)-1]
                    return Q(v[2*j+side],GRID)
                else:
                    # LP tail degree cells preserve the five birth groups
                    # only after rounding down to their actual minimum.
                    i=max(i for i,g in enumerate(GROUPS) if g<=d)
                    return Q(bins[i][side],GRID)
            x,y=S.price(LP,lambda d,k:cv(d,k,0),lambda d,k:cv(d,k,1))
            x,y=max(x,l),min(y,u)
            assert x<=y
            print(f'  sided [{fmt(floor(x*10**15))}, {fmt(ceil(y*10**15))}] width≤{fmt(ceil((y-x)*10**15))}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--caps',nargs=3,type=int,default=(2,3,2),metavar=('APEX','FORK','TIP'))
    p.add_argument('--mother',type=int,default=3)
    p.add_argument('--budget',type=int,default=850000)
    p.add_argument('--diagnostic',action='store_true')
    p.add_argument('--barrier',action='store_true')
    p.add_argument('--only32',action='store_true')
    p.add_argument('--forecast',action='store_true')
    p.add_argument('--sided',action='store_true')
    p.add_argument('--tip32',action='store_true',help='spend tip memory only while fork degree is at most three')
    a=p.parse_args()
    if not all(d==0 or 2<=d<=9 for d in (*a.caps,a.mother)):
        p.error('caps must be zero or in [2,9]')
    lo,hi=Path('b3f18_lo.bin'),Path('b3f18_hi.bin')
    B.authenticate(B.NEON/'all-contrast'/'scan.py',B.SCAN_SHA)
    B.authenticate(lo,B.PAYLOAD_SHA['lower'])
    B.authenticate(hi,B.PAYLOAD_SHA['upper'])
    names,domain=prepare(B.pairs(lo,hi),a.diagnostic or a.barrier or a.only32)
    m=birth_floors()
    solve(names,domain,(*a.caps,a.caps[-1]),a.mother,m,a.budget,a.barrier,a.only32,a.forecast,a.sided,a.tip32)


if __name__=='__main__':
    main()
