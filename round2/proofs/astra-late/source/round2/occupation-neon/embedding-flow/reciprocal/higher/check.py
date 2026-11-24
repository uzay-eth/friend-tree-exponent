"""Bounded exact cross-check of original physical hV coefficients."""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import prod

from hv import Reader, Role, H, NB, bernstein, contrast, moment_controls, opener


def direct(c, d, k):
    # Independent numeric path: original sparse Cert embedding at both
    # physical degrees D and D+1, no symbolic shift or denominator.
    a,u,w=x=tuple((D-K-1,K,1) for D,K in zip(d,k))
    if k[1]>c.K:
        return {},{}
    def mixed(xs,fn):
        return sum(prod(wt for _,wt in row)*int(fn(*(v for v,_ in row)))
                   for row in product(*(c.embed(s) for s in xs)))
    def hv(z):
        return mixed(z,lambda a,b,w:c.hvtab[c.par[a],c.hv[b],c.child[w]] if b in c.hv else 0)
    full=hv(x)
    up=lambda s:(s[0]+1,s[1],s[2])
    bud=lambda s:(s[0]-1,s[1]+1,s[2])
    live={}
    t={'discount':-c.lam*full}
    routes=((0,'mother',1,a in c.p4 and u in c.p4),
            (0,'children',k[0]-1,a in c.ybu and u in c.ybw),
            (1,'children',k[1]-1,u in c.yt),(2,'children',k[2],u in c.p4 and w in c.p4))
    for i,s in enumerate(x):
        new=x[:i]+(up(s),)+x[i+1:]
        dp=hv(new)-full
        t[i,'leaves']=s[0]*dp
        t[i,'internal']=sum((Q(1,d[j]) for j in range(3) if abs(i-j)==1),Q(0))*dp
        for j,side,n,route in routes:
            if i!=j or not n: continue
            if route:
                assert d[i]<H
                t[i,side,'centered']=n*Q(int((c.rpar if side=='mother' else c.rch)[c.state[s]]),c.Sr)*dp
            else:
                live[i,side]=(n,Q(dp,c.S))
        if not s[0]: continue
        v,b,z=x[:i]+(bud(s),)+x[i+1:]
        birth=Q(s[0],d[i])
        t[i,'promotion']=birth*(hv((v,b,z))-full)
        nz=0
        if i==0 and v in c.ybu and b in c.ybw:
            nz=mixed((z,),lambda q:c.ybtab[c.ybu[v],c.ybw[b],c.child[NB],c.child[q]])
        if i==1 and b in c.yt:
            nz=mixed((v,z),lambda p,q:c.yttab[c.par[p],c.yt[b],c.tp(q,NB)])
        if i==2 and b in c.p4 and z in c.p4:
            nz=mixed((v,),lambda p:c.p4tab[c.par[p],c.p4[b],c.p4[z],c.child[NB]])
        t[i,'newborn']=birth*nz
    G=lambda p,w:int(c.g[c.state[p],c.state[w]])
    for is_apex,neighbor,donor in [(True,w,0),(False,a,2)]:
        if d[1]<H:
            center=Q(int((c.rpar if is_apex else c.rch)[c.state[u]]),c.Sr)
            fee=(mixed((up(u),neighbor),G)-mixed((u,neighbor),G) if is_apex else
                 mixed((neighbor,up(u)),G)-mixed((neighbor,u),G))
        else:
            center=Q(int((c.trpar if is_apex else c.trch)[k[1]]),c.Sr)
            fee=Q(0)
            for state,p in c.embed(u):
                if state[0]=='T2':
                    v=('T1',k[1],1)
                    if is_apex:
                        fee+=p*sum(q*Q(G(('inf',),z)-G(v,z),H) for z,q in c.embed(neighbor))
                    else:
                        fee+=p*sum(q*Q(G(z,('inf',))-G(z,v),H) for z,q in c.embed(neighbor))
                if state[0]=='T3':
                    fee+=p*sum(q*int(c.yp[k[1],c.child[z]] if is_apex else c.yw[c.par[z],k[1]])
                               for z,q in c.embed(neighbor))
        t['g→hV','from apex' if is_apex else 'from tip']=(Q(1,d[donor])-center)*fee
    return {j:Q(v,c.S) for j,v in t.items()},live


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lower')
    p.add_argument('upper')
    p.add_argument('--authenticate',action='store_true')
    a=p.parse_args()
    readers=[Reader(opener(path,side,authenticate=a.authenticate))
             for path,side in [(a.lower,'lower'),(a.upper,'upper')]]
    cells=[((5,4,5),(3,2,2)),((17,17,17),(2,2,1)),((4,4,4),(1,2,1)),
           ((None,4,4),(2,2,1)),((4,None,4),(2,2,1)),((4,4,None),(2,2,1)),
           ((None,None,4),(2,2,1)),((None,4,None),(19,2,18)),
           ((None,None,None),(2,2,1)),((None,None,None),(1,17,0)),
           ((None,None,None),(19,18,18)),((18-1,None,18-1),(1,1,0)),
           ((3,3,2),(2,2,1))]
    counts={'rows':0,'terms':0,'ports':0,'extract':0,'controls':0,'zeros':0,'moment':0}
    nonzero=set()
    pole=None
    for ds,k in cells:
        roles=tuple(Role(d,j) for d,j in zip(ds,k))
        both=[r.row(roles) for r in readers]
        paired=contrast(*readers,roles,Q(7,5))
        degrees=tuple(7 if d is None else 0 for d in ds)
        bounds=[(Q(0),min(Q(1),Q(H,j+1))) for j in k]
        for row in [*both,paired]:
            assert all(row.base.degree(y)<=deg for y,deg in zip(row.base.gens,degrees))
            slopes=tuple(6 if d is None else 0 for d in ds)
            assert all(f.degree(y)<=deg for _,f in row.ports.values() for y,deg in zip(row.base.gens,slopes))
        # One tensor conversion per cell; verify both the exact power and
        # Bernstein evaluations at an interior rational control point.
        B=bernstein(paired.base,bounds,degrees)
        BD=bernstein(paired.denominator,bounds,degrees)
        assert min(BD.values())>0
        from math import comb
        v=Q(2,5)
        weights={j:prod(comb(n,i)*v**i*(1-v)**(n-i) for i,n in zip(j,degrees)) for j in B}
        point=tuple(v*h for _,h in bounds)
        assert sum(B[j]*w for j,w in weights.items())==Q(paired.base.eval(point))
        assert sum(BD[j]*w for j,w in weights.items())==Q(paired.denominator.eval(point))
        counts['controls']+=len(B)
        if ds==(None,4,4):
            f=paired.ports[0,'mother'][1]
            pole=Q(f.eval((-H,0,0)))
            assert pole
        if ds in [(None,4,4),(None,None,None)] and k==(2,2,1):
            axis=tuple(8 if D is None else 0 for D in ds)
            mass,fs=moment_controls(paired,3,bounds)
            W=Q(paired.denominator.eval(point))
            def moment(coeff):
                return sum(x*prod(comb(n,i)*v**i*(1-v)**(n-i) for i,n in zip(j,axis))
                           for j,x in coeff.items())/W
            assert moment(mass)==Q(paired.base.eval(point))/W
            for j,(m,g) in fs.items():
                assert moment(g)==Q(paired.ports[j][1].eval(point))/W
                counts['moment']+=1
            counts['moment']+=1
        for n in (18,19,53,10**9+7) if None in ds else (18,):
            d=tuple(max(n+i*(n!=18),k[i]+1) if D is None else D for i,D in enumerate(ds))
            point=tuple(Q(H,D) if v is None else Q(0) for v,D in zip(ds,d))
            for r,row in zip(readers,both):
                t,ports=direct(r.c,d,k)
                den=Q(row.denominator.eval(point))
                def ev(f): return Q(f.eval(point))/den
                assert {j:ev(f) for j,f in row.terms.items() if ev(f)}=={j:v for j,v in t.items() if v},(d,k)
                assert {j:(n,ev(f)) for j,(n,f) in row.ports.items()}==ports
                nonzero.update(j for j,z in t.items() if z)
                live={j:Q(m,11+i) for i,(j,(m,_)) in enumerate(ports.items())}
                assert row.at(d,live)==sum(t.values())+sum(live[j]*f for j,(_,f) in ports.items())
                if max(d)<H-1:
                    x=tuple((di-ki-1,ki,1) for di,ki in zip(d,k))
                    old,extra=r.raw.row('hV',x)
                    assert ev(row.base)==r.raw.sign*old/r.c.S
                    assert ports=={(i,'mother' if par else 'children'):(m,Q(r.raw.sign*w,r.c.S))
                                   for i,par,m,w in extra}
                    counts['extract']+=1
                counts['rows']+=1
                counts['terms']+=len(t)
                counts['ports']+=len(ports)
                counts['zeros']+=not t
            assert paired.at(d,live)==both[0].at(d,live)+Q(7,5)*both[1].at(d,live)
    assert all(k in nonzero for k in [(0,'newborn'),(1,'newborn'),(2,'newborn'),(0,'mother','centered'),
                                     (0,'children','centered'),(1,'children','centered'),(2,'children','centered')])
    print('exact fixtures',counts,'nonzero branches',len(nonzero),'genuine cubic pole numerator at y0=-18',pole)


if __name__=='__main__':
    main()
