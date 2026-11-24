"""Tagged Palm/adjacency Bellman bound: off-motif branches stay alive."""
import argparse
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb
from pathlib import Path
from rows import Rows, Cert, SPEC, label, up, bud

LO, HI=Q('.565673570322'),Q('.565775')
GRID=10**22


def down(x):
    return Q(x*GRID//1,GRID)


def tail_power(E,f,x,entry):
    """A positive lower polynomial face, in t=entry/D; zeros fail closed."""
    k,h=x[0][1:]
    n=7
    ts=[Q(entry,entry+j) for j in range(n+1)]

    def vals(d):
        t=Q(entry,d)
        z=label(x,0,(d-k-h,k,h))
        a,cols,ex=E.controls(f,z)
        # Polynomial endpoint controls kept separately (min is not polynomial).
        return [(E.c.D+Q(E.c.D,d))**3*(a+sum(v for _,v in r)+ex*p)
                for r in product(*cols) for p in (0,Q(1,2))]

    bases=[]
    for i,t in enumerate(ts):
        cp=[Q(1)]
        for j,u in enumerate(ts):
            if i==j:continue
            z=[Q(0)]*(1+len(cp))
            for l,a in enumerate(cp):
                z[l]-=a*u/(t-u)
                z[l+1]+=a/(t-u)
            cp=z
        bases.append(cp)
    pols=[[sum(z*b[i] for z,b in zip(ys,bases)) for i in range(n+1)]
          for ys in zip(*(vals(entry+j) for j in range(n+1)))]
    for d in (entry+11,entry+29):
        assert [sum(v*Q(entry,d)**i for i,v in enumerate(p)) for p in pols]==vals(d)
    # Lower-bound each control after factoring its own first nonzero power.
    q=0;mins=[]
    for p in pols:
        a=next((i for i,v in enumerate(p) if v),None)
        if a is None:return Q(0),0
        m=n-a
        low=min(sum(p[a+l]*Q(comb(j,l),comb(m,l)) for l in range(j+1)) for j in range(m+1))
        if low<=0:return Q(0),0
        q=max(q,a)
        mins.append(low)
    return min(mins)/(E.c.D+Q(E.c.D,entry))**3,q


def solve(E,f,degree,arity,tip,tail=False,start=18,input_cap=Q(1,2)):
    topo=SPEC[f]
    n=len(topo)
    count={}
    reach=[0,0]

    @cache
    def face(x):
        c,a=tail_power(E,f,x,degree+1)
        if not c:return Q(0)
        reach[0]+=1
        reach[1]=max(reach[1],a)
        # Stop on any tracked non-apex gain or any tracked promotion.
        # The remaining apex's actual gain is at most D at every age.
        hazard=1+sum(sum(z)+1 for z in x[1:])
        return down(c/(HI+hazard+a))

    @cache
    def V(x):
        m=len(x)
        d=sum(x[0]); k=x[0][1]
        if k>arity or any(sum(z)>tip for z in x[1:]):return Q(0)
        paid = m==n and k>=5 and all(z!=(1,0,1) for z in x[1:])
        if d>degree:
            return face(x) if tail and paid and d==degree+1 else Q(0)
        p=topo[:m]
        c,cols,ext=E.controls(f,x) if paid and d>=start else (Q(0),None,Q(0))
        if cols is not None:
            count[d]=count.get(d,0)+1
        const=c;den=HI
        todo=[]
        for j,z in enumerate(x):
            L,k,h=z; D=sum(z)
            ns=[i for i in range(m) if p[i]==j]+([p[j]] if j else [])
            g=L+sum(Q(1,sum(x[i])) for i in ns)
            b=Q(L,D)
            gain=V(label(x,j,up(z)))
            den+=g+b
            const+=g*gain
            if b:
                bz=label(x,j,bud(z))
                const+=b*V(bz)
                # The existing prefix survives. A distinct adjacent promotion
                # creates one new embedding with this prefix as owner.
                if m<n and j==topo[m]:const+=b*V(bz+((1,0,1),))
            if cols is None:
                nc=k-p.count(j)
                col=[(0,0),(Q(nc,2),0)] if nc else [(0,0)]
            else:col=cols[j]
            todo.append((gain,col))
        # Each independent two-corner control changes its minimizing corner
        # at most once as the trial value decreases from the first policy value.
        todo.append((V(label(x,0,up(x[0]))),[(0,0),(input_cap,ext*input_cap)]))
        v=Q(0)
        for _ in range(1+len(todo)):
            picked=[min(cols,key=lambda pair: pair[0]*(gain-v)+pair[1]) for gain,cols in todo]
            new=(const+sum(a*gain+w for (gain,_),(a,w) in zip(todo,picked)))/(den+sum(a for a,_ in picked))
            if new==v:break
            v=new
        # Direct dual check; floor at zero is valid by accepted nonnegative rows.
        v=max(Q(0),down(v))
        if v:
            assert den*v <= const+sum(min(a*(gain-v)+w for a,w in col) for gain,col in todo)
        return v

    out=LO*V(((1,0,1),))/(2 if f=='Yt' else 1)
    print(f,'pi lower',float(out),str(out),'states',V.cache_info().currsize,
          'paid descriptors',sum(count.values()),'tail tested/positive/max power',
          (face.cache_info().currsize,*reach),flush=True)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('blob',type=Path)
    p.add_argument('--family',choices=('hV','Yt','P4'),action='append')
    p.add_argument('--degree',type=int,default=40)
    p.add_argument('--arity',type=int,default=6)
    p.add_argument('--tip',type=int,default=3)
    p.add_argument('--start',type=int,default=18)
    p.add_argument('--tail',action='store_true')
    p.add_argument('--tail-only',action='store_true')
    p.add_argument('--palm',action='store_true',help='group the true capped newborn Palm labels')
    p.add_argument('--palm-group',type=int,action='append',choices=(3,5,10,20))
    a=p.parse_args()
    assert 18<=a.start<=a.degree and a.tip>=2 and 5<=a.arity<=17
    c=Cert(a.blob)
    assert (c.D,c.K,c.S,c.lam)==(18,17,10**11,Q('.565775' if c.mode==1 else '.565673'))
    e=Rows(c)
    def price(f):
        fn=lambda m: solve(e,f,a.degree,a.arity,a.tip,a.tail or a.tail_only,
                           a.degree+1 if a.tail_only else a.start,Q(1,m))
        v=fn(2)
        if not a.palm and not a.palm_group:return v
        # Directed floors on disjoint physical birth classes from next/replay.
        # All younger/new degrees retain the exact parent upper bound.
        weights=((3,Q('.16627050')),(5,Q('.13980026')),
                 (10,Q('.08556734')),(20,Q('.02650492')))
        total=v+sum(w*max(0,fn(i)-v) for i,w in weights if not a.palm_group or i in a.palm_group)
        print(f,'physical Palm-group price',float(total),str(total),flush=True)
        return total
    v=sum(map(price,a.family or ('hV',)))
    side='upper' if c.mode==1 else 'lower'
    den=Q('1.367169995210' if c.mode==1 else '1.356833021311')
    print(side,'increment',float(v/den),flush=True)
    if not a.tail and not a.tail_only and a.degree<=64:
        base=Q('.565771874464' if c.mode==1 else '.565676634491')
        y=(base+e.sign*v/den)*10**12
        ticks=e.sign*((e.sign*y)//1)
        print('cumulative outward',f'{ticks//10**12}.{ticks%10**12:012d}',flush=True)
