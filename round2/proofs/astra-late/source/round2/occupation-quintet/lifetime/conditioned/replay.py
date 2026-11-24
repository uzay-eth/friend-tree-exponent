#!/usr/bin/env python3
"""Original f lifetime with a live mother and latest real child."""
import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path
from itertools import product

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('quintet_life',HERE.parent/'check.py')
base=module_from_spec(spec);spec.loader.exec_module(base)
G=base.GRID


def rnd(n,d,up):
    return -(-n//d) if up else n//d


def coefficients(s,d,c,up):
    a,b=s.row_ticks(d,c); a+=base.GOLD[s.side]['cap_f']
    assert min(a,a+b*Q(c,2))>=0
    return tuple(rnd(t.numerator*G,t.denominator,up) for t in (a,b))


def broad(s,H,C,keep,up):
    """Same unconditioned grid solve as the original, retaining early slabs."""
    nxt=base.tail_suffix(s,H+1,C,upper=up)
    discount=base.LO if up else base.HI
    q,p=discount.denominator,discount.numerator
    saved={}
    for d in range(H,1,-1):
        cur=[0]*(C+2); w=2*d*q
        for c in range(min(d,C),0,-1):
            a,b=coefficients(s,d,c,up); L=d-c
            n=(a+L*nxt[c])*w+2*L*q*cur[c+1]
            k=2*d*p+2*L*q+L*w
            dx=c*d*q
            tips=(rnd(n,k,up),rnd(n+dx*(b+nxt[c]),k+dx,up))
            cur[c]=(max if up else min)(tips)
        if d<=keep: saved[d]=cur
        nxt=cur
    if (H,C)==(12000,30) and not up:
        assert saved[2][1]=={'lower':525907178733179898,'upper':991668031155825523}[s.side]
    return saved


def extremum(n,k,ports,up):
    # Port slopes are next value, plus reward slope for focal receive. The
    # optimal quotient activates the rates in slope order until no improvement.
    for z,r in sorted(ports,reverse=up):
        if (z*k>n if up else z*k<n):
            n+=z*r; k+=r
        else: break
    return rnd(n,k,up)


def conditional(s,W,D,K,M,T,up):
    # Mother exact 2..M; negative m denotes an unobserved degree >= -m.
    # Exact m advances monotonically and may enter the >=M+1 bucket.
    moms=tuple(dict.fromkeys((-3,-(M+1),*range(M,1,-1))))
    tips=[(0,0)]; tips += ([(T+1,0)]+[(j,i) for j in range(T,1,-1) for i in range(j,0,-1)]) if T else []
    def cells(c): return tips[:1] if c==1 or not T else tips[1:]
    newborn=(2,1) if T else (0,0)
    discount=base.LO if up else base.HI; pc,qc=discount.numerator,discount.denominator
    lookup={}
    checks=0; states=0
    def read(d,c,m,t):
        if d>D or c>K: return W[d][c]
        return lookup[d,c,m,t]

    for d in range(D,1,-1):
        for c in range(min(K,d),0,-1):
            L=d-c; a,z=coefficients(s,d,c,up)
            for m in moms:
                for t in cells(c):
                    j,i=t; live=m>0; child=i>0
                    Wg=read(d+1,c,m,t); Wb=read(d,c+1,m,newborn) if L else 0
                    den0=lcm(d,abs(m),j or 1,2)
                    den=den0*qc
                    D1=den//d
                    known=(den//m if live else 0)+(den//j if child else 0)
                    nx=den*L+known
                    n=a*den+z*known+nx*Wg+L*D1*Wb
                    k=den0*pc+nx+L*D1
                    rem=(c-1-int(j>0))* (den//2)+(0 if live else den//(-m))
                    if j and not child: rem+=den//j
                    ports=[(z+Wg,rem)]
                    if live:
                        u=read(d,c,m+1 if m<M else -(M+1),t)
                        n+=D1*u;k+=D1
                        # The mother already has the focal and her own parent;
                        # at birth both are cores. Her other slots count at most 1.
                        ports.append((u,(2*m-3)*(den//2)))
                    if child:
                        u=read(d,c,m,(j+1,i) if j<T else (T+1,0))
                        gain=(j-i)*den+D1
                        n+=gain*u;k+=gain
                        ports.append((u,(i-1)*(den//2)))
                        if j>i:
                            b=(j-i)*(den//j)
                            n+=b*read(d,c,m,(j,i+1));k+=b
                    v=extremum(n,k,ports,up)
                    if (d in (2,3,180,2000,D) and c<=5):
                        for bits in product((0,1),repeat=len(ports)):
                            nn=n+sum(z*r*b for (z,r),b in zip(ports,bits))
                            kk=k+sum(r*b for (_,r),b in zip(ports,bits))
                            assert (v*kk>=nn if up else v*kk<=nn)
                            checks+=1
                    # Both W and the conditional quotient bound the same
                    # physical continuation, so include the unconditioned fallback.
                    lookup[d,c,m,t]=(min if up else max)(v,W[d][c]);states+=1
        # d+1 is no longer needed once this entire slab was computed.
        if d<D:
            for c in range(1,min(K,d+1)+1):
                for m in moms:
                    for t in cells(c):
                        del lookup[d+1,c,m,t]
    # All nonfounder mothers partition into degree 2 and >=3 at the actual birth.
    roots=tuple(Q(read(2,1,m,(0,0)),G) for m in (2,-3))
    return roots,states,checks


def birth_mix(roots,up):
    # One degree-2 mother birth per nonfounder at most: the first exit from
    # its virgin residence, if by promotion rather than leaf gain.
    probabilities=(Q(1,2)/(base.HI+2),Q(1,2)/(base.LO+Q(3,2)))
    return (max if up else min)(p*roots[0]+(1-p)*roots[1] for p in probabilities)


def decimal(value,up):
    n=rnd(value.numerator*10**9,value.denominator,up)
    return f'{n//10**9}.{n%10**9:09d}'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--degree',type=int,default=360)
    p.add_argument('--arity',type=int,default=18)
    p.add_argument('--mother',type=int,default=2)
    p.add_argument('--tip',type=int,default=3)
    p.add_argument('--side',choices=('lower','upper'),action='append')
    p.add_argument('--bound',choices=('lower','upper'),action='append')
    p.add_argument('--preauthenticated',action='store_true')
    a=p.parse_args()
    assert 2<=a.degree<6000 and 1<=a.arity<=30 and a.mother>=2 and (a.tip==0 or a.tip>=2)
    for side in a.side or ('lower','upper'):
        with base.Snapshot(base.BASE/f'quintet_{side}.bin',side,verify=not a.preauthenticated) as s:
            for bound in a.bound or ('lower','upper'):
                up=bound=='upper'; H,C=(6000,64) if up else (12000,30)
                W=broad(s,H,C,a.degree+1,up)
                roots,n,checks=conditional(s,W,a.degree,a.arity,a.mother,a.tip,up)
                v=birth_mix(roots,up)
                mul=base.HI if up else base.LO
                overflow=65000*s.scale/(1+base.LO)**C if up else 0
                old=mul*Q(W[2][1],G)+overflow; new=mul*v+overflow
                improvement=(old-new) if up else (new-old)
                assert improvement>=0
                print(side,bound,'states',n,'corners',checks,'roots',roots,'mixed lifetime ticks',v,
                      'old',decimal(old,up),'uniform',decimal(mul*(max if up else min)(roots)+overflow,up),
                      'new',decimal(new,up),'improvement ticks',decimal(improvement,False),flush=True)
            s.check_unchanged()


if __name__=='__main__': main()
