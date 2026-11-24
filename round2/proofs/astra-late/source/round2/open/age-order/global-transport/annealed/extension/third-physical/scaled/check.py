"""Uniform birth-scale locked/physical cap and local debit checks."""
from fractions import Fraction as Q
from functools import cache
from math import factorial
from itertools import product
from pathlib import Path
from runpy import run_path

S=10**9

def root(k,x):
    """Lower fixed-grid kth root of x, for x>=1."""
    lo,hi=S,(x+1)*S
    while hi-lo>1:
        j=(lo+hi)//2
        if j**k<=S**k*x: lo=j
        else: hi=j
    return Q(lo,S)

@cache
def gamma_lower(k):
    # Gamma(1/k), split at 1. Right rectangles on [1,8].
    low=sum((-1)**j*Q(k,j*k+1)/factorial(j) for j in range(12))
    d=128
    step=Q(1,d)
    e=S
    mult=S*(1-step+step**2/2-step**3/6)//1
    for i in range(1,8*d+1):
        e=e*mult//S
        if i>d:
            # floor root via integer root of i*d^(k-1), divided by d
            low+=Q(e,S*i)*root(k,i*d**(k-1))/d
    return low

@cache
def cap(k,Y):
    if k>=32:
        return 2*Q(1,(k-1)*Y)
    # y^(r-1) = 1 / (y^(k-1))^(1/k).
    v=Q(1)/(root(k,Y**(k-1))*gamma_lower(k))
    return Q(-((-v*S)//1),S)

def chains(Y,M):
    a,base,fee,ft=map(Q,('.006','.004','.004','.008'))
    U=cap(2,Y);g=M-U; D=32
    assert g>0
    def profile(d,c,r):
        l=d-c; s=(l+r)/(d*(d+1)); t=Q(1,2*d); V=cap(d+1,Y)
        return s,g*(d*s+l*t)-l*(max(0,t-s)*U+max(0,s-t)*V)
    z=Q(2,D*Y)
    # Large cores pay all their own adverse directions, not just excess.
    # With l>d/2 the adverse leaf direction is nearly absent.
    rich=Q(D,D+1)*(g-z)/2+g/4-U/(2*(D+1))
    sparse=min(Q(0),(Q(D,D+1)*(g-z)+(g-U)/2)/2)
    assert min(rich,(a-fee)*D/2+sparse)>2*a+base
    assert z/(D+1)<base
    assert Q(D,D+1)*(g-z)>0
    end=set()
    for d in range(2,D):
        assert profile(d,0,Q(0))[1]>0
        for c in range(1,d+1):
            if d==c==2: continue
            f=ft if d==2 else fee
            V=cap(d+1,Y)
            for k in range(c+1):
                for r in (Q(k,2),Q(c,2)):
                    s,C=profile(d,c,r)
                    if C+a*(c-2)-base-k*f-(c-k)*max(0,s*V-base)<0:
                        return False,('site',d,c,k,r)
            for r in (Q(1,2),Q(c,2)):
                end.add((Q(1,d),profile(d,c,r)[0],V,f))
    for r in (Q(0),Q(1,D)):
        end.add((r,None,z/(D+1),fee))
    X=6*10**8; V=10**9
    def ceil(t): return -((-t)//1)
    recv=ceil(V*cap(3,Y)); credit=2*g*V//1
    E=[]
    for r,s,u,f in end:
        E.append((X*r//6,ceil(X*r/6),None if s is None else (X*s//1,ceil(X*s)),
                  ceil(V*u),X*V*f//1))
    @cache
    def port(lo,hi,e):
        _,_,s,u,f=e
        return (max(hi*recv,X*u) if s is None else max(0,(hi-s[0])*recv,(s[1]-lo)*u))-f
    for x,y in product(E,repeat=2):
        for length in (1,2,3):
            if length==1:
                ch=((x[0]+y[0],x[1]+y[1]),)
            else:
                ch=((X//12+x[0],X//12+x[1]),)+((X//6,X//6),)*(length-2)+((X//12+y[0],X//12+y[1]),)
            cost=port(*ch[0],x)+port(*ch[-1],y)-credit*sum(t for t,_ in ch)
            cost+=recv*sum(max(v-lo,hi-u) for (lo,hi),(u,v) in zip(ch,ch[1:]))
            if cost>base*X*V: return False,('chain',length,x,y,float(Q(cost,X*V)))
    return True,len(end)

if __name__=='__main__':
    from argparse import ArgumentParser
    p=ArgumentParser(description=__doc__)
    p.add_argument('--finite',action='store_true')
    p.add_argument('--start',type=int,default=20)
    p.add_argument('--stop',type=int,default=39)
    args=p.parse_args()
    if args.finite:
        B=run_path(str(Path(__file__).resolve().parents[1]/'band/check.py'))
        assert B['serial'](20)>Q('.2124')
        for b in range(args.start,args.stop):
            m=Q('.2124') if b==20 else Q('.211')
            l=Q(1,3)-Q(1,2400)
            ok,kind=B['ledger'](b,m,l,Q(1,3))
            assert ok,(b,kind)
        print(f'PASS finite older labels [{args.start},{args.stop}) on [799/2400,1/3]')
    else:
        assert 2*3**32>13*2**33 and Q(33,32)*Q(3,2)<2
        ok,info=chains(13,Q('.211'))
        assert ok,info
        assert info==735
        print('PASS unbounded degree/label ledger at y=(b+1)x>=13; 735 endpoint types')
