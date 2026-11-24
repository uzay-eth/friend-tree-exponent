"""Rational continuation and local ledger at x=1/3 (no old-tree sweep)."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import factorial, prod
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path

S=10**10
u=Q(2,3)

def up(v):
    return -((-v*S)//1)

@cache
def root_bound(k):
    lo,hi=S,2*S
    while hi-lo>1:
        j=(hi+lo)//2
        if 2*j**k>=3*S**k: hi=j
        else: lo=j
    return Q(hi,S)

@cache
def beta(n,k,prom=False):
    """Upper beta survival, divided by u^(n-1), in grid units."""
    r=Q(1,k)
    # sin(pi r)/(pi r) <= 1-z^2/6+z^4/120, decreasing on 0<=z<=pi/2.
    z=Q(314159,100000)*r
    normal=r*(1-z*z/6+z**4/120)*factorial(n-1)
    if prom:
        normal/=r*prod(i-r for i in range(1,n))
        a=1+r
    else:
        normal/=prod(i-r for i in range(1,n+1))
        a=r
    q=n-a
    coeff=Q(1)
    value=Q(1)
    for j in range(1,25):
        coeff*=u*(j-a)/j
        value+=coeff*q/(q+j)
    if not prom:
        value+=coeff*q/(q+25)*u/(1-u)
    return up(normal*root_bound(k)*(1 if prom else u)*value)

def tips(b):
    m=b+1
    D=m+14
    row={(k,j):beta(D,k) for k in range(2,D) for j in range(k)}
    for n in range(D-1,m-1,-1):
        nxt={}
        for k in range(2,n):
            cap=beta(n,k)
            prom=beta(n+1,k,True)
            for j in range(k):
                gain=row[k+1,j+1]
                hold=row[k,j]
                loss=row[k,j-1] if j else 0
                # Retain all actual siblings, and bound all other incoming
                # sources at the parent by their reciprocal-degree cap 1/2.
                y=prom+j*loss
                z=k*(j+1)
                small=y+z*gain+(k*n-z-j-1)*hold
                zz=k*(k+j+1)
                large=2*y+zz*gain+(2*k*n-zz-2*(j+1))*hold
                val=max(2*small,large)
                nxt[k,j]=min(cap,(val+3*k*n-1)//(3*k*n))
        row=nxt
    return {k:Q(v,S) for k,v in row.items()}

def ledger(b,M,fee,ft,a=Q(0),base=Q(0)):
    U=tips(b)
    g=M-U[2,0]
    if g<=0: return False,'serial'
    def profile(d,c,r):
        l=d-c
        s=(l+r)/(d*(d+1)); t=Q(1,2*d)
        C=g*(d*s+l*t)-l*max(0,t-s)*U[2,0]
        if l: C-=l*max(0,s-t)*U[d+1,l-1]
        return s,C
    if all(C>=c*s*U[d+1,d-c]+(d*s+Q(d-c,2*d))/10000
           for d in range(2,b-1) for c in range(d+1) for r in (Q(0),Q(c,2))
           for s,C in (profile(d,c,r),)):
        return True,'sitewise reserve'
    end=[]
    for d in range(2,b-1):
        if profile(d,0,Q(0))[1]<=0: return False,('star',d)
        for c in range(1,d+1):
            if d==c==2: continue
            F=ft if d==2 else fee
            for k,r in ((k,r) for k in range(c+1) for r in (Q(k,2),Q(c,2))):
                s,C=profile(d,c,r)
                if C+a*(c-2)-base-k*F-(c-k)*max(0,s*U[d+1,d-c]-base)<0:
                    return False,('site',d,c,k,r)
            for r in (Q(1,2),Q(c,2)):
                s,_=profile(d,c,r)
                end.append((Q(1,d),s,U[d+1,d-c],F))
    # Deduplicate exact endpoint types; integer enclosure accelerates pair sweep.
    X=6*10**8
    Y=10**9
    def ceil(z): return -((-z)//1)
    end=list({(Q(r*X,6)//1,ceil(r*X/6),s*X//1,ceil(s*X),ceil(v*Y),f*X*Y//1)
              for r,s,v,f in end})
    v=ceil(Y*U[3,0]); credit=2*g*Y//1
    @cache
    def port(a,b,e):
        _,_,s,t,w,f=e
        return max(0,(b-s)*v,(t-a)*w)-f
    for length in (1,2,3):
        for x,y in product(end,repeat=2):
            if length==1:
                chain=((x[0]+y[0],x[1]+y[1]),)
            else:
                chain=((X//12+x[0],X//12+x[1]),)+((X//6,X//6),)*(length-2)+((X//12+y[0],X//12+y[1]),)
            debt=port(*chain[0],x)+port(*chain[-1],y)
            debt+=v*sum(max(t-a,b-s) for (a,b),(s,t) in zip(chain,chain[1:]))
            debt-=credit*sum(a for a,_ in chain)
            if debt>base*X*Y: return False,('chain',length,x,y,Q(debt,X*Y))
    return True,'chain reserve'

if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--start',type=int,default=21)
    parser.add_argument('--stop',type=int,default=40)
    parser.add_argument('--ledger-only',action='store_true')
    args=parser.parse_args()
    def atan_part(r,n):
        return sum((-1)**i*r**(2*i+1)/(2*i+1) for i in range(n))
    assert 16*atan_part(Q(1,5),8)-4*atan_part(Q(1,239),3)>Q(314159,100000)
    joint=run_path(str(Path(__file__).resolve().parents[1]/'first-sender/joint/check.py'))
    M=Q(211,1000)
    if not args.ledger_only:
        assert joint['coefficient'](21,768,56)>M
        print('PASS pinned serial coefficient > .211 for b>=21 at x=1/3',flush=True)
    U=tips(21)
    assert U[2,0]<Q('.158') and U[3,0]<Q('.070')
    assert U[3,1]<U[3,0] and U[4,2]<U[4,0]
    print('PASS physical tip caps at b=21: (2,0) < .158, (3,0) < .070; strict stock gains',flush=True)
    assert 21 <= args.start <= args.stop <= 40
    for b in range(args.start,args.stop):
        ok,where=ledger(b,M,*map(Q,('.004','.008','.004','.004')))
        assert ok,(b,where)
        print(f'PASS b={b}: physical continuation and {where}',flush=True)
