"""Exact threshold-box replay with a fixed-label serial lower."""
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from pathlib import Path
from runpy import run_path
from types import FunctionType

HERE=Path(__file__).resolve().parent
parent=run_path(str(HERE.parent/'check.py'))
S=parent['S']
up=parent['up']
factorial,prod=parent['factorial'],parent['prod']

def cdiv(a,b):
    return -((-a)//b)

@cache
def power(q,k):
    """Ceiling of q^(1-1/k) on the integer grid."""
    a,b=q.numerator,q.denominator
    L,H=0,S
    while H-L>1:
        y=(L+H)//2
        if y**k*b**(k-1)>=S**k*a**(k-1): H=y
        else: L=y
    return Q(H,S)

@cache
def norm(n,k,prom):
    r=Q(1,k)
    z=Q(314159265,10**8)*r
    # Alternating upper, with the negative term of order six included.
    sn=sum((-1)**j*z**(2*j)/factorial(2*j+1) for j in range(5))
    return r*sn*factorial(n-1)/(r*prod(i-r for i in range(1,n)) if prom
                                else prod(i-r for i in range(1,n+1)))

@cache
def beta(n,k,prom,q):
    r=Q(1,k); a=r+int(prom); B=n-a
    p=v=Q(1)
    for j in range(1,25):
        p*=q*(j-a)/j
        v+=p*B/(B+j)
    if not prom:
        v+=p*B/(B+25)*q/(1-q)
    return norm(n,k,prom)*v

@cache
def tips(b,x,y):
    q=1-x; low=1-y
    m=b+1; D=m+14
    def pending(n,k): return up(power(q,k)*beta(n,k,False,q))
    row={(k,j):pending(D,k) for k in range(2,D) for j in range(k)}
    for n in range(D-1,m-1,-1):
        nxt={}
        for k in range(2,n):
            cap=pending(n,k)
            # q^(1-r)(1-qv)^r increases if q<=1-r. Otherwise bound separately.
            prom=up(power(q,k)*beta(n+1,k,True,q if q<=1-Q(1,k) else low))
            for j in range(k):
                gain,hold=row[k+1,j+1],row[k,j]
                loss=row[k,j-1] if j else 0
                z=k*(j+1); zz=k*(k+j+1)
                a=2*(z*gain+(k*n-z-j-1)*hold)
                h=zz*gain+(2*k*n-zz-2*(j+1))*hold
                hn,hd=q.numerator,q.denominator
                val=2*prom*hd+hn*(2*j*loss+max(a,h))
                nxt[k,j]=min(cap,cdiv(val,2*k*n*hd))
        row=nxt
    return {k:Q(v,S) for k,v in row.items()}

def ledger(b,M,x,y):
    f=FunctionType(parent['ledger'].__code__,dict(parent,tips=lambda b:tips(b,x,y)))
    return f(b,M,*map(Q,('.004','.008','.004','.004')))

def serial(b,n=320,d=24):
    """Exact marginal CDF and exact shifted-endpoint survival at fixed b."""
    tmax=8*d
    surv=[(1+Q(i,d*b))**(-b) for i in range(tmax+1)]
    weights=[S*(a-z)//1 for a,z in zip(surv,surv[1:])]
    fam=[(1+Q(i,d*b)) for i in range(1,tmax)]
    total=0
    # Group identical f values from the right endpoints of the uniform bins.
    from collections import Counter
    bins=Counter((b-1)*(n-j)//j for j in range(1,n+1))
    for f,count in bins.items():
        C=Q(b,b-1)*sum(Q(1,k) for k in range(f+2,f+b+1))
        if C>=8: continue
        prev=0; dist=[]
        for k,v in enumerate(fam,1):
            mass=(S*(1-v**(-f-1)))//1
            dist.append((mass-prev, up(C*v)))
            prev=mass
        for i in range(1,tmax):
            value=0
            for k,(mass,charge) in enumerate(dist[:i-1],1):
                p=S-cdiv(charge*d,i-k)
                if p<=0: break
                value+=mass*p
            total+=count*weights[i]*value
    return Q(total,n*S**3)

if __name__=='__main__':
    from argparse import ArgumentParser
    p=ArgumentParser(description=__doc__)
    p.add_argument('--start',type=int,default=20)
    p.add_argument('--stop',type=int,default=40)
    args=p.parse_args()
    assert 20<=args.start<=args.stop<=40
    def pin(path,expected):
        data=path.read_bytes()
        assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()==expected
    pin(HERE.parent/'check.py','e39f8bfdf66f19e34a59247fa3649ccaf884358b')
    pin(HERE.parents[1]/'first-sender/joint/check.py','d3b4ad5bc87760a843712d9af43af76450a07a37')
    def atan(v,n):
        return sum((-1)**j*v**(2*j+1)/(2*j+1) for j in range(n))
    assert 16*atan(Q(1,5),8)-4*atan(Q(1,239),3)>Q(314159265,10**8)
    assert serial(20)>Q('.2124')
    print('PASS b=20 serial coefficient > .2124',flush=True)
    for b in range(args.start,args.stop):
        M=Q('.2124') if b==20 else Q('.211')
        cuts=(Q(1,3),Q(41,120),Q(7,20),Q(3,8),Q(2,5),Q(1,2)) if b<=21 else (Q(1,3),Q(1,2))
        for x,y in zip(cuts,cuts[1:]):
            ok,kind=ledger(b,M,x,y)
            assert ok,(b,x,y,kind)
        print(f'PASS b={b} over [1/3,1/2]',flush=True)
