"""Outward limiting stock/degree comparison, then the frozen signed ledger."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from math import factorial
from pathlib import Path
from runpy import run_path
from types import FunctionType

HERE=Path(__file__).resolve().parent
FROZEN=HERE.parent/'annealed/extension/third-physical/scaled'
S=10**20
D=32
RATE=D+1


def ceil(p,q):
    return -((-p)//q)


def exp_floor(x):
    """Downward exp(x) for 0<=x<=32, argument on grid; also upper."""
    assert 0<=x<=32*S
    term, hi = S,S
    lo, upper = S,S
    for j in range(1,161):
        term=term*x//(j*S)
        hi=ceil(hi*x,j*S)
        lo+=term
        upper+=hi
    return lo,upper+ceil(hi*x,161*S-x)


def exp_neg(x):
    a,b=exp_floor(x)
    return S*S//b,ceil(S*S,a)


def root_floor(k,x):
    """Integer floor of kth root; floating seed is not trusted."""
    z=max(1,int(pow(x,1/k))) if x.bit_length()<1024 else 1<<(ceil(x.bit_length(),k))
    # Use a Newton step only if the initial guess is further wrong than a few units.
    if z**k>x:
        while True:
            t=((k-1)*z+x//z**(k-1))//k
            if t>=z:
                return z
            z=t
    if (z+1)**k>x:
        return z
    hi=2*z+2
    assert hi**k>x
    while hi-z>1:
        mid=(hi+z)//2
        if mid**k<=x:
            z=mid
        else:
            hi=mid
    return z


def moment(k,x):
    return root_floor(k, x*S**(k-1))


def low_integral(k,x,ex,z):
    # Integral_0^x t^(1/k-1) exp(-t)dt, Kummer positive expansion.
    term=total=k*S
    for j in range(1,161):
        term=term*x*k//(S*(j*k+1))
        if not term:
            break
        total+=term
    return ex*z*total//S**2


@cache
def gamma(k):
    x=32*S
    z=moment(k,x)
    ex,hi=exp_neg(x)
    lower=low_integral(k,x,ex,z)
    # Exact upper by using an upper grid recurrence for the positive expansion,
    # and t^(1/k-1)<=1 in the remaining tail.
    term=total=k*S
    for j in range(1,161):
        term=ceil(term*x*k,S*(j*k+1))
        total+=term
    total+=ceil(term*x*k, S*(161*k+1)-x*k)
    return lower,ceil(hi*(z+1)*total,S**2)+hi


def tails(k,x,ex,hi):
    z=moment(k,x)
    low,high=gamma(k)
    pending=min(S,S-low_integral(k,x,ex,z)*S//high)
    promo=min(S,pending+ceil(k*hi*(z+1),low))
    return pending,promo


def poisson(x,terms):
    _,p=exp_neg(x)
    w=[p]
    low,_=exp_neg(x)
    lw=[low]
    for j in range(1,terms+1):
        w.append(ceil(w[-1]*x,S*j))
        lw.append(lw[-1]*x//(S*j))
    return w,S-sum(lw)


def cap_box(lo,hi,N=256):
    """Upper limiting normalized caps on closed y box; T=2, degree stop D."""
    step=2*S//N
    assert step*N==2*S
    delta,_=exp_floor(step)
    ws,tail=poisson(RATE*step,12)
    ys=[lo*S//1]
    for i in range(N):
        ys.append(ys[-1]*delta//S)
    # For pointwise time and threshold monotonicity of each gamma survival,
    # take the left endpoint of the time slab and the lower y endpoint.
    row=None
    offsets={k:k*(k-1)//2-1 for k in range(2,D+1)}
    size=D*(D+1)//2-1
    for i in range(N,-1,-1):
        v=ys[i]
        ex,high=exp_neg(v)
        ends=[None,None]+[tails(k,v,ex,high) for k in range(2,D+1)]
        locks=[ends[k][0] for k in range(2,D+1) for j in range(k)]
        if i==N:
            row=locks
            continue
        accu=[w*ws[0]+tail*S for w in row]
        cur=row[:offsets[D]]+locks[offsets[D]:]
        for w in ws[1:]:
            nxt=locks.copy()
            for k in range(2,D):
                base=offsets[k]
                gg=offsets[k+1]+1
                p=ends[k][1]
                den=2*k*RATE
                for j in range(k):
                    h=cur[base+j]
                    dg=cur[gg+j]-h
                    gain=(k*(k+j+1) if dg>=0 else 2*k*(j+1))*dg
                    loss=2*j*(cur[base+j-1]-h) if j else 0
                    nxt[base+j]=ceil(h*den+gain+loss+2*(p-h),den)
            for s in range(size):
                accu[s]+=w*nxt[s]
            cur=nxt
        row=[min(l,ceil(v,S)) for l,v in zip(locks,accu)]
    # Relative to lo, each normalized payment e^y T_a(y*exp(t))
    # increases by at most sqrt(hi/lo), by the shifted gamma integral.
    _,sc=exp_floor(ceil(lo.numerator*S,lo.denominator))
    scale=ceil(sc*(1+root_floor(2,ceil(hi.numerator*lo.denominator*S*S,
                                      hi.denominator*lo.numerator))),S)
    return {(k,j):Q(ceil(10**7*scale*row[offsets[k]+j],S*S)+1,10**7)
            for k in range(2,D+1) for j in range(k)}


def pinned(name,blob):
    path=FROZEN/name
    data=path.read_bytes()
    assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()==blob
    return run_path(str(path))


def ledger(U,y,M,stock=False):
    f=pinned('check.py','85e9ce571c1b6a842e6d3aa8faa71dbda147b2df')
    assert y>=1 and y<=Q(5,2)
    assert 2*3**32>5*2**33 and Q(33,32)*Q(3,2)<2
    def cap(k,Y):
        if k>=D:
            return f['cap'](k,Y)
        return max(U[k,j] for j in range(k))
    if not stock:
        return FunctionType(f['chains'].__code__,dict(f,cap=cap))(y,M)
    # Preserve every signed term, tariff, corner test and high-degree bound;
    # just select the actual stock at each occurrence of an adverse cap.
    from inspect import getsource
    src=getsource(f['chains'])
    src=src.replace('V=cap(d+1,Y)\n        return', 'V=Ustock(d+1,l-1) if l else 0\n        return')
    src=src.replace('V=cap(d+1,Y)\n            for', 'V=Ustock(d+1,d-c)\n            for')
    src=src.replace('recv=ceil(V*cap(3,Y))','recv=ceil(V*Ustock(3,0))')
    assert 'V=cap(d+1,Y)' not in src and 'recv=ceil(V*cap(3,Y))' not in src
    def ust(k,j):
        return U[k,j] if k<D else f['cap'](k,y)
    ns=dict(f,cap=cap,Ustock=ust)
    exec(src,ns)
    return ns['chains'](y,M)


@cache
def serial_replay():
    # Parameterize only the two ends of the already frozen serial proof.
    from inspect import getsource
    f=pinned('near_two.py','56cb515405295df48198b405fc552a423e7a3bed')
    src=getsource(f['lower']).replace('def lower():','def lower(lo,hi):')
    src=src.replace('Q(8)','4*hi').replace('Q(8**i,factorial(i))','(4*hi)**i/factorial(i)').replace('Q(199,100)','lo')
    assert 'Q(8)' not in src and 'Q(8**i' not in src and 'Q(199,100)' not in src
    exec(src,f)
    return f['lower']


def certify(lo,hi,M,N,sub=Q(1,100)):
    U=cap_box(lo,hi,N)
    result=ledger(U,lo,M,stock=True)
    assert result==(True,945),(lo,hi,result)
    ll=lo
    minserial=Q(1)
    while ll<hi:
        hh=min(hi,ll+sub)
        minserial=min(minserial,serial_replay()(ll,hh))
        ll=hh
    assert minserial>M,(lo,hi,float(minserial),M)
    print('PASS',lo,hi,'M',M,'serial >',Q(minserial*10**7//1,10**7),
          'U20',U[2,0],'U30',U[3,0],'U31',U[3,1],result,flush=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('--lo',type=Q,default=Q('1.99'))
    p.add_argument('--hi',type=Q,default=Q('2'))
    p.add_argument('--steps',type=int,default=256)
    p.add_argument('--caps-only',action='store_true')
    p.add_argument('--ignore-stock',action='store_true')
    p.add_argument('--suite',action='store_true')
    p.add_argument('--M',type=Q,default=Q('.462'))
    a=p.parse_args()
    if a.suite:
        cuts=list(map(Q,('1.90','1.94','1.99','2.05','2.16','2.28','2.40','2.50')))
        serial=list(map(Q,('.462','.462','.4618','.4601','.4585','.4570','.4558')))
        for lo,hi,M in zip(cuts,cuts[1:],serial):
            certify(lo,hi,M,a.steps)
        raise SystemExit(0)
    U=cap_box(a.lo,a.hi,a.steps)
    print('CAPS',a.lo,a.hi,a.steps,[(k,j,float(U[k,j])) for k,j in [(2,0),(3,0),(3,1),(4,0),(4,2)]],flush=True)
    if not a.caps_only:
        result=ledger(U,a.lo,a.M,stock=not a.ignore_stock)
        print('LEDGER',result,flush=True)
        assert result[0],result
