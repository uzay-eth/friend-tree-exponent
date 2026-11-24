"""Verify the fringe transform, newborn derivatives, and star tail with fractions."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import ceil, factorial, floor

A = F(1,2)
BETA = F(1,16)
S = F(113123,200000)
UNIT = 10**12

def newborn(tag=2):
    return (tag,1,())

def plus(tree):
    t,l,cc=tree
    return t,l+1,cc

@cache
def info(tree):
    _,l,cc=tree
    return 1+l+len(cc),1+l+sum(info(y)[1] for y in cc)

def mark(tree):
    return F(1,tree[0])

def b(tree):
    return F(tree[1],info(tree)[0])

@cache
def physical(tree,parent=F(0)):
    t,l,cc=tree
    degree=info(tree)[0]
    out=defaultdict(F)
    gain=l+parent+sum(F(1,info(y)[0]) for y in cc)
    if gain:
        out[plus(tree)]+=gain
    if l:
        ch=newborn(2 if degree<16 else 16)
        out[t,l-1,tuple(sorted((*cc,ch)))]+=F(l,degree)
    for i,y in enumerate(cc):
        for yy,v in physical(y,F(1,degree)):
            nxt=t,l,tuple(sorted(cc[:i]+(yy,)+cc[i+1:]))
            out[nxt]+=v
    return tuple(out.items())

def L(tree,f):
    return sum(rate*(f(nxt)-f(tree)) for nxt,rate in physical(tree))

def D(tree,f):
    return sum(f(ch) for ch in tree[2])

def T(tree,f,z=F(1)):
    return f(tree)+z*sum(T(ch,f,z) for ch in tree[2])

def R(tree,f,constant=True):
    aa=A if constant else mark(tree)
    p=F(1,info(tree)[0])
    value=L(tree,f)-aa*f(tree)
    for y in tree[2]:
        ay=A if constant else mark(y)
        value+=(ay-p)*f(y)+p*f(plus(y))
    return value

@cache
def derivative(tree,k,kind):
    if not k:
        return b(tree)
    f=lambda t:derivative(t,k-1,kind)
    return L(tree,f) if kind=='intrinsic' else R(tree,f)

def check_transform():
    def f(t):
        deg,n=info(t)
        return F(n*(t[1]+1),deg+2)+mark(t)+F(len(t[2]),3)
    starting={newborn(),newborn(16),(2,14,(newborn(16),)),(16,0,(newborn(),))}
    trees=set(starting)
    for n in range(3):
        trees|={y for t in trees for y,_ in physical(t)}
    count=0
    for t in trees:
        deg,n=info(t)
        assert sum(v for _,v in physical(t))==n-F(1,deg)
        for z in (F(1),F(4,5)):
            for const in (False,True):
                aa=lambda t: A if const else mark(t)
                def P(t):
                    p=F(1,info(t)[0])
                    nj=newborn(2 if info(t)[0]<16 else 16)
                    cp=sum((aa(y)-p)*f(y)+p*f(plus(y)) for y in t[2])
                    return L(t,f)-aa(t)*f(t)+z*(cp+b(t)*f(nj))
                # In the constant gauge erased birth marks can still be observed
                # as auxiliary metadata, so the direct newborn uses its real tag.
                assert T(t,P,z)==L(t,lambda y:T(y,f,z))-aa(t)*f(t)
                count+=1
    print('weighted full-generator intertwining identities:',count)

def tail_bounds(lam,K,N):
    assert 2<=K<=N
    mass=F(1)/(3+2*lam)
    partial=F(0)
    for j in range(2,N):
        if j>=K:
            partial+=mass
        mass*=F(j*j)/(j*j+(2+lam)*j+lam)
    gamma=lam+2/(lam+2)
    numer=(N+gamma)*mass
    return (partial+numer/(1+lam),
            partial+numer/(1+lam-lam*lam/((lam+2)*N*N)))

def show(interval):
    lo,hi=interval
    return f'[{floor(UNIT*lo)}/{UNIT},{ceil(UNIT*hi)}/{UNIT}]'

if __name__=='__main__':
    check_transform()
    for kind in ('intrinsic','constant'):
        vals=[derivative(newborn(),k,kind) for k in range(6)]
        print(kind,'derivatives 0 through 5',','.join(map(str,vals)))
        for p in (3,5,7):
            k=p-2
            unique=F(factorial(k)*(p-1),p)
            assert (vals[k]-unique).denominator%p and not vals[k].denominator%p
    for a in (A,BETA):
        lam=S+a
        for K in (2,16,64,256):
            interval=tail_bounds(lam,K,384)
            print('mark',a,'first-promotion mass from',K,show(interval))
        # Error after the constant high-degree continuation in the one-child state.
        qp=F(49,125)
        d=F(4,25)
        hmin=F(471910343,10**9)
        kappa=F(891876000,721498129)
        alpha=1/hmin
        wj=1/(BETA+d)
        force=1+alpha*BETA*wj/(lam-qp)
        error=(force*(2/lam+1/lam**2)
               +alpha*kappa*(1+BETA/(lam-qp))*wj*(1+1/lam)/(lam-qp))
        for K in (64,256,1024):
            er=error*tail_bounds(lam,K,max(K,384))[1]/K
            print('mark',a,'root tail continuation at',K,'absolute error <',ceil(UNIT*er),'/',UNIT)
    theta=F(51,50)*(1-F('0.12884'))
    assert F('4.062')**90*theta**1079<1
    for n in (100,180):
        error=F(51,50)*theta**n/(1-theta)
        print('common-row steps',n,'absolute error <',ceil(UNIT*error),'/',UNIT)
