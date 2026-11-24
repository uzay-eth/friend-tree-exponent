"""Exact physical adverse-tip bounds at x=1/3; no mixed-moment inference."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from math import factorial, prod
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
EXT = HERE.parents[2]/'global-transport/annealed/extension'
S = 10**10
U = Q(2,3)


def grid(x):
    return -(-x*S//1)


@cache
def power(a):
    # Upper rational on (2/3)^a, including negative a.
    target = U**a.numerator
    q = a.denominator
    lo, hi = 0, 4*S
    assert Q(hi,S)**q >= target
    while hi-lo>1:
        m=(lo+hi)//2
        if Q(m,S)**q >= target:
            hi=m
        else:
            lo=m
    return Q(hi,S)


@cache
def beta(n,a):
    """Upper in ticks on Beta(a,n-a) survival divided by U^(n-1)."""
    i=a//1
    r=a-i
    if a==1:
        return S
    assert 0<r<1 and 0<=i<=2 and n>a
    t=Q(314159,100000)*min(r,1-r)
    sine=min(r,1-r)*(1-t*t/6+t**4/120)
    factor=factorial(n-1)*sine / (prod(r+j for j in range(i))*prod(j-r for j in range(1,n-i+1)))
    q=n-a
    coeff=value=Q(1)
    for j in range(1,25):
        coeff*=U*(j-a)/j
        value+=coeff*q/(q+j)
    value+=abs(coeff)*q/(q+25)*U/(1-U)
    assert value>0
    return grid(factor*power(1-a)*value)


@cache
def after(n,k):
    """Tag has just promoted: fresh leaf, permanent parent of degree k."""
    r=Q(1,k)
    a=1+r
    # Exact next-step probabilities; all other arrows leave the cap at a.
    value=(a*beta(n+1,a+1)+(n-a-1)*beta(n+1,a)
           +Q(beta(n+1,a-Q(1,2))+beta(n+1,1+Q(1,k+1)),2))
    return min(beta(n,a), -(-(U*value/n)//1))


def tips(b, refined=True, post=True):
    first=b+1
    end=first+14
    row={(k,j):beta(end,Q(1,k)) for k in range(2,end) for j in range(k)}
    for n in range(end-1,first-1,-1):
        nxt={}
        for k in range(2,n):
            cap=beta(n,Q(1,k))
            prom=after(n+1,k) if post else beta(n+1,1+Q(1,k))
            for j in range(k):
                gain,hold=row[k+1,j+1],row[k,j]
                loss=row[k,j-1] if j else 0
                rlo=j+1
                rhi=Q(k+j+1,2) if refined else n-Q(j+1,k)
                def v(r):
                    return U*(r*gain+(n-r-Q(j+1,k))*hold+Q(prom+j*loss,k))/n
                nxt[k,j]=min(cap,-(-max(v(rlo),v(rhi))//1))
        row=nxt
    return {k:Q(v,S) for k,v in row.items()}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--serial',action='store_true',help='also rerun the pinned serial quadrature')
    args=parser.parse_args()
    p=lambda r,n: sum((-1)**i*r**(2*i+1)/(2*i+1) for i in range(n))
    assert 16*p(Q(1,5),8)-4*p(Q(1,239),3)>Q(314159,100000)
    for name,pin in (
        ('third-physical/README.md','27c6b6c7a4aa559e166835f3adc17fd189c989bb'),
        ('third-physical/check.py','e39f8bfdf66f19e34a59247fa3649ccaf884358b'),
        ('first-sender/joint/check.py','d3b4ad5bc87760a843712d9af43af76450a07a37'),
        ('first-sender/check.py','5680c8e937ebafdec09c2e193e6b36fe47ac84b1'),
    ):
        raw=(EXT/name).read_bytes()
        assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==pin
    original=run_path(str(EXT/'third-physical/check.py'))
    if args.serial:
        joint=run_path(str(EXT/'first-sender/joint/check.py'))
        serial=joint['coefficient'](21,768,56)
        assert serial==Q(13530815625626578442946382009,64000000000000000000000000000)>Q(211,1000)
        print('PASS serial lower coefficient:',serial,flush=True)
    old=original['tips'](21)
    new=tips(21)
    for k,j in ((2,0),(3,0),(3,1),(4,0),(4,2),(5,3)):
        assert new[k,j] < old[k,j]
        print((k,j), 'prior',old[k,j], 'refined',new[k,j], 'absolute upper',new[k,j]*U**21,flush=True)
    # This substitution tests the same local ledger, with the new upper.
    original['ledger'].__globals__['tips']=lambda b: new
    trial=original['ledger'](21,Q(211,1000),*map(Q,('.004','.008','.004','.004')))
    assert trial[0],trial
    print('PASS b=21 refined physical tip upper and signed ledger:',trial[1],flush=True)
    # Stars bypass skeleton summation: directly check the same reserve.
    g=Q(211,1000)-new[2,0]
    for d in range(2,20):
        s,t=Q(1,d+1),Q(1,2*d)
        credit=g*(d*s+d*t)-d*max(0,s-t)*new[d+1,d-1]
        assert credit>=Q(3,250)
    assert trial[1]=='chain reserve'
    print('PASS 3/250 normalized numerator reserve including stars; gap >=',Q(3,250)*U**21/(20*21),flush=True)


if __name__=='__main__':
    main()
