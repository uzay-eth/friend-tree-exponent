"""Independent harmonic, directed stock, and exact chain-seam checks."""

from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
from runpy import run_path

S = 10**13
E = sum(Q((-1)**j, factorial(j)) for j in range(18))


def stock(p, degree):
    # Strip e^-1, start at (0,1), and generate from the harmonic equation.
    a = [Q(0), Q(1)]
    for j in range(1, degree):
        a.append(((j+1+p)*a[-1]-a[-2])/j)
    return a


def lower(lo, hi, D=40, K=96):
    after = [E*x*S//1 for x in stock(lo,D)]
    for d in range(D-1,1,-1):
        R = 0
        for k in range(K-1,1,-1):
            gain = k-1+Q(1,d)
            R = ((Q(S,k)+gain*R)/(hi+Q(1,k)+gain))//1
        row = [R*after[1]//S]
        for j in range(1,d+1):
            row.append(((j*after[j+1]+Q(j,d)*row[-1])/(hi+j+Q(j,d)))//1)
        after = row
    return Q(after[1],S)


def harmonic():
    rows = 0
    for p in (Q(0),Q(1,1024),Q(1,7),Q(1,2),Q(473,512),Q(1)):
        a = stock(p,81)
        rising = [Q(0),Q(1)]
        for j in range(1,81):
            rising.append(rising[-1]*(1+p/j))
        for j in range(1,81):
            assert a[j]==sum(rising[j-k]/factorial(k) for k in range(j))
            assert a[j+1]>a[j]>a[j-1]
            rows += 1
    assert lower(Q(0),Q(0))>Q(78,100)
    print(f'PASS {rows} exact harmonic/convolution rows; p=0 polynomial extension > .78',flush=True)


def chain_seam(src):
    lo,hi = Q(473,512),Q(474,512)
    M = lower(lo,hi)
    U,L = src['side_bounds'](lo,hi,M)
    g = M-U(2)
    B,F,FT = Q('.008'),Q('.004'),Q('.0131')
    ends = []
    for d in range(2,10):
        for c in range(1,d+1):
            if d==c==2:
                continue
            for r in (Q(1,2),Q(c,2)):
                s = (d-c+r)/(d*(d+1))
                ends.append((Q(1,d),L(d+1),s,U(d+1),FT if d==2 else F))
    for d in range(10,31):
        ends.append((Q(1,d),L(d+1),None,U(d+1)/(d+1),F))
    for r in (Q(0),Q(1,31)):
        ends.append((r,0,None,1/(lo*32**2),F))

    def port(s,e):
        r,Lp,t,U_p,fee = e
        f = max(0,U(3)-Lp)
        if t is None:
            return max(s*f,U_p)-fee
        return max(0,(s-t)*f,(t-s)*max(0,U_p-L(3)))-fee

    smallest = None
    for n,x,y in product((1,2,3),ends,ends):
        r,t = x[0],y[0]
        line = ((r+t)/6,) if n==1 else ((Q(1,2)+r)/6,)+(Q(1,6),)*(n-2)+((Q(1,2)+t)/6,)
        cost = port(line[0],x)+port(line[-1],y)
        cost += max(0,U(3)-L(3))*sum(abs(a-b) for a,b in zip(line,line[1:]))
        residual = B-cost+2*g*sum(line)
        smallest = residual if smallest is None else min(smallest,residual)
    assert smallest>0
    print(f'PASS {3*len(ends)**2} exact chain pairs on the threshold box; slack > {smallest*S//1}/{S}',flush=True)


if __name__ == '__main__':
    harmonic()
    writer = Path(__file__).resolve().parents[1]/'open/age-order/hub-ranking/fractional-middle/check.py'
    src = run_path(str(writer))
    minimum = Q(1)
    for j in range(256):
        lo,hi = Q(j,256),Q(j+1,256)
        m = lower(lo,hi)
        assert m==src['middle'](lo,hi)>Q(11,20)
        minimum = min(minimum,m)
    print(f'PASS independent recomputation of 256 lower boxes; minimum {minimum}',flush=True)
    chain_seam(src)
