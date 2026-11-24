"""Birth-scale sharp-tail enclosures using continuous physical clocks."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from math import factorial
from pathlib import Path
from runpy import run_path
from types import FunctionType, SimpleNamespace

HERE = Path(__file__).resolve().parent
S = 10**9


def pinned(path, pin):
    data = path.read_bytes()
    assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()==pin, path
    return run_path(str(path))


FR = pinned(HERE.parents[2]/'hub-ranking/mean-ledger/fractional/check.py',
            '76bbdfc337be47cf8a0a91138ee5d191264017ab')
HUB = pinned(HERE.parents[1]/'hostile/check.py','b5b63d953c0ae51a6d145fa3c00452c4493c1b95')


def up(q):
    return -((-S*q)//1)


@cache
def exp_bounds(x):
    """Rational grid enclosing exp(-x), by scaling then squaring."""
    assert x >= 0
    if x >= 30:
        return 0,1
    k = 0
    while x > Q(1,4):
        x /= 2
        k += 1
    lo = sum((-x)**i/factorial(i) for i in range(12))
    hi = lo+x**12/factorial(12)
    T = S**2
    a,b = lo*T//1,-((-hi*T)//1)
    for _ in range(k):
        a,b = a*a//T,(b*b+T-1)//T
    return a//S,(b+S-1)//S


@cache
def root_upper(x,k):
    a,b = 0, S*max(1,-((-x)//1))
    while b-a>1:
        i = (a+b)//2
        if Q(i,S)**k >= x:
            b = i
        else:
            a = i
    return Q(b,S)


@cache
def gamma_normal_lower(k):
    # Integral of t^(1/k) exp(-t) by a concave chord on each cell.
    # Integrate the chord exactly against exp(-t); use bounded atoms.
    q = Q(1,32)
    n = 320
    # Integrate the chord times the minimum density on each cell.
    L = [max(0,root_upper(j*q,k)-Q(1,S)) for j in range(n+1)]
    # The density's linear lower e^(-left)*(1-t) improves the quadrature.
    return sum((L[j]*(q/2-q*q/6)+L[j+1]*(q/2-q*q/3))
               *exp_bounds(j*q)[0]/S for j in range(n))


def solve(z, upper, dcut=24, rcut=32, steps=128, scale=Q(33,32)):
    """Reachable CT returns; log threshold grid grows by `scale` per bin."""
    h = scale-1
    # Every elapsed log-time h multiplies the threshold by e^h. Use its
    # upper envelope in lower payoffs and its lower envelope in upper ones.
    stretch = 1+h+h*h/(2*(1-h)) if not upper else 1+h+h*h/2
    tick = z*S//1 if upper else -((-z*S)//1)
    zz = [Q(tick,S)]
    for _ in range(steps):
        tick = tick*stretch//1 if upper else -((-tick*stretch)//1)
        zz.append(Q(tick,S))
    assert zz[0]>0

    def decreasing_upper(row):
        for i in range(steps-1,-1,-1):
            row[i] = max(row[i],row[i+1])
        return row

    @cache
    def kernel(rate):
        # Integral event masses by bin with directed rounding. Accumulate
        # survival bounds without repeatedly calling the exponential.
        a,b = exp_bounds(rate*h)
        l,u = S,S
        out = []
        for _ in range(steps):
            out.append(((l*(S-b))//S,(u*(S-a)+S-1)//S))
            l,u = l*a//S,(u*b+S-1)//S
        return out

    def race(branches, floor=False):
        rate = sum(r for r,_ in branches)
        weights = kernel(rate)
        row = [0]*(steps+1)
        if not floor:
            row[-1] = min(S,up(Q(3,2)/zz[-1]))
        for i in range(steps-1,-1,-1):
            payoff = 0
            # On the upper side extend the terminal value as a constant
            # for all later times; on the lower side discard those paths.
            for r,f in branches:
                if not r:
                    continue
                s = 0
                if floor:
                    for t in range(steps-i):
                        s += weights[t][0]*f[i+t+1]
                else:
                    # Upper rectangles on all finite cells and the tail.
                    mass = 0
                    for t in range(steps-i):
                        low,high = weights[t]
                        mass += low
                        s += high*f[i+t]
                    s += max(0,S-mass)*f[-1]
                payoff += r*s
            v = Q(payoff)/(S*rate)
            row[i] = max(0,min(S,v//1 if floor else -((-v)//1)))
        if not floor:
            decreasing_upper(row)
        else:
            for i in range(1,steps+1):
                row[i] = min(row[i],row[i-1])
        return row

    def controlled(fixed, moving, lo, hi, floor, boundary, extra=None):
        # On each time cell, bound a continuation payoff by its appropriate
        # monotone endpoint. Integrate the scalar differential comparison
        # with constant payoff. The endpoint rates give a valid barrier:
        # use the worst linear drift on the interval enclosed by the two
        # candidate trajectories (the constant-payoff one-state CTMDP).
        row = [0]*(steps+1)
        row[-1] = boundary[-1]
        for i in range(steps-1,-1,-1):
            j = i+floor
            F = sum(r*f[j] for r,f in fixed)
            R = sum(r for r,_ in fixed)
            levels = [(R+r,F+r*moving[j]) for r in (lo,hi)]
            if extra is not None:
                rate,f = extra
                levels += [(t+rate,p+rate*f[j]) for t,p in levels]
            w = row[i+1]
            # Explicit and implicit Euler bracket a linear segment whose
            # drift is outward for every admissible rate/branch choice.
            if floor:
                vals = [min(w+h*(p-r*w),(w+h*p)/(1+h*r)) for r,p in levels]
                z = max(0,min(vals)//1)
                row[i] = max(boundary[i],z)
            else:
                vals = [max(w+h*(p-r*w),(w+h*p)/(1+h*r)) for r,p in levels]
                z = min(S,-((-max(vals))//1))
                row[i] = min(boundary[i],max(w,z))
        if floor:
            for i in range(1,steps+1):
                row[i] = min(row[i],row[i-1])
        return row if floor else decreasing_upper(row)

    @cache
    def gamma(k, kind):
        norm = gamma_normal_lower(k)
        r = Q(1,k)
        # For shape 1+r, Jensen on the shifted Exp(1) survival integral.
        if kind == 'plain':
            return decreasing_upper([min(S,up(Q(exp_bounds(x)[1],S)*root_upper(x+1,k)/norm)) for x in zz])
        if kind == 'gain':
            # Recurrence Q(2+r,x)=Q(1+r,x)+e^-x*x^(1+r)/Gamma(2+r).
            return decreasing_upper([min(S,v+up(Q(exp_bounds(x)[1],S)*x*root_upper(x,k)/(norm*(1+r))))
                                     for x,v in zip(zz,gamma(k,'plain'))])
        # After promoting the fresh neighbour, its shape is <=1;
        # the exponential is an immediate uniform upper.
        assert kind == 'loss'
        return [exp_bounds(x)[1] for x in zz]

    @cache
    def promoted(k):
        a = 1+Q(1,k)
        v = race(((a,gamma(k,'gain')),(Q(1,2),gamma(k,'loss')),
                  (Q(1,2),gamma(k+1,'plain'))))
        return [min(x,y) for x,y in zip(v,gamma(k,'plain'))]

    @cache
    def U(k,j):
        hub = up(HUB['upper'](k,j)) if k < rcut else up(Q(1,j+1))
        norm = gamma_normal_lower(k)
        cap = [min(hub,up((1+Q(1,k))/x),
                   up(Q(exp_bounds(x)[1],S)*root_upper(x,k)/(k*norm*x))) for x in zz]
        if k>=rcut:
            return decreasing_upper(cap)
        prom = promoted(k)
        rise = U(k+1,j+1)
        fall = U(k,j-1) if j else rise
        return controlled(((Q(1,k),prom),(Q(j,k),fall)),rise,Q(j+1),Q(k+j+1,2),False,cap)

    @cache
    def poisson(j):
        row = []
        for x in zz:
            if x>=30:
                row.append(0)
                continue
            p,c = Q(1),Q(1)
            acc = Q(0)
            terms = [c]
            for i in range(1,j):
                p *= x/i
                c += p
                terms.append(c)
            for k in range(j):
                acc += terms[j-k-1]/factorial(k)
            row.append(Q(367879441,S)*exp_bounds(x)[0]*acc//1)
        return row

    @cache
    def ret(d,k):
        if k>=rcut:
            return [0]*(steps+1)
        pay = V(d+1,1)
        if upper:  # U and lower never mix in a physical row.
            raise AssertionError
        # Additional gains of the mark have rate at most d and restart.
        return controlled(((Q(1,k),pay),),ret(d,k+1),Q(0),k-1+Q(1,d),True,
                          [0]*(steps+1),extra=(d,pay))

    @cache
    def V(d,j):
        if d>=dcut:
            return poisson(j)
        if j==0:
            return ret(d,2)
        rise = V(d+1,j+1)
        return controlled(((Q(j,d),V(d,j-1)),),rise,Q(j),Q(d),True,poisson(j))

    @cache
    def L(k):
        return ret(1,k)[0] if k<rcut else 0

    @cache
    def Host(k):
        if k>=rcut:
            return V(2,1)
        return controlled(((Q(1)+Q(1,k),V(3,2)),(Q(1,2),V(2,0))),Host(k+1),
                          Q(0),k-Q(1,2),True,V(2,1))

    if upper:
        return U
    return V,L,Host


def evaluate(lo,hi,err=Q(1,1000),dcut=24,rcut=32,steps=128,scale=Q(33,32)):
    assert sum(Q((-1)**j,factorial(j)) for j in range(14))>Q(367879441,S)
    U0 = solve((1-err)*lo,True,dcut,rcut,steps,scale)
    V,L0,Host = solve((1+err)*hi,False,dcut,rcut,steps,scale)
    M = Q(V(2,1)[0],S)-err

    @cache
    def U(k,j):
        if k>=rcut:
            z = (1-err)*lo
            cap = Q(exp_bounds(z)[1],S)*root_upper(z,k)/(z*(k-1))
            return min(Q(1,j+1),cap+err)
        return min(Q(1,j+1),Q(U0(k,j)[0],S)+err)

    @cache
    def L(k):
        return max(0,Q(L0(k),S)-err)

    # The frozen fractional checker supplies a purely algebraic ledger.
    # lo=0 invokes only the stock cap on high endpoints, valid for sharp tails.
    B = FR['Bounds']
    bonuses = {k:Q(Host(k)[0]-V(2,1)[0],S) for k in range(2,min(rcut,8))}
    v = SimpleNamespace(stop=max(48,rcut),lo=Q(0),g=M-U(2,0),bonuses=bonuses,U=U,L=L)
    v.profile = lambda d,c,r:B.profile.__wrapped__(v,d,c,r)
    print(f'kernels on [{lo},{hi}]: M={M}, U(2,0)={U(2,0)}, L(2)={L(2)}',flush=True)
    for t in (('.011','.011','.007','.020','.002'),('.011','.011','.007','.020','.003'),('.010','.011','.007','.020','.002'),
              ('.011','.012','.007','.020','.001'),('.011','.012','.008','.020','.001'),
              ('.012','.012','.008','.022','.003'),('.015','.015','.007','.022','.005')):
        a,b,fee,ft,tr = map(Q,t)
        out = B.ledger(v,a,b,fee,ft,160,tr)
        if out[0]:
            print('PASS signed ledger',t,flush=True)
            # Chernoff: each relative deviation has probability <= exp(-m e²/(2(1+e))).
            target = 1
            while exp_bounds(Q(target))[1] > S*err:
                target += 1
            m = -((-(2*(1+err)*target/err**2))//1)
            print(f'PASS probability tails for all m >= {m}',flush=True)
            return t
        print('FAIL tariff',t,out[1],flush=True)
    return None


def finite():
    """Lower finite seam via the frozen interval kernel and b=20 tariffs."""
    A = pinned(HERE.parent/'b20.py','ea2f9e71a1700badc175c71f788aefbded70c0b3')
    B = A['A']
    # Verify constants needed by the imported analytical caps and Poisson floor.
    atan = lambda r,n:sum((-1)**i*r**(2*i+1)/(2*i+1) for i in range(n))
    assert 16*atan(Q(1,5),8)-4*atan(Q(1,239),3)>Q(314159,100000)
    assert sum(Q((-1)**i,factorial(i)) for i in range(10))>B['M']
    for b in range(22,25):
        lo,hi = Q(11,5*(b+1)),Q(1,10)
        if lo>=hi:
            continue
        cells = B['certify'](b,lo,hi,14,10)[0]
        print(f'PASS b={b}, {cells} interval box(es), y >= 11/5 up to x=1/10',flush=True)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__)
    p.add_argument('lo',type=Q,nargs='?')
    p.add_argument('hi',type=Q,nargs='?')
    p.add_argument('--finite',action='store_true')
    p.add_argument('--err',type=Q,default=Q(1,1000))
    p.add_argument('--dcut',type=int,default=24)
    p.add_argument('--rcut',type=int,default=32)
    p.add_argument('--steps',type=int,default=128)
    p.add_argument('--scale',type=Q,default=Q(33,32))
    a = p.parse_args()
    if a.finite:
        finite()
    else:
        assert 0<a.lo<=a.hi and 0<a.err<1 and a.dcut>=4 and a.rcut>=4 and a.steps>=1
        assert 1<a.scale<=Q(33,32)
        assert evaluate(a.lo,a.hi,a.err,a.dcut,a.rcut,a.steps,a.scale)
