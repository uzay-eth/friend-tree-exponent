"""Outward rational threshold boxes for physical signed pre-placement order."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from math import factorial, prod
from pathlib import Path
from runpy import run_path
from types import FunctionType

EXT = Path(__file__).resolve().parents[1] / 'annealed/extension'
for name, pin in (
    ('third-physical/check.py', 'e39f8bfdf66f19e34a59247fa3649ccaf884358b'),
    ('first-sender/joint/check.py', 'd3b4ad5bc87760a843712d9af43af76450a07a37'),
    ('first-sender/check.py', '5680c8e937ebafdec09c2e193e6b36fe47ac84b1'),
):
    data = (EXT / name).read_bytes()
    assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest() == pin, name
PHYSICAL = run_path(str(EXT / 'third-physical/check.py'))
JOINT = run_path(str(EXT / 'first-sender/joint/check.py'))
S = PHYSICAL['S']
M = Q(36787, 100000)


@cache
def root_bound(k, q):
    lo, hi = S, 2*S
    assert q >= Q(1, 2)
    while hi-lo > 1:
        j = (hi+lo)//2
        if q*j**k >= S**k:
            hi = j
        else:
            lo = j
    return Q(hi, S)


@cache
def beta(n, k, qlo, qhi, promoted):
    r = Q(1, k)
    z = Q(314159, 100000)*r
    normal = r*(1-z*z/6+z**4/120)*factorial(n-1)
    if promoted:
        normal /= r*prod(i-r for i in range(1,n))
        a, q = 1+r, qlo
    else:
        normal /= prod(i-r for i in range(1,n+1))
        a, q = r, qhi
    B = n-a
    coeff = value = Q(1)
    for j in range(1,49):
        coeff *= q*(j-a)/j
        value += coeff*B/(B+j)
    if not promoted:
        value += coeff*B/(B+49)*q/(1-q)
    return PHYSICAL['up'](normal*root_bound(k,q)*(1 if promoted else q)*value)


def tips(b, xlo, xhi, depth=14):
    qlo, qhi = 1-xhi, 1-xlo
    m, D = b+1, b+1+depth
    row = {(k,j):beta(D,k,qlo,qhi,False) for k in range(2,D) for j in range(k)}
    for n in range(D-1,m-1,-1):
        nxt = {}
        for k in range(2,n):
            cap = beta(n,k,qlo,qhi,False)
            prom = beta(n+1,k,qlo,qhi,True)
            for j in range(k):
                gain, hold = row[k+1,j+1], row[k,j]
                loss = row[k,j-1] if j else 0
                y = prom+j*loss
                z, zz = k*(j+1), k*(k+j+1)
                small = y+z*gain+(k*n-z-j-1)*hold
                large = 2*y+zz*gain+(2*k*n-zz-2*(j+1))*hold
                val = qhi.numerator*max(2*small,large)
                den = qhi.denominator*2*k*n
                nxt[k,j] = min(cap,(val+den-1)//den)
        row = nxt
    beta.cache_clear()
    return {key:Q(v,S) for key,v in row.items()}


def ledger(b, U, fee=Q('.004'), ft=Q('.008'), a=Q('.004'), base=Q('.004')):
    source = PHYSICAL['ledger']
    test = FunctionType(source.__code__, dict(source.__globals__, tips=lambda _:U))
    return test(b,M,fee,ft,a,base)


def favourable(b, xhi, depth):
    """Lower pending-tip tail by gaining a leaf before the finite horizon."""
    q = 1-xhi
    m, D = b+1, b+1+depth
    row = [0]*D
    for n in range(D-1,m-1,-1):
        nxt = [0]*n
        for k in range(2,n):
            gain = int(M*S)+(k*n-1)*row[k]
            # Promotion is disjoint from parent gain; g in [0,k].
            assert k*k <= k*n-1
            gain += min(0,k*k*(row[k+1]-row[k]))
            nxt[k] = q.numerator*gain//(q.denominator*k*n)
        row = nxt
    return [Q(v,S) for v in row]


def improvement(b, U, L):
    """Signed certificate by local baseline/excess, retaining opposite credit."""
    from itertools import combinations_with_replacement, product

    g = M-U[2,0]
    assert g > 0

    def profile(d,c,r):
        l = d-c
        s,t = (l+r)/(d*(d+1)),Q(1,2*d)
        ans = (d*s+l*t)*g
        if l:
            ans -= l*(max(0,t-s)*max(0,U[2,0]-L[d+1])
                      +max(0,s-t)*max(0,U[d+1,l-1]-L[2]))
        return s,ans

    def site(a, base, fee, ft, transfer):
        if any(profile(d,0,Q(0))[1] <= 0 for d in range(2,b-1)):
            return ('star',)
        for d in range(2,b-1):
            for c in range(1,d+1):
                if d==c==2:
                    continue
                F = ft if d==2 else fee
                adverse = U[d+1,d-c]
                debit = transfer if d>=8 else 0
                for k in range(c+1):
                    def reserve(r):
                        s,C = profile(d,c,r)
                        return s,C+a*(c-2)-base-k*F-(c-k)*debit
                    if all(R >= (c-k)*max(0,s*adverse-base)
                           for s,R in (reserve(Q(k,2)),reserve(Q(c,2)))):
                        continue
                    if d>8:
                        return ('coarse site',d,c,k)
                    # In the high bin retain its finite reciprocal interval
                    # and minimum opposite coefficient, but drop its payout.
                    for es in combinations_with_replacement((0,*range(2,min(8,b-1))),c-k):
                        for rhi in (Q(1,b-2),Q(1,8)):
                            r = Q(k,2)+sum(Q(1,e) if e else rhi for e in es)
                            s,R = reserve(r)
                            R += es.count(0)*transfer if d<=3 else 0
                            excess = sum(max(0,max(0,s-Q(1,d*(e or b-2)*((e or b-2)+1)))
                                             *max(0,adverse-(L[e+1] if e else 0))-base) for e in es)
                            if R < excess:
                                return ('site',d,c,k,es,rhi,R-excess)
        return None

    def chain(base,fee,ft):
        # Integer outward enclosures of the signed chain debit.
        X,Y = 6*10**8,10**9
        ceil = lambda z:-((-z)//1)
        end = set()
        for d in range(2,b-1):
            for c in range(1,d+1):
                if d==c==2:
                    continue
                for r in (Q(1,2),Q(c,2)):
                    s,_ = profile(d,c,r)
                    end.add((d,Q(X,6*d)//1,ceil(Q(X,6*d)),s*X//1,ceil(s*X),
                             ceil(Y*max(0,U[3,0]-L[d+1])),
                             ceil(Y*max(0,U[d+1,d-c]-L[3])),
                             (ft if d==2 else fee)*X*Y//1))
        v = ceil(Y*max(0,U[3,0]-L[3]))
        credit = 2*g*Y//1

        @cache
        def port(a,b,e):
            _,_,_,s,t,w,y,f = e
            return max(0,(b-s)*w,(t-a)*y)-f

        for length in (1,2,3):
            for x,y in product(end,repeat=2):
                if x[0]+y[0]+length>b-1:
                    continue
                if length==1:
                    ss = ((x[1]+y[1],x[2]+y[2]),)
                else:
                    ss = ((X//12+x[1],X//12+x[2]),)+((X//6,X//6),)*(length-2)+((X//12+y[1],X//12+y[2]),)
                debt = port(*ss[0],x)+port(*ss[-1],y)
                debt += v*sum(max(t-a,b-s) for (a,b),(s,t) in zip(ss,ss[1:]))
                debt -= credit*sum(a for a,_ in ss)
                if debt > base*X*Y:
                    return ('chain',length,x,y,Q(debt,X*Y)-base)
        return None

    return site,chain


def obstruction():
    # Twenty old vertices: nine length-two spokes plus one ordinary leaf.
    x = Q(1,4)
    U, L = tips(21,x,x), favourable(21,x,24)
    hub, arm, tip, leaf = Q(11,220), Q(11,60), Q(1,4), Q(1,20)
    A = 10*hub+18*arm+9*tip+leaf
    # Exact coefficients at both sides of all edges; opposite payout is
    # determined by the actual parent degree in that parallel completion.
    debt = (9*(arm-hub)*max(0,U[3,1]-L[11])
            +9*(tip-arm)*max(0,U[2,0]-L[3])
            +abs(leaf-hub)*max(0,U[2,0]-L[11]))
    score = A*(M-U[2,0])-debt
    need = U[2,0]+debt/A
    assert A*(Q(211,1000)-U[2,0])-debt < 0 < score
    print(f'CONTROL b=21, x=1/4, 9 two-edge spokes and 1 leaf: score={score}; serial coefficient needed > {need}',flush=True)


def certify(b, lo, hi, depth, splits):
    U = tips(b,lo,hi,depth)
    ok, where = ledger(b,U)
    if ok:
        return 1, where
    site, chain = improvement(b,U,favourable(b,hi,24))
    for fee,ft in ((Q('.004'),Q('.008')),(Q('.005'),Q('.009'))):
        where = site(Q('.004'),Q('.004'),fee,ft,Q('.002')) or chain(Q('.004'),fee,ft)
        if where is None:
            return 1, f'paired reserve; bills {fee}, {ft}'
    assert splits>0,(b,lo,hi,where)
    mid = (lo+hi)/2
    left,_ = certify(b,lo,mid,depth,splits-1)
    right,_ = certify(b,mid,hi,depth,splits-1)
    return left+right, 'partitioned'


def analytic(xlo):
    # Joint tail debit; the sixth-root bound also covers every b>=B.
    for B in range(21,101):
        lo,hi = 0,20000
        ratio = (1-xlo)/(B*xlo)
        assert ratio <= 2**6
        while hi-lo > 1:
            mid = (lo+hi)//2
            if ratio <= Q(mid,10000)**6:
                hi = mid
            else:
                lo = mid
        r = Q(hi,10000)
        if M*M*Q(157,50)*(B+Q(1,2))*xlo > (Q(4,3)+r/2)**2*(1-xlo):
            JOINT['corner'](B,xlo,M,r)
            return B,r
    raise AssertionError(xlo)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--lo',type=Q,default=Q(1,10))
    parser.add_argument('--hi',type=Q,default=Q(1,3))
    parser.add_argument('--start',type=int,default=21)
    parser.add_argument('--stop',type=int)
    parser.add_argument('--depth',type=int,default=14)
    parser.add_argument('--splits',type=int,default=8)
    args = parser.parse_args()
    assert 0 < args.lo <= args.hi <= Q(1,2)
    assert 20 <= args.start and (args.stop is None or args.start <= args.stop)
    assert args.depth >= 1 and args.splits >= 0
    atan = lambda r,n:sum((-1)**i*r**(2*i+1)/(2*i+1) for i in range(n))
    assert 16*atan(Q(1,5),8)-4*atan(Q(1,239),3)>Q(314159,100000)
    gamma = sum((-1)**k*Q(3,3*k+1)/factorial(k) for k in range(6))
    assert Q(22,7)<Q(1773,1000)**2
    assert Q(1773,1000)/gamma*Q(42,41)<Q(3,4)
    assert sum(Q((-1)**k,factorial(k)) for k in range(10)) > M
    print('PASS Poisson serial constant 1/e > 36787/100000',flush=True)
    B,r = analytic(args.lo)
    print(f'PASS analytic b>={B}, x>={args.lo}, root <= {r}',flush=True)
    for b in range(args.start, min(args.stop or B,B)):
        cells, where = certify(b,args.lo,args.hi,args.depth,args.splits)
        print(f'PASS b={b}: {cells} interval box(es); {where}',flush=True)
    obstruction()


if __name__ == '__main__':
    main()
