#!/usr/bin/env python3
"""Exact whole-tree generation polynomials and fixed-discount prime tests."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import gcd, isqrt, lcm

EDGE = (
    Q(191976090122981726641370633087, 339373453939523193250000000000),
    Q(24172191945911554486185483421, 42724375658564528593750000000),
)
MEMORY = (
    Q(2586244278310822397274177587149, 4571936719109803018500000000000),
    Q(6574835699364880591816419305857, 11621030179129551777500000000000),
)
TRIAL = Q(99, 175)


def v(x, p):
    assert x
    x = Q(x)
    a, b, k = x.numerator, x.denominator, 0
    while a % p == 0:
        a //= p
        k += 1
    while b % p == 0:
        b //= p
        k -= 1
    return k


def mod(x, p):
    x = Q(x)
    return x.numerator * pow(x.denominator, -1, p) % p


def bracket():
    assert EDGE[0] < MEMORY[0] < MEMORY[1] < EDGE[1]
    for lo, hi in (EDGE, MEMORY):
        assert Q(56, 99) < lo < TRIAL < hi < Q(43, 76)
        assert 99 * 43 - 56 * 76 == 1
        assert all(not lo <= Q(a, b) <= hi for b in range(1,175)
                   for a in range(1,b))
        hits = []
        for a in range(1,8):
            for b in range(-7,8):
                for c in range(-7,8):
                    d = b*b-4*a*c
                    if gcd(a,b,c) != 1 or d<=0 or isqrt(d)**2 == d:
                        continue
                    # All vertices for this small-height search lie outside the bracket.
                    assert not lo <= Q(-b,2*a) <= hi
                    f = lambda x: (a*x+b)*x+c
                    if f(lo)*f(hi) <= 0:
                        hits.append((max(a,abs(b),abs(c)),a,b,c))
        assert min(hits) == (7,3,-7,3)
        hits = []
        for d in range(1,isqrt(57200)+1):
            for n in range(d+1,57200//d+1):
                if gcd(n,d)==1 and lo*n*d <= n*n-d*d <= hi*n*d:
                    hits.append((n*d,n*n-d*d))
        assert min(hits) == (57200,32361)
    print('EDGE and post-EDGE minimums: rational 99/175; quadratic 3s²−7s+3; square case 32361/57200')


def flatten(tree):
    adj, depths = [], []

    def visit(t, parent, depth):
        i = len(adj)
        adj.append([] if parent is None else [parent])
        depths.append(depth)
        for child in t:
            adj[i].append(visit(child,i,depth+1))
        return i

    visit(tree,None,0)
    return adj, depths


@cache
def row(tree):
    adj, depths = flatten(tree)
    D = len(adj[0])+1
    rates = defaultdict(Q)
    for donor, nb in enumerate(adj):
        for target in nb:
            rates[target] += Q(1,len(nb)+(donor==0))
    assert sum(rates.values()) == len(adj)-Q(1,D)
    edges=defaultdict(Q)
    for target,r in rates.items():
        def encode(i,parent):
            kids=[encode(j,i) for j in adj[i] if j!=parent]
            if i==target:
                kids.append(())
            return tuple(sorted(kids))
        depth = depths[target] if target and len(adj[target])==1 else 0
        edges[encode(0,None),depth]+=r
    return D,tuple(edges.items())


def add(a,b,scale=1):
    if len(a)<len(b):
        a.extend([Q(0)]*(len(b)-len(a)))
    for i,x in enumerate(b):
        a[i]+=x*scale


def mul(a,b):
    c=[]
    for i,x in enumerate(a):
        add(c,[Q(0)]*i+b,x)
    return c


def evaluate(poly,s):
    ans=Q(0)
    for coef in reversed(poly):
        ans=ans*s+coef
    return ans


def primitive(p):
    h=lcm(*(x.denominator for x in p))
    ints=[int(h*x) for x in p]
    g=gcd(*ints)
    return [x//g for x in ints],Q(g,h)


def divide(a,b):
    a=list(a)
    while a and not a[-1]:
        a.pop()
    q=[Q(0)]*max(0,len(a)-len(b)+1)
    while len(a)>=len(b):
        k=len(a)-len(b)
        q[k]=a[-1]/b[-1]
        for j,x in enumerate(b):
            a[k+j]-=q[k]*x
        while a and not a[-1]:
            a.pop()
    return q,a


def reduce_holdings(num, den):
    # Denominator roots are known rational holding poles; this routine
    # verifies reduction with a monic Euclidean gcd over Q[s].
    a,b=den,num
    while b:
        b=[x/b[-1] for x in b]
        _,r=divide(a,b)
        a,b=b,r
    p,r=divide(num,a)
    assert not r
    d,r=divide(den,a)
    assert not r and d[-1]==1
    return p,d


def polynomials(last):
    # P_n(s,z) counts discounted births through size n, including initial
    # root. At each size use a common polynomial denominator in s.
    flow={((),):[Q(1)]}
    old=[[Q(1)]]
    den=[Q(1)]
    for n in range(2,last+1):
        Ds={row(t)[0] for t in flow}
        factors={d:[n-Q(1,d),Q(1)] for d in Ds}
        def product(ds):
            a=[Q(1)]
            for d in sorted(ds):
                a=mul(a,factors[d])
            return a
        complements={d:product(Ds-{d}) for d in Ds}
        step=product(Ds)
        den=mul(den,step)
        nxt=defaultdict(list)
        delta=[[] for _ in range(n)]
        for t,poly in flow.items():
            d,edges=row(t)
            mass=mul(poly,complements[d])
            for (y,dep),r in edges:
                add(nxt[y],mass,r)
                if dep:
                    add(delta[dep],mass,r)
        old=[mul(p,step) for p in old]+[[] for _ in range(len(delta)-len(old))]
        for i,p in enumerate(delta):
            add(old[i],p)
        yield n,den,delta,old
        flow=nxt


def literal(last,s):
    # Unmerged vertex-ID walks; independent of the canonical branching code.
    flow=[([(-1,1),(0,)], [0,1], Q(1))]
    for n in range(2,last+1):
        out=[Q(0)]*n
        nxt=[]
        for a,depth,weight in flow:
            mass=weight/(s+n-Q(1,len(a[0])))
            for u,nb in enumerate(a):
                for t in nb:
                    if t<0:
                        continue
                    r=mass/len(nb)
                    if t and len(a[t])==1:
                        out[depth[t]]+=r
                    if n<last:
                        nxt.append((a[:t]+[a[t]+(n,)]+a[t+1:]+[(t,)],depth+[depth[t]+1],r))
        yield out
        flow=nxt


def generation():
    literal_rows=list(literal(6,TRIAL))
    for n,den,diff,total in polynomials(8):
        D=evaluate(den,TRIAL)
        inc=[evaluate(p,TRIAL)/D for p in diff]
        coeff=[evaluate(p,TRIAL)/D for p in total]
        assert coeff[0]==1
        if n<=6:
            assert inc==literal_rows[n-2]
        if n in (5,7):
            for a in (TRIAL,Q(2,n),Q(1,n**2),Q(1,n)):
                e=-v(a,n)
                values=[evaluate(poly,a)/evaluate(den,a) for poly in diff]
                u=mod(n**e*a,n)
                if e==1 and u==1:
                    w=v(n*(a+n)-1,n)
                    power=(n-2)-w
                    factor=1
                else:
                    power=(n-1)*e-1
                    factor=1 if e>1 else u*pow(u-1,-1,n)%n
                scaled=[mod(c/Q(n)**power,n) for c in values]
                assert scaled==[0,(-factor)%n]+[0]*(len(inc)-2)
                if a==TRIAL:
                    print(f'physical cap {n}→{n+1}: p={n}, e={e}, scaled generation increment {scaled}')
        if n in (5,6,7,8):
            # Actual successive depth coefficients, not a Laurent coefficient.
            g=total[2].copy()
            add(g,total[1],-1)
            raw_deg=len(g)-1
            g,d=reduce_holdings(g,den)
            g,content=primitive(g)
            deg=len(g)-1
            hom=evaluate(g,TRIAL)*TRIAL.denominator**deg
            assert hom.denominator==1
            gain=1
            for p in (5,7):
                rev=tuple(x%p for x in reversed(g))
                lead_lost=next(i for i,x in enumerate(rev) if x)
                value=evaluate(g,TRIAL)
                assert value*content/evaluate(d,TRIAL)==coeff[2]-coeff[1]
                assert not any(x%p for x in g[deg-lead_lost+1:])
                h=v(hom,p)
                gain*=p**h
                print(f'size {n+1} depth 2−1: p={p} degree={deg} (raw {raw_deg}) leading zeros={lead_lost}, '
                      f'v(homogeneous g)={h}, v(depth gap)={v(coeff[2]-coeff[1],p)}')
            floor=gain*abs(content/evaluate(d,TRIAL))/TRIAL.denominator**deg
            assert 0<floor<abs(coeff[2]-coeff[1])
            print(f'  finite nonzero-form bound improves by {gain}; true gap / bound > 10^{digits_ratio(abs(coeff[2]-coeff[1])/floor)}')
            quadratic_resultant(g)


def digits_ratio(x):
    k=0
    while x>10:
        x/=10
        k+=1
    return k


def quadratic_resultant(g):
    # a is the independent height-seven witness, 3a²−7a+3=0.
    # Reduce g(a)=x+y*a and take its exact field norm.
    x=y=Q(0)
    for coef in reversed(g):
        x,y=coef-y,x+Q(7,3)*y
    d=len(g)-1
    result=3**d*(x*x+Q(7,3)*x*y+y*y)
    assert result and result.denominator==1
    # Homogeneities of both leading coefficients share infinity mod 3.
    shared_infinity=int(g[-1]%3==0)
    assert v(result,3)>=shared_infinity
    print(f'  height-7 quadratic: v3(resultant)={v(result,3)} '
          f'(common infinity test gives at least {shared_infinity})')


def counterexample():
    # Infinite 1x1 generation kernel: exact crossing remains our rational
    # trial, while physical-style depth gaps can keep a p-unit obstruction
    # at every finite size. Each retained holding chain has positive rates.
    p=5
    def full(s):
        return Q(TRIAL+2)/(s+2)
    def killed(s,n):
        return (1-Q(1,p**n))*full(s)
    for n in range(1,13):
        k=killed(TRIAL,n)
        assert v(k,p)==-n
        for m in (1,2,3):
            d=k**(m+1)-k**m
            assert v(d,p)==-n*(m+1)
            assert mod(p**(n*(m+1))*d,p)==(-1)**(m+1)%p
        assert killed(TRIAL,n+1)-k==Q(p-1,p**(n+1))
    assert full(TRIAL)==1
    print('exact rational crossing: every finite depth difference has a p-unit after scaling; size limit is not p-adic')

    # A second nested cutoff converges in BOTH places, but to distinct
    # values. All depth gaps have a nonzero limiting residue as well.
    theta=Q(2,3)
    t=lambda n: theta+(1-theta)*Q(p**n,p**n+1) if n else Q(0)
    for n in range(1,13):
        ts=[t(i)-t(i-1) for i in range(1,n+1)]
        assert min(ts)>0 and sum(ts)==t(n) and v(t(n)-theta,p)==n
        assert t(n)*full(TRIAL)==t(n) and 1-t(n)==(1-theta)/Q(p**n+1)
        assert mod(1-t(n)*full(TRIAL),p)==2
        for m in (1,2,3):
            assert mod(t(n)**(m+1)-t(n)**m,p) == 3*4**m%p != 0
    print('both-place convergence: real K(99/175)=1, 5-adic K(99/175)=2/3; every limiting depth gap is a 5-unit')


if __name__=='__main__':
    bracket()
    generation()
    counterexample()
