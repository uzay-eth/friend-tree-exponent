"""Independent exact boundary, signed debit and cutoff checks."""

from fractions import Fraction as Q
from itertools import product
from math import factorial


def harmonic():
    rows=0
    for x in (Q(0),Q(1,9),Q(3,13),Q(6,9),Q(31,32)):
        pmf=[Q(1)]
        for m in range(2,25):
            def step(p):
                return [a*(1-x)+b*x for a,b in zip(p+[0],[0]+p)]
            pmf=step(pmf)
            next_pmf=step(pmf)
            def cdf(v,r):
                return sum(v[:max(0,r)],Q(0))
            H=lambda v,j:sum((cdf(v,j-k)/factorial(k) for k in range(j+2)),Q(0))
            assert cdf(pmf,-1)==cdf(pmf,0)==0 and cdf(pmf,m)==1
            for j in range(1,m):
                terms=[(k*(cdf(next_pmf,j-k+1)-cdf(next_pmf,j-k))
                        -(cdf(next_pmf,j-k)-cdf(next_pmf,j-k-1)))/factorial(k)
                       for k in range(j+3)]
                assert terms[j]==cdf(next_pmf,1)/factorial(j-1)>0
                assert not any(terms[j+1:]) and sum(terms)==0
                assert m*H(pmf,j)==j*H(next_pmf,j+1)+(m-j-1)*H(next_pmf,j)+H(next_pmf,j-1)
                rows+=1
    print(f'PASS {rows} binomial PMF harmonic rows, including x=0 and the k=j edge')


def weights():
    rows=tight=0
    for n in range(2,8):
        for parents in product(*(range(i) for i in range(1,n))):
            adj=[set() for _ in range(n)]
            for i,a in enumerate(parents,1):
                adj[i].add(a); adj[a].add(i)
            d=list(map(len,adj))
            r=[sum(Q(1,d[a]) for a in v) for v in adj]
            w=[r[i]/(d[i]+1) for i in range(n)]
            A=sum(w)
            assert sum(r)==n and A>=1
            X=E=Q(0)
            for u,v in ((u,v) for u in range(n) for v in adj[u] if u<v):
                forward=r[u]*(r[v]-Q(1,d[u]*(d[u]+1)))
                reverse=r[v]*(r[u]-Q(1,d[v]*(d[v]+1)))
                diff=forward-reverse
                assert diff==w[v]/d[v]-w[u]/d[u]
                E+=abs(diff)
                if diff and d[u if diff<0 else v]==1:
                    X+=abs(diff)
            assert E<=A and 3*X<=A
            tight+=3*X==A
            rows+=1
    print(f'PASS {rows} exact paired charges, {tight} saturating the one-third bound')


def cutoff():
    gamma=sum(Q(3*(-1)**k,factorial(k)*(3*k+1))*Q(3,2)**(3*k+1) for k in range(8))
    assert gamma>Q(13,5)
    # Alternating arctangent bounds, with one extra positive term for the upper.
    pi=16*(Q(1,5)-Q(1,3*5**3)+Q(1,5**6)-Q(1,7*5**7))-Q(4,239)
    assert 6*pi*Q(231,1000)**2>1 and 36*Q(3031,10000)**3>1
    result=Q(0)
    for s,pre,w in ((Q(1,2),Q(231,1000),Q(4,3)),(Q(2,3),Q(3031,10000)/Q(13,5),Q(2,3))):
        D=Q(7,9)*s*(1-(1+s)/6)
        assert s>(s+1)*D/6
        result+=w*pre*(1-D/6)
    residual=sum(Q((-1)**k,factorial(k)) for k in range(12))-result
    assert residual>Q(19,10000)
    assert Q(6,20+6)==Q(3,13) and Q(6,21+6)==Q(2,9)
    print(f'PASS exact cutoff surplus > {residual} (>19/10000), including b=3')


if __name__=='__main__':
    harmonic()
    weights()
    cutoff()
