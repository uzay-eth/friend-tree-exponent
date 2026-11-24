"""Add an exact bare-core corrector; reuse pinned f/g floor receipts, no cross-scan."""
from fractions import Fraction as Q
from pathlib import Path
from collections import defaultdict
import json
import re

HERE=Path(__file__).resolve().parent
E=Q(7,10**7)
C0={'lower':Q(184124847315,10**11),'upper':Q(195820675182,10**11)}
M={'lower':Q(25227936004759,1700000000000000000000),
   'upper':Q(5046841434259,340000000000000000000)}


def physical_check():
    def graph(p):
        g=[[] for _ in p]
        for j,i in enumerate(p[1:],1): g[j].append(i);g[i].append(j)
        return g
    def B(p):
        g=graph(p)
        return sum(len(g[i])>=2 and all(len(g[j])>1 for j in g[i]) for i in range(1,len(p)))
    generation=[(-1,0,0)]; checked=0; births=deaths=0
    for _ in range(4):
        nxt=[]
        for p in generation:
            g=graph(p); rate=defaultdict(Q)
            for js in g:
                for j in js: rate[j]+=Q(1,len(js))
            exact=sum(r*(B(p+(j,))-B(p)) for j,r in rate.items())
            plus=minus=0
            for i in range(1,len(p)):
                if len(g[i])<2: continue
                L=sum(len(g[j])==1 for j in g[i])
                if L==1: plus+=Q(1,len(g[i]))
                if not L: minus+=sum(Q(1,len(g[j])) for j in g[i])
            assert exact==plus-minus
            births+=bool(plus); deaths+=bool(minus); checked+=1
            nxt.extend(p+(j,) for j in rate)
        generation=nxt
    assert births and deaths
    print(f'bare-core correction physical arrows: {checked} trees; birth exercised {births}, death exercised {deaths}')


def main():
    physical_check()
    assert 0<E<Q(113465364604,4*10**11)
    data=json.loads((HERE/'charged-f.json').read_text())
    receipt=(HERE/'combined.log').read_text()
    for name in ('lower','upper'):
        d=data[name]; s=d['sign']; lam=Q(d['lambda']); m=M[name]
        val={tuple(row['state']):min(map(Q,row['signed_corners'])) for row in d['rows']}
        a,b=val[1,0,1],val[0,1,1]
        ap,bp=a-E/2,b+E*lam
        A=lambda k:ap-(ap-bp)*k
        assert ap-bp>=m and -A(5)>=6*m
        old=lambda k:a-(a-b)*k
        count=0
        for (L,k,h),margin in val.items():
            D=L+k+h
            if D>=18: continue
            worst=E*lam if L==0 else -E/Q(D) if L==1 else 0
            assert margin+worst>=A(k)
            count+=1
        # Previous exact tail-only minima are reported as surplus over old(k).
        for k in range(5):
            t=Q(re.search(rf'^{name} tail coverage {k} min surplus (\S+)$',receipt,re.M)[1])
            # Covers all tail degrees, including capped tails. Even if L=1,
            # the new bare term costs at most E/18; discard positive L=0 correction.
            assert t+old(k)-A(k)-E/18>=m*(k+1)
        # For k>=5 the same fallback starts at k=5; include the E/18 tail loss,
        # and E/(k+2) bounds the worst inner loss (at L=1).
        assert -A(5)-E/18>=6*m and -A(5)>=E/7
        rate=lam+s*(bp+m)/(C0[name]+(E if s<0 else 0))
        print(name,'inner states',count,'a_prime',ap,'b_prime',bp,'endpoint',rate,
              f'{float(rate):.16f}','additional inward shift',float(s*(rate-(lam+s*(b+m)/C0[name]))))


if __name__=='__main__':main()
