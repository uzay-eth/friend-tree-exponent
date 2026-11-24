"""Exact upper replay of the one-family plus conditional-Markov certificate."""
from fractions import Fraction as Q
from pathlib import Path
import runpy

base = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'check.py'))
floor, ceil = (base[k] for k in ('floor', 'ceil'))
S = 10**8


def residual_lower(b, f):
    # C=b*sum_{j>f} 1/[(b+j)(j+1)] exactly by telescoping.
    return Q(b,b-1)*sum((Q(1,k) for k in range(f+2,f+b+1)),Q())


def certificate_ceiling(b=21, n=512, d=24, end=12):
    # q gives F=floor((b-1)(1-q)/q); t has exact Pareto survival.
    weights = []
    old = Q(1)
    for i in range(1,d*end+1):
        survival = (1+Q(i,d*b))**(-b)
        weights.append(ceil(S*(old-survival)))
        old = survival
    total_hit = total_miss = 0
    for j in range(1,n):
        lo = (b-1)*(n-j-1)//(j+1)
        hi = (b-1)*(n-j)//j
        C = floor(S*residual_lower(b,hi))
        # E reciprocal times expected non-live stock, conditional on F,t.
        # Bound the non-live stock factor at hi and reciprocal at lo.
        M = min(S,ceil(S*(1+Q(hi,2*(b-1)))*min(Q(1,2),Q(2,max(1,lo)))))
        # The uncapped W gives exactly the same payoff below the current t.
        family = []
        p = 0
        for k in range(1,d*end+1):
            cdf = ceil(S*(1-(1+Q(k,d*b))**(-hi-1)))
            family.append(cdf-p)
            p = cdf
        for i,weight in enumerate(weights,1):
            miss_pay = max(0,S-C*d//i)
            total_miss += weight*M*miss_pay
            hit = 0
            for k,mass in enumerate(family[:i],0):
                # Left endpoint family loss and a lower rounded penalty.
                loss = C*(b*d+k)//(b*(i-k))
                if loss>=S:
                    break
                hit += mass*(S-loss)
            total_hit += weight*hit
    return Q(total_hit,n*S**3),Q(total_miss,n*S**3),Q(1,n)+old


if __name__ == '__main__':
    for f in (0,1,5,50):
        z = Q(1,2)
        for h in range(f):
            z *= 1+Q(1,20+h)
        assert z+Q(1,2) == 1+Q(f,40)
    a,m,t = certificate_ceiling()
    assert a<Q(232,1000) and m<Q(38,1000) and t<Q(21,10000)
    assert a+m+t<Q(272,1000)
    # The chosen debit, not the physical adverse tail, exceeds .28.
    assert Q(2,3)**6<Q(2,21)
    assert Q(28,100)**2*Q(22,7)*Q(43,4)<Q(5,3)**2
    print('PASS exact upper averaging: hit < .232, miss credit < .038, discarded bins < .0021')
    print('PASS one-family certificate < .272 versus chosen beta debit > .28 at (21,1/3)')
