#!/usr/bin/env python3
"""Check local endpoint failures and the smooth seam without importing a solver."""
from fractions import Fraction as F
from uncapped import COARSE, analytic_tail, check, exact_inputs


def main():
    d=F(17,50)
    n=80
    a=F(11299,20000)
    P=exact_inputs(n,d,COARSE)
    check(n,d,a,P,COARSE)
    try:
        check(n,d,a,P,COARSE,'wrong-sender')
    except ValueError as e:
        assert 'edge corner' in str(e)
        print('rejected the wrong grandchild sender rate')
    else: raise AssertionError('wrong-sender escaped')

    M=160
    P=exact_inputs(M,d)
    C=analytic_tail(d)
    def edge(i): return P[i] if i<M else C/i
    rows=0
    for i in [2,3,5,15,M-1,M,M+1,2*M,10000]:
        for j in [2,3,15,M-1,M,M+1,2*M,10000]:
            q=lambda i: 1/(1+d*i)
            delta=lambda i: (q(i)-q(i+1))/d
            for u in (d,d+F(1,2)):
                for r in (F(1,j),i-2+F(1,j)):
                    # Unfactored coefficient of a grandchild score.
                    gain=(r*(edge(i+1)-edge(i))*edge(j)
                          +edge(i)*(edge(j+1)-edge(j))/i)
                    emitted=u*((q(i)+r*delta(i))*edge(j)+delta(j)/i)
                    assert gain+emitted<=u*edge(i)*edge(j)
                    rows+=1
    print('exact uncut two-level drift probes:',rows)

    # An actual lower normalization bound (no outward rounding upward).
    # A common denominator for these acyclic rows is unnecessary: keeping
    # each floor on the integer grid gives a certified lower table.
    from math import floor
    g=F('.53')-d
    tab={c:P[M]/g for c in range(1,M+1)}
    for i in range(M-1,1,-1):
        old,tab=tab,{}
        for c in range(i,0,-1):
            l=i-c; b=F(l,i)
            charge=b*(tab[c+1]+P[i]) if l else 0
            val=max((r*old[c]+charge)/(g+r+b) for r in (l,l+F(c,2)))
            tab[c]=F(floor(10**9*val),10**9)
    assert tab[1]>1
    print('fixed-edge relaxation at rho=.53 requires virgin >=',tab[1])


if __name__=='__main__': main()
