"""Cold arithmetic: retry clock, owner budgets, and founder interpolation."""

from fractions import Fraction as Q
from math import factorial


def main():
    slack = Q(1,200000)
    waits = 40+20*31+1
    age = 128
    stock = 2**32
    assert Q(19,20)*age>121
    assert Q(1,2)*age==64
    for families,exponent in ((20,121),(2,64)):
        assert Q(2*families*stock,2**exponent)+Q(8,stock)<Q(1,2**27)
    assert waits+age+1==790

    def eta(p):
        return Q(1,2**27)+790*p

    original = eta(Q(1,10**9))
    assert 6*original<slack and 4*original<slack
    print(f'PASS 661 mean retry clock; original maximum loss {original}; 4η and 6η fit')

    for n in range(1,4):
        # Internal bare-core weights are at most 1/3 and their s at most 1/6.
        owner = Q(n,3)+2*Q(1,3)+(n-1)*Q(1,6)
        assert owner<=2
    extended = eta(Q(3,2*10**9))
    assert 2*extended<4*extended<slack
    print(f'PASS refined chain budget 2η; extended retained slack {slack-4*extended}')

    # Taylor tail from 1/5! is at most a geometric series of ratio 1/6.
    e_upper = Q(11,4)
    assert sum(Q(1,factorial(k)) for k in range(5))+Q(1,100)<e_upper
    a = 1-Q(108,305)
    b = Q(387,1000)  # Published physical mean ceiling, strictly above its exact cap.
    ratio = min((2/e_upper)/a,(3/e_upper**2)/b)
    assert ratio>1
    print(f'PASS founder endpoint ratios; log concavity gives uniform ratio > {ratio}')


if __name__=='__main__':
    main()
