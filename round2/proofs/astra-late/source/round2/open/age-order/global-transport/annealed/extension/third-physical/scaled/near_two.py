"""Rational serial lower in the large-label size window near y=2."""
from fractions import Fraction as Q
from math import factorial


def exp_bounds(x, n=256):
    term = v = Q(1)
    for k in range(1, n+1):
        term *= x/k
        v += term
    assert x < n+2
    return v, v + term*x/(n+1)/(1-x/(n+2))


def lower():
    # A=4, lower log(4) by the positive atanh series.
    t = 4*sum(Q(1,3)**(2*i+1)/(2*i+1) for i in range(16))
    J, rate, S = 40, 40, 10**14
    _, e1 = exp_bounds(Q(1))
    _, e8 = exp_bounds(Q(8))
    _, et = exp_bounds(rate*(t+Q(1,10**15)))
    assert Q(1,10**15) > 4*Q(1,3)**33/(33*(1-Q(1,9)))
    G = [sum(Q(8**i,factorial(i)) for i in range(j))/e8 for j in range(J)]
    row = [sum(G[j-k]/factorial(k) for k in range(j))*S//e1 for j in range(J)] + [0]
    p, v = 1/et, row[1]/et
    for i in range(1, 161):
        row = [0] + [((rate*(j+1)-j*(j+1)-j)*row[j]
                     +j*(j+1)*row[j+1]+j*row[j-1])//(rate*(j+1))
                    for j in range(1,J)] + [0]
        p *= rate*t/i
        v += p*row[1]
    eyl, _ = exp_bounds(Q(199,100))
    return v*eyl/S


if __name__ == '__main__':
    assert lower() > Q(462,1000)
    print('PASS limiting serial coefficient > .462 uniformly on 1.99 <= y <= 2')
