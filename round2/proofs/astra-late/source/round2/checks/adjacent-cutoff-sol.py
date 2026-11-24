#!/usr/bin/env python3
"""Independent exact budget checks and a negative actual rooted score."""
from fractions import Fraction as Q
from math import factorial


def main():
    d, t, m = 1 << 12, 1 << 24, 1 << 56
    T, R = d*t, 1 << 40
    assert Q(5, 2)**64 > 2*m and Q(5, 2)**d > 8
    assert sum(Q(1, factorial(i)) for i in range(4)) == Q(8, 3)  # e > 8/3
    assert Q(7, 8)*Q(5, 8)/16 > Q(1, 32)
    runs = t//(2*d+128)
    assert runs == 2016
    lag = d*d*(Q(31, 32)**runs + Q(2, m+1))
    assert lag < Q(1, 2**28)

    # Include root in the depth-six incident-degree budget, and all
    # depth-R branching and leaf stocks in the killed-size budget.
    sub_exp = 640 + 60 + 4
    assert sub_exp < 800 and T*4*2 <= R-5
    assert R-5 > sub_exp+30
    E = 8*(60*R+81)
    a = E+60*(R+1)+10
    b = a+2*T+32
    u = 32+2*(b+36)
    assert b-a-2*T == 32 and 2*(b+36)-u == -32

    # Verify the Bernoulli--Poisson leading bound analytically in the note;
    # here check that every chosen numerical parameter is in its domain.
    assert b+36 < u
    eps = Q(1, 2**30) + Q(2, 2**32)
    beta = Q(113, 200)*(1-Q(1133, 2000))/2
    age = Q(113, 200)-Q(1, 4)
    assert 40*(Q(5, 4)/age)*(2/age) < 32**2
    assert Q(5, 4)**4 > 2
    # Geometric path multiplier, as opposed to the published coarse one.
    pre = sum(2**(10*j) for j in range(6)) + 1
    pref = pre*(Q(512, 2**128)+Q(32, 2**80))
    assert pref < Q(1, 2**22)
    assert 60*R+41+9+30 < 2**57
    assert 60*R+41+5+30 < E//8
    error = 6*(pref + Q(2, 2**30) + Q(1, d) + lag + 2*eps)
    assert beta == Q(97971, 800000)
    assert error < Q('0.001467') and beta-error > Q('0.120996')
    print(f'PASS: {runs} exact late-hit blocks; cutoff/Euler and stationary losses')

    k = R+1
    # J <= 2^j: 96(k+1) <= 2^49 and +6 fits in one further bit.
    j = 50+(k+1)*(E+1)
    assert 96*(k+1) < 2**49 and j < 2**92
    # Spatial prefactor: C*6*32*(H+C)/h_* <= 2^975.
    assert 6*32*(2**68+2**453) <= 2**462
    assert beta/8 > Q(1, 2**7)
    constraints = (515*2**348, 20*2**348*(982+2*k), 25*(71+j))
    assert max(constraints) < 2**394 < 2**396
    print('PASS: centered fixed-size transfer; three symbolic exponent tests, no huge integer expansion')

    # Two ordinary edges of a path, plus two inert slots on its left
    # endpoint. At step length 1/5, the first-jump tie between size four
    # and a degree-three middle hit must use the declared hit flag.
    internal = [Q(1, 3), Q(1, 2), Q(1, 2), Q(1)]
    step = Q(1, 5)
    middle_hit = step * (internal[0]+internal[-1])
    assert middle_hit == Q(4, 15) and 1-step*sum(internal) == Q(8, 15)
    assert middle_hit != step*2  # discarding inert slots changes the law
    print('PASS: exact Euler hit/kill tie and inert-slot control')

    # Old core path u--o--w with one leaf on u, none on o, and J-1 on w.
    # Place tips at o and w. Parent o's degree remains three up to its
    # first gain; the total gain rate there is at most two. The tag's
    # promotion rate 1/3 thus beats the first gain with probability 1/7.
    J = 256
    s_o = (Q(1, 2)+Q(1, J))/6
    s_w = (J-1+Q(1, 2))/(J*(J+1))
    contrast_lower = Q(1, 14)-Q(2, J)
    assert s_o > s_w and contrast_lower > 0
    signed_edge_upper = (s_w-s_o)*contrast_lower
    assert signed_edge_upper < 0
    print(f'PASS: finite tree has actual paired edge upper bound {signed_edge_upper} < 0')


if __name__ == '__main__':
    main()
