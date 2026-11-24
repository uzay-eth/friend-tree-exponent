"""Bounded independent first-sender quadrature; no physical-tree enumeration."""

from fractions import Fraction as Q
from math import factorial


S = 10**10


def down(x):
    return x.numerator//x.denominator


def up(x):
    return -down(-x)


def integral_lower(b, bins, rate):
    # Using the exact first-sender floor at b is uniform in later labels:
    # floor((b'-1)z)/b' >= floor((b-1)z)/b need not hold, so use the
    # continuous lower floor bound here, then round the log conservatively.
    score = []
    for j in range(1, bins+1):
        g = Q((b-1)*(bins-j), j)-1
        if g <= 0:
            break
        h = b/g
        if h >= 12:
            break
        r = h/(2+h)
        log_hi = 2*(sum(r**k/k for k in range(1,16,2)) + r**17/(17*(1-r*r)))
        score.append(up(S*(h+log_hi)))
    z = Q(1, rate)
    # Odd Taylor bound for per-bin survival factor. The bin mass is the
    # difference factor times left survival, each bounded separately.
    factor = down(S*sum((-z)**k/factorial(k) for k in range(6)))
    mass = down(S*sum((-1)**(k+1)*z**k/factorial(k) for k in range(1,5)))
    w = S
    integral = 0
    for i in range(12*rate):
        if i:
            integral += w*mass//S * sum(max(0,S-(a*rate+i-1)//i) for a in score)
        w = w*factor//S
    return Q(integral, bins*S*S)


def ceiling(b, bins=256, rate=32):
    # For this nonuniform obstruction, use the exact floor at this label.
    score = [0]
    for j in range(1, bins):
        f = (b-1)*(bins-j)//j
        h = Q(b,f) if f else Q(1000)
        score.append(down(S*(h+2*h/(2+h))))
    integral = 0
    last = Q(1)
    for i in range(1, 12*rate+1):
        tail = (1+Q(i,b*rate))**(-b)
        integral += up(S*(last-tail))*sum(max(0,S-a*rate//i) for a in score)
        last = tail
    return Q(integral,bins*S*S)+last


def main():
    gamma = sum((-1)**j*Q(3,3*j+1)/factorial(j) for j in range(6))
    assert gamma == Q(141077,58240)
    assert Q(22,7) < Q(1773,1000)**2
    assert Q(1773,1000)*Q(42,41) < Q(3,4)*gamma
    for b, x, bins, rate in ((20,Q(151,250),1024,128),
                             (21,Q(591,1000),1024,128),
                             (30,Q(1,2),2048,256)):
        # A finer rational outer sixth root, without floating-point tests.
        ratio = (1-x)/(b*x)
        k = next(k for k in range(1,1000) if Q(k,1000)**6 >= ratio)
        lower = integral_lower(b,bins,rate)
        debit2 = (Q(4,3)+Q(k,2000))**2*(1-x)/(Q(157,50)*(b+Q(1,2))*x)
        assert lower**2 > debit2
        # Print only exact outward rounded margins.
        print(f'PASS b={b}, x={x}: c >= {down(10**6*lower)}/10^6; '
              f'squared margin >= {down(10**6*(lower**2-debit2))}/10^6')
    assert ceiling(21) < Q(18,100)
    assert Q(3,5)**6 < Q(1,21)
    assert Q(19,100)**2*Q(22,7)*Q(43,2) < (Q(4,3)+Q(3,10))**2
    print('PASS exact-F ceiling < 0.18 versus chosen debit > 0.19 at b=21, x=1/2')
    for b in (20,21,30):
        survival = Q(1)
        for k in range(80):
            assert survival == Q(b-1,b+k-1)
            survival *= 1-Q(1,b+k)
    assert max(range(2,20))+1 == 20 and max(range(2,21))+1 == 21
    assert list(range(21,30)) == [21,22,23,24,25,26,27,28,29]
    print('PASS first-sender survival and older-label splice ledger')


if __name__ == '__main__':
    main()
