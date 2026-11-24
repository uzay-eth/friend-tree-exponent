"""Rational rectangle certificate; no simulation or physical-tree sweep."""
from fractions import Fraction as Q
from math import factorial

SCALE = 10**8


def floor(q):
    return q.numerator // q.denominator


def ceil(q):
    return -floor(-q)


def log_upper(h):
    r = h/(2+h)
    return 2*(sum(r**(2*j+1)/(2*j+1) for j in range(8))
              + r**17/(17*(1-r*r)))


def coefficient(B, nq=1024, d=128):
    step = Q(1, d)
    decay = floor(SCALE*(1-step+step**2/2-step**3/6))
    lost = floor(SCALE*(step-step**2/2))
    weights = []
    exp_lo = SCALE
    for _ in range(12*d):
        weights.append(exp_lo*lost//SCALE)
        exp_lo = exp_lo*decay//SCALE
    charges = []
    for j in range(1, nq):
        # Uniform q-bin ((j-1)/nq,j/nq] and Pareto z=(1-q)/q.
        f = (B-1)*Q(nq-j,j)-1
        if f <= 0:
            break
        h = B/f
        if h >= 12:
            break
        charges.append(ceil(SCALE*(h+log_upper(h))))
    subtotal = 0
    for i, weight in enumerate(weights):
        if i:
            subtotal += weight*sum(max(0,SCALE-a*d//i-1) for a in charges)
    return Q(subtotal, nq*SCALE**2)


def corner(B, x, c, root):
    assert (1-x)/(B*x) <= root**6
    theta = Q(3,4)*root
    assert c*c*Q(157,50)*(B+Q(1,2))*x > (Q(4,3)+Q(2,3)*theta)**2*(1-x)


def first_sender_law():
    for b in (2, 20, 21, 30):
        survival = Q(1)
        for k in range(30):
            assert survival == Q(b-1,b+k-1)
            survival *= 1-Q(1,b+k)


def relaxation_ceiling(B, nq=256, d=32):
    """Upper rectangles for this particular hazard/Markov certificate."""
    charges = [0]
    for j in range(1, nq):
        # The actual F is no larger than the unrounded Pareto variable.
        h = B*Q(j,(B-1)*(nq-j))
        charges.append(floor(SCALE*(h+2*h/(2+h))))
    subtotal = 0
    prev = Q(1)
    for i in range(1, 12*d+1):
        survival = (1+Q(i,d*B))**(-B)
        weight = ceil(SCALE*(prev-survival))
        subtotal += weight*sum(max(0,SCALE-a*d//i) for a in charges)
        prev = survival
    return Q(subtotal,nq*SCALE**2)+prev


if __name__ == '__main__':
    first_sender_law()
    gamma_third = sum((-1)**k*Q(3,3*k+1)/factorial(k) for k in range(6))
    assert Q(22,7)<Q(1773,1000)**2
    assert Q(1773,1000)/gamma_third*Q(42,41)<Q(3,4)
    for B, x, c, mesh in [(20,Q(151,250),Q(1635,10000),(1024,128)),
                           (21,Q(591,1000),Q(1639,10000),(1024,128)),
                           (30,Q(1,2),Q(167,1000),(2048,256))]:
        assert coefficient(B,*mesh)>c
        corner(B,x,c,Q(57,100))
        print(f'PASS b>={B}, x>={x}: middle coefficient > {c}')
    assert relaxation_ceiling(21)<Q(18,100)
    # The actual normalized all-tree debit used by this comparison at the
    # corner exceeds .19 even using root >= 3/5 and pi <= 22/7.
    assert Q(3,5)**6<Q(1,21)
    assert Q(19,100)**2*Q(22,7)*Q(43,2) < (Q(4,3)+Q(1,2)*Q(3,5))**2
    print('PASS the chosen scalar certificate cannot cover b=21, x=1/2')
