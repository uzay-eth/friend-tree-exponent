"""Uniform sub-half corner using the exact first-family marginal CDF."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import runpy

base = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'check.py'))
floor, ceil, log_upper, corner = (base[k] for k in ('floor','ceil','log_upper','corner'))
S = 10**9


def coefficient(B, n=448, d=36):
    step = Q(1,d)
    decay = floor(S*(1-step+step**2/2-step**3/6))
    lost = floor(S*(step-step**2/2))
    weights = []
    e = S
    for _ in range(10*d):
        weights.append(e*lost//S)
        e = e*decay//S
    # log(1+y) >= 2y/(2+y); at the q endpoint take the exact integer F_B.
    eps = Q(1,256)
    decay = ceil(S*(1-eps+eps**2/2))
    exp_hi = [S]
    for _ in range(20*256):
        exp_hi.append((exp_hi[-1]*decay+S-1)//S)
    total = 0
    for j in range(1,n):
        f = (B-1)*(n-j)//j
        h = Q(n-j,j)
        g = max(Q(0),B*h-1)
        # For b>=B+1 use the minimum of the initial shifted floor / b
        # and its limiting value. The zero branch decreases to the seam.
        shift = min(h,(B*h+Q(1,2))/(B+1))
        C = max((1+Q(1,B+f+1))*log_upper(Q(2*B,2*f+3)),
                (1+1/(B+2+g))*log_upper(1/shift))
        C = ceil(S*C)
        distr = []
        old = 0
        for k in range(1,10*d):
            a = Q(k,d)
            # (F_b+1)/b >= (F_B+1)/B need not hold across integer jumps.
            # Both (F_b+1)/b > (1-1/b)*(1-q)/q and F_b >= f do hold.
            # For b=B use the exact exponent; for b>=B+1 use the continuous floor.
            v0 = (f+1)*2*a/(2*B+a)
            v1 = Q(B*(n-j), (B+1)*j)*2*a/(2+a/(B+1))
            v = min(v0,v1)
            cdf = S-exp_hi[min(floor(v*256),20*256)]
            distr.append((cdf-old,ceil(C*(1+a/B))))
            old = cdf
        for i in range(1,10*d):
            score = 0
            for k,(mass,charge) in enumerate(distr,1):
                if k>=i:
                    break
                p = S-(charge*d+i-k-1)//(i-k)
                if p<=0:
                    break
                score += mass*p
            total += weights[i]*score
    return Q(total,n*S**3)


if __name__ == '__main__':
    gamma = sum((-1)**k*Q(3,3*k+1)/factorial(k) for k in range(6))
    assert Q(22,7)<Q(1773,1000)**2
    assert Q(1773,1000)/gamma*Q(42,41)<Q(3,4)
    for B,x,c,root,mesh in [(20,Q(49,100),Q(2084,10000),Q(6111,10000),(512,40)),
                             (21,Q(12,25),Q(2087,10000),Q(612,1000),(448,36)),
                             (40,Q(1,3),Q(2088,10000),Q(61,100),(256,24))]:
        assert coefficient(B,*mesh)>c
        corner(B,x,c,root)
        print(f'PASS b>={B}, x>={x}: coefficient > {c}')
