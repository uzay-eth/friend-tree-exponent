"""Outward rational rectangles for the red genealogy comparison."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import runpy

base = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'check.py'))
floor, ceil, log_upper, corner = (base[k] for k in ('floor','ceil','log_upper','corner'))
S = 10**8


def coefficient(B, n=256, d=16):
    step = Q(1,d)
    decay = floor(S*(1-step+step**2/2-step**3/6))
    lost = floor(S*(step-step**2/2))
    weights = []
    e = S
    for _ in range(10*d):
        weights.append(e*lost//S)
        e = e*decay//S

    # Lookup for exp(-r), always rounded up, at r=k/64.
    step = Q(1,64)
    decay = ceil(S*(1-step+step**2/2))
    exp_hi = [S]
    for _ in range(16*64):
        exp_hi.append((exp_hi[-1]*decay+S-1)//S)
    total = 0
    for j in range(1,n):
        f = (B-1)*Q(n-j,j)-1
        if f <= 0:
            break
        C = ceil(S*log_upper(B/(f+Q(1,2))))
        if C >= 10*S:
            continue
        cdf = [0]* (10*d)
        charge = [0]* (10*d)
        for k in range(1,10*d):
            a = Q(k,d)
            r = (f/B)*a/(1+a/B)
            cdf[k] = S-exp_hi[min(floor(64*r),16*64)]
            charge[k] = ceil(C*(1+a/B))
        assert all(a<=b for a,b in zip(cdf,cdf[1:]))
        assert all(a<=b for a,b in zip(charge,charge[1:]))
        for i in range(1,10*d):
            score = 0
            for k in range(1,i):
                p = max(0,S-(charge[k]*d+i-k-1)//(i-k))
                if not p:
                    break
                score += (cdf[k]-cdf[k-1])*p
            total += weights[i]*score
    return Q(total,n*S**3)


def family_controls():
    for f in range(8):
        for s in range(1,16):
            for theta,w in [(Q(1,2),Q(1,3)),(Q(7,8),Q(2,5))]:
                z = theta*(1-w)
                joint = comb(s+f,f)*(1-theta)**(f+1)*theta**s*s*(1-w)**(s-1)
                density = (f+1)*theta*(1-theta)**(f+1)/(1-z)**(f+2)
                conditional = comb(s+f,f+1)*(1-z)**(f+2)*z**(s-1)
                assert joint == density*conditional
                # The atom is exactly the survival of the uncapped law at the cap.
                assert (1+theta/(1-theta))**(-f-1) == (1-theta)**(f+1)
    print('PASS exact red-family conditional weights and cap atom')


if __name__ == '__main__':
    family_controls()
    for b, x, c, root, mesh in [(20,Q(1,2),Q(205,1000),Q(61,100),(512,40)),
                                (21,Q(1,2),Q(2028,10000),Q(61,100),(400,32))]:
        assert coefficient(b,*mesh)>c
        corner(b,x,c,root)
        print(f'PASS b>={b}, x>={x}: refined coefficient > {c}')
