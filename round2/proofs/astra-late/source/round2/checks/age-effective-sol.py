"""Exact arithmetic for an explicit nineteen-label spectral margin."""

from fractions import Fraction as F
from pathlib import Path
from runpy import run_path


l, u, rho, A = F(113,200), F(1133,2000), F(51,100), F(20,3)
b = (l-rho)/2


def constants():
    e = min(leaf_path(j) for j in range(2,21))
    g = e*(l-rho)/(19*A)
    a = g/64
    L, L2 = A/(g*b**2), A/(g*b**3)
    r = min(b, a/(8*L), a/(2*(L2+8*L**2/a)), F(1,4))
    eta = min(r/2, a*r*r/(131072*L))
    M = max(2, (8/a+16/(a*eta))/g, 131072/(g*a*r*r))
    h = g*g*(l-rho)**2/(A*(u+2)*(u+3))
    assert rho+b < l and u+2 < 3 and F(1,24)*(1+F(9,50)) <= l-rho
    assert g/(2*(u+2)**2) >= a
    assert 8*L*r <= a and (L2+8*L**2/a)*r <= a/2
    assert eta*L <= a*r*r/131072 and eta <= r/2
    return g, r, eta, M, h


def leaf_path(j):
    p = F(1)
    for k in range(2,j+1):
        p *= F(k-1) / (u+F(1,2)+k-1+F(k-1,k))
    return p/j


def lower_power(x):
    k = 0
    while F(1,2**k) > x:
        k += 1
    return k


def main():
    g,r,eta,M,h = constants()
    q = F(9,25)
    sigma_lo = rho+b
    B1 = A+64*A/b
    H = 2*B1+2
    residue = 64/(g*g)
    Cv = A+2*H*M*(A/b)**2+4*B1**2/H
    CJ = A+A*Cv/b+A*residue/(l-rho)
    lead = A*residue/(l-rho)**2
    CP = 1+CJ/b+lead
    h_hi = lead*(A/l+2/(l-q))
    C = 1+CP*(A/sigma_lo+2/(sigma_lo-q))+h_hi
    # Uniform phase loss from the exponential label-2 submeasure:
    # (1+x)^.5*(1+(1+x)^.5) <= 4 for x in [0,1/9].
    assert F(4,3)*(1+F(4,3)) < 4
    assert F(1,6*4*9) >= F(1,64*4096)
    # Coefficients for the quadratic bound on the core-degree moment.
    assert (40, F(35,4), F(1,4)) == (40,40*F(1,4)-F(40,32),F(40,32)-1)
    # Integer roots must be enclosed exactly; verify the floor and the
    # factor-two loss used in the clock estimate on this finite grid.
    for n in range(8,513):
        m = max(k for k in range(n+1) if k**20 <= n**19)
        assert m**20 <= n**19 < (m+1)**20 and (2*m)**20 >= n**19
    f = Path(__file__).resolve().parents[1]/'open/age-order/effective/threshold.py'
    w = run_path(str(f))
    w['self_check']()
    assert C <= 2**453 and h_hi <= 2**68 and h >= F(1,2**60)
    assert w['BETA']/8 > F(1,2**7)
    for k,D in ((0,2),(3,7),(9,200)):
        J = 96*(k+1)*(D+1)**(k+1)+6
        j = (J-1).bit_length()
        s = max(515*2**348,20*2**348*(982+2*k),25*(71+j))
        rows = w['envelopes'](F(2**453),F(1,2**348),F(1,2**60),F(2**68),k,D)
        # Verify using exact coefficient/budget bit lengths without forming
        # the large powers that the actual cutoff can entail.
        assert s >= w['cutoff_power'](rows)
    for name,x in (('g_*',g),('r',r),('eta',eta),('h_*',h)):
        k = lower_power(x)
        assert F(1,2**k) <= x and (not k or x < F(1,2**(k-1)))
        print(f'{name} >= 2^-{k}')
    print(f'inverse bound M <= 2^{lower_power(1/M)}; 505 clock floor controls')
    print(f'C <= 2^{lower_power(1/C)}; h(T3) <= 2^{lower_power(1/h_hi)}')
    print('symbolic cutoff formula: 3 radius/cap controls')


if __name__ == '__main__':
    main()
