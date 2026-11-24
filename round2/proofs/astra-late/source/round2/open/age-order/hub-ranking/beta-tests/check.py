"""Exact constants, kernel identities, and scope controls for binomial tests."""
from fractions import Fraction as Q
from math import comb, prod
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
MELLIN = run_path(str(HERE.parent/'mellin/check.py'))


def rise(x,j):
    return prod(x+i for i in range(j))


def test(k,r,z):
    return sum((comb(k,j)*z**j*(1-z)**(k-j) for j in range(r,k+1)), Q(0))


def kernel(k,r):
    return [r*comb(k,r)*(-1)**j*comb(k-r,j) for j in range(k-r+1)]


def integral(k,r,z):
    return sum(a*z**(r+j)/(r+j) for j,a in enumerate(kernel(k,r)))


def upper_family():
    c = Q(28,729)
    log_upper = sum(Q((-1)**(j-1), j*2**j) for j in range(1,8))
    # Uniform lower middle constant, using the monotone finite b factors.
    assert Q(2,9)*(1-log_upper)*Q(2,3)**3 >= c
    # For l=k-r+1 and b>=3, (l+b)/b is largest at b=3.
    # r>=400(l+5) implies r-1>=400(l+3).
    assert 4*Q(8,25)/(400*3) < c*c
    assert 4*Q(8,25)*5/(3*(1650-1)) < c*c
    assert 4*Q(8,25)*6/(3*(1800-1)) < c*c
    print('PASS uniform upper-family constants, including one and two misses')


def fixed_hit_family():
    # Positive endpoint comparison, not signed moments subtracted.
    for r in range(3,9):
        p=Q(r)
        middle=MELLIN['middle'](p,p)
        t2,t3=MELLIN['positive_tip'](p,p)
        debit=(4*t2+2*t3)/3
        kmax=10*r
        # The loss uses a separate unsigned Yule upper for the middle.
        target = Q((kmax-r)*r*(r+2))/(middle-debit)
        n = 1000*(target//1000+1)
        loss=Q((kmax-r)*r*(r+2),n+r+2)
        assert middle-loss>debit
        print(f'PASS fixed hit r={r}: {r}<=k<={kmax}, n>={n}, margin {float(middle-loss-debit):.8g}')


def obstruction():
    # [a,b] is a negative tail rectangle; [c,d] is a translated positive
    # rectangle of the same mass. They are genuine differences of law tails.
    a,b,c,d = Q(1,4),Q(3,8),Q(3,8),Q(1,2)
    base = Q(1,800)
    w = Q(1,2**21)
    # Relabelled two-point comparison, padded with a common atom at zero.
    # U mass w at a,d, V mass w at b,c. Add Beta(1,21) with mass base
    # to U, thereby meeting the sharpened all-power floor at n=20.
    # The right translated rectangle dominates the left for every p>=1.
    assert d-c == b-a and c==b and a<Q(1,3)<b
    # This positive tail supplies 420*base times the gamma ratio in (1).
    assert 420*base>Q(1,40)
    # At target its contribution is strictly less than the negative atom.
    assert base*(1-Q(1,3))**21 < w
    assert base+2*w<1 and base>=Q(371,12500)/(20*21)
    # Even the individual birth-degree beta upper caps hold.
    assert 2*w<=(1-b)**21
    assert 2*w+base*(1-a)**21<=(1-a)**20
    assert w+base*(1-d)**21<=(1-d)**20
    # Concrete beta test concentrating around x=1/3.
    k,r=600,201
    gap=test(k,r,d)-test(k,r,c)+test(k,r,a)-test(k,r,b)
    assert gap < -Q(9,10)
    profile=base*rise(Q(k-r+1),21)/rise(Q(k+1),21)
    assert profile+w*gap<0
    print('PASS all-power obstruction at x=1/3 and (k,r)=(600,201)')


def transfer():
    for k in range(1,13):
        for r in range(1,k+1):
            assert integral(k,r,Q(1))==1
            for z in (Q(1,3),Q(1,2),Q(4,5)):
                assert integral(k,r,z)==test(k,r,z)
    # Concrete one-sided sandwich at 1/3 with honest binomial errors.
    k=90
    x=Q(1,3)
    rl,rh=26,35
    assert test(k,rh,x)<Q(1,5) and 1-test(k,rl,x)<Q(1,5)
    print('PASS 234 kernel/mixture comparisons and rational threshold errors at 1/3')


if __name__=='__main__':
    upper_family()
    fixed_hit_family()
    obstruction()
    transfer()
