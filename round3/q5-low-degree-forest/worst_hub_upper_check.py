"""Small exact algebra and numerical stress checks for the worst-hub upper proof."""
from fractions import Fraction as F
from math import exp, expm1, factorial, log, log1p


def stock_and_constants():
    # Boundary values solve the finite birth-death Dirichlet equation.
    checked = 0
    for k in range(2, 41):
        increments = [F(1, factorial(j)) for j in range(1, k+1)]
        s = sum(increments)
        h = [F(0)]
        for inc in increments:
            h.append(h[-1] + inc/s)
        assert h[0] == 0 and h[k] == 1
        for i in range(1, k):
            assert i*(h[i+1]-h[i]) + F(i, i+1)*(h[i-1]-h[i]) == 0
            checked += 1
    theta = F(3,4)
    exp_upper = sum(theta**i/factorial(i) for i in range(6))
    exp_upper += theta**6/factorial(6)/(1-theta/7)
    assert exp_upper < F(17,8)
    assert (exp_upper-1)/2 < F(565713, 10**6)
    assert sum(F(1, factorial(i)) for i in range(5)) > F(8,3)
    # For B=3/8, an old-count-only estimate at d=1/2 undershoots
    # the required cost, although the optimized fresh cost exceeds B.
    mu_hi = F(565737, 10**6)
    assert mu_hi/(1+mu_hi) < F(3,8)
    assert 1+F(1,2)+F(1,8) > F(8,5)
    print('exact stock hitting-equation checks:', checked)
    print('exact concentration slack and raw-count negative control: passed')


def kl(x,p):
    if x == 0:
        return float('inf')
    if p >= 1:
        return -log(x)
    return p*log(p/x)+(1-p)*(log1p(-p)-log1p(-x))


def grid_cost(d, nstep):
    raw_min = cap_min = 1.
    clipped = 0
    for j in range(nstep+1):
        x = (1-d)*j/nstep
        p = x+d if j != nstep else 1.
        direct = -expm1(-kl(x,p))
        raw_min = min(raw_min, direct)
        # Compute the capped Chernoff objective directly from the Poisson MGF.
        if j == nstep:
            cost = 1-x
        else:
            tcap = 0.75/(1-p)
            opt = float('inf') if x == 0 else log(p/x)+log1p(-x)-log1p(-p)
            t = min(tcap, opt)
            clipped += opt > tcap
            cost = -x*expm1(t*(1-p)) - (1-x)*expm1(-t*p)
        cap_min = min(cap_min, cost)
    return raw_min, cap_min, clipped


def numerical_stress():
    active = 0
    gap = 0.
    for j in range(1,100):
        a,b,count = grid_cost(j/100, 2048)
        gap = max(gap, abs(a-b))
        active += count
    assert active > 100 and gap < 2e-12
    print(f'numerical envelope grid: {99*2049} cases, {active} clipped, max minimum-cost difference {gap:.3g}')
    for budget in (0.1, 0.375, 0.9):
        lo,hi = 0.,1.
        for _ in range(28):
            mid = (lo+hi)/2
            cost,_,_ = grid_cost(mid, 512)
            if cost < budget:
                lo = mid
            else:
                hi = mid
        print(f'numerical D* on the finite grid at B={budget:g}: {lo:.5f}')


if __name__ == '__main__':
    stock_and_constants()
    numerical_stress()
