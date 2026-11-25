"""Independent rational certificates for the parameters in the pinned HIGH-P.md."""
from fractions import Fraction as Q

GRID = 1 << 180


def down(x):
    return Q((x * GRID).__floor__(), GRID)


def up(x):
    return -down(-x)


def exp_pos(z):
    """Enclose exp(z) for rational z>=0 by positive Taylor terms and squaring."""
    assert z >= 0
    squares = 0
    while z > Q(1, 2):
        z /= 2
        squares += 1
    total = term = Q(1)
    for j in range(1, 33):
        term *= z / j
        total += term
    remainder = term * z / 33 / (1 - z / 34)
    low, high = down(total), up(total + remainder)
    for _ in range(squares):
        low, high = down(low * low), up(high * high)
    return low, high


def neg_log_1m(z):
    """Enclose -log(1-z) for rational 0<z<1."""
    assert 0 < z < 1
    total = sum(z**j / j for j in range(1, 41))
    return total, total + z**41 / (41 * (1 - z))


def audit(name, p, a, d, eta, theta, i):
    rho = (1-p)/p
    assert rho/(1-i) <= a and i*p > (1-p)**2
    xlow, xhigh = exp_pos(i/rho - Q(4,5))
    assert xlow > 100 and (1-d)*xlow >= 1
    _, eta_log_high = neg_log_1m(eta)
    c_low = eta - (1-eta)*eta_log_high
    assert c_low > 0
    b1lo = c_low*d*xlow
    assert b1lo > 1+a
    attenuation = rho*(1-d)/((1-eta)*d)
    gamma_low = 1/exp_pos(attenuation)[1]
    log_delta_high = neg_log_1m(d)[1]
    main_hi = exp_pos(attenuation + (1+a)*log_delta_high)[1]/(1-theta)
    denominator = 2 + (i + Q(17,10)*rho)**2
    bf_low = theta**2 * gamma_low**2 * xlow*(1-i)*(1-d)/(2*(1+a)**2*denominator)
    root2_low = Q(1414213,10**6)
    assert root2_low**2 < 2
    slope_low = 1-a-rho/root2_low
    assert slope_low > 0 and bf_low*slope_low > a

    # B*exp(-B) is decreasing above one. Cap the *whole product*, since
    # its original prefactor x grows with B. The second tail has fixed prefactor.
    first_arg = min(b1lo, Q(200))
    second_arg = min(bf_low, Q(200))
    tail1 = Q(101,100)*first_arg/((1-i)*c_low*d*exp_pos(first_arg)[0])
    tail2 = Q(101,100)/((1-i)*(1-eta)*d*exp_pos(second_arg)[0])
    c2 = main_hi + tail1 + tail2
    c1 = (1-a)/(1-i)
    required = max(c1,c2)*(1+a)/(a*(1-a))
    assert tail1 > 0 and tail2 > 0 and required < 1/rho

    if name == 'headline':
        assert xlow > Q('21773.76') and xhigh < Q('21773.78')
        assert b1lo > Q('70.76') and bf_low > Q('7.3516') and bf_low*slope_low > Q('4.0495')
        assert main_hi < Q('2.3274') and tail1 < Q('1.62e-26') and tail2 < Q('0.0129')
        assert c2 < Q('2.3403') and c2 < c1 and required == 14
        print('p=187/200: C2 < 2.3403; C=12/5; required=14 < 187/13')
    elif name == 'alternative':
        assert xlow > 40134 and tail1 < Q('1.3e-70') and tail2 < Q('4.2e-13')
        assert c2 < Q('2.223') and required < Q('12.97')
        print('p=19/20: positive T1 < 1.3e-70, positive T2 < 4.2e-13; required < 12.97 < 19')
    else:
        assert b1lo > 10**7 and bf_low > 10**7
        assert tail1 < Q('1e-75') and tail2 < Q('1e-75')
        print(f'p={p}: both negative-exponential arguments have magnitude > 10^7;'
              ' both tail upper bounds strictly between 0 and 10^-75; criterion passes')


def audit_decimal_failure():
    from decimal import Context, Decimal, ROUND_CEILING, Underflow, localcontext

    # This probe models only the documented arithmetic operations of the helper.
    with localcontext(Context(prec=60, rounding=ROUND_CEILING)) as ctx:
        val = Decimal(-10**7).exp()
        assert ctx.flags[Underflow] and val == 0
        ctx.prec = 50
        padded = (+val) + Decimal('1e-45')*abs(+val)
        assert padded == 0
    with localcontext(Context(prec=55, rounding=ROUND_CEILING)):
        pos = Decimal(5)/19
    with localcontext(Context(prec=28)):
        shifted = Q(-pos) + Q(5,19)
    assert shifted == Q(7,19*10**28)
    print('arithmetic probe: padded upper exp(-10^7) = 0;'
          ' default unary negation exceeds -5/19 by 7/(19*10^28)')


if __name__ == '__main__':
    audit('headline',Q(187,200),Q(2,5),Q(6,25),Q(4,25),Q(9,50),Q(3,4))
    audit('alternative',Q(19,20),Q(2,5),Q(1,5),Q(1,5),Q(1,5),Q(3,5))
    for eps,q in [(Q(1,5),Q(1,1000)),(Q(1,10),Q(1,10000)),(Q(1,20),Q(1,100000))]:
        audit('near-one',1-q,(1+eps)*q,eps/16,Q(1,4),eps/16,eps/4)
    audit_decimal_failure()
