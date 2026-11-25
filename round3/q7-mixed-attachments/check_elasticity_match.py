"""Outward bounds for the new weighted stopping inequality at p=927/1000."""
from fractions import Fraction as Q

from review_high_p_bounds import exp_pos, neg_log_1m


def certificate():
    p, a, i, delta, eta, theta = map(Q, ['.927', '.44', '.795', '.25', '.11', '.26'])
    rho = (1 - p) / p
    assert a >= rho / (1 - i) and i * p > (1 - p)**2
    x_lo, x_hi = exp_pos(i / rho - Q(4, 5))
    assert x_lo > 100 and (1 - delta) * x_lo >= 1
    chernoff = eta - (1 - eta) * neg_log_1m(eta)[1]
    b1 = chernoff * delta * x_lo
    atten = rho * (1 - delta) / ((1 - eta) * delta)
    gamma = 1 / exp_pos(atten)[1]
    denom = 2 + (i + Q(17, 10) * rho)**2
    bf = theta**2 * gamma**2 * x_lo * (1 - i) * (1 - delta) / (2 * (1 + a)**2 * denom)
    assert b1 > 17 and bf > Q('5.60')
    assert bf * (1 - a - rho / Q('1.4142')) > a and Q('1.4142')**2 < 2
    main = exp_pos(atten + (1 + a) * neg_log_1m(delta)[1])[1] / (1 - theta)
    tail1 = Q('1.01') * x_hi / ((1 - i) * exp_pos(b1)[0])
    tail2 = Q('1.01') / ((1 - i) * (1 - eta) * delta * exp_pos(bf)[0])
    upper = max((1 - a) / (1 - i), main + tail1 + tail2)
    pi_hi = Q(22, 7)
    sin_lo = 1 - ((Q(1, 2) - a) * pi_hi)**2 / 2
    required = upper * pi_hi * (1 + a) / sin_lo
    assert upper < Q('2.751') and 0 < tail1 and 0 < tail2
    assert required < Q('12.68') < 1 / rho
    print('p=927/1000: B1 > 17, B2*f0 > 5.60, C < 2.751; both tail bounds positive')
    print('weighted stopping cost < 12.68 < 927/73; uniform upper interval starts at 927/1000')


if __name__ == '__main__':
    certificate()
