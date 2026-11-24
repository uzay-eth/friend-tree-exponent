"""Integer upper bounds on a first differential send through one old edge."""

from fractions import Fraction
from math import factorial


SCALE = 10**12
LAST_DEGREE = 255
LAST_GAP = 300
LAST_INFECTED = 24


def ceil_div(num, den):
    return -(-num // den)


def tick_upper(q):
    return ceil_div(q.numerator * SCALE, q.denominator)


def displayed(q, places=12):
    unit = 10**places
    n = ceil_div(q.numerator * unit, q.denominator)
    return f'{n // unit}.{n % unit:0{places}d}'


def tagged_bound(s, d, initial_infected=1):
    a, b = s.numerator, s.denominator
    boundary = tick_upper(1 / (1 + s*d))
    next_gap = [boundary] * (LAST_INFECTED+2)
    for gap in range(LAST_GAP, 0, -1):
        current = [boundary] * (LAST_INFECTED+2)
        scale_time = d * (d+gap)
        new_contact = b * d * gap
        tagged_contact = b * gap
        for infected in range(LAST_INFECTED, initial_infected-1, -1):
            new_neighbour = b * scale_time * (gap+infected)
            denom = a*scale_time + tagged_contact + new_neighbour + new_contact
            num = (tagged_contact*SCALE + new_neighbour*next_gap[infected]
                   + new_contact*current[infected+1])
            current[infected] = min(boundary, ceil_div(num, denom))
        next_gap = current
    return Fraction(next_gap[initial_infected], SCALE)


def check_large_degrees(s, proposed):
    # (3/d)^s / [s(1-s)] bounds (d-1) times the discounted hit.
    x = Fraction(3, LAST_DEGREE+1)
    assert x**s.numerator < (proposed*s*(1-s))**s.denominator


def newborn_depth_error(s, q, depth):
    # Split E[min(N_depth(t), c_root(t))] at a deterministic time.
    # The split can be chosen from ordinary integers with no numerical search.
    t = Fraction(2*depth, 5)
    assert s*t < depth+1
    exp_upper = 1/(1+s*t/256)**256
    early = t**(depth+1)/(factorial(depth)*(depth+1-s*t))
    late = t/s + 1/s**2
    return q**depth * exp_upper * (early+late) / (2*s)


def main():
    trials = [('0.39', '0.986'), ('0.5', '0.710'), ('0.565', '0.605')]
    for discount, ceiling in trials:
        s, proposed = Fraction(discount), Fraction(ceiling)
        worst, degree = max(((d-1)*tagged_bound(s, d), d)
                            for d in range(2, LAST_DEGREE+1))
        check_large_degrees(s, proposed)
        assert worst < proposed
        print(f's={s} small-degree max <= {displayed(worst)} at d={degree}; '
              f'all-degree bound < {proposed}')
    s, q = Fraction('0.565'), Fraction('0.605')
    for depth in (16, 32, 40):
        print(f'virgin intrinsic Vb depth {depth}: '
              f'absolute cut error <= {displayed(newborn_depth_error(s,q,depth), 16)}')


if __name__ == '__main__':
    main()
