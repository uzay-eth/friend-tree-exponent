#!/usr/bin/env python3
"""Exact Laurent coefficients for the physical planted-fringe resolvents."""

from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from math import factorial, lcm


# Canonical core fringe: (birth label, real leaf count, sorted core children).
# Label zero represents the unmarked geometry. Every real core has a parent
# slot; the top's slot is silent. It contributes to degree but not size.
def virgin(label=0):
    return (label, 1, ())


def degree(s):
    return 1 + s[1] + len(s[2])


def b(s):
    return Q(s[1], degree(s))


def h(s):
    return Q(1, 2) if s[1:] == (1, ()) else Q(0)


def gain(s):
    return s[0], s[1] + 1, s[2]


def replace(s, index, child):
    children = s[2]
    return s[0], s[1], tuple(sorted(children[:index] + (child,) + children[index + 1 :]))


@cache
def size(s):
    return 1 + s[1] + sum(map(size, s[2]))


@cache
def total_b(s):
    return b(s) + sum(map(total_b, s[2]))


@cache
def intrinsic_jumps(s, incoming=0):
    a, leaves, children = s
    d = degree(s)
    p = Q(1, d)
    rate = leaves + incoming + sum(Q(1, degree(y)) for y in children)
    jumps = [(gain(s), rate)] if rate else []
    if leaves:
        y = virgin(min(d, 20) if a else 0)
        jumps.append(((a, leaves - 1, tuple(sorted(children + (y,)))), b(s)))
    for j, child in enumerate(children):
        for nxt, rate in intrinsic_jumps(child, p):
            jumps.append((replace(s, j, nxt), rate))
    return tuple(jumps)


def row(s, mode):
    p = Q(1, degree(s))
    intrinsic = intrinsic_jumps(s)
    holding = size(s) - p
    assert sum(rate for _, rate in intrinsic) == holding
    yield from intrinsic
    if mode == 'intrinsic':
        yield s, -holding
        return
    if mode == 'signed':
        yield s, -holding
        for child in s[2]:
            yield child, -p
            yield gain(child), p
        return
    yield s, -holding - Q(1, s[0])
    for child in s[2]:
        unshifted = Q(1, child[0]) - p
        assert unshifted >= 0
        if unshifted:
            yield child, unshifted
        yield gain(child), p
    if mode == 'trigger' and degree(s) > 2 and s[1]:
        yield virgin(min(degree(s), 20)), b(s)


def series(mode, count, start, rewards):
    scale = lcm(*range(1, max(20, count + 1) + 1))
    weights = {start: 1}
    result = [[] for _ in rewards]
    def scaled(value):
        assert scale % value.denominator == 0
        return value.numerator * (scale // value.denominator)

    for n in range(count):
        for values, reward in zip(result, rewards):
            values.append(Q(sum(weight * scaled(reward(s)) for s, weight in weights.items()), scale ** (n + 1)))
        if n + 1 == count:
            break
        if n + 2 == count and mode == 'trigger' and rewards == [h]:
            # Intrinsic jumps and shifted copies never land on a virgin.
            def last(s):
                p = Q(1, degree(s))
                return ((b(s) * Q(1, 2) if degree(s) > 2 else Q(0))
                        - (size(s) - p + Q(1, s[0])) * h(s)
                        + sum((Q(1, y[0]) - p) * h(y) for y in s[2]))
            result[0].append(Q(sum(w * scaled(last(s)) for s, w in weights.items()), scale ** (n + 1)))
            break
        next_weights = defaultdict(int)
        for s, weight in weights.items():
            for nxt, rate in row(s, mode):
                next_weights[nxt] += weight * scaled(rate)
        weights = {s: w for s, w in next_weights.items() if w}
        intrinsic_jumps.cache_clear()
        size.cache_clear()
        total_b.cache_clear()
    return result


def vp(q, p):
    assert q
    n, d = q.numerator, q.denominator
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    while d % p == 0:
        d //= p
        v -= 1
    return v


def multiply(a, b):
    return [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(len(a))]


def main():
    u, = series('signed', 14, virgin(), [b])
    root, total = series('intrinsic', 14, virgin(), [b, total_b])
    # U=G/(1+G), with G=V(Tb). Compare all 14 signed coefficients.
    assert [x + y for x, y in zip(u, multiply([Q(0)] + u, [Q(0)] + total)[1:])] == total
    print('Signed generator and independent total-promotion transfer: 14 exact coefficients agree')
    for k, value in enumerate(u, 1):
        print(f'u_{k} = {value}')
    a, = series('trigger', 13, virgin(2), [h])
    for k, value in enumerate(a, 1):
        print(f'a_{k} = {value}')
    for p in (5, 7, 11, 13):
        star = Q(factorial(p - 2) * (p - 1), p)
        x, y = u[p - 2], a[p - 1]
        assert vp(x, p) == vp(y, p) == -1
        assert (root[p - 2] - star).denominator % p
        assert (x - root[p - 2]).denominator % p
        assert (y - star / 2).denominator % p
        print(f'p={p}: v_p(u_{p-1})=v_p(a_{p})=-1; both star-subtracted remainders p-integral')
    for label, values in [('u_1..u_14', u), ('a_1..a_13', a)]:
        primes = (2, 3, 5, 7, 11, 13, 17, 19)
        supports = [tuple(p for p in primes if q.denominator % p == 0) for q in values]
        for q, support in zip(values, supports):
            den = q.denominator
            for p in support:
                den //= p ** (-vp(q, p))
            assert den == 1
        print(f'{label} denominator prime supports: {supports}')


if __name__ == '__main__':
    main()
