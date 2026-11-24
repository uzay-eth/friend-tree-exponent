#!/usr/bin/env python3
"""Independent exact controls for the free-tree left-tail audit."""

from collections import Counter
from fractions import Fraction as Q


def arrows(edges):
    degree = Counter(v for edge in edges for v in edge)
    return degree, {(x, y): Q(1, degree[x])
                    for a, b in edges for x, y in ((a, b), (b, a))}


def coupled_event(low, high, source, target, paired):
    _, a = arrows(low)
    _, b = arrows(high)
    common = min(a[source, target], b[source, target])
    if paired:
        assert common > 0
        v = max(max(edge) for edge in high + low) + 1
        low.append((target, v))
        high.append((target, v))
        return v
    surplus = a[source, target] - common
    assert surplus > 0
    low.append((target, max(max(edge) for edge in high + low) + 1))
    return surplus


def delayed_account():
    low = [(0, j) for j in range(1, 65)]
    high = low + [(0, 65)]
    x = coupled_event(low, high, 1, 0, True)
    child = coupled_event(low, high, 0, x, True)
    # The branch promotes and grows before its first mismatch.
    coupled_event(low, high, child, x, True)
    coupled_event(low, high, x, child, True)
    for _ in range(80):
        coupled_event(low, high, child, x, True)
    a, ra = arrows(low)
    b, rb = arrows(high)
    assert a[x] == b[x] == 83 and a[child] == b[child] == 2
    assert {u for u in a if a[u] != b[u]} == {0}
    excess = sum(ra[e] - rb[e] for e in ra if e[0] == 0)
    private_promotion = rb[0, 65]
    assert excess == private_promotion == Q(1, 66) < Q(1, 64)
    assert coupled_event(low, high, 0, x, False) == Q(1, 65 * 66)
    a, _ = arrows(low)
    b, _ = arrows(high)
    assert a[x] == b[x] + 1 and b[x] - 1 == 82
    print('PASS delayed seed account: pre-mismatch degree 83, existing core descendant, gap one on infection')


def rate(edges):
    d, a = arrows(edges)
    return sum(r for (_, y), r in a.items() if d[y] == 1)


def hazards():
    for n in (4, 5, 9, 31, 241):
        path = [(i, i + 1) for i in range(n - 1)]
        star = [(0, i) for i in range(1, n)]
        assert rate(path) == rate(star) == 1
        assert rate(path + [(1, n)]) == Q(7, 6)
        if n >= 5:
            assert rate(path + [(2, n)]) == Q(4, 3)
        # Actual first promotion of a size-n star: one endpoint of
        # the two-core tree has n-2 leaves, the other one has one leaf.
        assert rate(star + [(1, n)]) == Q(n - 2, n - 1) + Q(1, 2)
    # At size K+2 three leaves already imply the gap when C<K.
    # The one-unit looser author gate is harmless; below K+2 a path
    # can have C=K-1 and have rate one.
    K = 13
    assert rate([(i, i + 1) for i in range(K)]) == 1
    n = K + 2
    assert rate([(i, i + 1) for i in range(n - 2)] + [(1, n - 1)]) == Q(7, 6)
    survival = Q(1)
    for j in range(3, 513):
        assert survival == Q(2, j - 1)
        survival *= Q(j - 1, j)
    print('PASS bare paths, sharp 7/6 size-gated example, 5 post-star hazards, exact star survival to size 512')


def budgets():
    # At mu=3/5, take d=a**5 to make the split integral exact.
    for a in (2, 17, 10**12):
        d = a**5
        early = (a**2 - 1) / (Q(2, 5) * d)
        late = Q(5, 3 * a**3)
        assert (early + late) * a**3 == Q(25, 6) - Q(5, 2 * a**2)
    # Use rates from a large star; leaf terms pay a total of one,
    # and the high-degree centre term is d**(1-2*mu) = 1/a.
    a = 5
    d = a**5
    _, rates = arrows([(0, i) for i in range(1, d + 1)])
    full = sum(r * (Q(1, a**6) if y == 0 else 1) for (_, y), r in rates.items())
    assert full == 1 + Q(1, a) < 2
    # Worst possible passage value from a relative overshoot ≤1/4.
    assert 1 / (1 + Q(1, 4)) == Q(4, 5) > Q(3, 4)
    assert 1 - Q(1, 2) / Q(3, 4) == Q(1, 3)
    r = Q(13, 12)
    assert 1 < r < Q(7, 6) < 4 and Q(7, 6) / (Q(7, 6) - r) == 14
    print('PASS rational d^-mu split through degree 10^60; variance, overshoot and hitting-time budgets')


if __name__ == '__main__':
    delayed_account()
    hazards()
    budgets()
