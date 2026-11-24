"""Independent rational audit; no imports from the writer's checker."""

from collections import defaultdict
from fractions import Fraction as F
from functools import cache


def add(g, t):
    return tuple(ns + (len(g),) if v == t else ns for v, ns in enumerate(g)) + ((t,),)


EDGE = ((1,), (0,))
T3 = add(EDGE, 0)


@cache
def arrows(g):
    rates = defaultdict(F)
    for ns in g:
        for v in ns:
            rates[v] += F(1, len(ns))
    assert sum(rates.values()) == len(g)
    return tuple(rates[v] / len(g) for v in range(len(g)))


@cache
def stay(j, m, a, B):
    # Solve the Dirichlet problem by back substitution on differences,
    # without computing the upward-first-passage ratios.
    f = lambda i: min(F(1), F(2, i + 1))
    if j >= m:
        return 1 - f(j)
    if not a and not j:
        return F(0)
    low = 0 if a else 1
    coeffs = [F(1)]
    for i in range(low, m):
        p = coeffs[-1]
        p0 = coeffs[-2] if i > low else F(0)
        D = max(2, i + B)
        coeffs.append(p + ((i + 1) * p - i * p0) / (D * (i + a)))
    v = [c * (1 - f(m)) / coeffs[-1] for c in coeffs]
    assert v == sorted(v)
    for i in range(low, m):
        k = i - low
        prev = v[k - 1] if k else 0
        for d in (max(2, i + B), max(2, i + B) + 1, max(2, i + B) + 9):
            drift = (i + a) * (v[k + 1] - v[k]) + (i * (prev - v[k]) - v[k]) / d
            assert drift >= 0
            if d == max(2, i + B):
                assert drift == 0
    return max(v[j - low], 1 - f(j))


@cache
def core(d, l, cap=12):
    # Forward mass propagation to terminal survival outcomes.
    states = {l: F(1)}
    while d < cap:
        next_states = defaultdict(F)
        for stock in range(max(states, default=0), 0, -1):
            p = states.get(stock, 0) / (d + 1)
            if p:
                states[stock - 1] = states.get(stock - 1, 0) + p
                next_states[stock + 1] += d * p
        d += 1
        states = next_states
    return sum((p * (1 - F(1, d)**k) for k, p in states.items()), F(0))


@cache
def hub(g, v):
    n = tuple(map(len, g))
    d = n[v]
    stocks = lambda w: sum(n[y] == 1 for y in g[w])
    if d == 1:
        (u,) = g[v]
        j = stocks(u) - 1
        return F(0), 1 - stay(j, 32, 1, n[u] - j)
    l = stocks(v)
    if l:
        return core(d, l), F(1)
    return F(0), min(F(1), sum(1 - stay(stocks(u), 32, 0, n[u] - stocks(u)) for u in g[v]))


def rooted(g, v, p=-1):
    return tuple(sorted(rooted(g, w, v) for w in g[v] if w != p))


@cache
def diff(g, u, v, h):
    if rooted(g, u) == rooted(g, v):
        return F(0), F(0)
    a, b = hub(g, u)
    c, d = hub(g, v)
    lo, hi = a - d, b - c
    if h:
        children = [diff(add(g, t), u, v, h - 1) for t in range(len(g))]
        lo = max(lo, sum(p * a for p, (a, _) in zip(arrows(g), children)))
        hi = min(hi, sum(p * b for p, (_, b) in zip(arrows(g), children)))
    assert lo <= hi
    return lo, hi


def adjacent(g, h):
    # Enumerate two genuine jumps. Pair only after obtaining their full
    # normalized weights, without using the old-edge formula to compute them.
    n = len(g)
    probs = {(i, j): p * q for i, p in enumerate(arrows(g))
             for j, q in enumerate(arrows(add(g, i)))}
    assert sum(probs.values()) == 1
    lo = hi = F(0)
    for (i, j), p in probs.items():
        if j < i or j == i:
            continue
        c = p - probs.get((j, i), 0)
        r = [n * x for x in arrows(g)]
        if j == n:
            expect = r[i] / (len(g[i]) + 1)
        elif j in g[i]:
            expect = r[j] / (len(g[j]) * (len(g[j]) + 1)) - r[i] / (len(g[i]) * (len(g[i]) + 1))
        else:
            expect = F(0)
        assert c * n * (n + 1) == expect
        if c:
            bounds = diff(add(add(g, i), j), n, n + 1, h)
            a, b = sorted(c * x for x in bounds)
            lo += a
            hi += b
    return lo, hi


def exact_gaps(h):
    ctr = diff(T3, 0, 1, h)
    return ((ctr[0] / 2, ctr[1] / 2), adjacent(EDGE, h), adjacent(T3, h))


def newborn_classes(g):
    """Coefficients by the final shape and endpoint orbit, before bounding."""
    n = len(g)
    out = defaultdict(F)
    for i, p in enumerate(arrows(g)):
        a = add(g, i)
        for j, q in enumerate(arrows(a)):
            b = add(a, j)
            x, y = rooted(b, n), rooted(b, n + 1)
            if x < y:
                out[x, y] += p * q
            elif x > y:
                out[y, x] -= p * q
    return sorted(abs(c) for c in out.values() if c)


def individual(n):
    law = {EDGE: F(1)}
    for _ in range(2, n):
        law = {add(g, t): w*p for g, w in law.items() for t, p in enumerate(arrows(g))}
    assert sum(law.values()) == 1
    for v in range(5):
        bounds = [(w, hub(g, v)) for g, w in law.items()]
        yield tuple(sum(w * ab[k] for w, ab in bounds) for k in (0, 1))


def report(x, upper):
    unit = 10**12
    i = -((-x * unit) // 1) if upper else (x * unit) // 1
    return f'{"-" if i < 0 else ""}{abs(i)//unit}.{abs(i)%unit:012d}'


if __name__ == '__main__':
    assert stay(1, 4, 1, 1) == F(108, 305)
    assert (stay(1, 4, 1, 1) - F(1, 4)) / 2 == F(127, 2440)
    assert core(2, 2, 4) == F(437, 512)
    assert stay(0, 4, 0, 1) == 0 < stay(0, 4, 1, 1)
    assert diff(EDGE, 0, 1, 0) == diff(T3, 1, 2, 0) == (0, 0)
    path = add(add(T3, 1), 3)
    assert set(path[2]) != set(path[4]) and diff(path, 2, 4, 0) == (0, 0)
    assert newborn_classes(EDGE) == [F(1, 6)]
    assert newborn_classes(T3) == [F(1, 72), F(1, 24), F(1, 18)]
    print('exact small certificates: 108/305, 127/2440, 437/512; swap controls: edge, cherry, reflected P5')
    print('nonzero orbit coefficient magnitudes: T2 [1/6]; T3 [1/72, 1/24, 1/18]')
    points = list(individual(7))
    assert points[0] == points[1]
    for i, (lo, hi) in enumerate(points, 1):
        print(f'p{i}, size 7 ultimate bounds: [{report(lo, False)}, {report(hi, True)}]')
    for h in (0, 4):
        for i, (lo, hi) in enumerate(exact_gaps(h), 2):
            assert hi > 0 and ((lo > 0) if h or i == 2 else (lo < 0))
            print(f'p{i}-p{i+1}, h={h}: [{report(lo, False)}, {report(hi, True)}]')
