"""Independent cold checker for physical-trigger-hy4 (own code, no author imports).

Small-tree trigger identities, signs, printed drift/correction arithmetic and
deliberate-error controls. Scope: the scalar defined from K_s at a fixed s;
no assertion about the value or arithmetic of mu.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import ceil

S = F('0.565615')


def nb(i):
    return (i, 1, ())


@lru_cache(None)
def degree(t):
    return 1 + t[1] + len(t[2])


@lru_cache(None)
def size(t):
    return 1 + t[1] + sum(size(c) for c in t[2])


def mark(i):
    return F(1, i) if i else F(0)


def b(t):
    return F(t[1], degree(t))


def virgin(t):
    return t[1] == 1 and not t[2]


@lru_cache(None)
def moves(t):
    """Intrinsic plus old-child copy rates; parent slot kept at 1/degree."""
    out = defaultdict(F)
    p = F(1, degree(t))
    gain = t[1] + sum(F(1, degree(c)) for c in t[2])
    if gain:
        out[(t[0], t[1] + 1, t[2])] += gain
    if t[1]:
        out[(t[0], t[1] - 1, tuple(sorted((*t[2], nb(min(degree(t), 20))))))] += F(t[1], degree(t))
    for i, c in enumerate(t[2]):
        out[c] += mark(c[0]) - p
        out[(c[0], c[1] + 1, c[2])] += p
    return tuple(out.items())


def states(cap):
    seen = {nb(i) for i in range(2, cap + 1)}
    todo = list(seen)
    while todo:
        t = todo.pop()
        for u, r in moves(t):
            if r and size(u) <= cap and u not in seen:
                seen.add(u)
                todo.append(u)
    return sorted(seen)


def solve(M, rhs):
    n = len(M)
    a = [list(row) + [v] for row, v in zip(M, rhs)]
    for i in range(n):
        q = next(j for j in range(i, n) if a[j][i])
        a[i], a[q] = a[q], a[i]
        d = a[i][i]
        a[i] = [x / d for x in a[i]]
        for j in range(n):
            if j != i and a[j][i]:
                f = a[j][i]
                a[j] = [x - f * y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def system(s, cap, mode='true', column=0):
    st = states(cap)
    labels = range(2, cap + 1)
    ix = {t: i for i, t in enumerate(st)}
    n = len(st)
    base, killed, trig = [], [], []
    for t in st:
        den = s + mark(t[0]) + size(t) - F(1, degree(t))
        row = [F(0)] * n
        row[ix[t]] = den
        for u, r in moves(t):
            if u in ix:
                row[ix[u]] -= r
        j = min(degree(t), 20)
        assert (b(t) if degree(t) == 2 else 0) == F(virgin(t), 2)
        red = row[:]
        if mode != 'no-nontrigger' and degree(t) > 2 and b(t):
            red[ix[nb(j)]] -= b(t)
        h = F(degree(t) == 2, 2) if mode == 'all-deg2' else F(virgin(t), 2)
        base.append(row)
        killed.append(red)
        trig.append(h)
    A0 = [[row[j] for j in range(n)] for row in base]
    ls = list(labels)
    U = [solve(A0, [F(ls[i] == min(degree(t), 20)) * b(t) for t in st]) for i in range(len(ls))]
    K = [[U[i][ix[nb(j)]] for j in ls] for i in range(len(ls))]
    k = [K[i][column] for i in range(len(ls))]
    C = [[F(i == j) - (K[i][j] if j else F(0)) for j in range(len(k))] for i in range(len(k))]
    v = solve(C, k)
    f = v[0]
    actual = solve([[row[j] for j in range(n)] for row in killed], trig)[ix[nb(2)]]
    return st, ix, K, k, v, f, actual


def identity(s, cap, **kw):
    st, ix, K, k, v, f, actual = system(s, cap, **kw)
    assert f == actual, (s, cap, f, actual)
    for i in range(len(k)):
        lhs = sum(K[i][j] * v[j] for j in range(len(v)))
        assert lhs - v[i] == (f - 1) * k[i]
    return f


def signs(cap):
    vals = [(s, identity(s, cap)) for s in (F('0.4'), F('0.6'), F(1), F(2))]
    assert all(a > c for (_, a), (_, c) in zip(vals, vals[1:]))
    K = system(F('0.6'), cap)[2]
    z = [F(1)] * (len(K) - 1)
    for _ in range(60):
        z = [sum(K[i][j] * z[j] for j in range(len(z))) for i in range(len(z))]
        m = max(z)
        z = [x / m for x in z]
    r = max(sum(K[i][j] * z[j] for j in range(len(z))) / z[i] for i in range(len(z)))
    assert (r > 1) == (dict(vals)[F('0.6')] > 1)
    print('cap', cap, 'F decreasing on grid; sign(F-1) matches sign(r(K)-1), F(.6)=', float(dict(vals)[F('0.6')]))


def negatives(cap):
    true = identity(S, cap)
    for kw in ({'mode': 'all-deg2'}, {'mode': 'no-nontrigger'}, {'column': 1}):
        st, ix, K, k, v, f, actual = system(S, cap, **kw)
        assert f != true or actual != f, kw
    assert F(1, 2 * (S + 2)) < true
    print('cap', cap, 'deliberate errors rejected: all-degree-2 trigger, dropped non-trigger birth, wrong renewal column')


def printed_math():
    for d in range(2, 21):
        for c in range(1, d + 1):
            L = d - c
            gain = F(c, d * (d + 1))
            promo = F(L, d) * (-F(1, d))
            drift = F(L) * gain + promo
            assert drift >= 0 if c >= 2 else drift == -F(L, d * d * (d + 1))
            eta = F(1, d * (d + 1)) if c == 1 else F(0)
            assert (eta == F(1, 420)) == (d == 20 and c == 1)
            assert (eta == 0) == (c >= 2)
    u = F('0.5665')
    g2 = F('0.119101711051')
    assert ceil(10**12 * (1 - g2 / (2 * (u + F(1, 20) + F(3, 2))))) == 971863522077
    s0 = F('0.56564')
    g20 = F('0.117936293836')
    assert ceil(10**12 * (1 + F('0.006') - g20 / (2 * (s0 + F(1, 20) + F(3, 2))))) == 978127513699
    print('c>=2 drift >= 0, c=1 eta=1/420, and both complement contractions verified exactly')


def root_lower(s, uniform):
    """Independent label-20 first-event stencil; unit 1e12, outward floor."""
    D, unit = 20, 10**12
    lam = F(s) + F(1, 2)
    n, q = lam.numerator, lam.denominator
    flat = [[0] * (D + 2) for _ in range(D + 1)]
    for c in range(1, D + 1):
        eta = F(1, D * (D + 1)) if (uniform or c == 1) else F(0)
        flat[D][c] = (F(D - c, D) * unit / (lam + eta)).__floor__()
    for d in range(D - 1, 1, -1):
        for c in range(d, 0, -1):
            L = d - c
            b, g, G = 2 * q * L, 2 * d * q * L, 2 * d * q * L + d * q * (c - 1)
            den, Den = 2 * d * n + b + g, 2 * d * n + b + G
            source = b * flat[d][c + 1]
            v, V = source + g * flat[d + 1][c], source + G * flat[d + 1][c]
            flat[d][c] = min(v // den, V // Den)
    return F(flat[2][1], unit)


def child_message():
    for s, uniform, sharp in (('0.56564', '0.022244698369', '0.022279887682'),
                              ('0.5665', '0.022179480005', '0.022214534253')):
        assert root_lower(s, True) == F(uniform)
        assert root_lower(s, False) == F(sharp)
        assert F(uniform) <= F(sharp)
    print('child message exact: uniform', '0.022244698369/.56564', '0.022179480005/.5665',
          'sharp', '0.022214534253/.5665', 'all reproduced')


if __name__ == '__main__':
    for cap in (3, 4):
        print('cap', cap, 'finite trigger identity exact, F(S)=', float(identity(S, cap)))
        signs(cap)
        negatives(cap)
    printed_math()
    child_message()
    print('SCOPE: scalar defined from K_s at fixed s; no claim about mu')
    print('PHYSICAL TRIGGER COLD CHECK: PASS')
