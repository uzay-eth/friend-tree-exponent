"""Independent cold checker for response-transcendence-hy4.

Scope: the scalar response as a function of the Laplace parameter s only.
No statement about the arithmetic of the critical root is made or tested.
All tree and series arithmetic here is written from scratch.
"""
from fractions import Fraction as F
from collections import defaultdict
from functools import lru_cache
from math import factorial

A = F(1, 2)


def nb(tag=2):
    return (tag, 1, ())


def plus(t):
    return t[0], t[1] + 1, t[2]


@lru_cache(None)
def size(t):
    return 1 + t[1] + sum(size(c) for c in t[2])


@lru_cache(None)
def deg(t):
    return 1 + t[1] + len(t[2])


def mark(t):
    return F(1, t[0])


def reward(t):
    return F(t[1], deg(t))


def reward_wrong(t):
    return F(t[1], size(t))


@lru_cache(None)
def intrinsic(t, parent=F(0)):
    out = defaultdict(F)
    gain = t[1] + parent + sum(F(1, deg(c)) for c in t[2])
    if gain:
        out[plus(t)] += gain
    if t[1]:
        out[(t[0], t[1] - 1, tuple(sorted((*t[2], nb(2 if deg(t) < 16 else 16)))))] += F(t[1], deg(t))
    for i, y in enumerate(t[2]):
        for yy, r in intrinsic(y, F(1, deg(t))):
            out[(t[0], t[1], tuple(sorted(t[2][:i] + (yy,) + t[2][i + 1:])))] += r
    return tuple(out.items())


def residual(f, t, constant=True):
    aa = A if constant else mark(t)
    p = F(1, deg(t))
    val = -aa * f(t)
    for y, r in intrinsic(t):
        val += r * (f(y) - f(t))
    for y in t[2]:
        ay = A if constant else mark(y)
        val += (ay - p) * f(y) + p * f(plus(y))
    return val


def derivatives(kind, k):
    def g0(t):
        return reward(t)
    fs = [g0]
    for _ in range(k):
        prev = fs[-1]
        if kind == 'intrinsic':
            fs.append(lambda t, p=prev: sum(r * (p(y) - p(t)) for y, r in intrinsic(t)))
        else:
            fs.append(lambda t, p=prev: residual(p, t))
    return fs[k](nb())


def paths(k, rwd=reward):
    cur = {nb(): F(1)}
    for _ in range(k):
        nxt = defaultdict(F)
        for t, w in cur.items():
            for y, r in intrinsic(t):
                nxt[y] += w * r
        cur = nxt
    out = defaultdict(F)
    for t, w in cur.items():
        out[(size(t), deg(t), len(t[2]))] += w * rwd(t)
    return out


def prime_step():
    for p in (3, 5, 7):
        k = p - 2
        star = F(factorial(k) * (p - 1), p)
        by = paths(k)
        hits = [(st, v) for st, v in by.items() if v.denominator % p == 0]
        assert len(hits) == 1 and hits[0][0] == (p, p, 0) and hits[0][1] == star
        assert all(st[1] < p for st in by if st != (p, p, 0))
        for kind in ('intrinsic', 'constant'):
            d = derivatives(kind, k)
            assert d.denominator % p == 0
            assert (d - star).denominator % p != 0
        print('p', p, 'unique terminal star', star, 'and d_{p-2} has p in reduced denominator')
    wrong = paths(3, rwd=reward_wrong)
    assert any(v.denominator % 5 == 0 for st, v in wrong.items() if st != (5, 5, 0))
    print('negative control: degree->size reward breaks uniqueness at p=5')


def series_step():
    # Toy of the formal argument: P(z,H)=H-z-zH^2 = 0, Catalan coefficients.
    H = [F(1)]
    for n in range(1, 10):
        trunc = H[:]
        prod = defaultdict(F)
        for i, a in enumerate(trunc):
            for j, b in enumerate(trunc):
                prod[i + j] += a * b
        pred = prod[n - 1]
        H.append(pred)
    assert H[:7] == [F(1), F(1), F(2), F(5), F(14), F(42), F(132)]
    assert all(h.denominator == 1 for h in H)
    planted = H[:]
    planted[7] += F(1, 7)
    trunc = planted[:7]
    prod = defaultdict(F)
    for i, a in enumerate(trunc):
        for j, b in enumerate(trunc):
            prod[i + j] += a * b
    forced = prod[6]
    assert forced == H[7] and planted[7] != forced
    print('formal recurrence c*d_n = -[z^(n+r)]P(z,H_<n) reproduces Catalan; planted new prime 7 is excluded')


if __name__ == '__main__':
    prime_step()
    series_step()
    print('SCOPE: response as a function of s only; nothing claimed about the critical root')
    print('RESPONSE TRANSCENDENCE COLD CHECK: PASS')
