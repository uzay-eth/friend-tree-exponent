"""Independent exact controls for physical placements and power-order scope."""
from collections import Counter
from fractions import Fraction as Q
from math import factorial, prod


def attach(t, u):
    return tuple(sorted(t + ((u, len(t) + 1),)))


def rates(t):
    n = len(t) + 1
    neighbors = [[] for _ in range(n)]
    for u, v in t:
        neighbors[u].append(v)
        neighbors[v].append(u)
    d = list(map(len, neighbors))
    return d, [sum((Q(1, d[v]) for v in a), Q(0)) for a in neighbors]


def endpoint(t, tag):
    # Erase labels, preserving the marked vertex and its entire future state.
    adj = [[] for _ in range(len(t) + 1)]
    for u, v in t:
        adj[u].append(v)
        adj[v].append(u)

    def root(u, parent):
        return tuple(sorted(root(v, u) for v in adj[u] if v != parent))

    return root(tag, -1)


def identity(t):
    n = len(t) + 1
    d, r = rates(t)
    w = [r[u] / (d[u] + 1) for u in range(n)]
    s = [w[u] / d[u] for u in range(n)]
    direct, reduced = Counter(), Counter()

    def add(dst, done, first, second, mass):
        dst[endpoint(done, first)] += mass
        dst[endpoint(done, second)] -= mass

    for u in range(n):
        once = attach(t, u)
        for v, rate in enumerate(rates(once)[1]):
            add(direct, attach(once, v), n, n + 1, r[u] * rate)
        add(reduced, attach(once, n), n, n + 1, w[u])
    for u, v in t:
        add(reduced, attach(attach(t, u), v), n, n + 1, s[v] - s[u])
    assert direct == reduced
    a = sum(w)
    assert a >= 1 and sum(abs(s[u] - s[v]) for u, v in t) <= a
    x = sum(max(0, s[v] - s[u]) if d[u] == 1 else 0 for u, v in t)
    x += sum(max(0, s[u] - s[v]) if d[v] == 1 else 0 for u, v in t)
    assert 3*x <= a


def power_controls():
    three = ((0, 1), (1, 2))
    d, r = rates(three)
    assert r == [Q(1, 2), Q(2), Q(1, 2)] and d == [1, 2, 1]
    assert endpoint(three, 0) == endpoint(three, 2) != endpoint(three, 1)
    for p in (Q(1), Q(129,128), Q(2), Q(3), Q(7), Q(21,2)):
        for j in range(2, 22):
            phi = lambda k: prod(1 + p/i for i in range(1, k-1)) if k >= 2 else Q(0)
            drift = j*(phi(j+1)-phi(j)) + phi(j-1)-phi(j)-p*phi(j)
            assert phi(2) == 1 and drift >= 0
    for p in range(1, 15):
        beta = lambda a: prod((a+k)/(k+1) for k in range(p))
        assert beta(Q(1,2)) <= Q(1,2)
        for k in range(2, 12):
            a = Q(1,k)
            assert beta(1+a)/(1+k*p) == beta(a)
        normalization = Q(factorial(p), prod(3+k for k in range(p)))
        assert normalization/4 == Q(1, 2*(p+1)*(p+2))
        gaps = [Q(factorial(p), 40*n*(n+1)*prod(n+2+k for k in range(p))) for n in range(2,12)]
        def t(n):
            return Q(factorial(p), 40*(p+1)*prod(n+k for k in range(p+1)))
        assert sum(gaps) == t(2)-t(12)


def scope_controls():
    # Equal mass at 0 and 1 beats a sure 1/4 in powers but loses at zero.
    for p in range(1, 20):
        assert Q(1,2) > Q(1,4)**p
    assert Q(1,2) < 1  # P(older > 0) < P(younger > 0)
    # Equal mass at 1/8 and 1 beats a sure 1/2 in powers and ties at zero,
    # but loses at threshold 1/4 (for all real p>=1, (1/2)^p<=1/2).
    for p in range(1, 20):
        assert Q(1,2)+Q(1,2)*Q(1,8)**p > Q(1,2)**p
    assert Q(1,2) < 1  # older vs younger at threshold 1/4
    # Strengthen the example: the younger law has mass 1/4 at zero and
    # 3/4 at 1/2. Hubhood is now strictly ordered, but the tails still cross.
    for p in range(1, 20):
        assert Q(1,2)+Q(1,2)*Q(1,8)**p > Q(3,4)*Q(1,2)**p
    assert 1 > Q(3,4) > Q(1,2)

    def moment(law, p):
        return sum(w * x**p for x, w in law)

    old = ((Q(1,8), Q(1,2)), (Q(1), Q(1,2)))
    new = ((Q(0), Q(1,4)), (Q(1,2), Q(3,4)))
    mixture = {1: Q(7,11), 2: Q(2,9), 7: Q(13,12)}
    f = lambda x: sum(w * x**p for p, w in mixture.items())
    gap = sum(w*f(x) for x, w in old)-sum(w*f(x) for x, w in new)
    assert gap == sum(w*(moment(old,p)-moment(new,p)) for p,w in mixture.items()) > 0


if __name__ == '__main__':
    trees = [((0, 1),)]
    count = 0
    for _ in range(5):
        for tree in trees:
            identity(tree)
            count += 1
        trees = [attach(tree, u) for tree in trees for u in range(len(tree)+1)]
    print(f'PASS {count} physical two-birth identities on full marked-shape coefficients')
    power_controls()
    scope_controls()
    print('PASS founder stock generator, affine caps, normalized telescope and order scope')
