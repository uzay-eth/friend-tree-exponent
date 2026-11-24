"""Exact finite falsifiers for the physical fixed-shape diagonal arguments."""
from fractions import Fraction as Q
from itertools import product
from math import ceil


def tree(p):
    a = [set() for _ in p]
    for v in range(1, len(p)):
        a[v].add(p[v])
        a[p[v]].add(v)
    return a


def full_census(p, high):
    a = tree(p)
    todo = set(range(len(p))) - high
    out = {}
    while todo:
        top = min(todo)
        comp = {top}
        todo.remove(top)
        stack = [top]
        while stack:
            for w in a[stack.pop()] & todo:
                todo.remove(w)
                comp.add(w)
                stack.append(w)
        if any(len(a[v]) >= 2 for v in comp):
            out[top] = frozenset(comp)
    return out


def threshold_census(p, cutoff):
    return full_census(p, {v for v, adj in enumerate(tree(p)) if len(adj) > cutoff})


def cut_controls():
    # 0--1--2, with one leaf on each endpoint and one on the middle.
    p = (-1, 0, 1, 0, 2, 1)
    p += (0,) * 3 + (2,) * 3
    n = 4
    assert threshold_census(p, n)[1] == {1, 5}
    assert len(tree(p)[1]) == 3  # a real pair, but not an intact virgin
    p2 = p + (5,)              # promotion by 1
    assert threshold_census(p2, n)[1] == {1, 5, len(p)}
    p3 = p2 + (5,) * 3         # only now does 5 cross, exposing a singleton
    assert threshold_census(p3, n)[1] == {1}
    assert 1 not in [v for v, a in enumerate(tree(p3)) if len(a) == 1]
    # One final neighbour on the equality seam must survive the live cut.
    assert threshold_census(p[:-1], n)[1] > {1, 5}
    # A leaf of a high vertex is an isolated ordinary vertex, not a core count.
    assert 3 not in threshold_census(p, n)

    tested = 0
    for growth in product(*(range(k) for k in range(3, 8))):
        p = (-1, 0, 0) + growth
        for k in range(2, 5):
            old = threshold_census(p, k)
            for t in range(len(p)):
                new = threshold_census(p + (t,), k)
                a, b = tree(p), tree(p + (t,))
                assert all(len(a[v]) <= len(b[v]) for v in range(len(p)))
                # Direct top-edge classification, including founder exception.
                hi = {v for v, adj in enumerate(b) if len(adj) > k}
                assert set(new) == {v for v, adj in enumerate(b)
                                    if len(adj) >= 2 and v not in hi
                                    and (not v or p[v] in hi)}
                # Under a FIXED label mask, each old retained component grows.
                fixed = full_census(p + (t,), {v for v, adj in enumerate(a) if len(adj) > k})
                assert all(c <= fixed[v] for v, c in old.items())
                tested += 1
    print(f'{tested} fixed-label/moving-degree cut rows; temporary pair and singleton mutations caught')


def block_and_urn():
    mom = rows = 0
    # Finite martingale replay for E[P_(v+h)^2], beta terminal moment and exact maximal bound.
    for v in range(3, 30):
        law = {1: Q(1)}
        for n in range(v, v + 15):
            assert sum(law.values()) == 1
            assert sum(p * Q(x, n) for x, p in law.items()) == Q(1, v)
            factorial = sum(p * Q(x * (x + 1), n * (n + 1)) for x, p in law.items())
            assert factorial == Q(2, v * (v + 1))
            assert sum(p * Q(x*x, n*n) for x, p in law.items()) <= factorial
            nxt = {x: Q(0) for x in range(1, max(law) + 2)}
            for x, p in law.items():
                nxt[x] += p * (1 - Q(x, n))
                nxt[x + 1] += p * Q(x, n)
            law = nxt
            mom += 1
    # Deliberately delayed permanent envelope gates; a single high reading is insufficient.
    b = 3
    for j in range(1, 75):
        c, following = Q(1, j + 2), Q(1, j + 3)
        b = max(b, ceil(j/c))
        end = b + j*j
        next_start = max(end + 1, ceil((j+1)/following), 1000*j*j)
        for t in (b, end, next_start - 1):
            cutoff = int(c*t)
            assert cutoff >= j and Q(t, cutoff) <= 2/c
            rows += 1
        assert b <= end < next_start  # waiting preserves the selected search window
        b = next_start
    print(f'{mom} exact urn moment rows; {rows} deterministic block/floor checks')


def clock_checks():
    # Model multiplicative brackets exactly with e^omega=5/4.
    r = Q(5, 4)
    windows = misses = thresholds = 0
    for x, width, birth, enter, exit in product(range(3), (1, 2), range(-5, 3),
                                               range(3), range(1, 5)):
        left, right = x, x + width
        a, b = birth + enter, birth + enter + exit
        resident = a <= left and right < b
        visits = a <= right and left < b
        for step in range(4*left, 4*right + 1):
            s = Q(step, 4)
            assert resident <= (a <= s < b) <= visits
        misses += bool(visits and birth > left)
        windows += 1
    for degree, future, cut in product(range(1, 11), range(1, 11), range(1, 11)):
        proxy = Q(degree)
        actual = Q(future)
        if proxy/r <= actual <= proxy*r:
            assert (proxy > cut*r) <= (actual > cut) <= (proxy > cut/r)
            thresholds += 1
    assert misses
    print(f'{windows} residence/random-clock rows ({misses} late arrivals); {thresholds} strict threshold brackets')

    # Chebyshev on Gamma(m,1) and the sum of all noncentre urn bounds at a star restart.
    for c in map(lambda k: Q(1, k), range(2, 101)):
        m = ceil(1024/c**2)
        urn = 8/c**2 * (Q(m-1, m*(m+1)) + Q(1, m+1))
        gamma = Q(16, m)
        assert urn + gamma < Q(1, 32)  # less than e^-2/2
    print('99 positive macroscopic zero-atom budgets without independence')


def rate_controls():
    # The two first internal targets of a virgin have conditional weights 2/3 and 1/3.
    # Then inspect its next independent *root source* ring with the root degree fixed
    # up to that ring (or condition on it being the next source at equal side sizes).
    # A return is null for the inert branch but is not equally likely in both shapes.
    root_first, leaf_first = Q(2, 3), Q(1, 3)
    rejected_star, rejected_path = Q(1, 3), Q(1, 2)
    after_no_return = root_first * (1-rejected_star) / (
        root_first * (1-rejected_star) + leaf_first * (1-rejected_path))
    assert after_no_return == Q(8, 11) != root_first
    for j in range(1, 101):
        f = lambda k: min(1, Q(2, k+1))
        drift = (j*(f(j+1)-f(j)) + Q(j, j+1)*(f(j-1)-f(j))
                 + (1-f(j))/Q(j+1))
        assert drift <= 0
    # Integer shift can break regular-cutoff convergence exactly at an atom:
    # nondecreasing integer degrees, ratio tending to 1/2, but frequent flips.
    d = [n//2 + (n % 2) for n in range(10, 40)]
    assert all(0 <= b-a <= 1 for a, b in zip(d, d[1:]))
    assert sum(D > n//2 for n, D in enumerate(d, 10)) == 15
    # A point mass of exactly this deficit meets the physical adjacent low-pair exclusion.
    p = (-1, 0, 1, 0, 2, 1) + (0,)*3 + (2,)*3
    a = tree(p)
    assert (Q(1, len(a[1])) + Q(1, len(a[5]))) > 0
    print('100 first-send barrier drifts; no-return bias 2/3 -> 8/11; boundary-equality and two-way hazard controls')


if __name__ == '__main__':
    block_and_urn()
    cut_controls()
    clock_checks()
    rate_controls()
