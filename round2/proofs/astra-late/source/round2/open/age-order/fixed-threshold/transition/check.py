"""Exact controls for the two seed likelihoods and the signed leading sector."""
from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path


uniform = run_path(str(Path(__file__).resolve().parents[1] / 'uniform' / 'check.py'))
grow = uniform['grow']


def moment(a, m, k):
    return prod((a + i) / (m + i) for i in range(k))


def rates(tree):
    return tuple(sum(Q(1, len(tree[w])) for w in nbrs) for nbrs in tree)


def likelihood(tree, depth):
    # Both seeds have the same locked deficit D_0-R_0=1/2.
    m, a = len(tree), len(tree[0]) - Q(1, 2)
    count = 0
    for word in product((False, True), repeat=depth):
        states = [(tree, Q(1), Q(1), Q(1))]
        R, B = a, m - a
        expected = Q(1)
        for step in word:
            expected *= (R if step else B) / (R + B)
            R += step
            B += not step
            nxt = []
            for t, ref, like, physical in states:
                for new, q, ell, p, red in uniform['next_rows'](t, Q(1, 2)):
                    if red == step:
                        assert ref * q * like * ell == physical * p
                        nxt.append((new, ref * q, like * ell, physical * p))
                        count += 1
            states = nxt
        assert sum(ref for _, ref, _, _ in states) == expected
    return count


def endpoints(tree):
    # A common finite-horizon test after one further physical birth.
    m = len(tree)
    totals = [Q(0), Q(0)]
    for w, r in enumerate(rates(tree)):
        for j, tag in enumerate((m - 2, m - 1)):
            totals[j] += r / m * (len(tree[tag]) + (tag == w) >= 2)
    return totals


def two_births(tree):
    n = len(tree)
    r = rates(tree)
    direct = Q(0)
    for u in range(n):
        once = grow(tree, u)
        for w, second in enumerate(rates(once)):
            p, q = endpoints(grow(once, w))
            direct += r[u] * second * (p - q)

    leading = remainder = Q(0)
    for u, nbrs in enumerate(tree):
        d = len(nbrs)
        once = grow(tree, u)
        chain = grow(once, n)
        middle, tip = endpoints(chain)
        charge = r[u] / (d + 1)
        if d == 1:
            assert rates(once)[n] == Q(1, 2) and rates(chain)[n] == Q(3, 2)
            leading += charge * middle
            remainder -= charge * tip
        else:
            assert rates(chain)[n] == 1 + Q(1, d + 1) <= Q(4, 3)
            remainder += charge * (middle - tip)
        assert rates(chain)[n + 1] <= Q(1, 2)
        for w in nbrs:
            separate = grow(once, w)
            assert max(rates(separate)[n:]) <= Q(1, 2)
            p, q = endpoints(separate)
            remainder -= r[u] / (d * (d + 1)) * (p - q)
    assert direct == leading + remainder
    assert leading > 0
    return remainder != 0


if __name__ == '__main__':
    leaf = ((1,), (0, 2), (1, 3), (2,))
    middle = grow(leaf, 0)
    rows = sum(likelihood(tree, 3) for tree in (leaf, middle))
    print(f'PASS {rows} exact likelihood prefix checks and all depth-3 words for both seeds')

    count = 0
    for v in (4, 5, 19, 1001):
        for k in range(9):
            assert moment(Q(3, 2), v + 1, k) == 2 * v * moment(Q(1, 2), v, k + 1)
            count += 1
    print(f'PASS {count} exact beta tilt moment identities')

    trees = [((1,), (0,))]
    count = nonzero = 0
    for n in range(3, 6):
        trees = [grow(tree, u) for tree in trees for u in range(n - 1)]
        for tree in trees:
            nonzero += two_births(tree)
            count += 1
    assert nonzero > 0
    print(f'PASS {count} exact signed two-birth rows, {nonzero} with nonzero remainder')
