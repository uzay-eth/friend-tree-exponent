"""Independent alternating-relation and exact coupling checks (no Fable imports)."""

from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from itertools import product


@cache
def tree(par):
    adj = [set() for _ in range(len(par) + 1)]
    for i, p in enumerate(par, 1):
        adj[i].add(p)
        adj[p].add(i)
    return tuple(map(frozenset, adj))


@cache
def witness(A, B, a, b, pa=-1, pb=-1, even=True):
    ca, cb = sorted(A[a] - {pa}), sorted(B[b] - {pb})
    small, large = (cb, ca) if even else (ca, cb)

    def assign(i, remaining):
        if i == len(small):
            return ()
        for q in remaining:
            u, v = (q, small[i]) if even else (small[i], q)
            left = witness(A, B, u, v, a, b, not even)
            if left is not None:
                right = assign(i+1, [w for w in remaining if w != q])
                if right is not None:
                    return left + right
        return None

    kids = assign(0, large)
    return None if kids is None else ((a, b, even),) + kids


def rates(A):
    r = [Q(0) for _ in A]
    for row in A:
        for j in row:
            r[j] += Q(1, len(row))
    assert sum(r) == len(A)
    return r


def coupled(par, other, a, b, M):
    A, B = tree(par), tree(other)
    ra, rb = rates(A), rates(B)
    ca, cb = ra.copy(), rb.copy()
    joint = defaultdict(Q)
    for u, v, even in M:
        assert (ra[u] >= rb[v]) if even else (ra[u] <= rb[v])
        shared = min(ra[u], rb[v])
        joint[u, v] += shared
        ca[u] -= shared
        cb[v] -= shared
    rho = sum(ca)
    assert rho == sum(cb)
    if rho:
        for i, j in product(range(len(A)), repeat=2):
            joint[i, j] += ca[i]*cb[j]/rho
    for i in range(len(A)):
        assert sum(w for (j, _), w in joint.items() if j == i) == ra[i]
        assert sum(w for (_, j), w in joint.items() if j == i) == rb[i]
    for (i, j), mass in joint.items():
        if mass:
            assert witness(tree(par+(i,)), tree(other+(j,)), a, b) is not None
    return sum(w > 0 for w in joint.values()), rho == 0


def main():
    pairs = transitions = zeros = 0
    for n in range(2, 5):
        histories = list(product(*(range(i) for i in range(1, n))))
        for par, other in product(histories, repeat=2):
            for a, b in product(range(n), repeat=2):
                M = witness(tree(par), tree(other), a, b)
                if M is not None:
                    jumps, zero = coupled(par, other, a, b, M)
                    pairs += 1
                    transitions += jumps
                    zeros += zero
        witness.cache_clear()
        tree.cache_clear()
    print(f'{pairs:,} comparable rooted history pairs; {transitions:,} positive joint sites; {zeros:,} zero-residual rows')

    # On the prescribed all-append event after birth on a leaf, the old
    # parent stays at degree 2 and every other neighbour is a leaf.
    for v in range(3, 13):
        initial = (0,) * (v-2) + (1,)
        p = initial
        for n in range(v, v+12):
            assert rates(tree(p))[v-1] == n-v+Q(1, 2)
            p += (v-1,)
    print('120 exact leaf-target append hazards with a degree-two parent')


if __name__ == '__main__':
    main()
