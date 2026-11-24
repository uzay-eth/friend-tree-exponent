"""Independent exact labelled replay of the adjacent-birth flux."""

from fractions import Fraction as Q
from functools import cache
from itertools import product


@cache
def neighbours(tree):
    rows = [[] for _ in range(len(tree) + 1)]
    for child, parent in enumerate(tree, 1):
        rows[parent].append(child)
        rows[child].append(parent)
    return tuple(map(tuple, rows))


@cache
def rates(tree):
    rows = neighbours(tree)
    out = [Q(0)] * len(rows)
    # Initiating-arrow implementation, not the two-placement formula.
    for row in rows:
        for target in row:
            out[target] += Q(1, len(row))
    return tuple(out)


@cache
def harmonic(tree, label, mode, bound, stop):
    d = len(neighbours(tree)[label])
    size = len(tree) + 1
    if mode == "hit" and d >= bound:
        return Q(1)
    if size == stop:
        return Q(mode == "frac" and 3 * d > bound * stop)
    return sum(r * harmonic(tree + (u,), label, mode, bound, stop)
               for u, r in enumerate(rates(tree))) / size


def compare(tree, h):
    ds = tuple(map(len, neighbours(tree)))
    rs = rates(tree)
    n = len(ds)
    brute = sum(r * s * (h(tree + (u, v), n) - h(tree + (u, v), n + 1))
                for u, r in enumerate(rs)
                for v, s in enumerate(rates(tree + (u,))))

    def hd(u, v):
        final = tree + (u, v)
        return h(final, n) - h(final, n + 1)

    def owner(u):
        return u if ds[u] > 1 else neighbours(tree)[u][0]

    core_score = [Q(0)] * n
    used_edge = used_chain = False
    for u, row in enumerate(neighbours(tree)):
        # For each oriented edge there is exactly one negative delta summand.
        coeff = rs[u] / (ds[u] * (ds[u] + 1))
        for v in row:
            term = -coeff * hd(u, v)
            core_score[owner(u)] += term
            used_edge |= term != 0
        term = rs[u] / (ds[u] + 1) * hd(u, n)
        core_score[owner(u)] += term
        used_chain |= term != 0
    assert sum(core_score) == brute
    assert all(abs(x) <= 3 for x in core_score)
    return used_edge, used_chain


def main():
    counts = [0, 0, 0]
    for n in range(3, 6):
        for tree in product(*(range(v) for v in range(1, n))):
            for mode, bound in (("hit", 2), ("hit", 3), ("hit", 4), ("frac", 1)):
                h = lambda g, x: harmonic(g, x, mode, bound, n + 4)
                edge, chain = compare(tree, h)
                counts[0] += 1
                counts[1] += edge
                counts[2] += chain
            harmonic.cache_clear()
            neighbours.cache_clear()
            rates.cache_clear()
    print(f"{counts[0]:,} exact state/endpoint comparisons; nonzero directed-edge terms in "
          f"{counts[1]:,}, nonzero chain terms in {counts[2]:,}")
    assert counts[1] and counts[2]


if __name__ == "__main__":
    main()
