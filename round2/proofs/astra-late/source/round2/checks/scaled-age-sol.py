"""Exact controls for the scaled adjacent-age audit (standard library only)."""

from collections import defaultdict
from fractions import Fraction as Q


def trees(n):
    if n == 2:
        yield (frozenset({1}), frozenset({0}))
        return
    for old in trees(n - 1):
        for parent in range(n - 1):
            yield attach(old, parent)


def attach(g, parent):
    out = list(g)
    out[parent] = out[parent] | {len(g)}
    return tuple(out) + (frozenset({parent}),)


def counts(g):
    d = tuple(map(len, g))
    return d, tuple(sum(d[j] >= 2 for j in adj) for adj in g)


def physical(g):
    for i, nbrs in enumerate(g):
        for j in nbrs:
            yield i, j, Q(1, len(nbrs))


def endpoint_controls():
    nchecks = 0
    for n in range(2, 8):
        for g in trees(n):
            d, c = counts(g)
            drift_d, drift_c = [Q(0)] * n, [Q(0)] * n
            for i, j, rate in physical(g):
                h = attach(g, j)
                dd, cc = counts(h)
                for v in range(n):
                    drift_d[v] += rate * (dd[v] - d[v])
                    drift_c[v] += rate * (cc[v] - c[v])
                    if cc[v] > c[v]:
                        # A leaf j becomes core, which only its lone neighbour
                        # v can initiate. The new vertex is itself still a leaf.
                        assert v == i and d[j] == 1 and cc[v] == c[v] + 1
            for v in range(n):
                assert 0 <= d[v] - drift_d[v] <= c[v]
                assert drift_c[v] == Q(d[v] - c[v], d[v]) <= 1
                nchecks += 1
    return nchecks


def cut_controls():
    cuts = proposals = tilts = 0
    for v in range(3, 8):
        for birth in trees(v):
            mark = v - 1
            (u,) = birth[mark]
            padded = attach(birth, u)
            p1 = next(r for i, j, r in physical(birth) if (i, j) == (mark, u)) / v
            p2 = next(r for i, j, r in physical(padded) if (i, j) == (u, mark)) / (v + 1)
            assert p1 * p2 >= Q(1, v * (v + 1)**2)
            g = attach(padded, mark)
            mark_side = {mark, v + 1}
            exterior = set(range(v + 2)) - mark_side
            assert padded[v] == {u} and g[v] == {u}
            assert len(mark_side) == 2 and len(exterior) == v
            ds = tuple(map(len, g))
            for side, root, other in ((mark_side, mark, u), (exterior, u, mark)):
                expected = defaultdict(Q)
                for i, j, r in physical(g):
                    if i in side and j in side:
                        expected[i, j] += r
                    elif i in side:
                        assert (i, j) == (root, other) and r == Q(1, ds[root])
                n = len(side)
                assert sum(expected.values()) == n - Q(1, ds[root])
                # Independently enumerate slot/choice products, including
                # proposals that grow only the upper process.
                for top in (n, n + 1, n + 9):
                    law = defaultdict(Q)
                    law[None] = Q(top - n, top)
                    for i in side:
                        choices = sorted(g[i] & side)
                        if i == root:
                            choices.append(None)   # inert slot, not a deleted slot
                        assert len(choices) == ds[i]
                        for j in choices:
                            law[None if j is None else (i, j)] += Q(1, top * len(choices))
                    assert sum(law.values()) == 1
                    assert top * law[None] == top - n + Q(1, ds[root])
                    assert {k: top * r for k, r in law.items() if k is not None} == expected
                    proposals += 1
                    for theta, decay in ((Q(1, 3), Q(3, 5)), (Q(2), Q(1, 4))):
                        q = 1 / (1 + theta * decay)
                        h_now = q**top
                        h_jump = q**(top + 1)
                        # Feynman--Kac/Doob cancellation of the time drift
                        # uses every proposal, including dummy or cut rejection.
                        time_derivative = top * (1 - q) * h_now
                        assert time_derivative + top * (h_jump - h_now) == 0
                        rates = {k: top * r * h_jump / h_now for k, r in law.items()}
                        assert sum(rates.values()) == top * q
                        assert rates[None] == q * (top - n + Q(1, ds[root]))
                        assert all(rates[k] == q * r for k, r in expected.items())
                        tilts += 1
            # A genuine no-cross path is reweighted: even its first holding
            # interval includes the two hazards absent from the inert race.
            a, b = Q(1, ds[mark]), Q(1, ds[u])
            internal = Q(v + 2) - a - b
            joint = internal / (internal + a + b)
            separate = (2 - a) / 2 * (v - b) / v
            assert 0 < separate < joint < 1
            cuts += 1
    return cuts, proposals, tilts


def ratio_controls():
    nchecks = 0
    # Use a positive eta well inside the strict feasible ratio condition.
    for x in (Q(1, 10), Q(1, 3), Q(1, 2), Q(7, 10), Q(19, 20)):
        s = 1 - x
        theta = 1 / s - 1
        assert theta == x / (1 - x)
        assert s - 1 + theta * s == 0   # nonlog terms in the minimized cost
        for ep in (s / 100, s / 1000):
            eta = ep / 2
            denom = 1 - x * (1 + eta)
            assert denom > 0
            z = (1 + ep) * x * (s + ep) / denom
            assert z / ((1 + eta) * z + s + ep) > x
            for w in (z, 2 * z, 10 * z):
                assert w / ((1 + eta) * w + s + ep) > x
                nchecks += 1
    return nchecks


if __name__ == '__main__':
    endpoints = endpoint_controls()
    cuts, proposals, tilts = cut_controls()
    ratios = ratio_controls()
    print(f'{endpoints} endpoint drift controls; {cuts} exact two-jump cuts; '
          f'{proposals} proposal kernels; {tilts} tilted kernels; {ratios} strict ratio controls')
