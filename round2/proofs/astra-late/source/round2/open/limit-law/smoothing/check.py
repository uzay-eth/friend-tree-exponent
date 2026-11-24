"""Exact multi-slot embedded likelihood and Dirichlet concentration checks."""

from fractions import Fraction as Q


def rates(adj, slots):
    degree = [len(row) + s for row, s in zip(adj, slots)]
    internal = sum(Q(1, degree[i]) for i, row in enumerate(adj) for _ in row)
    killed = sum(Q(s, d) for s, d in zip(slots, degree))
    return degree, internal, killed


def main():
    # A central anchor with two cut slots and a single ordinary leaf.
    adj = [[1], [0]]
    weight = Q(1)
    for n in range(2, 34):
        slots = [2] + [0] * (n - 1)
        degree, lam, killed = rates(adj, slots)
        assert lam + killed == n and killed == Q(2, n + 1)
        # Every genuine arrow density acquires the same wait factor.
        for i, row in enumerate(adj):
            for _ in row:
                physical = Q(1, degree[i])
                assert physical / (n * physical / lam) == 1 - killed / n
        weight *= lam / n
        adj[0].append(n)
        adj.append([0])
    # For this prescribed star sequence lambda/n=(n-1)(n+2)/(n(n+1)).
    assert weight == Q(35, 99)
    print("32 multi-slot star states: genuine rates and telescoping terminal product agree")

    for m in range(3, 12):
        for n in (2, 5, 100, 1000):
            second = Q(n * (n + 1), (m * n) * (m * n + 1))
            assert m * second - Q(1, m) == (1 - Q(1, m)) / (m * n + 1)
    print("36 exact symmetric-Dirichlet concentration fixtures")


if __name__ == "__main__":
    main()
