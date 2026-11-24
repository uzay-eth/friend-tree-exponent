"""Exact finite checks for initiating-vertex flux and the reciprocal drift."""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import product


def check_flux():
    rows = 0
    for n in range(3, 9):
        for parents in product(*(range(v) for v in range(1, n))):
            adj = [[] for _ in range(n)]
            for v, p in enumerate(parents, 1):
                adj[v].append(p)
                adj[p].append(v)
            degree = list(map(len, adj))
            direct = {
                (x, degree[x] + 1): sum(Q(1, degree[z]) for z in adj[x])
                for x in range(n)
            }
            split = defaultdict(Q)
            for o in range(n):
                if degree[o] < 2:
                    continue
                for x in adj[o]:
                    split[x, degree[x] + 1] += Q(1, degree[o])
                split[o, degree[o] + 1] += sum(degree[x] == 1 for x in adj[o])
            assert direct == split and sum(split.values()) == n
            rows += 1
    print(f"Initiator split matches exact target-location weights on {rows:,} trees")


def check_drift():
    def f(j):
        return Q(1, j + 2)

    def generator(f, j, immigration):
        return (j + immigration) * (f(j + 1) - f(j)) + Q(j, 2) * (f(j - 1) - f(j))

    for j in range(1001):
        drift = -(j + 1) * generator(f, j, 1)
        assert drift == Q(j * j + j + 2, 2 * (j + 2) * (j + 3)) >= Q(1, 6)
        # Check the two moment equations for a family without immigration.
        assert generator(lambda x: x, j, 0) == Q(j, 2)
        assert generator(lambda x: x * x, j, 0) == j * j + Q(3 * j, 2)
    print("Reciprocal drift and two family moment equations checked at 1,001 stocks")


if __name__ == "__main__":
    check_flux()
    check_drift()
