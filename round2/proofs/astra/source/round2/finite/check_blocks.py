"""Exact finite checks for the one-port argument and the free-cut obstruction."""

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement, product
from math import ceil, prod
from typing import NamedTuple


class Side(NamedTuple):
    # A negative parent marks the endpoint of the cut. The port counts in its degree.
    parent: tuple
    leaves: tuple

    def gain(self, vertex):
        leaves = list(self.leaves)
        leaves[vertex] += 1
        return Side(self.parent, tuple(leaves))

    def promote(self, vertex):
        leaves = list(self.leaves)
        assert leaves[vertex] > 0
        leaves[vertex] -= 1
        leaves.append(1)
        return Side(self.parent + (vertex,), tuple(leaves))


@cache
def intrinsic(side):
    neighbors = [[] for _ in side.parent]
    for vertex, parent in enumerate(side.parent[1:], 1):
        neighbors[parent].append(vertex)
        neighbors[vertex].append(parent)
    degrees = [len(items) + side.leaves[i] + (i == 0)
               for i, items in enumerate(neighbors)]
    rates = {}
    for i in range(len(side.parent)):
        rates[i, 'gain'] = side.leaves[i] + sum(F(1, degrees[k]) for k in neighbors[i])
        if side.leaves[i]:
            rates[i, 'bud'] = F(side.leaves[i], degrees[i])
    return rates, F(1, degrees[0])


@cache
def full_rates(sides):
    # Build the complete vertex tree, including individual leaves. Enumerate all
    # ordered source-neighbor choices independently of the one-port formula.
    adj = defaultdict(list)

    def edge(a, b):
        adj[a].append(b)
        adj[b].append(a)

    edge((0, 0), (1, 0))
    for branch, side in enumerate(sides):
        for v, p in enumerate(side.parent):
            if p >= 0:
                edge((branch, p), (branch, v))
            for leaf in range(side.leaves[v]):
                edge((branch, v), (branch, v, leaf))

    rates = defaultdict(F)
    for neighbors in adj.values():
        weight = F(1, len(neighbors))
        for target in neighbors:
            branch, v = target[:2]
            rates[branch, v, 'bud' if len(target) == 3 else 'gain'] += weight
    assert sum(rates.values()) == len(adj)
    return dict(rates)


def check_path_densities():
    totals = Counter()
    inside = Side((-1,), (1,))
    exteriors = (Side((-1,), (1,)),
                 Side((-1, 0), (0, 2)),
                 Side((-1, 0, 0), (1, 1, 2)))
    durations = (F(2, 7), F(4, 5), F(3, 8), F(5, 9), F(1, 2))

    def visit(sides, depth, coefficient, common, gains, direct_survival, split_survival):
        local = [intrinsic(side) for side in sides]
        r, a = local[0]
        out, b = local[1]
        direct = full_rates(sides)
        projected = {(branch, i, event): rate for branch, (rates, _) in enumerate(local)
                     for (i, event), rate in rates.items()}
        projected[0, 0, 'gain'] += b
        projected[1, 0, 'gain'] += a
        assert projected == direct
        dt = durations[depth]
        direct_survival += dt * sum(direct.values())
        split_survival += dt * (sum(r.values()) + sum(out.values()) + a + b)
        assert direct_survival == split_survival

        # Compare coefficients of the same joint path density including its last
        # outside arrow. The survival exponents were checked independently above.
        full_product = prod(left + right for left, right in gains)
        by_last = prod(left for left, _ in gains)
        for j, (_, last) in enumerate(gains):
            by_last += (last * prod(left + right for left, right in gains[:j])
                        * prod(left for left, _ in gains[j+1:]))
        assert coefficient == common * full_product == common * by_last
        totals['paths'] += 1
        if len(gains) >= 2 and len({b for _, b in gains}) > 1:
            totals['paths with changing donor between root gains'] += 1
        if depth == 4:
            return
        for (branch, vertex, event), rate in direct.items():
            new_sides = list(sides)
            side = sides[branch]
            new_sides[branch] = side.gain(vertex) if event == 'gain' else side.promote(vertex)
            root_gain = (branch, vertex, event) == (0, 0, 'gain')
            new_gains = gains + ((r[0, 'gain'], b),) if root_gain else gains
            visit(tuple(new_sides), depth + 1, coefficient * rate,
                  common if root_gain else common * projected[branch, vertex, event],
                  new_gains, direct_survival, split_survival)

    for outside in exteriors:
        visit((inside, outside), 0, F(1), F(1), (), F(0), F(0))
    assert totals['paths with changing donor between root gains'] > 0
    return totals


def check_profiles():
    checks = 0
    for degrees in combinations_with_replacement(range(2, 9), 7):
        b_path = [F(1, degree) for degree in degrees]
        for grid in (2, 4, 10):
            error = F(0)
            for i, b in enumerate(b_path):
                left = i - i % 3
                beta = F(ceil(grid * b_path[left]), grid)
                assert beta >= b
                error += (beta - b) / 3
                assert error <= F(i+1, 3*grid) + F(1, 2)
                checks += 1
    return checks


def check_likelihood():
    taylor = 0
    for d, B, beta, K, extra in product((3, 4, 8, 16), (2, 5, 97),
                                         (F(0), F(1, 7), F(1, 2)),
                                         (1, 2, 8, 100), (F(0), F(1, 3), F(9, 5))):
        r, gamma = K + extra, F(1, B)
        u = (r + beta)/(r + gamma)
        residual = (r + gamma) * (u**d - 1) - d * (beta - gamma)
        bound = F(d*(d-1), 8*K) * F(3, 2)**(d-2)
        assert 0 <= residual <= bound
        taylor += 1

    counts = 0
    for K, M in ((1, 3), (2, 4), (3, 5), (5, 8)):
        for changes in product((1, -1), repeat=10):
            degree = 2
            leaves = 1
            low_gains = rare_promotions = 0
            for change in changes:
                if change == -1 and not leaves:
                    break
                if change == 1:
                    low_gains += leaves < K
                    degree += 1
                elif degree >= M and leaves <= K:
                    rare_promotions += 1
                leaves += change
                assert low_gains <= M + K + rare_promotions
                counts += 1
    return taylor, counts


def check_donor():
    def pgf(death, decay, z):
        zero = death*(1-decay)/(1-death*decay)
        q = (1-death)*decay/(1-death*decay)
        if q == 0:
            return zero
        return zero + (1-zero)*q*z/(1-(1-q)*z)

    checks = 0
    for b, x, y, z in product((2, 3, 9, 100), (F(1), F(2, 3), F(1, 101)),
                              (F(1), F(1, 5), F(5, 9)), (F(0), F(1, 2), F(9, 10))):
        death = F(1, b)
        assert pgf(death, x, pgf(death, y, z)) == pgf(death, x*y, z)
        assert pgf(death, F(1), z) == z
        assert pgf(death, F(0), z) == death
        checks += 1
    return checks


def check_parameters():
    for j in range(4, 65):
        gamma = F(1, j)
        c_d = F(j*(j-1), 8)*F(3, 2)**(j-2)
        z = (1+F(j, 2))**j
        K = max(1, ceil(j*c_d))
        M = K+1+ceil(j*(z-1)*K)
        V = j*gamma + c_d/K + (z-1)*K/M
        loss = F(1, 2*j) + F(1, 2*j) + V/j + F(2, j) + gamma + F(j-1, j*j)
        assert loss <= F(6, j) + F(2, j*j)
    return 61


def check_free_cut():
    # Direct transition rates into initially unmarked vertices, marks held fixed.
    def drift(graph, marked):
        return sum(F(1, len(ns)) for ns in graph.values() for y in ns
                   if y not in marked and len(graph[y]) == 1)

    checks = 0
    for K in range(1, 33):
        m = 2*K
        g = {i: [] for i in range(2*m+1)}
        for y in range(1, m+1):
            g[0].append(y)
            g[y] = [0, m+y]
            g[m+y] = [y]
        marked = set(range(m+1))
        original = drift(g, marked)
        for retained in range(K):
            cut = {v: ns.copy() for v, ns in g.items()}
            cut[0] = list(range(1, retained+1))
            for y in range(retained+1, m+1):
                cut[y].remove(0)
            difference = drift(cut, marked) - original
            assert difference == F(m-retained, 2)
            assert difference/(m+1) > F(1, 4)
            checks += 1
    return checks


if __name__ == '__main__':
    for label, value in check_path_densities().items():
        print(f'{label} = {value}')
    print('rounded monotone profile prefixes =', check_profiles())
    taylor, counts = check_likelihood()
    print('Renyi Taylor rows =', taylor)
    print('low-leaf gain-count prefixes =', counts)
    print('birth-death PGF semigroup rows =', check_donor())
    print('explicit parameter scales j=4..64 =', check_parameters())
    print('free-cut exact drift rows =', check_free_cut())
    print('all asserted rational equalities and inequalities passed')
