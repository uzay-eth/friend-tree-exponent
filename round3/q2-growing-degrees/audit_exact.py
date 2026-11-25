"""Independent finite controls for the frozen-mean audit and its Gaussian addendum."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import comb


def gamma_control():
    eta = F(2, 5)
    x = eta / 4
    g1, g2 = F(1099, 1000), F(89199, 1000)
    assert g1 < 1 + x and g2 > (1 - x) * 99
    excess = g1 / (g1 + g2) - (1 + eta / 2) / 100
    assert excess == F(964, 5643625)
    # The rational identities are polynomial certificates in eta.
    # With 0<eta<1/2 both right sides are positive.
    for j in range(1, 200):
        eta = F(j, 400)
        x = eta / 5
        assert (1 + eta / 2) * (1 - x) - (1 + x) == eta * (1 - eta) / 10
        assert (1 - x) - (1 - eta / 2) * (1 + x) == eta * (1 + eta) / 10
    print(f"strict gamma good-event counterexample: excess {excess}")
    print("199 rational substitutions in the two repaired comparison identities")


@lru_cache(maxsize=None)
def neutral(u, n, j, k):
    if u == n:
        return F(j == k)
    prob = F(j, u)
    return (1 - prob) * neutral(u + 1, n, j, k) + prob * neutral(u + 1, n, j + 1, k)


def choose(n, d):
    return comb(n, d) if 0 <= d <= n else 0


def column_check():
    checks = 0
    for n in range(5, 18):
        for u in range(3, n + 1):
            for k in range(2, n):
                # Independently remove the target rank, then draw u-2 ranks.
                total = comb(n - 2, u - 2)
                sample = [F(choose(k - 1, j - 1) * choose(n - k - 1, u - j - 1), total)
                          for j in range(1, u)]
                assert sum(sample) == 1
                for j, p in enumerate(sample, 1):
                    assert neutral(u, n, j, k) == F(u - 1, n - 1) * p
                stat = sum(j * p for j, p in enumerate(sample))
                var = sum((j - stat)**2 * p for j, p in enumerate(sample))
                assert var == F((u - 2) * (k - 1) * (n - k - 1) * (n - u),
                                (n - 2)**2 * (n - 3))
                checks += 1
    print(f"{checks} complete exact columns and variances, including deterministic terminal columns")


def seed(parents):
    rows = [set() for _ in range(len(parents) + 1)]
    for v, w in enumerate(parents, 1):
        rows[v].add(w)
        rows[w].add(v)
    return tuple(map(frozenset, rows))


def attach(tree, v):
    return tree[:v] + (tree[v] | {len(tree)},) + tree[v + 1:] + (frozenset({v}),)


def rates(tree):
    return [sum((F(1, len(tree[w])) for w in row), F()) for row in tree]


def one_case(tree, k, n):
    leaves = nodes = marked_tests = 0

    def walk(t, R, V, start_score):
        nonlocal leaves, nodes, marked_tests
        u = len(t)
        d = [len(row) for row in t]
        if u == n:
            leaves += 1
            M = d.count(k) - start_score - R
            return M, M * M - V
        nodes += 1
        a = rates(t)
        e = [deg - weight for deg, weight in zip(d, a)]
        choices = [F(weight, u) for weight in a]
        totals = Counter()
        for degree, prob in zip(d, choices):
            totals[degree] += prob

        def weight(j):
            return neutral(u + 1, n, j, k)

        gradient = {j: weight(j + 1) - weight(j) for j in set(d) if j > 1}
        gradient[1] = weight(2)
        raw = sum(totals[j] * gradient[j] for j in totals)
        bias = totals[1] * weight(2) - sum(
            e[i] * gradient[d[i]] / u for i in range(u) if d[i] >= 2)
        bracket = sum(p * (gradient[j] - raw)**2 for j, p in totals.items())
        score = sum(neutral(u, n, j, k) for j in d if j >= 2)
        means = [F(), F()]
        observed_drift = [F()] * u
        for v, prob in enumerate(choices):
            child = attach(t, v)
            after = sum(weight(len(row)) for row in child if len(row) >= 2)
            assert after - score == gradient[d[v]] - raw + bias
            # Check the proposed weighted degree-slice generator directly.
            a1 = rates(child)
            for root in range(u):
                if d[root] >= 2:
                    ep = len(child[root]) - a1[root]
                    test = ep * (len(child[root]) == k) - e[root] * (d[root] == k)
                    observed_drift[root] += prob * u * test
            terms = walk(child, R + bias, V + bracket, start_score)
            for j in range(2):
                means[j] += prob * terms[j]
        for v in range(u):
            if d[v] < 2:
                continue
            leaf_count = sum(d[w] == 1 for w in t[v])
            h = F(leaf_count, 2 * d[v]) + sum(
                a[w] / (d[w] * (d[w] + 1)) for w in t[v] if d[w] >= 2)
            predicted = a[v] * e[v] * ((d[v] == k - 1) - (d[v] == k)) + h * (d[v] == k)
            assert observed_drift[v] == predicted
            marked_tests += 1
        return tuple(means)

    start_score = sum(neutral(len(tree), n, len(row), k) for row in tree if len(row) >= 2)
    assert walk(tree, F(), F(), start_score) == (0, 0)
    return leaves, nodes, marked_tests


if __name__ == "__main__":
    gamma_control()
    column_check()
    cases = [one_case(seed(parents), k, 9)
             for parents in ((0, 0, 0, 0, 0), (0, 1, 2, 3, 4), (0, 0, 1, 2, 1))
             for k in (3, 4, 7)]
    leaves, nodes, sliced = map(sum, zip(*cases))
    print(f"{len(cases)} physical neutral-kernel cases: {leaves} paths, {nodes} internal histories")
    print(f"compensated first and second moments exact; {sliced} tagged slice-generator checks")
