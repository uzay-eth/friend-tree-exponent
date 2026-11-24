"""Exact controls for an explicit every-old-tree comparison."""
from fractions import Fraction as Q
from math import comb


def grow(g, u):
    return tuple(row + (len(g),) if k == u else row for k, row in enumerate(g)) + ((u,),)


def graph(edges):
    g = [[] for _ in range(1 + max(max(edge) for edge in edges))]
    for u, v in edges:
        g[u].append(v)
        g[v].append(u)
    return tuple(map(tuple, g))


def rates(g):
    return [sum(Q(1, len(g[w])) for w in row) for row in g]


def check_stock():
    g = graph([(0, 1), (0, 2), (0, 3), (1, 4), (4, 5), (2, 6), (3, 7)])
    r = rates(g)
    w = [q / (len(row) + 1) for q, row in zip(r, g)]
    s = [q / len(row) for q, row in zip(w, g)]
    assert s[1] == Q(5, 36) > s[0] == Q(1, 8)
    assert sum(abs(s[u] - s[v]) for u, row in enumerate(g) for v in row if u < v) <= sum(w)
    parallel = grow(grow(g, 0), 1)
    assert [u for u in parallel[1] if len(parallel[u]) == 1] == [9]
    for n in (8, 64, 1024):
        old = graph([(0, j) for j in range(1, n)])
        middle = grow(grow(old, 0), n)
        r = rates(middle)
        d = len(middle[0])
        red = 1 + Q(1, d)
        assert r[n] == red and r[n+1] == Q(1, 2)
        assert r[0] / (len(middle) - red) > 1 - Q(2, n)
        assert r[0] * Q(1, d * (d + 1)) <= Q(1, 3)


def compositions(s, length):
    if length == 1:
        yield (s,)
        return
    for k in range(s + 1):
        for rest in compositions(s - k, length - 1):
            yield (k,) + rest


def check_compositions():
    cases = 0
    for f in range(1, 6):
        for s in range(13):
            count = comb(s + f, f)
            observed = [0] * s
            for word in compositions(s, f + 1):
                for j in range(s - word[-1], s):
                    observed[j] += 1
            assert observed == [comb(j + f, f) for j in range(s)]
            assert sum((Q(h, j + 1) for j, h in enumerate(observed)), Q(0)) <= Q(count, f)
            cases += 1
    return cases


def check_budgets():
    h = 512
    m = z = 2**80
    b = m - Q(3, 2)
    r = (m + 63) // 64
    q = Q(7, 8)
    # Envelope second moments: sum (j+5)^2 q^j and add first-block cost.
    total = q * (1 + q) / (1 - q)**3 + 10*q / (1 - q)**2 + 25 / (1 - q) + 6
    assert q**4 > Q(1, 2) and total < 2**11
    assert Q(15, 16)**8 > Q(1, 2)
    assert Q(5, 3) * r / b < Q(1, 16)
    assert Q(5, 3) * r / b + Q(9 + 2**11, h*h) < Q(1, 4)
    assert Q(2**13 * h*h, 2**40) + Q((2**14 + 2**7) * h*h, z) < 1
    assert Q(4, 3) * b / (b - 1) * Q(5, 3) < 3  # uses pi^2 < 10
    assert Q(3, 4 * 3**6) > Q(1, 2**10)
    assert 2**90 >= z and Q(2**12, 2**45) < 1


if __name__ == '__main__':
    check_stock()
    print('PASS nonzero zero-stock adverse edge and three high-degree star controls')
    print(f'PASS {check_compositions()} conditional geometric composition bounds')
    check_budgets()
    print('PASS exact dyadic, log-error, Jensen and beta separation budgets for Y=2^90')
