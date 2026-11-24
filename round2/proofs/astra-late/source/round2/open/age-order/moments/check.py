#!/usr/bin/env python3
"""Independent rational tests of the finite double-placement moment law."""
from functools import lru_cache
from fractions import Fraction as Q
from itertools import product
from math import comb, prod


def add(t, u):
    return t + (u,)


@lru_cache(None)
def rates(t):
    adj = [[] for _ in range(len(t) + 1)]
    for u, v in enumerate(t, 1):
        adj[u].append(v)
        adj[v].append(u)
    deg = tuple(map(len, adj))
    return adj, deg, tuple(sum((Q(1, deg[v]) for v in nn), Q()) for nn in adj)


@lru_cache(None)
def terminal(t, mark, alpha, depth):
    _, deg, r = rates(t)
    n = len(deg)
    if not depth:
        return Q(deg[mark], n)**alpha
    return sum(r[u]/n * terminal(add(t, u), mark, alpha, depth-1) for u in range(n))


def check(t, alpha):
    adj, d, r = rates(t)
    n = len(d)

    def diff(s):
        return terminal(s, n, alpha, 1) - terminal(s, n+1, alpha, 1)

    direct = Q()
    for u in range(n):
        s = add(t, u)
        for w, r_w in enumerate(rates(s)[2]):
            direct += r[u]*r_w*diff(add(s, w))

    chain = sum(r[u]/(d[u]+1)*diff(add(add(t, u), n)) for u in range(n))
    delta = [Q(1, k*(k+1)) for k in d]
    edge = sum((r[w]*delta[w]-r[u]*delta[u])*diff(add(add(t, u), w))
               for u in range(n) for w in adj[u] if u < w)
    assert direct == chain + edge
    return edge != 0


def main():
    rows = nonzero_edges = 0
    for n in range(3, 6):
        for t in product(*(range(i) for i in range(1, n))):
            for alpha in (1, 2, 3):
                nonzero_edges += check(t, alpha)
                rows += 1
    assert nonzero_edges
    print(f'PASS: {rows} exact signed moment comparisons; {nonzero_edges} nonzero edge terms')

    # Beta power integrals by polynomial integration against the density,
    # versus its rising-factorial moment and the adjacent survival identity.
    rows = 0
    for m in range(3, 16):
        for a in range(1, 5):
            moment = (m-1)*(m-2)*sum(Q((-1)**i*comb(m-3, i), a+i+2)
                                      for i in range(m-2))
            assert moment == Q(prod(range(2, a+2)), prod(range(m, m+a)))
            for x in (Q(1, 9), Q(1, 2), Q(8, 9)):
                b = lambda j: (1-x)**(j-2)*(1+(j-2)*x)
                assert b(m)-b(m+1) == (m-1)*x*x*(1-x)**(m-2)
                rows += 1
    print(f'PASS: {rows} exact Beta moment/tail rows at actual promotion sizes')


if __name__ == '__main__':
    main()
