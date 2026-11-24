from collections import defaultdict
from fractions import Fraction as F


def levels(last):
    trees = {(0,): F(1)}
    for n in range(2, last+1):
        marginals = [defaultdict(F) for _ in range(n)]
        leaves = [F(0) for _ in range(n)]
        next_trees = {}
        for parents, mass in trees.items():
            adj = [[] for _ in range(n)]
            for v, w in enumerate(parents, 1):
                adj[v].append(w)
                adj[w].append(v)
            d = list(map(len, adj))
            for u, k in enumerate(d):
                marginals[u][k] += mass
                leaves[u] += mass*max(0, sum(d[v] == 1 for v in adj[u])-1)
            if n < last:
                rates = [sum((F(1, d[v]) for v in row), F(0)) for row in adj]
                assert sum(rates) == n
                for u, rate in enumerate(rates):
                    next_trees[parents+(u,)] = mass * rate / n
        assert all(sum(row.values()) == 1 for row in marginals)
        assert sum(sum(d*p for d,p in row.items()) for row in marginals) == 2*(n-1)
        assert marginals[0] == marginals[1]
        yield n, marginals, leaves
        trees = next_trees


if __name__ == '__main__':
    for n, rows, leaves in levels(8):
        if n == 2:
            continue
        print(f'n={n}; degree order 1..{n-1}')
        for u, row in enumerate(rows, 1):
            print(f"  u={u}: {' '.join(str(row[k]) for k in range(1,n))}"
                  f" ; E[D/n]={sum(k*p for k,p in row.items())/n}"
                  f" ; lower bound on E[Z]={leaves[u-1]/n}")
        for a,b in zip(rows, rows[1:]):
            assert all(sum(a[k] for k in range(j,n)) >= sum(b[k] for k in range(j,n))
                       for j in range(1,n))
    print('Total probability, handshake, initial symmetry and adjacent stochastic order checks passed.')
    assert leaves[1]/8 > sum(k*p for k,p in rows[2].items())/8
    print('Exact cross-time brackets certify E[Z_2] > E[Z_3].')
