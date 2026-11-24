"""Check branch flux and the no-interior-cut rule by labelled arrow enumeration."""

from fractions import Fraction
from itertools import product


def check():
    cases = asymmetric = 0
    for n in range(4, 8):
        for parents in product(*(range(v) for v in range(1, n))):
            adj = [set() for _ in range(n)]
            for v, p in enumerate(parents, 1):
                adj[v].add(p)
                adj[p].add(v)
            for v, p in enumerate(parents, 1):
                if min(len(adj[v]), len(adj[p])) < 2:
                    continue
                branch, queue = {v}, [v]
                for x in queue:
                    next_nodes = adj[x] - branch - {p}
                    branch.update(next_nodes)
                    queue.extend(next_nodes)
                cross = {(x, y) for x in branch for y in adj[x] - branch}
                assert cross == {(v, p)}
                rate = sum(Fraction(1, len(neighbors)) for neighbors in adj
                           for y in neighbors if y in branch)
                margin = Fraction(1, len(adj[p])) - Fraction(1, len(adj[v]))
                assert rate == len(branch) + margin
                assert abs(margin) <= Fraction(1, 2)
                asymmetric += margin != 0
                assert all(len(adj[x]) <= len(branch) for x in branch)
                assert all(len(adj[x]) == len(adj[x] & branch) + (x == v)
                           for x in branch)
                cases += 1
    assert asymmetric > 0
    print(f'{cases} physical core-edge branches; {asymmetric} unequal-flux cases; all PASS')


if __name__ == '__main__':
    check()
