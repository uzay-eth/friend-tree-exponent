"""Exact finite-horizon endpoint controls for the adjacent-label identity."""

from fractions import Fraction as Q
from functools import cache
from itertools import product


@cache
def graph(par):
    adj = [set() for _ in range(len(par)+1)]
    for v, u in enumerate(par, 1):
        adj[v].add(u)
        adj[u].add(v)
    degrees = tuple(map(len, adj))
    r = tuple(sum(Q(1, degrees[u]) for u in ns) for ns in adj)
    return adj, degrees, r


@cache
def hit(par, v, d, h):
    """Exact discrete harmonic recursion: hit degree d in h further births."""
    _, degrees, r = graph(par)
    if degrees[v] >= d:
        return Q(1)
    if h == 0 or degrees[v]+h < d:
        return Q(0)
    return sum(w*hit(par+(u,), v, d, h-1) for u,w in enumerate(r))/len(r)


def verify(par, d, h):
    adj, degree, r = graph(par)
    n = len(r)

    def gap(a, b):
        state = par+(a,b)
        return hit(state,n,d,h)-hit(state,n+1,d,h)

    raw = Q(0)
    for a, rate in enumerate(r):
        _, _, nxt = graph(par+(a,))
        raw += sum(rate*w*gap(a,b) for b,w in enumerate(nxt))
    edge = Q(0)
    by_core = {a: Q(0) for a in range(n) if degree[a]>=2}
    costs = dict(by_core)
    weight = [r[a]/(degree[a]*(degree[a]+1)) for a in range(n)]
    for a in range(n):
        for b in range(a+1,n):
            assert gap(a,b) == -gap(b,a)
            if b in adj[a]:
                edge += (weight[b]-weight[a])*gap(a,b)
        assert gap(a,a)==0
        root = a if degree[a]>=2 else next(iter(adj[a]))
        w = r[a]/(degree[a]+1)
        by_core[root] += w*gap(a,n) - weight[a]*sum(gap(a,b) for b in adj[a])
        costs[root] += w+degree[a]*weight[a]
    serial = sum(r[a]*gap(a,n)/(degree[a]+1) for a in range(n))
    assert raw == edge+serial == sum(by_core.values())
    assert all(abs(by_core[a])<=costs[a]<=3 for a in by_core)
    return raw,edge,serial


if __name__ == '__main__':
    count = 0
    orientations = set()
    for n in range(3,7):
        for par in product(*(range(i) for i in range(1,n))):
            for d,h in ((2,0),(2,1),(3,2)):
                raw,edge,serial = verify(par,d,h)
                if edge:
                    assert raw != -edge+serial
                    orientations.add(edge>0)
                count += 1
        hit.cache_clear()
        graph.cache_clear()
    # T3 is a fixture where the serial contribution already has nonzero
    # mass at zero future horizon: middle is degree 2, tip is degree 1.
    assert verify((0,0),2,0) == (Q(7,6),Q(0),Q(7,6))
    # At threshold two the serial middle is already a success. Comparing a
    # hit at time zero with the true hub probability needs a real survival
    # lower bound, rather than declaring hubhood from present degree.
    print(f'{count} exact harmonic rows passed; nonzero edge signs={sorted(orientations)}')
    print('T3 zero-horizon threshold-2: total=serial=7/6, old-target sector=0')
