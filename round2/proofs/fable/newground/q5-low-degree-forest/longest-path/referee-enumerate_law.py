from collections import Counter
from fractions import Fraction
from math import lcm


def direct(max_n):
    """Enumerate every chronological labelled attachment history from one edge."""
    weights = {n: Counter() for n in range(2, max_n + 1)}
    den = {2: 1}
    scales = {n: lcm(*range(1, n)) for n in range(2, max_n)}
    for n, scale in scales.items():
        den[n + 1] = n * scale * den[n]
    seen_states = Counter()
    adj = [[1], [0]]

    def visit(w):
        n = len(adj)
        seen_states[n] += 1
        unseen = {v for v, a in enumerate(adj) if len(a) == 2}
        largest = 0
        while unseen:
            todo = [unseen.pop()]
            count = 0
            while todo:
                v = todo.pop()
                count += 1
                for x in adj[v]:
                    if x in unseen:
                        unseen.remove(x)
                        todo.append(x)
            largest = max(largest, count)
        weights[n][largest] += w
        if n < max_n:
            scale = scales[n]
            for v, a in enumerate(adj[:n]):
                # Sum the donor-neighbour pair probabilities for recipient v.
                child_weight = w * sum(scale // len(adj[x]) for x in a)
                a.append(n)
                adj.append([v])
                visit(child_weight)
                adj.pop()
                a.pop()

    visit(1)
    for n, row in weights.items():
        assert sum(row.values()) == den[n]
    return ({n: {k: Fraction(w, den[n]) for k, w in row.items()}
             for n, row in weights.items()}, seen_states)


def shape_oracle(max_n):
    """Independent Fraction transition by donor then recipient, quotient by shape.

    Longest path is measured by trying all edge-avoiding walks from each vertex.
    """
    def canonical(adj):
        def code(v, p):
            return '(' + ''.join(sorted(code(x, v) for x in adj[v] if x != p)) + ')'
        return min(code(v, -1) for v in range(len(adj)))

    def longest(adj):
        def walk(v, p):
            if len(adj[v]) != 2:
                return 0
            return 1 + max(walk(x, v) for x in adj[v] if x != p)
        return max(walk(v, -1) for v in range(len(adj)))

    shapes = {'(())': (Fraction(1), [[1], [0]])}
    laws = {}
    shape_counts = {}
    for n in range(2, max_n + 1):
        row = Counter()
        nxt = {}
        for p, adj in shapes.values():
            row[longest(adj)] += p
            if n == max_n:
                continue
            for donor_targets in adj:
                w = p / (n * len(donor_targets))
                for v in donor_targets:
                    child = [a.copy() for a in adj] + [[v]]
                    child[v].append(n)
                    s = canonical(child)
                    if s in nxt:
                        nxt[s] = (nxt[s][0] + w, nxt[s][1])
                    else:
                        nxt[s] = (w, child)
        laws[n] = dict(row)
        shape_counts[n] = len(shapes)
        shapes = nxt
    return laws, shape_counts


if __name__ == '__main__':
    from math import factorial

    laws, histories = direct(10)
    oracle, shapes = shape_oracle(10)
    assert laws == oracle
    for n, row in laws.items():
        assert histories[n] == factorial(n - 1)
        if n >= 3:
            # Only the n-path can have n-2 degree-two vertices.
            # At size >=3 its two leaf-extension rates sum to 1.
            assert row[n - 2] == Fraction(2, factorial(n - 1))
        print(n, histories[n], shapes[n], ' '.join(f'{k}:{row[k]}' for k in sorted(row)))
    print('PASS: labelled enumeration, shape oracle and full-path probabilities agree')
