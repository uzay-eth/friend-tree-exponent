"""Independent exact arrow-tree audit of the fixed-interior change of law."""
from fractions import Fraction as Q
from itertools import product


def tree(parents):
    adj = tuple(() for _ in range(len(parents) + 1))
    out = list(adj)
    for i, p in enumerate(parents, 1):
        out[i] += (p,)
        out[p] += (i,)
    return tuple(out)


def birth(adj, target):
    return adj[:target] + (adj[target] + (len(adj),),) + adj[target + 1:] + ((target,),)


def arrows(adj):
    for sender, neighbours in enumerate(adj):
        for target in neighbours:
            yield sender, target, Q(1, len(neighbours))


def potentials(adj, tag, a):
    R = len(adj[tag]) - 1 + a
    B = len(adj) - R
    delta = R - sum(Q(1, len(adj[x])) for x in adj[tag])
    return R, B, delta


def descend(adj, tag, a, depth, pm=Q(1), qm=Q(1), likelihood=Q(1), count=0, out=0):
    assert pm == qm * likelihood
    if depth == 0:
        return 1
    n = len(adj)
    R, B, delta = potentials(adj, tag, a)
    assert delta >= 0 and delta <= 1 + count
    assert B == tag + 1 - a + out
    listed = list(arrows(adj))
    to_tag = sum(w for _, y, w in listed if y == tag)
    outside = sum(w for _, y, w in listed if y != tag)
    assert (to_tag, outside) == (R - delta, B + delta)
    assert sum(w for x, _, w in listed if x == tag) == 1
    for theta in (Q(1, 7), Q(1, 2), Q(6, 7)):
        signed = theta * (-delta / R) + (1 - theta) * delta / B
        assert signed == delta * (R / n - theta) * n / (R * B)
    # Keep outside winners as directed arrows even when histories coalesce.
    paths = descend(birth(adj, tag), tag, a, depth - 1, pm * to_tag / n,
                    qm * R / n, likelihood * (1 - delta / R), count, out)
    qout = Q(0)
    for x, y, w in listed:
        if y == tag:
            continue
        qstep = B / n * w / outside
        qout += qstep
        paths += descend(birth(adj, y), tag, a, depth - 1, pm * w / n,
                         qm * qstep, likelihood * (1 + delta / B), count + (x == tag), out + 1)
    assert qout == B / n
    return paths


def plants():
    adj = tree((0, 1, 2))  # leaf placement: new tag 3 has parent 2 of post-degree 2
    assert potentials(adj, 3, Q(1, 2)) == (Q(1, 2), Q(7, 2), 0)
    assert potentials(birth(adj, 2), 3, Q(1, 2)) == (Q(1, 2), Q(9, 2), Q(1, 6))
    after_gain = birth(adj, 3)
    assert potentials(birth(after_gain, 4), 3, Q(1, 2)) == (Q(3, 2), Q(9, 2), Q(1, 2))
    # High-degree recipient: almost all birth mass, but only order 1/n weighted mass.
    for n in (4, 8, 32, 128):
        star = tree((0,) * (n - 1))
        raw = sum(w / n for _, y, w in arrows(star) if y == 0)
        assert raw == Q(n - 1, n)
        assert raw / len(birth(star, 0)[0]) == Q(n - 1, n * n)
        assert sum(w / n for _, y, w in arrows(star) if y != 0) == Q(1, n)


if __name__ == '__main__':
    paths = placements = 0
    for n in range(2, 6):
        for parents in product(*(range(i) for i in range(1, n))):
            old = tree(parents)
            rates = [sum(w for _, y, w in arrows(old) if y == x) for x in range(n)]
            cores = [x for x in range(n) if len(old[x]) >= 2]
            assert sum(rates[x] / (len(old[x]) + 1) for x in cores) <= len(cores)
            if n > 2:
                assert sum(rates[x] for x in range(n) if x not in cores) == sum(
                    Q(sum(len(old[y]) == 1 for y in old[x]), len(old[x])) for x in cores)
            for x in range(n):
                post = birth(old, x)
                a = Q(1, len(post[x]))
                paths += descend(post, n, a, 2)
                placements += 1
    print(f'PASS {placements} tagged placements and {paths} two-step arrow histories; exact dual likelihood and signed drift')
    plants()
    print('PASS tag→old parent, tag→new stock leaf, and four high-degree star plants')
