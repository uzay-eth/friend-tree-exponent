from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache


@lru_cache(None)
def adj(parents):
    out = [[1], [0]] + [[] for _ in parents]
    for v, u in enumerate(parents, 2):
        out[u].append(v)
        out[v].append(u)
    return out


states = {((), None, None): F(1)}
joint = defaultdict(F)
for n in range(2, 9):
    nxt, exits, survival = defaultdict(F), defaultdict(F), F(0)
    for (tree, j, c4), prob in states.items():
        deg = list(map(len, adj(tree)))
        c = sum(d >= 2 for d in deg)
        if n == 4:
            c4 = c
        if n == 5:
            joint[c4, c] += prob
        if n >= 3:
            assert (j is None) == (c == 1)
            survival += prob if j is None else 0
            if j is not None:
                exits[j] += prob
        if n == 8:
            continue
        # Aggregate sender choices by the recipient, preserving exact rates.
        for recv, senders in enumerate(adj(tree)):
            child = (*tree, recv)
            q = sum(F(1, n * deg[u]) for u in senders)
            new_j = n if n >= 3 and j is None and deg[recv] == 1 else j
            if j != new_j:
                new_a = adj(child)
                cores = sorted((len(a), u) for u, a in enumerate(new_a) if len(a) >= 2)
                assert [d for d, _ in cores] == [2, n - 1]
                assert cores[0][1] in new_a[cores[1][1]]
            nxt[child, new_j, c4] += q * prob
    assert sum(states.values()) == 1
    if n >= 3:
        assert survival == F(2, n - 1)
        assert exits == {j: F(2, j * (j - 1)) for j in range(3, n)}
        print(f"n={n}: histories={len(states)} survival={survival}; exits={list(exits.values())}")
    states = nxt
assert list(v for _, v in sorted(joint.items())) == [F(1, 2), F(1, 6), F(1, 4), F(1, 12)]
print("(C4,C5):", sorted(joint.items()))

# Counterexample to the unmodified triple flag. At the start, the middle
# vertex has two core neighbours; a neighbour sends into it before it rings.
p = (1, 2, 3, 4, 5)  # path of seven vertices, triple 2,3,4
def send(tree, source, target):
    a = adj(tree)
    assert target in a[source] and len(a[target]) >= 2  # no promotion
    return (*tree, target)

p = send(p, 4, 3)   # first half, z sends to w
assert adj(p)[3] == [2, 4, 7]
p = send(p, 3, 4)   # first half, middle mark 5/12 selects z, not u
assert len(adj(p)[2]) == 2
p = send(p, 2, 1)   # second half, base mark 1/4 when leaves are ranked first
assert sum(len(a) >= 2 for a in adj(p)) == 5
print("Unmodified triple flag can fire with no promotion (positive-length mark intervals).")

lo, hi = F(565, 1000), F(5665, 10000)
assert 2*lo > 1 and hi < 1
assert 1/F(9, 8) + 1/F(16) + 1/F(144, 7) == 1
for s, cap in [(F(11, 10), 1), (F(2), F(4, 3)), (F(41, 20), F(4, 3))]:
    print(f"s={s}: max size order={F(9,8)*hi*s} < {cap}; inverse power={16*s} < 40")
    assert F(9, 8) * hi * s < cap and 16*s < 40
