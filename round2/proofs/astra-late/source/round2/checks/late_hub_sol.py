"""Exact finite controls for the deterministic late-label hub proof."""

from collections import defaultdict
from fractions import Fraction as Q


def arrows(parents):
    n = len(parents) + 1
    adj = [set() for _ in range(n)]
    for v, u in enumerate(parents, 1):
        adj[u].add(v)
        adj[v].add(u)
    return adj, {(s, t): Q(1, len(adj[s])) for s in range(n) for t in adj[s]}


def promotion_tiny_time(parents, target, remove_tag=False):
    # Report p'(0),p''(0) for promotion of the newly born leaf. Compute
    # from the rate once attached, allowing one jump of the suppressed chain.
    placed = parents+(target,)
    adj, moves = arrows(placed)
    tag = len(placed)
    assert adj[tag] == {target}
    hazard = moves.pop((target, tag))
    if remove_tag:
        moves.pop((tag, target))
    drift = Q(0)
    for (_, t), r in moves.items():
        _, nxt = arrows(placed+(t,))
        drift += r*(nxt[target, tag]-hazard)
    return hazard, drift-hazard**2


def check_birth_flux():
    # Physical tree: 0 has three leaf neighbours and one core neighbour 3;
    # vertex 3 also has two leaves. Retain labels in the one-jump law.
    pre = (0, 0, 0, 3, 3, 0)
    adj, moves = arrows(pre)
    flat = defaultdict(Q)
    for (_, t), rate in moves.items():
        flat[pre+(t,)] += rate
    split, wrong = defaultdict(Q), Q(0)
    for s, ns in enumerate(adj):
        if len(ns) == 1:
            continue
        for t in ns:
            split[pre+(t,)] += Q(1,len(ns))
            if len(adj[t]) == 1:
                wrong += Q(1,len(ns))  # counting immediate old-leaf promotions
        split[pre+(s,)] += sum(len(adj[t])==1 for t in ns)
    assert split == flat and sum(flat.values()) == 7
    assert wrong == Q(17,12) > 0  # newly born tags have zero mass at delay zero
    # The same absolute weighted small-time bound covers arbitrary leaf stock.
    for k in range(2, 51):
        star = (0,)*k
        slope, _ = promotion_tiny_time(star, 0)
        assert k*slope == Q(k,k+1) < 1
    return wrong


def pgfs(x, z):
    # x=e^(t/2); one-family transform, and immigrant transform.
    den = 1+2*(x-1)*(1-z)
    return 1-x*(1-z)/den, 1/den


def check_pgf():
    rows = 0
    for x in (Q(1), Q(5,4), Q(3,2), Q(2), Q(8)):
        for z in (Q(0), Q(1,3), Q(2,3), Q(1)):
            g, imm = pgfs(x,z)
            f = 1-g
            den = 1+2*(x-1)*(1-z)
            df = x*(1-z)/(2*den)-((x*(1-z))/den)**2
            assert df == (x/2)*(1-z)*(1-2*(1-z))/den**2
            assert -f*imm == -x*(1-z)/den**2
            for y in (Q(1), Q(3,2), Q(3)):
                gy, iy = pgfs(y,z)
                gxgy, ixgy = pgfs(x,gy)
                assert (gxgy, ixgy*iy)==pgfs(x*y,z)
                rows += 1
            assert 0 <= g <= 1 and 0 < imm <= 1
    assert pgfs(Q(2), Q(0))[1] == Q(1,3)  # no immigrants would give 1
    return rows


if __name__ == '__main__':
    wrong = check_birth_flux()
    # Start from T3; place at an old leaf. The new tag has a degree-two
    # parent, not an immediately promoted tag. Before its first promotion
    # only two channels can grow the parent: the tag and the old centre.
    physical = promotion_tiny_time((0,0),1)
    dropped = promotion_tiny_time((0,0),1,remove_tag=True)
    assert physical == (Q(1,2), -Q(1,2))
    assert dropped == (Q(1,2), -Q(1,3))
    print(f'new-label zero-age flux=0; wrong old-target flux={wrong}; star weighted density at zero <1')
    print(f'T3 leaf placement p\'(0),p\'\'(0): {physical}; deleting tag initiation: {dropped}')
    print(f'{check_pgf()} exact family/immigration semigroup tests passed')
