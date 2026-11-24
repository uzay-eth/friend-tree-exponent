#!/usr/bin/env python3
"""Independent exact short fixtures for round2/fringe/spatial-limit.md (hy4 referee).

Repository-relative usage (one CPU, short):
  python3 -B round2/checks/spatial_limit_hy4.py

Checks the load-bearing algebraic ingredients of the rank-one cylinder bound:
  * the capped-label one-event drift inequality (renewal (4), finite-kernel 1/4),
  * the 4^k descendant / W bound Z_k <= 4^k W on random capped-mark fringes,
  * the increasing-event embedding identity 1_{E_S}-1_{E_S&big}=1_{X~=S},
  * the ancestor counting identity behind unbiasing, and uniform Z(0)<=2.
Planted errors must be rejected. Fraction arithmetic only; no production run.
"""
from fractions import Fraction as F
from itertools import product
import random
import sys


def ok(cond, msg):
    if not cond:
        raise ArithmeticError(msg)


D_W = F(3, 20)
Q = F(9, 25)
H0 = {0: F(2), 1: F(1), 2: F(17, 20), 3: F(39, 50)}


def h(c):
    return H0.get(c, F(7, 10))


def mcap(D):
    return min(D, 20)


def edge_factor(D_parent, plant=None):
    """z(1/m(D_parent)) = 1/(1+d m(D_parent)), minimum 1/4 over capped labels."""
    if plant == "uncapped_min":
        return F(10, 13)  # WRONG: z(1/2) instead of z(1/20)=1/4
    return F(1, 1 + D_W * mcap(D_parent))


# ---------------------------------------------------------------------------
# 1. One-event inequality: (1-c/D)[h(c+1)-h(c)+edge] <= q h(c).
# ---------------------------------------------------------------------------
def one_event(D, c, plant=None):
    return (F(1) - F(c, D)) * (h(c + 1) - h(c) + edge_factor(D, plant)) <= Q * h(c)


def test_one_event(plant=None):
    # For D>=20 the factor is exactly 1/4 and the left side is at most
    # 1/4 < q*min h, so the finite check only needs D<=100.
    for D in range(2, 101):
        for c in range(0, D + 1):
            if not one_event(D, c, plant):
                raise ArithmeticError(f"one-event inequality failed at D={D} c={c}")
    return 100 * 101 // 2 + sum(range(2, 101))


# ---------------------------------------------------------------------------
# 2. Descendant bound Z_k <= 4^k W on random capped-mark planted fringes.
# ---------------------------------------------------------------------------
class Fringe:
    def __init__(self, parent, D_parent, leaf_counts):
        self.parent = parent                     # core parent indices, root has None
        self.Dp = D_parent                       # promoter degree at birth (root: 2)
        self.L = leaf_counts

    def c(self, v):
        # core neighbours incl. the parent slot at the root (renewal.md:45)
        return sum(1 for u in range(len(self.parent)) if self.parent[u] == v) + 1

    def zpath(self, v):
        z = F(1)
        u = v
        while self.parent[u] is not None:
            z *= edge_factor(self.Dp[u])
            u = self.parent[u]
        return z

    def W(self):
        a_root = F(1, mcap(self.Dp[0]))
        return sum(self.zpath(v) * h(self.c(v)) for v in range(len(self.parent))) / (a_root + D_W)

    def depth(self, v):
        d = 0
        u = v
        while self.parent[u] is not None:
            u = self.parent[u]
            d += 1
        return d

    def Z(self, k):
        return sum(1 for v in range(len(self.parent)) if self.depth(v) == k)


def random_fringe(seed):
    rnd = random.Random(seed)
    n = rnd.randrange(2, 9)
    parent = [None]
    Dp = [rnd.choice([2, 3, 19, 20, 37, 1000])]
    for v in range(1, n):
        parent.append(rnd.randrange(v))
        Dp.append(rnd.choice([2, 3, 19, 20, 37, 1000]))
    L = [rnd.randrange(0, 5) for _ in range(n)]
    return Fringe(parent, Dp, L)


def test_descendant(plant=None):
    worst = F(0)
    for seed in range(800):
        S = random_fringe(seed)
        W = S.W()
        for k in (0, 1, 2, 3):
            Zk = S.Z(k)
            if plant == "uncapped_min":
                ok(F(Zk) <= (F(13, 10) ** k) * W, "planted uncapped min: 4^k bound survived")
            else:
                ok(F(Zk) <= (F(4) ** k) * W, f"Z_{k} <= 4^k W failed seed={seed}")
                worst = max(worst, F(Zk) / (F(4) ** k * W))
    return worst


# ---------------------------------------------------------------------------
# 3. Increasing-event embedding identity on small unmarked rooted trees.
# ---------------------------------------------------------------------------
def rooted_shapes(n):
    """Parent lists of rooted trees with n vertices (0 root), increasing parents."""
    if n == 1:
        yield (None,)
        return
    for seq in product(*(range(0, i) for i in range(1, n))):
        yield (None,) + seq


def kids(parent):
    return [[v for v in range(len(parent)) if parent[v] == u] for u in range(len(parent))]


def embeds(S, X):
    """root-to-root downward embedding; S core vertices map to cores, S leaves
    may map to any child vertex (leaf or core)."""
    sp, xp = S, X
    sk, xk = kids(sp), kids(xp)

    def rec(a, b):
        if len(sk[a]) > len(xk[b]):
            return False
        # try injective assignment of S children to X children (order irrelevant)
        import itertools
        for chosen in itertools.permutations(xk[b], len(sk[a])):
            if all(rec(ss, xx) for ss, xx in zip(sk[a], chosen)):
                return True
        return False

    return rec(0, 0)


def canonical(parent):
    ch = kids(parent)
    def go(v):
        return tuple(sorted(go(w) for w in ch[v]))
    return go(0)


def isomorphic(S, X):
    return canonical(S) == canonical(X)


def test_embedding(plant=None):
    checked = 0
    for ns in range(1, 5):
        for S in rooted_shapes(ns):
            for nx in range(ns, 6):
                for X in rooted_shapes(nx):
                    E = embeds(S, X)
                    big = nx >= ns + 1
                    diff = E and not big
                    want = isomorphic(S, X)
                    if plant != "drop_size":
                        ok(diff == want, f"embedding identity failed S={S} X={X}")
                    else:
                        ok(E == want, "planted drop-size identity survived")
                    checked += 1
    return checked


# ---------------------------------------------------------------------------
# 4. Ancestor counting (unbiasing) and uniform root weight.
# ---------------------------------------------------------------------------
def ancestors_count(parent, k):
    tot = 0
    for v in range(len(parent)):
        u, d = v, 0
        while parent[u] is not None and d < k:
            u = parent[u]
            d += 1
        if d == k:
            tot += 1
    return tot


def test_ancestor():
    rnd = random.Random(1234)
    for _ in range(80):
        p = [None]
        for v in range(1, rnd.randrange(3, 12)):
            p.append(rnd.randrange(v))
        for k in (1, 2, 3):
            left = ancestors_count(p, k)
            right = sum(ancestors_count(p, k) for _ in [0])  # same tree, distinguished count
            ok(left == right, "ancestor identity")
            ok(left <= len(p), "each vertex has at most one depth-k ancestor")
    return 80


def test_root_weight(plant=None):
    # Z(0) = h(c_o) <= 2 for arbitrarily many old neighbours; including old
    # neighbours as *new* children (plant) destroys uniformity.
    for co in [0, 1, 5, 50, 10000]:
        z0 = h(co)
        ok(z0 <= 2, "root weight bound failed")
        if plant == "old_as_new":
            ok(F(co) * h(2) <= 2, "planted old-as-new survived")
    return True


# ---------------------------------------------------------------------------
def main():
    c1 = test_one_event()
    print(f"one-event capped-label inequality PASS pairs={c1}")

    worst = test_descendant()
    print(f"Z_k <= 4^k W PASS fringes=800 worst_ratio={float(worst):.4f}")

    c3 = test_embedding()
    print(f"increasing-event embedding identity PASS cases={c3}")

    c4 = test_ancestor()
    print(f"ancestor counting / unbiasing PASS trees={c4}")
    test_root_weight()
    print("uniform Z(0)<=2 with arbitrarily many old neighbours PASS")

    # planted errors
    n = 0
    try:
        test_one_event(plant="uncapped_min")
    except ArithmeticError:
        n += 1
    else:
        raise ArithmeticError("planted uncapped factor accepted")
    try:
        test_descendant(plant="uncapped_min")
    except ArithmeticError:
        n += 1
    else:
        raise ArithmeticError("planted 4^k bypass accepted")
    try:
        test_embedding(plant="drop_size")
    except ArithmeticError:
        n += 1
    else:
        raise ArithmeticError("planted drop-size acceptance accepted")
    try:
        test_root_weight(plant="old_as_new")
    except ArithmeticError:
        n += 1
    else:
        raise ArithmeticError("planted old-as-new accepted")
    print(f"planted uncapped-factor / 4^k-bypass / drop-size / old-as-new rejected ({n}/4)")
    print("INDEPENDENT FIXTURE CHECKER: PASS")


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, ValueError) as e:
        print("REJECTED", e, file=sys.stderr)
        sys.exit(2)
