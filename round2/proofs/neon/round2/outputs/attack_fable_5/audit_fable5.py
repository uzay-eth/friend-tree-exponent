#!/usr/bin/env python3
"""
Independent Gemini-family audit of Fable Claim 5 (ft2-fable@88e13f3):
Stochastic domination between adjacent labels via alternating-depth matching coupling,
exact transport certificates, strictness on [0, 1), and Beta(1/2, v-1/2) upper/lower tails.

Completely from-scratch pure Python (exact Fractions / integers); shares zero code with Fable.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
import math
import sys


# ---------------------------------------------------------------------------
# 1. Rooted tree representation & alternating-depth preorder >=
# ---------------------------------------------------------------------------
# A rooted tree is an immutable sorted tuple of child rooted trees; () is a leaf.

def tree_size(t):
    return 1 + sum(tree_size(c) for c in t)


def adj_to_rooted(adj, root, parent=-1):
    return tuple(sorted(adj_to_rooted(adj, w, root) for w in adj[root] if w != parent))


def rooted_to_adj(t):
    """Return (adj_list, depths) with root = 0."""
    adj = []
    depths = []

    def build(node, par, d):
        u = len(adj)
        adj.append([] if par < 0 else [par])
        depths.append(d)
        if par >= 0:
            adj[par].append(u)
        for c in node:
            build(c, u, d + 1)

    build(t, -1, 0)
    return adj, depths


def all_rooted_trees(max_n):
    """Generate all rooted trees up to size max_n."""
    by_n = {1: [()]}
    for n in range(2, max_n + 1):
        res = set()
        # partitions of n - 1 into sizes of child subtrees
        def gen(rem, max_sz, cur):
            if rem == 0:
                res.add(tuple(sorted(cur)))
                return
            for s in range(min(rem, max_sz), 0, -1):
                for child in by_n[s]:
                    cur.append(child)
                    gen(rem - s, s, cur)
                    cur.pop()
        gen(n - 1, n - 1, [])
        by_n[n] = sorted(res)
    return by_n


def bipartite_matching_exists(left, right, ok_fn):
    """Return an explicit injective mapping dict {i: j} or None."""
    nl, nr = len(left), len(right)
    if nl > nr:
        return None
    match_r = [-1] * nr

    def dfs(u, vis):
        for v in range(nr):
            if not vis[v] and ok_fn(left[u], right[v]):
                vis[v] = True
                if match_r[v] == -1 or dfs(match_r[v], vis):
                    match_r[v] = u
                    return True
        return False

    for u in range(nl):
        if not dfs(u, [False] * nr):
            return None
    inv = {}
    for v, u in enumerate(match_r):
        if u != -1:
            inv[u] = v
    return inv


@lru_cache(None)
def ge_even(a, b):
    """E(a, b): injection from children(b) into children(a) with O(phi(b'), b')."""
    return bipartite_matching_exists(b, a, lambda bc, ac: le_odd(ac, bc)) is not None


@lru_cache(None)
def le_odd(a, b):
    """O(a, b): injection from children(a) into children(b) with E(a', psi(a'))."""
    return bipartite_matching_exists(a, b, lambda ac, bc: ge_even(ac, bc)) is not None


def extract_matching(a, b):
    """
    Construct an explicit matched subtree M between rooted adjacency representations
    of a and b when ge_even(a, b) holds.
    Returns (adjA, depA, adjB, depB, matched_pairs) where matched_pairs is a list
    of (uA, uB, parity) with parity = 0 (even) or 1 (odd).
    """
    adjA, depA = rooted_to_adj(a)
    adjB, depB = rooted_to_adj(b)

    # Also store the rooted subtree at each vertex id for matching
    subA = [None] * len(adjA)
    subB = [None] * len(adjB)
    chA = [[] for _ in range(len(adjA))]
    chB = [[] for _ in range(len(adjB))]

    def fill(adj, sub, ch, u, par):
        cs = [w for w in adj[u] if w != par]
        ch[u] = cs
        for w in cs:
            fill(adj, sub, ch, w, u)
        sub[u] = tuple(sorted(sub[w] for w in cs))

    fill(adjA, subA, chA, 0, -1)
    fill(adjB, subB, chB, 0, -1)

    pairs = []

    def rec(uA, uB, parity):
        pairs.append((uA, uB, parity))
        if parity == 0:
            # inject chB[uB] into chA[uA] with le_odd(subA[a'], subB[b'])
            mp = bipartite_matching_exists(
                chB[uB], chA[uA], lambda vb, va: le_odd(subA[va], subB[vb])
            )
            assert mp is not None
            for idx_b, idx_a in mp.items():
                rec(chA[uA][idx_a], chB[uB][idx_b], 1)
        else:
            # inject chA[uA] into chB[uB] with ge_even(subA[a'], subB[b'])
            mp = bipartite_matching_exists(
                chA[uA], chB[uB], lambda va, vb: ge_even(subA[va], subB[vb])
            )
            assert mp is not None
            for idx_a, idx_b in mp.items():
                rec(chA[uA][idx_a], chB[uB][idx_b], 0)

    rec(0, 0, 0)
    return adjA, depA, adjB, depB, pairs


def attachment_rates(adj):
    deg = [len(nb) for nb in adj]
    return [sum(F(1, deg[w]) for w in adj[v]) for v in range(len(adj))]


# ---------------------------------------------------------------------------
# 2. Check Part (a) & (b): Preorder, One-Step Coupling, Corollary 1
# ---------------------------------------------------------------------------
def check_parts_ab():
    by_n = all_rooted_trees(7)
    counts = {n: len(by_n[n]) for n in range(1, 8)}
    # OEIS A000081: 1, 1, 2, 4, 9, 20, 48
    assert counts == {1: 1, 2: 1, 3: 2, 4: 4, 5: 9, 6: 20, 7: 48}, f"Bad counts: {counts}"
    print(f"[PART A] Rooted tree census n=1..7 verified (OEIS A000081): {counts}")

    # Reflexivity and transitivity on all trees of size <= 6 (37 shapes -> 50,653 triples)
    small = [t for n in range(1, 7) for t in by_n[n]]
    for t in small:
        assert ge_even(t, t), "Reflexivity failed"
    trans_checks = 0
    for a in small:
        for b in small:
            if not ge_even(a, b):
                continue
            for c in small:
                if ge_even(b, c):
                    assert ge_even(a, c), "Transitivity failed!"
                    trans_checks += 1
    print(f"[PART A] Reflexivity (37 shapes) and transitivity ({trans_checks} comparable triples across mixed sizes n<=6) PASSED")

    # Exhaustive one-step coupling lemma check on ALL comparable equal-size pairs for n = 3..6
    # and a broad sample for n = 7
    coupled_pairs_checked = 0
    coupled_outcomes_checked = 0
    for n in range(3, 7):
        trees = by_n[n]
        for a in trees:
            for b in trees:
                if not ge_even(a, b):
                    continue
                adjA, depA, adjB, depB, pairs = extract_matching(a, b)
                rA = attachment_rates(adjA)
                rB = attachment_rates(adjB)
                assert sum(rA) == n and sum(rB) == n

                # Check degree alternation and rate alternation on every matched pair
                matchedA = set()
                matchedB = set()
                common_mass = F(0)
                resA = [F(0)] * n
                resB = [F(0)] * n
                for uA, uB, par in pairs:
                    matchedA.add(uA)
                    matchedB.add(uB)
                    if par == 0:
                        assert len(adjA[uA]) >= len(adjB[uB])
                        assert rA[uA] >= rB[uB]
                        common_mass += rB[uB] / n
                        resA[uA] = (rA[uA] - rB[uB]) / n
                    else:
                        assert len(adjA[uA]) <= len(adjB[uB])
                        assert rA[uA] <= rB[uB]
                        common_mass += rA[uA] / n
                        resB[uB] = (rB[uB] - rA[uA]) / n

                for uA in range(n):
                    if uA not in matchedA:
                        resA[uA] = rA[uA] / n
                for uB in range(n):
                    if uB not in matchedB:
                        resB[uB] = rB[uB] / n

                rhoA = sum(resA)
                rhoB = sum(resB)
                assert rhoA == rhoB == 1 - common_mass and rhoA >= 0

                # Check exact marginals of the one-step coupling
                margA = [F(0)] * n
                margB = [F(0)] * n
                for uA, uB, par in pairs:
                    p = min(rA[uA], rB[uB]) / n
                    margA[uA] += p
                    margB[uB] += p
                    # Verify invariance after common step
                    if p > 0:
                        nA = [list(x) for x in adjA] + [[uA]]
                        nA[uA].append(n)
                        nB = [list(x) for x in adjB] + [[uB]]
                        nB[uB].append(n)
                        assert ge_even(adj_to_rooted(nA, 0), adj_to_rooted(nB, 0))
                        coupled_outcomes_checked += 1
                if rhoA > 0:
                    for uA in range(n):
                        if resA[uA] == 0:
                            continue
                        for uB in range(n):
                            if resB[uB] == 0:
                                continue
                            p = resA[uA] * resB[uB] / rhoA
                            margA[uA] += p
                            margB[uB] += p
                            # Verify invariance after residual step
                            nA = [list(x) for x in adjA] + [[uA]]
                            nA[uA].append(n)
                            nB = [list(x) for x in adjB] + [[uB]]
                            nB[uB].append(n)
                            assert ge_even(adj_to_rooted(nA, 0), adj_to_rooted(nB, 0))
                            coupled_outcomes_checked += 1
                assert margA == [x / n for x in rA]
                assert margB == [x / n for x in rB]
                coupled_pairs_checked += 1

    print(
        f"[PART A] Exact one-step coupling verified on ALL {coupled_pairs_checked} comparable pairs "
        f"(n=3..6), {coupled_outcomes_checked} coupled transitions: degrees, rates, marginals, and invariance exact"
    )

    # Planted defect test: flip parity requirement (same-parity injection fails rate monotonicity!)
    # Take 4-path rooted at interior vs 4-path rooted at leaf -> they satisfy ge_even,
    # while swapping them fails ge_even and breaks rate monotonicity.
    p4_center = ((), ((),))
    p4_leaf = ((((),),),)
    assert ge_even(p4_center, p4_leaf)
    assert not ge_even(p4_leaf, p4_center)
    print("[PART A] Planted direction reversal caught cleanly (p4_center >= p4_leaf, reverse false)")

    # Part (b): Corollary 1 on ALL trees of sizes 3..7
    cor1_checks = 0
    for n in range(3, 8):
        for t in by_n[n]:
            adj, _ = rooted_to_adj(t)
            leaves = [v for v in range(n) if len(adj[v]) == 1]
            has_leaf_nb = [v for v in range(n) if any(len(adj[w]) == 1 for w in adj[v])]
            for u in has_leaf_nb:
                ru = adj_to_rooted(adj, u)
                for w in leaves:
                    rw = adj_to_rooted(adj, w)
                    assert ge_even(ru, rw)
                    cor1_checks += 1
    print(f"[PART B] Corollary 1 (vertex with leaf neighbour >= every leaf) verified on {cor1_checks} (tree,u,w) triples (n=3..7)")


# ---------------------------------------------------------------------------
# 3. Independent Exact Enumeration + Exact Max-Flow Transport (Part c)
# ---------------------------------------------------------------------------
def exact_maxflow(supply_dict, demand_dict, rel_fn):
    """
    Exact integer Dinic max-flow between two Fraction distributions using common denominator.
    """
    S = [k for k, v in supply_dict.items() if v > 0]
    T = [k for k, v in demand_dict.items() if v > 0]
    den = 1
    for v in list(supply_dict.values()) + list(demand_dict.values()):
        den = math.lcm(den, v.denominator)

    ns, nt = len(S), len(T)
    src, snk = ns + nt, ns + nt + 1
    N = snk + 1
    adj = [[] for _ in range(N)]
    to = []
    cap = []

    def add_edge(u, v, c):
        adj[u].append(len(to))
        to.append(v)
        cap.append(c)
        adj[v].append(len(to))
        to.append(u)
        cap.append(0)

    total_sup = 0
    for i, s in enumerate(S):
        c = int(supply_dict[s] * den)
        total_sup += c
        add_edge(src, i, c)
    for j, t in enumerate(T):
        c = int(demand_dict[t] * den)
        add_edge(ns + j, snk, c)

    INF = total_sup + 1
    arcs = 0
    for i, s in enumerate(S):
        for j, t in enumerate(T):
            if rel_fn(s, t):
                add_edge(i, ns + j, INF)
                arcs += 1

    flow = 0
    while True:
        lvl = [-1] * N
        lvl[src] = 0
        q = deque([src])
        while q:
            u = q.popleft()
            for e in adj[u]:
                if cap[e] > 0 and lvl[to[e]] < 0:
                    lvl[to[e]] = lvl[u] + 1
                    q.append(to[e])
        if lvl[snk] < 0:
            break
        it = [0] * N

        def dfs(u, pushed):
            if u == snk:
                return pushed
            while it[u] < len(adj[u]):
                e = adj[u][it[u]]
                v = to[e]
                if cap[e] > 0 and lvl[v] == lvl[u] + 1:
                    tr = dfs(v, min(pushed, cap[e]))
                    if tr > 0:
                        cap[e] -= tr
                        cap[e ^ 1] += tr
                        return tr
                it[u] += 1
            return 0

        while True:
            p = dfs(src, INF)
            if p == 0:
                break
            flow += p

    return F(flow, den), arcs


def check_part_c():
    # Enumerate exact labeled histories by aggregating identical unrooted labeled?
    # Better: aggregate identical rooted/unrooted distributions step-by-step independently
    # using tuple-based canonical unrooted codes (center-based, completely independent of Fable's bitsets!)
    def unrooted_canon(adj):
        n = len(adj)
        if n == 1:
            return ((),)
        deg = [len(x) for x in adj]
        q = deque([v for v in range(n) if deg[v] <= 1])
        rem = n
        while rem > 2:
            k = len(q)
            rem -= k
            for _ in range(k):
                u = q.popleft()
                for w in adj[u]:
                    deg[w] -= 1
                    if deg[w] == 1:
                        q.append(w)
        centers = list(q)
        if len(centers) == 1:
            return (adj_to_rooted(adj, centers[0]),)
        c1, c2 = centers
        r1 = adj_to_rooted(adj, c1, c2)
        r2 = adj_to_rooted(adj, c2, c1)
        return tuple(sorted([r1, r2]))

    def unrooted_to_adj(uc):
        if len(uc) == 1:
            adj, _ = rooted_to_adj(uc[0])
            return adj
        adj1, _ = rooted_to_adj(uc[0])
        adj2, _ = rooted_to_adj(uc[1])
        off = len(adj1)
        adj = [list(x) for x in adj1] + [[w + off for w in x] for x in adj2]
        adj[0].append(off)
        adj[off].append(0)
        return adj

    # First verify our unrooted propagator against brute-force labeled history enumeration through N = 8!
    brute = [(F(1), ((1,), (0,)))]
    for n in range(2, 8):
        nxt = []
        for p, adj in brute:
            r = attachment_rates(adj)
            for w in range(n):
                nadj = [list(x) for x in adj] + [[w]]
                nadj[w].append(n)
                nxt.append((p * r[w] / n, tuple(tuple(x) for x in nadj)))
        brute = nxt
    assert len(brute) == math.factorial(7)
    assert sum(p for p, _ in brute) == 1

    brute_q = [dict() for _ in range(9)]  # 1-indexed labels at N=8
    for p, adj in brute:
        for u in range(1, 9):
            s = adj_to_rooted(adj, u - 1)
            brute_q[u][s] = brute_q[u].get(s, F(0)) + p

    # Check all pairs (u, u+1) at N = 8 from brute force
    assert brute_q[1] == brute_q[2]
    for u in range(1, 8):
        fl, arcs = exact_maxflow(brute_q[u], brute_q[u + 1], ge_even)
        assert fl == 1, f"Brute force N=8 pair ({u},{u+1}) failed: {fl}"
    print("[PART C] From-scratch brute-force 7! = 5,040 histories at N=8: Z_1=Z_2 and ALL 7 pairs (u,u+1) max-flow = 1 EXACT")

    # Now propagate unrooted laws up to N = 12 independently and compare EXACT numerators
    # against Fable's child-transport-all.txt & lead-transport_sparse_2_11.log!
    fable_expected = {
        (2, 3): (0, 1),
        (3, 4): (0, 2),
        (4, 5): (0, 28),
        (5, 6): (0, 1152),
        (6, 7): (6480, 309264),
        (7, 8): (2890560, 102063360),
        (8, 9): (10495556400, 279352830000),
        (9, 10): (76217022504000, 1766482425204000),
        (10, 11): (1815415088183726400, 37997911398273321600),
        (11, 12): (46506935955629227276800, 913672344783897871401600),
    }

    U = {2: {unrooted_canon([[1], [0]]): F(1)}}
    den = 1
    for n in range(2, 12):
        # Compute U[n+1] and newborn Q_{n+1} at size n+1
        nxt_U = {}
        Q_new = {}
        for uc, prob in U[n].items():
            adj = unrooted_to_adj(uc)
            r = attachment_rates(adj)
            for w in range(n):
                p_step = prob * r[w] / n
                nadj = [list(x) for x in adj] + [[w]]
                nadj[w].append(n)
                uc2 = unrooted_canon(nadj)
                nxt_U[uc2] = nxt_U.get(uc2, F(0)) + p_step
                rnew = adj_to_rooted(nadj, n)
                Q_new[rnew] = Q_new.get(rnew, F(0)) + p_step
        assert sum(nxt_U.values()) == 1 and sum(Q_new.values()) == 1
        U[n + 1] = nxt_U

        # Compute older label Q_u at size N = n+1 (where u = n)
        u = n
        N = n + 1
        den *= n * math.lcm(*range(1, n)) if n >= 2 else 1
        if u == 2:
            Q_old = {
                adj_to_rooted([[1, 2], [0], [0]], 0): F(1, 2),
                adj_to_rooted([[1, 2], [0], [0]], 1): F(1, 2),
            }
        else:
            Q_old = {}
            for rt, prob in Q_prev_new.items():
                adj, _ = rooted_to_adj(rt)
                r = attachment_rates(adj)
                for w in range(u):
                    nadj = [list(x) for x in adj] + [[w]]
                    nadj[w].append(u)
                    rt2 = adj_to_rooted(nadj, 0)
                    Q_old[rt2] = Q_old.get(rt2, F(0)) + prob * r[w] / u
        assert sum(Q_old.values()) == 1
        if N == 8:
            assert Q_old == brute_q[7] and Q_new == brute_q[8], "Unrooted propagator mismatch with brute force at N=8!"

        # Signed difference d = Q_old - Q_new
        keys = set(Q_old) | set(Q_new)
        leaf_pos = {}
        leaf_neg = {}
        univ = F(0)
        for k in keys:
            d = Q_old.get(k, F(0)) - Q_new.get(k, F(0))
            if d == 0:
                continue
            if len(k) >= 2:
                # Must be positive and have a leaf child ()
                assert d > 0 and () in k
                univ += d
            elif d > 0:
                leaf_pos[k] = d
            else:
                leaf_neg[k] = -d
        assert sum(leaf_pos.values()) + univ == sum(leaf_neg.values())

        # Compare exact integer numerators on Fable's denominator Den[N]
        want_num = sum(leaf_pos.values()) * den
        univ_num = univ * den
        assert want_num.denominator == 1 and univ_num.denominator == 1
        assert (want_num.numerator, univ_num.numerator) == fable_expected[(u, u + 1)], (
            f"Mismatch at ({u},{u+1}): got {(want_num.numerator, univ_num.numerator)} vs {fable_expected[(u,u+1)]}"
        )

        # Solve direct exact bipartite max-flow for N <= 10 (fast in pure Python!)
        if N <= 10:
            fl, _ = exact_maxflow(leaf_pos, leaf_neg, ge_even)
            assert fl == sum(leaf_pos.values()), f"Leaf transport failed at ({u},{u+1})"
            # Planted defect: testing the FALSE reverse domination Z_{u+1} >= Z_u must fail (< 1)
            fl_rev, _ = exact_maxflow(Q_new, Q_old, ge_even)
            assert fl_rev < 1, f"False reverse domination Z_{u+1}>=Z_u unexpectedly passed at ({u},{u+1})!"

        Q_prev_new = Q_new

    print(
        "[PART C] Independent exact unrooted-law propagator matches ALL 10 pairs (2,3)..(11,12) "
        "exact integer numerators (leaf_pos and universal) in Fable's logs, and independent direct "
        "max-flow passes (2,3)..(9,10) (with reversed-order planted defect rejected for (6,7)..(9,10))"
    )

    # Check the 5-vertex conditional obstruction from child-REPORT.md §2 / COUPLING.md §6
    # Seed: root 0 with 3 children: two leaves (1,2) and one child (3) with leaf (4)
    seed_adj = [[1, 2, 3], [0], [0], [0, 4], [3]]
    q6 = {}
    q7_old = {}
    q7_new = {}
    r5 = attachment_rates(seed_adj)
    for w1 in range(5):
        p1 = r5[w1] / 5
        a6 = [list(x) for x in seed_adj] + [[w1]]
        a6[w1].append(5)
        r6 = attachment_rates(a6)
        for w2 in range(6):
            p12 = p1 * r6[w2] / 6
            a7 = [list(x) for x in a6] + [[w2]]
            a7[w2].append(6)
            s_old = adj_to_rooted(a7, 5)
            s_new = adj_to_rooted(a7, 6)
            q7_old[s_old] = q7_old.get(s_old, F(0)) + p12
            q7_new[s_new] = q7_new.get(s_new, F(0)) + p12
    fl_cond, _ = exact_maxflow(q7_old, q7_new, ge_even)
    assert fl_cond == F(2159, 2160), f"Unexpected full conditional flow: {fl_cond}"
    # Also check cancelled residual flow (matches Fable's 1210/1220 on common denominator 21600)
    pos7, neg7 = {}, {}
    for k in set(q7_old) | set(q7_new):
        d = q7_old.get(k, F(0)) - q7_new.get(k, F(0))
        if d > 0:
            pos7[k] = d
        elif d < 0:
            neg7[k] = -d
    fl_res, _ = exact_maxflow(pos7, neg7, ge_even)
    tot_pos = sum(pos7.values())
    assert tot_pos == F(1220, 21600) and fl_res == F(1210, 21600), (tot_pos, fl_res)
    assert (tot_pos - fl_res) == F(1, 2160)
    print(
        f"[PART C] Conditional 5-vertex seed obstruction independently reproduced: "
        f"full flow = {fl_cond} (deficit 1/2160), residual flow = {int(fl_res*21600)}/{int(tot_pos*21600)} on denom 21600"
    )


# ---------------------------------------------------------------------------
# 4. Check Part (d) (Strictness) and Part (e) (Beta(1/2, v-1/2) & exact tail)
# ---------------------------------------------------------------------------
def check_parts_de():
    by_n = all_rooted_trees(7)

    # (d1) Verify drift bounds for D/n (supermartingale) and (L-1)/n (submartingale) on all rooted trees n=3..7
    drift_checks = 0
    for n in range(3, 8):
        for t in by_n[n]:
            adj, _ = rooted_to_adj(t)
            r = attachment_rates(adj)
            D = len(adj[0])
            L = sum(1 for w in adj[0] if len(adj[w]) == 1)
            # Exact expected next D/(n+1) and (L-1)/(n+1)
            e_next_D = F(0)
            e_next_L = F(0)
            for w in range(n):
                pw = r[w] / n
                nD = D + (1 if w == 0 else 0)
                if w == 0:
                    nL = L + 1
                elif w in adj[0] and len(adj[w]) == 1:
                    nL = L - 1
                else:
                    nL = L
                e_next_D += pw * F(nD, n + 1)
                e_next_L += pw * F(nL - 1, n + 1)
            assert e_next_D <= F(D, n)
            assert e_next_L >= F(L - 1, n)
            drift_checks += 1
    print(f"[PART D] Exact supermartingale D/n and submartingale (L-1)/n verified on all {drift_checks} rooted trees (n=3..7)")

    # (d2) Verify promotion supermartingale f(k) = 1/max(k-1, 1) for all 3 <= k < d <= 60
    promo_checks = 0
    for d in range(4, 61):
        for k in range(3, d):
            # Exact worst-case drift with birth rate r(p) >= k+1 > k (even with r(p) = k)
            drift_k = F(k, 1) * (F(1, k) - F(1, k - 1)) + F(k, d) * (F(1, k - 2) - F(1, k - 1)) + F(1, d) * (1 - F(1, k - 1))
            assert drift_k <= 0
            # Verify algebraic identity in child-REPORT.md line 23
            rhs = (-1 + (F(k, k - 2) + (k - 2)) / d) / (k - 1)
            assert drift_k == rhs and rhs <= 0
            promo_checks += 1
    print(f"[PART D] Promotion supermartingale f(k)=1/max(k-1,1) exact drift <= 0 verified on all {promo_checks} (d,k) pairs")

    # (d3) Verify strict residual positivity on universal pairs (s, t)
    res_checks = 0
    for n in range(3, 7):
        univ_s = [t for t in by_n[n] if len(t) >= 2 and () in t]
        leaf_t = [t for t in by_n[n] if len(t) == 1]
        for s in univ_s:
            adjA, _ = rooted_to_adj(s)
            rA = attachment_rates(adjA)
            # pick a leaf child a1 of root 0 in s
            a1 = next(w for w in adjA[0] if len(adjA[w]) == 1)
            for t in leaf_t:
                adjB, _ = rooted_to_adj(t)
                rB = attachment_rates(adjB)
                b1 = adjB[0][0]
                assert rA[0] > rB[0], "Root A residual not positive!"
                assert rB[b1] > rA[a1], "Parent-of-B residual not positive!"
                res_checks += 1
    print(f"[PART D] Strict positivity of both forced residual targets verified on all {res_checks} universal pairs (n=3..6)")

    # (e) Verify r_T(v) <= D_v - 1/2 for every vertex v in a tree of size >= 3 with at least one neighbour of degree >= 2
    polya_checks = 0
    for n in range(3, 8):
        for t in by_n[n]:
            adj, _ = rooted_to_adj(t)
            r = attachment_rates(adj)
            for v in range(n):
                if any(len(adj[w]) >= 2 for w in adj[v]):
                    assert r[v] <= len(adj[v]) - F(1, 2)
                    polya_checks += 1
    # Check lower-bound product identity prod_{j=v}^{m-1} (j - v + 1/2)/j = Gamma(v)/sqrt(pi) * Gamma(m-v+1/2)/Gamma(m)
    for v in range(3, 12):
        m = 500
        prod = 1.0
        for j in range(v, m):
            prod *= (j - v + 0.5) / j
        exact_gamma = (math.gamma(v) / math.sqrt(math.pi)) * math.exp(
            math.lgamma(m - v + 0.5) - math.lgamma(m)
        )
        assert abs(prod - exact_gamma) / exact_gamma < 1e-12
    print(f"[PART E] Polya rate bound r(v) <= D_v - 1/2 verified on {polya_checks} vertex instances (n=3..7), and lower-bound Gamma product matches to <1e-12 for v=3..11")


if __name__ == "__main__":
    check_parts_ab()
    check_part_c()
    check_parts_de()
    print("ALL INDEPENDENT GEMINI-FAMILY AUDITS OF FABLE CLAIM 5 PASSED.")
