#!/usr/bin/env python3
"""Independent exact-rational and numerical verification for Neon-2 Claim 7 (density-lp)."""

from fractions import Fraction
import math
import random

# Certified exact rational bracket for mu from Neon-2 Claim 6 (commit b04cc50)
MU_LO = Fraction(2825181729593748861, 4994114995300000000)
MU_HI = Fraction(439124781075521909, 776183857175000000)


def check_exponents():
    """Verify exact rational inequalities for gamma, beta, xi_*, and 1/(1 - xi_*)."""
    # gamma = 1 + 2/(2*mu - 1) is decreasing in mu
    gamma_hi = 1 + Fraction(2, 1) / (2 * MU_LO - 1)
    gamma_lo = 1 + Fraction(2, 1) / (2 * MU_HI - 1)
    assert gamma_hi < Fraction(33, 2), f"gamma_hi={float(gamma_hi)}"
    assert gamma_lo > 16, f"gamma_lo={float(gamma_lo)}"

    # beta = 1 + 4/(2*mu - 1) is decreasing in mu
    beta_hi = 1 + Fraction(4, 1) / (2 * MU_LO - 1)
    beta_lo = 1 + Fraction(4, 1) / (2 * MU_HI - 1)

    # xi_* = (2*mu - 1) / (2*mu^2 + 5*mu - 1)
    # Derivative wrt mu: [2(2mu^2+5mu-1) - (2mu-1)(4mu+5)] / denom^2
    # = [4mu^2+10mu-2 - (8mu^2+6mu-5)] / denom^2 = (-4mu^2+4mu+3)/denom^2 > 0 on (0, 1)!
    for mu in (MU_LO, MU_HI):
        deriv_num = -4 * mu * mu + 4 * mu + 3
        assert deriv_num > 0
    xi_lo = (2 * MU_LO - 1) / (2 * MU_LO * MU_LO + 5 * MU_LO - 1)
    xi_hi = (2 * MU_HI - 1) / (2 * MU_HI * MU_HI + 5 * MU_HI - 1)
    assert xi_lo > Fraction(1, 19), f"xi_lo={float(xi_lo)}"
    p_lo = Fraction(1, 1) / (1 - xi_lo)
    assert p_lo > Fraction(19, 18), f"p_lo={float(p_lo)}"
    # Verify ratio vs Astra's 1/80000
    assert xi_lo * 80000 > 4200

    # Verify exact Corollary 2 balancing identity:
    # (beta + 1/mu) * xi_* == 1/mu, and (1 - beta/(beta + 1/mu)) == xi_*
    for mu in (MU_LO, MU_HI):
        beta = 1 + Fraction(4, 1) / (2 * mu - 1)
        d = Fraction(1, 1) / mu
        xi = (2 * mu - 1) / (2 * mu * mu + 5 * mu - 1)
        assert xi == d / (beta + d)
        assert 1 - beta / (beta + d) == xi

    # Verify non-monotonicity of F(delta) = delta * (log(1/delta))^gamma near delta = 1/2
    gam = float(gamma_lo)
    f_50 = 0.50 * (math.log(1.0 / 0.50) ** gam)
    f_49 = 0.49 * (math.log(1.0 / 0.49) ** gam)
    assert f_50 < f_49 < 1.0, "F(delta) must be decreasing near 0.5 (motivating Corollary 2 padding)"

    # Verify exact 7-stage planting factors in Lemma 2:
    # exp powers: Stage1(2) + Stage2(1) + Stage3(2) + Stage4(2) + Stage5(3) + Stage6(3) + Stage7(4) = 17
    # denominators: Stage1(2) * Stage2(6) * Stage4(2) * Stage6(2) = 48
    exp_sum = 2 + 1 + 2 + 2 + 3 + 3 + 4
    denom_prod = 2 * 6 * 2 * 2
    assert exp_sum == 17 and denom_prod == 48
    # Planted defect control: old e^-15/24 exceeds true Type-2 e^-17/48 by 2*e^2
    # (and printed (1 - e^-1)*e^-15/24 exceeds exact Type-1 e^-16/24 by factor e - 1 > 1.7)
    assert math.exp(-15) / 24 > math.exp(-17) / 48 * 14.0
    assert (1.0 - math.exp(-1)) * math.exp(-15) / 24 > math.exp(-16) / 24 * 1.7

    print(
        f"Exponents & planting OK: gamma in [{float(gamma_lo):.6f}, {float(gamma_hi):.6f}] < 33/2; "
        f"beta in [{float(beta_lo):.6f}, {float(beta_hi):.6f}]; "
        f"xi_* in [{float(xi_lo):.6f}, {float(xi_hi):.6f}] > 1/19 ({float(xi_lo*80000):.1f}x Astra); "
        f"1/(1-xi_*) >= {float(p_lo):.6f} > 19/18; 7-stage p_* factor e^-{exp_sum}/{denom_prod}; "
        f"F(0.50)={f_50:.6f} < F(0.49)={f_49:.6f}."
    )


def check_supermartingale_exact():
    """Exact rational check of Lemma 3 supermartingale (Gh)(l) < 0 for l = 5..2000 + planted defects."""
    def gh(l, c_num):
        h_l = min(Fraction(1, 1), Fraction(c_num, l))
        h_p = min(Fraction(1, 1), Fraction(c_num, l + 1))
        h_m = min(Fraction(1, 1), Fraction(c_num, l - 1))
        alarm_rate = Fraction(2, l + 1)
        return l * (h_p - h_l) + alarm_rate * (1 - h_l) + (1 - alarm_rate) * (h_m - h_l)

    max_val = None
    for l in range(5, 2001):
        val = gh(l, 4)
        # Verify exact algebraic formula in README
        formula = (
            -Fraction(2, l + 1)
            + Fraction(4, l * (l - 1))
            - Fraction(8, l * (l + 1))
            - Fraction(8, l * (l + 1) * (l - 1))
        )
        assert val == formula, f"formula mismatch at l={l}"
        assert val == -Fraction(2 * (l + 2), l * (l + 1)), f"simplified formula mismatch at l={l}"
        assert val < 0, f"(Gh)({l}) = {val} >= 0"
        if max_val is None or val > max_val:
            max_val = val

    # Check K = 24 barrier bound: 3 * (4/24) == 1/2
    assert 3 * Fraction(4, 24) == Fraction(1, 2)

    # Planted defect 1: c_num = 1 fails for all l >= 2 (leading drift +(2 - c)/(l+1) = +1/(l+1) > 0)
    assert gh(10, 1) > 0, "planted defect c=1 must fail!"
    # Planted defect 2: uncoupled bound (dropping negative -2c/(l(l+1)) coupling) with c=3 fails at l=4
    uncoupled_c3_l4 = -Fraction(3 - 2, 4 + 1) + Fraction(3, 4 * 3)
    assert uncoupled_c3_l4 > 0, "planted defect uncoupled c=3 at l=4 must fail!"
    # Planted defect 3: c_num = 4 with alarm rate 5/(l+1) fails (leading drift +1/(l+1) > 0)
    l_t = 20
    h_l, h_p, h_m = Fraction(4, l_t), Fraction(4, l_t + 1), Fraction(4, l_t - 1)
    bad_alarm = l_t * (h_p - h_l) + Fraction(5, l_t + 1) * (1 - h_l) + (1 - Fraction(5, l_t + 1)) * (h_m - h_l)
    assert bad_alarm > 0, "planted defect alarm rate 5/(l+1) must fail!"
    assert Fraction(2, 24 + 2) < Fraction(2, 24 + 1)

    print(
        f"Supermartingale exact OK: (Gh)(l) < 0 for all l=5..2000 "
        f"(at l=5: {gh(5,4)}, max on 5..2000: {float(max_val):.6e}); "
        f"3 planted defects caught (c=1, uncoupled c=3 at l=4, alarm 5/(l+1))."
    )


def check_azuma_gamma_bounds():
    """Verify a_{j+1} >= (j+mu)^mu, a_n > n^mu / 2, and tail sum <= 2/(2*mu-1) * n_0^(1-2*mu)."""
    for mu in (float(MU_LO), float(MU_HI)):
        for j in range(1, 500):
            a_jp1 = math.exp(math.lgamma(j + 1 + mu) - math.lgamma(j + 1))
            assert a_jp1 >= (j + mu) ** mu - 1e-12
            assert a_jp1 <= (j + 1) ** mu + 1e-12
        for n0 in range(3, 100):
            a_n0 = math.exp(math.lgamma(n0 + mu) - math.lgamma(n0))
            assert a_n0 > 0.5 * (n0 ** mu)
            # Compute sum_{j=n0}^{20000} a_{j+1}^{-2} + exact integral tail from 20000
            partial = sum(
                math.exp(-2.0 * (math.lgamma(j + 1 + mu) - math.lgamma(j + 1)))
                for j in range(n0, 20000)
            )
            tail = (20000 - 1 + mu) ** (1.0 - 2.0 * mu) / (2.0 * mu - 1.0)
            rhs = 2.0 / (2.0 * mu - 1.0) * (n0 ** (1.0 - 2.0 * mu))
            assert partial + tail <= rhs
    print("Azuma/Gamma bounds OK: verified Gautschi lower/upper bounds and tail sums for n0=3..99.")


def pack_disjoint_gadgets(n_cores, edges, leaves_per_core):
    """Implement Lemma 2(1) disjoint packing on a core tree and verify exact disjointness."""
    adj = {u: [] for u in range(n_cores)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    assert n_cores >= 1
    if n_cores == 1:
        assert leaves_per_core[0] >= 2
        return [((0,), "T1")]

    c_deg = {u: len(adj[u]) for u in range(n_cores)}
    for u in range(n_cores):
        assert c_deg[u] >= 1
        if c_deg[u] == 1:
            assert leaves_per_core[u] >= 1

    S1 = {u for u in range(n_cores) if c_deg[u] <= 2 and leaves_per_core[u] >= 1}
    Sge3 = {u for u in range(n_cores) if c_deg[u] >= 3}
    M = S1 | Sge3

    blocks = [((u,), "T1") for u in sorted(S1)]
    visited = set(M)
    bare_paths = 0
    for u in range(n_cores):
        if u in visited:
            continue
        # BFS/DFS component of K \ M
        comp = []
        stack = [u]
        visited.add(u)
        while stack:
            x = stack.pop()
            comp.append(x)
            for nb in adj[x]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        bare_paths += 1
        # Order comp as a path
        comp_set = set(comp)
        if len(comp) == 1:
            path = comp
        else:
            ends = [x for x in comp if sum(1 for nb in adj[x] if nb in comp_set) == 1]
            assert len(ends) == 2
            path = [ends[0]]
            prev = None
            cur = ends[0]
            while True:
                nxts = [nb for nb in adj[cur] if nb in comp_set and nb != prev]
                if not nxts:
                    break
                prev, cur = cur, nxts[0]
                path.append(cur)
            assert len(path) == len(comp)
        for i in range(0, len(path) - 2, 3):
            u0, w0, x0 = path[i], path[i + 1], path[i + 2]
            # Verify Type-2 invariants strictly
            assert c_deg[w0] == 2 and leaves_per_core[w0] == 0
            assert set(adj[w0]) == {u0, x0}
            assert c_deg[u0] == 2 and leaves_per_core[u0] == 0
            blocks.append(((u0, w0, x0), "T2"))

    # Verify pairwise vertex-disjointness
    seen = set()
    for verts, kind in blocks:
        for v in verts:
            assert v not in seen, "OVERLAPPING BLOCKS!"
            seen.add(v)
    assert bare_paths <= len(M) - 1
    assert len(blocks) * 6 >= n_cores + 8
    return blocks


def all_labeled_trees(n):
    """Generate all labeled trees on 0..n-1 via Prufer sequences (for small n)."""
    if n == 1:
        yield []
        return
    if n == 2:
        yield [(0, 1)]
        return
    from itertools import product
    for seq in product(range(n), repeat=n - 2):
        deg = [1] * n
        for v in seq:
            deg[v] += 1
        edges = []
        for v in seq:
            for u in range(n):
                if deg[u] == 1:
                    edges.append((u, v))
                    deg[u] -= 1
                    deg[v] -= 1
                    break
        rem = [u for u in range(n) if deg[u] == 1]
        edges.append((rem[0], rem[1]))
        yield edges


def check_combinatorial_packing():
    """Verify Lemma 2(1) exhaustively on C=1 and labeled C=2..7 + extremal binary & caterpillar trees + random trees."""
    checked = 0
    # 0. C = 1 base case
    b1 = pack_disjoint_gadgets(1, [], [2])
    assert len(b1) == 1 and 6 * len(b1) >= 1
    checked += 1

    # 1. Exhaustive labeled core trees for C = 2..7 with minimal leaves (0 on non-leaves of core tree!)
    for C in range(2, 8):
        for edges in all_labeled_trees(C):
            deg = [0] * C
            for u, v in edges:
                deg[u] += 1
                deg[v] += 1
            # Worst case: ONLY c=1 cores have a leaf (L=1); all c>=2 cores have L=0
            leaves = [1 if deg[u] == 1 else 0 for u in range(C)]
            pack_disjoint_gadgets(C, edges, leaves)
            checked += 1

    # 2. Extremal sharpness check: full binary trees with n1 leaves, every edge subdivided by 2 bare vertices
    for depth in range(1, 7):
        # Build binary tree with 2^depth leaves -> 2^(depth+1)-1 vertices, 2^(depth+1)-2 edges
        # Wait: root has degree 2, make root degree 3 by attaching 3 binary subtrees!
        edges = []
        nxt = 1
        # Build rooted binary tree of given depth
        def add_sub(parent, d):
            nonlocal nxt
            cur = nxt
            nxt += 1
            # Subdivide edge (parent, cur) by 2 bare vertices
            b1 = nxt; nxt += 1
            b2 = nxt; nxt += 1
            edges.extend([(parent, b1), (b1, b2), (b2, cur)])
            if d > 0:
                add_sub(cur, d - 1)
                add_sub(cur, d - 1)
        for _ in range(3):
            add_sub(0, depth - 1)
        C = nxt
        deg = [0] * C
        for u, v in edges:
            deg[u] += 1
            deg[v] += 1
        leaves = [1 if deg[u] == 1 else 0 for u in range(C)]
        blocks = pack_disjoint_gadgets(C, edges, leaves)
        assert 6 * len(blocks) == C + 8, f"Expected exact sharpness 6R == C+8, got 6*{len(blocks)} vs {C+8}"
        checked += 1

    # 2b. Extremal subdivided cubic caterpillars (spine of k degree-3 vertices, each edge subdivided twice)
    for k in range(1, 7):
        # k spine vertices 0..k-1, each adjacent to leaf-end(s) so every spine vertex has core degree 3
        edges = []
        nxt = k
        def add_sub_edge(u, v):
            nonlocal nxt
            b1 = nxt; nxt += 1
            b2 = nxt; nxt += 1
            edges.extend([(u, b1), (b1, b2), (b2, v)])
        for i in range(k - 1):
            add_sub_edge(i, i + 1)
        for i in range(k):
            needed = 3 - (0 if k == 1 else (1 if i in (0, k - 1) else 2))
            for _ in range(needed):
                end = nxt; nxt += 1
                add_sub_edge(i, end)
        C = nxt
        deg = [0] * C
        for u, v in edges:
            deg[u] += 1
            deg[v] += 1
        leaves = [1 if deg[u] == 1 else 0 for u in range(C)]
        blocks = pack_disjoint_gadgets(C, edges, leaves)
        assert 6 * len(blocks) == C + 8, f"Expected exact sharpness 6R == C+8 on cubic caterpillar, got 6*{len(blocks)} vs {C+8}"
        checked += 1

    # 3. Random trees + long paths + caterpillars up to C=500
    rng = random.Random(20261124)
    for C in (10, 50, 200, 500):
        # Long path
        edges = [(i, i + 1) for i in range(C - 1)]
        leaves = [1 if i in (0, C - 1) else 0 for i in range(C)]
        blocks = pack_disjoint_gadgets(C, edges, leaves)
        # Planted defect control: Type-1 only on long path has only 2 blocks << C/6!
        t1_only = sum(1 for _, k in blocks if k == "T1")
        if C >= 50:
            assert t1_only == 2 and t1_only < C / 6, "planted defect control: Type-1 only must fail on path"
        checked += 1
        # Random Prufer trees
        for _ in range(5):
            seq = [rng.randrange(C) for _ in range(C - 2)]
            deg = [1] * C
            for v in seq:
                deg[v] += 1
            import heapq
            h = [u for u in range(C) if deg[u] == 1]
            heapq.heapify(h)
            edges = []
            for v in seq:
                u = heapq.heappop(h)
                edges.append((u, v))
                deg[u] -= 1
                deg[v] -= 1
                if deg[v] == 1:
                    heapq.heappush(h, v)
            edges.append((heapq.heappop(h), heapq.heappop(h)))
            deg_f = [0] * C
            for u, v in edges:
                deg_f[u] += 1
                deg_f[v] += 1
            leaves = [1 if deg_f[u] == 1 else (1 if rng.random() < 0.2 else 0) for u in range(C)]
            pack_disjoint_gadgets(C, edges, leaves)
            checked += 1

    print(
        f"Combinatorial packing OK: verified R >= (C+8)/6 with 100% disjointness on {checked} trees "
        f"(exhaustive labeled C=2..7, exact sharpness 6R=C+8 on 6 subdivided binary trees, "
        f"paths/random trees up to C=500; Type-1-only planted defect caught)."
    )


def check_polar_fubini_and_dirichlet():
    """Verify Lemma 4: Polar-Fubini bound -psi'' >= mu(1-mu)/16 and integral bound <= 32*pi/(mu(1-mu)M)."""
    mu = float(MU_LO)
    p = 1.0 / (1.0 - mu)
    bound_d2 = mu * (1.0 - mu) / 16.0
    bound_area = 32.0 * math.pi / (mu * (1.0 - mu))

    triples = [
        (1.0, 1.0, 1.0),
        (1.0, 0.5, 0.1),
        (1.0, 1.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1e-4, 0.0),
        (1.0, 1e-8, 0.0),
        (1.0, 1e-4, 1e-8),
        (1.0, 1e-8, 1e-16),
        (1.0, 0.999999, 1e-12),
    ]
    min_d2 = float("inf")
    max_slice_density = 0.0

    n_theta = 720
    for x1, x2, x3 in triples:
        S = 1.0 + (x2 ** p) + (x3 ** p)
        u1s, u2s, u3s = 1.0 / S, (x2 ** p) / S, (x3 ** p) / S
        g_max = S ** (1.0 - mu)
        # Verify Holder maximum on random points
        for a, b in ((0.3, 0.3), (0.01, 0.01), (0.98, 0.01), (0.01, 0.98)):
            c = 1.0 - a - b
            assert x1 * (a ** mu) + x2 * (b ** mu) + x3 * (c ** mu) <= g_max + 1e-12

        # Check all rays theta
        for k in range(n_theta):
            th = (2.0 * math.pi * (k + 0.37)) / n_theta
            v1 = math.cos(th)
            v2 = math.sin(th)
            v3 = -v1 - v2
            # Compute T_theta = max t > 0 with u* + t v in Delta_2
            t_max = float("inf")
            enters = True
            for ui, vi in ((u1s, v1), (u2s, v2), (u3s, v3)):
                if ui == 0.0 and vi <= 0.0:
                    enters = False
                    break
                if vi < 0.0:
                    t_max = min(t_max, -ui / vi)
            if not enters or t_max <= 1e-15:
                continue
            for frac in (1e-4, 0.01, 0.1, 0.5, 0.9, 0.999):
                s = frac * t_max
                u1, u2, u3 = u1s + s * v1, u2s + s * v2, u3s + s * v3
                assert u1 > 0 and u2 > 0 and u3 > 0
                neg_d2 = mu * (1.0 - mu) * (
                    x1 * (u1 ** (mu - 2.0)) * v1 * v1
                    + (x2 * (u2 ** (mu - 2.0)) * v2 * v2 if x2 > 0 else 0.0)
                    + (x3 * (u3 ** (mu - 2.0)) * v3 * v3 if x3 > 0 else 0.0)
                )
                assert neg_d2 >= bound_d2 - 1e-12, f"fail at {(x1,x2,x3)}, th={th}, s={s}: {neg_d2} < {bound_d2}"
                if neg_d2 < min_d2:
                    min_d2 = neg_d2

        # Check actual 1D pushforward density (numerical integral of t / (-psi'(t)) dtheta) at y = g_max - eps
        for eps in (1e-5, 1e-3, 0.05, 0.3):
            y = g_max - eps
            if y <= 0:
                continue
            total = 0.0
            dth = 2.0 * math.pi / n_theta
            for k in range(n_theta):
                th = (2.0 * math.pi * (k + 0.37)) / n_theta
                v1, v2 = math.cos(th), math.sin(th)
                v3 = -v1 - v2
                t_max = float("inf")
                enters = True
                for ui, vi in ((u1s, v1), (u2s, v2), (u3s, v3)):
                    if ui == 0.0 and vi <= 0.0:
                        enters = False
                        break
                    if vi < 0.0:
                        t_max = min(t_max, -ui / vi)
                if not enters or t_max <= 1e-15:
                    continue
                def psi(t):
                    return (
                        x1 * ((u1s + t * v1) ** mu)
                        + (x2 * ((u2s + t * v2) ** mu) if x2 > 0 else 0.0)
                        + (x3 * ((u3s + t * v3) ** mu) if x3 > 0 else 0.0)
                    )
                if psi(t_max * (1.0 - 1e-12)) >= y:
                    continue
                lo, hi = 0.0, t_max * (1.0 - 1e-12)
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if psi(mid) > y:
                        lo = mid
                    else:
                        hi = mid
                t = 0.5 * (lo + hi)
                u1, u2, u3 = u1s + t * v1, u2s + t * v2, u3s + t * v3
                neg_d1 = -mu * (
                    x1 * (u1 ** (mu - 1.0)) * v1
                    + (x2 * (u2 ** (mu - 1.0)) * v2 if x2 > 0 else 0.0)
                    + (x3 * (u3 ** (mu - 1.0)) * v3 if x3 > 0 else 0.0)
                )
                assert t / neg_d1 <= 16.0 / (mu * (1.0 - mu)) + 1e-9
                total += (t / neg_d1) * dth
            assert total <= bound_area
            if total > max_slice_density:
                max_slice_density = total

    # Verify exact 2-cut rate cancellation for s_i in {1, 2} and Dirichlet-3 density bound
    for m in range(2, 12):
        for s_i in (1, 2):
            for D in range(s_i + 1, m + s_i + 3):
                lam = Fraction(m - 1, 1) + Fraction(D - s_i, D)
                kap = Fraction(s_i, D)
                assert lam + kap == Fraction(m, 1)

    for N in (6, 15, 40):
        for n1 in (1, 2, N // 3):
            for n2 in (1, 2, N // 3):
                n3 = N - n1 - n2
                if n3 < 1:
                    continue
                # Peak of Dirichlet density at mode ((n1-1)/(N-3), (n2-1)/(N-3))
                u1 = (n1 - 1) / (N - 3) if n1 > 1 else 0.0
                u2 = (n2 - 1) / (N - 3) if n2 > 1 else 0.0
                u3 = 1.0 - u1 - u2
                log_f = (
                    math.lgamma(N) - math.lgamma(n1) - math.lgamma(n2) - math.lgamma(n3)
                    + ((n1 - 1) * math.log(u1) if n1 > 1 else 0.0)
                    + ((n2 - 1) * math.log(u2) if n2 > 1 else 0.0)
                    + ((n3 - 1) * math.log(u3) if n3 > 1 else 0.0)
                )
                assert math.exp(log_f) <= (N - 1) * (N - 2) + 1e-9

    print(
        f"Polar-Fubini & Dirichlet OK: min(-psi'')={min_d2:.6f} >= mu(1-mu)/16={bound_d2:.6f}; "
        f"max slice density={max_slice_density:.6f} <= 32*pi/(mu(1-mu))={bound_area:.3f}; "
        f"2-slot cancellation and Dir(n1,n2,n3) <= (N-1)(N-2) verified."
    )


if __name__ == "__main__":
    check_exponents()
    check_supermartingale_exact()
    check_azuma_gamma_bounds()
    check_combinatorial_packing()
    check_polar_fubini_and_dirichlet()
    print("ALL DENSITY-LP CHECKS PASSED.")
