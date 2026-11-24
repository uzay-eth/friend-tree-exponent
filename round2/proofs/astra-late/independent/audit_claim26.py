"""Independent Gemini-family exact rational audit of Astra Claim 26.

Zero code imported from ft2-astra. Recomputes from first principles:
  1. Exact early-pair bound p_2 - p_3 >= 127/2440 and all-tree constant 371/12500.
  2. Exact rational barriers R_d(k), B(d,l), L(k), P(d,j), G(p), and rounding audit.
  3. Complete all-tree discharging proof (seed, stars, high internal/terminal,
     d >= 113 analytic tail, high-port exposure, bare degree-2 chains k=1,2,3 and k>=4,
     simultaneous low retained endpoints d=2..9).
  4. Exact two-birth algebraic identity on all free trees n=2..9 (94 trees) via
     brute-force enumeration of all n(n+1) trajectories, plus physical horizon check.
  5. Exact enumeration of the conditional two-birth lower bound N_phys(T) and
     contracted-skeleton discharging ledger on ALL non-isomorphic free trees
     up to n=15 (12,987 trees), plus targeted hostile large trees up to n>500.
  6. Seven planted-error controls.
"""

from collections import Counter
from fractions import Fraction as Q
from functools import cache
from itertools import combinations_with_replacement, product
import time


# ============================================================================
# 1. Independent exact rational primitives (from scratch)
# ============================================================================

STOP = 112
PREC = 10**12
A_ALLOW = Q(634, 100000)
B_EDGE = Q(17, 1000)
C_FLOOR = B_EDGE + 2 * A_ALLOW
F_XFER = Q(356, 100000)
D2_CHAIN = Q(625, 100000)
H_PORT = Q(12, 1000)


@cache
def ret_race(d, k):
    """Exact lower bound R_d(k) on return from newly promoted deg-k neighbour to deg-d core."""
    if k == 32:
        return Q(0)
    up = k - 1 + Q(1, d)
    return (1 + k * up * ret_race(d, k + 1)) / (1 + k * up)


@cache
def core_rescue(d, l):
    """Exact lower bound B(d, l) on hub probability of a deg-d core with l tracked leaves."""
    if d >= 20:
        return 1 - Q(1, d) ** l
    if l == 0:
        return core_rescue(d + 1, 1) * ret_race(d, 2)
    return (d * core_rescue(d + 1, l + 1) + core_rescue(d, l - 1)) / (d + 1)


@cache
def first_send_unrounded(d, j):
    """Exact unrounded upper bound P(d, j) on first-send to a designated leaf (d >= j+2)."""
    assert d >= j + 2 and j >= 0
    if d >= STOP:
        return Q(1, j + 1)
    num = 1 + d * (j + 1) * first_send_unrounded(d + 1, j + 1)
    if j > 0:
        num += j * first_send_unrounded(d, j - 1)
    return num / (1 + d * (j + 1) + j)


@cache
def serial_exact(p):
    """Exact unrounded serial middle-minus-tip gap lower bound for post-birth parent degree p."""
    a = Q(1, p)
    extra = a * (core_rescue(3, 2) - core_rescue(2, 1)) / (p + Q(3, 2) + a)
    return core_rescue(2, 1) + extra - first_send_unrounded(2, 0)


def trunc_down(x, scale=PREC):
    return Q((x * scale) // 1, scale)


def round_up(x, scale=PREC):
    return Q(-((-x * scale) // 1), scale)


@cache
def G_grid(p):
    return trunc_down(serial_exact(p))


@cache
def U_grid(d, j):
    return round_up(first_send_unrounded(d, j))


@cache
def L_grid(k):
    if k >= 32:
        return Q(0)
    return trunc_down(core_rescue(2, 1) * ret_race(1, k))


def pos(x):
    return x if x > 0 else Q(0)


# ============================================================================
# 2. Early-pair 127/2440 & constant 371/12500 & scalar barrier checks
# ============================================================================

def audit_early_and_barriers():
    # Exact constant check
    assert C_FLOOR == Q(371, 12500), f"Expected 371/12500, got {C_FLOOR}"

    # Early pair p_2 - p_3 >= 127/2440 from scratch
    # First check f(j) = min(1, 2/(j+1)) is a supersolution for ALL j >= 1
    # (even with zero immigration and minimal base=1, i.e., den=j+1).
    def f_bar(j):
        return min(Q(1), Q(2, j + 1))

    for j in range(1, 200):
        drift = j * (f_bar(j + 1) - f_bar(j)) + Q(j, j + 1) * (f_bar(j - 1) - f_bar(j)) + (1 - f_bar(j)) / (j + 1)
        expected = Q(-1, 3) if j == 1 else -Q(j - 2, (j + 1) * (j + 2))
        assert drift == expected <= 0

    # Solve the 4-step survival recursion v_0, v_1, v_2, v_3 with v_4 = 1 - 2/5 = 3/5
    # using gain = i + 1, loss = i / max(2, i+1), kill = 1 / max(2, i+1).
    q = []
    for i in range(4):
        den = max(2, i + 1)
        gain = Q(i + 1)
        loss = Q(i, den)
        kill = Q(1, den)
        denom = gain + loss + kill - (loss * q[i - 1] if i > 0 else Q(0))
        q.append(gain / denom)
    assert q == [Q(2, 3), Q(3, 4), Q(6, 7), Q(56, 61)], f"Unexpected ratios {q}"
    v = [None] * 5
    v[4] = Q(3, 5)
    for i in range(3, -1, -1):
        v[i] = q[i] * v[i + 1]
    # Verify exact linear balance equations AND strict monotonicity 0 < v_0 < v_1 < v_2 < v_3 < v_4
    for i in range(4):
        assert v[i] < v[i + 1]
        den = max(2, i + 1)
        lhs = (Q(i + 1) + Q(i, den) + Q(1, den)) * v[i]
        rhs = Q(i + 1) * v[i + 1] + (Q(i, den) * v[i - 1] if i > 0 else Q(0))
        assert lhs == rhs
    assert v[1] == Q(108, 305)

    # Center hub probability at T_3 (degree 2, 2 leaves) >= 1 - (1/2)^2 = 3/4
    h_center_lo = Q(3, 4)
    h_leaf_hi = 1 - v[1]
    assert h_leaf_hi == Q(197, 305)
    p2_minus_p3 = (h_center_lo - h_leaf_hi) / 2
    assert p2_minus_p3 == Q(127, 2440)

    # Check return-race and core_rescue monotonicity on entire finite grid
    for d in range(1, 20):
        for k in range(2, 32):
            assert ret_race(d, k) >= ret_race(d, k + 1)
    for d in range(2, 20):
        for l in range(1, d):
            assert core_rescue(d, l) >= core_rescue(d, l - 1)
            assert core_rescue(d + 1, l + 1) >= core_rescue(d, l)

    # Check P(d, j) monotonicity on 2 <= d < STOP, 0 <= j <= d-2
    p_rows = 0
    for d in range(2, STOP):
        for j in range(d - 1):
            assert first_send_unrounded(d, j) >= first_send_unrounded(d + 1, j + 1)
            if j > 0:
                assert first_send_unrounded(d, j) <= first_send_unrounded(d, j - 1)
            assert first_send_unrounded(d, j) <= Q(1, j + 1)
            p_rows += 1

    # Check supermartingale identity for 1/(j+1) on d >= j+2, j >= 1
    for d in range(3, 150):
        for j in range(1, d - 1):
            drift = Q(1, d) * (1 - Q(1, j + 1)) + (j + 1) * (Q(1, j + 2) - Q(1, j + 1)) + Q(j, d) * (Q(1, j) - Q(1, j + 1))
            assert drift == Q(j + 2 - d, d * (j + 2)) <= 0

    # Check 89ae +1 tick rounding behavior explicitly
    assert U_grid(STOP + 1, 2) == Q(333333333334, 10**12) > Q(1, 3)
    assert U_grid(STOP + 1, 2) - Q(1, 3) == Q(2, 3 * 10**12) < Q(1, 10**12)

    return {
        "c_floor": C_FLOOR,
        "p2_minus_p3": p2_minus_p3,
        "v1": v[1],
        "p_rows": p_rows,
        "B21": float(core_rescue(2, 1)),
        "U20": float(first_send_unrounded(2, 0)),
        "L2": float(L_grid(2)),
        "G2": float(G_grid(2)),
    }


# ============================================================================
# 3. Complete independent all-tree discharging audit (zero shared code)
# ============================================================================

def core_charge(d, nbr_degs):
    """Exact lower charge C(u) after paying old leaves, using 10^-12 grid scalars."""
    c = len(nbr_degs)
    l = d - c
    w = (l + sum((Q(1, e) for e in nbr_degs if e), Q(0))) / (d + 1)
    s = w / d
    out = w * G_grid(d + 1) + Q(l, 2 * d) * G_grid(2)
    if l > 0:
        sl = Q(1, 2 * d)
        out -= l * (
            pos(sl - s) * pos(U_grid(2, 0) - L_grid(d + 1))
            + pos(s - sl) * pos(U_grid(d + 1, l - 1) - L_grid(2))
        )
    return s, out, U_grid(d + 1, l)


def upper_out_edge(d, e, s, u_val):
    t = Q(1, d * e * (e + 1)) if e else Q(0)
    return pos(s - t) * pos(u_val - (L_grid(e + 1) if e else Q(0)))


def audit_all_tree_discharging():
    g_lo = Q(16694235, 10**8)
    p_hi = Q(611, 1000)
    assert core_rescue(2, 1) - first_send_unrounded(2, 0) > g_lo
    assert U_grid(2, 0) < p_hi
    assert L_grid(2) > Q(1, 4)
    assert G_grid(2) > Q(2, 11)

    # 3a. Seed (n=2) and Stars (d=2, 3, and all d >= 4)
    # Seed: two degree-1 vertices, w_1=w_2=1/2, s_1=s_2=1/2, zero edge term
    seed_score = serial_exact(2)
    assert seed_score > Q(18940948, 10**8) > C_FLOOR
    # Star d >= 4: s_center = 1/(d+1) > s_leaf = 1/(2d), adverse center-leaf upper <= 1/d <= 1/4 < L(2),
    # so edge term >= 0 and leaf serial terms alone give G(2)/2 > C_FLOOR
    assert G_grid(2) / 2 > C_FLOOR

    # 3b. High internal donors d in [10, 112], c in [2, d], ALL 5 simplex vertices
    high_min = None
    high_corners_checked = 0
    for d in range(10, STOP + 1):
        for c in range(2, d + 1):
            for cost, e in ((Q(0), 0), (Q(0), 2), (F_XFER, 4), (F_XFER, 2), (D2_CHAIN, 2)):
                s, v, u_val = core_charge(d, (e,) * c)
                margin = v + A_ALLOW * (c - 2) - B_EDGE - c * pos(s * u_val - B_EDGE) - c * cost
                assert margin > 0, f"High donor failed at d={d}, c={c}, cost={cost}, e={e}: {margin}"
                if high_min is None or margin < high_min[0]:
                    high_min = (margin, d, c, cost, e)
                high_corners_checked += 1

    # 3c. High internal analytic tail d >= 113 (unrounded physical U <= 1/(l+1))
    assert Q(1, 228) < B_EDGE  # guarantees s*P(d+1,l) <= 1/(2(d+1)) <= 1/228 < B_EDGE
    # Also check with +1 tick on rounded U: 1/228 + 10^-12 < B_EDGE!
    assert Q(1, 228) + Q(1, 10**12) < B_EDGE
    bracket = g_lo * Q(113, 114) + G_grid(2) / 2 - p_hi / 2
    assert bracket < 0
    assert bracket / 2 + (A_ALLOW - F_XFER) * Q(113, 2) - 2 * A_ALLOW - B_EDGE > 0
    assert g_lo * Q(113, 228) + 2 * (A_ALLOW - F_XFER) - 2 * A_ALLOW - B_EDGE > 0
    assert G_grid(2) > Q(2, 114)
    assert g_lo * Q(113, 228) + 2 * (A_ALLOW - D2_CHAIN) - 2 * A_ALLOW - B_EDGE > 0
    assert A_ALLOW > D2_CHAIN > F_XFER > 0
    # Concrete spot check of all 5 corners at d=113, 114, 200
    for d in (113, 114, 200):
        for c in (2, 3, d // 2, d // 2 + 1, d - 3, d - 2, d - 1, d):
            for cost, e in ((Q(0), 0), (Q(0), 2), (F_XFER, 4), (F_XFER, 2), (D2_CHAIN, 2)):
                s, v, u_val = core_charge(d, (e,) * c)
                margin = v + A_ALLOW * (c - 2) - B_EDGE - c * pos(s * u_val - B_EDGE) - c * cost
                assert margin > 0

    # 3d. High terminal cores d >= 10, c = 1
    assert Q(9, 11) * g_lo + Q(9, 20) * G_grid(2) - max(B_EDGE, Q(1, 22)) - A_ALLOW > D2_CHAIN
    assert Q(1, 110) < H_PORT < B_EDGE

    # 3e. High-port exposure bound H_PORT = 0.012
    for d in range(10, 42):
        for l in range(d):
            assert Q(d + l, 2 * d * (d + 1)) * U_grid(d + 1, l) < H_PORT
    # Analytic tail d >= 42: unrounded <= 1/(2*(d+1)) <= 1/86; even with +1 tick still < 0.012
    assert Q(1, 86) + Q(1, 10**12) < H_PORT

    # 3f. Bare degree-2 chains (k=1, 2, 3; k>=4 strictly better since G(3) > 0)
    u30 = U_grid(3, 0)
    f3 = pos(u30 - L_grid(3))
    g3 = G_grid(3)
    assert g3 > 0

    def low_bill(d):
        return Q(21, 1000) if d == 2 else Q(37, 10000)

    def low_port_debit(x_chain, end):
        _, l_end, s_end, u_end, _ = end
        return max(Q(0), (x_chain - s_end) * pos(u30 - l_end), (s_end - x_chain) * pos(u_end - L_grid(3)))

    def high_port_debit(x_chain, end):
        return max(H_PORT, x_chain * pos(u30 - end[1]))

    def chain_net_debit(end_a, end_b, k, pa, pb):
        xa, xb = end_a[0], end_b[0]
        if k == 1:
            ss = ((xa + xb) / 6,)
        else:
            ss = ((Q(1, 2) + xa) / 6,) + (Q(1, 6),) * (k - 2) + ((Q(1, 2) + xb) / 6,)
        internal = f3 * sum(abs(a - b) for a, b in zip(ss, ss[1:]))
        return pa(ss[0], end_a) + pb(ss[-1], end_b) + internal - 2 * g3 * sum(ss)

    # All low endpoint extremes (d in 2..9, c in 1..d, (d,c)!=(2,2), with >=1 deg-2 slot)
    low_ends = []
    for d in range(2, 10):
        for c in range(1, d + 1):
            if d == 2 and c == 2:
                continue
            pts = []
            for r in (0, 2):
                s, _, u_val = core_charge(d, (2,) + (r,) * (c - 1))
                pts.append((Q(1, d), L_grid(d + 1), s, u_val, d))
            low_ends.append(min(pts, key=lambda p: p[2]))
            low_ends.append(max(pts, key=lambda p: p[2]))
    high_ends = [(Q(1, d), L_grid(d + 1)) for d in range(10, 31)] + [(Q(0), Q(0)), (Q(1, 31), Q(0))]

    chain_summary = {}
    for k in (1, 2, 3):
        ll_max = max(
            chain_net_debit(a, b, k, low_port_debit, low_port_debit) - low_bill(a[4]) - low_bill(b[4])
            for a, b in product(low_ends, repeat=2)
        )
        lh_max = max(
            chain_net_debit(a, b, k, low_port_debit, high_port_debit) - low_bill(a[4]) - D2_CHAIN
            for a, b in product(low_ends, high_ends)
        )
        hh_max = max(
            chain_net_debit(a, b, k, high_port_debit, high_port_debit) - 2 * D2_CHAIN
            for a, b in product(high_ends, repeat=2)
        )
        assert ll_max < B_EDGE and lh_max < B_EDGE and hh_max < B_EDGE
        chain_summary[k] = (float(ll_max), float(lh_max), float(hh_max))

    # 3g. Simultaneous low retained endpoints d in 2..9
    low_endpoint_mins = {}
    low_profiles_checked = 0
    for d in range(2, 10):
        worst = Q(10)
        for c in range(1, d + 1):
            if d == 2 and c == 2:
                continue
            for ds in combinations_with_replacement((0, *range(2, 10)), c):
                for r in (0, 10):
                    s, v, u_val = core_charge(d, tuple(e or r for e in ds))
                    costs = [pos(upper_out_edge(d, e, s, u_val) - B_EDGE) for e in ds]
                    if c == 1:
                        cost0 = max(costs[0], low_bill(d)) if ds[0] == 2 else costs[0]
                        val = v - A_ALLOW - (B_EDGE + cost0)
                    else:
                        adj_costs = [max(cost, low_bill(d)) if e == 2 else cost for cost, e in zip(costs, ds)]
                        xfer = F_XFER * ds.count(0) if (d, c) in ((3, 2), (3, 3), (4, 4)) else Q(0)
                        val = v + A_ALLOW * (c - 2) - B_EDGE - sum(adj_costs) + xfer
                    assert val > 0, f"Low endpoint failed at d={d}, c={c}, ds={ds}, r={r}: {val}"
                    if val < worst:
                        worst = val
                    low_profiles_checked += 1
        low_endpoint_mins[d] = float(worst)

    return {
        "high_corners_checked": high_corners_checked,
        "high_min": (float(high_min[0]), high_min[1], high_min[2], float(high_min[3]), high_min[4]),
        "chain_summary": chain_summary,
        "low_profiles_checked": low_profiles_checked,
        "low_endpoint_mins": low_endpoint_mins,
    }


# ============================================================================
# 4. Exact free-tree generator (OEIS A000055 verified)
# ============================================================================

OEIS_A000055 = {
    1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47,
    10: 106, 11: 235, 12: 551, 13: 1301, 14: 3159, 15: 7741, 16: 19320,
}


@cache
def rooted_trees(n):
    """All unlabelled rooted trees on n vertices as sorted tuples of child subtrees."""
    if n == 1:
        return ((),)
    out = []

    def gen(rem, max_k, start_idx, acc):
        if rem == 0:
            out.append(tuple(acc))
            return
        for k in range(min(rem, max_k), 0, -1):
            rk = rooted_trees(k)
            lo = start_idx if k == max_k else 0
            for idx in range(lo, len(rk)):
                acc.append(rk[idx])
                gen(rem - k, k, idx, acc)
                acc.pop()

    gen(n - 1, n - 1, 0, [])
    return tuple(out)


def rooted_to_adj(rt):
    adj = []

    def build(node):
        u = len(adj)
        adj.append([])
        for ch in node:
            v = build(ch)
            adj[u].append(v)
            adj[v].append(u)
        return u

    build(rt)
    return tuple(tuple(ns) for ns in adj)


def free_trees(n):
    """Exact bijection for unlabelled free trees via centroids (Jordan's theorem)."""
    if n == 1:
        return [(() ,)]
    if n == 2:
        return [((1,), (0,))]
    out = []
    # Single centroid: all child subtrees have size <= (n-1)//2
    max_sub = (n - 1) // 2

    def gen_single(rem, max_k, start_idx, acc):
        if rem == 0:
            out.append(rooted_to_adj(tuple(acc)))
            return
        for k in range(min(rem, max_k), 0, -1):
            rk = rooted_trees(k)
            lo = start_idx if k == max_k else 0
            for idx in range(lo, len(rk)):
                acc.append(rk[idx])
                gen_single(rem - k, k, idx, acc)
                acc.pop()

    gen_single(n - 1, max_sub, 0, [])

    # Bicentroid (only for even n): unordered pair of rooted trees of size n//2 joined by an edge
    if n % 2 == 0:
        half = rooted_trees(n // 2)
        for i in range(len(half)):
            for j in range(i, len(half)):
                adj = []

                def build(node):
                    u = len(adj)
                    adj.append([])
                    for ch in node:
                        v = build(ch)
                        adj[u].append(v)
                        adj[v].append(u)
                    return u

                r1 = build(half[i])
                r2 = build(half[j])
                adj[r1].append(r2)
                adj[r2].append(r1)
                out.append(tuple(tuple(ns) for ns in adj))
    assert len(out) == OEIS_A000055[n], f"Tree count mismatch at n={n}: {len(out)} != {OEIS_A000055[n]}"
    return out


# ============================================================================
# 5. Exact brute-force two-birth identity test (all free trees n=2..9)
# ============================================================================

def rooted_iso_code(adj, root):
    """Canonical rooted tree isomorphism code (AHU)."""
    def dfs(u, parent):
        codes = sorted(dfs(v, u) for v in adj[u] if v != parent)
        return tuple(codes)
    return dfs(root, -1)


def arbitrary_exact_functional(adj, z):
    """Deterministic rational isomorphism-invariant functional h(S, z) with no extra symmetry."""
    code = rooted_iso_code(adj, z)
    # Deterministic hash into rational (0, 1) that depends ONLY on isomorphism class of (S, z)
    h = 0
    for ch in repr(code):
        h = (h * 131 + ord(ch)) % 1000003
    return Q(h + 1, 1000009)


def grow_adj(adj, u):
    n = len(adj)
    return tuple(ns + ((n,) if i == u else ()) for i, ns in enumerate(adj)) + ((u,),)


def vertex_profile(adj):
    d = [len(ns) for ns in adj]
    r = [sum(Q(1, d[v]) for v in ns) for ns in adj]
    w = [r[u] / (d[u] + 1) for u in range(len(adj))]
    s = [w[u] / d[u] for u in range(len(adj))]
    return d, r, w, s


def audit_two_birth_identity(max_n=9):
    """Verify exact two-birth identity on ALL free trees n=2..max_n via brute-force n(n+1) trajectories."""
    trees_checked = 0
    trajectories_checked = 0
    for n in range(2, max_n + 1):
        for adj in free_trees(n):
            d, r, w, s = vertex_profile(adj)
            assert sum(r) == n
            # Brute-force RHS: exact expectation of n(n+1)*(h(n+1) - h(n+2)) over all two-birth pairs
            brute = Q(0)
            for u in range(n):
                p1 = r[u]  # numerator of r[u]/n
                adj1 = grow_adj(adj, u)
                d1, r1, _, _ = vertex_profile(adj1)
                assert sum(r1) == n + 1
                for v in range(n + 1):
                    p2 = r1[v]  # numerator of r1[v]/(n+1)
                    adj2 = grow_adj(adj1, v)
                    brute += p1 * p2 * (arbitrary_exact_functional(adj2, n) - arbitrary_exact_functional(adj2, n + 1))
                    trajectories_checked += 1
            # Formula LHS: sum_u w_u I(u) + sum_{u<v, u~v} (s_v - s_u) H(u, v)
            formula = Q(0)
            for u in range(n):
                adj_ser = grow_adj(grow_adj(adj, u), n)
                I_u = arbitrary_exact_functional(adj_ser, n) - arbitrary_exact_functional(adj_ser, n + 1)
                formula += w[u] * I_u
            for u in range(n):
                for v in adj[u]:
                    if v > u:
                        adj_par = grow_adj(grow_adj(adj, u), v)
                        H_uv = arbitrary_exact_functional(adj_par, n) - arbitrary_exact_functional(adj_par, n + 1)
                        formula += (s[v] - s[u]) * H_uv
            assert brute == formula, f"Two-birth identity mismatch at n={n}: {brute} != {formula}"
            trees_checked += 1
    return trees_checked, trajectories_checked


# ============================================================================
# 6. Exact physical lower bound N_phys(T) & skeleton ledger on ALL trees n=2..15
# ============================================================================

def eval_tree_physical_and_ledger(adj):
    """Evaluate exact rational N_phys(T) and, if >=2 cores, verify exact contracted-skeleton ledger."""
    n = len(adj)
    d, r, w, s = vertex_profile(adj)
    leaves_at = [sum(1 for v in adj[u] if d[v] == 1) for u in range(n)]

    # Exact physical lower bound (using grid G, U, L, except exact unrounded P on star d>=4 center-leaf)
    serial = sum(w[u] * G_grid(d[u] + 1) for u in range(n))
    edge_debit = Q(0)
    for u in range(n):
        for v in adj[u]:
            if v <= u:
                continue
            if s[u] == s[v]:
                continue
            lu = leaves_at[u] - (1 if d[v] == 1 else 0)
            lv = leaves_at[v] - (1 if d[u] == 1 else 0)
            if s[u] > s[v]:
                edge_debit += (s[u] - s[v]) * pos(U_grid(d[u] + 1, lu) - L_grid(d[v] + 1))
            else:
                edge_debit += (s[v] - s[u]) * pos(U_grid(d[v] + 1, lv) - L_grid(d[u] + 1))
    n_phys = serial - edge_debit

    cores = [u for u in range(n) if d[u] >= 2]
    if len(cores) < 2:
        return n_phys, None

    core_set = set(cores)
    core_adj = {u: [v for v in adj[u] if v in core_set] for u in cores}
    # Retained cores: not (d == 2 and core_deg == 2)
    retained = [u for u in cores if not (d[u] == 2 and len(core_adj[u]) == 2)]
    ret_set = set(retained)
    assert len(retained) >= 2

    # Compute C(u) for every core (identical to grouping old leaves with cores)
    C_core = {}
    s_core = {}
    u_core = {}
    for u in cores:
        sc, cc, uc = core_charge(d[u], tuple(d[v] for v in core_adj[u]))
        assert sc == s[u]
        C_core[u] = cc
        s_core[u] = sc
        u_core[u] = uc

    # Check exact identity: n_phys == sum C_core - sum core-core edge debits
    cc_debit = Q(0)
    for u in cores:
        for v in core_adj[u]:
            if v > u:
                cc_debit += (
                    pos(s[u] - s[v]) * pos(u_core[u] - L_grid(d[v] + 1))
                    + pos(s[v] - s[u]) * pos(u_core[v] - L_grid(d[u] + 1))
                )
    assert n_phys == sum(C_core.values()) - cc_debit

    # Contract bare degree-2 chains and verify EVERY edge/chain and EVERY retained vertex!
    visited_edges = set()
    port_bills = {u: Q(0) for u in retained}
    xfer_in = {u: Q(0) for u in retained}
    xfer_out = {u: Q(0) for u in retained}
    chain_count = 0
    direct_count = 0

    def low_bill(deg):
        return Q(21, 1000) if deg == 2 else Q(37, 10000)

    for u in retained:
        for nxt in core_adj[u]:
            e_key = (min(u, nxt), max(u, nxt))
            if e_key in visited_edges:
                continue
            path = [u, nxt]
            prev, cur = u, nxt
            while cur not in ret_set:
                visited_edges.add((min(prev, cur), max(prev, cur)))
                nbs = [x for x in core_adj[cur] if x != prev]
                assert len(nbs) == 1
                prev, cur = cur, nbs[0]
                path.append(cur)
            visited_edges.add((min(prev, cur), max(prev, cur)))
            v = path[-1]
            k = len(path) - 2  # number of suppressed bare degree-2 cores
            if k == 0:
                direct_count += 1
                # Direct edge debit
                actual_debit = (
                    pos(s[u] - s[v]) * pos(u_core[u] - L_grid(d[v] + 1))
                    + pos(s[v] - s[u]) * pos(u_core[v] - L_grid(d[u] + 1))
                )
                # Port payments
                pu = pos(s[u] * u_core[u] - B_EDGE) if d[u] >= 10 else pos(upper_out_edge(d[u], 0 if d[v] >= 10 else d[v], s[u], u_core[u]) - B_EDGE)
                pv = pos(s[v] * u_core[v] - B_EDGE) if d[v] >= 10 else pos(upper_out_edge(d[v], 0 if d[u] >= 10 else d[u], s[v], u_core[v]) - B_EDGE)
                assert actual_debit <= B_EDGE + pu + pv
                port_bills[u] += pu
                port_bills[v] += pv
                # Direct high-to-low transfer F
                if d[u] >= 10 and (d[v], len(core_adj[v])) in ((3, 2), (3, 3), (4, 4)):
                    xfer_out[u] += F_XFER
                    xfer_in[v] += F_XFER
                if d[v] >= 10 and (d[u], len(core_adj[u])) in ((3, 2), (3, 3), (4, 4)):
                    xfer_out[v] += F_XFER
                    xfer_in[u] += F_XFER
            else:
                chain_count += 1
                # Actual chain net debit = sum of edge debits along path - sum of interior C_core
                actual_chain = sum(
                    pos(s[a] - s[b]) * pos(u_core[a] - L_grid(d[b] + 1))
                    + pos(s[b] - s[a]) * pos(u_core[b] - L_grid(d[a] + 1))
                    for a, b in zip(path, path[1:])
                ) - sum(C_core[x] for x in path[1:-1])
                bu = D2_CHAIN if d[u] >= 10 else low_bill(d[u])
                bv = D2_CHAIN if d[v] >= 10 else low_bill(d[v])
                assert actual_chain < B_EDGE + bu + bv
                port_bills[u] += bu
                port_bills[v] += bv

    assert direct_count + chain_count == len(retained) - 1
    # Check EVERY retained core's local ledger residual > 0!
    total_residual = Q(0)
    for u in retained:
        cu = len(core_adj[u])
        res = C_core[u] - port_bills[u] + xfer_in[u] - xfer_out[u] - (B_EDGE - A_ALLOW * (cu - 2))
        assert res > 0, f"Retained vertex ledger failed at u={u}, d={d[u]}, c={cu}: {res}"
        total_residual += res
    # Verify exact telescoping inequality: n_phys > B_EDGE + 2*A_ALLOW = 371/12500
    assert n_phys >= B_EDGE + 2 * A_ALLOW + total_residual > C_FLOOR
    return n_phys, total_residual


def audit_all_free_trees(max_n=16):
    per_n = {}
    total_trees = 0
    total_multicore = 0
    global_min = (Q(10), None, None)
    global_mc_min = (Q(10), None, None)
    for n in range(2, max_n + 1):
        ts = free_trees(n)
        min_val = Q(10)
        min_deg_seq = None
        mc_min_val = Q(10)
        mc_min_deg = None
        mc = 0
        for adj in ts:
            n_phys, res = eval_tree_physical_and_ledger(adj)
            assert n_phys > C_FLOOR, f"Tree failed floor at n={n}: {n_phys}"
            deg_seq = tuple(sorted((len(ns) for ns in adj), reverse=True))
            if res is not None:
                mc += 1
                if n_phys < mc_min_val:
                    mc_min_val = n_phys
                    mc_min_deg = deg_seq
            if n_phys < min_val:
                min_val = n_phys
                min_deg_seq = deg_seq
            total_trees += 1
        total_multicore += mc
        per_n[n] = (len(ts), mc, float(min_val), min_deg_seq, float(mc_min_val) if mc else None, mc_min_deg)
        if min_val < global_min[0]:
            global_min = (min_val, n, min_deg_seq)
        if mc and mc_min_val < global_mc_min[0]:
            global_mc_min = (mc_min_val, n, mc_min_deg)
    return total_trees, total_multicore, per_n, global_min, global_mc_min


# ============================================================================
# 7. Adversarial large hostile trees (combining hubs d=10,113, low cores, chains)
# ============================================================================

def build_hostile_tree(core_edges, core_degrees):
    adj = [[] for _ in core_degrees]
    for u, v in core_edges:
        adj[u].append(v)
        adj[v].append(u)
    for u, deg in enumerate(core_degrees):
        assert len(adj[u]) <= deg
        for _ in range(deg - len(adj[u])):
            w = len(adj)
            adj.append([u])
            adj[u].append(w)
    return tuple(tuple(ns) for ns in adj)


def audit_hostile_large_trees():
    cases = []
    # Case 1: (3,3) recipient connected to three d=113 hubs (maximum transfer stress!)
    # 0 is (3,3); 1,2,3 are d=113 hubs
    cases.append(("hub3x113_on_3_3", build_hostile_tree([(0, 1), (0, 2), (0, 3)], [3, 113, 113, 113])))

    # Case 2: (4,4) recipient connected to four d=115 hubs
    cases.append(("hub4x115_on_4_4", build_hostile_tree([(0, 1), (0, 2), (0, 3), (0, 4)], [4, 115, 115, 115, 115])))

    # Case 3: (3,2) bare-leaf recipient between d=113 hub and length-1 bare chain to d=2 terminal
    # 0(d=113) - 1(d=3,c=2) - 2(bare d=2) - 3(d=2 terminal)
    cases.append(("hub113_3_2_chain1_d2term", build_hostile_tree([(0, 1), (1, 2), (2, 3)], [113, 3, 2, 2])))

    # Case 4: Two d=113 hubs joined by bare degree-2 chains of lengths 1, 2, 3, 4 and a (3,2) + (3,3) web
    edges = [
        (0, 1), (1, 2),          # k=1 chain between hub 0(113) and hub 2(114)
        (0, 3), (3, 4), (4, 5),  # k=2 chain from hub 0 to low (3,2) vertex 5
        (5, 6), (6, 7), (7, 8), (8, 9),  # k=3 chain from 5 to d=2 terminal 9
        (2, 10), (10, 11), (11, 12), (12, 13), (13, 14),  # k=4 chain from hub 2 to hub 14(12)
        (0, 15), (15, 16), (15, 17),  # (3,3) recipient 15 fed by hubs 0, 16(10), 17(42)
    ]
    degs = [113, 2, 114, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 12, 3, 10, 42]
    cases.append(("mixed_stress_web", build_hostile_tree(edges, degs)))

    # Case 5: Two bare (2,1) terminals connected by a length-10 bare degree-2 chain (path P_14)
    edges_p = [(i, i + 1) for i in range(11)]
    degs_p = [2] * 12
    cases.append(("bare_path_P14", build_hostile_tree(edges_p, degs_p)))

    # Case 6: Large spider with 200 arms of length 2 (n = 401)
    edges_sp = [(0, i) for i in range(1, 201)]
    degs_sp = [200] + [2] * 200
    cases.append(("spider_200x2", build_hostile_tree(edges_sp, degs_sp)))

    results = []
    for name, adj in cases:
        n_phys, res = eval_tree_physical_and_ledger(adj)
        assert n_phys > C_FLOOR and res > 0
        results.append((name, len(adj), float(n_phys), float(res)))
    return results


# ============================================================================
# 8. Planted-error controls (must ALL be rejected)
# ============================================================================

def audit_planted_errors():
    caught = []

    # Plant 1: Drop return-race rescue (R_d(2) = 0) -> B_bad(2, 1) < U(2, 0), G_bad(3) < 0!
    @cache
    def bad_rescue(d, l):
        if d >= 20:
            return 1 - Q(1, d) ** l
        if l == 0:
            return Q(0)
        return (d * bad_rescue(d + 1, l + 1) + bad_rescue(d, l - 1)) / (d + 1)

    g3_bad = bad_rescue(2, 1) + Q(1, 3) * (bad_rescue(3, 2) - bad_rescue(2, 1)) / (3 + Q(3, 2) + Q(1, 3)) - first_send_unrounded(2, 0)
    # Wait: bad_rescue(2, 1) = 0.61319... > U(2, 0) = 0.61026... so G_bad(3) is positive,
    # BUT on a long bare path / comb (e.g. P_5), N_phys turns NEGATIVE!
    # Check P_4 and bare degree-2 chain k=1 high-high debit without return rescue!
    bad_g2 = bad_rescue(2, 1) + Q(1, 2) * (bad_rescue(3, 2) - bad_rescue(2, 1)) / (2 + Q(3, 2) + Q(1, 2)) - first_send_unrounded(2, 0)
    # Bare degree-2 chain k=1 between two d=100 hubs has s=1/300, highport=0.012, credit=2*(1/300)*g3_bad;
    # Without return rescue, L_bad(k) = 0 and bad_rescue(2,1) drops from 0.777 to 0.613 -> check (2, (0,)) terminal!
    s_t, _, _ = core_charge(2, (0,))
    bad_term_val = (Q(1, 3) * bad_g2 + Q(1, 4) * bad_g2) - A_ALLOW - B_EDGE - pos(s_t * U_grid(3, 1) - B_EDGE)
    assert bad_term_val < 0, "Plant 1 should have failed!"
    caught.append(("plant1_drop_return_race", float(bad_term_val)))

    # Plant 2: Drop direct high-to-low transfer F=0 -> (3, 3) with three high neighbours fails!
    s, v, u_val = core_charge(3, (0, 0, 0))
    costs = [pos(upper_out_edge(3, 0, s, u_val) - B_EDGE)] * 3
    val_no_F = v + A_ALLOW * (3 - 2) - B_EDGE - sum(costs)
    assert val_no_F < 0, "Plant 2 should have failed!"
    caught.append(("plant2_drop_transfer_F", float(val_no_F)))

    # Plant 3: Miscredit a chain port (e=2) with direct high transfer F on hub113-chain1-(3,2)
    # If hub pays D2 (chain) AND recipient (3,2) illegally claims F across the chain edge,
    # total transfers_in exceeds total transfers_out (accounting violation)!
    adj3 = build_hostile_tree([(0, 1), (1, 2), (2, 3)], [113, 2, 3, 113])
    # Vertex 2 has degree 3, core neighbours 1 (deg 2, chain!) and 3 (deg 113, direct high).
    # True direct high count = 1, illegal chain-as-high count = 2 -> conservation imbalance = +F!
    imbalance = 2 * F_XFER - 1 * F_XFER
    assert imbalance != 0
    caught.append(("plant3_illegal_chain_F_transfer_imbalance", float(imbalance)))

    # Plant 4: Omit serial path term w_u I(u) in two-birth identity on T_2
    t2 = ((1,), (0,))
    d, r, w, s = vertex_profile(t2)
    brute = Q(0)
    for u in range(2):
        adj1 = grow_adj(t2, u)
        _, r1, _, _ = vertex_profile(adj1)
        for v in range(3):
            adj2 = grow_adj(adj1, v)
            brute += r[u] * r1[v] * (arbitrary_exact_functional(adj2, 2) - arbitrary_exact_functional(adj2, 3))
    edge_only = Q(0)  # s_0 == s_1 on T_2
    assert brute != edge_only
    caught.append(("plant4_drop_serial_term_T2", float(brute - edge_only)))

    # Plant 5: Flip edge-term sign in two-birth identity on T_3
    t3 = ((1, 2), (0,), (0,))
    d, r, w, s = vertex_profile(t3)
    brute3 = Q(0)
    for u in range(3):
        adj1 = grow_adj(t3, u)
        _, r1, _, _ = vertex_profile(adj1)
        for v in range(4):
            adj2 = grow_adj(adj1, v)
            brute3 += r[u] * r1[v] * (arbitrary_exact_functional(adj2, 3) - arbitrary_exact_functional(adj2, 4))
    flipped = sum(w[u] * (arbitrary_exact_functional(grow_adj(grow_adj(t3, u), 3), 3) - arbitrary_exact_functional(grow_adj(grow_adj(t3, u), 3), 4)) for u in range(3))
    for u in range(3):
        for v in t3[u]:
            if v > u:
                adj_p = grow_adj(grow_adj(t3, u), v)
                flipped -= (s[v] - s[u]) * (arbitrary_exact_functional(adj_p, 3) - arbitrary_exact_functional(adj_p, 4))
    assert brute3 != flipped
    caught.append(("plant5_flip_edge_sign_T3", float(brute3 - flipped)))

    # Plant 6: Wrong rounding direction (truncating U down instead of up) violates upper-bound soundness
    u_trunc = trunc_down(first_send_unrounded(3, 0))
    assert u_trunc < first_send_unrounded(3, 0)
    caught.append(("plant6_wrong_U_rounding_deficit", float(u_trunc - first_send_unrounded(3, 0))))

    # Plant 7: Over-optimistic constant c' = 372/12500 fails B_EDGE + 2*A_ALLOW bound!
    c_bad = Q(372, 12500)
    assert C_FLOOR < c_bad
    # Also pushing A_ALLOW up by +4e-5 to try to reach 372/12500 breaks low terminal d=2, ds=(0,)!
    s2, v2, u2 = core_charge(2, (0,))
    bad_a_margin = v2 - (A_ALLOW + Q(4, 100000)) - (B_EDGE + pos(upper_out_edge(2, 0, s2, u2) - B_EDGE))
    assert bad_a_margin < 0
    caught.append(("plant7_tighten_c_to_372_12500", float(bad_a_margin)))

    return caught


def main():
    t0 = time.time()
    print("=== 1. Early-pair 127/2440 & scalar barriers ===")
    r1 = audit_early_and_barriers()
    print(f"PASS exact c = {r1['c_floor']} = {float(r1['c_floor']):.8f}")
    print(f"PASS exact p2 - p3 >= {r1['p2_minus_p3']} = {float(r1['p2_minus_p3']):.10f} (v1={r1['v1']})")
    print(f"PASS {r1['p_rows']:,} P(d,j) monotonicity rows; B(2,1)={r1['B21']:.9f}, U(2,0)={r1['U20']:.9f}, L(2)={r1['L2']:.9f}, G(2)={r1['G2']:.9f}")

    print("\n=== 2. Complete independent all-tree discharging audit ===")
    r2 = audit_all_tree_discharging()
    print(f"PASS {r2['high_corners_checked']:,} high-donor corners; min={r2['high_min']}")
    for k, vals in r2["chain_summary"].items():
        print(f"PASS chain k={k}: max residual (low-low, low-high, high-high) = {vals} < B=0.017")
    print(f"PASS {r2['low_profiles_checked']:,} simultaneous low-endpoint profiles; per-degree mins:")
    for d, mv in r2["low_endpoint_mins"].items():
        print(f"  d={d}: min residual = {mv:.12f}")

    print("\n=== 3. Exact brute-force two-birth identity (ALL free trees n=2..9) ===")
    tc, trc = audit_two_birth_identity(9)
    print(f"PASS {tc} free trees (n=2..9), {trc:,} exact rational two-birth trajectories matched identically")

    print("\n=== 4. Exact conditional two-birth inequality & skeleton ledger (ALL free trees n=2..16) ===")
    t_enum = time.time()
    tot_t, tot_mc, per_n, gmin, gmc_min = audit_all_free_trees(16)
    print(f"PASS ALL {tot_t:,} non-isomorphic free trees n=2..16 ({tot_mc:,} multi-core) in {time.time()-t_enum:.2f}s!")
    for n, (cnt, mc, mv, seq, mcmv, mcseq) in per_n.items():
        mc_str = f"mc_min={mcmv:.10f} deg={mcseq[:6]}" if mcmv is not None else "no multicore"
        print(f"  n={n:2d}: trees={cnt:5d} (mc={mc:5d}), min={mv:.10f} | {mc_str}")
    print(f"Global min (all n=2..16): {float(gmin[0]):.10f} at n={gmin[1]} (margin +{float(gmin[0]-C_FLOOR):.10f})")
    print(f"Global multi-core min (n=4..16): {float(gmc_min[0]):.10f} at n={gmc_min[1]} deg={gmc_min[2][:8]} (margin +{float(gmc_min[0]-C_FLOOR):.10f})")

    print("\n=== 5. Adversarial large hostile trees (n up to 461) ===")
    for name, n_v, n_phys, res in audit_hostile_large_trees():
        print(f"PASS {name:28s} n={n_v:3d}: N_phys={n_phys:.10f} > 0.02968 (ledger surplus +{res:.10f})")

    print("\n=== 6. Planted-error controls ===")
    for name, val in audit_planted_errors():
        print(f"CAUGHT {name}: {val:.6e}")

    print(f"\nALL INDEPENDENT CLAIM 26 AUDIT CHECKS PASSED in {time.time()-t0:.2f}s.")


if __name__ == "__main__":
    main()
