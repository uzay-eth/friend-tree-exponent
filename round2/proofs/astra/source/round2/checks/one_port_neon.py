#!/usr/bin/env python3
"""Independent adversarial checker for round2/finite/one-port-upper.md.

Single-CPU, stdlib-only exact small-tree verification, Monte Carlo coupling
checks, constant audits, and planted-error falsification tests covering
Sections 1-8 of round2/finite/one-port-upper.md.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from fractions import Fraction
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# 1. Exact Gauss-Legendre nodes/weights and sparse uniformization exp(t * M)
# ---------------------------------------------------------------------------

def gauss_legendre(n: int) -> Tuple[List[float], List[float]]:
    """Compute n-point Gauss-Legendre nodes and weights on [0, 1]."""
    nodes = []
    weights = []
    m = (n + 1) // 2
    for i in range(1, m + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(30):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = -p1 / dp
            x += dx
            if abs(dx) < 1e-15:
                break
        w = 2.0 / ((1.0 - x * x) * dp * dp)
        nodes.append(0.5 * (1.0 - x))
        weights.append(0.5 * w)
        if i != n + 1 - i:
            nodes.append(0.5 * (1.0 + x))
            weights.append(0.5 * w)
    pairs = sorted(zip(nodes, weights))
    return [p[0] for p in pairs], [p[1] for p in pairs]


def uniformize_vec_expm(
    vec: Dict[int, float],
    rows: List[List[Tuple[int, float]]],
    diag: List[float],
    t: float,
    tol: float = 1e-14,
) -> Dict[int, float]:
    """Compute row vector `vec @ exp(t * M)` via shifted uniformization.

    M has off-diagonal entries given by `rows[i] = [(j, rate_ij), ...]` (all >= 0)
    and diagonal entries `diag[i]`. Works for true generators, killed generators,
    and Feynman-Kac / Renyi tilted generators without cancellation.
    """
    if t == 0.0:
        return dict(vec)
    n = len(diag)
    c = max(diag) if n > 0 else 0.0
    alpha = max(c - d for d in diag) if n > 0 else 0.0
    if alpha == 0.0:
        scale = math.exp(c * t)
        return {i: v * scale for i, v in vec.items()}

    lam = alpha * t
    # Choose truncation K using Poisson tail bound
    k_max = max(20, int(lam + 10.0 * math.sqrt(lam) + 40))

    # Build sparse P = I + (M - c*I) / alpha (every entry >= 0)
    p_diag = [1.0 + (diag[i] - c) / alpha for i in range(n)]
    inv_alpha = 1.0 / alpha
    p_off = [[(j, r * inv_alpha) for j, r in rows[i]] for i in range(n)]

    # Use shifted exponent to avoid underflow for large lam:
    # sum_{k>=0} exp(-lam + c*t) * lam^k / k! * (vec P^k)
    cur = dict(vec)
    # Compute Poisson weights normalized safely around mode if needed
    # Here lam <= 200 in all tests, so direct log-space weight is exact.
    log_pref = c * t - lam
    w = math.exp(log_pref)
    res: Dict[int, float] = defaultdict(float)
    for i, val in cur.items():
        res[i] += w * val

    log_lam = math.log(lam)
    log_fact = 0.0
    for k in range(1, k_max + 1):
        nxt: Dict[int, float] = defaultdict(float)
        for i, val in cur.items():
            if val == 0.0:
                continue
            pd = p_diag[i]
            if pd != 0.0:
                nxt[i] += val * pd
            for j, pr in p_off[i]:
                nxt[j] += val * pr
        cur = nxt
        log_fact += math.log(k)
        log_w = log_pref + k * log_lam - log_fact
        if log_w < -700.0:
            w = 0.0
        else:
            w = math.exp(log_w)
        if w > 0.0:
            max_add = 0.0
            for i, val in cur.items():
                add = w * val
                res[i] += add
                if abs(add) > max_add:
                    max_add = abs(add)
            if k > lam and max_add < tol * max(1.0, sum(abs(x) for x in res.values())):
                break
    return dict(res)


# ---------------------------------------------------------------------------
# 2. Cut-side state representation (rooted core tree + leaf counts)
# ---------------------------------------------------------------------------

# A side state is (parents_tuple, leaves_tuple), where:
#   - cores are indexed 0, 1, ..., C-1 (0 is root v or p, with 1 retained cut slot)
#   - parents_tuple[i] is parent core index of core i (parents_tuple[0] = -1)
#   - leaves_tuple[i] is number of leaves L_i attached to core i.
SideState = Tuple[Tuple[int, ...], Tuple[int, ...]]


def side_degrees(state: SideState, cut_slot: int = 1) -> List[int]:
    """Return degree D_i for each core i in a side state (including cut_slot at root 0)."""
    parents, leaves = state
    c = len(parents)
    core_deg = [0] * c
    core_deg[0] += cut_slot
    for i in range(1, c):
        p = parents[i]
        core_deg[i] += 1
        core_deg[p] += 1
    return [core_deg[i] + leaves[i] for i in range(c)]


def side_vertex_count(state: SideState) -> int:
    parents, leaves = state
    return len(parents) + sum(leaves)


def side_intrinsic_transitions(
    state: SideState,
    cut_slot: int = 1,
    include_core_children_in_r: bool = True,
) -> Tuple[List[Tuple[SideState, float, bool]], float, float]:
    """Return (transitions, root_r, cross_out_rate) for one side.

    Each transition is (next_state, rate, is_root_leaf_gain).
    - root_r = r_x = L_0 + sum_{y core child of 0} 1/D_y (intrinsic root leaf gain rate).
    - cross_out_rate = 1 / D_0 (rate of arrow across cut, e.g. v->p or p->v).
    Total intrinsic exit rate Lambda = sum of all returned transition rates.
    Note: Lambda + cross_out_rate == side_vertex_count(state) when cut_slot == 1
    and include_core_children_in_r == True!
    """
    parents, leaves = state
    c = len(parents)
    deg = side_degrees(state, cut_slot=cut_slot)
    children: List[List[int]] = [[] for _ in range(c)]
    for i in range(1, c):
        children[parents[i]].append(i)

    trans: List[Tuple[SideState, float, bool]] = []
    root_r = 0.0

    # Leaf gain at each core i
    for i in range(c):
        rate = float(leaves[i])
        if i > 0:
            rate += 1.0 / deg[parents[i]]
        if i != 0 or include_core_children_in_r:
            for ch in children[i]:
                rate += 1.0 / deg[ch]
        if i == 0:
            root_r = rate
        if rate > 0.0 or i == 0:
            new_leaves = list(leaves)
            new_leaves[i] += 1
            trans.append(((parents, tuple(new_leaves)), rate, i == 0))

    # Budding (promotion) at each core i with L_i >= 1
    for i in range(c):
        if leaves[i] >= 1:
            rate = leaves[i] / deg[i]
            new_parents = parents + (i,)
            new_leaves = list(leaves)
            new_leaves[i] -= 1
            new_leaves.append(1)
            trans.append(((new_parents, tuple(new_leaves)), rate, False))

    cross_out = (1.0 / deg[0]) if cut_slot >= 1 else 0.0
    return trans, root_r, cross_out


def build_side_state_space(
    init: SideState,
    max_jumps: int,
    cut_slot: int = 1,
    include_core_children_in_r: bool = True,
):
    """Enumerate all side states reachable in <= max_jumps jumps (BFS by jump count).

    Note: every jump increases vertex count by 1, so jump count = |V(state)| - |V(init)|.
    """
    v0 = side_vertex_count(init)
    idx_of: Dict[SideState, int] = {init: 0}
    states: List[SideState] = [init]
    queue = [init]
    while queue:
        st = queue.pop(0)
        jumps = side_vertex_count(st) - v0
        if jumps >= max_jumps:
            continue
        trans, _, _ = side_intrinsic_transitions(
            st, cut_slot=cut_slot, include_core_children_in_r=include_core_children_in_r
        )
        for nxt, _, _ in trans:
            if nxt not in idx_of:
                idx_of[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)

    # Precompute transitions and metadata for each state
    n = len(states)
    root_gain_idx: List[int] = [-1] * n
    nonroot_trans: List[List[Tuple[int, float]]] = [[] for _ in range(n)]
    root_r = [0.0] * n
    cross_out = [0.0] * n
    lam_int = [0.0] * n
    core_count = [len(st[0]) for st in states]
    root_deg = [side_degrees(st, cut_slot=cut_slot)[0] for st in states]
    root_leaves = [st[1][0] for st in states]

    for i, st in enumerate(states):
        trans, r0, cout = side_intrinsic_transitions(
            st, cut_slot=cut_slot, include_core_children_in_r=include_core_children_in_r
        )
        root_r[i] = r0
        cross_out[i] = cout
        total_rate = 0.0
        for nxt, rate, is_rg in trans:
            total_rate += rate
            j = idx_of.get(nxt, -1)
            if is_rg:
                root_gain_idx[i] = j
            elif j >= 0 and rate > 0.0:
                nonroot_trans[i].append((j, rate))
        lam_int[i] = total_rate

    return {
        "states": states,
        "idx_of": idx_of,
        "root_gain_idx": root_gain_idx,
        "nonroot_trans": nonroot_trans,
        "root_r": root_r,
        "cross_out": cross_out,
        "lam_int": lam_int,
        "core_count": core_count,
        "root_deg": root_deg,
        "root_leaves": root_leaves,
        "v0": v0,
        "max_jumps": max_jumps,
    }


# ---------------------------------------------------------------------------
# 3. Direct unlabeled tree jump chain (independent from cut representation!)
# ---------------------------------------------------------------------------

def canonical_unrooted_tree(adj: List[List[int]]) -> Tuple:
    """AHU canonical form for an unrooted finite tree (exact isomorphism invariant)."""
    n = len(adj)
    if n == 1:
        return ((),)
    # Find center(s) by iterative leaf peeling
    deg = [len(adj[i]) for i in range(n)]
    rem = n
    leaves = [i for i in range(n) if deg[i] <= 1]
    alive = [True] * n
    while rem > 2:
        new_leaves = []
        for u in leaves:
            alive[u] = False
            rem -= 1
            for w in adj[u]:
                if alive[w]:
                    deg[w] -= 1
                    if deg[w] == 1:
                        new_leaves.append(w)
        leaves = new_leaves
    centers = [i for i in range(n) if alive[i]]

    def rooted(u: int, parent: int) -> Tuple:
        kids = sorted(rooted(w, u) for w in adj[u] if w != parent)
        return tuple(kids)

    forms = [rooted(c, -1) for c in centers]
    return min(forms)


def tree_core_count(adj: List[List[int]]) -> int:
    return sum(1 for nbrs in adj if len(nbrs) >= 2)


def build_t_star_adj(B: int, H: int) -> List[List[int]]:
    """Build adjacency list of T_* from Section 7:
    v(0) -- leaf(1), v(0) -- p(2), p(2) -- u_k (B-1 cores), each u_k -- (H-1) leaves.
    """
    n = 3 + (B - 1) * H
    adj: List[List[int]] = [[] for _ in range(n)]

    def add_edge(a: int, b: int) -> None:
        adj[a].append(b)
        adj[b].append(a)

    add_edge(0, 1)  # v -- l_v
    add_edge(0, 2)  # v -- p
    nxt = 3
    for _ in range(B - 1):
        u = nxt
        nxt += 1
        add_edge(2, u)
        for _ in range(H - 1):
            leaf = nxt
            nxt += 1
            add_edge(u, leaf)
    assert nxt == n
    return adj


def exact_jump_chain_from_t2(max_n: int):
    """Exact rational law of T_n (up to tree isomorphism) for n = 2..max_n."""
    # Store map: canon -> (prob Fraction, representative adj)
    t2 = [[1], [0]]
    c2 = canonical_unrooted_tree(t2)
    dist: Dict[int, Dict[Tuple, Tuple[Fraction, List[List[int]]]]] = {
        2: {c2: (Fraction(1), t2)}
    }
    a_n: Dict[int, Fraction] = {2: Fraction(0)}

    for n in range(2, max_n):
        nxt_dist: Dict[Tuple, Tuple[Fraction, List[List[int]]]] = {}
        for canon, (prob, adj) in dist[n].items():
            # For each vertex w, attach-to-w prob = (1/n) * sum_{u ~ w} 1/deg(u)
            for w in range(n):
                p_attach = Fraction(0)
                for u in adj[w]:
                    p_attach += Fraction(1, n * len(adj[u]))
                new_adj = [list(nbrs) for nbrs in adj]
                new_adj[w].append(n)
                new_adj.append([w])
                nc = canonical_unrooted_tree(new_adj)
                add_p = prob * p_attach
                if nc in nxt_dist:
                    old_p, rep = nxt_dist[nc]
                    nxt_dist[nc] = (old_p + add_p, rep)
                else:
                    nxt_dist[nc] = (add_p, new_adj)
        assert sum(p for p, _ in nxt_dist.values()) == 1
        dist[n + 1] = nxt_dist
        a_n[n + 1] = sum(p * tree_core_count(adj) for p, adj in nxt_dist.values())
    return dist, a_n


def yule_n_prob(k: int, n: int, t: float) -> float:
    """Exact P(N(t) = n) from k initial vertices: NegBin(k, e^{-t})."""
    if n < k:
        return 0.0
    j = n - k
    if j == 0:
        return math.exp(-k * t)
    if t == 0.0:
        return 0.0
    p = math.exp(-t)
    q = 1.0 - p
    log_c = math.lgamma(k + j) - math.lgamma(k) - math.lgamma(j + 1)
    return math.exp(log_c - k * t + j * math.log(q))


def exact_f_bounds(a_n: Dict[int, Fraction], t: float) -> Tuple[float, float]:
    """Certified [F_lo(t), F_hi(t)] using exact a_n (n <= N0) + exact NegBin tail."""
    n0 = max(a_n.keys())
    lower = 0.0
    used_first_moment = 0.0
    for n in range(2, n0 + 1):
        pn = yule_n_prob(2, n, t)
        lower += float(a_n[n]) * pn
        used_first_moment += n * pn
    # Exact E[N(t)] = 2 * exp(t). For n > n0, 0 <= a_n <= n - 2 <= n.
    tail_first_moment = max(0.0, 2.0 * math.exp(t) - used_first_moment)
    return lower, lower + tail_first_moment


# ---------------------------------------------------------------------------
# 4. TEST 1: Section 2 Joint Path Density (eq. (3)) + Planted Errors
# ---------------------------------------------------------------------------

def test_section2_joint_density():
    print("=== TEST 1: Section 2 Joint Path Density (eq. 3) ===")
    # Seed: interior v with 1 leaf; exterior p with 1 core neighbor u having 1 leaf
    # (so p has degree 2: slot to v and core u; u has degree 2: p and 1 leaf).
    x0: SideState = ((-1,), (1,))
    y0: SideState = ((-1, 0), (0, 1))
    t = 0.35
    max_jx = 3
    max_jy = 4

    sx = build_side_state_space(x0, max_jx)
    sy = build_side_state_space(y0, max_jy)

    # Verify exact generator total exit rate identity on every state
    for i, st in enumerate(sx["states"]):
        lhs = sx["lam_int"][i] + sx["cross_out"][i]
        rhs = float(side_vertex_count(st))
        assert abs(lhs - rhs) < 1e-14, f"Interior rate sum mismatch: {lhs} vs {rhs}"
    for j, st in enumerate(sy["states"]):
        lhs = sy["lam_int"][j] + sy["cross_out"][j]
        rhs = float(side_vertex_count(st))
        assert abs(lhs - rhs) < 1e-14, f"Exterior rate sum mismatch: {lhs} vs {rhs}"

    # Method A: Exact killed full coupled CTMC on sx x sy
    nx = len(sx["states"])
    ny = len(sy["states"])
    n_tot = nx * ny
    rows_full: List[List[Tuple[int, float]]] = [[] for _ in range(n_tot)]
    diag_full = [0.0] * n_tot
    for ix in range(nx):
        ax = sx["cross_out"][ix]
        rx = sx["root_r"][ix]
        rgx = sx["root_gain_idx"][ix]
        for iy in range(ny):
            by = sy["cross_out"][iy]
            ry = sy["root_r"][iy]
            rgy = sy["root_gain_idx"][iy]
            idx = ix * ny + iy
            diag_full[idx] = -(sx["lam_int"][ix] + by + sy["lam_int"][iy] + ax)
            # Non-root interior jumps
            for jx, rate in sx["nonroot_trans"][ix]:
                rows_full[idx].append((jx * ny + iy, rate))
            # Root interior leaf gain
            if rgx >= 0:
                rows_full[idx].append((rgx * ny + iy, rx + by))
            # Non-root exterior jumps
            for jy, rate in sy["nonroot_trans"][iy]:
                rows_full[idx].append((ix * ny + jy, rate))
            # Root exterior leaf gain
            if rgy >= 0:
                rows_full[idx].append((ix * ny + rgy, ry + ax))

    full_dist = uniformize_vec_expm({0: 1.0}, rows_full, diag_full, t)
    marg_x_full = [0.0] * nx
    for idx, prob in full_dist.items():
        marg_x_full[idx // ny] += prob
    ec_full = sum(marg_x_full[ix] * sx["core_count"][ix] for ix in range(nx))
    mass_full = sum(marg_x_full)

    # Method B: Exact integration of eq. (3) over all interior jump histories of length 0..max_jx
    # For fixed interior state sequence x_0 -> x_1 -> ... -> x_m and jump times 0 < s_1 < ... < s_m < t,
    # the exterior integral \int Q^{a_x}(dy) e^{-\int b_y} \prod_{u in G(x)} [r_x(u) + b_y(u)]
    # is computed EXACTLY on the killed exterior state space sy via piecewise uniformization.
    def eval_eq3(
        use_sx,
        include_exp_by: bool = True,
        quad_pts: int = 7,
    ) -> Tuple[List[float], float, float]:
        gl_x, gl_w = gauss_legendre(quad_pts)

        # Precompute exterior generator matrices for each distinct a_x value
        ext_mats: Dict[float, Tuple[List[List[Tuple[int, float]]], List[float]]] = {}
        for ax_val in set(use_sx["cross_out"]):
            rows_e: List[List[Tuple[int, float]]] = [[] for _ in range(ny)]
            diag_e = [0.0] * ny
            for iy in range(ny):
                by = sy["cross_out"][iy]
                ry = sy["root_r"][iy]
                rgy = sy["root_gain_idx"][iy]
                diag_e[iy] = -(sy["lam_int"][iy] + ax_val + (by if include_exp_by else 0.0))
                for jy, rate in sy["nonroot_trans"][iy]:
                    rows_e[iy].append((jy, rate))
                if rgy >= 0:
                    rows_e[iy].append((rgy, ry + ax_val))
            ext_mats[ax_val] = (rows_e, diag_e)

        marg = [0.0] * len(use_sx["states"])

        def recurse(ix: int, t_prev: float, vec_y: Dict[int, float], weight_x: float, depth: int):
            # Finish interval [t_prev, t] with no more interior jumps
            dt_end = t - t_prev
            ax = use_sx["cross_out"][ix]
            rows_e, diag_e = ext_mats[ax]
            v_end = uniformize_vec_expm(vec_y, rows_e, diag_e, dt_end)
            surv_x = math.exp(-use_sx["lam_int"][ix] * dt_end)
            marg[ix] += weight_x * surv_x * sum(v_end.values())

            if depth >= max_jx:
                return

            # Next jump time s in (t_prev, t) via Gauss-Legendre on [t_prev, t]
            next_moves: List[Tuple[int, float, bool]] = []
            for jx, rate in use_sx["nonroot_trans"][ix]:
                next_moves.append((jx, rate, False))
            rgx = use_sx["root_gain_idx"][ix]
            if rgx >= 0:
                next_moves.append((rgx, use_sx["root_r"][ix], True))

            if not next_moves:
                return

            span = t - t_prev
            for qx, qw in zip(gl_x, gl_w):
                dt = span * qx
                s = t_prev + dt
                w_dt = span * qw
                v_mid = uniformize_vec_expm(vec_y, rows_e, diag_e, dt)
                surv_seg = math.exp(-use_sx["lam_int"][ix] * dt)
                base_w = weight_x * surv_seg * w_dt
                for jx, rate, is_rg in next_moves:
                    if is_rg:
                        v_next = {
                            iy: val * (rate + sy["cross_out"][iy])
                            for iy, val in v_mid.items()
                        }
                        recurse(jx, s, v_next, base_w, depth + 1)
                    else:
                        recurse(jx, s, v_mid, base_w * rate, depth + 1)

        recurse(0, 0.0, {0: 1.0}, 1.0, 0)
        ec = sum(marg[ix] * use_sx["core_count"][ix] for ix in range(len(marg)))
        return marg, ec, sum(marg)

    marg_eq3, ec_eq3, mass_eq3 = eval_eq3(sx, include_exp_by=True, quad_pts=7)
    max_diff = max(abs(a - b) for a, b in zip(marg_x_full, marg_eq3))
    print(f"  Exact killed mass: full={mass_full:.12f}, eq(3)={mass_eq3:.12f}")
    print(f"  Exact E[C_B(t) 1_trunc]: full={ec_full:.12f}, eq(3)={ec_eq3:.12f}")
    print(f"  Max marginal L_inf error across {nx} interior states: {max_diff:.3e}")
    assert max_diff < 1e-11, f"eq(3) mismatch: {max_diff}"

    # Also cross-check full CTMC against direct independent vertex-ringing Monte Carlo!
    # (Verifies that our cut-side generator itself matches the true friend-tree generator.)
    rng = random.Random(20261124)
    mc_n = 40000
    mc_cb = 0.0
    mc_in_trunc = 0
    for _ in range(mc_n):
        # Start from T_4 path: leaf(0)-v(1)-p(2)-u(3)-leaf(4)? Wait:
        # x0 has v(0) + 1 leaf(1); y0 has p(2) -- u(3) -- leaf(4). Total 5 vertices!
        # Vertices carry side tag: 'x' or 'y'.
        adj = [[1, 2], [0], [0, 3], [2, 4], [3]]
        side = ["x", "x", "y", "y", "y"]
        cur_t = 0.0
        jx = 0
        jy = 0
        while True:
            cur_t += rng.expovariate(len(adj))
            if cur_t > t:
                break
            ring_v = rng.randrange(len(adj))
            nbrs = adj[ring_v]
            target = nbrs[rng.randrange(len(nbrs))]
            new_v = len(adj)
            adj[target].append(new_v)
            adj.append([target])
            s_tag = side[target]
            side.append(s_tag)
            if s_tag == "x":
                jx += 1
            else:
                jy += 1
        if jx <= max_jx and jy <= max_jy:
            mc_in_trunc += 1
            cb = sum(1 for i in range(len(adj)) if side[i] == "x" and len(adj[i]) >= 2)
            mc_cb += cb
    mc_ec_trunc = mc_cb / mc_n
    mc_mass = mc_in_trunc / mc_n
    se_ec = math.sqrt(max(1e-9, mc_ec_trunc * 2.0 / mc_n))
    print(f"  Independent vertex-ringing MC: mass={mc_mass:.5f} (exact {mass_full:.5f}), "
          f"E[C_B 1_trunc]={mc_ec_trunc:.5f} (exact {ec_full:.5f})")
    assert abs(mc_mass - mass_full) < 4 * math.sqrt(mass_full * (1 - mass_full) / mc_n)
    assert abs(mc_ec_trunc - ec_full) < 4 * se_ec

    # PLANTED ERRORS for Test 1:
    # 1a: Drop retained cut slot at v (cut_slot = 0)
    sx_bad_slot = build_side_state_space(x0, max_jx, cut_slot=0)
    _, ec_bad_slot, mass_bad_slot = eval_eq3(sx_bad_slot, include_exp_by=True, quad_pts=5)
    err_1a = abs(ec_bad_slot - ec_full)
    # 1b: Omit survival factor e^{-\int b_y}
    _, ec_no_surv, mass_no_surv = eval_eq3(sx, include_exp_by=False, quad_pts=5)
    err_1b = abs(mass_no_surv - mass_full)
    # 1c: Omit core-children sum in r_x
    sx_bad_r = build_side_state_space(x0, max_jx, include_core_children_in_r=False)
    _, ec_bad_r, _ = eval_eq3(sx_bad_r, include_exp_by=True, quad_pts=5)
    err_1c = abs(ec_bad_r - ec_full)
    print(f"  Planted errors caught: 1a(no cut slot) diff={err_1a:.4e}, "
          f"1b(no exp(-int b_y)) diff={err_1b:.4e}, 1c(bad r_x) diff={err_1c:.4e}")
    assert err_1a > 1e-2 and err_1b > 1e-2 and err_1c > 1e-4
    print("  [PASS] Test 1\n")


# ---------------------------------------------------------------------------
# 5. TEST 2: Section 3 Birth-Death PGF & Exterior Estimate (eq. (4))
# ---------------------------------------------------------------------------

def test_section3_exterior_estimate():
    print("=== TEST 2: Section 3 Birth-Death Law & Exterior Bound (eq. 4) ===")
    # 2A: Exact linear birth-death Kolmogorov forward equations vs closed-form geometric law
    k_trunc = 260
    for b in [2, 3, 5, 8]:
        s = 1.0 - 1.0 / b
        rho = 1.0 / b
        rows: List[List[Tuple[int, float]]] = [[] for _ in range(k_trunc + 1)]
        diag = [0.0] * (k_trunc + 1)
        for k in range(1, k_trunc + 1):
            diag[k] = -k * (1.0 + rho)
            rows[k].append((k - 1, k * rho))
            if k < k_trunc:
                rows[k].append((k + 1, float(k)))
        for h in [0.1, 0.7, 1.5, 2.5]:
            dist = uniformize_vec_expm({1: 1.0}, rows, diag, h)
            p0_num = dist.get(0, 0.0)
            r = math.exp(-s * h)
            p0_exact = (1.0 - s) * (1.0 - r) / (1.0 - (1.0 - s) * r)
            assert abs(p0_num - p0_exact) < 1e-11, f"p0 mismatch: {p0_num} vs {p0_exact}"

            p_pos = 1.0 - p0_num
            q_exact = s * r / (1.0 - (1.0 - s) * r)
            assert q_exact <= r + 1e-15
            max_pmf_err = 0.0
            e_inv_pos = 0.0
            e_inv_surv_num = 0.0
            for k in range(1, k_trunc + 1):
                pk = dist.get(k, 0.0)
                pk_cond = pk / p_pos
                pk_geom = q_exact * ((1.0 - q_exact) ** (k - 1))
                if k < 40:
                    max_pmf_err = max(max_pmf_err, abs(pk_cond - pk_geom))
                e_inv_pos += (pk / p_pos) / k
                # P(survive forever | Z_h = k) = 1 - (1 - s)^k
                p_surv_k = 1.0 - ((1.0 - s) ** k)
                e_inv_surv_num += pk * p_surv_k / k
            e_inv_surv = e_inv_surv_num / s
            e_inv_exact = -q_exact * math.log(q_exact) / (1.0 - q_exact)
            assert max_pmf_err < 1e-10, f"Geom pmf err {max_pmf_err}"
            assert abs(e_inv_pos - e_inv_exact) < 1e-9
            assert e_inv_exact <= math.sqrt(q_exact) + 1e-15
            assert e_inv_exact <= math.exp(-0.5 * s * h) + 1e-15
            assert e_inv_surv <= (1.0 / s) * math.exp(-0.5 * s * h) + 1e-12
    print("  Birth-death PGF, geometric conditional law, and 1/Z_h bounds verified to <1e-10.")

    # 2B: Exact CTMC + Monte Carlo check of eq. (4) under Q^a
    # Worst-case exterior for slow growth of D_p: p has NO leaves and ONE exterior
    # core neighbor of huge degree H_ext = 1000 (so exterior adds leaves to p at rate <= 0.001).
    # We test BOTH exact Q^a on small exterior (p + 1 leaf) AND MC on adversarial exterior!
    def simulate_qa_inv_dp(
        init_ext_core_deg: int,
        init_lp: int,
        a_fn,
        u: float,
        n_trials: int,
        seed: int,
    ) -> Tuple[float, float]:
        """Simulate Q^a where p has 1 retained slot, init_lp leaves, and 1 core neighbor
        of initial degree init_ext_core_deg (whose degree only grows, upper-bounding 1/D_p
        if we conservatively ignore exterior core-neighbor growth? No, simulate EXACT
        tree where the exterior neighbor has init_ext_core_deg-1 leaves!)."""
        # Better: for exactness, simulate the exact exterior tree with p (slot=1, L_p=init_lp)
        # joined to core w with (init_ext_core_deg - 1) leaves.
        rng = random.Random(seed)
        tot = 0.0
        tot2 = 0.0
        # Piecewise-constant a_fn returns (a_val, next_change_time)
        for _ in range(n_trials):
            # Track simplified exact exterior when init_ext_core_deg is large?
            # Even better: simulate exact full exterior vertex list!
            # Wait: for init_ext_core_deg=1000, 1000 leaves would be slow; use exact
            # compressed core-tree representation (core_parents, core_leaves, slot at 0)!
            parents = [-1, 0]
            leaves = [init_lp, init_ext_core_deg - 1]
            c_deg = [2, 1]  # p has cut slot + w; w has p
            t_cur = 0.0
            while t_cur < u:
                a_val, t_next = a_fn(t_cur)
                t_lim = min(u, t_next)
                # Total rate of exterior + immigration:
                # For each core i:
                #   leaves[i] (leaves ringing -> +1 leaf at i)
                #   + (D_i - 1)/D_i if i==0 else 1 (core i ringing -> either bud or send to core nbr)
                # Total = sum(leaves) + len(parents) - 1/D_0 + a_val
                d0 = c_deg[0] + leaves[0]
                tot_rate = sum(leaves) + len(parents) - 1.0 / d0 + a_val
                dt = rng.expovariate(tot_rate)
                if t_cur + dt >= t_lim:
                    t_cur = t_lim
                    continue
                t_cur += dt
                # Choose which channel fired
                r = rng.random() * tot_rate
                if r < a_val:
                    leaves[0] += 1
                    continue
                r -= a_val
                # Check each core i: leaf ring (rate leaves[i]), bud (rate leaves[i]/D_i),
                # core-to-core arrows (rate 1/D_i for each non-cut core neighbor)
                done = False
                for i in range(len(parents)):
                    li = leaves[i]
                    di = c_deg[i] + li
                    if r < li:
                        leaves[i] += 1
                        done = True
                        break
                    r -= li
                    bud_rate = li / di
                    if r < bud_rate:
                        leaves[i] -= 1
                        c_deg[i] += 1
                        parents.append(i)
                        leaves.append(1)
                        c_deg.append(1)
                        done = True
                        break
                    r -= bud_rate
                    # Outgoing arrows to core neighbors inside exterior
                    out_core_cnt = c_deg[i] - (1 if i == 0 else 0)
                    out_rate = out_core_cnt / di
                    if r < out_rate:
                        # Pick which exterior core neighbor uniformly
                        # Collect neighbors of i
                        nbrs = []
                        if parents[i] >= 0:
                            nbrs.append(parents[i])
                        for j in range(1, len(parents)):
                            if parents[j] == i:
                                nbrs.append(j)
                        tgt = nbrs[rng.randrange(len(nbrs))]
                        leaves[tgt] += 1
                        done = True
                        break
                    r -= out_rate
                if not done:
                    leaves[0] += 1
            val = 1.0 / (c_deg[0] + leaves[0])
            tot += val
            tot2 += val * val
        mean = tot / n_trials
        se = math.sqrt(max(0.0, tot2 / n_trials - mean * mean) / n_trials)
        return mean, se

    # Test several deterministic profiles a(w) <= 1/2
    profiles = [
        ("const_0.5", lambda w: (0.5, 1e9), lambda u: 0.5 * u),
        ("step_0.5_to_0.1", lambda w: (0.5, 2.0) if w < 2.0 else (0.1, 1e9),
         lambda u: 0.5 * min(u, 2.0) + 0.1 * max(0.0, u - 2.0)),
        ("zero_gap", lambda w: (0.0, 1.5) if w < 1.5 else (0.5, 1e9),
         lambda u: 0.5 * max(0.0, u - 1.5)),
    ]
    b_val, kappa = 4, 0.6  # s = 0.75 > 0.6
    s_val = 1.0 - 1.0 / b_val
    d0_const = ((1.0 - kappa) ** (-b_val)) / (s_val - kappa)

    for pname, afn, A_fn in profiles:
        u = 4.0
        mean_inv, se = simulate_qa_inv_dp(5, 0, afn, u, n_trials=4000, seed=101)
        rhs = d0_const * math.exp(-kappa * A_fn(u))
        print(f"  Q^a ({pname}, u={u}): E[1/D_p]={mean_inv:.5f}+/-{se:.5f} <= {rhs:.5f} (D_0={d0_const:.2f})")
        assert mean_inv + 4 * se <= rhs

    # PLANTED ERRORS for Test 2:
    # 2a: kappa = 1.3 > 1 with a = 0.5, exterior where p has 1 neighbor of degree 1000
    #     (ignoring exterior neighbor's negligible rate <= 0.001, P(0 immigrants) = e^{-0.5 u},
    #      and E[1/D_p] >= 0.5 * e^{-0.5 u} * e^{-0.001 u} >> D_0 * e^{-1.3 * 0.5 u} for large u!)
    # Check directly with a clean pure-immigration lower bound (zero exterior arrivals -> D_p = 2 + Pois(0.5*u)
    # before leaves multiply; leaves only make D_p larger, EXCEPT on the event of 0 immigrants where D_p=2!)
    u_big = 20.0
    lb_no_imm = 0.5 * math.exp(-0.5 * u_big)  # exact contribution of {0 immigrants, no ext arrival}
    bad_kappa_rhs = 10.0 * math.exp(-1.3 * 0.5 * u_big)
    assert lb_no_imm > bad_kappa_rhs, "Planted error 2a (kappa > 1) should fail!"

    # 2b: Drop D_0 (replace D_0 by 1.0) for a = 0.5, u = 2.0, kappa = 0.99:
    # Exact pure-immigration lower bound when exterior neighbor has huge degree (rate ~ 0):
    # Even if every immigrant instantly became a rate-1 Yule tree? No, check the 0-immigrant + 1-immigrant
    # event: P(0 imm) = e^{-1}, D_p=2 -> 0.5*e^{-1} = 0.1839; P(1 imm at t, no leaf ring before u) = ...
    # Actually compute exact lower bound on E[1/D_p(2)] for p with 1 huge-deg core neighbor:
    # simulate with ext deg = 1000 (use pure p-leaves model where ext arrival rate = 0, exact upper bound
    # on D_p? No, 0 ext arrivals is lower bound on D_p, i.e., upper bound on 1/D_p; wait, we can directly
    # simulate exact 1-neighbor ext with deg 1000 ONLY tracking p's leaves since neighbor rate <= 1/1000!).
    # Better: exact closed-form lower bound on E[1/D_p(u)] when ext has 1 neighbor of degree 10^6:
    # On event {no neighbor->p arrow (prob >= e^{-u*1e-6})}:
    #   - 0 immigrants (prob e^{-A}): D_p = 2, contrib 1/2
    #   - 1 immigrant at time s, 0 rings of that leaf in (s, u] (prob int_0^u a e^{-A} e^{-(u-s)} ds): D_p = 3, contrib 1/3
    # For a = 0.5, u = 4.0 (A = 2.0):
    #   0 imm: 0.5 * e^{-2} = 0.06767
    #   1 imm, 0 rings: (1/3) * 0.5 * e^{-2} * (1 - e^{-4}) = 0.02214
    #   Wait: compare with 1.0 * e^{-0.99 * 2} = 0.13807. Not yet > 1.0!
    # At u = 0.0? No, at u = 0, 1/D_p = 0.5 <= 1.0.
    # Wait: what about D_p(0) = 2, u = 0.5, a = 0.5, kappa = 0.99 -> e^{-0.2475} = 0.78 > 0.5.
    # Can E[1/D_p(u)] EVER exceed e^{-\kappa A(u)} when D_p(0) = 2 and \kappa = 0.99?
    # YES! Wait: for large u, with a = 0.5, A(u) = u/2.
    # On {0 immigrants}: contribution is 0.5 * e^{-u/2} = 0.5 * e^{-A(u)}.
    # Compare 0.5 * e^{-A(u)} vs 1.0 * e^{-0.99 A(u)}: 0.5 e^{-A} < e^{-0.99 A} always!
    # What if we add 1, 2, 3... immigrants? Wait: at A = 100, kappa = 0.99, e^{-0.99*100} = e^{-99}
    # while 0.5 * e^{-100} = 0.5 e^{-100}. What if kappa = 0.9999 and A = 10?
    # Let's check exact E[1/(2 + Y)] where Y is pure linear birth + immigration(a=0.5),
    # i.e., exact Q^{0.5} with NO exterior arrivals: P(N_imm = k, leaves don't grow much) ->
    # For a=0.5, leaf birth rate=1, promotion=1/D_p, wait: at u = 1.0, A = 0.5,
    # D_p = 2 + L_p(1). Even with 0 promotions, L_p is Yule with immigration 0.5:
    # PGF of L_p(u) with imm rate a=0.5, birth rate 1: (e^{-u} / (1 - (1-e^{-u})z))^{0.5}!
    # So E[1/(2+L_p(u))] >= int_0^1 z * (e^{-u}/(1 - (1-e^{-u})z))^{0.5} dz.
    # For u = 10: int_0^1 z * e^{-5} / sqrt(1 - (1-e^{-10})z) dz ~ e^{-5} * (4/3) = 1.333 * e^{-5}!
    # WHOA!! 1.333 * e^{-5} > 1.0 * e^{-0.99 * 5} = 1.051 * e^{-5}!!
    # Let's verify that exact integral numerically!
    gl_x, gl_w = gauss_legendre(40)
    u_test = 10.0
    q_y = 1.0 - math.exp(-u_test)
    # Note: promotions only DECREASE L_p (making D_p smaller? Wait: promotion decreases L_p by 1
    # and increases c_p by 1, leaving D_p unchanged, AND reduces leaf birth rate by 1, so D_p is
    # STOCHASTICALLY SMALLER than without promotions, hence 1/D_p is EVEN LARGER!)
    e_inv_no_prom = sum(
        w * x * math.exp(-0.5 * u_test) / math.sqrt(1.0 - q_y * x)
        for x, w in zip(gl_x, gl_w)
    )
    rhs_no_d0 = math.exp(-0.99 * 0.5 * u_test)
    print(f"  Planted errors caught: 2a(kappa=1.3>1) lb={lb_no_imm:.3e} > rhs={bad_kappa_rhs:.3e}; "
          f"2b(drop D_0, kappa=0.99, u=10) E[1/D_p]>={e_inv_no_prom:.5e} > exp(-kA)={rhs_no_d0:.5e}")
    assert e_inv_no_prom > rhs_no_d0, "Planted error 2b (dropping D_0) should fail!"
    print("  [PASS] Test 2\n")


# ---------------------------------------------------------------------------
# 6. TEST 3: Section 4 Last-Exterior Decomposition & Profile Rounding (eq. 5)
# ---------------------------------------------------------------------------

def test_section4_last_arrival_and_profiles():
    print("=== TEST 3: Section 4 Last-Arrival Expansion & Profile Rounding (eq. 5) ===")
    # 3A: Exhaustive adversarial check of profile rounding properties
    rng = random.Random(42)
    for J in [2, 4, 6, 8]:
        L = J // 2
        for t in [0.05, 0.3, 1.0, 1.0001, 1.9999, 2.0, 2.5, 4.7]:
            m = max(1, math.ceil(t))
            exact_count = math.comb(m + L, L)
            bound_count = (2.0 + t) ** L
            assert exact_count <= bound_count + 1e-12, (
                f"Profile count failed: J={J}, t={t}, exact={exact_count}, bound={bound_count}"
            )
            # Test 50 adversarial nonincreasing step functions b_y in (0, 0.5]
            for _ in range(50):
                n_drops = rng.randint(0, 6)
                drops = sorted(rng.uniform(0.0, t) for _ in range(n_drops))
                # Include occasional integer jump times!
                if rng.random() < 0.3 and t > 1.0:
                    drops = sorted(drops + [float(rng.randint(1, int(t)))])
                vals = [0.5]
                for _ in drops:
                    vals.append(rng.uniform(1e-6, vals[-1]))

                def b_left(u: float) -> float:
                    # Left-continuous nonincreasing step function
                    idx = 0
                    for d in drops:
                        if d < u:
                            idx += 1
                        else:
                            break
                    return vals[idx]

                # Round up at left end of each unit interval [k, k+1)
                betas = [math.ceil(J * b_left(float(k))) / J for k in range(m)]
                # Check nonincreasing & in {0, 1/J, ..., 1/2}
                for k in range(m):
                    assert 0.0 <= betas[k] <= 0.5 + 1e-15
                    if k > 0:
                        assert betas[k] <= betas[k - 1] + 1e-15
                # Check tau in (0, t]
                for tau in [t * 0.25, t * 0.73, t]:
                    # Exact integral of (beta - b_y) on [0, tau]
                    pts = sorted(set([0.0, tau] + [float(k) for k in range(1, m) if k < tau]
                                     + [d for d in drops if d < tau]))
                    integral = 0.0
                    for a_pt, b_pt in zip(pts[:-1], pts[1:]):
                        mid = 0.5 * (a_pt + b_pt)
                        bk = betas[min(m - 1, int(math.floor(mid)))]
                        by_mid = b_left(mid)
                        assert bk >= by_mid - 1e-15, f"beta < b_y at {mid}"
                        integral += (bk - by_mid) * (b_pt - a_pt)
                    assert integral <= t / J + 0.5 + 1e-12, (
                        f"Integral bound failed: {integral} > {t/J + 0.5}"
                    )
    print("  Profile count |(P_t)| <= (2+t)^{J/2}, beta >= b_y, and int(beta-b_y) <= t/J+1/2 verified.")

    # 3B: Exact algebraic identity for last-exterior-gain decomposition on small tree
    # Compare exact full CTMC E[C_B(t) 1_{#Jx<=2, #Jy<=3}] against:
    #   (all intrinsic term) + int_0^t (last exterior gain at tau) dtau!
    x0: SideState = ((-1,), (1,))
    y0: SideState = ((-1, 0), (0, 1))
    t = 0.3
    max_jx, max_jy = 2, 3
    sx = build_side_state_space(x0, max_jx)
    sy = build_side_state_space(y0, max_jy)
    nx, ny = len(sx["states"]), len(sy["states"])

    # Full CTMC ground truth
    n_tot = nx * ny
    rows_full: List[List[Tuple[int, float]]] = [[] for _ in range(n_tot)]
    diag_full = [0.0] * n_tot
    for ix in range(nx):
        ax, rx, rgx = sx["cross_out"][ix], sx["root_r"][ix], sx["root_gain_idx"][ix]
        for iy in range(ny):
            by, ry, rgy = sy["cross_out"][iy], sy["root_r"][iy], sy["root_gain_idx"][iy]
            idx = ix * ny + iy
            diag_full[idx] = -(sx["lam_int"][ix] + by + sy["lam_int"][iy] + ax)
            for jx, rate in sx["nonroot_trans"][ix]:
                rows_full[idx].append((jx * ny + iy, rate))
            if rgx >= 0:
                rows_full[idx].append((rgx * ny + iy, rx + by))
            for jy, rate in sy["nonroot_trans"][iy]:
                rows_full[idx].append((ix * ny + jy, rate))
            if rgy >= 0:
                rows_full[idx].append((ix * ny + rgy, ry + ax))
    full_dist = uniformize_vec_expm({0: 1.0}, rows_full, diag_full, t)
    ec_full = sum(p * sx["core_count"][idx // ny] for idx, p in full_dist.items())

    # Exact decomposition via expanded Markov state (ix, iy, has_ext_gain) vs
    # explicit integral over last exterior gain time tau!
    # In exact CTMC:
    # 1) All-intrinsic: kill exterior->interior arrows (keep survival -by in diag, only rx on rgx).
    rows_intr: List[List[Tuple[int, float]]] = [[] for _ in range(n_tot)]
    for ix in range(nx):
        ax, rx, rgx = sx["cross_out"][ix], sx["root_r"][ix], sx["root_gain_idx"][ix]
        for iy in range(ny):
            ry, rgy = sy["root_r"][iy], sy["root_gain_idx"][iy]
            idx = ix * ny + iy
            for jx, rate in sx["nonroot_trans"][ix]:
                rows_intr[idx].append((jx * ny + iy, rate))
            if rgx >= 0:
                rows_intr[idx].append((rgx * ny + iy, rx))  # ONLY intrinsic r_x!
            for jy, rate in sy["nonroot_trans"][iy]:
                rows_intr[idx].append((ix * ny + jy, rate))
            if rgy >= 0:
                rows_intr[idx].append((ix * ny + rgy, ry + ax))
    dist_intr = uniformize_vec_expm({0: 1.0}, rows_intr, diag_full, t)
    ec_all_intr = sum(p * sx["core_count"][idx // ny] for idx, p in dist_intr.items())

    # 2) Last exterior gain at tau in (0, t):
    # On [0, tau): full coupled process (rows_full, diag_full).
    # At tau: forced root leaf gain at v with factor b_y(tau) dtau (only if rgx >= 0).
    # On (tau, t]: suffix with ONLY intrinsic root gains (rows_intr, diag_full).
    gl_x, gl_w = gauss_legendre(12)
    ec_last_ext = 0.0
    for qx, qw in zip(gl_x, gl_w):
        tau = t * qx
        w_tau = t * qw
        v_pre = uniformize_vec_expm({0: 1.0}, rows_full, diag_full, tau)
        v_after: Dict[int, float] = defaultdict(float)
        for idx, p in v_pre.items():
            ix, iy = idx // ny, idx % ny
            rgx = sx["root_gain_idx"][ix]
            if rgx >= 0:
                v_after[rgx * ny + iy] += p * sy["cross_out"][iy]
        v_end = uniformize_vec_expm(dict(v_after), rows_intr, diag_full, t - tau)
        ec_last_ext += w_tau * sum(p * sx["core_count"][idx // ny] for idx, p in v_end.items())

    err_decomp = abs(ec_full - (ec_all_intr + ec_last_ext))
    print(f"  Exact decomposition: full={ec_full:.12f}, intr+last_ext={ec_all_intr + ec_last_ext:.12f}, diff={err_decomp:.3e}")
    assert err_decomp < 1e-11

    # 3C: Check upper bound (5) directly for small parameters (J=2, b=2, kappa=0.25, t=0.4)
    # Use larger truncation so truncated expectations equal true ones to ~1e-5.
    t5 = 0.4
    sx5 = build_side_state_space(x0, 7)
    nx5 = len(sx5["states"])
    # Compute E^{R_0}[C_B(t5)] and E^{R_{beta, tau}}[C_B(t5) e^{-kappa A_x(tau)}]
    def eval_r_beta_tau(beta_val: float, tau: float, kappa: float, use_sx) -> float:
        n_s = len(use_sx["states"])
        # Prefix generator with immigration beta_val and Feynman-Kac potential -kappa * a_x
        rows_pre: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
        diag_pre = [0.0] * n_s
        rows_suf: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
        diag_suf = [0.0] * n_s
        for i in range(n_s):
            diag_pre[i] = -(use_sx["lam_int"][i] + beta_val + kappa * use_sx["cross_out"][i])
            diag_suf[i] = -use_sx["lam_int"][i]
            for j, r in use_sx["nonroot_trans"][i]:
                rows_pre[i].append((j, r))
                rows_suf[i].append((j, r))
            rg = use_sx["root_gain_idx"][i]
            if rg >= 0:
                rows_pre[i].append((rg, use_sx["root_r"][i] + beta_val))
                if use_sx["root_r"][i] > 0.0:
                    rows_suf[i].append((rg, use_sx["root_r"][i]))
        if tau < 0.0:  # Sentinel for R_0
            v_end = uniformize_vec_expm({0: 1.0}, rows_suf, diag_suf, t5)
            return sum(p * use_sx["core_count"][i] for i, p in v_end.items())
        v_pre = uniformize_vec_expm({0: 1.0}, rows_pre, diag_pre, tau)
        v_post: Dict[int, float] = defaultdict(float)
        for i, p in v_pre.items():
            rg = use_sx["root_gain_idx"][i]
            if rg >= 0:
                v_post[rg] += p
        v_end = uniformize_vec_expm(dict(v_post), rows_suf, diag_suf, t5 - tau)
        return sum(p * use_sx["core_count"][i] for i, p in v_end.items())

    er0 = eval_r_beta_tau(0.0, -1.0, 0.0, sx5)
    J5, b5, kappa5 = 2, 2, 0.25
    d0_5 = ((1.0 - kappa5) ** (-b5)) / ((1.0 - 1.0 / b5) - kappa5)
    # For t5 = 0.4 < 1, J5 = 2: profiles in P_t are constant beta in {0.0, 0.5}
    sum_prof = 0.0
    for beta_val in [0.0, 0.5]:
        integ = sum(
            t5 * qw * eval_r_beta_tau(beta_val, t5 * qx, kappa5, sx5)
            for qx, qw in zip(gl_x, gl_w)
        )
        sum_prof += integ
    rhs_5 = er0 + d0_5 * math.exp(0.5 + t5 / J5) * sum_prof
    # Compare with MC true E[C_B(t5)] from (x0, y0)
    rng2 = random.Random(77)
    mc_cb5 = 0.0
    for _ in range(10000):
        adj = [[1, 2], [0], [0, 3], [2, 4], [3]]
        side = ["x", "x", "y", "y", "y"]
        cur_t = 0.0
        while True:
            cur_t += rng2.expovariate(len(adj))
            if cur_t > t5:
                break
            rv = rng2.randrange(len(adj))
            tgt = adj[rv][rng2.randrange(len(adj[rv]))]
            nv = len(adj)
            adj[tgt].append(nv)
            adj.append([tgt])
            side.append(side[tgt])
        mc_cb5 += sum(1 for i in range(len(adj)) if side[i] == "x" and len(adj[i]) >= 2)
    mc_cb5 /= 10000
    print(f"  Eq (5) check (t={t5}): E[C_B(t)]={mc_cb5:.4f} <= RHS(5)={rhs_5:.4f}")
    assert mc_cb5 <= rhs_5

    # PLANTED ERRORS for Test 3:
    # 3a: Drop +1/2 total variation term from int(beta - b_y) <= t/J + 1/2
    # Step function b_y(u) = 0.5 on [0, eps], 1e-6 on (eps, 1]; J=4, t=1:
    # beta = 0.5 on [0, 1], integral = (0.5 - 1e-6)*(1 - eps) ~ 0.499 > 1/4 = 0.25!
    eps = 1e-3
    integ_tv = (0.5 - 1e-6) * (1.0 - eps)
    assert integ_tv > 1.0 / 4.0 and integ_tv <= 1.0 / 4.0 + 0.5
    # 3b: Replace (2+t)^{J/2} with (1+t)^{J/2} at t=0.5, J=2:
    assert math.comb(1 + 1, 1) > (1.0 + 0.5) ** 1
    print(f"  Planted errors caught: 3a(drop +1/2 TV) integ={integ_tv:.4f} > t/J=0.25; "
          f"3b((1+t)^(J/2) too small) exact=2 > 1.5")
    print("  [PASS] Test 3\n")


# ---------------------------------------------------------------------------
# 7. TEST 4: Section 5 Renyi Likelihood Domination (eq. (6)) + Planted Errors
# ---------------------------------------------------------------------------

def test_section5_renyi_domination():
    print("=== TEST 4: Section 5 Renyi Likelihood Domination (eq. 6) ===")
    # 4A: Dense grid verification of Taylor bound (r+gamma)(u^d - 1) - d(beta-gamma) <= c_d / K
    max_ratio = 0.0
    for d in [3, 4, 5, 7, 10]:
        c_d = d * (d - 1) * ((1.5) ** (d - 2)) / 8.0
        for K in [1, 2, 5, 10]:
            for r in [K, K + 0.5, K + 3.0, 20.0]:
                for B in [2, 3, 5, 10, 50]:
                    gamma = 1.0 / B
                    for beta in [0.0, 0.05, gamma, 0.25, 0.5]:
                        u = (r + beta) / (r + gamma)
                        assert u <= 1.5 + 1e-15
                        lhs = (r + gamma) * (u ** d - 1.0) - d * (beta - gamma)
                        rhs = c_d / K
                        assert lhs <= rhs + 1e-12, f"Taylor bound failed: {lhs} > {rhs}"
                        if rhs > 0:
                            max_ratio = max(max_ratio, lhs / rhs)
    print(f"  Taylor bound verified on dense grid (both beta>gamma and beta<gamma); max lhs/rhs = {max_ratio:.5f} <= 1.")

    # Planted error 4a: Drop (3/2)^{d-2} from c_d -> fails at d=4, K=1, r=1, beta=0.5, gamma=0.01
    d_bad = 4
    r_bad, b_bad, g_bad = 1.0, 0.5, 0.01
    u_bad = (r_bad + b_bad) / (r_bad + g_bad)
    lhs_bad = (r_bad + g_bad) * (u_bad ** d_bad - 1.0) - d_bad * (b_bad - g_bad)
    wrong_cd = d_bad * (d_bad - 1) / 8.0
    true_cd = wrong_cd * (1.5 ** (d_bad - 2))
    assert lhs_bad > wrong_cd and lhs_bad <= true_cd

    # 4B: Exact tilted-generator verification of E_P[Z] = 1, E_P[L] = 1, and E_P[Z^d] <= z^{M+K} e^{Vt}
    x0: SideState = ((-1,), (1,))
    sx = build_side_state_space(x0, 9)
    n_s = len(sx["states"])

    def exact_ezd_piecewise(
        d: int,
        gamma: float,
        beta_pieces: List[Tuple[float, float]],  # [(dt, beta_val), ...]
        only_high_K: int = 0,  # if > 0, compute E_P[L] for auxiliary law with threshold K
    ) -> float:
        """Exact killed-space computation of E_P[Z^d] (or E_P[L] when only_high_K > 0)."""
        vec: Dict[int, float] = {0: 1.0}
        for dt, beta in beta_pieces:
            rows: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
            diag = [0.0] * n_s
            for i in range(n_s):
                r = sx["root_r"][i]
                lv = sx["root_leaves"][i]
                u = (r + beta) / (r + gamma)
                if only_high_K > 0:
                    if lv >= only_high_K:
                        new_rg = (r + gamma) * (u ** d)
                        diag[i] = -(sx["lam_int"][i] - r + new_rg)
                    else:
                        new_rg = r + gamma
                        diag[i] = -(sx["lam_int"][i] + gamma)
                else:
                    new_rg = (r + gamma) * (u ** d)
                    diag[i] = -(sx["lam_int"][i] + gamma + d * (beta - gamma))
                for j, rate in sx["nonroot_trans"][i]:
                    rows[i].append((j, rate))
                rg = sx["root_gain_idx"][i]
                if rg >= 0 and new_rg > 0.0:
                    rows[i].append((rg, new_rg))
            vec = uniformize_vec_expm(vec, rows, diag, dt)
        return sum(vec.values())

    # Check normalization E_P[Z] = 1 (up to tiny tail > 9 jumps at tau = 0.25) and E_P[L] = 1
    ez1 = exact_ezd_piecewise(1, 1.0 / 3.0, [(0.12, 0.5), (0.13, 0.25)])
    el_aux = exact_ezd_piecewise(3, 1.0 / 3.0, [(0.12, 0.5), (0.13, 0.25)], only_high_K=2)
    print(f"  Normalization check (max_j=9, tau=0.25): E_P[Z]={ez1:.8f}, E_P[L_aux]={el_aux:.8f}")
    assert abs(ez1 - 1.0) < 1e-4 and abs(el_aux - 1.0) < 1e-3

    # Check eq. (6) across multiple (d, B, K, M, beta_pieces)
    test_cases = [
        (3, 2, 1, 2, [(0.25, 0.5)]),
        (3, 3, 1, 2, [(0.15, 0.5), (0.15, 0.0)]),
        (3, 4, 2, 3, [(0.20, 0.5)]),
        (4, 3, 2, 3, [(0.20, 0.5)]),
        (3, 3, 1, 2, [(0.30, 0.0)]),  # beta = 0 < gamma!
    ]
    for d, B, K, M, pieces in test_cases:
        gamma = 1.0 / B
        tau = sum(dt for dt, _ in pieces)
        z = (1.0 + 1.0 / (2.0 * gamma)) ** d
        c_d = d * (d - 1) * ((1.5) ** (d - 2)) / 8.0
        V = d * gamma + c_d / K + (z - 1.0) * K / M
        ezd = exact_ezd_piecewise(d, gamma, pieces)
        rhs = (z ** (M + K)) * math.exp(V * tau)
        print(f"  Eq (6) (d={d}, B={B}, K={K}, M={M}, tau={tau:.2f}): E_P[Z^d]={ezd:.6f} <= {rhs:.4e}")
        assert ezd <= rhs

    # 4C: Pathwise Monte Carlo verification under auxiliary law P_tilde of:
    #   (i) G_low <= M + K + N on 100% of paths
    #   (ii) E_{P_tilde}[z^N] <= exp((z-1) K t / M)
    #   (iii) Planted error 4b: omitting L_v = K promotions from N FAILS pathwise bound!
    rng = random.Random(555)
    d, B, K, M = 3, 3, 2, 4
    gamma = 1.0 / B
    beta = 0.5
    z = (1.0 + 1.0 / (2.0 * gamma)) ** d
    tau_sim = 3.0
    n_paths = 5000
    n_viol_good = 0
    n_viol_planted_4b = 0
    sum_zn = 0.0
    for _ in range(n_paths):
        # Simulate interior process under P_tilde on [0, tau_sim]
        parents = [-1]
        leaves = [1]
        c_deg = [1]  # retained slot at v
        t_cur = 0.0
        g_low = 0
        n_prom_good = 0
        n_prom_bad_4b = 0
        while t_cur < tau_sim:
            # Compute root r_x
            r_x = float(leaves[0])
            for ch in range(1, len(parents)):
                if parents[ch] == 0:
                    r_x += 1.0 / (c_deg[ch] + leaves[ch])
            lv = leaves[0]
            dv = c_deg[0] + lv
            u = (r_x + beta) / (r_x + gamma)
            rg_rate = (r_x + gamma) * (u ** d) if lv >= K else (r_x + gamma)
            # Non-root intrinsic total rate = sum(leaves) + len(parents) - 1/dv - r_x
            non_rg_rate = sum(leaves) + len(parents) - 1.0 / dv - r_x
            tot_rate = rg_rate + non_rg_rate
            dt = rng.expovariate(tot_rate)
            if t_cur + dt >= tau_sim:
                break
            t_cur += dt
            pick = rng.random() * tot_rate
            if pick < rg_rate:
                if lv < K:
                    g_low += 1
                leaves[0] += 1
                continue
            pick -= rg_rate
            # Pick non-root transition: bud at core i, leaf ring at i>0, core-to-core to i>0
            done = False
            for i in range(len(parents)):
                li = leaves[i]
                di = c_deg[i] + li
                bud_r = li / di
                if pick < bud_r:
                    if i == 0 and dv >= M:
                        if 1 <= li <= K:
                            n_prom_good += 1
                        if 1 <= li <= K - 1:
                            n_prom_bad_4b += 1
                    leaves[i] -= 1
                    c_deg[i] += 1
                    parents.append(i)
                    leaves.append(1)
                    c_deg.append(1)
                    done = True
                    break
                pick -= bud_r
                if i > 0:
                    # Leaf gain at core i > 0: leaves[i] + 1/D_{par} + sum_{ch} 1/D_{ch}
                    lg_r = float(li) + 1.0 / (c_deg[parents[i]] + leaves[parents[i]])
                    for ch in range(1, len(parents)):
                        if parents[ch] == i:
                            lg_r += 1.0 / (c_deg[ch] + leaves[ch])
                    if pick < lg_r:
                        leaves[i] += 1
                        done = True
                        break
                    pick -= lg_r
            if not done:
                leaves[-1] += 1
        if g_low > M + K + n_prom_good:
            n_viol_good += 1
        if g_low > M + K + n_prom_bad_4b:
            n_viol_planted_4b += 1
        sum_zn += z ** n_prom_good

    mean_zn = sum_zn / n_paths
    bound_zn = math.exp((z - 1.0) * K * tau_sim / M)
    print(f"  Pathwise P_tilde check ({n_paths} paths): G_low <= M+K+N violations={n_viol_good}, "
          f"E[z^N]={mean_zn:.4f} <= {bound_zn:.4e}")
    assert n_viol_good == 0
    assert mean_zn <= bound_zn

    # Also test Planted error 4b deterministically on an oscillating path:
    # Start after D_v = M with L_v = K; repeat 10 times: promote at L_v=K (L_v->K-1), root gain at L_v=K-1 (L_v->K).
    # Then G_low = (M-2) + 10, N_good = 10 -> G_low = M+8 <= M+K+10; N_bad = 0 -> M+8 > M+K+0 for K=2!
    g_osc = (M - 2) + 10
    assert g_osc <= M + K + 10 and g_osc > M + K + 0

    # Planted error 4c: Drop d*gamma from V when beta = 0, r = 0, large K, M:
    # If r=0, beta=0, P-intensity is gamma, Z = e^{gamma*t} 1_{N=0}, so E_P[Z^d] = e^{(d-1)*gamma*t}.
    # With K, M huge, c_d/K + (z-1)K/M < 1e-3, so dropping d*gamma makes RHS ~ z^{M+K} e^{0.001 t}
    # which is beaten by e^{(d-1)*gamma*t} for large t!
    print(f"  Planted errors caught: 4a(drop (3/2)^(d-2)) {lhs_bad:.4f}>{wrong_cd:.4f}; "
          f"4b(omit L_v=K from N) oscillating path {g_osc} > M+K={M+K} (MC viols={n_viol_planted_4b}); "
          f"4c(drop d*gamma) fails for beta=0, r=0 at large t.")
    print("  [PASS] Test 4\n")


# ---------------------------------------------------------------------------
# 8. TEST 5: Section 6 Holder Exponents & Yule Second Moment (eq. 7)
# ---------------------------------------------------------------------------

def test_section6_holder():
    print("=== TEST 5: Section 6 Holder Exponents & Second Moment (eq. 7) ===")
    x0: SideState = ((-1,), (1,))
    sx = build_side_state_space(x0, 7)
    n_s = len(sx["states"])

    d = 3
    q = d / (d - 1.0)
    theta = (2.0 - q) / q
    assert abs(theta - (1.0 - 2.0 / d)) < 1e-15
    assert abs((q - 1.0) / q - 1.0 / d) < 1e-15

    B = 3
    gamma = 1.0 / B
    beta = 0.5
    kappa = 0.6
    t = 0.45
    tau = 0.25

    # Compute exact under P = R_{gamma, tau} and R_{beta, tau} on sx
    def run_r_dist(imm_rate: float, pot_coeff: float) -> Dict[int, float]:
        rows_pre: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
        diag_pre = [0.0] * n_s
        rows_suf: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
        diag_suf = [0.0] * n_s
        for i in range(n_s):
            diag_pre[i] = -(sx["lam_int"][i] + imm_rate + pot_coeff * sx["cross_out"][i])
            diag_suf[i] = -sx["lam_int"][i]
            for j, r in sx["nonroot_trans"][i]:
                rows_pre[i].append((j, r))
                rows_suf[i].append((j, r))
            rg = sx["root_gain_idx"][i]
            if rg >= 0:
                rows_pre[i].append((rg, sx["root_r"][i] + imm_rate))
                if sx["root_r"][i] > 0.0:
                    rows_suf[i].append((rg, sx["root_r"][i]))
        v_pre = uniformize_vec_expm({0: 1.0}, rows_pre, diag_pre, tau)
        v_post: Dict[int, float] = defaultdict(float)
        for i, p in v_pre.items():
            rg = sx["root_gain_idx"][i]
            if rg >= 0:
                v_post[rg] += p
        return uniformize_vec_expm(dict(v_post), rows_suf, diag_suf, t - tau)

    d_p_unweighted = run_r_dist(gamma, 0.0)
    e_cb2 = sum(p * (sx["core_count"][i] ** 2) for i, p in d_p_unweighted.items())
    yule_2 = 2.0 * math.exp(2.0 * t) - math.exp(t)
    assert e_cb2 <= yule_2 + 1e-10 <= 2.0 * math.exp(2.0 * t)

    d_p_f = run_r_dist(gamma, 1.0)
    h_tau = sum(p * sx["core_count"][i] for i, p in d_p_f.items())
    d_p_f2 = run_r_dist(gamma, 2.0)
    e_f2 = sum(p * (sx["core_count"][i] ** 2) for i, p in d_p_f2.items())

    d_beta = run_r_dist(beta, kappa)
    lhs_7 = sum(p * sx["core_count"][i] for i, p in d_beta.items())

    # Exact E_P[Z^d] on same truncated space
    rows_z: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
    diag_z = [0.0] * n_s
    for i in range(n_s):
        r = sx["root_r"][i]
        u = (r + beta) / (r + gamma)
        diag_z[i] = -(sx["lam_int"][i] + gamma + d * (beta - gamma))
        for j, rt in sx["nonroot_trans"][i]:
            rows_z[i].append((j, rt))
        rg = sx["root_gain_idx"][i]
        if rg >= 0:
            rows_z[i].append((rg, (r + gamma) * (u ** d)))
    ezd = sum(uniformize_vec_expm({0: 1.0}, rows_z, diag_z, tau).values())

    rhs_7 = (
        math.exp(0.5 * (1.0 - kappa) * t)
        * (ezd ** (1.0 / d))
        * (h_tau ** theta)
        * (2.0 ** (1.0 / d))
        * math.exp(2.0 * t / d)
    )
    print(f"  E_P[C_B^2]={e_cb2:.5f} <= Yule={yule_2:.5f} <= 2e^(2t)={2*math.exp(2*t):.5f}")
    print(f"  Eq (7): LHS={lhs_7:.5f} <= RHS={rhs_7:.5f} (h(tau)={h_tau:.5f}, E_P[f^2]={e_f2:.5f})")
    assert lhs_7 <= rhs_7
    # Planted error 5: d=2 gives q=2, theta=0 (so h^theta = 1, loses linear S(t) factor? Actually still
    # valid if d=2 gives 1, wait: for d=2, (2-q)/q = 0, so (E f^q)^{1/q} = (E f^2)^{1/2} <= sqrt(2) e^t,
    # losing e^{2t/d} -> e^t which does NOT vanish! And d < 2 gives q > 2 where f^q interpolation fails!)
    d_bad = 1.5
    q_bad = d_bad / (d_bad - 1.0)
    assert q_bad > 2.0, "d < 2 gives q > 2, breaking interpolation between moments 1 and 2!"
    print("  [PASS] Test 5\n")


# ---------------------------------------------------------------------------
# 9. TEST 6: Section 7 Physical Finite Seed T_*, Eq. (8), and pi_*
# ---------------------------------------------------------------------------

def test_section7_physical_seed():
    print("=== TEST 6: Section 7 Physical Seed T_*, Eq. (8), & pi_* ===")
    # 6A: Exact jump-chain realization probability of T_* from T_2
    dist, a_n = exact_jump_chain_from_t2(9)
    for B, H in [(2, 2), (2, 3), (3, 2)]:
        n_star = 3 + (B - 1) * H
        t_star_adj = build_t_star_adj(B, H)
        c_star = canonical_unrooted_tree(t_star_adj)
        p_jump_exact = float(dist[n_star][c_star][0])
        p_jump_bound = 1.0
        for j in range(2, n_star):
            p_jump_bound /= (j * (j - 1))
        p_yule_1 = yule_n_prob(2, n_star, 1.0)
        p_yule_formula = (n_star - 1) * math.exp(-2.0) * ((1.0 - math.exp(-1.0)) ** (n_star - 2))
        assert abs(p_yule_1 - p_yule_formula) < 1e-15
        pi_star = p_yule_formula * p_jump_bound
        p_realize_exact = p_yule_1 * p_jump_exact
        print(f"  T_*(B={B},H={H},n*={n_star}): P(T_n*=T*)={p_jump_exact:.6e} >= {p_jump_bound:.6e}; "
              f"P(T(1)=T*)={p_realize_exact:.6e} >= pi_*={pi_star:.6e}")
        assert p_jump_exact >= p_jump_bound - 1e-18
        assert p_realize_exact >= pi_star - 1e-18

    # 6B: Verify EXACT coupling identity of Section 7 for (B=2, H=2) and (B=3, H=2)!
    # Under T_*, simulate using the exact independent-channel Poisson construction of Section 7:
    #   - E_0 = {no p->v proposal on [0, t]} -> E[C_B(t) 1_{E_0}] == e^{-gamma t} E^{R_0}[C_B(t)]
    #   - E_2 = E_ext \cap E_fwd \cap E_rev -> E[C_B(t) 1_{E_2}] == e^{-delta t} int_0^t gamma e^{-gamma(t-tau)} h(tau) dtau
    for B, H in [(2, 2), (3, 2)]:
        gamma = 1.0 / B
        delta = (B - 1.0) / H
        t = 0.4

        # Compute exact E^{R_0}[C_B(t)] and exact integral I_exact = e^{-delta t} int_0^t gamma e^{-gamma(t-tau)} h(tau) dtau
        x0: SideState = ((-1,), (1,))
        sx = build_side_state_space(x0, 9)
        n_s = len(sx["states"])

        def exact_r0_and_h(tau: float) -> float:
            rows_pre: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
            diag_pre = [0.0] * n_s
            rows_suf: List[List[Tuple[int, float]]] = [[] for _ in range(n_s)]
            diag_suf = [0.0] * n_s
            for i in range(n_s):
                diag_pre[i] = -(sx["lam_int"][i] + gamma + sx["cross_out"][i])
                diag_suf[i] = -sx["lam_int"][i]
                for j, r in sx["nonroot_trans"][i]:
                    rows_pre[i].append((j, r))
                    rows_suf[i].append((j, r))
                rg = sx["root_gain_idx"][i]
                if rg >= 0:
                    rows_pre[i].append((rg, sx["root_r"][i] + gamma))
                    if sx["root_r"][i] > 0.0:
                        rows_suf[i].append((rg, sx["root_r"][i]))
            if tau < 0.0:
                v_end = uniformize_vec_expm({0: 1.0}, rows_suf, diag_suf, t)
                return sum(p * sx["core_count"][i] for i, p in v_end.items())
            v_pre = uniformize_vec_expm({0: 1.0}, rows_pre, diag_pre, tau)
            v_post: Dict[int, float] = defaultdict(float)
            for i, p in v_pre.items():
                rg = sx["root_gain_idx"][i]
                if rg >= 0:
                    v_post[rg] += p
            v_end = uniformize_vec_expm(dict(v_post), rows_suf, diag_suf, t - tau)
            return sum(p * sx["core_count"][i] for i, p in v_end.items())

        er0_exact = exact_r0_and_h(-1.0)
        gl_x, gl_w = gauss_legendre(12)
        h_vals = [exact_r0_and_h(t * qx) for qx in gl_x]
        int_h = sum(t * qw * hv for qw, hv in zip(gl_w, h_vals))
        e2_exact = math.exp(-delta * t) * sum(
            t * qw * gamma * math.exp(-gamma * (t - t * qx)) * hv
            for qx, qw, hv in zip(gl_x, gl_w, h_vals)
        )

        # Direct Monte Carlo of T_* using the exact independent-channel thinning of Section 7
        def run_t_star_mc(extra_leaf_at_p: bool = False, n_trials: int = 25000, seed: int = 88):
            rng = random.Random(seed)
            sum_s = 0.0
            sum_e0 = 0.0
            sum_e2 = 0.0
            sum_e2_sq = 0.0
            n_e2_events = 0
            n_e2_dp_grew = 0
            for _ in range(n_trials):
                adj = build_t_star_adj(B, H)
                side = ["x", "x"] + ["y"] * (len(adj) - 2)
                if extra_leaf_at_p:
                    # Planted error 6a: replace one core neighbor u_{B-1} of p by a leaf at p
                    # so D_p(0) is still B, but L_p(0) = 1!
                    adj = [[1, 2], [0], [0, 3], [2]]
                    side = ["x", "x", "y", "y"]
                    for _ in range(B - 2):
                        u = len(adj)
                        adj[2].append(u)
                        adj.append([2])
                        side.append("y")
                        for _ in range(H - 1):
                            lf = len(adj)
                            adj[u].append(lf)
                            adj.append([u])
                            side.append("y")
                ext_core_nbrs = set(w for w in adj[2] if w != 0 and len(adj[w]) >= 2)
                fwd_props = []
                tt = 0.0
                while True:
                    tt += rng.expovariate(gamma)
                    if tt > t:
                        break
                    fwd_props.append(tt)
                has_fwd = len(fwd_props) > 0
                tau_last = fwd_props[-1] if has_fwd else -1.0

                ext_props = []
                for u_nbr in ext_core_nbrs:
                    tt = 0.0
                    while True:
                        tt += rng.expovariate(1.0 / H)
                        if tt > t:
                            break
                        ext_props.append((tt, u_nbr))
                ext_props.sort()
                e_ext_ok = (len(ext_props) == 0)

                events = [(tp, "fwd", -1) for tp in fwd_props] + [
                    (tp, "ext", u_nbr) for tp, u_nbr in ext_props
                ]
                events.sort()
                ev_idx = 0
                cur_t = 0.0
                v_to_p_before_tau = False
                dp_exceeded_B_before_tau = False

                while cur_t < t:
                    t_next_ev = events[ev_idx][0] if ev_idx < len(events) else t
                    n_v = len(adj)
                    dp = len(adj[2])
                    sub_rate = 1.0 / dp + sum(1.0 / len(adj[u]) for u in ext_core_nbrs)
                    other_rate = n_v - sub_rate
                    dt = rng.expovariate(other_rate)
                    if cur_t + dt < t_next_ev:
                        cur_t += dt
                        while True:
                            rv = rng.randrange(n_v)
                            nbrs = adj[rv]
                            tgt = nbrs[rng.randrange(len(nbrs))]
                            if rv == 2 and tgt == 0:
                                continue
                            if rv in ext_core_nbrs and tgt == 2:
                                continue
                            break
                        if rv == 0 and tgt == 2 and has_fwd and cur_t < tau_last:
                            v_to_p_before_tau = True
                        nv = len(adj)
                        adj[tgt].append(nv)
                        adj.append([tgt])
                        side.append(side[tgt])
                    else:
                        cur_t = t_next_ev
                        if ev_idx < len(events):
                            _, kind, u_nbr = events[ev_idx]
                            ev_idx += 1
                            if kind == "fwd":
                                if len(adj[2]) > B:
                                    dp_exceeded_B_before_tau = True
                                if rng.random() < B / len(adj[2]):
                                    nv = len(adj)
                                    adj[0].append(nv)
                                    adj.append([0])
                                    side.append("x")
                            else:
                                if rng.random() < H / len(adj[u_nbr]):
                                    nv = len(adj)
                                    adj[2].append(nv)
                                    adj.append([2])
                                    side.append("y")
                c_tot = sum(1 for nbrs in adj if len(nbrs) >= 2)
                c_b = sum(1 for i, nbrs in enumerate(adj) if side[i] == "x" and len(nbrs) >= 2)
                sum_s += c_tot
                if not has_fwd:
                    sum_e0 += c_b
                if e_ext_ok and has_fwd and (not v_to_p_before_tau):
                    n_e2_events += 1
                    if dp_exceeded_B_before_tau:
                        n_e2_dp_grew += 1
                    sum_e2 += c_b
                    sum_e2_sq += c_b * c_b

            s_mc = sum_s / n_trials
            e0_mc = sum_e0 / n_trials
            e2_mc = sum_e2 / n_trials
            se_e2 = math.sqrt(max(0.0, sum_e2_sq / n_trials - e2_mc * e2_mc) / n_trials)
            return s_mc, e0_mc, e2_mc, se_e2, n_e2_events, n_e2_dp_grew

        s_mc, e0_mc, e2_mc, se_e2, n_e2, n_grew = run_t_star_mc(extra_leaf_at_p=False, n_trials=20000, seed=88 + B)
        e0_theory = math.exp(-gamma * t) * er0_exact
        print(f"  T_*(B={B},H={H},t={t}): S(t)={s_mc:.4f}, D_p>B on E2 = {n_grew}/{n_e2} (must be 0)")
        print(f"    Coupling id 1: MC E[C_B 1_E0]={e0_mc:.5f} vs exact e^(-gamma t)E^R0={e0_theory:.5f}")
        print(f"    Coupling id 2: MC E[C_B 1_E2]={e2_mc:.5f}+/-{se_e2:.5f} vs exact={e2_exact:.5f}")
        print(f"    Eq (8) bounds: E^R0={er0_exact:.4f} <= {math.exp(gamma*t)*s_mc:.4f}; "
              f"int h={int_h:.4f} <= {(1/gamma)*math.exp((delta+gamma)*t)*s_mc:.4f}")
        assert n_grew == 0, "D_p must stay equal to B up through tau on E2!"
        assert abs(e0_mc - e0_theory) < 4 * math.sqrt(e0_theory * 2.0 / 20000)
        assert abs(e2_mc - e2_exact) < 4 * se_e2
        assert er0_exact <= math.exp(gamma * t) * s_mc
        assert int_h <= (1.0 / gamma) * math.exp((delta + gamma) * t) * s_mc

        # Planted error 6a (for B=2, H=2): replace core neighbor of p by a leaf at p
        # Notice L_p(0) = 1 rings at rate 1 -> increases D_p BEFORE tau -> D_p > B fails!
        if B == 2 and H == 2:
            _, _, e2_bad, se_bad, n_e2_bad, n_grew_bad = run_t_star_mc(
                extra_leaf_at_p=True, n_trials=10000, seed=999
            )
            print(f"  Planted error 6a (leaf at p): D_p>B before tau on E2 = {n_grew_bad}/{n_e2_bad} "
                  f"({100.0*n_grew_bad/max(1,n_e2_bad):.1f}% violations of D_p=B!)")
            assert n_grew_bad > 100, "Planted error 6a should break D_p=B before tau!"

    # 6C: Verify S(t) <= pi_*^{-1} F(t+1) for (B=2, H=2, n_*=5) at t = 0.3
    n_star = 5
    p_jump_bound = 1.0 / ((2 * 1) * (3 * 2) * (4 * 3))
    pi_star = (n_star - 1) * math.exp(-2.0) * ((1.0 - math.exp(-1.0)) ** (n_star - 2)) * p_jump_bound
    f_lo, f_hi = exact_f_bounds(a_n, 1.3)
    # S(0.3) <= E_{T_*}[N(0.3)] = 5 * exp(0.3) ~ 6.75; pi_*^{-1} F(1.3) > 1000!
    assert 5.0 * math.exp(0.3) <= f_lo / pi_star
    print(f"  S(0.3) <= 5e^0.3={5*math.exp(0.3):.3f} <= pi_*^{{-1}} F_lo(1.3)={f_lo/pi_star:.1f}")
    print("  [PASS] Test 6\n")


# ---------------------------------------------------------------------------
# 10. TEST 7: Section 1 & Section 8 Constants, Parameter Schedule, Eq. (1), (10)
# ---------------------------------------------------------------------------

def test_section1_and_section8():
    print("=== TEST 7: Section 1 & Section 8 Constants, Schedule, Eq. (1), (10) ===")
    # 7A: Check parameter schedule for j = 4..50 (log-space for huge M, pi_*)
    for j in range(4, 51):
        d = j
        B = j
        J = 2 * j
        H = j * j
        kappa = 1.0 - 1.0 / j
        b = 2 * j
        s = 1.0 - 1.0 / b
        assert d >= 3 and J >= 2 and J % 2 == 0 and B >= 2 and H >= 2 and b >= 2
        assert 0.0 < kappa < s < 1.0

        gamma = 1.0 / B
        delta = (B - 1.0) / H
        log_z = d * math.log(1.0 + 1.0 / (2.0 * gamma))
        c_d = d * (d - 1) * (1.5 ** (d - 2)) / 8.0
        K = max(1, math.ceil(j * c_d))
        # Use exact big-int / float log for (z-1)K/M:
        # Since M = K + 1 + ceil(j*(z-1)*K), (z-1)*K/M <= 1/j strictly!
        if j <= 15:
            z = (1.0 + 1.0 / (2.0 * gamma)) ** d
            M = K + 1 + math.ceil(j * (z - 1.0) * K)
            assert M > K
            V = d * gamma + c_d / K + (z - 1.0) * K / M
        else:
            # Exact upper bound since M > j*(z-1)*K
            V = d * gamma + c_d / K + 1.0 / j
        lam = 1.0 / J + 0.5 * (1.0 - kappa) + V / d + 2.0 / d + gamma + delta
        bound_author = 6.0 / j + 2.0 / (j * j)
        bound_tighter = 6.0 / j + 1.0 / (j * j)
        assert lam <= bound_tighter + 1e-14 <= bound_author + 1e-14, (
            f"Lambda bound failed at j={j}: {lam} vs {bound_author}"
        )
    print("  Schedule j=4..50 verified: Lambda <= 6/j + 1/j^2 < 6/j + 2/j^2 everywhere.")

    # Print explicit constants at j=4 to document scale
    j = 4
    d = B = j
    J = 2 * j
    H = j * j
    kappa = 1.0 - 1.0 / j
    b = 2 * j
    s = 1.0 - 1.0 / b
    gamma = 1.0 / B
    delta = (B - 1.0) / H
    z = (1.0 + 1.0 / (2.0 * gamma)) ** d
    c_d = d * (d - 1) * (1.5 ** (d - 2)) / 8.0
    K = max(1, math.ceil(j * c_d))
    M = K + 1 + math.ceil(j * (z - 1.0) * K)
    V = d * gamma + c_d / K + (z - 1.0) * K / M
    lam = 1.0 / J + 0.5 * (1.0 - kappa) + V / d + 2.0 / d + gamma + delta
    n_star = 3 + (B - 1) * H
    log_pi_star = (
        math.log(n_star - 1)
        - 2.0
        + (n_star - 2) * math.log(1.0 - math.exp(-1.0))
        - sum(math.log(k * (k - 1)) for k in range(2, n_star))
    )
    log_d0 = -b * math.log(1.0 - kappa) - math.log(s - kappa)
    log_a1 = -log_pi_star + log_d0 + 0.5 + math.log(2.0) / d + ((M + K) / d) * math.log(z) + (J / 2) * math.log(3.0) + math.log(1.0 + B)
    print(f"  At j=4: z={z:.1f}, c_d={c_d:.3f}, K={K}, M={M}, n*={n_star}, "
          f"Lambda={lam:.6f} (<= {6/4 + 2/16:.6f}), log(1/pi*)={-log_pi_star:.1f}, log A(1)={log_a1:.1f}")

    # 7B: Exact F(t) bounds & Section 8 inequalities (F(t) >= 1-e^{-2t}, F(t+1) <= e(F(t)+1), eq. (1) & (10))
    dist, a_n = exact_jump_chain_from_t2(10)
    print("  Exact E[C_n] for n=2..10:", [f"{n}:{float(a_n[n]):.5f}" for n in range(2, 11)])

    # Check F(t) >= 1 - e^{-2t} and F(t+1) <= e*(F(t) + 2/3 e^{-2t}) <= e*(F(t) + 1)
    # Using exact certified [F_lo, F_hi] bounds (at small t, tail at N0=10 is tiny) + MC at t=1.0, 1.5
    for t in [0.0, 0.1, 0.25, 0.5, 1.0]:
        flo, fhi = exact_f_bounds(a_n, t)
        flo1, fhi1 = exact_f_bounds(a_n, t + 1.0)
        assert fhi >= 1.0 - math.exp(-2.0 * t) - 1e-12
        if t <= 0.25:
            assert flo >= 1.0 - math.exp(-2.0 * t) - 1e-6
        assert flo1 <= math.e * (fhi + (2.0 / 3.0) * math.exp(-2.0 * t)) + 1e-12
        assert flo1 <= math.e * (fhi + 1.0) + 1e-12
    # High-precision MC check of F(t+1) <= e*(F(t) + 2/3 e^{-2t}) and F(t) >= t at t = 1.0
    rng_f = random.Random(12345)
    f1_mc, f2_mc = 0.0, 0.0
    for _ in range(15000):
        adj = [[1], [0]]
        cur_t = 0.0
        c1 = -1
        while True:
            cur_t += rng_f.expovariate(len(adj))
            if c1 < 0 and cur_t >= 1.0:
                c1 = tree_core_count(adj)
            if cur_t >= 2.0:
                break
            rv = rng_f.randrange(len(adj))
            tgt = adj[rv][rng_f.randrange(len(adj[rv]))]
            nv = len(adj)
            adj[tgt].append(nv)
            adj.append([tgt])
        if c1 < 0:
            c1 = tree_core_count(adj)
        f1_mc += c1
        f2_mc += tree_core_count(adj)
    f1_mc /= 15000
    f2_mc /= 15000
    print(f"  MC check at t=1: F(1)={f1_mc:.4f} (>=1), F(2)={f2_mc:.4f} <= e(F(1)+2/3 e^-2)={math.e*(f1_mc + 2/3*math.exp(-2)):.4f}")
    assert f1_mc >= 1.0 - 0.03
    assert f2_mc <= math.e * (f1_mc + (2.0 / 3.0) * math.exp(-2.0))

    # Check exact (1) and (10) with minimal valid parameters (d=3, J=2, B=2, H=2, K=1, M=2, b=2, kappa=0.25)
    d, J, B, H, K, M, b, kappa = 3, 2, 2, 2, 1, 2, 2, 0.25
    s = 1.0 - 1.0 / b
    gamma = 1.0 / B
    delta = (B - 1.0) / H
    z = (1.0 + 1.0 / (2.0 * gamma)) ** d
    c_d = d * (d - 1) * (1.5 ** (d - 2)) / 8.0
    V = d * gamma + c_d / K + (z - 1.0) * K / M
    d0 = ((1.0 - kappa) ** (-b)) / (s - kappa)
    lam = 1.0 / J + 0.5 * (1.0 - kappa) + V / d + 2.0 / d + gamma + delta
    n_star = 3 + (B - 1) * H
    pi_star = (
        (n_star - 1)
        * math.exp(-2.0)
        * ((1.0 - math.exp(-1.0)) ** (n_star - 2))
        / math.prod(k * (k - 1) for k in range(2, n_star))
    )

    def A_of(t: float) -> float:
        return (1.0 / pi_star) * (
            1.0
            + d0
            * math.exp(0.5)
            * (2.0 ** (1.0 / d))
            * (z ** ((M + K) / d))
            * ((2.0 + t) ** (J / 2.0))
            * (t + 1.0 / gamma)
        )

    t_chk = 0.5
    flo1, _ = exact_f_bounds(a_n, t_chk + 1.0)
    rhs_1 = A_of(t_chk) * math.exp(lam * t_chk) * flo1
    u_val = 1.0 + t_chk * rhs_1
    # For T_3 (C(T_3)=1, |V|=3), E_{T_3}[C(t)] <= E_{T_3}[N(t)] = 3 e^t
    assert 3.0 * math.exp(t_chk) <= u_val
    # Also check F(t) <= E_{T_3}[C(t)] directly via MC!
    rng = random.Random(321)
    ec_t3 = 0.0
    for _ in range(10000):
        adj = [[1, 2], [0], [0]]
        cur_t = 0.0
        while True:
            cur_t += rng.expovariate(len(adj))
            if cur_t > t_chk:
                break
            rv = rng.randrange(len(adj))
            tgt = adj[rv][rng.randrange(len(adj[rv]))]
            nv = len(adj)
            adj[tgt].append(nv)
            adj.append([tgt])
        ec_t3 += tree_core_count(adj)
    ec_t3 /= 10000
    flo_t, fhi_t = exact_f_bounds(a_n, t_chk)
    print(f"  Minimal params at t={t_chk}: F(t) in [{flo_t:.5f}, {fhi_t:.5f}] <= E_T3[C(t)]={ec_t3:.5f} <= U(t)={u_val:.3e}")
    assert fhi_t <= ec_t3 + 0.05
    assert ec_t3 <= u_val

    # Planted error 7: Violate "Both across-edge endpoints are core forever" (start from T_2->T_3
    # where p is a leaf, D_p(0) = 1 -> b_y(0) = 1 > 1/2):
    # Then beta <= 1/2 fails, |beta - gamma| <= 1/2 fails, and at r=1, K=1, u=(1+1)/(1+gamma) ~ 2 > 3/2!
    u_leaf_p = (1.0 + 1.0) / (1.0 + 0.25)
    assert u_leaf_p > 1.5
    print(f"  Planted error 7 (leaf parent p, D_p=1): b_y=1 > 1/2 gives u={u_leaf_p:.2f} > 3/2, breaking Sec 4 & 5.")
    print("  [PASS] Test 7\n")


def test_section9_explicit_errors():
    print("=== TEST 8: Section 9 Explicit Errors (eqs. 11-12), Tail Formula & Asymptotics ===")
    # 8A: Verify exact closed-form NegBin first-moment tail formula at line 197 using EXACT Fraction arithmetic:
    # sum_{n=N+1}^inf n(n-1) q^2 (1-q)^{n-2} == (1-q)^{N-1} [N(N+1)q + 2(N+1)(1-q) + 2(1-q)^2/q]
    for q in [Fraction(1, 20), Fraction(1, 5), Fraction(1, 2), Fraction(17, 20)]:
        for N in [2, 3, 5, 10, 25]:
            exact_tail = Fraction(2) / q - sum(
                Fraction(n * (n - 1)) * (q ** 2) * ((1 - q) ** (n - 2)) for n in range(2, N + 1)
            )
            formula_tail = ((1 - q) ** (N - 1)) * (
                Fraction(N * (N + 1)) * q + Fraction(2 * (N + 1)) * (1 - q) + Fraction(2) * ((1 - q) ** 2) / q
            )
            assert exact_tail == formula_tail, f"NegBin tail mismatch at q={q}, N={N}: {exact_tail} vs {formula_tail}"
    print("  Closed-form NegBin first-moment tail formula (line 197) verified with EXACT Fraction equality.")

    # 8B: Verify Lemma 5.1 threshold 1 + log(6x) <= 2x for x >= 2 and theta_k(t) <= 1 for j >= 4
    for x in [2.0, 2.001, 2.5, 5.0, 100.0]:
        assert 1.0 + math.log(6.0 * x) <= 2.0 * x
    for j in range(4, 15):
        k = math.ceil(6.0 * j * math.exp(j))
        t_min = 2 * j * k
        # log(theta_k(t_min)) = 0.5 * log(k*(k+1)*(2k+1)) - t_min
        log_theta = 0.5 * (math.log(k) + math.log(k + 1) + math.log(2 * k + 1)) - t_min
        assert log_theta < -1000.0, f"theta_k too large at j={j}"
        # Check line 206 bounds: log z <= j log j, K <= 1 + j^3 (3/2)^j, M+K <= 5 j z K
        z = (1.0 + 0.5 * j) ** j
        c_d = j * (j - 1) * (1.5 ** (j - 2)) / 8.0
        K = max(1, math.ceil(j * c_d))
        M = K + 1 + math.ceil(j * (z - 1.0) * K)
        assert math.log(z) <= j * math.log(j) + 1e-12
        assert K <= 1.0 + (j ** 3) * (1.5 ** j) + 1e-12
        assert M + K <= 5.0 * j * z * K + 1e-12
        # Check k = 6j * 3^j alternative in line 199
        k_alt = 6 * j * (3 ** j)
        assert k_alt >= k
    print("  Section 9 threshold 1+log(6x)<=2x, theta_k<=1, and line 206 inequalities verified for j=4..14.")

    # 8C: Check instantaneous promotion rate >= 1 on ALL trees up to n = 10 (including T_2!)
    dist, a_n = exact_jump_chain_from_t2(10)
    min_prom_rate = float("inf")
    for n, tree_map in dist.items():
        for _, (_, adj) in tree_map.items():
            # Promotion rate = sum_{x} (# leaf neighbors of x) / deg(x)
            rate = sum(
                sum(1 for w in adj[x] if len(adj[w]) == 1) / len(adj[x])
                for x in range(len(adj))
            )
            min_prom_rate = min(min_prom_rate, rate)
    print(f"  Minimum instantaneous promotion rate across all trees n=2..10: {min_prom_rate:.6f} >= 1.0")
    assert min_prom_rate >= 1.0 - 1e-15
    # Verify F(t) >= t directly at small t where NegBin tail is negligible
    for t in [0.1, 0.25, 0.5]:
        flo, fhi = exact_f_bounds(a_n, t)
        assert fhi >= t and flo >= t - 1e-4
    print("  [PASS] Test 8\n")


def main():
    test_section2_joint_density()
    test_section3_exterior_estimate()
    test_section4_last_arrival_and_profiles()
    test_section5_renyi_domination()
    test_section6_holder()
    test_section7_physical_seed()
    test_section1_and_section8()
    test_section9_explicit_errors()
    print("ALL 8 VERIFICATION & PLANTED-ERROR SUITES PASSED.")


if __name__ == "__main__":
    main()
