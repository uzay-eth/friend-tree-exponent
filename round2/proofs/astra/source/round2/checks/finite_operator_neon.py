#!/usr/bin/env python3
"""Independent cold-audit checker for round2/finite/finite-operator.md.

Single-CPU, stdlib-only independent reconstruction of all six exact rational
branching matrices, Perron-Frobenius / Collatz-Wielandt certified spectral
brackets, characteristic polynomial checks, certificate cross-check against
round2/outputs/finite-matrices.json, matrix mutation tests, and direct
Monte Carlo / coupling invariant verification for Sections 3-5.
"""

from __future__ import annotations

import json
import math
import random
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# 1. Independent canonical unrooted attributed tree representation
# ---------------------------------------------------------------------------

# An active component is represented canonically by an unrooted AHU tree key
# with vertex attributes (D_u, L_u). A concrete instance is (attrs, adj) where:
#   attrs: List[Tuple[int, int]]  -- (D_u, L_u) for u = 0..k-1
#   adj:   List[List[int]]        -- adjacency list of active edges
HUB_KEY = ("HUB",)


def canonical_active_tree(attrs: Sequence[Tuple[int, int]], adj: Sequence[Sequence[int]]) -> Tuple:
    """Center-peeling AHU canonical form for an unrooted vertex-attributed tree."""
    n = len(attrs)
    if n == 1:
        return (attrs[0], ())
    deg = [len(adj[i]) for i in range(n)]
    rem = n
    alive = [True] * n
    leaves = [i for i in range(n) if deg[i] <= 1]
    while rem > 2:
        nxt = []
        for u in leaves:
            alive[u] = False
            rem -= 1
            for w in adj[u]:
                if alive[w]:
                    deg[w] -= 1
                    if deg[w] == 1:
                        nxt.append(w)
        leaves = nxt
    centers = [i for i in range(n) if alive[i]]

    def encode_rooted(u: int, parent: int) -> Tuple:
        kids = sorted(encode_rooted(w, u) for w in adj[u] if w != parent)
        return (attrs[u], tuple(kids))

    return min(encode_rooted(c, -1) for c in centers)


def connected_components_without(
    attrs: Sequence[Tuple[int, int]],
    adj: Sequence[Sequence[int]],
    drop_vertex: int = -1,
    drop_edge: Optional[Tuple[int, int]] = None,
) -> List[Tuple[List[Tuple[int, int]], List[List[int]]]]:
    """Return connected components after deleting a vertex or an edge."""
    n = len(attrs)
    visited = [False] * n
    if drop_vertex >= 0:
        visited[drop_vertex] = True
    bad_e = None
    if drop_edge is not None:
        a, b = drop_edge
        bad_e = {(a, b), (b, a)}

    comps = []
    for start in range(n):
        if visited[start]:
            continue
        nodes = []
        q = deque([start])
        visited[start] = True
        while q:
            u = q.popleft()
            nodes.append(u)
            for w in adj[u]:
                if visited[w]:
                    continue
                if bad_e is not None and (u, w) in bad_e:
                    continue
                visited[w] = True
                q.append(w)
        remap = {u: idx for idx, u in enumerate(nodes)}
        c_attrs = [attrs[u] for u in nodes]
        c_adj: List[List[int]] = [[] for _ in nodes]
        for u in nodes:
            iu = remap[u]
            for w in adj[u]:
                if w in remap and (bad_e is None or (u, w) not in bad_e):
                    c_adj[iu].append(remap[w])
        comps.append((c_attrs, c_adj))
    return comps


# ---------------------------------------------------------------------------
# 2. Independent exact rational generator builder (with optional mutations!)
# ---------------------------------------------------------------------------

def build_branching_operator(
    M: int,
    B: int,
    *,
    mutate_cut_unique_outcomes: bool = False,
    mutate_drop_slots: bool = False,
    mutate_slot_send_rate: Fraction = Fraction(0),
    mutate_hub_newborn_leaves: int = 1,
):
    """Independently build exact rational generator G for (M, B) from Section 1-2.

    Optional mutation flags implement planted errors to test sensitivity.
    """
    assert M >= 3 and B >= 1
    s = (B + (M - 1) - 1) // (M - 1)  # ceil(B / (M - 1))

    # Founder: singleton (D=2, L=2) with no active edges
    founder_attrs = [(2, 2)]
    founder_adj: List[List[int]] = [[]]
    founder_key = canonical_active_tree(founder_attrs, founder_adj)

    key_to_idx: Dict[Tuple, int] = {founder_key: 0}
    idx_to_key: List[Tuple] = [founder_key]
    idx_to_rep: List[Optional[Tuple[List[Tuple[int, int]], List[List[int]]]]] = [
        (founder_attrs, founder_adj)
    ]
    core_counts: List[int] = [1]

    def register(key: Tuple, rep) -> int:
        if key not in key_to_idx:
            idx = len(idx_to_key)
            key_to_idx[key] = idx
            idx_to_key.append(key)
            idx_to_rep.append(rep)
            core_counts.append(1 if key == HUB_KEY else len(rep[0]))
            queue.append(idx)
        return key_to_idx[key]

    rows: List[Dict[int, Fraction]] = []
    queue = deque([0])
    processed = 0

    while processed < len(idx_to_key):
        idx = processed
        processed += 1
        key = idx_to_key[idx]
        row: Dict[int, Fraction] = defaultdict(Fraction)

        if key == HUB_KEY:
            # Hub stays alive and creates one singleton (2, mutate_hub_newborn_leaves) at rate 1
            nb_attrs = [(2, mutate_hub_newborn_leaves)]
            nb_adj: List[List[int]] = [[]]
            nb_idx = register(canonical_active_tree(nb_attrs, nb_adj), (nb_attrs, nb_adj))
            # nu_H = 1 (cancels -1_{H=H}), nu_nb = 1
            row[nb_idx] += Fraction(1)
            rows.append(dict(row))
            continue

        rep = idx_to_rep[idx]
        assert rep is not None
        attrs, adj = rep
        n = len(attrs)

        # Verify state invariants on every reachable active type
        if not mutate_drop_slots:
            assert 1 <= n <= B
            for u in range(n):
                du, lu = attrs[u]
                au = len(adj[u])
                assert 2 <= du < M, f"Invalid degree {du}"
                assert lu >= 0, f"Negative leaves {lu}"
                assert du - lu - au >= 0, f"Negative slots: D={du}, L={lu}, a={au}"

        # Transition 1: Gain at each vertex u
        for u in range(n):
            du, lu = attrs[u]
            au = len(adj[u])
            slots = du - lu - au
            rate = Fraction(lu)
            for w in adj[u]:
                rate += Fraction(1, attrs[w][0])
            if mutate_slot_send_rate != 0:
                rate += slots * mutate_slot_send_rate
            if rate == 0:
                continue

            new_du = du + 1
            new_lu = lu + 1
            row[idx] -= rate
            if new_du >= M:
                # Replace u by hub and release connected components of A \ {u}
                h_idx = register(HUB_KEY, None)
                row[h_idx] += rate
                if mutate_drop_slots:
                    # Planted error: decrement D_w of former neighbors of u
                    mod_attrs = list(attrs)
                    for w in adj[u]:
                        dw, lw = mod_attrs[w]
                        mod_attrs[w] = (max(2, dw - 1), lw)
                    comps = connected_components_without(mod_attrs, adj, drop_vertex=u)
                else:
                    comps = connected_components_without(attrs, adj, drop_vertex=u)
                for ca, cadj in comps:
                    ck = canonical_active_tree(ca, cadj)
                    c_idx = register(ck, (ca, cadj))
                    row[c_idx] += rate
            else:
                new_attrs = list(attrs)
                new_attrs[u] = (new_du, new_lu)
                nk = canonical_active_tree(new_attrs, adj)
                n_idx = register(nk, (new_attrs, [list(nbrs) for nbrs in adj]))
                row[n_idx] += rate

        # Transition 2: Promotion at each vertex u with L_u >= 1
        for u in range(n):
            du, lu = attrs[u]
            if lu <= 0:
                continue
            rate = Fraction(lu, du)
            row[idx] -= rate

            # Build post-promotion tree of size n + 1
            new_attrs = list(attrs)
            new_attrs[u] = (du, lu - 1)
            new_attrs.append((2, 1))
            new_adj = [list(nbrs) for nbrs in adj]
            new_adj[u].append(n)
            new_adj.append([u])

            if n + 1 <= B:
                nk = canonical_active_tree(new_attrs, new_adj)
                n_idx = register(nk, (new_attrs, new_adj))
                row[n_idx] += rate
            else:
                # Size B + 1: enumerate ALL legal physical edges (u_e < v_e) whose two sides have >= s cores
                legal_outcomes = []
                for a in range(n + 1):
                    for b in new_adj[a]:
                        if a < b:
                            if mutate_drop_slots:
                                cut_attrs = list(new_attrs)
                                cut_attrs[a] = (max(2, cut_attrs[a][0] - 1), cut_attrs[a][1])
                                cut_attrs[b] = (max(2, cut_attrs[b][0] - 1), cut_attrs[b][1])
                                comps = connected_components_without(cut_attrs, new_adj, drop_edge=(a, b))
                            else:
                                comps = connected_components_without(new_attrs, new_adj, drop_edge=(a, b))
                            assert len(comps) == 2
                            if len(comps[0][0]) >= s and len(comps[1][0]) >= s:
                                k0 = canonical_active_tree(comps[0][0], comps[0][1])
                                k1 = canonical_active_tree(comps[1][0], comps[1][1])
                                legal_outcomes.append((k0, comps[0], k1, comps[1]))
                assert len(legal_outcomes) >= 1, "Centroid theorem violated: no legal cut!"

                if mutate_cut_unique_outcomes:
                    # Planted error: deduplicate isomorphic unordered outcomes before averaging
                    seen = {}
                    for k0, c0, k1, c1 in legal_outcomes:
                        pair_key = (k0, k1) if k0 <= k1 else (k1, k0)
                        if pair_key not in seen:
                            seen[pair_key] = (k0, c0, k1, c1)
                    chosen_list = list(seen.values())
                else:
                    chosen_list = legal_outcomes

                w_edge = rate / Fraction(len(chosen_list))
                for k0, c0, k1, c1 in chosen_list:
                    i0 = register(k0, c0)
                    i1 = register(k1, c1)
                    row[i0] += w_edge
                    row[i1] += w_edge

        # Prune zero entries
        rows.append({j: v for j, v in row.items() if v != 0})

    # Verify off-diagonal nonnegativity and core-count birth rate <= C_tilde
    if not mutate_drop_slots and mutate_slot_send_rate == 0:
        for i, row in enumerate(rows):
            for j, val in row.items():
                if j != i:
                    assert val >= 0, f"Negative off-diagonal G[{i},{j}]={val}"
            gc_i = sum(val * core_counts[j] for j, val in row.items())
            assert gc_i <= Fraction(core_counts[i]), (
                f"Core birth rate exceeds C_tilde at type {i}: (Gc)_i={gc_i} > c_i={core_counts[i]}"
            )

    return {
        "M": M,
        "B": B,
        "s": s,
        "rows": rows,
        "keys": idx_to_key,
        "reps": idx_to_rep,
        "core_counts": core_counts,
    }


# ---------------------------------------------------------------------------
# 3. Tarjan SCCs, exact Collatz-Wielandt Perron certification, exact det()
# ---------------------------------------------------------------------------

def tarjan_sccs(n: int, rows: Sequence[Dict[int, Fraction]]) -> List[List[int]]:
    """Strongly connected components of off-diagonal positive graph of G."""
    index = 0
    stack: List[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    lowlink = [0] * n
    sccs: List[List[int]] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w, val in rows[v].items():
            if w == v or val <= 0:
                continue
            if indices[w] == -1:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif on_stack[w]:
                lowlink[v] = min(lowlink[v], indices[w])
        if lowlink[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            sccs.append(sorted(comp))

    for v in range(n):
        if indices[v] == -1:
            strongconnect(v)
    return sccs


def certify_scc_spectral_abscissa(
    rows: Sequence[Dict[int, Fraction]],
    scc: Sequence[int],
) -> Tuple[Fraction, Fraction, List[int]]:
    """Compute certified rational interval [lo, hi] for Perron eigenvalue of G_{SS}
    using exact Collatz-Wielandt quotients min_i (G_SS w)_i / w_i <= lambda <= max_i (G_SS w)_i / w_i.
    """
    k = len(scc)
    if k == 1:
        u = scc[0]
        val = rows[u].get(u, Fraction(0))
        return val, val, [1]

    pos = {u: idx for idx, u in enumerate(scc)}
    sub_exact: List[List[Tuple[int, Fraction]]] = [[] for _ in range(k)]
    sub_float: List[List[Tuple[int, float]]] = [[] for _ in range(k)]
    min_diag = 0.0
    for i, u in enumerate(scc):
        for v, val in rows[u].items():
            if v in pos:
                j = pos[v]
                sub_exact[i].append((j, val))
                fv = float(val)
                sub_float[i].append((j, fv))
                if i == j and fv < min_diag:
                    min_diag = fv

    # Shifted nonnegative matrix A = G_SS + shift * I (use shift > -min_diag + 1)
    shift = math.ceil(-min_diag) + 1.0
    # Power iteration in float, then refine with Rayleigh-quotient / shifted inverse iteration?
    # Even simpler: power iteration on A = G_SS + shift*I, or if slow, do exact Fraction
    # shifted inverse iteration? Wait: with float power iteration + inverse iteration (Gaussian elim)
    # on a <= 55x55 matrix we get 15-16 digits in milliseconds!
    w = [1.0] * k
    for _ in range(300):
        nw = [shift * w[i] for i in range(k)]
        for i in range(k):
            for j, fv in sub_float[i]:
                nw[i] += fv * w[j]
        norm = max(nw)
        w = [x / norm for x in nw]

    # Estimate lambda via Rayleigh-like quotient
    aw = [0.0] * k
    for i in range(k):
        s_val = 0.0
        for j, fv in sub_float[i]:
            s_val += fv * w[j]
        aw[i] = s_val
    lam_est = sum(aw[i] * w[i] for i in range(k)) / sum(w[i] * w[i] for i in range(k))

    # 3 steps of exact-ish shifted inverse iteration in high-precision float (with tiny nudge)
    # to nail eigenvector to full float64 precision (1e-15), then scale to 10^15 integers!
    def solve_linear_float(mat: List[List[float]], rhs: List[float]) -> List[float]:
        n = len(rhs)
        a = [row[:] + [rhs[i]] for i, row in enumerate(mat)]
        for col in range(n):
            piv = max(range(col, n), key=lambda r: abs(a[r][col]))
            a[col], a[piv] = a[piv], a[col]
            pv = a[col][col]
            if abs(pv) < 1e-18:
                return rhs
            inv_pv = 1.0 / pv
            for j in range(col, n + 1):
                a[col][j] *= inv_pv
            for r in range(n):
                if r != col and a[r][col] != 0.0:
                    fac = a[r][col]
                    for j in range(col, n + 1):
                        a[r][j] -= fac * a[col][j]
        return [a[i][n] for i in range(n)]

    for step in range(4):
        mu = lam_est + (1e-12 if step == 0 else 1e-15)
        mat = [[0.0] * k for _ in range(k)]
        for i in range(k):
            mat[i][i] = -mu
            for j, fv in sub_float[i]:
                mat[i][j] += fv
        y = solve_linear_float(mat, w)
        if all(x < 0 for x in y):
            y = [-x for x in y]
        if all(x > 0 for x in y):
            norm = max(y)
            w = [x / norm for x in y]
            aw = [sum(fv * w[j] for j, fv in sub_float[i]) for i in range(k)]
            lam_est = sum(aw[i] * w[i] for i in range(k)) / sum(w[i] * w[i] for i in range(k))

    # Refine from float64 eigenvector using exact integer power iteration on shifted L*G_SS + shift*I
    # until Collatz-Wielandt gap < 1e-17!
    den_lcm = 1
    for i in range(k):
        for _, val in sub_exact[i]:
            den_lcm = math.lcm(den_lcm, val.denominator)
    int_rows = [[(j, int(val * den_lcm)) for j, val in sub_exact[i]] for i in range(k)]
    min_int_diag = min((v for i in range(k) for j, v in int_rows[i] if i == j), default=0)
    int_shift = max(0, -min_int_diag) + 1

    scale = 10**25
    w_int = [max(1, int(round(x * scale))) for x in w]
    for it in range(5000):
        # Compute y_i = (L*G_SS w)_i and shifted_y_i = y_i + int_shift * w_i
        y_vals = [sum(v * w_int[j] for j, v in int_rows[i]) for i in range(k)]
        if it % 8 == 0:
            ratios = [Fraction(y_vals[i], den_lcm * w_int[i]) for i in range(k)]
            lo_r, hi_r = min(ratios), max(ratios)
            if hi_r - lo_r < Fraction(1, 10**17):
                return lo_r, hi_r, w_int
        nxt = [y_vals[i] + int_shift * w_int[i] for i in range(k)]
        mx = max(nxt)
        w_int = [max(1, (val * scale + mx // 2) // mx) for val in nxt]
    ratios = [Fraction(sum(v * w_int[j] for j, v in int_rows[i]), den_lcm * w_int[i]) for i in range(k)]
    return min(ratios), max(ratios), w_int


def exact_det_int(mat: List[List[int]]) -> int:
    """Exact integer determinant via fraction-free Bareiss algorithm."""
    n = len(mat)
    a = [row[:] for row in mat]
    sign = 1
    prev_pivot = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap_r = -1
            for r in range(k + 1, n):
                if a[r][k] != 0:
                    swap_r = r
                    break
            if swap_r == -1:
                return 0
            a[k], a[swap_r] = a[swap_r], a[k]
            sign = -sign
        piv = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * piv - a[i][k] * a[k][j]) // prev_pivot
        prev_pivot = piv
    return sign * a[n - 1][n - 1]


def verify_charpoly(
    rows: Sequence[Dict[int, Fraction]],
    scc: Sequence[int],
    scale_L: int,
    coeffs_desc: Sequence[int],
) -> None:
    """Verify det(x I - L G_{SS}) == sum coeffs_desc[r] x^{k-r} at k+1 integer points."""
    k = len(scc)
    assert len(coeffs_desc) == k + 1
    pos = {u: idx for idx, u in enumerate(scc)}
    int_mat = [[0] * k for _ in range(k)]
    for i, u in enumerate(scc):
        for v, val in rows[u].items():
            if v in pos:
                scaled = val * scale_L
                assert scaled.denominator == 1, f"Scale L={scale_L} does not clear {val}"
                int_mat[i][pos[v]] = scaled.numerator

    # Evaluate at k + 1 distinct integers x = -k//2 .. -k//2 + k
    x0 = -(k // 2)
    for x in range(x0, x0 + k + 1):
        m_x = [[-int_mat[i][j] + (x if i == j else 0) for j in range(k)] for i in range(k)]
        det_val = exact_det_int(m_x)
        poly_val = 0
        for c in coeffs_desc:
            poly_val = poly_val * x + c
        assert det_val == poly_val, f"Charpoly mismatch at x={x}: det={det_val} vs poly={poly_val}"


# ---------------------------------------------------------------------------
# 4. TEST A: Recompute all 6 exact matrices + compare to draft & certificate
# ---------------------------------------------------------------------------

TABLE_EXPECTED = [
    # (M, B, types, classes, max_cls_size, lo_str, hi_str)
    (3, 1, 4, 3, 2, "0.61803398874989", "0.61803398874990"),
    (3, 3, 9, 5, 5, "0.67833737318518", "0.67833737318519"),
    (3, 9, 27, 11, 17, "0.67928074796982", "0.67928074796983"),
    (4, 1, 8, 5, 4, "0.55612325604897", "0.55612325604898"),
    (4, 2, 23, 8, 16, "0.61485063481498", "0.61485063481499"),
    (4, 3, 68, 14, 55, "0.62755802142835", "0.62755802142836"),
]


def test_all_six_matrices_and_cert():
    print("=== TEST A: Independent Recomputation of All 6 Exact Matrices ===")
    recomputed = {}

    for M, B, exp_types, exp_classes, exp_max_sz, lo_s, hi_s in TABLE_EXPECTED:
        op = build_branching_operator(M, B)
        n = len(op["rows"])
        sccs = tarjan_sccs(n, op["rows"])
        assert n == exp_types, f"({M},{B}) types={n} != {exp_types}"
        assert len(sccs) == exp_classes, f"({M},{B}) classes={len(sccs)} != {exp_classes}"

        # Certify every SCC
        scc_bounds = []
        for scc in sccs:
            lo, hi, _ = certify_scc_spectral_abscissa(op["rows"], scc)
            scc_bounds.append((lo, hi, scc))
        scc_bounds.sort(key=lambda item: item[0], reverse=True)
        win_lo, win_hi, win_scc = scc_bounds[0]
        assert len(win_scc) == exp_max_sz, f"({M},{B}) max class size={len(win_scc)} != {exp_max_sz}"
        # Dominant class strictly beats all other classes!
        for other_lo, other_hi, _ in scc_bounds[1:]:
            assert win_lo > other_hi, f"({M},{B}) dominant class not strictly separated!"

        exp_lo = Fraction(lo_s)
        exp_hi = Fraction(hi_s)
        # Our certified [win_lo, win_hi] MUST lie inside [exp_lo, exp_hi] (or overlap? Better:
        # since win_hi - win_lo < 1e-14, check exp_lo <= win_lo <= win_hi <= exp_hi!)
        gap = float(win_hi - win_lo)
        print(f"  (M={M}, B={B}): types={n}, classes={len(sccs)}, win_sz={len(win_scc)}, "
              f"lam in [{float(win_lo):.15f}, {float(win_hi):.15f}] (gap={gap:.2e})")
        assert exp_lo <= win_lo and win_hi <= exp_hi, (
            f"({M},{B}) certified [{float(win_lo):.15f}, {float(win_hi):.15f}] "
            f"not inside draft [{lo_s}, {hi_s}]!"
        )
        recomputed[(M, B)] = (op, sccs, win_scc, win_lo, win_hi)

    # Check explicit matrix for (3, 1) (lines 105-113)
    op31, _, win31, _, _ = recomputed[(3, 1)]
    k_founder = canonical_active_tree([(2, 2)], [[]])
    k_hub = HUB_KEY
    k_s1 = canonical_active_tree([(2, 1)], [[]])
    k_s0 = canonical_active_tree([(2, 0)], [[]])
    order_keys = [k_founder, k_hub, k_s1, k_s0]
    kmap = {k: i for i, k in enumerate(op31["keys"])}
    perm = [kmap[k] for k in order_keys]
    mat31 = [[op31["rows"][perm[i]].get(perm[j], Fraction(0)) for j in range(4)] for i in range(4)]
    expected_31 = [
        [Fraction(-3), Fraction(2), Fraction(2), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1), Fraction(-1), Fraction(1, 2)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
    ]
    assert mat31 == expected_31, f"(3,1) matrix mismatch: {mat31}"
    print("  Exact (3,1) 4x4 matrix matches lines 107-112 identically.")

    # Check explicit characteristic polynomials in line 114!
    # (3,1): x^2 + x - 1 = 0 with x = lambda_{3,1} (scale L=1)
    verify_charpoly(op31["rows"], win31, 1, [1, 1, -1])
    # (3,3): x^5 + 17 x^4 + 94 x^3 + 162 x^2 - 132 x - 416 = 0 with x = 2 * lambda_{3,3}
    op33, _, win33, _, _ = recomputed[(3, 3)]
    verify_charpoly(op33["rows"], win33, 2, [1, 17, 94, 162, -132, -416])
    # (4,1): x^4 + 15 x^3 + 62 x^2 + 12 x - 270 = 0 with x = 3 * lambda_{4,1}
    op41, _, win41, _, _ = recomputed[(4, 1)]
    verify_charpoly(op41["rows"], win41, 3, [1, 15, 62, 12, -270])
    print("  Exact charpolys for (3,1), (3,3), (4,1) in line 114 verified by exact Bareiss det().")

    # Cross-check round2/outputs/finite-matrices.json completely independently!
    cert_path = Path(__file__).resolve().parents[1] / "outputs" / "finite-matrices.json"
    with cert_path.open("r", encoding="utf-8") as f:
        cert = json.load(f)
    # Inspect structure flexibly
    entries = cert["cases"] if isinstance(cert, dict) and "cases" in cert else (
        cert["matrices"] if isinstance(cert, dict) and "matrices" in cert else cert
    )
    assert isinstance(entries, list) and len(entries) == 6, f"Unexpected cert format: {type(entries)}"

    def parse_cert_type(enc) -> Tuple:
        if enc == "hub":
            return HUB_KEY
        # Nested triple [D, L, [child_enc, ...]]
        attrs: List[Tuple[int, int]] = []
        adj: List[List[int]] = []

        def dfs(node) -> int:
            u = len(attrs)
            attrs.append((int(node[0]), int(node[1])))
            adj.append([])
            for ch in node[2]:
                w = dfs(ch)
                adj[u].append(w)
                adj[w].append(u)
            return u

        dfs(enc)
        return canonical_active_tree(attrs, adj)

    for entry in entries:
        M, B = int(entry["M"]), int(entry["B"])
        op, _, win_scc, _, _ = recomputed[(M, B)]
        cert_types = entry["types"]
        cert_rows = entry["rows"]
        assert len(cert_types) == len(op["rows"])
        our_idx = {k: i for i, k in enumerate(op["keys"])}
        perm_c_to_o = [our_idx[parse_cert_type(enc)] for enc in cert_types]
        assert len(set(perm_c_to_o)) == len(op["rows"]), "Non-bijective type map with cert!"
        for ci, row_triples in enumerate(cert_rows):
            oi = perm_c_to_o[ci]
            cert_map = {perm_c_to_o[int(cj)]: Fraction(int(num), int(den)) for cj, num, den in row_triples}
            assert cert_map == op["rows"][oi], f"Cert row mismatch at ({M},{B}) row {ci}"
        # Verify recorded dominant charpoly in cert (char_scale / char_coefficients)
        L_scale = int(entry["char_scale"])
        coeffs_desc = [int(c) for c in entry["char_coefficients"]]
        assert coeffs_desc[0] == 1
        verify_charpoly(op["rows"], win_scc, L_scale, coeffs_desc)

        # Verify opposite endpoint signs & positive derivative of charpoly on [L*lo, L*hi]
        br_den = int(entry["bracket_denominator"])
        br_lo = Fraction(int(entry["lower_numerator"]), br_den)
        br_hi = Fraction(int(entry["upper_numerator"]), br_den)
        x_lo = L_scale * br_lo
        x_hi = L_scale * br_hi

        def eval_poly(cs: Sequence[int], x: Fraction) -> Fraction:
            acc = Fraction(0)
            for c in cs:
                acc = acc * x + c
            return acc

        p_lo = eval_poly(coeffs_desc, x_lo)
        p_hi = eval_poly(coeffs_desc, x_hi)
        assert p_lo < 0 and p_hi > 0, f"Charpoly endpoint signs failed at ({M},{B}): {p_lo}, {p_hi}"
        # Derivative P'(x) on [x_lo, x_hi]: use exact rational interval arithmetic!
        deg = len(coeffs_desc) - 1
        d_coeffs = [coeffs_desc[r] * (deg - r) for r in range(deg)]
        # Since x_lo > 0, interval bound on P'(x) = sum d_k x^k:
        # for each term c*x^p: if c >= 0 min is c*x_lo^p; if c < 0 min is c*x_hi^p.
        d_min = Fraction(0)
        for r, c in enumerate(d_coeffs):
            pwr = deg - 1 - r
            if c >= 0:
                d_min += c * (x_lo ** pwr)
            else:
                d_min += c * (x_hi ** pwr)
        assert d_min > 0, f"Charpoly derivative not certified positive at ({M},{B}): d_min={d_min}"

        # Verify recorded positive integer Collatz-Wielandt witnesses in EVERY class (gap < 1e-16)
        for cls_obj in entry["classes"]:
            c_idx = [perm_c_to_o[int(i)] for i in cls_obj["indices"]]
            w_vec = [int(w) for w in cls_obj["weights"]]
            assert all(w > 0 for w in w_vec)
            pos_c = {u: i for i, u in enumerate(c_idx)}
            ratios = []
            for i, u in enumerate(c_idx):
                num = sum(val * w_vec[pos_c[v]] for v, val in op["rows"][u].items() if v in pos_c)
                ratios.append(num / Fraction(w_vec[i]))
            c_gap = float(max(ratios) - min(ratios))
            assert c_gap < 1e-16, f"Cert witness gap >= 1e-16 at ({M},{B}): {c_gap}"
    print("  All 6 matrices, charpolys, endpoint signs, P'>0, and <1e-16 witnesses in finite-matrices.json verified 100%!")
    print("  [PASS] Test A\n")
    return recomputed


# ---------------------------------------------------------------------------
# 5. TEST B: Matrix Mutations & Planted Errors
# ---------------------------------------------------------------------------

def test_matrix_mutations(recomputed):
    print("=== TEST B: Matrix Mutations & Planted-Error Sensitivity ===")
    # Mutation 1: Deduplicate isomorphic outcomes before averaging cut edges!
    # For (M=3, B=9) or (M=4, B=3), check if physical-edge averaging differs from unique-outcome averaging!
    found_cut_diff = False
    for M, B in [(3, 9), (4, 3)]:
        op_true, _, _, lo_true, hi_true = recomputed[(M, B)]
        op_mut = build_branching_operator(M, B, mutate_cut_unique_outcomes=True)
        sccs_m = tarjan_sccs(len(op_mut["rows"]), op_mut["rows"])
        best_m = max(
            (certify_scc_spectral_abscissa(op_mut["rows"], sc) for sc in sccs_m),
            key=lambda x: x[0],
        )
        if best_m[1] < lo_true or best_m[0] > hi_true:
            found_cut_diff = True
            print(f"  Mutation 1 (unique-outcome cut avg) at ({M},{B}): "
                  f"lam_mut={float(best_m[0]):.10f} != lam_true={float(lo_true):.10f}")
    assert found_cut_diff, "Expected cut-edge averaging mutation to change lambda!"

    # Mutation 2: Drop inactive slots on hub removal / cut (decrement neighbor D_w)
    op_m2 = build_branching_operator(4, 2, mutate_drop_slots=True)
    sccs_m2 = tarjan_sccs(len(op_m2["rows"]), op_m2["rows"])
    best_m2 = max(
        (certify_scc_spectral_abscissa(op_m2["rows"], sc) for sc in sccs_m2),
        key=lambda x: x[0],
    )
    lo_42 = recomputed[(4, 2)][3]
    print(f"  Mutation 2 (drop inactive slots) at (4,2): "
          f"lam_mut={float(best_m2[0]):.10f} != lam_true={float(lo_42):.10f}")
    assert abs(float(best_m2[0]) - float(lo_42)) > 1e-3

    # Mutation 3: Inactive slots send leaves at rate 1/M instead of 0
    op_m3 = build_branching_operator(3, 3, mutate_slot_send_rate=Fraction(1, 3))
    sccs_m3 = tarjan_sccs(len(op_m3["rows"]), op_m3["rows"])
    best_m3 = max(
        (certify_scc_spectral_abscissa(op_m3["rows"], sc) for sc in sccs_m3),
        key=lambda x: x[0],
    )
    lo_33 = recomputed[(3, 3)][3]
    print(f"  Mutation 3 (inactive slots send at 1/M) at (3,3): "
          f"lam_mut={float(best_m3[0]):.10f} != lam_true={float(lo_33):.10f}")
    assert abs(float(best_m3[0]) - float(lo_33)) > 1e-3

    # Mutation 4: Hub creates (2, 2) instead of (2, 1) (forgetting hub slot on newborn)
    op_m4 = build_branching_operator(3, 1, mutate_hub_newborn_leaves=2)
    sccs_m4 = tarjan_sccs(len(op_m4["rows"]), op_m4["rows"])
    best_m4 = max(
        (certify_scc_spectral_abscissa(op_m4["rows"], sc) for sc in sccs_m4),
        key=lambda x: x[0],
    )
    lo_31 = recomputed[(3, 1)][3]
    print(f"  Mutation 4 (hub newborn (2,2) not (2,1)) at (3,1): "
          f"lam_mut={float(best_m4[0]):.10f} != lam_true={float(lo_31):.10f}")
    assert abs(float(best_m4[0]) - float(lo_31)) > 0.1

    # Mutation 5: Why dropping unreachable classes is essential in eq. (13)!
    # Check if there exists a valid (2 <= D < M, L >= 0, D-L-a >= 0) component NOT reachable
    # from the founder, e.g. at M=4, B=2: two active cores both with (D=3, L=2, a=1)!
    # Wait: is (3,2)-(3,2) reachable? Let's check!
    # At M=3, B=1: (2,2), (2,1), (2,0), hub are all 4.
    # At M=3, B=2 (s=1): consider a 2-core component (2,1)-(2,1): reachable?
    # Wait: founder (2,2) promotes -> (2,1)-(2,1)! So that's reachable.
    # What about (M=4, B=1)? Valid singletons have 2<=D<=3, 0<=L<=D:
    # (2,2),(2,1),(2,0),(3,3),(3,2),(3,1),(3,0) + hub = 8 types (all reachable!).
    # What about (M=4, B=2)? How many valid trees of size <= 2 exist vs 23 reachable?
    # Size 1: 7 types. Size 2: pairs of (D, L) with 2<=D<=3, 0<=L<=D-1:
    # valid vertex states with a=1 are (2,1),(2,0),(3,2),(3,1),(3,0) -> 5 states -> C(5+1,2)=15 pairs!
    # Total valid types = 1(hub) + 7 + 15 = 23 (all reachable!).
    # What about (M=3, B=3)? Size 1: 3 types ((2,2),(2,1),(2,0)).
    # Size 2: vertex states with a=1 are (2,1),(2,0) -> 3 pairs!
    # Size 3 (max active deg <= 2 -> path v0-v1-v2):
    #   middle v1 has a=2, D=2 -> MUST be (2,0)!
    #   ends v0, v2 have a=1 -> in {(2,1),(2,0)} -> 3 paths!
    # Total = 1(hub) + 3 + 3 + 3 = 10 valid types, vs ONLY 9 REACHABLE TYPES!!
    # WHOA!! Which type is UNREACHABLE at (M=3, B=3)?! Let's find it!
    # Probably path (2,1)-(2,0)-(2,1)? Wait: founder (2,2) -> promote -> (2,1)-(2,1) -> promote -> (2,1)-(2,0)-(2,1)!
    # Then what IS unreachable at (M=3, B=3)?! (2,0)-(2,0) or (2,0)-(2,0)-(2,0)?
    # Check right now!
    op33 = recomputed[(3, 3)][0]
    k_all_valid_33 = []
    # Size 1
    for l in [0, 1, 2]:
        k_all_valid_33.append((f"1:({l})", canonical_active_tree([(2, l)], [[]])))
    # Size 2
    for l1 in [0, 1]:
        for l2 in range(l1, 2):
            k_all_valid_33.append((f"2:({l1},{l2})", canonical_active_tree([(2, l1), (2, l2)], [[1], [0]])))
    # Size 3 path
    for l1 in [0, 1]:
        for l2 in range(l1, 2):
            k_all_valid_33.append((f"3:({l1},0,{l2})", canonical_active_tree([(2, l1), (2, 0), (2, l2)], [[1], [0, 2], [1]])))
    reachable_set = set(op33["keys"])
    unreachable = [name for name, k in k_all_valid_33 if k not in reachable_set]
    print(f"  Reachability check at (3,3): unreachable valid types = {unreachable} (count={len(unreachable)}, reachable=9/10)")
    assert len(unreachable) == 1

    # Also test adding a spurious unreachable absorbing/supercritical type (e.g. a hypothetical
    # unreachable component or entry perturbation by 1e-6) is caught by certificate check!
    op31_rows = [dict(r) for r in recomputed[(3, 1)][0]["rows"]]
    win31 = recomputed[(3, 1)][2]
    op31_rows[win31[0]][win31[0]] = op31_rows[win31[0]].get(win31[0], Fraction(0)) + Fraction(1, 10**6)
    lo_pert, hi_pert, _ = certify_scc_spectral_abscissa(op31_rows, win31)
    assert lo_pert > Fraction("0.61803398874990"), "1e-6 matrix perturbation must break bracket!"
    print(f"  Mutation 5 (1e-6 diagonal perturbation on (3,1)): lam={float(lo_pert):.14f} exits bracket!")
    print("  [PASS] Test B\n")


# ---------------------------------------------------------------------------
# 6. TEST C: Direct Annotated Branching & Section 3-4 Coupling Verification
# ---------------------------------------------------------------------------

def test_sections3_4_5_coupling_and_bounds(recomputed):
    print("=== TEST C: Sections 3-5 Potential Invariant, Exact Coupling & Bounds ===")

    # C1: Simulate annotated auxiliary process from T_3, verify:
    #   (i) Section 3 invariant s*J + P <= C_tilde - 1 at EVERY jump
    #   (ii) E[C_tilde(t)] matches exact matrix exponential (e^{t G} c)_0!
    def simulate_annotated_aux(M: int, B: int, t_end: float, rng: random.Random):
        """Simulate full annotated auxiliary branching process from founder T_3.

        Returns (c_tilde, J, P, state_snapshot) where state_snapshot can be lifted to a true tree!
        """
        s = (B + (M - 1) - 1) // (M - 1)
        # Vertices 0..C-1:
        #   is_hub[u]: bool
        #   D[u], L[u]: int (valid for nonhubs)
        #   comp[u]: active component id (or -1 if hub)
        #   edges: list of (u, v, kind) where kind in {"active", "hub", "seam"}
        is_hub = [False]
        D = [2]
        L = [2]
        active_adj: List[List[int]] = [[]]
        all_adj: List[List[int]] = [[]]
        seams = set()  # pairs (min(u,v), max(u,v))
        hub_edges = set()

        def check_sec3_invariant():
            c_tot = len(is_hub)
            # Compute active components and P
            vis = [False] * c_tot
            P = 0
            for u in range(c_tot):
                if is_hub[u] or vis[u]:
                    continue
                q = deque([u])
                vis[u] = True
                sz = 0
                while q:
                    x = q.popleft()
                    sz += 1
                    for w in active_adj[x]:
                        if not vis[w]:
                            vis[w] = True
                            q.append(w)
                assert 1 <= sz <= B
                if sz > s:
                    P += sz - s
            J = len(seams)
            assert s * J + P <= c_tot - 1, f"Sec 3 invariant failed: {s}*{J}+{P} > {c_tot-1}"
            # Check slot count == inactive incident edges for every nonhub
            for u in range(c_tot):
                if not is_hub[u]:
                    inact = len(all_adj[u]) - len(active_adj[u])
                    assert D[u] - L[u] - len(active_adj[u]) == inact

        check_sec3_invariant()
        t_cur = 0.0
        while t_cur < t_end:
            c_tot = len(is_hub)
            # Compute total event rate across all vertices
            ev_rates = []
            tot_rate = 0.0
            for u in range(c_tot):
                if is_hub[u]:
                    r = 1.0
                    ev_rates.append(("hub_birth", u, r))
                    tot_rate += r
                else:
                    g_rate = float(L[u]) + sum(1.0 / D[w] for w in active_adj[u])
                    if g_rate > 0.0:
                        ev_rates.append(("gain", u, g_rate))
                        tot_rate += g_rate
                    if L[u] > 0:
                        p_rate = L[u] / D[u]
                        ev_rates.append(("prom", u, p_rate))
                        tot_rate += p_rate
            dt = rng.expovariate(tot_rate)
            if t_cur + dt >= t_end:
                break
            t_cur += dt
            pick = rng.random() * tot_rate
            kind, u = ev_rates[-1][0], ev_rates[-1][1]
            for k_ev, u_ev, r_ev in ev_rates:
                if pick < r_ev:
                    kind, u = k_ev, u_ev
                    break
                pick -= r_ev

            if kind == "hub_birth":
                nb = len(is_hub)
                is_hub.append(False)
                D.append(2)
                L.append(1)
                active_adj.append([])
                all_adj.append([u])
                all_adj[u].append(nb)
                hub_edges.add((min(u, nb), max(u, nb)))
            elif kind == "gain":
                D[u] += 1
                L[u] += 1
                if D[u] >= M:
                    is_hub[u] = True
                    for w in list(active_adj[u]):
                        active_adj[w].remove(u)
                        hub_edges.add((min(u, w), max(u, w)))
                    active_adj[u].clear()
            else:  # prom
                L[u] -= 1
                nb = len(is_hub)
                is_hub.append(False)
                D.append(2)
                L.append(1)
                active_adj.append([u])
                active_adj[u].append(nb)
                all_adj.append([u])
                all_adj[u].append(nb)
                # Find active component containing u
                comp_nodes = []
                vis = {u}
                q = deque([u])
                while q:
                    x = q.popleft()
                    comp_nodes.append(x)
                    for w in active_adj[x]:
                        if w not in vis:
                            vis.add(w)
                            q.append(w)
                if len(comp_nodes) == B + 1:
                    # Find all legal edges
                    comp_set = set(comp_nodes)
                    legal_edges = []
                    for x in comp_nodes:
                        for y in active_adj[x]:
                            if x < y:
                                # Size of x's side without (x, y)
                                seen = {x}
                                qq = deque([x])
                                while qq:
                                    z = qq.popleft()
                                    for w in active_adj[z]:
                                        if (z == x and w == y) or (z == y and w == x):
                                            continue
                                        if w not in seen:
                                            seen.add(w)
                                            qq.append(w)
                                sz1 = len(seen)
                                sz2 = len(comp_nodes) - sz1
                                if sz1 >= s and sz2 >= s:
                                    legal_edges.append((x, y))
                    assert len(legal_edges) >= 1
                    cx, cy = legal_edges[rng.randrange(len(legal_edges))]
                    active_adj[cx].remove(cy)
                    active_adj[cy].remove(cx)
                    seams.add((min(cx, cy), max(cx, cy)))

            check_sec3_invariant()

        return {
            "is_hub": is_hub,
            "D": D,
            "L": L,
            "active_adj": active_adj,
            "all_adj": all_adj,
            "seams": seams,
            "hub_edges": hub_edges,
        }

    # Verify (e^{t G} c)_0 matches MC E[C_tilde(t)] AND Section 3 holds 100% of jumps!
    def exact_etg_c0(op, t: float) -> float:
        n = len(op["rows"])
        rows_f = [[(j, float(v)) for j, v in r.items() if j != i] for i, r in enumerate(op["rows"])]
        diag_f = [float(op["rows"][i].get(i, Fraction(0))) for i in range(n)]
        # Uniformize e_0 @ exp(t G)
        c_shift = max(diag_f)
        alpha = max(c_shift - d for d in diag_f)
        if alpha == 0.0:
            return float(op["core_counts"][0]) * math.exp(c_shift * t)
        lam = alpha * t
        k_max = max(25, int(lam + 10 * math.sqrt(lam) + 40))
        p_diag = [1.0 + (diag_f[i] - c_shift) / alpha for i in range(n)]
        p_off = [[(j, r / alpha) for j, r in rows_f[i]] for i in range(n)]
        cur = [0.0] * n
        cur[0] = 1.0
        log_pref = c_shift * t - lam
        w = math.exp(log_pref)
        res = [w * x for x in cur]
        log_lam = math.log(lam)
        log_fact = 0.0
        for k in range(1, k_max + 1):
            nxt = [0.0] * n
            for i in range(n):
                if cur[i] == 0.0:
                    continue
                nxt[i] += cur[i] * p_diag[i]
                for j, pr in p_off[i]:
                    nxt[j] += cur[i] * pr
            cur = nxt
            log_fact += math.log(k)
            w = math.exp(log_pref + k * log_lam - log_fact)
            for i in range(n):
                res[i] += w * cur[i]
        return sum(res[i] * op["core_counts"][i] for i in range(n))

    rng = random.Random(2026)
    for M, B, t_chk in [(3, 3, 1.0), (4, 2, 1.0)]:
        op = recomputed[(M, B)][0]
        exact_mean = exact_etg_c0(op, t_chk)
        n_mc = 4000
        vals = []
        for _ in range(n_mc):
            snap = simulate_annotated_aux(M, B, t_chk, rng)
            vals.append(len(snap["is_hub"]))
        mc_mean = sum(vals) / n_mc
        mc_se = math.sqrt(sum((v - mc_mean) ** 2 for v in vals) / (n_mc * (n_mc - 1)))
        print(f"  Sec 3 & exact E[C_tilde({t_chk})] at ({M},{B}): "
              f"MC={mc_mean:.4f}+/-{mc_se:.4f} vs (e^{{tG}}c)_0={exact_mean:.4f} (0 invariant violations)")
        assert abs(mc_mean - exact_mean) < 4 * mc_se

    # C2: EXACT Section 4 Coupling Simulation (starting from a seam-rich, hub-rich reachable state!)
    # Implement EVERY rule of Section 4 literally and assert ALL clean-pair invariants after EVERY jump,
    # verify E[I(t)] / C_0 <= d_{M, B}(t), and verify a planted coupling bug BREAKS invariants!
    def run_section4_coupling_window(
        M: int,
        B: int,
        init_snap: dict,
        t_win: float,
        rng: random.Random,
        *,
        plant_bug_ignore_hub_recipient: bool = False,
    ) -> Tuple[int, int, int]:
        """Run exact Section 4 coupling on [0, t_win] from a true lift of init_snap.
        Returns (C_true(t_win), C_aux(t_win), I(t_win)).
        """
        s = (B + (M - 1) - 1) // (M - 1)
        c0 = len(init_snap["is_hub"])

        # True tree cores: indexed by integer IDs 0..c0-1 initially (and new IDs for unpaired true cores)
        # Aux cores: indexed by integer IDs 0..c0-1 initially (and new IDs for unpaired aux cores)
        # Paired vertices share the SAME integer ID! Unpaired cores get distinct fresh IDs!
        next_id = c0
        true_cores = set(range(c0))
        aux_cores = set(range(c0))
        clean = set(range(c0))  # IDs that are currently clean pairs (in both true_cores and aux_cores)

        true_adj: Dict[int, List[int]] = {u: list(init_snap["all_adj"][u]) for u in range(c0)}
        true_L: Dict[int, int] = {}
        aux_hub: Dict[int, bool] = {u: init_snap["is_hub"][u] for u in range(c0)}
        aux_D: Dict[int, int] = {u: init_snap["D"][u] for u in range(c0)}
        aux_L: Dict[int, int] = {u: init_snap["L"][u] for u in range(c0)}
        aux_active_adj: Dict[int, List[int]] = {u: list(init_snap["active_adj"][u]) for u in range(c0)}
        aux_all_adj: Dict[int, List[int]] = {u: list(init_snap["all_adj"][u]) for u in range(c0)}
        seams = set(init_snap["seams"])
        hub_edges = set(init_snap["hub_edges"])

        # True lift (Section 4, line 42):
        # nonhub -> stored L_u; hub -> max(0, M - c_u) leaves
        for u in range(c0):
            if aux_hub[u]:
                cu = len(true_adj[u])
                true_L[u] = max(0, M - cu)
            else:
                true_L[u] = aux_L[u]

        def true_D(u: int) -> int:
            return len(true_adj[u]) + true_L[u]

        def check_clean_invariants():
            for u in clean:
                assert u in true_cores and u in aux_cores
                # No clean vertex has an unmatched core neighbor!
                assert set(true_adj[u]) == set(aux_all_adj[u]), (
                    f"Clean vertex {u} core neighbors differ: true={true_adj[u]} aux={aux_all_adj[u]}"
                )
                if aux_hub[u]:
                    assert true_D(u) >= M, f"Clean hub {u} has true_D={true_D(u)} < {M}"
                else:
                    assert true_D(u) == aux_D[u] < M, f"Clean low {u} D mismatch: {true_D(u)} vs {aux_D[u]}"
                    assert true_L[u] == aux_L[u], f"Clean low {u} L mismatch: {true_L[u]} vs {aux_L[u]}"
            i_true = len(true_cores - clean)
            i_aux = len(aux_cores - clean)
            assert abs(len(true_cores) - len(aux_cores)) <= i_true + i_aux

        def maybe_size_cut(start_u: int):
            if aux_hub[start_u]:
                return
            comp = []
            vis = {start_u}
            q = deque([start_u])
            while q:
                x = q.popleft()
                comp.append(x)
                for w in aux_active_adj[x]:
                    if w not in vis:
                        vis.add(w)
                        q.append(w)
            if len(comp) == B + 1:
                legal = []
                for x in comp:
                    for y in aux_active_adj[x]:
                        if x < y:
                            seen = {x}
                            qq = deque([x])
                            while qq:
                                z = qq.popleft()
                                for w in aux_active_adj[z]:
                                    if (z == x and w == y) or (z == y and w == x):
                                        continue
                                    if w not in seen:
                                        seen.add(w)
                                        qq.append(w)
                            if len(seen) >= s and len(comp) - len(seen) >= s:
                                legal.append((x, y))
                assert legal
                cx, cy = legal[rng.randrange(len(legal))]
                aux_active_adj[cx].remove(cy)
                aux_active_adj[cy].remove(cx)
                seams.add((min(cx, cy), max(cx, cy)))

        def promote_aux_hub(u: int):
            aux_hub[u] = True
            for w in list(aux_active_adj[u]):
                aux_active_adj[w].remove(u)
                hub_edges.add((min(u, w), max(u, w)))
            aux_active_adj[u].clear()

        check_clean_invariants()
        t_cur = 0.0
        while t_cur < t_win:
            # Build exact list of independent Poisson channels under Section 4 coupling
            events = []
            tot_rate = 0.0
            # 1. Clean pairs
            for u in clean:
                if aux_hub[u]:
                    # Clean hub: common rate-1 ring + true-only leaf births (rate true_L[u])
                    events.append(("clean_hub_ring", u, 1.0))
                    tot_rate += 1.0
                    if true_L[u] > 0:
                        events.append(("clean_hub_leaf", u, float(true_L[u])))
                        tot_rate += float(true_L[u])
                else:
                    # Clean low: synchronized leaf births (rate L_u) + common rate-1 core ring
                    if aux_L[u] > 0:
                        events.append(("clean_low_leaf", u, float(aux_L[u])))
                        tot_rate += float(aux_L[u])
                    events.append(("clean_low_ring", u, 1.0))
                    tot_rate += 1.0
            # 2. Dirty true cores (individual rates)
            for u in (true_cores - clean):
                if true_L[u] > 0:
                    events.append(("dirty_true_leaf", u, float(true_L[u])))
                    tot_rate += float(true_L[u])
                events.append(("dirty_true_ring", u, 1.0))
                tot_rate += 1.0
            # 3. Dirty aux cores (individual rates)
            for u in (aux_cores - clean):
                if aux_hub[u]:
                    events.append(("dirty_aux_hub", u, 1.0))
                    tot_rate += 1.0
                else:
                    if aux_L[u] > 0:
                        events.append(("dirty_aux_leaf", u, float(aux_L[u])))
                        tot_rate += float(aux_L[u])
                    # Active sends + promotions (total rate (L_u + a_u)/D_u)
                    r_core = (aux_L[u] + len(aux_active_adj[u])) / aux_D[u]
                    if r_core > 0.0:
                        events.append(("dirty_aux_ring", u, r_core))
                        tot_rate += r_core

            dt = rng.expovariate(tot_rate)
            if t_cur + dt >= t_win:
                break
            t_cur += dt
            pick = rng.random() * tot_rate
            ev_kind, u = events[-1][0], events[-1][1]
            for k_ev, u_ev, r_ev in events:
                if pick < r_ev:
                    ev_kind, u = k_ev, u_ev
                    break
                pick -= r_ev

            if ev_kind == "clean_low_leaf":
                true_L[u] += 1
                aux_L[u] += 1
                aux_D[u] += 1
                if aux_D[u] >= M:
                    promote_aux_hub(u)
            elif ev_kind == "clean_low_ring":
                du = aux_D[u]
                lu = aux_L[u]
                idx_choice = rng.randrange(du)
                if idx_choice < lu:
                    # Synchronized promotion!
                    true_L[u] -= 1
                    aux_L[u] -= 1
                    nb = next_id
                    next_id += 1
                    true_cores.add(nb)
                    aux_cores.add(nb)
                    clean.add(nb)
                    true_adj[nb] = [u]
                    true_adj[u].append(nb)
                    true_L[nb] = 1
                    aux_hub[nb] = False
                    aux_D[nb] = 2
                    aux_L[nb] = 1
                    aux_active_adj[nb] = [u]
                    aux_active_adj[u].append(nb)
                    aux_all_adj[nb] = [u]
                    aux_all_adj[u].append(nb)
                    maybe_size_cut(u)
                else:
                    # Pick paired core neighbor w
                    w = true_adj[u][idx_choice - lu]
                    e_key = (min(u, w), max(u, w))
                    # True ALWAYS sends a leaf to w
                    true_L[w] += 1
                    if w in aux_active_adj[u]:
                        # Active edge: aux ALSO sends leaf to w
                        aux_L[w] += 1
                        aux_D[w] += 1
                        if aux_D[w] >= M:
                            promote_aux_hub(w)
                    elif e_key in seams:
                        # Seam: aux sends nothing; if w is clean low, mark w dirty!
                        if w in clean and (not aux_hub[w]):
                            clean.discard(w)
                    else:
                        # Non-seam inactive edge -> goes to hub in aux genealogy
                        assert aux_hub[w]
            elif ev_kind == "clean_hub_leaf":
                true_L[u] += 1
            elif ev_kind == "clean_hub_ring":
                du = true_D(u)
                lu = true_L[u]
                idx_choice = rng.randrange(du)
                if idx_choice < lu:
                    # True picks leaf -> paired newborn!
                    true_L[u] -= 1
                    nb = next_id
                    next_id += 1
                    true_cores.add(nb)
                    aux_cores.add(nb)
                    clean.add(nb)
                    true_adj[nb] = [u]
                    true_adj[u].append(nb)
                    true_L[nb] = 1
                    aux_hub[nb] = False
                    aux_D[nb] = 2
                    aux_L[nb] = 1
                    aux_active_adj[nb] = []
                    aux_all_adj[nb] = [u]
                    aux_all_adj[u].append(nb)
                    hub_edges.add((min(u, nb), max(u, nb)))
                else:
                    # True picks core neighbor w -> unmatched aux newborn + true leaf at w
                    w = true_adj[u][idx_choice - lu]
                    true_L[w] += 1
                    nb = next_id
                    next_id += 1
                    aux_cores.add(nb)
                    aux_hub[nb] = False
                    aux_D[nb] = 2
                    aux_L[nb] = 1
                    aux_active_adj[nb] = []
                    aux_all_adj[nb] = [u]
                    aux_all_adj[u].append(nb)
                    hub_edges.add((min(u, nb), max(u, nb)))
                    # Mark source u dirty, newborn nb is already not in clean, and mark recipient w if clean!
                    clean.discard(u)
                    if not plant_bug_ignore_hub_recipient:
                        clean.discard(w)
            elif ev_kind == "dirty_true_leaf":
                true_L[u] += 1
            elif ev_kind == "dirty_true_ring":
                du = true_D(u)
                lu = true_L[u]
                idx_choice = rng.randrange(du)
                if idx_choice < lu:
                    true_L[u] -= 1
                    nb = next_id
                    next_id += 1
                    true_cores.add(nb)
                    true_adj[nb] = [u]
                    true_adj[u].append(nb)
                    true_L[nb] = 1
                else:
                    w = true_adj[u][idx_choice - lu]
                    true_L[w] += 1
                    clean.discard(w)
            elif ev_kind == "dirty_aux_hub":
                nb = next_id
                next_id += 1
                aux_cores.add(nb)
                aux_hub[nb] = False
                aux_D[nb] = 2
                aux_L[nb] = 1
                aux_active_adj[nb] = []
                aux_all_adj[nb] = [u]
                aux_all_adj[u].append(nb)
                hub_edges.add((min(u, nb), max(u, nb)))
            elif ev_kind == "dirty_aux_leaf":
                aux_L[u] += 1
                aux_D[u] += 1
                if aux_D[u] >= M:
                    promote_aux_hub(u)
            elif ev_kind == "dirty_aux_ring":
                lu = aux_L[u]
                anbrs = aux_active_adj[u]
                idx_choice = rng.randrange(lu + len(anbrs))
                if idx_choice < lu:
                    aux_L[u] -= 1
                    nb = next_id
                    next_id += 1
                    aux_cores.add(nb)
                    aux_hub[nb] = False
                    aux_D[nb] = 2
                    aux_L[nb] = 1
                    aux_active_adj[nb] = [u]
                    aux_active_adj[u].append(nb)
                    aux_all_adj[nb] = [u]
                    aux_all_adj[u].append(nb)
                    maybe_size_cut(u)
                else:
                    w = anbrs[idx_choice - lu]
                    aux_L[w] += 1
                    aux_D[w] += 1
                    if aux_D[w] >= M:
                        promote_aux_hub(w)
                    clean.discard(w)

            check_clean_invariants()

        i_end = len(true_cores - clean) + len(aux_cores - clean)
        return len(true_cores), len(aux_cores), i_end

    # Generate a reachable seam-rich & hub-rich annotated state for (M=4, B=3) and (M=6, B=4)
    for M, B in [(4, 3), (6, 4)]:
        s = (B + (M - 1) - 1) // (M - 1)
        # Run until we get a reachable state with >= 1 hub and >= 1 seam
        snap = None
        for seed_try in range(200):
            cand = simulate_annotated_aux(M, B, 2.5, random.Random(1000 + seed_try))
            if any(cand["is_hub"]) and len(cand["seams"]) >= 1:
                snap = cand
                break
        assert snap is not None
        c0 = len(snap["is_hub"])
        t_win = 0.5
        d_mb = (10.0 / M + 2.0 / s) * (math.exp(2.0 * t_win) - math.exp(t_win))
        n_trials = 1000
        sum_i = 0.0
        sum_ct = 0.0
        sum_ca = 0.0
        rng_c = random.Random(777 + M)
        for _ in range(n_trials):
            ct, ca, ie = run_section4_coupling_window(M, B, snap, t_win, rng_c)
            sum_ct += ct
            sum_ca += ca
            sum_i += ie
        ei_norm = (sum_i / n_trials) / c0
        diff_norm = abs(sum_ca - sum_ct) / (n_trials * c0)
        print(f"  Sec 4 coupling (M={M}, B={B}, s={s}, C0={c0}, hubs={sum(snap['is_hub'])}, "
              f"seams={len(snap['seams'])}, t={t_win}): "
              f"|E C_aux - E C_true|/C0={diff_norm:.4f} <= E[I]/C0={ei_norm:.4f} <= d_MB={d_mb:.4f}")
        assert ei_norm <= d_mb

        # Planted error on Section 4 coupling: ignore marking clean recipient dirty on hub core-hit
        if M == 4 and B == 3:
            caught = False
            for seed_bad in range(500):
                try:
                    run_section4_coupling_window(
                        M, B, snap, 1.5, random.Random(5000 + seed_bad),
                        plant_bug_ignore_hub_recipient=True,
                    )
                except AssertionError:
                    caught = True
                    break
            assert caught, "Planted coupling error (unmarked hub recipient) must violate invariants!"
            print("  Planted coupling error (omit hub recipient mark) caught by invariant assertion!")

    # C3: Section 5 analytic constant checks (lines 85-90)
    for M in [24**3, 24**3 + 1, 100000, 10**7]:
        B = M * M
        s = (B + (M - 1) - 1) // (M - 1)
        t = math.log(M) / 3.0
        d = (10.0 / M + 2.0 / s) * (math.exp(2.0 * t) - math.exp(t))
        assert s >= M
        assert d <= 12.0 * (M ** (-1.0 / 3.0)) + 1e-15
        assert d <= 0.5 + 1e-15
    # Verify (U + d)/(L - d) <= 3 * U / L for all U, L >= 1 and 0 <= d <= 1/2
    for U in [1.0, 2.5, 100.0]:
        for L in [1.0, 1.5, 50.0]:
            for d in [0.0, 0.25, 0.5]:
                assert (U + d) / (L - d) <= 3.0 * (U / L) + 1e-12
    print("  Section 5 scaling (B=M^2, t=log(M)/3, M>=24^3 -> d<=1/2, loss<=log(3)/t) verified.")
    print("  [PASS] Test C\n")


def main():
    recomputed = test_all_six_matrices_and_cert()
    test_matrix_mutations(recomputed)
    test_sections3_4_5_coupling_and_bounds(recomputed)
    print("ALL FINITE-OPERATOR VERIFICATION, MUTATION & COUPLING SUITES PASSED.")


if __name__ == "__main__":
    main()
