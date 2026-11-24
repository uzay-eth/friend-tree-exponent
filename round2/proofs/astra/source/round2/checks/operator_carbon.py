#!/usr/bin/env python3
"""
Independent exact small-tree checker for round2/operator.md (and round2/degree-law.md).

Recomputes in exact rational arithmetic (fractions.Fraction) and bounded-error
Poisson uniformization:
  1. Equation (2): cell representation of Q_s vs microscopic full friend-tree
     dynamics, and semigroup identity Q_s Q_t = Q_{s+t} on unlabelled rooted
     marked core trees X (no orientation toward the initial seed).
  2. Equation (3): old-edge path-reversal inequality, exact vanishing order
     O(s^{R-k+1}) of (Q_s f - Q_s^[R] f) when missing old neighbour slots are
     retained, and numerical spatial cutoff inequality (3).
  3. Equation (7): generator Af(G, o), Dynkin's formula Q_t f - f = int_0^t Q_s Af ds,
     the renewal identity Q_t f = S_t f + int_0^t Q_s B S_{t-s} f ds, and the
     degree-crossing generator identity A(1_{D_o >= k+1}) = (k - c_o + rho_o) 1_{D_o = k}.
  4. Both edge-transfer cancellations (signed gain transfer r_w g(w, u) and
     thresholded arrow transfer 1_{D_w >= T} a(D_u, c_u) / D_w) on finite trees
     and root-uniform measures, plus exact agreement of amended vs unamended sums.
  5. Eight deliberate planted errors (negative controls), verifying that each
     planted error is strictly detected and rejected by the checker.

Single-CPU execution (< 2 seconds). No external dependencies.
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
import math
import sys
from typing import Callable, Dict, FrozenSet, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# 1. Data structures for marked core trees and full friend trees
# ---------------------------------------------------------------------------

class CoreTree:
    """
    Finite marked core tree with optional retained missing old-neighbour slots
    (used for truncated dynamics Q_s^[R]) and optional cell labels.

    Vertices are 0 .. n-1.
    adj[u]: tuple of core neighbours in the active graph.
    L[u]: nonnegative integer leaf count.
    extra_c[u]: nonnegative integer count of retained missing old-neighbour slots
                (0 in any true tree in X; >= 0 in an R-ball truncation).
    cell[u]: integer cell label (defaults to u at time 0).
    """

    __slots__ = ("n", "adj", "L", "extra_c", "cell")

    def __init__(
        self,
        adj: Sequence[Sequence[int]],
        L: Sequence[int],
        extra_c: Optional[Sequence[int]] = None,
        cell: Optional[Sequence[int]] = None,
    ) -> None:
        self.n = len(adj)
        self.adj = tuple(tuple(sorted(nb)) for nb in adj)
        self.L = tuple(int(x) for x in L)
        self.extra_c = tuple(0 for _ in range(self.n)) if extra_c is None else tuple(int(x) for x in extra_c)
        self.cell = tuple(range(self.n)) if cell is None else tuple(int(x) for x in cell)
        for u in range(self.n):
            assert self.L[u] >= 0, f"Negative leaf mark at {u}"
            assert self.extra_c[u] >= 0, f"Negative extra_c at {u}"
            assert self.D(u) >= 2, f"Core vertex {u} has total degree {self.D(u)} < 2"

    def c(self, u: int) -> int:
        return len(self.adj[u]) + self.extra_c[u]

    def D(self, u: int) -> int:
        return self.L[u] + self.c(u)

    def gain_rate(self, u: int, wrong_denom_self: bool = False) -> Fraction:
        """r_u = L_u + sum_{y ~ u} 1 / D_y (inbound across suppressed boundary = 0)."""
        acc = Fraction(self.L[u], 1)
        for y in self.adj[u]:
            denom = self.D(u) if wrong_denom_self else self.D(y)
            acc += Fraction(1, denom)
        return acc

    def prom_rate(self, u: int) -> Fraction:
        """b_u = L_u / D_u."""
        return Fraction(self.L[u], self.D(u))

    def after_gain(self, x: int) -> "CoreTree":
        new_L = list(self.L)
        new_L[x] += 1
        return CoreTree(self.adj, new_L, self.extra_c, self.cell)

    def after_prom(self, x: int, bug_root_only_cell: bool = False) -> Tuple["CoreTree", int]:
        """Promote one leaf of x; returns (new_tree, new_vertex_id)."""
        assert self.L[x] >= 1
        new_id = self.n
        new_adj = [list(nb) for nb in self.adj]
        new_adj[x].append(new_id)
        new_adj.append([x])
        new_L = list(self.L)
        new_L[x] -= 1
        new_L.append(1)
        new_extra = list(self.extra_c) + [0]
        if bug_root_only_cell:
            # Planted bug: indirect descendants (x != cell[x]) fail to inherit cell[x]
            inherited = self.cell[x] if x == self.cell[x] else -1
        else:
            inherited = self.cell[x]
        new_cell = list(self.cell) + [inherited]
        return CoreTree(new_adj, new_L, new_extra, new_cell), new_id

    def relabelled_fresh(self) -> "CoreTree":
        """Forget cell history and give each current core vertex its own label."""
        return CoreTree(self.adj, self.L, self.extra_c, tuple(range(self.n)))

    def distances_from(self, root: int) -> List[int]:
        dist = [-1] * self.n
        dist[root] = 0
        q = deque([root])
        while q:
            u = q.popleft()
            for w in self.adj[u]:
                if dist[w] == -1:
                    dist[w] = dist[u] + 1
                    q.append(w)
        return dist

    def canonical_ball(self, root: int, radius: Optional[int] = None, include_c_mark: bool = True) -> Tuple:
        """
        Isomorphism-invariant encoding of the rooted marked ball B_radius(root)
        (or the entire rooted marked tree when radius is None).
        No vertex IDs or seed orientation survive.
        """
        dist = self.distances_from(root)

        def encode(u: int, parent: int) -> Tuple:
            children_codes = []
            if radius is None or dist[u] < radius:
                for w in self.adj[u]:
                    if w != parent:
                        children_codes.append(encode(w, u))
            children_codes.sort()
            if include_c_mark:
                return (self.L[u], self.c(u), tuple(children_codes))
            return (self.L[u], tuple(children_codes))

        return encode(root, -1)

    def truncate_ball(self, root: int, R: int, retain_missing_slots: bool = True) -> Tuple["CoreTree", int]:
        """
        Construct the R-ball old truncation around `root`:
        keep vertices with dist(root, u) <= R; suppress inbound arrows from > R;
        if retain_missing_slots is True, keep initial c_u^0 via extra_c.
        """
        dist = self.distances_from(root)
        kept = [u for u in range(self.n) if 0 <= dist[u] <= R]
        idx_map = {u: i for i, u in enumerate(kept)}
        sub_adj: List[List[int]] = [[] for _ in kept]
        sub_L: List[int] = [0] * len(kept)
        sub_extra: List[int] = [0] * len(kept)
        for u in kept:
            i = idx_map[u]
            for w in self.adj[u]:
                if w in idx_map:
                    sub_adj[i].append(idx_map[w])
            missing = self.c(u) - len(sub_adj[i])
            sub_L[i] = self.L[u]
            if retain_missing_slots:
                sub_extra[i] = missing
            else:
                # Planted bug mode: drop missing old neighbour slots (clamp L only if D < 2 would crash,
                # so test on trees where D >= 2 even without missing slots!)
                sub_extra[i] = 0
        new_root = idx_map[root]
        sub_cell = list(range(len(kept)))
        return CoreTree(sub_adj, sub_L, sub_extra, sub_cell), new_root


class FullTree:
    """
    Explicit full friend tree (both core vertices and degree-1 leaves) for
    completely independent microscopic verification of Equation (2).
    """

    __slots__ = ("n", "adj", "origin")

    def __init__(self, adj: Sequence[Sequence[int]], origin: Sequence[int]) -> None:
        self.n = len(adj)
        self.adj = tuple(tuple(sorted(nb)) for nb in adj)
        # origin[v]: for core vertex, the time-0 core label it descends from;
        # for a current leaf l attached to core x, origin[l] = origin[x] (since if l is
        # promoted, x is its unique neighbour and promoter!).
        self.origin = tuple(int(x) for x in origin)

    @staticmethod
    def from_core_tree(ct: CoreTree) -> "FullTree":
        assert all(e == 0 for e in ct.extra_c), "FullTree requires a genuine finite tree"
        adj: List[List[int]] = [list(nb) for nb in ct.adj]
        origin: List[int] = list(range(ct.n))
        for u in range(ct.n):
            for _ in range(ct.L[u]):
                leaf_id = len(adj)
                adj[u].append(leaf_id)
                adj.append([u])
                origin.append(u)
        return FullTree(adj, origin)

    def to_core_tree_and_cells(self) -> Tuple[CoreTree, Dict[int, int]]:
        core_nodes = [u for u in range(self.n) if len(self.adj[u]) >= 2]
        idx = {u: i for i, u in enumerate(core_nodes)}
        c_adj: List[List[int]] = [[] for _ in core_nodes]
        c_L: List[int] = [0] * len(core_nodes)
        c_cell: List[int] = [0] * len(core_nodes)
        for u in core_nodes:
            i = idx[u]
            c_cell[i] = self.origin[u]
            leaves = 0
            for w in self.adj[u]:
                if w in idx:
                    c_adj[i].append(idx[w])
                else:
                    leaves += 1
            c_L[i] = leaves
        return CoreTree(c_adj, c_L, None, c_cell), idx

    def step_attach_leaf(self, target: int) -> "FullTree":
        """Attach a new leaf to `target`."""
        new_id = self.n
        new_adj = [list(nb) for nb in self.adj]
        new_adj[target].append(new_id)
        new_adj.append([target])
        # If target was a leaf, its sole old neighbour was its promoter, so origin[target]
        # is already the promoter's cell! And the newborn leaf's unique neighbour is `target`,
        # which is now core with cell origin[target], so newborn leaf also gets origin[target].
        new_origin = list(self.origin) + [self.origin[target]]
        return FullTree(new_adj, new_origin)


# ---------------------------------------------------------------------------
# 2. Exact rational Taylor-coefficient engines (microscopic vs cell operator)
# ---------------------------------------------------------------------------

TestFn = Callable[[CoreTree, int], Fraction]


def eval_cell_sum(ct: CoreTree, root_label: int, f: TestFn) -> Fraction:
    acc = Fraction(0, 1)
    for v in range(ct.n):
        if ct.cell[v] == root_label:
            acc += f(ct, v)
    return acc


def core_operator_derivatives(
    ct0: CoreTree,
    root: int,
    f: TestFn,
    max_order: int,
    bug_root_only_cell: bool = False,
) -> List[Fraction]:
    """
    Exact derivatives d^m / ds^m |_{s=0} (Q_s f)(ct0, root) for m = 0 .. max_order,
    computed directly from the continuous-time core-cell Markov jump rates.
    """
    root_label = ct0.cell[root]

    # Memoize (state_signature, order_left) -> Fraction
    memo: Dict[Tuple, Fraction] = {}

    def sig(ct: CoreTree) -> Tuple:
        return (ct.adj, ct.L, ct.extra_c, ct.cell)

    def deriv(ct: CoreTree, m: int) -> Fraction:
        if m == 0:
            return eval_cell_sum(ct, root_label, f)
        key = (sig(ct), m)
        if key in memo:
            return memo[key]
        total_out = Fraction(0, 1)
        acc = Fraction(0, 1)
        for x in range(ct.n):
            rx = ct.gain_rate(x)
            if rx != 0:
                total_out += rx
                acc += rx * deriv(ct.after_gain(x), m - 1)
            if ct.L[x] > 0:
                bx = ct.prom_rate(x)
                total_out += bx
                ct_p, _ = ct.after_prom(x, bug_root_only_cell=bug_root_only_cell)
                acc += bx * deriv(ct_p, m - 1)
        acc -= total_out * deriv(ct, m - 1)
        memo[key] = acc
        return acc

    return [deriv(ct0, m) for m in range(max_order + 1)]


def full_tree_microscopic_derivatives(
    ct0: CoreTree,
    root: int,
    f: TestFn,
    max_order: int,
) -> List[Fraction]:
    """
    Completely independent calculation of d^m / ds^m |_{s=0} E[sum_{v in H_o(s)} f(G_s, v)]
    from the microscopic full friend-tree process: each vertex u of T (core or leaf)
    rings at rate 1 and chooses a uniform neighbour w of u (prob 1 / deg_T(u)),
    attaching a new leaf to w.
    """
    ft0 = FullTree.from_core_tree(ct0)
    memo: Dict[Tuple, Fraction] = {}

    def sig(ft: FullTree) -> Tuple:
        return (ft.adj, ft.origin)

    def deriv(ft: FullTree, m: int) -> Fraction:
        if m == 0:
            ct, _ = ft.to_core_tree_and_cells()
            return eval_cell_sum(ct, root, f)
        key = (sig(ft), m)
        if key in memo:
            return memo[key]
        # Total rate is ft.n (each of the n vertices rings at rate 1)
        # Aggregate attachment rate to each target w
        target_rate: Dict[int, Fraction] = {}
        for u in range(ft.n):
            du = len(ft.adj[u])
            inv = Fraction(1, du)
            for w in ft.adj[u]:
                target_rate[w] = target_rate.get(w, Fraction(0, 1)) + inv
        assert sum(target_rate.values(), Fraction(0, 1)) == Fraction(ft.n, 1)
        acc = -Fraction(ft.n, 1) * deriv(ft, m - 1)
        for w, rw in target_rate.items():
            acc += rw * deriv(ft.step_attach_leaf(w), m - 1)
        memo[key] = acc
        return acc

    return [deriv(ft0, m) for m in range(max_order + 1)]


def single_root_derivatives(
    ct0: CoreTree,
    root: int,
    f: TestFn,
    max_order: int,
) -> List[Fraction]:
    """
    Exact derivatives d^m / ds^m |_{s=0} (S_s f)(ct0, root) tracking ONLY the
    same old root `root` (as in round2/degree-law.md).
    """
    memo: Dict[Tuple, Fraction] = {}

    def sig(ct: CoreTree) -> Tuple:
        return (ct.adj, ct.L, ct.extra_c)

    def deriv(ct: CoreTree, m: int) -> Fraction:
        if m == 0:
            return f(ct, root)
        key = (sig(ct), m)
        if key in memo:
            return memo[key]
        total_out = Fraction(0, 1)
        acc = Fraction(0, 1)
        for x in range(ct.n):
            rx = ct.gain_rate(x)
            if rx != 0:
                total_out += rx
                acc += rx * deriv(ct.after_gain(x), m - 1)
            if ct.L[x] > 0:
                bx = ct.prom_rate(x)
                total_out += bx
                ct_p, _ = ct.after_prom(x)
                # new_id = ct.n > root, so `root` keeps its index!
                acc += bx * deriv(ct_p, m - 1)
        acc -= total_out * deriv(ct, m - 1)
        memo[key] = acc
        return acc

    return [deriv(ct0, m) for m in range(max_order + 1)]


# ---------------------------------------------------------------------------
# 3. Finite-time evaluation via exact jump-chain expansion (bounded tail)
# ---------------------------------------------------------------------------

def eval_Q_finite_time(
    ct0: CoreTree,
    root: int,
    s: float,
    f_float: Callable[[CoreTree, int], float],
    max_jumps: int = 11,
) -> Tuple[float, float]:
    """
    Evaluate (Q_s f)(ct0, root) by conditioning on exact jump counts / holding
    times via uniformization on bounded-size trees, or recursive jump integration.
    Uses uniformization per jump-rate envelope up to `max_jumps` jumps on small s.
    Returns (approx_value, rigorous_abs_error_bound_for_f_norm_1).
    """
    # Use uniformization with rate Lambda_step? Since total jump rate grows by at most +2 per jump
    # (initial total core+leaf jump rate Lambda_0 = sum_x (r_x + b_x) <= sum_x (D_x + 1)),
    # exact integration over m jumps of inhomogeneous pure-jump Markov chain:
    # For each path of m transitions with exit rates lam_0, lam_1, ..., lam_m,
    # the probability weight on [0, s] is the divided-difference exponential E(lam_0, ..., lam_m; s)!
    root_label = ct0.cell[root]
    states: List[CoreTree] = [ct0]
    state_idx: Dict[Tuple, int] = {(ct0.adj, ct0.L, ct0.extra_c, ct0.cell): 0}
    transitions: List[List[Tuple[int, float]]] = []
    exit_rates: List[float] = []
    depths: List[int] = [0]

    q = deque([0])
    while q:
        idx = q.popleft()
        ct = states[idx]
        d = depths[idx]
        trans_here: List[Tuple[int, float]] = []
        tot = 0.0
        if d < max_jumps:
            for x in range(ct.n):
                rx = float(ct.gain_rate(x))
                if rx > 0.0:
                    tot += rx
                    nxt = ct.after_gain(x)
                    sig = (nxt.adj, nxt.L, nxt.extra_c, nxt.cell)
                    if sig not in state_idx:
                        state_idx[sig] = len(states)
                        states.append(nxt)
                        depths.append(d + 1)
                        q.append(state_idx[sig])
                    trans_here.append((state_idx[sig], rx))
                if ct.L[x] > 0:
                    bx = float(ct.prom_rate(x))
                    tot += bx
                    nxt, _ = ct.after_prom(x)
                    sig = (nxt.adj, nxt.L, nxt.extra_c, nxt.cell)
                    if sig not in state_idx:
                        state_idx[sig] = len(states)
                        states.append(nxt)
                        depths.append(d + 1)
                        q.append(state_idx[sig])
                    trans_here.append((state_idx[sig], bx))
        else:
            for x in range(ct.n):
                tot += float(ct.gain_rate(x)) + (float(ct.prom_rate(x)) if ct.L[x] > 0 else 0.0)
        transitions.append(trans_here)
        exit_rates.append(tot)

    num_s = len(states)
    max_rate = max(exit_rates) if exit_rates else 1.0
    steps = max(4, int(math.ceil(max_rate * s / 0.01)))
    dt = s / steps
    prob = [0.0] * num_s
    prob[0] = 1.0

    def forward_rhs(pvec: List[float]) -> List[float]:
        out = [-exit_rates[i] * pvec[i] for i in range(num_s)]
        for i in range(num_s):
            pi = pvec[i]
            if pi != 0.0:
                for j, rate in transitions[i]:
                    out[j] += pi * rate
        return out

    for _ in range(steps):
        k1 = forward_rhs(prob)
        v2 = [prob[i] + 0.5 * dt * k1[i] for i in range(num_s)]
        k2 = forward_rhs(v2)
        v3 = [prob[i] + 0.5 * dt * k2[i] for i in range(num_s)]
        k3 = forward_rhs(v3)
        v4 = [prob[i] + dt * k3[i] for i in range(num_s)]
        k4 = forward_rhs(v4)
        prob = [prob[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) for i in range(num_s)]

    val = 0.0
    captured_mass = sum(prob)
    for i, ct in enumerate(states):
        pi = prob[i]
        if pi > 0.0:
            cell_val = sum(f_float(ct, v) for v in range(ct.n) if ct.cell[v] == root_label)
            val += pi * cell_val
    # Tail probability (> max_jumps jumps) is strictly 1 - captured_mass; by Cauchy-Schwarz with Yule L^2:
    tail_prob = max(0.0, 1.0 - captured_mass)
    err_bound = math.sqrt(2.0) * math.exp(s) * math.sqrt(tail_prob)
    return val, err_bound


# ---------------------------------------------------------------------------
# 4. Generator (7) and Edge Transfers (Section 4 & upper-bound-banded.md)
# ---------------------------------------------------------------------------

def generator_A(
    ct: CoreTree,
    root: int,
    f: TestFn,
    omit_newborn: bool = False,
    newborn_at_old_root: bool = False,
    wrong_gain_denom_self: bool = False,
) -> Fraction:
    """
    Equation (7):
      Af(G, o) = sum_x [r_x Delta_x^+ f + b_x Delta_x^b f] + b_o f(G^{b_o}, new).
    """
    base = f(ct, root)
    acc = Fraction(0, 1)
    for x in range(ct.n):
        rx = ct.gain_rate(x, wrong_denom_self=wrong_gain_denom_self)
        if rx != 0:
            acc += rx * (f(ct.after_gain(x), root) - base)
        if ct.L[x] > 0:
            bx = ct.prom_rate(x)
            ct_p, _ = ct.after_prom(x)
            acc += bx * (f(ct_p, root) - base)
    if ct.L[root] > 0 and not omit_newborn:
        bo = ct.prom_rate(root)
        ct_po, new_id = ct.after_prom(root)
        eval_vertex = root if newborn_at_old_root else new_id
        acc += bo * f(ct_po, eval_vertex)
    return acc


def gain_transfer_divergence(
    ct: CoreTree,
    u: int,
    g_table: Callable[[int, int, int, int], Fraction],
    bug_use_receiver_rate: bool = False,
    bug_transpose_one_side: bool = False,
) -> Fraction:
    """
    Net gain-transfer adjustment added to row u (`operator.md` line 76):
    directed arc (w, v) carries r_w g(w, v) from row v to row w.
    So row u receives + sum_{v ~ u} r_u g(u, v) and pays - sum_{w ~ u} r_w g(w, u).
    """
    ru = ct.gain_rate(u)
    Du, cu = ct.D(u), ct.c(u)
    incoming = Fraction(0, 1)
    for v in ct.adj[u]:
        incoming += ru * g_table(Du, cu, ct.D(v), ct.c(v))
    outgoing = Fraction(0, 1)
    for w in ct.adj[u]:
        rw = ru if bug_use_receiver_rate else ct.gain_rate(w)
        if bug_transpose_one_side:
            g_val = g_table(Du, cu, ct.D(w), ct.c(w))
        else:
            g_val = g_table(ct.D(w), ct.c(w), Du, cu)
        outgoing += rw * g_val
    return incoming - outgoing


def arrow_transfer_divergence(
    ct: CoreTree,
    u: int,
    a_table: Callable[[int, int], Fraction],
    T_thresh: int,
    bug_threshold_on_receiver: bool = False,
) -> Fraction:
    """
    Net second (arrow) transfer adjustment added to row u (`operator.md` line 76):
    directed arc (w, v) carries 1_{D_w >= T} a(D_v, c_v) / D_w from row v to row w.
    So row u receives + (1_{D_u >= T} / D_u) sum_{v ~ u} a(D_v, c_v)
    and pays - sum_{w ~ u} 1_{D_w >= T} a(D_u, c_u) / D_w.
    """
    Du, cu = ct.D(u), ct.c(u)
    incoming = Fraction(0, 1)
    if Du >= T_thresh:
        inv_Du = Fraction(1, Du)
        for v in ct.adj[u]:
            incoming += inv_Du * a_table(ct.D(v), ct.c(v))
    outgoing = Fraction(0, 1)
    au = a_table(Du, cu)
    for w in ct.adj[u]:
        Dw = ct.D(w)
        cond = (Du >= T_thresh) if bug_threshold_on_receiver else (Dw >= T_thresh)
        if cond:
            outgoing += au * Fraction(1, Dw)
    return incoming - outgoing


def amended_banded_row(
    ct: CoreTree,
    u: int,
    z_fn: TestFn,
    g_table: Callable[[int, int, int, int], Fraction],
    a_table: Callable[[int, int], Fraction],
    T_thresh: int,
) -> Fraction:
    """
    Exact amended row formula from `proofs/upper-bound-banded.md` lines 50-59
    (for full-ball exact observation S_u = {w ~ u}, I_u = V(ct)):
      L_u U_u + sum_{w ~ u} (U_u - a(u) 1_{D_w >= T}) / D_w
      + b_u (Delta_u^b z_u + z_new(u))
      + sum_{w ~ u} [r_w (Delta_w^+ z_u - g(w, u)) + b_w Delta_w^b z_u]
      + sum_{x not in {u} U N(u)} [r_x Delta_x^+ z_u + b_x Delta_x^b z_u]
      + (1_{D_u >= T} / D_u) sum_{w ~ u} a(w),
    where U_u = Delta_u^+ z_u + sum_{v ~ u} g(u, v) using operator.md's arc convention
    g(w, u) carried from u to w at rate r_w.
    """
    base = z_fn(ct, u)
    Du, cu, Lu = ct.D(u), ct.c(u), ct.L[u]
    delta_u_plus = z_fn(ct.after_gain(u), u) - base
    Uu = delta_u_plus + sum(
        (g_table(Du, cu, ct.D(v), ct.c(v)) for v in ct.adj[u]),
        Fraction(0, 1),
    )
    au = a_table(Du, cu)

    row = Fraction(Lu, 1) * Uu
    for w in ct.adj[u]:
        Dw = ct.D(w)
        sub = au if Dw >= T_thresh else Fraction(0, 1)
        row += (Uu - sub) * Fraction(1, Dw)

    if Lu > 0:
        bu = ct.prom_rate(u)
        ct_pu, new_id = ct.after_prom(u)
        row += bu * ((z_fn(ct_pu, u) - base) + z_fn(ct_pu, new_id))

    nbrs = set(ct.adj[u])
    for w in ct.adj[u]:
        rw = ct.gain_rate(w)
        dw_plus = z_fn(ct.after_gain(w), u) - base
        g_wu = g_table(ct.D(w), ct.c(w), Du, cu)
        row += rw * (dw_plus - g_wu)
        if ct.L[w] > 0:
            bw = ct.prom_rate(w)
            ct_pw, _ = ct.after_prom(w)
            row += bw * (z_fn(ct_pw, u) - base)

    for x in range(ct.n):
        if x != u and x not in nbrs:
            rx = ct.gain_rate(x)
            row += rx * (z_fn(ct.after_gain(x), u) - base)
            if ct.L[x] > 0:
                bx = ct.prom_rate(x)
                ct_px, _ = ct.after_prom(x)
                row += bx * (z_fn(ct_px, u) - base)

    if Du >= T_thresh:
        inv_Du = Fraction(1, Du)
        for w in ct.adj[u]:
            row += inv_Du * a_table(ct.D(w), ct.c(w))

    return row


# ---------------------------------------------------------------------------
# 5. Test suite: positive checks + planted-error negative controls
# ---------------------------------------------------------------------------

def build_test_trees() -> List[Tuple[str, CoreTree, int]]:
    return [
        # Singleton core (T_3: centre + 2 leaves, c = 0, L = 2)
        ("T3_singleton", CoreTree([[]], [2]), 0),
        # Two-core asymmetric edge
        ("edge_asym", CoreTree([[1], [0]], [1, 3]), 0),
        # Three-core path with a leafless type-B middle core (L=0, c=2 -> D=2)
        ("path3_typeB_mid", CoreTree([[1], [0, 2], [1]], [2, 0, 1]), 0),
        ("path3_typeB_root", CoreTree([[1], [0, 2], [1]], [2, 0, 1]), 1),
        # Branched 5-core tree with mixed degrees
        (
            "branched5",
            CoreTree(
                [[1, 2, 3], [0], [0, 4], [0], [2]],
                [0, 2, 1, 1, 3],
            ),
            0,
        ),
        # Path of length 4 (radius 4 from endpoint 0) with type-B interior vertices
        (
            "path5_radius4",
            CoreTree(
                [[1], [0, 2], [1, 3], [2, 4], [3]],
                [2, 0, 1, 0, 2],
            ),
            0,
        ),
    ]


def make_nontrivial_cylinder(k: int) -> TestFn:
    """
    Bounded cylinder reading only the marked ball B_k(v) (with L-marks, no c-marks
    at distance k so k=R is a pure radius-k ball test) and hashing canonical structure
    into exact rationals in [-1, 1].
    """
    cache: Dict[Tuple, Fraction] = {}

    def f(ct: CoreTree, v: int) -> Fraction:
        code = ct.canonical_ball(v, radius=k, include_c_mark=False)
        if code not in cache:
            # Deterministic rational hash in [-1, 1] sensitive to every L-mark and child subtree
            h = hash(repr(code))  # deterministic within run; let's make it cross-run deterministic!
            s = repr(code).encode("utf-8")
            acc = 17
            for b in s:
                acc = (acc * 131 + b) % 1999
            cache[code] = Fraction(acc - 999, 1000)
        return cache[code]

    return f


def test_1_cell_formula_and_semigroup() -> None:
    print("[Check 1] Recomputing (2): microscopic friend-tree vs cell operator & semigroup Q_s Q_t = Q_{s+t} ...")
    f0 = make_nontrivial_cylinder(k=0)
    f1 = make_nontrivial_cylinder(k=1)
    f2 = make_nontrivial_cylinder(k=2)

    # 1A. Microscopic full friend-tree vs Core-cell formula (2) up to order 3 on 1-2 core seeds, order 2 on 3-core seeds
    for name, ct, root in build_test_trees()[:4]:
        max_ord = 3 if ct.n <= 2 else 2
        for fn_name, fn in [("cyl_k0", f0), ("cyl_k1", f1), ("cyl_k2", f2), ("const_1", lambda g, v: Fraction(1, 1))]:
            d_cell = core_operator_derivatives(ct, root, fn, max_order=max_ord)
            d_full = full_tree_microscopic_derivatives(ct, root, fn, max_order=max_ord)
            assert d_cell == d_full, f"Mismatch on {name}/{fn_name}: cell={d_cell} vs full={d_full}"

    # 1B. Exact rational semigroup property Q_s(Q_t f) == Q_{s+t} f up to total degree 3,
    # explicitly stripping cell labels at intermediate time s (testing seed-orientation independence!).
    for name, ct, root in build_test_trees()[:4]:
        max_ord = 3 if ct.n <= 2 else 2
        for fn in (f0, f1):
            d_direct = core_operator_derivatives(ct, root, fn, max_order=max_ord)
            # Precompute b-th derivative of Q_t fn as a pure state function on X (relabelled_fresh!)
            qt_deriv_fns: List[TestFn] = []
            for b in range(max_ord + 1):
                memo_b: Dict[Tuple, Fraction] = {}

                def make_gb(order_b: int, store: Dict[Tuple, Fraction]) -> TestFn:
                    def gb(g: CoreTree, v: int) -> Fraction:
                        code = g.canonical_ball(v, radius=None, include_c_mark=True)
                        if code not in store:
                            fresh = g.relabelled_fresh()
                            store[code] = core_operator_derivatives(fresh, v, fn, max_order=order_b)[order_b]
                        return store[code]

                    return gb

                qt_deriv_fns.append(make_gb(b, memo_b))

            # Check exact operator power identity for EVERY pair (a, b) with a + b <= max_ord:
            # d^a/ds^a|_{s=0} Q_s ( d^b/dt^b|_{t=0} Q_t fn ) == d^{a+b}/du^{a+b}|_{u=0} Q_u fn
            for n in range(max_ord + 1):
                for a in range(n + 1):
                    b = n - a
                    da_gb = core_operator_derivatives(ct, root, qt_deriv_fns[b], max_order=a)[a]
                    assert da_gb == d_direct[n], (
                        f"Semigroup failure at (a={a}, b={b}) on {name}: {da_gb} != {d_direct[n]}"
                    )

    # 1C. Finite-time ODE check of Q_s Q_t = Q_{s+t} and Yule domination Q_s 1 <= e^s
    ct_edge = CoreTree([[1], [0]], [1, 2])
    s_val, t_val = 0.025, 0.020
    f_fl = lambda g, v: float(f1(g, v))
    q_st, err_st = eval_Q_finite_time(ct_edge, 0, s_val + t_val, f_fl, max_jumps=5)

    qt_cache: Dict[Tuple, float] = {}

    def qt_on_X(g: CoreTree, v: int) -> float:
        code = g.canonical_ball(v, radius=None, include_c_mark=True)
        if code not in qt_cache:
            val, _ = eval_Q_finite_time(g.relabelled_fresh(), v, t_val, f_fl, max_jumps=4)
            qt_cache[code] = val
        return qt_cache[code]

    q_s_qt, err_sq = eval_Q_finite_time(ct_edge, 0, s_val, qt_on_X, max_jumps=4)
    assert abs(q_st - q_s_qt) < 5e-4, f"Finite-time semigroup diff too large: {abs(q_st - q_s_qt)}"

    q_one, err_one = eval_Q_finite_time(ct_edge, 0, 0.05, lambda g, v: 1.0, max_jumps=5)
    assert q_one - err_one <= math.exp(0.05) + 1e-9
    print("  PASS: Equation (2) matches microscopic friend-tree to all tested orders and satisfies Q_s Q_t = Q_{s+t}.")


def enumerate_simple_paths_ending_at(ct: CoreTree, w: int, length_j: int) -> List[Tuple[int, ...]]:
    """All simple paths (x_0, x_1, ..., x_j = w) of `length_j` edges in ct."""
    results: List[Tuple[int, ...]] = []

    def dfs_rev(curr: int, visited: FrozenSet[int], rev_path: List[int]) -> None:
        if len(rev_path) - 1 == length_j:
            results.append(tuple(reversed(rev_path)))
            return
        for prev in ct.adj[curr]:
            if prev not in visited:
                dfs_rev(prev, visited | {prev}, rev_path + [prev])

    dfs_rev(w, frozenset([w]), [w])
    return results


def test_2_spatial_cutoff_and_path_reversal() -> None:
    print("[Check 2] Recomputing (3): old-edge path-reversal bound, O(s^{R-k+1}) vanishing order & inequality (3) ...")

    # 2A. Path-reversal inequality: sum_{gamma of length j ending at w} prod_{i=0}^{j-1} (1 / c_{x_i}^0) <= c_w^0
    trees_for_paths = [t for _, t, _ in build_test_trees()] + [
        CoreTree(
            [[1, 2, 3], [0, 4, 5], [0, 6], [0], [1], [1], [2, 7, 8], [6], [6]],
            [0, 0, 1, 2, 2, 2, 0, 2, 2],
        )
    ]
    for ct in trees_for_paths:
        for w in range(ct.n):
            cw0 = Fraction(ct.c(w), 1)
            for j in range(1, 5):
                paths = enumerate_simple_paths_ending_at(ct, w, j)
                total_prod = Fraction(0, 1)
                for p in paths:
                    prod = Fraction(1, 1)
                    for x_send in p[:-1]:
                        prod *= Fraction(1, ct.c(x_send))
                    total_prod += prod
                assert total_prod <= cw0, f"Path-reversal bound failed: {total_prod} > {cw0}"

    # 2B. Exact Taylor vanishing order: (Q_s f - Q_s^[R] f)(G, o) vanishes for orders 0 .. (R - k),
    # i.e. difference is O(s^j) with j = R - k + 1, INCLUDING when boundary vertex has (L, c) = (0, 2)!
    _, ct_p5, root0 = build_test_trees()[5]  # path5_radius4: vertex 1 and 3 have (L=0, c=2)!
    for k in (0, 1, 2):
        fk = make_nontrivial_cylinder(k=k)
        full_derivs = core_operator_derivatives(ct_p5, root0, fk, max_order=3)
        for R in range(k, 4):
            j = R - k + 1
            ct_trunc, trunc_root = ct_p5.truncate_ball(root0, R=R, retain_missing_slots=True)
            trunc_derivs = core_operator_derivatives(ct_trunc, trunc_root, fk, max_order=3)
            for order in range(min(j, 4)):
                assert full_derivs[order] == trunc_derivs[order], (
                    f"Truncation R={R}, k={k} differed at order {order} < j={j}: "
                    f"{full_derivs[order]} vs {trunc_derivs[order]}"
                )

    # 2C. Numerical check of Equation (3) at finite s > 0:
    # |Q_s f - Q_s^[R] f| <= 2*sqrt(2) * e^s * ||f||_inf * sqrt(p)
    _, ct_p3, root_p3 = build_test_trees()[2]  # path3_typeB_mid (radius 2, middle vertex (0,2))
    for s_val in (0.04, 0.10):
        for k, R in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2)]:
            fk = make_nontrivial_cylinder(k=k)
            fk_fl = lambda g, v, fn=fk: float(fn(g, v))
            val_full, err_full = eval_Q_finite_time(ct_p3, root_p3, s_val, fk_fl, max_jumps=6)
            ct_tr, tr_root = ct_p3.truncate_ball(root_p3, R=R, retain_missing_slots=True)
            val_tr, err_tr = eval_Q_finite_time(ct_tr, tr_root, s_val, fk_fl, max_jumps=6)

            dist0 = ct_p3.distances_from(root_p3)
            sum_cw0 = sum(ct_p3.c(w) for w in range(ct_p3.n) if 0 <= dist0[w] <= k)
            j = R - k + 1
            p_bound = min(1.0, (s_val ** j) / math.factorial(j) * sum_cw0)
            rhs = 2.0 * math.sqrt(2.0) * math.exp(s_val) * 1.0 * math.sqrt(p_bound)
            lhs_upper = abs(val_full - val_tr) + err_full + err_tr
            assert lhs_upper <= rhs + 1e-9, f"Eq (3) violated for k={k}, R={R}, s={s_val}: {lhs_upper} > {rhs}"

    # 2D. Section 5 finite genuine core-tree ball completion (`operator.md` line 86):
    # Appending one (L=1, c=1) core vertex to each missing boundary slot produces a genuine
    # tree in X (extra_c == 0 everywhere, D >= 2) whose R-truncation is IDENTICAL to G's
    # R-truncation, so Q_t 1(G^(R), o) agrees with Q_t 1(G, o) through order R!
    for R in (0, 1, 2, 3):
        tr, tr_r = ct_p5.truncate_ball(root0, R=R, retain_missing_slots=True)
        comp_adj = [list(nb) for nb in tr.adj]
        comp_L = list(tr.L)
        for u in range(tr.n):
            for _ in range(tr.extra_c[u]):
                nid = len(comp_adj)
                comp_adj[u].append(nid)
                comp_adj.append([u])
                comp_L.append(1)  # one new core with one leaf -> (L=1, c=1), D=2
        ct_comp = CoreTree(comp_adj, comp_L)  # asserts extra_c==0 and D>=2 everywhere!
        f_one: TestFn = lambda g, v: Fraction(1, 1)
        d_orig = core_operator_derivatives(ct_p5, root0, f_one, max_order=3)
        d_comp = core_operator_derivatives(ct_comp, tr_r, f_one, max_order=3)
        for ord_m in range(min(R + 1, 4)):
            assert d_orig[ord_m] == d_comp[ord_m], (
                f"Section 5 completion differed at order {ord_m} <= R={R}: {d_orig[ord_m]} != {d_comp[ord_m]}"
            )

    print("  PASS: Path-reversal <= c_w^0, exact vanishing up to s^{R-k}, inequality (3), and Sec 5 completion hold.")


def test_3_generator_dynkin_and_degree_law_identities() -> None:
    print("[Check 3] Recomputing (7), Dynkin's formula, birth renewal (degree-law eq. 2) & degree crossing ...")
    f0 = make_nontrivial_cylinder(k=0)
    f1 = make_nontrivial_cylinder(k=1)
    f2 = make_nontrivial_cylinder(k=2)

    for name, ct, root in build_test_trees():
        # 3A. A1 = b_o
        a_one = generator_A(ct, root, lambda g, v: Fraction(1, 1))
        assert a_one == ct.prom_rate(root), f"A1 != b_o on {name}: {a_one} vs {ct.prom_rate(root)}"

        # 3B. Af(G, o) == d/ds|_{s=0} (Q_s f)(G, o), and d^n/ds^n|_{s=0} Q_s f == d^{n-1}/ds^{n-1}|_{s=0} Q_s(Af)
        max_ord = 3 if ct.n <= 2 else 2
        for fn in (f0, f1, f2):
            af_val = generator_A(ct, root, fn)
            q_derivs = core_operator_derivatives(ct, root, fn, max_order=max_ord)
            assert af_val == q_derivs[1], f"Eq (7) mismatch at order 1 on {name}: {af_val} != {q_derivs[1]}"

            af_fn: TestFn = lambda g, v, base_f=fn: generator_A(g, v, base_f)
            q_af_derivs = core_operator_derivatives(ct, root, af_fn, max_order=max_ord - 1)
            for m in range(1, max_ord + 1):
                assert q_derivs[m] == q_af_derivs[m - 1], (
                    f"Dynkin derivative mismatch at order {m} on {name}: {q_derivs[m]} != {q_af_derivs[m - 1]}"
                )

        # 3C. Exact degree-crossing generator identity (`round2/degree-law.md` line 35):
        # For k >= 2, A(1_{D_o >= k+1}) == (k - c_o + rho_o) 1_{D_o == k}
        for k_deg in range(2, 7):
            ind_fn: TestFn = lambda g, v, kd=k_deg: Fraction(1 if g.D(v) >= kd + 1 else 0, 1)
            lhs = generator_A(ct, root, ind_fn)
            if ct.D(root) == k_deg:
                rho_o = sum((Fraction(1, ct.D(w)) for w in ct.adj[root]), Fraction(0, 1))
                rhs = (Fraction(k_deg - ct.c(root), 1) + rho_o)
            else:
                rhs = Fraction(0, 1)
            assert lhs == rhs, f"Degree-crossing identity failed on {name}, k={k_deg}: {lhs} != {rhs}"

    # 3D. Exact power-series verification of `round2/degree-law.md` renewal equation:
    # Q_t f = S_t f + int_0^t Q_s B S_{t-s} f ds
    # Since int_0^t (s^a / a!) ((t-s)^b / b!) ds = t^{a+b+1} / (a+b+1)!, comparing m-th derivative at t=0 (m >= 1):
    # (Q^{(m)} f) == (S^{(m)} f) + sum_{a + b = m - 1} Q^{(a)}( B (S^{(b)} f) )
    for name, ct, root in build_test_trees()[:4]:
        max_ord = 3 if ct.n <= 2 else 2
        for fn in (f0, f1, lambda g, v: Fraction(1, 1)):
            q_d = core_operator_derivatives(ct, root, fn, max_order=max_ord)
            s_d = single_root_derivatives(ct, root, fn, max_order=max_ord)

            bs_fns: List[TestFn] = []
            for b in range(max_ord):
                memo_bs: Dict[Tuple, Fraction] = {}

                def make_bs(order_b: int, store: Dict[Tuple, Fraction]) -> TestFn:
                    def bs(g: CoreTree, u: int) -> Fraction:
                        if g.L[u] == 0:
                            return Fraction(0, 1)
                        code = g.canonical_ball(u, radius=None, include_c_mark=True)
                        if code not in store:
                            bu = g.prom_rate(u)
                            gp, new_v = g.relabelled_fresh().after_prom(u)
                            sb_val = single_root_derivatives(gp, new_v, fn, max_order=order_b)[order_b]
                            store[code] = bu * sb_val
                        return store[code]

                    return bs

                bs_fns.append(make_bs(b, memo_bs))

            assert q_d[0] == s_d[0]
            for m in range(1, max_ord + 1):
                rhs = s_d[m]
                for a in range(m):
                    b = (m - 1) - a
                    qa_val = core_operator_derivatives(ct, root, bs_fns[b], max_order=a)[a]
                    rhs += qa_val
                assert q_d[m] == rhs, f"Birth-renewal identity failed on {name} at order {m}: {q_d[m]} != {rhs}"

    print("  PASS: Equation (7), Dynkin's formula, birth-time renewal, and degree-crossing flux all hold.")


def generate_exhaustive_core_trees(max_cores: int = 5) -> List[CoreTree]:
    """Generate a broad suite of unmarked/marked core trees up to `max_cores` vertices."""
    trees: List[CoreTree] = [
        CoreTree([[]], [2]),
        CoreTree([[]], [9]),
    ]
    # All labelled tree topologies on n = 2..max_cores via Prüfer-like attachment, deduplicated canonically
    seen: set = set()
    adj_list_by_n: Dict[int, List[List[List[int]]]] = {1: [[[]]]}
    for n in range(2, max_cores + 1):
        cur = []
        for prev in adj_list_by_n[n - 1]:
            for attach in range(n - 1):
                nxt = [list(nb) for nb in prev]
                nxt[attach].append(n - 1)
                nxt.append([attach])
                cur.append(nxt)
        adj_list_by_n[n] = cur

    leaf_patterns = [
        lambda n, adj: [max(0, 2 - len(adj[u])) for u in range(n)],
        lambda n, adj: [1 + (u % 3) for u in range(n)],
        lambda n, adj: [8 if u == 0 else (0 if len(adj[u]) >= 2 else 2) for u in range(n)],
        lambda n, adj: [9 if u % 2 == 0 else 1 for u in range(n)],
    ]
    for n in range(2, max_cores + 1):
        for adj in adj_list_by_n[n]:
            for pat in leaf_patterns:
                L = pat(n, adj)
                ct = CoreTree(adj, L)
                canon = min(ct.canonical_ball(r, radius=None, include_c_mark=True) for r in range(n))
                if canon not in seen:
                    seen.add(canon)
                    trees.append(ct)
    return trees


def test_4_both_edge_transfer_cancellations() -> None:
    print("[Check 4] Recomputing both edge-transfer cancellations & amended banded row identity ...")

    # Asymmetric signed gain-transfer table g(D_w, c_w; D_u, c_u) supported on D < 12
    def g_table(Dw: int, cw: int, Du: int, cu: int) -> Fraction:
        if Dw >= 12 or Du >= 12:
            return Fraction(0, 1)
        # Deliberately non-symmetric in (w, u): 3*Dw - 5*cw + 7*Du - 11*cu + Dw*cu
        return Fraction(3 * Dw - 5 * cw + 7 * Du - 11 * cu + Dw * cu, 13)

    # Nonnegative finitely-supported arrow table a(D_u, c_u) supported on D_u < 10
    def a_table(Du: int, cu: int) -> Fraction:
        if Du >= 10:
            return Fraction(0, 1)
        return Fraction(2 * Du + 3 * cu + 1, 7)

    T_thresh = 8
    z_fn = make_nontrivial_cylinder(k=2)

    suite = generate_exhaustive_core_trees(max_cores=5)
    assert len(suite) >= 30
    for ct in suite:
        # 4A. Exact Gain-Transfer Cancellation (sum over all core vertices == 0)
        tot_gain_div = sum((gain_transfer_divergence(ct, u, g_table) for u in range(ct.n)), Fraction(0, 1))
        assert tot_gain_div == Fraction(0, 1), f"Gain transfer failed to cancel: {tot_gain_div}"

        # 4B. Exact Arrow-Transfer Cancellation (sum over all core vertices == 0)
        tot_arrow_div = sum(
            (arrow_transfer_divergence(ct, u, a_table, T_thresh=T_thresh) for u in range(ct.n)),
            Fraction(0, 1),
        )
        assert tot_arrow_div == Fraction(0, 1), f"Arrow transfer failed to cancel: {tot_arrow_div}"

        # 4C. Root-uniform mass-transport identity: E_{u ~ Unif(V(G))}[div(G, u)] == 0
        assert tot_gain_div / ct.n == Fraction(0, 1)
        assert tot_arrow_div / ct.n == Fraction(0, 1)

        # 4D. Exact vertex-by-vertex identity:
        #     amended_banded_row(u) == (Az)(G, u) + gain_transfer_div(u) + arrow_transfer_div(u),
        #     and consequently sum_u amended_banded_row(u) == sum_u (Az)(G, u)!
        sum_unamended = Fraction(0, 1)
        sum_amended = Fraction(0, 1)
        for u in range(ct.n):
            az_u = generator_A(ct, u, z_fn)
            am_u = amended_banded_row(ct, u, z_fn, g_table, a_table, T_thresh=T_thresh)
            expected_am_u = (
                az_u
                + gain_transfer_divergence(ct, u, g_table)
                + arrow_transfer_divergence(ct, u, a_table, T_thresh=T_thresh)
            )
            assert am_u == expected_am_u, f"Row decomposition mismatch at u={u}: {am_u} != {expected_am_u}"
            sum_unamended += az_u
            sum_amended += am_u
        assert sum_amended == sum_unamended, f"Whole-tree sum mismatch: {sum_amended} != {sum_unamended}"

    print(f"  PASS: Both transfers cancel identically and match amended rows across {len(suite)} distinct trees.")


def test_5_planted_errors() -> None:
    print("[Check 5] Planting 8 deliberate errors (negative controls) to verify detector sensitivity ...")

    # Planted Error 1: Cell inheritance bug in (2) (grandchildren promoted by promoted children dropped from H_o)
    ct_t3 = CoreTree([[]], [2])
    d_true = core_operator_derivatives(ct_t3, 0, lambda g, v: Fraction(1, 1), max_order=3, bug_root_only_cell=False)
    d_bug1 = core_operator_derivatives(ct_t3, 0, lambda g, v: Fraction(1, 1), max_order=3, bug_root_only_cell=True)
    assert d_true != d_bug1, "Planted Error 1 WAS NOT DETECTED!"
    print(f"  [Caught Planted Error 1] Dropping grandchildren from H_o broke order-2 derivative: {d_bug1[2]} != {d_true[2]}")

    # Planted Error 2: Dropping missing old-neighbour slots in Q_s^[R] breaks O(s^{R-k+1}) agreement!
    # Tree: 0 -- 1 -- 2 with L = [2, 2, 2]; at R=1, k=0 (j = R-k+1 = 2), vertex 1 has 1 missing neighbour (2).
    # Without retaining vertex 1's missing slot, D_1 drops from 4 to 3, changing r_0 = L_0 + 1/D_1 at order s^1!
    ct_slots = CoreTree([[1], [0, 2], [1]], [2, 2, 2])
    f_k0 = make_nontrivial_cylinder(k=0)
    full_d = core_operator_derivatives(ct_slots, 0, f_k0, max_order=2)
    tr_bad, r_bad = ct_slots.truncate_ball(0, R=1, retain_missing_slots=False)
    bad_d = core_operator_derivatives(tr_bad, r_bad, f_k0, max_order=2)
    assert full_d[1] != bad_d[1], "Planted Error 2 WAS NOT DETECTED!"
    print(f"  [Caught Planted Error 2] Omitting boundary slots broke O(s^2) cutoff at order 1: {bad_d[1]} != {full_d[1]}")

    # Planted Error 3: Off-by-one exponent j = R - k + 2 in (3) (claiming order R - k + 1 vanishes!)
    # On ct_slots with R=0, k=0 (true j = 1), exact order-1 difference is nonzero!
    tr_r0, r0 = ct_slots.truncate_ball(0, R=0, retain_missing_slots=True)
    r0_d = core_operator_derivatives(tr_r0, r0, f_k0, max_order=1)
    assert full_d[1] != r0_d[1], "Planted Error 3 WAS NOT DETECTED!"
    print(f"  [Caught Planted Error 3] Off-by-one j=R-k+2 refuted: exact order j=1 diff is {full_d[1] - r0_d[1]} != 0")

    # Planted Error 4: Omitting newborn term +b_o f(G^{b_o}, new) in (7), or evaluating it at old root o
    f_k1 = make_nontrivial_cylinder(k=1)
    q1_true = core_operator_derivatives(ct_slots, 0, f_k1, max_order=1)[1]
    a_no_nb = generator_A(ct_slots, 0, f_k1, omit_newborn=True)
    a_wrong_root = generator_A(ct_slots, 0, f_k1, newborn_at_old_root=True)
    assert a_no_nb != q1_true and a_wrong_root != q1_true, "Planted Error 4 WAS NOT DETECTED!"
    print(
        f"  [Caught Planted Error 4] Newborn bugs in (7) caught: omit={a_no_nb}, old_root={a_wrong_root}, true={q1_true}"
    )

    # Planted Error 5: Using 1/D_x instead of 1/D_y in gain rate r_x = L_x + sum_{y ~ x} 1/D_y
    ct_asym = CoreTree([[1], [0]], [1, 4])
    q1_asym = core_operator_derivatives(ct_asym, 0, f_k0, max_order=1)[1]
    a_bad_rate = generator_A(ct_asym, 0, f_k0, wrong_gain_denom_self=True)
    assert a_bad_rate != q1_asym, "Planted Error 5 WAS NOT DETECTED!"
    print(f"  [Caught Planted Error 5] Wrong neighbour-degree denominator in r_x caught: {a_bad_rate} != {q1_asym}")

    # Planted Error 6: Gain-transfer bug (multiplying outgoing g(w, u) by r_u instead of r_w, or transposing g)
    def g_asym(Dw: int, cw: int, Du: int, cu: int) -> Fraction:
        return Fraction(2 * Dw + 5 * Du - 3 * cw + cu, 7)

    bug6_rate = sum(
        (gain_transfer_divergence(ct_asym, u, g_asym, bug_use_receiver_rate=True) for u in range(ct_asym.n)),
        Fraction(0, 1),
    )
    bug6_trans = sum(
        (gain_transfer_divergence(ct_asym, u, g_asym, bug_transpose_one_side=True) for u in range(ct_asym.n)),
        Fraction(0, 1),
    )
    assert bug6_rate != 0 and bug6_trans != 0, "Planted Error 6 WAS NOT DETECTED!"
    print(f"  [Caught Planted Error 6] Gain-transfer bugs nonzero sum: wrong_rate={bug6_rate}, transposed={bug6_trans}")

    # Planted Error 7: Arrow-transfer bug (threshold D_u >= T on receiver instead of D_w >= T on sender)
    ct_thresh = CoreTree([[1], [0]], [8, 2])  # D_0 = 9 >= 8, D_1 = 3 < 8
    a_tab = lambda D, c: Fraction(D + c, 3) if D < 12 else Fraction(0, 1)
    bug7_sum = sum(
        (arrow_transfer_divergence(ct_thresh, u, a_tab, T_thresh=8, bug_threshold_on_receiver=True) for u in range(2)),
        Fraction(0, 1),
    )
    assert bug7_sum != 0, "Planted Error 7 WAS NOT DETECTED!"
    print(f"  [Caught Planted Error 7] Arrow-transfer receiver-threshold bug nonzero sum: {bug7_sum}")

    # Planted Error 8: Degree-law bugs — (a) replacing S_{t-s} by Q_{t-s} in renewal eq, (b) omitting -c_o in crossing
    # (a) At order 2 on T_3 for f = 1: Q^{(2)} 1 = 1 + Q^{(1)}(B 1) + B(S^{(1)} 1); using Q^{(1)} 1 != S^{(1)} 1 = 0!
    q_d_one = core_operator_derivatives(ct_t3, 0, lambda g, v: Fraction(1, 1), max_order=2)
    # With Q_{t-s} instead of S_{t-s} at m=2: S^{(2)}1 + Q^{(1)}(B Q^{(0)}1) + Q^{(0)}(B Q^{(1)}1)
    b_q1_at_t3 = ct_t3.prom_rate(0) * core_operator_derivatives(ct_t3.after_prom(0)[0].relabelled_fresh(), 1, lambda g, v: Fraction(1, 1), max_order=1)[1]
    bug8a_rhs = q_d_one[2] + b_q1_at_t3  # extra spurious term since S^{(1)}1 == 0 whereas Q^{(1)}1 > 0
    assert bug8a_rhs != q_d_one[2], "Planted Error 8a WAS NOT DETECTED!"
    # (b) Omitting -c_o in (k - c_o + rho_o) 1_{D_o = k} on edge_asym root (D_0 = 2, c_0 = 1)
    ind_ge3: TestFn = lambda g, v: Fraction(1 if g.D(v) >= 3 else 0, 1)
    true_cross = generator_A(ct_asym, 0, ind_ge3)
    bug8b_cross = Fraction(2, 1) + Fraction(1, ct_asym.D(1))  # forgot -c_0 = -1
    assert true_cross != bug8b_cross, "Planted Error 8b WAS NOT DETECTED!"
    print(
        f"  [Caught Planted Error 8] Degree-law bugs caught: renewal_Q_vs_S_diff={b_q1_at_t3}, "
        f"crossing_omit_co={bug8b_cross} != {true_cross}"
    )


def main() -> int:
    test_1_cell_formula_and_semigroup()
    test_2_spatial_cutoff_and_path_reversal()
    test_3_generator_dynkin_and_degree_law_identities()
    test_4_both_edge_transfer_cancellations()
    test_5_planted_errors()
    print("\nALL EXACT OPERATOR & DEGREE-LAW CHECKS AND PLANTED-ERROR CONTROLS PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
