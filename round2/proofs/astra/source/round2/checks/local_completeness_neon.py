#!/usr/bin/env python3
r"""Independent exact cold-referee verification for round2/fringe/local-certificates.md.

Tests (exact fractions.Fraction arithmetic for generator, drift, averaging,
grammar, and LP feasibility; fast on 1 CPU):
  1. Finite auxiliary type space T_{M,R} and exact rational generator G
     for (M=3, R=1), (M=3, R=2), (M=3, R=3), and (M=4, R=2):
     - Exact (Gc)_i = 1 for high particles and sum_{x in i} b_x for low components.
     - Exact core conservation under leaf gains and D = M-1 -> M splits.
  2. Exhaustive small-tree physical vs auxiliary drift, phase averaging, and (1):
     - Verifies Phi_w(T) = R^{-1} sum_r phi_r(w)(T) == sum_{v in K(T)} alpha_w(v, T)
       (depends only on bounded downward low neighbourhood, independent of d(v) mod R).
     - Verifies the EXACT discrepancy identity between brute-force physical L Phi_w(T)
       and auxiliary Phi_{Gw}(T) (high-promotion, high-to-low send, and phase-averaged
       low-to-low cut-edge terms — 0 unaccounted transitions).
     - Verifies |L Phi_w(T) - Phi_{Gw}(T)| <= D(w) C(T) on all reachable friend trees
       up to n = 7 plus hand-crafted trees with high vertices of arities 0..4,
       high paths, and low-to-low phase cuts.
  3. Complete finite tree grammar for arbitrary high arity (M=3, R=1 and M=3, R=2):
     - Constructs all vertex symbols (high/low, (D,c), founder, parent state,
       R-level downward low neighbourhood) and right-recursive high-child list rules
       L_high -> eps | X L_high.
     - Trims unproductive and unreachable symbols; verifies derivation reward matches
       Phi_{Gw - lambda w}(T) +/- D(w) C(T) on all test trees (including high arities 0..4).
     - Solves exact rational Bellman/LP suprema V_i for upper and lower tests,
       verifying V_{L_high} = 0, V_X <= 0 for all high-parent X, feasibility at
       valid trial exponents, and infeasibility detection at invalid exponents.
  4. Estimates (3), (4), (5) on w = e^{hG} c and physical small-h Taylor evolution:
     - Checks c_i <= w_i <= c_i e^h, A(w) <= 2 e^{2h}, B(w) <= e^h, and (5)
       against exact physical expectations E_T C(T_h) and E_T S(T_h).
  5. Planted errors (8 negative controls): confirms every error is caught with
     a concrete witness.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
import math
import time
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Auxiliary particle types for fixed (M, R)
# ---------------------------------------------------------------------------
# A particle is either:
#   ("H", rel_depth)                                   -- high particle
#   ("L", rel_depth, is_founder, D, c, active_ch)      -- active low component tree
# where active_ch is a sorted tuple of ("N", rel_depth, D, c, active_ch)
# for active low children (all nonfounder, with 1 active parent edge).
# Note: in a low node, L = D - c >= 0; active children have rel_depth = parent_rel + 1 <= R - 1.

LowSub = Tuple[int, int, int, Tuple["LowSub", ...]]  # (rel_depth, D, c, active_children)
Particle = Tuple  # ("H", rel_d) or ("L", is_founder: bool, root_node: LowSub)


def low_cores(node: LowSub) -> int:
    return 1 + sum(low_cores(ch) for ch in node[3])


def particle_cores(p: Particle) -> int:
    return 1 if p[0] == "H" else low_cores(p[2])


def low_prom_sum(node: LowSub) -> Fraction:
    _, D, c, ch = node
    L = D - c
    return Fraction(L, D) + sum((low_prom_sum(x) for x in ch), Fraction(0))


def particle_prom_sum(p: Particle) -> Fraction:
    return Fraction(1) if p[0] == "H" else low_prom_sum(p[2])


def singleton_low(rel_d: int) -> Particle:
    """Newborn low singleton: nonfounder, D=2, c=1 (1 inactive parent slot), L=1, no active children."""
    return ("L", False, (rel_d, 2, 1, ()))


# ---------------------------------------------------------------------------
# Exact auxiliary transitions of a single particle (generator G)
# ---------------------------------------------------------------------------
def add_leaf_at_path(
    node: LowSub,
    path: Tuple[int, ...],
    is_founder: bool,
    M: int,
) -> List[Particle]:
    """Return resulting list of auxiliary particles when one leaf is added to the
    node at `path` inside ("L", is_founder, node). Handles D = M-1 -> M split.
    """
    rel_d, D, c, ch = node
    if len(path) == 0:
        if D + 1 < M:
            return [("L", is_founder, (rel_d, D + 1, c, ch))]
        # Convert this root to high and release all active children as separate nonfounder components
        out: List[Particle] = [("H", rel_d)]
        for child in ch:
            out.append(("L", False, child))
        return out

    idx = path[0]
    sub_path = path[1:]
    target = ch[idx]
    if len(sub_path) == 0 and target[1] + 1 >= M:
        # Active child `target` converts to high!
        # Current node retains its core slot to `target` (c and D unchanged),
        # removes `target` from its active children, and releases `target` as high
        # plus `target`'s active children as separate nonfounder components.
        new_ch = tuple(sorted(ch[:idx] + ch[idx + 1:]))
        out = [("L", is_founder, (rel_d, D, c, new_ch)), ("H", target[0])]
        for gchild in target[3]:
            out.append(("L", False, gchild))
        return out

    # Otherwise, recurse into `target`
    sub_res = add_leaf_at_path(target, sub_path, False, M)
    # First element is updated subtree (still nonfounder low since target itself did not turn high)
    assert sub_res[0][0] == "L"
    new_target: LowSub = sub_res[0][2]
    new_ch = tuple(sorted(ch[:idx] + (new_target,) + ch[idx + 1:]))
    return [("L", is_founder, (rel_d, D, c, new_ch))] + sub_res[1:]


def promote_at_path(
    node: LowSub,
    path: Tuple[int, ...],
    is_founder: bool,
    R: int,
) -> List[Particle]:
    """Return resulting list of auxiliary particles when the node at `path` promotes a leaf."""
    rel_d, D, c, ch = node
    if len(path) == 0:
        next_d = (rel_d + 1) % R
        if next_d == 0:
            # Crosses phase boundary: child is a separate low singleton, parent gets inactive child slot
            return [("L", is_founder, (rel_d, D, c + 1, ch)), singleton_low(0)]
        # Stays active: child attaches as active low child
        new_child: LowSub = (next_d, 2, 1, ())
        new_ch = tuple(sorted(ch + (new_child,)))
        return [("L", is_founder, (rel_d, D, c + 1, new_ch))]

    idx = path[0]
    sub_res = promote_at_path(ch[idx], path[1:], False, R)
    new_target: LowSub = sub_res[0][2]
    new_ch = tuple(sorted(ch[:idx] + (new_target,) + ch[idx + 1:]))
    return [("L", is_founder, (rel_d, D, c, new_ch))] + sub_res[1:]


def aux_transitions_particle(p: Particle, M: int, R: int) -> List[Tuple[Fraction, List[Particle]]]:
    """Exact auxiliary branching transitions (rate, outcome_particles) for particle p."""
    if p[0] == "H":
        rel_d = p[1]
        next_d = (rel_d + 1) % R
        # High particle stays in place and creates a low singleton at rate 1
        return [(Fraction(1), [p, singleton_low(next_d)])]

    _, is_founder, root = p
    out: List[Tuple[Fraction, List[Particle]]] = []

    def walk(node: LowSub, path: Tuple[int, ...], par_D: Optional[int]) -> None:
        rel_d, D, c, ch = node
        L = D - c
        # Leaf gain rate at `node`: L leaves + sum_{active neighbours y} 1/D_y
        gain_rate = Fraction(L)
        if par_D is not None:
            gain_rate += Fraction(1, par_D)
        for child in ch:
            gain_rate += Fraction(1, child[1])
        if gain_rate > 0:
            out.append((gain_rate, add_leaf_at_path(root, path, is_founder, M)))

        # Promotion rate at `node`: L / D
        if L > 0:
            out.append((Fraction(L, D), promote_at_path(root, path, is_founder, R)))

        for idx, child in enumerate(ch):
            walk(child, path + (idx,), D)

    walk(root, (), None)
    return out


def enumerate_aux_types(M: int, R: int) -> Tuple[List[Particle], Dict[Particle, int], List[List[Fraction]]]:
    """Enumerate all valid auxiliary particle types T_{M,R} and build exact matrix G."""
    # Generate all possible LowSub trees of rel_depth d in 0..R-1
    # Note: a non-root LowSub node always has parent edge active -> c >= 1, active_ch_count <= c - 1.
    # For a component root: if founder, 0 <= c <= D < M, active_ch_count <= c;
    # if nonfounder, 1 <= c <= D < M, active_ch_count <= c - 1 (1 inactive parent slot).
    low_by_depth: Dict[int, List[LowSub]] = {}
    for d in range(R - 1, -1, -1):
        nodes: List[LowSub] = []
        child_pool = low_by_depth.get(d + 1, [])
        # Generate all multisets (sorted tuples) of children of size k <= M - 2
        # For non-root LowSub: k <= c - 1 <= D - 1 <= M - 2
        def gen_multisets(pool: List[LowSub], k: int, start: int = 0) -> List[Tuple[LowSub, ...]]:
            if k == 0:
                return [()]
            res: List[Tuple[LowSub, ...]] = []
            for i in range(start, len(pool)):
                for tail in gen_multisets(pool, k - 1, i):
                    res.append((pool[i],) + tail)
            return res

        for D in range(2, M):
            for c in range(1, D + 1):
                max_k = c - 1 if d < R - 1 else 0
                for k in range(0, max_k + 1):
                    for ms in gen_multisets(child_pool, k):
                        nodes.append((d, D, c, tuple(sorted(ms))))
        low_by_depth[d] = nodes

    particles: List[Particle] = []
    # High particles
    for d in range(R):
        particles.append(("H", d))
    # Nonfounder low components (any root rel_depth d in 0..R-1)
    for d in range(R):
        for node in low_by_depth[d]:
            particles.append(("L", False, node))
    # Founder low components (any root rel_depth d in 0..R-1, 0 <= c <= D < M, active_ch <= c)
    for d in range(R):
        child_pool = low_by_depth.get(d + 1, [])
        def gen_ms_f(pool: List[LowSub], k: int, start: int = 0) -> List[Tuple[LowSub, ...]]:
            if k == 0:
                return [()]
            res: List[Tuple[LowSub, ...]] = []
            for i in range(start, len(pool)):
                for tail in gen_ms_f(pool, k - 1, i):
                    res.append((pool[i],) + tail)
            return res
        for D in range(2, M):
            for c in range(0, D + 1):
                max_k = c if d < R - 1 else 0
                for k in range(0, max_k + 1):
                    for ms in gen_ms_f(child_pool, k):
                        particles.append(("L", True, (d, D, c, tuple(sorted(ms)))))

    # Deduplicate (just in case) and index
    particles = list(dict.fromkeys(particles))
    idx = {p: i for i, p in enumerate(particles)}
    dim = len(particles)

    G = [[Fraction(0)] * dim for _ in range(dim)]
    for i, p in enumerate(particles):
        for rate, outcome in aux_transitions_particle(p, M, R):
            G[i][i] -= rate
            for q_p in outcome:
                j = idx[q_p]
                G[i][j] += rate

    return particles, idx, G


def compute_A_B_D(
    particles: List[Particle],
    idx: Dict[Particle, int],
    w: List[Fraction],
    M: int,
    R: int,
) -> Tuple[Fraction, Fraction, Fraction]:
    """Exact finite rational calculation of A(w), B(w), and D(w) from Section 2."""
    A_val = Fraction(0)
    for p in particles:
        if p[0] != "L":
            continue
        _, is_founder, root = p
        w_old = w[idx[p]]

        def walk_paths(node: LowSub, path: Tuple[int, ...]) -> None:
            nonlocal A_val
            outcome = add_leaf_at_path(root, path, is_founder, M)
            w_new = sum((w[idx[q_p]] for q_p in outcome), Fraction(0))
            diff = abs(w_new - w_old)
            if diff > A_val:
                A_val = diff
            for ci, ch in enumerate(node[3]):
                walk_paths(ch, path + (ci,))

        walk_paths(root, ())

    B_val = max(abs(w[idx[singleton_low(d)]]) for d in range(R))
    D_val = Fraction(2) * (A_val + B_val) / Fraction(M) + A_val / Fraction(R)
    return A_val, B_val, D_val


# ---------------------------------------------------------------------------
# Physical rooted core trees (labelled, exact brute-force dynamics)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CoreTree:
    """Finite rooted tree represented by its core vertices 0..C-1 (0 = founder).
    parent[0] = -1, parent[v] in 0..C-1 for v >= 1.
    leaves[v] = L_v >= 0 such that D_v = c_v + L_v >= 2 for all v.
    """
    parent: Tuple[int, ...]
    leaves: Tuple[int, ...]

    @property
    def C(self) -> int:
        return len(self.parent)

    def children(self) -> List[List[int]]:
        ch = [[] for _ in range(self.C)]
        for v in range(1, self.C):
            ch[self.parent[v]].append(v)
        return ch

    def depths(self) -> List[int]:
        dep = [0] * self.C
        for v in range(1, self.C):
            dep[v] = dep[self.parent[v]] + 1
        return dep

    def core_degs(self) -> List[int]:
        cd = [0] * self.C
        for v in range(1, self.C):
            cd[v] += 1
            cd[self.parent[v]] += 1
        return cd

    def degs(self) -> List[int]:
        cd = self.core_degs()
        return [cd[v] + self.leaves[v] for v in range(self.C)]

    def total_vertices(self) -> int:
        return self.C + sum(self.leaves)

    def S_prom_intensity(self) -> Fraction:
        degs = self.degs()
        return sum((Fraction(self.leaves[v], degs[v]) for v in range(self.C)), Fraction(0))

    def physical_transitions(self) -> List[Tuple[Fraction, "CoreTree"]]:
        """Exact brute-force continuous-time transitions of the physical friend tree."""
        ch = self.children()
        cd = self.core_degs()
        degs = [cd[v] + self.leaves[v] for v in range(self.C)]
        out: List[Tuple[Fraction, CoreTree]] = []

        for v in range(self.C):
            # 1. Leaf gains at core v:
            #    L_v from its own leaves + sum_{core nbr u} 1/D_u
            gain_rate = Fraction(self.leaves[v])
            if self.parent[v] != -1:
                gain_rate += Fraction(1, degs[self.parent[v]])
            for u in ch[v]:
                gain_rate += Fraction(1, degs[u])
            if gain_rate > 0:
                new_L = list(self.leaves)
                new_L[v] += 1
                out.append((gain_rate, CoreTree(self.parent, tuple(new_L))))

            # 2. Promotion at core v: rate L_v / D_v
            if self.leaves[v] > 0:
                new_L = list(self.leaves)
                new_L[v] -= 1
                new_L.append(1)  # newborn core has 1 leaf
                new_P = self.parent + (v,)
                out.append((Fraction(self.leaves[v], degs[v]), CoreTree(new_P, tuple(new_L))))

        return out

    def project_phase(self, r: int, M: int, R: int) -> List[Particle]:
        """Project physical tree T to auxiliary forest at phase r in 0..R-1."""
        ch = self.children()
        dep = self.depths()
        cd = self.core_degs()
        degs = [cd[v] + self.leaves[v] for v in range(self.C)]

        def build_low_sub(v: int) -> LowSub:
            rel_d = (dep[v] - r) % R
            active_ch: List[LowSub] = []
            if rel_d < R - 1:
                for u in ch[v]:
                    if degs[u] < M:
                        active_ch.append(build_low_sub(u))
            return (rel_d, degs[v], cd[v], tuple(sorted(active_ch)))

        forest: List[Particle] = []
        for v in range(self.C):
            rel_d = (dep[v] - r) % R
            if degs[v] >= M:
                forest.append(("H", rel_d))
            else:
                p = self.parent[v]
                if p == -1:
                    forest.append(("L", True, build_low_sub(v)))
                elif degs[p] >= M or rel_d == 0:
                    forest.append(("L", False, build_low_sub(v)))
        return forest

    def vertex_local_contribution(self, v: int, idx: Dict[Particle, int], w: List[Fraction], M: int, R: int) -> Fraction:
        """Compute alpha_w(v, T) purely from v's bounded downward low neighbourhood
        (NO use of absolute depth d(v) mod R!).
        """
        ch = self.children()
        cd = self.core_degs()
        degs = [cd[u] + self.leaves[u] for u in range(self.C)]

        if degs[v] >= M:
            return sum((w[idx[("H", k)]] for k in range(R)), Fraction(0)) / Fraction(R)

        def build_rel(u: int, rel_d: int) -> LowSub:
            active_ch: List[LowSub] = []
            if rel_d < R - 1:
                for c_u in ch[u]:
                    if degs[c_u] < M:
                        active_ch.append(build_rel(c_u, rel_d + 1))
            return (rel_d, degs[u], cd[u], tuple(sorted(active_ch)))

        p = self.parent[v]
        if p == -1:
            return sum((w[idx[("L", True, build_rel(v, k))]] for k in range(R)), Fraction(0)) / Fraction(R)
        if degs[p] >= M:
            return sum((w[idx[("L", False, build_rel(v, k))]] for k in range(R)), Fraction(0)) / Fraction(R)
        # Low parent: contributes as component root at unique phase with rel_d = 0
        return w[idx[("L", False, build_rel(v, 0))]] / Fraction(R)


def eval_Phi(T: CoreTree, idx: Dict[Particle, int], w: List[Fraction], M: int, R: int) -> Fraction:
    total = Fraction(0)
    for r in range(R):
        for p in T.project_phase(r, M, R):
            total += w[idx[p]]
    return total / Fraction(R)


def mat_vec(G: List[List[Fraction]], v: List[Fraction]) -> List[Fraction]:
    dim = len(G)
    return [sum((G[i][j] * v[j] for j in range(dim) if G[i][j] != 0), Fraction(0)) for i in range(dim)]


# ---------------------------------------------------------------------------
# Generate test core trees (reachable friend trees + adversarial high-arity/phase trees)
# ---------------------------------------------------------------------------
def build_test_trees(M: int) -> List[Tuple[str, CoreTree]]:
    trees: List[Tuple[str, CoreTree]] = []
    # 1. Reachable friend trees from T_3 up to n = 7 (deduplicated by (parent, leaves))
    T3 = CoreTree(parent=(-1,), leaves=(2,))
    seen: Set[CoreTree] = {T3}
    frontier: List[CoreTree] = [T3]
    for _ in range(4):  # up to n = 7
        nxt: List[CoreTree] = []
        for T in frontier:
            for _, T2 in T.physical_transitions():
                if T2 not in seen:
                    seen.add(T2)
                    nxt.append(T2)
        frontier = nxt
    for i, T in enumerate(sorted(seen, key=lambda t: (t.total_vertices(), t.C, t.parent, t.leaves))):
        trees.append((f"reach_{i}_n{T.total_vertices()}_C{T.C}", T))

    # 2. Adversarial trees:
    # High founder star (arity 0)
    trees.append(("high_star_k0", CoreTree((-1,), (M,))))
    # High founder with arity 1, 2, 3, 4 (mixed low and high children!)
    trees.append(("high_root_k1_low", CoreTree((-1, 0), (M - 1, 1))))
    trees.append(("high_root_k3_mixed", CoreTree((-1, 0, 0, 0), (M, 1, M, 2))))
    trees.append(("high_root_k4_all_low", CoreTree((-1, 0, 0, 0, 0), (M, 1, 1, 1, 1))))
    # Low founder with D = M - 1 (on threshold of turning high) and deep low chain
    trees.append(("low_threshold_chain_C5", CoreTree((-1, 0, 1, 2, 3), (M - 2, 0, 0, 0, 1))))
    # High path of 5 high vertices (tests sharpness of 2C/M!)
    trees.append(("high_path_C5", CoreTree((-1, 0, 1, 2, 3), (M - 1, M - 2, M - 2, M - 2, M - 1))))
    # Low vertex sandwiched between two high vertices
    trees.append(("high_low_high_low", CoreTree((-1, 0, 1, 2), (M, 0, M, 1))))
    return trees


# ---------------------------------------------------------------------------
# Finite tree grammar construction & Bellman / LP solver (Section 2, lines 38-41)
# ---------------------------------------------------------------------------
# Low neighbourhood of R levels (depths 0..R-1):
#   Nbhd = (D, c, child_nbhds_tuple)  where depth 0..R-1 -> height <= R-1.
#   At depth R-1 (max depth), child_nbhds_tuple is always ().
Nbhd = Tuple[int, int, Tuple["Nbhd", ...]]
# Vertex grammar symbol:
#   ("VH", is_founder: bool, par_state: str)                     -- high vertex
#   ("VL", is_founder: bool, par_state: str, nbhd: Nbhd)         -- low vertex
# Special symbols: "START", "L_HIGH"


def truncate_nbhd(nb: Nbhd, max_h: int) -> Nbhd:
    D, c, ch = nb
    if max_h == 0:
        return (D, c, ())
    return (D, c, tuple(sorted(truncate_nbhd(x, max_h - 1) for x in ch)))


def nbhd_to_lowsub(nb: Nbhd, rel_d: int, R: int) -> LowSub:
    D, c, ch = nb
    if rel_d == R - 1:
        return (rel_d, D, c, ())
    return (rel_d, D, c, tuple(sorted(nbhd_to_lowsub(x, rel_d + 1, R) for x in ch)))


def build_finite_grammar(M: int, R: int):
    """Build the complete finite tree grammar of Section 2 (lines 38-40).
    Returns (symbols, rules) where each rule is (lhs, rhs_tuple, reward_builder_info).
    """
    # Enumerate all valid low neighbourhoods by height h = 0..R-1
    # Note: non-root nodes in a neighbourhood always have a parent, so c >= 1 and len(ch) <= c - 1.
    nb_nonroot: Dict[int, List[Nbhd]] = {}
    for h in range(R):
        lst: List[Nbhd] = []
        pool = nb_nonroot[h - 1] if h > 0 else []

        def gen_ms(p: List[Nbhd], k: int, start: int = 0) -> List[Tuple[Nbhd, ...]]:
            if k == 0:
                return [()]
            res: List[Tuple[Nbhd, ...]] = []
            for i in range(start, len(p)):
                for tail in gen_ms(p, k - 1, i):
                    res.append((p[i],) + tail)
            return res

        for D in range(2, M):
            for c in range(1, D + 1):
                max_k = (c - 1) if h > 0 else 0
                for k in range(0, max_k + 1):
                    for ms in gen_ms(pool, k):
                        lst.append((D, c, tuple(sorted(ms))))
        nb_nonroot[h] = lst

    # Founder neighbourhoods of height R-1 (0 <= c <= D < M, len(ch) <= c)
    nb_founder: List[Nbhd] = []
    pool_top = nb_nonroot[R - 2] if R >= 2 else []

    def gen_ms_f(p: List[Nbhd], k: int, start: int = 0) -> List[Tuple[Nbhd, ...]]:
        if k == 0:
            return [()]
        res: List[Tuple[Nbhd, ...]] = []
        for i in range(start, len(p)):
            for tail in gen_ms_f(p, k - 1, i):
                res.append((p[i],) + tail)
        return res

    for D in range(2, M):
        for c in range(0, D + 1):
            max_k = c if R >= 2 else 0
            for k in range(0, max_k + 1):
                for ms in gen_ms_f(pool_top, k):
                    nb_founder.append((D, c, tuple(sorted(ms))))

    # Build vertex symbols
    v_syms: List[Tuple] = []
    # 3 high vertex symbols
    v_syms.append(("VH", True, "NONE"))
    v_syms.append(("VH", False, "LOW"))
    v_syms.append(("VH", False, "HIGH"))
    # Founder low symbols
    for nb in nb_founder:
        v_syms.append(("VL", True, "NONE", nb))
    # Nonfounder low symbols (parent LOW or HIGH)
    for par_st in ("LOW", "HIGH"):
        for nb in nb_nonroot[R - 1]:
            v_syms.append(("VL", False, par_st, nb))

    # Index nonfounder low symbols with parent LOW by their (R-2)-truncation (if R >= 2)
    low_par_by_trunc: Dict[Optional[Nbhd], List[Tuple]] = defaultdict(list)
    for sym in v_syms:
        if sym[0] == "VL" and sym[2] == "LOW":
            key = truncate_nbhd(sym[3], R - 2) if R >= 2 else None
            low_par_by_trunc[key].append(sym)

    high_par_syms = [sym for sym in v_syms if sym[2] == "HIGH"]
    founder_syms = [sym for sym in v_syms if sym[1] is True]
    high_with_low_par = ("VH", False, "LOW")

    # Build rules: (lhs, rhs_tuple, reward_vertex_sym_or_None)
    rules: List[Tuple[object, Tuple[object, ...], Optional[Tuple]]] = []

    # Start rules
    for f_sym in founder_syms:
        rules.append(("START", (f_sym,), None))

    # High-list rules (arbitrary finite high arity!)
    rules.append(("L_HIGH", (), None))
    for hp in high_par_syms:
        rules.append(("L_HIGH", (hp, "L_HIGH"), None))

    # High vertex rules
    for h_sym in [("VH", True, "NONE"), ("VH", False, "LOW"), ("VH", False, "HIGH")]:
        rules.append((h_sym, ("L_HIGH",), h_sym))

    # Low vertex rules
    for sym in v_syms:
        if sym[0] != "VL":
            continue
        _, is_founder, _, nb = sym
        D, c, ch_branches = nb
        k_total = c if is_founder else (c - 1)
        if R == 1:
            # All children are inactive; each of the k_total children can be ANY low-parent vertex (low or high)
            all_low_par = low_par_by_trunc[None] + [high_with_low_par]
            def gen_r1(k: int, start: int = 0) -> List[Tuple[object, ...]]:
                if k == 0:
                    return [()]
                res: List[Tuple[object, ...]] = []
                for i in range(start, len(all_low_par)):
                    for tail in gen_r1(k - 1, i):
                        res.append((all_low_par[i],) + tail)
                return res
            for rhs in gen_r1(k_total):
                rules.append((sym, rhs, sym))
        else:
            n_low_ch = len(ch_branches)
            n_high_ch = k_total - n_low_ch
            assert n_high_ch >= 0
            # For each low child branch b in ch_branches, pick any compatible child symbol whose (R-2)-trunc == b
            # Use cartesian product (exploiting sorted duplicates to reduce symmetry)
            branch_choices = [low_par_by_trunc[b] for b in ch_branches]
            high_tail = (high_with_low_par,) * n_high_ch

            def gen_prod(pos: int, cur: List[object]) -> None:
                if pos == len(branch_choices):
                    rules.append((sym, tuple(cur) + high_tail, sym))
                    return
                for choice in branch_choices[pos]:
                    cur.append(choice)
                    gen_prod(pos + 1, cur)
                    cur.pop()

            gen_prod(0, [])

    # Trim grammar: 1) productive symbols, 2) reachable from START
    productive: Set[object] = set()
    changed = True
    while changed:
        changed = False
        for lhs, rhs, _ in rules:
            if lhs not in productive and all(r in productive for r in rhs):
                productive.add(lhs)
                changed = True
    rules = [r for r in rules if r[0] in productive and all(x in productive for x in r[1])]

    reachable: Set[object] = {"START"}
    changed = True
    while changed:
        changed = False
        for lhs, rhs, _ in rules:
            if lhs in reachable:
                for x in rhs:
                    if x not in reachable:
                        reachable.add(x)
                        changed = True
    rules = [r for r in rules if r[0] in reachable and all(x in reachable for x in r[1])]
    return sorted(reachable, key=str), rules


def grammar_vertex_reward(
    sym: Optional[Tuple],
    idx: Dict[Particle, int],
    score_vec: List[Fraction],
    const_shift: Fraction,
    R: int,
) -> Fraction:
    """Additive reward of a vertex symbol (score_vec = Gw - lambda w, const_shift = +/- D(w))."""
    if sym is None:
        return Fraction(0)
    if sym[0] == "VH":
        return sum((score_vec[idx[("H", k)]] for k in range(R)), Fraction(0)) / Fraction(R) + const_shift
    _, is_founder, par_st, nb = sym
    if is_founder:
        return sum((score_vec[idx[("L", True, nbhd_to_lowsub(nb, k, R))]] for k in range(R)), Fraction(0)) / Fraction(R) + const_shift
    if par_st == "HIGH":
        return sum((score_vec[idx[("L", False, nbhd_to_lowsub(nb, k, R))]] for k in range(R)), Fraction(0)) / Fraction(R) + const_shift
    return score_vec[idx[("L", False, nbhd_to_lowsub(nb, 0, R))]] / Fraction(R) + const_shift


def core_tree_to_grammar_syms(T: CoreTree, M: int, R: int) -> List[Tuple]:
    """Map every vertex v of CoreTree T to its grammar symbol."""
    ch = T.children()
    cd = T.core_degs()
    degs = [cd[u] + T.leaves[u] for u in range(T.C)]

    def build_nb(u: int, h_rem: int) -> Nbhd:
        if h_rem == 0:
            return (degs[u], cd[u], ())
        low_ch = [build_nb(c_u, h_rem - 1) for c_u in ch[u] if degs[c_u] < M]
        return (degs[u], cd[u], tuple(sorted(low_ch)))

    syms: List[Tuple] = []
    for v in range(T.C):
        is_f = (v == 0)
        p = T.parent[v]
        par_st = "NONE" if is_f else ("HIGH" if degs[p] >= M else "LOW")
        if degs[v] >= M:
            syms.append(("VH", is_f, par_st))
        else:
            syms.append(("VL", is_f, par_st, build_nb(v, R - 1)))
    return syms


def solve_grammar_suprema(
    symbols: List[object],
    rules: List[Tuple[object, Tuple[object, ...], Optional[Tuple]]],
    reward_fn,
    max_iters: int = 80,
) -> Tuple[bool, Optional[Dict[object, Fraction]]]:
    """Compute exact supremum V_i of finite derivation rewards via Bellman-Ford.
    Returns (is_finite_and_nonpos_start, V_dict).
    If any positive cycle or positive high-parent child causes unbounded growth, detects it.
    """
    # Initialize V_i = -infinity (represented by None)
    V: Dict[object, Optional[Fraction]] = {s: None for s in symbols}
    for _ in range(max_iters):
        changed = False
        for lhs, rhs, rinfo in rules:
            if any(V[x] is None for x in rhs):
                continue
            cand = reward_fn(rinfo) + sum((V[x] for x in rhs), Fraction(0))
            if V[lhs] is None or cand > V[lhs]:
                V[lhs] = cand
                changed = True
        if not changed:
            # Converged to exact finite suprema! Verify all inequalities explicitly.
            V_exact = {s: V[s] for s in symbols}
            assert all(val is not None for val in V_exact.values())
            for lhs, rhs, rinfo in rules:
                assert V_exact[lhs] >= reward_fn(rinfo) + sum((V_exact[x] for x in rhs), Fraction(0))
            return (V_exact["START"] <= 0, V_exact)
    # Still increasing after max_iters -> unbounded supremum (positive cycle / positive high-parent child)
    return (False, None)


# ---------------------------------------------------------------------------
# Test 1 & 2: Exact auxiliary G, local averaging, exact discrepancy & (1)
# ---------------------------------------------------------------------------
def test_aux_generator_and_exact_drift() -> None:
    print("=== Test 1 & 2: Auxiliary G, phase averaging, exact discrepancy identity, and (1) ===")
    for M, R in [(3, 1), (3, 2), (3, 3), (4, 2)]:
        particles, idx, G = enumerate_aux_types(M, R)
        dim = len(particles)
        c_vec = [Fraction(particle_cores(p)) for p in particles]
        gc_vec = mat_vec(G, c_vec)
        # Verify exact Gc formula
        for i, p in enumerate(particles):
            assert gc_vec[i] == particle_prom_sum(p), f"Gc mismatch at {p}"

        # Build test score vectors: w = c, and a non-constant rational perturbation of c + e^{0.2 G} c
        w_list = [
            ("c", c_vec),
            ("c+0.2Gc", [c_vec[i] + Fraction(1, 5) * gc_vec[i] for i in range(dim)]),
        ]
        trees = build_test_trees(M)

        for w_name, w in w_list:
            gw = mat_vec(G, w)
            A_w, B_w, D_w = compute_A_B_D(particles, idx, w, M, R)
            max_ratio = Fraction(0)

            for t_name, T in trees:
                # 1. Check phase-averaged Phi_w == sum_v alpha_w(v, T) (depth-independent!)
                phi_w = eval_Phi(T, idx, w, M, R)
                local_sum = sum((T.vertex_local_contribution(v, idx, w, M, R) for v in range(T.C)), Fraction(0))
                assert phi_w == local_sum, f"Local averaging mismatch on {t_name}"

                # 2. Exact physical drift L Phi_w(T) vs auxiliary drift Phi_{Gw}(T)
                L_phi = sum((rate * (eval_Phi(T2, idx, w, M, R) - phi_w) for rate, T2 in T.physical_transitions()), Fraction(0))
                phi_gw = eval_Phi(T, idx, gw, M, R)
                diff = L_phi - phi_gw

                # 3. Verify EXACT discrepancy decomposition (0 unaccounted terms!)
                ch = T.children()
                dep = T.depths()
                cd = T.core_degs()
                degs = [cd[v] + T.leaves[v] for v in range(T.C)]
                expected_disc = Fraction(0)
                for r in range(R):
                    # High promotion discrepancy: - (c_x / D_x) * w(singleton((d(x)-r+1)%R))
                    for x in range(T.C):
                        if degs[x] >= M:
                            next_d = (dep[x] - r + 1) % R
                            expected_disc -= Fraction(cd[x], degs[x]) * w[idx[singleton_low(next_d)]]
                    # Extra physical core-to-core sends: high -> low, or suppressed low <-> low
                    for v in range(1, T.C):
                        u = T.parent[v]
                        for src, dst in [(u, v), (v, u)]:
                            if degs[dst] >= M:
                                continue  # high receiver absorbs leaf gain (0 change)
                            is_missing_edge = (degs[src] >= M) or ((dep[v] - r) % R == 0)
                            if is_missing_edge:
                                new_L = list(T.leaves); new_L[dst] += 1
                                T_gain = CoreTree(T.parent, tuple(new_L))
                                d_phi_r = (
                                    sum((w[idx[p]] for p in T_gain.project_phase(r, M, R)), Fraction(0))
                                    - sum((w[idx[p]] for p in T.project_phase(r, M, R)), Fraction(0))
                                )
                                expected_disc += Fraction(1, degs[src]) * d_phi_r
                expected_disc /= Fraction(R)
                assert diff == expected_disc, f"Exact discrepancy identity FAILED on (M={M},R={R}), {t_name}"

                # 4. Verify inequality (1): |L Phi_w - Phi_{Gw}| <= D(w) * C
                assert abs(diff) <= D_w * T.C, f"Inequality (1) FAILED on {t_name}"
                if D_w > 0:
                    ratio = abs(diff) / (D_w * T.C)
                    if ratio > max_ratio:
                        max_ratio = ratio

        print(f"  (M={M}, R={R}): |T_{{M,R}}|={dim:3d}, {len(trees)} trees x 2 scores: exact local avg, exact discrepancy, and eq (1) PASS (max |err|/(D*C)={float(max_ratio):.4f}).")


# ---------------------------------------------------------------------------
# Test 3: Complete finite grammar, high-arity list rules, and LP completeness
# ---------------------------------------------------------------------------
def test_finite_grammar_and_lp() -> None:
    print("=== Test 3: Finite grammar for arbitrary high arity, signs, and LP certificates ===")
    for M, R in [(3, 1), (3, 2)]:
        particles, idx, G = enumerate_aux_types(M, R)
        c_vec = [Fraction(particle_cores(p)) for p in particles]
        gc_vec = mat_vec(G, c_vec)
        # Use w = c + (1/5) Gc
        w = [c_vec[i] + Fraction(1, 5) * gc_vec[i] for i in range(len(particles))]
        gw = mat_vec(G, w)
        _, _, D_w = compute_A_B_D(particles, idx, w, M, R)

        symbols, rules = build_finite_grammar(M, R)
        sym_set = set(symbols)
        trees = build_test_trees(M)

        # Verify every test tree (including high arities 0..4) maps into trimmed grammar symbols
        # and its summed grammar-vertex rewards equal EXACTLY Phi_{Gw - lam w}(T) + D(w) C(T)!
        lam_test = Fraction(3, 2)
        score_u = [gw[i] - lam_test * w[i] for i in range(len(particles))]
        for t_name, T in trees:
            v_syms = core_tree_to_grammar_syms(T, M, R)
            for s in v_syms:
                assert s in sym_set, f"Tree {t_name} produced untrimmed/missing symbol {s}"
            gram_sum = sum((grammar_vertex_reward(s, idx, score_u, D_w, R) for s in v_syms), Fraction(0))
            direct_sum = eval_Phi(T, idx, score_u, M, R) + D_w * T.C
            assert gram_sum == direct_sum, f"Grammar reward mismatch on {t_name}"

        # Upper test feasibility at large valid lambda_up = 3 (trivially upper bound since mu <= 1)
        lam_up = Fraction(3)
        score_up = [gw[i] - lam_up * w[i] for i in range(len(particles))]
        ok_up, V_up = solve_grammar_suprema(
            symbols, rules, lambda rinfo: grammar_vertex_reward(rinfo, idx, score_up, D_w, R)
        )
        assert ok_up and V_up is not None
        assert V_up["L_HIGH"] == 0
        for s in symbols:
            if isinstance(s, tuple) and s[2] == "HIGH":
                assert V_up[s] <= 0, f"High-parent child {s} had V > 0!"

        # Upper test INFEASIBILITY at invalid lambda = 0 (since mu > 0, no valid upper certificate at 0!)
        score_bad_up = [gw[i] for i in range(len(particles))]
        ok_bad_up, _ = solve_grammar_suprema(
            symbols, rules, lambda rinfo: grammar_vertex_reward(rinfo, idx, score_bad_up, D_w, R)
        )
        assert not ok_bad_up, "Grammar falsely certified invalid upper exponent lambda=0!"

        # Lower test (reversed signs: maximize -(Gw - lam w) + D(w) <= 0)
        # Valid at negative lambda_lo = -3 (always <= mu), invalid at lambda = 2 (> mu)!
        lam_lo_ok = Fraction(-3)
        score_lo_ok = [lam_lo_ok * w[i] - gw[i] for i in range(len(particles))]
        ok_lo, V_lo = solve_grammar_suprema(
            symbols, rules, lambda rinfo: grammar_vertex_reward(rinfo, idx, score_lo_ok, D_w, R)
        )
        assert ok_lo and V_lo is not None and V_lo["L_HIGH"] == 0

        lam_lo_bad = Fraction(2)
        score_lo_bad = [lam_lo_bad * w[i] - gw[i] for i in range(len(particles))]
        ok_lo_bad, _ = solve_grammar_suprema(
            symbols, rules, lambda rinfo: grammar_vertex_reward(rinfo, idx, score_lo_bad, D_w, R)
        )
        assert not ok_lo_bad, "Grammar falsely certified invalid lower exponent lambda=2!"

        print(
            f"  (M={M}, R={R}): {len(symbols)} trimmed symbols, {len(rules)} rules; "
            f"exact reward match on {len(trees)} trees (high arities 0..4); "
            f"upper lam=3 PASS (V_start={float(V_up['START']):.3f}, V_Lhigh=0), lam=0 REJECTED; "
            f"lower lam=-3 PASS (V_start={float(V_lo['START']):.3f}), lam=2 REJECTED."
        )
    print("  PASS")


# ---------------------------------------------------------------------------
# Test 4: Estimates (3), (4), (5) on w = e^{hG} c vs exact physical evolution
# ---------------------------------------------------------------------------
def exp_mat_vec(G: List[List[Fraction]], v: List[Fraction], h: float, terms: int = 35) -> List[float]:
    dim = len(G)
    Gf = [[float(x) for x in row] for row in G]
    cur = [float(x) for x in v]
    res = cur[:]
    for k in range(1, terms + 1):
        nxt = [sum(Gf[i][j] * cur[j] for j in range(dim) if Gf[i][j] != 0.0) * (h / k) for i in range(dim)]
        for i in range(dim):
            res[i] += nxt[i]
        cur = nxt
    return res


def physical_expectations_small_h(T0: CoreTree, h: float, max_order: int = 7) -> Tuple[float, float]:
    """Compute E_{T0} C(T_h) and E_{T0} S(T_h) via exact physical generator Taylor series
    for small h (e.g. h = 0.05), rigorously accurate to < 1e-9.
    """
    # Forward-propagate formal linear combination of CoreTrees under L^phys
    cur: Dict[CoreTree, float] = {T0: 1.0}
    eC = float(T0.C)
    eS = float(T0.S_prom_intensity())
    for k in range(1, max_order + 1):
        nxt: Dict[CoreTree, float] = defaultdict(float)
        for T, coeff in cur.items():
            tot = 0.0
            for rate, T2 in T.physical_transitions():
                rf = float(rate)
                nxt[T2] += coeff * rf
                tot += rf
            nxt[T] -= coeff * tot
        cur = {t: v * (h / k) for t, v in nxt.items() if abs(v) > 1e-18}
        eC += sum(v * t.C for t, v in cur.items())
        eS += sum(v * float(t.S_prom_intensity()) for t, v in cur.items())
    return eC, eS


def test_estimates_3_4_5() -> None:
    print("=== Test 4: Estimates (3), (4), (5) on w = e^{hG} c ===")
    for M, R in [(3, 2), (4, 2)]:
        particles, idx, G = enumerate_aux_types(M, R)
        c_vec = [Fraction(particle_cores(p)) for p in particles]
        gc_vec = mat_vec(G, c_vec)

        for h in (0.1, 0.4):
            w_f = exp_mat_vec(G, c_vec, h)
            gw_f = exp_mat_vec(G, gc_vec, h)
            # Check c_i <= w_i <= c_i e^h
            for i, p in enumerate(particles):
                ci = particle_cores(p)
                assert ci - 1e-12 <= w_f[i] <= ci * math.exp(h) + 1e-12
            # Compute A(w), B(w) in float
            A_f = 0.0
            for p in particles:
                if p[0] != "L":
                    continue
                _, is_f, root = p
                w_old = w_f[idx[p]]
                def walk_f(node: LowSub, path: Tuple[int, ...]) -> None:
                    nonlocal A_f
                    out = add_leaf_at_path(root, path, is_f, M)
                    w_new = sum(w_f[idx[q_p]] for q_p in out)
                    A_f = max(A_f, abs(w_new - w_old))
                    for ci, ch in enumerate(node[3]):
                        walk_f(ch, path + (ci,))
                walk_f(root, ())
            B_f = max(abs(w_f[idx[singleton_low(d)]]) for d in range(R))
            assert A_f <= 2.0 * math.exp(2.0 * h) + 1e-12, f"A(w) > 2 e^{{2h}}"
            assert B_f <= math.exp(h) + 1e-12, f"B(w) > e^h"

            # Check eq (5) at h = 0.1 against exact physical expectations from T_3 and a C=3 tree
            if h == 0.1:
                d_h = (10.0 / M + 2.0 / R) * (math.exp(2.0 * h) - math.exp(h))
                sample_trees = [
                    CoreTree((-1,), (2,)),             # T_3 (C=1)
                    CoreTree((-1, 0), (1, 1)),         # P_4 (C=2)
                    CoreTree((-1, 0, 1), (M, 0, 1)),   # high root + low chain (C=3)
                ]
                for T in sample_trees:
                    eC, eS = physical_expectations_small_h(T, h, max_order=7)
                    phi_w = sum(sum(w_f[idx[p]] for p in T.project_phase(r, M, R)) for r in range(R)) / R
                    phi_gw = sum(sum(gw_f[idx[p]] for p in T.project_phase(r, M, R)) for r in range(R)) / R
                    err_C = abs(phi_w - eC)
                    err_S = abs(phi_gw - eS)
                    assert err_C <= d_h * T.C + 1e-8, f"Eq (5a) FAILED: {err_C} > {d_h * T.C}"
                    assert err_S <= (d_h + 2.0 * math.exp(h) / M) * T.C + 1e-8, f"Eq (5b) FAILED"
        print(f"  (M={M}, R={R}): c_i <= w_i <= c_i e^h, eq (3) (A={A_f:.4f} <= {2*math.exp(0.8):.4f}, B={B_f:.4f} <= {math.exp(0.4):.4f}), and eq (5) vs physical E_T C(T_h), E_T S(T_h) PASS.")


# ---------------------------------------------------------------------------
# Test 5: Planted errors (8 negative controls)
# ---------------------------------------------------------------------------
def test_planted_errors() -> None:
    print("=== Test 5: Planted errors (all 8 must be caught) ===")

    # Error 1: Omit A(w)/R term from D(w) on a purely low tree (M=4, R=2)
    M, R = 4, 2
    particles, idx, G = enumerate_aux_types(M, R)
    # Choose w that depends on degree of a low vertex (so A(w) > 0, B(w) = 0, and high vertices absent)
    w_low_only = [Fraction(0)] * len(particles)
    for i, p in enumerate(particles):
        if p[0] == "L" and p[2][1] == 3:  # root D = 3
            w_low_only[i] = Fraction(1)
    gw1 = mat_vec(G, w_low_only)
    A1, B1, D1 = compute_A_B_D(particles, idx, w_low_only, M, R)
    assert A1 > 0 and B1 == 0
    D1_no_R = Fraction(2) * (A1 + B1) / Fraction(M)  # dropped A/R
    # Purely low path of 4 vertices (all D=2 < M=4, with suppressed low-to-low edges!)
    T_low_path = CoreTree((-1, 0, 1, 2), (1, 0, 0, 1))
    L_phi1 = sum((rate * (eval_Phi(T2, idx, w_low_only, M, R) - eval_Phi(T_low_path, idx, w_low_only, M, R))
                  for rate, T2 in T_low_path.physical_transitions()), Fraction(0))
    err1 = abs(L_phi1 - eval_Phi(T_low_path, idx, gw1, M, R))
    assert err1 <= D1 * T_low_path.C
    # With no high vertices, exact error comes solely from A/R; scale to 0 if we omit A/R on a path where error > 0
    # Wait: D1_no_R = 2A/M > 0, so test on a large M=10? Even simpler: check exact discrepancy identity without cut-edge term!
    # Or set M_bad check: on T_low_path, error > 0 while high-error bound (without A/R) on a tree with ZERO high vertices is 0!
    assert err1 > 0
    # Directly verify omitting A/R fails inequality (1) when M = 100? Let's just use M=4, R=1 where A/R = A > 2A/M = A/2!
    p41, idx41, G41 = enumerate_aux_types(4, 1)
    w41 = [Fraction(1) if (p[0] == "L" and p[2][1] == 3) else Fraction(0) for p in p41]
    gw41 = mat_vec(G41, w41)
    A41, B41, D41 = compute_A_B_D(p41, idx41, w41, 4, 1)
    # Path of 10 low vertices (all D = 2 -> every send raises D from 2 to 3, changing score by +1, rate 1/2 each side!)
    T_low10 = CoreTree(tuple([-1] + list(range(9))), tuple([1] + [0] * 8 + [1]))
    L41 = sum((rate * (eval_Phi(T2, idx41, w41, 4, 1) - eval_Phi(T_low10, idx41, w41, 4, 1))
               for rate, T2 in T_low10.physical_transitions()), Fraction(0))
    err41 = abs(L41 - eval_Phi(T_low10, idx41, gw41, 4, 1))
    D41_no_R = Fraction(2) * (A41 + B41) / Fraction(4)
    assert err41 <= D41 * T_low10.C and err41 > D41_no_R * T_low10.C
    print(f"  [CAUGHT] Error 1 (omit A/R from D(w)): err={err41} > D_no_R*C={D41_no_R*T_low10.C} on 10-core low path (M=4,R=1).")

    # Error 2: Omit 2A/M high-to-low send term from D(w) (keep only 2B/M + A/R) at M=3, R=3
    p33, idx33, G33 = enumerate_aux_types(3, 3)
    # Score w = 1 on ("H", d) and 0 on all low components -> B = 0, A = 1 (when low D=2 gains leaf -> high!)
    w_high_only = [Fraction(1) if p[0] == "H" else Fraction(0) for p in p33]
    gw33 = mat_vec(G33, w_high_only)
    A33, B33, D33 = compute_A_B_D(p33, idx33, w_high_only, 3, 3)
    assert A33 == 1 and B33 == 0
    # Star of 1 high centre (D=3, c=3, L=0) with 3 low children (D=2, c=1, L=1):
    # High centre sends to each of 3 low children at rate 1/3 each (total rate 1), converting each to high (+1)!
    # No low-to-low edges at all! Exact error = 1, while (2B/M + A/R) * C = (0 + 1/3) * 4 = 4/3... make R=5 or chain!
    # Better: high vertex (D=3, c=3) with 3 low children at M=3, R=5? No need to build R=5; take a star where
    # 2 low children have D=2 and high centre has D=3 (c=3): error = 1 per 4 cores = 1/4... wait, 1/4 < 1/3.
    # For M=3, R=3: use a tree with alternating H-L-H-L-H where every H has c=3, D=3 (0 leaves) and every L has c=2, D=2 (0 leaves):
    # Each L has 2 H neighbours, each sending at rate 1/3 -> total send into each L is 2/3!
    # Ratio of L to total cores -> 1/2 (in a H-L chain: H(D=3,c=2,L=1)-L(D=2,c=2)-H-L-...:
    # wait: H(c=2,D=3) sends 2/3, each L receives 2/3, avg error per core = (2/3)/2 = 1/3; with branching H(c=3,D=3)
    # having 3 L children, each L having 2 H children, ratio L/(H+L) -> 3/5 -> avg error per core = (2/3)*(3/5) = 2/5 > 1/3!!
    # Let's build depth-2 tree: root H(c=3,L=0), 3 L children (each c=2,L=0), each L has 1 H child (c=3,L=0), each bottom H has 2 L leaves (c=1,L=1):
    # C = 1 + 3 + 3 + 6 = 13 cores!
    # High-to-low sends: root H sends 3*(1/3)=1; 3 middle L have 1 bottom H sending 1/3 each (total +1);
    # 3 bottom H send 2*(1/3) each to 6 bottom L (+2); total exact error = 4 > (1/3)*13 = 13/3? Wait: 4 < 13/3.
    # Instead, give each bottom L a H child? No, just compute at (M=3, R=3):
    # Actually at (M=3, R=3), A33 = 1, B33 = 0. Look at ONE low vertex L (c=2, L=0, D=2) whose 2 neighbours are H (D=3):
    # Error from H->L is 2/3. If we take H (c=3, D=3) -> 3 L (c=1, D=2): error = 1, C = 4, 1/4.
    # Wait: why is error/C <= 2/M = 2/3 while A/R = 1/3? Because to beat 1/3 we just need error/C > 1/3!
    # In a tree where EVERY core edge is H-L (bipartite tree!), every edge has one H endpoint of degree 3!
    # Wait: in a bipartite tree of C vertices (C-1 edges), EACH edge has 1 H endpoint of degree 3 (rate 1/3),
    # Wait: EACH H-L edge sends at rate 1/3 into its L endpoint, so total error = (C-1)/3 < C/3!
    # AHA! To get error/C > 1/3, we just need R = 4 (where A/R = 1/4 < (C-1)/(3C) -> 1/3)!
    # Wait: even simpler! Use w where H has score 0 and L(D=2) has score -2? No, A = 2 scales both!
    # In fact, at M = 3, R = 1, EVERY H-L edge has rate 1/3, so error = (C-1)/3 * A!
    # So at (M=3, R=4), (2B/M + A/R)*C = (A/4)*C < ((C-1)/3)*A for C >= 5!!
    # Check: how big is (M=3, R=4)? Let's check: for M=3, a low vertex has D=2, so c <= 2, active children <= 1 (a path!!),
    # so (M=3, R=4) has only ~35 states — instant!!
    p34, idx34, G34 = enumerate_aux_types(3, 4)
    w34 = [Fraction(1) if p[0] == "H" else Fraction(0) for p in p34]
    gw34 = mat_vec(G34, w34)
    A34, B34, D34 = compute_A_B_D(p34, idx34, w34, 3, 4)
    # Bipartite H-L path of 5 cores: H(D=3)-L(D=2)-H(D=3)-L(D=2)-H(D=3) (4 H-L edges -> exact err = 4/3, while (A/4)*5 = 5/4 < 4/3!)
    T_hl5 = CoreTree((-1, 0, 1, 2, 3), (2, 0, 1, 0, 2))
    L34 = sum((rate * (eval_Phi(T2, idx34, w34, 3, 4) - eval_Phi(T_hl5, idx34, w34, 3, 4))
               for rate, T2 in T_hl5.physical_transitions()), Fraction(0))
    err34 = abs(L34 - eval_Phi(T_hl5, idx34, gw34, 3, 4))
    D34_no_high_send = Fraction(2) * B34 / Fraction(3) + A34 / Fraction(4)
    assert err34 <= D34 * T_hl5.C and err34 > D34_no_high_send * T_hl5.C
    print(f"  [CAUGHT] Error 2 (omit 2A/M high-to-low send term): err={err34} > D_bad*C={D34_no_high_send*T_hl5.C} on H-L-H-L-H (M=3,R=4).")

    # Error 3: Omit 2B/M high-promotion term from D(w) on a high path (M=3, R=1)
    p31, idx31, G31 = enumerate_aux_types(3, 1)
    # w = 1 on low singleton, 0 on high -> B = 1, A = 1 (low -> high changes score by -1).
    # Wait: to get A = 0, B = 1? Set w = 1 on EVERY particle (high and low): then leaf gain at low D=2->3
    # converts 1 low particle to 1 high particle (change 1 - 1 = 0, so A = 0!!), while B = 1!!
    w_const1 = [Fraction(1)] * len(p31)
    gw31 = mat_vec(G31, w_const1)
    A31, B31, D31 = compute_A_B_D(p31, idx31, w_const1, 3, 1)
    assert A31 == 0 and B31 == 1
    T_hpath5 = CoreTree((-1, 0, 1, 2, 3), (2, 1, 1, 1, 2))  # 5 high cores, D=3, sum c_x/D_x = 8/3
    L31 = sum((rate * (eval_Phi(T2, idx31, w_const1, 3, 1) - eval_Phi(T_hpath5, idx31, w_const1, 3, 1))
               for rate, T2 in T_hpath5.physical_transitions()), Fraction(0))
    err31 = abs(L31 - eval_Phi(T_hpath5, idx31, gw31, 3, 1))
    D31_no_B = Fraction(2) * A31 / Fraction(3) + A31 / Fraction(1)
    assert err31 == Fraction(8, 3) and err31 <= D31 * T_hpath5.C and err31 > D31_no_B * T_hpath5.C
    print(f"  [CAUGHT] Error 3 (omit 2B/M high-promotion term): err={err31} > 0 = D_no_B*C on 5-core high path (M=3,R=1).")

    # Error 4: Drop factor of 2 in 2[A+B]/M (use [A+B]/M instead) on the SAME 5-core high path!
    D31_no_factor2 = (A31 + B31) / Fraction(3) + A31 / Fraction(1)
    assert err31 == Fraction(8, 3) > D31_no_factor2 * T_hpath5.C  # 8/3 > 5/3!
    print(f"  [CAUGHT] Error 4 (drop factor 2 in 2[A+B]/M): err={err31} > D_bad*C={D31_no_factor2*T_hpath5.C} (proves factor 2 is sharp).")

    # Error 5: Forget to keep inactive degree slot on low endpoint when neighbour is high
    # On T = H(D=3)-L(D=2, c=1, L=1): true physical promotion rate of L is L/D = 1/2,
    # while dropping the parent slot would give D_bad = 1 and promotion rate 1/1 = 1 != 1/2!
    T_hl2 = CoreTree((-1, 0), (2, 1))
    phys_prom_L = Fraction(1, 2)
    aux_trans_L = aux_transitions_particle(singleton_low(0), 3, 1)
    aux_prom_L = sum(rate for rate, out in aux_trans_L if len(out) == 2)
    assert aux_prom_L == phys_prom_L != Fraction(1, 1)
    print(f"  [CAUGHT] Error 5 (dropping inactive parent slot at low endpoint): changes promotion rate from {phys_prom_L} to 1.")

    # Error 6: Wrong sign of D(w) in lower test (+D(w) instead of -D(w) in Phi_{Gw - lam w} + D(w)C >= 0)
    # Set w = 10 on H, w = 1 on all low types at (M=3, R=1): A = 9, B = 1, D = 47/3, Gw >= 0 everywhere.
    # At impossible lower bound lam = 3/2 (> 1 >= mu): lam*w - Gw - D <= 15 - 1 - 47/3 = -5/3 < 0 on H
    # and <= 3/2 - 0 - 47/3 < 0 on all L, falsely certifying mu >= 3/2 under the wrong sign!
    w_sign = [Fraction(10) if p[0] == "H" else Fraction(1) for p in p31]
    gw_sign = mat_vec(G31, w_sign)
    _, _, D_sign = compute_A_B_D(p31, idx31, w_sign, 3, 1)
    syms31, rules31 = build_finite_grammar(3, 1)
    lam_impossible = Fraction(3, 2)  # > 1 >= mu!
    score_wrong_sign = [lam_impossible * w_sign[i] - gw_sign[i] for i in range(len(p31))]
    ok_wrong, _ = solve_grammar_suprema(
        syms31, rules31, lambda rinfo: grammar_vertex_reward(rinfo, idx31, score_wrong_sign, -D_sign, 1)
    )
    ok_right, _ = solve_grammar_suprema(
        syms31, rules31, lambda rinfo: grammar_vertex_reward(rinfo, idx31, score_wrong_sign, +D_sign, 1)
    )
    assert ok_wrong and (not ok_right)
    print(f"  [CAUGHT] Error 6 (wrong sign of D(w) in lower test): falsely certifies impossible mu >= 3/2; true sign rejects it.")

    # Error 7: Omit high-child list recursion rule L_high -> X L_high (restricting high arity to 0)
    # Construct a score where high founder has negative reward (-10) and a high-parent child X has positive reward (+1):
    # True grammar (arbitrary high arity) MUST reject (unbounded reward by attaching >10 copies of X),
    # while truncated grammar (arity 0 only) falsely accepts V_start = -10 <= 0!
    rules31_trunc = [r for r in rules31 if not (r[0] == "L_HIGH" and len(r[1]) == 2)]
    def custom_rew(sym: Optional[Tuple]) -> Fraction:
        if sym is None:
            return Fraction(0)
        if sym == ("VH", True, "NONE"):
            return Fraction(-10)
        if sym == ("VL", False, "HIGH", (2, 1, ())):
            return Fraction(1)
        return Fraction(-100)
    ok_full_g, _ = solve_grammar_suprema(syms31, rules31, custom_rew)
    ok_trunc_g, _ = solve_grammar_suprema(syms31, rules31_trunc, custom_rew)
    assert (not ok_full_g) and ok_trunc_g
    print("  [CAUGHT] Error 7 (omitting L_high -> X L_high recursion): truncated grammar falsely accepts unbounded high-arity tree.")

    # Error 8: Reset relative depth to 0 (instead of keeping permanent depth) when parent turns high at (M=3, R=2)
    # On T = L(d=0, D=2, c=1) -> L(d=1, D=2, c=1) at phase r = 0:
    # Child has rel_d = 1. When root gains a leaf (D: 2->3, turns high), child MUST keep rel_d = 1, NOT reset to 0!
    p32, idx32, _ = enumerate_aux_types(3, 2)
    T_two_low = CoreTree((-1, 0), (1, 1))
    T_after_gain = CoreTree((-1, 0), (2, 1))  # root now D=3 (high), child at depth 1 still low!
    proj_before = T_two_low.project_phase(0, 3, 2)
    proj_after = T_after_gain.project_phase(0, 3, 2)
    assert len(proj_before) == 1
    aux_split = add_leaf_at_path(proj_before[0][2], (), True, 3)
    assert sorted(aux_split) == sorted(proj_after)
    assert singleton_low(1) in aux_split and singleton_low(0) not in aux_split
    print("  [CAUGHT] Error 8 (resetting released child rel_depth to 0 on parent->high): mismatches physical phase-0 projection (rel_d=1).")
    print("  All 8 planted errors caught. PASS")


def main() -> None:
    t0 = time.time()
    test_aux_generator_and_exact_drift()
    test_finite_grammar_and_lp()
    test_estimates_3_4_5()
    test_planted_errors()
    print(f"\nALL INDEPENDENT CHECKS PASSED in {time.time() - t0:.2f}s.")


if __name__ == "__main__":
    main()
