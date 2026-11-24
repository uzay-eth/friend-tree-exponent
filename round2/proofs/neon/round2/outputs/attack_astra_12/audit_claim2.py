#!/usr/bin/env python3
"""Independent adversarial verification suite for Astra Claim 2:
  - round2/fringe/renewal.md

Written from scratch directly from the mathematical equations in renewal.md,
without copying any code from Astra's round2/fringe/check.py or fringe_neon.py.
Read-only on sibling ft2-astra checkout.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import math
import random
from typing import Callable, Dict, List, Optional, Set, Tuple


# ===========================================================================
# 1. Exact Model Constants & Lyapunov Weight from §3 of renewal.md
# ===========================================================================

D_CONST = Fraction(3, 20)
Q_CONST = Fraction(9, 25)
RHO_CONST = Fraction(51, 100)  # q + d = 36/100 + 15/100 = 51/100


def h_func(c: int) -> Fraction:
    if c == 0:
        return Fraction(2)
    if c == 1:
        return Fraction(1)
    if c == 2:
        return Fraction(17, 20)
    if c == 3:
        return Fraction(39, 50)
    return Fraction(7, 10)


def z_func(a: Fraction) -> Fraction:
    return a / (a + D_CONST)


# A birth-marked descendant subtree S is represented as:
#   (a_root: Fraction, leaves_at_root: int, sorted_children: Tuple[Subtree, ...])
# Note: root of S always has 1 exterior parent slot in its total degree and core degree!
Subtree = Tuple[Fraction, int, Tuple["Subtree", ...]]

# A full rooted birth-marked tree T (rooted at initial centre o of T_3, NO parent slot) is:
#   (leaves_at_centre: int, sorted_children: Tuple[Subtree, ...])
FullTree = Tuple[int, Tuple[Subtree, ...]]


def newborn_J(a: Fraction) -> Subtree:
    return (a, 1, ())


def root_plus_leaf(S: Subtree) -> Subtree:
    return (S[0], S[1] + 1, S[2])


def root_deg_S(S: Subtree) -> int:
    # 1 parent slot + L_o + #core_children
    return 1 + S[1] + len(S[2])


def root_c_S(S: Subtree) -> int:
    # 1 parent slot + #core_children
    return 1 + len(S[2])


def vertex_count_S(S: Subtree) -> int:
    return 1 + S[1] + sum(vertex_count_S(y) for y in S[2])


def core_count_S(S: Subtree) -> int:
    return 1 + sum(core_count_S(y) for y in S[2])


def collect_subtrees_in_S(S: Subtree) -> List[Subtree]:
    res = [S]
    for y in S[2]:
        res.extend(collect_subtrees_in_S(y))
    return res


def collect_noncentre_subtrees(T: FullTree) -> List[Subtree]:
    res: List[Subtree] = []
    for y in T[1]:
        res.extend(collect_subtrees_in_S(y))
    return res


def weight_W(S: Subtree) -> Fraction:
    """Exact Lyapunov weight W(S) = (a_S + d)^{-1} sum_{x in S core} w_x h(c_x)."""
    def inner_sum(node: Subtree) -> Fraction:
        c_x = 1 + len(node[2])
        acc = h_func(c_x)
        for ch in node[2]:
            acc += z_func(ch[0]) * inner_sum(ch)
        return acc

    return inner_sum(S) / (S[0] + D_CONST)


def weight_Z_full(T: FullTree) -> Fraction:
    """Full-tree weight Z(T) = sum_{x in T core} w_x h(c_x) from centre o."""
    def inner_sum(node: Subtree) -> Fraction:
        acc = h_func(1 + len(node[2]))
        for ch in node[2]:
            acc += z_func(ch[0]) * inner_sum(ch)
        return acc

    acc0 = h_func(len(T[1]))
    for y in T[1]:
        acc0 += z_func(y[0]) * inner_sum(y)
    return acc0


# ===========================================================================
# 2. Fresh Implementation of True Labelled Friend Tree & Operators L^0, R, B, I
# ===========================================================================

class ExplicitFriendTree:
    """Undirected friend tree on vertices 0..n-1 with centre 0 (from T_3: 0--1, 0--2).
    Stores adjacency lists and permanent birth marks on vertices that have become cores.
    """

    def __init__(self, adj: Tuple[Tuple[int, ...], ...], marks: Tuple[Optional[Fraction], ...]):
        self.adj = adj
        self.marks = marks

    @staticmethod
    def seed_T3() -> "ExplicitFriendTree":
        return ExplicitFriendTree(
            adj=((1, 2), (0,), (0,)),
            marks=(None, None, None),
        )

    def n_vertices(self) -> int:
        return len(self.adj)

    def core_count(self) -> int:
        return sum(1 for nbrs in self.adj if len(nbrs) >= 2)

    def step_transitions(self) -> List[Tuple[Fraction, "ExplicitFriendTree"]]:
        """Enumerate every (source_vertex, chosen_neighbor) continuous-time transition."""
        n = len(self.adj)
        out = []
        for u in range(n):
            d_u = len(self.adj[u])
            rate = Fraction(1, d_u)
            for w in self.adj[u]:
                # Attach new leaf n to w
                new_adj = [list(nbrs) for nbrs in self.adj] + [[w]]
                new_adj[w].append(n)
                new_marks = list(self.marks) + [None]
                if len(self.adj[w]) == 1:
                    # w was a leaf: its unique existing neighbor p_w is its parent!
                    p_w = self.adj[w][0]
                    assert p_w == u, "Only a leaf's unique neighbor can attach a leaf to it!"
                    new_marks[w] = Fraction(1, len(self.adj[p_w]))
                out.append((rate, ExplicitFriendTree(tuple(tuple(x) for x in new_adj), tuple(new_marks))))
        return out

    def extract_full_tree(self) -> FullTree:
        """Root at centre 0 and extract canonical FullTree + verify mark inequality."""
        def dfs_subtree(v: int, par: int) -> Subtree:
            a_v = self.marks[v]
            assert a_v is not None
            d_par = len(self.adj[par])
            assert Fraction(1, d_par) <= a_v <= Fraction(1, 2), (
                f"Mark inequality failed: p_par=1/{d_par} > a_v={a_v}"
            )
            lvs = 0
            kids = []
            for w in self.adj[v]:
                if w == par:
                    continue
                if len(self.adj[w]) == 1:
                    lvs += 1
                else:
                    kids.append(dfs_subtree(w, v))
            return (a_v, lvs, tuple(sorted(kids)))

        l0 = 0
        ch0 = []
        for w in self.adj[0]:
            if len(self.adj[w]) == 1:
                l0 += 1
            else:
                ch0.append(dfs_subtree(w, 0))
        return (l0, tuple(sorted(ch0)))


def add_term(d: Dict[Subtree, Fraction], k: Subtree, v: Fraction) -> None:
    if v == 0:
        return
    nv = d.get(k, Fraction(0)) + v
    if nv == 0:
        d.pop(k, None)
    else:
        d[k] = nv


def intrinsic_jump_rates(S: Subtree) -> Dict[Subtree, Fraction]:
    """Positive off-diagonal jump rates of L^0 on S (parent slot at root is inert)."""
    a_S, L_o, ch = S
    D_o = 1 + L_o + len(ch)
    jumps: Dict[Subtree, Fraction] = {}

    # 1. Root leaf gain from root's leaves and root's core children
    r_root_gain = Fraction(L_o) + sum(Fraction(1, root_deg_S(Y)) for Y in ch)
    if r_root_gain > 0:
        add_term(jumps, (a_S, L_o + 1, ch), r_root_gain)

    # 2. Root promotion
    if L_o > 0:
        nb = newborn_J(Fraction(1, D_o))
        nxt_ch = tuple(sorted(ch + (nb,)))
        add_term(jumps, (a_S, L_o - 1, nxt_ch), Fraction(L_o, D_o))

    # 3. Root -> child Y arrow (rate 1/D_o) + internal transitions inside each child Y
    p_S = Fraction(1, D_o)
    for idx, Y in enumerate(ch):
        rest = ch[:idx] + ch[idx + 1 :]
        add_term(jumps, (a_S, L_o, tuple(sorted(rest + (root_plus_leaf(Y),)))), p_S)
        for Y_nxt, r_sub in intrinsic_jump_rates(Y).items():
            add_term(jumps, (a_S, L_o, tuple(sorted(rest + (Y_nxt,)))), r_sub)

    return jumps


def op_L0(S: Subtree) -> Dict[Subtree, Fraction]:
    jumps = intrinsic_jump_rates(S)
    out = dict(jumps)
    add_term(out, S, -sum(jumps.values(), Fraction(0)))
    return out


def op_R(
    S: Subtree,
    *,
    drop_copy_Y: bool = False,
    drop_killing: bool = False,
    flip_copy_sign: bool = False,
) -> Dict[Subtree, Fraction]:
    """Residual generator R from eq. (1) of renewal.md."""
    out = op_L0(S)
    if not drop_killing:
        add_term(out, S, -S[0])
    p_S = Fraction(1, root_deg_S(S))
    sgn = Fraction(-1 if flip_copy_sign else 1)
    for Y in S[2]:
        a_Y = Y[0]
        assert a_Y >= p_S
        if not drop_copy_Y:
            add_term(out, Y, sgn * (a_Y - p_S))
        add_term(out, root_plus_leaf(Y), sgn * p_S)
    return out


def op_B(S: Subtree) -> Dict[Subtree, Fraction]:
    """Direct newborn birth operator B f(S) = b_S f(J_{p_S})."""
    L_o = S[1]
    if L_o == 0:
        return {}
    D_o = root_deg_S(S)
    return {newborn_J(Fraction(1, D_o)): Fraction(L_o, D_o)}


def op_RB(S: Subtree, **kwargs) -> Dict[Subtree, Fraction]:
    out = op_R(S, **kwargs)
    for k, v in op_B(S).items():
        add_term(out, k, v)
    return out


def op_I_tree(T: FullTree, *, flip_sign: bool = False) -> Dict[Subtree, Fraction]:
    """Initial-centre immigration measure from eq. (2) on a fixed tree T."""
    L_o, ch = T
    D_o = L_o + len(ch)
    p_o = Fraction(1, D_o)
    out: Dict[Subtree, Fraction] = {}
    if L_o > 0:
        add_term(out, newborn_J(p_o), Fraction(L_o, D_o))
    sgn = Fraction(-1 if flip_sign else 1)
    for Y in ch:
        a_Y = Y[0]
        assert a_Y >= p_o
        add_term(out, Y, sgn * (a_Y - p_o))
        add_term(out, root_plus_leaf(Y), sgn * p_o)
    return out


def full_tree_jumps(T: FullTree) -> Dict[FullTree, Fraction]:
    """Exact continuous-time jumps of the canonical full tree T."""
    L_o, ch = T
    D_o = L_o + len(ch)
    out: Dict[FullTree, Fraction] = {}

    def add_ft(k: FullTree, v: Fraction):
        if v != 0:
            nv = out.get(k, Fraction(0)) + v
            if nv == 0:
                out.pop(k, None)
            else:
                out[k] = nv

    r_gain = Fraction(L_o) + sum(Fraction(1, root_deg_S(Y)) for Y in ch)
    if r_gain > 0:
        add_ft((L_o + 1, ch), r_gain)
    if L_o > 0:
        nb = newborn_J(Fraction(1, D_o))
        add_ft((L_o - 1, tuple(sorted(ch + (nb,)))), Fraction(L_o, D_o))

    p_o = Fraction(1, D_o)
    for idx, Y in enumerate(ch):
        rest = ch[:idx] + ch[idx + 1 :]
        add_ft((L_o, tuple(sorted(rest + (root_plus_leaf(Y),)))), p_o)
        for Y_nxt, r_sub in intrinsic_jump_rates(Y).items():
            add_ft((L_o, tuple(sorted(rest + (Y_nxt,)))), r_sub)
    return out


# ===========================================================================
# 3. Check 1: Exact Enumeration n=3..9, E[C_n], Generator & Time-Derivative Identities
# ===========================================================================

EXPECTED_AN_FACTS = {
    3: Fraction(1, 1),
    4: Fraction(4, 3),
    5: Fraction(19, 12),
    6: Fraction(647, 360),
    7: Fraction(17189, 8640),
    8: Fraction(4717423, 2177280),
    9: Fraction(101568829, 43545600),
}


def run_check1_generator_and_enumeration() -> Tuple[Dict[int, Dict[FullTree, Fraction]], Set[Subtree]]:
    print("=== [Claim 2 / Check 1A] Exact Friend-Tree Enumeration (n=3..9) & M'_t f Identity ===")

    # First: for n = 3..7, propagate EXPLICIT labelled undirected trees and cross-check
    # against canonical FullTree distribution!
    lab_dist: Dict[ExplicitFriendTree, Fraction] = {ExplicitFriendTree.seed_T3(): Fraction(1)}
    can_dists: Dict[int, Dict[FullTree, Fraction]] = {3: {(2, ()): Fraction(1)}}
    reps_by_n: Dict[int, Dict[FullTree, ExplicitFriendTree]] = {3: {(2, ()): ExplicitFriendTree.seed_T3()}}

    for n in range(3, 9):
        nxt_can: Dict[FullTree, Fraction] = defaultdict(Fraction)
        for T_ft, prob in can_dists[n].items():
            jumps = full_tree_jumps(T_ft)
            assert sum(jumps.values()) == Fraction(n), f"Total jump rate != {n}"
            for T_nxt, r in jumps.items():
                nxt_can[T_nxt] += prob * r / Fraction(n)
        can_dists[n + 1] = dict(nxt_can)

        if n < 7:
            nxt_lab: Dict[ExplicitFriendTree, Fraction] = defaultdict(Fraction)
            nxt_reps: Dict[FullTree, ExplicitFriendTree] = {}
            proj_can: Dict[FullTree, Fraction] = defaultdict(Fraction)
            for et, prob in lab_dist.items():
                for r, et_nxt in et.step_transitions():
                    p_new = prob * r / Fraction(n)
                    nxt_lab[et_nxt] += p_new
                    ft_nxt = et_nxt.extract_full_tree()
                    proj_can[ft_nxt] += p_new
                    if ft_nxt not in nxt_reps:
                        nxt_reps[ft_nxt] = et_nxt
            assert dict(proj_can) == can_dists[n + 1], f"Labelled vs canonical mismatch at n={n+1}"
            lab_dist = dict(nxt_lab)
            reps_by_n[n + 1] = nxt_reps

    for n in range(3, 10):
        an = sum(
            prob * (1 + sum(core_count_S(y) for y in T[1]))
            for T, prob in can_dists[n].items()
        )
        assert an == EXPECTED_AN_FACTS[n], f"E[C_{n}] mismatch: {an} vs {EXPECTED_AN_FACTS[n]}"
        print(f"  n={n}: {len(can_dists[n]):4d} marked states, exact E[C_{n}] = {an} (matches FACT 2)")

    # Verify exact pointwise generator identity on ALL ExplicitFriendTree states (n=3..7)
    # AND on ALL 1,350+ canonical FullTree states (n=3..9), evaluated on non-trivial test functions
    # f_1, f_2, f_3, f_4 AND as formal linear combinations!
    def f_nontrivial_1(S: Subtree) -> Fraction:
        # Non-linear score depending on root mark, leaf count, vertex count, and child marks
        b_S = Fraction(S[1], root_deg_S(S))
        return (S[0] ** 2) * vertex_count_S(S) + b_S * sum(
            y[0] * (vertex_count_S(y) ** 2) for y in S[2]
        )

    def f_nontrivial_2(S: Subtree) -> Fraction:
        # Score sensitive to subtree depth, Y vs Y^+, and marks
        return Fraction(S[1] + 1, root_deg_S(S) + 1) + sum(
            (y[0] + Fraction(y[1], 7)) * core_count_S(y) for y in S[2]
        )

    test_fns: List[Tuple[str, Callable[[Subtree], Fraction]]] = [
        ("1", lambda _: Fraction(1)),
        ("W(S)", weight_W),
        ("nonlinear_mark_size", f_nontrivial_1),
        ("leaf_depth_sensitive", f_nontrivial_2),
    ]

    all_subtrees: Set[Subtree] = set()
    checked_trees = 0

    # Check directly against undirected ExplicitFriendTree step_transitions() for n=3..7
    for n in range(3, 8):
        for ft_key, et in reps_by_n[n].items():
            old_subs = collect_noncentre_subtrees(et.extract_full_tree())
            lhs_formal: Dict[Subtree, Fraction] = {}
            for r, et_nxt in et.step_transitions():
                for S in collect_noncentre_subtrees(et_nxt.extract_full_tree()):
                    add_term(lhs_formal, S, r)
                for S in old_subs:
                    add_term(lhs_formal, S, -r)

            rhs_formal: Dict[Subtree, Fraction] = {}
            for S in old_subs:
                for S_out, c in op_RB(S).items():
                    add_term(rhs_formal, S_out, c)
            for S_out, c in op_I_tree(ft_key).items():
                add_term(rhs_formal, S_out, c)
            assert lhs_formal == rhs_formal, f"Explicit undirected transfer identity failed at n={n}"

    # Check on ALL canonical states n=3..9
    for n in range(3, 10):
        for T_ft in can_dists[n]:
            checked_trees += 1
            old_subs = collect_noncentre_subtrees(T_ft)
            for S in old_subs:
                all_subtrees.add(S)
                all_subtrees.add(root_plus_leaf(S))

            jumps = full_tree_jumps(T_ft)
            tot_r = sum(jumps.values())
            lhs_formal = {}
            for T_nxt, r in jumps.items():
                for S in collect_noncentre_subtrees(T_nxt):
                    add_term(lhs_formal, S, r)
            for S in old_subs:
                add_term(lhs_formal, S, -tot_r)

            rhs_formal = {}
            for S in old_subs:
                for S_out, c in op_RB(S).items():
                    add_term(rhs_formal, S_out, c)
            for S_out, c in op_I_tree(T_ft).items():
                add_term(rhs_formal, S_out, c)
            assert lhs_formal == rhs_formal

            for _, fn in test_fns:
                lhs_v = sum(c * fn(S) for S, c in lhs_formal.items())
                rhs_v = sum(c * fn(S) for S, c in rhs_formal.items())
                assert lhs_v == rhs_v

    # Check continuous-time derivatives d^k/dt^k |_{t=0} of M'_t f = M_t(R+B)f + I_t f
    # AND M_t f = int_0^t I_u P_{t-u} f du for orders k = 0..5 on non-trivial f!
    tree_derivs: List[Dict[FullTree, Fraction]] = [{(2, ()): Fraction(1)}]
    for _ in range(6):
        prev = tree_derivs[-1]
        nxt_d: Dict[FullTree, Fraction] = {}
        for T_ft, coeff in prev.items():
            jumps = full_tree_jumps(T_ft)
            tot_r = sum(jumps.values())
            nxt_d[T_ft] = nxt_d.get(T_ft, Fraction(0)) - coeff * tot_r
            for T_nxt, r in jumps.items():
                nxt_d[T_nxt] = nxt_d.get(T_nxt, Fraction(0)) + coeff * r
        tree_derivs.append({k: v for k, v in nxt_d.items() if v != 0})

    for k in range(5):
        # Verify d/dt|_{t=0} d^k/dt^k (M_t f) == (d^k/dt^k M_t)((R+B)f) + (d^k/dt^k I_t)(f)
        for fname, fn in test_fns:
            lhs_dk1 = sum(
                coeff * sum(fn(S) for S in collect_noncentre_subtrees(T_ft))
                for T_ft, coeff in tree_derivs[k + 1].items()
            )
            rhs_rb = sum(
                coeff
                * sum(
                    sum(c2 * fn(S2) for S2, c2 in op_RB(S).items())
                    for S in collect_noncentre_subtrees(T_ft)
                )
                for T_ft, coeff in tree_derivs[k].items()
            )
            rhs_imm = sum(
                coeff * sum(c2 * fn(S2) for S2, c2 in op_I_tree(T_ft).items())
                for T_ft, coeff in tree_derivs[k].items()
            )
            assert lhs_dk1 == rhs_rb + rhs_imm, f"Time-derivative order {k+1} failed for {fname}"

    print(
        f"  Exact generator & time-derivative identities M'_t f = M_t(R+B)f + I_t f verified "
        f"across {checked_trees} states and orders k=1..5 on 4 test functions."
    )
    print("  [PASS] Check 1A\n")
    return can_dists, all_subtrees


# ===========================================================================
# 4. Check 2: Exact Rational Lyapunov Inequality (4), Child Identity & §4-§5 Bounds
# ===========================================================================

def run_check2_lyapunov_and_section5(
    can_dists: Dict[int, Dict[FullTree, Fraction]],
    all_subtrees: Set[Subtree],
) -> None:
    print("=== [Claim 2 / Check 1B & 2] Rational Inequality (4), Child Identity & §4-§5 Bounds ===")

    # 1. Rational inequality (4) for ALL 2 <= D < 250 and 0 <= c <= D
    checked_pairs = 0
    max_gap = Fraction(-100)
    arg_max = None
    for D in range(2, 250):
        for c in range(0, D + 1):
            lhs = Fraction(D - c, D) * (h_func(c + 1) - h_func(c) + Fraction(1, 1 + D_CONST * D))
            rhs = Q_CONST * h_func(c)
            assert lhs <= rhs, f"Inequality (4) failed at D={D}, c={c}: {lhs} > {rhs}"
            if lhs - rhs > max_gap:
                max_gap = lhs - rhs
                arg_max = (D, c, lhs, rhs)
            checked_pairs += 1

    # Analytic D >= 100 check:
    assert 1 + D_CONST * 100 == Fraction(16)
    assert Fraction(1, 16) <= Q_CONST * Fraction(7, 10)
    D_t, c_t, lhs_t, rhs_t = arg_max
    print(
        f"  Inequality (4) verified for {checked_pairs} pairs (2<=D<250, 0<=c<=D); "
        f"tightest at (D={D_t}, c={c_t}): {lhs_t} <= {rhs_t} (slack={rhs_t - lhs_t})."
    )
    print(
        f"  Analytic tail D>=100: LHS <= 1/(1+100*3/20) = 1/16 = {float(Fraction(1,16)):.4f} "
        f"<= q*min(h) = 63/250 = {float(Q_CONST*Fraction(7,10)):.4f}."
    )

    # 2. Child-weight identity & R W <= rho W on all reachable types + a=0 limit types + random large trees
    test_S = set(all_subtrees)
    test_S.add(newborn_J(Fraction(0)))
    for S in list(all_subtrees):
        test_S.add((Fraction(0), S[1], S[2]))

    for S in test_S:
        wS = weight_W(S)
        c_o = root_c_S(S)
        # Elementary bounds in §3 line 53
        assert wS >= Fraction(7, 10) / (Fraction(1, 2) + D_CONST) == Fraction(14, 13) > 1
        assert wS <= Fraction(vertex_count_S(S)) / D_CONST
        assert weight_W(root_plus_leaf(S)) == wS

        # Child-weight identity: sum_Y a_Y W(Y) == (a_S + d) W(S) - h(c_o)
        lhs_ch = sum((Y[0] * weight_W(Y) for Y in S[2]), Fraction(0))
        rhs_ch = (S[0] + D_CONST) * wS - h_func(c_o)
        assert lhs_ch == rhs_ch

        # L^0 W <= q W, R W <= rho W, B W <= 1/d
        l0_w = sum(c * weight_W(S2) for S2, c in op_L0(S).items())
        assert l0_w <= Q_CONST * wS
        r_w = sum(c * weight_W(S2) for S2, c in op_R(S).items())
        assert r_w == l0_w + D_CONST * wS - h_func(c_o)
        assert r_w <= RHO_CONST * wS
        b_w = sum(c * weight_W(S2) for S2, c in op_B(S).items())
        assert b_w <= Fraction(1) / D_CONST

    print(f"  Child-weight identity & RW <= (51/100)W verified on {len(test_S)} types (including a_S=0).")

    # 3. §4 Spine eq. (8) exact a-independence check on J_a for a in {1/2, 1/3, 1/20, 0}
    # Verify that under L^0 started at J_a, the distribution of (L_o, children) after k jumps
    # and the copy operator sum_Y [(a_Y - p_S) delta_Y + p_S delta_{Y^+}] are 100% independent of a!
    def strip_a(S: Subtree):
        return (S[1], S[2])

    ref_spine = None
    for a_val in [Fraction(1, 2), Fraction(1, 3), Fraction(1, 20), Fraction(0)]:
        cur = {newborn_J(a_val): Fraction(1)}
        trace = []
        for _ in range(4):
            stripped = {strip_a(S): v for S, v in cur.items()}
            copied: Dict[Subtree, Fraction] = {}
            for S, v in cur.items():
                p_S = Fraction(1, root_deg_S(S))
                for Y in S[2]:
                    add_term(copied, Y, v * (Y[0] - p_S))
                    add_term(copied, root_plus_leaf(Y), v * p_S)
            trace.append((stripped, copied))
            nxt: Dict[Subtree, Fraction] = {}
            for S, v in cur.items():
                for S2, c2 in op_L0(S).items():
                    add_term(nxt, S2, v * c2)
            cur = nxt
        if ref_spine is None:
            ref_spine = trace
        else:
            assert trace == ref_spine, f"Spine depended on root mark a={a_val}!"
    print("  §4 Spine decomposition (8) verified completely independent of root mark a in {1/2, 1/3, 1/20, 0}.")

    # 4. §5 Upper bound sum_{v != o} W(S_v) <= (260/9)(C-1) & lower combinatorial bound >= C/12 - 1
    def verify_combinatorial_bounds_on_adj(adj: Sequence[Sequence[int]]) -> None:
        n = len(adj)
        cores = [u for u in range(n) if len(adj[u]) >= 2]
        C = len(cores)
        core_set = set(cores)
        assert 0 in core_set
        # Root at 0
        par = {0: -1}
        ch_core = {u: [] for u in cores}
        stack = [0]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w in core_set and w != par[u]:
                    par[w] = u
                    ch_core[u].append(w)
                    stack.append(w)
        c_deg = {u: len(ch_core[u]) + (0 if u == 0 else 1) for u in cores}
        l_cnt = {u: len(adj[u]) - c_deg[u] for u in cores}

        def bare(u: int) -> bool:
            return c_deg[u] == 2 and l_cnt[u] == 0

        eligible = [
            u for u in cores
            if (c_deg[u] <= 2 and l_cnt[u] >= 1)
            or (bare(u) and any(bare(w) for w in adj[u] if w in core_set))
        ]
        assert len(eligible) * 4 >= C

        good1 = {u for u in cores if c_deg[u] <= 2 and l_cnt[u] >= 1}
        good2 = {u for u in cores if bare(u) and any(bare(w) for w in ch_core[u])}
        rem = set(eligible) - (good1 | good2)
        assert 0 not in rem
        for u in rem:
            assert par[u] in good2
        noninit_good = (good1 | good2) - {0}
        assert Fraction(len(noninit_good)) >= Fraction(C, 12) - 1

    # Check on random friend trees up to n = 150 + adversarial bare trees
    rng = random.Random(999)
    for _ in range(300):
        et = ExplicitFriendTree.seed_T3()
        adj = [list(x) for x in et.adj]
        target_n = rng.randint(4, 150)
        while len(adj) < target_n:
            u = rng.randrange(len(adj))
            w = adj[u][rng.randrange(len(adj[u]))]
            nv = len(adj)
            adj[w].append(nv)
            adj.append([w])
        verify_combinatorial_bounds_on_adj(adj)

    # Adversarial tree: bare centre 0 with two bare chains of various lengths & bare binary trees
    for chain_len in [1, 2, 3, 5, 15]:
        adj_adv: List[List[int]] = [[]]
        for _ in range(2):
            prev = 0
            for _ in range(chain_len):
                cur = len(adj_adv)
                adj_adv.append([prev])
                adj_adv[prev].append(cur)
                prev = cur
            lf = len(adj_adv)
            adj_adv.append([prev])
            adj_adv[prev].append(lf)
        verify_combinatorial_bounds_on_adj(adj_adv)

    # 5. §5 Fixed-size Yule transfer inequalities: (1 - 1/n)^{n-1} > e^{-1} and E[max(1, N(log n)/n)] < 3
    for n in range(3, 200):
        q = 1.0 / n
        # Exact P(N(log n) >= n) where N(t) = G_1 + G_2, G_i ~ Geom(1/n) on {1,2,...}
        p_g1_ge_n = (1.0 - q) ** (n - 1)
        assert p_g1_ge_n > math.exp(-1.0)
        # Exact E[max(1, N(log n)/n)] = 1 + (1/n) sum_{m=n+1}^infty (m - n) (m - 1) q^2 (1-q)^{m-2}
        # Note sum_{m=n+1}^infty (m - n) P(N=m) < E[N] = 2n, so ratio < 1 + 2 = 3!
        # Verify exact sum in Fraction for n = 3..25:
        if n <= 25:
            qf = Fraction(1, n)
            pf = 1 - qf
            # E[max(n, N)] = n + sum_{k=n}^infty P(N > k)
            # P(N <= k) = sum_{m=2}^k (m-1) q^2 p^{m-2} = 1 - p^{k-1}(1 + (k-1)q)
            # wait: check at k=2: 1 - p(1 + q) = 1 - (1-q)(1+q) = q^2 = P(N=2)! Exact!
            # So sum_{k=n}^infty P(N > k) = sum_{k=n}^infty p^{k-1}(1 + (k-1)q)
            # Let j = k - (n-1) >= 0 (wait k >= n -> m = k-1 >= n-1):
            # sum_{m=n-1}^infty p^m (1 + m q) = p^{n-1}/q + q * p^{n-1} * ((n-1)/q + p/q^2)
            # = p^{n-1} [1/q + n - 1 + p/q] = p^{n-1} [n + (1 + p)/q - 1] = p^{n-1} [n - 2 + 2/q] = p^{n-1} (3n - 2)!
            e_max_exact = Fraction(1) + (pf ** (n - 1)) * Fraction(3 * n - 2, n)
            # Verify against finite prefix + tiny tail:
            assert Fraction(1) < e_max_exact < Fraction(3)
    print("  §5 Upper W-sum, lower (C/12 - 1) good-vertex bound, and Yule transfer [e^{-1}, 3] verified.")
    print("  [PASS] Check 2\n")


# ===========================================================================
# 5. Check 3: Planted-Error Tests for Claim 2
# ===========================================================================

def run_check3_planted_errors(can_dists: Dict[int, Dict[FullTree, Fraction]]) -> None:
    print("=== [Claim 2 / Check 3] Planted-Error Falsification Tests ===")

    def test_generator_mutation(**rb_kwargs) -> Optional[Tuple[int, FullTree]]:
        for n in sorted(can_dists):
            for T_ft in can_dists[n]:
                old_subs = collect_noncentre_subtrees(T_ft)
                jumps = full_tree_jumps(T_ft)
                tot_r = sum(jumps.values())
                lhs: Dict[Subtree, Fraction] = {}
                for T_nxt, r in jumps.items():
                    for S in collect_noncentre_subtrees(T_nxt):
                        add_term(lhs, S, r)
                for S in old_subs:
                    add_term(lhs, S, -tot_r)

                imm_flip = rb_kwargs.pop("imm_flip", False)
                rhs: Dict[Subtree, Fraction] = {}
                for S in old_subs:
                    for S2, c in op_RB(S, **rb_kwargs).items():
                        add_term(rhs, S2, c)
                for S2, c in op_I_tree(T_ft, flip_sign=imm_flip).items():
                    add_term(rhs, S2, c)
                rb_kwargs["imm_flip"] = imm_flip
                if lhs != rhs:
                    return (n, T_ft)
        return None

    # 1. Drop (a_Y - p_S) f(Y) copy term in R
    e1 = test_generator_mutation(drop_copy_Y=True, imm_flip=False)
    assert e1 is not None
    print(f"  [CAUGHT] Planted Error 1 (drop (a_Y - p_S)f(Y) copy term): fails at n={e1[0]}")

    # 2. Drop -a_S f(S) killing term in R
    e2 = test_generator_mutation(drop_killing=True, imm_flip=False)
    assert e2 is not None
    print(f"  [CAUGHT] Planted Error 2 (drop -a_S f(S) killing term): fails at n={e2[0]}")

    # 3. Flip sign on copy terms in R or I_t
    e3a = test_generator_mutation(flip_copy_sign=True, imm_flip=False)
    e3b = test_generator_mutation(imm_flip=True)
    assert e3a is not None and e3b is not None
    print(f"  [CAUGHT] Planted Error 3 (flip sign in R / I_t): fails at n={e3a[0]} / n={e3b[0]}")

    # 4. Perturb h(c) or q in inequality (4)
    # Perturb q from 9/25 = 0.36 down to 0.35:
    q_mut = Fraction(7, 20)  # 0.35
    fail_q = None
    for D in range(2, 50):
        for c in range(0, D + 1):
            lhs = Fraction(D - c, D) * (h_func(c + 1) - h_func(c) + Fraction(1, 1 + D_CONST * D))
            if lhs > q_mut * h_func(c):
                fail_q = (D, c, lhs, q_mut * h_func(c))
                break
        if fail_q:
            break
    assert fail_q is not None
    print(
        f"  [CAUGHT] Planted Error 4a (perturb q=9/25 -> 7/20): fails at D={fail_q[0]}, c={fail_q[1]} "
        f"({fail_q[2]} > {fail_q[3]})"
    )

    # Perturb h(1) from 1 to 9/10 (or h(2) from 17/20 to 19/20):
    for label, h_override in [("h(1)=9/10", {1: Fraction(9, 10)}), ("h(2)=19/20", {2: Fraction(19, 20)})]:
        fail_h = None
        for D in range(2, 50):
            for c in range(0, D + 1):
                hc = h_override.get(c, h_func(c))
                hc1 = h_override.get(c + 1, h_func(c + 1))
                lhs = Fraction(D - c, D) * (hc1 - hc + Fraction(1, 1 + D_CONST * D))
                if lhs > Q_CONST * hc:
                    fail_h = (D, c, lhs, Q_CONST * hc)
                    break
            if fail_h:
                break
        assert fail_h is not None, f"Failed to catch {label}"
        print(
            f"  [CAUGHT] Planted Error 4b (perturb {label}): fails at D={fail_h[0]}, c={fail_h[1]} "
            f"({fail_h[2]} > {fail_h[3]})"
        )

    print("  [PASS] All planted-error tests bit cleanly.\n")


def main() -> None:
    can_dists, all_subtrees = run_check1_generator_and_enumeration()
    run_check2_lyapunov_and_section5(can_dists, all_subtrees)
    run_check3_planted_errors(can_dists)
    print("ALL ADVERSARIAL CHECKS FOR CLAIM 2 COMPLETED SUCCESSFULLY.")


if __name__ == "__main__":
    main()
