#!/usr/bin/env python3
r"""Independent exact cold-referee verification for round2/fringe/renewal.md.

Checks (all exact rational arithmetic using fractions.Fraction):
  1. Lyapunov rational inequality (4) for all 0 <= c <= D <= 200, D >= 2,
     plus analytical tail bound for D >= 100.
  2. Exhaustive exact small rooted marked growth histories from T_3:
     - Exhaustive enumeration of ALL individual directed-edge labelled growth
       histories from T_3 to n = 6 (verifying exact tree-law pushforward and
       mark inequality along every individual growth trajectory).
     - Exhaustive enumeration of ALL distinct reachable canonical marked rooted
       trees up to n = 9, cross-checking the brute-force vertex-by-vertex
       labelled physical transitions against the canonical generator.
     - Exact generator transfer identity:
         (L^phys Phi_f)(T) == sum_{v != o} ((R + B) f)(S_v) + I(T; f)
       verified as an exact formal linear combination of canonical types
       (simultaneously proving equality for every test function f).
  3. Exact multistep convolution identity M_t f = int_0^t I_u P_{t-u} f du:
     - Exact Taylor coefficients at t = 0 up to order k = 5 as formal linear
       combinations of types (testing iterated R+B, copy births, direct births,
       killing, and true initial-centre immigration).
     - Exact copy-separation / one-spine decomposition (8) of R_t up to order
       k = 5 at a = 1/2, 1/3, 1/10, and the limit mark a = 0, plus direct
       verification that L^0 and Gamma on J_a are completely independent of a.
     - Exact Volterra / first-direct-birth generation decomposition of P_t:
       (R + B)^k == R^k + sum_{j=0}^{k-1} R^j B (R + B)^{k-1-j} up to k = 5
       at a = 1/2, 1/3, and a = 0.
  4. Exact residual & full-tree Lyapunov weight identities and inequalities:
     - Child-weight identity sum_Y a_Y W(Y) == (a_S + d) W(S) - h(c_o)
       on all reachable types S and auxiliary a_S = 0 copies.
     - L^0 W <= q W, R W == L^0 W + d W - h(c_o) <= rho W - 7/10,
       W >= 14/13 > 1, W <= |S|/d, B W <= 1/d (including a_S = 0).
     - Full-tree L^phys Z <= q Z and exact immigration weight identity
       sum_{Y at o} ((a_Y - p_o) W(Y) + p_o W(Y^+)) == Z - h(c_o).
  5. Exact count comparisons and discrete de-Poissonization inequality:
     - Upper comparison sum_{v != o} W(S_v) <= (260/9)(C - 1) on all reachable
       trees up to n = 9 and adversarial deep trees.
     - Lower comparison tree combinatorics: |E| >= C/4, o not in B_rem,
       |p^{-1}(w)| <= 1 (w != o), |p^{-1}(o)| <= 2, |(G_1 U G_2) \ {o}| >= C/12 - 1
       on all reachable trees up to n = 9 and adversarial bare-chain/comb trees.
     - Exact discrete chain E[C_n] for n = 3..9 and promotion bound <= C_n / n.
  6. Deliberate corruptions (negative controls): verifies that every broken
     variant is immediately caught with a concrete counterexample.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
import time
from typing import Dict, List, Tuple, Set, Optional

# ---------------------------------------------------------------------------
# Parameters from renewal.md Section 3
# ---------------------------------------------------------------------------
D_PARAM = Fraction(3, 20)   # d = 3/20
Q_PARAM = Fraction(9, 25)   # q = 9/25
RHO_PARAM = Q_PARAM + D_PARAM  # rho = 51/100


def h_val(c: int) -> Fraction:
    if c == 0:
        return Fraction(2)
    if c == 1:
        return Fraction(1)
    if c == 2:
        return Fraction(17, 20)
    if c == 3:
        return Fraction(39, 50)
    return Fraction(7, 10)


def z_val(a: Fraction) -> Fraction:
    return a / (a + D_PARAM)


# ---------------------------------------------------------------------------
# Canonical type representations
# ---------------------------------------------------------------------------
TypeS = Tuple[Fraction, int, Tuple["TypeS", ...]]
TreeT = Tuple[int, Tuple[TypeS, ...]]
LinComb = Dict[TypeS, Fraction]


def J_type(a: Fraction) -> TypeS:
    """Newborn core J_a: 1 leaf, 0 core children, 1 parent slot, root mark a."""
    return (a, 1, ())


def add_leaf_root(S: TypeS) -> TypeS:
    """Return S^+ (S with a new leaf attached to its root)."""
    a_S, L_o, ch = S
    return (a_S, L_o + 1, ch)


def deg_S(S: TypeS) -> int:
    return 1 + S[1] + len(S[2])


def core_deg_S(S: TypeS) -> int:
    return 1 + len(S[2])


def size_S(S: TypeS) -> int:
    """Number of vertices in S (root + leaves + descendants, excluding parent slot)."""
    _, L_o, ch = S
    return 1 + L_o + sum(size_S(c) for c in ch)


def core_count_S(S: TypeS) -> int:
    return 1 + sum(core_count_S(c) for c in S[2])


def size_T(T: TreeT) -> int:
    L_o, ch = T
    return 1 + L_o + sum(size_S(c) for c in ch)


def core_count_T(T: TreeT) -> int:
    return 1 + sum(core_count_S(c) for c in T[1])


def all_noninit_subtrees_S(S: TypeS) -> List[TypeS]:
    out = [S]
    for c in S[2]:
        out.extend(all_noninit_subtrees_S(c))
    return out


def all_noninit_subtrees_T(T: TreeT) -> List[TypeS]:
    out: List[TypeS] = []
    for c in T[1]:
        out.extend(all_noninit_subtrees_S(c))
    return out


def lin_add(dst: LinComb, key: TypeS, coeff: Fraction) -> None:
    if coeff == 0:
        return
    val = dst.get(key, Fraction(0)) + coeff
    if val == 0:
        dst.pop(key, None)
    else:
        dst[key] = val


def lin_clean(d: Dict) -> Dict:
    return {k: v for k, v in d.items() if v != 0}


# ---------------------------------------------------------------------------
# Independent brute-force labelled rooted tree representation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LabelledTree:
    """Explicit labelled tree rooted at vertex 0 (initial centre).

    Vertices are 0..n-1.
    parent[0] = -1; parent[v] in 0..n-1 for v > 0.
    birth_mark[v] = Fraction(1, D_parent_at_promotion) for non-root core v,
                    None for root 0 and for leaves.
    """
    n: int
    parent: Tuple[int, ...]
    birth_mark: Tuple[Optional[Fraction], ...]

    def adj(self) -> List[List[int]]:
        a = [[] for _ in range(self.n)]
        for v in range(1, self.n):
            p = self.parent[v]
            a[p].append(v)
            a[v].append(p)
        return a

    def children(self) -> List[List[int]]:
        ch = [[] for _ in range(self.n)]
        for v in range(1, self.n):
            ch[self.parent[v]].append(v)
        return ch

    def degrees(self) -> List[int]:
        deg = [0] * self.n
        for v in range(1, self.n):
            deg[v] += 1
            deg[self.parent[v]] += 1
        return deg

    def attach_leaf(self, target: int) -> "LabelledTree":
        deg = self.degrees()
        new_parent = self.parent + (target,)
        marks = list(self.birth_mark) + [None]
        if deg[target] == 1:
            u = self.parent[target]
            assert u >= 0 and deg[u] >= 2
            marks[target] = Fraction(1, deg[u])
        return LabelledTree(self.n + 1, new_parent, tuple(marks))

    def physical_transitions(self) -> List[Tuple[Fraction, "LabelledTree"]]:
        adj = self.adj()
        deg = [len(adj[w]) for w in range(self.n)]
        out: List[Tuple[Fraction, LabelledTree]] = []
        for w in range(self.n):
            rate = Fraction(1, deg[w])
            for z in adj[w]:
                out.append((rate, self.attach_leaf(z)))
        return out

    def canonical_S(self, v: int, ch: List[List[int]], deg: List[int]) -> TypeS:
        assert v != 0 and deg[v] >= 2
        a_v = self.birth_mark[v]
        assert a_v is not None
        leaves = 0
        core_ch: List[TypeS] = []
        for c in ch[v]:
            if deg[c] == 1:
                leaves += 1
            else:
                core_ch.append(self.canonical_S(c, ch, deg))
        core_ch.sort()
        return (a_v, leaves, tuple(core_ch))

    def canonical_T(self) -> TreeT:
        ch = self.children()
        deg = self.degrees()
        assert deg[0] >= 2
        leaves = 0
        core_ch: List[TypeS] = []
        for c in ch[0]:
            if deg[c] == 1:
                leaves += 1
            else:
                core_ch.append(self.canonical_S(c, ch, deg))
        core_ch.sort()
        return (leaves, tuple(core_ch))

    def noninit_subtrees(self) -> List[TypeS]:
        ch = self.children()
        deg = self.degrees()
        return [self.canonical_S(v, ch, deg) for v in range(1, self.n) if deg[v] >= 2]


def T3_labelled() -> LabelledTree:
    return LabelledTree(n=3, parent=(-1, 0, 0), birth_mark=(None, None, None))


def T3_canonical() -> TreeT:
    return (2, ())


# ---------------------------------------------------------------------------
# Canonical generators: L^0, R, B, I, and canonical physical generator
# ---------------------------------------------------------------------------
def intrinsic_transitions_S(S: TypeS) -> Dict[TypeS, Fraction]:
    a_S, L_o, ch = S
    D_o = 1 + L_o + len(ch)
    out: Dict[TypeS, Fraction] = {}

    gain_rate_o = Fraction(L_o) + sum(Fraction(1, deg_S(Y)) for Y in ch)
    if gain_rate_o > 0:
        S_plus = (a_S, L_o + 1, ch)
        out[S_plus] = out.get(S_plus, Fraction(0)) + gain_rate_o

    if L_o > 0:
        prom_rate_o = Fraction(L_o, D_o)
        new_child = J_type(Fraction(1, D_o))
        new_ch = tuple(sorted(ch + (new_child,)))
        S_prom = (a_S, L_o - 1, new_ch)
        out[S_prom] = out.get(S_prom, Fraction(0)) + prom_rate_o

    for i, Y in enumerate(ch):
        prefix = ch[:i]
        suffix = ch[i + 1:]

        Y_plus = add_leaf_root(Y)
        ch_plus = tuple(sorted(prefix + (Y_plus,) + suffix))
        S_yplus = (a_S, L_o, ch_plus)
        out[S_yplus] = out.get(S_yplus, Fraction(0)) + Fraction(1, D_o)

        for Y_next, rate in intrinsic_transitions_S(Y).items():
            ch_next = tuple(sorted(prefix + (Y_next,) + suffix))
            S_next = (a_S, L_o, ch_next)
            out[S_next] = out.get(S_next, Fraction(0)) + rate

    return lin_clean(out)


def apply_L0_basis(S: TypeS) -> LinComb:
    trans = intrinsic_transitions_S(S)
    res: LinComb = {}
    total = Fraction(0)
    for S_next, rate in trans.items():
        lin_add(res, S_next, rate)
        total += rate
    lin_add(res, S, -total)
    return res


def apply_Gamma_basis(S: TypeS) -> LinComb:
    _, L_o, ch = S
    p_S = Fraction(1, 1 + L_o + len(ch))
    res: LinComb = {}
    for Y in ch:
        a_Y = Y[0]
        assert a_Y >= p_S, f"Mark inequality violated: a_Y={a_Y} < p_S={p_S}"
        lin_add(res, Y, a_Y - p_S)
        lin_add(res, add_leaf_root(Y), p_S)
    return res


def apply_R_basis(S: TypeS) -> LinComb:
    res = apply_L0_basis(S)
    lin_add(res, S, -S[0])
    for Y, coeff in apply_Gamma_basis(S).items():
        lin_add(res, Y, coeff)
    return res


def apply_B_basis(S: TypeS) -> LinComb:
    _, L_o, ch = S
    D_o = 1 + L_o + len(ch)
    res: LinComb = {}
    if L_o > 0:
        lin_add(res, J_type(Fraction(1, D_o)), Fraction(L_o, D_o))
    return res


def apply_RB_basis(S: TypeS) -> LinComb:
    res = apply_R_basis(S)
    for k, v in apply_B_basis(S).items():
        lin_add(res, k, v)
    return res


def apply_I_basis(T: TreeT) -> LinComb:
    L_o, ch = T
    D_o = L_o + len(ch)
    p_o = Fraction(1, D_o)
    b_o = Fraction(L_o, D_o)
    res: LinComb = {}
    if L_o > 0:
        lin_add(res, J_type(p_o), b_o)
    for Y in ch:
        a_Y = Y[0]
        assert a_Y >= p_o, f"Root mark inequality violated: a_Y={a_Y} < p_o={p_o}"
        lin_add(res, Y, a_Y - p_o)
        lin_add(res, add_leaf_root(Y), p_o)
    return res


def canonical_physical_transitions_T(T: TreeT) -> Dict[TreeT, Fraction]:
    L_o, ch = T
    D_o = L_o + len(ch)
    out: Dict[TreeT, Fraction] = {}

    gain_rate_o = Fraction(L_o) + sum(Fraction(1, deg_S(Y)) for Y in ch)
    if gain_rate_o > 0:
        T_plus = (L_o + 1, ch)
        out[T_plus] = out.get(T_plus, Fraction(0)) + gain_rate_o

    if L_o > 0:
        new_ch = tuple(sorted(ch + (J_type(Fraction(1, D_o)),)))
        T_prom = (L_o - 1, new_ch)
        out[T_prom] = out.get(T_prom, Fraction(0)) + Fraction(L_o, D_o)

    for i, Y in enumerate(ch):
        prefix = ch[:i]
        suffix = ch[i + 1:]
        ch_plus = tuple(sorted(prefix + (add_leaf_root(Y),) + suffix))
        T_yplus = (L_o, ch_plus)
        out[T_yplus] = out.get(T_yplus, Fraction(0)) + Fraction(1, D_o)
        for Y_next, rate in intrinsic_transitions_S(Y).items():
            ch_next = tuple(sorted(prefix + (Y_next,) + suffix))
            T_next = (L_o, ch_next)
            out[T_next] = out.get(T_next, Fraction(0)) + rate

    return lin_clean(out)


# ---------------------------------------------------------------------------
# Lyapunov weights W(S) and Z(T)
# ---------------------------------------------------------------------------
def W_tilde_S(S: TypeS) -> Fraction:
    _, _, ch = S
    total = h_val(1 + len(ch))
    for Y in ch:
        total += z_val(Y[0]) * W_tilde_S(Y)
    return total


def W_S(S: TypeS) -> Fraction:
    return W_tilde_S(S) / (S[0] + D_PARAM)


def Z_T(T: TreeT) -> Fraction:
    _, ch = T
    total = h_val(len(ch))
    for Y in ch:
        total += z_val(Y[0]) * W_tilde_S(Y)
    return total


def eval_lincomb(lc: LinComb, fn) -> Fraction:
    return sum((coeff * fn(S) for S, coeff in lc.items()), Fraction(0))


def apply_op_to_lincomb(op_basis, lc: LinComb) -> LinComb:
    out: LinComb = {}
    for S, coeff in lc.items():
        for S_next, c2 in op_basis(S).items():
            lin_add(out, S_next, coeff * c2)
    return out


def apply_Lphys_to_tree_comb(tcomb: Dict[TreeT, Fraction]) -> Dict[TreeT, Fraction]:
    out: Dict[TreeT, Fraction] = {}
    for T, coeff in tcomb.items():
        trans = canonical_physical_transitions_T(T)
        total_rate = Fraction(0)
        for T_next, rate in trans.items():
            v = out.get(T_next, Fraction(0)) + coeff * rate
            if v == 0:
                out.pop(T_next, None)
            else:
                out[T_next] = v
            total_rate += rate
        v0 = out.get(T, Fraction(0)) - coeff * total_rate
        if v0 == 0:
            out.pop(T, None)
        else:
            out[T] = v0
    return out


def check_mark_ineq_S(S: TypeS) -> None:
    p_S = Fraction(1, deg_S(S))
    for Y in S[2]:
        assert Y[0] >= p_S, f"Internal mark inequality failed: p_S={p_S}, a_Y={Y[0]}"
        assert Y[0] <= Fraction(1, 2), f"Mark > 1/2: {Y[0]}"
        check_mark_ineq_S(Y)


def check_mark_ineq_T(T: TreeT) -> None:
    p_o = Fraction(1, T[0] + len(T[1]))
    for Y in T[1]:
        assert Y[0] >= p_o, f"Root mark inequality failed: p_o={p_o}, a_Y={Y[0]}"
        assert Y[0] <= Fraction(1, 2), f"Mark > 1/2: {Y[0]}"
        check_mark_ineq_S(Y)


# ---------------------------------------------------------------------------
# Test 1: Lyapunov rational inequality (4)
# ---------------------------------------------------------------------------
def test_inequality_4() -> None:
    print("=== Test 1: Exact rational check of Lyapunov inequality (4) ===")
    worst_diff = Fraction(-100)
    worst_pair = None
    count = 0
    for D in range(2, 201):
        for c in range(0, D + 1):
            lhs = (1 - Fraction(c, D)) * (h_val(c + 1) - h_val(c) + Fraction(1, 1 + D_PARAM * D))
            rhs = Q_PARAM * h_val(c)
            diff = lhs - rhs
            assert diff <= 0, f"Inequality (4) FAILED at c={c}, D={D}: lhs={lhs}, rhs={rhs}"
            if diff > worst_diff:
                worst_diff = diff
                worst_pair = (c, D, lhs, rhs)
            count += 1
    assert 1 + D_PARAM * 100 == Fraction(16)
    assert Fraction(1, 16) < Q_PARAM * Fraction(7, 10)
    assert RHO_PARAM == Fraction(51, 100)
    assert Fraction(7, 10) / (Fraction(1, 2) + D_PARAM) == Fraction(14, 13) > 1
    c, D, lhs, rhs = worst_pair
    print(f"  Checked {count} pairs (D=2..200, 0<=c<=D): 0 failures.")
    print(f"  Tightest pair: c={c}, D={D}: lhs={lhs} ({float(lhs):.6f}) <= rhs={rhs} ({float(rhs):.6f}), margin={-worst_diff} ({float(-worst_diff):.6f})")
    print(f"  D>=100 tail bound: 1/16 = {float(Fraction(1,16)):.4f} < q*7/10 = {Q_PARAM*Fraction(7,10)} ({float(Q_PARAM*Fraction(7,10)):.4f}). PASS")


# ---------------------------------------------------------------------------
# Test 2a: Exhaustive individual directed-edge labelled growth histories (n=3..7)
# ---------------------------------------------------------------------------
def test_exhaustive_individual_histories(max_n: int = 7) -> Dict[int, Dict[TreeT, Fraction]]:
    print(f"=== Test 2a: Exhaustive individual directed-edge labelled growth histories (n=3..{max_n}) ===")
    # Each history is a sequence of directed edges (w -> z) chosen at each step.
    # For n = 3..7, number of directed-edge histories is 1, 4, 24, 192, 1920.
    histories: List[Tuple[LabelledTree, Fraction]] = [(T3_labelled(), Fraction(1))]
    hist_dists: Dict[int, Dict[TreeT, Fraction]] = {}
    for n in range(3, max_n + 1):
        dist: Dict[TreeT, Fraction] = {}
        for T_lab, prob in histories:
            T_can = T_lab.canonical_T()
            check_mark_ineq_T(T_can)
            dist[T_can] = dist.get(T_can, Fraction(0)) + prob
        assert sum(dist.values()) == Fraction(1)
        hist_dists[n] = dist
        print(f"  n={n}: {len(histories)} individual directed-edge histories -> {len(dist)} canonical marked trees, sum(p)=1.")
        if n < max_n:
            nxt_hist: List[Tuple[LabelledTree, Fraction]] = []
            for T_lab, prob in histories:
                for rate, T_next in T_lab.physical_transitions():
                    nxt_hist.append((T_next, prob * rate / Fraction(n)))
            histories = nxt_hist
    print("  PASS")
    return hist_dists


# ---------------------------------------------------------------------------
# Test 2b: Exhaustive exact reachable marked trees (n = 3..9) & transfer identity
# ---------------------------------------------------------------------------
def enumerate_reachable_trees(max_n: int) -> Tuple[Dict[int, Dict[TreeT, LabelledTree]], Set[TypeS]]:
    by_n: Dict[int, Dict[TreeT, LabelledTree]] = {3: {T3_canonical(): T3_labelled()}}
    all_types: Set[TypeS] = set()

    for n in range(3, max_n):
        nxt: Dict[TreeT, LabelledTree] = {}
        for T_can, T_lab in by_n[n].items():
            bf_rates: Dict[TreeT, Fraction] = {}
            for rate, T_lab_next in T_lab.physical_transitions():
                T_can_next = T_lab_next.canonical_T()
                bf_rates[T_can_next] = bf_rates.get(T_can_next, Fraction(0)) + rate
                if T_can_next not in nxt:
                    nxt[T_can_next] = T_lab_next
            can_rates = canonical_physical_transitions_T(T_can)
            assert bf_rates == can_rates, f"Mismatch bf vs canonical at T={T_can}"
            assert sum(bf_rates.values()) == Fraction(n), "Total physical rate != n"
        by_n[n + 1] = nxt

    for n in range(3, max_n + 1):
        for T_can in by_n[n]:
            check_mark_ineq_T(T_can)
            for S in all_noninit_subtrees_T(T_can):
                all_types.add(S)
                all_types.add(add_leaf_root(S))

    return by_n, all_types


def test_exhaustive_transfer_and_weights(by_n: Dict[int, Dict[TreeT, LabelledTree]], all_types: Set[TypeS]) -> None:
    max_n = max(by_n)
    print(f"=== Test 2b: Exhaustive brute-force transfer identity (n = 3..{max_n}) ===")
    total_trees = 0
    for n, trees in sorted(by_n.items()):
        for T_can, T_lab in trees.items():
            total_trees += 1
            old_subtrees = T_lab.noninit_subtrees()
            assert sorted(old_subtrees) == sorted(all_noninit_subtrees_T(T_can))
            lhs: LinComb = {}
            for rate, T_lab_next in T_lab.physical_transitions():
                new_subtrees = T_lab_next.noninit_subtrees()
                for S in new_subtrees:
                    lin_add(lhs, S, rate)
                for S in old_subtrees:
                    lin_add(lhs, S, -rate)

            rhs: LinComb = {}
            for S in old_subtrees:
                for S_out, coeff in apply_RB_basis(S).items():
                    lin_add(rhs, S_out, coeff)
            for S_out, coeff in apply_I_basis(T_can).items():
                lin_add(rhs, S_out, coeff)

            assert lhs == rhs, f"Transfer identity FAILED at n={n}, T={T_can}"
        print(f"  n={n}: {len(trees)} distinct marked rooted trees checked: 0 failures.")
    print(f"  Total marked rooted trees verified (exact formal basis): {total_trees}. PASS")

    print("=== Test 3: Lyapunov identities & inequalities on all reachable types (+ a=0) ===")
    aux_types: Set[TypeS] = set(all_types)
    aux_types.add(J_type(Fraction(0)))
    for S in list(all_types):
        aux_types.add((Fraction(0), S[1], S[2]))

    max_l0_ratio = Fraction(-100)
    max_r_ratio = Fraction(-100)
    min_W = Fraction(1000)
    max_W_over_size = Fraction(0)
    max_BW = Fraction(0)

    for S in aux_types:
        check_mark_ineq_S(S)
        wS = W_S(S)
        assert wS >= Fraction(14, 13) > 1
        assert wS <= Fraction(size_S(S)) / D_PARAM
        if wS < min_W:
            min_W = wS
        r_ws = wS / size_S(S)
        if r_ws > max_W_over_size:
            max_W_over_size = r_ws

        child_sum = sum((Y[0] * W_S(Y) for Y in S[2]), Fraction(0))
        expected_child_sum = (S[0] + D_PARAM) * wS - h_val(core_deg_S(S))
        assert child_sum == expected_child_sum
        assert W_S(add_leaf_root(S)) == wS

        l0_w = eval_lincomb(apply_L0_basis(S), W_S)
        assert l0_w <= Q_PARAM * wS
        if l0_w / wS > max_l0_ratio:
            max_l0_ratio = l0_w / wS

        r_w = eval_lincomb(apply_R_basis(S), W_S)
        assert r_w == l0_w + D_PARAM * wS - h_val(core_deg_S(S))
        assert r_w <= RHO_PARAM * wS - Fraction(7, 10)
        if r_w / wS > max_r_ratio:
            max_r_ratio = r_w / wS

        b_w = eval_lincomb(apply_B_basis(S), W_S)
        assert b_w <= Fraction(1) / D_PARAM
        if b_w > max_BW:
            max_BW = b_w

    max_z_ratio = Fraction(-100)
    for n, trees in sorted(by_n.items()):
        for T_can in trees:
            zT = Z_T(T_can)
            trans = canonical_physical_transitions_T(T_can)
            lz = sum((rate * (Z_T(T_next) - zT) for T_next, rate in trans.items()), Fraction(0))
            assert lz <= Q_PARAM * zT
            if lz / zT > max_z_ratio:
                max_z_ratio = lz / zT

            L_o, ch = T_can
            p_o = Fraction(1, L_o + len(ch))
            imm_ch_w = sum(((Y[0] - p_o) * W_S(Y) + p_o * W_S(add_leaf_root(Y)) for Y in ch), Fraction(0))
            assert imm_ch_w == zT - h_val(len(ch))

            imm_w = eval_lincomb(apply_I_basis(T_can), W_S)
            assert imm_w <= Fraction(1) / D_PARAM + zT - h_val(len(ch))

    print(f"  Checked {len(aux_types)} types (reachable + a=0 auxiliary) and {total_trees} full trees: 0 failures.")
    print(f"  min W = {min_W} ({float(min_W):.4f} >= 14/13), max W/|S| = {max_W_over_size} (<= 1/d={Fraction(1)/D_PARAM}), max BW = {max_BW} (<= 1/d)")
    print(f"  max (L^0 W)/W = {float(max_l0_ratio):.6f} (<= q={float(Q_PARAM)}), max (RW)/W = {float(max_r_ratio):.6f} (<= rho={float(RHO_PARAM)}), max (L^phys Z)/Z = {float(max_z_ratio):.6f} (<= q)")
    print("  PASS")


# ---------------------------------------------------------------------------
# Test 4: Exact multistep Taylor check of eq (3), eq (8), and Volterra eq
# ---------------------------------------------------------------------------
def test_multistep_identities() -> None:
    print("=== Test 4: Exact multistep Taylor check of eq (3), eq (8), and Volterra eq (orders k=1..5) ===")
    max_order = 5
    mu_list: List[Dict[TreeT, Fraction]] = [{T3_canonical(): Fraction(1)}]
    for _ in range(max_order):
        mu_list.append(apply_Lphys_to_tree_comb(mu_list[-1]))

    imm_list: List[LinComb] = []
    for j in range(max_order):
        lc: LinComb = {}
        for T, coeff in mu_list[j].items():
            for S, c2 in apply_I_basis(T).items():
                lin_add(lc, S, coeff * c2)
        imm_list.append(lc)

    for k in range(1, max_order + 1):
        lhs: LinComb = {}
        for T, coeff in mu_list[k].items():
            for S in all_noninit_subtrees_T(T):
                lin_add(lhs, S, coeff)

        rhs: LinComb = {}
        for j in range(k):
            term = dict(imm_list[j])
            for _ in range(k - 1 - j):
                term = apply_op_to_lincomb(apply_RB_basis, term)
            for S, coeff in term.items():
                lin_add(rhs, S, coeff)

        assert lhs == rhs, f"Multistep eq (3) Taylor order k={k} FAILED!"
        print(f"  eq (3) order k={k}: {len(lhs)} non-zero basis types match identically.")

    # Copy-separation eq (8) AND Volterra decomposition P_t = R_t + int_0^t R_s B P_{t-s} ds
    test_marks = [Fraction(1, 2), Fraction(1, 3), Fraction(1, 10), Fraction(0)]

    def strip_root_mark(S: TypeS) -> Tuple[int, Tuple[TypeS, ...]]:
        return (S[1], S[2])

    ref_p0_stripped = None
    ref_gamma_stripped = None

    for a in test_marks:
        S0 = J_type(a)
        R_pow: List[LinComb] = [{S0: Fraction(1)}]
        RB_pow: List[LinComb] = [{S0: Fraction(1)}]
        K_pow: List[LinComb] = [{S0: Fraction(1)}]
        P0_pow: List[LinComb] = [{S0: Fraction(1)}]

        def op_L0_minus_a(S: TypeS, _a: Fraction = a) -> LinComb:
            res = apply_L0_basis(S)
            lin_add(res, S, -_a)
            return res

        for k in range(1, max_order + 1):
            R_pow.append(apply_op_to_lincomb(apply_R_basis, R_pow[-1]))
            RB_pow.append(apply_op_to_lincomb(apply_RB_basis, RB_pow[-1]))
            K_pow.append(apply_op_to_lincomb(op_L0_minus_a, K_pow[-1]))
            P0_pow.append(apply_op_to_lincomb(apply_L0_basis, P0_pow[-1]))

            # eq (8): R^k = (L^0 - a)^k + sum_{j=0}^{k-1} R^{k-1-j} Gamma (L^0 - a)^j
            rhs_k: LinComb = dict(K_pow[k])
            for j in range(k):
                step = apply_op_to_lincomb(apply_Gamma_basis, K_pow[j])
                for _ in range(k - 1 - j):
                    step = apply_op_to_lincomb(apply_R_basis, step)
                for S, coeff in step.items():
                    lin_add(rhs_k, S, coeff)
            assert R_pow[k] == rhs_k, f"eq (8) failed at a={a}, k={k}"

            # Volterra: (R+B)^k = R^k + sum_{j=0}^{k-1} (R+B)^{k-1-j} B R^j (in measure action)
            vol_k: LinComb = dict(R_pow[k])
            for j in range(k):
                step = apply_op_to_lincomb(apply_B_basis, R_pow[j])
                for _ in range(k - 1 - j):
                    step = apply_op_to_lincomb(apply_RB_basis, step)
                for S, coeff in step.items():
                    lin_add(vol_k, S, coeff)
            assert RB_pow[k] == vol_k, f"Volterra generation eq failed at a={a}, k={k}"

        p0_stripped = [
            {strip_root_mark(S): c for S, c in P0_pow[k].items()}
            for k in range(max_order + 1)
        ]
        gamma_p0 = [
            apply_op_to_lincomb(apply_Gamma_basis, P0_pow[k])
            for k in range(max_order + 1)
        ]
        if ref_p0_stripped is None:
            ref_p0_stripped = p0_stripped
            ref_gamma_stripped = gamma_p0
        else:
            assert p0_stripped == ref_p0_stripped
            assert gamma_p0 == ref_gamma_stripped

    print(f"  eq (8) copy-separation, exact a-independence, and Volterra P_t = R_t + int R B P verified for a in {test_marks}, orders k=1..{max_order}. PASS")


# ---------------------------------------------------------------------------
# Test 5: Count comparisons & discrete de-Poissonization bound
# ---------------------------------------------------------------------------
def check_tree_count_comparisons(T_lab: LabelledTree, label: str = "") -> None:
    deg = T_lab.degrees()
    ch = T_lab.children()
    cores = [v for v in range(T_lab.n) if deg[v] >= 2]
    C = len(cores)
    assert C >= 1

    if all(T_lab.birth_mark[v] is not None for v in cores if v != 0):
        subtrees = T_lab.noninit_subtrees()
        sum_W = sum((W_S(S) for S in subtrees), Fraction(0))
        bound_W = Fraction(260, 9) * (C - 1)
        assert sum_W <= bound_W, f"Upper W bound failed on {label}: {sum_W} > {bound_W}"

    core_set = set(cores)
    adj = T_lab.adj()
    core_nbrs = {v: [u for u in adj[v] if u in core_set] for v in cores}
    leaf_cnt = {v: deg[v] - len(core_nbrs[v]) for v in cores}
    core_ch = {v: [u for u in ch[v] if u in core_set] for v in cores}

    def is_bare(v: int) -> bool:
        return len(core_nbrs[v]) == 2 and leaf_cnt[v] == 0

    eligible = {
        v for v in cores
        if (len(core_nbrs[v]) <= 2 and leaf_cnt[v] >= 1)
        or (is_bare(v) and any(is_bare(u) for u in core_nbrs[v]))
    }
    assert Fraction(len(eligible)) >= Fraction(C, 4), f"Eligibility < C/4 on {label}"

    G1 = {v for v in cores if len(core_nbrs[v]) <= 2 and leaf_cnt[v] >= 1}
    G2 = {v for v in cores if is_bare(v) and any(is_bare(u) for u in core_ch[v])}
    assert G1.isdisjoint(G2)
    B_rem = eligible - (G1 | G2)

    assert 0 not in B_rem, f"Root o in B_rem on {label}"
    fibres: Dict[int, List[int]] = defaultdict(list)
    for u in B_rem:
        p = T_lab.parent[u]
        assert p in G2, f"Parent of B_rem vertex not in G2 on {label}"
        fibres[p].append(u)
    for w, children_in_B in fibres.items():
        if w == 0:
            assert len(children_in_B) <= 2
        else:
            assert len(children_in_B) <= 1

    good_noninit = (G1 | G2) - {0}
    assert Fraction(len(good_noninit)) >= Fraction(C, 12) - 1


def make_adversarial_trees() -> List[Tuple[str, LabelledTree]]:
    trees: List[Tuple[str, LabelledTree]] = []

    parent = (-1, 0, 0, 1, 2, 3, 4)
    marks = (None, Fraction(1, 2), Fraction(1, 2), Fraction(1, 2), Fraction(1, 2), None, None)
    trees.append(("bare_root_2_bare_children_C5", LabelledTree(7, parent, marks)))

    par = [-1]
    mk: List[Optional[Fraction]] = [None]
    prev = 0
    for _ in range(10):
        idx = len(par); par.append(prev); mk.append(Fraction(1, 2)); prev = idx
    par.append(prev); mk.append(None)
    prev = 0
    for _ in range(10):
        idx = len(par); par.append(prev); mk.append(Fraction(1, 2)); prev = idx
    par.append(prev); mk.append(None)
    trees.append(("long_bare_chain_C21", LabelledTree(len(par), tuple(par), tuple(mk))))

    par = [-1]
    mk = [None]
    queue = [(0, 0)]
    while queue:
        u, dep = queue.pop(0)
        if dep < 4:
            for _ in range(2):
                b1 = len(par); par.append(u); mk.append(Fraction(1, 2))
                b2 = len(par); par.append(b1); mk.append(Fraction(1, 2))
                nxt = len(par); par.append(b2); mk.append(Fraction(1, 2))
                queue.append((nxt, dep + 1))
        else:
            par.append(u); mk.append(None)
    trees.append(("subdivided_binary_tree", LabelledTree(len(par), tuple(par), tuple(mk))))

    par = [-1, 0]
    mk = [None, None]
    prev = 0
    for _ in range(30):
        idx = len(par); par.append(prev); mk.append(Fraction(1, 2)); prev = idx
    par.append(prev); mk.append(None)
    trees.append(("deep_path_C31", LabelledTree(len(par), tuple(par), tuple(mk))))

    return trees


def test_count_comparisons_and_discrete_bound(
    by_n: Dict[int, Dict[TreeT, LabelledTree]],
    hist_dists: Dict[int, Dict[TreeT, Fraction]],
) -> None:
    print("=== Test 5: Count comparisons & discrete de-Poissonization bound ===")
    max_n = max(by_n)
    for n, trees in sorted(by_n.items()):
        for T_can, T_lab in trees.items():
            check_tree_count_comparisons(T_lab, f"reachable_n{n}")

    for name, T_lab in make_adversarial_trees():
        check_tree_count_comparisons(T_lab, name)
    print(f"  Upper W bound & lower eligibility/good-vertex bounds verified on all reachable n<={max_n} trees and 4 adversarial trees.")

    dist: Dict[TreeT, Fraction] = {T3_canonical(): Fraction(1)}
    a_vals: Dict[int, Fraction] = {3: Fraction(1)}
    for n in range(3, max_n):
        if n in hist_dists:
            assert dist == hist_dists[n], f"Markov chain dist != exhaustive history dist at n={n}"
        nxt_dist: Dict[TreeT, Fraction] = {}
        for T, prob in dist.items():
            C_T = core_count_T(T)
            trans = canonical_physical_transitions_T(T)
            prom_prob = Fraction(0)
            for T_next, rate in trans.items():
                step_p = rate / Fraction(n)
                nxt_dist[T_next] = nxt_dist.get(T_next, Fraction(0)) + prob * step_p
                if core_count_T(T_next) == C_T + 1:
                    prom_prob += step_p
            assert prom_prob <= Fraction(C_T, n)
        assert sum(nxt_dist.values()) == Fraction(1)
        a_next = sum((p * core_count_T(T) for T, p in nxt_dist.items()), Fraction(0))
        assert a_vals[n] <= a_next <= a_vals[n] * Fraction(n + 1, n)
        a_vals[n + 1] = a_next
        dist = nxt_dist
    if max_n in hist_dists:
        assert dist == hist_dists[max_n]

    print("  Exact a_n = E[C_n]:", ", ".join(f"a_{n}={float(a_vals[n]):.5f}" for n in range(3, max_n + 1)))
    print(f"  Verified Markov dist == exhaustive individual-history dist (n<=7), a_n <= a_{{n+1}} <= ((n+1)/n) a_n, and prom_prob <= C_n/n (n=3..{max_n-1}). PASS")


# ---------------------------------------------------------------------------
# Test 6: Deliberate corruptions (negative controls)
# ---------------------------------------------------------------------------
def test_deliberate_corruptions(by_n: Dict[int, Dict[TreeT, LabelledTree]]) -> None:
    print("=== Test 6: Deliberate corruptions (all must be caught) ===")

    def find_transfer_failure(rb_fn, imm_fn, use_current_parent_mark: bool = False) -> Optional[Tuple[int, TreeT]]:
        for n in sorted(by_n):
            for T_can, T_lab in by_n[n].items():
                if use_current_parent_mark:
                    def extract_cur(tl: LabelledTree) -> List[TypeS]:
                        deg = tl.degrees()
                        ch = tl.children()
                        def rec(v: int) -> TypeS:
                            a_cur = Fraction(1, deg[tl.parent[v]])
                            leaves = sum(1 for c in ch[v] if deg[c] == 1)
                            core_ch = tuple(sorted(rec(c) for c in ch[v] if deg[c] >= 2))
                            return (a_cur, leaves, core_ch)
                        return [rec(v) for v in range(1, tl.n) if deg[v] >= 2]
                    old_sub = extract_cur(T_lab)
                    lhs: LinComb = {}
                    for rate, T_next in T_lab.physical_transitions():
                        for S in extract_cur(T_next):
                            lin_add(lhs, S, rate)
                        for S in old_sub:
                            lin_add(lhs, S, -rate)
                    deg = T_lab.degrees(); ch = T_lab.children()
                    def rec_t(v: int) -> TypeS:
                        a_cur = Fraction(1, deg[T_lab.parent[v]])
                        leaves = sum(1 for c in ch[v] if deg[c] == 1)
                        core_ch = tuple(sorted(rec_t(c) for c in ch[v] if deg[c] >= 2))
                        return (a_cur, leaves, core_ch)
                    T_cur = (
                        sum(1 for c in ch[0] if deg[c] == 1),
                        tuple(sorted(rec_t(c) for c in ch[0] if deg[c] >= 2)),
                    )
                    rhs: LinComb = {}
                    for S in old_sub:
                        for S_out, c in rb_fn(S).items():
                            lin_add(rhs, S_out, c)
                    for S_out, c in imm_fn(T_cur).items():
                        lin_add(rhs, S_out, c)
                    if lhs != rhs:
                        return (n, T_can)
                else:
                    old_sub = T_lab.noninit_subtrees()
                    lhs = {}
                    for rate, T_next in T_lab.physical_transitions():
                        for S in T_next.noninit_subtrees():
                            lin_add(lhs, S, rate)
                        for S in old_sub:
                            lin_add(lhs, S, -rate)
                    rhs = {}
                    for S in old_sub:
                        for S_out, c in rb_fn(S).items():
                            lin_add(rhs, S_out, c)
                    for S_out, c in imm_fn(T_can).items():
                        lin_add(rhs, S_out, c)
                    if lhs != rhs:
                        return (n, T_can)
        return None

    # Corruption 1: Non-permanent current 1/D_parent mark
    w1 = find_transfer_failure(apply_RB_basis, apply_I_basis, use_current_parent_mark=True)
    assert w1 is not None
    print(f"  [CAUGHT] Corruption 1 (non-permanent current 1/D_parent mark): fails at n={w1[0]}")

    # Corruption 2: Drop (a_Y - p_S) f(Y) copy term in R
    def rb_no_copy_Y(S: TypeS) -> LinComb:
        res = apply_L0_basis(S)
        lin_add(res, S, -S[0])
        p_S = Fraction(1, deg_S(S))
        for Y in S[2]:
            lin_add(res, add_leaf_root(Y), p_S)
        for k, v in apply_B_basis(S).items():
            lin_add(res, k, v)
        return res

    w2 = find_transfer_failure(rb_no_copy_Y, apply_I_basis)
    assert w2 is not None
    print(f"  [CAUGHT] Corruption 2 (drop (a_Y-p_S)f(Y) copy in R): fails at n={w2[0]}")

    # Corruption 3: Replace Y^+ with Y in R copy term
    def rb_no_plus(S: TypeS) -> LinComb:
        res = apply_L0_basis(S)
        lin_add(res, S, -S[0])
        for Y in S[2]:
            lin_add(res, Y, Y[0])
        for k, v in apply_B_basis(S).items():
            lin_add(res, k, v)
        return res

    w3 = find_transfer_failure(rb_no_plus, apply_I_basis)
    assert w3 is not None
    print(f"  [CAUGHT] Corruption 3 (replace Y^+ by Y in R): fails at n={w3[0]}")

    # Corruption 4: Drop (a_Y - p_o) f(Y) from initial-centre immigration I
    def imm_no_copy_Y(T: TreeT) -> LinComb:
        L_o, ch = T
        D_o = L_o + len(ch)
        p_o = Fraction(1, D_o)
        res: LinComb = {}
        if L_o > 0:
            lin_add(res, J_type(p_o), Fraction(L_o, D_o))
        for Y in ch:
            lin_add(res, add_leaf_root(Y), p_o)
        return res

    w4 = find_transfer_failure(apply_RB_basis, imm_no_copy_Y)
    assert w4 is not None
    print(f"  [CAUGHT] Corruption 4 (drop (a_Y-p_o)f(Y) in immigration I): fails at n={w4[0]}")

    # Corruption 5: Omit direct birth B
    w5 = find_transfer_failure(apply_R_basis, apply_I_basis)
    assert w5 is not None
    print(f"  [CAUGHT] Corruption 5 (omit direct birth B): fails at n={w5[0]}")

    # Corruption 6: Omit parent slot in D_o
    def rb_no_parent_slot(S: TypeS) -> LinComb:
        res = apply_R_basis(S)
        _, L_o, ch = S
        D_bad = L_o + len(ch)
        if L_o > 0 and D_bad > 0:
            lin_add(res, J_type(Fraction(1, D_bad)), Fraction(L_o, D_bad))
        return res

    w6 = find_transfer_failure(rb_no_parent_slot, apply_I_basis)
    assert w6 is not None
    print(f"  [CAUGHT] Corruption 6 (omit parent slot in D_o): fails at n={w6[0]}")

    # Corruption 7: Replace killing rate a_S by constant 1/2
    def rb_const_kill(S: TypeS) -> LinComb:
        res = apply_L0_basis(S)
        lin_add(res, S, -Fraction(1, 2))
        for k, v in apply_Gamma_basis(S).items():
            lin_add(res, k, v)
        for k, v in apply_B_basis(S).items():
            lin_add(res, k, v)
        return res

    w7 = find_transfer_failure(rb_const_kill, apply_I_basis)
    assert w7 is not None
    print(f"  [CAUGHT] Corruption 7 (constant killing rate 1/2 instead of a_S): fails at n={w7[0]}")

    # Corruption 8: Tighten q = 9/25 to q' = 359/1000
    q_bad = Fraction(359, 1000)
    caught_q = None
    for D in range(2, 50):
        for c in range(0, D + 1):
            lhs = (1 - Fraction(c, D)) * (h_val(c + 1) - h_val(c) + Fraction(1, 1 + D_PARAM * D))
            if lhs > q_bad * h_val(c):
                caught_q = (c, D, lhs, q_bad * h_val(c))
                break
        if caught_q:
            break
    assert caught_q is not None
    print(f"  [CAUGHT] Corruption 8 (tighten q to 359/1000): fails at c={caught_q[0]}, D={caught_q[1]} ({caught_q[2]} > {caught_q[3]})")

    # Corruption 9: Perturb h(2) from 17/20 to 16/20
    def h_bad(c: int) -> Fraction:
        return Fraction(16, 20) if c == 2 else h_val(c)

    caught_h = None
    for D in range(2, 50):
        for c in range(0, D + 1):
            lhs = (1 - Fraction(c, D)) * (h_bad(c + 1) - h_bad(c) + Fraction(1, 1 + D_PARAM * D))
            if lhs > Q_PARAM * h_bad(c):
                caught_h = (c, D, lhs, Q_PARAM * h_bad(c))
                break
        if caught_h:
            break
    assert caught_h is not None
    print(f"  [CAUGHT] Corruption 9 (perturb h(2) to 16/20): fails at c={caught_h[0]}, D={caught_h[1]} ({caught_h[2]} > {caught_h[3]})")

    # Verify Good kind 2 root leaf-gain rate in L^0
    S_bare_with_bare_child: TypeS = (Fraction(1, 2), 0, ((Fraction(1, 2), 0, (J_type(Fraction(1, 2)),)),))
    trans_good = intrinsic_transitions_S(S_bare_with_bare_child)
    assert trans_good.get(add_leaf_root(S_bare_with_bare_child), Fraction(0)) == Fraction(1, 2)
    print("  [VERIFIED] Good kind 2 (bare vertex with bare child): root leaf-gain rate in L^0 is exactly 1/2.")
    print("  All 9 deliberate corruptions caught + Good kind 2 verified. PASS")


def main() -> None:
    t0 = time.time()
    test_inequality_4()
    hist_dists = test_exhaustive_individual_histories(max_n=7)
    by_n, all_types = enumerate_reachable_trees(max_n=9)
    test_exhaustive_transfer_and_weights(by_n, all_types)
    test_multistep_identities()
    test_count_comparisons_and_discrete_bound(by_n, hist_dists)
    test_deliberate_corruptions(by_n)
    dt = time.time() - t0
    print(f"\nALL INDEPENDENT CHECKS PASSED in {dt:.2f}s.")


if __name__ == "__main__":
    main()
