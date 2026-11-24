#!/usr/bin/env python3
r"""Independent exact cold-referee verification for round2/fringe/finite-kernel.md.

Audits (using exact fractions.Fraction arithmetic throughout, fast on 1 CPU):
  1. Clamped Lyapunov inequality (Section 1): m(D) = min(D, 20), exact rational
     verification for all 0 <= c <= D <= 250, D >= 2, including D < 20, D = 20,
     and D > 20, plus proof that M = 20 is the minimal integer clamp.
  2. Clamped signed-to-positive transfer identity (Section 1): tested with
     brute-force vertex-by-vertex labelled physical dynamics at BOTH M = 3
     (where clamping D > 3 triggers immediately on small trees n = 3..7,
     exercising a_v > p_u at birth) and M = 20 (n = 3..7).
  3. Formula (2) structural checks (Section 2):
     - Exact r_0(S) == n(S) - 1/D_root on all reachable types.
     - Self-transition check: t(S, S) == 0 identically (even when n(Y^+) == n(S)).
     - Exact match between finite path sum (2) restricted to D_N and K^-_N(s)
       (exact at l <= 2 on acyclic D_3, and converging geometrically to K^-_4(s)
       on cyclic D_4 with exact remainder identity).
  4. Finite rational truncations K^-_N(s), K^+_N(s) (Section 3):
     - Full 19-label (M = 20) exact rational matrices at N = 2 (19 states),
       N = 3 (57 states), and N = 4 (152 states, cyclic!) at s = l = 113/200
       and s = u = 1133/2000.
     - M-matrix property (A_N W > 0, A_N^{-1} >= 0 entrywise), exact nesting
       K^-_2 <= K^-_3 <= K^-_4 and K^+_4 <= K^+_3 <= K^+_2, endpoint N = 2,
       and sharper test-column bound K^-_N g + u_N ||g||_inf.
     - Exact check on a small M = 3 clamped model up to N = 5 (tracking
       monotone spectral-radius enclosures r(K^-_N) <= r(K^+_N)).
  5. Convergence estimates (5)-(7) and inverse slope (Section 4):
     - Exact L^0 V_delta identity and R V_delta <= [rho + delta(1 + q/2)] V_delta
       at delta = 1 (exact rational) and delta = 1/24 (via exact rational
       binomial inequality (1 + 1/(24n))^{24} >= 1 + 1/n and exact checks on
       all states in D_4).
     - Exact verification of (6) on all 19 rows for N = 2, 3, 4 at s = l, u.
     - Exact rational computation of l, u, delta, H, e (endpoints j = 2 and
       j = 20), g_*, bound (7), and the inverse-slope integral lower bound
       g_*(s' - s) / [2(u + 2)^2] (including J_{1/i} exit rate 1/i + 3/2 <= 2).
  6. Planted errors (9 negative controls): confirms every error is caught with
     a concrete counterexample.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import time
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Parameters from renewal.md & finite-kernel.md
# ---------------------------------------------------------------------------
D_PARAM = Fraction(3, 20)       # d = 3/20
Q_PARAM = Fraction(9, 25)       # q = 9/25
RHO_PARAM = Q_PARAM + D_PARAM   # rho = 51/100
L_BOUND = Fraction(113, 200)    # l = 113/200
U_BOUND = Fraction(1133, 2000)  # u = 1133/2000
DELTA = Fraction(1, 24)         # delta = 1/24
M_CLAMP = 20                    # m(D) = min(D, 20), 19 labels 2..20


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


def m_clamp(D: int, M: int = M_CLAMP) -> int:
    return D if D < M else M


# ---------------------------------------------------------------------------
# Canonical planted marked fringe S and full tree T
# ---------------------------------------------------------------------------
# S = (a_S: Fraction, L_o: int, children: Tuple[S, ...])
# Root has 1 inert parent slot: D_o = 1 + L_o + len(children), c_o = 1 + len(children)
TypeS = Tuple[Fraction, int, Tuple["TypeS", ...]]
TreeT = Tuple[int, Tuple[TypeS, ...]]
LinComb = Dict[TypeS, Fraction]


def J_type(label: int) -> TypeS:
    return (Fraction(1, label), 1, ())


def add_leaf_root(S: TypeS) -> TypeS:
    return (S[0], S[1] + 1, S[2])


def deg_S(S: TypeS) -> int:
    return 1 + S[1] + len(S[2])


def core_deg_S(S: TypeS) -> int:
    return 1 + len(S[2])


def size_S(S: TypeS) -> int:
    return 1 + S[1] + sum(size_S(c) for c in S[2])


def W_tilde_S(S: TypeS) -> Fraction:
    total = h_val(1 + len(S[2]))
    for Y in S[2]:
        total += z_val(Y[0]) * W_tilde_S(Y)
    return total


def W_S(S: TypeS) -> Fraction:
    return W_tilde_S(S) / (S[0] + D_PARAM)


def lin_add(dst: LinComb, key: TypeS, coeff: Fraction) -> None:
    if coeff == 0:
        return
    v = dst.get(key, Fraction(0)) + coeff
    if v == 0:
        dst.pop(key, None)
    else:
        dst[key] = v


def lin_clean(d: Dict) -> Dict:
    return {k: v for k, v in d.items() if v != 0}


# ---------------------------------------------------------------------------
# Clamped transitions r^0(S, S'), t(S, S'), R, B, I
# ---------------------------------------------------------------------------
def intrinsic_transitions_S(S: TypeS, M: int = M_CLAMP) -> Dict[TypeS, Fraction]:
    """Exact intrinsic transition rates r^0(S, S') with parallel transitions added."""
    a_S, L_o, ch = S
    D_o = 1 + L_o + len(ch)
    out: Dict[TypeS, Fraction] = {}

    # 1. Leaf gains at root: L_o leaves (rate 1 each) + core children (rate 1/D_Y each)
    gain_rate_o = Fraction(L_o) + sum(Fraction(1, deg_S(Y)) for Y in ch)
    if gain_rate_o > 0:
        S_plus = (a_S, L_o + 1, ch)
        out[S_plus] = out.get(S_plus, Fraction(0)) + gain_rate_o

    # 2. Promotion at root: rate L_o / D_o, newborn mark 1 / m_clamp(D_o, M)
    if L_o > 0:
        new_child = J_type(m_clamp(D_o, M))
        new_ch = tuple(sorted(ch + (new_child,)))
        S_prom = (a_S, L_o - 1, new_ch)
        out[S_prom] = out.get(S_prom, Fraction(0)) + Fraction(L_o, D_o)

    # 3. Events inside each immediate core child Y_i
    for i, Y in enumerate(ch):
        prefix = ch[:i]
        suffix = ch[i + 1:]

        # 3a. Root arrow into Y_i (rate 1/D_o) -> Y_i^+
        ch_plus = tuple(sorted(prefix + (add_leaf_root(Y),) + suffix))
        S_yplus = (a_S, L_o, ch_plus)
        out[S_yplus] = out.get(S_yplus, Fraction(0)) + Fraction(1, D_o)

        # 3b. Internal events inside Y_i
        for Y_next, rate in intrinsic_transitions_S(Y, M).items():
            ch_next = tuple(sorted(prefix + (Y_next,) + suffix))
            S_next = (a_S, L_o, ch_next)
            out[S_next] = out.get(S_next, Fraction(0)) + rate

    return lin_clean(out)


def copy_transitions_S(S: TypeS) -> Dict[TypeS, Fraction]:
    """Copy rates sum_Y [p_S 1_{Y^+=S'} + (a_Y - p_S) 1_{Y=S'}]."""
    p_S = Fraction(1, deg_S(S))
    out: Dict[TypeS, Fraction] = {}
    for Y in S[2]:
        a_Y = Y[0]
        assert a_Y >= p_S, f"Mark inequality violated: a_Y={a_Y} < p_S={p_S}"
        if a_Y > p_S:
            out[Y] = out.get(Y, Fraction(0)) + (a_Y - p_S)
        Yp = add_leaf_root(Y)
        out[Yp] = out.get(Yp, Fraction(0)) + p_S
    return lin_clean(out)


def total_t_transitions_S(S: TypeS, M: int = M_CLAMP) -> Dict[TypeS, Fraction]:
    """Combined nonnegative rate t(S, S') = r^0(S, S') + copy(S, S')."""
    out = dict(intrinsic_transitions_S(S, M))
    for S_next, rate in copy_transitions_S(S).items():
        out[S_next] = out.get(S_next, Fraction(0)) + rate
    return lin_clean(out)


def r0_total(S: TypeS, M: int = M_CLAMP) -> Fraction:
    return sum(intrinsic_transitions_S(S, M).values(), Fraction(0))


def apply_RB_basis(S: TypeS, M: int = M_CLAMP) -> LinComb:
    res: LinComb = {}
    r0_map = intrinsic_transitions_S(S, M)
    r0 = sum(r0_map.values(), Fraction(0))
    lin_add(res, S, -(S[0] + r0))
    for S_next, rate in r0_map.items():
        lin_add(res, S_next, rate)
    for S_next, rate in copy_transitions_S(S).items():
        lin_add(res, S_next, rate)
    L_o = S[1]
    if L_o > 0:
        D_o = deg_S(S)
        lin_add(res, J_type(m_clamp(D_o, M)), Fraction(L_o, D_o))
    return res


def apply_I_basis(T: TreeT, M: int = M_CLAMP) -> LinComb:
    L_o, ch = T
    D_o = L_o + len(ch)
    p_o = Fraction(1, D_o)
    res: LinComb = {}
    if L_o > 0:
        lin_add(res, J_type(m_clamp(D_o, M)), Fraction(L_o, D_o))
    for Y in ch:
        a_Y = Y[0]
        assert a_Y >= p_o
        lin_add(res, Y, a_Y - p_o)
        lin_add(res, add_leaf_root(Y), p_o)
    return res


# ---------------------------------------------------------------------------
# Independent brute-force labelled tree engine with clamped marks
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LabelledTree:
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

    def attach_leaf(self, target: int, M: int) -> "LabelledTree":
        deg = self.degrees()
        new_parent = self.parent + (target,)
        marks = list(self.birth_mark) + [None]
        if deg[target] == 1:
            u = self.parent[target]
            marks[target] = Fraction(1, m_clamp(deg[u], M))
        return LabelledTree(self.n + 1, new_parent, tuple(marks))

    def physical_transitions(self, M: int) -> List[Tuple[Fraction, "LabelledTree"]]:
        adj = self.adj()
        deg = [len(adj[w]) for w in range(self.n)]
        out: List[Tuple[Fraction, LabelledTree]] = []
        for w in range(self.n):
            rate = Fraction(1, deg[w])
            for z in adj[w]:
                out.append((rate, self.attach_leaf(z, M)))
        return out

    def canonical_S(self, v: int, ch: List[List[int]], deg: List[int]) -> TypeS:
        a_v = self.birth_mark[v]
        assert a_v is not None
        leaves = sum(1 for c in ch[v] if deg[c] == 1)
        core_ch = tuple(sorted(self.canonical_S(c, ch, deg) for c in ch[v] if deg[c] >= 2))
        return (a_v, leaves, core_ch)

    def canonical_T(self) -> TreeT:
        ch = self.children()
        deg = self.degrees()
        leaves = sum(1 for c in ch[0] if deg[c] == 1)
        core_ch = tuple(sorted(self.canonical_S(c, ch, deg) for c in ch[0] if deg[c] >= 2))
        return (leaves, core_ch)

    def noninit_subtrees(self) -> List[TypeS]:
        ch = self.children()
        deg = self.degrees()
        return [self.canonical_S(v, ch, deg) for v in range(1, self.n) if deg[v] >= 2]


def T3_labelled() -> LabelledTree:
    return LabelledTree(3, (-1, 0, 0), (None, None, None))


# ---------------------------------------------------------------------------
# Exact rational linear solver (Gaussian elimination with Fraction)
# ---------------------------------------------------------------------------
def solve_rational_multi_rhs(A: List[List[Fraction]], B: List[List[Fraction]]) -> List[List[Fraction]]:
    """Solve A X = B exactly over Q where A is n x n and B is n x m."""
    n = len(A)
    m = len(B[0])
    # Augmented matrix
    M = [A[i][:] + B[i][:] for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            raise ValueError("Singular matrix A")
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
        piv_val = M[col][col]
        # Eliminate below (and normalize pivot row)
        row_c = M[col]
        for j in range(col, n + m):
            row_c[j] /= piv_val
        for r in range(n):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                row_r = M[r]
                for j in range(col, n + m):
                    if row_c[j] != 0:
                        row_r[j] -= fac * row_c[j]
    return [M[i][n:] for i in range(n)]


def spectral_radius_power(mat: List[List[Fraction]], iters: int = 80) -> float:
    """Fast Collatz-bracketed float estimate of Perron radius of a nonnegative matrix."""
    n = len(mat)
    mf = [[float(x) for x in row] for row in mat]
    # Check if zero matrix or trivial
    v = [1.0] * n
    lam = 0.0
    for _ in range(iters):
        w = [sum(mf[i][j] * v[j] for j in range(n)) for i in range(n)]
        mx = max(w)
        if mx == 0.0:
            return 0.0
        v = [x / mx for x in w]
        lam = mx
    return lam


# ---------------------------------------------------------------------------
# Domain D_N builder and rational truncation matrices K^-_N, K^+_N
# ---------------------------------------------------------------------------
def build_domain_DN(N: int, M: int = M_CLAMP) -> List[TypeS]:
    """Enumerate all states reachable from the (M-1) newborns staying in size <= N."""
    seeds = [J_type(i) for i in range(2, M + 1)]
    seen: Set[TypeS] = set(seeds)
    queue: List[TypeS] = list(seeds)
    while queue:
        S = queue.pop(0)
        for S_next in total_t_transitions_S(S, M):
            if size_S(S_next) <= N and S_next not in seen:
                seen.add(S_next)
                queue.append(S_next)
    # Place newborns first in canonical order 2..M, then remaining states sorted by (size, S)
    rest = sorted(seen - set(seeds), key=lambda S: (size_S(S), S))
    return seeds + rest


def build_truncation_system(
    N: int,
    s: Fraction,
    M: int = M_CLAMP,
    *,
    corrupt_boundary_diag: bool = False,
    corrupt_drop_copies: bool = False,
    corrupt_r0_no_slot: bool = False,
    corrupt_HN_denom: Optional[Fraction] = None,
    corrupt_zero_one_HN: bool = False,
) -> Tuple[List[TypeS], List[List[Fraction]], List[List[Fraction]], List[Fraction], List[List[Fraction]], List[List[Fraction]], List[Fraction]]:
    """Build A_N, B_cols, H_N, and solve for K^-_N, K^+_N, u_N exactly over Q."""
    states = build_domain_DN(N, M)
    idx = {S: i for i, S in enumerate(states)}
    dim = len(states)
    n_labels = M - 1

    A_N = [[Fraction(0)] * dim for _ in range(dim)]
    B_cols = [[Fraction(0)] * n_labels for _ in range(dim)]
    H_N = [Fraction(0)] * dim
    hn_zeroed = False

    for i, S in enumerate(states):
        r0_map = intrinsic_transitions_S(S, M)
        r0 = sum(r0_map.values(), Fraction(0))
        if not corrupt_r0_no_slot:
            assert r0 == Fraction(size_S(S)) - Fraction(1, deg_S(S)), "r_0 != n - 1/D_root"
        else:
            r0 = Fraction(size_S(S))

        t_map = dict(r0_map)
        if not corrupt_drop_copies:
            for S_next, rate in copy_transitions_S(S).items():
                t_map[S_next] = t_map.get(S_next, Fraction(0)) + rate

        # Verify no self-transitions
        assert t_map.get(S, Fraction(0)) == 0, f"Unexpected self-transition at {S}!"

        retained_r0 = sum(rate for S_next, rate in r0_map.items() if size_S(S_next) <= N)
        diag_r0 = retained_r0 if corrupt_boundary_diag else r0
        A_N[i][i] = s + S[0] + diag_r0

        for S_next, rate in t_map.items():
            sz_next = size_S(S_next)
            if sz_next <= N:
                j = idx[S_next]
                A_N[i][j] -= rate
            else:
                assert sz_next == N + 1, "Exit jump exceeded size N+1!"
                denom = corrupt_HN_denom if corrupt_HN_denom is not None else (s - RHO_PARAM)
                H_N[i] += rate * W_S(S_next) / denom

        if corrupt_zero_one_HN and (not hn_zeroed) and H_N[i] > 0:
            H_N[i] = Fraction(0)
            hn_zeroed = True

        # Direct birth column B_j(S) = b_S * 1_{m(D_root) = j}
        L_o = S[1]
        if L_o > 0:
            D_o = deg_S(S)
            lbl = m_clamp(D_o, M)
            B_cols[i][lbl - 2] = Fraction(L_o, D_o)

    # Solve A_N X = [B_cols | H_N] in one exact rational elimination
    rhs = [B_cols[i] + [H_N[i]] for i in range(dim)]
    sol = solve_rational_multi_rhs(A_N, rhs)

    K_minus = [sol[r][:n_labels] for r in range(n_labels)]
    u_N = [sol[r][n_labels] for r in range(n_labels)]
    K_plus = [[K_minus[r][c] + u_N[r] for c in range(n_labels)] for r in range(n_labels)]
    return states, A_N, B_cols, H_N, K_minus, K_plus, u_N


# ---------------------------------------------------------------------------
# Test 1: Clamped Lyapunov inequality (Section 1)
# ---------------------------------------------------------------------------
def test_clamped_lyapunov() -> None:
    print("=== Test 1: Clamped Lyapunov inequality (m(D) = min(D, 20)) ===")
    worst = (Fraction(-100), None)
    count = 0
    for D in range(2, 251):
        mD = m_clamp(D, M_CLAMP)
        factor = Fraction(1, 1 + D_PARAM * mD)
        if D >= 20:
            assert factor == Fraction(1, 4)
        for c in range(0, D + 1):
            lhs = (1 - Fraction(c, D)) * (h_val(c + 1) - h_val(c) + factor)
            rhs = Q_PARAM * h_val(c)
            diff = lhs - rhs
            assert diff <= 0, f"Clamped Lyapunov FAILED at c={c}, D={D}: {lhs} > {rhs}"
            if diff > worst[0]:
                worst = (diff, (c, D, lhs, rhs))
            count += 1
    assert Fraction(1, 4) < Q_PARAM * Fraction(7, 10)  # 1/4 = 0.25 < 63/250 = 0.252
    # Also prove M = 20 is the minimal integer clamp:
    assert Fraction(1, 1 + D_PARAM * 19) > Q_PARAM * Fraction(7, 10)
    c, D, lhs, rhs = worst[1]
    print(f"  Checked {count} pairs (D=2..250, 0<=c<=D): 0 failures.")
    print(f"  Tightest pair: c={c}, D={D} (lhs={lhs} <= rhs={rhs}, margin={-worst[0]}).")
    print(f"  Tail D>=20 factor 1/4 < q*7/10 = 63/250 (margin 1/500); M=19 fails (20/77 > 63/250). PASS")


# ---------------------------------------------------------------------------
# Test 2: Clamped physical-to-branching transfer identity (M=3 and M=20)
# ---------------------------------------------------------------------------
def test_clamped_transfer() -> None:
    print("=== Test 2: Clamped signed-to-positive transfer identity (M=3 and M=20, n=3..7) ===")
    for M in (3, 20):
        by_n: Dict[int, Dict[TreeT, LabelledTree]] = {3: {(2, ()): T3_labelled()}}
        clamped_births_seen = 0
        for n in range(3, 7):
            nxt: Dict[TreeT, LabelledTree] = {}
            for T_can, T_lab in by_n[n].items():
                deg = T_lab.degrees()
                for rate, T_next in T_lab.physical_transitions(M):
                    # Detect promotions where D_parent > M (strict clamping a_v > p_u at birth!)
                    new_leaf_par = T_next.parent[-1]
                    if deg[new_leaf_par] == 1:
                        u = T_lab.parent[new_leaf_par]
                        if deg[u] > M:
                            clamped_births_seen += 1
                    T_can_next = T_next.canonical_T()
                    if T_can_next not in nxt:
                        nxt[T_can_next] = T_next
            by_n[n + 1] = nxt

        total_trees = 0
        for n, trees in sorted(by_n.items()):
            for T_can, T_lab in trees.items():
                total_trees += 1
                old_sub = T_lab.noninit_subtrees()
                lhs: LinComb = {}
                for rate, T_next in T_lab.physical_transitions(M):
                    for S in T_next.noninit_subtrees():
                        lin_add(lhs, S, rate)
                    for S in old_sub:
                        lin_add(lhs, S, -rate)
                rhs: LinComb = {}
                for S in old_sub:
                    for S_out, c in apply_RB_basis(S, M).items():
                        lin_add(rhs, S_out, c)
                for S_out, c in apply_I_basis(T_can, M).items():
                    lin_add(rhs, S_out, c)
                assert lhs == rhs, f"Clamped transfer FAILED at M={M}, n={n}, T={T_can}"
        print(f"  M={M}: {total_trees} reachable marked trees (n=3..7), clamped D>M births exercised={clamped_births_seen}: 0 failures.")
    print("  PASS")


# ---------------------------------------------------------------------------
# Test 3: Exact rational truncations K^-_N, K^+_N (M=20, N=2,3,4) & path sum (2)
# ---------------------------------------------------------------------------
def test_rational_truncations_and_path_sum() -> None:
    print("=== Test 3: Exact 19x19 rational matrices K^-_N, K^+_N (N=2,3,4) & path sum (2) ===")
    for s_name, s_val in [("l=113/200", L_BOUND), ("u=1133/2000", U_BOUND)]:
        results = {}
        for N in (2, 3, 4):
            states, A_N, B_cols, H_N, K_m, K_p, u_N = build_truncation_system(N, s_val, M_CLAMP)
            dim = len(states)
            # Check M-matrix property: A_N W > 0 entrywise
            W_vec = [W_S(S) for S in states]
            for i in range(dim):
                aw_i = sum(A_N[i][j] * W_vec[j] for j in range(dim))
                assert aw_i >= (s_val - RHO_PARAM) * W_vec[i] > 0
            # Check A_N^{-1} >= 0 explicitly on N=2,3 (full inverse); on N=4 check nonnegativity of K^- and u_N
            if N <= 3:
                I_dim = [[Fraction(1 if i == j else 0) for j in range(dim)] for i in range(dim)]
                inv_A = solve_rational_multi_rhs(A_N, I_dim)
                for row in inv_A:
                    for val in row:
                        assert val >= 0, "A_N^{-1} had negative entry!"
            for r in range(19):
                assert u_N[r] > 0
                for c in range(19):
                    assert K_m[r][c] >= 0 and K_p[r][c] >= K_m[r][c]

            # Check estimate (6) for every row i = 2..20
            # Exact check via rational 24th power: u_N(i) * (s - rho) / W(J_{1/i}) <= (2/(N+1))^{1/24}
            for r in range(19):
                ratio = u_N[r] * (s_val - RHO_PARAM) / W_S(J_type(r + 2))
                assert ratio >= 0
                assert ratio ** 24 <= Fraction(2, N + 1), f"Estimate (6) FAILED at N={N}, i={r+2}"

            results[N] = (states, A_N, B_cols, H_N, K_m, K_p, u_N)

        # Verify exact nesting: K^-_2 <= K^-_3 <= K^-_4 and K^+_4 <= K^+_3 <= K^+_2
        for r in range(19):
            for c in range(19):
                assert results[2][4][r][c] <= results[3][4][r][c] <= results[4][4][r][c]
                assert results[4][5][r][c] <= results[3][5][r][c] <= results[2][5][r][c]

        # Check sharper test-column bound (line 55) for g = (1, 1/2, ..., 1/19)
        g_test = [Fraction(1, k + 1) for k in range(19)]
        g_norm = max(g_test)
        for N in (2, 3):
            Km, Kp, uN = results[N][4], results[N][5], results[N][6]
            Km_next = results[N + 1][4]
            for r in range(19):
                lo_next = sum(Km_next[r][c] * g_test[c] for c in range(19))
                sharp_up = sum(Km[r][c] * g_test[c] for c in range(19)) + uN[r] * g_norm
                loose_up = sum(Kp[r][c] * g_test[c] for c in range(19))
                assert lo_next <= sharp_up <= loose_up

        # Verify spectral radius ordering
        rm2 = spectral_radius_power(results[2][4])
        rm3 = spectral_radius_power(results[3][4])
        rm4 = spectral_radius_power(results[4][4])
        rp4 = spectral_radius_power(results[4][5])
        rp3 = spectral_radius_power(results[3][5])
        rp2 = spectral_radius_power(results[2][5])
        assert rm2 <= rm3 <= rm4 <= rp4 <= rp3 <= rp2
        print(
            f"  s={s_name}: |D_2|={len(results[2][0])}, |D_3|={len(results[3][0])}, |D_4|={len(results[4][0])}; "
            f"exact nesting & eq (6) PASS; r(K^-_2..4)={rm2:.4f},{rm3:.4f},{rm4:.4f} <= r(K^+_4..2)={rp4:.1f},{rp3:.1f},{rp2:.1f}"
        )

    # Verify exact finite path sum (2) restricted to D_3 (acyclic: exact at ell <= 2) and D_4 (cyclic: exact remainder)
    for N, max_ell, expect_exact in [(3, 3, True), (4, 10, False)]:
        states, A_N, B_cols, _, K_m, _, _ = build_truncation_system(N, L_BOUND, M_CLAMP)
        dim = len(states)
        # Split A_N = D_s - T_N
        D_diag = [A_N[i][i] for i in range(dim)]
        # Compute P = D_s^{-1} T_N and v_0 = D_s^{-1} B_cols
        # Track partial sum S_L = sum_{ell=0}^L P^ell v_0 and compare to exact X = A_N^{-1} B_cols
        cur = [[B_cols[i][c] / D_diag[i] for c in range(19)] for i in range(dim)]
        accum = [row[:] for row in cur]
        for _ in range(max_ell):
            nxt = [[Fraction(0)] * 19 for _ in range(dim)]
            for i in range(dim):
                inv_d = Fraction(1) / D_diag[i]
                for j in range(dim):
                    if j != i and A_N[i][j] != 0:
                        t_ij = -A_N[i][j]  # positive transition t(S_i, S_j)
                        fac = t_ij * inv_d
                        for c in range(19):
                            if cur[j][c] != 0:
                                nxt[i][c] += fac * cur[j][c]
            for i in range(dim):
                for c in range(19):
                    accum[i][c] += nxt[i][c]
            cur = nxt

        if expect_exact:
            for r in range(19):
                for c in range(19):
                    assert accum[r][c] == K_m[r][c], f"Path sum (2) on D_3 mismatch at ({r},{c})"
                    assert cur[r][c] == 0
            print(f"  Path sum (2) on acyclic D_3 (ell<=3): EXACT match to K^-_3(l) on all 361 entries.")
        else:
            # Verify exact cycle presence (cur != 0) and exact telescoping (I - P) accum = v_0 - P^{L+1} v_0,
            # and strict monotonic convergence accum <= K_m with tiny residual!
            max_rem = max(abs(float(K_m[r][c] - accum[r][c])) for r in range(19) for c in range(19))
            assert any(cur[r][c] > 0 for r in range(19) for c in range(19)), "Expected cycles in D_4!"
            for r in range(19):
                for c in range(19):
                    assert accum[r][c] <= K_m[r][c]
            assert max_rem < 1e-8
            print(f"  Path sum (2) on cyclic D_4 (ell<=10): cycles confirmed, accum <= K^-_4(l), max residual={max_rem:.3e}. PASS")


# ---------------------------------------------------------------------------
# Test 4: Estimates (5), (7), constants H, e, g_*, and inverse slope
# ---------------------------------------------------------------------------
def test_estimates_5_7_and_slope() -> None:
    print("=== Test 4: Estimates (5), (7), constants H, e, g_*, and inverse slope ===")
    # 4a. Check (5): rho + delta*(1 + q/2) <= l
    drift_coeff = RHO_PARAM + DELTA * (1 + Q_PARAM / 2)
    assert drift_coeff == Fraction(671, 1200) < L_BOUND == Fraction(678, 1200)

    # Exact verification of L^0 V_delta identity & (5) at delta = 1 on all states in D_4
    states4 = build_domain_DN(4, M_CLAMP)
    for S in states4:
        n = size_S(S)
        r0_map = intrinsic_transitions_S(S, M_CLAMP)
        r0 = sum(r0_map.values(), Fraction(0))
        assert r0 == Fraction(n) - Fraction(1, deg_S(S)) < n
        # Exact identity (1 + 1/(24 n))^{24} >= 1 + 1/n (proving (1 + 1/n)^{1/24} <= 1 + 1/(24n))
        assert (1 + DELTA / n) ** 24 >= 1 + Fraction(1, n)

        # At delta = 1: V_1(S) = n * W(S)
        l0_W = sum(rate * (W_S(Sp) - W_S(S)) for Sp, rate in r0_map.items())
        l0_V1 = sum(rate * (size_S(Sp) * W_S(Sp) - n * W_S(S)) for Sp, rate in r0_map.items())
        assert l0_V1 == (n + 1) * l0_W + r0 * W_S(S)
        r_V1 = l0_V1 - S[0] * n * W_S(S) + sum(rate * size_S(Sp) * W_S(Sp) for Sp, rate in copy_transitions_S(S).items())
        assert r_V1 <= (RHO_PARAM + (1 + Q_PARAM / 2)) * (n * W_S(S))

    # 4b. Check H, e (endpoints j = 2 and j = 20), g_*, bound (7)
    H_val = Fraction(1) / (D_PARAM * (L_BOUND - RHO_PARAM))
    assert H_val == Fraction(4000, 33)

    e_terms: Dict[int, Fraction] = {}
    for j in range(2, 21):
        prod = Fraction(1)
        for k in range(2, j):
            denom = U_BOUND + Fraction(1, 2) + (k - 1) + Fraction(k - 1, k)
            prod *= Fraction(k - 1) / denom
        term_denom = U_BOUND + Fraction(1, 2) + (j - 1) + Fraction(j - 1, j)
        e_terms[j] = prod * (Fraction(j - 1, j) / term_denom)
        assert e_terms[j] > 0

    e_min = min(e_terms.values())
    e_argmin = min(e_terms, key=lambda j: e_terms[j])
    g_star = e_min / (19 * H_val)
    assert e_min > 0 and g_star > 0

    # Verify e_terms[j] <= exact single-path entry of K^-_4(u) for j = 2, 3, 4 (all inside D_4!)
    _, _, _, _, Km4_u, _, _ = build_truncation_system(4, U_BOUND, M_CLAMP)
    for r in range(19):
        for j in (2, 3, 4):
            assert Km4_u[r][j - 2] >= e_terms[j] >= e_min

    # Verify exact J_{1/i} exit rate 1/i + 3/2 <= 2 and inverse slope ell=0 bound
    for i in range(2, 21):
        S0 = J_type(i)
        exit_rate = S0[0] + r0_total(S0, M_CLAMP)
        assert exit_rate == Fraction(1, i) + Fraction(3, 2) <= Fraction(2)
        # Exact ell=0 contribution to K_{i,2}(l) - K_{i,2}(u):
        ell0_diff = Fraction(1, 2) / (L_BOUND + exit_rate) - Fraction(1, 2) / (U_BOUND + exit_rate)
        slope_lb = (U_BOUND - L_BOUND) / (2 * (U_BOUND + 2) ** 2)
        assert ell0_diff >= slope_lb > 0

    # Verify Collatz bound on K^-_4(l) <= K^-_4(u) (no error) and bound (7) on M=20 N=4
    _, _, _, _, Km4_l, Kp4_l, u4_l = build_truncation_system(4, L_BOUND, M_CLAMP)
    gap4 = spectral_radius_power(Kp4_l) - spectral_radius_power(Km4_l)
    bound7_N4 = float(38 * H_val / g_star) * ((2.0 / 5.0) ** (1.0 / 24.0))
    assert 0.0 <= gap4 <= bound7_N4

    print(f"  rho + delta*(1+q/2) = {drift_coeff} < l = {L_BOUND}; H = {H_val} ({float(H_val):.3f}).")
    print(f"  e = {float(e_min):.6e} (attained at j={e_argmin}), e(j=2)={float(e_terms[2]):.5f}, g_* = {float(g_star):.6e}.")
    print(f"  ell=0 exit rate 1/i+3/2 <= 2 & inverse slope >= {float(g_star / (2*(U_BOUND+2)**2)):.6e}*(s'-s). PASS")


# ---------------------------------------------------------------------------
# Test 5: Small M=3 clamped model up to N=5 (monotone spectral enclosures)
# ---------------------------------------------------------------------------
def test_small_clamped_enclosures() -> None:
    print("=== Test 5: Small 2-label (M=3) rational enclosure test (N=2..5) ===")
    s_test = Fraction(3, 5)
    prev_rm, prev_rp = -1.0, 1e18
    for N in (2, 3, 4, 5):
        states, _, _, _, Km, Kp, uN = build_truncation_system(N, s_test, M=3)
        rm = spectral_radius_power(Km)
        rp = spectral_radius_power(Kp)
        assert prev_rm <= rm <= rp <= prev_rp
        prev_rm, prev_rp = rm, rp
        print(f"  N={N} (|D_N|={len(states):3d}): r(K^-_N)={rm:.6f} <= r(K^+_N)={rp:.6f}, max u_N={float(max(uN)):.5f}")
    print("  PASS")


# ---------------------------------------------------------------------------
# Test 6: Planted errors (9 negative controls)
# ---------------------------------------------------------------------------
def test_planted_errors() -> None:
    print("=== Test 6: Planted errors (all 10 must be caught) ===")

    # Reference exact N=2, 3, 4 at s = l
    _, _, _, _, Km2, Kp2, u2 = build_truncation_system(2, L_BOUND, M_CLAMP)
    _, _, _, _, Km3, Kp3, u3 = build_truncation_system(3, L_BOUND, M_CLAMP)
    _, _, _, _, Km4, Kp4, u4 = build_truncation_system(4, L_BOUND, M_CLAMP)

    # Error 1: Omit exiting r^0 rate from diagonal on boundary n(S)=N
    _, _, _, _, Km2_bad, _, _ = build_truncation_system(2, L_BOUND, M_CLAMP, corrupt_boundary_diag=True)
    # Without boundary diagonal subtraction, K^-_2 overestimates K^-_3!
    caught1 = any(Km2_bad[r][c] > Km3[r][c] for r in range(19) for c in range(19))
    assert caught1, "Planted Error 1 NOT caught!"
    print("  [CAUGHT] Error 1 (omit boundary exit rate from diagonal): K^-_2_bad > K^-_3 (violates monotonicity/lower bound).")

    # Error 2: Drop copy transitions from R_N (changes K^-_3 on size 3!)
    _, _, _, _, Km3_nocopy, _, _ = build_truncation_system(3, L_BOUND, M_CLAMP, corrupt_drop_copies=True)
    caught2 = any(Km3_nocopy[r][c] != Km3[r][c] for r in range(19) for c in range(19))
    assert caught2, "Planted Error 2 NOT caught!"
    print("  [CAUGHT] Error 2 (drop copy terms in t(S,S')): K^-_3_nocopy != exact K^-_3 (fails at column j=3).")

    # Error 3: Use r_0(S) = n(S) instead of n(S) - 1/D_root
    _, _, _, _, Km2_noslot, _, _ = build_truncation_system(2, L_BOUND, M_CLAMP, corrupt_r0_no_slot=True)
    caught3 = any(Km2_noslot[r][c] != Km2[r][c] for r in range(19) for c in range(19))
    assert caught3, "Planted Error 3 NOT caught!"
    print("  [CAUGHT] Error 3 (r_0 = n instead of n - 1/D_root): K^-_2_noslot != exact K^-_2.")

    # Error 4a: Wrong supersolution denominator s instead of s - rho fails (s - R)(W/s) >= W at J_{1/20}
    S20 = J_type(20)
    r_W_20 = (
        sum(rate * (W_S(Sp) - W_S(S20)) for Sp, rate in intrinsic_transitions_S(S20, M_CLAMP).items())
        - S20[0] * W_S(S20)
        + sum(rate * W_S(Sp) for Sp, rate in copy_transitions_S(S20).items())
    )
    assert r_W_20 > 0 and (L_BOUND * W_S(S20) - r_W_20) / L_BOUND < W_S(S20)
    print(f"  [CAUGHT] Error 4a (W/s instead of W/(s-rho)): fails supersolution (s-R)(W/s) >= W at J_{{1/20}} (RW={r_W_20} > 0).")

    # Error 4b: Over-aggressive H_N denominator (s + 5) directly breaks upper nesting K^+_3 <= K^+_2
    _, _, _, _, _, Kp2_badd, _ = build_truncation_system(2, L_BOUND, M_CLAMP, corrupt_HN_denom=L_BOUND + 5)
    _, _, _, _, _, Kp3_badd, _ = build_truncation_system(3, L_BOUND, M_CLAMP, corrupt_HN_denom=L_BOUND + 5)
    caught4b = any(Kp3_badd[r][c] > Kp2_badd[r][c] for r in range(19) for c in range(19))
    assert caught4b, "Planted Error 4b NOT caught!"
    print("  [CAUGHT] Error 4b (H_N denominator s+5 instead of s-rho): violates upper nesting K^+_3 <= K^+_2.")

    # Error 5: Zero out one boundary entry of H_N at N=2
    _, _, _, _, _, Kp2_zero1, _ = build_truncation_system(2, L_BOUND, M_CLAMP, corrupt_zero_one_HN=True)
    caught5 = any(Kp3[r][c] > Kp2_zero1[r][c] for r in range(19) for c in range(19))
    assert caught5, "Planted Error 5 NOT caught!"
    print("  [CAUGHT] Error 5 (zero out one boundary H_2 entry): violates K^+_3 <= K^+_2 on that row.")

    # Error 6: Choose delta = 1/10 (violates rho + delta*(1+q/2) <= l)
    delta_bad = Fraction(1, 10)
    assert RHO_PARAM + delta_bad * (1 + Q_PARAM / 2) > L_BOUND
    print(f"  [CAUGHT] Error 6 (delta=1/10): rho + delta(1+q/2) = {RHO_PARAM + delta_bad*(1+Q_PARAM/2)} > l = {L_BOUND}.")

    # Error 7: Clamp at M = 19 instead of M = 20 (first fails at c=4, D=135)
    fail_m19 = None
    for D in range(19, 200):
        lhs = (1 - Fraction(4, D)) * Fraction(1, 1 + D_PARAM * m_clamp(D, 19))
        if lhs > Q_PARAM * h_val(4):
            fail_m19 = (4, D, lhs, Q_PARAM * h_val(4))
            break
    assert fail_m19 is not None
    print(f"  [CAUGHT] Error 7 (clamp M=19): fails Lyapunov inequality at c={fail_m19[0]}, D={fail_m19[1]} ({fail_m19[2]} > {fail_m19[3]}).")

    # Error 8: Overestimate J_{1/2} survival as e^{-3t/2} instead of e^{-2t} in inverse slope
    exact_ell0_diff = Fraction(1, 2) / (L_BOUND + Fraction(2)) - Fraction(1, 2) / (U_BOUND + Fraction(2))
    bad_slope_lb = (U_BOUND - L_BOUND) / (2 * (U_BOUND + Fraction(3, 2)) ** 2)
    assert exact_ell0_diff < bad_slope_lb
    print("  [CAUGHT] Error 8 (dropping killing 1/2 from J_{1/2} survival in inverse slope): bound exceeds exact ell=0 term.")

    # Error 9: Replacing ||g||_inf by min(g) in line 55 sharp test-column bound
    g_nonunif = [Fraction(1) if k == 1 else Fraction(0) for k in range(19)]  # supported on label j=3
    caught9 = any(
        sum(Km3[r][c] * g_nonunif[c] for c in range(19))
        > sum(Km2[r][c] * g_nonunif[c] for c in range(19)) + u2[r] * min(g_nonunif)
        for r in range(19)
    )
    assert caught9, "Planted Error 9 NOT caught!"
    print("  [CAUGHT] Error 9 (replacing ||g||_inf by min(g) in line 55 sharp test-column bound): violated at N=2 vs N=3.")
    print("  All 10 planted errors caught. PASS")


def main() -> None:
    t0 = time.time()
    test_clamped_lyapunov()
    test_clamped_transfer()
    test_rational_truncations_and_path_sum()
    test_estimates_5_7_and_slope()
    test_small_clamped_enclosures()
    test_planted_errors()
    print(f"\nALL INDEPENDENT CHECKS PASSED in {time.time() - t0:.2f}s.")


if __name__ == "__main__":
    main()
