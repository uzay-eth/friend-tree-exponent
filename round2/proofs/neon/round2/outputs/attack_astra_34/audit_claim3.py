#!/usr/bin/env python3
"""
Independent adversarial verification & planted-error test suite for Astra Claim 3:
  - round2/fringe/asymptotics.md
  - round2/degree-law.md
  - dependencies in round2/fringe/renewal.md

Written completely from scratch (no imports from ft2-astra/).
"""

from __future__ import annotations

import cmath
from collections import defaultdict
from fractions import Fraction
import math
from typing import Callable, Dict, List, Sequence, Tuple


# ============================================================================
# 1. Exact small-n enumeration of the discrete friend-tree chain (n = 2..10)
# ============================================================================

def tree_canonical_key(n: int, edges: Tuple[Tuple[int, int], ...]) -> Tuple:
    """AHU canonical isomorphism signature of an unlabelled free tree."""
    if n == 1:
        return ()
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    deg = [len(adj[u]) for u in range(n)]
    rem = n
    leaves = [u for u in range(n) if deg[u] <= 1]
    while rem > 2:
        nxt = []
        for u in leaves:
            rem -= 1
            for v in adj[u]:
                deg[v] -= 1
                if deg[v] == 1:
                    nxt.append(v)
        leaves = nxt

    def enc(u: int, p: int) -> Tuple:
        ch = [enc(v, u) for v in adj[u] if v != p]
        ch.sort()
        return tuple(ch)

    return min(enc(c, -1) for c in leaves)


def exact_friend_tree_chain(n_max: int = 10):
    """
    Enumerate exact rational distributions of T_n for n = 2..n_max started from T_2.
    Returns:
      dists[n]: dict[edges -> Fraction probability]
      E_ge[(n, k)]: exact E[N_{>=k}(n)]
      E_eq[(n, k)]: exact E[N_{=k}(n)]
    """
    dists: Dict[int, Dict[Tuple[Tuple[int, int], ...], Fraction]] = {}
    cur: Dict[Tuple, Tuple[Tuple[Tuple[int, int], ...], Fraction]] = {
        tree_canonical_key(2, ((0, 1),)): (((0, 1),), Fraction(1, 1))
    }

    E_ge: Dict[Tuple[int, int], Fraction] = {}
    E_eq: Dict[Tuple[int, int], Fraction] = {}

    for n in range(2, n_max + 1):
        dists[n] = {edges: prob for edges, prob in cur.values()}
        for k in range(1, n + 1):
            s_ge = Fraction(0, 1)
            s_eq = Fraction(0, 1)
            for edges, prob in dists[n].items():
                deg = [0] * n
                for u, v in edges:
                    deg[u] += 1
                    deg[v] += 1
                s_ge += prob * sum(1 for d in deg if d >= k)
                s_eq += prob * sum(1 for d in deg if d == k)
            E_ge[(n, k)] = s_ge
            E_eq[(n, k)] = s_eq

        if n < n_max:
            nxt_prob: Dict[Tuple, Fraction] = defaultdict(Fraction)
            nxt_rep: Dict[Tuple, Tuple[Tuple[int, int], ...]] = {}
            inv_n = Fraction(1, n)
            for edges, prob in dists[n].items():
                adj: List[List[int]] = [[] for _ in range(n)]
                for u, v in edges:
                    adj[u].append(v)
                    adj[v].append(u)
                # Attachment probability at target w: (1/n) * sum_{u ~ w} (1 / D_u)
                w_prob = [Fraction(0, 1)] * n
                for u in range(n):
                    share = inv_n * Fraction(1, len(adj[u]))
                    for w in adj[u]:
                        w_prob[w] += share
                assert sum(w_prob) == Fraction(1, 1)
                for w in range(n):
                    new_edges = edges + ((w, n),)
                    key = tree_canonical_key(n + 1, new_edges)
                    if key not in nxt_rep:
                        nxt_rep[key] = new_edges
                    nxt_prob[key] += prob * w_prob[w]
            cur = {key: (nxt_rep[key], p) for key, p in nxt_prob.items()}

    return dists, E_ge, E_eq


# ============================================================================
# 2. Exact Yule Laplace identity & T_2 <-> T_3 Exp(2) convolution check
# ============================================================================

def verify_yule_and_t2_t3_convolution(
    dists, E_ge, n_max: int = 10,
    plant_wrong_yule_weight: bool = False,
    plant_wrong_t2_t3_factor: bool = False,
) -> Tuple[bool, str]:
    """
    Check:
      (a) Exact CTMC evolution of T_2 vs Yule mixture q^2 sum_{n>=2} (n-1)(1-q)^{n-2} a_{n,k}:
          In continuous time from T_2, holding time at size n is Exp(n), independent of jump chain.
          Solve exact Kolmogorov forward ODEs for P(N(t)=n) for n=2..n_max via Taylor series
          d^m/dt^m|_{t=0} P(N(t)=n) and verify exact rational identity with
          q^2 (n-1) (1-q)^{n-2} expanded in t (where q = e^{-t})!
      (b) Exact Laplace transform ratio between T_2 and T_3:
          L{P_{T_2}(N(t)=n); z} / L{P_{T_3}(N(t)=n); z} == 2 / (2 + z) for every n >= 3.
    """
    # (a) Exact derivatives at t=0 of P_{T_2}(N(t)=n) from birth process generator:
    # p'_2 = -2 p_2; p'_n = (n-1) p_{n-1} - n p_n (n >= 3).
    max_ord = 8
    # dp[m][n] = d^m/dt^m|_{t=0} P(N(t)=n)
    dp = [[Fraction(0, 1)] * (n_max + 1) for _ in range(max_ord + 1)]
    dp[0][2] = Fraction(1, 1)
    for m in range(1, max_ord + 1):
        dp[m][2] = -2 * dp[m - 1][2]
        for n in range(3, n_max + 1):
            dp[m][n] = (n - 1) * dp[m - 1][n - 1] - n * dp[m - 1][n]

    # Compare with d^m/dt^m|_{t=0} [ (n-1) e^{-2t} (1 - e^{-t})^{n-2} ]
    # Expand (1 - e^{-t})^{n-2} = sum_{j=0}^{n-2} (-1)^j C(n-2, j) e^{-jt}
    # => e^{-2t}(1-e^{-t})^{n-2} = sum_{j=0}^{n-2} (-1)^j C(n-2, j) e^{-(j+2)t}
    # m-th derivative at t=0 is sum_{j=0}^{n-2} (-1)^j C(n-2, j) (-(j+2))^m.
    for n in range(2, n_max + 1):
        coeff = Fraction(n if plant_wrong_yule_weight else (n - 1), 1)
        for m in range(max_ord + 1):
            formula_deriv = Fraction(0, 1)
            for j in range(n - 1):
                c_nj = math.comb(n - 2, j)
                sign = -1 if (j % 2 == 1) else 1
                formula_deriv += coeff * sign * c_nj * ((-(j + 2)) ** m)
            if formula_deriv != dp[m][n]:
                return False, f"Yule Taylor mismatch at n={n}, order m={m}: {formula_deriv} != {dp[m][n]}"

    # Also verify E_{T_2}[N_{>=k}(t)] Taylor derivatives up to order max_ord match
    # direct sum_n dp[m][n] * a_{n,k} vs formula derivative for all k=2..n_max:
    for k in range(2, n_max + 1):
        for m in range(max_ord + 1):
            lhs = sum(dp[m][n] * E_ge.get((n, k), Fraction(0, 1)) for n in range(2, n_max + 1))
            rhs = Fraction(0, 1)
            for n in range(2, n_max + 1):
                w = Fraction(n if plant_wrong_yule_weight else (n - 1), 1)
                d_yule = sum(
                    (-1 if j % 2 else 1) * math.comb(n - 2, j) * ((-(j + 2)) ** m)
                    for j in range(n - 1)
                )
                rhs += w * d_yule * E_ge.get((n, k), Fraction(0, 1))
            if lhs != rhs:
                return False, f"E_T2 N_{{>={k}}} derivative mismatch at m={m}"

    # (b) Exact Laplace transform of P_{T_2}(N(t)=n) vs P_{T_3}(N(t)=n):
    # Holding times from T_2: Exp(2), Exp(3), ..., Exp(n) ->
    # L{P_{T_2}(N(t)=n); z} = (1/(z+n)) * prod_{j=2}^{n-1} (j / (z + j))
    # From T_3: Exp(3), ..., Exp(n) ->
    # L{P_{T_3}(N(t)=n); z} = (1/(z+n)) * prod_{j=3}^{n-1} (j / (z + j))
    # Also verify via Beta integral int_0^1 q^{z-1} (n-1) q^2 (1-q)^{n-2} dq!
    for z in (Fraction(5657, 10000), Fraction(113, 200), Fraction(1, 2), Fraction(1, 1)):
        factor = Fraction(1, 1 + z) if plant_wrong_t2_t3_factor else Fraction(2, 2 + z)
        for n in range(3, n_max + 1):
            # Beta integral for T_2: (n-1) * (n-2)! / prod_{j=2}^n (z + j)
            denom2 = Fraction(1, 1)
            for j in range(2, n + 1):
                denom2 *= (z + j)
            lap2 = Fraction(math.factorial(n - 1), 1) / denom2

            # Beta integral for T_3: ((n-1)(n-2)/2) * (n-3)! / prod_{j=3}^n (z + j)
            denom3 = Fraction(1, 1)
            for j in range(3, n + 1):
                denom3 *= (z + j)
            lap3 = Fraction(math.factorial(n - 1), 2) / denom3

            if lap2 != factor * lap3:
                return False, f"T2/T3 Laplace ratio mismatch at z={z}, n={n}: {lap2} != {factor}*{lap3}"

    return True, "Exact Yule Taylor identity (to order 8) and T_2/T_3 factor 2/(2+z) verified."


# ============================================================================
# 3. Degree-crossing identity (degree-law.md eq. 35) exact small-n check
# ============================================================================

def verify_crossing_identity(
    dists, E_ge, n_max: int = 10,
    drop_rho: bool = False,
    wrong_k_minus_1: bool = False,
) -> Tuple[bool, str]:
    """
    Check the exact continuous-time derivative identity behind degree-law.md eq. (35):
      d/dt E[N_{>=k+1}(t)] = E[ sum_{u: D_u = k} (k - c_u + rho_u) ]
    On a fixed size-n tree T, the instantaneous rate at which N_{>=k+1} increases by 1
    is sum_{u: D_u = k} (L_u + rho_u) = sum_{u: D_u = k} (k - c_u + rho_u).
    In discrete time, n * (E[N_{>=k+1}(n+1)] - E[N_{>=k+1}(n)]) MUST equal
    E_{T_n}[ sum_{u: D_u = k} (k - c_u + rho_u) ] for EVERY n >= 3 and EVERY k >= 2!
    (Because holding time at size n is Exp(n), so d/dt|_{size n} = n * E[Delta | T_n]!)
    """
    for n in range(3, n_max):
        for k in range(2, n + 1):
            lhs = Fraction(n, 1) * (E_ge.get((n + 1, k + 1), Fraction(0, 1)) - E_ge.get((n, k + 1), Fraction(0, 1)))
            rhs = Fraction(0, 1)
            for edges, prob in dists[n].items():
                adj: List[List[int]] = [[] for _ in range(n)]
                for u, v in edges:
                    adj[u].append(v)
                    adj[v].append(u)
                deg = [len(adj[u]) for u in range(n)]
                tree_sum = Fraction(0, 1)
                for u in range(n):
                    if deg[u] == k:
                        core_nbrs = [w for w in adj[u] if deg[w] >= 2]
                        c_u = len(core_nbrs)
                        rho_u = Fraction(0, 1) if drop_rho else sum(
                            Fraction(1, deg[w]) for w in core_nbrs
                        )
                        base = Fraction(k - 1, 1) if wrong_k_minus_1 else Fraction(k - c_u, 1)
                        tree_sum += base + rho_u
                rhs += prob * tree_sum
            if lhs != rhs:
                return False, f"Crossing identity failed at n={n}, k={k}: lhs={lhs} != rhs={rhs}"

    # Also check continuous-time derivative d/dt E_{T_2}[N_{>=k+1}(t)] at t = 0.35
    # truncated at n_max=10 (exact up to n_max-1 via holding times:
    # d/dt sum_{n=2}^N P(N(t)=n) a_{n,k+1} = sum_{n=2}^{N-1} n P(N(t)=n) (a_{n+1,k+1} - a_{n,k+1}) - N P(N(t)=N) a_{N,k+1})
    return True, f"Exact crossing identity n*(E[N_{{>=k+1}}(n+1)] - E[N_{{>=k+1}}(n)]) == E[sum_{{D_u=k}} (k - c_u + rho_u)] verified for all 3<=n<={n_max-1}, 2<=k<=n."


# ============================================================================
# 4. Resolvent / Laplace transform formula (3) of asymptotics.md at small t
# ============================================================================

class MarkedTree:
    """Actual tree from T_3 rooted at initial centre 0 with permanent marks a_v."""
    def __init__(self, n: int, edges: Tuple[Tuple[int, int], ...], marks: Dict[int, Fraction]):
        self.n = n
        self.edges = edges
        self.marks = marks
        self.adj: List[List[int]] = [[] for _ in range(n)]
        for u, v in edges:
            self.adj[u].append(v)
            self.adj[v].append(u)
        self.parent = [-1] * n
        self.children: List[List[int]] = [[] for _ in range(n)]
        stack = [0]
        seen = {0}
        while stack:
            u = stack.pop()
            for v in self.adj[u]:
                if v not in seen:
                    seen.add(v)
                    self.parent[v] = u
                    self.children[u].append(v)
                    stack.append(v)

    def deg(self, u: int) -> int:
        return len(self.adj[u])

    def is_core(self, u: int) -> bool:
        return len(self.adj[u]) >= 2

    def canonical_planted_subtree(self, v: int) -> Tuple:
        """Canonical representation of planted type S_v for v != 0 (core)."""
        assert v != 0 and self.is_core(v)
        def enc(u: int) -> Tuple:
            if not self.is_core(u):
                return ("L",)
            ch = sorted(enc(w) for w in self.children[u])
            mk = self.marks[u]
            return ("C", mk.numerator, mk.denominator, tuple(ch))
        return enc(v)


def planted_transitions_and_rates(S: Tuple) -> Tuple[
    int, Fraction, Fraction, List[Tuple[Tuple, Fraction]], List[Tuple[Tuple, Fraction]], Tuple
]:
    """
    Given a canonical planted type S = ("C", num, den, children_tuple),
    compute:
      - D_root (including 1 parent slot!)
      - a_S (root mark)
      - b_S = L_root / D_root
      - intrinsic L^0 transitions: list of (S_next, rate)
      - residual copy births: list of (Y_or_Yplus, rate)
      - direct B birth type J_{p_S}
    """
    assert S[0] == "C"
    a_S = Fraction(S[1], S[2])
    ch = list(S[3])
    n_leaves = sum(1 for c in ch if c == ("L",))
    core_children = [c for c in ch if c != ("L",)]
    D_root = len(ch) + 1  # +1 for exterior parent slot!
    p_S = Fraction(1, D_root)
    b_S = Fraction(n_leaves, D_root)

    def add_leaf_at_root(node: Tuple) -> Tuple:
        new_ch = sorted(list(node[3]) + [("L",)])
        return ("C", node[1], node[2], tuple(new_ch))

    def intrinsic_recursive(node: Tuple, is_root: bool) -> List[Tuple[Tuple, Fraction]]:
        """
        Return list of (new_node, rate) for all intrinsic L^0 transitions strictly inside `node`,
        with exterior parent slot inert (every node has 1 parent slot!).
        """
        res: List[Tuple[Tuple, Fraction]] = []
        nch = list(node[3])
        d_node = len(nch) + 1  # always has a parent slot (core parent or exterior slot)
        inv_d = Fraction(1, d_node)
        leaf_indices = [i for i, c in enumerate(nch) if c == ("L",)]
        core_indices = [i for i, c in enumerate(nch) if c != ("L",)]

        # 1. Leaves of `node` ring (each at rate 1) -> add leaf to `node`
        if leaf_indices:
            res.append((add_leaf_at_root(node), Fraction(len(leaf_indices), 1)))

        # 2. `node` rings (rate 1):
        #    - picks a leaf (rate len(leaf_indices)/d_node): promotes leaf to J_{1/d_node}
        if leaf_indices:
            idx0 = leaf_indices[0]
            promoted = ("C", 1, d_node, (("L",),))
            new_ch = nch[:idx0] + [promoted] + nch[idx0 + 1:]
            new_ch.sort()
            res.append((("C", node[1], node[2], tuple(new_ch)), Fraction(len(leaf_indices), d_node)))

        #    - picks a core child i (rate 1/d_node): adds leaf to child i
        for i in core_indices:
            new_ch = nch[:i] + [add_leaf_at_root(nch[i])] + nch[i + 1:]
            new_ch.sort()
            res.append((("C", node[1], node[2], tuple(new_ch)), inv_d))

        # 3. Core child i sends to `node` (rate 1/D_child) -> adds leaf to `node`
        for i in core_indices:
            d_child = len(nch[i][3]) + 1
            res.append((add_leaf_at_root(node), Fraction(1, d_child)))

        # 4. Internal transitions strictly inside core child i:
        for i in core_indices:
            for child_next, r in intrinsic_recursive(nch[i], is_root=False):
                new_ch = nch[:i] + [child_next] + nch[i + 1:]
                new_ch.sort()
                res.append((("C", node[1], node[2], tuple(new_ch)), r))

        return res

    l0_trans = intrinsic_recursive(S, is_root=True)
    # Aggregate identical successor states
    l0_agg: Dict[Tuple, Fraction] = defaultdict(Fraction)
    for nxt, r in l0_trans:
        l0_agg[nxt] += r

    copy_births: Dict[Tuple, Fraction] = defaultdict(Fraction)
    for Y in core_children:
        a_Y = Fraction(Y[1], Y[2])
        assert a_Y >= p_S, f"Mark inequality violated: a_Y={a_Y} < p_S={p_S}"
        if a_Y > p_S:
            copy_births[Y] += (a_Y - p_S)
        copy_births[add_leaf_at_root(Y)] += p_S

    J_pS = ("C", p_S.numerator, p_S.denominator, (("L",),))
    return D_root, a_S, b_S, list(l0_agg.items()), list(copy_births.items()), J_pS


def verify_small_t_resolvent_expansion() -> Tuple[bool, str]:
    """
    Directly verifies the small-t Taylor expansion (orders m = 0, 1, 2, 3, 4) of
      M_t f_k  vs  int_0^t I_u P_{t-u} f_k du
    AND the generation decomposition
      P_t f_k = R_t f_k + (R * B * R)_t f_k + (R * B * R * B * R)_t f_k + ...
    whose Laplace transform IS eq. (3) of asymptotics.md:
      L{M_t f_k}(z) = I_hat_z(G_z f_k + U_z (I - K_z)^{-1} j G_z f_k)!
    We verify exact rational equality of all Taylor coefficients t^m/m! for m = 0..4
    and k = 2, 3, 4, 5!
    """
    max_ord = 4

    # 1. Compute exact d^m/dt^m|_{t=0} M_t f_k directly from the physical T_3 CTMC!
    t3 = MarkedTree(3, ((0, 1), (0, 2)), {})

    def phys_sig(t: MarkedTree) -> Tuple:
        def enc(u: int) -> Tuple:
            if not t.is_core(u):
                return ("L",)
            mk = t.marks.get(u, Fraction(-1, 1))
            ch = sorted(enc(w) for w in t.children[u])
            return ("C", mk.numerator, mk.denominator, tuple(ch))
        return enc(0)

    def phys_step_list(t: MarkedTree) -> List[Tuple[MarkedTree, Fraction]]:
        out = []
        for u in range(t.n):
            du = t.deg(u)
            rate = Fraction(1, du)
            for w in t.adj[u]:
                nm = dict(t.marks)
                if not t.is_core(w):
                    nm[w] = Fraction(1, du)
                out.append((MarkedTree(t.n + 1, t.edges + ((w, t.n),), nm), rate))
        return out

    def phys_eval_M(t: MarkedTree, k: int) -> Fraction:
        return Fraction(sum(1 for v in range(1, t.n) if t.is_core(v) and t.deg(v) >= k), 1)

    memo_phys: Dict[Tuple, Fraction] = {}
    def phys_M_deriv(t: MarkedTree, m: int, k: int) -> Fraction:
        if m == 0:
            return phys_eval_M(t, k)
        key = (phys_sig(t), m, k)
        if key in memo_phys:
            return memo_phys[key]
        acc = -Fraction(t.n, 1) * phys_M_deriv(t, m - 1, k)
        for nxt, r in phys_step_list(t):
            acc += r * phys_M_deriv(nxt, m - 1, k)
        memo_phys[key] = acc
        return acc

    # 2. Compute exact d^m/dt^m|_{t=0} R_t f_k(S) and d^m/dt^m|_{t=0} (R_t B_op g)(S)
    # and generation terms P_t^{(g)} f_k(S) via Leibniz convolution:
    # (A * C)^{(m)}(0) = sum_{j=0}^{m-1} C(m-1, j) A^{(j)}(0) C^{(m-1-j)}(0) for int_0^t A(u) C(t-u) du!
    # Check: d/dt int_0^t A(u) C(t-u) du = A(t)C(0) + int_0^t A(u) C'(t-u) du
    # By induction, (A * C)^{(m)}(0) = sum_{j=0}^{m-1} A^{(j)}(0) C^{(m-1-j)}(0) * 1?
    # WAIT! Check binomial coefficient for Taylor coefficients vs derivatives!!
    # If A(t) = sum a_j t^j/j! and C(t) = sum c_r t^r/r!,
    # int_0^t u^j (t-u)^r du = j! r! / (j+r+1)! * t^{j+r+1}.
    # So the coefficient of t^m / m! in (A * C)(t) is:
    # sum_{j+r = m-1} a_j c_r  (NO binomial coefficient — the factorials cancel!!).
    # Let's verify: (1 * 1)(t) = t -> m=1: a_0 c_0 = 1. (t * 1)(t) = t^2/2 -> m=2: a_1 c_0 = 1. EXACT!

    memo_R: Dict[Tuple, Fraction] = {}
    def R_deriv(S: Tuple, m: int, fn: Callable[[Tuple], Fraction], fn_tag: Tuple) -> Fraction:
        """Exact m-th derivative at t=0 of (R_t fn)(S)."""
        if m == 0:
            return fn(S)
        key = (S, m, fn_tag)
        if key in memo_R:
            return memo_R[key]
        D_root, a_S, b_S, l0_trans, copy_births, _ = planted_transitions_and_rates(S)
        tot_l0 = sum((r for _, r in l0_trans), Fraction(0, 1))
        acc = -(tot_l0 + a_S) * R_deriv(S, m - 1, fn, fn_tag)
        for nxt, r in l0_trans:
            acc += r * R_deriv(nxt, m - 1, fn, fn_tag)
        for child_t, r in copy_births:
            acc += r * R_deriv(child_t, m - 1, fn, fn_tag)
        memo_R[key] = acc
        return acc

    def make_f_k(k: int) -> Callable[[Tuple], Fraction]:
        return lambda S: Fraction(1 if (len(S[3]) + 1) >= k else 0, 1)

    # Generation-g derivatives of P_t^{(g)} f_k(S) at t=0:
    # P^{(0)} = R_t f_k
    # P^{(g)} = int_0^t R_u [b_S P^{(g-1)}_{t-u} f_k(J_{p_S})] du
    # So (P^{(g)} f_k)^{(m)}(S) = sum_{j=0}^{m-1} (R^{(j)} [S' -> b_{S'} * (P^{(g-1)} f_k)^{(m-1-j)}(J_{p_{S'}})])(S)!
    # Notice: for g > m, (P^{(g)})^{(m)}(0) == 0 because it has g convolutions!
    memo_Pg: Dict[Tuple, Fraction] = {}
    def Pg_deriv(S: Tuple, m: int, g: int, k: int) -> Fraction:
        if g == 0:
            return R_deriv(S, m, make_f_k(k), ("fk", k))
        if m < g:
            return Fraction(0, 1)
        key = (S, m, g, k)
        if key in memo_Pg:
            return memo_Pg[key]
        acc = Fraction(0, 1)
        for j in range(g - 1, m):
            # r = m - 1 - j goes from 0 to m - g; actually j is order on R_u, r = m-1-j is order on P^{(g-1)}
            pass
        for j in range(0, m):
            r_ord = m - 1 - j
            if r_ord < g - 1:
                continue
            def B_of_Pg(S_arg: Tuple, rr=r_ord, gg=g - 1, kk=k) -> Fraction:
                D_r, _, b_r, _, _, J_p = planted_transitions_and_rates(S_arg)
                return b_r * Pg_deriv(J_p, rr, gg, kk)
            acc += R_deriv(S, j, B_of_Pg, ("B_Pg", r_ord, g - 1, k))
        memo_Pg[key] = acc
        return acc

    # Now compute d^m/dt^m|_{t=0} [ int_0^t I_u (sum_{g=0}^m P^{(g)}_{t-u} f_k) du ]
    # = sum_{j=0}^{m-1} I^{(j)} ( sum_{g=0}^{m-1-j} (P^{(g)} f_k)^{(m-1-j)} )!
    memo_I: Dict[Tuple, Fraction] = {}
    def I_eval(t: MarkedTree, score: Callable[[Tuple], Fraction]) -> Fraction:
        d0 = t.deg(0)
        p0 = Fraction(1, d0)
        leaves0 = [w for w in t.children[0] if not t.is_core(w)]
        cores0 = [w for w in t.children[0] if t.is_core(w)]
        b0 = Fraction(len(leaves0), d0)
        J_p0 = ("C", p0.numerator, p0.denominator, (("L",),))
        res = b0 * score(J_p0)
        for y in cores0:
            Y = t.canonical_planted_subtree(y)
            a_Y = t.marks[y]
            Y_plus = ("C", Y[1], Y[2], tuple(sorted(list(Y[3]) + [("L",)])))
            res += (a_Y - p0) * score(Y) + p0 * score(Y_plus)
        return res

    def I_deriv(t: MarkedTree, m: int, score: Callable[[Tuple], Fraction], tag: Tuple) -> Fraction:
        if m == 0:
            return I_eval(t, score)
        key = (phys_sig(t), m, tag)
        if key in memo_I:
            return memo_I[key]
        acc = -Fraction(t.n, 1) * I_deriv(t, m - 1, score, tag)
        for nxt, r in phys_step_list(t):
            acc += r * I_deriv(nxt, m - 1, score, tag)
        memo_I[key] = acc
        return acc

    for k in (2, 3, 4, 5):
        for m in range(max_ord + 1):
            lhs = phys_M_deriv(t3, m, k)
            rhs = Fraction(0, 1)
            for j in range(m):
                r_ord = m - 1 - j
                def P_sum_r(S_arg: Tuple, rr=r_ord, kk=k) -> Fraction:
                    return sum((Pg_deriv(S_arg, rr, g, kk) for g in range(rr + 1)), Fraction(0, 1))
                rhs += I_deriv(t3, j, P_sum_r, ("P_sum", r_ord, k))
            if lhs != rhs:
                return False, f"Resolvent/generation expansion mismatch at k={k}, order m={m}: lhs={lhs} != rhs={rhs}"

    # Verify eq. (3) Laplace-transform identity explicitly as Laurent series in 1/z (large z):
    # Since L{t^m / m!}(z) = z^{-(m+1)}, exact equality of ALL coefficients m=0..4 of
    # M_t f_k and I * (R + R*B*R + R*B*R*B*R + ...) PROVES exact agreement of the
    # asymptotic series in 1/z of LHS and RHS of eq. (3) up to O(z^{-6})!
    return True, "Exact Taylor/generation expansion of asymptotics.md eq. (3) verified for k=2,3,4,5 up to order t^4/4! (matching O(z^{-5}))."


# ============================================================================
# 5. Fixed-size Yule Tauberian transfer (asymptotics.md lines 62-75)
# ============================================================================

def verify_tauberian_transfer(
    plant_gamma_mu_plus_1: bool = False,
) -> Tuple[bool, str]:
    """
    Test asymptotics.md lines 62-75 on synthetic sequences
      a_n = c * n^mu + B + D * n^{-alpha}
    (nondecreasing for n >= 3, a_2 = 0) across multiple (mu, c, B, D, alpha).
    Directly evaluates:
      1. Continuous-time Yule sum E(t) = q^2 sum_{n>=2} (n-1)(1-q)^{n-2} a_n,
         checking E(t) / e^{mu t} -> c * Gamma(mu + 2) (NOT Gamma(mu + 1)!).
      2. Measure nu_x Laplace transform int e^{-sy} d nu_x(y) -> d^* s^{-(mu+2)}
         with 1 - q = e^{-s/x}, checking the (1-q)^2 factor.
      3. Monotonicity sandwich bounds recovering c = d^* / Gamma(mu + 2).
    """
    cases = [
        (0.565723, 1.42, -0.3, 0.2, 0.5),
        (0.71, 0.85, 0.5, -0.1, 0.3),
        (0.45, 2.10, 0.0, 0.4, 0.8),
    ]
    gamma_shift = 1.0 if plant_gamma_mu_plus_1 else 2.0

    for mu, c, B, D, alpha in cases:
        def a_seq(n: int) -> float:
            if n <= 2:
                return 0.0
            return c * (n ** mu) + B + D * (n ** (-alpha))

        # Check monotonicity for n >= 3
        for n in range(3, 500):
            assert a_seq(n + 1) >= a_seq(n)

        d_expected = c * math.gamma(mu + gamma_shift)

        # 1. Direct Yule series at t = 6.5 (q = e^{-6.5} ~ 0.001503):
        t = 6.5
        q = math.exp(-t)
        n_cut = int(30.0 / q)  # (1-q)^{n_cut} < e^{-30}
        # Compute E(t) accurately via log-space or recurrence
        val = 0.0
        one_minus_q = 1.0 - q
        pw = 1.0
        for n in range(2, n_cut + 1):
            val += (n - 1) * pw * a_seq(n)
            pw *= one_minus_q
        val *= q * q
        ratio = val * math.exp(-mu * t)
        rel_err = abs(ratio - d_expected) / d_expected
        if rel_err > 1e-2:
            return False, f"Yule limit mismatch for mu={mu}: E(t)e^(-mu t)={ratio:.6f} vs c*Gamma(mu+{gamma_shift})={d_expected:.6f} (rel_err={rel_err:.4f})"

        # 2. Exact algebraic check of nu_x Laplace transform vs Yule series at 1-q = e^{-s/x}:
        x = 500.0
        s = 1.3
        q_s = 1.0 - math.exp(-s / x)
        n_cut2 = int(28.0 * x / s)
        lap_nu = 0.0
        yule_sum = 0.0
        step_e = math.exp(-s / x)
        pw_e = step_e * step_e  # e^{-2s/x}
        pw_1mq = 1.0            # (1-q)^{n-2}
        for n in range(2, n_cut2 + 1):
            an = a_seq(n)
            lap_nu += (n - 1) * an * pw_e
            yule_sum += (n - 1) * an * pw_1mq
            pw_e *= step_e
            pw_1mq *= step_e
        lap_nu *= x ** (-(mu + 2.0))
        yule_E = (q_s ** 2) * yule_sum
        exact_identity_rhs = (x ** (-(mu + 2.0))) * ((1.0 - q_s) ** 2 / (q_s ** 2)) * yule_E
        if abs(lap_nu - exact_identity_rhs) > 1e-12 * max(1.0, abs(lap_nu)):
            return False, "Exact algebraic identity between nu_x Laplace transform and Yule series failed!"

        target_lap = (c * math.gamma(mu + 2.0)) * (s ** (-(mu + 2.0)))
        if abs(lap_nu - target_lap) / target_lap > 1.5e-2:
            return False, f"nu_x Laplace limit mismatch: {lap_nu} vs {target_lap}"

        # 3. Sandwich recovery of c from d_true = c * Gamma(mu + 2):
        d_true = c * math.gamma(mu + 2.0)
        c_recovered = d_true / math.gamma(mu + gamma_shift)
        n0 = 15000
        h = 0.02
        m_lo = math.ceil((1.0 - h) * n0)
        m_hi = int((1.0 + h) * n0)
        upper = (sum((m - 1) * a_seq(m) for m in range(n0, m_hi + 1)) * (n0 ** (-(mu + 2.0)))) / (
            sum(m - 1 for m in range(n0, m_hi + 1)) / (n0 * n0)
        )
        lower = (sum((m - 1) * a_seq(m) for m in range(m_lo, n0 + 1)) * (n0 ** (-(mu + 2.0)))) / (
            sum(m - 1 for m in range(m_lo, n0 + 1)) / (n0 * n0)
        )
        if not (lower <= a_seq(n0) / (n0 ** mu) <= upper):
            return False, "Sandwich inequality failed!"
        if not (lower * 0.985 <= c_recovered <= upper * 1.015):
            return False, f"Recovered c={c_recovered:.6f} outside sandwich [{lower:.6f}, {upper:.6f}]"

    return True, "Fixed-size Yule Tauberian transfer verified on 3 synthetic power-law+correction sequences; confirms c_k = d_k^* / Gamma(mu+2)."


# ============================================================================
# 6. Analytic pole / residue / aperiodicity + degree-law leaf drift checks
# ============================================================================

def verify_spectral_pole_and_aperiodicity(
    plant_wrong_residue_sign: bool = False,
    plant_lattice_kernel: bool = False,
) -> Tuple[bool, str]:
    """
    Tests Section 2 of asymptotics.md on an explicit 2-mark renewal kernel
    (continuous density vs lattice planted error), checking:
      - unique real mu with r(K_mu)=1,
      - lambda'(mu) = ell(K'_mu g) < 0,
      - residue formula (4) vs contour integral around z = mu,
      - invertibility of I - K_{mu + i tau} for tau != 0.
    """
    a_rates = [0.5, 1.0 / 3.0]
    M = [[1.4, 0.6], [1.6, 0.7]]
    b_rates = [0.7, 1.1]

    def K_mat(z: complex):
        out = [[0j, 0j], [0j, 0j]]
        for i in range(2):
            for j in range(2):
                if plant_lattice_kernel:
                    # Pure lattice at t = 1 => K_{mu + 2 pi i} == K_mu!
                    out[i][j] = M[i][j] * math.exp(-a_rates[i] - 0.1 * j) * cmath.exp(-z)
                else:
                    out[i][j] = M[i][j] / (z + a_rates[i] + b_rates[j])
        return out

    def K_deriv(s: float):
        out = [[0.0, 0.0], [0.0, 0.0]]
        for i in range(2):
            for j in range(2):
                if plant_lattice_kernel:
                    out[i][j] = -M[i][j] * math.exp(-a_rates[i] - 0.1 * j) * math.exp(-s)
                else:
                    out[i][j] = -M[i][j] / ((s + a_rates[i] + b_rates[j]) ** 2)
        return out

    def perron(s: float):
        km = [[x.real for x in row] for row in K_mat(s)]
        tr = km[0][0] + km[1][1]
        det = km[0][0] * km[1][1] - km[0][1] * km[1][0]
        lam = 0.5 * (tr + math.sqrt(tr * tr - 4.0 * det))
        # Right eigenvector (K - lam I) g = 0
        g = [km[0][1], lam - km[0][0]]
        s_g = sum(g)
        g = [x / s_g for x in g]
        # Left eigenvector ell (K - lam I) = 0
        ell = [km[1][0], lam - km[0][0]]
        dot = ell[0] * g[0] + ell[1] * g[1]
        ell = [x / dot for x in ell]
        return lam, g, ell

    lo, hi = 0.0, 5.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if perron(mid)[0] > 1.0:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    lam, g, ell = perron(mu)
    assert abs(lam - 1.0) < 1e-13

    kp = K_deriv(mu)
    ell_kp_g = sum(ell[i] * kp[i][j] * g[j] for i in range(2) for j in range(2))
    eps = 1e-7
    fd = (perron(mu + eps)[0] - perron(mu - eps)[0]) / (2 * eps)
    if abs(ell_kp_g - fd) > 1e-8 or not (ell_kp_g < 0):
        return False, f"Eigenvalue derivative mismatch: {ell_kp_g} vs {fd}"

    def F_val(z: complex) -> complex:
        km = K_mat(z)
        a00, a01 = 1.0 - km[0][0], -km[0][1]
        a10, a11 = -km[1][0], 1.0 - km[1][1]
        det = a00 * a11 - a01 * a10
        jgf = [1.0 / (z + a_rates[0] + 0.9), 1.0 / (z + a_rates[1] + 0.9)]
        inv_jgf = [(a11 * jgf[0] - a01 * jgf[1]) / det, (-a10 * jgf[0] + a00 * jgf[1]) / det]
        iu = [0.8 / ((z + 0.3) * (z + a_rates[0] + 0.5)), 0.5 / ((z + 0.3) * (z + a_rates[1] + 0.5))]
        igf = 0.6 / ((z + 0.3) * (z + 1.2))
        return igf + iu[0] * inv_jgf[0] + iu[1] * inv_jgf[1]

    jgf_mu = [1.0 / (mu + a_rates[0] + 0.9), 1.0 / (mu + a_rates[1] + 0.9)]
    iu_mu = [0.8 / ((mu + 0.3) * (mu + a_rates[0] + 0.5)), 0.5 / ((mu + 0.3) * (mu + a_rates[1] + 0.5))]
    denom = ell_kp_g if plant_wrong_residue_sign else (-ell_kp_g)
    d_formula = ((iu_mu[0] * g[0] + iu_mu[1] * g[1]) * (ell[0] * jgf_mu[0] + ell[1] * jgf_mu[1])) / denom

    # Contour integral around z = mu
    rad = 1e-5
    n_pts = 64
    res_num = 0j
    for idx in range(n_pts):
        ang = 2.0 * math.pi * idx / n_pts
        dz = cmath.exp(1j * ang)
        res_num += F_val(mu + rad * dz) * rad * dz / n_pts

    if abs(res_num - d_formula) > 1e-10:
        return False, f"Residue mismatch: contour={res_num} vs formula={d_formula}"

    # Aperiodicity check on Re z = mu, tau != 0:
    for tau in (0.5, 1.0, math.pi, 2.0 * math.pi, 4.0 * math.pi):
        km = K_mat(complex(mu, tau))
        det = (1.0 - km[0][0]) * (1.0 - km[1][1]) - km[0][1] * km[1][0]
        if abs(det) < 1e-8:
            return False, f"Peripheral pole detected at tau={tau}: |det(I-K)|={abs(det)}"

    return True, "Spectral pole, lambda'(mu)=ell(K'_mu g)<0, residue (4), and aperiodicity verified."


def verify_degree_law_bounds() -> Tuple[bool, str]:
    """
    Checks two specific quantitative bounds in degree-law.md:
      1. Yule L^2 moment from Y(0)=2: E[(e^{-t} Y(t))^2] == 6 - 2 e^{-t} exactly.
      2. Leaf drift bound E[L(t) - 1] >= e^t (ell - 1) even when L(t) can visit 0
         (tested via exact Markov CTMC ODE integration for worst-case rho=0, c(0)=1, L(0)=3!).
    """
    # 1. Exact check of Yule second moment from k=2 using P(Y(t)=n) = (n-1) q^2 (1-q)^{n-2}:
    for t in (0.2, 0.7, 1.5, 3.0):
        q = math.exp(-t)
        n_max = int(35.0 / q) + 20
        m2 = 0.0
        pw = 1.0
        for n in range(2, n_max):
            m2 += (n * n) * (n - 1) * (q * q) * pw
            pw *= (1.0 - q)
        exact = 6.0 - 2.0 * q
        if abs(m2 * (q * q) - exact) > 1e-10:
            return False, f"Yule L^2 moment mismatch at t={t}: {m2*q*q} vs {exact}"

    # 2. Exact CTMC evolution of (L, c) starting at (L=3, c=1) with WORST-CASE rho = 0!
    # Transitions: (L, c) -> (L+1, c) at rate L; (L, c) -> (L-1, c+1) at rate L/(L+c) (0 if L=0).
    # Notice L=0 IS reachable in 3 promotions! Check E[L(t) - 1] >= e^t (3 - 1) = 2 e^t
    # AND check P(L(t) == 0) > 0 to confirm L=0 really is visited!
    L_cap = 90
    c_cap = 25
    # State index (L, c) for 0 <= L <= L_cap, 1 <= c <= c_cap
    def idx(L, c):
        return L * c_cap + (c - 1)

    N_st = (L_cap + 1) * c_cap
    prob = [0.0] * N_st
    prob[idx(3, 1)] = 1.0
    T_end = 1.2
    steps = 1200
    dt = T_end / steps

    def rhs(p):
        dp = [0.0] * N_st
        for L in range(L_cap + 1):
            for c in range(1, c_cap + 1):
                i = idx(L, c)
                pi = p[i]
                if pi == 0.0:
                    continue
                if L < L_cap and L > 0:
                    r_up = float(L)
                    dp[i] -= r_up * pi
                    dp[idx(L + 1, c)] += r_up * pi
                if L > 0 and c < c_cap:
                    r_dn = float(L) / float(L + c)
                    dp[i] -= r_dn * pi
                    dp[idx(L - 1, c + 1)] += r_dn * pi
        return dp

    for step in range(1, steps + 1):
        k1 = rhs(prob)
        p2 = [prob[i] + 0.5 * dt * k1[i] for i in range(N_st)]
        k2 = rhs(p2)
        p3 = [prob[i] + 0.5 * dt * k2[i] for i in range(N_st)]
        k3 = rhs(p3)
        p4 = [prob[i] + dt * k3[i] for i in range(N_st)]
        k4 = rhs(p4)
        prob = [prob[i] + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(N_st)]

    p_zero = sum(prob[idx(0, c)] for c in range(1, c_cap + 1))
    p_cap = sum(prob[idx(L_cap, c)] for c in range(1, c_cap + 1))
    mean_L = sum(L * prob[idx(L, c)] for L in range(L_cap + 1) for c in range(1, c_cap + 1))
    bound = 2.0 * math.exp(T_end)
    assert p_zero > 1e-4, f"L=0 should be visited with positive probability, got {p_zero}"
    assert p_cap < 1e-6, f"Truncation leaked mass {p_cap}"
    if not (mean_L - 1.0 >= bound - 1e-5):
        return False, f"Leaf drift bound failed: E[L(t)-1]={mean_L - 1.0} < {bound}"

    return True, f"Yule L^2 moment 6-2e^(-t) exact; worst-case rho=0 leaf drift at t={T_end}: E[L(t)-1]={mean_L-1:.6f} >= 2e^t={bound:.6f} (with P(L=0)={p_zero:.5f} > 0)."


# ============================================================================
# Main runner
# ============================================================================

def main():
    print("=== INDEPENDENT ADVERSARIAL AUDIT OF ASTRA CLAIM 3 ===")
    dists, E_ge, E_eq = exact_friend_tree_chain(n_max=10)

    print("\n[1] Exact small-n enumeration (n = 3..10):")
    print("  n | E[C_n] = E[N_{>=2}]          | E[N_{=2}]                    | E[N_{>=3}]")
    for n in range(3, 11):
        c_n = E_ge[(n, 2)]
        eq2 = E_eq[(n, 2)]
        ge3 = E_ge[(n, 3)]
        assert c_n == sum(E_eq[(n, k)] for k in range(2, n + 1))
        assert E_ge[(n, 2)] - E_ge[(n, 3)] == eq2
        print(f" {n:2d} | {str(c_n):28s} | {str(eq2):28s} | {ge3}")

    # Verify against Astra FACTS #2 values for n=3..9!
    expected_astra = [
        Fraction(1, 1), Fraction(4, 3), Fraction(19, 12), Fraction(647, 360),
        Fraction(17189, 8640), Fraction(4717423, 2177280), Fraction(101568829, 43545600)
    ]
    for i, exp_v in enumerate(expected_astra, start=3):
        assert E_ge[(i, 2)] == exp_v, f"Disagreement with Astra FACT at n={i}: {E_ge[(i,2)]} != {exp_v}"
    print("  PASS: Exact E[C_n] matches n=3..9 and exact telescoping sum_{k>=2} E[N_{=k}(n)] = E[C_n] holds.")

    checks = [
        ("Yule Laplace & T_2<->T_3 factor 2/(2+mu)", lambda: verify_yule_and_t2_t3_convolution(dists, E_ge, 10)),
        ("Crossing identity (degree-law.md eq. 35)", lambda: verify_crossing_identity(dists, E_ge, 10)),
        ("Small-t resolvent / generation expansion (eq. 3)", verify_small_t_resolvent_expansion),
        ("Fixed-size Yule Tauberian transfer (Gamma(mu+2))", verify_tauberian_transfer),
        ("Spectral pole, residue (4), aperiodicity", verify_spectral_pole_and_aperiodicity),
        ("Degree-law Yule L^2 & worst-case leaf drift", verify_degree_law_bounds),
    ]

    print("\n[2] Core verification checks:")
    for name, fn in checks:
        ok, msg = fn()
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}: {msg}")
        assert ok, f"Core check failed: {name}: {msg}"

    planted = [
        ("Planted Error 1: Gamma(mu+1) instead of Gamma(mu+2)",
         lambda: verify_tauberian_transfer(plant_gamma_mu_plus_1=True)),
        ("Planted Error 2: Wrong T_2<->T_3 factor 1/(1+mu) instead of 2/(2+mu)",
         lambda: verify_yule_and_t2_t3_convolution(dists, E_ge, 10, plant_wrong_t2_t3_factor=True)),
        ("Planted Error 3: Drop rho_o in crossing identity",
         lambda: verify_crossing_identity(dists, E_ge, 10, drop_rho=True)),
        ("Planted Error 4: Replace k-c_o by k-1 in crossing identity",
         lambda: verify_crossing_identity(dists, E_ge, 10, wrong_k_minus_1=True)),
        ("Planted Error 5: Replace (n-1) by n in Yule Laplace identity",
         lambda: verify_yule_and_t2_t3_convolution(dists, E_ge, 10, plant_wrong_yule_weight=True)),
        ("Planted Error 6: Wrong residue sign +ell(K'_mu g)",
         lambda: verify_spectral_pole_and_aperiodicity(plant_wrong_residue_sign=True)),
        ("Planted Error 7: Lattice kernel (peripheral pole at mu+2pi i)",
         lambda: verify_spectral_pole_and_aperiodicity(plant_lattice_kernel=True)),
    ]

    print("\n[3] Planted-error negative controls (all MUST be caught = rejected):")
    for name, fn in planted:
        ok, msg = fn()
        assert not ok, f"Planted error was NOT caught: {name}"
        print(f"  [CAUGHT] {name} -> {msg}")

    print("\nALL CLAIM 3 CHECKS PASSED & ALL 7 PLANTED ERRORS CAUGHT.")


if __name__ == "__main__":
    main()
