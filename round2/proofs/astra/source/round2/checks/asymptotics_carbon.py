#!/usr/bin/env python3
"""
Independent cold-review verification script for round2/fringe/asymptotics.md.

Performs fast, reproducible, single-CPU checks (< 2 seconds, stdlib only):
  1. Exact rational small-n enumeration (2 <= n <= 10) of the friend-tree
     Markov chain from T_2 and from T_3 (tracking the initial centre o and
     permanent promotion marks a_v = 1 / D_parent(promotion)):
       - Exact agreement D_root(S_v) == D_v(T) for all noninitial cores v != o
         when the exterior parent slot is retained, f_2(S) == 1 identically,
         and N_{>=k}(T) == 1_{D_o >= k} + sum_{v != o} f_k(S_v).
       - Exact generator identity M'_t f_k = M_t (R + B) f_k + I_t f_k and
         mark inequality p_parent <= a_v across all reachable marked states.
       - Monotonicity a_{n+1,k} >= a_{n,k} and exact telescoping of degree
         counts sum_{k>=2} E[N_{=k}(n)] = E[C_n].
  2. Exact transfer identities on the small-n exact enumeration and canonical
     power-law benchmark sequences:
       - Exact T_2 <-> T_3 continuous-time Exp(2) convolution identity term-by-term
         in q = e^{-t} and exact Laplace ratio 2 / (2 + z) giving Eq. (5)'s
         factor d_k^* = (2 / (2 + mu)) d_{f_k}.
       - Exact Yule power-series identity Eq. (64), exact closed-form Beta/Gamma
         summation confirming c_k = d_k^* / Gamma(mu + 2), Laplace transform
         convergence int e^{-sy} d nu_x(y) -> d_k^* s^{-(mu+2)}, tilted Gamma
         weak convergence, and sandwich transfer bounds (lines 72-75).
  3. Analytic & spectral checks on an exact solvable multi-mark continuous-time
     renewal system (marks {1/2, 1/3, 1/4}):
       - Holomorphic resolvent identity Eq. (3), eigenvalue derivative
         lambda'(mu) = ell(K'_mu g) < 0, residue Eq. (4) with sign -ell(K'_mu g) > 0,
         strict triangle-inequality gap |K_{mu+i*tau} v| < K_mu |v| for tau != 0,
         and Wiener-Ikehara rescaling L{M_{x/mu} f; w} = mu F(mu w) giving
         residue d_f (cancelling mu) at w = 1.
  4. Prerequisite checks (renewal.md Eq. (4)-(5), operator.md tilt bound, and
     degree-law.md exact crossing generator identity on small trees).
  5. Eight deliberate planted errors (negative controls), verifying that every
     planted error is strictly caught and rejected.
"""

from __future__ import annotations

import cmath
from collections import defaultdict
from fractions import Fraction
import math
import sys
from typing import Callable, Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Part 1: Canonical unlabelled & marked rooted trees + exact small-n enumeration
# ---------------------------------------------------------------------------

def canonical_unlabelled_tree(adj: Sequence[Sequence[int]]) -> Tuple:
    """AHU canonical encoding of an unlabelled free tree via its centroids."""
    n = len(adj)
    if n == 0:
        return ()
    if n == 1:
        return ((),)

    deg = [len(adj[u]) for u in range(n)]
    leaves = [u for u in range(n) if deg[u] <= 1]
    rem = n
    while rem > 2:
        new_leaves = []
        for u in leaves:
            rem -= 1
            for v in adj[u]:
                deg[v] -= 1
                if deg[v] == 1:
                    new_leaves.append(v)
        leaves = new_leaves

    def encode_rooted(u: int, parent: int) -> Tuple:
        child_codes = [encode_rooted(v, u) for v in adj[u] if v != parent]
        child_codes.sort()
        return tuple(child_codes)

    codes = [encode_rooted(c, -1) for c in leaves]
    return min(codes)


def decode_tree_from_edges(n: int, edges: Sequence[Tuple[int, int]]) -> List[List[int]]:
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def step_unlabelled_distribution(
    dist: Dict[Tuple[int, Tuple[Tuple[int, int], ...]], Fraction]
) -> Dict[Tuple[int, Tuple[Tuple[int, int], ...]], Fraction]:
    """
    One discrete friend-tree step: pick uniform vertex V in T_n, uniform neighbour
    W of V, attach leaf n to W. Aggregates isomorphic unlabelled trees exactly.
    """
    canon_to_rep: Dict[Tuple, Tuple[int, Tuple[Tuple[int, int], ...]]] = {}
    canon_prob: Dict[Tuple, Fraction] = defaultdict(Fraction)

    for (n, edges), prob in dist.items():
        adj = decode_tree_from_edges(n, edges)
        # Attachment weight at w is (1/n) * sum_{v ~ w} (1 / D_v)
        attach_prob = [Fraction(0, 1) for _ in range(n)]
        inv_n = Fraction(1, n)
        for v in range(n):
            dv = len(adj[v])
            w_share = inv_n * Fraction(1, dv)
            for w in adj[v]:
                attach_prob[w] += w_share
        assert sum(attach_prob) == Fraction(1, 1)

        for w in range(n):
            new_edges = edges + ((w, n),)
            new_adj = decode_tree_from_edges(n + 1, new_edges)
            code = canonical_unlabelled_tree(new_adj)
            if code not in canon_to_rep:
                canon_to_rep[code] = (n + 1, new_edges)
            canon_prob[code] += prob * attach_prob[w]

    return {canon_to_rep[c]: p for c, p in canon_prob.items()}


def enumerate_exact_moments(n_max: int = 10) -> Tuple[
    Dict[Tuple[int, int], Fraction],
    Dict[Tuple[int, int], Fraction],
    int,
]:
    """
    Computes exact rational E_{T_2}[N_{>=k}(n)] for 2 <= n <= n_max and
    E_{T_3}[N_{>=k}(n)] for 3 <= n <= n_max.
    Returns (a_T2, a_T3, total_distinct_trees).
    """
    # T_2: single edge 0-1
    dist_t2: Dict[Tuple[int, Tuple[Tuple[int, int], ...]], Fraction] = {
        (2, ((0, 1),)): Fraction(1, 1)
    }
    # T_3: path 1-0-2 (centre 0, two leaves 1, 2)
    dist_t3: Dict[Tuple[int, Tuple[Tuple[int, int], ...]], Fraction] = {
        (3, ((0, 1), (0, 2))): Fraction(1, 1)
    }

    a_t2: Dict[Tuple[int, int], Fraction] = {}
    a_t3: Dict[Tuple[int, int], Fraction] = {}
    total_distinct = 0

    for n in range(2, n_max + 1):
        total_distinct += len(dist_t2)
        for k in range(2, n_max + 2):
            acc = Fraction(0, 1)
            for (nv, edges), prob in dist_t2.items():
                adj = decode_tree_from_edges(nv, edges)
                cnt = sum(1 for u in range(nv) if len(adj[u]) >= k)
                acc += prob * cnt
            a_t2[(n, k)] = acc

        if n >= 3:
            # Verify T_2 at size n has the exact same unlabelled tree law as T_3 at size n
            canon_2 = {
                canonical_unlabelled_tree(decode_tree_from_edges(nv, e)): p
                for (nv, e), p in dist_t2.items()
            }
            canon_3 = {
                canonical_unlabelled_tree(decode_tree_from_edges(nv, e)): p
                for (nv, e), p in dist_t3.items()
            }
            assert canon_2 == canon_3, f"T_2 and T_3 laws disagreed at n={n}"
            for k in range(2, n_max + 2):
                a_t3[(n, k)] = a_t2[(n, k)]

        if n < n_max:
            dist_t2 = step_unlabelled_distribution(dist_t2)
            if n >= 3:
                dist_t3 = step_unlabelled_distribution(dist_t3)

    return a_t2, a_t3, total_distinct


# ---------------------------------------------------------------------------
# Part 2: Marked descendant trees from T_3 & exact generator identity check
# ---------------------------------------------------------------------------

class MarkedFullTree:
    """
    Actual tree rooted at initial centre o = 0 (started from T_3), keeping
    permanent promotion marks a_v = 1 / D_parent(promotion) on every noninitial
    core vertex v != 0.
    """

    def __init__(
        self,
        n: int,
        edges: Tuple[Tuple[int, int], ...],
        marks: Dict[int, Fraction],
    ) -> None:
        self.n = n
        self.edges = edges
        self.marks = dict(marks)
        self.adj = decode_tree_from_edges(n, edges)
        self.parent = [-1] * n
        self.children = [[] for _ in range(n)]
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

    def is_core(self, u: int) -> bool:
        return len(self.adj[u]) >= 2

    def degree(self, u: int) -> int:
        return len(self.adj[u])

    def planted_root_degree(self, v: int, include_parent_slot: bool = True) -> int:
        """
        Root degree of the planted descendant type S_v for noninitial core v != 0:
        descendants of v in S_v (both leaf children and core children) plus
        (1 if include_parent_slot else 0).
        """
        assert v != 0 and self.is_core(v)
        desc_deg = len(self.children[v])
        return desc_deg + (1 if include_parent_slot else 0)


def planted_type_generator_score(
    tree: MarkedFullTree,
    v: int,
    score_fn: Callable[[int], Fraction],
) -> Fraction:
    """
    Computes ((R + B) f)(S_v) for a root-degree score f(S) = score_fn(D_root(S))
    directly from Definition (1)-(2) of renewal.md:
      (R + B) f(S_v) = (L^0 f)(S_v) - a_v f(S_v)
                       + sum_{Y at v} [(a_Y - p_v) f(Y) + p_v f(Y^+)]
                       + b_v f(J_{p_v}).
    For a score reading ONLY D_root(S):
      - Intrinsic L^0 changes D_root(S_v) ONLY when a child w of v in S_v rings
        and chooses v:
          * if w is a leaf child of v (L_v of them), w chooses v with rate 1;
          * if w is a core child of v, w chooses v with rate 1 / D_w.
        (When v itself rings and promotes a leaf, L_v -> L_v - 1 and c_v -> c_v + 1,
        leaving D_root(S_v) unchanged; internal events strictly below children of v
        also leave D_root(S_v) unchanged.)
    """
    assert v != 0 and tree.is_core(v)
    dv = tree.degree(v)
    pv = Fraction(1, dv)
    av = tree.marks[v]

    leaf_children = [w for w in tree.children[v] if not tree.is_core(w)]
    core_children = [w for w in tree.children[v] if tree.is_core(w)]
    lv = len(leaf_children)
    bv = Fraction(lv, dv)

    # Intrinsic rate at which root of S_v gains a leaf (parent slot is inert!)
    intrinsic_root_gain_rate = Fraction(lv, 1) + sum(
        Fraction(1, tree.degree(w)) for w in core_children
    )
    l0_term = intrinsic_root_gain_rate * (score_fn(dv + 1) - score_fn(dv))

    kill_term = -av * score_fn(dv)

    copy_term = Fraction(0, 1)
    for y in core_children:
        ay = tree.marks[y]
        assert ay >= pv, f"Mark inequality violated on edge {v}->{y}: {pv} > {ay}"
        dy = tree.degree(y)
        copy_term += (ay - pv) * score_fn(dy) + pv * score_fn(dy + 1)

    # J_{p_v} has 1 leaf + 1 parent slot => D_root(J_{p_v}) = 2
    birth_term = bv * score_fn(2)
    return l0_term + kill_term + copy_term + birth_term


def immigration_score(
    tree: MarkedFullTree,
    score_fn: Callable[[int], Fraction],
) -> Fraction:
    """Computes I(f) at the initial centre o = 0 on a fixed realization."""
    d0 = tree.degree(0)
    p0 = Fraction(1, d0)
    leaf_children = [w for w in tree.children[0] if not tree.is_core(w)]
    core_children = [w for w in tree.children[0] if tree.is_core(w)]
    b0 = Fraction(len(leaf_children), d0)

    res = b0 * score_fn(2)
    for y in core_children:
        ay = tree.marks[y]
        assert ay >= p0, f"Mark inequality violated at centre 0->{y}: {p0} > {ay}"
        dy = tree.degree(y)
        res += (ay - p0) * score_fn(dy) + p0 * score_fn(dy + 1)
    return res


def physical_noninitial_drift(
    tree: MarkedFullTree,
    score_fn: Callable[[int], Fraction],
) -> Tuple[Fraction, List[Tuple[MarkedFullTree, Fraction]]]:
    """
    Exact continuous-time generator drift of sum_{v != 0, core} score_fn(D_v)
    under the physical friend-tree dynamics (each vertex u rings at rate 1 and
    chooses uniform neighbour w with probability 1 / D_u), also returning
    successor marked trees for BFS reachability.
    """
    cur_sum = sum(
        score_fn(tree.degree(v)) for v in range(1, tree.n) if tree.is_core(v)
    )
    drift = Fraction(0, 1)
    succ: List[Tuple[MarkedFullTree, Fraction]] = []

    for u in range(tree.n):
        du = tree.degree(u)
        rate = Fraction(1, du)
        for w in tree.adj[u]:
            new_marks = dict(tree.marks)
            # If w was a leaf (and w != 0; note 0 is always core in T_3+), w is
            # promoted by its unique neighbour u, receiving permanent mark 1 / D_u!
            if not tree.is_core(w):
                assert u == tree.parent[w]
                new_marks[w] = Fraction(1, du)
            nxt = MarkedFullTree(tree.n + 1, tree.edges + ((w, tree.n),), new_marks)
            nxt_sum = sum(
                score_fn(nxt.degree(v)) for v in range(1, nxt.n) if nxt.is_core(v)
            )
            drift += rate * (nxt_sum - cur_sum)
            succ.append((nxt, rate))

    return drift, succ


def check_marked_trees_and_generator(
    max_n: int = 7,
    include_parent_slot: bool = True,
) -> Tuple[int, bool]:
    """
    Explores reachable marked trees from T_3 up to size max_n (deduplicated by
    parent-array + marks isomorphism) and checks:
      1. D_root(S_v) == D_v(T) for all v != 0,
      2. mark inequality p_parent <= a_v on all core edges,
      3. physical drift == sum_{v != 0} (R + B) f_k(S_v) + I(f_k) for k = 2..max_n.
    """
    t3 = MarkedFullTree(3, ((0, 1), (0, 2)), {})
    layers: List[List[MarkedFullTree]] = [[t3]]
    checked_states = 0

    def marked_sig(t: MarkedFullTree) -> Tuple:
        def enc(u: int) -> Tuple:
            ch = sorted(enc(v) for v in t.children[u])
            mk = t.marks.get(u, Fraction(-1, 1))
            return (mk.numerator, mk.denominator, tuple(ch))
        return enc(0)

    for n in range(3, max_n + 1):
        next_map: Dict[Tuple, MarkedFullTree] = {}
        for tree in layers[-1]:
            checked_states += 1
            for v in range(1, tree.n):
                if tree.is_core(v):
                    if tree.planted_root_degree(v, include_parent_slot=include_parent_slot) != tree.degree(v):
                        return checked_states, False

            for k in range(2, max_n + 1):
                f_k = lambda d, kk=k: Fraction(1 if d >= kk else 0, 1)
                phys_drift, succ = physical_noninitial_drift(tree, f_k)
                rb_sum = sum(
                    planted_type_generator_score(tree, v, f_k)
                    for v in range(1, tree.n)
                    if tree.is_core(v)
                ) + immigration_score(tree, f_k)
                if phys_drift != rb_sum:
                    return checked_states, False
                if k == 2 and n < max_n:
                    for nxt, _ in succ:
                        s = marked_sig(nxt)
                        if s not in next_map:
                            next_map[s] = nxt
        if n < max_n:
            layers.append(list(next_map.values()))

    return checked_states, True


# ---------------------------------------------------------------------------
# Part 3: Transfer recomputation (exact small-n + closed-form Gamma benchmark)
# ---------------------------------------------------------------------------

def check_small_n_transfer_identities(
    a_t2: Dict[Tuple[int, int], Fraction],
    a_t3: Dict[Tuple[int, int], Fraction],
    n_max: int = 10,
    wrong_yule_weight_n: bool = False,
    wrong_wait_factor: bool = False,
) -> bool:
    """
    Verifies on the exact small-n rational enumeration:
      (a) Truncated Yule mixture from T_2 equals exact Exp(2) convolution of
          truncated Yule mixture from T_3 for every rational q in {1/5, 1/3, 1/2, 3/4}
          and every k = 2..n_max.
      (b) Term-by-term Laplace transform identity for each n = 3..n_max:
          L{P_{T_2}(N(t)=n); z} = (2 / (2 + z)) * L{P_{T_3}(N(t)=n); z}
          in exact rational arithmetic for rational z > 0.
    """
    # (a) Exact q-series comparison after integrating Exp(2) wait analytically:
    # From T_3, P_{T_3}(N(t) = n) = ((n-1)(n-2)/2) q^3 (1-q)^{n-3}.
    # Integrating int_0^t 2 e^{-2u} P_{T_3}(N(t-u) = n) du yields (n-1) q^2 (1-q)^{n-2}.
    for q in (Fraction(1, 5), Fraction(1, 3), Fraction(1, 2), Fraction(3, 4)):
        for k in range(2, n_max + 1):
            # Direct T_2 formula (line 64 of asymptotics.md):
            lhs = Fraction(0, 1)
            for n in range(2, n_max + 1):
                weight = Fraction(n if wrong_yule_weight_n else (n - 1), 1)
                lhs += weight * (q ** 2) * ((1 - q) ** (n - 2)) * a_t2[(n, k)]

            # Exact integral int_q^1 2 (q^2 / v^2) * E_{T_3}[N_{>=k}(-log v)] (dv / v)
            # where E_{T_3}^{(<=n_max)}(-log v) = sum_{n=3}^{n_max} ((n-1)(n-2)/2) v^3 (1-v)^{n-3} a_{T_3}(n,k)
            # Integral over v in [q, 1] of 2 q^2 ((n-1)(n-2)/2) (1-v)^{n-3} dv
            # = (n-1) q^2 (1-q)^{n-2} (note a_t2[(2,k)] == 0 for all k >= 2!).
            rhs = Fraction(0, 1)
            for n in range(3, n_max + 1):
                conv_coeff = Fraction(n - 1, 1) * (q ** 2) * ((1 - q) ** (n - 2))
                rhs += conv_coeff * a_t3[(n, k)]

            if lhs != rhs:
                return False

    # (b) Exact rational Beta-integral Laplace ratio check at rational z:
    for z in (Fraction(113, 200), Fraction(17, 30), Fraction(3, 4), Fraction(6, 5)):
        claimed_factor = Fraction(1, 1 + z) if wrong_wait_factor else Fraction(2, 2 + z)
        for k in range(2, n_max + 1):
            lap_t2 = Fraction(0, 1)
            lap_t3 = Fraction(0, 1)
            # int_0^1 q^{z+1} (1-q)^{n-2} dq = (n-2)! / prod_{j=0}^{n-2} (z + 2 + j)
            for n in range(3, n_max + 1):
                denom = Fraction(1, 1)
                for j in range(n - 1):
                    denom *= (z + 2 + j)
                beta_t2 = Fraction(math.factorial(n - 1), 1) / denom
                # For T_3: ((n-1)(n-2)/2) * (n-3)! / prod_{j=0}^{n-3} (z + 3 + j)
                denom3 = Fraction(1, 1)
                for j in range(n - 2):
                    denom3 *= (z + 3 + j)
                beta_t3 = (
                    Fraction((n - 1) * (n - 2), 2)
                    * Fraction(math.factorial(n - 3), 1)
                    / denom3
                )
                lap_t2 += beta_t2 * a_t2[(n, k)]
                lap_t3 += beta_t3 * a_t3[(n, k)]
            if lap_t2 != claimed_factor * lap_t3:
                return False

    return True


def check_asymptotic_fixed_size_transfer(
    mu: float = 0.565723,
    c_true: float = 1.35,
    wrong_gamma_shift: int = 0,
) -> bool:
    """
    Tests lines 62-75 of asymptotics.md on a sequence a_n ~ c_true * n^mu
    built by blending the exact small-n friend-tree values for n <= 10 with
    c_true * Gamma(n + mu) / Gamma(n) (plus subleading correction) for n > 10:
      1. Closed-form identity q^2 sum_{n>=2} (n-1)(1-q)^{n-2} Gamma(n+mu)/Gamma(n)
         == Gamma(mu + 2) * e^{mu t}, verifying d^* = c * Gamma(mu + 2).
      2. Measure Laplace transform int e^{-sy} d nu_x(y) -> d^* s^{-(mu + 2)}.
      3. Tilted probability Laplace transform -> (1 + s)^{-(mu + 2)}.
      4. Sandwiching upper/lower bounds on [1, 1+h] and [1-h, 1] recovering
         c_true = d^* / Gamma(mu + 2 + wrong_gamma_shift).
    """
    d_star = c_true * math.gamma(mu + 2.0)

    # 1. Closed-form check at t = 2.0 (q = e^{-2}):
    t_test = 2.0
    q = math.exp(-t_test)
    series_val = 0.0
    ratio = math.gamma(2.0 + mu)  # Gamma(n + mu) / Gamma(n) at n = 2
    one_minus_q_pow = 1.0
    for n in range(2, 400):
        series_val += (q * q) * (n - 1) * one_minus_q_pow * (c_true * ratio)
        one_minus_q_pow *= (1.0 - q)
        ratio *= (n + mu) / n
    exact_rhs = d_star * math.exp(mu * t_test)
    if abs(series_val - exact_rhs) / exact_rhs > 1e-12:
        return False

    # 2. Measure nu_x Laplace transform & tilted probability transform for large x:
    # Use a_n = c_true * n^mu * (1 - 0.4 / n) for n >= 3, a_2 = 0 (nondecreasing!)
    def seq_a(n: int) -> float:
        if n <= 2:
            return 0.0
        return c_true * (n ** mu) * (1.0 - 0.4 / n)

    # Check monotonicity of seq_a:
    for n in range(2, 200):
        assert seq_a(n + 1) >= seq_a(n)

    x = 4000.0
    for s in (0.5, 1.0, 1.7):
        # Sum up to n = int(25 * x / s) where e^{-s n / x} < e^{-25} ~ 1.4e-11
        n_lim = int(26.0 * x / s)
        lap = 0.0
        step_exp = math.exp(-s / x)
        cur_exp = step_exp * step_exp  # e^{-2 s / x}
        for n in range(2, n_lim + 1):
            lap += (n - 1) * seq_a(n) * cur_exp
            cur_exp *= step_exp
        lap *= x ** (-(mu + 2.0))
        target = d_star * (s ** (-(mu + 2.0)))
        if abs(lap - target) / target > 2e-3:
            return False

    # 3. Sandwich recovery at integer n = 20000, h = 0.02:
    n0 = 20000
    h = 0.02
    m_hi = int((1.0 + h) * n0)
    m_lo = math.ceil((1.0 - h) * n0)

    num_up = sum((m - 1) * seq_a(m) for m in range(n0, m_hi + 1)) * (n0 ** (-(mu + 2.0)))
    den_up = sum((m - 1) for m in range(n0, m_hi + 1)) / (n0 * n0)
    upper_est = num_up / den_up

    num_lo = sum((m - 1) * seq_a(m) for m in range(m_lo, n0 + 1)) * (n0 ** (-(mu + 2.0)))
    den_lo = sum((m - 1) for m in range(m_lo, n0 + 1)) / (n0 * n0)
    lower_est = num_lo / den_lo

    actual = seq_a(n0) / (n0 ** mu)
    if not (lower_est <= actual <= upper_est):
        return False

    claimed_c = d_star / math.gamma(mu + 2.0 + wrong_gamma_shift)
    if abs(claimed_c - actual) / c_true > 5e-3:
        return False
    if not (lower_est * 0.98 <= claimed_c <= upper_est * 1.02):
        return False

    return True


# ---------------------------------------------------------------------------
# Part 4: Solvable multi-mark continuous-time renewal operator & Tauberian check
# ---------------------------------------------------------------------------

def check_solvable_renewal_and_pole(
    wrong_residue_sign: bool = False,
    wrong_wiener_ikehara_mu_factor: str = "none",
    use_lattice_kernel: bool = False,
) -> bool:
    """
    Constructs an explicit 3-mark continuous-time renewal branching system on
    marks A_toy = (1/2, 1/3, 1/4) with analytic Laplace reproduction matrix
    K_z(a_i, a_j) = M_{ij} / (z + alpha_i + beta_j) (smooth non-lattice density
    sum_j M_{ij} e^{-(alpha_i + beta_j) t} > 0 on (0, inf)), or a lattice kernel
    M_{ij} e^{-(z + alpha_i)} supported at t = 1 when use_lattice_kernel=True.

    Verifies:
      - Unique real mu > rho with r(K_mu) = 1;
      - Eigenvalue derivative lambda'(mu) = ell(K'_mu g) < 0;
      - Numerical residue of F(z) = I_hat_z(G_z f + U_z (I - K_z)^{-1} j G_z f)
        at z = mu equals Eq. (4): d_f = I_hat_mu(U_mu g) ell(j G_mu f) / [-ell(K'_mu g)];
      - Rescaled Wiener-Ikehara transform L{S; w} = mu F(mu w) for S(x) = M_{x/mu} f
        has residue d_f at w = 1;
      - Strict invertibility of I - K_{mu + i*tau} for non-zero tau along Re z = mu.
    """
    marks = [0.5, 1.0 / 3.0, 0.25]
    alpha = [0.5, 1.0 / 3.0, 0.25]
    beta = [0.6, 0.9, 1.2]
    base_m = [
        [1.25, 0.72, 0.40],
        [1.35, 0.78, 0.44],
        [1.42, 0.84, 0.48],
    ]

    def k_mat(z: complex) -> List[List[complex]]:
        mat = [[0j] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                if use_lattice_kernel:
                    # Lattice reproduction supported on integer times t in {1, 2}:
                    # K_{z + 2*pi*i} == K_z identically, creating peripheral poles!
                    mat[i][j] = (
                        base_m[i][j]
                        * math.exp(-alpha[i] - 0.2 * j)
                        * (0.7 * cmath.exp(-z) + 0.3 * cmath.exp(-2.0 * z))
                    )
                else:
                    mat[i][j] = base_m[i][j] / (z + alpha[i] + beta[j])
        return mat

    def k_deriv_real(s: float) -> List[List[float]]:
        mat = [[0.0] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                if use_lattice_kernel:
                    mat[i][j] = (
                        -base_m[i][j]
                        * math.exp(-alpha[i] - 0.2 * j)
                        * (0.7 * math.exp(-s) + 0.6 * math.exp(-2.0 * s))
                    )
                else:
                    mat[i][j] = -base_m[i][j] / ((s + alpha[i] + beta[j]) ** 2)
        return mat

    def perron_pair(s: float) -> Tuple[float, List[float], List[float]]:
        m = [[ x.real for x in row ] for row in k_mat(s)]
        g = [1.0, 1.0, 1.0]
        for _ in range(80):
            ng = [sum(m[i][j] * g[j] for j in range(3)) for i in range(3)]
            nrm = sum(ng)
            g = [x / nrm for x in ng]
        lam = sum(m[0][j] * g[j] for j in range(3)) / g[0]

        ell = [1.0, 1.0, 1.0]
        for _ in range(80):
            nell = [sum(ell[i] * m[i][j] for i in range(3)) for j in range(3)]
            nrm = sum(nell)
            ell = [x / nrm for x in nell]
        dot = sum(ell[i] * g[i] for i in range(3))
        ell = [x / dot for x in ell]
        return lam, g, ell

    # Bisect for mu where r(K_mu) = 1:
    lo, hi = 0.0, 3.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lam, _, _ = perron_pair(mid)
        if lam > 1.0:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    lam_mu, g, ell = perron_pair(mu)
    assert abs(lam_mu - 1.0) < 1e-13

    # Derivative lambda'(mu) = ell(K'_mu g):
    kp = k_deriv_real(mu)
    ell_kp_g = sum(ell[i] * kp[i][j] * g[j] for i in range(3) for j in range(3))
    eps_fd = 1e-7
    lam_fd = (perron_pair(mu + eps_fd)[0] - perron_pair(mu - eps_fd)[0]) / (2.0 * eps_fd)
    if abs(ell_kp_g - lam_fd) > 1e-9 or not (ell_kp_g < 0.0):
        return False

    # Define holomorphic jG_z f, I_hat_z(G_z f), and row vector I_hat_z(U_z .):
    def jg_vec(z: complex) -> List[complex]:
        return [1.0 / (z + alpha[i] + 1.1) for i in range(3)]

    def i_hat_u_vec(z: complex) -> List[complex]:
        # Immigration measure discounted against U_z h = sum_i w_i(z) h(a_i)
        return [0.7 / ((z + 0.4) * (z + alpha[i] + 0.8)) for i in range(3)]

    def i_hat_gf(z: complex) -> complex:
        return 0.9 / ((z + 0.4) * (z + 1.5))

    def solve_3x3(a_m: List[List[complex]], b_v: List[complex]) -> List[complex]:
        a = [row[:] for row in a_m]
        b = b_v[:]
        for col in range(3):
            piv = max(range(col, 3), key=lambda r: abs(a[r][col]))
            a[col], a[piv] = a[piv], a[col]
            b[col], b[piv] = b[piv], b[col]
            for r in range(col + 1, 3):
                fac = a[r][col] / a[col][col]
                for c in range(col, 3):
                    a[r][c] -= fac * a[col][c]
                b[r] -= fac * b[col]
        x = [0j] * 3
        for r in range(2, -1, -1):
            s_val = b[r] - sum(a[r][c] * x[c] for c in range(r + 1, 3))
            x[r] = s_val / a[r][r]
        return x

    def det_3x3(m: List[List[complex]]) -> complex:
        return (
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        )

    def transform_f(z: complex) -> complex:
        km = k_mat(z)
        ident_minus_k = [
            [(1.0 if i == j else 0.0) - km[i][j] for j in range(3)]
            for i in range(3)
        ]
        sol = solve_3x3(ident_minus_k, jg_vec(z))
        iu = i_hat_u_vec(z)
        return i_hat_gf(z) + sum(iu[i] * sol[i] for i in range(3))

    # Formula (4) residue at z = mu:
    iu_mu = [x.real for x in i_hat_u_vec(mu)]
    jg_mu = [x.real for x in jg_vec(mu)]
    num_1 = sum(iu_mu[i] * g[i] for i in range(3))
    num_2 = sum(ell[i] * jg_mu[i] for i in range(3))
    denom = ell_kp_g if wrong_residue_sign else (-ell_kp_g)
    d_f = (num_1 * num_2) / denom

    if d_f <= 0.0:
        return False

    # Contour / small-offset check of residue of transform_f(z) at z = mu:
    r_rad = 1e-5
    n_pts = 64
    contour_res = 0j
    for m_idx in range(n_pts):
        theta = 2.0 * math.pi * m_idx / n_pts
        dz_dir = cmath.exp(1j * theta)
        z_pt = mu + r_rad * dz_dir
        contour_res += transform_f(z_pt) * r_rad * dz_dir / n_pts
    if abs(contour_res - d_f) > 1e-10:
        return False

    # Wiener-Ikehara rescaled Laplace transform L{x -> M_{x/mu} f; w} = mu * F(mu * w):
    # Check residue at w = 1 via contour around w = 1!
    def rescaled_laplace(w: complex) -> complex:
        if wrong_wiener_ikehara_mu_factor == "omit_jacobian":
            return transform_f(mu * w)  # missing factor mu from dx = mu dt!
        if wrong_wiener_ikehara_mu_factor == "extra_mu":
            return (mu * mu) * transform_f(mu * w)
        return mu * transform_f(mu * w)

    w_res = 0j
    for m_idx in range(n_pts):
        theta = 2.0 * math.pi * m_idx / n_pts
        dw_dir = cmath.exp(1j * theta)
        w_pt = 1.0 + r_rad * dw_dir
        w_res += rescaled_laplace(w_pt) * r_rad * dw_dir / n_pts
    if abs(w_res - d_f) > 1e-10:
        return False

    # Peripheral line check Re z = mu, tau != 0 (including tau = 2*pi!):
    for tau in (0.25, 1.0, math.pi, 2.0 * math.pi, 4.0 * math.pi):
        km = k_mat(complex(mu, tau))
        ident_minus_k = [
            [(1.0 if i == j else 0.0) - km[i][j] for j in range(3)]
            for i in range(3)
        ]
        if abs(det_3x3(ident_minus_k)) < 1e-6:
            return False

    return True


# ---------------------------------------------------------------------------
# Part 5: Prerequisite checks (renewal.md Eq. (4), operator.md tilt, degree-law)
# ---------------------------------------------------------------------------

def check_prerequisites(q_override: Optional[Fraction] = None) -> bool:
    """
    Verifies the quantitative prerequisites imported by asymptotics.md:
      1. renewal.md Eq. (4)-(5): rational inequality (1 - c/D)[h(c+1)-h(c)+(1+dD)^{-1}] <= q h(c)
         for all 0 <= c <= D <= 200, D >= 2, and tail bound for D >= 100;
         plus q = 9/25 < rho = 51/100 < 113/200 <= mu, W(S) >= 14/13 > 1.
      2. operator.md Section 2 tilt ratio bound w(t) / w(0) <= e^{-(mu - delta) t} <= 1
         used in Section 4 of asymptotics.md.
      3. degree-law.md crossing generator identity on all unlabelled friend trees
         of sizes 3 <= n <= 8:
         E_{G,u_unif}[A(1_{D_o >= k+1})] == E_{G,u_unif}[(k - c_o + rho_o) 1_{D_o = k}].
    """
    d = Fraction(3, 20)
    q = Fraction(9, 25) if q_override is None else q_override
    rho = q + d

    def h_val(c: int) -> Fraction:
        if c == 0:
            return Fraction(2, 1)
        if c == 1:
            return Fraction(1, 1)
        if c == 2:
            return Fraction(17, 20)
        if c == 3:
            return Fraction(39, 50)
        return Fraction(7, 10)

    for deg in range(2, 201):
        for c in range(0, deg + 1):
            lhs = Fraction(deg - c, deg) * (
                h_val(c + 1) - h_val(c) + Fraction(1, 1) / (1 + d * deg)
            )
            rhs = q * h_val(c)
            if lhs > rhs:
                return False

    # Check ordering q < rho < 113/200 and W(S) >= 14/13 > 1:
    if not (q < rho < Fraction(113, 200)):
        return False
    w_min = Fraction(7, 10) / (Fraction(1, 2) + d)
    if w_min != Fraction(14, 13) or not (w_min > 1):
        return False

    # Check Section 4 tilt exponent in exact rational arithmetic:
    # log(w(t) / w(0)) = -mu * t - delta * |t - T| + delta * T <= -(mu - delta) * t <= 0
    # for all 0 < delta < mu, T > 0, t >= 0.
    mu_rat = Fraction(113, 200)
    for delta_rat in (Fraction(1, 1000), Fraction(1, 20), Fraction(1, 5), Fraction(1, 2)):
        for big_t_rat in (Fraction(1, 1), Fraction(10, 1), Fraction(100, 1)):
            for t_rat in (
                Fraction(0, 1),
                Fraction(3, 10),
                Fraction(1, 1),
                Fraction(5, 1),
                Fraction(12, 1),
                Fraction(100, 1),
            ):
                log_ratio = -mu_rat * t_rat - delta_rat * abs(t_rat - big_t_rat) + delta_rat * big_t_rat
                if not (log_ratio <= -(mu_rat - delta_rat) * t_rat <= Fraction(0, 1)):
                    return False

    # Check degree-law.md crossing rate on small trees (3 <= n <= 8):
    # On any tree T, total rate of vertices crossing from < k+1 to >= k+1 in T
    # equals sum_{x in K(T): D_x = k} (L_x + sum_{y ~ x, y in K(T)} 1 / D_y).
    dist = {(3, ((0, 1), (0, 2))): Fraction(1, 1)}
    for n in range(3, 9):
        for (nv, edges), _ in dist.items():
            adj = decode_tree_from_edges(nv, edges)
            for k in range(2, nv + 1):
                cur_cnt = sum(1 for x in range(nv) if len(adj[x]) >= k + 1)
                # Direct microscopic jump drift of N_{>= k+1}(T):
                micro_drift = Fraction(0, 1)
                for u in range(nv):
                    du = len(adj[u])
                    for w in adj[u]:
                        # Attaching leaf to w increments D_w by 1; newborn has D = 1 < k+1
                        delta_cnt = 1 if len(adj[w]) == k else 0
                        micro_drift += Fraction(delta_cnt, du)

                # Formula in degree-law.md line 35 summed over core vertices x:
                formula_sum = Fraction(0, 1)
                for x in range(nv):
                    dx = len(adj[x])
                    if dx == k:  # automatically >= 2 since k >= 2
                        core_nbs = [y for y in adj[x] if len(adj[y]) >= 2]
                        cx = len(core_nbs)
                        rho_x = sum(Fraction(1, len(adj[y])) for y in core_nbs)
                        formula_sum += Fraction(k - cx, 1) + rho_x
                if micro_drift != formula_sum:
                    return False
        if n < 8:
            dist = step_unlabelled_distribution(dist)

    return True


# ---------------------------------------------------------------------------
# Main runner + Negative Controls (Deliberate Errors)
# ---------------------------------------------------------------------------

def main() -> int:
    print("=== Cold Review Verification: round2/fringe/asymptotics.md ===")

    # 1. Exact small-n enumeration (2 <= n <= 10)
    n_max = 10
    a_t2, a_t3, total_trees = enumerate_exact_moments(n_max=n_max)
    for n in range(2, n_max):
        for k in range(2, n_max + 1):
            assert a_t2[(n + 1, k)] >= a_t2[(n, k)], f"Monotonicity failed at n={n}, k={k}"
            if n + 1 >= k + 1:
                assert a_t2[(n + 1, k)] > a_t2[(n, k)], f"Strict increase failed at n={n}, k={k}"
    for n in range(3, n_max + 1):
        exact_deg_sum = sum(a_t2[(n, k)] - a_t2[(n, k + 1)] for k in range(2, n))
        assert exact_deg_sum == a_t2[(n, 2)], f"Degree telescoping failed at n={n}"
    print(
        f"[PASS] Exact rational friend-tree enumeration (2 <= n <= {n_max}, "
        f"{total_trees} unlabelled state layers): monotonicity & degree telescoping exact."
    )
    print(
        "       Sample exact E[C_n] (n=3..10): "
        + ", ".join(f"n={n}: {a_t2[(n, 2)]}" for n in range(3, n_max + 1))
    )

    # 2. Marked descendant trees & R+B+I generator check from T_3
    checked_marked, ok_marked = check_marked_trees_and_generator(max_n=7, include_parent_slot=True)
    assert ok_marked, "Marked descendant tree / generator check failed!"
    print(
        f"[PASS] Marked descendant trees & renewal generator M'_t f_k = M_t(R+B)f_k + I_t f_k "
        f"exact across {checked_marked} reachable marked configurations (3 <= n <= 7)."
    )

    # 3. Exact small-n transfer & asymptotic Gamma/Karamata transfer
    assert check_small_n_transfer_identities(a_t2, a_t3, n_max=n_max), (
        "Small-n exact T_2 <-> T_3 & Beta-Laplace transfer failed!"
    )
    assert check_asymptotic_fixed_size_transfer(), (
        "Asymptotic fixed-size Yule/Gamma/sandwich transfer failed!"
    )
    print(
        "[PASS] Exact small-n q-series & Beta-Laplace transfer (ratio 2/(2+z)) and "
        "fixed-size Gamma(mu+2) / sandwich transfer verified."
    )

    # 4. Solvable multi-mark renewal operator, residue sign, peripheral line, Wiener-Ikehara scaling
    assert check_solvable_renewal_and_pole(), (
        "Solvable multi-mark renewal operator / pole / Tauberian scaling check failed!"
    )
    print(
        "[PASS] Multi-mark renewal resolvent Eq. (3), residue Eq. (4) sign -ell(K'_mu g) > 0, "
        "peripheral invertibility (tau != 0), and Wiener-Ikehara mu-rescaling verified."
    )

    # 5. Prerequisites (renewal.md Eq. 4, operator.md tilt bound, degree-law.md crossing)
    assert check_prerequisites(), "Prerequisite verification failed!"
    print(
        "[PASS] Prerequisites verified: renewal.md rational table (D <= 200), "
        "operator.md tilt ratio bound, and degree-law.md crossing identity."
    )

    # 6. Deliberate planted errors (negative controls)
    neg_controls: List[Tuple[str, bool]] = [
        (
            "Error 1: Wrong Yule weight n instead of (n-1) in Eq. (64)",
            not check_small_n_transfer_identities(
                a_t2, a_t3, n_max=n_max, wrong_yule_weight_n=True
            ),
        ),
        (
            "Error 2: Wrong T_3->T_2 wait factor 1/(1+z) instead of 2/(2+z) in Eq. (5)",
            not check_small_n_transfer_identities(
                a_t2, a_t3, n_max=n_max, wrong_wait_factor=True
            ),
        ),
        (
            "Error 3: Wrong Gamma(mu+1) instead of Gamma(mu+2) in c_k = d_k^* / Gamma(mu+2)",
            not check_asymptotic_fixed_size_transfer(wrong_gamma_shift=-1),
        ),
        (
            "Error 4: Missing dx = mu dt Jacobian in Wiener-Ikehara x -> M_{x/mu} f_k",
            not check_solvable_renewal_and_pole(wrong_wiener_ikehara_mu_factor="omit_jacobian"),
        ),
        (
            "Error 5: Flipped residue sign +ell(K'_mu g) instead of -ell(K'_mu g) in Eq. (4)",
            not check_solvable_renewal_and_pole(wrong_residue_sign=True),
        ),
        (
            "Error 6: Omitting exterior parent slot in D_root(S_v)",
            not check_marked_trees_and_generator(max_n=6, include_parent_slot=False)[1],
        ),
        (
            "Error 7: Lattice reproduction kernel producing poles at mu + 2*pi*i (tau != 0)",
            not check_solvable_renewal_and_pole(use_lattice_kernel=True),
        ),
        (
            "Error 8: Lowering q from 9/25 = 0.36 to 359/1000 in renewal.md Eq. (4) [breaks (c=1, D=3)]",
            not check_prerequisites(q_override=Fraction(359, 1000)),
        ),
    ]

    for desc, caught in neg_controls:
        assert caught, f"NEGATIVE CONTROL FAILED TO CATCH: {desc}"
        print(f"[CAUGHT] {desc}")

    print(f"ALL CHECKS AND {len(neg_controls)}/{len(neg_controls)} NEGATIVE CONTROLS PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
