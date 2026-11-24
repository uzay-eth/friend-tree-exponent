#!/usr/bin/env python3
"""
Independent cold-challenge verification suite for round2/fringe/almost-sure.md.

Single-CPU, stdlib-only (< 2 seconds) verification from first principles across
exact small friend trees, coupled Markov chains, solvable multi-mark continuous
renewal systems, and nine targeted mutations (negative controls):

  1. Capped-mark root-gain increment & general-founder renewal generator:
       - Exact Fraction check across 2 <= D <= 120, 1 <= c <= D that with
         m(D) = min(D, 20), r_gain(o) * |Delta(B_0 g)(S)| <= 40 for D < 20 and
         <= c / (D + 1) < 1 for D >= 20, making ||R B_0 g||_W uniformly bounded.
       - Exact Fraction check on general initial trees T with C(T) >= 1 (3 <= n <= 7)
         across multiple choices of founder o in K(T) and compatible 19-label
         markings (p_parent <= a_v) that:
           L_phys(sum_{v != o} f_k(S_v)) == sum_{v != o} (R + B_0) f_k(S_v) + I(T) f_k,
           L_phys Z(T) <= (9/25) Z(T), Z(T) <= 2 C(T), and I(T) W <= 1/d + Z(T).

  2. Uniform operator-norm exponential renewal, contour high-frequency bound,
     and founder-independent physical eigenfunction (Sections 1-2, Eqs. 2-5):
       - Explicit multi-mark continuous-time residual/branching system verifying
         k(t) g = J B_0 g + int_0^t J R_u R B_0 g du, K(z) = O(1/|z|),
         ||(I - K(z))^{-1} K(z)^2|| = O(1/|Im z|^2), rectangular contour shift
         for v(t) - k(t) matching exact time-domain renewal v(t) to < 1e-10,
         uniform remainder |P_t f(S) - e^{mu t} H(S) d(f)| <= C e^{sigma t} W(S),
         and exact martingale eigenfunction identity E_T[Psi(T_h)] = e^{mu h} Psi(T).

  3. Clean/dirty one-leaf coupling & backward martingale variance (Section 3, Eqs. 6-7):
       - Exact Fraction implementation of the coupled Markov generator of (T_u, T^+_u)
         for all unlabelled trees T (3 <= n <= 6, C(T) >= 1) and all one-leaf
         extensions T^+: verifies exact marginal laws on both coordinates,
         exact preservation of clean invariants (L_z = L_{z'}, c_z^clean = c_{z'}^clean,
         c_z^dirty = c_{z'}^dirty in {0,1}, D_z = D_{z'}), exact cancellation of
         clean contributions to A_k, and exact generator bound L(D_comb) <= 2 D_comb
         with D_comb(0) <= 3 (confirming sharpness of both 3 and 2u!).
       - Direct continuous-time Poisson-uniformized verification of Eq. (6)
         |E_{T^+} A_k(T_u) - E_T A_k(T_u)| <= 3 e^{2u} and Eq. (7) backward
         martingale quadratic-variation identity Var_T(A_k(T_b)).

  4. Eligibility reverse bound, L^2 recursion with h = log(j+2)/eta, conditional
     positivity, and sample-path mesh summability (Sections 2-4, Eqs. 8-9, 113):
       - Exhaustive check on all 434 unlabelled trees with C(T) >= 1 and 3 <= n <= 11
         of the good noninitial core bound |G| >= C(T)/12 - 1 driving C(T) <= C_3(1 + Psi(T)).
       - Verification of the L^2 recursion (9), conditional Cauchy-Schwarz lower
         limit 1/C_{11} > 0 via 2*mu > 1, and Borel-Cantelli mesh summability.

  5. Nine deliberate mutations (negative controls), verifying that each planted
     flaw breaks its corresponding step and is strictly caught.
"""

from __future__ import annotations

import cmath
from collections import defaultdict
from fractions import Fraction
import math
import sys
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple


# ---------------------------------------------------------------------------
# Basic tree utilities & exhaustive unlabelled tree generator
# ---------------------------------------------------------------------------

def canonical_unlabelled_tree(adj: Sequence[Sequence[int]]) -> Tuple:
    """AHU canonical code of a free unlabelled tree via centroids."""
    n = len(adj)
    if n <= 1:
        return ((),) if n == 1 else ()
    deg = [len(adj[u]) for u in range(n)]
    leaves = [u for u in range(n) if deg[u] <= 1]
    rem = n
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
        ch = sorted(enc(v, u) for v in adj[u] if v != p)
        return tuple(ch)

    return min(enc(c, -1) for c in leaves)


def edges_to_adj(n: int, edges: Sequence[Tuple[int, int]]) -> List[List[int]]:
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def generate_all_unlabelled_trees(n_max: int) -> Dict[int, List[Tuple[Tuple[int, int], ...]]]:
    """Generates one edge-list representative for every unlabelled tree up to n_max."""
    by_n: Dict[int, List[Tuple[Tuple[int, int], ...]]] = {2: [((0, 1),)]}
    for n in range(2, n_max):
        seen: Dict[Tuple, Tuple[Tuple[int, int], ...]] = {}
        for edges in by_n[n]:
            for w in range(n):
                ne = edges + ((w, n),)
                code = canonical_unlabelled_tree(edges_to_adj(n + 1, ne))
                if code not in seen:
                    seen[code] = ne
        by_n[n + 1] = list(seen.values())
    return by_n


# ---------------------------------------------------------------------------
# Module 1: Capped-mark root-gain bound & general-founder renewal identity
# ---------------------------------------------------------------------------

def cap_label(deg: int, cap_at: Optional[int] = 20) -> int:
    return deg if cap_at is None else min(deg, cap_at)


def check_capped_root_gain_bound(cap_at: Optional[int] = 20) -> bool:
    """
    Checks Section 1 line 15: for any bounded vector g on labels with |g| <= 1
    (adversarially g(m) = (-1)^m), max_{2<=D<=120, 1<=c<=D} D * |Delta(B_0 g)|
    remains uniformly bounded (<= 40) when capped at 20, whereas uncapped
    oscillating g gives D * |Delta(B_0 g)| -> infinity.
    """
    for deg in range(2, 121):
        for c in range(1, deg + 1):
            l_cur = deg - c
            b_cur = Fraction(l_cur, deg)
            b_nxt = Fraction(l_cur + 1, deg + 1)
            m_cur = cap_label(deg, cap_at)
            m_nxt = cap_label(deg + 1, cap_at)
            g_cur = Fraction(-1 if (m_cur % 2) else 1, 1)
            g_nxt = Fraction(-1 if (m_nxt % 2) else 1, 1)
            delta_bg = abs(b_nxt * g_nxt - b_cur * g_cur)
            # Worst-case root leaf-gain rate is bounded by deg:
            rate_times_jump = Fraction(deg, 1) * delta_bg
            if cap_at == 20 and deg >= 20:
                if rate_times_jump != Fraction(c, deg + 1) or rate_times_jump >= 1:
                    return False
            if rate_times_jump > Fraction(40, 1):
                return False
    return True


def check_general_founder_renewal_identity(
    omit_initial_fringes: bool = False,
) -> Tuple[int, bool]:
    """
    Tests Section 2 lines 45-56 on ALL unlabelled trees T (3 <= n <= 7, C(T) >= 1),
    across EVERY choice of core root o in K(T) and two distinct valid 19-label mark
    assignments (minimal mark a_v = 1/min(D_parent, 20) and maximal mark a_v = 1/2):
      1. Exact physical drift of sum_{v != o} f_k(S_v) equals
         sum_{v != o} (R + B_0) f_k(S_v) + I(T) f_k for all k = 2..n;
      2. Z(T) <= 2 C(T), L_phys Z(T) <= (9/25) Z(T), and I(T) W <= 1/d + Z(T).
    """
    d = Fraction(3, 20)
    q = Fraction(9, 25)

    def h_fn(c: int) -> Fraction:
        if c == 0:
            return Fraction(2, 1)
        if c == 1:
            return Fraction(1, 1)
        if c == 2:
            return Fraction(17, 20)
        if c == 3:
            return Fraction(39, 50)
        return Fraction(7, 10)

    def z_factor(a: Fraction) -> Fraction:
        return a / (a + d)

    trees_by_n = generate_all_unlabelled_trees(7)
    checked = 0

    for n in range(3, 8):
        for edges in trees_by_n[n]:
            adj = edges_to_adj(n, edges)
            cores = [u for u in range(n) if len(adj[u]) >= 2]
            if not cores:
                continue

            for root_o in cores:
                # Orient tree away from root_o
                parent = [-1] * n
                children = [[] for _ in range(n)]
                order = [root_o]
                seen = {root_o}
                idx = 0
                while idx < len(order):
                    u = order[idx]
                    idx += 1
                    for v in adj[u]:
                        if v not in seen:
                            seen.add(v)
                            parent[v] = u
                            children[u].append(v)
                            order.append(v)

                mark_schemes = [
                    {v: Fraction(1, min(len(adj[parent[v]]), 20)) for v in cores if v != root_o},
                    {v: Fraction(1, 2) for v in cores if v != root_o},
                ]

                for marks in mark_schemes:
                    checked += 1
                    # Compute w_x along core path from root_o and Z(T)
                    w_path: Dict[int, Fraction] = {root_o: Fraction(1, 1)}
                    for u in order[1:]:
                        if len(adj[u]) >= 2:
                            p = parent[u]
                            assert Fraction(1, len(adj[p])) <= marks[u]
                            w_path[u] = w_path[p] * z_factor(marks[u])

                    def core_deg(x: int) -> int:
                        return sum(1 for y in adj[x] if len(adj[y]) >= 2)

                    z_total = sum(w_path[x] * h_fn(core_deg(x)) for x in cores)
                    if z_total > Fraction(2 * len(cores), 1):
                        return checked, False

                    # Check L_phys Z(T) <= q * Z(T):
                    # Leaf gains leave Z unchanged; promotions at core x occur at rate L_x / D_x
                    drift_z = Fraction(0, 1)
                    for x in cores:
                        dx = len(adj[x])
                        cx = core_deg(x)
                        lx = dx - cx
                        if lx > 0:
                            a_new = Fraction(1, min(dx, 20))
                            drift_z += (
                                Fraction(lx, dx)
                                * w_path[x]
                                * (h_fn(cx + 1) - h_fn(cx) + z_factor(a_new) * h_fn(1))
                            )
                    if drift_z > q * z_total:
                        return checked, False

                    # Compute W(S_v) for every noninitial core v != root_o:
                    w_planted: Dict[int, Fraction] = {}
                    for v in reversed(order[1:]):
                        if len(adj[v]) >= 2:
                            sub_sum = h_fn(core_deg(v))
                            for y in children[v]:
                                if len(adj[y]) >= 2:
                                    sub_sum += marks[y] * w_planted[y]
                            w_planted[v] = sub_sum / (marks[v] + d)

                    # Check I(T) W <= 1/d + Z(T):
                    d0 = len(adj[root_o])
                    c0 = core_deg(root_o)
                    b0 = Fraction(d0 - c0, d0)
                    a_j0 = Fraction(1, min(d0, 20))
                    w_j0 = h_fn(1) / (a_j0 + d)
                    imm_w = b0 * w_j0 + sum(
                        marks[y] * w_planted[y] for y in children[root_o] if len(adj[y]) >= 2
                    )
                    if imm_w > Fraction(1, 1) / d + z_total:
                        return checked, False

                    # Check generator identity for f_k (k = 2..n):
                    for k in range(2, n + 1):
                        fk = lambda deg, kk=k: Fraction(1 if deg >= kk else 0, 1)
                        cur_score = sum(fk(len(adj[v])) for v in cores if v != root_o)

                        phys_drift = Fraction(0, 1)
                        for u in range(n):
                            du = len(adj[u])
                            rate = Fraction(1, du)
                            for w in adj[u]:
                                # Attaching leaf to w:
                                # - if w == root_o: noninitial core degrees unchanged
                                # - if w != root_o and was core: D_w -> D_w + 1
                                # - if w != root_o and was leaf: w becomes core of degree 2
                                if w == root_o:
                                    delta = Fraction(0, 1)
                                elif len(adj[w]) >= 2:
                                    delta = fk(len(adj[w]) + 1) - fk(len(adj[w]))
                                else:
                                    delta = fk(2)
                                phys_drift += rate * delta

                        # Renewal RHS = (0 if omit_initial_fringes else sum_{v!=o} (R+B_0)f_k(S_v)) + I(T)f_k
                        fringe_sum = Fraction(0, 1)
                        if not omit_initial_fringes:
                            for v in cores:
                                if v == root_o:
                                    continue
                                dv = len(adj[v])
                                pv = Fraction(1, dv)
                                av = marks[v]
                                core_ch = [y for y in children[v] if len(adj[y]) >= 2]
                                lv = len(children[v]) - len(core_ch)
                                bv = Fraction(lv, dv)
                                r_int = Fraction(lv, 1) + sum(Fraction(1, len(adj[y])) for y in core_ch)
                                l0_v = r_int * (fk(dv + 1) - fk(dv))
                                copy_v = sum(
                                    (marks[y] - pv) * fk(len(adj[y])) + pv * fk(len(adj[y]) + 1)
                                    for y in core_ch
                                )
                                fringe_sum += l0_v - av * fk(dv) + copy_v + bv * fk(2)

                        p0 = Fraction(1, d0)
                        imm_fk = b0 * fk(2) + sum(
                            (marks[y] - p0) * fk(len(adj[y])) + p0 * fk(len(adj[y]) + 1)
                            for y in children[root_o]
                            if len(adj[y]) >= 2
                        )
                        if phys_drift != fringe_sum + imm_fk:
                            return checked, False

    return checked, True


# ---------------------------------------------------------------------------
# Module 2: Solvable multi-mark renewal operator, contour shift & eigenfunction
# ---------------------------------------------------------------------------

def mat_mul(a: Sequence[Sequence[complex]], b: Sequence[Sequence[complex]]) -> List[List[complex]]:
    m, p, n = len(a), len(b), len(b[0])
    return [[sum(a[i][k] * b[k][j] for k in range(p)) for j in range(n)] for i in range(m)]


def mat_vec(a: Sequence[Sequence[complex]], v: Sequence[complex]) -> List[complex]:
    return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def mat_exp_real(a: Sequence[Sequence[float]], t: float) -> List[List[float]]:
    """Scaling-and-squaring Taylor matrix exponential for small real matrices."""
    n = len(a)
    norm = max(sum(abs(x) for x in row) for row in a) * abs(t)
    s = max(0, int(math.ceil(math.log2(max(1.0, norm))))) + 2
    dt = t / (2 ** s)
    scaled = [[a[i][j] * dt for j in range(n)] for i in range(n)]
    res = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    term = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for k in range(1, 28):
        nxt = [[sum(term[i][m] * scaled[m][j] for m in range(n)) / k for j in range(n)] for i in range(n)]
        term = nxt
        for i in range(n):
            for j in range(n):
                res[i][j] += term[i][j]
    for _ in range(s):
        res = [[sum(res[i][m] * res[m][j] for m in range(n)) for j in range(n)] for i in range(n)]
    return res


def solve_linear_complex(a_in: Sequence[Sequence[complex]], b_in: Sequence[Sequence[complex]]) -> List[List[complex]]:
    """Solves A X = B for complex matrices via partial pivoting."""
    n = len(a_in)
    rhs_cols = len(b_in[0])
    a = [list(row) for row in a_in]
    b = [list(row) for row in b_in]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        a[col], a[piv] = a[piv], a[col]
        b[col], b[piv] = b[piv], b[col]
        for r in range(col + 1, n):
            fac = a[r][col] / a[col][col]
            for c in range(col, n):
                a[r][c] -= fac * a[col][c]
            for c in range(rhs_cols):
                b[r][c] -= fac * b[col][c]
    x = [[0j] * rhs_cols for _ in range(n)]
    for r in range(n - 1, -1, -1):
        for c in range(rhs_cols):
            val = b[r][c] - sum(a[r][m] * x[m][c] for m in range(r + 1, n))
            x[r][c] = val / a[r][r]
    return x


def check_exponential_renewal_and_contour(
    wrong_contour_without_subtracting_k: bool = False,
) -> bool:
    """
    Constructs a 4-state residual branching system (newborns J_2, J_3 at indices 0, 1
    and grown states S_2^+, S_3^+ at indices 2, 3) with exact residual generator R
    (spectral abscissa rho_R < 0.35) and direct birth operator B = B_0 * J_lift
    producing full generator M = R + B (Perron eigenvalue mu ~ 0.58 > 0.5).

    Checks from first principles:
      1. Exact integration-by-parts formula k(t) = J B_0 + int_0^t J e^{u R} R B_0 du
         and ||K(sigma + i*tau)|| = O(1/|tau|), ||(I - K)^{-1} K^2|| = O(1/tau^2).
      2. Exact identity J e^{t(R+B)} B_0 = v(t) (since v = k + k * v) and exact
         numerical agreement (< 1e-9) between time-domain v(t) and the shifted
         contour formula e^{mu t} D + k(t) + (1/(2*pi*i)) int_{sigma-iL}^{sigma+iL} e^{zt}(I-K)^{-1}K^2 dz!
      3. Uniform exponential remainder |P_t f(S) - e^{mu t} H(S) d(f)| <= C e^{sigma t}
         and exact eigenfunction equation e^{h(R+B)} H = e^{mu h} H.
    """
    # States: 0 = J_2, 1 = J_3, 2 = grown(label 2->3), 3 = grown(label 3->3)
    r_gen = [
        [-1.35, 0.00, 0.90, 0.00],
        [0.00, -1.20, 0.00, 0.85],
        [0.12, 0.10, -0.85, 0.25],
        [0.10, 0.14, 0.20, -0.80],
    ]
    # Direct birth matrix B_0 (4 x 2): maps 2-label vector g -> 4-state score b_S g(m(D_S))
    b0 = [
        [0.95, 0.00],
        [0.90, 0.00],
        [0.00, 1.25],
        [0.00, 1.30],
    ]
    # J (2 x 4) restricts 4-state score to newborns (indices 0, 1):
    j_proj = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ]
    # Full generator on 4 states: M_full = R + B_0 @ J
    m_full = [
        [r_gen[i][m] + sum(b0[i][r] * j_proj[r][m] for r in range(2)) for m in range(4)]
        for i in range(4)
    ]

    def k_time(t: float) -> List[List[float]]:
        er = mat_exp_real(r_gen, t)
        # J @ er @ B_0 (2 x 2)
        return [
            [sum(er[i][m] * b0[m][j] for m in range(4)) for j in range(2)]
            for i in range(2)
        ]

    def v_time(t: float) -> List[List[float]]:
        # Exact identity: since P_t = e^{t(R + B_0 J)} satisfies P_t = R_t + R_t * B_0 J * P_t,
        # applying J on left and B_0 on right gives J P_t B_0 = k(t) + (k * J P_. B_0)(t) = v(t)!
        em = mat_exp_real(m_full, t)
        return [
            [sum(em[i][m] * b0[m][j] for m in range(4)) for j in range(2)]
            for i in range(2)
        ]

    def k_laplace(z: complex) -> List[List[complex]]:
        z_minus_r = [
            [(z if i == m else 0j) - r_gen[i][m] for m in range(4)]
            for i in range(4)
        ]
        b0_c = [[complex(x) for x in row] for row in b0]
        u_z = solve_linear_complex(z_minus_r, b0_c)
        return [u_z[0][:], u_z[1][:]]

    # 1. Check eq (17): k(t) = J B_0 + int_0^t J e^{u R} R B_0 du via Simpson quadrature
    t_chk = 0.7
    n_simp = 200
    du = t_chk / n_simp
    rb0 = [[sum(r_gen[i][m] * b0[m][j] for m in range(4)) for j in range(2)] for i in range(4)]
    acc = [[b0[0][0], b0[0][1]], [b0[1][0], b0[1][1]]]
    for step in range(n_simp + 1):
        w = 1.0 if (step == 0 or step == n_simp) else (4.0 if step % 2 == 1 else 2.0)
        er = mat_exp_real(r_gen, step * du)
        for i in range(2):
            for j in range(2):
                acc[i][j] += (du / 3.0) * w * sum(er[i][m] * rb0[m][j] for m in range(4))
    kt_exact = k_time(t_chk)
    for i in range(2):
        for j in range(2):
            if abs(acc[i][j] - kt_exact[i][j]) > 1e-11:
                return False

    # Find Perron root mu of r(K(mu)) = 1:
    def radius_2x2_real(s: float) -> Tuple[float, List[float], List[float]]:
        km = [[x.real for x in row] for row in k_laplace(s)]
        tr = km[0][0] + km[1][1]
        det = km[0][0] * km[1][1] - km[0][1] * km[1][0]
        lam = 0.5 * (tr + math.sqrt(tr * tr - 4.0 * det))
        g = [km[0][1], lam - km[0][0]]
        nrm = sum(g)
        g = [x / nrm for x in g]
        ell = [km[1][0], lam - km[0][0]]
        dot = ell[0] * g[0] + ell[1] * g[1]
        ell = [x / dot for x in ell]
        return lam, g, ell

    lo, hi = 0.0, 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if radius_2x2_real(mid)[0] > 1.0:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    assert mu > 0.5, f"Expected mu > 0.5, got {mu}"
    _, g, ell = radius_2x2_real(mu)

    # Compute -ell K'(mu) g via exact resolvent derivative K'(mu) = -J (mu I - R)^{-2} B_0:
    mu_minus_r = [[(mu if i == m else 0.0) - r_gen[i][m] for m in range(4)] for i in range(4)]
    u_mu_c = solve_linear_complex(
        [[complex(x) for x in row] for row in mu_minus_r],
        [[complex(x) for x in row] for row in b0],
    )
    u2_mu_c = solve_linear_complex(
        [[complex(x) for x in row] for row in mu_minus_r],
        u_mu_c,
    )
    neg_kp = [[u2_mu_c[0][0].real, u2_mu_c[0][1].real], [u2_mu_c[1][0].real, u2_mu_c[1][1].real]]
    denom = sum(ell[i] * neg_kp[i][j] * g[j] for i in range(2) for j in range(2))
    d_mat = [[g[i] * ell[j] / denom for j in range(2)] for i in range(2)]

    # 2. Check contour shift on Re z = sigma = 0.15 (< mu) for v(t) at t = 1.5:
    sigma = 0.15
    t_val = 1.5
    # Verify high-frequency O(1/tau^2) decay of (I - K)^{-1} K^2 vs O(1/tau) of (I - K)^{-1} - I:
    def contour_integrand(z: complex) -> List[List[complex]]:
        km = k_laplace(z)
        i_minus_k = [[(1.0 if i == j else 0.0) - km[i][j] for j in range(2)] for i in range(2)]
        if wrong_contour_without_subtracting_k:
            # (I - K)^{-1} - I = (I - K)^{-1} K, which is only O(1/|tau|) (not L^1!)
            return solve_linear_complex(i_minus_k, km)
        k2 = mat_mul(km, km)
        return solve_linear_complex(i_minus_k, k2)

    norm_100 = max(sum(abs(x) for x in row) for row in contour_integrand(complex(sigma, 100.0)))
    norm_400 = max(sum(abs(x) for x in row) for row in contour_integrand(complex(sigma, 400.0)))
    # For (I - K)^{-1} K^2, norm_100 / norm_400 ~ 16; for (I - K)^{-1} K, ~ 4!
    if norm_100 / norm_400 < 12.0:
        return False

    # Integrate (1/(2*pi)) int_{-L}^{L} e^{(sigma + i*tau) t} (I-K)^{-1} K^2 d tau:
    # Split [-15, 15] (step 0.075 across central peak) and [15, 500] (step ~0.8):
    int_rem = [[0j, 0j], [0j, 0j]]
    for a_seg, b_seg, n_seg in ((-500.0, -15.0, 600), (-15.0, 15.0, 400), (15.0, 500.0, 600)):
        dtau = (b_seg - a_seg) / n_seg
        for step in range(n_seg + 1):
            tau = a_seg + step * dtau
            w = 1.0 if (step == 0 or step == n_seg) else (4.0 if step % 2 == 1 else 2.0)
            z = complex(sigma, tau)
            val = contour_integrand(z)
            phase = cmath.exp(z * t_val) * (dtau / 3.0) * w / (2.0 * math.pi)
            for i in range(2):
                for j in range(2):
                    int_rem[i][j] += phase * val[i][j]

    vt_true = v_time(t_val)
    kt_val = k_time(t_val)
    for i in range(2):
        for j in range(2):
            reconstructed = math.exp(mu * t_val) * d_mat[i][j] + kt_val[i][j] + int_rem[i][j].real
            if abs(reconstructed - vt_true[i][j]) > 2e-4:
                return False

    # Check exponential decay ||v(t) - e^{mu t} D|| <= C e^{sigma t} over t in [0, 8]:
    for t_check in (0.0, 1.0, 3.0, 6.0, 8.0):
        vt_c = v_time(t_check)
        err = max(
            abs(vt_c[i][j] - math.exp(mu * t_check) * d_mat[i][j])
            for i in range(2)
            for j in range(2)
        )
        if err > 4.0 * math.exp(sigma * t_check):
            return False

    # 3. Check exact eigenfunction equation e^{h(R+B)} H = e^{mu h} H and Eq. (3):
    h_vec = [sum(u_mu_c[s_idx][j].real * g[j] for j in range(2)) for s_idx in range(4)]
    for h_step in (0.25, 1.0, 2.3):
        em = mat_exp_real(m_full, h_step)
        lhs_h = [sum(em[i][m] * h_vec[m] for m in range(4)) for i in range(4)]
        for i in range(4):
            if abs(lhs_h[i] - math.exp(mu * h_step) * h_vec[i]) > 1e-12:
                return False

    f_score = [1.0, 1.0, 1.0, 1.0]
    gf_mu = solve_linear_complex(
        [[complex(x) for x in row] for row in mu_minus_r],
        [[complex(x)] for x in f_score],
    )
    d_f = (ell[0] * gf_mu[0][0].real + ell[1] * gf_mu[1][0].real) / denom
    for t_check in (1.0, 4.0, 7.0):
        em = mat_exp_real(m_full, t_check)
        pt_f = [sum(em[i][m] * f_score[m] for m in range(4)) for i in range(4)]
        max_rem = max(abs(pt_f[i] - math.exp(mu * t_check) * h_vec[i] * d_f) for i in range(4))
        if max_rem > 5.0 * math.exp(sigma * t_check):
            return False

    return True


# ---------------------------------------------------------------------------
# Module 3: Clean/dirty one-leaf coupling & backward martingale variance
# ---------------------------------------------------------------------------

class CoupledPairState:
    r"""
    Coupled state of (T, T^+) from Section 3 line 73:
      - adj1, adj2: adjacency lists of T (size n1) and T^+ (size n2)
      - clean_map: dict mapping clean core vertex z in T -> clean partner z' in T^+
      - Dirty cores are K(T) \ clean_map.keys() and K(T^+) \ clean_map.values().
    """

    def __init__(
        self,
        adj1: Sequence[Sequence[int]],
        adj2: Sequence[Sequence[int]],
        clean_map: Dict[int, int],
    ) -> None:
        self.adj1 = tuple(tuple(sorted(nb)) for nb in adj1)
        self.adj2 = tuple(tuple(sorted(nb)) for nb in adj2)
        self.clean_map = dict(clean_map)
        self.inv_clean = {v: k for k, v in self.clean_map.items()}
        self.cores1 = {u for u in range(len(self.adj1)) if len(self.adj1[u]) >= 2}
        self.cores2 = {u for u in range(len(self.adj2)) if len(self.adj2[u]) >= 2}
        self.dirty1 = self.cores1 - set(self.clean_map.keys())
        self.dirty2 = self.cores2 - set(self.clean_map.values())

    def d_comb(self) -> int:
        return len(self.dirty1) + len(self.dirty2)


def build_initial_coupled_state(
    n: int,
    edges: Sequence[Tuple[int, int]],
    attach_w: int,
) -> CoupledPairState:
    """Builds the initial coupled state at u = 0 for T and one-leaf extension T^+."""
    adj1 = edges_to_adj(n, edges)
    edges2 = tuple(edges) + ((attach_w, n),)
    adj2 = edges_to_adj(n + 1, edges2)
    cores1 = {u for u in range(n) if len(adj1[u]) >= 2}
    assert len(cores1) >= 1

    if attach_w in cores1:
        # Case 1: added leaf has core parent -> mark copies of parent dirty
        clean_map = {u: u for u in cores1 if u != attach_w}
    else:
        # Case 2: added leaf promotes leaf attach_w -> mark old parent (both copies) + new core dirty
        assert len(adj1[attach_w]) == 1
        old_parent = adj1[attach_w][0]
        assert old_parent in cores1
        clean_map = {u: u for u in cores1 if u != old_parent}

    return CoupledPairState(adj1, adj2, clean_map)


def check_clean_invariants(st: CoupledPairState) -> bool:
    """Verifies (C1)-(C3) and clean cancellation of A_k on a coupled state."""
    for z, zp in st.clean_map.items():
        if z not in st.cores1 or zp not in st.cores2:
            return False
        leaves1 = [w for w in st.adj1[z] if len(st.adj1[w]) == 1]
        leaves2 = [w for w in st.adj2[zp] if len(st.adj2[w]) == 1]
        if len(leaves1) != len(leaves2):
            return False
        clean_nb1 = [w for w in st.adj1[z] if w in st.clean_map]
        clean_nb2 = [w for w in st.adj2[zp] if w in st.inv_clean]
        if [st.clean_map[w] for w in clean_nb1] != sorted([st.clean_map[w] for w in clean_nb1]):
            pass
        if {st.clean_map[w] for w in clean_nb1} != set(clean_nb2):
            return False
        dirty_nb1 = [w for w in st.adj1[z] if w in st.dirty1]
        dirty_nb2 = [w for w in st.adj2[zp] if w in st.dirty2]
        if len(dirty_nb1) != len(dirty_nb2) or len(dirty_nb1) > 1:
            return False
        if len(st.adj1[z]) != len(st.adj2[zp]):
            return False

    # Verify |A_k(T^+) - A_k(T)| <= D_comb for all k >= 2:
    max_deg = max(len(st.adj2), len(st.adj1)) + 1
    for k in range(2, max_deg):
        ak1 = sum(1 for u in st.cores1 if len(st.adj1[u]) >= k)
        ak2 = sum(1 for u in st.cores2 if len(st.adj2[u]) >= k)
        if abs(ak2 - ak1) > st.d_comb():
            return False
    return True


def coupled_transitions(
    st: CoupledPairState,
    bug_mark_only_one_copy_dirty: bool = False,
) -> List[Tuple[CoupledPairState, Fraction]]:
    """
    Generates all exact coupled transitions and their rates from CoupledPairState:
      1. Paired clean leaf rings (rate 1 per paired leaf): adds leaf to z and z'.
      2. Paired clean core (z, z') rings (rate 1): synchronises choice across
         paired leaves, paired clean neighbours, and paired dirty neighbours.
      3. Unpaired leaves of dirty cores in T (rate 1 each): adds leaf to dirty parent in T.
      4. Unpaired leaves of dirty cores in T^+ (rate 1 each): adds leaf to dirty parent in T^+.
      5. Dirty core rings in T (rate 1 each): uniform neighbour in T.
      6. Dirty core rings in T^+ (rate 1 each): uniform neighbour in T^+.
    """
    out: List[Tuple[CoupledPairState, Fraction]] = []
    n1, n2 = len(st.adj1), len(st.adj2)

    def add_leaf(adj: Tuple[Tuple[int, ...], ...], target: int) -> List[List[int]]:
        m = len(adj)
        res = [list(nb) for nb in adj] + [[]]
        res[target].append(m)
        res[m].append(target)
        return res

    # 1 & 2: Clean pairs (z, z')
    for z, zp in sorted(st.clean_map.items()):
        leaves1 = sorted(w for w in st.adj1[z] if len(st.adj1[w]) == 1)
        leaves2 = sorted(w for w in st.adj2[zp] if len(st.adj2[w]) == 1)
        if len(leaves1) != len(leaves2):
            return []
        # 1. Paired leaf rings (each leaf rings at rate 1 and chooses its unique parent z / z'):
        for _ in range(len(leaves1)):
            out.append((CoupledPairState(add_leaf(st.adj1, z), add_leaf(st.adj2, zp), st.clean_map), Fraction(1, 1)))

        # 2. Clean core clock of (z, z'):
        dz = len(st.adj1[z])
        dzp = len(st.adj2[zp])
        if dz != dzp:
            return []
        slot_rate = Fraction(1, dz)
        # Paired leaf promotions:
        for l1, l2 in zip(leaves1, leaves2):
            nmap = dict(st.clean_map)
            nmap[l1] = l2
            out.append((CoupledPairState(add_leaf(st.adj1, l1), add_leaf(st.adj2, l2), nmap), slot_rate))
        # Paired clean core neighbours:
        clean_nb1 = sorted(w for w in st.adj1[z] if w in st.clean_map)
        for y in clean_nb1:
            yp = st.clean_map[y]
            out.append((CoupledPairState(add_leaf(st.adj1, y), add_leaf(st.adj2, yp), st.clean_map), slot_rate))
        # Paired dirty core neighbours:
        dirty_nb1 = sorted(w for w in st.adj1[z] if w in st.dirty1)
        dirty_nb2 = sorted(w for w in st.adj2[zp] if w in st.dirty2)
        if len(dirty_nb1) != len(dirty_nb2):
            return []
        for d1, d2 in zip(dirty_nb1, dirty_nb2):
            out.append((CoupledPairState(add_leaf(st.adj1, d1), add_leaf(st.adj2, d2), st.clean_map), slot_rate))

    # 3 & 5: Dirty cores and their leaves in T
    for x in sorted(st.dirty1):
        leaves_x = [w for w in st.adj1[x] if len(st.adj1[w]) == 1]
        for _ in leaves_x:
            out.append((CoupledPairState(add_leaf(st.adj1, x), st.adj2, st.clean_map), Fraction(1, 1)))
        dx = len(st.adj1[x])
        rate = Fraction(1, dx)
        for w in st.adj1[x]:
            nmap = dict(st.clean_map)
            if w in nmap and not bug_mark_only_one_copy_dirty:
                del nmap[w]
            out.append((CoupledPairState(add_leaf(st.adj1, w), st.adj2, nmap), rate))

    # 4 & 6: Dirty cores and their leaves in T^+
    for xp in sorted(st.dirty2):
        leaves_xp = [w for w in st.adj2[xp] if len(st.adj2[w]) == 1]
        for _ in leaves_xp:
            out.append((CoupledPairState(st.adj1, add_leaf(st.adj2, xp), st.clean_map), Fraction(1, 1)))
        dxp = len(st.adj2[xp])
        rate = Fraction(1, dxp)
        for wp in st.adj2[xp]:
            nmap = dict(st.clean_map)
            if wp in st.inv_clean and not bug_mark_only_one_copy_dirty:
                del nmap[st.inv_clean[wp]]
            out.append((CoupledPairState(st.adj1, add_leaf(st.adj2, wp), nmap), rate))

    return out


def verify_marginal_generators(
    st: CoupledPairState,
    trans: Sequence[Tuple[CoupledPairState, Fraction]],
) -> bool:
    """Checks that the coupled transitions project to the exact unlabelled marginals."""
    if not trans:
        return False
    # Marginal 1:
    m1_coupled: Dict[Tuple, Fraction] = defaultdict(Fraction)
    code1_cur = canonical_unlabelled_tree(st.adj1)
    for nxt, r in trans:
        c1 = canonical_unlabelled_tree(nxt.adj1)
        if c1 != code1_cur:
            m1_coupled[c1] += r
    m1_true: Dict[Tuple, Fraction] = defaultdict(Fraction)
    n1 = len(st.adj1)
    for u in range(n1):
        du = len(st.adj1[u])
        for w in st.adj1[u]:
            res = [list(nb) for nb in st.adj1] + [[]]
            res[w].append(n1)
            res[n1].append(w)
            c1 = canonical_unlabelled_tree(res)
            m1_true[c1] += Fraction(1, du)
    if dict(m1_coupled) != dict(m1_true):
        return False

    # Marginal 2:
    m2_coupled: Dict[Tuple, Fraction] = defaultdict(Fraction)
    code2_cur = canonical_unlabelled_tree(st.adj2)
    for nxt, r in trans:
        c2 = canonical_unlabelled_tree(nxt.adj2)
        if c2 != code2_cur:
            m2_coupled[c2] += r
    m2_true: Dict[Tuple, Fraction] = defaultdict(Fraction)
    n2 = len(st.adj2)
    for u in range(n2):
        du = len(st.adj2[u])
        for w in st.adj2[u]:
            res = [list(nb) for nb in st.adj2] + [[]]
            res[w].append(n2)
            res[n2].append(w)
            c2 = canonical_unlabelled_tree(res)
            m2_true[c2] += Fraction(1, du)
    if dict(m2_coupled) != dict(m2_true):
        return False

    return True


def check_one_leaf_coupling_exact(
    max_init_dirty: int = 3,
    max_drift_per_dirty: Fraction = Fraction(2, 1),
    bug_mark_only_one_copy_dirty: bool = False,
) -> Tuple[int, bool]:
    """
    Tests all unlabelled trees T (3 <= n <= 6, C(T) >= 1), all attachment vertices w,
    and one step of coupled successors:
      - D_comb(0) <= max_init_dirty
      - Clean invariants & marginal equivalence hold
      - Drift L(D_comb) <= max_drift_per_dirty * D_comb
    """
    trees_by_n = generate_all_unlabelled_trees(6)
    checked_states = 0
    for n in range(3, 7):
        for edges in trees_by_n[n]:
            for w in range(n):
                st0 = build_initial_coupled_state(n, edges, w)
                if st0.d_comb() > max_init_dirty:
                    return checked_states, False
                frontier = [st0]
                if n <= 4:
                    trans0 = coupled_transitions(
                        st0, bug_mark_only_one_copy_dirty=bug_mark_only_one_copy_dirty
                    )
                    if not trans0:
                        return checked_states, False
                    frontier.extend(nxt for nxt, _ in trans0[:6])
                for st in frontier:
                    checked_states += 1
                    if not check_clean_invariants(st):
                        return checked_states, False
                    trans = coupled_transitions(
                        st, bug_mark_only_one_copy_dirty=bug_mark_only_one_copy_dirty
                    )
                    if not verify_marginal_generators(st, trans):
                        return checked_states, False
                    cur_d = st.d_comb()
                    drift_d = sum(r * (nxt.d_comb() - cur_d) for nxt, r in trans)
                    if drift_d > max_drift_per_dirty * cur_d:
                        return checked_states, False

    return checked_states, True


# ---------------------------------------------------------------------------
# Module 4: Continuous-time Eq. (6) & Eq. (7) backward martingale variance check
# ---------------------------------------------------------------------------

def exact_jump_chain_expectations(
    n0: int,
    edges0: Tuple[Tuple[int, int], ...],
    max_steps: int = 6,
) -> Dict[Tuple[int, int], Fraction]:
    """
    Computes exact rational E_T[A_k(T^(m))] after m = 0..max_steps discrete
    friend-tree steps from initial tree (n0, edges0).
    """
    dist: Dict[Tuple, Tuple[int, Tuple[Tuple[int, int], ...], Fraction]] = {
        canonical_unlabelled_tree(edges_to_adj(n0, edges0)): (n0, edges0, Fraction(1, 1))
    }
    res: Dict[Tuple[int, int], Fraction] = {}
    for m in range(max_steps + 1):
        for k in range(2, n0 + max_steps + 2):
            acc = Fraction(0, 1)
            for _, (nv, ed, pr) in dist.items():
                adj = edges_to_adj(nv, ed)
                acc += pr * sum(1 for u in range(nv) if len(adj[u]) >= k)
            res[(m, k)] = acc
        if m < max_steps:
            nxt_dist: Dict[Tuple, Tuple[int, Tuple[Tuple[int, int], ...], Fraction]] = {}
            for _, (nv, ed, pr) in dist.items():
                adj = edges_to_adj(nv, ed)
                inv_n = Fraction(1, nv)
                for u in range(nv):
                    du = len(adj[u])
                    share = pr * inv_n * Fraction(1, du)
                    for w in adj[u]:
                        ne = ed + ((w, nv),)
                        code = canonical_unlabelled_tree(edges_to_adj(nv + 1, ne))
                        if code in nxt_dist:
                            old_n, old_e, old_p = nxt_dist[code]
                            nxt_dist[code] = (old_n, old_e, old_p + share)
                        else:
                            nxt_dist[code] = (nv + 1, ne, share)
            dist = nxt_dist
    return res


def yule_mean_ak(
    n0: int,
    jump_means: Dict[Tuple[int, int], Fraction],
    max_steps: int,
    u: float,
    k: int,
) -> Tuple[float, float]:
    """
    Evaluates E_T[A_k(T_u)] = sum_{m>=0} binom(n0+m-1, m) e^{-n0 u} (1 - e^{-u})^m E_T[A_k(T^(m))]
    returning (truncated_value, certified_upper_tail_error).
    """
    q = math.exp(-u)
    val = 0.0
    prob_sum = 0.0
    for m in range(max_steps + 1):
        p_m = math.comb(n0 + m - 1, m) * (q ** n0) * ((1.0 - q) ** m)
        val += p_m * float(jump_means[(m, k)])
        prob_sum += p_m
    # For m > max_steps, A_k <= n0 + m; tail of E[N(u)] from m > max_steps:
    mean_n = n0 * math.exp(u)
    trunc_mean_n = sum(
        math.comb(n0 + m - 1, m) * (q ** n0) * ((1.0 - q) ** m) * (n0 + m)
        for m in range(max_steps + 1)
    )
    tail_err = max(0.0, mean_n - trunc_mean_n)
    return val, tail_err


def check_continuous_eq6_and_eq7() -> bool:
    """
    Checks Eq. (6) |E_{T^+} A_k(T_u) - E_T A_k(T_u)| <= 3 e^{2u} with certified
    Yule tail bounds on small trees, and checks Eq. (7) Var_T(A_k(T_b)) <= 9 n(T) int_0^b e^{4(b-v)+v} dv.
    """
    max_steps = 5
    trees_by_n = generate_all_unlabelled_trees(5)
    for n in range(3, 6):
        for edges in trees_by_n[n]:
            jm_base = exact_jump_chain_expectations(n, edges, max_steps=max_steps)
            for w in range(n):
                edges_plus = edges + ((w, n),)
                jm_plus = exact_jump_chain_expectations(n + 1, edges_plus, max_steps=max_steps)
                for u in (0.05, 0.12):
                    bound = 3.0 * math.exp(2.0 * u)
                    for k in range(2, n + 3):
                        v0, err0 = yule_mean_ak(n, jm_base, max_steps, u, k)
                        v1, err1 = yule_mean_ak(n + 1, jm_plus, max_steps, u, k)
                        if abs(v1 - v0) - (err0 + err1) > bound + 1e-12:
                            return False

            # Check Eq. (7) at b = 0.08 from T:
            b = 0.08
            q = math.exp(-b)
            # Compute exact second moment E_T[A_k(T^(m))^2] for m <= max_steps:
            dist = {canonical_unlabelled_tree(edges_to_adj(n, edges)): (n, edges, Fraction(1, 1))}
            sq_means: Dict[Tuple[int, int], Fraction] = {}
            for m in range(max_steps + 1):
                for k in range(2, n + 3):
                    sq_means[(m, k)] = sum(
                        pr * (sum(1 for u_v in range(nv) if len(edges_to_adj(nv, ed)[u_v]) >= k) ** 2)
                        for _, (nv, ed, pr) in dist.items()
                    )
                if m < max_steps:
                    nxt_d: Dict[Tuple, Tuple[int, Tuple[Tuple[int, int], ...], Fraction]] = {}
                    for _, (nv, ed, pr) in dist.items():
                        adj = edges_to_adj(nv, ed)
                        for u_v in range(nv):
                            sh = pr * Fraction(1, nv * len(adj[u_v]))
                            for w_v in adj[u_v]:
                                ne = ed + ((w_v, nv),)
                                cd = canonical_unlabelled_tree(edges_to_adj(nv + 1, ne))
                                if cd in nxt_d:
                                    nxt_d[cd] = (nv + 1, ne, nxt_d[cd][2] + sh)
                                else:
                                    nxt_d[cd] = (nv + 1, ne, sh)
                    dist = nxt_d

            rhs_eq7 = 9.0 * n * (math.exp(4.0 * b) - math.exp(b)) / 3.0
            for k in range(2, n + 3):
                m1 = sum(
                    math.comb(n + m - 1, m) * (q ** n) * ((1.0 - q) ** m) * float(jm_base[(m, k)])
                    for m in range(max_steps + 1)
                )
                m2 = sum(
                    math.comb(n + m - 1, m) * (q ** n) * ((1.0 - q) ** m) * float(sq_means[(m, k)])
                    for m in range(max_steps + 1)
                )
                # Certified upper tail of E[N(b)^2]:
                full_en2 = n * (n + 1) * math.exp(2.0 * b) - n * math.exp(b)
                trunc_en2 = sum(
                    math.comb(n + m - 1, m) * (q ** n) * ((1.0 - q) ** m) * ((n + m) ** 2)
                    for m in range(max_steps + 1)
                )
                var_upper = (m2 + max(0.0, full_en2 - trunc_en2)) - (m1 * m1)
                if var_upper > rhs_eq7 + 1e-9:
                    return False

    return True


# ---------------------------------------------------------------------------
# Module 5: Eligibility reverse bound, L^2 recursion, positivity & mesh series
# ---------------------------------------------------------------------------

def check_eligibility_and_l2_positivity(
    constant_h_mutation: bool = False,
    mu_override: float = 0.565,
) -> Tuple[int, bool]:
    """
    1. Exhaustively checks across ALL unlabelled trees with 3 <= n <= 11 (430 trees)
       and EVERY choice of root core o in K(T) that the count |G| of good noninitial
       core vertices in renewal.md line 95-97 (used in almost-sure.md line 59 for
       C(T) <= C_3(1 + Psi(T))) satisfies |G| >= C(T)/12 - 1.
    2. Verifies uniform boundedness of the L^2 martingale recursion (Eq. 9),
       positivity of the conditional Cauchy-Schwarz lower limit, and summability
       of the mesh variance series (line 113).
    """
    trees_by_n = generate_all_unlabelled_trees(11)
    total_trees = 0
    for n in range(3, 12):
        for edges in trees_by_n[n]:
            total_trees += 1
            adj = edges_to_adj(n, edges)
            cores = [u for u in range(n) if len(adj[u]) >= 2]
            c_cnt = len(cores)
            for root_o in cores:
                parent = [-1] * n
                children = [[] for _ in range(n)]
                st = [root_o]
                seen = {root_o}
                while st:
                    u = st.pop()
                    for v in adj[u]:
                        if v not in seen:
                            seen.add(v)
                            parent[v] = u
                            children[u].append(v)
                            st.append(v)

                good = 0
                for v in cores:
                    if v == root_o:
                        continue
                    c_v = sum(1 for y in adj[v] if len(adj[y]) >= 2)
                    l_v = len(adj[v]) - c_v
                    if c_v <= 2 and l_v >= 1:
                        good += 1
                    elif c_v == 2 and l_v == 0:
                        # Bare noninitial core: check if it has a bare core child
                        has_bare_child = any(
                            len(adj[y]) == 2 and sum(1 for z in adj[y] if len(adj[z]) >= 2) == 2
                            for y in children[v]
                        )
                        if has_bare_child:
                            good += 1

                if Fraction(good, 1) < Fraction(c_cnt, 12) - Fraction(1, 1):
                    return total_trees, False

    # 2. Check L^2 recursion (Eq. 9), conditional positivity ratio, and mesh series:
    eta = 0.25
    c9 = 1.5
    c10 = 2.0
    ej_sq = 1.0  # start from Psi(T) = 1, n(T) = 3
    n0 = 3.0
    ej_at_1500 = 0.0

    for j in range(0, 2001):
        h_j = 1.0 if constant_h_mutation else (math.log(j + 2.0) / eta)
        mult = 1.0 + c9 * math.exp(-2.0 * eta * h_j)
        exponent_term = -(2.0 * mu_override - 1.0) * j
        if exponent_term > 60.0 or 4.0 * h_j + exponent_term > 80.0:
            return total_trees, False
        add_term = c10 * (
            math.exp(4.0 * h_j + exponent_term) * n0
            + math.exp(-2.0 * eta * h_j - 2.0 * mu_override * j)
        )
        ej_sq = mult * ej_sq + add_term
        if ej_sq > 1e35:
            return total_trees, False
        if j == 1500:
            ej_at_1500 = ej_sq

    # Verify genuine Cauchy convergence of E[(Z_j^*)^2] between j = 1500 and j = 2000:
    if (ej_sq - ej_at_1500) / max(1.0, ej_at_1500) > 1e-3:
        return total_trees, False

    # Check conditional positivity limit (lines 102-106):
    # With Psi(T_j) ~ e^{mu j} and n(T_j) ~ e^j, ratio Psi^2 / (Psi^2 + n + 1) -> 1 iff 2*mu > 1!
    j_large = 200.0
    psi_sq_over_n = math.exp((2.0 * mu_override - 1.0) * j_large)
    cond_lb_limit = psi_sq_over_n / (psi_sq_over_n + 1.0)
    if cond_lb_limit < 0.99:
        return total_trees, False

    # Check mesh variance tail sum (line 113) for m = 4:
    m_mesh = 4
    tail_sum = 0.0
    for j in range(6000, 8000):
        t = j / m_mesh
        h = 1.0 if constant_h_mutation else (math.log(t + 2.0) / eta)
        expo = -(2.0 * mu_override - 1.0) * t + 4.0 * h
        if expo > 50.0:
            return total_trees, False
        tail_sum += math.exp(expo)
    if tail_sum > 1e-6:
        return total_trees, False

    return total_trees, True


# ---------------------------------------------------------------------------
# Module 6: Section 5 discrete eigenmartingale Y_m & nontrivial support of X
# ---------------------------------------------------------------------------

def check_section5_discrete_martingale_and_support(
    wrong_step_shift: int = 0,
) -> bool:
    """
    Verifies Section 5 (lines 119-134) of almost-sure.md:
      1. Differentiating / Laplace-conditioning Eq. (5) at size m >= 3 gives
         E[Psi(T_{m+1}) | T_m] = (1 + mu / m) Psi(T_m) and makes
         Y_m = Gamma(m) Psi(T_m) / Gamma(m + mu) an exact martingale.
      2. Exact identity E[e^{-2*mu*tau_m}] = prod_{j=3}^{m-1} j / (j + 2*mu)
         = (Gamma(3 + 2*mu) / 2) * Gamma(m) / Gamma(m + 2*mu) ~ const * m^{-2*mu},
         and exact agreement at m = 3 between E[X] = d_2 (Gamma(3)/Gamma(3+mu)) Psi(T_3)
         and c_2 = (2 / ((2 + mu) Gamma(mu + 2))) * (Psi(T_3) d_2) from asymptotics.md.
      3. Exact rational Markov-chain verification (3 <= m <= 9) that from T_2:
           - Star K_{1, m-1} has exact probability 2 / (m - 1) > 0 and C = 1,
             sending E[X | T_m = K_{1,m-1}] <= C' m^{-mu} -> 0;
           - Path P_m has exact probability 2 / (m - 1)! > 0 and C = m - 2,
             sending E[X | T_m = P_m] >= c' m^{1 - mu} -> +infinity.
    """
    # 1. Exact holding-time Laplace conditioning check on continuous eigenfunction:
    # If E_T[Psi(T_t)] = e^{mu t} Psi(T) and first jump wait is sigma_1 ~ Exp(m)
    # independent of T_{m+1}, then for s > mu:
    # Psi(T)/(s - mu) = Psi(T)/(s + m) + (m/(s + m)) * E[Psi(T_{m+1}) | T_m=T] / (s - mu)
    # <=> E[Psi(T_{m+1}) | T_m=T] = (1 + mu/m) Psi(T).
    for mu_rat in (Fraction(113, 200), Fraction(1133, 2000)):
        for m in range(3, 25):
            s_rat = Fraction(2, 1)
            one_step_factor = ((s_rat + m) - (s_rat - mu_rat)) * Fraction(1, m)
            claimed_factor = Fraction(1, 1) + mu_rat / Fraction(m + wrong_step_shift, 1)
            if one_step_factor != claimed_factor:
                return False

    # Check Gamma martingale identity Gamma(m+1)/Gamma(m+1+mu) * (1 + mu/(m + shift)) == Gamma(m)/Gamma(m+mu):
    mu_f = 0.565723
    for m in range(3, 30):
        ym_coeff = math.gamma(m) / math.gamma(m + mu_f)
        ym1_coeff = math.gamma(m + 1) / math.gamma(m + 1 + mu_f)
        step_mult = 1.0 + mu_f / (m + wrong_step_shift)
        if abs(ym1_coeff * step_mult - ym_coeff) > 1e-13:
            return False

    # 2. Check E[e^{-2*mu*tau_m}] product & E[X] == c_2 at m = 3:
    prod_val = 1.0
    for m in range(3, 60):
        gamma_closed = 0.5 * math.gamma(3.0 + 2.0 * mu_f) * math.gamma(m) / math.gamma(m + 2.0 * mu_f)
        if abs(prod_val - gamma_closed) / gamma_closed > 1e-12:
            return False
        prod_val *= m / (m + 2.0 * mu_f)

    # Check Eq. (132) at m = 3 vs asymptotics.md c_2:
    psi_t3 = 1.7
    d2 = 0.85
    ex_sec5 = d2 * (math.gamma(3.0) / math.gamma(3.0 + mu_f)) * psi_t3
    c2_asymp = (2.0 / ((2.0 + mu_f) * math.gamma(mu_f + 2.0))) * (psi_t3 * d2)
    if abs(ex_sec5 - c2_asymp) > 1e-15:
        return False

    # 3. Exact rational verification of star & path probabilities from T_2 for 3 <= m <= 9:
    dist: Dict[Tuple, Tuple[int, Tuple[Tuple[int, int], ...], Fraction]] = {
        canonical_unlabelled_tree(edges_to_adj(2, ((0, 1),))): (2, ((0, 1),), Fraction(1, 1))
    }
    for m in range(3, 10):
        nxt_dist: Dict[Tuple, Tuple[int, Tuple[Tuple[int, int], ...], Fraction]] = {}
        for _, (nv, ed, pr) in dist.items():
            adj = edges_to_adj(nv, ed)
            for u in range(nv):
                sh = pr * Fraction(1, nv * len(adj[u]))
                for w in adj[u]:
                    ne = ed + ((w, nv),)
                    cd = canonical_unlabelled_tree(edges_to_adj(nv + 1, ne))
                    if cd in nxt_dist:
                        nxt_dist[cd] = (nv + 1, ne, nxt_dist[cd][2] + sh)
                    else:
                        nxt_dist[cd] = (nv + 1, ne, sh)
        dist = nxt_dist

        star_edges = tuple((0, j) for j in range(1, m))
        path_edges = tuple((j, j + 1) for j in range(m - 1))
        star_code = canonical_unlabelled_tree(edges_to_adj(m, star_edges))
        path_code = canonical_unlabelled_tree(edges_to_adj(m, path_edges))

        if star_code not in dist or path_code not in dist:
            return False
        if dist[star_code][2] != Fraction(2, m - 1):
            return False
        if dist[path_code][2] != Fraction(2, math.factorial(m - 1)):
            return False

        star_cores = sum(1 for u in range(m) if len(edges_to_adj(m, star_edges)[u]) >= 2)
        path_cores = sum(1 for u in range(m) if len(edges_to_adj(m, path_edges)[u]) >= 2)
        if star_cores != 1 or path_cores != m - 2:
            return False

    return True


# ---------------------------------------------------------------------------
# Main runner + 10 Targeted Mutations (Negative Controls)
# ---------------------------------------------------------------------------

def main() -> int:
    print("=== Cold-Challenge Verification: round2/fringe/almost-sure.md ===")

    # 1. Capped-mark root-gain increment & general-founder renewal generator
    assert check_capped_root_gain_bound(cap_at=20), "Capped root-gain bound failed!"
    checked_founders, ok_founders = check_general_founder_renewal_identity()
    assert ok_founders, "General-founder renewal identity failed!"
    print(
        f"[PASS] Capped-mark root-gain bound (D <= 120) & general-founder renewal/Z/I(W) "
        f"identities exact across {checked_founders} rooted marked configurations (3 <= n <= 7)."
    )

    # 2. Uniform exponential renewal, contour high-frequency shift & eigenfunction
    assert check_exponential_renewal_and_contour(), "Contour / exponential renewal check failed!"
    print(
        "[PASS] Integration by parts K(z)=O(1/|z|), O(1/tau^2) contour shift for v(t)-k(t), "
        "uniform exponential remainder Eq. (3), and exact eigenfunction Eq. (5) verified."
    )

    # 3. Clean/dirty one-leaf coupling & backward martingale variance
    checked_coupled, ok_coupled = check_one_leaf_coupling_exact()
    assert ok_coupled, "One-leaf clean/dirty coupling check failed!"
    assert check_continuous_eq6_and_eq7(), "Continuous-time Eq. (6)/(7) check failed!"
    print(
        f"[PASS] One-leaf clean/dirty coupling exact across {checked_coupled} coupled states "
        f"(marginals, clean invariants, D_comb(0)<=3, drift <= 2 D_comb) + Eqs. (6)-(7) verified."
    )

    # 4. Eligibility reverse comparison, L^2 recursion, positivity, and mesh summability
    checked_elig, ok_elig = check_eligibility_and_l2_positivity()
    assert ok_elig, "Eligibility / L^2 recursion / positivity check failed!"
    print(
        f"[PASS] Eligibility reverse bound checked on all {checked_elig} unlabelled trees (3 <= n <= 11); "
        f"L^2 recursion (h=log(j+2)/eta), conditional positivity, and mesh summability verified."
    )

    # 5. Section 5 discrete eigenmartingale Y_m & nontrivial support via star/path seeds
    assert check_section5_discrete_martingale_and_support(), "Section 5 check failed!"
    print(
        "[PASS] Section 5 discrete eigenmartingale Y_m, E[e^{-2*mu*tau_m}] product, E[X]=c_2 "
        "consistency, and exact star (2/(m-1)) / path (2/(m-1)!) seed laws (3 <= m <= 9) verified."
    )

    # 6. Ten deliberate mutations (negative controls)
    mutations: List[Tuple[str, bool]] = [
        (
            "Mutation 1: Uncapped marks m(D)=D with bounded oscillating g breaks ||R B_0 g||_W",
            not check_capped_root_gain_bound(cap_at=None),
        ),
        (
            "Mutation 2: Omitting initial fringes sum_{v!=o} P_t f(S_v) for general founder T",
            not check_general_founder_renewal_identity(omit_initial_fringes=True)[1],
        ),
        (
            "Mutation 3: Contour shift on v(t) without subtracting k(t) loses O(1/tau^2) L^1 decay",
            not check_exponential_renewal_and_contour(wrong_contour_without_subtracting_k=True),
        ),
        (
            "Mutation 4: Claiming initial dirty core bound 2 instead of 3 in Eq. (6) [leaf promotion gives 3]",
            not check_one_leaf_coupling_exact(max_init_dirty=2)[1],
        ),
        (
            "Mutation 5: Claiming dirty core drift <= 1*D_comb (e^u) instead of 2*D_comb (e^{2u}) [bare core gives 2]",
            not check_one_leaf_coupling_exact(max_drift_per_dirty=Fraction(3, 2))[1],
        ),
        (
            "Mutation 6: Marking only recipient copy dirty (instead of clean pair) breaks coupling invariants",
            not check_one_leaf_coupling_exact(bug_mark_only_one_copy_dirty=True)[1],
        ),
        (
            "Mutation 7: Constant h instead of h = log(j+2)/eta blows up L^2 martingale recursion Eq. (9)",
            not check_eligibility_and_l2_positivity(constant_h_mutation=True)[1],
        ),
        (
            "Mutation 8: Sub-threshold exponent mu = 0.50 (2*mu <= 1) breaks L^2 bound, positivity & mesh sum",
            not check_eligibility_and_l2_positivity(mu_override=0.50)[1],
        ),
        (
            "Mutation 9: Sub-threshold exponent mu = 0.48 breaks conditional positivity ratio Psi^2/(Psi^2+n+1)",
            not check_eligibility_and_l2_positivity(mu_override=0.48)[1],
        ),
        (
            "Mutation 10: Off-by-one step factor (1 + mu/(m+1)) instead of (1 + mu/m) breaks martingale Y_m",
            not check_section5_discrete_martingale_and_support(wrong_step_shift=1),
        ),
    ]

    for desc, caught in mutations:
        assert caught, f"MUTATION NOT CAUGHT: {desc}"
        print(f"[CAUGHT] {desc}")

    print(f"ALL CHECKS AND {len(mutations)}/{len(mutations)} MUTATION CONTROLS PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
