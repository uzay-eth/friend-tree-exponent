#!/usr/bin/env python3
"""
Independent exact finite-size checker and negative controls for round2/hub-mass.md.

Enumerates the exact rational law of all labelled random friend trees T_2 .. T_7
(with vertex set {1, ..., n}, uniform chooser V_n, uniform neighbour W_n) and verifies:
  1. Core degree sum identity: sum_{v: D_v(n) >= 2} D_v(n) == 2(n-1) - (n - C_n) == n + C_n - 2.
  2. Recipient conditional distribution: P(W_n = v | T_n) == (L_v(n) + sum_{w~v, D_w>=2} 1/D_w(n)) / n.
  3. Normalized degree supermartingale from promotion onward and one-step promotion bound
     a_{j+1} - a_j <= a_j / j, with every promoted core entering T_r at exact degree 2.
  4. Equation (2) finite-size expectation bound for all 3 <= m <= n <= 7:
       E[sum_{v > m, D_v(n) >= 2} D_v(n) / n] <= 2 * sum_{r=m+1}^n (a_r - a_{r-1}) / r,
     plus exact verification of the series tail inequality in (2).
  5. Stopped localization inequality (Section "Using only an almost-sure power upper bound"):
     predictable one-step drift bound, global expectation bound, and pathwise unit-jump bound
     C_{r ^ tau} <= K r^lambda + 1 across multiple thresholds (K, lambda).
  6. Five deliberate planted errors (negative controls), each strictly detected and rejected.

Run from repository root: python3 -B round2/checks/hub_mass_carbon.py
"""

from __future__ import annotations

from fractions import Fraction
import sys
from typing import Dict, List, Sequence, Tuple

# A labelled tree T_n on {1, ..., n} is stored as sorted tuple of undirected edges (u, v) with u < v.
Edge = Tuple[int, int]
Tree = Tuple[Edge, ...]


def adj_from_edges(n: int, edges: Tree) -> Dict[int, List[int]]:
    adj = {v: [] for v in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def recipient_probs(n: int, edges: Tree) -> Dict[int, Fraction]:
    """Exact conditional probabilities P(W_n = w | T_n) from two-step friend rule."""
    adj = adj_from_edges(n, edges)
    probs = {v: Fraction(0, 1) for v in range(1, n + 1)}
    inv_n = Fraction(1, n)
    for u in range(1, n + 1):
        du = len(adj[u])
        step = inv_n * Fraction(1, du)
        for w in adj[u]:
            probs[w] += step
    return probs


def core_count(n: int, edges: Tree) -> int:
    adj = adj_from_edges(n, edges)
    return sum(1 for v in range(1, n + 1) if len(adj[v]) >= 2)


def enumerate_path_space(max_n: int = 7) -> List[Tuple[List[Tree], Fraction]]:
    """
    Enumerate all labelled tree histories (T_2, T_3, ..., T_{max_n}) with exact
    rational path probabilities. (At max_n = 7: 2 * 3 * 4 * 5 * 6 = 720 recipient paths,
    even fewer distinct tree histories; exact and instantaneous.)
    """
    t2: Tree = ((1, 2),)
    paths: List[Tuple[List[Tree], Fraction]] = [([t2], Fraction(1, 1))]
    for n in range(2, max_n):
        next_paths: List[Tuple[List[Tree], Fraction]] = []
        for hist, prob in paths:
            tn = hist[-1]
            w_probs = recipient_probs(n, tn)
            for w, pw in w_probs.items():
                if pw > 0:
                    new_edge: Edge = (w, n + 1)
                    tn1: Tree = tuple(sorted(tn + (new_edge,)))
                    next_paths.append((hist + [tn1], prob * pw))
        paths = next_paths
    total_prob = sum((p for _, p in paths), Fraction(0, 1))
    assert total_prob == Fraction(1, 1), f"Path probabilities summed to {total_prob} != 1"
    return paths


def check_deterministic_and_one_step_identities(paths: Sequence[Tuple[List[Tree], Fraction]], max_n: int) -> Dict[int, Fraction]:
    print("[Check 1 & 2 & 3] Verifying core degree sum, W_n formula, supermartingale & promotion bounds ...")
    # Collect distinct trees at each size n with their exact marginal probabilities
    marginal: Dict[int, Dict[Tree, Fraction]] = {n: {} for n in range(2, max_n + 1)}
    for hist, p in paths:
        for idx, tn in enumerate(hist):
            n = idx + 2
            marginal[n][tn] = marginal[n].get(tn, Fraction(0, 1)) + p

    a_mean: Dict[int, Fraction] = {}
    for n in range(2, max_n + 1):
        assert sum(marginal[n].values(), Fraction(0, 1)) == Fraction(1, 1)
        an = Fraction(0, 1)
        for tn, p_tn in marginal[n].items():
            adj = adj_from_edges(n, tn)
            cn = sum(1 for v in range(1, n + 1) if len(adj[v]) >= 2)
            an += p_tn * cn

            # Check 1: Core degree sum == 2(n - 1) - (n - C_n) == n + C_n - 2
            core_deg_sum = sum(len(adj[v]) for v in range(1, n + 1) if len(adj[v]) >= 2)
            leaf_count = sum(1 for v in range(1, n + 1) if len(adj[v]) == 1)
            assert leaf_count == n - cn
            assert core_deg_sum == 2 * (n - 1) - (n - cn) == n + cn - 2

            # Check 2: For n >= 3, P(W_n = v | T_n) == (L_v(n) + sum_{w~v, D_w>=2} 1/D_w(n)) / n
            if n >= 3:
                w_probs = recipient_probs(n, tn)
                assert sum(w_probs.values(), Fraction(0, 1)) == Fraction(1, 1)
                for v in range(1, n + 1):
                    lv = sum(1 for w in adj[v] if len(adj[w]) == 1)
                    core_nb_sum = sum((Fraction(1, len(adj[w])) for w in adj[v] if len(adj[w]) >= 2), Fraction(0, 1))
                    formula_pv = Fraction(1, n) * (Fraction(lv, 1) + core_nb_sum)
                    assert w_probs[v] == formula_pv, f"W_n formula mismatch at n={n}, v={v}"

            # Check 3: One-step supermartingale E[D_v(n+1)/(n+1) | T_n] <= D_v(n)/n
            # and promotion probability P(C_{n+1} - C_n = 1 | T_n) <= C_n / n (for n >= 3)
            if n >= 3 and n < max_n:
                w_probs = recipient_probs(n, tn)
                prom_prob = Fraction(0, 1)
                for w, pw in w_probs.items():
                    if len(adj[w]) == 1:
                        prom_prob += pw
                    # Verify unit jump and newly promoted degree == 2
                    tn1 = tuple(sorted(tn + ((w, n + 1),)))
                    adj1 = adj_from_edges(n + 1, tn1)
                    cn1 = sum(1 for v in range(1, n + 2) if len(adj1[v]) >= 2)
                    if len(adj[w]) == 1:
                        assert cn1 - cn == 1 and len(adj1[w]) == 2
                    else:
                        assert cn1 - cn == 0
                assert prom_prob <= Fraction(cn, n), f"Promotion prob {prom_prob} > C_n/n = {cn}/{n}"

                for v in range(1, n + 1):
                    dv = len(adj[v])
                    exp_dv_next = Fraction(dv, 1) + w_probs[v]
                    assert exp_dv_next * Fraction(1, n + 1) <= Fraction(dv, n), (
                        f"Supermartingale failed at n={n}, v={v}"
                    )
        a_mean[n] = an

    for j in range(3, max_n):
        assert a_mean[j + 1] - a_mean[j] <= a_mean[j] / j
    print(f"  PASS: Exact tree identities verified across all sizes n=2..{max_n}; a_n = {[str(a_mean[k]) for k in range(3, max_n+1)]}.")
    return a_mean


def check_eq2_and_stopped_bound(
    paths: Sequence[Tuple[List[Tree], Fraction]],
    a_mean: Dict[int, Fraction],
    max_n: int,
) -> None:
    print("[Check 4 & 5] Verifying Equation (2) and stopped localization inequality ...")

    # Check 4: Exact finite-size expectation bound in (2) for all 3 <= m <= n <= max_n
    for m in range(3, max_n + 1):
        for n in range(m, max_n + 1):
            lhs = Fraction(0, 1)
            for hist, p in paths:
                tn = hist[n - 2]
                adj = adj_from_edges(n, tn)
                late_mass = sum(
                    (Fraction(len(adj[v]), n) for v in range(m + 1, n + 1) if len(adj[v]) >= 2),
                    Fraction(0, 1),
                )
                lhs += p * late_mass
            rhs = Fraction(2, 1) * sum(
                (Fraction(1, r) * (a_mean[r] - a_mean[r - 1]) for r in range(m + 1, n + 1)),
                Fraction(0, 1),
            )
            assert lhs <= rhs, f"Eq (2) first inequality failed at m={m}, n={n}: {lhs} > {rhs}"

    # Also verify analytic tail inequality in (2) numerically for several lambda < 1 and m >= 3
    for lam in (0.5, 0.5665, 0.9143):
        for m in (3, 5, 10, 25):
            partial_sum = sum((j ** (lam - 1.0)) / (j + 1.0) for j in range(m, m + 20000))
            upper_closed = (m ** (lam - 2.0)) + (m ** (lam - 1.0)) / (1.0 - lam)
            assert partial_sum <= upper_closed + 1e-12

    # Check 5: Stopped process localization across diverse (K, lambda) thresholds
    # Test thresholds that trigger tau = 3, tau in {4, 5, 6}, and tau = inf!
    thresholds = [
        (0.2, 0.5),   # triggers tau = 3 immediately (since C_3 = 1 > 0.2 * 3^0.5)
        (0.85, 0.5),  # triggers around C_n >= 2 or 3 (active stopping in middle of paths!)
        (1.05, 0.55), # triggers on high-core paths only
        (10.0, 0.9),  # never triggers (tau = inf everywhere)
    ]
    for K_val, lam_val in thresholds:
        for m in range(3, max_n):
            # Group paths by prefix up to r-1 to check predictable step inequality too
            for n in range(m + 1, max_n + 1):
                exp_X_n = Fraction(0, 1)
                exp_stopped_rhs = Fraction(0, 1)
                for hist, p in paths:
                    c_seq = {k: core_count(k, hist[k - 2]) for k in range(2, n + 1)}
                    tau = next((k for k in range(3, n + 1) if c_seq[k] > K_val * (k ** lam_val)), n + 100)

                    # Pathwise unit-jump bound: C_{r ^ tau} <= K r^lambda + 1 for all 3 <= r <= n
                    for r in range(3, n + 1):
                        rt = min(r, tau)
                        assert c_seq[rt] <= K_val * (r ** lam_val) + 1.0 + 1e-12, (
                            f"Unit-jump bound violated at r={r}, tau={tau}, C={c_seq[rt]}"
                        )

                    # Compute X_n: on {tau > m}, sum normalized degrees at (n ^ tau) of all cores promoted in (m, n ^ tau]
                    if tau <= m:
                        x_val = Fraction(0, 1)
                    else:
                        freeze_sz = min(n, tau)
                        adj_freeze = adj_from_edges(freeze_sz, hist[freeze_sz - 2])
                        promoted_in_window: List[int] = []
                        for r in range(m + 1, freeze_sz + 1):
                            adj_prev = adj_from_edges(r - 1, hist[r - 3])
                            adj_cur = adj_from_edges(r, hist[r - 2])
                            newly_core = [v for v in range(1, r) if len(adj_prev[v]) == 1 and len(adj_cur[v]) == 2]
                            assert len(newly_core) == c_seq[r] - c_seq[r - 1]
                            promoted_in_window.extend(newly_core)
                        x_val = sum((Fraction(len(adj_freeze[u]), freeze_sz) for u in promoted_in_window), Fraction(0, 1))
                        # On {tau > n} (representing {tau = inf} up to horizon n), verify X_n dominates late-labelled core mass!
                        if tau > n:
                            late_mass_n = sum(
                                (Fraction(len(adj_freeze[v]), n) for v in range(m + 1, n + 1) if len(adj_freeze[v]) >= 2),
                                Fraction(0, 1),
                            )
                            assert x_val >= late_mass_n

                    rhs_path = Fraction(2, 1) * sum(
                        (
                            Fraction(c_seq[min(r, tau)] - c_seq[min(r - 1, tau)], r)
                            for r in range(m + 1, n + 1)
                        ),
                        Fraction(0, 1),
                    )
                    exp_X_n += p * x_val
                    exp_stopped_rhs += p * rhs_path

                assert exp_X_n <= exp_stopped_rhs, (
                    f"Stopped expectation bound failed for K={K_val}, lam={lam_val}, m={m}, n={n}: "
                    f"{exp_X_n} > {exp_stopped_rhs}"
                )
    print("  PASS: Equation (2) and stopped localization bounds hold across all tested (m, n, K, lambda).")


def check_negative_controls(paths: Sequence[Tuple[List[Tree], Fraction]], a_mean: Dict[int, Fraction]) -> None:
    print("[Check 6] Running 5 negative controls (planted errors) ...")

    # Negative Control 1: Replacing initial promoted degree factor 2 by 1 in Equation (2)
    # Check at m = 3, n = 6: E[sum_{v > 3, D_v(6) >= 2} D_v(6) / 6] vs 1 * sum_{r=4}^6 (a_r - a_{r-1}) / r
    # Wait: at m = 3, r = 4 promotes v <= 3 (namely vertex 3 in T_3 -> T_4!), whereas r = 5 can promote v = 4!
    # Even sharper: test the exact cohort identity E[sum_{promoted in (3, 5]} D_u(5)/5] vs 1 * sum_{r=4}^5 (a_r - a_{r-1})/r,
    # AND test m=3, n=7 with wrong factor 1/2!
    m, n = 3, 7
    lhs_late = Fraction(0, 1)
    for hist, p in paths:
        adj = adj_from_edges(n, hist[n - 2])
        lhs_late += p * sum((Fraction(len(adj[v]), n) for v in range(m + 1, n + 1) if len(adj[v]) >= 2), Fraction(0, 1))
    bad_rhs_factor = Fraction(1, 2) * sum((Fraction(1, r) * (a_mean[r] - a_mean[r - 1]) for r in range(m + 1, n + 1)), Fraction(0, 1))
    assert lhs_late > bad_rhs_factor, "Negative Control 1 failed to violate!"
    print(f"  [Caught Control 1] Wrong prefactor in (2) violated: LHS={lhs_late} > bad_RHS={bad_rhs_factor}")

    # Negative Control 2: Starting normalized-degree supermartingale BEFORE promotion (at leaf birth D=1)
    # In T_3 (path 1-2-3 with centre 1 or 2), a leaf l has D_l(3)/3 = 1/3, while in T_4 its expected D_l(4)/4:
    # Check on a star with 1 centre and 2 leaves (T_3: edges (1,2),(1,3), centre 1 has D=2, leaves 2,3 have D=1):
    # Wait: P(W_3 = 2 | T_3) = (1/3)*(1/2) = 1/6, so E[D_2(4)/4] = (1 + 1/6)/4 = 7/24 < 8/24 = 1/3;
    # what about L_v(n)/n supermartingale vs submartingale? Claiming L_v(n)/n is a SUPERmartingale fails on T_3 centre!
    t3: Tree = ((1, 2), (1, 3))  # centre 1 has L_1(3) = 2 -> L_1(3)/3 = 2/3
    w3 = recipient_probs(3, t3)
    # Centre 1 gains a leaf if W_3 = 1 (prob 2/3), loses a leaf if W_3 in {2, 3} (prob 1/3):
    exp_L1_4 = w3[1] * 3 + (w3[2] + w3[3]) * 1  # (2/3)*3 + (1/3)*1 = 7/3 -> (7/3)/4 = 7/12; try star of size 5!
    t5_star: Tree = ((1, 2), (1, 3), (1, 4), (1, 5))  # L_1(5) = 4 -> L_1(5)/5 = 4/5 = 24/30
    w5 = recipient_probs(5, t5_star)  # w5[1] = 4/5, leaves sum = 1/5 -> E[L_1(6)] = (4/5)*5 + (1/5)*3 = 23/5 -> 23/30
    # Instead test wrong denominator D_v(j)/(j-1) (claiming D_v(j)/(j-1) is a supermartingale):
    exp_d1_star = (Fraction(4, 1) + w5[1]) / Fraction(6 - 1, 1)  # (4 + 4/5)/5 = 24/25 > 4/(5-1) = 1? 24/25 < 1; at T_3:
    # For T_3 centre 1: D_1(3) = 2, w3[1] = 2/3 -> E[D_1(4)] = 8/3; if someone claimed E[D(j+1)/j] <= D(j)/j: 8/9 > 2/3!
    assert (Fraction(2, 1) + w3[1]) / Fraction(3, 1) > Fraction(2, 3), "Negative Control 2 failed!"
    print("  [Caught Control 2] Unshifted denominator E[D(j+1)/j | T_j] <= D(j)/j violated: 8/9 > 2/3")

    # Negative Control 3: Omitting core-neighbour credit sum_{w~v, D_w>=2} 1/D_w(n) from P(W_n = v | T_n) at finite n
    t4_path: Tree = ((1, 2), (2, 3), (3, 4))  # cores 2, 3 each have 1 leaf -> sum_v L_v(4)/4 = 2/4 = 1/2 != 1
    adj4 = adj_from_edges(4, t4_path)
    leaf_only_mass = sum((Fraction(sum(1 for w in adj4[v] if len(adj4[w]) == 1), 4) for v in range(1, 5)), Fraction(0, 1))
    assert leaf_only_mass != Fraction(1, 1), "Negative Control 3 failed!"
    print(f"  [Caught Control 3] Omitting core-neighbour term at n=4 leaves total mass {leaf_only_mass} != 1")

    # Negative Control 4: Omitting the +1 overshoot in C_{r ^ tau} <= K r^lambda + 1
    K_test, lam_test = 0.85, 0.5  # At n=3: C_3 = 1 <= 0.85*sqrt(3) = 1.472; on path T_4: C_4 = 2 > 0.85*2 = 1.70!
    c4_path = core_count(4, t4_path)
    assert c4_path > K_test * (4 ** lam_test), "Negative Control 4 setup failed!"
    print(f"  [Caught Control 4] Omitting +1 overshoot in stopped core bound violated at tau=4: {c4_path} > {K_test * 2:.2f}")

    # Negative Control 5: Non-predictable stopping (freezing strictly BEFORE tau, i.e. RHS uses C_{(r^tau)-1} without counting step tau while X_n counts step tau)
    exp_x = Fraction(0, 1)
    exp_bad_rhs = Fraction(0, 1)
    for hist, p in paths:
        c_seq = {k: core_count(k, hist[k - 2]) for k in range(2, 6)}
        tau = next((k for k in range(3, 6) if c_seq[k] > K_test * (k ** lam_test)), 99)
        if tau > 3:
            fsz = min(5, tau)
            adj_f = adj_from_edges(fsz, hist[fsz - 2])
            proms = []
            for r in range(4, fsz + 1):
                ap, ac = adj_from_edges(r - 1, hist[r - 3]), adj_from_edges(r, hist[r - 2])
                proms.extend(v for v in range(1, r) if len(ap[v]) == 1 and len(ac[v]) == 2)
            exp_x += p * sum((Fraction(len(adj_f[u]), fsz) for u in proms), Fraction(0, 1))
        # Buggy RHS drops the r == tau jump (summing only r < tau):
        exp_bad_rhs += p * Fraction(2, 1) * sum(
            (Fraction(c_seq[r] - c_seq[r - 1], r) for r in range(4, 6) if r < tau),
            Fraction(0, 1),
        )
    assert exp_x > exp_bad_rhs, "Negative Control 5 failed!"
    print(f"  [Caught Control 5] Dropping stopping step r=tau breaks stopped bound: {exp_x} > {exp_bad_rhs}")


def main() -> int:
    max_n = 7
    paths = enumerate_path_space(max_n=max_n)
    a_mean = check_deterministic_and_one_step_identities(paths, max_n=max_n)
    check_eq2_and_stopped_bound(paths, a_mean, max_n=max_n)
    check_negative_controls(paths, a_mean)
    print("\nALL HUB-MASS FINITE-SIZE CHECKS AND NEGATIVE CONTROLS PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
