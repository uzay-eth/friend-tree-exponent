#!/usr/bin/env python3
"""Exact enumeration of IR, DAN, and PAN trees on labelled histories from L_3 up to N=8,
plus high-precision verification of analytical roots mu_DAN and mu_PAN.
"""
from collections import defaultdict
from fractions import Fraction
import math
import numpy as np
from scipy import integrate, optimize


def canonical_degree_state(adj):
    """Canonical unlabeled tree key via AHU tree isomorphism algorithm."""
    n = len(adj)
    deg = [len(adj[u]) for u in range(n)]
    # Find center(s) of tree
    rem = set(range(n))
    cur_deg = deg[:]
    leaves = [u for u in range(n) if cur_deg[u] <= 1]
    while len(rem) > 2:
        new_leaves = []
        for u in leaves:
            rem.remove(u)
            for v in adj[u]:
                if v in rem:
                    cur_deg[v] -= 1
                    if cur_deg[v] == 1:
                        new_leaves.append(v)
        leaves = new_leaves
    centers = sorted(rem)

    def encode(u, parent):
        children = [encode(v, u) for v in adj[u] if v != parent]
        children.sort()
        return "(" + "".join(children) + ")"

    codes = [encode(c, -1) for c in centers]
    return min(codes)


def transition_probs(adj, model):
    """Return list of (attach_target_vertex, Fraction_prob) for adding vertex n."""
    n = len(adj)
    deg = [len(adj[u]) for u in range(n)]
    leaves_of = [[] for _ in range(n)]
    for u in range(n):
        for v in adj[u]:
            if deg[v] == 1:
                leaves_of[u].append(v)

    probs = defaultdict(Fraction)
    if model == "IR":
        inv_n = Fraction(1, n)
        for v in range(n):
            p = inv_n / deg[v]
            for w in adj[v]:
                probs[w] += p
    elif model == "DAN":
        inv_n = Fraction(1, n)
        for v in range(n):
            if deg[v] == 1:
                # leaf -> unique neighbor
                w = adj[v][0]
                probs[w] += inv_n
            elif len(leaves_of[v]) > 0:
                # rank-1 -> uniform among its leaf neighbors
                p = inv_n / len(leaves_of[v])
                for w in leaves_of[v]:
                    probs[w] += p
            else:
                # nucleus -> attach directly to v
                probs[v] += inv_n
    elif model == "PAN":
        non_nuc = [v for v in range(n) if deg[v] == 1 or len(leaves_of[v]) > 0]
        inv_d = Fraction(1, len(non_nuc))
        for v in non_nuc:
            if deg[v] == 1:
                w = adj[v][0]
                probs[w] += inv_d
            else:
                p = inv_d / len(leaves_of[v])
                for w in leaves_of[v]:
                    probs[w] += p
    else:
        raise ValueError(model)
    assert sum(probs.values()) == 1
    return probs


def enumerate_model(model, max_n=8):
    # Initial condition: 3-path L_3: 0-1-2
    init_adj = ((1,), (0, 2), (1,))
    init_key = canonical_degree_state(init_adj)
    states = {init_key: (init_adj, Fraction(1, 1))}

    results = {}
    for n in range(3, max_n + 1):
        c_dist = defaultdict(Fraction)
        q_line = Fraction(0, 1)
        for key, (adj, prob) in states.items():
            c = sum(1 for nbrs in adj if len(nbrs) >= 2)
            c_dist[c] += prob
            if c == n - 2 and sum(1 for nbrs in adj if len(nbrs) == 1) == 2:
                q_line += prob
        mean_c = sum(Fraction(c) * p for c, p in c_dist.items())
        mean_c2 = sum(Fraction(c * c) * p for c, p in c_dist.items())
        var_c = mean_c2 - mean_c * mean_c
        results[n] = {
            "num_shapes": len(states),
            "mean_c": mean_c,
            "var_c": var_c,
            "p_star": c_dist[1],
            "p_line": q_line,
            "c_dist": dict(sorted(c_dist.items())),
        }
        if n == max_n:
            break
        next_states = {}
        for key, (adj, prob) in states.items():
            tprobs = transition_probs(adj, model)
            for w, tp in tprobs.items():
                # Attach new vertex n to w
                new_adj = [list(nbrs) for nbrs in adj]
                new_adj[w].append(n)
                new_adj.append([w])
                new_adj_t = tuple(tuple(x) for x in new_adj)
                nkey = canonical_degree_state(new_adj_t)
                if nkey in next_states:
                    old_adj, old_p = next_states[nkey]
                    next_states[nkey] = (old_adj, old_p + prob * tp)
                else:
                    next_states[nkey] = (new_adj_t, prob * tp)
        states = next_states
    return results


def verify_roots():
    def F(mu):
        val, _ = integrate.quad(lambda x: (x ** mu) * math.exp(x - 1.0), 0.0, 1.0, epsabs=1e-14, epsrel=1e-14)
        return val

    mu_dan = optimize.brentq(lambda m: 1.0 - m * m - F(m), 0.1, 0.99, xtol=1e-15)
    mu_pan = optimize.brentq(lambda m: 1.0 - m - F(m), 0.1, 0.99, xtol=1e-15)
    print(f"Analytical roots:")
    print(f"  mu_DAN = {mu_dan:.12f} (paper: 0.771192)")
    print(f"  mu_PAN = {mu_pan:.12f} (paper: 0.549735)")
    print()


def main():
    verify_roots()
    for model in ["IR", "DAN", "PAN"]:
        res = enumerate_model(model, max_n=8)
        print(f"=== Exact enumeration for {model} (starting from L_3) ===")
        for n, d in res.items():
            mc = d["mean_c"]
            vc = d["var_c"]
            print(
                f"  N={n} (shapes={d['num_shapes']:2d}): "
                f"E[C]={mc} = {float(mc):.12f}, "
                f"Var[C]={float(vc):.12f}, "
                f"S_N={d['p_star']}, Q_N={d['p_line']}"
            )
        # Check known exact formulas
        for n in range(3, 9):
            assert res[n]["p_star"] == Fraction(2, n - 1), f"S_N mismatch for {model} N={n}"
        if model == "IR":
            assert res[4]["mean_c"] == Fraction(4, 3)
            assert res[5]["mean_c"] == Fraction(19, 12)
            assert res[6]["mean_c"] == Fraction(647, 360)
            print("  [OK] Matches prompt exact IR values E[C_4]=4/3, E[C_5]=19/12, E[C_6]=647/360")
        elif model == "DAN":
            for n in range(5, 9):
                # Q_5=1/6, Q_{k+1} = (2/k) Q_k
                expected_q = Fraction(1, 6)
                for k in range(5, n):
                    expected_q *= Fraction(2, k)
                assert res[n]["p_line"] == expected_q, f"DAN Q_{n}: {res[n]['p_line']} != {expected_q}"
            print("  [OK] Matches exact DAN S_N=2/(N-1) and Q_N=Eq.(A13) for N=3..8")
        elif model == "PAN":
            for n in range(4, 9):
                expected_q = Fraction(1, 3 * (2 ** (n - 4)))
                assert res[n]["p_line"] == expected_q, f"PAN Q_{n}: {res[n]['p_line']} != {expected_q}"
            print("  [OK] Matches exact PAN S_N=2/(N-1) and Q_N=Eq.(A11) for N=3..8")
        print()


if __name__ == "__main__":
    main()
