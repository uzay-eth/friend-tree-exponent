#!/usr/bin/env python3
"""Compare the apex score with the published whole-tree motif decomposition."""

from fractions import Fraction as F
from itertools import combinations

from small_fringe import Apex, NEW, intrinsic, load_reference
from motif_identity import build, decomposition, tree_states


def complete(tree, parent_degree):
    par = build((parent_degree - 1, [tree]))
    # Leaves of the artificial founder precede the one child.
    focal = parent_degree
    degrees, cores, children, states = tree_states(par)
    assert degrees[0] == parent_degree and par[focal] == 0
    return par, cores, children, states, focal


def reference_apex(cert, tree, degree):
    par, cores, children, states, focal = complete(tree, degree)
    apex = list(cores)
    apex.extend(par[w] for w in cores if par[w] >= 0)
    apex.extend(par[par[w]] for w in cores
                if par[w] >= 0 and par[par[w]] >= 0 and states[par[w]] in cert.hu)
    for v in cores:
        if states[v] in cert.cu:
            apex.extend(v for _ in combinations(children[v], 2))
    seen = decomposition(cert, par)
    assert len(apex) == len(seen)
    sign = 1 if cert.mode == 1 else -1
    return sum(sign * (cert.lam * f - drift) for root, (_, drift, f)
               in zip(apex, seen) if root == focal)


def main():
    states = [NEW]
    frontier = {NEW}
    for _ in range(5):
        frontier = {y for x in frontier for y in intrinsic(x)}
        states.extend(sorted(frontier))
    for mode in (1, 2):
        cert = load_reference(mode)
        apex = Apex(cert)
        hits = 0
        mismatches = 0
        for state in states:
            v0, v1 = apex.score(state, F()), apex.score(state, F(1,2))
            for d in (2, 3, 7, 29):
                ref = reference_apex(cert, state, d)
                assert ref == v0 + F(2, d) * (v1 - v0)
                assert ref == apex.score(state, F(1,d))
                hits += 1
                mismatches += ref != v0
        assert mismatches > 0
        print(f"mode {mode}: {hits} complete-tree apex comparisons agree, "
              f"silencing the outside parent changes {mismatches} scores.")


if __name__ == "__main__":
    main()
