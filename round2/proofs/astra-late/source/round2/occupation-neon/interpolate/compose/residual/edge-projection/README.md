> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../README.md#claim-21) for exact scope and citations.

# Remaining D18 first-body and edge slack after the edge correction

This is a **candidate residual dual**, pending cold/clean review. No new
exponent endpoint is composed here. Both original payload hashes are
authenticated, as are the edge producer and its frozen dependencies.

## Charge partition

Use the stronger disjoint edge witnesses at the original two trials
in [`interpolate/edge`](../../../edge/README.md), with its exact Bellman
tables `z,V` and margin `m=s z(J)−δM`. For a nonfounder census `x`, let
`Q_x=q(x)+t(k−1)` be the same conservative allocated potential upper.
For each existing oriented edge `e=(p,y)` whose parent is a nonfounder,
write `G_e` for its original signed charged edge row at the old `λ`, and
`G_res` for the original inner or high-parent reserve subtracted by the
edge producer. The disjoint remainders are

```
R_f = F + δ(M−Q_x) + (A−s)z + b V(post-parent,J),
R_g = G_e − G_res + (A_edge−s)V.
```

The newborn edge term goes in `R_f` **once**; the four transitions for
`R_g` act on the old surviving edge only. Regions outside the edge cut
cannot enter it by a gain or a promotion. On those edges the adjoint
term vanishes. The existing case checks guarantee the remainder is
nonnegative after removal of exactly the producer's two `g` reserves.
The all-degree census tail pays the earlier verified arity-only balance;
no new edge is retained beyond it. The exact `R_f` corner check spans all
2,079 finite census rows per side. On the priced `g` mask the replay
reconstructs its signed affine row and verifies all four corners anew.

With the same founder payment, the signed new drift dominates
`m(C−1) + original reserves + untouched six higher families + ΣR_f + ΣR_g`.
The objective pays `R_f` only through degree 18, and `R_g` only when both
endpoint degrees are at most 4. Those edges are all inner and lose the
entire uniform inner reserve before being booked. No high-parent reward,
higher-body charge, or gap from the conservative `−δΦ` allocation is
booked again. This projection is separate from fresh higher-mode work.

## Shared donor constraints and all-degree dual

Use the canonical pressure-memory implementation (`pressure-memory/marginal/replay.py`, historical reference)
at the repository root (single copy and SHA pin). Keep exact endpoint
census bins through degree 4 so the edge row's coefficients are exact.
If `u,v` are the exterior reciprocal inputs of the parent and child of
an edge, its two **total** receive rates are `u+1/D_child` and
`v+1/D_parent`. The incoming edge columns sum the latter against the
vertex `ρ` census moment, and outgoing columns sum the former against
the `kρ` census moment. The identical `ρ` enters the `f` objective and
the vertex gain in the balance equations. On exact child bins, the
same edge rate also gives `w=(L+v+1/D_parent)/D_child` for the pressure
source; the quadratic receiver hull and the other endpoint marginals
stay shared. This ties first moments; it does not claim to reconstruct
the full conditional distribution of adjacent exterior inputs.

The LP still includes the original tail mass, mass-transport balances,
reciprocal hulls through reciprocal zero, and separate excess-stub and
excess-upward-input budgets for unbounded arity. Outside the paid mask
the coefficient is zero, not an extrapolated signed charge. Rounded
dual multipliers are repaired by exact rational column inequalities for
**all** columns, including zero-objective and saturated tail columns.
That is the global validity of the lower prices, despite the small
mask. With memory switched off only the 20 `B` balance rows are removed;
all 9,674 columns have the same remaining coefficients and objectives.

## Measured lower prices per physical core

| Witness | Constraints | Unused `f` | Unused `g` alone | Joint `f+g` |
|---|---|---:|---:|---:|
| Lower | Shared endpoint flow | `8.97731984e−8` | `0` | `9.14435182e−8` |
| Lower | Flow + pressure memory | `1.016956395e−7` | `0` | `1.030961604e−7` |
| Upper | Shared endpoint flow | `5.10735310e−8` | `0` | `5.13398900e−8` |
| Upper | Flow + pressure memory | `5.98826281e−8` | `0` | `5.99909484e−8` |

The isolated edge objective provides only the nonnegativity floor on
this mask after exact dual repair. A joint objective has a real
additional floor relative to the displayed `f`-only certificate. These
are comparisons of certified lower bounds, not claims about optimizing
over the exact stationary law. The tiny upper-side joint increment
should not be extrapolated to higher degrees or larger motifs.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 15 /usr/bin/time -f 'wall=%e cpu=%U+%S' round2/.work/venv/bin/python -Bu round2/occupation-neon/interpolate/compose/residual/edge-projection/replay.py b3f18_hi.bin b3f18_lo.bin --exact 4
```

[Raw exact replay](projection.log): 199.48 s elapsed, 136.34+2.91 CPU s.
The matched-column assertion added afterward was checked separately:
9,674 identical columns after filtering out 20 memory rows. All paths
resolve relative to the checkout; no second copy of the memory source.
