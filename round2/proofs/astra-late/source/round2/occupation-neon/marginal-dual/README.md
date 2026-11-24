> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Stationary marginal dual for the D18 pair

**The radius-one stationary LP detects a very small positive paid contrast, but does not separate the full ratio or certify eight digits.** At `s=.56572365`, `t=1.24`, and degree cap eight, the signed combination of the paid `f` and `g` rows lies in `[9.067898e−10, 1.0004244e−9]` per core. Its two repaired duals lose at most `1.07e−15` to rounding. The unpaid charge, including all higher motifs and degree overflow, is not comparably small: a valid outer interval for the **whole** contrast is `[-1.227103894847e−4, 9.1283724734e−5]`. A sign on the paid contrast cannot be promoted to a full ratio sign. This is a static marginal calculation, with no old-flow ancestor-age pricing.

## Target and the sign of a separator

For the two production potentials `P=Φ−`, `Q=Φ+` (scores per core when averaged), write `a=.565673`, `b=.565775`, `δ−=LP−aP`, `δ+=bQ−LQ`. The whole positive drift scores satisfy

```
πδ−=(μ−a)πP,   πδ+=(b−μ)πQ,
D_t=δ−−tδ+ = L(P+tQ) − (aP+t bQ),
πD_t / π(P+tQ) = μ − [aπP+t bπQ]/π(P+tQ).
```

The generator identity holds for bounded per-core score allocations by the normalized finite-tree drift and common first-moment limit, as in the [production occupation argument](../renewal/README.md). For a proof at a fixed trial `s`, the test `(L−s)H` has zero mean for any included score `H`. In particular, using `H=P+tQ` and taking `D_t−[(s−a)P−t(b−s)Q]` as the objective just reproduces an equality. It cannot separate a feasible trial. Instead add local stationary rows that constrain how a pattern distribution can simultaneously realize **the actual full contrast and the two potentials**, or use those rows to refute trial feasibility directly. An uncertified LP infeasibility is not an exponent sign.

There is a calibration issue before any eight-digit solve. From the existing [potential comparison](../potential-ratio/README.md), `.999≤πQ/πP≤1.12`. **If** `.56572365` is critical, the whole residual ratio must be in `[.880685769926, .987355417736]` (outward). A whole ratio equal to `1.24` would instead imply `μ∈[.565729439064,.565732300737]`; even a perfectly known `1.24` leaves about `2.86e−6` with that potential interval. The earlier inner `f` proportionality and the `1.24` contrast used here apply to selected charges, not the sum over the production witnesses. A near-critical search should work on the *product* `r=(πδ−/πδ+)(πQ/πP)` or an equivalent direct separator. Since `μ=(a+b r)/(1+r)`, at `.56572365` the local derivative is about `2.59e−5` per unit of `r`; a width `10⁻⁸` needs an odds width about `3.9e−4`. Equivalently, near a centered `t`, the whole signed defect should be controlled on the per-core scale `10⁻⁸·π(P+tQ)`; at `t=1.24` the published `C0` ceilings make this at most `4.27e−8` even for a full `10⁻⁸` radius, or about `2.14e−8` for a half-width. The unpriced term in this pilot is orders of magnitude above either scale.

## Model and exact reconstruction

At trial `s`, introduce nonnegative measures for every degree/child-count type `u=(D,k)` through `D≤N` and an explicit overflow `*`. Radius zero has vertex mass `x_u`; radius one has the directed parent–child mass `y_uv`, counted once per edge / per core. Permanent birth marks and age are not needed by this instantaneous generator. The orientation distinguishes parent and child, so asserting `y_uv=y_vu` would be wrong. Instead re-rooting on the reversed incident dart uses the **same** measure by exchanging its two oriented slots; with these conventions visibility in the two orientations is automatic. The founder has vanishing normalized mass. Impose

* `Σ x=1`, `Σ_parent y_parent,u=x_u` for every type including overflow, `Σ_child y_u,child=k x_u` for each finite type, and the exact total edge normalization implied by these. The unknown overflow child count is represented by its outgoing edge mass, not by a fictitious terminal count.
* Both included-endpoint reciprocals are exact when finite. If a donor is overflow, introduce the joint moment of its reciprocal in that edge cell, between zero and `y/(N+1)`. The same columns are used in vertex and edge balances and in any affine charged reward.
* Balance every finite vertex test and every finite–finite oriented edge test using exact leaf, degree-gain, promotion, and creation rates. Newly promoted vertices enter the virgin state and their oriented edge enters the **post-promotion** parent cell. Birth intensity is exactly `s` in the vertex balances; the overflow contributes between zero and its mass to this intensity. Degrees never decrease, so no source returns from overflow. Retain the age/Yule upper bound on `x_*` by domination of the tracking degree gain by `D`; no time marginal is imposed beyond this valid cut.
* For existing edge gains from the third and subsequent cores use nonnegative joint input fluxes, with each flux capped by the number of other slots times `y/2`. On summing over the `k` children of a given parent, the outside donor input is `k` times total incoming reciprocal minus the included child reciprocals. On summing over the single parent of a child, it equals that child's incoming reciprocal from all its children. These are exact projections of the missing labelled radius-two patterns, not independently chosen gain intervals in each equality. They are weaker than full radius-two projectivity: a flux is a moment, not a complete three-site marginal.

The objective is the original pinned signed `f` contrast for finite core states plus the pinned inner `g` contrast for finite–finite edges. The low-cap replay uses constant inner `f` rows, but the compiler retains common affine `ρ` when finite tail rows are enabled. Both objectives are obtained from [the pinned operands](../rows/README.md) and use neither independently minimized slacks nor newly assigned signs. A direct jump fixture checks every generator equality against reevaluation of vertex and edge counts, including post-birth parent routing and overflow reciprocal terms. Its oriented cycle is only an algebraic device to exercise all vertices with a parent, not a claim that cyclic states occur in a tree or satisfy the stationary balances.

### Exact dual and overflow payment

For an LP `min c·x`, `Ax=b`, `Gx≤h`, `0≤x≤u`, round the proposed row multipliers to rationals and enforce `z≤0` for inequalities. A certified lower bound is

```
y·b + z·h − Σ_j u_j [(Aᵀy+Gᵀz−c)_j]_+ .
```

No solver feasibility tolerance is used in this check; HiGHS only proposes the multipliers. Each variable cap follows from probability, edge count, or physical rate. The dual `l1` norm also bounds amplification of a uniform additive right-hand-side error; for `N=8` its two paid-contrast values are below `1.57e−6`. The largest multiplier is below `8.94e−9` for those particular contrast solves, so a known bound on total row defect can also be charged in `l1` using that multiplier ceiling. These small numbers apply to the paid objectives and the stated row normalization only, not to arbitrary higher-motif duals.

| Cap | Variables / equalities / inequalities | Required trial overflow mass | Paid contrast interval |
|---:|---:|---:|---:|
| 6 | `949 / 397 / 595` | `[.2897378131102621,.3216507853218084]` | `[8.602715e−10,9.646528e−10]` |
| 8 | `2965 / 1150 / 1888` | `[.2439164361343310,.2649166111381036]` | `[9.067898e−10,1.0004244e−9]` |

The trial is conditional: these apply if the physical exponent equals `s`. The overflow lower bound is not a measured physical statistic independent of that assumption. The bound for the omitted score itself uses an all-tail accounting identity: since production charges are nonnegative and the paid terms appear exactly once, the unpaid `δ−` is at most `(s−a) C0− − paid_lower−`; the unpaid `δ+` is at most `(b−s) C0+ − paid_lower+`. At `N=8` these ceilings are `9.12827243096e−5` and `9.8960722802e−5`. This includes overflow and higher charges even when all of their mass is concentrated in the adverse part of the contrast. It is an upper budget, **not** evidence that the actual uncertainty is that wide or that a sharper hierarchy cannot succeed.

## Separator to try next

Use columns for ordered, oriented connected patterns rather than full balls: take each of the eight D18 motif supports, its one-donor extension for every gain, and the preimage obtained by deleting the newly promoted terminal vertex for each birth. Include role-labelled existing-star/path and mixed-degree cells at the D18 transition seams. Impose `y_u=Σ` of projected embeddings with the correct falling factorial when extra same-vertex arms are distinguished. Reverse re-rooting on shared oriented darts and quotient only actual pattern automorphisms, not distinct donor roles. Generate both the `D_t` and `P+tQ` coefficients on the **same** Bernstein controls before introducing tail caps. Grow the test basis using the largest exact missing reduced-cost contributions, especially the older `hV`, `P4`, `Yt`, and `Yb` classes identified in the production diagnostic (`round2/occupation-neon/rows/diagnostic/README.md`, historical reference). The pilot shows that enumerating more low-degree edge cells alone improves an already small term.

For a tail embedding keep joint positive moment columns for the Bernstein controls and for the actual donor reciprocal, with mass projections and valid moment inequalities in each role. A checker must either bound the residual of the full dual objective on each remaining mixed/overflow cell or use a valid total-mass envelope times its **paired** residual; paying independent upper and lower production budgets is too wide at the present scale. To exclude an entire trial interval rather than isolated points, reuse a frozen rational separating dual throughout it. With a fixed variable-cap vector and affine coefficient/RHS dependence on `s`, the repaired lower bound above is concave in `s` (sum of minima of affine functions), so strict positive endpoint values exclude that whole interval. For an infeasibility certificate with no observable objective, use the same check with `c=0`; a positive lower bound is then a contradiction. If a constraint has a nonlinear RHS (the optional Yule cut), freeze a valid interval-wide RHS before that endpoint check. No such separator at eight-digit scale was obtained here.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 taskset -c 13 nice -n 12 timeout 70s \
  ./round2/.venv/bin/python -Bu \
  round2/occupation-neon/marginal-dual/pilot.py --caps 6 8 --full-budget
```

The replay completed for both caps, including direct generator fixtures and rational dual repair. HiGHS via SciPy `1.14.1` is only used to find candidate duals. The right-hand side bounds displayed above are rebuilt in rational arithmetic. No run at cap 18 or for the proposed higher-pattern separator is claimed.
