> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Paired Neon production slack occupation

**The paid inner rows are almost proportional in occupation, at different scales for `f` and `g`.** For the two pinned D18 Neon witnesses, on nonfounder inner `f` rows the uncancelled share of the weighted sum is at most **0.015162%** at `t=1.2415`; on nonfounder inner–inner `g` edges it is at most **0.098775%** at `t=1`. The finite `f` tail does not share this uniform proportionality. These are occupation inequalities for the named summands; they are not new full-drift or exponent bounds.

Write `a=δ−` for the signed slack of the lower witness, `b=δ+` for the signed slack of the upper witness, and define, on a given paid class,

```
H− = π[(a−t b)+],   H+ = π[(t b−a)+],
S = π[a+t b],      O = π[min(a,t b)],
R = (H−+H+)/S = 1−2 O/S.
```

All these are per core, with **edge multiplicity** for `g`. At `t=1`, `S` is literally `π[δ−+δ+]`; at other `t` it is the sum in the stated pricing units. For any real opposite-sign pair of prices `(c,−d)` with `c>0,d≥0`, use `t=d/c` and scale the result by `c`; reverse signs by swapping which hinge is used. The rational CLI accepts a common `--t`; without it the trial values differ by family. Only use distinct multipliers if the pricing itself permits distinct prices per family.

## Pinned rows and observed local structure

The replay reconstructs both slacks from the [portable raw pinned operands](../rows/README.md), not from published minima. It reads both witnesses with their own signs, `λ` and recharges. At a nonfounder `f` state `(D,k)`, the production inner slack is constant as a function of `ρ`; the finite-tail slack is affine in **the same** `ρ∈[0,(k+1)/2]` for both witnesses. `g` uses the actual opposite endpoint as an included donor (`1/D_child` or `1/D_parent`). Its external donations are centered in the pinned inner–inner charge and passed to higher rows; there is no extra `ρ` in this charge. The tracked pair occupation below still includes the real external rates.

* **Mass-carrying `f`:** newborn, first gain, and second gain ratios are `1.2412775`, `1.2412014`, `1.2415874`. The bare `(D,k)=(2,1)` ratio is `1.2404711`. Pearson correlation across the eight sampled young rows is `0.99999983` and across all 152 inner rows `0.99992988` (unweighted grid summaries, not estimates under `π`). The inner minimum `1.1822093` is at the fully bare `(17,16)`; the maximum `1.2420935` is at `(14,0)`. Stable occupation puts much less weight on the bare extreme than a row envelope would.
* **Mass-carrying `g`:** the edges from the nonfounder `(L,k,h)=(0,1,1)`, `(1,1,1)`, `(1,2,1)` to a virgin have ratios `0.99910786`, `1.0002876`, `1.0003383`. Across 20,672 nonfounder inner pairs the ratio spans `0.9812146` to `1.0182753`. The unweighted Pearson correlation is only `0.00745`: the two deviations around the common ~15 ppb slack are not themselves proportional. The *paid levels*, rather than their centered deviations, are what cancel.
* **Shared-rate tail extreme:** at `(22,4)` the ratios at the two common corners are `11.33133` and `0.098139883`. At `t=1.2415`, the paired weighted sum has minimum **39.505938 ppm**; summing individually minimized terms gives only **7.8026445 ppm**. This loses about 80% of the true pointwise floor in that row. The choice of least favorable `ρ` is different for the two slacks. Suppressing the coupling also raises the occupation upper bounds on the two hinges over all finite `f` rows from `0.675246, 0.836032` ppm to `0.693570, 0.946423` ppm, respectively. The high bare rows have even more extreme pointwise ratios, so it would be unsafe to extrapolate an inner multiplicative envelope to the tails.

## Certified occupation results

The hinges below are in **ppb per core**, rounded outwards. Shares are dimensionless (multiply by 100 for percent).

| Paid class | `t` | `H−` ppb | `H+` ppb | Uncancelled share `R` with signed bisection |
|---|---:|---:|---:|---:|
| `f`, `D≤5,k≤1` | 1.2415 | `[.030146,.042061]` | `[.308534,.456855]` | `[.00011375,.00014909]` |
| `f`, `D≤17` | 1.2415 | `[.150636,.233595]` | `[.363728,.511249]` | `[.00011986,.00015162]` |
| `f`, `18≤D≤64,k≤4` | 1.2415 | `[.005960,314.951410]` | `[0,407.939908]` | `[.00000484,.35696247]` |
| `f`, `D≤64` | 1.2415 | `[.157051,675.245493]` | `[.363844,836.031096]` | `[.00007988,.12675999]` |
| `g`, inner parent and virgin child | 1 | `[.000947,.001777]` | `[.001121,.001841]` | `[.00049586,.00066362]` |
| `g`, both degrees `≤7` | 1 | `[.001951,.002990]` | `[.001773,.002720]` | `[.00050152,.00059680]` |
| `g`, both inner | 1 | `[.004974,.007506]` | `[.004614,.007170]` | `[.00079511,.00098775]` |

For comparison, with a **single** multiplier `t=1`, the inner `f` share is `[.10761641,.10773415]`; with a single `t=1.2415`, the inner `g` share is `[.10753055,.10785918]`. Thus a common calibration cannot be inferred from the two separate sharp results. To price the same multiplier across disjoint `f` and `g` charges, add the absolute hinge bounds per class before forming any joint ratio. Do not add nested classes or reuse these slacks if already spent. The first-body tail beyond `D=64`, mixed/tail edges, and all higher motifs are outside this calculation.

## Proof of the inexpensive bounds

**One core.** Use the [physical stable newborn formula](../renewal/README.md#stronger-integrate-the-physical-root-at-every-age), `π[r]=μ E_β∫exp(−μ v)r_v dv`, with the age of the tracked newborn and its full exterior. The tracked core has leaf gain at `L+ρ` and promotion at `L/D`. In decreasing order of `(D,k)`, with value zero beyond the maximum supported degree, upper/lower comparison potentials solve, with directed rounding, the maximum/minimum of

```
(r(ρ) + (L+ρ) V(D+1,k) + (L/D) V(D,k+1)) / (s+L+ρ+L/D).
```

On each affine segment of the reward this ratio is fractional linear in the *same* `ρ`. For a hinge upper, endpoints of the whole rate interval suffice by convexity after cross multiplication; for a lower include any interior zero of the contrast. For overlap, include the same kink and treat the minimum of two affines as concave (endpoints suffice for its iterative **minimum**, but all kinks are needed for its **maximum**). Use `LO≤μ≤HI` as pinned in `renewal/check.py`; for nonnegative rewards multiply the lower integral at `HI` by `LO`, and the upper integral at `LO` by `HI`. These are Dynkin inequalities against adaptive exterior input, without freezing the physical rates. Lower values are built from nonnegative rewards and nonnegative successor values.

**One edge.** From the same stable newborn, track each promotion at age `v`. It creates a child with `(D,k)=(2,0)`, while the parent's `k` increases by one and its `D` is unchanged. For each such child count its future edge reward without resetting the parent's age or cloning its state. Tonelli for nonnegative rewards and promotion compensation give the outer tracked-parent toll `(L/D) V_pair(post-parent, newborn)` at age `v`. The pair integral is at discount `μ` as well; there is only **one** leading factor `μ` (from the parent's age mixture). In a pair with parent `(D,k)` and child `(E,j)`, the gain rates are `L+1/E+ρ_p` and `l+1/D+ρ_w`, with `0≤ρ_p≤k/2`, `0≤ρ_w≤j/2`; promotions have rates `L/D` and `l/E`. Monotonicity makes degree exit permanent, and both age and childcount increase; a triangular Bellman sweep uses the four rate corners for an exact extremum over the containing rectangle. Physical conditional inputs may be correlated; this rectangle is only an outer control set. The `--prefix` option also computes a nonnegative lower sum using the positive [fringe residence prefix](../renewal/README.md#what-the-19-mark-rank-one-formulas-give); it takes the maximum of the two lower bounds, never their sum. Founder edges are excluded.

**Ratio.** First combine the nonnegative interval bounds using `H−+H++2O=S`. For a tighter bound test `π[|a−tb|−r(a+tb)]` directly. This signed test is convex piecewise affine in the common rate in an `f` row; the upper still uses endpoints and the lower checks kinks. Choose both discount endpoints locally in each signed Bellman step (the ratio in the discount is monotone for a fixed numerator, even if negative). Do not use positive-reward discount monotonicity on the signed test. The factor `μ` is positive and unnecessary for the sign, and a signed pair value can be fed through the outer parent comparator by pointwise order with `L/D≥0`. Directed rounding is at `10^−21`; the optional eight bisection steps per side only tighten endpoints proved by comparator signs.

## Replay

```sh
taskset -c 13 nice -n 19 timeout 50s python3 -Bu round2/occupation-neon/ratio-rows/check.py --tight-ratio
taskset -c 13 nice -n 19 timeout 45s python3 -Bu round2/occupation-neon/ratio-rows/check.py --t 1 --only f-inner g-inner --prefix --tight-ratio
taskset -c 13 nice -n 19 timeout 40s python3 -Bu round2/occupation-neon/ratio-rows/check.py --t 1.2415 --only f-inner g-inner --tight-ratio
```
