> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Signed stationary pressure on a distinguished adjacent donor

**Result:** the matched finite-tree census finds cheap adjacent donors worth several `10⁻⁶` per core in each witness's *correctly signed* positive components. This is a material candidate, not a stationary estimate. At a positive coefficient `.003–.004`, an endpoint move of `10⁻⁶` needs roughly `.00034–.00046` of certified weighted donor incidence/core; a `10⁻⁵` move needs ten times this. Knowing only *total* incoming gain is insufficient: even plugging the exact sampled total back in leaves widths `4.1e−6–1.13e−5` per core across the target families. Keeping the actual adjacent donors of degree at most four as well shrinks those diagnostic widths to `2.7e−8–7.9e−7`. No new stationary increment is booked.

## Owner, sign, and the exact flux relation

Fix `M ∈ {hV,P4,Yt}`, occurrence census `x`, and role `i`. Its charged row is `a(x)+Σ_i(c_i^p r_i^p+c_i^c r_i^c)` in the normalized signed Neon units. The `r` sums are **external** reciprocal sends. Start with the exact charged row (not raw `Δ_L Φ`): for a routed receive, subtract the transferred derivative from the gain coefficient before applying the sign. Expose only one directed child edge `v→i`; its weight is `c_i^c(x)/D_v`. The parent version is identical except that re-rooting the count can change the highest owner. The augmented connected object has one more vertex, distinguished receiver and donor, and retains the old occurrence's symmetry convention; no extra term from promoting at the donor is implied. For any common support test, summing over old occurrences and then donors equals summing over donors and their incident compatible old occurrences. The replay verifies this reversal and reconstructs every retained original charged row exactly on the sampled trees.

Let `π_M` be the stationary incidence measure; `I_M` is the unique post-promotion birth source (see [lineage](../../lineage/README.md)). On a fixed incidence the tracked generator `A_M` changes the degree of role `i` with rate `g_i=g_i^0+r_i^p+r_i^c`, where `g_i^0=L_i+Σ_inside 1/D_j`, or promotes it at rate `L_i/D_i`. For a bounded projected test `H`, exact discounted Dynkin gives

```
π_M(A_M H) = μ π_M(H) − I_M(H).
```

In particular take `H=1[D_i>q] W(x_−i)`. Since own promotion preserves `D_i`, and `W` does not depend on that role's state, exactly

```
π_M[ W 1[D_i=q] (r_i^p+r_i^c) ]
  = μ π_M(H) − I_M(H)
    − π_M[W 1[D_i=q] g_i^0] − π_M[A_transverse H].
```

`A_transverse` carries every other gain and promotion. If `W` also uses `k_i` or a projected support filter, its own promotion belongs to that seam too. Apply this identity on the *geometric* lineage, without killing it when the Neon coefficient or chosen support turns off. A compact aperture needs all its side seams; omitting the birth or seam terms is invalid. To isolate child pressure at the top slot also subtract the actual parent-edge weighted mass. At a non-top slot the parent is already inside `M`, so that subtraction is unnecessary. The stationary relation is exact, and eliminating its left-hand variable in one common linear program preserves all signed cancellations. It does not give a number without the right-hand occupations and the joins of endpoint states across apertures. Replacing them by separately optimized rates discards the useful constraint.

For a cell with incidence mass `m`, if `n_i` outside arrows in a role are allowed, the integrated pressure there is at most `n_i m/2`. Retain cheap donor **counts and reciprocal mass** on each directed role; the missing mass has bound `(n_i m−n_cheap)/5` for child donors if all `D≤4` were accounted for. Counts and both directions of the same edge must use a single joint edge measure. Given one total-pressure equation in a cell, subtract the cheap contributions before distributing the remainder into these shared boxes. Sorting the coefficients and filling the lowest (highest) first gives the exact lower (upper) endpoint of that relaxation. This works for signed coefficients and cannot spend the same total independently on positive and negative corners. For an upper stationarity constraint use the same columns and reverse the interval endpoint, not the edge orientation. A reachable local swap with external parent/child degrees `(2,4)` vs `(4,2)` has the same combined receive `3/4` but changes a `.001` child-only price by `1/4000` per incidence. This is an identifiability obstruction for the projected total, not a claimed stationary error.

For one-sided lower credit it is safer to subtract the independent GMP-style directional floors first. Each positive child slope leaves `c/D_v≥0`; each negative one leaves `|c|(1/2−1/D_v)≥0`. Thus one may discard unresolved arrows or keep a provably supported subset of these **excesses**. Do not book positive signed terms from the uncentered sum and discard negative ones. Cheap donors help the positive slope most; for a negative slope, degree two makes the excess zero. Shared-parent minimization, previously priced full rows, and the floor itself are not extra donor credit.

## Matched finite census and required mass

Four size-16,000 trees at seed `122419`, 546 cores. For common geometrical support use nonfounder top owners, `size>12`, `cores≥5`, and no virgin in the old `M` itself. The **adjoined** donor may be virgin. Bins: `13–17,k≥3`; `18–64,k=3,4`; `18–64,k≥5`; `65+,k≥5`; and the complement. Compare by max with any full-row price on these same occurrences (including older, old-flow, and high-tail), unless a disjoint reward split is proved. `lineage` supplies the exact identity and has no booked increment. `high-parent` pays `g` on hub-to-virgin edges and `hC` on pairs of virgins; neither is this objective. In particular do not rebook the newborn `hV` value assigned to the promoting old `g` row. This replay transposes only the live higher-row affine term of an **already present** occurrence. A separate `g` or one-edge lift is not an additional price for the same higher-row excess unless owners and masks are reconciled.

| Positive component (fixed descriptor) | coefficient `c` | sampled price/core | weighted mass for inward `10⁻⁶` | weighted mass for inward `10⁻⁵` |
|---|---:|---:|---:|---:|
| Upper `Yt`, high apex, tip child send | `.00380353` | `3.4831e−6` | `.00035945` | `.00359447` |
| Upper `P4`, apex `D=81,k=2`, first middle child send | `.00319338` | `2.9243e−6` | `.00042813` | `.00428127` |
| Upper `hV`, apex `D=14,k=4`, apex child sends | `.00102264` | `2.3412e−6` | `.00133690` | `.01336902` |
| Lower `P4`, high apex, terminal child send | `.00141521` | `1.2960e−6` | `.00095875` | `.00958748` |
| Lower `P4`, apex `D=21,k=2`, terminal child send | `.00207079` | `1.2642e−6` | `.00065523` | `.00655226` |

Thresholds use `(desired endpoint move × supplied ceiling)/c`; they are requirements for that fixed coefficient *as the sole extra component*. If certifying only degree-two senders, multiply each mass target by two to get required incidences/core; for the whole `D≤4` class, four times the target suffices without resolving degrees. The strongest sampled upper `Yt` donor was degree two and supplied the entire positive tip component. It has a **negative** coefficient on the lower witness at that same descriptor. Upper `P4`'s strongest first-middle contribution was also degree two and changes sign on the lower witness. Upper `hV` in the `13–17` host bin had three qualifying cheap-strong edges and `2.3412e−6` from those; the lower witness instead had no cheap-strong positive `hV` component in that bin. Across the two `P4` regions the matching lower cheap-strong source contributes about `3.527e−6` in this sample. These are targeted joint masses, not an unconditional abundance of all degree-two cores.

For `hV,P4,Yt`, the exact sums of signed **child** pressure in the sample (positive contribution minus negative magnitude) are respectively `+8.442e−6,−1.020e−6,+3.588e−6` upper and `−7.171e−6,+2.979e−6,−3.995e−6` lower. These omit the affine constants and parent terms and need not be nonnegative. The oracle diagnostic below bounds the full external affine pressure by using progressively more *true sample inputs*; it is not a certified stationary enclosure:

| Side / family | Independent corner width/core | Exact total oracle | Exact total + cheap edge oracle |
|---|---:|---:|---:|
| Upper `hV` | `2.233e−5` | `8.168e−6` | `7.909e−7` |
| Upper `P4` | `1.219e−5` | `1.133e−5` | `3.694e−7` |
| Upper `Yt` | `4.522e−6` | `4.499e−6` | `4.200e−8` |
| Lower `hV` | `2.035e−5` | `7.563e−6` | `7.051e−7` |
| Lower `P4` | `1.164e−5` | `1.059e−5` | `4.273e−7` |
| Lower `Yt` | `4.142e−6` | `4.127e−6` | `2.739e−8` |

A child-rate-weighted draw is useful **after restricting to the joined old-motif aperture**. The sampler enumerates that exact finite-tree sum; equivalently draw a donor from its reciprocal-weighted slot measure and evaluate the signed joined coefficient. To turn this into a certified endpoint change, certify either the low-degree crossing/occupation of the enlarged motif from its true birth source, or the aperture residual with the joint constraints above. Sampling more trees does not supply such a guarantee. For signed remainder tolerance `ε`, a missing weighted mass of unknown sign with coefficient bounded by `C` costs at most `C·mass`; e.g. `ε=10⁻⁸` at `C=.004` calls for at most `2.5×10⁻⁶` unresolved reciprocal incidence/core unless cancellations are certified. No claim is made that the common current stationary occupations attain the sample ranks.

## Replay

Using the original hash-pinned mmap payloads and existing Astra helper layout, run for each `BLOB`:

```sh
taskset -c "$CPU" nice -n 17 timeout 30s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  "$PYTHON" -Bu "$SCRIPT" "$BLOB" --verify-hash
```

`SCRIPT` is this directory's `replay.py`. Both commands exited 0 (`11.90 s` upper, `5.41 s` lower; under `144,000 KiB` RSS). They verify the payload hashes, three exact projected crossing/seam controls, two signed joint-box controls, and exact receiver/donor and original-row equalities on all four trees per blob. This does not re-run the full GMP verifier; those pinned obligations are documented in [the higher-order replay](../README.md#production-obligations-and-provenance).
