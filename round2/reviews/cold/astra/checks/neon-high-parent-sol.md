# Neon D18 cold check: size-20 plus high-parent

**CORRECT in scope.** The conservative outward candidate is **`[0.565676634491, 0.565771874464]`** on the two original D18 witnesses, conditional on the accepted full B GMP verifications, the stationary mean and physical stable-age identities, and the incoming physical bracket `[0.565673570322, 0.565775]`. This is distinct from the already-reviewed Claim 21 release (`occupation-neon/release/README.md`), `[0.565675793084, 0.565772420378]`.

## Charge ledger and overlap

The accounting key is the **current charged row and its physical embedding**, owned at its top core. A transition can occur in several terms of the decomposition without making their final charged rows the same. Conversely, changing the occupation argument, the observer, or the time of entry does not create a new charged row.

| Booked term | Present-time support and disjointness |
|---|---|
| All-time first body | Only `f` of nonfounders; the newborn vertex and edge potential terms assigned to `f` at promotion stay there. The remainder transferred to `g` is not counted again as `f`. |
| Size-12 higher | Top-owned `hV,hC,P4,Yt,Yb,S3` on whole planted fringes of size at most 12. The exterior parent is excluded; the root and all descendants and leaves are included. A planted degree-18 root alone forces at least 18 vertices, so no high-parent pair is in this class. |
| Cumulative size-20 higher | Same six families, size 13–20, three or four cores, current degrees at most 16 so every gain successor stays inner. Disjoint from size-12 by size and from high-parent `hC` both by the seam and by core count: the latter has at least four children plus the owner. |
| Actual inner-edge floor in the size-20 solve | The minimum is applied only to actual root→child `g` with both degrees below 18 in the retained solve (size at most 20, at most four cores). This part can pay before size 13; it has no competing `g` in the size-12 or `f` ledger. In particular, it does not use the all-edge shortcut that reimburses missing edges out of `f` surplus. |
| High-parent joint payment | `g` with a degree-18–191 parent and virgin child, plus `hC` over unordered pairs of virgins at the same degree and arity 4–10. The edge arity is 1–10. A virgin is exactly `(L,k,h)=(1,0,1)`: a bare degree-two child is not virgin. High-parent `g` is mixed, so it is outside the inner-edge floor even on small fringes. `g` and `hC` are different rows even when a paid pair contains two paid edges. |

The lower-price optimizations do not split these rows further: add only the joint high-parent result; degree anchors 18, 32, 64 and 128 are slices of its support. The size-20 joint physical reward is the result reproduced here. Residual renewal and physical birth-Palm weighting price the same reward, so one can take the better valid price, never add them. Likewise the young-state deficits in the all-degree **potential** ceiling are not signed numerator credits. Allocation of potential to a bottom or other core for the ceiling does not reassign ownership of a charged drift row.

Nearby competing ledgers were checked for the following seams:

- The arms pilot pays only higher rows on size greater than 12 and at most four cores. Its mixed/tail higher rows are disjoint from high-parent `g` and its at-least-five-core `hC`, but the all-inner suffix overlaps the size-20 higher charge. The arms aggregate cannot be added.
- The young-arm fringe pilot pays selected, observed children in states `(2,0)` and `(3,0)` in `(D,k)` coordinates. Its observed virgin edges and its all-virgin pairs where root arity is at least four overlap high-parent `g,hC` in the common degree and arity ranges. The historical observation gate does not prevent that overlap. Pairs with a degree-three child, mixed edges to degree-three children, and `S3` are different from the selected high-parent rows, but these components are not a separately booked price. In particular, young-arm `S3` at the high host cannot meet the size-20 higher degree guard. The companion large-crown fringe allowance excludes `g,hC,S3`; that alone does not separate its other higher rows from the size-20 suffix. Its parent-split bookkeeping reconstructs the same charged row rather than inventing a new charge at the payer.
- The diagnostic mixed-edge pilot in first virgin residence overlaps high-parent `g`. The follow-on residence prices a large-parent `g` edge only after the child either gains to `(3,0)` or buds to bare `(2,1)`. Both are outside the candidate: the parent is tail, the child is no longer virgin, and the family is still `g`. The degree-three branch can overlap the separately observed degree-three young-arm objective. No diagnostic price is booked here.
- The ratio continuation prices a selected `P4` with at least five physical cores and at least 13 vertices at its owner, so it does not intersect either higher size ledger or the new `g,hC`. Its separate witness contrast and price are not inputs to this bracket.

## Physical and numerical checks

For size 20, the suffix replay (`occupation-neon/renewal/next/replay.py`) checks the finite crown subsolution with the same exterior send `r` in the signed charge and in the root's gain successor. The intrinsic rate sum is `n−1/D`: the only missing send from the crown is the root's outgoing contact to its parent. The incoming parent send adds `r`; the actual internal donors are retained. Both cross-multiplied endpoint inequalities cover `r∈[0,1/i]` for each permanent birth-mark group, since the live parent's degree cannot decrease. Crossing the size or core cap is absorbing for this payoff. The all-inner higher guard checks the complete fringe, whereas `g` only checks its actual two endpoints.

The grouped newborn matrices count direct promotions by source and target mark. The physical law is birth-Palm `β`, with no killing; the residual law uses normalized left Perron `ell` and a killed original particle. Their discount is the upper physical prior. For each law the initial column minima are lower bounds, and subsequent simultaneous iterates put unassigned probability on the column minimum. The payoff puts that remaining probability on the smallest admissible value, not on a conditional or renormalized distribution. In the optional dual, `v≥0` makes `βA≤β` imply `α≤βf` after multiplying `α+v=f+Av` by `β`. The completed cap-20 run uses physical `β` and multiplies the discounted newborn value by the lower physical prior exactly once. I inspected this construction and replayed its exact suffix fractions; I did not independently reimplement the full cap-20 geometry table.

For the high-parent charge I compared the specialized formulas to the separate general physical charge extractor (`occupation-neon/fringe/charges.py`) from the original blobs at all 1,740 `(D,k)` settings per witness. This checks the pinned GMP `g` source terms at lines 1218–1255 and tail `g` transfer at lines 1280–1294 against the exact physical Bernstein contraction: `f→g` uses the physical first-body transfer, while the parent’s outgoing `g` transfer uses `3t²(1−t)(F₀−F₁)/18+t³Yref`, with `t=18/D`. The centered transfer is retained and the external donor sum multiplies the untransferred derivative. The `hC` pair takes both incident `g` transfers. On this tail-hub/virgin-pair face there is no Level-4 outgoing transfer or booked Level-4 newborn. The free donor counts are `k` for the edge and `k−1` for the pair; their signed affine lower envelopes agree exactly. The author replay also passed its concrete whole-star drift checks.

The companion fixture (`neon_high_parent_sol.py`) enumerates directed physical arrows on a labelled asymmetric degree-18 crown, with a real degree-seven exterior parent and virgin, bare and degree-three arms. Independently rebuilding the projected state after each arrow recovers gain, promotion and virgin loss. In general the projected gain rate is between `L+j/2` and `L+(k+1)/2`; promotion is `L/D`; virgin loss is `j(3/2+1/D)`. The last term includes the arm's leaf send, own promotion and parent send. Arrows leaving virgins toward the hub are gains of the hub, not losses of virgin status. The fixture solves a separate integer-grid backward recurrence, tests both exact cross-multiplied gain extremes, and combines `j` edges with `binom(j,2)` unordered pairs. The joint value exceeds the sum of separately minimized values on both sides. Degree and arity only increase, so the zero cutoff continuation is valid. Discounted Dynkin applies even with adapted untracked arms, and the physical stable-age identity supplies one factor of `LO`. No product law for arms or their parent is used.

## Exact composition

Use signs `(+1,−1)` on `ℒΦ−λΦ`. The all-time `f` replay (`occupation-neon/renewal/README.md#the-founder-universal-c0-and-the-average-denominator`) and size-12 replay (`occupation-neon/higher/README.md`) supply the first two conservative decimal floors. The cap-20 replay supplies the exact suffix fractions; the fixture computes high-parent jointly and does not use the high-parent `PAID` input for the cumulative amount.

| Term | Lower | Upper |
|---|---:|---:|
| Original trial | `.565673` | `.565775` |
| Charged `f` floor | `2.512846393e−6` | `2.025752646e−6` |
| Size-12 higher floor | `1.276902849680e−6` | `1.501029458674e−6` |
| Size-20 suffix | `58961819710386489/10^23` | `46635843154908967/10^23` |
| Cumulative before high-parent | `437936743978386489/10^23` | `399314053622308967/10^23` |
| High-parent joint | `55203067615510485475692014463/10^35` | `139999781774155931356847859909/(5·10^35)` |
| Fixed all-degree ceiling | `1.356833021311` | `1.367169995210` |

The mean signed inequality is `sign·(μ−λ)·barΦ ≥ paid`. Divide the summed reward by the one ceiling for that side, then add or subtract it at the original trial. This gives the exact unrounded endpoints

```text
76752873706217324196974475692014463/135683302131100000000000000000000000
386753165449918989299233643152140091/683584997605000000000000000000000000
```

Floor the first and ceil the second at `10^−12`. Neither the outward size-20 endpoint nor the already-reviewed release endpoint is substituted for a trial. The independently computed fractions are for these conservative input floors, not replacements for the release's separate exact fractions.

## Evidence

Replay from the repository root, with `PYTHON`, `D18_LOWER` and `D18_UPPER` set as in the release instructions (`occupation-neon/high-parent/README.md#replay`):

```sh
sha256sum -c round2/occupation-neon/release/inputs.txt
sha256sum "$D18_LOWER" "$D18_UPPER"
OPENBLAS_NUM_THREADS=1 "$PYTHON" -Bu round2/occupation-neon/renewal/next/replay.py --cap 20 "$D18_LOWER" "$D18_UPPER"
OPENBLAS_NUM_THREADS=1 "$PYTHON" -Bu round2/occupation-neon/high-parent/replay.py "$D18_LOWER" "$D18_UPPER"
OPENBLAS_NUM_THREADS=1 "$PYTHON" -Bu round2/checks/neon_high_parent_sol.py "$D18_LOWER" "$D18_UPPER"
```

All ten release manifest entries returned `OK`. Payload hashes were `9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec` and `1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41`. The cap-20 replay printed `0.565676227639` and `0.565772079266` with the exact suffix and cumulative fractions above. Both the high-parent replay and the independent fixture printed `0.565676634491` and `0.565771874464`; the fixture reported `1740 exact g/hC comparisons` on each side and passed the physical census and charge-support checks. Full GMP B and the all-degree ceiling scans were inherited, not rerun here.
