> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Staged hV → P4/Yt correction over the frozen edge witness

Matched no-memory exact replay gives outward **`[0.565677997980320, 0.565770504801779]`**. That improves the reviewed edge enclosure by `1.72942e−8` lower and `8.47817e−9` upper. The requested `5e−6` on each side was not reached: the upward-rounded Bellman comparison leaves shortfalls of `6.72920e−6` and `6.77705e−6` in signed margin at the respective `5e−6` targets, even granting a rounded-up ceiling on the entire frozen reserve *as booked*. This is a barrier for the tested finite box/cut hierarchy with its relaxed rates, not for the physical exponent or larger projections.

## Source, support and proof

Frozen edge files (SHA-256; none modified): `README.md` `ff7803d6d39078d14fe47d464be77c677cdcad381e06c4b1caa054cd0ef5cda6`; `replay.py` `e90651922acb6b3a427872cf3801f4b9097d2b83c7abe0ab75add6e60c78350f`; `regression.py` `7a5979e12e1c15b20995d65c228140aba277b750424423062f820ee7564ccb5a`. The frozen edge replay embeds the dependency and original payload hashes, and uses release commit `80fe02660b9f0a0d4e7fdc8240b0815615ff5896`. The hierarchy authenticates the edge replay, its source dependencies, and each original payload before using its certified rows. No pressure-memory is used.

For an extant connected ordered motif, retain every participating vertex's physical `(D,k)`, with unordered tips canonicalized for `Yt`. The real gains use leaf stock, actual `1/D` from each internal contact, and one exterior reciprocal sum at each vertex. The latter has upper bound one half times the number of external incidences. The charge from the production tables is affine in these same sums; coefficients are checked equal across exterior roles at a common receiver on the paid support. Minimize the shared-rate fractional-linear Bellman quotient by sorting these affine slopes. After downward rounding independently check the cross-multiplied inequality by minimizing the affine residual over the entire box. On the upward solve, a checked corner yields an upper comparison on that same finite monotone DAG. A cut violation cannot return. Only the reward, not the state or its continuation, is gated by charge support.

Each genuine new `P4` or `Yt` is assigned to the unique existing `hV` obtained by deleting its newborn terminal tip. The boundary value reads the promoter's **post-promotion** `k` and the virgin newborn; each old `hV` contributes one extension at its bottom and one at its middle, at their own promotion rates. Similarly a new `hV` is assigned to its old edge and a new edge to its old parent census. Every pre-existing incidence still changes normally on the same event. The singleton normalization gives zero newborn singleton value. Founder-owned corrections are excluded in all layers. Domains at a birth are nested into the predecessor layer. These are finite-tree generator inequalities, with stationary normalization only for the frozen reserve and denominator in the final composition.

**Disjoint mask.** Let `r` be the number of participating vertices. The complete planted apex fringe has at least `ΣD−r+1` ordinary vertices: count the tagged cores, their leaves and the distinct immediate off-motif children. It has at least `1+Σk` cores. Pay the three genuine charges only if the first bound exceeds 12 and either `Σk≥4` or a tagged `D≥17`. This excludes the frozen small-fringe block and the interior size-13–20 suffix (`cores≤4`, `D≤16`). High-parent pays `g,hC`, neither of which the higher layers spend. The original `g` reservations and `f` accounting stay exactly those of the frozen edge replay. `--replace` admits all retained descriptors, but drops the entire size-12 and suffix reserves (including their joint edge credit); it retains only high-parent. It is an overlapping alternative, never additive with the first mask.

The absolute per-core norm adds the frozen census/edge norm and the three motif sup norms, with a multiplier `k_middle_cap` for `Yt`: a bottom vertex identifies one path in `hV/P4`, and at most that many other tips in `Yt`. At the improved trials these bounds are `<.000707168` and `<.000207852`. The signed denominator adjustment uses this motif sup only on the lower side, because the added motif values are nonnegative. Both exact generators passed 128 full-tree potential-jump checks, including positive `P4` and `Yt` births; minimum surpluses over the aggregate signed local margin were `5.135075426e−6` and `4.264437811e−6`.

## Density and the target shortfall

The lower broad no-memory run admitted `79,663` `hV`, `164,582` `P4` and `74,389` `Yt` paid contexts, including `33,376`, `51,532` and `23,224` with apex arity ≥5. It evaluated `22,632` live four-body birth sources and `3,528` live `hV` birth sources. These are row counts, not physical mass estimates. At the `5e−6` target, `hV` alone added `2.51127e−8` to the lower margin; adding the two four-body families bought only `3.88866e−10` more. The upper matched run added `1.36214e−8` and a further `6.50402e−12`. Giving the upper solve all retained small descriptors (larger cut, `--replace`) added only `3.76339e−8` to its baseline margin, against loss of about `1.97e−6` of booked reserve. Thus the small result is not explained just by the disjoint gate. The next substantive revision would have to resolve the donor state or pressure at the particular paying receiver, or improve the exterior-rate relaxation; increasing precision on these same rows cannot bridge the shortfall. Any such comparison must keep the reserve mask or replace it explicitly.

## Replay

From the checkout, serialize the runs:

```sh
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1
P=./round2/.venv/bin/python
R=round2/occupation-neon/interpolate/edge/higher/replay.py
nice -n 10 taskset -c 7 "$P" -Bu "$R" b3f18_lo.bin --config 32 7 7 4 7 4 24 6 --full-only --s .565677930574
nice -n 10 taskset -c 7 "$P" -Bu "$R" b3f18_hi.bin --full-only --s .565770634792
nice -n 10 taskset -c 7 "$P" -Bu "$R" b3f18_lo.bin --config 32 7 7 4 7 4 24 6 --robust-upper
nice -n 10 taskset -c 7 "$P" -Bu "$R" b3f18_hi.bin --robust-upper
nice -n 10 taskset -c 7 "$P" -Bu "$R" b3f18_hi.bin --config 32 7 7 4 7 4 24 6 --full-only --replace --robust-upper
```

Completed respectively in 546.06 s, 91.81 s, 421.59 s, 146.06 s, and 699.02 s. The frozen edge was independently audited; this new layer has not yet received a separate cold review.
