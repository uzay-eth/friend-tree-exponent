> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Nonvirgin spines at high-degree Neon D18 apexes

An exact first-arm replay pays **`3.824730245e−14` / `1.749310756e−14` per core** (rounded down, lower/upper witnesses). The inward endpoint moves are at least `2.818865833e−14` / `1.279512249e−14`. Added to the [high-parent bracket](../../high-parent/README.md), that gives outward **`[0.5656766344910281, 0.5657718744639873]`**. The gain is about `4.02e−10` of the original `.000102` gap, and changes neither endpoint at 12 decimal places. It is too small to pursue by increasing this arm-size cap alone: going from seven to eight increases the joint payments by less than 0.12% on either witness.

## Support and composition

* One actual nonfounder, **the charged apex**, with `65≤D≤100` and `5≤k≤10`. Families `hV`, `P4`, `Yt` only, supported production rows only. No participating descendant motif vertex is virgin `(L,k,h)=(1,0,1)`. The apex has maximum degree among participating vertices. The full planted fringe of the unique first child on the motif has at most eight ordinary vertices, excluding its parent. Sibling arms and the apex's own parent are unrestricted. An embedding is counted once, including each unordered `Yt` pair.
* The size-12 and size-13–20 higher charges cannot reach this apex: its planted fringe alone has at least `D` vertices. The four-core arms suffix is also disjoint because `k≥5`. The high-parent add-on prices `g,hC`, distinct charge families. The signed bare P4 correction has apex degree two. The moderate-old and old-flow main bands end at degree 64. For any shared-age or tail certificate of the same high-apex charges, use the maximum of certified prices on the overlapping objective, not their sum.
* Use the inherited [all-degree ceilings](../../denom/README.md) `1.356833021311` / `1.367169995210`. The proof prices the full per-occurrence charge, including the actual apex-parent input, and makes no second claim for a parent correction. The replay uses the two previously authenticated D18 blobs and checks layout, witness sides, cutoff and exponents; it does not rescan payload hashes or ceilings.

## Two physical comparisons

For a fixed apex `(D,k)` and a retained first arm `X`, condition on the complete physical filtration. Put `L=D−k−1`, `b=L/D`, and let `u` be the sum of the actual reciprocals of the apex's parent and the `k−1` other arms. Thus `0≤u≤k/2`. Apex gain has rate `L+1/D_X+u`; apex promotion has rate `b`; the apex sends into `X` at rate `1/D`. Intrinsic transitions of `X` use all actual internal arrows, with the outgoing `1/D_X` accounted for only in the apex gain. Neither the live input nor its evolution is assumed independent of `X`. Dummy labels collapse geometries only; no capped birth-mark rate is used as a live rate.

The pinned signed charge restricted to those embeddings is `a+c u`. For this tail apex each outside donor multiplies the **same** apex degree-gain derivative. For other motif vertices, the exact outside donors are already inside `X`. The per-family clamped affine toll `max(0,a+c u)` equals the true nonnegative toll on physical inputs; the clamp matters on the relaxed interval. In the joint solve take the hinge after summing the three affine charges. Backward comparison minimizes the one-step discounted Bellman quotient over the two endpoints and any interior hinge, since it is fractional-linear on each side. Crucially, use the same `u` in the charge, apex gain and holding denominator. All arithmetic is rational, rounded down on the `10⁻²¹` grid. Overshoots have zero continuation.

The resulting value at a virgin arm is a conditional price at its actual birth. Sum these prices over **all** apex promotions. Conditional birth prices telescope under birth compensation (each first-arm reward belongs to exactly one birth); no independence between birth times and the eventual successful arms is needed. The outer `(D,k)` recursion therefore adds the arm's entry price to the continuation for each promotion, and compares its apex gain over `[L,L+(k+1)/2]`. This proves a uniform lower discounted lifetime payment from a physical newborn at `HI=.565771874464`. The [physical stable-age law](../../renewal/next/README.md) multiplies the conditional discounted payment by `LO=.565676634491`. That law uses the actual physical birth-Palm distribution; the lower bound holds for every admitted birth environment, so no estimated mark weights are needed. It is an occupation floor for the same normalized per-core mean, not a separate residual-resolvent credit.

## Replay

From the Astra checkout (same environment for each invocation):

```sh
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1
P=./round2/.venv/bin/python
S=round2/occupation-neon/fringe/old-hub
nice -n 10 taskset -c 7 timeout 25s "$P" -Bu "$S/check.py" b3f18_lo.bin b3f18_hi.bin
nice -n 10 taskset -c 7 timeout 95s "$P" -Bu "$S/pilot.py" b3f18_lo.bin b3f18_hi.bin --arm 7 --detail
nice -n 10 taskset -c 7 timeout 100s "$P" -Bu "$S/pilot.py" b3f18_lo.bin b3f18_hi.bin --arm 8
```

The full-crown checker matched 135 exact charge and hazard fixtures **per side**, with nonzero examples of all three families. It uses a separate per-occurrence physical charge formula and projects actual full-tree arrows. It caught an earlier implementation that counted the common apex coefficient twice; the recorded results are from the corrected replay. The arm-7 diagnostic ran in 53.86 s (79 arm geometries, 74,655 Bellman rows per side); arm 8 ran in 58.47 s (194 geometries, 183,330 rows per side).

Arm-7 breakdown, rounded down per core:

| Lower-bound objective | Lower witness | Upper witness |
|---|---:|---:|
| `hV` | `2.23251581e−15` | `1.710466759e−14` |
| `P4` | `3.584019451e−14` | `2.5729066e−16` |
| `Yt` | `0` | `5.887166e−17` |
| Joint | `3.823166376e−14` | `1.747389944e−14` |

All cap variants and the sum of separate family prices are competing lower bounds, not additional credit on top of the joint result. A zero lower price is not a zero physical occupation claim. An independent reviewer has not replayed the derivation.
