> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Original Gold2 lower: unused FY and post-FY edge credit

**Book exactly 19,240.878325014863 additional ticks/core** alongside the existing disjoint lower `FY ∪ g18` and a valid singleton pressure replacement. This falls **10,759.12 ticks/core short** of the 30k target. It is not a new endpoint or a ceiling on the remaining edge account. The certified edge total with the old three floors is **467,078.011494398 ticks/core**; no additional small-payer FY is claimed here.

## Original cap after stripping all FY

Let `G` be the **charged** outward original physical lower row on a single edge, with `σ=−1`. It includes singleton transfers, debits to H3 at M3, and the original newborn H3 at the post-promotion mark. The uncharged `L g − λg` alone is not the relevant capped row. For exact-G2 non-M3 endpoint `x` and physical payer `y`, let `j=Δ_x g(x,y)` and `A` be the original 20-knot spline. The ordered FY atom at distinct donor `z` is

```
max_k σ(ξ_k j − A(ξ_k)) − σ(j/D_z − A(1/D_z)).
```

This is nonnegative; adding every such atom to `G` leaves the original signed primal spline terms and replaces each external FY contribution by its knot maximum. Write the result `G*`. The cap claim is `G* ≤ 792`, hence `792−G` dominates every union of distinct FY atoms **plus** the residual `792−G*`.

At a mixed edge use `t=252/D_tail`. The finite endpoint jump is `(1−t)j_INF+t j_knot`; the convex knot maximum is dominated by the corresponding mixture of knot maxima. The primal spline at the tail reciprocal is affine on `[0,1/252]`, as is its singleton transfer, and projection preserves the original H3 newborn weight. In the signed adjacent-tail send, use `t max(σj_INF,σj_knot)/252`, including the finite endpoint's own `251→252` successor. The tail coordinate contributes `tγ Δ_signed g` and `tβ Δ_signed g` in the knot column, with `γ=receive/(D_tail+1)` and `β=(D_tail−c_tail)/D_tail`. These are exactly the physical ranges used in the original C++ verifier; overflow clamps only the weight index, and its promotion jump vanishes. At an overflow finite endpoint the two-way interior send is still covered, and there is no FY spline there. Thus the lifted knot/INF proof covers `G*` at mixed edges, without moving the signed maximum outside the expectation in the wrong direction. Both G2-tail endpoints imply a zero charged row and zero FY: the entire 24×24 corner face is zero, including all clamped tail classes. The checker checks 4,624 charged-row identities and cap inequalities, including 1,056 comparisons to the lifted mixed-tail envelope.

## Disjoint union and units

All FY masks below use physical distinct `(host, payer, donor)` roles, exact-G2 non-M3 host and `max(D_host,D_payer)>18`. All payers have `D>5` and hence are outside 2f74’s account. Power-of-two scale is `S=1099511627776` ticks per potential unit. The FY replay returns normalized units after one birth normalization and one division by `S`; these are multiplied back by `S` exactly once. Post-FY routines already return ticks/core.

| Additive region | Exact ticks/core | Rounded ticks/core |
|---|---:|---:|
| Donor `2..5`, payer `35..125,c≤3`, host `D≤251,c≤12` | `58162652610249882843115293/6201700363750000000000` | 9,378.500927 |
| Donor `2..5`, payer `6..34,c≤8`, host `D≤251,c≤18`; at least one of `D_host>150`, `c_host>10`, `c_payer>7` | `7127543113684417780424361/1240340072750000000000` | 5,746.442665 |
| Donor `D≥6`, payer `6..40,c≤5`, host `D≤251,c≤11` (donor absorbs at 9) | `25026790798544028667583523/6201700363750000000000` | 4,035.472424 |
| Non-FY edge residual | `6316495315166525735295397/78502536250000000000000` | 80.462309 |

The middle line excludes the entire previously booked large-payer mask within its coverage. The three new FY masks are mutually disjoint by payer range or donor range. The original degree-18 full edge credit cannot overlap an FY with paying-edge maximum above 18. The remaining residual line pays `792−G*`, so it leaves every FY atom available even when a path is geographically in the same tree. Its finite cut pays both endpoints `D≤25` or one endpoint `D≤3` and both `D<252`, always with `max(D)>18`; its mixed cut pays one endpoint `D≥252`, other `D≤10`; its double-tail cut pays only edges with parent already `D≥252` at registration once the child also reaches 252. These three subregions are edge-disjoint at payment. The mixed reward retains only `t` times the nonnegative knot residual; the INF remainder is discarded. The high-coordinate continuation is `t` times a clamped table; after dividing the inequality by `t`, `γ` multiplies the decay of the continuation directly. No affine-reward width scaling from either disputed H3 sweep is used.

The exact additional total is `119326162107114155945908957363/6201700363750000000000000` ticks/core, equivalently `119326162107114155945908957363/6818841661925773803520000000000000000` potential units/core (about `1.74995e−8`). These are stationary credits, not exponent displacements; a global margin must use the corresponding Gold2 lower comparison denominator and the selected nonsubstituted credit ledger. Smaller pilot masks overlap the final masks and are not added.

## Physical replay

All bounds use original Gold2 exact `HI` as discount, exact `LO` as birth normalization (except the very small double-tail charge, which safely widens both endpoints). The FY sweeps register both mother/child orders and both sibling roles, pay only at actual marks and evaluate exact adverse corners; a donor in the degree-9 absorbing class has reciprocal in `[0,1/9]` and retains its physical identity. The fixture and pricing runs authenticate the full original lower SHA-256, the pinned receipt and source `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`. Neither Gold1 numerical inputs nor the H3 payment replay are used.

```sh
R=round2/occupation-quintet-gold2/true-drift/lower/check.py
python3 -IBu "$R" --mode cap
python3 -IBu "$R" --mode wide-payer --host 251 --stock 12 --payer 125 --pstock 3 --donor 5
python3 -IBu "$R" --mode wide-host --host 251 --stock 18 --payer 34 --pstock 8 --donor 5
python3 -IBu "$R" --mode old-donor --host 251 --stock 11 --payer 40 --pstock 5 --donor 9
python3 -IBu "$R" --mode post-fy
python3 -IBu "$R" --mode union
```

The completed exact sweeps checked 13,997,830, 13,995,828 and 16,703,820 FY path states, respectively; the post-FY finite sweep checked 209,155 pairs. The cap and all pricing runs ended `PASS original Gold2 lower pins and exact comparisons`. The union command reproduces both exact totals above.
