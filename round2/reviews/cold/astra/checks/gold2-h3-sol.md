# Cold audit: Gold2 H3 and two FY floors at `99/175`

**HOLDS**, taking the original Gold2 production verifier and comparability certificate at its pinned receipt. The original lower witness plus exactly these three credits gives

```math
\mu\ \ge\ \frac{126905728962161948969119589499374909399}
                  {224328301684105796116575000000000000000}
       \,>\,\frac{99}{175}.
```

The lower bound is about `0.5657143035873726`; the strict gap is about `1.78731e-8`. This uses no singleton rebate, upper H3 price, or Gold1 numerical input.

## Physical row and original caps

The inputs are original Gold2 at source `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, with `S=2^40`, Gold2 lower outward sign `−1`, and lower payload SHA-256 `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`. The large-payer replay rehashed the lower payload. The threshold replay checked the certificate hash, detached source HEAD, and stable file identity. This audit reuses the pinned production caps rather than rerunning the exhaustive C++ verifier.

The H3 occurrence is the unoriented set underlying a chronological nonfounder mother → middle → child path, with caps `79/8, 5/4, 18/3` and mother degree at least `19` at payment. This picks at most one occurrence from any given three-label set: the middle has just one parent. Equal child labels do not identify distinct children. Extra children are consistently included in `c` and in the external receive sum. Both endpoints are exact G3 within this mask (`88/10`), including their one-step successor marks. A one-step exit from the **occupation** mask can still have nonzero **production** H3 weight and is evaluated accordingly in `Rows.triple`.

Against pinned `verify_quintet.cpp:630–719`, the H3 outward row `H` consists of signed own gain/promotion, the nominal subtraction, and two transfers from the incident `g` rows. For a middle gain from one of the two path neighbors, the transferred term uses that donor's `1/D` less the middle's `r1`, times the gain jump of the edge to the other neighbor. The C++ source envelopes certify an upper bound on these terms; the physical row uses the terms themselves. At any path vertex, an external gain is retained at `r1` where that exterior neighbor extends the pre-event motif to supported P4 or S3; otherwise the true gain is retained. The extension support tests depend on internal marks only: P4 uses `24/7` and internal sum `≤32`, S3 uses `20/7` and requires at least three arms. At a promotion the transmitting H3 row keeps the new P4/S3 weight, using post-promotion arity. P4/S3 receive the paired gain transfers only.

For each H3 occurrence use its **original** local middle cap `K_m=max(0,cap_h3_m[m])`, not the global H3 maximum `408` in place of that local cap. The original verifier sums `binom(c_m,2) K_m` at the middle (`1252–1259`). Since `H ≤ K_m`, the selected payment `R=max(0,−H)` gives `H+R=max(H,0)≤K_m`; unselected occurrences are unchanged. Thus the original positive aggregation is valid with the credit subtracted, including middles whose original signed cap was negative. No higher-row envelope overhang is booked.

## Occupation and birth normalization

An augmented physical generator acts on the sum of nonfounder singleton, ordered edge, and ordered path continuations. Its tracked gain at label `u` has rate `D_u−c_u + Σ 1/D_v` over core neighbors; promotion has rate `(D_u−c_u)/D_u`, leaves `D_u` fixed and creates `(2,1)`. Exterior events with no tracked mark update contribute zero to a continuation jump. At a promotion the old registered objects advance and any admissible new objects are seeded independently. For H3, a mother promotion seeds `(r⁺,J)` and preserves the updated singleton `r⁺`; a middle promotion seeds `(r,h⁺,J)` and preserves the updated prefix. A promotion of another vertex in the path changes its mark only. The edge table's newborn tip is in the third slot; the mother's first edge is read with the mother's post-promotion mark. Saturation `c=D` gives zero promotion hazard. Monotone upper cutoffs make exits permanent; the lower mother-degree payment gate is **not** a kill gate. A nonfounder singleton at every true core birth is the sole unallocated source.

`H3.layer` bounds the unknown external sums by `[0,(c−tracked_degree)/2]` and couples the same sum to both physical gain and the affine row reward. Its directed rounding and exact cross-multiplied corner residual certify the Bellman inequality at every point in those boxes. Taking the larger of subsolutions for zero and `−H` is safe for the clipped payoff: at the current state both tests use the same fixed successor values, and `max(0,−H)` dominates either branch pointwise. All stored continuations are nonnegative, so moving the discount from original `HI` down to `μ` preserves the inequality.

For clarity, let `A` be the sum of continuations over currently registered objects, in ticks, and `B1` the physical promotion rate at the observation root. Summing the Bellman inequalities gives `R_total ≥ μ A − L A + V(J) Σ B1`. This includes all independent prefix seed terms in `L A`, with just `V(J)` per true birth added back. The root law is the **physical** uniform-core limit `π` of `fringe/spatial-limit.md`, with the live exterior and chronological ancestry. Equivalently, to justify applying this balance without postulating an invariant product law, take finite-tree expectations in the displayed inequality, multiply by `exp(−μt)`, integrate to `T`, and divide by `T`. The weighted `μ A − L A` term is just the endpoint difference divided by `T`. It vanishes: the number of objects with nonzero continuation per core is uniformly bounded, and `exp(−μt) E C(t)` has a finite positive limit by the spatial mean expansion. The bounded local reward has that same limit factor times `πR` (one may allocate it at its middle). For the source, `E Σ B1 = d E C(t)/dt`, so the same integration proves `π B1=μ`. Hence `πR ≥ μ V(J) ≥ LO V(J)`. This is also the birth-Palm normalization with newborn law `πB/μ`, where the post-birth trajectory follows the same born core with the live physical neighbors. No independence among marks, clocks, or ages is used.

## FY ownership and sign

At an exact-G2 non-M3 host each ordered distinct `(payer,donor)` contributes one gap between the knot maximum of `σ(x Δg−A_h(x))` and the value at `x=1/D_donor`, `σ=−1`. Piecewise linearity gives a nonnegative gap. The C++ edge bound uses the knot maximum (or its convex corner envelope); replacing it by the actual value and crediting the gap stays within the original edge cap. The paired spline terms cancel when the ordered incidences are summed; H3 uses only the separate M3 transfer slots. These claims use the original 20 knots, not any Gold1 grid. The two masks are disjoint by payer degree:

| Credit | Paid mask | Ticks/core |
|---|---|---:|
| H3 | Chronological path above | `699435.940046906` |
| Small FY | Non-M3 exact-G2 host `19..251`, `c≤17`; payer `D≤5`; any distinct donor | `395261.4675692066` |
| Large FY | Non-M3 exact-G2 host `D≤150,c≤10`; payer `6..34,c≤7`; distinct donor `2..5`; `max(D_host,D_payer)>18` | `52341.70613199517` |

Both FY passes register the two mother/child role orders and both ordered sibling roles. They seed continuations at actual promotion marks. The small branch's degree-12 absorbed donor cannot leave the class; its ongoing mark changes are self-transitions of the tracked state and its reciprocal stays in `[0,1/12]`. Minimizing the reward over that interval's endpoints and interior original knots is valid even if that clock and receive rate are correlated. Omitting an already-absorbed prospective mother's registration source only lowers the floor. Both FY discounts use decimal endpoints directed outward from the exact Gold2 enclosure.

With `M=π(SΦ allocated per core)` and `m_min≤M≤m_max`, the corrected ledger gives `(μ−λ)M ≥ credit−19217`. This follows also by the finite-tree weighted integration above, with the existing bounded per-core allocation for `SΦ`. The right side is positive, so use `m_max=1808603193693`, **not** `m_min`. At `λ=7071421/12500000`, the exact required credit is `97537456765729/87500000` ticks; the exact three-floor surplus is `28066074690162200126495624365793/868238050925000000000000000` ticks (`32325.322139776385`).

## Replay evidence

All commands below ran from `[local path]` and exited successfully:

* `sha256sum -c round2/occupation-quintet-gold2/post-fy/threshold.sha256` — nine `OK` entries.
* `python3 -IBu round2/checks/gold2-h3-sol.py` — 926 literal directed arrows, 1,003 newly seeded registration objects across hypothetical promotions, 1,127 continuation cutoff exits, 11 H3 row equalities (all negative), 22 incoming and 27 outgoing transfers, 56 nonzero newborn terms, 42 ordered terminal controls, and 2,240 FY knot checks. A planted recipient-denominator substitution in the donor transfer was rejected by all 11 H3 comparisons. Sparse reader touched 1,727 words.
* `nice -n 19 python3 -IBu round2/occupation-quintet-gold2/post-fy/replay.py --kind threshold --side lower` — 271,350 states, 187,105 paid rows; normalized `1084423030952180046412413909/1704710415481443450880000000000000`.
* `nice -n 19 python3 -IBu round2/occupation-quintet-gold2/edges/fy.py --host 251 --stock 17 --donor 12 --absorb` — 3,585,120 states; normalized `39526146756920663734311/109951162777600000000000000000`; `PASS exact original Gold2 knot and physical Bellman replay`.
* `nice -n 19 python3 -IBu round2/occupation-quintet-gold2/payer/price.py --host 150 --stock 10 --receiver 34 --rstock 7` — 3,946,320 states, 516,312 exact reward labels; normalized `1046834122639903380009/21990232555520000000000000000`; `PASS pinned lower bytes and exact Bellman inequalities`.
