> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# High-stock occupation of raw higher-row slack

The matched short replay gives additional stationary slack lower bounds of about `3.27336e−7` (lower witness) and `3.43615e−7` (upper witness). Added **only** to first-body plus higher size `≤12`, these give the outward bracket **`[.565675822028,.565772451624]`**. The installment from genuinely unbounded stock is much smaller than the finite `size>12` installment.

| Witness | With zero at `D≥18` | With unbounded high-stock payment | Unbounded gain in stationary reward |
|---|---:|---:|---:|
| Lower | `3.273310331636669e−7` | `3.273357538588808e−7` | `4.720695214e−12` |
| Upper | `3.436094971804806e−7` | `3.436147907300660e−7` | `5.293549585e−12` |

Displayed reward decimals are approximations to the exact fractions from the replay. The final transfer uses the original witness reward units and the frozen `F_DEN` ceiling from `higher/extract.py`, not an interpolated new certificate. The two inward endpoint moves versus `.565675597658,.565772677869` are approximately `2.24370e−7` and `2.26246e−7`.

## Proof and scope

Only top-owned `hV,hC,P4,Yt,Yb,S3` are paid, at **physical fringe size greater than 12**, on at most four cores. A two-core fringe can be visited but owns none of these higher rows. Interior rows are imported from the exact decomposer and cross-checked against a Bernstein-extension formula. The definitive objective is the nonnegative original signed higher-row slack; there is no signed unknown terminal payoff. The proof takes the two physical exterior-rate endpoints jointly for the reward and the root gain before rounding the Bellman row down. Other gains and promotions are actual intrinsic moves, including both reciprocal directions between each pair of neighbors. At each first promotion the full actual post-promotion fringe is carried forward if still admitted; discarded paths give zero, without resetting the sibling or its age. This is the embedded-clock integral of the joint physical gain/Palm law, not a product of child histories conditional on a realized root path.

The extension evaluates the physical Bernstein motif weights at exact degrees, including gain and promotion successors. For transferred edge terms where the receiver is in the degree tail, it uses the actual production controls `3t²(1−t)(F0−F1)/18+t³Yref`, with the other arm Bernstein-mixed as required. Transfers from smaller rows are included once inside the chosen higher row, without booking any smaller motif reward. On explicit rows both exterior endpoint rewards are checked nonnegative. The two witness payloads are the existing production inputs from the verification record (`round2/attacks/neon1/outputs/production-final/session.log`, historical reference); this replay does not stream their multi-gigabyte hashes again.

After an explicit cap exit, the unbounded positive payment is used if exactly one core has `D≥18`, the others have degree at most five, and there are three or four cores. Keep the joint shape and all other stocks fixed until the first change outside the large vertex's own gain channel, with zero payoff thereafter. Allow the large stock's gains indefinitely. For either exterior-rate endpoint `r=0,1/2`, put `t=18/D`; the raw slack times `(18+t)^3` is a polynomial of degree at most seven on this one-coordinate face. Indeed the old/promoted weight is cubic, the one gain replacement uses `t'=18t/(18+t)`, degree gain rates have at worst a simple `1/t` that cancels against the gain difference, and reciprocal/promotion rates are affine in `t`. Incoming edge transfers on tail receivers are cubic and the center and cutoffs are fixed on this face. The checker interpolates with exact rationals, independently checks at `D=27,43`, factors a common `t^a` if needed, and lower bounds the remaining polynomial by its Bernstein coefficients on `[0,1]`. A failed positive coefficient test supplies zero, not a signed tail.

For a positive face with lower slack `c t^a`, the candidate `c t^a/(HI+κ+a)` is a lower continuation. `κ` bounds all exit intensities (including promotions and actual rate caps for sends affected by the large degree); the large vertex's gain rate is at most `D`. Its drift on `D^(−a)` costs at most `a` in relative units. Exterior root gains are accounted for by the same uniform hazard or, if at the large coordinate, also fit within the physical gain bound. Thus the pointwise comparison holds for arbitrary adapted parent input, with unrestricted Yule growth. Cores and size are monotone, so these face rewards cannot overlap the existing size-12 calculation. The stable-age transfer multiplies the lower uniform birth integral at discount `HI` by `LO`. Bounded-local approximations of the nonnegative face toll increase to it by restricting the large degree and taking limits; no unbounded raw observable is required in a finite-radius stationary convergence step.

**Do not add** this result to `renewal/next` higher/three/four/both. 4e4f identified overlapping `size>12` higher rows there. An independently verified *g-only* numerator is disjoint, but no such number is included here. For competing proposals on the same higher reward, combine by a verified pointwise maximum, or compare lower bounds on the same objective. No 19-mark continuation is added here.

## Replay

For each `SIDE=lo,hi` run from this directory:

```sh
taskset -c 7 nice -n 5 timeout 32s /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M' ./round2/.venv/bin/python -Bu replay.py b3f18_${SIDE}.bin --cap 20 --degree 17 --tail --unbounded
```

Omit `--unbounded` for the matched zero-boundary control. Both production runs visited 42,135 contexts, paid 9,679 interior contexts and admitted 19 positive high-stock boundary payments. CPU and peak RSS: lower with tail 6.42 s / 213,452 KiB, upper 6.32 s / 213,332 KiB; zero-boundary lower 8.06 s / 210,864 KiB, upper 8.59 s / 210,712 KiB. Earlier exploratory pilots used different caps and are not combined with these runs. All commands were niced and pinned to CPU7.
