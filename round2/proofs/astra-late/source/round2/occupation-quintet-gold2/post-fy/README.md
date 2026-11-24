> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Original Gold2 H3 threshold release input

The frozen lower-witness payment is **699435.940046906 ticks/core**, with exact normalized floor

```
1084423030952180046412413909 / 1704710415481443450880000000000000
```

in potential units/core (`S=2^40` ticks/unit). This is above `6e−7`; `699000` ticks/core is a conservative integer input to the combined release comparison. This is a stationary reward, not itself an exponent increment. The exact completed run is transcribed in [`threshold-lower.historical.log`](threshold-lower.historical.log). Its command was:

```sh
taskset -c 6 nice -n 19 timeout 105s /usr/bin/time -f 'wall=%e rss=%M' \
  python3 -Bu round2/occupation-quintet-gold2/post-fy/replay.py \
  --kind h3 --host 5 --stock 4 --receiver 79 --rstock 8 --donor 18 --tstock 3 --young 18
```

The release replay alias `--kind threshold --side lower` fixes these parameters, runs literal checks, and asserts the historical exact rational. Use only the original authenticated Gold2 pair and receipt in `round2/.work/replay-inputs`, clean detached revision `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`. No Gold1 numerical input enters the test. The exact prior endpoints `LO=35083855339718511/62017003637500000`, `HI=127993535256830043/226242120050000000` supply birth multiplier and discount. The lower payload uses the original trial `56571368/100000000` and outward sign −1.

The exact tick fraction is `1084423030952180046412413909/1550425090937500000000`. For the release ledger with zero singleton rebate, the required credit at `s=99/175` is `19217 + (99/175−56571368/100000000)·1808603193693 = 97537456765729/87500000` ticks/core, using the original Gold2 error and upper comparison denominator. With separate verified FY inputs of at least `395000` and `52000` ticks/core, the rounded total `1146000` exceeds that requirement by exactly `2737543234271/87500000` ticks/core. This arithmetic is conditional on those separate FY receipts and their mutual disjointness; this lane supplies the H3 input.

The same frozen mask on the upper witness returned **532790.0299087872 ticks/core**, exact `3304204122287698944634328103/6201700363750000000000`, or `3304204122287698944634328103/6818841661925773803520000000000000` in potential units/core. [`threshold-upper.log`](threshold-upper.log) is the captured stdout/stderr from `taskset -c 6 nice -n 19 timeout 130s /usr/bin/time -f 'wall=%e rss=%M exit=%x' python3 -Bu round2/occupation-quintet-gold2/post-fy/replay.py --kind threshold --side upper`: exit 0, 271350 states, 201436 paid comparisons, CPU `56.340769s`, wall `90.36s`, peak RSS `235588 KiB`. Verify code dependencies and both output logs with `sha256sum -c round2/occupation-quintet-gold2/post-fy/threshold.sha256` from `.`.

## Paid physical mask and production ownership

Pay one original H3 row per physical mother `r` → middle `h` → child `t` chain. The mother is a nonfounder, `19≤D_r≤79, c_r≤8`; middle `D_h≤5, c_h≤4`; child `D_t≤18, c_t≤3`. Actual labels include each core's own parent in `c`. All ages and permitted additional births are included, with the original degrees and arities at the time of payment. The degree ranges make the mother strictly older in the sense of the degree gate, not a permanent birth-degree mark. Founder-rooted chains and two-child arms earn nothing in this mask.

For signed original H3 row `H`, pay exactly `R=max(0,−H)` within the mask. At a middle vertex, `g→H3` transfers are one per removed tip, with actual `(1/D_tip−r1_middle) Δg_retained`. Their sender retains `r1_middle Δg` and subtracts the transferred portion. The original `H3→P4/S3` terms are subtracted once whenever the receiving motif is supported at the **pre-gain** labels, and all other external reciprocal contributions remain physical. The newborn P4/S3 term stays in the transmitting H3 row, computed at the actual post-promotion arity. In particular neither the envelope overhang nor a higher family cap is booked as additional reward. For a nonnegative production cap `K` on that occurrence, `H+R≤K`, so the unchanged host ledger covers the remaining outward charge. The old lower witness H3 family maximum is 408 ticks; no per-occurrence cap tightening is used here.

**FY disjointness:** original FY atoms are only at a non-M3 endpoint of a paying `g` edge. Every incoming `g→H3` transfer here occurs at the H3 middle, which is in M3; it is an original inter-level gain allocation and contains no FY term. H3 payments use only the already allocated H3 row. Singleton payments and edge cap residuals are distinct row slots. Do not add any diagnostic number or overlapping H3 reward. No standalone endpoint update or singleton rebate is taken here.

## Lifetime proof and verification

Use three registration levels: a nonfounder mother at actual birth, each actual mother→middle edge at its birth, and each actual child of that middle. An extra birth advances the original prefix or complete row to the post-event arity and, when relevant, seeds a new marked prefix independently. Mother and child are live throughout; their current `1/D` is used in both the payoff and transitions. Other exterior contacts enter the corresponding gain hazards as live reciprocal sums. Where the charge has a slope on that sum, the same variable enters the reward. The directed integer Bellman comparison is checked over the full exterior intervals. For a clipped payoff crossing zero, the solver takes the better certified subsolution for zero and the uncut affine reward; both are pointwise below the positive part. No arbitrary box discrepancy is itself a payment.

The physical gain rate is `D−c+Σ 1/D_neighbor`; promotion rate is `(D−c)/D` with old `D` fixed, old `c` incremented, and actual newborn `(2,1)`. Reverse degree/arity makes every within-mask update triangular. Summing the checked `HI` Bellman inequalities over physical registrations and applying stable-age balance gives `π[R]≥LO·V_mother(2,1)` in ticks/core. Payment is unique by genealogy even with repeated arm labels. The completed threshold lower solve checked 271350 path states, with 187105 states using the nonzero-piece reward comparison; CPU `52.486688s`, wall `59.39s`, peak RSS `232028 KiB`.

Separate completed fixture replays with `--kind g --fixtures --limit 12 --wing 32 --arm 4` passed on each witness: 37 H3 rows, 40 edge rows, 64 H3 outgoing transfers, 74 edge outgoing transfers, four FY-gap identities, and eight unique branching-path registrations. The small edge solve removes **all** FY gaps before pricing, and is optional comparison only; the signed threshold release input above is solely H3.

The Gold2 mother-high diagnostic from 835c guided the choice of mask but does not enter the certificate. The threshold input and its historical evidence are frozen at `79/8,5/4,18/3`. Wider masks are comparisons and do not replace this receipt.

## Independent audit and exact release-coordinate arithmetic

The d4d3 cold audit verdict was **HOLDS** ([audit](../../../../../../reviews/cold/astra/checks/gold2-h3-sol.md), [replay](../../checks/gold2-h3-sol.py)). It freshly replayed the two separate original Gold2 FY floors, supplying `39526146756920663734311/100000000000000000` and `1046834122639903380009/20000000000000000` ticks/core. The first pays at exact-G2 non-M3 host `D=19…251,c≤17`, payer `D≤5`, and any distinct donor (absorbed at `D=12`); the second pays at exact-G2 non-M3 nonfounder host `D≤150,c≤10`, payer `D=6…34,c≤7`, distinct donor `D=2…5`, and `max(D_host,D_payer)>18`. Their paying-edge degree ranges are disjoint.

With just these two floors and the frozen H3 floor, the exact total is `142271857779342357237253374909399/124034007275000000000000000` ticks/core, and the margin over the `99/175` requirement is `28066074690162200126495624365793/868238050925000000000000000` ticks/core (about `32325.32214`). With zero singleton rebate, the resulting original-Gold2 lower rational is `126905728962161948969119589499374909399/224328301684105796116575000000000000000`, about `0.5657143035873726`. This arithmetic adds no other H3 credit and does not itself authorize release or change the board.
