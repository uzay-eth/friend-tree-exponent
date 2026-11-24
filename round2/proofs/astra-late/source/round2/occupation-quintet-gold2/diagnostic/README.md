> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# GOLD2 physical census for the `99/175` target

On the original GOLD2 **lower** witness the largest targeted higher-row pool is H3: **2.571 / 2.550 ppm per core** pooled / held out on nonfounder H3 with middle arity 2–4 and maximum tip degree 19–1998, excluding the old strict-chain mask. The chronological mother-high orientation alone accounts for **2.268 / 2.236**. The broad P5+B5+S4 pool also clears 2 ppm: **3.222 / 3.464**. These are negative charged rows on preserved finite trees, not certified gains toward the roughly `5.62e−7` additional normalized numerator specified for the current stationary `f+host` ledger.

| Lower-witness physical slice | Pooled | Held out |
|---|---:|---:|
| All H3 | 3.682 | 3.693 |
| Priority H3, nonfounder and outside old strict mask | 2.571 | 2.550 |
| Priority H3, middle degree ≤5 | 2.452 | 2.428 |
| Priority H3, mother-high | 2.268 | 2.236 |
| Compact mother-high 65/18/5, mother arity ≤7, child arity ≤3 | 0.892 | 0.892 |
| P5+B5+S4 | 3.222 | 3.464 |
| B5 alone | 1.730 | 1.676 |
| S4 alone | 1.102 | 1.398 |
| Edge negative charged row | 1.597 | 1.573 |
| Edge FY knot gap outside diagnostic masks | 0.927 | 0.907 |
| Edge uncovered non-FY gain gap | 0.094 | 0.093 |
| Edge cap slack after subtracting FY and uncovered box gaps | 0.290 | 0.289 |

All table entries are `10^-6` normalized potential units per core; nested rows and edge subaccounts are not additive. For the broad H3 priority cut the lower witness's degree bands 19–37, 38–179, and 180–1998 carry **0.863, 1.122, 0.586** ppm pooled respectively. At mother-high with child degree 2–5, those bands carry **0.451, 0.428, 0.238**. With child degree 6–11 the first two bands add **0.176, 0.187**. A compact-only proof leaves most of the signal outside its mask; keeping middle degree ≤5 retains 95% while allowing a broader range at the tips. The post-FY edge residual is mainly on endpoint pairs with degrees 2–5 and 180–1998, max arity 5–8 (**0.180 / 0.179**); the edge physical negative reward is largest at older mothers with degrees 38–179, arity 5–8, and children of degree 2–5 and arity 2–4 (**0.347 / 0.328**). The upper witness is tabulated in `report.log`; its priority H3 is 3.992 / 3.861 and its quintet pool is 3.457 / 3.651.

For a stationary follow-up, the higher-row payment must be tied to the original charged H3 or quintet row and the **retained local host cap**. A measured negative row can support a cap rebate where a proof certifies its stationary occupation; it does not itself price that occupation. In particular the large H3 pool is mainly on low-degree middle hosts already present in the `f+host` LP, so disjointness cannot be argued from middle degree or from the mother being outside the LP degree cap. Keep each actual negative-row increment distinct from any already-booked adverse-box replacement or host saving. A larger-window H3 pricing and separate B5/S4 row pricing are the two measured pools above 2 ppm; the edge gaps are smaller alternatives.

## Inputs and exact replay

- Snapshot `round2/.work/replay-inputs/`, detached source `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`, receipt SHA-256 `5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284`.
- Two read-only 6,344,902,168-byte blobs. Lower SHA-256 `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`; upper `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`. Blob digests are inherited from the cold audit's authentication log. Replay rechecks the source, JSON receipt and pre/post blob identity/stat. Exact trial points `56571368/10^8` and `56573686/10^8`; `DF=4000`, G2 cutoff 252, scale `2^40`.
- Same eight age-20 inputs `numeric4-census-ct-48125000{0,1,2,3,4,5,6,7}.tsv` (896,887 pooled cores), and separate age-20 `numeric4-root-tree.tsv` (235,709 cores). The logs give per-tree counts and SHA-256. No simulation.
- `bash round2/occupation-quintet-gold2/diagnostic/replay.sh` exited 0. Lower/upper census wall time 14.70/17.69 s. Across 18 witness/tree comparisons the largest absolute difference between charged, raw, and event-generator sums was `1.564e−12` normalized potential units **before** dividing by core count.
- The separate sparse rational oracle tested four literal ordinary-tree fixtures per witness, including singleton and G2 seams, all eight signed transfer channels, and 12 actual ordered FY slots per witness. Another 700 exact signed knot evaluations per witness cover all 20 GOLD2 knots (denominator 10080), interior interpolation and ordered cancellation. Fixtures compare exact charged/raw totals and literal direct generator before comparing with the census. See `exact.log`.

No GOLD1 payment was used, no endpoint was changed, and the census provides no stationary lower bound.
