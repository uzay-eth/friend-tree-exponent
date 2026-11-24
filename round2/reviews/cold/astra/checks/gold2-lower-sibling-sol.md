# Cold audit: original Gold2 lower compact sibling floor

**HOLDS**, sibling component only: `16363762149897743918379249111/65029541606195200000000`
payload ticks/core = **251635.8219006546**, about `2.28861e-7` in potential
units/core. The audit does not validate the upstream or box amounts or
claim restoration of `99/175` on a composed ledger.

## Physical reward and two-child check

The paid path is a nonfounder mother `h` → middle `m` → child `y`, plus a
distinguished actual other child `s` of `h`. At payment, `19≤D_h≤79,c_h≤8`,
`D_m≤5,2≤c_m≤4`, `D_y≤5,c_y≤3`; registration requires `c_h≥3` including her
own parent. The sibling is live through `D_s=2,3` with all possible arities,
then remains a real known neighbor with reciprocal in `[0,1/4]`. Its own
transitions in that abstract state do not change the test function. The
first three roles are killed only at irreversible caps; below the mother's
entry degree the comparison continues with reward zero.

Write the original lower outward row as `P=A+Σ b_i ρ_i` in payload ticks,
using `σ=-1`, the original nominal trial, actual labelled gain and promotion
terms, original transfers, and newborn extensions at post-promotion arity.
The independent code uses the authenticated physical reader, not the price
or cohort solver. Let `B` be the row at the adverse corner for all exterior
slots (`c_h−1,c_m−2,c_y−1`), `b` the mother slope, and
`z_s=max(0,b/2)−b/D_s≥0`. Pay `F(z_s)`, where
`F(z)=[z−B]_+−[−B]_+`. For overflow use the endpoint giving the smallest
`z_s`. All other slots remain adverse **in the reward**; their actual rates
are still enclosed in the generator. `F` is convex, nonnegative, and vanishes
at zero, hence `Σ F(z_s) ≤ F(Σ z_s) ≤ [−P]_+−[−B]_+` for all siblings of
the same path. The last inequality follows because replacing the remaining
adverse slots by their actual incidences cannot increase the row. This
preserves the original cap and does not pay the same H3 negative row twice
when two siblings exist. In particular the sibling credit is not disjoint
from an unrestricted full `[-P]_+` price on overlapping paths.

The audit recomputed 2,048 labelled physical rows with two actual sibling
children, including 1,200 where both one-child rewards are positive and 416
where the sum is strictly below the joint replacement. These controls include
exact and overflow sibling degrees and uncovered/covered exterior incidences.

## Chronology, units, and replay

Exactly three birth orders are possible: `hmys`, `hmsy`, `hsmy`. A real
promotion both continues an existing prefix at post-promotion arity and
registers a longer prefix at the designated parent, adding a newborn
`(2,1)`. The auditor enumerated all 720 increasing-parent trees on seven
vertices and checked that registration matches direct genealogical
enumeration, with 197 incidences per order and no duplicates.

For the scaled Bellman inequality, the reward is `floor(2^20 F(z_s))`
times the common rate denominator. Gain ports carry the **unshifted**
successor only: the robust reward has no remaining slope in these ports.
The exact residual is checked over the full rate box at every state. Using
the original upper prior endpoint as discount and original lower prior
endpoint as birth multiplier, the three birth-order fractions match
`compact.log` exactly. To report normalized potential units divide payload
ticks by `2^40`.

```sh
taskset -c 5 nice -n 19 timeout 100s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /usr/bin/time -f 'wall=%e user=%U rss=%M exit=%x' \
  round2/.work/venv/bin/python -IBu round2/checks/gold2-lower-sibling-sol.py
```

Observed: `HOLDS sibling only 16363762149897743918379249111/65029541606195200000000 ticks/core 251635.8219006546`; `wall=24.15 user=23.20 rss=314900 exit=0`.
The replay fresh-hashes the lower payload against
`41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`, checks
the authenticated original receipt and immutable file identity, and does
not open the upper payload. Reviewed SHA-256s: `compact.log`
`c3d7d979fe82c4ac626cbdbef5f461411e966abd4e53ab35529a9cd5bb75c02a`;
`price.py` `589e5d8d14f91621d0f5ad6bfecfdc00e4018a05bb3100416d0de8b22a4587f1`;
`higher/cohort.py` `b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c`.
