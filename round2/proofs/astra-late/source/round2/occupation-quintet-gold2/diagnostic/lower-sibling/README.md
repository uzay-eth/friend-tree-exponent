> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# GOLD2 lower: actual maternal sibling and upstream H3 recovery

The exact floor for the combined non-box payment is **278,435.854688312 ticks/core**, about `2.53236e−7` normalized potential units/core. It replaces the corrected `0.4214974850656785`-tick 19–79 floor; the corrected `14,132.789764062283`-tick 80–100 seam remains disjoint. No original cap is paid. No endpoint is asserted by this receipt.

| Payment | Exact ticks/core | Approximate |
| --- | --- | ---: |
| Actual distinguished maternal sibling, three orders | `16363762149897743918379249111/65029541606195200000000` | 251635.8219006546 |
| Unique actual upstream of the mother | `1742793847212361294561574577/65029541606195200000000` | 26800.032787657383 |
| Combined non-box | `18106555997110105212940823688/65029541606195200000000` | 278435.854688312 |

## Pointwise budget and registration

The paid row is an original lower signed H3 charge on an actual nonfounder chronological path `h→m→y`, with current mother `D=19…79,c≤8`, middle `D≤5,c≤4`, and young `D≤5,c≤3`. The optional upstream payment starts only at actual nonfounder upstream births. The additional child is a child of `h` distinct from `m`. No condition on elapsed age enters the mask.

Write original outward charge `P=A+Σ b_jρ_j`. Let `B` be its upper bound obtained by maximizing each uncovered physical port over `[0,1/2]`; it is an adverse row value, not the production cap. For each real extra child and for the unique upstream, set `z_i=max(0,b_h/2)−b_h/D_i`. These reductions are nonnegative, and `P≤B−Σz_i` with unpriced ports still adverse. Pay each retained donor `f(z_i)=[z_i−B]_+−[-B]_+`. Convexity with `f(0)=0` implies `Σ f(z_i)≤f(Σz_i)≤[-P]_+`. This bound is on the sum over all extra children and the upstream, even if several contribute to the same H3 row. Thus the original nonnegative production cap still covers `P+Σf(z_i)≤max(P,0)`; the host savings are not paid here. The separate box diagnostic in the log is not in the reported sum.

Each distinguished sibling is registered exactly once across `hmys`, `hmsy`, `hsmy`; a different sibling is a different incidence. `uhmy` registers at most one upstream. The labelled registration oracle checks 26 unique sibling incidences and two upstream paths. Each auxiliary live donor is tracked through degree 3; at degree 4 it enters a permanent inverse-degree interval `[0,1/4]` retaining the same physical identity. Its marked neighbor's gain uses the same conservative interval for that donor; after upstream overflow no further marked children are seeded from that prefix. The tracked mother, middle and young pay only before their irreversible degree/arity exits.

## Exact arithmetic and reproduction

The lower payload is authenticated by the pinned Gold2 reader on the same open inode. `price.py` reuses the unaffected constant-reward sweep from the pinned cohort source. The defective `kind='physical'` shifted slope is never used: the pointwise robust subreward is determined before the recurrence and rounded toward zero on a nonnegative `GRID=2^20` grid in **raw witness ticks**. Clock denominators are exact: `den=lcm(HI.denominator,2,…,80)`, discount `HI`, birth factor `LO` once. The summed prefix Bellman inequalities cancel actual registration source terms, leaving the nonfounder birth value `LO·V(2,1)/GRID` in stationary ticks/core. For potential units divide by `S=2^40` exactly once.

The literal fixture with two actual sibling donors and one actual upstream uses `h=(38,4), m=(2,2), y=(2,1), s1=(2,1), s2=(3,1), u=(2,2)`. It rebuilds the signed original charge from actual gain, promotion, newborn and centering transfers, and checks the sum of donor payments against `[-P]_+`. With `s1` marked, the mother's fixed gain rate is `35` and its remaining exterior gain is `5/6`; both are integral on the replay denominator. The exact `s1` subreward is `18,546,192` raw ticks, so the integer numerator for the constant reward is `den·GRID·18,546,192`. The replay checks this conversion and both reward unit scales explicitly.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 10 nice -n 19 timeout 100s \
  /usr/bin/time -f 'wall=%e user=%U rss=%M exit=%x' \
  round2/.work/venv/bin/python -Bu \
  round2/occupation-quintet-gold2/diagnostic/lower-sibling/price.py --kind all --young 5
```

[Replay](compact.log): `wall=31.26 user=30.22 rss=274252 exit=0`. It checked 42 literal labelled H3 rows before the multi-donor fixture, 265,518 sibling terminal states and 207,900 upstream terminal states. The receipt SHA-256 is `c3d7d979fe82c4ac626cbdbef5f461411e966abd4e53ab35529a9cd5bb75c02a`.
