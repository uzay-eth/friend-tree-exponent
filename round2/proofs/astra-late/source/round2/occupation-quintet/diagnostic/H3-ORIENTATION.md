> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Original `b04cc50` H3 orientation diagnostic

Replay: `bash round2/occupation-quintet/diagnostic/h3_replay.sh`. The two witness passes and the report completed with exit 0. For all nine preserved age-20 trees, H3 counts match the prior physical census exactly; signed total, negative full-row reward and uncovered gaps agree within `2e-11` in normalized total units. Eight pooled trees contain 896,887 cores; the held-out tree contains 235,709. No new trees were generated.

Numbers below are lower/upper witness pairs, in `10^-6` normalized potential units per core. The reward is `max(0, -charged H3 row)`, clipped after the full row. A port's uncovered box gap overlaps that row; it is not an additional credit.

## Where the negative rows are

With middle arity `c=2–4`, high tip degree `19–1998`, and founders and the old masks excluded, the pooled negative reward splits as **4.958/6.235** (high tip is the mother), **0.576/0.605** (high tip is a child, low tip is the mother) and **0.223/0.281** (two children). Outside this priority slice, the full H3 census includes **0.824/0.721** in tied tips, **0.077/0.083** in the old strict chain, and **0.008/0.006** in the high-middle mask; the high-middle mask has no geometric overlap with the priority slice.

A compact single-orientation cut clears 2 on both witnesses and both cohorts: **mother degree 19–65, middle degree ≤5 and arity 2–4, child degree ≤18**, with founders and the old strict chain excluded. This covers 44,884 pooled rows and 11,986 held-out rows. Reward is **2.075/2.205 pooled** and **2.058/2.163 held out**. Degree caps of 100 at the mother and 11 at the child, with the same middle and exclusions, give **2.219/2.341 pooled** and **2.212/2.303 held out**.

Within the 65/18/5 cut, exact middle `(D,c)=(2,2)` supplies **1.610/1.702**; exact child `(D,c)=(2,1)` supplies only **0.500/0.617**. These stock slices can overlap and neither measures elapsed age. The leading lower-witness arity triples `(mother,middle,child)` are `(4,2,1)`, `(5,2,1)`, `(3,2,1)`, `(4,2,2)`, `(5,2,2)`, and `(6,2,1)`. Mother degree bands 19–24, 25–37 and 38–65 contribute **0.753/0.758**, **0.598/0.623**, and **0.724/0.824** respectively. In the first two bands the uncovered gap is numerically almost the full negative reward; in the third it is only **0.312/0.346**.

## Ports and seam

The 65/18/5 uncovered gap is concentrated at the high mother's receiving port: her other children contribute **1.061/1.139**, her own mother **0.545/0.458**, and the low child's exterior children **0.057/0.129**. There is no uncovered middle receiver. The actual arrows at the high mother's receiver are much smaller than those gaps: signed contributions from her own mother are **+0.002/−0.013**, and from her other children **+0.112/−0.209**. At high-mother arities 3–8 the upstream degree exceeds 65 in 49–68% of rows; an upstream degree above 65 and an other-child reciprocal sum at least `(c−2)/4` co-occur in 30–40%. These are occurrence frequencies, not bounds on rates or birth-flux measurements.

For the clean mixed slice with low tip degree ≤11 and no middle cut, exact high degree 1999 occurs 5 times pooled and once held out; exact 2000 occurs only once, held out. All are mother-high. Their held-out negative rewards are `0.00000307/0.00002222` at 1999 and `0.00000253/0.00000009` at 2000, in the same per-million units. On the lower witness the two observed incoming tip-transfer corrections into the middle are `+0.04550,+0.01367` at 1999 and `−0.01021,−0.00496` at 2000; these are different sampled rows. The neighboring 2001–2011 band has 95 pooled rows (90 mother-high, five two-child) and nine held-out (all mother-high). The seam counts are too sparse to infer an occupation or transition law.

`h3-report.log` contains the stock, cap, arity, and signed port tables; `h3-run.log` and `h3-report-time.log` record the replay timings. These are diagnostics on the original witnesses. They do not certify a stationary floor, change an endpoint, or establish disjoint ownership with other payments.
