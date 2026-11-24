> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Gold2 stationary frontier

This composes the original two Gold2 witnesses with the D24 host envelopes,
the coupled singleton replacement and one disjoint set of edge and H3 floors
per sign. The corrected base receipt gives
**`0.565714262991 < μ < 0.565735974375`**. The independently checked
live-sibling revision gives **`0.565714467904 < μ < 0.565735974375`**.
All ledger values in the receipts are exact ticks/core, with `2^40` ticks
per potential unit. `241/426` lies strictly inside both intervals.
The original production bracket was `[0.565713486333, 0.565736986679]`
(outward). The earlier H3-derived stationary interval used invalid credits
and is not a benchmark for improvement.

## Physical budget

Each mask below is evaluated at payment time. `D` is core degree and `c`
includes the mother. A payer and a donor are distinct physical neighbors.

| Slot | Lower | Upper | Frozen sources |
|---|---|---|---|
| Whole `g` row, D18 | Both endpoints `D≤18` | Same | [Lower](../../edges/edge18.log), [upper](../../edges/upper/ceiling8.log) (its D18 line) |
| `g`, small-payer FY slot | Host `D=19…251,c≤17`; payer `D≤5`; any distinct donor | Host `D=19…251,c≤22`; payer `D≤5`; any distinct donor | [Lower](../../edges/absorb12.log), [upper](../../edges/upper/absorb15.log) |
| `g`, wide-payer FY slot | Host `D≤150,c≤10`; payer `D=6…34,c≤7`; donor `D=2…5` | Host `D≤251,c≤12`; payer `D=6…35,c≤6`; donor `D=2…5` | [Lower](../../payer/README.md), [upper](../../payer/upper/README.md) |
| Whole `g` row, double tail | — | Parent was `D≥252` at actual child birth; pay once per extant edge from the child's first `D=252` | [True drift](../../true-drift/README.md), exact transform in its `check.py` |
| H3, threshold | Mother `D=19…79,c≤8`; middle `D≤5,c≤4`; child `D≤18,c≤3` | Same | [Lower audit](../../../checks/gold2-lower-h3-alt-sol.md), [upper repair](../../diagnostic/upper-correction/threshold-upper.log) |
| H3, lower seam | Mother `D=80…100,c≤8`; same middle; child `D≤5,c≤3` | — | [Lower audit](../../../checks/gold2-lower-h3-alt-sol.md) |
| H3, high mother and middle `(2,2)` | — | Mother `D=80…480,c≤12`; middle `D=c=2`; child `D≤12,c≤3` | [Mother-high replay](replay-upper-h3.py) |
| H3, high mother and middle `3…5` | — | Mother `D=80…120,c≤10`; middle `D=3…5,c≤4`; child `D≤7,c≤3` | [Corrected increment replay](replay-upper-increments.py) |
| H3, higher mother stock | — | Mother `D=19…79,c=9…12`; middle `D≤5,c≤4`; child `D≤7,c≤3` | [Corrected increment replay](replay-upper-increments.py) |

Both FY slots require an exact G2 host outside M3 and a paying edge with
maximum endpoint degree `>18`. The payer cut at 5 separates the two FY
rewards; D18 spends neither. The double tail is already beyond exact G2 at
both endpoints. Its zero original charged edge row leaves the entire edge cap
on the specified history subset. This is the only payment retained from the
true-drift edge recipe; the other three prices there are alternatives to the
chosen upper edge floors. The signed `3999→4000` entry-flux bound is not a
production credit.

The H3 rows follow actual nonfounder mother→middle→child paths. Their four
upper masks separate first by mother degree at 80, then by mother arity at 9
below 80 or by middle degree at 3 above it. Each pays `[−q]_+` on the original
signed charged H3 row, preserving its original transfers. Thus `q+[−q]_+`
is `max(q,0)`, still covered by the original nonnegative cap at the **middle**
(family cap 408/233). D24's host envelope retains these caps, so its savings
do not spend the H3 payment a second time.

The historical threshold and upper increment amounts are superseded. Their
affine solvers divided the full exterior reward change by a rate width scaled
by the common denominator `d`, omitting the corresponding `d` in the reward.
For physical width `R` and outward row slope `b`, the scaled rate width is
`d·R`; its continuation shift must be the directed rounding of
`d·(−b·R)/(d·R)=−b`. The separate cohort replay instead flipped the sign of
the rounded slope. The corrections use the frozen source in memory and keep
the original production certificate and host envelopes intact.

### Separate lower live-sibling calculation

On the lower witness, replace the corrected `19…79` threshold term with
two payments on the narrower mask: nonfounder mother `D=19…79,c≤8`, middle
`D≤5,c≤4`, and child `D≤5,c≤3`. Track a real maternal sibling or the unique
nonfounder upstream vertex through degree three and retain an absorbing
reciprocal band thereafter. The three chronological sibling registrations
count distinct physical siblings once. The lower seam at mother `D=80…100`
remains disjoint. The new [producer log](../../diagnostic/lower-sibling/compact.log)
records exact floors of about `251635.821901` ticks/core for siblings and
`26800.032788` for upstream. No box floor is booked.

For an H3 occurrence let `B` be its outward all-port maximum and `b` the
mother's outward receive slope. Each distinct tracked exterior donor `i`
recovers `z_i=max(0,b/2)−b/D_i≥0`. For `f(z)=[z−B]_+−[−B]_+`, convexity
and `f(0)=0` imply `Σ f(z_i)≤f(Σ z_i)`. Since the physical outward row
`q≤B−Σ z_i`, the summed payments are at most `[−q]_+`, retaining the
original cap on every occurrence. This argument allows joint payment of
sibling and upstream without paying for the same port twice.

One original singleton `f` reward is placed on each core, with newborn `g`
still charged to `f`, newborn H3 to `g`, and centered receive transfers
unchanged. The marked lifetime floor (mother degree through five, two newest
donors through degree three) **replaces** the previous singleton floor. Its
value is placed once on each actual birth at the matching mother mark;
`(0,0)` is the conservative overflow mark. The D24 host tables retain their
original full degree and arity overflow. They bound the original potential
and higher caps independently of the singleton replacement.

## Exact comparison

Let `d=+1/−1` for lower/upper, `E` be the original error, `K` the sum of
the disjoint edge and H3 floors, `X=πΦ>0` the original potential mean,
and `H−h` the D24 host savings. The original budget gives

```
d(μ−λ) X ≥ K−E + πF + π(H−h).
```

The scalar control uses the same singleton, edge and H3 floors and original
`M_max`, without host savings. For a trial `s` with `a=d(s−λ)≥0`, the common
LP prices the signed objective `b·V(mark)+H−h−a·u`, where `b` is the
actual promotion rate and `u` is the D24 upper comparison for `Φ`. The
birth-degree tags for `D=2…12` impose both bounds on actual marked birth flux; each
nonnegative tag slack has its own mass bound `1/LO` in the exact repair.
Every LP column is repriced against the rationally rounded dual. If `z` is
the resulting lower price in `π/μ` units, the valid price in `π` units is
`min(LO·z,HI·z)` even when `z<0`.

The uniform margin `m` gives `d(μ−s)X≥m`; when `m>0`, use the D24 upper
bound for `X` once to turn it into an inward endpoint. A restricted trial
uses `[LO,s]` below or `[s,HI]` above, only under the contrary hypothesis,
including equality. A positive margin excludes that hypothesis; a negative
margin does not prove a sign or a ceiling on a different proof. Restricted
margins are never divided to offset the trial.

| Receipt | Scalar control, approx. | Restricted `241/426` margins, ticks/core (lower / upper) |
|---|---|---|
| Corrected base | `[0.5657141084694073, 0.5657362035344400]` | `−18264451.378257 / −11254455.943124` |
| Lower live-sibling replacement | `[0.5657142624199445, 0.5657362035344400]` | `−17986015.945067 / −11254455.943124` |

The base's restricted `99/175` margin is `−30858.044834` ticks/core. The
live-sibling calculation changes that margin by exactly the net change in
the lower H3 budget: it uses the same roots, cells, prior and objective.

## Replay

From the Astra root, using the existing SciPy environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 15 nice -n 19 \
  round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/host/frontier/check.py
```

The default pins and consumes the original-witness singleton root exports
and rebuilds the rational repairs. `--replay-singleton` recomputes those
roots from the same read-only inodes it SHA-authenticates. No Claim 10
marginals from the root-export experiment enter this proof.

The base and live-sibling replay hashes are
`9f9402eb89f711e5c12729a3fcfe50e31ee2e3ad9ee3e862a9dcfd55dbf3b03b`
and `7a6ddda08b83574f6b745f11b3e5445ba80a60f47d7ef75c625b093a56cad20f`.

Add `--live-sibling` to replay `receipt-live-sibling.json`; it binds the
unchanged upper result to the exact base receipt by SHA-256. It gives
`+247577.388356574` ticks/core at the restricted `99/175` trial. The
[cold review](../../../../../../../reviews/cold/astra/checks/gold2-frontier-sol.md) independently
reproduced both compact payments at every Bellman corner, the revised
lower trial, and the unchanged upper endpoint.
