> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# D20 signed singleton and existing-edge corrector

**Producer result, pending cold review:** `[0.565681980098700, 0.565766097351300]`.
Compared with the reproduced D20 renewal bracket
`[0.565681716887067, 0.565766354357632]`, the inward moves are
**`2.63211632e−7` lower and `2.57006331e−7` upper**, before outward
display rounding. The `1e−6`-per-side target was not reached. The entire
D20 size-16 joint payment is retained on both sides.

| Same size-16 reserve and edge cut | Lower endpoint | Upper endpoint |
|---|---:|---:|
| D20 renewal | `.565681716887067` | `.565766354357632` |
| Edge + singleton birth-Palm, divide by D20 ceiling | `.565681717514894` | `.565766353731946` |
| Signed denominator-coupled birth-Palm | `.565681979741300` | `.565766097764600` |
| Signed stationary incidence LP | `.565681955022224` | `.565766123513244` |
| Signed pointwise corrector + reserve | **`.565681980098700`** | **`.565766097351300`** |

The edge-only increment over renewal is about `6.3e−10` per side; the
larger improvement comes from jointly pricing the D20 potential deficit
and signed singleton charge. Using only the all-stock caps at each exact
census can give a result worse than renewal. The `--exact-den` computation
uses the D20 `Charge` interface for the exact census envelope, checks its
all-stock caps against the pinned export, and inherits the exported
all-age ceilings `1.356756788233927 / 1.365141582117995` (outward display).
No D18 numeric inputs or duplicate GMP verifier are used.

## Exact rows and newborn compensation

Write `x=(D,k)`, `J=(2,0)`, `x⁺=(D+1,k)`, `xᵖ=(D,k+1)`, `b=L/D`.
For a nonfounder edge `p→w`, track both censuses until either exceeds
`D=16` or `k=5` in the principal run. Each endpoint receives the other
endpoint's exact reciprocal; its other incoming pressure is in a
physical range `[0,k_parent/2]` or `[0,k_child/2]`. A single pressure
variable enters both the signed charge and the corresponding physical
gain rate. Compute `E` by downward Bellman comparison:

```
(s − A_existing) E ≤ eligible D20 g.
```

For relaxed affine boxes the code uses the positive part and checks
the corners plus the affine zero/box intersections; it agrees with
the nonnegative production charge on physical inputs. Degree and child
count are monotone at the stopped boundary. Existing-edge propagation
contains both endpoint promotions but no newborn source.

The singleton table then solves

```
(s − A_single) W ≤ f + b E(xᵖ,J) + δ (M − Q(x)).
```

Here `δ=|s−λ|` at the pointwise trial (or the trial price used for the
stationary objective), and `Q(x)` is the exact D20 census potential
majorant including the export's zero-sum recharge. `M` is its all-stock
cap. The same total incoming pressure changes the singleton reward and
the gain rate. Beyond the finite census export use the exact D20 stock
floors and degree-independent recharged stock caps, with overflow at
the last bounded arity class. The overflow recharge cannot increase.

The compensation is for a **real promotion**: the old core stays with
its post-promotion `k+1`, and the newborn starts at `J`. The new edge
value belongs to the promoting singleton at precisely this one event.
Singleton birth has zero corrector value after recentering `W` by
`W(J)`. Neither operation replaces the newborn component already
assigned to the production `f` charge. Edges from the founder are
omitted, as are founder singleton terms; the potential normalization
and drift inequality in the core limit discard this boundary. Positivity
of the corrected signed potential follows from the original D20 lower
potential bound and the checked finite correction norm.

For the pointwise test add `Σ(W−W(J))+ΣE` to the **signed** production
potential. The normalized at-trial lower drift is at least
`s W(J)−δM` plus the untouched nonnegative reserve. The two exact trial
offsets are `19800987/10^13` and `19026487/10^13`. The combined margins
are strictly positive. The stationary birth-Palm version uses the
D20 bracket's upper discount and lower prefactor; the LP version uses
stationary identities directly. These are alternative certificates.

## Reserve accounting and stationary comparison

In `preserve` mode charge `g` only when the two endpoint censuses imply
planted parent fringe size `>16` or at least five cores. The lower size
count is `2+L_parent+L_child+2(k_parent−1+k_child)`, and the core count
is at least `1+k_parent+k_child`. A zero-paid edge may have positive
continuation to eligible rows; this is not immediate payment. The full
size-16 `g`/higher reserve is replayed unchanged (downward display
`5.7189520e−8 / 2.4708513e−8`). In `replace` mode the full size-16 joint
price is discarded, the six higher families are repriced together, and
the edge corrector may pay all edges. Upper cap-13 and lower cap-16
pointwise comparisons showed only small advantages from replacement;
the principal result uses `preserve`. No mode adds separate prices to
the full joint reserve.

The LP tracks exact censuses to degree 40, an aggregate infinite-degree
tail, true promotion sources and one shared parent-to-child incidence
table with counts and reciprocal moments at both ends. It prices
`f + bE(xᵖ,J) + δ(M−Q)` jointly. The genuine physical stationary law
satisfies the age rows `π[(μ−A)h]=μh(J)`, `πb=μ`, and total receive
`πρ=1−μ`. Tail arity and receive overflow have explicit nonnegative
columns. The common `μ` is relaxed to convex endpoint lifts; the true
law is feasible. Every objective coefficient for a grouped tail is
bounded below over its entire degree/arity cell. In each of the four
column classes (occurrence, incidence, excess stubs, excess receive),
the total physical mass is at most one, by probability normalization,
one oriented edge per nonfounder, `πk=1`, and `πρ≤1`. Subtracting the
largest positive exact reduced-cost defect of each class therefore
certifies the rational dual proposed by HiGHS. Final repair was below
`4.3e−14` lower and `3.6e−14` upper in the objective's original units.

At renewal plus `1e−6` the signed LP certificates give margins
`−9.55318678e−7` and `−9.70422493e−7`. These negative **lower bounds**
quantify the shortfall of this certificate; they are not upper bounds
on the physical result or a proof that no stronger correction works.

## Replay and evidence

Run from `.` for `S=lower` and `S=upper`:

```sh
taskset -c 6 nice -n 19 timeout 200s /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M exit=%x' \
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 round2/.work/venv/bin/python -Bu \
  round2/occupation-neon20/corrector/replay.py --side "$S" --mode preserve \
  --edge-degree 16 --arity 5 --stationary 40 --fixtures --exact-den
```

The authenticated lower/upper stationary runs exited zero (`153.38` /
`101.68` wall seconds; `82.74` / `68.72` CPU seconds). Each had 5,200
edge rows, 2,079 authenticated singleton rows, and a 1,109-row /
11,394-column LP. Later guard-only verification reruns (same flags with
`--stationary 40` omitted) reproduced both exact pointwise endpoints
and checked 81 full labelled D20 production drifts per sign, 114
nonzero newborn edge sources, and 77 fixtures rejecting omission of
those birth terms. Both exited zero (`54.52` / `43.02` wall seconds).
The full payload digests and all four D20 row/source hashes are checked
on each run. No external cold review is claimed in this producer pass.
