# Independent audit: stationary P4 band and unpaid reserve

**Correct in scope.** The tagged lifetime recursion prices the original D18 signed P4 charges. The `unpaid-*` objectives are disjoint from the frozen EDGE reserve and the audited `R_f/R_g` memory payment. Combining the supplied wide-band floors with the exact post-EDGE memory endpoints gives **`[0.565678089158861, 0.565770444554472]`** outward, conditional on the clean pin of that composition and the producer's wide-band Bellman receipt. The companion rechecks the full wide charge census and its sign corners, but does not repeat the wide Bellman run.

## What the band counts

Register each nonfounder apex `A` at its actual promotion into a degree-two core. The path is oriented down the core-child tree; its four census labels are `(D,k)` and the corresponding D18 labels are `(D−k−1,k,1)`. The actual row is charged once at `A` for each labelled `A→B→C→D` in the band. A positive contrast uses the two *already signed* raw D18 rewards at their respective original `λ`, evaluated on the same live trajectory; it is not a third additive reserve.

- **Gain and score:** The known arrival pressure is ordinary leaf count plus the reciprocals of adjacent path tags. Only the exterior mother at `A` and the off-path core children at `B` are omitted. Those same exterior arrivals multiply the reward's free coefficients and change the recipients' degrees. The two incoming centered transfers are from the three-core chains obtained by deleting `A` or `D`. No replacement at either end or bare-apex/`C` assumption is made.
- **Matched births:** At a tip promotion, one physical arrow at total rate `L/D` continues the old prefix with its updated arity and emits one extended prefix with a virgin last tag. At `B`, the child is a new side donor for the old prefix but not for the new extended prefix. An already registered `ABC` or `ABCD` acquires each later `B` child as a side donor. Multiplicity is additive even though all tags see the same physical arrows. Saturating a short prefix at its tip's arity cap only ends its ability to emit; it does not delete the longer tags it already emitted.
- **Hidden state and exits:** An exact side donor has gain pressure between `L+1/D_B` and `L+1/D_B+k/2`, promotion rate `L/D`, and sends at `1/d`. The mother gains at least `1/D_A`; her own core parent contributes at most `1/2` and her remaining neighbors at most `m−2`, giving the claimed `1/D_A+m−3/2` upper bound. This uses a nonfounder mother; the founder and its direct children do not contribute in the normalized stable-age limit. Hidden promotions and arrivals that leave the retained tuple unchanged have zero drift on that tuple, while subsequent rates still lie in their conditional boxes. Degree and arity cannot decrease, so a tag outside either cap cannot become eligible, and an overflow donor never returns to exact tracking.
- **Stationary reduction:** On each row the endpoint drift is affine in the optional arrival rates and the discount, with each rate used consistently in reward and continuation. The signed drift check is the direct robust inequality; directed rounding of each intercept/slope and continuation preserves it because rates are nonnegative. The projected DAG progresses by degree, arity, extension or irreversible overflow. No Markov assumption on the projected state is needed. Applying the inherited physical birth-Palm stable-age formula to each live `A` family multiplies its discounted lifetime by the true `μ`. Taking both discount-factor endpoints in the final product handles signed values; bounding over all entry-mother bins covers an arbitrary physical birth mixture. The independent wide scan finds both raw signed rewards nonnegative at all affine port corners. Dropping other charge rows or replacing the `unpaid-*` lower endpoint by zero if negative is therefore safe. The same operation is not justified for the signed contrast.

## Exactly which P4 charge is unpaid

For a full path put `j=k_B−1`. The four tags with their immediate ordinary leaves contain exactly `Σ(D−k)` vertices. The `j` additional, disjoint core-child branches each contain at least two vertices, so the full planted fringe rooted at `A` has size at least `M=Σ(D−k)+2j`, excluding the exterior mother, and at least `4+j` cores. The predicate `M>20 or (M>12 and j>0)` avoids both paid higher-family masks in the all-seven EDGE ledger: full `A` fringe size at most 12, and size 13–20 with three/four cores. This is a sufficient test; a descriptor may fail the test while a deeper or larger realization is in fact unpaid. Tests at the exact boundaries include:

| Path `(D,k)`, with minimal side branches | Actual size | Cores | Kept |
|---|---:|---:|---|
| `(3,1) (4,2) (3,1) (4,0)` | 12 | 5 | No |
| `(4,1) (5,3) (3,1) (5,0)` | 16 | 6 | Yes |
| `(8,1) (6,1) (5,1) (4,0)` | 20 | 4 | No |
| `(8,1) (6,1) (5,1) (5,0)` | 21 | 4 | Yes |

Putting the same P4 inside an ancestor's fringe does not change its owner. The frozen uniform inner-edge reserve, high-parent `g/hC` payment, and the subsequent `R_f/R_g` objectives spend other charged rows, even when they refer to the same tree, arrow, or edge. See the all-seven reserve (`occupation-neon/interpolate/compose/README.md#accounting`), [EDGE audit](neon-edge-sol.md), and [memory audit](edge-memory-sol.md). The complete band, the unpaid subset, and the nested small/large band solves are not additive. The tracked single target has `M=11`, so `single=True` has zero unpaid value even when an actual realization of that descriptor happens to grow beyond the reserve.

## Checked arithmetic and replay

Reported wide run: tagged degree caps `(12,8,6,16)`, `k_B≤3`, donor/mother caps `2/4`; 14,850 supported descriptors, 12,346 unpaid, 196,460 live states and 2,357,520 signed drift inequalities. The lower/upper unpaid reward floors are respectively `8830984420263341/(2·10^23)` and `6793831153142289/(2·10^23)`. Add each floor to its own signed stationary numerator, carrying forward the post-EDGE memory result once. The potential is unchanged, so use the existing ceilings `1357493815758092773/10^18` and `273436004214812983/(2·10^17)`. Equivalently, move the *exact* lower/upper memory endpoints inward by the new floor divided by its ceiling, about `3.2526794e−8` and `2.4846147e−8`. Exact endpoints before rounding:

```text
137070954632127182140563254618897/242312646112819559980500000000000
6574835410627056583269136805857/11621030179129551777500000000000
```

From `[local path]`:

```sh
OPENBLAS_NUM_THREADS=1 nice -n 15 taskset -c 3 timeout 80s \
  round2/.venv/bin/python -Bu round2/checks/p4-band-sol.py \
  b3f18_lo.bin b3f18_hi.bin --wide-rows
```

The companion (`p4-band-sol.py`) authenticated both full payloads and five source pins, including `path_band.py` at `3bb442ee42b93fa431ae00e8b9cbf0128fbb013c6f2cdca143a6ab49d3c2f642`. It routes raw labelled arrows through full eight-family potential differences and rebuilds each P4's transfers from the two induced three-core chains, separately checking the total P4 allocation and the total whole-tree drift on each witness. It passed 48 signed P4 comparisons and 176 nonzero omitted-term controls. It reuses the prior labelled projection oracle for 192 projected rows and 1,396 merged successors, and checks ten arrows from real apex birth through nonbare `A/C`, late side births and irreversible exit. The wider static scan reproduced both reported descriptor counts and passed 118,800 exact sign-corner checks.

With `single=True`, donor/mother `2/3`, original tagged degree and arity caps, the same discount bracket, and `θ=1.2415`, the four common objectives agreed exactly when the reference grid was aligned to `10^18`; with the original `10^16` grid all endpoint differences were below `10^−15`. The replay printed 255 live states and 3,060 signed drift checks for that comparison. Its contrast straddles zero under this particular small-cap mother box; the stronger single-cell sign under other caps or constrained Palm weights remains a separate result.
