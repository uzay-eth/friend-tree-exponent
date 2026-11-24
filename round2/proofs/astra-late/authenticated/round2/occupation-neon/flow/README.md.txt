# Joint reciprocal-edge flow for the D18 occupation dual

**Result.** Pricing the all-degree first-body numerators of the two pinned Neon D18 potentials with one shared endpoint measure improves both normalized census bounds. With census cutoff 48 and the existing [all-degree denominator ceilings](../denom/README.md), the resulting interval is **`[0.565674991393, 0.565773407322]`** from first-body rewards alone, or **`[0.565675932482, 0.565772309419]`** with the already certified, disjoint higher-motif rewards. The latter interval moves the [previous size-12 higher-motif bracket](../higher/README.md) inward by `3.34824e−7` and `3.68450e−7`. This is an occupation improvement for the original potentials, not a new GMP certificate or a change to the denominator proof.

All prices in the next table are directed **lower** bounds, in units of `10⁻⁶` per core. The comparison column uses the same truncated census, exact normalizations and scalar reciprocal budget, with only the joint endpoint constraints removed.

| Census cut | D18 side | Previous all-time Bellman | Normalized census, scalar flow | Plus joint endpoint flow |
|---:|---|---:|---:|---:|
| 36 | Lower | 2.512846393249 | 2.627859376212 | **2.651316728138** |
| 36 | Upper | 2.025752646833 | 2.112435943813 | **2.137203454752** |
| 48 | Lower | 2.512846393249 | 2.681149834029 | **2.701989116205** |
| 48 | Upper | 2.025752646833 | 2.157431269536 | **2.177461932639** |

## 1. Exact stationarity and edge reversal

Under the oriented stationary uniform-core law `π`, a typical core has one parent, `k` core children, `L` leaf neighbours, `D=L+k+1`, `b=L/D`, and receives core-to-core arrows at rate `ρ=Σ_core-neighbours 1/D_y`. Here mass transport is taken on the core tree, with the ancestral orientation retained. For bounded `f` (or nonnegative `f` by truncation),

```text
π(k)=1,                 π(b)=μ,
π(f(x) ρ(x)) = π(Σ_(y core adjacent to x) f(y)/D_x),
π(ρ)=π((k+1)/D)=1−μ.
```

For the first and third identities, reverse the finite-tree sums exactly, then use the [oriented mean limit](../../fringe/spatial-limit.md). The founder is a single vanishing fraction for vertex tests; its terms with an edge multiplicity are controlled as follows. There is uniform integrability of the number of core children in mean-rooted trees: for any `θ>0` with `a=e^θ−1<μ`, the generator on the sum of `e^(θk)` has drift at most `a` times that sum plus the mean core-birth rate, itself at most `E C(t)`. Each promotion increases one old `k` by one and creates a new `k=0` core. The resulting bound is `O(e^(μt))` by the mean estimate and `a<μ`; dividing by `E C(t)≍e^(μt)` gives a uniform exponential moment on large times. The same applies to incoming tests dominated by `(k+1)/2`. Thus reversal survives the limit without losing a high-arity edge mass. The birth balance `π(b)=μ` is the [physical cell-generator / stable-age identity](../renewal/README.md), and the last equality follows by setting `f=1`. None of these identities asserts `ρ=(k+1)/D` at a vertex or on a degree class.

Orient each physical edge from parent to child and let `E` be its joint endpoint measure per core. On an edge let `p=1/D_parent`, `q=1/D_child`. Its row count measure is `π(k·1_cell)`, and its column count is `π(1_cell)`. Its weighted row measure for `p` is `π(k/D·1_cell)=π((1-b-1/D)·1_cell)` for cells supporting outgoing edges. Its weighted column measure for `q` is `π(1/D·1_cell)`. Mass reversal with `f=1_cell` is the additional equation

```text
π(ρ·1_cell) = Σ_(edge parent in cell) E[q] + Σ_(edge child in cell) E[p].
```

The same edge column contains **both** endpoint reciprocals. In particular a high-arity state's incoming rate may exceed one, but it draws from the same child-stub supply that must fund the receivers' down-traffic elsewhere. Four independent marginal rates on an edge cannot enforce this. Nor is the law of a parent under edge sampling `π`; it is biased by `k`.

## 2. The projected birth-occupation LP

The physical stable-age occupation of the tracked root uses discount/birth rate `μ` and starts from `(D,k)=(2,0)`. At a real state it gains at `g=L+ρ` and promotes at `b`. For a census test `V`, the exact linear occupation rows are

```text
π[ μ(V−V(2,0)) + g(V−V(D+1,k)) + b(V−V(D,k+1)) ] = 0,
π(1)=π(k)=1,      π(b−μ)=0,      π(ρ+μ−1)=0.
```

These follow from the physical stable-age formula and discounted Dynkin, or from normalized birth occupation. They remain valid under history-dependent incoming rates; there is no transition law assumption on the two-coordinate census. For a truncated test, a gain out of the retained census increases its state to the overflow bin in the same column as the rate `ρ`. Integrability follows from bounded tests, the finite cutoff for degree-changing transitions and the core-child moment above. The LP uses `μ∈[.565673570322,.565775]` at the two endpoints in each occupation column and lets those columns mix independently. A real occupation lifts with the **same** convex coefficients for `μ` everywhere, so this further relaxation is outward. At retained states, the reward and all constraints depend affinely on `ρ∈[0,(k+1)/2]`; they use the same column interpolation. The LP does not minimize the price at one `ρ` and propagate the exit at another.

For the endpoint constraints group degree into bins ending at `2,3,4,5,6,8,11,14,17,21,27,36,48` (as applicable) and one unbounded bin, and arity into `0`, `1..3`, `≥4`. For each possible ordered pair of cells, make columns at the corners of their reciprocal intervals. Every actual edge lifts to a convex combination of these columns preserving both reciprocal coordinates. In particular the unbounded degree bin has reciprocal interval `[0,1/(cutoff+1)]`, not a fixed reciprocal. The edge measure and occupation measure both normalize to one, but are distinct measures.

The census overflow is **not killed**: above the cutoff keep `k=0..11` and one saturated `k≥12` bin. At fixed overflow arity the degree moves no longer change the test, and `b=1−(k+1)/D` stays coupled with `p=1/D` in that column. In the saturated bin both gains and promotions leave the test unchanged. Keep `b,p≥0`, `p≤1/(cutoff+1)`, `b+13p≤1`, the exact aggregate `kp=1-b-p`, an extra nonnegative stub column for `k−12`, and a separate nonnegative column for its whole incoming `ρ`. The latter columns have total masses at most one, by `π(k)=1` and `π(ρ)=1−μ<1`, respectively. We do not require the relaxed saturated row to recover the individual `k` of each endpoint; we retain both weighted and unweighted totals for its bin. This avoids capping real multi-arm traffic or pretending `k=12` throughout the tail.

The objective is the accepted [charged first-body price](../rows/README.md): exact signed affine values for each of the 2,079 nonfounder descriptors through degree 64 per side, and the certified all-stock floors for `k≤4` above that. In the aggregated degree tail use, per `k`, the minimum of the remaining exact row corners and the corresponding all-stock floor, zero where there is no floor. This is a conservative contribution on the whole degree domain. The old `tracked` baseline uses the extended exact rows through 64 and the same all-stock floors. Other nonnegative GMP motif rewards are omitted from the dual. The optional final interval adds a downward-rounded floor for just the six separately priced motifs: recover that floor by multiplying the published directed endpoint improvement by the exact denominator used in the [size-12 replay source](../higher/extract.py), then round down to `10⁻¹¹`. No first-body or edge surplus is counted twice. Dividing by the existing `1.356833021311` and `1.367169995210` denominator ceilings and rounding the two exponent endpoints outward is the usual signed-drift transfer.

## 3. Exact dual verification and exit diagnostics

`check.py` uses HiGHS only to propose the equality-row multipliers. Round them to rational multiples of `10⁻¹¹`. Independently check every exact column and collect the largest **positive** reduced-cost defect separately for occupation, edge, overflow-stub and overflow-incoming columns. Subtract these from the displayed dual objective (each class has total mass at most one). Equivalently, add one slack inequality per class and use the corresponding negative defect as its dual multiplier. The checker asserts every repaired column inequality with `Fraction`; the certified bound does not use the floating primal or dual objective. The tested maximum **sum** of repairs for the two cutoff-48 prices is below `4×10⁻⁸` in the table's micro-units. There is no total-unimodularity claim about the coefficient matrix: the exact equalities here are from unimodular mass transport, and their columns involve reciprocals.

At census cutoff 36, direct objectives `−ρ·1_exit` give the following directed **upper** bounds on the receive part of the gain-exit flux per core. The degree-36 rows are the actual last retained rows, so these fluxes drive gains into the unbounded census tail. They are unconditional integrated fluxes, not pointwise or conditionally uniform bounds on `ρ` at an individual core.

| Receiver class | Scalar flow | Joint endpoint flow |
|---|---:|---:|
| `18≤D≤36, k≥4` | `.023811781865` | **`.023562301124`** |
| `D=36, k≥4` | `.001041254316` | **`.001032995717`** |
| `D=36, k≥8` | `.000022185491` | **`.000021468472`** |

The reciprocal gain is bounded jointly where the exit contribution is charged; replacing the pointwise Bellman envelope by any of these flux numbers would be invalid. The static LP still relaxes occupations beyond the cutoff and has no tracking of the time correlations of particular arms.

## 4. Exact bad-model controls and bounded replay

On a six-arm top of degree eight with three children each at degrees two and four, the upward total is `9/4` and the downward total `3/4`. Filling all six upward slots with degree-two donors is impossible with those child margins. Testing just scalar reversal misses the wrong endpoint allocation; the fixture also tests the weighted identity restricted to several receiver degree classes. On a 15-arm core of degree 32 with parent degree 16 and five degree-two, ten degree-eight children, actual receive total is `61/16`, not `15/2+1/16`.

Even two individually margin-feasible endpoint pairings cannot be mixed by direction. With parents of degrees 16 and 32 and children of degrees 2 and 8, compare the two bijections between slots. Maximize upward contribution to the degree-32 parent and downward contribution to the degree-two child. Separately the maxima sum to `1/2+1/16=9/16`, while the maximum from a **single** joint edge measure is `17/32`. Convex mixtures cannot exceed that joint maximum. Thus keeping the same margins in two edge tables is not equivalent to a shared joint table.

From repository root (SciPy is available in the existing interpreter):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 45s taskset -c 7 nice -n 19 /usr/bin/time -f 'wall=%e rss=%M exit=%x' round2/.work/venv/bin/python -Bu round2/occupation-neon/flow/check.py --diagnostic
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 35s taskset -c 13 nice -n 19 /usr/bin/time -f 'wall=%e rss=%M exit=%x' round2/.work/venv/bin/python -Bu round2/occupation-neon/flow/check.py --limit 48
```

Both passed: cutoff 36 used 837 rows and 5,106 joint columns (40.82 s wall, 127,556 KiB); cutoff 48 used 1,360 rows and 7,644 joint columns (18.80 s wall, 155,996 KiB). Earlier cutoff-64 searches hit their bounded time limits without a certificate. Only the sparse charged-row package is read, with its pinned manifest and exact operand-to-corner checks; the raw multi-gigabyte payloads and full GMP inequalities are prior certified inputs and are not reread here.
