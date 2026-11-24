> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# D18 cross-ratio: paired stopping and a realizable obstruction

The lower and upper D18 witnesses have a **signed killed-prefix comparison** from an exact paired stop. A pointwise ratio comparison on all physical trees, however, cannot improve either endpoint of the incoming `[.565674722414,.565773666188]` bracket. Two reachable, repeatable families refute it even with the founder effect removed, every motif included, and all sends evaluated at their actual physical rates.

## Exact whole-tree comparison

Let `a=.565673`, `b=.565775`, `δ−=ℒΦ−−aΦ−`, `δ+=bΦ+−ℒΦ+`, and `Dθ=δ−−θδ+`. The cross-ratio uses limiting expected charge per core (the same physical normalization for both witnesses):

```
(mu−a)/(b−mu) = (πδ−/πδ+)(πΦ+/πΦ−).
```

The incoming bracket corresponds to odds `0.0171764605502 … 75.4725463559`. The identity can also be obtained by normalized Laplace averaging and the common first-moment limit; it does not require treating a pointwise ratio as an observable with that mean.

`replay.py` contracts the two pinned D18 tables on the **same** physical Bernstein embedding and enumerates `f,g,hV,hC,P4,Yt,Yb,S3` with the physical multiplicity. At a receiver, it sums all gain increments before multiplying by the actual `L+Σ1/D_donor`. At a promotion it counts each pre-existing motif's change plus the new occurrences at `L/D`; the appended core is a single virgin. No signed corner, common-parameter relaxation, separate absolute error, or unphysical independent donor rates are used in these whole-tree computations. For 19 fixtures (including tail/cap boundaries), both six-cell obstruction witnesses and two actual parent ports, a separate full-potential-before/after evaluation gives identical rational drift. It uses the reviewed D18 payloads (`round2/attacks/neon1/outputs/production-final/clone/round2/proof-b3f-motif4.md`, historical reference); this replay does not hash the multigigabyte files again.

### Bulk counterexamples, without a founder artifact

Start from a spine with a founder holding one leaf, `m` interior cells, and a terminal core holding one leaf. In one train every interior cell is bare and has just the spine child (degree 2). In the other every interior cell has one leaf, the spine child, and one bare hub child with 16 core children. The hub has degree 17, and each of its 16 children holds one leaf and has degree 2. These geometries are reachable from `T3`: grow the target as an ordinary finite rooted tree by attaching vertices in ancestor order, starting with the founder's leaf and its first spine vertex. Every prescribed attachment has positive rate, even when its target is a leaf that must promote. For example, each nonfounder hub can be promoted when its spine parent has degree 4; a legal attached permanent mark assignment requires no independent reciprocal choices. D18 itself is mark-blind and reads the live degrees.

For `m≥5`, all four totals `δ−,δ+,Φ−,Φ+` are affine in `m`. An occurrence spans at most three spine edges. In its drift the one-step send or newborn extension depends on at most one additional spine edge. Equivalently, for a gain stencil sum over each connected occurrence with one adjacent donor, and for a newborn stencil remove the newborn from its occurrence to see the old connected stencil. Their spine span is at most four. Counts of interior translates therefore grow exactly linearly; counts touching either end are constant once the ends cannot interact through one of these stencils. The replay checks the same increments from five to six and from six to seven. The finite and limiting ratios are:

| Interior cell | Six-cell normalized odds | Bulk slack ratio | Bulk normalized odds |
|---|---:|---:|---:|
| Bare spine | 0.000490856735040 | 0.000161834168… | 0.000162716708495 |
| Spine with bare 16-child hub | 130.626641394 | 121.722… | 125.714454236 |

Normalized odds here means `(δ−Φ+)/(δ+Φ−)` **on the specified tree**, rather than a π average. The exact `δ−,δ+,Φ−,Φ+` per inserted cell are printed by the replay. These examples defeat any improvement based just on universal extrema of the uncorrected complete-tree ratio (even normalized by the same tree's potential ratio), or a comparison of uncorrected slack with a uniformly negligible per-core boundary remainder. They do not refute a physical occupation bound or a comparison allowing a generator correction.

## Coupled stop: keep the actual exterior rate shared

Allocate a charged occurrence to its topmost vertex, including `f` and `g`. The [charged row extractor](../higher/extract.py) is used only in the all-inner zone with current degrees at most 16; transition targets still lie in the inner range. The entire descendant crown rooted at the focal vertex is observed. Every free donor below the crown root is paid at its real reciprocal degree. All rows with a free donor above the root use the **same** `p=1/d` for the actual parent; other charges that involve that parent are owned by their own apex, not copied here. The replay cross-checks the sum of all apex charges against whole-tree drift on an all-inner fixture with all eight classes. It also verifies complete trees for two physical parent ports `d=2,8`. At `θ=6/5` one crown's pooled signed minimum improves over minimizing each signed row separately by exactly `243953811/4000000000000 = .00006098845275` (this is an instantaneous charge difference, not an exponent improvement).

At crown size at most `N` the paired reward is `q(z,p)=q₀(z)+p q₁(z)`, formed by subtracting the two exact rows before taking extrema. The crown's real intrinsic jump rates are `r(z,y)` with sum `r₀`; the parent-to-focal send goes to `Gz`, the focal leaf-gain state. In a tracked bin `2≤d≤6`, the parent gain goes to `(z,d+1)` with actual predictable total rate `t∈[1/D_focal,d]`. The lower bound is the focal's own send to its parent; the upper bound follows because its `d` neighbours each send toward it at rate at most one. Parent gains and crown jumps have distinct receivers. A parent promotion keeps `d` fixed, and is an unobserved self-transition. At `d≥7` collapse the parent bin, using `p∈[0,1/7]`; `d` is nondecreasing. A fully free comparison uses `p∈[0,1/2]` throughout and no tracked parent jumps.

For common discount `s∈[.565674722414,.565773666188]`, use the lower and upper Bellman inequalities on

```
[q₀ + p q₁ + Σ_y r(z,y)V(y,d) + p V(Gz,d) + t V(z,d+1)]
 / [s + r₀ + p + t].
```

This expression is fractional linear in every varied scalar separately, with strictly positive denominator; iterating endpoint extrema is valid for signed values as well. Solve backwards in crown size and parent degree, round lower down and upper up, and set both boundary values to zero at first size exit. The resulting generator inequalities hold even for adaptive exterior evolution, by discounted stopped Dynkin. There is no independent draw of a parent's degree for a motif, generator term, or witness. At birth pool across the exact degree bins and the overflow bin, without a birth-law assumption.

At `N=8` the exact replay gives these rounded displays (the stored comparisons are rational):

| θ | Free exterior prefix | Tracked exterior prefix |
|---|---|---|
| 1/4 | `[1.49324e−6, 5.61792e−6]` | `[2.52453e−6, 3.81132e−6]` |
| 1 | `[−3.87012e−6, 2.60584e−6]` | `[−1.99128e−6, 2.18057e−7]` |
| 6/5 | `[−5.33789e−6, 1.82602e−6]` | `[−3.21368e−6, −6.72263e−7]` |
| 5/4 | `[−5.70643e−6, 1.63171e−6]` | `[−3.51973e−6, −8.94397e−7]` |

The negative upper bound at `6/5` is an improvement from tracking genuine outside transitions: the free-exterior comparison crosses zero. Only motif charges owned by the focal apex before the stop are paid. The bound is for the discounted prefix of a physical newborn, not the infinite lifetime or a bound on `πδ−/πδ+`.

### What a global transfer still needs

To turn the stop into an unlimited ratio, replace each zero terminal by lower/upper prices for the **same paired signed continuation**. A continuation interval valid at the first size exit bounds the signed terminal `E[e^{−sτ} V_exit]` in the same Bellman inequalities. A uniform upper terminal price `≤6.72e−7` would preserve the negative cap-8 comparison at `θ=6/5` because the discounted exit mass is at most one; the corresponding sufficient lower terminal price at `1/4` is `≥−2.52e−6`. These are sufficient seam targets, not measured errors. One can instead use a joint coboundary with zero stationary mean. The proof must cover physical tail rows, the higher motifs crossing the stop, and the conditional exterior at exit. Killed sign cannot be extrapolated by positivity of the two separate slacks: `Dθ` is signed, and the physical bulk examples show why universal ratio floors and ceilings cannot supply the required missing control. No new exponent endpoint is claimed.

## Reproduce

```sh
taskset -c 13 nice -n 10 /usr/bin/time -f 'cpu=%U+%S wall=%e' timeout 20s \
  ./round2/.venv/bin/python -Bu \
  round2/occupation-neon/ratio-transfer/replay.py b3f18_lo.bin b3f18_hi.bin
```
