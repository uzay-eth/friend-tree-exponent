# D20 signed corrector — cold audit

**Holds: `[0.565681980098700, 0.565766097351300]`.** The fixed-trial pointwise inequalities have strictly positive rational margins after keeping the full joint size-16 reserve. The stationary incidence LP also holds and gives the weaker outward interval `[0.565681955022224, 0.565766123513244]`. These are alternative uses of the same signed charges. The audit pins the local `corrector/` implementation by hash; `8790` did not resolve as a local Git object.

## Authenticated inputs and denominator

The two 11,190,881,056-byte originals were rehashed in both the producer and audit replays. The source is `ft2-neon@51f4cfadac675fdf1cef84d44935eebb5456aabb`: its proof, binary format, GMP verifier source and two logs match the pinned Git blobs. Both logs tie the same payload hash to `PASS (total_fails = 0)` for all eight families, the comparison and 34 concrete/boundary trees. The full GMP result is inherited; the audit does not rerun the exhaustive verifier. The authentic discounts are `λ_lo=.565680`, `λ_hi=.565768`; sign `σ=+1` lower and `−1` upper. The exact scales are `S=10^11`, `Sr=10^9`.

| Side | Payload SHA-256 |
|---|---|
| lower | `eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469` |
| upper | `55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2` |

The full comparison majorizes the original **potential**, not the signed charge. It includes all eight families: bottom charges for `g/hV/P4/Yb` and the proper apex/middle multiplicities for `hC/Yt/S3`, in addition to `f`. The audit replays the independent D20 occupation routines (`neon20-occupation-sol.py`) for all 404 cap groups, all 2,079 finite nonfounder census envelopes and singleton affine rows, and the five all-stock reward floors. The recharge `t(k−h)` in ticks sums to zero on the whole oriented core tree; its nonpositive `t` controls true arity overflow beyond the last stock class. The full cap `M` bounds the recharged envelope everywhere; at exact censuses the sharper bound is `Q(x)`. A nonnegative birth-Palm deficit payment at discount HI with prefactor LO yields the stationary upper denominator `C=M−deficit`:

| Side | `t` (ticks) | `S M` | Exact `C` |
|---|---:|---:|---|
| lower | `−200000000` | `181661912793` | `16959459852924085304227707/12500000000000000000000000` |
| upper | `−1000000000` | `183051705651` | `17064269776474933055949939/12500000000000000000000000` |

The pointwise inequality uses `M` together with the census deficit; subtracting the separately priced deficit from `M` there as well would double-count it. `C` is the full denominator for the stationary direct quotient and residual-to-increment step. The positivity check uses the original GMP **lower** comparison `c0` and the uniform norm bound `max|W−W(J)|+max E<c0` (at most one retained edge per nonfounder child).

## Corrector and actual reserve support

For an existing nonfounder mother–child edge, each vertex's gain uses the other endpoint's true reciprocal, its leaf stock and its exterior receive variable. The latter is the very same variable in the affine signed `g` reward. With the cutoff `D≤16, k≤5` at both ends, gains and promotions move forward, and a crossing still contributes its exit rate but has zero continuation. The old edge process includes promotions of its endpoints, without any newborn edge entrance. Testing the corners of the pressure rectangle plus intersections of its sides with the zero-charge line suffices: the cross-multiplied residual is affine on each piece of the hinge subdivision. The positive-part replacement is valid using the inherited nonnegativity of the **physical** charge. The independent checker clips each rectangle side using its two signed endpoint values, and checks interior rates as well.

A physical arrow into a leaf turns that leaf into a core of degree 2 and arity 0 (`J`). Its mother **keeps her degree** and raises arity by one. The new correction `E((D,k+1),J)` is placed once in her singleton equation at rate `L/D`. All pre-existing incident edge values also change under that same promotion. The recentered newborn singleton correction `W(J)−W(J)` is zero. This does not remove or re-credit the newborn **production potential** already allocated to `f`. We retain the founder's actual arrows, exclude her corrector and founder-apex edges, and allow her promotion to create a new nonfounder with zero initial singleton correction. A founder-origin boundary term in spatial normalization is bounded per founder and vanishes per core.

The reserve is for the **actual complete planted fringe** of the apex, with ordinary size `≤16` and at most four cores. Its independent finite-fringe replay retains both the old `g` row and all six higher families jointly, with one physical exterior receive in the reward and transition. Each off-edge core child requires itself and at least one further vertex. Thus for a parent–child edge its apex fringe has size at least `2+L_parent+L_child+2(k_parent−1+k_child)` and at least `1+k_parent+k_child` cores. The corrector pays `g` only when one of these **lower bounds** excludes the reserve. The converse is false: the audit finds actual out-of-reserve fringes that are still unpaid by this sufficient mask. This is conservative by GMP nonnegativity. A masked-out edge may still have positive continuation to paid edges. That continuation is not an immediate second payment. Exact joint reserves are `2859476002446237549/50000000000000000000000000` lower and `30885641464886320713/1250000000000000000000000000` upper.

## Stronger fixed-trial route

Put `s=λ+σδ` and `B=Σ_nonfounder (W−W(J))+Σ_nonfounder edges E`. From the Bellman tests, the signed drift of `Φ+σB` at `s` is bounded below by `sW(J)−δM` per nonfounder plus the original unspent reserve charge and the founder boundary. Indeed the old signed drift allows `f` and masked `g` to pay for `(s−G)B`; the `δ(M−Q)` term offsets `−δΦ` with the zero-sum-recharged comparison. Above the finite census use the all-stock reward floors and conservative recharged caps. There `W` is degree-independent, the unsaturated promotion rate is in `[max(0,1−(k+1)/65),1]`, and on saturated arity a promotion does not change the old vertex's stored value. There are no edge births priced by `E` above this singleton census.

To take the exponent comparison, normalize the signed drift by discounted expected core occupation as the discount decreases to `μ`, as in the prior D20 occupation audit. This avoids assuming convergence of an unbounded pointwise generator: integration by parts for `Φ+σB`, whose positive lower and finite upper bounds are proportional to core count, yields `σ(μ−s)` times its normalized mean up to a vanishing initial term. The separate bounded reserve gets its uniform birth-Palm lower payment; the founder's contribution vanishes. A positive combined margin excludes the wrong sign of `μ−s` on either side. In particular this step does not require the reserve to be positive on every physical census.

| Side | Exact `δ` | Exact positive margin `sW(J)−δM+reserve` |
|---|---|---|
| lower | `19800987/10^13` | `1329501250839720657049503/125000000000000000000000000000000` |
| upper | `19026487/10^13` | `22036192591132498427015907/2000000000000000000000000000000000` |

The trial rationals are already on the 15-decimal display grid. The audit computes both margins by a fixed-trial backward pass, without the producer's trial search.

## Stationary LP route

The LP keeps exact `(D,k)` for `D≤40` and a grouped infinite-degree tail. Its physical lift uses the same `μ` mixture in the occurrence columns; total probability, total child stubs and total edges each normalize to one. The newborn forcing in the age rows is `μh(J)`; promotions at saturated tail arity leave the old class unchanged. Tail degree gains leave the aggregate class unchanged, and the finite gain seam clips the target arity. For unsaturated tail arity, `b=1−(k+1)p` with `p=1/D`; at saturation `0≤b≤1−9p`. These are inside the column hull, with nonnegative excess child count and excess child receive accounting for the omitted higher arity. The child and parent incidence moments share one edge: each edge supplies the child's reciprocal to the parent receive marginal and the parent's reciprocal to the child receive marginal. Its two reciprocals lie in the respective coarse degree intervals; a convex lift into their tensor endpoints preserves both moments. Mass transport supplies the source and normalization identities `πb=μ`, `πρ=1−μ`, `πk=1`. No equality of conditional mean products is imposed.

The joint objective is `f+bE(post,J)+δ(M−Q)`, with true shared receive at exact censuses. The existing-edge Bellman inequality at HI also holds at `μ≤HI` since `E≥0`. Applying it to edge occupation gives the lower payment `π[bE(post,J)]` from the true stationary edge entrance. At a grouped tail the LP takes the objective minimum over all covered finite censuses and all-stock classes; dropping the nonnegative edge-birth term there is safe (in this run its cutoff is already below the tail). The audit captures the HiGHS proposal, recalculates the rational objective on every column, pairs the rounded dual with the exact sparse rows, and subtracts the maximum positive defect of each of the four groups. Each group mass is at most one: occurrences and edges by normalization, excess stubs by `πk`, excess receive by `πρ≤1`. This certifies an objective on the physical lift regardless of floating point errors. The signed residual at trial `δ0` bounds `(σ(μ−λ)−δ0)πΦ` from below; a positive residual buys at least `residual/C`. This verifies the reported LP interval with the full denominator. Both LP replays retain the same joint reserve as the pointwise route.

## Reproduction

From the project root run for each corresponding `side` and `payload` (`b3f20_lo.bin`, `b3f20_hi.bin`):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 timeout 250 round2/.work/venv/bin/python -Bu round2/occupation-neon20/corrector/replay.py --side "$side" --mode preserve --edge-degree 16 --arity 5 --stationary 40 --fixtures --exact-den
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 timeout 240 round2/.work/venv/bin/python -Bu round2/checks/neon20-corrector-sol.py "$side" "$payload"
```

Both producer replays passed, printing the pointwise and LP intervals above, 5,200 edge states, 81 exact production drift fixtures, and the 1,109-row / 11,394-column stationary model per sign. The independent replays also passed both signs, with the exact denominators, reserves and positive margins displayed above. Each sign checked 5,200 edge states (4,739 eligible, 461 unpaid with positive continuation), 3,914 retained reserve states (1,259 size exits and 9,060 core exits), and 99 labelled ordinary-arrow fixtures with full signed production drift equal to the sum of the eight charged families. Per sign, 275 nonzero newborn terms were seen; using a pre-promotion mother or deleting the newborn correction failed 95 fixtures each. The fixtures contained 40 reserve-eligible edge incidences (three at exact size 16), 22 census-certified outside incidences (one at size 16 with five cores), and eight actually outside but not paid by the census mask. The LP audit printed outward integer ticks `565681955022224` and `565766123513244` and matched the exact grouped repairs.
