# Independent review: absolute continuity of the Claim 6 limit

**Verdict: accept.** With the positive planted limits from `51a2a3e`, the proposed derivation proves absolute continuity. It also proves that the density is positive almost everywhere on `(0,∞)` using the already established full support (`open/limit-law/README.md`) and the argument in §4. No conditioning on future survival is used to assert independence.

## 1. A uniformly attainable fast edge

In a free tree with at least three vertices the cores induce a nonempty tree. Choose an endpoint `x` of that tree (or its only vertex). It has at least one leaf and at most one core neighbour. For one time unit insist on no ring of `x` and at least two rings of one specified leaf. That leaf cannot be promoted by another source, so its rings each append a leaf to `x`. After the silent unit, the first ring of `x` picks a leaf with conditional probability at least `3/4`: until that ring the stock cannot decrease, and its number of core neighbours is unchanged. At promotion time, `x` retains at least two leaves and the promoted vertex `y` has its new leaf. The probability of this setup is at least `q = (3/4)e^(-1)(1-2/e)`, independent of the seed. The two-vertex seed first reaches size three almost surely; the singleton is outside the model.

For each anchor track a family of currently genuine leaf neighbours, starting with one. Each tracked leaf sends to its anchor at rate one; give the newly born leaf its own fresh clocks and include it in that family if its parent particle is still tracked at the birth. Give each tracked leaf an independent marked death-proposal clock of rate `1/2`. On a proposal kill the lower-process particle regardless of the mark; in the physical tree promote that leaf only if the mark falls below `2/D_anchor`. Use residual physical directed-arrow channels for all other choices from the anchor, including the edge to the other anchor. A real leaf can only be promoted by its own anchor, so no unaccounted physical jump kills a tracked particle. The rate of each tracked physical promotion is exactly `1/D_anchor`, and all ordinary physical directed-arrow rates together are unchanged. The families are disjoint and their lower-process clocks are independent; acceptance is state-dependent in the physical tree and is **not** used to define the lower death. Fresh clocks may be assigned at discovery by the Markov construction; enlarging the tree via these channels is nonexplosive by the physical size rate.

Each lower process is linear birth–death with individual rates `1` and `1/2`. For completeness, `e^(-t/2)L_t` is an `L²`-bounded mean-one martingale: its compensated bracket has expected integral `(3/2)∫e^(-t/2)dt=3`. Its limit is positive with probability at least `1/4` by Cauchy–Schwarz. Independence gives a probability at least `1/16` for positive limits on both sides. On this event both anchor reciprocals have finite time integral: use the positive eventual coefficient of `e^(t/2)` for the tail and `D≥2` on bounded intervals.

Let `E` be the physical path event that some adjacent cores have both reciprocal integrals finite from their first time as adjacent cores. This is measurable in the increasing union of finite-time histories; throwing away a finite interval cannot change finiteness. At each integer time after size three, the conditional probability of `E` is at least `q/16` by the construction above, applied to the entire realized finite state. Lévy upward convergence gives `P(E)=1`. For each fixed label pair the cross count has compensator equal to the integral of the two directed rates when the edge exists. Finite compensator implies finite count on that event: stop when the compensator reaches integer `k`, use the bounded expectation of the stopped count, and take the union over `k`. This step does not require a finite unconditional expected compensator. Take a countable intersection over label pairs, then a countable union over candidate edges and deterministic integer cut times after the final cross. At the chosen cut both endpoints must be cores, so each half has at least two ordinary vertices.

## 2. Exact no-cross weight, with the slots retained

Fix the finite state and a core edge at a deterministic cut time. Run the two independent one-slot inert half paths with product law `Q`. Each internal arrow has its original physical rate, including arrows sent by the roots, since their slot counts toward degree. Use independent rate-one marked candidate clocks for the two suppressed crossing arrows and accept each at `1/D_root`, until the first accepted cross; one may continue the physical process separately after that cross. This constructs the correct physical law up to the first cross by equality of rates. Consequently the physical sublaw with no further cross has on the two-half path space the exact weight

```text
w = exp(-∫₀∞ [1/D_left(t) + 1/D_right(t)] dt) ≤ 1
```

relative to `Q`. The conditional Poisson calculation here is given the **unmasked inert paths**. The normalized no-cross path law generally is not a product. On `w>0`, the sublaw and `Q` restricted to that set are equivalent. The finite test gives a concrete unequal tilt: split a four-vertex path at its middle; an internal arrow growing the left root versus its leaf gives subsequent cross rate `5/6` versus `1` on an interval with no internal jump. Such prefixes can have equal internal rates and still get unequal weights. For the absolute-continuity argument only sublaw domination is required; equivalence on `w>0` is used below for positivity.

## 3. Smoothing the ratio

Under `Q`, condition on the second path and the first half's first genuine successor and entire continuation measured **from that jump**, without its absolute time. The race at the initial state has rate `λ=m−1/D_root`; its first genuine holding time is independent of the successor and, given the successor, independent of the restarted path. This follows by sampling exponential waiting time and winner separately, followed by a new Markov path. If `A,B` are the first restarted size and core limits, and `C,D` are the second half's size and core limits, all four are finite and strictly positive by `51a2a3e` (`open/forest-age/inert/README.md`). With `r=e^(-T)` of density `λr^(λ−1)` on `(0,1)`, the combined limit is

```text
f(r) = (B r^μ + D)/(A r + C)^μ,
f'(r) = μ (BC r^(μ−1) − AD)/(Ar+C)^(μ+1).
```

The numerator is strictly decreasing for `0<μ<1`, so there is at most one critical point. Off that point, cover by countably many intervals with locally Lipschitz inverse; the critical singleton has zero probability. This proves the pushforward under `Q` is absolutely continuous. Under no cross, core counts add exactly with the slot convention, ordinary sizes also add, and the combined ratio is the physical fixed-size limit since each genuine physical jump adds one vertex. A fixed cut time cancels in the ratio. For each null set, domination eliminates it on each fixed cut and state, and the countable cover in §1 eliminates it for the original law. In particular the proof handles a genuine interior turning point; a claim of global monotonicity would be false.

## 4. Extension: the law is equivalent to Lebesgue measure on `(0,∞)`

Call a planted path good if its root reciprocal degree has finite integral. This property, and the intrinsic ratio `b/a^μ` of its core and size limits, are invariant under deleting its first holding time and restarting at the genuine successor. Write `x₁,x₂` for the two intrinsic ratios. The combined ratio has form `p^μ x₁ + (1-p)^μ x₂` with `0<p<1`, so lies between `min(x₁,x₂)` and `2^(1−μ)max(x₁,x₂)`.

Fix `0<l<h<∞`. Full support gives positive physical probability of `X<l`, and of `X>2^(1−μ)h`. On each event apply the countable cover by fast edges and integer cut states, **retaining the fast-root conditions** in the subevent. Some fixed cut state and one side therefore give an inert seed with positive probability of a good path of ratio `<l`, by domination; possibly a different cut state and side give an inert seed with positive probability of a good path of ratio `>h`. Only the marginal existence is needed. Join the roots of those two seeds by making the two slots a real edge. The resulting finite tree is reachable up to labels from the edge seed with positive probability: prune to an edge and reverse, prescribing each attachment to a vertex with a neighbour that can send to it. Start the cut at this deterministic **size** hitting time, using the Markov property (or use a positive-probability finite-state cylinder at deterministic time).

For this joined state, sample each genuine successor and its restarted continuation under the product law. With positive probability both are good with their respective ratios `<l` and `>h`, independently of the two initial waits `T₁,T₂`. Conditional on such continuations, `r=exp(T₂−T₁)` has a strictly positive density on all of `(0,∞)` because the two independent exponential waits have positive rates. The common time factor cancels in the ratio, giving the same formula `f(r)` as in §3, now over `(0,∞)`. The two endpoint limits are the two intrinsic ratios, hence the range covers `(l,h)`. Off the critical point, change of variables gives positive probability of hitting **every** positive-Lebesgue-measure subset of `(l,h)`: the image of the regular pieces covers the interval except at most one value and the density of `r` is everywhere positive. Goodness gives `w>0` for every finite pair of waits; reweighting by `w` preserves this positive probability. Multiplying by the positive probability of reaching the chosen finite seed transfers it to the physical law. Exhausting `(0,∞)` by compact intervals proves the claim. Full support alone would not imply this; the ability to rejoin planted seeds and vary both waiting times is needed here.

## 5. Exact finite checks

```text
$ taskset -c 8 nice -n 10 timeout 8s python3 -B round2/checks/limit_density_sol.py
PASS: 872 labelled seeds, 3130 promoted-leaf cuts and exact arrow rates
COUNTEREXAMPLE: equal-prior path prefixes have cross hazards 5/6 and 1
PASS: exact two-branch, monotone, and endpoint-critical smoothing fixtures
```

The checks use independent rational arithmetic. The infinite-time statements are analytic, not simulation claims.
