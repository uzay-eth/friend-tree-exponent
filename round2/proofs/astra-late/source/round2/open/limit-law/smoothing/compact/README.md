> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../../README.md#claim-17) for exact scope and citations.

# Faster adapted cuts and a quantitative density estimate

For the edge-seeded physical core amplitude `X`, write `f_X` for its density. The method below gives:

- On each `[l,h]⊂(0,∞)`, `P(X∈A)≤C_l |A|^(4/243)` for Borel `A⊂[l,h]` of measure at most one. In particular `f_X∈L^p([l,h])` for `1≤p<243/239`.
- On the whole line, `P(X∈A)≤C |A|^(4/247)` for `|A|≤1`, and `f_X∈L^p` for `1≤p<247/243` (in particular `p=65/64`). This improves the previous [set exponent `1/80000`](../mixture/README.md) by more than three orders of magnitude.

This uses genuine early cuts at deterministic times in the original physical process. The key is to retry at the actual failure time, with a fast two-clock race for setup; the terminal no-cross weight is used only at the history-selected cut. The set estimate gives a uniform Hölder modulus for the distribution function, not a pointwise density modulus. Section 4 gives a multi-side kernel estimate and explicitly identifies the remaining mixing issue for continuity.

## 1. A physical edge with large stock, quickly and uniformly

Start at the three-vertex tree, which the edge-seeded chain reaches deterministically in shape. At each attempt start choose a core-tree endpoint `x` (including the lone-core case) and one of its leaf neighbours `ℓ`, using the present tree. Race their first rate-one clocks. If `x` wins, declare failure immediately, after processing the ring. Otherwise the ring at `ℓ` adds a leaf at `x`. Wait for the first subsequent ring of `x`, and declare failure if it selects a core neighbour. Before this ring the stock cannot fall and `x` still has at most one core neighbour; hence conditional success probability at the ring is at least `2/3`. At success the chosen leaf becomes the second anchor, with one new leaf, and `x` has at least one other leaf. No imposed fixed silent interval is needed. The first race duration has rate two and the identity of its winner is independent of that duration. The second wait is fresh rate one; its winning edge selection is made in the current state. All claims hold at every finite physical restart state with a core.

At this successful setup time start two disjoint lower leaf populations from one each, using the [independent physical embedding](../../density/README.md) (equivalently, see the [audited embedding](../../../../../../../../reviews/cold/astra/checks/limit-density-sol.md#1-a-uniformly-attainable-fast-edge)). Birth and death rates per tracked leaf are `1,1/2`; on a death proposal the track is erased even if its physical promotion was rejected. If a tracked family becomes extinct, restart immediately; otherwise never restart. The restart rules use only present information in the enlarged graphical filtration, and the construction at each restart has its original physical rates. At setup, the two lower family laws are independent, with their own drivers; the physical tree after additionally conditioning on one or both families surviving is not claimed to be independent of them.

Let `D` be the duration on failure, and put `r=1/24`. For one birth–death family with extinction time `T_0`, the [scalar law](../mixture/README.md#1-a-usable-edge-at-an-exponential-rate) gives `P(T_0<∞)=1/2` and `P(t<T_0<∞)≤e^(−t/2)/2`. For the minimum of the two extinction times on failure, the defective exponential moment is at most

```
F(r) = 3/4 + r/(1/2−r) = 37/44.
```

Indeed its defective tail is bounded by the sum of the two finite-extinction tails. Integrate the tail against `r e^(rt)` and add the failure probability `3/4`. Denoting `a_r=2/(2−r)` and `b_r=1/(1−r)`, we have, uniformly conditional on the restart history,

```
E[e^(rD); fail] ≤ a_r/2 + (a_r b_r/2)[1/3+(2/3)F(r)]
                = 11736/11891 =: ρ < 1.
```

The bracket is increasing as the setup success probability decreases, since `F(r)<1`. The bound uses the fresh family law conditional at the promotion time, not independence from the past degree path. The defective exponential moment of elapsed setup time on a never-failed attempt is at most `a_r b_r`. Thus multiplying conditional failure budgets by induction and summing over attempts gives `P(no never-failed attempt)=0` and `E e^(rS)≤a_r b_r/(1−ρ)`, where `S` is the setup completion time of the first never-failed attempt. `S` is **not** used as a stopping time or a cut selector.

At deterministic time `t`, take `K=ceil(e^(t/24))`. Except with probability `O(e^(−t/30))`, `S≤4t/5`. For any attempt reaching setup at real time `s≤4t/5`, the probability that either family has count in `[1,K)` at age `t−s` is at most `2K e^(−t/10)` conditional at setup, by the geometric law of the surviving count at a finite age. Sum over attempts, without a survival conditioning: the expected number of attempts *reached* is at most `1/(1−ρ)`. This is valid also for continuous setup times, by the strong Markov construction and the deterministic `t` in the age `t−s`. The lower processes of failed trials can be continued with their own laws if needed for the finite-age bound; only the never-failed attempt is used to certify real stock at time `t`. Thus except with probability `O(e^(−t/30))` the current tree has a core edge with at least `K` real leaves on each side incident to its respective endpoint. Select one such edge using only the current tree. The [first-send barrier](../../../age-order/early/README.md#a-no-send-bound-without-an-age-cutoff) gives a probability at most `4/(K+1)` of any future cross, conditional on the actual state and chosen edge. Therefore the physical no-cross retained sublaw at time `t` loses mass at most `C e^(−t/30)`. This is a quantitative bound on bad terminal cuts, not a retrospective selection of a permanent edge.

## 2. Beta kernel averaged over genuine cut sizes

Conditional on the finite time-`t` history with an eligible edge, use the [exact cut identity](../README.md#2-exact-dirichlet-image-including-the-core-limits). Under its unmasked reference, given the two embedded paths, the total amplitude is `g(P)=aP^μ+b(1−P)^μ`, with `P~Beta(n_1,n_2)`, the two ordinary cut sizes as beta parameters, and terminal shape weight `0≤L≤1`. We use only the upper bound by the reference kernel, never conditional independence under a shape-and-clock-dependent selector. At least one coefficient is positive on the weighted sublaw.

For integer `n_1,n_2≥2`, beta density is `(n_1+n_2−1)` times the binomial point likelihood with `n_1−1` successes and `n_2−1` failures, maximized at the empirical frequency. Stirling there gives its supremum at most `C (n_1+n_2)/sqrt(min(n_1,n_2)−1)` with an absolute constant. Both cut components contain the endpoint and its `K` real leaves, so the bound is at most `C |T_t|/sqrt(K)`. This is not a replacement of the physical cut-size law by a deterministic size cap.

If `A⊂[l,∞)` and the image intersects `A`, `a+b≥l`. The concave map has `−g''≥μ(1−μ)(a+b)`; separate the region of derivative magnitude at most `h` and the at most two monotone inverse branches outside it. Optimizing in `h` gives `Leb(g^−1(A))≤4 sqrt(|A|/[μ(1−μ)l])` for arbitrary Borel `A`. Averaging the product of the beta density bound and this sublevel bound over the actual physical history, using `E|T_t|=3e^t`, yields

```
P(X∈A) ≤ C e^(−t/30) + C_μ e^((47/48)t) sqrt(|A|/l).
```

There is no expectation conditional on future absence of crossings, no size cutoff and no sum over retrospectively successful edges.

## 3. Exponents and full density

Set `c=1/30`, `d=47/48`, `δ=|A|`. For the compact bound choose `t=log(1/δ)/[2(c+d)]`; then both powers are `ξ=c/[2(c+d)]=4/243`. For the global bound, split off `[0,l]` using the [left-tail estimate](../../left/README.md#matching-left-tail-order-and-the-negative-moment-threshold) `P(X≤l)≤C l^(1/μ)`. Choose `l=δ^(μξ)` and `t=(ξ/c) log(1/δ)` with `ξ=c/[2(c+d)+μc]≥4/247`. Null sets follow by outer approximation (absolute continuity is also known). For either set bound, on `{f_X>H}` in the relevant domain, `H Leb{f_X>H}≤P(X∈{f_X>H})≤C Leb{f_X>H}^ξ` at large `H`. Layer cake yields the stated open `p` ranges. In particular the proof does not identify a finite density supremum or rule out spikes.

## 4. Multi-side test and the remaining continuity input

A concrete six-side estimate uses three pairs. Conditional on the embedded paths and pair sums `s_j`, independent within-pair beta variables split each pair, by the gamma representation of the Dirichlet law. Suppose each paired coefficient sum `A_j>0`; let `M_j` bound its beta density. The sublevel proof gives `P(pair image∈E)≤C_μ M_j |E|^(1/2)/(A_j s_j^μ)^(1/2)`. For a probability density with set bound `B sqrt(|E|)`, the superlevel volume is at most `min(1/H,B²/H²)`; integrating at power `3/2` gives `||f||_(3/2)≤C B^(2/3)`. Thus the threefold convolution has a bounded continuous representative with bound

```
C_μ ∏_(j=1)^3 M_j^(2/3) (A_j s_j^μ)^(-1/3).
```

This follows from Young's inequality and continuity of translations in `L^(3/2)`. For fixed shapes, integration in the pair sums is legitimate: writing `N_j` for the pair shape totals and `N=ΣN_j`, the factor is `Γ(N)/Γ(N−μ) · ∏ Γ(N_j−μ/3)/Γ(N_j)`, finite since `N_j≥4`. The same argument gives a bounded/continuous **fixed-cut sublaw** if its history/shape average of the displayed bound is finite, using dominated convergence; an arbitrary future selection mask allows only domination, not continuity of its masked density. The estimate is conditional and does not assume uniform negative moments of real cut coefficients or any bound on the density of the discarded law. Balanced parameters and equal coefficients already yield order `1/N` concentration at the stationary maximum of the Dirichlet image as sizes `N` grow ([symmetric obstruction](../README.md#4-arbitrarily-large-permanent-cutsets-and-the-summation-obstruction)). So uniform fixed-cut bounds independent of sizes cannot solve this step.

The [sharp Q3 birth-mass density](../../../age-order/eventual-domination/README.md) concerns a stationary one-birth intensity. Its large-mass tail does not by itself bound the inverse pair coefficients or their mixing law under genuine early cuts. Even in the [quenched marked-count law](../../../age-order/spectrum/quenched/README.md), the normalized mark intensity removes the random factor `X`; the annealed tree law instead size-biases it. Substituting either intensity for a conditional law of cut coefficients would ignore the shared physical future. The power estimate above bypasses this issue; boundedness or continuity of the full density on compact sets remains open here.

## Bounded replay

```text
taskset -c 8 nice -n 10 timeout 8s python3 -IB round2/open/limit-law/smoothing/compact/check.py
PASS defective retry moment <= 11736/11891; missing-mass rate 1/30
PASS 612 beta profile comparisons, including unbalanced cuts
PASS compact/global exponents, three-pair convolution, symmetric concentration scales
```

These are finite arithmetic controls; the density estimates follow from the analytic retry and cut arguments.
