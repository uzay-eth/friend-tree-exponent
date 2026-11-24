> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# Adjacent labels at positive birth-scale thresholds

For every finite `K`, there is a label `N(K)` such that for `n≥N(K)` and **every** `0≤y≤K`,

```math
 \Pr(Z_{n+1}>y/n)>\Pr(Z_{n+2}>y/n).
```

In fact both the single-birth and the adjacent-difference scaled tails converge uniformly on such intervals. For fixed unscaled `0<x<1` a separate argument below gives the exact log tail `v⁻¹ log Pr(Z_v>x)→log(1−x)`, uniformly on compact subintervals of `(0,1)`. That log limit does not order adjacent tails at fixed `x`.

## 1. Positive continuous-time mass has no atoms, even in an infinite exterior

Fix any locally finite initial marked tree and focal vertex (possibly a leaf). Write `D(t)`, `L(t)`, `c(t)` for its degree, leaf stock, and number of core neighbours. The degree limit `W=lim e⁻ᵗ D(t)` exists by the nonnegative supermartingale property: the gain rate is at most `D`. For a core, `c(t)` can grow only by its own rate-one initiations; more generally the initially degree-one endpoint can have one initial neighbour become a core without such a promotion, or one new leaf added at its own first promotion. The bound `c(t)≤c(0)+1+Pois_clock(t)` suffices. In particular `c(t)=O(1+t)` a.s., and `e⁻ᵗ L(t)→W` on `{W>0}`.

Fix `w>0` and integer `t`; write `z=eᵗw` and `L=L(t)`, `D=D(t)`. Given the history, the remaining-time mass `eᵗW` is stochastically bounded above by `Gamma(D,1)` and below by `Gamma(L,1)−E`, with `E∼Exp(1)` and no independence needed for the latter pair. These are the [degree and tagged-leaf comparisons](../../../../../../astra/source/round2/fringe/hub-spectrum.md#2-a-vertexs-limiting-mass-as-seen-at-time-t); their constructions use only local directed channels and apply in the locally finite exterior. The upper comparison uses a Yule chain; the lower keeps phantom Yule descendants of each lost tracked leaf and dominates the lost stock by starting a family at every focal ring. Thus, for `L≥1` and `a>0`,

```math
 \begin{aligned}
 \Pr(W=w\mid\mathcal F_t)
 &\le \Pr(G_L\le z+a)+e^{-a}-\Pr(G_D<z)\\ &\le \|f_L\|_\infty(a+c(t))+e^{-a}.
 \end{aligned}
```

For the second line use the convolution representation `G_D=G_L+G_c'` with independent gamma increments **only for these marginal gamma laws**, and bound the cdf shift by `E G_c'=c(t)` times the supremum gamma density. This density is `O(1/√L)` by its maximum at `L−1` and Stirling (include `L=1` separately). On `{W=w}`, eventually `L≥(w/2)eᵗ` and `c(t)=O(t)`. Take `a=t`; the right side tends to zero on that event. But `Pr(W=w|F_t)→1_{W=w}` a.s. by upward convergence. Therefore `Pr(W=w)=0`. No cutoff at infinity, independence of hub futures, or finite global population is needed. The statement allows an atom at zero and makes no density claim.

## 2. Uniform convergence and sign on the full compact scale

Use the [late-birth decomposition](../late/README.md) for the single-birth score and the [audited oriented charges](../../../../../../../reviews/cold/astra/checks/adjacent-age-sol.md) for the signed two-birth score. With endpoint test `W>y` in the locally finite tree after each insertion, denote the scores by `F_y` and `G_y`. The pointwise bounds are `0≤F_y≤7`, `|G_y|≤3` (also the sum of absolute oriented and chain charges before inserting the endpoint differences is at most 3). Put

```math
 A(y)=\mathbb E[X]\pi(F_y),\qquad B(y)=\mathbb E[X]\pi(G_y).
```

For `y>0`, the endpoint event is continuous under local convergence even when the threshold varies toward `y`. In detail, for `e⁻ᵗD(t)` the absolute remaining integrated drift in expectation is at most `(c(0)+t+2)e⁻ᵗ`, by the bound on `c` (also at a leaf). Its martingale has expected remaining quadratic variation at most `D(0)e⁻ᵗ` since the degree gain rate is bounded by `D(t)` and `ED(t)≤eᵗD(0)`. This gives local uniform `L¹` approximation of `W` by finite-time endpoint reads. The [old-cell cutoff](../../../../../../astra/source/round2/operator.md#1-cells-and-evaluation-of-the-operator) gives continuity of those reads in law, also after each finite insertion. Section 1 eliminates the positive-threshold atom. In a finite tree of post-insertion size `m`, the actual endpoint is `mZ=(m/S)W`, with `S∼Gamma(m,1)`. The factor `m/S` converges in probability to one uniformly over initial shapes; the threshold sandwich does not require its independence from `W`.

At `y=0`, the endpoint test is hubhood. Its continuity follows from the hitting and stock-charge squeeze in the [adjacent audit](../../../../../../../reviews/cold/astra/checks/adjacent-age-sol.md#continuity-for-both-newborn-endpoint-types); that argument includes the initially degree-two middle. For any nonnegative sequence tending to zero, the endpoint probabilities lie between the zero-threshold hub test and every fixed small positive-threshold test eventually; monotone convergence in the limiting input as the positive threshold decreases to zero completes the squeeze. This also applies when the finite-tree endpoint is `mZ`. Thus the integrated scores are continuous on `[0,∞)` including right continuity at zero. For a convergent rooted sequence, the finitely many placement types and their coefficients eventually agree locally, even if their multiplicities are large at that particular root.

Rooted empirical weak convergence, followed by core-count uniform integrability and the score bounds, now gives for **every** sequence `y_n→y≥0`

```math
 \begin{aligned}
 n^{1-\mu}\Pr(Z_{n+1}>y_n/n)&\to A(y),\\ n^{2-\mu}\{\Pr(Z_{n+1}>y_n/n)-\Pr(Z_{n+2}>y_n/n)\}&\to B(y).
 \end{aligned}
```

For the varying test, apply the sequential continuity just proved to a coupling of the rooted empirical measures with `π` realizing their almost-sure local weak convergence (or use the extended portmanteau theorem on this Polish space). The denominator in the signed adjacent identity is `n(n+1)` and the finite-tree multiplier uses `m=n+2`; neither change affects the displayed limits. By subsequence compactness and continuity of `A,B`, both limits are uniform for `y∈[0,K]`.

Telescope adjacent differences from `n` to `⌊rn⌋` at the **same** threshold `y/n`. Uniform convergence on `[0,ry]` and Riemann integration give, for `r>1`,

```math
 A(y)-r^{\mu-1}A(ry)=\int_1^r s^{\mu-2}B(sy)\,ds.
```

Since `A` is nonincreasing, dividing by `r−1` and sending `r↓1` gives `B(y)≥(1−μ)A(y)` at every `y≥0`. At zero there is equality. For any finite `y`, `A(y)>0`: under the stationary local law the core-initiation weight for choosing a leaf is `π(L/D)=μ>0`. After such a placement, finitely many choices promote the marked tip and grow its leaf stock to any prescribed integer in one bounded time interval with positive probability. The tagged gamma lower comparison from there gives positive probability of mass above `y` by taking the stock sufficiently large. Consequently `inf_[0,K] B≥(1−μ)A(K)>0`. Uniform convergence proves the claimed strict comparison. It also implies `A` is continuously differentiable on `(0,∞)` by the same integral identity and `B(y)=(1−μ)A(y)−yA'(y)`.

## 3. Fixed unscaled thresholds: the exact exponential rate

Let `v` denote the marked birth label, and start the clock just after its birth. The upper bound `(1−x)^(v−1)` is the degree-one Pólya comparison in the [parent note](../README.md). The following matching lower bound on the exponential scale uses a physical two-vertex cut and the **actual** size limit. It avoids an independence claim about physical mass and physical size.

### A split with uniform exterior control

Write `u` for the marked leaf's birth parent. Require the marked leaf to initiate the next jump (probability `1/v`), placing another leaf on `u`, and then `u` to select the marked leaf at the following jump (probability at least `1/(v+1)^2`). This uniform polynomial-cost setup produces a two-vertex side (mark and its new leaf) and an exterior of `v` vertices, where `u` retains a leaf. Both cut endpoints are cores in the full graph. Restart the clock here. In the **inert** pair suppress the two cross arrows, keep the exterior slot in each boundary degree, and retain independent side drivers. The conditional probability of no cross forever is `exp[−∫(1/D_mark+1/D_u)dt]` on the inert paths, because each cut proposal can be given an independent Poisson driver with predictable thinning. Either side has size gain rate `N−1/D_boundary`, so its normalized size limit exists by the nonnegative supermartingale property. The lower bound uses events making both reciprocal integrals bounded; a blanket no-cross assumption would not suffice.

Embed the exterior inert process into a Yule proposal clock of total rate `N_+` with unit jumps, initial value `v`. At each proposal choose a uniform slot out of `N_+`; the `N` occupied slots represent exterior initiators. If occupied, use the physical uniform neighbour choice, including the exterior slot. On a rejected cross or dummy proposal increment only `N_+`; on acceptance increment both. Proposal holding times are independent of the marks and embedded exterior evolution. Under the normalized exponential tilt with weight `e^{-\theta S_+}`, `\theta>0`, the Yule rate is `N_+ q(t)` where

```math
 q(t)=(1+\theta e^{-t})^{-1}.
```

Indeed the conditional Laplace transform at time `t` is `(1+θe⁻ᵗ)^−N_+(t)`; its jump ratio is `q(t)`. Thus the exterior physical rates are multiplied by the same deterministic `q(t)`; in the clock `∫_0^t q(s)ds` they are ordinary, and `q(t)≥(1+θ)⁻¹`. On survival of the embedded tracked family `J`, `∫_0^∞ dt/(1+J(∫_0^t q(s)ds))` is finite a.s. Its tilted distribution gives, for any fixed `θ`, a uniform positive probability (over all admissible exteriors) of this integral being at most some finite `M_θ`: the dominated family has the same rates in the time-changed clock. Under the tilt `S_+` is `Gamma(v,1+θ)`, so intersecting the bounded-reciprocal event with `|S_+/v−1/(1+θ)|≤ε` still has probability bounded below for all large `v`. This intersection uses a union bound, not independence.

### Marked side: a tail with rate as close to one as desired

For any `z>0`, choose stock `k_v=ceil((1+η)ze^h v)` at time `h`, where `η,h>0` are fixed. On the inert marked side silence the focal rate-one clock until `h` (probability `e⁻ʰ`). The remaining clocks grow a pure star; its Yule leaf count at `h` has geometric tail `(1−e⁻ʰ)^(k_v−1)`. The outward reciprocal integral up to `h` is at most `h`. After `h` use the ordinary inert dynamics. At stock at least `k_v`, track just `k_v` leaves. With conditional probability tending to one uniformly in such states, its limiting mass from the cut is at least `zv`: use `e⁻ʰ(G_{k_v}−E)` and the gamma lower concentration bound. Also the outgoing reciprocal integral from `h` is at most `1` with probability tending to one: the [stock envelope](../README.md#tagged-leaf-bounds-used-for-the-limiting-comparisons) bounds failure by `O(1/k_v)` and on success gives `∫_h^∞dt/D_mark(t)≤2/k_v`. Consequently the mark-side event of mass at least `zv` with total outward reciprocal at most `h+1` has log probability per `v` at least

```math
 -(1+η)z\, e^h[-\log(1-e^{-h})].
```

The prefactor multiplying `z` decreases to one by taking `h→∞` and then `η↓0` (after the `v` limit). Crucially, on the marked side its ordinary size mass is at most `(1+o(1))` times its endpoint mass in this lower bound, as follows. Up to `h` its size is `k+1` on a star of `k` leaves. After `h` the size is dominated by a Yule process started at `k+1`, whereas the endpoint is bounded below by the gamma-minus-loss comparison from `k`. Thus on the star event, with conditional probability `1−o(1)` uniformly for `k≥k_v`, the cut-time mass `S_{
m mark}` satisfies `S_{
m mark}≤(1+η)W_{
m mark}` by simultaneous gamma concentration and the loss tail. No independence of these two estimates is needed.

### Optimize the actual ratio

Under no cross the two ordinary size masses add, and the physical fraction is `W_mark/(S_mark+S_ext)`. This identity uses the common restart clock, so the durations of the two setup jumps cancel. Set `s=(1+θ)⁻¹∈(0,1)`. Undoing the exterior tilt on the band `|S_+/v−s|≤ε` costs log probability per `v` at most `I(s)+θε`, where `I(s)=s−1−log s`: indeed the inverse likelihood is `(1+θ)^−v exp(θ S_+)`. There `S_ext≤(s+ε)v`. Choose `z` so that `z/[ (1+η)z+s+ε]>x`. The independent inert side events and bounded no-cross factor give, after the polynomial setup cost and sending the auxiliary errors to zero,

```math
 \liminf_{v\to\infty} v^{-1}\log\Pr(Z_v>x)
 \ge -I(s)-\frac{x}{1-x}s.
```

The minimum of `I(s)+xs/(1−x)` is at `s=1−x`, where the value is `−log(1−x)`. Together with the upper cap this proves the asserted log limit, even uniformly over the conditional birth tree because all lower bounds were uniform in that input. Uniformity over a compact interval of `x` follows by finite grid monotonicity and continuity of `log(1−x)`. The tilt is on an auxiliary upper size, while the success event explicitly bounds the two **actual** size masses. Neither the rate nor the compact birth-scale comparison determines the sign when `nx→∞`.
