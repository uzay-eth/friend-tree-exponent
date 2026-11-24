> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports. See the [release index](../../../../../../../README.md#claim-24) for exact scope and citations.

# Large labels, uniformly up to the upper edge

For each `0<p<1` the [interior equivalent](../README.md) extends uniformly to **all `p≤x<1`**, with the same constant:

```math
 \Pr(Z_v>x)=(1+o(1))\frac{\kappa}{\sqrt{x}}
 v^{\mu-3/2}(1-x)^{v-1/2},\qquad
 \frac{\Pr(Z_{v+1}>x)}{(1-x)\Pr(Z_v>x)}=1+o(1).\qquad\text{(A)}
```

The deterministic errors are uniform for `x∈[p,1)`. Thus one label cutoff works for every such `x`, even arbitrarily close to one. There is no claim at small `x`. For every leaf-parent birth state, without any local compactness requirement, its physical tail is between `K_-` and `1` times `Beta(1/2,v−1/2)` survival on this range, for all sufficiently large `v` depending only on `p`. One choice is `K_-=exp(−π_math²/3)/2`. The comparison alone, with the adjacent beta ratio and the core-placement estimate, gives the sign sufficiently near one.

Credit: Fable 5(c) gives finite full domination by rooted transport (see [audit](../../../../../attacks/fable5/flow/README.md)); 5(e) gives the Pólya comparison. We use the exact likelihood sandwich and limiting old-side process of the [reviewed interior theorem](../README.md), but change the clock and prove new estimates uniform as the auxiliary endpoint tends to one. The [fixed-seed near-one Palm law](../../../near-one-constant/README.md) has a different order of limits; the connection is in §5.

## 1. Group by outside jumps

Use a physical leaf-parent birth at size `v`, with Pólya red mass `a=1/2` and blue mass `b=v−1/2`. Conditional on the reference endpoint `Θ=θ∈[p,1)`, set `u=1−θ`. Before outside jump `f` (indexed from zero), the block contains `G_f` tag choices then one outside choice. The `G_f` are independent geometrics on `{0,1,…}` with survival `θ^k`. At the start of the block let `R_f=a+Σ_(h<f)G_h`, `B_f=b+f`, and `Δ_f=R_f−r_f`, the excess of reference red mass over the actual tag gain rate. It is constant through the tag choices, nonnegative and nondecreasing over outside jumps; `Δ_0=0`. Put `τ_n` immediately after `n` outside jumps and `L_n` for the likelihood at `τ_n`. Then

```math
 \log L_n=\sum_{f<n}\left[\log(1+\Delta_f/B_f)
           +\sum_{j<G_f}\log(1-\Delta_f/(R_f+j))\right].\qquad\text{(1)}
```

For the outside arrow, use its actual rate divided by `B_f+Δ_f`, evaluated after this block's tag choices. This preserves the physical conditional outside law. The start of block `f+1` has red mass `R_{f+1}=R_f+G_f`. Under the reference, `τ_n` is almost surely finite. Under the physical law, from any state blue mass is at least `b`, and away from outside jumps the hazard is at least `b/j` at size `j`; the product of survival factors vanishes. Inductively every `τ_n` is finite under the physical law too. The fixed-prefix upper sandwich therefore extends by summing over the disjoint countably many histories at `τ_n`. For the lower sandwich use the dominated measure of density `L_∞` from the reviewed proof; its value is the almost-sure limit along this subsequence. In particular, writing `Q_v(x)=Pr(Θ>x)`, we have

```math
 \mathbb E_{\rm ref}[L_\infty\mid\Theta>x]
 \le \Pr_T(Z>x)/Q_v(x)
 \le \mathbb E_{\rm ref}[L_n\mid\Theta>x].\qquad\text{(2)}
```

Both inequalities allow an inserted prefix indicator measurable at `τ_n` (and at all larger cutoffs). For the upper inequality the future comparison starts from the same `R,B` at that realized history, using the locked deficit at birth.

## 2. Uniform likelihood approximation

### Bounds in the outside clock

The number of tag-origin outside arrows before block `f` is dominated by `J_f`, a sum of independent Bernoullis of parameters `1/(b+h)` for `0≤h<f`. This holds even conditionally on the entire sequence of waits, by proposing the tag-origin decision with a fresh uniform at each outside jump. The old parent's contribution to `Δ_f` is at most `1/2`, and every new core neighbour costs at most one and comes from a tag-origin arrow. Consequently `Δ_f≤1/2+J_f`. The first increase of the old parent degree has rate at most two until it occurs, so a union bound over outside steps gives the sharper expectation estimate

```math
 \mathbb E_\theta\Delta_f\le 2\sum_{h<f}(b+h)^{-1}.\qquad\text{(3)}
```

In the first `m` outside choices, the probability of any change to `Δ` is at most `3m/b` (parent gain or tag-origin arrow). Both statements are uniform over initial shapes and waits. They retain the incoming slot to the tag at the old parent; a parent's choice of that slot causes a tag attachment, not a change to the parent degree.

Let `m=ceil(εv)` with fixed `ε>0`. With probability at least `1−C/H²` the following hold for all `f≥m`, uniformly in `θ∈[p,1)` (the probability bound does not need to depend on `ε≤1`):

```
|u R_f−θ f| ≤ H sqrt(v) (1+f/v)^(3/4),
Δ_f ≤ 1+H(1+f/v)^(1/8).
```

The first uses `E[uG_f]=θ`, `Var(uG_f)=θ≤1`, the maximal variance inequality on dyadic blocks of scale `v`, and the harmless `ua` initial term. The second follows from `E[J_f^2]≤C(1+log(1+f/v))²` and a dyadic union bound. For fixed `ε,H,p` and large `v` these imply `R_f≥pf/(2u)≥pf/2`, `Δ_f/R_f≤1/2`, and

```math
 e_f:=\log\frac{uR_f}{\theta f},\qquad |e_f|\le C_{p,\varepsilon} H v^{-1/2}(1+f/v)^{-1/4}.\qquad\text{(4)}
```

### Cancellation across a block

On this event replace (1) over blocks `m≤f<n` by `Σ Δ_f[B_f⁻¹−log(R_{f+1}/R_f)]`. The outside linearization costs at most `Σ (Δ_f/B_f)²`. The red linearization and the error of replacing the discrete harmonic sum by the log together cost at most

```math
 C\sum_{f=m}^{n-1}\Delta_f(1+\Delta_f)(R_f^{-1}-R_{f+1}^{-1}).\qquad\text{(5)}
```

Indeed `Σ_(j<G) (R+j)⁻² ≤ (1+1/R)(1/R−1/(R+G))`, by telescoping `1/[(R+j)(R+j+1)]`, and the integral approximation of the harmonic sum has error at most `1/R−1/(R+G)`. These estimates also cover `G=0` and are independent of how large a wait is. Sum (5) by dyadic blocks: on `f` of scale `2^k v` take the upper bound on `Δ_f(1+Δ_f)` at the block end and telescope the reciprocal increments. Both errors are at most `C_{p,ε}(1+H)^2/v`, uniformly in `n`, including `n→∞`.

Now `log(R_{f+1}/R_f)=log((f+1)/f)+e_{f+1}−e_f`. Summation by parts uses only the monotone increments of `Δ_f`. On each dyadic block the positive variation is at most its upper endpoint, so the errors from (4) are summable, with factor `2^{-k/8}`. The boundary terms obey the same bound. Thus the error of this second replacement is at most `C_{p,ε}(1+H)^2/√v`, uniformly over `n≥m`. Define the explicit path functional

```math
 d_f=\log(1+1/f)-1/(b+f),\quad
 W_{v,n}=\exp(-\textstyle\sum_{1\le f<n}\Delta_f d_f),\quad W_v=W_{v,\infty}.
```

Here `0≤d_f≤b/[f(b+f)]` since `b≥1`. With probability at least `1−3m/b−C/H²` the initial `m` blocks have `Δ_f=0` and all the preceding bounds apply; there `log L_n−log W_{v,n}` is uniformly small for every `n≥m`. In particular on those paths the finite-block log limit exists and gives the same estimate for `L_∞` versus `W_v`. For this conclusion only almost every `θ` in the lower sandwich is necessary. Taking `v→∞`, then `H→∞` and `ε↓0`, yields

```math
 \sup_{T,\;p\le\theta<1}
 \mathbb E_\theta|L_\infty\wedge1-W_v|=o(1)\qquad\text{(6)}
```

with the supremum for the almost-sure likelihood version understood in the essential sense. For fixed `M`, `n=⌊Mv⌋`, the analogous approximation holds in `L¹` for `L_n` without truncation: `W_{v,n}≤1`, and `L_n≤exp[n(1+J_n)/b]` has uniformly bounded second moment at fixed `M` (the Bernoulli parameters sum to `O_M(1)`). In particular no uniform bound on the untruncated infinite product is asserted. Finally (3) gives

```math
 \sup_{T,\theta}\mathbb E_\theta\sum_{f\ge Mv}\Delta_f d_f
 \le C(1+\log M)/M\quad(M\ge2).\qquad\text{(7)}
```

These estimates do not require a limiting birth environment. Jensen in (3) gives `liminf_v inf_(T,θ) E_θ W_v ≥ exp(−π_math²/3)`: for `f/v→s>0`, the harmonic sum tends to `log(1+s)` and `v d_f→1/[s(1+s)]`. The omitted small interval has upper bound `Cε` from `Σ_(h<f)(b+h)⁻¹≤f/b`; the upper tail is (7). Thus the upper limit of the mean cost is at most `2∫₀∞ log(1+s)/(s(1+s)) ds=π_math²/3`, using `t=log(1+s)`. Equations (2) and (6) prove the claimed universal `K_-`. The physical upper constant one follows from the locked urn at birth.

## 3. Local limit all the way to θ=1

An auxiliary continuous-time embedding makes the outside law local. At the start of block `f`, execute its prescribed tag choices without advancing this clock. Given the prescribed word and current state, draw a fresh `Exp(B_f+Δ_f)` wait independently of that step's outside winner, then take the winner. The successive rates are adapted. Tag choices change the tag degree, but not the total outside rate. With `T_f` the clock reading just after `f` outside jumps, the compensated clock has maximal squared deviation `O(1/v)` on `f≤Mv`, because the sum of the conditional waiting variances is at most `Σ_(h<Mv)(b+h)^(-2)`. Also, the expected difference between `Σ_(h<f)(B_h+Δ_h)^(-1)` and `Σ_(h<f)B_h^(-1)` is `O_M(1/v)` by (3), uniformly on that window. Therefore `T_⌊vs⌋→log(1+s)` in probability uniformly on `0≤s≤M`, including for sequences `θ_v→1`.

In this clock use independent physical drivers on the old inert side and the inert immigrant branches, and one independent rate-one tag Poisson clock. At its rings, the real tag uses its current degree and the exact physical fractions of leaf and core slots; the ideal process always starts a virgin. The independent genuine outside arrows inside each active side have their original rates. A boundary root's arrow to the tag is excluded from this clock, but the associated exterior slot stays in its degree. Thus until the first tag-to-core ring the old side and immigrants coincide with their ideal counterparts; instantaneous red choices only change the tag stock. On the first `εv` outside steps, the probability of a tag-origin ring is at most `εv/b+O(1/v)`. On the rest of an `Mv` window the tag degree is at least `pεv/2` with probability tending to one by the geometric count bound, and its core degree is at most `1+J_⌊Mv⌋`, tight uniformly. Conditional on the waits and past outside arrows, and on the pointwise degree bounds, the probability of tag-to-core at a step is at most the ratio of core degree to tag degree times `1/(b+f)`. Stopping at a fixed bound on the core degree, the sum of these probabilities is `O_{M,p,ε}(1/v)`; remove the stop using tightness. Let `ε↓0`. This proves the agreement in probability on each finite outside window, without coupling global discrete denominators in different trees.

For locally convergent initial old environments, the ideal finite-time paths at the boundary converge by the [continuous-time old-cell cutoff](../../../../../../../../astra/source/round2/operator.md#1-cells-and-evaluation-of-the-operator): retain the missing boundary slots and use ordered old-edge proposals with rates bounded by inverse initial core degrees. The bound for a path of length `h` by time `T` is `T^h/h!` times the inverse-degree product; sum by reversed simple random walk. The extra exterior slot at the old parent needs no incoming proposal. Couple the independent immigrations and their independent virgin continuations exactly. The auxiliary clock convergence and this local coupling give the claimed outside-clock limit of `Δ`. The first `εv` steps cause no singularity because the probability of any change of `Δ` there is `O(ε)`, independent of the tag waits.

For every sequence of leaf-birth environments converging locally to an old environment, (3), (7), and the finite-window convergence now show, for all choices of `θ_v≥p`,

```math
 \mathbb E_{\theta_v} W_v\longrightarrow
 g=\mathbb E\exp\left[-\int_0^\infty \frac{\delta(t)}{e^t-1}\,dt\right].\qquad\text{(8)}
```

Indeed `v d_⌊vs⌋→1/[s(1+s)]` uniformly on each clipped window `ε≤s≤M`, and the latter times `ds` becomes `dt/(e^t−1)`. The function `δ(t)` is the old-parent reciprocal deficit plus the deficits of the independent immigrants as in the reviewed proof. It is nondecreasing with tight total variation on bounded time windows. A uniform clock displacement of size `h` sandwiches its value between its translates by `±h`; on the clipped window their weighted integral differs by at most a constant times `h` times its variation on a slightly enlarged window. This justifies the weighted Riemann passage in probability. The resulting `g` is the same as in the interior proof. The expectation convergence is sequentially uniform over `θ∈[p,1)` for each convergent sequence of environments. Mix over the beta endpoint conditioned to exceed any `x_v≥p`, use (6) for the lower bound, and use fixed `M` then (7) for the upper bound. The physical tail divided by `Q_v(x_v)` tends to `g`, with no restriction on how fast `x_v` approaches one.

## 4. Annealed tail and one cutoff

For `α∈(0,1/2]`, `v≥3`, the beta tail at `x≥p` is bounded by

```math
 \frac{\Gamma(v)}{\Gamma(\alpha)\Gamma(v-\alpha+1)}
 (1-x)^{v-\alpha}
 \le\Pr(\mathop{\mathrm{Beta}}\nolimits (\alpha,v-\alpha)>x)
 \le\frac{\Gamma(v)}{\Gamma(\alpha)\Gamma(v-\alpha+1)}
 (1-x)^{v-\alpha}x^{\alpha-1}.\qquad\text{(9)}
```

For `α=1/2` the upper approximation has uniform relative error `O_p((1-x)/v)`: under the measure proportional to `h^{b-1}` on `0≤h≤1-x`, the mean of `(1-x)−h` is `(1-x)/(b+1)`, and `t^(-1/2)` is Lipschitz for `t≥p`. The gamma ratio gives `Q_v(x)∼x^(-1/2)v^(-1/2)(1-x)^(v-1/2)/√π_math` uniformly. For core placements `α=1/d≤1/3`; the beta upper is at most `C_p α v^(-2/3)(1-x)^(v-1/2)`. As in the interior argument, summing the original physical target weights uses `r_o/(D_o+1)≤1`, so this sector is `O_p(v^-1/6)` relative to the leaf scale even at the upper edge. For the leaf sector average the bounded score `(L_o/D_o)g`, using sequential uniformity in (8), extended portmanteau on the empirical core laws, and uniform integrability of `C_{v-1}/v^μ` from the existing local law. This yields (A) with the same `κ=E[X]π((L_o/D_o)g)/√π_math`. Positivity follows also from the universal comparison and `π(L_o/D_o)=μ>0`. Since `1-x≤1-p`, the uniform relative ratio implies one cutoff for every `x∈[p,1)`.

## 5. Two-limit interpretation

At *every* sequence `ε_v↓0`, even faster than any exponential in `v`, (A) gives

```math
 \Pr(Z_v>1-\varepsilon_v)\sim \kappa v^{\mu-3/2}\varepsilon_v^{v-1/2}.
```

Let `A_v=lim_(ε↓0) ε^(-(v-1/2)) Pr(Z_v>1−ε)`, which exists by the fixed-label prefactor theorem. Because the relative error in (A) is uniform, letting `x↑1` for each sufficiently large `v` gives `A_v∼κ v^{μ-3/2}`. Likewise along locally convergent leaf-birth seeds, the fixed-seed amplitude `Y_T` of the near-one theorem satisfies

```math
 \Gamma(v+1/2)\,\mathbb E[Y_T^{-(v-1/2)}]\longrightarrow g.
```

This follows by taking the fixed-`v` endpoint limit in the ratio to the beta survival before using the uniform relative convergence of §3. Thus the inverse moments normalize to the same constant in both orders.

### The departure and overshoot in the two-limit regime

Let `G` be the physical number of tag attachments before the first outside jump, so the physical pre-departure size is `K=v+G`. For any `ε_v↓0` the following limits are joint, conditional on `Z_v>1−ε_v`:

```math
 \left(\varepsilon_v(K-v),\; v\left[1-\frac{1-Z_v}{\varepsilon_v}\right]\right)
 \ \Longrightarrow\ \bigl(\mathrm{Exp}(1),\mathrm{Exp}(1)\bigr)
 \quad\text{with independent coordinates}.\qquad\text{(10)}
```

The centering is essential: `ε_v K` adds the deterministic shift `ε_v v`, which may tend to zero, a positive constant, or infinity. To prove (10), insert `H={G_0≥k}` into (2). Given `θ` and `H`, the first geometric wait is a fresh geometric plus `k` by memorylessness; all other waits are unchanged. Both dyadic bounds and the local limit in §§2–3 hold under this conditional law, uniformly for `uk` bounded, since the centered sum acquires only the bounded shift `uk`. The outside-index domination and the early-window bound are unchanged. The probability of the inserted event is exactly `θ^k`; the conditional likelihood has the same limit `g`. For `x=1−ε`, `k=ceil(z/ε)`, `z≥0`, all auxiliary `θ>x` obey `uk≤z+ε`. Under the conditional beta law, `u/ε→1` uniformly as `v→∞`: on `0≤u≤ε` the density is proportional to `u^{b-1}(1-u)^{-1/2}` and the probability of `u≤(1−η)ε` is at most `p^{-1/2}(1−η)^b` for `x≥p`. Hence the mixture factor `θ^k` tends to `e^{-z}`. The core sector with the inserted event is at most its original tail, so the weighted local averaging remains valid. For the joint survival with overshoot greater than `y≥0`, apply the same calculation with endpoint `x'=1−ε(1−y/v)` and **the same** `k`. There `uk` is still bounded, `u/ε→1`, and the ratio of physical tails at `x'` and `x` tends to `e^{-y}` by (A). This proves the product survival in (10). Rounding the offset by a single jump has no effect when `ε↓0`.

The inverse-tilted continuation amplitude in the fixed-seed theorem also has a large-label limit. Along any locally convergent sequence of leaf-birth seeds, the argument for the first coordinate of (10) holds at each seed, by (8) and positivity of `g`. At fixed `v`, letting `ε↓0` gives the law of `Y_v^* P_v`, where `Y_v^*` has the `Y^{-(v-1/2)}`-tilted amplitude law and `P_v` is its independent standard Pareto factor of index `v−1/2`, as in the near-one theorem. The fixed-`v` shift `εv` vanishes. A diagonal choice of `ε_v↓0` approximating these fixed-`v` distributions against bounded Lipschitz tests, together with the arbitrary-sequence first-coordinate limit, gives `Y_v^* P_v⇒Exp(1)`. Since `P_v→1` in probability, it follows that **`Y_v^*⇒Exp(1)`**. This also holds for the fixed-label mixture weighted by the physical leaf-placement prefactors, using the annealed version of (10). The conclusion is about the inverse-tilted amplitudes; the unweighted continuation laws are different.

## Bounded exact controls

```
python3 -IB round2/open/age-order/fixed-threshold/uniform/upper-edge/check.py
PASS 5580 exact stopped prefixes: block likelihood and reciprocal deficit
PASS 20 reciprocal-telescoping controls, including zero and long waits
PASS scaled geometric moments and prefix shifts at u=1/2, 1/10, 10^-24
```

The rational fixtures verify the physical likelihood, the deficit monotonicity, and the geometric scaling identities. The asymptotic claims use the estimates in the proof, not a numerical tail fit.
