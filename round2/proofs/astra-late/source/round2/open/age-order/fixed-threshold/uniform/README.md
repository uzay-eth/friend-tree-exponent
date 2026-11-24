> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports. See the [release index](../../../../../../README.md#claim-24) for exact scope and citations.

# Fixed interior thresholds: physical prefactor and adjacent sign

For the original edge seed, let `Z_v` be the limiting degree share of label `v`, and let `μ`, `X`, `π` be the core exponent, the limit of `C_n/n^μ`, and the limiting uniform-core environment. Using the [existing deterministic-size local law and uniform integrability](../../late/README.md), we obtain, uniformly on every compact subinterval of `0<x<1`,

```math
 \Pr(Z_v>x)\sim \frac{\kappa}{\sqrt{x}}\,v^{\mu-3/2}(1-x)^{v-1/2},
 \qquad
 \frac{\Pr(Z_{v+1}>x)}{\Pr(Z_v>x)}\longrightarrow 1-x .
```

Here `κ>0` is given below and is independent of `x`. This proves eventual strict adjacent order at every fixed interior threshold (uniformly on compact subintervals). It does not give an effective cutoff or order on the entire open interval at once. **Credit:** Fable Claim 5(c) gives finite-label full stochastic domination by alternating rooted transport (two receipts through 19→20 at `88e13f3`, and a later single receipt for 20→21; see the [audit](../../../../attacks/fable5/flow/README.md)). Claim 5(e) supplies the Pólya upper comparison and fixed-label endpoint ordering. The [fixed-seed near-one constant](../../near-one-constant/README.md) uses the exact first-defect hazard; here we transport the reference conditioned on an interior endpoint and control the physical likelihood along that entire path.

## 1. Exact reference and physical likelihood

At the marked birth, let `d≥2` be the birth parent's post-placement degree and set `a=1/d`, `b=v−a`, `A₀=1−a`. In any continuation with total size `j`, marked degree `D`, and attachment rate `r=Σ_(y~tag) 1/D_y`, put

```
R=D−A₀,       B=j−R,       Δ=R−r=A−A₀≥0,
A=Σ_(y~tag)(1−1/D_y).
```

The reference chooses the tag with probability `R/j`, and outside with probability `B/j`; after an outside choice use the **physical** conditional outside arrow law (total weight `B+Δ`). The tag choice attaches a leaf to the tag, regardless of the winner. Under the reference the tag indicators are a fractional Pólya sequence with endpoint `Θ~Beta(a,b)`. Given `Θ=θ`, their values are iid Bernoulli(`θ`); choose the outside arrows as above with independent fresh drivers. This disintegration follows by integrating any finite Bernoulli word against the beta density. Define the factors from state just before each jump by

```math
 \ell_{\rm tag}=1-\Delta/R,\qquad \ell_{\rm out}=1+\Delta/B .
```

For every finite history, their product `L` is the exact physical likelihood. If `F_i` is the number of outside choices before step `i` (`i=0` at birth), then `B_i=b+F_i`. The number `C_i` of outside arrows initiated by the tag bounds `Δ_i≤1+C_i`: changes in the old parent's reciprocal cost at most one, and each other core neighbour has been created by one tag initiation. On the full physical law the fraction has an almost-sure limit: `D/j` is a nonnegative supermartingale. Under the reference it tends to `Θ`.

Here and below finite lists of outside arrows retain all multiplicities. The likelihood does not replace the physical outside winner by an unweighted exterior choice. The nonnegative finite-history likelihood is a martingale under the unconditioned reference and converges there almost surely to `L_∞`. The measure with density `L_∞` is dominated by the physical law: on each finite-history cylinder use Fatou, and extend the inequality to the generated σ-algebra by the monotone class theorem. Both endpoints are the same path functional (limiting degree share). In particular, for any prefix event `H`,

```math
 \Pr_T(H,\ Z>x)\ge \int_x^1 f_{a,b}(\theta)\,\mathbb E_\theta[1_H L_\infty]\,d\theta. \qquad\text{(1)}
```

Equality of the infinite-path measures is not needed.

## 2. Global localization of the change of measure

Number the outside steps by `f=0,1,…`, so `B=b+f`. Given the Bernoulli word, the conditional chance that the outside arrow is from the tag is at most `1/(b+f)`; the tag's total initiation rate is one. Couple all these events below independent Bernoulli variables with these parameters. This uses a fresh uniform for the tag-origin decision even though the underlying state changes. Write `J_f` for the sum of the dominating variables before outside step `f`.

Fix `0<p<q<1`. For every fixed `k`, the independent Bernoulli envelope gives `E[J_f^k]≤C_k(1+log(1+f/v))^k`. Below we use dyadic ranges of the jump index, with `F_i≤i`, to get a simultaneous bound on `C_i`. For the likelihood one must keep its signed drift: separately bounding the sums of tag and outside log factors on an infinite window would lose their cancellation.

Let `S_i=R_i−a` be the tag count and define `z_i=R_i/(v+i)`. Conditional on a state, for the reference with fixed parameter `θ`,

```math
 \mathbb E_\theta[\ell_i-1\mid\mathcal F_i]
 =\Delta_i\left(\frac{1-\theta}{B_i}-\frac{\theta}{R_i}\right)
 =\frac{\Delta_i(z_i-\theta)}{(v+i) z_i(1-z_i)} . \qquad\text{(2)}
```

For rigorous localization after `i≥Mv`, use the stopping time at which any of the following bounds fails (fixed `H≥1`, `i≥Mv`, `M≥4/p`):

```
|S_i−θ i| ≤ H sqrt(v) (1+i/v)^(3/4),
C_i ≤ H(1+i/v)^(1/8),
z_i ≥ p/4,       1−z_i ≥ (1−q)/4.
```

Under the reference the probability of stopping is bounded by `C/H²`, uniformly in `v,M,θ`, for large `v` at fixed `H`. The first bound uses the maximal variance inequality on dyadic index blocks up to `2^k v` and summable `2^(−k/2)`. The second uses the outside-index moment bound with exponent `1/8` (since `F_i≤i`); the last two follow from the first at large `v`. On the stopped path, the sums of the absolute conditional drifts in (2), and of the conditional squares `E_θ[(ℓ_i−1)^2|F_i]`, starting at `Mv`, are bounded respectively by

```math
 C_{p,q}\left(H M^{-7/8}+H^2 v^{-1/2}M^{-1/8}\right),
 \qquad C_{p,q} H^2 v^{-1}M^{-3/4}. \qquad\text{(3)}
```

Indeed `|z_i−θ|≤(v+H√v(1+i/v)^{3/4})/(v+i)` and `Δ_i≤1+C_i`. On a stopped path the variance and drift bounds also control log likelihood (for large `v` each factor differs from one by less than `1/2`). Thus the log increment after `Mv` vanishes in probability in the iterated limit `v→∞`, `M→∞`, then `H→∞`. This controls the tail in probability. The following separate device controls expectations without assuming any moment of the exponential of a dominating product.

### Upper comparison after a finite prefix

Under the unconditioned reference, conditional on a finite history, its endpoint has law `Beta(R_i,B_i)`. The **physical** future conditional tail is bounded above by that same beta law, by coupling future red counts with the locked initial deficit `A₀`. In particular, for any event `H` measurable by step `i`,

```math
 \Pr_T(H,\ Z>x)\le \mathbb E_{\rm ref}[1_H L_i;\Theta>x]. \qquad\text{(4)}
```

Equations (1) and (4) give opposite sides without an infinite-time uniform-integrability assumption. All finite factors are positive because the tag's actual attachment rate is positive. On the tail without a stop the square and drift estimates give convergence of its log; equivalently one can first bound all finite log increments after `Mv` by the same stopped martingale maximal bound, then use the almost-sure limit `L_∞`. The conditional reference is only needed for almost every `θ` under the beta mixture in (1); the finite-horizon and stopped bounds are uniform in the compact interval.

## 3. Finite windows at time `t=log(1+(1−θ)i/v)`

For the leading sector, the mark is born onto a leaf; `a=1/2`. Call that promoted parent `u`, and remove the tag with its future descendants at the old edge to `u`, keeping a degree slot at `u`. In the limit the old side runs the inert physical dynamics in the above time `t`. Define `P(t)` to be an independent rate-one Poisson process. Its arrivals start fresh independent inert virgin branches, each rooted at a former tag leaf with its own exterior degree slot. In this limit

```math
 \delta(t)=\tfrac12-1/D_u(t)+\sum_{s\in P,\;s\le t}(1-1/D_s(t)),
 \quad
 \mathcal L=\exp\left[-\int_0^\infty \frac{\delta(t)}{e^t-1}\,dt\right]. \qquad\text{(5)}
```

The integral is finite almost surely: `δ` is zero until the first change of `D_u` or the first immigration, and `δ≤1+P(t)`. The old side boundary degree starts at 2 and has gain intensity at most 2 before it changes, so `E[δ(t)]≤2t` (the old term costs at most `1/2`, and the immigration terms at most `P(t)`). By Jensen, uniformly in the initial environment, `exp(−π_math²/3)≤E[ℒ]≤1`, since `∫₀∞ t/(e^t−1)dt=π_math²/6`. This is a potential along the actual limiting outside shape process, not a formal substitution into the fixed-seed inverse moment.

The constant can be computed as a product of two expectations. Write `k(t)=1/(e^t−1)`, and let `D_vir(a)` be the boundary degree of an inert virgin at age `a`. Then the Poisson exponential formula gives

```math
 g=\mathbb E_{\rm old}\exp\!\left[-\int_0^\infty (\tfrac12-1/D_u(t))k(t)\,dt\right]
 \exp\!\left[-\int_0^\infty
 \left(1-\mathbb E_{\rm vir}\exp\!\left[-\int_0^\infty (1-1/D_{\rm vir}(a))k(s+a)\,da\right]\right)ds\right].
```

The second factor is universal; its exponent is finite by `1−e^(−z)≤z` and `∫₀∞ t k(t)dt<∞`. No distribution of a typical (unweighted) birth parent is substituted for the old-side law.

Here are details of convergence, including the singularity at the start. Along any sequence of leaf-placement environments converging locally, and any `θ_v∈[p,q]` tending to `θ`, take `i=⌊vs⌋`. Conditional on `θ_v`, uniformly on finite `s` intervals, `R_i/v→θ s`, `B_i/v→1+(1−θ)s` in probability. Outside conditional arrow probabilities are the physical weights divided by `B_i+Δ_i`, with `Δ_i=O_P(1)` on finite `s` intervals by the outside-index bound. In `s` time each fixed local arrow of weight `w` therefore has limiting intensity `(1−θ)w/[1+(1−θ)s]`, which is just `w` in `t` time. The initial old side has no incoming tag arrow in the limit: on `s≥ε>0` there are `O_P(1)` core neighbours of the tag and its degree is of order `v`, so a tag ring goes to its leaf stock with probability tending to one. The total limiting tag-outside clock has rate one in `t` and hence starts the independent virgin branches. Their internal rates and the old-side rates retain exterior slots. No reference tag gain alters these branches.

For explicit justification of the local process convergence, first keep a fixed ball in the initial old core, with boundary slots, and stop its cells and the tag immigration cluster at a fixed population bound. Given the Bernoulli schedule, use independent uniform proposal slots for the finite bounded list of local jump channels. The factor `1/(B_i+Δ_i)` differs from `1/B_i` by `O(v^-2)` on the stopped window, so replacing it in these local slots costs `O(v^-1)` in probability. The remaining proposal can be assigned to outside changes that do not enter the tracked ball. Now the Bernoulli outside counting measure divided by `v` converges to `(1−θ)ds`, and the one-step no-local-jump products converge to the exponentials of the limiting hazards. Iteration over the bounded jump types gives the joint clock law; incoming tag-to-core arrows are negligible away from zero, and near zero all tag-origin arrows cost `O(ε)`. Remove the population stop using Yule and finite immigration domination (conditional expected increases are bounded by current tracked population plus a finite boundary rate, divided by `b`). For removal of the spatial cutoff, couple two outside evolutions with different boundary treatments by using the common local proposal denominators `B_i`: the cost of substituting these denominators on a fixed stopped ball is `O(v^-1)`, and use underlying proposal slots bounded by `1/(b c_z^0)` on old edges. A specified ordered path of length `h` has probability at most `(i/b)^h/h!` times the product of `1/c_z^0` for its senders in an `i`-step window. Sum by the reversed-walk bound of the [old-cell construction](../../../../../../../astra/source/round2/operator.md#1-cells-and-evaluation-of-the-operator). The same bound covers proposals entering a fixed ball, uniformly in the approximating environments. This justifies the finite-window joint convergence of `Δ` and its jump times. In the first `εv` steps, the chance of any change to `Δ` is `O(ε)`, uniformly in all leaf-placement environments: while unchanged, the tag-origin rate is at most one and the parent gain rate is at most two, so each outside step has chance at most `3/b` to change it. The corresponding limiting bound is the same order. Let `ε↓0` after taking the finite-window limit.

On `ε≤s≤M`, expand the logs of the factors and use the martingale difference bound for the Bernoulli draws; the quadratic error is `o_P(1)` after stopping at a fixed bound for `Δ`. The sum converges to

```math
 \int_{\varepsilon}^M \delta(t(s))
 \left(\frac{1-\theta}{1+(1-\theta)s}-\frac1s\right)ds
 =-\int_{t(\varepsilon)}^{t(M)} \frac{\delta(t)}{e^t-1}\,dt .
```

Consequently `L_⌊Mv⌋` converges in law to the truncation of (5). It also converges in expectation for each fixed `M`: on the first `Mv` steps the product is at most `exp(M(1+J_⌊Mv⌋) v/b)`, by discarding the tag factors. The dominating `J` has all fixed exponential moments bounded uniformly in `v` at fixed `M`. This gives uniform integrability for the finite-window calculation.

Combining the two sides avoids any infinite-time uniform-integrability assumption. For the upper bound use (4) at `i=⌊Mv⌋` and the finite-window expectations. For the lower bound use (1), (3), and bounded continuous truncations of `L_∞`: first send `v→∞`, then `M→∞`, then `H→∞` in the stopping bound. Formula (5) is at most one; the finite-window limiting likelihood decreases to it, allowing bounded convergence. The beta endpoint conditioned on `Θ>x_v` converges in probability to `x`; outside a compact neighbourhood of `x` its conditional mass is exponentially small. The finite-window expectation bound is uniform even on that discarded portion, since its outside-index domination does not depend on `θ`. These arguments give, along every locally convergent sequence of such birth environments and every `x_v→x∈(0,1)`,

```math
 \frac{\Pr_T(Z>x_v)}{\Pr(\mathrm{Beta}(1/2,v-1/2)>x_v)}\to \mathbb E_{\rm old}[\mathcal L]=:g(\text{old environment}) . \qquad\text{(6)}
```

In particular the likelihood argument is a transport under the actual rare event. The first outside jump offset after birth is tight, rather than of order `v`. The exact reference word gives `Pr_θ(first outside offset≥k)=θ^k`. Conditional on the first `k` choices all being tag choices, the remaining word is still iid and no outside object has changed, while all finite-window counts start at offsets at most `k`. The same finite-window and stopped-likelihood arguments with that prefix event inserted in (1) and (4) give the geometric survival limit `x^k` under the physical conditional law. Section 4 transfers this to the unconditional birth law as well. The fixed-label near-one Palm scaling `εK` cannot be applied with `ε` fixed as `v→∞`; the offset, not the absolute pre-size, is the informative departure variable here.

There is also a nonasymptotic bound on the discarded offsets. On a perfect prefix of length `k`, `Δ=0` and `L_k=1` exactly. Using (4), for any `x<q<1` the joint physical tail with this prefix is at most `q^k Pr(Θ>x)+Pr(Θ>q)`. Averaging over birth states and using the estimates in Section 4 bounds the conditional offset tail, for `x` in a fixed compact interval below `q`, by `C(q^k+e^{-cv})`, uniformly in `k` for all sufficiently large `v`. This uses the actual outside hazard at every prefix.

## 4. Negligible core-parent placements and summation

If the birth is onto a core, `a=1/d≤1/3`, so its locked beta prefactor is smaller. A uniform bound retaining the rate factor `a` is essential: total core-target weight for a single tree can be large. For `x` in a compact interior interval and `0<a≤1/3`, integrating the density and using `1/Γ(a)≤C a` gives

```math
 \Pr(\mathrm{Beta}(a,v-a)>x)\le C_{p,q}\,a v^{-2/3}(1-x)^{v-1/2}.
```

Since `r_o/(D_o+1)≤1`, average over core targets to bound their contribution by `O(v^{μ-5/3}(1-x)^{v-1/2})`, smaller by `v^(-1/6)` than the leaf sector. For leaf targets the birth initiator is their unique core neighbour `o`; conditional on `T_{v-1}`, their combined weight at that core is `L_o/[(v-1)D_o]`. For the old environment in (5), replace one of those leaves by the degree-two parent `u`, with one edge to `o` and one retained exterior slot (the edge to the removed tag). Equal leaf choices are isomorphic. The per-core score before birth is therefore `(L_o/D_o)g`. It is bounded by one and continuous in local topology: finite-window cutoff as above, followed by `E∫_T^∞ (1+P(t))/(e^t−1)dt→0`; at `L_o=0` give score zero. For convergent input environments the sequential asymptotic (6) applies, and every normalized finite score is bounded by one by the locked cap. Extended portmanteau on the empirical core laws, followed by uniform integrability of `C_{v-1}/v^μ`, gives the limit mean score `E[X] π((L_o/D_o)g)`. In particular `exp(−π_math²/3) μ E[X]≤E[X] π((L_o/D_o)g)≤μ E[X]`.

Finally the locked beta survival is uniformly asymptotic to `x^(-1/2)v^(-1/2)(1-x)^(v-1/2)/sqrt(π_math)` on compact interior intervals (endpoint integration and the gamma ratio), with

```math
 \kappa=\frac{\mathbb E[X]}{\sqrt{\pi_{\rm math}}}\,\pi\big((L_o/D_o)g\big).
```

The leaf-birth weight has order `v^(μ−1)` and the beta survival contributes `v^(−1/2)`, giving `v^(μ−3/2)`. The uniform adjacent ratio follows by taking the quotient at the same threshold. The same quotient at thresholds `x+y/v` and `x` gives `Pr(v(Z_v−x)>y | Z_v>x)→exp(−y/(1−x))` for `y≥0`. Inserting the fixed outside-offset prefix in both thresholds gives the joint limit: this exponential overshoot is independent of the geometric first-outside offset. Neither the finite-window comparison nor the ratio here is uniform up to either endpoint.

## Exact controls

```
python3 -IB round2/open/age-order/fixed-threshold/uniform/check.py
PASS 1968 physical/ref arrow descendants and all depth-4 Pólya words
PASS 72 signed drift and clock identities
```

The fixtures check the physical/reference likelihood at every prefix, with leaf and core birth parents, and the sign and time change of the drift. They do not estimate the prefactor or establish the limit numerically.
