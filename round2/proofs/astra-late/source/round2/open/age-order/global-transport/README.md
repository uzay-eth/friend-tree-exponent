> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../README.md#claim-26) for exact scope and citations.

# Finite-uniform pre-placement transport

**Result.** There is an absolute finite `Y` such that for **every** pre-placement tree of size `n≥2`, the next two births are strictly ordered at every `0<x<1` with `(n+1)x≥Y`. This does not use a typical-tree or leaf-density hypothesis, and allows thresholds arbitrarily close to one. The proof below also gives an annealed sufficient test on the remaining thresholds. We do not obtain its sign throughout the remaining range or a pre-placement counterexample.

## Signed comparison

Put `m=n+2`, `w_u=r_u/(D_u+1)`, `s_u=w_u/D_u`, `A=Σ w_u`. From the [exact two-birth identity](../adjacent/README.md), the numerator of the conditional gap, with denominator `n(n+1)`, is

```
N_T = Σ_u w_u (h_M,u - h_L,u) + Σ_{edges u,v} (s_v-s_u)(h_u|v - h_v|u).
```

Here `M,L` denote the serial middle and tip; `u|v` the parallel leaf placed at `u`. All endpoints are read at size `m` with the same threshold. Since `Σ_edges |s_v−s_u|≤A`, a uniform lower `B` on middle tails and uniform upper `U` on leaf tails imply `N_T≥A(B−2U)`. In every completion each leaf has a permanent core parent; take `U=Q_{m,1/2}(x)`, with `Q_{m,a}` the survival of `Beta(a,m−a)`, by the locked urn. Every serial middle, even on a bare old path, has a core parent and one fresh leaf.

## A uniform lower bound for the middle

In the serial completion use the [exact locked likelihood](../fixed-threshold/uniform/README.md) with `a=1+1/d∈[1,3/2]`, where `d≥2` is the starting degree of the middle's old neighbour, and `b=m−a`. Before outside jump `f`, write `G_f` for the number of red choices since the last outside jump, `R_f=a+Σ_(j<f)G_j`, `B_f=b+f`, and `Δ_f` for the physical attachment deficit relative to the locked red mass. Conditional on the beta endpoint `θ`, the `G_f` are iid geometrics with survival `θ^k`. Use the actual conditional outside-arrow law throughout. At outside horizon `N`, the log likelihood is

```
log L_N = Σ_(f<N) [log(1+Δ_f/B_f) + Σ_(j<G_f) log(1−Δ_f/(R_f+j))].
```

The deficit starts at zero and grows only at outside jumps. It is bounded by `1/2+J_f`, where `J_f` dominates the number of tag-origin arrows and is a sum of independent Bernoullis of means `1/(b+h)`, `0≤h<f`. This domination holds conditional on the whole red schedule, by proposing a tag-origin arrow with a fresh uniform. The contribution from the **old neighbour** alone, `Δ_f^old`, has conditional expectation at most `Σ_(h<f) 1/[3(b+h)]` even on a large-degree input: at old-neighbour degree `D` its gain weight is at most `D`, and the change to the deficit on a gain is `1/[D(D+1)]`. Do not require that the old neighbour never gains in the early window; that event has no useful uniform bound for large `D`.

**Initial blocks.** For `f≥1`, put `H_f=Σ_(j<G_f)1/(R_f+j)`. Conditional on `S=Σ_(j≤f)G_j=s`, the counts are uniform weak compositions of `s`. The probability that the last block includes the red index `j` is `binom(j+f,f)/binom(s+f,f)` for `0≤j<s`. Using `a≥1`, sum to get

```
E[H_f | S=s] ≤ Σ_(0≤j<s) binom(j+f,f)/[(1+j)binom(s+f,f)] ≤ 1/f.
```

The last inequality is the hockey-stick sum after dividing the binomial by `1+j`. Let `F=ceil(ε b)` with `0<ε<1/2`. With probability at least `1−F/b`, no tag-origin arrow occurs before block `F`. On that event `Δ_f=Δ_f^old≤1/2` for `f<F`. Dropping positive blue increments, the initial negative log is bounded there by `2 Σ_(1≤f<F)Δ_f^old H_f`; its expectation with the indicator of this event is at most `2F/(3b)`, by conditioning on all waits before applying the composition bound. At `f=0` the likelihood is one. Thus, by Markov's inequality, the probability that the initial log is below `−1` is at most `5F/(3b)` uniformly, even as `θ→1`.

**Remaining blocks.** Here are details of the bound used, to make the quantifiers on `θ` and tree shape explicit. For fixed `ε,H`, let `bθ` grow, take `u=1−θ` and `z_f=1+f/b`. With probability at least `1−C/H²` simultaneously on `f≥F`,

```
|uR_f−θf| ≤ 2+H sqrt(bθ) z_f^(3/4),   Δ_f ≤ 1+H z_f^(1/8).
```

For the first bound apply the maximal variance inequality on dyadic blocks starting at scale `b`, using `Var(uG_f)=θ`. For the second use the envelope's second moment `C(1+log(1+f/b))²` at the right endpoints, and sum the dyadic probabilities. On this event `R_f≥θf/(2u)` and `Δ_f/R_f≤1/2` for large `bθ`, uniformly in `θ`. Replace the log by `Σ Δ_f [1/B_f−log(R_(f+1)/R_f)]`: the blue error is at most `Σ (Δ_f/B_f)^2` and the red error plus integral approximation at most `C Σ Δ_f(1+Δ_f)(1/R_f−1/R_(f+1))`. On each dyadic block telescope the last factor after bounding the weight by its block maximum. Summing yields `O_(ε,H)(1/(bθ))` uniformly in the end time. Next write `log(R_(f+1)/R_f)=log(1+1/f)+e_(f+1)−e_f` where `|e_f|=|log(uR_f/(θf))|≤C_(ε,H)(bθ)^(-1/2) z_f^(-1/4)`. Monotonicity of `Δ` gives a summation-by-parts error `O_(ε,H)((bθ)^(-1/2))`, including the boundary terms: its variation on each dyadic block is bounded by its endpoint, and the block series has factor `2^(-1/8)`. This is the outside-block cancellation of the [upper-edge proof](../fixed-threshold/uniform/upper-edge/README.md), with the actual variance `θ` retained.

The replacement log is `−Σ Δ_f d_f`, where `0≤d_f=log(1+1/f)−1/(b+f)≤b/[f(b+f)]`. Conditional on the waits, `EΔ_f≤(4/3) Σ_(h<f)1/(b+h)` by the old-neighbour and tag-origin bounds. Accordingly `EΣ_(f≥1) Δ_f d_f≤C_0` for an absolute constant: sum for `f≤b` using `EΔ_f≤(4/3) f/b`, then group dyadically for `f>b`, with series `O(Σ (1+j)2^(-j))`. No local convergence of the old tree is required.

Choose `ε` small and `H` large so the two bad-event probabilities above sum to at most `1/4` for large `b`; then use the event `Σ Δ_f d_f≤4C_0` (probability at least `3/4`). For `bθ` sufficiently large the remaining-block error is at most one. With probability at least `1/2` all finite likelihoods with `N≥F` therefore exceed `exp(−2−4C_0)`. The infinite likelihood converges almost surely under the unconditioned reference, hence for beta-almost every `θ`; the outside stopping times tend to infinity on those paths. The [dominated likelihood measure](../fixed-threshold/uniform/README.md) gives the uniform bound

```math
 h_M(x)\ge c\,Q_{m,a}(x)\ge c(1-x)^{m-1}\qquad (bx\ge Y_0),
\quad c=\tfrac12\exp(-2-4C_0)>0.
```

In the last inequality beta survival increases with `a` at fixed total mass (the density ratio is increasing). The proof only uses the lower likelihood sandwich, so does not assume uniform integrability of the infinite product. Finally `Q_{m,1/2}(x)≤C(mx)^(-1/2)(1-x)^{m-1/2}` by replacing `t^(-1/2)` by `x^(-1/2)` in the integral and bounding the gamma ratio. Hence `2U/B≤(2C/c)(mx)^(-1/2)`. Increase `Y_0` to compensate for `b` versus `m−1` and to make this ratio less than one. This proves the stated conditional comparison. The constants are not optimized; no numerical value of `Y` is asserted.

The distinction from the [effective annealed route](../fixed-threshold/effective-uniform/README.md) is that **all** serial coefficients buy the lower bound, not just those from first attachment to an old leaf. The initial-block reciprocal estimate is what keeps this valid when the middle's old neighbour has unbounded starting degree. In return this comparison only signs the upper range and does not provide a normalized adjacent equivalent.

## Lower thresholds: an annealed first-promotion test

One can retain the exact edge coefficients without demanding pointwise order of every transfer. For each edge with `s_v>s_u`, the only leaf needing an upper bound is `v|u`; the positive `u|v` may be retained if a lower bound is available. Likewise if `s_u>s_v` use `u|v`. For any such leaf (or serial tip), call its parent's post-placement degree `d` and retain its `j` other leaf neighbours. Until its first promotion the parent receives rate-one return immigration from the tag, and the tracked stock increases from its own leaf births at rate equal to its stock. If the first promotion is at size `k`, the physical continuation tail is at most `Q_{k,3/2}(x)`.

For a short common elapsed time `t`, the leaf upper bound can be taken as

```math
 U_t=\min\{Q_{m,1/2}(x),\ \mathbb E[1_{\tau\le t}Q_{N(\tau),3/2}(x)+1_{\tau>t}F(D(t),J(t))Q_{N(t)+1,3/2}(x)]\},
```

where `τ` is its first promotion, `N` the actual total size, and `F(D,J)` any valid no-age first-send upper bound such as `min(1,2/(J+1))` from the [leaf barrier](../early/README.md). `J` may use only tracked other leaves. This follows by conditioning at the earlier of `τ,t`; for a pending leaf the eventual send probability is bounded by `F` and the beta cap decreases with the promotion size. The product in the expectation keeps age, size and parent stock together. It does **not** multiply an averaged stock payoff by an independent size-discount. At `t=0` it already allows a cheap exact degree/stock test without exposing a single new tree. For example, if `B_{M,u}` is any middle lower bound, a sufficient test at deterministic old size is

```
E Σ_u w_u B_M,u > E [Σ_u w_u U_tip,u + Σ_edges |s_u−s_v| U_adverse(u,v)].
```

This is an inequality of expectations over the actual old-tree law. It allows negative rows on individual trees and uses only one adverse endpoint per ordered coefficient. At zero the beta factors equal one; then small-stock adverse endpoints remain expensive and this test has not been shown positive for every old size. In particular neither a negative oriented transfer nor a post-placement reversal decides the whole-tree criterion. The upper-range theorem reduces any counterexample at large size to a bounded birth-scale threshold, but the known annealed positivity there is eventual and does not fill the finite unresolved label interval.

There is a no-enumeration version of the early-horizon upper charge. For a pending leaf starting at `(m,d,j)`, let `F=min(1,2/(j+1))`. The cumulative probability of promotion by size `m+k` is bounded by

```math
 C_k=\min\left\{F,\ 1-\prod_{i=0}^{k-1}\left(1-\frac1{(m+i)d}\right)\right\},\qquad C_0=0.
```

Indeed, conditional on all previous sizes and no send, its next-jump promotion hazard is at most `1/[(m+i)d]`, since the parent degree cannot decrease. The lifetime constraint `F` is independent of this hazard bound. With `q_k=Q_{m+k,3/2}(x)`, for any small integer `K≥1` a valid upper endpoint is

```
min(Q_(m,1/2)(x), Σ_(k=1)^K (C_k−C_(k−1)) q_k + (F−C_K) q_(K+1)).
```

This follows by summation by parts against the decreasing sequence of beta tails: assign the maximal allowed cumulative send mass to the earliest sizes first. It does not require those maxima to be jointly attainable in the tree. It improves (strictly at `0<x<1` if `C_1<F`) on the undiscounted early charge `F q_1` even at `K=1`, costs only scalar operations per placement, and can be combined by taking the minimum with the elapsed-time bound. At `x=0` all the `q_k` are one and it gives exactly `F`, identifying the limitation of age discount alone at hubhood.

## Replay

`python3 -IB round2/open/age-order/global-transport/check.py` checks the initial-block composition bound and the signed upper charge on small rational fixtures. It does not estimate `Y` or enumerate physical early births.
