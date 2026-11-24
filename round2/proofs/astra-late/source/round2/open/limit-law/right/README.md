> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../README.md#claim-17) for exact scope and citations.

# The right-tail rate exists: a finite-block variational formula

Put `a=1−μ`, `p=1/a`, and let `Φ` be the continuous-time eigenfunction in actual core-count normalization. For the edge-seeded fixed-size limit `X`,

```text
I = lim_(x→∞) [−log P(X>x)] / x^p ∈ (0,∞).
```

The following characterization uses only the finite embedded chain and the deterministic eigenfunction. Write `H_m(b)=P(Φ(T_m)>bm)` and use `−log 0=∞`. For any fixed `B` with `Φ(T_m)≤Bm` for all `m`,

```text
I = lim_(δ↓0) liminf_(m→∞) inf_(δ≤b≤B) [−log H_m(b)]/(m b^p).
```

In particular, along block sequences with `−log H_m(b)≤Dm` for a fixed `D` and `b≥δ>0`, repetition costs `−log H_m(b)+o(m)` per block and preserves the eigenfunction density up to a relative `o(1)` loss. The linear-cost hypothesis makes the root-source trimming negligible relative to the rare block probability; the density floor makes the bounded joining loss negligible relative to `bm` and the overhead negligible on the `m b^p` rate scale. No such uniform rate comparison is asserted for superlinear costs or for densities tending to zero. In the variational proof `δ` is held fixed through the block limit before being sent to zero. Neither a numeric value nor the existence of an optimizer at a fixed density is asserted.

Inputs: the [uniform finite-time free one-leaf gradient, embedded martingale, and right-tail order](../README.md). This proof does not use the degree-sensitive gradient or the left-tail prefactor.

## 1. Deterministic cost of joining free components

Take a finite scaffold and `k` disjoint free trees of size at least two, and join each tree to the scaffold by one edge. Colors are assigned to these `k+1` components; all future vertices inherit the color of their attachment target, so only the original edges cross colors. Delete those edges to view the evolving components as free trees. There is a uniform `G` such that, on every deterministic initial configuration of this form,

```text
|Φ(joined) − Σ_components Φ(free component)| ≤ G k.
```

To see this, let `f_u` on colored states be the sum of the disconnected free means of the core counts after duration `u`. The disconnected generator governs its backward equation. At a physical source with internal degree `d≥1` and `c` cross edges, the total deficit of internal-arrow rates compared to the disconnected generator is `c/(d+c)`; the total extra cross rate is the same. Sum over the endpoints to get a total absolute rate difference at most `4k`, including both deficits and extra arrows. Each arrow adds a leaf in exactly one component, changing `f_u` by at most `L_f e^(μu)` by the finite-time free gradient. Apply the time-dependent stopped jump formula to `f_(t−s)` along the joined process, then remove stops (finite-time size has a Yule second-moment bound and the residual integrand is uniformly bounded). The difference of the two expected free-component count sums at time `t` is bounded by `4k L_f ∫₀ᵗ e^(μ(t−s)) ds`. The physical count differs from the sum over cut components by at most `2k`, since only original cut endpoints can differ in core status. Divide by `e^(μt)` and let `t→∞`, giving `G=4L_f/μ`. This estimate does not require branch independence under the physical law or a sign for the effect of a join.

## 2. Scaffold and simultaneous repetition under the embedded law

For `k≥1`, start by growing a degree-three-capped scaffold to size `3k`. The [capped construction](../README.md#3-matching-lower-right-tail-rate-keep-degree-at-most-three) costs at most `12^(3k−3)` in reciprocal probability. Choose `k` distinct scaffold anchors measurably. In a first stage give each anchor a new leaf once, in any order, using any of its original scaffold neighbors as source (choose one per target). Source degrees during the stage are at most four. With `j` unused anchors the chance of such a next arrow is at least `j/(16k)`, since the size is at most `4k` and the arrows for distinct targets are disjoint. In a second stage promote each new leaf once by its anchor in any order, with probability at least `j/(20k)` when `j` such leaves remain. Thus the full setup costs at most `exp(c_0 k)` in reciprocal probability (use `k!≥(k/e)^k`). There are now `k` identifiable two-vertex blocks each with one external root slot, and an unchanged scaffold of size `3k` after removing the connecting edges. All these estimates refer to genuine steps of the free embedded chain.

Compare ordered directed-arrow histories after this setup to an auxiliary disconnected urn. At every step choose a color with probability proportional to its current ordinary size, and within that color sample the free directed arrow; a birth adds one to that color. Thus the urn schedule is independent of the prepared independent block histories, because sizes depend only on the schedule. Demand no scaffold draw for `L=k(m−2)` steps and exactly `m−2` births per block. The no-scaffold probability is

```text
∏_(i=0)^(L−1) (2k+i)/(5k+i) ≥ exp[−c_1 k(1+log m)].
```

This follows from `log(1+3k/(2k+i))≤3k/(2k+i)`. Conditional on no scaffold, the schedule has the ordinary `k`-color Pólya law with initial weight two each, since the displayed factor depends only on the total steps. The weight of a composition `(l_1,…,l_k)` with sum `L` is proportional to `∏(l_i+1)` by summing the exchangeable sequences. The balanced composition maximizes it by the arithmetic–geometric mean inequality. Its probability is at least the reciprocal of the number of compositions, itself at most `(em)^k` (or use the binomial coefficient `binom(L+k−1,k−1)`). This is a balance event at a **fixed total size**, with no trimming of oversized blocks.

For an internal arrow, the ratio of its physical to urn probability is one unless the source is a block root, in which case the ratio is `d/(d+1)≥1/2`, with `d` its current free degree. On the above schedule, no arrow touching the scaffold is used. In each free block let `R_m` count steps whose source is its specified original root (including the step from size two). Conditional on the past, the source is uniform at size `j`, so `R_m` is a sum of independent Bernoulli variables of probabilities `1/j`, `2≤j<m`. For `s_m=ceil(m/sqrt(log m))` and any fixed `D`, `P(R_m>s_m)≤exp(−Dm)` for all sufficiently large `m`, by the Chernoff bound `(e∑1/j /s_m)^s_m`. This is an **unconditional** trimming bound; root-source frequencies may be strongly biased inside the rare block event. In particular if `−log H_m(b)≤Dm`, trimming at `s_m` retains at least `H_m(b)/2` for all sufficiently large `m` (use `D+1` in the preceding estimate).

Condition on any successful setup. Summing the likelihood comparison over the prescribed schedules and trimmed free histories gives a physical event at size `n=k(m+3)` with probability at least

```text
exp(−k J_m(b)),   J_m(b)=−log H_m(b) + c_2(1+log m) + s_m log 2,
```

where the universal `c_2` also pays for setup and the factors of two per block. Histories here include source/target choices, so distinct ordered histories are disjoint even if the final shapes coincide. On this event, nonnegativity of the scaffold amplitude and §1 give `Φ(T_n)≥k(bm−G)`. No time-change, independent survival mask, or source-bias assumption is used.

## 3. From one good block to the entire right tail

Define `f(x)=−log P(X>x)` and `I_- = liminf f(x)/x^p`, `I_+ = limsup f(x)/x^p`. Both are positive and finite by the earlier order bound. For any sufficiently large fixed `m,b` to which the trimming bound in §2 applies and with `bm>G`, write `b_*=(bm−G)/(m+3)`. Given any `ε>0`, for each `x` choose the least integer `k` with `b_*[k(m+3)]^a≥(1+ε)x`. The embedded martingale mean at that size is `Y_n=γ_n Φ(T_n)≥(1+ε/2)x` for large `x`. By conditional variance, the probability of `X>x` on each such history tends to one uniformly as `k→∞`. Thus

```text
I_+ ≤ (1+ε)^p J_m(b)/[(m+3)b_*^p].
```

Let `ε↓0` here. For every fixed `δ>0` and any sequence with `δ≤b≤B` and bounded block quotients, the trimming hypothesis holds with a fixed `D`; the overhead is `o(m)` and `b_*/b→1`. This proves `I_+≤V(δ)` with `V(δ)` the inner liminf/inf in the formula (trivial if infinite). In particular every finite block used here is first fixed while `x→∞`; sending `m→∞` happens afterward.

For the other inequality, fix `0<ε<1/4`, and take a sequence `x→∞` attaining `I_-`. At the deterministic size `m=ceil(Mx^p)`, if `Φ(T_m)≤(1−2ε)x m^μ`, then `Y_m≤(1−ε)x` for large `x`. The uniform one-sided conditional concentration of the embedded limit bounds the contribution of this event to the tail by

```text
exp(−c ε² x² m^(2μ−1)) ≤ exp(−c ε² M^(2μ−1) x^p),
```

using `2+p(2μ−1)=p`. Choose `M` finite with this exponent constant strictly greater than `I_-`. The rest of the tail is covered by `H_m(b)` at `b=(1−2ε)x/m^a→(1−2ε)M^(−a)>0`. Along this subsequence `H_m(b)≥P(X>x)(1−o(1))`. Therefore `V(δ)≤I_-/(1−2ε)^p` for any `δ` smaller than the limiting `b`. Taking `δ↓0`, then `ε↓0`, gives `I_+≤lim_(δ↓0)V(δ)≤I_-`, as required.

## 4. Conditional growth localization and where clock bias enters

The extraction argument works along the full tail once the rate exists. More explicitly, for `m=ceil(Mx^p)` and every fixed `ε>0`, choosing `M` large enough yields

```text
P(|X−Y_m|>εx | X>x) ≤ exp(−c_ε x^p)
```

for large `x`, by the two-sided unconditional estimate obtained from conditional concentration and division by `P(X>x)`. The rate limit also gives `P(X>(1+ε)x | X>x)=exp(−[I((1+ε)^p−1)+o(1)]x^p)`. In particular, conditional on `X>x`, `Y_m/x` and `Φ(T_m)/(x m^μ)` lie within any prescribed fixed tolerance of one with exponentially high probability for sufficiently large fixed `M`. This is a conditional statement on the eigenfunction, not on a limiting empirical tree or a unique favorite growth history.

The embedded martingale can also be evaluated at the first **size** when `Y_n` exceeds a threshold, on the event that such a size is reached, by conditioning separately on each value of that stopping size. Its residual bound then uses that actual size, not a deterministic proxy size. A first hit is not the terminal event `X>x`; conversely the proof above avoids any first-hit approximation. A physical-time jump count or a first physical-time crossing also must not be identified with an unconditioned urn balance event. Independently sampled holding times at a deterministic size do not justify treating holding times stopped by a random size or selected by a physical-time event as unconditioned. The likelihood comparison and both tail bounds here use only genuine arrows and deterministic final sizes.

## 5. Deterministic bounds on the rate and exact controls

If `L` bounds the absolute one-leaf change in `Φ`, the conditional range length of a centered embedded increment is at most `2L γ_(j+1)`. Hoeffding's lemma therefore strengthens the one-sided bound by using `L² Σ_(j≥m)γ_(j+1)^2` as its variance proxy. Since that sum without `L²` is asymptotic to `m^(1−2μ)/(2μ−1)`, taking `m∼t x^p` yields

```text
I ≥ (2μ−1)/(2L²) · max_(0<t<B^(−p)) t^(2μ−1) (1−B t^a)^2.
```

The maximizing `t` satisfies `B t^a=2μ−1`. For an upper bound let `b_cap=liminf_m min Φ(S)/m`, with the minimum over size-`m` free trees of maximum degree three. This number is positive by the reverse comparison with core count. The cap survival cost at most `log 12` per step and the variational formula give `I≤log(12)/b_cap^p` by using each fixed `0<b<b_cap` and sending `b↑b_cap`.

```text
$ taskset -c 8 nice -n 10 timeout 8s python3 -B round2/open/limit-law/right/check.py
PASS: 108 exact urn compositions; scaffold avoidance and balance
PASS: exact directed-arrow likelihood and 77 cut-source generator rows
COUNTEREXAMPLE: stopping by physical time changes the clock probability from 7/8 to 8/9
```

The final fixture uses departure at size four from `T₃`: the first wait is independent of whether that departure happens, but it is truncated if departure is additionally required by a time cutoff. The asymptotic rate is proved by the inequalities above rather than estimated by the fixtures.
