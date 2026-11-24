> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../README.md#claim-17) for exact scope and citations.

# Matching left-tail order and the negative-moment threshold

For the Claim 6 limit from one edge, there exist `0<c<C<∞` such that for sufficiently small `x>0`

```
c x^(1/μ) ≤ P(X≤x) ≤ C x^(1/μ).
```

Consequently `E X^(−p)<∞` exactly for `0<p<1/μ`. The proof gives two-sided order, not an asymptotic equivalent, and makes no non-atomicity claim away from zero.

## 1. A degree-sensitive gradient

Use the continuous-time normalization `Φ(S)=lim exp(−μt) E_S C_t` and let `Z_S` be the limit of its continuous-time eigenmartingale from free seed `S`. The [free attachment partition](../README.md) and [conditional infection bound](../../../../../../../reviews/cold/astra/influence/check/README.md) use the same full-envelope clocks as here.

For an extra leaf at a site initially of degree `d≥2`, let `τ` be the first surplus send from the infected site to **any** of its common neighbours; this includes neighbours appended after time zero. Before `τ`, only the attachment vertex has a discrepancy among paired sites. Its private-neighbour count `k` starts at one. Before `τ` the only additions to `k` are sends from its private neighbours, at combined rate at most `k`; paired neighbours have equal degrees and send to the site at common rates. Hence `E[k_t;τ>t]≤exp(t)`. The total surplus intensity is at most `min(1,k/d)` since the lower degree never decreases. The compensation formula then gives, for the fixed exponent `0<μ<1`,

```
E exp(−μτ) ≤ ∫₀∞ exp(−μt) min(1,exp(t)/d) dt
           ≤ [1/μ+1/(1−μ)] d^(−μ).
```

This controls all retained infections *other than the seed*: conditional on the first surplus send, the total discounted mass of retained first infections thereafter is bounded by a constant even if the surplus went to an accounted vertex rather than a retained one. Explicitly the first surplus may infect at most one eligible old neighbour, and thereafter the seed spends at most rate one on new retained infections. Every retained nonseed vertex starts at its own first infection with the usual conditional `.84` offspring bound. This gives a factor `(1+1/μ)/(1−.84)` times the displayed estimate. No restart with a one-unit degree gap is made at the seed.

Before `τ`, compare the seed's new branch accounts in pairs instead of paying both full cell bounds from time zero. A paired account can first disagree only by a surplus send across its anchor edge; this is impossible before `τ`. Private accounts can first promote before `τ` via the higher-copy anchor choosing a private leaf. This intensity is also bounded by `min(1,k/d)`. After `τ`, charge first promotions of still unpromoted accounts as in the free gradient; at most one per unit time on each side. **Already promoted paired accounts before `τ`** cannot simply be charged as a single new core at `τ`: they may have grown large. Instead keep each such account matched until its own first mismatch. An account cannot disagree until a surplus send through its boundary edge. If it is already promoted, the first discrepancy is one extra leaf at its paired root, and the global finite-time one-leaf comparison still works in the surrounding pair as follows. Reroute that first infection into the retained recursion, using its exact-one private leaf and its unique incoming infector. Count it and later retained descendants using the same `.84` budget *on edges present at that infection*. Its newly attached branches are paid by the spatial bound, as usual. For an account first infected while still a leaf the identical argument starts with a singleton, whose core indicators differ at most one.

One precise accounting implements this by delaying *all* common-account cuts at the initial seed: every common neighbour of the initial seed, whether pre-existing or born later, may be retained at its first infection. Only private branches at the seed are charged as branch accounts directly. All other retained first infections follow the original freeze-at-infection protocol. A clean common branch at the seed contributes zero until it is infected from the seed, and all its surviving clean portions enter the usual recursive partition then. The rate-one post-`τ` seed allowance covers the newly eligible neighbours, with no multiplication by their number. For the seed private branches, first-promotion intensity before `τ` is at most `min(1,k/d)`, after `τ` at most one. The retained singleton contribution vanishes in the normalized limit by the already proved `.525<μ` comparison; this modified protocol still has seed discounted mass at most `1/.525` and the same nonseed contraction. Therefore

```
|Φ(S+_v)−Φ(S)| ≤ D d_v^(−μ)  for every vertex v,
```

for a uniform finite `D` (at `d=1` use the existing free gradient). This accounts for paired branches even when their first promotion precedes the first discrepancy. Their possibly large populations are never estimated from the time-`τ` branch count.

## 2. Conditional variance and a uniform negative-moment start

Write `r_v=Σ_(u∼v)1/d_u` for the physical arrival rate at a vertex and `b(S)=Σ_core L_v/d_v` for total promotions. A finite free tree with a core has no leaf-to-leaf edge. Thus `Σ_leaf r_v=b(S)≤C(S)` and `r_v≤d_v` on cores. Since `2μ>1`,

```
Σ_v r_v [Φ(S+_v)−Φ(S)]² ≤ 2 D² C(S).
```

The continuous-time predictable quadratic variation of `exp(−μt)Φ(S_t)` is the integral of this jump-rate expression against `exp(−2μt)`. The spatial cell bound gives `E_S C_t≤C_Q C(S)exp(μt)`, uniformly in all finite seeds with a core. Localization followed by martingale isometry therefore gives

```
Var(Z_S) ≤ V C(S),             E Z_S=Φ(S).
```

The identical inequality applies conditionally from a stopping time. In particular the many leaves in a star do not pay the uniform jump bound at total rate of order its size.

Only the existing reverse comparison for sufficiently large core counts is needed: take `a>0` and `K₀` with `Φ(S)≥a C(S)` when `C(S)≥K₀`. At every seed with at least `K≥K₀` cores, Chebyshev gives a small uniform probability of falling below `a K/2`. Iterating that estimate at deterministic times would be invalid because the failures refer to the same future limit. Instead iterate actual downward passages of the martingale.

## 3. Stopped downward passages

For a process starting from any seed with `C≥K`, take `a>0` with `Φ≥aC` on all such seeds. Since `C_t` never decreases, the eigenmartingale `M_t=exp(−μt)Φ(S_t)` is always at least `a K exp(−μt)`. For `y>0` below `a K/2`, let `σ_j` be the first passage below `2^(-j) aK` (set `σ_0=0`). On a finite passage at level `j≥1`, the value may overshoot by at most `D exp(−μσ_j)` from the attachment bound, and continuous drift cannot overshoot. Choosing `K≥4D/a` guarantees that the value at passage is at least `3/4` of the level: indeed before a passage jump its absolute normalized size is at most the pre-jump value times `D/(aK)`, and after the jump the level is at least the post-jump value. Equivalently the post-jump loss is at most `D/(aK)` times the post-jump value, using the same lower comparison there. Hence the next level is separated by at least one third of the martingale's current value. Conditional on a passage state the expected future squared fluctuation, divided by the squared present value, is at most `V/(a² K)` by §2 and the common time discount. Doob's maximal inequality for the centered martingale (or stopping at the next passage and using conditional isometry) bounds the conditional probability of the next passage by `η_K=16V/(a²K)`. The same bound applies to the first passage from `σ_0` since `M_0≥aK`. There is no skipping of levels for the chosen `K`.

If `Z_S≤y`, convergence forces a passage below each level strictly above `y`. The strong Markov iteration gives, uniformly at seeds with `C≥K`,

```
P(Z_S≤y) ≤ const(K) y^p
```

for any preassigned finite `p` by taking `K` large enough that `η_K<2^(−p)`. In particular choose such a `K` with a uniform bounded negative moment of order `r` for any fixed required `r` by using `p>r`. This is conditional control on an actual hitting state, not independence of its cell descendants.

## 4. Star departure is the unique rate-one bottleneck needed here

Work first in continuous time from `T₃`. Until the first promotion the tree is a star, with total promotion rate exactly one: the centre's clock always selects a leaf. Its first departure time `τ_*` is exponential of rate one. In discrete size the exact counterpart is `P(T_n is a star)=2/(n−1)`, since at size `j` the centre is selected as source with probability `1/j` and no other source promotes. The time and exit state need not be independent: a later departure typically sees a larger star.

For *every* tree with at least one core the promotion rate is at least one. When there are two or more cores, the core tree has at least two endpoints, each with at least one real leaf and hence promotion rate at least one half. For the strengthened bound, fix the large threshold `K` from §3. Whenever `2≤C<K` and `n≥K+3`, there are at least three ordinary leaves. Either some core-tree endpoint has at least two leaves, giving total rate at least `2/3+1/2=7/6`, or at least one nonendpoint has a leaf. Then its core degree is at most `K−2` and its promotion rate at least `1/K`, in addition to the two endpoints. (If there are more than two endpoints, three endpoints alone give rate at least `3/2`.) Thus the rate is uniformly at least `1+δ`, `δ=min(1/6,1/K)`, throughout that region.

Let `σ` be the first time after `τ_*` when both size is at least `K+3` and core count is at least `K`. The **additional wait** `D_*=σ−τ_*` has uniformly bounded conditional exponential moment of some order `r>1` given the full departure state/history. To see this, first wait for size `K+3` if needed. The holding times have rates equal to the current size, and the departure already has size at least four, so their finite sum has a uniform moment of any order `r<4`. Then, if there are fewer than `K` cores, at most `K−2` further promotions are needed. Each conditional waiting time has survival bounded by `exp(−(1+δ)t)` using the predictable intensity lower bound, even as the size and shape change. Conditional iteration of the moment bound gives the claim for `1<r<1+δ`. No independence between those waiting times or their exit shapes is assumed. This estimate includes all alternative bounded-core slow histories after the star, not just a prescribed deterministic history.

By the strong Markov property,

```
Z_(T₃) = exp(−μτ_*) Y,
Y = exp(−μ D_*) Z_from(σ),
E[Y^(−r/μ) | F_(τ_*)] ≤ M_r < ∞.
```

For the last inequality use §3 at `σ` with exponent `r/μ` and the conditional moment of `D_*`. Choosing `K` is noncircular: first choose it large enough in §3 for, say, negative moment order `2/μ`, then take `r<1+δ≤7/6<2`. Thus

```
P(Z_(T₃)≤x | F_(τ_*)) ≤ min(1, M_r x^(r/μ) exp(r τ_*)).
```

Integrate against the genuine rate-one exponential law of `τ_*`. Splitting at `t=log(1/x)/μ` for `0<x<1` gives a uniform upper bound `(1+M_r/(r−1)) x^(1/μ)`. The argument uses uniform conditional negative moments after departure; simply factoring an independent post-star limit would not be valid. It also shows why looking for another exceptionally cheap slow history cannot reduce the left-tail power: all histories up to the large-core threshold are covered, and the residual limit thereafter has the higher negative moment.

## 5. Return to fixed size

As in Claim 6, `Z_(T₃)=W_size^μ X` with positive `W_size` independent of `X`. This relation uses the overall core-count normalization `Φ=d₂Ψ`. Pick any finite `B>0` with `P(W_size≤B)>0`. Then

```
P(X≤x) P(W_size≤B) ≤ P(Z_(T₃)≤B^μ x) ≤ const · x^(1/μ).
```

The other direction is the exact long-star event and conditional Markov inequality from the [parent proof](../README.md): choose `n` of order `x^(−1/μ)` large enough that the conditional mean at the star is at most `x/2`, giving conditional probability at least `1/2` of `X≤x`. Tail integration now gives the stated negative-moment threshold, including divergence at equality. The independent continuous size factor smooths `Z_(T₃)` but does not prove that `X` has no positive atoms. Nor do the bounds imply a limiting prefactor for the left tail.

## Short control

```
taskset -c 8 nice -n 10 timeout 8s python3 -B round2/open/limit-law/left/check.py
```

The exact directed-arrow check covers 5,912 labelled trees, including 4,028 cases of the strict bounded-core hazard, and verifies star survival through size 79. The four-vertex path has total promotion rate one despite having two cores: it is a counterexample to dropping the size gate from §4. The high-degree gradient and the stopping-time iteration are analytical.
