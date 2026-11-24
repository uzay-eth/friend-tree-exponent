> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** supporting Astra work; the cited inputs and any recorded review are scoped in the linked claim index and historical record; no blanket independent review is asserted. See the [release index](../../../../../README.md) for exact scope and citations.

# Critical values: degree-prime signs do not specialize

**Result.** The `02d0bd8` physical enclosure raises the least possible rational denominator to **175**; the least quadratic height and least rational square-discriminant denominator do not improve. More substantially, even the exact first degree-prime congruences of Neon-2 Claim 4, together with a positive rational chronological expansion and a simple unique crossing, cannot exclude **any** subclass of algebraic points in the enclosure using only those properties. Below is a countermodel for every prescribed `a∈(0,1)` that reproduces those congruences for each scalar separately, and a separate example retaining the **entire** physical expansion at infinity. Neither is an alternative growth law for the physical tree. The arithmetic nature of `μ` remains open.

## 1. What the sharper enclosure rules out

Use the two exact endpoints in the [EDGE audit](../../../../../../../reviews/cold/astra/checks/neon-edge-sol.md), conditional on its inherited certificates, at `02d0bd87c64b53a2bda2004a4af65de9142dda13`. The outward interval is `[0.565677980686114,0.565770513279953]`. In the exact interval the first rational by denominator is `99/175`; the previous `43/76` is above the upper edge. The Farey certificate is particularly short: `56/99 < lo ≤ hi < 43/76`, the two outside fractions have determinant one, and their mediant `99/175` is inside. Thus every rational inside has denominator at least `76+99`.

With naive primitive polynomial height `max(|a|,|b|,|c|)`, the least quadratic height remains seven: the smaller root of `3s²−7s+3` is inside. The rational square-discriminant parameterization in Long's work is `s=u/v−v/u`, `u/v>1` reduced; `gcd(u²−v²,uv)=1`. The least denominator remains 57,200 with `32361/57200`, while Long's `2603/4602` is below the interval. The exact tests use monotonicity of the parameterization, and for quadratics take the extrema of `−a s²−b s` at the endpoints and any interior vertex. The square-discriminant class remains dense; these denominator cutoffs do not exclude the whole class. Long's original inert first-promotion height law distinguishes the square-discriminant class, but its coefficients are not the coefficients of the physical critical equation.

## 2. Why the large-degree prime does not evaluate at positive algebraic discount

In Neon-2 Claim 4 (`../../../../../ft2-board/neon2.md#claim-4--transcendence-over-mathbb-cs-of-usj-the-physical-degree-two-trigger-ts--fs-and-all-361-entries-k_ijs-track-2b`, historical reference), the degree-`p` star is the first place a denominator prime can enter a particular coefficient of the **large-`s` Poincaré expansion**. For the physical trigger, the new term is at `s^(-p)` and equals `(p−2)!(p−1)/(2p)` modulo `p`-integral terms. One further action after reaching the star is essential because the star itself earns no trigger reward; it sends an untriggered direct newborn into label 20. This controls a formal coefficient, not a Laurent pole at the positive value or a congruence for the value of the full continuation. The zero-gauge response instead sees the promotion reward at the star itself and its first term is at `s^(-(p−1))`, without the factor `1/2`. The stronger entrywise statements likewise concern different coefficients of different functions; they give no determinant congruence for `I−K(a)`.

There is a direct change of valuation already in the true chronological summand. At a good prime over `p` for algebraic `a` (unramified and `a`, the fixed marks, and `1/2` integral), the degree-`p` star has holding factor

```text
1/(a + a_i + p − 1/p) = p/[p(a+a_i+p)−1].
```

It has valuation **+1** in the normalization with `v(p)=1`. The singular direct emission `(p−1)/p` cancels this factor, giving an integral product and, modulo the prime, residue `1`. The upper- and lower-size rows elsewhere can bring unrelated nonunit denominators. Thus the star coefficient computed after first expanding the holding factors at infinity cannot be carried to a fixed algebraic `a` by reading the same sign or valuation. For a pole at infinity the formal `z=1/s` grading singles out the first permitted path. At a fixed ordinary discount every chronological length can contribute, and no adelic convergence or unique minimum valuation is established. Continued fractions of the real root do not repair the order of these operations.

## 3. Countermodel with the exact first-prime scalar congruences

Fix **any** real `a∈(0,1)`; it can in particular be algebraic of an arbitrarily prescribed degree (where that degree has points in the interval), a rational square-discriminant point, or the explicit quadratic above. Define two seed responses, indexed by `(k_p,b_p)`, with primes `p≥5`:

```text
trigger:     k_p=p,    b_p=(p−2)!(p−1)/(2p),      f₀(s)=1/[2(s+2)];
zero gauge:  k_p=p−1,  b_p=(p−2)!(p−1)/p,         f₀(s)=1/[2(s+2)];
S(s) = f₀(s) + Σ_p b_p/(s+p!)^k_p .
```

These two models are separate and only match the stated congruences. Each term is a positive serial holding path with rational rates/weights (equivalently its inverse Laplace density is a positive multiple of `t^(k_p−1)e^(−p!t)`). For example the summand is bounded by `2^(-p)` on `s≥0`: `b_p≤p!` and `(p!)^(k_p−1)≥2^p` for `p≥5`. Hence `0<S(a)<1`. Let `B=1+a∈(1,2)`, `u₀=1−S(a)`, and recursively set

```text
d_n = floor(B u_(n−1)) ∈ {0,1},   u_n = B u_(n−1)−d_n ∈ [0,1),
H_a(s) = S(s)+Σ_(n≥1) d_n/(s+1)^n .
```

This specifies a sequence of **rational** weights even for algebraic `a`; `a` encodes the digits, not a real-valued rate in any individual history. The identity `u₀−Σ_(n≤N)d_n B^(-n)=B^(-N)u_N` gives `H_a(a)=1` exactly, with no convergence rate assumption beyond `B>1`. On `Re s>0` the series converge locally uniformly, as do derivatives, and the nonzero positive immediate term implies `H_a'(s)<0` for real `s>0`. Hence the crossing is simple and unique. This is an explicit construction for a fixed positive algebraic trial point; possible ties in the floor are resolved by the displayed mathematical convention, not by a numerical approximation.

The large-`s` expansion is a complete Poincaré expansion. For each order `m`, include the finitely many primes with `k_p≤m+1`; the remaining summands are `O(s^(-m-1))` for real `s→∞`. Indeed one can bound each by `s^(-m-1)b_p/(p!)^(k_p-m-1)` when `k_p≥m+1`, and that constant sum converges after removing finitely many terms. The digits term is an analytic power series in `1/(s+1)` at infinity and has **integer** expansion coefficients in `1/s`. In coefficient `s^(-k_p)`, primes `q>p` are absent, and the earlier spikes have denominators prime to `p`: expanding integer shifts `q!` creates no new denominators. `f₀` is also `p`-integral. Thus the residue of the first nonintegral term at `p` is exactly the displayed physical `b_p`, not just an example with arbitrary denominator primes. The same finite-prime argument as Claim 4 proves function transcendence. Nevertheless the critical point is `a`. This combines the relevant scalar coefficient signs with positivity and a rational chronological realization; it does not preserve the entire physical generator or joint relations among all 361 entries.

### Even the complete expansion does not suffice without a uniqueness class

On the right half-plane where the genuine physical trigger `F(s)` is defined, for any fixed positive algebraic `a` in that region put

```text
G_a(s)=F(s)+(1−F(a)) exp(a−s).
```

Then `G_a(a)=1` and **every** Poincaré coefficient agrees with the physical one. To additionally preserve complete monotonicity and a single simple crossing, choose `a` strictly above the certified upper endpoint (then `1−F(a)>0` by monotonicity). The extra term is the Laplace transform of a positive atom at elapsed time one; it is flat on the large positive ray. For a point inside the certified interval this formula still matches the entire expansion and crosses the value one, but it need not preserve positivity when the correction is negative. The atom model is not a claim about a rational action-chain realization. Together the two models separate two missing hypotheses: arithmetic control on the physical continuation, and uniqueness/special-value rigidity under any proposed analytic resummation. A convergent analytic germ at infinity is not available from these Poincaré facts, and general Borel uniqueness, even if proved in an admissible class, would identify a function, not prove an arithmetic value theorem for it.

## 4. A quantitative target for genuine approximants

The existing chronological obstruction in the [parent note](../README.md) uses the positive physical omission and a surviving small holding pole. It shows that the *primitive numerator value* at `μ` grows exponentially for that approximant choice; it should not be generalized to resummed finite-state resolvents. Long's full-environment primitive height and signed noncancellation results likewise apply to his specified retained-history families. His proved nonsummability, or a nonzero numerator, is not a small-value estimate. Geometric decay in generation after an infinite occupation integral has been supplied is not an arithmetic rate in finite size.

For the actual cap-20 size-killed trigger `F_N`, use `D_N(s)=det(H(sI−Rhat_N))`, with `H=lcm(1,…,max(N,20))`, and `P_N(s)=D_N(s)(1−F_N(s))∈Z[s]` from the one-column determinant lemma. Let `a` have degree `d`, primitive minimal polynomial leading coefficient `b`, and conjugates `a_1=a,…,a_d`. Put `R_j=max(1,|a_j|)` and let `M` be the retained dimension. The row-sum argument gives coefficient ℓ¹ norm ≤ `A_N=[H(4N+4)]^M` for **each** determinant, so for `P_N` (the augmented determinant) the same bound applies. If `P_N(a)≠0`, the resultant has absolute value at least one, yielding

```text
|P_N(a)| ≥ b^(-M) A_N^(-(d−1)) (Π_(j>1) R_j)^(-M).
```

At a hypothetical `a=μ`, nonzero follows from the positive omitted path and invertibility of the killed resolvent. A sufficient small-value budget for this family would be

```text
|1−F_N(a)| < b^(-M) A_N^(-(d−1)) (Π_(j>1)R_j)^(-M) / |D_N(a)|.
```

For a uniform crude sufficient target one may replace `|D_N(a)|` by `A_N R_1^M`. This gives precision of exponential order in `M(log H+log N)` even at fixed small degree/height. The proved upper bound on the unscaled size tail at criticality is `O(N^(-.11764))`. The single omitted path to degree `N+1` followed by a reset to `J_20` gives an explicit lower bound, from the same leaf-gain product as in the parent note:

```text
1−F_N(μ) ≥ 1/[(N+1)^2 (N+2)(N+3)]     (N≥20).
```

There is a stronger obstruction for the **raw** cleared determinant, even with unrestricted lengths. The [physical semigroup comparison at `s*=0.5641`](../../../../../../astra/source/round2/finite/gauge/physical-trigger.md#4-positive-chronological-series-and-uniform-error) dominates the killed positive semigroup, so each eigenvalue `λ` of `Rhat_N` has `Re λ≤s*`. With `δ=lo−s*>0`, taken from the exact lower endpoint, `|μ−λ|≥δ`. The determinant is positive at `μ` by pairing complex conjugate eigenvalues, and `[Hδ]^M≤D_N(μ)`. Already `Hδ>2` for `N≥20` and `M≥N−1` by the reachable stars of sizes `2` through `N`. Thus

```text
P_N(μ) ≥ (Hδ)^M / [(N+1)^2 (N+2)(N+3)]  →  ∞.
```

This uses no algebraic hypothesis on `μ`. It gives no lower bound after polynomial common-factor and content cancellation; the factor `H^M` in the raw determinant makes that distinction essential. A concrete next step would be to bound the denominator of the reduced scalar rational function using only poles that remain visible from the trigger source and reward. Such a bound could extend the obstruction to primitive size-killed numerators; it would not itself exclude algebraic values. For an exclusion one still needs a genuine sequence of small nonzero integer forms or a special-value theorem. No infinite algebraic class is excluded here.

## Exact replay

```sh
taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/open/arithmetic/critical-value/check.py
```

The replay independently enumerates the bracket exclusions, verifies the two scalar congruences through prime 31, certifies a finite greedy prefix at the rational `99/175`, and checks cancellation of the degree-23 reciprocal in an actual star holding factor. The infinite countermodels and resultant estimates are analytical.
