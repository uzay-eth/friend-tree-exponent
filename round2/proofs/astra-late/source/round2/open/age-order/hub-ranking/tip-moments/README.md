> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# Physical upper bounds for the adverse tip moment

Starting in any finite completion with a marked leaf and parent of degree `d≥2`, let `W=lim exp(−t) D_tag(t)` in rate-one vertex time. Then the following upper bounds are uniform in the entire exterior, including zero other leaf siblings:

| exponent | `d=2` | `d≥3` |
|---|---:|---:|
| `E W / Γ(2)` | `< .387` | `< .265` |
| `E W² / Γ(3)` | `< .302` | `< .185` |

An explicit convergent recursion below covers **every real `p≥1`**. A simple closed-form bound for the entire range is `(1+1/d)_p/[Γ(p+1)(1+pd)]`. No finite-label certificates or hubhood estimates are needed.

## Physical cap after the first promotion

Define `A_p(a)=Γ(a+p)/(Γ(a)Γ(p+1))`. At first promotion the tag has degree two, one fresh leaf, and its permanent parent has degree `k` unchanged by the promotion. For subsequent degree `D` the tag gain rate is at most `D−1+1/k`; its initial effective seed is therefore **`a=1+1/k`**, not `2−1/k`. This comparison is a pure birth chain on `a+integer` with gain rate equal to its current mass, coupled by common gains plus nonnegative residual rates. Its normalized limit has Gamma shape `a` and rate one: the finite-time number of increments has negative-binomial law with probability generating function `[exp(−t)/(1−(1−exp(−t))s)]^a`, whose rescaled transform tends to the Gamma transform. Using the original parent slot forever is valid even when most old leaves later promote. The comparison requires no independence between the marked degree and the population.

For a further gain in the bound, restart at the first tag ring `T` or the affine cap's first gain `V`, whichever occurs first. Use for the cap an independent rate-`a` candidate before this race ends, accepting any actual tag gain by thinning it at the actual predictable tag gain rate divided by `a`. Candidate times are independent of the tag's own rate-one clock before this race: they replace only incoming arrows; allocate actual incoming senders by their rates on a hit, and use residual arrows for all other physical targets. This directed-channel split preserves all physical rates. At `V<T`, the actual tag degree is at most three, so the effective mass is at most `a+1`. At `T<V`, the degree is still two and its only two choices are:
> * **Fresh leaf** (probability `1/2`): its promotion makes that neighbour a permanent degree-at-least-two slot, reducing the effective seed to at most `a−1/2`.
> * **Old parent** (probability `1/2`): its degree rises to at least `k+1`, reducing the effective seed to at most `1+1/(k+1)`.

Other growth of the parent can only help the upper estimate. The discounted first-race identity is `A_p(a)=[a A_p(a+1)+A_p(a)]/(p+a+1)`. Replacing the second term by the two improved seeds yields

```
B_p(k) = A_p(a) − [2 A_p(a)−A_p(a−1/2)−A_p(1+1/(k+1))]/[2(p+a+1)].
```

Both caps apply conditional on the entire physical promotion history. Keeping only the cap race up to its first event avoids an invalid assertion that all later separated outgoing clocks are independent under an endpoint tilt. `B_p(k)` is nonincreasing in `k`: write it as `A_p(a)` times `1−1/(a+p+1)` plus `[(A_p(a−1/2)+A_p(1+1/(k+1)))/2]/(a+p+1)`. In the form `M+(1−1/(a+p+1))(A_p(a)−M)` it is increasing in each of `A_p(a)`, `M`, and the mixing weight on the larger value when these are varied in order. More explicitly, the weight on `A_p(a)` increases with `a` while both values increase as `k` decreases, and `A_p(a)≥M`. The monotonicity of `A_p` follows from Gamma additivity or the strictly increasing digamma function.

## Physical discount before first promotion

Until promotion, the tag's own rate-one clock always appends a leaf at the parent, increasing its degree by one. The true promotion intensity is `1/k` at parent degree `k`; every other parent gain also increases `k` and cannot be suppressed physically when simulating the chain. They may be ignored in an **upper** comparison via a decreasing test. Let `τ` be first promotion. On `τ=∞`, `W=0`; otherwise apply the Markov property and the evolution from the actual promotion size/state to bound the moment by

```
E[W^p]/Γ(p+1) ≤ E[exp(−pτ) B_p(D_parent(τ)); τ<∞].
```

This is bounded by the value `V_p(d)` of the one-dimensional race with `k→k+1` at rate one and absorption at rate `1/k` paying `B_p(k)`. Its recursion and a rigorous finite truncation at integer `K>d` are

```
V_p(k) = [B_p(k)/k + V_p(k+1)]/(p+1+1/k),
V_p(K) ≤ B_p(K)/(1+pK).
```

To check the comparison, realize the one-dimensional race by an independent Poisson degree-increase clock and a unit exponential cumulative-hazard threshold. Started at `k`, its absorption discount is no greater than that for the constant hazard `1/k`, and its payoff is at most `B_p(k)`. This proves the terminal bound for the race, not a coupling of promotion-timed payoffs between two real trees. The recursion then shows `V_p(k)≥V_p(k+1)`: use `V_p(k+1)≤ B_p(k+1)/(1+p(k+1))≤ B_p(k)/(1+pk)`. Hence using the discounted test `V_p` in the **physical** chain, each additional parent gain contributes a nonpositive term. Stopped Dynkin with bounded tests (mark absorption with its payoff) gives the desired upper, removing finite horizons by positivity. Equivalently a finite backwards recursion with terminal upper has the same monotonicity, and at a physical hit of `K` the remaining absorption time is stochastically no shorter than an exponential of rate `1/K`. Boundary uncertainty in the race propagates down by factors at most `1/(p+1)≤1/2`. Thus the recursion retains the real Poisson return clock without needing a false pathwise comparison of the actual and surrogate rewards at differently timed promotions.

Replacing `B_p(k)` by `A_p(1+1/k)` gives the easier upper `V^0_p`, also nonincreasing in parent degree. It is decreasing in the real exponent `p≥1`: each term in its convergent expansion is a product of `1/(p+1+1/j)` ending with `A_p(1+1/k)/[1+(p+1)k]`. For `a=1+1/k∈[1,3/2]`, log differentiation of `A_p(a)` gives `ψ(p+a)−ψ(p+1)≤(a−1)ψ'(p+1)`. The sum and integral split gives `ψ'(p+1)≤1/(p+1)^2+1/(p+1)`; half of this is at most `1/(p+3/2)` for all `p≥1` by multiplying out. This is at most the logarithmic loss from the final denominator. Thus the unrefined uniform ceilings on `[1,∞)` are `.414,.284` for the two parent classes, and `.325,.199` on `[2,∞)`. The sharper pointwise recursion with `B_p` is valid without using monotonicity in `p`.

If `Z` denotes the limiting share in the physical jump chain from common size `m`, then `W=S Z` with the whole-tree Yule amplitude `S∼Gamma(m,1)` independent of the embedded tree path. This independence holds under the original physical law from the fixed completion, not conditional on a clock-defined promotion time. Therefore `E Z^p = E W^p/(m)_p` at a fixed completion, and the bound compares endpoint moments without any additional moment factor depending on the old target.

## Replay

`nice -n 10 timeout 12s python3 -IB round2/open/age-order/hub-ranking/tip-moments/check.py` passed all eight rational moment ceilings and 64 explicit physical arrow fixtures. The cutoff `K=64` is used for the numerical ceilings only; the analytical recurrence works for any finite cutoff above the starting degree.
