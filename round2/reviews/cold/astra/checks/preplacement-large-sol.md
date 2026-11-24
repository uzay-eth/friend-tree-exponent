# Cold audit: explicit conditional large-threshold comparison

**Holds with `Y=2^90`.** For every fixed old tree `T_n` with `n≥2`, the conditional gap is strictly positive at every `0<x<1` with `(n+1)x≥2^90`. This combines the outside-block likelihood in 3a49's quantitative proof (`open/age-order/fixed-threshold/effective-uniform/README.md`) with ee0c's all-middle signed comparison (`open/age-order/global-transport/README.md`). No stationary flux, local law, or printed annealed tower is used. The 3a49 early no-change estimate is appropriate for its degree-two-parent sector only; it cannot be used for all old trees.

## Every adverse endpoint, even without stock

Use the common final size `m=n+2`. With `D_u` the **prebirth** degree and `r_u` the prebirth attachment rate, write `w_u=r_u/(D_u+1)`, `s_u=w_u/D_u`, and `A=Σ_u w_u>0`. The serial coefficient is `w_u` because only the external parent (now of degree `D_u+1`) sends to the first newborn. The parallel edge coefficient on `h_(u|v)−h_(v|u)` is `s_v−s_u`: the change in the sender's reciprocal degree is exactly `1/[D_u(D_u+1)]`, regardless of whether either endpoint is a core. Summing its absolute value over old edges is at most `Σ_u D_u s_u=A`. Each edge can cost at most `|s_v−s_u| U`, not twice that: one endpoint enters positively and the other is at most `U`. Each serial tip costs at most `w_u U`. Thus a common middle lower bound `B` gives numerator at least `A(B−2U)` before dividing by `n(n+1)`.

All adverse endpoints are leaves with permanent parent of post-placement degree at least two. Consequently `U=Q_(m,1/2)(x)` by the locked urn, even if the parent has **zero other leaves**. No stock or age cap is needed. For a nonzero zero-stock example, take a core of degree three joined to three degree-two cores; extend one of those arms by one more degree-two core, then terminate each arm with a leaf. Along the extended arm the degree-two core adjacent to the centre has `s=5/36`, above the centre's `s=1/8`. Its parallel tip is therefore adverse, yet it has zero other leaf siblings. On a star the serial tip also has zero other leaf siblings at its degree-two parent, which is the serial middle.

## Lower bound valid for every serial middle

The actual starting red mass is `a=1+1/d∈[1,3/2]` with `d=D_u+1≥2`, and `b=m−a`. Write `u=1−θ` in the likelihood calculation only. Let `G_f` be iid geometric waits (number of tag jumps before outside jump `f`) at tilted reference endpoint `θ`, `R_f=a+Σ_(h<f)G_h`, and `Δ_f=R_f−r_f`. The outside conditional arrow probability uses the actual rate divided by `b+f+Δ_f`. The exact finite-block log is `Σ_f [log(1+Δ_f/(b+f))+Σ_(j<G_f)log(1−Δ_f/(R_f+j))]`. Tag jumps do not change `Δ`; outside jumps cannot decrease it. Every tag-origin outside send has combined rate one, including sends to a core or to the old parent. Proposal uniforms dominate the count by independent Bernoullis `J_f` of means `1/(b+h)`, even conditional on the whole geometric schedule. At old-parent degree `D`, a gain of that parent has weight at most `D` and reciprocal cost `1/[D(D+1)]`; the expected increment to `Δ^old` conditional on the waits and past is at most `1/[3(b+f)]`. Hence `Δ≤1/2+J`, and `E[Δ_f | waits]≤(4/3)Σ_(h<f)1/(b+h)`. On an old star, putting the middle directly at the hub gives old-parent gain probability tending to one at the first outside jump, but its reciprocal increment is order `d^-2`: stopping on *any* old-parent change gives the wrong uniformity.

**Early blocks, quantitative.** Take `r=ceil(m/64)`, `H=512`, `z=mθ`, and suppose `z≥2^80`. The probability of any tag-origin outside arrow before `r` is at most `r/b`. On its complement `Δ_f=Δ_f^old≤1/2`, so the negative log through those blocks is at most `2Σ_(1≤f<r) Δ_f^old H_f` with `H_f=Σ_(j<G_f)1/(R_f+j)`. Given `s=Σ_(h≤f)G_h`, weak compositions are uniform. That last block includes red index `j` with probability `binom(j+f,f)/binom(s+f,f)`; using `a≥1`, the expected harmonic sum is at most `1/f`, by `binom(j+f,f)/(j+1)=binom(j+f,f−1)/f` and hockey-stick summation. Given **all** waits, the old-parent expectation is bounded by `f/(3b)`. Thus early log below `−1` costs at most `5r/(3b)<1/16`, without any independence between old-parent gains and the waits. Block zero incurs zero loss.

**Later blocks, explicit budgets.** On `[r,m)` require `|uΣG_h−θ f|≤H√z` and `J_f≤H`; on `[2^j m,2^(j+1)m)` require the bounds multiplied respectively by `2^(3j/4)` and `2^(j/8)`. Use maximal variance `Var(uG)=θ` for the first sequence; failure sums to less than `9/H²`. At each dyadic right endpoint the Bernoulli envelope mean is at most `j+4`, so its second moment is at most `(j+5)^2`; on the first interval its second moment is at most six. Summing with `2^(-1/4)<7/8` gives envelope failure less than `2^11/H²`. Integer endpoint exclusions are harmless: use the next interval's pointwise bound at its left endpoint. These bounds give `R_f≥θ f/(2u)` and `Δ_f/R_f≤1/2`. With `e_f=log(uR_f/(θ f))`, we may use `256H/√z` on the first interval and `4H·2^(-j/4)/√z` on dyadic interval `j`, including `ua≤3/2`.

The blue linearization error is at most `2^7H²/m` by summing `Δ_f²/(b+f)²`. The red linearization plus harmonic-to-log error is at most `5Σ_f Δ_f(1+Δ_f)(R_f^-1−R_(f+1)^-1)` by reciprocal telescoping; summing per block with `R_r≥z/128` and `R_(2^j m)≥2^j z/2` gives at most `2^14H²/z`. Finally replace `log(R_(f+1)/R_f)` by `log(1+1/f)+e_(f+1)−e_f`. On each block, summation by parts costs at most three times the block bound on `Δ` times the maximum `|e|`, since `Δ` is monotone. Using `Δ≤2H` on the first block and `Δ≤2H·2^(j/8)` thereafter, the total is at most `2^13H²/√z` (`2^(-1/8)<15/16`). These three errors sum to less than one at `z≥2^80`; all bounds are uniform in the terminal block. This is the 3a49 late-window calculation, with the early argument replaced, and does not require a bound on the duration of an individual wait as `θ↑1`.

The cost of the replacement sum over all `f≥1` is nonnegative and has mean less than three. Indeed `d_f=log(1+1/f)−1/(b+f)` is nonnegative and at most `b/[f(b+f)]`; harmonic summation and a decreasing right-endpoint Riemann sum give

```
E Σ Δ_f d_f ≤ (4/3)[b/(b−1)] ∫₀∞ log(1+s)/[s(1+s)] ds
            = (4/3)[b/(b−1)] π²/6 < 3.
```

The good event for the early log and the dyadic bounds has probability at least `3/4`. On it `L_∞≥exp(−2−Σ_f Δ_f d_f)`. Conditional Jensen then gives `E_θ L_∞ ≥ (3/4)e^(−6)>2^-10` (use `e<3`). Only beta-almost every `θ` is required: the nonnegative reference martingale converges under the mixture, and its limit defines a dominated physical submeasure by Fatou on cylinders. This invokes no infinite-time UI or independence of the good event and cost. Integrating over `θ>x` yields `h_M≥2^-10 Q_(m,a)(x)` whenever `mx≥2^80`.

## Locked tail separation and cutoff

The fixed-total-size beta survival increases in red mass by monotone density ratio, so `B=2^-10 (1−x)^(m−1)` works uniformly, including in the limit of arbitrarily large `d` where `a↓1`. For `b'=m−1/2`, bound the leaf integral by its maximum factor `x^-1/2` to get `U≤ Γ(m)/[Γ(1/2)Γ(m+1/2)] x^-1/2 (1−x)^(m−1/2)`. If `G~Gamma(m−1/2,1)`, Jensen gives `Γ(m)/Γ(m−1/2)=E√G≤√(m−1/2)`. Dividing by `m−1/2` and using `Γ(1/2)>1` gives `U≤2(mx)^(-1/2)(1−x)^(m−1/2)`. Thus `2U/B≤2^12 (mx)^(-1/2)√(1−x)<1` for `mx≥2^90`. Since `m=n+2`, the claimed cutoff in `(n+1)x` follows. The loss in the beta power at the high-degree middle is real; bounding it uniformly by the degree-two `a=3/2` tail would have no justification. This only controls the conditional large-threshold range. No assertion is made about every old tree at bounded birth-scale thresholds, or any annealed cutoff there.

## Replay

```text
taskset -c 13 nice -n 10 timeout 12s python3 -IB round2/checks/preplacement-large-sol.py
PASS nonzero zero-stock adverse edge and three high-degree star controls
PASS 65 conditional geometric composition bounds
PASS exact dyadic, log-error, Jensen and beta separation budgets for Y=2^90
```
