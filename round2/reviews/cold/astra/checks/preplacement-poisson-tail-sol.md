# Cold audit: Poisson–beta threshold minorant

**Holds** at proof (`open/age-order/hub-ranking/poisson-tail/README.md`) `1a436ac3a8703327dea6f9949dad334905b3c7af` and checker `0452a2890979bc03f3f59fb80a2fcd0d83079622`. The comparison holds conditionally on every old tree of size `n≥2` when `0<x<1` and `(n+1)x/(1−x)≥6`, with the stated strict surplus at the cutoff. The minorant separately holds for every `0≤x<1`, without a threshold or degree cutoff.

## Harmonic boundary and the strict tail

Use the binomial CDF with index `r−1`: it is identically zero for `r≤0`, and the saturated CDF at `r=m` is one. The urn equation therefore extends to negative indices (all terms zero) and holds at zero and saturation. In the Poisson convolution, rewrite the discrepancy using `r=j−k` and the forward difference as in the proof. **Keep `k=j` on the next-size gain term**, even though the old term at that index is zero. At `k=j`, the surviving summand is `jπ_j B_(m+1,1)`, which cancels the last nonzero negative term from `k=j−1`. Only `k≥j+1` is identically zero. The identity also works at `j=m−1`; the hold coefficient is then zero, not negative. The independent checker builds the binomial PMFs by Bernoulli convolution and checks the edge term explicitly, including `x=0`.

In the physical chain the gain and loss channels are disjoint. Relative to the harmonic test, extra gain probability and the reduction in loss probability separately produce nonnegative drift since the next-size stock score increases with stock. At stock zero the lower test is zero and all next values are nonnegative. The test is bounded between zero and one. For `x<u<1` and `Z<u`, the actual stock fraction has limsup at most `Z`, so eventually lies below `u` by a fixed positive margin on each such sample path. The binomial CDF bound then forces the test at `u` to tend to zero there. Bounded reverse Fatou bounds the initial test at `u` by `P(Z≥u)`. As `u` decreases to `x` from above, these events increase to `{Z>x}`, and the finite beta sum is continuous from the right also at `x=0`. No assumption about atoms of the physical `Z` is used. For the auxiliary chain the stock fraction is a bounded supermartingale and the harmonic is a bounded martingale. At continuity points of the limiting fraction, its harmonic converges to the tail indicator: below use the binomial bound; above retain each fixed Poisson prefix and then exhaust its mass. This identifies the law, including the missing Poisson mass at zero.

## Debit and all-size margin

The locked beta comparison can be seen directly in size time: at a fresh leaf with permanent parent degree at least `k`, the degree-gain rate at later degree `D` is at most `D−1+1/k`. Couple integer increments under a common uniform with an affine urn starting from `1/k` at total size `m`. The limiting fraction is dominated pathwise, hence the beta tail bounds even the strict event. A serial tip has parent degree two. For a parallel pair only the larger-`s` endpoint is adverse, from the exact coefficient `(s_v−s_u)` for the test of leaf `u` minus leaf `v`. At a core-side endpoint the old recipient's degree rises to at least three. An old leaf recipient needs the degree-two cap, and around a degree-`d` core with `L≥1` old leaves their aggregate expensive coefficient is at most `A_l(d−1)/(d+1)`, with `A_l=L/(2d)`. Three times this bound is at most `A_l+w_core`, since `w_core≥2d A_l/(d+1)`. Summing gives `X≤A/3`; the edge seed has zero parallel coefficient. The total edge coefficient is at most `ΣD_us_u=A`, and `A≥1` follows from `Σr_u=n` and `D_u+1≤n`. These are deterministic, not annealed statements. With a single tracked leaf on each serial middle the normalized credit is at least `1/e`; the resulting debit allocation is `4/3` at shape `1/2` and `2/3` at shape `1/3`.

For `r∈{1/2,1/3}`, substituting `z=x+qu` in the beta integral factors out `q^(B−b)x^(−s)=(q/x)^s`, where `B=b+1−r=b+s`. The remaining normalization is `Γ(B+r)/(Γ(r)Γ(B+1))`, bounded by `B^(−s)/Γ(r)` by applying Jensen to the `r`-moment of `Gamma(B,1)`. The Taylor remainder after the quadratic has negative sign for nonnegative arguments, since the third derivative is negative. The exact beta moments then yield the displayed linear correction. Uniformity really starts at `b=3`: `B≥7/2` gives `B/(B+1)≥7/9`, and `v=q/(Bx)≤1/6` follows from the cutoff since `B≥b`. The bracket is positive, and `v^s(1−d_sv)` increases through `1/6` by its derivative. The lower Taylor polynomial for the gamma integral and rational arctangent bounds give the announced constants. The independently recomputed surplus is strictly above `13901663/7005398400`, itself above `19/10000`. At `x=1` the strict tails are all zero and the claim uses `x<1`.

For the unconditional splice the cutoff threshold at older label `b` is `6/(b+6)`, decreasing in `b`. Paired early receipts through 19 then give `[3/13,1)`; adding the single `20→21` receipt gives `[2/9,1)`. This uses the inherited early receipt scope, not a treewise early assertion. Founders tie in law.

## Replay

```text
nice -n 10 timeout 60s python3 -IBu round2/open/age-order/hub-ranking/poisson-tail/check.py
PASS 2436 exact harmonic rows, 9022 tracked physical states, 104964 directed event checks, 24 saturated-star equalities
PASS normalized debit <= 8322209/22744800; normalized surplus > 19/10000 whenever (n+1)x/(1-x) >= 6

python3 -IB round2/checks/preplacement-poisson-tail-sol.py
PASS 1380 binomial PMF harmonic rows, including x=0 and the k=j edge
PASS 873 exact paired charges, 0 saturating the one-third bound
PASS exact cutoff surplus > 13901663/7005398400 (>19/10000), including b=3
```
