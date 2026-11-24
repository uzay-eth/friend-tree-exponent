> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports. See the [release index](../../../../../../README.md#claim-24) for exact scope and citations.

# An effective late cutoff for every threshold

For the original edge-seeded physical process, the estimates below give the following explicit (very loose) cutoff, using the same certified spectral constants as the [hub-only cutoff](../../effective/cylinder/README.md):

```math
n\ge 2^{2^{2^{123}}}\quad\Longrightarrow\quad
\Pr(Z_{n+1}>x)>\Pr(Z_{n+2}>x)\quad(0\le x<1).
```

At `x=1` the tails are equal. This is an upper bound on a starting label, not a claim about the first starting label or the remaining finite range. The wide-threshold argument needs only a degree-one stationary flux estimate at deterministic size; the small-threshold bridge uses explicit finite cylinders, without evaluating any cylinder kernel or doing exhaustive enumeration. Credit: Fable Claim 5 gives finite-label domination and the locked Pólya upper comparison; the [transition proof](../transition/README.md) identifies the same leading chain in the exact signed numerator.

## 1. Quantitative physical lower comparison

Consider the chain middle after first placing a tag on an old leaf (now its degree-two parent) and then placing a leaf on the tag. The final size is `m=n+2`. Use the [outside-count likelihood](../uniform/upper-edge/README.md), now with red `a=3/2` and blue `b=m−3/2≥m/2`. The argument there through the exact likelihood sandwich applies verbatim to these masses. At the reference endpoint `θ`, set `u=1−θ`, `z=mθ`. The scaled geometric waits have `E[uG]=θ`, `Var(uG)=θ`. At outside count `f`, again `Δ_f≤1/2+J_f`, where `J_f` is dominated by independent Bernoulli proposals of sum of means `Σ_(h<f)(b+h)^−1`; also `EΔ_f≤2Σ_(h<f)(b+h)^−1`. Extra initial tag leaves make no contribution to the deficit and do not change these bounds.

Here are explicit log-error budgets, so the lower comparison does not rely on a qualitative uniform limit. Set `ε=1/64`, `H=512`, `r=ceil(m/64)`, and use one first block of outside indices `[r,m)`, then the dyadic blocks `[2^j m,2^(j+1)m)` for `j≥0`. On the first block require the centered sum `uΣ_(h<f)G_h−θ f` to have absolute value at most `H√z` and `J_f≤H`; on dyadic block `j` require the bounds `H√z·2^(3j/4)` and `H·2^(j/8)`. Kolmogorov's maximal inequality up to each block's right endpoint gives total count failure at most `9/H²`. For the Bernoulli envelope at the first endpoint the mean is at most two; at each dyadic right endpoint it is at most `j+4`, by integral comparison. Its second moment is at most `(j+5)^2`, so the envelope failure is at most `2^11/H²` using `2^(-1/4)<7/8` and the geometric second-moment sum, including first-block failure. All these bounds are uniform in `θ`. The probability that `Δ` changes before outside count `r` is at most `3r/b≤1/8` for `m≥192`.

If `z≥2^78`, these events imply `R_f≥θ f/(2u)` and `Δ_f/R_f≤1/2` for all `f≥r`. In the first block the log discrepancy `e_f=log(uR_f/(θ f))` is bounded by `256H/√z`; in dyadic block `j` it is at most `4H·2^(−j/4)/√z` (the initial offset `ua≤3/2` is included). In the endpoint of the first block use the dyadic bound at `m`; at any other block endpoint use the next bound, which is smaller. Use [the block likelihood expansion](../uniform/upper-edge/README.md#cancellation-across-a-block): after linearizing, replace `log(R_(f+1)/R_f)` by `log(1+1/f)` by summation by parts. Monotonicity of `Δ` bounds the cost on a block by three times its upper bound times the supremum of `|e|` at the endpoints and within the block. Using the next endpoint bounds where needed gives total cost at most `2^13 H²/√z`. For the outside linearization, `Σ Δ_f²/(b+f)²≤2^7 H²/m`, by the same blocks and `Δ_f≤2H·2^(j/8)` in dyadic block `j`. For the red linearization and harmonic-sum replacement the per-block cost is at most `5` times the upper bound on `Δ(1+Δ)` divided by its first `R`. This follows from the reciprocal telescoping inequality in the upper-edge proof with `R≥1`. Summing gives at most `2^14 H²/z`, using `R_r≥z/128` and `R_(2^j m)≥2^j z/2`. The combined log error is less than one at `z≥2^78`. These estimates apply simultaneously to every finite end index and hence to `L_∞` in the lower sandwich.

Write `W=exp(−Σ_(f≥1)Δ_f d_f)` with `0≤d_f=log(1+1/f)−1/(b+f)≤b/[f(b+f)]`. The harmonic sum obeys `Σ_(h<f)(b+h)^−1≤log(1+f/(b−1))≤[b/(b−1)]log(1+f/b)`. Thus its nonnegative mean cost is at most `[2b/(b−1)]∫₀∞ log(1+s)/[s(1+s)] ds <4` for `b≥192`, using the decreasing integrand to bound the right-endpoint Riemann sum and the integral `π_math²/6`. The good event for the log error has probability at least `6/7`, by the early-window and the two dyadic error bounds. On this event `L_∞≥e^-1 W`. Conditional Jensen uses the mean cost on the good event at most `4/(6/7)`, without assuming independence. Hence `E_θ L_∞≥(6/7)e^(-1−14/3)>2^-10` (use `e<3`). Integrating over all `θ>x` gives

```math
 \Pr_T(Z>x)\ge 2^{-10} Q_{m,3/2}(x),\qquad mx\ge 2^{78},\quad 0<x<1.\qquad\text{(1)}
```

All constants are uniform over exterior trees, including atypically large degrees, and over `θ` close to one. No stationary limit goes into (1).

## 2. The exact signed comparison for `nx≥K`

Let `a_n=EC_n` and `F_n=EΣ_core L_o/D_o`. The stationary identity is `π(L/D)=μ≥113/200`. To obtain an effective deterministic-size floor for `F_n/a_n`, on the admitted balls of the [hub cylinder](../../effective/cylinder/README.md) use `f=L/D` and use zero otherwise. Its stationary loss is at most the excluded mass in that proof (`<2^-21`), its support and norm are smaller than those of the hub lower cylinder, and the same de-Poissonization error budgets apply. In particular `F_n≥a_n/2` for `n≥2^(2^396)` (two errors each at most `2^-7`).

Let `Q_{m,a}(x)` be the survival of `Beta(a,m−a)`. It increases in `a` by the density likelihood ratio. For `x>0`, comparing the densities for `a=4/3` and `a=3/2` on `t≥x` gives

```math
 Q_{m,4/3}(x)\le 4(mx)^{-1/6} Q_{m,3/2}(x).\qquad\text{(2)}
```

For `s=m−3/2`, Jensen applied to the `5/6` power of a gamma variable of shape `s+1/6` gives `Γ(s+1)/Γ(s+1/6)≤(s+1/6)^(5/6)`. Divide by `s` to get `Γ(s)/Γ(s+1/6)≤2m^(-1/6)` for `m≥4` (use `s≥5m/8` and `s+1/6≤m`). For the other coefficient, `Γ(3/2)≤1` by Jensen and `Γ(4/3)≥1/2`: on `[0,1]` integrate `t^(1/3)(1−t)` and on `[1,∞)` integrate `e^(-t)`. The density's ratio in `t` is `((1−t)/t)^(1/6)≤x^(−1/6)`.

The exact two-birth signed numerator from `T_n` uses the endpoint tails at the *same* threshold on the common size `m=n+2`. Its leading middle terms have coefficients summing in expectation to `F_n/2` (old leaf, then attach to its newborn). Every other endpoint is a tip with red mass `≤1/2` or a middle with an old core parent and red mass `≤4/3`. Their upper physical comparisons are the corresponding locked beta laws. The total absolute coefficient budget is at most `6a_n` counting both endpoints, before the denominator `n(n+1)`. This yields

```math
 n(n+1)(p_{n+1}(x)-p_{n+2}(x))
 \ge a_n Q_{m,3/2}(x)[2^{-12}-24(mx)^{-1/6}].\qquad\text{(3)}
```

In particular `K=2^102` suffices whenever `nx≥K`, since `24<2^5`. No adjacent one-label approximation is used and there is no upper restriction on `x<1`.

## 3. Bounded `y=nx`: quantitative stationary margin

Let `G_y` be the stationary signed score with test `W>y`, and `F_y` the single-birth score, from the [birth-scale ordering](../../scaled/README.md). There `πG_y≥(1−μ)πF_y`. For an explicit uniform positive lower bound choose `d*=2^118≥1024K`. On an old leaf placement, demand that the post-placement degree-two parent send to the tag before its own degree grows and within one time unit. One can require its first ring within one and suppress gains at the parent before that time by thinning at total rate two, independently of its own rate-one clock; then it chooses the tag with probability `1/2`. Integrating the first ring density against the silence probability gives `(1−e^-3)/6>1/32`. On promotion the tag has one leaf. Silence the tag for the next unit and run the Yule family from that leaf for one unit. Its chance of reaching at least `d*` is `(1−e^-1)^(d*−1)≥2^-d*`; tag silence costs at least `1/4`. At that time the remaining continuous mass is at least `Gamma(d*,1)−Exp(1)` in the coupling of the [hub mass bound](../../../../../../../astra/source/round2/fringe/hub-spectrum.md#2-a-vertexs-limiting-mass-as-seen-at-time-t). Chebyshev at `d*/2` and the exponential tail at `d*/4` give a chance at least `1/2` of remaining mass at least `d*/4`, without assuming their independence. Since elapsed time is at most two and `e²<9`, the original mass then exceeds `K`. Old-leaf placements have stationary flux `μ`. Using `μ≥113/200` and `1−μ≥1−1133/2000`, we may take the common stationary margin

```math
 \pi G_y\ge \beta_K:=2^{-A},\quad 0\le y\le K,\qquad A=2^{119}.\qquad\text{(4)}
```

The following cylinder bounds aim at errors `2^-A` absolutely; their exponents are not intended to be sharp.

## 4. Low mass and a positive-threshold modulus

All bounds in this section are uniform in the exterior and for post-insertion marked degrees one and two. Put `B=A+100`, `d=2^B`, `M=2^(4B)`, `h=16B` and

```
T = d · (2d+32B) · 128B,    y0 = 2^(-4(T+h)),    t = 64(T+h+B),
δ = 2^(-4t),                C = 2^(8B),            U = 2^(32B).
```

The [finite-late-hit argument](../../effective/cylinder/README.md#1-a-uniform-modulus-for-a-finite-next-gain) generalizes directly: use its first two stages of length `d`, then a birth–death family of rates `1,1/2` for `32B` units to reach stock `M`. Its mean is `e^(16B)` and Paley–Zygmund at half its mean still gives at least `1/16`. Each block has success floor `1/32`; subsequent send from stock `M` costs at most `2/(M+1)`. Use `128B` blocks per gap and fewer than `d²` (gap, source) cases. A finite degree-`d` hit later than `T` has probability at most `3·2^(-2B)`. At a hit by `T`, the fresh leaf has a family dominating a linear birth–death process of birth rate one and death rate `1/d` per member, since the true degree never decreases. Its probability of being below `d` at age `h` is at most `1/d+2d e^(-(1−1/d)h)`: the linear birth–death generating function gives an atom of size at most `1/d` at zero and, conditional on nonzero, a geometric law with parameter `(1−1/d)e^(-(1−1/d)h)/(1−d^(-1)e^(-(1−1/d)h))`. From stock at least `d` the gamma-minus-exponential comparison gives remaining mass at least `d/4` except with probability `4/d+e^(-d/4)`. Therefore `Pr(0<W≤2y0)≤2^(-B+4)`. This also bounds the error incurred by using the hub event in place of `W>y` for `y≤2y0`.

Here is a separate interval bound for `y≥y0`. At observation time `t`, with conditional probability lost at most `(t+3)/C`, the tag's core degree is at most `C`; it is at most its initial degree plus one plus the number of its rate-one rings. If its leaf count is below `U`, the conditional probability of `W≥y0/2` is at most `2(U+C)e^-t/y0`, by the gamma upper mean. Otherwise, conditional on the time-`t` history, the gamma upper and gamma-minus-exponential lower comparisons bound the conditional mass of any interval of original length at most `16δ` in `[y0/2,∞)` by

```
 (16δ e^t + C + B + 1) / sqrt(U−1) + e^(-B).
```

Indeed the gamma density with integer shape `L≥U` is bounded by `1/sqrt(L−1)` (evaluate at the mode and use Stirling's factorial lower bound), the gamma with shape `D` can be obtained by independently adding `Gamma(D−L,1)` for this comparison of marginal CDFs, and the exponential subtraction is paid at `B`. The bound also covers closed intervals. Since `T+h+B<2^(3B)`, the sum of all three displayed costs is below `2^(-B+5)`. Initial state bounds on degrees of any other vertex are unnecessary.

## 5. An explicit finite cylinder for all bounded y

Use a new observation time `S=64(t+B)`. For either inserted endpoint, estimate `W` by `W_S=e^-S D(S)`. The [birth-scale martingale tail](../../../../../../../../reviews/cold/astra/checks/adjacent-age-sol.md#positive-tails-a-birth-scale-extension) gives `E|W−W_S|≤√2 e^(-S/2)+(S+4)e^-S`, so `Pr(|W−W_S|>δ)<2^-B`. For a finite size `m=n+2`, write `mZ=(m/Y)W`, `Y~Gamma(m,1)`, without assuming independence. On `|Y/m−1|≤δ/(K+1)`, the threshold `yY/n` for testing `W` differs from `y` by at most `2δ` for `0≤y≤K` and `n≥4(K+1)/δ`; the extra factor `m/n` is included. The failure probability is at most `(K+1)^2/(m δ²)<2^-B` at the stated cutoff. Thus for `y≥y0`, compare with finite-time reads at `y±4δ`, and use the interval bound of §4 to pay their joint stationary gap. For `y≤y0`, on that same gamma event the threshold is at most `2y0`. Use the degree-`d` hit by time `T`: its upper error is the late-hit bound and its lower error is the fresh-family estimate at threshold `2y0` (plus the size-factor failure bound). These are pointwise sandwiches on finite trees, not just on stationary samples.

Use radius `R=2^(5B)` in the old core, read radius `k=R+1`, core-degree cap `J_c=2^(6B)`, and total-degree cap `D=2^E` with `E=2^(7B)`. As in the [weighted exclusion bound](../../effective/cylinder/README.md#3-weighted-loss-of-every-excluded-cylinder), the excluded stationary mass is at most

```
 (1+R) J_c^(R−1) [32·2^(-E/8) + 512·2^(-floor(J_c/8))] < 2^-B.
```

On admitted balls the queried post-insertion cells are based within distance five. Their initial *core* degrees have sum at most `2^12 (J_c+2)^7`, independently of the total-degree cap. The ordered-path error by time `S` is thus `<2^-B`: apply the bound with length `R−5`, use `(R−5)!≥((R−5)/e)^(R−5)` and `R−5≥8S`. Retain all boundary degree slots as in the hub cylinder. Kill at size `B_{
m kill}=2^{E+10BR+4S+4B}`; Yule domination gives error `<2^-B` from killing even after the insertions. Rational Euler uniformization with `2^{4B} (B_{
m kill} S)^2` steps (round the step count up to an integer) gives error `<2^-B` as in the hub cylinder. For the positive readout one can use dyadic threshold cells of width at most `δ` and rational bounds on `e^-S`, so that the readout threshold changes by at most an additional `2δ` (approximate it even at degree `B_kill`). The lower and upper reads still differ only within the `16δ` interval already budgeted. For the hub-style hit use the hit flag. The same support, time and arithmetic budgets work at every threshold, with no need to tabulate the threshold grid.

Insert the endpoint bounds in the exact signed placement formula, respecting the sign of each coefficient, and clip below at `−3`. Use `−3` on excluded balls. The resulting cylinder `L_y` satisfies `L_y≤G_y` in the stationary infinite law and bounds the deterministic-size finite score from below, with the finite size gamma error accounted for in the endpoint bounds. Its stationary loss is less than `2^(-B+12)` by the absolute endpoint-weight budget of six, the excluded-ball cost of six, and §§4–5. The same cylinder may be used for all `n` beyond the claimed cutoff: pay the uniform gamma error `2^-B` up front. Consequently `πL_y≥2^-A−2^(-B+12)` uniformly for `0≤y≤K`.

For the implemented rational readout, center the two tests at `y±5δ` before rounding, rather than using the unrounded `y±4δ` tests above. Round the threshold within `δ` and approximate `e^-S` with error at most `δ/B_kill`. The resulting `2δ` error leaves `3δ` of inward slack for the `2δ` size-factor shift and `δ` martingale error. On `|W−W_S|≤δ`, the stationary readouts can disagree only for `W∈[y−8δ,y+8δ]`, exactly the `16δ` interval budget in §4.

Apply the certified [signed de-Poissonization](../../effective/README.md#3-return-to-the-deterministic-original-labels) to `L_y`, with `M=6`, radius `k`, and jump bound `96(k+1)(D+1)^(k+1)+6`. The centered cylinder has constant part `−3−πL_y`, of magnitude at most six. The shell at depth `R+1` is read only for incidence to depth-`R` vertices, not for its marks. Changes of admitted summands can therefore be charged to the recipient or the promotion site within the radius-`R` ball on the admitted side of the jump. The reversed paths used in the jump bound stay within capped balls. The core lower bound and the two errors, at most `2^(-A−4)` each, follow at `N=2^s`, `s=2^(16A)`. Since `A=2^119`, the exact height is `log₂ log₂ N = 16A = 2^123`, hence `N=2^(2^(2^123))`. The sufficient tests are

```
s ≥ 515·2^348,
s ≥ 20·2^348 (982+2k+A+8),
s ≥ 25 (71+ceil(log2 J)+A+8),
ceil(log2 J) ≤ 16+5B+(k+1)(E+1).
```

All error envelopes decrease thereafter. The normalized expected signed score is strictly positive on `0≤nx≤K`, and (3) covers `nx≥K`. This supplies a numerical tower rather than only a procedure that terminates after finding acceptable stationary integrals. It assumes the same certified spectral constants as the earlier hub-only numerical result.

## Replay

```
python3 -IB round2/open/age-order/fixed-threshold/effective-uniform/check.py
PASS rational likelihood, stationary flux, and beta separation; K=2^102
PASS symbolic height: N=2^(2^(2^123)); smaller tower rejected by the spectral floor
PASS symbolic cylinder and de-Poissonization budgets; A=2^119, log2(log2(N))=2^123
```

The checker reads the stated cutoff in this README, strips two powers of two symbolically, and checks the remaining exponent against `16A`. As a regression control, the smaller tower `2^(2^123)` fails both that identity and the necessary floor for these budget tests, `log₂ N≥515·2^348`. The checker does not expand the cutoff, enumerate trees, or evaluate stationary integrals.
