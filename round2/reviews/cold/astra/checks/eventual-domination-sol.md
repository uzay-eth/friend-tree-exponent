# Cold review: Claim 24, eventual domination at all thresholds

**HOLDS as a synthesis of the frozen inputs.** Reviewed Claim 24 (`open/age-order/eventual-domination/README.md`), blob `a9957a051ae53b56aaf035cf40d831acce6f3f26`. I did not author or re-prove the three regime estimates. This review checks their stated scopes, matching constants, and all deductions in the synthesis. The common cutoff is existential. No finitely many unresolved pairs are settled here.

## Inputs, scales and one common cutoff

| Pin | Input actually used |
|---|---|
| `25a2a2f` | `scaled/README.md` §2, blob `13c7c9b9…`: uniform compact limits **including 0** for labels `n+1,n+2` at one threshold `z/n`. Both profiles continuous, `B≥(1−μ)A>0`, and `B(0)=(1−μ)A(0)`. Also supplies the profile derivative identity only on `(0,∞)`. |
| `389eeff` | `fixed-threshold/uniform/README.md` §4, blob `810719e6…`: compact interior only, with `κ=E[X] π((L_o/D_o)g)/√π_math>0`. |
| `4492e16` | `fixed-threshold/uniform/upper-edge/README.md` (A) and §4, blob `d862bf91…`: for each fixed `a>0`, the **relative** tail and ratio errors are uniform on `[a,1)`, with no lower bound on `1−x`. |
| `4492e16` | `fixed-threshold/transition/README.md` (1), (5)–(7), blob `be4e5789…`: the same one-label tail and a **signed** gap equivalent along every `x→0, vx→∞`, no minimum growth rate for `vx`. |
| `4233b1a` | `spectrum/README.md` §§1,3 and `moments/README.md`: the positive stationary birth intensity, its density on the open positive axis, and the fixed-power coefficient. The `82864bc` fixed-label endpoint theorem is not used to extend either large-label range. |

For `n=v−1`, put `z=ny/v`, not `z=y` at the same finite label. Exactly,

```text
v^(1−μ) p_v(y/v) = (v/n)^(1−μ) [n^(1−μ) p_(n+1)(z/n)],
v^(2−μ) d_v(y/v) = (v/n)^(2−μ) [n^(2−μ) (p_(n+1)−p_(n+2))(z/n)].
```

Here the second threshold is also `z/n`; using `y/(v+1)` for the later label is a different difference. Since `z=y−y/v`, these equations transfer uniformly on each compact down to zero. The frozen signed denominator is `n(n+1)`, not a same-size square. In the leading transition term the single birth has size `v`, red `1/2`, blue `v−1/2`; the chain middle has size `v+1`, red `3/2`, **the same blue mass**. Its reference density divided by the single-birth density is `2vθ`. The division by `2(n+1)=2v` in the physical ratio leaves the conditional beta mean `E[Θ | Θ>x] = x+O(1/v)`, which is asymptotic to `x` only because `vx→∞`. The two annealed scores both tend, after division by `v^μ`, to `E[X] π((L_o/D_o)g)`. The common leading flux has prebirth factor `1/(v−1)` and postbirth beta survival of order `v^(−1/2)x^(−1/2)(1−x)^(v−1/2)`. The ratio `v/(v−1)` can be absorbed into the scalar prefactor. **The power `v−1/2` cannot be replaced by `n−1/2`**: that multiplies the equivalent by `(1−x)^(−1)`, already wrong at a fixed interior threshold and unbounded near one.

If no common cutoff existed, choose strictly increasing failing labels and one failing threshold per label. From those thresholds extract a subsequence with `vx` bounded or tending to infinity. The first case contradicts the positive minimum of `B` on a compact. In the second case, extract a further subsequence on which `x` converges in `[0,1]`. Limit zero uses the signed input; a positive limit (including one) eventually lies in `[a,1)` for some fixed `a>0`. Writing the total upper-edge ratio as `(1−x)(1+e_v(x))`, the gap divided by `x p_v` is `1−(1−x)e_v(x)/x`, uniformly close to one for `x≥a`. This same countersequence reasoning proves the gap equivalent along **every** `vx→∞` sequence; there is no analogous deduction by subtracting one-label equivalents when `x→0` since `e_v/x` could diverge. The compact proof covers zero, where the gaps compare probabilities of positive mass. State space is `[0,1]`; at one both strict tails are zero. Eventual strict ordering is pointwise at each threshold with a common label cutoff, not a uniform positive lower bound on the unscaled gap. Early certification through `(19,20)` plus a single-program `(20,21)` agrees with the cited Fable audit.

## Matching, density, and uniform equivalents

Given arbitrary `y_j→∞`, choose increasing `v_j` so the two compact approximations at `y_j` each have relative error at most `1/j`; this is possible separately for each finite `y_j` since both profiles are strictly positive. Impose also `v_j≥j(1+y_j)^2`. Then `x_j=y_j/v_j→0`, `v_j x_j=y_j→∞`, and the log correction `(v_j−1/2)log(1−y_j/v_j)+y_j` is `O((y_j+y_j²)/v_j)=o(1)`. The scaled one-label equivalent contributes `κ y_j^(−1/2)e^(−y_j)`; multiplying the gap equivalent by `v_j^(2−μ)` adds exactly `v_j x_j=y_j`, giving `κ y_j^(1/2)e^(−y_j)`. Thus these are profile equivalents along every divergent `y_j`. This diagonal does not assert uniform compact errors on an independently prescribed growing window. Do not differentiate an equivalent: the frozen identity gives `−A'=[B−(1−μ)A]/y`, and `(1−μ)A/B∼(1−μ)/y`, so `−A'∼κ y^(−1/2)e^(−y)`.

On bounded `y=vx`, `log E_v=(v−1/2) log(1−y/v)+y=O_K(1/v)` and `E_v(0)=1`. The positive compact minima convert the absolute compact limits into relative limits. Along every `y→∞`, the large-profile formulas cancel `e^y` in `E_v`; the exact factor `(1−x)^(v−1/2)` survives. The two denominators in (39)–(40) then match the one-label and signed equivalents, respectively. Any hypothetical failure of uniform relative convergence has a subsequence in one of these cases. This also checks probabilities smaller than any exponential near one. Dividing gives a uniform relative **gap fraction** `d_v/p_v ∼ B(vx)/(v A(vx))`, including at zero, where it is `∼(1−μ)/v`. It is positive, and the original tail ratio is strictly below one. It is strictly above zero since the equivalent is positive at all sufficiently late labels for every `x<1`. If desired, a uniform relative equivalent for the **total tail ratio** itself is `p_(v+1)/p_v ∼ 1−x`: quotient the two one-label composites to obtain the extra factor `((v+1)/v)^(μ−1) e^x A(vx+x)/A(vx)`, which tends to one on bounded `vx` by continuity (`x→0`) and on divergent `vx` by the profile equivalent and `0≤x≤1`. This ratio alone does not resolve the degenerate sign at small `x`, and near one it should not be recovered by subtracting a relative gap estimate from one.

## Powers and stationary intensity

Layer cake uses strict tails with the same power at both labels: `E Z^α=α∫_0^1 x^(α−1)P(Z>x)dx` for every `α>0`. Thus the cutoff gives strict moment order for all positive powers simultaneously, even arbitrarily close to zero; it does not include `α=0`. More generally the moment of the later label is positive, so the total moment ratio is strictly between zero and one. The fixed-power `b_α` is the **unnormalized** coefficient, not a conditional moment given hubhood. For large `α` split its integral at fixed `Y>0`. The bounded-interval part before multiplication by `α` is `O(Y^α/α)`, negligible relative to `Γ(α−1/2)` by Stirling; the missing bounded interval in the gamma integral is negligible by the same argument (with exponent `α−3/2`). Tail bounds at arbitrary relative accuracy then give `b_α∼ακΓ(α−1/2)`.

For unbounded powers, the uniform relative error of each composite can be pulled through layer cake because its comparison integrand is positive and the error is independent of `α`. Write `a=α−1/2` for the single and `a=α+1/2` for the gap, and `b=v+1/2` for both. The weights after removing the large-`y` profile are precisely `x^(a−1)(1−x)^(b−1)`. For `T~Beta(a,b)` and `a>1`,

```text
E[(vT)^(-1)] = (a+b−1)/(v(a−1)) = 1/v + b/[v(a−1)] → 0.
```

Thus `Pr(vT≤Y)≤Y E[(vT)^(-1)]` for every fixed `Y`, and `E[(vT)^(-1/2)]≤sqrt(E[(vT)^(-1)])→0`. The single-profile ratio is bounded near zero and tends to one at infinity. The gap-profile ratio can grow like `y^(−1/2)` near zero, but is bounded by `C(1+y^(−1/2))`, so the first two inequalities control its full small-`y` contribution. On `vT>Y` take `Y` large enough for any prescribed relative error. This proves the two beta equivalents with no relation between `α` and `v` and no uniform-in-power assumption from the fixed-power note.

With `a=α−1/2`, `b=v+1/2`, the gap fraction is asymptotic to `B(a+1,b)/B(a,b)=a/(a+b)=(α−1/2)/(α+v)`. As a separate check on the **total** later-to-earlier moment ratio, apply the one-label equivalent at `v+1` using the same `α`; the exact shift is `B(a,b+1)/B(a,b)=b/(a+b)`. This gives a relative equivalent `(v+1/2)/(α+v)` even when `α/v→∞` and the total ratio vanishes. The two beta fractions add to one; neither small fraction should be read by subtracting the other's asymptotic equivalent from one. These are signed gap and positive ratio conclusions for the original law, not for independently sampled outcomes of the two labels.

The measure `ν` is finite on `(0,∞)`, of mass `A(0)/E[X]`. On each compact positive interval the fundamental theorem gives `ν((s,t])=(A(s)−A(t))/E[X]=∫_s^t (−A')/E[X]`; exhaustion shows that it has **no singular part on `(0,∞)`**. It excludes zero by definition. Continuity of `A` and `B(0)=(1−μ)A(0)` imply `yA'(y)→0` at zero, and the density is integrable there; they do not imply a bounded density or an equivalent at zero. They also say nothing about the mass of zero-endpoint births in an unrestricted intensity. At infinity the density equivalent gives the square-root blow-up of the Laplace integral at real `t=1` from below: split at fixed `Y`, bound the tail by relative multiples of `y^(−1/2)e^(-(1−t)y)`, and substitute `z=(1−t)y`. The compact part stays bounded. The integral is finite for every real `t<1` and diverges for `t≥1`, with leading coefficient `κ√π_math/E[X]`; this is a singularity of the transform, not of the measure.

## Evidence

The companion exact fixture (`eventual-domination-sol.py`) checks the prebirth/postbirth beta tilt, both moment beta ratios and the reciprocal beta mean. Rational endpoint mutations catch shifting the upper-edge power, losing the sign through cancellation of one-label errors, and dropping the quadratic condition from the matching diagonal. These are algebraic controls, not proofs of any of the three frozen probabilistic estimates.

```text
python3 -IB round2/checks/eventual-domination-sol.py
PASS 42 exact two-size placement/tilt checks
PASS 20 exact two-power-weight beta ratio and inverse-mean checks
PASS 10 rational overlap and endpoint controls
```
