# Sharp translation and regularity of the Q4 density from Fable Claim 6

This is a **consequence of Fable Claim 6**, not a replacement for Astra's independent [fixed-cut](../planted/README.md) or [retry-integral](../retry-integral/README.md) proofs. Use Fable `27a34b2aabbb8005160f4de00140abd1ef842562`, `round2/newground/law-of-X/density/REPORT.md` (Git blob `3195997ac2f153c0e1e52da4da5e99af38905564`) with its companion repair (blob `6a6170f33c51af39b636aff7497bee37b0be8b70`); see [Astra's outgoing audit](../../../../attacks/fable6-density/README.md). The bare-triple silence repair and Fable's original exponent bracket for its moment budgets are part of the input. No Gold2 or arithmetic sensitivity result is used.

Set `α=1/μ−1∈(0,1)`, `p*=1/(1−α)`, `A=κ/μ>0`. Extend the continuous density `f` by zero on the negative axis. For `h→0`, the full-line first differences satisfy

| Range | Equivalent for `||f(·+h)−f||_p` |
|---|---|
| `1≤p<p*` | `|h| ||f′||_p` |
| `p=p*` | `Aα |h| (log(1/|h|))^(1/p*)` |
| `p*<p<∞` | `A C_(α,p) |h|^(α+1/p)` |
| `p=∞` | `A |h|^α` |

Here `C_(α,p)^p = ∫_R |(t+1)_+^α−t_+^α|^p dt = 1/(1+αp)+∫₀∞[(t+1)^α−t^α]^p dt`. In the first row the constant is the norm of the **whole** derivative, rather than a constant determined by the cusp alone.

For all `1≤p≤∞` write `β_p=α+1/p` with `1/∞=0`. For the usual inhomogeneous Besov spaces of positive smoothness and `1≤q≤∞`, the exact result is

```
f ∈ B^s_(p,q)(R)  iff  0<s<β_p  or  (s=β_p and q=∞).
```

In particular **`f∈B¹_(p*,∞)` despite the logarithm in the first differences**, while `f∉W^{1,p*}`. For `1≤p<∞`, in the positive-order Sobolev–Slobodeckij convention (weak derivatives at integers), `f∈W^{s,p}` iff `s<β_p`. For `p=∞`, the endpoint holds in the Slobodeckij/Hölder convention `W^{s,∞}=C^{0,s}` for `0<s<1`, so `0<s≤α` is exactly the range. We do not use that notation for Bessel-potential spaces.

## 1. Inputs at zero and infinity

Fable §§5, 9 give `f(x)~A x^α` and `f′(x)~Aα x^(α−1)` as `x↓0`. These are two distinct proved equivalents. In §9 the scaled exit densities `v_j=f_(V_j)` have a uniform global `C^{1,1/20}` bound, and a uniform exponential right-tail probability bound upgraded to `|v_j(z)|+|v′_j(z)|≤C exp(−cz)` for `z≥1` (enlarge the constant from its stated large-`z` range). Their star mixture is

```
f(x) = Σ_(j≥3) 2/[j(j−1)] · j^μ v_j(j^μ x).
```

The derivative series converges on every compact positive interval by the same exponential bound. For `x≥1`, use `j^μ x≥(3^μ x+j^μ)/2` to factor exponential decay in `x` and a summable exponential in `j^μ`, even for the derivative weights `j^(2μ)`. Hence `|f(x)|+|f′(x)|≤C′ exp(−c′x)` at infinity. In particular no derivative tail control is inferred merely from the law's exponential moments or interior smoothness. At zero the derivative is locally integrable, `f(0)=0`, and integration of the classical derivative across zero shows that the zero extension is absolutely continuous on compact intervals and its weak derivative is `1_(x>0) f′(x)`, with no atomic term. It belongs to strong `L^p` precisely for `p<p*`, by the derivative equivalent and the checked infinity bound.

## 2. First-difference equivalents

By translation invariance of the norm take `h>0`. For `p<p*`, the difference quotient is the average of `f′(·+u)` over `0<u<h`. Strong continuity of translations in finite `L^p` gives convergence to `f′` in `L^p` and the first row. The derivative norm is nonzero since `f` is a nonzero density vanishing on a half-line.

For `p>p*`, scale `x=ht`. On any bounded set in `t`, `h^(−α) Δ_h f(ht)` converges uniformly to `A[(t+1)_+^α−t_+^α]` by the density equivalent. In `R<t<ε/h` for `R≥1` and fixed small `ε`, the derivative bound and its integral give `|h^(−α) Δ_h f(ht)|≤ C t^(α−1)` (use `x+h≤2ε` for small `h`). Its `p`th power is integrable when `p>p*`; for `p=∞` this bound tends to zero as `R→∞`. On `x≥ε`, the differences are `O(h)` in `L^p` and `L^∞` by the derivative bound on that half-line. After normalization this part vanishes since `α+1/p<1` (or `α<1` in supremum). The region `t<−1` contributes zero. These bounds prove convergence of the scaled norms to the model norms in both cases. The model increment has supremum one: on `−1≤t≤0` it increases to one, and on `t>0` concavity makes it decreasing.

At `p=p*`, the region `−h≤x≤Rh` has `p`th integral `O_R(h^(αp+1))=O_R(h^p)`, and `x≥ε` gives `O_ε(h^p)`. In the intermediate region `Rh≤x≤ε`, for any fixed small `η>0` choose `ε` so the relative derivative error is at most `η` up to `2ε`. Since the derivative model decreases, integration over `[x,x+h]` bounds the increment divided by `Aα h x^(α−1)` between `(1−η)(1+1/R)^(α−1)` and `1+η`. Now `p(α−1)=−1` and `∫_(Rh)^ε dx/x=log(1/h)+O_(R,ε)(1)`. First send `h→0`, then take `R→∞` and `η→0` in the bounds. This gives the coefficient of the critical logarithm, not just its order.

## 3. Higher differences and the Besov endpoint

Use fourth forward differences `Δ_h^4=(τ_h−I)^4`. A bound at `p≤p*` of the required order cannot be obtained just by taking higher differences of the first-difference estimate. Instead reuse the **four-pair** construction of the pinned Fable §§1–4 **for the scaled exit laws**. In each retained piece, condition on the history, marked paths and pair totals. The predecessor mask is shape-only on that piece, so the four within-pair Beta variables are independent there conditional on this data. With Fable's notation, `T_i` is pair mass, `B=min(1,a_1,…,a_8)`, `P=min_i T_i`, and `λ=j^μ P^μ B`. Each pair is a scaled Beta image; its coefficient sum is at least `λ` and the Beta shape parameters are the fixed `d,d`. The [angular translation kernel bound](../planted/README.md#1-a-derivative-kernel-in-the-physical-split-coordinate) (a deterministic Beta-image bound) gives total variation of its first difference at most `C_d min(1,√(|h|/λ))`. The remainder is a translation after conditioning and costs nothing in total variation.

For convolution, `Δ_h^4(K_1*…*K_4)=(Δ_h K_1)*…*(Δ_h K_4)` exactly. Young's inequality for finite signed measures therefore bounds the fourth-difference total variation by `C_d min(1,|h|² λ^(−2))`. Apply Fable's pre-mask moment and probability estimates **at scale power 2**, just as for its weighted exit density bound: Beta pair-total inverse moments first, then Hölder with `9/8,16,144/7` for size, inverse coefficient and disjointification indicator. The actual size exponent `(9/8)2μ<4/3`, inverse order `32<40`, and exponential miss in core count pay for the sum of pieces uniformly over `j`. Explicitly, the weighted sum of `λ^(−2)` over the retained pieces is uniformly finite; this does not require that the shape laws factor after the predecessor mask. The pieces account for mass one, so also use the trivial `2^4` bound. Thus

```
sup_j ||Δ_h^4 v_j||_1 ≤ C min(1,|h|²).
```

The star-mixture scaling in `L¹` cancels the Jacobian. Split at `j=|h|^(−1/μ)`; since `2μ>1`,

```
||Δ_h^4 f||_1 ≤ C Σ_(j≥3) j^(−2) min(1,|h|² j^(2μ))
               ≤ C |h|^(1/μ) = C |h|^(1+α),     0<|h|≤1.
```

The supremum first-difference bound gives `||Δ_h^4 f||_∞≤C|h|^α`. Interpolating these **higher-difference** bounds yields `||Δ_h^4 f||_p≤C|h|^(α+1/p)` for every `p` in the table. For the reverse bound, on `−4h<x<−3h` with `h>0` only the rightmost term is nonzero, exactly `f(x+4h)`. The cusp equivalent gives a lower bound `c h^(α+1/p)` on the norm for finite `p`, and `c h^α` for `p=∞`. The same bounds therefore hold for the fourth modulus `sup_(|u|≤h)||Δ_u^4 f||_p`.

Since `β_p<2`, the standard fourth-difference Besov norm for `0<s<4` reduces its small-scale test to integrability of `t^(β_p−s)` in `L^q(dt/t)` (or boundedness if `q=∞`). Higher smoothness is excluded by the ordinary nestedness of Besov spaces, already below order four. For finite `p` and noninteger `s>0`, use `W^{s,p}=B^s_(p,p)` in the Slobodeckij convention. At integer one the weak derivative test in §1 applies, including nonmembership at `p=p*`; higher integers are excluded by nesting into any lower noninteger order strictly above `β_p`. For `p=∞` and `0<s<1` use the Hölder convention and the first difference (or the equivalent Besov norm there); integer orders cannot hold because they would imply Lipschitz continuity. These claims are for the **zero extension** on the whole line, not an assertion about pointwise asymptotics of higher classical derivatives on the positive half-line.

## 4. An observable inverse for μ

With the actual limit law, use its continuous density and measure `ω(h)=sup_x|f(x+h)−f(x)|`. Then `lim_(h↓0) log ω(h)/log h=α`, and `μ=1/(1+α)`. Alternatively the CDF gives `lim log P(X≤x)/log x=1/μ` and is determined without a choice of density representative. For a finite `p`, the log-slope of fourth differences is `α+1/p` even when the first-difference slope is clipped at one. These are identifiability statements for the exact limiting law; they do not provide a finite-sample, finite-tree, or arithmetic estimate of the exponent.
