> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../../README.md#claim-17) for exact scope and citations.

# Q4: fractional smoothing of the full law with the planted future weight

**The full density satisfies** `||f(·+h)−f||₁ ≤ C |h|^(13/278)` for `|h|≤1`, hence belongs to `W^{s,1}(R)` for every `0<s<13/278`. The proof keeps the exact future no-cross likelihood on deterministic physical cuts. Local boundedness, continuity and unimodality remain open here.

## 1. A derivative kernel in the physical split coordinate

Use the two-side cuts and the retry budget in the [audited angular argument](../../../../../../../../reviews/cold/astra/checks/limit-iteration-sol.md). Conditional on the actual cut history at time `t`, and on the two embedded side paths in the independent rate-size reference, the rescaled physical amplitude on the no-cross sublaw is

```
g(p) = a p^μ + b(1−p)^μ,     p ~ Beta(n₁,n₂),
S = a+b > 0,                N = n₁+n₂.
```

The initial leaf gate gives `n₁,n₂≥2`. The exact terminal likelihood `L` depends on embedded paths only. All estimates below hold for fixed paths; `L` is retained when taking their average. Use `p=sin²θ`, with angular density `w` on `(0,π/2)`, and `F(θ)=g(sin²θ)`. The angular density is unimodal, vanishes at the endpoints, and has both supremum and total variation bounded by `C√N` by the [angular Beta bound](../iteration/README.md#1-a-sharper-two-side-kernel). These facts concern the *input* density only.

For `0<ε≤ε₀` and `H=c_μ S ε` with a sufficiently small fixed `c_μ>0`, discard the two angular endpoint strips of length `ε`, and discard `|F'|≤H` on the rest. Their total angular length is `≤C_μ ε`. Here are the geometric details (including zero coefficients):

* On `F<S/4`, `g(1/2)≥S/2` implies `|g'|≥S/2` by concavity and the secant to `1/2`. After endpoint trimming `sin(2θ)≥C ε`, so this region misses the derivative strip if `c_μ` is small enough.
* On `F≥S/4` and inside the derivative strip, use
  `F''=−4μ(1−μ)F+(2μ−1)(cot θ−tan θ)F'`.
  On the trimmed interval the last term is at most `C |F'|/ε`, so `F''≤−c'_μ S`. The high region is an interval, and the derivative can cross each strip boundary only downward there. The strip is thus an interval (possibly empty) of length at most `2H/(c'_μ S)`. It cannot leave and re-enter through the low region because the high region itself is an interval and the derivative-strip points there form one interval.

For the good part of the angular interval, the pushforward has an integrable BV density, with variation `≤C_μ√N/H`. To see the point that makes this quantitative, the inverse Jacobian `1/|F'|` has bounded variation `O(1/H)` on the good intervals, **without any uniform lower bound on the critical-point curvature on the whole angular domain**. In fact `F''` has at most two zeros in the open domain. With `k=2μ∈(1,2)`, `x=tan²θ`, and `c=k−1>0`, its sign is the sign of

```
Q(x) = a x^(μ−1)(c−x) + b(cx−1).
```

If `a=0`, this is affine. Otherwise

```
Q'(x)  = bc − a[(1−μ)c x^(μ−2) + μ x^(μ−1)],
Q''(x) = a(1−μ) x^(μ−3)[(2−μ)c + μ x] > 0.
```

A strictly convex function has at most two zeros. Cut the good region at those zeros; on each remaining interval both `F` and `|F'|` are monotone, with `|F'|≥H`. The inverse Jacobian varies by at most `1/H` per interval. The product of this Jacobian with the input density has variation `≤C√N/H` by the product inequality for BV, including boundary jumps when extending by zero. Changing from `θ` to `F(θ)` preserves variation on each monotone branch, and adding the finitely many contributions gives the claim.

Let `K` be the fixed-path probability kernel of `F(θ)` and `τ_h` translation by `h`. The discarded kernel has mass `≤C√N ε`; translation costs at most twice this mass. For the good kernel, apply `||τ_h v−v||₁≤|h|·TV(v)` to its BV density. Optimizing with `ε=√(|h|/S)` gives, for all real `h`,

```
||τ_h K−K||_TV ≤ C_μ √N √(|h|/S).
```

Here TV denotes the total variation norm of a signed measure (the `L¹` norm for densities). For too large a choice of `ε` use the trivial bound two and enlarge the constant. The estimate remains valid if one coefficient vanishes, and does not require a negative moment of either individual coefficient.

## 2. Averaging with the exact future likelihood

Select the edge measurably at deterministic `t` using only the actual leaf stocks and the threshold from the angular retry proof. Write `ν_t` for the sublaw on histories with an eligible edge and **no subsequent crossing** of that edge. The discarded mass, including ineligible histories, is at most `C exp(−ct)` with `c=13/252`. Requiring no crossing at one cut time is an exact shape likelihood for that reference, not a restriction on the Beta clocks. There is no subtraction of earlier attempted no-cross events inside the integrand.

In the reference conditional kernel bound, integrate over shape paths with the exact terminal factor `L`, then over the ordinary physical cut-time histories. Since `0<X=g(P)≤S` on that sublaw, and the [physical left-tail bound](../../left/README.md) gives `E X^(−1)<∞`,

```
∫ √N S^(−1/2) L dQ d(history)
  ≤ E_phys[√N X^(−1/2); retained]
  ≤ (E_phys |T_t|)^(1/2) (E_phys X^(−1))^(1/2)
  ≤ C exp(t/2).
```

In the first line the clock integral can be included because `S` is shape-only, and the full physical `X` appears only *after* transferring the exact weighted sublaw. Dropping `L` before this transfer would leave an unjustified planted inverse moment. The unresolved future paths need no density or kernel bound: both their law and its translate have the same missing mass. It follows for the actual full density `f` that

```
||τ_h f−f||₁ ≤ C exp(−ct) + C exp(t/2)√|h|.
```

Take `t=log(1/|h|)/(1+2c)`. The exact common exponent is `c/(1+2c)=13/278`.

Hence `f∈W^{s,1}(R)` for all `0<s<13/278`. By definition the fractional seminorm is the integral of `||τ_h f−f||₁/|h|^(1+s)` over `h∈R`; use the modulus near zero and `||τ_h f−f||₁≤2` at infinity. One also gets an `L^p` version: for `1<p<r<278/265`, interpolate this translate against its uniform `L^r` bound `2||f||_r`. Then `f∈W^{s,p}` for

```
0<s<(13/278)·(1/p−1/r)/(1−1/r).
```

This follows from the translation form of the `W^{s,p}` seminorm applied to the interpolated translate norm. Letting `r` approach `278/265` gives `s<1/p−265/278`, positive for every `p` strictly inside the audited strong range. There is no endpoint `W^{13/278,1}` claim.

## 3. Why a clock randomization alone does not produce boundedness

An independent randomization of cut time can average the preceding exact retained/rejected sublaws, but for any single deterministic-time cut kernel the **two positive coefficients** already furnish a quantitative obstruction to a bounded *conditional* kernel, even with its correct shape likelihood. At the unique interior maximum `x_*=g(p_*)`, the Beta density is strictly positive at `p_*` and `g''(p_*)<0`. Taylor and the inverse-branch formula give

```
K([x_*−δ,x_*]) ~ 2 w_P(p_*) √(2/|g''(p_*)|) √δ.
```

For each fixed pair of side paths with `L>0` the exact weighted kernel has the same asymptotic multiplied by `L`; randomizing a common elapsed hold after the cut leaves the ratio unchanged. This is not a lower bound on the original mixture at a common `x_*`: those maxima depend on the shape paths. It identifies a genuine additional input needed for this route: quantitative control on the future-weighted distribution of that maximum (for example, with extra bounds, a uniform integrable inverse square root of the distance to a test location). Merely choosing a time independently of future paths does not change the angular fold conditional on cut history and embedded paths. The rate for the translation modulus by itself is also insufficient to infer local boundedness or unimodality: integrable compactly truncated spikes `|x−x₀|^(−β)` for `0<β<1−13/278` satisfy this translation rate, and disjoint translates and sums preserve it. This is an obstruction to the inference, not evidence for spikes or multiple modes in the physical law.

A finite adapted downcrossing (for example, the first jump at which a present-state reciprocal passes a fixed threshold) can be fed into the no-cross formula by the strong Markov property, also with independent randomization of the rule. But it cannot yield the **full** law as a single no-cross sublaw. At any selected directed slot `v→w` with finite degree `D_v=d`, consider the next ring of `v` and those of its `d` present neighbours. With probability `1/(d+1)` the first is `v`'s; its degree has then stayed `d` because only an existing neighbour can cause its next gain. The chance it chooses the selected edge is then `1/d`. Thus the conditional future-crossing probability is at least `1/[d(d+1)]` from any such history, regardless of external growth between those waits. This applies after conditioning on an independent choice of the cut rule as well. A retrospective downcrossing selector based on the full future amplitude has a different issue: a gate on `X=g(P)` explicitly restricts the Beta variable even when the paths and their likelihood are held fixed. For example, gating on `x_*−δ≤g(P)≤x_*` and renormalizing preserves the fold and has acceptance of order `√δ`; the conditional density must have supremum at least `1/δ`. One cannot treat that extra mask as shape-only. The argument above bypasses both issues for fractional regularity by bounding the total variation of the unresolved future paths at each deterministic cut.

## Finite controls

`nice -n 10 timeout 12s python3 -IB round2/open/limit-law/smoothing/planted/check.py` tests the angular turning polynomial against independently differentiated monomials, including one-coefficient degeneracies, and the exact interpolation and fold exponents. These are algebra checks; the variation and weighted-measure arguments are analytic.
