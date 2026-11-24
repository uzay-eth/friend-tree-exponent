# Independent check: large-cutoff top-component age law

**Verdict: the proposed law holds for the ordinary decorated component rooted at its highest surviving core, with the external parent removed.** The retained root's original degree may still count the one exterior slot, as it does in the inert model. Do not include the degree or mark of the actual high parent (nor the top's permanent promoter-degree mark), or assert this approximation under sampling proportional to component size. All limits below take the original size to infinity at each fixed cutoff first, and then send the cutoff to infinity. In fact the epsilon-fraction discard is unnecessary: a birth after `sqrt(log N)` has enough slack for both the no-input coupling and the descendant-size cutoff, uniformly on a bounded parent-age window.

## Measure and target

Use Claim 10's one-ended marked core law `π`, with parent toward the end. In the fixed-threshold count formula (`low-forest/README.md`), count each component containing a surviving core at its unique top. By mass transport along the oriented edge, the same (unnormalized) test mass can be counted at its high parent. For `|f|≤1` on finite top-rooted ordinary trees, suppress the founder exception, and write

```
M_N(f) = π[ 1{D_o>N} Σ_(core children v of o) 1{D_v≤N} f(F_N(v)) ].
```

`F_N(v)` is the downward component after cutting vertices by their original degrees and restoring their surviving ordinary leaves. All tests may depend on `N`; the error bounds below hold simultaneously for this whole unit ball. Tests on ordinary types may equivalently retain degree labels reconstructed internally plus the single exterior root slot whenever there are no further cuts. The count formula is valid for signed tests because `π c_o<∞`. This measure does not include singleton ordinary leaves; if singleton surviving cores are also excluded their mass is `o(M_N(1))` by the conclusion, since an inert virgin has two ordinary vertices and never loses them.

Let `J` be the inert virgin (one leaf, one parent slot), write `P_u^0 f(J)` for the physical *intrinsic* planted law (not the positive residual branching family), and put `L=log N`. The stronger unnormalized conclusion is

```
sup_|f|≤1 | (N^μ/L) M_N(f) − E_β[W^μ] (1/L)∫_0^L P_u^0 f(J) du | → 0.       (1)
```

Here `W=lim e^(−t)D_t` for the focal newborn grown in its real environment from `β`, the stable-age birth law on full environments. In particular `M_N(1)∼E_β[W^μ] N^(−μ) log N` with positive coefficient. Dividing by this mass gives total-variation approximation by the inert law at an independent uniform age on `[0,log N]`; (1) gives full sum-of-absolute-masses control, without a factor `1/2` convention.

## 1. Stable age, old/new high event and integration tails

Apply the stable-age identity (`degree-law.md`) `π=μ∫e^(−μt)β S_t dt` to the nonnegative parts of the test, using monotone truncation for the unbounded child count. Under `β` the focal core starts with no core children. If `B_t` is its number of promotions to time `t`, then `c_t=1+B_t` and `B_t` is dominated by a Poisson count of rate one. With `t=L+y`, the scaled integral is `μ∫_(−L)^∞ e^(−μy)` times the expectation of the summand divided by `L`. Write `H_N(y)={D_t>N}`, and for fixed bounded `y` let

```
b=√L,    E_N(y)={e^(−b)D_b>e^(−y)}.
```

For large `N` on a bounded `y` interval, `0<b<t`. For every `y` with `P_β(W=e^(−y))=0`, both event indicators converge almost surely to `1{W>e^(−y)}` by the same supermartingale limit, even though the two observations come from the same path. Only countably many `y` are exceptional. The replacement error divided by `L` is at most

```
(E_β[B_t²]/L²)^(1/2) P_β(H_N(y) △ E_N(y))^(1/2).
```

This tends to zero pointwise almost everywhere and in the bounded-window integral by domination. One does **not** have to prove an asymptotic law for the number of promotions conditional on `H_N(y)` to do this replacement.

For the original integral outside the window there is a uniform integrable majorant. The Yule upper comparison for degree gives `E_β D_t²≤6e^(2t)` and hence `P(H_N(y))^(1/2)≤√6 e^y`. For `−L≤y<0` and `L≥1`, `(E_β B_t²)^1/2/L≤(t+1)/L≤2`, so the scaled integrand in absolute value after leaving off `μ` is at most `2√6 e^((1−μ)y)`. For `y≥0` use `(t+1)/L≤2+y`, giving `(2+y)e^(−μy)`. Set the negative integrand to zero below `−L`. The proposed limit coefficient's integrand has the same tail integrability, since `E W²≤6` by Fatou. It is enough to work with `E_N` only on each bounded window; no uniform relative error in that event at large negative `y` is needed.

## 2. Conditional child law after replacing the high event

Fix `A<∞` and `|y|≤A`. Discard all births before `b`; their expected contribution divided by `L` is at most `b/L`. For any retained birth time `τ∈[b,t]`, the event `E_N(y)` belongs to the full history at `τ`, and on it

```
D_vparent(r) ≥ e^(b−y)  for every r≥τ
```

by monotonicity of degree. `vparent` here is the focal high-parent candidate, not the promoted child. On `E_N` it need not yet have degree above `N` at observation; the high event was already replaced, and the summand to compare is always the **downward** clipped child fringe ignoring the parent. This formula for the summand is defined on all configurations, and agrees with the full component whenever `H_N` holds and `D_v≤N`.

Condition on the entire history at `τ` just after the promotion. The child has one leaf and an exterior parent slot, regardless of the parent's history. Drive it and an inert `J` by identical fresh physical drivers until the first incoming parent arrow; the child's own choices of its parent slot, which add leaves outside the branch, are nonmutating in the branch and stay in the physical law. Only incoming parent arrows can break agreement. On the measurable event above, the predictable incoming intensity is at most `e^(y−b)`, even through arbitrary feedback from the child to the parent. Thus the conditional failure probability through observation is at most `(t−τ)e^(y−b)`. Neither the no-hit event nor any further event involving the parent's future is conditioned on in the expectation. In particular one must not claim that the endpoint is inert *conditional on no incoming hit*; one compares unconditional marginals of a coupling and charges its failure probability.

The physical branch ordinary size `Z_u` after birth has gain rate `Z+1/D_parent−1/D_child` (all internal ordinary initiators counted once, minus the arrow out, plus the arrow in). It is bounded above by the rate `Z+1/2` and starts at 2, so `E[Z_u | F_τ]≤(5/2)e^u`. Whenever `Z_u≤N`, no vertex of this whole branch can have original degree `>N`: the top has at most `Z_u−1` internal neighbors and one parent, the rest have only internal neighbors. Thus the missing/tested-downward-root and further-cut error conditional on birth is bounded by `(5/2)e^u/N`. For all `τ≥b`, `u=t−τ≤L+y−b`; the bound is at most `(5/2)e^(y−b)`. In particular a separate removal of the oldest `εL` ages would be valid but is not needed. Multiplying the probabilities by 2 handles arbitrary signed tests in the unit ball. The expected count of retained opportunities is at most `t` (or use the compensator below), so their scaled total error is `O_A((L+1)e^(A−b))=o(1)` after division by `L`. These arguments use only the strong Markov property at actual births; siblings need not be independent.

## 3. Compensating birth times without selection bias

With the endpoint test now replaced by the deterministic function `g(u)=P_u^0 f(J)`, its retained expectation is exactly

```
E_β[1_E Σ_(b<τ≤t) g(t−τ)]
 = E_β[1_E ∫_b^t (1−c_r/D_r) g(t−r) dr].
```

`1_E` is known at `b`, and the pre-jump promotion rate is predictable; deterministic endpoints cause no additional atom. This step would not be valid with `1_H` inside the integrand in place of `1_E`. On `E`, the error of replacing the rate by 1 before dividing by `L` is at most

```
e^(y−b) ∫_b^t E_β c_r dr ≤ e^(y−b) (t+t²/2),
```

using the focal Poisson domination. No independence of `c_r` and `E` is used. The remaining age integral is `P(E)∫_0^(t−b)g(u)du`, and changing its endpoint to `L` costs at most `b+|y|`. Hence over each bounded window, uniformly on the whole unit ball, the scaled law is `(1/L)∫_0^L g(u)du` times `μ∫ e^(−μy)P(E_N(y))dy`, up to a vanishing error. The bounded-window limit, followed by the tail bounds in Section 1, proves (1) because

```
μ ∫_R e^(−μy) P(W>e^(−y)) dy = E W^μ.
```

There is no correction to the coefficient from the long-time promotion deficit; the argument estimates it directly on the early-high measurable event.

## 4. Logarithmic volume and scope

For the inert `J` write `Z^0(u)` for ordinary branch size omitting the exterior slot. The rate bounds `Z−1/2≤rate≤Z+1/2` permit whole-path stochastic comparison with the two linear pure birth chains started at 2, by common marked proposals at equal states. Both comparison chains satisfy `e^(−u)Z_±(u)→W_±` with finite strictly positive gamma limits (shape `2±1/2`), as in the ordinary descendant-size argument (`fringe/other-marginals.md`). In particular `(log Z^0(u))/u→1` almost surely in the coupled law, or in probability from marginal comparison alone. Represent the independent uniform age as `u=LV` with `V∼Unif[0,1]` independent of the whole inert path. For almost every `V>0`, `log Z^0(LV)/L→V`. Push forward the total-variation approximation (even though the pushforward map varies with `N`) to obtain `log |F_N| / log N ⇒ Unif[0,1]` under component-count normalization. The size here is ordinary volume, not number of cores; a core-count exponent would require separate input.

The relative count asymptotic and total-variation statement are for threshold limits *after* the fixed-cutoff local limit. This argument gives no diagonal approximation for a threshold depending on finite observed size. If stubs to additional high descendants are retained as marks, they have vanishing mass under the same cutoff comparison; the single parent stub is always present at a top, but its actual external degree is necessarily `>N` and cannot be compared to an inert mark with finite degree. No arbitrary-age joint-cohort independence is asserted.

## Independent physical fixture

`timeout 8s taskset -c 7 nice -n 10 python3 -B independent.py` gave `1930 physical core-edge branches; 1192 unequal-flux cases; all PASS`. The fixture enumerates directed arrows on labelled finite trees, checks the sign and direction of the two boundary fluxes, and checks exactly which branch vertex has an external slot. It checks the local identities; the total-variation and domination arguments above are analytic.
