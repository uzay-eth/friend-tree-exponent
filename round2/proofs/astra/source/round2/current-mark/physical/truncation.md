# Positive size truncation of the cap-20 kernel

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 13: Claude for transfer bookkeeping; GPT for the uncapped bound and root transfer; inputs include physical jump geometry and the Round 1 exponent and prior. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Take the residual of [the current-mark construction](README.md) and kill each particle on its own first intrinsic jump from size `N` to `N+1`. Size is the number of real vertices in its shape. The other residual transitions never increase lineage size: `Y` and `Z` are strict subtrees, `Y⁺` has at most the size of its source, and a dormant switch preserves size. The direct reward `B₀` is collected at the active source. In particular killing a particle at the boundary does not kill its previously emitted relatives.

Let `R_N` denote the killed first-moment generator, and `K_N(s)=J(s−R_N)⁻¹B₀`. The diagonal retains the full intrinsic rate `n−1/D`, including at size `N`, plus active killing `a` or dormant total departure `u`. Separate intrinsic transitions, active emissions, and dormant activations to get a nonnegative off-diagonal matrix. For rational `s` in the convergence domain below, `K_N(s)` is rational. Its finite inverse exists there: the killed expected active occupation is bounded by the unlimited one, and the killed dormant occupation is bounded by active occupation because each dormant has a uniformly positive chance of switching and surviving for a bounded interval. More precisely, for any dormant `(l,u;S)` of size at most `N`, require activation by elapsed time `1/2`, no intrinsic jump through elapsed time `2`, and no death before or after activation through elapsed time `2`. This event has a positive lower bound depending only on `N` and the finite label set. The resulting original is active throughout `[t+1,t+2]`; its discounted occupation bounds a constant times the dormant count at `t` after taking expectations. Finite discounted occupation for the entire finite matrix implies that its resolvent is the Laplace integral. By the common-family coupling and nonexplosion,

```
0 ≤ K_N(s) ≤ K_{N+1}(s) ≤ K_c(s),     K_N(s) ↑ K_c(s).
```

The exact error is a first-exit contribution: for `g≥0`, sum `exp(−s τ) U_s(g;X)` over boundary particles `X` at their exit times `τ`, in the killed family with exits recorded, where `U_s(g)=(s−R_c)⁻¹B₀g`. In the coupled unlimited family, each lost contribution belongs to the subtree of the first boundary crossing on its lineage. No copy or activation term is lost except through such a recorded exit. This formula follows by conditioning at the recorded exits and Tonelli.

## A qualitative uniform rate

Fix any real discount interval `[s_min,s_max]` with `s_min>s₀>0`, where `s₀<μ` is supplied by [gap.md §6](gap.md). For all discounts there, the gap gives

```
0 ≤ U_s(g;X) ≤ C ||g||∞ W_c(X) ≤ C' ||g||∞ n(X).
```

Late exits, with `τ>T`, contribute at most the unlimited discounted reward after time `T`, which is bounded uniformly on newborns by `C'' ||g||∞ exp(−(s_min−s₀)T)`. To bound early exits, retain every boundary particle as an inert record at its new size, and give every live particle and record weight `n²`. The total weight has drift at most `32` times itself. Indeed the intrinsic size-increase drift of each live particle is at most `n(2n+1)≤(5/2)n²`. Active child-copy terms have total cost at most `(9/8)n²`: use `n(Y⁺)≤3n(Y)/2`, the summed copy coefficient `q(D)≤1/2` per child, and disjoint subtrees. For root-drop emissions, `r≤D≤19` when the drop is nonzero, and the child sizes have square sum at most `n²`. Grandchild-drop emissions cost at most `n²/2` by the same argument and `p≤1/2`. Switches preserve weight and deaths decrease it. Exit recording preserves precisely the weight payment of an intrinsic size increase. The stopped drift bound and Fatou yield expected total weight `≤4 exp(32T)` from any newborn. Every early exit has `n=N+1`, so its expected weight `n` summed over records is at most `4 exp(32T)/(N+1)`. The nonnegative discount does not increase this bound.

Take `T=log(N+1)/(32+s_min−s₀)`. The resulting uniform entrywise (or bounded nonnegative action) error is

```
0 ≤ (K_c(s)−K_N(s))g ≤ C_* ||g||∞ (N+1)^(-η) 1,
η = (s_min−s₀)/(32+s_min−s₀) > 0.
```

This establishes a power-law size-tail bound. The proof does not give numerical values for `s₀` or `C_*`, so the bound does not yet turn a rational killed result into a numerical upper enclosure for the infinite matrix.

## What can be certified with finite blocks alone

For any rational `s>s₀` at which `r(K_c(s))>1`, some finite killed solve certifies `r(K_N(s))>1`. This follows from entrywise convergence and continuity of the finite matrix radius. A rational positive action `g` with `K_N(s)g>g` then exists: start with the positive Perron vector for a sufficiently large cutoff (every label pair is connected by a finite star path), and approximate inside the strict inequalities. Conversely, any such action for a killed kernel at a positive discount certifies `s<μ`: if `s≥μ`, its finite rewards are bounded by those of the full kernel, contradicting `r(K_c(s))≤1`. A killed result below one does not give the upper endpoint; an effective tail bound or a separate positive full-operator upper witness is needed.
