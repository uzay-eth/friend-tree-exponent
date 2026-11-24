# Cold review: conditional upper-range comparison

**Holds at 256; the same estimates also certify 232.** Reviewed sharp-large (`open/age-order/global-transport/sharp-large/README.md`) blob `0c81b9c15fb56a4efefe3feecd243fa70359811a`. This is a pointwise statement about *pre-placement* trees: for every old tree at size `n`, the next two labels have strict conditional tail order for `0<x<1` and `(n+1)x≥232`. Averaging gives that same upper-range conclusion for the original law at every eligible label. It does not give the unconditional all-label claim at lower thresholds, does not compare arbitrary post-placement endpoints, and says nothing at hubhood `x=0`. At `x=1` both tails are zero.

## Sender colours and conditional losses

Only the serial middle's fresh leaf starts red. Colour follows the initiator, even when the target has the opposite colour; the mark is permanently blue. The probability of a word is the Pólya word probability since the next colour rate depends only on the colour counts. Conditional on that word the within-colour sender and neighbour kernels are unchanged, as can be verified first on finite words then by the sequential construction with a beta parameter. In particular no conditional change to actual old-parent arrows is being made.

If a red sender is a good leaf, its new vertex is good. If a bad red sender adds a good leaf it only helps. Thus with `B,R` the populations and `F` the bad red fraction, the red conditional drift is at most zero, while the blue drift is at most `1/(BR)`. Indeed only the mark can promote a leaf adjacent to itself; that event loses at most one and occurs only on a blue step. One may instead permanently mark a lost red lineage as bad, in which case the red drift is exactly zero and the blue drift is nonnegative; the upper charge remains the same. The finite expectation argument needs neither a limit of `F` nor independence from red abundance. Under a fixed tilted word law the red share tends to `θ`; with `F≤1`, replacing the share in its product with `θ` costs `o(1)` in expectation. Fatou then gives the displayed limit bound on `E[(θ−Z)_+]` for beta-almost every `θ`, using convergence of the physical degree share. The optional full-loss identity (2) includes red waits before blue jump zero; summing their generating functions gives `u/(1−θs−ut)`, consistent with the printed log integral.

Given the entire colour word, the event of avoiding the blue mark during the first `k` blue selections has probability `(b−1)/(b+k−1)` with `b=n+1`. This follows stepwise even under previous avoidance: the mark is one of exactly `b+f` blue senders. Therefore the word law is unchanged, including the red string before the first blue choice and the string just before blue jump `k`. No red losses occur during the protected prefix. Beyond it the fresh choice drivers retain their conditional kernels. For `f≥1`, the reciprocal charge uses `R_f=1+Σ_(i≤f)G_i`, not the start-of-block count, giving `E[1/R_f]=u(1−u^f)/(θ f)`. The tail integral from `k−1` gives `E[(θ−Z)_+|E]≤(u/b)log(1+b/(k−1))`. The bounds hold uniformly in the old tree, including no other leaf at the parent.

## Adverse weights and numerical bound

Use the exact signed identity at a common threshold in the common final size. Put `s_u=r_u/[D_u(D_u+1)]` and `A=Σ_u D_u s_u`; then `Σ_edges |s_u−s_v|≤A`. Count only the adverse endpoint per edge, and separately pay the serial tips, for a total upper budget `2AU`. A leaf with no other siblings at its core parent still has locked red mass at most `1/2`; monotone reciprocal deficit gives the beta upper tail without stock assumptions. The checker includes a nonzero adverse edge on which the charged parallel leaf has no siblings. The gamma ratio in the upper beta integral follows from `Γ(q+1/2)/Γ(q)≤√q` with `q=m−1/2`, so the printed upper bound `U≤(1−x)^b/√(π b x)` follows.

The printed floors at `b≥256` are conservative and yield `π·256·c²>4`, with `c` the lower coefficient after the shift. To tighten using the same argument, let `k=ceil(11b/5)+1` and `t=5/6`. The prefix payoff is at least `(b−1)/(16b/5+1)`; the loss integral is at most `(u/b)log(16/11)`. For `θ≥x+t(1−x)/b`, the gap is at least `tu/b`, and the protected Markov success floor is `1−log(16/11)/t`. Its beta tail integral contributes the factor `(1−t/b)^b(1−x)^b`. By the logarithm series,

```
(1−t/b)^b ≥ exp(−t−t²/[2(b−t)]) ≥ exp(−t) [1−t²/(2(b−t))].
```

Both `b`-dependent factors in the resulting coefficient lower bound are increasing. Evaluate at `b=232` using odd alternating partial sums for the lower `exp(−t)` and the upper `log(1+5/11)`, and `π>157/50`. The rational replay gives `(157/50)·232·c²>4`, sufficient at `bx≥232`. This is only a modest optimization of prefix length and endpoint shift, not a lower bound on the best cutoff and not a claim of a single-digit cutoff.

## Independent bounded replay

```text
taskset -c 8 nice -n 10 timeout 8s python3 -IB round2/checks/preplacement-256-sol.py
PASS 102424 masked colored histories with reciprocal drift and exact prefix mass
PASS 12 geometric reciprocal brackets and nonzero adverse zero-sibling endpoint
PASS original 256 and tightened 232 cutoffs with rational bounds
```

The infinite-time implications are analytic; the checker only enumerates depth-four words and arrow histories and checks finite arithmetic.
