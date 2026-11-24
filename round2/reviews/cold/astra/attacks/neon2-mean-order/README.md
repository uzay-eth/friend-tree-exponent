# Outgoing audit — Neon-2 Claim 8 (mean degree shares)

**Verdict: HOLDS.** At `834bb960235f0d0afdc263da8af510f7c452997b`, the all-tree conditional mean bound and the edge-seeded founder bound both hold with their published constants. The result is strict order of the **means**, with equality between the founders.

## Source and replay

From the Astra root, run:

```sh
python3 -IB round2/attacks/neon2-mean-order/check.py
```

The checker reads `git archive <commit> round2/neon2/mean-order` from the Neon-2 object store (`[local path]` by default, overridable by the sole positional argument). It checks the archive's commit header and all three Git blob IDs. The archive is opened in memory, with no extraction, checkout, or changes to Neon-2. The certificate is executed from the archived bytes in isolated Python with and without optimization, and the two outputs agree byte-for-byte. Both pass 13,299 table rows, 24 certificate trees and 10 stock controls. The archived `AUDIT.md` is an earlier review labelled Claim 9: its eight-control count and `assert`/non-integer caveats concern the pre-hardening script. They do not describe the pinned executable. The pinned verifier checks positivity rather than pinning the printed numerical margins for arbitrary replacement tables; this audit pins the fractions separately.

## All-tree proof and founder setup

- **Continuous-time coupling.** At size `m`, the `Exp(m)` holding time is independent of the uniformly chosen caller and its neighbour mark. The holding-time sequence is therefore independent of the jump chain. Since `r(v)≤D_v`, the discrete degree share is a bounded supermartingale from the tag's birth. In continuous time, the size martingale `e^{-t}|T_t|` is L²-bounded and has mean `N` from a size-`N` start. Its limit `W` is independent of the limiting degree share `Z_v`. Multiplying the bounded degree share by the size martingale gives L¹ convergence to `W Z_v`, and the Markov property gives the closed martingale `e^{-t} M(T_t,v)`, with `M(T,v)=|T| E_T[Z_v]`. This argument does not require the projected local coordinates to be Markov. The elementary lower boundary follows from gain at rate at least `L` and promotions at rate `L/D≤1`: `(L−1)_+` has generator at least itself. The upper boundary follows from `r(v)≤D_v`. Localization at tree size is dominated on each fixed time horizon by the Yule size there.
- **Signed two-birth reduction.** The serial weight at `x` is `w_x=r_x/(D_x+1)`. Twin leaves have equal rooted laws and nonadjacent placements cancel in reverse order. On an old edge `xy`, attaching at `x` first reduces the next rate into `y` by `1/[D_x(D_x+1)]`, so the edge coefficient for the payoff at `x` minus the payoff at `y` is `s_y−s_x`, where `s_x=w_x/D_x`. Dropping the two nonnegative cross terms leaves a debit `s_x` per incident parallel tip at `x`. Thus the normalized bracket at `x` is its serial difference minus the average of its parallel tip payoffs. `Σ_x w_x≥1` follows from `Σ_x r_x=n` and `D_x+1≤n`. The two births contribute `n(n+1)` to the denominator and converting `M` to mean share contributes `n+2`.
- **All-tree local bounds.** For `v,w`, a core gains at rate in `[L,L+(D−L)/2]`, and leaf promotion occurs at `L/D`. For `U`, the pending leaf itself feeds its parent at rate one; the nonincreasing table makes this minimum rate the worst upper-drift case. Its large-parent boundary uses the predictable promotion hazard `≤1/d`, giving `E[e^{-τ}]≤1/(d+1)`, even with an atom at `τ=∞`. For `phi2`, adjacent cores remain cores and their leaf sets are disjoint, so no single birth changes both tracked pairs. Each contributes exactly `1/D` toward the other's gain, in addition to the live-leaf and external-core rates. For `phi1` and `phi1_x`, distinguished-tip promotion removes exactly one leaf from the middle and hands off to the adjacent-core row with the tip at `(2,1)`. Tracking the parent `x` adds a separate degree increment, with rate in `[1/D_a,D_x]`; the upper endpoint includes leaf callers, even when `x` was an old leaf. At parent exit use the untracked `phi1`; at middle exit use the one-vertex lower and pending-leaf upper. The generator drift is affine in the uncertain rates, so the box corners suffice even if the true rates are correlated. At an infinite stopping time the relevant degrees are bounded and the terminal share is zero. Bounded discounted local values and the closed martingale justify stopping, including for signed differences.
- **Founders.** The edge seed is invariant under label swap. At size three, founder 2 is the center with probability `1/2`; otherwise it and newborn 3 are the symmetric leaves and the difference is zero. Thus the mean gap is exactly `(M(center)−M(leaf))/6`: the additional factor `1/3` comes from the common size. The two-birth theorem then applies starting at `n=2`, including the comparison of labels 3 and 4. For arbitrary seeds the conditional claim means starting from that tree, or the corresponding transition kernel given a reached tree.

## Exact margins and independent attacks

The two completions at an old vertex of degree `D` give a parent of degree `D+1` in both cases. The serial middle has `(D_a,L_a)=(2,1)` with a distinct core parent even when `D=1`. The exact pointwise margin against **each** parallel tip is bounded below as follows:

| Old degree | Serial difference lower | Parallel tip upper | Margin lower |
|---|---:|---:|---:|
| 1 | `95490419/200000000` | `367808441/1000000000` | `54821827/500000000` |
| 2 | `49433051/125000000` | `70025291/250000000` | `28840811/250000000` |
| ≥3 | `834963/3125000` | `14205589/62500000` | `2493671/62500000` |

The last row uses the generic serial bound and monotonicity of the clipped `U` bound. Its margin is smallest and exceeds `1/26` by `1167723/812500000`. The center lower `v(2,2)=1187553139/1000000000` minus `U(2)`, divided by six, gives `409872349/3000000000`, above `2/15` by `9872349/3000000000`.

The independent replay decodes Prüfer words and uses literal caller-to-neighbour marks. It checks 3,160 mark pairs on 25 trees, comparing **coefficient vectors** indexed by rooted final tree, rather than a scalar test score. It also verifies the dropped cross-term vector and the founder row. For the physical drift it reconstructs the graph after each literal arrow and recomputes the tagged coordinates, using the pinned tables only as values. The resulting 581 concrete rows cover all six tables and 11 exit/promotion categories. These rows are regressions; the all-tree argument uses the certified corner inequalities. In addition to the 10 source controls, the independent replay rejects 10/10 audit plants: seven one-tick cutoff defects, one bool tick, a reversed signed edge coefficient, and an inflated physical promotion payoff.

## Transfer to sharp tails

The two-birth coefficient identity and the **parent-degree-tracked, promotion-aware** local recursion transfer directly to a strict tail payoff. For the tail, keep the exact signed edge coefficient where helpful; the mean proof's endpoint split discards useful cancellation. The state must also include the current size `s`. Replace the mean's discounted CT drift by a size-indexed one-step inequality: for a lower table at coordinate `c`, require `s F_s(c) ≤ Σ_e r_e F_{s+1}(c_e)`, including size-advancing births that leave `c` unchanged; reverse the sign for an upper table. The physical rates and the promotion handoffs are the same. A finite size horizon makes this an ordinary backward comparison, with independent tail lower/upper envelopes there (e.g. the Poisson live-leaf minorant (`open/age-order/hub-ranking/poisson-tail/README.md`) and the locked beta majorant). The mean tables and their degree-cutoff payoffs do not give those envelopes.

There is already a concrete one-step improvement to the serial-middle Poisson lower payoff from this transfer. At completion size `m=n+2`, let `d=D_T(x)+1` be the old parent's new degree, `z∈[0,1)` the threshold, and `q=1−z`. The middle has one tracked live leaf, gains at exact rate `1+1/d`, and loses its tracked leaf at exact rate `1/2`. The baseline Poisson row uses gain rate one and loss rate one. Applying the lower envelope at size `m+1` to the true first step therefore gives

```math
 \Pr(Z_{\rm middle}>z\mid T_{x-a-b})\ge
 H_{m,1}(z)+\frac{H_{m+1,2}(z)-H_{m+1,1}(z)}{md}
             +\frac{H_{m+1,1}(z)}{2m}
 =e^{-1}q^{m-1}\left[1+\frac{z}{d}+\frac{q}{m}\left(\frac1d+\frac12\right)\right].
```

Here `H_{m+1,2}−H_{m+1,1}=e⁻¹(q^m+m z q^{m−1})`, by cancellation of the one-stock beta component. The replay independently evaluates the first-step expectation from literal arrows and binomial tails in 745 rational rows. This is a conditional strict-tail credit at the correct power of `q` and uses no moment inequality. It does not settle the signed tail gap at thresholds where the adverse tips are still too expensive. Promotion-aware upper rows and retained edge cancellation are the corresponding places to look for a sharper debit. No extended tail band or stochastic ordering is claimed by this audit.
