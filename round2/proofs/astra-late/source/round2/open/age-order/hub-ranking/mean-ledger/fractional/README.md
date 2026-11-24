> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../../README.md#claim-27) for exact scope and citations.

# Stock-sensitive fractional moment ledger

**Every finite old tree has a strictly positive conditional adjacent moment gap for every real `0<p≤473/512`.** If the old tree has at least two cores, the normalized numerator is at least `.033`. The edge seed and stars sign separately. The [fractional-middle proof](../../fractional-middle/README.md) supplies `[473/512,1]`; no tiny-power perturbation is required to reach that seam.

Use `E W^p/Γ(1+p)` at common completion size `m=n+2`, as in the [real ledger](../real/README.md). The gap in degree-share moments is the normalized numerator times `Γ(1+p)/[n(n+1)(m)_p]`. The normalization uses the true size amplitude independent of the jump chain from a fixed completed tree, not after a return or promotion race.

## Keep the actual adverse leaf stock

Use `M=middle(lo,hi)` from the fractional-middle checker, with its harmonic stock boundary and downward recurrence. Use its rational Euler-product bound for the two-race postpromotion cap from [tip moments](../../tip-moments/README.md) at the true parent degree `k`. Write `B(k)` for the upward cap, also bounded by `1+1/k` on `0≤p≤1`. Indeed the unrefined `A_p(1+1/k)` increases to `A_1(1+1/k)=1+1/k` by the Euler product.

While the marked tip is a leaf, initially tag `j` of its other sibling leaves. The monitored stock gains have rate at least `j+1`, including the marked leaf's own arrow into its parent. Its parent's true degree `k` rises by one at any such gain. Tracked sibling promotions decrease stock at rate `j/k`; promotion of the marked leaf has rate `1/k`. On an extra parent gain, track the fresh leaf too, so it has exactly the same successor as a stock gain. At cutoff degree 48 in the first 36 cells and 32 thereafter, the remaining upper is

```
(1+1/k) min(1/(j+1), 1/(1+lo*k)).
```

The first factor inside the minimum is the physical all-time first-send barrier, retaining a permanent nonleaf slot; the second is the Laplace bound from promotion hazard at most `1/k` thereafter. In both cases the reward conditional on promotion is bounded by `1+1/k` since the parent degree cannot decrease. No independence or fractional interpolation of hub and mean probabilities is needed.

Below cutoff the discounted supersolution is obtained by unwinding

```
[B(k) + k(j+1) U(k+1,j+1) + j U(k,j-1)]/[1+k(j+1)+j+lo*k]
```

with upward rounding, increasing the value further if necessary to at least `U(k+1,j+1)`. The replay checks gain monotonicity and stock monotonicity on every finite row. The first makes extra physical gains nonpositive in the stopped drift; using `lo` rather than `p` weakens the discount safely. The marked tip on a core edge whose old parent has total degree `d` and leaf stock `l` has shifted type `(d+1,l)`; for the core side of an old-leaf edge it has `(d+1,l−1)` because that old leaf became a core. The adverse serial tip has `(2,0)`. Favourable parallel tips use the true discounted return lower at `hi`, truncated at parent degree 32, with postpromotion payoff `M`. The return argument is as in the integer ledger and does not require the adverse tip's stock.

## Low-exponent host credit and simultaneous payment

Below `1/4`, retain the old parent of the serial middle during one preliminary race. At its post-placement degree `k`, the middle has two incoming gain sources of rates `1` and `1/k`, the loss of its fresh leaf has rate `1/2`, and growth of the old parent has some predictable rate in `[0,k−1/2]`. The channels are disjoint; a middle-to-old-parent send is included in parent growth. On gain pay the harmonic-recursion lower `V(3,2)`; on loss pay the actual return-state lower `V(2,0)`. On parent growth continue the same race at `k+1`; at cutoff 64 use the unboosted `M`. Writing `x=(1+1/k)V(3,2)+V(2,0)/2` and `y=hi+3/2+1/k`, take the smaller of `x/y` and `(x+(k−1/2)H(k+1))/(y+k−1/2)` and round down. These two endpoint tests are sufficient because the discounted drift at the lower test is affine in the unknown rate. Taking the maximum with `M` is safe: the baseline is a separately available physical lower bound at every state of this race. The lower is `H(d+1)` at the serial host; on old leaves use `H(2)`, and on every suppressed site use `H(3)`. The favourable parallel-tip lower still uses `M` only.

Use `g=M−U(2,0)`, the old-leaf-adjusted credit from the [integer ledger](../README.md) with actual shifted uppers and any host bonus, and suppress bare degree-two chains. A retained endpoint pays baseline `b` and receives signed census credit `a(c−2)`. Use `a=b=.015`, per-chain bills `.007` and `.022` (the latter for degree-two terminals) in the first 36 cells. All later cells use `a=b=.011` and bills `.007` and `.020`. Pay excess above `b` only for direct ports. In the first 36 cells, every endpoint of degree at least 15 budgets `.005` on every direct port, and each direct neighbour of degree at most 3 actually receives it. Other budgets may be discarded. The transfer thus never increases the summed credit.

For the finite coarse check bound each directed debit by `s U(d+1,l)`. At each fixed number of chain ports `k`, use concavity in receive `r∈[k/2,c/2]`. If the maximum directed debit is at most `b`, the debit excess vanishes and only extreme port counts are needed. When a coarse corner fails, keep exact direct neighbour degrees through 14 for centres of degree at most four, or through six otherwise. Use their guaranteed opposite coefficient `1/[d e(e+1)]` and favourable lower `L(e+1)`. The remaining slots use a reciprocal interval up to `1/(cut+1)` with opposing coefficient and favourable lower set to zero. In particular an unresolved slot of a low transfer recipient always has degree at least 15. At a fixed multiset the resulting residual is concave in total receive.

For `d≥D` (`D=96` in the first 36 cells, 64 thereafter), write `x=l/d`, `T=(1+1/(D+1))/(D+1)` and `z=(D/(D+1))g+(g−U(2,0))/2`. All direct excess has vanished by `s U(d+1,l)≤T/2`. The core-side old-leaf debit is at most `T` in total by `U(d+1,l−1)≤(1+1/(d+1))/l`. If `x≤1/2`, credit before core allowances is at least `min(z,0)/2−T`. The remaining allowance per port is at least `a−.007−transfer`, and there are at least `D/2` ports. If `x>1/2`, the other old-leaf debit is at most `U(2,0)/(4(D+1))`, since `w≥l/(d+1)` and `(1/2−l/(d+1))_+≤1/[2(d+1)]`. Its serial stock credit alone is at least `(D/(D+1))g/2+g/4`; subtract that debit and `T`. Both inequalities pay `2a+b` in all cells.

The chain tests keep actual stock and the paired difference at low endpoints; from degree 10 through the degree just below the tip cutoff, keep exact reciprocals and the worst outgoing exposure across stocks. Beyond that use reciprocal endpoints `0,1/cutoff`, `L=0`, and the better of the stock and discount exposure ceilings at the first degree of the tail. Convexity in the ranged coordinates justifies corner checks. Test lengths one, two and three after both bills. Extra interior sites are identical with coefficient `1/6`, so add positive serial credit and no additional debit. Summing on the retained tree gives `b+2a≥.033`. For a large star the centre-side leaf debit is bounded by `T`; centre credit alone exceeds it. Small stars and the seed are checked separately.

## Replay

`python3 -IB round2/open/age-order/hub-ranking/mean-ledger/fractional/check.py` checks 473 cells of width `1/512` (allow several minutes). The first cell uses the zero endpoint in numerical enclosures for positive `p`, not a claim about `Z^0`. To shard the replay, use `--start` and `--stop` with half-open integer bounds. The imported lower checker is `fractional-middle/check.py` at blob `d8111c15376ad27f2fe0a583cabf8192f925ec88`, with proof blob `5921ba52b80a3cfc8ccba1b16b6a77f66d7cfcfa`.
