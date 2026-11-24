> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# Conditional threshold order down to one tenth

Write `b` for the older label of the comparison `b→b+1`, and condition on any old tree of size `b−1`. **For every `b≥21` and every `1/10≤x<1`, the conditional tail gap is strictly positive.** In particular this proves the quarter threshold on every such tree. For `21≤b≤60`, the interval ledger below covers `[1/10,1/3]`; for `b≥61` one analytic corner covers `[1/10,1)`. The [Poisson birth-scale theorem](../../hub-ranking/poisson-tail/README.md) covers `x≥1/3` already for `b≥21` (`bx/(1−x)≥21/2>6`), closing the upper seam independently of any one-third receipt.

For the **unconditional** all-label conclusion on `[1/10,1)`, combine this with [Fable's two-program receipts](../../../../attacks/fable5/README.md) for older labels `2,…,19` and its **single-program** exact receipt for `20→21`; the same audit gives strictness for those receipts for `0<x<1`. Labels 1 and 2 have equal laws. This certificate does not claim a conditional `20→21` result on the whole band.

## Threshold kernels

At size `n` divide each sharp tail by `q^(n−1)` where `q=1−x`; this is its **own** normalization at the actual threshold. The [Poisson stock lower](../../hub-ranking/poisson-tail/README.md) with one tracked leaf gives a normalized serial middle lower `M=36787/100000<1/e`, uniformly in `n,x` and the old parent degree. For the signed placement numerator, common completion has size `b+1`, and the positive factor back to the probability gap is `q^b/[b(b−1)]`. The replay therefore compares sharp tails directly.

For adverse pending tips use the [actual-stock Bellman bound](../annealed/extension/third-physical/README.md). On a closed interval `[x_lo,x_hi]`, let `q_lo=1−x_hi`, `q_hi=1−x_lo`. Its two locked beta caps have first shape `r=1/k` or `1+r`, respectively. In the normalized series of that proof the pending cap is increasing in `q` (factor `q^(1−r)` and nonnegative series), and the promoted cap is decreasing (factor `q^(−r)` and positive function represented by the series with nonpositive nonconstant coefficients). Evaluate the pending cap at `q_hi`, with 48 terms and its positive geometric remainder, and the promoted cap at `q_lo`, with 48 terms and no remainder. Gamma, root and final grid rounding are outward as in the source. At each backward step multiply the two candidate rate sums by `q_hi/n`: all entries and rates are nonnegative. Minimize with the independently valid pending beta cap. Horizon is 14 sizes beyond completion. This safely discards any correlation of the threshold dependencies within a row.

For favourable pending tips, use a separate 24-step lower. Upon promotion pay the same `M` because the mark has just acquired a fresh leaf. Until then the promotion rate is `1/k`, the parent growth rate lies in `[0,k]`, and other births leave `k` fixed. At every backward step use `q_lo/n` and the smaller of the two endpoint rate sums; use zero at the finite horizon and round down. Sibling losses need no separate state for this favourable lower. These rates are disjoint except that all events leaving `k` unchanged and not promoting the mark are grouped as holds.

## Signed ledger

Use the authenticated third-physical local ledger first: it tests sitewise reserve or suppresses bare degree-two chains. If it fails, improve each oriented edge with the favourable payout from the *opposite* side before taking its positive part, keeping the coefficient difference signed. In the improvement, a core with degree `d`, `c` core neighbours and receive `r` uses `s=(d−c+r)/[d(d+1)]`. On an old-leaf edge the debit slopes are `(U(2,0)−L(d+1))_+` and `(U(d+1,d−c−1)−L(2))_+`. Use the same census credit and baseline `.004` as the old ledger. A core of degree at least eight pays `.002` per direct port; a direct neighbour of degree at most three takes it. Unused transfer is discarded, and suppressed ports are excluded from both. The chain bills are `.004/.008` or `.005/.009` (second is for a terminal degree-two core).

Coarse site checks use receive extremes and discard favourable payout. For a failing small site keep the multiset of direct degrees `2,…,7` and one high bin `8,…,b−2`. The high bin has reciprocal endpoints `1/(b−2),1/8`, credits the matching transfer where eligible, and keeps a minimum opposite coefficient `1/[d(b−2)(b−1)]`; its favourable payout may be dropped. For a known neighbour of degree `e`, its coefficient is at least `1/[de(e+1)]` and its favourable payout at least `L(e+1)`. At each multiset the balance is concave in total receive, so endpoints suffice. Old stars are tested separately. For chains, test lengths 1–3 and both endpoint receives using outward integer bounds; longer chains add nonnegative credit at coefficient `1/6` and no new internal debit. To skip an endpoint pair with degrees `d,e`, note that a length-`h` chain and both endpoints use at least `h+2+(d−1)+(e−1)=h+d+e` distinct old vertices by acyclicity. Baseline and census accounting leave at least `.012` in a multi-core suppressed ledger. No old-tree enumeration is used. If a threshold box fails both tariffs and the original ledger, bisect and require both closed halves to pass.

The analytic `b≥61` corner imports the [paired degree-two/degree-three beta debit](../annealed/extension/first-sender/joint/README.md) and checks the same strict rational inequality, with `M=36787/100000` and a sixth-root upper `727/1000` at `(61,1/10)`. The normalized debit decreases in both `b` and `x`; no monotonicity of the physical signed gap is assumed.

## Replay

From the Astra root, run `python3 -IB round2/open/age-order/global-transport/quarter/check.py`. For parallel replay or time budgets use the disjoint label shards below (the stop is exclusive; the default analytic corner clips the final shard at 61):

```sh
R=round2/open/age-order/global-transport/quarter/check.py
python3 -IB "$R" --start 21 --stop 32
python3 -IB "$R" --start 32 --stop 38
python3 -IB "$R" --start 38 --stop 40
python3 -IB "$R" --start 40 --stop 47
python3 -IB "$R" --start 47 --stop 51
python3 -IB "$R" --start 51
```

The checker authenticates its three imported programs via Git blob hashes. The separately audited Poisson input proof and checker have blobs `1a436ac3a8703327dea6f9949dad334905b3c7af` and `0452a2890979bc03f3f59fb80a2fcd0d83079622`; replay the latter with `python3 -IB round2/open/age-order/hub-ranking/poisson-tail/check.py`. The fixed quarter control is a single feasible tree (`b=21`, nine two-edge spokes and one ordinary leaf). Even exact coefficients with the new opposite payouts give a negative *certificate score* at the former `.211` serial coefficient; the same score is strictly positive with the Poisson coefficient. It is not a physical counterexample.
