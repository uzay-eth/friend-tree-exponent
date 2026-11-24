> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# Hubhood: recycling a lost leaf in the signed pre-placement score

At `x=0`, every finite **length-two spider** with at least two arms has a strictly positive whole-tree pre-placement gap. The same certificate signs selected paths, comblike trees and length-three spiders. The general all-tree sign remains open; no negative physical score was found. Only degree/stock recurrences and the signed two-birth coefficients are used here. The [upper-threshold proof](../README.md) is separate.

## 1. A tracked leaf is not necessarily lost forever

Write `B(d,l)` for a lower bound on a marked core's hub probability with physical degree exactly `d` and `l` tracked adjacent leaves, `0≤l≤d−1`. The restart at stock zero below is only called after a last **actual** tracked leaf is promoted by the mark, with its new degree exactly two. It is not an assertion for general bare cores.

At that promotion let `k` be the degree of the newly promoted neighbour. Until the mark's next gain its degree is `d`, so the new neighbour sends a leaf back at rate `1/k`. The neighbour's gain rate is at most `k−1+1/d`: all other `k−1` neighbours contribute at most one each. A gain elsewhere at the mark counts as success just like the specified return, since it too supplies a fresh leaf. A lower probability of getting that success by the time `k` reaches 32 is `R_d(2)`, with

```
R_d(32)=0
R_d(k)=[1+k(k−1+1/d) R_d(k+1)]/[1+k(k−1+1/d)].
```

This is a bounded subsolution: `R_d(k)` decreases in `k`, so the maximal competing neighbour growth is pessimistic. At degree cutoff 32 failure is charged even if there would later be a return. Optional stopping works either via the race in each `k` or via the bounded subsolution; when neither exit has occurred the specified send rate remains at least `1/31`.

At positive stock use

```
B(d,l)=[d B(d+1,l+1)+B(d,l−1)]/(d+1)
B(d,0)=R_d(2) B(d+1,1)     (only at the last actual promotion)
B(20,l)=1−20^(−l).
```

Tracked leaf births have rate `l` and tracked promotions have rate `l/d`. All other gains at the mark give fresh leaves and are added to stock, so their drift is nonnegative by the verified inequality `B(d+1,l+1)≥B(d,l)`. These rates never require removing a fictitious tracked leaf; retaining the fresh degree-two neighbour at the **actual** last loss is essential. At cutoff the standard tracked birth–death comparison supplies the terminal payoff. The monitored transitions have almost-sure finite exits: each positive-stock stage has rate at least one until a tracked jump or degree gain, and a last-loss stage uses the stopped return race. There are only finitely many monitored stages before cutoff or failure, since each successful return increases the degree and stock losses at fixed degree decrease stock.

This raises the uniform serial middle lower bound from the old one-shot `0.613194634` to `0.777202011` (decimals rounded for display). The [frozen scalar first-send bound](../../preplacement-transfer/README.md#2-uniform-serial-gap-including-arbitrary-physical-exterior-input) gives the tip upper `P(2,0)=0.610259660`. Thus the uniform serial-gap lower bound rises from about `0.002935` to `0.166942351`.

Two additional signed credits are safe:

* For a leaf whose parent starts at degree `k≤32`, the same `R_1(k)` race gives a hub lower bound `R_1(k) B(2,1)`. Before it is promoted its degree stays one, so its return immigration into its parent is included in `k−1+1`. This buys credit for the **favourable parallel tip**. Above cutoff use zero.
* In a serial completion let the old neighbour of the middle have degree `p`, and put `a=1/p`. While `p` and the two-vertex middle/leaf configuration remain unchanged, the old neighbour sends to the middle at rate `a`; each such send earns `B(3,2)−B(2,1)`. The tracked gain and tracked loss rates are `1` and `1/2`, and the old neighbour's degree grows at rate at most `p−1+1/2` including any arrow from the middle. Use competing-rate upper `p+3/2+a`, conservatively. Stop on any of those events and retain the original continuation bound (the rescue bound after loss). The tracked gain and loss terms balance by the recursion; old-neighbour growth leaves the baseline payoff `B(2,1)` unchanged. Consequently the gain on top of the baseline is at least `a[B(3,2)−B(2,1)]/(p+3/2+a)`. Write `g_p` for the resulting serial-gap lower bound after subtracting `P(2,0)`.

## 2. Exact paired edge debit and an infinite family

Use the notation of the [parent signed formula](../README.md): `w_u=r_u/(D_u+1)`, `s_u=w_u/D_u`, `A=Σ w_u`. On an edge oriented so that `s_u>s_v`, the only adverse tag is the leaf at `u` in the parallel completion. Bound the edge debit by `(s_u−s_v)[U(u|v)−L(v|u)]_+`, using the frozen first-send upper and the return-race lower above. Add the signed serial credit `Σ w_u g_(D_u+1)`. This bound is on the numerator before division by `n(n+1)`.

For a spider with `k≥32` two-edge arms let `h=1/k`. At the intermediate arm vertex, `w=(1+h)/3` and `s=w/2`; at the old outer leaf both weights are `1/4`, and at the centre `s=1/[2(k+1)]`. There are two edge debits per arm. On the centre edge the adverse leaf has parent degree three and one other leaf; dropping both the centre's coefficient subtraction and its favourable leaf gives debit at most `(w/2)P(3,1)`. On the outer edge the adverse leaf has a degree-two parent with no other leaf; the favourable leaf has a degree-three parent. The remaining debit is at most `(1/4−w/2)[P(2,0)−R_1(3)B(2,1)]`. Drop the centre's positive serial credit. The lower bound per arm becomes

```
w*g_3 + g_2/4 − (w/2)*P(3,1)
    − (1/4−w/2)*(P(2,0)−R_1(3)*B(2,1)).
```

It is affine in `h`; the two exact endpoint tests at `0` and `1/32` prove it exceeds `0.006` per arm throughout. The finite range `2≤k<32` is checked by the same degree/stock bound, with no continuation-tree enumeration. This is a whole-tree sign, even on a sequence with an arbitrarily crowded bare centre. With the favourable tip discarded and without the initial old-neighbour credit the same simpler bound fails already on the 32-arm fixture. A failure of that lower bound is not a negative physical score.

## 3. Scope and replay

```
python3 -IB round2/open/age-order/global-transport/hostile/check.py
```

The rational replay checks the scalar comparison inequalities, the two endpoint tests and the small spider cases, five paths, four simple combs, and 36 crowded combs with 1–3-edge teeth. All 36 crowded combs sign; the weakest tested bound is `N_T/A > .042428843` (two-edge teeth, spine length eight, 16 teeth per spine vertex). These computations do not bound general trees. The remaining problem is to certify the signed local-profile sum on every tree, or find a profile where a **physical upper bound** on that sum is negative. The retry certificate works for hubhood and does not supply a bound uniform at positive birth-scale thresholds; strict positivity extends by right continuity only to a tree-dependent neighbourhood of zero.
