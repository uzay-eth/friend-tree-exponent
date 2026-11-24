# Signed zero-mark renewal: finite signed solver and an all-degree TV tail

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 14: Claude for transfer/scalar identities; GPT for full contraction and sign transfer; inputs include Claim 12 and the cap-48 continuation certificate. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


At a fixed discount `s` write `C h(S)=p_o Σ_Y h(Y)`, with `Y` ranging over immediate core children. The signed residual equation is `(s−L⁰)U=b+C ΔU`. Thus `H=ΔV_s b+T_s H`, `T_s=ΔV_s C`. This directory does not use a positive-gauge weight or infer a root of the full susceptibility from a killed value.

The later [six-count pair certificate](pair/README.md) bounds the unlimited norm uniformly by `0.860052` for `s ≥ 565658/10^6`, including the low-degree sector left open by the earlier tail below. The [physical sign theorem](critical.md) combines this with the permanent-mark trigger. See the [local audit](pair-audit/README.md) and cold referee report for their scopes.

## Exact finite controls

Kill the intrinsic Green evolution at its first jump above total size `N` and take both `U(S)` and `U(S⁺)` in this same domain. A child needed by `C` always has `size(Y)+1≤size(S)≤N`, so this differentiation never requires an outside value on a queried child. A descending-size Green sweep computes the entire rational child-occupation measure and reward on each state. Subtract the two exact measures before taking absolute values. For a proposed `H=x`, if `w−|T|w>0`, outward error in coordinate `i` is at most `w_i max_j |x−T x−f|_j/(w−|T|w)_j`, with `f=ΔV_s b`. All arithmetic for this final check is rational; only the proposal uses floating point. For `U(J)` use the positive occupation row with these signed intervals. Independently enumerate individual directed arrows and check the resulting score against the full signed generator, including both `+p` and `−p` child entries. A hand formula at size 3 verifies the loop through the bare parent. Deleting the plain negative copy gives an independent full residual of `0.414082642416` at size 8 and is rejected.

At `s=.565728`:

| N | Intrinsic shapes | Difference rows | Max TV (outward) | Shaped TV (outward) | Certified killed `U(J)` |
|--:|--:|--:|--:|--:|:--|
| 8 | 199 | 84 | .196865546961 | .173493738979 | [.538120420139,.538120420140] |
| 11 | 3,046 | 1,204 | .219457022399 | .195064771531 | [.595207568706,.595207568707] |
| 12 | 7,812 | 3,046 | .226994907188 | — | [.609998610225,.609998610226] |

In particular these exact row sums rule out a radius-`≥1` counterexample on the tested domains. They are not estimates of the unlimited radius. Size 13 did not finish within the short slot. The killed signed equation has no positivity comparison in `N`, even when each finite `T` is contractive. One must certify the face of the unlimited operator before using `U(J)=1`.

### Rank-one split

With `c` the actual number of core children, `p c=1−b−p`. Separate `T=A+u e_J^*` with `u=T1` and `Ah=T(h−h(J))`; in the killed calculation compute `u` by the row sum (the term from intrinsic lifetime is not `Δ(1/s)` at a killed boundary). The exact denominator of the rank-one update is `1−[(I−A)^{-1}u](J)`. At size 8 the denominator is in `[1.068280718340,1.068280718341]`; the residual replay gives the same `U(J)` enclosure. Centering at `J` increases the raw max TV from `.196865546961` to `.358453046768`, so it is not a TV improvement in this norm.

## Unlimited high-degree rows via common-child matching

The given global infection count `3e^t−t−2` is too large to integrate directly against `e^{-st}` at this `s`, and degree discrepancy alone does not control a TV toll that can tell arbitrary whole child shapes apart. Instead stop on the **first unmatched root selection**. Start with one private leaf at the root in A. Until the stop, couple the common children *as entire subtrees* and pair all common leaves. Independently drive leaf-to-root arrows; common ones add a paired leaf, while a private leaf adds another private leaf. Other sources' degrees are equal on both sides, so their sends to the root can be paired exactly. Split the root's individual nearest-neighbor channels: rate `1/(D+δ)` for each synchronized selection; its surplus on each common neighbor goes to B. Here `D` is B's current degree, and `δ` is the number of extra private leaves in A. An A-only selection of a private leaf is also a stop (it would promote it). The root's inert slot sends no intrinsic selection. Both the synchronized promotions and sends preserve the invariant. Before stopping, `δ` is the rate-one Yule population from 1, stopped on a competing hazard (its expected value before the stop is at most `e^t`). In particular no descendant can differ before this stop. This only uses locality up to a stopping time, not independence conditional on a future realized history.

Write `d=D(0)`, and `τ` for the first unmatched selection. Its predictable pre-stop rate is

```
(D−1)[1/D−1/(D+δ)] + δ/(D+δ) ≤ 2δ/d.
```

The first term counts all common neighbors except the inert slot; omitting the private-leaf term is a planted error. It follows that `P(τ≤t)≤2(e^t−1)/d`. For a bounded signed child toll with `||h||∞≤1`, before the stop the entire signed occupation cancels except the reciprocal discrepancy on common children, of magnitude at most `δ/d`. After the stop both total rates are bounded by 1, so their difference is bounded by 2. Pointwise in time the resulting norm is at most `min(2,5e^t/d)`. This is a bound on the TV of the true, **unlimited** signed occupation measure after integrating; it is not an estimate of a coupled difference for a single selected toll only. Integrating at `t*=log(2d/5)` gives for `d≥3` and `0<s<1`:

```
||T_s(S,·)||TV ≤ 2 (5/(2d))^s/[s(1−s)] − 5/[d(1−s)].
```

The original integral is nonincreasing in `d`, so checking its value at a degree threshold is uniform over larger degrees and over all child forests. At `.565728`, outward bounds at thresholds 32, 64, 128, 256, 1024 are `1.568112915`, `1.123282746`, `.790943220`, `.550469436`, `.260824113`. The rational replay rounds the exponent down only in the positive term. This resolves the high-degree sector in the unweighted norm but does not bound the low-degree, arbitrarily complex child forests by a constant below one. Conditioning on a first root discrepancy only gives the crude bound 2 thereafter; a low-degree certificate needs a finer return/occupation estimate in that sector or a shaped power bound. The high-degree bound alone supplies no full root sign.

```sh
taskset -c 8 timeout 12s python3 -Bu round2/influence/zero/solver.py --cutoff 8 --split
taskset -c 8 timeout 12s python3 -Bu round2/influence/zero/solver.py --cutoff 11
taskset -c 8 timeout 12s python3 -Bu round2/influence/zero/solver.py --cutoff 12 --unweighted
taskset -c 8 timeout 12s python3 -B round2/influence/zero/independent.py 3
taskset -c 8 timeout 12s python3 -B round2/influence/zero/independent.py 8
taskset -c 8 timeout 12s python3 -B round2/influence/zero/independent.py 8 --plant
taskset -c 8 timeout 12s python3 -B round2/influence/zero/obstruction.py
```
