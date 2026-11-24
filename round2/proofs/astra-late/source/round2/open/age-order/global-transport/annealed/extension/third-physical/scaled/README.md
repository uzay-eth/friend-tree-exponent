> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../../README.md#claim-26) for exact scope and citations.

# A uniform birth-scale cutoff of 5 with the Poisson-stock lemma

Using the new Poisson-stock serial lower bound, the conditional comparison holds on **every old tree for `b≥20` and `0<x<1` with `(b+1)x≥5`**. The finite physical boxes below the symbolic seam prove the same conditional order for **all `b≥20` on `[1/8,1)`**. Fable's paired receipts then give the strict annealed order on this interval for all adjacent pairs starting with `2→3`, with the founder tie. In particular, the quarter threshold needs no separate finite replay: `21/4>5` already at `b=20`.

The input is the Poisson-stock subharmonic bound from the separate `hub-ranking/poisson-tail/` proof: at serial middle stock one and common size `m=b+1`, `H(m,1;x)=e^(-1)B(m,1;x)=e^(-1)q^b`. The exact ledger uses `M=36787/100000`, verified below `e^(-1)` by its odd Taylor partial sum. The symbolic run changes only `M` and the entrance cutoff in the ledger described below. The large-degree bound uses `5^(1/32)<3/2` (as does the original run with 13). The finite run uses the pinned physical box checker unchanged, covers `b=20,…,38` from `1/8` to `5/(b+1)`, and uses rational cuts spaced at most `1/200` apart. At `b=39` the symbolic seam equals `1/8`.

At `x=1/10` the symbolic certificate starts at `b=49`. For `b=20`, the frozen 14-step continuation and fixed fees fail the degree-two terminal check even on the single point `1/10`: at `d=2,c=1,r=0` its credit is `(M−U)·7/12−U/12`, below the required `.008` plus the direct-edge excess. This is a certificate limitation; the remaining all-label tenth-threshold seam is not closed here.

### Poisson-bound replay

```
nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/scaled/poisson.py
PASS symbolic ledger at y>=5

nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/scaled/poisson.py --finite
# PASS for each b=20,…,38 from 1/8 through 5/(b+1)

nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/scaled/poisson.py --obstruction
LIMIT b=20, x=1/10: degree-two terminal misses reserve with U(2,0)=286974353/1000000000
```

## Baseline certificate without the Poisson-stock lemma

On every old tree, for **all older labels `b≥21`** and **all `0<x<1` with `(b+1)x≥13`**, the older of the two next births has a strictly larger conditional tail. This is a symbolic size range, not a sweep through a finite label cap. The accompanying finite-label replay extends the all-label interval to **`[799/2400,1)`** using only Fable's paired receipts. This does not prove the whole range `x≥1/10`: the symbolic result covers that threshold for `b≥129`, but does not cover the finite sizes below that seam.

## 1. Uniform entrance price

For a pending leaf at common size `m=b+1`, with parent degree `k`, put `r=1/k`, `q=1−x` and `y=mx`. The frozen locked beta bound, normalized by `q^(m−1)`, is no larger than

```
 Γ(m)/(Γ(m+1−r) Γ(r)) · (q/x)^(1−r)
 ≤ [qm/(m−r)]^(1−r) y^(r−1)/Γ(r).
```

For the first inequality, bound the beta-integral factor `z^(r−1)` by `x^(r−1)` over `z≥x`. For the second, use `Γ(m)/Γ(m−r)≤(m−r)^r` by Jensen on a gamma variable and the Gamma recurrence. The factor in square brackets is at most one whenever `y≥r`. Consequently for any `Y≥1` the normalized price is at most `Y^(r−1)/Γ(r)` throughout `y≥Y`, for **every** allowed `m`. No stationarity or limiting approximation in `m` is needed. The replay builds constant upward prices for `Y=13`. For `k<32` it uses lower rational roots and a lower bound for `Γ(1/k)`: integrate the odd alternating Taylor polynomial on `[0,1]`, then use downward exponential values and right-end rectangles for the decreasing density on `[1,8]`. For `k≥32`, `Γ(1/k)≥k−1` by integrating `1−t` on `[0,1]`, and `13^(1/k)<3/2`, so `2/[(k−1)13]` suffices.

### A physical postpromotion price, and where its use stops here

At a promotion taking the size to `N`, condition on the actual degree `k` of the permanent parent. The mark has one fresh leaf and the locked shape is `1+r`, `r=1/k`. For any threshold, the tail conditional on this history is at most

```
 (N−r)^r q^(N−1−r) / Γ(1+r)
      · [ x^r + r x^(r−1) q/(N−r) ].
```

Indeed, in the beta integral bound `z^r≤x^r+r x^(r−1)(z−x)` and integrate; the mean overshoot under the remaining density is `q/(N−r)`. The Gamma ratio is bounded by Jensen on shape `N−r`. The `k` here is the degree **at the actual promotion**, not the entrance degree. A pending stock state `(n,k,j)` has promotion weight `1/k`, sibling loss weight `j/k`, and gain weight in `[j+1,(k+j+1)/2]`; successor stocks are exactly as in the [physical Bellman proof](../README.md). Thus the bound can be inserted after promotion at the actual `N` without averaging stock independently of the size or the degree. At a symbolic finite horizon one can always take the minimum with the uniform entrance bound at the terminal state. The signed symbolic cutoff below uses the entrance bound itself; it does **not** claim a further numerical improvement from this new analytic promotion bound. The separate finite-label strip does use the existing stock-sensitive physical continuation with its sharper integral upper.

## 2. Unbounded degree ledger

Use the [reviewed serial lower](../README.md) `.211 q^b` for `b≥21`. For each normalized adverse price substitute the uniform entrance cap above at the appropriate post-placement parent degree. Writing `U=cap(2,13)`, `g=.211−U`, use the frozen signed credit `C` on the core and its `l` leaves; use `cap(d+1,13)` at core parents (this symbolic cap ignores their stock). Retain the frozen skeleton construction but take `a=.006`, baseline `.004`, and chain fees `.004` or `.008` at degree-two terminals. For each retained degree below 32 check the two receive extremes for **each integer chain-port count**; for chains of lengths 1–3 check the paired boundaries and internal debit at the receive corners. Both reductions are as in the frozen ledger. For high endpoints use reciprocal `[0,1/32]` and upper own direction `z/33`, where `z=2/(32·13)`. The chain-side direction is bounded simply by its full coefficient times `cap(3,13)`.

There is no degree cap in the theorem. At `d≥32`, all core-side costs, including outgoing direct edges and the core direction on leaf edges, cost at most `wz` together, and the direct excess vanishes since `sz≤z/33<.004`. After paying the core-side cost, for `l≤d/2` the remaining stock credit is bounded below by

```
(g−z)w + (g−U)l/(2d) ≥ (l/d) [(32/33)(g−z)+(g−U)/2].
```

For this range one may need a negative stock floor: take `min(0, [(32/33)(g−z)+(g−U)/2]/2)`. The replay verifies that even with that floor the available `(a−fee)d/2` pays `2a+.004`. For `l>d/2`, since `t−s≤1/[2d(d+1)]` whenever positive, the total adverse leaf direction is at most `U/[2(d+1)]`. The strict stock credit `(32/33)(g−z)/2+g/4` pays that toll and the same baseline and core-degree costs without needing the port allowance. High stars can instead use their centre credit after the core-side cost alone (`l=d` makes `t−s≤0`). These tests extend the ledger to every physical degree. Summing the skeleton gives positive normalized numerator at least `2a+.004=.016` in the nonstar case. Small stars sign in the explicit local grid.

## 3. Finite splice and remaining range

The new finite check imports the [band checker](../band/README.md) and adds a single **closed** threshold box `[799/2400,1/3]` for older labels `20,…,38`. It uses `.2124` at label 20 and `.211` elsewhere, exactly as before. For `b≥39` and `x≥799/2400`, `(b+1)x≥40·799/2400>13`, so the symbolic theorem is available. The frozen continuous band supplies the rest of `[1/3,1)` at the finite labels. With paired early receipts this proves the stated all-label interval, with founder equality. The numerical gain in the common threshold is small; the substantive gain is the uniform birth-scale cutoff.

At `Y=12`, the same unmodified symbolic ledger fails its degree-two terminal check. This is a limitation of the current normalized caps and fees, not a negative physical gap. Extending to `x=1/10` would need a lower-scale conditional argument for a finite range (in particular `b=20,…,128` after the early receipts), or a substantially lower symbolic cutoff. No finite-label claim on that range is inferred from the analytic promotion price.

## Replay

```
nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/scaled/check.py
PASS unbounded degree/label ledger at y=(b+1)x>=13; 735 endpoint types

nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/scaled/check.py --finite
PASS finite older labels [20,39) on [799/2400,1/3]
```
