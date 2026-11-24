> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../README.md#claim-26) for exact scope and citations.

# Physical continuation closes the finite gap at one third

At **`x=1/3`**, for every older label **`21≤b≤39`** the older birth has strictly greater conditional tail probability, **on every old tree**. The old tree has size `b−1`, and both continuations start at size `b+1`. The [earlier analytic bound](../first-sender/joint/README.md) covers `b≥40`; Fable's [paired receipts and single `20→21` receipt](../../../../../../attacks/fable5/README.md) cover the earlier labels unconditionally. Thus the strict annealed order at `x=1/3` now covers every adjacent pair starting with `2→3`; the first two labels have equal laws. The new certificate is at a single threshold. It uses no implication from moment order.

## 1. Actual stock on the pending leaf

Write `q=2/3` and normalize a tail at current size `n` by dividing by `q^(n−1)`. While the tagged tip is a leaf, let `k` be its parent's actual degree and `j` the actual number of other leaf neighbours of that parent. The next jump, in rate-one-vertex time with total rate `n`, has the following **disjoint** effects:

| Event | Rate weight | Successor |
|---|---:|---|
| Tag promotion | `1/k` | restart at size `n+1` with the physical degree-two tag and old parent still degree `k` |
| Other leaf promotion | `j/k` | `(n+1,k,j−1)` |
| Parent gains a fresh leaf | `g` | `(n+1,k+1,j+1)` |
| Other birth | `n−(j+1)/k−g` | `(n+1,k,j)` |

Here `j+1≤g≤j+1+(k−1−j)/2`: the tag and each sibling are leaf senders of rate one, and the other `k−1−j` neighbours are cores. Keeping the exact siblings makes the displayed successor states valid; there is no stocking approximation and no independence assumption. In particular the parent's own send to a sibling changes that sibling's degree and creates a leaf **at the sibling**, not at the parent. Thanks to 42b5 for the tighter core-neighbour bound on `g`.

The locked degree urn at a pending state gives `Beta(1/k,n−1/k)` as an upper tail; at the actual promotion the upper is `Beta(1+1/k,n+1−1−1/k)`. Indeed the reciprocal contributions of the core parent and any subsequent other core neighbours can only fall, so locking the current offset and coupling an affine Pólya urn dominates the attachment fraction. The shift is `1+1/k` after promotion because the tag then has one fresh leaf as well as its old parent. This argument applies at each physical stopping state; in particular it does not equate clock-discounted moments to threshold tails.

Start the Bellman calculation at `b+15` with the pending beta cap. Work back to `n=b+1`, using the upward normalized promotion cap `P` at size `n+1`; write `G,L,H` for the upward normalized rows at gain, sibling loss and hold respectively. The next upper payment is

```
min(pending_beta_cap, (q/n) max_g[(P+j L)/k + g G + (n−(j+1)/k−g) H]).
```

The maximum is over the **two** endpoints of the physical interval for `g`. Every row is a conditional upper bound by finite backward induction; no monotonicity in `j`, `k` or `n` of a rounded array is required. The factor `q` converts between the two consecutive size normalizations. The strict stock reductions in the replay illustrate the benefit beyond the bare-state cap. At `b=21` even the bare caps are below `.158` for `(2,0)` and `.070` for `(3,0)`.

## 2. Rational beta normalization

For `a=r` or `1+r`, `r=1/k`, and size `n`, put `B=n−a`. Substitute `1-z=qv` in the survival integral. After dividing by `q^(n−1)` the factor is

```
[Γ(n)/(Γ(a)Γ(B+1))] q^(1−a) · Σ_(i≥0) c_i B/(B+i),
c_0=1,  c_i=c_(i−1) q(i−a)/i.
```

For `a=r` the prefactor Gamma ratio is `(n−1)! sin(πr)/(π∏_(i=1)^n(i−r))`. For `a=1+r` replace the denominator product by `r∏_(i=1)^(n−1)(i−r)`. These are the Gamma recurrence and reflection formulas. Bound `sin(πr)/π` above by `r[1−z²/6+z⁴/120]` at `z=3.14159r`; the alternating Taylor bound applies through `π/2`, and the bounding polynomial is decreasing there. Machin's formula and rational alternating arctangent bounds verify the stated lower value of π. The factor `q^(−r)` is rounded upward by integer root comparison. The pending series has nonnegative coefficients; after term 24 the remaining terms cost at most `c_24 B/(B+25) · q/(1−q)`. The promotion series has all negative nonconstant coefficients, so truncating at 24 gives an upper bound. The resulting beta values and the Bellman recursion are rounded upward on a `10^10` grid.

## 3. Local signed ledger

Take the uniform serial lower coefficient **`M=.211`** from a finer, outward rational replay of the already reviewed [one-family lower](../first-sender/joint/README.md). Its proof is unchanged; the calculation is `coefficient(21,768,56)` at frozen checker blob `d3b4ad5bc87760a843712d9af43af76450a07a37`. Only the **adverse** side of the comparison is replaced here.

Use `U(k,j)` for the Bellman row at common completion size. The serial tip has `(2,0)`. At an old core with degree `d`, `c` core neighbours, `l=d−c` leaves and core reciprocal receive `r`, let `w=(l+r)/(d+1)`, `s=w/d`, and let `t=1/(2d)` be the old leaf coefficient. Put `g=M−U(2,0)`. Local credit after all incident old-leaf edges is at least

```
C = g(w+l t) − l(t−s)_+ U(2,0) − l(s−t)_+ U(d+1,l−1).
```

Drop the last term at `l=0`. In parallel placement on an old-leaf edge, the old leaf becomes a core in the completion; this is why its parent's other stock becomes `l−1` on the core side. The old-leaf side has no siblings. At a core–core edge the endpoint cap is `U(d+1,l)` and an oriented adverse debit is at most `s U(d+1,l)`. For `25≤b≤39` the endpoint receive checks `r=0,c/2` give `C≥c s U(d+1,l)+(w+l t)/10000`. Concavity in receive closes the entire interval; summing retains `A/10000>0`, with `A=Σw` over all old vertices.

For `21≤b≤24`, suppress maximal chains of bare degree-two old cores. Each retained core reserves a baseline `.004`, receives `.004(c−2)`, and pays a bill on each of its `k` chain ports (`.004` except `.008` at a degree-two terminal). For direct ports pay only the excess `(s U(d+1,l)−.004)_+` because an edge has at most one adverse direction and its baseline is needed only once. Thus the checked retained credit is

```
C + .004(c−2) − .004 − k bill − (c−k)(s U(d+1,l)−.004)_+ ≥ 0.
```

Check **all** integer `0≤k≤c` at receive endpoints `k/2,c/2`. All degrees on old trees are at most `b−2`, so no asymptotic degree cap is used. Stars are checked separately by `C>0`. On a suppressed chain, a site has serial credit `2sg`, and the bare-chain adverse cap is `U(3,0)`. Internal edge costs are bounded by `U(3,0)|s−s'|`. At a boundary edge use the actual signed coefficient and the endpoint caps, charging the maximum of zero and the two directed debits. At a retained endpoint the adjacent chain port contributes exactly `1/2` to receive; the remaining receives range between zero and `(c−1)/2`. For fixed endpoint degree and stock the net chain cost is convex in the ranged receives, so the two endpoints on each side suffice. Lengths 1–3 suffice since additional sites have the same `s=1/6` and positive credit, with no extra internal cost. After both bills every chain fits the single `.004` baseline. The skeleton on retained cores is a tree with `Σ(c−2)=−2`; summing leaves a normalized numerator of at least `.012`. The factor converting it to a conditional probability gap is `(2/3)^b/[b(b−1)]`.

## Replay

```
nice -n 10 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/check.py --start 21 --stop 21
PASS pinned serial coefficient > .211 for b>=21 at x=1/3
PASS physical tip caps at b=21: (2,0) < .158, (3,0) < .070; strict stock gains

nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/check.py --ledger-only
PASS physical tip caps at b=21: (2,0) < .158, (3,0) < .070; strict stock gains
# PASS for each b=21..24 (chain reserve) and b=25..39 (sitewise reserve)
```

Local integer profiles and endpoint pairs are the only finite lists; no old trees or sampled tails are used. The `--start/--stop` arguments bound a half-open label range; skipping the serial quadrature assumes the first replay has already passed.
