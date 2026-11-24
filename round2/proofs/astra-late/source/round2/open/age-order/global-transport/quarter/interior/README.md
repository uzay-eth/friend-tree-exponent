> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# Two interior birth-scale bands from physical clocks

Write `m=b+1`, where `b` is the older of the next two births, and `y=mx`. The signed sharp-tail comparison is strictly positive on **every** old tree on the following closed boxes:

| Birth-scale interval | Common completion size |
|---|---:|
| `[1/10,21/200]` | `m ≥ 14,014,000` |
| `[1,1001/1000]` | `m ≥ 240,002,400,000` |

There is also a finite-size extension for `b=22,23,24`: `11/5 ≤ y ≤ (b+1)/10`, adjoining the frozen tenth-threshold result. Neither band is an all-size or all-label claim. In particular the lower band does not join the tiny birth window.

## True-amplitude kernels

From a fixed restart with size `n`, the rate-one-vertex clock gives `W=S_n Z`, with `S_n~Gamma(n,1)` independent of the embedded path. The [Poisson-stock lower](../../../hub-ranking/poisson-tail/README.md) becomes a Poisson mixture of `Gamma(j-k,1)` tails for `0≤k<j` and zero for larger `k`; this follows by beta–gamma algebra at the fixed restart. The analogous locked beta upper gives `Gamma(r,1)` at a pending tip (`r=1/k` for parent degree `k`), and `Gamma(1+r,1)` just after physical promotion. The latter also follows from the affine pure-birth coupling. At physical stopping times restart the comparison relative to the current time, and multiply the *required amplitude* by `exp(elapsed time)`. No independence of the original size amplitude and a clock-defined stopping event is used.

For shape `1+r`, the survival upper is `e^-z (1+z)^r/Gamma(1+r)` by concavity inside the exponential overshoot integral. Shape `2+r` follows by gamma recurrence. For shape `r≤1`, the upper `e^-z z^(r−1)/Gamma(r)` follows by monotonicity of the density's power factor. At high degrees `Gamma(r)≥1/[r(1+r)]≥k−1` by integrating `1−t` on `[0,1]`. Elsewhere a rational lower for `Gamma(1+r)` integrates concave chords against the tangent lower for the exponential density on `[0,10]`. The upper includes the [two-race postpromotion improvement](../../../hub-ranking/tip-moments/README.md): candidate gain at rate `1+1/k` pays shape `2+1/k`, the tag's own ring splits into two rate-`1/2` channels paying at most shape `1` and `1+1/(k+1)`. The three channels are independent up to the first candidate or tag ring; this is a bound at the actual promotion history, not a moment-to-tail implication.

Before promotion retain actual sibling stock with rates `1/k` for promotion, `j/k` for sibling promotion and `[j+1,(k+j+1)/2]` for parent growth. Use the first-promotion barrier `1/(j+1)` as a further upper, along with the locked cap. The lower tracks stock: with marked degree `d` and positive tracked stock `j`, use gain rates in `[j,d]` and loss rate `j/d`. Upon last loss, the just-promoted former leaf starts at degree two. A specified return arrow has rate `1/k` from that neighbour, its growth rate is at most `k−1+1/d`, and other mark gains are disjoint and at most `d`; a return or other mark gain starts a new one-leaf row. At return cutoff use zero. At the positive-stock degree cutoff, and as a fallback at lower degrees, use the nonnegative Poisson–gamma payoff. The favourable parallel pending tip is the same return row with `d=1`. The serial host bonus keeps the original old parent through one race: gain at `1+1/k`, loss at `1/2`, host growth in `[0,k−1/2]`. On loss use the physical zero-stock return row. On host growth continue; at cutoff use the unboosted serial row. All these comparisons retain the appropriate actual degree and do not suppress an adverse physical arrow.

## Time mesh and translation to shares

Use log-time step `h`, evaluating at a rational outward enclosure of `z exp(ih)`. An upper row is made decreasing by lifting earlier entries, and a lower row by lowering later entries: each operation uses monotonicity of the underlying threshold event. In each slab take future branch payoffs at the early end for an upper bound and late end for a lower bound. With fixed branch payoffs the backward-time equation has affine drift `P−R v` at each rate corner. If `w` is a boundary value, the two candidate endpoints `w+h(P−Rw)` and `(w+hP)/(1+hR)` bound the linear barrier: at either end of the segment its slope is outward from every allowed affine drift when the maximum (upper) or minimum (lower) is taken over both endpoint forms and all corners. A separate valid bound can replace an endpoint, e.g. the locked cap or Poisson fallback. The postpromotion fixed-rate race instead uses outward exponential masses and monotone rectangular payoffs. Lower horizon discards unresolved paths except for the stated fallback; upper horizon uses the physical locked/stock bound. Root, exponential and all grid rounding use rationals.

From the common completion only, `S_m` has mean `m` and gamma Chernoff tails. For `0<e<1`, either one-sided relative deviation by `e` costs at most `exp(-m e²/[2(1+e)])` (integrate the logarithm in the usual gamma Chernoff exponent). Subtract this error from the `W>(1+e)y` lower and add it to the `W>(1−e)y` upper. This is a union bound, including at the strict threshold, and does not require independence between `W` and `S_m`. On a closed `y`-box use the right threshold for lower kernels and the left threshold for uppers. The replay uses `e=.001` and `e=.00001` in the first and second boxes respectively and ensures deviation cost at most `e` in each.

## Signed tree ledger and replay

The [authenticated algebraic ledger](../../../hub-ranking/mean-ledger/fractional/README.md) uses these sharp-tail prices in place of endpoint moments; only its local signed bookkeeping is imported. The high-degree branch uses only the stock caps `1/(j+1)`, which are still valid. The two tariff rows are `(a,b,fee,terminal fee,transfer)=(.011,.011,.007,.020,.003)` and `(.011,.011,.007,.020,.002)`, with degree split 160 and the high chain envelope starting at 48. Thus on each box the unnormalized signed sharp-tail numerator is at least `.033` for trees with at least two cores, before division by `(m−2)(m−1)`. Stars and the seed have strictly positive separate checks.

```sh
R=round2/open/age-order/global-transport/quarter/interior/check.py
python3 -IB "$R" 1/10 21/200 --steps 430 --scale 129/128 --dcut 12 --rcut 20
python3 -IB "$R" 1 1001/1000 --err 1/100000 --steps 240 --scale 129/128 --dcut 12 --rcut 20
python3 -IB "$R" --finite
```

The finite replay uses the frozen quarter checker, with `2,4,7` threshold boxes at labels `22,23,24`. Trying the same left endpoint `y=11/5` at `b=25` failed under the frozen tariffs at `[11/130,5633/66560]`: an outward length-one chain bill exceeded baseline by `4877580950041/30000000000000000`. This is a tested certificate failure, not a physical counterexample or a claim that a sharper subdivision or local allowance cannot work. Likewise the size cutoffs in the table are sufficient only; a size-ratio lower-size transfer remains open.
