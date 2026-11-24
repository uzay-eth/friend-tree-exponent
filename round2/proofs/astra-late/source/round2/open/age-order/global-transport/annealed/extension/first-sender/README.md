> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../README.md#claim-26) for exact scope and citations.

# First-sender averaging: a common threshold of 0.591

Let `b=n+1` be the older label in the pair `b→b+1`; the pre-placement tree has size `n`. For `0<x<1` the following are **strict pre-placement bounds on every old tree**, and hence also hold annealed:

| Older label `b` | Sufficient threshold |
|---:|---:|
| `b≥20` | `x≥151/250 = 0.604` |
| `b≥21` | `x≥591/1000 = 0.591` |
| `b≥30` | `x≥1/2` |

The [two Fable receipts](../../../../../../attacks/fable5/flow/README.md) cover the older labels `2,…,19` at **all** thresholds; the later, authenticated **single** receipt covers `b=20`. Strictness on `0<x<1` follows from the positive core-to-leaf transport mass as in the [audit](../../../../../../attacks/fable5/README.md). Thus `0.591≤x<1` is an all-label interval using the single receipt, and `0.604≤x<1` works with the paired receipts alone. Founders 1 and 2 have equal laws. At `x=1` all tails vanish. At `x=1/2` the uncovered pairs are exactly the older labels **21 through 29 for these certificates** (20 through 29 if requiring two Fable receipts). These are sufficient cutoffs, not identified optimal thresholds. The proof below certifies a whole interval by monotonicity of its bounds, not by assuming monotonicity of the tail *gap* in `x`.

## Replace the fixed protection window by a disjoint first-sender decomposition

At a serial middle completion the common size is `m=b+1`. Colour its one fresh leaf red and all other vertices blue. Each birth inherits its initiator's colour. As in the [protected colour argument](../../../sharp-large/README.md), the physical sender colour word has the exact `Beta(1,b)` mixture; conditional on endpoint `θ` the colours are iid. Write `u=1−θ`. Let `F` be the index, starting from zero, of the first blue choice of the mark as sender. Given the *entire* colour word the blue sender hazard at index `f` is `1/(b+f)`, until the first choice. Consequently `F` is independent of the word and satisfies

```
P(F≥k) = (b−1)/(b+k−1),    k=0,1,... .
```

Crucially, a fixed prefix avoided the loss at its last protected step; conditioning on `F=f` instead forces the loss opportunity at `f`. Let `R_j` count red vertices just before blue choice `j`, after the corresponding red string. No red leaf family is lost before `f`, at most one red leaf is lost at `f`, and subsequent blue choices select the mark with probability `1/(b+j)`. Good reds reproduce good reds; thus on red steps the expected bad fraction cannot increase. Conditional on `F=f` and `θ`, the loss bound by the same finite-time bad-fraction estimate and Fatou is

```
E[(θ−Z)_+ | θ,F=f] ≤ θ E[1/R_f | θ]
                      + θ Σ_(j>f) E[1/R_j | θ]/(b+j).
```

After step `f` this conditioning exposes no future within-colour drivers. The geometric word identity `θ E[1/R_j | θ]=u(1−u^j)/j` for `j≥1` remains valid. For `f≥1`, dropping `1−u^j` and integrating the decreasing remaining summand from `f` gives the simple charge

```
E[(θ−Z)_+ | θ,F=f] ≤ (u/b) a_b(f),
a_b(f)=b/f+log(1+b/f),   a_b(0)=∞.
```

Here the convention at zero means we use a lower success probability of zero, not that the physical loss diverges. For `θ>x`, put `t=b(θ−x)/(1−θ)`. Conditional on the beta endpoint exceeding `x`, the survival of `t` is exactly `(1+t/b)^(-b)`, independent of `x`. Markov at each `(F,t)`, retaining the entire positive part, yields

```
h_M(x) ≥ (1−x)^b E[(1−a_b(F)/t)_+].                     (1)
```

This averages all disjoint first-sender classes and all endpoint shifts. It does not select a single protection window or estimate a physical tail from finite simulations.

## A monotone coefficient and one rational corner per band

Generate `F=floor((b−1)(1−q)/q)` from `q` uniform on `(0,1)`. For `b≥B`, set `g_B(q)=(B−1)(1−q)/q−1`. When positive, `F/b≥g_B(q)/B` since `g_b(q)/b` is increasing in `b`. Hence `a_b(F)≤h+log(1+h)` with `h=B/g_B(q)`. Set the success floor to zero elsewhere. Also `(1+t/b)^(-b)≥exp(−t)` so a mean-one exponential can replace `t` in the increasing payoff of (1). This gives a uniform coefficient `c_B` by integrating over independent uniform and exponential variables. The replay lower bounds this integral with rectangles: use the right `q` endpoint for the charge, the left `t` endpoint for the payoff, and a rational lower bound on exponential bin mass. Drop boxes with `h≥12` and the `t>12` slice. For `log(1+h)` use the first eight terms of the positive series in `r=h/(2+h)`, plus the geometric upper remainder `2 r^17/[17(1−r²)]`. The exponential bin masses use `e^-z≥1−z+z²/2−z³/6` and `1−e^-z≥z−z²/2` at the grid spacing `z`. All score rounding is outwards.

On the adverse side use the [one-third paired-edge debit](../README.md), retaining the factors `b+1/2` in the high beta cap and `(1−x)^(1/6)` in the core-to-leaf density comparison:

```
U=(1−x)^b sqrt((1−x)/(π (b+1/2)x)),
Q_(m,1/3)(x) ≤ θ_* U,
θ_* = (3/4) ((1−x)/(bx))^(1/6)   (in the tested bands).
```

The `3/4` bound holds for `b≥20`: use `sqrt(π)<1773/1000` and the Gamma(1/3) integral lower bound from the parent proof. The signed two-birth numerator divided by its positive serial weight is at least `(1−x)^b c_B−U(4/3+2θ_*/3)`. The normalized cost decreases in `x` and `b` and the coefficient floor is valid for **all** `b≥B`. The replay checks the three rational corners with a common outward sixth-root bound `57/100`. There is no enumeration of old-tree laws and no assumption about an annealed degree profile.

## Where this leaves the nine pairs

At `b=21,x=1/2` the exact integral for the **chosen first-sender charge** `a_b(f)` in (1) is less than `0.18`, even before the lower-replay relaxations. For this upper check use `F≤(b−1)(1−q)/q`, `log(1+h)≥2h/(2+h)`, exact rational `(1+t/b)^(-b)` bin masses, and upper rather than lower payoff rectangles; the omitted high-`t` payoff is bounded by its tail mass. The chosen normalized paired beta debit at that corner exceeds `0.19`, even evaluating with exact roots and π. So just refining the first-sender quadrature cannot settle the first missing pair. This obstruction is to the stated charge-and-debit certificate only: it does not bound the real serial tail above or the real gap below. In particular the bound `a_b(f)` forgets both the factor `1−u^j` and the chance that the forced sender hits a good red leaf. Further progress at these labels would need sharper physical loss/debit information (for example retaining where that first arrow goes), or additional actual transport receipts. None of the band edges here supplies an all-threshold receipt for any of those nine pairs.

## Replay

```
nice -n 10 python3 -IB round2/open/age-order/global-transport/annealed/extension/first-sender/check.py
PASS b>=20, x>=151/250: middle coefficient > 327/2000
PASS b>=21, x>=591/1000: middle coefficient > 1639/10000
PASS b>=30, x>=1/2: middle coefficient > 167/1000
PASS the chosen scalar certificate cannot cover b=21, x=1/2
```
