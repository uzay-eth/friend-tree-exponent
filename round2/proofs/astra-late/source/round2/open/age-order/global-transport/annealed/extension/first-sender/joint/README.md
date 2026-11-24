> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../../README.md#claim-26) for exact scope and citations.

# Sub-half cutoffs and the joint red-family kernel

Strict adjacent tail order holds on every pre-placement tree in the following bands (`b` is the older birth label, and `x<1`):

| Older label | Threshold |
|---:|---:|
| `b≥20` | `x≥0.49` |
| `b≥21` | `x≥0.48` |
| `b≥40` | `x≥1/3` |

With [Fable's paired receipts and strictness](../../../../../../../attacks/fable5/README.md), the first row gives the all-label interval `[0.49,1)`. Including the single `20→21` receipt gives `[0.48,1)`. Founders are tied. At `1/3` the pairs with older label `21,…,39` remain unproved by these certificates (include `20` if using only paired receipts). The numerical advance here uses the first-family marginal; the joint kernel below is exact but has not been evaluated to certify a further band.

## 1. Sharpen the frozen one-family proof

Use the [frozen half-threshold construction](../refined/README.md), including its dummy charge on a miss. Its post-first-sender loss coefficient can be tightened by keeping the `j+1` in the denominator. For first sender `F=f`,

```
b Σ_(j>f) 1/[(b+j)(j+1)]
 ≤ (1+1/(b+f+1)) b Σ_(j>f) 1/[(b+j+1)(j+1)]
 ≤ (1+1/(b+f+1)) log(1+b/(f+3/2)).
```

The second line uses the midpoint integral of the convex function `1/[v(b+v)]`, starting at `f+3/2`. This bound works also at `f=0`. It replaces `C` in the payoff `(1−C(1+a/b)/(t−a))_+` on `a<t` from the frozen proof.

For a band `b≥B`, use a common uniform quantile `q`, put `h=(1−q)/q` and `f_b=floor((b−1)h)`. At `b=B` use `f_B` exactly. For `b≥B+1`, use `f_b≥g_b=max(0,(b−1)h−1)`, so `b+f_b+1≥B+2+g_(B+1)`. For the scaled denominator use the untruncated bound `(f_b+3/2)/b≥h+(1/2−h)/b≥min(h,(Bh+1/2)/(B+1))`. Taking the worse of the resulting charge and the exact-`B` charge covers every `b≥B`, including bins on both sides of `h=1/2`. Drop the last `q` bin, whose right edge is `h=0`.

The uncapped family-size survival is `(1+a/b)^(-(f_b+1))`. The inequality `log(1+y)≥2y/(2+y)` yields the exponential upper survival exponents

```
b=B:     v0=(f_B+1)·2a/(2B+a),
b≥B+1:   v1=[B/(B+1)]·h·2a/(2+a/(B+1)).
```

The second bound uses `(f_b+1)/b≥(1−1/b)h`. Use `min(v0,v1)`, replace `t` by the same mean-one exponential as in the frozen proof, and sum outward rational rectangles. The unchanged paired beta debit gives each interval by monotonicity from its corner.

## 2. Exact multi-hit update, including misses

At any checkpoint just before a tagged blue sender, let `p` be the absolute sender-word probability of uncharged reds (initially `θ`), let `r` be their current number, and let all of them be good leaves. The checkpoint may retain the actual marked degree `D` and local arrows. Each live vertex is a distinct family founder. The limiting shares *within that live colour* have law `Dirichlet(1,…,1)` conditional on the checkpoint: future live steps select uniformly among the live vertices and a newborn inherits the sender's family. For any `k<r` distinct founders chosen independently of future drivers, their shares and the aggregate remainder are jointly `Dirichlet(1^k,r−k)`. In particular distinct chosen shares have covariance `−1/[r²(r+1)]`. Independence or repeated Beta draws at fixed `r` would miss this dependence.

The actual arrow has these branches:

* **Miss**, probability `1−r/D`: expose the actual non-live neighbour if the subsequent geometry is needed. Delete no family, and leave `p,r` alone.
* **Hit**, probability `r/D`: the lost leaf is a uniformly chosen live founder. Expose its asymptotic fraction `V` with `Beta(1,r−1)` for `r≥2` and `V=1` for `r=1`. Mark its entire family as charged. Continue with `p'=p(1−V)`, `r'=r−1`. Every uncharged vertex is still a good leaf. The normalized surviving family shares are independent of `V` by Dirichlet neutrality.

Iterating this at the first `K` blue selections of the mark gives a joint disintegration even when the next victim is a **new** leaf. After conditioning on the exposed aggregate percentages and the finite arrow history, relabel each current live vertex as a founder at the next checkpoint. To see why this is valid, conditional on the previously exposed family shares the live versus charged word is iid, but the live within-colour sender driver is fresh and uniform; it was not used to sample any previously exposed share within the live class. The new founder distribution is consequently Dirichlet with unit counts at the actual new census, not the Dirichlet vector from the first census. Between checkpoints every live sender produces a new leaf at the mark. Only marked sends can spoil one. This establishes the invariant needed by the branches. On a miss the potential surviving shares are not conditioned on the destination of a *future* arrow. Replacing `r/D` by a fixed probability independent of the exposed stock is not justified.

## 3. A residual bound valid at each finite checkpoint

After the last exposed tagged arrow, say at blue index `f`, let `p,r` be the live parameters, and keep the original blue probability `u=1−θ`. If `r=0` or `p≤x`, use a zero lower success floor. Otherwise the live count before blue `f+k` dominates (in fact equals) `r+S_k`, where `S_k` counts live choices in `k` independent strings after deleting the already charged red choices. Its summands are geometric with survival parameter `p/(p+u)`. This remains true given the finite past and exposed aggregate shares; it does not require an independent degree model. The bad-fraction estimate gives

```
E[(p−Z)_+ | checkpoint] ≤ p Σ_(k≥1) E[1/(r+S_k)]/(b+f+k).
```

Put `λ=p/(p+u)`. The reciprocal is `∫₀¹ v^(r−1)[(1−λ)/(1−λv)]^k dv`. For `k≥2`, it is at most `u/[p(k−1)]`, from the case `r=1`; at every `k` it is at most `1/r`. After summing through `k=L≥1`, the remaining mean-deficit charge is at most `u·log(1+(b+f+1)/L)/(b+f+1)` plus `u/[L(b+f+1+L)]`, by bounding the first term and the remaining decreasing sum by an integral. Markov with slack `p−x` converts this into a lower success payoff; if the total is too large, use zero. In particular this supplies a convergent scalar error bound for the **residual** at each fixed positive-slack checkpoint, regardless of any error control chosen for the distribution of those checkpoints.

There is a precise obstruction to seeking a size-only uniform positive miss probability for the first few threats. On the physical event that the first `M` steps are red, there have been no blue sends, all red vertices are good, and `D=r+1` (the only other neighbour is the old parent). If the very next step is the first marked blue send, its conditional miss probability is exactly `1/(M+2)`. This event has positive probability for every finite `M` and any serial seed. Thus no uniform positive miss credit survives conditioning on the exact word, even at the first threat. More generally at that first threat, `r=R` and `R+1≤D≤R+1+f`: before it, only the `f` blue births can add non-red neighbours beyond the old parent. The miss probability therefore lies in `[1/(R+1),(f+1)/(R+f+1)]`. These bounds are conditional on the census and can be integrated against the negative-binomial law, but the upper end cannot be credited as a lower probability. This is not an upper bound on the complete multi-hit approach. Making it quantitative requires averaging the stock-dependent hit/miss kernel on reachable checkpoints, or proving a valid comparison for the joint stock/degree law; the marginal rectangles in §1 do not do that.

The first miss can be included without solving the degree law if one uses the lower end of this bracket. For any fixed payoff depending only on the post-arrow `(p,r)` and nondecreasing in both arguments (the residual success floor is one), use the genuine miss branch with weight `1/(R+1)` conditional on `R`, and the conservative dummy charge on the remainder. This follows by thinning genuine misses with probability `D/[(D−R)(R+1)]≤1`, using an independent uniform at the checkpoint. For example the resulting retained miss mass conditional on `θ,f` is exactly

```
∫₀¹ v [u/(1−θv)]^(f+1) dv.
```

For a weighted payoff the same factor must stay inside the sum over `R`; the retained miss census is tilted by `1/(R+1)`. This gives a scalar way to credit real misses at the first threat. Iterating with only this scalar census needs additional control, because the physical degree after a real hit or miss enters the next threat's branch probability.

## Replay

```
nice -n 10 python3 -IB round2/open/age-order/global-transport/annealed/extension/first-sender/joint/check.py
PASS b>=20, x>=49/100: coefficient > 521/2500
PASS b>=21, x>=12/25: coefficient > 2087/10000
PASS b>=40, x>=1/3: coefficient > 261/1250
```
