> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# An effective common near-zero window for all labels

Use the corrected [all-threshold late cutoff](../../fixed-threshold/effective-uniform/README.md) `N=2^(2^(2^123))` in **old-size** indexing, and the all-tree hub-gap numerator `b=0.02968` from the combined discharge. Define the following integers and threshold (no expansion of `N` is needed):

```
k = 10 + 2 log2 N,     d = 2^k = 1024 N^2,
h = 16 k,              T = d (2d+h) 256 k,
epsilon = 2^(-2(T+h)) / (2(N+2)).
```

Then `Pr(Z_v>x)>Pr(Z_(v+1)>x)` for **every** `v≥2` and **all** `0≤x≤epsilon`. The two seed labels have identical laws. The implication uses the cited late cutoff and the hub certificate as inputs; it does not verify either of them anew. The new ingredient is a finite-size, birth-state-uniform loss bound using the true common size amplitude.

## 1. Uniform small *positive* continuous-time mass

Start the physical unit-clock process in a finite tree with a tag of degree one or two. Let `W=lim exp(−t)D_tag(t)`. The estimates below work for integer `k≥10`, `d=2^k` and the corresponding `h,T` above. A positive `W` forces a finite degree-`d` hit. Control a **finite late hit**, not the event of no hit.

Before a next gain, the tag degree is fixed below `d`. For each current neighbour `u`, run consecutive blocks of length `2d+h` until any gain at the tag or until `u` has `M=d^4` other leaves. The [finite-next-gain comparison](../../effective/cylinder/README.md#1-a-uniform-modulus-for-a-finite-next-gain) supplies conditional success probability at least `1/32` in each block. For clarity, if `u` is still a leaf, allow `d` units for its rate-one ring to send to the tag or for `u` to become a core. Then allow `d` units for the tag's channel of rate at least `1/d` to give `u` an other leaf, if it does not already have one. In the last `h` units track a birth–death family from one such leaf, with per-member rates `1,1/2`. Its mean is `exp(h/2)>2M`, and its second moment is at most four times the square of its mean, so Paley–Zygmund gives probability at least `1/16` of reaching `M` in that stage. All stages stop early on a tag gain and the floors multiply to more than `1/32`.

Once `u` has the other stock, the all-time first-send barrier bounds a later send from `u` before the next tag gain by `2/(M+1)`. Thus a specified neighbour's probability of supplying the next gain after `256k` blocks, but at finite time, is bounded by `d^-8+2/(d^4+1)`: use `(31/32)^32≤exp(−1)<1/2`. Take a union bound over at most `d` candidates at each of fewer than `d` successive gain stages, using the conditional bound at each stage's stopping time. A finite hit later than `T` requires one such long completed gap, so it costs less than `d^-6+2/d^2`.

At the finite hit, follow its fresh leaf for `h` units by dominating from below with a linear birth–death process having per-member birth rate one and death rate `1/d`. This is valid on the restart history since the tag's physical degree never decreases. Its probability of stock below `d` is at most `1/d+2d exp(−(1−1/d)h)≤1/d+2/d^7`. Indeed its extinction-by-age probability is at most `1/d` and the conditional positive count is geometric with parameter at most `2 exp(−(1−1/d)h)`; use the union bound for its first `d` possible counts. This is a finite-age stock comparison; its normalized birth–death limit is not used as the degree mass.

At stock at least `d`, the [tagged-leaf bound](../../../../../../../astra/source/round2/fringe/hub-spectrum.md#2-a-vertexs-limiting-mass-as-seen-at-time-t) on the true remaining degree mass is `Gamma(L,1)−Exp(1)`, with `L≥d` and no independence required. The probability this bound is at most `d/4` is at most `8/d`: compare the gamma marginal to shape `d`, apply Chebyshev at `d/2`, and use Markov at `d/4` on the exponential. If the hit was by `T`, the successful lower mass at the original starting clock is greater than `(d/4) exp(−(T+h))`. Conditioning at both stops and adding the errors gives, uniformly in the initial tree,

```
P(0<W<= exp(−(T+h)))
 <= d^(-6)+2/d^2+1/d+2/d^7+8/d < 10/d.
```

The normalization is always relative to the original start. In particular the estimate does not replace the limiting positive degree mass by a survival indicator of the auxiliary family.

## 2. True amplitude and conditioning

Start just after the marked birth at size `m`, conditional on its full birth tree, or at size three for the older seed label. In the unit-clock construction at size `j`, the holding time has law `Exp(j)` independent of the winning arrow; the size rates depend only on `j`. Equivalently generate the entire sequence of holding times independently of the entire embedded graph sequence. Hence `S=lim exp(−t)N(t)` is `Gamma(m,1)` conditional on the starting tree, and is independent of `Z` **under that conditioning**. For example `N(t)` is the sum of `m` independent Yule family counts, or its one-particle geometric law gives exponential limit and convolution gives the same gamma law. The identity is `W=S Z` on the actual physical path. It does not make `S` independent of `W`, or of quantities observed at a subsequent age or degree hit.

For `m≤N+2` the gamma mean gives `Pr(S≤2(N+2))≥1/2`. If `0<Z≤epsilon` and `S≤2(N+2)`, then `0<W≤2^(-2(T+h))≤exp(−(T+h))`, using `exp(1)<4`. Independence here is only the size-clock/embedded-chain independence from the deterministic-size start. Thus uniformly in the birth tree,

```
Pr(0<Z≤epsilon | birth tree) < 20/d = (5/256)/N^2.
```

This bound does not condition on the mark ultimately being promoted. In particular it does not silently assign a birth-size gamma law to a process restarted at a degree hitting time.

## 3. Transfer the hub gap

At a common threshold `x`, write `F_v(x)=Pr(0<Z_v≤x)`. The adjacent tail gap is its hub gap minus `F_v(x)` plus `F_(v+1)(x)`. It suffices to pay the older label's loss once. For an old size `n<N` the supplied numerator gives a hub gap at least `.02968/[n(n+1)]≥.02968/N^2` once in the general certificate range. The low-mass loss above is strictly smaller (`5/256<.02968`). The [small-base hub gaps](../../early/README.md) for pairs `(2,3)` and `(3,4)` exceed `5/256` already; for `(2,3)` start at size three, where label two has degree one or two. The same loss estimate therefore pays the bases. For `n≥N`, invoke the corrected late theorem, which already covers the full interval `0≤x<1`. This proves the common near-zero window with no estimates of how fast an individual label approaches its atom at zero.

Replay the compact arithmetic with `python3 -IB round2/open/age-order/global-transport/near-zero/check.py`. This checks the numerical margins; the clock and branching comparisons are analytic.
