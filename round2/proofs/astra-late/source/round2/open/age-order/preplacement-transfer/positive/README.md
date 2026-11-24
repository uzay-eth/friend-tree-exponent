> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# A near-one comparison and the limits of random-rank ordering

For every fixed `v≥3`, the original unconditional tails satisfy `Pr(Z_v>x)>Pr(Z_(v+1)>x)` for all `x` sufficiently close to one. The comparison below is actually uniform over the birth tree of the earlier vertex. It does not cross the intermediate-threshold gap, and it does not resolve the all-label comparison at zero. No finite-horizon tree enumeration is used.

## 1. An explicit comparison next to one

Write `ε=1−x`. Conditional on any size-`v` birth tree, the new leaf has a parent of degree at most `v−1`; its promotion at the very next jump has probability at least `1/[v(v−1)]`. At size `v+1`, on this event the marked vertex has one leaf and hence at least `n−v` leaves at size `n` if every intervening birth attaches to it. The conditional probability of this pure-append segment up to size `n≥v+1` is at least

```math
 \prod_{j=v+1}^{n-1}\frac{j-v}{j}
 =\frac{v!}{(n-1)(n-2)\cdots(n-v)}\ge v!\,n^{-v}.
```

The count deficit `v` in this estimate includes the focal vertex itself, so it cannot simply be replaced by `v−1`. At the terminal size the [tagged-leaf comparison](../../../../../../../astra/source/round2/fringe/hub-spectrum.md) gives `E[W|T_n]≥n−v−1`, using the gamma stock and its mean-one loss. With the actual Yule size amplitude `S`, `W=S Z_v`, `W≤S`, and `S∼Gamma(n,1)`. Thus `E(S−W)≤v+1`. No independence between `S` and `W` is used. For `n=ceil(8(v+1)/ε)` the gamma variance bound gives `Pr(S<n/2)≤4/n` and Markov gives `Pr(S−W≥εn/2)≤2(v+1)/(ε n)`; the sum is below `1/2`. It follows that

```math
 \Pr(Z_v>1-\epsilon\mid T_v)\ge
 l_v\epsilon^v,\qquad l_v=\frac{v!}{2v(v-1)[9(v+1)]^v}. \qquad\text{(1)}
```

For label `v+1`, first promotion has size `m≥v+2`. The stopped shifted urn in [fixed-threshold](../../fixed-threshold/README.md) gives a `Beta(3/2,m−3/2)` cap, since the marked vertex has a permanent core parent after promotion. For fixed red parameter, the tail decreases with the other parameter (couple its gamma denominator via independent increments). Consequently, without any estimate of the promotion probability,

```math
 \Pr(Z_{v+1}>1-\epsilon)\le
 u_v\epsilon^{v+1/2},\qquad
 u_v=\frac{\Gamma(v+2)}{\Gamma(3/2)\Gamma(v+3/2)}. \qquad\text{(2)}
```

Indeed integrate the beta density using `t^(1/2)≤1`. The half-power separation survives the correct append deficit in (1). In particular choose `0<ε<min(1,(l_v/u_v)^2)`. This proves strict adjacent comparison at these thresholds, including for arbitrarily large deterministic `v`, albeit with a window that shrinks with `v`. At the initial comparison, labels 1 and 2 agree by symmetry. For labels 2 and 3, there is a probability `1/2` that label 2 is the centre at size three. Starting there with two leaves gives (by the same append and gamma argument) a lower bound `const·ε` because `n−L_n=1` on the append path. The promotion cap for label 3 starts at size four, giving `O(ε^(5/2))` and therefore the same conclusion near one.

## 2. Exact interpolation by a random rank

Fix any threshold `x≥0` and denote the original deterministic-label probability by `p_j(x)`. Independent placement of the observation rank is important: a random rank chosen from future successful vertices is a different conditional law. For `m≥2`, `d≥1`, let `K_t=m+Bin(d,t)`, with the binomial draw independent of the tree. Then

```math
 -\frac{d}{dt}\mathbb E\,p_{K_t}(x)=
 d\sum_{j=0}^{d-1}\binom{d-1}{j} t^j(1-t)^{d-1-j}
       [p_{m+j}(x)-p_{m+j+1}(x)]. \qquad\text{(3)}
```

This follows by differentiating the binomial weights or by giving each of `d` rank increments an independent uniform activation time. The brackets are precisely the averages of the signed two-birth preplacement score (for `m≥3`) with its true normalization, rather than postbirth probabilities at a common exchangeable location. Positivity of every adjacent bracket implies monotonicity of every interpolation; conversely monotonicity for every starting size and any fixed `d` implies every bracket by taking the right derivative at zero. For a single interpolation over a block the reverse inference fails. For example the three degree-two Bernstein coefficients of `(t−1/2)^2+1/16` are `(5/16,−3/16,5/16)`: the polynomial is strictly positive despite the negative middle coefficient. Integrating those coefficients as adjacent differences after scaling by a positive constant can keep all four synthetic probabilities between zero and one. This is a precise obstruction to extracting all ranks from a single monotone binomial transform. It is not a tree counterexample.

## 3. The conditional signed martingale cannot preserve sign by degree order

Consider the signed terminal reward `1_(Z_a>x)−1_(Z_b>x)` on one common postbirth tree. Its conditional expectations form a bounded martingale. A deterministic nonnegative submartingale certificate based only on the two marked degrees and order of labels, with nonnegative value at `(1,1)` and limsup bounded above by the terminal reward along the future, would make the conditional gap nonnegative in every reachable postplacement tree with those two leaf degrees (by bounded reverse Fatou, assuming a bounded certificate). This is impossible even at threshold zero: the [large-star-with-an-arm fixture](../../README.md) has an adjacent older leaf with hub probability `<.031` and the younger leaf with hub probability `≥.125`. The same example rules out an order-preserving two-degree coupling from all such states; there are states where the younger leaf has the larger next gain rate even though its current degree is equal.

What is not ruled out is a certificate **before** either birth, taking expectation over the two placements and exploiting opposite signs in their weights, or a test of that signed score averaged over the law of the old tree. Formula (3) supplies no new sign without information on that average or an argument that works for every starting rank. The near-one result uses a different source of sign: the exact degree deficit created by immediate promotion limits the endpoint of the comparison urn for the later label. Together with the [compact birth-scale comparison](../../scaled/README.md) it leaves an intermediate region, including every fixed `x∈(0,1)` as `v→∞`, not decided by these arguments.
