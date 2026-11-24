> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../README.md#claim-18) for exact scope and citations.

# Birth labels and limiting hub mass in the original random friend tree

For the original edge seed, $`\Pr(Z_1>0)=\Pr(Z_2>0)>\Pr(Z_3>0)>\cdots`$. Every positive real moment of the limiting degree share has the same strict order. Both adjacent comparisons hold conditional on the tree before the next two births. [Original Question 3](https://arxiv.org/html/2403.20185v1#S7) · [Hubhood proof](hub-ranking/README.md) · [Positive moments](hub-ranking/all-powers/positive/README.md)

For **every $`1/10\le x<1`$**, the original seed has $`\Pr(Z_1>x)=\Pr(Z_2>x)>\Pr(Z_3>x)>\cdots`$. The comparison is conditional on every common pre-birth tree for older labels $`b\ge20`$. The early unconditional comparisons use Fable’s eighteen paired receipts. At smaller thresholds the next two births from any old tree of size $`n`$ are strictly ordered if $`(n+2)x\le2^{-66}`$ or $`(n+1)x/(1-x)\ge6`$. All sufficiently late adjacent labels are strictly ordered for every $`0\le x<1`$. [Threshold ledger](global-transport/quarter/README.md) · [Label 20](global-transport/quarter/b20.py) · [Cold audit](../../../../../../reviews/cold/astra/open/age-order/hub-ranking/beta-tests/physical-tail/poisson/quarter.md) · [Small birth-scale window](global-transport/birth-window/sharp/README.md) · [Poisson–beta comparison](hub-ranking/poisson-tail/README.md) · [Late labels](eventual-domination/README.md)

The earlier argument below compares each fixed label with all sufficiently late labels at every threshold. Its conditional reversal fixes a tree after both placements.

## What the adjacent swap does and does not cancel

In the original edge seed the labels 1 and 2 are exchangeable, including their entire infinite trajectories. Thus `Z_1` and `Z_2` have identical laws, and their hub probabilities coincide. Neither a uniform-age prior nor a conditioned parent is involved.

At size `n`, on a fixed old tree `T`, write `a(x)=Σ_{z∼x} 1/D_z`; the next target is `x` with probability `a(x)/n`. Suppose labels `u=n+1` and `v=n+2` both attach to vertices in `T`. Conditional on their **unordered pair of old targets with multiplicity**, the two resulting placements give the same unlabelled graph, with the positions of the two new labels interchanged. For distinct targets `x,y`, put `δ_x=1/[D_x(D_x+1)]`. If `x∼y`, the probabilities of the two orders before conditioning are

```
p_xy = a(x)(a(y)−δ_x)/(n(n+1)),
p_yx = a(y)(a(x)−δ_y)/(n(n+1)).
```

For nonadjacent targets both correction terms vanish. For a repeated target the two newborns are twins, hence conditionally exchangeable. On an adjacent edge the sign of `p_xy−p_yx` is the sign of `a(y)δ_y−a(x)δ_x`. It can go in either direction relative to degree: the two exact fixtures have `(D_x,D_y)=(3,1)` with `p_xy<p_yx`, and `(3,3)` with `p_xy>p_yx`. A common conditional law for the **locations** therefore holds on the nonadjacent or repeated-target sector, not in general. Once the physical locations are fixed the future law is independent of the age marks, but the orientation bias still multiplies the difference between the two location-specific endpoint laws. If instead the later newborn attaches to the earlier one, only one chronological order exists; swapping the labels in the unlabelled length-two extension is not another allowed two-birth trajectory. An argument for unconditional order needs both this sector and an averaged sign for the adjacent-target sector.

The finite checker also finds consecutive labels 5 and 6, both leaves at size six, with respective incoming intensities `1/3` and `1/2`. A conditional upward-rate order based on degree and immutable age alone does not survive the next step. This does not decide the ordering of their unconditional degrees or limits.

## Tagged leaf bounds used for the limiting comparisons

Use independent rate-one initiating clocks and let `N(t)` be the size. Started at size `n`, `e^(−t)N(t)→S` with `S∼Gamma(n,1)>0`. Physical waiting times depend on size alone, so the discrete tree path (including its limit) is independent of the sequence of size waiting times. **It does not follow** that the continuous-time local mass `W=lim e^(−t)D(t)=S Z` and `S` are independent.

At a core with `k` selected neighbouring leaves, tag them and their copies born by their own rate-one rings. Each surviving tag is promoted only by the core, at per-tag rate `1/D≤1/2`. Two useful couplings are:

1. Globally couple a birth-death population (per-member rates `1,1/2`) below the tag count. Its survival probability is at least `1/2` from one and `1−2^(−k)` from `k`. Survival implies unbounded physical degree, which implies a hub as follows. At degree first hitting `d`, a newly acquired leaf was just born; preserve it with an embedded birth-death population of rates `1,1/d` to obtain unbounded leaf count with failure probability at most `1/d`. At first leaf count `m`, the next coupling yields hub probability tending to one. Taking these stopping times to infinity proves the implication (cf. [hub-children](../../../../../astra/source/round2/fringe/hub-children.md)).
2. Couple the tags to a rate-one-per-member Yule population, keeping a phantom Yule family whenever a tag is removed. Let `Y_t` be this upper population and `B_t` the living physical tags. A removal at time `r≤t` contributes expected phantom count `e^(t−r)` at time `t`, conditionally at removal, and the predictable removal intensity is at most one. Thus `E[e^(−t)(Y_t−B_t)]≤1`. The Yule limit is `G_k∼Gamma(k,1)`. Since `B_t≤D_t`, Fatou gives `E[(G_k−W)_+]≤1`. In particular `E W≥k−1`, and `P(W=0)≤P(G_k<k/2)+2/k`, tending to zero as `k→∞`. This completes the first argument without assuming independence of the loss and the gamma limit.

For a convenient uniform reciprocal estimate, stop instead on the first hit of `m<k` tags. Ignore any one marked neighbour outside the stock, so tags are unsullied by its eventual promotion. While the number `B` of tags is at least `m+1`, its gain rate is `B` and its loss rate is at most one. The function `1/B` has drift at most `−1/(B+1)+1/[B(B−1)]`: an independent linear-in-time upper envelope will not control the cumulative hazard at infinite time, but a lower **immigration-killed Yule** stock does. Construct `U` with transitions `j→j+1` at rate `j`, `j→j−1` at rate one, from `k`, absorbed at zero. The same clock and nested gain rates give `U≤B`. Its martingale is

```
M_t = e^(−t)(U_t−1)  (before absorption),
E[(M_(t∧τ)−(k−1))²] ≤ k+1/2,
```

where `τ` is absorption. In fact `E U_t≤k e^t` and the quadratic variation is bounded by `∫ e^(−2t)(k e^t+1)dt`. Stop the martingale at the first violation of `U_t ≥ (k/2)e^t`; such a violation requires a downward excursion of size at least `k/2−1`. The probability of this event is at most `(k+1/2)/(k/2−1)^2`. Thus with this exception removed, the tag stock is at least `(k/2)e^t` for all `t≥0`.

## A real *conditional* reversal, at adjacent labels

At size `n=k+4`, prescribe the reachable labelled tree whose vertex 1 has neighbours 2 and labels `3,…,n−1`, while vertex 2 has neighbours 1 and `n`. Let `u=n−1` and `v=n`. Both are deterministic labels and both are still leaves; `u` is a leaf of the degree-`k+2` vertex 1, and `v` is a leaf of the degree-two vertex 2. This conditioning event has strictly positive probability (every specified attachment has positive rate). The `k` other leaves of vertex 1 supply the stock above. Before `u` promotes, its promotion intensity is exactly `1/D_1`; on the good stock trajectory its cumulative intensity is at most `2/k`. More formally, integrate the promotion intensity stopped on the first lower-envelope failure, then add that failure probability. It follows that

```
P(Z_u>0 | T_n) ≤ (k+1/2)/(k/2−1)^2 + 2/k.
```

For `v`, race vertex 2's first ring selecting `v` (rate `1/2`), `v`'s own ring (rate one), and input from vertex 1 to vertex 2 (at most `1/2`). Until one of these events vertex 2 keeps degree two. Using a dominating independent clock for the last input gives probability at least `1/4` that `v` is promoted first. Then its new leaf initiates a surviving birth-death population with conditional probability at least `1/2`. Hence `P(Z_v>0|T_n)≥1/8`. With `k=200` the older upper is `29851/980100 < .031`, while the newer lower is `.125`. This is an asymptotic conditional reversal (and rules out conditional stochastic dominance), not a finite-size crossing. It does not disprove the unconditional comparison because the other trees at size `n` may contribute the opposite sign. It also illustrates why comparing by deterministic labels is not the same as conditioning a new vertex on its parent or the parent's future hub status.

## Unconditional partial comparisons and beta/Yule conditioning

- **Hub probability tends to zero.** From size `v≥3`, a vertex `v` still a leaf has permanent parent. For any cutoff `m≥v`, the probability it promotes by size `m` cannot be bounded by a small unconditioned age hazard merely using `1/D_parent≤1`. Instead use the identity at size `v−1`: the newborn target distribution is `a(x)/(v−1)`. A parent of degree `<K` has `a(x)≤D_x<K`. Every core of degree `<K` contributes at most `K` to the total target weight; all leaves together have target weight `Σ_{core x}L_x/D_x≤C_{v−1}`. Therefore this low-parent-degree event has probability at most `(K+1)E C_{v−1}/(v−1)`. On a parent of degree `≥K` with `k≥K/2` leaf neighbours at birth, excluding the newborn, the stock argument bounds the child's ever promotion by `O(1/K)`. Parents with fewer leaves than half their degree contribute total target weight at most `2Σ_x c_x≤4C_{v−1}`: for a core `a(x)≤D_x`, and `D_x≤2c_x` in this sector. Thus, as `v→∞`, choose `K→∞` slowly with `K E C_{v−1}/(v−1)→0`. The requisite sublinear mean follows already from the paper's almost-sure bound by dominated convergence. With the repository's mean power estimate the same reasoning gives `P(Z_v>0)=O(v^{-(1−λ)/2})` for any bound `E C_n=O(n^λ)`, `λ<1`. This is a tail bound, not termwise monotonicity.

  **Sharper rate.** No degree cutoff is needed. For `k≥4`, the stock bound above satisfies `b(k)≤20/k`: `k/2−1≥k/4`, so the two terms are at most `(18+2)/k`. For `1≤k<4`, use the same larger bound trivially. In the sector `k≥D/2` (where `k` is the parent's leaf stock *before* the newborn attaches), `a(x)b(k)≤D·20/k≤40` for each possible core parent. In the complementary sector `k<D/2`, charge `a(x)≤D≤2c_x`, and charge actual leaf parents in aggregate by at most `C_n` as before. Consequently

  ```
  P(Z_(n+1)>0) ≤ 45 E C_n/n,   n≥3.
  ```

  This improves the power to `O(n^(λ−1))`; the bound applies to deterministic birth labels by averaging over all parent targets.
- **Fixed old labels have positive mass with positive probability.** For every deterministic `u`, force a finite positive-probability sequence giving `u` at least one leaf neighbour and degree at least two, and apply the tag survival bound. Therefore for each `u` its hub probability strictly exceeds that of all sufficiently late labels.
- **Quantitative beta upper cap.** A fixed vertex has attachment probability at most `D_n(u)/n`. Couple its degree to a two-colour Pólya count, started at its deterministic birth size `v≥3` with one red and `v−1` white. At equality use the common uniform to keep the degree below the red count; while strictly lower a one-step move cannot overtake. Thus `P(Z_v>x)≤(1−x)^(v−1)` for `0≤x≤1`. For either original label use `Beta(1,1)` at size 2. More generally after stopping at a specified state `(n,d)` the corresponding conditional cap is `Beta(d,n−d)`. In continuous time this beta law is the ratio of independent gamma limits for the **comparison urn**. It is not the exact law of the physical limit given a target or a parent. Nor is it an upper bound below one on `P(Z_v>0)`.

### Upgrade to simultaneous stochastic order for sufficiently late labels

Fix `u`. Prescribe a reachable state at size `m` in which it is a core with `k₀≥1` leaf neighbours. Put `B=m−k₀≥1` and let `π>0` be the probability of the prescribed history. Then insist that every following birth up to size `n` attaches to `u`. Conditional on this history, its leaf stock is `n−B`; the probability of reaching this pure-append state is at least

```
π ∏_(j=m)^(n−1) (j−B)/j
 = π (m−1)!/(m−B−1)! · (n−B−1)!/(n−1)!
 ≥ π (m−1)!/(m−B−1)! · n^(−B).
```

At that state the gamma-minus-loss coupling gives `E W_u≥n−B−1`; since `W_u=S Z_u≤S`, we have `E(S−W_u)≤B+1` without any independence assertion. Also `S∼Gamma(n,1)` implies `P(S<n/2)≤4/n`. Set `ε=1−x` and `n=max(m,16,ceil(8(B+1)/ε))`; the probability that either `S<n/2` or `S−W_u≥ε n/2` is at most `4/n+2(B+1)/(ε n)≤1/2`. On the complementary event `Z_u>x`. Since `n≤A_u/ε` for a finite explicit `A_u` (for example `m+17+8(B+1)`), this proves the **uniform polynomial lower tail**

```
P(Z_u>x) ≥ c_u (1−x)^B,    0<x<1,
c_u = π (m−1)! / [2 (m−B−1)! A_u^B] > 0.
```

For `x≥1/2`, the ratio of the beta upper cap for `v` to this lower bound is at most `2^{−(v−1−B)}/c_u` when `v−1≥B`. For `x≤1/2`, use `P(Z_u>x)≥P(Z_u>1/2)>0` and `P(Z_v>x)≤P(Z_v>0)→0`. These choices of `v` are uniform in the threshold. Thus for each fixed `u` there is a finite `v₀(u)` such that the claimed strict order holds for all `v≥v₀(u)`. This argument gives no adjacent-label bound as `u` itself grows.

### How to certify a genuine limiting crossing if one is found

There are finite-state certificates for hub probabilities on **both** sides. For a specified `d≥2`, the event “degree has hit `d` by size `n`” gives the lower bound `(1−1/d) P(D_n(u)≥d)`: at its first hit the new leaf supplies a surviving tag population as above. For an initially degree-one vertex the first hit of degree two is also an actual one-leaf core birth. For the upper bound, let `r_n(x)` be the number of leaf neighbours of `x` at size `n`, excluding `u`. At an arbitrary state, a jump in `D(u)` has intensity `Σ_{x∼u}1/D(x)`. For each neighbour with `r_n(x)=k≥4`, the tag-envelope argument gives an upper bound

```
b(k)=min(1,(k+1/2)/(k/2−1)^2+2/k)
```

on any arrow from `x` into `u` after size `n`; set `b(k)=1` otherwise. This is a per-neighbour bound for the *first* such arrow, by charge before the tag-envelope failure plus failure probability. A hub must gain again, and its next gain has an old neighbour as source. Consequently `E[min(1,Σ_{x∼u} b(r_n(x)))]` is a rigorous hub upper bound. On each nonhub path `D(u)` is eventually a finite constant (by the unbounded-degree implication) and each of its finitely many permanent neighbours is a hub (paper Theorem 3.2), so their leaf stocks diverge and this charge tends to zero. Bounded convergence gives a limsup no larger than the true hub probability; the charge always bounds it from above. The lower approximation converges when first `n`, then `d` tends to infinity. These bounds describe an effective comparison with no appeal to a finite degree crossing; they are not expected to be sharp on very small enumerations.

For a positive threshold there is also an easy quantitative check that avoids independent gamma assumptions. On `T_n`, `E[Z_u|T_n] ≥ (L_n(u)−1)/n` by the paper's bounded leaf submartingale, while `D_n(u)/n` is a supermartingale. In the beta domination coupling from `(n,D_n(u))`, let `B` be the beta limit above `Z_u`; then `E(B−Z_u|T_n)≤(D_n(u)−L_n(u)+1)/n`. For `0<x<x+ε<1`, this gives the interval

```
P(B>x+ε | T_n) − (D_n−L_n+1)/(n ε) ≤ P(Z_u>x | T_n) ≤ P(B>x | T_n).
```

For integer beta parameters and rational thresholds the beta tails are finite binomial sums. Such an inequality can certify a strict crossing if the actual intervals separate. The lower bound may be negative and can then be clipped at zero.

## Exact finite checks

From this directory:

```
taskset -c 13 nice -n 10 timeout 30s python3 -uB check.py --size 10
```

The program enumerates increasing labelled parent sequences, using a common exact integer denominator at each size. It independently reconstructs each directed-arrow transition through size seven. The swap fixture checks 40,318 ordered placements. Through size ten (362,880 trees at that size), all finite-degree tails were nonincreasing by label, with exact equality of the first two rows. As a check at size nine, `P(D_9≥2)` drops from `1487379956579/2194698240000` at labels 1 and 2, to `390030836579/1097349120000` at label 3 and `9311761/292626432` at label 8. Finite agreement is evidence only; even a future finite reversal would not on its own settle either limiting question. An unconditional proof needs an averaged comparison of the two **location-specific limiting laws** on the biased-edge sector, balanced against the earlier-newborn sector, or a separate uniformly accurate test of limiting probabilities with a strict crossed gap.
