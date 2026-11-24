> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 3 and 10; Gemini reviewed the base sequential law, with separate scopes for addenda. See the [release index](../../../../../README.md#claim-20) for exact scope and citations.

# Fixed component intensities: young shapes and old singletons

This follows up the fixed-shape question in the [Neon-2 discussion](../../../../../../../reviews/cold/astra/attacks/neon2-forest/primary/AUDIT.md), which identified why the leading total-variation approximation at scale `N^-mu log N` does not give fixed-shape intensities at scale `N^-mu`. The limit here is sequential: first `n→∞` at each fixed cutoff `N`, then `N→∞`. The observation counts each component with an original core once, rooted at its highest surviving core. Deletion uses **original** degrees.

Let `M_N(R)` be its limiting count per physical `C_n` for a fixed finite ordinary rooted shape `R`. Use the stable-age [parent birth law](../../../../../../astra/source/round2/degree-law.md) `β`, and write `W_o=lim_t e^-t D_o(t)` for the focal parent's mass from its own promotion. Put `A=E_β W_o^mu`. At each promotion time `τ` of the focal parent a new physical child `v_τ` starts as a virgin with one leaf; the exterior continues to send to it. Let `T^0` be the inert, planted virgin dynamics, with its one exterior slot counted in the degree and omitted from the observed shape. Then:

* For **each** fixed `|R|≥2`, $`\displaystyle \lim_{N\to\infty}N^\mu M_N(R)=A\,I(R),\qquad I(R)=\int_0^\infty \Pr(T^0_u=R)\,du.`$ (1)
  This is precisely the inert occupation coefficient at the fine scale. The two-vertex shape has coefficient `2A/3`. Every rooted tree with at least two vertices can be reached by a finite sequence of positive-rate additions, so all coefficients in (1) are positive. The finite coefficient is the inert hit probability times `1/(|R|−1/D_root)`; the hit probability is a finite rational recursion from the virgin through trees of size at most `|R|`.
* The **singleton core** instead has coefficient $`\displaystyle \lim_{N\to\infty}N^\mu M_N(\bullet) = E_\beta\sum_{\tau<\infty}\mathbf1\{D_{v_\tau}(\infty)<\infty\} \left(\min\{W_o,\inf_{w}W_w\}\right)^\mu.`$ (2)
  The infimum ranges over the final nonparent neighbours of the bounded child; it is over a nonempty finite set. `W_w=lim_t e^-t D_w(t)` uses the **parent's** time origin, not the neighbour's own birth time. This expectation is finite and strictly positive. The inert singleton intensity is zero. Formula (2) is a physical lifetime coefficient with live input, not an occupation of the inert virgin or a product of independent hub events.

By the same proof one may keep the relative parent-degree mark. For a bounded continuous test `f` of `log(D_parent/N)`, the coefficient in (1) becomes `A I(R) ∫_0^∞ f(e) mu e^-mu e de`. In (2), let `m_*` be the minimum in the summand and replace `m_*^mu` by `W_o^mu ∫_{log(W_o/m_*)}^∞ f(e) mu e^-mu e de` on `m_*>0` (zero otherwise). Thus the other cut neighbours can impose a positive lower shift on the singleton's log overshoot. For an unrooted observation, sum over its finitely many distinct top-rootings, without a multiplicity for equivalent roots.

## 1. Stable-age formula and a full-count bound

Put `L=log N`, `t=L+y`. The [top-edge formula](../README.md) gives, setting the integrand to zero when `t<0`,

```math
N^\mu M_N(R)=\int_{-L}^\infty \mu e^{-\mu y}
 E_\beta\!\left[\mathbf1_{D_o(t)>N}\sum_{\tau\le t}
 \mathbf1_{D_{v_\tau}(t)\le N,\,F_N(v_\tau;t)=R}\right]dy. \qquad\text{(3)}
```

The focal parent's birth state has no nonparent core neighbour: it was a leaf until promoted. Its one inherited exterior core neighbour remains live. **All** child branches in (3) start at actual promotions of the focal parent, including when it promotes its inherited original leaf. Inside an older child branch we cannot discard any of its previously existing neighbours, nor can we assume that a late cut preserves its virgin leaf.

Write `s=|R|`. Call a counted child *intact* if no immediate edge anywhere in its downward component leads to a removed descendant; otherwise call it *cut inside*. The intact branch then has ordinary size exactly `s`, omitting only the parent edge. There is a crucial bound valid even without knowing its full past. At age `u=t−τ` the intact event requires at most `s−2` genuine internal births in the child's branch. Independently of outside interactions, their predictable intensity is at least `3/2` at all times: at size `m≥2` it is exactly `m−1/D_top`, with `D_top≥2`, and incoming parent arrows only increase ordinary size. A constant-rate thinning gives the bound `b_s(u)=Pr(Pois(3u/2)≤s−2)`, or zero if `s=1`. Compensating promotions, whose rate is at most one, bounds the expected number of intact contributors with `t−τ>h` by `∫_h^∞ b_s(u)du`, uniformly in the exterior.

The **live holding law** at an already attained branch shape of size `s` can also be written exactly. Freeze all genuine changes inside the branch and suppress the single parent-to-top input; retain exterior evolution, in particular top-to-parent sends. If `D_p^fr(a)` is the parent's degree in that frozen system, then the conditional survival probability of the holding time for another `u` is

```
exp(−(s−1/D_top)u) · E_frozen exp(−∫_0^u da / D_p^fr(a)).
```

Use separate directed clock channels until the first genuine change. The suppressed internal channels have constant total rate before that change; the missing incoming channel uses a fresh integrated-hazard threshold, without conditioning the frozen exterior on its eventual success. In particular for a virgin the inert holding mean is `2/3`; the live mean is between `1/(3/2+1/D_p(0))` and `2/3`, since `D_p^fr` is nondecreasing. The asymptotic coefficient in (1) uses the vanishing input of *late* promotions. The frozen formula does not replace the original dynamics of early children in (2).

For a child with an internal cut, some vertex in its retained component has a high descendant neighbour. Choose one such boundary vertex `w` at a depth at most `s−1` below the top. Its own original degree at time `t` is at most `N`. More importantly, at most `s−1` of its current real leaf children can be ordinary leaves (they all survive). The number of its core children is at most its personal rate-one clock count since it became core. We will use the stronger fact that a high descendant had to be born early in **parent clock**, regardless of this leaf-stock observation. Every such high vertex is a core at time `t`, and belongs to the downward core tree of the focal parent. Its ancestry begins with the component top, a child promoted at `τ`. At the moment a lower descendant promotes its degree is two. From age `t-u` its degree is dominated by a Yule process from two, with second moment at most `6e^{2(t-u)}`. The expected number of all core descendants at time `u` grows at most at rate one per core, so their birth intensity is at most `e^u` in expectation. Consequently the expected number of high descendants promoted after parent age `b` satisfies

```math
 E_\beta\#\{v\text{ in the focal descendant core tree}: \tau_v>b,\ D_v(t)>N\}
   \le 6 e^{2y}\int_b^t e^{-u}du\le 6e^{2y-b}. \qquad\text{(4)}
```

Choose one removed descendant per cut-inside component. Distinct counted top children have disjoint branches, so this is an injection. In particular (4) controls **all** cut-inside components whose tops promote after `b`, at the fixed-shape scale without a factor `log N`. This also explains why a fixed-size component that looks young after deletion must not simply be assigned a recent promotion time.

## 2. Integrable envelopes, including the high positive-age tail

For `y≤0` the cut-inside count is at most the number of high proper descendants, with mean at most `6e^{2y}` by (4) with `b=0`. For `y>0`, refine this as follows. A removed descendant with ancestor first child promoted after `2y` contributes to (4) with `b=2y`, costing at most six; if `2y>t` no such descendant exists. The number of the other component tops is at most the focal parent's promotion count up to `min(t,2y)`, whose expectation is at most `2y`. This gives an integrable bound after multiplying by `e^-mu y` on the right.

For intact counts use the integrable age bound in §1 on the right. On the left, retain the high-parent indicator by Cauchy–Schwarz. On each unit age strip `t−τ∈[j,j+1)`, the qualifying count is bounded by the focal rate-one clocks on that strip, whose fourth moment is uniformly bounded. Its expected qualifying count is at most `sup_{u∈[j,j+1]} b_s(u)`. If `Z_j` is the qualifying count then `P(Z_j>0)≤E Z_j`, and Cauchy–Schwarz gives `E Z_j²≤(E Pois(1)^4)^{1/2}(E Z_j)^{1/2}`. The square roots of these second moments sum over `j` by exponential decay. Thus the entire intact count has uniformly bounded second moment. The second-moment degree bound `E D_o(t)^2≤6e^{2t}` now yields a bound `C_s e^y` for the left envelope. Both tails in (3) are integrable for `0<mu<1`, uniformly in `N`. These estimates also apply with any bounded parent mark inserted.

## 3. Recent promotions: proof of (1)

At each fixed `y`, first restrict intact components to child age `≤h`. All these promotions occur at parent ages `τ≥t−h`, tending to infinity. Choose a fixed early history time `b`, and replace the high-parent event by `E_b={e^-b D_o(b)>e^-y}`. For all but countably many `y`, first `N→∞` then `b→∞` makes the replacement error tend to zero: the high indicators converge to `1{W_o>e^-y}`, and the recent parent clock count has bounded second moment at fixed `h`. On `E_b`, the parent has degree `>e^{b-y}` at every later promotion; the conditional risk of an incoming arrow to the child before observation is at most `h e^{y-b}`. In the inert coupling the entire branch remains of size `≤N` at these bounded ages except with probability at most `2 e^h/N`; in the live branch use the uniform size intensity `m+1/2`, giving the corresponding `5e^h/(2N)` bound. The event used for conditioning is measurable before these promotions. Conditional expected future observations are taken in the full, unconditioned marginal of each continuation; failure of the coupling is paid separately.

The focal promotion rate is `1−c_o(r)/D_o(r)`. To replace it by one on the recent interval, one **cannot** simply use `c_o(r)≤1+Pois(r)` together with a *fixed* degree floor at `b` and then take `r→∞`. On paths with `W_o>e^-y` the ratio tends to zero almost surely because `D_o(r)∼W_o e^r` and `c_o(r)=O(r)`. On the discrepancy event where `E_b` holds but `W_o≤e^-y`, use the trivial bound one and pay `h P(E_b,W_o≤e^-y)`, which vanishes as `b→∞`. Thus compensation gives the recent coefficient `Pr(W_o>e^-y)∫_0^h Pr(T^0_u=R)du` in the limit. Send `h→∞` using the uniform age tail in §1. For a log-parent mark one can use the early value `y+log(e^-b D_o(b))` on `E_b`, with the same bounded-clock and dominated-continuity argument. Finally integrate `mu e^-mu y` by §2; the identity `∫ mu e^-mu y Pr(W_o>e^-y)dy=A` gives (1), and changing variables `e=y+log W_o` gives the marked version.

## 4. Early promotions: why only a singleton survives at this scale

Fix an early age cutoff `b`. There are finitely many early children almost surely, with count dominated by a `Pois(b)` clock. On the event of subexponential time growth of personal core degree (each core degree is bounded by one plus its own rate-one clock since promotion), the following physical dichotomy applies simultaneously to the countably many vertices and edges. A vertex of bounded degree eventually has a fixed finite neighbour set. It cannot have a bounded neighbour, since an edge between two eventually fixed degrees would continue to deliver arrows at a positive rate. Every unbounded vertex has positive continuous mass almost surely, by the tagged-family argument: at each first hit of a new high degree `d`, the fresh leaf can seed a birth–death family with survival `1−1/d`; survival gives diverging stock and positive mass by the gamma-minus-loss comparison at growing stock. Letting `d→∞` excludes unbounded degree with zero mass. These are the [physical tagged-family and continuity inputs](../../age-order/scaled/README.md); they apply also in the locally finite live environment under `β`.

If a fixed early child is unbounded but has positive mass, then at times `t=L+y` it is either above threshold, or its stock of ordinary leaves grows exponentially. In the latter case all of those leaves remain attached to it in its low component, so its component cannot have fixed finite size. This does not require a limit for its high indicator at equality. A bounded child has finitely many final neighbours, **all** of positive mass. At observation, if any of them other than its parent survives the cut, that neighbour brings a diverging stock of ordinary leaves, again ruling out a fixed finite size. If all are cut, the child is a surviving singleton. Its parent is included in the high requirement. For all `y` off the countable set of positive fixed-threshold atoms of the relevant masses (or for almost every `y` by Fubini along trajectories), its singleton indicator converges to

```
1{D_child(infinity)<infinity} · 1{min(W_o, W_w : w final nonparent neighbour)>e^(-y)}.
```

Dominated convergence applies to the finitely many early children by the focal clock bound. The remaining cut-inside count for children born after `b` costs at most `6e^{2y-b}` by (4), regardless of their age at observation or any omitted older neighbourhood. The limiting early-child sum is increasing in `b`; take `b→∞` at fixed `y`, then use the uniform envelopes of §2 for the `y` integral. This gives zero contribution to (1) and exactly (2) by Tonelli. In particular a surviving singleton has *no inert-only prediction* even in the sequential law at this finer scale. Component sizes `≥2` come from the observation-age layer `u=O(1)`, whereas singletons come from bounded **parent-clock ages** of the child's promotion, in the sense of the successive truncations above. The focal parent itself is observed at age `log N+O(1)` in the integrable-window sense of §2.

**Strict positivity of (2).** Use the first focal promotion to create a virgin child `v` with leaf `w`. Arrange that the focal parent has a large leaf stock before this promotion, then that `v` promotes `w` before an incoming arrow to `v`. This makes `v` a bare degree-two core with neighbours `o,w`; the promotion does not increase `D_v`. Arrange a large leaf stock at `w` while preventing any gain to `v`. These finite operations have positive conditional probability: silence the relevant core source clocks over bounded intervals except for the two prescribed promotions; each available real leaf can ring enough times to grow its own core's stock. Incoming arrows to a bare `v` can only be from `o,w`, so no unmentioned external arrow can spoil the freeze; external additions to those sources beyond the specified scripts do not decrease their leaf stock while they are silent. Before the first focal promotion the virgin parent's original leaf supplies the initial stock growth. The directed first-send barrier `min(1,2/(j+1))` with zero receiving-leaf immigration applies to each of the two arrows into the bare core, from stocks `j` as large as desired. The two sources each have a high-probability positive mass event as well, by the tagged gamma-minus-loss bound at large stock. A union bound gives positive conditional probability of no gain ever to `v` and positive masses at both neighbours; no independence of those events is required. On this event the corresponding summand in (2) is strictly positive.

## Exact finite checks

`python3 -B round2/open/forest-age/fixed-shapes/check.py` returned:

```text
84 positive exact rooted-shape occupations at sizes 2–7; pair 2/3
5 exact old-core singleton cutoffs; inert singleton impossible
```

The occupation recursion uses actual directed rates and subtracts the null slot selection from the holding rate. The singleton fixture checks deletion by original degrees on a bare core between two large-stock neighbours. It does not estimate the lifetime expectation in (2).
