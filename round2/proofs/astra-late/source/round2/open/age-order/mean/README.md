> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# Adjacent mean limiting shares

Write `m_v=E Z_v`, where `Z_v=lim D_v(n)/n`. The seed labels are exchangeable. This note gives strict ordering from label 2 through 8 and proves strict ordering between consecutive sufficiently late labels. There is no effective eventual cutoff here.

## Continuous-clock mass and terminal bounds

From a size-`n` tree let `M(T,v)=lim_(t→∞) e^(−t)D_v(t)`. Its mean `U(T,v)` is `n E[Z_v|T]`: at each jump the waiting time depends only on the size and is independent of the donor/target choices, so the limiting Yule size factor has mean `n` and is independent of the entire discrete path. In particular this step does not assert independence between the mass `M` and the Yule factor. For a core, the same scaled limit is obtained using its leaf count: its core degree grows only at promotions of its own leaves, at total rate at most one.

For a core with `l` leaves, `U≥max(0,l−1)` by integrating the discounted leaf drift: the excess births from core donors are nonnegative and the loss intensity is at most one. If `l≥1`, retain a single tagged leaf as a lower stock, and stop on its first birth or removal. It gives birth at rate 1, and can be killed for the lower comparison at constant rate `1/2`, since the host has degree at least two. Thus it gives two leaves, at time with discounted probability `∫exp(−(1+1+1/2)t)dt=2/5`; from those two the previous lower bound is at least 1. This gives the additional lower floor `2/5` whenever `l>0`.

For the upper bound, old core donors contribute at most `Σ_w 1/D_w(0)` to discounted births, by nondecreasing degrees. Each future core neighbour is created by promoting one of the tag's leaves, at some time `s`. It can send discounted future births at rate at most `1/2`, so its entire charge is at most `e^(−s)/2`; the same promotion costs one discounted leaf `e^(−s)`. Pairing these terms gives `U≤l+Σ_w 1/D_w(0)`; this is at most the ordinary degree upper bound `D_v(0)` (each donor contributes at most one). All intensity integrals are integrable: new core neighbours accrue at rate at most one, hence their expected count by time `t` is at most `t`.

For a leaf with parent of degree `p≥2`, suppress the parent clock until its independent first ring. Jensen gives expected reciprocal degree at age `t` at least `e^(−t)/p`, as the suppressed degree's mean is bounded by `p e^t`. Selecting the leaf at that ring yields a newly promoted core with one leaf, worth at least `2/5` before discount. Therefore `U≥(2/5)∫e^(−3t)dt/p=2/(15p)`. For the upper bound, while the tag is unpromoted the promotion intensity is at most `1/p`; its discounted selection probability is at most `1/(p+1)` by thinning. After selection the upper value is at most `1+1/p` by the core bound, so `U≤1/p`. To improve this, stop at the first parent degree increase or promotion. Until then promotion has constant intensity `1/p`, and the rate-one clock of the tag always increases the parent's degree. One can sample those two clocks independently using directed-channel thinning; other parent degree gains may preempt both. If the degree increases first, its continuation value is at most `1/(p+1)`. Ignoring the other degree gains can only increase the bound: at a preempting time the continuation of the two-clock comparison is at least `1/(p+1)`. Indeed that comparison has discounted value

```
[(1+1/p)/p + 1/(p+1)] / (2+1/p) = 1/p − 1/[(p+1)(2p+1)] ≥ 1/(p+1).
```

## Exact signed comparison

Use `r_u=Σ_(x~u)1/D_x`, `δ_u=1/[D_u(D_u+1)]`, and `g(S,v)=U(S,v)/|S|`. Define `I(u)` as the difference of `g` at middle and tip after extending old vertex `u` by a path of length two; define `H(u,v)` as the difference between the two marked tips after attaching one at each old vertex. Conditional on the original size-`n` tree, the difference of future mean shares between its next two births is exactly

```
[ Σ_u r_u/(D_u+1) I(u) + Σ_{unordered u,v edge} (r_v δ_v − r_u δ_u) H(u,v) ] / [n(n+1)].
```

For distinct old targets, the second rate is `r_v−δ_u` if adjacent and `r_v` otherwise. Reversing the order cancels nonedges; repeated targets make the marks exchangeable. The only unmatched history is second attachment to the first newborn. In the program, each signed terminal difference is bounded jointly in its final tree, then propagated by the original rational next-target law for up to three jumps. The propagated interval is intersected with the all-time interval at the current state. Unmarked pre-placement shapes are grouped with exact probabilities.

Replay from `.`:

```sh
timeout 45 taskset -c 13 nice -n 10 /usr/bin/time -f 'elapsed=%es' python3 -Bu round2/open/age-order/mean/check.py
```

The replay passed 99 direct-versus-paired finite-share comparisons and all six strict ordering checks (8.95 s elapsed). The exact lower endpoints, in order, were:

| Gap | Rational lower bound |
|---|---:|
| 2–3 | `16267847/119750400` |
| 3–4 | `29786243/3353011200` |
| 4–5 | `7673996221/3138418483200` |
| 5–6 | `3132747937103/2824576634880000` |
| 6–7 | `295170272247530747/560207699251200000000` |
| 7–8 | `7334665005036749572441/26554517193745981440000000` |

## Eventual ordering

The normalised mass `U` extends locally to the infinite marked-core environments as the limit of expected scaled leaf stock. For a tag initially core with `c` core neighbours, the expected absolute drift of the stock after the unit exponential rescaling is bounded by `e^(−t)(1+(c+t)/2)`, since the tag promotes neighbours at rate at most one and permanent core donors each have rate at most `1/2`. For a tag initially a leaf, stock is zero before promotion and can gain its first leaf at rate at most one; its initial core degree at promotion is one. These bounds give uniformly vanishing integrated tails near any fixed marked input. The finite-time expectation is local and continuous by the old-cell cutoff: couple the fixed-time histories, then pass from bounded observations to stock using domination by a Yule degree process of the fixed initial degree, with a finite second moment. Thus `U` is continuous. The same construction and bounds apply after the finitely many insertions at a fixed rooted input.

Use `U` instead of `g` in the one-birth score of [the late hub law](../late/README.md). The conditional flux is `n(n+1) E[Z_(n+1)|T_n]=Σ_o F_M(T_n,o)`. The core-initiation term is bounded by 1, while the `L_o` leaf initiations give a bounded contribution because post-placement parent degree is `D_o+1` and a newborn's `U≤const/(D_o+1)`. The birth score is continuous, since `L_o` and the finite list of choices are fixed in a local neighbourhood. Its limiting stationary mean is positive: the leaf lower bound above, weighted by the attachment rates, and the `C/16` counting inequality in the Fable Claim 3 audit give `Σ_o F_M≥C/120`. The core local law and uniform integrability then give

```
m_(n+1) ~ b n^(μ−2),    b=E[X] π(F_M)>0.
```

For two births use `U` in the signed comparison, so the share denominator has the extra factor `n+2`. Group the initial-vertex score over each core and its adjacent leaves as in [the adjacent hub law](../adjacent/README.md). The marked values `U` are bounded by their starting degrees (at most two), and `r_u≤D_u`, so the per-core score is bounded and continuous. The local law and uniform integrability again imply that `n^(3−μ)(m_(n+1)−m_(n+2))` tends to a finite `a`. Since `m_n→0`, summing these differences gives `n^(2−μ) m_(n+1)→a/(2−μ)`, identifying `a=(2−μ)b>0`. This proves eventual strict adjacent ordering.
