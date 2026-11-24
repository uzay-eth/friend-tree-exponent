> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../../README.md#claim-17) for exact scope and citations.

# Quantitative deterministic-cut mixture

The full edge-seeded law has a density in **`L^p` for `1<p<1/(1−ξ)`**, where one explicit admissible `ξ` is given below (in particular `ξ=1/80000` works). For Borel sets `A` of Lebesgue measure at most one,

```
P(X∈A) ≤ C |A|^ξ.
```

This uses only two sides, so it does not require uniform boundedness of the conditional multi-side density. Do not select the last crossing retrospectively. At a deterministic integer time select the first eligible edge using the current physical tree, bound the probability of any subsequent crossing, and apply the fixed-cut sublaw on the rest.

## 1. A usable edge at an exponential rate

Use the two-anchor finite setup of [the multicut proof](../README.md#4-arbitrarily-large-permanent-cutsets-and-the-summation-obstruction): in a fixed integer time `d`, with conditional probability at least `q>0` from *every* tree with a core, prepare adjacent cores each with a designated private leaf. One explicit setup starts at a core-tree endpoint, silences it for a unit while a leaf rings at least twice, then requires its first ring in the next unit to select a leaf. Retain one spare leaf at the original core and the new core's first leaf by silencing both anchors after the promotion until integer time `d=2`. Incoming sends cannot remove these leaves. The endpoint and the designations can be chosen by birth labels from the history.

For this script take `q=(3/4)e^(−3)(1−2/e)(1−1/e)`: silence the old anchor for the first unit, use the leaf's Poisson tail, then pay for a first ring within one unit and for choosing one of at least three leaves with at most one core neighbour. Conditional on that promotion, silence both anchors until the end of the second unit at cost at least `e^(−2)`.

At the prepared edge, start the two embedded leaf populations with birth and death rates per member `1,1/2`, using the [independent family construction](../../density/README.md). These are lower bounds on distinct leaf stocks. If either family goes extinct, end the failed attempt at the next integer (allow equality), and retry from the physical state there. If setup fails, retry immediately at its deterministic end. Families from a failed trial need not be used again: the fresh future clocks at the restart give the same conditional bounds. All restart times are stopping times on the enlarged graphical space; the unrestarted successful attempt is used only to prove an eligibility probability, never to condition the cut.

For a family from one, let `Z_u` be the linear birth–death count, `T_0` its extinction time, and `z=e^(−u/2)`. Then

```
P(T_0≤u) = (1−z)/(2−z),      P(u<T_0<∞)=z/[2(2−z)],
P(1≤Z_u<K) ≤ K z.
```

The first formula follows from the scalar backward birth–death equation; given `Z_u>0` the count is geometric on positive integers with parameter `z/(2−z)`. These facts do not condition the physical exterior on survival. The probability of an attempt that never fails is at least `q/4`. The **defective** waiting time to a failed attempt's restart has uniformly exponential tail, since beyond setup only the first family extinction can cause a restart and the two finite-extinction tails sum to at most `e^(−u/2)`. In particular for some `r>0, ρ<1`, uniformly over restart histories,

```
E[e^(rD) ; attempt fails | history] ≤ ρ,
```

where `D` is the duration on failure. For `0<r≤1/4` the increase from `r=0` is at most `r(d+1+4)e^((d+1)/4)` by splitting the exponential tail integral at `d+1`. Take `r=q/(56e^(3/4))` and `ρ=1−q/8`. Induction multiplies these defective bounds. Infinitely many failures have probability zero, and the start time of the first never-failed setup (including its setup duration) has an exponential tail: its exponential moment is at most `e^(rd)/(1−ρ)` by summing over the number of failures. This argument does not assert that this start time is a stopping time.

At integer `t`, except with probability `O(e^(−rt/2))` a never-failed pair started by `t/2`. For every attempt with successful setup at integer `s≤t/2` use the unconditional bound `P(1≤Z_(t−s)<K | history at s) ≤ K e^(−t/4)` on each branch; these bounds are summed over at most `t+1` possible times, without conditioning on which trial survives. With `K=ceil(e^(t/8))`, the chance of having no core edge at time `t` with at least `K` real leaves at both ends is therefore `O(e^(−c t))` for some `c>0`. Select any such edge measurably from the time-`t` tree. The bound applies to existence, not to the subsequent law of a particular retrospectively successful family.

The [first-send barrier](../../../age-order/early/README.md) gives chance at most `2/(K+1)` of each future cross direction, conditional on the actual selected state and edge. It uses only the sender's other ordinary leaves; core targets supply no assumed return immigration. Put `c=r/2`. The no-cross retained sublaw has missing mass at most `B e^(−ct)` with `B=8e^(2r)/q+40`: the never-failed setup time costs at most `8e^(2r)e^(−ct)/q`, the small-stock union is bounded by `4(t+1)e^(−t/8)`, and the future cross costs at most `4e^(−t/8)`. Use `c≤1/16` and `sup_(u≥0) u e^(−u/16)=16/e` to absorb the latter two. There is **no size truncation**. The future mask is never factored as independent of the clock law except for the fixed-cut terminal likelihood in the [Dirichlet identity](../README.md#2-exact-dirichlet-image-including-the-core-limits).

## 2. Uniform estimate on the retained sublaw

For a fixed time-`t` state choose the edge by the rule above, which is adapted. On its no-cross sublaw the ratio is `g(x)=a x^µ+b(1−x)^µ`, with `x~Beta(n_1,n_2)` independent of the embedded paths under the reference construction. Here `a,b≥0` and `a+b>0` on the weighted law. The reference density is bounded by `n_1+n_2−1`: it is that factor times a binomial point probability (integer sizes). For any `h>0`, the exceptional derivative region `{|g'|≤h}` has length at most `2h/[µ(1−µ)(a+b)]` since `−g''≥µ(1−µ)(a+b)` on `(0,1)`. Elsewhere there are at most two monotone inverse branches with derivative at most `1/h`. Thus for arbitrary Borel `A`, not only intervals,

```
Leb(g^(-1)(A)) ≤ 2h/[µ(1−µ)(a+b)] + 2|A|/h,
inf over h > 0 = 4 sqrt(|A|/[µ(1−µ)(a+b)]).
```

For the part of `A` above `v>0` we may impose `a+b≥v` since `g≤a+b` pointwise. With terminal likelihood `0≤L≤1`, the conditional retained sublaw is bounded by `4(|T_t|−1) sqrt(|A|/[µ(1−µ)v])`. Average **in the original physical law of the time-`t` history**, multiplying by the history-measurable eligibility indicator and dropping that indicator only for the upper bound. The result is `12 e^t sqrt(|A|/[µ(1−µ)v])`, since `E|T_t|=3e^t` from `T₃` by the [free Yule size clock](../../../../../../../astra/source/round2/fringe/almost-sure.md#3-second-moments-and-positivity). There is no cardinality factor or conditional size expectation given no future crossing. The [edge-seeded left-tail upper bound](../../left/README.md#matching-left-tail-order-and-the-negative-moment-threshold) pays `P(X≤v)≤D v^(1/µ)`; enlarge `D` to cover all `0<v≤1`. The edge-seeded embedded chain passes deterministically through shape `T₃`, so using the continuous clock from this shape does not change the law of `X`.

For `δ=|A|∈(0,1)`, balance all three errors with

```
ξ = c/[2(1+c)+µc],      t = floor((ξ/c) log(1/δ)),      v = δ^(µξ).
```

The three exponents are equal because `1/2−ξ/c−µξ/2=ξ`. One may take `C=B e^c+12/sqrt(µ(1−µ))+D` in the set bound. Numerically `c≈2.63051e−5`; since `0<µ<1`, the displayed exponent exceeds `1/80000`. For comparison, even reoptimizing the version with an `e^(2t)` size cap in the kernel gives the weaker exponent `c/[2(2+c)+µc]` at the same missing-mass rate. Null sets follow by outer approximation. This is a single integer-time choice per tolerance; equivalently the deterministic integer cut mesh yields sublaws with exponentially small loss and a resolution-dependent density estimate.

Writing `f` for the resulting full density and applying the set bound to `{f>H}` (whose measure is at most `1/H`), one gets `Leb{f>H}≤(C/H)^(1/(1−ξ))` for `H≥1`. Integrating layer cake proves `f∈L^p` for every `1<p<1/(1−ξ)` (the sublevel contribution is bounded by `∫ f`). Each retained sublaw is dominated by the full law and hence has density bounded above by `f` almost everywhere. In particular these deterministic-mesh approximants are uniformly bounded in such `L^p`, and even converge to `f` in `L^p`: split where `f≤H` and where `f>H`, use vanishing discarded mass on the first set and domination by `f^p` on the second. This does not assert a globally bounded or continuous full-law density.

## Conditioning note: a cut parent's degree

For the normalized degree **fraction** of a history-selected cut endpoint, work on its fixed-cut no-cross sublaw. Within side `i` let `b_i` be the limit of `D/N_i` along side-size epochs. On `L>0` it exists: on the physical sublaw `e^(−u)D(u)` converges by the [physical degree-supermartingale bound](../../../age-order/scaled/README.md#1-positive-continuous-time-mass-has-no-atoms-even-in-an-infinite-exterior), and under the reference clock `e^(−u)N_i(u)` has a positive limit; the resulting size-epoch ratio is a function of the embedded shape sequence alone. The whole-tree fraction is `b_i P_i`. Conditional on shapes with `b_i>0`, the unmasked pushforward has no atoms. A future selector can be included as a mask by domination (also across a countable list of candidate cuts); it cannot be treated as preserving conditional Beta clocks. This argument does not remove a zero atom, nor deduce atomlessness at random thresholds depending on the same holds. For the [fixed-shape forest parent mark](../../../forest-age/fixed-shapes/README.md), the normalized limiting positive mark is atomless for a different reason: integrate over the continuous stable-age coordinate first, with masses and the singleton gate held fixed. The [simultaneous threshold comparison](../../../forest-age/fixed-shapes/diagonal/README.md) likewise removes boundary equalities under age integration; no independence after conditioning on eventual hub selection is needed.
