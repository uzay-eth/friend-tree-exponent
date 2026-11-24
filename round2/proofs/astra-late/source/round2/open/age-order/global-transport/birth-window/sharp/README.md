> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# Q3: a sharper conditional birth-scale window

For every old tree of size `n≥2` and every `0≤(n+2)x≤2^-66`, the older of the next two labels has strictly larger conditional probability of limiting share exceeding `x`. The founder comparison `(2,3)` holds for `0≤3x≤2^-66`. This raises the explicit birth-scale cutoff from `2^-131400` to **`2^-66`**. It does not reach an order-one birth scale or connect to the large-`vx` region.

The improvement uses the exponential moment of the actual monitored time and an exact stock exit. Neither the old signed ledger nor its hub endpoint inputs need to be changed.

## Keep the return clock in the recurrence

Use the [original retry events](../../hostile/README.md), reaching degree 20 with tracked stock or failing. Fix `θ=1/2`. At a last loss with marked degree `d`, use a virtual return sender of degree `k`, starting at the actual sender degree. It moves to `k+1` at rate `b=k−1+1/d`, and its return channel has rate `1/k`. Stop at a return, a different physical gain at the mark, or virtual degree 32.

To realize this comparison, couple each real sender gain with a virtual degree gain and pad with extra proposals. Until exit the real sender degree is at most `k`, so its gain rate is at most `b`. Thin the sender's real return channel to rate `1/k`. Any discarded return or other physical gain also ends the wait successfully. Thus the virtual return has the same lower probability as the frozen recurrence, including when the actual sender grows much more slowly. A cutoff is charged as failure even if a later real send would have succeeded. Before promotion of a favourable leaf use the same construction with `d=1` and its actual initial parent degree (discard immediately at degree 32 or higher).

An upper transform for the duration of a whole return, including failure, is

```
R(32)=1,
R(k) = [1/k + b R(k+1)]/[1/k+b−θ].
```

All denominators are positive at the chosen `θ`. The checker verifies `R(k)≥R(k+1)≥1`. Extra successful exits only shorten the race; optional stopping with the exponential clock and boundary value one gives the bound (or continue the nominal race with its own drivers after an early real exit). The frozen success bound at `θ=0` is computed with terminal payoff **zero** at 32, whereas this time bound has payoff **one** there. The time estimate therefore does not condition on an eventual success.

## Exit from terminal stock instead of waiting to a fixed age

At degree 20, start the tracked birth–death comparison with birth rate `j` and death rate `j/20`. The physical degree never decreases, so the lower comparison can use the real leaf births and extra death proposals. Let `K=96`. This chain, stopped on `0` or `K`, reaches `K` with probability `(1−20^-j)/(1−20^-K)≥1−20^-j`. Thus replacing terminal survival by this finite hit costs no probability before deadline and mass errors are deducted. The physical stock at this hit is at least 96.

For the duration up to either exit, solve the exact equations with boundary values `V_0=V_K=1`:

```
j V_(j+1) + (j/20) V_(j−1) = (j+j/20−θ) V_j.
```

The rational elimination checks positive pivots and `V_j≥1`. Applying the bounded-state generator to `e^(θt)V_j` stopped at either exit, then using Fatou, proves the transform bound. For the founder, use two initial leaves, death rate `j/2`, and `θ=1/4`; the hit probability is at least `3/4`.

Back through the degree-20 retry use `V_j` as the terminal row (rounded upward). For positive stock `j` and degree `d<20`, write `g` for the next-degree continuation at stock `j+1` and `l` for the same-degree continuation at stock `j−1`. The transform upper is

```
max(g,[j g + (j/d)l]/[j+j/d−θ]).
```

The maximum ensures that any additional physical gain ending at `(d+1,j+1)` has nonpositive drift in the upper comparison. At stock zero after a last loss multiply the return transform by the next-degree value at stock one. This allows all clocks and success histories to be dependent: each new bound holds conditional on the entire state at the preceding stop. The initial serial host-credit race has exit intensity at least one, and its possible continuation states are `(2,0),(2,1),(3,1),(3,2)`, so multiply their largest upper by `1/(1−θ)`. The favourable initial return instead uses `R_1(2)` times the degree-two, stock-one upper; this bounds any initial sender degree below cutoff by monotonicity. The computed probability that either monitored path to the stock exit lasts past `T=46` is less than `4.492×10^-7`.

## Convert the stock exit to a birth-scale threshold

At the `K` hit use the [physical tagged-leaf mass bound](../../../../../../../../astra/source/round2/fringe/hub-spectrum.md#2-a-vertexs-limiting-mass-as-seen-at-time-t). Its remaining mass is at least `Gamma(K,1)−Exp(1)` in a coupling, without assuming independence. It exceeds 32 except with probability at most `e^48/2^96+e^-16`, using the gamma negative Laplace transform at one for the first term. Starting at completion size `m≥3`, write `W=S Z` with `S~Gamma(m,1)`. The probability of `S>16m` is at most `(2e^-8)^m≤8e^-24`. This is a union bound with the monitored events, not a factorization at the random hit time. On all the good events,

```
mZ > (32/16)e^-46 > 2^-66.
```

The combined loss against the original successful-path lower is less than `5.72×10^-7`. The old adverse hub upper bounds remain valid. The [frozen perturbation budgets](../README.md#2-spend-the-losses-in-the-signed-ledger) cost at most four times this loss at a retained vertex and six times on a short chain, both below `1/200000`. Additional homogeneous chain sites still give positive credit. The seed and stars retain their old positive surplus. For the founder the total lower-endpoint loss is below `.000083`, versus the centre-minus-leaf margin `127/1220`.

## Why this method stalls before practical thresholds

This proof first constructs a positive-mass event and then deducts **every** successful path after a common deadline, regardless of its stock or mass. Even the nominal initial favourable return from parent degree two supplies a basic obstruction to making that deduction arbitrarily cheap: with combined rate `1/2+2=5/2`, the event that the first exit is the return and is later than `T` has probability `e^(-(5/2)T)/5`. The nominal comparison can continue after that return with success probability bounded below by its strictly positive frozen retry value. This makes the late-success error at least a positive constant times `e^(-(5/2)T)` for that surrogate. Shrinking the numerical error budget therefore forces a logarithmically growing deadline and an exponentially shrinking mass floor. This is a limitation of a common-deadline certificate, not evidence for a negative physical score. Optimizing the time tilt alone does not remove it.

A route beyond this window is to propagate a threshold payoff `F(y e^t)` through each actual holding-time distribution and retain the stock at each exit, rather than discard all paths later than `T`. At moderate `y`, the upper bounds should also discount late adverse sends: the fixed hub uppers presently collect those sends in full. Such a signed threshold comparison needs new endpoint inequalities. In particular the [large-`vx` asymptotics](../../../fixed-threshold/transition/README.md) do not supply a constant comparison down to the window proved here.

## Replay

From the repository root run `python3 -IB round2/open/age-order/global-transport/birth-window/sharp/check.py`. All decisions use rational arithmetic, including the bounds on exponentials; printed decimals are for readability.
