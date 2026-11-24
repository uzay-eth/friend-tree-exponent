> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../README.md#claim-26) for exact scope and citations.

# Full stock payoff at small birth scale

For every old tree at size `n≥2`, the older of the next two labels has strictly larger conditional share tail for

```
0 ≤ y=(n+2)x ≤ 1/422500000000000000 ≈ 2.367×10^-18.
```

The founder comparison also holds on `0≤3x≤1/422500000000000000`. This improves the previous explicit cutoff `2^-66` by a factor greater than 174. It is not a practical-scale window and does not reach the separate `.10` seam.

## The entire Poisson payoff at each successful exit

Keep the frozen serial and favourable-parallel retry rules through degree 20. From degree 20 use the same lower linear stock comparison as in the [deadline proof](../sharp/README.md), now stopping at zero or at `K=16`; if stock already exceeds `K`, stop immediately and use just `K` leaves. The chance of reaching `K` from stock `j<K` is `(1−20^-j)/(1−20^-K)`, at least the frozen terminal credit `1−20^-j`. Thus the probability of the credited exit is at least the corresponding frozen serial-middle or favourable-parallel input. Let `τ` denote elapsed time on this path and `A` its stock-exit event. Failed paths pay zero.

At the exit with actual size `s`, apply the [full `H(s,K;x)` payoff](../../../hub-ranking/poisson-tail/README.md). Given the history there, start the fresh physical clocks. The remaining amplitude is `S'Z`, with `S'~Gamma(s,1)` independent of the embedded continuation that defines `Z`. Multiplying the beta mixture by this independent gamma converts it into the law of `X` with a Poisson(1) index `N`: if `N<K` use `Gamma(K−N,1)`, otherwise use zero. In particular the **conditional remaining amplitude** stochastically dominates `X` at every exit. This conversion holds conditional on the exit history; it makes no assertion that the initial size amplitude is independent of the exit or of `τ`.

For threshold `z`, the full conditional loss of the payoff is at most `P(X≤z)`. For a light replay integrate it using

```
P(X≤z) ≤ α + I sqrt(z),
α = P(N≥K),    I ≥ E[X^(-1/2); X>0].
```

The checker uses the factorial weights and outward bounds for `e^-1`; at gamma shape one it bounds `sqrt(pi)` by `9/5`, and for shape `r≥2` it uses Cauchy–Schwarz and the exact inverse first moment `1/(r−1)`. It verifies `α<1.909×10^-14` and `I<.268529`. At the initial common completion size `m`, the size amplitude satisfies `P(S>16m)≤8e^-24` for every `m≥3`. On `S≤16m`, `W>16y` suffices for `Z>x`. Apply the exit payoff at `z=16y e^τ` and union-bound the size exception once. The resulting loss against the frozen endpoint is bounded by

```
α + 4 I sqrt(y) E[e^(τ/2); A] + 8 e^-24.
```

This charges no common deadline and keeps the actual threshold shift even on late successful paths. For the founder use the two-stock comparison with death rate `j/2` and the `1/4` time tilt. Cauchy–Schwarz gives `E[X^-1/4;X>0]≤sqrt(I)≤1`; the checker uses one for this inverse moment and `2^-13` for `(16y)^(1/4)`.

## Return and stock transforms on successful paths

At terminal stock solve the finite tilted exit equations with values **zero at stock zero and one at `K`**, dividing out the earlier elapsed time. With death rate `q=j/d`, use

```
a_0=0,   a_j=j/(j+q−θ−q a_(j−1)),   V_j=a_j V_(j+1).
```

For a last-loss return use the virtual sender from the [deadline proof](../sharp/README.md#keep-the-return-clock-in-the-recurrence), now with `R_32=0` and

```
R_k=max(1,R_(k+1), [1/k+(k−1+1/d)R_(k+1)]/[1/k+k−1+1/d−θ]).
```

The floor one is essential: an extra physical return can turn a virtual failure into a successful early exit, so it need not decrease the success-weighted payoff. With that floor it has nonpositive drift against the upper. The other maximum keeps the row decreasing in sender degree, and even a slower sender is covered by the coupled virtual growth with padded proposals. Positive-stock rows use the same directed drift check as in the deadline proof, multiplied by the relevant succeeding transform. For the initial serial host race, the tracked birth and old-parent send total `1+1/k`, the tracked loss rate is `1/2`, and a host gain preempts to the base `(2,1)` value. Maximizing the linear-fractional bound over the host growth rate is bounded by the larger of its zero-growth value and that base value. The same expression depends linearly-fractionally on `1/k`, so test `k=2` and its limit at infinity. The favourable initial race instead uses its actual parent degree, with zero cost if the lower already vanishes at the cutoff.

With `θ=1/2` the uniform serial transform is less than `466.261`; the uniform favourable-parallel transform is less than `1487.026`. These are separate budgets. Let their payoff losses be `η_s,η_e` from the preceding formula. At a retained core the serial weight summed over itself and attached old leaves is `w+l/(2d)≤3/2`. The total edge weight charged there is also at most `w+l/(2d)`: the old-leaf edges use `l(s+t)`, and adverse direct charges owned there use at most `(d−l)s`. Chain boundary edges are charged only to the chain. A short chain through length three uses serial weight at most one and at most four edge weights of size at most `1/3`, since both endpoints of each are cores. The budget per retained core is thus `3(η_s+η_e)/2`, and per short chain `η_s+4η_e/3`. At the stated cutoff the exact replay verifies losses below `4.843×10^-6` and `4.048×10^-6`, respectively, both below the frozen `1/200000` surplus. Longer chains still add positive homogeneous serial credit; seed and stars have ample surplus.

## Obstruction and seam

The limiting degree-two terminal with one direct high-degree core neighbour has frozen local slack only `5.47868233×10^-6`. The direct adverse coefficient does not vanish there (`s=1/6`). Cashing out the Poisson payoff at stock one is especially costly: at a degree-20 stock-one exit its zero-threshold value is `1/e` against the frozen `19/20` survival credit. Even ignoring return time, the limiting stock-one Poisson payoff at `y=.01` loses `e^-1(1−e^(-.01))` from its own zero-threshold value. At the terminal's combined serial weight `1/2`, this illustrative loss exceeds `.00182`, more than 300 times the frozen slack. This is a comparison-payoff loss, not a physical upper or a counterexample. The stock-16 exit avoids the atom loss, but its inverse moment times the tilted return transform is still too large for the available debit budget at `.01`.

A practical bridge needs threshold-dependent adverse credits or new local transfers at near-tight profiles, in addition to retaining the full payoff rather than bounding its lower tail by a single inverse moment. The present certificate does not overlap the separate work targeting `y≥.10`.

## Replay

`python3 -IB round2/open/age-order/global-transport/birth-window/payoff/check.py` checks the frozen margins read-only, the successful-path transform, rational payoff and size errors, and the explicit terminal controls.
