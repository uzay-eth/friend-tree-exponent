> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Stationary upper prices for the two realized D18 extremes

At common multiplier `θ=1.2415`, the actual *four-bare-degree-two P4 block* contributes between `−4.314e−7` and zero per core to the averaged signed contrast. The *bare degree-17 host with sixteen virgin children* contributes between zero and `3.830e−18` per core (sum of all eight apex-owned classes there). These bounds remove the two concrete counterexamples from the pointwise-comparison argument. They do not price the whole complement; no exponent endpoint moves yet.

## Spine: use the true entry flux and price the block before taking a hinge

Use the exact physical charged rows from [the all-eight extractor](../../fringe/README.md) at the same donor rates in the two witnesses. The signed P4 row is `q(d,p)=a+b p+c/d` when all four participating nonfounders are bare degree-two cores; `d` is the live arm degree just below them and `p` the live parent reciprocal. There is one such path per apex. Starting at its birth, the focal can enter that face only by **four consecutive tip promotions**: a gain anywhere in the growing chain would give a future participating vertex degree above two, irreversibly. At size `n=2,…,5` the desired rate is `1/2`, the intrinsic total rate is `n−1/2`, and exterior input only raises this total. Therefore the discounted probability of entry at `s≥LO` is at most `∏_(n=2,…,5) [(.5)/(LO+n−.5)] < .000479205`.

At entry the arm is always virgin. On the face a gain at any of the four protected cores kills the block, exactly at `3+1/d+p`. The remaining tracked arm gain and promotion rates are `g∈[l+.5,l+(k+1)/2]` and `b_arm=l/d`. For `q`, `−q`, or either pointwise hinge, `check.py` bounds the discounted occupation from above using

```
 (s+3+1/d+p+g+b_arm) V(d,k)
    >= reward(d,p) + g V(d+1,k) + b_arm V(d,k+1).
```

For signed rewards test both discounts, each shared exterior corner, and both gain-rate corners. For nonnegative rewards use discount `LO`; a hinge kink is retained too. Above `d=16` use a constant terminal equal to the maximum nonnegative reward on the true reciprocal interval divided by `LO+3`, rounded up. This covers unbounded live-arm growth rather than killing on an omitted arm event. The block upper bound is `HI × entry_upper × max(0,V(virgin))`; capping the entry product at zero on negative `V` does not use an upper entry frequency with the wrong sign. When `V(virgin)<0`, the **entire block** cannot contribute positively, by the stopped Poisson correction. This is a bound on the signed mean of the block, not the expectation of the pointwise positive hinge.

| Common multiplier | Raw pointwise hinge caps (positive / negative) | Signed block caps (positive / negative) |
|---|---:|---:|
| `1` | `1.104e−8 / 3.421e−7` | `0 / 3.416e−7` |
| `1.2415` | `7.269e−9 / 4.316e−7` | `0 / 4.314e−7` |

This uses the physical stable-age identity on the focal apex for the upper average; it does not require a sampled founder-free tree or a universal positive incidence of the face.

## Hubs: bound arity, simultaneous young-arm multiplicity, and the actual bare host

For a true physical newborn let `N` count its core-child births and `V` the number of still-virgin such children at age `x`. Births have predictable rate at most one. Every virgin permanently exits at rate `3/2+1/D_host` (one leaf gain, its own promotion, and the actual host send). Attach independent rate-`3/2` private clocks to the baseline losses at every birth and bound births by a rate-one proposal process, keeping the actual birth as a subset chosen predictably. Thus one marked immigration–death process simultaneously dominates `N` and `V`, not an independence hypothesis conditional on the physical host path. In particular, for any `j,K≥0`, the Poisson factorial identity at age `x` gives

```
 E[binom(V,j) 1(N≥K) | age=x]
 <= (1-e^(−3x/2))^j / [j! (3/2)^j] × P(Poisson(x)≥max(0,K−j)).
```

Integrating at the independent exponential age of rate `LO` gives the exact finite alternating sum in the replay; the dominating coupled observable is increasing in age in distribution, so its expectation increases when the exponential age rate is decreased. One direct construction for this monotone-age assertion extends the comparison process backwards in time from zero stock, sampling at the same terminal time and admitting all extra earlier births. This avoids multiplying marginal arity and multiplicity estimates. At `K=16`, the first/second/third binomial-moment bounds are `8.006e−4, 4.178e−4, 1.454e−4`. These are reusable with degree/arity-dependent row envelopes; the third one counts the actual number of triples, not merely the probability of some triplet.

The exact whole-tree balance `π(k)=1` improves the arity-only upper bound: for `1≤j≤K`, time to `K−j` further births stochastically exceeds the corresponding sum of rate-one holding times, even with adaptive rates. Thus `π(k≥K)≤(1+LO)^−(K−j) π(k≥j)`; summing the resulting lower bounds on the earlier tails gives `π(k≥16)<.000434179`. A separate robust Poisson solve on actual `(D,k)` using `g=L+ρ`, `b=L/D`, and `0≤ρ≤(k+1)/2`, with degree-overflow value zero, gives `π(D=17,k=16)<1.067e−13` with **arbitrary arms**. A bare host of that exact degree cannot be approximated by arity alone.

For the counterexample host with precisely sixteen virgins, the stronger bound races the next birth (rate ≤1) against the first loss (rate ≥`3j/2` when `j` virgins have been born). Loss cannot be undone; otherwise the target cannot occur at arity sixteen. At entry the visit intensity is at least `24` from virgin exits, and it cannot recur. Its age-weighted mass is at most

```
 HI × Π_(j=0,…,15) [1/(LO+1+3j/2)] / (LO+24) < 1.448e−18.
```

The actual total apex contrast at this host is affine in its true parent reciprocal, with extrema `2.425665` (rounded up) and `2.645585` (rounded up) at `θ=1.2415`. Its upper product is `3.830e−18`. This includes every charge *owned by that apex*; it does not assume every nearby apex charge vanishes. The `P4` spine face and the exact hub class used here are disjoint as charged apex states.

## What is enough for a material endpoint gain

Let `E=πδ−−θπδ+`, `P=πΦ−`, `Q=πΦ+`. Use `E=(μ−a)P−θ(b−μ)Q`, certified denominator ceilings, the GMP denominator floors, and the prior independent lower slack payments (`4.410155511e−6`, `3.921982939e−6`). If `μ≤l`, then

```
 E ≤ (l−a) P_hi − θ max((b−l) Q_lo, existing_upper_slack_floor).
```

If `μ≥r`, use the reverse inequality with `max((r−a)P_lo, existing_lower_slack_floor)` and `Q_hi`. For a **`1e−5`** inward move from each endpoint of the high-parent bracket, after the two explicit exception budgets above, the allowed contrast of the *entire remaining* charge at `θ=1.2415` is bounded below by `−1.2736295e−5` for the lower move or above by `+2.9319160e−6` for the upper move (strict inequalities with these rounded-in targets suffice). These are targets for complementary signed bounds; they are not currently attained. With only nonnegativity of the remaining opposing charges the shared all-motif sum cap is about `1.39414e−4`. Both old endpoints remain feasible in that coarse accounting relaxation even with these exception caps and the independent slack floors. Those are feasible **accounts**, not physical tree laws.

The next high-leverage check is therefore a signed upper bound on the **positive hinge of the complement** at one common multiplier, not increasing the depth of the exact sixteen-virgin fixture. For high arity, contract the true mixed/tail and higher-motif coefficients against joint arm binomial bounds (or against the physical receiver flow), with a separate bounded stock/arity block for populated lower degrees. Keep the contrasting tail corrections explicit: the denominator cap and the nearly proportional small rows do not bound the signs of the remaining rows. The all-eight extractor can supply actual row operands; no full complementary envelope was replayed here.

## Exact replay

```sh
taskset -c 13 nice -n 10 /usr/bin/time -f 'cpu=%U+%S wall=%e' timeout 8s \
  ./round2/.venv/bin/python -Bu \
  round2/occupation-neon/ratio-transfer/averaged/check.py b3f18_lo.bin b3f18_hi.bin
```

Passed (`0.39+0.02` CPU seconds, `0.90` wall). Exact fractions check each rounded Poisson inequality and the algebraic endpoint targets; the stdout point decimals are descriptive. Inputs are the pinned existing payloads, audited production comparisons (`round2/attacks/neon1/outputs/production-final/oracle.log`, historical reference), [denominator ceilings](../../denom/README.md), [first-moment slack floors](../../fringe/README.md) and [high-parent bracket](../../high-parent/README.md). This replay does not repeat the payload scans or claim an independent review of its own derivation.
