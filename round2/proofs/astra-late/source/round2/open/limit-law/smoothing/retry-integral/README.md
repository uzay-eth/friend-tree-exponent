> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../../README.md#claim-17) for exact scope and citations.

# Q4: integrate the crossing cost over adapted setup times

For the full edge-seeded density `f`, this gives

```
||f(·+h) − f||₁ ≤ C |h|^(747/11494),          |h|≤1,
f ∈ W^{s,1}(R)                              for 0<s<747/11494.
```

The cutoff increases from [`13/278`](../planted/README.md) (`0.04676259…`) to `747/11494` (`0.06499043…`). The gain is `14561/798833` (`0.01822784…`). The main change is to integrate the *actual* crossing bound against finite-age stock distributions over the setup times of all reached retries. A second change is to keep waiting through harmless endpoint rings during setup. Local boundedness remains open; §4 gives a necessary weighted fold-potential bound for it.

## 1. Integrate over reached setups, without a stock threshold

Start the physical clock at the three-vertex tree. At every reached attempt start, use the [two-anchor race](../compact/README.md#1-a-physical-edge-with-large-stock-quickly-and-uniformly) to set up adjacent cores and then start the [two independent lower-stock families](../../density/README.md#some-edge-receives-only-finitely-many-crossing-arrows), each from one. Restart on a setup failure or the first observed extinction of either stock. The setup race is modified in §2. Write `R_j` for attempt start and `U_j` for its successful finite preparation time, setting either to infinity if it is not reached. Both are stopping times on the graphical space (including the lower-family proposals). All the conditional bounds concern the history observed at these times.

Suppose that uniformly at reached starts the defective exponential moment of the increment to the next start is at most `ρ<1` at some `0<r<1/2`, and the conditional moment of the finite preparation time *relative to that start*, counting only completed preparations, is at most `M<∞`. Then

```
Σ_j E[e^(r U_j); U_j<∞] ≤ M / (1−ρ).
```

Indeed the reached-start exponential moments are bounded by `ρ^(j−1)` by conditioning at the previous actual restart. Infinitely many failures have probability zero, and any reached attempt completes or fails its preparation in finite time. Let `S` be the preparation time of the never-failed attempt. This witness need not be a stopping time. Since it is one of the `U_j`, `P(S>t)≤C e^(−rt)`.

At deterministic time `t`, if any edge joins two cores with positive stocks of incident real leaves, select an edge maximizing the smaller stock (break ties by birth labels). Call the maximum `B_t` and set it to zero if there is no such edge. By the [first-send barrier](../../../age-order/early/README.md#a-no-send-bound-without-an-age-cutoff), the conditional chance of any subsequent crossing of the selected edge is at most `4/(1+B_t)` on `B_t≥1`, uniformly in the physical exterior. Thus the no-cross sublaw rejects at most `E[min(1,4/(1+B_t))]`, including the histories with no eligible edge. The selection itself uses only the ordinary physical time-`t` history.

On `S≤t`, the never-failed attempt has two lower stocks positive at age `t−S`; the minimum of these is at most `B_t`. Both reciprocals are nonnegative. Bound the expectation by summing their terms over **every reached setup** with `U_j≤t`, counting `1/(1+Z_i)` only when its lower stock `Z_i` is positive. For each setup use its unconditional lower process law after conditioning at `U_j`; there is no conditioning on being the never-failed attempt. For this marginal estimate one may continue a discarded stock with fresh drivers after its trial terminates, keeping the same linear birth–death law. Domination of physical stock is needed only on the never-failed trial, where nothing was restarted.

At age `u`, conditional on being positive a lower stock is geometric on `{1,2,…}` with parameter `q=e^(−u/2)/(2−e^(−u/2))` by the [scalar law](../mixture/README.md#1-a-usable-edge-at-an-exponential-rate). For `0<q<1`, the geometric law has `E[1/Z | Z>0]=−q log(q)/(1−q)`, with limit one at `q=1`. Thus

```
E[(1+Z_u)^(−1); Z_u>0] ≤ C (1+u)e^(−u/2) ≤ C_r e^(−ru).
```

Condition on each finite reached `U_j`, plug in `u=t−U_j`, and sum. The resulting missing mass is bounded by

```
P(S>t) + C_r e^(−rt) Σ_j E[e^(r U_j); U_j<∞] ≤ C e^(−rt).
```

This integration uses the entire preparation-time distribution and charges each reciprocal at its own age. It does not multiply probabilities of no-cross events at different deterministic times. Indeed, even for one fixed edge, the events of a *future* crossing at two successive times are nested; their intersection is the later event, so a product bound for such crossings is in general false.

## 2. Faster adapted setup, with an exact moment certificate

At the current endpoint `x` and chosen incident leaf, keep the initial rate-two race. If `x` wins, process its ring, declare failure and restart immediately. If the leaf wins, there are at least two leaves at `x` and at most one core neighbour. Now wait for the first ring of `x` that **picks a leaf**. Rings of `x` choosing its core neighbour add a leaf at the neighbour and do not consume any leaves at `x` or increase its core-neighbour count. Other vertices' rings cannot do either before this promotion. At each successive rate-one ring of `x`, the conditional chance of choosing a leaf is at least `2/3`. Equivalently, order the leaf channels first in the fresh uniform mark: the first mark in `[0,2/3]` certifies a promotion if none has occurred yet. This stochastically bounds the wait by `Exp(2/3)` conditional at the leaf win. At promotion `x` retains a true leaf and the promoted core has one; start the two lower stocks there with fresh conditional drivers.

On this preparation branch the next failure is the first extinction of the two stocks. Use the [defective tail](../iteration/README.md#2-actual-retry-cost-and-full-law-powers) `Σ_(j≥1) (j+1) 2^(−j−2) exp(−ju/2)` for its time on failure, with total defective probability `3/4`. For `0<r<1/2`, a rational upper envelope for the defective exponential moment is

```
F_r = 3/4 + r [ Σ_(j=1)^8 w_j/(j/2−r) + (3/4−Σ_(j=1)^8 w_j)/(9/2−r) ],
w_j = (j+1)/2^(j+2).
```

The total failure moment now has the uniform upper bound

```
ρ_r = [1 + ((2/3)/(2/3−r)) F_r] / (2−r).
```

The first term in brackets pays for the initial endpoint win; the other uses the conditional wait bound followed by the fresh defective extinction moment at promotion. These multiply by conditioning, including at real failure times, rather than by treating physical future masks as independent. The preparation moment alone has a finite upper bound `M=((2/3)/(2/3−r))/(2−r)`. At `r=747/10000`, exact rational arithmetic gives `ρ_r<99991/100000<1`. Even without the faster preparation, keeping the old `13/200` clock and just integrating yields `13/226` for the Sobolev cutoff.

## 3. Transfer to the full law

Let `ν_t` be the actual sublaw retained by the deterministic selector in §1 and by absence of future crossings. The only terminal mask is the [shape likelihood `L`](../README.md#1-the-holding-weights-cancel-exactly) for this selected edge; the setup retries are used in the discarded-mass bound only. Each side contains its endpoint and at least one leaf. For fixed embedded paths and time-`t` physical history, set `S_a=a+b` for the two split coefficients and `N=|T_t|`. The [frozen translation estimate](../planted/README.md#1-a-derivative-kernel-in-the-physical-split-coordinate) bounds the conditional total variation modulus by `C √N √(|h|/S_a)`. Its angular argument also gives the same bound with `|h|` replaced by the volume of a Borel test set for small-set probability, using the [global angular inverse estimate](../../../../../../../../reviews/cold/astra/checks/limit-iteration-sol.md#separate-viable-refinement-use-the-two-coefficients-jointly).

Keep `L` in the average. The exact reference-to-physical identity and `X≤S_a` on the retained sublaw bound the averaged coefficient by `E[√N X^(−1/2); retained]≤C e^(t/2)`, using `E X^(−1)<∞` and `E N=3e^t`. By §1 the actual discarded mass is `O(e^(−rt))`. Hence, with `δ=|h|`,

```
||τ_h f−f||₁ ≤ C [e^(−rt) + e^(t/2)√δ].
```

Balance at `t=log(1/δ)/(1+2r)`, giving the claimed modulus at `σ=r/(1+2r)=747/11494`. The same balance for Borel sets gives `P(X∈A)≤C |A|^σ` for `|A|≤1`, hence weak `L^{11494/10747}` and strong `L^p` below it. Integrating the translation estimate gives the stated `W^{s,1}` range; interpolation with the strong `L^p` estimates gives `W^{s,p}` for `1<p<11494/10747` and `0<s<1/p−10747/11494`. The fractional endpoint is not asserted.

## 4. A necessary fold-potential condition for local boundedness

Fix a deterministic cut rule and time. Write `η` for its history/path parameter and `dλ(η)` for its mixing measure **including** `L`. Restrict to the parameters with `a,b>0`. At the interior maximum `m(η)`, the Beta image folds, with nonzero input density and nonzero quadratic curvature. The local inverse-branch expansion gives measurable `c(η)>0` and `d(η)>0` such that the fixed-path density is at least `c(η)/√(m(η)−x)` throughout `m(η)−d(η)<x<m(η)`. For example, use half its positive one-sided asymptotic coefficient and take the first dyadic radius on which this bound holds, testing the continuous local branch density at rationals. The radii may depend on the parameters; no uniform curvature or Beta weight at the fold is assumed.

If the full density is essentially bounded by `H` on an open positive interval `I`, its retained subdensity is dominated by it. Tonelli then bounds the following potential by `H` for almost every `x` in `I`:

```
∫ c(η) 1{0<m(η)−x<d(η)} / √(m(η)−x) dλ(η).
```

The integrand for fixed `η`, taken as zero outside the indicated open interval, is nonnegative and lower semicontinuous in `x`, so Fatou makes the integral lower semicontinuous. The bound therefore holds at **every** `x` in `I`, and with the same `H` at every deterministic cut used to test the law. In particular on any restricted part with `c≥c₀>0` and `d≥d₀>0`, the weighted dyadic shell masses above each `x` must be summable against the weights `2^(k/2)/√d₀` (shell widths `d₀ 2^(−k)`, open outer boundary). A mere `1/2`-power interval estimate for the weighted locations is insufficient: a location density proportional to `(m−x₀)^(−1/2)` just above `x₀` has that estimate uniformly in the location, yet its potential at `x₀` diverges logarithmically. This is a necessary condition on genuine weighted fold locations, not a claim that the physical mixing law saturates it. Controlling the rejection mass by itself does not verify it.

## Replay

`nice -n 10 timeout 12s python3 -IB round2/open/limit-law/smoothing/retry-integral/check.py` checks the exact defective-clock envelope and rates, the finite-age geometric reciprocal bound on rational fixtures, and the dyadic critical-scale obstruction. The stopping-time, measure and Sobolev arguments are analytic.
