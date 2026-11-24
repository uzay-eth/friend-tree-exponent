> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 3 and 10; Gemini reviewed the base sequential law, with separate scopes for addenda. See the [release index](../../../../../README.md#claim-20) for exact scope and citations.

# Positive inert core limit and joint component-size exponents

An inert planted tree has one parent slot at its root, counted in that root's degree. Each ordinary vertex rings at rate one; a root choice of the slot makes no change, and the slot never sends. All other choices append one ordinary leaf to the chosen vertex. Write `C(S)` for the number of ordinary vertices of degree at least two **including the slot in the root's degree**, and `m(S)` for the ordinary vertex count. In particular the virgin `J` (one root and one real leaf) has `C(J)=1`, `m(J)=2`. We consider all finite planted seeds with at least two ordinary vertices. Let `μ` be the Claim 6 growth exponent; `1/2<μ<1`.

**Result.** For every such seed `S`, in physical continuous time,

```text
e^(−μt) C(S_t) → Z⁰_S ∈ (0,∞)    almost surely and in L²,
e^(−t) m(S_t)  → W⁰_S ∈ (0,∞)    almost surely.
```

Consequently under the limiting uniform-component law `ν_N` of [forest-age](../README.md), with `K_N` the number of surviving original cores and `V_N` the number of ordinary vertices,

```text
(log K_N / log N, log V_N / log N) ⇒ (μU, U),    U ∼ Uniform[0,1].
```

The order is first the large-tree limit at fixed cutoff, then the large cutoff. The same pair arises when cores are measured by degree in the observed component, provided one uses `log(max(1, count))` for that convention.

## 1. The inert mean estimate and martingale

The inputs are the uniform one-cell estimate at test `1` in [spatial-limit, §1](../../../../../../astra/source/round2/fringe/spatial-limit.md#1-start-in-one-old-cell) (from the finite-fringe renewal estimate of [Claim 6, §3, (3)](../../../../../../astra/source/round2/standalone/claim-6.md#3-the-resolvent-and-the-contour-estimate)), and the uniform planted one-leaf gradient of influence/gradient, §§1–3 (`round2/influence/gradient/README.md`, historical reference). We use the gradient for the **physical** inert process with its slot retained; neither a positive residual particle nor a free seed is substituted for it.

Here are details of the mean transfer. Replace the slot by a real adjacent auxiliary vertex bearing `M≥1` additional leaves. At fixed elapsed time `t` the chance of an incoming auxiliary send before `t` is at most `t/(M+1)`, since auxiliary degree never decreases. Until the first such send, the internal laws can be coupled to agree: sending to the auxiliary has the same internal effect as a null choice in the inert process. Give each initial core its own spatial label. For a fixed initial planted core `v`, its future cell stays on the planted side of the cut; the other endpoint is already a core and cannot transmit a cell label by a promotion. Both versions of the cell count have second moment at most `2e^(2t)`, by rate-one per-core promotion domination, so the mean difference at fixed `t` tends to zero by Cauchy–Schwarz. The uniform rate of convergence `O(e^(−(μ−σ)t))` for the normalized mean of each free cell (`σ<μ`) passes to the inert cell: for two times compare their normalized means to the *same* free coefficient with auxiliary stock `M`, then let `M→∞`. This makes the inert normalized means uniformly Cauchy with the same remainder after sending the second time to infinity. Summing over the finitely many original core labels gives a coefficient `Φ⁰` in the normalization of the actual core count and constants independent of the seed:

```text
0 ≤ Φ⁰(S) ≤ B C(S),
|E_S C(S_t) − e^(μt) Φ⁰(S)| ≤ B C(S)e^(σt),
|Φ⁰(S⁺) − Φ⁰(S)| ≤ L₀
```

for every one-leaf extension `S⁺`. In particular `Φ⁰(J)>0`: this can be seen directly from the [spatial immigration formula, (3)](../../../../../../astra/source/round2/fringe/spatial-limit.md#1-start-in-one-old-cell) in the auxiliary construction. The virgin-cell root has promotion rate `1/2` while both its own clock and its leaf's clock have remained silent and no auxiliary arrow has arrived. That event has probability at least `e^(−3u)` up to time `u` (use an auxiliary dominating Poisson proposal clock of rate at most one). The positive direct immigration term in the coefficient on this event is `d₂ H(J_{1/2})/2` times `e^(−μu)`, with `d₂>0`, `H(J_{1/2})>0` by Claim 6. Integrating gives the uniform lower bound `d₂ H(J_{1/2})/[2(μ+3)]`, which survives the auxiliary limit. This uses only the coefficient for a virgin, not a free-seed almost-sure limit.

The Markov property, the mean remainder and finite expected population on finite intervals imply `E_S Φ⁰(S_t)=e^(μt)Φ⁰(S)`: apply the mean estimate after time `t` and send the remaining horizon to infinity. Thus `M_t=e^(−μt)Φ⁰(S_t)` is a nonnegative martingale. Every changing jump is one leaf addition, of total rate

```text
λ(S) = m(S) − 1/D_root(S).
```

The omitted root-to-slot arrow is the only null choice. Using the uniform jump gradient and `E_S m(S_u) ≤ m(S)e^u`, the stopped jump isometry and removal of the size stop give

```text
E_S (M_t−M_0)^2 ≤ L₀² m(S) ∫₀ᵗ e^((1−2μ)u) du
                 ≤ L₀² m(S)/(2μ−1).
```

In detail, between changing jumps the finite-variation derivative of `M_t` cancels the predictable mean of its jumps by the martingale identity. Localize at bounded sizes, compensate the squared genuine jumps, and pass to the unstopped process by the integrable bound on their total expectation. In particular `M_t→M_∞` almost surely and in `L²`, with mean `Φ⁰(S)`. Cauchy–Schwarz gives the seed-dependent lower bound `P_S(M_∞>0) ≥ Φ⁰(S)^2/[Φ⁰(S)^2+L₀² m(S)/(2μ−1)]` whenever `Φ⁰(S)>0`. This bound by itself need not be uniform over seeds (ordinary size can be large at fixed core count). We next establish a different uniform bound.

For the identification with the count, use the **finite-time** form of the same planted gradient. It gives `|E_{S⁺} C(S_h)−E_S C(S_h)|≤B₁e^(μh)` (the singleton cost has exponent `.525<μ`, and the branch costs are bounded by `e^(μh)` as in the gradient proof). The backward conditional mean martingale of `C(S_t)` therefore has jumps bounded in absolute value by `B₁e^(μ(t−u))` at time `u`. For deterministic `0≤r≤t`, the backward martingale isometry on `[r,t]` implies

```text
E_S [e^(−μt)(C(S_t)−E[C(S_t)|F_r])]^2
  ≤ B₁² m(S) ∫ᵣᵗ e^((1−2μ)u) du.
```

The backward equation and stopped compensation are justified on each finite horizon by the size Yule bound and the finite-time gradient. With `r=t/2`, the right side is exponentially summable on each fixed positive mesh. To control the remainder of the conditional mean, `e^(−μr) C(S_r)` is bounded by `e^(ε r)` eventually almost surely for every `ε>0`: the mean bound gives `E C(S_j)≤B' C(S)e^(μj)` at integer times, so use Markov–Borel–Cantelli there and monotonicity between integers. Choose `ε<μ−σ`. Then

```text
e^(−μt) E[C(S_t)|F_r] = M_r + O(e^(−(μ−σ)(t−r)) e^(−μr)C(S_r)) → M_∞.
```

Chebyshev and Borel–Cantelli give convergence along each mesh `1/k`. Monotonicity of `C(S_t)` gives the lower and upper interpolation factors `e^(−μ/k)` and `e^(μ/k)`; send `k→∞`. For the `L²` conclusion, the backward bound with `r=0` and the mean estimate yield `E C(S_r)^2≤B₂(C(S)^2+m(S)) e^(2μr)` at any deterministic `r`. The conditional mean remainder now tends to zero also in `L²` with `r=t/2`, and the centered error has the bound above. This proves `Z⁰_S=M_∞` before using strict positivity.

## 2. A uniform chance to establish a positive inert branch

Let `p=P_J(Z⁰_J>0)>0` by the virgin coefficient and the `L²` argument above. We will find a deterministic lower bound `δ>0` for `P_S(Z⁰_S>0)` over **all** finite planted seeds with at least two vertices. No claim that the entire tree has a uniformly positive limit *size* is needed.

First, consider the stock `L` of ordinary leaf neighbours of any core `x`, with `D=L+c` including the slot where appropriate. Stock gains have predictable intensity at least `L` from its leaves, and stock losses (promotions) have intensity `L/D≤1`. These descriptions hold in the full physical history; root-to-slot choices are null. Start with stock `k≥3` and let `ζ` be the first time it drops below three. For `f(l)=1/l`, at `l≥3` the conditional drift before stopping is at most

```text
−1/(l+1) + (l/D)/(l(l−1)) ≤ −1/(l+1) + 1/((l+1)(l−1))
                            ≤ −1/[2(l+1)] ≤ −3/(8l).
```

Here `D≥l+1` because `x` has a parent slot, real or inert. Stopped Dynkin (first at finite sizes and time) gives

```text
P(ζ<∞) ≤ 2/k,           E ∫₀^ζ ds / L_s ≤ 8/(3k),
```

At crossing `L_ζ=2`, and the stopped reciprocal is nonnegative; monotone convergence gives the integral. The estimates hold conditionally from any realized stock `k` at a stopping time, uniformly over the entire rest of the tree.

If a new child of `x` has just been promoted into a virgin and `L_x≥k`, couple its future with a fresh inert virgin using the internal arrows until the first arrow from `x` into that child. The child's sends to `x` are genuine in the full tree and null in the comparison, with identical effects *within* the branch. The comparison's marginal law remains that of `J` from the promotion history; extend it with its own physical drivers after separation. Up to stock stop `ζ`, the rate of the incoming arrow is `1/D_x≤1/L_x`, so by compensation the probability of a separating arrow before `ζ` is at most `8/(3k)`. Including the event of a stock stop bounds the probability of **any** separation by `14/(3k)`. All these bounds concern the joint construction, without conditioning the inert continuation on no incoming arrows. The union bound, not a product of dependent probabilities, now gives probability at least `p−14/(3k)` of a positive virgin limit and no separation. On that event the planted tree contains every core of the tracked branch at all subsequent times; its positive limit is at least the branch limit times the fixed elapsed-time factor `e^(−μτ)` from the promotion time `τ`.

To produce this configuration uniformly, choose a deepest core `x` in the orientation away from the root. Its only nonleaf slot is its parent (inert if `x` is the root); it has at least one ordinary leaf. Choose an existing leaf `ℓ` of `x`. During the first unit require no rings of `x` and at least `k` rings of `ℓ`. The leaf stays a leaf in this unit, each of its rings adds one leaf at `x`, and no exterior action can promote a leaf of `x`. Thus `L_x≥k+1`, with exactly one nonleaf slot still. This event costs `e^(−1)P(Poisson(1)≥k)`, uniformly, by the independent clocks of these two vertices. In the second unit require a first ring of `x` and that its first choice be a leaf. Until that choice the stock is nondecreasing, so its chance is at least `1/2`; the time criterion has probability `1−e^(−1)`. This gives an event of probability at least

```text
a_k = e^(−1) P(Poisson(1)≥k) (1−e^(−1))/2 > 0.
```

The successful event is known at the selected promotion time, no later than time two. The promotion leaves `L_x≥k` and opens a virgin with its one true exterior parent edge. Fix `k` large enough that `14/(3k)≤p/2`. The conditional construction from that stopping time proves `P_S(Z⁰_S>0) ≥ a_k p/2 =: δ>0` for every seed. The uniform probability concerns the future limit, not a uniform convergence time or a bound independent of the virgin input `p`.

For an arbitrary initial seed, at deterministic integer `j`, positivity of the limit is unchanged by the factor `e^(−μj)` relating the original limit to that of the restarted state. Thus `P_S(Z⁰_S>0 | F_j)≥δ` by the Markov property. These conditional probabilities converge to `1_{Z⁰_S>0}` because the event is measurable in the increasing union of the histories. This proves strict positivity almost surely, without importing a free-seed logarithmic lower bound. In particular `log C(S_t)/t→μ` almost surely.

## 3. Ordinary size and the high-cutoff joint law

In a state of size `m`, the holding time until the next **change** is exponential of rate `m−1/D_root`, not of rate `m` and not in general a function of `m` alone. Given the current state, that wait is independent of the winning genuine arrow. One may construct the whole chain with independent unit exponentials `E_m`, independent of the randomness used to pick genuine successors, so that the elapsed time to reach size `n` is

```text
τ_n = Σ_{m=m₀}^{n−1} E_m/(m−1/D_root(S at size m)).
```

The actual waits are not independent of the chain of shapes. However, `D_root≥2` and the difference of a reciprocal holding rate from `1/m` is nonnegative and bounded by `1/m²` for `m≥2`. The expected sum of the clock corrections is finite, hence the sum is finite almost surely. Also `Σ (E_m−1)/m` converges almost surely by its summable variances, and the deterministic harmonic sum minus `log n` converges. Thus `τ_n−log n` has a finite real limit. The waiting gaps tend to zero almost surely (use `E_m/(m−1/2)` and exponential tails), so `e^(−t)m(S_t)→W⁰_S∈(0,∞)` between epochs as well. No Gamma law or independence of `W⁰_S` from the shape trajectory is needed. In particular, if desired along the genuine embedded size chain, `C(S at size n)/n^μ→Z⁰_S/(W⁰_S)^μ` almost surely.

Now use [forest-age, (1)](../README.md) as the total-variation transfer for full ordinary rooted shapes. In the comparison mixture, put the virgin at independent uniform age `u=U log N` on a common sample path with `U∼Uniform[0,1]`. For `U>0`, both log rates apply almost surely there, giving the stated joint limit for the slot-reconstructed count `C` and ordinary size `m`. The statistic reconstructs one extra root degree for the inert slot. In an actual cutoff component any surviving nonroot original core can lose core status in the ordinary observation only if it has a neighbour deleted by the cutoff; the top may additionally lose its former parent's edge. More generally, total variation on *unmarked shapes* alone need not preserve original labels. For the count of **original** cores, use the proof of (1), specifically [the discarded-branch error on the unit ball](../README.md#the-branch-after-a-promotion): its coupling simultaneously rejects any further deletion within the branch, at the same uniform error. On each successful coupling the only removed adjacency from a surviving vertex is the top's parent edge, and restoring one slot recovers original core status exactly. Thus the TV transfer proof applies verbatim with `K_N` retained as an observable, on the same sequential limit. Counting ordinary observed cores removes at most one on the good coupling; `log(max(1,C−1))−log C` is bounded in absolute value by `log 2`, so that convention also has the stated joint limit. In particular the limiting pair shares the **same** uniform variable, not two independent ages.
