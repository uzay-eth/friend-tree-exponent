> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../README.md#claim-17) for exact scope and citations.

# The left-tail constant from large-star departure

Let `X` be the fixed-size Claim 6 limit from an edge and put `α=1/μ`. In the physical normalization of the actual core count, let `V` be the positive continuous-time core limit `lim e^(−μt) C_t` of the inert one-slot virgin (a root with a parent slot and one real leaf). Let `V_0,V_1,…` be independent copies, independent of the rate-one Poisson epochs `s_1<s_2<…`. Then

```text
Y_∞ = V_0 + Σ_(i≥1) e^(−μ s_i) V_i,
P(X≤x) ~ 2 E[Y_∞^(−1/μ)] x^(1/μ)       as x↓0.
```

The constant is finite and strictly positive. The dependencies are the [uniform spatial mean remainder](../../../../../../astra/source/round2/fringe/spatial-limit.md#1-start-in-one-old-cell), the [positive inert virgin limit](../../forest-age/inert/README.md), and the degree-sensitive free gradient and large-core downward-passage estimate of the [left-tail proof](../left/README.md). In particular the negative-moment argument below is conditional on the degree-sensitive comparison, not on a hypothesized independence after departure.

## 1. Conditional law at the first departure

Let `K≥4` be the first size with two cores. In a star at size `j≥3`, the only way to leave is to use the hub as source, of probability `1/j`. Thus `P(K>j)=2/(j−1)` and `P(K=k)=2/[(k−2)(k−1)]`. Given the entire embedded history through departure at `k`, the unlabelled departure seed is always `S_k`: hub of degree `k−2`, joined to a degree-two promoted root with one leaf, with `k−3` other leaves on the hub. Let `X_k` be the continuation's limit with denominator the actual total size to the power `μ`, and let `Y_k` be its continuous-time core limit with time zero at departure. The Markov property gives the full conditional law through this seed, not just its conditional expectation. In the free tree the total size rate is its size; sampling the wait and the directed-arrow winner independently at each jump gives

```text
Y_k = W_k^μ X_k,     W_k ∼ Gamma(k,1),     W_k independent of X_k.
```

Here `W_k=lim e^(−t)N_t` by the free Yule clock. Thus `k^μ X_k` and `Y_k` have the same weak limit if either converges, since `W_k/k→1` in probability. No independence from the elapsed physical time before departure is needed.

## 2. The hub limit on bounded intervals

Decompose the physical tree into the hub, its ordinary leaf stock, and the disjoint planted branches rooted at its core neighbours. Initially there is one virgin branch. At a hub ring, choosing an ordinary leaf creates a new virgin; choosing a core neighbour adds a leaf to that branch root. In either direction the parent slot remains in that root's degree. A branch-root send into the hub creates a hub leaf but has no internal effect in the sending branch, exactly like a root-to-slot null choice. A hub leaf sending to the hub also creates another hub leaf. Neither kind of addition can spoil the coupling below.

Use the same rate-one hub Poisson clock in the physical tree and an ideal system that creates an independent inert virgin at every hub ring. Both start with the common initial virgin. At a successful immigration share fresh independent branch drivers, with common internal arrows and the same root-to-slot marks. Stop agreement at the first hub ring selecting a core neighbour in the physical tree. Prior to this stop the number of such neighbours is `1+P_{t−}`, where `P` counts hub rings, and the hub degree is at least `k−2`. Compensation, or summing conditional failure probabilities at the rings, gives

```text
P(failure by t) ≤ (t+t²/2)/(k−2).
```

The ideal system is run even on failure; its marginal law remains an independent Poisson immigration process with independent inert evolutions. The physical ordinary leaf stock need not converge. The physical core count on `[0,t]` agrees with the ideal count `C_∞(u)=1+C_0(u)+Σ_(s_i≤u) C_i(u−s_i)` up to that error. This includes the hub's one core and omits the artificial slots from the vertex counts.

## 3. Passing to the infinite-time amplitude

A bounded-time coupling by itself does not justify passage to the limiting eigen-amplitude. Here there is a uniform approximation from the two-core seeds. For each finite free seed with `C(S)=2`, write `M_t=e^(−μt)Φ(S_t)` where `Φ` is in core-count normalization. The spatial mean gives `|E_S C_t−e^(μt)Φ(S)|≤b C(S)e^(σt)` and `E_S C_t≤b C(S)e^(μt)` (sum the uniform bound over the initial cells). For the martingale tail, the [rate-weighted squared jump bound](../left/README.md#2-conditional-variance-and-a-uniform-negative-moment-start) is at most a constant times `C(S_t)`. Hence martingale isometry and the same mean bound yield

```text
E_S |Y_S−M_t|² ≤ b' e^(−μt).
```

The spatial remainder does not directly compare `Φ(S_t)` with `C_t` on the realized state. Instead approximate the martingale **with a bounded-time forecast**, which avoids needing a degree-sensitive gradient for the centered finite-time counts.

For fixed `t,h>0`, put `F_{t,h}^k=E[e^(−μ(t+h))C_(t+h)^k | state at t]`. Its `L¹` distance to `M_t^k` is at most `2b² e^(−(μ−σ)h)`. Couple the physical and ideal systems using their joint filtration at time `t`, which reveals only elapsed proposal clocks, then continue using the same construction. Given this filtration the conditional forecast for each marginal is its own Markov forecast: after agreement failure use separate fresh drivers, and before failure use fresh independent uniforms and the shared proposal streams. Thus the expected absolute difference of forecasts is at most

```text
e^(−μ(t+h)) E[ |C_(t+h)^k−C_∞(t+h)| ; failure by t+h ].
```

For fixed `t,h` this goes to zero by Cauchy–Schwarz. In both systems the core count starts at two and has promotion intensity at most its current count (including the ideal hub's rate one); its second moment is bounded by that of the two-founder rate-one Yule process at `t+h`, independently of `k`. The argument uses the full core/hub state for the physical forecast and only the currently present planted branches for the ideal one. In the coupling the enlarged history contains neither system's future randomness.

The inert virgin limit in `L²` gives `e^(−μu) C_0(u)→V_0` in `L¹` and a uniform bound on its discounted expected count. Campbell integration shows that the expected total of all terms from arrivals after `R` is bounded by a constant times `e^(−μR)`, both in the normalized count at time `u` and in the proposed limit series. For finitely many arrivals before `R`, integrate the individual `L¹` convergence against their arrival times (also using the uniform bound). The normalized ideal count therefore tends to `Y_∞` in `L¹`, and the positive series is finite almost surely. For any fixed `t`, the ideal forecast as `h→∞` tends in `L¹` to `E[Y_∞ | ideal history at t]`; the latter tends to `Y_∞` in `L¹` as `t→∞`. Combining the estimates in the order **first `k→∞` at fixed `t,h`, then `h→∞`, then `t→∞`** gives `Y_k⇒Y_∞` (indeed convergence in the distance tested by all 1-Lipschitz functions). This proves the required interchange without a bound on the number of initial ordinary leaves.

## 4. Higher inverse-moment bound, including the size normalization

Choose `K₀` large enough that the large-core estimate of [left-tail §3](../left/README.md#3-stopped-downward-passages) gives a uniform negative moment of order `2/μ` from any seed with at least `K₀` cores. From any `S_k`, before hitting `K₀` cores, first wait if needed until size `K₀+3`; the total size rates are at least four. At and above that size while `2≤C<K₀`, the [promotion hazard](../left/README.md#4-star-departure-is-the-unique-rate-one-bottleneck-needed-here) is at least `1+δ` for `δ=min(1/6,1/K₀)`. Each interval to a new promotion has conditional exponential moment of order `1<r<1+δ` uniformly bounded, with at most `K₀−2` such intervals. The conditional bound follows by integrating the waiting-time survival estimate, valid along nonpromotion jumps as well. Multiplying the moment bounds at the successive stopping times and then using the strong Markov property at the large-core hit gives

```text
sup_(k≥4) E Y_k^(−r/μ) < ∞.
```

For `D_k=k^μ X_k`, independence of the free Yule clock gives `E Y_k^(−r/μ)=E[(W_k/k)^(−r)] E D_k^(−r/μ)` (Tonelli), and the first factor is finite and at least one since `k≥4>r` and `E W_k=k`. Thus the bound holds for `D_k` as well. It implies UI for the inverse `α` moments of both families, and convergence of those moments to `E Y_∞^(−α)` by weak convergence and truncation. In particular using the minimum of the two initial core rates at a four-vertex departure as if it were uniformly above one would be incorrect: the size gate handles that case.

## 5. Summation at the star scale

Write `h=x^α` and `F_k(u)=P(D_k≤u)`, `F_∞(u)=P(Y_∞≤u)`. The exact mixture is

```text
(1/h) P(X≤x) = Σ_(k≥4) [2/(h(k−2)(k−1))] F_k((kh)^μ).
```

Identify each summand with the integral on `[kh,(k+1)h)` of its value divided by `h`. For almost every `u>0`, `k=floor(u/h)→∞` and `F_k((kh)^μ)→F_∞(u^μ)`, by weak convergence at continuity points and monotonicity. The resulting step integrand is bounded by `const·u^(−2) min(1,u^r)` using the uniform inverse moment; here `k≥4` ensures `kh` and `u` are uniformly comparable. This is integrable at both ends because `r>1`. Dominated convergence and Tonelli give

```text
lim (1/h)P(X≤x) = 2∫₀∞ F_∞(u^μ) du/u² = 2 E Y_∞^(−α).
```

The same dominated convergence, with a bounded continuous test of `kh`, identifies the departure size under the small-tail event: conditional on `X≤x`, the variable `x^α K` converges in law to the density

```text
g(u) = F_∞(u^μ) / [u² E Y_∞^(−α)],   u>0.
```

This avoids inverting the physical Yule smoothing transform. For another description of the constant, with `L(q)=E exp(−q V)`, the Poisson exponential formula gives

```text
E exp(−q Y_∞) = L(q) exp(−α∫₀^q [1−L(v)] dv/v).
```

The integral is finite at zero by `E V<∞`. Integrating this Laplace transform against `q^(α−1)/Γ(α)` gives `E Y_∞^(−α)`; this is a characterization, not a numerical evaluation.

## Exact finite controls

```text
$ taskset -c 8 nice -n 10 timeout 8s python3 -B round2/open/limit-law/prefactor/check.py
PASS: exact departure and survival weights through size 256
PASS: 3159 exact hub degree/core-neighbour states
PASS: exact fixed-size mixture with a two-atom toy amplitude
```

The two-atom fixture exercises the departure-mixture normalization and the factor two. The convergence of the physical eigen-amplitudes and the inverse-moment bound are analytical.
