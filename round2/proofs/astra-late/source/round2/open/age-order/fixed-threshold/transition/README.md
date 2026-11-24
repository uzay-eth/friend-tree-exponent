> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports. See the [release index](../../../../../../README.md#claim-24) for exact scope and citations.

# Adjacent tails between birth scale and fixed thresholds

For the original edge seed write `p_v(x)=Pr(Z_v>x)`. Along **every** sequence `x=x_v→0` with `vx→∞`,

```math
 \frac{p_v(x)-p_{v+1}(x)}{x p_v(x)}\longrightarrow 1.
 \qquad\text{(1)}
```

This includes `x=(log v)^a/v` for every `a>0`, `x=v^(-α)` for `0<α<1`, and arbitrarily slowly divergent `vx`. In uniform-window language, (1) is uniform on `l_v/v≤x≤u_v` whenever `l_v→∞` and `u_v→0`. The proof uses the exact signed two-birth identity rather than differencing two one-label equivalents. The leading term is the middle endpoint of a two-birth chain grown from an old leaf. Its locked reference gets one extra red unit, and the factor this buys is exactly `2vθ` at reference endpoint `θ`. The chain rate supplies the compensating factor `1/2`.

## 1. Two seeds, one tilted likelihood limit

Use a post-insertion size `m`, and mark either (i) a leaf freshly attached to an old leaf `u` (now of degree two), or (ii) the middle of a chain formed by attaching to `u` as in (i) and then to the mark. The mark in (ii) has the same old neighbour `u` and one additional fresh leaf. Put `a=1/2` in (i) and `a=3/2` in (ii), and use

```
R_i=a+S_i,    B_i=m-a+F_i,    Δ_i=R_i-r_i,
```

where `S_i,F_i` count tag and outside choices and `r_i` is the physical rate of gain at the mark. Then `Δ_0=0`, and `Δ_i` is nonnegative and nondecreasing: existing neighbours' reciprocal contributions can only decrease, and a newly added neighbour is a leaf. The [exact likelihood and two sandwiches](../uniform/README.md#1-exact-reference-and-physical-likelihood) apply for both seeds: given `Θ=θ`, the indicators are iid Bernoulli(`θ`), and the physical conditional outside arrow law has denominator `B_i+Δ_i`. The likelihood factors are `1−Δ_i/R_i` on a tag and `1+Δ_i/B_i` outside. The finite-history likelihood `L_i` gives an upper bound on the physical tail by integrating `E_θ L_i` against the beta law restricted to `Θ>x`; the integral of `E_θ L_∞` gives the lower bound. The upper bound uses the future locked Pólya domination, since `Δ_i≥0`; the lower bound uses Fatou on cylinders and the common path-functional endpoint. Neither sandwich requires infinite-horizon uniform integrability.

For either seed, at locally convergent old-side inputs, the claim is

```math
 \Pr_T(Z>x) / Q_{m,a}(x) \longrightarrow g,\qquad
 Q_{m,a}(x)=\Pr(\mathrm{Beta}(a,m-a)>x),\quad x\to0,\ mx\to\infty,
 \qquad\text{(2)}
```

where `g` is the **same** old-side and virgin-immigration functional as in [the fixed-interior proof, §3](../uniform/README.md). In detail, in the inert exterior clock let `P` be a rate-one Poisson process, independent of the inert old side with its slot at `u`; each arrival starts an independent inert virgin branch. Define

```math
 \delta(t)=\tfrac12-D_u(t)^{-1}+\sum_{s\in P,\ s\le t}(1-D_s(t)^{-1}),\qquad
 g=\mathbb E\exp\{-\int_0^\infty \delta(t)/(e^t-1)\,dt\}.
```

Here `D_s(t)` is the boundary degree of the branch born at time `s`. The functional is continuous in the old-side input by finite-time local construction and the integrable late bound `Eδ(t)≤1+t`. The initial bound `Eδ(t)≤2t` follows from the rate at most two of the first gain of `u` and rate one of immigration. In particular `exp(−π_math²/3)≤g≤1` by Jensen. Neither seed's initially present leaves are exterior branches at time zero.

### Finite windows with a vanishing Bernoulli parameter

Fix `θ=θ_m→0` with `mθ→∞` and use jump-index time `s=i/m`. On a bounded interval `0≤s≤M`, the maximal variance inequality gives `S_i/(mθ)→s` uniformly in probability, and `F_i/m→s`. The number of tag-origin outside arrows is bounded by an independent Bernoulli envelope `J_f` indexed by outside choices, with parameters `1/(m−a+f)`: use a fresh uniform at each outside choice for the tag-origin decision, whose rate weight is at most one. Thus `Δ_i≤1+J_{F_i}` and `E[J_f^k]≤C_k(1+log(1+f/m))^k` for fixed `k`. Before `u` first gains or the first tag-origin arrow, `Δ_i=0`. The conditional probability of either event by `εm` is `O(ε)`, since their combined outside weight before the first event is at most three and the outside denominator is at least `m−a`. This holds also for the chain middle.

Here are the local clock details. Conditional on the Bernoulli schedule, at each outside slot the old-edge and internal channels retain their physical weights; on a bounded local population and envelope stop, replacing their denominator by `m−a+F_i` costs `O_M(m^-1)` on any fixed ball. The outside counting measure divided by `m` tends to `ds`; sums of a fixed channel hazard converge to its integrated hazard in time `t=log(1+s)`. The maximum proposal probability tends to zero, so the no-jump products converge to the corresponding exponential survivals; iterating the next-jump law gives the joint stopped evolution. Cutoffs can be removed as follows. For a fixed ball the expected tracked population on `i≤Mm` is bounded by the discrete rate-one population domination with bounded boundary input, using per-slot denominator at least `m−a`. The tag cluster uses the same bound with rate-one immigration. For the exterior cutoff, along any ordered old path of length `h` the conditional proposal probability is bounded by `(Mm/(m−a))^h/h!` times the product of inverse old core degrees of the senders; the [reversed-path sum](../../../../../../../astra/source/round2/operator.md#1-cells-and-evaluation-of-the-operator) bounds the remaining sum by the old core degree at each readout cell. This also bounds the effect of an inbound boundary treatment. On `ε≤s≤M`, the tag has degree tending to infinity in probability while its nonleaf neighbours are bounded in probability. Its rate-one initiation therefore promotes a leaf with probability tending to one; the cumulative weight of tag sends to existing core neighbours on that interval tends to zero. On `s≤ε` the probability of any tag send is `O(ε)`. Promotions start independent virgins in the limiting clock, with the mark's outgoing slots suppressed thereafter except for further immigration. Thus the old-side and branch process producing `Δ_i` tends to the process for `δ(t)` above. The extra initial leaf in (ii) disappears from these clock and deficit calculations; it only changes the initial red count.

For clarity, on `ε≤s≤M`, stop at a fixed bound for `Δ_i` and where `R_i` is at least `ε mθ/2`. The sum of conditional squares of the likelihood increments is `O_{ε,M}(1/(mθ)+1/m)` there; this bounds both the log expansion error in probability and the martingale part of the linearized sum. Its drift is exactly

```math
 \Delta_i[(1-\theta)/B_i-\theta/R_i]
   =\Delta_i\frac{a+S_i-\theta(m+i)}{R_i B_i} .
 \qquad\text{(3)}
```

Consequently the drift sum tends to `∫ δ(log(1+s))[(1+s)^−1−s^−1]ds` over the same interval. The substitution `t=log(1+s)` yields the kernel `−1/(e^t−1)`. Remove the initial `ε` using the first-defect bound and `Eδ(t)≤2t`. This proves convergence of `L_floor(Mm)` in law to the truncated exponential. It also converges in expectation: by dropping the tag factors, `L_floor(Mm)≤exp[C_M(1+J_floor(Mm))]`, and the envelope has all fixed exponential moments uniformly in `m`. The latter bound holds for **all** `θ`.

### Late window: the controlling parameter is `mθ`, not `mθ²`

For `w=1+i/m` and `i≥Mm` (`M≥1`), stop at the first failure of

```
|S_i−θi| ≤ H sqrt(mθ) w^(3/4),       J_(F_i) ≤ H w^(1/8).
```

The union probability of failure is `≤C/H²`: on dyadic blocks apply the maximal variance inequality to the centred Bernoulli sums (variance at most `θi`), and the second-moment bound to the increasing `J_f` at the block's right endpoint. On the stopped window, for `θ≤1/2` and `mθ` sufficiently large at fixed `H`, one has `R_i≥c mθ w`, `B_i≥c m w`. In (3) bound the numerator before the factor `Δ_i` by `a+mθ+H sqrt(mθ)w^(3/4)`. Writing `η_i=ℓ_i−1` for the one-step likelihood increment, the resulting sums of predictable bounds are

```math
\begin{aligned}
 \sum_{i\ge Mm}|\mathbb E_\theta[\eta_i\mid\mathcal F_i]|
 &\le C\{H M^{-7/8}+H^2(m\theta)^{-1/2}M^{-1/8}\},\\ \sum_{i\ge Mm}\mathbb E_\theta[\eta_i^2\mid\mathcal F_i]
 &\le C H^2 (m\theta)^{-1} M^{-3/4}.
\end{aligned}
\qquad\text{(4)}
```

In the second line the one-step bound is `Δ_i²(θ/R_i²+(1−θ)/B_i²)`; in particular one must keep the Bernoulli factor `θ` on the red variance. On the stopped window `sup|η_i|<1/2` for large `mθ`; `|log(1+η_i)−η_i|≤2η_i²`. Summation, Markov's inequality and the maximal martingale variance bound control the log increment uniformly over all finite end times on the event of no stop. Take the large-`m` limit, then `M→∞`, then `H→∞`. This also controls the increment up to `L_∞` for beta-almost every `θ`, by the unconditioned reference's martingale convergence. The drift estimate keeps the two signs in (3) together; bounding the separate red and outside series would not give (4).

To integrate, under either beta law conditioned on `Θ>x`, the overshoot is `O_P(1/m)` and `Θ∈[x,2x]` except with probability tending to zero (see the endpoint estimates below). All the preceding convergence statements are sequential for `θ→0`, `mθ→∞`, and therefore uniform on `[x,2x]` along the input sequence. For the upper sandwich use the finite-window expectation and its bound on the discarded parameter portion, then let `M→∞`. For the lower sandwich restrict to `[x,2x]`, and use bounded continuous truncations: on the good event the limit likelihood is at least the finite-window value times `e^{-ε}` for arbitrary log tolerance `ε>0`; the discarded good-event complement costs at most the truncation level times its probability. Use (4) and send `M→∞`, `H→∞`, `ε↓0`. The limiting truncated exponential is at most one, so these lower and upper bounds coincide. This proves (2), including along locally converging sequences of old-side inputs. There is no required lower bound on `mθ` beyond divergence.

## 2. Beta bounds and the signed cancellation

In the range `x→0`, `y=mx→∞`, endpoint integration gives, for fixed `a>0`,

```math
 Q_{m,a}(x)\sim \frac{m^{a-1} x^{a-1}(1-x)^{m-a}}{\Gamma(a)} .
```

For the bounds below one can substitute `t=x+(1-x)z/m` in the density integral. The remaining factor is `(1+(1-x)z/(mx))^(a-1)(1-z/m)^(m-a-1)` on `0≤z≤m`. It is dominated by `(1+z)^k e^{-z/2}` for bounded `a`, large `m` and `mx≥1`, with a fixed `k` for the range in use. A matching lower bound on `0≤z≤1` is uniform in that range. This proves the equivalent, the overshoot bound `E[Θ−x|Θ>x]=O(1/m)`, and the uniform estimates

```math
 Q_{m,a}(x)\le C a\,y^{a-1}(1-x)^{m-a}\quad(0<a\le3/2),
```

using `1/Γ(a)≤C a` and the uniform gamma-ratio bound `Γ(m)/Γ(m-a)≤C m^a`. In particular `Q_{m,a}≤C a Q_{m,1/2}` for `a≤1/2`, with an extra factor `y^-1/6` when `a≤1/3`. Likewise `Q_{m,a}≤C y^-1/6 Q_{m,3/2}` when `a≤4/3`. At the two different sizes needed for cancellation the **densities** satisfy

```math
 f_{v+1,3/2}(t)=2vt\, f_{v,1/2}(t),\qquad
 Q_{v+1,3/2}(x)=2vx\,Q_{v,1/2}(x)(1+O(1/(vx))). \qquad\text{(5)}
```

Now start at `T_n`, `n=v−1`, and denote by `S_v` the expectation of the sum over old leaves `u` of `r_u g(T_n,u)`, where the old-side input is formed by promoting `u` and leaving a slot for its marked neighbour. Equivalently this is the sum over old cores `o` with weights `L_o/D_o` times the corresponding `g`. The existing deterministic-size local law and core-count uniform integrability give `S_v/v^μ→EX·π((L/D)g)>0`. Use (2) and the same empirical convergence with bounded per-core scores to replace `g` by the corresponding normalized physical probability; both sums equal `S_v+o(v^μ)`. In detail, the sequential convergence in (2) for every converging old-side input gives convergence of the bounded varying scores under weak convergence of the rooted empirical laws. Multiply by `C_n/v^μ` using uniform integrability. No speed of empirical convergence is used.

For the single birth, core targets have `a=1/(D_u+1)≤1/3` and `r_u a≤1`; their contribution before division by `n` is `O(E C_n·Q_{v,1/2}(x)(vx)^-1/6)`. Therefore

```math
 p_v(x)=\frac{Q_{v,1/2}(x)}{n}\{S_v+o(v^\mu)\}.
\qquad\text{(6)}
```

For the signed difference use the [exact two-birth identity](../../adjacent/README.md#exact-cancellation-in-two-births) with the common final-state endpoint test `Z>x`. In every old-edge term both endpoints are leaves with locked `a≤1/2`. The same is true for every chain tip. A chain middle with first attachment to an old **core** has locked `a=1+1/(D_u+1)≤4/3`, since it has one fresh leaf and one old neighbour. All these bounds are physical upper comparisons by monotone deficit. Since the absolute two-birth coefficients charge at most three per old core (count each endpoint separately for a bound of six), the omitted terms in the signed numerator are at most `O(E C_n· Q_{v+1,3/2}(x)(vx)^-1/6)` in absolute value. This bound allows either sign. The remaining chain has first attachment to an old leaf; its middle is exactly seed (ii), and its coefficient is `r_u/(D_u+1)=r_u/2`. Thus

```math
 p_v(x)-p_{v+1}(x)
 =\frac{Q_{v+1,3/2}(x)}{2n(n+1)}\{S_v+o(v^\mu)\} .
\qquad\text{(7)}
```

Dividing (7) by (6) gives (1) using (5) and `n+1=v`. In particular, the same `g` matters at **relative error `o(1)` on the signed numerator itself**. No estimate comparing adjacent one-label remainders to `x` has been used. Equation (6) also extends the fixed-interior one-label prefactor to this transition range with the same constant `κ`.

## Scope and verification

The signed ratio error consists of the two normalized tilted-score errors (vanishing, without a rate here), `O((vx)^-1/6)` from the leftover sectors, and `O((vx)^-1)` from the beta overshoot in (5). Equation (4) gives explicit late-window control but does not furnish a numerical starting label. For `vx` bounded use [Claim 18's compact-scale comparison](../../scaled/README.md); for fixed interior thresholds use the [fixed-interior result](../uniform/README.md). The three statements imply eventual strict adjacent order uniformly on `0≤x≤q` for every fixed `q<1`: a countersequence has a subsequence with either a positive interior limit, a zero limit with bounded `vx`, or a zero limit with `vx→∞`. This leaves thresholds tending to one outside this argument.

`python3 -IB round2/open/age-order/fixed-threshold/transition/check.py` passed 466 likelihood prefix checks (all depth-three words for each seed), 36 beta tilt moment identities, and 32 signed two-birth rows. The checks are algebraic controls; the moving-window assertion uses the analytic estimates above.
