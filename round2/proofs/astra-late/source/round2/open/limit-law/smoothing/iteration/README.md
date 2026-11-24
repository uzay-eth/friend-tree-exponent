> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../../../README.md#claim-17) for exact scope and citations.

# Iterating physical cuts: stronger full-law integrability and continuous sublaws

For the edge-seeded amplitude `X`, the full density is **strong `L^{291/278}`**. More precisely the global all-set exponent is `13/(278+13μ)` and the local exponent on positive compact sets is `13/278`. This improves the [previous full-law bound](../../../../../../../../reviews/cold/astra/checks/limit-compact-sol.md). For a fixed deterministic time a physical six-side cut with an additional shape-only gate gives a **bounded continuous subdensity**, with explicit discarded-mass and inverse-moment bounds. These estimates do not establish local boundedness or continuity of the full density.

## 1. A sharper two-side kernel

For `P~Beta(n₁,n₂)`, put `P=sin²θ`, `0<θ<π/2`, `N=n₁+n₂`, where both `n_i≥2` are the *actual* physical component sizes. The density of `θ` is bounded by `C√N`, uniformly even for unbalanced sizes. Indeed it is

```
2 Γ(N)/(Γ(n₁)Γ(n₂)) · (sin²θ)^(n₁−1/2) (cos²θ)^(n₂−1/2).
```

Take its maximum at `sin²θ=(n₁−1/2)/(N−1)`. Set `u=n₁−1/2`, `v=n₂−1/2`; Stirling bounds give `Γ(u+1/2) ≍ u^u e^(−u)`, similarly for `v`, and `Γ(N) ≍ (u+v)^(u+v+1/2)e^(−u−v)`. All powers except `√(N−1)` cancel. Unlike a bound on the density of `P`, no minimum-size factor is lost.

Let `F(θ)=a sin^(2μ)θ+b cos^(2μ)θ`, with nonnegative `a,b` and `a+b>0`. The following curvature formula is valid on the open interval, including cases with a zero coefficient:

```
F'' = −4μ(1−μ) F + (2μ−1)(cot θ−tan θ) F'.
```

For any Borel set `A⊂[l,∞)`, trim angular distance `ε` at both endpoints (cost `2ε` in length). On the remaining region `F≥l`, wherever `|F'|≤h=c_μ l ε`, the last term has magnitude at most `2μ(1−μ)l` if `c_μ` is sufficiently small. Thus `F''≤−2μ(1−μ)l` there. The region `F≥l` is an interval since the map is unimodal (`a p^μ+b(1−p)^μ` is concave). Within its trimmed interval the derivative strip has at most one interval of positive length: at either level `±h` the derivative crosses only downward. Its length is at most `2h/[2μ(1−μ)l]`. Off that strip the at-most-two inverse branches cost `2|A|/h`. Choosing `ε=√(|A|/l)` proves an angular length bound `C_μ√(|A|/l)`; if this choice is too large for trimming, the whole interval satisfies the same bound by increasing `C_μ`. The argument uses only unimodality, the displayed curvature identity, and the area change on monotone branches. It neither assumes a coefficient lower bound nor drops angles near zero at zero cost.

Consequently, on the adapted no-cross history split, the exact [shape-weighted Dirichlet law](../README.md) gives, after averaging `E√|T_t|≤√3 e^(t/2)`,

```
P(X∈A) ≤ C e^(−ct) + C_μ e^(t/2) √(|A|/l),       A⊂[l,∞).
```

The shape weight is at most one. There is no `E[√N | no cross]` and no inverse moment of an individual coefficient in this step.

## 2. Actual retry cost and full-law powers

Use the [fast two-anchor race](../compact/README.md), restarting on observable failure. With `w=e^(−u/2)`, the defective tail of the first extinction of the two independent lower stocks is `q+q²`, `q=w/[2(2−w)]`. Its coefficient at `w^j` is `(j+1)/2^(j+2)` for `j≥1`. At `r=13/200` sum the first eight exponential integrals and bound the rest by placing their remaining coefficient mass at the ninth decay rate. Including the defective atom mass in the expectation gives

```
F_r ≤ 5715644171206736442597875 / 7046497926887658264039168,
ρ ≤ [1 + (1+2F_r)/(3(1−r))]/(2−r)
  = 11951004048194066431207745075 / 11951906448537490802583997242 < 1.
```

Multiplying these defective moments at actual failures bounds the completion time `S` of the first never-failed setup by `P(S>u)≤C e^(−ru)`. This is an existence witness, not a stopping time used for conditioning. At time `t`, set `α=100/126`, `c=13/252`, `K=ceil(exp(ct))`. The late setup, small-positive lower stock, and any future cross are bounded respectively by `C exp(−rαt)`, `C K exp(−(1−α)t/2)`, and `4/(K+1)`, all `O(exp(−ct))`. The finite-age bound is summed over the attempts reached and uses the unconditional lower-stock distribution at setup, not conditioning on which attempt is never-failed. The final selector sees only the physical state at the deterministic `t`.

Balancing yields the local exponent `c/(1+2c)=13/278`. With the [full-law left-tail bound](../../left/README.md) `P(X≤l)≤C l^(1/μ)`, use `l=δ^(μξ)`, `t=(ξ/c)log(1/δ)`, `δ=|A|`, to obtain `ξ=13/(278+13μ)`. Testing density superlevel sets gives weak `L^{(278+13μ)/(265+13μ)}` and strong `L^p` strictly below it. Since `μ<1`, this includes `p=291/278`. It gives no strong claim at the μ-dependent weak endpoint.

## 3. Multiple physical successes at deterministic times

The following explicit path experiment works for any fixed `m≥2` (we use `m=6`). For each of `m−1` stages of two units, begin at the current terminal core, which has at most one core neighbour. Require all earlier anchors to be silent throughout the stage. In the first unit silence the active anchor and require two rings from a designated neighbouring leaf. In the second unit require the active anchor's first ring to occur by its end and to choose a leaf; then silence it and the newly promoted anchor for the rest of the stage. Each old anchor retains a leaf, and the new one starts with a leaf. External growth cannot remove these leaves. If `d=2(m−1)`, the conditional probability of completing the script from any finite state with a core is at least

```
q_m = [(3/4)e^(−3)(1−2/e)(1−1/e)]^(m−1) e^(−(m−1)(m−2)).
```

The last factor pays two units of silence per strictly earlier anchor at each stage. All wait/choice conditions are conditional on the physical history then available. At script completion embed `m` independent disjoint lower leaf families, one at each anchor, as in the two-side proof. Conditional probability all survive is `2^(−m)`; conditional probability of at least `j` surviving is `2^(−m) Σ_{i=j}^m binom(m,i)`. All surviving is sufficient for the path cut; a collection with intervening failed anchors is not asserted to give that cut. A failed script restarts at its end; a failed batch restarts at its first observed extinction. In a failed batch `P(t<U<∞)≤(m/2)e^(−t/2)` by union over the finite-extinction tails. If `p_m=q_m 2^(−m)` and

```
r_m = min(1/4, p_m / [2 e^(d/4)(d+2m)]),
```

the uniform defective restart moment is at most `1−p_m/2`: use the failure probability at zero and the tail integral for `e^(r_m D)−1`, splitting at `d`. Thus the probability the successful setup is later than `αt` is at most `C_m exp(−r_m αt)`. The geometric finite-age bound for each positive stock gives an additional `C_m K exp(−(1−α)t/2)`. With `K=ceil(exp(γt))` and `γ=r_m/(1+4r_m)`, `α=1/(1+4r_m)`, both are `O_m(e^(−γt))`. This controls all `m` anchors simultaneously at the **same deterministic time**, before using any terminal masks. The conditional bound on `m`-fold success is a product only for the auxiliary stock drivers at preparation, not for success events of different retrospective cuts.

On selection remove the `m−1` path edges, retaining slots; each component has exactly one anchor, with all its private leaves still present. To bound the extra requirement of **no future arrow from an anchor into any core** (including across a cut), first discard histories with maximum core degree above `Bt`, for a sufficiently large fixed `B`. The exponential moment sum for core degrees is at most `C e^{e t}` by the [personal-clock argument](../../../../../../../astra/source/round2/fringe/hub-spectrum.md#3-counting-small-limiting-masses); choosing `B>e+1` makes this error `O(e^(−t))`.

Here is a conditional stock estimate, also giving negative moments of real masses. Starting at a stopping time from stock `L≥K`, let `ζ` be the first drop below a fixed integer `k₀>4+2q`, for any chosen `q>0`. Stock birth intensity is at least `L` and stock loss intensity is at most one. For `h_q(l)=Γ(l−2−q)/Γ(l−2)`,

```
h_q(l+1)/h_q(l) = 1−q/(l−2),
h_q(l−1)/h_q(l) = 1+q/(l−3−q).
```

Hence `A h_q≤−q h_q` when `l≥4+2q`, and killing at `ζ` only improves that inequality. Stopped Dynkin and `h_q(l)≍_q l^(−q)` give

```
E[e^(qu) h_q(L_u); u<ζ | history] ≤ h_q(L_0),
P(ζ<∞ | history) ≤ h_q(L_0)/h_q(k₀−1).
```

For the second line use the nonincreasing *stopped* test, not the killed test. On `ζ=∞`, Fatou and `D_u≥L_u` give `E[M_anchor^(−q); ζ=∞ | history]≤C_q L_0^(−q)` for `M_anchor=lim e^(−u)D_u`; zero mass on this event is excluded. The same bound applies to the continuous ordinary size mass of the actual component on the no-cross sublaw, since that side contains all the anchor's ordinary leaves. It holds for any prescribed finite order if the stock gate is chosen accordingly; it is **not** an unconditional inverse moment of a fixed hub mark, which can be zero.

Use `q=2` for the discarded-mass estimate. Conditional on history the anchor's core count is at most `Bt+R_u`, where `R_u` is its future rate-one ring count up to elapsed `u`. By Cauchy–Schwarz, the expected intensity of a core send before `ζ` is at most `C (Bt+u+√u)e^(−u)/K`. Integrating and adding the probability of `ζ<∞` bounds the chance of failure of this no-core-send requirement by `C_m (1+t)/K` in total. These bounds are on actual clocks without independent hub-future assumptions. Thus the fixed deterministic-time cut with this additional mask retains mass at least `1−C_m(1+t)e^(−γt)`.

## 4. Negative core-amplitude moments and continuous retained six-side laws

Let `Z_i,W_i` be the continuous core and ordinary masses of side `i`, in the time scale from the cut. Plain lower bounds on `W_i` are not lower bounds on `Z_i`. For the latter use the following additional physical argument on the no-anchor-to-core-send sublaw. In each independent inert-side construction give the anchor a rate-one clock including its inert slot choices, and use those slot choices as the reinstated outward crossing channels. Assign to every anchor ring an independent inert-virgin continuation in advance. When a ring promotes a leaf, use that continuation inside the new branch until the first incoming arrow from its parent. Outgoing child-to-parent arrows do not affect the branch internally. On the retained sublaw all anchor rings promote, and no such later parent arrow arrives. Counts are truly included in that cut side, not in the entire tree regardless of cut.

Fix `z_*>0` with `p_* = P(Z_J^0≥z_*)` (inert virgin) arbitrarily close to one, using its [positive limit](../../../forest-age/inert/README.md). The first preassigned mark meeting this threshold is at time `T_*~Exp(p_*)` in the **unrestricted auxiliary clock construction**. On the retained sublaw `Z_i≥z_* e^(−μ T_*)`; do not condition the auxiliary clock on the retained event. Thus for any `q μ<p_*`, its unnormalized inverse moment there is at most `z_*^(−q) p_* /(p_*−q μ)`, uniformly over the cut history. The product bound across sides is valid for their independent *unconditioned* inert drivers with suppressed crossing clocks, integrating a subprobability by domination. It is not an assertion that several success probabilities factor after selecting a future hub.

In the rate-`n` reference of the exact cut theorem, the extra ban on **internal** sends from the anchor to existing cores is an embedded-arrow condition (retain the initiator when sampling the successor). The stock-floor gate also depends only on the embedded path. Banning cut sends was already accounted for by the factorized likelihood `L`. Consequently these sublaws still have independent Gamma clocks given their embedded paths, with `W_i~Gamma(n_i,1)` and `Z_i=a_i W_i^μ`. For each side the preceding weighted `Z_i` moment therefore gives the genuine weighted coefficient bound

```
∫ a_i^(−q) dν_i ≤ C_q Γ(n_i)/Γ(n_i−μ q),       μq<min(p_*,n_i),
```

where `ν_i` is the shape measure including the terminal gate and that side's no-cross weight. This is an unnormalized moment, with no inverse of a success probability. Likewise the weighted moment of the ordinary size share `P_i` at order `q<n_i` is at most `Γ(N)Γ(n_i−q)/(Γ(N−q)Γ(n_i))` by the Dirichlet clock law and total weight at most one. These controls are for the specified true gate; they cannot simply be used for an arbitrary earlier-failure mask as an independence statement.

Pair the six sides arbitrarily. Apply the [three-pair convolution estimate](../compact/README.md#4-multi-side-test-and-the-remaining-continuity-input) conditional on paths, using `A_j^(−1/3)≤a_{i(j)}^(−1/3)` for one side in each pair and `M_j≤N_j` for its within-pair Beta maximum. Integrate the pair sums first; their factor is `Γ(N)/Γ(N−μ)∏ Γ(N_j−μ/3)/Γ(N_j)`. Then integrate the shape weights by the inverse coefficient bounds, whose product is at most `C∏ N_j^(μ/3)` and cancels the inverse pair-size powers up to constants. The result is bounded by `C N^μ ∏ N_j^(2/3)≤C N^(2+μ)`. Each conditional convolution is continuous, and the integrated bound makes dominated convergence legitimate. Averaging over the ordinary physical cut-time history gives a deterministic-time continuous subdensity with sup norm `≤ C e^((2+μ)t)`, by Yule moments. Its missing mass is at most `C(1+t)e^(−γt)`. This proves a quantitative continuous exhaustion as an approximation in total variation; the sublaws need not be increasing.

## 5. What iteration does and does not prove

On a deterministic time mesh the separate failure bounds are summable, so almost every history eventually satisfies all the deterministic-cut success events along that mesh. No multiplication of final success probabilities is justified: setup failures are observed at stopping times, but the terminal masks for different deterministic cuts interrogate overlapping futures. In particular a first-success decomposition, subtracting earlier cut events, need not retain the Gamma clock law for each term. Domination still holds. It costs at most the entire current kernel bound, growing like `e^((2+μ)t)`, and mere control of the earlier failures in probability does not multiply that supremum by their probability. Interpolating a decomposed piece's small mass with this large supremum can give finite `L^p` close to one, not a summable supremum. The stronger global result here comes from the unbalanced two-side kernel, not from asserting conditional independence of terminal masks.

There is a concrete general obstruction to promoting this approximation control alone to local boundedness. For any `d,c>0`, take an integrable density spike `f(x)∝χ(x)|x−1|^{−a}` with `0<a<c/(c+d)` and a continuous compact cutoff `χ` equal to one near 1. The continuous subdensities capped at `e^(dt)` have missing mass `O(e^{−d(1−a)t/a})=O(e^{−ct})`, yet the full law is locally unbounded at 1. This is an obstruction to the *approximation argument*, not a claim that the physical law has a spike. And the balanced-Dirichlet concentration example in the parent note still precludes a uniform conditional bound for arbitrary sizes. What is missing for full continuity is a stronger estimate on the physical cut mixture or a summable pointwise bound for the unresolved sublaws (including their true dependence on prior masks).

Exact finite arithmetic controls: `nice -n 10 timeout 12s python3 -IB round2/open/limit-law/smoothing/iteration/check.py`. Passed 64 defective-tail coefficient checks, the exact retry and exponent checks, 50 angular curvature checks, 279 integer Beta peak checks and 456 inverse-stock generator checks. The probabilistic arguments and density estimates above are analytic.
