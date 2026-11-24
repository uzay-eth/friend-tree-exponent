> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# Exact near-one constant from the first outside attachment

Fix a finite physical friend tree `T` of size `n₀≥2` and an existing tagged vertex `v`. This includes conditioning on any reachable birth tree of a fixed label. In this note `Z` is the limit of the tag's degree divided by the **whole physical size**. Write

```text
R = Σ_(y~v) 1/D_y,      α = R,      β = n₀−R,
N = n−D_v,             A = Σ_(y~v, D_y≥2) (1−1/D_y).
```

In particular `α,β>0`. If `K` is the size **before** the first attachment outside the tag, then

```math
\mathbb P_T(Z>1-\varepsilon)\sim
\frac{\Gamma(n_0)}{\Gamma(\alpha)}\,\mathbb E[Y^{-\beta}]\,\varepsilon^\beta .
```

The positive finite amplitude `Y` is specified below from inert old branches and independent Poisson immigration, started from the **physical first-defect mixture**. Conditional on `Z>1−ε`, `εK` has limiting density `β u^(−β−1) P(Y<u) / E[Y^(−β)]`, `u>0`. Recording the first size *after* departure instead gives the same limit. We also identify the joint departure type and continuation amplitude.

**Credit.** Fable Claim 5(e), `88e13f36b39162f6ee60561a2b465a982b34af02` (the pinned local `ft2-fable` file `round2/newground/stochastic-domination/COUPLING.md`, §5), proved the Pólya upper comparison and all-attachments matching order for fixed labels. The [independent endpoint audit](../../../attacks/fable5/endpoint/README.md) obtained the fixed-restart exponent `n₀−R` by the monotone reciprocal deficit. The argument here obtains an asymptotic equivalent, not merely an order. Its architecture follows the [physical `X` left-prefactor proof](../../limit-law/prefactor/README.md); the present first-defect kernel, all-time deficit martingale and inverse-moment estimate are different and do not rely on that proof's core-growth martingale estimates.

## 1. Exact departure weights and the marked restart

In physical time each vertex initiates at rate one and picks a uniform neighbour; its choice adds a leaf at the recipient. The gain rate at the tag is `D_v−A`, and the outside rate is `N+A`. Before departure all new vertices attach to the tag. Therefore `N=n₀−D_v(T)` and `A=D_v(T)−R` are constant on this prefix, with outside rate exactly `β`. The total rate at size `j` is `j`. Thus for integers `m≥n₀`,

```math
S_m:=\mathbb P_T(K\ge m)=\prod_{j=n_0}^{m-1}\frac{j-\beta}{j}
=\frac{\Gamma(n_0)\Gamma(m-\beta)}{\Gamma(\alpha)\Gamma(m)},\qquad
w_m:=\mathbb P_T(K=m)=S_m\frac{\beta}{m}.
```

In particular departure is almost sure and `w_m ~ β C m^(−β−1)` with `C=Γ(n₀)/Γ(α)`.

Remove the tag and ignore its leaf neighbours. Each other neighbour roots an old branch, with an exterior slot retained in its degree. Before departure these old branches have their original fixed shape. For the *physical* outside arrow conditional on departure at pre-size `m`:

| Arrow type | Conditional weight | Post-defect branches |
|---|---|---|
| Directed genuine arrow `x→z` within an old branch | `1/(β D_x)` | Add a leaf at `z` in that branch. |
| Tag chooses any of its leaf neighbours | `(D_v−c)/(β D_v)` | Keep the old branches and add one virgin (promoted leaf with its new leaf). |
| Tag chooses an existing core neighbour `y` | `1/(β D_v)` for each `y` | Add a leaf at that old branch's root. |

Here `D_v=D_v(T)+m−n₀` and `c` is the fixed number of old core neighbours. No other arrows can produce a first outside attachment. In particular, old branch roots can send into the tag but those are *prefix* arrivals, not defect arrivals. Leaf arrows into the tag are also prefix arrivals. Summing genuine internal arrows gives `Σ_B (|B|−1/D_root)=N+A−1=β−1`. Consequently the limiting type `J` has weight `1/β` for a new virgin and the listed `1/(β D_x)` for each internal old arrow; tag-to-core types vanish in this limit. If there are no old branches (`β=1`), the defect is always a virgin. These probabilities follow from the directed-arrow winner law of each jump; they do not use an unweighted random choice among exterior recipients. For the shape state after departure, leaf identity among the tag stock is immaterial.

## 2. The limiting amplitude

An inert planted branch preserves one root slot: choosing the slot causes no internal change and it emits no input into the branch. Let `W(B)=lim e^(−t)|B_t|` from a finite planted seed `B` and let `W₀` denote its virgin law. Use the random initial branch list `ℬ_J` from the limiting defect type. Run all its branches independently given the type, independently of an additional rate-one Poisson process with epochs `s_i`, and at each epoch launch another independent virgin. Then

```math
Y=\sum_{B\in\mathcal B_J} W(B)+\sum_{i\ge1}e^{-s_i}W_i,
```

with `W_i` i.i.d. copies of `W₀`. Each inert size `k` has genuine rate `k−1/D_root∈[k−1,k]`, so its amplitude exists, is positive and has mean at most its starting size. For a direct proof, wait for each successive genuine size change with an independent unit exponential divided by its rate, using fresh independent randomness for the arrow winner given the state. The corrections to the Yule holding reciprocals `1/k` are nonnegative and at most `1/[k(k−1)]`; their total is finite a.s. The centred Yule sum has summable variances, and the harmonic sum minus log size converges. Hence the time to size `k` minus `log k` converges to a finite real value; the successive gaps vanish almost surely. The mean bound follows by Fatou and Yule domination. Equivalently use the [inert size limit](../../forest-age/inert/README.md#3-ordinary-size-and-the-high-cutoff-joint-law). Campbell's formula now gives a finite mean for the immigration sum. There is always at least one initial branch of size at least two, so `Y>0` a.s. Its law has no positive atoms: for such a branch the first genuine wait is a continuous exponential independent of the winning arrow and subsequent elapsed-time evolution. Its amplitude is a strictly positive independent continuation amplitude times the exponential of minus that wait. Given the initial type this branch is independent of the other branches and the immigration. A finite mixture retains atomlessness. The summation proof only needs continuity almost everywhere, however.

Let `T_m^+` be the actual post-defect seed conditional on `K=m`, of total size `m+1`. In the physical continuation reset time to zero, and let `D_v(t)` be the tag's degree and `N_m(t)=M_m(t)−D_v(t)` where `M_m(t)` is total size. The exact decomposition is

```text
N_m(t)=1+Σ_(branches behind current core neighbours) (size−1).
```

It is **not** the sum of branch sizes plus one. Use the same tag clock in the ideal and physical trees. A ring that chooses a leaf starts an inert virgin in the ideal tree, and its actual new branch can be coupled with the same internal drivers. At such a birth the chosen former leaf has no earlier core branch to preserve. Share genuine internal arrows and root-to-slot choices. A physical root-to-tag send adds a tag leaf and is null for the branch; it does not stop the comparison. Tag leaves never affect the rest except through later tag promotions, which get fresh drivers. The first tag ring into a core neighbour can stop agreement. Up to that stop the number of eligible core neighbours is at most `c+1+P(t−)` for the tag's rate-one clock `P`, and tag degree is at least `D_v(T)+m−n₀`. Over `[0,t]` the failure probability is bounded by `[(c+1)t+t²/2]/[D_v(T)+m−n₀]`. Also couple the initial list by matching the finite type probabilities, with mismatch `O(1/m)`. These bounds concern unconditioned comparison marginals: continue the ideal clocks even on failure. They give convergence of the fixed-time deficits (jointly with type) to the ideal deficit defined by the same `size−1` formula. This does not claim that the physical tag leaf stock has a finite limit.

## 3. Uniform calendar-time tails and inverse moments

For each seed in the departure mixture, at every later state `N_m` increases by one at rate `N_m+A_m` and never decreases. The count of core neighbours, hence `A_m`, is bounded by `c+1+P(t)`, and `N_m(0)=n₀−D_v(T)+1` is fixed. In particular `E N_m(t)≤C_T e^t` by the drift equation (or its stopped form and monotone bounds), uniformly in `m`. The integrating-factor decomposition is

```math
e^{-t}N_m(t)=N_m(0)+\int_0^t e^{-u}A_m(u)\,du+L_m(t),
\qquad
\langle L_m\rangle_t=\int_0^t e^{-2u}(N_m(u)+A_m(u))\,du.
```

The expected tail of the positive drift is at most `C_T(1+t)e^(−t)` and that of the quadratic variation is at most `C_T e^(−t)`. Thus the amplitude `Y_m=lim e^(−t)N_m(t)` exists almost surely and its `L¹` distance to the time-`t` value is at most `C_T[(1+t)e^(−t)+e^(−t/2)]` uniformly over the mixture. The same estimate applies in the ideal system: intrinsic growth in a size-`k` branch adds one to the deficit at rate `k−1/D_root`, and the immigration clock adds one at rate one. Their total is exactly the ideal `N+A`, with `A` summing `1−1/D_root` over branches. Its limit is the stated `Y`: use `L¹` convergence for each inert branch (Yule domination gives a uniform second moment of its rescaled size), Campbell for the tail of immigration after a fixed time, and the correction `e^(−t)(1−#branches)→0` in `L¹`. Let `m→∞` first on fixed windows and then let `t→∞` using the preceding uniform tails. This proves `(J_m,Y_m)⇒(J,Y)`.

The random total size in the physical continuation is Yule from `m+1`, with amplitude `G_m∼Gamma(m+1,1)`. In the same realization `Y_m=G_m(1−Z)`, so `V_m=m(1−Z)` conditional on departure at `m` converges jointly with `J_m` to `Y`. Only `G_m/m→1` in probability is needed; in particular no independence of deficit and size amplitudes is assumed.

To control rare small amplitudes uniformly, lock the reciprocal deficit at its **pre**-departure value `A₀`. It is nondecreasing at all jumps, and at departure `N` increases by one. Thus at post-size `m+1` the shifted deficit `N+A₀=β+1`. From then on its increase probability is at least `(N+A₀)/n`. Couple it above a fractional Pólya count with initial red `β+1` and white `m−β>0` using common uniforms. Therefore `1−Z` dominates `Beta(β+1,m−β)` at every possible post-defect seed, and for any fixed `β<p<β+1`,

```math
\mathbb E[V_m^{-p}]\le
\frac{\Gamma(\beta+1-p)}{\Gamma(\beta+1)}
\frac{\Gamma(m+1)}{m^p\Gamma(m+1-p)}.
```

The right side is uniformly bounded for `m≥n₀` by the gamma-ratio asymptotic and the positive finite initial values. This conditional bound includes the vanishing-mass tag-to-core types; it does not need a separate estimate for them in the tail asymptotic. Portmanteau gives `E Y^(−p)<∞`. In particular the constant is positive and finite and inverse moments of order `β` are uniformly integrable. The strict extra unit of outside weight comes from the first defect, not from a presumed positive jump of `A` (which can stay unchanged).

## 4. Summation and conditional departure law

The exact mixture for the strict tail is `Σ_m w_m P(V_m<εm)`. For a bounded continuous test `f(u,j,y)` of scaled pre-size, type and scaled continuation, the corresponding summands are `w_m E[f(εm,J_m,V_m); V_m<εm]`. After division by `ε^β` identify them with step integrals over `[εm,ε(m+1))`. The integrand at almost every `u>0` tends to

```math
\beta C u^{-\beta-1}\,\mathbb E[f(u,J,Y); Y<u].
```

Indeed the moving indicator can be squeezed at continuity points under joint weak convergence. The inverse moment bound dominates the integrand by `C_f u^(−β−1) min(1,u^p)`: on occupied steps `u` and `εm` are uniformly comparable. This is integrable. With `f=1`, Tonelli integrates `β u^(−β−1)` above `Y` to `Y^(−β)`. This cancels the **β** in the local departure rate and gives the claimed constant. For general `f`, divide the resulting integral by `C E[Y^(−β)]`. Equivalently, given the near-one event the limiting `(J,V_m)` has the law of `(J,Y)` tilted by `Y^(−β)/E[Y^(−β)]`, and the limiting `εK` is that tilted amplitude times an **independent** standard Pareto variable of index `β` (survival `r^(−β)` for `r≥1`). In particular `(1−Z)/ε=V_K/(εK)` has limiting CDF `u^β` on `[0,1]`, independent of the tilted amplitude and type. The type is therefore not generally distributed by the original first-defect kernel under the near-one conditioning.

## 5. Fixed labels from the initial edge

At birth of `v≥3`, `n₀=v` and the sole neighbour's post-birth degree is `d≥2`. Thus `α=1/d`, `β=v−1/d`. The formula above is conditional on the entire reachable birth tree and uses its actual old branch, not a branch chosen just from the degree `d`. There are finitely many birth histories at fixed `v`. The smallest exponent `β*=v−1/2` occurs with positive probability (attachment to a leaf); all `d>2` terms are negligible at that scale. Precisely,

```math
\lim_{\varepsilon\downarrow0}\varepsilon^{-\beta_*}\mathbb P(Z_v>1-\varepsilon)
=\frac{\Gamma(v)}{\Gamma(1/2)}\,
\mathbb E\!\big[\mathbf1_{\{d=2\}}\,\mathbb E[Y_{T_v}^{-\beta_*}\mid T_v]\big].
```

Here the integrand is taken to be zero on `d>2`. The unconditional departure law is the corresponding mixture of the conditional joint limits, weighting each `d=2` tree by its prior probability times its prefactor. For either initial label start from the size-two edge; `α=β=1`, there are no old core branches, and the initial defect is a virgin with probability one. The prefactor is then `E[Y^(−1)]`, with the same law for the two labels by symmetry. These are fixed-label, not uniform-in-label asymptotics.

For a numerical characterization only, let `ℒ_B(q)=E exp(−q W(B))`. The Laplace transform of `Y` is

```math
\mathbb E\prod_{B\in\mathcal B_J}\mathcal L_B(q)\cdot
\exp\!\left[-\int_0^q(1-\mathcal L_{\rm virgin}(z))\,\frac{dz}{z}\right].
```

The integral is finite at zero from the first-moment bound; integrating the transform against `q^(β−1)/Γ(β)` gives the inverse moment. No numerical constant is required for the proof.

## Exact replay

```text
taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/open/age-order/near-one-constant/check.py
PASS 2616 perfect prefixes and 25042 directed outside arrows: deficit, Palm weights and postdefect shift
PASS 480 exact hazard/telescoping steps and inverse-tilted Pareto normalization
```

The checker enumerates ordinary labelled trees and computes the rates of the genuine outside arrows, including returned-arrow exclusion, first-virgin creation, tag-to-core transitions and the unit shift. Calendar-time convergence and the prefactor are proved above, not inferred from the finite controls.
