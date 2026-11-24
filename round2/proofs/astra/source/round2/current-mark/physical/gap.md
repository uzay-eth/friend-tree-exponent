# Strict gap for the current-mark residual

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 13: Claude for transfer bookkeeping; GPT for the uncapped bound and root transfer; inputs include physical jump geometry and the Round 1 exponent and prior. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


For cap 20, §§4–6 prove a strict gap for the newborn residual, and an exponential bound on expected active count from arbitrary states. The proof is qualitative: the exponent below `μ` has not been assigned a numerical value. The binary extension is stated only where the needed drift and contraction checks are defined.

## 1. Dormant flux in the path weight

Choose `d_*>0`, put `z(D)=q(D)/(q(D)+d_*)`, and on an unmarked shape let

```
Z(S) = h(D_o,c_o) + z(D_o) Σ_Y Z(Y),
W_a(S) = Z(S)/(a+d_*),
W_[l,u](S) = [1/(l+d_*)−1/(u+d_*)] Z(S).
```

Here `c` includes the parent slot; the dormant weight is for raw emissions. Write `δ(D)=z(D)−z(D+1)` and `k(D)=δ(D)/d_*`. An intrinsic root gain decreases the `Z` weight of all old child subtrees by `δ(D)ΣZ(Y)`. For active root `a`, the corresponding emitted dormant weight is exactly `k(D)ΣZ(Y)`. Combining the two, the effective coefficient is

```
r δ(D) a/[d_*(a+d_*)] Σ Z(Y).
```

At an immediate child `Y`, its external gain `p` likewise changes the path factor out of it, and a grandchild dormant is emitted. One may bound the intrinsic path decrease everywhere by zero for an upper bound, but the remaining **positive** root and grandchild dormant terms must still be charged. An actual certificate needs to control these two terms, not simply reuse the permanent `W` inequality. On dormant states, the activation minus extra killing cancels exactly: `R_c W_[l,u]=L⁰ W_[l,u]+d_* W_[l,u]` for the displayed weight. In fact `(u−l)/(l+d_*)−u[1/(l+d_*)−1/(u+d_*)]=d_*[1/(l+d_*)−1/(u+d_*)]`.

A pointwise positive weight `W` satisfying `R_c W≤ρW`, growing at most polynomially in `n`, dominating `1_active` and with bounded expected centre immigration `I_t W=O(e^(βt))` for `β<μ`, suffices for §2. No such cap-20 or binary `ρ<μ` certificate is claimed in this file. The weight algebra above is an identity, not a proof of that inequality.

## 2. Critical root, conditional on a strict gap

Suppose those inequalities hold with `ρ<μ`. The finite positive matrix is

```
K_c(s) = J(s−R_c)⁻¹ B₀,          B₀g(a,S)=b g(q(D)),
```

with `B₀=0` on dormants. For `s>ρ` its entries are positive sums over explicitly described finite lineages, are continuous, and strictly decrease on each nonzero entry. Strict positivity on every label in the image of `q` follows from consecutive intrinsic star gains to a suitable output degree, then direct production; no copy or dormant is needed for this path. Resolvent bounds are inherited from `W`. Group the full discounted active occupation by direct generations. If `r(K_c(s))<1`, the newborn series converges, all full-state occupation is `O(W)`, and for `s>β` the physical immigration identity gives a finite physical discounted core count. If `r(K_c(s))≥1`, the positive virgin forcing gives divergence in every newborn row; the founder supplies label 2 with positive expected discounted mass. The known physical growth `E C_t≈exp(μt)` forces `r(K_c(μ))=1` and strict monotonicity gives exactly one solution. The permanent-kernel proof identifies its critical root with the same physical exponent, with no similarity or entrywise identity needed.

## 3. Unconditional comparison at and above the physical threshold

For `s>μ` the left side of the positive physical identity for active unit score is integrable. The founder directly emits `J_{q(2)}` over an interval with strictly positive expected discounted mass, so full active occupation from that newborn is finite (Tonelli). Every other newborn label is reached by a positive finite path from it, using star gains followed by direct production under the full process. Hence full and residual active occupation are finite in all newborn rows. Positive separation into direct generations gives

```
J(s−P_c)⁻¹ 1_active = Σ_{j≥0} K_c(s)^j J(s−R_c)⁻¹ 1_active .
```

On newborns the forcing dominates the strictly positive original-particle survival contribution; therefore `r(K_c(s))<1`. All entries at `μ` itself are finite and `r(K_c(μ))≤1`: fix `s₁>μ`, use `(K_c(s)^2)_{ii}≥K_{ij}(s) K_{ji}(s₁)` for `μ<s≤s₁`, and `r(K_c(s))²≥(K_c(s)^2)_{ii}`. Monotone convergence gives the entries at `μ`, then continuity of the finite matrix radius gives the inequality. This does not prove a uniform below-critical gap, nor exclude a physical singularity caused by the residual or immigration instead of the direct-birth matrix. It applies to both finite gauges independently of the Lyapunov attempt.

## 4. Physical shell asymptotic in the required weight

Write `A^p_t(a,S)` for the old full permanent-mark semigroup applied to a bounded mark-independent active terminal score `f`. In this paragraph `S` can also carry arbitrary **compatible** permanent labels on strict descendants, and the top root label is free. The published exponential remainder [almost-sure.md §1](../../fringe/almost-sure.md) gives

```
A^p_t(a,S) = exp(μt) H^p_a(S) ℓ(f) + O(exp(σt) W^p_a(S)), σ<μ,
```

uniformly on such seeds, with a common scalar functional `ℓ`. This extends to an auxiliary root label zero: its first-residence formula is `A^p_t(0,S)=A^p_t(a,S)+a∫₀ᵗ E^0_S A^p_{t-u}(a,S_u)du`. The original intrinsic drift `L⁰ W^p ≤ 0.36 W^p` absorbs the exponential remainder and the missing tail in the leading integral, taking `σ>0.51`. Arbitrary compatible labels, including labels chosen equal to **current** ones in the initial seed, are admissible here: the mark inequality remains valid as degrees rise and the proof uses only that inequality, not the past birth history. The auxiliary root label zero requires no inequality against a parent, since it has no parent send.

Let `M_t(S)=E^0_S Σ_v f(S_v(t))`, including the top core, using the same planted degree convention at each vertex. The permanent transform at top label zero gives exactly `M_t=T A^p_t(0,S)`: apply the physical backward jump identity on this isolated seed, with no uncancelled root kill. Take `D=I−E` on **unmarked** shapes. Choose the compatible initial permanent labels to be `q(D_parent)` on all strict descendants of `S`, and when subtracting `M_t(Y)` reuse the identical labels below `Y` but set its top label to zero. Then all terms strictly below the immediate children cancel. In particular

```
D M_t(S)= A^p_t(0,S)+Σ_Y [A^p_t(q(D_o),Y)−A^p_t(0,Y)].      (4)
```

Use the old weight with `d_*=0.15`, `h(c)=(1,0.85,0.78,0.7,...)` for `c≥1`, and let `Z_c` use the factors `z(D)` evaluated at current degree. With the chosen initial labels, `W^p_0(S)=Z_c(S)/d_*`. Moreover **bound the differences** in (4), rather than their summands independently: their magnitude in the exponential remainder is at most a constant times `a W^p_a(Y) exp(σt)` by the first-residence integral above, and their leading coefficient has the analogous bound. Since `Σ_Y a W^p_a(Y) ≤ Z_c(S)` for `a=q(D_o)`, the full result is

```
D M_t(S)=exp(μt) χ(S) ℓ(f)+O(exp(σt) W_c(0,S)),          (5)
```

and `χ` is bounded in that weight. For the fixed cap, the same conclusion also follows by summing the two remainders separately, using `a≥1/20`. No assumption about an isolated physical exponent or the size of a boundary cancellation on entire trajectories enters (4).

## 5. Current-mark reconstruction and analytic renewal

Let `A^c_t` be the new full response. Its top-zero transform also gives `T A^c_t(0,S)=M_t(S)`, now summing current labels in `T`, by the transfer and dormant-difference identities. The common original-residence forcing gives

```
A^c_t(a,S)=A^c_t(0,S)−a∫₀ᵗ e^(-au) E^0_S A^c_{t-u}(0,S_u)du,
A^c_t(0,S)=D M_t(S)+Σ_Y a_Y∫₀ᵗ e^(-a_Y u) E^0_Y A^c_{t-u}(0,Y_u)du,
```

where `a_Y=q(D_o)` at the initial shape. The resolvent convolution acts in `W_c(0,·)` (not in the permanent-label norm along later intrinsic trajectories). Current path factors decrease on gains. A promotion contributes exactly the old local inequality with `z(D)` evaluated at the birth degree. Therefore `L⁰W_c≤0.36 W_c` for the cap by the original capped weight verification. The convolution norm on functions bounded by `exp(σt) W_c(0,S)` is at most

```
θ_σ = max_{a∈image(q)} (a+d_*)/(σ+a−0.36) < 1   (σ>0.51).
```

Indeed `W_c(0,Y)=(a_Y+d_*)W_c(a_Y,Y)/d_*` and `Σ a_Y W_c(a_Y,Y)≤d_* W_c(0,S)`. This yields a unique solution in the indicated space after subtracting its leading term: the coefficient is obtained by applying the same contraction at `μ` to `χ`; the truncated convolution of this exponential differs from the infinite one by `O(exp(0.36t)W_c)`. The resulting Volterra solution agrees with the actual finite-time response. For example, on each size-killed domain small time-step uniqueness is immediate; for unrestricted shapes one can instead use a norm weighted by vertex count, with `E^0 n(S_u)≤e^u n(S)` and disjoint children. On bounded time intervals the child convolution has a uniform finite kernel bound in that norm, giving factorial iteration uniqueness. Both the obtained solution and the actual response belong to that class by §3 of the construction. Fixed positive `a` and dormant differences inherit the expansion. Thus for each bounded common terminal score,

```
P^c_t f(a,S)=exp(μt) H^c_a(S) ℓ(f)+O(exp(σt) W_c(a,S)).        (6)
```

The same on a dormant state uses the difference weight or any equivalent weight, since the label set is finite. The auxiliary top-zero size bound is the same as for active positive labels. The binary gauge follows by this proof if an old permanent **binary** bounded-observable expansion with a uniform compatible-seed weight is supplied, along with `L⁰ W_c≤q₀ W_c` and `q₀+d_*<μ` for that same weight (a constant-at-high-degree profile must check the jump at the binary seam). The cap-20 proof does not assume those binary inputs.

Now take the 19-column choice `f=B₀`, and write `F(s)=J(s−P_c)⁻¹B₀`. For real `s>μ` the direct-generation expansion gives `I+F(s)=(I−K_c(s))⁻¹`. By (6), `F(s)` extends meromorphically across `μ` with at most a simple pole.

To verify that the pole is present, use `f=B₀1≥0`. The permanent functional `ℓ(f)` is strictly positive from its positive newborn residual integral (each newborn initially has `b=1/2`). The physical identity on the childless newborn gives `A^c_t(0,J)=M_t(J)=A^p_t(0,J)`. The permanent auxiliary label-zero response dominates `A^p_t(a,J)` by its first-residence formula, and the old newborn leading coefficient is strictly positive. Thus the new leading coefficient for this score is strictly positive at `(0,J)`. To pass to every positive label, let `C` be the entire positive emission row on actives. This row is bounded from the weighted score space to the `W_c(0,·)` space: child and grandchild sums are bounded by the first two weighted path layers, whose factors are bounded below, and `r` is bounded whenever the root label drops. By (6), `exp(−μt) C P^c_t f` has a nonnegative limit `G_*`. The original-residence formula and the intrinsic drift bound give the leading coefficient at label `a` as `V_{μ+a}G_*`, including at `a=0`: the initial terminal term grows at most as `exp(0.36t)`, and convolution of the remainder is `O(exp(σt))` in the weighted norm for `σ>0.51`. The same nonnegative integrand yields `V_{μ+a}G_*(J)>0` for every `a` if it is positive at `a=0`. Consequently `F` has a pole. Since `K_c(μ)` exists by §3, `r(K_c(μ))=1`; otherwise `F(s)` would stay bounded.

To obtain the gap, `(I+F)⁻¹` is meromorphic near `μ` by the finite adjugate formula and invertibility for real `s>μ`. Its limit from the right exists and is finite by §3, so any pole is removable. Thus `K_c` continues analytically across `μ`. Positivity promotes analytic continuation into actual convergence to the left: choose a real `s_1>μ` close enough that the analytic disk around it reaches some `s_0<μ`. Taylor coefficients of `K_c(s_1-x)` are nonnegative entrywise and given by integrated chronological powers. Summing at positive `x=s_1-s_0` and using Tonelli proves finiteness of the original Laplace integrals at `s_0`. For all real `s>s_0`, actual kernel entries strictly decrease, from the explicit positive direct paths, so `μ` is its unique critical root in that interval.

## 6. Full-state active-count bound below criticality

Let `A_s(S)` denote the Laplace transform in (6) with `f=1_active`, and `j_s=J A_s`. Initially take `s>μ`. Splitting at the first direct production gives

```
A_s = V_s + (s−R_c)⁻¹ B₀ j_s,   V_s=(s−R_c)⁻¹ 1_active.
```

For the vector rate observable, `P_s B₀ = (s−R_c)⁻¹B₀ (I+F)` in Laplace notation. Hence

```
V_s=A_s−(P_s B₀)(I−K_c(s))j_s .
```

All quantities on the right are meromorphic near `μ` with bounds in the weighted supremum norm. The singular terms cancel: both full transforms have the common left functional `ℓ` and leading state coefficient `H^c`, and `F(I−K_c)=K_c` shows that `(I−K_c)j_s` is analytic there. For an explicit cancellation argument, write the residue of `j_s` as `g ℓ(1)` with `g=JH^c`. The residue of `F` is `g ℓ(B₀)`. Since `(I−K_c)F=K_c`, `(I−K_c)g` vanishes at `μ` and
`ℓ(B₀)[(I−K_c)j_s]_{s=μ}=ℓ(1)` follows from `F(I−K_c)j_s=K_c j_s` by taking residues and `K_c(μ)g=g`. Thus `V_s` is analytic across `μ` in the weighted norm. Apply the same positive Taylor proof to obtain `V_{s_0}(S)≤C W_c(S)` uniformly, with a possibly closer `s_0<μ`. Finally, after elapsed time `t`, each active particle survives as an active original for an additional time `u` with probability at least `exp(−u/2)`. Integrating over `0≤u≤1` in the occupation tail gives

```
R^c_t 1_active(S) ≤ C' exp(s_0 t) W_c(S).
```

For the fixed finite cap, the same bound holds for the combined active and raw dormant count. A dormant present at time `t` switches by time `t+1` with probability at least `min_(l,u) (u−l)(1−exp(−u))/u>0`; intrinsic jumps neither change the two labels nor remove the particle. After such a switch the resulting original is active throughout `[t+1,t+2]` with conditional probability at least `exp(−1)`. Each particle at `t` has a distinct resulting original. The discounted active occupation on that interval therefore bounds a fixed positive multiple of `exp(−s₀t)` times the expected dormant count. This bound is specific to a finite set of positive drops. The proof does not supply explicit `s₀` or a pointwise test for `R_cW_c≤s₀ W_c`.
