# Long's inert generations through Astra's 19-label kernel

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 16: GPT proof review; Claude simulations; inputs include Long’s generation sums, the capped compact kernel and critical tilt. See the [Astra ledger](../../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


**Result.** C. D. Long's actual discounted inert generation coefficients satisfy, locally uniformly for real `s` near the physical exponent `μ`,

```
A_m^∞(s) = c(s) ζ(s)^(-m) + O(R^(-m)),     R>sup ζ(s),
r_∞(s)=1/ζ(s),     c(s)>0,
```

where `ζ(μ)=1`, and `ζ` is the simple positive root of `det(I−K(s,z))` in the depth variable. In particular the critical coefficients converge exponentially to a strictly positive constant. The root and amplitude are real analytic; the ratio and its first discount derivative converge exponentially as well. The proof uses the original permanent marks, with no current-mark or signed-zero-gauge numerical inputs.

The coefficients and their pressure/crossing are Long's: see the pinned [introduction, `main:computation`](https://github.com/long-mathematics/friend-tree-growth/blob/5bd25a7d506cb64a8444b86334ab531d295f6c78/sections/original/introduction.tex#L28-L72) and [pressure theorem, `gc:crossing`](https://github.com/long-mathematics/friend-tree-growth/blob/5bd25a7d506cb64a8444b86334ab531d295f6c78/sections/original/pressure.tex#L1-L49). The earliest recorded statement in the compared histories is Long's [`090bbbe`](https://github.com/long-mathematics/friend-tree-growth/commit/090bbbe5faec83ec4124bd427eb5fda32c79bef1); see the source comparison for the scope of that chronology. The connection here is to the separate [Astra 19-label kernel](../../../fringe/finite-kernel.md), not a similarity, compactness claim, or finite-rank approximation for Long's full-environment operator.

## 1. Exact depth transfer, including the inert founder

For a planted fringe write `L⁰` for intrinsic evolution, `a` for its permanent root birth mark, `Q` for the two positive old-child copy terms and `B=B₀J_mark` for direct newborn production. Here `B₀g(S)=b_S g(1/min(D_o,20))` and `J_mark` evaluates a score at the 19 positive-mark newborns. Define

```
R_z = L⁰−a+z Q,                 P_z = R_z+z B,
T_z f = Σ_v z^depth(v) f(S_v).
```

Root the inert tree at its initial core. Its inert slot counts in `D` and in the core degree, sends nothing, and selection of the slot is null. Attach permanent birth marks below the root and attach the auxiliary mark **zero** to the root only. For any score `f`, the pointwise finite-row identity is

```
L⁰ T_z f = T_z P_z f.
```

At a physical child `Y`, the gain missing from its own `L⁰` is the send from its parent, with rate `p=1/D_parent`. The parent's pair of copy terms is `a_Y f(Y)+p[f(Y⁺)−f(Y)]`; its factor `z` puts both parts at the child's depth. The first part cancels the child's killing exactly. The new fringe at a promotion has factor `z` relative to its promoter. An intrinsic jump changes the shape read by existing summands, without advancing their depth: even an intrinsic promotion is **not** multiplied by `z`. It is the separate direct newborn emission that accounts for the newly created summand. The only kill with no cancelling parent term is the top kill, which is zero here. For an artificially nonzero top mark, `T_z P_z f=L⁰ T_z f−a f(S)` instead. No ordinary ancestral immigration or external input is to be added.

In particular for the childless inert newborn `J_0`, `T_z f(J)=f(J_0)`. Thus its mean weighted physical count is `(P_t^z 1)(J_0)`, with `P_t^z` the full weighted first-moment evolution. For a general initially planted marked seed it is `T_z P_t^z1`, with its initial old fringes included. If the top is excluded instead, from `J` one must use immigration of **both** kinds, `z Q+z B`, with the top degree including the slot; add the top unit score separately. This is equivalent by first-residence conditioning at mark zero. For an actual free founder the same immigration formula uses its degree **without** the slot and counts the free founder separately. No claim is made here about retagging a newborn in an arbitrary active exterior.

For justification on real `z≥0` on any bounded time interval, use the positive marked branching family and the usual stopped backward jump cancellation. The original linear size and path-weight estimates give finite first moments also with arbitrary fixed finite emission multiplier `z`. In physical sums one cannot crudely use exponential weight in the number of vertices at all horizons. Instead for any `r≥max(1,|z|)` use `X_t=Σ_v r^depth(v)` on physical cores: the drift is at most `r X_t`, since each core promotes at rate at most one. For its square, a new summand is at most `r X_t`; the square drift is at most `(2r+r²)X_t²` using the individual parent's weight. Localization and Cauchy–Schwarz with the ordinary Yule size moments control polynomial size factors in sums of local scores/derivatives; for uniform integrability one can use higher integer moments of `X_t`, by the same drift bound with a moment-dependent constant. This proves the finite-time identity. Alternatively, use the identity for `0≤z<1` and identify the analytic germs. Only the nonnegative half-line is needed for the initial expectation/series equality below.

## 2. Holomorphic resolvent and the actual coefficient series

Use `W` from the permanent construction, including at the auxiliary zero root. It satisfies

```
QW=(a+.15)W−h_weight(c_o) ≤ .65 W,
R_1 W≤ρ W,       ρ=.51,
R_r W ≤ [ρ+.65(r−1)_+]W       (r≥0).
```

In particular the intrinsic part is unchanged even on the auxiliary row. The mark zero is only a starting type: all children and marked newborn outputs still have positive marks. Expand the residual Green function in time-ordered copies about `L⁰−a`. The bounded copy operator on the weighted score space gives absolute domination by the positive series at `|z|`; the displayed drift gives its integral norm bound `1/[Re s−ρ−.65(|z|−1)_+]`. This proves joint holomorphy of `G(s,z)=(s−R_z)⁻¹` on

```
Re s > .51 + .65 (|z|−1)_+.
```

For example `Re s>.549` admits `|z|<1.06`. Holomorphic here refers to the bounded resolvent on the weighted score space, or equivalently to its absolutely converging time-ordered series; no strong continuity assumption on a supremum space is required. Set

```
K(s,z)=z J_mark G(s,z) B₀,
α= (G(s,z)1)(J_0),   w= z (G(s,z) B₀)(J_0),
v= J_mark G(s,z)1.
```

All are jointly holomorphic on the same region. For real `s` near `μ` and sufficiently small positive `z`, the positive direct-generation sum converges since `K(s,0)=0`. Tonelli, the depth transfer, and integration of each core's presence from its physical promotion time give the germ identity

```
M_s(z) := Σ_(m≥0) A_m^∞(s) z^m
       = s ∫₀∞ e^(−st) E⁰_J Σ_v z^depth(v) dt
       = s [α + w (I−K(s,z))⁻¹ v].                       (1)
```

The physical sum runs over cores present at time `t`. For example `M_s(0)=1`, since the zero-mark particle is neither killed nor copied at `z=0`. The right-hand side is a meromorphic continuation of the **physical generation germ**, including for discounts below `μ` for which the defining integral at `z=1` diverges. It is not a formula for the expectation integral past its convergence circle. Complex `s` near `μ` is also valid in the germ by domination and the identity theorem: the expected depth-`m` count at time `t` is at most `t^m/m!` by the rate-one conditional promotion bound, hence `|A_m^∞(s)|≤(Re s)^(-m)` by bounding the discounted birth measure at real discount `Re s`. This gives a uniform initial power-series disk on compact discount sets.

## 3. Simple pole and aperiodicity in **depth**

For real `s` and real `z>0` in the stated domain, `K` is entrywise strictly positive and has a series `Σ_(j≥1) z^j K_j(s)` with nonnegative coefficients. The coefficient `K_1(s)` is **strictly positive entrywise**: the original marked virgin can gain leaves until a requested output degree is reached and then emit a direct newborn at nonzero rate, all with no prior copy on that lineage. Thus the Perron eigenvalue `λ(s,z)` is simple at `(μ,1)` and

```
λ(μ,1)=1,    λ_z=ℓ K_z g>0,    λ_s=ℓ K_s g<0,    ℓg=1.
```

The last inequality follows by differentiating the Laplace weights in chronological time. Here subscripts on `K` in derivative formulas denote partial derivatives, not discounted kernels at distinct arguments. The analytic implicit function theorem gives `ζ(s)` with `λ(s,ζ(s))=1`, `ζ(μ)=1`, with `ζ(s)>0` for nearby real `s`. For `|z|<ζ(s)` entrywise comparison to `K(s,|z|)` rules out eigenvalue one. On `|z|=ζ(s)`, if `K(s,z)u=u` for a nonzero vector, the triangle inequalities and the positive left Perron vector at radius `ζ(s)` force equality throughout. Thus `|u|` is a positive Perron vector. The `j=1` terms imply that for each row `i`, all phases `z u_j` agree with `u_i`; the diagonal `j=1` term then forces `z/|z|=1`. This is depth aperiodicity, not a time Fourier assertion. Positivity of `K(s,z)` alone would not prove it without the precise powers; e.g. replacing every depth charge by `z²` introduces poles at both signs.

At the positive root the inverse matrix has simple pole with singular part `gℓ/[λ_z(s,ζ)(ζ−z)]`, taking the eigenvectors at `(s,ζ(s))`. Both `w g` and `ℓ v` are strictly positive, since the auxiliary virgin can directly emit and ordinary originals can survive positive time. Therefore (1) has the same pole, with positive coefficient in `1/(1−z/ζ)` equal to

```
c(s)= s (w g)(ℓ v) / [ζ(s) λ_z(s,ζ(s))].               (2)
```

In particular the physical radius is exactly `ζ(s)`, not the residual resolvent radius. Isolate the pole at `(μ,1)`. The determinant has no other zero on the compact disk `|z|≤1`, and its zero at one is simple. Choose `R>1` slightly larger with no other zeros on `|z|≤R` at `μ`, and still inside the holomorphic-resolvent domain. Compactness away from a small pole disk, and the simple zero/argument principle on that disk, give the same statement for complex `s` in a sufficiently small neighbourhood of `μ`, with `|ζ(s)|<R`. After subtraction of (2) times `1/(1−z/ζ)` the remainder is jointly holomorphic through `z=ζ(s)` and bounded on `|z|=R` uniformly on a smaller discount neighbourhood. Cauchy gives the asserted error `O(R^(-m))`. This proves the critical limit and the actual pressure formula locally (without using Long's general generation-limit theorem for existence in this neighbourhood). This argument gives an exponential rate but no numerical `R` or certified endpoint window for it.

There is also a quantitative angular control, independent of the simple-pole argument near the positive axis. Every input row has one-step and two-step depth paths to the degree-three label: for the former gain a leaf and then emit, and for the latter intrinsically promote in the virgin, take the plus copy of its core child, and emit from this degree-three copy. At criticality the respective contributions are at least

```
L₁ = (2/3)/[(μ+2)(μ+1/2+8/3)],
L₂ = (1/2)(1/2)(2/3)/[(μ+2)(μ+1/2+5/2)(μ+1/2+8/3)].
```

These are lower bounds using maximal starting mark `1/2`: the intrinsic total rates of the virgin, two-core chain with one leaf, and two-leaf star are `3/2,5/2,8/3`. After Perron conjugation at `(μ,1)`, each row has total weight one and contains the two contributions in the same output column, each of mass at least `d_* = min(L₁,L₂) min_i(g_3/g_i)>0`. At angle `φ`, remove `d_*` from each of these two terms and apply the triangle inequality to the rest. The matrix norm after conjugation is at most `1−2d_*[1−|cos(φ/2)|]`. This only controls angles bounded away from zero; by itself it is not a radius for the coefficient remainder.

## 4. Derivatives and comparisons

On a small real neighbourhood where `c(s)` stays bounded away from zero, write `θ=max ζ(s)/R<1`. Then the physical ratio satisfies

```
A_(m+1)^∞(s)/A_m^∞(s) = r_∞(s) + O(θ^m),
ζ'(s)=−λ_s/λ_z >0,
r_∞'(s)=λ_s/[ζ(s)^2 λ_z] <0.
```

Derivatives of the ratio also converge exponentially: use its analytic expansion on a slightly larger complex discount disk where `c` has no zero and `max |ζ|/R<1`, and apply Cauchy on a smaller disk. Equivalently `(∂_s log A_m^∞)/m → ∂_s log r_∞`, with the amplitude term `c'/[m c]` and an exponentially small remainder. Both partial derivatives admit positive-mass formulas at the critical point. A chronological path contributing to `K` uses one mandatory direct emission and `j−1` copies, with total discounted delay `t`. Weight the critical matrix path measures with `ℓ_i g_j` (here `i,j` are the input/output labels). Their mass is `ℓKg=1`; if `N` is the depth exponent, then `λ_z=E N≥1` and `−λ_s=E t>0`. Thus `−r_∞'(μ)=E t/E N`. This distinguishes physical depth increments from the number of direct generations in the matrix inverse. Both moments are finite by holomorphy and can be compared by nonnegative path sums.

An additional normalization check uses the physical time residue `Ψ⁰(J)=lim_(t→∞) e^(−μt)E⁰_J C_t` of the isolated newborn (see the [isolated-seed estimate](../../zero/critical.md)). At `z=1`, formula (1) and the simple time pole give `μ Ψ⁰(J)=μ(wg)(ℓv)/(−λ_s)`. Consequently

```
c(μ)=μ Ψ⁰(J) ζ'(μ)=−μ Ψ⁰(J) r_∞'(μ).
```

That residue is a cross-check, not an input to the depth pole. No assertion is made about Long's whole-environment operator radius away from criticality, global analyticity of `r_∞` on all `s>0`, or an efficient evaluation of the 19 kernel entries.

## Short exact replay

```
taskset -c 8 timeout 8s python3 -B round2/influence/generation/analytic/check.py
```

The fixture lists real arrows on labelled vertices, independently of the residual copy/birth split. It checks 154 planted and 437 free-founder states at `z=7/5`, with a shape- and permanent-mark-sensitive test score. It rejects omitting either depth charge, charging intrinsic changes, killing the inert top, cancelling by current rather than permanent child marks, and giving the free founder a parent slot. These are finite checks of the transfer and boundary conventions; the analytic proof gives the infinite coefficient result.
