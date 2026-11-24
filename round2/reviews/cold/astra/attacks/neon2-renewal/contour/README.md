# Neon-2 Claim 3, part 3: contour provenance and replacement proof

**As written: GAP in contour independence. With the replacements below: HOLDS under the listed inputs.** This qualifies Astra's earlier conditional HOLDS: that review took the cited estimates as inputs. Fable is right that the printed proof imports a contour remainder. The dependency is through first moments, so it neither makes the proof circular in the almost-sure limit nor refutes that limit. This note concerns part 3 only.

The source is `ft2-neon2@ac55a7b9fa9236e73a305c32fbf509ae3bf028a3:round2/neon2/claim14-proofs/README.md`, blob `d486bf4901a54f6a41c6818c9000feed5c980cf2`. The source paths and line references below are at that commit. The triggering report is Fable's `attacks/fable-on-neon2-3.md`, §§1–2.

| Use in the pinned proof | Provenance as written | What is needed |
|---|---|---|
| Step 3.3(2), lines 236–244: variance (13) | `open/limit-law/README.md` §1 bounds private branches using the free cell bound from `fringe/spatial-limit.md` (3). That equation inserts `fringe/almost-sure.md` (3), obtained by Fourier inversion and rectangular contour shift after subtracting the first kernel term. The gradient there is stated for the core count. | A summable bound for each **fixed** lag, for every cumulative count. |
| Theorem 2(3), line 130: affine lower comparison | `standalone/claim-6.md` §4 gives positive eligible-core values using local scripts. The value bounds are independent of contours, but the displayed route from these values to physical means goes through its §3 equation (3) and §4 equation (4), with the finite-matrix contour shift. | A lower bound on the limiting free mean amplitude, then a real Abelian limit. |
| Lines 7, 58, 76, 128, 146 and 235: signed normalization and slope | `influence/zero/critical.md` §1 gets the vanishing boundary term from the uniform cell bound, including its inert extension. §3 also cites the spatial cell limit for the time-domain amplitude. Both go back through spatial equation (3) to the same contour estimate. The finite-tree transfer identity and contraction do not use it. Section 2 already has a positive scalar Laplace decomposition; it does not require a contour to make that decomposition. | The inert virgin's sharp mean bound and positive limiting amplitude, then stopped Dynkin and a real Abelian limit. No two-sided quantitative inverse estimate is needed. |

There is one more place to bypass in part 3: line 245 invokes Theorem 1(2) to name the degree tail. Theorem 1 Route B cites the contour-derived spatial limit, and Route A needs a projected balance on the full bounded intrinsic domain that is not supplied just by the cylinder-domain statement in `operator.md` §4. Part 3 needs only the degree tails, not that full stationary fringe-law theorem. Section 5 below also supplies contour-free identification with the **unique** spatial eigenmeasure, so uniqueness need not be imported from the old remainder.

## 1. Run the positive scalar renewal first

Keep as inputs the physical mean exponent and its bracket (`proofs/existence.md`, `proofs/bracket-motifs.md` and the refined physical bracket used in Claim 14); the capped positive transfer, weight and critical Perron equation in `fringe/finite-kernel.md` and `standalone/claim-6.md` **§2**; and the degree-two trigger comparison and chronological bound of `finite/gauge/physical-trigger.md` §§1, 4. Also keep the unlimited signed occupation contraction in `influence/zero/pair/README.md` and the **reciprocal-stock separation bound alone** in `open/forest-age/inert/README.md` §2. The arithmetic and certificates of these inputs are inherited, not re-certified here. Write `q₀=9/25`, and choose the certified chronological exponent `ρ'<μ` for the untriggered semigroup. This bound controls a different positive process; it is not a spectral remainder for the full renewal.

Write `j=J_(1/2)`, `h_tr=b·1_(D_o=2)`, let `P_t` be the full positive mean semigroup and `Ĥ_t` the untriggered one. Use `B` for a bound in `Ĥ_t W≤B e^(ρ't)W`. Define

```text
p(u) = (Ĥ_u h_tr)(j),          K(du) = e^(-μu)p(u)du,
m = ∫ u K(du),               a(S) = ∫₀∞ e^(-μu)(Ĥ_u h_tr)(S)du.
```

The first residence of `j` contributes `p(u)≥(1/2)e^(-2u)` (intrinsic rate `3/2`, kill rate `1/2`). The positive comparison gives `F(μ)=∫K=1`; the chronological bound gives `0<m<∞` and `a(S)≤BW(S)/(μ−ρ')`. For a score `f` with `|f|≤B_f W` and `|f(S)|≤L_f |S|`, let `z_f(t)=e^(-μt)(Ĥ_t f)(j)`. It is continuous: if `V` is the sum of ordinary sizes over living auxiliary particles, positive jumps have size at most `V` and their rate times added size sums to at most `3V`. Indeed the intrinsic term costs at most `V`; the old-child copy arithmetics use disjoint child subtrees and total rate at most `1/2` per child, costing at most `3V/4` even with the extra leaf; direct untriggered births cost at most `V`. Thus stopped `G V²≤9V²`, giving a uniform second-moment bound on compact time intervals. Nonexplosion and absence of deterministic-time jumps give continuity by uniform integrability. The exponential envelope from `Ĥ_t W` then makes `z_f` directly Riemann integrable, including its positive and negative parts.

First-trigger decomposition on lineages gives `Z_f=z_f+K*Z_f`, with `Z_f(t)=e^(-μt)(P_t f)(j)`. Iteration gives the renewal solution: on every compact interval `Z_f` is bounded by the full weight bound in Claim 6 §2, and the remainder is bounded by that supremum times the probability that the sum of `n` independent `K`-increments is within the interval, which tends to zero. Real key renewal gives `Z_f→β_f=(∫₀∞z_f)/m`, and `Z_f` is bounded on the whole half-line. In particular `β_1>0` by the original particle's residence at `j`.

For a general marked seed use `q_S(u)=e^(-μu)(Ĥ_u h_tr)(S)`, a **discounted expected intensity**, not a probability. First-trigger decomposition gives the untriggered term plus `∫₀ᵗ q_S(u) Z_f(t-u)du`. After subtracting `a(S)β_f`, the portion `u≤t/2` is bounded by `a(S) sup_(v≥t/2)|Z_f(v)−β_f|`. The rest, including the omitted tail past `t`, is bounded by a constant times `e^(-(μ−ρ')t/2)W(S)`; the untriggered term itself has the exponential envelope. Hence

```text
| e^(-μt) P_t f(S) − β_f a(S) | ≤ ε_f(t) W(S),
sup_t ε_f(t)<∞,   ε_f(t)→0.                                  (A)
```

For a finite free tree, orient from a core founder `o`, initialize **all** nonfounders with compatible marks (all `1/2` suffice), and use Claim 6 §2 equation (1). Let

```text
A(T) = Σ_(v≠o) a(S_v) + ∫₀∞ e^(-μu) I_(T,u) a du.
```

The input weight estimates are `Σ_(v≠o)W(S_v)≤K₁C(T)` and `I_(T,u)W≤K_I C(T)e^(q₀u)`. The first sums the ancestor weights using the ratio `10/13`; the second is the founder drift bound from §2. Apply (A) to initial fringes; apply it under the immigration integral, splitting at `t/2` with discounted envelope `K_I C(T)e^(-(μ−q₀)u)`. The resulting errors are bounded by `C(T)` times a sum of `K₁ ε_f(t)`, a constant times `sup_(v≥t/2) ε_f(v)+e^(-(μ−q₀)t/2)`, and the omitted immigration tail `O(e^(-(μ−q₀)t))`. The founder contributes at most `e^(-μt)` for a cumulative count. Thus for `f_k=1_(D_o≥k)`, including `f_2=1`,

```text
e^(-μt) E_T N_(≥k)(T_t) = β_(f_k) A(T) + o_k(1) C(T)          (B)
```

uniformly in all free seeds. In particular the free mean bound for the resolvent follows already at this stage.

## 2. Isolated virgin and signed critical facts

In the **physical inert** virgin `J`, keep its slot in the degree and suppress its arrows. The same marked cancellation gives the inert core mean `m_J(t)=1+∫₀ᵗ I_(J,u)^0 P_(t-u)1 du`. There are no initial nonfounders. The top drift still uses the true core degree including the slot, and only births of real cores change the path weight, so `I_(J,u)^0 W≤1/d+2e^(q₀u)`, `d=3/20`. The cancelled terms are exactly the artificial nonfounder kills and missing incoming parent arrows; no kill or incoming arrow is included at the top. Finite-time stopped transfer uses the same polynomial size estimates as Claim 6 §2, with inert size dominated by ordinary Yule size. Applying (A) under this integral as in (B) gives

```text
e^(-μt)m_J(t) → Ψ⁰(J) = β_1 ∫₀∞ e^(-μu) I_(J,u)^0 a du > 0.
```

For strict positivity, with probability at least `e^(-2u)` neither original clock rang before `u`; then the direct term in `I_(J,u)^0 a` pays `(1/2)a(j)=1/2`. The normalized mean is bounded for all `t`. This precedes any signed identities and any physical pathwise limit.

Now use the inherited contraction to define the signed `U_s`, and write `𝒯` for the sum over all core fringes in a planted tree. The exact signed transfer (Track 3 (5)) gives `(s−L⁰) 𝒯 U_s=(1−U_s(J)) L⁰ C`: incoming arrows are assigned to parents and promotions add a virgin score. For real `s>μ` the boundary term in stopped Dynkin starting from `J` vanishes, since `|𝒯 U_s|≤||U_s||∞ C` and `m_J(t)=O(e^(μt))`. With `M_s(J)=∫₀∞e^(-st)m_J(t)dt`, integration by parts on the core drift gives

```text
U_s(J) = (1−U_s(J))(s M_s(J)−1),
s M_s(J)(1−U_s(J)) = 1.                                      (C)
```

For any bounded function `g(t)` tending to `g∞`, the real Abelian identity `r∫₀∞e^(-rt)g(t)dt→g∞` follows by the substitution `x=rt` and dominated convergence. Apply it to `e^(-μt)m_J(t)`. Then `M_s(J)→∞`, `(s−μ)M_s(J)→Ψ⁰(J)`, and (C), continuity and differentiability of the contracted resolvent imply `U_μ(J)=1` and `−U'_μ(J)=1/(μ Ψ⁰(J))`. This recovers all the critical facts used in part 3.

The pinned root-independence and coboundary identities can now be used unchanged: their equality on `s>1` uses Yule domination and extends by analyticity, with no inversion. For a proof confined to real parameters, these are real-analytic identities on `(ℓ,∞)`: the real Laplace integrals have locally norm-convergent power series (integrate the exponential series with a strictly positive discount margin), and the inverse at each real center has a locally convergent Neumann expansion of its perturbation. Equality on `(1,∞)` extends along the connected interval by uniqueness of real power series. For real `s>μ` the free boundary term vanishes by (B), and (C) gives `1−U_s(J)>0`. Theorem 2(3)'s occupation identity and (B) therefore give

```text
Φ(T) = lim_(s↓μ) Φ_s(T) = β_1 A(T)/Ψ⁰(J) ≥ 0.                (D)
```

For the inert martingale at arbitrary planted seeds, the pinned positive-probability reachability from `J` extends its mean bound to each fixed seed. The same stopped calculation proves `𝒯 U_s≥0` for `s>μ`; the limit is nonnegative, is harmonic at `μ` by transfer, and takes value one at `J`. The pinned grafted jump bound uses just the free root-independence and the bounded root toll.

## 3. Lower comparison without a finite-time lower estimate

Use a new symbol `H` for the **positive** permanent-mark value (unrelated to the signed gradient). Let `R_t` be the residual semigroup excluding all direct newborns, `B₀g(S)=b(S)g_(m(D_o))`, and `g>0` the critical right Perron vector of the finite positive kernel. Define `H=∫₀∞e^(-μt)R_t B₀g dt`, so `H(J_i)=g_i` and `H=O(W)` by the residual growth bound. Split the direct payment into nontriggers and triggers. The resolvent equation for `H`, iterated over nontrigger direct births, gives `H=g₂∫₀∞e^(-μt)Ĥ_t h_tr dt=g₂ a`. More explicitly, at the newborn rows the residual-to-nontrigger matrix is `K_μ` with column two zeroed; it has radius below one by the trigger comparison, so the renewal equation is solved by its convergent generation series. The same first residual row at any `S` has finite toll by `R_tW≤e^(0.51t)W`; hence the newborn identity extends to all `S` by that row's first-direct-birth decomposition. No contour coefficient is being identified here.

For `C≥2` the eligible-core count from Claim 6 §4 is at least `(C+5)/4`. A nonfounder with a real leaf and `c≤2` pays at least `e^(-(μ+1/2))(1−e^(-1)) min_i g_i/3` in `H`: keep the original residual particle alive up to time one, and use its first root ring in that unit. Before that ring the root still has a leaf and the same core degree; the conditional promotion probability is at least `1/3`. The direct integral is the compensator of the promotion toll, so dropping all other payments is valid. A bare vertex with a bare child obtains a leaf in one unit by silencing the root and the child's other neighbour and requiring the child's first ring to choose the root. The child's degree is two before its first ring on this event. Independent survival of the original particle and the first script in the next unit give another universal positive lower bound. A bare eligible vertex with no bare child can be charged to its bare parent with a bare child; at most two charges go to each such parent. After excluding the founder there are at least `C/12−1` eligible nonfounders covered by these scripts. They have `a≥c_*>0`, uniformly in marks and geometry, since `H=g₂ a`.

All terms defining `A(T)` are nonnegative. Thus `A(T)≥c_* (C(T)/12−1)_+` for `C≥2`, and `A(T)≥0` also for `C=1`. Equations (B) and (D) give `C≤K₄(1+Φ)` for some universal `K₄`. The printed two-unit leaf-stock script then upgrades this to strict linear comparability using harmonicity: from every seed one can reach `C≥ceil(2K₄)` in two units with a universal positive probability, and at that size `Φ≥1`. With `Φ≥c_min>0` this gives `a₀ C≤Φ≤b₀ C`. This argument bounds the **limit** of the normalized mean; it does not assert the printed uniform finite-time lower bound from an estimate of the amplitude alone.

## 4. Fixed-lag sampling and pathwise conclusion

Couple two physical trees differing by one attached leaf. Pair clean cores, leaves and incident slots, syncing their clocks and choices; mark both copies of the attachment core dirty, or both copies of its old parent and the newly promoted core if the attachment target was a leaf. A dirty core initiates at rate one and can create one dirty core or infect a clean pair; other rings do not grow the dirty set. Clean degree counts cancel, even when a clean core sends a leaf to a dirty recipient, because the clean incident slots remain matched. Size localization and Grönwall give expected combined dirty count `≤3e^(2u)` over elapsed `u`. Therefore the backward conditional mean of the terminal `N_(≥k)(T_t)` has one-birth jump at time `v` at most `3e^(2(t−v))`, uniformly in the state and `k≥2`.

On `[t−h,t]`, the localized backward martingale is purely jump with compensator, with total jump rate `|T_v|`. Its mean-square bound passes to the limit by finite-time Yule moments. Conditional jump isometry and `E[|T_(t−h+w)| | F_(t−h)]=|T_(t−h)|e^w` give the replacement for (13):

```text
Var(e^(-μt) N_(≥k)(T_t) | F_(t−h))
 ≤ 9 e^(-2μt)|T_(t−h)| ∫₀ʰ e^(4(h−w)+w) dw.
```

The expectation is `O_h(e^(-(2μ−1)t))`, summable on every fixed mesh. On the countable intersection over fixed integer lags, mesh denominators and degree thresholds, the normalized conditional-mean residual tends to zero by Chebyshev and Borel–Cantelli. From (B), (D) and affine comparison the normalized error between the conditional mean and `d_k e^(-μ(t−h))Φ(T_(t−h))` is at most `o_k(1)` in `h` times `K₄[e^(-μ(t−h))+e^(-μ(t−h))Φ(T_(t−h))]`, where `d_k=(β_(f_k)/β_1)Ψ⁰(J)`. Steps 3.1–3.2 give the positive martingale limit `M∞` by the bounded jumps, `2μ>1`, the stopped reciprocal-stock bound, and the uniform restart argument. First let `t` tend to infinity along the mesh and only then let `h→∞`. Interpolate the unscaled monotone cumulative counts with factors `e^(±μ/m)` for a mesh of denominator `m→∞`. This yields `e^(-μt)N_(≥k)(T_t)→d_k M∞` at all real times; exact counts follow by subtraction.

The remaining pinned discrete argument and leaf relocation require only these limits, linear comparability, `μ<1`, and the free jump bound, as checked in [the earlier pathwise review](../real.md). In particular the limit at size `n` for the normalized core count is `X=Ψ⁰(J)Y∞`; the Gamma normalization and Doob's maximum dominate `C_n/n^μ` by an `L²` variable. The same domination applies to cumulative counts, giving `E N_(≥k)(n)/n^μ→(β_(f_k)/β_1) E X` and hence the ratio `c_k/c₂` directly. The fixed-size relocation argument gives support `[0,∞)`.

## 5. Spatial degree identification without the old remainder

The spatial operator construction and existence of a probability eigenmeasure in `operator.md` §§1–2 use old-edge localization and logarithmic mean growth only. To obtain uniqueness without the contour, reuse the **pre-remainder** cell identity and drift bound of `fringe/spatial-limit.md` §1, equations (2) and preceding text. They use only positive transfer; even with infinitely many old cells the given cell has finitely many old neighbours, and its ordinary size has birth intensity at most its size plus half the old core degree. This gives the local polynomial moments for stopped transfer. The weighted immigration includes only new children of the base root and satisfies `I_u W≤1/d+2e^(q₀u)` uniformly in the exterior. There are no initial nonfounders in this cell. Consequently the same convolution split with (A) gives for the mean sum over its new members

```text
e^(-μt) M_t f → β_f ∫₀∞ e^(-μu) I_u a du,
```

uniformly in all starting environments for each score in §1. In particular `e^(-μt)Q_t1→h`, where `h=β_1∫₀∞e^(-μu)I_u a du` is bounded. It is positive at each input: a root with a leaf can promote with positive probability on a bounded interval, and a bare root can first be sent one leaf by a chosen old neighbour. The latter event has positive probability by the finite-neighbour ordered-proposal construction (silence arrows into that neighbour until it sends). Each direct newborn has `a(J_i)>0`, by immediate virgin residence.

For any bounded radius-`r` cylinder `f`, use the same depth-`r` re-rooting score `F_(r,f)` as in `fringe/spatial-limit.md` §2. It sums the readings at depth `r` in a marked planted fringe; the ancestor's exterior slot preserves boundary degree marks. It is bounded both by `4^r W||f||∞` and by `|S| ||f||∞`, so (A) applies. Each cell member of depth greater than `r` is counted once by its depth-`r` ancestor below the base root, and the omitted shallow members cost in expectation at most `||f||∞ Σ_(j=0)^r t^j/j!` by the rate-one promotion bound. Hence `e^(-μt)Q_t f→h β_(F_(r,f))/β_1` uniformly. Integrate the `f=1` convergence against any probability eigenmeasure of `Q_1` at `e^μ` along integer times to get `π(h)=1`. Integrate the cylinder convergence along the same times to get `π(f)=β_(F_(r,f))/β_1`. Cylinders determine the probability, proving uniqueness. For `f=1_(D_o≥k)` take `r=0`, so `π(D_o≥k)=β_(f_k)/β_1=c_k/c₂`. The independent local coboundary Cesàro identity also identifies this ratio with `𝓜_μ(ψ_k)`; the stationary fringe-law identification on **all** bounded scores is not needed for this proof.

Thus all three reported imports and the degree-tail identification have replacements in the dependency chain for part 3. The convergence estimate used for the full normalized means is `o(1)`, not the old exponential remainder. The lag in the variance proof is fixed before the time limit. No revision of the claimant's text or board entry is implied by this note.
