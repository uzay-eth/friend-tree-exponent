# Critical characterization of the signed zero gauge

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 14: Claude for transfer/scalar identities; GPT for full contraction and sign transfer; inputs include Claim 12 and the cap-48 continuation certificate. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Let `l=.565658`, `u=.565798`, and let `J` be a virgin with one leaf and one silent parent slot. For the **unlimited** bounded signed solution, on the entire real half-line `s≥l`,

```
U_s(J)>1  iff  s<μ;     U_s(J)=1  iff  s=μ;     U_s(J)<1  iff  s>μ.
```

On `[l,u]` one has `|s−μ| ≤ 471522 |1−U_s(J)|`, and at criticality `−U'_μ(J)>1/845`. These are conservative analytic constants. They use the prior physical bracket only for the quantitative positive-matrix estimates and interval constants below; they do not assume an endpoint sign for the signed evaluator. The bounds apply in particular at both endpoints of the prior bracket, and on any smaller target window inside `[l,u]`.

## 1. Analytic signed solution and physical transfer

On bounded scores on unmarked finite planted shapes, set `V_z=(z−L⁰)⁻¹`, given by the intrinsic Laplace integral. Let `E f(S)=Σ_Y f(Y)`, where `Y` ranges over immediate core children retaining their parent slots. Write `T f` for the sum of `f` over all core descendant fringes, including the top; on finite shapes, `T=(I−E)⁻¹`. Let `Δf(S)=f(S⁺)−f(S)` add a leaf **at the top**, and `C f=p_S Σ_Y f(Y)` with `p_S=1/D_o`, `b=L_o/D_o`. The intrinsic Markov semigroup gives norm holomorphy of `V_z` for `Re z>0` by integrating the exponential power series. Here `||Δ||≤2` and `||C||≤1`. Define

```
H_z = (I−Δ V_z C)⁻¹ Δ V_z b,
U_z = V_z(b+C H_z).
```

The [unlimited occupation certificate](pair/README.md) gives `||Δ V_s C||≤.8600511260` for real `s≥l`. Its proof bounds the integral of the pointwise absolute paired occupation, so the identical majorant also applies for complex `Re z≥l` by `|exp(−zt)|=exp(−Re z·t)`. In any case the real estimate and norm holomorphy already give an open complex neighbourhood of each real `s≥l`, including `l`. Thus `U` is holomorphic for `Re z>l` and near the boundary points in question. The inverse on the gradient space has norm `<7.146`. This is the unique bounded solution of `(z−L⁰)U=b+CΔU`, with domain understood as the intrinsic integral formula (equivalently the pointwise generator equation for bounded scores by stopped Dynkin).

Let `M_t(S)=E⁰_S C_t` count all real cores from an **isolated planted seed**: the top slot is counted in degrees and never sends. The relation between a bounded fringe score and its whole-tree sum is

```
L⁰(TU) = T(L⁰U+CΔU+b U(J)).
```

Indeed, the only gain missing from each non-top fringe's intrinsic evolution is a leaf sent by its parent. Charge it to the parent's `p_SΔU(Y)`. Every promotion adds the virgin fringe and is charged to its promoter. The top has no incoming send; there is no centre without a slot in this identity. Consequently

```
(s−L⁰)TU_s = (1−U_s(J)) T b.
```

For real `s>μ`, integrate by finite-time stopped Dynkin. The boundary term `exp(−st)E⁰_S TU_s(S_t)` vanishes since `|TU_s|≤||U_s|| C_t` and `M_t(S)≤C_Q C_0(S)exp(μt)` by the spatial cell estimate (including its inert-slot extension). On finite time intervals the derivative and jumps admit polynomial size bounds and Yule localization. Since `T b` is the total promotion rate, `M'_t=E⁰_S T b(S_t)`. Put `M_s(S)=∫₀∞ exp(−st)M_t(S)dt`. Integrating by parts gives `TU_s(S)=(1−U_s(J))(s M_s(S)−C_0(S))`. As `(TU_s)(J)=U_s(J)` and `C_0(J)=1`, this yields

```
s M_s(J) (1−U_s(J)) = 1,                 s>μ.             (1)
```

The physical expectation is used first to justify (1) in its convergence domain. At no point do we treat the defining expectation integral as convergent at or below `μ`.

## 2. Positive trigger and the two-sided sign

Use the [permanent 19-label kernel](../../fringe/finite-kernel.md) and its [degree-two trigger](../../finite/gauge/physical-trigger.md). Write `G_z=(z−R)⁻¹`, `K_z=J_mark G_z B₀`, `j=J_(1/2)`, `h(S)=b(S) 1_(D_o=2)` (one half on a virgin and zero elsewhere). Here `J_mark` takes all 19 newborn evaluations. Let `B_nt` be the original direct production of all newborns **except** those triggered at degree two, and `Rhat=R+B_nt`. Write `Ghat_z=(z−Rhat)⁻¹` for its positive occupation. If `k_z=K_z e_2` and `C_nt(z)=K_z−k_z e_2^T` (column 2 set to zero), its resolvent is given by the convergent sum over nontriggered generations using `(I−C_nt)⁻¹`. This is not `C` from §1. Put `F(z)=(Ghat_z h)(j)`. The trigger result supplies a positive vector `g` at `s_0=.56564`, with `max g=1`, `min g>.117`, and `C_nt(s_0)g≤.979 g`; it therefore holds for all real `s≥s_0` and in norm for complex `Re z≥s_0` by positivity. Thus `Ghat_z` is holomorphic near the real domain of interest. `F(s)` is strictly decreasing and has critical equation `F(s)=1` iff `r(K_s)=1` there, by the positive trigger comparison. In particular `F(μ)=1`. No other pole from the untriggered denominator can occur in the target window.

Let `I_t^J` be the mean **permanent-mark immigration from the isolated top of J**, excluding the top itself from the count of descendants. Its row has direct newborn immigration at rate `b` and the usual child copies at total rate `a_Y` per old child, with the top degree counting its inert slot. Let `Ihat_z=∫₀∞ exp(−zt)I_t^Jdt`. The usual physical drift-transfer and weight proof work with the inert slot: the missing exterior incoming arrow is never included; the top still contributes one core and needs no kill or score. For real `s>μ` positive renewal and separation at the trigger give

```
M_s(J) = A(s) + B(s)/(1−F(s)),
A(z)   = 1/z + Ihat_z Ghat_z 1,
B(z)   = (Ihat_z Ghat_z h) (Ghat_z 1)(j).
```

Here `1` is active unit score on the permanent family. `A,B` are holomorphic on the stated domain; both are strictly positive on the real domain (in `B`, direct immigration of `j` and both unchanged original residence tolls already have positive mass). Equality (1) on `s>μ` gives, by the analytic identity theorem along the connected real domain with open complex neighbourhoods,

```
s [(1−F(s)) A(s)+B(s)] (1−U_s(J)) = 1−F(s).              (2)
```

This is an identity of analytic functions with **no division** by `1−F`. The bracket in (2) cannot vanish at a real point where `F≠1`; at the sole remaining point `μ` it is `B(μ)>0`. Hence it is positive everywhere on the real half-line under consideration. At `l`, the analytic identity extends to a neighbourhood of the boundary, so the same nonvanishing argument gives the strict sign when `l≠μ`. An interval for the **full** `U_s(J)` entirely above one certifies a lower exponent endpoint; entirely below one certifies an upper endpoint. No such inference is made from a killed signed value alone.

## 3. Quantitative inverse factor from the zero-mark positive row

Adjoin an auxiliary **positive-gauge** top label zero, used only for the starting root. Its strict descendants and all direct newborns still carry the original nineteen permanent labels. Let `G_s^0` denote evaluation of the permanent residual resolvent at this auxiliary virgin. Define `a_s=G_s^0 1`, `w_s=G_s^0 B₀` (a row), `v_s=J_mark G_s 1`. The physical full positive branching score from this virgin equals `M_t(J)`: sum its full backward score over the physical core fringes, with root mark zero and permanent marks below, and cancel parent copy and child kill terms. There is no uncancelled top kill, the initial virgin has no strict descendants, and the terminal sum of unit scores is the physical count. This is the same stopped transfer as in §1 with the positive rows; the weights below and Yule moments justify it. Hence for `s>μ`

```
M_s(J)=a_s+w_s(I−K_s)⁻¹ v_s.
```

For this row the original weight `W(J_0)=1/d=20/3` is admissible: the exact identity `R W≤ρW−h_weight(c_o)≤ρW−.7` uses only the descendant mark inequality and remains true at top mark zero. Stopped drift integration gives `a_s≤P:=200/21`; likewise `w_s g≤P` since `B₀g≤1`. On ordinary marked newborns `v_s≤M 1` with `M=50/7`, as `W(J_a)≤5`. These bounds hold throughout `[l,u]` and are the reason to use the auxiliary positive-gauge row for effective constants rather than estimating the physical top immigration term separately.

Keep the degree-two trigger `F` of §2 and let `D=(I−C_nt)⁻¹`. Rank-one inversion of `I−K=I−C_nt−k e_2^T` gives the same two terms `M_s(J)=A_0+B_0^sc/(1−F)` with

```
A_0 = a_s+w_s D v_s,
B_0^sc = (w_s D k)(e_2^T D v_s).
```

This is an alternative analytic decomposition of the same isolated response; `(1−F)A_0+B_0^sc` equals the bracket in (2). Use `m=.117`, `δ=.021`, and the fixed comparison vector `g` of §2. Then `D v_s≤[M/(mδ)]g`. Also `K_sg≤1.006 g` by the Perron property at `s_0`, the radius bound there, and monotonicity of entries. Thus `D k g_2=g+D(K_s−I)g≤v g`, with `v=1+.006/δ`. It follows that

```
A_0 ≤ A_*:=P[1+M/(mδ)],
B_0^sc ≤ B_*:=P v M/(mδ).
```

The immediate trigger contribution `1/[2(s+2)]` and positivity of all remaining chronological paths give `−F'(s)≥1/[2(u+2)^2]` on `[l,u]`. For `s≤μ` the positive bracket in (2) is at most `B_*`; for `s≥μ` it is at most `A_*+B_*` because `0<F≤1`. These give the uniform inverse factor

```
|s−μ| ≤ 2u(u+2)^2 (A_*+B_*) |1−U_s(J)|.                  (3)
```

The coefficient in (3) is less than `471522`. For a sharper critical derivative use `Ψ⁰(J)=lim_(t→∞) exp(−μt)M_t(J)` from the [spatial cell limit](../../fringe/spatial-limit.md) on the isolated seed. The virgin starts with one core, so this is its cell amplitude. It is strictly positive by the positive permanent renewal, or by the trigger simple pole with `B(μ)>0`. Bounded convergence after the rescaling `t=x/(s−μ)` gives `lim_(s↓μ) (s−μ)M_s(J)=Ψ⁰(J)`. Equivalently `Ψ⁰(J)=B(μ)/(−F'(μ))`. Differentiating (1) through its analytic inverse, or (2) at its common zero, yields the exact physical inverse slope

```
−U'_μ(J) = 1/(μ Ψ⁰(J)) = −F'(μ)/(μ B(μ)).                 (4)
```

In terms of the 19-label Perron pole, let `g_μ` be the right vector with maximum one and `ℓ_μ g_μ=1` for the left vector. Put `λ'(μ)=ℓ_μ K'_μ g_μ`. The auxiliary positive row gives

```
Ψ⁰(J) = (w_μ g_μ)(ℓ_μ v_μ)/(−λ'(μ)).
```

To compare without estimating `λ'` separately, the positive trigger formula gives `−F'(μ)≥1/[2(u+2)^2]`; a sharper permanent estimate follows directly from the retained *first original residence* in **every** newborn row. Indeed, the virgin has no copies before its first event. Write `K_{ij}(s)=N_{ij}(s)/(s+a_i+3/2)`, where `N_{ij}` is the immediate reward plus intrinsic successor resolvent contributions, and is a Laplace transform of a nonnegative measure. Hence `−K'_{ij}(μ)≥K_{ij}(μ)/(μ+2)`, giving `−λ'(μ)≥1/(μ+2)`. The minimum of `g_μ` is `>.117` by the trigger proof's critical bound. Therefore `ℓ_μ 1≤1/m` and `Ψ⁰(J)≤(μ+2) PM/m`. Inserting this into (4) yields `−U'_μ(J)>1/845` by rational arithmetic. This a priori derivative bound is strictly at the root; the full-window inverse factor uses (2), not an assumption of derivative monotonicity.

Formula (4) also agrees with the isolated eigenfunction transform: `(I−E)Ψ⁰/Ψ⁰(J)` solves the critical bounded signed equation, and bounded uniqueness by contraction identifies it with `U_μ`. In particular the top gradient is exactly `ΔΨ⁰/Ψ⁰(J)`; an attached physical seed with an active parent is not interchangeable with this isolated expectation.

## 4. Evaluator controls

For a bounded approximate gradient `H̃` let `ε` bound the **full, unlimited** defect `||H̃−ΔV_s(b+C H̃)||∞`. Then the error in the scalar proposal `V_s(b+C H̃)(J)` is at most `ε (V_s C1)(J)/(1−.8600511260)`. Add any certified error from evaluating the latter intrinsic integrals. A killed defect alone cannot substitute for `ε`. Alternatively, for the Neumann iteration starting from zero, truncation after `m` gradient terms costs at most `.8600511260^m/[s²(1−.8600511260)]` in the newborn scalar: `||ΔV_s b||∞≤1/s` since `0≤V_s b≤1/s`, and `V_s C1≤1/s`.

If the final scalar enclosure is `[L_s,H_s]`, then `L_s>1` proves `s<μ`, and `H_s<1` proves `s>μ`. For `s∈[l,u]`, the enclosure also implies `|s−μ|≤471522 max(|1−L_s|,|1−H_s|)`. This is a residual-to-exponent control, not a successful new numerical bracket. The small finite signed replay below tests the generator signs and ordinary finite equation separately from the unlimited sign theorem.

```
taskset -c 8 timeout 5s python3 -B round2/influence/zero/critical_controls.py
taskset -c 8 timeout 8s python3 -B round2/influence/zero/independent.py 3
```
