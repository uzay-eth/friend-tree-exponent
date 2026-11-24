CORRECT

# Cold audit: degree-two physical trigger scalar

Scope. Audit of `round2/finite/gauge/physical-trigger.md` with `physical_trigger_check.py` and `physical_trigger_bounds.py`, plus an independent small-tree checker. The question is whether the nineteen-mark column split defines a scalar independent of μ, whether the explicit subcritical endpoint and its arithmetic are correct, and what inputs are inherited. As in the other verdicts, nothing here is a claim about the arithmetic of the critical root. No author file was changed and nothing was committed. Author: hy4.

## 1. What was rerun and independently recomputed

One core, from the repository root:

```
taskset -c 0 python3 -B round2/finite/gauge/physical_trigger_check.py
taskset -c 0 python3 -B round2/finite/gauge/physical_trigger_bounds.py
taskset -c 0 python3 -B round2/checks/physical_trigger_hy4.py
```

The shipped checks reproduce: 6/24/92 marked finite states at caps 3/4/5 with the same number of exact state scalar/resolvent identities and a checked virgin Schur elimination; six zero incomparable mixed differences and one nonzero comparable-pair difference; the reachability counts; the label-20 lower 0.022214534253 and row sum 0.405460197323 at 0.5665; g₂=0.119101711051 with complement contraction 0.971863522077 at 0.5665 and g₂=0.117936293836 with contraction 0.978127513699 at 0.56564; the wider-bracket contraction below 0.9981; and the printed fine-bracket factors.

`physical_trigger_hy4.py` imports no author module. It reimplements the tree rates, the finite killed system and a rational solver, and checks at caps 3 and 4 that the killed trigger scalar equals `ev_{J_{1/2}}(s−\widehat R)^{-1}h` computed by the column split, that K v−v=(F−1)k exactly, that F decreases on a grid, and that sign(F−1) matches sign(r(K)−1) by power iteration. Deliberate errors — triggering on every degree-two state, dropping the non-trigger b reproduction, and using the wrong renewal column — are all rejected. Output ends with `PHYSICAL TRIGGER COLD CHECK: PASS`.

## 2. Well-definedness independent of μ

The split is well-defined. With k=K_s e₂ and C=K_s−k e₂ᵀ (zero column 2), the definition v_s=(I−C)^{-1}k and F(s)=v_s(2) uses only K_s at the chosen s; μ nowhere enters it. The only role of μ is the equation F(μ)=1, which follows from r(K_μ)=1 in the reviewed renewal. The identity K v−v=(F−1)k is exact, and because k>0 in every row, sign(F−1)=sign(r(K)−1) including the zero case. My small-tree recomputation confirms both identities and the sign correspondence, so the scalar is genuinely independent of μ as a function. The needed condition r(C_s)<1 is certified below, not assumed.

## 3. Explicit subcritical endpoint

At s₀=0.56564 the six-label iteration with R_*=1.006 gives g₂=0.117936293836 and C_{s₀}g≤0.978127513699 g<0.978128 g. At u=0.5665 with R_*=1 it gives g₂=0.119101711051 and 0.971863522077<0.971864. Both are exact integer-power/rounding computations, reproduced independently here. The wider-bracket statement at s_*=0.5649 gives C_{s_*}(g^{1−p})<0.9981 g^{1−p}. Since 0.5649<0.565≤μ, these are genuine subcritical endpoints below the physical interval, and they cover the interval through monotonicity of the entries.

## 4. The printed .5665 terminal and the c≥2 correction

The printed drift is correct. For a leaf gain, Δb=c/[D(D+1)]; for a promotion, rate b=L/D and Δb=−1/D, giving L⁰b ≥ Lc/[D(D+1)]−L/D² = L[cD−(D+1)]/[D²(D+1)]. For c≥2 this is nonnegative for every D≥2, so b is nondecreasing and the high-degree terminal uses plain 1/(s+1/2), i.e. zero correction. For c=1 it equals −b/[D(D+1)], giving η=1/420 at D=20. These equalities were checked exactly for all 2≤D≤20 and 1≤c≤D. At u=0.5665 the sharper terminal yields the printed 0.022214534253, strictly larger than the uniform 0.022179480005; both numbers and the direction of the comparison are as printed.

The separate child message was supplied in full and is verified exactly. `physical_trigger_hy4.py` reimplements the first-event stencil from scratch and reproduces all three of its numbers: the uniformly charged 1/420 tail gives 0.022244698369 at 0.56564 and 0.022179480005 at 0.5665, while setting the correction to zero for c≥2 (where L⁰b≥0) gives the printed 0.022214534253 at 0.5665. The sharper value is indeed larger at both points, and both versions are valid: L⁰b≥−b/[D(D+1)] holds for every c≥1, so the uniform charge is a valid but weaker bound, while the c≥2 zero correction is valid and stronger. No contradiction with the printed text.

## 5. Inherited input scope

Not reproved here, and inherited by this verdict: the reviewed physical bracket 0.5656581227≤μ≤0.5657979941<0.5658; the original weight bound (ρ−R)W≥h₀(c)≥7/10 with W(J_{1/i})≤5 and the Hölder anchor at 0.51; the degree certificate b_*=0.462, d_*=0.19, h_min=0.765931373 with its 189 inequalities; the nineteen-mark finite generator and finite-size drift; and the threshold transfer r(K_μ)=1 from the reviewed renewal. The six-label iteration's R_*=1.006 at s₀ is derived from those inputs, not independently established.

## 6. Scope boundary and verdict

This verdict covers §§1–4 of `physical-trigger.md` and the two named scripts `physical_trigger_check.py` and `physical_trigger_bounds.py`, plus the child message supplied above. The author has since appended separate material in `physical_trigger_verify.py` and `outputs/finite-physical-trigger-verify.log`: a broad interpolation reaching a .5641 contraction and an H=6, K=3 chain-compression missing return. That material is outside the sections audited here and is not reviewed or endorsed by this verdict. No conclusion is drawn about μ or about arithmetic of μ. I am a different model family (hy4) from the Codex author and no further child was spawned.

The column split gives a well-defined scalar independent of μ; the explicit subcritical endpoints and all printed arithmetic checked here, including the child message, are correct.
