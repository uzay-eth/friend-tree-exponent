CORRECT

# Cold audit: all-depth star/copy continuation, 0.807760 and 4800·e^{0.5630t}

Scope. Audit of `recursive-excursions.md`, `first-promotion.md` and the independent replay of `excursion-evidence.log`, starting from the reviewed binary renewal and positive dormant split in `scalar-hy4.md`. The verdict is replay-confirmed and algebraically re-derived, with one code-level limitation stated in §3. No author file was changed and nothing was committed. The response-transcendence material is deliberately excluded; it has its own verdict. Author: hy4.

## 1. What was actually rerun

One core, sequentially, from the repository root:

```
bash round2/finite/gauge/excursion-replay.sh > round2/checks/excursions_hy4_replay.log 2>&1
taskset -c 0 python3 -B round2/checks/excursions_hy4.py
```

The replay executes `excursion_check.py`, all six `excursion_rows.py` cases, and `excursion_constants.py`. Every numerical line of the fresh log matches `excursion-evidence.log` exactly: 58 raw / 24 Schur states and 116 exact row-column equalities; 954 retained multicore states with contraction ≤ 0.2954199 at u=0.5665, ≤ 0.2970853 at the four-bin anchor 0.54001, ≤ 0.2955441 / 0.2955169 at the two Sabba points, 156 states with ≤ 0.2077579 at s=1; fixed-point weighted errors ≤ 7·10⁻¹⁴. Wall times were 0.71–4.18 s per case (CPU user 0.22–1.37 s), all inside the five-second wrapper; the recorded log's 0.57–0.84 s wall times reflect a less loaded host. The N=9 sweep therefore passes on one core.

## 2. Independently re-derived algebra

`round2/checks/excursions_hy4.py` imports no author module. It reimplements the star weights, the γ telescoping identity, both tail moments, the final rational bounds and negative controls with exact fractions:

- ω_{j,j}=b_j/(λ+j−1+b_j) and ω_{j,k+1}=ω_{j,k}k²/(k²+(2+λ)k+λ); for j=2 this reproduces q₂=1/(3+2λ) of `first-promotion.md`. Weights decrease in k.
- The identity (k+γ)q_k−(k+1+γ)q_{k+1}=(1+λ)q_k−λ²q_{k+1}/[(λ+2)k²] with γ=λ+2/(λ+2) holds exactly, and a_k=k(1−r_k) increases to 2+λ. A 2·10⁶-term direct float sum lies inside the exact moment intervals for λ=s+1/2, s+1/16 at u and for λ=1+1/2.
- Four fresh inputs, rounded outward exactly as in `excursion_constants.py`: k_A≥0.5225363, k_β≥0.0706297, d_A≥0.2995629 and positive dormant row <4.002 at 0.54001. Then 1−(1+d_A)k_β/(1−k_A)=38567596224187/47746370000000<0.807760 exactly.
- Hölder at s=0.5630 has exponent p=350/2649; the integer-power comparison 4.002³⁵⁰·0.807760²²⁹⁹<0.9980²⁶⁴⁹ holds, and 1+w/[(1−0.9980)h_min]=1+[0.974/(1/16+0.14)]/[(0.002)·0.501218471] has ceiling exactly 4800. The same script reproduces the 160·e^{0.564t} and Sabba-endpoint multipliers.

Small-tree kernels were checked two ways: the shipped finite reference (58/24 states, 116 exact equalities) was rerun, and the independent script confirms the star elimination is the positive-series sum over first promotions with no residual child copied during the star phase. Negative controls behave as required: the wrong ratio k²/(k²+(1+λ)k+λ) breaks the γ identity; swapped moment endpoints fail to enclose the direct sum; dropping (1+d_A) gives 1−k_β/(1−k_A)>0.807760; cap 0.9981 pushes the multiplier above 4800; anchor 4.1 breaks the 0.9980 interpolation.

Output of the independent checker:

```
fresh rows match recorded log exactly: 24 parsed quantities
all six cases: nonnegative ordered rows, contraction < 1, 954/156 state counts OK
star weights, gamma telescoping identity, moment enclosures, monotone a_k: PASS
negative controls (wrong ratio / swapped moment endpoints) rejected
four fresh inputs -> critical excess 38567596224187/47746370000000 < .807760
Holder at .5630 gives d_beta < 499/500 ; multiplier ceil = 4800
negative controls (drop (1+d_A), cap .9981, anchor 4.1) all fail as required
EXCURSIONS COLD CHECK: PASS
```

## 3. Limitation and unresolved source prerequisites

- The 954-state enumeration was **re-executed, not re-implemented**. State generation, the finite solve and the contraction sweep remain the authors' code; my independent checks cover the surrounding star/moment algebra and all final inequalities, plus the shipped finite reference. This is a replay confirmation, not a from-scratch second solver of the 954-state system.
- `recursive-excursions.md` §5 states the size-drift bound for δ=0.013 (basic) or δ=0.025 (four-bin) only for s≥0.565. The 0.54001 anchor uses the four-bin profile. The finite-N enclosure at that anchor is still valid because its boundary tables are certified at fixed degree cutoff, but the stated N→∞ convergence of that particular anchor is not established by §5. The interpolation uses the N=9 value; this is a scope gap, not an arithmetic failure.
- The following remain external prerequisites, already reviewed elsewhere and not reproved here: the binary row/Lyapunov tables (q=49/125, d=4/25, h_min=0.501218471, α, κ), the four-bin virtual tables (Q=2/5, DAMP=7/50), the positive-dormant potential identity W≤W+m v<[1+w/((1−d_β)h_min)]W from `dormant-return.md`, and the threshold transfer r(K_s)=1 at μ that converts the critical d_β formula into a statement at μ. The 0.5630 bound inherits all of these.
- The 0.5630 bound is a uniform upper bound on the **enlarged residual with identical direct newborns**, not a new upper bound on the physical exponent. Nothing here identifies μ or changes the open endpoint signs.

## 4. Verdict

The all-depth star/copy continuation passes its replay on one core; the four inputs are reproduced and the 0.807760 and 4800·e^{0.5630t} derivations are arithmetically correct given the reviewed dormant potential and binary certificates. Conditional on those stated prerequisites and on the code-level limitation above, both target numbers are supported.
