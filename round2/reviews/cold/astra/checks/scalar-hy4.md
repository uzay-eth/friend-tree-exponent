CORRECT

# Cold audit: two-label determinant, scalar first return and dormant gap

Scope. First-principles audit of `binary-and-return.md` (`finite/gauge/binary-and-return.md`),
`dormant-return.md` (`finite/gauge/dormant-return.md`) and their binary/virtual Lyapunov and scalar-return inputs
(`binary-row-bounds.md` (`finite/gauge/binary-row-bounds.md`), `binary_weight.py` (`finite/gauge/binary_weight.py`),
`virtual_weight.py` (`finite/gauge/virtual_weight.py`), `return_constants.py` (`finite/gauge/return_constants.py`),
`transfer_check.py` (`finite/gauge/transfer_check.py`), `dormant_check.py` (`finite/gauge/dormant_check.py`)).
Author: Codex. No author file was changed; no long enumeration was rerun. The later high-hub work and the response-
transcendence argument are outside this review and were not used.

## 1. Binary degree rows (independent recomputation)

I reimplemented the reverse degree recursion from `binary-and-return.md` §1 from scratch (marks A=1/2 below degree 16,
β=1/16 otherwise; d=4/25; q=49/125; cutoff 64; upward rounding at 10⁻⁹). My results reproduce the frozen log exactly:
h(2,1)=0.998871532<1, hmin=0.471910343, hmax=1.079138952; every one of the 2·D degree/core endpoint drift
inequalities holds and the size condition δ·q·hmax/(β+d)≤hmin holds with δ=13/1000. The four-bin virtual tables
independently reproduce required/allowed newborn credits and hmin=0.501218471, hmax=1.072734084, δ=1/40, with all
6,232 endpoint rows satisfied. So the residual is a genuine nonnegative supersolution with L⁰W≤0.4W and RW≤0.54W−h.

## 2. Determinant and scalar algebra

For K=[[k_A,k_β],[k_A+d_A,k_β+d_β]] I verified exactly with rationals:

  det(I₂−K) = (1−d_β)(1−F),  F = k_A + k_β(1+d_A)/(1−d_β),

and that r(K)=1 iff det(I₂−K)=0 at r=1. Since D_s=K−1·k_s has zero first row and second row (d_A,d_β), r(D_s)=d_β and
(I−D_s)⁻¹1=(1,(1+d_A)/(1−d_β)), giving precisely the stated series F=k_s(I−D_s)⁻¹1=k_sΣD_sⁿ1. The frozen lower bounds
k_A≥0.471314072924, k_β≥0.049408249826, d_A≥0.193161058005 at u=0.5665 are strictly positive (positive residual),
and my own Perron/power comparison reproduces d_β(μ)<0.8884949. With d_β<1 the scalar denominator is strictly
positive throughout the prior bracket, so F(s)=1 and det(I₂−K_s)=0 have the same unique root locally; uniqueness uses
the positive renewal transfer r(K_s)=1 at μ and strict decrease of the immediate common newborn return, both from the
referenced renewal/finite-kernel prerequisite and not reproved here.

## 3. Dormant induction on small trees and positive gap

I implemented the planted-fringe generator myself (leaf gains, core promotions with the degree-16 mark rule, child
residual arrows) and checked it coincides with the shipped process on the small tree set. On a 54 active / 27 dormant
finite killed system I verified: positive excess d_A>0, 0<d_β<1; the geometric dormant return d_β/(1−d_β); the scalar
return of the enlarged process equals k_A+k_β(1+d_A)/(1−d_β); and the full-occupation identities before and after the
dormant split agree. Using my own log-convexity integer-power comparisons I obtain d_β(0.56470)<0.99917 and
d_β(0.56480)<0.99268, hence comparison factors below 11,600 and 1,320 and an enlarged-residual growth bound with
exponent 0.56470 (respectively 0.56480), strictly below the established physical lower bound 0.565. This is a gap
for the *enlarged residual with identical direct newborns*, not a new upper bound on the physical exponent.

## 4. Which gap is independently supported

**0.56470 is the independently supported gap.** It follows from the short rational binary-row/virtual-Lyapunov
certificates and the dormant potential identity, all of which I recomputed here on one core in about one second.
The later **0.5630** comes from `recursive-excursions.md` (`finite/gauge/recursive-excursions.md`): a 954-state all-depth
numerical scheme giving copy contraction <0.296, d_β(μ)<0.807760 and d_β(0.5630)<0.9980, with multiplier 4,800. Its
four numerical inputs are read from the recorded rows in `excursion-evidence.log` (`finite/gauge/excursion-evidence.log`)
by `excursion_constants.py` (`finite/gauge/excursion_constants.py`); establishing them requires the long 954-state
enumeration, which I did not rerun. I ran `excursion_constants.py` against the recorded log: its four inputs match and
the 0.5630 → multiplier 4,800 algebra is correct, but the enclosures themselves were not reproduced. Within this scope the
0.5630 gap is therefore **conditional**, not independently supported; the 0.56470 gap is unconditional given the
reviewed binary certificates.

## 5. What is required to identify the physical exponent

The gap bound alone does **not** identify μ: the dormant construction changes the process, and the bound 0.56470<0.565
is a subcriticality certificate for that enlarged residual. The exponent is identified by (3) of `dormant-return.md`
evaluated at μ together with the original threshold transfer r(K_s)=1 ⟺ s=μ. For a two-sided numerical identification
one still needs a rigorous enclosure of the four true entries of the original two-label kernel near μ (equivalently
k_A,k_β,d_A,d_β), from a finite killed resolvent with explicit truncation error, giving **opposite signs** of
F(s)−1 (equivalently det(I₂−K_s)) at bracketing trial points. The present material supplies only lower bounds at u
and wide upper intervals from the separate nineteen-label/excursion scheme, and it explicitly states that the
Sabba endpoint signs remain open (`README.md` (`finite/gauge/README.md`)). So the determinant characterisation and
scalar uniqueness are established, but a two-sided value of μ is not produced by this work.

## 6. Controls and limitations

Commands, from the repository root on one core:

```
taskset -c 0 python3 -B round2/checks/scalar_hy4.py
```

Output: `SCALAR RETURN COLD AUDIT: PASS`. Planted errors rejected: shifting the mark threshold 16→17 makes
h(2,1)>1; choosing q+d≥0.565 destroys the subcritical gap; omitting the (1+d_A) factor or using 1+d_β in the scalar
denominator breaks det(I₂−K)=(1−d_β)(1−F). I did not rerun the 954-state excursion enumeration, the high-hub work or
the transcendence argument; those are separately scoped. No defect was found in the audited finite two-label
determinant, scalar first return, positive residual or strict gap.
