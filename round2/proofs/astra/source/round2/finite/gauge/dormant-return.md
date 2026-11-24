# Identical direct newborns and a positive residual with a strict gap

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The [all-depth numerical scheme](recursive-excursions.md) eliminates every one-core visit with four precomputed tail moments. Its 954-state calculation gives a copy contraction below $`0.296`$, a dormant return below $`0.807760`$ at criticality and a uniform enlarged-residual upper exponent $`0.5630`$. The calculations in section 3 below use the simpler two-coordinate bounds for an independent initial gap.

## 1. A rank-one birth operator

Put $`A=1/2`$, $`\beta=1/16`$. Start with the [binary fringe construction](binary-and-return.md). Its active particles have the original residual generator $`R`$. Add one kind of dormant particle. Its current state is a whole old fringe whose root has mark $`\beta`$. It runs the original intrinsic dynamics with a silent exterior parent, retaining its marks, leaves and old children at every intrinsic event. It dies at rate $`\beta`$ and switches to the active phase at rate $`A-\beta`$, preserving the whole then-current shape. It emits nothing while dormant.

In the new residual process, an active particle also creates a dormant $`J_\beta`$ at rate $`b_H=b\,\mathbf1_{\{D\ge16\}}`$. The full process adds the *same* direct newborn $`J_A`$ at rate $`b`$ from every active particle. A simultaneous direct and dormant birth at high degree, or separate clocks with these intensities, give the same first-moment generator. Thus the direct operator is genuinely rank one. The internal dormant state and the active copies with birth mark $`\beta`$ are part of this construction. The literal process using permanent mark $`1/2`$ on every old copy has different rates.

Let $`\mathcal R_0`$ be the generator with dormant activation but without new dormant births. For a pair of scores $`f=(f_{\rm act},f_{\rm dorm})`$ its actions are
```math
 (\mathcal R_0 f)_{\rm act}=R f_{\rm act},\qquad
 (\mathcal R_0 f)_{\rm dorm}=(L^0-A)f_{\rm dorm}+(A-\beta)f_{\rm act}(S_\beta).
```
Write $`J_D`$ for the initial dormant state and $`r_H=(b_H,0)`$. The new residual and direct operators are
```math
\widetilde R=\mathcal R_0+r_H\,{\rm ev}_{J_D},\qquad
 \widetilde B=(b,0)\,{\rm ev}_{J_A}. \qquad\text{(1)} 
```
Each positive operator in (1) has its stated direct branching interpretation. There is no cancellation between sampled particles.

## 2. Original-root survival pays exactly for the dormant particle

Take any nonnegative physical newborn test $`g=(g_A,g_\beta)`$. Under $`\mathcal R_0`$, an active particle has response $`U_sg`$ to the running reward $`b g_{a(D)}`$. A dormant starts with no reward and then pays
```math
 (A-\beta)(V_{s+A}^0 U_\beta)(J)
       =(K_sg)(\beta)-(K_sg)(A). \qquad\text{(2)}
```
Indeed $`U_x=V^0_{s+x}G`$ on a fixed intrinsic shape, where the original root's remaining copy forcing and direct reward $`G`$ are independent of $`x`$. The intrinsic resolvent identity proves (2) with the positive sign. In particular the expected discounted number of new dormant particles from one $`\mathcal R_0`$ family starting at $`J_D`$ is precisely $`d_\beta(s)`$.

Under $`\widetilde R`$ let $`h`$ be the discounted number of future direct births from $`J_D`$. Each physical newborn opportunity now pays one direct root, and high opportunities also start an independent dormant residual. Equation (2) gives $`h=d_A+d_\beta(1+h)`$. By positive summation over dormant generations, if $`d_\beta<1`$, the actual common newborn kernel in (1) is
```math
 \widetilde K_s={\rm ev}_{J_A}(s-\widetilde R)^{-1}(b,0)
 = k_A+k_\beta\frac{1+d_A}{1-d_\beta}. \qquad\text{(3)}
```
It equals one exactly when $`r(K_s)=1`$ in this domain, by comparison in [the scalar first return](binary-and-return.md). The [algebraicity obstruction](isolated-transfer.md) applies solely to the scalar function $`\widetilde K_s`$. It has no implication for the algebraicity of the critical root $`\mu`$. Positivity of the immediate direct contribution makes $`\widetilde K_s`$ strictly decreasing there.

The full active occupation of a mark-independent active observable also agrees with the old full process. From any fixed shape, split an original root of mark $`\beta`$ into an active copy with root killing rate $`A`$ and a dormant copy of that shape. In first moment an old potential root reward at intrinsic age $`t`$ is integrated against $`e^{-\beta t}=e^{-At}+\int_0^t(A-\beta)e^{-Av}e^{-\beta(t-v)}\,dv`$. The Markov property at $`v`$ justifies keeping the complete intrinsic shape at activation. This equality includes rewards produced by residual and direct offspring, whose forcing is independent of the original root mark. Apply the split at each high-degree direct birth and sum over generations by monotone convergence. The same replacement handles high-degree direct immigration from the initial centre. For old core-child immigrants the real originating binary types can be retained.

## 3. A uniform bound on the enlarged residual

On active states use the weight $`W`$ from the four-bin binary certificate. On dormant states use the same $`W`$ for that shape with root mark $`\beta`$. Then $`\mathcal R_0 W\le0.54W`$. The dormant component is at most $`(q-\beta)W`$ with $`q=0.4`$. The active component has the sharper bound $`RW\le0.54W-h_{\min}`$, $`h_{\min}=0.501218471`$. Thus, for $`s>0.54`$,
```math
v=(s-\mathcal R_0)^{-1}r_H\le W/h_{\min},\qquad v(J_D)=d_\beta(s).
```
The first inequality uses $`(s-\mathcal R_0)W\ge h_{\min}r_H`$. All resolvents here have nonnegative kernels and polynomial bounds on each finite time interval.

At the critical point the two Perron equations and the [verified binary rows](binary-row-bounds.md) give $`d_\beta(\mu)<0.888495`$. The positive intrinsic survival difference (2) gives $`d_\beta(0.54001)<5`$ from the controlled table. Hölder makes $`d_\beta`$ log convex. For $`s_*=0.56470<0.565\le\mu\le u=0.5665`$, the worst interpolation point is $`u`$. Consequently
```math
d_\beta(s_*)\le5^{180/2649}(0.888495)^{2469/2649}<0.99917.
```
The last inequality is a direct integer-power comparison. Taking $`s=0.56480`$ instead gives $`d_\beta(s)<0.99268`$.

Put $`w=W(J_D)=0.974/(\beta+0.14)`$, $`m=w/(1-d_\beta(s_*))`$ and $`\widehat W=W+m v`$. By (1) and the base resolvent equation,
```math
(s_*-\widetilde R)\widehat W
 =(s_*-\mathcal R_0)W+[m-w-m d_\beta(s_*)]r_H\ge0. 
```
Moreover $`W\le\widehat W<11600\,W`$. The original residual and dormant processes are nonexplosive by linear total-size bounds, and the extra dormant birth rate is at most the number of active particles. Alternatively the positive drift inequality gives the stopped first-moment comparison directly. Therefore every seeded active or dormant state satisfies
```math
\widetilde R_t W\le11600\,e^{0.56470t}W,\qquad 0.56470<0.565\le\mu. \qquad\text{(4)}
```
At $`s=0.56480`$ the same construction gives the smaller multiplier $`1320`$. This is a uniform gap for the enlarged residual with identical direct newborns. Independently, (3) at $`\mu`$ and the original threshold argument identify the physical exponent.

There is also an explicit finite first-exit comparison. Let $`\delta=0.02470`$ and $`n`$ be fringe size. The binary degree certificate gives $`R(n^\delta W)\le s_*n^\delta W`$, because $`\delta q h_{\max}/(\beta+d)\le h_{\min}`$. On dormant states the base bound follows from the intrinsic size increment and $`n\ge2`$, giving $`[q-\beta+\delta(1+q/2)]n^\delta W<s_*n^\delta W`$. In the preceding construction replace $`W`$ by $`n^\delta W`$ in $`\widehat W`$ and $`w`$ by $`2^\delta w`$, retaining the same $`v`$. This gives a positive supersolution dominating $`n^\delta W`$ and at most $`11600 n^\delta W`$. Active copies cannot exceed the source size, and activation preserves size. Killing at the first size $`N+1`$ therefore gives the factor $`(2/(N+1))^\delta`$ on the initial newborn's discounted exit-weight account, with multiplier $`11600 W(J_A)`$. Converting this exit weight to missing reproduction uses the finite bound on $`\widetilde K`$'s residual reward, for instance (4) and $`s>s_*`$. Thus (3) can equally be computed from one scalar finite killed resolvent with an explicit error.

At the initial $`J_A`$ there are no old children or high birth opportunities. Its intrinsic rate is $`3/2`$ and its killing rate is $`1/2`$. The first-event equation makes $`(s+2)\widetilde K_s`$ nonincreasing, since the running reward is $`1/2`$ and all subsequent discounted contributions are nonincreasing. Hence $`|\mu-s|\le(s+2)|\widetilde K_s-1|`$ wherever (3) converges and the intermediate real segment is in this domain.

## 4. More dormant labels and a check on the generator

For any original finite mark set, create a dormant copy of label $`a<A`$ when the actual birth has that label, kill it at $`a`$, and activate at $`A-a`$. Equation (2) gives the nonnegative dormant-to-dormant matrix obtained from the corresponding rows and columns of $`D=K-\mathbf1 k`$. Whenever its radius is below one the positive potential estimate also extends. Let $`r_j`$ be the physical reward giving births of dormant label $`j`$, $`v_j=(s-\mathcal R_0)^{-1}r_j`$, and $`w_i=W(J_{D,i})`$. Then $`v_j(J_{D,i})=D_{ij}`$. With $`m=(I-D_{\rm dormant})^{-1}w`$, the weight $`W+\sum_j m_j v_j`$ is a supersolution by the calculation in section 3. For the nineteen-label scheme the previously bounded $`r(D_{0.556})<1`$ therefore gives a rank-one direct generator with a residual bound at $`0.556`$, conditional on the nineteen-label input under independent review.

As a small algebraic check for the actual threshold 16, keep 54 active states and 27 dormant states, including both root marks, all star degrees through 16 and selected multi-core old shapes. Kill every excluded intrinsic or copy lineage. At $`s=2`$ the checker solves the finite rational systems both before and after the split. All 54 full active-occupation identities, the dormant return and the scalar common return agree exactly. [Finite generator replay](dormant_check.py) · [Integer bounds](return_constants.py) · Recorded output
