# Summing every one-core visit in the scalar regenerative kernel

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. Unknowns and the exact star row

This is an evaluation scheme for $`\widetilde K_s=k_A+k_\beta(1+d_A)/(1-d_\beta)`$ from [the rank-one construction](dormant-return.md). It computes $`k_j=K_s(A,j)`$ and $`d_j=K_s(\beta,j)-K_s(A,j)`$ using positive residence-time integrals. At every depth it removes all one-core states, including old-star and plus-star copies. The only remaining finite-system unknowns are active **multicore** trees of size at most $`N`$. A single additional backward sweep gives each dormant output.

Let $`S_j^a`$ be a star with root mark $`a`$, total degree $`j\ge2`$, one parent slot and $`j-1`$ leaves. Write $`\lambda=s+a`$, $`b_j=(j-1)/j`$. The first intrinsic promotion at degree $`k`$ gives the state $`V_k^a`$ with the same root, $`k-2`$ leaves and one virgin child of mark $`a(k)`$, where $`a(k)=A=1/2`$ below 16 and $`\beta=1/16`$ thereafter. Its discounted mass from $`S_j^a`$ obeys
```math
 \omega^\lambda_{j,j}=\frac{b_j}{\lambda+j-1+b_j},\qquad
 \omega^\lambda_{j,k+1}=\omega^\lambda_{j,k}\frac{k^2}{k^2+(2+\lambda)k+\lambda}.
```
No residual child can be copied during this star phase. The exact response satisfies
```math
U_a(S_j)=\sum_{k\ge j}\omega^\lambda_{j,k}\,[g(a(k))+U_a(V_k)]. \qquad\text{(1)}
```
The direct reward and first promotion have the same predictable intensity. Thus (1) applies afresh whenever a lineage first enters **any** star, regardless of how often the lineage previously made a copy or started a dormant excursion.

Keep the ordinary active first-event rows for all reachable multicore states of size at most $`N`$. A successor which is a retained multicore stays an unknown. A successor which is a star is substituted by (1). A multicore successor above $`N`$ receives a certified boundary interval. The intrinsic part of a multicore row cannot enter a star. Its root retains its original core children. Every retained intrinsic jump increases size by one. Copies may visit old multicore states or stars of arbitrary degrees. For a finite starting set and $`N`$, each incoming star degree $`j`$ in these rows is at most $`N`$. Positive series over all descendants, with absorption at the first multicore oversize on each lineage, justify substitution simultaneously at every copy depth.

## 2. Two scalar tail moments per residence rate

Use the controlled boundary in [the binary row bounds](binary-row-bounds.md), with degree boundary $`K=\max(32,N+1)`$. Put $`u^\pm=\beta H^\pm_{a,\beta}(2,1)`$ and $`t^\pm=g(\beta)+u^\pm`$. At a multi-core exit before degree $`K`$, apply the two controlled root tables, adding each old child's $`\ell H^\pm_{a,\ell}(D_Y,c_Y)`$ and the upper old-grandchild weight from that proof. The code intersects this upper with $`\alpha W`$. In particular it keeps the whole old child in $`V_k`$.

For every $`k\ge K`$, the root of $`V_k`$ has $`c=2`$, $`b=1-2/k`$, $`\eta=0`$, and the sole old child is $`J_\beta`$. The latter contributes $`u^\pm`$ with zero old-grandchild weight. Direct root reward plus the lower/upper postpromotion continuation in (1) is therefore bounded by
```math
 i_- - h_-/k\ \le\ g(\beta)+U_a(V_k)\ \le\ i_+ - h_+/k,
 \quad i_\pm=g(\beta)+u^\pm+t^\pm/\lambda,\quad
 h_-=2t^-/\lambda,\quad h_+=2t^+/(\lambda+1). \qquad\text{(2)}
```
The upper terminal uses $`\int e^{-\lambda t}\mathbb E b(t)\,dt\le1/\lambda-(1-b)/(\lambda+1)`$. The lower uses $`b/\lambda`$. The coefficients are computable from two-dimensional tables, with no unknown newborn response or derivative in them.

Only the moments $`T_{j,K}^\lambda=\sum_{k\ge K}\omega^\lambda_{j,k}`$ and $`M_{j,K}^\lambda=\sum_{k\ge K}\omega^\lambda_{j,k}/k`$ occur in (2). All starts $`j\le K`$ use the same pair of normalised moments, multiplied by $`\omega^\lambda_{j,K}`$. The two active root marks give four scalar moments in total for a specified $`s`$. They are independent of the test column and every old tree.

For exact enclosures, sum a rational prefix through $`P-1`$, with $`P\ge K`$. The rational two-sided telescoper for the remaining $`T`$ is (4) in [the first-promotion calculation](first-promotion.md). Set $`r_k=\omega^\lambda_{j,k+1}/\omega^\lambda_{j,k}`$ and $`a_k=k(1-r_k)`$. Direct differentiation shows that $`a_k`$ increases to $`2+\lambda`$. With $`b_0=2+\lambda`$ the numerator of its derivative is $`(b_0^2-\lambda)k^2+2b_0\lambda k+\lambda^2>0`$. Since $`\omega_k/k=(\omega_k-\omega_{k+1})/a_k`$,
```math
 \omega_P/(2+\lambda)\ \le\ M_{j,P}^\lambda\ \le\ \omega_P/a_P. \qquad\text{(3)}
```
Every interval is rounded outwards. In the lower terminal combine the lower $`T`$ with the upper $`M`$, and reverse them in the upper terminal. At fixed $`K`$ both normalised moment intervals converge with increasing $`P`$.

## 3. A contraction for successive copies

Write the finite positive multicore equations after substitution as
```math
u=h+Iu+Cu.
```
Here $`I`$ contains only ordinary retained intrinsic jumps. All copy terms, including copy visits integrated through (1), are in $`C`$. As $`I`$ strictly increases size it is nilpotent. Eliminate it by a descending-size sweep and iterate only
```math
u=(1-I)^{-1}h + Tu,\qquad T=(1-I)^{-1}C. \qquad\text{(4)}
```
For positive trial weights $`w`$, the exact finite contraction bound is $`\theta_N=\max_i(Tw)_i/w_i`$. Computing $`Tw`$ needs one sparse multiplication and one triangular sweep.

There is a uniform bound as well on the genuine decorated state space. Write $`C_{\rm old}`$ for the original copy forcing, use $`(s-R)W\ge(s-\rho)W+h_{\min}`$, and note $`C_{\rm old}W\le\kappa(A+d)W`$ from the old-child and plus-weight bounds. Hence for $`\theta=\kappa(A+d)/[\kappa(A+d)+s-\rho]<1`$,
```math
C_{\rm old}W\le\theta (s+a-L^0)W.
```
Stopping the original particle at its first boundary and integrating proves the same comparison for copies before that boundary. A copied star's discounted weight on returning to a retained multicore is at most its starting $`W`$, by the star form of $`(s-R)W\ge0`$. Thus the intrinsic sweep and star substitution preserve this contraction. For the four-bin weight the implementation may use the maximum of the admissible local bin scores in each projected state to bound boundary weights. It checks the resulting finite norm $`\theta_N<1`$ directly.

For a rational or rounded numerical trial $`v`$ form the **exact** residual $`r=(1-I)^{-1}(h+Iv+Cv-v)`$. Then the fixed point of (4) satisfies
```math
\|u-v\|_w\le \|r\|_w/(1-\theta_N). \qquad\text{(5)}
```
Use (5) separately for the lower and upper source vectors, and propagate both corrections to each newborn row. An inner solver using ordinary integer rounding therefore does not determine the direction of any final inequality.

## 4. A positive dormant row at all excursions

Let $`D(S)=U_\beta(S)-U_A(S)`$, with identical old descendants in both responses. This equals $`(A-\beta)V_{s+A}^0 U_\beta(S)=(A-\beta)V_{s+\beta}^0 U_A(S)`$. For a retained multicore state, the first expression gives
```math
 (s+A+r_0)D(S)=(A-\beta)U_\beta(S)+\sum_{S'}r^0(S,S')D(S').
```
It needs one additional descending-size pass using the already bounded active values. A boundary $`D(S)`$ can use the second expression, the silent $`H^\pm_{\beta,A}`$ tables and the same upper old-grandchild weight. All terms are nonnegative.

For a dormant star, put $`\omega_b=\omega^{s+\beta}_{j,k}`$ and $`\omega_a=\omega^{s+A}_{j,k}`$. The exact positive decomposition is
```math
D(S_j)=\sum_{k\ge j}\{(\omega_b-\omega_a)[g(a(k))+U_\beta(V_k)]+\omega_a D(V_k)\}. \qquad\text{(6)}
```
The first term describes activation before promotion, and the second its continuation after promotion. The moment bounds for the difference use outward differences of the two rate-specific $`T`$ and $`M`$ intervals. The difference measures themselves are positive.

In the last term of (6), for $`k\ge K`$ the silent table gives the lower $`(A-\beta)g(\beta)(1-2/k)/[(s+A)(s+\beta)]`$. Its upper is constant. In detail, set $`v=s+\beta`$, $`\lambda=s+A`$, $`c_A=\kappa\alpha/(\lambda-q)`$ and $`X=\beta W(J_\beta)`$. The upper is
```math
(A-\beta)\Big\{[(g(\beta)+c_A\chi(K))/\lambda+c_A\chi(K)/(v-q)]/v
                      +c_A X/(v-q)\Big\}. 
```
Thus (6) uses the same four moments, including the old-child contribution. Once the four row entries are known, all subsequent dormant excursions use the scalar geometric series with ratio $`d_\beta`$. The code does not subtract two independently bounded complete rows to obtain $`d_A`$ or $`d_\beta`$.

## 5. Convergence as the multicore size cutoff grows

With exact star moments the only analytic boundary occurs when a multicore exceeds $`N`$. Finite ordinary exits use the $`\alpha W`$ upper. At a large star promotion, (2) gives an upper for the continuation at most $`u^++t^+/\lambda`$, bounded by a computable constant times $`W(V_k)`$ since $`W\ge h_{\min}/(A+d)`$. Let $`C_g`$ be the larger of $`\alpha`$ and the two resulting constants. All boundary widths are therefore at most $`C_g W`$. The same statement for the base dormant row follows from the stated constant and its old-child weight term.

Take $`\delta=0.013`$ with the basic profile or $`\delta=0.025`$ with the four-bin profile for $`s\ge0.565`$. The size drift inequalities give $`R(n^\delta W)\le s n^\delta W`$. The same holds for the base process started dormant, by the dormant size estimate in [the rank-one proof](dormant-return.md). Stop separately at the first oversize multicore of each lineage. The discounted exit-weight sum is at most $`(2/(N+1))^\delta W(J)`$ from a newborn, also when a star has grown far past $`N`$ before promoting. This gives a vanishing size error for each true row. At a fixed $`N`$, moment endpoints, outward rounding units and the residual (5) can all be refined to any tolerance before increasing $`N`$. This yields an effective scheme without a fixed large-tree continuation error.

## 6. Actual entries and separate replays

At $`u=0.5665`$, $`N=9`$, $`K=32`$ and $`P=96`$, the system has 954 multicore states and $`\theta_N\le0.2954199`$. For both true columns the residual correction is at most $`5\cdot10^{-14}`$ times the trial weight. The following enclosures are rounded outwards.

| True quantity | Lower | Upper |
|---|---:|---:|
| $`K(A,A)`$ | 0.52253633 | 0.96098322 |
| $`K(A,\beta)`$ | 0.07062970 | 0.53930328 |
| $`K(\beta,A)`$ | 0.83424723 | 3.52494776 |
| $`K(\beta,\beta)`$ | 0.29256522 | 3.21915168 |
| $`d_A`$ | 0.29956293 | 2.69982162 |
| $`d_\beta`$ | 0.19410746 | 2.81022683 |

At $`\mu`$ the Perron equations reduce to $`d_\beta=1-(1+d_A)k_\beta/(1-k_A)`$. The three lower bounds at $`u`$ thus give **$`d_\beta(\mu)<0.807760`$**. At $`0.54001`$ using the four-bin weight, the positive dormant row has upper $`4.001176330206<4.002`$ and $`\theta_N\le0.2970853`$. Hölder for the positive return gives $`d_\beta(s)\le4.002^p(0.807760)^{1-p}`$ for $`p=(u-s)/(u-0.54001)`$ between the anchor and $`u`$, using monotonicity above $`\mu`$. At $`0.5630`$ this gives $`d_\beta<0.9980`$, yielding the uniform enlarged-residual bound $`4800e^{0.5630t}W`$ by the potential proof. At $`0.564`$ it gives $`160e^{0.564t}W`$. Throughout the prior bracket $`d_\beta<0.8844`$. At the two Sabba endpoints it is below $`0.8522`$ and $`0.8423`$, respectively. Keeping 100 terms of the geometric series for the scalar first return then costs at most $`3.643\cdot10^{-6}`$ and $`1.062\cdot10^{-6}`$, respectively. This last error concerns the number of excursions. The boundary intervals in the table still dominate the total numerical uncertainty, and the endpoint signs remain open.

The two-dimensional determinant additionally gives finite comparative bounds for the scalar itself. Put $`r=r(K_s)`$ and $`b_*=0.5225363+0.0706297`$. The second eigenvalue is $`k_A+k_\beta+d_\beta-r`$, so $`\det(I-K)=(1-d_\beta)(1-\widetilde K_s)`$ gives
```math
\widetilde K_s-1=(r-1)\left[1+\frac{r-k_A-k_\beta}{1-d_\beta}\right].
```
Use $`k_A+k_\beta\ge b_*`$, the positive dormant ceilings just obtained, and the earlier binary anchor $`r(0.54001)<5.435`$. Log convexity with the prior interval gives $`0.95919<r(0.565615)<1.05819`$ and $`0.94667<r(0.565809)<1.04515`$. The resulting **prior-assisted** scalar intervals are $`(0.8468,1.2414)`$ and $`(0.8090,1.1746)`$, respectively. These bounds quantify the remaining scalar uncertainty near the two endpoints.

At $`s=1`$, $`N=7`$, the end-to-end calculation gives $`0.37148574<\widetilde K_1<0.50924616`$ from the current entries alone. The positive dormant row is at most $`0.18539285`$ there and the multicore contraction is at most $`0.2077579`$. The implementation accepts a scalar test, a tighter moment cutoff, additional iteration sweeps and more integer digits.\n\nA separate finite reference retains all star degrees through 18, including direct high-label births, and multicore size through 5. Solving the full 58-state killed system and the independently assembled 24-state Schur system gives 116 exact row-column equalities. [Schur implementation](excursion_rows.py) · [Finite reference](excursion_check.py) · [Exact contraction consequences](excursion_constants.py) · [Replay log](excursion-evidence.log)
