# A scalar renewal at a degree-two promoter

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Keeping the original nineteen marks, the expected discounted count of first degree-two direct births from $`J_{1/2}`$ is a scalar $`F(s)`$ with $`F(\mu)=1`$. The remaining newborn generations contract at $`s=0.5641`$, below the original lower bound $`0.565`$. The reviewed narrower interval improves the contraction factor at $`0.56564`$ to $`0.976433602362`$. A single killed scalar solve gives convergent rational lower bounds. At degree cap six and core cap three the discarded continuations still carry at least $`0.746245301832`$ of the critical return. Contraction checks · Killed scalar

## 1. The split and its critical equation

Use the original nineteen permanent birth marks $`a_i=1/i`$, $`2\le i\le20`$, with $`m(D)=\min(D,20)`$. In the positive fringe process write $`R`$ for the residual generator with all direct newborns set aside. Each intrinsic promotion still changes its entire current tree. The additional direct offspring at a planted root has rate $`b=L/D`$ and state $`J_{1/m(D)}`$. At degree two a positive birth rate occurs precisely on a virgin $`J_{a_i}`$, where it equals $`1/2`$. The other degree-two possibility has one old core child and no leaves. Set
```math
 h(S)=\tfrac12\mathbf1_{\{S\text{ virgin}\}},\qquad
 \widehat R f(S)= Rf(S)+ b(S)\mathbf1_{D_o>2} f(J_{1/m(D_o)}).
```
The full generator is $`\widehat R+h\,\mathrm{ev}_{J_{1/2}}`$. Every trigger emits a genuine $`J_{1/2}`$, independently of the promoter's own mark. The source particle remains alive in its old state in the direct-production term. This is the mean branching construction on planted fringes. [Nineteen-mark generator](../../fringe/finite-kernel.md)

Write $`K=K_s`$, $`k=K e_2`$ and $`C=K-k e_2^\top`$. Thus $`C`$ has zero column 2 and all its other columns are those of $`K`$. Wherever $`r(C)<1`$,
```math
 v_s=(I-C)^{-1}k,\qquad F(s)=v_s(2)
   =\sum_{n\ge0}e_2^\top C_s^n k_s
   =\mathrm{ev}_{J_{1/2}}(s-\widehat R)^{-1}h . \qquad\text{(1)}
```
The summand counts $`n`$ untriggered newborn generations followed by a triggered birth along one lineage. Positivity gives the identity by summing those mutually classified lineage contributions. The vector $`v_s`$ is strictly positive and $`K v_s-v_s=(F(s)-1)k_s`$, with $`k_s>0`$ in every row. Positive matrix comparison gives the same sign for $`F(s)-1`$ and $`r(K_s)-1`$, including equality. Consequently $`F(\mu)=1`$. The positive immediate term $`1/[2(s+2)]`$ is strictly decreasing, as are the positive chronological weights, so this scalar equation has a unique real root throughout the half-line certified below.

## 2. An explicit gap

Normalise the Perron vector of $`K_s`$ by $`g_{20}=1`$ and write $`r=r(K_s)`$. Each newborn row is an intrinsic Laplace transform in the root killing parameter $`a_i`$, with the same nonnegative forcing. Thus $`g_i`$ increases with $`i`$. There is a useful refinement in the first holding time. The original virgin has no old copies, intrinsic rate $`3/2`$, and direct payment $`g_2/2`$. After an intrinsic jump, its new shape has the same mark-independent intrinsic and copy forcing for every reassignment of its root mark. Hence $`(s+a_i+3/2)(K_sg)_i`$ increases with $`i`$ as well. In particular
```math
 k_i\ge\frac1{2(s+a_i+3/2)},\qquad
 (s+a_i+3/2)g_i\le s+1/20+3/2,\qquad
 C_sg\le\left(r-\frac{g_2}{2(s+1/20+3/2)}\right)g. \qquad\text{(2)}
```

At the residual discount $`\rho=0.51`$ the original weight obeys $`(\rho-R)W\ge h_0(c)\ge7/10`$, with $`W(J_{1/i})=(1/i+3/20)^{-1}\le5`$. The stopped drift inequality gives $`K_\rho\mathbf1\le 50/7`$ even at equality of the discount and the growth bound. Hölder on chronological paths through successive direct generations gives log convexity of $`r(K_s)`$. With the physical interval $`0.5656581227\le\mu\le0.5657979941<0.5658`$, putting $`s_0=0.56564`$ gives $`r(K_{s_0})\le (50/7)^{4/1395}<1.006`$. [Original weight](../../fringe/renewal.md) · Physical bracket

Using only the original intrinsic root, one has $`K_s(2,20)\ge0.022214534253`$ and $`(K_s\mathbf1)_2\ge0.405460197323`$ for $`0.51<s\le0.5665`$. For the latter, the virgin starts with $`b=1/2`$ and $`L^0 b\ge-b/6`$. In fact a leaf gain increases $`b`$ by $`c/[D(D+1)]`$ and a promotion decreases it by $`1/D`$. The gain rate is at least $`D-c`$, so $`L^0b\ge0`$ for $`c\ge2`$ and $`L^0b\ge-b/[D(D+1)]`$ for $`c=1`$. Nonnegative incoming rates preserve these inequalities. For the former entry, at $`D\ge20`$ use $`b/(s+1/2+\eta(c))`$, where $`\eta(1)=1/420`$ and $`\eta(c)=0`$ for $`c\ge2`$. Once $`c\ge2`$ it stays there. At $`D<20`$ solve backwards in $`(D,c)`$ using the minimum over extra gain rates $`0,(c-1)/2`$, zero direct toll, and denominator $`s+1/2+D-c+\xi+b`$. Integer rounding is downwards. At $`s=0.5665`$ a uniform $`\eta=1/420`$ gives $`0.022179480005`$. At $`s_0`$ the uniform and the varying terminals give $`0.022244698369`$ and $`0.022279887682`$, respectively. [Directed entry computation](../../outputs/finite-physical-trigger-final-bounds.log)

Taking $`a=0.0222`$, $`B=0.405`$, the row at label two gives $`r g_2\ge B g_2+a(1-g_2)`$, hence $`g_2\ge a/(r-B+a)`$. At $`s_0`$ this already gives $`r(C_{s_0})<0.998`$ by (2). Section 3 independently improves the same comparison to
```math
 g_2>0.11793629,\qquad C_{s_0}g\le0.978128g. \qquad\text{(3)}
```
The vector in (3) works simultaneously at every $`s\ge s_0`$ by monotonicity of the entries. These constants concern the **untriggered** generations and use the displayed prior on the physical exponent. [Rational bounds](physical_trigger_bounds.py)

## 3. A lower bound using six monotone label groups

Take representatives $`I=\{2,3,4,6,10,20\}`$ and $`j(d)=\max\{i\in I:i\le m(d)\}`$. Fix a common upper discount $`u\ge s`$ and an upper bound $`r(K_s)\le R_*`$. We maintain lower bounds $`x_i\le g_i`$ and $`H_i(d,c)`$ on the true $`U_sg`$ of a fringe having these root degrees and root mark at most $`1/i`$. The latter drop all rewards attributable to already present children. Initially $`x_{20}=1`$ and all the other lower bounds are zero.

Write $`\lambda_i=u+1/i`$. For a rate payoff $`f(d,c)`$ denote by $`J_{i,p}f`$ its lower discounted root table, with extinction discount $`\lambda_i`$, gain rate $`d-c+\xi`$, and promotion rate $`b=(d-c)/d`$. The extra gain rate lies in $`[0,(c-1)/2+p]`$. The first-event ratio is fractional linear in $`\xi`$ and the smaller endpoint value is used. This comparison permits adaptive physical inputs. A high-degree terminal is specified below.

From the previous tables define
```math
 w_{ij}=\left[J_{i,1/j}\min\{H_j(d,c),H_j(d+1,c)\}\right](2,1).
```
A newly promoted child of a degree-$`d`$ root has mark $`1/m(d)\le1/j(d)`$, and any later incoming arrow from that root has rate at most $`1/j(d)`$. While the child is present the two copy rates sum to $`1/m(d)`$. The worse of the unshifted and shifted child values is a valid lower payment for both copies. Conditioning at the child's birth and integrating its expected subsequent copy rate therefore adds at least $`w_{i,j(d)}/m(d)`$ to that birth's toll. Distinct root-child births give additive payments. The next tables use
```math
 t_i(d)=x_{j(d)}+w_{i,j(d)}/m(d),\qquad H_i^{\rm new}=J_{i,0}(b t_i).
```
At $`d\ge20`$ the toll is constant. Set $`\eta(c)=1/420`$ at $`c=1`$ and zero otherwise. The lower extension is $`H_i(d,c)=t_i(20)b/(\lambda_i+\eta(c))`$. For the child table the corresponding high lower integral is $`t_j^{\rm old}(20)b/[(\lambda_j+\eta(c))(\lambda_i+\eta(c))]`$. This follows from the drift of $`b`$ above, $`b(Y^+)\ge b(Y)`$ and the persistence of $`c\ge2`$. Finally replace $`x_i`$ by the maximum of its old value, $`H_i^{\rm new}(2,1)/R_*`$, and the current bound at the preceding representative. Induction proves all lower bounds, including when $`s< u`$. Rounding every division down to $`10^{-12}`$ gives (3) in 24 rounds. At $`(u,R_*)=(0.5665,1)`$ the same computation gives $`g_2(\mu)>0.11910171`$ and $`C_\mu g\le0.971864g`$.

One can dispense with the refined bracket on $`\mu`$. Use only $`0.565\le\mu\le0.5665`$ and the decreasing degree score with residual growth $`b_*=0.462`$, $`d_*=0.19`$, and $`h_{\min}=0.765931373`$. Its 189 degree inequalities imply $`C_{b_*}\mathbf1\le K_{b_*}\mathbf1\le1/[(1/20+d_*)h_{\min}]<5.45`$. With $`s_*=0.5649`$, $`p=(\mu-s_*)/(\mu-b_*)\le16/1045`$, Hölder row by row gives
```math
 C_{s_*}(g^{1-p})\le 5.45^p (0.971864)^{1-p}g^{1-p}
 <0.9981\,g^{1-p}.
```
This is a uniform weighted-norm contraction at a known number below the **original** lower endpoint. The minimum component of its unit-maximum comparison vector exceeds $`0.11910171`$. Degree certificate · [Integer comparisons](physical_trigger_bounds.py)

### 3.1 A stronger comparison from the original interval

Using all nineteen representatives in the same recursion, sixteen rounds give $`g_2(\mu)\ge0.126364920098`$ at the upper discount $`u=0.5665`$ and hence $`C_\mu g\le0.970147668298\,g`$. At $`s_0`$ the refined upper on the physical exponent gives $`g_2(s_0)\ge0.125103707001`$ and $`C_{s_0}g\le0.976433602362\,g`$. Each replay used one core for under three seconds. [Critical computation](../../outputs/finite-physical-trigger-nineteen.log) · [Fine computation](../../outputs/finite-physical-trigger-nineteen-fine.log)

For the lower discount $`b_*=0.462`$ put $`w_i=(1+4/20)/(1+4/i)`$ and $`\bar w_d=w_{m(d)}\mathbf1_{d>2}`$. In the decreasing degree weight of Section 3 write $`\alpha=1/h_{\min}`$, $`q_*=0.272`$, $`z_d=1/(1+0.19m(d))`$, and $`\lambda_i=b_*+1/i`$. The residual potential of the payoff $`b\,\bar w_D`$ is at most $`\alpha W`$ by the stopped drift. Count the original root and all the copies it emits during its lifetime. Each of its promotions contributes direct payoff $`\bar w_d`$ and creates a child of mark $`\ell=1/m(d)`$. The child's two future copy payments per unit time are at most $`\ell\alpha W(Y)`$, since $`W(Y^+)\le W(Y)`$. Even with the live parent's input, the child's mean intrinsic weight grows at most as $`e^{q_*t}W(J_\ell)`$. Given its birth, the future discounted copy payments up to the original root's death are therefore at most $`\alpha\ell W(J_\ell)/(\lambda_i-q_*)=\alpha z_d/(\lambda_i-q_*)`$. The compensator of intrinsic promotions gives the direct payoff at the same rate $`b`$. All the original root's children are born during that residence.

Maximise the degree recursion of Section 3 with toll $`T_i(d)=\bar w_d+\alpha z_d/(\lambda_i-q_*)`$. At $`D\ge20`$ the toll is constant. Its upper continuation is
```math
 T_i(20)\left(\frac1{\lambda_i}-\frac{1-b}{\lambda_i+1}\right).
```
Indeed $`L^0(1-b)\ge-(1-b)`$ because the gain rate is at most $`D`$ and promotions increase $`1-b`$. Upward integer divisions give $`\max_i(C_{b_*}w)_i/w_i\le3.488723`$. To interpolate, $`C_u g\le C_\mu g\le0.970148g`$. For $`b_*\le s\le u`$ let $`p=(u-s)/(u-b_*)`$. Apply Hölder to the joint measure of elapsed time and terminal label in each newborn row, with exponents $`1/p`$ and $`1/(1-p)`$. The result is
```math
 C_s(w^p g^{1-p})
 \le 3.489^p(0.970148)^{1-p}w^p g^{1-p}.
```
These comparison vectors are normalised to have maximum one and minimum at least $`0.126364920098`$. The upper factor at each displayed point also works for all larger real discounts with that same point's vector. [Anchor replay](../../outputs/finite-physical-trigger-final-bounds.log) · Exact power comparisons

| $`s`$ | $`p`$ | Upper factor using only $`0.565\le\mu\le0.5665`$ |
|---:|---:|---:|
| 0.5641 | $`24/1045`$ | 0.999090 |
| 0.5642 | $`23/1045`$ | 0.997867 |
| 0.5649 | $`16/1045`$ | 0.989348 |
| 0.5650 | $`3/209`$ | 0.988137 |
| 0.56564 | $`43/5225`$ | 0.980421 |

## 4. Positive chronological series and uniform error

For a finite marked fringe set $`d_s(S)=s+a_S+n(S)-1/D_o`$. Let $`t(S,S')`$ be the intrinsic rate plus the two old-child copy rates as in the nineteen-label sum. Put $`\widehat t=t+b(S)\mathbf1_{D_o>2}\mathbf1_{S'=J_{1/m(D_o)}}`$. Then (1) has the single positive expansion
```math
 F(s)=\sum_{\ell\ge0}\sum_{S_0=J_{1/2},\ldots,S_\ell}
 \left(\prod_{j<\ell}\frac{\widehat t(S_j,S_{j+1})}{d_s(S_j)}\right)
 \frac{h(S_\ell)}{d_s(S_\ell)} . \qquad\text{(4)}
```
The diagonal uses intrinsic changes and root death. Independent particle productions contribute to the numerator, including a production back to a type already seen. Nonexplosion and monotone truncation of chronological lineages justify the sum. A virgin visit through an old-child copy also earns the terminal rate $`h`$, with its actual root mark. The defining renewal is the **next degree-two direct birth**. Absolute domination by the sum at $`0.5641`$ gives analyticity on $`\Re s>0.5641`$.

For the sharper bound (3), $`k_{s_0}=(r g-C_{s_0}g)/g_2`$. Telescoping gives, uniformly for $`s\ge s_0`$,
```math
 0\le F(s)-\sum_{n=0}^{m-1}e_2^\top C_s^n k_s
 \le\left(1+\frac{0.006}{1-0.978128}\right)(0.978128)^m .
```
This error counts the number of untriggered newborn generations. It does not truncate intrinsic jumps.

The same gap gives a uniform chronological residual growth bound for every starting fringe. At $`s_0`$ set $`w_i=W(J_{1/i})`$, $`z=(I-C)^{-1}w`$ and $`u_j=(s_0-R)^{-1}[b\mathbf1_{m(D)=j}]`$. Then
```math
 \widehat W=W+\sum_{j=3}^{20}z_j u_j,\qquad
 W\le\widehat W,\quad \widehat W(J_{1/i})=z_i,\quad
 (s_0-\widehat R)\widehat W=(s_0-R)W\ge h_0(c_o).
```
Here $`z_i\le 5g_i/[g_2(1-0.978128)]`$, and $`\sum_j u_j\le(10/7)W`$. Thus $`\widehat W`$ is uniformly comparable to $`W`$, and the stopped semigroup comparison gives $`\widehat R_t W\le \mathrm{const}\, e^{s_0t} W`$. Replacing the parameters by the wider-bracket comparison gives the same conclusion at $`s_*`$. [Original finite-size drift](../../fringe/finite-kernel.md)

### 4.1 A single finite rational solve with an error tending to zero

For this paragraph set $`\theta=0.976434`$, $`m_0=0.12510370`$, $`A_0=1+0.006/(1-\theta)<1.255`$. The nineteen-group computation and the telescope in Section 4 bound the remaining generation tail by $`A_0\theta^m`$. This is below $`5.986\cdot10^{-10}`$ for $`m=900`$ and $`3.111\cdot10^{-11}`$ for $`m=1024`$. For all $`s\ge s_0`$ the same telescope starting in row $`i`$ gives $`v_s(i)\le A_0 g_i/m_0`$.

Let $`F_s(S)=[(s-\widehat R)^{-1}h](S)`$. First-generation decomposition gives $`F_s=U_s y`$, where $`y_2=1`$, $`y_i=v_s(i)`$ for $`i>2`$, and $`U_s y=(s-R)^{-1}[b\,y_{m(D)}]`$. Take the degree weight $`W'`$ with $`d'=0.18`$, $`q'=0.268`$, and $`h'_{\min}=81522947/125000000`$. Its two endpoint drift inequalities give $`(s_0-R)W'\ge h'_{\min}`$, hence
```math
 F_s(S)\le \frac{A_0}{m_0 h'_{\min}}W'(S).
```
Put $`\delta=0.11764`$ and $`V=n^\delta W'`$. The degree-weight size inequality gives $`RV\le s_0 V`$: the remainder $`\delta(q'-1/D_o)W'/n`$ in its divided drift is paid by $`h'(D_o,c_o)`$, since $`\delta q'\max h'/(1/20+d')\le h'_{\min}`$ and $`q'+d'+\delta=s_0`$. All 189 degree pairs, the high-degree tail and this last comparison were replayed in rational arithmetic. Size inequality · Replay

To include every untriggered newborn, set $`w'_i=W'(J_{1/i})`$, $`z'=(I-C_{s_0})^{-1}w'`$ and
```math
 \widehat V=V+2^\delta\sum_{j=3}^{20}z'_j u_j .
```
Every newborn has size two, so $`\widehat V(J_{1/i})=2^\delta z'_i`$ and $`(s_0-\widehat R)\widehat V=(s_0-R)V\ge0`$. Also $`z'_2\le[(1/20+d')(1-\theta)]^{-1}`$, by bounding $`w'`$ by $`g/[m_0(1/20+d')]`$ with the sharper denominator $`g_2`$ in place of $`m_0`$ before evaluating row 2.

Kill a lineage upon its first jump to size $`N+1`$. Let $`F_N(s)`$ be the solution at $`J_{1/2}`$ of the resulting finite scalar system $`(s-\widehat R_N)f=h`$. Its nonnegative resolvent exists since the full system has a semigroup bound at $`0.5641<s_0`$. Only intrinsic jumps can cross this size boundary. A stopped discounted $`\widehat V`$ account dominates the sum of all $`V`$ at the exits, and $`W'=V/(N+1)^\delta`$ there. It follows uniformly for $`s\ge s_0`$, $`N\ge2`$, that
```math
 0\le F(s)-F_N(s)
 \le\frac{A_0}{m_0h'_{\min}(1/20+d')(1-\theta)}
       \left(\frac{2}{N+1}\right)^{0.11764}
 <2837\left(\frac{2}{N+1}\right)^{0.11764}. \qquad\text{(5)}
```
For rational $`s`$, $`F_N(s)`$ requires only rational linear algebra. This proves effective convergence. The modest cutoffs in Section 5 use their own finite-sector upper bounds to measure the actual loss. Scalar, matrix and Schur replay

### 4.2 Eliminating every successive virgin holding time

Write $`J_i=J_{1/i}`$ in this paragraph. From $`J_i`$ there is a leaf gain of rate one and a promotion of rate $`1/2`$ to the bare root $`T_i`$ with its one child $`J_2`$. It has no pre-existing core children to copy. Hence
```math
 F_s(J_i)=\frac{1/2+F_s(J_i^+)+F_s(T_i)/2}{s+1/i+3/2}. \qquad\text{(6)}
```
It can be substituted at every subsequent arrival to any $`J_i`$. With missing successors set to zero this is exact on every geometry-killed system containing the newborns. Exact rational elimination of *all* virgin rows leaves the same response at every nonvirgin state at physical size caps three, four and five. At $`s=3/5`$, $`N=3`$, the scalar return is $`19047995/91316849`$. [Independent whole-tree checker](physical_trigger_check.py) · All three results

The numerator in (6) decreases with $`s`$. At $`J_2`$, $`(s+2)F(s)`$ is therefore nonincreasing. As benchmarks for an independent numerical calculation, the reviewed physical bracket implies $`F(0.565615)>1.00001680`$ and $`F(0.565809)<0.99999572`$. These two benchmarks use the prior on $`\mu`$. [Integer comparisons](physical_trigger_bounds.py)

## 5. How many old marks does a scalar computation need?

Fix a finite initial rooted geometry and distinguish its vertices before identifying unordered trees. In a term of (4) an old mark enters the holding denominator only at the lineage root. In a copy numerator it can enter only at the selected child, which becomes the new lineage root. Intrinsic rates and all new marks depend only on geometry. An untriggered direct newborn discards the old tree. Thus the old marks read by any individual path lie on a single ancestral chain. The same statement holds with geometry-only killing and for positive uniformised finite powers whose self coefficient depends only on geometry and the root mark. Use a common set of vertex-labelled walks for all assignments, including walks of zero weight. Every mixed difference in two incomparable old marks vanishes term by term. The argument works for the terminal score $`1`$ as well as $`h`$.

For each internal child whose parent has current degree $`d`$ use $`1/m(d)`$ as its baseline mark. Hold the exterior root label fixed. Write $`S|_v`$ for the assignment that retains the original marks only on the path from the root to $`v`$, and $`S|_o`$ for the all-baseline internal assignment. For each of the above response functions $`G`$,
```math
 G(S)=G(S|_o)+\sum_{v\ne o}\left[G(S|_v)-G(S|_{\mathrm{parent}(v)})\right]. \qquad\text{(7)}
```
For one walk the nonzero increments lie on its chain of read marks and telescope. Summing proves (7). Grouping terms by vertices expresses the right side as the sum at core leaves minus the branch-prefix values, each subtracted with multiplicity one less than the number of core children. Consequently for a nonnegative response its supremum over $`k`$-core geometries is at most $`\max(1,k-1)`$ times its maximum on the retained assignments. Two comparable old marks can have a nonzero mixed difference. A two-step exact check gives $`2372276870/53447875134621`$ for one such pair. This reduction applies to both the nineteen-label and the virgin-trigger reward. The latter also reduces the number of right-hand sides to one. Exact mark tests

For a degree cap $`H\le20`$, the possible birth-degree labels $`j_1\le\cdots\le j_b`$ of the $`b`$ children of a degree-$`d`$ vertex satisfy $`i+1\le j_i\le d`$. Before its $`i`$th promotion the parent needs a leaf in addition to its $`i-1`$ old children and its parent slot. Conversely, promote in sorted order after gaining leaves to each required degree, then gain to the final degree. If the parent is bare its existing child can send it a leaf. Child subtrees can subsequently be grown with the other arrows silent. This proves sufficiency recursively, and replacing off-path birth-degree labels by $`d`$ preserves the condition.

The following are counts per fixed exterior root label, at **exactly** $`k`$ cores and maximum degree six. To obtain the first column sum multisets of child unmarked types. For the chain column choose at most one exceptional child. Among subtrees of size $`l`$, its possible pairs of label and exceptional subtree number $`(d-1)C_l-U_l`$, with $`C_l,U_l`$ the chain and unmarked counts. The remaining children are an unmarked multiset. For unrestricted marks, add one multiset of subtrees at each birth-degree label $`j=2,\ldots,d`$, retaining at most $`j-1`$ cumulative children. [Count recurrence](physical_trigger_check.py) · Exact counts

| $`k`$ | Unmarked | Assignments retained by (7) | All reachable assignments |
|---:|---:|---:|---:|
| 1 | 5 | 5 | 5 |
| 2 | 25 | 75 | 75 |
| 3 | 185 | 1,435 | 1,775 |
| 4 | 1,530 | 28,805 | 49,410 |
| 5 | 13,915 | 587,315 | 1,512,760 |
| 6 | 133,751 | 12,041,401 | 49,165,401 |

In a geometry-killed scalar calculation started from $`J_{1/2}`$, just the five root labels $`2,\ldots,6`$ occur at $`H=6`$. A continuation supplied by an external high-degree parent can require all nineteen exterior labels. These counts are for the *original nineteen-mark* physical trigger. A different rule used in the previous low-sector computation assigned a common mark to low births. [Previous sector](low-sector/whole-sector.md)

### 5.1 Measured finite sectors

Cap every degree at $`H`$ and the core count at $`k`$. In the ordinary positive finite system $`Q=\widehat R_{H,k}`$ let $`\Lambda=k(H-1)+2`$ and $`P=(\Lambda I+Q)/(\Lambda+s_0)`$. A tree has at most $`k(H-1)+1`$ vertices and $`a_o\le1/2`$, so this choice gives $`P\ge0`$. Calculate $`x=(\Lambda+s_0)^{-1}\sum_{j<m}P^j h`$ and $`z=P^m\mathbf1`$ with a common integer denominator. In the reduced coordinates each row uses (7) on its successors and may have negative coefficients. The iterates still equal the pointwise values in the positive full system. For $`l=\max(1,k-1)`$, $`\tau=l\max_{\rm retained}z<1`$ implies
```math
 0\le F_{H,k}(s_0)-x(J_{1/2})
 \le z(J_{1/2})\,l\max_{\rm retained}x/(1-\tau).
```
This follows by iterating the positive *full* $`P^m`$. When using all states take $`l=1`$. The table gives outward endpoints for this finite response at $`m=96`$. All commands had a five-second timeout and ran on one core. [Integer sector code](physical_trigger_sector.py)

| $`H,k`$ | Representation | Rows | Nonzero entries | Wall seconds | Peak KiB | Killed scalar at $`0.56564`$ |
|---|---|---:|---:|---:|---:|---|
| 4,3 | Ancestral chains | 504 | 2,496 | 0.32 together | 19,444 together | $`[0.230673533259,0.230676493649]`$ |
| 4,3 | Full | 549 | 2,688 | combined | combined | $`[0.230673533259,0.230675013442]`$ |
| 4,4 | Ancestral chains | 4,341 | 25,620 | 1.83 | 46,080 | $`[0.231149350815,0.231224757262]`$ |
| 4,4 | Full | 5,520 | 30,564 | 2.07 | 47,916 | $`[0.231149350815,0.231174478897]`$ |
| 6,3 | Ancestral chains | 7,575 | 44,815 | 3.03 | 62,204 | $`[0.251490678253,0.253754698168]`$ |

Both $`x`$ and $`z`$ reconstructed exactly on all 549 full states of the first pilot, giving 1,098 integer identities. Since $`s_0<\mu`$ and $`F(\mu)=1`$, the $`(6,3)`$ boundary discards at least $`1-0.253754698168=0.746245301832`$ of the actual critical return. The full $`(4,4)`$ pilot similarly loses at least $`0.768825521103`$. The scalar construction and mark reduction make the finite solve smaller, but resolving the sign near criticality still requires substantially better boundary continuation. Paired pilot · Four-core chains · Four-core full · Six-degree pilot
