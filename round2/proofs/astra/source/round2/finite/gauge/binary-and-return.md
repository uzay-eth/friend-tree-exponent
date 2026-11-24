# Two birth marks and a scalar first return

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. Two actual newborn states suffice

Plant a fringe as in the [positive renewal construction](../../fringe/renewal.md). At the instant of promotion assign mark $`A=1/2`$ if the parent's degree is below 16, and $`\beta=1/16`$ otherwise. This ensures $`p=1/D_{\rm parent}\le a_Y`$ throughout the child's life. Actual internal arrows retain their physical reciprocal degrees. The rest of the positive residual and direct-birth construction is unchanged. The resulting newborn kernel has two rows and two columns.

Here is a short rational certificate of a gap. Take $`d=4/25`$, $`z(a)=a/(a+d)`$ and
```math
 W(S)=\frac{1}{a_o+d}\sum_{x\ {\rm core}}w_xh(D_x,c_x),\qquad
 w_x=\prod_{o<y\le x}z(a_y),\qquad h(2,1)=1.
```
Choose either $`(K,q)=(16,2/5)`$ or $`(64,49/125)`$. For all $`D\ge K`$ take $`h(D,c)`$ constant at the upward multiple of $`10^{-9}`$ of $`z(\beta)/q`$. Work in decreasing $`D`$ and, within each $`D`$, decreasing $`c`$. At every pair $`2\le D<K`$, $`1\le c\le D`$, set $`h(D,c)`$ to the upward multiple of $`10^{-9}`$ of
```math
 \max_{e\in\{0,(c+1)/2+d\}}
 \frac{(D-c+e)h(D+1,c)+(1-c/D)[h(D,c+1)+z(a(D))]}
 {q+D-c+e+1-c/D}.
```
The term at $`c+1`$ vanishes when $`c=D`$. The resulting $`h(2,1)`$ is below one. Raise it to one and check both endpoints in every row. In the larger table $`0.471910343\le h\le1.079138952`$. The explicit minima, rationals, and row counts are in [binary_weight.py](binary_weight.py).

The intrinsic contribution at each vertex is $`(D-c+e)\Delta_Dh+(1-c/D)[\Delta_ch+z(a(D))]`$ times its path weight divided by $`a_o+d`$. Physical input from core neighbours is at most $`c/2`$. The residual copies give $`dW-h(D_o,c_o)`$ plus
$`p_o\sum_Y\Delta_Dh(Y)/(a_Y+d)`$. At an immediate child the latter is precisely an additional gain coefficient $`p_o(a_o+d)/a_Y\le1/2+d`$ in its weighted row. Thus the checked interval covers both generators, giving
```math
 L^0W\le qW,\qquad RW\le(q+d)W-h(D_o,c_o).
```
For an initial centre one may add $`h(D,0)=2`$. Leaf gains then have zero increment there and promotion has increment $`h(D,1)+z(a(D))-2\le0`$. Adding a leaf changes $`W`$ by a bounded factor, obtained by checking $`\max(1,h(D+1,c)/h(D,c))`$. Consequently centre immigration satisfies the weighted bound in the renewal construction as well. The two certificates give $`q+d=0.56`$ and $`q+d=0.552`$, respectively, both strictly below the proved lower bound $`0.565`$ for $`\mu`$.

All four newborn entries are positive. To obtain the high-degree label, the original newborn can take sufficiently many consecutive root gains and then promote. The threshold and transfer argument of [finite-kernel.md](../../fringe/finite-kernel.md) now gives
```math
 r(K_\mu^{(2)})=1,\qquad
 \mu=\text{largest real zero in }(0.552,\infty)
 \text{ of }\det(I_2-K_s^{(2)}).
```
The superscript counts the two possible marks. Both are genuine states. For the size cutoff take $`\delta=13/1000`$. Since $`W/n\le h_{\max}/(\beta+d)`$, the checker verifies $`\delta qh_{\max}/(\beta+d)\le h_{\min}`$. The usual size-weight inequality then gives $`R(n^\delta W)\le0.565\,n^\delta W`$. The two-label construction saves matrix dimension while giving a smaller certified size power than the nineteen-label construction.

## 2. A positive common part of all newborn rows

The following argument works for any finite birth-mark scheme with maximal mark $`A=1/2`$. For a bounded nonnegative test $`v`$, the starting root mark $`a`$ changes only its original particle's constant killing rate. Thus
```math
 (K_s v)(a)=\int_0^\infty e^{-(s+a)t}T_{s,v}(t)\,dt.
```
Here $`T_{s,v}(t)`$ includes the reward at time $`t`$ and the future discounted offspring of copies emitted at $`t`$. Its intrinsic law, including the marks at all descendant births, is independent of $`a`$. Put $`k_s=K_s(A,\cdot)`$ and $`D_s=K_s-\mathbf1 k_s`$. All entries of $`D_s`$ are nonnegative. They are Laplace transforms of nonnegative measures in total chronological time because the extra factor at the initial time is $`e^{-at}-e^{-At}\ge0`$, and the remaining copy delays have their ordinary positive laws. One implementation couples an additional hazard $`A-a`$ to the initial particle. The common contribution is produced before the additional killing time, and the excess comes from later survival. Subsequent residual families use the physical birth-mark rule.

Whenever $`r(D_s)<1`$, define
```math
 F(s)=k_s(I-D_s)^{-1}\mathbf1
     =\sum_{n\ge0}k_sD_s^n\mathbf1.
```
For $`v=(I-D_s)^{-1}\mathbf1>0`$ we have $`K_sv-v=(F(s)-1)\mathbf1`$, which gives equivalence of both strict signs and equality by positive comparison. At the critical point let $`g>0`$ be the Perron vector with maximum one. Its minimum is $`g_A`$, since each row is nonincreasing in the starting mark. The identity $`D_sg=r(K_s)(g-g_A\mathbf1)`$ shows $`r(D_\mu)\le1-g_A<1`$. Continuity gives a neighbourhood on which the scalar equation $`F(s)=1`$ detects the exponent. This qualitative statement only needs the positive renewal theorem.

For the binary scheme there is no matrix inversion beyond a scalar denominator. Write $`k_j=K_s(A,j)`$, $`d_j=K_s(\beta,j)-k_j`$. The first row of $`D_s`$ is zero, so $`r(D_s)=d_\beta`$ and
```math
 F(s)=k_A+k_\beta\frac{1+d_A}{1-d_\beta}.
```
The four entries still contain genuine fringe continuations.

## 3. A quantitative scalar expansion with the nineteen-mark kernel

For the nineteen-label scheme the separate degree-weight calculations and row enclosures give $`r(K_s)<1.02`$, $`g_A>0.12884`$ throughout $`[0.565,0.5665]`$. Thus $`D_sg\le\theta g`$ with $`\theta=1.02(1-0.12884)=0.8885832`$. These quantitative inputs come from the residual drafts under independent review. For a nonnegative term bound, $`\mathbf1\le g/g_A`$ and $`k_sg=r(K_s)g_A`$ give
```math
 0\le F(s)-\sum_{n=0}^{N-1}k_sD_s^n\mathbf1
 \le\frac{1.02\,\theta^N}{1-\theta}.
```
This is below $`6.79\,10^{-5}`$ at $`N=100`$ and $`5.34\,10^{-9}`$ at $`N=180`$. For any trial vector $`w`$ the same cancellation gives
```math
 |F(s)-k_sw|\le \frac{1.02}{1-\theta}
 \,\|\mathbf1-(I-D_s)w\|_\infty.
```
The residual comparison controls errors in the scalar equation once the true row actions have been enclosed.

The entries of $`D_s`$ are positive Laplace transforms, so $`r(D_s)`$ is log convex by Hölder on paths through $`n`$ generations and the $`n`$-th root limit. At $`s=0.4481`$ it is bounded by $`r(K_s)<4.062`$ from the separate anchor calculation. At $`s=0.556`$, log convexity and $`4.062^{90}\theta^{1079}<1`$ give $`r(D_{0.556})<1`$. Consequently the positive scalar series is analytic on $`\Re s>0.556`$ by domination on real compact subintervals. Its real critical root is unique there by positivity and strict decrease of the immediate common newborn contribution. This construction has one regeneration class and retains the original marks inside each excursion. Whether every physical birth can instead carry mark $`1/2`$ is considered in constant-gauge.md.

## 4. Four auxiliary birth-degree bounds give residual exponent 0.54

In the proof remember the bin $`i\in\{2,4,8,16\}`$ containing the parent's degree at each core birth. Use the largest bin threshold at most that degree. The true mark is $`a_i=A`$ for $`i=2,4,8`$ and $`a_{16}=\beta`$. These extra records have no effect on any clock, copy probability or birth label. Given the physical tree and the actual marks, any admissible assignment of auxiliary bins has the same projected dynamics. The projected newborn kernel therefore still has exactly two rows and two columns.

Use the weight of section 1 with $`d=7/50`$, $`q=2/5`$ and local score $`h_i(D,c)`$. Set the four newborn coefficients $`u_i=h_i(2,1)`$ to $`(1,963/1000,931/1000,974/1000)`$. A newborn contributes
$`\chi(D)=a_{i(D)}u_{i(D)}/(a_{i(D)}+d)`$ to its parent's local row. At each nonroot the ordinary input from its parent is at most $`1/i`$, and that from its children is at most $`(c-1)/2`$. If this vertex is an immediate child of the particle root, the residual plus-copy adds an effective input $`p(a_o+d)/a_i\le(A+d)/(i a_i)`$. The row inequality in section 1 therefore uses
```math
0\le e\le(c-1)/2+[1+(A+d)/a_i]/i,\qquad z(a(D))\ \hbox{replaced by }\chi(D).
```
The silent root lies within this bound, as does the time-dependent root input of an old child while it remains attached. At $`D\ge40`$ take the constant upward rounding of $`\chi(40)/q`$. The backward recursion gives $`h_i(2,1)\le u_i`$ before imposing the specified four values. Rechecking the 6,232 rational endpoint inequalities gives
```math
h_{\min}=0.501218471,\qquad h_{\max}=1.072734084,\qquad
 L^0W\le0.4W,\quad RW\le0.54W-h_{i(o)}(D_o,c_o). 
```
The initial physical centre has no birth degree. Use any fixed bin there, take $`h(D,0)=2`$ before its first promotion, and use its physical inputs from children thereafter. The checker verifies $`h_i(D,1)+\chi(D)\le2`$ for each bin. The root-to-child accounting and bounded leaf-gain ratio then give the centre estimate of section 1. Treat any newborn that starts from actual mark $`A`$ as originating in bin 2 if the earlier birth degree is unavailable. Children already present retain their genuine birth records.

For $`\delta=1/40`$, the rational test $`\delta qh_{\max}/(\beta+d)\le h_{\min}`$ gives $`R(n^\delta W)\le0.565\,n^\delta W`$ by the same size calculation. In particular the two-by-two determinant characterises $`\mu`$ throughout $`(0.54,\infty)`$, and the finite cutoff gives an explicit two-sided computation. The stronger result for the two-label scalar first return uses its true row differences and the prior bracket. [Auxiliary-bin checker](virtual_weight.py) · [Scalar denominator](binary-row-bounds.md)
