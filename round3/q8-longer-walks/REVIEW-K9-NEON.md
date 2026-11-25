> **Public checkpoint export note:** Edited export of `walk-astra/REVIEW-K9-NEON.md` (frozen original SHA-256 `c3bb25ad02ec277d16d12bd3a0e5bfc38abd35d6853798f892672a9a349181ff`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

WRONG — Lemma 1(2)'s displayed general polynomial. HOLDS — the corrected certificate, the all-tree nine-step drift, the stochastic theorems, and the specified stock obstructions.

# Version and verdicts

Astra reviewed Neon's NINE-STEP.md at SHA-256 `8d444e33dd93392276c316da805f2ef319b8b29c89b1ce076a694ef861594974`. The verifier was exported by a Git blob read at `f2e72034e5fee54fcef207e6d1ddae6c7ba33ebb:q8/verify_k9.py`. The SHA-256 of its bytes is `dd1b8bfa5c918a08f81df7ed23fb1ef6fe3bcaa0df9091d13a61675feed97d23`. Both exports are read-only, named by digest under the reviewer's board directory `pinned-k9/`. The named earlier proof files matched the frozen versions used in the seven-step review. The replay evidence and dependency hashes are in RESULT-K9-REVIEW.md.

| Scope | Verdict | Reason |
|---|---|---|
| Lemma 1(2), line 61, claimed exact polynomial | **WRONG** | The expression is for walk length $`2m+3`$. An actual degree-17 branch gives a discrepancy of $`(16/17)^5`$ at nine steps. The pinned program evaluates the intended degree-five polynomial. |
| Lemma 1(1), unrestricted $`J_0,m`$ | **UNSUPPORTED by the printed induction; HOLDS by reconstruction** | The induction stops at $`J_0\ge m`$ in line 77. The first-visit identity below proves the conclusion with constant one for all stated cutoffs and lengths. |
| Corrected analytic-to-certificate passage at $`k=9`$ | **HOLDS** | Killed branch walks give the verifier's recurrence, all branch classes are covered, and formal multiplication bounds their polynomial degrees. All 90 mixed Bernstein coefficients are positive. The general polynomial in line 61 cannot serve as that bridge. |
| Proposition 1 on every finite tree | **HOLDS** | Exact increment, first-visit bound, endpoint and Bernstein inequalities. For the specified weights the sharp drift bound is $`1+(800/801)X_v`$. |
| Theorem 1, positive attained maximum, uniform degree remainder and founder support | **HOLDS** | Analytic conditional floor, terminal-core attempt, aggregate second moment and simultaneous martingale tails, including labels born later. |
| Theorem 2, total share, leaves, color, hub neighbors and conditional attachment law | **HOLDS** | Bounded-time planting, prefix cut and upward convergence with a prefix-independent constant. |
| Theorem 3(1) at nine steps and Theorem 3(2) at odd $`k\ge11`$ | **HOLDS in the stated stock classes** | Explicit sequences of leaf-rich trees force the necessary inequalities. They impose no conclusion on existence or absence of linear hubs at these longer lengths. |

The underlying process uses a fresh uniform start and exactly the indicated non-lazy number of steps at every birth. Seed trees are fixed, finite and have at least two vertices.

## 1. Independent first-visit identity on every tree

For a vertex $`x`$ let $`g(x)`$ be the change in the focal score $`X_v`$ if the next vertex attaches at $`x`$. Thus $`0\le g\le1`$, $`g(v)=0`$, and $`g`$ is supported on odd distances from $`v`$. At distance one the increment is one on a leaf and $`s(y)`$ on a nonleaf $`y`$. The only other positive increment can occur at distance three, on a leaf child of an outward pure star of a split neighbor $`y`$. It is $`\gamma(y)`$. Growing a leaf child of $`y`$ changes L to shallow star and leaves the deep count fixed. Adding a new child directly at a shallow star also preserves its class. An attachment cannot make an already deep branch shallow. This checks the cases where the sign or the number of promotions could otherwise be wrong.

Let $`P`$ be the walk transition matrix and $`K_B`$ the substochastic matrix killed on first reaching $`v`$ in the component $`B`$ rooted at $`y\sim v`$. Degrees are the full-tree degrees, including the cut edge. Put $`f(x)=1/d_x`$. With one unit of starting mass per vertex,

```math
N=\sum_x (P^{2m+1}g)(x).
```

The start at $`v`$ contributes $`(P^{2m+1}g)(v)`$. From $`B`$, paths avoiding $`v`$ contribute $`\sum_{x\in B} d_x g(x)(K_B^{2m+1}f)(x)`$ by killed reversibility. A first visit to $`v`$ at an odd time leaves an even number of steps and hence zero final reward. The total starting mass first visiting $`v`$ at time $`2i+2`$ is

```math
\sum_{x\in B}K_B^{2i+1}(x,y)/d_y=(K_B^{2i+1}f)(y).
```

Conditional on that visit the remaining reward is $`(P^{2m-2i-1}g)(v)`$. Therefore there is an exact identity

```math
N=(P^{2m+1}g)(v)+
\sum_B\left[\sum_{x\in B}d_x g(x)(K_B^{2m+1}f)(x)
+\sum_{i=0}^{m-1}(K_B^{2i+1}f)(y)(P^{2m-2i-1}g)(v)\right].
```

All root-started rewards in this expression are at most one. A branch consisting of a single leaf has zero killed contribution at positive times, including in the two-vertex seed. The remaining terms are precisely the two kinds of within-branch reward and the $`m`$ source terms needed in $`\Psi`$. This proves the root elimination with additive constant one whenever the killed branch values satisfy the stated envelopes, for any $`J_0,m`$. Alternatively the missing high-degree induction in the manuscript follows at once from $`c_p\le j`$, $`a\le1-1/j`$ and $`H_p+c_p a\le1+j(1-1/j)=j`$.

## 2. The branch polynomial and the manuscript discrepancy

For a degree-only high row let $`t=1-1/j`$, $`a=j^{-1}\sum_{z\ne v:z\sim y}1/d_z`$ and $`w=js`$. Conditioning outward in the killed branch bounds its odd walk values at $`y`$ by $`h_0=a`$, $`h_{p+1}=t-a+ah_p`$. Pure deep rows at small degrees have the same envelope on the smaller interval $`0\le a\le (j-1)/(2j)`$. For $`m\ge2`$ direct summation gives

```math
\sum_{p=0}^{m-1} h_p+w h_m =
(m-1+w)t-\frac1j\sum_{i=1}^{m-2}(m-1-i+w)a^i
-s a^{m-1}+(1-w)a^m+w a^{m+1}. \qquad (1)
```

Indeed $`h_p=t-(1/j)\sum_{i=1}^{p-1}a^i-a^p+a^{p+1}`$ for $`p\ge1`$. For $`m=1`$ the polynomial is $`wt+(1-w)a+w a^2`$. Line 61 shifts $`m`$ to $`m+1`$ throughout. In an actual branch whose root has degree 17 and whose sixteen outward neighbors are leaves, $`s=0`$, $`a=t=16/17`$ and the recursion is exact. At $`m=4`$ its polynomial is $`t+t^2+t^3+t^4`$; the displayed polynomial adds $`t^5=1048576/1419857`$. This also explains why a verifier checking just the degree-five identity and the two intended endpoints passes.

For $`0\le w\le1`$ the derivative in (1), after zero coefficients are omitted, has negative terms only below the nonnegative top terms. Equivalently for positive $`a`$ divide by $`a^{m-2}`$: the negative lower powers then increase or are constant, and the top powers increase. Hence a maximum on any $`[0,b]`$ occurs at an endpoint. This includes cases where the derivative has no positive zero. At $`m=1`$ the polynomial is already nondecreasing. At zero and $`t`$ the endpoint values printed in line 62 are correct.

**A separate numerical typo.** In lines 99 and 103, $`2s(2,1)=2(2388-1946)/1000=221/250`$, whereas the manuscript prints $`111/125`$. Both the inequality $`2s\le1`$ and the stated margins hold. The other factors and all 30 displayed finite endpoint margins in that table agree with the intended recurrence.

## 3. Coverage of the continuous small-degree rows

Treat an ordinary leaf as a shallow star with $`\theta=0`$. For each shallow outward child $`z`$ use $`\theta_z=1-1/d_z`$ and let $`U_z(p)`$ represent $`d_z(K_B^{2p}f)(z)`$ before bounding. Its initial value is one and its update is the preceding odd value at $`y`$ plus $`\theta_z U_z(p-1)`$. For $`b`$ deep children let $`a_B=j^{-1}\sum_B1/d_z`$. A deep child with degree $`d_z`$ contributes at most $`(1-1/d_z)+h_{p-1}/d_z`$ before division by $`j`$, because the onward killed expectation is at most one. Thus upper envelopes propagate forward by

```math
h_0=a_B+j^{-1}\sum_{\rm shallow}(1-\theta_z),\quad U_z(0)=1,
```

```math
U_z(p)=h_{p-1}+\theta_z U_z(p-1),\qquad
h_p=b/j-a_B+a_B h_{p-1}+j^{-1}\sum_{\rm shallow}(1-\theta_z)U_z(p).
```

The branch bound is $`\sum_{p<m}h_p+js h_m+\gamma\sum_z\theta_z U_z(m)`$. For S the last term counts all $`d_z-1`$ promotable leaves. This recursion also works when several S stars have different degrees. For all admissible parameters the operations propagate upper bounds monotonically in the previous walk values.

There are three degree-two child types L, S, B and six unordered degree-three pairs LL, LS, SS, LB, SB, BB. L and LL are the two constant rows. B and BB are the small pure endpoint rows. S, LS and LB each have one parameter; SS and SB each have two. In particular SB includes the deep reciprocal parameter separately from the star parameter. The degree-only rows require no classification of their actual children. Real parameters in the certified rectangles contain all allowed integer degrees, including the limiting large-star end.

With $`\theta=(1+x)/2`$ or $`a_B=y/6`$ as applicable, the total polynomial degrees obey $`\deg U_z(p)\le p`$ and $`\deg h_p\le p+1`$ by induction. The complete nine-step residual thus has total degree at most five. This supplies an analytic reason that a $`(5,5)`$ tensor basis represents the entire function. The independent checker multiplies formal polynomials over the rationals and converts monomials directly to Bernstein coefficients, using $`\binom{k}{i}/\binom{5}{i}`$ in each variable. No fitted interpolation or off-grid sampling is needed for this check. The three six-coefficient lists match individually. Their minima in order are $`4037/256000`$, $`88873/562500`$ and $`103643/20250000`$. The two 36-coefficient tables have minima $`11486981/40500000`$ and $`3209143/20250000`$. All are strictly positive.

Together the first-visit identity, corrected endpoint calculation, the two all-leaf rows, these five mixed rows and the constant tail prove (G9) for every finite tree. There is no reliance here on extrapolation from small trees.

**Sharp bound for the given stock.** Across the 30 finite pure endpoints, the two constant rows and the mixed rows bounded using Bernstein coefficients, the largest upper ratio $`\Psi/\Phi`$ is $`1074337/1075680`$, at $`j=12,a=t`$. This is only $`7/95735520`$ below $`800/801`$. In the tail $`s=0`$, both endpoint values are below four and $`4/G_{17}=800/801`$. Consequently

```math
N\le1+(800/801)X_v.
```

Stars force the additive constant one. For sharpness of the slope take $`C`$ degree-$`j`$ nonleaf spokes each with $`j-1`$ outward leaves, and $`C^2-C`$ leaf neighbors at the focal root. For fixed $`j\ge17`$ the upgrade increment vanishes. The root even walk values tend to zero as $`C\to\infty`$, the branch odd walk values tend to $`t^{p+1}`$, and unrolling leaf returns gives $`N/C\to t+t^2+t^3+t^4`$. Equivalently in the first-visit identity the root reward at each fixed odd time tends to one on this leaf-rich family. Here $`X_v=C G_{17}`$. Taking $`j\to\infty`$ proves that any uniform affine drift for these same weights requires slope at least $`800/801`$, even with an arbitrary finite additive constant. This certifies an $`O(n^{800/801})`$ uniform degree remainder by the same proof. Optimality of the actual stochastic remainder is open.

## 4. Almost-sure and total-share assertions

**Convergence, positivity and the moving maximum.** Reversibility gives $`q_v\le d_v`$; five leaf choices give $`q_v\ge d_v(1-C_v/d_v)^5\ge d_v-5C_v`$. Hence $`d_v(n)/n`$ is a bounded supermartingale and converges for all fixed labels. Set $`Y=X+250/111`$ at the stated $`\beta=999/1000`$. From size $`s`$, the conditional expected $`C`$ is bounded by $`Y_sP_j`$ for $`P_s=1`$, $`P_{j+1}=(1+\beta/j)P_j`$. Differencing $`P_j/j`$ bounds the tail loss in the degree fraction by $`5000Y_s/s`$. This implies the $`d/(2s)`$ floor when $`d\ge41000(C+1)`$. The quadratic degree supermartingale bounds $`E[Z_v^2\mid\mathcal F_s]\le2(d/s)^2`$, giving conditional probability at least $`1/8`$ when this floor applies.

At every size at least three the nonleaf-induced tree has a terminal vertex, allowing its singleton case. At this vertex $`d\ge2,C\le1`$. Until its next direct hit or leaf promotion the hit weight is at least $`1/16`$. The promotion weight is at most $`d e_4\le1+4C\le5`$ by unrolling exact root returns. Other updates preserve the required bounds. The first such event occurs almost surely by the divergent lower hazard, and conditioning over its possible times gives hit probability at least $`1/81`$. At most 82000 consecutive such hits yield a qualifying history. The resulting uniform conditional lower probability of some hub at every starting size and upward convergence imply probability one. At the edge seed a specified founder star to size $`s`$ has probability $`1/(s-1)`$; the conditional expectation floor on it tends to one, proving the founder assertion for every $`\varepsilon>0`$.

The shifted stock has increments in $`[0,1]`$. For $`F(u)=u(u+1)`$, $`F(u+h)-F(u)\le2(u+1)h`$ on this interval. Aggregate over old vertices and add the newborn whose shifted stock is at most seven. Division by $`\prod(1+2\beta/j)`$ and addition of the remaining deterministic immigration tail give a nonnegative supermartingale since $`2\beta>1`$. Thus $`\max C_v=O(n^\beta)`$ almost surely. Each present label's predictable tail loss is now uniformly $`O(n^{\beta-1})`$, while the centered tail has squared jump budget less than $`1/n`$. The two-sided maximal Azuma bound at $`\sqrt{8\log n/n}`$, summed over $`v\le n`$, is at most $`2n^{-3}`$. Borel–Cantelli covers all late sizes. For unborn $`s>n`$ apply the estimate at birth, obtaining $`Z_s\le1/s+O(s^{\beta-1})`$ uniformly. This handles late rivals and forces $`Z_s\to0`$, so the positive supremum of the $`Z_s`$ is attained. It prevents exchanging coordinate convergence and a maximum without justification. The proof works with $`\beta=800/801`$ and shift $`1/\beta`$ for the remainders.

**Total share.** Fatou on normalized excess degrees gives $`\sum Z_v\le1`$. Write $`R_0=41000`$, $`M=32R_0`$ and $`J=R_0(M+1)`$. Starting at a bounded edge $`vy`$, stop on $`J`$ hits at $`v`$, a hit on $`y`$ or a current leaf of $`v`$, or size $`2n`$. While active the target weight of $`v`$ is at least $`(M+J)^{-4}M^{-5}`$ from the single alternating walk; bad weight is at most $`2M+J`$, including if $`y`$ is itself a leaf. Put guaranteed good marks of length $`a/n`$ at the bottom of a fresh uniform and all possible bad marks in the top $`b/n`$, with $`a=(M+J)^{-4}M^{-5}/2`$, $`b=2M+J`$. Exactly $`J`$ bottom marks and no top marks in $`n`$ trials has probability at least $`(a/2)^Je^{-2(a+b)}/J!`$ for all sufficiently large $`n`$. The marks may be sampled to the horizon even on early termination. Conditional on success by the stopping size, the rich floor is at least $`1/(4n)`$. This gives a deterministic $`p_M/n`$ per eligible endpoint; simultaneous success across endpoints is unnecessary.

Outside any fixed old prefix, denote the internal edges by $`E`$, the boundary edges with a nonleaf outside endpoint by $`B`$, low–low internal edges by $`E_{LL}`$, highs by $`H`$, and total rich-high degree by $`D_R`$. High nonrich degrees are below $`2R_0C_v`$. In $`\sum_H C_v`$, high–high edges contribute at most $`2|H|`$, and old neighbors at most $`B`$. Low nonleaves with at least two high neighbors contribute at most $`2|H|`$ by the bipartite forest bound. Each low nonleaf with just one high neighbor has another edge, costing at most $`2E_{LL}+B`$ in total. With $`M|H|\le2E+B`$ this yields

```math
E\le328002E_{LL}+2D_R+328001B.
```

This covers empty outside sets and arbitrary numbers of leaves attached to the prefix, which are excluded from $`B`$. At least $`2E_{LL}/M`$ distinct bounded-edge endpoints are disjoint from the rich highs. Summing the conditional floors gives $`n E[S_{>N}\mid\mathcal F_n]\ge\kappa(E-328001B)`$ for a fixed $`\kappa>0`$ independent of the old prefix size $`N`$. The exact identity $`E=n-1-\sum_{v\le N}d_v+|E(T_n[1,N])|`$ and $`B\le\sum_{v\le N}C_v=o(n)`$ give the stated prefix inequality by upward convergence. Increasing $`N`$ after countable intersection proves unit share.

Every leaf has one neighbor, giving leaf dominance by fixed-prefix sums of $`d_v-C_v`$. The excess-degree vectors have limiting total mass one and coordinate convergence, so they converge in $`\ell^1`$. Red and blue each have degree sum $`n-1`$; excess shares are consequently $`1-R`$ and $`R`$. Since the newborn copies the fresh start's color, the ordinary two-color urn has the stated beta law. On a hub the leaf fraction at that hub tends to one, so the conditional promotion lower bound $`(L_v/d_v)^5/n`$ forces infinitely many distinct promotions. Eventually their labels exceed any fixed cutoff and they have degree two and at most one nonleaf neighbor. The uniform terminal floor at those random histories and upward convergence for the existence of an adjacent hub past the cutoff yield infinitely many direct hub neighbors. Finally the five-excursion error bound gives coordinatewise attachment convergence, and unit total share gives total-variation convergence.

## 5. Obstructions and the equality boundary

Take $`D=C^2`$ with $`C`$ identical nonleaf branches. In every test the even moments at the focal root tend to zero for fixed walk horizon since the chance to leave its leaf excursions is $`O(m/C)`$. The leaf-promotion weight divided by $`C`$ tends to the sum of the first $`m`$ limiting odd branch moments. Prescribing a finite birth order has positive probability because alternating an odd walk along any incident edge gives the desired target positive probability. An arbitrary finite additive constant disappears upon division by $`C`$.

For $`k=9`$, take degree $`j\ge3`$ with all outward children deep and both exterior star sizes diverging. The zeroth odd branch moment tends to zero, the succeeding fixed odd moments to $`t=1-1/j`$ since a walk sent out two steps reaches the large leaf-stars. The leaf source tends to $`3t`$ and the hit weight tends to $`j-1`$. Hence for $`\beta\le1`$ the necessary degree-only rows give $`G_j/(j-1)-G_{j+1}/j\ge3/j^2`$. Dropping the nonnegative terminal weight yields $`G_3\ge6\sum_{j=3}^{16}j^{-2}=173672066749/86572886400>2`$. On a length-two spoke, the limiting source is $`15/16`$, the $`y`$ hit is $`1/16`$, and the outward leaf hit is $`31/32`$. Its transition L to S(2) has nonnegative increment by hypothesis. Thus the required coefficient exceeds one. This argument concerns scores using degree alone at *all* $`j\ge3`$ and the specified comparison at degree two.

For odd $`k=2m+1`$, a degree-two spoke with a single deep child has limiting source $`(m-1)/2`$ and hit weight one. The hit goes to degree three with one deep child. Thus $`\beta w_{2,1}\ge(m-1)/2`$ by the required comparison. A degree-two spoke with a single shallow star of diverging degree has vanishing source and hit weights. The combined promotion weight of its leaves tends to one: in the exact shallow update $`\theta\to1`$ and the finitely many $`h_p\to0`$. A direct hit at $`y`$ has increment bounded below by $`-1`$, so its vanishing rate suffices even without monotonicity for that transition. Thus $`\beta\ge w_{2,1}-1`$. Nonnegativity and the first test force $`\beta>0`$. Combining proves $`2\beta(1+\beta)\ge m-1`$, excluding $`\beta<1`$ at eleven and $`\beta\le1`$ from thirteen onwards. These tests exploit the common weight on L and every shallow star and the specified degree-three comparison. Coarser or broader interpretations about every local method or about hub nonexistence do not follow.

**Reviewer extension at $`k=11`$.** In the full class with weights $`w_{j,b}`$ by degree and number of deep children (degree-only cutoffs allowed), the equality case $`\beta=1`$ already fails if the two extra comparisons $`w_{4,1}\ge w_{3,1}`$ and $`w_{4,2}\ge w_{3,2}`$ hold. These are among the monotone degree increments used in §1. At equality the first two tests force $`w_{2,1}=2`$ and, retaining the hit increment in the deep test, $`w_{3,1}=2`$. Two deep children at degree three force $`w_{3,2}\ge 4(2/3)=8/3`$ because a hit at $`y`$ cannot decrease its weight. Now give a degree-three spoke one deep child and one shallow leaf-star, with both exterior scales diverging. The odd branch moments tend to zero at $`p=0`$ and $`1/3`$ at $`p\ge1`$. This follows directly by conditioning on which child was chosen: the deep child supplies even walk value one after at least two steps and the shallow child supplies zero. The limiting source is $`4/3`$, the $`y`$ hit weight is one, and the shallow-to-deep promotion weight is $`1+4/3=7/3`$ by the exact shallow update. The required row would give

```math
2=w_{3,1}\ \ge\ 4/3+(w_{4,1}-2)+(7/3)(w_{3,2}-2)\ \ge\ 26/9.
```

This contradiction has margin $`8/9`$. It uses extra monotonicity beyond the sparse assumptions of Theorem 3(2). It concerns this proof method alone.

The new exact calculations, commands and version pins are recorded separately in RESULT-K9-REVIEW.md. The stochastic conclusions use the analytic arguments above. This is a cross-family AI audit for human verification.
