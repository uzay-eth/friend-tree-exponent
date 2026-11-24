# Astra Claim 3 — Sharp mean counts and the limiting pointwise degree law

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 3: Claude and Gemini; GPT standalone audit; inputs include Round 1 growth, the marked renewal and spatial weights (Claims 2 and 4). See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Statement and input ledger

Start with a **fixed single undirected edge** $`T_2`$. At size $`n`$, choose a vertex uniformly, then choose one of its neighbours uniformly, and attach the new vertex as a leaf to that neighbour. There is no distinguished seed endpoint in the rule. Write $`D_x`$ for total degree, $`C_n`$ for the number of vertices with degree at least two, and $`N_{\ge k}(n)`$ and $`N_{=k}(n)`$ for the cumulative and exact degree counts. In particular $`N_{\ge2}(n)=C_n`$. All means are with respect to this growth law, not a random initial seed.

**Claim.** Let $`\mu`$ be the logarithmic mean exponent of $`C_n`$. For every fixed integer $`k\ge2`$ there are constants $`c_k>c_{k+1}>0`$, with $`c_k\downarrow0`$, such that

```math
\mathbb E N_{\ge k}(n)\sim c_k n^\mu,
\qquad
\lim_{n\to\infty}\frac{\mathbb E N_{=k}(n)}{\mathbb E C_n}=p_k:=\frac{c_k-c_{k+1}}{c_2}.
```

Here $`c_2=\lim_{n\to\infty}n^{-\mu}\mathbb E C_n>0`$, $`\sum_{k\ge2}p_k=1`$, and for some $`0<A<\infty`$,

```math
\sum_{j\ge k}p_j\sim A k^{-\mu},\qquad p_k\sim\mu A k^{-1-\mu}\quad(k\to\infty).
```

First take $`n\to\infty`$ at **fixed** $`k`$; only then let $`k\to\infty`$ in the resulting probability law. This is a ratio of expectations, not an expectation of a random ratio, an almost-sure empirical-degree limit, or an estimate uniform in $`n,k`$.

**Repository inputs and scope.** The one prior growth theorem used is [existence.md, §0, Theorem 1 and §7, “Proof of the continuous-time limit”](../../proofs/existence.md): for the edge-start rule above, both $`\log \mathbb E C_n/\log n`$ and $`t^{-1}\log\mathbb E_{T_2} C(t)`$ converge to the same $`\mu\in[0,1]`$, using rate-one clocks per vertex. The [motif bracket, §1, Theorem 1 and the remark after Corollary 2](../../proofs/bracket-motifs.md) bounds edge-start mean core counts by positive multiples of products with exponents $`113/200`$ and $`1133/2000`$ (from size three). With the growth theorem this gives $`113/200\le\mu\le1133/2000`$; below we only need $`51/100<\mu\le1`$. The finite rational check in [renewal.md, §3, (4)](../fringe/renewal.md) concerns every integer $`2\le D<100, 0\le c\le D`$ for the explicit function in §2 below; we include a complete reproducible specification of that arithmetic and an analytic bound for all larger $`D`$.

The compact-renewal construction is from `bfcbbe5`, [fringe/renewal.md, §§1–4](../fringe/renewal.md). Commit `8e4db79` supplies the [pole and fixed-size transfer, fringe/asymptotics.md, §§1–4](../fringe/asymptotics.md), the [local cell operator and tilt, operator.md, §§1–2](../operator.md), and the [birth-age and flux proof, degree-law.md](../degree-law.md). These arguments are derived below. The constant-factor part of Claim 2, the minimum or full spectral-radius characterization of the local operator, and sharper numerical brackets are not needed.

## 1. The original clock and the two seeds

Give every current vertex an independent rate-one Poisson clock, with a fresh clock for each new vertex. At a ring the vertex chooses a uniform neighbour. At size $`n`$, the next wait has law $`\mathop{\mathrm{Exp}}\nolimits (n)`$ and is independent of the uniform ringing vertex and its uniform neighbour choice. By the fresh-increment property applied at each birth, the sequence of waits is independent of the discrete jump chain. This is the standard embedding of [existence.md, §1.1, Lemma 1.1 and §1.3, Lemma 1.3](../../proofs/existence.md); it does not change the discrete attachment law. The size is the pure-birth chain with birth rate $`n`$ in state $`n`$, hence is nonexplosive. From $`r`$ initial vertices it is a sum of $`r`$ independent geometric variables on $`\{1,2,\ldots\}`$ with parameter $`e^{-t}`$: attribute births to the clock's family, or solve the pure-birth forward equations starting at one and superpose $`r`$ independent families. Thus for $`q_t=e^{-t}`$,

```math
\Pr_r(N(t)=r+j)=\binom{r+j-1}{j}q_t^r(1-q_t)^j,\quad
\mathbb E_rN(t)=re^t,\quad \mathbb E_r N(t)^2=r(r+1)e^{2t}-re^t.
```

All fixed-time moments are finite. These are the assertions of [existence.md, §1.3, Lemma 1.4](../../proofs/existence.md), also seen from the geometric laws. We use later that the size-only chain begun at $`r`$ satisfies $`e^{-t}N(t)`$ a martingale, by the conditional mean formula.

At the first jump, after a wait $`\mathop{\mathrm{Exp}}\nolimits (2)`$, the edge becomes the three-vertex path $`T_3`$. Its centre depends on which endpoint was chosen, but the unlabelled tree is deterministic and future growth from its centre has the same law in either case. From here until the fixed-size transfer we often start **at time zero from $`T_3`$**, rooting the physical tree at its centre $`o`$. That centre is already core and stays core. For a current core $`x`$, write $`L_x`$ for its leaf count and $`c_x`$ for its number of core neighbours, so $`D_x=L_x+c_x`$. Its degree rises at rate $`L_x+\sum_{y\sim x,\ y\text{ core}}1/D_y`$. It promotes a leaf at rate $`b_x=L_x/D_x`$: the leaf obtains a new leaf and becomes core, while the promoter's degree stays the same and $`(L_x,c_x)`$ becomes $`(L_x-1,c_x+1)`$. Each new core therefore starts at $`(L,c)=(1,1)`$, unlike the initial centre $`(2,0)`$. Degrees never decrease.

If $`F_r(t)=\mathbb E_{T_r}C(t)`$, then $`F_2(t)=\int_0^t2e^{-2u}F_3(t-u)\,du`$. Both means are nondecreasing, so $`F_2(t)\le F_3(t)\le F_2(t+1)/(1-e^{-2})`$. In particular $`t^{-1}\log F_3(t)\to\mu`$. The centre contributes exactly one to $`F_3(t)`$, so the same logarithmic limit holds for the mean number of noninitial cores.

## 2. Positive renewal on planted descendant types

At the promotion of a leaf $`v`$ by $`u`$, give $`v`$ the permanent birth mark $`a_v=1/D_u`$ at that instant. Because the promoter's degree never falls, at all later times $`p_u:=1/D_u\le a_v`$. For a noninitial core $`v`$, let $`S_v`$ be its descendant subtree, including leaves and birth marks. Its root keeps a **core parent slot**, counted in both $`D`$ and $`c`$ but not filled by a vertex of $`S_v`$. The auxiliary type space contains finite oriented descendants with these slots and the internal inequalities $`p_u\le a_v`$; the root mark also allows $`0`$. Put $`\mathcal A=\{1/2,1/3,\ldots\}\cup\{0\}`$, a compact subset of the real line, and let $`J_a`$ be the root of mark $`a`$ with a single leaf and the parent slot. The mark $`0`$ is an auxiliary limit, never a physical birth. Denote the root's $`L/D`$ and $`1/D`$ by $`b_S`$ and $`p_S`$. Let $`Y`$ run over all immediate core-child subtrees with multiplicity, and let $`Y^+`$ add a leaf at the root of $`Y`$.

In the intrinsic generator $`L^0`$, the slot sends nothing; the root's choice of that slot makes no internal change. The only input from outside $`S_v`$ missing in the physical evolution is a leaf gain at its root from the actual parent, at the parent's reciprocal degree. Hence in a sum of old scores the parent accounts for $`p_S[f(Y^+)-f(Y)]`$ per core child. A root promotion also creates the new summand $`f(J_{p_S})`$. Adding and subtracting one $`a_v f(S_v)`$ per noninitial vertex reorganizes the signed transfer into a positive branching system:

```math
Rf(S)=L^0f(S)-a_S f(S)+\sum_Y\big[(a_Y-p_S)f(Y)+p_Sf(Y^+)\big],
\qquad Bf(S)=b_S f(J_{p_S}).
```

A residual particle evolves intrinsically, dies at rate $`a_S`$, and makes additional copies of the indicated types at the nonnegative rates in $`R`$, retaining the source at copy times. The full process additionally makes direct $`J`$ births at rate $`b_S`$. Denote the two first-moment semigroups by $`R_t`$ and $`P_t`$. The direct birth can have its own clock; independence here only realizes the required **first-moment generator**, not the physical joint law. For $`M_t f=\mathbb E_{T_3}\sum_{v\ne o,\ v\text{ core}} f(S_v(t))`$, the only remaining row is immigration from the actual centre (which has no parent slot):

```math
I_t f=\mathbb E_{T_3}\left[b_o f(J_{p_o})+\sum_{Y\text{ at }o}((a_Y-p_o)f(Y)+p_o f(Y^+))\right].
```

Indeed each $`a_Y f(Y)`$ from a parent's row, including a centre row, cancels the corresponding death term in the child's row and leaves the directed gain. Intrinsic changes in each old descendant tree and new summands at promotion are then accounted for exactly. At $`T_3`$ the sum $`M_0 f`$ is empty. The convolution consequent to this generator identity is

```math
M_t f=\int_0^t I_u P_{t-u}f\,du.\qquad\text{(1)}
```

**Weight and justification of (1).** Put $`d=3/20`$, $`q_0=9/25`$, $`\rho=q_0+d=51/100`$, and let $`h(0)=2,h(1)=1,h(2)=17/20,h(3)=39/50,h(c)=7/10`$ for $`c\ge4`$. For each core $`x`$ in a planted type, $`w_x`$ is the product of $`a_y/(a_y+d)`$ along its path from the root (excluding the root). Define

```math
W(S)=\frac{1}{a_S+d}\sum_{x\text{ core}}w_x h(c_x).
```

Leaf gains leave $`W`$ unchanged. A promotion of a degree-$`D`$, core-degree-$`c`$ vertex has rate $`1-c/D`$ and changes its weighted contribution, after dividing out the old path factor, by $`h(c+1)-h(c)+1/(1+dD)`$. For all integer $`D\ge\max(2,c)`$,

```math
(1-c/D)[h(c+1)-h(c)+1/(1+dD)]\le q_0 h(c).\qquad\text{(2)}
```

For $`D\ge100`$ the left side is at most $`1/16<q_0(7/10)`$. The remaining arithmetic, including the ordinary centre case $`c=0`$, is exhausted by:

```python
from fractions import Fraction as F
def h(c):
    return [F(2), F(1), F(17,20), F(39,50)][c] if c < 4 else F(7,10)
for D in range(2,100):
    for c in range(D+1):
        assert F(D-c,D)*(h(c+1)-h(c)+1/(1+F(3,20)*D)) <= F(9,25)*h(c)
```

This gives $`L^0W\le q_0 W`$. Factoring the first edge on each path gives $`\sum_Y a_YW(Y)=(a_S+d)W(S)-h(c_{\rm root})`$. Also $`W(Y^+)=W(Y)`$, so $`RW\le\rho W-h(c_{\rm root})`$. Since $`1<W(S)\le |S|/d`$ and $`BW\le1/d`$, the stopped, exponentially discounted sum of particle weights gives $`R_tW\le e^{\rho t}W`$ and $`P_tW\le e^{(\rho+1/d)t}W`$. To see that removing the stopping is valid, stop at the $`m`$th event, then use the weight bound to bound the chance of ever reaching $`K`$ simultaneous particles before a fixed time by $`O(1/K)`$. On a path with bounded count, the death intensity is bounded by half that count, so there are only finitely many deaths and therefore only finitely many births (including direct births). Each individual has nonexplosive intrinsic growth, since its jump rate is at most its size. This excludes explosion and Fatou yields the unstopped bounds. The same promotion estimate for the weighted sum rooted at the actual centre, without a root denominator, gives mean at most $`2e^{q_0t}`$. The immigration copy terms then sum in weight to that sum minus the centre term; thus $`I_tW\le 1/d+2e^{q_0t}`$.

For $`|f|\le C W`$, the first-event equation at each finite type gives the backward derivative $`\partial_v P_vf=(R+B)P_vf`$. This use is pointwise: finitely many initial successor multisets, initial rate $`O(|S|+1)`$ and successor weights $`O(|S|+1)`$ give a derivative bounded by $`O_t((1+|S|)^2)`$ for $`0\le v\le t`$. Continuity and the derivative follow from the first-event integral and the time-local weight bound. In the physical jump formula apply the time-dependent sum with scores $`P_{t-u}f`$, for $`0\le u\le t`$. Stop at a finite size; the physical Yule moments dominate sums of the scores, derivatives and jumps on this interval when the stop is removed. The backward derivative cancels the $`(R+B)`$ row, leaving just the integrated immigration and proving (1). This is the unbounded-system justification of [renewal.md, §3, (5)–(6) and final paragraph](../fringe/renewal.md).

## 3. Compact reproduction, critical pole and no other pole on the line

For complex $`z`$ with $`\Re z>\rho`$, work with the weighted supremum space $`\|f\|_W=\sup_S|f(S)|/W(S)`$ and define

```math
G_zf(S)=\int_0^\infty e^{-zt}R_tf(S)\,dt,\quad
U_z\varphi=G_z[b_S\varphi(p_S)],\quad
jF(a)=F(J_a),\quad K_z=jU_z.
```

Also let $`\widehat I_z F=\int_0^\infty e^{-zt}I_tF\,dt`$ for $`\Re z>q_0`$. The integrals are componentwise; the weight estimates justify the Taylor expansions in operator norm by termwise integration against a slightly smaller real discount. In particular no strong continuity of $`R_t`$ in weighted supremum norm is required. We have $`|U_z\varphi|\le W\|\varphi\|_\infty/(\Re z-\rho)`$.

**Compactness and positivity.** In a family started at $`J_a`$, separate its original particle and the families started by its copies. The original particle's intrinsic evolution and the marks of copied descendants do not depend on $`a`$; its survival to time $`t`$ supplies the factor $`e^{-at}`$. Thus $`K_z\varphi(a)`$ is the integral with factor $`e^{-(z+a)t}`$ of the intrinsic expectation of

```math
b_S\varphi(p_S)+\sum_Y[(a_Y-p_S)U_z\varphi(Y)+p_SU_z\varphi(Y^+)].
```

For $`\|\varphi\|_\infty\le1`$, its absolute expectation is at most $`1+e^{q_0t}/(\Re z-\rho)`$, using the same path-weight cancellation and the initial value $`h(1)=1`$. Multiplying this bound by $`t e^{-\Re z\, t}|a-a'|`$ proves a uniform Lipschitz bound in the root mark, also at the accumulation point $`0`$. Arzelà–Ascoli on compact $`\mathcal A`$ gives compactness on $`C(\mathcal A)`$. The identical separation for $`jG_z f`$ with $`f=1`$ or a root-degree indicator uses $`f`$ as the original toll and $`G_z f`$ for the copies; $`f`$ on the intrinsic original particle does not depend on its root mark. Thus $`jG_z f`$ is continuous in $`a`$ and analytic in $`z`$ in this space.

For real $`s>\rho`$, the kernel has a positive reproduction measure in birth time and mark. Each row charges each $`1/j`$ ($`j\ge2`$): the original root can keep its first leaf, gain $`j-2`$ leaves before any promotion or death, and then remain in a state with direct birth rate $`(j-1)/j`$ over a positive-length interval. The continuous positive cone has nonempty interior, and this support property maps every nonzero nonnegative continuous function into the interior (continuity also rules out a function supported only at $`0`$). Strong Krein–Rutman therefore gives a simple positive spectral radius $`r(K_s)`$, a strictly positive eigenfunction with positive minimum, and a positive dual eigenmeasure. The radius decreases strictly: for $`s_1<s_2`$ the positive difference on the $`s_2`$ eigenfunction dominates a strictly positive multiple of it, since births occur at positive times, and iteration gives $`r(K_{s_1})>r(K_{s_2})`$. Continuity follows by bounding nearby $`K_s`$ on the fixed positive eigenfunction from above and below by $`(r\pm\varepsilon)`$ times that function; the weighted supremum comparison or iteration yields the corresponding radius bounds.

**Selecting $`\mu`$.** Group a full family by generations of direct $`J`$ births, retaining the residual occupation between them. For real $`s>\rho`$, the discounted score $`1`$ at newborns is $`\sum_{n\ge0}K_s^n jG_s1`$, with nonnegative terms. The forcing dominates $`1/(s+1/2)`$ by survival of the original particle. For $`r(K_s)<1`$, the series converges in norm, and from arbitrary types the discounted total is $`G_s1+U_s\sum_{n\ge0}K_s^n jG_s1=O_s(W)`$; (1) and the immigration bound then imply finiteness of the physical count transform. For $`r(K_s)\ge1`$, the forcing dominates a positive multiple of the positive eigenfunction, giving divergence at every newborn mark. Immigration directly from the centre has positive expected rate, so the physical transform diverges by Tonelli. Before its first own ring the centre has only leaves, so this rate is positive on any finite time interval. The logarithmic growth limit for $`M_t1=F_3(t)-1`$ places the transform's abscissa at $`\mu>\rho`$. Accordingly continuity gives $`r(K_\mu)=1`$, its unique real crossing on $`(\rho,\infty)`$. This is the part of [renewal.md, §4, (7)–(9)](../fringe/renewal.md) needed here.

For $`f=1`$ or $`f=f_k(S)=\mathbf1_{D_{\rm root}\ge k}`$, group its occupation in the same manner and use (1). For $`\Re z>\mu`$ the absolute series is dominated by the real-discount series at $`\Re z`$, so

```math
\mathscr M_f(z):=\int_0^\infty e^{-zt}M_t f\,dt
=\widehat I_z\left(G_z f+U_z(I-K_z)^{-1}jG_z f\right).\qquad\text{(3)}
```

**Residue.** At $`z=\mu`$ let $`g,\ell`$ be the positive eigenfunction and dual eigenmeasure, with $`\ell(g)=1`$. The isolated simple eigenvalue has a local analytic branch $`\lambda(z)`$: integrate the resolvent around a small circle enclosing only $`1`$, use the norm inverse series to obtain analytic projections, and note that a sufficiently small perturbation of a rank-one projection has rank one. On the complementary subspace the inverse of $`I-K_z`$ is analytic, since it exists there at $`\mu`$. Differentiate $`K_zg_z=\lambda(z)g_z`$ with $`\ell(g_z)=1`$ to obtain $`\lambda'(\mu)=\ell(K'_\mu g)<0`$, because differentiation multiplies the positive reproduction measure by minus its positive birth time. The pole in (3) is simple, with residue

```math
d_f=\frac{\widehat I_\mu(U_\mu g)\,\ell(jG_\mu f)}{-\ell(K'_\mu g)}.\qquad\text{(4)}
```

The first factor is positive by direct centre immigration and $`U_\mu g(J_a)=g(a)>0`$. For each fixed degree $`k\ge2`$, the original particle can reach exactly $`k`$ by leaf gains alone and then hold there for positive time while alive. Thus $`\ell(jG_\mu f_k)>0`$ and $`\ell(jG_\mu(f_k-f_{k+1}))>0`$. All factors are finite by the weight estimates. No Tauberian step will be applied to the nonmonotone exact-degree test.

**Nonlattice argument.** Fix $`\tau\ne0`$. If $`K_{\mu+i\tau}v=v`$ for a nonzero continuous complex $`v`$, then $`|v|\le K_\mu|v|`$. The measure $`\ell`$ has full support: for nonzero continuous $`f\ge0`$, $`\ell(f)=\ell(K_\mu f)>0`$ by strong positivity. Hence equality holds everywhere and $`|v|`$ is a positive multiple of $`g`$. Equality in the triangle inequality in each row would require all nonzero integrands to have the same phase. The initial $`J_a`$ has intrinsic jump rate $`3/2`$ (one leaf gain and a rate-$`1/2`$ promotion), death rate $`a`$, and positive direct birth rate $`1/2`$ before any intrinsic jump. It therefore supplies a mark-$`1/2`$ component with time density at least $`e^{-(\mu+3/2+a)t}/2`$ in the **positive time-and-mark reproduction measure discounted at the real rate $`\mu`$**. The complex kernel at $`\mu+i\tau`$ integrates against that measure with the Fourier factor $`e^{-i\tau t}`$. On the mark-$`1/2`$ component, $`e^{-i\tau t}v(1/2)`$ cannot have constant phase on an interval. This contradiction excludes eigenvalue 1, and the compact-operator Fredholm alternative gives invertibility there. Invertibility persists locally by the inverse series. Combining this with (4) shows that $`\mathscr M_f(z)-d_f/(z-\mu)`$ extends analytically near **every** point of $`\Re z=\mu`$. No uniform strip or estimate at large imaginary part is used.

## 4. Tauberian transfer, then fixed-size transfer

For each fixed $`k\ge2`$ the cumulative count of noninitial cores is nondecreasing, as the statistic counts all new cores that qualify and old degrees never fall. Apply the Laplace-function Wiener–Ikehara theorem to $`S(x)=M_{x/\mu} f_k`$. The transform is $`\mu\mathscr M_{f_k}(\mu z)`$ on $`\Re z>1`$; the scaling leaves residue $`d_{f_k}`$ at 1. The nonnegative nondecreasing $`S`$ is locally integrable and log-linearly slowly decreasing. Local analytic extension of the pole-subtracted transform implies the theorem's local pseudofunction boundary condition: after multiplication by a compactly supported smooth cutoff the boundary function is integrable and its Fourier transform vanishes at infinity. The theorem gives $`S(x)\sim d_{f_k}e^x`$, or $`M_t f_k\sim d_{f_k}e^{\mu t}`$. The initial centre contributes at most one, so the full $`T_3`$ mean has the same asymptotic. This applies also to $`f_2=1`$ on the planted core types.

The first birth from the edge is independent of the subsequent evolution. As $`N_{\ge k}=0`$ before that birth, convolution against its rate-two density, with the bounded normalized $`T_3`$ mean as integrand, yields by dominated convergence

```math
e^{-\mu t}\mathbb E_{T_2} N_{\ge k}(t)\ \longrightarrow\ d_k^*=\frac{2}{2+\mu}d_{f_k}.\qquad\text{(5)}
```

Here boundedness follows from convergence at infinity and finite-time Yule domination. Neither $`d_{f_k}`$ nor $`d_k^*`$ is the fixed-size amplitude.

Let $`a_{n,k}=\mathbb E N_{\ge k}(n)`$. Independence of size and jump chain gives the exact negative-binomial mixture

```math
\mathbb E_{T_2} N_{\ge k}(t)=q_t^2\sum_{n\ge2}(n-1)(1-q_t)^{n-2}a_{n,k}.
```

For $`x\to\infty`$, set $`\nu_x=x^{-(\mu+2)}\sum_{n\ge2}(n-1)a_{n,k}\delta_{n/x}`$. For any fixed $`s>0`$ put $`1-q_t=e^{-s/x}`$. The Laplace transform of $`\nu_x`$ at $`s`$ is $`e^{-2s/x}(xq_t)^{-(\mu+2)}`$ times the left side of (5) with its factor $`e^{-\mu t}`$. As $`xq_t\to s`$, it tends to $`d_k^* s^{-(\mu+2)}`$. This is the transform of the density $`d_k^* y^{\mu+1}/\Gamma(\mu+2)`$ on $`y\ge0`$ (by substitution in Euler's gamma integral). For interval convergence, tilt $`\nu_x`$ by $`e^{-y}`$ and normalize by its mass, which tends to $`d_k^*>0`$. The tilted probability transforms tend to $`(1+s)^{-(\mu+2)}`$, those of the gamma probability of shape $`\mu+2`$ and rate one. For instance, pushing forward by $`u=e^{-y}`$ gives convergence of every integer moment on $`[0,1]`$ including mass; polynomial approximation gives weak convergence there. Undoing this coordinate change and the tilt on a bounded interval with positive endpoints is valid since the gamma law assigns no endpoint mass. Thus

```math
\nu_x([a,b])\longrightarrow \frac{d_k^*}{\Gamma(\mu+2)}\int_a^b y^{\mu+1}\,dy\qquad(0<a<b<\infty).
```

For integer $`x=n`$ and $`0<h<1`$, monotonicity of $`a_{n,k}`$ bounds it below by the $`(m-1)`$-weighted average on $`[(1-h)n,n]`$ and above by the analogous average on $`[n,(1+h)n]`$. After division by $`n^\mu`$, the numerators in these ratios are $`\nu_n([1-h,1])`$ and $`\nu_n([1,1+h])`$, and the denominators are $`n^{-2}\sum(m-1)`$ on the respective size intervals. These tend to the corresponding integrals of $`y`$. Taking $`n\to\infty`$ and then $`h\downarrow0`$ gives

```math
c_k=\frac{d_k^*}{\Gamma(\mu+2)}=\frac{2d_{f_k}}{(2+\mu)\Gamma(\mu+2)}.\qquad\text{(6)}
```

Equation (4) gives $`c_k>c_{k+1}>0`$. For $`k\to\infty`$, at each finite time every residual family from a newborn has finitely many particles of finite degree. The weighted moment bound permits dominated convergence in its occupation, the time integral (since $`\mu>\rho`$), and then the finite measure $`\ell`$ (since $`W(J_a)\le 1/d`$ uniformly). So $`\ell(jG_\mu f_k)\downarrow0`$, and $`c_k\downarrow0`$. For fixed $`k`$, subtract the two cumulative-count asymptotics; their unequal leading coefficients imply $`\mathbb E N_{=k}(n)\sim(c_k-c_{k+1})n^\mu`$. The claimed ratio follows, and the sum of the $`p_k`$ is one by telescoping.

## 5. A local law with this degree marginal

To extract a pointwise tail, we need a probability law of the environment of a core, not just its degree marginal. The following is the portion of [operator.md, §§1–2, (2)–(5)](../operator.md) used here.

Consider finite or countably infinite locally finite rooted **core** trees $`G`$, marked by their leaf counts $`L_x\in\{0,1,\ldots\}`$, with $`D_x=L_x+c_x\ge2`$ and $`c_x`$ the number of neighbours in $`G`$. A singleton core is allowed. There are no permanent birth labels or seed orientations in this space. Its local topology is agreement of marked balls; one more radius determines boundary degrees when needed. Compatible balls determine a locally finite marked rooted tree. The resulting space is Polish: it is the closed space of consistent finite-ball records in the countable product of countable discrete spaces, with rooted isomorphism recovered by consistent finite ball isomorphisms. The ball cylinders determine probabilities.

**Cell dynamics.** At initialization give each core its own cell label; a newborn core inherits the promoter's cell label. Cell membership is only a device for tracing descendants of an old core, not the permanent numerical birth mark used in the renewal. A cell grows at rate at most its number of cores, so its core population from one is dominated by a Yule process from one. If $`H_o(t)`$ is the cell of initial root $`o`$, define

```math
Q_t f(G,o)=\mathbb E_{G,o}\sum_{v\in H_o(t)} f(G_t,v).
```

For infinite old trees construct the dynamics with independent proposals on each oriented old edge $`x\to y`$ of rate $`1/c_x(0)`$ (no clocks on nonexistent edges); accept with probability $`c_x(0)/D_x(t-)`$. Each old edge thus sends at the physical rate $`1/D_x`$. All other transitions act within cells: rate-$`L_x`$ leaf production at $`x`$, and rate-one core clocks using lengths $`1/D_x`$ for each internal neighbour and $`L_x/D_x`$ for promotion, with the old-edge choices omitted from these clocks. In a finite union of old cells the endogenous rate is at most total size and there are finitely many incoming old-edge proposals on a bounded interval, so there is no explosion there.

Compare truncations of the old ball of radius $`R`$ retaining all missing degree slots and suppressing incoming arrows across the cut. The marked old ball of radius $`R+1`$ determines the missing slots at the boundary of the $`R`$-ball. Only old-edge proposals communicate between cells, so a discrepancy reaching an old cell $`w`$ from outside must travel along a path with increasing proposal times. For a test of radius $`k\le R`$ about a member of the root cell, only cells based in the old $`k`$-ball can matter. Erase loops from an influence path to such a $`w`$, keeping the proposal times in order. Its length is still at least $`j=R-k+1`$, since it starts outside the $`R`$-ball. Count only its **last $`j`$ steps**, which give a simple directed path ending at $`w`$ with ordered proposals in $`[0,t]`$. Each such path has expected tuple count $`t^j/j!`$ times the product of the reciprocals of the old sender degrees in the core. Reversing paths identifies this product with the reversed simple random-walk path probability times $`c_w(0)/c_{\rm start}(0)\le c_w(0)`$. Summing gives failure probability at most

```math
\min\left(1,\frac{t^j}{j!}\sum_{w\in B_k(o)}c_w(0)\right).
```

The bound is summable in $`R`$ and controls every larger truncation with the same old proposals; thus on every fixed finite set of cells and bounded time interval the processes stabilize almost surely. This gives the local process with the usual future Poisson drivers, also at local stopping times such as a promotion. The expected output error for a bounded cylinder is at most its supremum norm times $`2\sqrt{2}e^t`$ times the square root of the failure bound, by Cauchy–Schwarz and the cell Yule second moments. A truncated output depends only on a finite marked ball with boundary degrees, so this proves continuity of $`Q_t f`$ for bounded cylinders. It also gives Borel kernels for the cell and for tracking a single root. Relabelling the cores at an intermediate time partitions the later cell (conditional expectation, without independence of cells), and hence $`Q_tQ_u=Q_{t+u}`$ and $`1\le Q_t1\le e^t`$.

**Existence and moments of the tilted law.** For ordinary growth from $`T_3`$ let $`m_t(f)=\mathbb E_{T_3}\sum_{x\text{ core}} f(G_t,x)`$; then $`m_t(1)=F_3(t)`$ and $`m_{t+h}=m_tQ_h`$. For $`0<\delta<\mu`$ and $`T>0`$ set $`w(t)=e^{-\mu t-\delta|t-T|}`$, $`Z=\int_0^\infty w(t)F_3(t)\,dt`$, and $`\pi_{\delta,T}=Z^{-1}\int_0^\infty w(t)m_t\,dt`$. The logarithmic growth limit implies $`Z<\infty`$. For fixed $`\delta`$, $`w(0)/Z\to0`$ as $`T\to\infty`$ since $`Z\ge e^{-\mu(T+1)-\delta}F_3(T)`$, using $`[T,T+1]`$. Choose a sequence $`\delta\downarrow0`$ and then $`T`$ large at each step so $`w(0)/Z\to0`$. For fixed $`h\ge0`$, the ratio $`w(u-h)/w(u)`$ for $`u\ge h`$ lies in $`[e^{(\mu-\delta)h},e^{(\mu+\delta)h}]`$. Changing variables in $`\pi_{\delta,T}Q_h`$ bounds its total variation difference from $`e^{\mu h}\pi_{\delta,T}`$ by an error tending to zero: the missing part on $`[0,h]`$ is $`O_h(w(0)/Z)`$ and the ratio error is uniform.

At each promotion the newborn has total degree two and core degree one; the initial $`T_3`$ centre also has total degree two and core degree zero. At age $`a`$, conditionally on the full environment at promotion, its expected degree is at most $`2e^a`$ by the bound on its gain rate, and its core degree is dominated by $`1+\mathop{\mathrm{Pois}}\nolimits (a)`$ since only its own promotions increase core degree. Thus the conditional moments $`D^\alpha`$ ($`0<\alpha<\min(1,\mu)`$) and $`e^{\theta c}`$ ($`0<\theta<\log(1+\mu)`$) are bounded by $`B e^{\kappa a}`$ with $`(B,\kappa)=(2^\alpha,\alpha)`$ and $`(e^\theta,e^\theta-1)`$ respectively. The total intensity of promotions is at most the core count; $`w(t+a)\le w(t)e^{-(\mu-\delta)a}`$. Compensating births and adding the initial centre therefore gives the tilted moment bound $`B(1+w(0)/Z)/(\mu-\delta-\kappa)`$ for $`\delta<\mu-\kappa`$.

For tightness of whole balls, the root degree is tight by the first moment estimate, and the two estimates with Cauchy–Schwarz give uniform convergence to zero of $`\pi_{\delta,T}(c_o\mathbf 1_{D_o>J})`$ as $`J\to\infty`$ along the chosen sequence (discarding finitely many terms). In each finite tree, exchanging the two endpoints in a sum of ordered paths preserves the root-uniform mean measure after summation. Consequently the expected number of paths of length $`i`$ from a root satisfying earlier degree cutoffs $`J_0,\ldots,J_{i-1}`$ and with final degree $`>J_i`$ is at most $`(J_1\cdots J_{i-1})\pi_{\delta,T}(c_o\mathbf1_{D_o>J_i})`$: in the reverse order there are at most $`c_o`$ choices at the bad endpoint, then at most the indicated cutoffs. Choose the cutoffs successively with summable failure probabilities. At every fixed radius the allowed marked balls form a finite set, and simultaneous cutoffs give a compact set in the local topology, so the laws are tight. Prokhorov's theorem supplies a weak subsequential probability limit $`\pi`$. For each fixed $`h`$, convergence against bounded cylinders and the continuity of their $`Q_h`$ images yield $`\pi Q_h=e^{\mu h}\pi`$. The same subsequence works for every $`h`$, since the total variation estimate held for each $`h`$ on the original sequence. Truncating the continuous root moment then proves $`\pi(e^{\theta c_o})<\infty`$.

**Identification of the degree marginal.** The residue (4) and the $`T_3`$ asymptotic in §4 give $`m_t(\mathbf1_{D_o\ge k})/F_3(t)\to d_{f_k}/d_{f_2}`$ for fixed $`k`$. This ratio equals $`c_k/c_2`$ by (6), where the edge-seed and fixed-size factors cancel; neither factor enters the $`T_3`$ local law. Under the tilt, bounded time intervals have mass tending to zero, since on each such interval $`w(t)\le w(0)e^{-(\mu-\delta)t}`$ and $`w(0)/Z\to0`$. Thus integrating the ratio against $`w(t)F_3(t)/Z`$ has the same limit. The degree indicator is a bounded continuous cylinder, so $`\pi(D_o\ge k)=c_k/c_2`$ and $`\pi(D_o=k)=p_k`$. No uniqueness of the full environment law is needed.

## 6. Exponential age, positive tail constant, and crossing flux

Let $`S_t f`$ track the same old root under physical evolution (one root, rather than a cell sum), and $`\mathcal B f(G,o)=b_o f(G^{b_o},v)`$ where $`v`$ is the promoted leaf and $`G^{b_o}`$ is the state just after promotion. In core-tree coordinates a promotion always creates the same marked attachment, irrespective of which indistinguishable leaf is chosen. The persistent old root contributes $`S_t f`$; each other member of the cell has one promotion time. Conditional expectation for its subsequent trajectory and compensation at the previous promoter give, for bounded nonnegative $`f`$,

```math
Q_t f=S_t f+\int_0^t Q_u\mathcal B S_{t-u} f\,du.
```

Indeed the conditional total birth rate is bounded by the Yule-dominated cell size, so the integral is finite. The stopped finite construction gives the compensation formula and the local construction gives the same predictable rates after removing the stop; alternatively verify on cylinders and extend by a monotone class. Integrate against $`\pi`$, divide by $`e^{\mu t}`$, and change variable $`a=t-u`$. The old-root term vanishes as $`t\to\infty`$ since $`\mu>0`$, yielding

```math
\pi f=\int_0^\infty e^{-\mu a}\pi\mathcal B S_a f\,da.\qquad\text{(7)}
```

At $`f=1`$ this gives $`\pi\mathcal B1=\mu`$, so $`\beta=\pi\mathcal B/\mu`$ is a probability of a newborn and its full environment. Formula (7) says to start from $`\beta`$ and observe the tracked root at an independently sampled $`\mathop{\mathrm{Exp}}\nolimits (\mu)`$ age. Its seed is a newborn of degree two with one leaf and its real promoter, **not** the isolated $`T_2`$ edge, the $`T_3`$ centre, or a planted particle with an inert slot.

Starting from $`\beta`$, the tracked degree $`D(t)`$ grows by single steps at rate at most $`D(t)`$. Its exponentially normalized value $`e^{-t}D(t)`$ is a nonnegative supermartingale, by conditional rate integration, and converges almost surely to $`W_\infty`$. Couple the degree below a Yule process $`Y`$ from two by supplying the difference of rates with extra jumps; this works for every exterior since $`D(t)\le Y(t)`$ is preserved at each jump. The squared mean of the Yule martingale $`e^{-t}Y(t)`$ is $`6-2e^{-t}`$, so Doob's $`L^2`$ maximal inequality gives $`\mathbb E[\sup_t(e^{-t}Y(t))^2]\le24`$. Hence $`A=\mathbb E_\beta W_\infty^\mu`$ is finite, and the normalized degrees are uniformly integrable. These bounds use only the newborn degree, not the exterior. The same estimate from any given finite starting degree gives uniform integrability with constant depending on that degree.

For strict positivity, silence the newborn's promotion clock for one unit and let its tagged initial leaf ring at least twice. This leaf remains a leaf throughout because only the root can promote it. One may implement the rate-$`L`$ leaf stream in the local construction by the individual rate-one leaf clocks, so this event has probability at least $`e^{-1}(1-2e^{-1})>0`$ uniformly over the exterior. Then $`L(1)\ge3`$. From that state leaf drift is $`L+\rho_o-b_o\ge L-1`$, where $`\rho_o=\sum_{w\sim o,\ w\text{ core}}1/D_w`$. Stopped compensation and the degree moment bound imply $`\mathbb E[L(1+t)-1\mid\mathcal F_1]\ge e^t(L(1)-1)`$; use the integrating factor on $`L-1`$. Consequently, by conditional uniform integrability (the degree at time one is finite), $`\mathbb E[e W_\infty\mid\mathcal F_1]\ge L(1)-1>0`$ on this event. Thus $`A>0`$.

With $`T_k=\inf\{t:D(t)\ge k\}`$, (7) gives $`\pi(D_o\ge k)=\mathbb E_\beta[e^{-\mu T_k};T_k<\infty]`$. When $`D`$ grows without bound, $`D(T_k)=k`$ and $`T_k\to\infty`$ by nonexplosion, so $`k e^{-T_k}\to W_\infty`$. For bounded degree, the stopped contribution is eventually zero and $`W_\infty=0`$. In either case the $`\mu`$th power is bounded by the integrable maximal variable just constructed, so dominated convergence yields $`k^\mu\pi(D_o\ge k)\to A`$.

**Pointwise flux.** Subtracting two cumulative asymptotics alone would not justify the pointwise tail. For $`k\ge2`$, a newborn of degree two cannot already lie above level $`k`$. The cell count of vertices with degree at least $`k+1`$ increases only on unit degree gains from exactly $`k`$. The instantaneous intensity per member is

```math
a_k(G,o)=\big(k-c_o+\rho_o\big)\mathbf 1_{D_o=k}.
```

It is a bounded continuous cylinder, at most $`k`$. With $`F_k=\mathbf1_{D_o\ge k+1}`$, compensation gives $`Q_tF_k-F_k=\int_0^tQ_u a_k\,du`$; the Yule first moment bounds the integral, even for an infinite exterior. Apply $`\pi Q_t=e^{\mu t}\pi`$ and divide by $`(e^{\mu t}-1)/\mu`$ at any $`t>0`$ to get

```math
\mu\pi(D_o\ge k+1)=\pi(a_k).\qquad\text{(8)}
```

All core neighbours have degree at least two, so $`0\le\rho_o\le c_o/2`$ and $`a_k\le k\mathbf1_{D_o=k}`$. Conversely for fixed $`0<\varepsilon<1`$, $`a_k\ge (1-\varepsilon)k\mathbf1_{D_o=k,\ c_o\le\varepsilon k}`$. With the exponential moment, $`\pi(c_o>\varepsilon k)\le\pi(e^{\theta c_o})e^{-\theta\varepsilon k}`$. Thus (8) bounds $`k^{1+\mu}\pi(D_o=k)`$ from below in the limit by $`\mu A`$, and from above by $`\mu A/(1-\varepsilon)`$. Let $`k\to\infty`$ first, then $`\varepsilon\downarrow0`$. This proves the claimed pointwise tail, with strictly positive constant $`\mu A`$. In particular $`A`$ is determined by the degree marginal, even if the extraction of $`\pi`$ is not taken to give a unique full environment.

## References and audit scope

- **Compact positive operators:** the strong Krein–Rutman result used in §3 is the case of a compact strongly positive operator on a Banach space with a solid generating cone: the spectral radius is positive and algebraically simple, with an interior eigenvector and positive dual eigenfunctional. See Lei Zhang, [*A generalized Krein–Rutman theorem*, Theorem 1.1](https://doi.org/10.3934/dcds.2020054), taking essential spectral radius zero. Here complexification is used for holomorphy; the ordered space is real $`C(\mathcal A)`$. The eigenfunctional is a finite positive measure by the Riesz representation theorem on compact Hausdorff spaces (Rudin, *Real and Complex Analysis*, Theorem 2.14). Arzelà–Ascoli is used in its compact-domain form for uniformly bounded equicontinuous families (Rudin, *Principles of Mathematical Analysis*, Theorem 7.25). The compact-operator Fredholm alternative used here is the Riesz–Schauder assertion that every nonzero spectral value of a compact operator is an eigenvalue (Conway, *A Course in Functional Analysis*, 2nd ed., Chapter II, §4, compact operators); the local inverse and contour argument are given in §3.
- **Tauberian step:** G. Debruyne and J. Vindas, [*Complex Tauberian theorems for Laplace transforms with local pseudofunction boundary behavior*, Theorem 1.1 and the paragraph following](https://arxiv.org/abs/1611.09765): for locally integrable, log-linearly slowly decreasing $`S`$ with transform convergent on $`\Re z>1`$, subtracting $`d/(z-1)`$ and having local pseudofunction boundary behavior on $`\Re z=1`$ yields $`S(x)=e^x(d+o(1))`$. The analytic-boundary sufficient condition and its application are in §4. The fixed-size transfer there derives interval convergence by polynomial approximation on $`[0,1]`$ (Weierstrass, Rudin, *Principles of Mathematical Analysis*, Theorem 7.26), rather than assuming an automatic transfer for an arbitrary sequence.
- **Probability tools:** the Poisson fresh-increment and compensation rules used here apply to predictable intensities, with localization and the displayed integrability bounds; see Brémaud, *Point Processes and Queues: Martingale Dynamics*, Chapter II, “Intensity and Compensators.” Tightness implies relative weak compactness on a Polish probability space (Billingsley, *Convergence of Probability Measures*, 2nd ed., §5, Prokhorov's theorem); weak convergence tests bounded continuous functions and bounded tests with null discontinuity set (ibid., §2, Portmanteau). The nonnegative supermartingale convergence theorem and Doob's $`L^2`$ maximal inequality apply to right-continuous adapted processes, the latter to a square-integrable martingale on each finite horizon and then by monotone convergence (Durrett, *Probability: Theory and Examples*, 5th ed., Chapter 4, §§4.2 and 4.4). In §6 the dominating martingale has the requisite uniform second-moment bound explicitly.
- **Original audited renewal (`bfcbbe5`):** Neon referee report, [independent code](../checks/fringe_neon.py), and recorded output. The report audits the positivity and mark transfer, weight and immigration inequalities, analytic compactness and critical threshold, plus the separate constant-factor claim. Its finite checks include formal generator agreement on 1,039 reachable marked trees, the rational weight inequality, and nine planted errors. The [author's arithmetic replay](../fringe/check.py) and output are separate from the independent checks.
- **Original audited asymptotics (`8e4db79`):** Carbon referee report, [checker](../checks/asymptotics_carbon.py), and output cover the transform, pole, rescaling, fixed-size sandwich and degree-law bridge, with exact finite count/mark checks and eight planted errors (including the wrong seed wait, Yule weight, gamma factor, and a lattice reproduction kernel). Cell operator review and degree-law review, with [checker](../checks/operator_carbon.py) and output, audit local cutoff, cell semigroup, birth-time mixture and crossing flux; negative controls distinguish tracking the same root from using its full cell, and detect omission of $`-c_o`$ in (8). The original commit's operator proof is needed here only through §§1–2; the current operator review also covers its later §5, which this proof does not invoke.

For the standalone draft, the author replayed `round2/checks/fringe_neon.py` and `round2/checks/asymptotics_carbon.py` on CPU 13 with `timeout 12s python3 -B`; both passed, including nine and eight planted controls respectively. The independent Sol review of this standalone proof found it correct as stated, taking the cited growth theorem and motif bracket as inputs. Sol ran a fresh standalone fixture on CPU 13 and completed the full `round2/checks/operator_carbon.py` replay on CPU 8: all four positive check groups passed, including the 79-tree edge-transfer check, and all eight planted defects were caught. The finite checks do not substitute for the analytic extensions, infinite-volume limits, or flux arguments above.
