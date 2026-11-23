> Produced by GPT-6.2 Sol (under Gemini 4 Carbon).

# An almost-sure lower bound from local drift certificates

Start the random friend tree with an edge. In the step from $`n`$ to $`n+1`$, choose a uniform vertex $`V_n`$, then a uniform neighbour $`W_n`$ of $`V_n`$, and attach a new leaf to $`W_n`$. Let $`\mathcal F_n`$ record the history through size $`n`$. Write $`K_n`$ for the core (vertices of degree at least two), $`C_n=|K_n|`$, and root the tree at the centre $`r`$ of $`T_3`$. For a core vertex $`u`$, write $`D_u=L_u+c_u`$ for its degree, with $`L_u`$ leaf neighbours and $`c_u`$ core neighbours.

## Theorem and certificate hypotheses

Let $`f_-`$ and $`f_+`$ be radius-$`1`$ local per-vertex tables of the kind in `REPORT.md` §§3.3–3.5: the state $`x_u`$ of a core vertex $`u`$ is a function of its own $`(L_u,c_u)`$, the current $`(L_v,c_v,D_v)`$ of its core neighbours $`v\sim u`$ (optionally distinguishing its parent in $`T_n`$ from its children), and optionally a static birth label determined by its creator's $`(L,c,D)`$ at the step $`u`$ was promoted; moreover, there is a finite cutoff $`H\ge 3`$ such that $`f_\pm(x_u)`$ does not change when $`u`$ gains a leaf at $`L_u\ge H`$, and any neighbour (or parent/creator) of degree at least $`H`$ has a permanently saturated descriptor in $`x_u`$. The two tables may use different state encodings and caps. Assume $`0\le f_-\le F`$, $`1\le f_+\le U<\infty`$, and $`a_0=f_-(x_r(T_3))>0`$, the initial positivity used in the report's mean lower bound. Assume each local drift certificate holds at every true choice of rates, with exponents $`0<\delta\le\lambda<2\delta`$. Then

```math
\boxed{\displaystyle\liminf_{n\to\infty}\frac{C_n}{n^\delta}>0\quad\text{almost surely}.}
```

In fact the whole-tree lower potential divided by $`n^\delta`$ has a limit in $`(0,\infty]`$ almost surely. Consequently $`\liminf_n(\log C_n)/\log n\ge\delta`$. (For radius-$`R`$ local tables in which a vertex of degree at least $`H`$ also blocks dependence on vertices beyond it, the same proof applies upon replacing pattern depth $`2`$ in §4 by depth $`R+1`$, since the expected number of core vertices at depth at most $`R`$ is $`O((\log n)^R)=o(n^\delta)`$.)

Here is the precise meaning of a local drift certificate in the proof. If $`x_{u,n}`$ is the state evaluated in the **full tree**, and $`y`$ denotes the old leaf promoted when the receiver is a leaf of $`u`$, the certified quantity for $`u`$ is

```math
n\,\mathbb E\!\left[f(x_{u,n+1})-f(x_{u,n})
       +\mathbf 1_{\{W_n\text{ is a leaf of }u\}} f(x_{y,n+1})\ \middle|\ \mathcal F_n\right]. \qquad (1)
```

The last term is zero in the absence of a promotion at $`u`$. The lower certificate bounds (1) below by $`\delta f_-(x_{u,n})`$, and the upper bounds it above by $`\lambda f_+(x_{u,n})`$. In the report, (1) is evaluated by listing events that change the state of $`u`$ and crediting the entire newborn contribution to its creator. Robust bounds on the event rates give (1) at the actual rates. We use the original creator labels and root throughout. In a parent-distinguished radius-$`1`$ variant, the child counts only track children below a finite degree cap and the parent's $`(L_p,c_p,D_p)`$ descriptor is permanently `P_top` above its own cap. Any coarsening of such a radius-$`1`$ state is allowed provided the table stays bounded and satisfies the degree-$`H`$ saturation condition above. Initial positivity and (1) are the hypotheses supplied by a certificate.

## 1. Receivers, jumps, and the conditional second moment

Every leaf has a core neighbour from size three onward. Given $`\mathcal F_n`$, a specified vertex $`w`$ is the receiver with probability $`n^{-1}\sum_{z\sim w}D_z^{-1}`$. Hence for a core vertex $`u`$ the rates (probabilities multiplied by $`n`$) of receiving a leaf and of promoting one of its leaves are, respectively,

```math
L_u+\rho_u,\qquad L_u/D_u,\qquad
\rho_u=\sum_{v\sim u,\,v\in K_n}1/D_v.
```

A leaf gain changes $`(L_u,c_u)`$ to $`(L_u+1,c_u)`$. A promotion changes it to $`(L_u-1,c_u+1)`$ at unchanged degree and adds one core vertex with $`(L,c)=(1,1)`$. Degrees never decrease. The probability of any promotion is $`S_n/n`$, where $`S_n=\sum_{u\in K_n}L_u/D_u\le C_n`$. The credit of a single vertex need not be bounded by one. The bound used here is

```math
\sum_{u\in K_n}\rho_u=\sum_{v\in K_n}c_v/D_v\le C_n. \qquad (2)
```

For the `certs9.py` lower table a leaf gain at $`u`$ can change its own entry and the entries at some of its old core neighbours. A class-$`\ell`$ pendant has only one such neighbour. Class B has two. A tracked “other” of degree $`d<D_{\rm cap}`$ has at most $`D_{\rm cap}-1`$ such neighbours. An “other” in `top` changes none of its neighbours' entries. An “other” stays “other” under gains and promotions, since its degree never falls and a promotion increases its core-neighbour count. A promotion at $`u`$ changes $`u`$'s entry, introduces one new entry, and can change the entry of at most one old core neighbour, namely when $`u`$ was a pendant. It could not have been class B, and a tracked “other” retains its degree. The absolute change of each summand lies in $`[0,F]`$, including the introduction of a summand. Thus the lower jump is bounded by $`B=F\max(3,D_{\rm cap})`$.

If a core receiver $`u`$ has $`L_u\ge L_{\max}`$, its own leaf coordinate is already flat. It is neither a tracked pendant nor class B, and it is in `top` whenever it has a core neighbour: then $`D_u\ge L_{\max}+1\ge D_{\rm cap}`$. Here $`L_{\max}>K_0`$ and $`L_{\max}+1\ge D_{\rm cap}`$ are the verifier's conditions. If $`c_u=0`$, there is no neighbour entry to change. This checks the claim also in the case $`L_{\max}=D_{\rm cap}-1`$. Consequently

```math
\begin{split}
\mathbb P(W_n\in K_n,\ L_{W_n}<L_{\max}\mid\mathcal F_n)
&\le (L_{\max}+1) C_n/n,\\
\mathbb E[(\Phi_{n+1}^- -\Phi_n^-)^2\mid\mathcal F_n]
&\le B^2(L_{\max}+2)C_n/n, \qquad (3)
\end{split}
```

where $`\Phi_n^\pm=\sum_{u\in K_n}f_\pm(x_{u,n})`$ uses the state for its own table. The first inequality follows by summing the receiver rates and applying (2). The second adds the promotion probability.

For both tables and all the variants in the theorem, choose one integer $`H\ge3`$ at least as large as both own-leaf cutoffs, all neighbour or creator degree caps, and $`K_0+2`$ for each pendant classification in use (for the two final tables of `REPORT.md` §3.5, $`H=60`$ suffices). A core vertex of degree at least $`H`$ is an “other” in any pendant/B classification. Its current descriptor in any adjacent state is permanently saturated, and any static creator label of an already existing neighbour is fixed. Thus a gain at a core vertex can affect old neighbour entries only if its degree before the step is below $`H`$. A promotion can affect old neighbour entries only if the creator's degree is below $`H`$; the degree of the creator is unchanged by the promotion. In particular even when the rooted table records the parent's $`L,c`$ separately, at most $`H-1`$ old neighbours can change in either step. Including the receiver or creator itself and a possible promoted leaf, $`H+1`$ bounds the number of changed summands. For a table with entries in $`[0,J]`$ one may use $`B_J=(H+1)J`$. A core receiver with $`L\ge H`$ changes no entries at all. These facts give the whole-tree bound (3) with $`B_J`$ and $`H+2`$ in place of $`B`$ and $`L_{\max}+2`$.

## 2. Restriction across a saturated parent edge

Suppose $`v`$ is promoted at a vertex $`p`$ whose degree is at least $`H`$ at the promotion. From that time on let $`R_n`$ consist of all vertices on the $`v`$-side of the edge $`\{p,v\}`$, and put $`I_n=K_n\cap R_n`$, $`N_n=|I_n|`$. Use $`I_n=K_n`$, $`N_n=C_n`$ for the whole-tree case. Define $`X_n=\sum_{u\in I_n} f_-(x_{u,n})`$ and $`Y_n=\sum_{u\in I_n} f_+(x_{u,n})`$. A promotion belongs to $`I_{n+1}\setminus I_n`$ exactly when its creator belonged to $`I_n`$. Therefore summing (1) over $`I_n`$ gives the exact drift of the restricted sum, and

```math
\mathbb E[X_{n+1}\mid\mathcal F_n]\ge(1+\delta/n)X_n,\qquad
\mathbb E[Y_{n+1}\mid\mathcal F_n]\le(1+\lambda/n)Y_n,\qquad
0\le X_n\le F N_n,\quad N_n\le Y_n\le U N_n. \qquad (4)
```

These are full-tree states, so $`v`$ retains its parent descriptor, and $`1/D_p`$ still contributes to $`v`$'s true receiver rate. Formula (4) allows that rate to depend on the exterior history. The boundary degree is used to control jumps: an exterior receiver changes no state inside the branch. The only exterior vertex that can be seen by an inside vertex is $`p`$, and only $`v`$ sees it. Inside $`v`$'s state this neighbour is permanently `top` or `P_top` where tracked. If the table uses a remembered creator, $`p`$ was already “other” when $`v`$ was born, so its label is permanently at the cap (or the static no-claim value when that convention applies). Exterior changes do not change the membership or the $`L,c`$ of any inside core vertex. These arguments also cover an exterior receiver that is a leaf of $`p`$.

For a restricted jump, the probability of an inside promotion is at most $`N_n/n`$. For gains of inside core vertices with $`L<H`$, the sum of the credits counted over $`I_n`$ is

```math
\sum_{u\in I_n}\rho_u
=\sum_{z\in I_n}\frac{|\{u\in I_n:u\sim z\}|}{D_z}+\frac1{D_p}\le N_n+1\le 2N_n. \qquad (5)
```

The boundary term is omitted in the whole-tree case, and $`N_n\ge1`$ in either case. Using the same jump sizes as above, (5) bounds the second moments of the increments simultaneously:

```math
\mathbb E[(\Delta X_n)^2\mid\mathcal F_n]\le K_- N_n/n,\qquad
\mathbb E[(\Delta Y_n)^2\mid\mathcal F_n]\le K_+ N_n/n,\quad
K_-=(H+3)B_F^2,\quad K_+=(H+3)B_U^2. \qquad (6)
```

## 3. The normalized martingale and a uniform survival probability

Put $`P_a(m,n)=\prod_{j=m}^{n-1}(1+a/j)`$, with $`P_a(m,m)=1`$. Uniformly in $`n\ge m\ge3`$ and for fixed $`a>0`$,

```math
e^{-a^2}(n/m)^a\le P_a(m,n)\le e^a(n/m)^a. \qquad (7)
```

Indeed $`x-x^2/2\le\log(1+x)\le x`$ for $`x\ge0`$ by integration. Comparing the harmonic sum with the integral gives $`0\le\sum_{j=m}^{n-1}1/j-\log(n/m)\le1/m`$ for $`n>m`$. Also $`\sum_{j=m}^\infty j^{-2}\le1/(m-1)`$. For fixed $`m`$, $`P_a(m,n)/n^a`$ has a positive finite limit: the absolute discrepancies between $`\log(1+a/j)`$ and $`a\log(1+1/j)`$ are bounded by $`(a^2+a)/(2j^2)`$, and $`\sum_{j=m}^{n-1}\log(1+1/j)=\log(n/m)`$.

Condition on the history at a starting size $`m`$. This may be an almost surely finite stopping size on the histories of interest, with the branch already specified. The following calculations can be done at each fixed $`m`$ and each history atom through that size, so they also apply conditionally at such a stopping size. From (4), iteration gives

```math
\mathbb E[N_j\mid\mathcal F_m]\le Y_m P_\lambda(m,j)\le U N_m P_\lambda(m,j). \qquad (8)
```

The normalized lower sum $`M_n^{(m)}=X_n/P_\delta(m,n)`$ satisfies

```math
M_n^{(m)}=X_m+A_n+Z_n,\qquad
Z_{j+1}-Z_j=\frac{X_{j+1}-\mathbb E[X_{j+1}\mid\mathcal F_j]}{P_\delta(m,j+1)}. \qquad (9)
```

Here $`A_m=Z_m=0`$, and each increment of $`A`$ is the nonnegative conditional mean increment of $`M^{(m)}`$. The numerator in (9) is also the centred increment of $`X`$. Its conditional variance is at most the uncentred conditional second moment in (6). Let $`q=2\delta-\lambda>0`$. By (6)–(8), the sum of martingale increment second moments, conditional on the starting history, is at most

```math
\begin{split}
K_- U N_m \sum_{j=m}^\infty\frac{P_\lambda(m,j)}{jP_\delta(m,j+1)^2}
&\le K_- U e^{\lambda+2\delta^2} N_m \sum_{j=m}^\infty \frac{(j/m)^{-q}}j\\
&\le D_- N_m,\qquad D_-=K_- U e^{\lambda+2\delta^2}(1+1/q). \qquad (10)
\end{split}
```

The last sum is bounded by its first term $`1/m`$ plus $`\int_m^\infty (t/m)^{-q}\,dt/t=1/q`$. This also proves $`\sum_{j\ge3} C_j/j^{1+2\delta}<\infty`$ almost surely by (8), (7), and finiteness of the expectation of the nonnegative sum.

For completeness, orthogonality of the centred increments shows that $`Z_n`$ is Cauchy in $`L^2`$. To see almost-sure convergence, for $`s<t`$ let the first index with $`|Z_k-Z_s|\ge b`$ be the crossing index, if any through $`t`$. On each event that this first index equals $`k`$, conditional Jensen gives $`\mathbb E[(Z_t-Z_s)^2\mid\mathcal F_k]\ge (Z_k-Z_s)^2`$. Summing gives probability of a crossing at most $`b^{-2}\mathbb E[(Z_t-Z_s)^2]`$. Take $`t\to\infty`$ and starting indices $`s_i`$ with remaining incremental second moments at most $`2^{-4i}`$. With $`b=2^{-i}`$, the probabilities are summable, so almost surely only finitely many crossings occur for these pairs $`(s_i,2^{-i})`$. The sequence is then almost surely Cauchy. Its limit $`Z_\infty`$ has conditional mean zero and conditional second moment at most $`D_-N_m`$. Since $`A_n`$ is nondecreasing, $`M_n^{(m)}`$ converges almost surely in $`[0,\infty]`$.

For the whole tree this gives $`\mathbb E[M_n^{(3)}]\ge a_0`$. For a second-moment bound, apply the same decomposition to $`Y_n/P_\lambda(3,n)=Y_3-A_n^++Z_n^+`$, with $`A_n^+`$ nondecreasing. Its martingale has bounded second moments by (10) with $`K_+`$ in place of $`K_-`$ and both exponents equal to $`\lambda`$. Thus $`0\le Y_n/P_\lambda(3,n)\le Y_3+Z_n^+`$ and $`\sup_n\mathbb E[(Y_n/P_\lambda(3,n))^2]<\infty`$. In particular,

```math
\mathbb E[(M_n^{(3)})^2]\le F^2\big(P_\lambda(3,n)/P_\delta(3,n)\big)^2\,\sup_j\mathbb E[(Y_j/P_\lambda(3,j))^2]
=O(n^{2(\lambda-\delta)}). \qquad (11)
```

The square-summability in (10) controls the martingale part. The lower predictable part can grow, and (11) gives a uniform bound when $`\delta=\lambda`$. In general the normalized lower sum can have unbounded second moments. For example, the report's lower table also certifies $`\delta=0.542`$. With that normalization its mean is bounded below by a positive constant times $`n^{0.001}`$, by the same table's $`0.543`$ drift and (7). Jensen's inequality then gives a divergent lower bound for its second moment.

If $`X_m=a>0`$, set $`J=a+Z_\infty`$. By Cauchy–Schwarz, conditionally on the starting history,

```math
a=\mathbb E[J\mid\mathcal F_m]\le a/2+
\big((a^2+D_-N_m)\,\mathbb P(J>a/2\mid\mathcal F_m)\big)^{1/2}.
```

On $`\{J>a/2\}`$, (9) gives a positive normalized limit, and (7) gives $`\liminf X_n/n^\delta>0`$. Thus as long as $`X_m\ge a_*>0`$ and $`N_m\le k_*<\infty`$, the conditional probability of this last event is at least

```math
p_* = \frac{a_*^2}{4(a_*^2+D_-k_*)}>0. \qquad (12)
```

The whole-tree calculation started at arbitrary $`m`$ instead gives the bound $`X_m^2/[4(X_m^2+D_-C_m)]`$. Initial positivity and bounded summands alone do not keep this bound uniformly positive at late times. A branch of bounded starting size supplies the needed uniform restart.

## 4. Reaching a positive state in such a branch

### A finite pattern with positive weight

By (4) for the whole tree and (7), $`\mathbb E[\Phi_n^-]\ge a_0 P_\delta(3,n)`$ grows at least as a positive constant times $`n^\delta`$. Every depth-one core vertex was promoted at $`r`$. At each size that promotion has probability $`L_r/(nD_r)\le1/n`$. Thus the expected number of core vertices of depth at most one is at most $`1+\sum_{j=3}^{n-1}1/j`$. If every state of positive weight occurring in a reachable history belonged to a vertex of depth at most one, $`\mathbb E[\Phi_n^-]`$ would be at most $`F(1+\sum_{j=3}^{n-1}1/j)`$ for all $`n`$. Since $`\delta>0`$, some finite reachable history to a size $`k_*`$ has a vertex of depth at least two with positive weight $`a_*`$. Fix that history, described by its successive receivers from $`T_3`$.

At a qualifying promotion the branch consists of its founder $`v`$ and one new leaf. First attach another leaf to $`v`$, then replay the fixed receivers inside the branch, with $`v`$ corresponding to $`r`$. The extra edge to $`p`$ changes the full $`L,c,D`$ data relative to the pattern only at the copy of $`r`$ (it has the same $`L`$, one extra core neighbour, and degree one higher). A vertex of pattern depth at least two sees neither that extra edge nor the copy of $`r`$ in its own or adjacent $`L,c,D`$ data. The creator's data at its birth agree too, since its creator has depth at least one. Thus the corresponding positive state has weight $`a_*`$, including when a creator label or a rooted parent slot is used. At the end the branch has $`k_*`$ vertices in total and $`X_m\ge a_*`$. This argument does not require that the newborn with a saturated parent itself have positive weight.

The same reasoning applies if the lower potential first becomes positive later. It suffices that $`\mathbb E[\Phi_s^-]>0`$ at some fixed $`s`$, since (4) then gives the mean lower bound from size $`s`$.

### Uniform probability of creating and replaying the branch

Use the following elementary rule for waiting for a step in a specified class. If until the first such step its conditional probability is at least $`b/n`$, for some $`b>0`$ fixed by the starting history, the step occurs almost surely. Indeed the probability of no occurrence through a finite range is bounded by the product of the factors $`1-b/n`$, which tends to zero by divergence of the harmonic series. If on those same histories a desired subclass has conditional probability at least a fraction $`h`$ of the total probability of the class, the first step is in the subclass with probability at least $`h`$: sum the one-step inequalities with the indicator of no previous occurrence. This reasoning also applies at finite stopping sizes by conditioning on their histories.

At any size $`m`$ choose a leaf and let $`u`$ be its core neighbour. While $`D_u<H`$ and $`L_u\ge1`$, wait for a gain at $`u`$ or a promotion at $`u`$. The total rate is at least one, and the gain rate $`L_u+\rho_u`$ is at least the promotion rate $`L_u/D_u`$. The next such step is a gain with conditional probability at least $`1/2`$. A run of at most $`H`$ such gains reaches $`D_u\ge H`$ with a leaf still present, with conditional probability at least $`2^{-H}`$. If $`u`$ already has degree at least $`H`$, this stage is empty. Once $`D_u\ge H`$ and $`u`$ has a leaf, wait for the first promotion at $`u`$. Until then $`L_u/D_u`$ is nondecreasing, because the only possible changes at $`u`$ are gains at fixed $`c_u`$. The promotion rate is therefore bounded below by its positive value at the start of this wait, so a qualifying promotion occurs almost surely on reaching this stage.

The same waiting argument verifies sibling spawning. If a vertex of degree at least $`H`$ starts with $`R`$ leaves, it performs at least $`R`$ promotions almost surely. After fewer than $`R`$ promotions it still has a leaf, as each promotion uses one and gains only add leaves. The promoted vertices give distinct children with disjoint descendant sides. In particular if $`C_m\le m^{0.6}`$, some vertex has at least $`(m-C_m)/C_m\ge m^{0.4}-1`$ leaves by pigeonhole, so it can supply arbitrarily many such children as $`m`$ grows through these sizes.

After the promotion, when the branch has $`k\ge2`$ vertices, any particular receiver inside has an inside neighbour of degree at most $`k`$. Its conditional rate is at least $`1/k`$. The total rate of receivers inside is at most $`k+1`$, because only the $`k`$ inside vertices and the single exterior parent can be the uniform first vertex $`V_n`$ for such a receiver. The inside tree and its size stay fixed until the next inside receiver. The waiting rule shows that its identity matches any prescribed target with conditional probability at least $`1/[k(k+1)]`$. Multiplying for the one extra leaf and the finite replay gives conditional probability at least

```math
p_{\mathrm{pattern}}=\prod_{k=2}^{k_*-1}\frac1{k(k+1)}>0 \qquad (13)
```

of reproducing the pattern at a finite stopping size. All bounds are uniform over the exterior history, including changes in $`D_p`$.

## 5. From a uniform chance to probability one

Let $`E=\{\liminf_n\Phi_n^-/n^\delta>0\}`$. For every deterministic $`m`$, the procedure in §4 can be started at $`m`$. Condition successively at the qualifying promotion and the end of the replay, then apply (12). Throughout the resulting branch $`X_n\le\Phi_n^-`$ by nonnegativity. We obtain, almost surely for every $`m\ge3`$,

```math
\mathbb P(E\mid\mathcal F_m)\ge 2^{-H}p_{\mathrm{pattern}}p_*>0. \qquad (14)
```

Here is a direct justification that (14) implies $`\mathbb P(E)=1`$. Put $`p_0=2^{-H}p_{\mathrm{pattern}}p_*`$. For every event $`A`$ in the algebra $`\bigcup_m\mathcal F_m`$, (14) gives $`\mathbb P(E\cap A)\ge p_0\mathbb P(A)`$. Every event in $`\sigma(\bigcup_m\mathcal F_m)`$ is approximable by this algebra in probability of symmetric difference. Indeed the class of approximable events contains the algebra, is closed under complements, and is closed under countable unions by first approximating a union by a finite subunion in probability and then approximating those finitely many events. Since $`E^c`$ lies in the generated sigma-field, approximation gives $`0\ge p_0\mathbb P(E^c)`$.

This uses one branch per conditional trial, so it applies even when trials share vertices or depend on the same history. It also shows there are almost surely infinitely many qualifying promotions somewhere in the tree. For each fixed cutoff size, the conditional probability at every later size of a qualifying promotion after the cutoff is at least $`2^{-H}`$ by §4, so the same algebra argument applies.

Finally, $`C_n\ge\Phi_n^-/F`$ and $`E`$ has probability one. Equation (9) at $`m=3`$, together with the positive finite product limit in §3, proves convergence of $`\Phi_n^-/n^\delta`$ in $`(0,\infty]`$. With the certificates at $`\delta=543/1000`$, $`\lambda=293/500`$, the theorem gives $`\liminf_n C_n/n^{0.543}>0`$ almost surely.
