# Astra Claim 6 — Almost-sure core size and degree frequencies

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 6: Claude and Gemini; GPT standalone audit; inputs include Round 1 almost-sure exponent and Claims 2, 3 and 5. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Start from a single edge. At size $`n`$, choose a vertex uniformly, choose one of its neighbours uniformly, and attach a new leaf to the neighbour. Write $`T_n`$ for the resulting chain, $`C_n`$ for the number of vertices of degree at least two (the **core**), and $`N_{\ge k}(n)`$ and $`N_{=k}(n)`$ for cumulative and exact degree counts. Let $`\mu`$ be the mean core growth exponent. Then there is a finite, square-integrable random variable $`X>0`$ almost surely such that, simultaneously for all integers $`k\ge2`$,

```math
\frac{C_n}{n^\mu}\longrightarrow X,\qquad
\frac{N_{\ge k}(n)}{C_n}\longrightarrow\frac{c_k}{c_2},\qquad
\frac{N_{=k}(n)}{C_n}\longrightarrow p_k=\frac{c_k-c_{k+1}}{c_2}
\quad\text{almost surely}.
```

Here $`\mathbb E N_{\ge k}(n)\sim c_k n^\mu`$, $`c_k>c_{k+1}>0`$, $`c_k\downarrow0`$, and $`\mathbb E X=c_2`$. The $`p_k`$ form a probability distribution; the previously established degree-law estimate gives $`p_k\sim\mu A k^{-1-\mu}`$ with $`0<A<\infty`$. For every $`\varepsilon,M>0`$, $`\mathbb P(X<\varepsilon)>0`$ and $`\mathbb P(X>M)>0`$. In particular, $`X`$ is not constant.

## 1. Inputs and time parametrization

The two earlier exponent inputs, for exactly the single-edge chain just defined, are:

1. [Mean exponent, Theorem 1 of `proofs/existence.md`](../../proofs/existence.md): $`\log \mathbb E C_n/\log n\to\mu`$. Its continuous-time proof (Sections 6–8) also gives $`\log \mathbb E_{T_2} C(t)/t\to\mu`$. [Finite-size product bounds, Theorem 1 of `proofs/bracket-motifs.md`](../../proofs/bracket-motifs.md) put the mean exponent in $`[113/200,1133/2000]`$. We need only $`51/100<\mu<1`$, in particular $`2\mu>1`$.
2. [Logarithmic sample-path exponent, Theorem 2 of `proofs/almost-sure-exponent.md`](../../proofs/almost-sure-exponent.md): $`\log C_n/\log n\to\mu`$ almost surely. This is the earlier logarithmic result, proved by blocks and finite-branch comparisons, not convergence at the normalization $`n^\mu`$. We use it only to prove strict positivity of the limiting martingale.

The spatial eigenmeasure inputs used **only for the large-$`k`$ statement** are specified in Section 7. All other model-specific steps are proved below.

In continuous time give each vertex a rate-one clock. At a ring it selects a uniform neighbour. Thus the arrow from $`x`$ to an actual neighbour $`y`$ has intensity $`1/D_x`$, where $`D_x`$ is the degree of the **source** $`x`$. Every ring adds one vertex. When the size is $`m`$, the next holding time is exponential of rate $`m`$, independent of the uniform source and its subsequent choice. Since $`m`$ grows deterministically by one at each jump, construct the whole continuous chain using a discrete jump chain and an independent sequence $`E_m\sim\mathop{\mathrm{Exp}}\nolimits (m)`$. From a seed of size $`m_0`$, the epoch at size $`n`$ is $`\tau_n=\sum_{m=m_0}^{n-1}E_m`$. The centered sum converges almost surely by the summable variances $`\sum m^{-2}`$, and the harmonic sum minus $`\log n`$ converges, so $`\tau_n-\log n`$ converges to a finite real random variable. In particular $`\tau_n\to\infty`$ and $`n e^{-\tau_n}\to W_{\rm size}>0`$. The size process is rate-one Yule from $`m_0`$: it is distributed at time $`t`$ as the sum of $`m_0`$ geometric variables on $`\{1,2,\ldots\}`$ of parameter $`e^{-t}`$, by the pure-birth forward equations. It has all finite-time polynomial moments and mean $`m_0e^t`$. Since $`\tau_{n+1}-\tau_n\to0`$, interpolation also gives $`e^{-t}|T(t)|\to W_{\rm size}`$. The geometric law yields $`W_{\rm size}\sim\mathop{\mathrm{Gamma}}\nolimits (m_0,1)`$ (shape, rate) by taking Laplace transforms. These are the clock and size facts of [`existence.md`, Lemmas 1.1–1.4](../../proofs/existence.md#1-continuous-time); we use the explicit independent construction above.

From size two the first birth produces the three-vertex path up to isomorphism. We start the continuous process at that path, denoted $`T_3`$, and use arbitrary finite seeds of size at least three in conditional arguments. The discrete chain from size three has the required law. For the mean logarithmic exponent, if $`F_2,F_3`$ are the mean core counts from $`T_2,T_3`$, respectively, then $`F_2(t)\le F_3(t)\le F_2(t+1)/(1-e^{-2})`$ by the independent first wait and monotonicity. Thus the mean logarithmic exponent from $`T_3`$ is still $`\mu`$.

## 2. Marked fringes and weighted renewal

Fix any deterministic finite physical tree $`T`$ with at least one core and orient it away from a core founder $`o`$. Write $`c_x`$ for core degree and $`L_x`$ for number of leaf neighbours, so at a core $`D_x=c_x+L_x`$. An arrow to a leaf makes it a degree-two core; this is a **promotion**. A core promotes at rate $`b_x=L_x/D_x\le1`$. Degrees and core status never decrease. For each old nonfounder core choose any mark in $`\{1/i:2\le i\le20\}`$ at least the reciprocal degree of its parent (all marks $`1/2`$ work). At a promotion by a core of degree $`D`$ give the promoted vertex the permanent mark $`a=1/m(D)`$, where $`m(D)=\min(D,20)`$. The inequality $`a\ge1/D_{\rm parent}`$ persists as the parent degree grows. Marks never change physical rates.

For a nonfounder core $`v`$ let $`S_v`$ be its downstream tree, with a single inert exterior parent slot at its root, counted in **both** total and core degree. The slot sends nothing, and a root choice of the slot does nothing in the fringe. All other arrows have their ordinary physical rates. Let $`L^0`$ be this intrinsic generator on finite planted marked types. For a fringe $`S`$ use the notation $`a_S,b_S,p_S`$ for the root mark, promotion rate and reciprocal degree. Sum over immediate core child fringes $`Y`$ with multiplicity. The type $`Y^+`$ adds one leaf to the root of $`Y`$. The newborn $`J_a`$ has one structural leaf and the inert parent slot.

### Positive families and their uniform bound

In the residual particle family let each type evolve intrinsically, kill it at rate $`a_S`$, and emit additional particles of types $`Y^+`$ and $`Y`$ at rates $`p_S`$ and $`a_Y-p_S`$, respectively. The parent remains present at a copy birth. Let $`R_t`$ be its first-moment semigroup on scores, with generator

```math
Rf(S)=L^0f(S)-a_S f(S)+\sum_Y\big[p_Sf(Y^+)+(a_Y-p_S)f(Y)\big].
```

Let $`Jf=(f(J_{1/i}))_{i=2}^{20}`$ and $`(B_0u)(S)=b_Su_{m(D_{\rm root})}`$. The full family also emits $`J_{1/m(D_{\rm root})}`$ at rate $`b_S`$; denote its first-moment semigroup by $`P_t`$, with generator $`R+B_0J`$. The intrinsic promotion in a particle and the additional direct birth may use separate clocks: the identities below are first-moment identities, not assertions of independence of physical core fringes. This is the construction of [`renewal.md` §§1–2](../fringe/renewal.md#1-marked-descendant-trees), with the capped marks of [`finite-kernel.md` §1](../fringe/finite-kernel.md#1-nineteen-newborn-labels).

Set $`d=3/20`$, $`q=9/25`$, $`\rho=q+d=51/100`$, $`z(a)=a/(a+d)`$, and

```math
h(0),h(1),h(2),h(3),h(c\ge4)=2,1,17/20,39/50,7/10.
```

For a fringe, give each core $`x`$ the product weight $`w_x`$ of $`z(a_y)`$ over the downward edges to $`x`$ (root weight one). Define $`W(S)=\sum_x w_x h(c_x)/(a_S+d)`$. For a full physical tree use $`Z(T)=\sum_x w_xh(c_x)`$ with the founder's true core degree and no slot. A leaf gain changes neither weight. At a promotion the parent changes $`c`$ to $`c+1`$ and the new core contributes relative weight $`(1+d m(D))^{-1}h(1)`$. The inequality used for every possible promoter, including the founder, is

```math
(1-c/D)\,[h(c+1)-h(c)+(1+d m(D))^{-1}]\le qh(c)
\qquad(D\ge2,\ 0\le c\le D).
```

For $`D\ge20`$ the left side is at most $`1/4<q(7/10)`$. For $`2\le D<20`$ the minimum of the right minus left over this indicated finite set is $`1/4350`$ at $`(D,c)=(3,1)`$ (substitute the five values of $`h`$). This is the capped version of [`renewal.md` §3, (4)–(6)](../fringe/renewal.md#3-a-subcritical-estimate-for-the-residual-process). Factoring the first edge gives $`\sum_Y a_Y W(Y)=(a_S+d)W(S)-h(c_{\rm root})`$. Since $`W(Y^+)=W(Y)`$, these identities give

```math
L^0W\le qW,\quad RW\le\rho W-h(c_{\rm root}),\quad
B_0JW\le1/d,\quad 1<W(S)\le |S|/d.
```

Stopping the families at finite event cutoffs gives the exponential supermartingale bound for their total weight with exponents $`\rho`$ and $`\rho+1/d`$, respectively. Since each particle carries weight at least one, the bound on the stopped weight bounds the probability of reaching a large particle count before any fixed time. With the particle count bounded, deaths have bounded total rate and the number of copy/direct births is bounded by deaths plus the maximal count. A finite set of intrinsic lineages has no explosion (each is dominated by vertex Yule growth). Thus neither family explodes, and Fatou gives $`R_tW\le e^{\rho t}W`$ and $`P_tW\le e^{(\rho+1/d)t}W`$. The same stopped promotion inequality on the physical tree gives $`\mathbb E_T Z(T(t))\le e^{qt}Z(T)`$, with $`Z(T)\le2C(T)`$.

### Physical first-moment identity, with all initial fringes

For the founder, which has no exterior slot, put

```math
\mathcal I(T)f=b_o f(J_{1/m(D_o)})+\sum_{Y\text{ at }o}[(a_Y-p_o)f(Y)+p_o f(Y^+)],\qquad
I_{T,u} f=\mathbb E_T\mathcal I(T(u))f.
```

An old nonfounder fringe omits exactly the exterior arrow from its parent into its root. Move the missing contribution $`p_o[f(Y^+)-f(Y)]`$ into the parent's sum (use the actual reciprocal degree at the relevant parent). The positive copy terms add $`a_Y f(Y)`$, cancelling the artificial kill at that child. The founder's contributions are precisely $`\mathcal I`$, and each promotion adds the new summand once. The weighted child sum at the founder is $`Z-h(c_o)`$, hence $`I_{T,u}W\le1/d+2C(T)e^{qu}`$. Applying the generator identity to the time-dependent test $`P_{t-u}f`$ gives, for $`\|f\|_W=\sup |f|/W<\infty`$,

```math
\mathbb E_T\sum_{v\ne o\text{ core at }t}f(S_v(t))
=\sum_{v\ne o\text{ core at }0}P_tf(S_v)+\int_0^t I_{T,u}(P_{t-u}f)\,du. \qquad\text{(1)}
```

For justification, for each fixed finite type the first-event integral equation for $`P_sf`$ has finitely many initial successor types and gives a continuous pointwise solution and the backward derivative. On $`0\le s\le t`$, a type of size $`n`$ has initial rate $`O(n)`$ and each one-event sum of successor weights is $`O(n+1)`$, giving the bound $`O_t((1+n)^2)\|f\|_W`$ on the derivative. In the stopped physical jump formula the sums, their time derivatives and absolute jump intensities are thus bounded by polynomials of the physical size, integrable on bounded time intervals by the Yule moments. Remove stopping by domination. This proves (1) for signed or complex scores as well. In [`renewal.md` §2, (2)–(3)](../fringe/renewal.md#2-a-positive-first-moment-system) the initial sum vanishes because the seed is $`T_3`$; it must not be deleted in (1).

## 3. The resolvent and the contour estimate

Use the weighted supremum space of fringe scores and the maximum norm on the nineteen newborn coordinates (with the induced row-sum norm for matrices). For $`\Re z>\rho`$ set

```math
G_z f(S)=\int_0^\infty e^{-zt}R_tf(S)\,dt,\quad
U_z=G_z B_0,\quad k(t)=JR_tB_0,\quad K(z)=JU_z=\int_0^\infty e^{-zt}k(t)\,dt.
```

The integrals on the weighted score space are defined **componentwise**. The bound $`\|G_z\|_W\le1/(\Re z-\rho)`$ is uniform in the seed. Integrating the exponential Taylor series about $`z_0`$ gives coefficient operator norms at most $`(\Re z_0-\rho)^{-j-1}`$ at order $`j`$, proving norm holomorphy without assuming strong continuity on the weighted supremum space. In particular $`K`$ is a holomorphic matrix there.

For real $`s>\rho`$ the entries of $`K(s)`$ are strictly positive: the first particle can survive and gain leaves up to any target degree $`2,\ldots,20`$ without promotion, then spend positive time at that degree with positive direct birth rate. Each entry is strictly decreasing in $`s`$. The positive-matrix Perron–Frobenius theorem (a finite strictly positive matrix has a simple positive spectral radius and positive left/right vectors) gives continuity and strict decrease of $`r(K(s))`$ (compare on the positive right eigenvector at the larger $`s`$). To identify the crossing, group the full family by direct birth generations. For its discounted time-integrated count started at the newborns, the result is $`\sum_{j\ge0}K(s)^j JG_s1`$. The forcing is bounded and strictly positive, by survival of the original residual particle. This sum is finite if $`r(K(s))<1`$ and infinite in every component otherwise, by comparison with the positive Perron vector. If finite, the discounted full-family count is bounded by a multiple of $`W`$ at every planted seed, by the first-direct-birth decomposition. Formula (1) then makes the discounted physical count from $`T_3`$ finite, since $`s>\rho>q`$. If the sum is infinite, the founder has positive direct immigration over a time interval, so (1) and Tonelli make the physical integral infinite. The previously established mean logarithmic exponent separates convergence for $`s>\mu`$ from divergence for $`\rho<s<\mu`$. Continuity gives $`r(K(\mu))=1`$, without using Claim 6 to identify the crossing. This recovers the particular determinant fact needed from [`2ccd472`, `finite-kernel.md` §1](../fringe/finite-kernel.md#1-nineteen-newborn-labels).

### High frequencies, including arbitrary root degrees

For $`\|u\|_\infty\le1`$, put $`F=B_0u`$. The rate times the absolute increment of $`F`$ under root leaf gains is bounded by 40 for degrees $`D<20`$. For $`D\ge20`$ the label in $`F`$ stays fixed and $`|\Delta F|\le c/[D(D+1)]`$, where $`c`$ is the root core degree, including the slot. The root gain rate is $`L+\sum_{Y}1/D_Y\le D`$, so this product is at most one. A root promotion changes $`b`$ by $`1/D`$, at rate at most one, and other intrinsic events contribute zero. Killing contributes at most $`1/2`$ and absolute copy contributions are at most $`\sum_Y a_Y W(Y)\le(13/20)W(S)`$. The predictable absolute jump intensity of the total score in a residual family is therefore bounded by $`42`$ times its total weight. This is a bound on nonzero score changes, not on the total intrinsic jump rate uniformly over the seed.

For a finite deterministic horizon the expected total variation of the sum of this score is finite, by $`R_tW\le e^{\rho t}W`$. Integrate the stopped jump formula and let the cutoff increase; the integrable total variation dominates the stopped score. This gives $`R_tF=F+\int_0^t R_u RF\,du`$ componentwise. Integration by parts, with $`\|RF\|_W\le42`$, now gives for $`\epsilon>0`$ and $`\Re z\ge\rho+\epsilon`$

```math
K(z)u=\frac{JB_0u+\int_0^\infty e^{-zt}JR_t R B_0u\,dt}{z},\qquad
\|K(z)\|\le\frac{1+280/\epsilon}{|z|}. \qquad\text{(2)}
```

Here $`JW\le1/d`$. The same integration before applying $`J`$ yields an $`O_\epsilon(1/|z|)`$ operator bound for $`U_z`$ from newborn vectors to weighted scores. No horizontal truncation constant grows with the initial degree. The degree cap is important: an arbitrary bounded score on uncapped degree labels could oscillate from $`D`$ to $`D+1`$.

### The only pole near the leading line

Choose positive Perron vectors $`g,\ell`$ of $`K(\mu)`$ with $`\ell g=1`$. The simple eigenvalue has analytic continuation $`\lambda(z)`$ near $`\mu`$: the resolvent integral around a small circle enclosing only the eigenvalue one gives an analytic rank-one projection. Differentiating a right vector normalized by $`\ell`$ yields $`\lambda'(\mu)=\ell K'(\mu)g<0`$, since $`-K'(\mu)=\int_0^\infty t e^{-\mu t}k(t)dt`$ and the kernel has positive mass at positive ages. Thus $`(I-K(z))^{-1}`$ has a simple pole at $`\mu`$ with residue $`D_*=(g\otimes\ell)/[-\ell K'(\mu)g]`$; the complementary inverse is analytic. For $`\Re z>\mu`$, $`|K(z)^j v|\le K(\Re z)^j|v|`$, so its Neumann inverse exists.

If $`K(\mu+i\tau)v=v`$ for $`\tau\ne0`$ and $`v\ne0`$, then $`|v|\le K(\mu)|v|`$. Multiplication by the strictly positive left vector forces equality in every coordinate, so $`|v|`$ is a positive multiple of $`g`$. In particular $`v_2\ne0`$. Every row's newborn label-2 contribution has density at least $`\tfrac12 e^{-(2+1/i)t}`$ at $`t>0`$ (the first particle has no clock ring and survives death). Its phases $`e^{-i\tau t}v_2`$ vary on every interval. This gives a strict triangle inequality for that contribution and hence $`|K(\mu+i\tau)v|<K(\mu)|v|`$, a contradiction. The positivity-of-modulus step is essential before using the label-2 density.

Take $`\epsilon=(\mu-\rho)/2`$ and a frequency threshold so that (2) bounds $`\|K(z)\|\le1/2`$ uniformly above that threshold in $`\Re z\ge\rho+\epsilon`$. On the bounded portion of the line $`\Re z=\mu`$ outside a small disk about $`\mu`$, invertibility extends to a uniform neighbourhood by compactness. Inside the disk the simple-pole representation excludes other singularities. Thus there exists a single $`\sigma\in(\rho+\epsilon,\mu)`$ such that in $`\Re z\ge\sigma`$ the only singularity of the inverse is $`\mu`$. The inverse is bounded on $`\Re z=\sigma`$, by compactness at bounded frequencies and its Neumann bound at high frequencies. These choices depend only on the nineteen-type kernel.

### Inversion and the uniform mean remainder

Write $`v(t)=\sum_{n\ge1}k^{*n}(t)`$, with matrix multiplication in time order. The kernel $`k`$ is continuous, and on $`[0,T]`$, if $`M_T=\sup\|k\|`$, the $`n`$th term has norm at most $`M_T^n t^{n-1}/(n-1)!`$. Thus $`w=v-k`$ is continuous with $`w(0)=0`$. At a fixed real $`s>\mu`$ the extension by zero of $`e^{-st}w(t)`$ is continuous and integrable on the real line, by the nonnegative matrix generation series. Its Fourier transform, evaluated at frequency $`\tau`$, is

```math
F_w(s+i\tau),\qquad F_w(z)=(I-K(z))^{-1}-I-K(z)=(I-K(z))^{-1}K(z)^2.
```

Uniformly over $`\sigma\le x\le s`$ at large $`|\tau|`$, $`\|F_w(x+i\tau)\|\le 2(1+280/\epsilon)^2/\tau^2`$ by (2). In particular the transforms on both vertical lines are integrable. Fourier inversion for a continuous $`L^1`$ function with $`L^1`$ transform gives the vertical integral on $`\Re z=s`$ for every $`t\ge0`$. Shift it on the rectangle of heights $`\pm L`$: the horizontal integrals have norm $`O((s-\sigma)e^{st}/L^2)`$ and tend to zero for each $`t`$. The residue is $`e^{\mu t}D_*`$; the remaining vertical integral has norm at most $`e^{\sigma t}(2\pi)^{-1}\int_{\mathbb R}\|F_w(\sigma+i\tau)\|\,d\tau`$. Adding back $`k(t)=O(e^{\rho t})`$ proves $`v(t)=e^{\mu t}D_*+O(e^{\sigma t})`$ for every $`t\ge0`$. Inverting the unsubtracted transform by an absolute-integrability argument would lack the second power of frequency decay.

The first-direct-birth equations are $`JP_tf=JR_tf+(v*JR_\cdot f)(t)`$ and $`P_tf=R_tf+\int_0^t R_u B_0 JP_{t-u}f\,du`$. Define

```math
d(f)=\frac{\ell JG_\mu f}{-\ell K'(\mu)g},\qquad H=U_\mu g.
```

In the newborn equation $`\|JR_u f\|\le(1/d)e^{\rho u}\|f\|_W`$. Convolving the error costs at most a constant times $`e^{\sigma t}/(\sigma-\rho)`$; extending the leading integral to infinity costs at most a constant times $`e^{\rho t}/(\mu-\rho)`$. Thus $`JP_t f=e^{\mu t}g d(f)+O(e^{\sigma t}\|f\|_W)`$. In the seed equation $`|(R_u B_0v)(S)|\le e^{\rho u}W(S)\|v\|_\infty`$, so exactly those two integrals give

```math
|P_tf(S)-e^{\mu t}H(S)d(f)|\le K_0 e^{\sigma t}W(S)\|f\|_W. \qquad\text{(3)}
```

In resolvent language the full Laplace resolvent on $`\Re z>\mu`$ is $`G_z+U_z(I-K(z))^{-1}JG_z`$, by the same first-birth decomposition and positive domination. Its simple-pole part is $`H\otimes d/(z-\mu)`$; (3) gives an analytic remainder on $`\Re z>\sigma`$ with weighted operator norm at most $`K_0/(\Re z-\sigma)`$. No inversion of an unbounded-seed generator in an unweighted norm is used. The contour is taken only for the finite newborn matrix. These are the uniform versions of [`almost-sure.md` §1, (2)–(3), `c98885f`](../fringe/almost-sure.md#1-the-finite-renewal-kernel-has-an-exponential-remainder).

## 4. The physical eigenfunction and its comparison to core count

Put

```math
\Psi(T)=\sum_{v\ne o\text{ core at }0}H(S_v)+\int_0^\infty e^{-\mu u}I_{T,u}H\,du.
```

For fixed $`x`$, its ancestral weights in the sum of fringe weights are bounded by $`(2/d)`$ times a geometric sequence of ratio at most $`z(1/2)=10/13`$. Thus $`\sum_{v\ne o} W(S_v)\le K_1 C(T)`$ with $`K_1=(2/d)/(1-10/13)`$. Since $`C(T)\ge1`$, the immigration bound from Section 2 becomes $`I_{T,u}W\le K_I C(T)e^{qu}`$ with $`K_I=1/d+2`$. Neither constant depends on seed degrees or leaf stock.

Insert (3) into (1) and extend the leading immigration integral to infinity. For a score $`f`$, the initial-fringe error is at most $`K_0K_1 C(T)e^{\sigma t}\|f\|_W`$. The error in the immigration convolution and the discarded leading tail are bounded, respectively, by

```math
\frac{K_0 K_I}{\sigma-q} C(T)e^{\sigma t}\|f\|_W,
\qquad
\frac{K_I\|H\|_W\|d(\cdot)\|}{\mu-q} C(T)e^{qt}\|f\|_W,
```

where $`\|d(\cdot)\|`$ is the norm of the coefficient functional on weighted scores. The denominators come from $`\int_0^t e^{-(\sigma-q)u}du`$ and $`e^{\mu t}\int_t^\infty e^{-(\mu-q)u}du=e^{qt}/(\mu-q)`$. For $`A_k(T)=N_{\ge k}(T)`$ and $`f_k(S)=\mathbf1_{\{D_{\rm root}\ge k\}}`$, we have $`\|f_k\|_W\le1`$. All nonfounder fringe root degrees are their true physical degrees; the founder contributes at most one, absorbed by $`C(T)e^{\sigma t}\ge1`$. Thus with $`d_k=d(f_k)`$,

```math
|\mathbb E_T A_k(T(t))-e^{\mu t}\Psi(T)d_k|\le K_2 C(T)e^{\sigma t}, \qquad\text{(4)}
```

with one constant for all finite seeds with a core and all $`k\ge2`$. Here $`d_2=d(1)>0`$ by survival of an original newborn, and $`0\le\Psi(T)\le K_3 C(T)`$ by the displayed bounds. In particular the true seed error is proportional to the number of **cores**, not the vertex count. This distinction is later needed when passing from conditional means to almost-sure convergence.

For the reverse comparison, call a core eligible if it has $`c\le2,L\ge1`$, or if it is bare ($`c=D=2`$) and has a bare neighbour. For $`C=1`$ the centre is eligible. For $`C\ge2`$, the number of core vertices with $`c\ge3`$ is at most the number with $`c=1`$ minus two (sum $`c-2`$ over the core tree). The isolated bare vertices with no bare neighbour use two distinct core edges each, so their number is at most $`(C-1)/2`$. All $`c=1`$ vertices are eligible; these bounds give at least $`(C+5)/4`$ eligible. This is the counting lemma in [`existence.md` §2, Lemma 2.1](../../proofs/existence.md#2-eligible-vertices).

A nonfounder eligible vertex of the first kind has a uniform lower bound on $`H`$. In the original residual particle require survival through time one, and let the first root clock ring occur in that unit. At that ring there is still at least one leaf and $`c\le2`$, so the probability of a promotion is at least $`1/3`$. The expected discounted root-promotion payment, with newborn coefficient at least $`\min_i g_i`$, is a lower bound for $`H`$ by compensating the intrinsic promotion clock (the direct-birth intensity is the same). It is at least $`e^{-(\mu+1/2)}(1-e^{-1})\min_i g_i/3`$. A bare nonfounder with a bare child has a lower bound too: silence the root and that child's other neighbour for a unit of time, require the child's first ring to select the root during that unit, and require independent residual survival. The child remains of degree two before its first choice because only its two neighbours can increase its degree. The root now has a leaf and can use the first bound in the next unit. These probabilities do not depend on unobserved degrees. The bare eligible vertices with no bare child can be charged to bare parents with a bare child; each such parent receives at most two charges, and a single founder may have to be excluded. Thus there are at least $`C/12-1`$ nonfounders covered by the two bounds. This is the reverse comparison of [`renewal.md` §5](../fringe/renewal.md#5-constant-factor-growth-of-the-ordinary-mean) with the same local scripts under the capped marks. It gives $`C(T)\le K_4(1+\Psi(T))`$. For $`T_3`$, $`\Psi(T_3)>0`$ because there is positive founder direct immigration and $`H(J_{1/i})=g_i>0`$.

By (4) at $`k=2`$, $`\Psi(T)=d_2^{-1}\lim_{t\to\infty}e^{-\mu t}\mathbb E_T C(T(t))`$, a physical quantity independent of the chosen founder and compatible marks. Conditioning at any fixed elapsed $`u`$, apply (4) from the random tree $`T(u)`$ and send the remaining horizon to infinity. The mean of the error is at most $`K_2e^{-(\mu-\sigma)t}\mathbb E_T C(T(u))`$, which tends to zero because that expectation is bounded by $`|T|e^u`$. Hence $`\mathbb E_T\Psi(T(u))=e^{\mu u}\Psi(T)`$. The process $`M_t=e^{-\mu t}\Psi(T(t))`$ is a nonnegative càdlàg martingale by the Markov property. This derives the eigenfunction of the physical process before taking a random normalized-count limit.

## 5. Concentration and convergence at all times

### A uniform one-leaf difference

Compare a tree with a core to the same tree with a new leaf. If the attachment target is a core, mark its two copies dirty. If it was a leaf, mark the two copies of its old parent and the newly promoted core on the larger side. At most three dirty cores exist initially, counted across both copies. Couple rings of clean matched cores and their paired structural leaves, matching the neighbour choices. A clean core has matched leaf stock and matched neighbours, including a corresponding slot for each adjacent dirty core. Dirty cores ring independently. A dirty ring can promote a private core (+1 dirty) or send into a clean pair (mark **both** recipients dirty, +2); paired leaf rings into a clean root are synchronized, and leaf rings into a dirty root cannot spread dirt. Clean paired promotions create clean pairs. These rules preserve the clean invariants: a new dirty core born at a dirty parent is attached only to that parent, and making a clean pair dirty changes the status of the same neighbour slots at every remaining clean core. The conditional drift of the combined dirty count is at most twice its current value. Stop at bounded sizes and use Yule domination and Grönwall. For elapsed $`u`$ the expected count is at most $`3e^{2u}`$. All clean degrees agree, so

```math
|\mathbb E_{T^+} A_k(T(u))-\mathbb E_T A_k(T(u))|\le3e^{2u}. \qquad\text{(5)}
```

In particular one must not use a bound starting at two dirty cores for a leaf promotion or charge only one new dirty core when a clean recipient is hit. The proof uses a coupling of the actual physical trees, not the auxiliary positive branching law.

### Quadratic variation and an $`L^2`$ bound uniform over seeds

For a fixed observation time $`b`$, the backward conditional mean of $`A_k(T(b))`$ is a square-integrable martingale. Its jump at time $`u`$ is a one-leaf difference at horizon $`b-u`$, so (5) bounds the squared jump by $`9e^{4(b-u)}`$. Between jumps the conditional mean is differentiable by the first-event equation; on a finite interval the martingale has no continuous martingale part. The square-integrable martingale isometry and the jump compensator give

```math
\mathop{\mathrm{Var}}\nolimits _T A_k(T(b))\le9|T|\int_0^b e^{4(b-u)+u}\,du. \qquad\text{(6)}
```

This compensation can first be stopped at bounded sizes and then extended using the finite-time Yule second moment. Let $`\eta=\mu-\sigma>0`$ and $`g_h(T)=d_2^{-1}e^{-\mu h}\mathbb E_T C(T(h))`$. Apply the same backward calculation to $`g_h(T(1))`$ (the conditional horizon at $`u\le1`$ is $`h+1-u`$). This gives variance at most $`K_5 e^{4h}|T|`$. By (4), $`|\Psi-g_h|\le K_6 e^{-\eta h}C`$. The physical core is dominated by rate-one Yule from $`C(T)`$ since each core promotes at rate at most one, so $`\mathbb E_T C(T(1))^2\le 2e^2 C(T)^2`$. Using the variance inequality for a sum gives

```math
\mathop{\mathrm{Var}}\nolimits _T\Psi(T(1))\le K_7(e^{4h}|T|+e^{-2\eta h}C(T)^2). \qquad\text{(7)}
```

At each integer $`j`$ choose $`h=\eta^{-1}\log(j+2)`$. The conditional variance of $`M_{j+1}`$ adds to the previous second moment. Use (7), the reverse core comparison, and $`\mathbb E_T |T(j)|=|T|e^j`$ to obtain

```math
\mathbb E_T M_{j+1}^2\le (1+K_8(j+2)^{-2})\mathbb E_T M_j^2
+K_9\big[|T|(j+2)^{4/\eta}e^{-(2\mu-1)j}+(j+2)^{-2}e^{-2\mu j}\big].
```

The products and additive sums are finite. Conditional Jensen between deterministic integers extends this bound to all $`t\ge0`$. The $`L^2`$ martingale convergence theorem for a càdlàg martingale bounded in $`L^2`$ gives almost-sure and $`L^2`$ convergence to $`M_\infty`$, with

```math
\mathbb E_T M_\infty=\Psi(T),\qquad \mathbb E_T M_\infty^2\le K_{10}(\Psi(T)^2+|T|+1), \qquad\text{(8)}
```

uniformly in deterministic seeds. A fixed $`h`$ in (7) would not make the multiplicative error in this proof summable; that does not show the actual martingale is unbounded for some other proof strategy.

### Strict positivity, simultaneous counts, and interpolation

On the process from $`T_3`$, the earlier logarithmic theorem and the Yule clock give $`\log C(T(j))/j\to\mu`$ and $`\log |T(j)|/j\to1`$. The core comparison implies $`\Psi(T(j))\to\infty`$ and $`|T(j)|/\Psi(T(j))^2\to0`$. Conditional on the history up to $`j`$, apply Cauchy–Schwarz to the nonnegative future martingale limit from $`T(j)`$ and (8). The conditional probability that $`M_\infty>0`$ is at least

```math
\frac{\Psi(T(j))^2}{K_{10}(\Psi(T(j))^2+|T(j)|+1)}.
```

The event is measurable in the sigma-field generated by the increasing histories. The upward convergence theorem for conditional expectations makes these conditional probabilities converge to its indicator; their lower limit is positive. Thus $`M_\infty>0`$ almost surely, with no assumption of positivity of the desired $`X`$.

Fix $`k\ge2`$ and a mesh of spacing $`1/m`$. At mesh time $`t`$ put $`h=\eta^{-1}\log(t+2)`$, $`u=t-h\ge0`$ for all sufficiently large $`t`$. Conditional on time $`u`$, (6) bounds the unconditional squared centered error at normalization $`e^{-\mu t}`$ by

```math
K_{11}|T_3|e^{-(2\mu-1)t+4h}.
```

This is summable on the mesh. Chebyshev and Borel–Cantelli (at countably many error thresholds) make the centered error tend to zero almost surely. By (4) the scaled conditional mean differs from $`d_k M_u`$ by at most $`K_2 e^{-\eta h}e^{-\mu u}C(T(u))`$. The reverse comparison bounds the last normalized count by $`K_4(e^{-\mu u}+M_u)`$, almost surely bounded as $`u\to\infty`$. Hence $`e^{-\mu t}A_k(T(t))\to d_k M_\infty`$ on the mesh. Intersect over all integer $`k,m`$. Since each cumulative count is nondecreasing, the values at two successive mesh points bound the rescaled intervening value with factors $`e^{\pm\mu/m}`$. Letting $`m\to\infty`$ proves convergence at **all** real times, including between and at random Yule epochs. Exact-degree counts follow by subtraction, not by monotonicity of exact-degree counts.

## 6. Fixed sizes, limiting coefficients, and support

Along the Yule epochs from $`T_3`$ we obtain

```math
\frac{A_k(T_n)}{n^\mu}\longrightarrow d_k W_{\rm size}^{-\mu}M_\infty,
\qquad X=d_2 W_{\rm size}^{-\mu}M_\infty.
```

This is a pathwise substitution, not a passage of expectation through the random denominator. The fixed-size limit $`X`$ is measurable from the discrete jump chain, and is therefore independent of $`W_{\rm size}`$. The continuous limit $`M_\infty`$ is **not** independent of $`W_{\rm size}`$: its conditional mean given $`W_{\rm size}`$ is $`W_{\rm size}^\mu\mathbb E X/d_2`$, which is nonconstant. Here $`0<\mathbb E X<\infty`$ follows either from the discrete martingale argument below or from independence and $`0<\mathbb E M_\infty<\infty`$ by Tonelli.

To control expectation and the support, the first holding time from size $`m`$ is independent of the successor. Laplace-integrating the eigenfunction identity at any real $`s>1`$ and splitting at that first time gives

```math
\frac{\Psi(T)}{s-\mu}=\frac{\Psi(T)}{s+m}+\frac{m}{s+m}\frac{\mathbb E[\Psi(T_{m+1})\mid T_m=T]}{s-\mu}.
```

Thus $`\mathbb E[\Psi(T_{m+1})\mid T_m]=(1+\mu/m)\Psi(T_m)`$ and $`Y_m=\Gamma(m)\Psi(T_m)/\Gamma(m+\mu)`$ is a discrete martingale. At the stopping time $`\tau_m`$, use conditional Jensen for the $`L^2`$-closed continuous martingale (or bounded stopping followed by Fatou). Independence of the jump chain and the epochs gives

```math
\left(\prod_{j=3}^{m-1}\frac{j}{j+2\mu}\right)\mathbb E\Psi(T_m)^2
=\mathbb E M_{\tau_m}^2\le\mathbb E M_\infty^2.
```

The product equals $`\Gamma(3+2\mu)\Gamma(m)/(\Gamma(3)\Gamma(m+2\mu))\asymp m^{-2\mu}`$ by the gamma recurrence and the asymptotic gamma ratio; the $`\Gamma(3)`$ records the initial size three. Thus $`Y_m`$ is bounded in $`L^2`$ and its limit is $`X/d_2`$ by the epoch identity. In particular, with the discrete history through size $`m`$,

```math
\mathbb E[X\mid T_3,\ldots,T_m]=d_2\frac{\Gamma(m)}{\Gamma(m+\mu)}\Psi(T_m).
```

The comparison $`C\le K_4(1+\Psi)`$ makes $`n^{-\mu} A_k(T_n)`$ uniformly bounded in $`L^2`$, so its almost-sure limit also identifies the leading mean coefficient:

```math
c_k=d_k\frac{\Gamma(3)}{\Gamma(3+\mu)}\Psi(T_3),\qquad c_2=\mathbb E X.
```

This also identifies the constants of [`8e4db79`, `asymptotics.md` §3, (1) and (5)](../fringe/asymptotics.md#3-returning-to-ordinary-time-and-to-a-fixed-number-of-vertices), whose conclusion for the identical single-edge model is $`\mathbb E N_{\ge k}(n)\sim c_kn^\mu`$ for each fixed $`k\ge2`$. We have obtained that conclusion directly without importing its Tauberian step. Positivity of $`d_k-d_{k+1}`$ follows by keeping an original residual newborn alive and obtaining exactly $`k-2`$ root leaf gains, then waiting without changes for an interval: the residual integral of the exact-degree indicator is strictly positive. Dominated convergence gives $`d_k\downarrow0`$ since every type has finite degree, $`0\le f_k\le1\le W`$, and $`JG_\mu W`$ is finite. This proves the ratio limits and normalization of $`p_k`$.

Stars and paths occur at every prescribed finite size with positive probability. From a size-$`j`$ star, sending from any leaf into the centre preserves the star, with probability $`(j-1)/j`$. From a size-$`j`$ path, sending into either endpoint preserves the path, with total probability $`1/j`$: the source next to each endpoint has degree two, even when both sources are the same vertex at $`j=3`$. A broken star or path cannot return to that shape, since core count and maximum degree never decrease. The two events therefore have probabilities $`2/(m-1)`$ and $`2/(m-1)!`$ at size $`m\ge3`$. The star has one core and the path has $`m-2`$. At a star the displayed conditional mean of $`X`$ is $`O(m^{-\mu})`$ by the upper comparison for $`\Psi`$, so for sufficiently large $`m`$ that conditional mean is smaller than $`\varepsilon`$. At a path it is at least a positive constant times $`m^{1-\mu}`$ for sufficiently large $`m`$, by the reverse comparison, so it exceeds $`M`$. Since each conditioning event has positive probability, these conditional mean inequalities force the stated two support properties. This use of $`\mu<1`$ and integrability is separate from the proof of strict positivity.

## 7. Identifying the large-degree tail of the limiting law

This identification imports two scoped spatial results. In a finite or countable, connected, locally finite rooted core tree with leaf counts $`L`$ and total degrees $`L+c\ge2`$, let $`Q_t f`$ be the expected sum of $`f`$ rerooted at the descendants of the old root's initial core cell. Here a new core inherits the cell of its promoter, and rates are the full physical rates, including across edges between old cells. [`operator.md` §§1–2, (2)–(5), `8e4db79`](../operator.md#2-a-probability-eigenmeasure-at-the-growth-rate) constructs this Markov growth and, for the model with logarithmic mean exponent $`\mu>0`$, obtains a probability $`\pi`$ with $`\pi Q_t=e^{\mu t}\pi`$ for all $`t\ge0`$ and $`\pi(e^{\theta c_o})<\infty`$ for $`0<\theta<\log(1+\mu)`$. Specifically this $`\pi`$ is a weak subsequential limit of the normalized mean measures $`\int_0^\infty w(t)\mathbb E_{T_3}\sum_{x\text{ core}}\delta_{(T(t),x)}dt`$, with $`w(t)=e^{-\mu t-\delta|t-T|}`$ and parameters chosen so $`\delta\downarrow0`$ and $`w(0)/\int w(t)\mathbb E C(T(t))dt\to0`$.

For this same cell process, [`degree-law.md`, (1)–(2), §Point probabilities, `8e4db79`](../degree-law.md) assumes $`\mu>0`$, a probability eigenmeasure satisfying the preceding identity for all $`t`$, and some finite exponential moment of $`c_o`$. It concludes $`\pi(D_o\ge k)\sim A k^{-\mu}`$ and $`\pi(D_o=k)\sim\mu A k^{-1-\mu}`$ with $`0<A<\infty`$. Its constant is $`A=\mathbb E_\beta W_{\rm deg}^{\mu}`$, with $`\beta`$ the eigenmeasure's promotion-biased newborn environment, and $`W_{\rm deg}=\lim_{t\to\infty} e^{-t}D_o(t)`$ for a tracked newborn in that environment. These results predate Claim 6; neither requires an almost-sure empirical local limit.

To see that this is our degree distribution, from (4) the ratio of cumulative mean counts from $`T_3`$ tends to $`d_k/d_2`$. In the normalized mean tilts defining $`\pi`$, the triangle inequality gives $`w(t)/w(0)\le e^{(-\mu+\delta)t}`$. With $`0<\delta<\mu`$ and $`Z_{\delta,T}=\int_0^\infty w(t)\mathbb E_{T_3} C(t)\,dt`$, the mass of each fixed interval $`[0,B]`$ is at most

```math
\frac{w(0)}{Z_{\delta,T}}\int_0^B \mathbb E_{T_3} C(t)\,dt\longrightarrow0.
```

Thus the same ratio holds in the tilt limit. The degree-threshold test is bounded and continuous in the topology of marked rooted balls, so $`\pi(D_o\ge k)=d_k/d_2`$ and $`\pi(D_o=k)=p_k`$. This proves the large-$`k`$ assertion about the distribution obtained **after** the fixed-$`k`$ limits, without claiming a uniform estimate at degrees growing with $`n`$.

## 8. Provenance, review scope, and replay

The primary source is [`c98885f`, `fringe/almost-sure.md`](../fringe/almost-sure.md), Sections 1–5. The prerequisite chain is `bfcbbe5` (positive fringe renewal), `8e4db79` (mean asymptotics, spatial eigenmeasure and degree law), and `2ccd472` (nineteen permanent marks). This proof needs no later numerical estimate of $`\mu`$ and no signed surrogate. The scope of the archived Gemini 4 Carbon **CORRECT** review is the five sections of the original Claim 6 argument; it takes the capped renewal construction as a dependency. Gemini 4 Neon **CORRECT** covers that finite-label determinant and rational truncations, not the random limit. The earlier Carbon mean-asymptotics review and degree-law review cover the respective leading-mean/transfer and eigenmeasure degree-tail arguments. These are reviews of the cited sources and fixtures, not a review of this standalone document.

No unresolved gap was found in the cited Claim 6 chain. The reconstruction supplies details where a shortcut would leave an incomplete argument: positivity of the modulus on the peripheral line, absolute inversion after subtracting $`k`$, every old nonfounder term in the conditional mean, and independence only for quantities that actually depend on disjoint parts of the clock construction.

The independent GPT-6.2 Sol Ultrafast cold review found the standalone draft **CORRECT** conditional on the stated earlier inputs; it is a same-family review. Its [exact fixture](../checks/standalone_6_sol.py) tested 323 rational star/crown endpoint rows through $`D=10^{60}`$, including fixed core count with widely separated leaf stocks, dense crowns, and opposite label scores across the cap seam. The largest tested capped gain variation was $`6489/190`$; the uncapped plant was rejected. Six separate exact rows checked exponent budgets for the review's rate extension, including $`\eta=10^{-9}`$, and rejected an oversampled mesh. The rate is derived in that review; the fixtures do not certify the physical spectrum or simulate asymptotic sample paths.

Commands run from the repository root with the Python standard library:

```sh
nice -n 10 taskset -c 8 timeout 35s python3 -Bu round2/checks/strong_limit_carbon.py
nice -n 10 taskset -c 8 timeout 20s python3 -Bu round2/checks/finite_kernel_neon.py
nice -n 10 taskset -c 13 timeout 20s python3 -Bu round2/checks/asymptotics_carbon.py
```

All three passed: the first reported 118 rooted marked configurations, 128 coupled states, 434 unlabelled eligibility fixtures, and 10/10 negative controls (archived output); the second verified the clamped weights and the 19-label rational truncations at sizes 2–4 and caught all 10 plants; the third reported 88 marked configurations and 8/8 negative controls. The first checker's contour experiment is floating-point quadrature on a **four-state, two-label toy model**. Its docstring mentions $`10^{-9}`$, but the actual reconstruction assertion uses $`2\times10^{-4}`$; the test is not a certified enclosure for the physical 19-label inverse or its spectral gap. The proof of the contour and all-seed uniformity is analytic in Sections 2–4. The convergence controls that vary the smoothing window or lower $`\mu`$ test failure of the stated estimates to yield the proof, not falsity of a limit in a different model.
