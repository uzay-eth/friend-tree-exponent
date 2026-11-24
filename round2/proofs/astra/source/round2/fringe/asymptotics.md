# The leading coefficient and the limiting degree frequencies in mean

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 3: Claude and Gemini; GPT standalone audit; inputs include Round 1 growth, the marked renewal and spatial weights (Claims 2 and 4). See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Statement submitted for independent review

Let $`N_{\ge k}(n)`$ count the vertices of degree at least $`k`$ when the tree has $`n`$ vertices, starting from a single edge. For each fixed $`k\ge2`$, the following argument gives a constant $`c_k>0`$ such that
```math
\mathbb E N_{\ge k}(n)\sim c_k n^\mu .
\qquad\text{(1)}
```
In particular $`c_2`$ is the leading coefficient for the mean core size. Write $`p_k=(c_k-c_{k+1})/c_2`$. Then $`p_k>0`$, $`\sum_{k\ge2}p_k=1`$, and
```math
\lim_{n\to\infty}\frac{\mathbb E N_{=k}(n)}{\mathbb E C_n}=p_k .
\qquad\text{(2)}
```
Together with the separate marked-core argument this gives $`p_k\sim\mu A k^{-1-\mu}`$ as $`k\to\infty`$, with $`0<A<\infty`$. Equations (1) and (2) concern means. The derivation uses the residual semigroup, birth operator and immigration measure defined in [the renewal construction](renewal.md).

## 1. Discounted occupation

Use the Banach space of scores on finite planted types with norm $`\|f\|_W=\sup_S |f(S)|/W(S)`$. Define $`f_k(S)=\mathbf1_{\{D_{\rm root}\ge k\}}`$. Write $`f=1`$ or $`f=f_k`$. For a complex discount $`z`$ with $`\Re z>\rho=51/100`$, put
```math
G_z f(S)=\int_0^\infty e^{-zt} R_t f(S)\,dt,\quad
\widehat I_z(F)=\int_0^\infty e^{-zt} I_t F\,dt,\quad
jF(a)=F(J_a).
```
Recall $`U_z h=G_z[b_S h(p_S)]`$ and $`K_z=jU_z`$. The bound $`R_tW\le e^{\rho t}W`$ gives norm-convergent Taylor expansions for $`G_z f`$ and $`U_z`$ about every such $`z`$, by integrating the scalar expansions termwise with a uniform exponential majorant. The bound on $`I_t W`$ gives the same conclusion for $`\widehat I_z`$ on $`\Re z>9/25`$. These definitions can all be made with componentwise integrals. For the scores used here $`jG_z f`$ is continuous in the birth mark. The first-particle calculation in Section 4 of the renewal construction applies with the source $`f`$, whose value along that particle's intrinsic path is independent of its own birth mark. The identical estimates with $`\Re z`$ in place of $`s`$ give compactness of $`K_z`$ and holomorphy with values in the bounded operators on $`C(\mathcal A)`$.

Let $`M_t f`$ count noninitial cores in the actual tree from $`T_3`$ as in the renewal construction. For $`\Re z>\mu`$ the integral is absolutely convergent and the decomposition by successive births gives
```math
\int_0^\infty e^{-zt}M_t f\,dt
=\widehat I_z\left(G_z f+U_z(I-K_z)^{-1} jG_z f\right).
\qquad\text{(3)}
```
Indeed the absolute values of all generation terms for $`z`$ are bounded by the corresponding terms at $`\Re z`$, for $`|f|`$. At that real argument $`r(K_{\Re z})<1`$, so the operator series converges in norm, also for $`z`$. Equation (3) then follows from the immigration convolution and Fubini.

## 2. The pole on the line of convergence

Let $`g`$ and $`\ell`$ be positive right and left eigenvectors of $`K_\mu`$ at eigenvalue one, normalised by $`\ell(g)=1`$. The dual eigenmeasure has full support. For a nonzero continuous $`v\ge0`$, strong positivity gives $`\ell(v)=\ell(K_\mu v)>0`$. Simplicity and compactness isolate this eigenvalue. Its analytic eigenvalue branch $`\lambda(z)`$ has
```math
\lambda'(\mu)=\ell(K'_\mu g)<0 .
```
The strict inequality follows from integrating the positive reproduction kernel times $`-t`$ over positive times. For completeness, take the spectral projection of $`K_z`$ on a small circle about 1. Its rank is one locally, by continuity in operator norm. The inverse of $`I-K_z`$ on the complementary projection is analytic nearby, and on the spectral line it is multiplication by $`1/(1-\lambda(z))`$. Differentiating the eigenvalue equation after normalising the eigenvector by $`\ell`$ gives the derivative above. The residue of the inverse is therefore $`g\otimes\ell/[-\ell(K'_\mu g)]`$.

Fix $`\tau\ne0`$. If $`K_{\mu+i\tau}v=v\ne0`$, then $`|v|\le K_\mu |v|`$. Applying $`\ell`$ and using its full support makes equality hold everywhere, so $`|v|`$ is a positive multiple of $`g`$. For any row $`a`$, the triangle inequality in the integral for $`K_{\mu+i\tau}v(a)`$ must attain equality. Each row contains positive mass for children with birth mark $`1/2`$ at times in every interval: take the original $`J_a`$, with no root or leaf event and no killing up to $`t`$. This gives the component at $`1/2`$ an absolutely continuous density bounded below on each compact interval in $`(0,\infty)`$. Since $`v(1/2)\ne0`$, the complex numbers $`e^{-i\tau t}v(1/2)`$ have varying phases on this component. Equality is impossible. Compactness hence gives invertibility of $`I-K_{\mu+i\tau}`$, which continues to a neighbourhood of this point by the ordinary inverse-series argument.

It follows that the transform in (3) minus $`d_f/(z-\mu)`$ extends analytically to a neighbourhood of every point of the line $`\Re z=\mu`$, with
```math
d_f=\frac{\widehat I_\mu(U_\mu g)\,\ell(jG_\mu f)}{-\ell(K'_\mu g)}.
\qquad\text{(4)}
```
The first factor in the numerator is strictly positive by direct immigration at the initial centre. The second is positive for every $`f_k`$, since the intrinsic original particle can acquire the required leaves before any promotion or death, and can then hold that degree for an interval. The same observation applies to $`f_k-f_{k+1}`$ supported at any single degree $`k`$.

## 3. Returning to ordinary time and to a fixed number of vertices

For $`f_k`$, the actual count $`M_t f_k`$ is nondecreasing in $`t`$. Each new core has degree two, and existing degrees never decrease. Apply Wiener–Ikehara to the function $`x\mapsto M_{x/\mu} f_k`$. Its Laplace transform is $`\mu`$ times (3) evaluated at $`\mu z`$, with residue $`d_{f_k}`$ at 1. Monotonicity implies the one-sided Tauberian condition, and local analytic extension of the pole-subtracted transform is stronger than the required boundary condition. Thus $`M_t f_k\sim d_{f_k} e^{\mu t}`$. The initial centre contributes at most one. To return to the initial edge, condition on its first birth, whose wait is exponential of rate two and independent of the subsequent jump chain. Convolution, using boundedness of the normalised $`T_3`$ mean, gives
```math
\mathbb E_{T_2}N_{\ge k}(t)\sim d_k^* e^{\mu t},\qquad
d_k^*=\frac{2}{2+\mu}d_{f_k}.
\qquad\text{(5)}
```
The Tauberian implication with local analytic extension follows, for example, from Theorem 1.1 and its subsequent paragraph in [Debruyne and Vindas](https://arxiv.org/abs/1611.09765).

For clarity, we give the fixed-size transfer. Put $`a_{n,k}=\mathbb E N_{\ge k}(n)`$ and $`q=e^{-t}`$. The holding times at size $`n`$ are exponential of rate $`n`$ and independent of the discrete jump chain; the size from two is the sum of two independent Yule populations from one. Thus the left side of (5) is
```math
q^2\sum_{n\ge2}(n-1)(1-q)^{n-2} a_{n,k}.
```
For $`x\to\infty`$ define the locally finite measure
```math
\nu_x=x^{-(\mu+2)}\sum_{n\ge2}(n-1)a_{n,k}\,\delta_{n/x}.
```
Setting $`1-q=e^{-s/x}`$ in the preceding power series gives $`\int e^{-sy}d\nu_x(y)\to d_k^*s^{-(\mu+2)}`$ for each $`s>0`$. The limit is the Laplace transform of the measure with density $`d_k^* y^{\mu+1}/\Gamma(\mu+2)`$ on $`[0,\infty)`$. To see convergence on bounded intervals, tilt $`\nu_x`$ by $`e^{-y}`$ and divide by its mass. The probability transforms converge to $`(1+s)^{-(\mu+2)}`$, the transform of a gamma law. The continuity theorem for Laplace transforms gives weak convergence. Undoing the tilt on any bounded interval with positive endpoints then gives the asserted interval convergence.

Take integer $`x=n`$. Monotonicity of $`a_{m,k}`$ bounds its value at $`n`$ from above by its $`(m-1)`$-weighted average over $`n\le m\le(1+h)n`$, and from below by that over $`(1-h)n\le m\le n`$. Each denominator divided by $`n^2`$ tends to the integral of $`y`$ over the corresponding interval. Use the convergence of $`\nu_n`$ on the same intervals and let $`h\downarrow0`$. This proves (1) with
```math
c_k=\frac{d_k^*}{\Gamma(\mu+2)} .
```
In particular $`c_2>0`$ and $`c_k-c_{k+1}>0`$ by (4). By dominated convergence, $`jG_\mu f_k(a)\downarrow0`$ for each $`a`$, since each residual particle has finite degree at finite time, the family is nonexplosive, and $`R_t W(J_a)\le e^{\rho t}W(J_a)`$. The forcing is uniformly bounded in $`k,a`$. Integrating over $`\ell`$ gives $`c_k\downarrow0`$, so the point masses in (2) sum to one.

## 4. The tail of the limiting law

The probability $`\pi`$ constructed in the [local operator proof](../operator.md) is a weak limit of the tilted root-uniform mean measures. Equation (5) and its $`T_3`$ version show that $`m_t(f_k)/m_t(1)\to d_{f_k}/d_1`$. Along the tilts in that construction the normalised mass of each bounded time interval tends to zero, because $`w(0)/Z\to0`$ and $`w(t)/w(0)`$ is bounded there. Since $`f_k`$ is a bounded continuous function of the rooted local marks, $`\pi(D_o\ge k)=c_k/c_2`$. The positive finite $`A`$ in the [degree-tail argument](../degree-law.md) is the $`\mu`$th moment of the limiting normalised degree of a newborn sampled from $`\pi`$ at a promotion. Applying the point-probability estimate there gives $`p_k\sim\mu A k^{-1-\mu}`$. The large-$`k`$ statement concerns the distribution obtained after the fixed-$`k`$ limit in (2).
