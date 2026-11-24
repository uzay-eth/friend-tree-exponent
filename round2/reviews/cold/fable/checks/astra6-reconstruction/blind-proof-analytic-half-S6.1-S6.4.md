# Independent proof of S6.1–S6.4 (capped construction)

Notation: use entrywise modulus for matrices, the maximum norm for vectors, and the induced maximum row-sum norm for matrices. Write $`n=19`$ in the kernel argument, and $`l_*=113/200, \ u_*=1133/2000`$ for the endpoints of the given bracket (to distinguish them from $`\ell`$). All proofs are from the supplied statements. Constants and their effectiveness are collected in Step 20.

## S6.1 — kernel

**1. GIVEN.** Use the capped generators with the actual, uncapped send rates, the stated nonexplosion and weight estimates, and the given fact $`r(K(\mu))=1`$, with $`l_*\le\mu\le u_*`$. In particular $`1<W`$, $`W(J_{1/i})\le 1/d`$, $`\|R_t\|_W\le e^{\rho t}`$, and $`\|P_t\|_W\le e^{(\rho+1/d)t}`$. All particle expectations below apply also to signed or complex scores by absolute domination.

**2. PROVED (absolute jump estimate).** Let $`F=B_0 g_0`$ with $`\|g_0\|_\infty\le1`$; then $`|F|\le1`$. Consider the score summed over the residual population. For one type with root degree $`D`$, only a root-degree increase or root promotion can change its intrinsic score. The rate of root-degree increases is $`L+\sum_Y 1/D_Y\le D`$. When the capped label stays the same, the jump of $`F`$ under a degree increase is at most $`1/(D+1)`$ since $`(L, D)\mapsto(L+1, D+1)`$. When the label changes, $`D\le19`$ and the jump is at most 2. Root promotions have total rate $`b\le1`$ and jump size at most 1. Death contributes at most $`1/2`$. All copy events together contribute at most $`\sum_Y a_Y`$, and
```math
\sum_Y a_Y\le (1/2+d)\sum_Y z(a_Y)
\le \frac{(1/2+d)^2}{7/10}\,W(S),
```
by retaining the immediate children in the weight. Consequently the **absolute** score-jump intensity per particle is at most $`42W(S)`$, and so is $`|RF(S)|`$. The potentially unbounded leaf-ring rate does not cause a problem even at enormous degree.

**3. STANDARD (jump compensation), PROVED (regularity and decay).** The usual rate compensation identity for a nonexplosive pure jump process applies: stop at finitely many events, use the conditional event rates, and remove the stops by integrability of the absolute jump sum. Here integrability follows from Step 2 and $`\int_0^t42R_sW(S)\,ds<\infty`$; the identity for the nonnegative absolute jump sum follows from monotone convergence. It gives $`R_tF=F+\int_0^t R_sRF\,ds`$ pointwise. Restriction gives the claimed identity for $`k`$. In particular $`k`$ is locally absolutely continuous, with $`\|k(t)\|\le e^{\rho t}/d`$ and $`\|k'(t)\|\le42 e^{\rho t}/d`$ a.e. (use the 19 coordinate tests on a common full-measure set and then linearity and the bound for arbitrary $`g_0`$). For $`\Re z\ge\rho+\varepsilon`$, the discounted derivative is integrable and $`e^{-zt}k(t)\to0`$. Integration by parts gives
```math
K(z)=\frac{k(0)+\int_0^\infty e^{-zt}k'(t)\,dt}{z}, \qquad
\|K(z)\|\le \frac{1+42/\varepsilon}{d|z|}.
```
Domination by $`e^{-\Re z\,t}e^{\rho t}t^j/d`$ proves analyticity and the derivative formulas in $`\Re z>\rho`$.

**4. PROVED (positivity and nonlattice).** From any newborn, let the original residual particle (no death) take only root leaf-ring growth steps until root degree $`j`$, then hold that type before evaluation at time $`t>0`$. The simplex of such jump times has strictly positive density, or for $`j=2`$ use the no-jump event. Hence every $`k_{ij}(t)>0`$ for every $`t>0`$. Integrating this path against $`e^{-st}`$ gives $`K_{ij}(s)\ge E:=e_*`$ for $`l_*\le s\le u_*`$: each no-other-event denominator is $`s+1/i+D-1/D`$, each growth numerator $`D-1`$, and the final score $`(j-1)/j`$. This is the product in S5.2 with denominators upper bounded using $`u_*+1/2`$. In addition $`k_{ii}(t)\ge\tfrac12 e^{-(1/i+3/2)t}`$ by no event. At $`z=\mu+i\tau`$, $`\tau\ne0`$, integration of this absolutely continuous component gives strict cancellation in each diagonal; more explicitly, with $`A_i=\mu+1/i+3/2`$,
```math
|K_{ii}(\mu+i\tau)|\le K_{ii}(\mu)-\tfrac12\left(\frac1{A_i}-\frac1{\sqrt{A_i^2+\tau^2}}\right).
```
All other entries satisfy $`|K_{ij}(\mu+i\tau)|\le K_{ij}(\mu)`$. Thus this argument excludes lattice effects using a positive density over an interval, rather than using only positivity of the time-integrated matrix.

**5. STANDARD (Perron–Frobenius), PROVED (simple pole).** The strictly positive finite matrix $`K(\mu)`$ has a simple eigenvalue 1 with positive right and left eigenvectors, and all others have modulus below 1. Normalize $`\max_i g_i=1`$, $`\ell g=1`$; since $`K_{ij}(\mu)\ge E`$, $`g_i\ge E`$, $`\ell\mathbf1\le1/E`$. In the weighted maximum norm $`\|x\|_g=\max|x_i|/g_i`$, Step 4 gives $`\|K(\mu+i\tau)\|_g<1`$ for $`\tau\ne0`$. Also $`\|K(x+i\tau)\|_g<1`$ for $`x>\mu`$ since even the diagonal real integrals strictly decrease. Let $`\alpha=-\ell K'(\mu)g`$; the diagonal density in Step 4 gives
```math
\alpha\ge \beta=\frac1{2(u_*+2)^2}>0.
```
Decompose $`\mathbb C^n=\mathbb Cg\oplus\ker\ell`$ with projection $`g\otimes\ell`$. The off-diagonal blocks of $`I-K(\mu)`$ vanish, its complementary block is invertible, and the scalar Schur complement for $`I-K(z)`$ is $`\alpha(z-\mu)+O((z-\mu)^2)`$. This is because its upper-left coordinate is $`\ell(I-K(z))g`$ and the product of the off-diagonal blocks is of quadratic order. Block inversion gives the simple pole and the stated residue $`D^{\rm ker}`$.

**6. PROVED (effective strip).** Here are explicit conservative constants; the bounds are deliberately not optimized. Put
```math
s_0=(\rho+l_*)/2, \qquad M_j=\frac{j!}{d(s_0-\rho)^{j+1}}\quad(j=0,1,2).
```
These bound $`\|K^{(j)}\|`$ on $`\Re z\ge s_0`$. The matrix $`\mathrm{diag}(g)^{-1}K(\mu)\mathrm{diag}(g)`$ is stochastic with entries at least $`E^2`$. Subtract the matrix with every entry $`E^2`$ and apply the induced operator on the quotient by the constant vectors (with the quotient of the maximum norm). The subtracted matrix acts by zero there and the remaining row sums are $`1-nE^2`$. Thus each non-Perron eigenvalue has modulus at most $`1-nE^2`$. If $`\Delta(z)=\det(I-K(z))`$, the Schur complement from Step 5 gives $`|\Delta'(\mu)|\ge L_0=\beta(nE^2)^{n-1}`$. By permutation expansion with the product rule, on the same half-plane,
```math
|\Delta''(z)|\le F_2=n! \left[n M_2(1+M_0)^{n-1}+n(n-1)M_1^2(1+M_0)^{n-2}\right].
```
Take $`r=\min\{l_*-s_0,1, L_0/F_2\}`$. Taylor's integral remainder along the line segment gives $`|\Delta(z)|\ge L_0|z-\mu|/2`$ when $`|z-\mu|\le r`$. For $`|\tau|\ge r/2`$, Step 4 gives a diagonal deficit of at least
```math
\delta_0=\frac1{2(u_*+2)}\left(1-\frac{u_*+2}{\sqrt{(u_*+2)^2+(r/2)^2}}\right)>0
```
in each row; here each factor $`1/A_i`$ and $`1-A_i/\sqrt{A_i^2+\tau^2}`$ is bounded below by its value at $`A_i=u_*+2, \ |\tau|=r/2`$. Thus $`\|K(\mu+i\tau)\|_g\le1-\delta_0`$. Put
```math
\eta=\min\{r/4, E\delta_0/(2M_1)\}, \qquad \sigma=\mu-\eta.
```
Within $`\sigma\le x\le\mu`$ and $`|\tau|\ge r/2`$, shifting in real part changes the matrix in the $`g`$-norm by at most $`\eta M_1/E\le\delta_0/2`$. When $`|\tau|<r/2`$ in this horizontal range, use the disc bound instead. Combined with Step 5 this proves the whole half-plane assertion. In particular the strip width is explicit in the supplied numerical bracket and constants, without needing a numerical search for complex zeros. If a rational line is desired, any rational line in $`[\mu-\eta, \mu-\eta/2]`$ works; S5.2's decaying enclosures and strict slope bracket $`\mu`$ to arbitrary accuracy by testing rational points strictly to either side.

**7. PROVED (transform and contour), STANDARD (Fourier inversion and residue theorem).** The bound on $`k`$ gives $`\|k^{*j}(t)\|\le d^{-j}e^{\rho t}t^{j-1}/(j-1)!`$. The series defining $`v`$ converges locally uniformly and $`v-k`$ vanishes at 0. On a real line $`c>\rho+1/d`$, absolute Laplace integrability permits interchanging the sum and integrals, and gives the transform of $`v-k`$ as
```math
F(z):=K(z)^2(I-K(z))^{-1}=(I-K(z))^{-1}-I-K(z).
```
Use $`c=\rho+2/d+1`$. The transform of $`v-k`$, unlike the full renewal transform, is $`O(|z|^{-2})`$ uniformly for $`\sigma\le\Re z\le c`$, by Step 3 and the Neumann series at large $`|z|`$. The discounted function $`e^{-ct}(v-k)(t)`$, extended as zero for negative times, is continuous and integrable and its Fourier transform is integrable. The Fourier inversion theorem therefore gives its inverse everywhere, including at 0. Apply the residue theorem to $`e^{zt}F(z)`$ on rectangles with vertical lines $`c, \sigma`$ and horizontal heights $`\pm R`$. For any fixed $`t\ge0`$ the horizontal integrals vanish as $`R\to\infty`$, bounded by $`(c-\sigma)e^{ct} O(R^{-2})`$. The residue is $`e^{\mu t}D^{\rm ker}`$. On $`\Re z=\sigma`$ the integral of the norm is finite (no poles on the line, high-frequency decay), and the shifted integral is bounded by $`e^{\sigma t}(2\pi)^{-1}\int_{\mathbb R}\|F(\sigma+i\tau)\|\,d\tau`$. Adding $`k(t)`$ proves the assertion for **all** $`t\ge0`$. No smoothness of $`k`$ beyond the local absolute continuity and exponential bound on its a.e. derivative proved above is needed; inversion is on the twice-decaying remainder, not on a conditionally integrable transform.

## S6.2 — branching first moments

**8. PROVED (Duhamel and renewal).** Separate particles whose ancestral line has no direct birth from the descendants of the first direct birth on each other ancestral line. In the full process the former particles have exactly the residual law; a direct birth does not change its source. Rate compensation for the direct births off this residual population and the independent evolution of a newborn give $`P_tf=R_tf+(b*JP_\cdot f)(t)`$. For $`f\ge0`$ the decomposition is justified by summing nonnegative contributions; for $`|f|\le W`$, apply it to positive and negative parts (or real and imaginary parts), with finite $`P_tW`$. Restriction gives $`x=h_f+k*x`$ with $`x(t)=JP_tf`$. Iterate; the remainder $`k^{*m}*x`$ tends to zero on every bounded time interval by the factorial convolution bound in Step 7 and local boundedness of $`x`$. This yields $`x=h_f+v*h_f`$ pointwise, including time zero.

**9. PROVED (uniform asymptotic).** Write $`A`$ for a bound on $`e^{-\sigma t}\|v(t)-e^{\mu t}D^{\rm ker}\|`$ from Step 7. Use $`\|h_f(t)\|\le e^{\rho t}/d`$, uniformly in $`f`$. The convolution of the error in $`v`$ with $`h_f`$ is at most $`A e^{\sigma t}/[d(\sigma-\rho)]`$. Extending the main convolution to infinity after factoring out $`e^{\mu t}`$ gives $`e^{\mu t}g d(f)`$; the added tail is at most $`\|D^{\rm ker}\|e^{\rho t}/[d(\mu-\rho)]`$. Including $`h_f`$, one may take the first error constant
```math
A_J=1/d + A/[d(\sigma-\rho)]+\|D^{\rm ker}\|/[d(\mu-\rho)].
```
We have $`|d(f)|\le M_d:=1/[d E\beta(\mu-\rho)]`$, and $`0\le H(S)\le M_H W(S)`$ with $`M_H=1/(\mu-\rho)`$ since $`\|g\|_\infty=1`$. In Duhamel, $`\|b(t)\|_{\infty\to W}\le e^{\rho t}`$, so convolution of the newborn error costs at most $`A_Je^{\sigma t}W(S)/(\sigma-\rho)`$, and the added main tail costs at most $`M_d e^{\rho t}W(S)/(\mu-\rho)`$. The intrinsic term costs at most $`e^{\rho t}W(S)`$. Take
```math
A_P=1+A_J/(\sigma-\rho)+M_d/(\mu-\rho).
```
These estimates do not restrict the mark, size, root degree, or score beyond the stated hypotheses.

## S6.3 — physical means

**10. PROVED (summing fringe weights).** Every link along a core path has $`z(a)\le\theta=z(1/2)=10/13<1`$, including at the cap. A fixed core $`x`$ contributes to its own fringe and that of each eligible non-founder ancestor, at most one at each upward distance $`j`$. The contributions sum to at most $`(2/d)\sum_{j\ge0}\theta^j`$. Hence $`\sum_{v\ne o}W(S_v)\le C_0 C(T)`$ with $`C_0=2/[d(1-\theta)]`$ (core indices throughout). In $`Z`$ every coefficient is at most one and $`h\le2`$, giving $`Z\le2C`$. A bare path is covered by the geometric factor even where every mark is $`1/2`$; bounding the number of ancestor/descendant pairs without that factor would fail. A star has just the founder core and $`Z=2`$, for any number of leaves; leaves are never counted in these sums. The depth estimate is also unaffected by large degree capped marks: $`z(1/20)=1/4`$, still below $`\theta`$.

**11. PROVED (full-tree drift).** A leaf addition to an existing core changes neither old core marks nor old core degrees *within the core tree*. A promotion by vertex $`x`$ changes $`Z`$ by
```math
 w_x\big[z(1/m(D_x))+h(c_x+1)-h(c_x)\big],
```
since the added core has $`c=1`$, and old descendant path products are unchanged. Its aggregate rate is $`(D_x-c_x)/D_x`$. For $`c_x=0`$ the bracket is negative. For $`1\le c\le4`$, at $`2\le D\le20, \ D\ge c`$, the corresponding maximal rate times bracket is as follows (for $`c\ge4`$, use the row $`c=4`$ since the increment of $`h`$ is zero and the factor decreases with $`c`$). The exact finite maxima may be checked by substituting each integer $`D=\max(2, c), \ldots,20`$ into $`(D-c)/D\,[20/(20+3D)+h(c+1)-h(c)]`$.

| $`c`$ | finite-domain maximum (at $`D`$) | bound for $`D>20`$ | $`q h(c)`$ |
|---|---|---|---|
| 1 | $`313/870`$ (3) | $`1/10`$ | $`9/25`$ |
| 2 | $`289/950`$ (6) | $`9/50`$ | $`153/500`$ |
| 3 | $`103/440`$ (8) | $`17/100`$ | $`351/1250`$ |
| $`\ge4`$ | $`140/583`$ (11) | $`1/4`$ | $`63/250`$ |

The bounds for $`D>20`$ follow from $`z=1/4`$ and the factor at most one. For a direct check of the finite maxima without enumerating: put $`s=h(c)-h(c+1)`$; the derivative in real $`D`$ has the sign of $`20(20c+6cD-3D^2)-cs(20+3D)^2`$, whose leading coefficient is negative and whose constant coefficient is positive for these four rows. There is a single positive derivative zero, between 3 and 4, 5 and 6, 7 and 8, and 10 and 11 respectively, by plugging in those endpoints. Comparing the two adjacent integers gives the entries. For $`c`$ larger than 20 there is only the capped domain. Summing gives the physical generator bound $`\mathcal L Z\le qZ`$. Jump compensation is integrable on bounded time intervals: each increment is bounded absolutely by 2, the promotion rate is at most the current size, and size is a rate-one-per-vertex Yule process. This process has finite moments of every fixed order at fixed times, e.g. the size started from one is geometric with parameter $`e^{-t}`$ by the elementary birth equations, and an initial size $`m`$ gives a sum of $`m`$ such variables. Gronwall yields the asserted $`\mathbb E_T Z(T_t)\le e^{qt}Z(T)`$.

**12. PROVED (founder influx).** Inserting a leaf at the root of a planted tree does not change any old mark, core count or $`c_x`$; $`W(Y^+)=W(Y)`$. Thus the founder copy contribution to $`mathcal I(T)W`$ (where $`\mathcal I`$ is the operator in the question) is exactly $`\sum_Y a_Y W(Y)=Z(T)-h(c_o)`$, by factorizing outer path weights at each core child. Direct immigration costs at most $`1/d`$. This proves $`I_{T, t}W\le 1/d+e^{qt}Z(T)\le A_I C(T)e^{qt}`$, with $`A_I=1/d+2`$ since $`C(T)\ge1`$.

**13. PROVED (arbitrary-start first-moment identity).** For a score $`F`$, put $`G_F(T)=\sum_{v\ne o\;\mathrm{core}}F(S_v)`$. At a physical event, record its changes in each old fringe plus any newly created fringe. Each old fringe's changes with both endpoints internal give $`L^0F`$. Its only outside change is receipt of a leaf from the parent at rate $`p_{\rm parent}`$, costing $`F(S_v^+)-F(S_v)`$. New fringes start at $`J`$ with the rounded parent label, at each parent's promotion rate. Summing $`(R+B)F`$ over the old non-founder fringes and adding $`\mathcal I(T)F`$ gives exactly this list: for each old fringe the $`-a_vF(S_v)`$ cancels the $`+a_vF(S_v)`$ from the parent copy terms after writing them as $`a_vF(S_v)+p_{\rm parent}(F(S_v^+)-F(S_v))`$. The founder is included only through $`\mathcal I`$, not in $`G_F`$. Therefore
```math
\mathcal L G_F=G_{(R+B)F}+\mathcal I(T)F.
```
Use this equality for $`F=P_{t-s}f`$, $`0\le s\le t`$, in the time-dependent jump formula. Here are regularity and integrability details. Conditioning on the first event of one branching particle gives, for each type, a holding term and a convolution against an exponential density with a finite linear combination of successor first moments. (At a birth use the sum of source and new-particle contributions.) The time functions are measurable and locally bounded by the given semigroup bounds. The exponential convolution against a locally bounded function is continuous, so each first moment is continuous; putting this back into the first-event equation gives the pointwise backward derivative $`\partial_u P_uf=(R+B)P_uf`$, with right derivative at zero. For physical size $`m`$, all old and one-jump successor fringes have size at most $`m+1`$, hence $`|P_{t-s}f|\le e^{(\rho+1/d)t}(m+1)/d`$ there. The per-particle total intrinsic/death/emission rates are $`O(m+1)`$ with an explicit bound $`3(m+1)`$. Thus the sum of time derivatives over fringes and the physical absolute jump intensities are at most a fixed polynomial in $`m+1`$ times $`e^{(\rho+1/d)t}/d`$; specifically $`12(m+1)^3 e^{(\rho+1/d)t}/d`$ bounds each. The Yule moment bound in Step 11 justifies interchanging the time integral and expectation and using jump compensation. The derivative cancels the $`G_{(R+B)F}`$ term, leaving
```math
\mathbb E_T G_f(T_t)=\sum_{v\ne o\;\mathrm{core}}P_tf(S_v)+\int_0^t I_{T, s}P_{t-s}f\,ds.
```
This proves the stated extension, without requiring the supplied identity to apply to any seed other than $`T_3`$.

**14. PROVED (uniform estimate after immigration).** Sum Step 9's errors over initial non-founder fringes using Step 10. In the integral put $`P_{t-s}f=e^{\mu(t-s)}H d(f)+E_{t-s}`$. Step 12 bounds the integrated error by $`A_P A_I C(T)e^{\sigma t}/(\sigma-q)`$. Extending the main term to $`s=\infty`$ costs at most $`A_I M_H M_d C(T)e^{qt}/(\mu-q)`$. This also proves convergence of the integral defining $`\Psi`$. Thus one constant for all scores and trees is
```math
A_G=A_P C_0+\frac{A_PA_I}{\sigma-q}+\frac{A_IM_HM_d}{\mu-q}.
```

**15. PROVED (actual degrees, all thresholds).** For $`k\ge2`$, let $`f_k(S)=\mathbf1_{\{D_{\rm root}\ge k\}}`$, using the true root degree including the parent slot, and not $`m(D_{\rm root})`$. This is bounded by $`W`$. A newborn original residual particle can reach degree $`k`$ by $`k-2`$ consecutive leaf-ring growths and can then hold for any positive interval, with strictly positive probability, so $`d_k>0`$. Quantitatively the integral from any newborn is at least $`L_k=(u_*+1/2+k-1/k)^{-1}\prod_{j=2}^{k-1}(j-1)/(u_*+1/2+j-1/j)`$ by the same path integral as Step 4, now with score 1. Since $`\alpha\le (\ell\mathbf1)M_1`$, this gives $`d_k\ge L_k/M_1`$. The founder contributes at most 1 to $`A_k`$, including $`k=2`$, so $`C_1=A_G+1`$ works simultaneously for all $`k`$ and also for Step 14. Both the founder's direct births and the label of the kernel are capped; this affects neither the actual physical degree nor the indicator just used.

## S6.4 — eigenfunction

**16. PROVED (upper comparison).** Positivity is immediate from the definition. Steps 9, 10 and 12 give $`\Psi(T)\le C_2 C(T)`$ with $`C_2=M_H[C_0+A_I/(\mu-q)]`$.

**17. PROVED (lower comparison, including bare paths).** A type whose root has $`c\le4`$ and at least one leaf has $`H\ge h_L=E/[5(u_*+3/2)]`$: retain only the original particle, and at time $`t`$ use the event that it has neither died nor had its root ring, of probability at least $`e^{-3t/2}`$. On this event $`c`$ is fixed, its leaves can only increase, and $`b\ge1/5`$. If instead that root is bare (no leaves), with $`c\le4`$, and it has a bare core child also with $`c\le4`$, then $`H\ge h_B=h_L/[4(u_*+1/2+4)]`$. Indeed, in the intrinsic process up to the first change of either vertex's degree or leaf count, the intended child-to-root send has constant rate $`1/c_{\rm child}\ge1/4`$. All changes at the pair while bare are incoming core sends at aggregate rate at most 4: together they have at most eight core neighbours counted with multiplicity, each send at rate at most $`1/2`$, and the external parent slot sends nothing. Discard emitted particles, and censor the original particle on death (rate at most $`1/2`$). The survival probability without a local change or death to time $`t`$ is at least $`e^{-(4+1/2)t}`$; this follows by rate compensation and the differential inequality for the survival probability, even with unrelated intrinsic events in between. The discounted probability of the designated send first is therefore at least $`1/[4(u_*+1/2+4)]`$, by integrating its rate against that survival probability. After the send use the first lower bound by the strong Markov property of the intrinsic pure jump process and independent death (the first local change is a stopping time and the exponential lifetime is memoryless). This argument does not need a lower bound on the send rate of a neighbour of huge degree.

The core-induced tree is a tree on $`C`$ vertices: on any path between two cores all intermediate vertices have degree at least two. Thus $`\sum c=2(C-1)`$. At most $`2(C-1)/5`$ cores have $`c\ge5`$; call these high. From the set of non-founder, non-high cores, each leaf-bearing core is already one of the first good types. Each bare core in that set has at least one core child (a bare non-founder has $`c=D\ge2`$); choose one child arbitrarily. Chosen children are distinct by uniqueness of the parent. Discard the parents whose chosen child is high, losing at most the number of high cores. If the chosen child is non-high and has leaves, charge its parent to that child's fringe; if the child is non-high and bare, charge to the parent's own fringe by the second bound. The remaining starting vertices (including those with leaves, charged to themselves) number at least $`C-1-4(C-1)/5=(C-1)/5`$. Each charged fringe is non-founder and can be used at most twice, once as its own fringe and once by its unique parent. Hence $`\sum_{v\ne o}H(S_v)\ge h_B (C-1)/10`$. In particular $`C_3=\max(1,10/h_B)`$ works. In $`T_3`$, there are only the founder and two leaves initially. Up to its first ring, the leaves can only increase, so $`b_o=1`$ and $`H(J_{1/m(D_o)})\ge h_L`$. The first ring is exponential of rate 1 independently of additions to the founder. Consequently $`\Psi(T_3)\ge h_L/(\mu+1)>0`$.

**18. PROVED (intrinsic definition and harmonicity).** Step 15 for $`k=2`$, after discount and taking the limit, gives $`d_2\Psi(T)=\lim_{u\to\infty}e^{-\mu u}\mathbb E_TC(T_u)`$. The degree-count dynamics on the right use no orientation and no marks; $`d_2>0`$ is a fixed scalar of the capped construction. This proves independence for every compatible starting orientation and marks. For any $`h\ge0`$, Step 15 bounds $`e^{-\mu u}\mathbb E_{T_h}C(T_u)`$ by $`(d_2 C_2+C_1)C(T_h)`$ uniformly for $`u\ge0`$; the latter is integrable by the Yule size bound. By dominated convergence and the physical Markov property,
```math
d_2\mathbb E_T\Psi(T_h)=\lim_{u\to\infty}e^{-\mu u}\mathbb E_TC(T_{h+u})=d_2 e^{\mu h}\Psi(T).
```

**19. STANDARD (martingale convergence).** Start in the stated domain $`C(T)\ge1`$ (in particular from $`T_3`$). The process never leaves this domain. Harmonicity and the Markov property give the martingale identity conditional on the physical history; the random variables are integrable by Step 16 and the finite expected Yule size. The process is nonnegative and càdlàg: there are finitely many tree jumps on compact time intervals and the time factor is continuous. The continuous-time nonnegative martingale convergence theorem applies (bounded expectations, right-continuous paths), giving an almost-sure finite limit. Strict positivity of that random limit is not part of S6.4 and has not been inferred here.

**20. PROVED (effectiveness of constants).** Bounds in Steps 2–7 are explicit. For clarity a completely bounded version of $`A`$ in Step 9 is as follows. Let $`T_0=(1+42/(s_0-\rho))/d`$, $`R_0=1+2T_0`$. For $`|\tau|\ge R_0`$ on the line $`\sigma+i\tau`$, $`\|F\|\le 2T_0^2/\tau^2`$. For $`|\tau|\le R_0`$, the inverse norm is at most
```math
U_0=\max\left\{\frac{2}{E\delta_0}, \frac{2n!(1+M_0)^{n-1}}{L_0\eta}\right\}.
```
The first bound is the geometric series in the weighted norm for $`|\tau|\ge r/2`$, converted using $`g_i\ge E`$; the second is the adjugate bound in the disc with distance at least $`\eta`$ from $`\mu`$. Thus one can take $`A=1/d+[2R_0 M_0^2U_0+4T_0^2/R_0]/(2\pi)`$. All subsequent error/comparison constants are displayed explicitly. Their occurrences of $`\mu`$ or $`\sigma`$ can use the certified bounds (e.g. $`\sigma-\rho\ge s_0-\rho`$), and $`\|D^{\rm ker}\|\le1/(E\beta)`$. Threshold coefficients $`d_k`$ and the functions $`H, \Psi`$ are defined by integrals, not by closed-form evaluations; Steps 9, 15–17 give explicit bounds, including strict positivity where required. $`g, \ell`$ are determined up to the stated normalizations by the finite positive matrix. The spectral argument needs no unnamed compactness constants.

## Adversarial simulations

**21. PROVED (algorithm law) and numerical evidence.** `tests.py` runs the physical process without simulating every leaf birth. A clock of rate one runs per core, from when it becomes a core. Between times at which a core's leaf count is queried or changed by a core event, its leaves independently generate additional leaves at total rate equal to their number; with $`L`$ leaves, the number of additional leaves after duration $`s`$ has the negative-binomial law for $`L`$ successes and success probability $`e^{-s}`$, from the sum of $`L`$ geometric Yule families. Core clocks and core-neighbour lists can be maintained independently of these latent draws until a ring needs a leaf count; a send to an existing core first updates its latent leaf count, then adds the new leaf. Promotion removes one leaf of the source and creates a new core with one leaf. Poisson memorylessness supplies the promoted vertex's new independent core clock. Hence this sampling is exact in law (up to numerical RNG arithmetic) and works even when the underlying number of leaves exceeds billions. The starting caterpillar is a 101-core path with $`1+(i\bmod 5)^2`$ leaves on backbone vertex $`i`$ (0-based); the broom has a 101-core path with $`100{,}000`$ leaves at one end, one at the other, and none along the handle. The bare paths have one leaf at each end and none inside. Stars use each power of ten from $`10`$ through $`10^6`$.

Command (in this directory): `taskset -c 12 nice -n 5 round2/.venv/bin/python -B tests.py --log tests.log`. `tests.log` gives the complete table and paired standard errors. An independent literal per-vertex-ring simulation on a three-core bare path (2,048 runs per implementation) had discrepancies of $`-0.39,-0.20,-0.56`$ standard errors at times 1, 2, 3. Over the stress family the estimates at time 10 range from about 0.50 for the longest bare path and 0.50 for the broom to 1.46–1.51 for stars; the caterpillar is about 1.17. All paired differences from time 8 to 10 are under two paired Monte Carlo standard errors. The longest bare paths and broom show a visible early transient, but increasing initial size or maximum degree shows no growing late-time drift in these tests. Different limiting heights are allowed. The finite sample and finite time grid cannot certify the uniform bound or its exponent.

## Statements proved

S6.1 in Steps 1–7; S6.2 in Steps 8–9; S6.3 in Steps 10–15; S6.4 in Steps 16–19. Constants in Step 20. The adversarial tests are separate evidence, with no role in the proofs.

## Statements repaired

No substantive repair required. In S6.3 the use of $`f_k`$ is made explicit as the indicator of the true planted root degree, and the founder must be restored separately. The eigenfunction and martingale statements are used on their stated domain $`C(T)\ge1`$.

## Gaps

None found in S6.1–S6.4, conditional on the allowed S2/S5 facts. The umbrella theorem S6 (including positivity of the martingale limit and the discrete-time limits) was outside this assignment.

## Where the proof would break

- Without the exponential weighted moment estimate for the residual process, neither the jump-derivative bound in time nor the uniform renewal tails would follow by this proof.
- Without a positive diagonal time density, positivity of the real matrix alone would not exclude nonreal poles on the Perron line. Without subtracting $`k`$, the integration-by-parts estimate alone does not give absolute vertical integrability.
- Without a uniform upper bound $`z(a)<1`$, a bare path could make the sum of ancestor-fringe contributions grow superlinearly. Without the actual send rate in the copy cancellation, the arbitrary-start physical generator identity would change.
- The lower comparison does not assume there are many leaves initially or that a bare core has a neighbour sending at a uniformly high rate without further information. It uses the core-tree degree count and tests the local send only between two bounded-degree bare cores.
