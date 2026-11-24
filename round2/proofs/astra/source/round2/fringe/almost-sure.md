# Exponential convergence in continuous time and almost-sure degree frequencies

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 6: Claude and Gemini; GPT standalone audit; inputs include Round 1 almost-sure exponent and Claims 2, 3 and 5. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Statement submitted for independent review

For the ordinary tree from one edge there is a finite random variable $`X>0`$ almost surely such that
```math
 C_n/n^\mu\longrightarrow X,\qquad
 N_{\ge k}(n)/C_n\longrightarrow c_k/c_2\quad\text{almost surely for every fixed }k\ge2 .
\qquad\text{(1)}
```
Here $`c_k`$ are the constants in [the mean asymptotic](asymptotics.md). In particular the large-$`k`$ asymptotic for the limiting point probabilities in that result also describes the almost-sure limiting core frequencies. This argument uses only the nineteen-mark renewal construction, its strict residual growth bound $`\rho<\mu`$, and the repository's almost-sure logarithmic growth result with $`\mu>1/2`$. Time in Sections 1–4 is continuous, with individual vertex clocks of rate one. The letter $`P_t`$ denotes the *branching first-moment* semigroup of the renewal construction; an ordinary tree expectation is always denoted $`\mathbb E_T`$.

## 1. The finite renewal kernel has an exponential remainder

Use the capped birth marks in [the finite kernel construction](finite-kernel.md). Put $`B_0g(S)=b_S g(m(D_{\rm root}))`$ for a column on the nineteen labels, and write $`J`$ for restriction to newborns. The time kernel is the matrix $`k(t)=JR_t B_0`$. Its entries are continuous, dominated by $`C e^{\rho t}`$, and its Laplace transform on $`\Re z>\rho`$ is $`K(z)`$. For a vector $`g`$ of maximum modulus one, $`\|R B_0g\|_W\le C`$ for an absolute finite $`C`$. Indeed changes of the old particle's score arise only from events at its root. A gain while $`D<20`$ has rate at most 20 and changes that score by at most 2. A gain at $`D\ge20`$ changes $`b`$ by $`c/[D(D+1)]`$ with $`c`$ the core degree, at rate at most $`D`$. A root promotion has rate at most one and bounded increment. The killing term is uniformly bounded. Absolute copy contributions are bounded by $`\sum_Y a_YW(Y)\le(a_S+d)W(S)`$. These facts, $`W\ge1`$ and the finite-time moment estimates give
```math
 k(t)g = JB_0g + \int_0^t JR_u R B_0g\,du .
```
For instance apply the jump formula to the sum of this bounded score in the residual family and stop before large population sizes. Its *nonzero* absolute jump contributions have total expected rate at most a constant times the expected sum of $`W`$, by the very bounds just given; for a high-degree root this uses the gain increment bound. This permits removal of stopping. Integration by parts gives $`K(z)=O(1/|z|)`$ uniformly on each half-plane $`\Re z\ge\rho+\varepsilon`$.

The Perron eigenvalue of $`K(\mu)`$ is simple and $`I-K(z)`$ is invertible on the rest of $`\Re z=\mu`$ by the nonlattice argument of [the mean asymptotic](asymptotics.md), with positive density of input label 2 at all positive times. It is also invertible to the right, by comparison with $`K(\Re z)`$. At large $`|\Im z|`$ the uniform $`O(1/|z|)`$ bound makes its inverse a convergent Neumann series. Compactness of the remaining finite line segment, isolation of the simple pole, and the open inverse condition now give $`\rho<\sigma<\mu`$ such that the only singularity of $`(I-K(z))^{-1}`$ on $`\Re z\ge\sigma`$ is that pole. Choose positive Perron vectors $`g,\ell`$ with $`\ell g=1`$, and put
```math
D=\frac{g\otimes\ell}{-\ell K'(\mu)g}.
```

Let $`v(t)=\sum_{n\ge1}k^{*n}(t)`$, where products in each convolution are ordinary matrix products in time order. On bounded intervals the series converges uniformly by the exponential-series bound for convolution powers of a bounded time kernel. Its Laplace transform for $`\Re z>\mu`$ is $`(I-K(z))^{-1}-I`$. For $`v-k`$, the transform is $`(I-K)^{-1}-I-K=(I-K)^{-1}K^2`$. This is $`O(1/|\Im z|^2)`$ on both the vertical lines through $`\sigma`$ and through any fixed $`s>\mu`$, and uniformly on high horizontal segments between them. On the latter vertical line discounted $`v-k`$ is integrable by the nonnegative generation series. Fourier inversion and a rectangular contour shift, using the integrable high-frequency bound, give
```math
 v(t)=e^{\mu t}D+O(e^{\sigma t}),\qquad t\ge0 . \qquad\text{(2)}
```
The first convolution term $`k(t)=O(e^{\rho t})`$ is absorbed in the remainder.

For any score $`f`$ with $`|f|\le W`$, put $`h_f(t)=JR_t f`$ and $`b(t)=R_t B_0`$. The matrix renewal equation gives $`JP_t f=h_f(t)+(v*h_f)(t)`$. Thus
```math
 JP_t f=e^{\mu t}g\,d(f)+O(e^{\sigma t}),\qquad
 d(f)=\frac{\ell\int_0^\infty e^{-\mu u}h_f(u)\,du}{-\ell K'(\mu)g}.
```
All constants in the error are uniform for $`|f|\le W`$: the tail from replacing the upper integration endpoint $`t`$ in the leading term by infinity is $`O(e^{\rho t})`$, and $`\sigma>\rho`$. For a general planted starting state the first-generation equation is $`P_t f=R_t f+(b*JP_\cdot f)(t)`$. The same convolution estimates yield
```math
 \big|P_t f(S)- e^{\mu t} H(S)d(f)\big|\le C e^{\sigma t}W(S),\qquad
 H(S)=\int_0^\infty e^{-\mu u}(R_u B_0g)(S)\,du . \qquad\text{(3)}
```

## 2. A physical eigenfunction valid for every finite initial tree

Start a finite tree $`T`$ with at least one core, rooted at an initial core $`o`$, with any permanent marks satisfying the parent mark inequality and taking one of the nineteen values. Ordinary trees grown from $`T_3`$ satisfy this. Define $`\mathcal I(T)f`$ as the *unaveraged* immigration expression (2) in the renewal proof, using the rounded birth label. The full tree has a root with no parent slot and its unnormalised weight $`Z(T)`$ includes $`h(c_o)`$. The sum of $`W(S_v)`$ over every noninitial core is at most $`C_0 C(T)`$ by the geometric path factor $`z(1/2)<1`$. Also $`Z(T)\le2 C(T)`$. The promotion bound gives $`\mathbb E_T Z(T_t)\le e^{qt}Z(T)`$. In particular
```math
 I_{T,t}W:=\mathbb E_T\mathcal I(T_t)W
 \le 1/d+e^{qt} Z(T).
```
With general initial fringes the physical immigration identity has the additional term $`\sum_{v\ne o}P_t f(S_v)`$. Apply (3) there and in the immigration integral. It gives
```math
 \left|\mathbb E_T \sum_{v\ne o} f(S_v(t))-e^{\mu t}\Psi(T)d(f)\right|
       \le C_1 C(T)e^{\sigma t},
 \qquad
 \Psi(T)=\sum_{v\ne o}H(S_v)+\int_0^\infty e^{-\mu u} I_{T,u}H\,du . \qquad\text{(4)}
```
The error is uniform in $`T`$ and in $`|f|\le W`$. The founder's contribution to any degree count is at most one. Write $`A_k(T)=N_{\ge k}(T)`$ and $`d_k=d(f_k)>0`$ with $`f_k=\mathbf1_{\{D_{\rm root}\ge k\}}`$ on the finite-mark space. Then (4) applies to $`\mathbb E_T A_k(T_t)`$ including the founder with an increased constant, and also to $`C=A_2`$.

The function $`\Psi`$ is nonnegative and $`\Psi(T)\le C_2 C(T)`$. The pointwise eligible-vertex comparison in the renewal proof with the rounded marks gives $`C(T)\le C_3(1+\Psi(T))`$. Indeed newborn right eigenvector values are bounded below and the two local scripts and counting bound there do not use an unrounded mark. Moreover $`\Psi(T_3)>0`$ by the positive immigration rate from the founder. Equation (4) for $`f=1`$ makes $`\Psi(T)=d_2^{-1}\lim_{u\to\infty}e^{-\mu u}\mathbb E_T C(T_u)`$ independently of the compatible initial orientation or marks. For actual finite elapsed $`h`$, apply the Markov property and (4) at $`T_h`$; its error divided by $`e^{\mu u}`$ has expectation tending to zero, since $`\mathbb E_T C(T_h)<\infty`$. This gives
```math
 \mathbb E_T\Psi(T_h)=e^{\mu h}\Psi(T).
\qquad\text{(5)}
```
Thus $`Z_t^*=e^{-\mu t}\Psi(T_t)`$ is a nonnegative martingale and has an almost-sure finite limit $`Z_\infty^*`$.

## 3. Second moments and positivity

We give the details since the continuous-time vertex population grows at rate one while the core grows at rate $`\mu`$. For any finite tree $`T`$ with a core and any one-leaf extension $`T^+`$ there is a coupling such that
```math
 |\mathbb E_{T^+}A_k(T_u)-\mathbb E_T A_k(T_u)|\le 3 e^{2u}\quad(k\ge2).
\qquad\text{(6)}
```
If the added leaf has a core parent, mark the copies of its parent dirty. If it promotes a leaf, mark both copies of that leaf's old parent and the new core. Other cores are clean and initially paired with matching leaf counts and corresponding core neighbours. Synchronise rings at clean paired sources (including their leaf gains and choice of neighbour). A dirty ring in either tree occurs at rate one per dirty core. It may create one dirty newborn or send a leaf to a clean recipient, marking that pair dirty. A leaf gain at a dirty core changes only its own leaf stock. A clean ring that promotes creates a clean pair, and sends into dirty vertices do not create new marks. Dirty vertices stay dirty. The expected combined number of dirty cores is therefore at most $`3e^{2u}`$ by Grönwall, with standard finite-time localization. All clean contributions to $`A_k`$ cancel.

Write $`n(T)`$ for vertex count. Conditioning to a fixed future time $`b`$ gives the backward martingale $`\mathbb E_{T_v}A_k(T_{b-v})`$ for $`0\le v\le b`$. The total physical rate at $`v`$ is $`n(T_v)`$ and the jump sizes are bounded by (6) with $`u=b-v`$. Martingale isometry, first stopped at finite sizes and then in $`L^2`$ using the finite-time Yule moments, yields
```math
 \mathop{\mathrm{Var}}\nolimits _T(A_k(T_b))
 \le 9\,n(T)\int_0^b e^{4(b-v)+v}\,dv . \qquad\text{(7)}
```
We will also use its one-unit counterpart for the score $`g_h(T)=d_2^{-1}e^{-\mu h}\mathbb E_T C(T_h)`$. The conditional expectation of $`g_h(T_1)`$ at epoch $`0\le v\le1`$ equals $`d_2^{-1}e^{-\mu h}\mathbb E_{T_v} C(T_{1-v+h})`$. Hence
```math
 \mathop{\mathrm{Var}}\nolimits _T(g_h(T_1))\le C_4 e^{4h}n(T).
```
Put $`\eta=\mu-\sigma>0`$. Equation (4) gives $`|\Psi-g_h|\le C_5e^{-\eta h} C`$. Since $`C(T_1)`$ is dominated by a Yule population from $`C(T)`$, $`\mathbb E_T C(T_1)^2\le C_6 C(T)^2`$. The variance inequality for a sum now gives, uniformly in $`T`$ and $`h\ge0`$,
```math
 \mathop{\mathrm{Var}}\nolimits _T\Psi(T_1)\le C_7\big(e^{4h} n(T)+e^{-2\eta h}C(T)^2\big).
\qquad\text{(8)}
```
For integer $`j\ge0`$ choose $`h=\eta^{-1}\log(j+2)`$ in (8). By (5), conditional martingale second moments on $`[j,j+1]`$ equal the previous square plus the conditional variance. Use $`C(T_j)^2\le C_8(1+\Psi(T_j)^2)`$ and $`\mathbb E_T n(T_j)=e^j n(T)`$ to obtain
```math
 \mathbb E_T(Z_{j+1}^*)^2\le (1+C_9(j+2)^{-2})\,\mathbb E_T(Z_j^*)^2+
 C_{10}\big[(j+2)^{4/\eta}e^{-(2\mu-1)j} n(T)+(j+2)^{-2}e^{-2\mu j}\big].
```
The product of the coefficients is finite and the additive series converges since $`2\mu>1`$. Thus
```math
 \mathbb E_T (Z_\infty^*)^2\le C_{11}(\Psi(T)^2+n(T)+1),\qquad
 \mathbb E_T Z_\infty^*=\Psi(T).
\qquad\text{(9)}
```
The integer-time $`L^2`$ bound controls the continuous martingale on all bounded intervals by (5) and conditional Jensen.

From the original $`T_3`$ process the existing almost-sure logarithmic exponent and the Yule size limit imply $`\log C(T_j)/j\to\mu`$ and $`\log n(T_j)/j\to1`$. Hence $`\Psi(T_j)\to\infty`$ and $`n(T_j)/\Psi(T_j)^2\to0`$ almost surely. Conditional on the history at integer time $`j`$, the probability of a positive martingale limit is bounded below by
```math
 \frac{\Psi(T_j)^2}{C_{11}(\Psi(T_j)^2+n(T_j)+1)}
```
by (9) and Cauchy–Schwarz. The conditional probabilities converge to the indicator of that event. The displayed lower limit is positive almost surely, giving $`Z_\infty^*>0`$ almost surely. [Almost-sure logarithmic exponent](../../proofs/almost-sure-exponent.md)

## 4. Core and degree counts on sample paths

Fix a mesh $`1/m`$ and consider times $`t=j/m\to\infty`$. Put $`h=\eta^{-1}\log(t+2)`$, $`u=t-h`$, discarding a finite prefix so that $`u>0`$. Conditional on the tree at $`u`$, the variance bound (7) and $`\mathbb E n(T_u)=n(T_0)e^u`$ give
```math
 \mathbb E\big[e^{-2\mu t}(A_k(T_t)-\mathbb E_{T_u} A_k(T_h))^2\big]
       \le C_{12} n(T_0)e^{-(2\mu-1)t+4h}.
```
It is summable on the mesh. Chebyshev and Borel–Cantelli send this centered term to zero almost surely. By (4) the conditional mean divided by $`e^{\mu t}`$ differs from $`d_k Z_u^*`$ by at most a constant times $`e^{-\eta h}e^{-\mu u}C(T_u)`$. The latter normalized count is bounded almost surely by the reverse comparison to $`1+\Psi`$ and martingale convergence. Thus $`e^{-\mu t} A_k(T_t)\to d_k Z_\infty^*`$ on the mesh. Intersect over all $`k\ge2`$ and all meshes $`1/m`$; monotonicity of $`A_k`$ sandwiches intervening times with factors $`e^{\pm\mu/m}`$. This gives the same convergence at all real times.

Starting the continuous clock from $`T_3`$, let $`\tau_n`$ be the epoch of size $`n`$. The Yule process satisfies $`ne^{-\tau_n}\to Y_\infty>0`$ almost surely, for example by the sum of three independent rate-one Yule limits. Therefore $`A_k(T_{\tau_n})/n^\mu\to d_k Z_\infty^*/Y_\infty^\mu`$ almost surely. Starting at the original edge produces the same discrete chain from size three. Finally $`d_k/d_2=c_k/c_2`$ by the mean limits in continuous time from $`T_3`$ in the two equivalent descriptions (one using rounded birth marks, the other using exact birth marks). This proves the stated implications conditional on the renewal construction and the published exponent result.

## 5. The limiting size has nontrivial support

The possible values of $`X`$ in (1) approach zero and are unbounded above. To see this, use (5) at an arbitrary tree of size $`m`$. Its total jump rate is $`m`$ and its one-event successors are finite in number. Conditioning on the first holding time, or differentiating at zero using $`\Psi(T)\le C_2 C(T)`$ and Yule moments, gives
```math
 \mathbb E[\Psi(T_{m+1})\mid T_m]=(1+\mu/m)\Psi(T_m).
```
Therefore $`Y_m=\Gamma(m)\Psi(T_m)/\Gamma(m+\mu)`$ is a nonnegative discrete martingale. It is bounded in $`L^2`$. Indeed in the continuous-time coupling $`\tau_m`$, the time to size $`m`$ from three, is independent of the entire jump chain, and
```math
 \mathbb E e^{-2\mu\tau_m}
 =\prod_{j=3}^{m-1}\frac{j}{j+2\mu}\asymp m^{-2\mu}.
```
By (9) and optional sampling of the uniformly integrable martingale, $`\mathbb E (Z_{\tau_m}^*)^2\le\mathbb E (Z_\infty^*)^2`$. Equivalently this follows by applying conditional Jensen at bounded stopping times and Fatou. Independence in the displayed identity gives the $`L^2`$ bound for $`Y_m`$. Since $`ne^{-\tau_n}\to Y_\infty^{\rm size}>0`$ and $`\Gamma(n+\mu)/\Gamma(n)\sim n^\mu`$, (1) and the construction in Section 4 give $`X=d_2\lim_m Y_m`$. In particular
```math
 \mathbb E[X\mid T_m]= d_2\,\frac{\Gamma(m)}{\Gamma(m+\mu)}\Psi(T_m).
```
A star of each size $`m`$ has positive probability under the rule and has only one core. The upper comparison $`\Psi\le C_2 C`$ makes its conditional mean tend to zero. A path of each size $`m`$ also has positive probability and has $`m-2`$ cores. The lower comparison makes its conditional mean grow at least as a positive constant times $`m^{1-\mu}`$ for all sufficiently large $`m`$. Thus $`\mathbb P(X<\varepsilon)>0`$ and $`\mathbb P(X>M)>0`$ for every $`\varepsilon,M>0`$. These statements strengthen the nontrivial-support condition in the first question of Section 7 of [Addario-Berry et al.](https://arxiv.org/abs/2403.20185).
