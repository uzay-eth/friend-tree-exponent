# A finite matrix of explicit series and rational enclosures

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 5: Gemini; Claude checked weights and the determinant identity; inputs include Claim 2 and the Round 1 bracket. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Construction submitted for independent review

Let $`m(D)=\min(D,20)`$. There are nineteen possible newborn labels $`2,\ldots,20`$. Define below a $`19`$-by-$`19`$ matrix $`K(s)`$ without reference to $`\mu`$ or to any limiting distribution of the original tree. Its entries are sums of explicit nonnegative rational functions of $`s`$. On the half-line $`s>51/100`$, $`\mu`$ is the largest real zero of
```math
\det(I_{19}-K(s))=0. \qquad\text{(1)}
```
It is the unique solution there of $`r(K(s))=1`$. A pair of finite rational calculations bounds this radius from below and above. Their convergence estimate is given in (7). The determinant has fixed dimension but its entries use arbitrarily large finite trees. An evaluation of (1) to prescribed numerical accuracy requires controlling those entries.

## 1. Nineteen newborn labels

Repeat the [positive renewal construction](renewal.md) with birth mark $`a_v=1/m(D_{\rm parent})`$, assigned when $`v`$ becomes a core. The current inverse parent degree is $`p=1/D_{\rm parent}`$, including at later times, so $`p\le a_v`$. Keep the copy rates $`p_S`$ and $`a_Y-p_S`$ in (1) of that proof unchanged in form. In a direct birth use $`J_{1/m(D_{\rm root})}`$ instead of $`J_{p_S}`$. Both the physical signed-to-positive identity and initial-centre immigration still hold. The rates $`p_S`$ in the signed difference $`p_S[f(Y^+)-f(Y)]`$ remain the actual reciprocal degrees.

For the same $`d=3/20`$, $`q=9/25`$, $`\rho=q+d`$ and $`h`$ as in the renewal proof, the child path factor at promotion is now $`1/(1+d\,m(D))`$. On $`D<20`$ the finite check is unchanged. On $`D\ge20`$ the factor is $`1/4`$, and $`1/4\le q(7/10)=63/250`$. Monotonicity of $`h`$ gives the same intrinsic bound $`L^0W\le qW`$, the same copy-weight identity, and $`RW\le\rho W`$. These estimates include the actual tree at its initial centre.

For $`s>\rho`$ count the expected discounted direct offspring, by their labels, in the residual family of a newborn with label $`i`$. Call the resulting matrix $`K(s)`$. The argument in Sections 2–4 of the renewal construction applies with this finite label set. Each matrix entry is positive since the original newborn can gain leaves until the requested output label is available, remain alive, and produce a direct offspring. All radii above one give divergent generations, and radii below one give finite total discounted occupation. Thus $`r(K(s))=1`$ precisely at $`\mu`$. Above $`\mu`$, all eigenvalues have modulus below one. At $`\mu`$ the matrix has eigenvalue one, giving the largest-real-zero specification (1).

## 2. A formula for the entries without an unknown function

A state $`S`$ is a finite planted marked fringe. Its number of vertices $`n(S)`$ counts its core and leaf vertices, omitting the parent slot. Write $`r^0(S,S')`$ for an intrinsic transition rate, with parallel transitions added. It can be computed directly by listing the finitely many internal arrows in the marked tree. An arrow from a leaf to its parent gives that parent another leaf at rate one. A core-to-core internal arrow $`x\to y`$ gives $`y`$ one leaf at rate $`1/D_x`$. A core with $`L`$ leaf children promotes one of them at total rate $`L/D`$, creating the appropriate mark. Each of these changes increases $`n(S)`$ by one. Write $`r_0(S)=\sum_{S'}r^0(S,S')=n(S)-1/D_{\rm root}`$, since only the root-to-parent slot is inert. Put
```math
t(S,S')=r^0(S,S')+
 \sum_{Y\ {\rm core\ child}}
 \big[p_S\mathbf1_{\{Y^+=S'\}}+(a_Y-p_S)\mathbf1_{\{Y=S'\}}\big].
```
The finite path formula is
```math
K_{ij}(s)=
\sum_{\ell=0}^{\infty}
\quad\sum_{\substack{S_0=J_{1/i},\,S_1,\ldots,S_\ell\\m(D_{\rm root}(S_\ell))=j}}
\left(\prod_{u=0}^{\ell-1}
 \frac{t(S_u,S_{u+1})}{s+a_{S_u}+r_0(S_u)}\right)
 \frac{b_{S_\ell}}{s+a_{S_\ell}+r_0(S_\ell)}. \qquad\text{(2)}
```
Only positive transitions occur in the path sum, hence each inner sum is finite. One way to verify it is to separate the original particle's next intrinsic event or death from all the independent copies it emits while waiting. The expected discounted number of such copies in state $`Y`$ equals its production rate divided by $`s+a_S+r_0(S)`$. The discounted intrinsic successor weights and direct-production reward have the same denominator. Iterating by labelled lineages gives (2) by monotone convergence. More formally, every direct production in a bounded time interval has finitely many intrinsic transitions and copy births along its unique ancestry by nonexplosion. Truncating that number of transitions exhausts the actual discounted production. Section 1 bounds each entry by $`W(J_{1/i})/(s-\rho)`$. This establishes absolute convergence using only the explicitly given rates.

## 3. Finite rational lower and upper matrices

Fix $`N\ge2`$. Include all planted types reachable from the nineteen newborns while their size is at most $`N`$, using the transitions with $`t(S,S')>0`$. Discard a particle at its first intrinsic transition to size $`N+1`$. Copies have size at most their source size. Let $`R_N`$ be the residual first-moment generator killed in this way, acting on column scores. Its diagonal subtraction is $`a_S+r_0(S)`$, also on the boundary, and its positive entries are $`t(S,S')`$ for retained successors, with any self-transition added there. The finite matrix $`A_N=sI-R_N`$ is invertible with nonnegative inverse, because the killed finite semigroup is bounded by the full semigroup in weight $`W`$ and $`s>\rho`$. For a rational $`s`$ all entries of the inverse are rational.

Write $`B_j(S)=b_S\mathbf1_{\{m(D_{\rm root})=j\}}`$. Let $`J`$ denote restriction of a column to the nineteen newborns and define
```math
K^-_N(s)_{\cdot j}=J A_N^{-1}B_j,\qquad
H_N(S)=\sum_{n(S')=N+1} t(S,S')\,\frac{W(S')}{s-\rho},\qquad
u_N=J A_N^{-1}H_N.
\qquad\text{(3)}
```
Set $`K^+_N(s)_{\cdot j}=K^-_N(s)_{\cdot j}+u_N`$. In (3) each exiting particle is assigned future reward $`0`$ for the lower comparison, and an upper bound on its expected direct offspring of the prescribed label for the upper comparison. Decomposition at first exits on each lineage gives, entry by entry,
```math
K^-_N(s)\le K(s)\le K^+_N(s),\qquad
r(K^-_N(s))\le r(K(s))\le r(K^+_N(s)).
\qquad\text{(4)}
```
For a fixed test column $`g\ge0`$, a sharper upper comparison adds $`u_N\|g\|_\infty`$ to $`K^-_Ng`$ instead of applying $`K^+_N`$ to $`g`$. Any other verified positive supersolution with smaller boundary values can replace $`W/(s-\rho)`$ for the corresponding source. Both matrix bounds in (4) are nested with $`N`$: revealing one more intrinsic transition and copy generation inside the larger domain replaces the previous bound on the boundary by an improved bound, by $`(s-R)W/(s-\rho)\ge W\ge B_j`$ and positivity of the killed resolvent.

## 4. A rate of convergence

For $`0<\delta\le1`$ put $`V_\delta(S)=n(S)^\delta W(S)`$. The intrinsic generator obeys
```math
L^0 V_\delta=(n+1)^\delta L^0W+
 \big((n+1)^\delta-n^\delta\big)r_0(S)W .
```
Each copied state has size at most $`n`$, and $`\sum_Y a_YW(Y)=(a_S+d)W(S)-h(c_{\rm root})`$. Using $`(1+1/n)^\delta\le1+\delta/n`$, $`r_0\le n`$ and $`n\ge2`$ gives
```math
R V_\delta\le\big[\rho+\delta(1+q/2)\big]V_\delta .
\qquad\text{(5)}
```
Choose $`\rho+\delta(1+q/2)\le s`$. Stop each residual lineage on its first exit at size $`N+1`$. The expected sum of discounted $`V_\delta`$ of these exiting particles is at most $`V_\delta(S_0)`$. For instance leave the total discounted weights of past exits in an absorbing account, apply (5) to the particles still alive in the finite domain and integrate the drift. Contributions from their births inside the domain are exactly those of the residual first-moment generator. The bound follows first at finite time and then by monotone convergence of the account. At an exit, $`W=(N+1)^{-\delta}V_\delta`$. Consequently
```math
0\le u_N(i)\le \left(\frac{2}{N+1}\right)^\delta
       \frac{W(J_{1/i})}{s-\rho}. \qquad\text{(6)}
```

For clarity, the rate can be passed to the spectral radii with constants independent of unknown spectral data. Use $`l=113/200`$, $`u=1133/2000`$, $`\delta=1/24`$, so that (5) is at most $`lV_\delta`$. On $`l\le s\le u`$ all single entries of $`K(s)`$ are at most $`H=1/[d(l-\rho)]`$. Here is a common strictly positive lower bound $`e`$ on the entries. In the intrinsic original particle use only root leaf gains, from degree two up to degree $`j`$, then produce directly at degree $`j`$. At each intermediate degree $`k`$, with no core children, $`L=k-1`$, $`r_0=k-1+(k-1)/k`$ and the root gain rate is $`k-1`$. Its death rate is at most $`1/2`$. The direct rate at the terminal state is $`(j-1)/j`$. Dropping all other paths gives
```math
e=\min_{2\le j\le20}\!
 \left[\prod_{k=2}^{j-1}
 \frac{k-1}{u+1/2+k-1+(k-1)/k}\right]
 \frac{(j-1)/j}{u+1/2+j-1+(j-1)/j}>0 .
```
Normalise a positive eigenvector of $`K(s)`$ by $`\max g=1`$. The radius is at most $`19H`$. Each row at a maximal input contributes at least $`e`$, so $`g_{\min}\ge g_*:=e/(19H)`$. Comparing (4) on $`g`$ and using positivity gives
```math
 0\le r(K^+_N(s))-r(K^-_N(s))
 \le \frac{38H}{g_*}\left(\frac{2}{N+1}\right)^{1/24}. \qquad\text{(7)}
```
Indeed every row error applied to $`g`$ is at most $`19\max_i u_N(i)`$, and separate Collatz bounds are valid even for reducible matrices and for a negative trial lower bound.

There is also a uniform inverse slope. For $`l\le s<s'\le u`$, evaluate $`K(s)-K(s')`$ on the positive Perron vector of $`K(s')`$ with maximum one. Direct production of label 2 in the original newborn before any intrinsic change or killing contributes at least
```math
g_*\int_0^\infty \tfrac12(e^{-st}-e^{-s't})e^{-2t}\,dt
 \ge \frac{g_*}{2(u+2)^2}(s'-s).
```
Every other contribution has the same sign. Since $`g\le1`$, comparison on $`g`$ shows that $`r(K(s))-r(K(s'))\ge g_*(s'-s)/[2(u+2)^2]`$. Thus (4), with an increasingly fine rational grid or bisection, gives rigorous lower and upper enclosures for $`\mu`$ which converge with an explicit power rate, using only rational finite matrices and real-root comparisons. Constants (6)–(7) are coarse. They establish effective convergence and must not be read as an eight-digit computation.
