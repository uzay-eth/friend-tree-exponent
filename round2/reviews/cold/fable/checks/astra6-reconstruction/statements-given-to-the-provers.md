# Random friend tree: definitions and statements

## Model and notation

Start at size two with an edge. At size $`n`$, choose a uniform vertex $`V`$, choose a uniform neighbour of $`V`$, and attach a new leaf to that neighbour. In continuous time every vertex rings at rate one and chooses its neighbour at its ring. Write $`T_n`$ for the discrete tree, $`T_t`$ when the subscript denotes elapsed time, and $`\mathbb E_T`$ for expectation in continuous time from the finite tree $`T`$. The seed at size three is denoted by $`T_3`$. Label the edge endpoints 1 and 2 and remaining vertices in order of arrival. At a vertex write $`D`$ for total degree, $`L`$ for the number of leaf neighbours and $`c`$ for the number of core neighbours; a core is a vertex of degree at least two. A promotion is a ring of a core choosing a leaf neighbour, raising that neighbour to degree two. Let $`C=N_{\ge2}`$ count cores, $`N_{\ge k}`$ count degree at least $`k`$, and $`N_{=k}=N_{\ge k}-N_{\ge k+1}`$. We use the deterministic exponent $`\mu`$ defined by $`\log\mathbb E C_n/\log n\to\mu`$. The known bounds and logarithmic limits are $`113/200\le\mu\le1133/2000`$, $`\log C_n/\log n\to\mu`$ almost surely, and $`t^{-1}\log\mathbb E_{T_3} C(T_t)\to\mu`$.

## S2. Marked fringes and renewal (including the capped version S5)

**S2.1 (states and rates).** Root the physical tree at the centre $`o`$ of $`T_3`$. At each promotion of $`v`$ by its parent $`u`$, assign the permanent mark $`a_v=1/D_u`$ using the degree at promotion. At later times $`1/D_u\le a_v`$. For a noninitial core $`v`$, its planted fringe $`S_v`$ is the rooted tree of $`v`$ and all its descendants, with all leaves and core marks, and an extra parent slot at its root counted in both $`D`$ and $`c`$ but not as a vertex. Auxiliary types are finite such marked planted trees: the root mark is in $`\mathcal A=\{1/2,1/3, \ldots\}\cup\{0\}`$ with its usual subspace topology, and each other core mark is a nonzero element of $`\mathcal A`$ at least the current reciprocal degree of its parent. The root mark zero is an auxiliary value. For a type $`S`$ use $`a_S`$ for its root mark, $`p_S=1/D_{\rm root}`$ and $`b_S=L_{\rm root}/D_{\rm root}`$. Write $`Y`$ for the planted subtree of an immediate core child of the root (counted with multiplicity), $`Y^+`$ for $`Y`$ with one extra root leaf and unchanged old marks, and $`J_a`$ for one core with mark $`a`$, one leaf and one parent slot.

The intrinsic process, with generator $`L^0`$, has the physical rates within the planted fringe. A leaf rings and gives its parent a new leaf at rate one; a core sends a new leaf to a specified core neighbour inside the fringe at rate $`1/D`$; and a core promotes a specified one of its own leaves at rate $`1/D`$ (total promotion rate $`L/D`$), with its newborn mark assigned as above. Every intrinsic change updates the whole type. The exterior parent slot contributes to the root's degree, but sends nothing and makes no internal change when selected by the root.

In the residual branching process each particle evolves intrinsically and dies at the constant rate $`a_S`$ given by its own root mark. For each current core child $`Y`$ of its root it emits an additional particle of type $`Y^+`$ at rate $`p_S`$ and of type $`Y`$ at rate $`a_Y-p_S`$. Source particles remain in place; emitted particles use the inherited marks in their types and evolve independently. Its first-moment semigroup $`R_t`$ maps a score $`f`$ to the expected sum of $`f`$ over particles alive at time $`t`$, starting from a single particle of the given type. Its generator is
```math
Rf(S)=L^0 f(S)-a_S f(S)+\sum_Y\big[(a_Y-p_S)f(Y)+p_S f(Y^+)\big].
```
The full branching process further emits $`J_{p_S}`$ at rate $`b_S`$, retaining the source particle; these are the direct births. Its added generator is $`Bf(S)=b_S f(J_{p_S})`$, and its first-moment semigroup is $`P_t`$, with generator $`R+B`$.

**S2.2 (weight and bounds).** Put $`d=3/20`$, $`q=9/25`$, $`\rho=q+d=51/100`$, $`z(a)=a/(a+d)`$ and
```math
h(0)=2, \quad h(1)=1, \quad h(2)=17/20, \quad h(3)=39/50, \quad h(c)=7/10\ (c\ge4).
```
Let $`w_x`$ be the product of $`z(a_y)`$ along the path from the root to $`x`$, excluding the root and including $`x`$, with $`w_{\rm root}=1`$. Set $`W(S)=(a_S+d)^{-1}\sum_{x\ {\rm core\ in}\ S} w_xh(c_x)`$. The sum includes the root; $`c_x`$ includes the root's parent slot when $`x`$ is the root. The corresponding unnormalised weight of a full physical tree is $`Z(T)=\sum_{x\ {\rm core}}w_x h(c_x)`$, where its founder has no parent slot. The branching processes are nonexplosive and satisfy
```math
1<\frac{7/10}{1/2+d}\le W(S)\le |S|/d, \quad BW\le 1/d, \quad L^0 W\le qW, \quad
RW\le\rho W-h(c_{\rm root}), \quad R_tW\le e^{\rho t}W, \quad P_tW\le e^{(\rho+1/d)t} W.
```
Here $`|S|`$ counts vertices of $`S`$ without the parent slot.

**S2.3 (physical first moment).** From $`T_3`$ write $`M_tf=\mathbb E_{T_3}\sum_{v\ne o, \ v\ {\rm core}} f(S_v(t))`$. Define
```math
I_t f=\mathbb E_{T_3}\! \left[b_o f(J_{p_o})+\sum_{Y\text{ at }o}\big((a_Y-p_o)f(Y)+p_o f(Y^+)\big)\right],
```
where $`p_o=1/D_o`$, $`b_o=L_o/D_o`$ are for the founder and the child subtrees are those at time $`t`$. Then $`I_t W\le 1/d+2e^{qt}`$ and, for scores bounded in absolute value by a constant times $`W`$,
```math
M_t f=\int_0^t I_u P_{t-u} f\,du.
```

**S2.4 (renewal characterization).** For $`s>\rho`$ and $`g\in C(\mathcal A)`$ let
```math
U_sg(S)=\int_0^\infty e^{-st}R_tb_{(\cdot)}g(p_{(\cdot)}) (`S`)\,dt, \qquad (K_s g)(a)=U_s g(J_a).
```
The bracket denotes the type score $`S'\mapsto b_{S'}g(p_{S'})`$. $`K_s`$ is a compact positive operator on $`C(\mathcal A)`$ with the supremum norm. For real $`s>51/100`$, $`r(K_s)=1`$ if and only if $`s=\mu`$. The mean core count from the edge is $`\mathbb E C_n\asymp n^\mu`$.

**S5.1 (nineteen labels and determinant).** Set $`m(D)=\min(D,20)`$. In the capped variant, assign new core mark $`1/m(D_{\rm parent})`$ both physically and intrinsically, use $`J_{1/m(D_{\rm root})}`$ for a direct birth, and make the same replacement for the founder's direct immigration. Leave $`p_S=1/D_{\rm root}`$ and the actual send and copy rates unchanged in form. The weight, semigroup bounds and physical first-moment identity above hold in this variant. Number newborn types by labels $`i=2, \ldots,20`$, with type $`J_{1/i}`$. The matrix $`K(s)`$ has entry $`K_{ij}(s)`$ equal to the expected total $`e^{-st}`$-discounted number of direct offspring of label $`j`$ at their birth times $`t`$, from the residual family begun at newborn label $`i`$. For $`s>51/100`$, $`r(K(s))=1`$ exactly at $`s=\mu`$, and $`\mu`$ is the largest real zero there of $`\det(I_{19}-K(s))`$.

**S5.2 (entries and enclosures).** Let $`n(S)=|S|`$, let $`r^0(S, S')`$ be the total intrinsic rate from $`S`$ to $`S'`$, and put $`r_0(S)=\sum_{S'}r^0(S, S')=n(S)-1/D_{\rm root}`$. Put
```math
t(S, S')=r^0(S, S')+\sum_{Y\text{ at the root}}\big[p_S\mathbf1_{Y^+=S'}+(a_Y-p_S)\mathbf1_{Y=S'}\big],
```
with the sum over core children. Then, with only positive transitions included and the empty product equal to one,
```math
K_{ij}(s)=\sum_{\ell=0}^\infty
\sum_{\substack{S_0=J_{1/i}, S_1, \dots, S_\ell\,;\ m(D_{\rm root}(S_\ell))=j}}
\left(\prod_{u=0}^{\ell-1}\frac{t(S_u, S_{u+1})}{s+a_{S_u}+r_0(S_u)}\right)
\frac{b_{S_\ell}}{s+a_{S_\ell}+r_0(S_\ell)}.
```
This is absolutely convergent for $`s>\rho`$. For integer $`N\ge2`$ let the retained types be those reachable from the newborns using positive transitions while $`n(S)\le N`$. On these types let $`R_N`$ be the residual generator with discarded successors at size $`N+1`$: use diagonal subtraction $`a_S+r_0(S)`$ and add $`t(S, S')`$ for retained successors. Set $`A_N=sI-R_N`$, let $`J`$ restrict a score to the nineteen newborns, and put
```math
B_j(S)=b_S\mathbf1_{m(D_{\rm root})=j}, \quad H_N(S)=\sum_{n(S')=N+1}t(S, S')\frac{W(S')}{s-\rho}, \quad u_N=JA_N^{-1}H_N,
```
```math
K_N^-(s)_{\cdot j}=JA_N^{-1}B_j, \qquad K_N^+(s)_{\cdot j}=K_N^-(s)_{\cdot j}+u_N.
```
$`A_N^{-1}`$ exists, is nonnegative and is rational for rational $`s`$. The nested entrywise bounds $`K_N^-\le K\le K_N^+`$ bound the radii in the same order. For $`0<\delta\le1`$ and $`s\ge\rho+\delta(1+q/2)`$ with $`s>\rho`$,
```math
0\le u_N(i)\le (2/(N+1))^\delta\,W(J_{1/i})/(s-\rho).
```
With $`l=113/200`$, $`u=1133/2000`$, $`H_*=1/[d(l-\rho)]`$, define
```math
e_* =\min_{2\le j\le20}\left[\prod_{k=2}^{j-1}\frac{k-1}{u+1/2+k-1+(k-1)/k}\right]\frac{(j-1)/j}{u+1/2+j-1+(j-1)/j}, \qquad g_*=e_*/(19H_*).
```
On $`l\le s<s'\le u`$, the bounds are
```math
0\le r(K_N^+(s))-r(K_N^-(s))\le \frac{38H_*}{g_*}(2/(N+1))^{1/24}, \qquad
r(K(s))-r(K(s'))\ge\frac{g_*}{2(u+2)^2}(s'-s),
```
with the width estimate also at $`s=u`$.

## S3. Sharp means and the limiting degree law

**S3.1 (constants and limits).** Use exact marks as in S2. For $`\Re z>\rho`$, set $`G_z f=\int_0^\infty e^{-zt}R_t f\,dt`$, $`\widehat I_z(F)=\int_0^\infty e^{-zt}I_t F\,dt`$, and $`(jF)(a)=F(J_a)`$. Extend the definitions of $`U_z, K_z`$ by the same formulas with complex $`z`$. Let $`g, \ell`$ be positive right and left eigenvectors of $`K_\mu`$ for eigenvalue one (with $`\ell`$ interpreted as a positive measure), normalised by $`\ell(g)=1`$, and let $`K'_\mu`$ be the derivative with respect to the discount. With $`f_k(S)=\mathbf1_{D_{\rm root}\ge k}`$, define
```math
d_f=\frac{\widehat I_\mu(U_\mu g)\,\ell(jG_\mu f)}{-\ell(K'_\mu g)}, \qquad d_k^*=\frac{2}{2+\mu}d_{f_k}, \qquad
c_k=\frac{d_k^*}{\Gamma(\mu+2)}, \qquad p_k=\frac{c_k-c_{k+1}}{c_2}.
```
Here $`\Gamma`$ is the gamma function. For each fixed $`k\ge2`$, $`\mathbb E_{T_3}N_{\ge k}(t)\sim d_{f_k}e^{\mu t}`$, $`\mathbb E_{T_2}N_{\ge k}(t)\sim d_k^* e^{\mu t}`$, and
```math
\mathbb E N_{\ge k}(n)\sim c_k n^\mu, \qquad c_k>c_{k+1}>0, \qquad
\lim_{n\to\infty}\frac{\mathbb E N_{=k}(n)}{\mathbb E C_n}=p_k.
```
In particular $`\mathbb E C_n\sim c_2 n^\mu`$, $`p_k>0`$, $`\sum_{k\ge2}p_k=1`$. The sequence of point probabilities just defined satisfies $`p_k\sim\mu A k^{-1-\mu}`$ as $`k\to\infty`$, with $`A`$ as below. The fixed-$`k`$ limit in $`n`$ precedes the large-$`k`$ limit.

**S3.2 (rooted core laws).** The spatial state space $`\mathcal X`$ consists of isomorphism classes of finite or countably infinite, connected, locally finite rooted core trees. Each core has a nonnegative integer mark $`L_x`$ giving its number of leaves, with $`D_x=L_x+c_x\ge2`$ where $`c_x`$ is its degree in the core. A one-vertex core is allowed, and the topology is agreement on finite marked balls. These are leaf-count marks, not the birth marks of S2. For $`x`$ in the core, a gain of a new leaf has rate $`L_x+\sum_{y\sim x\ {\rm in\ the\ core}}1/D_y`$ and a promotion has rate $`b_x=L_x/D_x`$. A promotion replaces $`(L_x, c_x)`$ by $`(L_x-1, c_x+1)`$ and attaches a new adjacent core with $`(L, c)=(1,1)`$.

Assign each old core a different label and let a newly promoted core inherit the label of its promoter. The cell $`H_o(t)`$ is the set of cores carrying the initial label of the root. Write $`(Q_t f)(G, o)=\mathbb E_{G, o}\sum_{v\in H_o(t)}f(G_t, v)`$ for bounded scores on rooted states. This defines a positive semigroup on $`C_b(\mathcal X)`$ (bounded continuous functions), with $`Q_t1\le e^t`$, independent of a boundary condition at infinity. There is a Borel probability $`\pi`$ with $`\pi Q_t=e^{\mu t}\pi`$ for every $`t\ge0`$. The spectral radius of $`Q_1`$ on $`C_b(\mathcal X)`$ is $`e^\mu`$, and every probability eigenmeasure of $`Q_1`$ at a positive eigenvalue has that eigenvalue. For $`0<\alpha<\min(1, \mu)`$ and $`0<\theta<\log(1+\mu)`$ the law $`\pi`$ has finite $`D_o^\alpha`$ and $`e^{\theta c_o}`$ moments, and $`\pi(b_o)=\mu`$.

To specify the probability used in S3, let $`m_t(f)=\mathbb E_{T_3}\sum_{x\ {\rm core}}f(G_t, x)`$ with $`F_3(t)=m_t(1)`$. For $`0<\delta<\mu`$ and $`T>0`$, let $`w(t)=\exp(-\mu t-\delta |t-T|)`$, $`Z^{\rm tilt}=\int_0^\infty w(t)F_3(t)\,dt`$, and $`\pi_{\delta, T}=(Z^{\rm tilt})^{-1}\int_0^\infty w(t)m_t\,dt`$. The law $`\pi`$ in S3 is a weak limit along a sequence with $`\delta\downarrow0`$ and $`w(0)/Z^{\rm tilt}\to0`$. Its degree masses are $`\pi(D_o\ge k)=c_k/c_2`$ and $`\pi(D_o=k)=p_k`$.

**S3.3 (newborn law and degree constant).** Define $`S_t f(G, o)=\mathbb E_{G, o} f(G_t, o)`$ for an old root kept fixed, and define the spatial birth operator $`B^{\rm sp} f(G, o)=b_o f(G^{b_o}, v)`$, with $`G^{b_o}`$ the rooted core's underlying tree after promotion at $`o`$ and $`v`$ the newborn root in the score evaluation. Define the newborn probability $`\beta=\pi B^{\rm sp}/\mu`$. The root under $`\beta`$ has one leaf and one core neighbour. For bounded nonnegative $`f`$,
```math
\pi f=\int_0^\infty e^{-\mu a}\pi B^{\rm sp}S_a f\,da.
```
For the fixed newborn root from $`\beta`$, put $`W^{\rm deg}=\lim_{t\to\infty} e^{-t} D(t)`$, the almost-sure finite limit of the nonnegative normalized degree supermartingale. Set $`A=\mathbb E_\beta[(W^{\rm deg})^\mu]`$; then $`0<A<\infty`$. For $`T_k=\inf\{t:D(t)\ge k\}`$,
```math
\pi(D_o\ge k)=\mathbb E_\beta[e^{-\mu T_k};T_k<\infty]\sim A k^{-\mu}, \qquad
\pi(D_o=k)\sim\mu A k^{-1-\mu}.
```
For $`k\ge2`$, with $`\rho_o=\sum_{y\sim o\ {\rm in\ the\ core}}1/D_y`$,
```math
\mu\pi(D_o\ge k+1)=\pi[(k-c_o+\rho_o)\mathbf1_{D_o=k}].
```

## S6. Almost-sure limits and expectation estimates

**S6 (theorem).** From the edge there is a finite random variable $`X>0`$ almost surely such that
```math
C_n/n^\mu\longrightarrow X, \qquad N_{\ge k}(n)/C_n\longrightarrow c_k/c_2, \qquad N_{=k}(n)/C_n\longrightarrow p_k
```
almost surely, simultaneously for every fixed $`k\ge2`$. For every $`\varepsilon, M>0`$, $`\mathbb P(X<\varepsilon)>0`$ and $`\mathbb P(X>M)>0`$. The limiting $`p_k`$ are as in S3, including their large-$`k`$ asymptotic.

The following statements use capped marks throughout. Use the weight $`W`$ and constants $`d, q, \rho, h`$ of S2.2. Put $`B_0 g(S)=b_S g(m(D_{\rm root}))`$, and let $`J`$ restrict to the nineteen newborns, with labels as in S5. Write $`k(t)=JR_t B_0`$ (matrix); its Laplace transform is $`K(z)`$ from S5. In this section $`g, \ell`$ are the positive Perron column and row of the capped $`K(\mu)`$ with $`\ell g=1`$. Define
```math
D^{\rm ker}=\frac{g\otimes\ell}{-\ell K'(\mu)g}, \qquad
v(t)=\sum_{j\ge1} k^{*j}(t), \qquad h_f(t)=JR_tf, \qquad b(t)=R_t B_0.
```
Here $`g\otimes\ell`$ is the matrix of their outer product and convolution uses matrix/operator products in time order, $`(a*b)(t)=\int_0^t a(u)b(t-u)du`$.

**Lemma S6.1: kernel.** The entries of $`k(t)`$ are continuous and bounded by $`Ce^{\rho t}`$. For $`\max_i|g_0(i)|\le1`$, $`\|R B_0g_0\|_W\le C`$, with $`\|F\|_W=\sup_S |F(S)|/W(S)`$, and
```math
k(t)g_0=JB_0g_0+\int_0^t JR_u R B_0g_0\,du.
```
For each $`\varepsilon>0`$, $`K(z)=O(1/|z|)`$ uniformly on $`\Re z\ge\rho+\varepsilon`$ at large $`|z|`$. There exists $`\rho<\sigma<\mu`$ such that on $`\Re z\ge\sigma`$ the only singularity of $`(I-K(z))^{-1}`$ is a simple pole at $`\mu`$ of residue $`D^{\rm ker}`$. For all $`t\ge0`$,
```math
v(t)=e^{\mu t}D^{\rm ker}+O(e^{\sigma t}).
```

**Lemma S6.2: branching first moments.** For every score $`|f|\le W`$, set
```math
d(f)=\frac{\ell\int_0^\infty e^{-\mu u}h_f(u)\,du}{-\ell K'(\mu)g}, \qquad H(S)=\int_0^\infty e^{-\mu u}(R_u B_0g)(S)\,du.
```
For $`t\ge0`$,
```math
JP_tf=h_f(t)+(v*h_f)(t)=e^{\mu t}g\,d(f)+O(e^{\sigma t}), \qquad P_t f=R_tf+(b*JP_\cdot f)(t),
```
```math
|P_t f(S)-e^{\mu t}H(S)d(f)|\le C e^{\sigma t}W(S).
```
The error constants are uniform in $`|f|\le W`$, and the latter bound holds for all planted starting states.

**Lemma S6.3: physical means, uniform in the starting tree.** Start at any finite tree $`T`$ with at least one core, rooted at a core $`o`$, with capped marks obeying the parent mark inequality. Let $`\mathcal I(T)F`$ be the expression inside the expectation defining $`I_tF`$ in S2.3, with the direct birth rounded as in S5, and put $`I_{T, t}F=\mathbb E_T\mathcal I(T_t)F`$. Then $`\sum_{v\ne o, \ v\ {\rm core}} W(S_v)\le C_0 C(T)`$, $`Z(T)\le 2C(T)`$, $`\mathbb E_T Z(T_t)\le e^{qt}Z(T)`$, and
```math
I_{T, t}W\le 1/d+e^{qt}Z(T).
```
The mean excluding the founder is $`\sum_{v\ne o, \ v\ {\rm core}}P_tf(S_v)+\int_0^t I_{T, u}P_{t-u}f\,du`$. Define
```math
\Psi(T)=\sum_{v\ne o, \ v\ {\rm core}}H(S_v)+\int_0^\infty e^{-\mu u}I_{T, u}H\,du.
```
For all $`t\ge0`$ and all $`|f|\le W`$,
```math
\left|\mathbb E_T\sum_{v\ne o, \ v\ {\rm core}}f(S_v(t))-e^{\mu t}\Psi(T)d(f)\right|\le C_1 C(T)e^{\sigma t}.
```
The error is uniform in $`T`$ and $`f`$. Put $`A_k(T)=N_{\ge k}(T)`$ and $`d_k=d(f_k)>0`$ in the capped construction. The same bound with an increased constant holds for $`\left|\mathbb E_T A_k(T_t)-e^{\mu t}\Psi(T)d_k\right|`$, including $`k=2`$ with $`A_2=C`$.

**Lemma S6.4: eigenfunction.** $`0\le\Psi(T)\le C_2 C(T)`$, $`C(T)\le C_3(1+\Psi(T))`$ and $`\Psi(T_3)>0`$. Moreover $`\Psi(T)=d_2^{-1}\lim_{u\to\infty}e^{-\mu u}\mathbb E_T C(T_u)`$ is independent of compatible orientation or marks, and $`\mathbb E_T\Psi(T_h)=e^{\mu h}\Psi(T)`$ for $`h\ge0`$. Consequently $`Z_t^*=e^{-\mu t}\Psi(T_t)`$ is a nonnegative martingale with an almost-sure finite limit $`Z_\infty^*`$.

**Lemma S6.5: one-leaf influence.** For any finite tree $`T`$ with a core and any one-leaf extension $`T^+`$ there is a coupling such that, for $`u\ge0`$ and $`k\ge2`$,
```math
|\mathbb E_{T^+} A_k(T_u)-\mathbb E_T A_k(T_u)|\le 3e^{2u}.
```
Under the coupling the clean cores are paired and contribute equally to $`A_k`$ in both copies, and the expected total number of dirty cores across both copies is at most $`3e^{2u}`$; a dirty core is one excluded from this clean pairing.

**Lemma S6.6: variance.** With $`n(T)`$ the vertex count, for $`b\ge0`$ and $`k\ge2`$,
```math
\mathrm{Var}_T(A_k(T_b))\le 9n(T)\int_0^b e^{4(b-v)+v}\,dv .
```
For $`g_h(T)=d_2^{-1}e^{-\mu h}\mathbb E_T C(T_h)`$ and $`\eta=\mu-\sigma>0`$, the following bounds hold uniformly in $`T`$ and $`h\ge0`$:
```math
\mathrm{Var}_T(g_h(T_1))\le C_4 e^{4h}n(T), \quad |\Psi(T)-g_h(T)|\le C_5e^{-\eta h}C(T), \quad
\mathbb E_T C(T_1)^2\le C_6 C(T)^2,
```
```math
\mathrm{Var}_T\Psi(T_1)\le C_7\big(e^{4h}n(T)+e^{-2\eta h}C(T)^2\big).
```

**Lemma S6.7: martingale second moments and positivity.** For integer $`j\ge0`$, $`C(T_j)^2\le C_8(1+\Psi(T_j)^2)`$ and $`\mathbb E_T n(T_j)=e^j n(T)`$. There are finite constants such that
```math
\mathbb E_T(Z_{j+1}^*)^2\le (1+C_9(j+2)^{-2})\mathbb E_T(Z_j^*)^2+C_{10}\big[(j+2)^{4/\eta}e^{-(2\mu-1)j} n(T)+(j+2)^{-2}e^{-2\mu j}\big],
```
```math
\mathbb E_T(Z_\infty^*)^2\le C_{11}(\Psi(T)^2+n(T)+1), \qquad \mathbb E_T Z_\infty^*=\Psi(T).
```
If $`\mathcal F_j`$ is the continuous-time history at time $`j`$, then
```math
\mathbb P_T(Z_\infty^*>0\mid\mathcal F_j)\ge
\frac{\Psi(T_j)^2}{C_{11}(\Psi(T_j)^2+n(T_j)+1)}.
```
In particular from $`T_3`$ the limit is strictly positive almost surely.

**Lemma S6.8: conditional count estimates.** On a fixed mesh $`t=j/m`$ with integer $`m\ge1`$, set $`h=\eta^{-1}\log(t+2)`$ and $`u=t-h>0`$. Then for $`k\ge2`$,
```math
\mathbb E\left[e^{-2\mu t}(A_k(T_t)-\mathbb E_{T_u}A_k(T_h))^2\right]\le C_{12} n(T_0)e^{-(2\mu-1)t+4h}.
```
The difference between $`e^{-\mu t}\mathbb E_{T_u} A_k(T_h)`$ and $`d_k Z_u^*`$ is bounded by a constant times $`e^{-\eta h}e^{-\mu u} C(T_u)`$. From $`T_3`$, $`e^{-\mu t}A_k(T_t)\to d_k Z_\infty^*`$ as $`t\to\infty`$ almost surely, simultaneously for all $`k\ge2`$. If $`\tau_n`$ is the epoch of size $`n`$ from $`T_3`$ and $`Y_\infty^{\rm size}=\lim_n ne^{-\tau_n}>0`$ is the Yule size limit, then $`A_k(T_{\tau_n})/n^\mu\to d_k Z_\infty^*/(Y_\infty^{\rm size})^\mu`$ almost surely. The deterministic ratios satisfy $`d_k/d_2=c_k/c_2`$.

**Lemma S6.9: discrete martingale and conditional mean.** Let $`T_m`$ denote the tree at size $`m\ge3`$ and $`\tau_m`$ its size-$`m`$ time from $`T_3`$. Then
```math
\mathbb E[\Psi(T_{m+1})\mid T_m]=(1+\mu/m)\Psi(T_m), \qquad
\mathbb E e^{-2\mu\tau_m}=\prod_{j=3}^{m-1}\frac{j}{j+2\mu}\asymp m^{-2\mu},
```
```math
\mathbb E(Z_{\tau_m}^*)^2\le\mathbb E(Z_\infty^*)^2, \qquad
\mathbb E[X\mid T_m]=d_2\frac{\Gamma(m)}{\Gamma(m+\mu)}\Psi(T_m).
```
The nonnegative discrete martingale $`Y_m=\Gamma(m)\Psi(T_m)/\Gamma(m+\mu)`$ is bounded in $`L^2`$ and $`X=d_2\lim_m Y_m`$. Conditional on the size-$`m`$ tree being a star, the conditional mean tends to zero as $`m\to\infty`$; on a path it grows at least as a positive constant times $`m^{1-\mu}`$ for all sufficiently large $`m`$.

## S7. Hub mass and next attachment

**S7 (theorem).** For every label $`v`$, let $`Z_v=\lim_{n\to\infty}D_v(n)/n`$, the known almost-sure fixed-label limit. Then $`\sum_{v\ge1}Z_v=1`$ almost surely. Write $`W_n`$ for the recipient label in the attachment from size $`n`$ to size $`n+1`$. Almost surely the conditional law of $`W_n`$ given $`T_n`$ converges in total variation to $`(Z_v)_{v\ge1}`$ and the conditional law of $`D_{W_n}(n)/n`$ converges weakly to $`\sum_v Z_v\delta_{Z_v}`$, with $`\delta_z`$ denoting unit mass at $`z`$. The mass conclusion holds using just $`C_n=O(n^\lambda)`$ almost surely for deterministic $`\lambda<1`$, in place of an expectation power bound.

**Lemma S7.1: normalized degree and expectation tail.** A core promoted at size $`r`$ has degree two at size $`r`$, and from promotion onwards its normalized degree $`D_v(j)/j`$ is a nonnegative supermartingale in the discrete tree filtration. Write $`a_j=\mathbb E C_j`$. If $`a_j\le K j^\lambda`$ for $`j\ge3`$ with $`\lambda<1`$, then $`a_{j+1}-a_j\le a_j/j`$ and for $`n\ge m\ge3`$,
```math
\mathbb E \sum_{\substack{v>m\\D_v(n)\ge2}}\frac{D_v(n)}n\le 2\sum_{r=m+1}^n\frac{a_r-a_{r-1}}r\le 2K\sum_{j=m}^\infty\frac{j^{\lambda-1}}{j+1}\le 2K\left(m^{\lambda-2}+\frac{m^{\lambda-1}}{1-\lambda}\right).
```
The degree sum over cores is $`n+C_n-2`$. The fixed-label leaf-neighbour limits are $`L_v(n)/n\to Z_v`$ almost surely. The conditional attachment probability is
```math
\mathbb P(W_n=v\mid T_n)=\frac{L_v(n)+\sum_{w\sim v, \ D_w(n)\ge2}1/D_w(n)}n.
```

**Lemma S7.2: localized tail.** For $`K>0`$ and deterministic $`\lambda<1`$ put $`\tau=\inf\{n\ge3:C_n>K n^\lambda\}`$ (infimum of the empty set is infinity). At size $`m`$ start a sum with no counted cores. On $`\{\tau>m\}`$ include the normalized degrees of cores promoted after $`m`$ through the stopping step inclusive, and freeze the sum at $`\tau`$; on $`\{\tau\le m\}`$ it is zero. For $`n\ge m\ge3`$ the expected stopped sum at size $`n`$ is at most
```math
2\mathbb E\sum_{r=m+1}^n\frac{C_{r\wedge\tau}-C_{(r-1)\wedge\tau}}r
=O_\lambda((K+1)m^{\lambda-1}).
```
The bound is uniform in $`n`$, and
```math
\mathbb E\left[\mathbf1_{\{\tau=\infty\}}\left(1-\sum_{v\le m}Z_v\right)\right]\le O_\lambda((K+1)m^{\lambda-1}).
```
