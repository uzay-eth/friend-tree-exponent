# Renewal at the birth of a non-leaf

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 2: Claude and Gemini, with a DeepSeek finite-generator check; inputs include Round 1 exponent and bracket. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Proposed result

In the undirected friend tree, write $`a_n=\mathbb E C_n`$. The construction below gives a compact positive next-generation operator $`K_s`$ on $`C(\mathcal A)`$, where $`\mathcal A=\{1/2,1/3,\ldots\}\cup\{0\}`$. Its kernel is defined independently of the growth exponent. For $`s>51/100`$, the equation $`r(K_s)=1`$ has exactly one real solution, $`s=\mu`$. It gives $`a_n\asymp n^\mu`$. A further argument for a limiting leading coefficient is given separately. An eight-digit numerical enclosure does not follow from the constants here.

The argument uses the repository's existence theorem and its lower bound $`\mu\ge113/200>51/100`$. All times below are continuous, with a rate-one clock at each vertex. [Existence theorem](../../proofs/existence.md) · [Lower bound](../../proofs/bracket-motifs.md)

## 1. Marked descendant trees

Root the actual tree at the centre of the three-vertex seed. When $`u`$ promotes a leaf $`v`$ to a core, put the permanent mark $`a_v=1/D_u`$ on $`v`$, using the degree at that instant. The degree of $`u`$ does not change at promotion and never decreases later, so every such edge has $`p_u:=1/D_u\le a_v`$. For a core $`v`$ other than the initial centre, $`S_v`$ consists of $`v`$ and all its descendants, including leaves and marks. Its root has one additional parent slot in its total degree and core degree. For any such finite type $`S`$ write $`a_S`$ for the root mark, $`b_S=L_o/D_o`$, $`p_S=1/D_o`$, and $`Y`$ for the types of its immediate core children, with multiplicities. $`Y^+`$ is $`Y`$ with a new leaf at its root. Write $`J_a`$ for a single core with one leaf, one parent slot and mark $`a`$. The value $`a=0`$ at the initial root of an auxiliary type serves as the limit of large parent birth degrees.

The intrinsic process on $`S`$, with generator $`L^0`$, leaves the exterior parent slot inert. All internal leaf gains, promotions, and core-to-core arrows have the ordinary rates. At the root its own clock may choose the parent slot with no internal change, and the slot sends no leaf back. Every intrinsic promotion adds the appropriate permanent mark to its new core. In particular the whole finite type, including its old subtrees, changes when a promotion occurs anywhere inside it. All types have finitely many internal transitions, with total rate at most their vertex count.

If $`f(S)`$ is summed over noninitial cores, the one input omitted by $`L^0`$ on $`S_v`$ is a gain at $`v`$ from its parent. Its contribution is $`p_{\mathrm{parent}}[f(S_v^+)-f(S_v)]`$. Moving this term into the parent's row gives, for old noninitial parents, the signed transfer $`\sum_Yp_S[f(Y^+)-f(Y)]`$. A promotion at the root of $`S`$ also creates a new summand $`f(J_{p_S})`$; charge it to that root at rate $`b_S`$. The contribution from the initial centre will be retained explicitly.

## 2. A positive first-moment system

On a particle of type $`S`$ run the following independent-particle branching process. Its type undergoes the intrinsic changes above. Kill the particle at rate $`a_S`$. For each current core child $`Y`$ produce an additional particle of type $`Y^+`$ at rate $`p_S`$, and one of type $`Y`$ at rate $`a_Y-p_S`$. The parent particle is retained at these births. These rates are nonnegative by the mark inequality. Call this process the residual process and write $`R_t`$ for its first-moment semigroup. Its generator on a score $`f`$ is

```math
R f=L^0 f-a_S f(S)+\sum_Y\big[(a_Y-p_S)f(Y)+p_S f(Y^+)\big]. \qquad\text{(1)}
```

The full process also produces $`J_{p_S}`$ at rate $`b_S`$, leaving the source particle in place. Write $`B f(S)=b_S f(J_{p_S})`$ and $`P_t`$ for the first-moment semigroup with generator $`R+B`$. The added birth and the intrinsic promotion may be driven separately because their *sum of first-moment changes* is the required generator. In particular the physical model is being compared in mean.

Let $`M_t f=\mathbb E\sum_{v\ne o} f(S_v(t))`$, summed over noninitial cores of the ordinary tree. A positive immigration measure from the initial centre $`o`$ is

```math
I_t f=\mathbb E\!\left[b_o f(J_{p_o})+\sum_{Y\text{ at }o}((a_Y-p_o)f(Y)+p_o f(Y^+))\right]. \qquad\text{(2)}
```

Summing (1), the terms $`a_Yf(Y)`$ from parents cancel the killing terms $`-a_vf(S_v)`$, including the boundary terms supplied by (2). What remains of the copy terms are the actual directed input transfers. The intrinsic evolution of every old summand and every new summand are now counted once. Thus $`M'_t f=M_t(R+B)f+I_t f`$. At $`T_3`$ there are no noninitial cores. The resulting identity is

```math
M_t f=\int_0^t I_u P_{t-u}f\,du. \qquad\text{(3)}
```

A justification of this identity for the unbounded system follows the estimates in the next section.

## 3. A subcritical estimate for the residual process

Put $`d=3/20`$, $`z(a)=a/(a+d)`$ and define

```math
h(0)=2,\quad h(1)=1,\quad h(2)=17/20,\quad h(3)=39/50,\quad h(c)=7/10\quad(c\ge4).
```

Give a core vertex $`x`$ in $`S`$ the factor $`w_x=\prod_{o<y\le x} z(a_y)`$ along its path from the root, with $`w_o=1`$. Set

```math
W(S)=\frac1{a_S+d}\sum_{x\in S\text{ core}} w_x h(c_x).
```

In this sum $`c_x`$ counts all core neighbours, including the parent slot at the root. A gain of a leaf does not change $`W`$. At a promotion by a vertex of degree $`D`$ and core degree $`c`$, the newborn receives factor $`1/(1+dD)`$ times its parent's factor. For all $`D\ge\max(2,c)`$ one has

```math
(1-c/D)\big[h(c+1)-h(c)+(1+dD)^{-1}\big]\le q h(c),\qquad q=9/25. \qquad\text{(4)}
```

For $`D\ge100`$ the left side is at most $`1/16`$, since $`h`$ is nonincreasing and $`1+dD\ge16`$. For $`D<100`$ (4) is a finite rational inequality. The included arithmetic checker tests all those pairs, including $`c=0`$ for the actual centre. This gives $`L^0 W\le qW`$. The child sum satisfies $`\sum_Y a_Y W(Y)=(a_S+d)W(S)-h(c_o)`$ by factoring the first edge on every path. Since $`W(Y^+)=W(Y)`$, (1) yields

```math
RW\le(q+d)W-h(c_o)\le\rho W,\qquad \rho=51/100. \qquad\text{(5)}
```

Also $`W(S)\ge(7/10)/(1/2+d)>1`$, $`W(S)\le |S|/d`$ for auxiliary types, and $`BW\le1/d`$. In particular $`R_tW\le e^{\rho t}W`$, and $`P_tW\le e^{(\rho+1/d)t}W`$. For these first-moment bounds, stop the branching population at finite event cutoffs and apply the nonnegative supermartingale bound for the sum of weights multiplied by the corresponding exponential discount. Its stopped mean bounds the probability that the particle count reaches a large level. While the count is bounded, deaths have bounded intensity and the number of copy births is bounded by the count plus the deaths. With finitely many particles born, each undergoes a nonexplosive intrinsic growth process. This rules out finite-time explosion. Fatou gives the asserted unstopped bounds.

For the actual full tree use $`Z=\sum_x w_x h(c_x)`$ from its centre. The same one-event inequality, now including $`c=0`$, gives $`\mathbb E Z(t)\le2 e^{qt}`$. In (2) the weighted sum over children is precisely $`Z-h(c_o)`$, and $`W(J_a)\le1/d`$. Consequently

```math
I_t W\le 1/d+2e^{qt}. \qquad\text{(6)}
```

For (3), first take $`0\le f\le C W`$. For every finite starting type the branching backward equation holds by conditioning on its first event. The initial total event rate is $`O(|S|)`$ and every one-event sum of successor weights is $`O(|S|+1)`$, so $`|\partial_s P_sf(S)|\le C_t(1+|S|)^2`$ on $`[0,t]`$. More explicitly the first-event integral, with finitely many initial transitions and the finite-time $`W`$ bound, gives continuity and the pointwise backward derivative. Apply the physical jump formula to the time-dependent scores $`P_{t-u}f`$, from $`u=0`$ to $`t`$. The original population is Yule with finite moments of every order on a bounded interval; these dominate the polynomial bounds on the sums of scores, derivatives and jump rates. The backward and forward generator terms cancel, leaving (3). This argument also applies to signed scores bounded by a multiple of $`W`$.

## 4. The compact reproduction operator

For real $`s>\rho`$ and $`g\in C(\mathcal A)`$ set

```math
U_s g(S)=\int_0^\infty e^{-st}R_t[b_{(\cdot)}g(p_{(\cdot)})](S)\,dt,\qquad (K_s g)(a)=U_sg(J_a). \qquad\text{(7)}
```

This gives the mean discounted offspring produced by a residual family when direct births of new $`J`$ particles are counted and set aside. Equation (5) gives $`|U_sg|\le W\|g\|_\infty/(s-\rho)`$. The kernel in (7) can be constructed by positive operations on finite types.

The only dependence on $`a`$ in the original particle $`J_a`$ is its constant death rate. In particular the law of its intrinsic tree and the marks of every child it ever copies are independent of $`a`$. Separate that particle from the independent residual families it emits. Its contribution to (7), together with the future contributions of those copies, is the integral against $`e^{-(s+a)t}`$ of the intrinsic expectation of

```math
b_Sg(p_S)+\sum_{Y\text{ at }S}[(a_Y-p_S)U_sg(Y)+p_S U_sg(Y^+)]. \qquad\text{(8)}
```

For $`\|g\|_\infty\le1`$ the absolute intrinsic expectation is at most $`1+e^{qt}/(s-\rho)`$, by the child-weight identity and $`L^0W\le qW`$, starting at a newborn. Therefore all the rows in (7) depend Lipschitz-continuously on $`a`$ in total variation, uniformly for $`s`$ in compact subintervals of $`(\rho,\infty)`$. Arzelà–Ascoli gives compactness on $`C(\mathcal A)`$. The exponential bound also gives norm continuity in $`s`$.

Every row of $`K_s`$ charges every $`1/j`$, $`j\ge2`$. The original residual particle can survive and its initial root can acquire leaves up to degree $`j`$, keeping its first leaf and experiencing no promotions, then produce a direct child. Each prescribed finite sequence has positive probability. Thus $`K_s`$ is strongly positive, since a nonzero continuous nonnegative $`g`$ is positive at some $`1/j`$, and the output is continuous on a compact set. The strong Krein–Rutman theorem gives a simple positive spectral radius with a continuous eigenfunction having positive minimum and a positive dual eigenmeasure. The radius decreases strictly with $`s`$: for $`s_2>s_1`$ the difference $`(K_{s_1}-K_{s_2})h_2`$ on the positive eigenfunction of $`K_{s_2}`$ is everywhere strictly positive, hence bounded below by a positive multiple of $`h_2`$. Simplicity or this same comparison gives continuity of the radius. [Positive compact operators](https://www.aimsciences.org/article/doi/10.3934/dcds.2020054)

For clarity, let $`g_s(S)=\int_0^\infty e^{-st}R_t1(S)\,dt`$ and $`u_s(S)=\int_0^\infty e^{-st}P_t1(S)\,dt`$, allowing infinity for the latter. Group by the generation of direct $`J`$ births. At newborns the sum is

```math
u_s(J_{(\cdot)})=\sum_{n\ge0} K_s^n[g_s(J_{(\cdot)})]. \qquad\text{(9)}
```

The forcing on newborns is bounded and at least $`(s+1/2)^{-1}`$ by survival of the first particle. It is continuous by the calculation (8), with the leading term changed from $`b_Sg(p_S)`$ to 1 and with $`g_s`$ on copied types. If $`r(K_s)<1`$, (9) converges in norm and $`u_s(S)=g_s(S)+U_s[u_s(J_{(\cdot)})](S)\le C_s W(S)`$. Inserting this in (3) and using (6) makes $`\int_0^\infty e^{-st}M_t1\,dt`$ finite. If $`r(K_s)\ge1`$, the positive forcing dominates a multiple of the positive eigenfunction, and (9) diverges at every newborn state. The centre has positive mean rate of direct $`J`$ immigration over a time interval, so the same integral of $`M_t1`$ diverges by (2),(3) and Tonelli. Now $`M_t1=\mathbb E_{T_3}C(t)-1`$. Its logarithmic rate is $`\mu>\rho`$, so its Laplace integral converges for $`s>\mu`$ and diverges for $`\rho<s<\mu`$. Strict monotonicity and continuity identify $`r(K_\mu)=1`$ and give uniqueness.

## 5. Constant-factor growth of the ordinary mean

Let $`g>0`$ satisfy $`K_\mu g=g`$ and put $`H=U_\mu g`$. This is bounded above by a multiple of $`W`$, and $`H(J_a)=g(a)`$. The resolvent identity for $`R_t`$ gives

```math
R_tH=e^{\mu t}H-\int_0^t R_s B(e^{\mu(t-s)}H)\,ds.
```

Decomposition by the first generation of direct births gives $`P_t H=R_tH+\int_0^tR_s B P_{t-s}H\,ds`$. The Volterra equation has at most one solution locally bounded in the weighted supremum norm, by (5), boundedness of $`B`$ in this norm and iteration of the integral inequality. Consequently $`P_tH=e^{\mu t}H`$. From (3) and (6),

```math
e^{-\mu t} M_tH=\int_0^t e^{-\mu s} I_s H\,ds\ \longrightarrow\ A_0\in(0,\infty). \qquad\text{(10)}
```

Positivity follows from $`g`$ having positive minimum and the initial centre promoting at a positive rate before its first promotion.

Here is the comparison to the actual count. Every core can contribute to at most one ancestral $`W`$ at each distance. Since every birth mark satisfies $`z(a)\le z(1/2)<1`$, summing the geometric series gives $`\sum_{v\ne o}H(S_v)\le C_1 C`$ on every tree. In the reverse direction, a noninitial vertex with $`c\le2`$ and at least one leaf has $`H(S_v)\ge c_1>0`$ uniformly. In the intrinsic original residual particle, insist it survives one unit (probability at least $`e^{-1/2}`$). Until its first root clock ring its leaf count cannot fall, so at that ring it promotes with probability at least $`1/3`$. The expected discounted direct-birth reward through that unit is at least $`e^{-(\mu+1/2)}(1-e^{-1})\min g/3`$. This use of the intrinsic promotion is legitimate since (7) integrates its identical predictable rate for direct births. A bare noninitial vertex ($`c=D=2`$) having a bare core *child* also has a uniform positive lower bound: silence that root and the bare child's other neighbour for one unit, and let the bare child's first ring select the root. This gives the root a leaf; also require survival of the residual particle. Apply the preceding bound in the next unit. Rings sent into the two silent vertices do not change the bare child's degree before its first choice.

At least $`C/4`$ vertices satisfy the eligibility lemma of the existence proof. Among them, a bare vertex with a bare neighbour but no bare child has a bare parent with a bare child. There are at most two such vertices per such parent (only the initial centre can have two core children while bare). The initial centre costs at most one from the resulting count. Thus the number of noninitial vertices of the two good kinds just described is at least $`C/12-1`$. This proves $`C\le C_2(1+\sum_{v\ne o}H(S_v))`$. Together with (10) it gives $`\mathbb E_{T_3}C(t)\asymp e^{\mu t}`$; the independent first wait of rate two gives the same comparison starting from the edge.

At fixed size $`n\ge3`$, $`a_j`$ is nondecreasing and $`a_m\le(m/n)a_n`$ for $`m\ge n`$ by the promotion probability bound $`C_j/j`$. Write $`N(t)`$ for the Yule size started from the edge, independent of the jump chain. At $`t=\log n`$, $`\Pr(N(t)\ge n)\ge e^{-1}`$ and $`\mathbb E N(t)=2n`$. Hence

```math
e^{-1}a_n\le\mathbb E a_{N(\log n)}\le a_n\mathbb E\max(1,N(\log n)/n)\le3a_n.
```

This transfers the constant-factor estimate to $`a_n`$.
