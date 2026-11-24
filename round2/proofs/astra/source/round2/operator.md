# A first-moment operator on rooted core trees

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../SOURCES.tsv)).
> **Review and dependencies:** Claim 3: Claude and Gemini; GPT standalone audit; inputs include Round 1 growth, the marked renewal and spatial weights (Claims 2 and 4). Claim 4: Claude and Gemini; inputs include Round 1 one-edge coupling, planting and growth. See the [Astra ledger](../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Statement under review

Let $`\mathcal X`$ consist of isomorphism classes of finite or countably infinite, connected, locally finite rooted trees of core vertices. Give each vertex a nonnegative integer mark $`L`$ for its number of leaves, with $`D=L+c\ge2`$, where $`c`$ is its degree in the core. Use the topology of agreement on marked balls. A one-vertex core with $`c=0`$ is allowed. No orientation towards the initial seed is part of a state.

The construction below specifies a bounded positive operator $`Q=Q_1`$ without using the unknown exponent. Consider positive eigenvalues of its action on finite measures which have a probability eigenmeasure. The proposed characterization is

```math
e^\mu=\min\{r>0:\ \exists\text{ a Borel probability }\pi\text{ on }\mathcal X,\quad \pi Q=r\pi\}. \qquad\text{(1)}
```

In fact a single probability in (1) satisfies $`\pi Q_s=e^{\mu s}\pi`$ for every $`s\ge0`$. For this law $`\int L_o/D_o\,d\pi=\mu`$. No choice of boundary condition at infinity enters $`Q_s`$. The spatial cutoff estimate below applies at fixed time. Equation (1) does not itself supply numerical digits, a simple eigenvalue, convergence to the eigenmeasure from every seed, or a positive right eigenfunction at the critical parameter.

## 1. Cells and evaluation of the operator

At time zero give each old core its own label. A newly promoted core inherits the label of the vertex which promoted it. The population with label $`o`$ is $`H_o(s)`$. Old labels partition the future core. Old core edges are the only edges between different labels, because every new core is born adjacent to its promoter. Define

```math
(Q_s f)(G,o)=\mathbb E_{G,o}\sum_{v\in H_o(s)} f(G_s,v). \qquad\text{(2)}
```

All rates are specified by the marked tree. A vertex $`x`$ gains a leaf at rate $`r_x=L_x+\sum_{y\sim x}1/D_y`$. It promotes one leaf at rate $`b_x=L_x/D_x`$, changes $`(L_x,c_x)`$ to $`(L_x-1,c_x+1)`$, and creates a new adjacent core in state $`(1,1)`$. A gain changes $`L_x`$ to $`L_x+1`$. Each member promotes at rate at most one. The size of one cell is dominated by a rate-one Yule population begun with one individual. Consequently $`Q_s1\le e^s`$ and $`\mathbb E|H_o(s)|^2\le 2e^{2s}`$. Relabelling at time $`s`$ and partitioning descendants gives $`Q_sQ_t=Q_{s+t}`$ by conditional expectation. Independence of the cells is unnecessary.

Here are details for infinite initial cores. Write $`c_x^0`$ for the initial core degree. For every oriented *old* edge $`x\to y`$, use an independent Poisson process of rate $`1/c_x^0`$. At a point accept with probability $`c_x^0/D_x`$, using a fresh uniform mark, and on acceptance give $`y`$ a leaf. These rates are admissible since $`D_x\ge c_x^0`$. All other transitions act within one cell. They use independent cell drivers: rate-$`L`$ leaf births, promotions, and sends along the edges internal to the cell. These latter sends may be implemented with a rate-one clock at each core and subintervals of lengths $`1/D`$ for each internal neighbour; allocate a further subinterval of length $`L/D`$ for promotion and leave the old-edge slots idle. This gives precisely the rates above. A finite union of cells is nonexplosive in bounded time. The number of incoming proposals from its finite old boundary is bounded by a finite-rate Poisson process, and all endogenous births are bounded in rate by the current population, counting the actual leaves.

Compute inside the old ball of radius $`R`$ about $`o`$ with all missing old neighbour slots retained and all inbound arrows across its boundary suppressed. Use the same drivers on the common cells whenever two such truncations are compared. If an exterior change has reached cell $`w`$ in the interior by time $`s`$, a time-ordered sequence of proposals has crossed the old path from the boundary to $`w`$. Indeed a cell depends on another cell only via an old-edge proposal. A repeated vertex in a traversal can be erased, keeping an increasing sequence of times. For a specified simple oriented path of $`j`$ edges, the expected number of ordered proposal tuples is

```math
\frac{s^j}{j!}\prod_{x\text{ sending on the path}}\frac1{c_x^0}.
```

For paths of length $`j`$ ending at $`w`$, the sum of these products is at most $`c_w^0`$. To see this, reverse each path, use $`c_x^0 P(x,y)=c_y^0 P(y,x)=1`$ for simple random walk, and use $`c_x^0\ge1`$ at the original start. The sum of probabilities of the reversed simple paths is at most one. For a singleton there is no path to consider.

If $`f`$ is bounded and reads a marked ball of radius $`k\le R`$, the final balls of the members of $`H_o(s)`$ can meet only cells whose bases are in the old $`B_k(o)`$. With $`j=R-k+1`$, put

```math
p=\min\left(1,\,\frac{s^j}{j!}\sum_{w\in B_k(o)}c_w^0\right).
```

A union bound gives failure probability at most $`p`$ for agreement on these cells throughout $`[0,s]`$. On failure use Cauchy–Schwarz separately for the two Yule-dominated cell sizes. Thus

```math
|Q_s f-Q_s^{[R]}f|(G,o)\le 2\sqrt2\,e^s\|f\|_\infty\sqrt p. \qquad\text{(3)}
```

The ordered-proposal event, which depends only on the prescribed old-edge drivers, simultaneously controls all larger balls. Its probabilities are summable in $`R`$ for fixed $`s,k,G,o`$. The dynamics therefore stabilizes almost surely on each fixed set of cells over each bounded time interval. This constructs the process and the Markov property using independent future driver increments. Formula (3) gives continuity of $`Q_sf`$ on $`\mathcal X`$ for each bounded cylinder $`f`$: the factor from $`B_k(o)`$ is fixed in a neighbourhood of a given input, and each truncated answer depends only on a finite marked old ball with the boundary degrees. These cylinders determine Borel probabilities. The same argument with Yule uniform integrability and local convergence of the coupled output gives continuity for every bounded continuous $`f`$.

## 2. A probability eigenmeasure at the growth rate

Start the finite chain at the centre and two leaves, at time zero. Write $`m_t(f)=\mathbb E\sum_{u\in K_t}f(G_t,u)`$ and $`F_3(t)=m_t(1)`$. Equation (2) gives $`m_{t+h}=m_tQ_h`$. The result for an initial edge in the repository gives $`\log F_3(t)/t\to\mu`$. Indeed the first birth from the edge has an independent exponential wait of rate two. If $`F`$ is the mean from the edge, monotonicity gives $`F(t)\le F_3(t)\le F(t+1)/(1-e^{-2})`$. [Growth from finite seeds](../proofs/existence.md)

For $`0<\delta<\mu`$ and $`T>0`$ let $`w(t)=\exp(-\mu t-\delta|t-T|)`$, $`Z=\int_0^\infty w(t)F_3(t)\,dt`$ and $`\pi_{\delta,T}=Z^{-1}\int_0^\infty w(t)m_t\,dt`$. The integral is finite by the logarithmic growth rate. For fixed $`\delta`$, $`w(0)/Z\to0`$ as $`T\to\infty`$: integration on $`[T,T+1]`$ gives $`Z\ge e^{-\mu(T+1)-\delta}F_3(T)`$, whereas $`w(0)=e^{-\delta T}`$. Take $`\delta\downarrow0`$, choosing $`T`$ so that $`w(0)/Z\to0`$. For fixed $`h`$, $`w(u-h)/w(u)`$ lies between $`e^{(\mu-\delta)h}`$ and $`e^{(\mu+\delta)h}`$ on $`u\ge h`$. The lost contribution on $`[0,h]`$ is $`O_h(w(0)/Z)`$ since $`F_3(t)\le e^t`$. Changing variables in the integral yields

```math
\|\pi_{\delta,T}Q_h-e^{\mu h}\pi_{\delta,T}\|_{\rm TV}\longrightarrow 0. \qquad\text{(4)}
```

The measures cannot lose mass at infinite degree. At the promotion time of a vertex its total degree is two. After an elapsed time $`a`$ its expected degree is at most $`2e^a`$, since its gain rate is at most its degree. Its core degree is dominated by $`1+\mathop{\mathrm{Pois}}\nolimits (a)`$ using its own promotion clock. Both estimates are conditional on the full tree at promotion. They also bound the initial centre. Hence, for $`0<\alpha<\min(1,\mu)`$ and $`0<\theta<\log(1+\mu)`$, the conditional means of $`D^\alpha`$ and $`e^{\theta c}`$ are bounded respectively by $`2^\alpha e^{\alpha a}`$ and $`e^\theta e^{(e^\theta-1)a}`$. The mean total promotion intensity at time $`t`$ is at most $`F_3(t)`$. If one of these moment bounds is $`Be^{\kappa a}`$, summing over births, integrating and using $`w(t+a)\le w(t)e^{-(\mu-\delta)a}`$ gives

```math
\pi_{\delta,T}(\text{moment})\le\frac{B(1+w(0)/Z)}{\mu-\delta-\kappa}\qquad(\delta<\mu-\kappa). \qquad\text{(5)}
```

The root-uniform measures obey the mass-transport identity, by interchanging the two vertex sums on every finite tree. Here are details of ball tightness. Put $`E_J=\{D_o>J\}`$. Equation (5) and Cauchy–Schwarz give $`\sup\pi_{\delta,T}(c_o\mathbf1_{E_J})\to0`$ along the chosen sequence after omitting finitely many terms. To make all marks in a fixed ball tight, first bound the root degree with probability $`1-\varepsilon`$. If the earlier vertices along a simple path of length $`i`$ have degrees bounded by chosen cutoffs $`J_0,\ldots,J_{i-1}`$, the expected number of such paths from the root ending at degree $`>J_i`$ is at most

```math
(J_1\cdots J_{i-1})\,\pi_{\delta,T}(c_o\mathbf1_{D_o>J_i}).
```

Reverse the path by mass transport. The first step out of the bad endpoint has $`c_o`$ choices and every further step has the relevant previous cutoff. Taking $`J_i`$ large in sequence bounds the failure probability on a path of each length $`i`$. There are only finitely many marked balls at fixed radius with these cutoffs. Choosing summable failure probabilities over all radii gives tightness in the local topology; this space of locally finite rooted, discretely marked trees is Polish (compatible balls in a Cauchy sequence specify the limiting rooted tree). Extract a weak limit $`\pi`$. For every fixed $`h`$, test (4) on bounded cylinders, use the continuity from (3) and then the determining property. The same subsequence works for every $`h`$, since (4) held on the entire sequence for each $`h`$. This gives $`\pi Q_h=e^{\mu h}\pi`$.

## 3. Minimality

For each fixed $`(G,o)`$ one has $`\liminf_{s\to\infty}s^{-1}\log Q_s1(G,o)\ge\mu`$. The needed input from the finite-seed argument is specific. Given $`M\ge2`$ and $`0<\epsilon<1`$, it supplies a fixed $`k`$ and for each sufficiently large $`s`$ a rooted $`k`$-vertex seed such that a copy attached by one cut edge to a parent of degree at least $`M`$ contributes, in conditional mean, at least

```math
\frac{e^{-\eta s}}{k}\big(F(s)-k-\theta_k(s)\big),\qquad \eta=7\log2/\log(1/\epsilon)+9/M,\quad \theta_k(s)=O_k(e^{-s}). \qquad\text{(6)}
```

This uses the weighted one-edge coupling and the choice among finitely many roots in Corollary 5.3 and Lemma 6.1 of the existence proof. For a locally finite exterior fix the elapsed time, the planted seed, and the attaching edge. Complete successively larger old balls to finite core trees, attaching a core with one leaf to each missing core-neighbour slot to preserve the boundary degrees. The finite one-edge estimate is uniform in these exteriors. The ordered-proposal estimate in Section 1 couples the finite completions to the prescribed initial tree on the planted cells over the fixed elapsed time, with disagreement probability tending to zero. Their cell counts are uniformly integrable by the Yule domination for the finitely many planted labels. Taking the limit of the expected count transfers (6). This argument can be applied at the end of the planting script, conditional on its realized tree. [One-edge estimate](../proofs/existence.md)

The cell of $`o`$ can plant such a copy with a fixed positive probability, uniform in the choice among the finitely many seeds once $`(G,o),M,k`$ are fixed. If $`o`$ starts without a leaf, choose one of its neighbours $`w`$, silence all neighbours of $`w`$ for one time unit (including $`o`$), and let $`w`$ ring with its first choice at $`o`$. The original degree of $`w`$ stays fixed until that choice; all these conditions have positive probability. If $`o`$ has a leaf, simply silence $`o`$ for this unit. In a further unit silence $`o`$ and require a tagged leaf of $`o`$ to ring at least $`M`$ times. This guarantees $`D_o\ge M`$ with $`c_o`$ still the given finite number. In the next unit require one ring of $`o`$ choosing a leaf and protect that leaf and its newborn until the unit ends. The promotion probability at that ring is bounded below by $`(M+1)/(M+1+c_o)`$. Silence $`o`$ and script the other $`k-2`$ births, one per unit in parent-before-child order, exactly as in Lemma 3.1 of the existence proof. The cells of all newly promoted cores here descend from $`o`$. Conditions on finitely many clocks have the same probabilities in every larger-ball construction. They give a delay $`h`$ and a probability $`p(G,o,M,k)>0`$ independent of $`s`$. Multiplying (6) by $`p`$ bounds $`Q_{s+h}1(G,o)`$ below. First send $`s\to\infty`$, then $`\epsilon\downarrow0`$ and $`M\to\infty`$.

If $`\sigma Q=r\sigma`$ for a probability $`\sigma`$, (2) gives $`\int Q_n1\,d\sigma=r^n`$ at every integer time $`n`$. The preceding pointwise lower estimate and Fatou's lemma applied to $`e^{-(\mu-a)n}Q_n1`$ for $`a>0`$ imply $`r\ge e^\mu`$. Section 2 supplies equality, giving (1).

## 4. Drift equations and the two edge transfers

For a bounded cylinder $`f`$ write $`\Delta_x^+f`$ and $`\Delta_x^bf`$ for the changes at the old root when $`x`$ gains or promotes. The operator on test functions at time zero is

```math
Af(G,o)=\sum_x [r_x\Delta_x^+f+b_x\Delta_x^bf]+b_o f(G^{b_o},\mathrm{new}). \qquad\text{(7)}
```

Only finitely many old vertices affect a cylinder on any fixed tree. The extra term gives the newborn to its promoter's cell. For a cylinder supported on balls with bounded total degrees throughout the ball, the difference terms in (7) are uniformly bounded. A gain or promotion cannot remove a large-degree old vertex from the ball or reduce its total degree. Including the constant cylinder gives $`A1=b_o`$. Dynkin's formula for the sum over a cell is $`Q_t f-f=\int_0^t Q_sAf\,ds`$ for these tests. One can first stop the cell at a size bound and use (7) on the finite truncations, then remove the size bound using the Yule moments and (3). Thus $`\pi(Af)=\mu\pi(f)`$, in particular $`\pi(b_o)=\mu`$.

For the gain transfer in the degree-band certificate, the directed arc $`(w,u)`$ carries $`r_w g(u,w)`$ from row $`u`$ to row $`w`$. The first argument of $`g`$ is the focal vertex whose row is debited. For the second transfer it carries $`\mathbf1_{D_w\ge T}a(D_u,c_u)/D_w`$ from row $`u`$ to row $`w`$. Their divergences integrate to zero under $`\pi`$ by the same mass transport as on a finite tree. The limit is unimodular: first pass bounded transports supported on bounded-degree balls to the weak limit, then exhaust the balls and use monotone convergence for nonnegative transports. For the signed tables use their finite support and (5). These transfers change how rows of (7) are bounded. They need no additional assumption about independence of the neighbouring states. [Degree-band transfers](../proofs/upper-bound-banded.md)

## 5. Spectral radius from the uniform one-edge estimate

The separate [one-edge upper bound](finite/one-port-upper.md), equations (1)–(2), upgrades (1) to an equality with the full spectral radius on $`C_b(\mathcal X)`$. For any finite starting core and any one old label $`o`$, all of its new members are in the disjoint downward branches founded at direct promotions of $`o`$. The original member stays in the cell. Its direct promotion intensity is at most one. At each such birth the conditional mean size of the branch over an additional $`u`$ time units is bounded by $`A_j(u)e^{\Lambda_j u}F(u+1)`$, uniformly over the entire physical exterior. This is the branch estimate, before the sum over initial cores in the cited proof. Therefore
```math
Q_t1(G,o)\le 1+t A_j(t)e^{\Lambda_j t}F(t+1).
\qquad\text{(8)}
```
The bound applies in the infinite space as well. Complete a finite old ball to a finite genuine core tree with the correct degrees on the ball, by appending one new core with one leaf to each unrepresented core-neighbour slot on its boundary. The value at the original root of $`Q_t1`$ for these finite completions tends to its infinite-tree value by (3). All the constants in (8) are independent of the completion.

Positivity gives $`\|Q_n\|=\|Q_n1\|_\infty`$ on $`C_b(\mathcal X)`$, and the semigroup gives $`Q_n=Q_1^n`$. For each fixed $`j`$, (8), $`\log F(t)/t\to\mu`$ and $`\log A_j(t)=O_j(\log(2+t))`$ give $`\limsup n^{-1}\log\|Q_1^n\|\le\mu+\Lambda_j`$. Let $`j\to\infty`$. Conversely the probability eigenmeasure in Section 2 gives $`\|Q_1^n\|\ge e^{\mu n}`$. The spectral-radius formula now yields $`r(Q_1)=e^\mu`$. A probability eigenmeasure of $`Q_1`$ at any positive eigenvalue must have this eigenvalue, by (1) and the spectral-radius upper bound.
