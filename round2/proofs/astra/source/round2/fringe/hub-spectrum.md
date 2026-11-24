# From fixed degrees to ranked hub masses

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 11: HY4 Flash for ranked/intermediate limits; GPT for the fixed-hub extensions; inputs include Claims 3, 6 and 7. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


Let $`X=\lim C_n/n^\mu`$, let $`p_k`$ be the limiting fraction of cores of degree $`k`$, and set $`P_k=\sum_{j\ge k}p_j`$. Write $`A=\lim_{k\to\infty}k^\mu P_k>0`$. These limits and $`p_k\sim\mu A k^{-1-\mu}`$ come from the [mean and sample-path results](almost-sure.md) and the [degree tail](../degree-law.md). Let $`Z_v=\lim_n D_v(n)/n`$ for each fixed label, with $`\sum_v Z_v=1`$ almost surely as in [the hub mass calculation](../hub-mass.md).

The proposed connection is the almost-sure limit
```math
 \lim_{x\downarrow0}x^\mu\#\{v:Z_v\ge x\}=A X.                 \qquad\text{(1)}
```
Equivalently the $`j`$th largest limiting degree fraction is asymptotic to $`(AX/j)^{1/\mu}`$. For every range $`k\to\infty`$, $`k/n\to0`$, the number of degrees at least $`k`$ is asymptotic to $`AX(n/k)^\mu`$. The probability-one event is common to all such choices. There is also a deterministic $`\gamma>0`$ such that
```math
 \sup_{2\le k\le n^\gamma} k^{1+\mu}
 \left|N_{=k}(n)/n^\mu-Xp_k\right|\longrightarrow0
 \quad\hbox{almost surely}.                                      \qquad\text{(2)}
```

## 1. A uniform estimate at moderate degrees

Work in continuous time from $`T_3`$. Use $`d_k=d(\mathbf1_{D\ge k})`$, the physical eigenfunction $`\Psi`$, and $`Z^*_t=e^{-\mu t}\Psi(T_t)`$ from (3)–(9) of [almost-sure.md](almost-sure.md). Thus $`P_k=d_k/d_2`$, $`\sup_k k^\mu d_k<\infty`$, $`\sup_k k^{1+\mu}(d_k-d_{k+1})<\infty`$, and $`Z^*_t\to Z^*_\infty>0`$. With $`U_k(t)=N_{\ge k}(T_t)`$, all three bounds below have constants independent of $`k`$. Conditional on the tree at $`u`$, for $`h\ge0`$,
```math
 \begin{aligned}
 |\mathbb E[U_k(u+h)\mid\mathcal F_u]-e^{\mu h}\Psi(T_u)d_k|
     &\le K C(T_u)e^{\sigma h},\quad \sigma<\mu,\\ \mathop{\mathrm{Var}}\nolimits (U_k(u+h)\mid\mathcal F_u)&\le K n(T_u)e^{4h},\\ \mathbb E n(T_u)&=3e^u .
 \end{aligned}                                                     \qquad\text{(3)}
```
Also $`e^{-\mu u}C(T_u)`$ stays bounded almost surely. The uniformity in (3) follows from the weighted supremum estimate for bounded fringe scores and the single-leaf coupling for all degree thresholds.

Put $`\eta=\mu-\sigma`$, $`a=(2\mu-1)/8`$, $`g=\eta(2\mu-1)/100`$, and $`b=2g`$. For each integer $`j`$ condition at the *same* time $`u=(1-a)j`$ for the whole interval $`[j,j+1]`$. Divide that interval into $`\lceil e^{bj}\rceil`$ equal subintervals. At a mesh point $`v`$ the centered variable $`e^{-\mu v}(U_k(v)-\mathbb E[U_k(v)\mid\mathcal F_u])`$ has second moment at most $`K'e^{-(2\mu-1)j+4aj}`$. Chebyshev's inequality, summed over mesh points and $`2\le k\le e^{g(j+1)}+1`$ with squared weight $`k^{2(1+\mu)}`$, gives a summable series since $`4a+b+(3+2\mu)g<2\mu-1`$. By Borel–Cantelli the maximum weighted centered error tends to zero. The first line of (3) handles the conditional means since $`(1+\mu)g<\eta a`$. Consequently, at all these mesh points simultaneously,
```math
 U_k(v)= e^{\mu v}\big[d_k Z^*_{(1-a)j}+o(k^{-1-\mu})\big],      \qquad\text{(4)}
```
where the little-oh is uniform over the indicated $`k`$. Each integer block has its own common $`Z^*_{(1-a)j}`$.

For $`t`$ between successive points $`v,w`$, use $`U_k(v)\le U_k(t)\le U_k(w)`$. For the exact-degree count use
```math
U_k(v)-U_{k+1}(w)\le N_{=k}(T_t)\le U_k(w)-U_{k+1}(v).
```
After (4), the extra error in the latter two leading coefficients, multiplied by $`k^{1+\mu}`$, is at most $`K k(w-v)Z^*_{(1-a)j}=o(1)`$ uniformly. This use of the same conditioning time avoids differencing two uncontrolled martingale errors. Replacing the common martingale value by its limit now costs $`o(1)`$ after scaling by $`k^\mu`$ for $`U_k`$ or $`k^{1+\mu}`$ for $`N_{=k}`$. This gives uniform convergence for both scores throughout $`2\le k\le e^{gt}`$.

At the size-$`n`$ epoch $`\tau_n`$, $`n e^{-\tau_n}\to Y>0`$ and $`X=d_2 Z^*_\infty/Y^\mu`$. Taking any $`0<\gamma<g`$ gives (2) and its version for $`N_{\ge k}`$ with weight $`k^\mu`$ and limit $`X P_k`$.

## 2. A vertex's limiting mass as seen at time $`t`$

For a vertex already present, put $`M_v=\lim_{h\to\infty}e^{-h}D_v(t+h)`$, conditionally on the tree at $`t`$. This limit is $`e^t Y Z_v`$. The same argument applies when $`t`$ is a size epoch. From a current degree $`D`$, $`M_v`$ is dominated by a gamma variable of shape $`D`$ and unit scale, by comparison with a pure-birth process of rate its own population size. This includes vertices which are still leaves.

For a current core having $`L=D-c`$ adjacent leaves, there is also a coupling with $`M_v\ge G_L-E`$, where conditionally $`G_L`$ has the gamma law of shape $`L`$ and unit scale and $`E`$ is exponential of mean one. Set $`G_0=0`$. Independence of these two variables is unnecessary. Tag all $`L`$ current leaves. Whenever a tagged leaf rings, it creates another tagged leaf at the focal core. A ring at the focal core may promote a tagged leaf, in which case it moves to an auxiliary population and starts a fresh independent rate-one Yule family there. Ignore all incoming untagged leaves. The total tagged population, counting its auxiliary members, has single births at rate its current size and loses no members. Thus its normalized limit has the stated law $`G_L`$. To bound the missing population, start an auxiliary Yule family at *every* rate-one focal ring, using that family's future in the missing population if a tag was promoted at that ring. The resulting dominating population has immigration rate one and unit birth rate per member. Its size plus one is a Yule process starting at one, so its normalized limit is $`E`$. The physical leaf count dominates the remaining real tags. The construction uses fresh independent family drivers at each focal ring.

For every fixed $`0<\epsilon<1/4`$ and integer threshold $`k`$ sufficiently large, these comparisons have the following conditional consequences with constants depending only on $`\epsilon`$:
```math
 \begin{aligned}
 D\le(1-\epsilon)k &\implies \Pr(M_v\ge k\mid\mathcal F_t)\le K e^{-c_\epsilon k},\\ D\ge(1+\epsilon)k,\ c\le\epsilon k/4
   &\implies \Pr(M_v< k\mid\mathcal F_t)\le K e^{-c_\epsilon k}.
 \end{aligned}                                                       \qquad\text{(5)}
```
For the second line take $`E\le\epsilon k/4`$, use $`L\ge(1+3\epsilon/4)k`$, and apply an exponential lower-tail bound to $`G_L`$. Both gamma estimates follow directly from their moment generating functions. The same calculations show that if $`c\le\epsilon D/4`$ the probability of $`|M_v/D-1|>\epsilon`$ is at most $`K e^{-c_\epsilon D}`$.

A vertex first born at a subsequent time $`r`$ starts at degree one. Its contribution to the number of future vertices with $`e^tYZ_v\ge k`$ is bounded in conditional mean by $`\exp[-k e^{r-t}]`$ at its birth. The conditional expected total birth rate at $`r`$ is $`n(T_t)e^{r-t}`$. Compensation gives the bound
```math
 n(T_t)\int_0^\infty e^h \exp(-k e^h)\,dh= n(T_t)e^{-k}/k .        \qquad\text{(6)}
```
These estimates use no independence between the futures of different vertices.

## 3. Counting small limiting masses

The maximum core degree at time $`T`$ is $`O(T)`$ almost surely as $`T\to\infty`$. Indeed, from its arrival every new vertex has at most one old neighbour and can promote another neighbour only on its own rate-one initiating clock. The number of its initiations over a time interval of length at most $`T`$ is bounded in law conditionally on arrival by $`\mathop{\mathrm{Poisson}}\nolimits (T)`$. Summing exponential moments over the initial vertices and expected births through $`T`$ bounds the expectation of the sum of $`e^{c_v(T)}`$ by $`K e^{eT}`$. Markov's inequality and Borel–Cantelli at integer times, then monotonicity, give the asserted bound.

Put $`Z_v^{\rm CT}=YZ_v`$. On any fixed mesh $`t=j/m`$ observe the physical tree at $`T=t+2\log(t+2)`$, so $`k=e^{T-t}=(t+2)^2`$. The preceding maximum bound gives $`c_v(T)\le\epsilon k/4`$ eventually, for each fixed $`\epsilon`$. On the event enforcing that bound, (5) and a union bound over old vertices give failure probability at most $`K\mathbb E n(T_T)e^{-c_\epsilon k}`$, where $`n(T_T)`$ means the population at observation time $`T`$. Equation (6) treats future vertices. The two probabilities are summable on the mesh. Integer cutoffs in the resulting degree counts can be taken by floor or ceiling with an error of at most one in the threshold. Almost surely, eventually,
```math
 U_{\lceil(1+\epsilon)k\rceil}(T)
 \le \#\{v: Z_v^{\rm CT}\ge e^{-t}\}
 \le U_{\lfloor(1-\epsilon)k\rfloor}(T).
```
The degree thresholds are polynomial in $`T`$, so Section 1 applies. Since $`e^{\mu T}k^{-\mu}=e^{\mu t}`$ and $`P_j\sim A j^{-\mu}`$, the limiting lower and upper ratios to $`e^{\mu t}`$ are squeezed between $`A d_2Z^*_\infty(1+\epsilon)^{-\mu}`$ and $`A d_2Z^*_\infty(1-\epsilon)^{-\mu}`$. Send rational $`\epsilon\downarrow0`$ and intersect all meshes. The count is nondecreasing in $`t`$, so consecutive mesh points bound its ratios everywhere with factors $`e^{\pm\mu/m}`$. This proves the full continuous limit. Substituting $`e^{-t}=Y x`$ and using $`d_2 Z^*_\infty/Y^\mu=X`$ gives (1). Generalized inversion of the monotone count gives the ranked statement.

## 4. All intermediate cumulative thresholds

For clarity take the threshold $`l_n=\lceil(\log(n+1))^2\rceil`$ at size $`n`$. At $`t=\tau_n`$ we have maximum core degree $`O(\log n)`$ almost surely. By the last estimate after (5), uniformly among all existing vertices with $`D_v\ge l_n/2`$, $`e^{\tau_n}YZ_v/D_v\to1`$ almost surely. Indeed for each fixed rational $`\epsilon>0`$ the failure probability restricted to $`\max c_v\le\epsilon l_n/8`$ is at most $`Kn\exp(-c_\epsilon l_n/2)`$, summable over $`n`$. For an existing vertex with $`D_v<l_n/2`$, the probability that $`e^{\tau_n}Y Z_v\ge l_n`$ obeys the same summable union bound. For vertices born after size $`n`$, apply (6) at $`\tau_n`$. Finally $`e^{\tau_n}Y/n\to1`$.

For thresholds $`k\ge 2l_n`$, the high-degree comparisons thus sandwich $`N_{\ge k}(n)`$ between the counts of limiting masses above $`(1+o(1))k/n`$ and $`(1-o(1))k/n`$, with uniform perturbations. Equation (1) gives $`N_{\ge k}(n)\sim AX(n/k)^\mu`$ throughout $`k\ge2l_n`$, $`k/n\to0`$. For smaller $`k\to\infty`$ Section 1 and the limiting cumulative tail give the same result. In particular the large-$`k`$ and large-$`n`$ statements hold on a single event for every intermediate sequence. As another consequence of (1), for $`s>0`$ the sum $`\sum_v Z_v^s`$ is finite exactly when $`s>\mu`$, and $`(s-\mu)\sum_v Z_v^s\to\mu AX`$ as $`s\downarrow\mu`$, by integrating the count in (1).
