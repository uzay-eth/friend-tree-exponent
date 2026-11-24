> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 3 and 10; Gemini reviewed the base sequential law, with separate scopes for addenda. See the [release index](../../../../README.md#claim-20) for exact scope and citations.

# A uniformly counted component below a large degree cutoff

Delete all vertices of the original friend tree whose **original** degree exceeds $`N`$. Among the remaining components containing at least one original nonleaf, choose a component uniformly and root it at its highest surviving core, in the orientation toward the founder. Keep its ordinary vertices and edges and forget the removed parent. An isolated surviving core counts as one component. This sampling rule answers a version of the question about young friend trees in §7(5) of Addario-Berry et al. [Original question](https://arxiv.org/html/2403.20185v1#S7)

Let $`C_n`$ count the nonleaves at size $`n`$. Throughout, take $`n\to\infty`$ at fixed $`N`$, then $`N\to\infty`$. Write $`d_N`$ for the limiting number of the specified components per $`C_n`$, and $`\nu_N`$ for their limiting probability law on finite rooted trees. Let $`J`$ be a root with one real leaf and one inert parent slot. Inside this planted tree every real vertex rings at rate one and chooses a neighbour uniformly. A root choice of the slot is null. Let $`P_a^0(J,\cdot)`$ return the ordinary rooted tree at elapsed time $`a`$, omitting the slot from the output. Write $`A`$ for the positive constant in the limiting cumulative degree tail $`\pi(D_o\ge k)\sim A k^{-\mu}`$. Then

```math
 d_N\sim A N^{-\mu}\log N,\qquad
 \left\|\nu_N-\frac{1}{\log N}\int_0^{\log N}P_a^0(J,\cdot)\,da\right\|_{\rm TV}\longrightarrow0.
 \qquad\text{(1)}
```

Thus an independent uniform age on $`[0,\log N]`$ specifies the mixture for the **entire rooted component**, even for tests that vary with $`N`$. For its ordinary vertex count $`V_N`$,

```math
 \frac{\log V_N}{\log N}\ \Longrightarrow\ {\rm Uniform}[0,1].
 \qquad\text{(2)}
```

The degree of the exterior parent and the top vertex's birth mark are omitted from the observations. Under uniform sampling of all ordinary vertices the limiting observation at fixed $`N`$ is an isolated vertex. Indeed at most $`(N+1)C_n`$ vertices can lie in the components counted here. Sampling a surviving core uniformly would instead weight a component by its number of surviving cores.

## The component measure at a boundary edge

We use three results for the original model. Claim 10 gives almost-sure convergence of the empirical neighbourhood of a uniformly sampled core to a deterministic law $`\pi`$ on locally finite decorated core trees with one end. The ancestral direction points toward that end. Its proof gives the convergence with finite ancestral information, and the downward fringe at a typical core is finite almost surely. The birth-time identity represents $`\pi`$ by the physical evolution of a newborn in its live environment, sampled from the law $`\beta`$ at a core birth and observed at an independent $`{\rm Exp}(\mu)`$ age. From $`\beta`$, the newborn’s full degree $`D(t)`$ satisfies $`e^{-t}D(t)\to W`$ almost surely with $`A=\mathbb E_\beta W^\mu>0`$. We use $`0<\mu<1`$. [Oriented empirical convergence](../../../../../astra/source/round2/fringe/spatial-limit.md) · [Finite downward fringes](../../../../../astra/source/round2/fringe/other-marginals.md) · [Birth law and degree-tail constant](../../../../../astra/source/round2/degree-law.md)

Every component counted here except a possible founder component has one top $`v`$, of original degree at most $`N`$, whose parent has degree greater than $`N`$. Its component $`F_N(v)`$ is the part of its downward fringe connected to $`v`$ after deletion. For $`|f|\le1`$ define the finite measure

```math
 M_N(f)=\mathbb E_{\pi}\left[
  {\bf1}_{D_o>N}\sum_{v\ {\rm core\ child\ of}\ o}
  {\bf1}_{D_v\le N} f(F_N(v))\right].
 \qquad\text{(3)}
```

Counting first at the top and then at its parent gives (3) by the mass-transport identity for the uniform-core law. The child count is integrable under $`\pi`$. Equivalently, one can count at the top alone, where the indicator is bounded by one. For each finite component type the event can be read in a finite neighbourhood by examining its boundary degrees. The founder exception costs at most $`1/C_n`$, which tends to zero almost surely. The total mass converges by counting top edges. Since downward fringes under $`\pi`$ are finite, the pointwise type limits account for all this mass. Truncation to a finite list of types therefore gives convergence in total variation. In particular $`d_N=M_N(1)`$ and $`\nu_N=M_N/d_N`$ whenever $`d_N>0`$. Positivity for all sufficiently large $`N`$ follows from (1).

## Removing the conditioning on the final parent degree

Start the prospective high parent at its own birth in a $`\beta`$-distributed environment. It has one real leaf, one core parent, and no core children. Let $`K_t`$ count its promotions up to $`t`$. Its core degree is $`1+K_t`$, and the predictable promotion rate is $`1-(1+K_{t-})/D(t-)`$. Its clock dominates $`K_t`$ by a rate-one Poisson count. Conditional on its environment the gain rate of $`D(t)`$ is at most $`D(t)`$. A rate-one-per-member Yule process from two thus gives $`\mathbb E_\beta D(t)^2\le 6e^{2t}`$.

Put $`L=\log N`$ and $`t=L+y`$. The birth-time identity applied to (3), first to nonnegative bounded truncations of the sum and then by monotone convergence, gives

```math
 \frac{N^\mu}{L}M_N(f)=
 \int_{-L}^\infty \mu e^{-\mu y}\,
 \mathbb E_\beta\left[
  {\bf1}_{D(t)>N}\frac1L
  \sum_{\tau\le t}{\bf1}_{D_{v_\tau}(t)\le N} f(F_N(v_\tau;t))
 \right]dy .
 \qquad\text{(4)}
```

Here $`\tau`$ ranges over the focal parent’s promotions and $`v_\tau`$ is its new core child. At such a promotion the child's own branch starts with two ordinary vertices. Formula (4) initially conditions on the parent's **final** degree. This can bias the law of the child's earlier arrows.

Fix $`y`$ with $`\Pr_\beta(W=e^{-y})=0`$, which excludes at most countably many $`y`$. Take $`b=\sqrt L`$, so $`0<b<t`$ for large $`N`$, and form the event

```math
 E_N(y)=\{e^{-b}D(b)>e^{-y}\}.
```

Both $`{\bf1}_{E_N(y)}`$ and $`{\bf1}_{D(t)>N}`$ converge almost surely to $`{\bf1}_{W>e^{-y}}`$. By Cauchy–Schwarz the error inside the integrand after replacing the latter indicator by $`{\bf1}_{E_N(y)}`$ is at most
$`\sqrt{t+t^2}\Pr_\beta(E_N(y)\mathbin{\triangle}\{D(t)>N\})^{1/2}/L`$.
It tends to zero uniformly over the test functions. Discarding promotions through $`b`$ costs at most $`b/L`$. The conditioning for the remaining births is now measurable at $`b`$.

## The branch after a promotion

On $`E_N(y)`$, at every retained promotion $`\tau>b`$ the parent’s degree is already greater than $`e^{b-y}`$ and stays above that threshold. Conditional on the history through that promotion, couple the child branch to $`J`$ using the same subsequent internal arrows. The first arrow from the parent into the branch can separate their shapes. Its predictable rate before separation is at most $`e^{y-b}`$, so the probability of separation before $`t`$ is at most $`t e^{y-b}`$. Upward arrows from the branch can still change the parent. Their internal effect agrees with a null slot selection, and the bound on the reverse rate holds after those changes.

Independently of this coupling, a physical branch of current ordinary size $`m`$ has total rate of increase
$`m-1/D_{\rm top}+1/D_{\rm parent}\le m+1/2`$.
With $`u=t-\tau`$, its conditional mean size is at most $`(5/2)e^u`$. If its size is at most $`N`$, every degree in the branch, including the top's exterior slot, is at most $`N`$. Hence the probability of any further degree deletion in the branch is at most $`(5/2)e^{y-b}`$ for these births. Summing the two errors over promotions and dividing by $`L`$ costs at most
$`2(t/L)(t+5/2)e^{y-b}`$ on the unit ball of test functions. This tends to zero at every fixed $`y`$.

Write $`h_f(u)=P_u^0 f(J)`$. Having removed the final conditioning, we can take conditional expectations at the individual births and then compensate the promotion clock. The resulting inner expectation is

```math
 \frac1L\mathbb E_\beta\left[{\bf1}_{E_N(y)}
  \int_b^t \left(1-\frac{1+K_{r-}}{D(r-)}\right)h_f(t-r)\,dr\right].
```

Replacing the promotion rate by one costs at most $`e^{y-b}(t+t^2/2)/L`$, since $`\mathbb E_\beta K_r\le r`$. Changing the resulting integration interval for the age from $`[0,t-b]`$ to $`[0,L]`$ costs at most $`(b+|y|)/L`$. Thus at each of the specified $`y`$, the inner expectation in (4) differs by $`o(1)`$, uniformly for $`|f|\le1`$, from

```math
 \Pr_\beta(W>e^{-y})\,\frac1L\int_0^L h_f(u)\,du .
 \qquad\text{(5)}
```

For $`y\le0`$ and $`t\ge0`$, the absolute physical integrand before the density factor is at most
$`\sqrt{t+t^2}\Pr(D(t)>N)^{1/2}/L\le\sqrt{12}\,e^y`$ when $`L\ge1`$.
For $`y>0`$ use $`\mathbb E K_t/L\le 1+y`$. The bounds for (5) follow from $`\mathbb E W^2\le6`$. Multiplication by $`e^{-\mu y}`$ gives integrable bounds at both ends because $`0<\mu<1`$. Extend the physical integrand by zero for $`y<-L`$. Dominated convergence for the uniform error, and the identity

```math
 \int_{\mathbb R}\mu e^{-\mu y}\Pr_\beta(W>e^{-y})\,dy
 =\mathbb E_\beta W^\mu ,
```

give (1), first as convergence of the unnormalized measures after multiplication by $`N^\mu/L`$, then after division by their masses.

## The volume

The ordinary size $`m`$ of the inert branch grows at rate $`m-1/D_{\rm top}`$, between $`m-1`$ and $`m`$. Its size is bounded on a common probability space below by one plus a Yule process from one, and above by a Yule process from two. Both comparison processes have finite positive limits after multiplication by $`e^{-u}`$, so the log size of the branch divided by its elapsed time $`u`$ tends to one almost surely. Let $`V`$ be uniform on $`[0,1]`$, independent of the branch. At $`u=LV`$ the log size divided by $`L`$ tends to $`V`$ almost surely. Total variation controls even an $`N`$-dependent observation, giving (2). These results retain the sequential order of limits. A threshold growing with the observed finite tree size requires an additional uniform estimate.

The independent derivation and directed-arrow replay check the boundary rate, the time at which the parent event becomes measurable, and the normalization by component count. [Independent check and replay](../../../../../../reviews/cold/astra/open/forest-age/independent.md)
