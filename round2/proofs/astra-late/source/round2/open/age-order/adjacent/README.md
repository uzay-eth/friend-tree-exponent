> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# Adjacent hub probabilities at late birth labels

Write $`D_i(n)`$ for the degree at size $`n`$ of vertex $`i`$ in the undirected random friend tree, $`Z_i=\lim_{n\to\infty}D_i(n)/n`$, and $`p_i=\Pr(Z_i>0)`$. Start the process with the vertices labelled 1 and 2 joined by an edge. With $`0<\mu<1`$ its core-growth exponent and $`c>0`$ the late-birth hub constant,
```math
 \lim_{n\to\infty}n^{2-\mu}(p_{n+1}-p_{n+2})=(1-\mu)c.
 \qquad\text{(1)}
```
In particular **the probability of being a hub is strictly decreasing between every two consecutive sufficiently late labels**. The argument leaves the first possible exceptional label and the ordering across all early labels undetermined. The two endpoints at labels 1 and 2 are exchangeable. This addresses the probability-of-hubhood part of Question 3 in Section 7. [Original question](https://arxiv.org/html/2403.20185v1#S7) · [First-order late-birth law](../late/README.md)

## Exact cancellation in two births

In a tree $`T`$ of size $`n`$, let $`D_u`$ be the degree of $`u`$, set $`r_u=\sum_{v\sim u}1/D_v`$, and let $`\delta_u=1/[D_u(D_u+1)]`$. A birth attaches a leaf to $`u`$ with probability $`r_u/n`$. For any reachable tree with a marked vertex $`z`$, let $`h(S,z)`$ be the conditional probability its degree is ultimately unbounded. For a finite starting tree this equals its probability of a positive limiting degree fraction. For distinct old vertices $`u,v`$, add a leaf at each, and write $`H(u,v)`$ for $`h`$ of the leaf at $`u`$ minus $`h`$ of the leaf at $`v`$ in that single final tree. Thus $`H(v,u)=-H(u,v)`$. Also let $`I(u)`$ be $`h`$ of the middle vertex minus $`h`$ of the tip after growing a path of two vertices from $`u`$. All terms use the physical future from the indicated final tree. [Tagged degree test](../late/README.md)

Given a first attachment at $`u`$, a distinct old vertex $`v`$ has second attachment rate $`r_v-\mathbf1_{\{u\sim v\}}\delta_u`$, before division by $`n+1`$. The new leaf sends only to $`u`$, and the degree of $`u`$ has increased by one. Reversing the order at distinct adjacent targets and combining their contributions leaves $`(r_v\delta_v-r_u\delta_u)H(u,v)`$ per unordered edge. The two orders at nonadjacent targets cancel. If both births use the same old target, the final tips are exchangeable. If the second birth attaches to the first newborn, the corresponding weight is $`r_u/(D_u+1)`$. These cases exhaust the two-birth law.

Charge each old vertex with the oriented sum
```math
 S(T,u)=\frac{r_u}{D_u+1}\,I(u)
       -r_u\delta_u\sum_{v\sim u} H(u,v).
```
For each old core $`o`$, of degree at least two, put $`G(T,o)=S(T,o)+\sum_{\ell\sim o,\,D_\ell=1}S(T,\ell)`$. At $`n\ge3`$ each old leaf belongs to exactly one of these sums. The tower property at the time of the second birth gives
```math
 n(n+1)(p_{n+1}-p_{n+2})
 =\mathbb E\sum_{o\text{ core in }T_n}G(T_n,o). \qquad\text{(2)}
```
The common final environment permits the signed cancellation, including when one of the old targets was a leaf. Each $`\lvert S(T,u)\rvert\le 2r_u/(D_u+1)`$. The core term is bounded by 2 since $`r_o\le D_o`$. An old leaf adjacent to $`o`$ contributes at most $`1/D_o`$, so $`\lvert G\rvert\le3`$ on all trees.

## Local continuity, then the sign

The future hub probability $`h`$ extends to locally finite marked trees. We recall why it is continuous at the marked trees above. For a tag starting at degree one or two, a lower test is the probability of hitting degree $`d\ge3`$ before time $`t`$, times $`1-1/d`$. At a new hit, a fresh leaf starts a birth–death family with per-leaf upward rate one and death rate $`1/d`$ that can be embedded below the actual stock. Survival forces unbounded degree. On sending $`t`$, then $`d`$, to infinity, the lower tests approach $`h`$. For an upper test at time $`t`$, bound the probability of any further gain by the sum, capped at one, of the probabilities of a first subsequent arrow from the current neighbours. The stock bound for a neighbour with $`k`$ other leaves tends to zero as $`k\to\infty`$. On the event of bounded degree the tag eventually has a fixed finite set of neighbours. Each such neighbour has unbounded degree and diverging leaf stock: two endpoints of eventually bounded degree would still have an eventual positive directed rate, and unbounded degree implies diverging leaf stock by the tagged-family estimate. The excess of the upper test tends to zero by bounded convergence. Both finite-history bounds can be localized to finitely many initial cells for fixed $`t`$. [Finite-history and stock estimates](../late/README.md) · [Local construction](../../../../../../astra/source/round2/operator.md) · [Independent continuity check](../../../../../../../reviews/cold/astra/checks/adjacent-age-check.md)

The finitely many two-birth insertions made at any one rooted state preserve local convergence. Thus $`G`$ is bounded and continuous. Let $`C_n`$ count cores and $`\pi_n`$ be their uniform empirical law. The local-limit theorem gives $`\pi_nG\to\pi G`$ almost surely. The core scaling theorem gives $`C_n/n^\mu\to X`$ with uniform integrability. Applying these inputs to (2) proves existence of a finite constant $`a=\mathbb E[X]\pi G`$ such that
```math
 p_{n+1}-p_{n+2}=(a+o(1))n^{\mu-2}.
```
The series of differences is absolutely summable because $`\mu<1`$. Since $`p_n\to0`$, summing from $`n`$ to infinity gives $`n^{1-\mu}p_{n+1}\to a/(1-\mu)`$. The error after summation is bounded by the supremum of the remaining $`o(1)`$ errors times the power tail. The separately proved $`p_n\sim c n^{\mu-1}`$, with $`c>0`$, identifies $`a=(1-\mu)c`$ and proves (1). [Core local limit](../../../../../../astra/source/round2/fringe/spatial-limit.md) · [Core scaling and integrability](../../../../../../astra/source/round2/standalone/claim-6.md) · [Late-birth coefficient](../late/README.md)

Two independently implemented GPT-family checks enumerate the original two births on labelled trees. One compares 456 exact harmonic rows, including both signs of the old-edge term. Dropping the path of length two already fails on the initial three-vertex tree. A second checks 128 state and endpoint comparisons, with a nonzero path contribution in all 128. Its review also derives comparisons at thresholds of order $`1/n`$. Both checks use finite-horizon endpoint tests for the identity. Formula (1) follows from the limiting argument above. [Harmonic rows and replay](../../../../../../../reviews/cold/astra/checks/adjacent-age-check.md) · [Second enumerator and threshold comparison](../../../../../../../reviews/cold/astra/checks/adjacent-age-sol.md)
