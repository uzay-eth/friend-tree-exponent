> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../README.md#claim-26) for exact scope and citations.

# The probability that a birth label becomes a hub

Start the undirected random friend tree with the edge joining labels 1 and 2. At each birth choose a vertex uniformly, choose one of its neighbours uniformly, and attach the new leaf to that neighbour. Write $`Z_i=\lim_{n\to\infty}D_i(n)/n`$, $`H_i=\{Z_i>0\}`$, and $`\mathcal F_n`$ for the history through size $`n`$. The answer to the hub-probability part of Question 3 in Section 7 of *Random Friend Trees* is
```math
 \Pr(H_1)=\Pr(H_2)>\Pr(H_3)>\Pr(H_4)>\cdots .
```
Every comparison between labels at least 3 already holds conditional on their common past before either is born. [Original question](https://arxiv.org/html/2403.20185v1#S7) · [Two cold reviews and the frozen proof](../global-transport/hostile/RELEASE.md)

## A gap conditional on the past

Put $`c=371/12500`$. For every $`n\ge2`$,
```math
 \Pr(H_{n+1}\mid\mathcal F_n)-\Pr(H_{n+2}\mid\mathcal F_n)
 \ \ge\ \frac{c}{n(n+1)}. \qquad\text{(1)}
```
Here the future includes both placements. For an old vertex $`u`$ write $`r_u=\sum_{z\sim u}1/D_z`$, $`w_u=r_u/(D_u+1)`$, and $`s_u=w_u/D_u`$. Let $`I(u)`$ be the difference in hub probability between the middle and tip of a newly grown two-vertex path at $`u`$. Let $`H(u,v)`$ be the difference between two new leaves, one at each of the old vertices $`u,v`$, measured from that common resulting tree. Pair the two orders of distinct old targets. A nonadjacent pair cancels because the first birth leaves the other target's rate unchanged. At repeated targets the two leaves are exchangeable. Multiplying the conditional gap by $`n(n+1)`$ gives exactly
```math
 \sum_u w_u I(u)+\sum_{\{u,v\}\text{ old edge}}(s_v-s_u)H(u,v).
```
Reversing the orientation of an old edge leaves its summand unchanged. Replacing the two leaves' future probabilities by the same value before multiplying would discard an actual signed term. [Two-birth identity](../adjacent/README.md)

The frozen proof bounds a path middle using its tracked leaves. When the last tracked leaf is promoted, the resulting core can send a fresh leaf back. Keeping that return leaves enough surplus to pay the adverse old edges. In the tree induced by old cores, contract each maximal path of internal degree-two vertices of total degree two. Each remaining edge spends one allowance $`b=0.017`$. At each retained endpoint, the proof pays the additional directed-edge or contracted-chain cost for each real incident port, using its leaf stock and the other endpoint's reciprocal degree. The one- and two-site chains require separate estimates. Starting at length three, each further internal site adds positive credit and no further debit. With $`K\ge2`$ retained cores, the tree has $`K-1`$ edges and $`\sum_u(\deg_{\rm core}(u)-2)=-2`$. The verified local inequalities with terminal allowance $`a=0.00634`$ therefore leave $`b+2a=c`$. This includes the two-terminal tree left by contracting an entire core path. [Complete edge and chain budgets](../global-transport/hostile/edge-dual/README.md) · [Simultaneous endpoint checks](../global-transport/hostile/interfaces/low-end.md)

The same numerical constant in (1) covers the seed and every star. The unrounded local first-send barrier with $`j`$ other leaves is at most $`1/(j+1)`$. Indeed, at $`j\ge1`$, substituting that expression into the comparison recurrence gives numerator minus denominator times the present value equal to $`1-D/(j+2)\le0`$ for $`D\ge j+2`$. At $`j=0`$ the value one bounds a probability. On a star of centre degree $`d\ge4`$, $`s_{\rm centre}=1/(d+1)>s_{\rm leaf}=1/(2d)`$. Its adverse direction compares the first-send cap at the centre after the two births, at most $`1/d`$, with the opposite leaf's rescue bound $`L(2)>1/4`$. The charge for every old edge is therefore zero. The leaves' serial contributions alone sum to $`G(2)/2\ge189409485129/(2\cdot10^{12})>c`$. On the edge seed and at star degrees 2 and 3, the exact signed scores, before division by $`n(n+1)`$, exceed $`0.18940948,\ 0.20181908,\ 0.22577249`$, respectively. The fixture below restores the total serial weight in those scores. [Scalar barriers and star comparison](../global-transport/hostile/discharge/README.md) · [Rounding convention](../global-transport/hostile/RELEASE.md#cold-reviews-and-rounding)

Taking expectations in (1) gives the consecutive comparisons starting with labels 3 and 4. The two founders have the same law by edge symmetry. At size three the third vertex is a leaf, and founder 2 has probability one half to be the centre. The centre has hub probability at least $`3/4`$, and the leaf has hub probability at most $`1-108/305`$. Hence $`\Pr(H_2)-\Pr(H_3)\ge127/2440`$. Fable's earlier exact programs certified a finite prefix of the ordering and its transport proof compared the full degree-share tails on early pairs. [Founder comparison](../early/README.md) · [Fable's finite hub-probability certificates](../../../../../../../reviews/cold/astra/attacks/fable3/extension/README.md) · [Fable's early transport result](../../../attacks/fable5/README.md)

## Later labels viewed from the same earlier tree

For $`2\le m\le i-1`$ and $`i<j`$, conditioning (1) on the earlier history $`\mathcal F_m`$ and summing gives
```math
 \Pr(H_i\mid\mathcal F_m)-\Pr(H_j\mid\mathcal F_m)
 \ \ge\ c\left(\frac1{i-1}-\frac1{j-1}\right). \qquad\text{(2)}
```
The positive lower bound is uniform over every such history. Letting $`j`$ increase and using nonnegativity yields $`\Pr(H_i\mid\mathcal F_m)\ge c/(i-1)`$. Thus the expected number of future hubs with labels $`i,\ldots,j`$, conditional on this same past, is at least $`c\sum_{v=i}^j1/(v-1)`$. The tower property proves these statements on the original common tree, including when the two future hub events are dependent. They also hold with $`m`$ replaced by an almost surely finite stopping time and deterministic positive offsets for the two future labels, by conditioning on the observed stopping size.

For the unconditional degree-share tails, strict hub gaps and right continuity yield an interval starting at zero for each finite pair. Claim 24 covers all thresholds beyond a single finite label. Taking the minimum over the remaining finitely many pairs gives one positive near-zero interval for *every* adjacent pair starting with label 2. The comparisons between arbitrary labels at every intermediate positive threshold remain open. [Late-label theorem](../eventual-domination/README.md)

## Replay

Run `python3 -IB round2/open/age-order/hub-ranking/check.py` from the repository root. It checks the literal seed and star scores, the scalar barrier, the early founder arithmetic and the reciprocal sum in (2). The frozen proof lists the 12 producer and review commands for the general core-tree bound. [Release instructions](../global-transport/hostile/RELEASE.md#archive-replay)
