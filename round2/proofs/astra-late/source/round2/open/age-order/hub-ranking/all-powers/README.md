> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../README.md#claim-27) for exact scope and citations.

# Moment order for every birth label

Start the undirected random friend tree from a single edge. Given the tree at each birth, select a vertex uniformly, select one of its neighbours uniformly, and attach the new vertex to that neighbour. Let $`Z_i`$ be the limiting fraction of total tree size contributed by the degree of label $`i`$, and let $`\mathcal F_n`$ be the tree history at size $`n`$. For every real $`p\ge1`$,
```math
 \mathbb E Z_1^p=\mathbb E Z_2^p>\mathbb E Z_3^p>\mathbb E Z_4^p>\cdots .
```
For two future labels the strict inequality holds on every common history before either birth. It requires no starting-label cutoff. Earlier results include the strict early tail comparisons of Fable Claim 5 and the eventual all-threshold theorem of Astra Claim 24. [Finite transport audit](../../../../attacks/fable5/README.md) · [Late-label theorem](../../eventual-domination/README.md)

## Conditional gap

For $`n\ge2,\ p\ge1`$, set $`c=9/5000`$. Then
```math
 \mathbb E[Z_{n+1}^p-Z_{n+2}^p\mid\mathcal F_n]
 \ \ge\ c\frac{\Gamma(p+1)\Gamma(n)}{\Gamma(n+p+2)}. \qquad\text{(1)}
```
The constant can be increased to $`1/40`$ throughout $`1\le p\le3`$ and to $`21/1000`$ for $`p\ge7`$.

After both placements the size is $`m=n+2`$. In the rate-one continuous-time construction, the true size amplitude $`S`$ has law $`\mathop{\mathrm{Gamma}}\nolimits (m,1)`$ and is independent of the subsequent embedded jump chain conditional on the completed tree. The physical degree amplitude $`W_i=S Z_i`$ therefore satisfies $`\mathbb E W_i^p=(m)_p\mathbb E Z_i^p`$ at that tree. Pairing the two orders of placement gives a signed numerator $`N_p`$ for the continuation test $`\mathbb E W_i^p/\Gamma(p+1)`$. The factor converting $`N_p`$ to the conditional gap is $`\Gamma(p+1)/[n(n+1)(n+2)_p]`$, equal to the gamma ratio in (1). [Two-birth identity](../../adjacent/README.md) · [Physical clock and moment proof](../mellin/README.md)

| Real powers | Certified lower bound for $`N_p`$ | Computation or estimate |
|---|---:|---|
| $`1\le p\le3`$ | $`1/40`$ | On 87 closed intervals each vertex pays its incident debits. On the other 169 the proof sums the debits along paths of degree-two cores before charging their ends. |
| $`3\le p\le7`$ | $`1892360505093/10^{15}`$ | 2,048 closed intervals compare the serial middle with both classes of adverse tips. |
| $`7\le p<\infty`$ | $`22890923/1058158080`$ | A continuous-time lower comparison for the serial middle gives the same alternating-sum bound for all these powers. |

In the last two rows the original lower bounds include a factor $`A=\sum_u r_u/(D_u+1)`$, where $`r_u=\sum_{z\sim u}1/D_z`$ and $`D_u`$ is the old degree. Every initiator sends at total rate one, so $`\sum_u r_u=n`$. Since $`D_u+1\le n`$, we have $`A\ge1`$, also on the seed edge. The three ranges give (1) on every finite tree. [Real powers from one to three](../mean-ledger/real/README.md) · [Real powers from three to infinity](../mellin/README.md)

## The founders

At size three the tree is a two-leaf star. The centre has two tracked leaves. The stock submartingale of the serial-middle proof, started at that time, gives $`\mathbb E W_{\rm centre}^p/\Gamma(p+1)\ge(2-1)_p/\Gamma(p+1)=1`$ for $`p\ge1`$. A leaf at this same tree retains its parent, whose degree is always at least two. Its attachment rate at future degree $`D`$ is at most $`D-1+1/2`$. The affine pure-birth comparison started at degree one gives $`\mathbb E W_{\rm leaf}^p/\Gamma(p+1)\le(1/2)_p/\Gamma(p+1)`$. This last ratio is the $`p`$-th moment of a $`\mathop{\mathrm{Beta}}\nolimits (1/2,1/2)`$ variable. It is at most its mean $`1/2`$.

Label 2 is the centre with probability $`1/2`$. Otherwise it is a leaf exchangeable with label 3. Dividing the normalized half-difference by the size-three conversion factor gives
```math
 \mathbb E[Z_2^p-Z_3^p]\ge\frac{\Gamma(p+1)}{4(3)_p}
   =\frac{1}{2(p+1)(p+2)}.
```
Exchanging the labels of the initial edge preserves the entire tree law and gives the tie between the founders. [Stock submartingale](../mellin/README.md) · [Physical comparison at a leaf](../tip-moments/README.md)

## A common earlier history and mixtures of powers

For $`2\le m\le i-1`$ and $`i<j`$, condition (1) on $`\mathcal F_m`$ for each intervening pair and sum. Euler's beta function satisfies $`\mathrm B(n,p+1)-\mathrm B(n+1,p+1)=(p+1)\Gamma(p+1)\Gamma(n)/\Gamma(n+p+2)`$. Hence
```math
 \mathbb E[Z_i^p-Z_j^p\mid\mathcal F_m]\ge
 \frac{c}{p+1}\{\mathrm B(i-1,p+1)-\mathrm B(j-1,p+1)\}. \qquad\text{(2)}
```
The same bound applies after an almost surely finite stopping size for two deterministic positive offsets, by partitioning on the observed size.

Let $`\rho`$ be any nonzero finite positive measure on $`[1,\infty)`$. The same label order holds for the expectations of $`\int z^p\,\rho(dp)`$, by conditional integration of (1), (2), or the founder bound. Taking atoms of mass $`t^k/k!`$ at $`k=1,2,\ldots`$ gives strict order for $`\mathbb E e^{tZ_i}`$ at every $`t>0`$. General first-order stochastic dominance at intermediate positive thresholds remains open. [Thresholds covered at every label](../../global-transport/annealed/extension/first-sender/joint/README.md) · [Exact zero-threshold ranking](../README.md)

Run the exact arithmetic check in this directory and the two interval checkers from the repository root.

~~~sh
python3 -IB round2/open/age-order/hub-ranking/all-powers/check.py
python3 -IB round2/open/age-order/hub-ranking/mean-ledger/real/check.py
python3 -IB round2/open/age-order/hub-ranking/mellin/check.py
~~~
