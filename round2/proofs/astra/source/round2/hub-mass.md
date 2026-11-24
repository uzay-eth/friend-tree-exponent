# The total normalized degree of the hubs

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../SOURCES.tsv)).
> **Review and dependencies:** Claim 7: Claude and Gemini; GPT standalone audit; inputs include sublinear core growth and fixed-label degree limits. See the [Astra ledger](../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Statement submitted for independent review

Let $`Z_v=\lim_n D_v(n)/n`$, whose almost-sure existence for every fixed vertex is given by Addario-Berry et al. Then
```math
 \sum_{v\ge1} Z_v=1\qquad\text{almost surely}.
\qquad\text{(1)}
```
The conditional distribution of the recipient label $`W_n`$ in one new attachment tends in total variation to the random probability vector $`(Z_v)_{v\ge1}`$. Its normalized degree $`D_{W_n}(n)/n`$ therefore has conditional law converging weakly, almost surely, to $`\sum_v Z_v\delta_{Z_v}`$. These are the mass and empirical-degree conjectures stated in Section 3 of [Addario-Berry et al.](https://doi.org/10.1002/rsa.70082).

## A bound in expectation

This is a consequence of an existing upper bound on core growth. It does not require a new value or characterization of $`\mu`$. Write $`a_n=\mathbb E C_n`$, and let $`\lambda<1`$, $`K<\infty`$ satisfy $`a_n\le K n^\lambda`$ for $`n\ge3`$, as in the [published repository bound](../proofs/bracket-motifs.md).

Given a tree at size $`j\ge3`$, a particular core of degree $`D`$ gains a leaf with probability at most $`D/j`$, since each of its neighbours contributes at most one to the numerator of the attachment probability. Its normalized degree $`D(j)/j`$ is therefore a nonnegative supermartingale from its promotion onward. A core promoted at size $`r`$ begins with degree two. For $`n\ge m\ge3`$, every core of label larger than $`m`$ present at size $`n`$ was promoted at some size $`r>m`$. At most one promotion occurs at each step. Conditioning on the tree just after each promotion bounds the expected degree mass at $`n`$ from late-labelled cores by
```math
 \mathbb E\sum_{\substack{v>m\\D_v(n)\ge2}}\frac{D_v(n)}n
 \le 2\sum_{r=m+1}^n\frac{a_r-a_{r-1}}r
 \le 2K \sum_{j=m}^\infty\frac{j^{\lambda-1}}{j+1}
 \le 2K\left(m^{\lambda-2}+\frac{m^{\lambda-1}}{1-\lambda}\right).
\qquad\text{(2)}
```
The second inequality uses the one-step bound $`a_{j+1}-a_j\le a_j/j`$ on the promotion probability.

At size $`n`$ the degree sum over cores is $`2(n-1)-(n-C_n)=n+C_n-2`$. Its expected value divided by $`n`$ tends to one. For each of the first $`m`$ labels, $`D_v(n)\mathbf1_{D_v(n)\ge2}/n\to Z_v`$ almost surely and is bounded by one. Dominated convergence and (2) give $`\mathbb E\sum_{v\le m}Z_v\ge1-O(m^{\lambda-1})`$. The reverse inequality $`\sum_v Z_v\le1`$ holds almost surely by the fixed-vertex leaf-neighbour limits in Proposition 3.1 of Addario-Berry et al., or by Fatou and $`C_n/n\to0`$ almost surely. Sending $`m\to\infty`$ proves (1).

For fixed $`v`$, $`\mathbb P(W_n=v\mid T_n)=n^{-1}(L_v(n)+\sum_{w\sim v,\ D_w(n)\ge2}1/D_w(n))`$. The second term after division by $`n`$ is bounded by $`C_n/n\to0`$ almost surely, and $`L_v(n)/n\to Z_v`$. Thus the conditional weights converge pointwise on the countable label space. Their sums and the sum of the limits equal one, implying total variation convergence (truncate the first $`m`$ labels and then send $`m\to\infty`$). For a bounded continuous function of the normalized degree, the first $`m`$ summands converge individually and the remaining weight tends to zero with $`m`$. This gives the weak convergence stated above.

## Using only an almost-sure power upper bound

The qualitative conclusion also follows from $`C_n=O(n^\lambda)`$ almost surely with deterministic $`\lambda<1`$, as supplied by the earlier low-degree theorem in the [preprint](https://arxiv.org/abs/2403.20185). Here is the localization, since an almost-sure power upper bound need not imply the expectation bound used in (2). Let $`\tau=\inf\{n\ge3:C_n>K n^\lambda\}`$. Start at size $`m`$ with no counted cores. On $`\{\tau>m\}`$, add each core promoted after $`m`$ up to and including the stopping step, and freeze the sum of their normalized degrees at $`\tau`$; set the sum to zero if $`\tau\le m`$. While running, old entries have nonpositive conditional drift and a promotion adds $`2/r`$ at size $`r`$. Hence the expected stopped sum is bounded by
```math
 2\,\mathbb E\sum_{r=m+1}^n\frac{C_{r\wedge\tau}-C_{(r-1)\wedge\tau}}{r}.
```
For $`r\ge3`$, $`C_{r\wedge\tau}\le K r^\lambda+1`$ by the unit jump bound. Partial summation makes the display $`O_\lambda((K+1)m^{\lambda-1})`$ uniformly in $`n`$. On $`\{\tau=\infty\}`$ the sum dominates the actual normalized degree mass of all cores labelled after $`m`$. There the total normalized core degree tends to one. Fatou gives
```math
 \mathbb E\left[\mathbf1_{\{\tau=\infty\}}
     \left(1-\sum_{v\le m}Z_v\right)\right]\le O_\lambda((K+1)m^{\lambda-1}).
```
Let $`m\to\infty`$, then $`K\to\infty`$. The events $`\{\tau=\infty\}`$ exhaust an event of probability one. This proves (1) even under the weaker upper-bound hypothesis.
