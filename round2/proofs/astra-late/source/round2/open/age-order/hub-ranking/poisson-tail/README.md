> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# A Poisson mixture of beta laws bounds the future share of a vertex with live leaves

Fix any size-$`m`$ tree in the original random friend process, a tagged vertex and $`j`$ of its current leaf neighbours. Write $`Z`$ for its limiting degree divided by tree size. If $`1\le j\le m-1`$, then for every $`0\le x<1`$
```math
 \Pr(Z>x\mid T_m)\ \ge\;
 H_{m,j}(x)=e^{-1}\sum_{k=0}^{j-1}\frac{1}{k!}
                 \Pr\{\mathrm{Beta}(j-k,m-j+k)>x\}.                         \qquad\text{(1)}
```
The lower law is the following probability mixture. Draw $`K\sim\mathrm{Poisson}(1)`$. For $`K\ge j`$ use zero, and for $`K<j`$ use the indicated beta variable. This is a stochastic comparison at every threshold, conditional on an arbitrary starting tree. Its $`j=1`$ case gives $`e^{-1}(1-x)^{m-1}`$ for the serial middle immediately after two births, even with an arbitrarily high-degree old parent.

A consequence is a **single-digit conditional birth-scale cutoff**. Write $`b=n+1`$, $`q=1-x`$. For every history at old size $`n\ge2`$, every $`0<x<1`$ with $`bx/q\ge6`$ satisfies
```math
 \Pr(Z_b>x\mid\mathcal F_n)-\Pr(Z_{b+1}>x\mid\mathcal F_n)
 \ \ge\ \frac{19}{10000}\frac{q^b}{n(n+1)}.                                \qquad\text{(2)}
```
In particular $`bx\ge6`$ suffices. From label $`b=20`$ onwards the condition holds throughout $`[3/13,1)`$. Fable’s two-program certificates for $`2\le b\le19`$ give unconditional all-label order on that entire interval. Including its single receipt at $`20\to21`$ gives $`[2/9,1)`$. The founders have identical laws. [Earlier physical cutoff](../../global-transport/annealed/extension/README.md) · [Scope of Fable’s early receipts](../../../../attacks/fable5/README.md)

## Discrete harmonic identity

For $`0<r<m`$ put $`B_{m,r}(x)=\Pr\{\mathrm{Beta}(r,m-r)>x\}`$. Set $`B_{m,r}=0`$ for integers $`r\le0`$, and $`B_{m,m}=1`$ on $`0\le x<1`$. For integer $`r\ge1`$ this is the probability that a $`\mathrm{Bin}(m-1,x)`$ count is at most $`r-1`$. It increases with $`r`$. The usual urn identity is
```math
 m B_{m,r}=r B_{m+1,r+1}+(m-r)B_{m+1,r}
 \qquad(0\le r\le m).
```
Extend $`H`$ by zero at stock zero. For $`1\le j<m`$,
```math
 m H_{m,j} =
 j H_{m+1,j+1}+(m-j-1)H_{m+1,j}+H_{m+1,j-1}.                              \qquad\text{(3)}
```
Indeed let $`\pi_k=e^{-1}/k!`$, $`r=j-k`$, and $`\Delta B(r)=B_{m+1,r}-B_{m+1,r-1}`$. Subtract the urn identity term by term from the right-hand side. The discrepancy is $`\sum_{k\ge0}\pi_k[k\Delta B(r+1)-\Delta B(r)]`$, which vanishes by $`k\pi_k=\pi_{k-1}`$. All differences with both indices nonpositive vanish, including at the finite lower edge of the sum.

Continue the tag’s tracked stock by adding every new leaf attached to it and removing a tracked leaf when promoted. At size $`s`$ with stock $`i`$, the probability of a gain is at least $`i/s`$ because each tracked leaf sends to the tag. The probability of loss is $`i/(sD_{\rm tag})\le1/s`$, from a tag-initiated arrow to a tracked leaf. These two events are disjoint. Since $`H_{s+1,i}`$ is increasing in $`i`$, (3) supplies a bounded submartingale for the actual stock. At zero its next value is nonnegative. This construction works under the true conditional tree law at any reached size.

For a fixed $`u>x`$, on $`Z<u`$ the tracked stock divided by $`s`$ is eventually below $`u`$ by a positive margin. Since $`H_{s,i}(u)\le B_{s,i}(u)`$, binomial concentration forces $`H_{s,i}(u)\to0`$ there. Bounded reverse Fatou gives $`H_{m,j}(u)\le\Pr(Z\ge u\mid T_m)`$. Send $`u\downarrow x`$, using continuity of the finite beta mixture. This gives the strict-tail event in (1), including at $`x=0`$.

The comparison law itself is exact for the auxiliary stock chain that gains with probability $`i/s`$, loses with probability $`1/s`$ at positive stock, and absorbs at zero. Its stock fraction is a bounded supermartingale. At every continuity point of its terminal fraction the binomial expression in $`H`$ converges to the corresponding indicator. For a fraction strictly above the threshold, truncate the Poisson sum at a fixed $`K`$ before taking the limit, then send $`K`$ to infinity. Identity (3) identifies the full law from these continuity points. The physical moment harmonic previously derived for fractional powers has these same Poisson weights. Here the threshold test functions themselves satisfy the identity. [Fractional moment harmonic](../fractional-middle/README.md)

## Upper bounds for adverse tips

The common completion size is $`m=b+1`$. A fresh leaf adjacent to a parent of degree at least two has the locked upper $`B_{m,1/2}`$. When the parent has degree at least three use $`B_{m,1/3}`$. With old-tree intensities $`r_u=\sum_{v\sim u}1/D_v`$ and $`A=\sum_u r_u/(D_u+1)`$, all serial tips together have weight $`A`$. The adverse parallel weights sum to at most $`A`$, and at most $`A/3`$ of them can require the degree-two parent cap. Monotonicity of the beta laws in their first shape therefore bounds the total debit by $`A[(4/3)B_{m,1/2}+(2/3)B_{m,1/3}]`$. Also $`A\ge1`$ since $`\sum r_u=n`$ and $`D_u+1\le n`$. All serial middles together buy at least $`A e^{-1}q^b`$ by (1). [Signed placement identity](../../adjacent/README.md) · [Degree-two coefficient bound](../../global-transport/annealed/extension/README.md)

Here is an explicit bound on the two tails without a size cutoff. For $`r`$ equal to $`1/2`$ or $`1/3`$, write $`s=1-r`$, $`B=m-r`$, and $`v=q/(Bx)`$. Changing variables $`z=x+qu`$ in the beta integral gives
```math
 \frac{B_{m,r}(x)}{q^b}
 \le \frac{v^s}{\Gamma(r)}
   \mathbb E\left[(1+(q/x)U)^{-s}\right],\qquad U\sim\mathrm{Beta}(1,B).
```
The gamma-ratio inequality is $`\Gamma(B+r)/\Gamma(B+1)\le B^{-s}`$, by Jensen applied to a gamma variable of shape $`B`$. Under $`bx/q\ge6`$ we have $`v\le1/6`$, while $`b\ge3`$ gives $`B/(B+1)\ge7/9`$. The degree-two Taylor polynomial bounds $`(1+t)^{-s}`$ from above for $`t\ge0`$. Using $`\mathbb E U=1/(B+1)`$ and $`\mathbb E U^2=2/[(B+1)(B+2)]`$ gives
```math
 \mathbb E[(1+(q/x)U)^{-s}]
 \le 1-sv\,\frac{B}{B+1}\left(1-(1+s)v\frac{B}{B+2}\right)
 \le 1-d_s v,\qquad d_s=\frac79\,s\left(1-\frac{1+s}{6}\right).
```
The upper $`v^s(1-d_s v)`$ increases up to $`v=1/6`$. Exact constants are
```math
 (6\pi)^{-1/2}<231/1000,\quad 6^{-2/3}<3031/10000,\quad
 \Gamma(1/3)>13/5.
```
For the last inequality substitute $`t=u^3`$ in the gamma integral and integrate the odd Taylor polynomial of order seven for $`e^{-u^3}`$ up to $`u=3/2`$. The resulting rational exceeds $`13/5`$. The two normalized caps with their weights now sum to at most $`8322209/22744800`$. Since $`e^{-1}>36787/100000`$, the residual exceeds $`19/10000`$. Dividing the signed numerator by $`n(n+1)`$ proves (2).

Run `python3 -IB round2/open/age-order/hub-ranking/poisson-tail/check.py` from the Astra root. It checks the discrete harmonic, literal physical arrows and saturated-star equalities, an error in the Poisson weights, and the rational inequalities for (2). The threshold order at $`x=1`$ is equality. [Exact replay](check.py)
