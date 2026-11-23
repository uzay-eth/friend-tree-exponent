> Produced by Claude Fable 5.2.

# The degree tail of the nucleus has the same exponent

In the random friend tree, $`T_2`$ is one edge. Given $`T_n`$, pick $`V_n`$ uniformly among the $`n`$ vertices, pick $`W_n`$ uniformly among the neighbours of $`V_n`$, and attach a new leaf to $`W_n`$. Write $`D_v(n)`$ for the degree of $`v`$ in $`T_n`$ and $`C_n = \#\{v : D_v(n) \ge 2\}`$ for the number of core vertices. For $`k \ge 2`$, write

```math
N_k(n) := \#\{v : D_v(n) \ge k\}, \qquad N_2(n) = C_n.
```

Krapivsky and Redner (J. Stat. Mech. 2017, eq. (1)) report numerically that the number of vertices of degree $`k`$ behaves like $`N^\mu k^{-(1+\mu)}`$. The theorem below gives a cumulative, logarithmic form of this observation. Uniformly in $`k`$, the expected number of vertices of degree at least $`k`$ in the $`n`$-vertex tree is bounded below by a constant times the expected core size at scale $`n/(8k)`$ and bounded above by an explicit expression in $`\mathbb E[C_m]`$ for $`m \ge 2n/k`$. Part (d) compares this expected count with the expected core size at scale $`n/k`$ on the logarithmic scale only.

**Theorem D.** Write $`a_m := \mathbb E[C_m]`$.

* **(a)** For all $`n \ge 3`$, $`k \ge 2`$ and every integer $`m_0`$ with $`3 \le m_0 \le n`$,

  ```math
  \mathbb E[N_k(n)] \le a_{m_0} + \frac{2n}{k}\sum_{m=m_0}^{n-1}\frac{a_{m+1} - a_m}{m+1}.
  ```

* **(b)** For all $`k \ge 2`$, $`m' \ge 3`$ and $`n \ge 8m'k`$,

  ```math
  \mathbb E[N_k(n)] \ge \frac{a_{m'}}{44160}.
  ```

* **(c)** Suppose $`a m^\delta \le a_m \le A m^\lambda`$ for all $`m \ge 3`$, with constants $`a, A > 0`$ and $`0 < \delta \le \lambda < 1`$. Then for all $`k \ge 2`$ and $`n \ge 24k`$,

  ```math
  \frac{a}{44160}\left(\frac{n}{12k}\right)^\delta
  \le \mathbb E[N_k(n)]
  \le A\left(2 + \frac{1}{1-\lambda}\right)\left(\frac{3n}{k}\right)^\lambda.
  ```

* **(d)** If $`\log a_m / \log m \to \mu \in (0,1)`$ as $`m \to \infty`$, then

  ```math
  \frac{\log \mathbb E[N_k(n)]}{\log(n/k)} \to \mu
  \qquad \text{as } n/k \to \infty,
  ```

  uniformly in $`k \ge 2`$.

**Corollary.** The [cross-charging bounds](cross-charging.md) give the following unconditional version. Gautschi's inequality for $`0 < s < 1`$ is

```math
n(n+1)^{s-1} < \frac{\Gamma(n+s)}{\Gamma(n)} < n^s.
```

By this inequality, Claim 2 gives $`0.44\,m^{0.55} \le a_m \le 3.974\,m^{0.58}`$ for all $`m \ge 3`$. The lower coefficient uses $`0.5022\,(3/4)^{0.45} \approx 0.4412`$. Hence for all $`k \ge 2`$ and $`n \ge 24k`$,

```math
\frac{0.44}{44160}\left(\frac{n}{12k}\right)^{0.55}
\le \mathbb E[N_k(n)] \le 33\left(\frac{n}{k}\right)^{0.58},
```

where $`3.974\cdot(2 + 1/0.42)\cdot 3^{0.58} \approx 32.93`$.

Using the sharper certified exponents of [Gemini 4 Neon's bounds](bracket-motifs.md) and the [existence theorem by GPT-6.2 Astra](existence.md) as inputs, part (d) gives the same logarithmic exponent $`\mu`$ for the cumulative degree tail, between $`0.5650`$ and $`0.5665`$. The theorem concerns cumulative counts. The corresponding pointwise degree law remains open. A pointwise power law with this cumulative exponent would give density exponent $`\nu = 1 + \mu`$, for which they report $`1.566`$.

## Proofs

Throughout, $`\mathcal F_j`$ is the sigma-field of $`T_3, \ldots, T_j`$. For a vertex $`w`$ of $`T_j`$,

```math
\mathbb P(W_j = w \mid \mathcal F_j) = \frac{1}{j}\sum_{u \sim w}\frac{1}{D_u(j)},
```

because $`V_j`$ is uniform on the $`j`$ vertices and, given $`V_j = u`$, $`W_j`$ is uniform on the $`D_u(j)`$ neighbours of $`u`$.

**Lemma 1 (degrees grow at most linearly in mean).** Let $`m \ge 3`$ and let $`v`$ be a vertex of $`T_m`$. For all $`n \ge m`$,

```math
\begin{split}
\mathbb E[D_v(n) \mid \mathcal F_m] &\le D_v(m)\frac{n}{m}, \\
\mathbb E[D_v(n)(D_v(n)+1) \mid \mathcal F_m]
&\le D_v(m)(D_v(m)+1)\frac{n(n+1)}{m(m+1)}.
\end{split}
```

*Proof.* The degree $`D_v`$ increases by one exactly when $`W_j = v`$, and $`\mathbb P(W_j = v \mid \mathcal F_j) \le D_v(j)/j`$ since each of the $`D_v(j)`$ neighbours contributes at most $`1/j`$. Thus $`\mathbb E[D_v(j+1) \mid \mathcal F_j] \le D_v(j)(j+1)/j`$, and iterating gives the first bound.

Put $`Q = D(D+1)`$. When $`W_j = v`$, $`Q`$ increases by $`2(D_v(j)+1)`$, so

```math
\mathbb E[Q(j+1) \mid \mathcal F_j]
\le Q(j) + 2(D_v(j)+1)\frac{D_v(j)}{j} = Q(j)\frac{j+2}{j}.
```

The product $`\prod_{j=m}^{n-1}(j+2)/j`$ equals $`n(n+1)/(m(m+1))`$. ∎

**Proof of (a).** A vertex of degree at least $`k \ge 2`$ in $`T_n`$ is a core vertex of $`T_n`$. It is either already core in $`T_{m_0}`$ (at most $`C_{m_0}`$ vertices) or becomes core at a step $`m`$ with $`m_0 \le m \le n-1`$, by receiving the new leaf while it is a leaf of $`T_m`$. At most one vertex becomes core at each step, and $`C_{m+1} - C_m`$ is the indicator of that event. Such a vertex $`v`$ has $`D_v(m+1) = 2`$, so by Lemma 1 and Markov's inequality,

```math
\mathbb P(D_v(n) \ge k \mid \mathcal F_{m+1}) \le \frac{2n}{(m+1)k}.
```

Taking expectations and summing over $`m`$ gives (a). ∎

**Lemma 2 (take-off).** Let $`m \ge 3`$ and let $`u`$ be a core vertex of $`T_m`$ with $`\ell := L_u(m) \ge 3`$ leaf neighbours and $`c_u(m) \le 2`$ core neighbours. Then for every $`n \ge m`$,

```math
\mathbb P\!\left(L_u(n) \ge \frac{(\ell-1)n}{2m} \mid \mathcal F_m\right)
\ge \frac{(\ell-1)^2}{4(\ell+2)(\ell+3)} \ge \frac{1}{30}.
```

*Proof.* The leaf count $`L_u`$ increases by one when $`W_j = u`$, which has conditional probability $`(L_u(j) + \rho_u(j))/j`$ with $`\rho_u \ge 0`$ the contribution of the core neighbours. It decreases by one when one of the leaves of $`u`$ receives the new leaf, which has conditional probability $`L_u(j)/(jD_u(j)) \le 1/j`$. Hence

```math
\begin{split}
\mathbb E[L_u(j+1) - L_u(j) \mid \mathcal F_j] &\ge \frac{L_u(j)-1}{j}, \\
\mathbb E[L_u(j+1) - 1 \mid \mathcal F_j] &\ge (L_u(j)-1)\frac{j+1}{j}.
\end{split}
```

Taking conditional expectations given $`\mathcal F_m`$ and iterating yields $`\mathbb E[L_u(n)-1 \mid \mathcal F_m] \ge (\ell-1)n/m`$. The factors are positive for either sign of $`L_u-1`$. It follows that $`\mathbb E[L_u(n) \mid \mathcal F_m] \ge (\ell-1)n/m`$.

By Lemma 1, using $`D_u(m) \le \ell+2`$ and $`(n+1)/(m+1) \le n/m`$,

```math
\begin{split}
\mathbb E[L_u(n)^2 \mid \mathcal F_m]
&\le \mathbb E[D_u(n)(D_u(n)+1) \mid \mathcal F_m] \\
&\le (\ell+2)(\ell+3)\frac{n(n+1)}{m(m+1)}
\le (\ell+2)(\ell+3)\left(\frac{n}{m}\right)^2.
\end{split}
```

By the Paley–Zygmund inequality for the nonnegative variable $`L_u(n)`$,

```math
\begin{split}
\mathbb P\!\left(L_u(n) \ge \frac{\mathbb E[L_u(n) \mid \mathcal F_m]}{2} \mid \mathcal F_m\right)
&\ge \frac{1}{4}\frac{\mathbb E[L_u(n) \mid \mathcal F_m]^2}{\mathbb E[L_u(n)^2 \mid \mathcal F_m]} \\
&\ge \frac{(\ell-1)^2}{4(\ell+2)(\ell+3)}.
\end{split}
```

The last expression is increasing for $`\ell \ge 3`$ and equals $`4/120 = 1/30`$ at $`\ell = 3`$. ∎

Call a core vertex $`u`$ **eligible** if $`c_u \le 2`$ and $`L_u \ge 1`$, or if $`(L_u,c_u) = (0,2)`$ and $`u`$ has a core neighbour $`w`$ with $`(L_w,c_w) = (0,2)`$.

**Lemma 3 (a quarter of the core is eligible).** A tree with $`C \ge 1`$ core vertices has at least $`C/4`$ eligible vertices. This is the counting lemma used in [GPT-6.2 Astra's proof of existence of the mean exponent, Lemma 2.1](existence.md#2-eligible-vertices).

*Proof.* If $`C = 1`$, the core vertex has $`c = 0`$ and $`L \ge 2`$. Let $`C \ge 2`$. The core vertices span a tree. Let $`n_1`$ be the number of core vertices with $`c = 1`$ (they have $`L \ge 1`$, so they are eligible), $`n_3`$ the number with $`c \ge 3`$, and $`I`$ the number with $`(L,c) = (0,2)`$ that have no neighbour of that kind. An ineligible vertex has $`c \ge 3`$ or is counted in $`I`$. The $`I`$ vertices are pairwise non-adjacent and each has exactly two core edges, so $`2I \le C-1`$. In a tree $`n_3 \le n_1-2`$. Thus the number $`E`$ of eligible vertices satisfies

```math
E \ge C-n_3-\frac{C-1}{2}, \qquad E \ge n_1 \ge n_3+2.
```

Adding gives $`2E \ge (C+1)/2+2`$. ∎

**Lemma 4 (activation).** Let $`m' \ge 3`$ and let $`u`$ be an eligible core vertex of $`T_{m'}`$. Then

```math
\mathbb P(L_u(8m') \ge 3 \text{ and } c_u(8m') \le 2 \mid \mathcal F_{m'})
\ge \frac{1}{368}.
```

*Proof.* The choices $`V_j`$ for $`j \ge m'`$ are independent, with $`V_j`$ uniform on the vertices of $`T_j`$ and independent of $`\mathcal F_j`$. Let $`G`$ be the event that $`V_j \ne u`$ for all $`m' \le j < 8m'`$. Then

```math
\mathbb P(G) = \prod_{j=m'}^{8m'-1}\left(1-\frac{1}{j}\right)
= \frac{m'-1}{8m'-1} \ge \frac{2}{23}.
```

On $`G`$, the vertex $`u`$ never chooses a neighbour, so no leaf of $`u`$ is ever promoted during $`[m',8m')`$. The leaf count $`L_u`$ is nondecreasing, $`c_u`$ stays fixed, and each leaf of $`u`$ stays a leaf of $`u`$. Conditionally on $`G`$, the choices $`V_j`$ in that interval are still independent, and $`V_j`$ is uniform on the $`j-1`$ vertices other than $`u`$.

**Case A.** Suppose $`c_u \le 2`$ and $`L_u \ge 1`$. Fix a leaf $`x`$ of $`u`$ in $`T_{m'}`$. Each $`j`$ with $`V_j = x`$ gives $`W_j = u`$, adding one more leaf at $`u`$. Given $`G`$, the probability that $`V_j = x`$ for some $`j`$ in $`[m',3m')`$ is

```math
1-\prod_{j=m'}^{3m'-1}\frac{j-2}{j-1}
= 1-\frac{m'-2}{3m'-2} \ge \frac{2}{3},
```

and the probability that $`V_j = x`$ for some $`j`$ in $`[3m',8m')`$ is $`1-(3m'-2)/(8m'-2) \ge 5/8`$. The two events are independent given $`G`$. On their intersection with $`G`$, $`L_u(8m') \ge 3`$ and $`c_u(8m') = c_u(m') \le 2`$. The probability is at least $`(2/23)(2/3)(5/8) = 5/138`$.

**Case B.** Suppose $`(L_u,c_u) = (0,2)`$ and $`u`$ has a core neighbour $`w`$ with $`(L_w,c_w) = (0,2)`$. Let $`z`$ be the other neighbour of $`w`$. Given $`G`$, consider the event $`H_1`$ that $`V_j \ne z`$ for all $`j`$ in $`[m',2m')`$ and $`V_j = w`$ for some $`j`$ in $`[m',2m')`$. Its conditional probability is

```math
\begin{split}
\prod_{j=m'}^{2m'-1}\frac{j-2}{j-1} - \prod_{j=m'}^{2m'-1}\frac{j-3}{j-1}
&= \frac{m'-2}{2m'-2}\left(1-\frac{m'-3}{2m'-3}\right) \\
&= \frac{(m'-2)m'}{(2m'-2)(2m'-3)} \ge \frac{1}{4} \qquad (m' \ge 3).
\end{split}
```

On $`G`$ and $`H_1`$, let $`\tau`$ be the first $`j \ge m'`$ with $`V_j = w`$. Before $`\tau`$, neither $`u`$ nor $`z`$ has been chosen as $`V_j`$, so $`w`$ has received no leaf and still has exactly the two neighbours $`u`$ and $`z`$. Hence $`W_\tau = u`$ with conditional probability $`1/2`$. In that event, $`u`$ has a leaf $`x`$ (the vertex created at step $`\tau`$), which stays a leaf on $`G`$.

As in Case A, given $`G`$ the events that $`V_j = x`$ for some $`j`$ in $`[2m',4m')`$ and for some $`j`$ in $`[4m',8m')`$ have probabilities $`1-(2m'-2)/(4m'-2) \ge 1/2`$ and $`1-(4m'-2)/(8m'-2) \ge 1/2`$, respectively. They are independent of each other and of what happened before $`2m'`$. On the intersection of the events in this construction, $`L_u(8m') \ge 3`$ and $`c_u(8m') = 2`$. The probability is at least $`(2/23)(1/4)(1/2)(1/4) = 1/368`$. ∎

**Proof of (b).** Let $`m' \ge 3`$, $`k \ge 2`$ and $`n \ge 8m'k`$. For an eligible vertex $`u`$ of $`T_{m'}`$, apply Lemma 4 and then Lemma 2 at time $`m = 8m'`$ with $`\ell \ge 3`$, so that $`(\ell-1)n/(2m) \ge n/(8m') \ge k`$. This gives

```math
\mathbb P(D_u(n) \ge k \mid \mathcal F_{m'})
\ge \mathbb P(L_u(n) \ge k \mid \mathcal F_{m'})
\ge \frac{1}{368}\cdot\frac{1}{30}.
```

By Lemma 3 there are at least $`C_{m'}/4`$ eligible vertices in $`T_{m'}`$. Hence

```math
\mathbb E[N_k(n)] \ge \frac{\mathbb E[C_{m'}]}{4\cdot368\cdot30} = \frac{a_{m'}}{44160}.
```

∎

**Proof of (c).** For the upper bound, take $`m_0 = \lceil 2n/k \rceil`$ in (a). Then $`3 \le m_0 \le n`$ because $`48 \le 2n/k \le n`$. Put $`b_j = 1/j`$. If $`m_0 = n`$ (which happens exactly when $`k = 2`$), the sum in (a) is empty. If $`m_0 < n`$, summation by parts gives

```math
\sum_{m=m_0}^{n-1}(a_{m+1}-a_m)b_{m+1}
= a_n b_n - a_{m_0}b_{m_0+1} + \sum_{j=m_0+1}^{n-1} a_j(b_j-b_{j+1}).
```

In both cases, since all $`a_j > 0`$ and $`b_j-b_{j+1} = 1/(j(j+1)) \le j^{-2}`$,

```math
\begin{split}
\sum_{m=m_0}^{n-1}(a_{m+1}-a_m)b_{m+1}
&\le a_n b_n + \sum_{j=m_0+1}^{n-1} a_j(b_j-b_{j+1}) \\
&\le A n^{\lambda-1} + A\sum_{j>m_0} j^{\lambda-2} \\
&\le A m_0^{\lambda-1}\left(1+\frac{1}{1-\lambda}\right),
\end{split}
```

where the last step uses $`n \ge m_0`$, $`\lambda < 1`$ and $`\sum_{j>m_0} j^{\lambda-2} \le m_0^{\lambda-1}/(1-\lambda)`$. Thus

```math
\begin{split}
\mathbb E[N_k(n)]
&\le A m_0^\lambda + \frac{2n}{k} A m_0^{\lambda-1}\left(1+\frac{1}{1-\lambda}\right) \\
&\le A\left(2+\frac{1}{1-\lambda}\right)m_0^\lambda,
\end{split}
```

using $`2n/k \le m_0`$. Also $`m_0 \le 2n/k+1 \le 3n/k`$. For the lower bound, take $`m' = \lfloor n/(8k) \rfloor \ge 3`$ in (b). Then

```math
m' \ge \frac{n}{8k}-1 \ge \frac{2}{3}\frac{n}{8k} = \frac{n}{12k}.
```

∎

**Proof of (d).** For every $`\varepsilon`$ with $`0 < \varepsilon < \min(\mu,1-\mu)`$, there are $`a_\varepsilon, A_\varepsilon > 0`$ with

```math
a_\varepsilon m^{\mu-\varepsilon} \le a_m \le A_\varepsilon m^{\mu+\varepsilon}
```

for all $`m \ge 3`$ ($`a_m \ge 1`$). Apply (c) with $`\delta = \mu-\varepsilon > 0`$ and $`\lambda = \mu+\varepsilon < 1`$, take logarithms, divide by $`\log(n/k)`$, let $`n/k \to \infty`$, then let $`\varepsilon \to 0`$. ∎

## Remarks

1. Parts (a) and (b) are statements about means. For each fixed $`k\ge2`$, the almost-sure comparability of successive cumulative counts in Theorem 4.6 of the published paper by Addario-Berry et al. (Theorem 4.8 in the preprint), iterated from the core count, also gives its almost-sure logarithmic exponent. The corresponding uniform result for growing $`k`$ remains open. [Random Friend Trees, Theorem 4.6](https://doi.org/10.1002/rsa.70082) · [Sample-path exponent of the core](almost-sure-exponent.md)
2. The constant $`1/44160`$ is far from optimal. Only its positivity matters for the exponent.
3. As a numerical illustration, in one run with $`C = 7.4 \times 10^6`$ at $`n = 1.0 \times 10^{12}`$ ([archived histogram](../simulations/outputs/hist_t27.txt)), the fractions of core vertices with degree at least $`k`$ were $`0.616, 0.338, 0.206, 0.131, 0.114`$ for $`k = 3, 6, 12, 24, 30`$, with local log-slopes of $`-0.64`$ between $`k = 16`$ and $`24`$ and $`-0.62`$ between $`24`$ and $`30`$. These are finite-size slopes from a single tree.
