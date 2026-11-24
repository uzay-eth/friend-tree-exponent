# Astra Claim 7 — Hub mass and the law of the attachment recipient

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 7: Claude and Gemini; GPT standalone audit; inputs include sublinear core growth and fixed-label degree limits. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Model and claim

Start from the single edge $`T_2`$ on labels 1 and 2. Given $`T_n`$, choose $`V_n`$ uniformly from its $`n`$ vertices, choose $`W_n`$ uniformly from the neighbours of $`V_n`$, and attach the new vertex $`n+1`$ to $`W_n`$ by one edge. Let $`D_v(n)`$ be the degree of label $`v`$ in $`T_n`$, extended by zero before birth. Write $`\mathcal F_n`$ for the full history through the construction of $`T_n`$, before the choices at size $`n`$. A **core** vertex has degree at least two. Let $`C_n`$ count the core vertices. A **promotion** occurs when a leaf first reaches degree two.

**Claim.** Simultaneously for all labels, $`D_v(n)/n`$ has an almost-sure limit $`Z_v`$, and almost surely

```math
\sum_{v\ge1}Z_v=1,\qquad
\big\|\mathcal L(W_n\mid T_n)-(Z_v)_{v\ge1}\big\|_{\mathrm{TV}}\longrightarrow0,
\qquad
\mathcal L\!\left(D_{W_n}(n)/n\mid T_n\right)
\ \Longrightarrow\ \sum_{v\ge1}Z_v\,\delta_{Z_v}.
```

Here total variation on the labels is half the $`\ell^1`$ distance. The last limit is weak convergence of random probability measures on $`[0,1]`$, equivalently convergence in the Prokhorov metric. The “empirical law” here is the conditional law for the next attachment at size $`n`$: it weights vertices by their recipient probabilities, not uniformly. It is neither an average over past times nor a claim that the realized sequence $`D_{W_n}(n)/n`$ converges almost surely.

## 1. One-step estimates and the core-degree identity

For $`v\le n`$, conditioning on the uniform source gives

```math
p_n(v):=\mathbb P(W_n=v\mid\mathcal F_n)
       =\mathbb P(W_n=v\mid T_n)
       =\frac1n\sum_{u\sim v}\frac1{D_u(n)}\le\frac{D_v(n)}n.
\qquad\text{(1)}
```

Set $`p_n(v)=0`$ for $`v>n`$. For a vertex already present, its degree increases by $`\mathbf1_{\{W_n=v\}}`$. Thus, starting at the deterministic size $`\max(2,v)`$, its normalized degree is a nonnegative supermartingale:

```math
\mathbb E\!\left[\frac{D_v(n+1)}{n+1}\mid\mathcal F_n\right]
=\frac{D_v(n)+p_n(v)}{n+1}\le\frac{D_v(n)}n.
\qquad\text{(2)}
```

The supermartingale convergence theorem gives the fixed-label limits, on a single probability-one event by countability. These are the degree limits in Addario-Berry et al.; we do not need to import their leaf-neighbour limits.

At the founder edge both vertices are leaves, so $`C_2=0`$ and $`C_3=1`$ deterministically. At every size $`j\ge3`$, the sole neighbour of a leaf is a core vertex. If $`L_u(j)`$ denotes the number of leaf neighbours of $`u`$, then

```math
\mathbb E[C_{j+1}-C_j\mid\mathcal F_j]
=\mathbb P(W_j\text{ is a leaf}\mid\mathcal F_j)
=\frac1j\sum_{u:\,D_u(j)\ge2}\frac{L_u(j)}{D_u(j)}\le\frac{C_j}{j}.
\qquad\text{(3)}
```

There is at most one promotion per step, since the new vertex is born a leaf and only the recipient's degree increases among the existing vertices. The last inequality in (3) is not used at $`j=2`$.

Define the normalized core degrees $`q_n(v)=D_v(n)\mathbf1_{\{D_v(n)\ge2\}}/n`$. The tree has total degree $`2(n-1)`$, and its $`n-C_n`$ leaves each contribute one. Hence

```math
\sum_v q_n(v)=\frac{n+C_n-2}{n},\qquad
0\le D_v(n)/n-q_n(v)\le 1/n.
\qquad\text{(4)}
```

In particular $`q_n(v)\to Z_v`$ for each fixed $`v`$. If $`C_n/n\to0`$, then the core degree mass divided by $`n`$ tends to one. The total degree divided by $`n`$ tends to two: the other unit comes from the leaves, each of which has normalized degree $`1/n`$ while it is a leaf.

## 2. Mean bound and tightness of late-labelled cores

Let $`a_n=\mathbb E C_n`$. The quantitative input is [Theorem 1 in §1 of `proofs/bracket-motifs.md`][motifs], at commit `34b6d52e050465d1d7e160edb61d0eda0d544ff2`. It gives, for $`n\ge3`$,

```math
a_n\le\frac{14{,}606{,}753{,}031}{2{,}989{,}767{,}809}
\prod_{j=3}^{n-1}\left(1+\frac{1133}{2000j}\right)
\le K n^\alpha,\qquad \alpha=1133/2000<1,
\qquad\text{(5)}
```

for some finite deterministic $`K`$; the last bound follows from $`\log(1+x)\le x`$ and the harmonic sum estimate. No mean-growth exponent needs to exist for this argument. Equation (5) gives $`C_n/n\to0`$ almost surely: for each fixed $`\varepsilon>0`$, Markov's inequality bounds $`\mathbb P(C_{2^k}/2^k>\varepsilon)`$ by $`(K/\varepsilon)2^{-k(1-\alpha)}`$ for $`k\ge2`$. Borel–Cantelli for a countable sequence of $`\varepsilon`$ decreasing to zero yields the dyadic limit. For $`2^k\le n<2^{k+1}`$, monotonicity gives $`C_n/n\le 2C_{2^{k+1}}/2^{k+1}`$.

Small core counts alone do not turn pointwise convergence into conservation of mass: a core with a late label could still carry a large degree. To control this, let $`\sigma_v=\inf\{r\ge3:D_v(r)\ge2\}`$, with $`\inf\varnothing=\infty`$. On $`\{\sigma_v=r\}`$, an event in $`\mathcal F_r`$, $`D_v(r)=2`$. For $`r\le n`$ and $`v\le r`$, iteration of (2) and conditioning on the history **after** the promotion give

```math
\mathbb E\!\left[\mathbf1_{\{\sigma_v=r\}}\frac{D_v(n)}n\mid\mathcal F_r\right]
\le\frac2r\mathbf1_{\{\sigma_v=r\}}.
\qquad\text{(6)}
```

This handles the random promotion time by conditioning on each possible size $`r`$; the estimate starts at $`T_r`$, after the step $`T_{r-1}\to T_r`$. For $`n\ge m\ge3`$, a label $`v>m`$ that is a core at size $`n`$ must have been promoted at some $`m<r\le n`$. At each such $`r`$, the number of these promotions is at most $`C_r-C_{r-1}`$. Summing (6) and using (3) and (5),

```math
\begin{aligned}
\mathbb E\sum_{v>m}q_n(v)
&\le 2\sum_{r=m+1}^n\frac{a_r-a_{r-1}}r
 \le 2K\sum_{j=m}^{\infty}\frac{j^{\alpha-1}}{j+1}\\
&\le 2K\left(m^{\alpha-2}+\frac{m^{\alpha-1}}{1-\alpha}\right)=:b_m.
\end{aligned}
\qquad\text{(7)}
```

The last bound separates the first term and integrates the decreasing function $`x^{\alpha-2}`$ for the remaining terms. In particular $`b_m\to0`$ uniformly as a bound over $`n\ge m`$.

By (4) and $`C_n/n\to0`$, every finite sum of the $`Z_v`$ is at most one almost surely; hence $`S:=\sum_v Z_v\le1`$. For fixed $`m`$, bounded convergence gives $`\mathbb E\sum_{v\le m}q_n(v)\to\mathbb E\sum_{v\le m}Z_v`$. The expected total in (4) tends to one by (5). Subtracting the tail bound (7) therefore gives $`\mathbb E\sum_{v\le m}Z_v\ge1-b_m`$. Monotone convergence as $`m\to\infty`$ gives $`\mathbb E S\ge1`$. Since $`S\le1`$ almost surely, $`S=1`$ almost surely.

## 3. Alternate proof of mass conservation from an almost-sure bound

The same mass conclusion follows from $`C_n=O(n^\beta)`$ almost surely for any deterministic $`0<\beta<1`$, without a bound on $`a_n`$. For this tree, Proposition 6.22 in both versions of Addario-Berry et al. gives $`n^{-\beta}C_n\to0`$ almost surely for $`\sqrt2-1/2<\beta<1`$. The argument below uses that pathwise input directly.

Fix $`K\ge1`$ and let $`\tau=\inf\{r\ge3:C_r>K r^\beta\}`$. At size $`m\ge3`$, start with no counted vertices. At step $`T_{r-1}\to T_r`$ for $`r>m`$, if $`\tau\ge r`$, update the normalized degrees of the already counted vertices to size $`r`$ and add any vertex promoted on that step with weight $`2/r`$. Call the resulting sum at size $`n`$ $`X_n`$. This includes the stopping step $`r=\tau`$ and freezes the entire sum thereafter, including its denominator. If $`\tau\le m`$, the sum stays zero. The active indicator $`\mathbf1_{\{\tau\ge r\}}`$ is $`\mathcal F_{r-1}`$-measurable, so the old contributions have nonpositive conditional drift by (2). The new contribution is exactly $`2(B_r-B_{r-1})/r`$, where $`B_r=C_{r\wedge\tau}`$. Consequently

```math
\mathbb E X_n\le 2\,\mathbb E\sum_{r=m+1}^n\frac{B_r-B_{r-1}}r.
\qquad\text{(8)}
```

For $`r<\tau`$, $`B_r\le Kr^\beta`$. If $`r\ge\tau\ge4`$, then $`B_r=C_\tau\le C_{\tau-1}+1\le K(\tau-1)^\beta+1\le Kr^\beta+1`$. If $`\tau=3`$, the same bound follows from $`C_3=1`$. This unit overshoot is why we can keep the stopping increment. For $`n>m`$, summation by parts now gives the pathwise bound

```math
\begin{aligned}
\sum_{r=m+1}^n\frac{B_r-B_{r-1}}r
&=\frac{B_n}{n}-\frac{B_m}{m+1}
 +\sum_{r=m+1}^{n-1}\frac{B_r}{r(r+1)}\\
&\le (K+1)\left(1+\frac1{1-\beta}\right)m^{\beta-1}.
\end{aligned}
```

On $`\{\tau=\infty\}`$, $`X_n\ge\sum_{v>m}q_n(v)`$ because all such cores were promoted after $`m`$. On that event, $`C_n/n\to0`$ and (4) implies $`\sum_{v>m}q_n(v)\to 1-\sum_{v\le m}Z_v\ge0`$. Fatou's lemma and (8) bound the expectation of $`\mathbf1_{\{\tau=\infty\}}(1-\sum_{v\le m}Z_v)`$ by $`2(K+1)(1+1/(1-\beta))m^{\beta-1}`$. First let $`m\to\infty`$ at fixed $`K`$, by bounded convergence on $`\{\tau=\infty\}`$, to obtain $`S=1`$ there almost surely. Then let integer $`K\to\infty`$: the events $`\{\tau=\infty\}`$ increase to a probability-one event under the almost-sure power bound. That bound also gives $`C_n/n\to0`$ almost surely, as needed below.

## 4. Conditional laws

For $`n\ge3`$, the core induces a tree: every internal vertex on a path between two core vertices also has degree at least two. Let $`c_v(n)`$ be the number of core neighbours of a core vertex $`v`$, and set it to zero off the core. Then $`\sum_v c_v(n)=2(C_n-1)`$. The leaf-source part of (1) is $`L_v(n)/n`$; it vanishes on leaves and equals $`(D_v(n)-c_v(n))/n`$ on cores. The core-source part of (1) is a nonnegative vector of total mass $`C_n/n`$, because each core source has probability $`1/n`$ and sends it all to its neighbours. Thus

```math
\|p_n-q_n\|_1\le\frac{C_n+2(C_n-1)}n\longrightarrow0\qquad\text{almost surely}.
\qquad\text{(9)}
```

In particular $`p_n(v)\to Z_v`$ for each label. Since both $`p_n`$ and $`Z`$ have total mass one, for each fixed $`m`$ the tail mass $`\sum_{v>m}p_n(v)`$ tends to $`\sum_{v>m}Z_v`$. The $`\ell^1`$ distance on the first $`m`$ labels tends to zero, and the distance on the tail is bounded by the sum of the two tail masses. Taking $`n\to\infty`$ and then $`m\to\infty`$ proves the total-variation limit. This is why there is no attachment dust: arbitrarily large fractions of recipient probability are eventually carried by finitely many labels of positive limiting mass, even though the leaves collectively carry one unit of normalized degree.

Finally set $`d_n(v)=D_v(n)/n`$. For every continuous $`f:[0,1]\to\mathbb R`$, the conditional expectation of $`f(d_n(W_n))`$ is $`\sum_v p_n(v)f(d_n(v))`$. Its difference from $`\sum_v Z_v f(Z_v)`$ in absolute value is at most

```math
\|f\|_\infty\|p_n-Z\|_1+\sum_v Z_v\,|f(d_n(v))-f(Z_v)|.
```

The first term tends to zero by total variation. The second tends to zero by dominated convergence against the fixed probability vector $`Z`$ on each path with simultaneous fixed-label convergence. This proves weak convergence on the same probability-one event. In particular, taking $`f(x)=x`$ gives the conditional mean limit $`\sum_v Z_v^2`$; bounded convergence gives $`\mathbb E[D_{W_n}(n)/n]\to\mathbb E\sum_v Z_v^2`$.

## References and verification

- Astra source: [`round2/hub-mass.md`][source] and the [motif growth bound][motifs], both at `34b6d52e050465d1d7e160edb61d0eda0d544ff2`.
- Louigi Addario-Berry, Simon Briend, Luc Devroye, Serte Donderwinkel, Céline Kerriou and Gábor Lugosi, *Random Friend Trees*, **Random Structures & Algorithms** 69(1), e70082 (2026), [DOI 10.1002/rsa.70082][doi] ([author-hosted copy of the journal version][journal]). In that version the fixed-label limits are Proposition 3.1 (§3, p. 4); the two conjectures appear in §3, p. 5, after the discussion of Proposition 3.1 and immediately before Proposition 3.7 (the conditional-mean result), respectively. In [arXiv:2403.20185v1, 29 March 2024][preprint], those results are Theorems 3.1 and 3.3, and the conjectures appear in the corresponding paragraphs in §3. The almost-sure core bound is Proposition 6.22 in both versions.
- The pinned [cold review][review] checks both proof routes in the source. The [finite checker][checker] enumerates histories through size seven, checks the degree and recipient identities and the late-core estimates, and rejects five planted errors. These finite checks do not establish the infinite-size limits.
- The standalone cold audit independently checks post-promotion conditioning and prefix-wise stopped drift, including crossing-step predictability, the frozen denominator and unit overshoot, and uses an exact double star to distinguish recipient, uniform-vertex and degree mass and test the label-$`\ell^1`$ bound.

[source]: ../hub-mass.md
[motifs]: ../../proofs/bracket-motifs.md#theorem-1-finite-size-bounds-for-the-mean-core-size
[review]: ../../../../../reviews/cold/astra/checks/hub-mass-carbon.md
[checker]: ../checks/hub_mass_carbon.py
[doi]: https://doi.org/10.1002/rsa.70082
[journal]: https://luc.devroye.org/Addario+Briend+Devroye+Donderwinkel+Kerriou+Lugosi-Random_friendship_trees-2026-.pdf
[preprint]: https://arxiv.org/html/2403.20185v1#S3
