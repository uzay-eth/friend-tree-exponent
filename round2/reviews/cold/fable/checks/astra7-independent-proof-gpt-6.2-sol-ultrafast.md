# Theorem 7: independent proof

**1. PROVED — Transition probabilities and drift.** Write $`\mathcal F_n=\sigma(T_2, \ldots, T_n)`$, and let $`p_n(v)=\mathbb P(W_n=v\mid\mathcal F_n)`$. The histories are recoverable from $`T_n`$ by deleting later labels, so this is also the conditional law given $`T_n`$. For a present vertex $`v`$, the sampling rule gives
```math
np_n(v)=\sum_{u\sim v}1/D_u(n), \qquad
0\le np_n(v)\le D_v(n), \qquad
0\le D_v(n)-np_n(v)\le C_n.
```
Indeed each leaf neighbour contributes zero to the difference, and each core neighbour contributes at most one. For absent vertices put $`p_n(v)=0`$, so these inequalities also hold there. For a present vertex, $`D_v(n+1)=D_v(n)+\mathbf 1_{\{W_n=v\}}`$. Thus, writing $`X_v(n)=D_v(n)/n`$,
```math
\mathbb E[X_v(n+1)\mid\mathcal F_n]-X_v(n)
=\frac{p_n(v)-X_v(n)}{n+1}\le0.
```

**2. STANDARD — Limit (a) and supermartingale (e).** For fixed $`v`$, on the deterministic index range $`n\ge \max(2, v)`$, $`X_v(n)`$ is an adapted nonnegative supermartingale bounded by one. The nonnegative supermartingale convergence theorem gives a finite almost sure limit $`Z_v`$; a countable intersection gives all the limits simultaneously. The one-step inequality holds throughout the lifetime of the vertex, in particular on $`\{\tau_v\le n\}`$ where $`\tau_v`$ is its first core time. If read on the clock shifted by the promotion, the same inequality holds on $`\{\tau_v<\infty\}`$ with the shifted filtration: split any conditional-expectation test event at shifted time $`k`$ into the events $`\{\tau_v=t\}`$, on each of which the event is in $`\mathcal F_{t+k}`$; then apply the deterministic-time inequality and sum. This proves (e).

**3. INPUT (I1) — Core estimates used.** Use only $`\mathbb E C_n=n^{\mu+o(1)}`$, $`\log C_n/\log n\to\mu`$ almost surely, and $`\mu<1`$ (from the supplied upper bound). In particular, with $`\alpha=(1+\mu)/2<1`$, $`\mathbb E C_n\le n^\alpha`$ eventually and $`C_n/n\to0`$ almost surely, since eventually $`C_n\le n^\alpha`$ there as well.

**4. PROVED — Immigration and the tail drift.** A tree here starts with one edge and gains one edge at every step. Counting edge ends gives $`\sum_v D_v(n)=2(n-1)`$. Its $`n-C_n`$ leaves account for one degree each; the core therefore accounts for $`n-2+C_n`$. For $`K\ge2`$ and $`n\ge K`$, let $`S_n^K`$ be the sum of degrees of core vertices with labels $`>K`$. Initially $`S_K^K=0`$. Let $`q_n=\mathbb P(W_n\text{ is a leaf}\mid\mathcal F_n)`$. The core count increases by exactly one if and only if $`W_n`$ is a leaf, so $`\mathbb E q_n=\mathbb E C_{n+1}-\mathbb E C_n`$. The new vertex adds nothing to $`S_n^K`$. The receiving vertex adds one if it was a late core vertex, and two if it was a late leaf. Hence
```math
\begin{aligned}
\mathbb E[S_{n+1}^K-S_n^K\mid\mathcal F_n]
&=\sum_{v>K, \,D_v(n)\ge2}p_n(v)+2\sum_{v>K, \,D_v(n)=1}p_n(v)\\ &\le S_n^K/n+2q_n,\\ \mathbb E[S_{n+1}^K/(n+1)-S_n^K/n\mid\mathcal F_n]&\le 2q_n/(n+1).
\end{aligned}
```

**5. PROVED — Vanishing expected tail bound.** Put $`c_n=\mathbb E C_n`$. Summing the expectation of the last inequality gives for all $`N\ge K`$
```math
\mathbb E[S_N^K/N]\le R_K:=2\sum_{n=K}^{\infty}\frac{c_{n+1}-c_n}{n+1}.
```
For $`M\ge K`$, shifting indices and collecting coefficients gives
```math
\sum_{n=K}^M\frac{c_{n+1}-c_n}{n+1}
=\frac{c_{M+1}}{M+1}-\frac{c_K}{K+1}+\sum_{j=K+1}^{M}\frac{c_j}{j(j+1)}.
```
The first term tends to zero by step 3. The last series converges by comparison with the convergent $`p`$-series of terms $`j^{-(2-\alpha)}`$. Since the increments of $`c_n`$ are nonnegative, $`0\le R_K\le2\sum_{j>K} c_j/[j(j+1)]\to0`$.

**6. STANDARD — No escaping mass, (b).** A fixed degree differs from its contribution to core degree by at most one. Therefore the finite sum of the first $`K`$ core shares tends to $`\sum_{v\le K}Z_v`$. The total core share is $`(n-2+C_n)/n\to1`$ by steps 3–4. Subtract to obtain
```math
S_n^K/n\longrightarrow 1-\sum_{v\le K} Z_v\ge0.
```
These assertions hold for every $`K`$ on the joint probability-one event from steps 2–3. Fatou's lemma for the nonnegative $`S_n^K/n`$ and step 5 bound the expected limit by $`R_K`$. The nonnegative random variable $`1-\sum_{v\ge1}Z_v`$ is bounded above by each such limit, so its expectation is zero. This proves (b).

**7. PROVED — Total variation, (c).** Step 1 and $`C_n/n\to0`$ give pointwise $`p_n(v)\to Z_v`$. The former sums to one by the sampling rule, and the latter by (b). On the joint probability-one event, for any $`K`$,
```math
\sum_v|p_n(v)-Z_v|\le\sum_{v\le K}|p_n(v)-Z_v|+\Bigl(1-\sum_{v\le K}p_n(v)\Bigr)+\Bigl(1-\sum_{v\le K}Z_v\Bigr).
```
The limit superior is at most $`2(1-\sum_{v\le K}Z_v)`$, which tends to zero. In particular the absolute probability difference on any set of labels, bounded by the sum on the left, tends to zero uniformly over sets.

**8. PROVED — Weak convergence, (d).** For a bounded continuous function $`g`$ on $`[0,1]`$, the conditional expectation is $`\sum_v p_n(v)g(X_v(n))`$. Subtract $`\sum_v Z_vg(Z_v)`$ by adding and subtracting $`\sum_v Z_vg(X_v(n))`$. The absolute value is at most
```math
\|g\|_\infty \sum_v|p_n(v)-Z_v|+\sum_{v\le K}Z_v|g(X_v(n))-g(Z_v)|+2\|g\|_\infty\sum_{v>K}Z_v.
```
The first two terms vanish for fixed $`K`$ by (a), (c) and continuity. The third vanishes as $`K\to\infty`$ by (b). The event used is independent of $`g`$, so this gives weak convergence to the stated measure.

**9. PROVED — A protected leaf-neighbour share.** For a present vertex $`v`$, write $`L_v(n)`$ for its number of leaf neighbours. If it receives the new leaf this count increases by one. If a leaf neighbour receives the new leaf, the count decreases by one; otherwise there is no change. These events are disjoint. A leaf neighbour $`u`$ of $`v`$ has $`p_n(u)=1/(nD_v(n))`$. Since each leaf neighbour contributes one to $`np_n(v)`$,
```math
\mathbb E[L_v(n+1)-L_v(n)\mid\mathcal F_n]
=p_n(v)-\frac{L_v(n)}{nD_v(n)}\ge\frac{L_v(n)}{n}-\frac1n,
```
where the last step is $`L_v(n)\le D_v(n)`$. Subtracting the denominators explicitly, this gives
```math
\mathbb E[(L_v(n+1)-1)/(n+1)\mid\mathcal F_n]\ge (L_v(n)-1)/n.
```

**10. STANDARD — Positive conditional chance of positive share.** For a deterministic $`m\ge\max(2, v)`$ and any $`A\in\mathcal F_m`$, repeated conditional expectation in step 9 yields for $`N\ge m`$
```math
\mathbb E[\mathbf 1_A X_v(N)]\ge\mathbb E[\mathbf 1_A(L_v(N)-1)/N]\ge\mathbb E[\mathbf 1_A(L_v(m)-1)/m].
```
Since $`0\le X_v(N)\le1`$ and (a) holds, bounded convergence gives the left limit $`\mathbb E[\mathbf 1_A Z_v]`$. In particular, if $`A\subseteq\{L_v(m)\ge2\}`$ has positive probability, then $`\mathbb E[\mathbf 1_A Z_v]\ge\mathbb P(A)/m>0`$, and $`\mathbb P(A\cap\{Z_v>0\})>0`$. In fact this argument with all such $`A`$ shows the conditional probability of positivity is at least $`1/m`$ on $`\{L_v(m)\ge2\}`$, since $`Z_v\le1`$.

**11. PROVED — Positivity for every label at arbitrarily late promotion, (f).** Fix $`v`$ and any integer $`t\ge\max(2, v)`$. Choose an initial label $`r\in\{1,2\}\setminus\{v\}`$. Require all steps from size 2 through size $`t-1`$ to attach to $`r`$, then require the steps at sizes $`t, t+1`$ to attach to $`v`$. Every present vertex has a neighbour of positive degree, so step 1 gives strictly positive probability to each prescribed choice; by multiplying conditional probabilities the finite event has positive probability. Through size $`t`$, vertex $`v`$ is a leaf. It is promoted exactly at size $`t+1`$; at size $`t+2`$ it has two new leaf neighbours. Apply step 10 at $`m=t+2`$ to the prescribed history to obtain $`\mathbb P(\tau_v=t+1, \ Z_v>0)>0`$. Moreover, on any history that promotes $`v`$ at size $`t+1`$, it already has one new leaf neighbour and has strictly positive conditional probability to receive the very next leaf, so the same test-event argument gives a positive conditional chance of $`Z_v>0`$ on every such history.

## Inputs used

- **I1:** $`\mathbb E C_n=n^{\mu+o(1)}`$, $`\log C_n/\log n\to\mu`$ almost surely, and $`\mu<1`$. No proof files were consulted.
- **I2:** basic conditional expectation (including its tower and test-event properties), nonnegative supermartingale convergence, Fatou's lemma, bounded convergence, and countable intersections of probability-one events.

## Gaps

None in (a)–(f).
