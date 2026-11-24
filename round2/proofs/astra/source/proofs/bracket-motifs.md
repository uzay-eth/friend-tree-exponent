
> **Author/model:** Gemini 4 Neon. **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../SOURCES.tsv)).
> **Review and dependencies:** Round 1 prerequisite. This release grades only the downstream uses; see the dependencies and scoped reviews in the [Astra ledger](../../../../LEDGER.md#astra).

> Produced by Gemini 4 Neon. The title and introduction have been edited to distinguish the numerical estimate from the exact-value question. Sections 2 to 9 have been reorganised for readability; the potential, the certificates and the checker are unchanged.

# Product bounds for the mean number of non-leaves

## 1. Model and Main Theorem

Consider the **random friend tree** (Krapivsky–Redner isotropic redirection model with $`r = 1`$, arXiv:1706.01514 and arXiv:2403.20185):
- Start at $`n = 2`$ from a single edge $`T_2`$ on two vertices.
- Given $`T_n`$, pick a vertex $`V_n \in V(T_n)`$ uniformly at random, pick a neighbour $`W_n`$ of $`V_n`$ uniformly at random, and attach a new leaf $`n+1`$ to $`W_n`$ to form $`T_{n+1}`$.
- Let $`C_n := \#\{v \in V(T_n) : \deg_{T_n}(v) \ge 2\}`$ be the number of non-leaves (**core** vertices) of $`T_n`$.

Krapivsky and Redner measured a growth exponent near $`0.566`$. Strogatz later asked whether $`17/30`$ could be its exact value. The following bounds on $`\mathbb{E}[C_n]`$ place the upper exponent below that fraction. [Krapivsky and Redner, §5.1](https://arxiv.org/pdf/1706.01514#page=12) · [Strogatz’s post](https://x.com/stevenstrogatz/status/1572287647715934210)

### Theorem 1 (Finite-size bounds for the mean core size)
For every integer $`n \ge 3`$,
```math
\frac{7{,}304{,}352{,}308}{8{,}763{,}936{,}965} \prod_{j=3}^{n-1}\left(1 + \frac{113}{200\, j}\right) \;\le\; \mathbb{E}[C_n] \;\le\; \frac{14{,}606{,}753{,}031}{2{,}989{,}767{,}809} \prod_{j=3}^{n-1}\left(1 + \frac{1133}{2000\, j}\right).
```
Consequently,
```math
0.5650 \;=\; \frac{113}{200} \;\le\; \liminf_{n\to\infty}\frac{\log \mathbb{E}[C_n]}{\log n} \;\le\; \limsup_{n\to\infty}\frac{\log \mathbb{E}[C_n]}{\log n} \;\le\; \frac{1133}{2000} \;=\; 0.5665 \;<\; \frac{17}{30} = 0.5666\overline{6}.
```
The upper exponent differs from $`17/30`$ by exactly $`17/30-1133/2000=1/6000`$.

### Corollary 2 (Almost-Sure Upper Bound Strictly Below $`17/30`$)
Almost surely,
```math
\limsup_{n\to\infty} \frac{\log C_n}{\log n} \;\le\; \frac{1133}{2000} \;=\; 0.5665 \;<\; \frac{17}{30}.
```
*Proof.* Fix any $`\varepsilon > 0`$ and put $`n_k = 2^k`$. By Theorem 1, $`\mathbb{E}[C_{n_k}] = O(n_k^{1133/2000})`$ since $`\prod_{j=3}^{n-1}(1 + \alpha/j) = \frac{2}{\Gamma(3+\alpha)}\frac{\Gamma(n+\alpha)}{\Gamma(n)} \sim \frac{2}{\Gamma(3+\alpha)} n^\alpha`$. By Markov's inequality,
```math
\mathbb{P}\!\left(C_{n_k} > n_k^{1133/2000 + \varepsilon}\right) \le \frac{\mathbb{E}[C_{n_k}]}{n_k^{1133/2000 + \varepsilon}} = O(2^{-k\varepsilon}),
```
which is summable over $`k \ge 2`$. By Borel–Cantelli, a.s. $`C_{n_k} \le n_k^{1133/2000 + \varepsilon}`$ for all sufficiently large $`k`$. Because $`C_n`$ is nondecreasing in $`n`$, for $`n_k \le n < n_{k+1} = 2 n_k`$ we have $`C_n \le C_{n_{k+1}} \le (2n)^{1133/2000+\varepsilon}`$, so $`\limsup_{n\to\infty} \frac{\log C_n}{\log n} \le \frac{1133}{2000} + \varepsilon`$ a.s. Intersecting over rational $`\varepsilon \downarrow 0`$ yields the claim. $`\square`$

*(Remark: Theorem 1 and Corollary 2 are completely self-contained and do not assume existence of $`\lim \frac{\log \mathbb{E}[C_n]}{\log n}`$. Combined with the independent existence and almost-sure convergence theorems of [`existence.md`](existence.md) and [`almost-sure-exponent.md`](almost-sure-exponent.md), which show that a deterministic $`\mu`$ exists with $`\mathbb{E}[C_n] = n^{\mu+o(1)}`$ and $`\frac{\log C_n}{\log n} \to \mu`$ a.s., Theorem 1 yields $`\mu \in \left[\frac{113}{200}, \frac{1133}{2000}\right] = [0.5650, 0.5665]`$.)*

---

## 2. Dynamics in Core Coordinates

### 2.1 Rooting, core and local states
At $`n = 3`$ the tree $`T_3`$ is the star $`K_{1,2}`$. Its centre $`r`$ is the unique vertex of degree $`2`$. For every $`n \ge 3`$ we root $`T_n`$ at $`r`$.

**Lemma 2.1 (structure, $`n \ge 3`$).**
1. $`\deg(r) \ge 2`$, so $`r \in K(T_n) := \{v : \deg_{T_n}(v) \ge 2\}`$.
2. Every non-root core vertex $`u`$ has a core parent $`p_u`$. Hence $`K(T_n)`$ is a connected rooted subtree containing $`r`$.
3. Every leaf has exactly one neighbour, its parent, and that parent is a core vertex.

*Proof.* Degrees never decrease and $`\deg_{T_3}(r) = 2`$, which gives (1). If $`u \ne r`$, its parent $`p_u`$ is either $`r`$ or has a parent of its own as well as the child $`u`$, so $`\deg(p_u) \ge 2`$. That gives (2). A leaf has degree $`1`$. It is not $`r`$, so its only neighbour is its parent, which is core by (2). $`\square`$

For $`u \in K := K(T)`$ write $`h_u := \mathbf{1}_{u \ne r}`$, $`L_u :=`$ the number of leaf children of $`u`$, $`k_u :=`$ the number of core children of $`u`$, and $`D_u := \deg(u) = L_u + k_u + h_u \ge 2`$. The triple $`s_u := (L_u, k_u, h_u)`$ is the **state** of $`u`$. Write $`\mathrm{ch}_K(u)`$ for the core children of $`u`$, $`\iota_u := 1/D_u`$, and
```math
\rho_u := \sum_{v \in N_K(u)} \iota_v = h_u\,\iota_{p_u} + \sum_{c \in \mathrm{ch}_K(u)} \iota_c, \qquad r_u := L_u + \rho_u, \qquad b_u := \frac{L_u}{D_u}.
```
Here $`N_K(u)`$ is the set of core neighbours of $`u`$; it has $`h_u + k_u`$ elements. Since every core degree is at least $`2`$,
```math
0 \le \rho_u \le \frac{h_u + k_u}{2}, \qquad L_u \le r_u \le D_u - \frac{h_u + k_u}{2} < D_u + 1 .
```
For a state $`s = (L, k, h)`$ put $`s^+ := (L+1, k, h)`$ and, when $`L \ge 1`$, $`s^{\mathrm{bud}} := (L-1, k+1, h)`$. Let $`\mathrm{nb} := (1, 0, 1)`$.

### 2.2 One step of the process
**Lemma 2.2 (exact one-step law).** Let $`T_n = T`$ with $`n \ge 3`$. Exactly one of the following happens at step $`n \to n+1`$.
1. **Leaf gain at a core vertex $`u`$**, with probability $`r_u/n`$. The new leaf becomes a child of $`u`$, so $`s_u \mapsto s_u^+`$ and no other state changes.
2. **Budding at a core vertex $`u`$ with $`L_u \ge 1`$**, with probability $`b_u/n`$. The new leaf is attached to one of the $`L_u`$ leaf children $`\ell`$ of $`u`$, each with probability $`1/(nD_u)`$. Then $`\ell`$ becomes a new core vertex in state $`\mathrm{nb}`$, $`s_u \mapsto s_u^{\mathrm{bud}}`$ (so $`D_u`$ is unchanged), and no other state changes.

The probabilities add up to one: $`\sum_{u \in K}(r_u + b_u) = n`$.

*Proof.* The new leaf is attached to $`w`$ with probability $`\frac{1}{n}\sum_{v \sim w} \iota_v`$. For core $`w`$, its leaf neighbours contribute $`1`$ each and its core neighbours contribute $`\rho_w`$, which gives $`r_w`$. A leaf $`\ell`$ has the single neighbour $`p_\ell`$ (Lemma 2.1(3)), which gives $`\iota_{p_\ell}`$. The state changes are read off directly. The total is $`\frac{1}{n}\sum_w \sum_{v \sim w}\iota_v = \frac{1}{n}\sum_v \deg(v)\,\iota_v = 1`$. $`\square`$

**Lemma 2.3 (monotonicity).** Along the process, $`h_u`$ is constant, $`D_u`$ and $`k_u`$ never decrease, and core vertices stay core. Hence $`C_n`$ is nondecreasing. For the sets $`H_u`$ and $`C_u`$ of §3 (each is defined by $`D < \text{const}`$, $`k \le \text{const}`$ and a lower bound on $`k`$, and $`H_u`$ also fixes $`h = 1`$):
- a vertex that leaves $`H_u`$ or $`C_u`$ never re-enters;
- a vertex enters $`H_u`$ only by budding from $`k = 0`$ to $`k = 1`$, and enters $`C_u`$ only by budding from $`k = 1`$ to $`k = 2`$.

*Proof.* Leaf gain raises $`D`$ by one and keeps $`k`$. Budding keeps $`D`$ and raises $`k`$ by one. Membership in $`H_u`$ needs $`k \ge 1`$, and membership in $`C_u`$ needs $`k \ge 2`$. The upper bounds can only be broken, never restored, and the lower bound on $`k`$ can only be crossed by budding, one unit at a time. $`\square`$

### 2.3 Scaled generator
For any functional $`\Phi`$ of rooted trees, Lemma 2.2 gives
```math
(\mathcal{L}\Phi)(T) := n\,\mathbb{E}\big[\Phi(T_{n+1}) - \Phi(T_n) \,\big|\, T_n = T\big] = \sum_{u \in K(T)} \Big[r_u\,\Delta^{L}_u \Phi(T) + b_u\,\Delta^{B}_u \Phi(T)\Big],
```
where $`\Delta^L_u\Phi`$ and $`\Delta^B_u\Phi`$ are the exact changes of $`\Phi`$ under leaf gain and under budding at $`u`$. The budding change $`\Delta^B_u\Phi`$ includes the terms of $`\Phi`$ that the newborn vertex creates.

---

## 3. The Motif Potential $`\Phi`$

### 3.1 Parameters
The potential uses eight integer parameters. The conditions below are the ones the argument uses, each tagged with where it is used.

| Parameter | Shipped value | Condition | Used in |
|---|---:|---|---|
| $`L_F`$ | 800 | $`D_0 \le L_F`$ | §5.2, f-cases (C) and (D): $`L = L_F`$ implies a tail state |
| $`K_F`$ | 150 | $`D_0 \le K_F + 1`$ | §5.2, f-case (E), and §6: an overflow state is a tail state |
| $`D_0`$ | 24 | $`3 \le D_0`$ | $`\mathrm{nb}`$ and $`T_3`$'s root are inner |
| $`K_T`$ | 12 | $`K_T + 1 \le K_F`$, $`K_T \ge 0`$ | clamp $`\kappa = \min(k, K_T+1)`$ |
| $`D_{hu}, K_{hu}`$ | 16, 8 | $`2 \le D_{hu} \le D_0`$, $`K_{hu} \ge 1`$ | middle vertices of $`h_V`$ are inner |
| $`D_{cu}, K_{cu}`$ | 12, 6 | $`2 \le D_{cu} \le D_0`$, $`K_{cu} \ge 2`$ | apex vertices of $`h_C`$ are inner |

The verifier tests only some of these conditions on the header. The shipped headers satisfy all of them, for example $`D_0 = 24 \le K_F + 1 = 151`$. See `bracket-motifs-code-map.md`, §3.

### 3.2 State sets
- A state $`(L, k, h)`$ is **inner** if $`D = L + k + h < D_0`$ and **tail** otherwise. $`\mathcal{I}`$ is the set of inner states with $`D \ge 2`$.
- **Tail endpoints:** $`\mathrm{INF}`$ (standing for $`D = \infty`$) and $`T_1(\kappa, h)`$ for $`\kappa \in \{0, \dots, K_T+1\}`$, $`h \in \{0, 1\}`$ (standing for $`D = D_0`$ with clamped child count $`\kappa`$).
- $`\mathcal{G} := \mathcal{I} \cup \{\mathrm{INF}\} \cup \{T_1(\kappa, h)\}`$. The **parent states** $`\mathcal{P} \subset \mathcal{G}`$ are $`\mathrm{INF}`$ and the inner or $`T_1`$ states with $`k \ge 1`$ (resp. $`\kappa \ge 1`$). The **child states** $`\mathcal{W} \subset \mathcal{G}`$ are $`\mathrm{INF}`$ and the inner or $`T_1`$ states with $`h = 1`$.
- **Middle and apex domains:**
  $`\displaystyle H_u := \{(L, k, 1) : 1 \le k \le K_{hu},\ L + k + 1 < D_{hu}\}, \qquad C_u := \{(L, k, h) : 2 \le k \le K_{cu},\ L + k + h < D_{cu}\}.`$
  Both are subsets of $`\mathcal{I}`$, because $`D_{hu}, D_{cu} \le D_0`$.
- **Vertex projection:** $`\pi_f(L, k, h) := \mathrm{VSTAR}_h`$ if $`k > K_F`$, and $`(\min(L, L_F), k, h)`$ otherwise.
- For $`v \in \mathcal{G}`$ put $`\mathrm{inv}(v) := 1/D_v`$ if $`v \in \mathcal{I}`$, $`\mathrm{inv}(\mathrm{INF}) := 0`$ and $`\mathrm{inv}(T_1(\kappa, h)) := 1/D_0`$.

### 3.3 Tail interpolation
For a core vertex $`v`$ define a probability vector $`\mu_v`$ on $`\mathcal{G}`$:
- if $`D_v < D_0`$, then $`\mu_v := \delta_{s_v}`$;
- if $`D_v \ge D_0`$, then with $`t_v := D_0/D_v \in (0, 1]`$ and $`\kappa_v := \min(k_v, K_T+1)`$,
  $`\displaystyle \mu_v := (1 - t_v)\,\delta_{\mathrm{INF}} + t_v\,\delta_{T_1(\kappa_v, h_v)} .`$

For a tail vertex, $`\iota_v = t_v/D_0`$. That is exactly the $`\mu_v`$-average of $`\mathrm{inv}`$. For an inner vertex, $`\iota_v = \mathrm{inv}(s_v)`$. So in all cases
```math
\iota_v = \sum_{a} \mu_v(a)\,\mathrm{inv}(a).
```
We refer to this identity as (3.1).

### 3.4 Tables and the potential
The certificate supplies rational tables $`f`$ on $`\mathrm{dom}(\pi_f)`$, $`g`$ on $`\mathcal{P} \times \mathcal{W}`$ with $`g[a, b] = 0`$ whenever $`a, b \notin \mathcal{I}`$, $`h_V`$ on $`\mathcal{P} \times H_u \times \mathcal{W}`$, and $`h_C`$ on $`C_u \times \mathcal{W} \times \mathcal{W}`$, symmetric in its last two arguments. These are evaluated on core vertices by
```math
g(p \to w) := \sum_{a, b} \mu_p(a)\mu_w(b)\,g[a, b], \qquad h_V(p \to u \to w) := \mathbf{1}_{s_u \in H_u}\sum_{a, b} \mu_p(a)\mu_w(b)\,h_V[a, s_u, b],
```
```math
h_C(u; \{w_1, w_2\}) := \mathbf{1}_{s_u \in C_u}\sum_{a, b} \mu_{w_1}(a)\mu_{w_2}(b)\,h_C[s_u, a, b],
```
and
```math
\Phi(T) := \sum_{u \in K} f(\pi_f(s_u)) + \sum_{p \to w} g(p \to w) + \sum_{p \to u \to w} h_V(p \to u \to w) + \sum_{u \in K}\ \sum_{\{w_1, w_2\} \subset \mathrm{ch}_K(u)} h_C(u; \{w_1, w_2\}).
```
The second sum runs over core edges and the third over vertical core paths of length two. Every summand is a **term**. A term $`\tau`$ has a **support** $`S(\tau)`$ of one, two or three core vertices, and its value depends only on the states of $`S(\tau)`$. If both ends of an edge are tail vertices, every pair $`(a, b)`$ in the sum has $`a, b \notin \mathcal{I}`$, so the term is $`0`$.

The certificate also supplies centring constants $`r_{\mathrm{par}}, r_{\mathrm{ch}}, x_{\mathrm{par}}, x_{\mathrm{ch}}, x_{cc}`$ on $`\mathcal{I}`$, and knot values for three piecewise-linear functions $`A_{pc}[s], A_{cp}[s], A_{cc}[s]`$ on $`[0, 1/2]`$ for each $`s \in \mathcal{I}`$. All three use the knots $`0 = \xi_0 < \xi_1 = 1/D_0 < \dots < \xi_{m-1} = 1/2`$.

---

## 4. The Exact Drift Identity

This section has no inequalities. It rewrites $`(\mathcal{L}\Phi)(T)`$ exactly as a sum of one **local expression** per term.

### 4.1 Event-wise form
For a term $`\tau`$ and $`x \in S(\tau)`$, let $`\Delta^L_x\tau`$ be the change of $`\tau`$ when $`s_x`$ is replaced by $`s_x^+`$. Let $`\Delta^B_x\tau`$ be the change when $`s_x`$ is replaced by $`s_x^{\mathrm{bud}}`$, with $`\Delta^B_x\tau := 0`$ if $`L_x = 0`$. These changes include the indicators $`\mathbf{1}_{H_u}`$ and $`\mathbf{1}_{C_u}`$: a path or cherry whose middle or apex leaves its domain jumps to $`0`$. For a core vertex $`x`$ with $`L_x \ge 1`$ let
```math
N_x := f(\mathrm{nb}) + g(s_x^{\mathrm{bud}} \to \mathrm{nb}) + h_x\,h_V(p_x \to s_x^{\mathrm{bud}} \to \mathrm{nb}) + \sum_{c \in \mathrm{ch}_K(x)} h_C(s_x^{\mathrm{bud}}; \{c, \mathrm{nb}\}),
```
where each motif is evaluated with the current states of $`p_x`$ and $`c`$.

**Identity (I0).**
```math
(\mathcal{L}\Phi)(T) = \sum_{\tau} \sum_{x \in S(\tau)} \Big[r_x\,\Delta^L_x\tau + b_x\,\Delta^B_x\tau\Big] + \sum_{x \in K} b_x\,N_x .
```

*Proof.* By Lemma 2.2, an event at $`x`$ changes only the state of $`x`$, plus the newborn in the case of budding. A term changes only if $`x`$ is in its support, and then it changes by $`\Delta^L_x\tau`$ or $`\Delta^B_x\tau`$. Leaf gain creates no new core vertex, so it creates no new term. Budding at $`x`$ creates the core vertex $`\mathrm{nb}`$, which has no core child. The new terms are therefore exactly the terms whose support contains $`\mathrm{nb}`$:
- $`f(\mathrm{nb})`$;
- the edge $`x \to \mathrm{nb}`$;
- the path $`p_x \to x \to \mathrm{nb}`$ when $`h_x = 1`$ ($`\mathrm{nb}`$ can only be the bottom of a path);
- the cherries $`\{c, \mathrm{nb}\}`$ at apex $`x`$, for $`c \in \mathrm{ch}_K(x)`$.

These are the four parts of $`N_x`$. Could a path or cherry with an *old* support become nonzero because its middle or apex enters $`H_u`$ or $`C_u`$? By Lemma 2.3 that happens only when $`x`$ itself buds from $`k = 0`$ (no earlier child, so no old path through $`x`$) or from $`k = 1`$ (one earlier child, so no old cherry at $`x`$). Multiply by the probabilities of Lemma 2.2 and by $`n`$. $`\square`$

### 4.2 Allocation of neighbour credit
In (I0) each rate is $`r_x = L_x + \sum_{v \in N_K(x)} \iota_v`$. Each product $`\iota_v\,\Delta^L_x\tau`$ is kept in the local expression of $`\tau`$, or moved to a neighbouring term, according to this table. "Kept" means it stays in $`\tau`$. A centring constant $`c`$ splits a product as $`\iota_v\Delta = c\Delta + (\iota_v - c)\Delta`$. A Fenchel–Young (FY) split with function $`A`$ writes it as $`\iota_v\Delta = A(\iota_v) + R`$ with $`R := \iota_v\Delta - A(\iota_v)`$.

| Term $`\tau`$ | Vertex $`x`$ | Neighbour $`v`$ of $`x`$ | Allocation of $`\iota_v\,\Delta^L_x\tau`$ |
|---|---|---|---|
| $`f(u)`$ | $`u`$ tail | any | kept |
| $`f(u)`$ | $`u`$ inner | $`p_u`$ | $`r_{\mathrm{par}}[s_u]\Delta`$ kept; $`(\iota_{p_u} - r_{\mathrm{par}}[s_u])\Delta`$ to $`E(p_u \to u)`$ |
| $`f(u)`$ | $`u`$ inner | $`c \in \mathrm{ch}_K(u)`$ | $`r_{\mathrm{ch}}[s_u]\Delta`$ kept; $`(\iota_c - r_{\mathrm{ch}}[s_u])\Delta`$ to $`E(u \to c)`$ |
| $`g(p \to w)`$ | $`x \in \{p, w\}`$ tail | any | kept |
| $`g(p \to w)`$ | $`p`$ inner | $`w`$ | kept |
| $`g(p \to w)`$ | $`p`$ inner | $`p_p`$, $`s_p \in H_u`$ | $`x_{\mathrm{par}}[s_p]\Delta`$ kept; $`(\iota_{p_p} - x_{\mathrm{par}}[s_p])\Delta`$ to $`P(p_p \to p \to w)`$ |
| $`g(p \to w)`$ | $`p`$ inner | $`p_p`$, $`s_p \notin H_u`$ | FY with $`A_{pc}[s_p]`$: $`R`$ kept; $`A_{pc}[s_p](\iota_{p_p})`$ to $`E(p_p \to p)`$ |
| $`g(p \to w)`$ | $`p`$ inner | sibling $`c \ne w`$, $`s_p \in C_u`$ | $`x_{cc}[s_p]\Delta`$ kept; $`(\iota_c - x_{cc}[s_p])\Delta`$ to $`Ch(p; \{w, c\})`$ |
| $`g(p \to w)`$ | $`p`$ inner | sibling $`c \ne w`$, $`s_p \notin C_u`$ | FY with $`A_{cc}[s_p]`$: $`R`$ kept; $`A_{cc}[s_p](\iota_c)`$ to $`E(p \to c)`$ |
| $`g(p \to w)`$ | $`w`$ inner | $`p`$ | kept |
| $`g(p \to w)`$ | $`w`$ inner | $`c \in \mathrm{ch}_K(w)`$, $`s_w \in H_u`$ | $`x_{\mathrm{ch}}[s_w]\Delta`$ kept; $`(\iota_c - x_{\mathrm{ch}}[s_w])\Delta`$ to $`P(p \to w \to c)`$ |
| $`g(p \to w)`$ | $`w`$ inner | $`c \in \mathrm{ch}_K(w)`$, $`s_w \notin H_u`$ | FY with $`A_{cp}[s_w]`$: $`R`$ kept; $`A_{cp}[s_w](\iota_c)`$ to $`E(w \to c)`$ |
| $`h_V`$, $`h_C`$ | any | any | kept |

Conditions such as "$`p_p`$ exists" ($`h_p = 1`$), "a sibling exists" ($`k_p \ge 2`$) and "a child exists" ($`k_w \ge 1`$) are implicit. When they fail, the row is empty.

### 4.3 The local expressions
With this allocation the local expressions are as follows. Write $`[\cdot]`$ for an indicator, and $`\Delta^L_x g`$, $`\Delta^B_x g`$ for the jumps of the term being described.

**Vertex term, $`u \in K`$.** Let $`\theta_u := h_u r_{\mathrm{par}}[s_u] + k_u r_{\mathrm{ch}}[s_u]`$ if $`u`$ is inner, and $`\theta_u := \rho_u`$ if $`u`$ is tail. Then
```math
V(u) := (L_u + \theta_u)\,\Delta^L_u f + b_u\Big(f(\pi_f(s_u^{\mathrm{bud}})) - f(\pi_f(s_u)) + f(\mathrm{nb}) + g(s_u^{\mathrm{bud}} \to \mathrm{nb})\Big).
```

**Edge term, $`p \to w`$.** Let $`\sigma_p`$ and $`\sigma_w`$ be the kept rates. If $`p`$ is inner, $`\sigma_p := L_p + \iota_w + [s_p \in H_u]\,x_{\mathrm{par}}[s_p] + [s_p \in C_u](k_p - 1)\,x_{cc}[s_p]`$. If $`w`$ is inner, $`\sigma_w := L_w + \iota_p + [s_w \in H_u]\,k_w x_{\mathrm{ch}}[s_w]`$. For a tail end, the kept rate is the full rate ($`\sigma_p := r_p`$, resp. $`\sigma_w := r_w`$). Then
```math
\begin{aligned}
E(p \to w) :=\;& [p\ \mathrm{inner}](\iota_w - r_{\mathrm{ch}}[s_p])\Delta^L_p f + [w\ \mathrm{inner}](\iota_p - r_{\mathrm{par}}[s_w])\Delta^L_w f \\
&+ \sigma_p\,\Delta^L_p g + b_p\big(\Delta^B_p g + [s_p^{\mathrm{bud}} \in C_u]\,h_C(s_p^{\mathrm{bud}}; \{w, \mathrm{nb}\})\big) \\
&+ \sigma_w\,\Delta^L_w g + b_w\big(\Delta^B_w g + [s_w^{\mathrm{bud}} \in H_u]\,h_V(p \to s_w^{\mathrm{bud}} \to \mathrm{nb})\big) \\
&+ [p\ \mathrm{inner}, h_p = 1, s_p \notin H_u]\big(\iota_{p_p}\Delta^L_p g - A_{pc}[s_p](\iota_{p_p})\big) \\
&+ [p\ \mathrm{inner}, s_p \notin C_u] \sum_{c \in \mathrm{ch}_K(p) \setminus \{w\}} \big(\iota_c\,\Delta^L_p g - A_{cc}[s_p](\iota_c)\big) \\
&+ [w\ \mathrm{inner}, s_w \notin H_u] \sum_{c \in \mathrm{ch}_K(w)} \big(\iota_c\,\Delta^L_w g - A_{cp}[s_w](\iota_c)\big) \\
&+ [w\ \mathrm{inner}, s_w \notin H_u]\,k_w\,A_{pc}[s_w](\iota_p) + [p\ \mathrm{inner}, h_p = 1, s_p \notin H_u]\,A_{cp}[s_p](\iota_w) \\
&+ [p\ \mathrm{inner}, s_p \notin C_u](k_p - 1)\,A_{cc}[s_p](\iota_w).
\end{aligned}
```
The last two lines are the FY pieces this edge **receives**. It gets one $`A_{pc}[s_w](\iota_p)`$ from each path $`p \to w \to c`$ ($`k_w`$ of them), one $`A_{cp}[s_p](\iota_w)`$ from the path $`p_p \to p \to w`$, and one $`A_{cc}[s_p](\iota_w)`$ from each sibling $`c \ne w`$ ($`k_p - 1`$ of them). The two budding terms in the second and third lines are the newborn cherries and paths of $`N_p`$ and $`N_w`$, booked on this edge.

**Vertical path, $`p \to u \to w`$ with $`s_u \in H_u`$.**
```math
P(p \to u \to w) := (\iota_p - x_{\mathrm{par}}[s_u])\,\Delta^L_u g(u \to w) + (\iota_w - x_{\mathrm{ch}}[s_u])\,\Delta^L_u g(p \to u) + \sum_{x \in \{p, u, w\}} \big(r_x\,\Delta^L_x h_V + b_x\,\Delta^B_x h_V\big).
```

**Cherry, $`u`$ with $`s_u \in C_u`$, $`\{w_1, w_2\} \subset \mathrm{ch}_K(u)`$.**
```math
Ch(u; \{w_1, w_2\}) := (\iota_{w_1} - x_{cc}[s_u])\,\Delta^L_u g(u \to w_2) + (\iota_{w_2} - x_{cc}[s_u])\,\Delta^L_u g(u \to w_1) + \sum_{x \in \{u, w_1, w_2\}} \big(r_x\,\Delta^L_x h_C + b_x\,\Delta^B_x h_C\big).
```

**Identity (I).** For every rooted tree $`T`$ arising at some $`n \ge 3`$,
```math
(\mathcal{L}\Phi)(T) = \sum_{u \in K} V(u) + \sum_{p \to w} E(p \to w) + \sum_{p \to u \to w,\ s_u \in H_u} P(p \to u \to w) + \sum_{u:\ s_u \in C_u}\ \sum_{\{w_1, w_2\}} Ch(u; \{w_1, w_2\}).
```

*Proof.* Start from (I0) and expand each $`r_x`$ as $`L_x`$ plus its neighbour credits. Every product $`\iota_v\,\Delta^L_x\tau`$ in (I0) is listed in exactly one row of the table in §4.2. Every row replaces it by pieces that add back to it: $`c\Delta + (\iota_v - c)\Delta = \iota_v\Delta`$, and $`A(\iota_v) + R = \iota_v\Delta`$. The pieces are placed as follows.
- $`V(u)`$ keeps $`L_u + h_u r_{\mathrm{par}} + k_u r_{\mathrm{ch}}`$ for inner $`u`$. The parent edge and each of the $`k_u`$ child edges receive one residual.
- The kept parts of $`g(p \to w)`$ are $`\sigma_p`$, $`\sigma_w`$ and the three $`R`$ sums.
- The residuals of the centring splits $`x_{\mathrm{par}}, x_{\mathrm{ch}}`$ (resp. $`x_{cc}`$) are the two source terms of $`P`$ (resp. $`Ch`$). Each path $`p \to u \to w`$ with $`s_u \in H_u`$ receives exactly one residual of each kind, one from edge $`u \to w`$ (neighbour $`p`$) and one from edge $`p \to u`$ (neighbour $`w`$). Each cherry receives one from each of its two edges.
- The $`A(\iota_v)`$ pieces go to the edges listed in the last two lines of $`E`$, with the multiplicities stated there.
- The budding parts $`b_x\Delta^B_x\tau`$ and the four parts of $`b_x N_x`$ are placed as follows: $`f(\mathrm{nb}) + g(s_x^{\mathrm{bud}} \to \mathrm{nb})`$ in $`V(x)`$, the newborn path in $`E(p_x \to x)`$, and each newborn cherry $`\{c, \mathrm{nb}\}`$ in $`E(x \to c)`$.

Nothing is counted twice or left out. $`\square`$

**Check against the certificate (independent of the verifier).** The [identity test](../verify/tree-drift/motif_identity.py) implements this section independently of the production checker. It computes $`(\mathcal{L}\Phi)(T)`$ by brute force over all $`n`$ attachment targets, computes all local expressions exactly, and compares. Identity (I) held exactly on every tree tested (§8).

---

## 5. Bounding the Local Expressions in Both Directions

Fix a direction. **UPPER** means $`\lambda = 1133/2000`$, $`\mathop{\mathrm{ext}} = \max`$ and $`\lessgtr`$ is $`\le`$. **LOWER** means $`\lambda = 113/200`$, $`\mathop{\mathrm{ext}} = \min`$ and $`\lessgtr`$ is $`\ge`$. The goal of this section is
```math
V(u) \lessgtr \lambda f(\pi_f(s_u)), \quad E(p \to w) \lessgtr \lambda g(p \to w), \quad P \lessgtr \lambda h_V, \quad Ch \lessgtr \lambda h_C \qquad \text{for every term of every tree.}
```
We refer to these four inequalities together as (5.0).
By Identity (I), (5.0) implies $`(\mathcal{L}\Phi)(T) \lessgtr \lambda\Phi(T)`$.

### 5.1 Tools
**Lemma 5.1 (piecewise-linear Fenchel–Young).** Let $`A`$ interpolate $`(\xi_j, a_j)_{j<m}`$ with $`\xi_0 = 0`$, $`\xi_1 = 1/D_0`$, $`\xi_{m-1} = 1/2`$. Put $`A^{*}_{\max}(y) := \max_j(\xi_j y - a_j)`$ and $`A^{*}_{\min}(y) := \min_j(\xi_j y - a_j)`$. Then:
1. For $`x \in [0, 1/2]`$ and all $`y`$: $`A^*_{\min}(y) \le xy - A(x) \le A^*_{\max}(y)`$.
2. $`A^*_{\max}`$ is convex and $`A^*_{\min}`$ is concave.
3. $`A(t/D_0) = (1-t)A(0) + tA(1/D_0)`$ for $`t \in [0, 1]`$.

*Proof.* Write $`x = (1-\alpha)\xi_i + \alpha\xi_{i+1}`$. Then $`xy - A(x) = (1-\alpha)(\xi_i y - a_i) + \alpha(\xi_{i+1}y - a_{i+1})`$, which proves (1). (2) holds because a max (min) of affine functions is convex (concave). (3) holds because $`[0, 1/D_0]`$ is one knot interval. $`\square`$

So for the FY remainders of §4.3, $`R = \iota_v\Delta - A(\iota_v)`$ satisfies $`R \lessgtr A^*_{\mathop{\mathrm{ext}}}(\Delta)`$, with $`\Delta`$ the relevant jump. $`A^*_{\mathop{\mathrm{ext}}}`$ is convex for UPPER and concave for LOWER.

**Lemma 5.2 (tail vertices).** Let $`v`$ be a tail vertex with $`\kappa := \kappa_v`$, $`h := h_v`$. For any $`F[\cdot]`$ on $`\mathcal{G}`$ write $`F(\mu_v) := \sum_a \mu_v(a)F[a]`$.
1. $`r_v\,\Delta^L_v F(\mu_v) = t_v\,\gamma_v\,(F[\mathrm{INF}] - F[T_1(\kappa, h)])`$ with $`\gamma_v := r_v/(D_v+1) \in [\gamma_{\min}(\kappa, h), 1]`$. Here $`\gamma_{\min}(\kappa, h) := (D_0 - \kappa - h)/(D_0 + 1)`$ for $`\kappa \le K_T`$, and $`\gamma_{\min}(K_T + 1, h) := 0`$.
2. If $`k_v \le K_T`$, then $`b_v\,\Delta^B_v F(\mu_v) = t_v\,b_v\,(F[T_1(\kappa+1, h)] - F[T_1(\kappa, h)])`$ with $`b_v \in [b_{\min}(\kappa, h), 1]`$, where $`b_{\min}(\kappa, h) := (D_0 - \kappa - h)/D_0`$. If $`k_v > K_T`$, then $`\Delta^B_v F(\mu_v) = 0`$.

*Proof.* (1) Leaf gain keeps $`\kappa`$ and $`h`$ and sends $`t_v \mapsto D_0/(D_v+1)`$, a change of $`-t_v/(D_v+1)`$. The bound $`\gamma_v < 1`$ comes from $`r_v < D_v + 1`$ (§2.1). For the lower bound, $`\gamma_v \ge L_v/(D_v+1) = (D_v - k_v - h)/(D_v + 1)`$. If $`k_v = \kappa \le K_T`$, this is nondecreasing in $`D_v`$ because $`k_v + h + 1 \ge 0`$, so it is at least its value at $`D_v = D_0`$. Otherwise $`0`$ is a valid bound. (2) Budding keeps $`D_v`$, hence $`t_v`$, and changes $`\kappa`$ only if $`k_v \le K_T`$. Also $`b_v = 1 - (k_v + h)/D_v \ge 1 - (\kappa + h)/D_0`$. $`\square`$

**Lemma 5.3 (quadratic credit).** For $`t \in [0, 1]`$ and reals $`d_0, d_1`$:
```math
\frac{t}{D_0}\min(d_0, d_1) \le \frac{t}{D_0}\big((1-t)d_0 + t d_1\big) \le \frac{t}{D_0}\max(d_0, d_1).
```
With a second parameter $`t'`$ and bi-affine $`d_{ij}`$, apply this for each fixed $`t'`$.

**Lemma 5.4 (corner reduction).** Let the tail vertices of a support carry parameters $`t = (t_1, \dots, t_q)`$, $`q \le 2`$. Let $`\omega_c(t) \ge 0`$, $`\sum_c \omega_c = 1`$, be the product weights $`\prod_i \mu_{v_i}`$ over the corners $`c \in \{\mathrm{INF}, T_1\}^q`$. Suppose a local expression satisfies $`X(t) \lessgtr \sum_c \omega_c(t)\,\widehat X(c)`$, while the term is exactly $`\tau(t) = \sum_c \omega_c(t)\,\tau[c]`$. If $`\widehat X(c) \lessgtr \lambda\,\tau[c]`$ at every corner, then $`X(t) \lessgtr \lambda\tau(t)`$. $`\square`$

The bound $`X(t) \lessgtr \sum_c \omega_c(t)\widehat X(c)`$ comes from the following facts. Each one holds in the UPPER form; the LOWER form is the same with $`\max \leftrightarrow \min`$ and $`\le \leftrightarrow \ge`$.
- (a) Multi-affine pieces equal their corner interpolation. The source terms $`(\iota_p - x)\Delta^L_u g(u \to w)`$ are of this kind: by (3.1), $`\iota_p`$ is affine in $`t_p`$, and $`\Delta^L_u g(u \to w)`$ is affine in $`t_w`$. So are the transfers $`(\iota_w - r_{\mathrm{ch}})\Delta^L_p f`$ and the received pieces $`A(\iota_v)`$ for tail $`v`$, by Lemma 5.1(3).
- (b) A rate $`r`$ that is independent of $`t`$, lies in a known interval, and multiplies a multi-affine jump $`Y(t)`$ satisfies $`rY(t) = \sum_c \omega_c\,rY(c) \le \sum_c \omega_c \max_{r' \in \{r_{\mathrm{lo}}, r_{\mathrm{hi}}\}} r'Y(c)`$.
- (c) $`A^*_{\max}(Y(t)) \le \sum_c \omega_c A^*_{\max}(Y(c))`$, by convexity (Jensen).
- (d) Tail rates (Lemma 5.2): $`t_v\gamma_v Y(t') \le t_v\sum_{c'}\omega_{c'}(t')\max_{\gamma \in \{\gamma_{\min}, 1\}}\gamma Y(c')`$. This puts weight at the $`T_1`$ corner of $`v`$ and $`0`$ at its $`\mathrm{INF}`$ corner. Budding at tail vertices works the same way with $`b`$.
- (e) Products of a tail credit $`\iota_v = t_v/D_0`$ with a jump that depends on $`t_v`$ (the only non-multi-affine pieces) are bounded by Lemma 5.3. That again puts weight only at the $`T_1`$ corner of $`v`$.

### 5.2 All cases, family by family
Each line below gives a case of the local expression, the bound used, and the resulting finite inequality. Each finite inequality is checked by the verifier (map in `bracket-motifs-code-map.md`). "Corner" means a table state in $`\mathcal{G}`$; for tail endpoints, $`\mathrm{inv}`$ replaces $`\iota`$ as in (3.1).

**Vertex family $`V(u) \lessgtr \lambda f`$.** Cases by $`\pi_f(s_u)`$:
- (A) **Inner** ($`D < D_0`$, so $`L < L_F`$). $`V(u)`$ is a fixed rational:
  $`\displaystyle (L + h r_{\mathrm{par}} + k r_{\mathrm{ch}})\Delta^L f + B(s) \lessgtr \lambda f(s), \qquad B(s) := [L \ge 1]\tfrac{L}{D}\big(f(s^{\mathrm{bud}}) - f(s) + f(\mathrm{nb}) + g(s^{\mathrm{bud}} \to \mathrm{nb})\big).`$
  When $`s^{\mathrm{bud}}`$ is inner, $`g(s^{\mathrm{bud}} \to \mathrm{nb})`$ is a table entry, since $`s^{\mathrm{bud}}`$ has $`D`$ unchanged.
- (B) **Tail with $`L < L_F`$ and $`k \le K_F`$.** Here $`\theta_u = \rho_u \in [0, (h+k)/2]`$, and every other quantity is a fixed rational. $`g(s^{\mathrm{bud}} \to \mathrm{nb})`$ is evaluated exactly at $`t = D_0/D`$, by the verifier's `eval_g_q`. The two inequalities are those with $`r \in \{L, L + (h+k)/2\}`$ in place of $`L + \theta_u`$.
- (C) **$`L = L_F`$ exactly, $`k \le K_F`$.** This is a tail state because $`D \ge L_F \ge D_0`$. $`\Delta^L f = 0`$, since $`\pi_f`$ clamps $`L_F + 1 \mapsto L_F`$. So $`V = \frac{L_F}{D}(f(\pi_f(L_F - 1, k+1, h)) - f + f(\mathrm{nb}) + g(\cdot \to \mathrm{nb}))`$, a single exact inequality.
- (D) **$`L \ge L_F + 1`$, $`k \le K_F`$.** $`\Delta^L f = 0`$, and $`f(\pi_f(s^{\mathrm{bud}})) = f(\pi_f(L_F, k+1, h))`$, which is $`\mathrm{VSTAR}_h`$ if $`k = K_F`$. The free quantities are $`b = L/D \in [\frac{L_F+1}{L_F+1+k+h}, 1]`$ and $`t = D_0/D \in (0, \frac{D_0}{L_F+1+k+h}]`$. The expression $`b\cdot(\text{const} + \text{affine}(t))`$ is bilinear in $`(b, t)`$. It is checked at the four corners of the rectangle, which contains every attainable $`(b, t)`$.
- (E) **$`k > K_F`$ ($`\mathrm{VSTAR}_h`$).** $`\Delta^L f = \Delta^B f = 0`$, and $`V = b(f(\mathrm{nb}) + g(s^{\mathrm{bud}} \to \mathrm{nb}))`$. Since $`D \ge K_F + 1 + h \ge D_0`$ (this is where the condition $`D_0 \le K_F + 1`$ is needed), $`s^{\mathrm{bud}}`$ is tail with $`\kappa = K_T + 1`$, and $`t \le D_0/(K_F + 1 + h)`$. The checked rectangle is $`b \in [0, 1]`$, $`t \in [0, D_0/(K_F+1+h)]`$.

In (A) and (B) the vertex $`\mathrm{nb}`$ is inner ($`D = 2 < D_0`$). In every case $`s^{\mathrm{bud}}`$ has the same degree as $`s`$.

**Edge family $`E(p \to w) \lessgtr \lambda g`$.** Corners are $`(a, b) \in \mathcal{P} \times \mathcal{W}`$.
- **Both ends tail** ($`a, b \notin \mathcal{I}`$). Then $`g \equiv 0`$. Events keep both ends tail, so every jump of $`g`$ is $`0`$. Tail ends transfer nothing and receive nothing, and the newborn indicators vanish because tail states lie outside $`H_u`$ and $`C_u`$. Hence $`E \equiv 0 = \lambda g`$. The verifier checks $`g[a, b] = 0`$ for these pairs and skips their rows. This case needs no inequality.
- **At least one inner end.** The following summands are evaluated at the corner $`(a, b)`$; each item names its source in §5.1.
  - Transfers $`(\mathrm{inv}(b) - r_{\mathrm{ch}}[a])\Delta^L f(a)`$ for inner $`a`$, and $`(\mathrm{inv}(a) - r_{\mathrm{par}}[b])\Delta^L f(b)`$ for inner $`b`$ (fact (a)).
  - Received pieces $`A_{\mathrm{tot}}[a, b] := [b \in \mathcal{I}\setminus H_u, k_b \ge 1]\,k_b A_{pc}[b](\mathrm{inv}(a)) + [a \in \mathcal{I}\setminus H_u, h_a = 1]\,A_{cp}[a](\mathrm{inv}(b)) + [a \in \mathcal{I} \setminus C_u, k_a \ge 2](k_a - 1)A_{cc}[a](\mathrm{inv}(b))`$ (fact (a) with Lemma 5.1(3)).
  - Kept $`p`$-rate for inner $`a`$: $`(L_a + [b \in \mathcal{I}]\mathrm{inv}(b) + [a \in H_u]x_{\mathrm{par}}[a] + [a \in C_u](k_a - 1)x_{cc}[a])\,\Delta_p g`$, where $`\Delta_p g := g[a^+, b] - g[a, b]`$. If $`D_a = D_0 - 1`$, then $`a^+ = T_1(\kappa(k_a), h_a)`$.
  - Budding at inner $`a`$: $`\frac{L_a}{D_a}(g[a^{\mathrm{bud}}, b] - g[a, b] + [a^{\mathrm{bud}} \in C_u]h_C[a^{\mathrm{bud}}, b, \mathrm{nb}])`$ (exact).
  - The same three items for inner $`b`$ with $`\Delta_w g := g[a, b^+] - g[a, b]`$, kept rate $`L_b + [a \in \mathcal{I}]\mathrm{inv}(a) + [b \in H_u]k_b x_{\mathrm{ch}}[b]`$, and newborn path $`[b^{\mathrm{bud}} \in H_u]h_V[a, b^{\mathrm{bud}}, \mathrm{nb}]`$.
  - FY remainders: $`A^*_{pc}[a](\Delta_p g)`$ if $`a \in \mathcal{I} \setminus H_u`$, $`h_a = 1`$; $`(k_a - 1)A^*_{cc}[a](\Delta_p g)`$ if $`a \in \mathcal{I}\setminus C_u`$, $`k_a \ge 2`$; $`k_b A^*_{cp}[b](\Delta_w g)`$ if $`b \in \mathcal{I}\setminus H_u`$, $`k_b \ge 1`$ (Lemma 5.1(1), fact (c)). The $`k_a - 1`$ sibling remainders all have the same argument $`\Delta^L_p g`$, so each is bounded by the same $`A^*`$.
  - Tail end $`a = T_1(\kappa, h)`$: $`\mathop{\mathrm{ext}}_{\gamma \in \{\gamma_{\min}(\kappa, h), 1\}}\gamma(g[\mathrm{INF}, b] - g[a, b])`$; if $`\kappa \le K_T`$, also $`\mathop{\mathrm{ext}}_{\beta \in \{b_{\min}, 1\}}\beta(g[T_1(\kappa+1, h), b] - g[a, b])`$; and the quadratic term $`\frac{1}{D_0}\mathop{\mathrm{ext}}(g[\mathrm{INF}, b^+] - g[\mathrm{INF}, b],\ g[a, b^+] - g[a, b])`$ for the credit $`\iota_p\,\Delta^L_w g`$ (facts (d), (e)). At $`a = \mathrm{INF}`$ all these summands are $`0`$.
  - Tail end $`b = T_1(\kappa, 1)`$: the mirror image, with quadratic term $`\frac{1}{D_0}\mathop{\mathrm{ext}}(g[a^+, \mathrm{INF}] - g[a, \mathrm{INF}],\ g[a^+, b] - g[a, b])`$.

  The check is $`\mathrm{LHS}_g(a, b) \lessgtr \lambda g[a, b]`$ for every pair with at least one inner end.

**Path family $`P \lessgtr \lambda h_V`$.** Corners are $`(a, s, b) \in \mathcal{P} \times H_u \times \mathcal{W}`$. The middle $`s`$ is inner with $`D_s + 1 \le D_{hu} \le D_0`$, so $`s^+`$ is a table state. If $`s^+ \notin H_u`$ or $`s^{\mathrm{bud}} \notin H_u`$, the corresponding value is $`0`$ (exit, Lemma 2.3).
- Sources: $`(\mathrm{inv}(a) - x_{\mathrm{par}}[s])(g[s^+, b] - g[s, b]) + (\mathrm{inv}(b) - x_{\mathrm{ch}}[s])(g[a, s^+] - g[a, s])`$ (fact (a)).
- Top $`a`$ inner: rate in $`[L_a + \mathrm{inv}(s),\ L_a + \mathrm{inv}(s) + (h_a + k_a - 1)/2]`$, since $`s`$ is one of its $`h_a + k_a`$ core neighbours (fact (b)); budding exact. Top $`a`$ tail: Lemma 5.2 together with the quadratic term $`\frac{1}{D_0}\mathop{\mathrm{ext}}(h_V[\mathrm{INF}, s^+, b] - h_V[\mathrm{INF}, s, b],\ h_V[a, s^+, b] - h_V[a, s, b])`$ for $`\iota_p\,\Delta^L_u h_V`$.
- Bottom $`b`$ inner: rate in $`[L_b + \mathrm{inv}(s),\ L_b + \mathrm{inv}(s) + k_b/2]`$. Bottom $`b`$ tail: as for the top, with quadratic term $`\frac{1}{D_0}\mathop{\mathrm{ext}}(h_V[a, s^+, \mathrm{INF}] - h_V[a, s, \mathrm{INF}],\ \cdot)`$.
- Middle: rate in $`[r_0, r_0 + (k_s - 1)/2]`$ with $`r_0 = L_s + [a \in \mathcal{I}]\mathrm{inv}(a) + [b \in \mathcal{I}]\mathrm{inv}(b)`$, since $`s`$ has $`k_s - 1`$ other core neighbours. Budding is exact.
- When both $`a`$ and $`b`$ are tail, the middle's credit $`(t_p + t_w)/D_0`$ multiplies a jump that is bi-affine in $`(t_p, t_w)`$. Lemma 5.3 is applied to each summand separately. The verifier adds the two quadratic terms independently at each corner, which is the correct bound.

**Cherry family $`Ch \lessgtr \lambda h_C`$.** Corners are $`(s, b_1, b_2) \in C_u \times \mathcal{W} \times \mathcal{W}`$, unordered, including $`b_1 = b_2`$. Two distinct children in the same table state give the corner $`b_1 = b_2`$, and the two child ends are still treated as separate vertices.
- Sources: $`(\mathrm{inv}(b_1) - x_{cc}[s])(g[s^+, b_2] - g[s, b_2]) + (\mathrm{inv}(b_2) - x_{cc}[s])(g[s^+, b_1] - g[s, b_1])`$.
- Each child end $`w_i`$, inner: rate in $`[L_i + \mathrm{inv}(s),\ L_i + \mathrm{inv}(s) + k_i/2]`$. Tail: Lemma 5.2 and the quadratic term in $`\iota_{w_i}\Delta^L_u h_C`$.
- Apex: rate in $`[r_0, r_0 + (h_s + k_s - 2)/2]`$ with $`r_0 = L_s + \sum_i [b_i \in \mathcal{I}]\mathrm{inv}(b_i)`$. $`s`$ may be the root ($`h_s = 0`$).

Applying Lemma 5.4 to each family gives (5.0) in both directions, given that every finite inequality holds. That is what the verifier checks, with 0 failures (§8).

---

## 6. Comparison of $`\Phi`$ with $`C`$ in Both Directions

Assign each edge term and each path term to its bottom vertex, and each cherry term to its apex. For a core vertex $`u`$, let $`W(u)`$ be $`f(\pi_f(s_u))`$ plus the terms assigned to $`u`$. Then $`\Phi(T) = \sum_{u \in K} W(u)`$. Vertex $`u`$ receives at most one edge term (none for the root), at most one path term, and $`\binom{k_u}{2}[s_u \in C_u] \le \binom{K_{cu}}{2}`$ cherry terms.

For $`w \in \mathcal{W}`$ and $`s \in C_u`$ put (UPPER with $`\mathop{\mathrm{ext}} = \min`$, LOWER with $`\max`$; the $`0`$ covers an absent term)
```math
g_{\mathrm{ext}}(w) := \mathop{\mathrm{ext}}\big(0, \mathop{\mathrm{ext}}_{a \in \mathcal{P}} g[a, w]\big), \quad h_{V,\mathrm{ext}}(w) := \mathop{\mathrm{ext}}\big(0, \mathop{\mathrm{ext}}_{a \in \mathcal{P},\, s' \in H_u} h_V[a, s', w]\big), \quad h_{C,\mathrm{ext}}(s) := \mathop{\mathrm{ext}}\big(0, \mathop{\mathrm{ext}}_{b_1, b_2} h_C[s, b_1, b_2]\big).
```
For a core vertex $`u`$ define its **candidate endpoints** $`E(u) \subset \mathcal{G}`$: $`E(u) := \{s_u\}`$ if $`D_u < D_0`$, and $`E(u) := \{\mathrm{INF}, T_1(\kappa_u, h_u)\}`$ if $`D_u \ge D_0`$. Then $`\mathrm{supp}\,\mu_u \subseteq E(u)`$. The inclusion can be strict: at $`D_u = D_0`$ we have $`t_u = 1`$, so $`\mathrm{INF}`$ has weight $`0`$, but it stays in $`E(u)`$. An evaluated term is a convex combination of table entries whose parent (resp. middle) argument lies in $`\mathcal{P}`$ (resp. $`H_u`$) and whose bottom argument lies in $`\mathrm{supp}\,\mu_u \subseteq E(u)`$. An average over a subset of $`E(u)`$ is at least the minimum over all of $`E(u)`$. So in the UPPER direction
```math
W(u) \ge f(\pi_f(s_u)) + h_u\big(\mathop{\mathrm{ext}}_{b \in E(u)} g_{\mathrm{ext}}(b) + \mathop{\mathrm{ext}}_{b \in E(u)} h_{V,\mathrm{ext}}(b)\big) + [s_u \in C_u]\tbinom{k_u}{2}h_{C,\mathrm{ext}}(s_u) =: \underline W(u),
```
and in the LOWER direction the same holds with $`\le`$ and $`\mathop{\mathrm{ext}} = \max`$; call that right side $`\overline W(u)`$. The verifier takes the extremum over both candidate endpoints for every tail state, including the corner of weight $`0`$. This gives a bound that is never sharper than the one over $`\mathrm{supp}\,\mu_u`$, and it is still valid. The right side depends only on $`\pi_f(s_u)`$. For $`L < L_F`$ and $`k \le K_F`$ this is the state itself. For $`L \ge L_F`$ and $`k \le K_F`$ the state is tail because $`D \ge L_F \ge D_0`$, and $`E(u) = \{\mathrm{INF}, T_1(\min(k, K_T+1), h)\}`$ does not depend on $`L`$. For $`\mathrm{VSTAR}_h`$ the state is tail because $`D \ge K_F + 1 \ge D_0`$, $`\kappa = K_T + 1`$, and $`s_u \notin C_u`$. So $`\underline W`$ and $`\overline W`$ range over finitely many values. Let $`c_0`$ be the minimum of $`\underline W`$ (UPPER) and $`C_0`$ the maximum of $`\overline W`$ (LOWER), over all values of $`\pi_f(s_u)`$. The verifier finds
```math
c_0 = \tfrac{2{,}989{,}767{,}809}{10^{10}} > 0, \qquad C_0 = \tfrac{17{,}527{,}873{,}930}{10^{10}} > 0,
```
and hence $`\Phi_{\mathrm{up}}(T) \ge c_0\,C(T)`$ and $`\Phi_{\mathrm{lo}}(T) \le C_0\,C(T)`$ for every tree with $`n \ge 3`$.

---

## 7. Completion of the Proof of Theorem 1

$`T_3`$ has a single core vertex, the root in state $`(2, 0, 0)`$, and no core edge. So $`\Phi(T_3) = f(2, 0, 0)`$. This equals $`14{,}606{,}753{,}031/10^{10}`$ for UPPER and $`14{,}608{,}704{,}616/10^{10}`$ for LOWER. Both are positive.

**UPPER.** By §5, $`(\mathcal{L}\Phi_{\mathrm{up}})(T) \le \lambda\Phi_{\mathrm{up}}(T)`$ for every $`T`$. Hence $`\mathbb{E}\Phi_{\mathrm{up}}(T_{n+1}) = \mathbb{E}\Phi_{\mathrm{up}}(T_n) + \frac{1}{n}\mathbb{E}(\mathcal{L}\Phi_{\mathrm{up}})(T_n) \le (1 + \frac{\lambda}{n})\mathbb{E}\Phi_{\mathrm{up}}(T_n)`$. Induction from $`n = 3`$ and §6 give
```math
\mathbb{E}[C_n] \le \frac{\mathbb{E}\Phi_{\mathrm{up}}(T_n)}{c_0} \le \frac{\Phi_{\mathrm{up}}(T_3)}{c_0}\prod_{j=3}^{n-1}\Big(1 + \frac{1133}{2000\,j}\Big), \qquad \frac{\Phi_{\mathrm{up}}(T_3)}{c_0} = \frac{14{,}606{,}753{,}031}{2{,}989{,}767{,}809}.
```
**LOWER.** By §5, $`(\mathcal{L}\Phi_{\mathrm{lo}})(T) \ge \lambda\Phi_{\mathrm{lo}}(T)`$. Hence $`\mathbb{E}\Phi_{\mathrm{lo}}(T_{n+1}) \ge (1 + \frac{\lambda}{n})\mathbb{E}\Phi_{\mathrm{lo}}(T_n)`$, and
```math
\mathbb{E}[C_n] \ge \frac{\mathbb{E}\Phi_{\mathrm{lo}}(T_n)}{C_0} \ge \frac{\Phi_{\mathrm{lo}}(T_3)}{C_0}\prod_{j=3}^{n-1}\Big(1 + \frac{113}{200\,j}\Big), \qquad \frac{\Phi_{\mathrm{lo}}(T_3)}{C_0} = \frac{14{,}608{,}704{,}616}{17{,}527{,}873{,}930} = \frac{7{,}304{,}352{,}308}{8{,}763{,}936{,}965}.
```
Neither induction step needs $`\Phi \ge 0`$: each step only takes expectations of a pointwise linear inequality. The exponent statements follow from $`\prod_{j=3}^{n-1}(1 + \alpha/j) = \frac{2\,\Gamma(n+\alpha)}{\Gamma(3+\alpha)\Gamma(n)} \sim \frac{2}{\Gamma(3+\alpha)}n^{\alpha}`$. $`\square`$

---

## 8. Verification Record

The certificates are produced by `verify/motif-bracket/cert_build.c` and checked by `verify/motif-bracket/cert_verify_gmp.c`. The proof depends only on the checker; the builder is untrusted. The binary format is given in §9. The correspondence between this text and the checker, line by line, is in `bracket-motifs-code-map.md`. That file also separates what the checker establishes for the two shipped headers from what it would establish for an arbitrary binary, and records the integer headroom.

Run on 2026-11-23 (GCC 13.3.0, GMP 6.3.0, x86_64 Linux):

The lossless reference copies can be decompressed with `xz -dc` and passed directly to the GMP verifier. Upper reference array · Lower reference array

| Item | UPPER `upper_5665.bin` | LOWER `lower_5650.bin` |
|---|---|---|
| sha256 of rebuilt file | `3a4aef888a5b567a8ada05f9546f98c979d361da9f9b33ba6e794eb91fb30abe` | `9e966fba9aee6eae3da659b693d3bb9c481c1cea66635b5218d0887f182a60d9` |
| f / g / hV / hC checks | 484,138 / 160,545 / 13,519,800 / 3,164,625 | same |
| failures | 0 | 0 |
| comparison constant | $`c_0 = 2{,}989{,}767{,}809/10^{10}`$ | $`C_0 = 17{,}527{,}873{,}930/10^{10}`$ |
| prefactor | $`14{,}606{,}753{,}031/2{,}989{,}767{,}809 \approx 4.885581`$ | $`7{,}304{,}352{,}308/8{,}763{,}936{,}965 \approx 0.833456`$ |
| brute-force trees in checker | 37 of 37 pass | 37 of 37 pass |
| UBSan build (`-fsanitize=undefined -fno-sanitize-recover=all`) | 0 runtime errors, same output | 0 runtime errors, same output |
| exponent override (`1132/2000`, resp. `114/200`) | rejected, exit 2 | rejected, exit 2 |
| single-entry perturbation of `f`, `g`, `hV`, `hC` | each rejected, exit 2 | each rejected, exit 2 |

The 37 checker trees are $`T_3`$, $`K_{1,3}`$, $`P_4`$, four constructed trees (degrees $`D_0 - 1`$ and $`D_0 + 5`$ next to a hub, a chain of $`L = 0`$ vertices, a root with $`k = K_F + 2`$, a root with $`L = L_F + 2`$), and 30 sampled friend trees with $`5 \le n \le 45`$. `verify/run_all.sh` applies the perturbation controls to the UPPER file only. The LOWER perturbation rows above come from this audit's own runs.

**Identity (I) and (5.0) on trees (audit script, independent of the checker).** The script reads the regenerated binaries with its own loader. On each tree it evaluates (i) $`(\mathcal{L}\Phi)(T)`$ by brute force, (ii) every local expression of §4.3 in exact `Fraction` arithmetic, (iii) inequality (5.0) for every term, and (iv) the comparison of §6.
- Random set: 172 trees with up to $`n = 254`$ vertices per certificate, from four seeds. Each seed gives $`T_3`$, $`K_{1,3}`$, $`P_4`$, 20 sampled friend trees with $`n < 160`$, and 20 constructed trees with $`n \le 260`$. Together these trees contain 182 tail vertices, 266 edges with only the parent in the tail, 64 with only the child in the tail, 46 tail–tail edges, 18 inner vertices with $`D = D_0 - 1`$, 205 path middles in $`H_u`$ (85 with a tail end), 49 cherry apexes in $`C_u`$ (16 with a tail end), all three FY splits, and clamped $`\kappa = K_T + 1`$.
- Boundary set: five trees covering $`L = L_F + 1`$ at the root, $`L = L_F + 2`$ below the root, $`L = L_F`$ with $`k = K_F`$ (budding into $`\mathrm{VSTAR}_1`$), $`\mathrm{VSTAR}_0`$ with $`k = K_F + 2`$, and $`\mathrm{VSTAR}_1`$ with $`k = K_F + 1`$. They have up to $`n = 1104`$ vertices.
- Result, both certificates: identity failures 0, failures of (5.0) 0, comparison failures 0 (4,460 local expressions in the random set and 1,838 in the boundary set, counted over both certificates).

This confirms the term inventory, the multiplicities and the newborn bookkeeping of §4 on concrete trees. It does not replace the finite check of §5, which is what proves (5.0) for all trees.

---

## 9. Binary Certificate Format

The digests in Section 8 identify the byte streams used for the reported constants. Compiler floating-point choices during the search can change the proposed array. On a different build, check the resulting array and its seed and core-comparison values in exact arithmetic.

Both certificates (`upper_5665.bin` and `lower_5650.bin`) are generated by `verify/motif-bracket/cert_build.c` in the binary format defined in `verify/motif-bracket/cert_format.h` (little-endian x86_64):
1. `CertHeader` (`120` bytes, `struct.unpack("<Qiiiiqqiiiiiiiiii9ii", hdr)`): `magic = 0x494E564433563031ULL`, `mode` (`1`=UPPER, `2`=LOWER), `P, Q`, `reserved0`, `S = 10^10`, `S_r = 10^6`, `LF=800, KF=150, D0=24, KT=12, DHU=16, KHU=8, DCU=12, KCU=6`, `Q_K=240`, `nknots=9`, `knot_num = {0, 10, 15, 25, 40, 60, 80, 100, 120}`, `reserved1`.
   - All integer weights $`W \in \text{int64\_t}`$ in `f, g, hV, hC, k_pc, k_cp, k_cc` represent exact rationals $`W/S`$.
   - All integer centering constants $`R \in \text{int64\_t}`$ in `rpar, rch, xpar, xch, xcc` represent exact rationals $`R/S_r`$.
   - Knot abscissae are $`\xi_j = \text{knot\_num}[j] / Q_K`$.
2. Exact deterministic ordering of states:
   - `INNER_N = 572`: nested loop `h = 1, 0`; `D = 2..D0-1`; `k = 0..D-h`; `L = D-h-k`.
   - `INF_ID = 572`.
   - `T1(kap, h)`: `573 + (h ? 0 : KT+2) + kap` for `kap = 0..KT+1`, `h = 1, 0` (`NG = 601`).
   - `pars` (`npars = 555`): all `u in 0..NG-1` in increasing order with `u == INF_ID` or `gK[u] >= 1`.
   - `chs` (`nchs = 290`): all `u in 0..NG-1` in increasing order with `u == INF_ID` or `gH[u] == 1`.
   - `hu_list` (`NHU_ACT = 84`): `k = 1..KHU`, `L = 0..DHU-1` with `L+k+1 < DHU`.
   - `cu_list` (`NCU_ACT = 75`): `h = 1, 0`, `k = 2..KCU`, `L = 0..DCU-1` with `L+k+h < DCU`.
   - `vidx(L, k, h)`: `(h ? 0 : (KF+1)*(LF+1)) + k*(LF+1) + L`; `VSTAR1 = 2*(KF+1)*(LF+1)`, `VSTAR0 = VSTAR1 + 1` (`NV = 241904`).
3. Payload arrays (exact `int64_t` in order):
   - `rpar[INNER_N]`, `rch[INNER_N]`, `xpar[INNER_N]`, `xch[INNER_N]`, `xcc[INNER_N]`
   - `k_pc[INNER_N][9]`, `k_cp[INNER_N][9]`, `k_cc[INNER_N][9]`
   - `f[NV]`
   - `g[NG * NG]` (row-major `p * NG + w`)
   - `hV[npars * NHU_ACT * nchs]` (row-major `(pi * NHU_ACT + ua) * nchs + wi`)
   - `hC[NCU_ACT * nchs * (nchs + 1) / 2]` (row-major `ca * pair_n + (i1*nchs - i1*(i1-1)/2 + (i2-i1))` for `0 <= i1 <= i2 < nchs`).
