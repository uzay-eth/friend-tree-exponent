> **Public checkpoint export note (Part 1 of 2 — Sections 1–3):** Edited export of `order-neon/PROOF.md` (frozen original SHA-256 `e85e9cebfbde5439fad9246da2fec2c6eb7f014b09248e960c20d016411ce191`). To stay within GitHub’s per-page browser math-rendering limit, the exported proof is split at the Section 3 / Section 4 boundary into [`PROOF.md`](PROOF.md) (Sections 1–3) and [`PROOF-PART-2.md`](PROOF-PART-2.md) (Sections 4–10) with navigation links added; no mathematical content, formula, label, erratum, or proof step was added or altered (all audit verdicts apply to the single frozen original source SHA-256). Formatting changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in these reformatted exports.** Read alongside [`REVIEW-Q3.md`](REVIEW-Q3.md) (which records two local wording corrections and the 11,307 accepted vs. 12,105 attempted finite count distinction) and [`REVIEW-IMPLEMENTATION.md`](REVIEW-IMPLEMENTATION.md).

# Complete Self-Contained Proof of Birth-Order Stochastic Dominance (Addario-Berry et al. §7 Question 3, Second Half)

**Author/lead:** Gemini 4 Neon (`order-neon`)  
**External proof dependencies:** **None** (100% self-contained; relies on zero Tier B / Tier C round-2 results and zero floating-point inequalities).

> **Navigation (Part 1 of 2 — Sections 1–3):** [§1 Model & Main Theorems](#1-exact-model-and-main-theorems) · [§2 Signed Identity](#2-exact-pre-placement-signed-identity) · [§3 Endpoint Laws](#3-exact-endpoint-comparison-laws) | **Continue to Part 2 ([`PROOF-PART-2.md`](PROOF-PART-2.md)):** [§4 All-m Bounds](PROOF-PART-2.md#4-exact-discrete-time-continuation-bounds-for-all-m-ge-m_0--200-cert_all_mcpp) · [§5 Finite-m Induction](PROOF-PART-2.md#5-exact-integer-n-backward-induction-for-4-le-m-le-199-cert_small_mcpp) · [§6 Tree Ledger](PROOF-PART-2.md#6-dimensionless-normalized-all-tree-signed-ledger) · [§7 High-Tau Tail](PROOF-PART-2.md#7-self-contained-analytic-high-tau-tail-tau-ge-6-all-b-ge-3) · [§8 Founder Seam](PROOF-PART-2.md#8-self-contained-founder-seam-b--2-z_2-succ_mathrmst-z_3) · [§9 Int128 Arithmetic](PROOF-PART-2.md#9-exact-outward-rounded-__int128-arithmetic-guarantees) · [§10 Replay & Evidence](PROOF-PART-2.md#10-replay-commands-and-exact-verification-evidence)

---

## 1. Exact Model and Main Theorems

### 1.1 The Random Friend Tree Process
Start from the two-vertex tree $`T_2`$ consisting of a single edge $`\{1, 2\}`$. For each integer $`n \ge 2`$, conditionally on $`T_n`$ on vertex set $`\{1, \dots, n\}`$ with degree function $`d_n(u) := \deg_{T_n}(u)`$, form $`T_{n+1}`$ by adding vertex $`n+1`$ attached to a random vertex $`W_n \in V(T_n)`$ sampled by picking $`U_n \in V(T_n)`$ uniformly at random and then picking $`W_n`$ uniformly among the neighbors of $`U_n`$ in $`T_n`$. Equivalently,

```math
\mathbb{P}(n+1 \sim u \mid T_n) = \frac{w_n(u)}{n}, \qquad w_n(u) := \sum_{v \sim u} \frac{1}{d_n(v)}.
```

Since $`\sum_{u \in V(T_n)} w_n(u) = \sum_{v \in V(T_n)} d_n(v)\cdot \frac{1}{d_n(v)} = n`$, this defines a probability distribution on $`V(T_n)`$.

For each vertex label $`u \ge 1`$, the sequence $`(d_N(u)/N)_{N \ge u}`$ is a nonnegative supermartingale (its conditional expectation at step $`N+1`$ is $`(d_N(u) + w_N(u)/N)/(N+1) \le d_N(u)/N`$ because $`w_N(u) \le d_N(u)`$), so the limiting degree share

```math
Z_u := \lim_{N \to \infty} \frac{d_N(u)}{N} \in [0, 1]
```

exists simultaneously for all $`u \ge 1`$ almost surely (Addario-Berry et al., *Random Structures & Algorithms* 69(1), e70082 (2026); arXiv:2403.20185v1 §3).

### 1.2 Main Results

**Theorem 1 (Universal Pre-Placement Conditional Stochastic Dominance).**  
*For every integer $`n \ge 2`$, every finite tree $`T_n`$ on $`n`$ vertices, and every threshold $`x \in [0, 1)`$,*

```math
\mathbb{P}(Z_{n+1} > x \mid T_n) > \mathbb{P}(Z_{n+2} > x \mid T_n).
```

**Theorem 2 (Complete Resolution of Addario-Berry et al. §7 Question 3).**  
*The limiting degree shares $`(Z_u)_{u \ge 1}`$ satisfy exact founder symmetry and strict first-order stochastic dominance from birth label $`2`$ onward:*

```math
Z_1 \stackrel{d}{=} Z_2 \succ_{\mathrm{st}} Z_3 \succ_{\mathrm{st}} Z_4 \succ_{\mathrm{st}} \cdots,
```

*that is, $`\mathbb{P}(Z_1 > x) = \mathbb{P}(Z_2 > x)`$ and $`\mathbb{P}(Z_u > x) > \mathbb{P}(Z_v > x)`$ for every $`2 \le u < v`$ and every $`x \in [0, 1)`$ (with trivial equality at $`x \ge 1`$, where both sides are $`0`$, and at $`x < 0`$, where both sides are $`1`$).*

*Proof of Theorem 2 from Theorem 1 and the Founder Seam (Section 8).*  
The graph automorphism swapping $`1 \leftrightarrow 2`$ in $`T_2 = \{1, 2\}`$ preserves the transition law, so $`Z_1 \stackrel{d}{=} Z_2`$. Section 8 proves $`\mathbb{P}(Z_2 > x) > \mathbb{P}(Z_3 > x)`$ for every $`x \in [0, 1)`$. For every $`b = n + 1 \ge 3`$ ($`n \ge 2`$) and every $`x \in [0, 1)`$, taking expectation over $`T_n`$ in Theorem 1 gives $`\mathbb{P}(Z_b > x) > \mathbb{P}(Z_{b+1} > x)`$. Transitivity gives $`\mathbb{P}(Z_u > x) > \mathbb{P}(Z_v > x)`$ for all $`2 \le u < v`$ and $`x \in [0, 1)`$. $`\square`$

---

## 2. Exact Pre-Placement Signed Identity

Fix $`n \ge 2`$ and an arbitrary finite tree $`T_n`$. Write $`b := n+1`$, $`m := n+2 \ge 4`$, $`d(u) := d_n(u)`$, $`w(u) := w_n(u) = \sum_{v \sim u} 1/d(v)`$, and define the **edge-drop weight**

```math
s(u) := \frac{w(u)}{d(u)\big(d(u)+1\big)} \qquad (u \in V(T_n)).
```

Let $`T_m(b \sim u, m \sim v)`$ be the tree formed from $`T_n`$ by attaching leaf $`b`$ to $`u`$ and leaf $`m`$ to $`v`$.

**Lemma 2.1 (Exact Pre-Placement Decomposition).**  
*For every finite tree $`T_n`$ ($`n \ge 2`$) and every $`x \in [0, 1)`$,*

```math
n(n+1)\Big(\mathbb{P}(Z_b > x \mid T_n) - \mathbb{P}(Z_m > x \mid T_n)\Big) = \sum_{u \in V(T_n)} d(u)s(u)\, \Delta_{\mathrm{ser}}(u) - \sum_{\{u, v\} \in E(T_n)} \Big[(s(u) - s(v))_+ \Delta_{\mathrm{par}}(v, u) + (s(v) - s(u))_+ \Delta_{\mathrm{par}}(u, v)\Big],
```

*where:*
1. *$`\Delta_{\mathrm{ser}}(u) := \mathbb{P}(Z_b > x \mid T_m(b \sim u, m \sim b)) - \mathbb{P}(Z_m > x \mid T_m(b \sim u, m \sim b))`$ is the serial gap at $`u`$: $`b`$ has degree $`2`$, $`1`$ leaf neighbor $`m`$, and $`1`$ host neighbor $`u`$ of degree $`d(u)+1`$, while $`m`$ is a leaf on host $`b`$ of degree $`2`$ with $`0`$ other leaf neighbors.*
2. *$`\Delta_{\mathrm{par}}(v, u) := \mathbb{P}(Z_m > x \mid T_m(b \sim v, m \sim u)) - \mathbb{P}(Z_b > x \mid T_m(b \sim v, m \sim u))`$ is the adverse parallel gap on $`\{b \sim v, m \sim u\}`$: $`m`$ is a leaf on host $`u`$ of degree $`d(u)+1`$, and $`b`$ is a leaf on host $`v`$ of degree $`d(v)+1`$.*

*Proof.* Condition on $`b \sim u`$ (probability $`w(u)/n`$). Vertex $`m`$ then attaches to:
- **$`u`$ (twin):** $`b`$ and $`m`$ are symmetric leaves on $`u`$ in $`T_m`$, so their marginal laws match and contribute $`0`$.
- **$`b`$ (serial):** $`U_{n+1}`$ must be $`u`$ (the unique neighbor of $`b`$) and $`W_{n+1} = b`$ (probability $`1/((n+1)(d(u)+1))`$). Multiplying by $`w(u)/n`$ yields $`\frac{w(u)}{n(n+1)(d(u)+1)} = \frac{d(u)s(u)}{n(n+1)}`$.
- **$`v \neq u`$ not adjacent to $`u`$:** attaching $`b`$ to $`u`$ does not change any neighbor degree of $`v`$, so $`\mathbb{P}(b \sim u, m \sim v \mid T_n) = \frac{w(u)w(v)}{n(n+1)} = \mathbb{P}(b \sim v, m \sim u \mid T_n)`$. Swapping $`b \leftrightarrow m`$ is a tree isomorphism between $`T_m(b \sim u, m \sim v)`$ and $`T_m(b \sim v, m \sim u)`$, so every non-adjacent pair cancels identically.
- **$`v \sim u`$ (an edge $`\{u, v\} \in E(T_n)`$):** attaching $`b`$ to $`u`$ increases $`d(u)`$ to $`d(u)+1`$, reducing $`w(v)`$ by $`\frac{1}{d(u)} - \frac{1}{d(u)+1} = \frac{1}{d(u)(d(u)+1)}`$. Hence

```math
\mathbb{P}(b \sim u, m \sim v \mid T_n) = \frac{w(u)w(v) - s(u)}{n(n+1)}, \qquad \mathbb{P}(b \sim v, m \sim u \mid T_n) = \frac{w(u)w(v) - s(v)}{n(n+1)}.
```

Pairing the two orientations cancels their common mass by the $`b \leftrightarrow m`$ isomorphism, leaving net probability $`\frac{(s(u) - s(v))_+}{n(n+1)}`$ on $`\{b \sim v, m \sim u\}`$ and $`\frac{(s(v) - s(u))_+}{n(n+1)}`$ on $`\{b \sim u, m \sim v\}`$. $`\square`$

**Leaf-stock invariant at completion $`T_m`$.** Let $`l(u)`$ be the number of leaf neighbors of $`u`$ in $`T_n`$. In the adverse configuration $`\{b \sim v, m \sim u\}`$ along edge $`\{u, v\} \in E(T_n)`$:
- If $`d(v) \ge 2`$, then $`v`$ remains a non-leaf neighbor of $`u`$ in $`T_m`$, so host $`u`$ has degree $`k = d(u)+1`$, at least one non-leaf neighbor ($`v`$), and $`j = l(u) \le d(u) - 1 = k - 2`$ other leaf neighbors besides $`m`$.
- If $`d(v) = 1`$, then $`b \sim v`$ turns $`v`$ into a degree-$`2`$ non-leaf neighbor of $`u`$ in $`T_m`$, so host $`u`$ has degree $`k = d(u)+1`$, at least one non-leaf neighbor ($`v`$), and $`j = l(u) - 1 \le d(u) - 1 = k - 2`$ other leaf neighbors besides $`m`$.
- In the serial tip $`\{b \sim u, m \sim b\}`$, host $`b`$ has degree $`k = 2`$, one non-leaf neighbor ($`u`$), and $`j = 0 = k - 2`$ other leaf neighbors besides $`m`$.

Thus **every adverse pending leaf in Lemma 2.1 satisfies the permanent invariant $`k \ge j + 2`$** (the host has a permanent non-leaf neighbor, since vertex degrees never decrease).

---

## 3. Exact Endpoint Comparison Laws

For real $`N > r > 0`$ and $`x \in [0, 1)`$, let

```math
B_{N, r}(x) := \mathbb{P}\big(\mathrm{Beta}(r, N-r) > x\big) = \frac{\Gamma(N)}{\Gamma(r)\Gamma(N-r)} \int_x^1 u^{r-1}(1-u)^{N-r-1}\,du,
```

with $`B_{N, r}(x) := 0`$ for $`r \le 0`$ and $`B_{N, N}(x) := 1`$. For integer $`r \in \{1, \dots, N-1\}`$, $`B_{N, r}(x) = \mathbb{P}(\mathrm{Bin}(N-1, x) \le r-1) = (1-x)^{N-1}\sum_{i=0}^{r-1}\binom{N-1}{i}\left(\frac{x}{1-x}\right)^i`$. Note that $`B_{N, r}(x)`$ is strictly decreasing in $`x \in [0, 1)`$, strictly decreasing in real $`N > r`$, and strictly increasing in $`r \in (0, N)`$.

**Lemma 3.1 (Poisson-Beta Submartingale Lower Bound).**  
*At any size $`N \ge 2`$, for any vertex $`v \in V(T_N)`$ with at least $`j \in \{1, \dots, N-1\}`$ leaf neighbors and any $`x \in [0, 1)`$,*

```math
\mathbb{P}(Z_v > x \mid T_N) \ge H_{N, j}(x) := e^{-1}\sum_{k=0}^{j-1}\frac{1}{k!} B_{N, j-k}(x).
```

*Proof.* Set $`H_{N, 0} \equiv 0`$. From the exact Pólya-urn identity $`N B_{N, r} = r B_{N+1, r+1} + (N-r) B_{N+1, r}`$ and $`k(e^{-1}/k!) = e^{-1}/(k-1)!`$, telescoping $`\Delta B(r) = B_{N+1, r} - B_{N+1, r-1}`$ gives the discrete harmonic identity for all $`1 \le j < N`$:

```math
N H_{N, j}(x) = j H_{N+1, j+1}(x) + (N - j - 1) H_{N+1, j}(x) + H_{N+1, j-1}(x).
```

Designate $`j`$ current leaf neighbors of $`v`$ as tracked leaves, add each subsequent leaf attached to $`v`$ to the tracked set, and remove a tracked leaf when a new vertex attaches to it. At tree size $`s`$ with tracked stock $`i \ge 1`$, the probability of attaching to $`v`$ is at least $`i/s`$ (each tracked leaf has degree $`1`$), the disjoint probability of attaching to a tracked leaf is $`i/(s d_s(v)) \le 1/s`$, and $`H_{s+1, i}`$ increases in $`i`$. Thus $`(H_{s, I_s}(u))_{s \ge N}`$ is a bounded submartingale for every $`u \in (x, 1)`$. On $`\{Z_v < u\}`$, $`I_s/s \le d_s(v)/s \to Z_v < u`$, so $`H_{s, I_s}(u) \le B_{s, I_s}(u) = \mathbb{P}(\mathrm{Bin}(s-1, u) < I_s) \to 0`$. Reverse Fatou gives $`H_{N, j}(u) \le \mathbb{P}(Z_v \ge u \mid T_N)`$, and letting $`u \downarrow x`$ by continuity of $`H_{N, j}`$ completes the proof. $`\square`$

**Lemma 3.2 (Locked Beta Upper Bound).**  
*At any size $`N \ge 2`$, for any vertex $`v \in V(T_N)`$ with current attachment weight $`r_0 = w_N(v) < N`$ and any $`x \in [0, 1)`$,*

```math
\mathbb{P}(Z_v > x \mid T_N) \le B_{N, r_0}(x).
```

*Proof.* Any attachment to a neighbor $`u \sim v`$ can only increase $`d_s(u)`$, hence can only decrease $`1/d_s(u)`$, while each attachment to $`v`$ adds a new leaf neighbor of weight $`1`$. Thus $`w_s(v) \le r_0 + (d_s(v) - d_N(v))`$ for all $`s \ge N`$, so $`(d_s(v) - d_N(v))_{s \ge N}`$ is stochastically dominated by a Pólya urn starting with $`(r_0, N-r_0)`$ balls, whose limiting fraction is $`\mathrm{Beta}(r_0, N-r_0)`$. $`\square`$

**Lemma 3.3 (Static Time-Independent Hubhood Upper Bound).**  
*Fix $`K_{\max} \ge 2`$. Define $`U_{\mathrm{hub}}(K_{\max}+1, j) := \frac{1}{j+1}`$ for $`0 \le j \le K_{\max}-1`$, $`U_{\mathrm{hub}}(K_{\max}+1, K_{\max}) := \min(1, \frac{2}{K_{\max}+1})`$, and recursively for $`k = K_{\max}`$ down to $`2`$ and $`0 \le j < k`$:*

```math
U_{\mathrm{hub}}(k, j) := \min\left(c(k, j),\; \max_{g \in \{j+1,\, (k+j+1)/2\}} \frac{\frac{1}{k}\cdot 1 + \frac{j}{k} U_{\mathrm{hub}}(k, j-1) + g U_{\mathrm{hub}}(k+1, j+1)}{\frac{1+j}{k} + g}\right),
```

*where $`U_{\mathrm{hub}}(k, -1) := 0`$ and $`c(k, j) := \frac{1}{j+1}`$ if $`j \le k - 2`$, while $`c(k, k-1) := \min\left(1, \frac{2}{k}\right)`$. Then for every $`N \ge 4`$, every $`x \in [0, 1)`$, and every leaf $`v \in V(T_N)`$ attached to a host $`w`$ of degree $`k \ge 2`$ with $`j \in \{0, \dots, k-1\}`$ other leaf neighbors,*

```math
\mathbb{P}(Z_v > x \mid T_N) \le \mathbb{P}(Z_v > 0 \mid T_N) \le \mathbb{P}(v \text{ ever receives a child} \mid T_N) \le U_{\mathrm{hub}}(k, j).
```

*Proof.* A leaf $`v`$ that never receives a child has $`d_N(v) \equiv 1`$, hence $`Z_v = 0`$. Track the state $`(k, j)`$ of host $`w`$ (degree $`k`$, $`j`$ other leaf neighbors besides $`v`$) until $`v`$ receives its first child (payoff $`1`$). Before $`v`$ receives a child, the only events changing $`(k, j)`$ or stopping are:
- attachment to $`v`$ (rate $`1/k`$, stops with payoff $`1`$);
- attachment to one of the $`j`$ other leaves of $`w`$ (rate $`j/k`$, moves to $`(k, j-1)`$);
- attachment to $`w`$ (rate $`g = (j+1) + \sum_{\text{non-leaf } z \sim w} 1/d(z) \in [j+1, (k+j+1)/2]`$, moves to $`(k+1, j+1)`$).

First, verify the analytic supermartingale caps $`c(k, j)`$:
1. **Permanent non-leaf invariant ($`j \le k - 2`$):** Since degrees never decrease, $`k \ge j + 2`$ persists on every transition. Test $`f_1(j) := \frac{1}{j+1}`$ with payoff $`1`$ at promotion. For $`j = 0`$, $`f_1(0) = 1`$ bounds all probabilities. For $`j \ge 1`$, the generator drift is

```math
g\left(\frac{1}{j+2} - \frac{1}{j+1}\right) + \frac{j}{k}\left(\frac{1}{j} - \frac{1}{j+1}\right) + \frac{1}{k}\left(1 - \frac{1}{j+1}\right) \le -\frac{1}{j+2} + \frac{1}{k(j+1)} + \frac{j}{k(j+1)} = -\frac{1}{j+2} + \frac{1}{k} \le 0
```

because $`g \ge j + 1`$ and $`k \ge j + 2`$. By optional stopping (the marked leaf $`v`$ itself triggers an event at rate $`1/k > 0`$, or host $`w`$ grows at rate $`\ge 1`$, so the chain exits every bounded region a.s.), $`\mathbb{P}(v \text{ ever promoted}) \le \frac{1}{j+1}`$.
2. **General case including $`j = k - 1`$ (no non-leaf neighbor):** Test $`f_2(j) := \min\left(1, \frac{2}{j+1}\right)`$. For $`j \in \{0, 1\}`$, $`f_2(j) = 1`$. For $`j \ge 2`$, using $`g \ge j + 1`$ and $`k \ge j + 1`$, the drift is

```math
g\left(\frac{2}{j+2} - \frac{2}{j+1}\right) + \frac{j}{k}\left(\frac{2}{j} - \frac{2}{j+1}\right) + \frac{1}{k}\left(1 - \frac{2}{j+1}\right) \le -\frac{2}{j+2} + \frac{j + 1}{k(j+1)} = -\frac{2}{j+2} + \frac{1}{k} \le -\frac{2}{j+2} + \frac{1}{j+1} = -\frac{j}{(j+1)(j+2)} < 0.
```

Hence $`\mathbb{P}(v \text{ ever promoted}) \le \min(1, \frac{2}{j+1})`$ unconditionally for all $`(k, j)`$.

Now backward induction on $`(K_{\max}+1 - k, j)`$ proves $`\mathbb{P}(v \text{ ever promoted} \mid k, j) \le U_{\mathrm{hub}}(k, j)`$: the base $`k = K_{\max}+1`$ is $`c(K_{\max}+1, j)`$, and for $`k \le K_{\max}`$, conditioning on the first jump of $`\{v, \text{other } j \text{ leaves}, w\}`$ gives an affine-in-$`g`$ mixture whose maximum over $`g \in [j+1, (k+j+1)/2]`$ is achieved at the two endpoints, and taking the minimum with the independently valid bound $`c(k, j)`$ preserves validity. Note that states with $`j \le k - 2`$ never transition to $`j' = k' - 1`$. $`\square`$


---

> **End of Part 1 (Sections 1–3).** **Continue to Part 2 (Sections 4–10):** [`PROOF-PART-2.md`](PROOF-PART-2.md) ([§4 `cert_all_m.cpp`](PROOF-PART-2.md#4-exact-discrete-time-continuation-bounds-for-all-m-ge-m_0--200-cert_all_mcpp) · [§5 `cert_small_m.cpp`](PROOF-PART-2.md#5-exact-integer-n-backward-induction-for-4-le-m-le-199-cert_small_mcpp) · [§6 Tree Ledger](PROOF-PART-2.md#6-dimensionless-normalized-all-tree-signed-ledger) · [§7 High-Tau Tail](PROOF-PART-2.md#7-self-contained-analytic-high-tau-tail-tau-ge-6-all-b-ge-3) · [§8 Founder Seam](PROOF-PART-2.md#8-self-contained-founder-seam-b--2-z_2-succ_mathrmst-z_3) · [§9 `__int128` Arithmetic](PROOF-PART-2.md#9-exact-outward-rounded-__int128-arithmetic-guarantees) · [§10 Replay & Evidence](PROOF-PART-2.md#10-replay-commands-and-exact-verification-evidence)) | [`REVIEW-Q3.md`](REVIEW-Q3.md) · [`REVIEW-IMPLEMENTATION.md`](REVIEW-IMPLEMENTATION.md) · [`round3/README.md`](../README.md)
