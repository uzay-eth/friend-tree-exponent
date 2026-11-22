> Produced by Gemini 4 Neon.

# Rigorous Refutation of Strogatz's 17/30 Exponent Conjecture and Sharp Bracket `[0.5650, 0.5665]`

## 1. Model and Main Theorem

Consider the **random friend tree** (Krapivsky–Redner isotropic redirection model with $`r = 1`$, arXiv:1706.01514 and arXiv:2403.20185):
- Start at $`n = 2`$ from a single edge $`T_2`$ on two vertices.
- Given $`T_n`$, pick a vertex $`V_n \in V(T_n)`$ uniformly at random, pick a neighbour $`W_n`$ of $`V_n`$ uniformly at random, and attach a new leaf $`n+1`$ to $`W_n`$ to form $`T_{n+1}`$.
- Let $`C_n := \#\{v \in V(T_n) : \deg_{T_n}(v) \ge 2\}`$ be the number of non-leaves (**core** vertices) of $`T_n`$.

Steven Strogatz (2022-09-20) highlighted the open problem of proving whether $`\mathbb{E}[C_n] = n^{\mu + o(1)}`$ with $`\mu = 17/30 = 0.5666\overline{6}`$. We rigorously **disprove** $`\mu = 17/30`$ and establish the sharpest known finite-$`n`$ bounds on $`\mathbb{E}[C_n]`$.

### Theorem 1 (Exact Finite-$`n`$ Bounds and Refutation of $`17/30`$)
For every integer $`n \ge 3`$,
```math
\frac{7{,}304{,}352{,}308}{8{,}763{,}936{,}965} \prod_{j=3}^{n-1}\left(1 + \frac{113}{200\, j}\right) \;\le\; \mathbb{E}[C_n] \;\le\; \frac{14{,}606{,}753{,}031}{2{,}989{,}767{,}809} \prod_{j=3}^{n-1}\left(1 + \frac{1133}{2000\, j}\right).
```
Consequently,
```math
0.5650 \;=\; \frac{113}{200} \;\le\; \liminf_{n\to\infty}\frac{\log \mathbb{E}[C_n]}{\log n} \;\le\; \limsup_{n\to\infty}\frac{\log \mathbb{E}[C_n]}{\log n} \;\le\; \frac{1133}{2000} \;=\; 0.5665 \;<\; \frac{17}{30} = 0.5666\overline{6}.
```
In particular, the conjectured exponent $`17/30`$ is strictly false (exceeding the upper bound by an exact margin of $`\frac{17}{30} - \frac{1133}{2000} = \frac{1}{6000}`$).

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

*(Remark: Theorem 1 and Corollary 2 are completely self-contained and do not assume existence of $`\lim \frac{\log \mathbb{E}[C_n]}{\log n}`$. Combined with the independent existence and almost-sure convergence theorems proved by the Astra lead—showing that a deterministic $`\mu`$ exists with $`\mathbb{E}[C_n] = n^{\mu+o(1)}`$ and $`\frac{\log C_n}{\log n} \to \mu`$ a.s.—Theorem 1 yields $`\mu \in \left[\frac{113}{200}, \frac{1133}{2000}\right] = [0.5650, 0.5665]`$.)*

---

## 2. Canonical Rooting and Exact Scaled Generator

At $`n = 3`$, $`T_3`$ is deterministically the $`3`$-vertex star $`K_{1,2}`$ with a unique degree-$`2`$ vertex $`r`$. For every $`n \ge 3`$, root $`T_n`$ at $`r`$.
Because vertex degrees never decrease, the core $`K(T_n) := \{v \in V(T_n) : \deg_{T_n}(v) \ge 2\}`$ is always a connected rooted subtree containing $`r`$. Every non-root core vertex $`u \in K(T_n) \setminus \{r\}`$ has a unique core parent $`p_u \in K(T_n)`$, and every leaf of $`T_n`$ is a child of a unique core vertex.

For each core vertex $`u \in K(T)`$, define its local triple $`(L_u, k_u, h_u)`$:
- $`h_u := 0`$ if $`u = r`$, and $`h_u := 1`$ if $`u \ne r`$ (parent indicator);
- $`L_u := \#\{\text{leaf children of } u\} \ge 0`$;
- $`k_u := \#\{\text{core children of } u\} \ge 0`$;
- $`D_u := \deg_T(u) = L_u + k_u + h_u \ge 2`$.

Let $`\vec{E}(K)`$ be the set of directed core edges $`(p \to w)`$ (parent $`p = p_w`$ to child $`w`$), and define the **credit** received by $`u \in K`$ from its core neighbours:
```math
\rho_u := \sum_{v \sim_K u} \frac{1}{D_v} = \frac{h_u}{D_{p_u}} + \sum_{w \in \mathrm{ch}_K(u)} \frac{1}{D_w}.
```
Note that since $`D_v \ge 2`$ for every $`v \in K`$, we always have
```math
0 \;\le\; \rho_u \;\le\; \frac{k_u + h_u}{2}, \qquad\text{and}\qquad L_u + \rho_u \;\le\; L_u + \frac{k_u + h_u}{2} \;=\; D_u - \frac{k_u + h_u}{2} \;<\; D_u + 1.
```

In one discrete step $`T_n \to T_{n+1}`$, the probability that the new leaf attaches to $`w \in V(T_n)`$ is $`\frac{1}{n}\sum_{v \sim w} \frac{1}{D_v}`$. Hence:
1. **Leaf gain at $`u \in K`$:** with probability $`\frac{L_u + \rho_u}{n}`$, a new leaf attaches to $`u`$, updating $`(L_u, k_u, h_u) \mapsto (L_u + 1, k_u, h_u)`$ ($`D_u \mapsto D_u + 1`$, no other vertex changes).
2. **Budding (leaf conversion) at $`u \in K`$ (when $`L_u \ge 1`$):** each of the $`L_u`$ leaf children of $`u`$ receives the new leaf with probability $`\frac{1}{n D_u}`$, so at total probability $`\frac{L_u}{n D_u}`$, one leaf child becomes a newborn core child $`\mathrm{nb}`$ with state $`(L, k, h) = (1, 0, 1)`$ ($`D = 2`$), while $`u`$ updates $`(L_u, k_u, h_u) \mapsto (L_u - 1, k_u + 1, h_u)`$ ($`D_u`$ **unchanged**).

For any functional $`\Phi`$ on finite rooted trees ($`n \ge 3`$), its **exact scaled drift** is therefore
```math
(\mathcal{L}\Phi)(T) := n\,\mathbb{E}\!\left[\Phi(T_{n+1}) - \Phi(T_n) \mid T_n = T\right] = \sum_{u \in K(T)} \left[(L_u + \rho_u)\,\Delta_{L_u}\Phi(T) + \frac{L_u}{D_u}\,\Delta_{\mathrm{bud}, u}\Phi(T)\right],
```
where $`\Delta_{L_u}\Phi(T)`$ and $`\Delta_{\mathrm{bud}, u}\Phi(T)`$ are the exact jumps under leaf gain and budding at $`u`$.

---

## 3. The Analytic $`1/D`$-Tail 3-Vertex + Cherry Potential $`\Phi(T)`$

Fix integer parameters satisfying the nesting conditions
```math
2 < D_0 \le L_F, \qquad 2 \le D_{cu} \le D_0, \qquad 2 \le D_{hu} \le D_0, \qquad D_0 \le K_F + 1, \qquad K_T + 1 \le K_F, \qquad K_{hu} \ge 1, \qquad K_{cu} \ge 2.
```
*(In both our upper and lower certificates: $`L_F = 800`$, $`K_F = 150`$, $`D_0 = 24`$, $`K_T = 12`$, $`D_{hu} = 16`$, $`K_{hu} = 8`$, $`D_{cu} = 12`$, $`K_{cu} = 6`$.)*

Define the finite state sets:
- **Vertex table states for $`f`$:** all $`(L, k, h) \in \{0, \dots, L_F\} \times \{0, \dots, K_F\} \times \{0, 1\}`$ with $`L + k + h \ge 2`$, plus two overflow states $`\mathrm{VSTAR}_0, \mathrm{VSTAR}_1`$ for $`k > K_F`$. Define the projection
  ```math
  \pi_f(L, k, h) := \begin{cases} \mathrm{VSTAR}_h & \text{if } k > K_F, \\ (\min(L, L_F), k, h) & \text{if } k \le K_F. \end{cases}
  ```
- **Inner and tail endpoint states for $`g, h_V, h_C`$:**
  - $`\mathcal{I} := \{(L, k, h) : h \in \{0, 1\},\, L, k \ge 0,\, 2 \le L + k + h < D_0\}`$ (**inner** states, exact degree $`D < D_0`$);
  - $`\mathrm{INF}`$ (single tail endpoint representing $`t = 0`$, i.e. $`D = \infty`$);
  - $`T_1(\kappa, h)`$ for $`\kappa \in \{0, \dots, K_T + 1\}`$ and $`h \in \{0, 1\}`$ (tail endpoints representing $`t = 1`$, i.e. $`D = D_0`$, with clamped child count $`\kappa(k) := \min(k, K_T + 1)`$).
  Set $`\mathcal{G} := \mathcal{I} \cup \{\mathrm{INF}\} \cup \{T_1(\kappa, h)\}`$. Let $`\mathcal{P} \subset \mathcal{G}`$ be valid **parent** states ($`\mathrm{INF}`$, plus all inner/T1 states with $`k \ge 1`$ or $`\kappa \ge 1`$) and $`\mathcal{W} \subset \mathcal{G}`$ be valid **child** states ($`\mathrm{INF}`$, plus all inner/T1 states with $`h = 1`$).
- **Middle-vertex domains for $`h_V`$ and $`h_C`$:**
  ```math
  H_u := \{(L, k, 1) \in \mathcal{I} : 1 \le k \le K_{hu},\, D < D_{hu}\}, \qquad C_u := \{(L, k, h) \in \mathcal{I} : 2 \le k \le K_{cu},\, D < D_{cu}\}.
  ```
  Because $`D_u`$ and $`k_u`$ never decrease under either transition, once a vertex exits $`H_u`$ (or $`C_u`$) it can **never re-enter**. A vertex can enter $`H_u`$ only when $`(L, 0, 1)`$ with $`D < D_{hu}`$ buds to $`(L-1, 1, 1) \in H_u`$; it can enter $`C_u`$ only when $`(L, 1, h)`$ with $`D < D_{cu}`$ buds to $`(L-1, 2, h) \in C_u`$.

### Analytic $`1/D`$-Tail Interpolation
For every core vertex $`v \in K(T)`$, define its formal probability weight vector $`\mu_v`$ on $`\mathcal{G}`$:
- If $`D_v < D_0`$ (inner): $`\mu_v := \delta_{(L_v, k_v, h_v)}`$.
- If $`D_v \ge D_0`$ (tail): let $`t_v := \frac{D_0}{D_v} \in (0, 1]`$ and $`\kappa_v := \min(k_v, K_T + 1)`$, and set
  ```math
  \mu_v := (1 - t_v)\,\delta_{\mathrm{INF}} + t_v\,\delta_{T_1(\kappa_v, h_v)}.
  ```
Given rational tables:
- $`f : \mathrm{dom}(\pi_f) \to \mathbb{Q}`$,
- $`g : \mathcal{P} \times \mathcal{W} \to \mathbb{Q}`$ with **$`g(p, w) = 0`$ whenever both $`p, w \notin \mathcal{I}`$** (tail-tail edges vanish),
- $`h_V : \mathcal{P} \times H_u \times \mathcal{W} \to \mathbb{Q}`$,
- $`h_C : C_u \times \mathcal{W} \times \mathcal{W} \to \mathbb{Q}`$ symmetric in its last two arguments ($`h_C(u; a, b) = h_C(u; b, a)`$),

extend $`g, h_V, h_C`$ to all core vertices via multilinear evaluation against $`\mu`$:
```math
g(p \to w) := \sum_{a, b \in \mathcal{G}} \mu_p(a)\,\mu_w(b)\,g[a, b],
```
```math
h_V(p \to u \to w) := \mathbf{1}_{u \in H_u} \sum_{a, b \in \mathcal{G}} \mu_p(a)\,\mu_w(b)\,h_V[a, u, b],
```
```math
h_C(u; \{w_1, w_2\}) := \mathbf{1}_{u \in C_u} \sum_{a, b \in \mathcal{G}} \mu_{w_1}(a)\,\mu_{w_2}(b)\,h_C[u, a, b].
```
Finally, define the **potential functional** on any rooted tree $`T`$ ($`n \ge 3`$):
```math
\Phi(T) := \sum_{u \in K(T)} f(\pi_f(L_u, k_u, h_u)) + \sum_{(p \to w) \in \vec{E}(K)} g(p \to w) + \sum_{\substack{p \to u \to w \\ \text{vertical 2-path}}} h_V(p \to u \to w) + \sum_{u \in K(T)} \sum_{\{w_1, w_2\} \subset \mathrm{ch}_K(u)} h_C(u; \{w_1, w_2\}).
```

---

## 4. Exact Algebraic Drift Decomposition and Endpoint Reduction

The certificate also stores auxiliary rational centering parameters $`r_{\mathrm{par}}[u], r_{\mathrm{ch}}[u], x_{\mathrm{par}}[u], x_{\mathrm{ch}}[u], x_{cc}[u] \in \mathbb{Q}`$ for $`u \in \mathcal{I}`$, rational knot abscissae
```math
0 = \xi_0 < \xi_1 = \frac{1}{D_0} < \xi_2 < \dots < \xi_{m-1} = \frac{1}{2},
```
and rational knot value vectors $`k_{pc}[u], k_{cp}[u], k_{cc}[u] \in \mathbb{Q}^m`$ for $`u \in \mathcal{I}`$ defining piecewise-linear functions $`A_{pc}[u], A_{cp}[u], A_{cc}[u] : [0, 1/2] \to \mathbb{R}`$ (rational-valued on $`[0, 1/2] \cap \mathbb{Q}`$) by linear interpolation between consecutive $`(\xi_j, a_j)`$.

### Lemma 1 (Exact Piecewise-Linear Fenchel–Young Tautology)
Let $`A : [0, 1/2] \to \mathbb{R}`$ be the piecewise-linear interpolant of $`(\xi_j, a_j)_{j=0}^{m-1}`$ with $`\xi_0 = 0, \xi_1 = 1/D_0, \xi_{m-1} = 1/2`$. Then:
1. $`A`$ is **exact affine** on $`[0, 1/D_0]`$: for every $`t \in [0, 1]`$, $`A(t/D_0) = (1-t) A(0) + t A(1/D_0)`$.
2. For all $`x \in [0, 1/2]`$ and all $`y \in \mathbb{R}`$,
   ```math
   A(x) + \min_{0 \le j < m}(\xi_j y - a_j) \;\le\; x y \;\le\; A(x) + \max_{0 \le j < m}(\xi_j y - a_j).
   ```
3. $`y \mapsto \max_j(\xi_j y - a_j)`$ is convex and $`y \mapsto \min_j(\xi_j y - a_j)`$ is concave.

*Proof.* Write $`x = (1-\alpha)\xi_i + \alpha\xi_{i+1}`$ with $`\alpha \in [0, 1]`$ on its containing subinterval $`[\xi_i, \xi_{i+1}]`$. Then $`A(x) = (1-\alpha)a_i + \alpha a_{i+1}`$, so $`xy - A(x) = (1-\alpha)(\xi_i y - a_i) + \alpha(\xi_{i+1} y - a_{i+1})`$, which lies between $`\min_j(\xi_j y - a_j)`$ and $`\max_j(\xi_j y - a_j)`$. Properties (1) and (3) are immediate. $`\square`$

### Lemma 2 (Exact Analytic $`1/D`$-Tail Leaf-Gain Identity)
For any tail vertex $`p \in K`$ with $`D_p \ge D_0`$ and any function $`F(\mu_p) = (1 - t_p) F[\mathrm{INF}] + t_p F[T_1(\kappa_p, h_p)]`$ with $`t_p = D_0/D_p`$:
1. Under leaf gain at $`p`$ ($`D_p \mapsto D_p + 1`$, $`\kappa_p`$ and $`h_p`$ unchanged),
   ```math
   (L_p + \rho_p)\,\Delta_{L_p} F(\mu_p) = t_p \cdot \gamma_p \cdot \big(F[\mathrm{INF}] - F[T_1(\kappa_p, h_p)]\big), \qquad \gamma_p := \frac{L_p + \rho_p}{D_p + 1}.
   ```
   Moreover, $`\gamma_p`$ always lies in the closed interval $`[\gamma_{\min}(\kappa_p, h_p),\, 1]`$, where
   ```math
   \gamma_{\min}(\kappa, h) := \begin{cases} \dfrac{\max(0, D_0 - \kappa - h)}{D_0 + 1} & \text{if } \kappa \le K_T, \\ 0 & \text{if } \kappa = K_T + 1. \end{cases}
   ```
2. Under budding at $`p`$ ($`D_p`$ unchanged, so $`t_p`$ unchanged; $`k_p \mapsto k_p + 1`$ at rate $`b_p = L_p/D_p`$):
   - If $`k_p > K_T`$ (so $`\kappa_p = K_T + 1`$), then $`\kappa_p`$ is unchanged and $`\Delta_{\mathrm{bud}, p} F(\mu_p) = 0`$.
   - If $`k_p \le K_T`$, then $`b_p\,\Delta_{\mathrm{bud}, p} F(\mu_p) = t_p \cdot b_p \cdot \big(F[T_1(\kappa_p + 1, h_p)] - F[T_1(\kappa_p, h_p)]\big)`$, with $`b_p \in [b_{\min}(\kappa_p, h_p),\, 1]`$, where $`b_{\min}(\kappa, h) := \frac{\max(0, D_0 - \kappa - h)}{D_0}`$.

*Proof.* (1) $`\Delta_{L_p} t_p = \frac{D_0}{D_p+1} - \frac{D_0}{D_p} = -\frac{t_p}{D_p+1}`$. For the upper bound on $`\gamma_p`$, $`L_p + \rho_p \le D_p - \frac{k_p+h_p}{2} < D_p + 1 \implies \gamma_p < 1`$. For the lower bound, $`\rho_p \ge 0`$ gives $`\gamma_p \ge \frac{L_p}{D_p+1} = 1 - \frac{k_p+h_p+1}{D_p+1}`$; when $`k_p \le K_T`$, this is nondecreasing in $`D_p \ge D_0`$, hence minimized at $`D_p = D_0`$. (2) Analogous using $`b_p = 1 - \frac{k_p+h_p}{D_p}`$. $`\square`$

### Lemma 3 (Quadratic Credit Bound on $`[0, 1]`$ and $`[0, 1]^2`$)
1. For all $`t \in [0, 1]`$ and $`d_0, d_1 \in \mathbb{R}`$,
   ```math
   t \cdot \frac{\min(d_0, d_1)}{D_0} \;\le\; \frac{t}{D_0}\big((1 - t) d_0 + t d_1\big) \;\le\; t \cdot \frac{\max(d_0, d_1)}{D_0}.
   ```
2. For all $`(t_1, t_2) \in [0, 1]^2`$ and bi-affine $`\Delta(t_1, t_2) = (1-t_1)(1-t_2)d_{00} + t_1(1-t_2)d_{10} + (1-t_1)t_2 d_{01} + t_1 t_2 d_{11}`$,
   ```math
   \frac{t_1}{D_0}\Delta(t_1, t_2) \;\le\; (1-t_2)t_1 \frac{\max(d_{00}, d_{10})}{D_0} + t_2 t_1 \frac{\max(d_{01}, d_{11})}{D_0},
   ```
   and symmetrically for $`\frac{t_2}{D_0}\Delta(t_1, t_2)`$ (with $`\max`$ replaced by $`\min`$ for the lower bound).

*Proof.* (1) Immediate since $`t/D_0 \ge 0`$ and $`(1-t)d_0 + t d_1 \in [\min(d_0, d_1), \max(d_0, d_1)]`$. (2) Write $`\Delta(t_1, t_2) = (1-t_2)((1-t_1)d_{00} + t_1 d_{10}) + t_2((1-t_1)d_{01} + t_1 d_{11})`$ and apply (1) to each bracket. $`\square`$

### Lemma 4 (Exact Tree-Level Decomposition of $`(\mathcal{L}\Phi)(T)`$)
For every finite rooted tree $`T`$ ($`n \ge 3`$), expand $`(\mathcal{L}\Phi)(T)`$ into vertex, edge, vertical 2-path, and cherry terms by:
1. **Vertex-to-edge credit transfer** (only for inner $`u \in \mathcal{I}`$): write
   ```math
   \rho_u \Delta_{L_u} f(u) = \big(h_u r_{\mathrm{par}}[u] + k_u r_{\mathrm{ch}}[u]\big)\Delta_{L_u} f(u) + h_u\left(\frac{1}{D_{p_u}} - r_{\mathrm{par}}[u]\right)\Delta_{L_u} f(u) + \sum_{w \in \mathrm{ch}_K(u)}\left(\frac{1}{D_w} - r_{\mathrm{ch}}[u]\right)\Delta_{L_u} f(u),
   ```
   assigning the first term to vertex $`u`$, the second to incoming edge $`(p_u \to u)`$, and each summand of the third to outgoing edge $`(u \to w)`$. (Notice: for $`D_u \ge D_0`$, no credit is transferred; vertex $`u`$ directly bounds $`\rho_u \in [0, (h_u + k_u)/2]`$.)
2. **Edge-to-3-vertex credit transfer** (only for inner $`u \in \mathcal{I}`$, because for tail $`u`$, $`\Delta_{L_u} g`$ is handled directly by Lemma 2):
   - When $`u`$ is the middle vertex of a vertical 2-path $`(p \to u \to w)`$:
     - If $`u \in H_u`$, split $`\frac{1}{D_p} \Delta_{L_u} g(u \to w) = x_{\mathrm{par}}[u]\Delta_{L_u} g(u \to w) + (\frac{1}{D_p} - x_{\mathrm{par}}[u])\Delta_{L_u} g(u \to w)`$ and $`\frac{1}{D_w}\Delta_{L_u} g(p \to u) = x_{\mathrm{ch}}[u]\Delta_{L_u} g(p \to u) + (\frac{1}{D_w} - x_{\mathrm{ch}}[u])\Delta_{L_u} g(p \to u)`$, keeping the $`x`$-terms on edges $`(u \to w)`$ and $`(p \to u)`$ and moving the residual source terms to $`h_V(p \to u \to w)`$.
     - If $`u \notin H_u`$ ($`h_u = 1, k_u \ge 1`$), decouple by Lemma 1:
       ```math
       \frac{1}{D_p}\Delta_{L_u} g(u \to w) \lessgtr A_{pc}[u](1/D_p) + A_{pc}^*(\Delta_{L_u} g(u \to w)), \qquad \frac{1}{D_w}\Delta_{L_u} g(p \to u) \lessgtr A_{cp}[u](1/D_w) + A_{cp}^*(\Delta_{L_u} g(p \to u)),
       ```
       assigning $`k_u A_{pc}[u](1/D_p)`$ and $`k_u A_{cp}^*(\Delta_{L_u} g(p \to u))`$ to edge $`(p \to u)`$, and $`A_{cp}[u](1/D_w)`$ and $`A_{pc}^*(\Delta_{L_u} g(u \to w))`$ to edge $`(u \to w)`$.
   - When $`u`$ is the apex of a cherry $`\{w_1, w_2\} \subset \mathrm{ch}_K(u)`$ ($`k_u \ge 2`$):
     - If $`u \in C_u`$, split with $`x_{cc}[u]`$: $`(k_u - 1) x_{cc}[u] \Delta_{L_u} g(u \to w)`$ stays on edge $`(u \to w)`$, and $`(\frac{1}{D_{w_1}} - x_{cc}[u])\Delta_{L_u} g(u \to w_2) + (\frac{1}{D_{w_2}} - x_{cc}[u])\Delta_{L_u} g(u \to w_1)`$ moves to $`h_C(u; \{w_1, w_2\})`$.
     - If $`u \notin C_u`$, decouple by Lemma 1 with $`A_{cc}[u]`$: edge $`(u \to w)`$ receives $`(k_u - 1)\big(A_{cc}[u](1/D_w) + A_{cc}^*(\Delta_{L_u} g(u \to w))\big)`$.
3. **Newborn terms on budding at $`u`$ ($`L_u \ge 1`$, rate $`L_u/D_u`$):**
   Let $`u' = (L_u - 1, k_u + 1, h_u)`$ and $`\mathrm{nb} = (1, 0, 1)`$.
   - $`+f(\mathrm{nb}) + g(u' \to \mathrm{nb})`$ is assigned to vertex $`u`$.
   - If $`h_u = 1`$ and $`u' \in H_u`$ (note $`D_u = D_{u'} < D_{hu} \le D_0`$, so $`u \in \mathcal{I}`$), the newborn vertical 2-path $`+h_V(p_u \to u' \to \mathrm{nb})`$ is assigned to incoming edge $`(p_u \to u)`$.
   - If $`u' \in C_u`$ (so $`D_u = D_{u'} < D_{cu} \le D_0`$, $`u \in \mathcal{I}`$), for each existing core child $`w`$ of $`u`$, the newborn cherry $`+h_C(u'; \{w, \mathrm{nb}\})`$ is assigned to outgoing edge $`(u \to w)`$.
   - No other $`h_V`$ or $`h_C`$ is created because $`\mathrm{nb}`$ has $`k_{\mathrm{nb}} = 0`$, existing vertical paths and cherries through $`u`$ simply transition from $`u`$ to $`u'`$ (or vanish if $`u'`$ exits $`H_u`$ or $`C_u`$), and every newborn vertical path $`(p_u \to u' \to \mathrm{nb})`$ or newborn cherry $`(u'; \{w, \mathrm{nb}\})`$ created whenever $`u' \in H_u`$ or $`u' \in C_u`$ (including both first entries $`k_u = 0 \to 1`$, $`k_u = 1 \to 2`$ and subsequent buddings while remaining in $`H_u`$ or $`C_u`$) is accounted for by the two edge terms above.
4. **External leaf-gain rates in $`h_V`$ and $`h_C`$:**
   In $`h_V(p \to u \to w)`$: for inner $`p`$, $`L_p + \rho_p \in [L_p + 1/D_u,\, L_p + 1/D_u + (h_p+k_p-1)/2]`$; for inner $`w`$, $`L_w + \rho_w \in [L_w + 1/D_u,\, L_w + 1/D_u + k_w/2]`$; for $`u`$, $`L_u + \rho_u \in [L_u + 1/D_p + 1/D_w,\, L_u + 1/D_p + 1/D_w + (k_u-1)/2]`$ (with $`1/D_p = t_p/D_0`$ or $`1/D_w = t_w/D_0`$ handled via Lemma 3 whenever $`p`$ or $`w`$ is tail!). Analogously in $`h_C(u; \{w_1, w_2\})`$, $`u`$ has external core neighbours $`h_u + k_u - 2`$.

Applying Lemmas 1–3 to every tail parameter $`t_p, t_w \in [0, 1]`$ reduces the infinite-tree drift inequality $`(\mathcal{L}\Phi)(T) \le \frac{P}{Q}\Phi(T)`$ (for UPPER; reversed $`\ge`$ for LOWER) to a **finite list** of exact rational inequalities. Write $`\operatorname{ext} = \max`$ and $`\lessgtr`$ for non-strict $`\le`$ in UPPER, and $`\operatorname{ext} = \min`$ and $`\lessgtr`$ for non-strict $`\ge`$ in LOWER, with $`\lambda = P/Q`$.

### 4.5 Explicit List of Finite Rational Inequalities

Let $`\mathrm{nb} = (1, 0, 1) \in \mathcal{I}`$. For any $`v \in \mathcal{G}`$, define $`\mathrm{inv}(v) = 1/D_v`$ for $`v \in \mathcal{I}`$, $`\mathrm{inv}(\mathrm{INF}) = 0`$, $`\mathrm{inv}(T_1(\kappa, h)) = 1/D_0`$; define $`v^+ := (L_v+1, k_v, h_v)`$ if $`D_v + 1 < D_0`$, $`T_1(\min(k_v, K_T+1), h_v)`$ if $`D_v = D_0 - 1`$, and $`\mathrm{INF}`$ if $`v \notin \mathcal{I}`$; define $`v^{\mathrm{bud}} := (L_v-1, k_v+1, h_v)`$ if $`v \in \mathcal{I}, L_v \ge 1`$, and $`T_1(\kappa_v+1, h_v)`$ if $`v = T_1(\kappa_v, h_v)`$ with $`\kappa_v \le K_T`$. For $`u \in \mathcal{I}`$, put $`\Delta_L f(u) := f(L_u+1, k_u, h_u) - f(L_u, k_u, h_u)`$ (with $`f`$ always evaluated via $`\pi_f`$).

1. **Vertex inequalities (`f`):**
   - **$`\mathrm{VSTAR}_h`$ ($`h \in \{0, 1\}`$, $`k > K_F`$, $`D \ge K_F + 1 + h \ge D_0`$):** for both $`b \in \{0, 1\}`$ and both $`g_{\mathrm{nb}} \in \{g[\mathrm{INF}, \mathrm{nb}],\, g((0, K_F+1, h) \to \mathrm{nb})\}`$,
     ```math
     b\big(f(\mathrm{nb}) + g_{\mathrm{nb}}\big) \;\lessgtr\; \lambda\, f(\mathrm{VSTAR}_h).
     ```
   - **Box boundary $`L = L_F`$ ($`0 \le k \le K_F, h \in \{0, 1\}`$):**
     - Exact $`L = L_F`$: $`\frac{L_F}{L_F+k+h}\big(f(\pi_f(L_F-1, k+1, h)) - f(L_F, k, h) + f(\mathrm{nb}) + g((L_F-1, k+1, h) \to \mathrm{nb})\big) \lessgtr \lambda\, f(L_F, k, h)`$.
     - Tail $`L \ge L_F + 1`$: for both $`b \in \left\{\frac{L_F+1}{L_F+1+k+h},\, 1\right\}`$ and both $`g_{\mathrm{nb}} \in \{g[\mathrm{INF}, \mathrm{nb}],\, g((L_F, k+1, h) \to \mathrm{nb})\}`$,
       ```math
       b\big(f(\pi_f(L_F, k+1, h)) - f(L_F, k, h) + f(\mathrm{nb}) + g_{\mathrm{nb}}\big) \;\lessgtr\; \lambda\, f(L_F, k, h).
       ```
   - **Inner box $`L < L_F, D = L+k+h \in [2, D_0-1]`$:** let $`u = (L, k, h) \in \mathcal{I}`$, $`B(u) := \mathbf{1}_{L \ge 1}\frac{L}{D}\big(f(L-1, k+1, h) - f(u) + f(\mathrm{nb}) + g((L-1, k+1, h) \to \mathrm{nb})\big)`$. Check
     ```math
     \big(L + h\,r_{\mathrm{par}}[u] + k\,r_{\mathrm{ch}}[u]\big)\Delta_L f(u) + B(u) \;\lessgtr\; \lambda\, f(u).
     ```
   - **Outer box $`L < L_F, k \le K_F, D \ge D_0`$:** for both $`r \in \{L,\, L + (h+k)/2\}`$,
     ```math
     r\,\big(f(L+1, k, h) - f(L, k, h)\big) + B(L, k, h) \;\lessgtr\; \lambda\, f(L, k, h).
     ```

2. **Directed edge inequalities (`g`):** First require $`g[p, w] = 0`$ for all $`p, w \notin \mathcal{I}`$. For every $`(p, w) \in \mathcal{P} \times \mathcal{W}`$ with at least one of $`p, w \in \mathcal{I}`$, let $`v = g[p, w]`$, $`\Delta_p g = g[p^+, w] - v`$, $`\Delta_w g = g[p, w^+] - v`$.
   - Precomputed PL additive terms $`A_{\mathrm{tot}}[p, w]`$:
     - If $`w \in \mathcal{I} \setminus H_u`$ with $`h_w = 1, k_w \ge 1`$: add $`k_w A_{pc}[w](\mathrm{inv}(p))`$.
     - If $`p \in \mathcal{I} \setminus H_u`$ with $`h_p = 1, k_p \ge 1`$: add $`A_{cp}[p](\mathrm{inv}(w))`$.
     - If $`p \in \mathcal{I} \setminus C_u`$ with $`k_p \ge 2`$: add $`(k_p - 1) A_{cc}[p](\mathrm{inv}(w))`$.
   - Parent ($`p`$) terms:
     - If $`p \in \mathcal{I}`$: add $`(\mathrm{inv}(w) - r_{\mathrm{ch}}[p])\Delta_L f(p) + \big(L_p + \mathbf{1}_{w \in \mathcal{I}}\mathrm{inv}(w) + \mathbf{1}_{p \in H_u}x_{\mathrm{par}}[p] + \mathbf{1}_{p \in C_u}(k_p-1)x_{cc}[p]\big)\Delta_p g`$, plus if $`L_p \ge 1`$: $`\frac{L_p}{D_p}\big(g[p^{\mathrm{bud}}, w] - v + \mathbf{1}_{p^{\mathrm{bud}} \in C_u} h_C[p^{\mathrm{bud}}, w, \mathrm{nb}]\big)`$.
     - If $`p = T_1(\kappa_p, h_p)`$: add $`\operatorname{ext}_{r \in \{\gamma_{\min}(\kappa_p, h_p), 1\}} r\,\Delta_p g`$, plus if $`\kappa_p \le K_T`$: $`\operatorname{ext}_{b \in \{b_{\min}(\kappa_p, h_p), 1\}} b\,(g[p^{\mathrm{bud}}, w] - v)`$, plus quadratic term $`\frac{1}{D_0}\operatorname{ext}\big(g[\mathrm{INF}, w^+] - g[\mathrm{INF}, w],\, \Delta_w g\big)`$.
     - If $`p = \mathrm{INF}`$: no $`p`$-transition terms.
   - Child ($`w`$) terms (symmetric):
     - If $`w \in \mathcal{I}`$: add $`(\mathrm{inv}(p) - r_{\mathrm{par}}[w])\Delta_L f(w) + \big(L_w + \mathbf{1}_{p \in \mathcal{I}}\mathrm{inv}(p) + \mathbf{1}_{w \in H_u} k_w x_{\mathrm{ch}}[w]\big)\Delta_w g`$, plus if $`L_w \ge 1`$: $`\frac{L_w}{D_w}\big(g[p, w^{\mathrm{bud}}] - v + \mathbf{1}_{w^{\mathrm{bud}} \in H_u} h_V[p, w^{\mathrm{bud}}, \mathrm{nb}]\big)`$.
     - If $`w = T_1(\kappa_w, 1)`$: add $`\operatorname{ext}_{r \in \{\gamma_{\min}(\kappa_w, 1), 1\}} r\,\Delta_w g`$, plus if $`\kappa_w \le K_T`$: $`\operatorname{ext}_{b \in \{b_{\min}(\kappa_w, 1), 1\}} b\,(g[p, w^{\mathrm{bud}}] - v)`$, plus quadratic term $`\frac{1}{D_0}\operatorname{ext}\big(g[p^+, \mathrm{INF}] - g[p, \mathrm{INF}],\, \Delta_p g\big)`$.
     - If $`w = \mathrm{INF}`$: no $`w`$-transition terms.
   - Fenchel–Young conjugate terms $`A^*(y) := \operatorname{ext}_{0 \le j < m}(\xi_j y - a_j)`$:
     - If $`p \in \mathcal{I} \setminus H_u`$ with $`h_p = 1`$: add $`A^*_{pc}[p](\Delta_p g)`$.
     - If $`p \in \mathcal{I} \setminus C_u`$ with $`k_p \ge 2`$: add $`(k_p - 1) A^*_{cc}[p](\Delta_p g)`$.
     - If $`w \in \mathcal{I} \setminus H_u`$ with $`k_w \ge 1`$: add $`k_w A^*_{cp}[w](\Delta_w g)`$.
   - Check total $`\mathrm{LHS}_g(p, w) \lessgtr \lambda\, g[p, w]`$.

3. **Vertical 2-path inequalities (`h_V`):** For every $`(p, u, w) \in \mathcal{P} \times H_u \times \mathcal{W}`$, let $`v = h_V[p, u, w]`$ (interpreting $`h_V[\cdot, u', \cdot] = 0`$ if $`u' \notin H_u`$):
   - Source: $`(\mathrm{inv}(p) - x_{\mathrm{par}}[u])(g[u^+, w] - g[u, w]) + (\mathrm{inv}(w) - x_{\mathrm{ch}}[u])(g[p, u^+] - g[p, u])`$.
   - $`p`$-terms: if $`p \in \mathcal{I}`$, add $`\operatorname{ext}_{r \in \{L_p + 1/D_u,\, L_p + 1/D_u + (h_p+k_p-1)/2\}} r\,(h_V[p^+, u, w] - v)`$ plus $`\mathbf{1}_{L_p \ge 1}\frac{L_p}{D_p}(h_V[p^{\mathrm{bud}}, u, w] - v)`$. If $`p = T_1(\kappa_p, h_p)`$, add $`\operatorname{ext}_{r \in \{\gamma_{\min}, 1\}} r(h_V[\mathrm{INF}, u, w] - v)`$, plus $`\mathbf{1}_{\kappa_p \le K_T}\operatorname{ext}_{b \in \{b_{\min}, 1\}} b(h_V[p^{\mathrm{bud}}, u, w] - v)`$, plus $`\frac{1}{D_0}\operatorname{ext}\big(h_V[\mathrm{INF}, u^+, w] - h_V[\mathrm{INF}, u, w],\, h_V[p, u^+, w] - v\big)`$.
   - $`w`$-terms: if $`w \in \mathcal{I}`$, add $`\operatorname{ext}_{r \in \{L_w + 1/D_u,\, L_w + 1/D_u + k_w/2\}} r\,(h_V[p, u, w^+] - v)`$ plus $`\mathbf{1}_{L_w \ge 1}\frac{L_w}{D_w}(h_V[p, u, w^{\mathrm{bud}}] - v)`$. If $`w = T_1(\kappa_w, 1)`$, add $`\operatorname{ext}_{r \in \{\gamma_{\min}, 1\}} r(h_V[p, u, \mathrm{INF}] - v)`$, plus $`\mathbf{1}_{\kappa_w \le K_T}\operatorname{ext}_{b \in \{b_{\min}, 1\}} b(h_V[p, u, w^{\mathrm{bud}}] - v)`$, plus $`\frac{1}{D_0}\operatorname{ext}\big(h_V[p, u^+, \mathrm{INF}] - h_V[p, u, \mathrm{INF}],\, h_V[p, u^+, w] - v\big)`$.
   - $`u`$-terms: add $`\operatorname{ext}_{r \in \{r_0,\, r_0 + (k_u - 1)/2\}} r\,(h_V[p, u^+, w] - v)`$ where $`r_0 = L_u + \mathbf{1}_{p \in \mathcal{I}}\mathrm{inv}(p) + \mathbf{1}_{w \in \mathcal{I}}\mathrm{inv}(w)`$, plus $`\mathbf{1}_{L_u \ge 1}\frac{L_u}{D_u}(h_V[p, u^{\mathrm{bud}}, w] - v)`$.
   - Check total $`\mathrm{LHS}_{h_V}(p, u, w) \lessgtr \lambda\, h_V[p, u, w]`$.

4. **Cherry inequalities (`h_C`):** For every $`u \in C_u`$ and unordered $`\{w_1, w_2\} \subset \mathcal{W}`$ (including $`w_1 = w_2`$), let $`v = h_C[u, w_1, w_2]`$ (symmetric in $`w_1, w_2`$, and $`0`$ when $`u' \notin C_u`$):
   - Source: $`(\mathrm{inv}(w_1) - x_{cc}[u])(g[u^+, w_2] - g[u, w_2]) + (\mathrm{inv}(w_2) - x_{cc}[u])(g[u^+, w_1] - g[u, w_1])`$.
   - $`w_1`$-terms (and symmetrically $`w_2`$-terms): if $`w_1 \in \mathcal{I}`$, add $`\operatorname{ext}_{r \in \{L_1 + 1/D_u,\, L_1 + 1/D_u + k_1/2\}} r(h_C[u, w_1^+, w_2] - v)`$ plus $`\mathbf{1}_{L_1 \ge 1}\frac{L_1}{D_1}(h_C[u, w_1^{\mathrm{bud}}, w_2] - v)`$. If $`w_1 = T_1(\kappa_1, 1)`$, add $`\operatorname{ext}_{r \in \{\gamma_{\min}, 1\}} r(h_C[u, \mathrm{INF}, w_2] - v)`$, plus $`\mathbf{1}_{\kappa_1 \le K_T}\operatorname{ext}_{b \in \{b_{\min}, 1\}} b(h_C[u, w_1^{\mathrm{bud}}, w_2] - v)`$, plus $`\frac{1}{D_0}\operatorname{ext}\big(h_C[u^+, \mathrm{INF}, w_2] - h_C[u, \mathrm{INF}, w_2],\, h_C[u^+, w_1, w_2] - v\big)`$.
   - $`u`$-terms: add $`\operatorname{ext}_{r \in \{r_0,\, r_0 + (h_u + k_u - 2)/2\}} r(h_C[u^+, w_1, w_2] - v)`$ where $`r_0 = L_u + \mathbf{1}_{w_1 \in \mathcal{I}}\mathrm{inv}(w_1) + \mathbf{1}_{w_2 \in \mathcal{I}}\mathrm{inv}(w_2)`$, plus $`\mathbf{1}_{L_u \ge 1}\frac{L_u}{D_u}(h_C[u^{\mathrm{bud}}, w_1, w_2] - v)`$.
   - Check total $`\mathrm{LHS}_{h_C}(u, w_1, w_2) \lessgtr \lambda\, h_C[u, w_1, w_2]`$.

### Lemma 5 (Global Pointwise Comparison $`\Phi(T)`$ vs. $`C(T)`$)
Assign each edge $`g(p \to w)`$ and each vertical 2-path $`h_V(p \to u \to w)`$ to its unique bottom vertex $`w`$ ($`h_w = 1`$), and assign each cherry $`h_C(u; \{w_1, w_2\})`$ to its middle vertex $`u \in C_u`$. Every non-root core vertex $`w`$ has a unique parent and at most one grandparent, hence receives at most one $`g`$ term and at most one $`h_V`$ term; every $`u \in C_u`$ receives $`\binom{k_u}{2}`$ cherry terms.
Define for $`w \in \mathcal{W}`$ and $`u \in C_u`$ (using $`\min`$ with $`0`$ for UPPER, $`\max`$ with $`0`$ for LOWER):
```math
g_{\mathrm{ext}}(w) := \operatorname*{ext}_{p \in \mathcal{P}} \big(0, g[p, w]\big), \qquad h_{V,\mathrm{ext}}(w) := \operatorname*{ext}_{p \in \mathcal{P},\, u \in H_u} \big(0, h_V[p, u, w]\big), \qquad h_{C,\mathrm{ext}}(u) := \operatorname*{ext}_{w_1, w_2 \in \mathcal{W}} \big(0, h_C[u, w_1, w_2]\big).
```
Since multilinear interpolation on $`[0, 1]`$ or $`[0, 1]^2`$ attains its extrema at the corners $`\{\mathrm{INF}, T_1\}`$, for every core vertex $`u \in K(T)`$ its grouped contribution $`W(u)`$ satisfies (with $`\ge`$ for UPPER, $`\le`$ for LOWER):
```math
W(u) \;\gtrless\; f(\pi_f(L_u, k_u, h_u)) + \mathbf{1}_{h_u=1}\big(g_{\mathrm{ext}}^{\text{eff}}(u) + h_{V,\mathrm{ext}}^{\text{eff}}(u)\big) + \mathbf{1}_{u \in C_u} \binom{k_u}{2} h_{C,\mathrm{ext}}(u),
```
where for tail $`u`$ ($`D_u \ge D_0`$), $`g_{\mathrm{ext}}^{\text{eff}}(u) = \operatorname{ext}(g_{\mathrm{ext}}(\mathrm{INF}), g_{\mathrm{ext}}(T_1(\kappa_u, 1)))`$ (and similarly for $`h_V`$).
Taking the finite minimum $`c_0 > 0`$ (UPPER) or finite maximum $`C_0 > 0`$ (LOWER) over all $`(L, k, h) \in [0..L_F]\times[0..K_F]\times\{0,1\}`$ ($`D \ge 2`$) and $`\mathrm{VSTAR}_0, \mathrm{VSTAR}_1`$ gives
```math
\Phi_{\mathrm{up}}(T) \ge c_0\, C(T), \qquad \Phi_{\mathrm{lo}}(T) \le C_0\, C(T) \qquad (\forall\, T \text{ with } n \ge 3).
```

### Completion of Proof of Theorem 1
Since $`\mathbb{E}[\Phi(T_{n+1}) \mid T_n] = \Phi(T_n) + \frac{1}{n}(\mathcal{L}\Phi)(T_n)`$, taking expectations and applying Lemma 4 gives
```math
\mathbb{E}[\Phi_{\mathrm{up}}(T_{n+1})] \le \left(1 + \frac{1133}{2000\, n}\right)\mathbb{E}[\Phi_{\mathrm{up}}(T_n)], \qquad \mathbb{E}[\Phi_{\mathrm{lo}}(T_{n+1})] \ge \left(1 + \frac{113}{200\, n}\right)\mathbb{E}[\Phi_{\mathrm{lo}}(T_n)].
```
Inducting from $`n = 3`$ (where $`T_3`$ is deterministic with $`\Phi(T_3) = f(2, 0, 0)`$) and applying Lemma 5 yields Theorem 1 with exact rational prefactors $`A_{\mathrm{up}} = \Phi_{\mathrm{up}}(T_3)/c_0`$ and $`A_{\mathrm{lo}} = \Phi_{\mathrm{lo}}(T_3)/C_0`$. $`\square`$

---

## 5. Exact Binary Certificate Format (for Independent Verifiers)

Both certificates (`upper_5665.bin` and `lower_5650.bin`) are generated deterministically by `cert/cert_build.c` in the binary format defined in `cert/cert_format.h` (little-endian x86_64):
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

---

## 6. Exact Verification Summary (Arbitrary-Precision GMP)

| Certificate | Target Exponent | Exact Checks (`f / g / hV / hC`) | Failures | Exact Prefactor $`A`$ | Runtime (Build / Verify) |
|---|---|---|---:|---|---|
| `upper_5665.bin` | $`1133/2000 = 0.5665 < 17/30`$ | `484,138 + 160,545 + 13,519,800 + 3,164,625 = 17,329,108` | **0** | $`14{,}606{,}753{,}031 / 2{,}989{,}767{,}809 \approx 4.885581`$ | 4.6 s / 2.6 s |
| `lower_5650.bin` | $`113/200 = 0.5650`$ | `484,138 + 160,545 + 13,519,800 + 3,164,625 = 17,329,108` | **0** | $`7{,}304{,}352{,}308 / 8{,}763{,}936{,}965 \approx 0.833456`$ | 5.3 s / 2.4 s |

Both certificates additionally pass:
- **37 exact-rational brute-force tree drift checks** (all small trees $`T_3, K_{1,3}, P_4`$, extreme adversarial trees hitting $`D=D_0-1`$, $`D=D_0`$, $`D > D_0`$, $`L=0`$ chains, $`k > K_F`$, $`L > L_F`$, and 30 random friend trees up to $`n = 45`$) with **0 failures**;
- **Negative controls** (stricter exponents `1132/2000` for UPPER and `114/200` for LOWER, and single-entry perturbations to each of `f`, `g`, `hV`, `hC`), all **rejected** (`exit=2`).
