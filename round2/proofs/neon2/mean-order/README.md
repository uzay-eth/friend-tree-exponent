# All-Label Mean Degree-Share Monotonicity in the Friend Tree (Open Question 3 — Mean Shares)

**Pinned claim commit:** `834bb960235f0d0afdc263da8af510f7c452997b`.  
**Review and dependencies:** GPT-family exact-table and full-proof audits; Claude independently reviewed the analytic skeleton. No inherited proof dependencies.

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-8.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-8.md).

**Author:** Neon-2 lead (Gemini 4 Neon).  
**Dependencies:** **NONE** — 100% self-contained from the definition of the friend-tree process (zero Astra dependencies, zero asymptotics, zero floating-point arithmetic).

---

## Headline Theorems

Let $`(T_n)_{n \ge 2}`$ be the discrete-time friend tree started from any finite tree $`T`$ of size $`n \ge 2`$ (and in particular the original edge-seeded tree $`T_2 = \{1, 2\}`$), where at step $`m \to m + 1`$, a uniform vertex $`u \in T_m`$ picks a uniform neighbor $`x \sim u`$ and attaches vertex $`m + 1`$ as a new leaf adjacent to $`x`$. For every fixed label $`v`$, $`D_v(m)/m \in [0, 1]`$ is a bounded nonnegative supermartingale (Lemma 1), so $`Z_v := \lim_{m \to \infty} D_v(m)/m \in [0, 1]`$ exists almost surely and in $`L^1`$.

### Theorem 1 (Strict Conditional Next-Two-Birth Mean-Share Ordering — Every Tree, Every Size $`n \ge 2`$)
For **every** integer $`n \ge 2`$ and **every** finite tree $`T`$ of size $`n`$,
```math
\mathbb E\big[Z_{n+1} - Z_{n+2} \mid T_n = T\big] \;\ge\; \frac{\delta_*}{n(n+1)(n+2)} \;>\; 0,
\qquad\text{(1)}
```
with the exact certified rational constant
```math
\delta_* \;=\; \frac{2493671}{62500000} \;=\; 0.039898736 \;>\; \frac{1}{26}.
```

### Theorem 2 (All-Label Mean Degree-Share Monotonicity)
For the original edge-seeded friend tree ($`T_2 = \{1, 2\}`$), writing $`m_u := \mathbb E[Z_u]`$,
```math
m_1 \;=\; m_2 \;>\; m_3 \;>\; m_4 \;>\; m_5 \;>\; \cdots,
\qquad\text{(2)}
```
with explicit certified lower bounds
```math
m_2 - m_3 \;\ge\; \frac{409872349}{3000000000} \;\approx\; 0.136624 \;>\; \frac{2}{15}, \qquad
m_u - m_{u+1} \;\ge\; \frac{2493671}{62500000\,(u-1)u(u+1)} \qquad (\forall\, u \ge 3).
```

---

## Complete Proof

### Lemma 1 (Yule-Process Embedding & Scale-Invariant Local Expectation $`M(T', v)`$)
Let $`(T_t)_{t \ge 0}`$ be the continuous-time friend tree started from any finite tree $`T'`$ of size $`N \ge 2`$ (each vertex $`u \in T_t`$ carries an independent rate-1 Poisson clock; on a ring, $`u`$ picks a uniform neighbor $`x \sim u`$ and attaches a new leaf to $`x`$). Then:
1. For every $`v \in T'`$, the discrete degree share $`D_v(m)/m`$ ($`m \ge N`$) is a bounded nonnegative supermartingale, hence converges a.s. and in $`L^1`$ to $`Z_v \in [0, 1]`$.
2. Defining the continuous-time local expected amplitude
   $`\displaystyle M(T', v) \;:=\; \lim_{t \to \infty} e^{-t}\,\mathbb E_{T'}[D_v(t)], \qquad\text{(3)}`$
   we have the exact scale-invariant identity
   $`\displaystyle M(T', v) \;=\; N\,\mathbb E_{T'}[Z_v] \;\ge\; 0, \qquad\text{(4)}`$
   and $`e^{-t} M(T_t, v) = \mathbb E_{T'}[W_N Z_v \mid \mathcal F_t]`$ is a closed (uniformly integrable) Doob martingale on $`[0, \infty]`$ with terminal value $`W_N Z_v`$ ($`W_N = \lim_{t\to\infty} e^{-t} |T_t| \sim \mathrm{Gamma}(N, 1)`$).
3. For every $`T'`$ and every vertex $`v \in T'`$ with degree $`D_v \ge 1`$ and leaf-neighbor count $`L_v \ge 0`$,
   $`\displaystyle (L_v - 1)_+ \;\le\; M(T', v) \;\le\; D_v. \qquad\text{(5)}`$

**Proof.**
At discrete size $`m`$, vertex $`v`$ gains a leaf with probability $`r_{T_m}(v)/m`$, where $`r_T(v) := \sum_{u \sim v} 1/D_T(u) \le D_T(v)`$. Hence $`\mathbb E[D_v(m+1) \mid T_m] = D_v(m) + r_{T_m}(v)/m \le \frac{m+1}{m} D_v(m)`$, so $`D_v(m)/m \in [0, 1]`$ is a bounded supermartingale and converges a.s. and in $`L^1`$ to $`Z_v`$.

In continuous time, the total jump rate at size $`m`$ is $`m`$, so $`N(t) := |T_t|`$ is a rate-1 Yule process started at $`N(0) = N`$, and the embedded discrete sequence $`(T_{\tau_m})_{m \ge N}`$ is identical in law to the discrete-time friend tree and independent of the jump times $`(\tau_m)_{m \ge N}`$. Standard Yule theory gives $`e^{-t} N(t) \to W_N \sim \mathrm{Gamma}(N, 1)`$ a.s. and in $`L^1`$ (in fact $`L^2`$), with $`\mathbb E[W_N] = N`$ and $`W_N \perp Z_v`$. Since $`0 \le e^{-t} D_v(t) = (D_v(t)/N(t))\cdot e^{-t} N(t) \le e^{-t} N(t)`$ and $`e^{-t} N(t) \to W_N`$ in $`L^1`$, generalized dominated convergence yields
```math
M(T', v) = \lim_{t\to\infty} e^{-t}\mathbb E_{T'}[D_v(t)] = \mathbb E_{T'}[W_N Z_v] = \mathbb E[W_N]\,\mathbb E_{T'}[Z_v] = N\,\mathbb E_{T'}[Z_v].
```
Moreover, by the Markov property at time $`t`$, $`\mathbb E_{T'}[W_N Z_v \mid \mathcal F_t] = \lim_{s\to\infty} e^{-s}\mathbb E_{T'}[D_v(s) \mid \mathcal F_t] = e^{-t} M(T_t, v)`$, so $`e^{-t} M(T_t, v)`$ is a closed Doob martingale.

For (5): the generator $`\mathcal L`$ satisfies $`\mathcal L(D_v) = r_T(v) \le D_v`$, so $`e^{-t} D_v(t)`$ is a nonnegative supermartingale (justified by stopping at the auxiliary size cutoff $`\tau_{K'} = \inf\{t : N(t) \ge K'\}`$ and passing $`K' \to \infty`$ by dominated convergence using the single integrable dominating variable $`N(t)`$ at each fixed $`t`$), giving $`M(T', v) \le D_v`$. For the lower bound, let $`g(v) := (L_v - 1)_+`$. When $`L_v = 0`$, $`\mathcal L g = 0 = g`$; when $`L_v = 1`$, $`g = 0`$ and leaf gain moves $`L_v \to 2`$ while promotion moves $`L_v \to 0`$ (change $`0`$), so $`\mathcal L g \ge 0 = g`$; when $`L_v \ge 2`$, leaf gain occurs at rate $`r_T(v) \ge L_v`$ and leaf promotion occurs at rate $`L_v / D_v \le 1`$, so $`\mathcal L g = r_T(v) - L_v / D_v \ge L_v - 1 = g`$. Thus $`e^{-t}(L_v(t) - 1)_+`$ is a submartingale (justified by the same size localization and domination by $`e^{-t} N(t) \in L^1`$) with $`(L_v(t) - 1)_+ \le D_v(t)`$, yielding $`(L_v - 1)_+ \le M(T', v)`$. $`\square`$

---

### Lemma 2 (Signed Two-Birth Identity & Endpoint-Split Reduction)
Fix any finite tree $`T`$ of size $`n \ge 2`$. For each $`x \in T`$, let $`T_{x-a-b}`$ be the tree of size $`n+2`$ obtained by attaching $`a`$ to $`x`$ and $`b`$ to $`a`$, and write $`\Delta_{\mathrm{path}}(T, x) := M(T_{x-a-b}, a) - M(T_{x-a-b}, b)`$. For each edge $`\{x, y\} \in E(T)`$, let $`T_{x, y}`$ be the tree of size $`n+2`$ obtained by attaching a leaf $`a_x`$ to $`x`$ and a leaf $`a_y`$ to $`y`$. Then:
```math
n(n+1)(n+2)\,\mathbb E[Z_{n+1} - Z_{n+2} \mid T_n = T]
\;=\; \sum_{x \in T} \frac{r_T(x)}{D_x + 1}\,\Delta_{\mathrm{path}}(T, x)
\;+\; \sum_{\{x, y\} \in E(T)} \left(\frac{r_T(y)}{D_y(D_y + 1)} - \frac{r_T(x)}{D_x(D_x + 1)}\right)\!\big(M(T_{x, y}, a_x) - M(T_{x, y}, a_y)\big).
\qquad\text{(6)}
```
Consequently,
```math
n(n+1)(n+2)\,\mathbb E[Z_{n+1} - Z_{n+2} \mid T_n = T]
\;\ge\; \sum_{x \in T} \frac{r_T(x)}{D_x + 1}\!\left[\Delta_{\mathrm{path}}(T, x) \;-\; \frac{1}{D_x}\sum_{y \sim x} M(T_{x, y}, a_x)\right],
\qquad\text{(7)}
```
and furthermore $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge 1`$.

**Proof.**
Write $`a = n + 1`$ and $`b = n + 2`$. At step $`n \to n + 1`$, $`a`$ attaches to $`x \in T`$ with probability $`r_T(x)/n`$, yielding $`T' = T \cup \{x - a\}`$ where $`D_{T'}(x) = D_x + 1`$, $`D_{T'}(a) = 1`$, and $`D_{T'}(z) = D_z`$ for all $`z \in T \setminus \{x\}`$. At step $`n + 1 \to n + 2`$, conditioning on $`T_{n+2}`$ and applying Lemma 1 eq. (4) ($`(n+2)\mathbb E[Z_a - Z_b \mid T_{n+2}] = M(T_{n+2}, a) - M(T_{n+2}, b)`$), four disjoint cases occur:
1. $`b`$ attaches to $`a`$ (only possible when $`x`$ rings and selects $`a`$, probability $`\frac{1}{(n+1)(D_x+1)}`$): contributes $`\frac{r_T(x)}{n(n+1)(D_x+1)} \Delta_{\mathrm{path}}(T, x)`$.
2. $`b`$ attaches to $`x`$: $`a`$ and $`b`$ are twin leaves of $`x`$, so $`M(T_{n+2}, a) - M(T_{n+2}, b) = 0`$ by graph symmetry.
3. $`b`$ attaches to $`y \in T \setminus \{x\}`$ with $`y \not\sim x`$: $`r_{T'}(y) = r_T(y)`$, so $`P(a \to x, b \to y) = \frac{r_T(x)r_T(y)}{n(n+1)} = P(a \to y, b \to x)`$ and the opposite assignments cancel identically.
4. $`b`$ attaches to $`y \sim x`$ in $`T`$: incrementing $`D_x \to D_x + 1`$ changes $`r_{T'}(y)`$ to $`r_T(y) - \frac{1}{D_x} + \frac{1}{D_x+1} = r_T(y) - \frac{1}{D_x(D_x+1)}`$ (while leaf $`a`$ is not adjacent to $`y`$). Pairing $`(a \to x, b \to y)`$ with $`(a \to y, b \to x)`$ gives the exact net signed probability
   $`\displaystyle \frac{r_T(x)\big(r_T(y) - \frac{1}{D_x(D_x+1)}\big) - r_T(y)\big(r_T(x) - \frac{1}{D_y(D_y+1)}\big)}{n(n+1)} \;=\; \frac{1}{n(n+1)}\left[\frac{r_T(y)}{D_y(D_y+1)} - \frac{r_T(x)}{D_x(D_x+1)}\right],`$
   proving (6).

For (7): since $`M \ge 0`$, expanding each unordered edge term in (6) and dropping the two nonnegative cross terms $`\frac{r_T(y)}{D_y(D_y+1)} M(T_{x,y}, a_x) \ge 0`$ and $`\frac{r_T(x)}{D_x(D_x+1)} M(T_{x,y}, a_y) \ge 0`$ yields
```math
\left(\frac{r_T(y)}{D_y(D_y+1)} - \frac{r_T(x)}{D_x(D_x+1)}\right)\!\big(M(T_{x,y}, a_x) - M(T_{x,y}, a_y)\big)
\;\ge\; -\frac{r_T(x)}{D_x(D_x+1)} M(T_{x,y}, a_x) \;-\; \frac{r_T(y)}{D_y(D_y+1)} M(T_{x,y}, a_y).
```
Regrouping over all vertices $`x \in T`$ (each $`x`$ has $`D_x`$ incident edges $`\{x, y\}`$) yields (7). Finally, in any tree $`T`$ of size $`n \ge 2`$, every vertex satisfies $`D_x \le n - 1`$, hence $`D_x + 1 \le n`$, and $`\sum_{x \in T} r_T(x) = \sum_{u \in T} \sum_{x \sim u} 1/D_u = n`$, so $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge \frac{1}{n}\sum_{x \in T} r_T(x) = 1`$. $`\square`$

---

### Lemma 3 (Six Exact Finite-State Local Martingales — Zero-Float Verification)
Fix integer cutoffs $`D_M = 40`$, $`K = 14`$, $`D_{x,\max} = 8`$, and scale $`S = 10^9`$. The standalone pure-`Fraction` script `cert_mean_order.py` constructs and independently verifies six integer-tick tables (divided by $`S`$) satisfying exact boundary and generator-drift inequalities at all extreme external rates:

1. **1-vertex core bounds $`v(D, L)`$ and $`w(D, L)`$ ($`2 \le D \le D_M, 0 \le L \le D`$):**
   For any core vertex $`u`$ in any tree $`T'`$, its degree $`D_u`$ and leaf count $`L_u`$ jump only via leaf gain $`(D, L) \to (D+1, L+1)`$ at rate $`r_u \in [L, L + (D - L)/2]`$ (since each of the $`D - L`$ core neighbors has degree $`\ge 2`$) and leaf promotion $`(D, L) \to (D, L-1)`$ at rate $`L/D`$. On $`\tau := \inf\{t \ge 0 : D_u(t) = D_M + 1\}`$, on $`\{\tau = \infty\}`$ we have $`D_u(t) \le D_M`$ for all $`t`$, hence $`Z_u = 0`$ a.s. and $`e^{-t} v, e^{-t} w \to 0`$. At $`D = D_M + 1`$, $`v(D_M+1, L) \le (L-1)_+ \le M(T_\tau, u) \le D_M + 1 \le w(D_M+1, L)`$ by Lemma 1 eq. (5). Since $`e^{-t} M(T_t, u)`$ is a closed Doob martingale and $`e^{-t} v(D_u(t), L_u(t))`$ (resp. $`e^{-t} w(D_u(t), L_u(t))`$) is a bounded submartingale (resp. supermartingale) on $`[0, \tau]`$, optional stopping gives
   $`\displaystyle v(D_u, L_u) \;\le\; M(T', u) \;\le\; w(D_u, L_u) \qquad (\forall\, 2 \le D_u \le D_M,\; 0 \le L_u \le D_u).`$
   In particular, $`v(2, 2) = \frac{1187553139}{1000000000}`$ and $`w(2, 1) = \frac{639419937}{500000000} \approx 1.27884`$.

2. **Leaf upper bound $`U(d)`$ ($`2 \le d \le D_M + 1`$):**
   Let $`a`$ be a leaf attached to a core vertex $`x`$ of degree $`d \ge 2`$ in any tree $`T'`$, and let $`\tau_a := \inf\{t : D_a(t) = 2\}`$. Since $`D_x(t-) \ge d`$ for all $`t`$, the promotion time $`\tau_a`$ has predictable intensity $`\mathbf 1_{\{t < \tau_a\}}/D_x(t-) \le \mathbf 1_{\{t < \tau_a\}}/d`$; thus the survival function $`H(t) := \mathbb P(t < \tau_a) = 1 - \int_0^t \mathbb E[\mathbf 1_{\{s < \tau_a\}}/D_x(s-)]\,ds`$ satisfies $`H'(t) \ge -H(t)/d`$ a.e., hence $`H(t) \ge e^{-t/d}`$ and $`\mathbb E[e^{-\tau_a}] = 1 - \int_0^\infty e^{-t}H(t)\,dt \le \frac{1}{d+1}`$. At $`t = \tau_a`$, $`a`$ is a degree-2 core with 1 leaf, so $`M(T_{\tau_a}, a) \le w(2, 1)`$ (while on $`\{\tau_a = \infty\}`$, $`Z_a = 0`$), giving the global bound $`M(T', a) \le \frac{w(2, 1)}{d+1}`$ for all $`d \ge 2`$. Stopping at $`\tau := \tau_a \wedge \inf\{t : D_x(t) = D_M + 1\}`$: while $`t < \tau`$, $`a`$ rings at rate $`1`$ incrementing $`D_x \to D_x + 1`$ (and any other neighbor of $`x`$ increments $`D_x`$ at nonnegative rate $`r_{\mathrm{other}} \ge 0`$), while $`x`$ promotes $`a`$ at exact rate $`1/D_x`$. Because `cert_mean_order.py` verifies $`U(d+1) \le U(d)`$ for all $`2 \le d \le D_M`$, $`U(D_M+1) \ge w(2, 1)/(D_M+2)`$, and $`(U(d+1) - U(d)) + \frac{1}{d}(w(2, 1) - U(d)) - U(d) \le 0`$, optional stopping yields
   $`\displaystyle M(T', a) \;\le\; U(\min(D_{T'}(x),\, D_M + 1)) \qquad (\forall\, T',\; \text{leaf } a \sim x \text{ with } D_{T'}(x) \ge 2),`$
   with $`U(2) = \frac{367808441}{1000000000}`$, $`U(3) = \frac{70025291}{250000000}`$, and $`U(4) = \frac{14205589}{62500000}`$ (and $`U(d) \le U(4)`$ for all $`d \ge 4`$).

3. **Phase-2 adjacent-core difference lower bound $`\phi_2(D_a, L_a; D_b, L_b)`$ ($`2 \le D_a, D_b \le K`$, $`0 \le L_a \le D_a - 1`$, $`0 \le L_b \le D_b - 1`$):**
   For any two adjacent core vertices $`a \sim b`$ in $`T'`$, let $`\tau := \inf\{t : D_a(t) = K + 1 \text{ or } D_b(t) = K + 1\}`$. On $`\{\tau = \infty\}`$, $`Z_a = Z_b = 0`$ a.s.; on $`\{\tau < \infty\}`$, $`\phi_2 \le v(D_a, L_a) - w(D_b, L_b) \le M(T_\tau, a) - M(T_\tau, b)`$. While $`t < \tau`$, no single jump alters both $`(D_a, L_a)`$ and $`(D_b, L_b)`$; $`a`$ gains a leaf at rate $`L_a + 1/D_b + r_{\mathrm{ext}}(a)`$ with $`r_{\mathrm{ext}}(a) \in [0, (D_a - L_a - 1)/2]`$, $`b`$ gains a leaf at rate $`L_b + 1/D_a + r_{\mathrm{ext}}(b)`$ with $`r_{\mathrm{ext}}(b) \in [0, (D_b - L_b - 1)/2]`$, and leaf promotions occur at exact rates $`L_a/D_a`$ and $`L_b/D_b`$. Verifying submartingale drift $`\ge 0`$ at all 4 extreme-rate corners on all 10,816 interior states gives $`\phi_2(D_a, L_a; D_b, L_b) \le M(T', a) - M(T', b)`$.

4. **Phase-1 core–leaf difference lower bound $`\phi_1(D_a, L_a)`$ ($`2 \le D_a \le K, 1 \le L_a \le D_a - 1`$):**
   For any core $`a`$ with at least one core neighbor and a distinguished leaf $`b \sim a`$, let $`\tau := \inf\{t : D_b(t) = 2 \text{ or } D_a(t) = K + 1\}`$. At $`D_a(\tau) = K + 1`$ (with $`b`$ still a leaf), $`M(a) - M(b) \ge v(K+1, L_a) - U(K+1) \ge \phi_1(K+1, L_a)`$; when $`b`$ promotes (exact rate $`1/D_a`$), $`(a, b)`$ enters Phase 2 with $`M(a) - M(b) \ge \phi_2(D_a, L_a - 1; 2, 1)`$; other leaves of $`a`$ promote at rate $`(L_a - 1)/D_a`$, and $`a`$ gains a leaf at rate $`r_a \in [L_a, L_a + (D_a - L_a)/2]`$. Verifying submartingale drift $`\ge 0`$ on all 91 interior states gives $`\Delta_{\mathrm{path}}(T, x) \ge \phi_1(2, 1) = 834963/3125000`$ for **all** $`x \in T`$.

5. **Phase-1 parent-degree-tracked lower bound $`\phi_{1,x}(D_x, D_a, L_a)`$ ($`2 \le D_x \le D_{x,\max}`$, $`2 \le D_a \le K`$, $`1 \le L_a \le D_a - 1`$):**
   Tracking also the core neighbor $`x \sim a`$ up to $`\tau := \inf\{t : D_b(t) = 2 \text{ or } D_a(t) = K + 1 \text{ or } D_x(t) = D_{x,\max} + 1\}`$: $`x`$ sends exact rate $`1/D_x`$ to $`a`$ (so $`a`$'s gain rate lies in $`[L_a + 1/D_x,\; L_a + 1/D_x + (D_a - L_a - 1)/2]`$), while $`x`$'s degree increments at rate $`r_x = 1/D_a + \sum_{u \sim x, u \ne a} 1/D_u \in [1/D_a,\; D_x]`$ (valid even when $`x`$'s other neighbors are leaves!). At $`D_x(\tau) = D_{x,\max} + 1`$, the boundary $`\phi_1(D_a, L_a)`$ holds by part (4). Verifying submartingale drift $`\ge 0`$ on all 637 interior states yields $`\Delta_{\mathrm{path}}(T, x) \ge \phi_{1,x}(D_T(x) + 1, 2, 1)`$ whenever $`D_T(x) + 1 \le D_{x,\max}`$.

---

### Proof of Theorems 1 and 2
In Lemma 2 eq. (7), for every vertex $`x \in T`$ and every neighbor $`y \sim x`$, the tree $`T_{x, y}`$ has $`D_{T_{x,y}}(x) = D_T(x) + 1 \ge 2`$ with $`a_x`$ a leaf of $`x`$, while $`T_{x-a-b}`$ has $`D_{T_{x-a-b}}(x) = D_T(x) + 1 \ge 2`$, $`D_a = 2, L_a = 1`$, and $`b`$ a leaf of $`a`$. Splitting by $`D_T(x)`$ and applying Lemma 3 (`cert_mean_order.py` exact output):
- **$`D_T(x) = 1`$ ($`D'_x = 2`$):** $`\Delta_{\mathrm{path}}(T, x) - M(T_{x,y}, a_x) \ge \phi_{1,x}(2, 2, 1) - U(2) = \frac{54821827}{500000000} \approx 0.109644`$.
- **$`D_T(x) = 2`$ ($`D'_x = 3`$):** $`\Delta_{\mathrm{path}}(T, x) - M(T_{x,y}, a_x) \ge \phi_{1,x}(3, 2, 1) - U(3) = \frac{28840811}{250000000} \approx 0.115363`$.
- **$`D_T(x) \ge 3`$ ($`D'_x \ge 4`$):** $`\Delta_{\mathrm{path}}(T, x) - M(T_{x,y}, a_x) \ge \phi_1(2, 1) - U(4) = \frac{2493671}{62500000} \approx 0.0398987`$.

Thus every bracketed summand in (7) is at least $`\delta_* := \frac{2493671}{62500000} > \frac{1}{26}`$. Since $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge 1`$ (Lemma 2), eq. (7) immediately proves Theorem 1.

For Theorem 2: on the edge seed $`T_2 = \{1, 2\}`$, $`m_1 = m_2`$ by symmetry. At step $`n = 2 \to 3`$, vertex $`3`$ attaches uniformly to vertex $`1`$ or $`2`$, producing the 3-vertex path $`T_3`$ where vertex $`2`$ is the center ($`D = 2, L = 2`$) with probability $`1/2`$ and a leaf with probability $`1/2`$, while vertex $`3`$ is always a leaf ($`D = 1`$, attached to the degree-2 center). Hence by Lemma 1 eq. (4) and Lemma 3,
```math
m_2 - m_3 \;=\; \frac{1}{2}\big(\mathbb E_{T_3}[Z_{\mathrm{center}}] - \mathbb E_{T_3}[Z_{\mathrm{leaf}}]\big)
\;=\; \frac{1}{6}\big(M(T_3, \mathrm{center}) - M(T_3, \mathrm{leaf})\big)
\;\ge\; \frac{v(2, 2) - U(2)}{6}
\;=\; \frac{409872349}{3000000000} \;>\; \frac{2}{15}.
```
Taking expectations of Theorem 1 at $`n = u - 1 \ge 2`$ gives $`m_u - m_{u+1} \ge \frac{\delta_*}{(u-1)u(u+1)} > 0`$ for all $`u \ge 3`$, completing the proof. $`\square`$

---

## Exact Replay Command

```sh
python3 -B round2/proofs/neon2/mean-order/cert_mean_order.py
```

Executes in ~3 seconds using zero floating-point arithmetic (including under `python3 -OB`, with explicit runtime checks and strict integer tick validation), checking all `13,299` rational martingale rows across the 6 tables, all boundary and nonincreasing monotonicity conditions, the exact signed two-birth identity and endpoint-split inequality on all `24` unlabelled trees of sizes $`n = 2..7`$, and `10/10` planted single-tick and non-integer negative controls.
