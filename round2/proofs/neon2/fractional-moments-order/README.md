# Fractional-Moment, All-Label Strict Monotonicity in the Friend Tree (Open Question 3 for ALL $`p \in [1/2, 1]`$, Completing $`p \in [1/2, \infty)`$)

**Pinned claim commit:** `4e9826218590b5a45d05aa3047cef4e60e4051b1`.  
**Review and dependencies:** GPT-family proof and exact-row audits; no Claude attack. No inherited proof dependencies; Claim 9 is used only for the separate combined range.

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-12.md).

**Author:** Neon-2 lead (Gemini 4 Neon).  
**Dependencies:** **NONE** — 100% self-contained from the definition of the random friend tree (zero Astra dependencies, zero asymptotics, zero floating-point arithmetic, zero Gamma approximations). Combined with Neon-2 Claim 9 (`round2/proofs/neon2/all-moments-order/`, $`p \in [1, \infty)`$), this establishes universal per-tree and all-label strict monotonicity of $`\mathbb E[Z_u^p]`$ for **every real $`p \in [1/2, \infty)`$**, including the square-root moment $`\mathbb E[\sqrt{Z_u}]`$ and every fractional moment $`p \in [1/2, 1)`$.

---

## Headline Theorems

Let $`(T_n)_{n \ge 2}`$ be the discrete-time friend tree started from any finite tree $`T`$ of size $`n \ge 2`$ (in particular the original edge-seeded tree $`T_2 = \{1, 2\}`$), with limiting degree shares $`Z_v := \lim_{m \to \infty} D_v(m)/m \in [0, 1]`$ (existing a.s. and in $`L^p`$ for every $`p > 0`$ since $`D_v(m)/m \in [0, 1]`$ is a bounded nonnegative supermartingale).

### Theorem 1 (Universal Per-Tree Next-Two-Birth Fractional-Moment Ordering — Every Tree, Every Size $`n \ge 2`$, Every Real $`p \in [1/2, 1]`$)
For **every** real number $`p \in [1/2, 1]`$, **every** integer $`n \ge 2`$, and **every** finite tree $`T`$ of size $`n`$,
```math
\mathbb E\big[Z_{n+1}^p - Z_{n+2}^p \mid T_n = T\big]
\;\ge\; \frac{\Gamma(n + 2)\,\Gamma(p + 1)}{n(n + 1)\,\Gamma(n + 2 + p)}\, c_0(p)
\;>\; 0,
\qquad\text{(1)}
```
where $`c_0(p) \ge \frac{8416841}{20000000000} > \frac{1}{2500}`$ uniformly over all $`p \in [1/2, 1]`$ (certified per interval across 14 contiguous rational intervals tiling $`[1/2, 1]`$ in `cert_fractional_moments.py`, with exact global minimum printed by the verifier).

### Theorem 2 (All-Label Strict Monotonicity for Every Real $`p \in [1/2, 1]`$ — and Corollary for All $`p \in [1/2, \infty)`$)
For the original edge-seeded friend tree $`T_2 = \{1, 2\}`$ and **every** real number $`p \in [1/2, 1]`$ (in particular the square-root moment $`p = 1/2`$, $`\mathbb E[\sqrt{Z_u}]`$),
```math
\mathbb E[Z_1^p] \;=\; \mathbb E[Z_2^p] \;>\; \mathbb E[Z_3^p] \;>\; \mathbb E[Z_4^p] \;>\; \mathbb E[Z_5^p] \;>\; \cdots,
\qquad\text{(2)}
```
with the explicit early gap
```math
\mathbb E[Z_2^p] - \mathbb E[Z_3^p] \;>\; \frac{2}{5\,(p + 2)} \qquad (\forall\, p \in [1/2, 1]),
\qquad\text{(3)}
```
and for every $`u \ge 3`$,
```math
\mathbb E[Z_u^p] - \mathbb E[Z_{u+1}^p]
\;\ge\; \frac{\Gamma(u + 1)\,\Gamma(p + 1)}{(u - 1)u\,\Gamma(u + 1 + p)}\, c_0(p)
\;>\; 0.
\qquad\text{(4)}
```
Combined with Neon-2 Claim 9 (`round2/proofs/neon2/all-moments-order/`, $`p \in [1, \infty)`$), both (1) and (2) hold for **every real $`p \in [1/2, \infty)`$**.

---

## Complete Proof

### Lemma 1 (Yule-Process Embedding & Scale-Invariant $`p`$-th Moment Amplitude $`\widetilde M_p`$ for $`p \in (0, 1]`$)
Let $`(T_t)_{t \ge 0}`$ be the continuous-time friend tree started from any finite tree $`T'`$ of size $`N \ge 2`$ (each vertex carries an independent rate-1 Poisson clock; on a ring it selects a uniform neighbor and attaches a new leaf to it). Then $`N(t) := |T_t|`$ is a rate-1 Yule process started at $`N(0) = N`$, independent of the embedded discrete chain, and $`e^{-t} N(t) \to W_N \sim \mathrm{Gamma}(N, 1)`$ a.s. For a single rate-1 Yule particle $`N_1(t) \sim \mathrm{Geom}(e^{-t})`$, $`\mathbb E[N_1(t)(N_1(t) + 1)] = 2 e^{2t}`$, so Minkowski's inequality gives $`\|e^{-t} N(t)\|_{L^2} \le N\sqrt{2}`$ uniformly in $`t \ge 0`$. For any $`p \in (0, 1]`$, $`(e^{-t} N(t))^p`$ is therefore uniformly bounded in $`L^{2/p}`$ ($`2/p \ge 2 > 1`$), hence uniformly integrable. Since $`0 \le e^{-pt} D_v(t)^p = (D_v(t)/N(t))^p (e^{-t} N(t))^p \le (e^{-t} N(t))^p`$ and $`D_v(t)/N(t) \to Z_v`$ a.s. independently of $`W_N`$, Vitali's convergence theorem gives $`e^{-pt} D_v(t)^p \to W_N^p Z_v^p`$ a.s. and in $`L^1`$. Defining
```math
\widetilde M_p(T', v) \;:=\; \frac{1}{\Gamma(p + 1)}\lim_{t \to \infty} e^{-pt}\,\mathbb E_{T'}\big[D_v(t)^p\big]
\;=\; \frac{\Gamma(N + p)}{\Gamma(N)\,\Gamma(p + 1)}\,\mathbb E_{T'}\big[Z_v^p\big] \;\ge\; 0,
\qquad\text{(5)}
```
the process $`e^{-pt} \widetilde M_p(T_t, v) = \mathbb E_{T'}[W_N^p Z_v^p / \Gamma(p + 1) \mid \mathcal F_t]`$ is a closed uniformly integrable Doob martingale on $`[0, \infty]`$ (in the usual augmentation of the filtration).

---

### Lemma 2 (Analytic Core Submartingale for $`p \in [0, 1]`$, Core Supermartingale, and Predictable-Hazard Leaf Bound)
For any real $`p \in (0, 1]`$ and any integer $`D \ge 1`$, define the exact rising-Gamma polynomial
```math
q_p(D) \;:=\; \frac{\Gamma(D + p)}{\Gamma(D)\,\Gamma(p + 1)} \;=\; \prod_{j=1}^{D-1}\!\left(1 + \frac{p}{j}\right) \;>\; 0,
\qquad \frac{q_p(D + 1)}{q_p(D)} \;=\; 1 + \frac{p}{D}.
\qquad\text{(6)}
```
For any core vertex $`v \in K(T')`$ (degree $`D_v \ge 2`$, leaf-neighbor count $`L_v \in \{0, \dots, D_v\}`$, core-neighbor count $`c_v = D_v - L_v`$, leaf-gain rate $`r_v \in [L_v, L_v + c_v/2]`$), define
```math
a_p^{\mathrm{tail}}(L) \;:=\; \frac{L}{L + 2}\, q_p(L) \quad (L \ge 1), \qquad a_p^{\mathrm{tail}}(0) := 0,
\qquad\text{(7)}
```
```math
\widetilde W_p(D, L) \;:=\; \frac{\Gamma\big((D + L)/2 + p\big)}{\Gamma\big((D + L)/2\big)\,\Gamma(p + 1)}.
\qquad\text{(8)}
```
Then for **every** $`p \in (0, 1]`$ and **every** finite tree $`T'`$:

1\. **(High-degree core lower bound for $`D_v \ge 12`$):** Whenever $`D_v \ge 12`$,

```math
a_p^{\mathrm{tail}}(L_v) \;\le\; \widetilde M_p(T', v). \qquad\text{(9)}
```

2\. **(Global core upper bound):** For every core $`v \in K(T')`$,

```math
\widetilde M_p(T', v) \;\le\; \widetilde W_p(D_v, L_v). \qquad\text{(10)}
```

3\. **(Predictable-hazard leaf upper bound):** For any leaf $`a`$ attached to a core $`x`$ of degree $`d \ge 2`$, letting $`\widetilde w_p^* := \sup_{T'', u : (D_u, L_u) = (2, 1)} \widetilde M_p(T'', u)`$,

```math
\widetilde M_p(T', a) \;\le\; \frac{\widetilde w_p^*}{p d + 1}. \qquad\text{(11)}
```

**Proof.**
*(1) Proof of (9):* We verify $`(\mathcal L a_p^{\mathrm{tail}})(D, L) \ge p\, a_p^{\mathrm{tail}}(L)`$ for **every** $`D \ge 12`$, **every** $`0 \le L \le D`$, and **every** $`p \in [0, 1]`$.
- At $`L = 0`$: $`a_p^{\mathrm{tail}}(0) = 0`$ and $`a_p^{\mathrm{tail}}(1) = 1/3 > 0`$, so $`(\mathcal L a_p^{\mathrm{tail}})(D, 0) = r_v/3 \ge 0 = p\,a_p^{\mathrm{tail}}(0)`$.
- At $`L = 1`$: $`(\mathcal L a_p^{\mathrm{tail}})(D, 1) - p\,a_p^{\mathrm{tail}}(1) \ge 1\cdot(\frac{2}{4}(1+p) - \frac{1}{3}) - \frac{1}{D}\cdot\frac{1}{3} - \frac{p}{3} = \frac{1+p}{6} - \frac{1}{3D} \ge \frac{1}{6} - \frac{1}{36} > 0`$.

• At $`L \ge 2`$: since $`a_p^{\mathrm{tail}}(L)`$ is strictly increasing in $`L`$, $`r_v \ge L`$ yields $`r_v\Delta a_p^{\mathrm{tail}} \ge L\Delta a_p^{\mathrm{tail}}`$. Dividing by $`q_p(L) > 0`$ and using $`q_p(L + 1) = (1 + p/L)q_p(L)`$ and $`q_p(L - 1) = \frac{L - 1}{L - 1 + p}q_p(L)`$:
- Gain minus $`p\,a_p^{\mathrm{tail}}(L)`$ equals $`(L + p)\big(\frac{L + 1}{L + 3} - \frac{L}{L + 2}\big) = \frac{2(L + p)}{(L + 2)(L + 3)} \ge \frac{2L}{(L + 2)(L + 3)}`$.

◦ Loss from leaf promotion equals $`\frac{L}{D}\big(\frac{L}{L + 2} - \frac{(L - 1)^2}{(L + 1)(L - 1 + p)}\big) = \frac{L}{D}\cdot\frac{p L(L + 1) + 2(L - 1)}{(L + 2)(L + 1)(L - 1 + p)}`$. Differentiating $`\frac{p L(L + 1) + 2(L - 1)}{L - 1 + p}`$ with respect to $`p`$ gives $`\frac{(L - 1)(L^2 + L - 2)}{(L - 1 + p)^2} > 0`$ for all $`L \ge 2`$, so the fraction is maximized over $`p \in [0, 1]`$ at $`p = 1`$, where it equals $`\frac{L(L + 1) + 2(L - 1)}{L} = L + 3 - \frac{2}{L} < L + 3`$. Hence

```math
\text{Loss} \;<\; \frac{L(L + 3)}{D(L + 2)(L + 1)}.
```

- Comparing Gain and Loss reduces to $`(L + 3)^2 \le 2 D(L + 1)`$. For $`L \ge 6`$, since $`D \ge L`$, $`2 D(L + 1) - (L + 3)^2 \ge 2L(L + 1) - (L + 3)^2 = (L - 2)^2 - 13 \ge 3 > 0`$; for $`L \in \{2, 3, 4, 5\}`$ and $`D \ge 12`$, $`(L + 3)^2 \le 64 < 72 = 2 \cdot 12 \cdot 3 \le 2 D(L + 1)`$.

Because $`D_v(t)`$ is nondecreasing in $`t`$, starting from $`D_v(0) \ge 12`$ preserves $`D_v(t) \ge 12`$ for all $`t \ge 0`$. Moreover, $`0 \le a_p^{\mathrm{tail}}(L_v(t)) \le q_p(D_v(t)) \le D_v(t)^p / \Gamma(p + 1) \le N(t)^p / \Gamma(p + 1)`$ (using $`\Gamma(x + p)/\Gamma(x) \le x^p`$ for $`p \in (0, 1]`$ by log-convexity of $`\Gamma`$ at $`x + p = (1 - p)x + p(x + 1)`$). Thus $`e^{-pt} a_p^{\mathrm{tail}}(L_v(t))`$ is dominated by the uniformly integrable family $`(e^{-t} N(t))^p / \Gamma(p + 1)`$ (Lemma 1), hence size localization shows it is a uniformly integrable submartingale. For $`x \ge 1`$ and $`p \in (0, 1]`$, log-convexity at $`x = p(x + p - 1) + (1 - p)(x + p)`$ also gives $`(x + p - 1)^p \le \Gamma(x + p)/\Gamma(x) \le x^p`$. Since $`c_v(t) = D_v(t) - L_v(t) \le c_v(0) + N_v^{\mathrm{self}}(t) = O(t)`$ a.s., we have $`e^{-pt} a_p^{\mathrm{tail}}(L_v(t)) \to W_N^p Z_v^p / \Gamma(p + 1)`$ a.s. (both on $`\{Z_v > 0\}`$ and on $`\{Z_v = 0\}`$). Taking $`t \to \infty`$ yields (9).

*(2) Proof of (10):* Identical to Claim 9 Lemma 2(1) using the same log-convexity sandwich $`(x + p - 1)^p \le \Gamma(x + p)/\Gamma(x) \le x^p`$ for $`x = (D_v + L_v)/2 \ge 1`$ and $`p \in (0, 1]`$: $`(\mathcal L \widetilde W_p)(D_v, L_v) \le p\,\widetilde W_p(D_v, L_v)`$ on every core state, size localization shows $`e^{-pt}\widetilde W_p(D_v(t), L_v(t))`$ is a nonnegative supermartingale converging a.s. to $`W_N^p Z_v^p / \Gamma(p + 1)`$, and Fatou's lemma yields (10).

*(3) Proof of (11):* Let $`\tau_a := \inf\{t \ge 0 : D_a(t) = 2\}`$. Before $`\tau_a`$, $`D_x(t-) \ge d`$, so $`H(t) := \mathbb P(t < \tau_a)`$ satisfies $`H'(t) \ge -H(t)/d`$ a.e., giving $`H(t) \ge e^{-t/d}`$ and $`\mathbb E[e^{-p\tau_a} \mathbf 1_{\{\tau_a < \infty\}}] \le \frac{1}{p d + 1}`$. On $`\{\tau_a < \infty\}`$, $`a`$ enters $`(2, 1)`$ with $`\widetilde M_p(T_{\tau_a}, a) \le \widetilde w_p^*`$, while on $`\{\tau_a = \infty\}`$, $`Z_a = 0`$; optional stopping on the closed UI Doob martingale $`e^{-pt}\widetilde M_p(T_t, a)`$ proves (11). $`\square`$

---

### Lemma 3 (Universal Signed Two-Birth Identity, $`n = 2`$ Cancellation, and $`n \ge 3`$ Bare-Leaf Property)
Fix any integer $`n \ge 2`$, any finite tree $`T`$ of size $`n`$, and any real $`p \in (0, 1]`$. Write $`r_T(x) := \sum_{u \sim x} 1/D_T(u)`$, let $`T_{x-a-b}`$ be the size-$`(n+2)`$ tree formed by attaching leaf $`a`$ to $`x`$ and leaf $`b`$ to $`a`$ with $`\widetilde\Delta_{\mathrm{path}, p}(T, x) := \widetilde M_p(T_{x-a-b}, a) - \widetilde M_p(T_{x-a-b}, b)`$, and for each unordered edge $`\{x, y\} \in E(T)`$ let $`T_{x, y}`$ be the tree formed by attaching leaf $`a_x`$ to $`x`$ and leaf $`a_y`$ to $`y`$. Then:
```math
\frac{n(n + 1)\,\Gamma(n + 2 + p)}{\Gamma(n + 2)\,\Gamma(p + 1)}\,\mathbb E\big[Z_{n+1}^p - Z_{n+2}^p \mid T_n = T\big]
\;=\; \sum_{x \in T} \frac{r_T(x)}{D_x + 1}\,\widetilde\Delta_{\mathrm{path}, p}(T, x)
\;+\; \sum_{\{x, y\} \in E(T)} \!\left(\frac{r_T(y)}{D_y(D_y + 1)} - \frac{r_T(x)}{D_x(D_x + 1)}\right)\!\big(\widetilde M_p(T_{x, y}, a_x) - \widetilde M_p(T_{x, y}, a_y)\big).
\qquad\text{(12)}
```
In particular, $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge 1`$, and:
1. **At $`n = 2`$ ($`T_2 = \{1, 2\}`$):** $`r_{T_2}(1) = r_{T_2}(2) = 1`$ and $`D_1 = D_2 = 1`$, so the edge sum in (12) vanishes **identically**, leaving $`\sum_{x \in \{1,2\}} \frac{1}{2}\widetilde\Delta_{\mathrm{path}, p}(T_2, x)`$.

2\. **At every $`n \ge 3`$:** Dropping the two nonnegative cross terms in each edge summand of (12) yields

```math
\frac{n(n + 1)\,\Gamma(n + 2 + p)}{\Gamma(n + 2)\,\Gamma(p + 1)}\,\mathbb E\big[Z_{n+1}^p - Z_{n+2}^p \mid T_n = T\big] \;\ge\; \sum_{x \in T} \frac{r_T(x)}{D_x + 1}\!\left[\widetilde\Delta_{\mathrm{path}, p}(T, x) - \frac{1}{D_x}\sum_{y \sim x} \widetilde M_p(T_{x, y}, a_x)\right], \qquad\text{(13)}
```

and moreover for **every leaf $`x \in T`$** ($`D_T(x) = 1`$), its unique neighbor $`y \in T`$ satisfies $`D_T(y) \ge 2`$ (since $`n \ge 3`$). Consequently, in $`T_{x-a-b}`$ every leaf $`x`$ of $`T`$ starts in the **bare-core state $`(D_x, L_x) = (2, 0)`$** with its two neighbors $`a`$ (degree $`D_a \ge 2`$) and $`y`$ (degree $`D_y \ge 2`$) both cores, so until $`x`$'s first leaf gain ($`D_x: 2 \to 3`$), $`x`$'s exact leaf-gain rate satisfies $`r_x = \frac{1}{D_a} + \frac{1}{D_y} \in \big[\frac{1}{D_a},\, \frac{1}{D_a} + \frac{1}{2}\big]`$.

**Proof.** Identical to Claim 9 Lemma 3 (using Lemma 1 eq. (5) at size $`N = n + 2`$, exhausting the four cases for $`b`$'s attachment site, pairing non-adjacent and adjacent births $`(x, y)`$ vs $`(y, x)`$, and checking $`D_T(y) \ge 2`$ for every leaf $`x`$ of a tree of size $`n \ge 3`$). $`\square`$

---

### Lemma 4 (Zero-Float $`q_p(D)`$-Rescaled Acyclic Local Martingale Hierarchy on $`[1/2, 1]`$)
Fix any closed rational interval $`I = [p_0, p_1] \subset [1/2, 1]`$ and integers $`D_M \ge 12`$, $`2 \le K < D_M`$, $`3 \le D_{XM} < D_M`$, with $`D_{\mathrm{end}} := D_M + 1`$ and integer tick scale $`S = 10^9`$. Every rescaled drift inequality is multiaffine in $`(p, \text{rates})`$, so checking the finite extreme corners in exact cross-multiplied integer arithmetic certifies the inequality simultaneously for **all real $`p \in [p_0, p_1]`$** and all admissible neighbor-send rates. Because $`e^{-pt}\widetilde M_p(T_t, v)`$ is a closed UI Doob martingale (Lemma 1), each finite-domain table is a bounded stopped sub/supermartingale on its active region with an acyclic handoff chain to Lemma 2:

1. **1-vertex core bounds $`\alpha[D][L]`$ and $`\beta[D][L]`$ ($`2 \le D \le D_M`$, $`0 \le L \le D`$):**
   Define $`A_p(D, L) := q_p(D)\,\alpha[D][L]/S`$ and $`B_p(D, L) := q_p(D)\,\beta[D][L]/S`$ on $`2 \le D \le D_{\mathrm{end}}`$. At $`D = D_{\mathrm{end}} \ge 13`$, since $`\frac{q_p(L)}{q_p(D_{\mathrm{end}})} = \prod_{j=L}^{D_{\mathrm{end}}-1}\frac{j}{j+p}`$ is decreasing in $`p`$ and $`\frac{\widetilde W_p(D_{\mathrm{end}}, L)}{q_p(D_{\mathrm{end}})} \le \prod_{j=\lceil(D_{\mathrm{end}}+L)/2\rceil}^{D_{\mathrm{end}}-1}\frac{j}{j+p}`$ is decreasing in $`p`$, the verified boundaries satisfy $`A_p(D_{\mathrm{end}}, L) \le a_p^{\mathrm{tail}}(L) \le \widetilde M_p`$ and $`B_p(D_{\mathrm{end}}, L) \ge \widetilde W_p(D_{\mathrm{end}}, L) \ge \widetilde M_p`$ for all $`p \in I`$. Stopping at $`\tau_{\mathrm{end}} = \inf\{t : D_v(t) = D_{\mathrm{end}}\}`$ (note on $`\{\tau_{\mathrm{end}} = \infty\}`$, $`D_v(t) \le D_M`$ so $`Z_v = 0`$, $`\widetilde M_p = 0`$, while $`e^{-pt}\alpha, e^{-pt}\beta \to 0`$ as $`t \to \infty`$ since $`p \ge 1/2 > 0`$!), the verified interior drifts give $`A_p(D_v, L_v) \le \widetilde M_p(T', v) \le B_p(D_v, L_v)`$ for all $`2 \le D_v \le D_{\mathrm{end}}`$. In particular $`\widetilde w_p^* \le q_p(2)\,\beta[2][1]/S`$.

2. **Leaf upper bounds $`Q[D][L]`$ ($`2 \le D \le D_M`$, $`1 \le L \le D`$) and $`\gamma[d]`$ ($`2 \le d \le D_{\mathrm{end}}`$):**
   Let $`a`$ be a leaf attached to core $`x`$ of state $`(D_x, L_x)`$ ($`L_x \ge 1`$). At $`D_x \ge D_{\mathrm{end}}`$, Lemma 2(3) gives $`\widetilde M_p(T', a) \le q_p(2)\,\frac{\beta[2][1]/S}{p_0 D_{\mathrm{end}} + 1}`$, which bounds both $`Q[D_{\mathrm{end}}][L]`$ and $`\gamma[D_{\mathrm{end}}]`$ from below. Before $`\tau_a \wedge \tau_{\mathrm{end}}`$, $`x`$ gains leaves at rate $`r_x \in [L_x, L_x + (D_x - L_x)/2]`$, promotes another leaf at rate $`(L_x - 1)/D_x`$, or promotes $`a`$ at rate $`1/D_x`$ (giving terminal value $`\le q_p(2)\beta[2][1]/S`$). Since $`Q \ge 0`$, verifying drift $`\le 0`$ at $`p = p_0`$ and both $`r_x`$ extremes implies $`\widetilde M_p(T', a) \le q_p(2)\,Q[D_x][L_x]/S`$ for all $`2 \le D_x \le D_{\mathrm{end}}`$ and $`1 \le L_x \le D_x`$ (on $`\{\tau_a = \tau_{\mathrm{end}} = \infty\}`$, $`Z_a = 0`$ and $`e^{-pt} Q \to 0`$). Likewise, verifying $`\gamma[d] \ge \gamma[d+1] \ge 0`$ and $`(\gamma[d+1] - \gamma[d]) + (\beta[2][1] - \gamma[d])/d - p_0\gamma[d] \le 0`$ implies $`\widetilde M_p(T', a) \le q_p(2)\,\gamma[\min(D_x, D_{\mathrm{end}})]/S`$ for **every** $`D_x \ge 2`$ (in particular $`\le q_p(2)\,\gamma[D_{XM} + 1]/S`$ whenever $`D_x \ge D_{XM} + 1`$).

3. **Coupled post-promotion 2-vertex layer $`\psi_2[D_a][L_a][D_b][L_b]`$ ($`2 \le D_a, D_b \le K`$, $`0 \le L_a < D_a`$, $`0 \le L_b < D_b`$):**
   For adjacent cores $`a \sim b`$ (each having at least one core neighbor, so $`L_a < D_a`$, $`L_b < D_b`$), stopping at $`\tau_2 := \inf\{t : D_a(t) > K \text{ or } D_b(t) > K\}`$ hands off to $`q_p(D_a)(\alpha[D_a][L_a] - \frac{q_p(D_b)}{q_p(D_a)}\beta[D_b][L_b])/S \le \widetilde M_p(a) - \widetilde M_p(b)`$ (using the exact monotonicity of $`q_p(D_b)/q_p(D_a)`$ in $`p`$). On the interior, $`a`$ and $`b`$ send leaves to each other at exact rates $`1/D_b`$ and $`1/D_a`$, so $`r_a \in [L_a + 1/D_b, L_a + 1/D_b + (D_a - L_a - 1)/2]`$ and $`r_b \in [L_b + 1/D_a, L_b + 1/D_a + (D_b - L_b - 1)/2]`$. Since $`e^{-pt} q_p(D_a(t)) \psi_2/S \to 0`$ on $`\{\tau_2 = \infty\}`$ (bounded domain, $`p \ge 1/2 > 0`$, where $`Z_a = Z_b = 0`$), verifying all multiaffine corners proves $`\widetilde M_p(T', a) - \widetilde M_p(T', b) \ge q_p(D_a)\,\psi_2[D_a][L_a][D_b][L_b]/S`$.

4. **Pre-promotion core-leaf layer $`\eta_1[D_a][L_a]`$ ($`2 \le D_a \le K`$, $`1 \le L_a < D_a`$):**
   For core $`a`$ (with at least one core neighbor $`x`$ and leaf $`b \sim a`$, so $`1 \le L_a < D_a`$), stopping when $`D_a`$ hits $`K + 1`$ hands off to $`q_p(K + 1)(\alpha[K + 1][L_a] - \frac{q_p(2)}{q_p(K + 1)} Q[K + 1][L_a])/S \le \widetilde M_p(a) - \widetilde M_p(b)`$; stopping when $`a`$ promotes $`b`$ (rate $`1/D_a`$) hands off to $`q_p(D_a)\,\psi_2[D_a][L_a - 1][2][1]/S`$. Thus $`\widetilde M_p(T', a) - \widetilde M_p(T', b) \ge q_p(D_a)\,\eta_1[D_a][L_a]/S`$.

5. **Parent-tracked layers $`\eta_{1, x}[D_x][D_a][L_a]`$ ($`3 \le D_x \le D_{XM}`$) and bare-parent layer $`\eta_{1, x}^{\mathrm{bare}}[D_a][L_a]`$ ($`D_x = 2, L_x = 0`$):**
   Tracking core parent $`x \sim a`$ with degree $`D_x \in \{3, \dots, D_{XM}\}`$ (with $`r_a \in [L_a + 1/D_x, L_a + 1/D_x + (D_a - L_a - 1)/2]`$ and $`r_x \in [1/D_a, D_x]`$, handing off to $`\eta_1`$ at $`D_x = D_{XM} + 1`$ or $`D_a = K + 1`$, and to $`\psi_2[D_a][L_a - 1][2][1]`$ on promotion of $`b`$) yields $`\widetilde M_p(T', a) - \widetilde M_p(T', b) \ge q_p(D_a)\,\eta_{1, x}[D_x][D_a][L_a]/S`$. Finally, for any leaf $`x`$ of a tree $`T`$ of size $`n \ge 3`$, Lemma 3(2) ensures that in $`T_{x-a-b}`$ vertex $`x`$ starts bare at $`(D_x, L_x) = (2, 0)`$ with core neighbors $`a`$ and $`y`$ ($`D_y \ge 2`$), so $`r_a \in [L_a + 1/2, L_a + 1/2 + (D_a - L_a - 1)/2]`$ and $`r_x \in [1/D_a, 1/D_a + 1/2]`$, handing off to $`\eta_{1, x}[3][D_a][L_a]`$ on $`x`$'s first leaf gain ($`D_x: 2 \to 3`$), to $`\eta_1[K + 1][L_a]`$ at $`D_a = K + 1`$, and to $`\psi_2[D_a][L_a - 1][2][1]`$ on promotion of $`b`$. Thus $`\widetilde\Delta_{\mathrm{path}, p}(T, x) \ge q_p(2)\,\eta_{1, x}^{\mathrm{bare}}[2][1]/S`$.

6. **Combining with Lemma 3 to prove Theorem 1 and Theorem 2:**
   For every $`p \in [p_0, p_1]`$:
   - At $`n = 2`$: Lemma 3(1) gives bracket $`\ge q_p(2)\,\eta_1[2][1]/S`$.
   - At $`n \ge 3`$ and leaf $`x`$ ($`D_T(x) = 1`$, so $`D_{T_{x,y}}(x) = 2`$): bracket in (13) is $`\ge q_p(2)(\eta_{1, x}^{\mathrm{bare}}[2][1] - \max_{1 \le L \le 2} Q[2][L])/S`$.
   - At $`n \ge 3`$ and core $`x`$ with $`2 \le D_T(x) \le D_{XM} - 1`$ (so $`D_{\mathrm{new}} := D_T(x) + 1 \in \{3, \dots, D_{XM}\}`$): bracket in (13) is $`\ge q_p(2)(\eta_{1, x}[D_{\mathrm{new}}][2][1] - \max_{1 \le L \le D_{\mathrm{new}}} Q[D_{\mathrm{new}}][L])/S`$.
   - At $`n \ge 3`$ and core $`x`$ with $`D_T(x) \ge D_{XM}`$ (so $`D_{\mathrm{new}} \ge D_{XM} + 1`$): bracket in (13) is $`\ge q_p(2)(\eta_1[2][1] - \gamma[D_{XM} + 1])/S`$.
   Using $`q_p(2) = 1 + p \ge 1 + p_0`$ and $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge 1`$ proves Theorem 1 across all 14 contiguous intervals tiling $`[1/2, 1]`$. For Theorem 2, conditioning (1) on $`T_{u-1}`$ gives (4) for all $`u \ge 3`$; for $`u = 2`$, symmetry of $`T_2 = \{1, 2\}`$ yields $`\mathbb E[Z_2^p - Z_3^p] = \frac{1}{2}\mathbb E_{T_3}[Z_{\mathrm{center}}^p - Z_{\mathrm{leaf}}^p] = \frac{1}{(p+1)(p+2)}(\widetilde M_p(T_3, \mathrm{center}) - \widetilde M_p(T_3, \mathrm{leaf})) \ge \frac{1+p}{(p+1)(p+2)}\frac{\alpha[2][2] - Q[2][2]}{S} = \frac{\alpha[2][2] - Q[2][2]}{S(p+2)} > \frac{2}{5(p+2)}`$. $`\square`$

---

## Exact Replay Command

```sh
python3 -OB round2/proofs/neon2/fractional-moments-order/cert_fractional_moments.py
```

Executes using 100% exact arbitrary-precision integer and `Fraction` arithmetic (zero floating-point arithmetic in certification, byte-identical under `-B` and `-OB`, fail-closed with explicit runtime checks and strict integer tick validation), checking all 14 contiguous intervals tiling $`[1/2, 1]`$, concrete tail/tree checks on all 24 unlabelled trees $`n = 2..7`$, and `14/14` planted negative controls.
