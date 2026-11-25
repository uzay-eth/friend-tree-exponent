> **Public checkpoint export note (Part 2 of 2 — Sections 4–10):** Continuation of the edited export of `order-neon/PROOF.md` (frozen original SHA-256 `e85e9cebfbde5439fad9246da2fec2c6eb7f014b09248e960c20d016411ce191`; Part 1 — Sections 1–3 — is in [`PROOF.md`](PROOF.md)). Split solely at the Section 3 / Section 4 boundary to stay within GitHub’s per-page browser math-rendering limit, with navigation links added and zero new mathematical material (all cross-family audit verdicts apply to the single frozen original source SHA-256 listed in [`../PROVENANCE.md`](../PROVENANCE.md)). Formatting changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. **Line numbers cited in errata and audits refer to the frozen original source bytes, not line numbers in these reformatted exports.** Read alongside [`REVIEW-Q3.md`](REVIEW-Q3.md) and [`REVIEW-IMPLEMENTATION.md`](REVIEW-IMPLEMENTATION.md).

# Complete Self-Contained Proof of Birth-Order Stochastic Dominance — Part 2 (Sections 4–10)

> **Navigation (Part 2 of 2 — Sections 4–10):** **Back to Part 1 ([`PROOF.md`](PROOF.md)):** [§1 Model & Main Theorems](PROOF.md#1-exact-model-and-main-theorems) · [§2 Signed Identity](PROOF.md#2-exact-pre-placement-signed-identity) · [§3 Endpoint Laws](PROOF.md#3-exact-endpoint-comparison-laws) | **This file:** [§4 All-m Bounds](#4-exact-discrete-time-continuation-bounds-for-all-m-ge-m_0--200-cert_all_mcpp) · [§5 Finite-m Induction](#5-exact-integer-n-backward-induction-for-4-le-m-le-199-cert_small_mcpp) · [§6 Tree Ledger](#6-dimensionless-normalized-all-tree-signed-ledger) · [§7 High-Tau Tail](#7-self-contained-analytic-high-tau-tail-tau-ge-6-all-b-ge-3) · [§8 Founder Seam](#8-self-contained-founder-seam-b--2-z_2-succ_mathrmst-z_3) · [§9 Int128 Arithmetic](#9-exact-outward-rounded-__int128-arithmetic-guarantees) · [§10 Replay & Evidence](#10-replay-commands-and-exact-verification-evidence)

---

## 4. Exact Discrete-Time Continuation Bounds for All $`m \ge m_0 = 200`$ (`cert_all_m.cpp`)

To close all $`m = n+2 \ge m_0 := 200`$ simultaneously without any continuous-time Gamma/$`S_m`$ Chernoff union bound, we work **directly in discrete tree size $`N \ge m`$** on **raw survival probabilities $`\mathbb{P}(Z_v > x \mid T_N)`$**, parameterized by the dimensionless birth-scale threshold

```math
\tau := (m - 1)\log\frac{1}{1 - x} \in [0, \infty) \qquad \Longleftrightarrow \qquad q := 1 - x = e^{-\tau/(m-1)}.
```

Fix a closed rational box $`[\tau_{\mathrm{lo}}, \tau_{\mathrm{hi}}] \subset [0, 6]`$, a uniform time step $`h := 1/128`$, and a finite horizon $`T = M_{\mathrm{steps}} h`$ with grid knots $`t_i = i h`$ ($`i = 0, \dots, M_{\mathrm{steps}}`$). Note that $`m_0 - 1 = 199 > 128 = 1/h`$.

### 4.1 Exact Raw Endpoint Enclosures Uniform in $`m \ge m_0`$

**Lemma 4.1 (Uniform Lower Enclosure of $`H_{N, j}(x)`$).**  
*Let $`t_{\mathrm{lo}}(N) := \log(N/m)`$. For any $`m \ge m_0`$, any $`\tau \le \tau_{\mathrm{hi}}`$, and any integer $`N \ge m`$ with $`t_{\mathrm{lo}}(N) \le t_{\mathrm{hi}}`$, set $`a_+ := \frac{m_0 e^{t_{\mathrm{hi}}} - 1}{m_0 - 1}`$. Then*

```math
H_{N, j}(x) \ge e^{-1} e^{-\tau_{\mathrm{hi}}(a_+ + 1/(m_0 - 1))} \sum_{k=0}^{j-1} \frac{1}{k!} \sum_{i=0}^{j-k-1} \frac{\tau_{\mathrm{hi}}^i}{i!} \prod_{s=0}^{i-1} \left(a_+ - \frac{s}{m_0 - 1}\right)_+.
```

*Proof.* Since $`H_{N, j}(x)`$ is decreasing in both $`x`$ and $`N`$, it suffices to bound $`H_{N, j}(x)`$ at $`\tau = \tau_{\mathrm{hi}}`$ and $`N \le m e^{t_{\mathrm{hi}}}`$. For $`N \le m e^{t_{\mathrm{hi}}}`$ and $`m \ge m_0`$,

```math
\frac{N - 1}{m - 1} \le \frac{m e^{t_{\mathrm{hi}}} - 1}{m - 1} = e^{t_{\mathrm{hi}}} + \frac{e^{t_{\mathrm{hi}}} - 1}{m - 1} \le a_+.
```

Evaluating the exact binomial tail $`B_{N, r}(x) = q^{N-1}\sum_{i=0}^{r-1}\binom{N-1}{i}(x/q)^i`$ at the single upper envelope $`N_* := 1 + (m-1)a_+ \ge N`$ (or noting $`B_{N, r}(x)`$ is decreasing in $`N`$, and for $`N = \lfloor m e^{t_{\mathrm{hi}}} \rfloor \le 1 + (m-1)a_+`$, we have $`N \le 1 + (m-1)a_+`$ while $`N + 1 > m e^{t_{\mathrm{hi}}}`$? More directly: $`B_{N, r}(x)`$ is decreasing in real $`N`$, so evaluate at real $`N_1 := \lceil 1 + (m-1)a_+ \rceil \ge N`$, which satisfies $`(N_1 - 1)/(m-1) \in [a_+, a_+ + \frac{1}{m-1}] \subset [a_+, a_+ + \frac{1}{m_0-1}]`$!). At this integer $`N_1 \ge N`$:
1. $`q^{N_1 - 1} = \exp\big(-\tau_{\mathrm{hi}}\frac{N_1 - 1}{m - 1}\big) \ge \exp\big(-\tau_{\mathrm{hi}}(a_+ + \frac{1}{m_0 - 1})\big)`$;
2. $`\frac{x}{q} = e^{\tau_{\mathrm{hi}}/(m-1)} - 1 \ge \frac{\tau_{\mathrm{hi}}}{m - 1}`$;
3. For each $`0 \le s \le i-1`$, $`\frac{N_1 - 1 - s}{m - 1} \ge a_+ - \frac{s}{m - 1} \ge \left(a_+ - \frac{s}{m_0 - 1}\right)_+`$.
Multiplying these inequalities yields the stated bound (`H_raw_lower_int`). $`\square`$

**Lemma 4.2 (Uniform Upper Enclosure of $`B_{N, r}(x)`$).**  
*Let $`t_{\mathrm{up}}(N) := \log\frac{N-1}{m-1}`$. For any $`m \ge m_0`$, any $`\tau \ge \tau_{\mathrm{lo}} > 0`$ with $`u_0 := \frac{\tau_{\mathrm{lo}}}{m_0 - 1} \le \frac{1}{2}`$, any $`r \in (0, 3)`$, and any integer $`N \ge m`$ with $`t_{\mathrm{up}}(N) \ge t_{\mathrm{lo}}`$, set $`L_{\mathrm{lo}} := \tau_{\mathrm{lo}} e^{t_{\mathrm{lo}}}`$ and $`B_0 := (m_0 - 1)e^{t_{\mathrm{lo}}} + 1 - r`$. Then substituting $`1 - u = q e^{-v/B}`$ ($`B := N_{\mathrm{lo}} - r`$ at real $`N_{\mathrm{lo}} := 1 + (m-1)e^{t_{\mathrm{lo}}} \le N`$) gives the exact identity*

```math
B_{N, r}(x) \le B_{N_{\mathrm{lo}}, r}(x_{\mathrm{lo}}) = e^{-L_{\mathrm{lo}}} \cdot \frac{\Gamma(B+r)}{\Gamma(B)B^r} \cdot \frac{1}{\Gamma(r)} \int_0^\infty \big(w + B(1 - e^{-v/B})\big)^{r-1} e^{-v}\,dv, \qquad w := B\frac{x_{\mathrm{lo}}}{1 - x_{\mathrm{lo}}},
```

*which is bounded above uniformly in $`m \ge m_0`$ by:*
1. ***For $`0 < r < 1`$:*** *$`\frac{\Gamma(B+r)}{\Gamma(B)B^r} \le 1`$ (Gautschi's inequality), $`B \ge B_0`$, and $`w > L_{\mathrm{lo}}`$. Splitting $`[0, \infty)`$ at $`V_{\mathrm{cut}} := \min(24, \lfloor B_0/2 \rfloor)`$ with $`\alpha := 1 - \frac{V_{\mathrm{cut}}}{2 B_0}`$ (using $`B(1 - e^{-v/B}) \ge \alpha v`$ on $`[0, V_{\mathrm{cut}}]`$ and $`\ge \alpha V_{\mathrm{cut}}`$ on $`[V_{\mathrm{cut}}, \infty)`$) gives*

```math
B_{N, r}(x) \le e^{-L_{\mathrm{lo}}} \left[ \alpha^{r-1} e^{L_{\mathrm{lo}}/\alpha} Q\!\left(r, \frac{L_{\mathrm{lo}}}{\alpha}\right) + \frac{(L_{\mathrm{lo}} + \alpha V_{\mathrm{cut}})^{r-1} e^{-V_{\mathrm{cut}}}}{\Gamma(r)} \right].
```

2. ***For $`r = k + a \in [1, 3)`$ ($`k \in \{1, 2\}, a \in [0, 1)`$):*** *$`B \le (m-1)e^{t_{\mathrm{lo}}}`$, so $`w \le w_{\mathrm{hi}} := L_{\mathrm{lo}}\frac{e^{u_0} - 1}{u_0}`$, $`B(1 - e^{-v/B}) \le v`$, and $`\frac{\Gamma(B+r)}{\Gamma(B)B^r} \le C_r(B_0)`$ where $`C_{1+a}(B_0) = 1 + \frac{a}{B_0}`$ and $`C_{2+a}(B_0) = (1 + \frac{1+a}{B_0})(1 + \frac{a}{B_0})`$, giving*

```math
B_{N, r}(x) \le e^{-L_{\mathrm{lo}}} C_r(B_0) \cdot e^{w_{\mathrm{hi}}} Q(r, w_{\mathrm{hi}}).
```

*Here $`Q(r, v) := \frac{1}{\Gamma(r)}\int_v^\infty t^{r-1}e^{-t}\,dt`$ is evaluated outward-rounded via odd-truncated alternating Taylor series for $`v \le 2`$ and the $`(2K+1)`$-st Stieltjes continued-fraction convergent ($`K = 45`$, guaranteed strict upper bound for $`0 < a < 1`$) for $`v > 2`$, extended to $`r \in [1, 3)`$ by the exact recurrence $`e^v Q(1+a, v) = e^v Q(a, v) + v^a/\Gamma(1+a)`$ (`B_raw_upper_int`).*

### 4.2 Exact Discrete-Time Affine Knot-Crossing Theorem

The key observation that eliminates all time-discretization error for **all $`m \ge m_0 = 200`$** is that the discrete time increments satisfy the universal one-sided inequalities

```math
N\Delta t_{\mathrm{lo}}(N) = N\log\!\left(1 + \frac{1}{N}\right) < 1, \qquad N\Delta t_{\mathrm{up}}(N) = N\log\!\left(1 + \frac{1}{N-1}\right) > 1 \qquad (\forall N \ge 2),
```

while both step sizes satisfy $`\Delta t(N) < \frac{1}{m_0 - 1} = \frac{1}{199} < \frac{1}{128} = h`$ so a single step $`N \to N+1`$ crosses **at most one grid knot $`t_{i+1}`$**.

**Lemma 4.3 (Affine Knot-Crossing Lemma).**  
*Let $`f(t, s)`$ be continuous and piecewise-linear on $`[0, T + h]`$ with nonpositive slope $`-\sigma_i(s) \le 0`$ on each slab $`[t_i, t_{i+1}]`$ (so $`f(t, s) = f_{i+1}(s) + \sigma_i(s)(t_{i+1} - t)`$ for $`t \in [t_i, t_{i+1}]`$), and consider a Markov transition at tree size $`N`$ from $`(t, s)`$ ($`t \in [t_i, t_{i+1})`$) to $`(t', S_{N+1})`$ with $`t' = t + \Delta t`$, $`0 < \Delta t < h`$, and transition probabilities $`p_{s \to s'} = r_{s \to s'}/N`$ ($`s' \neq s`$, $`\sum_{s' \neq s} r_{s \to s'} = r_{\mathrm{tot}} \le N`$).*
1. ***(Lower submartingale, $`N\Delta t < 1`$):*** *Suppose slab $`i`$ satisfies the two endpoint conditions*

```math
\text{(L-R)}\quad \sigma_i(s) \le \sum_{s' \neq s} r_{s \to s'}\big(f_{i+1}(s') - f_{i+1}(s)\big), \qquad \text{(L-L)}\quad \sigma_i(s) \le \frac{\sum_{s' \neq s} r_{s \to s'} f_i(s') - r_{\mathrm{tot}} f_{i+1}(s)}{1 + h r_{\mathrm{tot}}}.
```

*Then $`\mathbb{E}[f(t', S_{N+1})] - f(t, s) \ge 0`$ for EVERY $`t' \in [t_i, t_{i+1}]`$, AND for every knot-crossing step $`t' \in (t_{i+1}, t_{i+2}]`$ provided either (a) $`(i+1, s)`$ is also a continuing state satisfying (L-R) and (L-L) on slab $`i+1`$, or (b) $`(i+1, s)`$ is a stopping state satisfying the explicit knot condition*

```math
\text{(L-KNOT)}\quad -\sigma_{i+1}(s) + \sum_{s' \neq s} r_{s \to s'}\Big[f_{i+1}(s') - f_{i+1}(s) - \big(\sigma_{i+1}(s') - \sigma_{i+1}(s)\big)\delta_2\Big] \ge 0 \qquad (\delta_2 \in \{0, \Delta t_{\max}^{\mathrm{cor}}\}).
```

2. ***(Upper supermartingale, $`N\Delta t > 1`$):*** *Suppose slab $`i`$ satisfies*

```math
\text{(U-R)}\quad \sigma_i(s) \ge \sum_{s' \neq s} r_{s \to s'}\big(f_{i+1}(s') - f_{i+1}(s)\big), \qquad \text{(U-L)}\quad \sigma_i(s) \ge \frac{\sum_{s' \neq s} r_{s \to s'} f_i(s') - r_{\mathrm{tot}} f_{i+1}(s)}{1 + h r_{\mathrm{tot}}}.
```

*Then $`\mathbb{E}[f(t', S_{N+1})] - f(t, s) \le 0`$ for EVERY $`t' \in [t_i, t_{i+1}]`$, AND for every knot-crossing step $`t' \in (t_{i+1}, t_{i+2}]`$ whenever $`(i+1, s)`$ is a continuing state satisfying (U-R) and (U-L) on slab $`i+1`$. At the terminal horizon $`t_M = T`$, lifting $`f_M^-(s) := f_{\mathrm{stop}}(s) + \Delta t_{\max}(M-1)\big(\sum_{s' \neq s} r_{s \to s'}(f_{\mathrm{stop}}(s') - f_{\mathrm{stop}}(s))\big)_+`$ guarantees $`\mathbb{E}[f_{\mathrm{stop}}(S_{N+1})] \le f(t, s)`$ for every step crossing $`t_M`$.*

*Proof.* Write the exact one-step drift $`N\big(\mathbb{E}[f(t', S_{N+1})] - f(t, s)\big) = N\big(f(t', s) - f(t, s)\big) + \sum_{s' \neq s} r_{s \to s'}\big(f(t', s') - f(t', s)\big)`$.
- **Interior ($`t' \in [t_i, t_{i+1}]`$):** $`N(f(t', s) - f(t, s)) = -(N\Delta t)\sigma_i(s)`$ (which is $`\ge -\sigma_i(s)`$ if $`N\Delta t < 1`$, and $`\le -\sigma_i(s)`$ if $`N\Delta t > 1`$), while $`\sum_{s' \neq s} r_{s \to s'}(f(t', s') - f(t', s))`$ is affine in $`t' \in [t_i, t_{i+1}]`$, so its extrema occur at $`t' = t_{i+1}`$ (giving (L-R)/(U-R)) and $`t' = t_i`$ (using $`f_i(s) = f_{i+1}(s) + h\sigma_i(s)`$, giving (L-L)/(U-L)).
- **Knot crossing ($`t' \in (t_{i+1}, t_{i+2}]`$):** Write $`\delta_1 = t_{i+1} - t \ge 0`$, $`\delta_2 = t' - t_{i+1} > 0`$, $`\delta_1 + \delta_2 = \Delta t`$. Then

```math
N\big(\mathbb{E}[f(t', S_{N+1})] - f(t, s)\big) = -N\sigma_i(s)(\Delta t - \delta_2) - N\sigma_{i+1}(s)\delta_2 + \sum_{s' \neq s} r_{s \to s'}\big(f(t_{i+1} + \delta_2, s') - f(t_{i+1} + \delta_2, s)\big),
```

which is **affine in $`\delta_2 \in [0, \Delta t]`$**! At $`\delta_2 = 0`$, it reduces to the $`t' = t_{i+1}`$ condition on slab $`i`$ ((L-R)/(U-R)). At $`\delta_2 = \Delta t`$, it equals the interior drift on slab $`i+1`$ at $`t' = t_{i+1} + \Delta t \in [t_{i+1}, t_{i+2}]`$ (already certified if $`(i+1, s)`$ is continuing, or directly checked by (L-KNOT) if $`(i+1, s)`$ is a lower stopping state).
- **Terminal lift ($`t < t_M \le t'`$ in the upper bound):** For $`t \in [t_{M-1}, t_M)`$, $`f(t, s) \ge f_M^-(s)`$ and $`1/N < 1/(N-1) \le \Delta t_{\max}(M-1)`$, so $`\mathbb{E}[f_{\mathrm{stop}}(S_{N+1})] = f_{\mathrm{stop}}(s) + \frac{1}{N}\sum_{s' \neq s} r_{s \to s'}(f_{\mathrm{stop}}(s') - f_{\mathrm{stop}}(s)) \le f_M^-(s) \le f(t, s)`$. $`\square`$

### 4.3 Tracked Markov States (Exact DAG Ordering)
Both lower and upper state graphs are acyclic within each slab $`i`$, enabling exact computation of $`f_i(s')`$ before $`f_i(s)`$:
1. **Lower submartingale** (at $`\tau_{\mathrm{hi}}`$, returning $`M_{\mathrm{lo}} := \max(V_0(2,1), H_{\mathrm{raw\_lower}}(1, 0))`$, $`H_{\mathrm{host}}(p) := \max(H_0(p), M_{\mathrm{lo}})`$, $`L_{\mathrm{lo}}(k) := R_0(1, k)`$):
   - **Tracked positive stock $`V(d, j)`$** ($`2 \le d \le D_{\max}, 1 \le j < d`$): degree $`d`$, $`j`$ tracked leaves; gain rate $`j \to V(d+1, j+1)`$, loss rate $`j/d \to V(d, j-1)`$ (or $`R(d, 2)`$ if $`j = 1`$); immediate-stop option $`H_{\mathrm{raw\_lower}}(j, t_{i+1})`$ (mandatory at $`d = D_{\max}`$ and whenever (L-KNOT) fails).
   - **Zero-stock return $`R(d, k)`$** ($`1 \le d < D_{\max}, 2 \le k \le R_{\max}`$): degree $`d`$, $`0`$ leaves, one tracked neighbor of degree $`k`$; return rate $`1/k \to V(d+1, 1)`$, neighbor growth rate $`\le k - 1 + 1/d \to R(d, k+1)`$ (zero at $`k = R_{\max}+1`$).
   - **Serial host bonus $`H(p)`$** ($`2 \le p \le H_{\mathrm{cut}}`$, handing off to $`V(2, 1)`$ at $`p = H_{\mathrm{cut}}+1`$): degree-$`2`$ serial middle with $`1`$ leaf and $`1`$ host neighbor of degree $`p`$; gain rate $`1 + 1/p \to V(3, 2)`$, loss rate $`1/2 \to R(2, 2)`$, host growth rate $`r_{\mathrm{host}} \in [0, p - 1/2] \to H(p+1)`$ (both endpoints checked).
   - **Pathwise monotonicity:** At every knot $`i`$, `cert_all_m.cpp` explicitly asserts $`V_i(d+1, j+1) \ge V_i(d, j)`$, $`R_i(d, k) \ge R_i(d, k+1)`$, $`V_i(d+1, 1) \ge R_i(d, 2)`$, and $`H_i(p) \ge H_i(p+1)`$. By piecewise linearity, these hold for all real $`t \in [0, T+h]`$, justifying dropping extra gains and inflating competitor growth rates.
2. **Upper supermartingale** (at $`\tau_{\mathrm{lo}} > 0`$; at $`\tau_{\mathrm{lo}} = 0`$ uses Lemma 3.3 `U_hub`):
   - **Post-promotion 1-race $`P(k)`$** ($`2 \le k \le K_{\max}`$; stopped at $`k = K_{\max}+1`$ using $`t_{\max(0, i-1)}`$ so linear interpolation on $`[t_i, t_{i+1}]`$ stays $`\ge B_{\mathrm{raw\_upper}}(1 + 1/(K_{\max}+1), t_i)`$): degree-$`2`$ promoted vertex with $`1`$ leaf and $`1`$ host of degree $`k`$; gain rate $`1 + 1/k \to`$ stop with $`B_{\mathrm{raw\_upper}}(2 + 1/k, t_i)`$, leaf-loss rate $`1/2 \to`$ stop with $`B_{\mathrm{raw\_upper}}(1/2 + 1/k, t_i)`$, host-growth rate $`g_w \in \{1/2, k - 1/2\} \to P(k+1)`$ (both endpoints checked; monotone lift $`P_i(k) \ge P_i(k+1)`$ enforced).
   - **Pending adverse leaf $`U(k, j)`$** ($`2 \le k \le K_{\max}, 0 \le j < k`$; stopped at $`k = K_{\max}+1`$ with $`\min(B_{\mathrm{raw\_upper}}(1/(K_{\max}+1), t_{\max(0,i-1)}), U_{\mathrm{hub}})`$): promotion rate $`1/k \to P(k)`$, sibling loss rate $`j/k \to U(k, j-1)`$, host gain rate $`g \in \{j+1, (k+j+1)/2\} \to U(k+1, j+1)`$ (both endpoints checked). At $`t = 0`$, capped by $`\min(U_0(k, j), B_{\mathrm{raw\_upper}}(1/k, 0), U_{\mathrm{hub}}(k, j))`$.

---

## 5. Exact Integer-$`N`$ Backward Induction for $`4 \le m \le 199`$ (`cert_small_m.cpp`)

For every finite $`m = b + 1 \in \{4, \dots, 199\}`$, parameterize $`x \in [0, 1)`$ globally by $`y := \log\frac{1}{1 - x} \in [0, \infty)`$. By Section 7 (High-$`\tau`$ Tail), every pair $`(m, y)`$ with $`(m - 1)y \ge 6`$ is already certified, so only $`y \in [0, 2]`$ ($`m_{\max}(y) := \min(199, \lfloor 6/y \rfloor + 1)`$) remains.

On each exact rational box $`[y_{\mathrm{lo}}, y_{\mathrm{hi}}] \subset [0, 2]`$:
- If $`199 y_{\mathrm{hi}} \le 6`$: initialize all state arrays $`(V, R, H, U, P)`$ at $`N_{\mathrm{start}} = 200`$ using Section 4's exact `init_at_200` solver on $`[199 y_{\mathrm{lo}}, 199 y_{\mathrm{hi}}]`$.
- If $`199 y_{\mathrm{lo}} \ge 6`$: initialize directly at $`N_{\mathrm{start}} := \max(m_{\max} + 4, \min(900, \lceil 24/y_{\mathrm{lo}} \rceil + 2))`$ using exact fixed-$`N`$ `H_exact_N_lower` and `B_exact_N_upper`.
- Then perform **one literal exact integer-$`N`$ backward sweep** $`N = N_{\mathrm{start}} - 1, \dots, 4`$:
  - Every 1-step transition probability at tree size $`N`$ is an exact rational with denominator $`N d`$, $`N k d`$, or $`2 N k`$, rounded **down** for $`(V, R, H)`$ (taking $`\max`$ with exact $`H_{N, j}(x_{\mathrm{hi}})`$) and **up** for $`(P, U)`$ (taking $`\min`$ with exact $`B_{N, r}(x_{\mathrm{lo}})`$ and $`U_{\mathrm{hub}}`$).
  - For $`4 \le N < 60`$, `B_exact_N_upper` computes $`\frac{\Gamma(N)}{\Gamma(r)\Gamma(N-r)B^r}`$ in outward-rounded log-space via `iv_gamma_rat` and bounds $`\int_0^\infty (w + B(1-e^{-v/B}))^{r-1} e^{-v}\,dv`$ by exact outward-rounded monotone quadrature (decreasing integrand for $`r < 1`$ with exact analytic antiderivative on $`[0, 1/64]`$; increasing integrand bounded by $`(w+B)^{r-1}`$ for $`r > 1`$).
  - At every size $`N \in \{4, \dots, m_{\max}\}`$, the exact snapshot $`(V, R, H, U, P)`$ at size $`N`$ is checked against the exact finite-tree signed ledger for $`m = N`$ (maximum vertex degree $`d \le m - 3`$).

---

## 6. Dimensionless Normalized All-Tree Signed Ledger

Fix any box and normalize all lower bounds $`(H_{\mathrm{host}}(p), L_{\mathrm{lo}}(k))`$ downward and all upper bounds $`U_{\mathrm{hi}}(k, j)`$ upward by dividing by $`M_{\mathrm{lo}} > 0`$ (in integer scale $`L_{\mathrm{SC}} = 10^9`$), so $`\widetilde{M} = 1`$, $`\widetilde{H}(p) \ge 1`$, $`g_{\mathrm{base}} := 1 - \widetilde{U}(2, 0) > 0`$, and for $`k > K_{\max}+1`$, $`\widetilde{U}(k, j) := \min\big(\frac{1}{(j+1)M_{\mathrm{lo}}}, \widetilde{u}_{\mathrm{tail0}}\big)`$.

### 6.1 Tree-Potential Charging Identity
Let $`T_n`$ be any finite tree of size $`n \ge 2`$ ($`m = n+2`$).
- **Seed ($`n = 2`$, $`T_2`$ = edge):** $`s(1) = s(2) = 1/2`$, no edge debit; Lemma 2.1 reduces to $`\widetilde{H}(2) > \widetilde{U}(2, 0)`$.
- **Stars ($`n = d + 1`$, $`S_{d+1}`$):** center has degree $`d \ge 2`$, $`c = 0`$, $`l = d`$, $`r = 0`$, $`s = \frac{1}{d+1}`$, leaf weight $`t = \frac{1}{2d}`$; Lemma 2.1 reduces to $`C_{\mathrm{lo}}(d, 0, 0) > 0`$ (checked for all $`d \le D_{\mathrm{tail}} = 300`$, and for $`d > 300`$ by $`\frac{14}{10} g_{\mathrm{base}} > \frac{1}{2}\widetilde{U}(302, 300)`$, since $`\frac{d}{d+1} + \frac{1}{2} > 1.49 > 1.4`$ and $`\frac{d-1}{2(d+1)} < \frac{1}{2}`$).
- **Multi-core trees (not an edge or star):** Deleting all leaves yields the core tree $`C(T_n)`$ with $`|V(C(T_n))| \ge 2`$. Every vertex $`u \in C(T_n)`$ has $`d(u) \ge 2`$, core degree $`c(u) \ge 1`$, leaf count $`l(u) = d(u) - c(u)`$, and receive sum $`r(u) := \sum_{v \sim u,\, v \in C(T_n)} 1/d(v) \in (0, c(u)/2]`$. Contracting maximal paths of bare degree-$`2`$ vertices ($`d = c = 2, l = 0`$) into chain-edges yields a site tree $`\mathcal{T}_S`$ on $`|S| \ge 2`$ sites with $`|E(\mathcal{T}_S)| = |S| - 1`$ edges and $`\sum_{u \in S}(c(u) - 2) = -2`$.

Choose nonnegative tariffs $`(a, b_{\mathrm{base}}, \mathrm{fee}_{\mathrm{lo}}, \mathrm{fee}_{\mathrm{hi}}, \mathrm{fee}_t, \mathrm{xfer}, D_{\mathrm{hi}})`$. Define:
1. **Site credit (after paying all incident leaf edges):**

```math
C_{\mathrm{lo}}(d, c, r) := d s(\widetilde{H}(d+1) - \widetilde{U}(2,0)) + \frac{l}{2d}(\widetilde{H}(2) - \widetilde{U}(2,0)) - l(t - s)_+(\widetilde{U}(2,0) - \widetilde{L}(d+1))_+ - l(s - t)_+(\widetilde{U}(d+1, l-1) - \widetilde{L}(2))_+,
```

where $`s = \frac{l+r}{d(d+1)}`$ and $`t = \frac{1}{2d}`$.
2. **Direct site-site edges $`\{u, v\}`$:** At most one of $`s(u) - s(v)`$ and $`s(v) - s(u)`$ is positive, so the edge debit is bounded by $`b_{\mathrm{base}} + (J_{u \to v} - b_{\mathrm{base}})_+ + (J_{v \to u} - b_{\mathrm{base}})_+`$. Each direct edge between a high site ($`d \ge D_{\mathrm{hi}}`$) and a low site ($`e \le 4`$) transfers $`\mathrm{xfer}`$ from the high site to the low site (net sum $`0`$).
3. **Contracted chains:** Each chain port at site $`u`$ pays fee $`F_{\mathrm{ch}}(u)`$ ($`\mathrm{fee}_t`$ if $`d=2`$, $`\mathrm{fee}_{\mathrm{lo}}`$ if $`3 \le d < D_{\mathrm{hi}}`$, $`\mathrm{fee}_{\mathrm{hi}}`$ if $`d \ge D_{\mathrm{hi}}`$) from $`u`$ to the chain (net sum $`0`$). For every bare degree-$`2`$ chain $`v_1 - \dots - v_L`$ ($`L \ge 1`$) between sites $`x, y`$, its net debt (sum of internal edge debits + port debits minus internal serial credits minus $`F_{\mathrm{ch}}(x) + F_{\mathrm{ch}}(y)`$) is verified to be strictly $`< b_{\mathrm{base}}`$ for $`L \in \{1, 2, 3\}`$. For $`L \ge 4`$, every extra interior vertex $`v_i`$ ($`3 \le i \le L-2`$) has $`s(v_i) = 1/6`$ equal to both neighbors (zero extra edge debit) and positive serial credit $`\frac{1}{3}(\widetilde{H}(3) - \widetilde{U}(2, 0)) > 0`$, so $`\mathrm{debt}(L) \le \mathrm{debt}(3) < b_{\mathrm{base}}`$.
4. **Site balance:** Every site $`u \in S`$ satisfies $`R(u) := C_{\mathrm{lo}}(u) + a(c(u) - 2) - b_{\mathrm{base}} - \sum_{\text{chains}} F_{\mathrm{ch}}(u) \pm \mathrm{xfer} \ge \mathrm{exc}(u) := \sum_{\text{direct } v} (J_{u \to v} - b_{\mathrm{base}})_+`$.

Summing over all sites $`u \in S`$ and all edges of $`\mathcal{T}_S`$ (using $`\sum_{u \in S} (a(c(u) - 2) - b_{\mathrm{base}}) = -2a - (|S|-1)b_{\mathrm{base}} - b_{\mathrm{base}}`$) yields:

```math
\frac{n(n+1)}{M_{\mathrm{lo}}}\Big(\mathbb{P}(Z_b > x \mid T_n) - \mathbb{P}(Z_m > x \mid T_n)\Big) \ge \sum_{u \in S}\big(R(u) - \mathrm{exc}(u)\big) + \sum_{\text{chains } P}\big(b_{\mathrm{base}} - \mathrm{debt}(P)\big) + (2a + b_{\mathrm{base}}) \ge 2a + b_{\mathrm{base}} > 0.
```

### 6.2 Completeness of Finite Checks & Convexity/Concavity Reductions
- **Low sites ($`2 \le d < D_{\mathrm{hi}} = 10`$):** For each $`(d, c)`$ and chain count $`k \in \{0, \dots, c\}`$, `R(r) - exc(r)` is concave in $`r`$ (since $`-(t-s)_+`$, $`-(s-t)_+`$, and $`-(J_{u \to v} - b_{\mathrm{base}})_+`$ are all concave in $`s`$), so checking the two extreme values of $`r`$ suffices. If the coarse check fails, exact multisets of direct neighbor degrees in $`\{2, \dots, D_{\mathrm{hi}}-1\} \cup \{\ge D_{\mathrm{hi}}\}`$ are enumerated with exact receive denominator $`\mathrm{lcm}(2, \dots, 10) = 2520`$, using opposite coefficient lower bound $`s(v) \ge \frac{1}{d e(e+1)}`$ and favorable payout $`\widetilde{L}(e+1)`$.
- **High finite sites ($`D_{\mathrm{hi}} \le d \le D_{\mathrm{tail}} = 300`$):** If $`s_{\max}\widetilde{U}(d+1, l) \le b_{\mathrm{base}}`$, then $`\mathrm{exc} \equiv 0`$, so $`R`$ is jointly concave in $`(k, p_{\mathrm{pay}}, r)`$ and the 3 simplex vertices $`\{(0,0), (c,0), (0,c)\} \times \{r_{\min}, r_{\max}\}`$ are exact; otherwise we iterate over all integer $`k \in \{0, \dots, c\}`$, for each of which $`R - \mathrm{exc}`$ is jointly concave in $`(p_{\mathrm{pay}}, r)`$, checking $`\{0, c-k\} \times \{r_{\min}, r_{\max}\}`$.
- **Analytic high-degree site tail ($`d > D_{\mathrm{tail}} = 300`$ in `cert_all_m.cpp`):**
  1. $`s\widetilde{U}(d+1, l) \le \frac{\widetilde{U}(302, 0)}{302} \le b_{\mathrm{base}}`$, so $`\mathrm{exc} \equiv 0`$.
  2. Since $`s - t \le \frac{l-1}{2d(d+1)}`$, core-side leaf debit is at most $`T_{\mathrm{cl}} := \min\big(\frac{1}{2(D_{\mathrm{tail}}+2)M_{\mathrm{lo}}}, \frac{\widetilde{U}(302,0)}{2}\big)`$.
  3. Require $`a \ge \max(\mathrm{fee}_{\mathrm{hi}}, \mathrm{xfer})`$. If $`l \le d/2`$ ($`c \ge D_{\mathrm{tail}}/2`$), serial credit minus leaf-side debit is $`\ge \frac{l}{d} z_{\mathrm{lo}} \ge \min(0, z_{\mathrm{lo}}/2)`$ with $`z_{\mathrm{lo}} = \frac{D_{\mathrm{tail}}}{D_{\mathrm{tail}}+1}g_{\mathrm{base}} + \frac{g_{\mathrm{base}} - \widetilde{U}(2,0)}{2}`$, and net tariff adds $`\ge (a - \max(\mathrm{fee}_{\mathrm{hi}}, \mathrm{xfer}))\lfloor D_{\mathrm{tail}}/2 \rfloor - 2a - b_{\mathrm{base}}`$. If $`l > d/2`$, then $`2l \ge d+1`$ so $`(t-s)_+ = 0`$, while serial credit $`\ge (\frac{D_{\mathrm{tail}}}{2(D_{\mathrm{tail}}+1)} + \frac{1}{4})g_{\mathrm{base}}`$ pays $`T_{\mathrm{cl}} + 2a + b_{\mathrm{base}}`$ (we even subtract an extra conservative $`\frac{\widetilde{U}(2,0)}{2(D_{\mathrm{tail}}+2)}`$).
- **Bare degree-$`2`$ chains ($`L \in \{1, 2, 3\}`$):** Port debt is convex in $`s_x`$ (so extreme receives $`2r \in \{1, c\}`$ suffice) and convex in $`(q_x, q_y)`$ on $`q \in [0, \frac{1}{6(K_{\max}+1)}]`$ for the high tail $`d > K_{\max}`$ (so the two convex corners $`q \in \{0, \frac{1}{6(K_{\max}+1)}\}`$ with exact stock-coupled `own_cap` $`\max\big(\frac{\widetilde{u}_{\mathrm{tail0}}}{2(K_{\max}+2)}, \min(\frac{1}{4(K_{\max}+1)M_{\mathrm{lo}}}, \frac{\widetilde{u}_{\mathrm{tail0}}}{K_{\max}+2})\big)`$ suffice).
- **One-line universal shortcut:** Letting $`A := \sum_{u \in V(T_n)} d(u) s(u) > 0`$, we have $`\sum_{\{u,v\}} |s(u) - s(v)| \le A`$ and, for every non-leaf $`v`$ with $`l \ge 1`$ leaves, $`l(t - s(v))_+ \le \frac{l(d-1)}{2d(d+1)} < \frac{1}{3}\left(\frac{l}{d+1} + \frac{l}{2d}\right) \le \frac{1}{3}(d(v)s(v) + \frac{l}{2d})`$, so leaf-originating parallel weight is at most $`A/3`$. Hence if $`\widetilde{U}(2,0) \ge U_3^{\max} := \max_{k \ge 3, 0 \le j < k} \widetilde{U}(k, j)`$ and $`3 > 4\widetilde{U}(2,0) + 2 U_3^{\max}`$, then $`G(T_n) \ge A M_{\mathrm{lo}}(1 - \frac{4}{3}\widetilde{U}(2,0) - \frac{2}{3}U_3^{\max}) > 0`$ on every finite tree $`T_n`$.

---

## 7. Self-Contained Analytic High-$`\tau`$ Tail ($`\tau \ge 6`$, All $`b \ge 3`$)

**Theorem 7.1.** *For every $`b = n + 1 \ge 3`$ ($`m = b+1 \ge 4`$), every tree $`T_n`$, and every $`x \in (0, 1)`$ with $`\tau := b\log\frac{1}{1-x} \ge 6`$,*

```math
\mathbb{P}(Z_b > x \mid T_n) - \mathbb{P}(Z_{b+1} > x \mid T_n) > \frac{19}{10000}\frac{(1-x)^b}{n(n+1)} > 0.
```

*Proof.* Set $`q := 1 - x`$. Since $`b x/q = b(e^{\tau/b} - 1) \ge \tau \ge 6`$, every serial middle has $`\mathbb{P}(Z_b > x) \ge H_{m, 1}(x) = e^{-1} q^b`$ (Lemma 3.1), while adverse tips on hosts of degree $`\ge 2`$ and $`\ge 3`$ satisfy $`\mathbb{P}(Z_m > x) \le B_{m, 1/2}(x)`$ and $`\le B_{m, 1/3}(x)`$ respectively (Lemma 3.2, monotonic in shape). By the one-line inequality of Section 6.2 with $`A = \sum_u w(u)/(d(u)+1) \ge \frac{1}{n}\sum_u w(u) = 1`$,

```math
n(n+1)\big(\mathbb{P}(Z_b > x \mid T_n) - \mathbb{P}(Z_{b+1} > x \mid T_n)\big) \ge q^b\left(e^{-1} - \frac{4}{3}\frac{B_{m, 1/2}(x)}{q^b} - \frac{2}{3}\frac{B_{m, 1/3}(x)}{q^b}\right).
```

For $`r \in \{1/2, 1/3\}`$, put $`s := 1 - r`$, $`B := m - r = b + s \ge 7/2`$, and $`v := q/(B x) \le q/(b x) \le 1/6`$. Substituting $`u = x + q t`$ in $`B_{m, r}(x)`$ and using Gautschi's bound $`\frac{\Gamma(B+r)}{\Gamma(B+1)} \le B^{-s}`$ gives

```math
\frac{B_{m, r}(x)}{q^b} \le \frac{v^s}{\Gamma(r)} \mathbb{E}\!\left[\left(1 + \frac{q}{x} U\right)^{-s}\right], \qquad U \sim \mathrm{Beta}(1, B).
```

Since $`(1+z)^{-s} \le 1 - s z + \frac{s(s+1)}{2} z^2`$ for $`z \ge 0`$, $`\mathbb{E}[U] = \frac{1}{B+1}`$, $`\mathbb{E}[U^2] = \frac{2}{(B+1)(B+2)}`$, and $`\frac{B}{B+1} \ge \frac{7}{9}`$ (as $`B \ge 7/2`$),

```math
\mathbb{E}\!\left[\left(1 + \frac{q}{x} U\right)^{-s}\right] \le 1 - s v \frac{B}{B+1}\left(1 - (1+s)v\frac{B}{B+2}\right) \le 1 - d_s v, \qquad d_s := \frac{7}{9} s\left(1 - \frac{1+s}{6}\right).
```

On $`v \in (0, 1/6]`$, $`v^s(1 - d_s v)`$ is strictly increasing (its derivative is $`v^{s-1}(s - (s+1)d_s v) > 0`$), so it is bounded above by its value at $`v = 1/6`$. Exact outward-rounded `__int128` verification in `cert_all_m.cpp` (`verify_high_tau_tail_exact()`) proves

```math
e^{-1} - \frac{4}{3}\frac{6^{-1/2}}{\Gamma(1/2)}\left(1 - \frac{d_{1/2}}{6}\right) - \frac{2}{3}\frac{6^{-2/3}}{\Gamma(1/3)}\left(1 - \frac{d_{2/3}}{6}\right) > 0.0050397 > \frac{19}{10000}. \quad \square
```

---

## 8. Self-Contained Founder Seam ($`b = 2`$: $`Z_2 \succ_{\mathrm{st}} Z_3`$)

**Theorem 8.1.** *For every $`x \in [0, 1)`$, $`\mathbb{P}(Z_2 > x) > \mathbb{P}(Z_3 > x)`$.*

*Proof.* At $`N = 3`$, vertex $`3`$ attaches uniformly to $`1`$ or $`2`$. Conditioned on $`3 \sim 1`$ (probability $`1/2`$), $`T_3`$ is the path $`2 - 1 - 3`$, in which $`2`$ and $`3`$ are symmetric leaves, contributing $`0`$. Conditioned on $`3 \sim 2`$ (probability $`1/2`$), $`T_3`$ is the path $`1 - 2 - 3`$. Expanding one more step $`T_3 \to T_4`$:
- $`U_3 = 1`$ or $`U_3 = 3`$ (total probability $`2/3`$): vertex $`4`$ attaches to center $`2`$, giving the star $`S_4`$ centered at $`2`$ with leaves $`\{1, 3, 4\}`$. Center $`2`$ has degree $`3`$ and $`3`$ leaves (tracking any $`2`$ of them gives $`\mathbb{P}(Z_2 > x \mid S_4) \ge V_4(3, 2)`$ and also $`\ge H_{4, 3}(x) > H_{4, 1}(x)`$), while leaf $`3`$ is in state $`(k=3, j=2)`$ (so $`\mathbb{P}(Z_3 > x \mid S_4) \le U_4(3, 2)`$ and $`\le B_{4, 1/3}(x)`$).
- $`U_3 = 2, W_3 = 1`$ (probability $`1/6`$): $`T_4 = (4 - 1 - 2 - 3)`$, where $`2`$ is a serial middle of state $`H_4(2)`$ ($`\ge H_{4, 1}(x)`$) and $`3`$ is a serial tip of state $`U_4(2, 0)`$ ($`\le B_{4, 1/2}(x)`$).
- $`U_3 = 2, W_3 = 3`$ (probability $`1/6`$): $`T_4 = (1 - 2 - 3 - 4)`$, where $`2`$ and $`3`$ are symmetric centers, contributing $`0`$.

Hence for every $`x \in [0, 1)`$,

```math
\mathbb{P}(Z_2 > x) - \mathbb{P}(Z_3 > x) = \frac{1}{3}\Big(\mathbb{P}(Z_{\mathrm{center}(S_4)} > x) - \mathbb{P}(Z_{\mathrm{leaf}(S_4)} > x)\Big) + \frac{1}{12}\Big(\mathbb{P}(Z_{\mathrm{mid}(P_4)} > x) - \mathbb{P}(Z_{\mathrm{leaf}(P_4)} > x)\Big).
```

- For $`y = \log\frac{1}{1-x} \in [0, 2]`$: `cert_small_m.cpp` certifies BOTH $`V_4(3, 2) > U_4(3, 2)`$ (using the rigorous $`2/(j+1)`$ cap at $`j = k - 1`$ from Lemma 3.3) AND $`H_4(2) > U_4(2, 0)`$ across all 84 exact contiguous boxes covering $`[0, 2]`$.
- For $`y \ge 2`$: $`\tau = 3 y \ge 6`$, so Section 7 already proves $`H_{4, 1}(x) = e^{-1} q^3 > \frac{4}{3} B_{4, 1/2}(x) > B_{4, 1/2}(x) > B_{4, 1/3}(x)`$, making both terms strictly positive. $`\square`$

---

## 9. Exact Outward-Rounded `__int128` Arithmetic Guarantees

Every inequality in `cert_all_m.cpp` and `cert_small_m.cpp` is evaluated in exact 128-bit signed integer arithmetic (`__int128`) with explicit directed rounding (`div_fl` / `div_cl`) and **zero floating-point trust**:
- Internal transcendental scale $`S_2 = 10^{15}`$, grid probability scale $`S = 10^{12}`$, normalized ledger scale $`L_{\mathrm{SC}} = 10^9`$, chain-coordinate scale $`X = 6 \times 10^8`$.
- `iv_exp_pos(x)` ($`0 \le x \le 16`$, asserted): halving to $`u \le 1/8`$ + 15-term Taylor series + exact geometric remainder $`2 u^{16}/16!`$.
- `iv_exp_neg(x)` (any $`x \ge 0`$): chunked by $`4`$ via $`e^{-4}`$ (overflow-free for all $`x \ge 0`$), bounded by $`[0, 1/S_2]`$ for $`x \ge 42`$ (or $`45`$).
- `iv_log_ge1_raw(y)` ($`y \ge 1`$): power-of-2 reduction to $`[1, 2)`$ + 20-term arctanh series in $`w = (z-1)/(z+1) \in [0, 1/3]`$ + exact remainder $`\frac{9}{8} w^{41}/41`$.
- `iv_gamma_rat(p, q)` ($`0 < p < q`$): shift to $`x = 60 + p/q`$ using Binet's elementary Bernoulli bounds $`\frac{1}{12x} - \frac{1}{360x^3} < \log\Gamma(x) - (x-\frac{1}{2})\log x + x - \frac{1}{2}\log(2\pi) < \frac{1}{12x}`$ and outward-rounded log-space reduction.
- Contiguous interval coverage is explicitly asserted inside both binaries (`covered_up_to == 6 * Q` for $`\tau \in [0, 6]`$ and `covered_up_to == 2 * YQ` for $`y \in [0, 2]`$).

---

## 10. Replay Commands and Exact Verification Evidence

**Source commit reference:** local research commit `4edb12d84426bb5837ea679475ddecd6f5883989` (pinned public baseline `c250a1f2fdd53ced91b7b485e94ac8a16106b9c7`; `paper/` untouched).

### Exact One-Line Replay Command
```bash
cd round3/q3-birth-order && \
g++ -O3 -std=c++17 -Wall -Wextra cert_all_m.cpp -o ./cert_all_m && \
g++ -O3 -std=c++17 -Wall -Wextra cert_small_m.cpp -o ./cert_small_m && \
./cert_all_m && ./cert_small_m && \
python3 -B audit_primitives.py && \
python3 -B planted_defects.py && \
rm -f ./cert_all_m ./cert_small_m
```

### Real Output (Observed at Commit `4edb12d84426bb5837ea679475ddecd6f5883989`, Zero Compiler Warnings, Exit Code 0)
```text
[PRIMITIVES] Gamma(1/2) in [1.772453850880, 1.772453873162]
[PRIMITIVES] Gamma(1/3) in [2.678938534668, 2.678938568626]
[HIGH-TAU EXACT] PASS tau >= 6 (all b >= 3, all trees): inv_e=0.3678794 - debit_hi=0.3628397 = margin +0.0050397 > 0.0019
  PASS [0.00000, 0.00025] | M=0.77308 U20/M=0.7900 | EXACT PASS (tariff#0)
  PASS [0.00025, 0.00100] | M=0.76635 U20/M=0.7953 | EXACT PASS (tariff#0)
  ... (52 contiguous boxes covering [0, 6] exactly) ...
  PASS [5.75000, 6.00000] | M=0.00089 U20/M=0.7462 | EXACT PASS (tariff#0)

=== CERTIFICATE SUMMARY (ALL m >= 200, ALL x in [0, 1)) ===
Boxes PASSED: 52, FAILED: 0, exact contiguous tau-coverage: [0, 600000/100000]
  PASS y=[0.000000, 0.000001] (m=4..199, N_start=200)
  ... (84 contiguous boxes covering y = log(1/(1-x)) in [0, 2] exactly) ...
  PASS y=[1.753769, 2.000000] (m=4..4, N_start=16)

=== SMALL-m CERTIFICATE SUMMARY (ALL m in {4..199} + Founder Seam b=2, ALL x in [0, 1)) ===
Boxes PASSED: 84, FAILED: 0, total (box, m) checks: 12105, exact contiguous y-coverage: [0, 39800000/19900000]
HAS_SCIPY=False
Audit [all] PASSED: {'GAMMA': 9, 'HRAW': 150, 'BRAW': 180, 'HEX': 0, 'BEX': 0}, max_rel_excess_BEX=0.00000
Audit [small] PASSED: {'GAMMA': 0, 'HRAW': 0, 'BRAW': 0, 'HEX': 320, 'BEX': 567}, max_rel_excess_BEX=0.23413
DEFECT 1 (poisoned U20 +10%) CAUGHT: FAIL:lo_site_d2_c1
DEFECT 2 (tau >= 5 cutoff) CAUGHT: margin = -0.029445 < 0
DEFECT 3 (poisoned M_lo -15% at tau in [1,1.2]) CAUGHT: FAIL:hi_tail_AB
ALL 3 PLANTED DEFECTS CAUGHT LOUDLY!
```

### Non-Proof Exploratory Diagnostics (Phase 1 Falsification Probes)
Before constructing the proof above, we also ran falsification probes (strictly diagnostic, not used in the proof):
- `research_order/probe_preplacement_finite.py`: exact rational finite-size degree-tail check $`\mathbb{P}(D_{n+1}(N) \ge k \mid T_n) - \mathbb{P}(D_{n+2}(N) \ge k \mid T_n)`$ over all unlabelled trees of sizes $`n = 3, 4, 5, 6`$ up to $`N = 12`$ (`0` negative cases).
- `research_order/fable_hall_asymptotic.py`: exact rational linear system for $`\lim_{n\to\infty} \mathbb{E}[\mathrm{dem}_{n-c}(T_n)] / \mathbb{E}[\mathrm{sup}_{n-c}(T_n)]`$ in Fable's local Hall criterion (`c = 2..9`, ratios $`2/9 = 0.2222, 257/540 = 0.4759, \dots, 0.5609`$, all $`< 1`$).
- `research_order/probe_adv.cpp`, `probe_landscape.cpp`, `probe_falsify_sweep.cpp`: Monte Carlo sweeps over extremal trees across $`y = (n+2)x \in [0, 5]`$ (all positive).

---

> **End of Part 2 (Sections 4–10).** **Back to Part 1 (Sections 1–3):** [`PROOF.md`](PROOF.md) ([§1 Main Theorems](PROOF.md#1-exact-model-and-main-theorems) · [§2 Signed Identity](PROOF.md#2-exact-pre-placement-signed-identity) · [§3 Endpoint Laws](PROOF.md#3-exact-endpoint-comparison-laws)) | [`REVIEW-Q3.md`](REVIEW-Q3.md) · [`REVIEW-IMPLEMENTATION.md`](REVIEW-IMPLEMENTATION.md) · [`round3/README.md`](../README.md) · [`round3/PROVENANCE.md`](../PROVENANCE.md)
