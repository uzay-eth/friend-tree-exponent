# Certified Logarithmic Height and Diameter Bounds for the Random Friend Tree

**Pinned claim commit:** `93484ba92061166a5d4c9c8f8ccdf730654c2385`.  
**Review and dependencies:** GPT-family exact-table and proof audits; Claude independently checked consequences of the established rates. Uses Biggins’ speed theorem and standard Yule convergence, with no Astra inputs.

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-11.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-11.md).

**Author:** Neon-2 lead (Gemini 4 Neon).  
**Dependencies:** **ZERO Astra dependencies** — self-contained from the definition of the friend-tree process (upper bound uses only elementary discrete-time martingales and Borel–Cantelli; lower bound uses our explicit pathwise-coupled finite-offspring CMJ branching process plus Biggins' 1977 branching-random-walk speed theorem and standard rate-1 Yule convergence). Zero floating-point arithmetic.

---

## Headline Theorem

Let $`(T_n)_{n \ge 2}`$ be the edge-seeded random friend tree started from $`T_2 = \{1, 2\}`$, let $`H(T_n) := \max_{v \in T_n} \mathrm{dist}_{T_n}(1, v)`$ be its height rooted at vertex $`1`$, and let $`\mathrm{diam}(T_n) := \max_{u, v \in T_n} \mathrm{dist}_{T_n}(u, v)`$ be its diameter. Then almost surely:
```math
\frac{5000}{5067} \;\approx\; 0.986777
\;\le\; \liminf_{n \to \infty} \frac{H(T_n)}{\log n}
\;\le\; \limsup_{n \to \infty} \frac{H(T_n)}{\log n}
\;\le\; 1.462248,
\qquad\text{(1)}
```
```math
\frac{10000}{5067} \;\approx\; 1.973554
\;\le\; \liminf_{n \to \infty} \frac{\mathrm{diam}(T_n)}{\log n}
\;\le\; \limsup_{n \to \infty} \frac{\mathrm{diam}(T_n)}{\log n}
\;\le\; 2.924496.
\qquad\text{(2)}
```
*(Note: for comparison, the naive rate-1 Yule domination bound on core promotions gives $`\limsup H(T_n)/\log n \le e \approx 2.718282`$; our depth-tilted core supermartingale cuts this nearly in half to $`1.462248`$, and pairs it with a matching lower bound $`0.986777`$, pinning the logarithmic height constant within a factor of $`1.482`$.)*

Moreover, for **every** finite size $`n \ge 3`$ and every $`c > 1.462248`$, letting $`H_{\mathrm{core}}(T_n) := \max_{u \in K(T_n)} \mathrm{dist}_{T_n}(o, u)`$ (where $`o \in \{1, 2\}`$ is the unique core of $`T_3`$, so $`\mathrm{dist}_{T_n}(1, o) \in \{0, 1\}`$, every leaf lies at distance $`1`$ from some core, and by the triangle inequality $`H_{\mathrm{core}}(T_n) - 1 \le H(T_n) \le H_{\mathrm{core}}(T_n) + 2`$ and $`H_{\mathrm{core}}(T_n) \le \mathrm{diam}(T_n) \le 2 H_{\mathrm{core}}(T_n) + 2`$), we have the explicit non-asymptotic polynomial tail bound (`cert_height.py`):
```math
\mathbb P\big(H_{\mathrm{core}}(T_n) \ge c \log n\big)
\;\le\; \frac{\Phi_z(T_3)}{f_{\min}}\, n^{-(c \log z - \lambda)},
\qquad z = \frac{31}{10},\quad \lambda = \frac{165439}{100000},\quad
\frac{\Phi_z(T_3)}{f_{\min}} = \frac{861558416425}{171115606632} \le 5.034950.
\qquad\text{(3)}
```

---

## Complete Proof

### Part 1: Certified Upper Bound (Pure Discrete-Time Depth-Tilted Core Supermartingale)
Root every tree $`T`$ with $`C(T) \ge 1`$ at the unique initial core $`o`$ of $`T_3`$, and write $`\mathrm{depth}(u) := \mathrm{dist}_T(o, u)`$ for each core $`u \in K(T)`$. Fix $`z = 31/10`$, $`\lambda = 165439/100000`$, $`D_M = 100`$, and $`S = 10^{12}`$. Given an integer-tick table $`F[D][L] \in \mathbb Z_{>0}`$ ($`2 \le D \le D_M`$, $`0 \le L \le D`$), define $`f(D, L) := F[D][L]/S`$ for $`D \le D_M`$ and $`f(D, L) := 1`$ for all $`D \ge D_M + 1`$, and set
```math
\Phi_z(T) \;:=\; \sum_{u \in K(T)} z^{\mathrm{depth}(u)}\, f(D_u, L_u).
```
Under the continuous-time friend-tree generator $`\mathcal L`$ (equivalently $`n`$ times the discrete one-step drift at size $`n`$), every birth affects **only** the attachment vertex $`x`$ (and, if $`x`$ was a leaf attached to a core $`u`$, the pair $`(u, x)`$):
1. **Leaf gain at core $`u \in K(T)`$:** $`(D_u, L_u) \to (D_u + 1, L_u + 1)`$ with $`\mathrm{depth}(u)`$ unchanged, at exact rate $`r_u = L_u + \sum_{w \in K(T), w \sim u} 1/D_w \in [L_u, L_u + (D_u - L_u)/2]`$.
2. **Leaf promotion at core $`u \in K(T)`$:** $`(D_u, L_u) \to (D_u, L_u - 1)`$ at $`\mathrm{depth}(u)`$, creating one newborn core $`(2, 1)`$ at $`\mathrm{depth}(u) + 1`$, at exact rate $`b_u = L_u / D_u \in [0, 1]`$ (including $`C(T) = 1`$, $`T_3`$, where $`(2, 2) \to (2, 1)`$ plus a newborn $`(2, 1)`$ at depth $`1`$ at rate $`2/2 = 1`$).

Summing over $`u \in K(T)`$ decouples completely by vertex:
- For any core $`u`$ with $`D_u \ge D_M + 1`$: $`f(D_u, L_u) = f(D_u + 1, L_u + 1) = f(D_u, L_u - 1) = 1`$, so the per-vertex drift divided by $`z^{\mathrm{depth}(u)}`$ is $`\frac{L_u}{D_u} z\, f(2, 1) - \lambda \le z\, f(2, 1) - \lambda \le 0`$, certified by the single exact rational inequality $`z\, F[2][1] \le \lambda S`$.
- For every finite core state $`2 \le D \le D_M`$, $`0 \le L \le D`$: `cert_height.py` verifies $`F[D][L] > 0`$ and the exact rational drift inequality
  $`\displaystyle r\big(f(D + 1, L + 1) - f(D, L)\big) + \frac{L}{D}\big(f(D, L - 1) - f(D, L) + z\, f(2, 1)\big) \;\le\; \lambda\, f(D, L)`$
  (with the promotion term omitted when $`L = 0`$, where $`L/D = 0`$) at both rate extremes $`r \in \{L, L + (D - L)/2\}`$ across all `5,148` interior state rows ($`\sum_{D=2}^{100}(D+1) = 5{,}148`$; by affinity in $`r`$, this covers every tree $`T`$).

Therefore $`(\mathcal L \Phi_z)(T) \le \lambda \Phi_z(T)`$ **pointwise on every finite tree $`T`$ with $`C(T) \ge 1`$**. In discrete time ($`n \ge 3`$), $`\mathbb E[\Phi_z(T_{n+1}) \mid T_n] = \Phi_z(T_n) + \frac{1}{n}(\mathcal L \Phi_z)(T_n) \le (1 + \lambda/n)\Phi_z(T_n)`$, so induction from deterministic $`T_3`$ yields
```math
\mathbb E\big[\Phi_z(T_n)\big]
\;\le\; \Phi_z(T_3) \prod_{j=3}^{n-1}\!\left(1 + \frac{\lambda}{j}\right)
\;\le\; \Phi_z(T_3)\exp\!\left(\lambda\sum_{j=3}^{n-1}\frac{1}{j}\right)
\;\le\; \Phi_z(T_3)\, n^\lambda \qquad (\forall\, n \ge 3),
```
using $`1 + x \le e^x`$ and $`\sum_{j=3}^{n-1} 1/j \le \int_2^{n-1} dx/x = \log((n-1)/2) < \log n`$. Since $`f_{\min} := \min_{D, L} f(D, L) > 0`$ and $`\Phi_z(T_n) \ge f_{\min} z^{H_{\mathrm{core}}(T_n)}`$, Markov's inequality proves (3) for every $`n \ge 3`$.

For the almost-sure $`\limsup`$: let $`L_z := 2\sum_{k=0}^{25}\frac{1}{2k+1}\left(\frac{z-1}{z+1}\right)^{2k+1} < \log z`$ (all omitted terms are strictly positive) and $`c_{\mathrm{hi}} := \lambda / L_z > \lambda / \log z`$. For any $`\varepsilon > 0`$, setting $`c := c_{\mathrm{hi}} + \varepsilon`$ gives $`c \log z - \lambda > \varepsilon \log z > 0`$, so (3) summed along the dyadic subsequence $`n_k = 2^k`$ ($`k \ge 2`$) is summable: $`\sum_{k=2}^\infty \mathbb P(H_{\mathrm{core}}(T_{2^k}) \ge c\log(2^k)) < \infty`$. By Borel–Cantelli, eventually $`H_{\mathrm{core}}(T_{2^k}) < c\log(2^k)`$ for all $`k \ge k_0`$. Because $`H_{\mathrm{core}}(T_n)`$ is nondecreasing in $`n`$, every $`n \in [2^{k-1}, 2^k]`$ obeys
```math
\frac{H_{\mathrm{core}}(T_n)}{\log n}
\;\le\; \frac{H_{\mathrm{core}}(T_{2^k})}{\log(2^{k-1})}
\;<\; c\,\frac{k}{k - 1} \;\longrightarrow\; c.
```
Taking a countable sequence $`\varepsilon_m \downarrow 0`$ and using $`H(T_n) \le H_{\mathrm{core}}(T_n) + 2`$ and $`\mathrm{diam}(T_n) \le 2 H_{\mathrm{core}}(T_n) + 2`$ proves the upper bounds in (1) and (2). $`\square`$

---

### Part 2: Certified Lower Bound (Pathwise-Coupled Reflected CMJ Tree + Convex Global Certificate)

#### Step 2.1 (Pathwise Coupling to an Independent Finite-Offspring CMJ Branching Process)
Fix integers $`L_{\max} = 50`$ and $`J_{\max} = 50`$. In the continuous-time friend tree $`(T_t)_{t \ge 0}`$ (each vertex carrying an independent rate-1 Poisson clock), consider any newborn core vertex $`v`$ born at time $`s_v`$ in state $`(D_v, L_v) = (2, 1)`$ (one core parent $`c_v(0) = 1`$, one real leaf $`L_v(0) = 1`$). Define a local Markov chain $`(\tilde L_v(u), j_v(u))_{u \ge 0}`$ on $`\{0, 1, \dots, L_{\max}\} \times \{1, \dots, J_{\max} + 1\}`$, started at $`(1, 1)`$ at age $`u = 0`$ (clock time $`t = s_v + u`$), coupled pathwise to $`v`$ as follows:
- Track $`v`$'s own rate-1 Poisson clock ring times (in absolute global time) $`s_v < \sigma_{v, 1} < \sigma_{v, 2} < \dots`$, with independent uniform neighbor-selection marks $`U_{v, 1}, U_{v, 2}, \dots \sim \mathrm{Unif}(0, 1)`$ revealed only at the corresponding ring times. Set $`j_v(u) := 1 + N_v^{\mathrm{self}}(u)`$ (capped at $`J_{\max} + 1`$), where $`N_v^{\mathrm{self}}(u)`$ counts rings of $`v`$'s own clock in $`(s_v, s_v + u]`$.
- While $`\tilde L_v \ge 1`$ and $`j_v \le J_{\max}`$, maintain a designated subset $`\mathcal G_v(u)`$ of size $`\tilde L_v(u)`$ of the real leaves currently adjacent to $`v`$ (initially the single newborn leaf $`\ell_0`$).
  - Whenever any leaf in $`\mathcal G_v(u-)`$ rings, it attaches a new leaf $`w`$ to $`v`$: if $`\tilde L_v(u-) < L_{\max}`$, add $`w`$ to $`\mathcal G_v`$ (so $`\tilde L_v \to \tilde L_v + 1`$ at rate $`\tilde L_v`$); if $`\tilde L_v(u-) = L_{\max}`$, do not add $`w`$ to $`\mathcal G_v`$ ($`\tilde L_v`$ stays at $`L_{\max}`$). Any leaf attached to $`v`$ by a leaf outside $`\mathcal G_v`$ or by a core neighbor of $`v`$ is also ignored by $`\mathcal G_v`$.
  - At the $`j`$-th ring $`\sigma_{v, j}`$ of $`v`$'s own clock ($`1 \le j \le J_{\max}`$): in the true friend tree, $`v`$ promotes a uniform leaf with probability $`\frac{L_v}{L_v + c_v}`$ (say on the event $`U_{v, j} < \frac{L_v}{L_v + c_v}`$). Because $`L_v \ge |\mathcal G_v| = \tilde L_v`$ and $`c_v \le 1 + (j - 1) = j`$ (since only rings of $`v`$'s own clock can promote a leaf adjacent to $`v`$), we have
    $`\displaystyle \frac{L_v}{L_v + c_v} \;\ge\; \frac{\tilde L_v}{\tilde L_v + j}.`$
    Whenever $`U_{v, j} < \frac{\tilde L_v}{\tilde L_v + j}`$, the true process **also** promotes a leaf adjacent to $`v`$ into a genuine newborn $`(2, 1)`$ core child $`w`$ of $`v`$ at depth $`\mathrm{depth}(v) + 1`$ and birth age $`\tau = \sigma_{v, j} - s_v`$; we record $`w`$ as a coupled child of $`v`$ (and $`w`$ starts its own coupled process using only $`w`$'s own clock and $`w`$'s own newly generated leaves). Regardless of the outcome of $`U_{v, j}`$, we remove one leaf from $`\mathcal G_v`$ (removing the promoted leaf if it belonged to $`\mathcal G_v`$, or an arbitrary leaf of $`\mathcal G_v`$ otherwise), so $`\tilde L_v \to (\tilde L_v - 1)_+`$ and $`j_v \to j + 1`$.
- Once $`\tilde L_v = 0`$ or $`j_v = J_{\max} + 1`$, the coupled individual $`v`$ stops emitting children.

By induction over jumps, the pathwise invariants $`L_v(s_v + u) \ge \tilde L_v(u)`$ and $`c_v(s_v + u) \le j_v(u)`$ hold for all $`u \ge 0`$ until $`j_v`$ reaches $`J_{\max} + 1`$. Crucially, across all simultaneously active coupled core individuals $`v`$ at any time $`t`$, their driving vertex sets $`\{v\} \cup \mathcal G_v(t)`$ are **pairwise disjoint** (a leaf has only one core neighbor, and when $`v`$ promotes a leaf $`w`$ at $`\sigma_{v,j}`$, $`w`$ is removed from $`\mathcal G_v`$ if it belonged to $`\mathcal G_v`$ — otherwise an adapted leaf of $`\mathcal G_v`$ is removed — while $`w`$'s own active set $`\{w\} \cup \mathcal G_w`$ consists solely of $`w`$ and leaves newly attached to $`w`$ after $`\sigma_{v,j}`$). In the full natural filtration (with marks $`U_{v,j}`$ revealed only at ring times), at every jump time the future clocks on these pairwise disjoint vertex sets are fresh independent rate-1 Poisson processes: thus the projected coupled forest has history-independent, additive jump rates across active individuals — gain rate $`\tilde L_v \mathbf 1\{\tilde L_v < L_{\max}\}`$ and self-ring rate $`1`$ per active individual $`v`$, using a fresh $`\mathrm{Unif}(0,1)`$ mark to emit a newborn child in state $`(1, 1)`$ with probability $`\tilde L_v/(\tilde L_v + j_v)`$ while updating $`(\tilde L_v, j_v) \to ((\tilde L_v - 1)_+, j_v + 1)`$ — and is nonexplosive (embedded in the rate-1 Yule tree). Because these jump rates coincide identically with running independent copies of the local chain and starting an independent copy at each birth, the coupled descendants of $`v`$ have the **exact law of an i.i.d. Crump–Mode–Jagers (CMJ) branching random walk** embedded as a subtree of the true core tree $`K(T_t)`$ (preserving exact birth times and core depths!), with deterministic offspring cap $`0 \le N \le J_{\max}`$.

#### Step 2.2 (Exact Acyclic Laplace Transform $`\phi(s)`$ and Convex Global Certificate for $`\beta_* \le \beta_{\mathrm{hi}}`$)
For any $`s \ge 0`$, let $`\phi(s; L, j) := \mathbb E_{(L, j)}[\sum_{i=1}^N e^{-s \tau_i}]`$ be the expected discounted offspring count started from state $`(L, j)`$, with $`\phi(s) := \phi(s; 1, 1)`$, and let $`\psi(s; L, j) := -\frac{d}{ds}\phi(s; L, j) = \mathbb E_{(L, j)}[\sum_{i=1}^N \tau_i e^{-s \tau_i}] \ge 0`$. Because every offspring satisfies $`\tau_i \le \sigma_{J_{\max}} \sim \mathrm{Gamma}(J_{\max}, 1)`$ and $`N \le J_{\max}`$, $`\phi(s)`$ is smooth and finite on $`(-1, \infty)`$ and satisfies the exact finite-state first-step linear system (with boundary $`\phi(s; 0, j) = \phi(s; L, J_{\max} + 1) = 0`$):
```math
(s + L + 1)\,\phi(s; L, j) \;=\; L\,\phi(s; L + 1, j) + \frac{L}{L + j} + \phi(s; L - 1, j + 1) \qquad (1 \le L < L_{\max}),
```
```math
(s + 1)\,\phi(s; L_{\max}, j) \;=\; \frac{L_{\max}}{L_{\max} + j} + \phi(s; L_{\max} - 1, j + 1),
```
and differentiating in $`s`$ gives the exact system for $`\psi = -\phi'`$ (with boundary $`\psi(s; 0, j) = \psi(s; L, J_{\max} + 1) = 0`$):
```math
(s + L + 1)\,\psi(s; L, j) \;=\; \phi(s; L, j) + L\,\psi(s; L + 1, j) + \psi(s; L - 1, j + 1) \qquad (1 \le L < L_{\max}),
```
```math
(s + 1)\,\psi(s; L_{\max}, j) \;=\; \phi(s; L_{\max}, j) + \psi(s; L_{\max} - 1, j + 1).
```
Because all coefficients are nonnegative and the dependency graph is strictly acyclic ($`j`$ descending $`J_{\max} \to 1`$, then $`L`$ descending $`L_{\max} \to 1`$), any integer-tick tables $`u_{\mathrm{lo}}, u_{\mathrm{hi}}, v_{\mathrm{lo}}, v_{\mathrm{hi}}`$ (divided by $`S = 10^{12}`$) satisfying the $`\le / \ge`$ inequalities at every state rigorously bracket $`u_{\mathrm{lo}}/S \le \phi(s) \le u_{\mathrm{hi}}/S`$ and $`v_{\mathrm{lo}}/S \le -\phi'(s) \le v_{\mathrm{hi}}/S`$ by immediate backward induction.

Set $`\theta_1 = 865/1000`$, $`\theta_2 = 867/1000`$, and $`\beta_{\mathrm{hi}} = 10134/10000`$, and define
```math
h(\theta) \;:=\; \log \phi(\theta) + \beta_{\mathrm{hi}}\,\theta \qquad (\theta \ge 0).
```
By Hölder's inequality on $`\phi(\theta) = \mathbb E[\sum_{i=1}^N e^{-\theta \tau_i}]`$, $`\phi`$ is log-convex, so $`h`$ is **convex** on $`[0, \infty)`$ with derivative $`h'(\theta) = \beta_{\mathrm{hi}} - \psi(\theta)/\phi(\theta)`$. `cert_height.py` verifies in exact `Fraction` arithmetic:
1. **Supercriticality:** $`\mathbb E[N] = \phi(0) \ge 12.856 > 1`$ (and clearly $`\mathbb P(N \ge 2) > 0`$).
2. **Derivative bracket:** $`v_{\mathrm{lo}}(\theta_1) \ge \beta_{\mathrm{hi}}\,u_{\mathrm{hi}}(\theta_1)`$ (hence $`h'(\theta_1) \le 0`$) and $`v_{\mathrm{hi}}(\theta_2) \le \beta_{\mathrm{hi}}\,u_{\mathrm{lo}}(\theta_2)`$ (hence $`h'(\theta_2) \ge 0`$).
3. **Global positivity $`h(\theta) > 0`$ on $`\mathbb R`$** (using Taylor's theorem with Lagrange remainder $`e^{-x} = \sum_{k=0}^{32} \frac{(-x)^k}{k!} - \frac{e^{-\xi} x^{33}}{33!} < \sum_{k=0}^{32} \frac{(-x)^k}{k!}`$ for every $`x > 0`$):
   - For $`\theta \ge \theta_2`$: since $`h'(\theta_2) \ge 0`$ and $`h`$ is convex, $`h(\theta) \ge h(\theta_2) > 0`$ (verified via $`u_{\mathrm{lo}}(\theta_2)/S > \sum_{k=0}^{32} \frac{(-\beta_{\mathrm{hi}}\theta_2)^k}{k!} > e^{-\beta_{\mathrm{hi}}\theta_2}`$).
   - For $`\theta \in [0, \theta_1]`$: since $`h'(\theta_1) \le 0`$ and $`h`$ is convex, $`h(\theta) \ge h(\theta_1)`$. Moreover, for all $`\theta \in (-1, 0]`$, convexity gives $`h'(0) \le h'(\theta_1) \le 0`$, so $`h(\theta) \ge h(0) = \log\phi(0) > 0`$; for $`\theta \le -1`$, the positive-probability path of $`49`$ leaf gains before the first self-ring reaches $`L = L_{\max} = 50`$, where the next self-ring is an independent $`\mathrm{Exp}(1)`$ emitting with probability $`50/51`$, so $`\phi(\theta) = +\infty`$ and $`e^{\beta_{\mathrm{hi}}\theta}\phi(\theta) = +\infty > 1`$.
   - For $`\theta \in [\theta_1, \theta_2]`$: by convexity, the tangent line at $`\theta_1`$ lies below $`h`$, giving $`h(\theta) \ge h(\theta_1) + h'(\theta_1)(\theta_2 - \theta_1) \ge h(\theta_1) - \mathrm{slope}_{\mathrm{hi}}(\theta_2 - \theta_1)`$ where $`\mathrm{slope}_{\mathrm{hi}} := v_{\mathrm{hi}}(\theta_1)/u_{\mathrm{lo}}(\theta_1) - \beta_{\mathrm{hi}} \ge -h'(\theta_1) \ge 0`$. `cert_height.py` verifies $`u_{\mathrm{lo}}(\theta_1)/S > \sum_{k=0}^{32} \frac{(-x_{1,\mathrm{eff}})^k}{k!} > e^{-x_{1,\mathrm{eff}}}`$ with $`x_{1,\mathrm{eff}} = \beta_{\mathrm{hi}}\theta_1 - \mathrm{slope}_{\mathrm{hi}}(\theta_2 - \theta_1) > 0`$, proving $`h(\theta) > 0`$ simultaneously on $`[0, \theta_1]`$ and $`[\theta_1, \theta_2]`$.

Hence $`\inf_{\theta \in \mathbb R} e^{\beta_{\mathrm{hi}}\theta}\phi(\theta) = \inf_{\theta > 0} e^{\beta_{\mathrm{hi}}\theta}\phi(\theta) > 1`$. Because every child is born at or after $`\sigma_{v, 1} - s_v \sim \mathrm{Exp}(1) > 0`$ a.s., dominated convergence gives $`\phi(\theta) \to 0`$ as $`\theta \to +\infty`$, so $`\frac{-\log\phi(\theta)}{\theta} > 0`$ for large $`\theta`$ (hence $`\beta_* > 0`$). Conversely, on the event that $`v`$'s first self-ring precedes $`\ell_0`$'s first ring (probability $`1/2`$, with conditional ring age $`\mathrm{Exp}(2)`$) and $`U_{v, 1} < 1/2`$ (independent probability $`1/2`$), the first emission contributes $`\frac{1}{2} \cdot \frac{1}{2} \cdot \mathbb E[e^{-\theta\,\mathrm{Exp}(2)}] = \frac{1}{4}\cdot\frac{2}{\theta + 2} = \frac{1}{2(\theta + 2)}`$, giving $`\phi(\theta) \ge \frac{1}{2(\theta + 2)}`$ for all $`\theta > 0`$. Thus $`\frac{-\log\phi(\theta)}{\theta} \to 0`$ as $`\theta \to +\infty`$ while $`\frac{-\log\phi(\theta)}{\theta} \to -\infty`$ as $`\theta \downarrow 0`$ (since $`\phi(0) > 1`$), so the continuous function $`\theta \mapsto \frac{-\log\phi(\theta)}{\theta}`$ attains its positive maximum at a finite $`\theta_* \in (0, \infty)`$, where $`h(\theta_*) > 0`$ proves the strict inequality
```math
0 \;<\; \beta_* \;:=\; \sup_{\theta > 0} \frac{-\log \phi(\theta)}{\theta} \;<\; \beta_{\mathrm{hi}} \;=\; \frac{10134}{10000}.
\qquad\text{(4)}
```

#### Step 2.3 (Biggins Speed, Two-Branch Diameter, Independent Restarts, and Discrete-Time Conversion)
Let $`M_g := \min_{|u| = g} S(u)`$ be the earliest birth time in generation $`g`$ of the coupled supercritical CMJ branching random walk (with $`0 \le N \le J_{\max}`$, $`\mathbb E[N] \ge 12.856042 > 1`$, $`\phi(\theta) < \infty`$ on $`(-1, \infty)`$, and $`\mathbb P(N = 0) \ge \mathbb P(\sigma_{v,1} < \tau_{\ell_0},\, U_{v,1} \ge 1/2) = 1/4`$ so $`0 < \rho < 1`$). By Biggins' branching-random-walk speed theorem (Biggins 1977, *Chernoff's theorem in the branching random walk*, J. Appl. Probab. 14, 630–636, Theorem 2; stated directly in on-survival speed form in Biggins 2010, *Branching out*, Theorem 3.1; applied to displacements $`X_i = -\tau_i \le 0`$ so $`\max_{|u|=g} X(u) = -M_g`$ and $`\gamma = \sup\{a : \inf_{\theta \in \mathbb R} e^{-\theta a} \phi(\theta) \ge 1\} = -\beta_* > -\beta_{\mathrm{hi}}`$), on the survival event $`\mathcal S_{\mathrm{surv}}`$ (of probability $`\rho > 0`$),
```math
\lim_{g \to \infty} \frac{M_g}{g} \;=\; \beta_* \;<\; \beta_{\mathrm{hi}} \qquad \text{almost surely on } \mathcal S_{\mathrm{surv}}.
```
Consequently, on $`\mathcal S_{\mathrm{surv}}`$ (started at any finite birth time $`s < \infty`$ from a newborn $`(2, 1)`$ core $`v`$), taking a countable sequence $`\varepsilon_m \downarrow 0`$, for every large $`t`$ generation $`g = \lfloor (t - s)/(\beta_* + \varepsilon_m) \rfloor`$ has a vertex born by time $`t`$ at depth $`\ge g`$, giving $`\liminf_{t \to \infty} H_{\mathrm{core}}(T_t)/t \ge 1/\beta_* > 1/\beta_{\mathrm{hi}}`$. Moreover, let $`Y_g`$ count the individuals in generation $`g`$ whose own independent CMJ subtrees survive forever; before revealing reproduction from generation $`g`$ onward, $`Y_g \mid Z_g \sim \mathrm{Bin}(Z_g, \rho)`$. The event that the tree survives with at most one infinite line of descent implies $`Y_g = 1`$ for every $`g \ge 0`$, so its probability is bounded for every $`g`$ by
```math
\mathbb P(Y_g = 1) \;=\; \mathbb E\big[\mathbf 1_{\{Z_g \ge 1\}}\, Z_g \rho (1 - \rho)^{Z_g - 1}\big] \;\longrightarrow\; 0 \qquad (g \to \infty)
```
by dominated convergence, because the integrand is bounded by $`1`$ and converges to $`0`$ almost surely both on extinction (where eventually $`Z_g = 0`$) and on $`\mathcal S_{\mathrm{surv}}`$ (where $`Z_g \to \infty`$ a.s.). Hence a.s. on $`\mathcal S_{\mathrm{surv}}`$ there exists a generation $`g_1`$ containing two distinct vertices $`w_1 \ne w_2`$ whose independent CMJ subtrees both survive forever. Intersecting the Biggins speed event simultaneously over all countably many Ulam–Harris subtree roots ensures that both surviving subtrees rooted at $`w_1`$ and $`w_2`$ achieve speed $`1/\beta_*`$ (up to finite birth-time offsets), producing for every large $`t`$ vertices $`u_1(t), u_2(t)`$ in disjoint subtrees below their LCA at distance $`\ge 2(1/\beta_* - \varepsilon_m)t`$, giving $`\liminf_{t\to\infty}\mathrm{diam}(T_t)/t \ge 2/\beta_* > 2/\beta_{\mathrm{hi}}`$ a.s. on $`\mathcal S_{\mathrm{surv}}`$.

To upgrade from $`\mathcal S_{\mathrm{surv}}`$ (probability $`\rho > 0`$) to **probability 1**:
1. On the extinction event of any coupled trial, the total progeny is finite and every individual stops after at most $`J_{\max}`$ rings of its own rate-1 clock, so the exact extinction instant $`T_{\mathrm{ext}}`$ (first time zero active individuals remain in the trial) is an a.s. finite $`\mathcal F_t`$-stopping time.
2. In **every** finite tree $`T`$ with $`C(T) \ge 1`$, the total core promotion rate obeys $`(\mathcal L C)(T) = \sum_{u \in K(T)} L_u / D_u \ge 1`$ (if $`C(T) = 1`$, $`L_o/D_o = 1`$; if $`C(T) \ge 2`$, the core tree has at least two leaves $`u`$ with $`c_u = 1`$, each contributing $`L_u/D_u = 1 - 1/D_u \ge 1/2`$). Hence after any finite stopping time $`T_{\mathrm{ext}}`$, the waiting time until the next core promotion is stochastically dominated by $`\mathrm{Exp}(1)`$ and therefore finite a.s.
3. Starting an independent fresh coupled CMJ trial from the next newly promoted $`(2, 1)`$ core after each extinction succeeds in finitely many trials with probability $`1`$ (since $`(1 - \rho)^m \to 0`$), yielding $`\liminf_{t\to\infty} H_{\mathrm{core}}(T_t)/t \ge 1/\beta_{\mathrm{hi}}`$ and $`\liminf_{t\to\infty} \mathrm{diam}(T_t)/t \ge 2/\beta_{\mathrm{hi}}`$ **almost surely**.

Finally, in continuous time, $`N(t) = |T_t|`$ is a standard rate-1 Yule process started at $`N(0) = 3`$, so $`e^{-t} N(t) \to W_3 \sim \mathrm{Gamma}(3, 1) \in (0, \infty)`$ a.s., whence the jump times $`\tau_n := \inf\{t \ge 0 : |T_t| = n\}`$ satisfy $`\tau_n / \log n \to 1`$ a.s. Evaluating along $`t = \tau_n`$ (since $`T_{\tau_n} = T_n`$) and using $`H(T_n) \ge H_{\mathrm{core}}(T_n) - 1`$ proves the lower bounds in (1) and (2) with $`1/\beta_{\mathrm{hi}} = 5000/5067`$ and $`2/\beta_{\mathrm{hi}} = 10000/5067`$. $`\square`$

---

## Exact Replay Command

```sh
python3 -OB round2/proofs/neon2/height-diameter/cert_height.py
```

Executes in ~2 seconds using zero floating-point arithmetic in certification (byte-identical under `-B` and `-OB`, fail-closed with explicit runtime checks and strict integer tick validation), checking all `5,148` upper supermartingale state rows (`2` rate extremes each) + far-tail bound, all `7,500` CMJ state rows (`25,000` recurrence inequalities) at $`s \in \{0, \theta_1, \theta_2\}`$ + exact Taylor convexity certificate, all `23` unlabelled trees $`n = 3..7`$ (`59` rooted core configurations) against the direct physical arrow generator, and `11/11` planted negative controls.
