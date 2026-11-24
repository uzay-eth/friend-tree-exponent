# All-Real-Moment, All-Label Strict Monotonicity in the Friend Tree (Open Question 3 for ALL $`p \in [1, \infty)`$)

**Pinned claim commit:** `fa615a179f65b34501f87f204605f8c7a7174c31`.  
**Review and dependencies:** GPT-family table, hazard and tail audits; Claude rederived the analytic skeleton. No inherited proof dependencies.

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-9.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-9.md).

**Author:** Neon-2 lead (Gemini 4 Neon).  
**Dependencies:** **NONE** — 100% self-contained from the definition of the friend-tree process (zero Astra dependencies, zero asymptotics, zero floating-point arithmetic, zero Gamma approximations).

---

## Headline Theorems

Let $`(T_n)_{n \ge 2}`$ be the discrete-time friend tree started from any finite tree $`T`$ of size $`n \ge 2`$ (and in particular the original edge-seeded tree $`T_2 = \{1, 2\}`$), with limiting degree shares $`Z_v := \lim_{m \to \infty} D_v(m)/m \in [0, 1]`$ (which exist a.s. and in $`L^p`$ for every $`p \ge 1`$ since $`D_v(m)/m \in [0, 1]`$ is a bounded nonnegative supermartingale).

### Theorem 1 (Universal Per-Tree Next-Two-Birth $`p`$-th Moment Ordering — Every Tree, Every Size $`n \ge 2`$, Every Real $`p \in [1, \infty)`$)
For **every** real number $`p \in [1, \infty)`$, **every** integer $`n \ge 2`$, and **every** finite tree $`T`$ of size $`n`$,
```math
\mathbb E\big[Z_{n+1}^p - Z_{n+2}^p \mid T_n = T\big]
\;\ge\; \frac{\Gamma(n + 2)\,\Gamma(p + 1)}{n(n + 1)\,\Gamma(n + 2 + p)}\, c_*(p)
\;>\; 0,
\qquad\text{(1)}
```
where $`c_*(p) > 0`$ (defined as the maximum of the applicable certified lower bounds at interval junctions and $`p = 12`$, and the analytic expression below for $`p > 12`$) satisfies the explicit uniform lower bounds:
- $`c_*(p) \ge \frac{17903587}{2500000000} > \frac{1}{140}`$ for all $`p \in [1, 12]`$ (certified per interval across 19 contiguous rational intervals in `cert_all_moments.py`);
- $`c_*(p) \ge \frac{1}{3} - \frac{4\sqrt{p + 1}}{\sqrt{3}\,(2p + 1)} > \frac{1}{4000}`$ for all $`p \in [12, \infty)`$ (with $`c_*(p) \to \frac{1}{3}`$ as $`p \to \infty`$);
- in particular, $`c_*(p) > \frac{1}{4000}`$ uniformly over **all** real $`p \in [1, \infty)`$.

### Theorem 2 (All-Label, All-$`L^p`$-Norm Strict Monotonicity on $`T_2 = \{1, 2\}`$)
For the original edge-seeded friend tree $`T_2 = \{1, 2\}`$ and **every** real number $`p \in [1, \infty)`$ (in particular, every positive integer moment $`k \in \mathbb Z_{\ge 1}`$ and every $`L^p`$ norm $`\|Z_u\|_p = (\mathbb E[Z_u^p])^{1/p}`$),
```math
\mathbb E[Z_1^p] \;=\; \mathbb E[Z_2^p] \;>\; \mathbb E[Z_3^p] \;>\; \mathbb E[Z_4^p] \;>\; \mathbb E[Z_5^p] \;>\; \cdots,
\qquad\text{(2)}
```
with the closed-form zero-computer early gap
```math
\mathbb E[Z_2^p] - \mathbb E[Z_3^p]
\;\ge\; \frac{1}{(p + 1)(p + 2)}\!\left(1 - \frac{2\sqrt{p + 1}}{\sqrt{3}\,(2p + 1)}\right)
\;>\; \frac{2}{5\,(p + 1)(p + 2)}
\qquad (\forall\, p \in [1, \infty)),
\qquad\text{(3)}
```
and for every $`u \ge 3`$,
```math
\mathbb E[Z_u^p] - \mathbb E[Z_{u+1}^p]
\;\ge\; \frac{\Gamma(u + 1)\,\Gamma(p + 1)}{(u - 1)u\,\Gamma(u + 1 + p)}\, c_*(p)
\;>\; 0.
\qquad\text{(4)}
```

---

## Complete Proof

### Lemma 1 (Yule-Process Embedding & Scale-Invariant $`p`$-th Moment Amplitude $`\widetilde M_p`$)
Let $`(T_t)_{t \ge 0}`$ be the continuous-time friend tree started from any finite tree $`T'`$ of size $`N \ge 2`$ (each vertex carries an independent rate-1 Poisson clock; on a ring it chooses a uniform neighbor and attaches a new leaf to it). Then $`N(t) := |T_t|`$ is a rate-1 Yule process started at $`N(0) = N`$, independent of the embedded discrete chain $`(T_{\tau_m})_{m \ge N}`$, and $`e^{-t} N(t) \to W_N \sim \mathrm{Gamma}(N, 1)`$ a.s. For **every** real $`p \in [1, \infty)`$:
1. For a single rate-1 Yule particle $`N_1(t) \sim \mathrm{Geom}(e^{-t})`$ on $`\{1, 2, \dots\}`$ and any integer $`m \ge 1`$,
   $`\displaystyle \mathbb E\big[N_1(t)^{\overline{m}}\big] \;=\; \sum_{j=1}^\infty e^{-t}(1 - e^{-t})^{j-1}\frac{(j + m - 1)!}{(j - 1)!} \;=\; m!\, e^{mt},`$
   so $`\|e^{-t} N_1(t)\|_{L^m} \le (m!)^{1/m}`$ uniformly in $`t \ge 0`$. Since $`N(t)`$ is a sum of $`N`$ independent copies of $`N_1(t)`$, Minkowski's inequality gives $`\|e^{-t} N(t)\|_{L^m} \le N(m!)^{1/m}`$ uniformly in $`t \ge 0`$. Choosing $`m := \lceil p \rceil + 1 > p`$, the family $`\{(e^{-t} N(t))^p : t \ge 0\}`$ is bounded in $`L^{m/p}`$, hence uniformly integrable.
2. Since $`0 \le e^{-pt} D_v(t)^p = (D_v(t)/N(t))^p (e^{-t} N(t))^p \le (e^{-t} N(t))^p`$ and $`D_v(t)/N(t) \to Z_v`$ a.s. independently of $`W_N`$, Vitali's convergence theorem gives $`e^{-pt} D_v(t)^p \to W_N^p Z_v^p`$ a.s. and in $`L^1`$. Defining the normalized scale-invariant $`p`$-th moment amplitude
   $`\displaystyle \widetilde M_p(T', v) \;:=\; \frac{1}{\Gamma(p + 1)}\lim_{t \to \infty} e^{-pt}\,\mathbb E_{T'}\big[D_v(t)^p\big], \qquad\text{(5)}`$
   and using $`\mathbb E[W_N^p] = \Gamma(N + p)/\Gamma(N)`$, we obtain the exact identity
   $`\displaystyle \widetilde M_p(T', v) \;=\; \frac{\Gamma(N + p)}{\Gamma(N)\,\Gamma(p + 1)}\,\mathbb E_{T'}\big[Z_v^p\big] \;\ge\; 0, \qquad\text{(6)}`$
   and $`e^{-pt} \widetilde M_p(T_t, v) = \mathbb E_{T'}[W_N^p Z_v^p / \Gamma(p + 1) \mid \mathcal F_t]`$ is a closed uniformly integrable Doob martingale on $`[0, \infty]`$ (working in the usual augmentation of the filtration).

### Lemma 2 (Exact Rising-Gamma Submartingale, Supermartingale, and First-Ring Coupling)
For every real $`p \in [1, \infty)`$, define piecewise for any core vertex $`v \in K(T')`$ (degree $`D_v \ge 2`$, leaf-neighbor count $`L_v \in \{0, \dots, D_v\}`$, core-neighbor count $`c_v = D_v - L_v`$):
```math
\widetilde f_p(L) \;:=\; \begin{cases} 0, & L \in \{0, 1\}, \\ \dfrac{\Gamma(L - 1 + p)}{\Gamma(L - 1)\,\Gamma(p + 1)}, & L \ge 2, \end{cases}
\qquad
\widetilde W_p(D, L) \;:=\; \frac{\Gamma\big((D + L)/2 + p\big)}{\Gamma\big((D + L)/2\big)\,\Gamma(p + 1)}.
\qquad\text{(7)}
```
Then for **every** real $`p \in [1, \infty)`$, **every** finite tree $`T'`$, and **every** core vertex $`v \in K(T')`$:
1. **(Exact global sandwich):**
   $`\displaystyle \widetilde f_p(L_v) \;\le\; \widetilde M_p(T', v) \;\le\; \widetilde W_p(D_v, L_v). \qquad\text{(8)}`$
2. **(Exact three-line first-ring lower bound at $`L_v \ge 1`$):** Whenever $`v \in K(T')`$ has $`L_v \ge 1`$ (in particular $`(D_v, L_v) = (2, 1)`$),
   $`\displaystyle \widetilde M_p(T', v) \;\ge\; \frac{1}{3}. \qquad\text{(9)}`$
3. **(Exact predictable-hazard leaf upper bound):** For any leaf $`a`$ attached to a core vertex $`x`$ of degree $`d \ge 2`$ in $`T'`$, letting $`\widetilde w_p^* := \sup_{T'', u : (D_u, L_u) = (2, 1)} \widetilde M_p(T'', u)`$,
   $`\displaystyle \widetilde M_p(T', a) \;\le\; \frac{\widetilde w_p^*}{p d + 1} \;\le\; \frac{\widetilde W_p(2, 1)}{p d + 1} \;<\; \frac{2\sqrt{p + 1}}{\sqrt{3}\,(p d + 1)}. \qquad\text{(10)}`$

**Proof.**
Write $`x_v := (D_v + L_v)/2 = L_v + c_v/2 \ge 1`$. On any leaf gain at $`v`$ (rate $`r_v \in [L_v, x_v]`$), $`(D_v, L_v) \to (D_v + 1, L_v + 1)`$ and $`x_v \to x_v + 1`$; on any leaf promotion (rate $`L_v / D_v \le 1`$, occurring only when $`L_v \ge 1`$), $`(D_v, L_v) \to (D_v, L_v - 1)`$ and $`x_v \to x_v - 1/2 \ge 1/2`$.

*Upper bound in (8):* Since $`\Gamma(z + 1) = z\Gamma(z)`$, $`\widetilde W_p(D + 1, L + 1) - \widetilde W_p(D, L) = \frac{p}{x_v}\widetilde W_p(D, L)`$, so the gain term contributes $`r_v \frac{p}{x_v}\widetilde W_p \le p\,\widetilde W_p`$. Meanwhile $`y \mapsto \Gamma(y + p)/\Gamma(y)`$ is strictly increasing on $`(0, \infty)`$ (log-derivative $`\psi(y + p) - \psi(y) > 0`$), so promotion jumps are strictly negative. Hence $`(\mathcal L \widetilde W_p)(D_v, L_v) \le p\,\widetilde W_p(D_v, L_v)`$ on every core state. Because $`\widetilde W_p(D, L) \le \frac{\Gamma(D + p)}{\Gamma(D)\Gamma(p+1)} \le (D + p)^p / \Gamma(p+1)`$ (by log-convexity of $`\Gamma`$: writing $`p = m + \theta`$ with $`m \in \mathbb Z_{\ge 0}, \theta \in [0, 1)`$, $`\Gamma(y + \theta) \le y^\theta \Gamma(y)`$ gives $`\Gamma(x + p)/\Gamma(x) \le (x + m)^p \le (x + p)^p`$), auxiliary size localization at $`\tau_{K'} = \inf\{t : N(t) \ge K'\}`$ (dominated at each fixed $`t`$ by the single integrable variable $`(N(t) + p)^p`$) shows $`e^{-pt}\widetilde W_p(D_v(t), L_v(t))`$ is a nonnegative supermartingale. Furthermore, for $`y > 1 - \theta`$, log-convexity at $`y = \theta(y + \theta - 1) + (1 - \theta)(y + \theta)`$ yields $`\Gamma(y + \theta)/\Gamma(y) \ge (y + \theta - 1)^\theta`$, hence $`(x - 1)^p \le \Gamma(x + p)/\Gamma(x) \le (x + p)^p`$ for all $`x \ge 1`$. Since $`c_v(t) \le c_v(0) + N_v^{\mathrm{self}}(t) = O(t)`$ a.s. (only rings of $`v`$'s own rate-1 clock promote its leaves), on $`\{Z_v > 0\}`$ we have $`e^{-t} x_v(t) = e^{-t}(D_v(t) - c_v(t)/2) \to W_N Z_v`$ a.s., so $`e^{-pt}\widetilde W_p(D_v(t), L_v(t)) \to W_N^p Z_v^p / \Gamma(p + 1)`$ a.s. (and also to $`0`$ on $`\{Z_v = 0\}`$). Fatou's lemma gives $`\widetilde M_p(T', v) \le \widetilde W_p(D_v, L_v)`$.

*Lower bound in (8):* Check $`(\mathcal L \widetilde f_p)(L)`$ at every $`L \ge 0`$:
- $`L = 0`$: $`\mathcal L \widetilde f_p = 0 = p\,\widetilde f_p(0)`$.
- $`L = 1`$: $`\widetilde f_p(1) = \widetilde f_p(0) = 0`$ and $`\widetilde f_p(2) = 1`$, so $`\mathcal L \widetilde f_p = r_v \ge 1 > 0 = p\,\widetilde f_p(1)`$.
- $`L = 2`$: $`\widetilde f_p(2) = 1`$, $`\widetilde f_p(3) = 1 + p`$, $`\widetilde f_p(1) = 0`$; gain rate $`r_v \ge 2`$ and promotion rate $`2/D_v \le 1`$ give $`\mathcal L \widetilde f_p \ge 2p - 1 \ge p = p\,\widetilde f_p(2)`$ since $`p \ge 1`$.
- $`L \ge 3`$: gain jump is $`\frac{p}{L - 1}\widetilde f_p(L)`$ at rate $`r_v \ge L`$, and promotion jump is $`-\frac{p}{L + p - 2}\widetilde f_p(L)`$ at rate $`L/D_v \le 1`$, so
  $`\displaystyle (\mathcal L \widetilde f_p)(L) \;\ge\; p\!\left(\frac{L}{L - 1} - \frac{1}{L + p - 2}\right)\widetilde f_p(L) \;=\; p\!\left(1 + \frac{1}{L - 1} - \frac{1}{L + p - 2}\right)\widetilde f_p(L) \;\ge\; p\,\widetilde f_p(L)`$
  since $`L + p - 2 \ge L - 1`$ for all $`p \ge 1`$.
Because $`0 \le \widetilde f_p(L_v(t)) \le (N(t) + p)^p/\Gamma(p+1)`$ with $`e^{-pt}(N(t) + p)^p`$ uniformly integrable across all $`t \ge 0`$ (Lemma 1(1)), size localization shows $`e^{-pt}\widetilde f_p(L_v(t))`$ is a uniformly integrable submartingale converging a.s. (by the same sandwich since $`L_v(t) = D_v(t) - c_v(t)`$ with $`c_v(t) = O(t)`$) to $`W_N^p Z_v^p / \Gamma(p + 1)`$. Taking $`t \to \infty`$ gives $`\widetilde f_p(L_v) \le \widetilde M_p(T', v)`$.

*Proof of (9):* Let $`\sigma := \inf\{t > 0 : v\text{'s own Poisson clock rings at } t\} \sim \mathrm{Exp}(1)`$. Pick one initial leaf $`\ell_0 \sim v`$ at $`t = 0`$, and track the family $`\mathcal Y_t`$ of leaves adjacent to $`v`$ generated recursively by rings of leaves in $`\mathcal Y`$ before $`\sigma`$ ($`\mathcal Y_0 = \{\ell_0\}`$). Because no leaf adjacent to $`v`$ can be promoted before $`v`$'s own clock rings at $`\sigma`$, every ring of $`u \in \mathcal Y_{t-}`$ ($`t < \sigma`$) attaches a new leaf to $`v`$; extending $`|\mathcal Y_t|`$ after $`\sigma`$ by independent rate-1 leaf clocks yields a standard rate-1 Yule process $`(Y(t))_{t \ge 0}`$ started at $`Y(0) = 1`$, **strictly independent of $`\sigma \sim \mathrm{Exp}(1)`$**, with $`L_v(t) \ge Y(t)`$ for all $`t < \sigma`$. At $`t = \sigma`$, $`v`$ promotes at most one leaf (and no leaf clock ties a.s.), so $`L_v(\sigma) \ge Y(\sigma) - 1`$ a.s. Applying optional stopping to the closed UI Doob martingale $`e^{-pt}\widetilde M_p(T_t, v)`$ at $`\sigma < \infty`$ and using (8) and monotonicity of $`\widetilde f_p`$:
```math
\widetilde M_p(T', v) \;=\; \mathbb E_{T'}\big[e^{-p\sigma}\widetilde M_p(T_\sigma, v)\big]
\;\ge\; \mathbb E\big[e^{-p\sigma}\widetilde f_p(Y(\sigma) - 1)\big].
```
For fixed $`t > 0`$ with $`q := e^{-t} \in (0, 1)`$, $`Y(t) \sim \mathrm{Geom}(q)`$ on $`\{1, 2, \dots\}`$, and by the negative binomial series $`\sum_{j=0}^\infty \binom{p + j}{j} z^j = (1 - z)^{-(p + 1)}`$ ($`|z| < 1`$):
```math
\mathbb E\big[\widetilde f_p(Y(t) - 1)\big]
\;=\; \sum_{m=3}^\infty q(1 - q)^{m-1}\frac{\Gamma(m - 2 + p)}{\Gamma(m - 2)\,\Gamma(p + 1)}
\;=\; q(1 - q)^2 \sum_{j=0}^\infty \binom{p + j}{j}(1 - q)^j
\;=\; (1 - q)^2 e^{pt}.
```
Multiplying by $`e^{-pt}`$ cancels the exponential factor identically, and integrating over $`\sigma \sim \mathrm{Exp}(1)`$ ($`q = e^{-t}`$) gives
```math
\widetilde M_p(T', v) \;\ge\; \int_0^\infty e^{-t}(1 - e^{-t})^2\,dt \;=\; \int_0^1 (1 - q)^2\,dq \;=\; \frac{1}{3}.
```

*Proof of (10):* Let $`\tau_a := \inf\{t \ge 0 : D_a(t) = 2\}`$. Before $`\tau_a`$, $`D_x(t-) \ge d`$, so the survival function $`H(t) := \mathbb P(t < \tau_a) = 1 - \int_0^t \mathbb E[\mathbf 1_{\{s < \tau_a\}}/D_x(s-)]\,ds`$ satisfies $`H'(t) \ge -H(t)/d`$ a.e., hence $`H(t) \ge e^{-t/d}`$ and $`\mathbb E[e^{-p\tau_a}] = 1 - \int_0^\infty p e^{-pt} H(t)\,dt \le 1 - \frac{p}{p + 1/d} = \frac{1}{pd + 1}`$. On $`\{\tau_a < \infty\}`$, $`a`$ enters state $`(2, 1)`$ (so $`\widetilde M_p(T_{\tau_a}, a) \le \widetilde w_p^* \le \widetilde W_p(2, 1)`$), while on $`\{\tau_a = \infty\}`$, $`Z_a = 0`$. Optional stopping gives the first two inequalities in (10). Finally, log-convexity of $`\Gamma`$ gives $`\Gamma(p + 3/2)^2 \le \Gamma(p + 1)\Gamma(p + 2) = (p + 1)\Gamma(p + 1)^2`$, and $`\Gamma(3/2) = \sqrt{\pi}/2 > \sqrt{3}/2`$ (as $`\pi > 3`$), yielding $`\widetilde W_p(2, 1) = \frac{\Gamma(p + 3/2)}{\Gamma(3/2)\Gamma(p + 1)} < \frac{2\sqrt{p + 1}}{\sqrt{3}}`$, proving (10). $`\square`$

---

### Lemma 3 (Signed Two-Birth Identity & Endpoint Split for Every Real $`p \ge 1`$)
Fix any integer $`n \ge 2`$, any finite tree $`T`$ of size $`n`$, and any real $`p \in [1, \infty)`$. Write $`r_T(x) := \sum_{u \sim x} 1/D_T(u)`$ for the attachment weight of $`x \in T`$. For each $`x \in T`$, let $`T_{x-a-b}`$ be the size-$`(n+2)`$ tree formed by attaching leaf $`a`$ to $`x`$ and then leaf $`b`$ to $`a`$, and set $`\widetilde\Delta_{\mathrm{path}, p}(T, x) := \widetilde M_p(T_{x-a-b}, a) - \widetilde M_p(T_{x-a-b}, b)`$. For each unordered edge $`\{x, y\} \in E(T)`$, let $`T_{x, y}`$ be the size-$`(n+2)`$ tree formed by attaching a leaf $`a_x`$ to $`x`$ and a leaf $`a_y`$ to $`y`$. Then:
```math
\frac{n(n + 1)\,\Gamma(n + 2 + p)}{\Gamma(n + 2)\,\Gamma(p + 1)}\,\mathbb E\big[Z_{n+1}^p - Z_{n+2}^p \mid T_n = T\big]
\;=\; \sum_{x \in T} \frac{r_T(x)}{D_x + 1}\,\widetilde\Delta_{\mathrm{path}, p}(T, x)
\;+\; \sum_{\{x, y\} \in E(T)} \!\left(\frac{r_T(y)}{D_y(D_y + 1)} - \frac{r_T(x)}{D_x(D_x + 1)}\right)\!\big(\widetilde M_p(T_{x, y}, a_x) - \widetilde M_p(T_{x, y}, a_y)\big).
\qquad\text{(11a)}
```
Consequently,
```math
\frac{n(n + 1)\,\Gamma(n + 2 + p)}{\Gamma(n + 2)\,\Gamma(p + 1)}\,\mathbb E\big[Z_{n+1}^p - Z_{n+2}^p \mid T_n = T\big]
\;\ge\; \sum_{x \in T} \frac{r_T(x)}{D_x + 1}\!\left[\widetilde\Delta_{\mathrm{path}, p}(T, x) - \frac{1}{D_x}\sum_{y \sim x} \widetilde M_p(T_{x, y}, a_x)\right],
\qquad\text{(11)}
```
and furthermore $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge 1`$.

**Proof.**
Label the next two births $`a = n + 1`$ and $`b = n + 2`$. Conditioning on $`T_{n+2}`$ and applying Lemma 1 eq. (6) at size $`N = n + 2`$ gives $`\frac{\Gamma(n+2+p)}{\Gamma(n+2)\Gamma(p+1)}\mathbb E[Z_a^p - Z_b^p \mid T_{n+2}] = \widetilde M_p(T_{n+2}, a) - \widetilde M_p(T_{n+2}, b)`$. At step $`n \to n + 1`$, vertex $`a`$ attaches to $`x \in T`$ with probability $`r_T(x)/n`$, yielding $`T' = T \cup \{x - a\}`$ where $`D_{T'}(x) = D_x + 1`$, $`D_{T'}(a) = 1`$, and $`D_{T'}(z) = D_z`$ for $`z \in T \setminus \{x\}`$. At step $`n + 1 \to n + 2`$, four disjoint cases exhaust the attachment site of $`b`$:
1. $`b`$ attaches to $`a`$ (possible only when $`x`$ calls $`a`$, probability $`\frac{1}{(n+1)(D_x+1)}`$): yields $`T_{x-a-b}`$ with scaled weight $`\frac{r_T(x)}{D_x + 1}`$ and signed payoff $`\widetilde\Delta_{\mathrm{path}, p}(T, x)`$.
2. $`b`$ attaches to $`x`$: $`a`$ and $`b`$ are twin leaves of $`x`$, so $`\widetilde M_p(T_{n+2}, a) - \widetilde M_p(T_{n+2}, b) = 0`$ by rooted-tree isomorphism invariance.
3. $`b`$ attaches to $`y \in T \setminus \{x\}`$ with $`y \not\sim x`$: $`r_{T'}(y) = r_T(y)`$, so both birth orders $`(a \to x, b \to y)`$ and $`(a \to y, b \to x)`$ have identical scaled weight $`r_T(x)r_T(y)`$ and opposite signed payoffs on the same unlabelled tree, cancelling identically.
4. $`b`$ attaches to $`y \sim x`$ in $`T`$: raising $`D_x`$ to $`D_x + 1`$ changes $`r_{T'}(y)`$ by $`\frac{1}{D_x + 1} - \frac{1}{D_x} = -\frac{1}{D_x(D_x + 1)}`$. Pairing $`(a \to x, b \to y)`$ with $`(a \to y, b \to x)`$ yields the exact net scaled weight
   $`\displaystyle r_T(x)\!\left(r_T(y) - \frac{1}{D_x(D_x + 1)}\right) - r_T(y)\!\left(r_T(x) - \frac{1}{D_y(D_y + 1)}\right) \;=\; \frac{r_T(y)}{D_y(D_y + 1)} - \frac{r_T(x)}{D_x(D_x + 1)}`$
   multiplying $`\widetilde M_p(T_{x, y}, a_x) - \widetilde M_p(T_{x, y}, a_y)`$, proving (11a).

Because $`\widetilde M_p \ge 0`$, expanding each unordered edge term in (11a) and dropping the two nonnegative cross terms $`\frac{r_T(y)}{D_y(D_y+1)}\widetilde M_p(T_{x, y}, a_x) \ge 0`$ and $`\frac{r_T(x)}{D_x(D_x+1)}\widetilde M_p(T_{x, y}, a_y) \ge 0`$ gives
```math
\left(\frac{r_T(y)}{D_y(D_y + 1)} - \frac{r_T(x)}{D_x(D_x + 1)}\right)\!\big(\widetilde M_p(T_{x, y}, a_x) - \widetilde M_p(T_{x, y}, a_y)\big)
\;\ge\; -\frac{r_T(x)}{D_x(D_x + 1)}\widetilde M_p(T_{x, y}, a_x) - \frac{r_T(y)}{D_y(D_y + 1)}\widetilde M_p(T_{x, y}, a_y).
```
Regrouping by vertex $`x \in T`$ (which has $`D_x`$ incident edges $`\{x, y\}`$) proves (11). Finally, in every size-$`n`$ tree ($`n \ge 2`$), $`D_x + 1 \le n`$ and $`\sum_{x \in T} r_T(x) = \sum_{u \in T}\sum_{x \sim u} 1/D_u = n`$, so $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge \frac{1}{n}\sum_{x \in T} r_T(x) = 1`$. $`\square`$

---

### Lemma 4 (Analytic Proof for All Real $`p \ge 12`$ and Zero-Float $`q_p(D)`$-Rescaled Proof for All Real $`p \in [1, 12]`$)

1. **Analytic regime $`p \in [12, \infty)`$ (zero computer needed):**
   For every $`x \in T`$ and every neighbor $`y \sim x`$, $`a`$ starts at $`(2, 1)`$ in $`T_{x-a-b}`$, $`b`$ is a leaf on degree $`2`$, and $`a_x`$ is a leaf on degree $`D_T(x) + 1 \ge 2`$. By Lemma 2 eqs. (9)–(10),
   $`\displaystyle \widetilde\Delta_{\mathrm{path}, p}(T, x) - \widetilde M_p(T_{x, y}, a_x) \;=\; \widetilde M_p(T_{x-a-b}, a) - \widetilde M_p(T_{x-a-b}, b) - \widetilde M_p(T_{x, y}, a_x) \;>\; \frac{1}{3} - \frac{4\sqrt{p + 1}}{\sqrt{3}\,(2p + 1)}.`$
   Because $`\frac{p + 1}{(2p + 1)^2}`$ is strictly decreasing on $`[1, \infty)`$ (derivative $`-\frac{2p + 3}{(2p + 1)^3} < 0`$) and at $`p = 12`$ we have $`\frac{16 \cdot 13}{3 \cdot 25^2} = \frac{208}{1875} < \left(\frac{1}{3} - \frac{1}{4000}\right)^2 = \frac{15976009}{144000000}`$ (exact integer check $`15974400 < 15976009`$, and $`4 \cdot 12^2 - 44 \cdot 12 - 47 = +1 > 0`$), every bracketed factor in (11) exceeds $`\frac{1}{3} - \frac{4\sqrt{p+1}}{\sqrt{3}(2p+1)} > \frac{1}{4000}`$ for **all** real $`p \ge 12`$; weighting by $`\sum_{x \in T} \frac{r_T(x)}{D_x + 1} \ge 1`$ proves Theorem 1 on $`[12, \infty)`$.

2. **Compact regime $`p \in [1, 12]`$ (19 contiguous rational intervals, $`q_p(D)`$-rescaled exact local martingales):**
   For any integer $`D \ge 2`$ and real $`p > 0`$, define the exact rising-Gamma polynomial
   $`\displaystyle q_p(D) \;:=\; \frac{\Gamma(D + p)}{\Gamma(D)\,\Gamma(p + 1)} \;=\; \prod_{j=1}^{D-1}\!\left(1 + \frac{p}{j}\right) \;>\; 0, \qquad \frac{q_p(D + 1)}{q_p(D)} \;=\; 1 + \frac{p}{D}.`$
   Note that for integers $`D_a, D_b \ge 2`$, the ratio $`q_p(D_b)/q_p(D_a)`$ is constant ($`=1`$) if $`D_b = D_a`$, nonincreasing in $`p`$ if $`D_b < D_a`$ ($`\prod_{j=D_b}^{D_a-1} \frac{j}{j+p}`$), and nondecreasing in $`p`$ if $`D_b > D_a`$ ($`\prod_{j=D_a}^{D_b-1} (1 + \frac{p}{j})`$); thus its exact maximum on $`[p_0, p_1]`$ is the exact rational product $`q_{\mathrm{upper}}(D_a, D_b; p_0, p_1)`$ evaluated at $`p_0`$ (for $`D_b < D_a`$) or $`p_1`$ (for $`D_b > D_a`$).

   On each closed rational interval $`[p_0, p_1] \subset [1, 12]`$, `cert_all_moments.py` constructs and independently verifies $`p`$-independent integer-tick tables (divided by $`S = 10^9`$) defining $`q_p(D)`$-rescaled local processes. Fix any real $`p \in [p_0, p_1]`$. In every local comparison below, the process is stopped at its exit time $`\tau`$ (evaluating the prescribed boundary/handoff table value on the finite exit jump, and $`0`$ on $`\{\tau = \infty\}`$). Because the pre-exit and immediate post-exit coordinate domain is finite, the stopped local value is uniformly bounded in absolute value for fixed $`p`$ (even for signed tables). Since $`(T_t)_{t \ge 0}`$ is nonexplosive, localizing at finite tree sizes $`N(t) = K'`$ and applying bounded convergence justifies the stopped discounted drift inequality on $`[0, t \wedge \tau]`$ for every deterministic $`t > 0`$, and letting $`t \to \infty`$ by bounded convergence yields $`\Phi(0) \le \mathbb E[\mathbf 1_{\{\tau < \infty\}} e^{-p\tau}\Phi(\tau)]`$ (for submartingales; reversed for supermartingales). Separately, optional sampling on Lemma 1's closed UI Doob martingale $`e^{-pt}\widetilde M_p(T_t, v)`$ (or the difference of two such martingales) gives the identical representation with terminal value $`\widetilde M_p(T_\tau, v)`$ on $`\{\tau < \infty\}`$ and $`W_N^p Z_v^p/\Gamma(p+1)`$ on $`\{\tau = \infty\}`$. In every case below, on $`\{\tau = \infty\}`$ the tracked vertex degree(s) stay below their finite cutoff for all $`t \ge 0`$ while $`N(t) \to \infty`$, forcing $`Z_v = 0`$ a.s. (and both $`Z_a = Z_b = 0`$ for difference bounds). Hence verifying the finite-exit inequality and interior drift sign proves the corresponding global amplitude bound in the following acyclic order:

   - **Step A ($`\alpha, \beta`$ — 1-vertex core bounds $`v_p(D, L) := \alpha(D, L)\,q_p(D)`$ and $`w_p(D, L) := \beta(D, L)\,q_p(D)`$, $`2 \le D \le D_M`$, $`0 \le L \le D`$):**
     Stop at $`\tau := \inf\{t : D_v(t) = D_M + 1\}`$. While $`t < \tau`$, $`(D_v, L_v)`$ jumps to $`(D_v+1, L_v+1)`$ at rate $`r_v \in [L_v, L_v + c_v/2]`$ and to $`(D_v, L_v-1)`$ at rate $`L_v/D_v`$. Dividing $`(\mathcal L v_p - p v_p)`$ (resp. $`\mathcal L w_p - p w_p`$) by $`q_p(D) > 0`$ gives
     $`\displaystyle r\!\left(\!\left(1 + \frac{p}{D}\right)\alpha(D + 1, L + 1) - \alpha(D, L)\right) + \frac{L}{D}\big(\alpha(D, L - 1) - \alpha(D, L)\big) - p\,\alpha(D, L),`$
     which is **bilinear (separately affine) in $`(r, p)`$**; checking the 4 box corners $`\{L, L + c/2\} \times \{p_0, p_1\}`$ certifies the drift for all $`p \in [p_0, p_1]`$ and all admissible $`r`$. At $`D = D_M + 1`$, the lower ratio $`\widetilde f_p(L)/q_p(D_M+1) = \prod_{j=L-1}^{D_M} \frac{j}{j + p}`$ ($`L \ge 2`$; $`0`$ for $`L \le 1`$) is minimized at $`p_1`$, while for the upper ratio, setting $`m := \lceil (D_M + 1 + L)/2 \rceil \le D_M + 1`$ and using monotonicity of $`y \mapsto \Gamma(y + p)/\Gamma(y)`$ gives
     $`\displaystyle \frac{\widetilde W_p(D_M + 1, L)}{q_p(D_M + 1)} \;\le\; \frac{\Gamma(m + p)/\Gamma(m)}{\Gamma(D_M + 1 + p)/\Gamma(D_M + 1)} \;=\; \prod_{j=m}^{D_M} \frac{j}{j + p} \;\le\; \prod_{j=m}^{D_M} \frac{j}{j + p_0},`$
     an exact rational product with zero Gamma functions and zero $`\pi`$. The verifier checks all boundaries, $`\beta(D, L) \ge 0`$ everywhere, and all interior corners, proving $`\alpha(D_v, L_v)\,q_p(D_v) \le \widetilde M_p(T', v) \le \beta(D_v, L_v)\,q_p(D_v)`$.

   - **Step B ($`\gamma`$ — leaf upper bound $`U_p(d) := \gamma(d)\,q_p(2) = \gamma(d)(1 + p)`$, $`2 \le d \le D_M + 1`$):**
     For any leaf $`a`$ attached to a core $`x`$ of degree $`d \ge 2`$, Lemma 2(3) and Step A give the global bound $`\widetilde M_p(T', a) \le \frac{\beta(2, 1)(1+p)}{p d + 1}`$ (using $`\beta(2, 1) \ge 0`$). For $`d \ge D_M + 1`$, this is at most $`\frac{\beta(2, 1)(1+p)}{p_0(D_M + 1) + 1} \le \gamma(D_M + 1)(1+p)`$ (verified). For $`2 \le d \le D_M`$, stop at $`\tau := \tau_a \wedge \inf\{t : D_x(t) = D_M + 1\}`$; before exit, $`D_x`$ increments at physical rate $`r_x \in [1, d]`$ (leaf $`a`$ alone contributes rate $`1`$) and $`a`$ promotes to $`(2, 1)`$ (bounded by $`\beta(2, 1)(1+p)`$) at rate $`1/d`$, with no tie a.s. Dividing the supermartingale drift by $`1 + p > 0`$ gives $`r_x(\gamma(d+1) - \gamma(d)) + \frac{1}{d}(\beta(2, 1) - \gamma(d)) - p\gamma(d)`$. Because the verifier checks $`\gamma(d) \ge 0`$, $`\gamma(d+1) \le \gamma(d)`$ (so $`r_x = 1`$ is worst case), and drift $`\le 0`$ at both $`p \in \{p_0, p_1\}`$, optional stopping proves $`\widetilde M_p(T', a) \le \gamma(\min(D_{T'}(x), D_M + 1))(1 + p)`$.

   - **Step C (Mode-specific coupling bounds and certified initial margins):**
     - **`1v` mode** (9 intervals tiling $`[3, 12]`$, $`D_M = 45`$): Uses Steps A–B directly; at every $`x \in T`$ and $`y \sim x`$, $`\widetilde\Delta_{\mathrm{path}, p}(T, x) - \widetilde M_p(T_{x, y}, a_x) \ge (\alpha(2, 1) - 2\gamma(2))(1 + p) > 0`$.
     - **`vx_only` mode** (4 intervals tiling $`[8/5, 3]`$, $`D_M = 40, D_{x,\max} = 7`$): For a core $`a`$ with an adjacent tracked core $`x \sim a`$ ($`0 \le L_a \le D_a - 1`$), define $`v_{x, p}(D_x, D_a, L_a) := \mathrm{vx}(D_x, D_a, L_a)\,q_p(D_a)`$ and stop at $`\tau := \inf\{t : D_a(t) = D_M + 1 \text{ or } D_x(t) = D_{x,\max} + 1\}`$. Before exit, $`a`$ gains a leaf at rate $`r_a \in [L_a + 1/D_x, L_a + 1/D_x + (D_a - L_a - 1)/2]`$ (multiplying $`q_p(D_a)`$ by $`1 + p/D_a`$), $`x`$ gains a leaf at rate $`r_x = 1/D_a + \sum_{u \sim x, u \ne a} 1/D_u \in [1/D_a, D_x]`$ (leaving $`q_p(D_a)`$ unchanged), and $`a`$'s leaves promote at rate $`L_a/D_a`$ (disjoint from $`x`$, which is core). Dividing by $`q_p(D_a) > 0`$ yields a triaffine drift in $`(r_a, r_x, p)`$, checked at all 8 box corners, with boundary $`\mathrm{vx} \le \alpha`$ at both exit faces (valid by Step A). Splitting $`D_T(x) \in \{1, 2, \ge 3\}`$ (so $`D'_x = D_T(x) + 1 \in \{2, 3, \ge 4\}`$) and using $`\gamma(d) \le \gamma(4)`$ for $`d \ge 4`$ gives bracketed margin $`\ge \min(\mathrm{vx}(2, 2, 1) - 2\gamma(2),\, \mathrm{vx}(3, 2, 1) - \gamma(2) - \gamma(3),\, \alpha(2, 1) - \gamma(2) - \gamma(4))(1 + p) > 0`$.
     - **`full` mode** (6 intervals tiling $`[1, 8/5]`$, $`D_M = 40, K \in \{12, 14\}, D_{x,\max} \in \{7, 8\}`$):
       1. *Phase-2 adjacent-core difference* $`\phi_{2, p}(D_a, L_a; D_b, L_b) := \psi_2(D_a, L_a; D_b, L_b)\,q_p(D_a)`$ ($`2 \le D_a, D_b \le K`$, $`0 \le L_a \le D_a - 1`$, $`0 \le L_b \le D_b - 1`$): stop at $`\tau := \inf\{t : D_a(t) = K + 1 \text{ or } D_b(t) = K + 1\}`$. Before exit, no single birth affects both cores; gain rates satisfy $`r_a \in [L_a + 1/D_b, L_a + 1/D_b + (D_a - L_a - 1)/2]`$ (only $`a`$-gain multiplies $`q_p(D_a)`$ by $`1 + p/D_a`$) and $`r_b \in [L_b + 1/D_a, L_b + 1/D_a + (D_b - L_b - 1)/2]`$, with leaf promotions at $`L_a/D_a`$ and $`L_b/D_b`$. At exit, $`\psi_2 \le \alpha(D_a, L_a) - \beta(D_b, L_b)\,q_{\mathrm{upper}}(D_a, D_b; p_0, p_1)`$ (valid by $`\beta \ge 0`$ and Step A). Checking all 8 box corners of $`(r_a, r_b, p)`$ proves $`\phi_{2, p} \le \widetilde M_p(T', a) - \widetilde M_p(T', b)`$.
       2. *Phase-1 core–leaf difference* $`\phi_{1, p}(D_a, L_a) := \eta_1(D_a, L_a)\,q_p(D_a)`$ ($`2 \le D_a \le K, 1 \le L_a \le D_a - 1`$, core $`a`$ having a distinct core neighbor and a distinguished leaf $`b \sim a`$): stop at $`\tau := \inf\{t : D_b(t) = 2 \text{ or } D_a(t) = K + 1\}`$. Before exit, $`a`$ gains a leaf at rate $`r_a \in [L_a, L_a + (D_a - L_a)/2]`$ (multiplying $`q_p(D_a)`$ by $`1 + p/D_a`$), other leaves of $`a`$ promote at rate $`(L_a - 1)/D_a`$, and $`b`$ promotes at rate $`1/D_a`$ directly into Phase 2 at $`(D_a, L_a - 1; 2, 1)`$ (leaving $`q_p(D_a)`$ unchanged, handed off to $`\psi_2(D_a, L_a - 1; 2, 1)`$). At $`D_a = K + 1`$, $`\eta_1 \le \alpha(K+1, L_a) - \gamma(K+1)\,q_{\mathrm{upper}}(K+1, 2; p_0, p_1)`$ (valid by $`\gamma \ge 0`$ and Steps A–B). Checking all 4 box corners of $`(r_a, p)`$ proves $`\widetilde\Delta_{\mathrm{path}, p}(T, x) \ge \eta_1(2, 1)(1 + p)`$ for all $`x \in T`$.
       3. *Phase-1 parent-tracked difference* $`\phi_{1, x, p}(D_x, D_a, L_a) := \eta_{1, x}(D_x, D_a, L_a)\,q_p(D_a)`$ ($`2 \le D_x \le D_{x,\max}, 2 \le D_a \le K, 1 \le L_a \le D_a - 1`$): stop at $`\tau := \inf\{t : D_b(t) = 2 \text{ or } D_a(t) = K + 1 \text{ or } D_x(t) = D_{x,\max} + 1\}`$, refining $`r_a \in [L_a + 1/D_x, L_a + 1/D_x + (D_a - L_a - 1)/2]`$ and tracking $`x`$-gain at rate $`r_x \in [1/D_a, D_x]`$, with $`b`$-promotion handoff to $`\psi_2(D_a, L_a - 1; 2, 1)`$, $`D_x`$-exit handoff $`\eta_{1, x} \le \eta_1(D_a, L_a)`$, and $`D_a`$-exit bound $`\eta_{1, x} \le \alpha(D_a, L_a) - \gamma(D_a)\,q_{\mathrm{upper}}(D_a, 2; p_0, p_1)`$. Checking all 8 box corners of $`(r_a, r_x, p)`$ and splitting $`D_T(x) \in \{1, 2, \ge 3\}`$ gives bracketed margin $`\ge \min(\eta_{1, x}(2, 2, 1) - \gamma(2),\, \eta_{1, x}(3, 2, 1) - \gamma(3),\, \eta_1(2, 1) - \gamma(4))(1 + p) > 0`$.

   Across all 19 contiguous intervals (`115,964` exact rational rows), every interval achieves certified $`c_*(p) \ge \eta_*(I)(1 + p_0) \ge \frac{17903587}{2500000000} > \frac{1}{140}`$. Together with Lemma 3 and part (1), this proves Theorem 1 for every real $`p \in [1, \infty)`$. $`\square`$

---

### Proof of Theorem 2
On $`T_2 = \{1, 2\}`$, $`\mathbb E[Z_1^p] = \mathbb E[Z_2^p]`$ by symmetry. After step $`2 \to 3`$, $`T_3`$ is deterministically the 3-vertex path where vertex 2 is center $`(D, L) = (2, 2)`$ with probability $`1/2`$ and a leaf (symmetric to vertex 3) with probability $`1/2`$, while vertex 3 is always a leaf on degree 2. By Lemma 1 eq. (6) (with $`\frac{\Gamma(3)\Gamma(p+1)}{\Gamma(3+p)} = \frac{2}{(p+1)(p+2)}`$) and Lemma 2 eqs. (8) & (10):
```math
\mathbb E[Z_2^p] - \mathbb E[Z_3^p]
\;=\; \frac{1}{2}\cdot\frac{2}{(p + 1)(p + 2)}\big(\widetilde M_p(T_3, \text{center}) - \widetilde M_p(T_3, \text{leaf})\big)
\;\ge\; \frac{1}{(p + 1)(p + 2)}\!\left(1 - \frac{2\sqrt{p + 1}}{\sqrt{3}\,(2p + 1)}\right),
```
since $`\widetilde f_p(2) = 1`$. Because $`\frac{p+1}{(2p+1)^2}`$ is strictly decreasing on $`[1, \infty)`$, at $`p = 1`$ we have $`\left(\frac{2\sqrt{2}}{3\sqrt{3}}\right)^2 = \frac{8}{27} < \frac{9}{25} = \left(\frac{3}{5}\right)^2`$, yielding $`1 - \frac{2\sqrt{p+1}}{\sqrt{3}(2p+1)} > \frac{2}{5}`$ and proving (3). Taking expectation of Theorem 1 at $`n = u - 1 \ge 2`$ gives (4) for all $`u \ge 3`$. $`\square`$

---

## Certified 19-Interval Cover of $`[1, 12]`$ (`cert_all_moments.py`)

| # | Interval $`[p_0, p_1]`$ | Mode | $`(D_M, K, D_{x,\max})`$ | Verified Rows | Certified $`\eta_*`$ | Certified $`c_*(p) \ge \eta_*(1+p_0)`$ |
|---:|---|---|---:|---:|---:|---:|
| 1 | $`[1, 21/20]`$ | `full` | $`(40, 14, 8)`$ | 13,299 | `0.012222` | `0.024445` |
| 2 | $`[21/20, 11/10]`$ | `full` | $`(40, 14, 8)`$ | 13,299 | `0.016038` | `0.032879` |
| 3 | $`[11/10, 6/5]`$ | `full` | $`(40, 14, 8)`$ | 13,299 | `0.003778` | `0.007934` |
| 4 | $`[6/5, 13/10]`$ | `full` | $`(40, 14, 8)`$ | 13,299 | `0.010594` | `0.023308` |
| 5 | $`[13/10, 29/20]`$ | `full` | $`(40, 12, 7)`$ | 8,146 | `0.003550` | `0.008166` |
| 6 | $`[29/20, 8/5]`$ | `full` | $`(40, 12, 7)`$ | 8,146 | `0.011188` | `0.027411` |
| 7 | $`[8/5, 9/5]`$ | `vx_only` | $`(40, -, 7)`$ | 6,669 | `0.002754` | `0.007161` |
| 8 | $`[9/5, 21/10]`$ | `vx_only` | $`(40, -, 7)`$ | 6,669 | `0.004820` | `0.013496` |
| 9 | $`[21/10, 13/5]`$ | `vx_only` | $`(40, -, 7)`$ | 6,669 | `0.004980` | `0.015439` |
| 10 | $`[13/5, 3]`$ | `vx_only` | $`(40, -, 7)`$ | 6,669 | `0.018950` | `0.068220` |
| 11–19 | $`[3, 12]`$ (9 steps) | `1v` | $`(45, -, -)`$ | $`9 \times 2{,}200`$ | `0.002344..0.009708` | `0.009377..0.075945` |

*(Decimal entries in the table above are 6-place downward truncations of the exact rational margins printed by `cert_all_moments.py`.)*

Total certified exact rational rows: **`115,964`** (plus **`[12, \infty)`** analytically in Lemma 4(1)).

---

## Exact Replay Command

```sh
python3 -OB round2/proofs/neon2/all-moments-order/cert_all_moments.py
```

Executes in ~25 seconds using zero floating-point arithmetic (identical output under `-B` and `-OB`, fail-closed with explicit runtime checks and strict integer tick validation), checking all `115,964` rational rows across all 19 contiguous intervals tiling $`[1, 12]`$, the exact analytic inequalities for $`p \ge 12`$ and $`u = 2`$ vs $`3`$, a finite-tree verification of the exact signed two-birth identity and endpoint split on all `24` unlabelled trees $`n = 2..7`$ (supplementing Lemma 3's universal proof), and `14/14` planted single-tick, boundary, non-integer, and malformed-tiling negative controls.
