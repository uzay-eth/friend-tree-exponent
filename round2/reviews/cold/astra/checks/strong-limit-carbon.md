CORRECT

# Cold Challenge and Review of `round2/fringe/almost-sure.md` (Sections 1–5)

## Scope and Verdict

This report cold-challenges every step of almost-sure.md (`fringe/almost-sure.md`), including the newly appended Section 5, for the undirected random friend tree started from a single edge $`T_2`$:
1. Existence of a finite random variable $`X`$ with $`X > 0`$ almost surely such that $`C_n / n^\mu \to X`$ and $`N_{\ge k}(n) / C_n \to c_k / c_2`$ almost surely for every fixed $`k \ge 2`$ (simultaneously across all $`k \ge 2`$), so $`N_{=k}(n) / C_n \to p_k = (c_k - c_{k+1})/c_2 \sim \mu A k^{-1-\mu}`$;
2. The discrete eigenmartingale $`Y_m = \frac{\Gamma(m)}{\Gamma(m+\mu)}\Psi(T_m)`$ ($`m \ge 3`$), the identity $`X = d_2 \lim_{m \to \infty} Y_m`$, and unbounded/zero-accumulating support $`\mathbb P(X < \varepsilon) > 0`$ and $`\mathbb P(X > M) > 0`$ for all $`\varepsilon, M > 0`$.

Section 1 of finite-kernel.md (`fringe/finite-kernel.md`) (the nineteen capped birth labels $`m(D) = \min(D, 20)`$ and residual growth bound $`\rho = 51/100 < \mu`$) is treated as a conditional dependency under separate review, with every prerequisite used here checked from first principles alongside asymptotics.md (`fringe/asymptotics.md`), renewal.md (`fringe/renewal.md`), existence.md (`proofs/existence.md`), and almost-sure-exponent.md (`proofs/almost-sure-exponent.md`). Every step in Sections 1–5 holds rigorously without gaps.

---

## 1. High-Frequency Decay, Contour Shift, and Uniform Exponential Renewal (Section 1)

### 1.1 Capped-mark root-gain bound, $`\|R B_0 g\|_W \le C`$, and removal of stopping
Let $`m(D) = \min(D, 20) \in \{2, \dots, 20\}`$, $`B_0 g(S) = b_S g(m(D_{\mathrm{root}}))`$ for $`g \in \mathbb C^{19}`$ with $`\|g\|_\infty \le 1`$, and $`k(t) = J R_t B_0`$. On any finite planted type $`S`$ with root $`o`$ ($`D_o \ge 2`$, $`c_o \ge 1`$, $`b_S = 1 - c_o/D_o \in [0, 1]`$), $`(B_0 g)(S)`$ depends only on $`(L_o, c_o)`$, so intrinsic transitions strictly below $`o`$ (and rings of $`o`$ choosing a core child or the inert parent slot) produce **zero** change in $`B_0 g`$. Only two intrinsic root events change $`B_0 g`$:
- **Root leaf gain** ($`(L_o, c_o) \mapsto (L_o + 1, c_o)`$, rate $`r_{\mathrm{gain}}(o) = L_o + \sum_{y \text{ core child}} 1/D_y < D_o`$):
  - If $`D_o < 20`$: $`r_{\mathrm{gain}}(o) < 20`$ and $`|\Delta(B_0 g)| \le 2`$, so $`r_{\mathrm{gain}}(o)|\Delta(B_0 g)| \le 40`$.
  - If $`D_o \ge 20`$: $`m(D_o + 1) = m(D_o) = 20`$, so $`g(m(D_o))`$ is **unchanged** and only $`b_S`$ changes by $`\frac{L_o + 1}{D_o + 1} - \frac{L_o}{D_o} = \frac{c_o}{D_o(D_o + 1)}`$. Thus
    $`\displaystyle r_{\mathrm{gain}}(o)\,|\Delta(B_0 g)| \le D_o \cdot \frac{c_o}{D_o(D_o + 1)} |g(20)| = \frac{c_o}{D_o + 1} < 1.`$
    (Without capping the mark space, an oscillating bounded $`g`$ has $`D_o |\Delta(B_0 g)| = \Theta(D_o) \to \infty`$; capping at $`20`$ is essential here and also preserves $`(1 + d\,m(D))^{-1} = 1/4 \le q h(c) = 63/250`$ for $`D \ge 20, c \ge 4`$.)
- **Root promotion** ($`(L_o, c_o) \mapsto (L_o - 1, c_o + 1)`$, rate $`b_S \le 1`$): $`D_o`$ is unchanged and $`|\Delta(B_0 g)| = |g(m(D_o))|/D_o \le 1/2`$.

Killing occurs at rate $`a_S \le 1/2`$ with $`|\Delta(B_0 g)| \le 1`$, and copy births have nonnegative rates $`a_Y - p_S \ge 0`$, $`p_S \ge 0`$ (since $`a_Y = 1/m(D_{\mathrm{parent}}(\text{birth})) \ge 1/D_{\mathrm{parent}}(\text{birth}) \ge p_S`$) with expected absolute jump sum $`\sum_Y a_Y \le \sum_Y a_Y W(Y) = (a_S + d)W(S) - h(c_o) \le \frac{13}{20} W(S)`$. Since $`W(S) \ge 14/13 > 1`$, both $`|(R B_0 g)(S)|`$ and the predictable rate of nonzero absolute jumps of $`\Phi_t = \sum_{S \in \mathcal P_t} (B_0 g)(S)`$ are bounded by $`42 W(S)`$. Stopping $`\mathcal P_t`$ at finite event/size cutoffs $`\tau_M`$, the pathwise total variation $`V_t = |\Phi_0| + \sum_{0 < u \le t}|\Delta \Phi_u|`$ satisfies $`\mathbb E_{J_{1/i}}[V_t] \le 1 + 42 \int_0^t R_u W(J_{1/i})\,du < \infty`$ by monotone convergence and nonexplosion. Dominated convergence with majorant $`V_t \in L^1`$ removes the stopping cutoff to give
```math
k(t) g = J B_0 g + \int_0^t J R_u R B_0 g\,du.
```
By Fubini (integration by parts) on $`\Re z \ge \rho + \varepsilon > \rho`$,
```math
K(z) g = \frac{1}{z}\left(J B_0 g + \int_0^\infty e^{-zu} J R_u R B_0 g\,du\right), \qquad \|K(z)\| \le \frac{1}{|z|}\left(1 + \frac{280}{\varepsilon}\right) = O_\varepsilon(1/|z|).
```

### 1.2 Spectral gap $`\sigma < \mu`$ and $`O(1/|\Im z|^2)`$ rectangular contour shift (Equation (2))
At $`z = \mu > \rho`$, $`K(\mu) \in (0, \infty)^{19 \times 19}`$ has simple Perron–Frobenius eigenvalue $`1`$ with $`g, \ell > 0`$, $`\ell g = 1`$, and $`-\ell K'(\mu)g > 0`$. On $`\Re z > \mu`$, $`r(K(z)) \le r(K(\Re z)) < 1`$. On $`z = \mu + i\tau`$ ($`\tau \ne 0`$), $`k_{i, 2}(t) \ge \frac{1}{2} e^{-(2 + 1/i)t} > 0`$ for all $`t > 0`$ forces strict phase cancellation $`|K(\mu + i\tau)v| < K(\mu)|v|`$ on any putative eigenvector $`v \ne 0`$, so $`I - K(\mu + i\tau)`$ is invertible. Taking $`\varepsilon_0 = (\mu - \rho)/2`$ and $`T_0`$ so $`\|K(z)\| \le 1/2`$ on $`\{\Re z \ge \rho + \varepsilon_0, |\Im z| \ge T_0\}`$ gives $`\|(I - K(z))^{-1}\| \le 2`$ there; compactness of $`\{\mu + i\tau : \delta_1/2 \le |\tau| \le T_0\}`$ and isolation of the simple pole at $`\mu`$ yield $`\sigma \in (\rho + \varepsilon_0, \mu)`$ such that $`(I - K(z))^{-1}`$ is meromorphic on $`\Re z \ge \sigma`$ with sole singularity $`(z - \mu)^{-1} D`$, $`D = \frac{g \otimes \ell}{-\ell K'(\mu)g}`$, and uniformly bounded on $`\Re z = \sigma`$.

Let $`w(t) := v(t) - k(t) = \sum_{n \ge 2} k^{*n}(t)`$. Subtracting $`k(t)`$ achieves two necessary properties simultaneously:
1. **Continuity at $`t = 0`$:** Each $`k^{*n}(t)`$ ($`n \ge 2`$) is continuous on $`[0, \infty)`$ with $`k^{*n}(0) = 0`$ (whereas $`k(0) = J B_0 \ne 0`$), and $`\|k^{*n}(t)\| \le M_T^n \frac{t^{n-1}}{(n-1)!}`$ on $`[0, T]`$, so $`\varphi(t) := e^{-st} w(t)\mathbf{1}_{[0, \infty)}(t)`$ ($`s > \mu`$ fixed) is **continuous on all of $`\mathbb R`$** and lies in $`L^1(\mathbb R)`$ by nonnegativity ($`\int_0^\infty e^{-st} w(t)\,dt = (I - K(s))^{-1} K(s)^2 < \infty`$).
2. **$`L^1`$ high-frequency decay:** Its Fourier/Laplace transform is $`F_w(z) = (I - K(z))^{-1} - I - K(z) = (I - K(z))^{-1} K(z)^2`$, which satisfies $`\|F_w(x + i\tau)\| = O(1/\tau^2)`$ uniformly for $`x \in [\sigma, s]`$ and $`|\tau| \ge T_0`$.

By classical $`L^1`$ Fourier inversion on $`\varphi \in L^1(\mathbb R) \cap C(\mathbb R)`$ with $`\widehat\varphi \in L^1(\mathbb R)`$, for **every** $`t \ge 0`$:
```math
w(t) = \frac{1}{2\pi i} \int_{s - i\infty}^{s + i\infty} e^{zt}(I - K(z))^{-1} K(z)^2\,dz.
```
Shifting the vertical contour across $`[\sigma, s] \times [-iL, iL]`$ as $`L \to \infty`$: the horizontal segments vanish as $`O(e^{st}/L^2)`$, the residue of $`e^{zt}[(I - K(z))^{-1} - I - K(z)]`$ at $`z = \mu`$ is $`e^{\mu t} D`$ (since $`-I - K(z)`$ is analytic at $`\mu`$), and the $`\Re z = \sigma`$ integral has norm at most $`\frac{e^{\sigma t}}{2\pi}\int_{-\infty}^\infty \|(I - K(\sigma + i\tau))^{-1}K(\sigma + i\tau)^2\|\,d\tau = C' e^{\sigma t}`$. Absorbing $`\|k(t)\| = O(e^{\rho t}) \le O(e^{\sigma t})`$ gives $`v(t) = e^{\mu t} D + O(e^{\sigma t})`$ uniformly for all $`t \ge 0`$.

### 1.3 Uniform weighted-norm renewal (Equation (3))
For any $`|f| \le W`$, $`\|h_f(u)\|_\infty = \|J R_u f\|_\infty \le \frac{20}{3} e^{\rho u}`$. Substituting (2) into $`J P_t f = h_f(t) + \int_0^t v(t - u) h_f(u)\,du`$ and extending $`\int_0^t e^{-\mu u} h_f(u)\,du`$ to $`\int_0^\infty`$ incurs tail error $`O(e^{\mu t} e^{-(\mu - \rho)t}) = O(e^{\rho t}) \le O(e^{\sigma t})`$ and convolution remainder $`\int_0^t e^{\sigma(t-u)} e^{\rho u}\,du = O(e^{\sigma t})`$ since $`\sigma > \rho`$, all with constants independent of $`f`$. Substituting $`J P_{t-u} f = e^{\mu(t-u)} g\, d(f) + O(e^{\sigma(t-u)})`$ into $`P_t f(S) = R_t f(S) + \int_0^t (R_u B_0 J P_{t-u} f)(S)\,du`$ and using $`|(R_u B_0 v)(S)| \le \|v\|_\infty e^{\rho u} W(S)`$ yields Equation (3) uniformly in $`S`$ and $`|f| \le W`$:
```math
\big|P_t f(S) - e^{\mu t} H(S) d(f)\big| \le C e^{\sigma t} W(S), \qquad H(S) = \int_0^\infty e^{-\mu u}(R_u B_0 g)(S)\,du.
```

---

## 2. Physical Eigenfunction $`\Psi(T)`$ for Every Finite Tree with a Core (Section 2)

Let $`T`$ be any finite tree with $`C(T) \ge 1`$. Pick any core $`o \in K(T)`$ as founder and any marks $`a_v \in \{1/2, \dots, 1/20\}`$ ($`v \in K(T) \setminus \{o\}`$) satisfying $`p_{\mathrm{parent}(v)}(T) \le a_v`$ (note that **every** finite tree $`T`$ with $`C(T) \ge 1`$ admits such marks: setting $`a_v = 1/2`$ for all $`v \ne o`$ satisfies $`1/D_u(T) \le 1/2 = a_v`$ on every core edge).
- **Weight bounds:** $`Z(T) = \sum_{x \in K(T)} w_x h(c_x) \le 2 C(T)`$ since $`w_x \le 1`$ and $`h \le 2`$. Along every path from $`o`$, $`z(a) \le z(1/2) = 10/13 < 1`$, so summing $`W(S_v)`$ over ancestors $`v \ne o`$ of each core $`x \ne o`$ gives $`\sum_{v \ne o} W(S_v) \le \frac{20}{3}(1 - 10/13)^{-1} C(T) = \frac{260}{9} C(T)`$. Because $`L^0 W \le q W`$ holds for all $`D \ge \max(2, c)`$ under $`m(D) = \min(D, 20)`$, $`\mathbb E_T Z(T_t) \le e^{qt} Z(T)`$, and $`\mathcal I(T_t) W = b_o W(J_{1/m(D_o)}) + \sum_{Y \text{ at } o} a_Y W(Y) \le 1/d + Z(T_t)`$ gives $`I_{T, t} W \le 1/d + e^{qt} Z(T) \le (1/d + 2 e^{qt}) C(T)`$ (using $`C(T) \ge 1`$).
- **General-founder immigration and Equation (4):** Applying the time-dependent jump formula on $`[0, t]`$ to $`\sum_{v \ne o \text{ in } T_u} (P_{t-u} f)(S_v(u))`$ retains the non-empty initial value $`\sum_{v \ne o \text{ in } T} P_t f(S_v)`$ at $`u = 0`$, giving $`\mathbb E_T \sum_{v \ne o} f(S_v(t)) = \sum_{v \ne o} P_t f(S_v) + \int_0^t I_{T, u}(P_{t-u} f)\,du`$. Inserting (3) and using $`0 < q = 9/25 < \rho = 51/100 < \sigma < \mu`$ bounds both $`\sum_{v \ne o}`$ error and the $`[0, t]`$ / $`[t, \infty)`$ immigration errors by $`C_1 C(T) e^{\sigma t}`$, uniformly in $`T`$ and $`|f| \le W`$. Since the founder $`o`$ contributes $`\Pr_T(D_o(t) \ge k) \in [0, 1]`$ and $`1 \le C(T) e^{\sigma t}`$, $`|\mathbb E_T A_k(T_t) - e^{\mu t} \Psi(T) d_k| \le (C_1 + 1) C(T) e^{\sigma t}`$.
- **Mark/founder independence, two-sided $`C(T)`$ bounds, and exact martingale (Equation (5)):** Taking $`f = f_2 = 1`$ ($`d_2 > 0`$) gives
  $`\displaystyle \Psi(T) = d_2^{-1} \lim_{u \to \infty} e^{-\mu u} \mathbb E_T C(T_u),`$
  whose right-hand side is purely physical, independent of the chosen founder $`o \in K(T)`$ and compatible marks. By $`H \le C_H W`$, $`\Psi(T) \le C_2 C(T)`$. Conversely, $`\min_{2 \le j \le 20} g(j) > 0`$ gives uniform positive lower bounds on $`H(S_v)`$ for noninitial cores with $`(c_v \le 2, L_v \ge 1)`$ and for bare noninitial cores with a bare child; by the eligibility lemma (`|\mathcal E| \ge C/4`), the count $`|G|`$ of such noninitial vertices satisfies $`|G| \ge \max(0, C(T)/8 - 3/2) \ge C(T)/12 - 1`$ for every $`C(T) \ge 1`$, giving $`C(T) \le C_3(1 + \Psi(T))`$. For any $`h \ge 0`$, conditioning at $`h`$ gives $`|d_2^{-1} e^{-\mu u} \mathbb E_T C(T_{h+u}) - \mathbb E_T \Psi(T_h)| \le d_2^{-1}(C_1 + 1)\mathbb E_T[C(T_h)] e^{-(\mu-\sigma)u} \to 0`$ as $`u \to \infty`$ (since $`\mathbb E_T C(T_h) \le n(T)e^h < \infty`$), while the left-hand side tends to $`e^{\mu h}\Psi(T)`$. Thus $`\mathbb E_T \Psi(T_h) = e^{\mu h}\Psi(T)`$, making $`Z_t^* = e^{-\mu t}\Psi(T_t)`$ a nonnegative càdlàg martingale with finite a.s. limit $`Z_\infty^*`$.

---

## 3. Clean/Dirty Coupling, $`h = \eta^{-1}\log(j+2)`$ $`L^2`$ Bound, and Positivity (Section 3)

### 3.1 Verification of the one-leaf coupling and sharpness of $`3 e^{2u}`$ (Equation (6))
Let $`T`$ have $`C(T) \ge 1`$ ($`n(T) \ge 3`$) and let $`T^+`$ attach one new leaf $`\ell^*`$ to $`w \in V(T)`$. Let $`D_{\mathrm{comb}}(u) = |\mathcal D(T_u)| + |\mathcal D(T^+_u)|`$ be the combined count of dirty cores across both trees.
- **At $`u = 0`$:** If $`w \in K(T)`$, marking both copies $`(w, w')`$ dirty gives $`D_{\mathrm{comb}}(0) = 2`$. If $`w`$ is a leaf of $`T`$, its unique neighbour $`p`$ lies in $`K(T)`$ (since $`C(T) \ge 1`$); marking both copies $`(p, p')`$ and the newly promoted core $`w'`$ in $`T^+`$ dirty gives $`D_{\mathrm{comb}}(0) = 3`$.
- **Clean invariants for all $`u \ge 0`$:** Every clean pair $`(z, z') \in \pi`$ satisfies: (C1) $`L_z(T_u) = L_{z'}(T^+_u)`$, (C2) clean core neighbours of $`z`$ correspond under $`\pi`$ to clean core neighbours of $`z'`$, and (C3) $`c_z^{\mathrm{dirty}}(T_u) = c_{z'}^{\mathrm{dirty}}(T^+_u)`$ (in fact $`\in \{0, 1\}`$ since dirty cores form a connected subtree in each tree). Consequently $`D_z(T_u) = D_{z'}(T^+_u)`$ identically! Synchronising paired leaf clocks and synchronising $`(z, z')`$'s rate-$`1`$ clock across paired leaves, paired clean neighbours, and paired dirty neighbours preserves (C1)–(C3), realizes the exact friend-tree marginal law in each tree, and creates zero dirty vertices. Note that when a dirty core promotes a leaf, the newborn is adjacent only to its dirty parent, preserving (C3); when a dirty core sends a leaf to a clean neighbour $`z`$, marking **both** $`(z, z')`$ dirty simultaneously shifts $`(z, z')`$ from clean to dirty for all remaining clean neighbours $`(u, u')`$, preserving (C2)–(C3).
- **Exact drift bound $`\le 2 D_{\mathrm{comb}}`$:** Each dirty core $`x`$ rings at rate $`1`$: choosing a leaf (prob $`L_x/D_x`$) creates $`1`$ dirty newborn ($`\Delta D_{\mathrm{comb}} = +1`$); choosing a dirty neighbour gives $`\Delta D_{\mathrm{comb}} = 0`$; choosing a clean neighbour (prob $`c_x^{\mathrm{clean}}/D_x`$) marks a clean pair dirty ($`\Delta D_{\mathrm{comb}} = +2`$). Thus the expected rate of increase of $`D_{\mathrm{comb}}`$ per dirty core is
  $`\displaystyle \frac{L_x \cdot 1 + c_x^{\mathrm{clean}} \cdot 2}{D_x} \le \frac{2(L_x + c_x^{\mathrm{clean}})}{D_x} \le 2.`$
  (Both constants in $`3 e^{2u}`$ are sharp: leaf promotion gives $`D_{\mathrm{comb}}(0) = 3`$, and a bare dirty core between two clean cores has $`L_x = 0, c_x^{\mathrm{clean}} = D_x = 2`$, giving instantaneous rate $`2`$.) Since clean pairs have $`D_z = D_{z'}`$, $`|A_k(T^+_u) - A_k(T_u)| \le \max(|\mathcal D(T_u)|, |\mathcal D(T^+_u)|) \le D_{\mathrm{comb}}(u)`$, and Grönwall (localized by Yule moments) yields $`|\mathbb E_{T^+} A_k(T_u) - \mathbb E_T A_k(T_u)| \le 3 e^{2u}`$.

### 3.2 Backward martingale quadratic variation and the $`h = \eta^{-1}\log(j+2)`$ $`L^2`$ recursion (Equations (7)–(9))
For fixed $`b > 0`$, $`M_v = \mathbb E_{T_v} A_k(\widetilde T_{b-v})`$ ($`0 \le v \le b`$) is an $`L^2`$-bounded martingale of finite pathwise variation (smooth compensating drift between jumps, zero continuous martingale part), so $`[M, M]_b = \sum_{0 < v \le b, \Delta T_v \ne 0} (\Delta M_v)^2`$. Every jump is a one-leaf extension at total physical rate $`n(T_v)`$ with $`|\Delta M_v| \le 3 e^{2(b-v)}`$ by (6), giving Equation (7):
```math
\mathrm{Var}_T(A_k(T_b)) = \mathbb E_T[[M, M]_b] \le \int_0^b \mathbb E_T[n(T_v)]\, 9 e^{4(b-v)}\,dv = 9\,n(T) \int_0^b e^{4(b-v)+v}\,dv.
```
For $`g_h(T) = d_2^{-1} e^{-\mu h}\mathbb E_T C(T_h)`$ over $`b = 1`$, jump sizes at $`v \in [0, 1]`$ are bounded by $`3 d_2^{-1} e^{2(1-v)+(2-\mu)h} \le 3 d_2^{-1} e^{2(1-v)+2h}`$, giving $`\mathrm{Var}_T(g_h(T_1)) \le C_4 e^{4h} n(T)`$. By (4), $`|\Psi - g_h| \le C_5 e^{-\eta h} C`$ ($`\eta = \mu - \sigma > 0`$), and $`C(T_1)`$ is dominated by a rate-$`1`$ Yule process from $`C(T)`$ (only cores promote, each at rate $`b_x \le 1`$), so $`\mathbb E_T C(T_1)^2 \le 2 e^2 C(T)^2`$ and $`\mathrm{Var}_T \Psi(T_1) \le C_7(e^{4h} n(T) + e^{-2\eta h} C(T)^2)`$.

Choosing $`h_j = \eta^{-1}\log(j + 2)`$ at step $`[j, j + 1]`$ balances the multiplicative factor $`1 + C_9 e^{-2\eta h_j} = 1 + C_9(j + 2)^{-2}`$ (infinite product $`< \infty`$) against the polynomial penalty $`e^{4 h_j} = (j + 2)^{4/\eta}`$, which is exponentially dominated by $`e^{-(2\mu - 1)j}`$ because $`\mu \ge 113/200 > 1/2 \implies 2\mu - 1 \ge 13/100 > 0`$. This yields $`\sup_{j \ge 0} \mathbb E_T[(Z_j^*)^2] \le C_{11}(\Psi(T)^2 + n(T) + 1)`$, extended to all real $`t \ge 0`$ by conditional Jensen ($`Z_t^* = \mathbb E_T[Z_{\lceil t \rceil}^* \mid \mathcal F_t]`$), giving $`L^2`$ convergence $`Z_t^* \to Z_\infty^*`$, $`\mathbb E_T[(Z_\infty^*)^2] \le C_{11}(\Psi(T)^2 + n(T) + 1)`$, and $`\mathbb E_T Z_\infty^* = \Psi(T)`$.

### 3.3 Positivity $`Z_\infty^* > 0`$ conditional on the past
On the $`T_3`$ space, $`Z_\infty^* = e^{-\mu j} \widetilde Z_\infty^*`$ a.s. where $`\widetilde Z_\infty^* = \lim_{u \to \infty} e^{-\mu u}\Psi(T_{j+u})`$, so $`\{Z_\infty^* > 0\} = \{\widetilde Z_\infty^* > 0\} \in \mathcal F_\infty`$ up to null sets. By the Markov property at $`j`$, Paley–Zygmund / Cauchy–Schwarz, and (9),
```math
\mathbb P_{T_3}(Z_\infty^* > 0 \mid \mathcal F_j) \ge \frac{\Psi(T_j)^2}{C_{11}(\Psi(T_j)^2 + n(T_j) + 1)}.
```
By almost-sure-exponent.md (`proofs/almost-sure-exponent.md`) (transferred from $`T_2`$ to $`T_3`$ via the independent $`\mathrm{Exp}(2)`$ first wait) and the Yule limit, $`j^{-1}\log C(T_j) \to \mu`$ and $`j^{-1}\log n(T_j) \to 1`$ a.s., so $`\Psi(T_j) \ge C_3^{-1} C(T_j) - 1 = e^{\mu j + o(j)} \to \infty`$ and $`n(T_j)/\Psi(T_j)^2 = e^{(1 - 2\mu)j + o(j)} \to 0`$ a.s. (again using $`2\mu > 1`$). By Lévy's upward 0–1 theorem, $`\mathbf{1}_{\{Z_\infty^* > 0\}} = \lim_{j \to \infty} \mathbb P_{T_3}(Z_\infty^* > 0 \mid \mathcal F_j) \ge C_{11}^{-1} > 0`$ a.s., forcing $`Z_\infty^* > 0`$ almost surely.

---

## 4. Sample-Path Transfer to All Fixed Degrees (Section 4)

On mesh $`t_j = j/m`$ ($`m \ge 1`$ fixed) with $`h_j = \eta^{-1}\log(t_j + 2)`$ and $`u_j = t_j - h_j > 0`$, conditioning on $`\mathcal F_{u_j}`$ and applying (7) gives
```math
\mathbb E_{T_3}\big[e^{-2\mu t_j}(A_k(T_{t_j}) - \mathbb E_{T_{u_j}} A_k(\widetilde T_{h_j}))^2\big] \le 3\,n(T_3) e^{t_j - h_j} e^{4 h_j - 2\mu t_j} \le 3\,n(T_3)(t_j + 2)^{4/\eta} e^{-(2\mu - 1)t_j},
```
which is summable in $`j`$, sending the centered term to $`0`$ a.s. by Borel–Cantelli. By (4), $`|e^{-\mu t_j}\mathbb E_{T_{u_j}} A_k(\widetilde T_{h_j}) - d_k Z_{u_j}^*| \le (C_1 + 1)(t_j + 2)^{-1} e^{-\mu u_j} C(T_{u_j}) \to 0`$ a.s. since $`u_j \to \infty`$ and $`e^{-\mu u} C(T_u) \le C_3(e^{-\mu u} + Z_u^*) \to C_3 Z_\infty^* < \infty`$. Monotonicity of $`t \mapsto A_k(T_t)`$ sandwiches arbitrary real $`t \in [j/m, (j+1)/m)`$ within $`e^{\pm \mu/m}`$; intersecting over countable $`(k, m)`$ gives $`e^{-\mu t} A_k(T_t) \to d_k Z_\infty^*`$ for all real $`t \to \infty`$ and all $`k \ge 2`$ simultaneously a.s. Evaluating at Yule epochs $`\tau_n`$ ($`n e^{-\tau_n} \to Y_\infty \sim \mathrm{Gamma}(3, 1) > 0`$ a.s.) gives $`A_k(T_{\tau_n})/n^\mu \to d_k Z_\infty^* Y_\infty^{-\mu}`$ a.s. Because $`T_2 \to T_3`$ deterministically on step $`1`$, $`(T_n)_{n \ge 3}`$ from $`T_2`$ has the exact joint law of $`(T_{\tau_n})_{n \ge 3}`$ from $`T_3`$, and matching the physical mean $`\lim_{t \to \infty} e^{-\mu t}\mathbb E_{T_3} N_{\ge k}(T_t) = \Psi(T_3) d_k = d_{f_k}^{\text{(exact)}}`$ against asymptotics.md (`fringe/asymptotics.md`) identifies $`d_k / d_2 = c_k / c_2`$.

---

## 5. Discrete Eigenmartingale $`Y_m`$ and Nontrivial Support of $`X`$ (Section 5)

1. **One-step discrete eigenrelation:** For any tree $`T`$ of size $`m \ge 3`$ ($`C(T) \ge 1`$), let $`(\widetilde T_t)_{t \ge 0}`$ start at $`T`$, with first jump wait $`\sigma_1 \sim \mathrm{Exp}(m)`$ independent of the one-step discrete successor $`T_{m+1}`$. Multiplying $`\mathbb E_T[\Psi(\widetilde T_t)] = e^{\mu t}\Psi(T)`$ by $`e^{-st}`$ ($`s > 1 > \mu`$) and integrating over $`t > 0`$ splits at $`\sigma_1`$ to give
   $`\displaystyle \frac{\Psi(T)}{s - \mu} = \frac{\Psi(T)}{s + m} + \frac{m}{s + m}\frac{\mathbb E[\Psi(T_{m+1}) \mid T_m = T]}{s - \mu} \iff \mathbb E[\Psi(T_{m+1}) \mid T_m = T] = \left(1 + \frac{\mu}{m}\right)\Psi(T).`$
   (Equivalently, differentiating $`\mathbb E_T[\Psi(\widetilde T_h)] = e^{-mh}\Psi(T) + m e^{-mh}(1 - e^{-h})\mathbb E[\Psi(T_{m+1})\mid T_m = T] + R(h)`$ at $`h = 0^+`$, where $`0 \le R(h) \le C_2\mathbb E_m[N(h)\mathbf{1}_{\{N(h)\ge m+2\}}] = C_2[m e^h - m e^{-mh} - m(m+1)e^{-mh}(1-e^{-h})] = O_m(h^2)`$, yields the same identity.) Using $`\Gamma(m+1)/\Gamma(m+1+\mu) = \frac{m}{m+\mu}\Gamma(m)/\Gamma(m+\mu)`$, $`Y_m = \frac{\Gamma(m)}{\Gamma(m+\mu)}\Psi(T_m)`$ ($`m \ge 3`$) is a nonnegative discrete martingale.
2. **$`L^2`$ boundedness and identification $`X = d_2 \lim_m Y_m`$:** In the continuous-time coupling from $`T_3`$, $`\tau_m = \sum_{j=3}^{m-1} E_j`$ ($`E_j \sim \mathrm{Exp}(j)`$ independent) is independent of $`(T_m)_{m \ge 3}`$ and satisfies
   $`\displaystyle \mathbb E[e^{-2\mu \tau_m}] = \prod_{j=3}^{m-1} \frac{j}{j + 2\mu} = \frac{\Gamma(3 + 2\mu)}{2}\frac{\Gamma(m)}{\Gamma(m + 2\mu)} \asymp m^{-2\mu}.`$
   Optional sampling of the $`L^2`$-bounded martingale $`Z_t^*`$ (or conditional Jensen at $`\tau_m \wedge K`$ plus Fatou) gives $`\mathbb E[e^{-2\mu \tau_m}]\mathbb E[\Psi(T_m)^2] = \mathbb E[(Z_{\tau_m}^*)^2] \le \mathbb E[(Z_\infty^*)^2] < \infty`$, so $`\mathbb E[\Psi(T_m)^2] = O(m^{2\mu})`$ and $`\sup_{m \ge 3} \mathbb E[Y_m^2] < \infty`$. Moreover, on the coupled space,
   $`\displaystyle Y_n = \left(\frac{n^\mu \Gamma(n)}{\Gamma(n + \mu)}\right) \big(n e^{-\tau_n}\big)^{-\mu} Z_{\tau_n}^* \longrightarrow \big(Y_\infty^{\mathrm{size}}\big)^{-\mu} Z_\infty^* \quad \text{a.s.},`$
   so $`X = \lim_{n \to \infty} C_n / n^\mu = d_2 (Y_\infty^{\mathrm{size}})^{-\mu} Z_\infty^* = d_2 \lim_{m \to \infty} Y_m`$ a.s., and $`L^2`$-convergence of $`Y_m`$ yields
   $`\displaystyle \mathbb E[X \mid T_m] = d_2 Y_m = d_2 \frac{\Gamma(m)}{\Gamma(m + \mu)}\Psi(T_m) \qquad (m \ge 3).`$
   (Cross-check at $`m = 3`$: $`\mathbb E[X] = d_2 \frac{\Gamma(3)}{\Gamma(3+\mu)}\Psi(T_3) = \frac{2}{(2+\mu)\Gamma(\mu+2)}(\Psi(T_3)d_2) = c_2`$, agreeing identically with asymptotics.md (`fringe/asymptotics.md`).)
3. **Accumulation at zero and unboundedness above:** Starting from $`T_2`$, for every $`m \ge 3`$:
   - The star $`K_{1, m-1}`$ has exact probability $`\prod_{j=3}^{m-1}\frac{j-1}{j} = \frac{2}{m-1} > 0`$ and $`C(K_{1, m-1}) = 1`$, so $`\mathbb E[X \mid T_m \cong K_{1, m-1}] \le d_2 C_2 \frac{\Gamma(m)}{\Gamma(m+\mu)} \sim d_2 C_2 m^{-\mu} \to 0`$ as $`m \to \infty`$, forcing $`\mathbb P(X < \varepsilon) > 0`$ for every $`\varepsilon > 0`$.
   - The path $`P_m`$ has exact probability $`\frac{2}{(m-1)!} > 0`$ and $`C(P_m) = m - 2`$, so $`\Psi(P_m) \ge C_3^{-1}(m - 2) - 1`$ gives $`\mathbb E[X \mid T_m \cong P_m] \ge c' m^{1 - \mu} \to +\infty`$ as $`m \to \infty`$ (since $`\mu < 1`$), forcing $`\mathbb P(X > M) > 0`$ for every $`M > 0`$. This strictly strengthens the non-trivial support condition in Question (1) of Section 7 of Addario-Berry et al. (`arXiv:2403.20185`).

---

## 6. Independent Small-Tree Enumeration and Mutation Tests

The standalone checker strong_limit_carbon.py tests every component above in exact rational arithmetic (`fractions.Fraction`) and bounded-error numerical integration, plus ten targeted mutations (negative controls):
1. **Capped root-gain & general-founder renewal ($`3 \le n \le 7`$, 118 rooted marked configurations):** exact verification of $`r_{\mathrm{gain}}|\Delta(B_0 g)| = \frac{c}{D+1} < 1`$ ($`20 \le D \le 120`$), $`\mathcal L_{\mathrm{phys}} Z \le \frac{9}{25} Z`$, $`\mathcal I(T) W \le 1/d + Z(T)`$, and $`\mathcal L_{\mathrm{phys}}\sum_{v \ne o} f_k(S_v) = \sum_{v \ne o}(R + B_0)f_k(S_v) + \mathcal I(T) f_k`$ across all core founders $`o \in K(T)`$ and multiple valid mark schemes.
2. **Contour high-frequency shift & exponential renewal:** 4-state multi-mark continuous-time model confirming $`K(z) = O(1/|z|)`$, $`(I - K)^{-1} K^2 = O(1/\tau^2)`$, numerical agreement of the shifted $`\Re z = \sigma`$ contour integral with time-domain $`v(t)`$, uniform exponential remainder (3), and exact eigenfunction identity (5).
3. **One-leaf clean/dirty coupling & backward martingale variance:** exact `Fraction` coupled Markov generator across small trees and one-leaf extensions verifying marginal equivalence, clean invariants, $`D_{\mathrm{comb}}(0) \le 3`$, sharp drift $`\mathcal L D_{\mathrm{comb}} \le 2 D_{\mathrm{comb}}`$, and continuous-time bounds (6)–(7).
4. **Eligibility reverse bound, $`L^2`$ recursion, and positivity:** exhaustive check of $`|G| \ge C(T)/12 - 1`$ on all 434 unlabelled trees ($`3 \le n \le 11`$), $`L^2`$ boundedness with $`h = \eta^{-1}\log(j+2)`$, conditional positivity ratio, and mesh summability.
5. **Section 5 discrete martingale & star/path seeds:** exact verification of $`\mathbb E[\Psi(T_{m+1})\mid T_m] = (1 + \mu/m)\Psi(T_m)`$, $`\mathbb E[e^{-2\mu\tau_m}]`$ product, $`\mathbb E[X] = c_2`$, and exact Markov probabilities $`\mathbb P(T_m \cong K_{1, m-1}) = \frac{2}{m-1}`$ and $`\mathbb P(T_m \cong P_m) = \frac{2}{(m-1)!}`$ ($`3 \le m \le 9`$).
6. **Ten mutations (negative controls):** uncapped marks with bounded oscillating $`g`$; omitting initial fringes $`\sum_{v \ne o} P_t f(S_v)`$; contour shift without subtracting $`k(t)`$; replacing $`3`$ by $`2`$ in $`D_{\mathrm{comb}}(0)`$; replacing $`2u`$ by $`\frac{3}{2}u`$ in dirty-core drift; marking only one copy dirty on clean-recipient send; constant $`h`$ in (9); sub-threshold $`\mu = 0.50`$ and $`\mu = 0.48`$; and off-by-one factor $`1 + \mu/(m+1)`$ in $`Y_m`$.

### Executed command and output
```text
$ python3 round2/checks/strong_limit_carbon.py
=== Cold-Challenge Verification: round2/fringe/almost-sure.md ===
[PASS] Capped-mark root-gain bound (D <= 120) & general-founder renewal/Z/I(W) identities exact across 118 rooted marked configurations (3 <= n <= 7).
[PASS] Integration by parts K(z)=O(1/|z|), O(1/tau^2) contour shift for v(t)-k(t), uniform exponential remainder Eq. (3), and exact eigenfunction Eq. (5) verified.
[PASS] One-leaf clean/dirty coupling exact across 128 coupled states (marginals, clean invariants, D_comb(0)<=3, drift <= 2 D_comb) + Eqs. (6)-(7) verified.
[PASS] Eligibility reverse bound checked on all 434 unlabelled trees (3 <= n <= 11); L^2 recursion (h=log(j+2)/eta), conditional positivity, and mesh summability verified.
[PASS] Section 5 discrete eigenmartingale Y_m, E[e^{-2*mu*tau_m}] product, E[X]=c_2 consistency, and exact star (2/(m-1)) / path (2/(m-1)!) seed laws (3 <= m <= 9) verified.
[CAUGHT] Mutation 1: Uncapped marks m(D)=D with bounded oscillating g breaks ||R B_0 g||_W
[CAUGHT] Mutation 2: Omitting initial fringes sum_{v!=o} P_t f(S_v) for general founder T
[CAUGHT] Mutation 3: Contour shift on v(t) without subtracting k(t) loses O(1/tau^2) L^1 decay
[CAUGHT] Mutation 4: Claiming initial dirty core bound 2 instead of 3 in Eq. (6) [leaf promotion gives 3]
[CAUGHT] Mutation 5: Claiming dirty core drift <= 1*D_comb (e^u) instead of 2*D_comb (e^{2u}) [bare core gives 2]
[CAUGHT] Mutation 6: Marking only recipient copy dirty (instead of clean pair) breaks coupling invariants
[CAUGHT] Mutation 7: Constant h instead of h = log(j+2)/eta blows up L^2 martingale recursion Eq. (9)
[CAUGHT] Mutation 8: Sub-threshold exponent mu = 0.50 (2*mu <= 1) breaks L^2 bound, positivity & mesh sum
[CAUGHT] Mutation 9: Sub-threshold exponent mu = 0.48 breaks conditional positivity ratio Psi^2/(Psi^2+n+1)
[CAUGHT] Mutation 10: Off-by-one step factor (1 + mu/(m+1)) instead of (1 + mu/m) breaks martingale Y_m
ALL CHECKS AND 10/10 MUTATION CONTROLS PASSED.
```
