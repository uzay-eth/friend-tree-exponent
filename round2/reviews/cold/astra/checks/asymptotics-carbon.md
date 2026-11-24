CORRECT

# Cold Review of `round2/fringe/asymptotics.md`

## Scope and Verdict

This is an independent cold review of `round2/fringe/asymptotics.md` establishing, for the undirected random friend tree started from a single edge $`T_2`$:
1. $`\mathbb E N_{\ge k}(n) \sim c_k n^\mu`$ as $`n \to \infty`$ for every fixed $`k \ge 2`$, with $`c_k > c_{k+1} > 0`$, and in particular $`\mathbb E C_n = \mathbb E N_{\ge 2}(n) \sim c_2 n^\mu`$;
2. $`\lim_{n \to \infty} \frac{\mathbb E N_{=k}(n)}{\mathbb E C_n} = p_k := \frac{c_k - c_{k+1}}{c_2} > 0`$ with $`\sum_{k \ge 2} p_k = 1`$;
3. $`p_k = \pi(D_o = k) \sim \mu A k^{-1-\mu}`$ as $`k \to \infty`$ with $`0 < A < \infty`$ (for the fixed-$`k`$ limiting distribution $`(p_k)_{k \ge 2}`$).

Sections 1–4 of `round2/fringe/renewal.md` are treated as a conditional dependency under separate review, while every prerequisite imported from `round2/fringe/renewal.md`, `round2/operator.md`, and `round2/degree-law.md` was checked directly. All steps—the discounted Laplace transform, Banach-space holomorphy and compactness, complex Krein–Rutman perturbation, residue signs, absence of non-real poles on $`\Re z = \mu`$, Wiener–Ikehara rescaling, $`T_3 \to T_2`$ first-birth convolution, Yule/Gamma fixed-size transfer, dominated convergence $`c_k \downarrow 0`$, and the local-operator degree-law bridge—are **mathematically complete and correct**.

---

## 1. Discounted Occupation, Holomorphy, and Renewal Transform (Section 1)

### 1.1 Planted root degree vs. actual degree, and $`f_2 \equiv 1`$
In `round2/fringe/renewal.md` Section 1, the actual tree from $`T_3`$ is rooted at the initial centre $`o`$ (which has degree $`\ge 2`$ at $`T_3`$ and forever after). Every ancestor of any noninitial core $`v \ne o`$ promoted a child to reach $`v`$, so every vertex on the path from $`o`$ to $`v`$ is core. Consequently, all leaf neighbours of $`v`$ in the actual tree are descendants of $`v`$, and $`v`$'s core neighbours consist of its core children in $`S_v`$ together with its unique parent toward $`o`$. Retaining the single inert exterior parent slot in $`S_v`$ makes
```math
D_{\mathrm{root}}(S_v) = L_v + (\#\{\text{core children of } v\} + 1) = L_v + c_v = D_v(T)
```
for every noninitial core $`v \ne o`$. Moreover, every reachable finite planted type $`S`$ has $`D_{\mathrm{root}}(S) \ge 2`$ (a newborn $`J_a`$ has $`L_o = 1`$ and $`c_o = 1`$, so $`D_{\mathrm{root}}(J_a) = 2`$; intrinsic transitions never decrease degrees; copied child subtrees $`Y, Y^+`$ retain a parent slot at a core vertex and hence satisfy $`D_{\mathrm{root}}(Y) \ge 2`$ and $`D_{\mathrm{root}}(Y^+) \ge 3`$). Thus $`f_2(S) = \mathbf{1}_{\{D_{\mathrm{root}}(S) \ge 2\}} \equiv 1`$ identically on finite planted types, and on the actual tree at time $`t`$ from $`T_3`$,
```math
N_{\ge k}(t) = \mathbf{1}_{\{D_o(t) \ge k\}} + \sum_{v \ne o} f_k(S_v(t)), \qquad M_t f_k \le \mathbb E_{T_3} N_{\ge k}(t) = M_t f_k + \Pr_{T_3}(D_o(t) \ge k) \le M_t f_k + 1.
```

### 1.2 Componentwise definitions and norm-holomorphy on $`\Re z > \rho`$
Let $`\mathcal S`$ be the countable set of isomorphism classes of finite planted types with marks in $`\mathcal A = \{1/2, 1/3, \dots\} \cup \{0\}`$, and let $`\mathcal B_W = \{F : \mathcal S \to \mathbb C : \|F\|_W = \sup_S |F(S)|/W(S) < \infty\}`$. From `renewal.md` Section 3 ($`d = 3/20`$, $`q = 9/25`$, $`\rho = q + d = 51/100`$), every $`S \in \mathcal S`$ has $`0 \le a_S \le 1/2`$ and $`h(c_o) \ge 7/10`$, so
```math
W(S) \ge \frac{7/10}{1/2 + 3/20} = \frac{14}{13} > 1,
```
giving $`\|1\|_W \le 13/14`$ and $`\|f_k\|_W \le 13/14`$. Extending $`R_t`$ and $`I_t`$ complex-linearly, $`|R_t F(S)| \le R_t |F|(S) \le \|F\|_W e^{\rho t} W(S)`$ and $`|I_t F| \le \|F\|_W (1/d + 2 e^{qt})`$.

Because $`\mathcal B_W`$ is an $`\ell^\infty`$-weighted space where strong continuity in norm need not hold uniformly over unbounded $`|S|`$, defining $`G_z f(S)`$, $`(U_z h)(S)`$, and $`\widehat I_z(F)`$ via componentwise scalar integrals (where $`t \mapsto R_t f(S)`$ and $`t \mapsto I_t F`$ are continuous by the first-event renewal equation and Yule domination) avoids any Bochner-measurability issue. About any $`z_0`$ with $`\Re z_0 > \rho`$, setting
```math
c_m(S) = \int_0^\infty \frac{(-t)^m}{m!} e^{-z_0 t} R_t f(S)\,dt
```
gives $`\|c_m\|_W \le \|f\|_W (\Re z_0 - \rho)^{-(m+1)}`$, so $`\sum_{m \ge 0} (z - z_0)^m c_m`$ converges in $`\|\cdot\|_W`$-norm on $`|z - z_0| < \Re z_0 - \rho`$ and equals $`G_z f`$ pointwise by dominated convergence. Identical majorants yield norm-convergent Taylor series for $`U_z \in \mathcal L(C(\mathcal A), \mathcal B_W)`$ on $`\Re z > \rho = 51/100`$ and for $`\widehat I_z \in \mathcal B_W^*`$ on $`\Re z > q = 9/25`$.

### 1.3 Continuity of $`j G_z f`$ and compactness/holomorphy of $`K_z`$
Starting from $`J_a`$, the original particle's killing rate $`a`$ is constant and acts independently of its intrinsic evolution $`S^{\mathrm{int}}(t)`$ (whose generator $`L^0`$ and child promotion marks $`a_Y = 1/D_{\mathrm{root}}(\text{at promotion})`$ do not read the root mark $`a`$). For $`f = 1`$ or $`f = f_k`$ (or $`f_k - f_{k+1}`$), $`f(S^{\mathrm{int}}(t))`$ depends only on $`D_{\mathrm{root}}(S^{\mathrm{int}}(t))`$, which is independent of $`a`$. Decomposing the residual family of $`J_a`$ into the surviving original particle and the child subtrees copied while it lives, and applying Fubini on $`\Re z > \rho`$ (using $`\sum_{Y \text{ at } S} a_Y W(Y) = (a_S + d)W(S) - h(c_o)`$ and $`\mathbb E_{\mathrm{int}}[(a+d)W(S^{\mathrm{int}}(t))] \le e^{qt}`$ with $`q = 9/25 < 51/100 < \Re z`$), yields
```math
(G_z f)(J_a) = \int_0^\infty e^{-(z+a)t} \Phi_{z, f}(t)\,dt, \qquad
(K_z h)(a) = \int_0^\infty e^{-(z+a)t} \Psi_{z, h}(t)\,dt,
```
where $`\Phi_{z, f}(t)`$ and $`\Psi_{z, h}(t)`$ are independent of $`a \in \mathcal A`$ and satisfy $`|\Phi_{z, f}(t)| \le \|f\|_\infty\big(1 + \frac{e^{qt}}{\Re z - \rho}\big)`$ and $`|\Psi_{z, h}(t)| \le \|h\|_\infty\big(1 + \frac{e^{qt}}{\Re z - \rho}\big)`$. Since $`|e^{-a_1 t} - e^{-a_2 t}| \le |a_1 - a_2| t`$ for $`a_1, a_2 \in \mathcal A \subset [0, 1/2]`$, we obtain $`j G_z f \in C(\mathcal A)`$, and $`\{K_z h : \|h\|_\infty \le 1\}`$ is uniformly bounded and equi-Lipschitz on the compact metric space $`\mathcal A`$. By Arzelà–Ascoli, $`K_z`$ is compact on $`C(\mathcal A)`$, and $`K_z = j U_z`$ is holomorphic in operator norm on $`\{\Re z > \rho\}`$.

### 1.4 Generation expansion and Equation (3) on $`\Re z > \mu`$
For real $`\sigma = \Re z > \mu \ge 113/200 > \rho`$, `renewal.md` Section 4 gives $`r(K_\sigma) < r(K_\mu) = 1`$. By positivity of $`R_t`$ and $`b_S \ge 0`$, $`|(K_z h)(a)| \le (K_\sigma |h|)(a)`$ for every $`h \in C(\mathcal A)`$, and induction gives $`|(K_z^m h)(a)| \le (K_\sigma^m |h|)(a)`$, hence $`\|K_z^m\|_{\mathcal L(C(\mathcal A))} \le \|K_\sigma^m\|_{\mathcal L(C(\mathcal A))}`$. By Gelfand's formula $`\sum_{m \ge 0} \|K_\sigma^m\| < \infty`$, so $`(I - K_z)^{-1} = \sum_{m \ge 0} K_z^m`$ converges absolutely in $`\mathcal L(C(\mathcal A))`$. Grouping $`P_t f = \sum_{n \ge 0} P_t^{(n)} f`$ by the generation $`n`$ of direct $`B`$-births on the nonexplosive branching process, Tonelli at $`\sigma = \Re z`$ bounds the sum of absolute integrals by $`\widehat I_\sigma\big(G_\sigma |f| + U_\sigma(I - K_\sigma)^{-1} j G_\sigma |f|\big) \le C_\sigma \widehat I_\sigma(W) < \infty`$ (since $`\sigma > \mu > 9/25`$). Fubini and norm continuity of $`U_z`$ and $`\widehat I_z`$ then give Equation (3):
```math
\int_0^\infty e^{-zt} M_t f\,dt = \widehat I_z\big(G_z f + U_z(I - K_z)^{-1} j G_z f\big) \qquad (\Re z > \mu).
```

---

## 2. Spectral Perturbation, Residue Signs, and Peripheral Line (Section 2)

### 2.1 Complexification, Krein–Rutman simplicity, and full support of $`\ell`$
On $`C_{\mathbb R}(\mathcal A)`$, the cone $`P = \{v \ge 0\}`$ has nonempty interior $`\{v : \min_{\mathcal A} v > 0\}`$. Because $`\{1/j : j \ge 2\}`$ is dense in $`\mathcal A`$, any $`v \in P \setminus \{0\}`$ has $`v(1/j) > 0`$ for some $`j \ge 2`$. Starting from $`J_a`$ ($`a \in [0, 1/2]`$), keeping the root silent while its leaf rings $`j - 2`$ times in $`[0, 1]`$, keeping all $`j`$ vertices silent in $`[1, 2]`$, and surviving killing through $`t = 2`$ (probability $`\ge e^{-2a} \ge e^{-1}`$, uniform down to $`a = 0`$) gives $`p_{S(t)} = 1/j`$ and $`b_{S(t)} = (j-1)/j \ge 1/2`$ throughout $`[1, 2]`$. Thus $`(K_\mu v)(a) > 0`$ for all $`a \in \mathcal A`$, i.e., $`K_\mu(P \setminus \{0\}) \subset \mathrm{int}(P)`$.

By the strong Krein–Rutman theorem, $`1 = r(K_\mu)`$ is an algebraically simple eigenvalue on $`C_{\mathbb R}(\mathcal A)`$ with $`g \in \mathrm{int}(P)`$ and positive dual eigenmeasure $`\ell`$ normalized by $`\ell(g) = 1`$. Extending to the complexification $`C_{\mathbb C}(\mathcal A)`$ preserves compactness, spectral radius, geometric simplicity (real and imaginary parts of $`w \in \ker(K_\mu - I)`$ lie in $`\mathbb R g`$), and algebraic simplicity (if $`(K_\mu - I)w = g`$, taking real parts and applying $`\ell`$ yields $`0 = \ell(g) = 1`$). For any nonzero $`v \in P`$, $`K_\mu v \ge c g`$ with $`c = (\min K_\mu v)/(\max g) > 0`$, so $`\ell(v) = \ell(K_\mu v) \ge c > 0`$; in particular $`\mathrm{supp}(\ell) = \mathcal A`$, and any continuous $`w \ge 0`$ with $`\ell(w) = 0`$ vanishes identically on $`\mathcal A`$.

### 2.2 Re-derivation of $`\lambda'(\mu)`$ and the residue sign in Equation (4)
Isolating $`1 \in \sigma(K_\mu)`$ by a small circle $`\Gamma`$, the Riesz projection $`\Pi_z = \frac{1}{2\pi i}\oint_\Gamma (\zeta I - K_z)^{-1}\,d\zeta`$ is norm-holomorphic near $`z = \mu`$ with $`\mathrm{rank}(\Pi_z) = \mathrm{rank}(\Pi_\mu) = 1`$ and $`\Pi_\mu = g \otimes \ell`$. Normalizing $`\tilde g(z) = \Pi_z g / \ell(\Pi_z g)`$ so $`\ell(\tilde g(z)) \equiv 1`$, we have $`K_z \tilde g(z) = \lambda(z)\tilde g(z)`$ with $`\lambda(z) = \ell(K_z \tilde g(z))`$ holomorphic near $`\mu`$ and $`\lambda(\mu) = 1`$. Differentiating at $`z = \mu`$ gives
```math
K'_\mu g + K_\mu \tilde g'(\mu) = \lambda'(\mu) g + \tilde g'(\mu).
```
Applying $`\ell`$ cancels $`\ell(K_\mu \tilde g'(\mu)) = \ell(\tilde g'(\mu))`$ and leaves $`\lambda'(\mu) = \ell(K'_\mu g)`$. Since
```math
-(K'_\mu g)(a) = \int_0^\infty t e^{-\mu t} R_t[b_{(\cdot)} g(p_{(\cdot)})](J_a)\,dt > 0 \qquad \forall a \in \mathcal A
```
(strictly positive on every interval by survival of $`J_a`$ with $`\min g > 0`$), we have $`-\ell(K'_\mu g) > 0`$, i.e., $`\lambda'(\mu) < 0`$.

On $`\mathrm{ran}(I - \Pi_z)`$, $`(I - K_z)^{-1}(I - \Pi_z)`$ is holomorphic near $`\mu`$. On $`\mathrm{ran}(\Pi_z)`$,
```math
(I - K_z)^{-1} \Pi_z = \frac{1}{1 - \lambda(z)} \Pi_z = \frac{1}{-\lambda'(\mu)(z - \mu) + O((z-\mu)^2)} \Pi_z = \frac{g \otimes \ell}{[-\ell(K'_\mu g)](z - \mu)} + O(1).
```
Sign check: for real $`z = \sigma > \mu`$, $`z - \mu > 0`$ and $`\lambda(\sigma) < 1 \implies 1 - \lambda(\sigma) > 0`$, matching $`[-\ell(K'_\mu g)](\sigma - \mu) > 0`$. Since $`\widehat I_z`$, $`G_z f`$, $`U_z`$, and $`j G_z f`$ are all norm-holomorphic on $`\Re z > \rho`$ (containing $`\Re z = \mu`$), $`\int_0^\infty e^{-zt} M_t f\,dt - \frac{d_f}{z - \mu}`$ is holomorphic in a neighborhood of $`z = \mu`$ with
```math
d_f = \frac{\widehat I_\mu(U_\mu g)\,\ell(j G_\mu f)}{-\ell(K'_\mu g)}.
```

### 2.3 Strict positivity of $`d_f`$ for $`f_k`$ and $`f_k - f_{k+1}`$
- $`\widehat I_\mu(U_\mu g) \ge \int_0^1 e^{-\mu t} \mathbb E_{T_3}[b_o (U_\mu g)(J_{p_o})]\,dt = \int_0^1 e^{-\mu t} \mathbb E_{T_3}[b_o g(p_o)]\,dt > 0`$, because $`(U_\mu g)(J_a) = (K_\mu g)(a) = g(a) \ge \min g > 0`$ and keeping the centre $`o`$ silent on $`[0, 1]`$ (probability $`e^{-1}`$) maintains $`c_o = 0, L_o \ge 2 \implies b_o = 1`$.
- For any $`k \ge 2`$ and $`a \in \mathcal A`$, starting from $`J_a`$, keeping the root silent on $`[0, 1]`$ while its initial leaf rings **exact** $`k - 2`$ times and newborn leaves stay silent, then keeping all $`k`$ vertices silent on $`[1, 2]`$ while the particle survives killing through $`t = 2`$, yields $`D_{\mathrm{root}}(S(t)) = k`$ throughout $`[1, 2]`$ with positive probability. Hence $`j G_\mu(f_k - f_{k+1})(a) > 0`$ for all $`a \in \mathcal A`$, so $`\ell(j G_\mu(f_k - f_{k+1})) > 0`$. Thus $`d_{f_k} - d_{f_{k+1}} = d_{f_k - f_{k+1}} > 0`$ and $`d_{f_k} > 0`$.

### 2.4 Invertibility of $`I - K_{\mu + i\tau}`$ for $`\tau \ne 0`$
Fix $`\tau \ne 0`$. Because $`K_{\mu + i\tau}`$ is compact, non-invertibility of $`I - K_{\mu + i\tau}`$ requires a nonzero $`v \in C(\mathcal A)`$ with $`K_{\mu + i\tau} v = v`$. Pointwise $`|v| \le K_\mu |v|`$, so $`w := K_\mu |v| - |v| \in C_{\mathbb R}(\mathcal A)`$ satisfies $`w \ge 0`$ and $`\ell(w) = \ell(K_\mu |v|) - \ell(|v|) = 0`$. By strict positivity of $`\ell`$ on $`P \setminus \{0\}`$, $`w \equiv 0`$, so $`|v| = K_\mu |v| = c g`$ with $`c > 0`$; in particular $`v(1/2) \ne 0`$.

Write $`(K_z v)(a) = \sum_{j \ge 2} \int_0^\infty e^{-zt} \pi_{a, j}(t) v(1/j)\,dt`$, where $`\pi_{a, j}(t) = R_t[b_{(\cdot)}\mathbf{1}_{\{p_{(\cdot)} = 1/j\}}](J_a) \ge 0`$. On the event that $`J_a`$ experiences no root ring, no leaf ring, and no killing on $`[0, t]`$ (probability $`e^{-(2+a)t}`$), $`J_a`$ remains in its initial state with $`p = 1/2`$ and $`b = 1/2`$, so $`\pi_{a, 2}(t) \ge \frac{1}{2} e^{-(2+a)t} > 0`$ for all $`t > 0`$. Then
```math
|v(a)| = |(K_{\mu + i\tau} v)(a)| \le \left|\int_0^\infty e^{-\mu t} \pi_{a, 2}(t) e^{-i\tau t} v(1/2)\,dt\right| + \sum_{j \ge 3} \int_0^\infty e^{-\mu t} \pi_{a, j}(t)|v(1/j)|\,dt < (K_\mu |v|)(a) = |v(a)|,
```
where the strict inequality holds because $`\tau \ne 0`$, $`|v(1/2)| > 0`$, and $`e^{-\mu t}\pi_{a, 2}(t) > 0`$ on $`(0, \infty)`$ forces $`|\int_0^\infty e^{-\mu t}\pi_{a, 2}(t)e^{-i\tau t}\,dt| < \int_0^\infty e^{-\mu t}\pi_{a, 2}(t)\,dt`$. Contradiction. Hence $`I - K_{\mu + i\tau}`$ is invertible on $`C(\mathcal A)`$, and norm-holomorphy of $`z \mapsto K_z`$ extends $`(I - K_z)^{-1}`$ analytically to a neighborhood of $`\mu + i\tau`$ via Neumann series. Thus (3) minus $`d_f/(z - \mu)`$ extends analytically across every point of $`\Re z = \mu`$.

---

## 3. Tauberian Step, $`T_3 \to T_2`$ Return, and Fixed-Size Transfer (Section 3)

### 3.1 Wiener–Ikehara rescaling and Debruyne–Vindas citation
Pathwise from $`T_3`$, noninitial core vertices are never lost and their degrees never decrease, so $`t \mapsto M_t f_k`$ is nonnegative, locally bounded ($`\le 3 e^t`$), and nondecreasing on $`[0, \infty)`$. Put $`S(x) = M_{x/\mu} f_k \in L^1_{\mathrm{loc}}[0, \infty)`$, which is nondecreasing and hence log-linearly slowly decreasing in the sense of Debruyne and Vindas (`arXiv:1611.09765`, Eq. (1.1)). With $`t = x/\mu`$ ($`dx = \mu\,dt`$), for $`\Re w > 1`$:
```math
\mathcal L\{S; w\} = \int_0^\infty e^{-w x} M_{x/\mu} f_k\,dx = \mu \int_0^\infty e^{-(\mu w)t} M_t f_k\,dt.
```
Near $`w = 1`$ ($`z = \mu w`$ near $`\mu`$), the Jacobian $`\mu`$ and denominator $`\mu w - \mu = \mu(w - 1)`$ cancel:
```math
\mathcal L\{S; w\} - \frac{d_{f_k}}{w - 1} = \mu\left(\int_0^\infty e^{-(\mu w)t} M_t f_k\,dt - \frac{d_{f_k}}{\mu w - \mu}\right),
```
which extends analytically to a neighborhood of every point of $`\Re w = 1`$. Theorem 1.1 and the paragraph immediately following it in Debruyne and Vindas (`arXiv:1611.09765`, lines 31–79, via the Riemann–Lebesgue lemma on analytic/$`L^1_{\mathrm{loc}}`$ boundary extensions; see also Theorem 3.6) give $`S(x) \sim d_{f_k} e^x`$ pointwise as $`x \to \infty`$, hence $`M_t f_k = S(\mu t) \sim d_{f_k} e^{\mu t}`$ and $`\mathbb E_{T_3} N_{\ge k}(t) \sim d_{f_k} e^{\mu t}`$.

### 3.2 Return to the initial edge $`T_2`$ (Equation (5))
At $`T_2`$ both vertices have degree $`1 < 2 \le k`$, so $`N_{\ge k}(T_2) = 0`$. The first birth time $`\tau_1 \sim \mathrm{Exp}(2)`$ deterministically creates $`T_3`$ and is independent of the subsequent process, giving the exact convolution
```math
e^{-\mu t} \mathbb E_{T_2} N_{\ge k}(t) = \int_0^t 2 e^{-(2+\mu)u} \big(e^{-\mu(t-u)}\mathbb E_{T_3} N_{\ge k}(t - u)\big)\,du.
```
Since $`s \mapsto e^{-\mu s}\mathbb E_{T_3} N_{\ge k}(s)`$ is continuous on $`[0, \infty)`$ and converges to $`d_{f_k}`$, it is globally bounded on $`[0, \infty)`$. Dominated convergence against $`2 e^{-(2+\mu)u} \in L^1(0, \infty)`$ yields
```math
\lim_{t \to \infty} e^{-\mu t} \mathbb E_{T_2} N_{\ge k}(t) = d_{f_k} \int_0^\infty 2 e^{-(2+\mu)u}\,du = \frac{2}{2 + \mu} d_{f_k} = d_k^*.
```

### 3.3 Fixed-size transfer and $`\sum_{k \ge 2} p_k = 1`$ (Equations (64)–(75))
Starting from $`T_2`$, the Yule size $`N(t)`$ is independent of the jump chain $`(T_n)_{n \ge 2}`$ with $`\Pr_{T_2}(N(t) = n) = (n-1) q^2 (1 - q)^{n-2}`$ ($`q = e^{-t}`$, $`n \ge 2`$), so $`\mathbb E_{T_2} N_{\ge k}(t) = q^2 \sum_{n \ge 2} (n - 1)(1 - q)^{n-2} a_{n, k}`$. For $`\nu_x = x^{-(\mu + 2)} \sum_{n \ge 2} (n - 1) a_{n, k}\,\delta_{n/x}`$ and fixed $`s > 0`$, setting $`1 - q = e^{-s/x}`$ (so $`q = 1 - e^{-s/x} \sim s/x`$, $`t = -\log q \to \infty`$, and $`(1-q)^{n-2} = e^{2s/x} e^{-sn/x}`$) gives
```math
\int_{[0, \infty)} e^{-sy}\,d\nu_x(y) = e^{-2s/x} (x q)^{-(\mu + 2)} \big(e^{-\mu t} \mathbb E_{T_2} N_{\ge k}(t)\big) \xrightarrow{x \to \infty} d_k^* s^{-(\mu + 2)}.
```
Tilting $`\tilde\nu_x(dy) = e^{-y}\nu_x(dy)`$ has mass $`Z_x = \int e^{-y}\,d\nu_x(y) \to d_k^* > 0`$, and the probability measures $`P_x = Z_x^{-1}\tilde\nu_x`$ satisfy $`\int e^{-uy}\,dP_x(y) \to (1 + u)^{-(\mu + 2)}`$ for all $`u \ge 0`$, the Laplace transform of $`P_\infty(dy) = \frac{y^{\mu+1} e^{-y}}{\Gamma(\mu + 2)}\,dy`$. By Feller's continuity theorem, $`P_x \Rightarrow P_\infty`$. For any bounded interval $`[a, b]`$ with $`0 < a < b < \infty`$, $`y \mapsto e^y \mathbf{1}_{[a, b]}(y)`$ is bounded and continuous outside $`\{a, b\}`$, which has $`P_\infty`$-measure zero, so Portmanteau yields
```math
\nu_x([a, b]) = Z_x \int e^y \mathbf{1}_{[a, b]}(y)\,dP_x(y) \longrightarrow \frac{d_k^*}{\Gamma(\mu + 2)} \int_a^b y^{\mu + 1}\,dy.
```
Since $`N_{\ge k}(T_{m+1}) \ge N_{\ge k}(T_m)`$ pathwise, $`a_{m, k}`$ is nondecreasing in $`m`$. For $`0 < h < 1`$ and integer $`x = n`$:
```math
\frac{\nu_n([1 - h, 1])}{n^{-2}\sum_{(1-h)n \le m \le n} (m - 1)} \le \frac{a_{n, k}}{n^\mu} \le \frac{\nu_n([1, 1 + h])}{n^{-2}\sum_{n \le m \le (1+h)n} (m - 1)}.
```
Sending $`n \to \infty`$ gives $`\frac{d_k^*}{\Gamma(\mu+2)}\frac{\int_{1-h}^1 y^{\mu+1}\,dy}{\int_{1-h}^1 y\,dy} \le \liminf \frac{a_{n, k}}{n^\mu} \le \limsup \frac{a_{n, k}}{n^\mu} \le \frac{d_k^*}{\Gamma(\mu+2)}\frac{\int_1^{1+h} y^{\mu+1}\,dy}{\int_1^{1+h} y\,dy}`$, and letting $`h \downarrow 0`$ proves $`a_{n, k} \sim c_k n^\mu`$ with $`c_k = d_k^* / \Gamma(\mu + 2)`$.

Finally, at each finite $`t`$ the nonexplosive residual family from $`J_a`$ has finitely many particles, each of finite degree, dominated in score by $`\sum W \in L^1`$ ($`R_t W(J_a) \le e^{\rho t} W(J_a)`$), so $`R_t f_k(J_a) \downarrow 0`$ as $`k \to \infty`$. Dominated convergence against $`e^{-(\mu-\rho)t} W(J_a) \in L^1(0, \infty)`$ gives $`j G_\mu f_k(a) \downarrow 0`$ pointwise on $`\mathcal A`$; uniform boundedness $`0 \le j G_\mu f_k(a) \le \frac{W(J_a)}{\mu - \rho} \le \frac{20/3}{\mu - \rho}`$ and finiteness of $`\ell(\mathcal A)`$ give $`\ell(j G_\mu f_k) \downarrow 0`$, hence $`c_k \downarrow 0`$ and $`\sum_{k \ge 2} p_k = 1`$.

---

## 4. Tail of the Limiting Law and Prerequisite Checks (Section 4)

### 4.1 Bridge to `round2/operator.md` and `round2/degree-law.md`
In `round2/operator.md` Section 2, $`m_t(f) = \mathbb E_{T_3}\sum_{u \in K_t} f(G_t, u)`$, $`F_3(t) = m_t(1)`$, $`w(t) = \exp(-\mu t - \delta|t - T|)`$, $`Z = \int_0^\infty w(t) F_3(t)\,dt`$, and $`\pi`$ is a weak limit of $`\pi_{\delta, T} = Z^{-1}\int_0^\infty w(t) m_t\,dt`$ along $`\delta \downarrow 0`$ and $`T \to \infty`$ chosen so $`w(0)/Z \to 0`$.
- For $`\tilde f_k(G, o) = \mathbf{1}_{\{D_o \ge k\}}`$ on $`\mathcal X`$, $`m_t(\tilde f_k) = \mathbb E_{T_3} N_{\ge k}(t) \sim d_{f_k} e^{\mu t}`$ and $`m_t(1) = \mathbb E_{T_3} N_{\ge 2}(t) \sim d_{f_2} e^{\mu t} = d_1 e^{\mu t}`$ (since $`f_2 \equiv 1`$), so $`m_t(\tilde f_k)/m_t(1) \to d_{f_k}/d_1 = c_k / c_2`$.

• For $`0 < \delta < \mu`$, $`T > 0`$, and every $`t \ge 0`$, the reverse triangle inequality $`|t - T| \ge T - t`$ gives

```math
\frac{w(t)}{w(0)} = \exp\big(-\mu t + \delta(T - |t - T|)\big) \le e^{-(\mu - \delta)t} \le 1.
```

Hence for every fixed $`M > 0`$, $`Z^{-1}\int_0^M w(t) m_t(1)\,dt \le \frac{w(0)}{Z}\int_0^M 3 e^t\,dt \to 0`$, which forces $`\pi_{\delta, T}(\tilde f_k) \to c_k / c_2`$.
- Because agreement on the marked ball of radius $`1`$ fixes both $`L_o`$ and $`c_o`$, $`\tilde f_k`$ is bounded and clopen (continuous) on $`\mathcal X`$. Weak convergence $`\pi_{\delta_j, T_j} \Rightarrow \pi`$ therefore gives $`\pi(D_o \ge k) = c_k / c_2`$ and $`p_k = \frac{c_k - c_{k+1}}{c_2} = \pi(D_o = k)`$ (independent of which subsequential limit $`\pi`$ is taken).

### 4.2 Audit of prerequisites used
- **`renewal.md` Sections 1–4:** Rational inequality (4) $`(1 - c/D)[h(c+1) - h(c) + (1 + dD)^{-1}] \le q h(c)`$ with $`d = 3/20`$, $`q = 9/25`$ holds for all $`D \ge \max(2, c)`$ (the tightest finite case is $`(c, D) = (1, 3)`$, giving $`313/870 \le 9/25 = 313.2/870`$; for $`D \ge 100`$ the left side is at most $`1/(1 + 15) = 1/16 < 9/25 \cdot 7/10 = 63/250`$). Child-weight cancellation gives $`R W \le \rho W`$ with $`\rho = 51/100 < 113/200 \le \mu`$ and $`I_t W \le 20/3 + 2 e^{(9/25)t}`$.
- **`operator.md` Sections 1–2 and `degree-law.md`:** $`\pi`$ satisfies $`\pi Q_t = e^{\mu t}\pi`$ and $`\pi(e^{\theta c_o}) < \infty`$ for $`0 < \theta < \log(1 + \mu)`$. Decomposing $`Q_t f = S_t f + \int_0^t Q_s B S_{t-s} f\,ds`$ by unique promotion birth times within a cell gives $`\pi f = \int_0^\infty \mu e^{-\mu a} (\beta S_a f)\,da`$ with $`\beta = \pi B / \mu`$. Under $`\beta`$, $`e^{-t} D(t)`$ is a nonnegative supermartingale dominated by a rate-one Yule martingale $`e^{-t} Y(t)`$ from $`Y(0) = 2`$ ($`\mathbb E[\sup_t (e^{-t} Y(t))^2] \le 24`$), giving $`A = \mathbb E_\beta[W^\mu] < \infty`$; silencing the newborn's clock on $`[0, 1]`$ while its leaf rings at least twice (probability $`e^{-1}(1 - 2e^{-1}) > 0`$) yields $`L_o(1) \ge 3`$ and submartingale leaf drift $`r_o - b_o \ge L_o - 1`$, proving $`A > 0`$ and $`\pi(D_o \ge k) = \mathbb E_\beta[e^{-\mu T_k}; T_k < \infty] \sim A k^{-\mu}`$. Finally, compensating upward crossings of $`\{D_v \ge k + 1\}`$ in a cell (where newborns start at degree $`2 < k + 1`$ and $`g_k(G, o) = (k - c_o + \rho_o)\mathbf{1}_{\{D_o = k\}}`$ is a radius-$`2`$ cylinder bounded by $`k`$) yields $`\mu \pi(D_o \ge k + 1) = \pi[(k - c_o + \rho_o)\mathbf{1}_{\{D_o = k\}}]`$, which sandwiches $`p_k = \pi(D_o = k) \sim \mu A k^{-1-\mu}`$ using $`0 \le \rho_o \le c_o/2`$ and exponential core-degree tails $`\pi(c_o > \delta k) = o(k^{-1-\mu})`$. Note that `asymptotics.md` explicitly states at the end of Section 4 that the large-$`k`$ asymptotic applies to the fixed-$`k`$ limiting law $`(p_k)_{k \ge 2}`$.

---

## 5. Reproducible Enumeration, Transfer Check, and Deliberate Errors

The standalone script `round2/checks/asymptotics_carbon.py` recomputes in exact rational arithmetic (`fractions.Fraction`) and bounded-error analytic checks:
1. **Exact friend-tree enumeration ($`2 \le n \le 10`$, 200 unlabelled state layers):** computes exact $`a_{n, k} = \mathbb E_{T_2}[N_{\ge k}(n)]`$ and $`\mathbb E_{T_3}[N_{\ge k}(n)]`$, confirming $`a_{n+1, k} \ge a_{n, k}`$, strict increase for $`n + 1 \ge k + 1`$, and exact degree telescoping $`\sum_{k=2}^{n-1}(a_{n, k} - a_{n, k+1}) = a_{n,2}`$.
2. **Marked descendant trees & generator identity ($`3 \le n \le 7`$, 88 reachable marked states from $`T_3`$):** verifies $`D_{\mathrm{root}}(S_v) = D_v(T)`$, $`p_{\mathrm{parent}} \le a_v`$, and $`\mathcal L_{\mathrm{phys}}\big(\sum_{v \ne o} f_k(S_v)\big) = \sum_{v \ne o}(R + B)f_k(S_v) + I(f_k)`$ in exact rationals.
3. **Transfer recomputation:** verifies exact $`q`$-series $`\mathrm{Exp}(2)`$ convolution between $`T_3`$ and $`T_2`$, exact rational Beta-integral Laplace ratio $`\frac{2}{2 + z}`$ term-by-term on small-$`n`$ enumeration, exact closed-form summation $`q^2 \sum_{n \ge 2}(n-1)(1-q)^{n-2}\frac{\Gamma(n+\mu)}{\Gamma(n)} = \Gamma(\mu + 2)e^{\mu t}`$ confirming $`c_k = d_k^*/\Gamma(\mu + 2)`$, $`\nu_x`$ Laplace convergence to $`d_k^* s^{-(\mu+2)}`$, and the $`[1-h, 1] / [1, 1+h]`$ sandwich bounds.
4. **Multi-mark renewal & Tauberian pole checks:** verifies Equation (3), $`\lambda'(\mu) = \ell(K'_\mu g) < 0`$, contour residue match to Equation (4), peripheral invertibility on $`\Re z = \mu`$ ($`\tau \ne 0`$), and Wiener–Ikehara rescaling residue $`d_f`$ at $`w = 1`$.
5. **Eight deliberate planted errors (negative controls):** wrong Yule weight $`n`$ vs. $`n-1`$; wrong $`T_3 \to T_2`$ wait factor $`1/(1+z)`$ vs. $`2/(2+z)`$; wrong $`\Gamma(\mu+1)`$ vs. $`\Gamma(\mu+2)`$; omitted $`dx = \mu\,dt`$ Jacobian in $`M_{x/\mu} f_k`$; flipped residue sign $`+\ell(K'_\mu g)`$; omitted parent slot in $`D_{\mathrm{root}}(S_v)`$; lattice reproduction kernel creating poles at $`\mu + 2\pi i`$; and perturbing $`q`$ from $`9/25`$ to $`359/1000`$ in `renewal.md` Eq. (4).

### Executed command and output
```text
$ python3 round2/checks/asymptotics_carbon.py
=== Cold Review Verification: round2/fringe/asymptotics.md ===
[PASS] Exact rational friend-tree enumeration (2 <= n <= 10, 200 unlabelled state layers): monotonicity & degree telescoping exact.
       Sample exact E[C_n] (n=3..10): n=3: 1, n=4: 4/3, n=5: 19/12, n=6: 647/360, n=7: 17189/8640, n=8: 4717423/2177280, n=9: 101568829/43545600, n=10: 200684675729/80621568000
[PASS] Marked descendant trees & renewal generator M'_t f_k = M_t(R+B)f_k + I_t f_k exact across 88 reachable marked configurations (3 <= n <= 7).
[PASS] Exact small-n q-series & Beta-Laplace transfer (ratio 2/(2+z)) and fixed-size Gamma(mu+2) / sandwich transfer verified.
[PASS] Multi-mark renewal resolvent Eq. (3), residue Eq. (4) sign -ell(K'_mu g) > 0, peripheral invertibility (tau != 0), and Wiener-Ikehara mu-rescaling verified.
[PASS] Prerequisites verified: renewal.md rational table (D <= 200), operator.md tilt ratio bound, and degree-law.md crossing identity.
[CAUGHT] Error 1: Wrong Yule weight n instead of (n-1) in Eq. (64)
[CAUGHT] Error 2: Wrong T_3->T_2 wait factor 1/(1+z) instead of 2/(2+z) in Eq. (5)
[CAUGHT] Error 3: Wrong Gamma(mu+1) instead of Gamma(mu+2) in c_k = d_k^* / Gamma(mu+2)
[CAUGHT] Error 4: Missing dx = mu dt Jacobian in Wiener-Ikehara x -> M_{x/mu} f_k
[CAUGHT] Error 5: Flipped residue sign +ell(K'_mu g) instead of -ell(K'_mu g) in Eq. (4)
[CAUGHT] Error 6: Omitting exterior parent slot in D_root(S_v)
[CAUGHT] Error 7: Lattice reproduction kernel producing poles at mu + 2*pi*i (tau != 0)
[CAUGHT] Error 8: Lowering q from 9/25 = 0.36 to 359/1000 in renewal.md Eq. (4) [breaks (c=1, D=3)]
ALL CHECKS AND 8/8 NEGATIVE CONTROLS PASSED.
```
