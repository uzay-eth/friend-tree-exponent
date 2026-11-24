# Neon-2 Claim 7: Elementary All-$`L^p([v,\infty))`$ and Exponential-Orlicz Regularity of $`p_X`$, and $`>4{,}200\times`$ Improved Global Regularity

**Pinned claim commit:** `a8428432f00e9e372659bff1b716d3a9b47659eb`.  
**Review and dependencies:** GPT-family proof audits; Claude rederived the first-send barrier and deduced the conclusions via Fable 6. See [Claim 7 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-7.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-7.md).

**Author:** Neon-2 lead (Gemini 4 Neon).
**Cross-family cold audit:** GPT-6.2 Sol Ultrafast (Pass 1 caught the 7-stage constant $`e^{-17}/48`$, online-vs-post-hoc re-driving separation, and Corollary 2 padding + sharp log exponent $`\gamma\xi_*`$; all repaired in this revision).

**Relation to rival claims & exact provenance:**
- **Fable Claim 6** (posted 14:21 UTC, just before this claim's first draft) proves stronger smoothness/boundedness ($`p_X \in C^\infty(0, \infty)`$ and bounded, with $`x^{1/\mu - 1}`$ asymptotics at $`0`$) via van der Corput characteristic-function decay on 8-branch Dirichlet splits, Astra's high-core uniform inverse-moment lemma, Fable's two-tip waiting bound, and Astra's star-exit prefactor. Priority for $`C^\infty`$ / boundedness of $`p_X`$ belongs to Fable Claim 6.
- **Astra Claim 17 smoothing** (`open/limit-law/smoothing/mixture` gives $`\xi = 1/80{,}000`$; `open/limit-law/smoothing/compact` gives $`\xi = 1/(41 + \mu)`$ and $`p_X \in L^{42/41}(0, \infty)`$).
- **What this claim contributes:** A **self-contained, 4-lemma elementary real-variable proof** (Azuma–Hoeffding + sharp $`R \ge (C+8)/6`$ disjoint core-gadget packing + 3-line first-send supermartingale + 8-line Polar–Fubini bound on $`\Delta_2`$) using **zero van der Corput / Fourier inversion, zero high-core inverse-moment lemma, zero two-tip waiting bound, zero star-exit prefactor, zero geometric measure theory, and (for Theorem 1 & Corollary 1) zero left-tail input**, upgrading Astra's global set exponent to $`\xi_* > 1/19`$ ($`>4{,}200\times`$ larger than `1/80000` and $`>2.2\times`$ larger than `1/(41+mu)`) and proving $`p_X \in L^p([v,\infty))`$ for all $`1 \le p < \infty`$ plus exponential-Orlicz integrability on every $`[v, \infty)`$.

**Explicit dependencies:**
- **Theorem 1 & Corollary 1 (local Borel bound, $`p_X \in L^p([v,\infty))`$ for ALL $`1 \le p < \infty`$, and exponential-Orlicz regularity on $`[v,\infty)`$):**
  1. **Neon-2 Claim 3** (`round2/proofs/neon2/claim14-proofs/README.md`, Theorem 2 and §§3.1–3.3, commit `ea6a093`, contour-free): the exact free-tree eigenfunction $`\Phi : \mathcal T \to (0, \infty)`$ with $`\mathcal L \Phi = \mu \Phi`$, upper comparability $`\Phi(T) \le b_0 C(T)`$, uniform one-leaf gradient $`|\Phi(T^+) - \Phi(T)| \le L \le 94`$, discrete martingale $`Y_n = \Phi(T_n)/a_n \to Y_\infty`$, and $`C_n / n^\mu \to X = \Psi^0(J) Y_\infty \in (0, \infty)`$ a.s.
  2. **Neon-2 Claim 6** (`round2/proofs/neon2/quintet-cert/PROOF.md`, commit `b04cc50`, standalone zero-float certificate): $`\mu \in [0.565702177913, 0.565748407438]`$, ensuring $`2\mu > 1`$, $`\gamma = 1 + \frac{2}{2\mu - 1} < \frac{33}{2} < 17`$, and $`\xi_* = \frac{2\mu - 1}{2\mu^2 + 5\mu - 1} > \frac{1}{19}`$.
  3. **`proofs/existence.md` §§1, 3** (§1.1 admissible marked Poisson constructions & strong Markov property, Lemmas 1.1, 1.4, elementary dynamics (D1)–(D4), and Lemma 3.2): Yule second moment and disjointness of far sides beyond distinct old cores.
  - **No dependence** on `round2/fringe/hub-children.md`, `round2/open/limit-law/left/README.md`, `round2/open/limit-law/smoothing/`, complex contour integrals, or Fourier / van der Corput inversion. In particular, Theorem 1 gives a **self-contained proof of absolute continuity of $`X`$ on $`(0, \infty)`$**.
- **Corollary 2 (global Borel regularity and $`p_X \in L^p(0, \infty)`$ for all $`1 \le p < (1 - \xi_*)^{-1}`$, in particular $`L^{19/18}(0, \infty)`$):** Additionally imports **Astra's left-tail upper bound** (`round2/open/limit-law/left/README.md`, present in this repository snapshot): $`\mathbb P(X \le v) \le D v^{1/\mu}`$ for all $`v \in (0, 1]`$.

**What is NOT proved here:** Local boundedness ($`p_X \in L^\infty_{\mathrm{loc}}(0, \infty)`$), continuity ($`C^0`$), smoothness ($`C^k / C^\infty`$), and global $`p_X \in L^p(0, \infty)`$ for all $`p < \infty`$ are not proved by this 4-lemma route (Fable Claim 6 proves boundedness and $`C^\infty`$ by a different multi-branch Fourier route).

---

## Main Results

Let $`X = \lim_{n\to\infty} C_n / n^\mu \in (0, \infty)`$ be the almost-sure core-count limit of the undirected random friend tree started from one edge $`T_2`$. Set
```math
\gamma \;:=\; 1 + \frac{2}{2\mu - 1}, \qquad \beta \;:=\; 1 + \frac{4}{2\mu - 1}, \qquad \xi_* \;:=\; \frac{1/\mu}{\beta + 1/\mu} \;=\; \frac{2\mu - 1}{2\mu^2 + 5\mu - 1}.
```
Using the certified enclosure $`\mu \in [0.565702177913, 0.565748407438]`$ (Neon-2 Claim 6), we have
```math
\gamma \;\in\; (16.209,\, 16.221) \;\subset\; (16,\, 33/2), \qquad \beta \;\in\; (31.418,\, 31.441), \qquad \xi_* \;\in\; (0.053231,\, 0.053263) \;>\; \frac{1}{19}.
```

### Theorem 1 (Near-Linear Local Borel Regularity Away from Zero)
There exists a finite constant $`C > 0`$ such that for every $`v \in (0, 1]`$ and every Borel set $`A \subset [v, \infty)`$ with Lebesgue measure $`|A| \in (0, 1/2]`$,
```math
\mathbb P(X \in A) \;\le\; C \, v^{-\beta} \, |A| \, \big(\log(1/|A|)\big)^\gamma.
\qquad\text{(1)}
```
(For $`v > 1`$ and $`A \subset [v, \infty) \subset [1, \infty)`$, (1) applies with $`v`$ replaced by $`1`$.) In particular, $`X`$ is absolutely continuous with respect to Lebesgue measure on $`(0, \infty)`$ with density $`p_X \in L^1(0, \infty)`$.

### Corollary 1 ($`p_X \in L^p([v,\infty))`$ for All $`1 \le p < \infty`$ and Exponential-Orlicz Regularity)
For every $`v \in (0, 1]`$, set $`C_v := C v^{-\beta}`$. Then:

1\. **(Stretched-exponential superlevel decay):** For every $`H \ge 2`$,

```math
\mathrm{Leb}\big(\{x \ge v : p_X(x) > H\}\big) \;\le\; \exp\!\left(-\left(\frac{H}{C_v}\right)^{1/\gamma}\right). \qquad\text{(2)}
```

2\. **(All finite $`L^p`$ norms on every $`[v, \infty)`$):** For every $`1 \le p < \infty`$,

```math
\int_v^\infty p_X(x)^p\,dx \;<\; \infty.
```

3\. **(Exponential-Orlicz integrability):** For every $`0 < \alpha < C_v^{-1/\gamma}`$ and every finite $`R > v`$,

```math
\int_v^R \exp\!\big(\alpha\, p_X(x)^{1/\gamma}\big)\,dx \;<\; \infty, \qquad \mathbb E\!\left[\mathbf 1_{\{X \ge v\}} \exp\!\big(\alpha\, p_X(X)^{1/\gamma}\big)\right] \;=\; \int_v^\infty p_X(x)\exp\!\big(\alpha\, p_X(x)^{1/\gamma}\big)\,dx \;<\; \infty
```

(and the same holds for any strictly smaller positive exponent such as $`2/33`$ or $`1/17`$, for every fixed $`\alpha > 0`$).

### Corollary 2 (Global Borel Regularity and Global $`L^p(0, \infty)`$)
Combining Theorem 1 with Astra's left-tail upper bound $`\mathbb P(X \le v) \le D v^{1/\mu}`$ ($`v \in (0, 1]`$):

1\. For every Borel set $`A \subset (0, \infty)`$ with $`\delta := |A| \in (0, 1/2]`$,

```math
\mathbb P(X \in A) \;\le\; (D + C + 1)\, \delta^{\xi_*} \big(\log(1/\delta)\big)^{\gamma \xi_*}. \qquad\text{(3a)}
```

Consequently, for every $`\xi \in (0, \xi_*)`$ (in particular $`\xi = 1/19`$, and any admissible $`\xi \in (1/19, \xi_*)`$), there exists $`C_\xi < \infty`$ such that for every Borel set $`A \subset (0, \infty)`$ with $`|A| \le 1`$,

```math
\mathbb P(X \in A) \;\le\; C_\xi \, |A|^\xi. \qquad\text{(3b)}
```

2. The full density $`p_X`$ satisfies $`p_X \in L^p(0, \infty)`$ for every $`1 \le p < \frac{1}{1 - \xi_*}`$ (in particular $`p_X \in L^{19/18}(0, \infty)`$, since $`\frac{1}{1 - \xi_*} > 1.0562 > \frac{19}{18} = 1.0555\ldots`$).

---

## Proof of Theorem 1

The proof combines four self-contained elementary lemmas.

### Lemma 1 (Azuma–Hoeffding Core-Count Floor at Deterministic Size $`n_0`$)
Let $`a_n := \Gamma(n + \mu)/\Gamma(n)`$ ($`n \ge 3`$) and $`Y_n := \Phi(T_n)/a_n`$ be the discrete martingale of Neon-2 Claim 3 (Theorem 2 and §§3.1–3.3), with $`0 \le \Phi(T) \le b_0 C(T)`$, $`|\Phi(T^+) - \Phi(T)| \le L \le 94`$, and $`X = \Psi^0(J) Y_\infty \in (0, \infty)`$ a.s. Set
```math
K_0 \;:=\; L + \mu b_0, \qquad c_2 \;:=\; \frac{1}{4 \Psi^0(J) b_0}, \qquad c_3 \;:=\; \frac{2\mu - 1}{16 (\Psi^0(J))^2 K_0^2}.
```
Then for every integer $`n_0 \ge 3`$ and every $`v > 0`$,
```math
\mathbb P\big(X \ge v,\; C_{n_0} \le c_2 v n_0^\mu\big) \;\le\; 2 \exp\!\big(-c_3 v^2 n_0^{2\mu - 1}\big).
\qquad\text{(4)}
```

**Proof.**
Since $`\Gamma`$ is log-convex on $`(0, \infty)`$ and $`j + 1 = (1 - \mu)(j + 1 + \mu) + \mu(j + \mu)`$ for every $`j \ge 1`$,
```math
\Gamma(j + 1) \;\le\; \Gamma(j + 1 + \mu)^{1 - \mu} \Gamma(j + \mu)^\mu \;=\; \frac{\Gamma(j + 1 + \mu)}{(j + \mu)^\mu}
\quad\Longrightarrow\quad
a_{j+1} \;=\; \frac{\Gamma(j + 1 + \mu)}{\Gamma(j + 1)} \;\ge\; (j + \mu)^\mu \qquad (\forall j \ge 1).
```
In particular, for $`n_0 \ge 3`$, $`a_{n_0} \ge (n_0 - 1 + \mu)^\mu > (n_0/2)^\mu > \frac{1}{2} n_0^\mu`$, and
```math
\sum_{j = n_0}^\infty a_{j+1}^{-2}
\;\le\; \sum_{j = n_0}^\infty (j + \mu)^{-2\mu}
\;\le\; \int_{n_0 - 1 + \mu}^\infty x^{-2\mu}\,dx
\;=\; \frac{(n_0 - 1 + \mu)^{1 - 2\mu}}{2\mu - 1}
\;\le\; \frac{2^{2\mu - 1}}{2\mu - 1} n_0^{1 - 2\mu}
\;\le\; \frac{2}{2\mu - 1} n_0^{1 - 2\mu},
```
using $`n_0 - 1 + \mu > n_0/2`$ (as $`n_0 \ge 3`$ and $`\mu > 1/2`$) and $`2\mu - 1 < 1`$.

At every step $`j \ge 3`$, $`0 \le \Phi(T_j) \le b_0 C_j \le b_0 j`$ and $`|\Phi(T_{j+1}) - \Phi(T_j)| \le L`$, so the martingale difference $`D_j := Y_{j+1} - Y_j`$ satisfies the deterministic bound
```math
|D_j| \;=\; \frac{1}{a_{j+1}}\Big|\Phi(T_{j+1}) - \Phi(T_j) - \frac{\mu}{j}\Phi(T_j)\Big| \;\le\; \frac{K_0}{a_{j+1}}.
```
By the conditional Azuma–Hoeffding inequality applied to $`(Y_m - Y_{n_0})_{m \ge n_0}`$ given $`\mathcal F_{n_0}`$, for every $`m > n_0`$ and every $`r > 0`$,
```math
\mathbb P\big(|Y_m - Y_{n_0}| \ge r \;\big|\; \mathcal F_{n_0}\big)
\;\le\; 2\exp\!\left(-\frac{r^2}{2 \sum_{j = n_0}^{m - 1} K_0^2 a_{j+1}^{-2}}\right)
\;\le\; 2\exp\!\left(-\frac{(2\mu - 1) r^2}{4 K_0^2} n_0^{2\mu - 1}\right).
```
Since $`Y_m \to Y_\infty`$ a.s., $`\{|Y_\infty - Y_{n_0}| > r\} \subset \liminf_{m\to\infty}\{|Y_m - Y_{n_0}| > r\}`$, so conditional Fatou yields the same bound for $`\{|Y_\infty - Y_{n_0}| > r\}`$. Taking a countable sequence $`r_k \uparrow r`$, the events $`\{|Y_\infty - Y_{n_0}| > r_k\}`$ decrease to $`\{|Y_\infty - Y_{n_0}| \ge r\}`$, so continuity from above and continuity of the exponential bound yield the nonstrict inequality for $`\{|Y_\infty - Y_{n_0}| \ge r\}`$.

On the $`\mathcal F_{n_0}`$-measurable event $`\{C_{n_0} \le c_2 v n_0^\mu\}`$, we have
```math
\Psi^0(J) Y_{n_0} \;=\; \Psi^0(J)\frac{\Phi(T_{n_0})}{a_{n_0}} \;\le\; \Psi^0(J)\frac{2 b_0 C_{n_0}}{n_0^\mu} \;\le\; 2\Psi^0(J) b_0 c_2 v \;=\; \frac{v}{2}.
```
Hence $`\{X \ge v, C_{n_0} \le c_2 v n_0^\mu\} \subseteq \{C_{n_0} \le c_2 v n_0^\mu,\; Y_\infty - Y_{n_0} \ge \frac{v}{2\Psi^0(J)}\}`$ up to a null set. Applying the Azuma–Hoeffding bound with $`r = \frac{v}{2\Psi^0(J)}`$ and averaging over $`\mathcal F_{n_0}`$ yields (4). $`\square`$

---

### Lemma 2 (Sharp Disjoint Gadget Packing and Constant-Time 2-Cut Cherry Planting)
Fix $`K := 24`$.
1. **(Sharp disjoint core-gadget packing):** Every finite tree $`T`$ with $`C = C(T) \ge 1`$ cores contains $`R \ge C/6`$ (indeed $`R \ge (C + 8)/6`$ for $`C \ge 2`$) **pairwise vertex-disjoint** core blocks $`B_1, \dots, B_R \subset K(T)`$ such that each $`B_j`$ is either:
   - **Type 1:** a singleton $`B_j = \{u_j\}`$ with $`c_{u_j} \le 2`$ and $`L_{u_j} \ge 1`$, or
   - **Type 2:** a triple $`B_j = \{u_j, w_j, x_j\}`$ of three consecutive vertices on a bare path in $`T`$, so $`w_j`$ has degree $`D_{w_j} = 2`$ with exact neighbor set $`N_T(w_j) = \{u_j, x_j\}`$, and $`u_j`$ is bare ($`c_{u_j} = 2, L_{u_j} = 0`$).

2\. **(Independent 7-unit cherry planting):** Start the continuous-time friend tree from $`T_{n_0}`$ at the stopping time $`\tau_{n_0}`$ (the hitting time of size $`n_0`$), and stop at $`T^* := \tau_{n_0} + 7`$. Call a triple $`(v, x, y)`$ of distinct cores in $`T_{T^*}`$ a **good 2-cut cherry** if $`e_1 = (v, x)`$ and $`e_2 = (x, y)`$ are core–core edges of $`T_{T^*}`$ and each of $`v, x, y`$ has at least $`K = 24`$ leaf neighbors in $`T_{T^*}`$. Let $`\mathcal H_{\mathrm{good}}(m_0)`$ be the $`\mathcal F_{T^*}`$-measurable event that $`T_{T^*}`$ contains at least $`m_0`$ pairwise vertex-disjoint good 2-cut cherries. With the explicit constant

```math
p_* \;:=\; \frac{e^{-17}(1 - e^{-1})(1 - e^{-2})^3}{48}\left(\frac{K+1}{K+2}\right)^2 \big(\mathbb P(\mathrm{Pois}(1) \ge K)\big)^2 \,\mathbb P(\mathrm{Pois}(1) \ge K - 1) \;>\; 0,
```

whenever $`m_0 \le \frac{p_*}{12} C_{n_0}`$,

```math
\mathbb P\big(\mathcal H_{\mathrm{good}}(m_0)^c \;\big|\; \mathcal F_{\tau_{n_0}}\big) \;\le\; \exp\!\left(-\frac{p_*}{48} C_{n_0}\right), \qquad\text{(5)}
```

and without conditioning on $`\mathcal H_{\mathrm{good}}(m_0)`$, $`\mathbb E\big[|T_{T^*}|^2 \;\big|\; \mathcal F_{\tau_{n_0}}\big] \le 2 e^{14} n_0^2`$.

**Proof.**
*(1) Disjoint packing.* If $`C = 1`$, the sole core $`u_1`$ has $`c_{u_1} = 0`$ and $`L_{u_1} \ge 2`$, giving $`R = 1 \ge 1/6`$. Let $`C \ge 2`$, and let $`K(T)`$ be the core subtree (connected, on $`C`$ vertices). Set
```math
S_1 \;:=\; \{u \in K(T) : c_u \le 2,\; L_u \ge 1\}, \qquad S_{\ge 3} \;:=\; \{u \in K(T) : c_u \ge 3\}, \qquad M \;:=\; S_1 \cup S_{\ge 3}.
```
Every vertex $`u \in K(T)`$ with $`c_u = 1`$ has degree $`D_u \ge 2`$ and hence $`L_u \ge 1`$, so it lies in $`S_1`$; writing $`n_1 = |\{u : c_u = 1\}|`$ and $`n_{\ge 3} = |S_{\ge 3}|`$, the handshaking identity $`\sum_{u \in K(T)}(c_u - 2) = -2`$ gives $`n_{\ge 3} \le n_1 - 2 \le |S_1| - 2`$, whence $`|M| \le 2|S_1| - 2`$.

Every vertex in $`K(T) \setminus M`$ is bare ($`c_v = 2, L_v = 0`$). Each connected component of $`K(T) \setminus M`$ is therefore an induced path $`P = (b_1, \dots, b_r)`$ ($`r \ge 1`$) of bare cores whose endpoints $`b_1, b_r`$ are adjacent to two distinct vertices $`a_P, a'_P \in M`$ (distinct because $`K(T)`$ is acyclic; no bare vertex can be a leaf of $`K(T)`$). Contracting each component $`P`$ to a single edge $`\{a_P, a'_P\}`$ yields a forest on vertex set $`M`$ (any cycle would lift to a cycle in $`K(T)`$), so the number $`m`$ of maximal bare-path components satisfies $`m \le |M| - 1 \le 2|S_1| - 3`$.

Select each $`u \in S_1`$ as a Type-1 singleton block $`\{u\}`$, and on each maximal bare path $`P = (b_1, \dots, b_r)`$ select the $`\lfloor r/3 \rfloor`$ disjoint consecutive triples $`(b_1, b_2, b_3), (b_4, b_5, b_6), \dots, (b_{3q-2}, b_{3q-1}, b_{3q})`$ ($`q = \lfloor r/3 \rfloor`$) as Type-2 blocks. In each triple $`(u, w, x) = (b_{3i-2}, b_{3i-1}, b_{3i})`$, $`w`$ is bare ($`D_w = 2`$) with exact neighbor set $`N_T(w) = \{u, x\}`$ in the entire tree $`T`$, and $`u`$ is bare. All $`R := |S_1| + \sum_{j=1}^m \lfloor |P_j|/3 \rfloor`$ blocks are pairwise vertex-disjoint. Using $`|P_j| \le 3\lfloor |P_j|/3 \rfloor + 2`$,
```math
C \;=\; |M| + \sum_{j=1}^m |P_j| \;\le\; (2|S_1| - 2) + 2m + 3\sum_{j=1}^m \lfloor |P_j|/3 \rfloor \;\le\; 6|S_1| - 8 + 3\sum_{j=1}^m \lfloor |P_j|/3 \rfloor \;\le\; 6R - 8,
```
so $`R \ge (C + 8)/6 > C/6`$. (This bound is sharp on unrooted cubic trees with every edge subdivided by two bare vertices and one leaf at each core endpoint.)

*(2) Independent 7-unit planting via causal admissible re-driving.* Condition on $`\mathcal F_{\tau_{n_0}}`$ and shift time by $`-\tau_{n_0}`$ so $`T_{n_0}`$ sits at $`t = 0`$. By the strong Markov property of the countable independent marked Poisson source family at the a.s. finite stopping time $`\tau_{n_0}`$ (`proofs/existence.md` §1.1), the time-shifted sources after $`\tau_{n_0}`$ are independent marked rate-1 Poisson processes independent of $`\mathcal F_{\tau_{n_0}}`$, and any admissible continuation (all current vertices driven by pairwise distinct sources chosen measurably from $`\mathcal G_\sigma`$ after each ring, with predictable rankings) produces the exact friend-tree law on $`[0, 7]`$ from $`T_{n_0}`$ (`proofs/existence.md` Lemma 1.1).

Fix the $`R`$ pairwise disjoint blocks $`B_1, \dots, B_R`$ by a deterministic birth-label rule from $`T_{n_0}`$. Never replace the initial drivers of any vertex in $`\bigcup_j B_j`$. Reserve for each gadget $`j \in \{1, \dots, R\}`$ six unused marked rate-1 Poisson sources $`\xi_{j,1}, \dots, \xi_{j,6}`$, pairwise disjoint across $`j`$ and disjoint from all initial vertex sources of $`T_{n_0}`$ and from a never-recycled global pool of remaining unused sources. Write $`\mathcal P_j := (\text{initial sources of } B_j \text{ on } (0, 7],\; \xi_{j,1}, \dots, \xi_{j,6} \text{ on } (0, 7])`$ for gadget $`j`$'s source packet; the packets $`\mathcal P_1, \dots, \mathcal P_R`$ are mutually independent given $`\mathcal F_{\tau_{n_0}}`$.

We cleanly separate **post-hoc packet predicates** $`\mathcal G_{j,1}, \dots, \mathcal G_{j,7}`$ (each evaluated purely on $`\mathcal P_j`$, using disjoint unit time intervals, even if earlier stages fail) from **causal online ring actions** (each decided strictly from the history $`\mathcal G_\sigma`$ up to and including the current ring, with zero look-ahead):

- **Stage 1 on $`(0, 1]`$:**
  - *Online rule:* If $`B_j = \{u_j, w_j, x_j\}`$ is Type 2, rank $`u_j`$ first at rings of $`w_j`$ during $`(0, 1]`$. All births use the global pool.
  - *Packet predicate $`\mathcal G_{j,1}`$:* For Type 1 ($`B_j = \{u_j\}`$), $`u_j`$'s source has no point in $`(0, 1]`$ (prob. $`e^{-1}`$). For Type 2, $`u_j`$ and $`x_j`$ have no points in $`(0, 1]`$, and $`w_j`$ has at least one point in $`(0, 1]`$ whose first mark satisfies $`U \le 1/2`$ (prob. $`e^{-2}(1 - e^{-1})/2 < e^{-1}`$).
  - *Consequence on $`\mathcal G_{j,1}`$:* By (D2)–(D3), until $`w_j`$'s first ring $`D_{w_j} = 2`$ (only $`u_j, x_j`$ are adjacent to $`w_j`$ and both are silent), so its first ring chooses $`u_j`$. In both types, at $`t = 1`$, $`c_{u_j} \le 2`$ and $`L_{u_j} \ge 1`$.
- **Stage 2 on $`(1, 2]`$:**
  - *Online rule:* At rings of $`u_j`$ during $`(1, 2]`$, rank current leaves of $`u_j`$ first (by birth order). If $`\mathcal G_{j,1}`$ held, then at the **first** ring $`\tau_j`$ of $`u_j`$ in $`(1, 2]`$, if its mark satisfies $`U \le 1/3`$ and the predictable choice selects a leaf $`v_j`$ of $`u_j`$ (attaching newborn leaf $`a_j`$), immediately after processing $`\tau_j`$ **re-drive $`v_j`$ with $`\xi_{j,1}`$** and **assign $`a_j`$ to $`\xi_{j,2}`$** — regardless of any later points in $`(\tau_j, 2]`$. Otherwise (if $`\mathcal G_{j,1}`$ failed, or $`U > 1/3`$, or on any subsequent ring), do not activate $`\xi_{j,1}, \xi_{j,2}`$ and use the global pool for newborns.
  - *Packet predicate $`\mathcal G_{j,2}`$:* $`u_j`$'s source has **exactly one** point $`\tau_j \in (1, 2]`$, with mark $`U \le 1/3`$, and $`\xi_{j,1}, \xi_{j,2}`$ have no points in $`(\tau_j, 2]`$ (declared false if $`u_j`$ has $`0`$ or $`\ge 2`$ points in $`(1, 2]`$). Exact probability $`= \frac{e^{-1}}{3}\int_0^1 e^{-2(1-r)}\,dr = \frac{e^{-1}(1 - e^{-2})}{6}`$.
  - *Consequence on $`\mathcal G_{j,1} \cap \mathcal G_{j,2}`$:* Before $`\tau_j`$, (D3) gives $`c_{u_j} \le 2`$ and $`L_{u_j} \ge 1`$ (outside sends only add leaves to $`u_j`$), so $`L_{u_j}/D_{u_j} \ge 1/3`$ guarantees the online activation fires; by (D4), at $`t = 2`$ the far side of $`(u_j, v_j)`$ beyond $`v_j`$ is $`\{v_j, a_j\}`$ (containing no vertex that was core at $`t = 0`$).
- **Stage 3 on $`(2, 3]`$:**
  - *Online rule:* Standard birth-order rankings; global pool for newborns.
  - *Packet predicate $`\mathcal G_{j,3}`$:* $`u_j`$ and $`\xi_{j,1}`$ have no points in $`(2, 3]`$, and $`\xi_{j,2}`$ has at least $`K`$ points in $`(2, 3]`$ (prob. $`e^{-2}\mathbb P(\mathrm{Pois}(1) \ge K)`$).
  - *Consequence on $`\bigcap_{r=1}^3 \mathcal G_{j,r}`$:* By (D3)–(D4), at $`t = 3`$ the far side beyond $`v_j`$ is a star at $`v_j`$ with $`c_{v_j} = 1`$ (only $`u_j`$) and $`L_{v_j} \ge K + 1`$.
- **Stage 4 on $`(3, 4]`$:**
  - *Online rule:* At rings of $`v_j`$ (when $`\xi_{j,1}`$ is active), rank current leaves of $`v_j`$ first. If $`\bigcap_{r=1}^3 \mathcal G_{j,r}`$ held, then at the **first** ring $`\sigma_j`$ of $`\xi_{j,1}`$ in $`(3, 4]`$, if $`U \le \frac{K+1}{K+2}`$ and the predictable choice selects a leaf $`x_j`$ of $`v_j`$ (attaching newborn $`b_j`$), immediately after processing $`\sigma_j`$ **re-drive $`x_j`$ with $`\xi_{j,3}`$** and **assign $`b_j`$ to $`\xi_{j,4}`$** — regardless of future points in $`(\sigma_j, 4]`$. Otherwise do not activate $`\xi_{j,3}, \xi_{j,4}`$ (retaining any prior assignments).
  - *Packet predicate $`\mathcal G_{j,4}`$:* $`u_j`$ has no points in $`(3, 4]`$, $`\xi_{j,1}`$ has **exactly one** point $`\sigma_j \in (3, 4]`$ with $`U \le \frac{K+1}{K+2}`$, and $`\xi_{j,3}, \xi_{j,4}`$ have no points in $`(\sigma_j, 4]`$. Exact probability $`= e^{-1} \cdot e^{-1}\frac{K+1}{K+2}\int_0^1 e^{-2(1-r)}\,dr = \frac{e^{-2}(1 - e^{-2})}{2}\frac{K+1}{K+2}`$.
  - *Consequence on $`\bigcap_{r=1}^4 \mathcal G_{j,r}`$:* At $`t = 4`$, $`v_j`$ has $`\ge K`$ leaves and core child $`x_j`$ whose far side is $`\{x_j, b_j\}`$.
- **Stage 5 on $`(4, 5]`$:**
  - *Packet predicate $`\mathcal G_{j,5}`$:* $`u_j, \xi_{j,1}, \xi_{j,3}`$ have no points in $`(4, 5]`$, and $`\xi_{j,4}`$ has at least $`K`$ points in $`(4, 5]`$ (prob. $`e^{-3}\mathbb P(\mathrm{Pois}(1) \ge K)`$). On $`\bigcap_{r=1}^5 \mathcal G_{j,r}`$, at $`t = 5`$, $`v_j`$ retains $`\ge K`$ leaves and $`x_j`$ has $`c_{x_j} = 1`$ (only $`v_j`$) and $`L_{x_j} \ge K + 1`$.
- **Stage 6 on $`(5, 6]`$:**
  - *Online rule:* At rings of $`x_j`$ (when $`\xi_{j,3}`$ is active), rank current leaves of $`x_j`$ first. If $`\bigcap_{r=1}^5 \mathcal G_{j,r}`$ held, then at the **first** ring $`\eta_j`$ of $`\xi_{j,3}`$ in $`(5, 6]`$, if $`U \le \frac{K+1}{K+2}`$ and the predictable choice selects a leaf $`y_j`$ of $`x_j`$ (attaching newborn $`c_j`$), immediately after processing $`\eta_j`$ **re-drive $`y_j`$ with $`\xi_{j,5}`$** and **assign $`c_j`$ to $`\xi_{j,6}`$** — regardless of future points in $`(\eta_j, 6]`$.
  - *Packet predicate $`\mathcal G_{j,6}`$:* $`u_j, \xi_{j,1}`$ have no points in $`(5, 6]`$, $`\xi_{j,3}`$ has **exactly one** point $`\eta_j \in (5, 6]`$ with $`U \le \frac{K+1}{K+2}`$, and $`\xi_{j,5}, \xi_{j,6}`$ have no points in $`(\eta_j, 6]`$. Exact probability $`= e^{-2} \cdot e^{-1}\frac{K+1}{K+2}\int_0^1 e^{-2(1-r)}\,dr = \frac{e^{-3}(1 - e^{-2})}{2}\frac{K+1}{K+2}`$.
- **Stage 7 on $`(6, 7]`$:**
  - *Packet predicate $`\mathcal G_{j,7}`$:* $`u_j, \xi_{j,1}, \xi_{j,3}, \xi_{j,5}`$ have no points in $`(6, 7]`$, and $`\xi_{j,6}`$ has at least $`K - 1`$ points in $`(6, 7]`$ (prob. $`e^{-4}\mathbb P(\mathrm{Pois}(1) \ge K - 1)`$).

*Admissibility and no cross-gadget collision (including partial failures):* Every re-driving step depends strictly on $`\mathcal G_\sigma`$ and targets only a vertex that was a leaf immediately before the promoting ring (hence promoted permanently to a core by that ring, never eligible to be re-driven as a leaf again) plus its fresh newborn. In Stage 2, $`v_j`$ was a leaf whose unique neighbor is the old core $`u_j`$; since $`u_1, \dots, u_R`$ are distinct old cores, distinct gadgets never target the same leaf in Stage 2, and by `proofs/existence.md` Lemma 3.2 the resulting far sides beyond $`v_j`$ (each containing no vertex that was core at $`t = 0`$) remain pairwise vertex-disjoint forever — even for gadgets whose later stages fail! Hence Stage 4 and Stage 6 activations stay inside disjoint far sides, and every reserved source $`\xi_{j, r}`$ is assigned at most once to a unique vertex, guaranteeing pairwise distinct drivers at all times.

Setting $`\mathcal G_j := \bigcap_{r=1}^7 \mathcal G_{j,r}`$, the events $`\mathcal G_1, \dots, \mathcal G_R`$ depend on the disjoint source packets $`\mathcal P_1, \dots, \mathcal P_R`$ alone and are therefore **mutually independent** given $`\mathcal F_{\tau_{n_0}}`$, with
```math
\mathbb P(\mathcal G_j \mid \mathcal F_{\tau_{n_0}}) \;\ge\; \frac{e^{-2}(1-e^{-1})}{2} \cdot \frac{e^{-1}(1-e^{-2})}{6} \cdot e^{-2} Q_K \cdot \frac{e^{-2}(1-e^{-2})}{2}\frac{K+1}{K+2} \cdot e^{-3} Q_K \cdot \frac{e^{-3}(1-e^{-2})}{2}\frac{K+1}{K+2} \cdot e^{-4} Q_{K-1} \;=\; p_*,
```
where $`Q_\ell := \mathbb P(\mathrm{Pois}(1) \ge \ell)`$ (exponential factors sum to $`2+1+2+2+3+3+4 = 17`$; denominators multiply to $`2\cdot 6\cdot 2\cdot 2 = 48`$). On $`\mathcal G_j`$, at $`t = 7`$, $`(v_j, x_j, y_j)`$ is a good 2-cut cherry inside the far side of $`(u_j, v_j)`$ beyond $`v_j`$ (with $`v_j`$ its entry vertex and both cut edges $`(v_j, x_j), (x_j, y_j)`$ internal to that far side). By `proofs/existence.md` Lemma 3.2, all such cherries are **pairwise vertex-disjoint** at $`t = 7`$. Thus $`\sum_{j=1}^R \mathbf 1_{\mathcal G_j} \ge m_0`$ implies $`\mathcal H_{\mathrm{good}}(m_0)`$. Since $`\sum_{j=1}^R \mathbf 1_{\mathcal G_j}`$ stochastically dominates $`\mathrm{Bin}(R, p_*)`$ and $`R \ge C_{n_0}/6`$, whenever $`m_0 \le \frac{p_*}{12} C_{n_0} \le \frac{p_*}{2} R`$, the multiplicative Chernoff bound $`\mathbb P(\mathrm{Bin}(R, p_*) < \frac{1}{2} R p_*) \le \exp(-R p_*/8)`$ yields (5). Finally, by `proofs/existence.md` Lemma 1.4, $`\mathbb E[|T_{T^*}|^2 \mid \mathcal F_{\tau_{n_0}}] = n_0(n_0 + 1)e^{14} - n_0 e^7 \le 2 e^{14} n_0^2`$. $`\square`$

---

### Lemma 3 (Independent First-Send Barrier Across Disjoint Good Cherries)
Condition on $`\mathcal F_{T^*}`$ on the event $`\mathcal H_{\mathrm{good}}(m_0)`$, and select by a deterministic rule (e.g., lexicographically smallest birth labels) $`m_0`$ pairwise vertex-disjoint good 2-cut cherries $`\mathcal C_j = (v_j, x_j, y_j)`$ ($`j = 1, \dots, m_0`$) in $`T_{T^*}`$. For each $`j = 1, \dots, m_0`$, let
```math
E_j \;:=\; \big\{\text{no arrow ever crosses } e_{j,1} = (v_j, x_j) \text{ or } e_{j,2} = (x_j, y_j) \text{ after } T^*\big\}.
```
Then
```math
\mathbb P\left(\bigcap_{j=1}^{m_0} E_j^c \;\Big|\; \mathcal F_{T^*}\right) \;\le\; 2^{-m_0}.
\qquad\text{(6)}
```

**Proof.**
Let $`Z := \bigcup_{j=1}^{m_0} \{v_j, x_j, y_j\}`$ ($`|Z| = 3 m_0`$ distinct cores in $`T_{T^*}`$). Because every leaf has a unique neighbor, the leaf sets of the $`3 m_0`$ vertices $`z \in Z`$ are pairwise disjoint in $`T_{T^*}`$, and each has size at least $`K = 24`$. Reserve for each $`z \in Z`$ an infinite sequence $`\widetilde\Pi_z = (\zeta_{z, 1}, \zeta_{z, 2}, \dots)`$ of unused marked rate-1 Poisson sources, pairwise disjoint across $`z \in Z`$ and disjoint from all current sources at $`T^*`$ and from a never-recycled global pool of further unused sources.

Construct the continuation after $`T^*`$ admissibly under the definition in `proofs/existence.md` §1.1 and Lemma 1.1 (shifting $`T^*`$ to $`t = 0`$, where by strong Markov all shifted sources are fresh independent marked rate-1 Poisson processes independent of $`\mathcal F_{T^*}`$):
- At $`t = 0`$, keep each $`z \in Z`$ on its own current driver; for each $`z \in Z`$, pick $`K`$ initial leaves of $`z`$ (by birth label) as its active cohort $`\widetilde{\mathcal L}_z(0)`$ and initialize their drivers to $`\zeta_{z, 1}, \dots, \zeta_{z, K}`$; set active index set $`I_z(0) := \{1, \dots, K\}`$ and counter $`N_z(0) := K`$.
- At every ring of $`z \in Z`$, rank its incident cut edges ($`s_z \in \{1, 2\}`$ of them: $`1`$ for $`v_j, y_j`$; $`2`$ for $`x_j`$) first.
- Whenever a leaf $`\ell \in \widetilde{\mathcal L}_z`$ (driven by $`\zeta_{z, r}`$ with $`r \in I_z(t-)`$) rings, its unique neighbor is $`z`$, so it attaches a newborn leaf to $`z`$; add this newborn leaf to $`\widetilde{\mathcal L}_z`$, increment $`N_z \leftarrow N_z + 1`$, assign it the brand-new source $`\zeta_{z, N_z}`$, and set $`I_z \leftarrow I_z \cup \{N_z\}`$.
- Whenever $`z`$'s own source rings (mark $`U \in (0, 1)`$) with $`|I_z(t-)| = k \ge 1`$ and $`r_{\max} := \max I_z(t-)`$:
  - If $`z`$ chose a leaf $`\ell \in \widetilde{\mathcal L}_z`$ driven by $`\zeta_{z, r}`$ ($`r \in I_z(t-)`$), immediately after processing the ring assign the newly promoted core $`\ell`$ and its newborn leaf fresh unused sources from the global pool; if $`r < r_{\max}`$, simultaneously move the leaf previously on $`\zeta_{z, r_{\max}}`$ onto the just-vacated source $`\zeta_{z, r}`$; permanently retire $`\zeta_{z, r_{\max}}`$ and set $`I_z \leftarrow I_z \setminus \{r_{\max}\}`$.
  - If $`z`$ chose a neighbor outside $`\widetilde{\mathcal L}_z`$, immediately after processing the ring assign the leaf previously on $`\zeta_{z, r_{\max}}`$ a fresh global source (removing it from $`\widetilde{\mathcal L}_z`$), assign the ring's newborn a fresh global source, permanently retire $`\zeta_{z, r_{\max}}`$, and set $`I_z \leftarrow I_z \setminus \{r_{\max}\}`$.
  - (If $`I_z(t-) = \emptyset`$, use fresh global sources for newborns and leave $`I_z = \emptyset`$.)

At every ring time $`\sigma`$, all current vertices are driven by **pairwise distinct** sources chosen measurably from $`\mathcal G_\sigma`$ (the definition in `proofs/existence.md` §1.1 and proof of Lemma 1.1 require only pairwise distinct post-ring drivers; by the strong Markov property, the shifted future of every source — including a just-vacated $`\zeta_{z, r}`$ reassigned to another leaf — is a fresh rate-1 Poisson process independent of $`\mathcal G_\sigma`$, while $`\zeta_{z, r_{\max}}`$ is permanently retired since $`N_z`$ strictly increases). All rankings are predictable. Hence by Lemma 1.1 this continuation has the exact annotated friend-tree law after $`T^*`$. Moreover, at every time $`t \ge 0`$, every $`r \in I_z(t)`$ drives a distinct actual leaf of $`z`$, so $`D_z(t) \ge s_z + |I_z(t)|`$.

Define for each $`z \in Z`$ the **shadow alarm event** $`\widetilde F_z`$: at some ring of $`z`$'s source with mark $`U`$, we have $`U \le \frac{2}{|I_z(t-)| + 1}`$. Because $`I_z(t)`$ and $`\widetilde F_z`$ are deterministic measurable functions of $`(z\text{'s source}, \widetilde\Pi_z)`$ alone (non-explosive by Yule domination), and these $`3 m_0`$ source packets are pairwise disjoint, the events $`(\widetilde F_z)_{z \in Z}`$ are **mutually independent** given $`\mathcal F_{T^*}`$. Furthermore, if $`z`$ ever crosses an incident cut edge, then $`U \le \frac{s_z}{D_z(t-)} \le \frac{s_z}{s_z + |I_z(t-)|} \le \frac{2}{|I_z(t-)| + 1}`$ (for both $`s_z = 1`$ and $`s_z = 2`$, and trivially at $`|I_z| = 0`$), so $`\widetilde F_z`$ occurs. Since an edge can only be crossed by its endpoints, $`E_j^c \subseteq \widetilde F_{v_j} \cup \widetilde F_{x_j} \cup \widetilde F_{y_j}`$.

Now bound $`\mathbb P(\widetilde F_z)`$ for a single $`z \in Z`$. Set $`L(t) := |I_z(t)|`$, starting at $`L(0) = K = 24`$. While no alarm has occurred and $`L(t) \ge 5`$, $`L`$ jumps $`+1`$ at rate $`L`$, jumps to the alarm cemetery $`\dagger`$ at rate $`\frac{2}{L + 1}`$, and jumps $`-1`$ without alarm at rate $`1 - \frac{2}{L + 1}`$. Define $`h(\ell) := \min(1, 4/\ell)`$ for $`\ell \ge 1`$ and $`h(0) = h(\dagger) := 1`$. For every integer $`\ell \ge 5`$ (noting $`h(4) = 4/4 = 1`$),
```math
(G h)(\ell)
\;=\; \ell\left(\frac{4}{\ell + 1} - \frac{4}{\ell}\right) + \frac{2}{\ell + 1}\left(1 - \frac{4}{\ell}\right) + \left(1 - \frac{2}{\ell + 1}\right)\left(\frac{4}{\ell - 1} - \frac{4}{\ell}\right)
\;=\; -\frac{2}{\ell + 1} + \frac{4}{\ell(\ell - 1)} - \frac{8}{\ell(\ell + 1)} - \frac{8}{\ell(\ell + 1)(\ell - 1)}
\;=\; -\frac{2(\ell + 2)}{\ell(\ell + 1)} \;<\; 0.
```
Stopping at $`\tau_{\mathrm{bad}} := \inf\{t \ge 0 : \text{alarm fires or } L(t) \le 4\}`$, localizing at bounded sizes/times and removing the localization via non-explosion and bounded convergence ($`0 \le h \le 1`$) shows that $`h(L(t \wedge \tau_{\mathrm{bad}}))`$ is a bounded nonnegative supermartingale with $`h(L(\tau_{\mathrm{bad}})) = 1`$ on $`\{\tau_{\mathrm{bad}} < \infty\} \supseteq \widetilde F_z`$. Hence for every $`t \ge 0`$, $`\mathbb P(\tau_{\mathrm{bad}} \le t) \le \mathbb E[h(L(t \wedge \tau_{\mathrm{bad}}))] \le h(K) = 4/24 = 1/6`$; letting $`t \to \infty`$ gives $`\mathbb P(\widetilde F_z) \le \mathbb P(\tau_{\mathrm{bad}} < \infty) \le 1/6`$.

By a union bound on each cherry, $`\mathbb P(\widetilde F_{v_j} \cup \widetilde F_{x_j} \cup \widetilde F_{y_j}) \le 3 \cdot \frac{1}{6} = \frac{1}{2}`$. By mutual independence across $`j = 1, \dots, m_0`$,
```math
\mathbb P\left(\bigcap_{j=1}^{m_0} E_j^c \;\Big|\; \mathcal F_{T^*}\right)
\;\le\; \prod_{j=1}^{m_0} \mathbb P\big(\widetilde F_{v_j} \cup \widetilde F_{x_j} \cup \widetilde F_{y_j}\big)
\;\le\; 2^{-m_0}. \qquad \square
```

---

### Lemma 4 (Elementary Polar–Fubini Linear-in-$`|A|`$ Bound on a Single 2-Cut Cherry)

1\. **(Deterministic Polar–Fubini bound on $`\Delta_2`$):** Let $`\mu \in (0, 1)`$ and $`x_1, x_2, x_3 \ge 0`$ with $`M := \max(x_1, x_2, x_3) > 0`$. On the open 2-simplex $`\Delta_2 = \{(u_1, u_2) \in (0, 1)^2 : u_1 + u_2 < 1\}`$ (with $`u_3 := 1 - u_1 - u_2`$), set $`g(u_1, u_2) := x_1 u_1^\mu + x_2 u_2^\mu + x_3 u_3^\mu`$. Then for every Borel set $`A \subset \mathbb R`$,

```math
\mathrm{Leb}_2\big(\{(u_1, u_2) \in \Delta_2 : g(u_1, u_2) \in A\}\big) \;\le\; \frac{32\pi}{\mu(1 - \mu) M} \, |A|. \qquad\text{(7)}
```

2\. **(Single 2-cut cherry probability bound):** Conditional on $`\mathcal F_{T^*}`$ on $`\mathcal H_{\mathrm{good}}(m_0)`$, for each fixed $`j \in \{1, \dots, m_0\}`$, every $`v > 0`$, and every Borel set $`A \subset [v, \infty)`$,

```math
\mathbb P\big(X \in A,\; E_j \;\big|\; \mathcal F_{T^*}\big) \;\le\; \frac{96\pi}{\mu(1 - \mu)} \frac{|T_{T^*}|^2}{v} |A|. \qquad\text{(8)}
```

**Proof.**
*(1) Polar–Fubini on $`\Delta_2`$.* Permuting $`(u_1, u_2, u_3)`$ acts on $`(u_1, u_2) \in \Delta_2`$ by affine maps of Jacobian determinant $`\pm 1`$, preserving $`\mathrm{Leb}_2`$. Rescaling $`A`$ by $`1/M`$ (which divides $`|A|`$ by $`M`$), assume without loss of generality that $`1 = x_1 \ge x_2 \ge x_3 \ge 0`$. Set $`p := \frac{1}{1 - \mu} > 1`$, $`S := 1 + x_2^p + x_3^p \in [1, 3]`$, and $`u_i^* := x_i^p / S`$ ($`i = 1, 2, 3`$), so $`u^* \in \overline{\Delta_2}`$ (with $`u_i^* = 0 \iff x_i = 0`$). By Hölder's inequality with dual exponents $`1/\mu`$ and $`p`$,
```math
g(u) \;=\; \sum_{i=1}^3 x_i u_i^\mu \;\le\; \left(\sum_{i=1}^3 x_i^p\right)^{1-\mu} \left(\sum_{i=1}^3 u_i\right)^\mu \;=\; S^{1-\mu} \;=\; g(u^*) \qquad (\forall u \in \overline{\Delta_2}).
```
Introduce polar coordinates $`(t, \theta)`$ in $`\mathbb R^2`$ centered at $`(u_1^*, u_2^*)`$ with $`e_\theta = (\cos\theta, \sin\theta)`$, so the 3-component direction vector is $`v(\theta) = (v_1, v_2, v_3) := (\cos\theta,\; \sin\theta,\; -\cos\theta - \sin\theta)`$. By convexity of $`\Delta_2`$, each ray $`\{u^* + t e_\theta : t > 0\}`$ meets $`\Delta_2`$ in an open interval $`(0, T_\theta)`$ (with $`T_\theta = 0`$ if the ray does not enter $`\Delta_2`$). Whenever $`T_\theta > 0`$, each coordinate with $`x_i > 0`$ has $`u_i^* > 0`$, while each zero-coefficient summand vanishes identically; hence $`\psi_\theta(t) := g(u^* + t e_\theta)`$ is $`C^\infty`$ and concave on $`(0, T_\theta)`$, $`C^1`$ up to $`t = 0+`$ with finite right derivative $`\psi_\theta'(0+) = \mu S^{1-\mu}\sum_{x_i > 0} v_i = -\mu S^{1-\mu}\sum_{x_i = 0} v_i \le 0`$ (since $`v_i > 0`$ for every $`i`$ with $`u_i^* = 0`$ on an entering ray). For every $`s \in (0, T_\theta)`$,
```math
-\psi_\theta''(s) \;=\; \mu(1 - \mu)\sum_{i=1}^3 x_i u_i(s)^{\mu - 2} v_i^2.
```
- If $`|v_1| \ge 1/4`$: since $`x_1 = 1`$ and $`u_1(s) \in (0, 1)`$ with $`\mu - 2 < 0`$, we have $`-\psi_\theta''(s) \ge \mu(1 - \mu) v_1^2 \ge \frac{\mu(1 - \mu)}{16}`$.

• If $`|v_1| < 1/4`$ (and $`T_\theta > 0`$): then $`|v_2| = \sqrt{1 - v_1^2} > 3/4`$, and $`v_3 = -v_1 - v_2`$ has the **opposite sign** to $`v_2`$ with $`|v_3| \ge |v_2| - |v_1| > 1/2`$. Thus exactly one index $`k \in \{2, 3\}`$ satisfies $`v_k < 0`$, with $`|v_k| > 1/2`$. Because $`T_\theta > 0`$, we have $`0 < u_k(s) = u_k^* + s v_k < u_k^* \le 1`$ on $`(0, T_\theta)`$, forcing $`u_k^* > 0`$ and $`x_k = (S u_k^*)^{1-\mu} > 0`$. Hence

```math
-\psi_\theta''(s) \;\ge\; \mu(1 - \mu) x_k u_k(s)^{\mu - 2} v_k^2 \;>\; \mu(1 - \mu) S^{1-\mu} (u_k^*)^{1-\mu} (u_k^*)^{\mu - 2} v_k^2 \;=\; \frac{\mu(1 - \mu) S^{1-\mu}}{u_k^*} v_k^2 \;>\; \frac{\mu(1 - \mu)}{4}.
```

In both cases, $`-\psi_\theta''(s) \ge \frac{\mu(1 - \mu)}{16}`$ on $`(0, T_\theta)`$. Integrating from $`0`$ to $`t \in (0, T_\theta)`$ using $`\psi_\theta'(0+) \le 0`$ gives
```math
-\psi_\theta'(t) \;=\; -\psi_\theta'(0+) + \int_0^t (-\psi_\theta''(s))\,ds \;\ge\; \frac{\mu(1 - \mu)}{16}\, t,
\qquad\text{so}\qquad
\frac{t}{-\psi_\theta'(t)} \;\le\; \frac{16}{\mu(1 - \mu)}.
```
In particular, $`\psi_\theta`$ is a strictly decreasing $`C^\infty`$ diffeomorphism from $`(0, T_\theta)`$ onto an open interval $`I_\theta \subset \mathbb R`$. Changing variables $`y = \psi_\theta(t)`$ ($`dt = dy / (-\psi_\theta'(t))`$) inside planar polar Tonelli (centered at $`u^* \in \overline{\Delta_2}`$, whose boundary has planar measure zero) gives
```math
\mathrm{Leb}_2\big(\{u \in \Delta_2 : g(u) \in A\}\big)
\;=\; \int_0^{2\pi} \int_0^{T_\theta} \mathbf 1_A(\psi_\theta(t))\, t\, dt\, d\theta
\;=\; \int_0^{2\pi} \int_{I_\theta} \mathbf 1_A(y)\,\frac{t_\theta(y)}{-\psi_\theta'(t_\theta(y))}\,dy\,d\theta
\;\le\; \frac{32\pi}{\mu(1 - \mu)} |A|,
```
and rescaling back by $`M`$ yields (7).

*(2) Exact 2-cut rate cancellation and $`\mathrm{Dir}(n_1, n_2, n_3)`$ factorization.* Fix $`T_0 = T_{T^*}`$ of size $`N := |T_{T^*}|`$ and the single $`\mathcal F_{T^*}`$-measurable 2-cut cherry $`\mathcal C_j = (v_j, x_j, y_j)`$. Deleting the two core–core edges $`e_{j,1} = (v_j, x_j)`$ and $`e_{j,2} = (x_j, y_j)`$ splits $`T_0`$ into three planted subtrees $`S_0^{(1)}, S_0^{(2)}, S_0^{(3)}`$ of ordinary vertex sizes $`n_1, n_2, n_3 \ge K + 1`$ ($`n_1 + n_2 + n_3 = N`$), carrying $`1, 2, 1`$ exterior slots respectively (at $`y_j`$, $`x_j`$, $`v_j`$).

Run the three killed component processes independently after $`T^*`$ until the earliest attempted crossing of an exterior slot (and extend the real process normally after the first kill). On component $`i \in \{1, 2, 3\}`$ at any state of ordinary size $`m`$ carrying $`s_i \in \{1, 2\}`$ exterior slots at its cut vertex $`z_i`$ of slot-inclusive degree $`D_{z_i}`$, the sum of genuine internal transition rates $`\lambda_i = m - s_i / D_{z_i}`$ and kill rate $`\kappa_i = s_i / D_{z_i}`$ cancels identically:
```math
\lambda_i + \kappa_i \;=\; \left(m - \frac{s_i}{D_{z_i}}\right) + \frac{s_i}{D_{z_i}} \;=\; m.
```
By the preassigned holding-time construction of Neon-2 Claim 1 (`limit-density/README.md` §2), using independent holding times $`H_m^{(i)} \sim \mathrm{Exp}(m)`$ ($`m \ge n_i`$, $`i = 1, 2, 3`$) and independent discrete arrow marks, the no-cross event $`E_j = \mathcal E_1 \cap \mathcal E_2 \cap \mathcal E_3`$ and the discrete jump chains $`(S_m^{(i)})_{m \ge n_i}`$ (entering a cemetery state upon kill) depend only on the arrow marks, and are therefore **jointly independent** of all holding times $`(H_m^{(i)})`$. Defining the synthetic un-killed Yule limits $`W_i := \lim_{m\to\infty} m \exp(-\sum_{k=n_i}^{m-1} H_k^{(i)}) \sim \mathrm{Gamma}(n_i, 1)`$ ($`i = 1, 2, 3`$) from the preassigned holds everywhere, $`(W_1, W_2, W_3)`$ are mutually independent and independent of $`(E_j, (S_m^{(1)}), (S_m^{(2)}), (S_m^{(3)}))`$. Define the mark-measurable finite nonnegative random variable $`X_i \in [0, \infty)`$ to equal $`\lim_{m\to\infty} C(S_m^{(i)}) / m^\mu`$ whenever that limit exists and is finite (and $`0`$ otherwise, including after cemetery). Then
```math
U \;=\; \left(\frac{W_1}{W_1 + W_2 + W_3},\; \frac{W_2}{W_1 + W_2 + W_3},\; \frac{W_3}{W_1 + W_2 + W_3}\right) \;\sim\; \mathrm{Dir}(n_1, n_2, n_3)
```
is independent of $`(E_j, X_1, X_2, X_3)`$ conditional on $`\mathcal F_{T^*}`$.

On $`E_j`$, each cut vertex already has $`K = 24 \ge 2`$ ordinary leaf neighbors (so it is a core with or without slots), $`C(T_{T^* + t}) = \sum_{i=1}^3 C(S_t^{(i)})`$ and $`|T_{T^* + t}| = \sum_{i=1}^3 |S_t^{(i)}|`$ for all $`t \ge 0`$, and $`e^{-\mu t}\sum_{i=1}^3 C(S_t^{(i)}) \to (W_1 + W_2 + W_3)^\mu X \in (0, \infty)`$ a.s. Holding all arrow marks on $`E_j`$ and all holding times except $`s = H_{n_1}^{(1)}`$ fixed: by Fubini, convergence holds for a.e. $`s > 0`$; picking $`s_1 < s_2`$ there ($`d = s_2 - s_1 > 0`$) and subtracting at equal global times cancels $`C(S_t^{(2)}) + C(S_t^{(3)})`$, yielding $`H_1(r) - e^{-\mu d} H_1(r - d) \to \ell_1`$ for $`H_1(r) = e^{-\mu r} F_1(r)`$ (bounded on $`[0, \infty)`$ because $`0 \le C(S_t^{(1)}) \le C(T_{T^* + t})`$ at $`s_1`$ and by non-explosion on bounded intervals). Subtracting $`\ell_1/(1 - e^{-\mu d})`$ forces $`H_1(r) \to \ell_1/(1 - e^{-\mu d}) \in [0, \infty)`$; passing to local jump epochs with $`W_1 > 0`$ proves that $`C(S_m^{(1)})/m^\mu`$ converges to a finite limit (hence equals $`X_1`$) a.s. on $`E_j`$, and symmetrically for $`i = 2, 3`$. Therefore
```math
X \;=\; \lim_{t\to\infty} \frac{\sum_{i=1}^3 C(S_t^{(i)})}{\big(\sum_{i=1}^3 |S_t^{(i)}|\big)^\mu} \;=\; X_1 U_1^\mu + X_2 U_2^\mu + X_3 U_3^\mu \qquad \text{a.s. on } E_j.
```
The $`\mathrm{Dir}(n_1, n_2, n_3)`$ density on $`\Delta_2`$ satisfies the pointwise bound
```math
f_{(n_1, n_2, n_3)}(u_1, u_2) \;=\; (N - 1)(N - 2)\binom{N - 3}{n_1 - 1,\; n_2 - 1,\; n_3 - 1} u_1^{n_1 - 1} u_2^{n_2 - 1} u_3^{n_3 - 1} \;\le\; (N - 1)(N - 2) \;<\; N^2,
```
since the multinomial factor is a probability mass of $`\mathrm{Multinomial}(N - 3; u_1, u_2, u_3)`$ and hence at most $`1`$. For any fixed $`(x_1, x_2, x_3) \in [0, \infty)^3`$, whenever $`x_1 u_1^\mu + x_2 u_2^\mu + x_3 u_3^\mu \in A \subset [v, \infty)`$, we have $`v \le x_1 + x_2 + x_3 \le 3 M`$, so $`M = \max(x_1, x_2, x_3) \ge v/3`$ (otherwise the preimage is empty). Combining $`f_{(n_1, n_2, n_3)} < N^2`$ with (7) at $`M \ge v/3`$ and integrating over the independent $`(X_1, X_2, X_3, \mathbf 1_{E_j})`$ gives (8). $`\square`$

---

### Assembly of Theorem 1 and Proofs of Corollaries 1 & 2

**Proof of Theorem 1.**
Fix $`v \in (0, 1]`$ and a Borel set $`A \subset [v, \infty)`$ with $`\delta := |A| \in (0, 1/2]`$ (so $`\log(1/\delta) \ge \log 2`$). Set
```math
m_0 \;:=\; \lceil \log_2(1/\delta) \rceil \;\ge\; 1, \qquad
n_0 \;:=\; \max\!\left\{3,\; \left\lceil \left(\frac{\log(1/\delta)}{c_3 v^2}\right)^{\frac{1}{2\mu - 1}} \right\rceil,\; \left\lceil \left(\frac{48 \log(1/\delta)}{p_* c_2 v}\right)^{\frac{1}{\mu}} \right\rceil\right\}.
```
Observe that:
- $`2^{-m_0} \le \delta`$ and $`m_0 \le 1 + \frac{\log(1/\delta)}{\log 2} \le 3 \log(1/\delta)`$ (since $`\log(1/\delta) \ge \log 2`$);
- $`c_3 v^2 n_0^{2\mu - 1} \ge \log(1/\delta)`$ and $`\frac{p_* c_2}{48} v n_0^\mu \ge \log(1/\delta)`$;
- on $`\{C_{n_0} > c_2 v n_0^\mu\}`$, we have $`\frac{p_*}{12} C_{n_0} > 4 \cdot \frac{p_* c_2}{48} v n_0^\mu \ge 4 \log(1/\delta) > m_0`$, so Lemma 2(2) applies;
- since $`v \in (0, 1]`$, $`\log(1/\delta) \ge \log 2`$, and $`\frac{2}{2\mu - 1} > \frac{1}{2\mu - 1} > \frac{1}{\mu}`$ (as $`\mu \in (1/2, 1)`$), there is a constant $`C_1 < \infty`$ (independent of $`v, A`$) such that $`n_0 \le C_1 v^{-\frac{2}{2\mu - 1}} (\log(1/\delta))^{\frac{1}{2\mu - 1}}`$.

Decomposing $`\{X \in A\}`$ and applying Lemmas 1, 2, 3, and 4 (with a union bound over $`j = 1, \dots, m_0`$ in the last term, avoiding any disjointification of cherries):
```math
\begin{aligned}
\mathbb P(X \in A)
&\;\le\; \mathbb P\big(X \ge v,\; C_{n_0} \le c_2 v n_0^\mu\big)
\;+\; \mathbb P\big(C_{n_0} > c_2 v n_0^\mu,\; \mathcal H_{\mathrm{good}}(m_0)^c\big) \\ &\qquad +\; \mathbb P\left(\mathcal H_{\mathrm{good}}(m_0),\; \bigcap_{j=1}^{m_0} E_j^c\right)
\;+\; \sum_{j=1}^{m_0} \mathbb P\big(X \in A,\; \mathcal H_{\mathrm{good}}(m_0),\; E_j\big) \\ &\;\le\; 2\exp\!\big(-c_3 v^2 n_0^{2\mu - 1}\big)
\;+\; \exp\!\left(-\frac{p_* c_2}{48} v n_0^\mu\right)
\;+\; 2^{-m_0}
\;+\; m_0\,\frac{96\pi}{\mu(1-\mu)}\frac{\mathbb E[|T_{T^*}|^2]}{v} \delta \\ &\;\le\; 4\delta \;+\; \frac{576\pi e^{14} C_1^2}{\mu(1-\mu)}\, v^{-1 - \frac{4}{2\mu - 1}} \, \delta \, \big(\log(1/\delta)\big)^{1 + \frac{2}{2\mu - 1}} \\ &\;\le\; C \, v^{-\beta} \, |A| \, \big(\log(1/|A|)\big)^\gamma.
\end{aligned}
```
For any Lebesgue-null Borel set $`N \subset (0, \infty)`$ and any $`v \in (0, 1]`$, covering $`N \cap [v, \infty)`$ by open subsets of $`[v, \infty)`$ of positive measure $`\delta \downarrow 0`$ (padding within $`[v, \infty)`$ if needed so the cover has measure $`\delta \in (0, 1/2]`$) gives $`\mathbb P(X \in N \cap [v, \infty)) = 0`$; taking countable $`v \downarrow 0`$ with $`\mathbb P(X > 0) = 1`$ gives $`\mathbb P(X \in N) = 0`$, proving absolute continuity. $`\square`$

**Proof of Corollary 1.**
Fix $`v \in (0, 1]`$ and $`C_v := C v^{-\beta}`$. Choose a Borel nonnegative version of the density $`p_X`$ (finite a.e.). For any $`H \ge 2`$, let $`E_H := \{x \ge v : p_X(x) > H\}`$. Since $`\int_{E_H} p_X \le 1`$, we have $`|E_H| \le 1/H \le 1/2`$. If $`|E_H| > 0`$, Theorem 1 gives
```math
H |E_H| \;\le\; \mathbb P(X \in E_H) \;\le\; C_v |E_H| \big(\log(1/|E_H|)\big)^\gamma,
```
and dividing by $`|E_H|`$ yields $`\log(1/|E_H|) \ge (H/C_v)^{1/\gamma}`$, i.e., $`|E_H| \le \exp(-(H/C_v)^{1/\gamma})`$ (which holds trivially if $`|E_H| = 0`$). This proves (2).

For any $`1 \le p < \infty`$, the pointwise inequality $`y^p \le 2^{p-1} y + \int_2^\infty p H^{p-1} \mathbf 1_{\{y > H\}}\,dH`$ ($`y \ge 0`$) and Tonelli yield
```math
\int_v^\infty p_X(x)^p\,dx \;\le\; 2^{p-1} + \int_2^\infty p H^{p-1} \exp\!\left(-\left(\frac{H}{C_v}\right)^{1/\gamma}\right) dH \;<\; \infty.
```
Similarly, for any $`0 < \alpha < C_v^{-1/\gamma}`$ and any finite $`R > v`$:
- using $`\exp(\alpha y^{1/\gamma}) \le e^{\alpha 2^{1/\gamma}} + \int_2^\infty \frac{\alpha}{\gamma} H^{1/\gamma - 1} e^{\alpha H^{1/\gamma}} \mathbf 1_{\{y > H\}}\,dH`$ gives $`\int_v^R \exp(\alpha p_X(x)^{1/\gamma})\,dx < \infty`$;
- using $`y\exp(\alpha y^{1/\gamma}) \le e^{\alpha 2^{1/\gamma}} y + \int_2^\infty (1 + \frac{\alpha}{\gamma} H^{1/\gamma}) e^{\alpha H^{1/\gamma}} \mathbf 1_{\{y > H\}}\,dH`$ gives $`\int_v^\infty p_X(x)\exp(\alpha p_X(x)^{1/\gamma})\,dx < \infty`$.
(For any strictly smaller positive exponent $`\theta < 1/\gamma`$, such as $`2/33`$ or $`1/17`$, $`\alpha y^\theta`$ is eventually dominated by $`\frac{1}{2} C_v^{-1/\gamma} y^{1/\gamma}`$, so both integrals converge for every fixed $`\alpha > 0`$.) $`\square`$

**Proof of Corollary 2.**
Let $`A \subset (0, \infty)`$ be Borel with $`\delta := |A| \in (0, 1/2]`$, and write $`F(\delta) := \delta (\log(1/\delta))^\gamma > 0`$. Set
```math
v \;:=\; \min\!\Big\{1,\; F(\delta)^{\frac{1}{\beta + 1/\mu}}\Big\} \;\in\; (0, 1].
```
Because $`x \mapsto F(x)`$ is not monotone near $`x = 1/2`$, before applying Theorem 1 to $`A_v := A \cap [v, \infty)`$ (which has $`|A_v| \le \delta`$), pad $`A_v`$ to a Borel superset $`\widetilde A_v \subset [v, \infty)`$ of **exact** Lebesgue measure $`|\widetilde A_v| = \delta`$: since $`\mathrm{Leb}([v, \infty) \setminus A_v) = \infty`$, the continuous function $`t \mapsto \mathrm{Leb}([v, t] \setminus A_v)`$ attains $`\delta - |A_v| \ge 0`$ at some $`t_* \ge v`$, and $`\widetilde A_v := A_v \cup ([v, t_*] \setminus A_v) \subset [v, \infty)`$ satisfies $`A_v \subseteq \widetilde A_v`$ and $`|\widetilde A_v| = \delta \in (0, 1/2]`$.

Using Astra's left tail $`\mathbb P(X < v) \le D v^{1/\mu}`$ and Theorem 1 on $`\widetilde A_v`$:

• if $`v < 1`$, then $`v = F(\delta)^{\frac{1}{\beta + 1/\mu}}`$ satisfies $`v^{1/\mu} = v^{-\beta} F(\delta) = F(\delta)^{\xi_*} = \delta^{\xi_*} (\log(1/\delta))^{\gamma \xi_*}`$, giving

```math
\mathbb P(X \in A) \;\le\; \mathbb P(X < v) + \mathbb P(X \in \widetilde A_v) \;\le\; (D + C)\,\delta^{\xi_*} \big(\log(1/\delta)\big)^{\gamma \xi_*};
```

- if $`v = 1`$, then $`F(\delta) \ge 1`$, so $`F(\delta)^{\xi_*} \ge 1`$ and trivially $`\mathbb P(X \in A) \le 1 \le \delta^{\xi_*} (\log(1/\delta))^{\gamma \xi_*}`$.

This proves (3a). For every $`\xi \in (0, \xi_*)`$, absorbing $`(\log(1/\delta))^{\gamma \xi_*}`$ into $`\delta^{\xi_* - \xi}`$ (using absolute continuity for $`\delta = 0`$, and $`\mathbb P(X \in A) \le 1 \le 2^\xi |A|^\xi`$ for $`\delta \in (1/2, 1]`$) proves (3b). Finally, for $`H \ge 1`$, $`E_H := \{x > 0 : p_X(x) > H\}`$ has $`|E_H| \le 1/H \le 1`$, so $`H |E_H| \le C_\xi |E_H|^\xi`$ yields $`|E_H| \le (C_\xi / H)^{1/(1-\xi)}`$ (trivially if $`|E_H| = 0`$). Using $`y^p \le y + \int_1^\infty p H^{p-1} \mathbf 1_{\{y > H\}}\,dH`$ ($`y \ge 0, p \ge 1`$) and Tonelli yields $`p_X \in L^p(0, \infty)`$ for every $`1 \le p < \frac{1}{1 - \xi}`$, and hence for every $`1 \le p < \frac{1}{1 - \xi_*}`$. $`\square`$

---

## Verification and Replay

Replay command (from the repository root):
```sh
python3 -B round2/proofs/neon2/density-lp/verify_density_lp.py
```
The script checks:
1. **Exact rational exponent enclosures & balancing:** Using exact `Fraction` endpoints of Neon-2 Claim 6, verifies $`\gamma \in (16, 33/2)`$, $`\xi_* > 1/19`$ ($`>4{,}200\times`$ Astra's $`1/80{,}000`$), $`\frac{1}{1 - \xi_*} > 19/18`$, exact algebraic identity $`v^{1/\mu} = v^{-\beta} F(\delta) = F(\delta)^{\xi_*}`$, and non-monotonicity of $`F(\delta)`$ near $`\delta = 1/2`$ (justifying Corollary 2's measure-preserving padding).
2. **Exact rational planting constant & supermartingale drift:** Exact rational check of the 7-stage exponential sum $`17`$ and denominator $`48`$ (catching the rejected $`e^{-15}/24`$ defect) and exact rational drift $`(G h)(\ell) = -2(\ell+2)/(\ell(\ell+1)) < 0`$ for $`\ell = 5, \dots, 2000`$ (plus 3 planted defect controls: $`c = 1`$, uncoupled $`c = 3`$ at $`\ell = 4`$, and alarm rate $`5/(\ell+1)`$).
3. **Exact integer combinatorial packing:** Exhaustively verifies $`R = 1 \ge 1/6`$ at $`C = 1`$ and $`R \ge (C + 8)/6`$ with strict vertex-disjointness over all labelled core trees of sizes $`C = 2, \dots, 7`$ (with minimal leaf counts), exact sharpness $`6R = C + 8`$ on 6 subdivided binary trees and 6 subdivided cubic caterpillars, plus long paths and random trees up to $`C = 500`$ (catching the Type-1-only planted defect on paths).
4. **Exact 2-slot rate cancellation & floating-point numerical probes:** Exact `Fraction` check of $`\lambda_i + \kappa_i = m`$ for $`s_i \in \{1, 2\}`$; floating-point `lgamma` checks of Gautschi inequalities ($`j = 1 \dots 499`$) and tail sums ($`n_0 = 3 \dots 99`$) at both bracket endpoints, Dirichlet-3 mode peaks $`\le (N-1)(N-2)`$, and sampled-ray floating-point probes of $`-\psi_\theta''(s) \ge \mu(1 - \mu)/16`$ and slice integrals $`\le 32\pi/(\mu(1-\mu))`$ on 9 extreme/boundary triples $`(1, x_2, x_3)`$.
