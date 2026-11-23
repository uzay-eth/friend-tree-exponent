# Growth of the non-leaf population in the random friend tree

*Status: Preprint-level work, computer-assisted, not peer reviewed. How the results were produced and tested is described in the README ("Who did what") and in the review record below. The exact exponent $`\mu`$ remains open.*

---

## Executive Summary

In the random friend tree $`(T_n)_{n \ge 2}`$ (Krapivsky–Redner isotropic/complete redirection with $`r=1`$, starting from one edge $`T_2`$), let $`C_n`$ be the number of non-leaves (vertices of degree at least $`2`$). Astra proves that a single growth exponent $`\mu := \lim_{n \to \infty} \frac{\log \mathbb{E}[C_n]}{\log n}`$ exists and that $`\frac{\log C_n}{\log n} \to \mu`$ almost surely ([`proofs/existence.md`](proofs/existence.md), [`proofs/almost-sure-exponent.md`](proofs/almost-sure-exponent.md)). The existence theorem and [motif bounds](proofs/bracket-motifs.md) give $`113/200\le\mu\le1133/2000`$. The upper endpoint lies $`1/6000`$ below $`17/30`$, the fraction mentioned in [Strogatz’s 2022 post](https://x.com/stevenstrogatz/status/1572287647715934210).

In Krapivsky and Redner’s paper and in Strogatz’s thread, $`N_k`$ denotes the number of vertices of degree **exactly** $`k`$. Here $`N_k(n)`$ counts vertices of degree **at least** $`k`$, so $`N_2(n)=C_n`$. The exact-degree count is $`N_k(n)-N_{k+1}(n)\le C_n`$ for every $`k\ge2`$. Thus both finite-size mean upper bounds below apply to the expectation of each exact-degree count. The cumulative lower bounds give no matching lower bound for that difference, and its proposed pointwise scaling remains open. [Krapivsky and Redner, §5.1](https://arxiv.org/pdf/1706.01514#page=12) · [Strogatz’s thread](https://x.com/stevenstrogatz/status/1572287647715934210)

Two independently derived upper exponents fall strictly below $`\frac{17}{30}`$ for all $`n \ge 3`$: Neon's finite-$`n`$ bracket
```math
\frac{7304352308}{8763936965}\prod_{j=3}^{n-1}\!\left(1+\frac{113}{200j}\right) \le \mathbb{E}[C_n] \le \frac{14606753031}{2989767809}\prod_{j=3}^{n-1}\!\left(1+\frac{1133}{2000j}\right)
```
($`17{,}329{,}108`$ exact GMP checks per side) and Astra's independent upper bound
```math
\mathbb{E}[C_n] \le \frac{1539738114735}{3908512529}\prod_{j=3}^{n-1}\!\left(1+\frac{5657}{9983j}\right)
```
($`226{,}545{,}617`$ exact integer checks, minimum gap $`\frac{1}{299490}`$ below $`\frac{17}{30}`$; [`proofs/upper-bound-banded.md`](proofs/upper-bound-banded.md)).

Combining Fable's Theorem D(d) ([`proofs/degree-tail.md`](proofs/degree-tail.md)) with Astra's existence theorem and either interior exponent bracket proves the uniform cumulative mean degree-tail limit $`\lim_{n/k \to \infty} \frac{\log \mathbb{E}[N_k(n)]}{\log(n/k)} = \mu`$ over $`k \ge 2`$. Independently of Astra and Neon, Fable proves $`\liminf \frac{\log C_n}{\log n} \ge \frac{11}{20} = 0.55`$ a.s. ([`proofs/almost-sure-lower.md`](proofs/almost-sure-lower.md), [`proofs/cross-charging.md`](proofs/cross-charging.md)), while Carbon / Sol (using both earlier baseline tables in [`verify/baseline/`](verify/baseline/)) gives a loss-free polynomial-scale proof of $`\liminf C_n / n^{543/1000} > 0`$ a.s. and a.s. convergence of $`\Phi_n^-/n^{543/1000}`$ in $`(0,\infty]`$ ([`proofs/polynomial-scale-lower.md`](proofs/polynomial-scale-lower.md); as a standalone inequality, $`C_n/n^{0.543} \to \infty`$ a.s. already follows from Fable's $`0.55`$). All proofs were reviewed only by AI models across competing model families, every finite certificate is re-executed from scratch by [`verify/run_all.sh`](verify/run_all.sh), no human referee or formal proof assistant has audited any proof, and the exact value of $`\mu`$ remains open.

---

## 1. Model and Line-by-Line Definition Matching

### 1.1 The common discrete model

All five bodies of work (the earlier baseline and the four independent leads) study the identical discrete-time random tree process $`(T_n)_{n \ge 2}`$:

* Start from $`T_2`$, a single edge between two vertices ($`|V(T_2)| = 2`$).
* Given $`T_n`$ on $`n`$ vertices, choose $`V_n \in V(T_n)`$ uniformly at random, then choose $`W_n`$ uniformly at random among the neighbours of $`V_n`$ in $`T_n`$, and form $`T_{n+1}`$ by attaching a fresh vertex $`n+1`$ as a leaf adjacent to $`W_n`$.
* Let $`D_v(n) := \deg_{T_n}(v)`$ and $`C_n := \#\{v \in V(T_n) : D_v(n) \ge 2\}`$ (the *core* or *nucleus*). Thus $`C_2 = 0`$, $`T_3`$ has deterministic unlabelled shape $`P_3 = K_{1,2}`$ with a unique degree-$`2`$ centre ($`C_3 = 1`$), and $`C_n`$ is nondecreasing with $`C_{n+1} - C_n \in \{0,1\}`$.
* For $`k \ge 2`$, $`N_k(n) := \#\{v \in V(T_n) : D_v(n) \ge k\}`$, so $`N_2(n) = C_n`$.

For any $`n \ge 3`$ and core vertex $`u`$ with $`L_u`$ leaf neighbours and $`c_u`$ core neighbours ($`D_u = L_u + c_u \ge 2`$), the one-step attachment probabilities are identically
```math
\mathbb{P}(W_n = u \mid T_n) = \frac{1}{n}\!\left(L_u + \sum_{v \sim u,\, D_v \ge 2} \frac{1}{D_v}\right), \qquad
\mathbb{P}(W_n \text{ is a leaf adjacent to } u \mid T_n) = \frac{L_u}{n D_u},
```
because every leaf in a tree with $`n \ge 3`$ vertices has a unique core neighbour.

### 1.2 Line-by-line comparison across all five works

| Dimension | Earlier baseline (`[0.543, 0.586]`) | Astra (GPT-6.2 Astra) | Neon (Gemini 4 Neon) | Fable (Claude Fable 5.2) | Carbon (Gemini 4 Carbon / GPT-6.2 Sol) |
|---|---|---|---|---|---|
| **Starting tree & $`n`$** | $`T_2`$ one edge; $`n = \|V(T_n)\|`$; $`T_3`$ has deterministic shape $`K_{1,2}`$ ($`C_3=1`$). | $`T_2`$ one edge; $`n = \|V(T_n)\|`$; existence & a.s. proofs also use continuous-time Yule embedding $`T(t)`$ with rate-1 vertex clocks (jump chain $`\equiv (T_n)`$, independent of holding times). | $`T_2`$ one edge; $`n = \|V(T_n)\|`$; for $`n \ge 3`$, $`T_n`$ canonically rooted at the unique centre $`r`$ of $`T_3`$ (rooting adds no extra random choice). | $`T_2`$ one edge; $`n = \|V(T_n)\|`$; degree tail uses $`N_k(n) = \#\{v : D_v(n) \ge k\}`$ ($`N_2(n)=C_n`$). | $`T_2`$ one edge; $`n = \|V(T_n)\|`$; uses both earlier baseline tables (`certs9`, `certs6`). |
| **Definition of $`C_n`$** | $`\#\{v : D_v(n) \ge 2\}`$ ($`C_2=0, C_3=1`$). | Same ($`C(t)`$ in continuous time). | Same. | Same. | Same. |
| **Filtration / auxiliary state** | Lower (`certs9`): tree $`T_n`$ only. Upper (`certs6`): history-dependent birth creator degree $`P_u`$ ($`\mathcal{F}_n`$). | Existence / a.s.: continuous-time $`\mathcal{F}_t`$. Transfer / banded bounds: history-dependent birth creator degree, pooled lower-bound mark, distance-2 witness ($`\mathcal{F}_n`$). | Canonical root $`r`$ of $`T_3`$; potential depends only on rooted tree $`(T_n, r)`$ (Markov on rooted trees). | A.s. lower: tree $`T_n`$ (imports baseline `certs9`). Cross-charging: tree $`T_n`$ only (no history marks). Degree tail: tree $`T_n`$ only. | Tree $`T_n`$ for lower (`certs9`), birth creator degree $`P_u`$ ($`\mathcal{F}_n`$) for upper (`certs6`). |
| **Type of statement** | Two-sided $`\mathbb{E}[C_n]`$ bounds + a.s. $`\limsup \frac{\log C_n}{\log n}`$. | Existence of $`\mu = \lim \frac{\log \mathbb{E}[C_n]}{\log n} \in [0,1]`$; $`\frac{\log C_n}{\log n} \to \mu`$ a.s.; finite-$`n`$ $`\mathbb{E}[C_n]`$ bounds. | Theorem 1: two-sided finite-$`n`$ $`\mathbb{E}[C_n]`$ bounds. Corollary 2: $`\limsup \frac{\log C_n}{\log n} \le \frac{1133}{2000}`$ a.s. | Thm A/B: $`\liminf \frac{\log C_n}{\log n} \ge 0.543`$ a.s. Cross-charging: two-sided $`\mathbb{E}[C_n]`$ bounds + $`\liminf \frac{\log C_n}{\log n} \ge 0.55`$ a.s. Thm D: $`\mathbb{E}[N_k(n)]`$ bounds & log limit. | Asymptotic: $`\liminf C_n / n^{543/1000} > 0`$ a.s. and $`\Phi_n^-/n^{543/1000} \to L \in (0,\infty]`$ a.s. (no $`-\varepsilon`$ exponent loss). |
| **Exact bound form** | $`0.502\,\frac{\Gamma(n+0.543)}{\Gamma(n)} \le \mathbb{E}[C_n] \le 2.17\,\frac{\Gamma(n+0.586)}{\Gamma(n)}`$ (scripts also print product form). | Product $`c \prod_{j=3}^{n-1}(1+a/j)`$: transfer $`[563/1000, 568/1000]`$; banded upper $`\frac{1539738114735}{3908512529}\prod_{j=3}^{n-1}\!\left(1+\frac{5657}{9983j}\right)`$. | Exact rational product: $`\frac{7304352308}{8763936965}\prod_{j=3}^{n-1}\!\left(1+\frac{113}{200j}\right) \le \mathbb{E}[C_n] \le \frac{14606753031}{2989767809}\prod_{j=3}^{n-1}\!\left(1+\frac{1133}{2000j}\right)`$. | Cross-charging: $`0.5022\,\frac{\Gamma(n+0.55)}{\Gamma(n)} \le \mathbb{E}[C_n] \le 3.974\,\frac{\Gamma(n+0.58)}{\Gamma(n)}`$. Degree tail: explicit $`n^a`$ bounds in D(c), uniform log limit in D(d). | Asymptotic polynomial scale $`C_n / n^{543/1000}`$ (exact product $`P_a(m,n) = \prod_{j=m}^{n-1}(1+a/j)`$ inside proof). |
| **Range of $`n`$ (and $`k, m_0, m'`$)** | All $`n \ge 3`$. | Finite-$`n`$ bounds: all $`n \ge 3`$. Existence / a.s.: $`n \to \infty`$. | All $`n \ge 3`$. | Cross-charging: all $`n \ge 3`$. Degree tail: D(a) all $`n \ge 3`$, $`k \ge 2`$, $`3 \le m_0 \le n`$; D(b) all $`k \ge 2`$, $`m' \ge 3`$, $`n \ge 8m'k`$; D(c) all $`k \ge 2`$, $`n \ge 24k`$; D(d) $`n/k \to \infty`$. | Asymptotic ($`n \to \infty`$; intermediate drift bounds hold for all $`n \ge 3`$). |

### 1.3 Compatibility proofs and non-combinable boundaries

Four technical points govern how statements from different proofs can (and cannot) be combined:

1. **Exact product $`\leftrightarrow`$ $`\Gamma`$-ratio $`\leftrightarrow`$ power-law conversion.** For every $`a \in (0,1)`$ and every integer $`n \ge 3`$,
   ```math
   P_a(3,n) := \prod_{j=3}^{n-1}\!\left(1 + \frac{a}{j}\right) = \frac{2}{\Gamma(3+a)}\,\frac{\Gamma(n+a)}{\Gamma(n)},
   ```
   with the empty product $`P_a(3,3) = 1`$. Hence $`P_a(3,n) = \frac{2}{\Gamma(3+a)}\,n^a\,(1 + O(1/n))`$ as $`n \to \infty`$, and by Gautschi's inequality $`n(n+1)^{a-1} < \frac{\Gamma(n+a)}{\Gamma(n)} < n^a`$ ($`0 < a < 1`$, $`n \ge 3`$),
   ```math
   \frac{2}{\Gamma(3+a)}\!\left(\frac{3}{4}\right)^{1-a}\! n^a \;<\; P_a(3,n) \;<\; \frac{2}{\Gamma(3+a)}\,n^a \qquad (\forall n \ge 3).
   ```
   Applying this to Neon Theorem 1 (`113/200 = 0.5650`, `1133/2000 = 0.5665`) gives, for every integer $`m \ge 3`$,
   ```math
   0.41\,m^{0.5650} \;\le\; \mathbb{E}[C_m] \;\le\; 2.74\,m^{0.5665}
   ```
   (exact Gautschi lower constant $`\frac{7304352308}{8763936965}\cdot\frac{2}{\Gamma(3.565)}\cdot(3/4)^{0.435} = 0.41166674\ldots > 0.41`$; exact upper constant $`\frac{14606753031}{2989767809}\cdot\frac{2}{\Gamma(3.5665)} = 2.73020494\ldots < 2.74`$).

   Selected certified finite-$`n`$ bounds ($`n \ge 3`$) converted between product and $`\Gamma`$-ratio forms (decimals are 8-digit approximations):

   | Certificate | Side & Exponent $`a`$ | Exact / Reported Product Prefactor $`c_{\mathrm{prod}}`$ | Equivalent $`\Gamma`$-Ratio Prefactor $`c_\Gamma = \frac{2\,c_{\mathrm{prod}}}{\Gamma(3+a)}`$ | Exact Margin $`\frac{17}{30} - a`$ |
   |---|---|---|---|---|
   | Neon (`motif-bracket`) | Lower $`113/200 = 0.5650`$ | $`7304352308 / 8763936965 \approx 0.83345560`$ | $`\approx 0.46654564`$ | $`+1/600`$ |
   | Neon (`motif-bracket`) | Upper $`1133/2000 = 0.5665`$ | $`14606753031 / 2989767809 \approx 4.88558108`$ | $`\approx 2.73020494`$ | **$`+1/6000`$** |
   | Astra (`banded-upper`) | Upper $`5657/9983 \approx 0.56666333`$ | $`1539738114735 / 3908512529 \approx 393.94478`$ (indep. audit: $`393.94479`$) | $`\approx 220.10737`$ | **$`+1/299490`$** |
   | Astra (`transfer-bounds`) | Lower $`563/1000 = 0.5630`$ | $`15430450607 / 17611596558 \approx 0.87615286`$ | $`\approx 0.49155024`$ | $`+11/3000`$ |
   | Astra (`transfer-bounds`) | Upper $`568/1000 = 0.5680`$ | $`15361771189 / 640204837 \approx 23.99508767`$ | $`\approx 13.38654204`$ | $`-1/750`$ |
   | Fable (`cross-charging`) | Lower $`11/20 = 0.5500`$ | $`\approx 0.88219234`$ | $`\approx 0.50220838 \ge 0.5022`$ | $`+1/60`$ |
   | Fable (`cross-charging`) | Upper $`29/50 = 0.5800`$ | $`\approx 7.21981746`$ | $`\approx 3.97372002 \le 3.974`$ | $`-1/75`$ |
   | Baseline (`baseline`) | Lower $`543/1000 = 0.5430`$ | $`\approx 0.87551816`$ | $`\approx 0.50232622 \ge 0.502`$ | $`+71/3000`$ |
   | Baseline (`baseline`) | Upper $`293/500 = 0.5860`$ | $`\approx 3.95462373`$ | $`\approx 2.16187785 \le 2.17`$ | $`-29/1500`$ |

2. **History-dependent potentials, canonical rooting, and continuous-time embedding do not block combination.** In Neon, $`T_3`$ has deterministic shape $`K_{1,2}`$ with a unique degree-$`2`$ centre $`r`$, so rooting $`T_n`$ ($`n \ge 3`$) at $`r`$ requires no extra random choice and preserves the law of $`(T_n)_{n \ge 3}`$. In Astra's drift bounds and the baseline upper bound, the potential $`\Phi_n`$ is adapted to the full growth history $`\mathcal{F}_n`$, and $`\mathbb{E}[\Phi_{n+1} \mid \mathcal{F}_n] \le (1+\lambda/n)\Phi_n`$ implies $`\mathbb{E}[\Phi_n] \le \Phi_3 P_\lambda(3,n)`$ by the tower property; pointwise comparison $`\Phi_n \ge c_{\min} C_n`$ then yields an unconditional bound on $`\mathbb{E}[C_n]`$. In Astra's existence and sample-path proofs, assigning each vertex an independent rate-1 Poisson clock gives total jump rate $`|V(T(t))|`$ and a uniform ringing vertex independent of the holding time, so the jump chain of $`T(t)`$ has the exact law of $`(T_n)_{n \ge 2}`$, and de-Poissonization converts the continuous-time limit back to discrete $`n`$.

3. **Mean bounds do not imply existence of $`\mu`$ or almost-sure lower bounds on their own.** Neon Theorem 1 and Astra's drift certificates bound $`\mathbb{E}[C_n]`$. By Markov's inequality on dyadic times $`2^k`$, monotonicity of $`C_n`$, and Borel–Cantelli, an expectation upper bound $`\mathbb{E}[C_n] \le A n^\lambda`$ directly implies $`\limsup \frac{\log C_n}{\log n} \le \lambda`$ a.s. (Neon Corollary 2). It does **not** prove that $`\lim \frac{\log \mathbb{E}[C_n]}{\log n}`$ exists, nor does it give any almost-sure lower bound. Therefore:
   * Neon Theorem 1 bounds $`\mu`$ **only when combined with Astra's existence theorem** ([`proofs/existence.md`](proofs/existence.md)).
   * Fable Theorem D(d) requires $`\mu := \lim \frac{\log \mathbb{E}[C_m]}{\log m} \in (0,1)`$. Because Astra's existence theorem alone proves $`\mu \in [0,1]`$, Theorem D(d) combines Astra's existence theorem with **any** positive lower exponent and upper exponent strictly below $`1`$ on $`\mathbb{E}[C_m]`$ (e.g., Neon Theorem 1's $`[0.5650, 0.5665]`$, Astra's $`[0.563, 0.56666]`$, Fable's $`[0.55, 0.58]`$, or the baseline's $`[0.543, 0.586]`$).

4. **Hypotheses and scope of Carbon's polynomial-scale lower bound.** Carbon's proof ([`proofs/polynomial-scale-lower.md`](proofs/polynomial-scale-lower.md)) upgrades a drift lower certificate from $`\mathbb{E}[C_n] \ge c n^\delta`$ to $`\liminf C_n / n^\delta > 0`$ a.s. (and a.s. convergence of $`\Phi_n^-/n^\delta`$ in $`(0,\infty]`$) with **zero $`-\varepsilon`$ exponent loss**. Its hypotheses require **per-vertex** drift inequalities without inter-vertex transfers, **radius-1 locality** with degree-$`H`$ exterior saturation, an accompanying upper potential with $`f_+ \ge 1`$, and $`\lambda < 2\delta`$. None of the existing sharper lower certificates satisfy all of these hypotheses: Fable's `0.55` uses inter-vertex cross-charging transfers; Astra's `0.563` uses signed edge transfers and a distance-2 witness; Astra's `0.560` uses a distance-2 witness and pooled history marks; Neon's `0.5650` uses 2- and 3-vertex motif potentials with a non-saturating $`1/D`$ tail. Thus Carbon's theorem currently applies only to the baseline pair (`543/1000`, `293/500`). As a standalone inequality, $`\liminf C_n / n^{543/1000} > 0`$ a.s. is **already implied** both by Fable's cross-charging bound alone ($`\liminf \frac{\log C_n}{\log n} \ge 0.55`$ a.s., which itself is independent of Astra and Neon) and by Astra + Neon ($`\frac{\log C_n}{\log n} \to \mu \ge 0.5650`$ a.s.). Carbon's distinct contribution is the loss-free martingale method itself (yielding polynomial-scale growth at the certificate's own exponent $`\delta`$, plus a.s. convergence of $`\Phi_n^-/n^{543/1000}`$ in $`(0,\infty]`$).

---

## 2. The Combined Theorem (and What Is Not Proved)

### 2.1 Combined Theorem

Let $`(T_n)_{n \ge 2}`$ be the random friend tree, $`C_n = \#\{v \in V(T_n) : \deg_{T_n}(v) \ge 2\}`$, and $`N_k(n) = \#\{v \in V(T_n) : \deg_{T_n}(v) \ge k\}`$ for $`k \ge 2`$. Every part below states its exact proof dependency.

1. **Finite-$`n`$ mean bounds for the core size.**
   * *(a) Tightest two-sided exponent bracket — depends solely on Neon Theorem 1 ([`proofs/bracket-motifs.md`](proofs/bracket-motifs.md)):* For every integer $`n \ge 3`$,
     ```math
     \frac{7304352308}{8763936965}\prod_{j=3}^{n-1}\!\left(1 + \frac{113}{200j}\right)
     \;\le\; \mathbb{E}[C_n] \;\le\;
     \frac{14606753031}{2989767809}\prod_{j=3}^{n-1}\!\left(1 + \frac{1133}{2000j}\right).
     ```
     *(At small $`n`$, e.g. $`n=3`$, Fable's unrounded lower product prefactor $`0.882192\ldots`$ is larger than Neon's $`0.833455\ldots`$; Neon has the tightest two-sided exponent bracket.)* In particular, these exact product bounds imply $`\mathbb{E}[C_n] = \Omega(n^{0.5650})`$, $`\mathbb{E}[C_n] = O(n^{0.5665})`$, and the explicit power-law bounds $`0.41\,n^{0.5650} \le \mathbb{E}[C_n] \le 2.74\,n^{0.5665}`$ for all $`n \ge 3`$, so
     ```math
     \limsup_{n \to \infty}\frac{\log \mathbb{E}[C_n]}{\log n} \;\le\; \frac{1133}{2000} \;=\; \frac{17}{30} - \frac{1}{6000} \;<\; \frac{17}{30},
     \qquad\text{hence }\; \mathbb{E}[C_n] \neq n^{17/30 + o(1)}.
     ```
   * *(b) Independent upper bound below $`17/30`$ — depends solely on Astra's banded upper bound ([`proofs/upper-bound-banded.md`](proofs/upper-bound-banded.md)):* For every integer $`n \ge 3`$,
     ```math
     \mathbb{E}[C_n] \;\le\; \frac{1539738114735}{3908512529}\prod_{j=3}^{n-1}\!\left(1 + \frac{5657}{9983j}\right),
     \qquad \frac{5657}{9983} = \frac{17}{30} - \frac{1}{299490} < \frac{17}{30}.
     ```
     *(Independently, [`proofs/bounds-transfer.md`](proofs/bounds-transfer.md) proves $`\mathbb{E}[C_n] \ge \frac{15430450607}{17611596558}\prod_{j=3}^{n-1}(1+\frac{563}{1000j})`$ for all $`n \ge 3`$.)*

2. **Existence of a single exponent $`\mu`$ for the mean and almost surely.**
   * *(a) Mean exponent — depends solely on [`proofs/existence.md`](proofs/existence.md):* There exists a deterministic constant $`\mu \in [0,1]`$ such that
     ```math
     \lim_{n \to \infty}\frac{\log \mathbb{E}[C_n]}{\log n} = \mu, \qquad \text{i.e., } \mathbb{E}[C_n] = n^{\mu + o(1)}.
     ```
   * *(b) Almost-sure exponent — depends on [`proofs/existence.md`](proofs/existence.md) and [`proofs/almost-sure-exponent.md`](proofs/almost-sure-exponent.md):* Almost surely,
     ```math
     \lim_{n \to \infty}\frac{\log C_n}{\log n} = \mu.
     ```
   * *(c) Numerical location of $`\mu`$:*
     * *Tightest certified bracket — combines [`proofs/existence.md`](proofs/existence.md) with [`proofs/bracket-motifs.md`](proofs/bracket-motifs.md):*
       ```math
       \frac{113}{200} = 0.5650 \;\le\; \mu \;\le\; \frac{1133}{2000} = 0.5665 \;<\; \frac{17}{30} \qquad (\text{width } 3/2000 = 0.0015;\text{ certified minimum gap } 1/6000).
       ```
     * *Independent Astra-only bracket — combines [`proofs/existence.md`](proofs/existence.md), [`proofs/bounds-transfer.md`](proofs/bounds-transfer.md), and [`proofs/upper-bound-banded.md`](proofs/upper-bound-banded.md):*
       ```math
       \frac{563}{1000} = 0.5630 \;\le\; \mu \;\le\; \frac{5657}{9983} \approx 0.5666633 \;<\; \frac{17}{30} \qquad (\text{certified minimum gap } 1/299490).
       ```

3. **Cumulative mean degree tail.**
   * *(a) Uniform logarithmic limit — combines [`proofs/degree-tail.md`](proofs/degree-tail.md) (Theorem D(d)) with [`proofs/existence.md`](proofs/existence.md) and any positive lower exponent and upper exponent strictly below $`1`$ on $`\mathbb{E}[C_m]`$ (e.g., Neon's $`[0.5650, 0.5665]`$ or Fable's $`[0.55, 0.58]`$, supplying $`\mu \in (0,1)`$):*
     ```math
     \lim_{n/k \to \infty}\frac{\log \mathbb{E}[N_k(n)]}{\log(n/k)} = \mu \qquad \text{uniformly over all integers } k \ge 2.
     ```
   * *(b) Explicit finite-$`(n,k)`$ bounds — combines Theorem D(c) ([`proofs/degree-tail.md`](proofs/degree-tail.md)) with Neon Theorem 1 ([`proofs/bracket-motifs.md`](proofs/bracket-motifs.md)) via Gautschi's inequality ($`0.41\,m^{0.5650} \le \mathbb{E}[C_m] \le 2.74\,m^{0.5665}`$ for all $`m \ge 3`$):* For every integer $`k \ge 2`$ and every integer $`n \ge 24k`$,
     ```math
     \frac{0.41}{44160}\!\left(\frac{n}{12k}\right)^{0.5650}
     \;\le\; \mathbb{E}[N_k(n)] \;\le\;
     22\!\left(\frac{n}{k}\right)^{0.5665}
     ```
     *(Theorem D(c) yields upper constant $`A\bigl(2 + \frac{1}{1-\lambda}\bigr)3^\lambda = 21.988467\ldots < 22`$ with $`A = 2.74`$, and $`21.909862\ldots`$ using the sharper $`\Gamma`$-ratio constant $`A = 2.7302049\ldots`$, at $`\lambda = 0.5665`$.)*

4. **Almost-sure bounds not relying on Astra's existence/sample-path proofs.**
   * *(a) Upper almost-sure bound — Neon Corollary 2 ([`proofs/bracket-motifs.md`](proofs/bracket-motifs.md)):* $`\limsup_{n\to\infty} \frac{\log C_n}{\log n} \le \frac{1133}{2000}`$ a.s. (similarly $`\le \frac{5657}{9983}`$ from Astra's banded bound by the same dyadic Markov + Borel–Cantelli argument).
   * *(b) Logarithmic lower bound — Fable ([`proofs/almost-sure-lower.md`](proofs/almost-sure-lower.md), [`proofs/cross-charging.md`](proofs/cross-charging.md)):* $`\liminf_{n\to\infty} \frac{\log C_n}{\log n} \ge \frac{11}{20} = 0.55`$ a.s. (and $`\ge 0.543`$ using the baseline `certs9` table).
   * *(c) Loss-free polynomial-scale lower bound from the baseline tables — Carbon / Sol ([`proofs/polynomial-scale-lower.md`](proofs/polynomial-scale-lower.md)), using both tables in [`verify/baseline/`](verify/baseline/):*
     ```math
     \liminf_{n \to \infty}\frac{C_n}{n^{543/1000}} > 0 \qquad \text{almost surely,}
     ```
     and the normalized baseline lower potential $`\Phi_n^- / n^{543/1000}`$ converges almost surely in $`(0,\infty]`$. *(As a standalone inequality on $`C_n`$, $`\liminf C_n / n^{543/1000} > 0`$ a.s. already follows from 4(b), since $`0.55 > 0.543`$; Carbon's distinct content is eliminating the $`-\varepsilon`$ exponent loss at the certificate's own exponent and proving a.s. convergence of $`\Phi_n^-/n^{543/1000}`$ in $`(0,\infty]`$.)*

### 2.2 What is NOT proved

To prevent any overstatement, the following remain **unproved**:

* **Exact value or closed form of $`\mu`$.** Nothing proves whether $`\mu`$ is rational, algebraic, transcendental, or expressible by any known formula.
* **Simulation estimates are not proofs.** High-precision continuous-time Monte Carlo runs give an effective-slope estimate $`\approx 0.565723 \pm 0.000003`$ (sampling standard error; 912 runs to $`n \approx 4.3 \times 10^{12}`$ and 91 runs to $`n \approx 8.1 \times 10^{13}`$, $`2.22 \times 10^{10}`$ core firings; consistent with independent lazy simulation $`0.565724 \pm 0.000018`$ at fixed time $`t=27`$, $`\mathbb{E}[N(t)] \approx 1.6 \times 10^{12}`$). These heuristics sit inside $`[0.5650, 0.5665]`$ and $`\sim 300`$ sampling standard errors below $`17/30`$, but do not guarantee asymptotic convergence or accuracy for the true limit $`\mu`$.
* **Regular variation / polynomial prefactor at $`\mu`$.** Neither $`\lim \mathbb{E}[C_n]/n^\mu \in (0,\infty)`$ nor $`\liminf C_n / n^\mu > 0`$ nor convergence of $`C_n / \mathbb{E}[C_n]`$ is proved. Astra's theorems permit arbitrary subpolynomial factors $`n^{o(1)}`$ (such as powers of $`\log n`$).
* **Pointwise degree distribution and uniform growing-$`k`$ almost-sure degree tail.** Fable Theorem D proves a *cumulative* *mean* statement for $`\mathbb{E}[N_k(n)]`$. For each **fixed** $`k \ge 2`$, Theorem 4.8 of Addario-Berry et al. ([arXiv:2403.20185](https://arxiv.org/abs/2403.20185)) already yields $`N_k(n) \asymp C_n`$ a.s. (and $`N_2(n) = C_n`$ is immediate), so combining it with Astra's sample-path theorem gives $`\lim_{n\to\infty}\frac{\log N_k(n)}{\log(n/k)} = \mu`$ a.s. simultaneously for all fixed $`k \ge 2`$ via a countable intersection. What remains **unproved** is: (i) an almost-sure logarithmic limit **uniform over growing $`k`$** (with $`n/k \to \infty`$), (ii) a two-sided constant-factor comparison $`\mathbb{E}[N_k(n)] \asymp \mathbb{E}[C_{n/k}]`$, and (iii) the Krapivsky–Redner pointwise scaling $`\mathbb{E}[\#\{v : D_v(n) = k\}] \sim c\,n^\mu k^{-(1+\mu)}`$ or convergence of normalized degree counts.
* **Excluded constructions.** An alternate tagged-neighbour construction by Astra (`tag.cpp`), whose strict $`< 17/30`$ upper bound remained unverified (`GAP` verdict in [`reviews/review-excluded-tagged-fable.md`](reviews/review-excluded-tagged-fable.md)), is excluded and not used in any theorem above.

---

## 3. Short Account of Each Method

**Astra — Banded drift certificate with two cancellations** ([`proofs/existence.md`](proofs/existence.md), [`proofs/almost-sure-exponent.md`](proofs/almost-sure-exponent.md), [`proofs/bounds-hub.md`](proofs/bounds-hub.md), [`proofs/bounds-transfer.md`](proofs/bounds-transfer.md), [`proofs/upper-bound-banded.md`](proofs/upper-bound-banded.md)). Astra constructs a supermartingale-type vertex potential $`\Phi_n = \sum_{u \in \mathrm{core}(T_n)} f(x_u)`$ where each core vertex $`u`$'s label $`x_u`$ stores its degree $`D_u`$, core degree $`c_u`$, observed degrees and possibly capped core degrees of tracked neighbours, pooled lower-bound marks for untracked neighbours, and one distinguished distance-2 witness, with degree cutoffs $`h_c(D)`$ decreasing across 7 degree bands up to $`H = 2400`$ and an explicit affine tail $`f = A + B\,q`$ ($`q = L/D`$) for $`D \ge 2400`$. To tighten the per-vertex drift bound without expanding the state space combinatorially, Astra introduces two exact zero-sum cancellations in the summed drift: (1) *signed edge transfers* $`g(D,c; d,k)`$ (`63,309` entries) that add and subtract the same gain-rate term across observable core-core edges, and (2) nonnegative *one-way arrow transfers* $`a(D_u, c_u)/D_w`$ (`18` entries, bounded by $`a_*`$) subtracted from $`u`$'s row and added to core neighbour $`w`$'s row for every ordered core edge $`w \to u`$ with $`D_w \ge T = 8`$, which cancel identically on every realized tree independently of local states. Separately, the integer trial rewards used to upper-bound newborn contributions are verified to exceed the actual $`c=1`$ newborn weights in the completed table (minimum reserve `4,340,525` scale units). Replaying the compressed recipe in 128-bit integer arithmetic verifies all $`226{,}545{,}617`$ finite states plus the analytic $`D \ge 2400`$ tail at $`\lambda = 5657/9983 = 17/30 - 1/299490`$ (minimum integer drift slack `3`).

**Neon — Directed 3-vertex and cherry motif potential with analytic $`1/D`$ tail** ([`proofs/bracket-motifs.md`](proofs/bracket-motifs.md), [`verify/motif-bracket/`](verify/motif-bracket/)). Rather than tracking history marks or a high-dimensional neighbourhood vector at each vertex, Neon roots $`T_n`$ ($`n \ge 3`$) canonically at the unique centre $`r`$ of $`T_3`$ and writes $`\Phi(T_n)`$ as a sum of four local motif tables: a single-vertex potential $`f(L,k,h)`$ on core vertices ($`L \le 800`$, $`k \le 150`$, plus overflow states $`\mathrm{VSTAR}_0, \mathrm{VSTAR}_1`$), a directed parent–child core edge potential $`g`$, a 3-vertex directed vertical path potential $`h_V`$ on core triples $`p \to u \to w`$ (`84` active middle states), and a sibling cherry potential $`h_C`$ on core pairs sharing a parent (`75` active apex states). For $`g, h_V, h_C`$, endpoints of degree $`D < D_0 = 24`$ use exact inner states while tail endpoints ($`D \ge D_0`$, with tail child-count parameter $`K_T = 12`$) use exact multilinear $`D_0/D`$ interpolation between $`\mathrm{INF}`$ and $`T_1(\kappa, h)`$, eliminating large-degree truncation slack; unmodeled 3-vertex credit products are decoupled via piecewise-linear Fenchel–Young centring interpolants on `9` rational knots. Because a single leaf attachment alters motifs within a bounded graph radius (whose count can grow with degree), Neon proves an exact global tree-level drift decomposition (Lemma 4) bounding $`n\,\mathbb{E}[\Phi(T_{n+1}) - \Phi(T_n) \mid T_n] - \lambda\,\Phi(T_n)`$ by `17,329,108` explicit rational inequalities across all $`f, g, h_V, h_C`$ states per side. Certified in exact GMP integer arithmetic, the same architecture produces both the sharpest upper bound $`\lambda = 1133/2000 = 17/30 - 1/6000`$ and the sharpest lower bound $`\delta = 113/200 = 0.5650`$, with compact product constants (`0.8335` and `4.8856`) and ~20-second end-to-end verification.

**Fable — Cross-charging and the cumulative degree tail** ([`proofs/cross-charging.md`](proofs/cross-charging.md), [`proofs/degree-tail.md`](proofs/degree-tail.md)). To keep certificates tiny while breaking past the single-vertex bottleneck, Fable pairs a modest single-vertex potential $`f(L, c, \text{small-degree multiset})`$ ($`D_{\mathrm{cap}} = 4`$, $`c_1 = 3`$, $`L_{\max} = 60`$, $`c_{\max} = 12`$; only $`5{,}666`$ table entries) with two zero-sum cross-charging redistribution terms across core-core edges: a TOP-giver edge charge $`\kappa(\mathrm{cls}(u))/D_w`$ and a non-TOP middle-vertex charge $`\nu(\mathrm{cls}(w), \mathrm{cls}(u))/D_z`$ that move drift credit from high-degree hubs to low-degree neighbours without adding any history-dependent labels. This verifies $`\mathbb{E}[C_n] \in [0.5022\,\frac{\Gamma(n+0.55)}{\Gamma(n)},\, 3.974\,\frac{\Gamma(n+0.58)}{\Gamma(n)}]`$ in $`5{,}851`$ / $`5{,}853`$ exact rational inequalities (under $`1`$ second). Separately, Fable proves Theorem D on the cumulative mean degree tail $`\mathbb{E}[N_k(n)]`$ directly from the discrete dynamics: part (a) combines a conditional first-moment recursion for vertex degree with Markov's inequality (followed by summation by parts in part (c)), while part (b) activates eligible vertices already present in the core at time $`m'`$ and applies conditional first and second moments with the Paley–Zygmund inequality.

**Almost-sure arguments — Astra, Fable, and Carbon** ([`proofs/existence.md`](proofs/existence.md), [`proofs/almost-sure-exponent.md`](proofs/almost-sure-exponent.md), [`proofs/almost-sure-lower.md`](proofs/almost-sure-lower.md), [`proofs/polynomial-scale-lower.md`](proofs/polynomial-scale-lower.md)). Three distinct mechanisms bridge mean bounds to almost-sure statements, with different trade-offs between conclusion strength and hypotheses: (1) *Astra* embeds $`(T_n)`$ in continuous time via independent rate-1 clocks (Yule total size $`N(t)`$); existence plants rooted $`k`$-vertex seeds at at least $`C/4`$ eligible sites (needing no independence between plantings, using marginal first moments) and establishes $`\mu = \lim t^{-1}\log \mathbb{E}[C(t)]`$ by iterating an almost-supermultiplicative inequality with $`\eta \to 0`$ plus de-Poissonization, while the sample-path proof runs conditionally independent lower trials inside logarithmic blocks and applies a martingale-difference law across blocks to prove $`\frac{\log C_n}{\log n} \to \mu`$ a.s. (2) *Fable (Theorem B)* is a self-contained martingale transfer theorem: for any bounded, locally supported nonnegative potential $`\Phi_n = \sum_{u \in \mathrm{core}(T_n)} f(x_u)`$ satisfying hypotheses (H1)–(H5) (including a global drift lower bound $`\mathbb{E}[\Delta\Phi_n \mid \mathcal{F}_n] \ge \frac{\delta}{n}\Phi_n`$ and positive newborn weight $`\eta > 0`$), the stopped reciprocal supermartingale $`1/(\Phi_{n \wedge \tau}\, p_{n \wedge \tau})`$ with $`p_m = \prod_{j=\sigma}^{m-1}(1 - (\delta-\varepsilon)/j)`$—started at the a.s.-finite crossing time $`\sigma`$ guaranteed by recurring conversions—proves $`\liminf \frac{\log C_n}{\log n} \ge \delta`$ a.s., giving both $`0.543`$ (baseline) and $`0.55`$ (cross-charging). (3) *Carbon / Sol* eliminates the $`-\varepsilon`$ exponent loss at $`\delta = 543/1000`$: isolating a replayed finite pattern inside a branch behind a degree-saturated parent ($`D_p \ge H`$) where exterior moves cannot perturb radius-1 states, it bounds the conditional second moment $`\mathbb{E}[(\Delta X_j)^2 \mid \mathcal{F}_j] \le (H+3)B_F^2 N_j / j`$ via $`N_j \le Y_j`$ and the baseline upper certificate ($`\lambda = 0.586 < 2\delta = 1.086`$), proving $`L^2`$-boundedness of the **centered martingale part $`Z_n`$** in the Doob decomposition $`M_n = X_m + A_n + Z_n`$ (noting explicitly that the full normalized sum $`M_n`$ can have unbounded second moments and may converge to $`+\infty`$ in $`(0,\infty]`$); conditional Cauchy–Schwarz gives a uniform positive lower bound on $`\mathbb{P}(\liminf \Phi_n^-/n^\delta > 0 \mid \mathcal{F}_m)`$ at every starting size $`m`$, which lifts to probability $`1`$ by approximation in $`\sigma(\bigcup_m \mathcal{F}_m)`$ (trials need not be independent and can share history).

---

## 4. Clean-Clone Verification

All seven finite verification suites are packaged in [`verify/`](verify/) with pinned Python dependencies ([`verify/requirements.txt`](verify/requirements.txt): `numpy==2.1.3`, `mpmath==1.3.0`, `scipy==1.14.1`) and a single entry point [`verify/run_all.sh`](verify/run_all.sh).

### 4.1 Environment

* **Hardware:** AMD Ryzen 7 5800X (8 cores / 16 threads), 62 GiB RAM.
* **OS & Toolchain:** Ubuntu Linux (`6.8.0-149-generic` x86_64), `gcc`/`g++` `13.3.0` (with OpenMP and `-fsanitize=undefined`), `Python 3.12.3`, `libgmp-dev 2:6.3.0+dfsg-2ubuntu6.1`, `gzip 1.12`.

### 4.2 Packaging notes

Three packaging notes apply to this self-contained repository (no proof text or certificate was altered):
1. During initial clean-clone testing of the unmerged lead work, the baseline scripts (`certs9.py`, `find_upper_lin.py`) and Carbon's `--report-tables` option failed without `numpy` and `mpmath` (and Fable's `xcharge.py` requires `scipy`). All three are pinned in [`verify/requirements.txt`](verify/requirements.txt) and installed automatically into `verify/.venv` by [`verify/run_all.sh`](verify/run_all.sh).
2. Relative path lookups in three helper scripts (`verify/motif-bracket/verify_python_cross.py`, `verify/cross-charging/as_lower_check.py`, and `verify/polynomial-lower/verify_as_lower.py`) were adapted to the `verify/` directory layout so every check runs from a fresh clone with zero manual path setup.
3. On Astra's banded bound, the primary `check.cpp` (`prefactor = 1539738114735/3908512529 ≈ 393.944782`) and the independent `verify_banded.cpp` (`prefactor = 1539738113358/3908512428 ≈ 393.944791`) differ by $`10^{-5}`$ in the product prefactor due to conservative integer rounding; both certify the identical exponent `5657/9983`.

*(Note on historical references inside `proofs/`: references inside `proofs/polynomial-scale-lower.md`, `proofs/almost-sure-lower.md`, and `proofs/cross-charging.md` to "`REPORT.md` §§3.3–3.5" or "the earlier report" refer to the earlier `[0.543, 0.586]` baseline write-up rather than this repository's top-level `REPORT.md`. That earlier baseline write-up is not packaged separately here; its exact verification scripts live in [`verify/baseline/`](verify/baseline/) and its two cold reviews in [`reviews/review-baseline-gemini.md`](reviews/review-baseline-gemini.md) and [`reviews/review-baseline-glm.md`](reviews/review-baseline-glm.md).)*

### 4.3 Verification suites executed by `verify/run_all.sh`

| Suite | Verifier(s) | Exact checks & negative controls |
|---|---|---|
| `motif-bracket` (Neon) | `cert_build.c`, `cert_verify_gmp.c`, `verify_python_cross.py` | `17,329,108` GMP checks/side (`1133/2000`, `113/200`), 37 C trees/side, 47 Python trees, `n=3..21` exact `E[C_n]` cross-check, 6/6 negative controls rejected. |
| `banded-upper` (Astra) | `check.cpp` (UBSan) | `226,545,617` rows, `63,309` signed transfers, `18` arrows, $`D \ge 2400`$ tail, min drift slack `3`, `PASS upper 5657/9983`. |
| `banded-upper-audit` (indep.) | `verify_banded.cpp` (UBSan), `test_banded_audit.py` | Smoke `3/5`, `436` exhaustive histories, `50` adversarial trees, 14/14 negative controls, full `226,545,617`-row production `PASS upper 5657/9983`. |
| `transfer-bounds` (Astra) | `transfer.cpp`, `hub.cpp` (UBSan) | `9,639,742` finite states + $`D \ge 208`$ tail (`563/1000`, `568/1000`); hub `[56/100, 57/100]`. |
| `almost-sure-lower` (Fable) | `as_lower_check.py`, `certs9.py` | `1,205,767` exact rational checks at `543/1000`, hypotheses (H1)–(H5), negative control (`11/20`) rejected. |
| `cross-charging` (Fable) | `xcharge.py`, `xdump.py`, `xverify.py` | 200-tree drift identity, `5,853` upper (`29/50`) and `5,851` lower (`11/20`) exact checks, 2/2 negative controls (`579/1000`, `552/1000`) rejected. |
| `polynomial-lower` (Carbon/Sol) | `find_upper_lin.py`, `verify_as_lower.py`, `check_as_lower.py` | `319,086` exact baseline upper checks (`293/500`) + negative control (`29/50`), 102 trees, 84 saturated branches, 80 grafts, independent audit with 14/14 negative controls caught. |

---

## 5. Review and Attack Record

All cold reviews and attack reports preserved from the project live in [`reviews/`](reviews/) (indexed in [`reviews/README.md`](reviews/README.md)). Every proof included in Section 2 passed cold review by at least one model from a different model family than its author, and the combined report received two cold cross-family consistency checks ([`reviews/review-combined-claude.md`](reviews/review-combined-claude.md), [`reviews/review-combined-gpt.md`](reviews/review-combined-gpt.md), both `CONSISTENT` / `APPROVE`).

Two non-trivial `GAP` verdicts occurred during review:
* **Degree tail (`proofs/degree-tail.md`):** the first GPT-6.2 Sol cold review returned `GAP` because the summation-by-parts equality in part (c) failed at the admissible endpoint $`k = 2`$ ($`m_0 = n`$), alongside two overstated prose sentences. The proof was repaired and re-checked `CORRECT` ([`reviews/review-degree-tail-gpt.md`](reviews/review-degree-tail-gpt.md); second cold review [`reviews/review-degree-tail-gemini.md`](reviews/review-degree-tail-gemini.md)).
* **Excluded tagged-neighbour construction:** Claude Fable 5.2 returned `GAP` on an alternate upper-bound construction (`tag.cpp`) whose strict $`< 17/30`$ bound remained unverified ([`reviews/review-excluded-tagged-fable.md`](reviews/review-excluded-tagged-fable.md)). It is excluded from all claims.

### Stated limits of the audits

* **No human refereeing or formal proof-assistant verification.** All proofs and C/C++/Python verifiers have been audited by AI models only.
* **Tree-level vs. row-level coverage at high degree:** in the Opus 5.5 attack on Astra's banded bound ([`reviews/attack-upper-bound-banded-claude.md`](reviews/attack-upper-bound-banded-claude.md)), direct tree-level term-by-term verification mainly covered vertices with $`D < 1000`$ (testing $`D \ge 1000`$ via ~36 constructed hub histories), while the full `226,545,617`-row algebraic verification and analytic $`D \ge 2400`$ tail were verified end-to-end by both `check.cpp` and `verify_banded.cpp`. In the Opus 5.5 attack on Neon ([`reviews/attack-bracket-motifs-claude.md`](reviews/attack-bracket-motifs-claude.md)), direct tree-level checks covered all $`32{,}972`$ rooted trees up to $`n=14`$ plus thousands of targeted trees, with partial coverage of $`n > 1200`$ hub trees; the algebraic decomposition (Lemma 4) was audited analytically and all `17,329,108` inequalities per side were re-verified by four independent verifiers.

---

## 6. Prior Work and Literature Check (2026-11-22)

* **Saramäki & Kaski (2004)** (*Scale-free networks generated by random walkers*, *Physica A* 341, 80–86) introduced random-walk attachment models related to redirection.
* **Cannings & Jordan (2013)** (*Random walk attachment graphs*, [arXiv:1303.1052](https://arxiv.org/abs/1303.1052)) analyzed random-walk attachment of length $`\ell`$; for $`\ell = 1`$ (which coincides with the random friend tree) their results show $`C_n = o(n)`$ a.s. and leave the nucleus exponent open.
* **Krapivsky & Redner (2017)** (*Emergent Network Modularity*, [arXiv:1706.01514](https://arxiv.org/abs/1706.01514), *J. Stat. Mech.* 2017, 073405) introduced *complete (isotropic) redirection* ($`r = 1`$), observed numerically that the number of non-leaves scales as $`N^\mu`$ with $`\mu \approx 0.566`$ and degree counts $`N_k \sim N^\mu k^{-(1+\mu)}`$, and highlighted the puzzle of determining $`\mu`$. (See also **Krapivsky & Redner (2023 / 2024)**, *The Magic of Networks Grown by Redirection*, [arXiv:2305.10628](https://arxiv.org/abs/2305.10628), *Indian J. Phys.* 98, 3777.)
* **Steven Strogatz (2022)** publicized the question of whether the exponent equals $`17/30 = 0.5666\overline{6}`$ ([2022-09-20 post](https://x.com/stevenstrogatz/status/1572287647715934210)).
* **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou & Lugosi (2024 / 2026)** (*Random friend trees*, [arXiv:2403.20185](https://arxiv.org/abs/2403.20185); *Random Structures & Algorithms* 69(1), e70082, 2026, [DOI](https://doi.org/10.1002/rsa.70082)) coined the name *random friend tree* and proved the first non-trivial rigorous bounds on the nucleus size: Propositions 6.20 and 6.22 show $`n^{-\delta} C_n \to \infty`$ a.s. for every $`\delta < \frac{\sqrt{13}-3}{6} \approx 0.100925`$ and $`n^{-\lambda} C_n \to 0`$ a.s. for every $`\lambda > \sqrt{2} - \frac{1}{2} \approx 0.914214`$, giving
  ```math
  \frac{\sqrt{13}-3}{6} \;\le\; \liminf_{n\to\infty}\frac{\log C_n}{\log n} \;\le\; \limsup_{n\to\infty}\frac{\log C_n}{\log n} \;\le\; \sqrt{2} - \frac{1}{2} \qquad \text{almost surely}
  ```
  *(the arXiv abstract displays $`n^{9/10}`$ as an upper ceiling, whereas the written proof of Proposition 6.22 yields $`1 - \gamma/2 = \sqrt{2}-1/2 \approx 0.914214`$)*. Their Theorem 4.8 also proves $`c_k < \liminf_{n\to\infty} N_{k+1}(n)/N_k(n)`$ a.s. (with $`N_{k+1}(n)/N_k(n) \le 1`$ by nesting) for every fixed $`k \ge 2`$. Existence and exact value of $`\mu`$ were left open.
* **Hartle, Krapivsky, Redner & Zhang (2026)** (*Anomalous scaling in redirection networks*, [arXiv:2604.01540](https://arxiv.org/abs/2604.01540), *Phys. Rev. E* 114, 034314) cite the earlier $`\mu \approx 0.566`$ estimate, report $`\mu \approx 0.569`$ from their own finite-size isotropic-redirection fit over $`N = 10^2\text{--}10^6`$ (Fig. 6), and solve two *modified* local redirection models (DAN and PAN) analytically, noting that isotropic redirection itself still lacks a closed recurrence or exact exponent.

**Dated literature check (2026-11-22):** live arXiv API queries across `math.PR`, `math.CO`, `cond-mat.stat-mech`, and `physics.soc-ph` and web searches found no published paper or preprint excluding $`\mu = 17/30`$, proving existence of $`\mu`$, or determining the exact value of $`\mu`$.

---

## 7. Open Problems

1. **Tightening the certified bracket beyond $`[0.5650, 0.5665]`$.** Neon's motif-plus-analytic-tail architecture verifies in ~20 seconds at `17.3M` checks; extending its active motif domains or combining its analytic $`1/D`$ tail with signed transfers could plausibly narrow the certified interval substantially around the simulated $`\approx 0.565723`$.
2. **Sharp asymptotics ($`C_n / n^\mu`$) and polynomial-scale convergence.** Proving that $`\mathbb{E}[C_n] \asymp n^\mu`$ (or $`\mathbb{E}[C_n]/n^\mu \to c \in (0,\infty)`$) and that $`C_n / \mathbb{E}[C_n]`$ converges to a non-degenerate random variable $`W > 0`$ remains open.
3. **Exact value or operator characterization of $`\mu`$.** Whether $`\mu`$ can be characterized as the leading eigenvalue / Malthusian parameter of a tractable infinite-dimensional linear operator on a local weak limit, or has no closed form at all, is completely open.
4. **Uniform growing-$`k`$ almost-sure degree tail and pointwise degree distribution.** While fixed-$`k`$ almost-sure convergence $`\log N_k(n)/\log(n/k) \to \mu`$ follows from combining Astra's sample-path theorem with Addario-Berry et al. Theorem 4.8, proving an almost-sure limit uniform over growing $`k`$ (with $`n/k \to \infty`$), proving $`\mathbb{E}[N_k(n)] \asymp \mathbb{E}[C_{n/k}]`$, and establishing the pointwise Krapivsky–Redner law $`\mathbb{E}[\#\{v : D_v(n) = k\}] \sim c\,n^\mu k^{-(1+\mu)}`$ all remain open.
5. **Human refereeing and formalization.** A human probabilist has not yet audited the proofs, and none of the finite-arithmetic checkers has been formalized in Lean, Coq, or Isabelle.
