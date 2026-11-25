> **Public checkpoint export note:** Edited export of `degrees-neon/CONTINUATION_DISCRETE_MEAN.md` (frozen original SHA-256 `d3ba5f9147e37e133810ee001fb4945c335068e1469f76e714b9903eddb0e456`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.** **IMPORTANT:** Frozen line 92 of the original proof contains a false Beta-ratio bound (caught by the explicit counterexample recorded by the audit). Read together with [`ERRATUM_DISCRETE_MEAN.md`](ERRATUM_DISCRETE_MEAN.md) and [`REVIEW.md`](REVIEW.md), which repair the proof and verify Theorem A.

# Continuation: Pure Tier A Pointwise Discrete-Time Mean via an $`O(\log^* k)`$ Iterated-Log Tower and Effective Depoissonization

**Author:** Gemini 4 Neon (`degrees-neon`).
**Status (historical pre-audit line):** Post-snapshot continuation (separate from the preserved reciprocal-reviewed snapshot `degrees-neon/PROOF.md`, `CHECKPOINT.md`, `REVIEW.md`). Written as self-checked by Gemini 4 Neon before outside-family review; *subsequently cross-family audited in [`REVIEW.md`](REVIEW.md) and repaired in [`ERRATUM_DISCRETE_MEAN.md`](ERRATUM_DISCRETE_MEAN.md); still awaiting human refereeing.*
**Pinned baseline:** `uzay-eth/friend-tree-exponent` at `c250a1f2fdd53ced91b7b485e94ac8a16106b9c7`.
**New finite identity checker:** `verify_tier_a_discrete_mean.py` (internal finite checker, not exported; see [`audit_exact.py`](audit_exact.py) for the exported reviewer exact-rational controls).

---

## 0. Exact Claim, Scope, and Honest Lineage Note

### Main Results Proved Here (100% Pure Tier A Inputs Only)
Let $`(T_n)_{n \ge 2}`$ be the discrete random friend tree started from $`T_2`$ (or from any reachable finite seed $`\tau`$ of size $`m_0 = |\tau| \ge 2`$), and let $`a_n(k) := \mathbb E[N_{=k}(T_n)]`$ (resp. $`a_n^{(\tau)}(k) := \mathbb E_\tau[N_{=k}(T_n)]`$).

- **Theorem A (Pure Tier A Pointwise Discrete-Time Exact-Degree Mean — Full Intermediate Regime):**
  Using **only Tier A inputs** from the pinned manuscript (Theorems 1, 4, 5: $`\mu \in (1/2, 1)`$ with $`e^{1/3} - 1 < \mu`$, $`C_n/n^\mu \to X`$ a.s. and in $`L^2`$ with $`\mathbb E[X] = c_2 > 0`$, fixed-degree law $`N_{=d}(T_n)/C_n \to p_d`$ a.s. for each fixed $`d \ge 2`$ with $`p_d \sim \mu A d^{-(1+\mu)}`$ as $`d \to \infty`$, plus Neon's Tier A Lemmas 1 and 2 from `degrees-neon/PROOF.md` (pre-continuation snapshot, not exported; Tier A Lemmas 1–2 are reconstructed in `REVIEW.md` §2)), uniformly as $`k \to \infty`$ with $`k/n \to 0`$:

```math
\lim_{\substack{k \to \infty \\ k/n \to 0}} \frac{\mathbb E[N_{=k}(T_n)]}{\mu A c_2 \, n^\mu k^{-(1+\mu)}} = 1,
```

and for every reachable finite seed $`\tau`$ with $`C(\tau) \ge 1`$:

```math
\lim_{\substack{k \to \infty \\ k/n \to 0}} \frac{\mathbb E_\tau[N_{=k}(T_n)]}{\mu A \mathbb E_\tau[X] \, n^\mu k^{-(1+\mu)}} = 1 \qquad \left(\mathbb E_\tau[X] = \frac{\Gamma(|\tau|)}{\Gamma(|\tau|+\mu)} d_2 \Psi(\tau)\right).
```

**Zero Tier B inputs are used:** neither `hub-spectrum.md` §1 eq. (2) (small-power exact-degree law) nor `hub-spectrum.md` Claim 11a (intermediate cumulative law) is assumed.

- **Corollary B (Pure Tier A Upgrades of Cumulative Mean, Continuous-Time Mean, and Effective Depoissonization):**
  From Theorem A + Neon's Tier A Lemmas 1–2 alone (upgrading all mean statements of our earlier checkpoint from Tier A + Claim 11a to **pure Tier A**):
  1. *(Pure Tier A full-intermediate cumulative mean)* Uniformly as $`k \to \infty`$ with $`k/n \to 0`$:

```math
\mathbb E[N_{\ge k}(T_n)] \sim A c_2 \, n^\mu k^{-\mu}, \qquad \mathbb E_\tau[N_{\ge k}(T_n)] \sim A \mathbb E_\tau[X] \, n^\mu k^{-\mu}.
```

2. *(Pure Tier A continuous-time mean)* Uniformly as $`k \to \infty`$ with $`k e^{-t} \to 0`$:

```math
\mathbb E_{T_2}[N_{=k}(T(t))] \sim \mu A \Gamma(2+\mu) c_2 \, e^{\mu t} k^{-(1+\mu)}, \qquad \mathbb E_\tau[N_{=k}(T(t))] \sim \mu A d_2 \Psi(\tau) \, e^{\mu t} k^{-(1+\mu)}.
```

3. *(Exact Yule-to-Poisson conditioning & uniform high-degree depoissonization)* Conditional on the terminal Yule size $`Y_2 = \lim_{s\to\infty} e^{-s}|T(s)| \sim \mathrm{Gamma}(2, 1)`$, the size satisfies $`|T(t)| - 2 \mid \{Y_2 = y\} \sim \mathrm{Poisson}((e^t - 1)y)`$ identically, so $`\mathbb E_{T_2}[N_{=k}(T(t))] = \int_0^\infty y e^{-y} \Pi_k((e^t-1)y)\,dy`$ where $`\Pi_k(\lambda) := \sum_{m \ge 2} a_m(k) e^{-\lambda}\lambda^{m-2}/(m-2)!`$ is the exact Poissonization. Uniformly over all high degrees $`k \to \infty, k/n \to 0`$:

```math
\frac{a_n(k)}{\mathbb E_{T_2}[N_{=k}(T(\log n))]} \longrightarrow \frac{1}{\Gamma(2+\mu)}, \qquad \frac{a_n(k)}{\Pi_k(n)} \longrightarrow 1.
```

- **Proposition C (Exact Depoissonization Obstruction Without Transport):**
  Despite the exact one-step discrete regularity lower bound $`a_{n+1}(k) \ge (1 - k/n) a_n(k)`$ (for all $`n \ge 2, k \ge 3`$), **no black-box Tauberian or depoissonization theorem** using only (i) continuous-time fixed-$`r`$ Yule averages (Theorem 1 for all fixed $`r \ge 2`$), (ii) fixed-$`\delta`$ harmonic averages (Theorem 3(2)), (iii) exact Poissonization $`\Pi_k(n) \sim \mu A c_2 n^\mu k^{-(1+\mu)}`$ at high degrees $`k \gg \sqrt{n}`$, and (iv) $`a_{n+1}(k) \ge (1 - k/n) a_n(k)`$ can deduce $`a_n(k) \sim \mu A c_2 n^\mu k^{-(1+\mu)}`$: the explicit sequence $`b_n(k) := \mu A c_2 n^\mu k^{-(1+\mu)}(1 + \frac{1}{2}\sin(k \log n))`$ satisfies (i)–(iv) yet oscillates between $`0.5`$ and $`1.5`$ times the target. This explains structurally why a tree-level transport step (such as Theorem A's iterated-log tower) is required for high degrees $`k \gtrsim \sqrt{n}`$.

### Honest Attribution & Independence Scope
- This continuation was developed **after** reciprocal review of `degrees-astra/PROOF.md` (separate dense/sparse Tier-B-dependent lane, not exported).
- **What is shared conceptually:** Like Astra §5, step 2 of Theorem A uses a stopped likelihood comparison between the physical vertex-degree transition and a neutral Pólya urn.
- **What is strictly new and independent of Astra's proof:**
  1. **Zero Tier B dependencies:** Astra §§2–7 initializes at the Tier B small-power exact-degree profile (`hub-spectrum.md` §1 eq. (2)) up to $`m^\gamma`$ and uses $`O(\log k)`$ size-doubling steps with a multi-color urn and exact Vandermonde column weight $`h_s^{(d)}`$ plus $`\log s`$ core caps. Here we eliminate `hub-spectrum.md` completely (both §1 eq. (2) and Claim 11a!), initializing instead at a **finite set of fixed degrees** via pure Tier A Theorem 5.
  2. **$`O(\log^* k)`$ iterated-log tower & $`r`$-free core control:** Using Neon's Tier A exponential core-degree moment (Lemma 2), the factor $`s_j^\mu = r^\mu L_j^\mu`$ cancels the exact-degree scale $`r^\mu L_{j-1}^{-1}`$ uniformly in $`r = n/k`$, making the core-degree cap $`B_j = \lceil 30 \log L_{j-1} \rceil`$ completely independent of $`\log n`$ (even if $`n = e^{e^k}`$!). This allows super-exponential jumps $`L_j = \lceil (\log L_{j-1})^7 \rceil`$ in only $`J \le \log^* k`$ steps with a 1-vertex 2-color urn.
  3. **6-line Rényi-moment likelihood tail:** Replaces Bernstein–Freedman + separate drift control with a direct conditional $`\alpha`$-moment supermartingale bound on $`(dP'/dQ)^{\pm \lambda}`$ ($`\alpha \in \{1+\lambda, -\lambda\}`$).
- **What Theorem A does NOT prove:** Pointwise high-degree fluctuations, dense $`L^p`$/a.s. laws, critical mixed-Poisson laws, or sparse multiplicities (which are proved in Astra §§7–19 using multi-color urns and `hub-spectrum.md` §1 eq. (2)). Theorem A strictly targets the **pointwise discrete mean** (and its Tier A upgrades of cumulative/CT means and uniform high-degree depoissonization).

---

## 1. Exact Yule-to-Poisson Conditioning and the Depoissonization Obstruction (Proposition C)

### Lemma 1.1 (Exact Conditional Poissonization Given Terminal Yule Size)
For the continuous-time friend tree $`(T(t))_{t \ge 0}`$ started from $`T_2`$, let $`Y_2 := \lim_{s \to \infty} e^{-s}|T(s)| \sim \mathrm{Gamma}(2, 1)`$. For every $`t > 0`$, $`y > 0`$, and integer $`m \ge 2`$:

```math
\mathbb P_{T_2}\big(|T(t)| = m \;\big|\; Y_2 = y\big) = e^{-\lambda} \frac{\lambda^{m-2}}{(m-2)!}, \qquad \lambda := (e^t - 1) y.
```

Consequently, since $`Y_2`$ is independent of the embedded discrete chain $`(T_m)_{m \ge 2}`$,

```math
\mathbb E_{T_2}\big[N_{=k}(T(t)) \;\big|\; Y_2 = y\big] = \Pi_k\big((e^t - 1) y\big), \qquad \mathbb E_{T_2}\big[N_{=k}(T(t))\big] = \int_0^\infty y e^{-y} \, \Pi_k\big((e^t - 1) y\big) \, dy,
```

where $`\Pi_k(\lambda) := \sum_{m=2}^\infty a_m(k) e^{-\lambda} \frac{\lambda^{m-2}}{(m-2)!}`$.

**Proof.** For a pure-birth Yule process started from $`2`$, $`\mathbb P(|T(t)| = m) = (m-1) e^{-2t}(1 - e^{-t})^{m-2}`$. Given $`|T(t)| = m`$, the future scaled limit satisfies $`Y_2 = e^{-t} \lim_{u \to \infty} e^{-u}|T(t+u)| \sim \mathrm{Gamma}(m, e^t)`$ with density $`\frac{e^{mt}}{(m-1)!} y^{m-1} e^{-e^t y}`$. Dividing their product by the marginal $`\mathrm{Gamma}(2, 1)`$ density $`y e^{-y}`$ yields

```math
\frac{(m-1) e^{-2t} (1 - e^{-t})^{m-2} \cdot \frac{e^{mt}}{(m-1)!} y^{m-1} e^{-e^t y}}{y e^{-y}} = \frac{((e^t - 1) y)^{m-2}}{(m-2)!} e^{-(e^t - 1) y}.
```

Independence of $`(T_m)_{m \ge 2}`$ from the holding times (hence from $`(|T(t)|, Y_2)`$) yields both expectation identities. $`\square`$

### Proof of Proposition C (One-Step Regularity and High-Degree Depoissonization Obstruction)
1. *(Exact one-step discrete regularity)* For every $`n \ge 2`$ and $`k \ge 3`$, $`N_{=k}(T_{n+1}) - N_{=k}(T_n)`$ increases by $`1`$ iff a degree-$`(k-1)`$ vertex receives the new leaf and decreases by $`1`$ iff a degree-$`k`$ vertex receives it, so

```math
a_{n+1}(k) - a_n(k) = \frac{\mathbb E[R_{k-1}(T_n)] - \mathbb E[R_k(T_n)]}{n} \ge -\frac{\mathbb E[R_k(T_n)]}{n} \ge -\frac{k}{n} a_n(k),
```

since $`R_{k-1} \ge 0`$ and $`R_k(T_n) = \sum_{D_v = k}(k - c_v + \rho_v) \le \sum_{D_v = k}(k - c_v/2) \le k N_{=k}(T_n)`$.
- *Remark (Low-degree sufficiency of Poissonization):* If $`k \ll \sqrt{n}`$, then over $`|m - n| \le M \sqrt{n}`$, $`\prod_{j=n}^{n+M\sqrt{n}} (1 - k/j) = 1 - O(Mk/\sqrt{n}) \to 1`$, so the one-step lower bound $`a_{n+1}(k) \ge (1 - k/n) a_n(k)`$ immediately depoissonizes $`\Pi_k(\lambda) \sim \mu A c_2 \lambda^\mu k^{-(1+\mu)}`$ to $`a_n(k) \sim \mu A c_2 n^\mu k^{-(1+\mu)}`$ for $`k \ll \sqrt{n}`$.
2. *(High-degree obstruction)* Consider $`b_m(k) := \mu A c_2 m^\mu k^{-(1+\mu)} (1 + \frac{1}{2}\sin(k \log m)) \mathbf 1_{\{m > k\}}`$.
   - For $`m > k`$, $`\frac{d}{d\theta}\log(1 + \frac{1}{2}\sin\theta) = \frac{\frac{1}{2}\cos\theta}{1 + \frac{1}{2}\sin\theta} \in [-\frac{1}{\sqrt{3}}, \frac{1}{\sqrt{3}}]`$, so for all large $`k`$ and $`m/k`$, $`b_{m+1}(k) \ge (1 - \frac{2}{3} k/m) b_m(k) \ge (1 - k/m) b_m(k)`$.
   - For every fixed $`\delta > 0`$, $`\sum_{m=n}^{\lfloor (1+\delta)n \rfloor} \frac{1}{m} m^\mu \sin(k \log m) = \int_n^{(1+\delta)n} x^{\mu - 1} \sin(k \log x)\,dx + O(k n^{\mu - 1}) = O(n^\mu / k + k n^{\mu-1}) = o(n^\mu)`$, so $`b_m(k)`$ matches Theorem 3(2) for every fixed $`\delta > 0`$.
   - For every fixed $`r \ge 2`$, summing by parts on $`\mathbb E[|T^{(r)}(t)|^\mu \sin(k \log |T^{(r)}(t)|)]`$ gives $`O_r(e^{\mu t} / k) = o(e^{\mu t})`$ (since $`\sum_{m \ge r} |\mathbb P(|T^{(r)}(t)| = m+1) - \mathbb P(|T^{(r)}(t)| = m)| = O(e^{-t})`$ while the phase $`k \log m`$ oscillates with period $`O(e^t / k)`$), so $`b_m(k)`$ matches Theorem 1 for every fixed $`r \ge 2`$.
   - When $`k \gg \sqrt{n}`$, the Poisson distribution $`\mathrm{Poisson}(n)`$ has standard deviation $`\sqrt{n}`$ while $`k \log m`$ has local period $`2\pi n / k \ll \sqrt{n}`$, so summing by parts on $`\sum_m e^{-n} \frac{n^{m-2}}{(m-2)!} m^\mu \sin(k \log m)`$ gives $`O(n^\mu \cdot \frac{\sqrt{n}}{k}) = o(n^\mu)`$, hence $`\Pi_k^{(b)}(n) \sim \mu A c_2 n^\mu k^{-(1+\mu)}`$ too!
   - Yet $`b_n(k) / (\mu A c_2 n^\mu k^{-(1+\mu)}) = 1 + \frac{1}{2}\sin(k \log n) \in [0.5, 1.5]`$ does not converge to $`1`$. $`\square`$

---

## 2. Elementary 1-Vertex Neutral Urn & Rényi-Moment Lemmas

Let $`q_{s \to t}(a, l) := \mathbb P(\mathrm{Urn}(a, s - a; t) = l)`$ be the transition probability of a neutral 2-color Pólya urn from $`(a, s - a)`$ at total size $`s`$ to $`(l, t - l)`$ at total size $`t > s`$:

```math
q_{s \to t}(a, l) = \frac{\binom{l - 1}{a - 1} \binom{t - l - 1}{s - a - 1}}{\binom{t - 1}{s - 1}}, \qquad 1 \le a \le l, \quad 0 \le l - a \le t - s.
```

### Lemma 2.1 (Exact Beta-Binomial Tails & Hypergeometric Backward Expectation)
Let $`1 \le a < s < t`$.
1. *(Exact Beta-Binomial representation & two-sided multiplicative tails)* $`U_t \sim \mathrm{Urn}(a, s - a; t)`$ satisfies $`U_t - a \mid \Theta \sim \mathrm{Bin}(t - s, \Theta)`$ with $`\Theta \sim \mathrm{Beta}(a, s - a)`$. Consequently, if $`a/s \le 1/4`$ and $`\eta \in (0, 1/2)`$,

```math
\mathbb P\!\left(U_t \ge (1 + \eta)\frac{t}{s} a\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{64}\right), \qquad \mathbb P\!\left(U_t \le (1 - \eta)\frac{t}{s} a\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{64}\right).
```

2. *(Maximal ratio tail)* For $`U_u \sim \mathrm{Urn}(a, s - a; u)`$ ($`s \le u \le t`$), if $`a/s \le 1/6`$:

```math
\mathbb P\!\left(\max_{s \le u \le t} \frac{U_u}{u} \ge \frac{1}{2}\right) \le \left(\frac{4a}{s}\right)^a.
```

3. *(Exact Hypergeometric backward column & inverse-power expectation)* For every $`1 \le l < t`$,

```math
\frac{t - 1}{s - 1} q_{s \to t}(d, l) = \frac{\binom{l - 1}{d - 1} \binom{t - l - 1}{s - d - 1}}{\binom{t - 2}{s - 2}} = \mathbb P\big(D = d\big),
```

where $`D - 1 \sim \mathrm{Hypergeom}(t - 2, l - 1, s - 2)`$ has mean $`\mu_D := 1 + (s - 2)\frac{l - 1}{t - 2}`$ and variance $`\mathrm{Var}(D) \le \mu_D`$. For any interval $`I = [d_-, d_+] \subset [1, l]`$ such that $`d_- \ge \frac{1}{2}\mu_D`$ and $`\min(\mu_D - d_-, d_+ - \mu_D) \ge \gamma \mu_D`$ with $`\gamma \in (0, 1/2)`$:

```math
\mathbb P(D \notin I) \le 2 \exp\!\left(-\frac{\gamma^2 \mu_D}{4}\right),
```

and there is a universal constant $`C_\mu < \infty`$ (depending only on $`\mu`$) such that whenever $`\mu_D \ge 4`$ and $`\exp(-\gamma^2 \mu_D / 8) \le \mu_D^{-4}`$,

```math
\left|\frac{\sum_{d \in I} q_{s \to t}(d, l) \, s^\mu d^{-(1+\mu)}}{t^\mu l^{-(1+\mu)}} - 1\right| \le C_\mu \left(\frac{1}{\mu_D} + \frac{1}{s}\right).
```

**Proof.**
1. Direct Beta integration gives $`\mathbb E[\binom{t-s}{l-a}\Theta^{l-a}(1-\Theta)^{t-s-l+a}] = q_{s \to t}(a, l)`$. Write $`\Theta \stackrel{d}{=} G_1/(G_1 + G_2)`$ with independent $`G_1 \sim \mathrm{Gamma}(a, 1)`$ and $`G_2 \sim \mathrm{Gamma}(s - a, 1)`$. For $`\xi := \eta/4 \in (0, 1/8)`$ and $`a/s \le 1/4`$:
   - On $`\{G_1 < (1+\xi)a\} \cap \{G_2 > (1-\xi)(s-a)\}`$, we have $`\Theta < \frac{1+\xi}{1 - \xi + 2\xi(a/s)}\frac{a}{s} \le \frac{1+\xi}{1-\xi/2}\frac{a}{s} \le (1 + 2\xi)\frac{a}{s} = (1 + \frac{\eta}{2})\frac{a}{s}`$. By Gamma Chernoff ($`\mathbb P(G_1 \ge (1+\xi)a) \le e^{-a\xi^2/3}`$ and $`\mathbb P(G_2 \le (1-\xi)(s-a)) \le e^{-(s-a)\xi^2/2} \le e^{-a\xi^2/2}`$), $`\mathbb P(\Theta \ge (1+\frac{\eta}{2})\frac{a}{s}) \le 2 e^{-a\eta^2/48}`$. On $`\{\Theta < (1+\frac{\eta}{2})\frac{a}{s}\}`$, $`a + \mathrm{Bin}(t-s, \Theta)`$ is dominated by $`a + B`$ where $`B \sim \mathrm{Bin}(t-s, (1+\frac{\eta}{2})\frac{a}{s})`$ has mean $`\nu = (1+\frac{\eta}{2})\frac{t-s}{s}a \le \frac{3}{2}\frac{t}{s}a`$ and target excess $`\Delta := (1+\eta)\frac{t}{s}a - a - \nu = \frac{\eta}{2}(\frac{t}{s}+1)a \ge \frac{\eta}{2}\frac{t}{s}a`$. Standard Binomial Bernstein/Chernoff $`\mathbb P(B \ge \nu + \Delta) \le \exp(-\frac{\Delta^2}{2(\nu + \Delta/3)})`$ with $`\nu + \Delta/3 \le 2\frac{t}{s}a`$ gives $`\le \exp(-\frac{\eta^2}{16}\frac{t}{s}a) \le \exp(-\frac{\eta^2 a}{16})`$.
   - Symmetrically, on $`\{G_1 > (1-\xi)a\} \cap \{G_2 < (1+\xi)(s-a)\}`$, $`\Theta > \frac{1-\xi}{1+\xi}\frac{a}{s} \ge (1 - \frac{\eta}{2})\frac{a}{s}`$; on this event, $`a + \mathrm{Bin}(t-s, \Theta)`$ dominates $`a + \mathrm{Bin}(t-s, (1-\frac{\eta}{2})\frac{a}{s})`$ with mean $`\nu' + a = (1-\frac{\eta}{2})\frac{t}{s}a + \frac{\eta}{2}a`$, so the shortfall to $`(1-\eta)\frac{t}{s}a`$ is $`\Delta' = \frac{\eta}{2}(\frac{t}{s} + 1)a`$, and standard lower-tail Binomial Chernoff $`\mathbb P(B \le \nu' - \Delta') \le \exp(-\frac{(\Delta')^2}{2\nu'})`$ with $`\nu' \le \frac{t}{s}a`$ gives $`\le \exp(-\frac{\eta^2}{8}\frac{t}{s}a) \le \exp(-\frac{\eta^2 a}{8})`$.
   - Combining yields $`\le 3 \exp(-\eta^2 a / 64)`$ for both tails.
2. Under the neutral Pólya urn, $`U_u / u = \mathbb E[\Theta \mid \mathcal G_u]`$ is a nonnegative martingale, so $`(U_u/u)^a`$ is a nonnegative submartingale. By Doob's maximal inequality and Jensen's inequality:

```math
\mathbb P\!\left(\max_{s \le u \le t} \frac{U_u}{u} \ge \frac{1}{2}\right) \le 2^a \mathbb E[\Theta^a] = 2^a \frac{a(a+1)\cdots(2a-1)}{s(s+1)\cdots(s+a-1)} \le 2^a \left(\frac{2a}{s}\right)^a = \left(\frac{4a}{s}\right)^a.
```

3. Direct factorial cancellation gives $`\frac{t-1}{s-1} q_{s \to t}(d, l) = \mathbb P(D = d)`$ with $`D - 1 \sim \mathrm{Hypergeom}(t - 2, l - 1, s - 2)`$, mean $`\mu_D = 1 + (s-2)\frac{l-1}{t-2} = \frac{s}{t} l \cdot (1 + O(\frac{1}{\mu_D} + \frac{1}{s}))`$, and $`\mathrm{Var}(D) = (s-2)\frac{l-1}{t-2}\frac{t-l-1}{t-2}\frac{t-s}{t-3} \le \mu_D`$. By Hoeffding's convex-order comparison between hypergeometric and binomial sampling, $`D - 1`$ inherits the two-sided multiplicative Binomial Chernoff bounds at deviation $`\gamma \mu_D`$ (noting $`\mu_D - 1 \ge \frac{3}{4}\mu_D`$ for $`\mu_D \ge 4`$), giving $`\mathbb P(D \notin I) \le 2 \exp(-\gamma^2 \mu_D / 4)`$. For the inverse-power sum:

```math
\sum_{d \in I} q_{s \to t}(d, l) \, s^\mu d^{-(1+\mu)} = \frac{s - 1}{t - 1} s^\mu \, \mathbb E\big[D^{-(1+\mu)} \mathbf 1_{\{D \in I\}}\big].
```

- Lower bound: by convexity of $`x \mapsto x^{-(1+\mu)}`$ on $`(0, \infty)`$ and Jensen's inequality,

```math
\mathbb E\big[D^{-(1+\mu)} \mathbf 1_{\{D \in I\}}\big] \ge \mathbb E[D^{-(1+\mu)}] - \mathbb P(D \notin I) \ge \mu_D^{-(1+\mu)} - 2 e^{-\gamma^2 \mu_D / 4} \ge \mu_D^{-(1+\mu)} (1 - 2 \mu_D^{-2}).
```

- Upper bound: on $`\{D \in I\}`$, $`D \ge \frac{1}{2}\mu_D`$, so second-order Taylor expansion with Lagrange remainder gives

```math
D^{-(1+\mu)} \le \mu_D^{-(1+\mu)} - (1+\mu)\mu_D^{-(2+\mu)}(D - \mu_D) + C_2 \mu_D^{-(3+\mu)}(D - \mu_D)^2.
```

Multiplying by $`\mathbf 1_{\{D \in I\}}`$ and taking expectation, $`|\mathbb E[(D - \mu_D)\mathbf 1_{\{D \in I\}}]| = |\mathbb E[(D - \mu_D)\mathbf 1_{\{D \notin I\}}]| \le \sqrt{\mathrm{Var}(D)\mathbb P(D \notin I)} \le \sqrt{2\mu_D} e^{-\gamma^2 \mu_D / 8} \le 2 \mu_D^{-3}`$, while $`\mathbb E[(D-\mu_D)^2 \mathbf 1_{\{D \in I\}}] \le \mathrm{Var}(D) \le \mu_D`$, giving $`\mathbb E[D^{-(1+\mu)}\mathbf 1_{\{D \in I\}}] \le \mu_D^{-(1+\mu)}(1 + C_3 / \mu_D)`$.
- Finally, $`\frac{s-1}{t-1} s^\mu \mu_D^{-(1+\mu)} = \frac{s-1}{t-1} \frac{t}{s} \left(\frac{s l / t}{\mu_D}\right)^{1+\mu} t^\mu l^{-(1+\mu)} = \big(1 + O(\frac{1}{\mu_D} + \frac{1}{s})\big) t^\mu l^{-(1+\mu)}`$. $`\square`$

### Lemma 2.2 (6-Line Rényi-Moment Comparison for One Single Vertex)
Fix integers $`1 \le d < s < t`$, a core-degree cap $`B \ge 1`$, and a deterministic vertex $`v`$ in a tree $`T_s`$ satisfying $`d := D_v(s) < s/2`$ and $`c_v(s) \le B`$. Assume $`d \ge 4 B`$. Let $`\tau_{\mathrm{stop}} := \inf\{u \ge s : c_v(u) > B \text{ or } D_v(u) \ge u/2\} \wedge t`$.
Construct coupled laws $`(P', Q)`$ on histories from $`s`$ to $`t`$ conditional on $`\mathcal F_s`$:
- While $`u < \tau_{\mathrm{stop}}`$, let $`Z_u := \mathbf 1_{\{v \text{ receives leaf } u+1\}}`$; under $`P'`$, sample $`Z_u \sim \mathrm{Bernoulli}(p_v(u))`$ with physical probability $`p_v(u) = (D_v(u) - b_v(u))/u`$ ($`b_v(u) := c_v(u) - \rho_v(u) \in [0, c_v(u)]`$); under $`Q`$, sample $`Z_u \sim \mathrm{Bernoulli}(q_v(u))`$ with $`q_v(u) = D_v(u)/u`$; given $`Z_u = 1`$ attach leaf $`u+1`$ to $`v`$, and given $`Z_u = 0`$ sample the recipient $`w \ne v`$ from the identical conditional physical distribution $`\mathbb P_{\mathrm{phys}}(\cdot \mid \text{recipient} \ne v, T_u)`$.
- For $`u \ge \tau_{\mathrm{stop}}`$, evolve $`D_v(u)`$ alone under both $`P'`$ and $`Q`$ by $`Z_u \sim \mathrm{Bernoulli}(D_v(u)/u)`$.

Then under $`Q`$, $`(D_v(u))_{u=s}^t \stackrel{d}{=} \mathrm{Urn}(d, s - d; u)`$, while $`P'`$ agrees identically with the physical law $`P`$ up to $`\tau_{\mathrm{stop}}`$. Furthermore, for every real $`\alpha`$ with $`|\alpha| \ge 1`$ and $`|\alpha| \cdot \frac{B}{d} \le \frac{1}{4}`$ (and $`s \ge 4 d`$), setting $`S_{s, t} := \sum_{u=s}^{t-1} \frac{1}{u}`$:

```math
\mathbb E_Q\!\left[\left(\frac{dP'}{dQ}\right)^\alpha\right] \le \exp\!\left(\frac{8 \alpha^2 B^2}{d} S_{s, t}\right).
```

Consequently, for every $`\lambda \ge 1`$ with $`(1 + \lambda)\frac{B}{d} \le \frac{1}{4}`$ and every $`\delta > 0`$, writing $`\varepsilon_{\mathrm{Ren}} := \exp\big(-\lambda \delta + \frac{8(1+\lambda)^2 B^2}{d} S_{s, t}\big)`$, for every integer $`l`$:

```math
e^{-\delta} q_{s \to t}(d, l) - \varepsilon_{\mathrm{Ren}} - \mathbb P_P(\tau_{\mathrm{stop}} < t \mid \mathcal F_s) \le \mathbb P_P\big(D_v(t) = l \;\big|\; \mathcal F_s\big) \le e^{\delta} q_{s \to t}(d, l) + \varepsilon_{\mathrm{Ren}} + \mathbb P_P(\tau_{\mathrm{stop}} < t \mid \mathcal F_s).
```

**Proof.**
Before $`\tau_{\mathrm{stop}}`$, $`D_v(u) \ge d \ge 4 B > c_v(u)`$ and $`D_v(u) < u/2`$, so $`0 < p_v(u) \le q_v(u) < 1/2`$, and the conditional transitions given $`Z_u`$ agree under $`P'`$ and $`Q`$. Hence $`\frac{dP'}{dQ} = \prod_{u=s}^{t-1} R_u`$, where $`R_u = 1`$ for $`u \ge \tau_{\mathrm{stop}}`$ and $`R_u = (1 + x_1)^{Z_u}(1 + x_0)^{1 - Z_u}`$ for $`u < \tau_{\mathrm{stop}}`$, with $`x_1 := -\frac{b_v(u)}{D_v(u)}`$ and $`x_0 := +\frac{b_v(u)}{u - D_v(u)}`$. Since $`u - D_v(u) > u/2 \ge s/2 \ge 2d`$, we have $`\max(|x_1|, |x_0|) \le \frac{B}{d}`$, so $`|\alpha x_i| \le \frac{1}{4}`$.

For any real $`|\alpha| \ge 1`$ and $`|\alpha x| \le \frac{1}{4}`$, $`|\log(1+x) - x| \le \frac{2}{3}x^2`$, so $`y := \alpha \log(1+x) = \alpha x + \rho`$ satisfies $`|\rho| \le \frac{2}{3}|\alpha|x^2 \le \frac{2}{3}\alpha^2 x^2`$, $`|y| \le \frac{1}{3} < 1`$, and $`y^2 \le 3 \alpha^2 x^2`$; using $`e^y \le 1 + y + y^2`$ for $`|y| \le 1`$ yields $`(1+x)^\alpha \le 1 + \alpha x + 4 \alpha^2 x^2`$. Taking conditional expectation under $`Q`$ and using $`q x_1 + (1-q) x_0 = 0`$:

```math
\mathbb E_Q\big[R_u^\alpha \;\big|\; \mathcal G_u\big] \le 1 + 4\alpha^2 \big(q x_1^2 + (1-q) x_0^2\big) = 1 + \frac{4\alpha^2 b_v(u)^2}{D_v(u)(u - D_v(u))} \le 1 + \frac{8 \alpha^2 B^2}{d \, u} \le \exp\!\left(\frac{8 \alpha^2 B^2}{d \, u}\right).
```

Iterating by the tower property from $`u = t-1`$ down to $`s`$ gives $`\mathbb E_Q[(dP'/dQ)^\alpha] \le \exp(\frac{8\alpha^2 B^2}{d} S_{s, t})`$.

By Markov's inequality:
- With $`\alpha = 1 + \lambda`$: $`\mathbb P_{P'}(dP'/dQ > e^\delta) \le e^{-\lambda \delta} \mathbb E_{P'}[(dP'/dQ)^\lambda] = e^{-\lambda \delta} \mathbb E_Q[(dP'/dQ)^{1+\lambda}] \le \varepsilon_{\mathrm{Ren}}`$.
- With $`\alpha = -\lambda`$: $`\mathbb P_Q(dP'/dQ < e^{-\delta}) \le e^{-\lambda \delta} \mathbb E_Q[(dP'/dQ)^{-\lambda}] \le \varepsilon_{\mathrm{Ren}}`$.

Let $`E := \{D_v(t) = l\}`$. Since $`P`$ and $`P'`$ agree identically on $`\{\tau_{\mathrm{stop}} = t\}`$ (and $`\mathbb P_{P'}(\tau_{\mathrm{stop}} < t) = \mathbb P_P(\tau_{\mathrm{stop}} < t)`$):

```math
\mathbb P_P(E \mid \mathcal F_s) \le \mathbb P_{P'}(E \cap \{dP'/dQ \le e^\delta\}) + \mathbb P_{P'}(dP'/dQ > e^\delta) + \mathbb P_P(\tau_{\mathrm{stop}} < t) \le e^\delta q_{s \to t}(d, l) + \varepsilon_{\mathrm{Ren}} + \mathbb P_P(\tau_{\mathrm{stop}} < t),
```

```math
\mathbb P_P(E \mid \mathcal F_s) \ge \mathbb P_{P'}(E \cap \{dP'/dQ \ge e^{-\delta}\}) - \mathbb P_P(\tau_{\mathrm{stop}} < t) \ge e^{-\delta} \mathbb P_Q(E \cap \{dP'/dQ \ge e^{-\delta}\}) - \mathbb P_P(\tau_{\mathrm{stop}} < t) \ge e^{-\delta} q_{s \to t}(d, l) - \varepsilon_{\mathrm{Ren}} - \mathbb P_P(\tau_{\mathrm{stop}} < t). \qquad \square
```

---

## 3. Proof of Theorem A (Pure Tier A $`O(\log^* k)`$ Iterated-Log Tower)

We give the proof from $`T_2`$; for any reachable finite seed $`\tau`$ of size $`m_0`$ with $`C(\tau) \ge 1`$, conditioning our Tier A Lemmas 1(2) and 2(2) and fixed-degree Tier A Theorem 5 on the positive-probability event $`\{T_{m_0} = \tau\}`$ yields identical constants (scaled by $`\mathbb P_{T_2}(T_{m_0} = \tau)^{-1}`$) and replaces $`c_2 = \mathbb E[X]`$ by $`\mathbb E_\tau[X]`$ verbatim.

### Step 1 (Tier A inputs used)
From the pinned manuscript (Theorems 1, 4, 5) and Neon's Tier A Lemmas 1–2 (`degrees-neon/PROOF.md` (pre-continuation snapshot, not exported; Tier A Lemmas 1–2 are reconstructed in `REVIEW.md` §2)):
1. $`\mu \in (1/2, 1)`$, $`e^{1/3} - 1 < \mu`$, and $`\mathbb E[C_s] \le K_0 s^\mu`$ for all $`s \ge 2`$.
2. *(Neon Lemma 1(2), Tier A)* $`\mathbb E[N_{\ge d}(T_s)] \le K_{\mathrm{disc}} s^\mu d^{-\mu}`$ for all $`2 \le d \le s`$.
3. *(Neon Lemma 2(2), Tier A)* $`\mathbb E[\sum_{v \in \text{Core}(T_s)} e^{c_v(s)/3}] \le K'_{1/3} s^\mu`$ for all $`s \ge 2`$.
4. *(Tier A Theorems 4–5)* For each fixed integer $`d \ge 2`$, $`N_{=d}(T_s)/C_s \to p_d`$ a.s., and since $`0 \le N_{=d}(T_s)/s^\mu \le C_s/s^\mu`$ which converges in $`L^2`$ to $`X`$, dominated/UI convergence gives $`\lim_{s \to \infty} a_s(d) / s^\mu = c_2 p_d`$; moreover $`p_d \sim \mu A d^{-(1+\mu)}`$ as $`d \to \infty`$.

### Step 2 (Tower parameters)
Fix an arbitrary $`\varepsilon \in (0, 1/8)`$. Choose a large threshold $`M_1 = M_1(\varepsilon) \ge 10^{70}`$ (fixed depending only on $`\varepsilon`$ and the Tier A constants above) so that:
- For every integer $`d \ge \lfloor \frac{3}{4} M_1 \rfloor`$, $`\left|\frac{p_d}{\mu A d^{-(1+\mu)}} - 1\right| \le \varepsilon`$.
- For every $`x \ge \exp(M_1^{1/7})`$, we have $`\log x \ge 10^{10}`$ (so $`\sqrt{\log x} \ge 10^5`$), $`x \ge 2^{14} \lceil (\log x)^7 \rceil`$, and all asymptotic error inequalities below hold with room to spare.
- Defining $`M_0 := \lceil \exp(M_1^{1/7}) \rceil`$ (still a fixed constant depending only on $`\varepsilon`$), choose $`R_0 = R_0(\varepsilon) \ge 24`$ so that for every integer $`s \ge R_0 M_1`$ and every $`d`$ in the **finite** set $`[\lfloor \frac{3}{4} M_1 \rfloor, \lceil \frac{5}{4} M_0 \rceil]`$:

```math
\left|\frac{a_s(d)}{s^\mu} - c_2 p_d\right| \le \varepsilon c_2 p_d.
```

Now take any $`(n, k)`$ with $`k > M_0`$ and $`r := n / k \ge R_0`$. Define the super-exponentially shrinking degree scales and matching sizes:

```math
L_0 := k, \quad s_0 := n; \qquad L_j := \big\lceil (\log L_{j-1})^7 \big\rceil, \quad s_j := \lfloor r L_j \rfloor \quad (j = 1, 2, \dots, J),
```

stopping at the first $`J \ge 1`$ for which $`L_J \le M_0`$. Because $`L_{J-1} > M_0 \ge \exp(M_1^{1/7})`$, we automatically have:

```math
M_1 \le L_J \le M_0, \qquad L_{j-1} \ge 2^{14} L_j \quad \forall j \in \{1, \dots, J\}, \qquad s_J \ge R_0 M_1.
```

Define expanding relative band radii going backward from $`j = 0`$ to $`j = J`$:

```math
\rho_0 := 0, \qquad \rho_j := \sum_{i=1}^j L_i^{-1/4}, \qquad I_j := \big[(1 - \rho_j) L_j,\; (1 + \rho_j) L_j\big] \cap \mathbb Z.
```

Since $`L_{j-1} \ge 2^{14} L_j`$, the series $`\rho_j`$ is geometrically dominated by $`2 L_J^{-1/4} \le 2 M_1^{-1/4} < 1/4`$, so $`I_j \subset [\frac{3}{4} L_j, \frac{5}{4} L_j]`$ for every $`j = 0, \dots, J`$, while $`I_0 = \{k\}`$.

### Step 3 (One-step transfer from stage $`j`$ to stage $`j - 1`$)
Fix any $`j \in \{1, \dots, J\}`$, write $`s := s_j`$, $`t := s_{j-1}`$, and fix any target degree $`l \in I_{j-1}`$. Set:

```math
B_j := \lceil 30 \log L_{j-1} \rceil, \qquad \delta_j := (\log L_{j-1})^{-1/2} \le L_j^{-1/14}, \qquad \lambda_j := (\log L_{j-1})^3.
```

Note that $`s \ge 24 L_j \ge 4 d`$ for all $`d \in I_j`$, and since $`r L_j - 1 < s \le r L_j`$ and $`r L_{j-1} - 1 < t \le r L_{j-1}`$, we have $`S_{s, t} = \sum_{u=s}^{t-1} \frac{1}{u} \le \log(t/s) + 1/s \le \log L_{j-1}`$.

Decompose the count $`N_{=l}(T_t) = \sum_{v \in T_t} \mathbf 1_{\{D_v(t) = l\}}`$ into five exhaustive pieces according to vertex status at size $`s`$:

1. **In-band, low-core vertices ($`v \in \text{Core}(T_s)`$ with $`D_v(s) \in I_j`$ and $`c_v(s) \le B_j / 2`$):**
   For each such vertex $`v`$ with $`d := D_v(s) \in I_j`$, check the hypotheses of Lemma 2.2 with $`B = B_j`$, $`\lambda = \lambda_j`$, $`\delta = \delta_j`$:
   - $`(1 + \lambda_j) B_j / d \le 2 (\log L_{j-1})^3 \cdot \frac{31 \log L_{j-1}}{\frac{3}{4} (\log L_{j-1})^7} < 90 (\log L_{j-1})^{-3} < 1/4`$.
   - $`\frac{8(1+\lambda_j)^2 B_j^2}{d} S_{s,t} \le \frac{8 \cdot 4 (\log L_{j-1})^6 \cdot 961 (\log L_{j-1})^2}{\frac{3}{4} (\log L_{j-1})^7} \log L_{j-1} \le 42000 (\log L_{j-1})^2`$, whereas $`\lambda_j \delta_j = (\log L_{j-1})^{5/2}`$. Since $`\log L_{j-1} \ge 10^{10}`$ (hence $`\sqrt{\log L_{j-1}} \ge 10^5 > 84000`$), we have $`\lambda_j \delta_j \ge 2 \cdot 42000 (\log L_{j-1})^2`$, yielding:

```math
\varepsilon_{\mathrm{Ren}} \le \exp\!\left(-\frac{1}{2} (\log L_{j-1})^{5/2}\right) \le L_{j-1}^{-10}.
```

- Bound $`\mathbb P_P(\tau_{\mathrm{stop}} < t \mid \mathcal F_s)`$:
- *Core cap failure:* Starting from $`c_v(s) \le B_j/2`$, at each step $`u \in [s, t)`$ the conditional probability that $`c_v`$ increases by $`1`$ is $`L_v(u)/(u D_v(u)) \le 1/u`$. Taking MGF parameter $`\theta = 1`$ ($`e - 1 < 2`$):

```math
\mathbb P_P\!\left(c_v(t) - c_v(s) > \frac{B_j}{2} \;\Big|\; \mathcal F_s\right) \le e^{-15 \log L_{j-1}} \prod_{u=s}^{t-1} \left(1 + \frac{e - 1}{u}\right) \le e^{-15 \log L_{j-1} + 2 S_{s, t}} \le L_{j-1}^{-13}.
```

- *Ratio boundary $`D_v(u) \ge u/2`$:* Since $`p_v(u) \le D_v(u)/u`$, $`(D_v(u))_{u=s}^t`$ is stochastically dominated by a neutral urn $`U_u \sim \mathrm{Urn}(d, s - d; u)`$. Since $`d/s \le \frac{5}{4 r} \le \frac{1}{12}`$, Lemma 2.1(2) gives

```math
\mathbb P_P\!\left(\max_{s \le u \le t} \frac{D_v(u)}{u} \ge \frac{1}{2} \;\Big|\; \mathcal F_s\right) \le \left(\frac{4d}{s}\right)^d \le 3^{-3 L_j / 4} \le L_{j-1}^{-13}.
```

- Summing Lemma 2.2 over all in-band low-core vertices and taking expectation (using Lemma 1(2): $`\mathbb E[\#\{v \in T_s : D_v(s) \in I_j\}] \le \mathbb E[N_{\ge \lfloor 3 L_j/4 \rfloor}(T_s)] \le K_{\mathrm{disc}} s^\mu (3 L_j / 4)^{-\mu} \le 2 K_{\mathrm{disc}} r^\mu`$):
the total additive error from $`\varepsilon_{\mathrm{Ren}} + \mathbb P_P(\tau_{\mathrm{stop}} < t \mid \mathcal F_s)`$ is at most

```math
2 K_{\mathrm{disc}} r^\mu \cdot 3 L_{j-1}^{-10} = 6 K_{\mathrm{disc}} L_{j-1}^{-9} \cdot \big(r^\mu L_{j-1}^{-1}\big).
```

2. **High-core vertices at size $`s`$ ($`v \in \text{Core}(T_s)`$ with $`c_v(s) > B_j / 2`$):**
   By Lemma 2(2) ($`\theta = 1/3`$):

```math
\mathbb E\big[\#\{v \in \text{Core}(T_s) : c_v(s) > B_j / 2\}\big] \le K'_{1/3} s^\mu e^{-B_j / 6} \le K'_{1/3} r^\mu L_j^\mu L_{j-1}^{-5} \le K'_{1/3} L_{j-1}^{-3} \cdot \big(r^\mu L_{j-1}^{-1}\big),
```

since $`L_j^\mu \le L_{j-1}`$. This bounds BOTH the physical contribution of high-core vertices in the upper bound AND the omitted $`q_{s \to t}(d, l) \le 1`$ terms in the lower bound!

3. **Lower wrong-side core vertices at size $`s`$ ($`v \in \text{Core}(T_s)`$ with $`2 \le D_v(s) < (1 - \rho_j) L_j`$):**
   Conditional on $`\mathcal F_s`$, each such vertex is dominated by $`U_t \sim \mathrm{Urn}(d_*, s - d_*; t)`$ where $`d_* := \lfloor (1 - \rho_j) L_j \rfloor`$. Since $`l \ge (1 - \rho_{j-1}) L_{j-1} = (1 - \rho_j + L_j^{-1/4}) L_{j-1}`$ while $`\frac{t}{s} d_* \le (1 + O(1/L_j))(1 - \rho_j) L_{j-1}`$, we have $`l \ge (1 + \frac{1}{2} L_j^{-1/4}) \frac{t}{s} d_*`$. By Lemma 2.1(1) (with $`\eta = \frac{1}{2} L_j^{-1/4}`$ and $`d_* \ge \frac{1}{2} L_j`$):

```math
\mathbb P_P\big(D_v(t) \ge l \;\big|\; \mathcal F_s\big) \le 3 \exp\!\left(-\frac{L_j^{1/2}}{512}\right) = 3 \exp\!\left(-\frac{(\log L_{j-1})^{7/2}}{512}\right) \le L_{j-1}^{-10}.
```

Summing over all $`v \in \text{Core}(T_s)`$ ($`\mathbb E[C_s] \le K_0 s^\mu \le K_0 r^\mu L_j^\mu \le K_0 r^\mu L_{j-1}`$) gives expected contribution at most $`K_0 L_{j-1}^{-8} \cdot (r^\mu L_{j-1}^{-1})`$.

4. **Vertices promoted after size $`s`$ (at a unique step $`u \in \{s + 1, \dots, t\}`$):**
   On $`E_u := \{C_u - C_{u-1} = 1\}`$, the unique promoted vertex $`v_u`$ has $`D_{v_u}(u) = 2 \le d_*`$, so conditional on $`\mathcal F_u`$, $`D_{v_u}(t)`$ is dominated by $`\mathrm{Urn}(2, u - 2; t) \preceq \mathrm{Urn}(d_*, s - d_*; t)`$ (by coupling with common uniforms, since $`2 \le d_* \le U_u`$). Thus $`\mathbb P_P(D_{v_u}(t) \ge l \mid \mathcal F_u) \le L_{j-1}^{-10}`$. Summing over $`u = s+1, \dots, t`$ and taking expectation gives at most

```math
\mathbb E[C_t - C_s] \cdot L_{j-1}^{-10} \le K_0 t^\mu L_{j-1}^{-10} \le K_0 L_{j-1}^{-8} \cdot \big(r^\mu L_{j-1}^{-1}\big).
```

5. **Upper wrong-side core vertices at size $`s`$ ($`v \in \text{Core}(T_s)`$ with $`D_v(s) > (1 + \rho_j) L_j`$ and $`c_v(s) \le B_j / 2`$):**
   The event $`\{c_v(t) > B_j\}`$ has conditional probability $`\le L_{j-1}^{-13}`$ (as in Piece 1). On $`\mathcal C_v := \{\max_{s \le u \le t} c_v(u) \le B_j\}`$, $`p_v(u) \ge (D_v(u) - B_j)/u`$, so coupling with common uniforms gives $`D_v(u) - B_j \ge W_u`$ where $`W_u \sim \mathrm{Urn}(d^*, s - d^*; u)`$ starts at $`d^* := \lceil (1 + \rho_j) L_j \rceil - B_j \ge (1 + \rho_{j-1} + \frac{3}{4} L_j^{-1/4}) L_j`$. Hence $`l - B_j \le (1 + \rho_{j-1}) L_{j-1} \le (1 - \frac{1}{4} L_j^{-1/4}) \frac{t}{s} d^*`$, so Lemma 2.1(1) gives $`\mathbb P(D_v(t) \le l,\; \mathcal C_v \mid \mathcal F_s) \le 3 \exp(-L_j^{1/2}/2048) \le L_{j-1}^{-10}`$. Summing over all $`v`$ with $`D_v(s) > (1 + \rho_j) L_j`$ (expected count $`\le K_{\mathrm{disc}} r^\mu`$ by Lemma 1(2)) gives at most $`2 K_{\mathrm{disc}} L_{j-1}^{-9} \cdot (r^\mu L_{j-1}^{-1})`$.

**Combining Pieces 1–5:** There is a constant $`C_0`$ (depending only on Tier A constants) such that uniformly over all $`l \in I_{j-1}`$,

```math
e^{-\delta_j} \sum_{d \in I_j} q_{s_j \to s_{j-1}}(d, l) \, a_{s_j}(d) - C_0 L_{j-1}^{-3} r^\mu L_{j-1}^{-1} \le a_{s_{j-1}}(l) \le e^{\delta_j} \sum_{d \in I_j} q_{s_j \to s_{j-1}}(d, l) \, a_{s_j}(d) + C_0 L_{j-1}^{-3} r^\mu L_{j-1}^{-1}. \qquad (T)
```

### Step 4 (Telescoping the $`O(\log^* k)`$ tower)
For each $`j \in \{0, \dots, J\}`$, define

```math
m_j := \min_{d \in I_j} \frac{a_{s_j}(d)}{s_j^\mu d^{-(1+\mu)}}, \qquad M_j := \max_{d \in I_j} \frac{a_{s_j}(d)}{s_j^\mu d^{-(1+\mu)}}.
```

For any $`l \in I_{j-1}`$, $`\mu_D = 1 + (s_j - 2)\frac{l - 1}{s_{j-1} - 2} \in [(1 - \rho_{j-1} - O(L_j^{-1})) L_j,\; (1 + \rho_{j-1} + O(L_j^{-1})) L_j]`$, so the distance from $`\mu_D`$ to $`\partial I_j`$ is at least $`\frac{1}{2} L_j^{3/4} \ge \frac{1}{3} L_j^{-1/4} \mu_D`$. Applying Lemma 2.1(3) (with $`\gamma = \frac{1}{3} L_j^{-1/4}`$, so $`\gamma^2 \mu_D \ge \frac{1}{12} L_j^{1/2} \gg \log L_j`$) to (T) and noting $`r^\mu L_{j-1}^{-1} \le 4 s_{j-1}^\mu l^{-(1+\mu)}`$ on $`I_{j-1}`$ gives:

```math
M_{j-1} \le e^{\alpha_j} M_j + \beta_j, \qquad m_{j-1} \ge e^{-\alpha_j} m_j - \beta_j,
```

where $`\alpha_j := \delta_j + C_1 L_j^{-1} \le 2 L_j^{-1/14}`$ and $`\beta_j := 4 C_0 L_{j-1}^{-3}`$.

Unrolling from $`j = J`$ down to $`j = 0`$:

```math
M_0 \le e^{\sum_{j=1}^J \alpha_j} \left(M_J + \sum_{j=1}^J \beta_j\right), \qquad m_0 \ge e^{-\sum_{j=1}^J \alpha_j} m_J - \sum_{j=1}^J \beta_j.
```

Since $`L_{j-1} \ge 2^{14} L_j`$ and $`L_J \ge M_1`$, we have $`\sum_{j=1}^J \alpha_j \le 4 M_1^{-1/14} \le \varepsilon`$ and $`\sum_{j=1}^J \beta_j \le 8 C_0 M_1^{-3} \le \varepsilon \mu A c_2`$ (for $`M_1(\varepsilon)`$ sufficiently large). At the bottom stage $`j = J`$, our choice of $`M_1, M_0, R_0`$ gives $`(1 - \varepsilon)^2 \mu A c_2 \le m_J \le M_J \le (1 + \varepsilon)^2 \mu A c_2`$. Since $`I_0 = \{k\}`$ and $`s_0 = n`$, we have $`m_0 = M_0 = \frac{a_n(k)}{n^\mu k^{-(1+\mu)}}`$, hence

```math
\Big(e^{-\varepsilon}(1-\varepsilon)^2 - \varepsilon\Big) \mu A c_2 \le \frac{a_n(k)}{n^\mu k^{-(1+\mu)}} \le e^\varepsilon \Big((1+\varepsilon)^2 + \varepsilon\Big) \mu A c_2
```

for all $`k > M_0(\varepsilon)`$ and $`n/k \ge R_0(\varepsilon)`$. Sending $`\varepsilon \downarrow 0`$ completes the proof of Theorem A. $`\square`$

---

## 4. Proof of Corollary B (Pure Tier A Upgrades & Effective Depoissonization)

1. **Cumulative mean (Pure Tier A):** Fix any $`M > 1`$. For $`k \to \infty`$ with $`k/n \to 0`$, every $`j \in [k, \lfloor M k \rfloor)`$ satisfies $`j \to \infty`$ and $`j/n \le M k / n \to 0`$, so Theorem A gives uniformly

```math
\sum_{j=k}^{\lfloor M k \rfloor - 1} a_n(j) = (1 + o(1)) \mu A c_2 n^\mu \sum_{j=k}^{\lfloor M k \rfloor - 1} j^{-(1+\mu)} = (1 + o(1)) A c_2 n^\mu k^{-\mu} (1 - M^{-\mu}).
```

Meanwhile, Neon's Tier A Lemma 1(2) bounds the tail $`\mathbb E[N_{\ge \lfloor M k \rfloor}(T_n)] \le K_{\mathrm{disc}} n^\mu \lfloor M k \rfloor^{-\mu} \le 2 K_{\mathrm{disc}} M^{-\mu} n^\mu k^{-\mu}`$. Sending $`k, n/k \to \infty`$ and then $`M \to \infty`$ yields $`\mathbb E[N_{\ge k}(T_n)] \sim A c_2 n^\mu k^{-\mu}`$ (and identically from any reachable seed $`\tau`$ with $`C(\tau) \ge 1`$).

2. **Continuous-time mean (Pure Tier A):** For $`T_2`$ or any reachable seed $`\tau`$ of size $`m_0`$ with $`C(\tau) \ge 1`$, the continuous-time cumulative mean satisfies

```math
e^{-\mu t} k^\mu \mathbb E_\tau[N_{\ge k}(T(t))] = \mathbb E_\tau\!\left[\big(e^{-t}|T(t)|\big)^\mu \cdot \frac{\mathbb E_\tau[N_{\ge k}(T_{|T(t)|})]}{|T(t)|^\mu k^{-\mu}} \mathbf 1_{\{|T(t)| > k\}}\right].
```

By part (1), the inner ratio converges to $`A \mathbb E_\tau[X]`$ a.s. as $`k \to \infty, k e^{-t} \to 0`$ (since $`e^{-t}|T(t)| \to Y_{m_0} \sim \mathrm{Gamma}(m_0, 1) > 0`$ a.s., so $`k / |T(t)| \to 0`$ a.s.), and by Neon's Tier A Lemma 1(2) the integrand is bounded by $`K_{\mathrm{disc}, \tau} (e^{-t}|T(t)|)^\mu`$, which is uniformly integrable (its $`L^2`$ norm is bounded uniformly in $`t \ge 0`$ by the geometric moments of $`\mathrm{NegBin}(m_0, e^{-t})`$). Hence $`\mathbb E_\tau[N_{\ge k}(T(t))] \sim A \mathbb E_\tau[X] \frac{\Gamma(m_0+\mu)}{\Gamma(m_0)} e^{\mu t} k^{-\mu} = A d_2 \Psi(\tau) e^{\mu t} k^{-\mu}`$ (with $`d_2 \Psi(T_2) = \Gamma(2+\mu) c_2`$), purely from Tier A. Plugging this Pure Tier A cumulative mean into Steps 1–4 of Neon's Theorem 1 (`degrees-neon/PROOF.md` (pre-continuation snapshot, not exported; Tier A Lemmas 1–2 are reconstructed in `REVIEW.md` §2)) immediately upgrades Neon's continuous-time exact-degree mean to **pure Tier A**.

3. **Uniform high-degree depoissonization:** Comparing Theorem A with part (2) at $`t = \log n`$ immediately gives $`a_n(k) / \mathbb E_{T_2}[N_{=k}(T(\log n))] \to 1/\Gamma(2+\mu)`$ uniformly for $`k \to \infty, k/n \to 0`$. For the exact Poissonization $`\Pi_k(n) = \mathbb E[a_{2 + \mathrm{Pois}(n)}(k)]`$:
   - On $`\{|2 + \mathrm{Pois}(n) - n| \le n/2\}`$, $`m \in [n/2, 3n/2]`$ satisfies $`k/m \to 0`$ uniformly, so Theorem A gives $`a_m(k) = (1 + o(1)) \mu A c_2 m^\mu k^{-(1+\mu)}`$; since $`(2 + \mathrm{Pois}(n))/n \to 1`$ in $`L^2`$, the contribution of this event to $`\Pi_k(n) / (\mu A c_2 n^\mu k^{-(1+\mu)})`$ tends to $`1`$.
   - On $`\{|2 + \mathrm{Pois}(n) - n| > n/2\}`$, Poisson Chernoff gives $`\mathbb P(|2 + \mathrm{Pois}(n) - n| > n/2) \le 2 e^{-c n}`$, while deterministically $`a_m(k) \le m`$; since $`\mathbb E[(2+\mathrm{Pois}(n))\mathbf 1_{\{|2+\mathrm{Pois}(n)-n| > n/2\}}] \le C n e^{-c n/2} = o(n^\mu k^{-(1+\mu)})`$ (as $`k < n`$), the tail is negligible.
   Thus $`a_n(k) / \Pi_k(n) \to 1`$ uniformly over all $`k \to \infty, k/n \to 0`$. $`\square`$

---

## 5. Exact Finite Replay Evidence

Run `verify_tier_a_discrete_mean.py` (internal finite checker, not exported; see [`audit_exact.py`](audit_exact.py) for the exported reviewer exact-rational controls) (tests finite algebraic, probabilistic, and inequality identities only — not the asymptotic theorem).
