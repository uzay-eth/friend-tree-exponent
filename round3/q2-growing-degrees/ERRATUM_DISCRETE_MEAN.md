> **Public checkpoint export note:** Edited export of `degrees-neon/ERRATUM_DISCRETE_MEAN.md` (frozen original SHA-256 `b998a84ad6b5bcee412b8ea9cc95e9d90c3cc5af6e1f7cc801aa4b0007a56ef5`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

# Erratum to Frozen `CONTINUATION_DISCRETE_MEAN.md` (Pure Tier A Discrete-Time Mean)

**Author:** Gemini 4 Neon (`degrees-neon`).
**Target frozen file (unmodified):** [`CONTINUATION_DISCRETE_MEAN.md`](CONTINUATION_DISCRETE_MEAN.md), SHA-256 `d3ba5f9147e37e133810ee001fb4945c335068e1469f76e714b9903eddb0e456`.
**Triggering audit:** [`REVIEW.md`](REVIEW.md), SHA-256 `096fa7869f8a7d7d2ed477edf0ca2d2ed47557895f01f86728b37fb146c76ce2` (with the explicit counterexample recorded by the audit).
**Verification replay:** `verify_errata_and_gaussian.py` (author script, not exported; see independent exact-rational check [`audit_exact.py`](audit_exact.py)) (Task 1 suite).

---

## 1. Verdict & Scope of Erratum

- **Frozen line 92 (Lemma 2.1(1) upper Beta ratio bound):** **`WRONG` as written.** The explicit counterexample recorded by the audit ($`a=1, s=100, \eta=2/5, \xi=1/10, G_1=1099/1000, G_2=89199/1000`$) satisfies $`G_1 < (1+\xi)a`$ and $`G_2 > (1-\xi)(s-a)`$, yet $`G_1/(G_1+G_2) - (1+\eta/2)(a/s) = 964/5643625 > 0`$, because bounding $`a/s \le 1/4`$ in the denominator $`1 - \xi + 2\xi(a/s)`$ goes in the wrong direction when $`a/s`$ is small.
- **Two literal rounding inequalities in Theorem A Step 3:** `WRONG` as literal inequalities (`(\log L_{j-1})^{-1/2} \le L_j^{-1/14}` reverses the ceiling $`L_j = \lceil (\log L_{j-1})^7 \rceil`$; `d/s \le 5/(4r)` ignores the floor in $`s = \lfloor r L_j \rfloor`$).
- **All headline claims of `CONTINUATION_DISCRETE_MEAN.md` (Theorem A for $`T_2`$ and every fixed reachable seed $`\tau`$, Corollary B, Proposition C):** **`HOLD`** after substituting the repaired proofs and explicit constants below. No Tier B input is introduced; the fixed-degree handoff and uniformity as $`k \to \infty, n/k \to \infty`$ arbitrarily slowly remain intact.

---

## 2. Replacement Proof & Constants for Lemma 2.1(1) (Both Urn Tails)

Let $`1 \le a < s < t`$ with $`p := a/s \le 1/4`$ (hence $`s - a \ge 3 a`$), and let $`\eta \in (0, 1/2)`$. Set

```math
\xi := \frac{\eta}{5} \in \left(0, \frac{1}{10}\right).
```

Represent $`\Theta \sim \mathrm{Beta}(a, s - a)`$ as $`\Theta = G_1/(G_1 + G_2)`$ with independent $`G_1 \sim \mathrm{Gamma}(a, 1)`$ and $`G_2 \sim \mathrm{Gamma}(s - a, 1)`$.

1. **Upper tail $`\mathbb P(U_t \ge (1 + \eta)\frac{t}{s} a)`$:**
   - On the upper good event $`\{G_1 < (1 + \xi) a\} \cap \{G_2 > (1 - \xi)(s - a)\}`$, since $`p > 0`$ gives $`1 - \xi + 2\xi p > 1 - \xi`$:

```math
\Theta < \frac{(1 + \xi) a}{(1 + \xi) a + (1 - \xi)(s - a)} = \frac{1 + \xi}{1 - \xi + 2\xi p} \, p < \frac{1 + \xi}{1 - \xi} \, p < \left(1 + \frac{\eta}{2}\right) p,
```

using the exact polynomial identity $`(1 + \frac{\eta}{2})(1 - \frac{\eta}{5}) - (1 + \frac{\eta}{5}) = \frac{\eta(1 - \eta)}{10} > 0`$ for all $`\eta \in (0, 1)`$.
- For Gamma exponential tails at $`\xi \in (0, 1/10)`$:
- Upper Gamma rate: $`\xi - \log(1 + \xi) \ge \frac{\xi^2}{2(1 + \xi)} \ge \frac{\xi^2}{2.2} = \frac{\eta^2}{55}`$ (conservatively $`\ge \frac{\xi^2}{3} = \frac{\eta^2}{75}`$).
- Lower Gamma rate: $`-\xi - \log(1 - \xi) \ge \frac{\xi^2}{2} = \frac{\eta^2}{50}`$ (conservatively $`\ge \frac{\xi^2}{3} = \frac{\eta^2}{75}`$).
Since $`s - a \ge 3 a`$, $`\mathbb P(\Theta \ge (1 + \frac{\eta}{2}) p) \le 2 e^{-\eta^2 a / 55} \le 2 e^{-\eta^2 a / 75}`$.
- Conditional on $`\Theta < (1 + \frac{\eta}{2}) p`$, $`U_t \preceq a + B`$ where $`B \sim \mathrm{Bin}(t - s, (1 + \frac{\eta}{2}) p)`$ has mean $`\nu = (1 + \frac{\eta}{2})\frac{t - s}{s} a \le \frac{3}{2}\frac{t}{s} a`$ and excess $`\Delta := (1 + \eta)\frac{t}{s} a - a - \nu = \frac{\eta}{2}(\frac{t}{s} + 1) a \ge \frac{\eta}{2}\frac{t}{s} a`$. Binomial Bernstein gives $`\mathbb P(B \ge \nu + \Delta) \le \exp(-\frac{\Delta^2}{2(\nu + \Delta/3)}) \le e^{-\eta^2 a / 16}`$.

2. **Lower tail $`\mathbb P(U_t \le (1 - \eta)\frac{t}{s} a)`$:**
   - On the lower good event $`\{G_1 > (1 - \xi) a\} \cap \{G_2 < (1 + \xi)(s - a)\}`$, since $`1 + \xi - 2\xi p < 1 + \xi`$:

```math
\Theta > \frac{(1 - \xi) a}{(1 - \xi) a + (1 + \xi)(s - a)} = \frac{1 - \xi}{1 + \xi - 2\xi p} \, p > \frac{1 - \xi}{1 + \xi} \, p > \left(1 - \frac{\eta}{2}\right) p,
```

using the exact polynomial identity $`(1 - \frac{\eta}{5}) - (1 - \frac{\eta}{2})(1 + \frac{\eta}{5}) = \frac{\eta(1 + \eta)}{10} > 0`$ for all $`\eta > 0`$.
- By the exact same Gamma rates above (and $`s - a \ge 3a`$), $`\mathbb P(\Theta \le (1 - \frac{\eta}{2}) p) \le 2 e^{-\eta^2 a / 55} \le 2 e^{-\eta^2 a / 75}`$.
- Conditional on $`\Theta > (1 - \frac{\eta}{2}) p`$, $`U_t`$ dominates $`a + B'`$ where $`B' \sim \mathrm{Bin}(t - s, (1 - \frac{\eta}{2}) p)`$ has mean $`\nu' \le \frac{t}{s} a`$ and shortfall $`\Delta' = (a + \nu') - (1 - \eta)\frac{t}{s} a = \frac{\eta}{2}(\frac{t}{s} + 1) a \ge \frac{\eta}{2}\frac{t}{s} a`$. Lower Binomial Chernoff gives $`\mathbb P(B' \le \nu' - \Delta') \le e^{-(\Delta')^2 / (2\nu')} \le e^{-\eta^2 a / 8}`$ (and $`0`$ if $`(1 - \eta)\frac{t}{s} a < a`$).

3. **Resulting two-sided constants:**
   - **Conservative unified bound (matching Astra's audit):**

```math
\mathbb P\!\left(U_t \ge (1 + \eta)\frac{t}{s} a\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{100}\right), \qquad \mathbb P\!\left(U_t \le (1 - \eta)\frac{t}{s} a\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{100}\right).
```

- **Sharp unified bound (restoring denominator $`64`$, in fact $`55`$):**

```math
\mathbb P\!\left(U_t \ge (1 + \eta)\frac{t}{s} a\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{55}\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{64}\right), \qquad \mathbb P\!\left(U_t \le (1 - \eta)\frac{t}{s} a\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{55}\right) \le 3 \exp\!\left(-\frac{\eta^2 a}{64}\right).
```

- **Resulting off-band constants in Theorem A Step 3 Pieces 3–5 (using conservative denominator $`100`$):**
- Lower wrong-side & late promotions ($`a = d_* \ge L_j / 2`$, $`\eta = \frac{1}{2} L_j^{-1/4}`$): $`\le 3 \exp(-L_j^{1/2}/800) \le L_{j-1}^{-10}`$.
- Upper wrong-side ($`a = d^* \ge L_j / 2`$, $`\eta = \frac{1}{4} L_j^{-1/4}`$): $`\le 3 \exp(-L_j^{1/2}/3200) \le L_{j-1}^{-10}`$.
(Both hold with huge margin since $`L_j^{1/2} \ge (\log L_{j-1})^{7/2}`$ and $`\log L_{j-1} \ge 10^{10}`$.)

---

## 3. Literal Rounding Repairs & Scope Clarifications

1. **Ceiling inequality for $`\delta_j`$ (frozen line 165):** Replace $`\delta_j := (\log L_{j-1})^{-1/2} \le L_j^{-1/14}`$ with:

```math
\delta_j := (\log L_{j-1})^{-1/2} = \left(\frac{L_j}{(\log L_{j-1})^7}\right)^{1/14} L_j^{-1/14} \le 2^{1/14} L_j^{-1/14} \le \frac{3}{2} L_j^{-1/14},
```

since $`L_j = \lceil (\log L_{j-1})^7 \rceil \le 2(\log L_{j-1})^7`$. Hence $`\alpha_j := \delta_j + C_1 L_j^{-1} \le 2 L_j^{-1/14}`$ and $`\sum_{j=1}^J \alpha_j \le 4 M_1^{-1/14}`$ hold verbatim.
2. **Floor inequality for $`d/s`$ (frozen line 174):** Replace $`d/s \le \frac{5}{4 r} \le \frac{1}{12}`$ with:

```math
\frac{d}{s} \le \frac{\frac{5}{4} L_j}{r L_j - 1} = \frac{5}{4 r \big(1 - (r L_j)^{-1}\big)} < \frac{1}{18} < \frac{1}{12} \qquad (r \ge 24,\; L_j \ge M_1 \ge 10^{70}).
```

3. **Rate caveat on "effective depoissonization" and fixed-degree handoff:** Theorem A and Corollary B(3) prove **uniform asymptotic equivalence** ($`a_n(k)/\Pi_k(n) \to 1`$, $`a_n(k)/\mathbb E_{T_2}[N_{=k}(T(\log n))] \to 1/\Gamma(2+\mu)`$) as $`k \to \infty, n/k \to \infty`$ arbitrarily slowly, with **no computable all-regime remainder rate** (since the fixed-degree handoff at $`s_J \ge R_0 M_1`$ uses qualitative convergence $`a_s(d)/s^\mu \to c_2 p_d`$).
4. **Low-degree depoissonization remark in Proposition C:** For $`k = o(\sqrt{n})`$, the two-sided squeezing step applies the forward one-step bound $`a_j \ge (1 - o(1)) a_n`$ on $`[n, n + 2M\sqrt{n}]`$ inside $`\Pi_k(n + M\sqrt{n})`$ to establish both the sharp upper bound and a uniform envelope on $`[n/2, 2n]`$, and then the backward bound $`a_j \le (1 + o(1)) a_n`$ on $`[n - 2M\sqrt{n}, n]`$ inside $`\Pi_k(n - M\sqrt{n})`$ (with $`a_j \le j`$ on the exponentially small tail outside $`[n/2, 2n]`$), sending $`n \to \infty`$ then $`M \to \infty`$.

---

## 4. Replay Evidence

```text
$ python3 -B verify_errata_and_gaussian.py  # author worktree run
PASS [Task 1]: Explicit counterexample recorded by the audit confirmed [human attribution neutralized in export]; xi=eta/5 certificates (196 rational checks) & sharper 1/55 + rounding bounds verified
```
*(Finite rational certificates only; does not test asymptotic limits.)*
