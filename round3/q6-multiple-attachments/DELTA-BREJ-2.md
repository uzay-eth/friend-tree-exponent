> **Public checkpoint export note:** Edited export of `multiple-neon/DELTA-BREJ-2.md` (frozen original SHA-256 `4f5d8e2c39dc2a3d5179a8340653a76a0ad9e275d398a97198811526637cb64d`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

# Extended Delta Note 2 for Frozen `ADDENDUM-BREJ.md` (Integrating Astra's Independent Review)

**Frozen/preserved files (all untouched byte-for-byte):**
- `multiple-neon/PROOF.md` (`c7bd8636d8b36425d4f7cc328e9b50db0a36de67a31e9acf502d6a4fb679e698`)
- `multiple-neon/ADDENDUM-BREJ.md` (`a2a2f85e0738deb610fff128e19ff53d4b5a2f5566dee92f5e420bdb3b492583`, verifier commit `894bbbc5eb71db19fdf89b64be23dd8a3e836a09`)
- `multiple-neon/DELTA-BREJ.md` (`b6d1b7e6e3ec046ccac7ec099c618f0262a667924c2f190e9a6494498b3487b5`)
- `multiple-neon/ADDENDUM-BREJ-RATES.md` (`93d722635d7350dbefa1d73eeeb66351a78ce588815573122999c8d021591831`)

**Incoming independent review integrated here:** `multiple-astra/brej-review/REVIEW.md` (SHA-256 `d17dc108858dcf17bea66babf8218501e7cab8cad95c38afd19bc95c9c3e8980`) and `multiple-astra/brej-review/check.py` (SHA-256 `0cd2cd2f407af89db817bd2b7c491df7bea24e5da0ee99725f62074e8d715430`).
**Regression verifier:** `q6/verify_delta2_brej.py` (author script, not exported; see independent exact-rational audit checker [`check_delta2.py`](check_delta2.py)) (100% exact `fractions.Fraction` + `math.isqrt`, zero floats).

---

## 1. Correction of Frozen Theorem 7(4): Exact Maximum-Degree Limits, Degree Conventions, and Moving-Maximum Proof

### 1.1 Exact Defect in Frozen `ADDENDUM-BREJ.md` Theorem 7(4)
Frozen Theorem 7(4) wrote that *"every one of those four [comparators: `A-simple`, `A-multi`, `B-multi`, `B-simple-col`] has $`\max_v p_n(v) \to \max_v Z_v > 0`$, $`S_n \to \sum Z_v^2 > 0`$, and $`\max_v d_n(v)/n \to 2\max_v Z_v > 0`$ a.s."*

That common coefficient $`2\max_v Z_v`$ is **WRONG** for **`B-simple-col`** (and also differs between native multidegree and projected simple degree in **`A-multi`** and **`B-multi`**): under **`B-simple-col`**, a duplicate draw $`(v, v)`$ adds only **one** simple edge to $`v`$ (not two), so the limiting maximum normalized simple degree is **$`2z - z^2`$** (where $`z := \max_{v \ge 1} h_v \in (0, 1)`$), differing from $`2z`$ by $`z^2 > 0`$ almost surely.

### 1.2 Corrected Table Across All Six Studied $`m = 2`$ Conventions

Writing $`z := \max_{v \ge 1} Z_v \in (0, 1)`$ (resp. $`\max_{v \ge 1} h_v \in (0, 1)`$ in `B-simple-col`) for each comparator's own random limiting share vector (with $`\sum_{v \ge 1} Z_v = 1`$ a.s., hence $`z > 0`$ a.s.), and explicitly specifying each model's admissible seed class and both **native degree** $`d_n^{\mathrm{native}}(v)`$ and **projected simple degree** $`d_n^{\mathrm{simp}}(v)`$:

| Convention ($`m = 2`$) | Admissible seed class | Native degree | $`\lim_{n\to\infty} \max_{v \in [n]} \frac{d_n^{\mathrm{native}}(v)}{n}`$ a.s. | $`\lim_{n\to\infty} \max_{v \in [n]} \frac{d_n^{\mathrm{simp}}(v)}{n}`$ a.s. | $`\lim M_n`$ & $`\lim S_n`$ a.s. |
|---|---|---|---|---|---|
| **`A-simple`** | Connected simple, $`\delta \ge 2`$ | Simple degree | $`2z > 0`$ | $`2z > 0`$ | $`z > 0`$, $`\sum Z_v^2 > 0`$ |
| **`A-multi`** | Connected loopless, $`\delta \ge 2`$ (Theorem 3 domain) | Multidegree | $`2z > 0`$ | $`z > 0`$ | $`z > 0`$, $`\sum Z_v^2 > 0`$ |
| **`B-multi`** | Connected loopless, $`e_0 \ge 1`$ | Multidegree | $`2z > 0`$ | $`2z - z^2 > 0`$ | $`z > 0`$, $`\sum Z_v^2 > 0`$ |
| **`B-simple-col`** | Connected simple, $`e_0 \ge 1`$ | Simple degree | **$`2z - z^2 > 0`$** | **$`2z - z^2 > 0`$** | $`z > 0`$, $`\sum h_v^2 > 0`$ |
| **`B-rej-pair`** | Connected simple, $`n_0 \ge 2`$ | Simple degree | $`0`$ | $`0`$ | $`0`$, $`0`$ |
| **`B-rej-seq`** | Connected simple, $`n_0 \ge 2`$ | Simple degree | $`0`$ | $`0`$ | $`0`$, $`0`$ |

**Qualitative separation preserved:** Across every comparator (on its admissible seed class) and under **both** native and projected simple degree conventions, all four non-rejection comparators have strictly positive almost-sure limits for $`\max_v d_n(v)/n`$, $`M_n`$, and $`S_n`$, whereas both rejection rules (`B-rej-pair`, `B-rej-seq`) have $`\max_v d_n(v)/n \to 0`$, $`M_n \to 0`$, and $`S_n \to 0`$ a.s. while preserving $`L_n \to 1`$ a.s. and $`E_n^{\mathrm{simp}}/n \to 2`$.

### 1.3 Rigorous Moving-Maximum Proof (Handling Changing Maximizers)

Let $`d_n(v)`$ denote either the native degree or the projected simple degree of label $`v \in [n]`$ in any of the four comparators. Pad inactive labels $`v > k`$ by $`p_k(v) := 0`$ and $`X_k(v) := 0`$, where $`X_k(v) := d_{k+1}(v) - d_k(v) \in \{0, 1, 2\}`$ is the old-vertex degree increment at step $`k \to k+1`$, and $`\mu_k(v) := \mathbb{E}[X_k(v) \mid \mathcal{F}_k]`$. Every vertex $`v \in [n]`$ satisfies

```math
d_n(v) = d_v^{\mathrm{init}} + \sum_{k = n_0}^{n - 1} X_k(v), \qquad 0 \le d_v^{\mathrm{init}} \le B_0 := \max\!\big(2, \max_{u \in [n_0]} d_{n_0}(u)\big).
```

1. **Uniform Azuma–Hoeffding concentration across all $`n`$ labels:** For each $`v \in [n]`$, $`M_n(v) := \sum_{k=n_0}^{n-1} (X_k(v) - \mu_k(v))`$ is a martingale with bounded differences $`|X_k(v) - \mu_k(v)| \le 2`$ and at most $`n`$ increments. For any fixed $`\eta > 0`$, Azuma–Hoeffding and a union bound over $`v \in [n]`$ give

```math
\mathbb{P}\!\left(\max_{v \in [n]} \frac{|M_n(v)|}{n} > \eta\right) \le 2 n \exp\!\left(-\frac{\eta^2 n}{8}\right),
```

which is summable in $`n`$. By Borel–Cantelli along $`\eta = 1/r`$ ($`r \in \mathbb{N}`$),

```math
\max_{v \in [n]} \left|\frac{d_n(v)}{n} - \frac{1}{n}\sum_{k=n_0}^{n-1} \mu_k(v)\right| \le \frac{B_0}{n} + \max_{v \in [n]} \frac{|M_n(v)|}{n} \xrightarrow{\mathrm{a.s.}} 0. \qquad (\star)
```

2. **Uniform $`\ell^\infty`$ convergence of Cesàro means:** In all four comparators, Theorems 2 and 5 give $`\|p_k - Z\|_{\mathrm{TV}} \to 0`$ a.s., hence $`\sup_{v \ge 1} |p_k(v) - Z_v| \to 0`$ a.s. Evaluating $`\mu_k(v)`$:
   - **`B-simple-col` (native = simple) and `B-multi` (projected simple):** An old vertex $`v \le k`$ gains $`+1`$ simple edge iff at least one of the two i.i.d. draws hits $`v`$, so $`\mu_k(v) = 1 - (1 - p_k(v))^2 = g(p_k(v))`$ where $`g(t) := 2t - t^2`$ satisfies $`g(0) = 0`$ and $`|g'(t)| \le 2`$ on $`[0, 1]`$. Thus $`\sup_{v \ge 1} |\mu_k(v) - g(Z_v)| \le 2 \sup_{v \ge 1} |p_k(v) - Z_v| \to 0`$ a.s.
   - **`A-simple` (native = simple) and `A-multi` / `B-multi` (native multidegree):** $`\mu_k(v) = 2 p_k(v)`$, so $`\sup_{v \ge 1} |\mu_k(v) - 2 Z_v| \to 0`$ a.s.
   - **`A-multi` (projected simple, on Theorem 3's $`\delta \ge 2`$ seed domain):** Writing $`r_k(x, v) := k_{G_k}(x, v)/D_k(x)`$ (so $`p_k(v) = \frac{1}{k}\sum_{x=1}^k r_k(x, v)`$), averaging $`1 - (1 - r_k(x, v))^2 = 2 r_k(x, v) - r_k(x, v)^2`$ over uniform mediator $`x \in [k]`$ gives

```math
\mu_k^{\mathrm{simp}}(v) - p_k(v) = \frac{1}{k}\sum_{x=1}^k r_k(x, v)\big(1 - r_k(x, v)\big).
```

For every pure $`2`$-leaf $`x \in \mathcal{P}_k`$, $`r_k(x, v) \in \{0, 1\}`$, so the summand vanishes identically; for non-pure mediators it is in $`[0, 1/4]`$. Hence $`\sup_{v \ge 1} |\mu_k^{\mathrm{simp}}(v) - p_k(v)| \le 1 - P_k/k \to 0`$ a.s. by Theorem 3(1), giving $`\sup_{v \ge 1} |\mu_k^{\mathrm{simp}}(v) - Z_v| \to 0`$ a.s.
3. **Conclusion:** In each case, Cesàro averaging preserves $`\ell^\infty`$ convergence, so $`(\star)`$ yields $`\sup_{v \in [n]} |d_n(v)/n - \psi(Z_v)| \to 0`$ a.s. (with $`\psi(t) = 2t - t^2`$, $`2t`$, or $`t`$ respectively). Since $`y \mapsto \sup_v y_v`$ is $`1`$-Lipschitz in $`\ell^\infty`$ and $`\psi`$ is strictly increasing on $`[0, 1]`$ with $`\psi(0) = 0`$, $`\max_{v \in [n]} d_n(v)/n \to \psi(\max_{v \ge 1} Z_v) = \psi(z) > 0`$ a.s. $`\square`$

---

## 2. Integration of All Additional §2–4 Analytic Repairs from Astra's Review (Not Previously in `DELTA-BREJ.md`)

1. **Theorem 7, Lemma 7.1 & Step 3 — Global floor $`a = \frac{10}{49}`$ and direct relative debit bound $`R^- \le 2w`$:**
   - *Defect in frozen `ADDENDUM-BREJ.md`:* Used the finite-table minimum $`c_{\mathrm{fin}} = \frac{40866367}{200000000} = w(71, 0)`$ as a global lower bound $`c_{\mathrm{fin}} C_n \le \Phi(G_n)`$, which fails on the tail $`d \ge 339, \ell = 0`$ where $`w(d, 0) = a = \frac{10}{49} < c_{\mathrm{fin}}`$ (falsified directly by seed $`K_{340}`$, where all $`340`$ vertices sit at $`(339, 0)`$ with $`\Phi/C = a < c_{\mathrm{fin}}`$).
   - *Repair:* Replace $`c_{\mathrm{fin}}`$ by the true global floor $`a = \frac{10}{49}`$ throughout Step 3 ($`a\,\mathbb{E}[C_{n \wedge \sigma_N}] \le \mathbb{E}[\Phi(G_{n \wedge \sigma_N})]`$). In Lemma 7.1, replace the loose bound $`R^- \le \frac{2M}{c_{\mathrm{fin}}} w`$ by the exact universal bound $`0 \le R^-(d, \ell) = \frac{2\ell}{d}(w(d, \ell) - w(d, \ell - 1)) \le 2 w(d, \ell)`$ (valid for all $`d \ge 3, 0 \le \ell \le d`$ since $`w \ge 0`$ and $`\frac{2\ell}{d} \le 2`$). With $`g = \frac{1}{1-\delta}`$ and $`j = 1 - 5\delta`$, this gives

```math
R_{\mathrm{rej},\delta}(d, \ell) = g(R^+ - R^-) + (g - j) R^- \le \big(g \lambda_0 + 2(g - j)\big) w(d, \ell) =: \lambda_{\mathrm{tight}}(\delta)\, w(d, \ell),
```

with $`\lambda_{\mathrm{tight}}(1/20000) \approx 0.980649 < \lambda_{\mathrm{printed}}(1/20000) = \frac{39488432850800}{40047037208017} < \frac{99}{100}`$.
2. **Theorem 7, Lemma 7.1 — Cap vs off-diagonal $`\Delta_+ w`$ sign precision:**
   - Frozen text wrote that $`\Delta_+ w \ge 2.619 \times 10^{-5} > 0`$ across all $`57{,}624`$ states. Precisely: the $`336`$ cap states $`(d, d)`$ ($`3 \le d \le 338`$) satisfy $`\Delta_+ w(d, d) = 0`$; the strict bound $`\Delta_+ w > 2.619601 \times 10^{-5}`$ holds across all $`57{,}288`$ off-diagonal states ($`3 \le d \le 338, 0 \le \ell < d`$), while $`\Delta_+ w \ge 0`$ holds globally (which is all Lemma 7.2 requires).
3. **Proposition 1(2) — All-$`M`$ strict maximizer ordering (removing the $`M_n < 1 - 1/\sqrt{2}`$ threshold):**
   - For *any* nonuniform positive probability vector $`p`$ on $`n \ge 3`$ vertices and *any* maximizer $`i`$ ($`p_i = M`$), writing $`h(i) = 1`$ and $`h(j) = \frac{M}{1 - p_j}`$ ($`j \ne i`$), sequential sampling draws $`J \sim p`$ and includes $`i`$ with conditional probability $`h(J)`$, while tuple rejection tilts $`J`$ by $`(1 - p_J)/(1 - S)`$ with the same conditional second draw. Hence $`q_{\mathrm{seq}}(i) - q_{\mathrm{pair}}(i) = \mathrm{Cov}_p(h(J), p_J)/(1 - S) > 0`$ strictly (including when $`M`$ is tied, since $`(h(u) - h(v))(p_u - p_v) > 0`$ whenever $`p_u \ne p_v`$ and $`0`$ on ties). The maximizer multiplicative enclosure also sharpens on the lower side to $`\frac{M - S}{2(1 - M)} \le 1 - f_{\mathrm{seq}}(i) \le \frac{M - S}{2(1 - M)^2}`$.
4. **Theorem 6 Step 4 Item 4 & Theorem 7 Step 3 (cross-reference with `DELTA-BREJ.md`):**
   - In Theorem 6 Step 4 item 4, both our $`\mathbf{1}_{\{Y_N \le M\}}`$ guard (`DELTA-BREJ.md` §3) and Astra's deterministic bound $`Y_{n \wedge \sigma_M} \le \max(N, b_0 + \Delta_{\max})`$ (verified against Astra's $`K_{5, 95}`$ initialization witness where $`Y_{100} = 5 > 1 + 2 = 3`$) cleanly repair the initialization gap.
   - In Theorem 7 Step 3, replacing the single dyadic threshold $`2^{k(49/50+\varepsilon)}`$ by a strictly smaller exponent ($`49/50 + \varepsilon/2`$ in `DELTA-BREJ.md` §4, or $`49/50 + 3\varepsilon/4`$ in Astra §4) together with global floor $`a = 10/49`$ yields the simultaneous zero limit $`C_n/n^{49/50 + \varepsilon} \to 0`$ a.s.
   - In Theorem 7(2), the direct degree identity $`\sum_{u \in \mathcal{C}_n} d_n(u) = 2 e_n - 2|\mathcal{L}_n| - |\mathcal{S}_n^{(1)}| = 2n + 2 C_n + O_{\mathrm{seed}}(1)`$ gives $`\frac{1}{2e_n}\sum_{u \in \mathcal{C}_n} d_n(u) \to \frac{1}{2}`$ a.s. directly.

---

## 3. Strict Partition of Review Status (Verified Regime vs Unreviewed Rate Addenda)

| Package | Review Status | Notes |
|---|---|---|
| **Theorems 1–5** (frozen `PROOF.md`, `THEOREM5-ERRATUM-AND-STATUS.md`) | **Cross-family AI reviewed: HOLDS** (Astra / GPT; no human review yet) | Includes repaired Theorem 4 Step 2a, §6 AM–GM erratum (`2323519823/67500000`), and $`m=2`$ vs $`m\ge 3`$ contrast with Astra §17. |
| **Rejection regime (`ADDENDUM-BREJ.md` + `DELTA-BREJ.md` + `DELTA-BREJ-2.md`)** — §1 bipartiteness, corrected Proposition 1, Theorem 6 (all $`m \ge 2`$ delocalization), Theorem 7(1–3) ($`m=2`$ sublinear core & mesoscopic hubs), corrected Theorem 7(4) ($`2z - z^2`$ for `B-simple-col` & qualitative separation) | **Cross-family AI reviewed on frozen `ADDENDUM-BREJ.md` with repairs derived in Astra `brej-review/REVIEW.md` and integrated here** (no human review yet) | Astra verified every underlying mechanism (`HOLDS` for §1, Theorem 6 with init repair, Theorem 7(1–3) with floor/threshold repairs, and qualitative Theorem 7(4) with $`2z - z^2`$ correction). *(Historical pre-audit note: the text of `DELTA-BREJ-2.md` was written as an author integration of Astra's repairs + exact regression suite `q6/verify_delta2_brej.py`, and was subsequently audited `HOLDS` in [`REVIEW-DELTA-BREJ-2.md`](REVIEW-DELTA-BREJ-2.md); still awaiting human refereeing.)* |
| **Quantitative rate & one-step obstruction addendum (`ADDENDUM-BREJ-RATES.md`)** — Candidate Theorem 8 ($`\liminf (\log n) M_n \ge 1`$ / $`\ge 2`$ a.s.) and Candidate Proposition 2 ($`K_2`$-reachable balanced $`k`$-hub one-step obstruction) | **UNREVIEWED candidate** (author checks only in `q6/verify_delta_and_lower_bound.py`; not covered by Astra's `brej-review/REVIEW.md`; no human review) | Kept strictly separate from the reviewed regime above. |
