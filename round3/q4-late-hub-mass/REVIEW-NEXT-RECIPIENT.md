> **Public checkpoint export note:** Edited export of `late-mass-neon/REVIEW-NEXT-RECIPIENT-0cd50c2.md` (frozen original SHA-256 `123817cfbabdf7674518a802cbc225b8061269f14db83ef10c34f3b156cd471a`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

VALIDITY VERDICT: HOLDS (for frozen copies `NEXT-RECIPIENT.md` SHA-256 `0cd50c23647af388b9e84e8663d5e70d1a3b71a0023743a1d4aea1f7a7f3c87c` and `NEXT-RECIPIENT-AVERAGING.md` SHA-256 `f5bed4429f0f2bf5155754b9ab8b51b0a8f87a586d9aab0903ee6dd4ef27580c`, conditional on stated inputs; **HASH MISMATCH** reported on `NEXT-RECIPIENT.md` vs dispatch pin `646ca34a67defb142bfa8e56983b662f6ac48dd7319edc57a558780c5579dc33`)
INHERITED-DEPENDENCY VERDICT: CONDITIONAL (rests on Round 1/2 Tier-A/B/C baseline claims, `9be4767` Neon-reviewed Q4 notes, and two previously Neon-audited post-review inputs — `ADDENDUM-CUTOFF-THRESHOLD.md` eq. (6) and `CRITICAL-EXCESS-EVIDENCE.md` — with zero reliance on `NEXT-RECIPIENT-EVIDENCE.md`, `ADDENDUM-DIRICHLET.md`, `ADDENDUM-GAUSSIAN.md`, `CRITICAL-EXCESS-DRIFT.md`, or `Astra 20b`)

# Independent Cross-Family Audit of Astra's Next-Recipient Threshold (`NEXT-RECIPIENT.md` & `NEXT-RECIPIENT-AVERAGING.md`)

**Reviewer:** Neon lead (`late-mass-neon`, Gemini 4 Neon)  
**Frozen exact audit copies:** Verified against SHA-256 digests listed below prior to auditing.

---

## 1. Hash Pinning & Frozen-Copy Status (Mismatch Reported)

| File | Dispatch Pin SHA-256 | Observed & Frozen Copy SHA-256 | Status |
|---|---|---|---|
| `NEXT-RECIPIENT.md` | `646ca34a67defb142bfa8e56983b662f6ac48dd7319edc57a558780c5579dc33` | `0cd50c23647af388b9e84e8663d5e70d1a3b71a0023743a1d4aea1f7a7f3c87c` | **HASH MISMATCH** (file on board updated at `03:14` UTC after dispatch, 1 min after `NEXT-RECIPIENT-EVIDENCE.md` appeared at `03:13` UTC) |
| `NEXT-RECIPIENT-AVERAGING.md` | `f5bed4429f0f2bf5155754b9ab8b51b0a8f87a586d9aab0903ee6dd4ef27580c` | `f5bed4429f0f2bf5155754b9ab8b51b0a8f87a586d9aab0903ee6dd4ef27580c` | **EXACT MATCH** |

- Exact byte-for-byte copies of both board files were frozen in our scratch folder prior to auditing and verified against the observed hashes above.
- **No replay reliance:** Neither `NEXT-RECIPIENT-EVIDENCE.md` nor `check_next_recipient.py` was read, replayed, or treated as evidence. Every mathematical step—including the discounted split intervals `[0.1191026, 0.1310137]` and `[0.2523362, 0.3156958]` reported in line 66 of frozen `NEXT-RECIPIENT.md`—was audited and verified strictly by independent proof and zero-shared-code exact rational computation.
- Cross-family AI review (Gemini auditing GPT) is not human referee verification.

---

## 2. Independent Reconstruction of the Five Required Mechanisms

### 2.1 Exact Positive-Jump Identity (Eqs. (3)–(5)) and Pathwise Monotonicity in $`m`$

For any edge-seeded tree $`T_j`$ ($`j \ge m \ge 2`$), letting $`e_v = \#\{u \le m : u \sim v\}`$ be the number of early neighbors of $`v`$, the unnormalized late recipient rate is identically

```math
H_j := j \, p(m, j) = \sum_{v=1}^j \frac{d_v - e_v}{d_v} = j - \sum_{v=1}^j \frac{e_v}{d_v}.
```

When newborn $`j+1 > m`$ attaches to recipient $`w \in \{1, \dots, j\}`$:
1. The leading $`j`$ increases by $`+1`$;
2. Every $`v \ne w`$ keeps $`(d_v, e_v)`$ unchanged;
3. At recipient $`w`$, $`e_w`$ is unchanged (since $`j+1 > m`$) while $`d_w \mapsto d_w + 1`$, raising $`-\frac{e_w}{d_w}`$ by $`+\frac{e_w}{d_w(d_w+1)} =: g_j(w) \ge 0`$;
4. At newborn $`j+1`$, $`d_{j+1} = 1`$ and $`e_{j+1} = 1 - I_j`$, contributing $`-(1 - I_j) = I_j - 1`$.

Summing all four changes yields **eq. (3)** identically:

```math
H_{j+1} - H_j = I_j + g_j(w), \qquad g_j(w) = \frac{e_w}{d_w(d_w + 1)} \ge 0.
```

Dividing by $`j+1`$ and subtracting $`H_j/j = p(m, j)`$ yields $`p(m, j+1) - p(m, j) = \frac{I_j - p(m, j)}{j+1} + \frac{g_j(w_j)}{j+1}`$. Hence $`p(m, j)`$ is a bounded submartingale converging a.s. to $`R_m`$ (for fixed $`m`$), and summing over $`j \ge n`$ gives **eq. (4)**:

```math
R_m - p(m, n) = M_{m, n} + B_{m, n}, \qquad M_{m, n} = \sum_{j \ge n} \frac{I_j - p(m, j)}{j+1}, \qquad B_{m, n} = \sum_{j \ge n} \frac{g_j(w_j)}{j+1} \ge 0.
```

- **Exact pathwise monotonicity in $`m`$ at fixed $`n`$:** Fixing a single physical trajectory, increasing $`m \mapsto m'`$ ($`m \le m' \le n`$) monotonically enlarges $`\{1, \dots, m\}`$, so $`e_w(m) \le e_w(m')`$ for every $`w`$ and every $`j \ge n`$, whereas $`d_w(j)`$ and $`w_j`$ are independent of $`m`$. Thus $`g_j(w_j; m)`$ and $`B_{m,n}`$ are pathwise nondecreasing in $`m`$, while $`p(m, n) = \sum_{v > m} p_n(v)`$ is pathwise nonincreasing in $`m`$. Consequently $`B_{m,n} / \sqrt{p(m, n)/n}`$ is pathwise nondecreasing in $`m`$ wherever $`p(m, n) > 0`$.
- **Exact one-step mean (eq. (5)):** Splitting $`j \mathbb{E}_j[g_j(w_j)] = \sum_{w=1}^j r_w e_w / [d_w(d_w + 1)]`$ across all leaves of early cores ($`e_w = 1`$, $`r_w = 1/d_x`$, $`g_j(w) = 1/2`$, including early leaves since $`m \ge 2`$ forces every early leaf's unique neighbor to be an early core), early cores ($`e_x = d_x(m)`$ fixed for $`j \ge m`$), and crossing late cores ($`e_y = 1`$, while non-crossing late vertices have $`e = 0`$) recovers eq. (5) identically.
- **Exact telescoping charges (line 50):** An early vertex $`x \le m`$ has fixed $`e_x = d_x(m)`$, so its unweighted charges after any time $`t_0`$ telescope to $`\le d_x(m)/d_x(t_0) \le 1`$. A young vertex $`y > m`$ has $`e_y = 0`$ unless its parent $`x \le m`$ is early, in which case its promotion (degree $`1 \to 2`$) contributes $`1/2`$ and all subsequent gains telescope to $`1/2 - 1/d_y(\infty) \le 1/2`$ (total $`\le 1`$).

### 2.2 Uniform Square Bound per Old Host (Eq. (A)) and Post-Observation Host Protection

1. **Cell host counts & polynomial-in-$`n`$ concentration (`NEXT-RECIPIENT-AVERAGING.md` line 15):** Under `SMALL-POWERS.md`, cell host counts $`1 + B_x`$ are conditionally independent given $`\mathcal{F}_{T_m}`$ with uniformly bounded moments of every order. Expanding centered $`2p`$-th moments gives $`\mathbb{E}[(\sum_{i=1}^{C_m} \tilde{Y}_i)^{2p} \mid \mathcal{F}_{T_m}] \le K_p C_m^p`$, so Markov yields failure probability $`\le K_p C_m^{-p}`$. Under the polynomial cutoff hypothesis $`m \ge n^\rho`$ ($`\rho > 0`$) and the size gate $`C_m \ge m^{\mu/2} \ge n^{\rho \mu/2}`$, choosing $`p = \lceil 2 a / (\rho \mu) \rceil`$ gives $`O_a(n^{-a})`$ for arbitrary $`a > 0`$.
2. **Repeated-trial spoiled-branch envelope (lines 17–23):** Every old host $`x`$ has $`d_x \ge 2`$ from slot activation onward (zero prior crossings exist at slot activation). Whenever $`L_x^{\mathrm{young}} \ge 1`$ with no active trial, start $`J`$ at $`1`$ with rates $`(j, j/2)`$, ordering young leaves first ($`j/d_x \le \min(1, j/2)`$ when $`J = L_x^{\mathrm{young}}`$) and trial-born crossing endpoints last (disjoint from young leaves). Crucially, Astra bounds the number of **distinct spoiled branches**: every crossing created during a finite (extinct) trial is written off in full via the host's ring count in that trial, waiting periods ($`L_x^{\mathrm{young}} = 0`$) create zero crossings, and during the first surviving trial only crossings before stabilization $`\tau_{\mathrm{stab}}`$ (written off via rings up to $`\tau_{\mathrm{stab}}`$) or trial-born crossings hit by a ring in the upper strip $`\min(1, N(s-) e^{-s/4})`$ after $`\tau_{\mathrm{stab}}`$ can be spoiled.
   - *Independent verification of the exact linear BD$`(j, j/2)`$ law:* Solving $`\partial_s f = (f - 1)(f - 1/2)`$ with $`f(0) = z`$ gives $`p_0(s) = \frac{1 - e^{-s/2}}{2 - e^{-s/2}}`$ and $`J(s) \mid \{J(s) > 0\} \sim \mathrm{Geom}(a)`$ with $`a = \frac{e^{-s/2}}{2 - e^{-s/2}}`$ exactly, giving $`\mathbb{P}(s < \tau_{\mathrm{ext}} < \infty) = \frac{a}{2} \le a`$. Cauchy–Schwarz over integer intervals gives finite moments of all orders for finite-trial ring counts.
   - *Stabilization & upper strip:* $`\mathbb{P}(0 < J(i) < 4 e^{i/4}) \le 4 e^{-i/4}`$; by gambler's ruin on the $`(2/3, 1/3)`$ embedded walk, $`\mathbb{P}_j(\inf J \le \lfloor j/2 \rfloor) = 2^{-\lceil j/2 \rceil} \le 2^{-\lfloor j/2 \rfloor}`$; since $`2 > e^{1/4}`$, $`J(s) \ge e^{s/4}`$ for all $`s \ge I_*`$ with exponential tail on $`I_*`$. On $`[i, i+1)`$, upper-strip hits $`X_i \le N([i, i+1))`$ have $`\mathbb{P}(X_i \ge 1) \le \mathbb{E}[X_i] \le (i+1) e^{-i/4}`$, so Cauchy–Schwarz gives $`\|X_i\|_{L^p} = O_p((i+1)^{1/(2p)} e^{-i/(8p)})`$, summable in $`i`$ for every $`p \ge 1`$. Summing over $`\mathrm{Geom}(1/2)`$ failed trials + one surviving trial yields i.i.d. slot spoiled-branch envelopes with all moments, transferring the purely physical bound $`\sum_x (\text{spoiled}_x)^2 \le K_3 C_m`$ to the calendar virgin coupling.
3. **Marked calendar Poisson envelope & eq. (7) (lines 25–29):** For $`Y(u) = \int \mathbf{1}_{u+a > L-1} e^{L-u-a} Q(da) \in [0, e]`$, $`\int_0^\infty \mathbb{E}[Y(u)]\,du \le e`$, and convexity $`e^{z Y} - 1 \le \frac{e^{ze}-1}{e} Y`$ gives $`\mathbb{E}[\exp(z \sum Y(u_k))] \le \exp(e^{ze} - 1)`$ uniformly in $`L`$. Combined with own-label charges ($`\le e`$) and spoiled branches, this proves (A): $`\sum_x W_x^2 \le K_0 C_m`$ outside probability $`O_a(n^{-a})`$ on $`\{C_m \ge m^{\mu/2}\}`$. Exponential maximal bounds on $`\sum_{i \ge m} (E_i - 1)/i`$ give $`\mathbb{P}(\sup_{j \ge m}|T_j - T_m - \log(j/m)| > \varepsilon) = O_\varepsilon(e^{-c_\varepsilon m}) = O_a(n^{-a})`$ since $`m \ge n^\rho`$; Cauchy–Schwarz then yields eq. (7).
4. **Fixed-time classification & post-$`U`$ host protection (lines 31–50):**
   - Noncircular quantifier hierarchy $`\eta \to (D, D') \to J \to K \to U \to (R, M_R) \to n \to \infty`$; active ($`\ell \ge K(1+a) \ge 2`$) and passive (all neighbors core with $`\ge K(1+d_{\mathrm{tag}})`$ leaves) are disjoint.
   - At fixed $`U`$, using `operator.md` §1's oriented old-edge proposals (rate $`1/c_x^0`$), simple-path reversal bounds the annulus alarm by $`(D + D D') U^{R-1}/(R-1)!`$. Truncating to balls $`B_{2R+4}(o)`$ of max degree $`\le M_R`$ creates a dependency graph of max degree $`\le 1 + M_R^{2R+2}`$, colorable into finitely many conditionally independent classes where Hoeffding gives stretched-exponential tails in $`C_m`$ around bounded cylinder expectations converging via Claim 10's local LLN.
   - *Active protection (line 46):* Private young-leaf lower stock from $`\ell \ge K(1+a)`$ + Poisson bound $`N(s) \le \sqrt{K}(s+1)`$ + upper strip $`\min(1, 4(a + \sqrt{K}(s+1)) e^{-s}/\ell)`$ integrates to $`\le 4/K + 8/\sqrt{K}`$ (using $`\int_0^\infty (s+1)e^{-s}\,ds = 2`$), independent across active hosts. Successful actives promote a young leaf at every post-$`U`$ ring (hence belong to $`H_m`$ and never return to any crossing branch).
   - *Passive protection (line 48):* Every donor core $`z`$ adjacent to a retained passive tag has at most $`D_1 = \max(1, D, D')`$ passive targets (early donors have initial degree $`\le \max(D, D')`$; late donors have $`\le 1`$ early neighbor). Independent donor lower-stock protections bound failed passive tags by $`\eta C_m`$; successful passives stay bounded degree ($`\notin H_m`$), produce zero future crossings, and their existing late neighbors have degrees $`\ge K(1 + d_{\mathrm{tag}})`$, leaving remaining total reciprocal drop $`\le 1/K`$ per tag (weighted sum $`\le e(D+1)C_m / K`$).
   - Consequently $`||A| - H_m| = O(\eta C_m)`$, so $`|A|/C_m = h + O(\eta)`$; exceptional hosts (even selected using the future!) contribute $`O(\sqrt{\eta}\, C_m)`$ by (A) + Cauchy–Schwarz; active old-label charges contribute $`\le e D (D+1) C_m / K \to 0`$.

### 2.3 Age Averaging and the Discounted Observation-Time Split

1. **Pre-$`U`$ existing crossing interiors (lines 56–58):** Every crossing branch interior of size $`s`$ gains new vertices at rate $`\le s + 1`$, so $`s + 1`$ is dominated simultaneously across disjoint branches by independent Yule processes started at $`3`$ (via thinning). At bounded $`U`$, all but $`\eta C_m`$ of the $`O_U(C_m)`$ pre-$`U`$ interiors have size $`\le A_0`$. From this finite set of size-$`\le A_0`$ rooted states, fresh independent inert comparisons started at $`U`$ match successful active branches forever; their remaining reciprocal drop after age $`L - 1 - U \to \infty`$ has max expectation $`\varepsilon_{\max}(n) \to 0`$ over the finite state set, and Hoeffding over all size-$`\le A_0`$ interiors bounds the total contribution uniformly over the future-selected successful-active subset.
2. **Post-$`U`$ marked calendar averaging (line 60):** On the $`\mathcal{F}_{T_m + U}`$-measurable active list $`A`$, a rate-1 Poisson ring at $`u \ge U`$ with inert virgin reward $`Q(ds)`$ contributes $`Y_v(u) = \int_{[0,\infty)} \mathbf{1}_{u+s > L+v} e^{L-u-s} Q(ds)`$. By Tonelli,

```math
\int_U^\infty \mathbb{E}[Y_v(u)]\,du = e^{-v} \mathbb{E}\left[\int_{[0, \infty)} \min\big(1, e^{L+v-U-s}\big) Q(ds)\right] \xrightarrow{L - U \to \infty} e^{-v}\mathbb{E}[Q([0, \infty))] = e^{-v}\left[\frac{1}{2} + \left(\frac{1}{2} - b\right)\right]
```

by dominated convergence (no convergence rate needed). Uniform exponential moments give a.s. concentration across $`x \in A`$; removing failed actives costs $`O(\sqrt{\eta}\, C_m)`$ on both physical and ideal sides; sending $`\eta \downarrow 0`$ proves **eq. (6)** of `NEXT-RECIPIENT.md`.
3. **Discounted split before/after observation (`NEXT-RECIPIENT.md` line 66):** Setting $`v = 0`$, a post-promotion gain at age $`s > 0`$ on a branch born at $`u`$ occurs after $`L`$ ($`u + s > L`$) with weight $`e^{L - u - s}`$:
   - Promoted **after** $`L`$ ($`u > L`$): $`\int_L^\infty e^{L - u - s}\,du = e^{-s} \implies \alpha_{\mathrm{new}} = \int_{(0, \infty)} e^{-s} \Gamma(ds) = \mathbb{E}_J[\int_0^\infty e^{-s}(-d(1/D_{\mathrm{root}}(s)))]`$.
   - Promoted **before** $`L`$ ($`U \le u \le L`$, $`u > L - s`$): $`\int_{L - s}^L e^{L - u - s}\,du = 1 - e^{-s} \implies \alpha_{\mathrm{old}} = \int_{(0, \infty)} (1 - e^{-s}) \Gamma(ds) = (1/2 - b) - \alpha_{\mathrm{new}}`$.
   - Replacing $`L`$ by $`T_n \in [L - \varepsilon, L + \varepsilon]`$ only flips branches promoted in $`[L - \varepsilon, L + \varepsilon]`$, which have $`\mathrm{Poisson}(2\varepsilon)`$ count per productive host and cost $`O(\varepsilon C_m)`$.

### 2.4 Almost-Sure Conditional $`L^1`$ Lift (`NEXT-RECIPIENT.md` Lines 60–64)

Let $`Z_n := n B_{m,n} / C_m \ge 0`$ for deterministic $`m \ge n^\rho, n/m \to \infty`$:
1. Sandwiching $`T_n - T_m \in [L - \varepsilon, L + \varepsilon]`$ and $`n/(j+1) \in e^{\pm \varepsilon} e^{L - s_j}`$ in eq. (6) (the first jump $`j = n`$ occurs strictly after $`T_n`$) gives $`Z_n \to \kappa h`$ a.s.
2. Taking $`a = 6`$ in eq. (7) gives $`\mathbb{P}(Z_n > K, C_m \ge m^{\mu/2}) = O(n^{-6})`$. By Markov on $`\mathbf{1}_{G_n}\mathbb{P}(Z_n > K \mid \mathcal{F}_n)`$ ($`G_n = \{C_m \ge m^{\mu/2}\} \in \mathcal{F}_n`$) and Borel–Cantelli ($`\sum n^{4-6} < \infty`$), $`\mathbb{P}(Z_n > K \mid \mathcal{F}_n) \le n^{-4}`$ eventually a.s.
3. From eq. (4), $`B_{m,n} = (R_m - p) - M_{m,n}`$ with $`|R_m - p| \le 1`$ and $`\mathbb{E}[M_{m,n}^2 \mid \mathcal{F}_n] \le \frac{1}{4}\sum_{j \ge n} (j+1)^{-2} \le \frac{1}{4n} \le 1`$, so $`\mathbb{E}[B_{m,n}^2 \mid \mathcal{F}_n] \le 4`$ and $`\mathbb{E}[Z_n^2 \mid \mathcal{F}_n] \le 4 n^2`$ (for $`C_m \ge 1`$).
4. Conditional Cauchy–Schwarz bounds the tail $`\mathbb{E}[Z_n \mathbf{1}_{Z_n > K} \mid \mathcal{F}_n] \le (2n)(n^{-2}) = 2/n \to 0`$ a.s.
5. For the capped variables $`\tilde{Z}_n = \min(Z_n, K) \in [0, K]`$, the decreasing $`\mathcal{F}_\infty`$-measurable majorant $`D_N^* = \sup_{r \ge N} |\tilde{Z}_r - \kappa h| \downarrow 0`$ and Lévy's upward martingale convergence theorem give $`\mathbb{E}[|\tilde{Z}_n - \kappa h| \mid \mathcal{F}_n] \to 0`$ a.s. Together these prove $`\mathbb{E}[|Z_n - \kappa h| \mid \mathcal{F}_n] \to 0`$ a.s., hence eq. (2) (and identically for each nonnegative component of $`B_{m,n}`$).

### 2.5 Distinct Threshold $`m^{1+\mu}/n`$, Shift $`\kappa h \sqrt{\lambda X / c}`$, and Monotone Interpolation

1. **Exact distinct scaling & shift:** Since $`C_m \sim X m^\mu`$ and $`p(m, n) \sim q(m, n) \sim R_m \sim c X m^{\mu - 1}`$ a.s.,

```math
\frac{C_m / n}{\sqrt{p(m, n)/n}} = \frac{C_m}{\sqrt{n \, p(m, n)}} \sim \frac{X m^\mu}{\sqrt{n c X m^{\mu - 1}}} = \sqrt{\frac{X}{c}} \sqrt{\frac{m^{1+\mu}}{n}}.
```

Thus the critical boundary for the next-recipient center $`p(m, n)`$ is **$`m^{1+\mu} / n \to \lambda \in (0, \infty)`$** (without the $`(\log n)^2`$ factor of the observed excess-degree center $`q_n`$!). At $`m^{1+\mu}/n \to \lambda`$, $`m \asymp n^{1/(1+\mu)}`$ satisfies $`m \ge n^{1/2}`$ (since $`\mu < 1`$) and $`m = o(n)`$ (since $`\mu > 1/2`$), so eq. (2)'s conditional $`L^1`$ limit combines with the compensated-noise CLT $`M_{m,n}/\sqrt{p(m,n)/n} \mid \mathcal{F}_n \Rightarrow N(0, 1)`$ (`ADDENDUM-CUTOFF-THRESHOLD.md` eq. (6), valid for all deterministic $`m \to \infty, m = o(n)`$) via conditional Slutsky (no independence between $`M_{m,n}`$ and $`B_{m,n}`$) to give Kolmogorov convergence to $`N(\kappa h \sqrt{\lambda X/c}, 1)`$ and matching conditional mean convergence. In particular, the older excess-degree critical regime $`m^{1+\mu}(\log n)^2 / n \to \lambda`$ satisfies $`m^{1+\mu}/n \to 0`$, so $`p(m, n)`$ remains strictly centered $`N(0, 1)`$ there (matching our independent corollary in `REVIEW-CRITICAL-EXCESS-DRIFT-16f0c12.md`).
2. **Monotone interpolation to every deterministic diverging sublinear $`m = o(n)`$ (`NEXT-RECIPIENT.md` line 74):**
   - Fix the countable grid $`m_a(n) = \lfloor a n^{1/(1+\mu)} \rfloor`$ ($`a \in \mathbb{Q}_{>0}`$). Simultaneously on one probability-$`1`$ event $`\Omega_0`$, for all $`a \in \mathbb{Q}_{>0}`$, $`B_{m_a, n}/\sqrt{p(m_a, n)/n} \to \kappa h \sqrt{a^{1+\mu} X / c}`$ in conditional $`L^1(\mathcal{F}_n)`$, while for the arbitrary given $`m(n) \to \infty, m(n) = o(n)`$, $`p(m, n) > 0`$ eventually and $`M_{m, n}/\sqrt{p(m, n)/n} \mid \mathcal{F}_n \Rightarrow N(0, 1)`$ uniformly in Kolmogorov distance.
   - *Subcritical ($`m^{1+\mu}/n \to 0`$):* Eventually $`m \le m_a`$, so pathwise monotonicity gives $`0 \le \mathbb{E}[B_{m,n}/\sqrt{p(m,n)/n} \mid \mathcal{F}_n] \le \mathbb{E}[B_{m_a, n}/\sqrt{p(m_a, n)/n} \mid \mathcal{F}_n] \to \kappa h \sqrt{a^{1+\mu} X/c}`$. Letting $`a \downarrow 0`$ yields vanishing conditional $`L^1`$ bias (even for ultra-slow $`m`$, e.g. $`m = \log\log n`$, where $`m \ge n^\rho`$ fails!) and hence $`N(0, 1)`$.
   - *Supercritical ($`m^{1+\mu}/n \to \infty`$):* Eventually $`m \ge m_a`$, so $`B_{m,n}/\sqrt{p(m,n)/n} \ge B_{m_a, n}/\sqrt{p(m_a, n)/n}`$ escapes to $`+\infty`$ in conditional probability as $`a \uparrow \infty`$; tightness of $`M_{m,n}/\sqrt{p(m,n)/n}`$ yields escape of $`(R_m - p)/\sqrt{p/n}`$ to $`+\infty`$.
   - *Positive limsup ($`\limsup m^{1+\mu}/n > a^{1+\mu} > 0`$):* Along a deterministic subsequence $`m(n_k) \ge m_a(n_k)`$, conditional probability of $`B_{m,n_k}/\sqrt{p/n_k} \ge d(\omega) := \frac{1}{2}\kappa h \sqrt{a^{1+\mu} X/c} > 0`$ tends to $`1`$; uniform Kolmogorov convergence of $`M_{m,n_k}/\sqrt{p/n_k}`$ evaluates cleanly at the pathwise point $`-d(\omega)`$, giving $`\limsup \mathbb{P}((R_m - p)/\sqrt{p/n} \le 0 \mid \mathcal{F}_{n_k}) \le \Phi(-d(\omega)) < 1/2`$, excluding $`N(0, 1)`$.

---

## 3. Independent Zero-Shared-Code Computational Evidence (Including Independent Certification of the Discounted Split)

Without reading `NEXT-RECIPIENT-EVIDENCE.md` or running `check_next_recipient.py`, we implemented an independent exact rational script (`audit_next_recipient.py`, reviewer scratch script, not exported) checking:
1. Exact positive-jump identity (3), monotonicity in $`m`$, and 1-step mean identity (5) across all $`5,912`$ edge-seeded trees ($`n = 3..8`$), $`40,318`$ `(tree, m)` pairs, $`34,404`$ direct post-jump trees ($`n \le 7`$), and $`34,406`$ monotonicity comparisons (plus $`31,784`$ positive early-leaf witnesses proving omitting early leaves from $`L_x`$ in (5) fails).
2. Exact Kolmogorov forward equations, geometric conditional law $`a = e^{-s/2}/(2 - e^{-s/2})`$, extinction tail $`a/2 \le a`$, and gambler's ruin bound $`2^{-\lceil j/2 \rceil} \le 2^{-\lfloor j/2 \rfloor}`$ for the linear BD$`(j, j/2)`$ trial process.
3. **Independent zero-shared-code exact rational Bellman enclosure for $`\alpha_{\mathrm{new}}`$ and $`\alpha_{\mathrm{old}}`$:** Using integration by parts $`\alpha_{\mathrm{new}}(T) = \frac{1}{d_0} - \mathbb{E}_T[\int_0^\infty e^{-u} \frac{1}{D_{\mathrm{root}}(u)}\,du]`$, dominating $`D_{\mathrm{root}}(u)`$ by a pure-birth process $`Y(u)`$ of rate $`k - 1`$ (since the exterior slot never sends) gives the exact upper bound $`\alpha_{\mathrm{new}}(T) \le (d_0 - 1)(\sum_{k=d_0}^\infty k^{-2} - d_0^{-1})`$ (bounded rationally by splitting at $`K = d_0 + 60`$ and telescoping $`\sum_{k > K} k^{-2} < 2/(2K + 1)`$), while waiting for the first ring among the $`\ell`$ initial root leaves before the root's first ring yields the exact lower bound $`\alpha_{\mathrm{new}}(T) \ge \frac{\ell}{(\ell + 2) d_0 (d_0 + 1)}`$. Propagating through the exact discounted Bellman equation $`\alpha_{\mathrm{new}}(T) = \frac{1}{\lambda_T + 1}\sum_{T'} \mathrm{rate}(T \to T')[\mathbf{1}_{\mathrm{root\ gain}} \frac{1}{d_0(d_0+1)} + \alpha_{\mathrm{new}}(T')]`$ on our independent Dyck-encoded size-10 tree graph ($`1,204`$ states) **strictly certifies both of Astra's line-66 intervals from inside**:
   - $`\alpha_{\mathrm{new}} \in [0.1196058, 0.1305172] \subsetneq [0.1191026, 0.1310137]`$;
   - $`\alpha_{\mathrm{old}} = (1/2 - b) - \alpha_{\mathrm{new}} \in [0.2528328, 0.3151925] \subsetneq [0.2523362, 0.3156958]`$.

Command run (exit 0):
```text
python3 -uB audit_next_recipient.py  # reviewer scratch run
```
Actual output:
```text
[1] Positive-jump & eq.(5) identities PASSED: trees=5912, (tree,m)-pairs=40318, direct_jumps=34404, mono_checks=34406, planted_defect_early_leaf_witnesses=31784
[2] Linear BD(j, j/2) master equation, geometric law, and gambler's ruin PASSED
[3] Independent discounted-split enclosure PASSED (cap=10, states=1204):
    b             in [0.0652016, 0.1166501]
    alpha_new     in [0.1196058, 0.1305172]   (exact: [0.1196058664, 0.1305171074])
    alpha_old(sub)in [0.2528328, 0.3151925]
    alpha_old(dir)in [0.2466875, 0.3168936]   (exact: [0.2466875637, 0.3168935482])
    1-step hand identity verified: alpha_new(J) >= 1/15 = 0.066667
```

---

## 4. Inherited-Dependency Verdict & Provenance Ledger

**Inherited-Dependency Verdict: `CONDITIONAL`** — All deductions inside frozen `NEXT-RECIPIENT.md` (`0cd50c23...`) and `NEXT-RECIPIENT-AVERAGING.md` (`f5bed442...`) hold cleanly, conditional on the following exact inputs:

### 4.1 Round 1 / Round 2 Public Baseline Inputs (`round2/LEDGER.md`, pinned base `c250a1f2fdd53ced91b7b485e94ac8a16106b9c7`)
1. **Astra 6 (Tier A):** Core growth scaling $`C_n / n^\mu \to X \in (0, \infty)`$ a.s., with $`\mu \in (1/2, 1)`$.
2. **Astra 7 / Theorem 9 (Tier A):** Unit total limiting share $`\sum_{v \ge 1} Z_v = 1`$ and $`\mathrm{TV}(p_j, Z) \to 0`$ a.s.
3. **Astra 4 (`operator.md` §1, Tier B for cutoff error):** Old-cell oriented proposal construction (rate $`1/c_x^0`$) and simple-path reversal bound $`c_w^0 U^\ell / \ell!`$.
4. **Astra 10 (Tier B):** Empirical uniform-core local limit $`\pi_n \to \pi`$ a.s. on locally finite rooted core trees (bounded cylinder LLN and local tightness of balls).
5. **Physical degree dichotomy:** `Neon 2(i)` dichotomy (**Tier B**) for finite trees ($`D_v(\infty) = \infty \iff Z_v > 0`$ simultaneously across labels) and `late/README.md` §3 (`Astra 18b`, **Tier C** in ledger, independently cold-audited in Neon's `PROOF.md` §1.2) for the infinite locally finite environment dichotomy.

### 4.2 Previously Cross-Family Reviewed Q4 Notes (`9be4767` — `late-mass-neon/REVIEW.md` (prior snapshot review, not exported), `HOLDS`)
6. `PROOF.md`: $`R_m \sim c X m^{\mu - 1}`$ a.s. (also independently proved in Neon's `PROOF.md`).
7. `FINITE-TREE.md`: $`p(m, n) \sim q(m, n) \sim R_m \sim c X m^{\mu - 1}`$ a.s. for all diverging $`m = o(n)`$ (and $`\mathrm{TV}(p_j, Z) = O_\omega(j^{\mu-1})`$ used inside eq. (6) below).
8. `SMALL-POWERS.md`: Conditional-independent per-cell old-leaf promotion bounds $`B_x`$ with all moments uniform in $`L_x`$, and autonomous lower leaf-stock bound $`(\sqrt{e}/2)^\ell + 2 e^{-\ell/8}`$.
9. `POSITIVE-LABELS.md`: $`H_m / m^\mu \to X M_0 / \mu`$ a.s., hence $`H_m / C_m \to h = M_0 / \mu`$ a.s.

### 4.3 Previously Neon-Audited Post-Review Continuation Inputs ([`REVIEW-CRITICAL-EXCESS-DRIFT.md`](REVIEW-CRITICAL-EXCESS-DRIFT.md), `HOLDS`)
10. `ADDENDUM-CUTOFF-THRESHOLD.md` **eq. (6) only**: Compensated-noise conditional CLT $`\mathrm{law}(M_{m,n}/\sqrt{q_n/n} \mid \mathcal{F}_n) \to N(0, 1)`$ a.s. for every deterministic $`m \to \infty, m = o(n)`$ (independently audited line-by-line in our `16f0c12` review).
11. `CRITICAL-EXCESS-EVIDENCE.md`: Inert-virgin constant enclosure $`0.8833499 \le \kappa = 1 - b \le 0.9347984`$ (independently rebuilt with zero shared code).

### 4.4 Confirmed Non-Dependencies
- Does **not** depend on `CRITICAL-EXCESS-DRIFT.md`'s log-squared asymptotics, `ADDENDUM-DIRICHLET.md`, `ADDENDUM-GAUSSIAN.md`, `ADDENDUM-FIRST-HIT.md`, `ADDENDUM-CUTOFF-THRESHOLD.md` eqs. (5)/(7)/(9)/(10)/(11), `NEXT-RECIPIENT-EVIDENCE.md`, or `Astra 20b`.
