> **Public checkpoint export note:** Edited export of `forest-neon/REVIEW-WORST-HUB-UPPER.md` (frozen original SHA-256 `bd0b62d1c6c67896b9b84baf86254e47dfcc4bb27dd9ccba3acf1c861f4ce120`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

# Independent Cross-Audit of Astra's `WORST-HUB-UPPER.md`: Sharp Worst Single-Parent Log-Volume KS Limit $`D^*(B)`$ on $`\alpha \in (\frac{\mu}{1+\mu}, 1)`$

**Reviewer:** Neon lead (Gemini 4 Neon, XHigh).  
**Freeze-Before-Read Verification:** Both read-only target files were copied verbatim into an immutable local scratch snapshot **before reading**, with exact SHA-256 match (zero mismatch):

| Target file | Expected SHA-256 | Frozen-byte SHA-256 | Status |
| :--- | :--- | :--- | :--- |
| `WORST-HUB-UPPER.md` | `9b2b7b0a02f37757857ce80fe7d565685ef0653ed9dc90850fc3b8dff9c82cd4` | `9b2b7b0a02f37757857ce80fe7d565685ef0653ed9dc90850fc3b8dff9c82cd4` | **EXACT MATCH (zero mismatch)** |
| `worst_hub_upper_check.py` | `ae77f2b79c5547428c6906c859f70c5d87002f31c7b177003fb381a2c7e99fb5` | `ae77f2b79c5547428c6906c859f70c5d87002f31c7b177003fb381a2c7e99fb5` | **EXACT MATCH (zero mismatch)** |

**Scope & Separation:** Preserves all prior local commits (`0b36e06e`, `7c3c5ca5`, `d363c702`, `91a20378`, `50b0a875`, `3800e612`, `5444dca2`), rival files, and `paper/` completely untouched. An earlier replay of Astra's script verified finite stock identities and a numerical variational grid; below we independently audit the full analytic theorem. Automated cross-family AI audit, not human peer review.

---

## 1. Scoped Verdicts (Distinguishing Fixed Observed Forest from Joint Dynamic Laws)

| Component / Claim in `WORST-HUB-UPPER.md` | Verdict | Exact Scope & Dependencies |
| :--- | :--- | :--- |
| 1. Block-uniform degree sandwiches & eligible old-child observable $`R_u`$ (§1) | **HOLDS** | **Fixed-observed-forest view** for deterministic $`N(n) = n^{\alpha + o(1)}`$, fixed $`\alpha \in (\frac{\mu}{1+\mu}, 1)`$. Conditional on Tier-A Yule clock + Tier-B Claim 5a §3.1 tagged-leaf bounds & $`O(t)`$ max core arity. Margin $`w = \log D_f + \sqrt{m} = o(L)`$ handles arbitrary subpower $`n^{o(1)}`$ variation within integer blocks $`I_m = [e^m, e^{m+1})`$. |
| 2. Eligible old-child exponential census at $`\theta = 3/4`$ (§2, eqs. (2) and (2a)) | **HOLDS** | Exact birth–death lower stock chain $`J_v`$ ($`j \to j+1`$ at rate $`j`$, $`j \to j-1`$ at rate $`j/(j+1)`$) gives $`\mathbb P_1(\lim e^{-v} J_v > 0) \ge 1/(e - 1) > 1/2`$, hence $`\limsup q_m < 1/2`$ and $`\frac{1}{2}(e^{3/4} - 1) < 9/16 < \mu`$. Marginal per-parent coupling + stopping-time compensation requires zero inter-parent independence. |
| 3. Old-measurable gate, conditional Laplace transform, & unnormalized physical-to-ideal transfer (§3, eqs. (3)–(4)) | **HOLDS** | **Fixed-observed-forest view.** Gate $`G_{\mathrm{old}}(m) = \{A_m \le e^{(B + \delta_d)L}\}`$ is strictly $`\mathcal F_s`$-measurable (zero future conditioning on $`\tau_n`$ or $`H(\tau_n)`$). Unnormalized linear score $`M_u^{\mathrm{fav}} - p M_u \ge 0`$ handles empty rows ($`M_u = 0 \implies 0 \ge 0`$) and endpoints $`x \in \{0, 1\}`$ cleanly; $`\max_{\mathcal H_m}(2 S_u + B_u + 1) = o(L)`$ a.s. |
| 4. Clipped two-bin variational optimization (§4, eq. (5): $`\inf_{0 \le z \le 1-d} I_{3/4}(z, z+d) = 1 - e^{-J(d)}`$) | **HOLDS** | Exact second derivative $`H''_\theta(y) = e^\theta e^{-\theta/y}\theta[\theta(d+y) - 2dy]/y^4 > 0`$ on $`y \in (0, 1 - d]`$ (since $`2dy/(d+y) \le 1/2 < 3/4`$), non-clipping as $`y \downarrow 0`$ ($`y t_* \to 0`$), and $`z = 0`$ bound $`\frac{1-d}{d} J(d) \le \frac{2\log v}{v+1} < 2/e < 3/4`$ prove clipping never raises the KS optimum above $`1 - e^{-J(d)}`$. Raw $`\theta = \log(1+\mu)`$ genuinely fails (verified at $`B = 3/8, d = 1/2`$). |
| 5. Uniform block summability & exact match with Neon's lower bound Proposition D1 (eq. (1): $`\max_{u \in H(\tau_n)} d_{\mathrm{KS}}^{(1)}(u) \to D^*(B)`$ a.s.) | **HOLDS** | **Fixed-observed-forest view.** For every fixed $`\alpha \in (\frac{\mu}{1+\mu}, 1)`$ and deterministic $`N(n) = n^{\alpha + o(1)}`$, Astra's upper bound $`\limsup \le D^*(B)`$ matches Neon's Proposition D1 (`REVIEW-ASTRA-WORST-PARENT-KS.md` §3) lower bound $`\liminf \ge D^*(B)`$ identically, closing our previously open single-hub supercritical regime `(R1)` at fixed $`\alpha`$. |
| 6. Two-parent extremal constant $`\max_{i \ne j} d_{\mathrm{KS}}^{(2)}(i, j)`$, oscillating exponents $`\alpha_- < \alpha_+`$, sub-polynomial $`\alpha_- = 0`$ | **OUTSIDE SCOPE** | Explicitly excluded from eq. (1) by Astra line 14 (two-parent constant requires a separate two-coordinate $`(x, y)`$ variational proof; oscillating $`\alpha_- < \alpha_+`$ breaks $`L_+ - L_- = o(L)`$ within blocks; $`\alpha_- = 0`$ lies outside Engine A). |
| 7. Process-level joint dynamic law / whole-forest joint independence across hubs or epochs | **UNSUPPORTED** | Eq. (1) controls the quenched maximum over empirical parent rows inside the fixed observed tree $`T_n`$. No inter-hub or dynamic whole-forest independence law is proved or implied. |

---

## 2. Independent Reconstruction of the Four Audited Mechanisms

### 2.1 Block Setup and Eligible Old-Child Exponential Census (§§1–2, eqs. (2) & (2a))
On integer block $`I_m := \{n : e^m \le n < e^{m+1}\}`$ with $`L := \alpha m`$ and $`L_\pm := \max/\min_{n \in I_m} \log N(n) = L + o(L)`$, localizing on $`\{Y \in (e^{-C}, e^C)\}`$ brackets $`\tau_n \in [t_-, t_+] = [m - C, m + 1 + C]`$ for all $`n \in I_m`$. Setting deterministic $`s := t_+ - L_- + 4\log m`$, $`h_\pm := t_\pm - s = L + o(L)`$, $`D_f := 2 e^{L_+ - h_-} \ge m^4`$, and $`w := \log D_f + \sqrt{m}`$ yields $`w = o(L)`$, $`s = (1 - \alpha)m + o(m) > w`$, $`e^{L_- - h_+} = m^4`$, and $`D_f e^{-w} = e^{-\sqrt{m}} \to 0`$.
- **Degree sandwiches:** For $`D_u(s) < m^4/2`$, upper Yule domination gives $`\mathbb P(D_u(t_+) > e^{L_-} \mid \mathcal F_s) \le 2^{m^4/2} e^{-m^4/2} = e^{-((1-\log 2)/2)m^4}`$ (same super-exponential tail for post-$`s`$ births), so $`\bigcup_{n \in I_m} H(\tau_n) \subseteq \mathcal H_m := \{u \in K(T_s) : D_u(s) \ge m^4/2\}`$ eventually a.s. For $`D_u(s) > D_f`$ on $`\max c(s) \le K_*(m + O_C(1))`$, tagged-leaf lower domination (`Claim 5a` §3.1) gives $`e^{-h_-} D_u(t_-) \ge \frac{3}{4} D - \frac{1}{8} D = \frac{5}{8} D > \frac{1}{2} D_f = e^{L_+ - h_-}`$ outside probability $`O(e^{-m^4/128})`$, so every $`y \in K(T_s)`$ with $`D_y(s) > D_f`$ is deleted at **all** $`n \in I_m`$ and contributes **zero** to any parent's row. Thus pre-$`s`$ tops in $`u`$'s row are bounded by $`R_u := \#\{y \sim u \text{ core at } s : D_y(s) \le D_f\}`$.
- **Lower stock chain $`J`$ and $`\theta = 3/4`$ census:** A newly promoted child $`y`$ starts with $`L_y = 1`$ leaf and $`c_y \ge 1`$ core neighbors forever; at stock $`j \ge 1`$ its leaf gain rate is $`\ge j`$ and promotion loss rate is $`j/(j + c_y) \le j/(j + 1)`$. Pathwise order coupling via $`y`$'s leaf-birth clocks and $`y`$'s own initiating clock (putting leaf choices in $`[0, L_y/D_y)`$) dominates the birth–death chain $`J_v`$ ($`j \to j+1`$ at rate $`j`$, $`j \to j-1`$ at rate $`j/(j+1)`$, absorbed at $`0`$). Solving $`i(h(i+1) - h(i)) + \frac{i}{i+1}(h(i-1) - h(i)) = 0`$ gives $`h_k(1) = 1/\sum_{j=1}^k 1/j! \to 1/(e - 1)`$. Since $`e^{-v} J_v`$ is a nonnegative supermartingale with mean $`\ge k - 1`$ and second moment $`\le k(k+1)`$ started from $`k`$, Paley–Zygmund + strong Markov at $`T_k`$ gives $`\mathbb P_1(W_J > 0) \ge h_k(1)\frac{(k-1)^2}{k(k+1)} \to \frac{1}{e - 1} > \frac{1}{2}`$ (already $`> 1/2`$ at $`k = 20`$!). Since $`D_f e^{-w} = e^{-\sqrt{m}} \to 0`$, $`\limsup q_m := \limsup \mathbb P_1(J_w \le D_f) \le 1 - \frac{1}{e - 1} < \frac{1}{2}`$. For any single parent $`u`$ promoted at time $`v \le s`$, preassigning independent private $`J`$-drivers to $`u`$'s rate-$`1`$ unthinned proposal clock (disjoint newborn roots and leaf clocks, needed only marginally per parent!) gives: any accepted child promoted by $`s - w`$ with $`J_w > D_f`$ has $`D_y(s) \ge D_y(\sigma + w) \ge J_w > D_f`$ by degree monotonicity, hence does NOT enter $`R_u`$. Thinning $`(v, s - w]`$ by $`q_m`$ and keeping all proposals on $`(\max(v, s - w), s]`$ gives

```math
\mathbb E[e^{\theta R_u} \mid \mathcal F_v] \;\le\; e^\theta \exp\!\Big((e^\theta - 1)\big[q_m(s - w - v)_+ + \min(w, s - v)\big]\Big).
```

At $`\theta = 3/4`$ ($`q_m \le 1/2`$ eventually), $`q_m(e^{3/4} - 1) \le \frac{1}{2}(e^{3/4} - 1) < \frac{9}{16} < 0.565713 < \mu`$. Compensating promotion stopping times $`\tau_k \le s`$ against the deterministic function $`f_m(v)`$ (no predictability of $`R_u`$ required!) with $`\mathbb E C(T_v) \le C e^{\mu v}`$ yields $`\mathbb E A_m \le C_\theta \exp\{\mu s + (e^\theta - 1)w\}`$, and Markov + Borel–Cantelli over $`m`$ yields $`\limsup L^{-1}\log A_m \le B`$ a.s. (eq. (2), and analogously eq. (2a) conditional on $`\mathcal F_{v_0}`$).

### 2.2 Old-Measurable Gate, Conditional Transform, and Unnormalized Transfer (§3, eqs. (3)–(4))
Conditional on $`\mathcal F_s`$, for fixed threshold $`x \in [0, 1]`$ and candidate $`u \in \mathcal H_m`$, the monotone endpoint tests on $`(s, t_-]`$ ($`V_{\mathrm{ideal}}(t_- - \sigma) \le e^{x L_+}`$ upward; $`V_{\mathrm{ideal}}(t_+ - \sigma) > e^{x L_-}`$ downward) yield independent Poissons $`(X_u, Y_u)`$ with deterministic means $`L(z + o(1))`$ and $`L(1 - z + o(1))`$ ($`z = x`$ upward, $`z = 1 - x`$ downward) — including at endpoints $`x = 0`$ (upward: $`V(T^0_a) \ge 2 \implies \mathbb E X_u = 0`$) and $`x = 1`$ (downward: $`h_+ = L_- - 4\log m \implies \mathbb E X_u \le 2 L / m^4 = o(L)`$!).
For any $`t \ge 0`$ with $`t(1 - p) \le \theta`$, exponential Markov summed over the $`\mathcal F_s`$-measurable set $`\mathcal H_m`$ using $`\sum_{u \in \mathcal H_m} e^{t(1-p)R_u} \le A_m`$ gives eq. (4) with **zero future conditioning**. On the degree sandwiches, every deleted $`u \in H(\tau_n)`$ satisfies the exact unnormalized inequality (avoiding division by $`M_u`$ and covering empty rows $`M_u = 0`$, where downward deviation still obeys $`0 \ge 0`$!):

```math
M_u^{\mathrm{fav}} - p M_u \;\le\; (1 - p)(X_u + R_u) - p Y_u + (1 + p) S_u + (1 - p) B_u + p \;\le\; (1 - p)(X_u + R_u) - p Y_u + 2 S_u + B_u + 1.
```

On the $`\mathcal F_s`$-measurable gate $`\{A_m \le e^{(B + \delta)L}\}`$ (so $`|\mathcal H_m| \le A_m`$), $`S_u \preceq \mathrm{Pois}(O(m^{-2}))`$ and $`B_u \sim \mathrm{Pois}(1 + 2C)`$ have tails $`\exp(-c_\varepsilon L \log L)`$, giving $`\max_{u \in \mathcal H_m}(2 S_u + B_u + 1) = o(L)`$ a.s.

### 2.3 Clipped Two-Bin Variational Optimization & Endpoints (§4, eq. (5))
For $`d \in (0, 1)`$, $`z \in [0, 1 - d]`$, $`p = z + d`$, and $`\theta = 3/4`$:
- Unclipped $`\inf_{t \ge 0} [z e^{(1-p)t} + (1-z)e^{-pt}] = e^{-D_{\mathrm{KL}}(\mathrm{Ber}(p)\|\mathrm{Ber}(z))} \le e^{-J(d)} =: H_*`$ (at $`p = 1`$, constraint $`0 \cdot t \le \theta`$ is vacuous and infimum is $`1 - d`$).
- Where clipping $`(1-p)t_* > \theta`$ binds, put $`y = 1 - p \in (0, 1 - d]`$: $`H_\theta(y) = e^\theta[1 - d - y + (d + y)e^{-\theta/y}]`$ has exact second derivative

```math
H_\theta''(y) \;=\; \frac{\theta e^\theta e^{-\theta/y}}{y^4}\big[\theta(d + y) - 2 d y\big] \;>\; 0
```

since $`\frac{2dy}{d+y} \le \frac{d+y}{2} \le \frac{1}{2} < \frac{3}{4} = \theta`$. Hence $`H_\theta(y)`$ is **strictly convex** on every clipping interval, attaining its maximum at an endpoint.
- Near $`y \downarrow 0`$ ($`p \uparrow 1`$), $`y t_* = y \log(1/y) + O(y) \to 0 < \theta`$, so clipping **never** binds near $`y = 0`$.
- At interior switching points, $`H_\theta`$ equals the unclipped value $`\le H_*`$.
- At the sole remaining boundary $`y = 1 - d`$ ($`z = 0`$), $`H_\theta(1 - d) = e^{-\theta d / (1 - d)}`$. Testing symmetric $`z = (1 - d)/2`$ ($`v = (1+d)/(1-d) > 1`$) yields

```math
\frac{1 - d}{d}\, J(d) \;\le\; (1 - d)\log v \;=\; \frac{2\log v}{v + 1} \;<\; \frac{2\log v}{v} \;\le\; \frac{2}{e} \;\approx\; 0.7358 \;<\; \frac{3}{4} \;=\; \theta,
```

so $`e^{-\theta d/(1-d)} < H_*`$. This proves eq. (5): $`\inf_{0 \le z \le 1-d} I_{3/4}(z, z+d) = 1 - e^{-J(d)}`$! *(By contrast, raw arity $`\theta = \log(1 + \mu) \approx 0.4484 < \max \frac{2\log v}{v+1} \approx 0.5569`$ genuinely fails eq. (5), e.g. gap $`> 0.05`$ at $`d = 1/2`$.)*

### 2.4 Uniform Block Summability & Exact Match with Neon's Proposition D1
For any rational $`d \in (D^*(B), 1)`$: strict monotonicity of $`p \mapsto D_{\mathrm{KL}}(\mathrm{Ber}(p)\|\mathrm{Ber}(z_*))`$ on $`[z_*, 1]`$ (moving $`p_*`$ inward to $`p' = z_* + \frac{1}{2}(D^*(B) + d) \in (z_*, 1)`$) proves $`1 - e^{-J(d)} > B`$, say $`= B + 4\delta_d`$ ($`\delta_d > 0`$). Choosing a finite tilt $`t(z, d)`$ for each of the finitely many grid points $`x_j = j/q`$ (both directions with $`z + d \le 1`$) makes eq. (4) on the $`\mathcal F_s`$-gate $`\{A_m \le e^{(B+\delta_d)L}\}`$ decay as $`O_q(e^{-\frac{1}{2}\delta_d \alpha m})`$, **geometrically summable over integer blocks $`m`$** for every deterministic $`N(n) = n^{\alpha + o(1)}`$. Grid monotonicity adds $`\le 1/q`$; countable intersection over $`C, K_*, q \to \infty`$, rational $`d \downarrow D^*(B)`$ yields $`\limsup_n \max_{u \in H(\tau_n)} d_{\mathrm{KS}}^{(1)}(u) \le D^*(B)`$ a.s.

**Exact match with Neon's Proposition D1:** In `REVIEW-ASTRA-WORST-PARENT-KS.md` §3 (SHA-256 `8d876c8be0842903e4288143a480ecda3198fd3f3c46b6d0a4211abf96794ea8`), Proposition D1 proved $`\liminf_n \max_{u \in H(\tau_n)} d_{\mathrm{KS}}^{(1)}(u) \ge D^*(B)`$ a.s. for all $`0 < \alpha_- \le \alpha_+ < 1`$ with $`B = \mu(1/\alpha_+ - 1) \in (0, 1)`$, using the identical variational value $`D^*(B)`$ ($`1 - e^{-D_{\mathrm{KL}}} < B \iff D_{\mathrm{KL}} < \log\frac{1}{1-B}`$). Combining with Astra's sharp upper bound completes the **full fixed-exponent phase diagram for worst single-parent log-volume KS** across $`N(n) = n^{\alpha + o(1)}`$, $`N(n)/n \to 0`$, $`\alpha \in (0, 1]`$:
- $`\alpha \in (0, \frac{\mu}{1+\mu}]`$: $`\max_{u \in H(\tau_n)} d_{\mathrm{KS}}^{(1)}(u) \to 1`$ a.s. (Astra §5);
- $`\alpha \in (\frac{\mu}{1+\mu}, 1)`$: $`\max_{u \in H(\tau_n)} d_{\mathrm{KS}}^{(1)}(u) \to D^*(B) \in (B, 1)`$ a.s. (Neon D1 + Astra `WORST-HUB-UPPER.md`);
- $`\alpha = 1`$ ($`N/n \to 0`$): $`\max_{u \in H(\tau_n)} d_{\mathrm{KS}}^{(1)}(u) \to 0`$ a.s. (Neon D2).

*(Remaining open regimes: two-parent extremal constant $`\max_{i \ne j} d_{\mathrm{KS}}^{(2)}(i, j)`$ on $`\alpha \in (\frac{\mu}{1+\mu}, 1)`$, oscillating schedules $`\alpha_- < \alpha_+`$, sub-polynomial $`\alpha_- = 0`$, and whole-forest / process-level joint laws.)*

---

## 3. Independent Adversarial Verification Evidence

Script: `round3/q5-transition/adversarial_review_astra_worst_hub_upper.py` (pure Python standard library, single-threaded, zero shared code with `worst_hub_upper_check.py`).

Reviewer's local adversarial script (`adversarial_review_astra_worst_hub_upper.py`, not exported; see [`worst_hub_upper_check.py`](worst_hub_upper_check.py) for the exported author finite checker):
```sh
python3 -B adversarial_review_astra_worst_hub_upper.py  # reviewer worktree run
```
Actual output (exit code `0`):
```text
PASS [W1] (EXACT RATIONAL) Birth-death hitting recurrence exact (1225 eqs, k=2..50; k=20 rational P_1(W_J>0) > 1/2); (e^(3/4)-1)/2 < 9/16 < mu_lo; raw-arity B=3/8 control verified
PASS [W2] (EXACT FINITE TREE TRAVERSAL) Eligible R_0=1 < c_0=2 and unnormalized score bound verified across 4080 row tests (incl. empty row, bad sends, late births); planted no-S_bad and no-R_0 controls caught
PASS [W3] (EXACT ANALYTIC / NUMERIC) Eq. (5) verified via continuous opt across d=0.01..0.99 (all 99 clipped; max err=2.22e-16; 4950 convexity checks; max 2*log(v)/(v+1)=0.5569 < 2/e < 3/4); planted theta=0.4 curvature & theta=log(1+mu) eq.(5) failure (gap=0.0514) caught
PASS [W4] (NUMERIC) Neon D1 matches Astra upper bound to <1e-10 (a=0.4:D*=0.8496>B=0.849, a=0.5:D*=0.6138>B=0.566, a=0.7:D*=0.3667>B=0.242, a=0.9:D*=0.1795>B=0.063); subpower block scales verified
=== ALL 4 INDEPENDENT ADVERSARIAL CHECKS ON ASTRA WORST-HUB-UPPER PASSED ===
```
*(Honest scope: finite rational, tree-traversal, and continuous-optimization checks with planted negative controls verify the algebraic recurrences, transfer inequalities, convexity/endpoint bounds, and D1 match; the asymptotic almost-sure limit rests on the analytic proof audited in §2. Automated cross-family AI review is not human peer review.)*
