# Neon Rival Attack on Fable Claim 3 (Hub Probability by Label: $`p_1=p_2>p_3>\dots>p_8`$ & Scale $`\Theta(u^{\mu-1})`$)

**VERDICT: HOLDS**

**Scope of verdict:** Fable Claim 3 parts (a), (b), and (c) pinned to `ft2-fable@fd889cf` (`round2/newground/q3-hub-probability/REPORT.md` and `exact.{py,cpp}`). The 62,125-replicate Monte Carlo section is diagnostic (not part of the claim) and was not re-run.

**Attacker:** Neon lead (Gemini 4 Neon, independent of Fable's Claude family and Astra's GPT family). Alongside Astra's GPT attack (`attacks/astra-on-fable-3.md`), this completes three-family verification of Fable Claim 3.

---

## 1. Line-by-Line Audit of All Three Subparts

1. **Part (a) — Finite-state classifiers, symmetry reductions, and exact ordering $`p_1 = p_2 > p_3 > \dots > p_8`$:**
   - *Lower classifier $`B(l, c)`$:* Coupling $`L_v`$ below by the frozen-$`c`$ birth-death chain (rates $`k`$ up, $`k/(k+c)`$ down) yields the exact gambler's-ruin escape probability $`B(l, c) = (\sum_{k=0}^{l-1} a_k)/(\sum_{k=0}^\infty a_k)`$ with $`a_0=1, a_k=1/[(c+1)\cdots(c+k)]`$, and exact geometric tail bound $`\sum_{k=25}^\infty a_k < a_{24}/(c+24)`$. Passage from unbounded leaf stock to $`Z_v > 0`$ via $`Y_1`$ (rates $`k`$ up, $`1`$ down, $`0`$ absorbing) is valid: $`e^{-t}(Y_1(t)-1)`$ is a submartingale (even after absorption at $`0`$, where $`-e^{-t}`$ has positive drift $`+e^{-t}`$), dominated by pure Yule($`k`$) whose limit has second moment $`k(k+1)`$, so Paley–Zygmund gives $`\mathbb{P}(W_1 > 0) \ge (k-1)^2/[k(k+1)] \to 1`$.
   - *Upper classifiers:* For a marked leaf at core $`w`$ ($`c_w \ge 1`$), $`f(k) = 1/k`$ obeys the exact cancellation $`(1 - 1/k) + (k-1)(\frac{1}{k-1} - \frac{1}{k}) = 1`$, so promotion drift $`1/D_w \le 1/(k+1)`$ cancels birth drift $`-1/(k+1)`$. For a star ($`c_w = 0`$), $`f(k) = 1/(k-1)`$ ($`k \ge 3`$) gives promotion drift $`(k^2-3k+3)/[k(k-1)(k-2)] \le 1/(k-1)`$ for all $`k \ge 3`$. For a leafless core ($`l = 0`$), each neighbour $`w`$ has $`D_w \ge k+1`$ and $`f(k) = 1/(k-1)`$ gives selection+promotion drift $`(k^2-3k+4)/[(k-1)(k-2)D_w] \le 1/(k-1)`$ for all $`k \ge 4`$ and $`2/D_w \le 1/2`$ at $`k=3`$. The refined tridiagonal recursion $`F_c(k)`$ is decreasing and $`\le 1/k`$, hence a valid supersolution for any larger birth rate or larger core degree $`c_{\mathrm{true}} \ge c`$.
   - *Symmetry & analytic pairs:* Seed symmetry gives $`\mathrm{law}(Z_1) = \mathrm{law}(Z_2)`$. On $`T_3`$, $`B(2,0) = 2/e > 11/18`$ proves $`p_2 > p_3`$. At $`n = 4`$, $`\mathrm{law}(Z_3) - \mathrm{law}(Z_4) = \frac{1}{6}(I - E)`$, and one step on $`P_4`$ yields $`17/[16(e-1)] > 34/55 > 239/392`$ (at $`K=3`$, $`F_1(2)=22/49, F_1(1)=31/49`$), proving $`p_3 > p_4`$.
   - *Exact rooted-shape enumeration:* Independently coded from scratch in pure Python `fractions.Fraction` through $`n = 9`$ (`286` canonical rooted shapes) without reading Fable's `exact.py`/`exact.cpp`, confirming $`p_1 = p_2`$ identically and strictly positive adjacent lower bounds on all six pairs $`2\text{--}3`$ through $`7\text{--}8`$; ALSO independently compiled and replayed Fable's `exact.cpp` at $`m = 11`$ (all `3,628,800` histories, `1,842` shapes), reproducing the posted lower differences (`7-8: 66855530980254312786365850240 / 16724654043955200000000000000000 = 0.00399742...`). **Verified.**

2. **Part (b) — Universal scale bounds $`\frac{\mathbb{E}[C_{u-1}]}{32(e-1)(u-1)} \le p_u \le \frac{2\mathbb{E}[C_{u-1}]}{u-1}`$ ($`u \ge 4`$):**
   - *Upper bound:* Summing core attachment bounds $`(L_w + \rho_w)/(L_w+1) \le 1 + \rho_w`$ (or $`L_w/L_w = 1`$ for a star) plus leaf attachment rates $`L_w/D_w`$ gives $`\sum_{w \in K} (1 + (c_w + L_w)/D_w) = 2 C_{u-1}`$, hence $`p_u \le 2\mathbb{E}[C_{u-1}]/(u-1)`$.
   - *Lower bound:* Suppressing the newborn's parent $`w`$ until its independent $`\mathrm{Exp}(1)`$ first ring gives $`\mathbb{E}[D_w(s) \mid T_{u-1}, w] \le d e^s`$ ($`d = D_w^{\mathrm{pre}}+1`$), so Jensen yields promotion probability $`\ge \int_0^\infty e^{-s}(d e^s)^{-1} ds = 1/(2d)`$, followed by $`B(1,1) = 1/(e-1)`$. In any core tree on $`C \ge 2`$ vertices, $`J = \{w : c_w \le 2, L_w \ge 1\}`$ contains all $`A \ge B + 2`$ core leaves (so $`|J| > B`$), each contributing $`\ge 1/4`$; the remaining $`H = \{w : c_w = 2, L_w = 0\}`$ has $`|H| \ge C - 2|J|`$ and at least $`4|H| - 2C + 2 \ge 2C - 8|J|`$ directed internal $`H \to H`$ edges, each contributing $`\frac{1/2}{3} = \frac{1}{6}`$; optimizing $`\frac{|J|}{4} + \frac{\max(0, 2C - 8|J|)}{6}`$ at $`|J| = C/4`$ yields $`\ge C/16`$ (and a star $`C=1`$ gives $`\ge 2/3 \ge 1/16`$). Combining gives $`\frac{\mathbb{E}[C_{u-1}]}{32(e-1)(u-1)}`$. **Verified.**

3. **Part (c) — Degree dichotomy $`\{Z_v = 0\} = \{\sup_n D_v(n) < \infty\}`$ a.s.:**
   - Independent continuous-time proof agrees with Neon Claim 2(i) (which proved the same dichotomy in discrete time via a coupled Pólya urn supermartingale): stopping at leaf-present times spaced by $`\ge 1`$ shows any core with leaves at arbitrarily late times hits every $`K`$ a.s., contradicting bounded leaf stock on $`\{Z_v = 0\}`$, so every non-hub core eventually loses all leaves and receives no further attachments. **Verified.**

---

## 2. Independent Verification Script & Replay

Script: `ft2-neon/round2/outputs/attack_fable_3/audit_fable3.py` (log `audit_fable3.log`, stdlib only):
- `[Fable 3 Classifiers] PASS`: exact rational tridiagonal verification of $`F_c(1\dots 60)`$ monotonicity & supersolution for $`c = 1\dots 10`$, exact $`P_4`$ gap $`34/55 > 239/392`$, star & leafless core drifts for $`k \le 200`$, and $`B(l, c)`$ tail truncation bound.
- `[Fable 3 Scale Bounds] PASS`: verified $`\sum_{w \in K} r_w/(D_w+1) \ge C/16`$ (`min ratio = 0.4166`) and upper sum $`\le 2C`$ (`max ratio = 2.0000`) across `2,000` random friend trees + adversarial bare-core paths; caught planted defect claiming $`\ge C/2`$.
- `[Fable 3 Exact n=9] PASS`: independent pure-Python `Fraction` Markov-chain enumeration of all `286` rooted shapes at $`n = 9`$ proves $`p_1 = p_2`$ identically and strictly positive adjacent lower bounds for all six gaps (`2-3: +0.24180`, `3-4: +0.03763`, `4-5: +0.01753`, `5-6: +0.00881`, `6-7: +0.00414`, `7-8: +0.00121`).
- Independently compiled and replayed Fable's `exact.cpp` on `exact-11.in` (all `3,628,800` histories), reproducing the posted size-11 lower bounds cleanly.

**What is NOT proved by Claim 3 (correctly stated in the claim):** Monotonicity of $`p_u`$ for $`u \ge 8`$ and stochastic domination $`Z_u \succeq_{\mathrm{st}} Z_{u+1}`$ remain open.
