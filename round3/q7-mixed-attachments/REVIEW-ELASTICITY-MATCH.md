> **Public checkpoint export note:** Edited export of `mixed-neon/REVIEW-DB439AC.md` (frozen original SHA-256 `b913c377025ed2c39ba4a719277401a09145ab252df62a92e080bcedc7df6482`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

# Independent Audit of Astra Q7 Revision `db439ac` (`ELASTICITY-MATCH.md` & `PHYSICAL-REDIRECT.md`)

- **Reviewer:** Mixed Attachments — Neon lead (Gemini 4 Neon, XHigh)
- **Reviewed read-only board inputs (revision `db439acd6c534072f66cd23602a2748ec51ce672`), all four SHA-256 pins verified byte-for-byte:**
  - [`ELASTICITY-MATCH.md`](ELASTICITY-MATCH.md): `6ae3a3be86530a53831c922ae3e6cfaa194e4004d7d256968b12c24155aabe26`
  - `PHYSICAL-REDIRECT.md` (not exported): `e712714c2602ac1fe4187f609993493e64120aac3d601e01a09b1692c06b778a`
  - [`check_elasticity_match.py`](check_elasticity_match.py): `91c563cac8bcd58ee103bc0552043a9ea25c8e8d817666ce64b1a02f0ea76e75`
  - `physical_redirect.py` (not exported): `89d9d323b12a750c0232fccc7de8379909b0ed0200cca3f98d6f52055be44bf3`
- **Independent zero-shared-code verifier:** `audit_db439ac.py` (reviewer worktree script, not exported) — SHA-256 `a2a0c5d910d8632a49368f7b8bff2aaea5e74cc92a84f968d9bde4d7686fc80b` (commit `7d2dbd7` on `research/ft3-mixed-neon`). **No Astra Python file (`check_elasticity_match.py`, `physical_redirect.py`, or `review_high_p_bounds.py`) was imported or executed.** No human review.

---

## Executive Verdicts (Two Separate Parts)

1. **Audit 1 (`ELASTICITY-MATCH.md` `6ae3a3be...` & `check_elasticity_match.py` `91c563ca...`): HOLDS.**
   - The age-weighted stochastic-ordering stopping inequality $`p/q < C B(a)`$ ($`B(a) = \pi(1+a)/\sin(\pi a)`$) — including the early-interval $`[0, t_1]`$ piecewise majorant $`w(s)`$ and Chebyshev covariance bound — holds rigorously.
   - Its uniform transfer to **every fixed $`p \in [927/1000, 1)`$** ($`\beta(p) = 0`$, $`m(0) = 1`$, advancing Neon's earlier threshold $`187/200 = 0.935`$ down to $`927/1000 = 0.927`$) holds rigorously across all four inherited conditions of Neon's `HIGH-P.md` Proposition 1 (with the narrowed fresh-immigrant coupling scoped in `REVIEW-HIGH-P.md` / `HIGH-P-DELTA.md`) and all endpoints, independently re-certified in pure `Fraction` arithmetic (`B1_lo = 17.1074 > 17`, `L_lo = 5.6053 > 5.60`, `0 < T1_hi = 1.995e-3`, `0 < T2_hi = 8.145e-2`, `C_hi = 2.750065 < 2.751`, `C * B_hi = 12.671301 < 12.68 < 927/73`).
   - Also holds: exact deterministic scalar threshold $`k = q/y \implies \int k = 1`$ iff $`p/q \ge R`$ (and $`Y(t) \sim c_p t^{1 + a_-}`$ for $`p/q > R`$), and the finite-age physical ratio necessity $`\sup_{t > 0} y(t)k(t)/q > 5/R`$ if $`\beta(5/6) > 0`$.
   - **Scope boundary:** Establishes **no** monotonicity of $`\beta(p)`$, **no** unique transition point in $`(R/(1+R), 927/1000)`$, and **no** sign at $`p = 5/6`$.
2. **Audit 2 (`PHYSICAL-REDIRECT.md` `e712714c...` & `physical_redirect.py` `89d9d323...`): HOLDS.**
   - The exact positive generator decomposition $`(\mathcal{G}_{\mathrm{phys}} - L) F(x) = \sum_i (k^{-1} + r_i)\nabla_i F(x)`$, the backward identity Eq. (1), the resolvent identity Eq. (3), the tagged-donor two-sided bound Eq. (4) ($`\sigma = 1/28`$), and the localization of the $`3/5`$-moment criterion to the single grandchild occupation band $`[11 u, 164 u]`$ (Eq. (5) and pre-birth statistic $`\sum_{k \ge 2} \mathbb{E} Z_k`$) all hold rigorously.
   - The finite-tree grandchild lower bound Eq. (6) ($`\int_0^1 G_i(h, x)\,dh \ge c/[d_i(d_i+1)k(k+1)]`$ with $`c > 27/1000`$, yielding $`\ge 0.000081`$ by original age $`39`$ across all $`3{,}047`$ rooted shapes through size $`11`$, using zero lower bound on $`m(0)`$), the root-to-child upper bound Eq. (7) ($`P < 0.00960`$), and the unconditional ODE deficit bound Eq. (8) ($`M(t) \ge 1 - (1 + 4t)^{-1/4}`$, giving $`M(1024) > 7/8`$ at original age $`6{,}144`$, $`B_\downarrow(\infty) > 0.08476`$, and $`B_\downarrow((250^4-1)/4) > 0.08076`$ at original age $`3(250^4-1)/2`$, combining Eq. (7) with the previously audited $`b_0 \le 254913347315347/281474976710656`$) all hold rigorously.
   - **Scope boundary:** Makes **no** claim of $`\beta(5/6) > 0`$ (or $`\beta(5/6) = 0`$) and **no** claim of a unique transition; $`0 \le \beta(5/6) < 1/1000`$ remains unresolved.

---

## Audit 1 Reconstruction: Age-Weighted Stopping Inequality & Uniform Transfer to $`[927/1000, 1)`$

### 1.1 Age-Weighted Stopping Inequality and Early-Interval Covariance Step
In original-clock time, with $`Y = y - 1`$, $`v = Y' = p + k * v`$, and $`H(t) = t v(t) - (1 + a) Y(t)`$ ($`0 < a < 1`$), at a putative first positive zero $`T \in (0, \infty)`$:
1. Substituting $`r = T - s`$ into $`H(T) = -a p T + (k * H)(T) + \int_0^T s k(s) v(T - s)\,ds`$ (with $`(k * H)(T) \le 0`$) gives

```math
a p T \le \int_0^T (T - r) k(T - r)\,dY(r).
```

2. Since $`(Y(r)/r^{1+a})' = H(r)/r^{2+a} \le 0`$ on $`(0, T]`$, the CDF $`F_Y(r) = Y(r)/Y(T)`$ satisfies $`F_Y(r) \ge (r/T)^{1+a} = F_a(r)`$ on $`[0, T]`$, i.e. $`\mu_Y := dY(r)/Y(T) \preceq_{\mathrm{st}} \mu_a := (1+a)(r/T)^a dr/T`$. Hence every nondecreasing nonnegative test function $`\phi(r)`$ satisfies $`\int_0^T \phi\,d\mu_Y \le \int_0^T \phi\,d\mu_a`$ (monotone truncation handles integrable singularities at $`r = T`$).
   *(Counterexample confirming nondecreasingness is essential: for $`Y(r)/Y(T) = r/T \ge (r/T)^{1+a}`$, the decreasing test $`\psi(r) = r^{-a}`$ gives $`\int_0^1 u^{-a}\,du = 1/(1-a) > 1+a = \int_0^1 u^{-a}(1+a)u^a\,du`$, reversing the inequality.)*
3. Under the physical stopped bounds inherited from `HIGH-P.md` Proposition 1 Steps 1–2 ($`T > t_1(p)`$, $`Y(t_1) < p t_1/(1-i)`$, $`s k(s) \le q/p`$ on $`(0, t_1]`$, $`k(s) \le C_2 q/y(s)`$ on $`[t_1, T]`$), set $`C := \max(\frac{1-a}{1-i}, C_2)`$, $`A := q T^{1+a}/Y(T)`$, and

```math
w(s) := \begin{cases} C t_1^{-a}/(1-a), & 0 < s \le t_1,\\ C s^{-a}, & t_1 < s \le T. \end{cases}
```

- **Monotonicity:** Since $`\frac{1}{1-a} > 1`$, $`w(s)`$ drops downward at $`s = t_1`$ and decreases strictly on $`(t_1, T]`$, so $`r \mapsto w(T - r)`$ is nondecreasing on $`[0, T)`$.
- **Pointwise majorization:** On $`(0, t_1]`$, $`Y(T)/T^{1+a} \le Y(t_1)/t_1^{1+a} < p/[(1-i)t_1^a] \implies s k(s) \le q/p < A t_1^{-a}/(1-i) \le A w(s)`$; on $`(t_1, T]`$, $`y(s) = 1 + Y(s) > Y(T)(s/T)^{1+a}`$ gives strictly $`s k(s) < C A s^{-a} = A w(s)`$.
- **Chebyshev covariance on $`[0, t_1]`$:** Under $`S \sim \mathrm{Unif}(0, t_1)`$, both $`S^{-a}`$ and $`(T - S)^a`$ are strictly decreasing, so $`\mathrm{Cov}(S^{-a}, (T - S)^a) > 0`$. Since $`\mathbb{E}[S^{-a}] = t_1^{-a}/(1 - a)`$,

```math
\int_0^{t_1} w(s)(T - s)^a\,ds = C t_1\,\mathbb{E}[S^{-a}]\,\mathbb{E}[(T - S)^a] < C t_1\,\mathbb{E}[S^{-a}(T - S)^a] = \int_0^{t_1} C s^{-a}(T - s)^a\,ds.
```

- Adding $`\int_{t_1}^T C s^{-a}(T - s)^a\,ds`$ and evaluating $`\int_0^T s^{-a}(T - s)^a\,ds = T \Gamma(1-a)\Gamma(1+a) = T \pi a/\sin(\pi a)`$ cancels $`a q T`$ and yields $`p/q < C B(a)`$.

### 1.2 Verification of All Four Inherited Conditions, Endpoints, and Uniformity on $`[927/1000, 1)`$
With $`(p_0, a, i, \delta, \eta, \theta) = (927/1000, 11/25, 159/200, 1/4, 11/100, 13/50)`$ ($`\rho_0 = 73/927`$):
1. **Inherited Condition (1) ($`a \ge \rho_0/(1 - i)`$):** $`\rho_0/(1 - i) = 14600/38007 \approx 0.38414 < 0.44 = a`$; since $`\rho(p)`$ strictly decreases on $`[p_0, 1)`$, $`T > t_1(p)`$ holds for all $`p \in [p_0, 1)`$.
2. **Inherited Condition (2) ($`X > 100`$, $`(1-\delta)X \ge 1`$, $`i p_0 > q_0^2`$):** $`i p_0 = 0.736965 > 0.005329 = q_0^2`$, so $`\frac{d}{dp}\log t_1(p) = i/(1-p)^2 - 1/p > 0`$ and both $`t_1(p)`$ and $`X(p) = p t_1(p)`$ strictly increase on $`[p_0, 1)`$, with $`X(p_0) \approx 10{,}888 > 100`$ and $`(1 - \delta)X(p_0) > 8{,}166 \ge 1`$.
3. **Inherited Condition (3) ($`B_1 = c_\eta \delta X > 1 + a`$):** Certified $`B_1 \ge 17.1074 > 1.44 = 1 + a`$. Because $`B_1(p) \ge B_1 > 1 + a > 1`$, both $`r \mapsto r^{1+a} e^{-B_1(p) r}`$ ($`r \ge 1`$) and $`x \mapsto x e^{-c_\eta \delta x}`$ ($`x \ge X(p_0)`$) decrease strictly, transferring Tail 1 uniformly to $`[p_0, 1)`$.
4. **Inherited Condition (4) ($`1 - a - \rho_0/\sqrt{2} > 0`$ and $`L(1 - a - \rho_0/\sqrt{2}) > a`$):** Here $`L = B_2 f_0 = \frac{\theta^2 e^{-2b} X(1-i)(1-\delta)}{2(1+a)^2[2 + (i + 17\rho_0/10)^2]}`$ matches `HIGH-P.md` line 205 identically since $`X = p_0 t_1(p_0)`$. Certified $`L \ge 5.6053 > 5.60`$, $`1 - a - \rho_0/\sqrt{2} > 0.504 > 0`$, and $`L(1 - a - \rho_0/\sqrt{2}) > 2.82 > 0.44 = a`$. Since $`\rho(p)`$ decreases while $`X(p)`$ and $`e^{-2b(p)}`$ increase on $`[p_0, 1)`$, $`L(p)`$ increases and Tail 2 transfers uniformly to $`[p_0, 1)`$.
5. **Endpoint handling in $`C`$ and $`B(a)`$:** Both $`T_1 \le 1.995 \times 10^{-3}`$ and $`T_2 \le 8.145 \times 10^{-2}`$ are certified strictly positive using outward rational bounds (`x_hi / exp_lo(b1_lo)` and `1 / exp_lo(L_lo)`); $`C_{\mathrm{main}} = \exp(b + (1+a)(-\log(1-\delta)))/(1-\theta)`$ is outward-rounded up; $`C_1 = (1-a)/(1-i) = 112/41 \approx 2.7317 < C_2 \le 2.750065 < 2.751`$; and $`\pi < 22/7`$, $`\sin(\pi a) = \cos(\pi(1/2 - a)) \ge 1 - ((1/2 - a)\cdot 22/7)^2/2`$ give $`C B(a) \le 12.671301 < 12.68 < 927/73 = p_0/q_0 \le p/q`$ for every $`p \in [927/1000, 1)`$.

---

## Audit 2 Reconstruction: Grandchild Accounting & Remaining Occupation Condition at $`p = 5/6`$

### 2.1 Exact Generator Decomposition and Localization to $`[11 u, 164 u]`$
- **Generator identity:** In redirection-clock age $`s = t/6`$, on functions $`F(x)`$ of the children's degree list $`x = (d_1, \dots, d_{k-1})`$, all transitions except root-to-child sends (rate $`1/k`$) and grandchild-to-child sends (rate $`r_i`$) match the direct-aging reference generator $`L`$ (verified algebraically in exact `Fraction` arithmetic on all $`37`$ rooted shapes of size $`1 \le |T| \le 6`$ in `audit_db439ac.py`). Localization by population size (dominated by rate-$`6`$ Yule) yields Eqs. (1) and (3).
- **Exponential tail Eq. (2):** Splitting $`F_h(\varnothing) = \int_0^1 \exp(-5\int_0^h H(v, u)\,dv)\,du`$ at $`u_0 = \frac{1}{10} e^{-h/102}`$ gives $`u e^{v/100} \ge \frac{1}{10}`$ for all $`v \ge \frac{50}{51} h`$ (since $`-\frac{1}{102} + \frac{50}{5100} = 0`$), so $`5\int_0^h H \ge 5\cdot\frac{1}{10}\cdot\frac{h}{51} = \frac{h}{102}`$ and $`F_h(\varnothing) \le \frac{11}{10} e^{-h/102}`$.
- **Two-sided tagged-donor response Eq. (4) and Eq. (5) localization:**
  - At $`\sigma = 1/28`$, $`1 + 5/\sigma = 141 < 143.6895 = \sum_{j=0}^9 5^j/j! < e^5`$, so $`\widehat{g}_d(1/28) = (\log 141)/5 < 1`$, giving $`G_i(h, x) \le C_\sigma e^{h/28}/[d_i(d_i+1)k^2]`$; the independent rate-$`5`$ immigrant forest gives $`G_i(h, x) \le \frac{11}{10} h e^{-h/102}/[d_i(d_i+1)]`$.
  - Using $`r_i \le d_i - 1 \implies \sum_i (k^{-1} + r_i)/[d_i(d_i+1)] \le \sum_i d_i^{-1} \le a(x)`$, $`\mathbb{E}[T_{k+1}^{3/5}] \le (k/5)^{3/5}`$, and the root-jump compensator $`\mathbb{E}\int a(x_s) b(s, k_s)\,ds = \sum_{k \ge 1} \mathbb{E} b(T_{k+1}, k)`$:
    - $`h \le 11 u`$: summand exponent $`-2 + \frac{11}{28} + \frac{3}{5} = -1 - \frac{1}{140} < -1`$ (converges);
    - $`h \ge 164 u`$: summand exponent $`-\frac{164}{102} + \frac{3}{5} = -1 - \frac{2}{255} < -1`$ (converges);
    - Root sends on $`[11 u, 164 u]`$: $`5 \sum_i \nabla_i V \le a(x) V(x) \le C a(x)(1+u)/k`$ gives summand $`k^{-7/5}\log k`$ (converges).
  - By `ZERO-SCALE.md` §4, $`\int_0^\infty t^{3/5} g(t)\,dt < \infty \iff \beta(5/6) > 0`$, so Eq. (5) — equivalently $`\sum_{k \ge 2} \mathbb{E} Z_k < \infty`$ via the predictable compensator of $`N_{\mathrm{root}}`$ — is necessary and sufficient for $`\beta(5/6) > 0`$.

### 2.2 Finite-Tree Grandchild Bound (Eq. (6)), Root-Send Cap (Eq. (7)), and ODE Deficit (Eq. (8))
- **Eq. (6) & age-$`39`$ occupation:** Using $`1 - z \le H(v, z) \le 1 - z/[1 + (e^v - 1)(1-z)]`$ (whose upper integral is $`L = \log(1 + (e^h - 1)(1-z))`$), coupling $`C_{d+1} - C_d`$ with $`e^{-A_{d+1}} - e^{-A_d} \ge e^{-L/d}(A_d - A_{d+1})`$ and keeping $`\{N(v) = 0\}`$ yields Eq. (88–89), and Jensen under $`k(k+1) z^{k-1}(1-z)`$ (mean of $`1-z`$ is $`2/(k+2)`$, $`2a(x)/(k+2) \le 3`$ for $`k \ge 2`$) yields Eq. (6). Independently verified on a 160-rectangle pure-`Fraction` grid ($`c \ge 0.027697 > 27/1000`$) and an independent base-$`10^{18}`$ rooted-shape enumeration through size $`11`$ (confirming all $`3{,}047`$ shapes, internal rate sum $`6j - 1/k`$, half-unit holding retention finishing by rescaled time $`11/2 + 1 = 13/2 \leftrightarrow`$ original age $`39`$, and grandchild-only extra return $`\ge 0.00008330 > 0.000081`$, using zero lower bound on $`m(0)`$).
- **Eq. (7) ($`P < 0.00960`$):** Adding one fresh child multiplies the integrand of $`F_h(x)`$ by $`1 - H(h, z)`$, so $`5[V(x) - V(x \oplus 1)] \le \int_0^1 z^{k-1}\int_0^\infty 5H e^{-5\int H}\,dh\,dz = 1/k`$. Inserting into $`L V = -1/k`$ gives $`\sum_i \nabla_i V \le (a(x) - 5)/(25 k) \le a(x)\frac{k-1}{25 k(k+4)}`$, and the root compensator yields $`P \le \frac{1}{25}\sum_{k=1}^\infty \frac{k-1}{k^2(k+4)} = \frac{1}{25}(\frac{125}{192} - \frac{\pi^2}{24}) \approx 0.0095923 < 0.00960`$.
- **Eq. (8) ($`M(t) \ge 1 - (1+4t)^{-1/4}`$):** Before the first unit crossing ($`d(s) = 1 - M(s) > 0`$), $`U(s) \le \sum M(s)^j = 1/d(s)`$, so $`A(t) = 1 + 5\int_0^t d(s)^{-1}\,ds`$ obeys $`-d' \ge 1/A`$ and $`A' = 5/d \implies d'/d \le -A'/(5A) \implies d \le A^{-1/5} \implies A' \ge 5 A^{1/5} \implies (A^{4/5})' \ge 4 \implies A \ge (1+4t)^{5/4} \implies M(t) \ge 1 - (1+4t)^{-1/4}`$. Hence at original age $`6{,}144`$ ($`t = 1024`$), $`M(1024) \ge 1 - 4097^{-1/4} > 7/8`$; combining with Eq. (7) and the preserved audited bound $`b_0 \le 254913347315347/281474976710656`$ yields $`B_\downarrow(\infty) > 0.08476`$ and $`B_\downarrow((250^4-1)/4) > 0.08076`$ at original age $`3(250^4-1)/2`$.

---

## Independent Verification Replay Evidence

Command:
```bash
python3 -B audit_db439ac.py  # reviewer worktree run (see check_elasticity_match.py for exported scalar certificate)
```
Exact output (exit code 0):
```text
[PASS] All 4 read-only db439ac board file SHA-256 pins verified.
[PASS] Independent Audit 1 (ELASTICITY-MATCH.md at p0=927/1000) verified:
       B1_lo=17.1074 > 17, L_lo=5.6053 > 5.60, 0 < T1_hi=1.995e-03, 0 < T2_hi=8.145e-02, C_hi=2.750065 < 2.751, C*B_hi=12.671301 < 12.68 < 12.698630.
[PASS] Independent Audit 2 (PHYSICAL-REDIRECT.md at p=5/6) verified:
       e^5 partial=143.6895 > 141; margins 139/140 < 1 & 257/255 > 1; exact generator identity on 37/37 shapes;
       grad_c(160)=0.027697 > 0.027; 3047 shapes (size<=11); gc_age39_lo=0.00008330 > 0.000081;
       P_hi=0.0095923 < 0.00960; all-age gc > 0.084774 (>0.08476); finite gc > 0.080774 (>0.08076); M(1024) > 7/8.
ALL INDEPENDENT DB439AC AUDIT CHECKS PASSED.
```
