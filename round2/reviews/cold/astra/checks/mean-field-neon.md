CORRECT

# Cold Referee Review of `round2/closure/mean-field.md` and `round2/closure/mean_field.py`

- **Targets reviewed**:
  - `round2/closure/mean-field.md` (`closure/mean-field.md`) (45 lines; SHA-256 `27e74f51abe1d05a88076ad1f71a5a2446d1d72f41df1382d9d647813c9ad3eb`)
  - `round2/closure/mean_field.py` (`closure/mean_field.py`) (138 lines; SHA-256 `b233d18a9ca5048428b2a129f7bcd0c52336c031aab733070b42f0ccf4704de1`)
  - Replayed output: `round2/outputs/mean-field.txt` (`outputs/mean-field.txt`) (SHA-256 `f517fe0cfb7cc18355a745790296445a5d2a38d8c111632bd359e0f1b4967207`)
- **Background references checked**: `round2/operator.md` (`operator.md`), `proofs/bracket-motifs.md`.
- **Independent checker**: `round2/checks/mean_field_neon.py` (`mean_field_neon.py`) (stdlib-only, single CPU, exact formal power-series derivation, exact `Fraction` triangular grid vs directed rounding, cross-parameter root verification, direct stopped-CTMC Monte Carlo, pathwise `(L, c)` monotone coupling, and five negative controls).

## Scope and Verdict Summary

This review audits the mathematics and numerical certification of the two **heuristic mean-field closures** ($`\gamma_s(c) = 1 - s`$ and $`\gamma_s(c) = (1 - s)c/2`$) themselves—not claims about the true friend-tree exponent $`\mu`$. Both closure roots lie strictly outside the repository's certified interval $`[0.5650, 0.5665]`$ for $`\mu`$, and the numerical proximity of the per-edge closure root (`0.5688424819...`) to the rounded simulation fit `0.569` in Hartle et al. is treated strictly as a numerical coincidence with no supporting theoretical link.

Within this scope, every analytical step and numerical certification in `mean-field.md` and `mean_field.py` is **CORRECT**: the stationary mean identities, uniqueness via promotion-counting `(L, c)` coupling, inverse-degree polynomial recursion, Lagrange remainder and global error bound, exact rational directed outward rounding, opposite endpoint signs, and reproduced 11-digit brackets all hold without mathematical error or implementation bug.

---

## Point-by-Point Independent Re-Derivation

### 1. Stationary Mean Identities (Lines 3–7)
- **$`\pi(\rho) = \pi(c/D) = 1 - \mu`$**: From `round2/operator.md` (`operator.md`) (lines 11, 75, 77), $`\pi`$ is a unimodular probability eigenmeasure satisfying $`\pi(b_o) = \pi(L_o/D_o) = \mu`$. Since $`D_o = L_o + c_o`$, $`\pi(c/D) = 1 - \pi(b) = 1 - \mu`$. Applying mass transport to the nonnegative transport $`f(G, u, v) = \mathbf 1_{\{u \sim_K v\}}/D_v`$ sends total mass $`\rho_o = \sum_{v \sim_K o} 1/D_v`$ out of $`o`$ and receives $`\sum_{u \sim_K o} 1/D_o = c_o/D_o`$ at $`o`$, giving $`\pi(\rho) = \pi(c/D) = 1 - \mu`$.
- **$`\pi(c) = 2`$**: On every finite nonempty core tree $`K_t`$ with $`C_t \ge 1`$ cores, $`\sum_{u \in K_t} c_u = 2(C_t - 1)`$, so $`m_t(c_o) = 2 F_3(t) - 2`$ and $`\pi_{\delta, T}(c_o) = 2 - \frac{2}{Z}\int_0^\infty w(t)\,dt`$. Splitting at any fixed $`M_0 > 0`$: $`\frac{1}{Z}\int_0^{M_0} w(t)\,dt \le M_0 w(0)/Z \to 0`$ along the construction of `operator.md` Section 2, while $`\frac{1}{Z}\int_{M_0}^\infty w(t)\,dt \le 1/F_3(M_0) \to 0`$ as $`M_0 \to \infty`$, giving $`\pi_{\delta, T}(c_o) \to 2`$. Since $`\mu > 0`$, `operator.md` eq. (5) supplies a uniform exponential moment $`\sup \pi_{\delta, T}(e^{\theta c_o}) < \infty`$ ($`0 < \theta < \log(1 + \mu)`$), yielding uniform integrability of the local continuous functional $`c_o`$ and hence $`\pi(c) = 2`$.
- Consequently, at $`s = \mu`$, both $`\gamma_s(c) = 1 - s`$ and $`\gamma_s(c) = (1 - s)c/2`$ match the stationary mean $`\pi(\rho) = 1 - \mu`$.

### 2. Closure Equation and Uniqueness by Promotion Counting (Line 9)
- In the associated CMJ branching model, each core evolves on $`\{(D, c) : D \ge 2, 1 \le c \le D\}`$ with gain rate $`D - c + \gamma_s(c)`$ and promotion rate $`b = 1 - c/D \in [0, 1]`$. Total jump rate is linear in $`D`$ (non-explosive), and $`0 \le U_s(D, c) = \mathbb E_{(D, c)}[\sum_{k \ge 1} e^{-s \tau_k}] \le 1/s`$. For $`s > 1`$, $`U_s(2, 1) \le 1/s < 1`$, so any root of $`U_s(2, 1) = 1`$ lies in $`(0, 1]`$.
- **Pathwise $`(L, c)`$ coupling for $`0 < s_2 < s_1 \le 1`$**: Couple $`(L_1, c_1)`$ (parameter $`s_1`$, weaker influx) and $`(L_2, c_2)`$ (parameter $`s_2`$, stronger influx) starting from $`(1, 1)`$ via common gain and promotion Poisson point processes on $`(0, \infty)^2`$. Maintain $`L_1 \le L_2`$ and $`c_1 \le c_2`$:
  - *Gain marks*: $`c_1, c_2`$ do not move. At $`L_1 = L_2 = L`$ and $`c_1 \le c_2`$, since $`\gamma_{s_1}(c_1) \le \gamma_{s_2}(c_1) \le \gamma_{s_2}(c_2)`$ (using that $`\gamma_{s_2}(c)`$ is nondecreasing in $`c`$ for $`s_2 \le 1`$ in both models), process 1 cannot gain without process 2 also gaining.
  - *Promotion marks* ($`(L, c) \mapsto (L-1, c+1)`$ at rate $`L/(L+c)`$): At $`L_1 = L_2 = L`$ and $`c_1 \le c_2`$, $`L/(L+c_2) \le L/(L+c_1)`$, so process 2 cannot promote (decrementing $`L_2`$) without process 1 also promoting. At $`c_1 = c_2 = c`$ and $`L_1 \le L_2`$, $`L_1/(L_1+c) \le L_2/(L_2+c)`$, so process 1 cannot promote (incrementing $`c_1`$) without process 2 also promoting.
- Therefore $`c_1(t) \le c_2(t)`$ for all $`t \ge 0`$. Since promotion count up to $`t`$ equals $`c(t) - 1`$, every $`k`$-th promotion time satisfies $`\tau_k^{(2)} \le \tau_k^{(1)}`$. Combined with strictly smaller discount $`e^{-s_1 t} < e^{-s_2 t}`$ on $`\{\tau_1^{(2)} < \infty\}`$, $`U_{s_1}(2, 1) < U_{s_2}(2, 1)`$ strictly. Together with continuity of $`s \mapsto U_s(2, 1)`$ on $`(0, 1]`$ (uniform tail $`\int_T^\infty e^{-st} b(t)\,dt \le e^{-sT}/s`$), opposite endpoint signs certify existence and **uniqueness** of the root.

### 3. Inverse-Degree Polynomial Recursion and Exact Bellman Residual (Lines 21–40)
- Putting $`v = 1/s - U_s \in [0, 1/s]`$ yields $`(s - \mathcal L) v = c/D = c x`$ with $`x = 1/D`$ and $`1/(D+1) = x/(1+x)`$.

• Expanding $`(1 + x)^{-i} = \sum_{k \ge 0} (-1)^k \binom{i+k-1}{i-1} x^k`$ in $`(s - \mathcal L)\sum_{i=1}^m P_i(c) x^i = c x`$ reproduces every coefficient and sign in lines 28–31:

```math
[(s + j) - \Delta] P_j = \mathbf 1_{\{j=1\}} c - c \Delta P_{j-1} + \sum_{i=1}^{j-1} P_i \left[(-1)^{j+1-i}\binom{j}{i-1} + (\gamma_s(c) - c)(-1)^{j-i}\binom{j-1}{i-1}\right].
```

- Inductively $`\deg Q_j \le j`$, so $`\Delta^{j+1} Q_j \equiv 0`$ and $`P_j = \sum_{r=0}^j \Delta^r Q_j / (s+j)^{r+1}`$ is exact with rational coefficients for rational $`s`$.
- **Lagrange remainder bound**: For $`g(x) = (1+x)^{-i}`$ ($`i \ge 1`$) and any $`x \ge 0`$, $`g^{(k)}(\xi)/k! = (-1)^k \binom{i+k-1}{i-1}(1+\xi)^{-(i+k)}`$ has magnitude $`\le \binom{i+k-1}{i-1}`$ on $`\xi \in (0, x)`$. Taking $`k = m + 2 - i`$ (for the $`-x^{-1}`$ piece) and $`k = m + 1 - i`$ (for the $`-(\gamma_s(c)-c)`$ piece) together with the promotion remainder $`+c\Delta P_m(c) x^{m+1}`$ bounds the exact Bellman residual $`|(s-\mathcal L)v_m(1/D, c) - c/D| \le \epsilon`$ uniformly for all $`D \ge N`$ and $`1 \le c < C`$, matching lines 37–39 and `tail_error()` lines 75–84 identically.

### 4. Global Error Bound (Lines 41–45)
- Define $`\widehat v`$ by $`v_m(1/D, c)`$ on $`D \ge N, 1 \le c \le C`$, by $`0`$ on $`2 \le D < N, c = C`$, and by the exact triangular Bellman recursion on $`2 \le D < N, 1 \le c < C`$. Note that at $`(D, C-1)`$ with $`D \ge N`$, promotion moves to $`(D, C)`$ where $`\widehat v = v_m`$ agrees with $`\mathcal L v_m`$, so the interior residual bound $`|(s - \mathcal L)(v - \widehat v)| \le \epsilon \mathbf 1_{\{D \ge N\}}`$ holds everywhere on $`1 \le c < C`$.

• Since $`v - \widehat v`$ is uniformly bounded on $`\{1 \le c \le C\}`$, Dynkin's formula up to $`\tau_C = \inf\{t : c_t = C\}`$ has zero transversality remainder and gives

```math
|v(2, 1) - \widehat v(2, 1)| \le \frac{\epsilon}{s} + \left(\frac{1}{s} + B\right)\mathbb E_{(2,1)}[e^{-s \tau_C}],
```

where $`B = \sum_{i=1}^m |P_i(C)| N^{-i}`$ bounds $`|\widehat v(D_{\tau_C}, C)|`$ on $`D_{\tau_C} \ge N`$ (and $`0 \le B`$ on $`D_{\tau_C} < N`$).
- Because the promotion intensity $`1 - c/D \le 1`$, reaching $`c = C`$ requires $`C - 1`$ promotions whose arrival time $`\tau_C`$ stochastically dominates $`\mathrm{Gamma}(C - 1, 1)`$, yielding $`\mathbb E[e^{-s \tau_C}] \le (1 + s)^{-(C - 1)}`$ and confirming line 43.

### 5. Directed Outward Rounding and Rational Arithmetic in `mean_field.py`
- For rational $`s = p/q \in (0, 1)`$, `polynomials()`, `high_degree_values()`, and `tail_error()` stay purely in exact `Fraction`/`int` arithmetic (zero floats).
- In `enclosed(s, model)`:
  - `math.floor(x * SCALE)` and `-math.floor(-x * SCALE)` invoke exact `Fraction.__floor__` (arbitrary-precision integer division, valid for positive and negative rationals).
  - Clearing denominators by $`d \cdot Q_{\mathrm{gain}} > 0`$ ($`Q_{\mathrm{gain}} = q`$ or $`2q`$) yields exact integer weights `r_num > 0`, `b_num >= 0` (identically `0` when `c == d`), and `denominator > 0`, so `low_num // denominator` and `-(-hi_num // denominator)` strictly bracket `SCALE * \widehat v(d, c)` at every step; `lower[C] = upper[C] = 0` is reset after every degree swap.
  - Converting back to $`U_s = 1/s - v`$ correctly reverses `upper_next[1]` and `lower_next[1]` and subtracts/adds `error` (line 110).

---

## Independent Verification and Negative-Control Evidence

1. **Exact replay of author's script** from repository root (`python3 -B round2/closure/mean_field.py`, exit `0`): stdout matches `round2/outputs/mean-field.txt` (`outputs/mean-field.txt`) byte-for-byte (`diff -u` empty):
```text
constant: closure root in [0.57982073459, 0.57982073462].
  Lower endpoint reward minus one >= 42191/10^15.
  Upper endpoint one minus reward >= 28042/10^15.
per-edge: closure root in [0.56884248194, 0.56884248197].
  Lower endpoint reward minus one >= 29885/10^15.
  Upper endpoint one minus reward >= 41117/10^15.
```

2. **Independent checker run** from repository root (`python3 -B round2/checks/mean_field_neon.py`, exit `0`, single CPU):
```text
=== TEST 1: Finite-Tree Mass Transport & Core-Degree Identities ===
  Verified sum(c)=2(C-1) and mass transport sum(rho)=sum(c/D)=C-sum(b) on 46 trees (n=3..8).
  [PASS] Test 1

=== TEST 2: Independent Formal Series & Exact Bellman Residual ===
  constant: P_1..P_8 exact Fraction match; max |residual|=8.053e-12 <= eps=1.712e-11
  Negative control 1 (constant, sign flip): |residual|=8.311e-05 >> eps=1.712e-11 (caught!)
  per-edge: P_1..P_8 exact Fraction match; max |residual|=1.255e-12 <= eps=3.190e-12
  Negative control 1 (per-edge, sign flip): |residual|=8.376e-05 >> eps=3.190e-12 (caught!)
  [PASS] Test 2

=== TEST 3: Exact Fraction Solver vs Directed Rounding & Negative Controls ===
  constant (N=25,C=10,m=4): exact grid U=0.999586616397, enclosed=[0.95960084, 1.03957239] contains ref=0.999579383779
  Negative controls 2 & 3 (constant): bad ceil, swapped U=1/s-v, and omitted C=3 error [1.104969, 1.104969] != ref 0.999579 all caught!
  per-edge (N=25,C=10,m=4): exact grid U=0.997270845367, enclosed=[0.95826895, 1.03627274] contains ref=0.997265802868
  Negative controls 2 & 3 (per-edge): bad ceil, swapped U=1/s-v, and omitted C=3 error [1.098552, 1.098552] != ref 0.997266 all caught!
  [PASS] Test 3

=== TEST 4: Independent Reproduction & Cross-Parameter Bracket Certification ===
  constant (2400,88,8): [0.57982073459, 0.57982073462] -> lo-1 >= 42191/10^15, 1-hi >= 28042/10^15 (exact match to mean-field.txt!)
  constant (1500,78,7): independently certified strict signs: lo-1 >= 17384/10^15, 1-hi >= 3235/10^15
  per-edge (2400,88,8): [0.56884248194, 0.56884248197] -> lo-1 >= 29885/10^15, 1-hi >= 41117/10^15 (exact match to mean-field.txt!)
  per-edge (1500,78,7): independently certified strict signs: lo-1 >= 24903/10^15, 1-hi >= 36135/10^15
  [PASS] Test 4

=== TEST 5: Direct Stopped-CTMC MC & Pathwise (L, c) Monotone Coupling ===
  Direct MC constant (stopped at D=30): MC=0.80907+/-0.00666 vs exact=0.80795 (z=+0.17)
  Direct MC per-edge (stopped at D=30): MC=0.80135+/-0.00650 vs exact=0.80925 (z=-1.21)
  Coupled (L,c) constant (s1=0.65 > s2=0.50): L-crossings=0, c-crossings=0, U_4(s1)=0.7639 < U_4(s2)=0.9958
  Coupled (L,c) per-edge (s1=0.65 > s2=0.50): L-crossings=0, c-crossings=0, U_4(s1)=0.7334 < U_4(s2)=0.9492
  Negative control 5 (decreasing gamma=(1-s)/c): L-crossings=71 > 0 (caught!)
  [PASS] Test 5

ALL MEAN-FIELD CLOSURE VERIFICATION & NEGATIVE-CONTROL SUITES PASSED.
```

## Minor Non-Blocking Notes
- Line 9 of `mean-field.md` relies on $`\gamma_s(c)`$ being nondecreasing in $`c`$ for $`s \in (0, 1]`$ so that gain rates preserve $`L_1 \le L_2`$ when $`c_1 \le c_2`$ (confirmed by Negative Control 5, where a decreasing $`\gamma(c) = (1-s)/c`$ breaks $`L_1 \le L_2`$). Both closures satisfy this on $`(0, 1]`$, and no root can lie in $`(1, \infty)`$ since $`U_s \le 1/s < 1`$.
- No changes were made to any author file and no commits were created.
