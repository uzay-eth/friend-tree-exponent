CORRECT

# Cold Referee Review: `round2/fringe/finite-kernel.md`

- **Target:** `round2/fringe/finite-kernel.md` (`fringe/finite-kernel.md`) (unmodified; depends on `round2/fringe/renewal.md` (`fringe/renewal.md`), `proofs/existence.md`, `proofs/bracket-motifs.md`)
- **Independent checker:** `round2/checks/finite_kernel_neon.py` (`finite_kernel_neon.py`)
- **Verdict:** **CORRECT** — Both the 19-type determinant characterization $`\det(I_{19} - K(s)) = 0`$ with explicit combinatorial series (2) and the two monotone rational matrix truncations $`K^-_N(s) \le K(s) \le K^+_N(s)`$ with convergence estimates (5)–(7) and uniform inverse slope are rigorously proved. Every diagonal, self-transition, path-count, and endpoint case checks out, and all exact rational matrix checks and planted-error controls pass.

---

## 1. Nineteen Newborn Labels, Clamping Endpoints, and Determinant (1) (§1)

- **Clamping $`m(D) = \min(D, 20)`$ and transfer identity:** Assigning permanent birth mark $`a_v = 1/m(D_{\rm parent}) \in \{1/2, 1/3, \dots, 1/20\}`$ at promotion preserves $`p_u(t) = 1/D_u(t) \le 1/D_u(\tau_v) \le 1/m(D_u(\tau_v)) = a_v`$ for all $`t \ge \tau_v`$. Even when $`D_u(\tau_v) > 20`$ (so $`a_v - p_u > 0`$ already at birth), the killing term $`-a_v f(S_v)`$ cancels $`+a_v f(S_v)`$ from the parent's copy terms $`(a_v - p_u) f(S_v) + p_u f(S_v^+)`$, recovering the exact physical transfer $`p_u[f(S_v^+) - f(S_v)]`$ with actual reciprocal degrees $`p_u`$. Direct births produce $`J_{1/m(D_{\rm root})}`$, matching the physical birth mark at every degree including the endpoint $`D = 20`$ and tail $`D > 20`$.
- **Lyapunov bound on $`D < 20`$, $`D = 20`$, and $`D > 20`$:** With $`d = 3/20`$, $`q = 9/25`$, $`\rho = q + d = 51/100`$, and $`h`$ from `renewal.md` (`fringe/renewal.md`), the newborn path factor is $`(1 + d\,m(D))^{-1}`$.
  - For $`2 \le D < 20`$, $`m(D) = D`$ reproduces the already verified rational inequality (4) of `renewal.md` (tightest pair $`(c, D) = (1, 3)`$, margin $`1/4350`$).
  - For $`D \ge 20`$ (endpoint $`D = 20`$ and all $`D > 20`$), $`m(D) = 20`$ gives factor $`(1 + 20 \cdot 3/20)^{-1} = 1/4 = 0.25`$. Since $`h`$ is nonincreasing and $`1 - c/D \le 1`$, the LHS is at most $`1/4 < 63/250 = 0.252 = q \cdot (7/10) \le q h(c)`$ for every $`c \ge 0`$ (margin $`1/500`$). Notably, $`M = 20`$ is the minimal integer clamp for these constants ($`M = 19`$ gives $`20/77 > 63/250`$ and fails at $`c = 4, D = 135`$).
  - The child-weight identity $`\sum_Y a_Y W(Y) = (a_S + d) W(S) - h(c_{\rm root})`$, $`L^0 W \le q W`$, $`RW \le \rho W - h(c_{\rm root}) \le \rho W`$, and initial-centre estimate $`I_t W \le 1/d + 2 e^{qt}`$ all carry over unchanged.
- **Largest-real-zero specification (1):** Every entry $`K_{ij}(s)`$ ($`i, j \in \{2, \dots, 20\}`$, $`s > \rho`$) is strictly positive and strictly decreasing and continuous in $`s`$, so by finite-dimensional Perron–Frobenius $`r(K(s))`$ is a simple positive eigenvalue, continuous and strictly decreasing on $`(\rho, \infty)`$. The Laplace threshold argument from `renewal.md` §4 (now on 19 labels) identifies $`r(K(s)) = 1`$ uniquely at $`s = \mu \in [113/200, 1] \subset (\rho, \infty)`$. For $`s > \mu`$, every eigenvalue of $`K(s)`$ has modulus $`\le r(K(s)) < 1`$, so $`\det(I_{19} - K(s)) \ne 0`$; at $`s = \mu`$, $`1`$ is an eigenvalue, so $`\det(I_{19} - K(\mu)) = 0`$. Hence $`\mu`$ is the largest real zero of (1) on $`s > 51/100`$. (Subdominant real eigenvalues could cross $`1`$ for $`s \in (\rho, \mu)`$, so "largest real zero" is the exact statement needed.)

---

## 2. Explicit Series (2): Diagonal, Self-Transitions, and Path Count (§2)

- **Diagonal exit rate $`r_0(S)`$:** In any planted marked fringe $`S`$ with $`n(S)`$ vertices (root + leaves + descendants, omitting the parent slot), each leaf and non-root core vertex has total internal arrow rate $`1`$, while the root sends rate $`1/D_{\rm root}`$ to the inert parent slot and $`1 - 1/D_{\rm root}`$ internally. Adding parallel transitions preserves the sum, so $`r_0(S) = \sum_{S'} r^0(S, S') = n(S) - 1/D_{\rm root}`$ identically. Since $`S`$ has a root and at least one child ($`D_{\rm root} \ge 2`$ with 1 parent slot), $`n(S) \ge 2`$.
- **Self-transitions ($`t(S, S) \equiv 0`$):**
  - Every intrinsic transition satisfies $`n(S') = n(S) + 1 > n(S)`$, so $`r^0(S, S) = 0`$.
  - Every $`Y`$-copy has $`n(Y) \le n(S) - 1 < n(S)`$, so $`Y \ne S`$.
  - Every $`Y^+`$-copy has $`n(Y^+) = n(Y) + 1 \le n(S)`$, with equality $`n(Y^+) = n(S)`$ possible *only* when $`S`$ has a bare root ($`L_{\rm root}(S) = 0`$, $`k_{\rm root}(S) = 1`$). Even in that same-size case, $`Y^+`$ has $`L_{\rm root}(Y^+) = L_{\rm root}(Y) + 1 \ge 1 \ne 0 = L_{\rm root}(S)`$, so $`Y^+ \ne S`$.
  - Thus $`t(S, S) = 0`$ identically for every state $`S`$. Line 40's remark (*"with any self-transition added there"*) is harmlessly defensive and makes $`(R_N)_{S, S'} = t(S, S') - (a_S + r_0(S))\delta_{S, S'}`$ valid verbatim without relying on $`t(S, S) = 0`$.
  - Directed cycles of length $`\ge 2`$ *do* occur once $`N \ge 4`$ (e.g., $`J_{1/2} \xrightarrow{r^0} (1/2, 0, (J_{1/2},)) \xrightarrow{r^0} (1/2, 1, (J_{1/2},)) \xrightarrow{\text{copy } Y} J_{1/2}`$), so $`A_N`$ is genuinely non-triangular for $`N \ge 4`$.
- **Path count, multiplicities, and verification of (2):** Parallel intrinsic arrows are summed in $`r^0(S, S')`$, and $`\sum_{Y\text{ core child}}`$ runs with multiplicities (also adding $`p_S`$ and $`a_{Y_2} - p_S`$ whenever $`Y_1^+ = Y_2`$). Every state has finite out-degree in $`t(\cdot, \cdot)`$, and along any length-$`\ell`$ walk from $`J_{1/i}`$ one has $`n(S_u) \le 2 + u \le 2 + \ell`$, so each inner sum at fixed $`\ell`$ has finitely many positive terms. Grading residual lineages by the total number $`\ell`$ of intrinsic jumps and copy births along their unique ancestry (finite on every bounded time interval by nonexplosion), each waiting interval in state $`S`$ lasts $`\mathrm{Exp}(a_S + r_0(S))`$, yielding expected discounted direct reward $`b_S \mathbf 1_{\{m(D_{\rm root})=j\}}/(s + a_S + r_0(S))`$ and expected discounted $`(\ell+1)`$-step offspring count $`t(S, S')/(s + a_S + r_0(S))`$. Summing over $`\ell \ge 0`$ by monotone convergence gives (2), bounded absolutely by $`W(J_{1/i})/(s - \rho)`$.

---

## 3. Rational Truncations $`K^-_N(s)`$ and $`K^+_N(s)`$, Nesting, and Endpoints (§3)

- **Finiteness, boundary diagonal, and M-matrix invertibility:** For $`N \ge 2`$, the domain $`\mathcal D_N`$ of states reachable from the 19 newborns within size $`\le N`$ is finite (clamped marks take only 19 values). Since copies have size $`\le n(S)`$, only intrinsic transitions from boundary states $`n(S) = N`$ exit to size $`N + 1`$. Retaining the full subtraction $`a_S + r_0(S)`$ on the diagonal (including on the boundary $`n(S) = N`$) kills particles upon exit to size $`N + 1`$. Since $`R_N`$ is Metzler and $`(A_N W)(S) \ge ((s - R)W)(S) \ge (s - \rho) W(S) > 0`$ on $`\mathcal D_N`$ for $`s > \rho`$, $`A_N = sI - R_N`$ is a nonsingular M-matrix with entrywise nonnegative inverse $`A_N^{-1} \ge 0`$ (rational for $`s \in \mathbb Q`$).
- **First-exit enclosure (4):** Partitioning walks in (2) by whether they exit $`\mathcal D_N`$ gives on $`\mathcal D_N`$ the exact identity
  $`\displaystyle U_s^{(j)} = A_N^{-1} B_j + A_N^{-1} \Big[S \mapsto \sum_{n(S')=N+1} t(S, S') U_s^{(j)}(S')\Big].`$
  Inserting $`0 \le U_s^{(j)}(S') \le W(S')/(s - \rho)`$ yields $`K^-_N(s) \le K(s) \le K^+_N(s)`$ entrywise, and monotonicity of the spectral radius of nonnegative matrices gives $`r(K^-_N(s)) \le r(K(s)) \le r(K^+_N(s))`$.
- **Nesting in $`N`$ and sharper test-column bound (line 55):**
  - $`K^-_N \le K^-_{N+1}`$ holds since $`\mathcal D_N \subset \mathcal D_{N+1}`$.
  - For $`K^+_{N+1} \le K^+_N`$, let $`v_N = A_N^{-1}(B_j + H_N)`$ on $`\mathcal D_N`$ and $`v_N = W/(s-\rho)`$ outside $`\mathcal D_N`$. First, $`A_N(W/(s-\rho) - v_N) \ge W - B_j \ge 0`$ on $`\mathcal D_N`$, so $`v_N \le W/(s-\rho)`$ on $`\mathcal D_N`$. On $`\mathcal D_{N+1}`$, $`(A_{N+1} v_N)(S) = B_j(S) = (A_{N+1} v_{N+1})(S)`$ for $`n(S) \le N`$; for $`n(S) = N + 1`$, using $`v_N(S') \le W(S')/(s-\rho)`$ on copies $`S' \in \mathcal D_N`$ returning into the interior gives
    $`\displaystyle (A_{N+1} v_N)(S) \ge \frac{(s - R)W(S)}{s - \rho} + H_{N+1}(S) \ge W(S) + H_{N+1}(S) \ge B_j(S) + H_{N+1}(S) = (A_{N+1} v_{N+1})(S).`$
    Multiplying by $`A_{N+1}^{-1} \ge 0`$ gives $`v_{N+1} \le v_N`$, hence $`K^+_{N+1}(s) \le K^+_N(s)`$.
  - For any test vector $`g \ge 0`$, the total direct reward $`\sum_j g_j B_j(S) \le \|g\|_\infty W(S)`$ implies $`U_s g(S') \le \|g\|_\infty W(S')/(s - \rho)`$, so $`K(s) g \le K^-_N(s) g + u_N \|g\|_\infty`$ (avoiding the factor-of-19 inflation of $`K^+_N g`$).
- **Endpoint $`N = 2`$:** $`\mathcal D_2`$ consists solely of the 19 newborns ($`A_2`$ diagonal with entries $`s + 1/i + 3/2`$), $`K^-_2(s)`$ has only column $`j = 2`$ nonzero (reducible, $`r(K^-_2(s)) = \frac{1/2}{s+2}`$), and $`H_2(J_{1/i}) = \frac{941}{520} \frac{W(J_{1/i})}{s - \rho}`$, giving valid enclosures even at $`N = 2`$.

---

## 4. Estimates (5)–(7) and Uniform Inverse Slope (§4)

- **Estimate (5):** Since every intrinsic jump adds $`+1`$ to $`n(S)`$, $`L^0(n^\delta W) = (n+1)^\delta L^0 W + ((n+1)^\delta - n^\delta) r_0 W`$ identically. For $`0 < \delta \le 1`$, $`(1 + 1/n)^\delta \le 1 + \delta/n`$, $`r_0 < n`$, $`n \ge 2`$, and $`(n+1)^\delta > 0`$ give $`L^0 V_\delta \le [q + \delta(1 + q/2)] V_\delta`$. Since copies have size $`\le n`$ and $`\sum_Y a_Y W(Y) = (a_S + d) W(S) - h(c_{\rm root})`$, $`R V_\delta \le [\rho + \delta(1 + q/2)] V_\delta - n^\delta h(c_{\rm root}) \le [\rho + \delta(1 + q/2)] V_\delta`$.
- **Estimate (6):** For $`\rho + \delta(1 + q/2) \le s`$, restricting $`(s - R) V_\delta \ge 0`$ to $`\mathcal D_N`$ gives $`A_N V_\delta|_{\mathcal D_N} \ge \sum_{n(S')=N+1} t(S, S') V_\delta(S') = (N+1)^\delta (s - \rho) H_N`$. Multiplying by $`A_N^{-1} \ge 0`$ at $`S_0 = J_{1/i}`$ ($`n(J_{1/i}) = 2`$) yields $`0 \le u_N(i) \le (2/(N+1))^\delta W(J_{1/i})/(s - \rho)`$ (equivalently via the absorbing exit account in line 69).
- **Constants and estimate (7):** With $`l = 113/200`$, $`u = 1133/2000`$, $`\delta = 1/24`$:
  - $`\rho + \delta(1 + q/2) = 51/100 + 59/1200 = 671/1200 < 678/1200 = l`$.
  - On $`s \in [l, u]`$, $`K_{ij}(s) \le W(J_{1/i})/(s - \rho) \le 1/[d(l - \rho)] = H = 4000/33`$, so $`r(K(s)) \le 19 H`$.
  - Along the pure root-leaf-gain path from $`J_{1/i}`$ through degrees $`k = 2, \dots, j`$ ($`2 \le j \le 20`$), each intermediate state has $`L = k - 1`$, $`0`$ core children, root gain rate $`k - 1`$, $`r_0 = k - 1/k = k - 1 + (k - 1)/k`$, death rate $`1/i \le 1/2`$, and terminal direct rate $`(j-1)/j`$ into label $`m(j) = j`$ (including both endpoints $`j = 2`$ with empty product and $`j = 20`$ at the clamp). Thus $`K_{ij}(s) \ge e > 0`$ uniformly on $`[l, u]`$ (exact value $`e \approx 3.294 \times 10^{-4}`$, attained at $`j = 20`$).
  - For the Perron vector $`g > 0`$ of $`K(s)`$ with $`\max g = 1`$, $`r(K(s)) g_i \ge e \implies g_{\min} \ge g_* := e/(19 H) \approx 1.430 \times 10^{-7}`$.
  - On $`g`$, $`((K^+_N - K^-_N) g)_i = u_N(i)\sum_j g_j \le 19 \max_i u_N(i) \le 19 H (2/(N+1))^{1/24} =: \varepsilon_N`$, giving $`K^+_N g \le (r(K) + \varepsilon_N/g_*) g`$ and $`K^-_N g \ge (r(K) - \varepsilon_N/g_*) g`$. Separate Collatz–Wielandt bounds hold for any nonnegative matrix (even reducible $`K^-_N`$, by iterating $`M^k g`$ and taking $`k`$-th roots, and trivially when $`r(K) - \varepsilon_N/g_* \le 0`$ since $`r(K^-_N) \ge 0`$), yielding (7): $`0 \le r(K^+_N(s)) - r(K^-_N(s)) \le (38 H / g_*)(2/(N+1))^{1/24}`$.
- **Uniform inverse slope:** At $`S_0 = J_{1/i}`$, $`r_0(J_{1/i}) = 2 - 1/2 = 3/2`$ and $`a_{J_{1/i}} = 1/i \le 1/2`$, so total exit rate is $`a_{J_{1/i}} + r_0(J_{1/i}) = 1/i + 3/2 \le 2`$ (equality at $`i = 2`$). For $`l \le s < s' \le u`$, the $`\ell = 0`$ direct-production term into label $`j = 2`$ (rate $`b = 1/2`$) gives
  $`\displaystyle (K_{i,2}(s) - K_{i,2}(s')) g_2 \ge g_* \int_0^\infty \tfrac{1}{2}(e^{-st} - e^{-s't}) e^{-2t}\,dt = \frac{g_*(s' - s)}{2(s+2)(s'+2)} \ge \frac{g_*(s' - s)}{2(u+2)^2}.`$
  Every other path contribution in (2) is nonnegative for $`s < s'`$, and $`g_i \le 1`$, so $`K(s) g \ge (r(K(s')) + \frac{g_*(s'-s)}{2(u+2)^2}) g`$, whence $`r(K(s)) - r(K(s')) \ge g_*(s' - s)/[2(u+2)^2]`$ by Collatz–Wielandt. Combined with (7), rational grid/bisection enclosures converge at an explicit power rate $`O((N+1)^{-1/24})`$.

---

## 5. Mathematical Form of the Determinant vs. Elementary / Named-Function Closed Forms

It is important to distinguish precisely what (1)–(2) gives from an elementary or named-function closed form:

- **What (1)–(2) is:** A **fixed-dimension ($`19 \times 19`$) Fredholm/secular determinant equation** whose 361 entries are **explicit, absolutely convergent combinatorial walk series** (2) of positive rational functions of $`s`$ over the countable state space of finite planted marked trees, depending only on the local friend-tree transition rates (no reference to $`\mu`$ or any unknown limiting distribution), accompanied by **rigorous monotone rational matrix enclosures** $`K^-_N(s) \le K(s) \le K^+_N(s)`$ with explicit convergence rate $`O((N+1)^{-1/24})`$ and explicit inverse slope.
- **What (1)–(2) is not:** It is **not** an algebraic equation, nor an elementary closed form (a finite expression in rational functions, radicals, exponentials, logarithms, and trigonometric functions), nor a formula in standard named special functions ($`\Gamma`$, $`\zeta`$, hypergeometric, Bessel, Painlevé, etc.). The finite $`19 \times 19`$ determinant encapsulates an infinite-dimensional tree-walk resolvent inside its entries, summing over arbitrarily large finite trees; while it provides a self-contained analytic characterization and an effective (though with exponent $`1/24`$, computationally coarse) rational approximation procedure, it does not express $`\mu`$ in elementary or named-function closed form.

---

## 6. Independent Rational Matrix Test and Planted Errors (`finite_kernel_neon.py`)

Replay command from the repository root (1 CPU, exact `fractions.Fraction` arithmetic, wall time 1.13 s, exit code 0):

```text
python3 round2/checks/finite_kernel_neon.py
```

Verified checks:
1. **Clamped Lyapunov inequality:** All 31,623 exact rational pairs ($`2 \le D \le 250`$, $`0 \le c \le D`$) with $`m(D) = \min(D, 20)`$ pass (0 failures; tightest pair $`(c, D) = (1, 3)`$, margin $`1/4350`$; tail factor $`1/4 < 63/250`$; $`M = 19`$ fails).
2. **Clamped signed-to-positive transfer identity:** Exhaustive brute-force labelled physical check at $`M = 3`$ (79 reachable marked trees for $`n = 3 \dots 7`$, exercising 18 strictly clamped $`D > M`$ births with $`a_v > p_u`$ at birth) and $`M = 20`$ (88 trees, $`n = 3 \dots 7`$): 0 failures.
3. **Exact $`19 \times 19`$ rational truncations & path sum (2):** Built exact $`A_N, K^-_N, K^+_N, u_N`$ for $`M = 20`$ at $`N = 2`$ ($`19`$ states), $`N = 3`$ ($`57`$ states), and $`N = 4`$ ($`152`$ states, cyclic) at both $`s = l = 113/200`$ and $`s = u = 1133/2000`$. Verified $`r_0(S) = n(S) - 1/D_{\rm root}`$, $`t(S, S) = 0`$, $`A_N W \ge (s - \rho) W > 0`$, $`A_N^{-1} \ge 0`$, exact nesting $`K^-_2 \le K^-_3 \le K^-_4`$ and $`K^+_4 \le K^+_3 \le K^+_2`$ on all 361 entries, exact 24th-power rational verification of (6) on all 19 rows, sharper test-column bound (line 55), exact match of path sum (2) to $`K^-_3(l)`$ at $`\ell \le 3`$ on acyclic $`\mathcal D_3`$, and cyclic convergence on $`\mathcal D_4`$ ($`\ell \le 10`$ residual $`9.574 \times 10^{-10}`$).
4. **Estimates (5), (7), constants, and inverse slope:** Verified exact $`L^0 V_1`$ and $`R V_1`$ bounds and exact rational inequality $`(1 + 1/(24n))^{24} \ge 1 + 1/n`$ on all 152 states of $`\mathcal D_4`$; verified $`H = 4000/33`$, $`e \approx 3.294 \times 10^{-4}`$ (attained at $`j = 20`$), $`g_* \approx 1.430 \times 10^{-7}`$, bound (7) at $`N = 4`$, exact exit rate $`1/i + 3/2 \le 2`$ for $`J_{1/i}`$, and inverse slope $`\ge 1.086 \times 10^{-8}(s' - s)`$.
5. **Small 2-label ($`M = 3`$) monotone enclosures ($`N = 2 \dots 5`$, up to 44 states):** Verified monotone tightening $`r(K^-_2) = 0.1923 < \dots < r(K^-_5) = 0.3566 \le r(K^+_5) = 19.116 < \dots < r(K^+_2) = 29.157`$.
6. **Planted errors (all 10 caught):**
   - Error 1: omitting exiting $`r^0`$ rate from boundary diagonal ($`K^-_{2, \text{bad}} > K^-_3`$),
   - Error 2: dropping copy terms in $`t(S, S')`$ ($`K^-_3`$ mismatch at column $`j = 3`$),
   - Error 3: using $`r_0 = n`$ instead of $`n - 1/D_{\rm root}`$ ($`K^-_2`$ mismatch),
   - Error 4a: replacing $`W/(s-\rho)`$ by $`W/s`$ (fails supersolution $`(s-R)(W/s) \ge W`$ at $`J_{1/20}`$, where $`RW = 135/104 > 0`$),
   - Error 4b: using $`H_N`$ denominator $`s + 5`$ (violates upper nesting $`K^+_3 \le K^+_2`$),
   - Error 5: zeroing one boundary entry of $`H_2`$ (violates $`K^+_3 \le K^+_2`$),
   - Error 6: $`\delta = 1/10`$ ($`\rho + \delta(1+q/2) = 157/250 > l = 113/200`$),
   - Error 7: clamping at $`M = 19`$ instead of $`20`$ (fails Lyapunov inequality at $`c = 4, D = 135`$: $`524/2079 > 63/250`$),
   - Error 8: dropping killing rate $`1/2`$ from $`J_{1/2}`$ survival in inverse slope (exceeds exact $`\ell = 0`$ difference),
   - Error 9: replacing $`\|g\|_\infty`$ by $`\min(g)`$ in the line 55 test-column bound (violated at $`N = 2`$ vs. $`N = 3`$ on $`g = \mathbf 1_{\{j=3\}}`$).
