CORRECT

# Cold Referee Review of `round2/finite/one-port-upper.md`

- **Target reviewed**: `round2/finite/one-port-upper.md` (`finite/one-port-upper.md`) (206 lines, Sections 1–9, eqs. (1)–(12); SHA-256 `c80dc0b3d6513e03e7a4f1162bac050cf64be92714a5acc12184725acd9f619a`).
- **Background reference checked**: `proofs/existence.md` (Lemmas 1.1–1.4, 2.1, 3.1–3.2, 4.1, 5.1–5.2, Corollary 5.3, Lemma 6.1, Proposition 6.2, Section 7).
- **Independent checker**: `round2/checks/one_port_neon.py` (`one_port_neon.py`) (stdlib-only, single CPU, exact small-tree killed-CTMC uniformization, exact rational jump-chain enumeration, Monte Carlo coupling identities, and planted-error falsification).

## Verdict Summary

All claims, intermediate lemmas, couplings, constants, parameter schedules, explicit finite-time error bounds (11)–(12), closed-form negative-binomial tail bound, and asymptotic convergence rates in Sections 1–9 hold without mathematical error or logical gap. Every single step was independently re-derived and stress-tested for edge cases ($`t = 0`$, non-integer $`t`$, $`\beta < \gamma`$ vs $`\beta > \gamma`$, $`\beta = 0`$, $`a = 0`$, boundary promotions at $`L_v = K`$, integer jump times of $`b_y`$, and the $`T_2 \to T_3`$ distinction now explicitly stated on line 3).

---

## Section-by-Section Independent Re-Derivation

### 1. Setup, Parameters, and Schedule (Section 1, eqs. (1)–(2))
- **Preconditions**: The interval starts immediately after an existing core $`p`$ promotes one of its leaves $`v`$, so $`v`$ has degree $`2`$ (edge $`pv`$ plus one leaf), $`p`$ has degree $`\ge 2`$, and both endpoints remain core forever. Thus $`a_x(u) = 1/D_v(x_{u-}) \le 1/2`$ and $`b_y(u) = 1/D_p(y_{u-}) \le 1/2`$ hold permanently. (The separate founding-edge transition $`T_2 \to T_3`$, where $`p`$ has degree $`1`$ and $`b_y = 1 > 1/2`$, is explicitly excluded from (1) and handled in Section 8.)

• **Parameter constraints**: For integer $`j \ge 4`$ with $`d = B = j`$, $`J = 2j`$, $`H = j^2`$, $`\kappa = 1 - 1/j`$, $`b = 2j`$, $`K = \max(1, \lceil j c_d \rceil)`$, $`M = K + 1 + \lceil j(z-1)K \rceil`$:
- $`d \ge 4 \ge 3`$, $`J \ge 8`$ even, $`B, H \ge 4 \ge 2`$, $`K \ge 1`$, $`M > K`$, $`s = 1 - 1/(2j) > 1 - 1/j = \kappa > 0`$.
- $`V = d\gamma + c_d/K + (z-1)K/M \le 1 + 1/j + 1/j = 1 + 2/j`$, so $`V/d \le 1/j + 2/j^2`$.

◦ Summing $`\Lambda = 1/J + (1-\kappa)/2 + V/d + 2/d + \gamma + \delta`$:

```math
\Lambda \le \frac{1}{2j} + \frac{1}{2j} + \left(\frac{1}{j} + \frac{2}{j^2}\right) + \frac{2}{j} + \frac{1}{j} + \left(\frac{1}{j} - \frac{1}{j^2}\right) = \frac{6}{j} + \frac{1}{j^2} < \frac{6}{j} + \frac{2}{j^2}.
```

- At fixed $`j`$, $`A(t) = O_j((2+t)^{J/2}(t+\gamma^{-1}))`$, so $`\log A(t) = O_j(\log(2+t))`$.

### 2. Joint Path Density Across the Cut (Section 2, eq. (3))
- Projecting onto the cut-labeled state space $`(X, Y)`$ (rooted core trees on each side with newly promoted cores labeled by promotion order separately on each side, and leaf counts $`L_w`$ per core) is Markovian with zero symmetry/collision factors: all transitions from any state produce distinct next states, and no jump alters both $`X`$ and $`Y`$.
- Total intrinsic exit rate on the $`v`$-side is $`\Lambda_{\mathrm{int}}(x) = |V(x)| - a_x`$, and on the $`p`$-side is $`\Lambda_{\mathrm{ext}}(y) = |V(y)| - b_y`$, so the coupled exit rate is $`\Lambda_{\mathrm{int}}(x) + b_y + \Lambda_{\mathrm{ext}}(y) + a_x = |V(x)| + |V(y)|`$, matching the exact friend-tree generator.

• Merging $`m`$ interior jumps and $`n`$ exterior jumps over all $`\binom{m+n}{m}`$ interleavings factors the joint path measure on the product simplex into

```math
I(x)\,dx \; Q^{a_x}(dy) \; e^{-\int_0^t b_y(u)\,du} \prod_{u \in G(x)} [r_x(u) + b_y(u)],
```

where $`Q^{a_x}(dy)`$ is the normalized non-explosive prescribed-input exterior law (total birth rate $`\le |V(y)| + 1/2 < |V(y)| + 1`$). Tonelli over $`(m, n)`$ gives (3) on the full path space with no conditioning of the true exterior on $`x`$.
- **Exact small-tree check**: On the 5-vertex cut seed ($`x_0`$: $`v`$ + 1 leaf; $`y_0`$: $`p`$ + core neighbor $`u`$ + 1 leaf), integrating eq. (3) via Gauss-Legendre quadrature and exact piecewise uniformization of killed $`Q^{a_x}`$ agrees with the exact killed full coupled CTMC across all 16 interior states to $`L_\infty`$ error `3.181e-15`, and both agree with independent vertex-ringing Monte Carlo. Planting errors (dropping cut slot, omitting $`e^{-\int b_y}`$, omitting core children from $`r_x`$) breaks the identity by `6.0e-3` to `1.7e-1`.

### 3. Exterior Estimate for One Late Arrow (Section 3, eq. (4))
- Under $`Q^a`$ ($`0 \le a \le 1/2`$), after the first $`b`$ prescribed immigrants, $`D_p \ge b`$ permanently. Coupling each subsequent tracked leaf with an independent rate-$`1/b`$ death proposal (always removed from the tracked family, thinned with probability $`b/D_p \in (0, 1]`$ for true promotion) preserves exact marginal leaf-ringing ($`L_p`$) and promotion ($`L_p/D_p`$) rates while ensuring every tracked leaf is an actual leaf of $`p`$, and making the tracked families i.i.d. linear birth-death processes (birth $`1`$, death $`1/b`$, independent of arrival times).
- Survival probability is $`s = 1 - 1/b`$. In operational time $`A(u) = \int_0^u a(w)\,dw`$ (extending $`a`$ by $`1/2`$ past $`u`$), the operational arrival time $`W`$ of the first surviving family after the first $`b`$ immigrants has law $`\mathrm{Gamma}(b, 1) + \mathrm{Exp}(s)`$ independently, giving $`\mathbb E e^{\kappa W} = (1-\kappa)^{-b} s/(s-\kappa)`$ for $`0 < \kappa < s`$.
- Solving the backward PGF ODE $`\partial_h f = (f-1)(f-(1-s))`$, $`f(0, z) = z`$ gives $`\Pr(Z_h = 0) = \frac{(1-s)(1-r)}{1-(1-s)r}`$ and $`\mathbb E[z^{Z_h} \mid Z_h > 0] = \frac{qz}{1-(1-q)z}`$ with $`r = e^{-sh}`$ and $`q = \frac{sr}{1-(1-s)r} \le e^{-sh}`$.
- For $`G \sim \mathrm{Geom}(q)`$ on $`\mathbb Z_{\ge 1}`$, $`\mathbb E[1/G] = \frac{-q\log q}{1-q} \le \sqrt{q}`$ because $`g(x) = x - 1/x - 2\log x`$ has $`g(1) = 0`$ and $`g'(x) = (1 - 1/x)^2 \ge 0`$ for $`x = q^{-1/2} \ge 1`$.
- Using $`\mathbf 1_{\text{surv}} \le \mathbf 1_{\{Z_h > 0\}}`$ gives $`\mathbb E[1/Z_h \mid \text{surv}] \le s^{-1} e^{-sh/2}`$. Because $`a \le 1/2`$, physical age on $`\{W \le A(u)\}`$ satisfies $`h \ge 2(A(u) - W)`$, canceling the $`1/2`$ factor in the exponent: $`e^{-sh/2} \le e^{-s(A(u)-W)} \le e^{-\kappa(A(u)-W)}`$. Combining with $`1/D_p \le 1 < s^{-1} e^{\kappa(W - A(u))}`$ on $`\{W > A(u)\}`$ yields $`s^{-1} e^{-\kappa A(u)} \mathbb E e^{\kappa W} = D_0 e^{-\kappa A(u)}`$, matching (4) with $`D_0 = (1-\kappa)^{-b}/(s-\kappa)`$.
- **Planted-error check**: Choosing $`\kappa = 1.3 > 1`$ fails (4) asymptotically against the zero-immigrant lower bound $`\frac{1}{2} e^{-A(u)}`$; dropping $`D_0`$ (setting $`D_0 = 1`$) fails (4) at $`a = 1/2, u = 10, \kappa = 0.99`$ (`8.821e-3 > 7.083e-3`).

### 4. Last Exterior Arrival and Finitely Many Profiles (Section 4, eq. (5))
- Telescoping expansion $`\prod_{j=1}^k (r_j + b_j) = \prod_{j=1}^k r_j + \sum_{m=1}^k \big(\prod_{j<m}(r_j + b_j)\big) b_m \big(\prod_{j>m} r_j\big)`$ is exact; re-indexing the last exterior gain time as $`\tau \in (0, t)`$ partitions the jump simplex bijectively into prefix $`(0, \tau)`$, deterministic forced root gain at $`\tau`$, and suffix $`(\tau, t)`$.
- Rounding $`b_y(k) = 1/D_p(y_{k-})`$ up to multiples of $`1/J`$ on each unit interval $`[k, k+1)`$ ($`k = 0, \dots, m-1`$, $`m = \lceil t \rceil`$) gives a nonincreasing profile $`\beta \ge b_y`$ with:
  - $`|\mathcal P_t| = \binom{m + J/2}{J/2} \le (m + 1)^{J/2} \le (2 + t)^{J/2}`$ (since $`\lceil t \rceil < t + 1`$ for all $`t > 0`$, and $`2^{J/2}`$ covers $`t \in (0, 1]`$).
  - Integral error $`\int_0^\tau (\beta - b_y)\,du \le \frac{\tau}{J} + \sum_{k=0}^{m_\tau - 2}[b_y(k) - b_y(k+1)] + [b_y(m_\tau - 1) - b_y(\tau)] = \frac{\tau}{J} + b_y(0) - b_y(\tau) \le \frac{t}{J} + \frac{1}{2}`$.
- Summing over $`\beta \in \mathcal P_t`$ decouples the interior factor from $`y`$, leaving $`\int b_y(\tau)\,Q^{a_x}(dy) \le D_0 e^{-\kappa A_x(\tau)}`$ by (4), yielding (5).
- **Exact & planted-error checks**: Exact killed-CTMC last-arrival decomposition matches full CTMC to `6.661e-16`. Dropping $`+1/2`$ total variation fails on an early-drop step function (`0.4995 > 0.25`); replacing $`(2+t)^{J/2}`$ by $`(1+t)^{J/2}`$ fails at $`t = 0.5, J = 2`$ (`2 > 1.5`).

### 5. Rényi Likelihood Domination (Section 5, eq. (6))
- Since $`\gamma = 1/B > 0`$, $`R_{\beta, \tau} \ll P = R_{\gamma, \tau}`$ even when $`\beta = 0`$.
- On the predictable set $`E_{\ge} = \{L_v(x_{w-}) \ge K\}`$, $`r_x \ge K \ge 1`$, $`|\Delta| = |\beta - \gamma| \le 1/2`$, and $`u = 1 + \Delta/(r + \gamma) \in (0, 3/2]`$ for BOTH $`\beta \ge \gamma`$ and $`\beta < \gamma`$. Lagrange Taylor remainder for $`g(\xi) = (1+\xi)^d - 1 - d\xi`$ ($`d \ge 3`$) gives $`(r+\gamma)(u^d - 1) - d(\beta - \gamma) \le \frac{1}{2} d(d-1)(3/2)^{d-2} \frac{\Delta^2}{r+\gamma} \le \frac{d(d-1)(3/2)^{d-2}}{8K} = \frac{c_d}{K}`$.
- The auxiliary law $`\widetilde P`$ changes the root gain rate to $`(r+\gamma)u^d`$ only on $`E_{\ge}`$ (total rate bounded by $`C_d(|V| + 1)`$, non-explosive, equivalent to $`P`$, preserving all promotion intensities).
- On $`E_{<} = \{L_v < K\}`$: compensator is $`\le d\gamma t`$; per-jump factor $`u^d \le (1 + 1/(2\gamma))^d = z`$. Tracking $`\Phi(w) = \min(L_v(x_w), K)`$ after $`\sigma_M = \inf\{w: D_v \ge M\} \wedge \tau`$ gives at most $`M-2 \le M`$ gains before $`\sigma_M`$, and after $`\sigma_M`$ each low-leaf gain adds $`+1`$ to $`\Phi`$ while decreases occur ONLY at root promotions with $`1 \le L_v \le K`$ (count $`N`$), giving $`G_{\mathrm{low}} \le M + K + N`$ pathwise. Under $`\widetilde P`$, the predictable intensity of $`N`$ is $`\le K/M`$, so $`M_w = z^{N_w} \exp(-(z-1)\int_0^w \lambda_N)`$ is a nonnegative local supermartingale giving $`\mathbb E_{\widetilde P}[z^N] \le e^{(z-1)Kt/M}`$, proving (6).
- **Exact & planted-error checks**: Exact tilted-generator uniformization confirms $`\mathbb E_P[Z] \approx 1`$, $`\mathbb E_P[L] \approx 1`$, and (6) across $`\beta > \gamma`$ and $`\beta = 0 < \gamma`$; 5000 $`\widetilde P`$-trajectories gave `0` pathwise violations of $`G_{\mathrm{low}} \le M + K + N`$. Dropping $`(3/2)^{d-2}`$ from $`c_d`$ fails Taylor bound (`1.9436 > 1.5`); omitting $`L_v = K`$ promotions from $`N`$ fails pathwise on oscillating trajectories (`12 > 6`); dropping $`d\gamma`$ from $`V`$ fails at $`\beta = 0, r = 0`$ for large $`t`$.

### 6. Hölder Exponents and Second Moment (Section 6, eq. (7))
- With $`q = d/(d-1) \in (1, 2)`$, $`\theta = (2-q)/q = 1 - 2/d`$, $`(q-1)/q = 1/d`$, and $`A_x(\tau) \le \tau/2 \le t/2`$:
  - Under $`P = R_{\gamma, \tau}`$, $`C_B(0) = 1`$, forced leaf gain at $`\tau`$ does not change $`C_B`$, and total promotion rate is $`\sum L_w/D_w \le C_B`$, so $`\mathbb E_P[f^2] \le \mathbb E_P[C_B(t)^2] \le 2 e^{2t} - e^t \le 2 e^{2t}`$.
  - Hölder with $`(d, q)`$ then interpolation $`f^q = f^{2-q}(f^2)^{q-1}`$ with $`(1/(2-q), 1/(q-1))`$ yields (7) cleanly.

### 7. Physical Finite Seed $`T_*`$, Linear Bound, and $`\pi_*`$ (Section 7, eqs. (8)–(9))
- $`T_*`$ has $`n_* = 3 + (B-1)H`$ vertices, $`D_p(0) = B = \gamma^{-1}`$, $`L_p(0) = 0`$, and $`B-1`$ exterior core neighbors $`u_i`$ of $`p`$ with permanent degree $`D_{u_i} \ge H`$.
- Both inequalities in (8) were independently re-derived in TWO ways: (i) via the author's independent-channel Poisson thinning coupling (where excluding $`u_i \to p`$ proposals on $`[0, t]`$ (prob $`e^{-\delta t}`$), conditioning on the last $`p \to v`$ proposal $`\tau`$ with density $`\gamma e^{-\gamma(t-\tau)}`$ and unconditioned rate-$`\gamma`$ Poisson prefix, and requiring no $`v \to p`$ arrow on $`[0, \tau)`$ (cond. prob $`e^{-A_x(\tau)}`$) keeps $`D_p = B`$ identically through $`\tau`$ so every prefix proposal and the last one are accepted), and (ii) directly from the exact path density (3) of Section 2.
- Realization mass $`\pi_* = \Pr(N(1) = n_*) \prod_{j=2}^{n_*-1} \frac{1}{j(j-1)}`$ matches Lemma 1.3–1.4 of `proofs/existence.md` and yields $`S(t) \le \pi_*^{-1} F(t+1)`$ by the Markov property.
- Since $`h^\theta \le 1 + h`$ and $`S(t) \ge 1`$, (9) follows, and substituting (6)–(9) into (5) reproduces (1) and (2) with every single factor exact.
- **Exact & planted-error checks**: Exact jump-chain enumeration confirms $`\Pr(T_{n_*} \cong T_*) \ge \prod_{j=2}^{n_*-1} \frac{1}{j(j-1)}`$ for $`(B, H) \in \{(2,2),(2,3),(3,2)\}`$. Independent channel-thinned Monte Carlo verifies both coupling identities directly (`0/2697` and `0/1516` violations of $`D_p = B`$ before $`\tau`$ on $`E_2`$). Planting a leaf at $`p`$ ($`L_p(0)=1, D_p(0)=B`$) immediately breaks $`D_p = B`$ on `243/1614` (`15.1%`) of $`E_2`$ trajectories.

### 8. Uniform Upper Bound and Explicit Block-Time Errors (Sections 8–9, eqs. (10)–(12))
- For any finite tree $`T`$ with $`C(T) \ge 1`$, core ancestry partitions all cores born in $`(0, t]`$ into disjoint branches rooted at direct children of initial cores. Strong Markov at each direct promotion time, (1), and compensation with direct promotion intensity $`L_u/D_u \le 1`$ per initial core yield (10).
- Deterministic iteration from $`T_3`$ plus $`F(u) \le \mathbb E_{T_3} C(u)`$ gives $`\mu \le t^{-1}\log U(t)`$.
- From $`F'(u) \le F(u) + 2 e^{-2u}`$, multiplying by $`e^{-u}`$ and integrating over $`[t, t+1]`$ gives $`F(t+1) \le e(F(t) + \frac{2}{3} e^{-2t}) \le e(F(t) + 1)`$ (valid for all $`t \ge 0`$). Also $`F(t) \ge 1 - e^{-2t}`$ (and in fact $`F(t) \ge t`$ since instantaneous promotion rate $`\ge 1`$ on every tree, including $`T_2`$ where it equals $`2`$).

• **Section 9 checks**:
- Eq. (11) follows directly from (2), (10), $`1 \le c^{-1} F(t)`$, and $`F(t+1) \le e(1 + c^{-1}) F(t)`$ for $`t \ge 1`$.
- Eq. (12) matches Proposition 6.2 and Lemma 5.1 of `proofs/existence.md` with $`M_0 = j, \epsilon_0 = e^{-j}, k = \lceil 6je^j \rceil, h = k+1`$: for $`x = t/(jk) \ge 2`$, $`1 + \log(6x) \le 2x`$ holds strictly (at $`x=2`$: $`3.4849 \le 4`$, derivative $`2 - 1/x \ge 3/2 > 0`$); $`\theta_k(t) \le 1`$; $`F(t) \ge t \ge 2(k+1)`$ makes the Proposition 6.2 bracket $`\ge F(t)/2`$; and $`f(t) \le 1 + \log 2 < 2`$ for $`t \ge 1`$ gives $`f(t) - \mu \le \eta_j + \frac{2h + \log(8k/p_j)}{t}`$ (and $`\le 2`$ trivially).

◦ Exact algebraic derivation confirms the closed-form negative-binomial first-moment tail (line 197):

```math
\sum_{n=N+1}^\infty n(n-1) q^2 (1-q)^{n-2} = (1-q)^{N-1}\! \left[N(N+1)q + 2(N+1)(1-q) + 2(1-q)^2/q\right],
```

verified with exact `Fraction` equality across $`(q, N)`$.
- Asymptotic rate $`|t^{-1}\log F(t) - \mu| = O(\log\log t / \log t)`$: for $`j = \lfloor \log t / (4\log\log t) \rfloor`$, $`\log A_j(t) \le \exp(j\log j + O(j)) + O(j\log t) \le t^{1/4 + o(1)}`$, while $`2h + \log(8k/p_j) = \exp(2j + O(\log j)) = t^{o(1)}`$ and $`2jk = t^{o(1)} \le t`$.

---

## Independent Verification Run Output

Command executed from repository root (single CPU):
```bash
python3 round2/checks/one_port_neon.py
```

Observed output (exit code `0`):
```text
=== TEST 1: Section 2 Joint Path Density (eq. 3) ===
  Exact killed mass: full=0.941893863167, eq(3)=0.941893863167
  Exact E[C_B(t) 1_trunc]: full=1.102444709671, eq(3)=1.102444709671
  Max marginal L_inf error across 16 interior states: 3.181e-15
  Independent vertex-ringing MC: mass=0.94068 (exact 0.94189), E[C_B 1_trunc]=1.10157 (exact 1.10244)
  Planted errors caught: 1a(no cut slot) diff=1.1170e-01, 1b(no exp(-int b_y)) diff=1.7028e-01, 1c(bad r_x) diff=6.0341e-03
  [PASS] Test 1

=== TEST 2: Section 3 Birth-Death Law & Exterior Bound (eq. 4) ===
  Birth-death PGF, geometric conditional law, and 1/Z_h bounds verified to <1e-10.
  Q^a (const_0.5, u=4.0): E[1/D_p]=0.14078+/-0.00251 <= 78.43599 (D_0=260.42)
  Q^a (step_0.5_to_0.1, u=4.0): E[1/D_p]=0.18893+/-0.00305 <= 126.75840 (D_0=260.42)
  Q^a (zero_gap, u=4.0): E[1/D_p]=0.22556+/-0.00270 <= 123.01212 (D_0=260.42)
  Planted errors caught: 2a(kappa=1.3>1) lb=2.270e-05 > rhs=2.260e-05; 2b(drop D_0, kappa=0.99, u=10) E[1/D_p]>=8.82117e-03 > exp(-kA)=7.08341e-03
  [PASS] Test 2

=== TEST 3: Section 4 Last-Arrival Expansion & Profile Rounding (eq. 5) ===
  Profile count |(P_t)| <= (2+t)^{J/2}, beta >= b_y, and int(beta-b_y) <= t/J+1/2 verified.
  Exact decomposition: full=1.016154063631, intr+last_ext=1.016154063631, diff=6.661e-16
  Eq (5) check (t=0.4): E[C_B(t)]=1.2277 <= RHS(5)=15.1608
  Planted errors caught: 3a(drop +1/2 TV) integ=0.4995 > t/J=0.25; 3b((1+t)^(J/2) too small) exact=2 > 1.5
  [PASS] Test 3

=== TEST 4: Section 5 Renyi Likelihood Domination (eq. 6) ===
  Taylor bound verified on dense grid (both beta>gamma and beta<gamma); max lhs/rhs = 0.69684 <= 1.
  Normalization check (max_j=9, tau=0.25): E_P[Z]=0.99999710, E_P[L_aux]=0.99999707
  Eq (6) (d=3, B=2, K=1, M=2, tau=0.25): E_P[Z^d]=0.999997 <= 2.3675e+03
  Eq (6) (d=3, B=3, K=1, M=2, tau=0.30): E_P[Z^d]=1.042091 <= 6.4723e+04
  Eq (6) (d=3, B=4, K=2, M=3, tau=0.20): E_P[Z^d]=1.030532 <= 5.9755e+08
  Eq (6) (d=4, B=3, K=2, M=3, tau=0.20): E_P[Z^d]=1.025861 <= 2.6621e+10
  Eq (6) (d=3, B=3, K=1, M=2, tau=0.30): E_P[Z^d]=1.071582 <= 6.4723e+04
  Pathwise P_tilde check (5000 paths): G_low <= M+K+N violations=0, E[z^N]=8.1526 <= 3.3677e+09
  Planted errors caught: 4a(drop (3/2)^(d-2)) 1.9436>1.5000; 4b(omit L_v=K from N) oscillating path 12 > M+K=6 (MC viols=0); 4c(drop d*gamma) fails for beta=0, r=0 at large t.
  [PASS] Test 4

=== TEST 5: Section 6 Holder Exponents & Second Moment (eq. 7) ===
  E_P[C_B^2]=1.91626 <= Yule=3.35089 <= 2e^(2t)=4.91921
  Eq (7): LHS=1.18432 <= RHS=1.94769 (h(tau)=1.12931, E_P[f^2]=1.51337)
  [PASS] Test 5

=== TEST 6: Section 7 Physical Seed T_*, Eq. (8), & pi_* ===
  T_*(B=2,H=2,n*=5): P(T_n*=T*)=8.333333e-02 >= 6.944444e-03; P(T(1)=T*)=1.139435e-02 >= pi_*=9.495291e-04
  T_*(B=2,H=3,n*=6): P(T_n*=T*)=9.166667e-02 >= 3.472222e-04; P(T(1)=T*)=9.903578e-03 >= pi_*=3.751355e-05
  T_*(B=3,H=2,n*=7): P(T_n*=T*)=4.012346e-03 >= 1.157407e-05; P(T(1)=T*)=3.288215e-04 >= pi_*=9.485236e-07
  T_*(B=2,H=2,t=0.4): S(t)=3.4751, D_p>B on E2 = 0/2697 (must be 0)
    Coupling id 1: MC E[C_B 1_E0]=0.99575 vs exact e^(-gamma t)E^R0=0.99536
    Coupling id 2: MC E[C_B 1_E2]=0.16895+/-0.00332 vs exact=0.16925
    Eq (8) bounds: E^R0=1.2157 <= 4.2445; int h=0.4580 <= 10.3685
  Planted error 6a (leaf at p): D_p>B before tau on E2 = 243/1614 (15.1% violations of D_p=B!)
  T_*(B=3,H=2,t=0.4): S(t)=4.6990, D_p>B on E2 = 0/1516 (must be 0)
    Coupling id 1: MC E[C_B 1_E0]=1.07010 vs exact e^(-gamma t)E^R0=1.06398
    Coupling id 2: MC E[C_B 1_E2]=0.09730+/-0.00261 vs exact=0.09541
    Eq (8) bounds: E^R0=1.2157 <= 5.3693; int h=0.4573 <= 24.0301
  S(0.3) <= 5e^0.3=6.749 <= pi_*^{-1} F_lo(1.3)=1190.4
  [PASS] Test 6

=== TEST 7: Section 1 & Section 8 Constants, Schedule, Eq. (1), (10) ===
  Schedule j=4..50 verified: Lambda <= 6/j + 1/j^2 < 6/j + 2/j^2 everywhere.
  At j=4: z=81.0, c_d=3.375, K=14, M=4495, n*=51, Lambda=1.560059 (<= 1.625000), log(1/pi*)=313.6, log A(1)=5287.1
  Exact E[C_n] for n=2..10: ['2:0.00000', '3:1.00000', '4:1.33333', '5:1.58333', '6:1.79722', '7:1.98947', '8:2.16666', '9:2.33247', '10:2.48922']
  MC check at t=1: F(1)=1.4585 (>=1), F(2)=2.9182 <= e(F(1)+2/3 e^-2)=4.2100
  Minimal params at t=0.5: F(t) in [0.81536, 0.83251] <= E_T3[C(t)]=1.52740 <= U(t)=4.427e+06
  Planted error 7 (leaf parent p, D_p=1): b_y=1 > 1/2 gives u=1.60 > 3/2, breaking Sec 4 & 5.
  [PASS] Test 7

=== TEST 8: Section 9 Explicit Errors (eqs. 11-12), Tail Formula & Asymptotics ===
  Closed-form NegBin first-moment tail formula (line 197) verified with EXACT Fraction equality.
  Section 9 threshold 1+log(6x)<=2x, theta_k<=1, and line 206 inequalities verified for j=4..14.
  Minimum instantaneous promotion rate across all trees n=2..10: 1.000000 >= 1.0
  [PASS] Test 8

ALL 8 VERIFICATION & PLANTED-ERROR SUITES PASSED.
```

## Minor Non-Blocking Notes
- In Section 1 (line 36), $`\Lambda \le 6/j + 2/j^2`$ can even be tightened to $`\Lambda \le 6/j + 1/j^2`$ since $`\delta = (j-1)/j^2 = 1/j - 1/j^2`$.
- In Section 8 (line 155), $`F(t+1) \le e(F(t)+1)`$ holds for all $`t \ge 0`$ (indeed $`F(t+1) \le e(F(t) + \frac{2}{3} e^{-2t})`$); $`t \ge 1`$ is only needed in Section 9 to bound $`1 \le c^{-1} F(t)`$.
- No changes were made to the author's draft.
