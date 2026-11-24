VERDICT: HOLDS

# Adversarial Audit of Astra Claim 3 — Sharp Mean and Pointwise Degree Law

- **Target:** Astra Claim 3 (`ft2-astra@8e4db79`), covering `round2/fringe/asymptotics.md` and `round2/degree-law.md` (with dependencies `round2/fringe/renewal.md` and `round2/operator.md`).
- **Auditor:** Gemini 4 Neon (independent model family from GPT-6.2 Astra).
- **Independent test script:** `round2/outputs/attack_astra_34/audit_claim3.py` (written from scratch; `[local directory]` left strictly read-only).

---

## 1. Summary Verdict

Every mathematical step in `asymptotics.md` and `degree-law.md` holds under line-by-line adversarial inspection, and all identities and constants pass independent exact-rational and numerical recomputation (with 7/7 planted-error controls caught). No mathematical bug or logical gap was found; only two minor expositional compressions are noted below, both easily filled with zero change to the proof.

---

## 2. Line-by-Line Proof Audit

### 2.1 `round2/fringe/asymptotics.md`

• **§1 (Discounted occupation, lines 17–34):**
- From `renewal.md` (5)–(6), $`R_t W \le e^{\rho t} W`$ ($`\rho = 51/100`$) and $`I_t W \le 1/d + 2e^{qt}`$ ($`q = 9/25`$). Termwise integration of the scalar Taylor series of $`e^{-zt}`$ against the uniform exponential majorant establishes holomorphy of $`G_z \in \mathcal{L}(\mathcal{B}_W)`$ and $`U_z \in \mathcal{L}(C(\mathcal{A}), \mathcal{B}_W)`$ on $`\Re z > 51/100`$, and of $`\widehat{I}_z \in \mathcal{B}_W^*`$ on $`\Re z > 9/25`$.

◦ **Continuity of $`j G_z f_k`$ on $`\mathcal{A} = \{1/2, 1/3, \dots\} \cup \{0\}`$ (including at $`a = 0`$):** Conditioning on the original particle's intrinsic path $`S(t)`$ (independent of root mark $`a`$) and death time $`\mathrm{Exp}(a)`$ gives

```math
(j G_z f)(a) = \int_0^\infty e^{-(z+a)t}\,\mathbb{E}_J^0\! \left[f(S(t)) + \sum_{Y\text{ at }S(t)}\big((a_Y - p_{S(t)})G_z f(Y) + p_{S(t)} G_z f(Y^+)\big)\right]dt.
```

Since $`W(Y^+) = W(Y)`$, $`\sum_Y a_Y W(Y) \le (a_S + d)W(S)`$, and $`\mathbb{E}_J^0[(a_S + d)W(S(t))] \le e^{qt}`$ with $`q = 9/25 < \rho < \Re z`$, the bracketed intrinsic expectation is $`O_{\Re z}(e^{qt})`$. Thus $`a \mapsto j G_z f(a)`$ is Lipschitz on $`\mathcal{A}`$ (including at $`a=0`$). The same bound gives equicontinuity of $`K_z`$ on unit balls (compactness by Arzelà–Ascoli) and holomorphy of $`z \mapsto j G_z f \in C(\mathcal{A})`$ and $`z \mapsto K_z \in \mathcal{L}(C(\mathcal{A}))`$ (since $`C(\mathcal{A})`$ is a closed subspace of $`\ell^\infty(\mathcal{A})`$ and $`\|j F\|_\infty \le d^{-1}\|F\|_W`$).
- **Resolvent formula (3):** For $`\Re z > \mu`$, $`r(K_{\Re z}) < 1`$ (`renewal.md` §4), so generation-by-generation absolute bounds at $`\Re z`$ justify norm convergence of $`\sum_{n \ge 0} K_z^n = (I - K_z)^{-1}`$ and Fubini on $`M_t f = \int_0^t I_u P_{t-u} f\,du`$. Confirmed to exact order $`t^4/4!`$ ($`O(z^{-5})`$) below.

- **§2 (Pole on $`\Re z = \mu`$ and aperiodicity, lines 35–51):**
  - Simplicity of eigenvalue $`1`$ of $`K_\mu`$ follows from strong positivity + compactness (if $`(I - K_\mu)u = cg`$, applying $`\ell`$ forces $`c=0`$). Local rank-1 Riesz projection $`P_z`$, normalized eigenvector $`\ell(g(z)) \equiv 1`$, and differentiation of $`K_z g(z) = \lambda(z) g(z)`$ give $`\lambda'(\mu) = \ell(K'_\mu g) < 0`$ strictly, yielding the residue $`g \otimes \ell / [-\ell(K'_\mu g)]`$ with positive denominator and formula (4).
  - **Full support of $`\ell`$ and aperiodicity (line 43):** Any nonzero continuous $`w \ge 0`$ on $`\mathcal{A}`$ must be positive at some $`1/m`$ ($`m \ge 2`$), because $`1/m \to 0`$ forces $`w(0) = 0`$ otherwise; strong positivity then yields $`\min_{\mathcal{A}} K_\mu w > 0`$, so $`\ell(w) = \ell(K_\mu w) > 0`$. Hence $`K_{\mu+i\tau} v = v \ne 0`$ implies $`|v| \le K_\mu |v|`$ and $`\ell(K_\mu |v| - |v|) = 0`$, forcing $`K_\mu |v| = |v|`$ everywhere on $`\mathcal{A}`$, so $`|v| = c g > 0`$ (in particular $`v(1/2) \ne 0`$). For every row $`a`$, keeping the original $`J_a`$ alive with no leaf or root event on $`[0, t]`$ (probability $`e^{-(2+a)t} > 0`$) gives the $`1/2`$-child component an absolutely continuous density $`\ge \frac{1}{2} e^{-(\mu+2+a)t}`$ bounded below on every compact interval in $`(0, \infty)`$. Equality in the complex triangle inequality forces $`e^{-i\tau t} v(1/2)`$ to have constant phase for a.e. $`t > 0`$, impossible for $`\tau \ne 0`$. Compactness (Fredholm alternative) gives invertibility of $`I - K_{\mu+i\tau}`$ and local analytic extension.
  - Strict positivity of $`\widehat{I}_\mu(U_\mu g)`$ (initial centre keeps $`b_o=1`$ until first promotion) and of $`\ell(j G_\mu f_k)`$ and $`\ell(j G_\mu(f_k - f_{k+1}))`$ (original $`J_a`$ reaches degree $`k`$ before any promotion/death and holds it on an interval) holds for all $`k \ge 2`$.

• **§3 (Wiener–Ikehara and fixed-size Yule transfer, lines 52–77):**
- $`M_t f_k`$ is nonnegative and nondecreasing in $`t`$ ($`M_0 f_k = 0`$). Rescaling $`S(x) = M_{x/\mu} f_k`$ gives $`\int_0^\infty e^{-wx} S(x)\,dx = \mu F(\mu w)`$, whose residue at $`w = 1`$ is $`\mu \cdot d_{f_k}/\mu = d_{f_k}`$, with analytic extension of $`\mu F(\mu w) - d_{f_k}/(w-1)`$ to a neighborhood of every point on $`\Re w = 1`$. Every hypothesis of Wiener–Ikehara (e.g. Debruyne–Vindas arXiv:1611.09765 Thm 1.1) holds, yielding $`M_t f_k \sim d_{f_k} e^{\mu t}`$.
- Initial centre adds $`\le 1`$; $`\mathrm{Exp}(2)`$ wait from $`T_2`$ to $`T_3`$ and dominated convergence on $`e^{-\mu t} \mathbb{E}_{T_2} N_{\ge k}(t) = \int_0^t 2 e^{-(2+\mu)u} [e^{-\mu(t-u)}\mathbb{E}_{T_3} N_{\ge k}(t-u)]\,du`$ yield $`d_k^* = \frac{2}{2+\mu} d_{f_k}`$ (eq. (5)).

◦ **Fixed-size transfer (lines 62–75):** Setting $`1 - q = e^{-s/x}`$ in $`\nu_x = x^{-(\mu+2)}\sum_{n \ge 2}(n-1)a_{n, k}\delta_{n/x}`$ gives the exact identity

```math
\int_0^\infty e^{-sy}\,d\nu_x(y) = x^{-(\mu+2)} \frac{(1-q)^2}{q^2}\,\mathbb{E}_{T_2} N_{\ge k}(-\log q) \sim d_k^* (1-q)^2 (xq)^{-(\mu+2)} \longrightarrow d_k^* s^{-(\mu+2)}.
```

*(Minor expositional note: the factor $`(1-q)^2 = e^{-2s/x} \to 1`$ is omitted in prose on line 70, but its limit is $`1`$.)* Tilting by $`e^{-y}`$ gives probability transforms converging to $`(1+s)^{-(\mu+2)}`$, the $`\mathrm{Gamma}(\mu+2, 1)`$ transform; weak convergence and continuity of the limit density yield $`\nu_x([a, b]) \to \frac{d_k^*}{\Gamma(\mu+2)} \int_a^b y^{\mu+1}\,dy`$. The monotonicity sandwich on $`[1, 1+h]`$ and $`[1-h, 1]`$ gives $`c_k = d_k^* / \Gamma(\mu+2)`$ (confirmed $`\Gamma(\mu+2)`$, not $`\Gamma(\mu+1)`$).
- **Tightness of $`p_k`$ (line 76):** Nonexplosion of the residual family and $`R_t W(J_a) \le e^{\rho t} W(J_a)`$ justify DCT for $`j G_\mu f_k(a) \downarrow 0`$ with uniform bound $`\frac{1}{d(\mu-\rho)}`$; integrating against the finite positive measure $`\ell`$ gives $`c_k \downarrow 0`$ and $`\sum_{k \ge 2} p_k = 1`$.

- **§4 (Link to $`\pi`$, lines 78–80):**
  - Along the tilt $`w(t) = e^{-\mu t - \delta|t-T|}`$ of `operator.md` §2, $`w(t)/w(0) = e^{-\mu t - \delta(|t-T|-T)} \le e^{-(\mu-\delta)t} \le 1`$, so the normalized mass on any bounded $`[0, R]`$ is at most $`(w(0)/Z) \int_0^R e^t\,dt \to 0`$. Since $`m_t(f_k)/m_t(1) \to d_{f_k}/d_1 = c_k/c_2`$ and $`f_k \in C_b(\mathcal{X})`$ (locally constant on $`B_1(o)`$), every subsequential weak limit $`\pi`$ from `operator.md` §2 satisfies $`\pi(D_o \ge k) = c_k/c_2`$ and $`\pi(D_o = k) = p_k`$.

### 2.2 `round2/degree-law.md`

- **Birth times & eq. (2) (lines 11–22):** Splitting $`H_o(t)`$ into $`o`$ and vertices promoted at times $`\tau \in (0, t]`$ and compensating promotion clocks yields $`Q_t f = S_t f + \int_0^t Q_s B S_{t-s} f\,ds`$ (passes from finite trees to $`\mathcal{X}`$ by continuity of $`Q_s, S_s, B`$ on $`C_b(\mathcal{X})`$ and Yule domination). Applying $`\pi`$, multiplying by $`e^{-\mu t}`$, and sending $`t \to \infty`$ kills $`e^{-\mu t}\pi(S_t f)`$ since $`|\pi(S_t f)| \le \|f\|_\infty`$ and $`\mu > 0`$, yielding $`\pi f = \int_0^\infty e^{-\mu a}\pi B S_a f\,da`$, $`\pi B 1 = \mu`$, and newborn law $`\beta = \pi B / \mu`$.
- **Supermartingale & Yule $`L^2`$ bound (line 25):** $`D(t)`$ rises at rate $`L + \sum_{w\sim o} 1/D_w \le D`$, so $`e^{-t}D(t)`$ is a nonnegative supermartingale converging a.s. to $`W`$, coupled below a rate-1 Yule process $`Y(t)`$ from $`Y(0)=2`$. Its normalized second moment is $`\mathbb{E}[(e^{-t}Y(t))^2] = 6 - 2e^{-t} \le 6`$; Doob's $`L^2`$ maximal inequality gives $`\mathbb{E}[\sup_t (e^{-t}Y(t))^2] \le 24`$, dominating $`\sup_t (e^{-t}D(t))^\mu`$ (since $`\mu \le 1 < 2`$) and giving uniform integrability of $`e^{-t}D(t)`$.
- **Positivity $`A = \mathbb{E}_\beta W^\mu > 0`$ (line 27):** Silencing $`o`$ on $`[0, 1]`$ and requiring $`\ge 2`$ rings of its tagged leaf has uniform probability $`e^{-1}(1 - 2e^{-1}) > 0`$ and gives $`L(1) \ge 3`$ (with $`D(1) < \infty`$ a.s.). From any state, the leaf drift satisfies $`L + \rho - b = (L - 1) + \rho + (1 - b) \ge L - 1`$ **even when $`L = 0`$** (where $`b = 0`$ and $`\rho \ge 0 > -1`$). Hence $`e^{-t}(L(1+t)-1)`$ is a submartingale (integrable by Yule domination from $`D(1)`$), giving $`\mathbb{E}[e^{-t} D(1+t) \mid \mathcal{F}_1] \ge \mathbb{E}[e^{-t}(L(1+t)-1) \mid \mathcal{F}_1] \ge L(1) - 1 \ge 2`$. Conditional UI from the Yule $`L^2`$ bound from $`D(1)`$ yields $`e\,\mathbb{E}[W \mid \mathcal{F}_1] \ge 2`$ on $`\{L(1) \ge 3\}`$, hence $`\mathbb{E}_\beta[W] > 0`$ and $`A = \mathbb{E}_\beta[W^\mu] > 0`$. *(Minor expositional note: $`\ell \ge 2`$ already gives $`\ell - 1 \ge 1 > 0`$; $`\ell \ge 3`$ is stronger and equally valid.)*
- **Tail asymptotic (line 29):** By (2) and right-continuous monotonicity of $`D(t)`$ with unit jumps, $`\pi(D_o \ge k) = \mathbb{E}_\beta[e^{-\mu T_k}; T_k < \infty]`$. On $`\{\sup D = \infty\}`$, $`T_k \uparrow \infty`$ (nonexplosion via $`D \le Y`$) and $`k e^{-T_k} = D(T_k)e^{-T_k} \to W`$; on $`\{\sup D < \infty\}`$, $`(k e^{-T_k})^\mu \mathbf{1}_{T_k < \infty} = 0 = W^\mu`$ eventually. Since $`(k e^{-T_k})^\mu \mathbf{1}_{T_k < \infty} \le \sup_t (e^{-t} D(t))^\mu \in L^1`$, DCT gives $`\pi(D_o \ge k) \sim A k^{-\mu}`$.
- **Point probabilities (lines 31–37):** For $`k \ge 2`$, newborns have degree $`2 < k+1`$, and cell counts of $`\{D_u \ge k+1\}`$ jump by $`+1`$ solely on leaf gains at degree-$`k`$ vertices (intensity $`k - c_u + \rho_u \le k`$ per cell member, bounded by $`k|H_o(s)|`$). Compensation on $`[0, 1]`$ and integration against $`\pi`$ yield $`\mu \pi(D_o \ge k+1) = \pi[(k - c_o + \rho_o)\mathbf{1}_{D_o = k}]`$. Using $`0 \le \rho_o \le c_o/2`$ gives upper bound $`k \pi(D_o = k)`$ and lower bound $`(1-\delta)k[\pi(D_o=k) - \pi(c_o > \delta k)]`$; Markov + $`\pi(e^{\theta c_o}) < \infty`$ (`operator.md` (5)) makes $`\pi(c_o > \delta k) = O(e^{-\theta\delta k}) = o(k^{-1-\mu})`$. Taking $`k \to \infty`$ then $`\delta \downarrow 0`$ yields $`p_k = \pi(D_o = k) \sim \mu A k^{-1-\mu}`$.

---

## 3. Independent Recomputation & Planted-Error Evidence

Command run:
```bash
python3 -u -B round2/outputs/attack_astra_34/audit_claim3.py
```

Exact output:
```text
=== INDEPENDENT ADVERSARIAL AUDIT OF ASTRA CLAIM 3 ===

[1] Exact small-n enumeration (n = 3..10):
  n | E[C_n] = E[N_{>=2}]          | E[N_{=2}]                    | E[N_{>=3}]
  3 | 1                            | 1                            | 0
  4 | 4/3                          | 2/3                          | 2/3
  5 | 19/12                        | 2/3                          | 11/12
  6 | 647/360                      | 253/360                      | 197/180
  7 | 17189/8640                   | 9709/12960                   | 32149/25920
  8 | 4717423/2177280              | 1739489/2177280              | 1488967/1088640
  9 | 101568829/43545600           | 42281369/49766400            | 172193683/116121600
 10 | 200684675729/80621568000     | 5079690167/5643509760        | 33215693089/20901888000
  PASS: Exact E[C_n] matches n=3..9 and exact telescoping sum_{k>=2} E[N_{=k}(n)] = E[C_n] holds.

[2] Core verification checks:
  [PASS] Yule Laplace & T_2<->T_3 factor 2/(2+mu): Exact Yule Taylor identity (to order 8) and T_2/T_3 factor 2/(2+z) verified.
  [PASS] Crossing identity (degree-law.md eq. 35): Exact crossing identity n*(E[N_{>=k+1}(n+1)] - E[N_{>=k+1}(n)]) == E[sum_{D_u=k} (k - c_u + rho_u)] verified for all 3<=n<=9, 2<=k<=n.
  [PASS] Small-t resolvent / generation expansion (eq. 3): Exact Taylor/generation expansion of asymptotics.md eq. (3) verified for k=2,3,4,5 up to order t^4/4! (matching O(z^{-5})).
  [PASS] Fixed-size Yule Tauberian transfer (Gamma(mu+2)): Fixed-size Yule Tauberian transfer verified on 3 synthetic power-law+correction sequences; confirms c_k = d_k^* / Gamma(mu+2).
  [PASS] Spectral pole, residue (4), aperiodicity: Spectral pole, lambda'(mu)=ell(K'_mu g)<0, residue (4), and aperiodicity verified.
  [PASS] Degree-law Yule L^2 & worst-case leaf drift: Yule L^2 moment 6-2e^(-t) exact; worst-case rho=0 leaf drift at t=1.2: E[L(t)-1]=7.224508 >= 2e^t=6.640234 (with P(L=0)=0.00411 > 0).

[3] Planted-error negative controls (all MUST be caught = rejected):
  [CAUGHT] Planted Error 1: Gamma(mu+1) instead of Gamma(mu+2) -> Yule limit mismatch for mu=0.565723: E(t)e^(-mu t)=1.971838 vs c*Gamma(mu+1.0)=1.263965 (rel_err=0.5600)
  [CAUGHT] Planted Error 2: Wrong T_2<->T_3 factor 1/(1+mu) instead of 2/(2+mu) -> T2/T3 Laplace ratio mismatch at z=5657/10000, n=3: 200000000/914851649 != 10000/15657*10000/35657
  [CAUGHT] Planted Error 3: Drop rho_o in crossing identity -> Crossing identity failed at n=4, k=2: lhs=1 != rhs=2/3
  [CAUGHT] Planted Error 4: Replace k-c_o by k-1 in crossing identity -> Crossing identity failed at n=3, k=2: lhs=2 != rhs=1
  [CAUGHT] Planted Error 5: Replace (n-1) by n in Yule Laplace identity -> Yule Taylor mismatch at n=2, order m=0: 2 != 1
  [CAUGHT] Planted Error 6: Wrong residue sign +ell(K'_mu g) -> Residue mismatch: contour=(0.6484648347585015-4.2814127351915865e-17j) vs formula=-0.6484648347575035
  [CAUGHT] Planted Error 7: Lattice kernel (peripheral pole at mu+2pi i) -> Peripheral pole detected at tau=6.283185307179586: |det(I-K)|=2.4535904415291853e-16

ALL CLAIM 3 CHECKS PASSED & ALL 7 PLANTED ERRORS CAUGHT.
```

Also replayed Astra's own checker read-only (`python3 -u -B round2/checks/asymptotics_carbon.py`): exited 0 with `ALL CHECKS AND 8/8 NEGATIVE CONTROLS PASSED.`
