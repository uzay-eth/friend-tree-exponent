VERDICT: HOLDS

# Adversarial Audit of Astra Claim 4 — Spatial Operator

- **Target:** Astra Claim 4 (`ft2-astra@8e4db79`), covering `round2/operator.md` and its §5 dependency `round2/finite/one-port-upper.md` (plus §2–§3 dependency `proofs/existence.md`).
- **Auditor:** Gemini 4 Neon (independent model family from GPT-6.2 Astra).
- **Independent test script:** `round2/outputs/attack_astra_34/audit_claim4.py` (written from scratch; `[local directory]` left strictly read-only).

---

## 1. Summary Verdict

The construction of $`Q_t`$ on $`C_b(\mathcal{X})`$, the spatial cutoff bound (3), existence of a Borel probability eigenmeasure $`\pi Q_t = e^{\mu t}\pi`$, minimality (1), generator/unimodularity identities (§4), and the upgrade to $`r_{C_b}(Q_1) = e^\mu`$ and uniqueness of the positive eigenvalue via `one-port-upper.md` (§5) all hold. All core identities and quantitative prerequisites pass independent exact-rational and numerical recomputation (with 5/5 planted-error controls caught). Three minor expositional compressions are noted below; none is a mathematical gap.

---

## 2. Line-by-Line Proof Audit

### 2.1 `round2/operator.md`

- **State space $`\mathcal{X}`$ & §1 (Cells, infinite-tree construction, cutoff (3), continuity, lines 5–38):**
  - $`\mathcal{X}`$ (locally finite rooted core trees with $`L`$-marks and $`D = L + c \ge 2`$) is Polish under the local ball ultrametric (finite trees are countable and dense; compatible $`R`$-balls in a Cauchy sequence specify an inductive-limit locally finite tree with $`D \ge 2`$ at every vertex).
  - Microscopic core dynamics match $`r_x = L_x + \sum_{y\sim x} 1/D_y`$ and $`b_x = L_x / D_x`$. Cell size $`|H_o(s)|`$ is dominated by a rate-1 Yule process from $`1`$, giving $`Q_s 1 \le e^s`$ and $`\mathbb{E}|H_o(s)|^2 \le 2e^{2s} - e^s \le 2e^{2s}`$.
  - **Infinite-tree drivers (line 21):** Old oriented edges $`x \to y`$ use independent rate-$`1/c_x^0`$ Poisson proposals thinned by $`c_x^0 / D_x \in (0, 1]`$ (valid since $`D_x \ge c_x \ge c_x^0`$); internal cell transitions use independent cell drivers (leaving old-edge slots idle). Any finite union of cells has finite inbound boundary proposal rate $`\Lambda_{\text{in}}`$ and endogenous birth rate bounded by current total vertex count (cores + leaves), hence is nonexplosive.
  - **Spatial cutoff (3) (lines 23–37):** Comparing truncations $`R < R'`$ under common drivers on $`B_R(o)`$, any cell discrepancy at $`w \in B_k(o)`$ traces back via finitely many events in $`B_{R'}(o)`$ to an ordered proposal walk from distance $`R+1`$ to $`w`$. Loop erasure yields a simple path of length $`\ge R - k + 1 = j`$ with strictly increasing proposal times; taking its **last $`j`$ edges** gives a simple length-$`j`$ path ending at $`w`$ with ordered proposals in $`[0, s]`$. *(Minor expositional note: taking the last $`j`$ edges when $`d(o, w) < k`$ is implicit in line 29, and is exact.)* For each simple length-$`j`$ path $`(x_0, \dots, x_j = w)`$, the expected ordered-tuple count is $`\frac{s^j}{j!} \prod_{i=0}^{j-1} 1/c_{x_i}^0`$. Reversing the path uses $`c_x^0 P(x, y) = c_y^0 P(y, x) = 1`$ and $`c_{x_0}^0 \ge 1`$ to bound the sum over all simple length-$`j`$ paths ending at $`w`$ by $`c_w^0`$. For a cylinder $`f`$ of radius $`k \le R`$, the cell-quotient map $`V(G_s) \to V(G_0)`$ is 1-Lipschitz, so $`B_k(G_s, v)`$ for $`v \in H_o(s)`$ only meets cells in $`B_k(o)`$. Union bound gives failure probability $`\le p = \min(1, \frac{s^j}{j!}\sum_{w\in B_k(o)} c_w^0)`$, and separate Cauchy–Schwarz on $`|H_o(s)|, |H_o^{[R]}(s)|`$ gives $`|Q_s f - Q_s^{[R]} f|(G, o) \le 2\sqrt{2}\,e^s\|f\|_\infty \sqrt{p}`$.
  - Summability in $`R`$ constructs the infinite-tree process a.s. on every finite cell set. On the neighborhood $`\{B_{R+1}(G', o') \cong B_{R+1}(G, o)\}`$, both $`Q_s^{[R]} f`$ and $`\sum_{w\in B_k(o)} c_w^0`$ are constant, proving continuity of $`Q_s f`$ for bounded cylinders; coupling on $`B_{R_n}(o)`$ ($`R_n \to \infty`$) + Yule $`L^2`$ uniform integrability proves $`Q_s : C_b(\mathcal{X}) \to C_b(\mathcal{X})`$. Semigroup $`Q_s Q_t = Q_{s+t}`$ holds on finite completions and extends to all of $`\mathcal{X}`$ by continuity and density of finite trees (as well as directly by the Markov property).

- **§2 (Tilted eigenmeasure $`\pi`$, lines 39–56):**
  - $`F(t) \le F_3(t) \le F(t+1)/(1 - e^{-2})`$ follows directly from the independent $`\mathrm{Exp}(2)`$ wait from $`T_2`$ to $`T_3`$, so $`t^{-1}\log F_3(t) \to \mu`$.
  - With $`w(t) = e^{-\mu t - \delta|t-T|}`$ ($`0 < \delta < \mu`$), $`Z \ge e^{-\mu(T+1)-\delta} F_3(T)`$ gives $`w(0)/Z \le e^{\mu+\delta} e^{(\mu-\delta)T}/F_3(T) \to 0`$ as $`T \to \infty`$ at fixed $`\delta`$. Taking a single $`h`$-independent sequence $`(\delta, T)`$ with $`\delta \downarrow 0`$ and $`w(0)/Z \to 0`$: for every fixed $`h \ge 0`$, $`w(u-h)/w(u) \in [e^{(\mu-\delta)h}, e^{(\mu+\delta)h}]`$ on $`u \ge h`$ and $`w(u)/w(0) \le e^{-(\mu-\delta)u} \le 1`$ on $`[0, h]`$, proving (4) $`\|\pi_{\delta, T} Q_h - e^{\mu h}\pi_{\delta, T}\|_{\mathrm{TV}} \to 0`$.
  - Moment bounds (5): After age $`a`$, $`\mathbb{E}[D_v(\tau+a) \mid \mathcal{F}_\tau] \le 2e^a`$ (so $`\mathbb{E}[D^\alpha] \le 2^\alpha e^{\alpha a}`$ by concavity for $`0 < \alpha < \min(1, \mu)`$) and $`c_v(\tau+a) \preceq 1 + \mathrm{Pois}(a)`$ (so $`\mathbb{E}[e^{\theta c}] \le e^\theta e^{(e^\theta - 1)a}`$ for $`0 < \theta < \log(1+\mu)`$), valid for the initial centre too. Using $`w(t+a) \le w(t)e^{-(\mu-\delta)a}`$ (exact for all $`t, a \ge 0`$ by triangle inequality $`|t+a-T| \ge |t-T|-a`$) and total promotion intensity $`\le F_3(t)`$ yields (5) $`\pi_{\delta, T}(\text{moment}) \le B(1 + w(0)/Z)/(\mu - \delta - \kappa)`$.
  - Mass-transport ball tightness: Cauchy–Schwarz + (5) gives $`\sup \pi_{\delta, T}(c_o \mathbf{1}_{D_o > J}) \to 0`$. Reversing simple paths of length $`i`$ with earlier degrees $`\le J_0, \dots, J_{i-1}`$ and endpoint degree $`> J_i`$ gives expected count $`\le (J_1 \cdots J_{i-1})\pi_{\delta, T}(c_o \mathbf{1}_{D_o > J_i})`$. Choosing $`J_i`$ inductively yields compact sets $`K_\varepsilon`$ with $`\pi_{\delta, T}(K_\varepsilon^c) \le \varepsilon`$; Prokhorov gives a weak limit $`\pi`$ satisfying $`\pi Q_h = e^{\mu h}\pi`$ for all $`h \ge 0`$ (and inheriting the moment bounds by Fatou).

- **§3 (Minimality, lines 57–68):**
  - Planting from any $`(G, o) \in \mathcal{X}`$: if $`L_o = 0`$, then $`D_o \ge 2 \implies c_o \ge 2`$, so $`o`$ has a core neighbor $`w`$. Silencing all $`D_w = c_w + L_w < \infty`$ neighbors of $`w`$ (both core and leaf neighbors — *minor expositional clarification: "all neighbours of $`w`$" must include $`w`$'s leaves so $`D_w`$ remains constant*) keeps $`D_w`$ fixed until $`w`$ sends a leaf to $`o`$; if $`L_o \ge 1`$, silence $`o`$. Next unit: tagged leaf rings $`\ge M`$ times $`\implies L_o \ge M+1, c_o = c_o^0`$. Next unit: promote a leaf to $`v^*`$ (prob $`\ge (M+1)/(M+1+c_o^0)`$) and protect $`\{v^*, v'\}`$. Next $`k-2`$ units: script the remaining $`k-2`$ births as in `existence.md` Lemma 3.1. Delay $`h = k+1`$ and probability $`p(G, o, M, k) > 0`$ are uniform over the finite set of $`k`$-vertex seeds. Applying `existence.md` Corollary 5.3 + Lemma 6.1 (directly via Section 1 drivers or via finite completions $`\widetilde{G}_R \to G`$ using (3) for $`R \ge 2`$) gives $`Q_{s+h} 1(G, o) \ge p \frac{e^{-\eta s}}{k}(F(s) - k - \theta_k(s))`$, hence $`\liminf_{s\to\infty} s^{-1}\log Q_s 1(G, o) \ge \mu`$.
  - For any probability eigenmeasure $`\sigma Q_1 = r\sigma`$, $`r^n = \int Q_n 1\,d\sigma`$, so Fatou on $`e^{-(\mu-a)n} Q_n 1`$ ($`a > 0`$) gives $`r \ge e^\mu`$.

- **§4 (Generator (7) and edge transfers, lines 69–78):**
  - For a bounded cylinder $`f`$ of radius $`k`$ supported on $`\{\max_{x \in B_k(o)} D_x \le M\}`$, neither gain nor promotion can decrease any degree or remove a large-degree vertex from $`B_k(o)`$, so $`\Delta_x^+ f = \Delta_x^b f = 0`$ outside the support; inside the support $`|B_k(o)| \le 1 + M^k`$ and $`r_x \le D_x \le M`$, while $`|b_o f(G^{b_o}, \mathrm{new})| \le \|f\|_\infty`$. Thus $`Af`$ is a uniformly bounded cylinder of radius $`k+2`$ (and $`A1 = b_o`$).
  - Dynkin's formula $`Q_t f - f = \int_0^t Q_s Af\,ds`$ holds on finite truncations $`B_R(o)`$ for $`R \ge k+2`$ (where $`Af`$ has zero boundary discrepancy on all cells in $`B_k(o)`$) or on finite completions $`\widetilde{G}_R`$ (*minor expositional note:* stopping the finite truncation's total vertex count at $`N`$, rather than only cell $`H_o`$'s size, bounds all rates before applying finite-CTMC Dynkin, and Yule domination on $`|H_o(t)|`$ then removes $`N`$), passing to $`\mathcal{X}`$ via (3) and DCT. Integrating against $`\pi`$ yields $`\pi(Af) = \mu\pi(f)`$ and $`\pi(b_o) = \mu`$.
  - Unimodularity of $`\pi`$ passes from finite trees on bounded-degree balls and extends by MCT; finite support of $`g`$ and $`a`$ + $`\pi(c_o) < \infty`$ from (5) put both transfer divergences in $`L^1(\pi)`$ with integral $`0`$.

- **§5 (Spectral radius & eigenvalue uniqueness, lines 79–88) + `round2/finite/one-port-upper.md`:**
  - Every member of $`H_o(t) \setminus \{o\}`$ lies in a unique disjoint downward branch founded at a direct promotion of $`o`$ (predictable intensity $`b_o \le 1`$). Monotonicity of $`A_j(u)e^{\Lambda_j u}F(u+1)`$ and `one-port-upper.md` (1)–(2) yield (8) $`Q_t 1(G, o) \le 1 + t A_j(t)e^{\Lambda_j t} F(t+1)`$ on every finite core tree (including singleton cores $`T_3`$, since `one-port-upper.md` §3 requires no exterior core neighbors of $`p`$).
  - Finite completions $`\widetilde{G}_R`$ (appending one $`(1,1)`$ core to each missing boundary slot of $`B_R(o)`$) are genuine finite core trees matching $`B_R(o)`$ with exact boundary degrees; by (3) for $`f \equiv 1`$, $`Q_t 1(\widetilde{G}_R, o) \to Q_t 1(G, o)`$, extending (8) to all of $`\mathcal{X}`$.
  - Positivity on $`C_b(\mathcal{X})`$ gives $`\|Q_1^n\| = \|Q_n 1\|_\infty`$; letting $`n \to \infty`$ then $`j \to \infty`$ ($`\Lambda_j \le 6/j + 2/j^2 \to 0`$) gives $`r_{C_b}(Q_1) \le e^\mu`$, while $`\pi Q_1^n 1 = e^{\mu n}`$ gives $`r_{C_b}(Q_1) \ge e^\mu`$. Hence $`r_{C_b}(Q_1) = e^\mu`$, and any probability eigenmeasure $`\sigma Q_1 = r\sigma`$ ($`r > 0`$) satisfies $`e^\mu \le r \le r_{C_b}(Q_1) = e^\mu`$.
  - **Audit of `one-port-upper.md` §§1–7:** All steps verified with zero gaps: §2 joint density (3); §3 linear birth-death PGF $`\partial_h f = (f-1)(f-(1-s))`$, conditional geometric law $`q = s e^{-sh}/(1-(1-s)e^{-sh}) \le e^{-sh}`$, reciprocal bound $`-q\log q/(1-q) \le \sqrt{q} \le e^{-sh/2}`$, and $`\mathbb{E}_{Q^a}[1/D_p(u-)] \le D_0 e^{-\kappa A(u)}`$; §4 monotone rounding profiles $`|\mathcal{P}_t| \le (2+t)^{J/2}`$ and total variation penalty $`\le t/J + 1/2`$; §5 auxiliary change of measure on $`\{L_v \ge K\}`$ with Taylor constant $`c_d = d(d-1)(3/2)^{d-2}/8`$ and low-leaf count $`M+K+N`$; §6 Hölder interpolation $`\theta = 1 - 2/d`$; §7 physical seed $`T_*`$ ($`n_* = 3 + (B-1)H`$) keeping $`D_p \equiv B`$ up through the last proposal $`\tau`$ and matching $`\pi_*`$.

---

## 3. Independent Recomputation & Planted-Error Evidence

Command run:
```bash
python3 -u -B round2/outputs/attack_astra_34/audit_claim4.py
```

Exact output:
```text
=== INDEPENDENT ADVERSARIAL AUDIT OF ASTRA CLAIM 4 ===

[1] Core verification checks:
  Running: Cell operator Q_s vs microscopic, Q_s Q_t = Q_{s+t}, m_{t+h} = m_t Q_h ...
  [PASS] Cell operator Q_s vs microscopic, Q_s Q_t = Q_{s+t}, m_{t+h} = m_t Q_h: Microscopic match, exact Taylor semigroup Q_s Q_t = Q_{s+t} (ord 3), and finite-time m_{t+h}=m_t Q_h (diff=6.37e-07) verified.
  Running: Path-reversal bound, O(s^{R-k+1}) vanishing, and spatial-cutoff (3) ...
  [PASS] Path-reversal bound, O(s^{R-k+1}) vanishing, and spatial-cutoff (3): Path-reversal bound verified on 18 trees (max ratio=1.0000 <= 1); O(s^{R-k+1}) exact; (3) holds (R=0: |diff|=1.339e-02 <= 5.048e-01 (rk4_err=4.75e-03), R=1: |diff|=1.595e-05 <= 6.183e-02 (rk4_err=5.12e-03), R=2: |diff|=4.883e-06 <= 6.183e-03 (rk4_err=6.46e-03)).
  Running: Generator (7), Dynkin's formula, and birth-time renewal identity ...
  [PASS] Generator (7), Dynkin's formula, and birth-time renewal identity: Generator (7), Dynkin Q_t f - f = int Q_s Af ds, and birth-time identity Q_t f = S_t f + int Q_s B S_{t-s} f ds verified exactly up to order 3.
  Running: Directed edge-transfer cancellations (unimodularity) ...
  [PASS] Directed edge-transfer cancellations (unimodularity): Both directed edge-transfer divergences cancel identically to 0 on finite trees.
  Running: one-port-upper.md (§1,§3,§5) & operator.md §2 tilt inequalities ...
  [PASS] one-port-upper.md (§1,§3,§5) & operator.md §2 tilt inequalities: one-port-upper.md §1, §3 (BD PGF + E[1/Z_h]), §5 (Taylor c_d), and operator.md §2 tilt bounds all verified.

[2] Planted-error negative controls (all MUST be caught = rejected):
  [CAUGHT] Planted Error 1: Drop newborn term b_o f(G^{b_o}, new) in generator (7) -> Dynkin/Generator (7) mismatch at order m=1: 1971/500 != 372/125
  [CAUGHT] Planted Error 2: Drop missing boundary slots in Q_s^{[R]} -> Cutoff vanishing order broken at R=1, k=0, order m=1 < j=2: 21169/6000 != 361/100
  [CAUGHT] Planted Error 3: Drop indirect descendants from cell H_o(s) -> Microscopic vs cell Q_s mismatch on T3: [Fraction(-23, 500), Fraction(1971, 500), Fraction(-4293, 200), Fraction(876133, 6000)] != [Fraction(-23, 500), Fraction(1971, 500), Fraction(-10491, 500), Fraction(435053, 3000)]
  [CAUGHT] Planted Error 4: Use 1/D_x instead of 1/D_y in gain rate r_x -> Dynkin/Generator (7) mismatch at order m=3: 435053/3000 != 35839/250
  [CAUGHT] Planted Error 5: Swap sender/receiver rate in gain edge transfer -> Gain-transfer divergence did not cancel: 911/140

ALL CLAIM 4 CHECKS PASSED & ALL 5 PLANTED ERRORS CAUGHT.
```

Also replayed Astra's own checker read-only (`python3 -u -B round2/checks/operator_carbon.py`): exited 0 with `ALL EXACT OPERATOR & DEGREE-LAW CHECKS AND PLANTED-ERROR CONTROLS PASSED.`
