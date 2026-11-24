VERDICT: HOLDS

# Adversarial Audit of Astra Claim 2 — Renewal Characterization & Genuine Power Law (`ft2-astra@bfcbbe5`)

- **Auditor:** Independent adversarial child of Neon lead (Gemini 4 Carbon (High), distinct model family from GPT-6.2 Astra)
- **Target audited:** `round2/fringe/renewal.md` (with replay check of `round2/checks/fringe_neon.py`)
- **Fresh independent verification script:** `round2/outputs/attack_astra_12/audit_claim2.py`

---

## 1. Independent Recomputation & Exact Generator Verification

A fresh stdlib-only verifier (`audit_claim2.py`) was written directly from the definitions and equations in `round2/fringe/renewal.md` without copying code from Astra's scripts.

1. **Exact friend-tree enumeration ($`n = 3, \dots,9`$) and forward-backward generator identity:**
   - Propagated both explicit labelled undirected trees from $`T_3`$ ($`n = 3, \dots,7`$, verifying exact pushforward agreement) and canonical birth-marked rooted trees up to $`n = 9`$ ($`1, 2, 6, 18, 61, 208, 743`$ states; $`1{,}039`$ distinct marked rooted trees total).
   - Verified the exact rational mean core counts $`\mathbb E[C_n]`$ for $`n = 3, \dots,9`$:
     $`1, \; 4/3, \; 19/12, \; 647/360, \; 17189/8640, \; 4717423/2177280, \; 101568829/43545600`$, matching FACT 2 of `astra.md` identically.
   - Verified the exact generator identity $`(\mathcal L^{\mathrm{phys}} \Phi_f)(T) = \sum_{v\ne o}((R+B)f)(S_v) + I(T; f)`$ across all $`1{,}039`$ states as formal linear combinations of canonical types, and verified continuous-time derivatives $`\frac{d^k}{dt^k}\big|_{t=0}`$ of $`M'_t f = M_t(R+B)f + I_t f`$ for orders $`k = 1, \dots,5`$ on four test functions (including non-linear functions of $`a_S`$, $`b_S`$, subtree sizes, leaf counts, and child marks).
2. **Exact rational Lyapunov verification (§3):**
   - Verified inequality (4) in `fractions.Fraction` across all $`31{,}372`$ pairs with $`2 \le D < 250`$ and $`0 \le c \le D`$ (tightest pair $`(D, c)=(3,1)`$, where $`\frac{2}{3}(-\frac{3}{20} + \frac{20}{29}) = \frac{313}{870} \le \frac{9}{25} = \frac{313.2}{870}`$, slack $`\frac{1}{4350}`$), plus the analytic tail bound for $`D \ge 100`$: $`\mathrm{LHS} \le (1 + 100\cdot\frac{3}{20})^{-1} = \frac{1}{16} = 0.0625 < q\min h = \frac{63}{250} = 0.252`$.
   - Verified $`\sum_Y a_Y W(Y) = (a_S + d)W(S) - h(c_o)`$, $`W(Y^+) = W(Y)`$, $`L^0 W \le q W`$, $`R W = L^0 W + d W - h(c_o) \le \rho W`$ ($`\rho = 51/100`$), $`W(S) \ge 14/13 > 1`$, $`W(S) \le |S|/d`$, and $`B W \le 1/d`$ on all $`947`$ reachable and $`a_S = 0`$ limit types.
3. **Planted-error falsification tests:**
   - Drop $`(a_Y - p_S)f(Y)`$ copy term in $`R`$: caught at $`n = 6`$.
   - Drop $`-a_S f(S)`$ killing term in $`R`$: caught at $`n = 4`$.
   - Flip copy-term sign in $`R`$ or $`I_t`$: caught at $`n = 5`$ and $`n = 4`$.
   - Perturb $`q = 9/25 \to 7/20`$: caught at $`(D, c)=(3,1)`$ ($`313/870 > 7/20`$).
   - Perturb $`h(1) = 1 \to 9/10`$ or $`h(2) = 17/20 \to 19/20`$: caught at $`(D, c)=(2,1)`$ and $`(D, c)=(3,1)`$.

---

## 2. Line-by-Line Proof Audit of `round2/fringe/renewal.md`

- **§1–§2 (Mark inequality, telescoping transfer, and rigor of (3)):**
  - When $`u`$ promotes a leaf $`v`$, $`D_u \ge 2`$ is unchanged at promotion and nondecreasing thereafter, so $`p_u(t) = 1/D_u(t) \le 1/D_u(\text{birth of } v) = a_v \le 1/2`$ holds permanently.
  - Every internal transition inside $`S_v`$ ($`v\ne o`$) is counted in $`(L^0 f)(S_v)`$ for every ancestor $`o \prec w \preceq v`$, promotions at $`v`$ add $`b_v f(J_{p_v}) = (Bf)(S_v)`$, and parent-to-child leaf gains $`u \to x`$ ($`p_u[f(S_x^+) - f(S_x)]`$) rewrite via $`\pm a_x f(S_x)`$ into parent copy births $`(a_x - p_u)f(S_x) + p_u f(S_x^+)`$ (with nonnegative rates!) plus child killing $`-a_x f(S_x)`$, matching (1)–(2) with zero discrepancy.
  - In §3 line 59, for $`|f| \le C W`$, conditioning on the first jump of $`S`$ (finitely many transitions, total rate $`\le \frac{3}{2}|S| + \frac{3}{2}`$, one-step successor weight sum $`O(|S|+1)`$) proves continuity and continuous differentiability of $`s \mapsto P_s f(S)`$ with $`\partial_s P_s f = (R+B)P_s f`$ and $`|\partial_s P_s f(S)| \le C_t(1+|S|)^2`$. Applying Dynkin's jump formula to $`\sum_{v\ne o} (P_{t-u}f)(S_v(u))`$ on $`[0, t]`$ is fully integrable because the Yule vertex count $`N(t) \sim \mathrm{NegBin}(3, e^{-t})`$ has finite moments of all polynomial orders, rigorously establishing (3).
- **§3 (Elementary $`W`$ bounds, nonexplosion, and (6)):**
  - Since every core $`x`$ in an auxiliary type $`S`$ has $`c_x \ge 1`$ (the root of $`S`$ retains its parent slot!), $`7/10 \le h(c_x) \le h(1) = 1`$ throughout $`S`$ (whereas $`h(0)=2`$ occurs only at the full-tree centre $`o`$). With $`0 \le a_S \le 1/2`$ and $`z(a) \le 1`$, this gives $`W(S) \ge (7/10)/(1/2 + 3/20) = 14/13 > 1`$, $`W(S) \le \mathrm{cores}(S)/d \le |S|/d`$, and $`BW(S) = b_S/(p_S + d) < 1/d`$.
  - Stopping at particle-count and event cutoffs, $`W \ge 14/13 > 1`$ bounds alive particle count $`N_{\mathrm{alive}}`$ in probability via the discounted supermartingale; while $`N_{\mathrm{alive}} < K`$, death intensity is $`< K/2`$ so total deaths $`D_{\mathrm{deaths}}(t) < \infty`$ a.s., particle conservation $`B_{\mathrm{births}} = N_{\mathrm{alive}} - 1 + D_{\mathrm{deaths}}`$ bounds total births, and each born particle's $`L^0`$ vertex count is Yule-dominated, ruling out explosion. Fatou gives $`R_t W \le e^{\rho t}W`$ and $`P_t W \le e^{(\rho+1/d)t}W`$; on the full tree $`Z(T_3) = h(0) = 2`$ and $`\mathcal L^{\mathrm{phys}} Z \le q Z`$ give (6).
- **§4 (Compactness at $`a=0`$, strong positivity, Krein–Rutman, and dichotomy (9)):**
  - $`\mathcal A = \{1/2, 1/3, \dots\} \cup \{0\}`$ is compact. In $`J_a`$, $`L^0`$, $`b_S`$, $`p_S`$, and all child marks $`a_Y`$ are completely independent of the root mark $`a`$ (verified symbolically to order $`4`$ in `audit_claim2.py`); $`a`$ enters (8) solely through the outer killing factor $`e^{-(s+a)t}`$. For $`\|g\|_\infty \le 1`$, the $`a`$-independent expectation $`\Phi_{s, g}(t)`$ in (8) satisfies $`|\Phi_{s, g}(t)| \le 1 + \mathbb E^{\mathrm{int}}_{J_*}[\widetilde W(S_t)]/(s-\rho) \le 1 + e^{qt}/(s-\rho)`$ since $`\widetilde W(J_*) = h(1) = 1`$. Because $`s > \rho = q + d > q`$, $`\int_0^\infty t e^{-st}(1 + e^{qt}/(s-\rho))\,dt < \infty`$, giving a uniform Lipschitz bound in $`a \in \mathcal A`$ (including at $`a=0`$) across $`\|g\|_\infty \le 1`$; Arzelà–Ascoli yields compactness of $`K_s`$ on $`C(\mathcal A)`$.
  - If $`g \in C(\mathcal A)_+ \setminus \{0\}`$, continuity at $`0`$ implies $`g(1/j) > 0`$ for some $`j \ge 2`$ (even if $`g`$ was constructed to peak around $`0`$). From $`J_a`$, survival on $`[0,2]`$ ($`\ge e^{-1}`$), root silence on $`[0,1]`$ ($`e^{-1}`$), $`j-2`$ rings of the initial leaf with silent new leaves on $`[0,1]`$, and total silence on $`[1,2]`$ ($`e^{-j}`$) has positive $`a`$-uniform probability and holds $`p_S = 1/j, b_S = (j-1)/j \ge 1/2`$ throughout $`[1,2]`$, so $`K_s g \in \mathrm{int}(C(\mathcal A)_+)`$.
  - Strong Krein–Rutman gives simple $`r(K_s) > 0`$ with $`h_s \in \mathrm{int}(C(\mathcal A)_+)`$ and strictly positive dual eigenmeasure $`\nu_s`$; pairing $`(K_{s_1}-K_{s_2})h_{s_2} \ge \delta h_{s_2}`$ with $`\nu_{s_1}`$ proves strict decrease of $`r(K_s)`$, and norm continuity of $`s \mapsto K_s`$ gives continuity.
  - Tonelli's theorem justifies the generation expansion (9). Since $`g_s(J_{(\cdot)}) \ge (s+1/2)^{-1} > 0`$ and $`g_s(J_{(\cdot)}) \in C(\mathcal A)`$, $`r(K_s) < 1`$ implies norm convergence of (9), $`u_s \le C_s W`$, and $`\int_0^\infty e^{-st} M_t 1\,dt < \infty`$ by (3) and (6); whereas $`r(K_s) \ge 1`$ implies $`g_s(J_{(\cdot)}) \ge c h_s \implies u_s(J_a) = \infty`$ for all $`a \in \mathcal A \implies \int_0^\infty e^{-st} M_t 1\,dt = \infty`$. Since $`\mu \ge 113/200 > 51/100 = \rho`$, $`r(K_s) = 1`$ has the unique root $`s = \mu`$ on $`(\rho, \infty)`$.
- **§5 (Volterra eigenfunction, two-sided $`C`$ comparison, and fixed-size transfer):**
  - With $`K_\mu g = g > 0`$ and $`H = U_\mu g \ge 0`$, $`(BH)(S) = b_S H(J_{p_S}) = b_S g(p_S)`$, so Tonelli + $`R_t R_u = R_{t+u}`$ directly give $`R_t H = \int_t^\infty e^{-\mu(v-t)} R_v(BH)\,dv = e^{\mu t}H - \int_0^t R_s B(e^{\mu(t-s)}H)\,ds`$. First-direct-birth decomposition gives $`P_t H = R_t H + \int_0^t R_s B P_{t-s}H\,ds`$, and Grönwall in $`\|\cdot\|_W`$ (using $`\|R_s B\|_W \le d^{-1}e^{\rho s}`$) uniquely forces $`P_t H = e^{\mu t}H`$, giving $`e^{-\mu t}M_t H \uparrow A_0 \in (0, \infty)`$.
  - Upper bound: $`z(a) \le z(1/2) = 10/13 < 1`$ and each core $`x\ne o`$ has at most one ancestor $`v\ne o`$ at each distance $`k \ge 0`$, so $`\sum_{v\ne o} W(S_v) \le \frac{20}{3}(1 - \frac{10}{13})^{-1}(C-1) = \frac{260}{9}(C-1)`$ and $`\sum_{v\ne o} H(S_v) \le C_1 C`$.
  - Lower bound: For noninitial $`v \in G_1`$ ($`c_v \le 2, L_v \ge 1`$), root survival on $`[0,1]`$ ($`\ge e^{-1/2}`$) and $`\mathbb E^{\mathrm{int}}\int_0^1 b_{S_t}\,dt = \mathbb E^{\mathrm{int}} N_{\mathrm{prom}}(1) \ge (1-e^{-1})/3`$ give $`H(S_v) \ge e^{-(\mu+1/2)}(1-e^{-1})\min g / 3 = c_1 > 0`$. For noninitial $`v \in G_2`$ (bare $`c_v=D_v=2`$ with a bare child $`w`$), $`w`$'s other neighbor $`x`$ is a core child of $`w`$ inside $`S_v`$; silencing $`v`$ and $`x`$ on $`[0,1]`$ keeps $`D_w = 2`$ up to $`w`$'s first ring (incoming rings into silent $`x`$ do not reach $`w`$), so $`w`$'s first ring on $`[0,1]`$ gives $`v`$ a leaf with probability $`\frac{1}{2}(1-e^{-1})`$, and $`v`$'s silence on $`[0,1]`$ preserves $`c_o(S_1)=2, L_o(S_1)\ge 1`$ up to $`t=1`$, giving $`H(S_v) \ge \frac{1}{2}e^{-(\mu+5/2)}(1-e^{-1}) c_1 > 0`$. By `existence.md` Lemma 2.1, $`|\mathcal E| \ge C/4`$; every $`u \in \mathcal E \setminus (G_1 \cup G_2)`$ has $`u \ne o`$ and bare parent $`p(u) \in G_2`$, with $`|p^{-1}(w)| \le 1`$ ($`w\ne o`$) and $`|p^{-1}(o)| \le 2`$, so $`|\mathcal E| \le |G_1| + 3|G_2| \le 3|G_1 \cup G_2| \implies |(G_1 \cup G_2)\setminus\{o\}| \ge C/12 - 1`$, proving $`C \le C_2(1 + \sum_{v\ne o} H(S_v))`$ for all $`C \ge 1`$.
  - Yule transfer: $`\Pr(C_{j+1}-C_j = 1 \mid T_j) = j^{-1}\sum_{u\in K} L_u/D_u \le C_j/j`$ gives $`a_m \le (m/n)a_n`$ ($`m \ge n \ge 3`$). At $`t = \log n`$, $`\Pr(N(t) \ge n) \ge (1 - 1/n)^{n-1} > e^{-1}`$ and $`\mathbb E\max(1, N(t)/n) \le 1 + \mathbb E N(t)/n = 3`$, yielding $`e^{-1} a_n \le F(\log n) \le 3 a_n`$ and hence $`\mathbb E C_n \asymp n^\mu`$.

---

## 3. Exact Commands & Outputs

```text
$ python3 -B round2/outputs/attack_astra_12/audit_claim2.py
=== [Claim 2 / Check 1A] Exact Friend-Tree Enumeration (n=3..9) & M'_t f Identity ===
  n=3:    1 marked states, exact E[C_3] = 1 (matches FACT 2)
  n=4:    2 marked states, exact E[C_4] = 4/3 (matches FACT 2)
  n=5:    6 marked states, exact E[C_5] = 19/12 (matches FACT 2)
  n=6:   18 marked states, exact E[C_6] = 647/360 (matches FACT 2)
  n=7:   61 marked states, exact E[C_7] = 17189/8640 (matches FACT 2)
  n=8:  208 marked states, exact E[C_8] = 4717423/2177280 (matches FACT 2)
  n=9:  743 marked states, exact E[C_9] = 101568829/43545600 (matches FACT 2)
  Exact generator & time-derivative identities M'_t f = M_t(R+B)f + I_t f verified across 1039 states and orders k=1..5 on 4 test functions.
  [PASS] Check 1A

=== [Claim 2 / Check 1B & 2] Rational Inequality (4), Child Identity & §4-§5 Bounds ===
  Inequality (4) verified for 31372 pairs (2<=D<250, 0<=c<=D); tightest at (D=3, c=1): 313/870 <= 9/25 (slack=1/4350).
  Analytic tail D>=100: LHS <= 1/(1+100*3/20) = 1/16 = 0.0625 <= q*min(h) = 63/250 = 0.2520.
  Child-weight identity & RW <= (51/100)W verified on 947 types (including a_S=0).
  §4 Spine decomposition (8) verified completely independent of root mark a in {1/2, 1/3, 1/20, 0}.
  §5 Upper W-sum, lower (C/12 - 1) good-vertex bound, and Yule transfer [e^{-1}, 3] verified.
  [PASS] Check 2

=== [Claim 2 / Check 3] Planted-Error Falsification Tests ===
  [CAUGHT] Planted Error 1 (drop (a_Y - p_S)f(Y) copy term): fails at n=6
  [CAUGHT] Planted Error 2 (drop -a_S f(S) killing term): fails at n=4
  [CAUGHT] Planted Error 3 (flip sign in R / I_t): fails at n=5 / n=4
  [CAUGHT] Planted Error 4a (perturb q=9/25 -> 7/20): fails at D=3, c=1 (313/870 > 7/20)
  [CAUGHT] Planted Error 4b (perturb h(1)=9/10): fails at D=2, c=1 (187/520 > 81/250)
  [CAUGHT] Planted Error 4b (perturb h(2)=19/20): fails at D=3, c=1 (371/870 > 9/25)
  [PASS] All planted-error tests bit cleanly.

ALL ADVERSARIAL CHECKS FOR CLAIM 2 COMPLETED SUCCESSFULLY.
```

```text
$ python3 -B round2/checks/fringe_neon.py
...
ALL INDEPENDENT CHECKS PASSED in 2.27s.
```
