> Packaging note: internal branch and job-runner names on line 8 were redacted for public release; all mathematical content is unchanged.

VERDICT: CORRECT
BOUNDS: [0.543, 0.586] (exact verified prefactors: 0.502326 * Gamma(n+0.543)/Gamma(n) <= E[C_n] <= 2.161878 * Gamma(n+0.586)/Gamma(n) for all n >= 3)

# Independent Cold Audit of the Random Friend Tree Exponent Bounds

Audit conducted strictly cold on an isolated check branch without reading the baseline certificate scripts, the author's independent verifier folder, any other review/check directory, or any other git branch.

---

## 1. Model and Literature Check (`REPORT.md` §§1–2)

### 1.1 Model equivalence
Checked directly against the primary sources:
- **Krapivsky & Redner (2017)** (*Emergent Network Modularity*, J. Stat. Mech. 073405; arXiv:1706.01514 v2): Section 1 defines *isotropic redirection* (IR) at $`r = 1`$ as selecting a provisional target node uniformly at random and attaching the new node to a uniformly random neighbour of that target; Section 2 specifies that "Unless otherwise stated, we assume that the initial network is a dimer" ($`T_2`$).
- **Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou, Lugosi (2024 / 2026)** (*Random Friend Trees*, arXiv:2403.20185 v1; *Random Structures & Algorithms* 69(1), e70082, 2026): Section 1 ("The Model") starts from $`T_2`$ (a single edge joining vertices 1 and 2) and inductively picks $`V_n \in \{1,\dots,n\}`$ uniformly, $`W_n`$ uniformly among neighbours of $`V_n`$, and attaches vertex $`n+1`$ to $`W_n`$.

`REPORT.md` §1 matches both sources identically.

### 1.2 Verification of quoted claims and previous bounds
- **Krapivsky & Redner (2017):** Verified directly in both the arXiv v2 PDF and TeX source (`paper.tex`).
  - Eq. (1): $`N_k \sim N^\mu / k^{1+\mu}`$ ($`k \gg 1`$).
  - Section 5.1, eq. (32) and Fig. 8: $`\mathcal{N} \sim N^\mu`$, $`N_k \sim N^\mu`$, with quoted $`\mu \approx 0.566`$, slope range `0.5652` to `0.5673`, accuracy "within 0.001" (exact quote in `REPORT.md` lines 38–40 is verbatim).
  - Eqs. (34)–(35): the mean-field sum rule $`d\mathcal{N}/dN = (\mathcal{N}/N)\sum_{k\ge 2}\sum_{\ell<k} c_{k,\ell}\,\ell/k`$ and $`\mu = \sum c_{k,\ell}\,\ell/k`$, used non-rigorously to argue $`\mu < 1`$.
  - Section 7: calls determining $`\mu`$ "an outstanding theoretical challenge".
  - Confirmed `"17/30"` does not appear anywhere in `paper.tex`.
- **Addario-Berry et al. (arXiv:2403.20185 v1 and RSA 2026 published PDF):**
  - The main polynomial-bound result is **Theorem 4.7** in arXiv:2403.20185 and renumbered to **Theorem 4.5** in the published RSA 2026 PDF (`REPORT.md` line 60 accurately notes "Thm 4.7 in the arXiv numbering").
  - Both the arXiv version and the published 2026 RSA paper state in the Abstract ("at most $`n^{9/10}`$"), Introduction ("between $`n^{0.1}`$ and $`n^{0.9}`$"), Theorem 4.7 / 4.5 ("$`0.1 < \delta \le \lambda < 0.9`$"), and Open Question 1 ("$`0.1 < \delta < \lambda < 0.9`$") that $`\lambda < 0.9`$.
  - **Verification of the 0.9142 discrepancy:** Checked directly in both versions (Proposition 6.20, Lemma 6.21, Proposition 6.22 — same numbering in both arXiv and RSA 2026):
    - **Proposition 6.20 (lower):** $`\alpha = (\sqrt{13}-3)/2 \approx 0.302776`$ (positive root of $`x = (1-2x)/(1+x)`$), giving $`n^{-\delta} X_n^{\ge 2} \to \infty`$ a.s. for any $`0 < \delta < \alpha/3 = (\sqrt{13}-3)/6 \approx 0.100925 > 0.1`$.
    - **Lemma 6.21 & Proposition 6.22 (upper):** $`\gamma = 3 - 2\sqrt{2} \approx 0.171573`$ (positive root of $`x = (1-5x)/(1-x)`$, matching $`(1-\gamma)/(4\gamma) = 1/(1-\gamma)`$). Proposition 6.22 explicitly sets $`\delta_{\text{upper}} = 1 - \gamma/2 = \sqrt{2} - 1/2 \approx 0.914214`$ and proves $`n^{-\lambda} X_n^{\ge 2} \to 0`$ a.s. for $`\lambda > 1 - \gamma/2`$ ("Omitting the details of the coupling", verbatim in both versions).
    - Thus `REPORT.md` is **correct**: even in the published 2026 RSA version, the written proof in Proposition 6.22 yields only $`\lambda > \sqrt{2} - 1/2 \approx 0.914214`$, not $`\lambda < 0.9`$.
  - Minor rounding note (`REPORT.md` lines 12, 13, 71, 356): writing `[0.1009, 0.9142]` rounds both algebraic endpoints $`(\sqrt{13}-3)/6 \approx 0.100925`$ and $`\sqrt{2}-1/2 \approx 0.914214`$ to 4 decimal places; strictly, the published almost-sure upper bound holds for any $`\lambda > \sqrt{2}-1/2 \approx 0.914214`$.
- **Secondary literature spot-checked (`REPORT.md` §2):**
  - Cannings & Jordan (2013, arXiv:1303.1052): confirms proportion of degree-1 vertices $`\to 1`$ a.s. for 1-step random walk attachment.
  - Krapivsky & Redner (2024, arXiv:2305.10628): confirms $`\nu \approx 1.567`$ ($`C \sim N^{\nu-1}`$) and quote "currently appear to be beyond the reach of available techniques"; no `"17/30"`.
  - Hartle, Krapivsky, Redner, Zhang (2026, arXiv:2604.01540): confirms $`\mu \approx 0.569`$ for IR, exact DAN ($`1-\mu^2 = \int_0^1 x^\mu e^{x-1}dx`$, $`\mu \approx 0.771192`$) and PAN ($`1-\mu = \int_0^1 x^\mu e^{x-1}dx`$, $`\mu \approx 0.549735`$) equations, and quote "the fractions $`q_\ell`$ do not satisfy a closed recurrence"; no `"17/30"`.

### 1.3 Sources read and not read
- **Read:**
  - `https://arxiv.org/abs/1706.01514`, `https://arxiv.org/e-print/1706.01514`, `https://arxiv.org/pdf/1706.01514` (Krapivsky & Redner, 2017; read 2026-11-21 UTC).
  - `https://arxiv.org/abs/2403.20185`, `https://arxiv.org/e-print/2403.20185`, `https://arxiv.org/pdf/2403.20185` (Addario-Berry et al., arXiv v1; read 2026-11-21 UTC).
  - `https://luc.devroye.org/Addario+Briend+Devroye+Donderwinkel+Kerriou+Lugosi-Random_friendship_trees-2026-.pdf` (Addario-Berry et al., *Random Structures & Algorithms* 69:e70082, 2026; read 2026-11-21 UTC).
  - `https://arxiv.org/abs/1303.1052` (read 2026-11-21 UTC) and `https://arxiv.org/pdf/1303.1052` (Cannings & Jordan, 2013; read 2026-11-22 UTC).
  - `https://arxiv.org/abs/2305.10628` (read 2026-11-21 UTC) and `https://arxiv.org/pdf/2305.10628` (Krapivsky & Redner, 2024; read 2026-11-22 UTC).
  - `https://arxiv.org/abs/2604.01540` (read 2026-11-21 UTC) and `https://arxiv.org/pdf/2604.01540` (Hartle et al., 2026; read 2026-11-22 UTC).
- **Not read (with reasons):**
  - `https://x.com/stevenstrogatz/status/1572287647715934210`: not read — X blocks unauthenticated automated access (and local files outside this worktree were not opened to preserve cold-audit isolation).
  - `https://doi.org/10.1002/rsa.70082` (Wiley landing page, e.g. online date `2026-07-30`): not read — verified the published article via Luc Devroye's open-access copy of the Wiley PDF (`Random Structures & Algorithms, 2026; 69:e70082`, accepted 4 June 2026).
  - J. Stat. Mech. (2017) 073405 publisher page: not read — paywalled; verified via arXiv:1706.01514 v2 PDF and TeX source.
  - Background citations in §2 table (Saramäki & Kaski 2004; Levens et al. 2022; Zhao et al. 2026; Hartle & Krapivsky 2025/2026; Angel et al. 2025; Krapivsky 2026; Briend and Hartle theses): not read — secondary references outside the primary bounds audit.

---

## 2. Hand Audit of the Proof (`REPORT.md` §§3.1–3.5)

No mathematical error or unbridgeable gap was found. Below is the step-by-step verification and a list of minor expositional notes (by section and line) that do not affect validity.

### 2.1 Lemma 1: One-step dynamics (§3.2, lines 104–125)
- For $`n \ge 3`$, $`T_n`$ is a connected tree with $`\ge 3`$ vertices, so no two leaves can be adjacent; every leaf has a unique neighbour in the core $`K = \{u : D_u \ge 2\}`$.
- For any $`w \in V(T_n)`$, $`P(W_n = w \mid T_n) = \frac{1}{n}\sum_{v\sim w} \frac{1}{D_v}`$.
  - If $`w = u \in K`$, its $`L_u`$ leaf neighbours contribute $`1`$ each and its $`c_u`$ core neighbours $`v`$ contribute $`1/D_v`$, giving $`(L_u + \rho_u)/n`$; attaching leaf $`n+1`$ increments $`L_u \to L_u + 1`$ ($`c_u`$ unchanged, $`C_n`$ unchanged).
  - If $`w = \ell`$ is a leaf with neighbour $`u \in K`$, $`P(W_n = \ell \mid T_n) = 1/(n D_u)`$; summing over the $`L_u`$ disjoint leaf neighbours of $`u`$ gives $`L_u/(n D_u)`$. Converting $`\ell`$ makes $`\ell \in K`$ with $`(L_\ell, c_\ell) = (1, 1)`$ and shifts $`(L_u, c_u) \to (L_u - 1, c_u + 1)`$ with $`D_u`$ unchanged, incrementing $`C_n`$ by $`1`$.
- The probabilities sum to $`\frac{1}{n}(\sum_{u\in K} L_u + \sum_{v\in K} c_v/D_v + \sum_{u\in K} L_u/D_u) = \frac{1}{n}((n-|K|) + |K|) = 1`$.
- Hence $`E[C_{n+1} - C_n \mid T_n] = S_n / n`$. **Valid.**

### 2.2 Lemma 2: Local states, event list, and drift identity (§3.3, lines 127–225)
- **Classification & permanence of "other":** For $`K_0 \in \{1, 2\}`$, a degree-2 core neighbour of $`u`$ has $`(L_w, c_w) \in \{(1,1), (0,2)\}`$ (class-1 pendant or class B), so every "other" core neighbour has $`D_w \ge 3`$. Since $`c_w`$ and $`D_w`$ never decrease and $`L_w`$ decreases only upon conversion (which raises $`c_w`$ to $`\ge 2`$ while keeping $`D_w \ge 3`$), an "other" neighbour remains "other" forever. **Valid.**
- **Remembered creator & labels ($`P`$, $`Q`$):** Every core vertex $`u`$ except the centre of $`T_3`$ is created when its creator $`p`$ converts $`u`$ from a leaf. Immediately after conversion, $`p`$ has $`(L_p - 1, c_p + 1)`$ with $`D_p`$ unchanged:
  - $`p`$ is a class-$`(L_p-1)`$ pendant of $`u`$ iff $`c_p = 0`$ and $`1 \le L_p - 1 \le K_0`$ (note $`c_p = 0, D_p \ge 2 \implies L_p \ge 2`$);
  - $`p`$ is class B iff $`(L_p, c_p) = (1, 1)`$;
  - otherwise $`p`$ is "other" (with $`D_p \ge 3`$), so $`p`$ is *remembered* and stays an "other" neighbour of $`u`$ forever, guaranteeing $`o_u \ge 1`$ forever.
  - Because $`D_p`$ never decreases, $`1/D_p \le 1/P`$ whenever $`P = \min(D_p(\text{birth}), P_{\max}) \ge 4`$, so $`\rho_u \in [B_0, B_0 + (o-1)/3 + 1/P]`$ (and $`[B_0, B_0 + o/3]`$ for $`P = 3`$).
  - For dynamic $`Q`$, while $`Q < Q_{\max}`$ we have $`D_p = Q`$ exact and $`p`$ gains a leaf at rate $`L_p + \rho_p \le L_p + c_p/2 \le D_p = Q`$ (in fact $`\le Q - 1/2`$ since $`c_p \ge 1`$, so $`[0, Q]`$ is a valid enclosing interval).
- **Tracked degrees variant (final lower bound, lines 193–209):** Classifying "other" neighbours by current degree $`d \in \{3, \dots, D_{\text{cap}}-1\}`$ (counts $`n_d`$) and $`\ge D_{\text{cap}}`$ (count $`\text{top}`$) depends only on $`T_n`$. Credit satisfies $`\rho_u \in [B_0 + \sum_d n_d/d,\; B_0 + \sum_d n_d/d + \text{top}/D_{\text{cap}}]`$. Newcomers leaving a pendant or B class enter at exact degrees:
  - E3$`_{K_0}`$: $`(K_0, 1) \to (K_0+1, 1)`$, degree $`K_0+2`$ ($`= 4`$ for $`K_0=2`$);
  - E4$`_\ell`$ ($`\ell \ge 2`$): $`(\ell, 1) \to (\ell-1, 2)`$, degree $`\ell+1`$ ($`= 3`$ for $`\ell=2`$);
  - E5: $`(0, 2) \to (1, 2)`$, degree $`3`$.
  - E7$`_d`$: each of the $`n_d`$ distinct degree-$`d`$ "other" neighbours gains a leaf with rate $`L_w + \rho_w \in [0, d]`$ (since $`\rho_w \le c_w/2 \le c_w`$), moving $`n_d \to n_d - 1`$ and $`n_{d+1} \to n_{d+1}+1`$ (or `top` if $`d = D_{\text{cap}}-1`$). Summing over disjoint choices of $`W_n`$ gives total rate in $`[0, n_d d]`$.
- **Exhaustiveness & disjointness of events E1–E7:** Partitioning $`W_n \in V(T_n)`$ by distance from $`u`$: distance $`0`$ is E1; distance $`1`$ splits disjointly into leaves (E2), class-$`\ell`$ pendants (E3$`_\ell`$), class B (E5), and "other" (E6 or E7$`_d`$); distance $`2`$ affects $`x_u`$ only when $`W_n`$ is a leaf neighbour of a pendant of $`u`$ (E4$`_\ell`$; since $`T_n`$ is a tree, leaf sets of distinct neighbours are disjoint, class-B vertices have $`L=0`$, and leaf conversions at "other" neighbours leave both "other" status and degree unchanged). Distance $`\ge 3`$ and non-leaf distance-$`2`$ vertices do not alter $`x_u`$.
- **Drift identity (★★):** Holds by linearity of conditional expectation across $`u \in K`$ (simultaneous transitions at adjacent core vertices do not affect expectation) plus the disjoint newborn sum $`\sum_{u\in K} (L_u/(n D_u)) f(\text{nb}(x_u))`$.
- **Minor expositional note (§3.3, lines 164–165):** The credit formulas $`(o-1)/3 + 1/P`$ ($`P \ge 4`$) and $`1/Q + (o-1)/3`$ ($`Q < Q_{\max}`$) implicitly use $`o \ge 1`$. This invariant follows from lines 138–144 (remembered creators stay "other" forever) and is closed under all transitions, so unreachable grid rows with $`o = 0, P \ge 4`$ (or $`o = 0, Q < Q_{\max}`$) are outside the reachable state space. Stating "$`(o \ge 1)`$" explicitly on lines 164–165 would avoid ambiguity for readers implementing a verifier from the spec.

### 2.3 Lemma 3 & Theorem/Corollary (§3.4 lines 226–245; §3.1 lines 77–102)
- Each uncertain rate ($`\rho_u`$, $`R_B`$, and E6/E7$`_d`$ rates) enters a separate linear term of $`\mathcal{A} f(x_u)`$, so taking the min (lower) or max (upper) independently over each interval's endpoints validly lower- or upper-bounds $`\mathcal{A}_{\text{true}} f(x_u)`$.
- Because $`T_n`$ has at most $`n`$ core vertices and $`f`$ is bounded on all states, $`\Phi_n = \sum_{u\in K} f(x_u)`$ is finite and integrable for every $`n \ge 3`$.
- $`T_3`$ is deterministically the path $`P_3`$ with unique core vertex $`(L=2, c=0, m=0, b=0, \text{no claim})`$, so $`\Phi_3 = f(x_{\text{centre}})`$.
- Lower: $`f \ge 0`$ everywhere and $`\sup f < \infty`$ imply $`E[C_n] \ge E[\Phi_n]/\sup f \ge (\Phi_3/\sup f)\prod_{m=3}^{n-1}(1+\delta/m)`$.
- Upper: $`f \ge 1`$ everywhere implies $`E[C_n] \le E[\Phi_n] \le \Phi_3 \prod_{m=3}^{n-1}(1+\lambda/m)`$.
- Product identity $`\prod_{m=3}^{n-1}(1+a/m) = 2\Gamma(n+a)/(\Gamma(n)\Gamma(3+a))`$ is exact for all $`n \ge 3`$ (empty product $`= 1`$ at $`n=3`$).
- **Almost-sure upper bound (Corollary, lines 93–98):** Since $`C_n`$ is nondecreasing, Markov's inequality on $`n = 2^j`$ with $`E[C_{2^j}] = O(2^{0.586 j})`$ gives $`\sum_j P(C_{2^j} > 2^{j(0.586+\varepsilon)}) < \infty`$; first Borel–Cantelli + monotonicity $`C_n \le C_{2^j} \le 2^{0.586+\varepsilon} n^{0.586+\varepsilon}`$ for $`2^{j-1} < n \le 2^j`$ gives $`\limsup_n \log C_n/\log n \le 0.586+\varepsilon`$ a.s.
- **Minor expositional notes (§3.1):**
  - Line 8 rounds the lower prefactor conservatively to `0.50`, whereas Theorem §3.1 (line 81) gives `0.502` (exact verified value `0.502326...`, and even the chained lower bound $`2 \cdot 0.875 / \Gamma(3.543) = 0.502029... > 0.502`$ holds; for upper, $`2 \cdot 3.955/\Gamma(3.586) = 2.162084... < 2.17`$).
  - Line 97–98 leaves the standard countable intersection over $`\varepsilon_k = 1/k \to 0`$ implicit when concluding $`\limsup_n \log C_n/\log n \le 0.586`$ a.s.

### 2.4 Infinite states outside the finite table (§3.5, lines 246–262) — Most critical check
The global extension rules (lines 248–250) define $`f(L, c, \cdot) = f(L_{\max}, c, \cdot)`$ for $`L > L_{\max}`$ and $`f \equiv 0`$ (lower) or $`f \equiv G`$ (upper) for $`c > c_{\max}`$.
- **Consistency on the corner quadrant $`L > L_{\max}, c > c_{\max}`$:** The two extension rules commute: clamping $`L \mapsto L_{\max}`$ first and then applying $`c > c_{\max}`$ gives the same constant ($`0`$ or $`G`$) as applying $`c > c_{\max}`$ directly. (Minor expositional note on lines 249–250: specifying "$`L > L_{\max}, c \le c_{\max}`$" on line 249 makes the partition visibly disjoint.)
- **Region (b): $`L > L_{\max}`$ and $`0 \le c \le c_{\max}`$ (lines 254–258):**
  - **E1:** $`L + 1 > L_{\max}`$, so $`f(L+1, \cdot) - f(L, \cdot) = f(L_{\max}, \cdot) - f(L_{\max}, \cdot) = 0`$ identically.
  - **E2 target:** Since $`L \ge L_{\max} + 1`$, we have $`L - 1 \ge L_{\max}`$, so $`f(L-1, c+1, \dots) = f(L_{\max}, c+1, \dots)`$ (or $`0`$/$`G`$ if $`c = c_{\max}`$), completely independent of $`L`$.
  - **Newborn invariance (line 255):** For $`L \ge L_{\max} + 1`$, $`L - 1 \ge L_{\max} > K_0`$ (not a pendant), $`(L, c) \ne (1, 1)`$ (not class B), and $`D = L + c \ge L_{\max} + 1`$. Thus $`\text{nb}(x_u)`$ is independent of $`L`$ **provided** $`L_{\max} + 1 \ge \max(K_0 + 2, P_{\max})`$ (upper) or $`L_{\max} + 1 \ge \max(K_0 + 2, D_{\text{cap}})`$ (lower). For all certificates in the report this holds with wide margin ($`61 \ge 20`$ for upper; $`61 \ge 6`$ and $`51 \ge 6`$ for lower9; $`61 \ge 8`$ for lower7). *Minor expositional note (§3.5, line 255):* this side condition on $`L_{\max}`$ is unstated in line 255.
  - **Dependence on $`L/D`$ and $`1/D`$ (lines 257–258):** After taking the worst case over the uncertain rates (whose worst-case endpoints depend only on the signs of $`\Delta f`$, which are independent of $`L`$ for $`L > L_{\max}`$):
    - E4, E6, E7 are independent of $`L`$ and $`D`$;
    - E2 has coefficient $`C_{\text{E2}} \cdot (L/D)`$;
    - E3$`_\ell`$ has $`m_\ell \ell \Delta_{\text{E3}_\ell} + m_\ell \Delta_{\text{E3}_\ell} \cdot (1/D)`$;
    - E5 has $`\text{ext}(0, (b/2)\Delta_{\text{E5}}) + b \Delta_{\text{E5}} \cdot (1/D)`$.
    - Setting $`y = 1/D \in (0, 1/(L_{\max}+1+c)]`$, we have $`L/D = 1 - c y`$, so the worst-case drift is an **exact affine function of $`y`$** whose extrema over all $`L > L_{\max}`$ are attained at $`y = 0`$ and $`y = 1/(L_{\max}+1+c)`$ (our check `T1`). Moreover, relaxing $`(L/D, 1/D)`$ independently to the rectangle $`[(L_{\max}+1)/(L_{\max}+1+c), 1] \times [0, 1/(L_{\max}+1+c)]`$ (lines 257–258, our checks `T2` and `T3`) is strictly more conservative and therefore sufficient. Our independent verifier confirmed that **all three** tail bounds (`T1`, `T2`, `T3`) pass with zero failures on every table.
- **Region (c): $`c > c_{\max}`$ (lines 259–261):**
  - Because $`c`$ never decreases in any event E1–E7, every post-event state of $`u`$ still has $`c > c_{\max}`$, so $`f(\text{new}) - f(x_u) = 0`$ identically for E1–E7. Only the newborn creation from E2 remains: $`\mathcal{A} f(x_u) = (L/D) f(\text{nb}(x_u))`$ exactly.
  - **Lower:** $`f(x_u) = 0`$ and $`f(\text{nb}) \ge 0`$ everywhere $`\implies \mathcal{A} f(x_u) \ge 0 = \delta \cdot 0`$.
  - **Upper:** $`f(x_u) = G \ge 1`$ and $`0 \le L/D < 1`$, so $`\mathcal{A} f(x_u) \le f(\text{nb}(x_u))`$. Note that for $`c > c_{\max} = 10`$ and $`L \ge 1`$, $`D = L + c \ge 12`$, so only newborns with $`P \in \{12, \dots, 20\}`$ are actually reachable; our verifier checked that **both** $`\max_{P \ge 12} f(1,1,0,0,P) \approx 2.661580 \le \lambda G`$ **and** the stronger $`\max_{\text{all nb}} f(\text{nb}) \approx 2.919026 \le \lambda G \approx 2.921471`$ hold in exact arithmetic!
- **Global sup/inf:** Since the $`L > L_{\max}`$ tail copies values from $`L = L_{\max}`$ and $`c > c_{\max}`$ takes constant $`0`$ (lower) or $`G`$ (upper), checking $`f \ge 0`$ on the box guarantees $`f \ge 0`$ and $`\sup f = \max_{\text{box}} f`$ globally (lower), and checking $`f \ge 1`$ on the box plus $`G \ge 1`$ guarantees $`f \ge 1`$ globally (upper). **Completely airtight.**

---

## 3. Independent Exact-Rational Verifier & Exact Tree Enumeration

Written strictly from §§3.3–3.5 in pure Python stdlib (no third-party packages, no access to `certs*.py` or `independent/`):
- `checks/strogatz-gemini/verify_certs.py`
- `checks/strogatz-gemini/verify_exact_trees.py`

Tables were dumped to temporary storage via `attempts/x-strogatz/dump_tables.py` and read solely as text (`float.hex` parsed directly to exact dyadic integers and cross-checked against `Fraction.from_float(float.fromhex(...))`). Every state check uses 100% exact arbitrary-precision integer arithmetic (clearing exact common denominators) with periodic direct `fractions.Fraction` cross-checks (asserting exact equality on thousands of states) and closure assertions confirming no transition ever leaves the valid state set.

### 3.1 Exact-rational certificate results

| Table | Verified exponent | Valid box states | Tail states | Min slack (box / tail T1 / tail T3) | Constants | Result |
|---|---|---:|---:|---|---|---|
| `upper_table.txt` (60×10, $`K_0=1, P_{\max}=20`$) | $`\lambda = 293/500 = 0.586`$ | 245,564 | 4,026 | $`2.445\times 10^{-3}`$ / $`2.445\times 10^{-3}`$ / $`2.445\times 10^{-3}`$ | $`\min f = 1`$, $`G = 4.985445`$, $`\Phi_3 = 3.954624 \le 3.955`$, $`2\Phi_3/\Gamma(3.586) = 2.161878 \le 2.17`$; $`\max_{\text{all}} f(\text{nb}) = 2.919026 \le \lambda G = 2.921471`$ | **PASS** (0 fails) |
| `lower9b_table.txt` (60×10, $`K_0=2, D_{\text{cap}}=6`$) | $`\delta = 543/1000 = 0.543`$ | 1,186,319 | 19,448 | $`1.083\times 10^{-6}`$ (on $`f>0`$; 9 zero-$`f`$ states have slack $`0`$) / $`2.051\times 10^{-4}`$ / $`2.051\times 10^{-4}`$ | $`\min f = 0`$, $`\sup f = 0.639069`$, $`\Phi_3 = 0.559516`$, $`\Phi_3/\sup f = 0.875518 \ge 0.875`$, $`2(\Phi_3/\sup f)/\Gamma(3.543) = 0.502326 \ge 0.502`$ | **PASS** (0 fails) |
| `lower9_table.txt` (50×9, $`K_0=2, D_{\text{cap}}=6`$) | $`\delta = 27/50 = 0.540`$ | 583,431 | 11,440 | $`3.015\times 10^{-6}`$ ($`f>0`$; 8 zero-$`f`$ states) / $`6.430\times 10^{-4}`$ / $`6.430\times 10^{-4}`$ | $`\Phi_3/\sup f = 0.879102`$, Gamma prefactor $`0.506075`$ | **PASS** (0 fails) |
| `lower_table.txt` (60×10, $`K_0=2, Q_{\max}=8`$) | $`\delta = 107/200 = 0.535`$ | 279,125 | 4,576 | $`3.472\times 10^{-7}`$ ($`f>0`$; 9 zero-$`f`$ states) / 0 fails / 0 fails | $`\Phi_3 = 0.502007`$, $`\sup f = 0.577726`$ | **PASS** (0 fails) |

Note on `upper_table.txt` state count: the dumped file contains 314,028 rows (full grid), of which 68,464 are unreachable ($`L + c < 2`$, or $`o = 0`$ with $`P \ge 4`$). Our verifier explicitly verified completeness and closure of the 245,564 valid reachable-superset box states (plus 4,026 tail states).

### 3.2 Negative controls (all failed as required)
1. **Perturbed `lower9b_table.txt` (+0.01 at initial state `(2,0,0,0,0,0,0,0)`):** **FAIL** — caught at `(2,0,0,0,0,0,0,0)` with slack $`-0.035249`$.
2. **Perturbed `lower9b_table.txt` ($`+10^{-4}`$ at newborn state `(1,1,0,0,0,0,0,0)`):** **FAIL** — caught at `(1,1,0,0,0,0,0,0)` with slack $`-9.005\times 10^{-5}`$.
3. **Perturbed `upper_table.txt` ($`-0.05`$ at `(2,0,0,0,3)`):** **FAIL** — caught at `(2,0,0,0,3)` with slack $`-0.176855`$.
4. **Unperturbed `lower9b_table.txt` at stricter $`\delta = 544/1000`$:** **FAIL** — 1,186,311 box/tail failures.
5. **Unperturbed `upper_table.txt` at stricter $`\lambda = 58/100`$:** **FAIL** — 241,540 box/tail failures.
6. **Perturbed `upper_table.txt` with $`G`$ decreased by $`1.0`$:** **FAIL** — caught by `C_GT_CMAX_FAIL` ($`\max f(\text{nb}) > \lambda G`$).

### 3.3 Exact enumeration over all small trees & extreme-tree brute force (`verify_exact_trees.py`)
- **All unlabeled trees up to $`n = 11`$** (OEIS A000055 counts `1, 1, 2, 3, 6, 11, 23, 47, 106, 235` matched exactly; exact `Fraction` probabilities):
  - Verified Lemma 1 ($`E[C_{n+1} - C_n \mid T_n] = S_n / n`$) and Lemma 2 tracked-degrees drift identity (★★) with true rates **identically** on all 199 unlabeled trees ($`n = 3..10`$) using both an injective rational hash function and the actual final lower certificate `lower9b_table.txt`, plus rate-interval containment and $`n E[\Phi_{n+1}-\Phi_n \mid T_n] \ge 0.543 \Phi_n`$.
  - Confirmed exact $`E[C_{10}] = 200684675729 / 80621568000 \approx 2.48921822`$ (matching `REPORT.md` line 319) and $`E[C_{11}] = 446709364753753 / 169305292800000 \approx 2.63848435`$.
- **All static-label-$`P`$ history-equivalence classes up to $`n = 9`$** (`1, 1, 2, 4, 12, 36, 121, 416` distinct `(tree, P, creator)` classes):
  - Verified the remembered-creator invariant ($`p`$ stays "other", $`o \ge 1`$, $`D_p \ge P`$), rate-interval containment, Lemma 2 static-$`P`$ drift identity (★★) (injective hash + actual `upper_table.txt`), and $`n E[\Phi_{n+1}-\Phi_n \mid \text{state}] \le 0.586 \Phi_n`$ identically on all 176 history-states ($`n = 3..8`$).
- **Explicit extreme trees & random RFTs ($`n \le 94`$):**
  - Constructed extreme trees hitting `L > 60` (`maxL = 64`), `c > 10` (`maxc = 13`), both simultaneously (`L = 69, c = 12`), `P = 20`, `top >= 1`, and all neighbour classes `m1, m2, b, n3, n4, n5`, plus 5 random friend trees at $`n = 75`$. Exact brute-force enumeration over all $`(V, W)`$ moves matched (★★) identically and satisfied both final certificate inequalities on every tree.
  - Drift negative control: corrupting the E4$`_1`$ rate from $`m_1/2`$ to $`m_1/3`$ was immediately caught by the exact drift check.

---

## 4. Small-Size C Simulation Check (`verify_sim.py`)

Inspected `attempts/x-strogatz/sim.c` (confirmed exact implementation of the model with xoshiro256** and unbiased Lemire range sampling), compiled it unmodified into temporary storage, and ran small sizes only:
- **Unmodified `sim.c` at `NMAX=1000`, 20,000 runs (seed 42):**
  - $`n = 10`$: $`\hat{E}[C_{10}] = 2.49790 \pm 0.00815`$ (exact $`2.48922`$, $`z = +1.06`$), $`\hat{E}[S_{10}]/10 = 0.14930 \pm 0.00023`$ (exact $`0.14927`$, $`z = +0.16`$).
  - $`n = 1000`$: $`\hat{E}[C_{1000}] = 35.651 \pm 0.127`$, $`\hat{\mu}_{\text{eff}} = 0.56906 \pm 0.00032`$.
  - Checkpoint increments $`\hat{E}[C_{n_{k+1}}] - \hat{E}[C_{n_k}]`$ agree with the interpolated sum $`\sum_{m=n_k}^{n_{k+1}-1} \hat{E}[S_m]/m`$ across all intervals ($`|z| \le 1.14`$).
- **Unmodified `sim.c` at `NMAX=10000`, 2,000 runs (seed 20261121):**
  - $`n = 10000`$: $`\hat{E}[C_{10^4}] = 134.80 \pm 1.51`$, $`\hat{\mu}_{\text{eff}} = 0.56627 \pm 0.00054`$ (consistent with `REPORT.md` line 333: $`0.56648 \pm 0.00062`$); all checkpoint intervals within $`|z| \le 1.13`$.
- **Direct consecutive-step C check (50,000 runs, paired $`(C_{n+1}-C_n) - S_n/n`$):**
  - $`n = 10`$: $`\hat{E}[C_{11}-C_{10}] = 0.146860`$ vs $`\hat{E}[S_{10}]/10 = 0.149315`$ (paired diff $`-0.002455 \pm 0.001577`$, $`z = -1.56`$).
  - $`n = 20`$: $`0.111480`$ vs $`0.109978`$ ($`z = +1.07`$).
  - $`n = 50`$: $`0.075740`$ vs $`0.074271`$ ($`z = +1.25`$).
  - $`n = 100`$: $`0.055680`$ vs $`0.055143`$ ($`z = +0.53`$).
  - $`n = 500`$: $`0.026960`$ vs $`0.027461`$ ($`z = -0.69`$).
  - $`n = 1000`$: $`0.020900`$ vs $`0.020328`$ ($`z = +0.90`$).

All numerical results confirm $`E[C_{n+1}] - E[C_n] = E[S_n]/n`$.

---

## 5. Summary of Remaining Unchecked Items

None within the scope of the theorem. (As instructed, the long $`n = 10^7, 10^8`$ Monte Carlo simulation runs in §4 were not rerun, and secondary background references outside the primary bounds were not read.)
