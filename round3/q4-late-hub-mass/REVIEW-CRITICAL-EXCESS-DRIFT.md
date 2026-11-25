> **Public checkpoint export note:** Edited export of `late-mass-neon/REVIEW-CRITICAL-EXCESS-DRIFT-16f0c12.md` (frozen original SHA-256 `935d2c0a592435867778f6a83226617ee74e7d23997e74738f8af077a0e94bc7`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

VERDICT: HOLDS (scoped to `CRITICAL-EXCESS-DRIFT.md` and `CRITICAL-EXCESS-EVIDENCE.md` at `16f0c1202cee968fc361fd34519acad03e461c72`, together with the four cited upstream items (2), (3), (6), (8) of `ADDENDUM-CUTOFF-THRESHOLD.md` which were independently audited from scratch here; conditional on the listed inherited Round 1/2 baseline and `9be4767` reviewed inputs)

# Independent Cross-Family Audit of Astra's Critical-Scale Excess-Degree Drift (`16f0c12`)

**Reviewer:** Neon lead (`late-mass-neon`, Gemini 4 Neon)  
**Target commit:** `16f0c1202cee968fc361fd34519acad03e461c72` (`research/ft3-late-mass-astra`, read-only board export)  
**Pinned target files & verified SHA-256 digests:**
- [`CRITICAL-EXCESS-DRIFT.md`](CRITICAL-EXCESS-DRIFT.md) — `4c27a043e261f9c2f3e9db1dd529dfce70cce0c8f548e875039fb85e262fa4e8`
- `CRITICAL-EXCESS-EVIDENCE.md` (not exported) — `395997919be9a72d7b7e769a076430f44fb19b93cb51490c36b6a1b2fc21af02`
- Cited upstream post-review addendum audited for its four imported ingredients ((2), (3), (6), (8)): `ADDENDUM-CUTOFF-THRESHOLD.md` (not exported) — `ab84cf0e09826c43a56e020387ce7256849acc7c1ab36976f135d97f17405bb3`

**Method:** Line-by-line adversarial proof audit across four separate attack surfaces, direct audit of the four imported claims from the previously unreviewed `ADDENDUM-CUTOFF-THRESHOLD.md`, explicit center-specificity analysis across $`q_n`$, $`\pi_n`$, and $`p_n`$, derivation of independent finite-tree and limiting identities, and a zero-shared-code exact rational verification (distinct from replaying Astra's `check_critical_drift.py`). Cross-family AI review is not human referee verification.

---

## 1. Executive Summary

| Component | Statement Audited | Verdict | Notes |
|---|---|---|---|
| Pillar 1 (`CRITICAL-EXCESS-DRIFT.md` §§1–4) | Crossing count $`t/(C_m L) \to h = M_0/\mu`$ a.s. (2), cohort lower bound (4), and uniform reciprocal-degree average $`\sup_{n \le j \le n^2} \lvert \frac{1}{t}\sum_{y\text{ crossing at }n} 1/d_y(j) - b\rvert  \to 0`$ a.s. (5) | **HOLDS** | Simultaneous registration coupling, pre-registration bound $`d K_m`$, Yule domination on $`G`$, cut-tree interior disjointness, incoming-arrow bound $`O(C_m)`$, finite-age truncation $`S_0`$, and all-prefix Hoeffding over $`k \ge C_m`$ all verify cleanly. Two minor micro-step clarifications noted below (non-rectangular marked region and bare-tree gate measurability); neither is a gap. |
| Pillar 2 (`CRITICAL-EXCESS-DRIFT.md` §5 & `ADDENDUM-CUTOFF-THRESHOLD.md` (2), (3), (8)) | Exact decomposition $`R_m - q_n = M_n + V_n - 2 q_n / n`$ and finite bound $`\mathbb{E}[\sum_{j \ge n} (t_j - t + k_j)/[j(j+1)] \mid \mathcal{F}_n] \le 24 k_n / n`$ via $`\mathbb{E}[k_j \mid \mathcal{F}_n] \le 12 k_n`$ | **HOLDS** | Re-proved from scratch. Exact graph identity holds for every edge-seeded tree $`T_j`$ ($`j \ge 3, 2 \le m \le j`$); per-host old-leaf expected promotion cap verified to be $`\le 10.1771 < 11`$; exact Tonelli double sum $`= 1/n + 1/n`$ gives $`24 k_n / n`$. |
| Pillar 3 (`CRITICAL-EXCESS-DRIFT.md` §5 & `ADDENDUM-CUTOFF-THRESHOLD.md` (6)) | Conditional $`L^1`$ concentration $`\mathbb{E}[\lvert n V_n/t - (1-b)\rvert  \mid \mathcal{F}_n] \to 0`$ a.s. (3), martingale CLT $`M_n / \sqrt{q_n/n} \mid \mathcal{F}_n \Rightarrow N(0,1)`$ a.s. at critical scale, and Kolmogorov limit (1) with shift $`\delta(X) = (1-b)\frac{M_0}{1+\mu}\sqrt{\lambda X/c}`$ | **HOLDS** | Audited `ADDENDUM-CUTOFF-THRESHOLD.md` eq. (6) directly: its conditional characteristic-function proof applies to every deterministic $`m \to \infty, m = o(n)`$ (including the critical scale) and does not depend on `ADDENDUM-GAUSSIAN.md` or `ADDENDUM-DIRICHLET.md`. Bounded full-path majorant $`\sup_{r \ge N} B_r`$, conditional Slutsky without independence, and conditional mean limit all verify. |
| Pillar 4 (`CRITICAL-EXCESS-EVIDENCE.md`) | Exact rational enclosure $`b \in [0.0652016, 0.1166501]`$, $`1-b \in [0.8833499, 0.9347984]`$, and convergence to $`b`$ as size cap $`\to \infty`$ | **HOLDS** | Upper payoff $`d^{-\ell-1}`$, protection coupling, generator inequality $`(6-2j)/[(j+1)(j-1)(j-2)] \le 0`$, absorbing recurrence, union bound, backwards intersection, and cap $`\to \infty`$ convergence proof all hold. Independently reconstructed with zero shared code (Dyck encoding, source-pair kernel, Thomas tridiagonal solve) reproducing `[0.0652016, 0.1166501]` identically. |
| Center-Specificity | Scaling $`m^{1+\mu}(\log n)^2/n \to \lambda`$ and shift $`+\delta(X)`$ apply specifically to the observed excess-degree center $`q_n`$ | **HOLDS** | Verified directly via our independent finite-tree identities: at the exact same scale, the star-weight center $`\pi_n`$ has negative shift $`-\frac{b}{1-b}\delta(X) = -b\frac{M_0}{1+\mu}\sqrt{\lambda X/c}`$, while the next-recipient center $`p_n`$ has **zero shift** ($`N(0,1)`$). |

---

## 2. Detailed Adversarial Audit by Pillar

### 2.1 Pillar 1: Crossing-Edge Count (2) and Uniform Reciprocal Endpoint Average (5)

1. **Critical scaling and Borel–Cantelli summability:**
   - From $`m^{1+\mu}(\log n)^2 / n \to \lambda \in (0, \infty)`$ with $`\mu \in (1/2, 1)`$: $`\log m / \log n \to 1/(1+\mu) \in (1/2, 2/3)`$, $`L = \log(n/m) \sim \frac{\mu}{1+\mu}\log n`$, and $`m \sim (\lambda n / (\log n)^2)^{1/(1+\mu)} \to \infty`$ polynomially in $`n`$ with $`m = o(n)`$.
   - On the $`\mathcal{F}_m`$-measurable size gate $`C_m \ge m^{\mu/2}`$ (which holds eventually a.s. since $`C_m / m^\mu \to X > 0`$), $`C_m / L \ge n^{\frac{\mu}{2(1+\mu)} - o(1)} \gg n^{1/7}`$, so $`\exp(-a_B C_m / L)`$ and $`\exp(-c C_m / d)`$ ($`d = \max(3, \lceil \sqrt{L}\rceil)`$) are summable over **all** integers $`n \ge 3`$ even when $`m(n)`$ repeats values across consecutive $`n`$.
   - The pathwise bound $`K_m = O_\omega(C_m)`$ on early labels $`\le m`$ that **ever** become cores follows genuinely from `SMALL-POWERS.md` lines 12–23: the private per-host bound $`B_x = P_x + L \mathbf{1}_{\mathrm{exc}}`$ bounds all promotions of old leaves of $`x`$ (not only ultimate hubs), so $`K_m \le C_m + \sum_{x \text{ core at } m} B_x`$; monotonicity of $`K_m`$ in $`m`$ and regular variation of $`C_m`$ interpolate dyadic $`m = 2^j`$ to all $`m`$.

2. **Upper bound $`\limsup t/(C_m L) \le h`$ (`CRITICAL-EXCESS-DRIFT.md` lines 35–41):**
   - *Exact crossing-creation characterization:* Every vertex $`y > m`$ enters as a leaf attached to a single parent $`x`$. If $`x \le m`$ was a leaf when $`y`$ attached, $`x`$ becomes an early core while $`y`$ is a young leaf. While $`y`$ is a leaf, $`x`$ is its unique neighbor; hence a crossing core edge $`(x, y)`$ ($`x \le m < y`$, both cores) is created if and only if early core $`x`$ rings and promotes young leaf $`y`$.
   - *Pre-registration cost $`\le d K_m`$:* Vertices $`x \le m`$ that never become cores create zero crossings. Any $`x \in K_m`$ with $`d_x(m) \ge d`$ has zero young neighbors until its first post-$`m`$ degree gain, at which instant it registers immediately (0 pre-registration crossings). Any $`x \in K_m`$ with $`d_x(m) < d`$ maintains degree $`< d`$ strictly prior to registration, so it has $`< d`$ neighbors total and creates at most $`d - 1 \le d`$ crossing edges before registration. Total pre-registration crossings $`\le d K_m = O_\omega(C_m \sqrt{L}) = o_\omega(C_m L)`$.
   - *Simultaneous registration coupling (line 37):* At registration, $`x`$ just received a fresh young leaf and satisfies $`D_x(s) \ge d \ge 3`$ for all $`s \ge 0`$. The virtual linear birth-death chain $`J_x`$ (started at $`1`$, birth rate $`j`$, death rate $`j/d`$, 0 absorbing) stays $`\le L_x^{\mathrm{phys}}`$: virtual births drive physical leaf-source births (supplemented by residual rate $`L_x^{\mathrm{phys}} - J_x \ge 0`$); at an $`x`$-ring with physical leaves ordered first in $`[0, L_x^{\mathrm{phys}}/D_x)`$, when $`J_x = L_x^{\mathrm{phys}}`$ we have $`L_x^{\mathrm{phys}}/D_x = J_x/D_x \le \min(1, J_x/d)`$, so every physical leaf loss triggers a virtual death; when $`J_x < L_x^{\mathrm{phys}}`$, losing one leaf cannot invert the inequality; extra deaths at rate $`(J_x/d - 1)_+`$ bring total virtual death rate to $`\min(1, j/d) + (j/d - 1)_+ = j/d`$. Distinct hosts have disjoint leaf sets; any promoted old leaf takes a fresh core clock at promotion by memorylessness; independent registration tapes can be pre-sampled before $`T_m`$ and assigned in registration order at stopping times without altering any physical rate.
   - *No look-ahead / filtration bug (lines 39–41):* Survival is **not** used in the stopping rule. The compensated counting bound for post-registration source rings stops at the **adapted** stopping time $`\tau = \inf\{s \ge 0 : k(T_m + s) > B C_m\} \wedge (L + 1)`$ (where predictable intensity is $`\le B C_m`$ and total compensator $`\le B C_m (L + 1)`$), yielding $`|N_{\mathrm{post}}(\tau) - \Lambda_{\mathrm{post}}(\tau)| \le C_m`$ with failure $`\le 2\exp(-a_B C_m/(L+1))`$. Separately, padding the independent registration tapes to $`N_B = \lceil B C_m \rceil`$ ($`\mathcal{F}_{T_m}`$-measurable), the extinction count is $`\mathrm{Binomial}(N_B, 1/d)`$ (by gambler's ruin $`q = 1/d`$), exceeding $`2 N_B / d`$ with conditional probability $`\le \exp(-c_0 B C_m / d)`$. On the intersection with $`\{K_m \le B C_m, T_n - T_m \le L + 1\}`$, we have $`\tau = L + 1 \ge T_n - T_m`$, at most $`N_B`$ tapes are used, every surviving registered vertex belongs to $`H_m`$ by the physical degree dichotomy, and therefore pathwise $`\Lambda_{\mathrm{post}}(L + 1) \le (H_m + 2 N_B / d)(L + 1)`$. Dividing by $`C_m L`$ and using $`H_m / C_m \to h`$, $`d \sim \sqrt{L} \to \infty`$ gives $`\limsup t / (C_m L) \le h`$ a.s. Because the resulting inequality concerns a purely physical event, each fixed-$`n`$ coupling transfers its summable failure probability directly to the physical space for Borel–Cantelli.

3. **Selected cohort lower bound (4) (`CRITICAL-EXCESS-DRIFT.md` lines 45–51):**
   - *Yule domination at $`j_0 = \lfloor m e^{\varepsilon L} \rfloor`$:* Restarting continuous time at $`T_{j_0}`$, $`e^{-t} N(T_{j_0} + t) \to S \sim \mathrm{Gamma}(j_0, 1)`$, and for each fixed $`x \le m`$ the degree gain rate satisfies $`r_x(t) = \sum_{w \sim x} 1/D_w(t) \le D_x(t)`$, so $`S Z_x`$ is marginally stochastically dominated by $`\mathrm{Gamma}(d_x(j_0), 1)`$. On $`\{S \ge j_0/2\}`$ and $`d_x(j_0) < r j_0 / (4m)`$, $`\mathbb{P}(\mathrm{Gamma}(d_x(j_0), 1) > r j_0 / (2m) \mid \mathcal{F}_{j_0}) \le (2/e)^{r j_0/(4m)}`$ by exponential Markov at $`\theta = 1/2`$, while $`\mathbb{P}(S < j_0/2) \le (\sqrt{e}/2)^{j_0}`$. Since $`j_0 / m \sim n^{\varepsilon \mu/(1+\mu)}`$ grows polynomially in $`n`$, union-bounding over $`x \le m`$ is summable in $`n`$.
   - *Micro-step verification on $`h(r)`$ (line 27 & line 45):* Claim 22 states convergence on bounded marked rectangles $`[a, b] \times (u, \infty)`$ ($`a, u > 0`$), whereas $`\{(v, u) \in (0, 1] \times (0, \infty) : u > r v\}`$ reaches the origin. For the lower bound $`\liminf |G|/C_m \ge h(r) - \pi_{\mathrm{core}}(d > D)`$ used in line 45, approximating from inside by finite unions of rectangles $`(t_{i-1}, t_i] \times (r t_i, \infty)`$ with $`t_0 = \delta > 0`$ immediately yields $`\liminf \#\{x \le m : Z_x > r/m\}/C_m \ge \int_\delta^1 v^{\mu-1}\nu((r v, \infty))\,dv \to h(r)`$ as mesh $`\to 0`$ and $`\delta \downarrow 0`$ (and the matching upper bound follows from $`H_{\lfloor \delta m\rfloor} \sim \delta^\mu h C_m`$). Moreover, early leaves at $`m`$ have $`d_x(m) = 1 \le D`$ and are automatically included in $`\{x \le m : d_x(m) \le D\}`$ while $`\{x \le m : d_x(m) > D\}`$ consists strictly of cores ($`D \ge 2`$), so Claim 10's local LLN gives $`\#\{x \le m : d_x(m) > D\}/C_m \to \pi_{\mathrm{core}}(d > D)`$ a.s.
   - *Young-leaf fraction on $`G`$ and eq. (4):* For every $`x`$, source rings through size $`n^2`$ form a sum of independent $`\mathrm{Bernoulli}(1/j)`$ indicators with mean $`\le 1 + 2\log n`$; exponential Markov at $`\theta = 1`$ gives $`\mathbb{P}(R_x(n^2) > 8 \log n) \le e^{e-1} n^{2(e-1) - 8} \le e^2 n^{-4.5}`$, summable over $`n^2`$ labels and all $`n`$. Every non-young-leaf neighbor of $`x \in G`$ at size $`j \in \{j_0, \dots, n\}`$ is either one of $`\le D`$ old neighbors or a young neighbor promoted at a ring of $`x`$ ($`\le 8 \log n`$), while $`d_x(j) \ge r j_0 / (4m)`$ grows polynomially. Thus the young-leaf fraction is $`1 - o(1)`$ uniformly on $`G`$ and $`j \in \{j_0, \dots, n\}`$. Since $`G`$ is $`\mathcal{F}_{j_0}`$-measurable, discrete Freedman on $`\{|G| \le B C_m, C_m \ge m^{\mu/2}\}`$ with $`\sum_{j=j_0}^{n-1} 1/j = (1 - \varepsilon) L + O(1)`$ yields (4). Sending $`D \to \infty`$, $`r \downarrow 0`$, $`\varepsilon \downarrow 0`$ along countable sequences combines with the upper bound to prove both (2) and $`(t - t_G)/t \to 0`$ a.s.

4. **Uniform reciprocal endpoint average (5) (`CRITICAL-EXCESS-DRIFT.md` lines 55–69):**
   - *Cut-tree interior disjointness & inert coupling:* Every crossing $`(x, y)`$ has $`x \le m < y`$, and every vertex in the subtree below $`(x, y)`$ has label $`\ge y > m`$ and unique path to $`\{1, \dots, m\}`$ through $`(x, y)`$. Therefore distinct crossing subtrees are pairwise disjoint, and only $`x`$ can send an incoming arrow into $`y`$'s subtree. Up to the first incoming arrow $`x \to y`$, the subtree rooted at $`y`$ (with edge $`(y, x)`$ mapped to the exterior parent slot) couples identically to an independent inert virgin tape $`\mathrm{Tape}_i`$: internal arrows match identically, and a root send $`y \to x`$ attaches a leaf to $`x`$ outside $`y`$'s subtree while producing a null event in the inert copy.
   - *Incoming-arrow bound through $`n^2`$:* On $`\{|G| \le B C_m\}`$ and the observable gate $`\{\max_{x \in G} c_x^{\mathrm{young}}(j) \le M \log n\}`$ (micro-step clarification: gating directly on the $`\mathcal{F}_j`$-measurable crossing-child count $`c_x^{\mathrm{young}}(j) \le R_x(n^2)`$ requires no source-augmented history in $`\mathcal{F}_j`$), the conditional incoming-arrow probability at step $`j \in \{j_0, \dots, n^2 - 1\}`$ is $`\le B C_m \cdot (M \log n) \cdot \frac{4 m}{r j_0 j}`$. Summing over $`j \le n^2`$ gives compensator $`O_B(C_m (\log n)^2 m / (r j_0)) = o_B(C_m)`$; Freedman bounds actual incoming arrows by $`\le 2 C_m = o(t)`$ eventually a.s.
   - *Finite age $`S_0`$ and all-prefix Hoeffding (no independence of $`t`$ from tapes needed):* Crossings created in $`(\lfloor n e^{-2S_0} \rfloor, n]`$ have stopped compensator $`\le B C_m (2 S_0 + 1)`$, hence count $`O_{\omega, S_0}(C_m) = o_{\omega, S_0}(t)`$. Earlier crossings have age $`\ge S_0`$ at $`T_n`$ (hence at all $`T_j`$, $`n \le j \le n^2`$) since $`T_n - T_{\lfloor n e^{-2S_0}\rfloor} \to 2 S_0 > S_0`$ a.s. For any good (selected, age $`\ge S_0`$, unspoiled through $`n^2`$) crossing $`i`$, monotonicity of $`D_i^0`$ gives $`0 \le 1/d_{y_i}(j) - U_i \le V_i(S_0) = 1/D_i^0(S_0) - U_i`$ **uniformly** in $`n \le j \le n^2`$. For every bad crossing, $`|1/d_{y_i}(j) - U_i| \le 1/2`$. Summing $`V_i(S_0) \ge 0`$ over the entire prefix $`1 \le i \le t`$, we obtain

```math
\sup_{n \le j \le n^2} \left|\frac{1}{t}\sum_{i=1}^t \frac{1}{d_{y_i}(j)} - b\right| \le \left|\frac{1}{t}\sum_{i=1}^t U_i - b\right| + \frac{1}{t}\sum_{i=1}^t V_i(S_0) + \frac{|\mathcal{B}|}{2 t}.
```

Crucially, because the tapes are independent of $`\mathcal{F}_{T_m}`$ ($`C_m`$ known), Hoeffding for $`[0, 1/2]`$-valued i.i.d. variables union-bounded over **all** integer prefixes $`k \ge C_m`$ gives conditional failure $`\le \frac{2}{1 - e^{-8\eta^2}} e^{-8 C_m \eta^2}`$, summable in $`n`$ on $`\{C_m \ge m^{\mu/2}\}`$. Since $`t \ge C_m`$ eventually by (2), the prefix means apply at the **random, tape-dependent** index $`k = t`$ without requiring independence between $`t`$ and the tapes. Sending $`(D, r, \varepsilon)`$, then $`S_0 \to \infty`$ ($`\mathbb{E}[V_1(S_0)] \to 0`$ by bounded convergence, no rate needed), then $`\eta \downarrow 0`$ proves (5).

---

### 2.2 Pillar 2: Negligibility of New Crossings and Upstream `ADDENDUM-CUTOFF-THRESHOLD.md` Eqs. (2), (3), (8)

Because `ADDENDUM-CUTOFF-THRESHOLD.md` belongs to the previously unreviewed `b404108` continuation, we independently audited its imported equations (2), (3), and (8) line-by-line rather than treating them as black-box citations.

1. **Exact finite-tree identity (2) and drift decomposition (3):**
   - For any edge-seeded tree $`T_j`$ ($`j \ge 3`$) and any cutoff $`2 \le m \le j`$:
     - Early leaves ($`w \le m, d_w = 1`$) have their unique neighbor in $`\{1, \dots, m\}`$ (since $`m \ge 2`$ induces a connected tree), contributing $`0`$ to $`j p_j`$.
     - Early cores contribute $`a_j = \sum_{x \le m, d_x \ge 2} l_x / d_x \in [0, k_j]`$.
     - Late leaves contribute $`\#\{\text{late leaves adjacent to a late core}\}`$.
     - Each late core $`y > m`$ has at most one early neighbor (its parent, which must be an early core), precisely when $`y`$ is the late endpoint of one of the $`t_j`$ crossing core edges; hence late cores contribute $`k_L - \sum_{y \text{ crossing at } j} 1/d_y(j)`$, where $`k_L`$ is the number of late cores.
     - Double-counting half-edges of the late core forest (which has $`k_L`$ vertices and $`t_j`$ components, hence $`k_L - t_j`$ internal edges) yields $`A_j = \sum_{v > m} (d_v - 1) = k_L - t_j + \#\{\text{late leaves adjacent to a late core}\}`$.
     - Subtracting cancels both $`k_L`$ and late leaves identically, giving the exact finite-tree identity

```math
j p_j - A_j = t_j - \sum_{y \text{ crossing at } j} \frac{1}{d_y(j)} + a_j.
```

- Since $`A_{j+1} - A_j = I_j = \mathbf{1}\{\text{recipient at step } j > m\}`$, telescoping $`Y_{j+1} - Y_j = (I_j - p_j)/(j+1) + (p_j - Y_j)/(j+1)`$, using $`Y_N \to R_m`$ a.s. (from unit mass $`\sum Z_v = 1`$, baseline Theorem 9), conditional $`L^2`$ convergence of $`M_n = \sum_{j \ge n}(I_j - p_j)/(j+1)`$ ($`\mathbb{E}[M_n \mid \mathcal{F}_n] = 0`$), and $`q_n - Y_n = A_n/(n-2) - A_n/n = 2 q_n / n`$ gives $`R_m - q_n = M_n + V_n - 2 q_n / n`$ identically.

2. **Unconditional finite bound $`\mathbb{E}[k_j \mid \mathcal{F}_n] \le 12 k_n`$ and $`24 k_n / n`$ (`ADDENDUM-CUTOFF-THRESHOLD.md` eq. (8)):**
   - At any size $`n \ge 3`$ and $`2 \le m \le n`$, every early leaf ($`v \le m, d_v(n) = 1`$) is attached to one of the $`k = k_n \ge 1`$ early current cores (seed edge $`1`$–$`2`$ ensures $`k \ge 1`$).
   - Applying `SMALL-POWERS.md`'s autonomous lower leaf stock $`J`$ (started at each early host $`x`$'s full current leaf count $`L = L_{\mathrm{early}} + L_{\mathrm{late}}`$) covers all promotions of early leaves of $`x`$. The exponential Doob maximal bounds $`\mathbb{P}(\inf_s e^{-s} Y(s) \le L/2) \le (\sqrt{e}/2)^L`$ and $`\mathbb{P}(\sup_s e^{-s}(1+I(s)) \ge L/4) \le 2 e^{-L/8}`$, together with $`\int_0^\infty \min(1, 4 e^{-s})\,ds = 1 + \log 4`$, bound the expected number of promoted old leaves per host by

```math
1 + \log 4 + \sup_{L \ge 0} L\big[(\sqrt{e}/2)^L + 2 e^{-L/8}\big] \le 1 + \log 4 + \frac{1}{e(\log 2 - 1/2)} + \frac{16}{e} = 10.177024\ldots < 11
```

(and integer maximum over $`L \ge 0`$ is $`10.033364\ldots < 11`$).
- Adding the $`k`$ early cores already present gives $`\mathbb{E}[k_j \mid \mathcal{F}_n] \le \mathbb{E}[K_m \mid \mathcal{F}_n] \le 12 k`$ unconditionally for every $`j \ge n`$ and every tree $`T_n`$.
- A new crossing core edge can only arise when an early core selects a late leaf (promoting an early leaf gives that newly promoted early core only an early parent and a late **leaf**, adding zero crossing core edges); thus $`\mathbb{E}[t_{j+1} - t_j \mid \mathcal{F}_j] \le k_j / j`$ and $`\mathbb{E}[t_j - t \mid \mathcal{F}_n] \le 12 k \sum_{i=n}^{j-1} 1/i`$.
- By Tonelli's theorem on nonnegative terms:

```math
\sum_{j=n+1}^\infty \frac{1}{j(j+1)} \sum_{i=n}^{j-1} \frac{1}{i} = \sum_{i=n}^\infty \frac{1}{i} \sum_{j=i+1}^\infty \left(\frac{1}{j} - \frac{1}{j+1}\right) = \sum_{i=n}^\infty \frac{1}{i(i+1)} = \frac{1}{n},
```

and $`\sum_{j \ge n} \frac{1}{j(j+1)} = \frac{1}{n}`$. Combining yields $`\mathbb{E}[\sum_{j \ge n}(t_j - t + k_j)/[j(j+1)] \mid \mathcal{F}_n] \le 12 k / n + 12 k / n = 24 k / n`$ exactly.
- Splitting $`V_n`$ into old crossings ($`j \le n^2`$ and $`j > n^2`$, where the tail weight is $`\sum_{j > n^2} \frac{n}{j(j+1)} = \frac{n}{n^2+1} < \frac{1}{n}`$ and each integrand deviation is $`\le 1/2 \le 1`$) plus new crossings (each in $`[0, 1]`$) and $`a_j \in [0, k_j]`$ yields eq. (6) of `CRITICAL-EXCESS-DRIFT.md`: $`\mathbb{E}[|n V_n / t - (1-b)| \mid \mathcal{F}_n] \le \mathbb{E}[B_n \mid \mathcal{F}_n] + 1/n + 24 k_n / t`$.

---

### 2.3 Pillar 3: Conditional Concentration (3), Martingale CLT at Critical Scale, and Gaussian Centering (1)

1. **Direct audit of `ADDENDUM-CUTOFF-THRESHOLD.md` eq. (6) at the critical scale:**
   - Importantly, `CRITICAL-EXCESS-DRIFT.md` does **not** rely on `ADDENDUM-GAUSSIAN.md` (which only proved sufficiency under $`m^{1+\mu}(\log n)^2 = o(n)`$) or `ADDENDUM-DIRICHLET.md`; it uses eq. (6) of `ADDENDUM-CUTOFF-THRESHOLD.md`, which proves $`\mathrm{law}(M_n / \sqrt{q_n/n} \mid \mathcal{F}_n) \to N(0, 1)`$ a.s. for **every** deterministic $`m \to \infty, m = o(n)`$ (including the critical scale $`m^{1+\mu}(\log n)^2 / n \to \lambda`$).
   - We verified every step of that conditional characteristic-function proof:
     - By `FINITE-TREE.md`, $`\sup_{j \ge n} |p_j - R_m| \le \sup_{j \ge n} \mathrm{TV}(p_j, Z) = O_\omega(n^{\mu - 1})`$ (since $`\mu < 1`$), while $`q_n \sim R_m \sim c X m^{\mu - 1}`$, giving $`\sup_{j \ge n} |p_j / q_n - 1| \to 0`$, $`q_n \to 0`$, and $`n q_n \to \infty`$ a.s.
     - Dominating the exit indicator $`\mathbf{1}_{E_n(\varepsilon)}`$ by the decreasing full-path tail $`\mathbf{1}\{\exists r \ge N : E_r(\varepsilon)\}`$ and applying Lévy's upward martingale convergence theorem gives $`\mathbb{P}(E_n(\varepsilon) \mid \mathcal{F}_n) \to 0`$ a.s.
     - At the pre-jump stopping time $`\tau = \inf\{j \ge n : |p_j/q_n - 1| > \varepsilon\}`$ ($`\{\tau > j\} \in \mathcal{F}_j`$), replacing post-$`\tau`$ increments by independent centered $`\mathrm{Bernoulli}(q_n)`$ variables preserves the sum on $`E_n(\varepsilon)^c`$ and keeps $`|\tilde{p}_j - q_n| \le \varepsilon q_n`$ everywhere. Using the $`1`$-Lipschitz bound on $`p(1-p)`$, for $`q_n \le 1/2`$ and $`\varepsilon \le 1/4`$ we get $`|\tilde{v}_j - v_j| \le 2 \varepsilon v_j`$ for $`v_j = n(1-q_n)/(j+1)^2`$, with $`\sum_{j \ge n} v_j \to 1`$, $`\max_{j \ge n} v_j \le 1/n`$, and $`|\tilde{X}_j| \le 1/\sqrt{n q_n}`$.
     - Cubic Taylor expansion of $`e^{i u \tilde{X}_j}`$, telescoping against $`\exp(-u^2 v_j / 2)`$ with unit modulus of the past factor, passing to $`N \to \infty`$ by conditional $`L^2`$ convergence, and applying Lévy/Pólya (either via regular conditional distributions on the single probability-1 event independent of $`u`$, or on countable dense $`u`$ with tightness from eq. (4) $`\mathbb{E}[M_n^2 \mid \mathcal{F}_n] \le 2\mathbb{E}[R_m \mid \mathcal{F}_n]/n`$) yields $`\sup_z |\mathbb{P}(M_n / \sqrt{q_n/n} \le z \mid \mathcal{F}_n) - \Phi(z)| \to 0`$ a.s.

2. **Full-path majorant for (3) and conditional Slutsky without independence for (1):**
   - $`B_n \in [0, 1/2]`$ is purely $`\mathcal{F}_\infty`$-measurable and tends to $`0`$ a.s. by (5). Setting $`M_N^* = \sup_{r \ge N} B_r \in [0, 1/2]`$, Lévy's upward theorem gives $`\limsup_{n \to \infty} \mathbb{E}[B_n \mid \mathcal{F}_n] \le M_N^*`$ a.s. for each $`N`$, hence $`\mathbb{E}[B_n \mid \mathcal{F}_n] \to 0`$ a.s. Since $`k_n \le K_m = O_\omega(C_m)`$ and $`t \sim h C_m L`$, eq. (6) implies (3): $`\mathbb{E}[|n V_n / t - (1 - b)| \mid \mathcal{F}_n] \to 0`$ a.s.
   - Algebraic constant verification:

```math
\delta_n := (1 - b)\frac{t}{\sqrt{n q_n}} \sim (1 - b) \frac{(M_0/\mu) X m^\mu \cdot \frac{\mu}{1+\mu}\log n}{\sqrt{n c X m^{\mu-1}}} = (1 - b)\frac{M_0}{1+\mu}\sqrt{\frac{X}{c}}\sqrt{\frac{m^{1+\mu}(\log n)^2}{n}} \to (1 - b)\frac{M_0}{1+\mu}\sqrt{\frac{\lambda X}{c}} = \delta(X).
```

- Writing $`(R_m - q_n)/\sqrt{q_n/n} = M_n/\sqrt{q_n/n} + \delta_n + \Delta_n`$ with $`\delta_n \in \mathcal{F}_n`$ and $`\mathbb{E}[|\Delta_n| \mid \mathcal{F}_n] \le \frac{t}{\sqrt{nq_n}}\mathbb{E}[|n V_n/t - (1-b)| \mid \mathcal{F}_n] + 2\sqrt{q_n/n} \to 0`$ a.s.:
- Since $`\mathbb{E}[M_n \mid \mathcal{F}_n] = 0`$, the conditional mean claim $`(\mathbb{E}[R_m \mid \mathcal{F}_n] - q_n)/\sqrt{q_n/n} = \delta_n + \mathbb{E}[\Delta_n \mid \mathcal{F}_n] \to \delta(X)`$ a.s. (and $`\mathbb{E}[R_m \mid \mathcal{F}_n] - q_n = (1-b)t/n + o_\omega(\sqrt{q_n/n})`$) follows in one line.
- For uniform Kolmogorov convergence (1), the $`\pm \eta`$ event sandwich uses only Markov's inequality on $`\mathbb{E}[|\Delta_n| \mid \mathcal{F}_n]`$ and the $`1/\sqrt{2\pi}`$-Lipschitz continuity of $`\Phi`$, requiring **no independence** between $`M_n`$ and $`\Delta_n`$:

```math
\sup_{z \in \mathbb{R}} \left|\mathbb{P}\left(\frac{R_m - q_n}{\sqrt{q_n/n}} \le z \;\middle|\; \mathcal{F}_n\right) - \Phi(z - \delta_n)\right| \le D_n + \frac{\eta}{\sqrt{2\pi}} + \frac{\mathbb{E}[|\Delta_n| \mid \mathcal{F}_n]}{\eta},
```

which vanishes a.s. (and likewise with $`\delta(X)`$ in place of $`\delta_n`$ since $`|\delta_n - \delta(X)| \to 0`$ a.s.).

---

### 2.4 Pillar 4: Exact Multiplier Enclosure in `CRITICAL-EXCESS-EVIDENCE.md`

1. **Upper payoff $`d^{-\ell - 1}`$:** For any finite rooted inert tree with root degree $`d \ge 2`$ (including exterior slot) and $`\ell \ge 0`$ ordinary leaf neighbors, a linear birth-death chain from $`\ell`$ with birth rate $`j`$ and death rate $`j/d`$ couples below the root's leaf stock (using the exact Pillar-1 ordering, since choosing the exterior slot promotes no leaf). Extinction has probability $`d^{-\ell}`$ (trivially $`1`$ for $`\ell = 0`$), and survival forces $`D_{\mathrm{root}}(\infty) = \infty \implies 1/D_{\mathrm{root}}(\infty) = 0`$; on extinction $`1/D_{\mathrm{root}}(\infty) \le 1/d`$, giving $`\mathbb{E}_T[1/D_{\mathrm{root}}(\infty)] \le d^{-\ell - 1}`$.
2. **Lower payoff and protection recurrence:**
   - For each ordinary neighbor $`z`$ of the root (with $`L_z`$ degree-$`1`$ neighbors; the root has degree $`d \ge 2`$ and is excluded), tracking lower stock $`j`$ (birth rate $`j`$; at each rate-$`1`$ ring of $`z`$, kill with prob $`1/(j+1)`$ and otherwise decrement $`j`$; fail at $`0`$) catches every physical send $`z \to \mathrm{root}`$ because $`D_z \ge j + 1`$ until failure, and ordering the root first places $`[0, 1/D_z) \subseteq [0, 1/(j+1))`$.
   - Stopping at $`\{2, \dagger\}`$ with $`f(j) = 1/(j-1)`$ ($`j \ge 3`$) and $`f(2) = f(\dagger) = 1`$, the generator satisfies

```math
(\mathcal{L} f)(j) = j[f(j+1) - f(j)] + \frac{1 - f(j)}{j+1} + \frac{j}{j+1}[f(j-1) - f(j)] = \frac{6 - 2j}{(j+1)(j-1)(j-2)} \le 0 \quad (j \ge 3).
```

Nonexplosion (Yule domination) and bounded stopping give failure probability $`\le 1/(j - 1)`$ for all $`j \ge 2`$.
- Below $`N = 24`$, the finite absorbing chain exits $`\{1, \dots, N-1\}`$ a.s. and its exact hitting-before-failure probability $`p_j^{(N)}`$ satisfies $`(j+1) p_j^{(N)} = j p_{j+1}^{(N)} + \frac{j}{j+1} p_{j-1}^{(N)}`$; multiplying by $`1 - 1/(N-1)`$ gives $`g_j \le \mathbb{P}_j(z \text{ never sends to root})`$.
- Union-bounding over present ordinary root neighbors (no independence needed) gives $`\mathbb{P}_T(D_{\mathrm{root}}(\infty) = d) \ge \max(0, 1 - \sum_z (1 - g_{L_z}))`$, since no new root neighbor can arrive without an existing ordinary neighbor sending first and the inert slot never sends.
3. **Transition kernel, backwards intersection, and convergence as cap $`\to \infty`$:**
   - Total size-changing jump rate at ordinary size $`s \ge 2`$ is $`s - 1/d_0 > 0`$, with exact recipient rate $`\sum_{v \sim_{\mathrm{ord}} w} 1/d_v`$; dynamics lump cleanly onto rooted isomorphism classes; backwards intersection at each node is valid by the strong Markov property at the first jump.
   - As size cap $`k \to \infty`$: on $`\{D_{\mathrm{root}}(\infty) = \infty\}`$, both payoffs are $`\le 1/D_{\mathrm{root}}(\tau_k) \to 0`$. On $`\{D_{\mathrm{root}}(\infty) < \infty\}`$, every fixed ordinary neighbor $`z`$ of the root must have $`D_z(\infty) = \infty`$ a.s. (otherwise $`z`$ would retain constant positive send rate $`1/d_z`$ to the root forever and send infinitely often). For every label $`z`$ and every $`d \ge 2`$, the fresh-leaf birth-death comparison at $`z`$'s first hit of degree $`d`$ bounds $`\mathbb{P}(D_z(\infty) = \infty, L_z(t) \not\to \infty) \le 1/d`$; sending $`d \to \infty`$ and union-bounding over countably many labels yields $`L_z(\tau_k) \to \infty`$ and $`\ell(\tau_k) = 0`$ eventually a.s. on $`\{D_{\mathrm{root}}(\infty) < \infty\}`$. Thus both raw payoffs converge a.s. to $`1/D_{\mathrm{root}}(\infty)`$, and bounded convergence closes the interval as $`k \to \infty`$.

---

## 3. Center-Specificity Analysis and Independent Identities

### 3.1 Independent Exact Finite-Tree Center-Difference Identities and Critical-Scale Corollaries

From our first-principles derivation in §2.2, **every** edge-seeded tree $`T_n`$ ($`n \ge 3`$) and **every** cutoff $`2 \le m \le n`$ satisfies the three exact finite-tree identities linking the observed excess-degree center $`q_n = A_n/(n-2)`$, the star-weight center $`\pi_n = (A_n + t_n)/n`$, and the next-recipient center $`p_n = \mathbb{P}(\text{next recipient} > m \mid \mathcal{F}_n)`$:

```math
n(p_n - q_n) = t_n - h_n + a_n - 2 q_n, \qquad n(\pi_n - q_n) = t_n - 2 q_n, \qquad n(\pi_n - p_n) = h_n - a_n, \qquad (\text{I1--I3})
```

where $`h_n = \sum_{y \text{ crossing at } n} 1/d_y(n)`$ and $`a_n = \sum_{x \le m, d_x \ge 2} l_x / d_x`$, together with the exact one-step drift identity

```math
\mathbb{E}[Y_{n+1} - Y_n \mid \mathcal{F}_n] = \frac{t_n - h_n + a_n}{n(n+1)}. \qquad (\text{I4})
```

Combined with equations (2), (3), (5) of `CRITICAL-EXCESS-DRIFT.md` and eq. (6) of `ADDENDUM-CUTOFF-THRESHOLD.md`, these identities yield **three new limiting consequences** at the critical scale $`m^{1+\mu}(\log n)^2 / n \to \lambda \in (0, \infty)`$ that prove directly why the shift $`+\delta(X)`$ is specific to $`q_n`$:

1. **Observable center-separation limits:** Since $`h_n / t_n \to b`$ a.s. (eq. (5) at $`j = n`$) and $`0 \le a_n \le k_n = o_\omega(t_n)`$,

```math
\frac{n(p_n - q_n)}{t_n} \to 1 - b, \qquad \frac{n(\pi_n - p_n)}{t_n} \to b, \qquad \frac{p_n - q_n}{\sqrt{q_n/n}} \to \delta(X) \quad \text{a.s.}
```

2. **Star-weight center $`\pi_n`$ at the critical scale (different sign and coefficient):** Since $`R_m - \pi_n = M_n + (V_n - t_n/n)`$ and $`\pi_n \sim q_n`$ a.s.,

```math
\sup_{z \in \mathbb{R}} \left|\mathbb{P}\left(\frac{R_m - \pi_n}{\sqrt{\pi_n/n}} \le z \;\middle|\; \mathcal{F}_n\right) - \Phi\left(z + b \frac{M_0}{1+\mu}\sqrt{\frac{\lambda X}{c}}\right)\right| \to 0 \quad \text{a.s.},
```

i.e., the **same** scaling $`m^{1+\mu}(\log n)^2/n \to \lambda`$ is critical for $`\pi_n`$, but with **negative** shift $`-\frac{b}{1-b}\delta(X) = -b\frac{M_0}{1+\mu}\sqrt{\lambda X/c}`$ instead of $`+\delta(X)`$.
3. **Next-recipient center $`p_n`$ at the critical scale (zero shift — not critical for $`p_n`$!):** Since

```math
R_m - p_n - M_n = \sum_{j \ge n} \frac{1}{j(j+1)} \left[\sum_{y \text{ crossing at } n} \left(\frac{1}{d_y(n)} - \frac{1}{d_y(j)}\right) + \sum_{\text{new } y} \left(1 - \frac{1}{d_y(j)}\right) + (a_j - a_n)\right],
```

the uniform supremum (5) bounds $`\frac{1}{t_n}\sum_{y \text{ crossing at } n}(1/d_y(n) - 1/d_y(j)) \le 2 B_n`$ for all $`n \le j \le n^2`$, yielding $`\frac{n}{t_n}\mathbb{E}[|R_m - p_n - M_n| \mid \mathcal{F}_n] \le 2\mathbb{E}[B_n \mid \mathcal{F}_n] + 1/n + 25 k_n/t_n \to 0`$ a.s., hence $`\mathbb{E}[|R_m - p_n - M_n|/\sqrt{q_n/n} \mid \mathcal{F}_n] \to 0`$ a.s. and (using $`p_n \sim q_n`$ a.s.):

```math
\sup_{z \in \mathbb{R}} \left|\mathbb{P}\left(\frac{R_m - p_n}{\sqrt{p_n/n}} \le z \;\middle|\; \mathcal{F}_n\right) - \Phi(z)\right| \to 0 \quad \text{a.s.}
```

Thus the scaling $`m^{1+\mu}(\log n)^2 / n \to \lambda`$ and positive shift $`+\delta(X)`$ apply **strictly** to the observed excess-degree center $`q_n`$ (and $`Y_n = A_n/n = q_n - 2 q_n / n`$), NOT to $`\pi_n`$ (negative shift $`-b/(1-b)\cdot\delta(X)`$) or $`p_n`$ (centered $`N(0, 1)`$, zero shift).

### 3.2 Independent Analytical Identities for $`b = \mathbb{E}_J[1/D_{\mathrm{root}}(\infty)]`$ and the Protection Chain

1. **Exact 1-step and 2-step Bellman identities for the inert virgin $`J`$:** Enumerating the first two size-changing transitions of $`J`$ directly by hand gives

```math
b = \frac{2}{3} b(S_2) + \frac{1}{3} b(P_3) = \frac{1}{2} b(S_3) + \frac{7}{30} b(B_4) + \frac{1}{5} b(Y_4) + \frac{1}{15} b(P_4),
```

where $`S_k`$ is the $`k`$-leaf star at the root, $`P_k`$ is the path of $`k`$ vertices, $`B_4`$ is the size-4 broom (root with 1 leaf and 1 path-of-2), and $`Y_4`$ is the root attached to a degree-3 child with 2 leaves. Applying the raw upper payoff at size 3 immediately yields $`b \le \frac{2}{3}\cdot\frac{1}{27} + \frac{1}{3}\cdot\frac{1}{2} = \frac{31}{162} \approx 0.191358`$. On the lower side, one step from $`Y_4`$ produces $`Y_5`$ (root degree $`2`$, unique child with $`L_1 = 3`$ leaves) with probability $`\frac{1}{5}\cdot\frac{5}{7} = \frac{1}{7}`$, where the generator bound $`1 - 1/(3-1) = 1/2`$ yields the five-line hand enclosure

```math
\frac{1}{28} \le b \le \frac{31}{162} \qquad (0.035714 \le b \le 0.191358).
```

2. **Exact two-sided analytical sandwich for the killed protection chain:** Rewriting $`(j+1) w_j = j w_{j+1} + \frac{j}{j+1} w_{j-1}`$ ($`w_0 = 0, w_1 = 1`$) in forward differences $`\Delta_j = w_{j+1} - w_j`$ gives $`\Delta_j = \frac{1}{j+1}\Delta_{j-1} + \frac{1}{j(j+1)} w_j > 0`$. Hence $`w_j \uparrow w_\infty < \infty`$, and summing $`\Delta_r \ge w_j / [r(r+1)]`$ over $`r \ge j`$ yields $`w_\infty - w_j \ge w_j / j`$, i.e., $`g_j^{(\infty)} = w_j / w_\infty \le j/(j+1)`$. Together with Astra's generator bound, the exact infinite-horizon protection probability satisfies the sharp two-sided sandwich

```math
\frac{j - 2}{j - 1} \le g_j^{(\infty)} \le \frac{j}{j + 1} \qquad \text{for all } j \ge 2,
```

proving that Astra's asymptotic boundary value $`1 - 1/(N-1)`$ is sharp to first order ($`1/N + O(N^{-2})`$).

---

## 4. Independent Zero-Shared-Code Computational Evidence

Rather than merely replaying `check_critical_drift.py`, we wrote and executed an independent zero-shared-code exact rational checker (`audit_critical_16f0c12.py` (reviewer scratch script, not exported)) using:
- Exhaustive enumeration of all $`5,912`$ edge-seeded labeled trees ($`n = 3..8`$) and all $`40,318`$ `(tree, m)` pairs for (I1)–(I4) plus a planted-defect control (dropping $`a_n`$).
- Balanced-parenthesis (Dyck-word) canonical rooted tree encoding, direct `(initiating source, neighbor)` transition-rate enumeration, exact Thomas tridiagonal elimination for the absorbing protection BVP ($`N = 24, 100, 200`$), and bottom-up BFS dynamic programming through size cap $`10`$ (confirming OEIS A000081 counts $`1, 2, 4, 9, 20, 48, 115, 286, 719`$; total $`1,204`$ states), plus a planted-defect generator control ($`f(j) = 1/j`$ fails with positive generator).

Command run (exit 0):
```text
python3 -uB audit_critical_16f0c12.py  # reviewer scratch run
```
Actual output:
```text
[1] Finite-tree identities PASSED: trees=5912, (tree,m)-pairs=40318, 1-step-drifts=718, nonzero(p-q)=34192, nonzero(pi-p)=32730, planted_defect_caught=True
[2] Independent inert enclosure PASSED: states(2..10)=1204, N=24 cap=10 -> b in [0.0652016, 0.1166501] (1-b in [0.8833499, 0.9347984]), N=100 cap=10 -> b in [0.0653235, 0.1166501]
    exact float b(N=24,cap=10): [0.0652016469, 0.1166500430]
    exact float b(N=100,cap=10): [0.0653235644, 0.1166500430]
    exact 1-step/2-step/Y5 hand identities verified: b >= 1/28 = 0.035714, b <= 31/162 = 0.191358
[3] Analytic constants PASSED: eq.(8) separate_sup=10.177024 < 11, integer_max=10.033364 < 11; protection generator identity verified for j=3..200
```

Our zero-shared-code reconstruction certifies Astra's exact outward interval $`b \in [0.0652016, 0.1166501]`$ ($`1 - b \in [0.8833499, 0.9347984]`$) identically (and increasing the protection BVP cap from $`N = 24`$ to $`N = 100`$ at size cap $`10`$ tightens the lower endpoint slightly to $`0.0653235`$).

---

## 5. Complete Inherited Input List & Review-Boundary Ledger

### 5.1 Audited Directly in This Review (`HOLDS`)
1. `CRITICAL-EXCESS-DRIFT.md` (`16f0c12`, SHA-256 `4c27a043e261f9c2f3e9db1dd529dfce70cce0c8f548e875039fb85e262fa4e8`): full proof of (1)–(6).
2. `CRITICAL-EXCESS-EVIDENCE.md` (`16f0c12`, SHA-256 `395997919be9a72d7b7e769a076430f44fb19b93cb51490c36b6a1b2fc21af02`): payoffs, protection recurrence, convergence proof, and exact rational interval $`[0.0652016, 0.1166501]`$ (independently rebuilt without shared code).
3. `ADDENDUM-CUTOFF-THRESHOLD.md` (`b404108`, SHA-256 `ab84cf0e09826c43a56e020387ce7256849acc7c1ab36976f135d97f17405bb3`), **equations (2), (3), (4), (6), and (8) only**: exact finite-tree drift identity, martingale/drift decomposition, martingale conditional $`L^2`$ bound, conditional martingale CLT for all deterministic $`m \to \infty, m = o(n)`$, and $`\mathbb{E}[k_j \mid \mathcal{F}_n] \le 12 k_n \implies 24 k_n / n`$ bound.

### 5.2 Previously Cross-Family Reviewed by Neon (`9be4767` — `late-mass-neon/REVIEW.md` (prior snapshot review, not exported), `HOLDS`)
4. `PROOF.md`: $`R_m \sim c X m^{\mu - 1}`$ a.s., $`c = M_1/(1-\mu)`$ (also independently proved in Neon's `PROOF.md`).
5. `FINITE-TREE.md`: $`\mathrm{TV}(p_n, Z) = O_\omega(n^{\mu-1})`$, $`\mathrm{TV}(q_n, Z) = O_\omega(n^{\mu-1})`$, and $`q_n \sim R_m`$ a.s. for all diverging $`m = o(n)`$.
6. `SMALL-POWERS.md`: autonomous lower leaf-stock comparison and per-host old-leaf promotion bound giving $`K_m = O_\omega(C_m)`$ (and used in `ADDENDUM-CUTOFF-THRESHOLD.md` eq. (8)).
7. `POSITIVE-LABELS.md`: $`H_m / m^\mu \to X M_0 / \mu`$ a.s., hence $`H_m / C_m \to h = M_0 / \mu`$ a.s.

### 5.3 Inherited Round 1 / Round 2 Public Baseline Inputs (`round2/LEDGER.md` at pinned base `c250a1f2fdd53ced91b7b485e94ac8a16106b9c7`)
8. **Theorem 9 / Astra 7 (Tier A):** Almost-sure unit total limiting share $`\sum_{v \ge 1} Z_v = 1`$ (used in $`Y_N \to R_m`$).
9. **Astra 6 (Tier A):** Core growth scaling $`C_n / n^\mu \to X \in (0, \infty)`$ a.s., with $`\mu \in (1/2, 1)`$ (Round 1 / Gold2 certified).
10. **Astra 10 (Tier B):** Uniform empirical core local limit $`\pi_n \to \pi`$ a.s. and properness of the limiting core root-degree law ($`\pi_{\mathrm{core}}(d > D) \to 0`$).
11. **Astra 22 (Tier B):** Simultaneous marked rectangle theorem for $`(v/m, v Z_v)`$ with positive birth-mark intensity $`\nu`$ and $`M_0 = \nu((0, \infty)) > 0`$.
12. **Physical degree dichotomy (`Neon 2(i)` dichotomy Tier B; `late/README.md` §3 / `hub-children.md` under Astra 18b Tier C, independently cold-audited in Neon's `PROOF.md` §1.2):** Simultaneously for all labels $`v \ge 1`$, $`D_v(\infty) = \infty \iff Z_v > 0`$.
13. **Inert virgin definition (`forest-age/inert/README.md`):** Only the elementary Markov definition of the inert virgin $`J`$ (with exterior slot retained and changing rate $`s - 1/d_0`$) is used; `CRITICAL-EXCESS-DRIFT.md` does **not** import the Tier-C `Astra 20b` inert core-scaling or TV theorems.

### 5.4 Explicitly Unreviewed Continuation Material (NOT Depended Upon by `16f0c12`)
- `ADDENDUM-DIRICHLET.md`, `ADDENDUM-GAUSSIAN.md`, `ADDENDUM-FIRST-HIT.md`, `ADDENDUM-EVIDENCE.md`, `EXTREMES.md`, and the converse/star-weight/diagnostic sections (5), (7), (9), (10), (11) of `ADDENDUM-CUTOFF-THRESHOLD.md` are **outside** the scope of this `HOLDS` verdict and are not required for any claim in `CRITICAL-EXCESS-DRIFT.md` or `CRITICAL-EXCESS-EVIDENCE.md`.
