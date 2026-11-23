# Independent Referee Report


The editing of this report for publication shortened two third-party quotations. The mathematical assessments are the referee’s.
**R3 · Gemini 4 Neon · xhigh reasoning.** Independent cold report on §§2–4, 6–8, production note and bibliography at `526f350`.

**Manuscript:** *Growth of the non-leaf population in the random friend tree* (`paper/friend-tree-exponent.tex`, commit `526f350`)
**Target venue standard:** *Annals of Applied Probability* / *Electronic Journal of Probability*
**Recommendation:** **Minor revision** (`minor`)

---

## 1. Executive Summary and Independent Claim Summary

### 1.1 Main claims established by the manuscript
Let $`(T_n)_{n\ge 2}`$ be the random friend tree (one-step random-walk attachment / isotropic redirection tree): starting from a single edge $`T_2`$, given $`T_n`$, select a uniform vertex $`V_n \in V(T_n)`$, select a uniform neighbour $`W_n`$ of $`V_n`$, and attach a new leaf to $`W_n`$ to obtain $`T_{n+1}`$. Let $`C_n = N_{\ge 2}(n) = |\{v \in V(T_n) : \deg_{T_n}(v) \ge 2\}|`$ denote the number of non-leaves (core vertices), and put $`P_a(n) = \prod_{j=3}^{n-1}(1 + a/j)`$ with $`P_a(3) = 1`$.

The manuscript establishes three main results:

1. **Existence of a common deterministic growth exponent (Theorem 1.1, §2).** There exists a deterministic constant $`\mu \in [0, 1]`$ such that
   ```math
   \lim_{n\to\infty}\frac{\log \mathbb{E}[C_n]}{\log n} = \mu
   \qquad\text{and}\qquad
   \lim_{n\to\infty}\frac{\log C_n}{\log n} = \mu \quad \text{almost surely.}
   ```
   This proof is purely analytic (no computer assistance) and embeds the process in continuous time with rate-$`1`$ Poisson clocks per vertex, proving almost-supermultiplicativity of $`F(t) = \mathbb{E}[C(t)]`$ by planting arbitrary finite seeds at eligible core vertices, de-Poissonising to $`\mathbb{E}[C_n]`$, and running conditionally independent site trials over disjoint reserved clock sets to lift the rate to almost-sure convergence.
2. **Sharp finite-$`n`$ two-sided product bracket via motif potentials (Theorem 1.2, first display, §3).** For every integer $`n \ge 3`$,
   ```math
   \frac{7{,}304{,}352{,}308}{8{,}763{,}936{,}965}\,P_{113/200}(n)
   \;\le\; \mathbb{E}[C_n] \;\le\;
   \frac{14{,}606{,}753{,}031}{2{,}989{,}767{,}809}\,P_{1133/2000}(n).
   ```
   Combined with Theorem 1.1, this yields $`113/200 = 0.5650 \le \mu \le 1133/2000 = 0.5665`$, separating $`\mu`$ from $`17/30`$ by at least $`1/6000`$. Independently of Theorem 1.1, dyadic Markov and Borel–Cantelli yield $`\limsup_{n\to\infty} (\log C_n)/\log n \le 1133/2000`$ almost surely, and the upper product bounds $`\mathbb{E}[N_{=k}(n)]`$ for every $`k \ge 2`$.
3. **Independent finite-$`n`$ upper bound below $`17/30`$ via banded history-dependent potentials (Theorem 1.2, second display, §4).** For every integer $`n \ge 3`$,
   ```math
   \mathbb{E}[C_n] \;\le\; \frac{1{,}539{,}738{,}114{,}735}{3{,}908{,}512{,}529}\,P_{5657/9983}(n),
   ```
   where $`5657/9983 = 17/30 - 1/299490 < 17/30`$. This provides an independent confirmation of $`\limsup_{n\to\infty} (\log \mathbb{E}[C_n])/\log n < 17/30`$ and $`\limsup_{n\to\infty} (\log C_n)/\log n \le 5657/9983`$ almost surely.

In addition, §6 summarizes supplementary results (an almost-sure lower transfer theorem yielding $`\liminf (\log C_n)/\log n \ge 543/1000`$ and $`11/20`$ via cross-charging, a polynomial-scale lower bound $`\liminf C_n/n^{543/1000} > 0`$ a.s., and the uniform cumulative mean degree-tail limit $`\lim_{n/k\to\infty} \log \mathbb{E}[N_{\ge k}(n)] / \log(n/k) = \mu`$ over $`k\ge 2`$, plus fixed-$`k`$ a.s. convergence $`\log N_{\ge k}(n)/\log n \to \mu`$ via Addario-Berry et al., Theorem 4.8), and §7 presents discrete and continuous-time numerical measurements (explicitly kept separate from all proofs).

### 1.2 Overall mathematical verdict
All three headline claims—**existence and a.s. limit (Theorem 1.1)**, **sharp motif bracket $`[113/200, 1133/2000]`$ (Theorem 1.2)**, and **independent banded upper bound $`5657/9983 < 17/30`$ (Theorem 1.2)**—**hold**. Every step of the analytic proofs in §§2–4 (and the referenced supplementary proofs in §6) has been audited, all finite-certificate verifiers have been rebuilt with undefined-behaviour sanitizers and re-executed cleanly, and independent small-case exact-rational computations corroborate the generator identities, bounds, state counts, and simulation summaries. No proof step was found to be false or beyond straightforward self-contained repair.

### 1.3 Top three issues for revision
1. **Self-containment of the §2 coupling and simultaneous-sites construction (`paper/friend-tree-exponent.tex`, lines 102–111 and 133–138).** Three subtle admissibility points in the paper text (executing the triggering ring before switching in §2.2; proving that the re-driven promoted leaf $`v_i^*`$ in §2.4 Stage 3 cannot collide with any other site's reserved vertices and clarifying failure-reversion before Stage 3; defining the free reference histories and exterior tests unconditionally at every site so $`(Y_i)_{i=1}^m`$ are manifestly i.i.d.) are stated too tersely in the manuscript (though expanded in `proofs/existence.md` and `proofs/almost-sure-exponent.md`). They should be made explicit in §2 so the main article stands alone at *AAP/EJP* rigor.
2. **Bibliography, missing literature, and archival stability of the supplement (`paper/friend-tree-exponent.tex`, lines 1307–1331).** Reference `[AB]` (Addario-Berry et al.) is cited only as a 2024 arXiv preprint despite being published in *Random Structures & Algorithms* **69**(1), e70082 (2026); several key foundational and closely related works (notably Evans–Saramäki 2005, Krapivsky–Redner 2001 and 2024, Gabel–Krapivsky–Redner 2013/2014, Vázquez 2003, and continuous-time branching/urn literature such as Janson 2004 and Holmgren–Janson 2017) are omitted; and `[supplement]` points to an earlier GitHub revision (`c4653036b993`, which does not contain `paper/`) rather than a permanent archival DOI synchronized with the manuscript.
3. **Minor formal, notational, and verification-hygiene cleanups (`paper/friend-tree-exponent.tex`, lines 236, 631–659, 1128–1190, 1282–1306; `verify/motif-bracket/cert_verify_gmp.c`, line 284; `verify/run_all.sh`, line 49).** Fix the undefined symbol $`\pi^{\mathrm{end}}`$ on line 236 (should be $`\zeta_s`$); resolve symbol overloading across §§2–4 (e.g., $`H, h, b, A, T`$) and clarify the $`N_k`$ vs. $`N_{\ge k}`$ notation difference when linking to `proofs/degree-tail.md`; add the missing programmatic header check `D0 <= KF + 1` to `cert_verify_gmp.c` and enable UBSan on `cert_verify_gmp.c` in `verify/run_all.sh`; clarify the editorial boundary of §6 (results whose proofs reside only in the supplement); and state the human author's end-to-end verification responsibility in the production disclosure.

---

## 2. Step-by-Step Proof Audit of §§2–4 and §§6–8

Each logical step is anchored to line numbers and equation labels in `paper/friend-tree-exponent.tex` (with supplementary anchors where §6 relies on `proofs/*.md`) and graded using the four-level rubric: **holds**, **holds with a repair (supply proof)**, **not justified as written**, or **false (show failure)**.

### 2.1 Section 2: Growth from finite seeds (`paper/friend-tree-exponent.tex`, lines 63–161)

#### Continuous-time embedding, admissible couplings, and Yule moments (lines 63–79)
1. **Lines 64–65 (Poisson embedding, jump chain, non-explosion):**
   - *Claim:* Rate-$`1`$ Poisson clocks per vertex induce exponential holding times of rate $`j`$ at size $`j`$, independent of the discrete jump chain $`(T_j)_{j\ge 2}`$; the sum of holding times diverges a.s.
   - *Check:* Minimum of $`j`$ independent $`\mathrm{Exp}(1)`$ waits is $`\mathrm{Exp}(j)`$, independent of the uniform argmin $`V_j`$ and its uniform neighbour $`W_j`$. Writing holding times as $`e_j/j`$ with $`e_j \stackrel{\text{i.i.d.}}{\sim} \mathrm{Exp}(1)`$, $`\sum_{j\ge 2} (e_j-1)/j`$ has summable variances $`\sum j^{-2} < \infty`$ and converges in $`L^2`$ and a.s., while $`\sum_{j\ge 2} 1/j = \infty`$, yielding non-explosion a.s.
   - *Verdict:* **holds**.
2. **Lines 66–76 (Coupling admissibility and two-step choice):**
   - *Claim:* Conditional on the history after a ring or at a deterministic time, remaining waits are i.i.d. $`\mathrm{Exp}(1)`$; giving the correct conditional law of holding time, source, and neighbour at each ring determines the law on every bounded interval; predictable ranking and two-step choice with acceptance $`d/(d+1)`$ produce uniform neighbour choices.
   - *Check:* Memorylessness / strong Markov property of independent Poisson processes applies at ring times and finitely many deterministic stage boundaries; for two-step choice, each of the first $`d`$ neighbours has probability $`(1/d)\cdot d/(d+1) = 1/(d+1)`$ and the last has probability $`1 - d/(d+1) = 1/(d+1)`$.
   - *Verdict:* **holds**.
3. **Lines 78–79 (Yule families and first two moments):**
   - *Claim:* Starting from $`k`$ vertices, $`N(t)`$ is a sum of $`k`$ independent Yule families with $`p_j(t) = e^{-t}(1-e^{-t})^{j-1}`$ solving $`p_j' = (j-1)p_{j-1} - jp_j`$, $`p_1(0)=1`$, giving $`\mathbb{E}N(t) = ke^t`$ and $`\mathbb{E}N(t)^2 = k(k+1)e^{2t} - ke^t`$.
   - *Check:* Assigning each newborn's clock to the clan of the ringing source (not the attachment target) gives $`k`$ independent rate-$`1`$ pure-birth processes. Direct differentiation verifies $`p_j'(t) = e^{-t}(1-e^{-t})^{j-2}[-(1-e^{-t}) + (j-1)e^{-t}] = (j-1)p_{j-1}(t) - jp_j(t)`$. For geometric $`G`$ with $`q=e^{-t}`$, $`\mathbb{E}G = e^t`$ and $`\mathrm{Var}(G) = e^{2t}-e^t`$, so summing $`k`$ independent copies gives $`\mathbb{E}N(t)^2 = (ke^t)^2 + k(e^{2t}-e^t) = k(k+1)e^{2t}-ke^t`$.
   - *Verdict:* **holds**.

#### Lemma 2.1 (`lem:eligible`, lines 81–85)
4. **Line 84 ($`C=0`$ and $`C=1`$):**
   - *Claim:* Assertion is empty for $`C=0`$ and holds at the star centre for $`C=1`$.
   - *Check:* For $`C=1`$, the unique core vertex has $`c=0\le 2`$ and $`L = D \ge 2 \ge 1`$, hence $`E = 1 \ge 1/4`$.
   - *Verdict:* **holds**.
5. **Line 84 ($`C\ge 2`$: edge counting for isolated bare vertices and degree sum):**
   - *Claim:* The $`I`$ bare vertices ($`c=2, L=0`$) without bare neighbours satisfy $`2I \le C-1`$; the count $`n_{\ge 3}`$ of core vertices with $`c\ge 3`$ is at most the number $`n_1`$ of core-tree ends ($`c=1`$), each of which has a leaf; hence $`C \le E + I + E`$ and $`E \ge C/4`$.
   - *Check:* Core vertices induce a tree with $`C-1`$ edges. No two of the $`I`$ isolated bare vertices share an edge, so their $`2I`$ incident core edges are pairwise distinct, giving $`2I \le C-1`$. Summing $`c_u - 2`$ over the core tree gives $`-n_1 + \sum_{c_u\ge 3}(c_u - 2) = -2`$, so $`n_{\ge 3} \le n_1 - 2 < n_1`$. Every core end has $`c_u = 1`$ and $`D_u \ge 2`$, hence $`L_u \ge 1`$, so it is eligible ($`n_1 \le E`$). Partitioning core vertices gives $`C = E + I + n_{\ge 3} \le 2E + (C-1)/2`$, so $`E \ge (C+1)/4 > C/4`$.
   - *Verdict:* **holds**.

#### Lemma 2.2 (`lem:plant`, lines 86–98) and disjointness (line 99)
6. **Line 90 (Stage 1, unit 1):**
   - *Claim:* Keeping $`u`$ silent for one unit (and if $`u`$ is bare with bare neighbour $`w`$ and other neighbour $`z`$ of $`w`$, keeping $`z`$ silent and requiring $`w`$'s first ring to choose $`u`$) yields a tagged leaf at $`u`$ with $`c_u\le 2`$ and probability at least $`e^{-2}(1-e^{-1})/2`$.
   - *Check:* If $`u`$ has a leaf, $`u`$ staying silent (probability $`e^{-1} > e^{-2}(1-e^{-1})/2`$) prevents any promotion of its leaves and keeps $`c_u \le 2`$. If $`u`$ is bare, $`w`$ has neighbours $`\{u, z\}`$ only; while $`u, z`$ are silent, no leaf can attach to $`w`$, so at $`w`$'s first ring $`D_w = 2`$ and it chooses $`u`$ with probability $`1/2`$; the newborn leaf at $`u`$ cannot be promoted while $`u`$ stays silent. Independence of the three clocks gives $`e^{-2}(1-e^{-1})/2`$.
   - *Verdict:* **holds**.
7. **Line 90 (Stage 2, unit 2):**
   - *Claim:* Keeping $`u`$ silent and requiring $`\ge M`$ rings of the tagged leaf has conditional probability $`e^{-1}\mathbb{P}(\mathrm{Pois}(1)\ge M)`$ and yields $`L_u \ge M+1`$ and $`c_u \le 2`$.
   - *Check:* Distinct clocks give independence; while $`u`$ is silent the tagged leaf stays a leaf and every ring attaches a leaf to $`u`$, while exterior arrivals can only add leaves to $`u`$.
   - *Verdict:* **holds**.
8. **Line 92 (Stage 3, unit 3):**
   - *Claim:* Requiring one ring of $`u`$ with leaves ordered first and mark $`< M/(M+2)`$ promotes a leaf $`v`$; requiring silence of $`v`$ and its newborn until the end of the unit has conditional probability at least $`e^{-2}`$, giving stage factor $`e^{-3}M/(M+2)`$.
   - *Check:* Just before $`u`$'s ring at stopping time $`\tau \in (0, 1)`$, $`L_u \ge M+1`$ and $`c_u \le 2`$, so $`L_u/D_u \ge (M+1)/(M+3) > M/(M+2)`$; thus a mark below $`M/(M+2)`$ selects a leaf $`v`$ and attaches a leaf $`v'`$ to it, leaving $`D_u \ge M+1 > M`$. By the strong Markov property at $`\tau`$, the post-$`\tau`$ clocks of $`u, v, v'`$ are independent rate-$`1`$ Poisson processes, so the conditional silence probability of $`(v, v')`$ on $`(\tau, 1]`$ is $`e^{-2(1-\tau)} \ge e^{-2}`$. Multiplying by $`e^{-1}`$ (single ring of $`u`$) and $`M/(M+2)`$ yields the stated bound, and on this event the far side beyond $`v`$ is the single edge $`vv'`$.
   - *Verdict:* **holds**.
9. **Lines 92–97 (Stage 4, $`k-2`$ units, and no old cores):**
   - *Claim:* Silencing $`u`$ for $`k-2`$ units and adjoining the remaining $`k-2`$ vertices of $`(T, v)`$ one per unit has conditional probability at least $`(e^{-(k+1)}/(k+1))^{k-2}`$; the planted branch contains no vertex that was core at the start of the interval.
   - *Check:* In unit $`j \in \{1, \dots, k-2\}`$, the partial branch has $`j+1 \le k-1`$ vertices and is connected, so the desired recipient has an internal neighbour $`z`$ with at most $`j \le k-2`$ internal neighbours plus at most one exterior edge to $`u`$, hence $`D_z \le k-1 < k+1`$ (constant until $`z`$ rings because all neighbours of $`z`$ are silent). Requiring one ring of $`z`$ choosing the recipient (prob $`\ge e^{-1}/(k+1)`$), silence of the $`j \le k-2`$ other branch vertices ($`\ge e^{-(k-2)}`$ ), silence of $`u`$ ($`e^{-1}`$), and post-ring silence of the newborn ($`\ge e^{-1}`$) gives at least $`e^{-(k+1)}/(k+1)`$ per unit (and $`1`$ when $`k=2`$). Since $`v`$ was a leaf just before Stage 3, by degree monotonicity it was not core at time $`t`$, and all other branch vertices were born at or after Stage 3. (Note: `proofs/existence.md` line 122 uses the sharper valid denominator $`k`$ rather than $`k+1`$.)
   - *Verdict:* **holds**.
10. **Line 99 (Disjointness of planted branches):**
    - *Claim:* Branches planted at distinct old core vertices (or distinct cut edges leaving the old core subtree without containing old cores) remain disjoint.
    - *Check:* The old core set $`\mathcal{C}_t`$ induces a connected subtree $`K_t`$ at all times $`r\ge t`$. Rooting $`T(r)`$ at any $`o \in \mathcal{C}_t`$, every cut edge $`u_i w_i`$ with $`u_i \in \mathcal{C}_t`$ and $`w_i \notin \mathcal{C}_t`$ points away from $`o`$, and its far side is the subtree rooted at $`w_i`$. Two subtrees of distinct vertices $`w_1 \ne w_2`$ can intersect only if one root (say $`w_2`$) lies in the other's subtree, forcing its parent $`u_2 \in \mathcal{C}_t`$ to lie in the first far side—contradiction. (If $`w_1 = w_2`$, then $`w_1`$ has two parents $`u_1 \ne u_2`$, impossible in a tree.) Subsequent leaf additions preserve the rooted-tree orientation and disjointness.
    - *Verdict:* **holds**.

#### Section 2.2: Comparison across one cut edge (lines 101–119)
11. **Lines 102–111 and Eq. (1) (`eq:cut`, coupling across one cut edge):**
    - *Claim:* Coupling the planted branch on $`(t+h, t+h+s]`$ with a free copy of $`T`$ via two-step root acceptance $`d/(d+1)`$ and exterior silence of marks $`\le 1/M`$ at $`u`$ preserves the actual tree law and yields Eq. (1): conditional on $`\mathcal{F}_{t+h}`$ on the event of successful preparation,
      ```math
      \mathbb{E}[X \mid \mathcal{F}_{t+h}] \;\ge\; e^{-s/M}\,\mathbb{E}\!\left[C_T(s)\prod_j \frac{d_j}{d_j+1}\right].
      ```
    - *Check:* The probability computation is exact: given the full free reference history $`R`$, acceptance marks at root rings are independent $`\mathrm{Unif}[0,1]`$ (joint acceptance probability $`\prod_j d_j/(d_j+1)`$), independent of the rate-$`1/M`$ Poisson process of marks $`\le 1/M`$ at $`u`$ (silence probability $`e^{-s/M}`$). On agreement, $`D_u \ge M`$ prevents inward redirection and every reference core vertex is an actual core vertex on the far side of $`uv`$. However, the sentence "At a discrepancy continue the actual tree with its ordinary transition kernel and the reference with its own randomness" (lines 102–103) omits the brief precision that the triggering discrepant ring itself must execute its prescribed transition (redirecting to $`u`$ on a failed root acceptance, or attaching to $`v`$ on an inward mark at $`u`$) and that switching to uncoupled transitions and fresh newborn sources happens *immediately after* that triggering ring is processed (leaving the reference process untouched).
    - *Verdict:* **holds with a repair (supply proof)** — *Repair:* Let $`\rho`$ be the first time after $`t+h`$ of either (T1) a ring of $`u`$ with mark $`U \le 1/D_u`$ (ranking $`v`$ first) or (T2) a ring of $`v`$ with acceptance mark $`A > d/(d+1)`$. Process every ring in $`(t+h, \rho]`$—including the ring at $`\rho`$—by the coupled rules (so at $`\rho`$, (T1) attaches a leaf to $`v`$ and (T2) attaches a leaf to $`u`$, preserving the exact uniform conditional transition law at stopping time $`\rho`$). Immediately after processing $`\rho`$, assign every subsequent actual newborn (including the one born at $`\rho`$) an unused source outside the reference's reserved family and evaluate subsequent actual choices using an independent mark coordinate, while the reference continues on its own drivers. By strong Markov at $`\rho`$, this yields an admissible construction with the exact friend-tree law, and on agreement $`\rho > t+h+s`$.
12. **Lines 113–116 (Good reference histories and exceptional contribution):**
    - *Claim:* For $`A=8`$, $`0<\epsilon<1`$, $`b = 7/\log(1/\epsilon)`$, $`k \ge 6M/\epsilon`$, we have $`b > \epsilon`$, $`b\log(b/\epsilon)-b+\epsilon \ge 6`$, $`A\log A - A + 1 > 5`$; the exceptional probability is $`O_k(e^{-4s})`$ for $`s \ge s_0(k,M,\epsilon)`$ uniformly over $`k`$-vertex seeds $`T`$, and its contribution to $`\mathbb{E}C_T(s)`$ is $`O_k(e^{-s})`$.
    - *Check:* $`\epsilon\log(1/\epsilon) \le 1/e < 7 \implies b > \epsilon`$; $`b\log(b/\epsilon)-b+\epsilon = 7 + (b\log b - b) + \epsilon > 7 - 1 = 6`$; $`8\log 8 - 7 \approx 9.6355 > 5`$. Exponential Markov gives $`\mathbb{P}(\mathrm{Pois}(as) \ge xs) \le \exp(-s(x\log(x/a)-x+a))`$ for $`x>a`$, bounding all $`2k`$ clock violations by $`k e^{-6s} + k e^{-5s}`$. Ranking one original neighbour first at each original vertex (which exists since $`k\ge 2`$ and $`T`$ is connected), while all original degrees remain $`< M`$, each original ring with mark $`< 1/M`$ increases the sum of original degrees by $`1`$; hence avoiding degree $`\ge M`$ by $`\epsilon s`$ forces fewer than $`kM`$ marks from a $`\mathrm{Pois}(k\epsilon s/M)`$ variable with mean $`\lambda \ge 6s`$. Since $`\lambda \mapsto \mathbb{P}(\mathrm{Pois}(\lambda) \le kM)`$ is strictly decreasing in $`\lambda > 0`$, this tail is at most $`\mathbb{P}(\mathrm{Pois}(6s)\le kM) = O_{k,M}(e^{-6s}s^{kM}) \le e^{-4s}`$ for $`s \ge s_0(k,M,\epsilon)`$, independent of $`T`$. Cauchy–Schwarz with $`C_T \le N_T`$ and $`\mathbb{E}N_T(s)^2 \le k(k+1)e^{2s}`$ gives $`\mathbb{E}[C_T(s)\mathbf{1}_{\text{bad}}] \le \sqrt{k(k+1)(2k+1)}\,e^{-s}`$. Note that "good" and the root-ring weight depend only on the annotated path law, invariant under admissible rankings.
    - *Verdict:* **holds**.
13. **Line 118 (Root assignment and comparison weight lower bound):**
    - *Claim:* Assigning each good history to the least-numbered original vertex $`v`$ with $`D_v(\epsilon s)\ge M`$ gives comparison weight $`\ge e^{-\eta s}`$ with $`\eta = b\log 2 + (A+1)/M`$, and some fixed original $`v`$ accounts for at least $`(\mathbb{E}C_T(s) - O_k(e^{-s}))/k`$ of the good expectation.
    - *Check:* On $`[0, \epsilon s]`$, at most $`bs`$ rings of $`v`$, each factor $`d/(d+1) \ge 1/2`$; on $`(\epsilon s, s]`$, at most $`As`$ rings of $`v`$, each with $`d \ge M`$ so $`d/(d+1) = (1+1/d)^{-1} \ge e^{-1/d} \ge e^{-1/M}`$. Multiplying by $`e^{-s/M}`$ gives $`e^{-(b\log 2 + (A+1)/M)s} = e^{-\eta s}`$. Pigeonhole over the $`k`$ original vertices yields the stated $`v`$.
    - *Verdict:* **holds**.

#### Section 2.3: Iteration and the mean exponent (lines 120–129)
14. **Line 121 (Reduction from $`T_2`$ to $`k`$-vertex seed):**
    - *Claim:* Stopping at the first visit $`\tau_k`$ to size $`k`$ gives $`F(s) \le k + \max_{|T|=k}\mathbb{E}C_T(s)`$.
    - *Check:* On $`\{\tau_k > s\}`$, $`C(s) \le N(s) < k`$; on $`\{\tau_k \le s\}`$, $`C(s) \le C(\tau_k + s)`$ by monotonicity. Taking expectations and applying the strong Markov property at $`\tau_k`$ gives the inequality.
    - *Verdict:* **holds**.
15. **Lines 121–125 (Path functional $`X_u`$ and almost-supermultiplicative recursion):**
    - *Claim:* Defining $`X_u`$ from the actual tree path alone as the core count at $`t+h+s`$ on all far sides of edges $`uw`$ at $`t+h`$ containing no time-$`t`$ cores gives $`\sum_u X_u \le C(t+h+s)`$, and summing single-site continuations over eligible $`u`$ yields
      ```math
      F(t+h+s) \ge F(t)\,\frac{p}{4k}\,e^{-\eta s}\big(F(s) - k - O_k(e^{-s})\big).
      ```
    - *Check:* Distinct outward edges from $`\mathcal{C}_t`$ with no $`\mathcal{C}_t`$-vertex on the far side have pairwise disjoint far sides (Step 10), and each such edge has a unique endpoint in $`\mathcal{C}_t`$, so $`\sum_{u\in\mathcal{C}_t} X_u \le C(t+h+s)`$ pointwise. Because $`X_u`$ depends only on the tree path, $`\mathbb{E}[X_u \mid \mathcal{F}_t]`$ is identical across admissible continuations; evaluating it in the single-site continuation for $`u`$ with Lemma 2.2 and Eq. (1), summing over the $`\ge C(t)/4`$ eligible vertices (Lemma 2.1), and taking expectations gives the displayed bound (trivially valid at $`t=0`$ since $`F(0)=0`$).
    - *Verdict:* **holds**.
16. **Line 126 (Divergence $`F(t)\to\infty`$ and existence of $`\mu = \lim_{t\to\infty} t^{-1}\log F(t)`$):**
    - *Claim:* After the first birth, promotion intensity $`\sum_u L_u/D_u \ge 1`$, giving $`F(t)\to\infty`$; iterating $`F(t+h+s) \ge B_s F(t)`$ from $`t=1`$ gives $`\underline\mu \ge \overline\mu - \eta`$, and sending $`M\to\infty, \epsilon\downarrow 0`$ yields $`\mu \in [0,1]`$.
    - *Check:* If $`C=1`$, $`L/D=1`$; if $`C\ge 2`$, at least two core leaves ($`c=1`$) each contribute $`L/(L+1)\ge 1/2`$. Integrating the intensity (bounded by $`N(r)`$ with $`\mathbb{E}\int_0^t N(r)\,dr < \infty`$) gives $`F(t) \ge \mathbb{E}(t-\tau_1)_+ \to \infty`$. For $`s\ge s_1`$, $`F(s)-k-O_k(e^{-s})\ge F(s)/2 > 0`$, $`F(1) > 0`$, and monotonicity on $`[1+j(h+s), 1+(j+1)(h+s))`$ yields $`\liminf t^{-1}\log F(t) \ge (\log B_s)/(h+s)`$ regardless of the sign of $`\log B_s`$. Taking $`s_i\to\infty`$ attaining $`\overline\mu`$ for fixed $`(M,\epsilon,k)`$ gives $`\underline\mu \ge \overline\mu - \eta`$; sending $`\eta\to 0`$ yields $`\mu \in [0, 1]`$ (upper bound from $`F(t)\le 2e^t`$).
    - *Verdict:* **holds**.
17. **Line 128 (De-Poissonisation to discrete mean $`a_n = \mathbb{E}C_n`$):**
    - *Claim:* $`a_n \le e F(\log n)`$ and, for $`t=(1-\alpha)\log n`$ ($`0<\alpha<1`$), $`F(t) \le a_n + o(1)`$, whence $`(\log a_n)/\log n \to \mu`$.
    - *Check:* $`F(t) = \mathbb{E}[a_{N(t)}]`$ with $`a_n`$ nondecreasing. At $`t=\log n`$, $`\mathbb{P}(N(\log n)\ge n) \ge \mathbb{P}(G_1\ge n) = (1-1/n)^{n-1} \ge e^{-1}`$ for all $`n\ge 2`$, giving $`a_n \le eF(\log n)`$ and $`\limsup (\log a_n)/\log n \le \mu`$. For $`t=(1-\alpha)\log n`$, $`G_1+G_2 > n \implies \max(G_1,G_2) \ge \lfloor n/2\rfloor + 1`$, so $`\mathbb{P}(N(t)>n) \le 2(1-e^{-t})^{\lfloor n/2\rfloor} \le 2e^{e^{-t}}\exp(-n^\alpha/2)`$. Cauchy–Schwarz with $`\mathbb{E}N(t)^2 \le 6n^{2(1-\alpha)}`$ gives $`\mathbb{E}[N(t); N(t)>n] = O(n^{1-\alpha}e^{-n^\alpha/4}) = o(1)`$, so $`a_n \ge F(t) - o(1) \ge F(t)/2`$ for large $`n`$, yielding $`\liminf (\log a_n)/\log n \ge (1-\alpha)\mu`$, and letting $`\alpha\downarrow 0`$ completes the mean limit.
    - *Verdict:* **holds**.

#### Section 2.4: Independent contributions and the sample-path exponent (lines 130–161)
18. **Line 131 (Selecting at least $`C/6`$ disjoint sites):**
    - *Claim:* For $`C\ge 2`$, taking all $`A_0`$ sites with $`c\le 2, L\ge 1`$ and one end of each of $`J`$ disjoint bare triples yields at least $`C/6`$ sites; for $`C=1`$, take the centre.
    - *Check:* Core ends lie in $`A_0`$, so $`n_{\ge 3} \le A_0 - 2`$, giving at most $`2A_0 - 2 < 2A_0`$ non-bare core vertices. Suppressing maximal bare paths gives a tree on non-bare vertices with at most $`2A_0 - 3 < 2A_0 - 1`$ edges (maximal bare paths). Splitting each path into disjoint triples leaves at most $`2`$ bare vertices per path, so $`C \le 2A_0 + 3J + 2(2A_0 - 1) < 6A_0 + 3J`$, and $`A_0 + J \ge A_0 + J/2 \ge C/6`$.
    - *Verdict:* **holds**.
19. **Line 133 (Simultaneous disjoint site preparations with common probability $`p`$):**
    - *Claim:* Reserving disjoint clock sets per site, re-driving the Stage 3 promoted leaf and newborn from site-specific inputs, thinning by coins, and reverting to ordinary transitions on failures produces independent site preparation flags of common probability $`p > 0`$.
    - *Check:* The construction is sound, but two points essential for admissibility are compressed in line 133: (i) why re-driving the randomly chosen leaf $`v_i^*`$ in Stage 3 can never collide with another site's reserved vertices, and (ii) precisely how Stage 3–4 rules behave on a site that failed Stages 1–2 so that failure reversion neither breaks admissibility nor correlates flags across sites.
    - *Verdict:* **holds with a repair (supply proof)** — *Repair:* Let each site $`i`$ hold its reserved vertex clocks on $`(t, t+h+s]`$, unused tagged clocks $`\sigma^{(i)}_0,\dots,\sigma^{(i)}_k`$ and reference newborn clocks $`(\zeta^{(i)}_r)_{r\ge 1}`$, and independent uniform thinning coins. Define the four stage flags $`F_1^{(i)},\dots,F_4^{(i)}`$ purely in terms of site $`i`$'s reserved inputs (with A-site Stage 1 thinned by $`e^{-1}(1-e^{-1})/2`$, Stage 3 thinned by $`e^{-2\tau_i}`$, and Stage 4 unit $`j`$ thinned so each unit matches $`e^{-(k+1)}/(k+1)`$—or kept at the fixed seed-only probability $`\hat p_4 \ge p_4`$ common to all sites), and put $`\mathrm{Prep}_i = \bigcap_{m=1}^4 F_m^{(i)}`$. In the actual tree evolution, rule (R2) re-drives the vertex chosen by $`u_i`$'s first ring in unit 3 *only if* that vertex is a leaf of $`u_i`$ (and in particular reverts to standard rules if not, or if $`F_1^{(i)}\cap F_2^{(i)}`$ failed). Any leaf of $`u_i`$ has unique neighbour $`u_i`$, so it is neither a core vertex (hence disjoint from all triples) nor a leaf of any other A-site base $`u_j \ne u_i`$; moreover, even if $`v_i^*`$ equals $`u_i`$'s own Stage 2 tagged leaf, site $`i`$ no longer reads that source after time $`t+2`$. Thus re-driving never touches another site's reserved sources or future site-$`i`$ sources, and by disjointness of site inputs $`(\mathrm{Prep}_i)_{i=1}^m`$ are i.i.d. $`\mathrm{Bernoulli}(p)`$ (or $`\hat p \ge p`$). On $`\mathrm{Prep}_i`$, $`L_{u_i}/D_{u_i} > M/(M+2)`$ holds automatically before $`\tau_i`$, so $`v_i^*`$ is guaranteed to be a leaf of $`u_i`$ and site $`i`$ finishes Stage 4 with a copy of $`(T, v)`$ driven by $`\sigma_1^{(i)},\dots,\sigma_k^{(i)}`$ and $`D_{u_i} \ge M+1`$, regardless of exterior or failing-site arrivals.
20. **Lines 135–138 (Independent block contributions $`Y_i`$ and pointwise domination):**
    - *Claim:* Generating independent free reference histories, acceptance marks, and exterior tests per site yields i.i.d. copies of $`Y`$ with $`\sum_{i=1}^m Y_i \le C(t+h+s)`$ and $`\mathbb{E}Y \ge (p/k)e^{-\eta s}(F(s)-k-O_k(e^{-s}))`$.
    - *Check:* The wording "At each site that succeeds, generate..." is slightly informal if one wishes $`(Y_1, \dots, Y_m)`$ to be manifestly i.i.d. without conditioning on the number of successes, though mathematically equivalent to generating $`(R_i, E_i, \mathrm{Acc}_i)`$ at *every* site $`i = 1, \dots, m`$ from disjoint post-$`(t+h)`$ site randomness and setting $`Y_i = \mathbf{1}_{\mathrm{Prep}_i}\mathbf{1}_{G_v(R_i)}\mathbf{1}_{E_i}\mathbf{1}_{\mathrm{Acc}_i}C_{R_i}(s)`$. Since bases $`u_i`$ are distinct and planted branches are pairwise disjoint (Step 10), driving each prepared branch by $`R_i`$ up to its own first discrepancy is admissible simultaneously across all prepared sites and gives $`\sum_{i=1}^m Y_i \le C(t+h+s)`$ pointwise; the mean bound follows from Step 13 and Step 14.
    - *Verdict:* **holds with a repair (supply proof)** — *Repair:* Define $`(R_i, E_i, \mathrm{Acc}_i)`$ unconditionally for every $`i \in \{1,\dots,m\}`$ using site $`i`$'s reserved post-$`(t+h)`$ clocks (only coupling the actual branch of site $`i`$ to $`R_i`$ on the $`\mathcal{G}_{t+h}`$-measurable event that site $`i`$ succeeded in preparation). Then $`(\mathrm{Prep}_i, R_i, E_i, \mathrm{Acc}_i)_{i=1}^m`$ are i.i.d. and independent of $`\mathcal{F}_t`$, so $`Y_1, \dots, Y_m`$ are i.i.d. given $`\mathcal{F}_t`$ with the stated mean lower bound.
21. **Line 139 (Truncation and Chebyshev bound):**
    - *Claim:* Truncating $`Y`$ at $`R < \infty`$ so $`\mathbb{E}\min(Y, R) \ge \mathbb{E}Y/2`$ gives $`\mathbb{P}(C(t+h+s) < C\mathbb{E}Y/24 \mid \mathcal{F}_t) \le 16R^2/(m(\mathbb{E}Y)^2) \le 96R^2/(C(\mathbb{E}Y)^2)`$.
    - *Check:* $`\mathbb{E}Y \le \mathbb{E}N_T(s) = ke^s < \infty`$, so monotone convergence provides $`R < \infty`$ whenever $`\mathbb{E}Y > 0`$. With $`\bar Y_i = \min(Y_i, R)`$ and $`S = \sum_{i=1}^m \bar Y_i`$, $`\mathbb{E}S \ge m\mathbb{E}Y/2`$, $`\mathrm{Var}(S) \le mR^2`$, and $`C\mathbb{E}Y/24 \le m\mathbb{E}Y/4`$, so Chebyshev gives $`\mathbb{P}(S < m\mathbb{E}Y/4) \le mR^2/(m\mathbb{E}Y/4)^2 = 16R^2/(m(\mathbb{E}Y)^2) \le 96R^2/(C(\mathbb{E}Y)^2)`$.
    - *Verdict:* **holds**.
22. **Lines 141–145 ($`C(t)\to\infty`$ almost surely):**
    - *Claim:* Conditional promotion probability $`\ge 1/n`$ at size $`n\ge 3`$ implies $`\mathbb{P}(\text{no promotion on } [n_0, n_1)) \le (n_0-1)/(n_1-1) \to 0`$, giving $`C(t)\to\infty`$ a.s.
    - *Check:* Telescoping product $`\prod_{n=n_0}^{n_1-1}(1-1/n) = (n_0-1)/(n_1-1)`$, union over countable $`n_0 \ge 3`$, and $`N(t)\to\infty`$ a.s. (finite holding times) give $`C(t) = C_{N(t)} \to \infty`$ a.s.
    - *Verdict:* **holds**.
23. **Lines 147–158 (Block martingale, Kronecker's lemma, and a.s. lower bound):**
    - *Claim:* For a block with $`g = \mathbb{E}Y/24 > 1`$, $`q_j = \mathbb{P}(I_j = 1 \mid \mathcal{F}_{t_{j-1}}) \to 0`$ a.s., $`\sum (I_j - q_j)/j`$ is an $`L^2`$-bounded martingale converging a.s., Kronecker's lemma gives $`J^{-1}\sum_{j\le J} I_j \to 0`$ a.s., hence $`\liminf t^{-1}\log C(t) \ge (\log g)/(h+s)`$ a.s.; taking integer $`s\to\infty`$ and countable $`(M,\epsilon)`$ with $`\eta\downarrow 0`$ gives $`\liminf t^{-1}\log C(t) \ge \mu`$ a.s. (and directly for $`\mu=0`$).
    - *Check:* Every step is exact: $`|I_j-q_j| \le 1`$ gives $`\sum j^{-2} < \infty`$; $`C(t_{j-1}) \to \infty`$ a.s. implies $`q_j \to 0`$ a.s. and Cesàro mean $`\to 0`$; after the a.s.-finite first grid index $`j_0`$ with $`C(t_{j_0})\ge 1`$, monotonicity ensures failed blocks have nonnegative $`\log C`$ increment while successful blocks add $`\ge \log g > 0`$; monotonicity between grid points transfers the lower limit to continuous $`t`$. If $`\mu > 0`$, $`\eta < \mu`$ implies $`s^{-1}\log(\mathbb{E}Y/24) \ge \mu - \eta > 0`$ for large integer $`s`$; countable intersection over integer $`s`$, integer $`M\ge 2`$, and rational $`\epsilon \in (0,1)`$ preserves probability $`1`$.
    - *Verdict:* **holds**.
24. **Line 160 (Almost-sure upper bound and return to discrete time):**
    - *Claim:* Markov's inequality at integer times + Borel–Cantelli + monotonicity + countable $`\delta\downarrow 0`$ give $`\limsup t^{-1}\log C(t) \le \mu`$ a.s.; hitting times $`\tau_n`$ satisfy $`\tau_n/\log n \to 1`$ a.s., yielding Theorem 1.1.
    - *Check:* $`\mathbb{P}(C(j) > e^{(\mu+\delta)j}) \le F(j)e^{-(\mu+\delta)j} = e^{-j(\delta+o(1))}`$ is summable; $`\tau_n = \sum_{j=2}^{n-1} 1/j + \sum_{j=2}^{n-1}(e_j-1)/j = \log n + O(1)`$ a.s.; substituting $`t = \tau_n`$ into $`\log C(t)/t \to \mu`$ and using $`(C(\tau_n))_{n\ge 2} \stackrel{d}{=} (C_n)_{n\ge 2}`$ completes Theorem 1.1.
    - *Verdict:* **holds**.

---

### 2.2 Section 3: Motif bounds for the mean core size — Full Motif Drift Audit (`paper/friend-tree-exponent.tex`, lines 162–702)

#### Sections 3.1–3.2: States, transitions, and four weight tables (lines 173–258)
1. **Lines 175–195 (Core rooting, states, one-step rates, bounds):**
   - *Claim:* Rooting at the centre of $`T_3`$ keeps the root in $`K`$ and gives every non-root core vertex a core parent; every leaf has a core parent; gain at $`u`$ has probability $`r_u/n = (L_u+\rho_u)/n`$, budding at $`u`$ has probability $`b_u/n = L_u/(nD_u)`$ creating $`\nu=(1,0,1)`$ and preserving $`D_u`$; $`0\le\rho_u\le c_u/2`$ and $`L_u \le r_u \le D_u - c_u/2`$.
   - *Check:* Every attachment target is either a core vertex $`u`$ (leaf neighbours contribute $`L_u \cdot 1`$, core neighbours contribute $`\sum_{v\sim_K u} 1/D_v = \rho_u`$) or a leaf of $`u`$ (only neighbour $`u`$, probability $`1/(nD_u)`$, summed over $`L_u`$ leaves gives $`b_u/n`$). Total probability $`\sum_{u\in K}(r_u+b_u)/n = n^{-1}\sum_{v\in V} D_v/D_v = 1`$. Since core degrees are $`\ge 2`$, $`0\le \rho_u \le c_u/2`$ and $`L_u \le L_u + \rho_u \le L_u + c_u/2 = D_u - c_u/2`$.
   - *Verdict:* **holds**.
2. **Lines 199–258 (Tables, endpoint measure $`\zeta_s`$, domains $`\mathcal{P}, \mathcal{W}, H, J`$, knots $`\Xi`$):**
   - *Claim:* $`\sum_a \zeta_s(a) i(a) = 1/D`$ exactly; $`\mathcal{P}`$ and $`\mathcal{W}`$ cover all parent and child indices on core edges; $`\Phi`$ is well-defined; knots $`\Xi`$ have second knot $`1/D_0`$ and last knot $`1/2`$.
   - *Check:* For $`D\ge D_0`$, $`(1-t)\cdot 0 + t(1/D_0) = (D_0/D)(1/D_0) = 1/D`$. On any core edge $`p\to w`$, $`k_p \ge 1`$ (so $`\kappa_p \ge 1`$, support in $`\mathcal{P}`$) and $`h_w = 1`$ (support in $`\mathcal{W}`$). Note $`10/240 = 1/24 = 1/D_0`$ and $`120/240 = 1/2`$. **Expository typo:** line 236 writes "under the measures $`\pi^{\mathrm{end}}`$" instead of $`\zeta_s`$ (defined on line 212 and used in the immediately following display on lines 238–239).
   - *Verdict:* **holds** (minor notation typo at line 236).

#### Section 3.3 and Lemma 3.1 (`lem:decompose`): Exact generator decomposition (lines 260–387)
3. **Lines 262–278 (Event-wise generator and newborn list):**
   - *Claim:* $`\mathcal{L}\Phi = \sum_\tau G(\tau) + \sum_{x\in K} b_x\bigl(f(\nu)+g(x^b,\nu)+h_V(p_x,x^b,\nu)+\sum_{w\in B(x)}h_C(x^b,w,\nu)\bigr)`$, where entry into $`H`$ or $`J`$ never activates an old term and exits are covered by zero extension.
   - *Check:* Leaf gain changes only $`s_x \mapsto s_x^+`$, increasing $`D_x`$ at fixed $`(k_x, h_x)`$, so no old path or cherry can enter $`H`$ or $`J`$. Budding at $`x`$ increases $`k_x`$ by $`1`$ at fixed $`(D_x, h_x)`$ and creates $`\nu=(1,0,1)`$ with $`B(\nu)=\emptyset`$. Entry into $`H`$ ($`1\le k \le 8, D<16, h=1`$) requires budding from $`k_x=0`$ (no old core child, hence no old path through $`x`$). Entry into $`J`$ ($`2\le k \le 6, D<12`$) requires budding from $`k_x=1`$ (only one old core child, hence no old cherry at $`x`$). An existing path $`p\to u\to w`$ with $`s_u \notin H`$ already has $`h_u=1, k_u\ge 1`$, so $`s_u \notin H`$ means $`k_u>8`$ or $`D_u\ge 16`$, which persists under every event ($`G(h_V)=0`$). Likewise for an existing cherry with $`s_u \notin J`$, $`k_u \ge 2`$ implies $`k_u>6`$ or $`D_u\ge 12`$, persisting under every event ($`G(h_C)=0`$). Because $`\nu`$ has no core children, the four displayed terms exhaust all motifs containing $`\nu`$.
   - *Verdict:* **holds**.
4. **Lines 280–295 (Vertex-to-edge redistribution: $`V_u`$ in Eq. (2) and $`U_e`$):**
   - *Claim:* Splitting inner vertex gain credit $`\rho_u = h_u i_{p_u} + \sum_{w\in B(u)} i_w`$ into $`\theta_u = h_u r_{\mathrm{par}}(s_u) + k_u r_{\mathrm{ch}}(s_u)`$ kept in $`V_u`$ and transferring $`h_u(i_{p_u}-r_{\mathrm{par}}(s_u))\Delta_u f`$ to $`(p_u\to u)`$ and $`(i_w-r_{\mathrm{ch}}(s_u))\Delta_u f`$ to $`(u\to w)`$ preserves the sum; tail vertices keep $`\theta_u = \rho_u`$.
   - *Check:* Every non-root core vertex has exactly one parent edge ($`h_u = 1`$) and $`k_u = |B(u)|`$ child edges. Summing $`U_e`$ over all core edges $`e=(p\to w)`$ yields $`\sum_{p\in K}[s_p\in\mathcal{I}]\sum_{w\in B(p)}(i_w-r_{\mathrm{ch}}(s_p))\Delta_p f + \sum_{w\in K}[s_w\in\mathcal{I}]h_w(i_{p_w}-r_{\mathrm{par}}(s_w))\Delta_w f`$, which cancels $`\sum_u(\theta_u-\rho_u)\Delta_u f`$ identically. Newborn vertex $`f(\nu)`$ and edge $`g(u^b,\nu)`$ are booked once in $`V_u`$.
   - *Verdict:* **holds**.
5. **Lines 297–367 and Lemma 3.1 (`lem:decompose`, Eq. (3)–(6), edge/path/cherry redistribution):**
   - *Claim:* On each edge $`e=(p\to w)`$, splitting every gain credit $`i_z \Delta_x g(e)`$ ($`x\in\{p,w\}`$ inner, $`z\ne y`$ another core neighbour of $`x`$) by the table on lines 303–309 and booking newborn cherries $`b_p h_C(p^b,w,\nu)`$ and newborn paths $`b_w h_V(p,w^b,\nu)`$ on $`e`$ yields Eq. (6):
     ```math
     \mathcal{L}\Phi = \sum_u V_u + \sum_e E_e + \sum_{p\to u\to w:\, s_u\in H} P_{puw} + \sum_{u:\, s_u\in J}\sum_{\{w_1,w_2\}\subset B(u)} B_{u12}.
     ```
   - *Check (exhaustive orientation-by-orientation bookkeeping):*
     - *Tail endpoint $`x\in\{p,w\}`$:* Keeps full rate $`\sigma_x = r_x`$, transfers nothing, contributes $`0`$ to $`n_1,n_2,n_3`$ and $`R_e`$.
     - *Inner $`p`$, neighbour $`w`$:* Kept in $`\sigma_p`$.
     - *Inner $`p`$, parent $`z=p_p`$ ($`h_p=1`$):* If $`s_p\in H`$ (note $`s_p\in H \implies h_p=1`$), $`x_{\mathrm{par}}(s_p)d_p`$ is kept in $`\sigma_p d_p`$ and $`(i_{p_p}-x_{\mathrm{par}}(s_p))\Delta_p g(p,w)`$ goes to the unique path $`p_p\to p\to w`$ (first source term of Eq. (4), `eq:motif-path`). If $`s_p\in\mathcal{I}\setminus H`$ ($`n_1=1`$), $`i_{p_p}d_p - A_{pc}^p(i_{p_p})`$ stays in $`R_e`$ (row 1, line 342) and $`A_{pc}^p(i_{p_p})`$ moves to parent edge $`e'=(p_p\to p)`$; conversely $`e'=(p'\to w')`$ receives $`A_{pc}^{w'}(i_{p'})`$ once from each of its $`k_{w'}`$ child edges, matching $`n_3' A_{pc}^{w'}(i_{p'})`$ in $`I_{e'}`$ (line 335).
     - *Inner $`p`$, sibling $`z\in B(p)\setminus\{w\}`$ ($`k_p - 1`$ siblings):* If $`s_p\in J`$, $`(k_p-1)x_{cc}(s_p)d_p`$ is kept in $`\sigma_p d_p`$ and $`(i_z - x_{cc}(s_p))\Delta_p g(p,w)`$ goes to the unique unordered cherry $`(p;\{w,z\})`$ (Eq. (5), `eq:motif-cherry`, receiving one term from $`p\to w`$ and one from $`p\to z`$). If $`s_p\in\mathcal{I}\setminus J`$ ($`n_2 = k_p-1`$), $`i_z d_p - A_{cc}^p(i_z)`$ stays in $`R_e`$ (row 2, line 343) and $`A_{cc}^p(i_z)`$ moves to sibling edge $`e'=(p\to z)`$; conversely $`e'=(p\to z)`$ receives $`A_{cc}^p(i_z)`$ from each of its $`k_p - 1`$ sibling edges, matching $`n_2 A_{cc}^p(i_z)`$ in $`I_{e'}`$ (line 335).
     - *Inner $`w`$, parent $`p`$:* Kept in $`\sigma_w`$.
     - *Inner $`w`$, child $`z\in B(w)`$ ($`k_w`$ children):* If $`s_w\in H`$, $`k_w x_{\mathrm{ch}}(s_w)d_w`$ is kept in $`\sigma_w d_w`$ and $`(i_z - x_{\mathrm{ch}}(s_w))\Delta_w g(p,w)`$ goes to the unique path $`p\to w\to z`$ (second source term of Eq. (4), `eq:motif-path`). If $`s_w\in\mathcal{I}\setminus H`$ ($`n_3 = k_w`$), $`i_z d_w - A_{cp}^w(i_z)`$ stays in $`R_e`$ (row 3, line 344) and $`A_{cp}^w(i_z)`$ moves to child edge $`e'=(w\to z)`$; conversely $`e'=(p'\to w')`$ receives $`A_{cp}^{p'}(i_{w'})`$ from its unique parent edge iff $`h_{p'}=1`$ and $`s_{p'}\in\mathcal{I}\setminus H`$, matching $`n_1' A_{cp}^{p'}(i_{w'})`$ in $`I_{e'}`$ (line 335).
     - *Newborn paths and cherries:* Every newborn path $`b_x h_x h_V(p_x, x^b, \nu)`$ corresponds bijectively to the edge $`p_x\to x`$ (booked on $`w=x`$ in Eq. (3)); every newborn cherry $`b_x h_C(x^b, w, \nu)`$ corresponds bijectively to the edge $`x\to w`$ (booked on $`p=x`$ in Eq. (3)).
   - *Verdict:* **holds** (independently verified symbolically and numerically in exact `Fraction` arithmetic in Section 3.3 below).

#### Section 3.4 and Lemma 3.2 (`lem:motif-corners`): Reduction to finite inequalities (lines 389–587)
6. **Lines 397–408 and 439–441 (Lemma 3.2(1), multi-affine interpolation and rate intervals):**
   - *Claim:* Separate affinity gives $`d = \sum_c \omega_c d(c)`$; for $`r\in I=[l,h]`$, $`rd \preceq \sum_c \omega_c M_I(d(c))`$.
   - *Check:* Immediate from $`\omega_c \ge 0`$ and $`rd(c) \preceq M_I(d(c))`$ for either sign of $`d(c)`$.
   - *Verdict:* **holds**.
7. **Lines 409–414 and 441–446 (Lemma 3.2(2), Fenchel–Young remainders and received knot functions):**
   - *Claim:* For $`z\in[0,1/2]`$, $`zd - A(z) \preceq A^*(d) = M_{\xi\in\Xi}\{\xi d - A(\xi)\} \preceq \sum_c \omega_c A^*(d(c))`$; at a tail endpoint $`v`$, $`A(i_v) = (1-t)A(0) + tA(1/D_0)`$.
   - *Check:* On each knot interval $`[\xi_j, \xi_{j+1}]`$, $`zd - A(z)`$ is affine in $`z`$, hence bounded by its endpoint values; $`A^*`$ is convex (UPPER, $`M=\max`$) or concave (LOWER, $`M=\min`$), so Jensen's inequality gives $`A^*(\sum_c \omega_c d(c)) \preceq \sum_c \omega_c A^*(d(c))`$. Because $`[0, 1/D_0]`$ is the first knot interval of $`\Xi`$, $`A(t/D_0)`$ is exact linear interpolation in $`t`$.
   - *Verdict:* **holds**.
8. **Lines 415–426 and 448–459 (Lemma 3.2(3), tail endpoint gain and budding):**
   - *Claim:* At tail endpoint $`v`$ with $`q=q_{\kappa h}`$, $`r_v\Delta_v F`$ has corner bounds $`0`$ at $`\infty`$ and $`M_{[\gamma_0,1]}(F[\infty]-F[q])`$ at $`q`$; $`b_v\nabla_v F`$ has $`0`$ at $`\infty`$ and $`M_{[\beta_0,1]}(F[q_{\kappa+1,h}]-F[q])`$ at $`q`$ (for $`\kappa \le K_T`$, else $`0`$), with $`\gamma_0 = (D_0-\kappa-h)/(D_0+1)`$, $`\beta_0 = (D_0-\kappa-h)/D_0`$ for $`\kappa \le K_T`$ and $`\gamma_0=0`$ for $`\kappa = K_T+1`$.
   - *Check:* Under gain, $`t^+ - t = D_0/(D+1) - D_0/D = -t/(D+1)`$, so $`r_v \Delta_v F = t \cdot \frac{r_v}{D+1}(F[\infty]-F[q])`$. Since $`r_v \le D < D+1`$, $`r_v/(D+1) \le 1`$; for $`\kappa \le K_T`$, $`k_v = \kappa`$, so $`r_v/(D+1) \ge (D-\kappa-h)/(D+1) = 1 - (\kappa+h+1)/(D+1)`$, increasing in $`D \ge D_0`$ (and positive since $`D_0 - K_T - 1 = 11 > 0`$). For $`\kappa = K_T+1`$, $`r_v \ge 0`$. Under budding, $`t`$ is unchanged; for $`\kappa = K_T+1`$, $`\kappa`$ is unchanged so $`\nabla_v F = 0`$; for $`\kappa \le K_T`$, $`b_v = 1 - (\kappa+h)/D \in [\beta_0, 1]`$ (increasing in $`D\ge D_0`$). Interpolating any second tail endpoint before applying $`M_I`$ gives the stated corner bounds.
   - *Verdict:* **holds**.
9. **Lines 427–435 and 461–470, Eq. (7) (`eq:motif-quadratic`, Lemma 3.2(4), quadratic tail credit):**
   - *Claim:* For tail $`v`$ contributing $`i_v`$ to the gain rate of another vertex $`x`$ in motif $`F`$, $`i_v \Delta_x F`$ has corner bounds $`0`$ at $`v=\infty`$ and $`D_0^{-1}M\{d_0, d_1\}`$ at $`v=q`$, including when a second tail endpoint is present.
   - *Check:* $`i_v \Delta_x F = (t/D_0)((1-t)d_0 + td_1) \preceq (t/D_0) M\{d_0, d_1\}`$. With a second tail parameter $`t'`$, $`d_0(t'), d_1(t')`$ are affine in $`t'`$, so $`M\{d_0(t'), d_1(t')\}`$ is convex (UPPER) or concave (LOWER) in $`t'`$, hence bounded by its $`t'`$-linear interpolation; multiplying by $`t/D_0 \ge 0`$ assigns $`0`$ to $`(\infty, c')`$ and $`D_0^{-1}M\{d_0(c'), d_1(c')\}`$ to $`(q, c')`$. When two tail endpoints both contribute to the same centre's gain rate, applying the rule to each summand separately and adding is valid.
   - *Verdict:* **holds**.
10. **Lines 473–503 (Vertex inequalities, Cases 1–5):**
    - *Claim:* The five cases bound Eq. (2) by $`\lambda f(s)`$ for all core states $`(L,k,h)`$ with $`D\ge 2`$.
    - *Check:* Case 1 ($`D < D_0`$, inner): exact rational evaluation; $`u^b`$ has same degree $`D < D_0`$ and $`\nu`$ is inner. Case 2 ($`D\ge D_0, k\le K_F, L < L_F`$): $`\rho_u \in [0, c/2]`$ covered by the two rate endpoints; $`B_s`$ evaluated at exact $`t = D_0/D`$. Case 3 ($`k\le K_F, L=L_F`$): $`D \ge L_F = 800 \ge D_0`$, $`\Delta_u f = 0`$, exact $`B_s`$. Case 4 ($`k\le K_F, L\ge L_F+1`$): $`\Delta_u f = 0`$, post-bud state projects to $`\pi(L_F, k+1, h)`$; $`b = L/(L+k+h) \in [(L_F+1)/D_*, 1]`$ and $`t = D_0/D \in (0, D_0/D_*]`$ with $`D_* = L_F+1+k+h`$; $`B_s`$ is separately affine in $`(b, t)`$, so checking the four rectangle corners covers all $`(b, t)`$. Case 5 ($`k > K_F`$, $`\star_h`$): both vertex jumps vanish, $`D \ge K_F+1+h \ge 151 > D_0`$, $`\kappa(u^b) = K_T+1`$, separately affine on $`[0,1]\times[0, D_0/(K_F+1+h)]`$. All five cases match `verify/motif-bracket/cert_verify_gmp.c` lines 396–495 line for line.
    - *Verdict:* **holds**.
11. **Lines 505–538 (Edge inequalities):**
    - *Claim:* Checking $`g[\text{tail},\text{tail}]=0`$ and the listed corner sums on the remaining pairs in $`\mathcal{P}\times\mathcal{W}`$ implies $`E_e \preceq \lambda g(e)`$ for every core edge.
    - *Check:* Tail–tail pairs have $`g\equiv 0`$, zero jumps (degrees stay $`\ge D_0`$), zero $`U_e, I_e, R_e`$, and cannot enter $`H`$ or $`J`$ (since $`D_0 = 24 > 16, 12`$), so $`E_e = 0 = \lambda g(e)`$. For pairs with at least one inner endpoint, every summand matches Lemma 3.2(1)–(4): $`U_e, I_e`$ and inner budding terms (including newborn $`h_C, h_V`$) are affine in any tail endpoint; gain from $`D_0 - 1`$ lands at $`q_{\kappa h}`$ ($`t=1`$)—and if the opposite end is tail, $`g[q_{\kappa h}, \text{tail}]=0`$ matches the tail–tail block; $`R_e`$ is bounded via Lemma 3.2(2); inner kept gain omits any tail neighbour credit, which is supplied by Eq. (7) alongside the tail endpoint's Lemma 3.2(3) terms. Matches `cert_verify_gmp.c` lines 341–384 and 506–689 exactly.
    - *Verdict:* **holds**.
12. **Lines 540–586 (Path and cherry inequalities, and global drift conclusion):**
    - *Claim:* Checking the listed corner bounds on $`\mathcal{P}\times H\times\mathcal{W}`$ and unordered pairs in $`J\times\mathcal{W}^2`$ (with repetition) implies $`P_{puw}\preceq\lambda h_V`$ and $`B_{u12}\preceq\lambda h_C`$, hence $`\mathcal{L}\Phi \preceq \lambda\Phi`$.
    - *Check:* In both Eq. (4) and Eq. (5), each edge-source product pairs the degree credit of one endpoint with an edge jump depending only on the centre and the *other* endpoint, so it is separately affine even when both endpoints are tail. All five outside-core-neighbour counts $`N`$ in the table (lines 555–563) are exact: top of path $`h_p+k_p-1\ge 0`$, bottom of path $`k_w`$, middle of path $`k_u-1\ge 0`$ (since $`h_u=1`$ and parent $`p`$ + child $`w`$ are inside the motif), child of cherry $`k_{w_j}`$, apex of cherry $`h_u+k_u-2\ge 0`$ (valid for root $`h_u=0`$ and non-root $`h_u=1`$). Tail endpoints apply Lemma 3.2(3) and Eq. (7) independently (including both endpoints simultaneously and diagonal $`b_1=b_2`$ in cherries, since two distinct children can share the same table index). Exits from $`H`$ or $`J`$ evaluate to $`0`$. Matches `cert_verify_gmp.c` lines 692–960 exactly.
    - *Verdict:* **holds**.

#### Section 3.5, Lemma 3.3 (`lem:motif-certificate`), and Motif Proof of Theorem 1.2 (lines 588–701)
13. **Lines 590–618 (Per-vertex comparison with $`|K|`$):**
    - *Claim:* Grouping each edge and path to its bottom vertex and each cherry to its apex yields $`\Phi = \sum_{u\in K} W_u`$, bounded pointwise by a finite enumeration over $`0\le L\le L_F, 0\le k\le K_F, h\in\{0,1\}`$ ($`D\ge 2`$) and $`\star_0, \star_1`$.
    - *Check:* Every non-root vertex ($`h_u=1`$) has one incoming edge and at most one incoming path (covered by including $`0`$ in $`g_m, v_m`$), while the root ($`h_u=0`$) has neither; an apex $`u\in J`$ has precisely $`\binom{k_u}{2}`$ cherries (and $`0`$ if $`u\notin J`$). Taking $`S_u = \{\infty, q_{\kappa h}\}`$ for tail states (even when $`t=1`$ at $`D=D_0`$) is conservative and valid. Increasing $`L > L_F`$ changes neither $`\pi(s_u)`$ nor $`S_u`$ (since $`L_F \ge D_0 > 12`$), and $`k > K_F`$ maps to $`\star_h`$ with $`D \ge K_F+1+h > D_0 > 12`$ and $`\kappa = K_T+1`$. Matches `cert_verify_gmp.c` lines 963–1017.
    - *Verdict:* **holds**.
14. **Lines 620–671 (Lemma 3.3, finite verification):**
    - *Claim:* Exact arithmetic evaluation in `cert_verify_gmp.c` certifies the stated drift inequalities and per-vertex comparison bounds for $`\Phi_+`$ and $`\Phi_-`$; denominator $`Q\cdot 10^{10}\cdot 2D_1 D_u D_2(D_0+1)\cdot 10^6`$ clears every path/cherry denominator.
    - *Check:* All counts (`572` inner, `601` endpoint, `555` parent, `290` child, `84` path centres, `75` cherry centres, `241,904` vertex entries, `484,138 / 160,545 / 13,519,800 / 3,164,625` checks) and denominator clearing (noting $`D_1=D_0`$ or $`D_2=D_0`$ whenever a tail $`q`$ occurs, so $`D_0 \mid 2D_1D_uD_2(D_0+1)`$) are exact. Independently confirmed by compiling `cert_verify_gmp.c` with UBSan and re-running both certificates (Section 3.2). Note as a software-hardening observation (documented in `proofs/bracket-motifs-code-map.md` §3) that `cert_verify_gmp.c` does not programmatically assert `D0 <= KF + 1` on arbitrary headers, though the shipped headers satisfy $`24 \le 151`$ and print `Domains: f(800,150) D0=24(...)`.
    - *Verdict:* **holds**.
15. **Lines 673–701 (Proof of motif bounds in Theorem 1.2 and a.s. upper bound):**
    - *Claim:* Iterating $`\mathbb{E}\Phi_\pm(T_{j+1}) \preceq (1 + \lambda_\pm/j)\mathbb{E}\Phi_\pm(T_j)`$ from $`j=3`$ (where $`\Phi_\pm(T_3) = f_\pm(2,0,0)`$) gives the first pair of product bounds in Theorem 1.2; dyadic Markov + Borel–Cantelli + monotonicity yield $`\limsup (\log C_n)/\log n \le 1133/2000`$ a.s.
    - *Check:* Taking expectations of the pointwise linear inequality does not require $`\Phi \ge 0`$; $`1 + \lambda_\pm/j > 0`$ preserves direction; $`c_0 > 0, C_0 > 0`$ give $`\mathbb{E}C_n \le \mathbb{E}\Phi_+/c_0`$ and $`\mathbb{E}C_n \ge \mathbb{E}\Phi_-/C_0`$. Exact ratio reduction: $`14608704616 / 17527873930 = 7304352308 / 8763936965`$ ($`\gcd = 2`$), while $`\gcd(14606753031, 2989767809) = 1`$. For $`\alpha = 1133/2000`$ and $`\varepsilon > 0`$, $`\mathbb{P}(C_{2^j} > 2^{j(\alpha+\varepsilon)}) = O(2^{-j\varepsilon})`$ is summable; monotonicity on $`[2^j, 2^{j+1})`$ and intersecting over $`\varepsilon = 1/m`$ gives the a.s. upper limit.
    - *Verdict:* **holds**.

---

### 2.3 Section 4: An independent upper bound from degree bands (`paper/friend-tree-exponent.tex`, lines 703–1072)

#### Section 4.1 and Lemma 4.1 (`lem:banded-history`): States along a growth history (lines 715–813)
1. **Lines 717–770 (Rates, band table, individual/witness/pooled alphabets, tie-breaking):**
   - *Claim:* Gain rate $`r_u = L + \sum_{w\sim_K u} 1/D_w`$ and conversion rate $`b_u = L/D`$ partition all transitions with rates summing to $`n`$; band cutoffs $`h_c(D)`$ and core-degree caps $`K_1=20, K_2=7, K_{c\ge 3}=5`$ define finite multisets with at most one pooled mark $`> h_c`$.
   - *Check:* Every band cutoff $`h_c(D)`$ in lines 741–747 matches `verify/banded-upper/geometry.txt` and is nonincreasing in both $`c`$ and $`D`$. All 7 band state counts match independent combinatorial computation (Section 3.1).
   - *Verdict:* **holds**.
2. **Lines 772–813 (Lemma 4.1, history invariant):**
   - *Claim:* Along every realised history, recorded low degrees and unsaturated core degrees are exact, saturated core degrees and pooled marks are valid lower bounds, and every required witness is an identifiable distance-$`2`$ vertex with the recorded data.
   - *Check (all transition classes):*
     - *Base ($`n=3`$):* Centre has $`(2,0)`$ and empty multiset.
     - *Central gain ($`D\mapsto D+1`$):* Within a band, multiset unchanged ($`c=1`$ alphabet and witness condition are band-independent because $`h_1 = 620`$ is constant). Across a band boundary, $`h_c'`$ can only decrease, so old exact degrees either stay exact or encode as a valid lower-bound mark, old pooled marks $`b \ge h_c \ge h_c'`$ project to valid lower-bound marks $`\le b`$, and resetting non-largest marks to $`h_c'`$ preserves lower bounds. At $`D+1=H`$, multiset dropped.
     - *Central conversion ($`c\mapsto c+1, L>0`$):* Both $`h_c`$ and $`K_c`$ can only decrease, so old labels project validly; fresh newborn $`(2,1)`$ is inserted (encoded as pooled mark $`2`$ when $`h_{c+1}=2`$, or individual $`(2,1)`$ otherwise—never requiring a witness since $`k=1`$).
     - *Observed neighbour gain/conversion and exact witness gain/conversion:* Gain raises $`d\mapsto d+1`$ (entering first pool mark $`h_c`$ at $`d+1=h_c`$, preserving witness while $`d+1\le 14`$); conversion at $`k<K_c, k<d`$ raises $`k\mapsto k+1`$ and, if $`c=1, d\le 14, 2\le k+1\le 4`$, selects the newly promoted $`(2,1)`$ vertex (guaranteed distinct from $`u`$ because it was a leaf of $`w`$) as witness; saturated conversion leaves $`(d, K_c)`$ fixed; exact witness $`(e,j)`$ increments $`e`$ (pooling at $`6`$) or $`j`$; events at pooled neighbours/witnesses leave lower bounds valid by degree monotonicity.
     - *Newborn initialisation (lines 801–810):* Newborn has $`(2,1)`$ and sees creator with $`(D_u, c_u+1)`$. A witness is required iff $`D_u\le 14`$ and $`1\le c_u\le 3`$. In band $`[2,16)`$, $`h_1=620, h_2=240, h_3=45`$ are all $`\ge 6`$ (in fact $`h_1,h_2,h_3 \ge 7`$ across all bands) and $`K_1=20, K_2=7, K_3=5`$ are all $`\ge 5`$. Hence the first entry of the creator's *old* sorted multiset (breaking ties by vertex birth order) is a core neighbour of $`u`$ distinct from the newborn, gives its exact pair $`(e,j)`$ whenever $`e<6`$ (even at $`(e,j)=(5,5)`$ with $`K_3=5`$, since $`j\le e = 5`$), and certifies $`e\ge 6`$ otherwise. Tail creators ($`D_u\ge H`$) give pooled mark $`H`$.
   - *Verdict:* **holds**.

#### Section 4.2: The generator and the two transfers (lines 815–867)
3. **Lines 817–835 (Generator identity):**
   - *Claim:* With $`I_u = \{u\}\cup S_u\cup J_u`$ (disjoint by the tree property; $`S_u=J_u=\emptyset`$ in the tail), $`n\mathbb{E}[\Phi_{n+1}-\Phi_n\mid\mathcal{F}_n] = \sum_{u\in K}[\sum_{x\in I_u}(r_x\Delta_x^+ z_u + b_x\Delta_x^b z_u) + b_u z_{\mathrm{new}}(u)]`$.
   - *Check:* Every leaf has a unique core neighbour, and events outside $`I_u`$ leave $`z_u`$ unchanged by Lemma 4.1; interchanging finite sums over events and old vertices counts each old change and each newborn exactly once.
   - *Verdict:* **holds**.
4. **Lines 837–849 (Signed edge transfer $`g_{uw}`$ and mutual observability):**
   - *Claim:* Requiring $`D_u, D_w < H`$, $`D_w < h_{c_u}(D_u)`$, $`c_w < K_{c_u}`$, $`D_u < h_{c_w}(D_w)`$, $`c_u < K_{c_w}`$ for nonzero $`g_{uw}`$ ensures both endpoints read $`g_{uw}`$ on every history and yields exact cancellation between $`+r_u g_{wu}`$ in row $`u`$ and $`-r_w g_{uw}`$ in row $`w`$.
   - *Check:* Because $`D_u, c_u, D_w, c_w`$ never decrease and cutoffs/caps never increase along any history, if $`w`$ had ever been pooled or saturated at $`u`$ earlier, it would still satisfy $`D_w \ge h_{c_u}(D_u)`$ or $`c_w \ge K_{c_u}`$ now. Contrapositively, strict inequality guarantees exact, unsaturated observation at both ends on every realised history, independently of witnesses. Summing $`\sum_u r_u \sum_{w\sim_K u} g_{wu} - \sum_u \sum_{w\in S_u} r_w g_{uw}`$ cancels identically upon relabelling $`u\leftrightarrow w`$. Matches `verify/banded-upper/check.cpp` lines 196–201, 277–283, 288, 311.
   - *Verdict:* **holds**.
5. **Lines 851–866 (Nonnegative arrow transfer $`a(u)`$ and amended row $`R_u`$):**
   - *Claim:* Subtracting $`a(u)/D_w`$ from row $`u`$ and adding it to row $`w`$ for every ordered core edge with $`D_w\ge T=8`$ cancels identically (including tail senders) and gives $`\sum_u R_u = n\mathbb{E}[\Phi_{n+1}-\Phi_n\mid\mathcal{F}_n]`$.
   - *Check:* $`\sum_u \sum_{w\sim_K u} \mathbf{1}_{\{D_w\ge T\}} a(u)/D_w = \sum_w (\mathbf{1}_{\{D_w\ge T\}}/D_w)\sum_{u\sim_K w} a(u)`$. Expanding $`r_u U_u`$ into leaves and core neighbours produces the displayed formula for $`R_u`$.
   - *Verdict:* **holds**.

#### Section 4.3: Finite rows and the boundary condition (lines 868–965)
6. **Lines 870–895 (Term upper bound $`Q_z`$ and central events):**
   - *Claim:* $`Q_z(v, r_-, r_+) = \lceil\max\{r_-(v-z), r_+(v-z)\}\rceil`$ upper-bounds $`r(v-z)`$ for $`r\in[r_-, r_+]`$; grouping identical individual and pooled labels gives valid upper bounds for central gain (with arrow subtraction) and central conversion (using trial reward $`s(t)`$).
   - *Check:* For pooled mark $`b`$, true degree $`D_w \ge b`$ gives $`1/D_w \in [0, 1/b]`$; if $`b\ge T`$ then $`D_w\ge T`$ certainly holds, and if $`b<T`$ omitting $`-a(u)\mathbf{1}_{\{D_w\ge T\}} \le 0`$ only increases the upper bound (since $`a(u)\ge 0`$ and $`1/D_w > 0`$). For `check.cpp` line 30 (`up(n,d) = n/d + (n%d > 0)`), C++17 truncation toward zero makes `up(n,d)` exact $`\lceil n/d \rceil`$ for both $`n\ge 0`$ and $`n<0`$.
   - *Verdict:* **holds**.
7. **Lines 897–925 (Neighbour and witness rate intervals):**
   - *Claim:* The four bulleted rate intervals cover every exact/saturated near gain, witnessed near gain, near conversion, and exact witness gain/conversion; multiplying interval endpoints by group multiplicity before $`Q_z`$ is valid.
   - *Check:* Every representative of an identical-label group shares the same transition target and transfer $`g_{uw}`$, so summing rates across the group yields $`[pr_-, pr_+]`$.
     - Unwitnessed gain: $`r_w = d - k' + 1/D + \sum_{k'-1} 1/D_y \in [d-k'+1/D,\ d-(k'+1)/2+1/D]`$. At $`k<K_c`$, $`k'=k`$ is exact. At saturation $`k=K_c`$, $`k \le k' \le d`$ gives lower bound $`1/D`$ and upper bound decreasing in $`k'`$, maximized at $`k'=k`$.
     - Witnessed gain ($`k`$ exact since $`k\le 4 < K_1=20`$): replaces one of the $`k-1`$ other core credits by $`1/e`$ (exact) or $`[0, 1/6]`$ (pooled $`e\ge 6`$), leaving $`k-2`$ credits in $`[0, 1/2]`$.
     - Near conversion: $`(d-k)/d`$ when $`k<K_c, k<d`$; $`0`$ at $`k=d`$; zero weight change at saturation.
     - Exact witness $`(e,j)`$: gain rate in $`[e-j+1/d,\ e-j+1/d+(j-1)/2]`$, conversion rate $`(e-j)/e`$ when $`j<e`$; pooled witness has zero change.
     All four formulas match `check.cpp` lines 300–322 and `verify_banded.cpp` lines 881–908 identically.
   - *Verdict:* **holds**.
8. **Lines 926–940 (Sender charge $`F_u`$ and finite row inequality Eq. (8), `eq:banded-row`):**
   - *Claim:* $`F_u = \mathbf{1}_{\{D\ge T\}}\lceil D^{-1}\sum_t p_t a_{\max}(t)\rceil`$ upper-bounds the outgoing arrow sum, so Eq. (8) implies $`R_u \le \lambda z_u`$ once newborn consistency holds.
   - *Check:* Every true neighbour type consistent with label $`t`$ satisfies $`a(w) \le a_{\max}(t)`$; summing ceilings of upper bounds gives $`R_u \le F_u + \sum Q_z \le \lfloor \lambda z\rfloor \le \lambda z`$.
   - *Verdict:* **holds**.
9. **Lines 942–964 (Tail region $`D \ge H`$ and boundary rounding at $`D = H-1`$):**
   - *Claim:* For $`z = A + Bq`$ ($`q = L/D`$, $`A>0, B\ge 0`$), $`R_u \le (B+a_*)(1-q)+sq`$; checking $`B+a_* \le \lambda A`$ and $`s \le \lambda(A+B)`$ suffices for all $`D\ge H`$ (margins $`6235`$ and $`4755`$ after multiplying by $`9983`$); rounding the boundary target at $`D=H-1`$ up to $`A + \lceil B(L+1)/H\rceil`$ is valid.
   - *Check:* Gain changes $`q`$ by $`(L+1)/(D+1) - L/D = (1-q)/(D+1)`$; conversion (rate $`b_u = q`$) changes $`q`$ by $`-1/D`$, contributing $`-q/D`$. Since $`r_u \le L + c/2 \le D < D+1`$ and $`q\in[0,1]`$, $`r_u(1-q)/(D+1) - q/D \le 1-q`$. Signed transfers and $`a(u)`$ vanish; outgoing arrow charge is $`\le a_* c/D = a_*(1-q)`$; newborn contribution is $`\le q s(H)`$. Both sides are affine in $`q\in[0,1]`$, so checking $`q=0, 1`$ suffices. Exact integer margins verified in Section 3.1: $`5657\cdot 1380739398522 - 9983(383976882193+398437500000) = 6235 > 0`$ and $`5657(1380739398522+383976882193) - 9983\cdot 10^{12} = 4755 > 0`$. At $`D=H-1`$, every rate multiplying $`z^+`$ is nonnegative, so rounding $`z^+`$ up preserves the inequality direction.
   - *Verdict:* **holds**.

#### Section 4.4, Lemma 4.2 (`lem:banded-certificate`), and Proposition 4.3 (lines 966–1072)
10. **Lines 978–1018 (Lemma 4.2, DAG ordering, residual monotonicity, closing recursion):**
    - *Claim:* Reverse order in $`(D, c, \text{lex multiset})`$ computes every non-newborn target before use; the integer residual is nondecreasing in $`z`$; checking $`s(t) \ge z(2,1;t)`$ for all $`c=1`$ labels (including $`t=H`$) closes the recursion.
    - *Check:* Central gains raise $`D`$; central conversions raise $`c`$ at fixed $`D`$; near gains, near conversions, and witness events all strictly advance one label (including witness replacement by $`(2,1)`$ on near conversion because $`k`$ increases first in the lexicographic order, and pool entry creates only the first pool mark $`h_c`$, never demoting an existing mark $`> h_c`$), hence strictly advance the sorted multiset. Increasing $`z`$ increases $`\lfloor\lambda z\rfloor`$ and decreases every $`v-z`$, and since all rate endpoints are nonnegative, each $`-Q_z`$ is nondecreasing in $`z`$. Because every newborn rate $`b_u \ge 0`$ (and tail $`q \ge 0`$), $`z(2,1;t) \le s(t)`$ lets actual newborn weights replace trial rewards everywhere.
    - *Verdict:* **holds**.
11. **Lines 1020–1038 (Dual verifier outputs):**
    - *Claim:* `check.cpp` and independent `verify_banded.cpp` verify all $`226{,}545{,}617`$ rows with the reported slacks, newborn reserves, weight ranges, and initial weights.
    - *Check:* Independently compiled both with `-fsanitize=undefined -fno-sanitize-recover=all` and re-executed on `geometry.txt` and `recipe.txt.gz`; every digit reproduced identically (Section 3.2).
    - *Verdict:* **holds**.
12. **Lines 1040–1071 (Proposition 4.3 and a.s. upper bound):**
    - *Claim:* Summing $`R_u \le \lambda z_u`$ gives $`\mathbb{E}[\Phi_{n+1}\mid\mathcal{F}_n] \le (1+\lambda/n)\Phi_n`$, yielding the second upper product bound in Theorem 1.2 and $`\limsup (\log C_n)/\log n \le 5657/9983`$ a.s.
    - *Check:* Tower property + iteration from $`\Phi_3 = z_{(2,0)}`$ + pointwise lower bound $`m C_n \le \Phi_n`$ ($`m = \min(\min_{\text{fin}} z, A) > 0`$) give the product bound for all $`n \ge 3`$. Dyadic Markov + Borel–Cantelli + monotonicity give the a.s. upper limit.
    - *Verdict:* **holds**.

---

### 2.4 Section 6: Related bounds and degree counts (`paper/friend-tree-exponent.tex`, lines 1128–1190)

Because §6 states three results whose main proofs reside in the supplement plus one inline deduction, every statement is checked both for faithfulness to the supplement and for step-by-step soundness of the underlying proof.

1. **Lines 1135–1149 (Lower drift on individual paths; `proofs/almost-sure-lower.md` Theorems A–B and `proofs/cross-charging.md` Lemmas 2–3, §§6–7):**
   - *Claim:* Under (H1)–(H5), $`\liminf (\log C_n)/\log n \ge \delta`$ a.s.; applies to the $`543/1000`$ baseline table and the $`11/20`$ cross-charged certificate (independently of Theorem 1.1).
   - *Check:* In `proofs/almost-sure-lower.md`, Lemma 1 ($`C_n\to\infty`$ a.s. by pigeonhole + conditional Borel–Cantelli), Lemma 2 ($`\sup\Phi_n=\infty`$ a.s. via bounded-increment martingale dichotomy and positive newborn weight $`\eta`$), Lemma 3 (reciprocal drift $`\mathbb{E}_n[1/\Phi_{n+1}] \le (1/\Phi_n)(1-(\delta-\varepsilon)/n)`$ on $`\{\Phi_n \ge K_0\}`$ using $`\Delta^2 \le B(\Delta+2\Delta^-)`$ and union bound over $`\le \Phi_n/f_{\min}`$ positive-weight vertices), and Theorem B (stopped supermartingale $`G_n = 1/(\Phi_{n\wedge\tau}p_{n\wedge\tau})`$ from crossing time $`\sigma`$, maximal inequality $`\mathbb{P}(\tau<\infty)\le K_0/K_1`$) all hold step by step. In `proofs/cross-charging.md`, the two exact cancellations ($`\kappa`$ and $`\nu`$), interval/tail/beyond checks, and verification of (H1)–(H5) all hold. Note that `paper/friend-tree-exponent.tex` lines 1136–1142 faithfully state all five hypotheses (boundedness $`f \le F`$ and existence of $`f_{\min} > 0`$ follow automatically from finitely many nonnegative values plus positive newborn weight).
   - *Verdict:* **holds**.
2. **Lines 1150–1170 (Growth on the polynomial scale; `proofs/polynomial-scale-lower.md` §§1–5):**
   - *Claim:* Combining baseline lower ($`\delta=543/1000`$) and upper ($`\lambda=293/500`$) tables (with $`\lambda < 2\delta`$) yields $`\liminf C_n/n^\delta > 0`$ and $`\lim \Phi_n^-/n^\delta \in (0,\infty]`$ a.s.
   - *Check:* In `proofs/polynomial-scale-lower.md`, §1 (jump size and conditional second-moment bound $`O(C_n/n)`$ via $`\sum \rho_u \le C_n`$), §2 (restriction across a degree-$`H`$ saturated parent edge, preserving full-tree per-vertex inequalities (1) because baseline tables use no inter-vertex transfers and $`p`$ stays permanently saturated), §3 (product bounds $`e^{-a^2}(n/m)^a \le P_a(m,n) \le e^a(n/m)^a`$, conditional $`L^2`$-boundedness of the centered martingale part $`Z_n`$ using $`q = 2\delta-\lambda > 0`$, and conditional Cauchy–Schwarz lower bound $`p_* > 0`$ on $`\{J > a/2\}`$), §4 (existence of a reachable depth-$`\ge 2`$ positive-weight pattern since expected depth-$`\le 1`$ cores are $`O(\log n)`$, uniform probability $`\ge 2^{-H}`$ of reaching saturation + a.s. promotion + uniform replay probability $`p_{\mathrm{pattern}} > 0`$), and §5 (0–1 lifting via algebra approximation) all hold. The paper accurately qualifies that $`C_n/n^\delta \to \infty`$ already follows from $`11/20 > \delta`$ and that the distinct content is the loss-free convergence of $`\Phi_n^-/n^\delta`$.
   - *Verdict:* **holds**.
3. **Lines 1172–1186 (Cumulative mean degree tail; `proofs/degree-tail.md` Theorem D(a)–(d)):**
   - *Claim:* The two displayed bounds on $`\mathbb{E}N_{\ge k}(n)`$ imply $`\lim_{n/k\to\infty} (\log \mathbb{E}N_{\ge k}(n))/\log(n/k) = \mu`$ uniformly over $`k\ge 2`$ given Theorems 1.1 and 1.2.
   - *Check:* In `proofs/degree-tail.md`, Lemma 1 (linear first-moment and quadratic $`D(D+1)`$ recursion), Part (a) (cohort decomposition + Markov), Lemma 2 (take-off via $`\mathbb{E}_j[L_u(j+1)-1] \ge (L_u(j)-1)(j+1)/j`$ and Paley–Zygmund $`\ge 1/30`$), Lemma 3 ($`\ge C/4`$ eligible vertices), Lemma 4 (activation probability $`\ge 1/368`$; note for Case B, $`\frac{m'(m'-2)}{(2m'-2)(2m'-3)} \ge \frac{1}{4}`$ holds for all $`m'\ge 3`$ since $`4m'(m'-2)-2(m'-1)(2m'-3) = 2m'-6 \ge 0`$), Part (b) ($`4\cdot 368\cdot 30 = 44160`$), Part (c) (summation by parts properly separating $`k=2`$, $`m_0=n`$), and Part (d) (uniform logarithmic limit for $`\mu\in(0,1)`$) all hold.
   - *Verdict:* **holds**.
4. **Lines 1187–1189 (Inline proof: fixed-$`k`$ a.s. limit for $`N_{\ge k}(n)`$):**
   - *Claim:* Iterating `[AB, Theorem 4.8]` from $`N_{\ge 2}(n)=C_n`$ gives $`(\log N_{\ge k}(n))/\log n \to \mu`$ a.s. for each fixed $`k\ge 2`$.
   - *Check:* Theorem 4.8 of Addario-Berry et al. gives deterministic constants $`c_j > 0`$ with $`c_j < \liminf_{n\to\infty} N_{\ge j+1}(n)/N_{\ge j}(n) \le 1`$ a.s. for every $`j\ge 2`$. Multiplying for $`j=2,\dots,k-1`$ gives $`\prod_{j=2}^{k-1} c_j \le \liminf N_{\ge k}(n)/C_n \le \limsup N_{\ge k}(n)/C_n \le 1`$ a.s.; taking logs and dividing by $`\log n`$ yields $`(\log N_{\ge k}(n) - \log C_n)/\log n \to 0`$ a.s., so Theorem 1.1 implies $`(\log N_{\ge k}(n))/\log n \to \mu`$ a.s.
   - *Verdict:* **holds**.

---

### 2.5 Sections 7–8 and Production Disclosure (`paper/friend-tree-exponent.tex`, lines 1191–1306)

- **Section 7 (Numerical measurements, lines 1191–1271):**
  - The identity $`d\mathbb{E}[C(t)]/dt = \mathbb{E}[S(t)]`$ (where $`S(t) = \sum_{u\in K} L_u/D_u`$) and $`d\log\mathbb{E}[C(t)]/d\log\mathbb{E}[N(t)] = \mathbb{E}[S(t)]/\mathbb{E}[C(t)]`$ for the $`T_3`$-started continuous-time process ($`\mathbb{E}[N(t)]=3e^t`$) hold rigorously (bounded intensity $`S(t) \le N(t)`$).
  - Every reported fit, jackknife standard error, firing count, slope, and mean population in Table 1 and the prose was independently recomputed from the raw archived files (`simulations/fixed-size/core_counts.tsv.gz`, `simulations/outputs/ct_t28_seedA.txt.gz`, `ct_t28_seedB.txt.gz`, `ct_t31_seedC.txt.gz`, `simulations/route-exact/replay.txt`) and matches to the last digit (Section 3.5).
  - The Gaussian approximation thresholds for large Poisson intensities (`10^9` in `simulations/sim/friend_ct.cpp` line 83, `10^{14}` in `simulations/route-exact/sim_lazy_exact.c` line 105) are verified in source and candidly disclosed. (Minor technical note: `sim_lazy_exact.c` lines 40–44 uses rejection-free 128-bit multiplication for uniform integers, introducing negligible bias $`< D/2^{64}`$, whereas `friend_ct.cpp` lines 48–57 uses exact rejection Lemire sampling.)
  - The separation between simulations and proofs is completely respected ("The simulations provide no bound on the distance from these estimates to $`\mu`$").
- **Section 8 (Open problems, lines 1272–1280):**
  - Accurately delineates what remains open: exact value of $`\mu`$, sharp asymptotics $`\mathbb{E}[C_n] \sim c n^\mu`$ and convergence of $`C_n/\mathbb{E}[C_n]`$, uniform growing-$`k`$ sample-path degree tail, and pointwise scaling of $`N_{=k}(n)`$.
- **Production disclosure (`How this work was produced`, lines 1282–1306):**
  - Unusually thorough and transparent, accurately matching the authorship headers across `README.md`, `REPORT.md`, and `proofs/`.
  - *Editorial note for AAP/EJP:* The closing sentence "No human referee has yet read the proofs" should be complemented upon submission/revision by an explicit statement of the human author's personal verification and scholarly responsibility for all proofs and code, in line with IMS / journal guidelines on AI-assisted research.

---

## 3. Independent Computations, Small-Case Checks, and Certificate Audits

All computations below were performed independently without modifying the assigned snapshot.

### 3.1 Exact rational $`\mathbb{E}[C_n]`$, Theorem 1.2 arithmetic, and §4 band state counts
1. **Exact arithmetic checks on Theorem 1.2 and §4 tail:**
   - Exponent margins below $`17/30`$: $`17/30 - 1133/2000 = 1/6000`$; $`17/30 - 5657/9983 = 1/299490`$.
   - Reduced prefactors: $`\gcd(14608704616, 17527873930) = 2`$ giving $`7304352308/8763936965`$; $`\gcd(14606753031, 2989767809) = 1`$; $`\gcd(1539738114735, 3908512529) = 1`$.
   - §4 tail margins (`paper/friend-tree-exponent.tex`, lines 959–961): with $`s=10^{12}`$, $`a_*=398437500000`$, $`A=1380739398522`$, $`B=383976882193`$,
     ```math
     5657 A - 9983(B + a_*) = 6235, \qquad 5657(A + B) - 9983 s = 4755,
     ```
     matching the paper verbatim.
2. **Independent exact enumeration of $`\mathbb{E}[C_n]`$ for $`n = 3, \dots, 9`$** (pure Python `fractions.Fraction` over canonical unlabeled trees starting from $`T_2`$), compared against `verify/motif-bracket/exact_rationals_n21.txt` (all five columns—unlabeled tree count, $`\mathbb{E}[C_n]`$, $`\mathbb{E}[S_n]`$, $`\mathbb{E}[S_n]/\mathbb{E}[C_n]`$, $`\mathrm{Var}(C_n)`$—match identically) and all three Theorem 1.2 product bounds (verified for all $`n = 3, \dots, 21`$):

| $`n`$ | Trees | Exact $`\mathbb{E}[C_n]`$ | Decimal | Matches `exact_rationals_n21.txt` | Motif LB | Motif UB | Banded UB |
|---:|---:|---|---:|:---:|---:|---:|---:|
| 3 | 1 | $`1`$ | `1.000000` | Exact (all 5 cols) | `0.833456` | `4.885581` | `393.944782` |
| 4 | 2 | $`4/3`$ | `1.333333` | Exact (all 5 cols) | `0.990423` | `5.808142` | `468.356347` |
| 5 | 3 | $`19/12`$ | `1.583333` | Exact (all 5 cols) | `1.130320` | `6.630720` | `534.707798` |
| 6 | 6 | $`647/360`$ | `1.797222` | Exact (all 5 cols) | `1.258046` | `7.382012` | `595.308295` |
| 7 | 11 | $`17189/8640`$ | `1.989468` | Exact (all 5 cols) | `1.376512` | `8.079064` | `651.535648` |
| 8 | 23 | $`4717423/2177280`$ | `2.166659` | Exact (all 5 cols) | `1.487616` | `8.732961` | `704.279739` |
| 9 | 47 | $`101568829/43545600`$ | `2.332471` | Exact (all 5 cols) | `1.592679` | `9.351395` | `754.165656` |
| 10 | 106 | (from file) | `2.489218` | Verified bounds | `1.692665` | `9.939775` | `801.640653` |
| 21 | 2,144,505 | (from file) | `3.881185` | Verified bounds | `2.590654` | `15.229864` | `1228.433228` |

3. **Independent combinatorial count of all §4 band states and `recipe.txt.gz` integrity:**
   - SHA-256 of `verify/banded-upper/recipe.txt.gz`: `584bb93c2797138ae3de54cb85c7d8e487bb7d7cc5320c021b5a64ec0cf0f60a` (exact match with `paper/friend-tree-exponent.tex`, lines 974–975).
   - Record counts in `recipe.txt.gz`: `scale = 1000000000000`, `threshold = 8`, `12,697` `reward` entries, `63,309` `edge` entries, `18` `arrow` entries (max `398437500000` at `(2, 2)`), matching lines 971 and 959.
   - Pure stars-and-bars combinatorial counts from `verify/banded-upper/geometry.txt` (`c=1` label count `12,697`):
     - $`[2,16)`$: `53,201,381`; $`[16,32)`$: `23,199,024`; $`[32,80)`$: `22,401,552`; $`[80,200)`$: `17,550,960`; $`[200,500)`$: `22,446,000`; $`[500,1000)`$: `21,291,500`; $`[1000,2400)`$: `66,455,200`; total = **`226,545,617`** (exact match on every row of the table on lines 741–747).

### 3.2 Full UBSan replay of all headline certificate verifiers
1. **Motif GMP verifier (`verify/motif-bracket/cert_build.c` + `cert_verify_gmp.c` compiled with `-O3 -fopenmp -fsanitize=undefined -fno-sanitize-recover=all`):**
   - Rebuilt binary SHA-256 digests:
     - UPPER (`1133/2000`): `3a4aef888a5b567a8ada05f9546f98c979d361da9f9b33ba6e794eb91fb30abe`
     - LOWER (`113/200`): `9e966fba9aee6eae3da659b693d3bb9c481c1cea66635b5218d0887f182a60d9`
     Both match `proofs/bracket-motifs.md` line 372 identically.
   - UPPER check: `f: 484,138`, `g: 160,545`, `hV: 13,519,800`, `hC: 3,164,625`, total **`17,329,108`** exact GMP checks, `0` failures, `0` UB traps; `Phi(T3) = 14606753031/10^{10}`, `min W(u) = 2989767809/10^{10}`, exact prefactor `14606753031/2989767809`; 37/37 built-in exact brute-force trees PASS.
   - LOWER check: identical `17,329,108` checks, `0` failures, `0` UB traps; `Phi(T3) = 14608704616/10^{10}`, `max W(u) = 17527873930/10^{10}`, exact prefactor `7304352308/8763936965`; 37/37 built-in exact brute-force trees PASS.
   - All 6 negative controls (`--override-pq 1132 2000`, `--override-pq 114 200`, `--perturb f|g|hV|hC`) rejected with exit code `2` (for example, `--perturb g` triggered 1 failure on `g(p=276 -> w=0)`, `--perturb hV` 1 failure on `hV(p=1, u=(0,1), w=0)`, `--perturb hC` 1 failure on `hC(u=(0,2,1), w1=0, w2=0)`).
2. **Primary banded verifier (`verify/banded-upper/check.cpp` compiled with `-std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all`):**
   - `verified rows 226545617, signed transfers 63309, minimum integer drift slack 3`
   - `scale 1000000000000, minimum newborn reserve 4340525, weight range [3908512529, 1764716280715], centre 1539738114735`
   - `PASS upper 5657/9983` (0 UB traps, exact match with `paper/friend-tree-exponent.tex`, lines 1024–1029).
3. **Independent banded verifier (`verify/banded-upper/independent/verify_banded.cpp`, same UBSan flags):**
   - Verified all 7 bands (`226,545,617` rows), `min finite slack = 0`, `newborn consistency: min(trial - actual) = 4341454` (worst $`c=1`$ index 286: actual `1015179606101` vs. trial `1015183947555`), `global [3908512428, 1764716280715]`, root weight `1539738113358`, `T3 exact scaled drift / root weight = 872508277036/1539738113358 = 0.566660180369`, `PASS upper 5657/9983` (exact match with `paper/friend-tree-exponent.tex`, lines 1034–1036).

### 3.3 Independent small-case motif decomposition checks (pure Python `Fraction`)
An independent pure-Python implementation of `paper/friend-tree-exponent.tex` §§3.2–3.3 (reading the binary certificate arrays directly, evaluating $`\Phi_\pm(T)`$, computing brute-force $`\mathcal{L}\Phi_\pm(T)`$ over all attachment targets, and separately evaluating every local expression $`V_u, E_e, P_{puw}, B_{u12}`$ from Eqs. (2)–(5)) confirmed `L_bf == V_sum + E_sum + P_sum + B_sum` identically (zero rational discrepancy) and verified every individual local term inequality ($`V_u \preceq \lambda f`$, $`E_e \preceq \lambda g`$, $`P_{puw} \preceq \lambda h_V`$, $`B_{u12} \preceq \lambda h_C`$) on all 9 test trees covering $`T_3, K_{1,3}, P_4`$, both $`P_5`$ rootings, active cherries, active paths + cherries, and the $`D_0=24`$ tail seam:

| Tree | Dir | Exact / Decimal $`\Phi`$ | Exact / Decimal $`\mathcal{L}\Phi`$ | Margin ($`\lambda\Phi-\mathcal{L}\Phi`$ or rev.) | Decomposition $`(V, E, P, B)`$ |
|---|:---:|---|---|---:|---|
| $`T_3`$ | UP | $`\frac{14606753031}{10^{10}}`$ (`1.46067530`) | $`\frac{1654715823}{2\cdot 10^9}`$ (`0.82735791`) | $`\frac{2292954123}{2\cdot 10^{13}}`$ (`1.1465e-4`) | `(0.827358, 0, 0, 0)` |
| $`T_3`$ | LO | $`\frac{1826088077}{125\cdot 10^7}`$ (`1.46087046`) | $`\frac{1650907733}{2\cdot 10^9}`$ (`0.82545387`) | $`\frac{3878481}{625\cdot 10^8}`$ (`6.2056e-5`) | `(0.825454, 0, 0, 0)` |
| $`K_{1,3}`$ | UP | $`\frac{3882153933}{25\cdot 10^8}`$ (`1.55286157`) | $`\frac{4397906807}{5\cdot 10^9}`$ (`0.87958136`) | $`\frac{573599089}{5\cdot 10^{12}}`$ (`1.1472e-4`) | `(0.879581, 0, 0, 0)` |
| $`K_{1,3}`$ | LO | $`\frac{7767119667}{5\cdot 10^9}`$ (`1.55342393`) | $`\frac{1755494601}{2\cdot 10^9}`$ (`0.87774730`) | $`\frac{62778129}{10^{12}}`$ (`6.2778e-5`) | `(0.877747, 0, 0, 0)` |
| $`P_4`$ | UP | $`\frac{2629575843}{125\cdot 10^7}`$ (`2.10366067`) | $`\frac{23831857331}{2\cdot 10^{10}}`$ (`1.19159287`) | $`\frac{639187}{4882812500}`$ (`1.3091e-4`) | `(1.161190, 0.030403, 0, 0)` |
| $`P_4`$ | LO | $`\frac{4202434769}{2\cdot 10^9}`$ (`2.10121738`) | $`\frac{5936306411}{5\cdot 10^9}`$ (`1.18726128`) | $`\frac{29383983}{4\cdot 10^{11}}`$ (`7.3460e-5`) | `(1.158509, 0.028753, 0, 0)` |
| $`P_5`$ (centre) | UP / LO | `2.64555699` / `2.63634602` | `1.49856547` / `1.48961746` | `1.4256e-4` / `8.1955e-5` | `B = -0.005479 / -0.007857` |
| $`P_5`$ (endcore) | UP / LO | `2.64482651` / `2.63618958` | `1.49815315` / `1.48952797` | `1.4107e-4` / `8.0856e-5` | `P = -0.003488 / -0.005416` |
| Cherry ($`n=7`$) | UP / LO | `3.71370583` / `3.69034130` | `2.10360588` / `2.09569110` | `2.0848e-4` / `1.0648e-2` | Both `P` and `B` active |
| Path+Cherry ($`n=8`$) | UP / LO | `4.24400782` / `4.20268563` | `2.40151772` / `2.39350022` | `2.7127e-3` / `1.8983e-2` | Both `P` and `B` active |
| Star $`D=24`$ | UP / LO | `1.66263735` / `1.66298661` | `0.94177145` / `0.93964997` | `1.1261e-4` / `6.2533e-5` | Tail seam $`D=D_0`$ |
| Edge $`(23,24)`$ | UP / LO | `3.27241774` / `3.27137752` | `1.85290879` / `1.84851853` | `9.1586e-4` / `1.9023e-4` | Inner–tail edge (`E<0`) |

### 3.4 Independent small-case banded history drift checks ($`n = 3, 4, 5, 6`$)
Using `verify/tree-drift/read_weights.cpp` (compiled with UBSan) to extract verified weights from `check.cpp` and evaluating the exact rational one-step generator over all $`1 + 3 + 12 + 60 = 76`$ attachment histories from $`T_3`$ through $`n = 6`$ (covering 43 distinct finite states):

- Extracted small-state weights: $`z_{(2,0)} = 1{,}539{,}738{,}114{,}735`$, $`z_{(3,0)} = 1{,}637{,}065{,}521{,}366`$, $`z_{(2,1;(2,1))} = 1{,}108{,}795{,}790{,}246`$, $`z_{(2,1;(2,2;2,1))} = 1{,}120{,}441{,}477{,}455`$.
- Exact hand-verified $`T_3`$ drift:
  ```math
  \mathcal{L}\Phi(T_3) = 2 z_{(3,0)} + 2 z_{(2,1;(2,1))} - 3 z_{(2,0)} = 872{,}508{,}279{,}019,
  ```
  with $`9983\,\mathcal{L}\Phi(T_3) = 8{,}710{,}250{,}149{,}446{,}677 < 8{,}710{,}298{,}515{,}055{,}895 = 5657\,\Phi(T_3)`$, ratio $`0.566660181 < 5657/9983 \approx 0.566663328`$.
- At $`n=4`$: $`K_{1,3}`$ (`path=(0,)`): $`\Phi = 1{,}637{,}065{,}521{,}366`$, $`\mathcal{L}\Phi = 927{,}660{,}174{,}901`$, ratio $`0.566660383`$; $`P_4`$ (`path=(1,)` and `(2,)`): $`\Phi = 2{,}217{,}591{,}580{,}492`$, $`\mathcal{L}\Phi = 1{,}256{,}622{,}945{,}165`$, ratio $`0.566661127`$.
- Across all $`12`$ histories at $`n=5`$ and $`60`$ histories at $`n=6`$: `0` failures, maximum ratio $`\mathcal{L}\Phi/\Phi = 0.566661576`$ ($`n=5`$) and $`0.566661369`$ ($`n=6`$), both strictly below $`5657/9983`$.

### 3.5 Independent verification of §7 simulation summaries
Checksum verification (`sha256sum -c simulations/SHA256SUMS`) passed on all 13 archived simulation files. Independent Python recomputation directly from the raw gzipped trajectory tables matched every number reported in §7:

- Fixed-size 2,400 runs (`simulations/fixed-size/core_counts.tsv.gz`): count-unit nonlinear least-squares fit $`0.5660272 \pm 0.0002417`$, log-unit linear fit $`0.5675481 \pm 0.0004706`$, adjacent slopes $`0.569271 \pm 0.003482`$ ($`100\to 316`$) and $`0.566350 \pm 0.000374`$ ($`316228\to 10^6`$).
- Continuous-time 912 runs (`ct_t28_seedA/B`, `ct_t31_seedC`): $`t=4,\dots,12`$ fits vs. $`3e^t`$ give $`0.5655751 \pm 0.0004760`$ (count units) and $`0.5676083 \pm 0.0008194`$ (log units); $`t=31`$ (91 runs) gives mean size $`8.1167\times 10^{13}`$ and slope $`0.5657203 \pm 0.0000038`$; pooled $`t>21`$ gives $`22{,}223{,}928{,}834`$ firings, $`A/F = 0.5657230 \pm 0.0000032`$, and $`Cv/F = 0.5657253 \pm 0.0000032`$.
- Separate 48-run replay (`simulations/route-exact/replay.txt`): $`0.5657242 \pm 0.0000175`$ at $`t=27`$.

---

## 4. Which Finite Checks Imply Which Claims

A computer-assisted proof should separate **proof-bearing finite checks** (required logically by a theorem) from **corroborative diagnostic checks** (which test code against independent formulas but do not enter the proof).

| Claim in paper | Proof-bearing analytic ingredients | Proof-bearing finite checks (trusted code + data) | Untrusted / diagnostic only (NOT needed for proof) |
|---|---|---|---|
| **Theorem 1.1** (existence of $`\mu`$ & $`\log C_n/\log n\to\mu`$ a.s.) | Entire §2 (`paper/friend-tree-exponent.tex`, lines 63–161) | **None** (purely analytic) | `verify/tree-drift/existence_counterexample.py` (illustrates why unweighted free comparison fails) |
| **Theorem 1.2, Motif bracket** ($`113/200`$ lower, $`1133/2000`$ upper) + a.s. $`\limsup \le 1133/2000`$ | §3.1–3.5: Lemmas 3.1–3.2, five vertex cases, edge/path/cherry corner reductions, vertex-grouped core comparison, plus header constraint $`D_0 \le K_F + 1`$ | `verify/motif-bracket/cert_verify_gmp.c` checks `[1/5]`–`[5/5]` (`17,329,108` GMP checks/side, $`g[\text{tail},\text{tail}]=0`$, $`c_0, C_0, \Phi(T_3) > 0`$) on the two binary tables, with `Domains` confirmed ($`D_0=24\le K_F+1=151`$) and zero signed overflow (UBSan) | `cert_build.c` (search/proposal only); `[6/6]` 37 brute-force trees; `--override-pq` and `--perturb` negative controls; `verify_python_cross.py`; `exact_rationals_n21.txt`; `motif_identity.py` |
| **Theorem 1.2, Banded upper bound** ($`5657/9983`$) + a.s. $`\limsup \le 5657/9983`$ | §4.1–4.4: Lemma 4.1, two exact cancellations, finite row bounds, $`D\ge H`$ affine bound, DAG ordering, newborn consistency | **Either** `verify/banded-upper/check.cpp` **or** independent `verify/banded-upper/independent/verify_banded.cpp` (with UBSan) on `geometry.txt` and `recipe.txt.gz`: verifies observability of all `63,309` edges, `18` arrows, `226,545,617` finite rows, newborn reserve $`\ge 0`$, and the 2 tail inequalities | Floating-point proposals inside `check.cpp`/`verify_banded.cpp` (only select initial trial $`z`$; exact 128-bit slack decides acceptance); `test_banded_audit.py` (smoke/histories/negative controls); `verify/tree-drift/banded_histories.py` |
| **Combined $`113/200 \le \mu \le 1133/2000 < 17/30`$** | Theorem 1.1 + Theorem 1.2 (motif) | Same as motif bracket above | All simulations in §7 |
| **§6 Supplementary bounds** ($`543/1000`$ a.s. & polynomial scale, $`[11/20, 29/50]`$ cross-charging, $`\mathbb{E}N_{\ge k}(n)`$ degree tail) | `proofs/almost-sure-lower.md`, `proofs/cross-charging.md`, `proofs/polynomial-scale-lower.md`, `proofs/degree-tail.md` | Baseline `certs9.py` / `find_upper_lin.py` (`as_lower_check.py`), cross-charging `xverify.py` (degree tail Theorem D is analytic given any interior bracket $`\mu\in(0,1)`$) | `xcharge.py` LP solver & 200-tree test; `verify/polynomial-lower/check_as_lower.py` tree/branch tests; transfer/hub suites (`563/1000`, `568/1000`, `56/100`, `57/100`) are not cited in the paper's theorems |

---

## 5. Independent Verification of Bibliography and Literature Coverage

Every bibliography entry in `paper/friend-tree-exponent.tex` (lines 1307–1331) was fetched independently from primary sources (arXiv TeX sources, Crossref API metadata, X oEmbed API, and GitHub commit API).

### 5.1 Entry-by-entry verification

1. **`[SK]` J. Saramäki and K. Kaski (2004).**
   - *Cited data:* *Physica A* **341** (2004), 80–86; `arXiv:cond-mat/0404088`.
   - *Primary-source verification:* Crossref (`DOI:10.1016/j.physa.2004.04.110`) confirms exact authors, title ("Scale-free networks generated by random walkers"), journal (*Physica A: Statistical Mechanics and its Applications*), volume 341, pages 80–86, October 2004. arXiv source `cond-mat/0404088` (`saramaki_randomwalk_rev1.tex`, lines 117–135) defines the growth rule where a new vertex initiates a random walk of length $`l`$ from a uniform vertex to mark $`m`$ attachment targets; setting $`m=1, l=1, m_0=2`$ recovers the random friend tree.
   - *Assertion check (line 29):* Accurate. (Note that Saramäki–Kaski heuristically claimed BA scaling $`\gamma = 3`$ even for $`l=1`$, which `[CJ]` later refuted for $`m=1, l=1`$.)
   - *Access limitation:* Verified from arXiv TeX source and Crossref metadata; Elsevier publisher PDF behind paywall.

2. **`[CJ]` C. Cannings and J. Jordan (2013).**
   - *Cited data:* *Electronic Communications in Probability* **18** (2013), paper 77, 1–5; `doi:10.1214/ECP.v18-2518`.
   - *Primary-source verification:* Crossref confirms authors, title ("Random walk attachment graphs"), journal (*Electronic Communications in Probability*), volume 18, pages 1–5, 2013, DOI `10.1214/ecp.v18-2518`; paper no. 77 corroborated via bibliographic index records (`arXiv:1303.1052`). arXiv TeX source (`1303.1052`, Section 3, Theorem 1, `\label{l1}`, line 102) states verbatim: *"When $`\ell = 1`$, as $`n\to\infty`$, $`L_n \to 1`$, almost surely."*
   - *Assertion check (line 29):* Accurate.
   - *Access limitation:* Direct HTTP request to `doi.org/10.1214/ECP.v18-2518` (Project Euclid) was blocked by an Incapsula anti-bot challenge (not bypassed); full text was verified from `arXiv:1303.1052v2` TeX source and Crossref metadata.

3. **`[KR]` P. L. Krapivsky and S. Redner (2017).**
   - *Cited data:* *Journal of Statistical Mechanics: Theory and Experiment* (2017), 073405; `arXiv:1706.01514`.
   - *Primary-source verification:* Crossref (`DOI:10.1088/1742-5468/aa7a3f`) confirms authors, title ("Emergent network modularity"), journal, year 2017, issue 7, article 073405. arXiv source (`1706.01514v2`, `paper.tex`) confirms every single detail:
     - Section 1 (line 105) coins *isotropic redirection (IR)*; lines 131–140 define the *nucleus* ($`\mathcal{N} = \sum_{k\ge 2} N_k`$, with $`N_k`$ counting degree **exactly** $`k`$) and propose $`N_k \sim N^\mu / k^{1+\mu}`$.
     - **Section 5.1** (*Sublinear Growth*, lines 633–670) states eq. (core) $`\mathcal{N}\sim N^\mu, N_k\sim N^\mu`$ with $`\mu\approx 0.566`$, explicitly reports slopes *"between 0.5652 and 0.5673"* (line 650), and eq. (ck) $`c_k = N_k/\mathcal{N} \sim k^{-(1+\mu)}`$.
   - *Assertion check (lines 29, 1193, 1278):* Exact match on every section number, notation, and numerical figure.
   - *Access limitation:* IOP publisher PDF behind paywall; verified from arXiv v2 TeX source and Crossref.

4. **`[Strogatz]` S. Strogatz (20 September 2022).**
   - *Cited data:* Post on the exponent in isotropic redirection, 20 September 2022, `https://x.com/stevenstrogatz/status/1572287647715934210`.
   - *Primary-source verification:* X oEmbed confirms the author, date, conversation with Redner and reported exponent of 0.566. The pertinent sentence is quoted below.
   - *Assertion check (line 29):* Date, author, Redner attribution, and $`17/30`$ problem confirmed. Minor wording nuance: Strogatz wrote *"It's an open problem to prove that 17/30 is really the exponent"* (slightly more affirmative than "asked whether $`17/30`$ could be the exponent").
   - *Access limitation:* The attached image (`pic.twitter.com/EqQoeRP4GD`) and any thread replies could not be fetched without X authentication/JS rendering, so the statement on line 29 that $`N_k`$ in Strogatz's post/image denotes exact degree $`k`$ could not be directly verified from the image itself (though it matches `[KR]`, from which Redner shared the problem).

5. **`[AB]` L. Addario-Berry, S. Briend, L. Devroye, S. Donderwinkel, C. Kerriou, and G. Lugosi.**
   - *Cited data:* `Random friend trees (2024). arXiv:2403.20185`.
   - *Bibliographic defect found:* Crossref confirms this paper has been published in ***Random Structures & Algorithms*** **69** (2026), no. 1, e70082, `DOI:10.1002/rsa.70082` (published online July 30, 2026). The manuscript should cite the published journal version alongside the arXiv identifier (notably, `REPORT.md` line 243 already records the RSA citation).
   - *Assertion check (lines 31 and 1188, verified in `2403.20185v1`, `template.tex`):*
     - **Proposition 6.20** (line 1153, `\label{prop:urns}`): $`\alpha = (\sqrt{13}-3)/2`$, for $`0 < \delta < \alpha/3 = (\sqrt{13}-3)/6 \approx 0.100925`$, $`n^{-\delta} X_n^{\ge 2} \to \infty`$ a.s. Exact match!
     - **Proposition 6.22** (line 1205, `\label{lem:upperbound_nonleaves}`): $`\gamma = 3 - 2\sqrt{2}`$, $`\delta = 1 - \gamma/2 = \sqrt{2} - 1/2 \approx 0.914214`$, for $`\delta < \lambda < 1`$, $`n^{-\lambda} X_n^{\ge 2} \to 0`$ a.s. Exact match!
     - **Theorem 4.8** (line 409, `\label{thm:all_finite_deg_large}`): there exists a positive sequence $`(c_k)_{k\ge 2}`$ such that $`c_k < \liminf_{n\to\infty} X_n^{\ge k+1}/X_n^{\ge k} \le 1`$ a.s. Exact match!
   - *Access limitation:* Wiley published RSA PDF is paywalled; verified against `arXiv:2403.20185v1` TeX source and Crossref. The authors should verify whether proposition/theorem numbering changed between `2403.20185v1` and the published RSA 2026 article.

6. **`[HKRZ]` H. Hartle, P. L. Krapivsky, S. Redner, and Y. Zhang (2026).**
   - *Cited data:* *Anomalous scaling in redirection networks*, *Physical Review E* **114** (2026), 034314; `arXiv:2604.01540`.
   - *Primary-source verification:* Crossref (`DOI:10.1103/9lmp-pydp`) confirms exact authors, title, *Physical Review E* vol. 114, issue 3, article 034314, September 2026. arXiv source (`2604.01540v1`, `paper.tex`) confirms:
     - Abstract and Section 2–3 introduce leaf-redirecting variants (DAN, PAN, VAN($`\omega`$)) admitting closed leaf-degree recurrences (Eqs. `M-ell`, `q-ell`).
     - **Figure 6** (`\label{fig:n1_growth}`, lines 271–278) caption gives the mean number of core vertices for the three models over 2,400 independent realizations per point, sizes $`N=10^2`$ to $`N=10^6`$, and the isotropic-redirection fit of 0.569. Weighting is unspecified; note also line 514 confirms $`L_3`$ (the 3-vertex path) is their canonical initial condition.
   - *Assertion check (lines 31, 1194–1197, 1230–1231):* Exact match.
   - *Access limitation:* APS publisher PDF behind paywall; verified from arXiv v1 TeX source and Crossref.

7. **`[supplement]` C. Uzay, GitHub revision `c4653036b993bf81799c51b60cdd18a500142147`.**
   - *Primary-source verification:* GitHub API confirms commit `c4653036b993bf81799c51b60cdd18a500142147` exists (HTTP 200, 168 tree entries). Comparing all git blob hashes at `c4653036b993` against the assigned snapshot `526f350` shows that every file under `proofs/`, `verify/`, and `simulations/` (except `simulations/README.md`) is byte-for-byte identical, while `paper/` itself was added only after `c4653036b993`.
   - *Archival issue:* Citing a GitHub commit earlier than the manuscript's own revision without a permanent DOI (e.g., Zenodo or journal-hosted supplementary archive) does not meet archival permanence standards for *AAP/EJP*.

### 5.2 Missing relevant work
For a probability journal such as *Annals of Applied Probability* or *Electronic Journal of Probability*, the 7-entry bibliography is unduly sparse. The following verified works directly belong in the introduction and open-problems discussion:

1. **T. S. Evans and J. P. Saramäki (2005)**, *Scale-free networks from self-organization*, *Phys. Rev. E* **72**(2), 026138 (`DOI:10.1103/PhysRevE.72.026138`) — directly analyzed random-walk attachment including length-$`1`$ walks; cited by both `[CJ]` and `[AB]`.
2. **P. L. Krapivsky and S. Redner (2001)**, *Organization of growing random networks*, *Phys. Rev. E* **63**(6), 066123 (`DOI:10.1103/PhysRevE.63.066123`) — foundational paper introducing redirection in growing random trees.
3. **A. Vázquez (2003)**, *Growing network with local rules: Preferential attachment, clustering hierarchy, and degree correlations*, *Phys. Rev. E* **67**(5), 056104 (`DOI:10.1103/PhysRevE.67.056104`) — early random-walk / local-rule network growth model preceding `[SK]`.
4. **A. Gabel, P. L. Krapivsky, and S. Redner (2013, 2014)**, *Highly dispersed networks by enhanced redirection*, *Phys. Rev. E* **88**(5), 050802 (`DOI:10.1103/PhysRevE.88.050802`) and *J. Stat. Mech.* (2014), P04009 (`DOI:10.1088/1742-5468/2014/04/P04009`) — established the $`N_k \sim N^\mu / k^{1+\mu}`$ sublinear-nucleus phenomenology and lack of self-averaging cited directly in `[KR]` Section 5.1.
5. **P. L. Krapivsky and S. Redner (2024)**, *The magic of networks grown by redirection*, *Indian J. Phys.* **98**, 3777 (`DOI:10.1007/s12648-023-02874-x`, `arXiv:2305.10628`) — survey highlighting isotropic redirection and open problems on $`\mu`$ (already noted in `REPORT.md`, omitted from `paper/friend-tree-exponent.tex`).
6. **S. Janson (2004)**, *Functional limit theorems for multitype branching processes and generalized Pólya urns*, *Stochastic Process. Appl.* **110**(2), 177–245 (`DOI:10.1016/j.spa.2003.12.002`) and **C. Holmgren and S. Janson (2017)**, *Fringe trees, Crump–Mode–Jagers branching processes and $`m`$-ary search trees*, *Probab. Surveys* **14**, 53–154 (`DOI:10.1214/16-PS272`) — standard probability literature on continuous-time branching-process embeddings of growing random trees (providing mathematical context for §2 and the open problem of $`C_n/\mathbb{E}[C_n]`$ convergence in §8).

---

## 6. Exposition, Notation, and Presentation Issues

1. **Undefined notation typo (line 236):** the phrase "under the measures $`\pi^{\mathrm{end}}`$" (TeX `\pi^{\rm end}`) should read "under the measures $`\zeta_s`$".
2. **Symbol overloading across §§2–4:**
   - $`H`$ denotes the path-centre domain in §3 (line 226) and the degree cutoff $`2400`$ in §4 (line 731).
   - $`h`$ denotes preparation time $`k+1`$ in §2 (line 86), parent indicator $`h_u\in\{0,1\}`$ in §3 (line 177), and neighbour cutoffs $`h_c`$ in §4 (line 732).
   - $`b`$ denotes $`7/\log(1/\epsilon)`$ in §2 (line 113), budding rate $`b_u`$ in §§3–4, and pooled lower-bound mark in §4 (line 761).
   - $`T`$ denotes a tree throughout §§1–3, and the arrow threshold $`T = 8`$ in §4 (line 853).
   - $`A`$ denotes $`8`$ in §2 (line 113), knot functions $`A_{pc}, A_{cp}, A_{cc}`$ in §3 (line 251), and the affine tail intercept in §4 (line 943).
   A brief notation warning at the start of §§3–4 (or renaming the arrow threshold $`T=8`$ to e.g. $`D_{\mathrm{arr}}=8`$ and $`H`$ in §3 to $`\mathcal{H}`$) would greatly improve readability.
3. **Awkward phrasing in §4.3 (line 948):** "Thus the rate times the mean increment of $`q`$, summed over the two events" should read "Thus the scaled drift of $`q`$, summed over the two events" (each event has a deterministic increment of $`q`$, weighted by its rate).
4. **Cross-reference hygiene (`paper/friend-tree-exponent.tex`):** Seven `\label`s (`lem:plant`, `lem:decompose`, `sec:motif-reduction`, `lem:banded-certificate`, `sec:checks`, `sec:related`, `sec:sim`) are defined but never referenced; in §3.4 line 585, citing Lemma 3.1 (`\ref{lem:decompose}`) alongside `\eqref{eq:motif-sum}`, and in §4.4 Proposition 4.3 citing Lemma 4.2 (`\ref{lem:banded-certificate}`), would tighten internal signposting.
5. **Supplement notation mismatch warning for §6:** While the paper cleanly distinguishes $`N_{\ge k}(n)`$ from $`N_{=k}(n)`$, the cited supplement `proofs/degree-tail.md` writes $`N_k(n)`$ for the cumulative count $`N_{\ge k}(n)`$. Adding a parenthetical sentence in §6 prevents confusion for readers inspecting the supplement.
6. **Verifier hardening (`verify/motif-bracket/cert_verify_gmp.c` and `verify/run_all.sh`):** Add `if (H.D0 > H.KF + 1)` to the structural header checks in `cert_verify_gmp.c` (line 284) so the checker is sound on arbitrary binary inputs without relying on manual inspection of `Domains`, and add `-fsanitize=undefined -fno-sanitize-recover=all` to the `cert_verify_gmp.c` compile line in `verify/run_all.sh` (line 49) to match the banded verifiers.
7. **Section 6 scope for journal publication:** Section 6 states three substantial theorems (almost-sure lower transfer, polynomial-scale lower bound, and cumulative degree tail) whose proofs appear only in the external supplement. For *AAP/EJP*, the authors should either frame §6 concisely as a discussion of complementary results in the archival supplement (without presenting external results as proved in the article) or incorporate self-contained proofs of any results meant to be part of the published paper (in particular, Theorem D on the cumulative degree tail is short and self-contained enough to include as an appendix if desired).

---

## 7. Recommendation

**Recommendation: Minor revision (`minor`).**

This manuscript resolves two longstanding open questions on the random friend tree / isotropic redirection model—existence and almost-sure convergence of the non-leaf growth exponent $`\mu`$, and rigorous separation of $`\mu`$ from $`17/30`$—with clean, creative analytic arguments and exceptionally well-engineered, independently corroborated exact-arithmetic certificates. Both headline constructions pass full line-by-line proof auditing and clean UBSan exact-arithmetic replay. Addressing the three top revision points above (making the §2 coupling/re-driving details fully self-contained in the main text, updating the bibliography and archiving the supplement at a permanent DOI synchronized with the manuscript, and applying the minor notational/verifier/disclosure cleanups) will bring the manuscript fully up to publication standard for *Annals of Applied Probability* or *Electronic Journal of Probability*.
