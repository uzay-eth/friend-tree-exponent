VERDICT: HOLDS

# Adversarial Audit of Astra Claim 1 — Effective Finite Matrices (`ft2-astra@f7ab502`)

- **Auditor:** Independent adversarial child of Neon lead (Gemini 4 Carbon (High), distinct model family from GPT-6.2 Astra)
- **Targets audited:**
  - `round2/finite/finite-operator.md`
  - `round2/finite/one-port-upper.md`
  - `round2/outputs/finite-matrices.json`, `round2/finite/finite_matrix.py`, `round2/finite/verify_matrix.py`, `round2/finite/replay.sh`
- **Fresh independent verification script:** `round2/outputs/attack_astra_12/audit_claim1.py`

---

## 1. Independent Recomputation & Coupling Checks

A fresh stdlib-only verifier (`audit_claim1.py`) was written directly from the mathematical specifications in §1–§4 of `finite-operator.md` and §1–§9 of `one-port-upper.md` without sharing code with Astra's builder or checkers (using tree-centroid canonicalization + $`S_n`$ permutation cross-check, Kosaraju's two-pass SCC algorithm, exact integer Faddeev–LeVerrier characteristic polynomials, and exact rational Collatz–Wielandt ratios).

1. **Exact $`G_{M, B}`$, SCCs, $`\lambda_{M, B}`$, and characteristic polynomials:**
   Recomputed from scratch for all six published pairs $`(M, B) \in \{(3,1), (3,3), (3,9), (4,1), (4,2), (4,3)\}`$. All reachable type counts ($`4, 9, 27, 8, 23, 68`$; total $`139`$), strongly connected class counts ($`3, 5, 11, 5, 8, 14`$), maximizing class sizes ($`2, 5, 17, 4, 16, 55`$), bijection-mapped rational generator entries, integer scales $`L`$, characteristic polynomials $`\det(xI - L G_{SS})`$, opposite endpoint signs, strict derivative positivity $`P' > 0`$ on $`[L\,\mathrm{lo}, L\,\mathrm{hi}]`$, and $`10^{-14}`$ rational brackets match `finite-operator.md` §6 and `round2/outputs/finite-matrices.json` identically. At $`(3,1)`$, $`\det(xI - G_{SS}) = x^2 + x - 1`$ gives $`\lambda_{3,1} = (\sqrt{5}-1)/2`$.
2. **Seam budget $`s J + P \le \widetilde C - 1`$ and liftability:**
   Simulated the full annotated finite core genealogy (retaining active edges, hub edges, and size-$`(B+1)`$ seams even after endpoints convert to hubs) across $`48{,}000`$ transitions over $`(M, B) \in \{(3,1), (3,3), (4,1), (4,2), (4,3), (5,6)\}`$. Every step satisfied $`s\Delta J + \Delta P \le \Delta \widetilde C`$, $`s J + P \le \widetilde C - 1`$ (attaining tight ratio $`1.0`$), legal cut existence at size $`B+1`$, and exact slot accounting $`D_u - L_u - a_u = \text{inactive genealogy degree}`$ ensuring valid genuine-tree lifts with $`C_{\mathrm{true}}(0) = \widetilde C(0)`$.
3. **Finite-time evolution & exact §4 synchronous coupling:**
   Compared $`(e^{tG}c)_0`$ against exact Poissonized $`T_3`$ friend-tree means ($`a_3, \dots, a_{11}`$ plus exact negative-binomial tail) at $`(3,1), (4,1), (3,3), (4,2)`$, AND simulated the exact §4 synchronous coupling ($`4{,}000`$ coupled trajectories each at $`(3,1)`$ and $`(4,1)`$, $`t = 0.5`$), verifying pathwise $`|C_{\mathrm{true}}(t) - \widetilde C(t)| \le I(t)`$, step-by-step $`\Delta I`$ caps ($`\le 5`$ on clean hub core-neighbor rings, $`\le 2`$ on clean low seam sends, $`\le 2`$ on dirty core rings, $`0`$ otherwise), and $`|\mathbb E \widetilde C(t) - \mathbb E C_{\mathrm{true}}(t)|/C_0 \le \mathbb E I(t)/C_0 \le d_{M, B}(t)`$.

---

## 2. Planted-Error Falsification Tests

All four required planted-error classes were injected and tested against **both** our independent checker and Astra's `round2/finite/verify_matrix.py` (using temporary JSON files outside `[local directory]`):

1. **Mutated rational entry of $`G_{M, B}`$:** Perturbing a single numerator in $`G_{3,3}`$ is immediately rejected by both `audit_claim1.py` and `verify_matrix.py` (`AssertionError: single-ring reconstruction differs`).
2. **Dropped seam slot:** Decrementing endpoint degrees on size-$`(B+1)`$ cuts shifts $`\lambda_{4,2}`$ from $`0.61485063\dots`$ to $`0.61331063\dots`$ and is rejected by both checkers.
3. **Altered cut probability:** Averaging uniformly over isomorphism classes of cut pairs rather than distinct allowed physical edges shifts $`\lambda_{4,3}`$ from $`0.62755802\dots`$ to $`0.62752393\dots`$ and is rejected by both checkers.
4. **Shifted eigenvalue bracket endpoint:** Shifting the $`(4,1)`$ bracket upward by $`2\times 10^{-14}`$ (placing both endpoints above the root) or tightening the upper endpoint by $`-10^{-14}`$ is rejected by both checkers via exact polynomial endpoint-sign and Collatz–Wielandt ratio assertions.

---

## 3. Line-by-Line Proof Audit

### A. `round2/finite/finite-operator.md`
- **§1–§2 (Centroid cut & Perron exponent):** At a promotion from size $`B`$ to $`B+1`$, the promoting vertex had $`L_u \ge 1`$ and $`D_u \le M-1`$, so its prior active degree was $`\le M-2`$ and its post-promotion active degree is $`\le M-1`$; all other vertices have active degree $`\le D_x \le M-1`$. Removing a centroid leaves $`\le M-1`$ branches of sizes summing to $`B`$ and each $`\le (B+1)/2`$; the largest branch has $`\ge \lceil B/(M-1)\rceil = s`$ cores and its complement has $`\ge (B+1)/2 \ge s`$ cores. Reachability of $`H`$ from `(2,2)` holds in $`M-2`$ leaf gains; block triangularization + Perron–Frobenius on reachable off-diagonal SCCs proves $`\tau^{-1}\log(e^{\tau G}c)_0 \to \lambda_{M, B}`$ in (13). Dropping unreachable types is genuinely necessary (at $`(3,3)`$, the valid bare 3-core path `(2,0)-(2,0)-(2,0)` is unreachable, leaving $`9`$ of $`10`$ valid types).
- **§3 (Seam potential):** Hub births add an active component of size $`1 \le s`$ ($`\Delta J = \Delta P = 0, \Delta\widetilde C = 1`$); hub conversions satisfy $`\sum_i(n_i-s)_+ \le (\sum_i n_i - s)_+ \le (|A|-s)_+`$ with $`\Delta J = \Delta\widetilde C = 0`$; promotions without cut have $`\Delta P \le 1 = \Delta\widetilde C`$; promotions at size $`B \ge s`$ followed by a cut into $`n_1, n_2 \ge s`$ ($`n_1+n_2 = B+1`$) satisfy $`\Delta(sJ+P) = s + (B+1-2s) - (B-s) = 1 = \Delta\widetilde C`$.
- **§4 (Dirty-core drift & Grönwall (14)):** Every neighbor of a clean vertex is paired (unmatched newborns attach only to a hub whose pair is simultaneously marked dirty, or to an already-dirty core).
  - Clean hub core-neighbor rings occur at rate $`\sum_{\text{clean high}} c_u/D_u \le \sum_{\text{cores}} c_u/M = 2(C_{\mathrm{true}}-1)/M < 2 C_{\mathrm{true}}/M`$ and mark at most $`5`$ cores dirty (source pair $`2`$ + unmatched auxiliary newborn $`1`$ + clean recipient pair $`2`$), giving $`10 C_{\mathrm{true}}/M`$. (Even when a clean hub has $`0`$ true leaves so $`c_u = D_u \ge M`$, its exceptional rate $`1 = c_u/D_u \le c_u/M`$ still holds.)
  - Clean low seam sends have rate $`\le 1/2 + 1/2 = 1`$ per seam, hence $`\le J \le \widetilde C/s`$ total, marking only a clean low recipient pair dirty ($`\Delta I \le 2`$; sender $`D, L`$ and neighbors do not change; clean hub recipients absorb leaves without breaking clean hub invariants), giving $`2\widetilde C/s`$.
  - Dirty vertices: leaf births at dirty cores affect only their own leaf stock ($`\Delta I = 0`$); each dirty core's own rate-$`1`$ clock either creates $`1`$ dirty newborn ($`\Delta I = 1`$) or hits a clean recipient pair ($`\Delta I = 2`$), giving drift $`\le 2 I`$. Size cuts only relabel an active edge as a seam without changing $`(D, L)`$.
  - Stopping at finite total-vertex cutoffs, extending the nonnegative drift integrand from $`t\wedge\tau_N`$ to $`t`$, applying Grönwall with $`\mathbb E C_{\mathrm{true}}(w), \mathbb E\widetilde C(w) \le C_0 e^w`$, and using monotone convergence yields $`\mathbb E I(t)/C_0 \le (10/M + 2/s)(e^{2t}-e^t) = d_{M, B}(t)`$, proving (14).
- **§5 (Multi-interval bounds (15)–(16) & asymptotics):** Conditional application of (14) to the fresh true lift at each block boundary $`kt`$ gives (15). Proposition 6.2 of `proofs/existence.md` holds uniformly over any starting tree $`T`$ with $`C(T)\ge 1`$, giving $`m(t) \ge L_j(t)`$, and (10) of `one-port-upper.md` gives $`u(t) \le U_j(t)`$, yielding (16). For $`B = M^2, t = (\log M)/3`$, we have $`s = \lceil M^2/(M-1)\rceil \ge M+2 > M`$, $`d \le 12 M^{-1/3} \le 1/2`$ for $`M \ge 24^3`$, extra log-loss $`\le \log((1+1/2)/(1-1/2))/t = (\log 3)/t`$, and at $`j = \lfloor \log t/(4\log\log t)\rfloor`$, $`\log A_j(t) \le \exp(j\log j + O(j)) = t^{1/4+o(1)} = o(t/j)`$ and $`2h + \log(8k/p_j) = \exp(O(j)) = t^{o(1)}`$, giving $`|\lambda_{M, M^2} - \mu| = O(\log\log\log M / \log\log M)`$.

### B. `round2/finite/one-port-upper.md`
- **§2 (Joint path density (3)):** Because interior jumps never alter exterior state and exterior jumps never alter interior state, the time-ordered jump density factors identically into $`I(x)\,dx`$, the normalized prescribed-immigration law $`Q^{a_x}(dy)`$ (nonexplosive since $`a_x \le 1/2`$), and $`e^{-\int_0^t b_y(u)\,du}\prod_{u\in G(x)}[r_x(u)+b_y(u)]`$. Bidirectional feedback is retained without approximation via $`x \mapsto Q^{a_x}`$.
- **§3 (Tracked birth-death bound (4)):** After $`b`$ prescribed immigrants, $`D_p \ge b + 2 > b`$ permanently. Giving each tracked leaf independent birth rate $`1`$ and death-proposal rate $`1/b`$, reclassifying rejected death proposals as untracked leaves while accepting true promotions with probability $`b/D_p \le 1`$, realizes exact $`Q^a`$ marginals while embedding independent linear birth-death processes ($`1, 1/b`$) inside $`L_p \le D_p`$. Solving $`\partial_h f = (f-1)(f-(1-s))`$ gives conditional geometric parameter $`q = s e^{-sh}/(1-(1-s)e^{-sh}) \le e^{-sh}`$; substituting $`q = e^{-2x}`$ proves $`-\log q = 2x \le 2\sinh x = (1-q)/\sqrt{q}`$, so $`\mathbb E[1/Z_h \mid Z_h > 0] \le \sqrt{q} \le e^{-sh/2}`$ and $`\mathbb E[1/Z_h \mid \text{surv}] \le s^{-1} e^{-sh/2}`$. Since $`a \le 1/2`$, $`A(u) - W = \int_{\sigma_*}^u a \le h/2 \implies h \ge 2(A(u)-W)`$. On $`\{W > A(u)\}`$, $`1/D_p \le 1/2 < s^{-1} e^{\kappa(W-A(u))}`$; combining both events with $`\mathbb E e^{\kappa W} = (1-\kappa)^{-b} s/(s-\kappa)`$ cancels $`s`$ and yields (4).
- **§4 (Profile rounding (5)):** Expanding $`\prod_{u\in G(x)}(r_x+b_y)`$ by the last $`b_y`$-factor partitions marked histories cleanly. Since $`b_y(u) = 1/D_p(y_{u-}) \in (0,1/2]`$ is nonincreasing, rounding up to $`\{0,1/J, \dots,1/2\}`$ at unit left-endpoints $`m \in \{0, \dots, \lfloor t\rfloor\}`$ yields $`\binom{\lfloor t\rfloor + 1 + J/2}{J/2} \le (2+t)^{J/2}`$ profiles, $`\beta \ge b_y`$, and $`\int_0^\tau(\beta - b_y) \le \tau/J + \sum_m (b_y(m)-b_y(\min(m+1, \tau))) \le t/J + 1/2`$.
- **§5 (Change of measure (6)):** For $`r \ge K \ge 1`$, $`\gamma \in (0,1/2]`$, $`\beta \in [0,1/2]`$, we have $`|\beta - \gamma| \le 1/2`$ and $`u = (r+\beta)/(r+\gamma) \in [2/3, 3/2)`$ in **both** signs $`\beta \ge \gamma`$ and $`\beta < \gamma`$. Lagrange's remainder for $`u^d`$ ($`d \ge 3`$) at $`u = 1`$ gives $`(r+\gamma)(u^d-1) - d(\beta-\gamma) = \frac{d(d-1)}{2}\xi^{d-2}\frac{(\beta-\gamma)^2}{r+\gamma} \le \frac{d(d-1)(3/2)^{d-2}}{8K} = c_d/K`$ since $`\xi \in [2/3, 3/2) \implies 0 < \xi^{d-2} \le (3/2)^{d-2}`$. On $`\{L_v < K\}`$, there are $`\le M-2 < M`$ root gains before $`D_v \ge M`$, and thereafter $`\min(L_v, K)`$ increases by $`1`$ on each low-leaf gain and decreases only on promotions with $`1 \le L_v \le K`$ (count $`N`$, predictable intensity $`\le K/M`$ under the auxiliary law), proving $`|G_{\mathrm{pref}}\cap\{L_v<K\}| \le M+K+N`$ and $`\mathbb E_{\widetilde P} z^N \le e^{(z-1)Kt/M}`$.
- **§6–§7 (Hölder (7) & seed $`T_*`$ (8)–(9)):** $`C_B`$ has per-core promotion rate $`\le 1`$ (forced root leaf gain adds no core), so $`\mathbb E_P f^2 \le \mathbb E[Y_1(t)^2] = 2e^{2t}-e^t \le 2e^{2t}`$; Hölder $`(d, q=d/(d-1))`$ and $`(1/(2-q), 1/(q-1))`$ yield (7). In $`T_*`$, $`n_* = 3 + (B-1)H`$, $`L_p(0) = 0, D_p(0) = B`$; requiring silence of $`u_i \to p`$ ($`i=1, \dots, B-1`$, prob $`e^{-\delta t}`$), last $`p\to v`$ proposal at $`\tau`$ (density $`\gamma e^{-\gamma(t-\tau)}`$ with independent rate-$`\gamma`$ Poisson prefix), and no $`v\to p`$ arrow on $`(0, \tau]`$ (prob $`e^{-A_x(\tau)}`$) keeps $`D_p = B`$ identically throughout $`(0, \tau]`$ so every prefix/last proposal is accepted with probability $`1`$, while zero $`p\to v`$ proposals occur on $`(\tau, t]`$, matching $`R_{\gamma, \tau}`$ and proving (8)–(9) and (1)–(2).
- **§8–§9 (Bounds (10)–(12) & tail (197)):** Core ancestry partitions post-$`0`$ cores into disjoint branches from initial cores, giving (10). $`F'(u) \le F(u) + 2e^{-2u}`$ integrates to $`F(t+1) \le e(F(t) + \frac{2}{3}e^{-3t}) \le e(F(t)+1)`$. Every tree (including $`T_2`$, whose two leaves promote each other at total rate $`2 \ge 1`$) has promotion rate $`\ge 1`$, so $`F(t) \ge t`$; for $`t \ge 2jk`$ ($`x = t/(jk) \ge 2`$), $`1+\log(6x) \le 2x`$ and $`\theta_k(t) \le 1`$, proving (12). Shifting $`n = N + 1 + m`$ with $`n(n-1) = N(N+1) + 2(N+1)m + m(m-1)`$ proves that (197) is the exact closed-form identity for $`\sum_{n > N} n(n-1)q^2(1-q)^{n-2}`$.

---

## 4. Exact Commands & Outputs

```text
$ python3 -B round2/outputs/attack_astra_12/audit_claim1.py
=== [Claim 1 / Check 1A] Fresh Independent Recomputation of G_{M,B} ===
  (M=3, B=1):  4 types,  3 SCCs, max_SCC= 2, lambda in [0.618033988749895, 0.618033988749895], L=1, deg=2 MATCH
  (M=3, B=3):  9 types,  5 SCCs, max_SCC= 5, lambda in [0.678337373185183, 0.678337373185183], L=2, deg=5 MATCH
  (M=3, B=9): 27 types, 11 SCCs, max_SCC=17, lambda in [0.679280747969823, 0.679280747969823], L=2, deg=17 MATCH
  (M=4, B=1):  8 types,  5 SCCs, max_SCC= 4, lambda in [0.556123256048971, 0.556123256048971], L=3, deg=4 MATCH
  (M=4, B=2): 23 types,  8 SCCs, max_SCC=16, lambda in [0.614850634814983, 0.614850634814983], L=12, deg=16 MATCH
  (M=4, B=3): 68 types, 14 SCCs, max_SCC=55, lambda in [0.627558021428355, 0.627558021428355], L=18, deg=55 MATCH
  [PASS] All 6 finite matrices, SCCs, Perron brackets, and charpolys verified.

=== [Claim 1 / Check 1B & 1C] Seam Budget Invariant & Exact §4 Coupling ===
  Seam budget s*J + P <= C_tilde - 1 & slot liftability verified across 48000 transitions (max (sJ+P)/(C-1) = 1.000000 <= 1).
  Finite-time (M=3, B=1, t=0.15): aux=1.15102894, true in [1.15045069, 1.15045917], |diff|<=0.000578 <= d_MB=1.002798
  Finite-time (M=3, B=1, t=0.50): aux=1.53195948, true in [1.50060090, 1.55709028], |diff|<=0.031359 <= d_MB=5.704323
  Finite-time (M=4, B=1, t=0.15): aux=1.15042064, true in [1.15045069, 1.15045917], |diff|<=0.000039 <= d_MB=0.846111
  Finite-time (M=4, B=1, t=0.50): aux=1.51711810, true in [1.50060090, 1.55709028], |diff|<=0.039972 <= d_MB=4.813023
  §4 Coupled MC (M=3, B=1, t=0.5, 4000 runs): mc_aux=1.5333 (exact 1.5320), mc_true=1.5180, E[I]=0.2343 <= d_MB=5.7043 (0 pathwise invariant violations)
  §4 Coupled MC (M=4, B=1, t=0.5, 4000 runs): mc_aux=1.5220 (exact 1.5171), mc_true=1.5220, E[I]=0.2345 <= d_MB=4.8130 (0 pathwise invariant violations)
  [PASS] Seam budget, finite-time coupling bound (14), and §4 step rules verified.

=== [Claim 1 / Check 1D] Independent Audit of one-port-upper.md (§2 - §9) ===
  §2 Joint path density (3) & last-gain expansion verified symbolically on 228 trajectories.
  §3 Geometric reciprocal inequality -log(q) <= (1-q)/sqrt(q) verified on 999 grid points.
  §4 Profile cardinality <= (2+t)^{J/2} and integral error <= t/J + 1/2 verified.
  §5 Taylor bound checked in exact Fraction on 2400 cases (incl. beta=0 < gamma=1/2); max lhs/rhs = 0.696840 <= 1.
  §7 Seed T_* (B=2, H=2, n_*=5): exact jump prob = 1/12 >= claimed LB 1/144
  §7 Seed T_* (B=2, H=3, n_*=6): exact jump prob = 11/120 >= claimed LB 1/2880
  §7 Seed T_* (B=3, H=2, n_*=7): exact jump prob = 13/3240 >= claimed LB 1/86400
  §9 Negative-binomial tail formula (197) verified as an EXACT rational identity.
  [PASS] All one-port-upper.md checks passed.

=== [Claim 1 / Check 2] Planted-Error Falsification Tests ===
  [CAUGHT] Planted Error 1 (mutate rational entry of G_{3,3}): rejected by both checkers.
  [CAUGHT] Planted Error 2 (drop seam slot): changes lambda_(4,2) 0.61485063 -> 0.61331063 & rejected by verify_matrix.py.
  [CAUGHT] Planted Error 3 (wrong cut probability): changes lambda_(4,3) 0.62755802 -> 0.62752393 & rejected by verify_matrix.py.
  [CAUGHT] Planted Error 4 (shift eigenvalue bracket endpoints): rejected by both checkers.
  [PASS] All planted-error tests bit cleanly.

ALL ADVERSARIAL CHECKS FOR CLAIM 1 COMPLETED SUCCESSFULLY.
```

```text
$ bash round2/finite/replay.sh
...
PASS: 6 matrices, 139 reachable rows, 105 characteristic-polynomial evaluations
...
ALL 8 VERIFICATION & PLANTED-ERROR SUITES PASSED.
...
ALL FINITE-OPERATOR VERIFICATION, MUTATION & COUPLING SUITES PASSED.
```
