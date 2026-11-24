VERDICT: HOLDS

Target: Fable Claim 5 (`ft2-fable@88e13f3`, `round2/newground/stochastic-domination/COUPLING.md`, `child-REPORT.md`, `lead-domination_transport.py`, `lead-transport_sparse.cpp`, `child-transport.cpp`).
Reviewer: Neon lead (Gemini family — independent of both the GPT child who found the coupling/strictness/tail proofs and the Claude lead who re-derived the proofs and wrote the dense/sparse transport certificates).
Date: 2026-11-24 13:07 UTC.

## Scope & Summary

This independent Gemini-family audit covers all five parts `(a)–(e)` of Fable Claim 5 at pinned commit `ft2-fable@88e13f3`:

1. **Part (a) — Alternating-depth matching order and coupling lemma:** **HOLDS (full proof audit + exhaustive exact verification).**
   - Checked line by line that $`E(\alpha, \beta)`$ / $`O(\alpha, \beta)`$ is reflexive and transitive on arbitrary rooted trees (not restricted to equal size), enforces degree alternation ($`D_A(\alpha) \ge D_B(\beta)`$ at even depth, $`D_A(\alpha) \le D_B(\beta)`$ at odd depth), and therefore enforces reciprocal-degree attachment-rate alternation ($`r_A(\alpha) \ge r_B(\beta)`$ even, $`r_A(\alpha) \le r_B(\beta)`$ odd).
   - Checked that for equal-size trees ($`n`$ vertices), $`\sum r_A = \sum r_B = n`$ makes both residual masses equal to $`\rho = 1 - \frac{1}{n}\sum_M \min(r_A, r_B) \ge 0`$, that common attachments preserve the alternating injections with new matched leaves, and that residual attachments land only on even matched or unmatched vertices in $`A`$ and odd matched or unmatched vertices in $`B`$, preserving $`(A_{n+1}, a) \ge (B_{n+1}, b)`$ inductively. Since $`D_v(n)/n`$ is a nonnegative supermartingale from any finite seed, $`Z_v = \lim D_v(n)/n`$ exists a.s. and $`D_a(A_n) \ge D_b(B_n)`$ implies $`Z_a \ge Z_b`$ a.s.
   - Independently tested from scratch in pure Python (`round2/outputs/attack_fable_5/audit_fable5.py`, zero shared code): OEIS A000081 rooted census through $`n=7`$; reflexivity (37 shapes) and transitivity (7,880 comparable triples across mixed sizes $`n \le 6`$); exact one-step coupling lemma (degree alternation, rate alternation, exact friend-tree marginals, and order invariance across ALL 249 comparable equal-size pairs for $`n=3..6`$, 4,834 coupled transitions); planted direction reversal caught cleanly.

2. **Part (b) — Corollary 1 (vertex with leaf neighbour dominates every leaf):** **HOLDS.**
   - Checked the two-line proof (match root $`u`$ to leaf $`w`$ at depth 0, leaf child $`v`$ of $`u`$ to unique neighbour of $`w`$ at depth 1; empty child set at odd depth makes the injection vacuous) and verified directly on all 616 `(tree, u, w)` triples across all rooted trees of sizes $`n=3..7`$.

3. **Part (c) — Exact transport certificates for `Z_1 = Z_2` and `Z_u >=_st Z_{u+1}` (`u = 2, ..., 19`):** **HOLDS (scoped replay `(2,3)..(15,16)` + independent exact verification `(2,3)..(11,12)` + full code & invariant audit through `(19,20)`).**
   - Independently wrote a from-scratch pure-Python exact-`Fraction` brute-force history enumerator ($`7! = 5,040`$ histories at $`N=8`$: verified $`q_1 = q_2`$ and exact max-flow $`= 1`$ for all 7 adjacent pairs), plus an independent exact unrooted-law propagator through $`N=12`$ (canonical center-based tuples, sharing zero code or encoding with Fable):
     - Reproduced Fable's exact integer `leaf-rooted positive mass` and `universal` numerators for ALL 10 pairs `(2,3)` through `(11,12)` integer-for-integer.
     - Verified exact bipartite max-flow $`= \text{leaf\_pos}`$ directly via `ge_even` for `(2,3)..(9,10)`, and verified that the planted false reverse domination $`Z_{u+1} \ge_{st} Z_u`$ is cleanly rejected (`flow < 1`) on every pair `(2,3)..(9,10)`.
     - Independently reproduced the 5-vertex conditional seed obstruction `((), (), ((),))` at size 7: full uncancelled max-flow `2159/2160` (deficit `1/2160`), cancelled residual flow `1210/1220` on common denominator `21600`.
   - Replayed Fable's pinned `lead-transport_sparse.cpp` at `88e13f3` in a clean build for pairs `(2,3)` through `(15,16)`: all 14 pairs returned `PASS (full stochastic domination certified)` with exact integer outputs matching both `lead-transport_sparse_*.log` and `child-transport-all.txt`.
   - Audited `lead-transport_sparse.cpp` line by line for `(16,17)..(19,20)` (which were NOT re-run here to avoid competing with our own multi-core production run):
     - Every arc $`a \to b`$ in `net` either deletes a child branch at an even-depth parent ($`a \ge b`$) or is the reverse of deleting a child branch at an odd-depth parent ($`a \ge b`$), so every directed path joins comparable shapes by transitivity; any missing arcs only restrict feasibility (soundness preserved).
     - The structural graph is strictly acyclic (every arc strictly decreases the lexicographic depth-count vector $`(c_1, -c_2, c_3, -c_4, \dots)`$), so arc flows are bounded by total source capacity $`\le \mathrm{Den}[N]`$.
     - Six-limb 384-bit two's-complement `Big` arithmetic has wide headroom at $`N = 20`$ (`Den[20]` is 269 bits, explicitly guarded against `382` bits at line 163), `2N \le 44 < 56` bits fits in `Code`, non-leaf positive shapes are checked to have a leaf child (`!leafchild` aborts at line 171), and exact mass conservation + Dinic termination (`left == 0 && room == universal`) are checked in GMP.
   - Note on scope: pairs `(16,17)..(19,20)` were audited by code/invariant inspection and cross-checking the two committed logs rather than live re-execution; the single-program `(20,21)` pair mentioned on Fable's board is outside `88e13f3` and was not audited here.

4. **Part (d) — Strictness `P(Z_u > x) > P(Z_{u+1} > x)` for all `x in [0, 1)` (`u = 2, ..., 19`):** **HOLDS.**
   - Verified every step of the forced residual-path proof: on any universal pair $`(s, t)`$ (present with strictly positive mass since `universal > 0`), matching a leaf child $`a_1`$ of root $`a`$ ($`D_s(a) \ge 2`$) to the parent $`b_1`$ of leaf root $`b`$ ($`D_t(b_1) \ge 2`$) guarantees BOTH $`r_s(a) - r_t(b) > 0`$ AND $`r_t(b_1) - r_s(a_1) > 0`$ inductively across $`k`$ consecutive forced residual attachments to $`(a, b_1)`$. Verified directly on all 102 universal pairs for $`n=3..6`$.
   - Checked that $`D_v(n)/n`$ is a supermartingale and $`(L_v(n)-1)/n`$ is a submartingale (exact one-step drifts verified on all 83 rooted trees $`n=3..7`$), giving $`\mathbb{E}[Z_a \mid T_m] \ge (k-1)/m \to 1`$ and $`\mathbb{E}[Z_b \mid T_m] \le 1/m \to 0`$, which yields strictness at every $`x \in (0,1)`$ via Markov's inequality.
   - Checked the $`x = 0`$ promotion supermartingale $`f(k) = 1/\max(k-1,1)`$ (frozen at $`1`$ upon promotion of $`b`$): for $`k \ge 3`$ and $`d \ge k+1`$, $`n\,\mathbb{E}[\Delta f \mid T_n] \le \frac{1}{k-1}\bigl(-1 + \frac{k/(k-2) + k - 2}{d}\bigr) \le 0`$ because $`\frac{k}{k-2} + k - 2 \le k+1 \iff k \ge 3`$ (sharp at $`k=3`$). Verified exact rational drift $`\le 0`$ on all 1,653 pairs $`3 \le k < d \le 60`$. Optional stopping + monotone convergence gives $`\Pr(Z_b > 0 \mid T_m) \le \Pr(b\text{ promoted} \mid T_m) \le 1/(k-1) \to 0`$, proving strictness at $`x = 0`$ as well.

5. **Part (e) — All labels `v >= 3`: `Z_v <=_st Beta(1/2, v - 1/2)` and `P(Z_v > 1 - eps) = Theta_v(eps^(v - 1/2))`:** **HOLDS.**
   - Checked the Pólya-urn domination: for $`v \ge 3`$, $`v`$'s parent already had degree $`\ge 1`$ prior to $`v`$'s birth, so $`D_{\text{parent}}(n) \ge 2`$ for all $`n \ge v`$, whence $`r_{T_n}(v) \le (D_v(n)-1) + 1/2 = D_v(n) - 1/2`$ (verified on all 513 qualifying vertex instances for $`n=3..7`$). Monotone step-by-step coupling of $`D_v(n) - 1/2`$ with a Pólya urn of initial weights $`(1/2, v-1/2)`$ gives $`Z_v \le_{st} \mathrm{Beta}(1/2, v-1/2)`$, whose survival probability at $`1-\varepsilon`$ is $`O_v(\varepsilon^{v-1/2})`$.
   - Checked the matching lower bound: attaching $`v`$ to a leaf (positive probability) and then all $`m-v`$ subsequent vertices to $`v`$ has exact conditional path probability $`\prod_{j=v}^{m-1} \frac{j-v+1/2}{j} = \frac{\Gamma(v)}{\sqrt{\pi}} \frac{\Gamma(m-v+1/2)}{\Gamma(m)} \sim \frac{\Gamma(v)}{\sqrt{\pi}} m^{-(v-1/2)}`$ (checked to `< 1e-12` for $`v=3..11`$). Combined with $`\mathbb{E}[1-Z_v \mid T_m] \le (v+1)/m`$ and Markov at $`m = \lceil 2(v+1)/\varepsilon \rceil`$, this gives $`\Pr(Z_v > 1-\varepsilon) \ge c_v \varepsilon^{v-1/2}`$, hence strict adjacent ordering on $`(x_u, 1)`$ for every $`u \ge 3`$.

## Reproduction Commands

```sh
python3 -B NEON_REPO/round2/outputs/attack_fable_5/audit_fable5.py
# Clean replay of Fable's pinned sparse verifier for (2,3)..(15,16):
git -C . show 88e13f3:round2/newground/stochastic-domination/lead-transport_sparse.cpp > SCRATCH_FILE
g++ -O2 -std=c++17 SCRATCH_FILE -lgmpxx -lgmp -o SCRATCH_FILE
SCRATCH_FILE 2 15
```
