CORRECT

# Cold Referee Review of `round2/finite/finite-operator.md`

- **Target reviewed**: `round2/finite/finite-operator.md` (`finite/finite-operator.md`) (118 lines, Sections 1–6, eqs. (13)–(16); SHA-256 `7a1dc1d9e2e28c4d0aca56a3d2ac6b8c885d065299e9b912af9d06e10f5b7967`).
- **Certificate cross-checked**: `round2/outputs/finite-matrices.json` (`outputs/finite-matrices.json`) (SHA-256 `5cbb167d459dff83ec2004cfe72910485944e6e42073e51f9b89f4eb9099d7dd`).
- **Upstream theorems referenced**: `round2/finite/one-port-upper.md` (`finite/one-port-upper.md`) (separately audited as `CORRECT` in [`round2/checks/one-port-neon.md`](one-port-neon.md)), `proofs/existence.md` (Proposition 6.2), and `proofs/bracket-motifs.md`.
- **Independent checker**: `round2/checks/finite_operator_neon.py` (`finite_operator_neon.py`) (stdlib-only, single CPU, independent center-peeling AHU canonicalization and exact `Fraction` generator construction, Tarjan SCCs, exact Collatz–Wielandt certification to $`<10^{-17}`$, exact Bareiss determinants, five matrix mutations, and direct Section 3–4 annotated branching & coupling simulations).

## Verdict Summary

The auxiliary branching construction, finite rational Metzler operator, spectral characterization (13), seam-potential accounting, one-window coupling bound (14), iterated squeeze (15)–(16), convergence $`\lambda_{M, M^2} \to \mu`$ at rate $`O(\log\log\log M / \log\log M)`$, and all six published exact matrices, algebraic polynomials, and 14-digit certified brackets in Section 6 are mathematically sound and verified without error or gap.

---

## Section-by-Section Independent Re-Derivation

### 1. The Auxiliary Branching Law and Centroid Cut Existence (Section 1)
- For $`M \ge 3`$ and $`B \ge 1`$, $`s = \lceil B/(M-1) \rceil \in [1, B]`$ (since $`M - 1 \ge 2`$).
- **State invariants**: Every reachable active core satisfies $`2 \le D_u < M`$, $`L_u \ge 0`$, and $`D_u - L_u - a_u \ge 0`$ (verified inductively across founder $`(2, 2, a=0)`$, hub newborns $`(2, 1, a=0)`$, promotion newborns $`(2, 1, a=1)`$, gains, hub replacements, and size cuts).
- **Legal cut existence**: Prior to a promotion at $`u`$, $`L_u \ge 1`$ implies $`a_u \le D_u - L_u \le M - 2`$, so after attaching the newborn core $`(2, 1)`$ every vertex of the oversized $`(B+1)`$-core tree still has active degree $`\le M - 1`$. Taking a centroid $`v_c`$ (all branches $`\le (B+1)/2`$), pigeonhole across its $`\le M - 1`$ branches gives a largest branch $`A_1`$ with $`|A_1| \ge \lceil B/(M-1) \rceil = s`$, and the complement $`A_2`$ has $`|A_2| = (B+1) - |A_1| \ge (B+1)/2 \ge |A_1| \ge s`$. Since $`s \ge 1`$, both sides lie in $`[s, B + 1 - s] \subseteq [1, B]`$.
- **Core birth rate**: Each hub births 1 core at rate $`1`$, each active core promotes at rate $`L_u / D_u \le 1`$, and gains/hub replacements/cuts preserve total cores $`\widetilde C`$. Hence $`\widetilde C`$ has birth rate $`\le \widetilde C`$ and $`\mathbb E \widetilde C(t) \le \widetilde C(0) e^t`$.

### 2. Finite Rational Operator and Accessible Classes (Section 2, eq. (13))
- For the finite-type multi-type branching generator $`G_{ij} = \sum_e r_{i, e}(\mathbb E_e \nu_j - \mathbf 1_{\{i=j\}})`$, entries are rational and off-diagonal entries are nonnegative (Metzler matrix), with $`G_{ij} > 0`$ ($`i \ne j`$) iff some positive-rate event from type $`i`$ can produce type $`j`$. Thus reachability in the branching law coincides with reachability in the off-diagonal directed graph of $`G`$.
- Conditioning on the first event gives $`\frac{d}{dt} M(t) = G M(t)`$, so $`\mathbb E_i \widetilde C(t) = (e^{tG} c)_i`$ with $`c_i \ge 1`$.
- Upper bound $`\limsup t^{-1}\log(e^{tG} c)_0 \le \lambda_{M, B}`$ follows from block-triangularizing $`G`$ by reachable SCCs and bounding the Jordan blocks of $`e^{tG}`$.
- Lower bound $`\liminf t^{-1}\log(e^{tG} c)_0 \ge \lambda_{M, B}`$ follows because $`(e^G)_{0, i} > 0`$ for every reachable $`i \in S^*`$ in a maximizing SCC $`S^*`$, dropping cross-class nonnegative off-diagonal entries only reduces $`e^{sG}`$ entrywise, and $`c_{S^*} \ge \varepsilon v_{S^*}`$ for the strictly positive Perron eigenvector $`v_{S^*}`$ of the irreducible Metzler block $`G_{S^* S^*}`$ (or trivially for a singleton SCC).
- Each $`\lambda_{M, B}`$ is a real algebraic number (dominant root of $`\det(\lambda I - G_{S^* S^*}) = 0`$). Dropping unreachable types is strictly required: e.g. at $`(M, B) = (3, 3)`$, the syntactically valid 3-core path `(2,0)-(2,0)-(2,0)` is unreachable from the founder (`9` reachable vs `10` valid types).

### 3. Seam Counting Potential (Section 3)
- For $`\Phi = sJ + P - (\widetilde C - 1)`$ with $`P = \sum_{\text{active } A} (|A| - s)_+`$:
  - Founding singleton $`T_3`$: $`J = P = 0`$, $`\widetilde C = 1 \implies \Phi = 0`$.
  - Hub birth: $`\Delta \widetilde C = +1`$, $`\Delta J = 0`$, $`\Delta P = (1 - s)_+ = 0 \implies \Delta \Phi = -1 \le 0`$.
  - Gain not hitting $`M`$: $`\Delta \Phi = 0`$.
  - Hub replacement (component of size $`n`$ splits into $`n_1, \dots, n_r \ge 1`$ with $`\sum n_i = n - 1`$): $`\Delta \widetilde C = \Delta J = 0`$, and superadditivity $`\sum (n_i - s)_+ \le (\sum n_i - s)_+ = (n - 1 - s)_+ \le (n - s)_+`$ gives $`\Delta P \le 0 \implies \Delta \Phi \le 0`$.
  - Promotion without cut: $`\Delta \widetilde C = +1`$, $`\Delta J = 0`$, $`\Delta P \in \{0, 1\} \implies \Delta \Phi \le 0`$.
  - Promotion with size cut ($`B \to B+1 \to (n_1, n_2)`$ with $`n_1, n_2 \ge s`$, noting $`B \ge s`$): $`\Delta \widetilde C = +1`$, $`\Delta J = +1`$, $`\Delta P = (n_1 - s + n_2 - s) - (B - s) = 1 - s \implies \Delta \Phi = s + (1 - s) - 1 = 0`$.
- Hence $`sJ + P \le \widetilde C - 1`$ (so $`J \le \widetilde C / s`$) holds pathwise at every time along every reachable annotated history and any continuation thereof.

### 4. Uniform Comparison on One Time Interval (Section 4, eq. (14))
- **True lift**: By induction on the annotated auxiliary history, $`D_u - L_u - a_u`$ equals the exact number of inactive core-tree edges incident to nonhub $`u`$. Restoring all inactive edges (both hub edges and seams) and attaching $`L_u`$ leaves to nonhubs and $`\max(0, M - c_u)`$ leaves to hubs produces a genuine finite tree with exact core set of size $`C_0 = \widetilde C(0)`$, matching $`(D_u, L_u)`$ at every nonhub and $`D_u \ge M`$ at every hub.
- **Coupling marginals and dirty-count drift**: Every event class was verified individually:
  - *Clean low pair*: Synchronized leaf births (rate $`L_u`$) and common rate-1 core ring match both true and auxiliary marginals (inactive-slot choices are discarded as null aux events). Active sends to clean low neighbors update both; sends to clean hubs leave them clean ($`D_{\mathrm{true}} \ge M`$ preserved); sends to dirty neighbors need no new marks (even if a dirty active neighbor turned into a hub in aux); seam sends to clean low neighbors mark both copies of the recipient dirty ($`+2`$ to $`I`$, total rate across both orientations of all $`J`$ seams at most $`2 \cdot \frac{1}{2} J = J`$, giving drift $`\le 2J \le 2\widetilde C/s`$).
  - *Clean hub pair*: Common rate-1 ring synchronized between true core ring and aux hub birth; true leaf births need no aux action. True leaf choice creates a paired clean newborn; true core-neighbor choice creates an unmatched aux newborn, sends a true leaf to the neighbor, and marks source pair ($`+2`$), newborn ($`+1`$), and clean recipient pair ($`\le +2`$), adding at most $`5`$ to $`I`$. By handshaking on the connected core tree, $`\sum_{\text{clean hub}} c_u / D_{u, \mathrm{true}} \le 2(C_{\mathrm{true}} - 1)/M \le 2 C_{\mathrm{true}}/M`$, giving drift $`\le 10 C_{\mathrm{true}}/M`$.
  - *Dirty vertices*: Individual true/aux rates preserve marginals; each dirty true core has core ring rate $`1`$ and each dirty aux core has promotion+active-send rate $`\le 1`$ (or hub birth rate $`1`$), adding at most $`+2`$ to $`I`$ per event, giving drift $`\le 2(I_{\mathrm{true}} + I_{\mathrm{aux}}) = 2I`$. Unmatched births always leave their parent dirty, so no clean vertex ever acquires an unmatched neighbor. Size cuts alter only inactive-slot labels and preserve all clean invariants.
- **Grönwall**: Integrating $`\mathbb E I(t) \le \int_0^t \big(2 \mathbb E I(w) + C_0(10/M + 2/s) e^w\big)\,dw`$ with $`I(0) = 0`$ yields $`\mathbb E I(t) \le C_0(10/M + 2/s)(e^{2t} - e^t)`$. Since $`|C_{\mathrm{true}}(t) - \widetilde C(t)| = |I_{\mathrm{true}}(t) - I_{\mathrm{aux}}(t)| \le I(t)`$, eq. (14) holds uniformly across all reachable starting states, including existing seams.

### 5. Multi-Interval Squeeze and Convergence Rate (Section 5, eqs. (15)–(16))
- Conditioning on the full annotated state at times $`nt`$ ($`n = 0, 1, \dots`$) and applying (14) against the uniform bounds $`m(t) \le \mathbb E_{T_{\mathrm{lift}}} C(t)/C(T_{\mathrm{lift}}) \le u(t)`$ iterates to $`(m(t) - d_{M, B}(t))^N \le \mathbb E \widetilde C(Nt) \le (u(t) + d_{M, B}(t))^N`$, yielding (15) as $`N \to \infty`$.
- **Uniform finite replacements**:
  - Upper: $`u(t) \le U_j(t)`$ by `one-port-upper.md` eq. (10).
  - Lower: `proofs/existence.md` Proposition 6.2 (line 330) holds uniformly for every finite tree $`T`$ with $`C(T) \ge 1`$ (including $`C(T) = 1`$, where the star center satisfies (E1) and $`|\mathcal E| = 1 \ge C/4`$). For $`r = t - h \ge 2jk`$, `one-port-upper.md` Section 9 verifies $`r \ge s_0`$ and $`F(r) - k - \theta_k(r) \ge F(r)/2`$, giving $`\mathbb E_T C(t)/C(T) \ge \frac{p_j}{4k} e^{-\eta_j r} \frac{F(r)}{2} = \frac{p_j}{8k} e^{-\eta_j r} F(r)`$, and trivially $`\ge 1`$. Hence $`m(t) \ge L_j(t)`$ and (16) follows whenever $`d < L_j(t)`$.
- **Asymptotic choice $`B = M^2`$, $`t = (\log M)/3`$**:
  - $`s = \lceil M^2/(M-1) \rceil \ge M + 1 > M`$, so $`d \le \frac{12}{M} e^{2t} = 12 M^{-1/3} \le 1/2`$ for $`M \ge 24^3`$ (where $`t \ge \log 24 > 1`$).
  - Integrating $`(e^{-u} F(u))' \le 2 e^{-3u}`$ over $`[r, t]`$ yields $`F(t) \le e^h(F(r) + \frac{2}{3} e^{-2r}) \le e^h(F(r) + 1) \le 2 e^h F(r)`$ (since $`F(r) \ge r \ge 8 > 1`$), hence $`L_j(t) \ge \frac{p_j}{16k} e^{-h - \eta_j t} F(t)`$, matching lines 87–88 exactly.
  - With $`U_j, L_j \ge 1`$ and $`d \le 1/2`$, $`\frac{U_j + d}{L_j - d} \le \frac{1+1/2}{1-1/2}\frac{U_j}{L_j} = 3\frac{U_j}{L_j}`$, adding at most $`\log 3 / t`$.
  - Canceling $`\log F(t)/t`$ against eq. (11) of `one-port-upper.md` and taking $`j = \lfloor \log t / (4 \log\log t) \rfloor`$ yields $`|\lambda_{M, M^2} - \mu| = O(\log\log t / \log t) = O(\log\log\log M / \log\log M)`$.

### 6. Independent Matrix Recomputation, Certificate Check, and Mutations (Section 6)
- All six matrices `(3,1), (3,3), (3,9), (4,1), (4,2), (4,3)` were reconstructed from scratch using center-peeling AHU canonicalization and exact `Fraction` arithmetic, independently of `finite_matrix.py` and `verify_matrix.py`.
- Every count (`types`, `classes`, `maximizing class size`), the explicit $`4 \times 4`$ matrix for `(3,1)`, the explicit characteristic polynomials for `(3,1)`, `(3,3)`, and `(4,1)`, all six entries of `round2/outputs/finite-matrices.json` (row-by-row identical up to type permutation), all six characteristic polynomials (up to degree $`55`$, checked at $`|S|+1`$ integers via exact Bareiss determinants), opposite endpoint signs, positive rational interval derivatives $`P' > 0`$ on the reported brackets, and $`< 10^{-16}`$ Collatz–Wielandt spreads were independently verified.
- Five concrete matrix mutations (deduplicating isomorphic cut outcomes before averaging at `(4,3)`, dropping inactive slots on cuts/hub removals at `(4,2)`, giving inactive slots send rate $`1/M`$ at `(3,3)`, changing hub newborns from `(2,1)` to `(2,2)` at `(3,1)`, and perturbing a single diagonal entry of `(3,1)` by $`10^{-6}`$) all measurably shifted $`\lambda_{M, B}`$ and broke certification.
- The Section 6 remark on `proofs/bracket-motifs.md` is accurate: the certified interval is $`\mu \in [0.5650, 0.5665]`$, bracketed on opposite sides by $`\lambda_{4,1} \approx 0.55612 < 0.5650`$ and $`\lambda_{3,1} \approx 0.61803 > 0.5665`$.

---

## Independent Verification Run Output

Command executed from repository root (single CPU):
```bash
python3 round2/checks/finite_operator_neon.py
```

Observed output (exit code `0`):
```text
=== TEST A: Independent Recomputation of All 6 Exact Matrices ===
  (M=3, B=1): types=4, classes=3, win_sz=2, lam in [0.618033988749895, 0.618033988749895] (gap=5.19e-24)
  (M=3, B=3): types=9, classes=5, win_sz=5, lam in [0.678337373185183, 0.678337373185183] (gap=3.09e-19)
  (M=3, B=9): types=27, classes=11, win_sz=17, lam in [0.679280747969823, 0.679280747969823] (gap=3.34e-18)
  (M=4, B=1): types=8, classes=5, win_sz=4, lam in [0.556123256048971, 0.556123256048971] (gap=4.24e-19)
  (M=4, B=2): types=23, classes=8, win_sz=16, lam in [0.614850634814983, 0.614850634814983] (gap=8.39e-19)
  (M=4, B=3): types=68, classes=14, win_sz=55, lam in [0.627558021428355, 0.627558021428355] (gap=4.51e-18)
  Exact (3,1) 4x4 matrix matches lines 107-112 identically.
  Exact charpolys for (3,1), (3,3), (4,1) in line 114 verified by exact Bareiss det().
  All 6 matrices, charpolys, endpoint signs, P'>0, and <1e-16 witnesses in finite-matrices.json verified 100%!
  [PASS] Test A

=== TEST B: Matrix Mutations & Planted-Error Sensitivity ===
  Mutation 1 (unique-outcome cut avg) at (4,3): lam_mut=0.6275239305 != lam_true=0.6275580214
  Mutation 2 (drop inactive slots) at (4,2): lam_mut=0.6125469871 != lam_true=0.6148506348
  Mutation 3 (inactive slots send at 1/M) at (3,3): lam_mut=0.7594363999 != lam_true=0.6783373732
  Mutation 4 (hub newborn (2,2) not (2,1)) at (3,1): lam_mut=0.8136065026 != lam_true=0.6180339887
  Reachability check at (3,3): unreachable valid types = ['3:(0,0,0)'] (count=1, reachable=9/10)
  Mutation 5 (1e-6 diagonal perturbation on (3,1)): lam=0.61803471235678 exits bracket!
  [PASS] Test B

=== TEST C: Sections 3-5 Potential Invariant, Exact Coupling & Bounds ===
  Sec 3 & exact E[C_tilde(1.0)] at (3,3): MC=2.2925+/-0.0225 vs (e^{tG}c)_0=2.2531 (0 invariant violations)
  Sec 3 & exact E[C_tilde(1.0)] at (4,2): MC=2.1582+/-0.0195 vs (e^{tG}c)_0=2.1649 (0 invariant violations)
  Sec 4 coupling (M=4, B=3, s=1, C0=8, hubs=3, seams=1, t=0.5): |E C_aux - E C_true|/C0=0.1113 <= E[I]/C0=0.6442 <= d_MB=4.8130
  Planted coupling error (omit hub recipient mark) caught by invariant assertion!
  Sec 4 coupling (M=6, B=4, s=1, C0=8, hubs=3, seams=2, t=0.5): |E C_aux - E C_true|/C0=0.0674 <= E[I]/C0=0.4189 <= d_MB=3.9217
  Section 5 scaling (B=M^2, t=log(M)/3, M>=24^3 -> d<=1/2, loss<=log(3)/t) verified.
  [PASS] Test C

ALL FINITE-OPERATOR VERIFICATION, MUTATION & COUPLING SUITES PASSED.
```

## Minor Non-Blocking Note
- In Section 5 (line 90), "Combine this with (10)–(12)" refers to equations (10)–(12) of `round2/finite/one-port-upper.md` (with `finite-operator.md` continuing numbering at (13)–(16)). Everything matches seamlessly when read together.
- No changes were made to any author draft or certificate file.
