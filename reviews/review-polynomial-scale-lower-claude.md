VERDICT: CORRECT

# Cold referee audit: `AS_LOWER_PROOF.md` (almost-sure lower bound liminf C_n / n^0.543 > 0)

* **Referee:** Claude Opus 5.5 (`claude-opus-5-5`, Claude family). The proof was written by a GPT-family model.
* **Date:** 2026-11-22.
* **Read:** `attempts/strogatz-carbon/sol_as_lower/AS_LOWER_PROOF.md` in full, and `attempts/x-strogatz/REPORT.md`
  §1 (model) and §§3.1–3.5. I did **not** open `verify_as_lower.py`, the `certs*.py` scripts, or any other
  review.
* **Re-read:** later on 2026-11-22 I re-read the revised `AS_LOWER_PROOF.md`, which applies my wording fixes
  M1–M4, to record their status. See §0, item 2, and §8.
* **Own script:** `checks/strogatz-carbon-as-lower-claude/check_as_lower.py`, written from scratch. It uses
  exact `fractions` arithmetic and ran on one core (`taskset -c 0`). Its full output is in §9.

## 0. What the verdict covers

The proof is a deduction. **If** the two certificates of REPORT §3.5 hold, **then** liminf_n C_n/n^δ > 0 almost
surely, with δ = 543/1000. The certificates are:

* the lower table `certs9.py … 543/1000 60 10 2 6`, with δ = 0.543;
* the upper table `find_upper_lin.py 60 10 1 20 … 293/500`, with λ = 0.586.

I checked that deduction line by line, and I found **no mathematical error and no missing step that affects
the headline result**.

Two things are outside the verdict:

1. **The certificate tables were not re-verified.** The brief limits me to the proof and REPORT. I accepted
   REPORT §3.5's statements as hypotheses: the exact-rational checks (a)–(c), Φ_3 = 0.55952, sup f_- = 0.63907,
   sup f_+ = 4.8045, and f_+ ≥ 1. I did independently check, by brute force, the bridge from those certificates
   to the proof's hypothesis (1): REPORT's Lemma 2 event list at true rates, and the rate intervals (§9 A).
2. **The general form of the theorem was loosely worded in the version I first audited.** It covered "the
   rooted version" and "any shell that drops information". Read one natural way, that theorem was not proved
   as worded (issue M1). The δ = 0.543 corollary never depended on that wording: both of its tables are
   explicitly defined in REPORT and are radius-1 local, which is the property the proof needs.
   **Status after the revision:** the theorem now states radius-1 locality and degree-H saturation as explicit
   hypotheses, and with those hypotheses every step is correct. **M1–M4 are resolved.** The revision adds one
   new side remark about radius-R tables, which overstates slightly (nit N1 in §8). It does not touch the
   theorem as stated or the corollary.

The minor issues are listed in §8. None of them changes the conclusion.

---

## 1. The setup and the theorem statement

* The model matches REPORT §1: a uniform V_n, then a uniform neighbour W_n, then a leaf attached to W_n.
  C_n = |{deg ≥ 2}|. The tree is rooted at the centre r of T_3, which stays core forever.
* Creator = parent. A promoted leaf y is a non-root vertex, so its unique neighbour u (its creator) is its
  parent. Parent relations never change. This justifies "original creator labels and root throughout".
* Hypothesis (1) is a per-vertex drift: n·E[f(x_{u,n+1}) − f(x_{u,n}) + 1{W_n leaf of u} f(x_{y,n+1}) | F_n].
  This is exactly the per-vertex summand of REPORT (★★), provided the event list of Lemma 2 is complete and
  exact. Every leaf has exactly one core neighbour, so the newborn indicators are disjoint across u. I checked
  this identity by brute force (A: 5,408 vertex comparisons, 0 mismatches) for:
  * the certs9 state with (K0, D_cap) ∈ {(1,4), (2,6), (3,5), (2,4), (2,5)};
  * the upper state (L, c, m1, b, P) with P_max ∈ {3, 6, 20}.

  Three planted event-list errors were all detected: a wrong E4 rate, a misclassified newborn, and dropped E7
  events. So REPORT's comparison lemma ("min over the rate intervals of 𝒜f ≥ δf") does imply (1) ≥ δ f_- at
  the true rates. The true rates lie in REPORT's intervals in every case tested (credit, E3/E4 exact, E5, E7).
* **Exponent hypotheses.** The first version assumed 1/2 < δ ≤ λ < 2δ, but only q = 2δ − λ > 0 and δ > 0 are
  used (nit M2). The revision assumes 0 < δ ≤ λ < 2δ, so M2 is resolved. The condition δ ≤ λ is automatic:
  c·n^δ ≤ E Φ^- ≤ F·E C_n ≤ F·U·C·n^λ.
  For the application, δ = 0.543, λ = 0.586, 2δ = 1.086, q = 0.5. ✓
* The other hypotheses hold for the final tables (per REPORT): 0 ≤ f_- ≤ F = 0.63907, 1 ≤ f_+ ≤ U < ∞, and
  a_0 = Φ_3 = 0.55952 > 0. ✓

## 2. §1: receivers, jumps, second moment

* **Receiver law.** P(W = w | F_n) = n^{-1} Σ_{z∼w} 1/D_z. The receiving rate of a core vertex u is L_u + ρ_u.
  The rate at which u promotes one of its leaves is L_u/D_u. ✓ Exact on 170 random trees (A); the planted
  preferential-attachment law was rejected on 170/170 trees.
* **Promotion.** (L,c) → (L−1, c+1) at the same degree, and the newborn has (L,c) = (1,1). The promotion
  probability is S_n/n ≤ C_n/n. ✓
* **(2).** Σ_{u∈K} ρ_u = Σ_{v∈K} c_v/D_v ≤ C_n, obtained by exchanging the double sum. ✓ (A)
* **Changed summands for certs9.**
  * A gain at u changes u's own entry, plus the entries of neighbours that track u: one if u is a pendant, two
    if u is in class B, at most D_cap − 1 if u is a tracked "other" (since c_u ≤ d), and none if u is in `top`.
    Total at most max(3, D_cap).
  * A promotion changes at most 3 summands: the creator, the newborn, and at most one old neighbour, and only
    when the creator was a pendant. Class B is impossible because L > 0, and an "other" keeps its degree.
  * So B = F·max(3, D_cap). ✓ (A: J1, J2)
* **The L ≥ L_max claim.** Assume L_max > K0 and L_max + 1 ≥ D_cap. Then a gain at a core u with L_u ≥ L_max
  leaves the own entry fixed (flat in L, and c unchanged) and leaves every neighbour entry fixed (u is in `top`).
  ✓ (A: J3, 1,667 tests).
  Negative control: when L_max + 1 < D_cap, 152 violations were found, so the verifier's condition really is
  needed.
* **(3).** P(W ∈ K, L_W < L_max) ≤ Σ_{L_u<L_max}(L_u + ρ_u)/n ≤ L_max·C_n/n ≤ (L_max+1)·C_n/n. Adding the
  promotion probability C_n/n gives E[(ΔΦ^-)²] ≤ B²(L_max+2)C_n/n. ✓ (A: 850 exact tests; the largest
  LHS/RHS ratio observed was 0.04)
* **The uniform H.**
  * Take H ≥ the own-leaf cutoffs, all degree and label caps, and K0 + 2. A vertex of degree ≥ H is "other" in
    every classification: pendants have degree ≤ K0 + 1 and class B has degree 2. Its descriptor is saturated:
    `top`, or Q = Q_max, or a static P.
  * So only receivers or creators of degree < H can change neighbour entries, at most H − 1 of them. Including
    the receiver or creator and the newborn, at most H + 1 summands change, giving B_J = (H+1)J.
  * A core receiver with L ≥ H changes nothing.
  * ✓ (A: J4a–J4d, and J4U for the upper table.)
  * For the final pair of tables, H = max(60, 6, 4, 60, 20, 3) = 60 works. The first version never named a
    value (nit M4). The revision states "H = 60 suffices" in §1, so M4 is resolved.

## 3. §2: restriction across a saturated parent edge

* **Exact drift of the restricted sums.** R_n is the v-side of {p, v}, and I_n = K_n ∩ R_n. The newborn lies
  in I_{n+1} \ I_n exactly when its creator lies in I_n, because the leaves of p are on p's side and the leaves
  of inside vertices are inside. Summing (1) over u ∈ I_n therefore gives the exact drift of X and of Y. This
  identity holds even without saturation: (1) counts every change to an inside summand, whatever causes it.
  So (4) follows from the per-vertex certificates evaluated at full-tree states and true rates, and the 1/D_p
  term in ρ_v is covered by the certificate's interval.
  ✓ (B: 46,050 exact outcome-level checks; omitting the newborn breaks the identity on 6,420 outcomes.)
* **Saturation is what controls the jumps.** Let D_p ≥ H ≥ D_cap. Then an exterior receiver changes no inside
  state:
  * the only exterior vertex any inside vertex sees is p, and only v sees it;
  * p stays in `top`, or P = P_max;
  * a promotion at p keeps D_p fixed.
  ✓ (B: 33,073 exterior outcomes, including p and the leaves of p, checked against both the raw certs9 and the
  raw upper states.)
  Negative control: with an unsaturated p (D_p < D_cap), 231 exterior outcomes changed an inside state.
* **(5).** Σ_{u∈I} ρ_u = Σ_{z∈I} |N(z) ∩ I|/D_z + 1/D_p ≤ N + 1 ≤ 2N, using N ≥ 1 because v is core. ✓ (B)
* **(6).** The inside change probability is at most [(H−1)N + 2N + N]/n = (H+2)N/n ≤ (H+3)N/n, and the jump is
  at most B_F, so E[(ΔX)²] ≤ (H+3)B_F² N/n. The same holds for Y with B_U. ✓ (B: 1,076 exact checks; the
  largest ratio was 0.010.)
* Remark: Y's second moment is never needed for the main theorem. Only Y's restricted drift is used, in (8),
  and that drift needs no saturation.

## 4. §3: normalized martingale and uniform survival probability

* **(7).**
  * Upper bound: log(1+x) ≤ x and Σ_{j=m}^{n−1} 1/j ≤ log(n/m) + 1/m, so the upper constant can be e^{a/m} ≤ e^a.
  * Lower bound: log(1+x) ≥ x − x²/2 and Σ_{j≥m} j^{-2} ≤ 1/(m−1) ≤ 1/2, so the lower constant can be
    e^{−a²/4} ≥ e^{−a²}.
  * The discrepancy |log(1+a/j) − a·log(1+1/j)| ≤ (a² + a)/(2j²) is summable, so P_a(m,n)/n^a has a positive
    finite limit.
  * ✓ (E: 3,008 numeric cases, the Γ identity, and the harmonic and discrepancy bounds. Removing either
    constant factor fails, as expected.)
* **(8).** Iterate (4) for Y. Since Y ≤ U·n, integrability is trivial. ✓
* **(9).** M_{j+1} − M_j = (X_{j+1} − (1 + δ/j)X_j)/P_δ(m, j+1). Splitting off the conditional mean gives
  A-increments (E[X_{j+1}|F_j] − (1+δ/j)X_j)/P_δ ≥ 0 and martingale increments as stated. The conditional
  variance is at most the raw second moment. ✓
* **(10).** P_λ(m,j) ≤ e^λ (j/m)^λ and P_δ(m,j+1) ≥ P_δ(m,j) ≥ e^{−δ²}(j/m)^δ. This gives the summand bound
  e^{λ+2δ²}(j/m)^{−q}/j. For a decreasing summand, the sum is at most the first term plus the integral, which is
  at most 1/m + 1/q. ✓
  * (E: the numeric sums are 1.88–2.26, against a bound of 6.48–7.56.)
  * Negative control: at q = 0 the sum grows like log J.
* **Convergence.** The a.s. convergence of Z is the standard maximal-inequality / Borel–Cantelli argument; the
  proof's version is correct. Z_∞ has conditional mean 0 and conditional second moment ≤ D_- N_m (L² limit).
  Since A is nondecreasing and M ≥ 0, M converges in [0, ∞]. ✓
* **(11) and the δ = 0.542 remark.** Both are correct. Neither is needed for the theorem.
* **(12).**
  * E[J 1{J ≤ a/2}] ≤ a/2, which still holds when J < 0.
  * Cauchy–Schwarz on {J > a/2} gives P(J > a/2) ≥ a²/(4·E J²) ≥ a²/(4(a² + D_- N_m)).
  * On {J > a/2}, M_∞ ≥ X_m + Z_∞ = J > a/2, so liminf X_n/n^δ ≥ (a/2)·e^{−δ²}·m^{−δ} > 0.
  * The bound increases in a and decreases in N_m, which gives p_*.
  * ✓ (E: 20,000 random laws. The planted "≥ 1/4" version without the variance term fails.)
* The proof notes that the whole-tree restart bound is not uniform, and that this is why the bounded branch is
  used. That remark is correct and is the right motivation.

## 5. §4: reaching a positive state in a branch

* **Pattern existence.**
  * E Φ_n^- ≥ a_0 P_δ(3,n) ≍ n^δ.
  * Every depth-one core vertex was promoted at r, with probability L_r/(nD_r) ≤ 1/n per step. So the expected
    number of core vertices at depth ≤ 1 is at most 1 + Σ_{j=3}^{n−1} 1/j = O(log n).
  * If no reachable history had a depth-≥2 core vertex of positive weight, then E Φ_n^- = O(log n), a
    contradiction. ✓ (D: exact enumeration of all 20,160 histories to n = 9 confirms the depth-≤1 increment
    E[L_r/D_r]/n and the bound.)
* **Replay invariance.**
  * The copy of r has the same L and one more core neighbour.
  * A vertex at pattern depth ≥ 2 sees only its own (L, c) and its neighbours' (L, c, D) data (depth ≥ 1),
    plus its creator's birth data (depth ≥ 1). All of these agree between pattern and copy.
  * ✓ (C: 14,310 comparisons of raw certs9 states for four (K0, D_cap) pairs and of raw upper states, made after
    every replay step with exterior noise at p and its leaves: 0 mismatches.)
  * Depth-1 states do differ, 6,676 times, so "depth ≥ 2" is genuinely needed.
  * Negative control: a planted **radius-2** table (own state plus the parent's full state) gives 4,084
    mismatches at depth ≥ 2. The locality hypothesis is therefore essential (issue M1).
* **The waiting rule.** If the class probability is ≥ b/n, the class step happens a.s., because
  Π(1 − b/n) → 0. If the subclass makes up at least a fraction h of the class, the first class step falls in the
  subclass with probability ≥ h, by summing over the first-occurrence time. ✓
* **Stage 1.** The gain rate L_u + ρ_u is at least L_u ≥ 1, and at least the promotion rate L_u/D_u. So each
  next step at u is a gain with conditional probability ≥ 1/2. At most H − D_u gains are needed, giving
  ≥ 2^{−H}. Gains keep L_u ≥ 1. After that, until the first promotion at u, L_u/D_u is nondecreasing (gains at
  fixed c), so the qualifying promotion happens a.s. ✓ (A checks the rate inequality exactly; F checks the
  empirical frequency heuristically.)
* **Replay rates.**
  * Every inside vertex has an inside neighbour of degree ≤ k: D_v ≤ (k−1) + 1. So its receiving probability
    is ≥ 1/(nk).
  * P(W ∈ R) ≤ P(V ∈ R ∪ {p}) = (k+1)/n.
  * So each prescribed target is hit with probability ≥ 1/(k(k+1)), uniformly in the exterior history. This
    gives (13). ✓ (B: exact, 1,307 checks.)
  * At the end, N ≤ k_* and X ≥ a_* (f_- ≥ 0). ✓ (C: the branch has exactly k_* vertices.)

## 6. §5: from a uniform chance to probability one

* **(14).** Condition at three stopping sizes: stage 1 succeeds (≥ 2^{−H}), the qualifying promotion happens
  (a.s.), and the replay ends (≥ p_pattern; each inside step happens a.s.). Then apply (12) at the end of the
  replay, conditionally on that history atom. Because X_n ≤ Φ_n^-, this gives P(E | F_m) ≥ p_0 > 0 almost surely
  for every m. ✓
* **The 0–1 step.** P(E ∩ A) ≥ p_0 P(A) holds on the algebra ∪_m F_m. It extends by approximation to
  σ(∪F_m) ∋ E^c, which gives P(E^c) = 0. This is equivalent to Lévy's 0–1 law, and the approximation lemma's
  sketch (closure under complements and under countable unions via finite subunions) is correct. ✓
* **Conclusion.** C_n ≥ Φ_n^-/F. Also M_n^{(3)} = Φ_n^-/P_δ(3,n) converges a.s. in [0, ∞] and
  P_δ(3,n)/n^δ → c ∈ (0, ∞). So Φ_n^-/n^δ has an a.s. limit, and since P(E) = 1 that limit lies in (0, ∞]. ✓
  With δ = 0.543 and λ = 0.586 this gives **liminf C_n/n^{0.543} > 0 a.s.**, hence
  liminf log C_n/log n ≥ 0.543 a.s.

## 7. The final application

The proof needs these structural properties of the two concrete tables, and every one holds:

* **certs9:** the state is a function of the current tree, built from own (L, c) and neighbours' (L, c, D).
  That is radius 1.
* **Upper table:** own (L, c), neighbours' (L, c), and the creator's (L, c, D) at birth, which is static.
  Also radius 1.
* Both tables are bounded, flat in L beyond L_max = 60, and constant for c > c_max = 10.
* The verifier conditions hold: L_max = 60 > K0 = 2 and L_max + 1 = 61 ≥ D_cap = 6.
* A finite H exists (H = 60).

## 8. Issues found (none affects the δ = 0.543 conclusion)

### Status after the revision (re-read 2026-11-22)

| issue | status | what the revised `AS_LOWER_PROOF.md` now says |
|---|---|---|
| M1 | **resolved** | The theorem paragraph defines the class explicitly: x_u is a function of (i) own (L_u, c_u); (ii) the current (L_v, c_v, D_v) of the core neighbours, optionally with the parent distinguished; (iii) optionally a static birth label set by the creator's (L, c, D) at promotion. A cutoff H ≥ 3 makes f_± flat under gains at L_u ≥ H and makes the descriptor of any neighbour, parent or creator of degree ≥ H permanently saturated. The parent variant records only the parent's (L_p, c_p, D_p) with `P_top` above a cap, and "any shell" became "any coarsening of such a radius-1 state". |
| M2 | **resolved** | The exponents are now 0 < δ ≤ λ < 2δ. |
| M3 | **resolved** | `p_{\mathrm{pattern}}` in (13) and (14), `$\{p,v\}$`, `$\log(1+a/j)$`, `$\bigcup_m\mathcal F_m$`. A grep finds no remaining `\$`. |
| M4 | **resolved** | §1 now says "for the two final tables of `REPORT.md` §3.5, H = 60 suffices". |
| M5 | unchanged | This is a dependence on the certificates, not a flaw. |
| N1 | new nit | The radius-R side remark; see below. |

**I checked the new M1 hypothesis against every step that uses it.** It is exactly what each step needs:

* **§1 jump counts.** A gain or promotion at u changes only u's own (L, c), plus the newborn. So only u and
  the states that read u's data can change, and those are at most c_u ≤ D_u − L_u neighbours. When D_u ≥ H their
  descriptors of u are saturated, and the birth labels are static. That gives at most H + 1 changed summands,
  and no change at all when L_u ≥ H.
* **§2 exterior invariance.** The only exterior vertex visible inside is p. Its descriptor in v is saturated
  because D_p ≥ H, and v's birth label is static.
* **§4 replay.** A vertex at pattern depth ≥ 2 reads only data at depth ≥ 1, and that data is unchanged in the
  copy. The radius-2 negative control in §9 C (4,084 mismatches) shows that this locality is really needed.
* **The final tables satisfy the hypothesis.** certs9 depends on own (L, c) and neighbours' (L, c, D) and
  saturates to `top` at D ≥ 6. The upper table depends on own (L, c), neighbours' (L, c), and a static P built
  from the creator's (L, c, D) at birth. Both are flat for L ≥ 60, and H = 60 works.

**N1 (new nit, side remark only).** The revised theorem paragraph adds: "For radius-R local tables, the same
proof applies verbatim upon replacing pattern depth 2 in §4 by depth R+1, since the expected number of core
vertices at depth at most R is O((log n)^R)."

* The counting is right. The depth-(j+1) core vertices are promoted at depth-j vertices, at rate at most
  E[#depth-j core]/n, so induction gives O((log n)^R). The replay step is also right at depth R + 1.
* **"Verbatim" is not right for R ≥ 2.** §2's exterior invariance uses the fact that "the only exterior vertex
  seen by an inside vertex is p, and only v sees it". For R ≥ 2, v (and inside vertices up to distance R − 1
  from v) also read data of p itself and of p's other neighbours. D_p ≥ H saturates only p's descriptor, not
  what lies behind p.
  * Example (a hand argument, not run): v's state records p's full certs9 state. That state starts with L_p, and
    L_p changes whenever p gains a leaf, even though D_p ≥ H.
  * So the saturation hypothesis must be strengthened: everything seen through a vertex of degree ≥ H must be
    hidden, as the `P_top` convention does. The changed-summand bound H + 1 also becomes O(H^R).
* My own earlier M1 suggestion ("for radius-2 tables use pattern depth ≥ 3") had the same omission; it covered
  only the replay step.
* **Fix:** drop "verbatim", or add the strengthened saturation hypothesis to the remark.
* This affects neither the radius-1 theorem as stated nor the δ = 0.543 corollary.

### Issues as reported in the first audit (kept for the record)

* **M1 (statement scope; should be fixed in the text).** The theorem covers "local tables of the kind in
  REPORT §§3.3–3.5, including the rooted version which records the parent's state…", and says "any shell that
  drops information is allowed".
  * REPORT §§3.3–3.5 never defines a rooted version.
  * §4's replay claim ("a vertex of pattern depth at least two sees neither that extra edge nor the copy of r")
    holds only for tables that are **radius-1 local**: own and neighbours' (L, c, D), plus static creator birth
    data.
  * If "parent's state" meant the parent's full local state, the claim fails at depth 2. My negative control
    shows this (4,084 mismatches).
  * Elsewhere ("records the parent's L,c separately", "P_top above its own cap") the proof seems to intend the
    radius-1 reading, under which everything is correct.
  * **Fix:** state locality and the degree-H saturation as explicit hypotheses, or restrict the theorem to the
    named tables. Alternatively, for radius-2 tables use pattern depth ≥ 3: E[#core at depth ≤ 2] = O(log² n),
    which is still o(n^δ).
  * "Any shell" is fine if it means "any function of such a state".
* **M2 (nit).** The hypothesis 1/2 < δ is never used. The proof needs only δ > 0 and λ < 2δ, as the proof itself
  later says.
* **M3 (typos).** `p_{mathrm{pattern}}` in (13) is missing a backslash. There are stray `\$` in "\$\{p,v\}$",
  "\$\log(1+a/j)$" and "\$\bigcup_m".
* **M4 (nit).** H is never given a value for the final certificates; H = 60 works. The constants in (3), (6) and
  (10) are loose, which is harmless.
* **M5 (dependence, not a flaw).** The result is conditional on REPORT's computer-assisted certificates, which
  this audit did not re-run. The bridge from those certificates to hypothesis (1) is checked in §9 A.

## 9. Independent check script and its exact output

Script: `checks/strogatz-carbon-as-lower-claude/check_as_lower.py`
(sha256 `a822ad16720ccd89e05bf15d311194e86a15e545ae2e4b0258900779f4bf3a28`).
It covers:

* **A:** one-step facts on 170 random trees, some grown with forced high-degree receivers, all reachable;
* **B:** branch restriction, 60 branches, all outcomes enumerated at every intermediate tree;
* **C:** replay invariance, 120 random patterns with k_* ≤ 28;
* **D:** exact enumeration of all histories up to n = 9;
* **E:** the analytic inequalities;
* **F:** heuristic Monte Carlo.

There are 14 planted-error negative controls, and all must be detected.

Command, run from the repository root (seed 20261122, wall time 4m13s):

```
cd checks/strogatz-carbon-as-lower-claude
taskset -c 0 env OMP_NUM_THREADS=1 ../../.venv/bin/python3 -I check_as_lower.py > check_as_lower.log
```

Output (`check_as_lower.log`, verbatim):

```
seed = 20261122

=== PART A: one-step facts on random trees (exact) ===
[PASS] receiver law P(W=w)=n^-1 sum_(z~w) 1/D_z (backward = forward enumeration), core rate L+rho, leaf rate L/D: 170 trees, bad=0
[NEG-CONTROL DETECTED] planted preferential-attachment law D_w/(2(n-1)) differs from true law on 170/170 trees
[PASS] identity (2) sum_u rho_u = sum_v c_v/D_v <= C_n and S_n <= C_n: bad=0
[PASS] stage-1: gain rate L+rho >= promotion rate L/D for every core vertex: bad=0
[PASS] per-vertex drift (1) by brute force == REPORT event-list formula at true rates (certs9 (K0,Dcap) in [(1, 4), (2, 6), (3, 5), (2, 4), (2, 5)]; upper P_max in [3, 6, 20]): 5408 vertex comparisons, mismatches=0
[PASS] certs9 true rates inside REPORT intervals: credit 3380 checks bad=0; E3/E4/E5/E7 3038 checks bad=0
[PASS] upper (L,c,m1,b,P) true rates inside REPORT intervals: credit 2028 checks bad=0; E3/E5 717 checks bad=0
[NEG-CONTROL DETECTED] planted event-list error 'E4rate' changes the formula: mismatching vertices=630
[NEG-CONTROL DETECTED] planted event-list error 'nb' changes the formula: mismatching vertices=1324
[NEG-CONTROL DETECTED] planted event-list error 'noE7' changes the formula: mismatching vertices=566
[NEG-CONTROL DETECTED] planted event-list error 'U-E4rate' changes the formula: mismatching vertices=378
[PASS] J1 gain: changed entries within receiver+neighbours, count <= max(3,Dcap): tests=3380 bad=0
[PASS] J2 promotion: <=3 changed summands, old neighbour changes only if creator was pendant: tests=22275 bad=0
[PASS] J3 gain at core w with L_w >= L_max (L_max>K0, L_max+1>=Dcap): no truncated entry changes: tests=1667 bad=0
[PASS] J4a gain at core w with D_w >= H: no neighbour entry changes: tests=1796 bad=0
[PASS] J4b gain at core w with L_w >= H: no truncated entry changes: tests=1429 bad=0
[PASS] J4c promotion at creator with D >= H: no old neighbour entry changes: tests=20272 bad=0
[PASS] J4d every step: <= H+1 changed truncated summands: tests=25655 bad=0
[PASS] J4U upper table: receiver/creator of degree >= H changes no neighbour entry: tests=10840 bad=0
[NEG-CONTROL DETECTED] planted L_max+1 < D_cap (K0=1,Dcap=6,L_max=2): gains at L>=L_max that change a truncated entry = 152
[PASS] second-moment bound (3) E[(dPhi)^2] <= B^2(L_max+2)C/n, B=F max(3,Dcap), plus P(W core, L_W<L_max) <= (L_max+1)C/n: tests=850 bad=0; max observed ratio LHS/RHS = 0.0398

=== PART B: restriction across a saturated parent edge (section 2), exact ===
[PASS] saturated parent (D_p >= H=6 >= D_cap=5, P_max=5): every exterior receiver (incl. p and leaves of p) leaves all raw inside certs9 and upper states unchanged: exterior outcomes=33073 bad=0
[NEG-CONTROL DETECTED] planted unsaturated parent (D_p < D_cap): exterior outcomes that change an inside state = 231
[PASS] restricted bookkeeping (4): I_(n+1) = I_n + newborn iff creator in I_n; X_(n+1)-X_n = sum over I_n of summand changes + newborn, for X and Y, every outcome; E[dX] = n^-1 sum_(u in I) (1)_u: checks=46050 bad=0
[NEG-CONTROL DETECTED] planted bookkeeping error (newborn omitted) changes dX on 6420 outcomes
[PASS] (5): sum_(u in I) rho_u = sum_(z in I) |N(z) cap I|/D_z + 1/D_p <= N+1 <= 2N: checks=1307 bad=0
[PASS] (6): P(X changes) <= (H+2)N/n, E[dX^2] <= (H+3)B_F^2 N/n, E[dY^2] <= (H+3)B_U^2 N/n: checks=1076 bad=0; max ratio for X = 0.01013
[PASS] replay rates: each inside vertex has P(W=x) >= 1/(nk), P(W in branch) <= (k+1)/n: checks=1307 bad=0

=== PART C: replay of a pattern inside a saturated branch (section 4) ===
[PASS] depth>=2 pattern vertices: raw certs9 states ((K0,Dcap) in [(1, 4), (2, 6), (2, 5), (3, 7)]) and upper (L,c,m1,b,P<= 8) states identical in pattern and replayed branch, after every replay step with exterior noise at p: comparisons=14310 mismatches=0
[PASS] copy of r: same L, one extra core neighbour: bad=0
[PASS] branch has exactly k_* vertices at the end: bad=0
[INFO] depth-1 vertices whose certs9 state differs (expected, the proof only uses depth>=2): 6676; upper-state differences at depth 1: 2599
[NEG-CONTROL DETECTED] depth-1 states do change in the copy, so the restriction to depth>=2 is needed
[NEG-CONTROL DETECTED] planted radius-2 table (own state + parent's full state): depth>=2 mismatches = 4084

=== PART D: exact enumeration of all receiver histories from T_3 up to n = 9 ===
[INFO] histories enumerated to n=9: 20160
[INFO] n=3: E[C_n]=1.000000  E[#core depth<=1]=1.000000  bound 1+H=1.000000
[INFO] n=4: E[C_n]=1.333333  E[#core depth<=1]=1.333333  bound 1+H=1.333333
[INFO] n=5: E[C_n]=1.583333  E[#core depth<=1]=1.541667  bound 1+H=1.583333
[INFO] n=6: E[C_n]=1.797222  E[#core depth<=1]=1.697222  bound 1+H=1.783333
[INFO] n=7: E[C_n]=1.989468  E[#core depth<=1]=1.822897  bound 1+H=1.950000
[INFO] n=8: E[C_n]=2.166659  E[#core depth<=1]=1.929072  bound 1+H=2.092857
[INFO] n=9: E[C_n]=2.332471  E[#core depth<=1]=2.021397  bound 1+H=2.217857
[PASS] sum rule E[C_(n+1)] - E[C_n] = E[S_n]/n exactly, n=3..8
[PASS] E[#core of depth<=1] increments by E[L_r/D_r]/n <= 1/n and is <= 1+sum_(j=3)^(n-1) 1/j
[PASS] some reachable history has a core vertex of depth >= 2 (E[C] > E[depth<=1 core])

=== PART E: analytic inequalities (7), (10), (12), (13) ===
[PASS] (7) e^(-a^2)(n/m)^a <= P_a(m,n) <= e^a (n/m)^a, a in [0.1,3], m>=3: 3008 cases, bad=0
[NEG-CONTROL DETECTED] planted (7) without the e^a factor (P <= (n/m)^a) fails in 1474 cases
[NEG-CONTROL DETECTED] planted (7) without the e^(-a^2) factor (P >= (n/m)^a) fails in 1120 cases
[PASS] harmonic 0 <= sum 1/j - log(n/m) <= 1/m and |log(1+a/j)-a log(1+1/j)| <= (a^2+a)/(2j^2): bad=0
[PASS] P_a(3,n) = 2 Gamma(n+a)/(Gamma(n) Gamma(3+a)) (so P_a(3,n)/n^a has a positive finite limit)
[INFO] (10) m=3: partial sum to J=300000 = 1.86225, + tail bound 0.02049 = 1.88274 <= e^(lam+2del^2)(1/m+1/q) = 7.56094
[INFO] (10) m=10: partial sum to J=300000 = 1.94914, + tail bound 0.03742 = 1.98655 <= e^(lam+2del^2)(1/m+1/q) = 6.80485
[INFO] (10) m=100: partial sum to J=300000 = 1.95960, + tail bound 0.11832 = 2.07792 <= e^(lam+2del^2)(1/m+1/q) = 6.51321
[INFO] (10) m=1000: partial sum to J=300000 = 1.88415, + tail bound 0.37417 = 2.25832 <= e^(lam+2del^2)(1/m+1/q) = 6.48405
[PASS] (10) per-term bound and total bound with delta=0.543, lambda=0.586, q=0.5
[NEG-CONTROL DETECTED] planted lambda = 2 delta (q=0): the sum in (10) grows ~log J: J=1000: 5.186, J=10000: 7.276, J=100000: 9.367, J=1000000: 11.458, J=1999999: 12.087
[PASS] (12) for any J with E J = a > 0: P(J > a/2) >= a^2/(4 E J^2) on 20000 random laws: bad=0
[NEG-CONTROL DETECTED] planted variance-free version P(J > a/2) >= 1/4 fails on 275 laws
[INFO] (13) k_*=4: p_pattern = 1/72 ~ 1.389e-02 > 0
[INFO] (13) k_*=6: p_pattern = 1/43200 ~ 2.315e-05 > 0
[INFO] (13) k_*=10: p_pattern = 1/658409472000 ~ 1.519e-12 > 0

=== PART F: Monte Carlo (heuristic only, proves nothing) ===
[INFO] n=1000: C_n/n^0.543 over 12 runs: min=0.211 median=0.940 max=1.386
[INFO] n=10000: C_n/n^0.543 over 12 runs: min=0.209 median=0.861 max=1.292
[INFO] n=100000: C_n/n^0.543 over 12 runs: min=0.220 median=0.894 max=1.384
[INFO] n=200000: C_n/n^0.543 over 12 runs: min=0.224 median=0.898 max=1.416
[INFO] stage 1 (H=8): reached D_u >= H before any promotion at u in 142/150 trials; mean of the proof's lower bound 2^-(H-D_u) = 0.8688
[PASS] stage-1 empirical success frequency is above the proof's bound

=== SUMMARY ===
failures: 0
negative controls missed: 0
ALL CHECKS PASSED
```

Summary of the first run: **32 PASS lines, 0 FAIL, 14 negative controls, all detected.**

A second run with a different seed (`check_as_lower.py 7 > check_as_lower_seed7.log`) gave the same result.
Key lines, verbatim:

```
[PASS] per-vertex drift (1) by brute force == REPORT event-list formula at true rates (certs9 (K0,Dcap) in [(1, 4), (2, 6), (3, 5), (2, 4), (2, 5)]; upper P_max in [3, 6, 20]): 5280 vertex comparisons, mismatches=0
[PASS] saturated parent (D_p >= H=6 >= D_cap=5, P_max=5): every exterior receiver (incl. p and leaves of p) leaves all raw inside certs9 and upper states unchanged: exterior outcomes=35006 bad=0
[PASS] depth>=2 pattern vertices: raw certs9 states ((K0,Dcap) in [(1, 4), (2, 6), (2, 5), (3, 7)]) and upper (L,c,m1,b,P<= 8) states identical in pattern and replayed branch, after every replay step with exterior noise at p: comparisons=16835 mismatches=0
[NEG-CONTROL DETECTED] planted radius-2 table (own state + parent's full state): depth>=2 mismatches = 5417
failures: 0
negative controls missed: 0
ALL CHECKS PASSED
```

Summary of the second run: 32 PASS lines, 0 FAIL, 14 negative controls, all detected.

What the script does **not** do:

* It does not check the 1.2-million-state certificate tables.
* The tree checks are exhaustive per tree (every receiver), but they cover random finite trees; they are not a
  proof for all trees. Each claim they test is also proved by hand above.
* Part F is a simulation only.
