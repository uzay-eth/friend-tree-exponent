VERDICT: CORRECT

The proof of Theorem B (the general transfer principle from (H1)–(H5) to the almost-sure lower bound `liminf log C_n / log n >= delta`) is mathematically sound, every inequality and constant checks out, and hypotheses (H1)–(H5) genuinely hold for the `certs9.py` tracked-degree certificate (`Lmax=60, cmax=10, K0=2, Dcap=6`, `delta = 543/1000`) with `b0 = 6` and `R = 126`. Only three minor expository omissions/nitpicks were found (all trivially repaired in one line without affecting any statement or constant), and our independent brute-force exact rational verification over 120 random friend trees (`n <= 150`, 20,516 `(V, W)` moves) plus extremal corner-case trees confirms (H1)–(H5) and catches all planted wrong constants. Assuming the previously independently re-verified exact box-and-tail drift check (H2) of `certs9.py` (also spot-checked by brute force on all our trees), Theorem A follows.

---

## 1. Line-by-line check of Theorem B

### Setting and hypotheses (H1)–(H5)
- **One-step dynamics (i)–(ii):** Valid for `n >= 3`. In any tree on `n >= 3` vertices, no two leaves are adjacent, so every leaf has a unique core neighbour. From `P(W_n = w | T_n) = (1/n) sum_{v ~ w} 1/D_v`: a core vertex `u` is picked with probability `(L_u + rho_u)/n` (effect `L_u -> L_u + 1`); a leaf `l` adjacent to core vertex `u` is picked with probability `1/(n D_u)` (creating core vertex `l` with `(1, 1)` and changing `u` to `(L_u - 1, c_u + 1)`), giving total conversion probability at `u` equal to `L_u/(n D_u) <= 1/n`. Every core neighbour `v` has `D_v >= 2`, so `rho_u <= c_u/2` and `L_u + rho_u <= D_u`.
- **Well-definedness of constants:** (H4) guarantees that a conversion (which occurs with positive probability at step `n = 3`) produces a core vertex `w` with `f(x_w) >= eta > 0`, so the finite range of `f` in (H1) contains at least one positive value; hence `fmin > 0` and `F >= fmin > 0` exist, and `B := (b0 + 1) F > 0` (as `b0 >= 0`).
- **Step bound `|Phi_{n+1} - Phi_n| <= B` in (H3):** Valid. Since `core(T_n) subset core(T_{n+1})` and `|core(T_{n+1}) \ core(T_n)| <= 1`, `Phi_{n+1} - Phi_n` is a sum of at most `b0` differences `f(x_v^{(n+1)}) - f(x_v^{(n)})` (each in `[-F, F]` since `0 <= f <= F`) plus at most one newborn term in `[0, F]`. Thus `|Phi_{n+1} - Phi_n| <= (b0 + 1) F = B`.

### Lemma 1 (conversions recur)
- **Validity:** Valid up to a trivial indexing nitpick for `M = 1`.
- **Line-by-line check:**
  - `C_n` is nondecreasing since degrees never decrease (`C_3 = 1`, `C_{n+1} - C_n in {0, 1}`).
  - *Minor nitpick:* Line 61 writes "Fix an integer `M >= 1` and let `n >= M^2 + M`." When `M = 1`, `M^2 + M = 2`, whereas the next sentence relies on `n >= 3` (which the text itself notes parenthetically: "(for `n >= 3` every leaf is adjacent to a core vertex)"). **Repair:** write `n >= max(3, M^2 + M)` (or restrict to `M >= 2`, since `{sup C_n <= 1} subset {sup C_n <= 2}`). Theorem survives completely unchanged.
  - For `n >= max(3, M^2 + M)` on `{C_n <= M}`, we have `1 <= C_n <= M` and `n - C_n >= M^2` leaves, each adjacent to one of the `C_n <= M` core vertices. By pigeonhole, some core vertex `u` has `L_u >= (n - C_n)/C_n >= M`, while `c_u <= C_n - 1 <= M - 1`.
  - Since `L/(L + c) = 1/(1 + c/L)` is nondecreasing in `L` and nonincreasing in `c`, `L_u / (n D_u) >= M / (n (2M - 1)) >= 1/(2n)` for all `M >= 1`.
  - On `{C_n <= M for all n}`, the sum over `n >= 3` of `P(conversion at step n | F_n)` is at least `sum_{n >= max(3, M^2 + M)} 1/(2n) = infinity`. Because `{conversion at step n} in F_{n+1}`, the conditional Borel–Cantelli lemma (Durrett, 5th ed., Thm 4.3.4: `{B_n i.o.} = {sum P(B_n | F_{n-1}) = infinity}` a.s.) implies that a.s. on `{sup_n C_n <= M}`, conversions occur infinitely often, contradicting `sup_n C_n <= M`. Hence `P(sup_n C_n <= M) = 0` for all `M`, so `C_n -> infinity` a.s.

### Lemma 2 (`Phi` is unbounded)
- **Validity:** Valid as written.
- **Line-by-line check:**
  - **Submartingale & Doob decomposition:** `0 <= Phi_n <= n F` is bounded (hence integrable) and `F_n`-measurable. By (H1)–(H2), `E_n[Delta_n] >= (delta/n) Phi_n >= 0`, so `(Phi_n)_{n >= 3}` is a submartingale. Its Doob decomposition `Phi_n = Phi_3 + M_n + A_n` with `A_n = sum_{k=3}^{n-1} E_k[Delta_k]` and `M_n = sum_{k=3}^{n-1} (Delta_k - E_k[Delta_k])` satisfies `A_n` nondecreasing, `A_n >= delta * sum_{k=3}^{n-1} Phi_k / k`, and `(M_n)` a martingale with increments `|M_{k+1} - M_k| = |Delta_k - E_k[Delta_k]| <= |Delta_k| + E_k[|Delta_k|] <= 2B` by (H3).
  - **Use of bounded-increment martingale dichotomy:** On `E := {sup_n Phi_n < infinity}`, since `Phi_3 >= 0` and `A_n >= 0`, we have `M_n = Phi_n - Phi_3 - A_n <= Phi_n <= sup_m Phi_m < infinity`, so `limsup_n M_n < infinity`. Durrett (5th ed.), Thm 4.3.1 states that a martingale with uniformly bounded increments a.s. either converges to a finite limit or has `limsup = +infinity` and `liminf = -infinity`. Therefore, a.s. on `E`, `M_n` converges to a finite limit `M_inf` (in particular, `(M_n)` is bounded on `E` a.s.).
  - **Convergence of `A_n` and `Phi_n`, and the step `Phi_inf = 0`:** A.s. on `E`, `A_n = Phi_n - Phi_3 - M_n` is bounded above; being nondecreasing, `A_n` converges to a finite limit `A_inf`. Hence `Phi_n = Phi_3 + M_n + A_n` also converges to a finite limit `Phi_inf := Phi_3 + M_inf + A_inf >= 0`. Crucially, because `Phi_n` is already known to converge (not merely bounded), `sum_{k=3}^infinity Phi_k / k <= A_inf / delta < infinity` together with `delta > 0` and `sum 1/k = infinity` forces `Phi_inf = 0` (otherwise `Phi_k >= Phi_inf/2` eventually, causing the series to diverge).
  - **Contradiction:** By Lemma 1, conversions occur at infinitely many steps `n_1 < n_2 < ...` a.s.; by (H4) and `f >= 0`, `Phi_{n_i + 1} >= eta > 0` for every `i`, contradicting `Phi_n -> 0`. Hence `P(E) = 0`.

### Lemma 3 (reciprocal supermartingale once `Phi` is large)
- **Validity:** Valid as written.
- **Line-by-line check:**
  - **Algebraic identity & lower bound on `Phi + Delta`:** On `{Phi >= 2B}`, `|Delta| <= B` implies `Phi + Delta >= Phi - B >= Phi/2 > 0` (so `1/Phi_{n+1}` is well-defined and bounded by `2/Phi <= 1/B`). Direct expansion confirms
    `1/Phi - Delta/Phi^2 + Delta^2 / (Phi^2 (Phi + Delta)) = (Phi(Phi + Delta) - Delta(Phi + Delta) + Delta^2) / (Phi^2 (Phi + Delta)) = 1/(Phi + Delta)`.
    Since `Delta^2 >= 0` and `Phi + Delta >= Phi/2 > 0`, the third term is at most `2 Delta^2 / Phi^3`.
  - **Bound `Delta^2 <= B(Delta + 2 Delta^-)`:** Valid. Since `|Delta| <= B`, `Delta^2 <= B |Delta|`, and checking signs (`Delta >= 0` gives `Delta + 0 = |Delta|`; `Delta < 0` gives `Delta + 2(-Delta) = -Delta = |Delta|`) confirms `|Delta| = Delta + 2 Delta^-`. Substituting into the previous line and taking `E_n[.]` gives equation (1):
    `E_n[1/Phi_{n+1}] <= 1/Phi - (E_n[Delta]/Phi^2)(1 - 2B/Phi) + 4B E_n[Delta^-]/Phi^3`.
  - **Applying (H2):** Since `Phi >= 2B`, `1 - 2B/Phi >= 0`, so `-(1 - 2B/Phi)/Phi^2 <= 0`. Combining with `E_n[Delta] >= (delta/n) Phi` gives `-(E_n[Delta]/Phi^2)(1 - 2B/Phi) <= -(delta/(n Phi))(1 - 2B/Phi)`, with no upper bound on `E_n[Delta]` required.
  - **Union bound for `P_n(Delta < 0)`:** If `Delta < 0`, since any newborn contributes `f >= 0`, at least one pre-existing core vertex `v in core(T_n)` must have `f(x_v^{(n+1)}) < f(x_v^{(n)})`, which forces `f(x_v^{(n)}) > 0` (as `f >= 0`) and `f(x_v)` to change. Every vertex in `S_n^+ := {v in core(T_n) : f(x_v^{(n)}) > 0}` contributes at least `fmin`, so `|S_n^+| <= Phi / fmin`. Because `(T_n)` is Markov and `x_v` is a function of `T_n` alone, `P_n(. ) = P(. | T_n)`, so by (H5) each `v in S_n^+` changes `f(x_v)` with `P_n`-probability at most `R/n`. Union bound gives `P_n(Delta < 0) <= |S_n^+| (R/n) <= R Phi / (n fmin)`. Since `Delta^- <= B * 1_{Delta < 0}`, `E_n[Delta^-] <= B R Phi / (n fmin)`.
  - **Constants `kappa` and `K0`:** Substituting both bounds into (1) gives
    `1/Phi - delta/(n Phi) + 2 B delta / (n Phi^2) + 4 B^2 R / (n fmin Phi^2) = (1/Phi)[1 - (1/n)(delta - kappa/Phi)]`
    with `kappa = 2 B delta + 4 B^2 R / fmin > 0`. On `{Phi >= K0}` with `K0 = max(2B, kappa/eps)`, both `Phi >= 2B` and `kappa/Phi <= eps` hold, yielding `E_n[1/Phi_{n+1}] <= (1/Phi)(1 - (delta - eps)/n)`. Exact.

### Final stopping-time argument (Proof of Theorem B)
- **Validity:** Mathematically valid; two routine details (`Phi_tau > 0` and integrability/shift-measurability of `G_{sigma+k}`) are left implicit and are easily supplied.
- **Line-by-line check:**
  - **Stopping times `sigma, tau`:** `sigma = inf{n >= 3 : Phi_n >= K1}` and `tau = inf{n >= sigma : Phi_n < K0}` are `(F_n)`-stopping times; `sigma < infinity` a.s. by Lemma 2, and since `K1 > K0`, `Phi_sigma >= K1 > K0` implies `tau >= sigma + 1` on `{sigma < infinity}`.
  - **Well-definedness of `G_n` (implicit point 1):** For `j >= 3` and `0 < delta - eps < 1`, each factor `1 - (delta - eps)/j` lies in `(2/3, 1)`, so `p_m in (0, 1]` (with empty product `p_sigma = 1`). Moreover, `Phi_{min(n, tau)}` is strictly positive: for `sigma <= n < tau` we have `Phi_n >= K0 >= 2B > 0`, and on `{tau < infinity}` we have `tau - 1 >= sigma` with `Phi_{tau - 1} >= K0 >= 2B`, so by (H3) `Phi_tau >= Phi_{tau - 1} - B >= K0 - B >= B > 0`. (The written proof omits stating `Phi_tau > 0` explicitly before writing `G_tau = 1/(Phi_tau p_tau) > 1/K0`, though Lemma 3 already proves `Phi_{n+1} >= Phi_n/2 > 0` on `{Phi_n >= 2B}`.)
  - **Supermartingale property of `(G_{sigma+k})_{k >= 0}` w.r.t. `(F_{sigma+k})_{k >= 0}` (implicit point 2):**
    - *Integrability:* For fixed `k >= 0`, `p_{min(sigma+k, tau)} >= (2/3)^k` and `Phi_{min(sigma+k, tau)} >= B`, so `0 < G_{sigma+k} <= 1/(B (2/3)^k)` is uniformly bounded on `{sigma < infinity}`, hence integrable (setting `G_{sigma+k} = 0` on the null event `{sigma = infinity}`).
    - *Supermartingale inequality at random start `sigma`:* Let `A in F_{sigma+k}`. For each `s >= 3`, `A_s := A intersect {sigma = s} in F_{s+k}`. Writing `n = s + k`, on `{sigma = s}` the factors `p_n, p_{n+1}` are deterministic constants (depending only on `s, n`) with `p_{n+1} = p_n (1 - (delta - eps)/n)`, and `{tau <= n}, {tau > n} in F_n`. On `A_s intersect {tau <= n}`, `G_{n+1} = G_n`. On `A_s intersect {tau > n}`, `Phi_n >= K0` and `min(n+1, tau) = n + 1`, so Lemma 3 gives `E[G_{n+1} 1_{A_s, tau > n}] = E[E_n[1/Phi_{n+1}] / p_{n+1} * 1_{A_s, tau > n}] <= E[G_n 1_{A_s, tau > n}]`. Summing over `s >= 3` by monotone convergence (since `G >= 0`) gives `E[G_{sigma+k+1} 1_A] <= E[G_{sigma+k} 1_A]`. Finally `G_sigma = 1/Phi_sigma <= 1/K1`.
  - **Consequence (a) — maximal inequality:** For any nonnegative supermartingale `(Y_k)_{k >= 0}` and `lambda > 0`, optional stopping at `T ^ N` (`T = inf{k : Y_k >= lambda}`) gives `P(exists k : Y_k >= lambda) <= E[Y_0]/lambda`, and taking `c -> lambda^-` gives `P(sup_k Y_k >= lambda) <= E[Y_0]/lambda`. With `lambda = 1/K0` and `E[G_sigma] <= 1/K1`, this yields `P(sup_k G_{sigma+k} >= 1/K0) <= K0/K1`. On `{tau < infinity}`, `0 < Phi_tau < K0` and `0 < p_tau <= 1` give `G_tau = 1/(Phi_tau p_tau) > 1/K0` (attained at finite `k = tau - sigma`), so `P(tau < infinity) <= K0/K1`. Exact.
  - **Consequence (b) — convergence and bound on `p_n`:** By Doob's supermartingale convergence theorem, `G_{sigma+k}` converges a.s. to an a.s. finite limit, so `G* := sup_{n >= sigma} G_n < infinity` a.s. (and `G* >= G_sigma > 0`). On `{tau = infinity}` (intersected with the probability-1 event `{sigma < infinity, G* < infinity}`), `Phi_n = 1/(G_n p_n) >= 1/(G* p_n)` for all `n >= sigma`. Using `1 - x <= e^{-x}` for each factor of `p_n`,
    `p_n = prod_{j=sigma}^{n-1} (1 - (delta - eps)/j) <= exp(-(delta - eps) sum_{j=sigma}^{n-1} 1/j)`.
    Because `t |-> 1/t` is decreasing on `(0, infinity)`, `sum_{j=sigma}^{n-1} 1/j >= integral_sigma^n (1/t) dt = log(n / sigma)`, and multiplying by `-(delta - eps) < 0` yields `p_n <= (sigma / n)^{delta - eps}`. Both inequalities go in the right direction. Therefore `Phi_n >= (n/sigma)^{delta - eps} / G*` for all `n >= sigma`, so `log Phi_n / log n >= (delta - eps) - ((delta - eps) log sigma + log G*)/log n`, giving `liminf log Phi_n / log n >= delta - eps` on `{tau = infinity}`.
  - **Conclusion:** Since `{liminf log Phi_n / log n >= delta - eps}` does not depend on `K1`, letting `K1 -> infinity` gives `P(liminf log Phi_n / log n >= delta - eps) = 1`. By (H1), `Phi_n <= F C_n`, so `C_n >= Phi_n / F` gives `liminf log C_n / log n >= delta - eps` a.s. Intersecting over a countable sequence `eps_m -> 0` yields `liminf log C_n / log n >= delta` a.s.

---

## 2. Verification of (H3), (H4), (H5) for `certs9.py` (`Lmax=60, cmax=10, K0=2, Dcap=6`)

### Analytic audit of Section 3
- **(H3) with `b0 = 6`:**
  - *Type (i) (`W_n = w` core):* Only `(L_w, c_w) -> (L_w + 1, c_w)` changes. A vertex of degree `D_w >= 6` has class `top` before and after, so another core vertex `v` changes state `x_v` only if `v ~ w` and `D_w <= 5`, in which case `c_w <= D_w <= 5`. Hence at most `1 + 5 = 6` pre-existing core vertices change state and `0` are created. This bound `b0 = 6` is sharp: a vertex `w` with `(L_w, c_w) = (0, 5)` (class `O_5`) changes to `(1, 5)` (`top`) when `W_n = w`, changing `x_w` and all 5 of its core neighbours (`6` states and `6` values of `f`).
  - *Type (ii) (conversion at `u`):* `(L_u, c_u) -> (L_u - 1, c_u + 1)` (with `D_u` unchanged) and leaf `l` becomes newborn core vertex `w` with `(1, 1)`. If `c_u >= 2`, then since `L_u >= 1`, `u` was `O_{D_u}` or `top`, and afterwards `c_u' >= 3` with same `D_u`, so its class is unchanged (only `x_u` changes). If `c_u = 1`, `u` has only one pre-existing core neighbour `v` (at most 2 pre-existing states change). *Minor expository omission:* lines 182–184 do not explicitly mention the case `c_u = 0` (e.g. a star `S_n`), where `u` has no pre-existing core neighbours at all so only `x_u` changes (1 existing state, trivially `<= 2`). In all subcases at most 2 pre-existing states change and 1 core vertex is created.
  - Thus `b0 = 6` and `<= 1` created core vertex hold unconditionally.
- **(H4):** Every newborn `w` has `(L_w, c_w) = (1, 1)` and a single core neighbour `u`, whose post-conversion class is one of the 7 classes `{P_1, P_2, B, O_3, O_4, O_5, top}`. Checking all 7 exact rational values of `f` on the rebuilt table confirms all 7 are strictly positive with `eta = min = 0.359043 > 0`.
- **(H5) with `R = 126`:** If `f(x_v) > 0`, then `c_v <= cmax = 10` (as `f = 0` for `c > 10`). From the case analysis of (H3), `f(x_v)` can change only via: (a) `W_n = v` with `L_v < 60` (since for `L_v >= 60`, `f(L_v + 1, .) = f(60, .) = f(L_v, .)`), prob `(L_v + rho_v)/n <= (59 + 5)/n < (60 + 5)/n`; (b) conversion at `v`, prob `L_v/(n D_v) <= 1/n`; (c) `W_n = w` for a core neighbour `w` of `v` with `D_w <= 5`, prob `(L_w + rho_w)/n <= D_w/n <= 5/n` for each of `<= 10` neighbours; (d) conversion at a core neighbour `w` of `v`, prob `<= 1/n` for each of `<= 10` neighbours. Summing gives `<= (60 + 5 + 1 + 10*5 + 10*1)/n = 126/n`.

### Independent brute-force verification program
We wrote an independent verification script `checks/fable-as-lower-gemini/check_hypotheses.py` that:
1. Rebuilds the exact rational table from `certs9.Model(60, 10, 2, 6)` and `certs9.lower_table(M, 0.5433)`.
2. Simulates 120 random friend trees (`n` spanning `3..150`: one for every `n in 3..30`, 40 with `n in 31..100`, 32 with `n in 101..149`, and 20 at `n = 150`; 20,516 directed `(V, W)` moves total), plus 8 constructed extremal/corner-case trees (`T_3`, `S_4`, `S_70`, `S_200`, tight `b0 = 6`, `L >= 60`, high `R = 102`, and `c > 10` with a dead `(0, 2)` all-`top` vertex).
3. Computes `x_v = (L, c, m_1, m_2, b, n_3, n_4, n_5)` directly from each tree's adjacency list without using `certs9.py`'s transition rules, enumerates all `2(n - 1)` moves `(V, W)` with exact rational probability `1/(n deg(V))` (cross-checked grouped vs literal ungrouped on the first 10 trees with exact equality), builds `T_{n+1}` from scratch, and checks (H3), (H4), (H5), and (H2) in exact `Fraction` arithmetic.
4. Runs 6 negative controls with deliberately wrong parameters (`b0 = 2`, `b0 = 5`, `R = 1`, `R = 100`, `delta = 0.600`, and omitting the `L >= 60` cap on `S_200`), confirming every one is caught.

### Command run and real output (single core)
Command:
```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 taskset -c 0 \
  attempts/strogatz-fable/.venv/bin/python \
  checks/fable-as-lower-gemini/check_hypotheses.py \
  2>&1 | tee checks/fable-as-lower-gemini/run.log
```

Output (exit code 0, wall time 38.8s on 1 core):
```text
=== 1. Building certs9 lower table (Lmax=60, cmax=10, K0=2, Dcap=6, build=0.5433) ===
Table built in 34.7s: 1186328 entries, Perron rho = 1.000033
(H1) sup f = 0.639069, fmin = 3.611209e-03, zero entries = 18
     zero entries breakdown: 9 with D<2, 9 with D>=2
       dead D>=2 zero state: (0, 2, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 3, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 4, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 5, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 6, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 7, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 8, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 9, (0, 0), 0, (0, 0, 0))
       dead D>=2 zero state: (0, 10, (0, 0), 0, (0, 0, 0))
(H4) All 7 newborn states have f > 0; eta = 0.359043
Constants: B = 7*sup_f = 4.473482, R = 126, kappa = 2.7930e+06

=== 2. Brute-force exact check on random friend trees (n in [3, 150]) ===
Checked 120 random friend trees (20516 directed (V,W) moves) in 0.7s:
  Grouped vs literal ungrouped (V,W) cross-check on first 10 trees: EXACT MATCH
  max L observed = 147 (43 trees had L >= 60); max c observed = 6 (0 trees had c > 10)
  (H3) max pre-existing state changes = 5 (limit 6), max f changes = 5, max new core = 1 (limit 1), max |Delta Phi| = 0.346477 (limit B = 4.473482)
  (H4) min newborn f observed = 0.359043 (theoretical eta = 0.359043)
  (H5) max n*P(relevant state change, c<=10) = 65.700000 (limit 126), max n*P(f change, f>0) = 65.700000
       witness for max n*P(rel): n=74, v=0, state=(56, 5, (2, 1), 2, (0, 0, 0)), n*P=65.700000
  (H2) min exact drift ratio (n*E[Delta]/Phi) = 0.543316 (max = 0.587254, target delta = 0.543000); min exact slack = 1.808332e-04
  ALL CHECKS PASSED ON RANDOM TREES: True

=== 3. Brute-force exact check on extremal / corner-case trees ===
  [T_3 (path on 3)] n=3, C_n=1: max_state_chg=1, max_f_chg=1, max_new=1, max_nP_rel=3.0000, max_nP_f=3.0000, drift_ratio=0.543323, pass=True
  [Star S_4 (c=0, L=3)] n=4, C_n=1: max_state_chg=1, max_f_chg=1, max_new=1, max_nP_rel=4.0000, max_nP_f=4.0000, drift_ratio=0.543321, pass=True
  [Star S_70 (c=0, L=69 >= 60)] n=70, C_n=1: max_state_chg=1, max_f_chg=1, max_new=1, max_nP_rel=1.0000, max_nP_f=1.0000, drift_ratio=0.543321, pass=True
  [Tight b0=6: (L=0,c=5) + 5 P_1 nbrs] n=11, C_n=6: max_state_chg=6, max_f_chg=6, max_new=1, max_nP_rel=11.0000, max_nP_f=11.0000, drift_ratio=0.569217, pass=True
  [Cap L>=60: L=62, c=2] n=69, C_n=3: max_state_chg=2, max_f_chg=2, max_new=1, max_nP_rel=5.5000, max_nP_f=5.5000, drift_ratio=0.545103, pass=True
  [High R: L=59, c=10 + ten O_5 nbrs] n=110, C_n=11: max_state_chg=2, max_f_chg=2, max_new=1, max_nP_rel=102.0000, max_nP_f=102.0000, drift_ratio=0.543712, pass=True
  [High c>10 & dead (L=0,c=2,top)] n=75, C_n=14: max_state_chg=3, max_f_chg=2, max_new=1, max_nP_rel=6.5769, max_nP_f=6.5769, drift_ratio=0.559891, pass=True

=== 4. Negative controls (planted wrong claims must be caught) ===
  [PASS] Planted wrong b0=2 CAUGHT on random tree #5 (n=8): move W=1 changed 3 pre-existing core states ([1, 3, 4])
  [PASS] Planted wrong b0=5 CAUGHT on extremal tree (n=11): move W=0 changed 6 pre-existing core states ([0, 1, 2, 3, 4, 5]) -> b0=6 is sharp!
  [PASS] Planted wrong R=1 CAUGHT on random tree #0 (n=3): vertex v=0 (state=(2, 0, (0, 0), 0, (0, 0, 0))) had n*P = 3 = 3.0000 > 1
  [PASS] Planted wrong R=100 CAUGHT on extremal high-R tree (n=110): vertex v=0 (state=(59, 10, (0, 0), 0, (0, 0, 10))) had n*P = 102.0000 > 100
  [PASS] Planted wrong delta=0.600 CAUGHT on random tree #0 (n=3): n*E[Delta]/Phi = 0.543323 < 0.600000
  [PASS] Star S_200 (L=199 >= 60): raw state-change n*P = 200 > 126, but relevant/f-change n*P = 1.0000 (f-change = 1.0000) <= 126 (confirms L>=60 cap is essential for H5!)

ALL VERIFICATION CHECKS AND NEGATIVE CONTROLS PASSED in 38.8s total.
```

---

## 3. Brute-force exact rational check of (H2)

We rebuilt the exact `certs9.py` table (`certs9.Model(60, 10, 2, 6)`, `certs9.lower_table(M, 0.5433)`) in 34.7s on 1 core (1,186,328 entries; Perron eigenvalue `1.000033`), obtaining the exact same `sup f = 0.639069`, `fmin = 3.611209e-03`, `18` zero entries, `eta = 0.359043`, `B = 4.473482`, `R = 126`, and `kappa = 2.7930e+06` reported in `attempts/strogatz-fable/outputs/as_lower_check.log`.

On all 120 random friend trees (`n in [3, 150]`) and all 8 extremal/corner-case trees (including `L >= 60`, `c > 10`, and dead `L = 0` states), brute-force enumeration over all `(V, W)` moves in exact `fractions.Fraction` arithmetic verified
`n * E[Phi_{n+1} - Phi_n | T_n] >= (543/1000) * Phi_n`
with zero failures:
- Minimum exact drift ratio `n * E[Delta] / Phi` across the 120 random friend trees: `0.543316 >= 0.543000` (minimum exact slack `1.808332e-04 > 0`).
- Minimum exact drift ratio across the extremal trees: `0.543321` (on `S_4` and `S_70`, and `0.559891` on the `c > 10` + dead-state tree).
- Negative control: testing the same table against a planted wrong `delta = 600/1000` immediately fails on `T_3` (`0.543323 < 0.600000`).

(We also inspected `certs9.verify` lines 155–199, including the `L > Lmax` tail endpoints and `c > cmax` handling, and found the encoding of the worst-case drift lower bound consistent with `attempts/x-strogatz/REPORT.md` §§3.3–3.5.)

---

## 4. Does Theorem A follow? Fairness of the literature comparison

- **Does Theorem A follow?** Yes. Theorem B is valid, hypotheses (H1), (H3), (H4), and (H5) hold for the `certs9.py` potential with `delta = 543/1000` (verified both analytically and computationally above), and (H2) holds by the exact rational certificate verification of `certs9.py` (previously independently re-verified twice, re-run in `as_lower_check.log`, and brute-force spot-checked on all 128 trees in our script). Combining Theorem B with the `certs9.py` certificate yields `liminf_{n -> infinity} log C_n / log n >= 0.543` almost surely (Theorem A). Combining Theorem A with the cited almost-sure upper bound `limsup log C_n / log n <= 0.586` from `attempts/x-strogatz/REPORT.md` §3.1 indeed gives `n^(0.543 - eps) <= C_n <= n^(0.586 + eps)` a.s. for every `eps > 0` and all large `n`.
- **Fairness of the literature comparison:** Fair and accurate, with two minor contextual qualifications worth noting:
  1. In Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou, and Lugosi, *Random friend trees* (arXiv:2403.20185 / *Random Structures & Algorithms* 2026), Proposition 6.20 indeed proves that for `alpha = (sqrt(13) - 3)/2` and any `0 < delta < alpha/3 = (sqrt(13) - 3)/6 = 0.1009252...`, almost surely `n^{-delta} X_{>= 2}^{(n)} -> infinity`. Citing `(sqrt(13) - 3)/6` (note: `as_lower.md` writes `= 0.1009` instead of `≈ 0.1009`) as the best published almost-sure lower exponent for `C_n = X_{>= 2}^{(n)}` is accurate and fair.
  2. For complete context, Addario-Berry et al.'s Theorem 4.7 proves an almost-sure polynomial lower bound `n^delta` (with `delta > 0.1`) for `X_n^{>= k}` for *every* fixed `k >= 2` by a purely human proof, whereas Theorem A improves the exponent specifically for `C_n = X_n^{>= 2}` (vertices of degree `>= 2`) and relies on a computer-checked rational certificate of 1,205,767 inequalities. Neither point undercuts the claim made in `as_lower.md`, which explicitly states that Theorem A is about `C_n` (degree `>= 2`), explicitly cites Proposition 6.20, and transparently describes the computer-checked certificate.
