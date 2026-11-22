VERDICT: CORRECT

Theorem D (parts (a)–(d)) and the Corollary follow from the written proofs (taking Claim 2 of `RESULT.md` as given for the Corollary). Every lemma (Lemmas 1–4), every probabilistic conditioning and independence argument, every product identity and lower bound, every numerical constant (`1/30`, `1/4`, `1/368`, `44160`, `0.44`, `33`), and Gautschi's inequality as applied are mathematically sound and corroborated by independent exact and Monte Carlo numerical tests. Only two trivial editorial/boundary slips were found—the summation-by-parts equality on line 129 of Part (c) fails strictly at the empty-sum boundary `k = 2` (`m0 = n`), where `<=` holds and the bound is trivial; and Part (d) (lines 136–137) omits `eps < mu` alongside `eps < 1 - mu`—both admit immediate one-line repairs with no impact on any result or constant.

---

## 1. Line-by-Line Audit of `attempts/strogatz-fable/proofs/degree_tail.md`

### Preliminaries (lines 3–7, 43–45)
* **Transition probability `P(W_j = w | F_j) = (1/j) sum_{u ~ w} 1/D_u(j)`**: **VALID.**
  In `T_j` (a tree on `j >= 2` vertices, so `D_u(j) >= 1` for all `u`), `V_j` is chosen uniformly from the `j` vertices and `W_j` uniformly from the `D_{V_j}(j)` neighbours of `V_j`. Summing over `u ~ w` yields `P(W_j = w | F_j) = (1/j) sum_{u ~ w} 1/D_u(j)`. Because the random friend tree `(T_j)` is Markov with persistent vertex labels `1, ..., j`, this holds conditionally on `F_j = sigma(T_3, ..., T_j)`.

### Lemma 1 (lines 47–54)
* **First-moment one-step bound and iteration (`E[D_v(n) | F_m] <= D_v(m) * n/m`)**: **VALID.**
  At step `j`, `D_v(j+1) - D_v(j) = 1_{W_j = v}` because only `W_j` receives the new leaf and degrees never decrease. Each of the `D_v(j)` neighbours `u` of `v` has `1/D_u(j) <= 1`, so `P(W_j = v | F_j) <= D_v(j)/j`, giving `E[D_v(j+1) | F_j] <= D_v(j)(j+1)/j`. Since `(j+1)/j > 0`, iterating via the tower property over `j = m, ..., n-1` gives `prod_{j=m}^{n-1} (j+1)/j = n/m` (including empty product `1` at `n = m`). (Equality is achieved at the center of a star.)
* **Second-moment one-step bound (`Q = D(D+1)`)**: **VALID.**
  When `W_j = v`, `Q(j+1) - Q(j) = (D+1)(D+2) - D(D+1) = 2(D+1)`, and `0` otherwise. Using `P(W_j = v | F_j) <= D_v(j)/j` yields
  `E[Q(j+1) | F_j] <= Q(j) + 2(D_v(j)+1) D_v(j)/j = Q(j)(j+2)/j`.
* **Product formula (`prod_{j=m}^{n-1} (j+2)/j = n(n+1)/(m(m+1))`)**: **VALID.**
  For `n = m` both sides equal `1`; for `n > m`, numerator `(m+2)...(n+1) = (n+1)!/(m+1)!` divided by denominator `m...(n-1) = (n-1)!/(m-1)!` equals `n(n+1)/(m(m+1))`. Iterating with positive factors `(j+2)/j` preserves the inequality.

### Theorem D, Part (a) (lines 15–17, 56–60)
* **Accounting for every vertex of degree `>= k`**: **VALID.**
  Since `k >= 2`, every vertex `v` of `T_n` with `D_v(n) >= k` is a core vertex of `T_n` (`D_v(n) >= 2`). Either `v` is already core in `T_{m0}` (contributing at most `C_{m0}` vertices), or `v` is not core in `T_{m0}` (either a leaf of `T_{m0}` or a vertex created at some step `j >= m0`, which is born as a leaf of degree `1`). Since degrees increase by at most `1` per step and only for `W_m`, any such `v` has a unique promotion step `m` with `m0 <= m <= n-1` where `D_v(m) = 1` and `D_v(m+1) = 2` (i.e. `W_m = v` is a leaf of `T_m`). At each step `m`, at most one vertex is promoted, so `C_{m+1} - C_m in {0, 1}` is the exact indicator of promotion at step `m`. For `m0 = n`, the sum is empty and `N_k(n) <= C_n = C_{m0}` holds directly.
* **Legitimacy of conditioning on `F_{m+1}` and Markov/Lemma 1 application**: **VALID.**
  Both `C_{m+1} - C_m` and the promoted vertex `v = W_m` (the unique neighbour of the new vertex `m+1` in `T_{m+1}`, equivalently the unique vertex of `T_m` whose degree grew from `1` to `2`) are `F_{m+1}`-measurable. On `{C_{m+1} - C_m = 1}`, `D_v(m+1) = 2`, `m+1 >= m0 + 1 >= 4 >= 3`, and `m+1 <= n`, so Lemma 1 applies conditionally on `F_{m+1}` (future evolution is Markov in `T_{m+1}`):
  `P(D_v(n) >= k | F_{m+1}) <= E[D_v(n) | F_{m+1}]/k <= 2n/((m+1)k)`.
  Multiplying by the `F_{m+1}`-measurable nonnegative indicator `C_{m+1} - C_m`, taking expectations, and summing over `m = m0, ..., n-1` yields `E[N_k(n)] <= a_{m0} + (2n/k) sum_{m=m0}^{n-1} (a_{m+1} - a_m)/(m+1)`.

### Lemma 2 (lines 62–77)
* **Drift of `L_u`**: **VALID.**
  At step `j >= m`, `u` is core (`D_u(j) >= D_u(m) >= 3`). `L_u` increases by `1` iff `W_j = u`, which has conditional probability `(L_u(j) + rho_u(j))/j` with `rho_u(j) = sum_{w ~ u, D_w(j)>=2} 1/D_w(j) >= 0`. `L_u` decreases by `1` iff `W_j` is one of the `L_u(j)` leaf neighbours `x` of `u`; each such `x` has `u` as its unique neighbour, so `P(W_j = x | F_j) = 1/(j D_u(j))`, and summing over disjoint events gives `L_u(j)/(j D_u(j)) <= 1/j` (valid even when `L_u(j) = 0`, since `0 <= 1/j`). All other choices of `W_j` leave `L_u` unchanged. Hence
  `E[L_u(j+1) - L_u(j) | F_j] = (L_u(j) + rho_u(j))/j - L_u(j)/(j D_u(j)) >= (L_u(j) - 1)/j`,
  which holds even when `L_u(j) = 0` (`rho_u(j)/j >= 0 > -1/j`).
* **Iteration when `L_u - 1` can be negative**: **VALID.**
  Rewriting gives `E[L_u(j+1) - 1 | F_j] >= (L_u(j) - 1)(j+1)/j`. Monotonicity of conditional expectation (`A >= B => E[A | F_m] >= E[B | F_m]`) holds regardless of sign, and multiplying by the positive deterministic factors `(j+1)/j > 0` preserves `>=` regardless of sign. Iterating via the tower property from `m` to `n` yields `E[L_u(n) - 1 | F_m] >= (l - 1) n/m`, so `E[L_u(n) | F_m] >= 1 + (l - 1)n/m > (l - 1)n/m`.
* **Second-moment bound**: **VALID.**
  Since `0 <= L_u(n) <= D_u(n)`, we have `L_u(n)^2 <= D_u(n)(D_u(n)+1)`. By Lemma 1, `D_u(m) = l + c_u(m) <= l + 2`, and `(n+1)/(m+1) <= n/m` (equivalent to `m <= n` for `n >= m >= 1`), we get
  `E[L_u(n)^2 | F_m] <= (l+2)(l+3) n(n+1)/(m(m+1)) <= (l+2)(l+3) (n/m)^2`.
* **Paley-Zygmund inequality**: **VALID.**
  `L_u(n)` is nonnegative, bounded (`<= n`), and has positive conditional mean `E[L_u(n) | F_m] >= (l-1)n/m >= 2 > 0`. Applying Paley-Zygmund with `theta = 1/2` gives
  `P(L_u(n) >= E[L_u(n) | F_m]/2 | F_m) >= (1/4) E[L_u(n) | F_m]^2 / E[L_u(n)^2 | F_m] >= (l-1)^2 / (4(l+2)(l+3))`.
  Because `E[L_u(n) | F_m]/2 >= (l-1)n/(2m)`, the event `{L_u(n) >= E[L_u(n) | F_m]/2}` is contained in `{L_u(n) >= (l-1)n/(2m)}`, yielding the stated inequality.
* **Monotonicity and constant `1/30`**: **VALID.**
  Writing `(l-1)^2 / (4(l+2)(l+3)) = (1/4) * (1 - 3/(l+2)) * (1 - 4/(l+3))`, both factors are positive and strictly increasing for `l >= 3`. At `l = 3`, the value is `2^2 / (4 * 5 * 6) = 4/120 = 1/30`.

### Lemma 3 (lines 79–90)
* **Case `C = 1`**: **VALID.**
  If a tree has exactly one core vertex, all other vertices are leaves, so the core vertex has `c = 0` core neighbours and `L = D >= 2 >= 1` leaf neighbours. Hence it satisfies `(c <= 2 and L >= 1)`, giving `E = 1 >= 1/4 = C/4`.
* **Counting argument for `C >= 2`**: **VALID.**
  Any path in a tree between two core vertices passes only through vertices of degree `>= 2` (hence core), so the core vertices induce a connected subtree `T_core` of `C` vertices and `C - 1` edges; every core vertex has `c >= 1`.
  - Vertices with `c = 1` have `L = D - 1 >= 1`, so all `n_1` of them are eligible (`E >= n_1`).
  - Vertices with `c = 2` are eligible if `L >= 1` or if `(L, c) = (0, 2)` and adjacent to another `(0, 2)` vertex. Thus the non-eligible vertices are precisely `{c >= 3}` (size `n_3`) disjoint union `I` (vertices with `(L, c) = (0, 2)` having no `(0, 2)` neighbour), giving `E = C - n_3 - I`.
  - The `I` vertices form an independent set in `T_core`, each incident to `2` edges of `T_core`; no edge of `T_core` has both endpoints in `I`, so the incident edge sets are pairwise disjoint, giving `2I <= |E(T_core)| = C - 1`.
  - In any tree on `C >= 2` vertices, `sum (c_v - 2) = -2`, so `n_1 = 2 + sum_{c >= 3}(c - 2) >= n_3 + 2`, i.e. `n_3 <= n_1 - 2`.
  - Adding `E >= C - n_3 - (C-1)/2 = (C+1)/2 - n_3` and `E >= n_1 >= n_3 + 2` cancels `n_3` and gives `2E >= (C+1)/2 + 2`, i.e., `E >= (C+5)/4 > C/4`. (Furthermore, the constant `1/4` is asymptotically sharp on subdivided 3-regular trees, where `C = 4k-5` and `E = k = (C+5)/4`.)

### Lemma 4 (lines 92–119)
* **Independence and law of `V_j`**: **VALID.**
  With vertex labels `{1, ..., j}` for `T_j`, each `V_j` (`j >= m'`) is drawn uniformly from `{1, ..., j}` independently of the entire prior history `H_j = sigma(T_2, V_2, W_2, ..., V_{j-1}, W_{j-1}) supseteq F_j`. Hence `(V_j)_{j >= m'}` are mutually independent, `V_j ~ Unif({1, ..., j})`, and jointly independent of `F_{m'}`.
* **Event `G` and what it guarantees**: **VALID.**
  Conditioned on `F_{m'}`, `u in {1, ..., m'}` is fixed. For `G = {V_j != u for all m' <= j < 8m'}`,
  `P(G | F_{m'}) = prod_{j=m'}^{8m'-1} (j-1)/j = (m'-1)/(8m'-1)`.
  Since `(m'-1)/(8m'-1) = 1/8 - (7/8)/(8m'-1)` is strictly increasing in `m' >= 1`, its minimum for `m' >= 3` is `(3-1)/(24-1) = 2/23`.
  On `G`, no leaf `x` of `u` (existing or newly created during `[m', 8m')`) can ever receive a leaf during `[m', 8m')`, because `u` is `x`'s unique neighbour so `W_j = x` requires `V_j = u`. Thus every leaf of `u` stays a leaf of `u`, `L_u` never decreases, and `c_u` cannot increase (nor decrease, since degrees are non-decreasing), so `c_u(8m') = c_u(m')`. Because `G` is a product event in the independent coordinates `(V_j)_{j=m'}^{8m'-1}`, conditioning on `G` preserves mutual independence of `(V_j)_{j=m'}^{8m'-1}` and makes `V_j ~ Unif({1, ..., j} \ {u})`.
* **Case A exact products and lower bounds (`2/3`, `5/8`, `5/138`)**: **VALID.**
  Fix an `F_{m'}`-measurable leaf `x` of `u` in `T_{m'}`. On `G`, `x` remains a leaf of `u`, so every `j` with `V_j = x` forces `W_j = u` and adds a new leaf to `u`. Given `(F_{m'}, G)`:
  - `P(exists j in [m', 3m') : V_j = x | F_{m'}, G) = 1 - prod_{j=m'}^{3m'-1} (j-2)/(j-1) = 1 - (m'-2)/(3m'-2) = 2m'/(3m'-2) > 2/3`.
  - `P(exists j in [3m', 8m') : V_j = x | F_{m'}, G) = 1 - prod_{j=3m'}^{8m'-1} (j-2)/(j-1) = 1 - (3m'-2)/(8m'-2) = 5m'/(8m'-2) > 5/8`.
  The windows `[m', 3m')` and `[3m', 8m')` are disjoint, so the two events are independent given `(F_{m'}, G)`. On their intersection with `G`, `L_u(8m') >= 1 + 1 + 1 = 3` and `c_u(8m') = c_u(m') <= 2`. Product bound: `(2/23)(2/3)(5/8) = 5/138 > 1/368`.
* **Case B exact product for `H_1` and lower bound `1/4`**: **VALID.**
  Here `u, w, z` are three distinct vertices of `T_{m'}` (`w` has exact neighbours `u` and `z` since `(L_w, c_w) = (0, 2)`). Given `(F_{m'}, G)`, `H_1 = {forall j in [m', 2m'), V_j != z} \ {forall j in [m', 2m'), V_j notin {z, w}}` has probability
  `prod_{j=m'}^{2m'-1} (j-2)/(j-1) - prod_{j=m'}^{2m'-1} (j-3)/(j-1)`.
  Since `(j-3)/(j-1) = ((j-3)/(j-2)) * ((j-2)/(j-1))`, the two telescoping products equal `(m'-2)/(2m'-2)` and `((m'-2)/(2m'-2)) * ((m'-3)/(2m'-3))` (valid for all `m' >= 3`, including `m' = 3` where the second product is `0`). Subtracting gives
  `((m'-2)/(2m'-2)) * (1 - (m'-3)/(2m'-3)) = (m'-2)m' / ((2m'-2)(2m'-3))`.
  Furthermore `4(m'-2)m' - (2m'-2)(2m'-3) = 2m' - 6 >= 0` for all `m' >= 3`, proving `>= 1/4` (with exact equality at `m' = 3` and limit `1/4` as `m' -> infinity`).
* **Case B: why `w` still has degree 2 at time `tau`**: **VALID.**
  On `G cap H_1`, `tau = min{j >= m' : V_j = w}` satisfies `m' <= tau < 2m'`. By induction on `j in [m', tau]`, if `w` has only neighbours `u` and `z` in `T_j` (`j < tau`), then `W_j = w` requires `V_j in {u, z}`, which is excluded by `G` (`V_j != u`) and `H_1` (`V_j != z`). Hence no leaf attaches to `w` before `tau`, so `w` still has exact neighbours `{u, z}` in `T_tau`.
* **Case B: the factor `1/2`**: **VALID.**
  The events `G`, `H_1`, and the time `tau` depend solely on `(F_{m'}, (V_j)_{j=m'}^{8m'-1})`. Given `(F_{m'}, (V_j)_{j=m'}^{8m'-1})` on `G cap H_1`, the uniform choice of `W_tau` among the two neighbours `{u, z}` of `V_tau = w` in `T_tau` uses fresh randomness independent of the entire sequence `(V_j)` and earlier `W_j`, so `P(W_tau = u | F_{m'}, G, H_1) = 1/2`.
* **Case B: identity of the new leaf `x`, independence of later hits, and bounds `1/2`, `1/368`**: **VALID.**
  On `G cap H_1 cap {W_tau = u}`, the new vertex attached to `u` at step `tau` has label `x = tau + 1 in {m'+1, ..., 2m'}` (so `x != u`, `x` is already in `T_{2m'}` even if `tau = 2m'-1`, and `x` remains a leaf of `u` on `G` through time `8m'`). Conditioned on `F_{m'}` and `G`, `G` splits as `G_{<2m'} cap G_{[2m',4m')} cap G_{[4m',8m')}` across the mutually independent blocks `H_{2m'} = sigma(T_2, V_2, W_2, ..., V_{2m'-1}, W_{2m'-1})`, `(V_j)_{j=2m'}^{4m'-1}`, and `(V_j)_{j=4m'}^{8m'-1}`. On each atom `{tau = t} cap H_1 cap {W_t = u} in H_{2m'}` (`m' <= t <= 2m'-1`), `x = t+1` is fixed, so the events `A_1 = {exists j in [2m', 4m') : V_j = x}` and `A_2 = {exists j in [4m', 8m') : V_j = x}` are conditionally independent with probabilities
  `p_1 = 1 - (2m'-2)/(4m'-2) = 2m'/(4m'-2) > 1/2`, `p_2 = 1 - (4m'-2)/(8m'-2) = 4m'/(8m'-2) > 1/2`.
  Because `p_1` and `p_2` do not depend on `t`, summing over `t = m', ..., 2m'-1` yields exact factorization `P(H_1 cap {W_tau = u} cap A_1 cap A_2 | F_{m'}, G) = P(H_1 cap {W_tau = u} | F_{m'}, G) * p_1 * p_2`. (This is the precise sense of the text's slightly compressed phrase "independent of each other and of what happened before `2m'`".) On the intersection, `u` gains `x` at `tau` plus at least two more leaves in `[2m', 4m')` and `[4m', 8m')`, so `L_u(8m') >= 3` and `c_u(8m') = 2`. Product bound: `(2/23)(1/4)(1/2)(1/2)(1/2) = 2/(23*32) = 1/368`.

### Theorem D, Part (b) (lines 19–21, 121–124)
* **Use of Lemma 4 and Lemma 2 at time `m = 8m'`**: **VALID.**
  Given `m' >= 3`, `k >= 2`, `n >= 8m'k`: we have `m = 8m' >= 24 >= 3` and `n >= 8m'k >= 16m' > m`. Let `A_u = {L_u(8m') >= 3 and c_u(8m') <= 2} in F_{8m'}`. For any eligible vertex `u` of `T_{m'}`, Lemma 4 gives `P(A_u | F_{m'}) >= 1/368`. On `A_u`, Lemma 2 at `m = 8m'` with `l = L_u(8m') >= 3` satisfies `(l-1)n/(2m) >= 2n/(16m') = n/(8m') >= k`, so `P(L_u(n) >= k | F_{8m'}) >= 1/30`. By the tower property (`F_{m'} subseteq F_{8m'}`) and `D_u(n) >= L_u(n)`, `P(D_u(n) >= k | F_{m'}) >= (1/368)(1/30)`.
* **Summing over eligible vertices and constant `44160`**: **VALID.**
  For `m' >= 3`, `C_{m'} >= C_3 = 1` a.s., so Lemma 3 gives at least `C_{m'}/4` eligible vertices in `T_{m'}`. Summing `1_{D_u(n) >= k}` over eligible vertices of `T_{m'}` and taking expectations gives `E[N_k(n)] >= E[C_{m'}]/(4 * 368 * 30) = a_{m'}/44160` (`4 * 368 * 30 = 120 * 368 = 44160`).

### Theorem D, Part (c) (lines 23–26, 126–134)
* **Choice of `m0 = ceil(2n/k)`**: **VALID.**
  For `k >= 2` and `n >= 24k`, we have `48 <= 2n/k <= n`, so `m0 = ceil(2n/k)` is an integer satisfying `3 <= 48 <= m0 <= n`.
* **Summation by parts on line 129**: **NEEDS REPAIR (boundary slip at `k = 2`, where `m0 = n`).**
  Line 129 writes:
  `sum_{m=m0}^{n-1} (a_{m+1} - a_m) b_{m+1} = a_n b_n - a_{m0} b_{m0+1} + sum_{j=m0+1}^{n-1} a_j (b_j - b_{j+1})`.
  - For `m0 <= n - 1` (which covers all `k >= 3`, since `n >= 24k >= 72 => ceil(2n/k) <= 2n/3 + 1 < n`): this identity is **VALID** (expanding `sum_{j=m0+1}^n a_j b_j - sum_{j=m0}^{n-1} a_j b_{j+1}`).
  - For `m0 = n` (which occurs whenever `k = 2`, explicitly included in `k >= 2`): the LHS is the empty sum `0`, while the RHS evaluates to `a_n b_n - a_n b_{n+1} = a_n / (n(n+1)) > 0`. Thus the `=` sign is strictly false at `k = 2`, although `LHS <= RHS` still holds.
  - **Repair**: Either note that for `k = 2` (`m0 = n`) the sum in (a) is empty (`= 0 <= A m0^(lambda-1)(1 + 1/(1-lambda))`, or directly `E[N_2(n)] = a_n <= A n^lambda <= A(2 + 1/(1-lambda))(3n/2)^lambda`) and apply summation by parts for `m0 <= n - 1`, or replace `=` on line 129 with `<=`.
* **Upper-bound inequalities and constant `A(2 + 1/(1-lambda))(3n/k)^lambda` (lines 130–133)**: **VALID.**
  - `a_n b_n <= A n^(lambda-1) <= A m0^(lambda-1)` (since `m0 <= n` and `lambda - 1 < 0`); `-a_{m0} b_{m0+1} < 0` is safely dropped.
  - `a_j (b_j - b_{j+1}) = a_j / (j(j+1)) <= A j^lambda / j^2 = A j^(lambda-2)`.
  - Since `x -> x^(lambda-2)` is strictly decreasing on `(0, infinity)` (`lambda < 1`), `sum_{j > m0} j^(lambda-2) < int_{m0}^infinity x^(lambda-2) dx = m0^(lambda-1)/(1-lambda)`.
  - Combining with `a_{m0} <= A m0^lambda` and `2n/k <= m0`:
    `E[N_k(n)] <= A m0^lambda + m0 * A m0^(lambda-1)(1 + 1/(1-lambda)) = A(2 + 1/(1-lambda)) m0^lambda`.
  - Since `n/k >= 24 >= 1`, `m0 = ceil(2n/k) <= 2n/k + 1 <= 3n/k`, giving `m0^lambda <= (3n/k)^lambda`.
* **Lower-bound choice of `m'` and constant `(a/44160)(n/(12k))^delta` (lines 133–134)**: **VALID.**
  Take `m' = floor(n/(8k))`. Since `n >= 24k`, `n/(8k) >= 3`, so `m' >= 3` and `8m'k <= n`, satisfying (b). Also `m' > n/(8k) - 1 >= n/(8k) - (1/3)(n/(8k)) = (2/3)(n/(8k)) = n/(12k)`. Since `delta > 0`, `a_{m'}/44160 >= (a/44160)(n/(12k))^delta`.

### Theorem D, Part (d) (lines 28–29, 136–138)
* **Proof of (d)**: **NEEDS REPAIR (minor wording: restrict `eps < mu` as well as `eps < 1 - mu`).**
  Line 136–137 says "For every `eps > 0` ... Apply (c) with `delta = mu - eps` and `lambda = mu + eps < 1`". To invoke (c), one needs `0 < delta <= lambda < 1`, so one should state `0 < eps < min(mu, 1 - mu)`. With that trivial wording fix, everything else is **VALID**: since `C_m >= C_3 = 1` a.s. for all `m >= 3`, `1 <= a_m <= m` for all `m >= 3`, so `log a_m / log m -> mu` implies positive constants `a_eps, A_eps` exist with `a_eps m^(mu-eps) <= a_m <= A_eps m^(mu+eps)` for all `m >= 3`. Applying (c) for `n/k >= 24`, taking logarithms, and dividing by `log(n/k)` yields error terms `O_eps(1/log(n/k))` that depend only on `eps` and `n/k` (not on `k` separately), proving uniform convergence in `k >= 2` as `n/k -> infinity` followed by `eps -> 0`.

### Corollary (lines 31–36)
* **Gautschi's inequality `n(n+1)^(s-1) < Gamma(n+s)/Gamma(n) < n^s` for `0 < s < 1`**: **VALID.**
  Standard Gautschi's inequality states `n^(1-s) < Gamma(n+1)/Gamma(n+s) < (n+1)^(1-s)` for `n > 0` and `0 < s < 1`. Substituting `Gamma(n+1) = n Gamma(n)`, taking reciprocals, and multiplying by `n` gives `n(n+1)^(s-1) < Gamma(n+s)/Gamma(n) < n^s`.
* **Numerical constants `0.44` and `33`**: **VALID.**
  - Lower bound: for `m >= 3`, `0.5022 * Gamma(m+0.55)/Gamma(m) > 0.5022 * m * (m+1)^(-0.45) = 0.5022 * (m/(m+1))^0.45 * m^0.55 >= 0.5022 * (3/4)^0.45 * m^0.55 = 0.4412190721... * m^0.55 >= 0.44 m^0.55`. (Matches text's `0.4412`.)
  - Upper bound: `3.974 * Gamma(m+0.58)/Gamma(m) < 3.974 * m^0.58` for `m >= 3`. Applying Theorem D(c) with `a = 0.44, A = 3.974, delta = 0.55, lambda = 0.58` gives lower bound `(0.44/44160)(n/(12k))^0.55` and upper prefactor `3.974 * (2 + 1/0.42) * 3^0.58 = 32.9250712252... <= 33`. (Matches text's `32.93`.)

---

## 2. Numerical Tests (Single CPU Core, with Planted Errors)

All scripts were written in Python 3 under `checks/fable-degree-tail-gemini/` and executed on a single CPU core (`taskset -c 0` with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`) using `attempts/strogatz-fable/.venv/bin/python`.

### Commands Run
```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 taskset -c 0 bash -c '
PY=attempts/strogatz-fable/.venv/bin/python
DIR=checks/fable-degree-tail-gemini
echo "########## TEST 1: LEMMA 3 ##########"
$PY $DIR/test1_lemma3.py
echo ""
echo "########## TEST 2: LEMMA 1 ##########"
$PY $DIR/test2_lemma1.py
echo ""
echo "########## TEST 3: LEMMA 2 ##########"
$PY $DIR/test3_lemma2.py
echo ""
echo "########## TEST 4: LEMMA 4 ##########"
$PY $DIR/test4_lemma4.py
'

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 taskset -c 0 \
  attempts/strogatz-fable/.venv/bin/python \
  checks/fable-degree-tail-gemini/test4b_caseB_subevent.py

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 taskset -c 0 \
  attempts/strogatz-fable/.venv/bin/python \
  checks/fable-degree-tail-gemini/test5_parts_ab_cor.py
```

### Test (i): Lemma 3 on all trees with `3 <= n <= 13` (`test1_lemma3.py`)
* Exhaustively enumerated all **2,286** non-isomorphic unrooted trees on `3 <= n <= 13` vertices (counts per `n` matched OEIS A000055: `1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301`; every tree with `n >= 3` is reachable by the random friend tree process).
* Every tree satisfied `E >= C/4`, core-spans-a-tree, `2I <= C - 1`, `n_3 <= n_1 - 2`, and `E >= (C+5)/4` (for `C >= 2`).
* **Tight cases**:
  - Smallest ratio `E/C` for `n <= 13`: **`3/7 = 0.428571`**, achieved by 14 trees with `C = 7, E = 3` (subdivided `K_{1,3}` core: one `c=3` center, three `(0,2)` subdivisions, three `(L>=1, c=1)` core leaves, starting at `n = 10` with 3 leaves).
  - Smallest additive margin `E - C/4`: **`3/4`** (the 11 stars `K_{1, n-1}` for `n = 3..13`, where `C = 1, E = 1`).
  - Exactly 39 trees with `C >= 2` achieve equality `E = (C+5)/4` (`25` with `C = 3, E = 2` and `14` with `C = 7, E = 3`). Extending the subdivided-3-regular-tree family confirms `E/C = (C+5)/(4C) -> 1/4` (`k = 50 => n = 245, C = 195, E = 50, E/C = 0.256410`).
* **Planted error**: Omitting the Case-B eligibility clause (`(0,2)` with a `(0,2)` neighbour) caused **7 trees** on `n <= 13` to fail `E >= C/4` (first failure: path `P_11` with `C = 9, E_broken = 2 < 2.25`).

Real output excerpt:
```text
Total non-isomorphic trees checked (3 <= n <= 13): 2286 (matches sum of OEIS A000055)
ALL 2286 trees satisfy E >= C/4, core-tree property, 2I <= C-1, n3 <= n1-2, and E >= (C+5)/4 (for C>=2).
Minimum ratio E/C over 3 <= n <= 13: 3/7 (0.428571), achieved by 14 tree(s):
  n=10, C=7, E=3, n1=3, n3=1, I=3, core (L,c)=((0, 2), (0, 2), (0, 2), (0, 3), (1, 1), (1, 1), (1, 1))
...
Minimum margin E - C/4 over 3 <= n <= 13: 3/4 (0.750000), achieved by 11 tree(s):
  n=3, C=1, E=1, deg_seq=(2, 1, 1), core (L,c)=((2, 0),)
...
Total trees with C >= 2 achieving exact proof bound E = (C+5)/4: 39
  Counts by C: {3: 25, 7: 14}
=== Asymptotically sharp subdivided-binary-tree family (E = (C+5)/4 -> C/4) ===
  k= 2: n=  5, C=  3, E= 2, E/C=0.666667 (exact (C+5)/(4C))
  k= 3: n= 10, C=  7, E= 3, E/C=0.428571 (exact (C+5)/(4C))
  k= 4: n= 15, C= 11, E= 4, E/C=0.363636 (exact (C+5)/(4C))
  k= 5: n= 20, C= 15, E= 5, E/C=0.333333 (exact (C+5)/(4C))
  k=10: n= 45, C= 35, E=10, E/C=0.285714 (exact (C+5)/(4C))
  k=50: n=245, C=195, E=50, E/C=0.256410 (exact (C+5)/(4C))
=== Planted deliberate error: drop Case-B eligibility clause ===
Planted error detected! 7 trees with 3 <= n <= 13 FAIL E >= C/4 when Case-B clause is omitted.
  FAIL example: n=11, C=9, E_broken=2 < C/4=2.25, deg_seq=(2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1)
```

### Test (ii): Lemma 1 exact one-step enumeration (`test2_lemma1.py`)
* Checked both product identities exactly with `fractions.Fraction` for all `3 <= m <= n <= 64`.
* Checked exact one-step enumeration over all `(V_j, W_j)` pairs on **98 trees** (simulated random friend trees plus stars and paths across `j in {3, 4, 5, 10, 20, 50, 100, 300, 500}`, totaling **7,904 vertices**): **0 violations** of `E[D_v(j+1) | T_j] <= D_v(j)(j+1)/j` or `E[Q_v(j+1) | T_j] <= Q_v(j)(j+2)/j`, with exact tightness `1.0000000000` attained on star centers.
* **Planted error**: Replacing `(j+1)/j` by `(j+0.9)/j` and `(j+2)/j` by `(j+1.9)/j` produced immediate violations on both the star (`j=20`, ratio `1.004785` / `1.004566`) and simulated trees (`j=10`: `D:1, Q:1`; `j=50`: `D:1, Q:0`; `j=300`: `D:5, Q:2`).

Real output excerpt:
```text
Exact product formulas verified for all 3 <= m <= n <= 64.
  j=  3: 12 trees (   36 vertices), violations=(D:0, Q:0), max LHS/RHS=(D:1.0000000000, Q:1.0000000000)
  ...
  j=500:  7 trees ( 3500 vertices), violations=(D:0, Q:0), max LHS/RHS=(D:1.0000000000, Q:1.0000000000)
PASS: Checked 98 trees (7904 total vertices), 0 violations. Tightness max LHS/RHS = (1.0000000000, 1.0000000000) [exact 1.0 on star center].
=== Planted deliberate error: tighten factors (j+1)/j -> (j+0.9)/j and (j+2)/j -> (j+1.9)/j ===
  On star j=20: violations=(D:1, Q:1), max LHS/RHS=(D:1.004785, Q:1.004566)
  On simulated tree j= 10: violations=(D:1, Q:1), max LHS/RHS=(D:1.009174, Q:1.008403)
  On simulated tree j= 50: violations=(D:1, Q:0), max LHS/RHS=(D:1.000926, Q:0.999889)
  On simulated tree j=300: violations=(D:5, Q:2), max LHS/RHS=(D:1.000270, Q:1.000207)
```

### Test (iii): Lemma 2 Monte Carlo (`test3_lemma2.py`)
* Note on realizable sizes: `c = 0, l = 3` uniquely determines the 4-vertex star `K_{1,3}` (`m = 4`), so we tested `c = 0, l = 3` at `m = 4` and `c = 0, l = 4` at `m = 5`; `c = 2, l = 3` requires `m >= 8`.
* Checked exact one-step drift `E[dL_u | F_j] >= (L_u - 1)/j` including `L_u = 0` (`+0.200000 >= -0.200000`).
* Across **38** `(config, m, l, c, n/m)` combinations covering `c in {0, 1, 2}`, diverse surroundings (giant star hubs minimizing `rho_u`, paths, binary trees, asymmetric hubs), `m` from `4/5` to `200`, and `n/m` from `2` to `50` (`600` to `2500` trials each): every case exceeded `1/30 = 0.0333`, with smallest empirical `P_hat = 0.8117 +- 0.0160` (`c1_star, m=200, n/m=50`, ~`24.4x` the bound), and all verified `E[L] >= (l-1)n/m` and `E[L^2] <= (l+2)(l+3)(n/m)^2`.
* **Planted errors**:
  - (A) Violating `c_u <= 2` (`c = 20, l = 3, m = 50, n/m = 20`): empirical `E[L^2] = 41314.6 > 12000.0 = (l+2)(l+3)(n/m)^2` (`3.44x`, caught).
  - (B) Corrupted leaf-hit transition (`c = 2, m = 25, n/m = 20`): `P_hat = 0.0000 +- 0.0000 < 1/30` (caught).

Real output excerpt:
```text
config                m  l  c  n/m trials    P_hat      SE  E[L]/lo E[L2]/hi
----------------------------------------------------------------------------
c0_star               4  3  0    2   2500   0.9688  0.0035    1.280    0.247
c0_star               4  3  0   50   2500   0.8820  0.0065    1.218    0.234
c0_star               5  4  0   50   2500   0.9280  0.0052    1.135    0.306
c1_star               6  3  1   50   2500   0.9280  0.0052    1.424    0.321
c1_path               6  3  1   50   2500   0.9252  0.0053    1.415    0.318
c1_binary            25  3  1   20   1500   0.8733  0.0086    1.362    0.337
c1_star             200  3  1   50    600   0.8117  0.0160    1.164    0.262
c2_two_stars          8  3  2   50   2500   0.9440  0.0046    1.617    0.418
c2_paths              8  3  2   50   2500   0.9564  0.0041    1.639    0.425
c2_asym_star_leaf    25  3  2   20   1500   0.9093  0.0074    1.462    0.375
c2_two_stars        200  3  2   50    600   0.8250  0.0155    1.252    0.304
c2_paths            200  3  2   50    600   0.9383  0.0098    1.700    0.505
PASS: All 38 Lemma 2 cases exceed 1/30 = 0.0333; smallest empirical P_hat = 0.8117 (about 24.4x the bound).
=== Planted deliberate errors for Test (iii) ===
  Planted error A (violate c<=2: c=20, m=50, n/m=20): E[L^2] = 41314.6 > (l+2)(l+3)(n/m)^2 = 12000.0 (ratio = 3.44) -> CAUGHT!
  Planted error B (broken leaf-hit rule, c=2, m=25, n/m=20): P_hat = 0.0000 +- 0.0000 < 1/30 = 0.0333 -> CAUGHT!
```

### Test (iv): Lemma 4 exact products & Monte Carlo (`test4_lemma4.py` and `test4b_caseB_subevent.py`)
* **Exact products (`m' = 3..200`, `Fraction` arithmetic)**: Every product matched its closed form identically and satisfied every stated lower bound (`2/23`, `2/3`, `5/8`, `1/4`, `1/2`, `5/138`, `1/368`).
* **Monte Carlo activation probability**: Note that the smallest tree containing a Case B vertex is `P_6` (`m' = 6`, since two adjacent `(0,2)` core vertices require two other core neighbours plus at least two leaves). Tested Case A at `m' in {3, 4, 10, 50}` and Case B at `m' in {6, 10, 50}` across paths, stars, worst-case `(L=1, c=2)` giant-star neighbours, and double-star bridges (`3000–5000` trials each). Empirical activation probabilities ranged from `0.1640 +- 0.0068` (`60.4x` above `1/368`) to `0.6996 +- 0.0065` (`257.5x` above `1/368`).
* **Bonus exactness check of Case B sub-event (`test4b_caseB_subevent.py`, `m' = 10`, `40,000` trials)**: Verified directly that `deg(w) == 2` at `tau` had `0` mismatches, the sub-event implied activation with `0` mismatches, and every conditional factor (`P(G)`, `P(H_1|G)`, `P(W_tau=u|G,H_1) = 0.49918 +- 0.01434 vs 0.50000`, `P(A_1 & A_2|G,H_1,W_tau=u) = 0.26194 +- 0.01785 vs 0.26991`, and `P(all) = 0.00398 +- 0.00031 vs 0.00402`) agreed with exact theory (`|z| <= 1.10`).
* **Planted errors**:
  - (1) Tightening `H_1`'s bound `1/4 -> 0.26`: `190` failures in `m' = 3..200` (first at `m' = 3`, value `0.25`; last at `m' = 200`, value `0.250623`).
  - (2) Off-by-one upper index in `P(G)` (`8m'` inclusive instead of `8m'-1`) at `m' = 3`: gives `1/12 = 0.083333 < 2/23 = 0.086957`.
  - (3) Running activation on an ineligible `c_u = 3` vertex (`m' = 10`): `P_hat = 0.000000 < 1/368`.

Real output excerpt:
```text
  P(G)      >= 2/23  (0.086957): min over m'=3..200 is 2/23 (0.086957) at m'=3
  Case A W1 >= 2/3   (0.666667): min over m'=3..200 is 200/299 (0.668896) at m'=200
  Case A W2 >= 5/8   (0.625000): min over m'=3..200 is 500/799 (0.625782) at m'=200
  Case B H1 >= 1/4   (0.250000): min over m'=3..200 is 1/4 (0.250000) at m'=3
  Case B W1 >= 1/2   (0.500000): min over m'=3..200 is 200/399 (0.501253) at m'=200
  Case B W2 >= 1/2   (0.500000): min over m'=3..200 is 400/799 (0.500626) at m'=200
  Case A prod >= 5/138 (0.036232): min is 0.050819 at m'=3
  Case B prod >= 1/368 (0.002717): min is 0.003557 at m'=3

Case                             m' trials    P_hat      SE  ratio to 1/368
----------------------------------------------------------------------------
CaseA m'=3  P_3 center            3   5000   0.6996  0.0065          257.5x
CaseA m'=4  K_{1,3} center        4   5000   0.6888  0.0065          253.5x
CaseA m'=4  P_4 core(1,1)         4   5000   0.4758  0.0071          175.1x
CaseA m'=10 path core(1,1)       10   4000   0.4963  0.0079          182.6x
CaseA m'=10 worst (L=1,c=2)      10   4000   0.2465  0.0068           90.7x
CaseA m'=50 path core(1,1)       50   3000   0.4910  0.0091          180.7x
CaseA m'=50 worst (L=1,c=2)      50   3000   0.2363  0.0078           87.0x
CaseB m'=6  P_6                   6   5000   0.2198  0.0059           80.9x
CaseB m'=10 P_10 mid             10   4000   0.2300  0.0067           84.6x
CaseB m'=10 double-star br.      10   4000   0.1855  0.0061           68.3x
CaseB m'=50 P_50 mid             50   3000   0.2247  0.0076           82.7x
CaseB m'=50 double-star br.      50   3000   0.1640  0.0068           60.4x

=== Lemma 4 Case B sub-event exactness check (m'=10, 40000 trials) ===
  P(G)                    : emp=0.11315 +- 0.00158  exact=0.11392  z=-0.49
  P(H1 | G)               : emp=0.26867 +- 0.00659  exact=0.26144  z=+1.10
  P(W_tau=u | G, H1)      : emp=0.49918 +- 0.01434  exact=0.50000  z=-0.06
  P(A1&A2 | G,H1,W_tau=u) : emp=0.26194 +- 0.01785  exact=0.26991  z=-0.45
  P(all)                  : emp=0.00398 +- 0.00031  exact=0.00402  z=-0.14
  deg(w)==2 at tau mismatches: 0; activation implication mismatches: 0
```

### Test (v): Parts (a), (b), (c), and Corollary on simulated data up to `n = 2e5` (`test5_parts_ab_cor.py`)
* Verified Gautschi's inequality via `math.lgamma` for `s in {0.55, 0.58}` across `n = 1..200,000`, and Corollary constants (`0.4412190721 > 0.44`, `32.9250712252 < 33`).
* Simulated **120 independent runs** up to `n = 200,000` (15.83s on 1 core), tracking `C_m` for every `m = 3..200,000` and `N_k(n)` at `n in {20,000, 50,000, 200,000}` for `k in {2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000}`.
* Every `k` at every checkpoint satisfied the Part (a) upper bound (with `m0 = ceil(2n/k)`, ratio `UB/E[Nk]` ranging cleanly from `1.00x` at `k = 2` to `4.65x` at `k = 5000`), the Part (b) lower bound (with `m' = floor(n/(8k))`, `E[Nk]/LB ~ 1.1e5`), and both Corollary bounds.
* Confirmed numerically the line 129 summation-by-parts boundary slip at `k = 2` (`m0 = n = 200000`): `LHS = 0.0`, `line-129 RHS = 1.8778656107e-08 != 0`, whereas at `k = 3, 10` `|RHS - LHS| <= 2e-18`.
* **Planted errors**:
  - (A1) Omitting `a_{m0}` in Part (a): caught at 3 values of `k` (e.g. `k = 2`: `E[Nk] = 751.15 > 0.00`).
  - (A2) Scaling Part (a) UB by `0.2x`: caught at `8/12` values of `k` (e.g. `k = 2`: `751.15 > 150.23`).
  - (B) Replacing Part (b)'s `1/44160` factor by `10` (`10 * a_{m'} <= E[Nk]`): caught at `k = 10` (`175.58 < 635.75`).

Real output excerpt (`n = 200,000`):
```text
    k      m0     m' |     E[Nk]     SE |  Part(a)UB  UB/ENk |  Part(b)LB    ENk/LB |     CorLB      CorUB
----------------------------------------------------------------------------------------------------------
    2  200000  12500 |   751.150 32.144 |    751.150    1.00 |   0.003565  210730.8 |  0.001428   26212.83
    3  133334   8333 |   462.517 19.702 |    722.394    1.56 |   0.002839  162909.2 |  0.001143   20719.58
    5   80000   5000 |   293.075 12.460 |    637.741    2.18 |   0.002126  137866.2 |  0.000863   15406.67
   10   40000   2500 |   175.583  7.560 |    499.063    2.84 |   0.001440  121962.4 |  0.000589   10306.50
   20   20000   1250 |   111.700  4.838 |    371.701    3.33 |   0.000968  115339.2 |  0.000403    6894.68
   50    8000    500 |    63.100  2.700 |    241.156    3.82 |   0.000578  109238.7 |  0.000243    4052.37
  100    4000    250 |    41.642  1.824 |    170.125    4.09 |   0.000381  109403.8 |  0.000166    2710.89
  200    2000    125 |    28.217  1.233 |    118.559    4.20 |   0.000254  111254.3 |  0.000113    1813.48
  500     800     50 |    16.658  0.716 |     72.618    4.36 |   0.000152  109252.3 |  0.000069    1065.88
 1000     400     25 |    10.850  0.490 |     49.826    4.59 |   0.000100  108279.3 |  0.000047     713.04
 2000     200     12 |     7.275  0.321 |     33.821    4.65 |   0.000064  113054.8 |  0.000032     477.00
 5000      80      5 |     4.317  0.178 |     20.077    4.65 |   0.000037  117911.8 |  0.000019     280.36

=== Boundary check: Part (c) line 129 summation-by-parts identity at k = 2 (m0 = n) ===
  k= 2 (m0=200000, n=200000): LHS=0.0000000000e+00, line-129 RHS=1.8778656107e-08, diff=+1.8778656107e-08
  k= 3 (m0=133334, n=200000): LHS=9.5726460929e-04, line-129 RHS=9.5726460929e-04, diff=-7.5894152074e-19
  k=10 (m0= 40000, n=200000): LHS=4.9390758817e-03, line-129 RHS=4.9390758817e-03, diff=+1.7347234760e-18

=== Planted deliberate errors for Test (v) ===
  Planted error A1 (omit a_{m0} in Part (a)): caught at 3 values of k, e.g. k=2: E[Nk]=751.15 > UB_broken=0.00
  Planted error A2 (wrong factor 0.2x in Part (a)): caught at 8/12 values of k, e.g. k=2: E[Nk]=751.15 > UB_broken=150.23
  Note: at k=10, E[Nk]=175.58 vs a_{m'}=63.58 (ratio=2.762); testing planted error B (a_{m'}*10 <= E[Nk]): 175.58 < 635.75 -> CAUGHT!
```

---

## 3. Conclusion

**Theorem D (parts (a)–(d)) and the Corollary follow.** All core lemmas and arguments are completely valid, and the two minor boundary/wording slips (line 129 `=` at `k = 2` where `m0 = n`, and `eps < mu` in part (d)) are trivial to repair without modifying any statement or constant.
