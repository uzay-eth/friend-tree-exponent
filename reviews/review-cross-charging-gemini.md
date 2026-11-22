VERDICT: CORRECT

**Short justification.** The proof in `attempts/strogatz-fable/proofs/xcharge.md` is mathematically sound: Lemma 1's event list, transition rates, and exact cancellations of the $`\kappa`$- and $`\nu`$-multipliers hold on every tree $`T_n`$ ($`n \ge 3`$); the interval enclosures of Section 4 are valid; and Lemmas 2 and 3 correctly bound every box state, the affine-in-$`y = 1/D`$ tail $`L > L_{\max}`$, and the beyond regime $`c > c_{\max}`$ to yield the product (and hence Gamma-ratio) bounds on $`\mathbb{E}[C_n]`$. Only minor notational clarifications are needed in `proofs/xcharge.md` (explicitly noting how conversion targets transition across the $`c = c_1`$ and $`c = c_{\max}`$ boundaries and omitting conversion when $`L = 0`$, plus the implicit non-degeneracy $`\sup f > 0`$ and $`\delta \ge -3`$ in Lemma 3). An independent exact rational verifier built solely from `proofs/xcharge.md` confirms that `upper_4_3_60_12.txt` (at $`\lambda = 29/50`$: 5,853 inequalities) and `lower_4_3_60_12.txt` (at $`\delta = 11/20`$: 5,851 inequalities) satisfy every sign constraint even without clipping and pass every inequality with zero failures, while failing cleanly on tighter targets and single-entry/multiplier perturbations. Direct tree-level tests bypassing row bookkeeping confirm both drift inequalities in exact rational arithmetic across all 434 unlabelled trees ($`3 \le n \le 11`$), 350 random friend trees ($`n \le 400`$), and 124 hand-made extreme trees. Finally, high-precision Gamma evaluation confirms that `0.5022` and `3.974` are safe-direction roundings. The claimed bound follows for all $`n \ge 3`$.

---

## 1. Line-by-Line Audit of `attempts/strogatz-fable/proofs/xcharge.md`

### Section 1 (Model and notation, lines 11–23) — **VALID**
- For $`T_2`$ (one edge), attaching a leaf to $`W`$ yields deterministically the 3-vertex path $`T_3`$, which has $`C_3 = 1`$ core vertex with $`(L, c, S) = (2, 0, \emptyset)`$.
- For every tree $`T_n`$ with $`n \ge 3`$, every leaf is adjacent to a core vertex (degree $`\ge 2`$), and the core is a nonempty connected subtree (if $`C_n = 1`$, the unique core vertex has $`c_u = 0, L_u = n - 1 \ge 2`$; if $`C_n \ge 2`$, every core vertex has $`c_u \ge 1, D_u = L_u + c_u \ge 2`$).
- Transition completeness and rates: choosing $`V`$ uniformly ($`1/n`$) and $`W`$ uniformly among $`V`$'s neighbours ($`1/D_V`$) lands on a core vertex $`w`$ with $`n \cdot \Pr(W = w) = L_w + \sum_{z \in N(w)} 1/D_z`$ (gain at $`w`$: $`L_w \to L_w + 1`$) and on a leaf of core vertex $`w`$ (when $`L_w \ge 1`$) with $`n \cdot \Pr = L_w/D_w`$ (conversion at $`w`$: $`(L_w, c_w) \to (L_w - 1, c_w + 1)`$ plus a newborn $`(1, 1)`$). Summing over all core vertices gives $`\sum_w L_w + \sum_z c_z/D_z + \sum_w L_w/D_w = (\#\text{leaves}) + (\#\text{core}) = n`$, so the event list is complete and exact.
- Both $`D_u = L_u + c_u`$ and $`c_u`$ are nondecreasing under all transitions.

### Section 2 (States, lines 25–48) — **VALID** (one minor wording note: **NEEDS REPAIR**)
- Parameters $`D_{\text{cap}} \ge 3`$, $`c_1 \ge 1`$, $`L_{\max} \ge D_{\text{cap}}`$, $`c_{\max} \ge c_1`$: valid. Since $`D_{\text{cap}} \ge 3`$ and $`c_1 \ge 1`$, the newborn class $`(1, 1)`$ is always non-TOP ($`1 + 1 = 2 < D_{\text{cap}}`$ and $`1 \le c_1`$).
- TOP absorption and minimum degree: a vertex with $`c_w \ge 1`$ is TOP iff $`D_w \ge D_{\text{cap}}`$ or $`c_w \ge c_1 + 1`$. Since $`D_w`$ and $`c_w`$ never decrease, TOP is absorbing, and every TOP vertex has $`D_w \ge \min(D_{\text{cap}}, c_1 + 1) =: d_T \ge 2`$. **VALID**.
- **NEEDS REPAIR (minor wording, line 40):** Line 29 defines $`\mathrm{cls}(w)`$ only for core vertices with $`c_w \ge 1`$, while line 40 states *"so is every vertex with $`L_u \ge L_{\max}`$"* without excepting $`c_u = 0`$. **Repair:** clarify *"every vertex with $`c_u \ge 1`$ and $`L_u \ge L_{\max}`$"* (or extend $`\mathrm{cls}`$ to $`c_u = 0`$ by the same rule, which makes $`L_u \ge L_{\max} \ge D_{\text{cap}}`$ TOP; either way a $`c_u = 0`$ vertex has no core neighbours, so its class never enters any neighbour's state or multiplier sum).

### Section 3 (The drift identity, Lemma 1, lines 49–109) — **VALID** (notational boundary clarifications: **NEEDS REPAIR**)
- **Events and rates ($`\kappa = \nu = 0`$):** **VALID**, except for the minor boundary-notation repair below.
  - A gain at $`w`$ (rate $`L_w + \rho(w)`$) changes $`x_w`$ by $`L_w \to L_w + 1`$ (giving $`Df`$, which is $`0`$ for $`L_w \ge L_{\max}`$ since $`f(L, \cdot) = f(L_{\max}, \cdot)`$ for all $`L > L_{\max}`$) and changes $`x_u`$ only for detailed neighbours $`u \in N(w)`$ when $`\mathrm{cls}(w) = a = (l, k) \ne \mathrm{TOP}`$, shifting one copy of $`a`$ in $`S_u`$ to $`a^+ = \mathrm{cls}(l + 1, k)`$ at rate $`l + 1/D_u + \rho'(w \setminus u)`$.
  - A conversion at $`w`$ ($`L_w \ge 1`$, rate $`L_w / D_w`$) creates a newborn in state $`(1, 1, \{\mathrm{cls}(L_w - 1, c_w + 1)\})`$ (always detailed; for coarse/beyond $`w`$, $`c_w + 1 > c_1`$ so the newborn sees TOP), transitions $`w`$ from $`(L_w, c_w)`$ to $`(L_w - 1, c_w + 1)`$ (adding $`(1, 1)`$ to $`S_w`$ if $`c_w + 1 \le c_1`$), and changes one copy of $`a = (l, k)`$ to $`a^- = \mathrm{cls}(l - 1, k + 1)`$ in $`S_u`$ for detailed neighbours $`u \in N(w)`$ when $`a \ne \mathrm{TOP}`$.
  - Regrouping by the vertex whose $`f`$-value changes reproduces $`\sum_{u \text{ core}} \mathrm{Row}(u)\big|_{\kappa=\nu=0}`$.
  - **NEEDS REPAIR (notational clarity at conversion boundaries, lines 62, 76, 88):**
    1. Line 62 (and line 88) writes $`f(L - 1, c + 1, S + \{(1, 1)\})`$ for any detailed vertex, whereas when $`c = c_1`$ the post-conversion state is coarse $`(L - 1, c_1 + 1)`$ (if $`c_1 < c_{\max}`$) or beyond $`G`$ (if $`c_1 = c_{\max}`$).
    2. Line 76 writes $`f(L - 1, c + 1)`$ for a coarse vertex, whereas when $`c = c_{\max}`$ the post-conversion state is beyond $`G`$.
    3. Lines 62 and 76 should explicitly note "(absent if $`L = 0`$)" so $`f(-1, \dots)`$ is never evaluated (parallel to line 69's note for $`l = 0`$).
    4. **Repair:** Define $`x_u^{\mathrm{conv}}`$ as the Section 2 state of vertex $`u`$ after $`(L, c) \to (L - 1, c + 1)`$ (detailed $`(L - 1, c + 1, S + \{(1, 1)\})`$ if $`c < c_1`$, coarse $`(L - 1, c + 1)`$ if $`c_1 \le c < c_{\max}`$, and $`G`$ if $`c \ge c_{\max}`$), and state that the conversion term is $`0`$ when $`L = 0`$. (Also clarify in lines 60–61 that $`\sum_{a \in S, a \ne \mathrm{TOP}}`$ sums with multiplicity over the multiset $`S`$, and $`S - a`$, $`S - \mathrm{TOP}`$ denote multiset deletion of one copy.)
- **$`\kappa`$-term cancellation:** **VALID**.
  - Subtracted across all core vertices $`u`$: $`\sum_{u: c_u \ge 1} \rho_{\mathrm{top}}(u) \kappa(\mathrm{cls}(u)) = \sum_{(u, w) \text{ adj. core},\, \mathrm{cls}(w) = \mathrm{TOP}} \kappa(\mathrm{cls}(u))/D_w`$ (line 61 for detailed $`u`$ with $`c_u \ge 1`$, zero for $`c_u = 0`$ since $`\rho_{\mathrm{top}} = 0`$; line 75 for coarse/beyond $`u`$ where $`\mathrm{cls}(u) = \mathrm{TOP}`$).
  - Added across all TOP vertices $`w`$: $`\sum_{w: c_w \ge 1,\, \mathrm{cls}(w) = \mathrm{TOP}} \frac{1}{D_w} \sum_{v \in N(w)} \kappa(\mathrm{cls}(v))`$ (line 67 for detailed TOP $`w`$; line 77 for coarse/beyond $`w`$; empty sum $`0`$ if $`c_w = 0`$).
  - These two sums over ordered adjacent core pairs $`(u, w)`$ with $`\mathrm{cls}(w) = \mathrm{TOP}`$ cancel identically, including TOP–TOP edges, detailed-TOP receivers, coarse, and beyond vertices.
- **$`\nu`$-term cancellation:** **VALID**.
  - Consider ordered 2-edge paths $`(u, w, z)`$ of core vertices with $`u, z \in N(w)`$, $`u \ne z`$, and $`\mathrm{cls}(w) = a \in A_2`$ (forcing $`w`$ to be detailed and non-TOP, with $`c_w \ge 2`$).
  - Summing $`\nu(a, \mathrm{cls}(u))/D_z`$ over $`(u, z)`$ for fixed $`w`$ gives $`\sum_{z \in N(w)} \frac{1}{D_z} N(a; S_w - \mathrm{cls}(z))`$. Splitting by $`\mathrm{cls}(z) \ne \mathrm{TOP}`$ and $`\mathrm{cls}(z) = \mathrm{TOP}`$ gives line 60 plus the $`N`$-part of line 61 in $`\mathrm{Row}(w)`$ (using that $`N(a; S_w - \mathrm{TOP})`$ is identical for every TOP neighbour $`z`$ of $`w`$, and indeed equals $`N(a; S_w)`$ since $`\nu(a, \mathrm{TOP}) = 0`$). For $`\mathrm{cls}(w) \notin A_2`$ (including detailed TOP, coarse, beyond, $`c_w \le 1`$), $`N(\mathrm{cls}(w); \cdot) = 0`$ so nothing is added.
  - Regrouping the same triple sum by $`u`$ gives $`\sum_u \sum_{w \in N(u), \mathrm{cls}(w) = a \in A_2} \nu(a, \mathrm{cls}(u)) \rho'(w \setminus u)`$. This is subtracted in line 65 for detailed non-TOP $`u`$; for detailed TOP $`u`$, line 65 subtracts $`\nu(a, \mathrm{TOP}) = 0`$; for coarse/beyond $`u`$, nothing is subtracted, matching $`\nu(a, \mathrm{TOP}) = 0`$; for non-TOP $`w`$ with $`k = 1`$, $`\rho'(w \setminus u) = 0`$ and $`\nu(a, \cdot) = 0`$.
  - Cancellation is exact on every tree $`T_n`$ ($`n \ge 3`$).

### Section 4 (Intervals, lines 111–120) — **VALID**
- Detailed vertex $`u = (L, c, S)`$: every term of $`\mathrm{Row}(u)`$ except $`\rho_{\mathrm{top}}(u)`$ and $`\sum_{w \in N(u), \mathrm{cls}(w)=a} \rho'(w \setminus u)`$ (for $`a = (l, k) \in A_2`$, i.e., $`k \ge 2`$) is determined solely by $`(L, c, S)`$ (note that for $`k = 1`$, $`w`$ has no core neighbour other than $`u`$, so $`\rho'(w \setminus u) = 0`$).
  - Each of the $`n_{\mathrm{top}}`$ TOP neighbours has $`D_w \ge d_T`$, so $`0 \le \rho_{\mathrm{top}}(u) \le n_{\mathrm{top}} / d_T`$.
  - Each of the $`m_a`$ neighbours of class $`a = (l, k)`$ ($`k \ge 2`$) has $`k - 1`$ other core neighbours, each of degree $`\ge 2`$, so $`0 \le \sum_{w: \mathrm{cls}(w)=a} \rho'(w \setminus u) \le m_a (k - 1)/2`$.
- Coarse or beyond vertex $`u = (L, c)`$ ($`c > c_1`$) with $`n_{\mathrm{top}} \in \{0, \dots, c\}`$ TOP neighbours:
  - $`\rho_{\mathrm{top}}(u) \in [0, n_{\mathrm{top}}/d_T]`$.
  - $`\rho(u) - \rho_{\mathrm{top}}(u) \in [0, (c - n_{\mathrm{top}})/2]`$ (each of the $`c - n_{\mathrm{top}}`$ non-TOP core neighbours has degree $`\ge 2`$).
  - Every neighbour $`v \in N(u)`$ has $`c_v \ge 1`$, so $`\sum_{v \in N(u)} \kappa(\mathrm{cls}\, v) \in [c \min_b \kappa(b),\, c \max_b \kappa(b)]`$.

### Section 5 (Certificates, Lemma 2 and Lemma 3, lines 121–164)
- **$`\max(0, \cdot)`$ / $`\min(0, \cdot)`$ treatment:** **VALID**. For any $`t \in [0, h]`$ with $`h \ge 0`$ and real $`e`$, $`\min(0, h e) \le t e \le \max(0, h e)`$.
- **Coarse rows (upper and lower):** **VALID**.
  - Upper: $`(\rho - \rho_{\mathrm{top}}) Df + \rho_{\mathrm{top}}(Df - \kappa(\mathrm{TOP})) \le (c - n_{\mathrm{top}}) \frac{1}{2} Df^+ + n_{\mathrm{top}} \frac{1}{d_T}(Df - \kappa(\mathrm{TOP}))^+ \le \max(0, \frac{c}{2} Df, \frac{c}{d_T}(Df - \kappa(\mathrm{TOP})))`$ since $`(c - n_{\mathrm{top}}) + n_{\mathrm{top}} = c`$ and $`0 \le n_{\mathrm{top}} \le c`$. (In fact, since $`\kappa(\mathrm{TOP}) \ge 0`$ and $`d_T \ge 2`$, the third argument is always $`\le \max(0, \frac{c}{2} Df)`$, so keeping it is completely valid and harmless.) Payments satisfy $`\frac{1}{D} \sum_{v \in N(u)} \kappa(\mathrm{cls}\, v) \le \frac{c}{D} K_p`$ since $`K_p = \max(0, \max_b \kappa(b)) \ge \max_b \kappa(b)`$.
  - Lower: since $`\kappa(\mathrm{TOP}) = 0`$, the uncertain credit term collapses to $`\rho(u) Df`$. Because $`d_T \ge 2`$, $`\rho(u) \in [0, c/2]`$, so $`\rho(u) Df \ge \min(0, \frac{c}{2} Df)`$. Since $`\kappa(b) \ge 0`$ for all $`b`$, the coarse payment term $`\frac{1}{D}\sum_v \kappa(\mathrm{cls}\, v) \ge 0`$ can be safely omitted (for detailed TOP vertices, $`(1/D)\sum_{b' \in S}\kappa(b')`$ is exact and stays in $`\mathrm{main}(x)`$).
- **Tail $`L > L_{\max}`$ (upper and lower):** **VALID**.
  - For every integer $`L \ge L_{\max} + 1`$, $`L - 1 \ge L_{\max}`$, so $`f(L, \cdot) = f(L + 1, \cdot) = f(L_{\max}, \cdot)`$ ($`Df = 0`$) and the post-conversion $`f`$-value $`f(L - 1, c + 1, \dots)`$ equals its value at $`L_{\max}`$ (or $`f(G)`$ if $`c = c_{\max}`$), independent of $`L`$.
  - Since $`L + c \ge L_{\max} + 1 \ge D_{\text{cap}} + 1 > D_{\text{cap}}`$, the vertex itself (if $`c \ge 1`$) is TOP, $`N(\mathrm{TOP}; \cdot) = 0`$, $`\nu(a, \mathrm{TOP}) = 0`$, and after conversion $`(L - 1, c + 1)`$ also has degree $`L + c > D_{\text{cap}}`$ and $`c + 1 \ge 1`$, so the newborn always sees TOP.
  - Setting $`y = 1/D = 1/(L + c) \in (0, 1/(L_{\max} + 1 + c)]`$, the only $`L`$-dependent rates remaining are $`L/D = 1 - c y`$ and $`1/D = y`$ (the $`L \cdot Df`$ term vanishes identically because $`Df = 0`$); every uncertain bound $`\mathrm{hi}_j`$ and coefficient $`e_j`$ is independent of $`y`$, and the RHS ($`\lambda f`$ or $`\delta f`$) is constant in $`y`$.
  - Therefore, the LHS of (U) and (L) is strictly an affine function $`A + B y`$ on $`(0, 1/(L_{\max} + 1 + c)]`$. An affine inequality on $`(0, y_{\max}]`$ holds iff it holds at the endpoints $`y = 0`$ and $`y = y_{\max} = 1/(L_{\max} + 1 + c)`$. (For $`c = 0`$, $`A + B y`$ is constant in $`y`$, so both endpoints agree and still suffice.)
- **Beyond states ($`c > c_{\max}`$):** **VALID**.
  - Upper: $`Df = 0`$, $`f(x_u) = f(x_u^{\mathrm{conv}}) = G`$, so $`\mathrm{Row}(u) = (L/D) f(1, 1, \{\mathrm{TOP}\}) - \rho_{\mathrm{top}}(u) \kappa(\mathrm{TOP}) + (1/D)\sum_{v \in N(u)} \kappa(\mathrm{cls}\, v)`$. Using $`\kappa(\mathrm{TOP}) \ge 0`$, $`\rho_{\mathrm{top}}(u) \ge 0`$, and $`\kappa(\mathrm{cls}\, v) \le K_p`$ gives $`\mathrm{Row}(u) \le (L/D) f(1, 1, \{\mathrm{TOP}\}) + (c/D) K_p \le \max(f(1, 1, \{\mathrm{TOP}\}), K_p) \le \lambda G`$ (since $`L, c \ge 0`$ and $`L/D + c/D = 1`$; this covers $`L = 0`$ as well).
  - Lower: $`f(G) = 0`$, $`Df = 0`$, $`\kappa(\mathrm{TOP}) = 0`$, $`f(1, 1, \{\mathrm{TOP}\}) \ge 0`$, and $`\kappa(b) \ge 0`$ for all $`b`$, so $`\mathrm{Row}(u) = (L/D) f(1, 1, \{\mathrm{TOP}\}) + (1/D)\sum_v \kappa(\mathrm{cls}\, v) \ge 0 = \delta f(G)`$ automatically.
- **Sign constraints on $`\kappa`$:** **VALID**.
  - Upper requires $`\kappa(\mathrm{TOP}) \ge 0`$ (needed to drop $`-\rho_{\mathrm{top}}(u)\kappa(\mathrm{TOP}) \le 0`$ in beyond rows; note also $`K_p = \max(0, \max_b \kappa(b)) = \max_b \kappa(b) \ge 0`$ automatically) and allows non-TOP $`\kappa(b)`$ of either sign.
  - Lower requires $`\kappa(b) \ge 0`$ for all $`b`$ (needed to drop coarse and beyond payment terms $`\ge 0`$) and $`\kappa(\mathrm{TOP}) = 0`$ (needed to drop $`-\rho_{\mathrm{top}}(u)\kappa(\mathrm{TOP})`$ in coarse and beyond rows).
- **Final product formula and Gamma constants:** **VALID** (minor implicit-hypothesis note for Lemma 3: **NEEDS REPAIR**).
  - Upper (Lemma 2): $`f \ge 1, G \ge 1`$ and $`f(1, 1, \{\mathrm{TOP}\}) \le \lambda G`$ force $`\lambda \ge 1/G > 0`$, so $`1 + \lambda/m > 0`$ for all $`m \ge 3`$. By Lemma 1 and (U) + tail + beyond, $`\mathbb{E}[\Phi_{m+1}] \le (1 + \lambda/m)\mathbb{E}[\Phi_m]`$ for all $`m \ge 3`$. Inducting from $`\Phi_3 = f(x_0)`$ and using $`C_n \le \Phi_n`$ (from $`f \ge 1, G \ge 1`$) yields $`\mathbb{E}[C_n] \le f(x_0) \prod_{m=3}^{n-1}(1 + \lambda/m)`$.
  - Lower (Lemma 3): **NEEDS REPAIR (minor implicit hypotheses):** Lemma 3 divides by $`\sup f`$ and multiplies $`\mathbb{E}[\Phi_{m+1}] \ge (1 + \delta/m)\mathbb{E}[\Phi_m]`$ inductively, which requires $`\sup f > 0`$ (i.e., $`f \not\equiv 0`$) and $`\delta \ge -3`$ (so $`1 + \delta/m \ge 0`$). **Repair:** state explicitly "$`f \not\equiv 0`$ (equivalently $`\sup f > 0`$) and $`\delta \ge -3`$ (in particular $`\delta > 0`$)." Both hold trivially for our certificate ($`\sup f \approx 1.7699 > 0`$, $`\delta = 11/20 > 0`$). Then $`\Phi_n \le (\sup f) C_n`$ gives $`\mathbb{E}[C_n] \ge (f(x_0)/\sup f)\prod_{m=3}^{n-1}(1 + \delta/m)`$.
  - Product-to-Gamma identity (exact for all $`n \ge 3`$, $`\alpha > -3`$, using $`\Gamma(z+1) = z\Gamma(z)`$ and $`\Gamma(3) = 2`$):
    ```math
    \prod_{m=3}^{n-1}\left(1 + \frac{\alpha}{m}\right) = \prod_{m=3}^{n-1}\frac{m+\alpha}{m} = \frac{\Gamma(n+\alpha)/\Gamma(3+\alpha)}{\Gamma(n)/\Gamma(3)} = \frac{2}{\Gamma(3+\alpha)}\,\frac{\Gamma(n+\alpha)}{\Gamma(n)}.
    ```

---

## 2. Independent Exact Rational Verification & Negative Controls

Script: `verify_certs.py` (written solely from `proofs/xcharge.md` and the text format in `xdump.py`; does not import `xcharge.py` or `xverify.py`; parses every C99 hex float directly to an exact `fractions.Fraction` and cross-checks against `Fraction.from_float(float.fromhex(...))`).

### Summary of results
- **Completeness & raw sign constraints (no clipping needed):** Both `upper_4_3_60_12.txt` and `lower_4_3_60_12.txt` contain every expected box state (5,116 detailed + 549 coarse = 5,665 entries, unit `starts = 0..60`, all 5 non-TOP classes + TOP, 6 $`\kappa`$ values, 15 $`\nu`$ values, plus $`G`$ and $`K_p`$ in upper) and satisfy **all** sign constraints with **0 violations even without clipping** (for upper, the recomputed $`K_p = \max(0, \max_b \kappa(b))`$ matches the file's `Kp` entry identically).
- **Upper certificate at $`\lambda = 29/50`$:** **PASS** — **5,853** inequalities checked (`5116 det_box + 549 coarse_box + 168 det_tail + 18 coarse_tail + 2 beyond`), **0 failed**, smallest slack $`9.999999999992 \times 10^{-5}`$ at `('det_box', 0, 2, (5, 5))` (beyond slacks: $`8.25845 \times 10^{-4}`$ for $`f(1,1,\{\mathrm{TOP}\}) \le \lambda G`$ and $`2.50462`$ for $`K_p \le \lambda G`$).
- **Lower certificate at $`\delta = 11/20`$:** **PASS** — **5,851** inequalities checked (`5116 det_box + 549 coarse_box + 168 det_tail + 18 coarse_tail`), **0 failed**, smallest slack $`0`$ (exact) at the 11 states with $`L = 0`$ and all-TOP/coarse neighbours where $`f(x) = 0`$ (`(0,2,(5,5))`, `(0,3,(5,5,5))`, and coarse `(0,4)..(0,12)`), and smallest strictly positive slack $`3.742778289501 \times 10^{-5}`$ at `('det_box', 0, 3, (2, 5, 5))`. Every state with $`c \le 1`$ has $`f \ge 1 > 0`$, guaranteeing $`\Phi(T_n) > 0`$ on every tree $`T_n`$ ($`n \ge 3`$).
- **Negative controls — tighter targets:**
  - Upper at $`\lambda = 579/1000`$: **FAIL** (5,668 / 5,853 failed; min slack $`-7.4326 \times 10^{-3}`$).
  - Lower at $`\delta = 552/1000`$: **FAIL** (5,747 / 5,851 failed; min slack $`-2.6548 \times 10^{-3}`$).
- **Negative controls — single-entry and single-multiplier perturbations:** Every table or multiplier perturbation tested on upper (`f(2,0,-) - 1/1000`, `f(0,2,(5,5)) ± 1/1000`, `kappa(1) + 1/1000`, `nu(1,0) + 1/1000`) and lower (`f(2,0,-) + 1/1000`, `kappa(0) - 1/1000`, `nu(1,0) + 1/100`) caused the verifier to fail (with `nu(1,0) + 1/1000` in lower still within slack, turning into a failure at `+1/100`).

### Command and real output (Items 2 and 4)
```text
$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 attempts/strogatz-fable/.venv/bin/python checks/fable-xcharge-gemini/verify_certs.py
=== SIGN CONSTRAINTS (RAW, NO CLIPPING) ===
UPPER: 0 violations (checked 5665 f-entries, G, 6 kappas)
LOWER: 0 violations (checked 5665 f-entries, G, 6 kappas)

=== PRIMARY VERIFICATION (NO CLIPPING) ===
UPPER @ 29/50: ok=True, checked=5853 {'det_box': 5116, 'coarse_box': 549, 'det_tail': 168, 'coarse_tail': 18, 'beyond': 2}, failed=0
  min_slack = 9.999999999992e-05 at ('det_box', 0, 2, (5, 5))
  Kp computed=2.285277145679, Kp_file=2.285277145679, exact_equal=True
  beyond [f(1,1,{TOP}) <= lam*G]: slack = 8.258450023400e-04
  beyond [Kp <= lam*G]: slack = 2.504623867889e+00
LOWER @ 11/20: ok=True, checked=5851 {'det_box': 5116, 'coarse_box': 549, 'det_tail': 168, 'coarse_tail': 18, 'beyond': 0}, failed=0
  min_slack = 0.000000000000e+00 at ('det_box', 0, 2, (5, 5))
  zero_slack count = 11: [('det_box', 0, 2, (5, 5)), ('det_box', 0, 3, (5, 5, 5)), ('coarse_box', 0, 4), ('coarse_box', 0, 5), ('coarse_box', 0, 6), ('coarse_box', 0, 7), ('coarse_box', 0, 8), ('coarse_box', 0, 9), ('coarse_box', 0, 10), ('coarse_box', 0, 11), ('coarse_box', 0, 12)]
  min_pos_slack = 3.742778289501e-05 at ('det_box', 0, 3, (2, 5, 5))

LOWER zero f states: 2 detailed [(0, 2, (5, 5)), (0, 3, (5, 5, 5))], 9 coarse [(0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (0, 10), (0, 11), (0, 12)]
LOWER min f among c <= 1 states: 1.000000000000000 (> 0 guarantees Phi(T) > 0 for every tree T_n, n >= 3)

=== NEGATIVE CONTROLS: TIGHTER TARGETS ===
UPPER @ 579/1000: ok=False, checked=5853, failed=5668, min_slack=-7.432605021085e-03 at ('det_tail', 'y=0', 3, (0, 0, 4))
LOWER @ 552/1000: ok=False, checked=5851, failed=5747, min_slack=-2.654830237718e-03 at ('det_box', 60, 0, ())

=== NEGATIVE CONTROLS: SINGLE-ENTRY / MULTIPLIER PERTURBATIONS ===
  UPPER f(2,0,-) -= 1/1000: ok=False, sign_viols=0, failed=1, min_slack=-2.858018253534e-03 at ('det_box', 2, 0, ())
  UPPER f(0,2,(5,5)) += 1/1000 (sign-safe, tight state): ok=False, sign_viols=0, failed=3, min_slack=-2.373784852613e-03 at ('det_box', 0, 2, (2, 5))
  UPPER f(0,2,(5,5)) -= 1/1000 (violates f >= 1): ok=False, sign_viols=1, failed=1, min_slack=-9.800000000001e-04 at ('det_box', 0, 2, (5, 5))
  UPPER kappa(1) += 1/1000: ok=False, sign_viols=0, failed=36, min_slack=-4.092742764908e-04 at ('coarse_box', 0, 4)
  UPPER nu(1,0) += 1/1000: ok=False, sign_viols=0, failed=6, min_slack=-7.579049039732e-04 at ('det_box', 0, 2, (0, 0))
  LOWER f(2,0,-) += 1/1000: ok=False, sign_viols=0, failed=1, min_slack=-2.769309700379e-03 at ('det_box', 2, 0, ())
  LOWER kappa(0) -= 1/1000: ok=False, sign_viols=0, failed=3, min_slack=-2.206344630419e-04 at ('det_box', 1, 3, (0, 0, 0))
  LOWER nu(1,0) += 1/1000: ok=True, sign_viols=0, failed=0, min_slack=0.000000000000e+00 at ('det_box', 0, 2, (5, 5))
  LOWER nu(1,0) += 1/100: ok=False, sign_viols=0, failed=1, min_slack=-4.417592691267e-03 at ('det_box', 1, 1, (1,))

=== CONSTANTS CHECK (ITEM 4, 50-digit mpmath) ===
UPPER: f(2,0,-) = 7.21981746463390017254369 (exact hex: 0x1.ce117d456c4dap+2)
UPPER: inf f = 1.000000000000000, sup f = 8.258450023393921, G = 8.258450023393804
UPPER: Gamma(3.58) = 3.633782668569140425579337
UPPER: f(x0)*2/Gamma(3.58) = 3.973720017480747047034255
UPPER: <= 3.974 ? True (margin = +0.000279982519253)
UPPER at n=3: 3.974 * Gamma(3.58)/Gamma(3) = 7.22032616244688 >= E[C_3]=1
LOWER: f(2,0,-) = 1.561380599241928068821039 (exact hex: 0x1.8fb6a3925a5bdp+0)
LOWER: inf f = 0.000000000000000, sup f = 1.769886825145300379702462 at ('f', (60, 0, ()))
LOWER: f(x0)/sup f = 0.882192339678976471783771
LOWER: Gamma(3.55) = 3.513252144693200644551331
LOWER: (f(x0)/sup f)*2/Gamma(3.55) = 0.5022083831993305882786317
LOWER: >= 0.5022 ? True (margin = +8.38319933059e-6)
LOWER at n=3: 0.5022 * Gamma(3.55)/Gamma(3) = 0.882177613532463 <= E[C_3]=1
```

---

## 3. Direct Tree-Level Verification (Bypassing Row Bookkeeping)

Script: `verify_trees.py`. Computes $`\Phi(T) = \sum_{u \text{ core}} f(x_u)`$ directly from Section 2's state definition, enumerates all directed moves $`(V, W)`$ with exact rational weight $`n \cdot \Pr(V, W) = 1/\deg(V)`$ (verifying $`\sum_{V, W} n\Pr(V, W) = n`$), and checks $`n \mathbb{E}[\Phi(T') - \Phi(T)] \le (29/50)\Phi(T)`$ (upper) and $`\ge (11/20)\Phi(T)`$ (lower) in exact `Fraction` arithmetic. The local-delta evaluator was cross-checked for exact rational equality against a full from-scratch $`O(n^2)`$ brute-force tree reconstruction on **all 434** unlabelled trees ($`3 \le n \le 11`$), **15** random trees ($`n \le 80`$), and **77** extreme trees ($`n \le 120`$).

### Summary of smallest relative slacks (0 violations across all 908 trees)

Here relative slack is reported both as exponent slack ($`\lambda - \mathrm{drift}/\Phi`$ for upper; $`\mathrm{drift}/\Phi - \delta`$ for lower) and normalized by $`\lambda`$ or $`\delta`$.

| Family | Trees | Upper viols | Upper min $`\lambda - \mathrm{drift}/\Phi`$ (min $`/ \lambda`$) | Lower viols | Lower min $`\mathrm{drift}/\Phi - \delta`$ (min $`/ \delta`$) |
|---|---:|---:|---|---:|---|
| (i) All unlabelled trees $`3 \le n \le 11`$ | 434 | 0 | $`1.000000 \times 10^{-4}`$ ($`1.724 \times 10^{-4}`$) at `n=8#12` | 0 | $`5.000000 \times 10^{-4}`$ ($`9.091 \times 10^{-4}`$) at `n=8#12` |
| (ii) Random friend trees ($`12 \le n \le 400`$) | 350 | 0 | $`1.000000 \times 10^{-4}`$ ($`1.724 \times 10^{-4}`$) at `rand#89(n=48)` | 0 | $`5.000000 \times 10^{-4}`$ ($`9.091 \times 10^{-4}`$) at `rand#73(n=38)` |
| (iii-a) Stars (incl. $`L > 60`$, up to $`n = 400`$) | 12 | 0 | $`1.000000 \times 10^{-4}`$ ($`1.724 \times 10^{-4}`$) at `Star(n=59,L=58)` | 0 | $`5.000000 \times 10^{-4}`$ ($`9.091 \times 10^{-4}`$) at `Star(n=60,L=59)` |
| (iii-b) Paths (up to $`n = 350`$) | 9 | 0 | $`1.000000 \times 10^{-4}`$ ($`1.724 \times 10^{-4}`$) at `Path(n=4)` | 0 | $`5.000000 \times 10^{-4}`$ ($`9.091 \times 10^{-4}`$) at `Path(n=4)` |
| (iii-c) Brooms & Double Brooms ($`L`$ up to 300) | 67 | 0 | $`1.000000 \times 10^{-4}`$ ($`1.724 \times 10^{-4}`$) at `Broom(core=2,L_end=2)` | 0 | $`5.000000 \times 10^{-4}`$ ($`9.091 \times 10^{-4}`$) at `Broom(core=2,L_end=5)` |
| (iii-d) Caterpillars with $`L = 0`$ in middle | 6 | 0 | $`4.722768 \times 10^{-4}`$ ($`8.143 \times 10^{-4}`$) at `Cat_L0_c2_len50_e2` | 0 | $`1.846648 \times 10^{-2}`$ ($`3.358 \times 10^{-2}`$) at `Cat_L0_c3_len8_e80` |
| (iii-e) Vertices with $`c > 12`$ (incl. $`(0,2)`$/$`(0,3)`$ arms, depth-2) | 19 | 0 | $`4.822043 \times 10^{-3}`$ ($`8.314 \times 10^{-3}`$) at `CoarseHub(c=12,L=0,arms=(0,3))` | 0 | $`5.646957 \times 10^{-4}`$ ($`1.027 \times 10^{-3}`$) at `HighC(c=14,L=0,arms=20)` |
| (iii-f) Vertices with $`L > 60`$ (tail, up to $`L = 200`$) | 11 | 0 | $`3.650018 \times 10^{-3}`$ ($`6.293 \times 10^{-3}`$) at `HighL(L=100,c=2,nbs=(1,1),(1,1))` | 0 | $`5.630672 \times 10^{-4}`$ ($`1.024 \times 10^{-3}`$) at `HighL(L=61,c=1,nb=(1,1))` |

Notice how the tightest tree-level slacks match the LP solve margins ($`0.58 - 0.5799 = 10^{-4}`$ for upper and $`0.5505 - 0.55 = 5 \times 10^{-4}`$ for lower) to 9–11 significant figures.

### Command and real output
```text
$ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 attempts/strogatz-fable/.venv/bin/python checks/fable-xcharge-gemini/verify_trees.py
Enumerated all 434 unlabelled trees with 3 <= n <= 11 in 0.06s
Family: (i) All unlabelled trees n<=11 (count=434, bruteforce_checked=True)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 9.999999979811e-05 [= 1.724137927554e-04 * lam] at ('n=8#12', 8, 0.5799000000002019)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 4.999999999968e-04 [= 9.090909090851e-04 * del] at ('n=8#12', 8, 0.5504999999999968)
Cross-checked fast vs brute-force on 15 random trees (n<=80): EXACT MATCH.
Family: (ii) 350 random friend trees (12 <= n <= 400) (count=350, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 9.999999999626e-05 [= 1.724137930970e-04 * lam] at ('rand#89(n=48)', 48, 0.5799000000000037)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 4.999999999957e-04 [= 9.090909090831e-04 * del] at ('rand#73(n=38)', 38, 0.5504999999999957)
Family: (iii-a) Stars (incl. L > 60) (count=12, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 9.999999999670e-05 [= 1.724137930978e-04 * lam] at ('Star(n=59,L=58)', 59, 0.5799000000000033)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 4.999999999946e-04 [= 9.090909090811e-04 * del] at ('Star(n=60,L=59)', 60, 0.5504999999999947)
Family: (iii-b) Paths (count=9, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 9.999999999914e-05 [= 1.724137931020e-04 * lam] at ('Path(n=4)', 4, 0.5799000000000009)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 4.999999999996e-04 [= 9.090909090901e-04 * del] at ('Path(n=4)', 4, 0.5504999999999995)
Family: (iii-c) Brooms & Double Brooms (count=67, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 9.999999999898e-05 [= 1.724137931017e-04 * lam] at ('Broom(core=2,L_end=2)', 5, 0.579900000000001)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 4.999999999995e-04 [= 9.090909090900e-04 * del] at ('Broom(core=2,L_end=5)', 8, 0.5504999999999995)
Family: (iii-d) Caterpillars with L=0 in middle (count=6, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 4.722767595840e-04 [= 8.142702751448e-04 * lam] at ('Cat_L0_c2_len50_e2', 54, 0.579527723240416)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 1.846647995436e-02 [= 3.357541809884e-02 * del] at ('Cat_L0_c3_len8_e80', 180, 0.5684664799543597)
Family: (iii-e) Vertices with c > 12 (beyond G) (count=19, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 4.822043286228e-03 [= 8.313867734876e-03 * lam] at ('CoarseHub(c=12,L=0,arms=(0,3))', 61, 0.575177956713772)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 5.646956995779e-04 [= 1.026719453778e-03 * del] at ('HighC(c=14,L=0,arms=20)', 295, 0.5505646956995779)
Family: (iii-f) Vertices with L > 60 (tail) (count=11, bruteforce_checked=False)
  UPPER (lam=29/50): viols=0, min_rel_slack(lam - drift/Phi) = 3.650017748350e-03 [= 6.293134048879e-03 * lam] at ('HighL(L=100,c=2,nbs=(1,1),(1,1))', 105, 0.5763499822516504)
  LOWER (del=11/20): viols=0, min_rel_slack(drift/Phi - del) = 5.630672119746e-04 [= 1.023758567227e-03 * del] at ('HighL(L=61,c=1,nb=(1,1))', 64, 0.5505630672119746)
Cross-checked fast vs brute-force on 77 extreme trees (n<=120): EXACT MATCH.
Total elapsed: 2.23s
```

---

## 4. Constants Check

Using 50-digit `mpmath` arithmetic on the exact rational table entries (output shown in Section 2):

- **Upper constant ($`\lambda = 29/50 = 0.58`$):**
  - $`f(x_0) = f(2, 0, -) = \texttt{0x1.ce117d456c4dap+2} = 7.21981746463390017254369\ldots`$
  - $`\inf f = 1`$ (satisfying $`f \ge 1, G \ge 1`$); $`\sup f = 8.258450023393921\ldots`$
  - $`\Gamma(3 + 0.58) = 3.633782668569140425579337\ldots`$
  - $`f(x_0) \cdot \dfrac{2}{\Gamma(3.58)} = 3.973720017480747047034255\ldots < 3.974`$ (safe rounding **up**, margin $`+2.7998 \times 10^{-4}`$).
- **Lower constant ($`\delta = 11/20 = 0.55`$):**
  - $`f(x_0) = f(2, 0, -) = \texttt{0x1.8fb6a3925a5bdp+0} = 1.561380599241928068821039\ldots`$
  - $`\inf f = 0`$; $`\sup f = 1.769886825145300379702462\ldots`$, attained at $`f(60, 0, -)`$ (note $`f(G) = 0 \le \sup f`$).
  - $`f(x_0)/\sup f = 0.882192339678976471783771\ldots`$
  - $`\Gamma(3 + 0.55) = 3.513252144693200644551331\ldots`$
  - $`\dfrac{f(x_0)}{\sup f} \cdot \dfrac{2}{\Gamma(3.55)} = 0.5022083831993305882786317\ldots > 0.5022`$ (safe rounding **down**, margin $`+8.3832 \times 10^{-6}`$).
- **Base case $`n = 3`$ sanity check:** $`0.5022 \cdot \Gamma(3.55)/\Gamma(3) \approx 0.882178 \le 1 = \mathbb{E}[C_3] \le 7.220326 \approx 3.974 \cdot \Gamma(3.58)/\Gamma(3)`$.

---

## 5. Conclusion

**The claim follows.** Combining Lemma 1, Section 4, Lemmas 2 and 3, the independently verified exact rational certificates `upper_4_3_60_12.txt` (at $`\lambda = 29/50`$) and `lower_4_3_60_12.txt` (at $`\delta = 11/20`$), and the safe-direction Gamma constant bounds establishes rigorously that for all $`n \ge 3`$,
```math
0.5022\,\frac{\Gamma(n + 0.55)}{\Gamma(n)} \;\le\; \mathbb{E}[C_n] \;\le\; 3.974\,\frac{\Gamma(n + 0.58)}{\Gamma(n)}.
```
