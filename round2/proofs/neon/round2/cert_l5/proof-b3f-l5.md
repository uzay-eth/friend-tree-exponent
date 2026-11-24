# Theorem 2: Exact Pure-GMP Drift Certificates for B3-FullH + Endpoint-Only Coarse Level-5 (`version=2`)

**Authors:** Gemini 4 Neon (Neon lead); GPT-6.2 Sol (implementation). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Claim pin:** `52b6b6522b4d36b015420481bed37311231af769`.

**Review:** GPT-6.2 Sol cold proof and code audit plus exact external-witness replay; no rival-family attack before close. **Dependencies:** the [Level-4 drift theorem](../proof-b3f-motif4.md), the existence theorem for the exponent, and the version-2 GMP checker for historical endpoints. See [Neon Claim 4](../../README.md) for review and release scope. Source paths are relative to the archived lane root (`round2/proofs/neon`).

This document states and proves the mathematical soundness of the `version=2` certificate extension implemented in `round2/cert_l5/` (`cert_b3f_l5_format.h`, `cert_b3f_l5_build.c`, `cert_b3f_l5_verify_gmp.c`, `l5_verify_rows.h`, `tree_identity_l5_check.py`). It extends Theorem 1 of `round2/proof-b3f-motif4.md` by adding a 5-vertex rooted core motif potential $`q_5`$ whose support allows coarse-degree endpoints while keeping every interior vertex strictly exact.

---

## 1. Supported 5-Motifs, Projection, and Cutoff Compatibility

Let $`T`$ be any finite rooted tree with core set $`K(T) = \{u \in V(T) : \deg_T(u) \ge 2\}`$, rooted at a fixed core vertex (so $`h_u = 0`$ at the root and $`h_u = 1`$ for every non-root core vertex, $`D_u = L_u + k_u + h_u`$). Let $`\Phi^{(1)}(T)`$ denote the `version=1` eight-table B3-FullH + Level-4 potential defined in `round2/proof-b3f-motif4.md` with parameters $`(D_0, K_T, L_F, K_F, D_{P4}, K_{P4}, D_{Yt}, K_{Yt}, D_{Yb}, K_{Yb}, D_{S3}, K_{S3})`$.

Fix integer Level-5 parameters:
- exact cutoff $`D_c \ge 4`$,
- coarse endpoint bucket count $`B \in \{1, 2\}`$,
- receiver selection bitmask $`\mathrm{mask} \in \{1, 2, 3\}`$ (bit 0 = L4 interior receivers, bit 1 = L4 endpoint receivers; production uses $`\mathrm{mask} = 3`$).

**Cutoff compatibility (enforced by `b3f_l5_validate_header`):**
```math
D_c \le \min(D_0, D_{P4}, D_{Yt}, D_{Yb}, D_{S3}), \qquad (B = 2 \implies D_c < D_0), \qquad \min(K_{P4}, K_{Yt}, K_{Yb}, K_{S3}) \ge D_c - 1.
```

**Lemma 1 (Exact-state completeness in Level-4).** *Under cutoff compatibility, every core state $`(L, k, h)`$ with $`D = L + k + h < D_c`$ and $`h \in \{0, 1\}`$ satisfies $`D < \min(D_{P4}, D_{Yt}, D_{Yb}, D_{S3})`$ and $`k \le D - h \le D_c - 1 \le \min(K_{P4}, K_{Yt}, K_{Yb}, K_{S3})`$. In particular, every exact state (including apex $`(L, k, h) = (0, D_c - 1, 0)`$) is present in every Level-4 table whose arity condition on $`k`$ it meets.*

**Supported 5-motifs.** Let $`\mathcal{M}_5(T)`$ be the set of all 5-vertex connected induced subtrees $`M`$ of the core tree $`K(T)`$ (each rooted at its unique topmost vertex in $`T`$, i.e., closest to the root of $`T`$). For $`v \in V(M)`$, let $`\deg_M(v)`$ denote its degree *within* $`M`$. Call $`v`$ an **interior vertex** of $`M`$ if $`\deg_M(v) \ge 2`$ and an **endpoint** of $`M`$ if $`\deg_M(v) = 1`$.

Define the supported set
```math
\mathcal{M}_5^{\mathrm{supp}}(T) := \{M \in \mathcal{M}_5(T) : D_v < D_c \text{ for every } v \in V(M) \text{ with } \deg_M(v) \ge 2\}.
```
Thus **every interior vertex of a supported 5-motif is exact ($`D_v < D_c`$)**, while an endpoint $`v`$ ($`\deg_M(v) = 1`$) may have arbitrary degree $`D_v \ge 2`$.

For $`M \in \mathcal{M}_5^{\mathrm{supp}}(T)`$, define the **coarse-endpoint projection** $`\mathrm{proj}(M)`$ by replacing each vertex state $`s_v = (L_v, k_v, h_v)`$ as follows:
- If $`D_v < D_c`$: retain the exact inner state $`s_v`$.
- If $`D_v \ge D_c`$ (which can occur only at an endpoint $`\deg_M(v) = 1`$):
  - for $`B = 1`$: map $`s_v`$ to the single coarse bucket $`\mathrm{BIG}`$;
  - for $`B = 2`$: map $`s_v`$ to $`\mathrm{BIG\_INNER}`$ if $`D_c \le D_v < D_0`$, and to $`\mathrm{TAIL}`$ if $`D_v \ge D_0`$.
  At the apex of $`M`$ ($`v = 0`$, when $`\deg_M(0) = 1`$), coarse buckets collapse both physical values $`h_0 \in \{0, 1\}`$ into one apex bucket state.

The `version=2` potential is
```math
\Phi^{(2)}(T) := \Phi^{(1)}(T) + \sum_{M \in \mathcal{M}_5^{\mathrm{supp}}(T)} q_5(\mathrm{proj}(M)),
```
where $`q_5`$ is stored as `nL5` signed integers scaled by $`S`$ (`I5 = S * q5`) in canonical rooted-shape / apex-$`h`$ / multiset colex rank order. Because $`q_5`$ depends only on exact states with $`D < D_c`$ and $`B \in \{1, 2\}`$ coarse endpoint buckets, **`nL5` is completely independent of $`D_0`$**.

---

## 2. Exact Generator Decomposition and Zero-Sum L4 $`\to`$ L5 Transfer

Let $`\mathcal{L}\Phi^{(2)}(T) = n\,\mathbb{E}[\Phi^{(2)}(T_{n+1}) - \Phi^{(2)}(T_n) \mid T_n = T]`$ be the continuous-time generator action. By linearity,
```math
\mathcal{L}\Phi^{(2)}(T) = \mathcal{L}\Phi^{(1)}(T) + \sum_{M \in \mathcal{M}_5^{\mathrm{supp}}(T)} \mathcal{L}_{\mathrm{own}} q_5(M) + \sum_{M' \in \mathcal{M}_5^{\mathrm{supp}}(T_{\mathrm{after}}) \setminus \mathcal{M}_5^{\mathrm{supp}}(T)} \mathrm{rate}(M')\, q_5(\mathrm{proj}(M')).
```

### 2.1 Newborn 5-Motif Ownership (Zero Large-Promoter Newborns)

A new 5-motif $`M' \in \mathcal{M}_5^{\mathrm{supp}}`$ is created at a single step of the friend tree **if and only if** a leaf $`\ell`$ attached to an existing core vertex $`u`$ is promoted (turning $`\ell`$ into a new core vertex of state $`v_{\mathrm{nb}} = (1, 0, 1)`$ with degree $`2 < D_c`$) and $`M' = m \cup \{\ell\}`$ for a 4-vertex connected core subtree $`m \subset K(T)`$ containing $`u`$.

Crucially, in $`M' = m \cup \{\ell\}`$, vertex $`u`$ is adjacent both to $`\ell`$ and to at least one vertex of $`m`$, so **$`\deg_{M'}(u) \ge 2`$ ($`u`$ is an interior vertex of $`M'`$)**. Therefore:
1. If $`D_u \ge D_c`$ (including all tail states $`D_u \ge D_0`$ and all overflow states $`k_u > K_T`$), then $`M' \notin \mathcal{M}_5^{\mathrm{supp}}`$ — **no large promoter can ever create a supported 5-motif!**
2. If $`D_u < D_c`$, then every interior vertex of $`m`$ is also an interior vertex of $`M'`$, so $`M' \in \mathcal{M}_5^{\mathrm{supp}}`$ if and only if every interior vertex of $`m`$ has degree $`< D_c`$ (i.e., $`m`$ is a supported 4-motif under `l5_project`). Because $`u`$ has $`D_u < D_c`$ and $`L_u \ge 1`$ before promotion, Lemma 1 guarantees that $`m`$ is represented in one of the four active Level-4 tables (`P4, Ytop, Ybot, S3`), with every spectator $`w \in V(m) \setminus \{u\}`$ either exact ($`D_w < D_c`$) or embedded as a convex combination of controls in its coarse bucket fiber (`BIG_INNER` for $`D_c \le D_w < D_0`$, `TAIL` for $`D_w \ge D_0`$).
3. Because $`q_5(\mathrm{proj}(m_{\mathrm{bud}(u)} \cup \{v_{\mathrm{nb}}\}))`$ depends on each coarse spectator $`w`$ **only through its bucket**, it is identical across every control in $`w`$'s Bernstein-3 embedding (including `INF` and overflow $`k_w > K_T`$). Adding $`q_5(\mathrm{proj}(m_{\mathrm{bud}(u)} \cup \{v_{\mathrm{nb}}\}))`$ directly to the promotion delta at promoter $`u`$ (rate $`L_u / D_u`$) inside each active L4 row equation therefore accounts for **every** newborn supported 5-motif with exact weight $`1`$ (since convex weights sum to $`1`$).

### 2.2 Exact Zero-Sum L4 $`\to`$ L5 Edge-Rate Transfer

Consider any active 4-motif $`m`$ in `P4, Ytop, Ybot, S3` that is L5-supported (all interior vertices of $`m`$ have $`D < D_c`$) and any vertex $`u \in V(m)`$ that is exact ($`D_u < D_c`$) and selected by `mask` (`l5_selected`).

In `version=1`, each core neighbour $`z \in N_K(u) \setminus V(m)`$ outside $`m`$ contributes leaf-arrival rate $`1/D_z \in [0, 1/2]`$ to $`u`$, multiplying the L4 jump $`\Delta \Phi_4(m, u) := \Phi_4(m_{\mathrm{up}(u)}) - \Phi_4(m)`$. In `version=2`, split this rate identically as:
```math
\frac{1}{D_z} = \min\!\left(\frac{1}{D_z}, \frac{1}{D_c}\right) + \mathbf{1}_{\{D_z < D_c\}} \left(\frac{1}{D_z} - \frac{1}{D_c}\right).
```
- The first term lies in $`[0, 1/D_c]`$ for **every** external core neighbour $`z`$ (both $`D_z < D_c`$ and $`D_z \ge D_c`$). Replacing $`u`$'s external rate interval $`[0, c_{\mathrm{ext}}/2]`$ by $`[0, c_{\mathrm{ext}}/D_c]`$ in the L4 row equation covers the entire first term.
- For each external core neighbour $`z`$ with $`D_z < D_c`$, attaching $`z`$ as a leaf of $`u`$ forms a unique 5-motif $`M = m \cup \{z\} \in \mathcal{M}_5^{\mathrm{supp}}(T)`$ in which $`z`$ is an exact endpoint ($`\deg_M(z) = 1`$, $`D_z < D_c`$) adjacent to $`u`$. Moving the exact positive residual term
  $`\displaystyle \left(\frac{1}{D_z} - \frac{1}{D_c}\right) \Delta \Phi_4(m, u)`$
  from the L4 equation for $`m`$ into the L5 equation for $`M`$ is an **exact zero-sum identity on every tree $`T`$**.

Because a coarse endpoint spectator $`w \in V(m) \setminus \{u\}`$ has physical L4 contribution given by a convex combination $`\sum_{\gamma \in \mathcal{F}(w)} \alpha_\gamma \Delta I_4(m[w \mapsto \gamma], u)`$ with the **same** nonnegative weights $`\alpha_\gamma`$ ($`\sum \alpha_\gamma = 1`$) before and after the arrival at $`u`$ (since $`w`$'s state does not change), for upper bounds (`mode=1`) and lower bounds (`mode=2`) respectively we have the exact inequalities
```math
\min_{\vec\gamma \in \prod \mathcal{F}} \Delta I_4(m[\vec\gamma], u) \;\le\; S\,\Delta \Phi_4(m, u) \;\le\; \max_{\vec\gamma \in \prod \mathcal{F}} \Delta I_4(m[\vec\gamma], u),
```
where $`\mathcal{F}`$ is the exact fiber of L4 controls compatible with each coarse spectator's role (apex: both $`h \in \{0, 1\}`$; child: $`h = 1`$) and bucket (`BIG_INNER`: inner states $`D_c \le D < D_0`$; `TAIL`: `INF` and all $`T_1, T_2, T_3`$ including overflow $`k = K_T + 1`$; unsplit `BIG`: union). Multiplying this direction-selected exact GMP extremum by $`(1/D_z - 1/D_c) > 0`$ gives the exact L5 incoming transfer term verified in `l5_verify_rows.h`.

### 2.3 Exact Level-5 Own-Drift Bounds

For each vertex $`i \in V(M)`$ of a supported 5-motif $`M \in \mathcal{M}_5^{\mathrm{supp}}(T)`$:
1. **Exact vertex ($`D_i < D_c`$):**
   - *Leaf arrival at $`i`$:* transitions $`s_i`$ to $`(L_i + 1, k_i, h_i)`$ if $`D_i + 1 < D_c`$; to the coarse bucket $`\mathrm{BIG}`$ / $`\mathrm{BIG\_INNER}`$ if $`D_i + 1 = D_c`$ and $`\deg_M(i) = 1`$; and out of support (value $`0`$) if $`D_i + 1 = D_c`$ and $`\deg_M(i) \ge 2`$. Let $`\Delta_i = q_5(\mathrm{after}) - q_5(\mathrm{before})`$. The arrival rate at $`i`$ is
     $`\displaystyle r_i = L_i + \sum_{\substack{j \in N_M(i) \\ D_j < D_c}} \frac{1}{D_j} + \sum_{\substack{j \in N_M(i) \\ D_j \ge D_c}} \frac{1}{D_j} + \sum_{w \in N_K(i) \setminus V(M)} \frac{1}{D_w}.`$
     The first two sums are exact; each coarse motif neighbour $`j`$ has $`1/D_j \in [0, 1/D_c]`$ (`BIG`/`BIG_INNER`) or $`[0, 1/D_0]`$ (`TAIL`); each of the $`c_{\mathrm{ext},i} = k_i + h_i - \deg_M(i)`$ external core neighbours has $`1/D_w \in [0, 1/2]`$. Selecting the adverse endpoint of these intervals according to the sign of $`\Delta_i`$ is exact.
   - *Leaf promotion at $`i`$ ($`L_i \ge 1`$):* occurs at exact rate $`L_i / D_i`$, leaves $`D_i < D_c`$ unchanged, and transitions $`s_i`$ to $`(L_i - 1, k_i + 1, h_i)`$.
2. **Coarse `TAIL` endpoint ($`D_i \ge D_0`$, $`\deg_M(i) = 1`$):** Degree only increases on leaf arrival and is invariant on promotion; hence `TAIL` stays `TAIL` with zero change.
3. **Coarse `BIG_INNER` endpoint ($`D_c \le D_i \le D_0 - 1`$, $`\deg_M(i) = 1`$, $`B = 2`$):**
   - Promotion leaves $`D_i`$ unchanged (stays `BIG_INNER`).
   - Leaf arrival leaves $`i`$ in `BIG_INNER` if $`D_i \le D_0 - 2`$, and crosses to `TAIL` with jump $`\Delta_i^{\mathrm{cross}} = q_5(\mathrm{TAIL}) - q_5(\mathrm{BIG\_INNER})`$ if $`D_i = D_0 - 1`$.
   - At $`D_i = D_0 - 1`$, every core neighbour of $`i`$ has degree $`\ge 2`$, so the leaf-arrival rate at $`i`$ satisfies
     $`\displaystyle 0 \;\le\; r_i = L_i + \sum_{w \sim_K i} \frac{1}{D_w} \;\le\; L_i + \frac{k_i + h_i}{2} \;\le\; L_i + k_i + h_i = D_0 - 1.`$
   - Therefore, uniformly over all $`D_i \in [D_c, D_0 - 1]`$, the exact own-change contribution of a `BIG_INNER` endpoint lies between $`\min(0, (D_0 - 1)\Delta_i^{\mathrm{cross}})`$ and $`\max(0, (D_0 - 1)\Delta_i^{\mathrm{cross}})`$.

Each L5 row inequality is verified in pure GMP integer arithmetic over the common positive denominator $`Q \cdot \mathrm{lcm}(1, \dots, D_0)`$.

---

## 3. Exact Interior-Anchored Comparison Bounds (`c0 * C(T) <= Phi^(2)(T) <= C0 * C(T)`)

For every $`M \in \mathcal{M}_5^{\mathrm{supp}}(T)`$, the set of interior vertices of $`M`$ is nonempty and connected. Define its unique **charge vertex** $`u^*(M) \in V(M)`$ as:
- the apex $`0`$ of $`M`$ if $`\deg_M(0) \ge 2`$;
- the unique child $`1`$ of apex $`0`$ in $`M`$ if $`\deg_M(0) = 1`$ (in which case vertex $`1`$ necessarily has $`\deg_M(1) \ge 2`$ and $`h_{u^*} = 1`$).

In both cases, $`u^*(M)`$ is an **exact interior vertex** of $`M`$, so $`D_{u^*} < D_c`$, and every other interior vertex $`v \in V(M)`$ satisfies $`D_v < D_c`$ and hence has at most $`D_c - 2`$ core children below it.

**Proposition 1 (Per-state, per-shape exact count bound).** *Let $`u \in K(T)`$ have exact state $`g = (L, k, h)`$ with $`D = L + k + h < D_c`$, and let $`S`$ be one of the 9 rooted 5-vertex tree shapes. Define recursively for any rooted tree shape $`R`$ and root child budget $`m \ge 0`$:*
```math
N(R, m) := \begin{cases} 1, & |R| = 1, \\ 0, & \text{root degree of } R > m, \\ \displaystyle\left(\prod_{j=1}^r \binom{m - \sum_{t<j} \alpha_t}{\alpha_j}\right) \prod_{c \in \mathrm{ch}(R)} N(c, D_c - 2), & \text{otherwise}, \end{cases}
```
*where the root's child subtrees partition into $`r`$ distinct isomorphism classes of multiplicities $`\alpha_1, \dots, \alpha_r`$. Then the number of motifs $`M \in \mathcal{M}_5^{\mathrm{supp}}(T)`$ of shape $`S`$ charged to $`u`$ is at most:*
- $`N(S, k)`$ *if the root of $`S`$ has $`\ge 2`$ children;*
- $`h \cdot N(S_{\mathrm{child}}, k)`$ *if the root of $`S`$ has $`1`$ child (with child shape $`S_{\mathrm{child}}`$), since the apex of $`M`$ must be $`u`$'s unique parent (requiring $`h = 1`$).*

*Furthermore, summing over all 9 shapes $`S`$ never exceeds the Lagrange-inversion outer cap*
```math
B(D) := [x^4](1 + F(x))^D = \frac{D}{4}\binom{D - 1 + 4(D_c - 2)}{3}, \qquad F(x) = x(1 + F(x))^{D_c - 2}.
```

*Proof.* Every motif $`M`$ of shape $`S`$ charged to $`u`$ maps the root of $`S`$ (or the unique child of the root when $`S`$ is unary, with the root fixed to $`u`$'s parent) to $`u`$. Distinct unordered child-subtree classes of multiplicity $`\alpha_j`$ choose $`\alpha_j`$ distinct core children from the remaining budget of $`u`$'s $`k`$ core children; every non-leaf child of $`u`$ in $`M`$ is an interior vertex of $`M`$, hence has $`h = 1`$ and $`D < D_c`$, so it has at most $`D_c - 2`$ core children available to embed its subtree. By induction on $`|R|`$, the count of embeddings is at most $`N(R, m)`$. For the outer cap, every supported motif charged to $`u`$ embeds into the full rooted ambient tree with $`D`$ branches at $`u`$ and $`D_c - 2`$ forward branches at every other vertex, whose number of 5-vertex connected subsets containing the root is $`[x^4](1 + F)^D = B(D)`$ by Lagrange inversion. $`\square`$

For every exact state $`g`$ ($`D < D_c`$) and every 5-vertex shape $`S`$, the verifier scans all rows of $`q_5`$ charged to $`(g, S)`$ to compute
```math
m_{g,S} := \min\!\bigl(0, \min_{\text{rows }(g,S)} I_5\bigr) \le 0, \qquad M_{g,S} := \max\!\bigl(0, \max_{\text{rows }(g,S)} I_5\bigr) \ge 0.
```
Adding $`\sum_S b_{g,S} m_{g,S}`$ to the lower per-vertex charge $`\mathrm{loz}_g`$ and $`\sum_S b_{g,S} M_{g,S}`$ to the upper per-vertex charge $`\mathrm{hiz}_g`$ in exact `mpz_t` arithmetic **before** applying the zero-sum edge recharge $`\sum_{u \in K(T)} \frac{1}{4}(k_u - h_u) = 0`$ yields exact integers $`4 S c_0`$ and $`S C_0`$ such that
```math
c_0\, C(T) \;\le\; \Phi^{(2)}(T) \;\le\; C_0\, C(T) \qquad \text{for every finite tree } T \text{ with } C(T) \ge 1.
```

---

## 4. Main Certificate Theorem (Theorem 2)

**Theorem 2.** *If a `version=2` certificate file passes `cert_b3f_l5_verify_gmp` with `total_fails = 0` (covering all 8 Level-1–4 tables, all `nL5` Level-5 rows, header/size validation, and the mode-required strict comparison positivity $`S C_0 > 0, S\Phi^{(2)}(T_3) > 0`$ for both modes plus $`4 S c_0 > 0`$ for upper mode; both pinned production payloads satisfy all three strictly):*
- *In LOWER mode (`mode=2`),* $`\mathcal{L}\Phi^{(2)}(T) \ge \frac{P}{Q}\Phi^{(2)}(T)`$ for every finite tree $`T`$ with $`C(T) \ge 1`$, hence for all $`n \ge 3`$.
  $`\displaystyle \mathbb{E}[C_n] \;\ge\; \frac{\Phi^{(2)}(T_3)}{C_0} \prod_{j=3}^{n-1} \left(1 + \frac{P/Q}{j}\right) \qquad \text{and} \qquad \mu \;\ge\; \frac{P}{Q}.`$
- *In UPPER mode (`mode=1`),* $`\mathcal{L}\Phi^{(2)}(T) \le \frac{P}{Q}\Phi^{(2)}(T)`$ for every finite tree $`T`$ with $`C(T) \ge 1`$, hence for all $`n \ge 3`$.
  $`\displaystyle \mathbb{E}[C_n] \;\le\; \frac{\Phi^{(2)}(T_3)}{c_0} \prod_{j=3}^{n-1} \left(1 + \frac{P/Q}{j}\right) \qquad \text{and} \qquad \mu \;\le\; \frac{P}{Q}.`$

*Proof.* By Sections 2.1–2.3, every event on every tree $`T`$ with $`C(T) \ge 1`$ decomposes identically into verified Level-1–4 row terms, zero-sum L4 $`\to`$ L5 edge transfers, exact-promoter L5 newborns, and verified Level-5 own-drift terms; thus `total_fails = 0` implies the global pointwise drift inequality $`\mathcal{L}\Phi^{(2)}(T) \ge \frac{P}{Q}\Phi^{(2)}(T)`$ (lower) or $`\le \frac{P}{Q}\Phi^{(2)}(T)`$ (upper). Since $`\mathbb{E}[\Phi^{(2)}(T_{n+1}) \mid T_n] = \Phi^{(2)}(T_n) + \frac{1}{n}\mathcal{L}\Phi^{(2)}(T_n)`$, induction from $`T_3`$ combined with Section 3's uniform comparison $`c_0 C(T) \le \Phi^{(2)}(T) \le C_0 C(T)`$ yields both finite-size product bounds and the asymptotic exponent inequalities via $`\mu = \lim \frac{\log \mathbb{E}[C_n]}{\log n}`$. $`\square`$

---

## 5. Verified Benchmarks (Self-Contained Pure-GMP `version=2`)

1. **`D0=10, Dc=6, B=2` (`round2/cert_l5/REPORT.md`):**
   - Certified $`\mu \in [5656316/10000000,\; 5658254/10000000] = [0.5656316,\; 0.5658254]`$ (width `0.0001938`, `-24.6%` vs pure-GMP `D0=10` without L5), `2,127,088` L5 rows, `total_fails = 0` under both standard and ASan+UBSan pure-GMP verifiers, `34/34` GMP trees, `38/38` Python `Fraction` trees, and all structural/override/planted negative controls rejected cleanly.
2. **`D0=12, Dc=7, B=2` (`round2/cert_l5/REPORT.md`):**
   - Certified $`\mu \in [5656682/10000000,\; 5657836/10000000] = [0.5656682,\; 0.5657836]`$ (width `0.0001154`, `-13.6%` vs historical no-L5 float width), `14,215,993` L5 rows, `total_fails = 0`, `34/34` GMP trees, `38/38` Python `Fraction` trees, `4.85M / 7.24M` inward override rejections.
3. **Production `D0=20, Dc=7, B=2` (`B3-FullH-20 + Coarse-L5` — Corollary 4 below):** Pure-GMP certified `[5656929/10000000, 5657570/10000000] = [0.5656929, 0.5657570]` (width `0.0000641`, `-27.2%` narrower than `B3-FullH-20-L4` at the exact same `D0=20` for `+1.6%` payload size).

---

## 6. Corollary 4: Production `B3-FullH-20 + Coarse-L5 (Dc=7, B=2)` Enclosure & Product Bounds

With parameters `(D0, KT, LF, KF) = (20, 19, 16000, 200)`, `P4=(12,9)`, `Yt=(11,8)`, `Yb=(12,9)`, `S3=(10,8)`, `(Dc, buckets, mask) = (7, 2, 3)`, `S = 10^{11}`, `S_r = 10^9`, `Q = 10^7`:

- **Exact binary layout (`11,370,156,264` bytes per file):**
  - Header: `96 + 24 = 120` bytes (`version=2`, `NG=523, INNER_N=396, NT=42, npars=481, nchs=253, NHU=228, NCU=432`).
  - Auxiliary rate/transfer entries: `2*396 + 4*42 + 10,626 + 10,101 = 21,687` `int64_t` values.
  - Stored potential entries across the 9 tables (`1,421,247,831` total `int64_t` values; `1,421,095,995` reported verifier state rows, i.e., `1,421,095,991` active drift and equality checks excluding the 4 `D < 2` padding cells in `f`, with `g[INF, INF] == 0` verified by exact equality):
    - `f`: `6,432,404`; `g`: `273,529` (`121,693` active `(par, ch)` pairs checked); `hV`: `27,746,004`; `hC`: `13,880,592`;
    - `qP4`: `354,856,788`; `qYt`: `540,925,385`; `qYb`: `331,822,656`; `qS3`: `131,094,480`; `q5` (L5): `14,215,993`.
  - Total `int64_t` entries: `1,421,247,831 + 21,687 = 1,421,269,518` (`120 + 8 * 1,421,269,518 = 11,370,156,264` bytes).

- **Pure-GMP verified certificates (`round2/outputs/cert_b3f_20_l5_{lo,hi}.log`):**
  - **Lower certificate** (`SHA-256 9936d424fb3d87bc9193c79093870f1dd23772c16dd8b27a66233bc42a5a8474`):
    - Exact exponent `P/Q = 5656929/10000000 = 0.5656929`, `total_fails = 0` across all 9 tables (`1,421,095,995` reported rows / `1,421,095,991` valid/equality checks), `34/34` GMP concrete trees, `38/38` independent Python `Fraction` trees.
    - Exact comparison integers (matched identically by GMP and Python): `4*S*c0 = 76840591160`, `S*C0 = 184481683938`, `S*Phi^(2)(T_3) = 146042051850`.
  - **Upper certificate** (`SHA-256 cd4553129c758b2c8f74e3f1686146a613e8b4f01f4f720a6a5b86da09b8d317`):
    - Exact exponent `P/Q = 5657570/10000000 = 565757/1000000 = 0.5657570`, `total_fails = 0` across all 9 tables (`1,421,095,995` reported rows / `1,421,095,991` valid/equality checks), `34/34` GMP concrete trees, `38/38` independent Python `Fraction` trees.
    - Exact GMP comparison integers: `4*S*c0 = 88296060680` (matched identically by Python), `S*C0 = 200835990548` (Python joint tail-child extremum gives `200701592409 <= 200835990548`), `S*Phi^(2)(T_3) = 146044250013`.

- **Corollary 4.** *For the undirected random friend tree started from one edge, with zero external prior bracket on $`\mu`$:*
  $`\displaystyle \frac{5656929}{10000000} \;\le\; \mu \;\le\; \frac{565757}{1000000}, \qquad \text{i.e., } \mu \in [0.5656929,\; 0.5657570] \quad (\text{width } 0.0000641).`$
  *Moreover, for every integer $`n \ge 3`$,*
  $`\displaystyle \frac{8113447325}{10248982441} \prod_{j=3}^{n-1} \left(1 + \frac{5656929/10000000}{j}\right) \;\le\; \mathbb{E}[C_n] \;\le\; \frac{146044250013}{22074015170} \prod_{j=3}^{n-1} \left(1 + \frac{565757/1000000}{j}\right),`$
  *where $`\frac{8113447325}{10248982441} = \frac{146042051850}{184481683938}`$ and $`\frac{146044250013}{22074015170} = \frac{4 \cdot 146044250013}{88296060680}`$.*
