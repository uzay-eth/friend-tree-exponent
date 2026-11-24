# Rigorous Unoriented Five-Vertex (Diameter-4) Exponent Certificate (`FT2_U5V1`)

**Pinned claim commit:** `b04cc50c1e391e0420eb7b006ac781d484c68543` (gold1), `b435b4c824e2ba5835f58b85e67a27de7d6a0d49` (gold2).  
**Review and dependencies:** GPT and Claude proof audits and separate historical replays of gold1 and gold2 using the same authored numerical verifier. See [Claim 6 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-6.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-6.md).

**Author:** Gemini 4 Neon (Neon-2).

This note proves that the standalone zero-floating-point verifier (`verify_quintet.cpp`) together with the exact rational tensor-Bernstein tail verifier (`verify_tail_bernstein.py`) certifies rigorous two-sided enclosures $`[\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$ for the core-growth exponent $`\mu`$ of the random friend tree — **without orientation, without any founder boundary lemma, and without occupation-measure post-processing**.

---

## 1. Unoriented Core Tree and Exact Physical Generator

For any finite tree $`T`$ with $`|V(T)| \ge 3`$, let $`K(T) := \{u \in V(T) : \deg_T(u) \ge 2\}`$ be its core set, with $`C(T) := |K(T)| \ge 1`$. Because $`T`$ is a tree, the induced subgraph $`T_K := T[K(T)]`$ is a connected tree on $`C(T)`$ vertices with $`C(T) - 1`$ edges. Each core vertex $`u \in K(T)`$ carries its **unoriented** physical state $`(D_u, c_u)`$ where:
- $`D_u := \deg_T(u) \ge 2`$ is its total degree in $`T`$,
- $`c_u := \deg_{T_K}(u) \in \{0, 1, \dots, D_u\}`$ is its core degree ($`c_u = 0 \iff C(T) = 1`$, and $`c_u \ge 1`$ whenever $`C(T) \ge 2`$),
- $`L_u := D_u - c_u \ge 0`$ is its leaf count, and $`\iota_u := 1/D_u`$.

In the continuous-time friend-tree embedding (each vertex independently chooses a uniformly random neighbor at rate $`1`$ and attaches a fresh leaf to that chosen neighbor), the **exact** core-tree generator on any function $`\Phi`$ of the unoriented core tree $`(T_K, (D_u)_{u \in K(T)})`$ is:
```math
(\mathcal L \Phi)(T) \;=\; \sum_{u \in K(T)} \Big(r_u \, \Delta_u^L \Phi(T) \;+\; b_u \, \Delta_u^B \Phi(T)\Big),
```
where:
- $`r_u := L_u + \sum_{v \sim_{T_K} u} \iota_v`$ is the leaf-gain rate at $`u`$ (each of $`u`$'s $`L_u`$ leaves chooses $`u`$ at rate $`1`$, and each core neighbor $`v \sim_{T_K} u`$ chooses $`u`$ at rate $`1/D_v = \iota_v`$), changing $`(D_u, c_u) \mapsto (D_u + 1, c_u)`$ and leaving $`T_K`$ unchanged;
- $`b_u := L_u / D_u`$ is the leaf-promotion rate at $`u`$ ($`u`$ chooses one of its $`L_u`$ leaves at rate $`L_u / D_u`$, attaching a fresh leaf to that leaf so it becomes a new core leaf $`J`$ with state $`(2, 1)`$ adjacent to $`u`$, while $`u`$ changes from $`(D_u, c_u)`$ to $`(D_u, c_u + 1)`$; zero if $`L_u = 0`$).

Notice that this identity holds **globally at every core vertex of every finite tree $`T`$ with $`C(T) \ge 1`$** (including $`C(T) = 1`$, where $`r_u = D_u`$ and $`b_u = 1`$, and $`C(T) = 2`$), with no root/founder exception and no dummy vertex.

---

## 2. Eight Unoriented Motif Shapes and $`c`$-Capped Endpoint Interpolation

### 2.1. Endpoint spaces with $`c`$-capped inner overflow
At each endpoint level $`j \in \{2, 3, 4, 5\}`$ with integer cutoffs $`(d_j, K_j)`$ ($`2 \le K_j < d_j`$), define the finite corner set $`G_j`$ consisting of four disjoint kinds:
1. **Exact inner** (`kind = 0`): pairs $`(D, c)`$ with $`2 \le D < d_j`$ and $`1 \le c \le \min(D, K_j)`$, reciprocal $`\mathrm{inv} = 1/D`$.
2. **Exact-degree inner overflow** (`kind = 3`): one state $`\mathrm{OVF}_j(D)`$ for each $`K_j + 1 \le D < d_j`$ representing all states $`(D, c)`$ with $`c \ge K_j + 1`$, exact reciprocal $`\mathrm{inv} = 1/D`$.
3. **One-knot degree tail** (`kind = 1`): classes $`T_j(\kappa)`$ for $`\kappa = 1, \dots, K_j + 1`$, representing $`D \ge d_j`$ with $`\min(c, K_j + 1) = \kappa`$, reciprocal $`\mathrm{inv} = 1/d_j`$.
4. **Infinity corner** (`kind = 2`): $`\mathrm{INF}_j`$, reciprocal $`\mathrm{inv} = 0`$.

For any physical core vertex $`u`$ with $`c_u \ge 1`$, its deterministic level-$`j`$ interpolation distribution $`\mu_u^{(j)}`$ on $`G_j`$ is:
```math
\mu_u^{(j)} := \begin{cases}
\delta_{(D_u, c_u)} & \text{if } D_u < d_j \text{ and } c_u \le K_j, \\ \delta_{\mathrm{OVF}_j(D_u)} & \text{if } D_u < d_j \text{ and } c_u > K_j, \\ \frac{d_j}{D_u}\,\delta_{T_j(\min(c_u, K_j+1))} + \big(1 - \frac{d_j}{D_u}\big)\,\delta_{\mathrm{INF}_j} & \text{if } D_u \ge d_j.
\end{cases}
```
Observe that $`\mathbb E_{\mu_u^{(j)}}[\mathrm{inv}] = 1/D_u = \iota_u`$ **identically** for every $`(D_u, c_u)`$.

Whenever $`d_{j+1} \le d_j`$ and $`K_{j+1} \le K_j`$, define the exact stochastic projection $`\Pi_{j \to j+1} : G_j \to \mathcal P(G_{j+1})`$ by applying the exact same rule to the corner label $`(D, \kappa)`$ of $`G_j`$ (mapping $`\mathrm{INF}_j \mapsto \delta_{\mathrm{INF}_{j+1}}`$). Because degree scaling composes ($`\frac{d_j}{D}\cdot\frac{d_{j+1}}{d_j} = \frac{d_{j+1}}{D}`$) and clamping composes ($`\min(\min(c, K_j+1), K_{j+1}+1) = \min(c, K_{j+1}+1)`$), we have the exact identity $`\mu_u^{(j)} \Pi_{j \to j+1} = \mu_u^{(j+1)}`$ for every $`(D_u, c_u)`$.

### 2.2. Eight unoriented motifs and nested downward-closed supports
A quantized potential $`\Phi(T)`$ in integer ticks $`S = 2^{40}`$ is the sum over all unoriented embeddings into $`T_K`$ of **only eight shapes**, covering every connected tree up to 5 vertices (diameter 4):
1. **Singleton** $`f(D_u, c_u)`$: exact table for $`2 \le D_u < D_F`$, analytic polynomial in $`(q, x) = (L_u/D_u, 1/D_u)`$ for $`D_u \ge D_F`$.
2. **Edge** $`g(\{u, v\})`$: symmetric bilinear expectation over $`\mu_u^{(2)} \otimes \mu_v^{(2)}`$, with $`g(\mathrm{INF}_2, \mathrm{INF}_2) = 0`$.
3. **3-path** $`H_3(m; \{a, b\})`$: internal middle $`m \in M_3`$ ($`D_m < d_{3m}, 2 \le c_m \le K_{3m}`$); endpoints $`\{a, b\}`$ evaluated over $`\mu_a^{(3)} \otimes \mu_b^{(3)}`$.
4. **4-path** $`P_4(a - b - c - d)`$: palindrome-symmetric, internal edge $`(b, c)`$ with $`D_b, D_c < d_{4i}`$, $`2 \le c_b, c_c \le K_{4i}`$, $`D_b + D_c \le S_{P4}`$; endpoints $`(a, d)`$ over $`\mu_a^{(4)} \otimes \mu_d^{(4)}`$.
5. **3-star** $`S_3(o; \{a, b, c\})`$: hub $`o \in M_{S3}`$ ($`D_o < d_{S3}, 3 \le c_o \le K_{S3}`$); endpoints $`\{a, b, c\}`$ over $`(\mu^{(4)})^{\otimes 3}`$.
6. **5-path** $`P_5(a - b - c - d - e)`$: palindrome-symmetric, internal 3-path $`(b, c, d)`$ with $`D_b, D_c, D_d < d_{5i}`$, $`2 \le c_b, c_c, c_d \le K_{5i}`$, $`D_b + D_c \le S_{P4}`$, $`D_c + D_d \le S_{P4}`$, $`D_b + D_c + D_d \le S_{P5}`$; endpoints $`(a, e)`$ over $`\mu_a^{(5)} \otimes \mu_e^{(5)}`$.
7. **5-broom** $`B_5(o; \{a, b\}; c - d)`$: hub $`o`$ ($`D_o < d_{5o}, 3 \le c_o \le K_{5o}`$) adjacent to $`\{a, b\}`$ and arm $`c`$ ($`D_c < d_{5c}, 2 \le c_c \le K_{5c}`$) with $`D_o + D_c \le S_{B5}`$; endpoints $`\{a, b\}`$ and $`d`$ over $`(\mu^{(5)})^{\otimes 3}`$.
8. **4-star** $`S_4(o; \{a, b, c, d\})`$: hub $`o \in M_{S4}`$ ($`D_o < d_{S4}, 4 \le c_o \le K_{S4}`$); endpoints $`\{a, b, c, d\}`$ over $`(\mu^{(5)})^{\otimes 4}`$.

The header validation (`q5_validate_header`) enforces strict nesting of all internal and endpoint cutoffs:
```math
M_3 \subseteq G_2^{\mathrm{ex}}, \quad M_{4P} \subseteq (M_3 \cap G_3^{\mathrm{ex}})^2, \quad M_{S3} \subseteq M_3, \quad M_{5P} \subseteq (M_{4v} \cap G_4^{\mathrm{ex}})^3, \quad M_{5B} \subseteq (M_{4v} \cap M_{S3}) \times (M_{4v} \cap G_4^{\mathrm{ex}}), \quad M_{S4} \subseteq M_{S3}.
```
**Downward closure & exact newborn bijection:** Every core event increases either $`D_v`$ (gain) or $`c_v`$ (promotion), so an existing motif embedding can never enter a higher support (it can only stay in support or drop to $`0`$). Conversely, when a vertex $`u`$ promotes and spawns $`J = (2, 1)`$, deleting $`J`$ from any newly created supported $`k`$-vertex motif ($`2 \le k \le 5`$) yields a **unique** pre-existing supported $`(k-1)`$-vertex motif (because of the strict nesting above), and its newborn weight is allocated cleanly to that parent motif's promotion term.

---

## 3. Exact Transfers, Envelopes, and Local Corner Inequalities

Let $`\sigma = +1`$ (UPPER) or $`\sigma = -1`$ (LOWER). Every local inequality bounds $`\sigma S (\text{row drift} - \lambda \cdot \text{row weight})`$ from above in integer ticks.

1. **Paired inter-level transfers (exact cancellation):**
   - $`f \to g`$: At any exact $`G_2`$ vertex $`u`$ ($`D_u < d_2, 1 \le c_u \le K_2`$), replace $`u`$'s neighbor gain rate $`\sum_{v \sim u} \iota_v`$ in the $`f`$-row by $`c_u r_1(u)`$, transferring $`(\iota_v - r_1(u))\Delta_u^L f`$ to the incident edge $`g(u, v)`$.
   - $`g \to H_3`$ & Fenchel–Young (FY): At an exact $`G_2`$ vertex $`u`$ with $`c_u \ge 2`$:
     - If $`u \in M_3`$: replace the external donor rate $`\sum_{w \sim u, w \ne v} \iota_w`$ in $`g(u, v)`$ by $`(c_u - 1) r_1(u)`$, transferring $`(\iota_w - r_1(u))\Delta_u^L g(u, v)`$ to $`H_3(u; \{v, w\})`$.
     - If $`u \notin M_3`$: transfer a piecewise-linear function $`A_u(\iota_w)`$ (on rational knot grid $`0 = \xi_0 < \xi_1 = 1/d_2 < \dots < \xi_{m-1} = 1/2`$, denominator $`Q_K`$) from edge $`g(u, v)`$ to edge $`g(u, w)`$ for every $`w \ne v`$. Summing over all ordered pairs of incident edges at $`u`$ cancels identically. Thus $`g(u, v)`$ receives $`(c_u - 1) A_u(\iota_v)`$ (since $`A_u`$ is affine on $`[\xi_0, \xi_1] = [0, 1/d_2]`$, its expectation over $`\mu_v^{(2)}`$ is **exact** even when $`D_v \ge d_2`$) and keeps, in signed units,
       $`\displaystyle \sum_{w \ne v} \sigma\big(\iota_w \Delta_u^L g(u, v) - A_u(\iota_w)\big) \;\le\; (c_u - 1) \max_{0 \le k < m} \sigma\big(\xi_k \Delta_u^L g(u, v) - A_u(\xi_k)\big),`$
       because each segment $`[\xi_k, \xi_{k+1}]`$ is affine in $`\iota_w`$, and $`z \mapsto \max_k \sigma(\xi_k z - A_u(\xi_k))`$ is convex in $`z = \Delta_u^L g(u, v) = \mathbb E_{V \sim \mu_v^{(2)}}[\Delta_u^L g(u, V)]`$ so Jensen's inequality passes directly to the corners $`V \in G_2`$.
   - $`H_3 \to P_4, S_3`$ and $`P_4, S_3 \to P_5, B_5, S_4`$: Whenever adjoining an external core neighbor $`w`$ of an exact-inner motif vertex $`x`$ produces a supported motif of the next size, retain $`r_1(x)`$ per external donor in the lower-level row and transfer $`(\iota_w - r_1(x))\Delta_x^L h`$ to the unique higher-level motif. Because every higher-level motif has full $`G_{j+1}`$ tail endpoints, the support condition depends **only** on the internal vertices, never on the degree of $`w`$ or the other endpoints. Uncovered external donors at exact-inner vertices have $`\iota_w \in [0, 1/2]`$, bounded by extremizing over $`[0, \frac{1}{2} n_{\mathrm{ext}}]`$.
   All transfers occur in exact cancelling pairs, so summing over all embeddings on any tree $`T`$ preserves $`(\mathcal L \Phi)(T)`$ identically.

2. **Exact 1D, 2D, and 3D projection envelopes:**
   Each cross-level source change $`d_{\mathrm{fine}}`$ on $`G_j^k`$ ($`k \in \{1, 2, 3\}`$) is enclosed by stored integer arrays $`(m, r)`$ on $`G_{j+1}^k`$ with $`r \ge 0`$ verified over **every** fine state/pair/triple in exact `__int128` arithmetic (clearing integer denominators):
   $`\displaystyle \big| d_{\mathrm{fine}}(s) - (\Pi^{\otimes k} m)(s) \big| \;\le\; (\Pi^{\otimes k} r)(s).`$
   By independence of the donor coordinate from the remaining coordinates and convexity of $`|\cdot|`$,
   $`\displaystyle \sigma(\iota_w - r_1)\,\mathbb E[d_{\mathrm{fine}}] \;\le\; \mathbb E\big[\sigma(\mathrm{inv}_{j+1}(w) - r_1)\,m + |\mathrm{inv}_{j+1}(w) - r_1|\,r\big].`$

3. **Exact-degree overflow inner endpoints (`kind = 3`), one-knot tail endpoints (`kind = 1`), and double-tail / `INF` edges:**
   - For `kind = 3` ($`K_j < D_u < d_j`$, $`c_u \ge K_j + 1`$): $`\mathrm{inv} = 1/D_u`$ is exact; promotion leaves the state unchanged ($`\Delta_u^B = 0`$); its adjacent motif neighbor contributes exact rate $`\iota_{\mathrm{nbr}}`$ (when inner) or tail quadratic credit (charged once at the tail when the adjacent vertex is a tail); and its remaining gain rate $`L_u + \sum_{z \ne \mathrm{nbr}} \iota_z \le (D_u - c_u) + \frac{c_u - 1}{2} = D_u - \frac{c_u + 1}{2} \le D_u - \frac{K_j + 2}{2}`$ lies in $`[0, (2 D_u - K_j - 2)/2]`$.
   - For `kind = 1` ($`D_u \ge d_j`$, $`t = d_j / D_u`$): writing $`H_1, H_0`$ for $`T_j(\kappa)`$ and $`\mathrm{INF}_j`$, leaf gain contributes $`t \gamma (H_0 - H_1)`$ with $`\gamma = (L_u + \sum \iota_z)/(D_u + 1) \in [\gamma_{\min}, 1]`$, where $`\gamma_{\min} = \max(0, (d_j - \kappa + \iota_{\mathrm{nbr,lo}})/(d_j + 1))`$ for exact $`\kappa \le K_j`$ (and $`0`$ for overflow $`\kappa = K_j + 1`$); promotion contributes $`t b (H_1^{\mathrm{bud}} - H_1)`$ with $`b \in [\max(0, d_j - \kappa)/d_j, 1]`$ for $`\kappa \le K_j`$ (and zero weight change for overflow $`\kappa = K_j + 1`$, since promotion leaves the clamped class unchanged, contributing $`0`$ even when $`b > 0`$); and whenever the adjacent vertex $`y`$ is inner (`kind = 0` or `kind = 3`), $`y`$'s leaf gain from $`u`$ satisfies $`\sigma \frac{t}{d_j}((1-t)d_0 + t d_1) \le \frac{t}{d_j}\max(\sigma d_0, \sigma d_1)`$ (absorbing both the $`t`$ and $`1-t`$ parts of $`\Delta_y^L`$ into the $`t`$-row). When several distinct tail endpoints attach to the same internal hub, the same linear inequality applies independently to each coordinate after taking expectation over the others.
   - **Double-tail and `INF` corners in $`g(u, v)`$:** Only $`g(u, v)`$ has two adjacent endpoints that can simultaneously be physical tails ($`D_u, D_v \ge d_2`$, with $`t_u = d_2/D_u, t_v = d_2/D_v`$). Expanding $`g(u, v)`$ over $`\mu_u^{(2)} \otimes \mu_v^{(2)}`$ and using $`g(\mathrm{INF}_2, \mathrm{INF}_2) = 0`$: the $`(\mathrm{INF}_2, \mathrm{INF}_2)`$ corner carries zero value, zero incoming transfers, and zero residual; on $`(\mathrm{INF}_2, y)`$ with $`y`$ inner, $`y`$ still evolves under its own leaves and external neighbors while $`u`$'s credit $`\iota_u`$ is already accounted for in $`(T_u, y)`$. When both $`u`$ and $`v`$ have $`D_u, D_v \ge d_2`$, the cross-donor rates $`\iota_v`$ and $`\iota_u`$ are **not** dropped: they are fully contained inside $`r_u`$ and $`r_v`$ (hence inside $`\gamma_u, \gamma_v \in [(d_2 - \kappa)/(d_2 + 1), 1]`$ for exact classes $`\kappa \le K_2`$ — taking $`\iota_{\mathrm{nbr,lo}} = 0`$ since $`\iota > 0`$ — and inside $`[0, 1]`$ for overflow class $`\kappa = K_2 + 1`$), and averaging $`t_u \gamma_u (g(\mathrm{INF}_2, \cdot) - g(T_u, \cdot))`$ over $`\mu_v^{(2)}`$ assigns exact weights $`t_u t_v`$ to row $`(T_u, T_v)`$ and $`t_u(1 - t_v)`$ to row $`(T_u, \mathrm{INF}_2)`$ (and symmetrically for $`v`$).

4. **Analytic singleton tail ($`D \ge D_F`$) and seam ($`D = D_F - 1 \to D_F`$):**
   Writing $`q = L/D \in [0, 1]`$ and $`x = 1/D \in [0, 1/D_F]`$, the unsigned analytic potential in ticks is $`w(q, x) = f_{p0} + f_{p1} q`$ (UPPER) or $`w(q, x) = f_{p0} q - f_{p1} q^2 - f_{p2} q x`$ (LOWER). At the seam $`D = D_F - 1`$, `verify_quintet.cpp` evaluates the exact rational jump $`w((D_F - c)/D_F, 1/D_F) - f(D_F - 1, c)`$ directly in exact 128-bit arithmetic. For $`D \ge D_F`$, no supported higher motif can contain the promoted vertex as the internal neighbor of $`J`$, so promotion adds only $`f(2,1) + g((D, c+1), J)`$ with post-promotion class $`\min(c+1, K_2+1) \in \{1, \dots, K_2+1\}`$. With $`C_{\mathrm{nb}} = f(2,1) + g(\mathrm{INF}_2, J)`$ and $`d_g = d_2(\mathrm{ext}_{\kappa \in \{1,\dots,K_2+1\}} g(T_2(\kappa), J) - g(\mathrm{INF}_2, J))`$ ($`\max`$ for UPPER, $`\min`$ for LOWER), multiplying the signed residual at either normalized gain-rate endpoint $`\rho \in \{q, (1+q)/2\}`$ by $`(1+x)^2`$ yields an exact polynomial $`P_\rho(q, x)`$ enclosed on $`[0, 1] \times [0, 1/D_F]`$ by exact rational 2D tensor Bernstein coefficients across 16 equal $`q`$-panels (`verify_tail_bernstein.py`). Dividing a panel's Bernstein upper bound $`M_j`$ by $`(1+x)^2 \in [1, (1 + 1/D_F)^2]`$ gives the rigorous upper bound $`M_j`$ if $`M_j \ge 0`$ and $`M_j / (1 + 1/D_F)^2`$ if $`M_j < 0`$.

---

## 4. Per-Core Aggregation, Global Comparability, and Exponent Deduction

Let $`e_f := \max(\mathrm{cap}_{f,\mathrm{in}}, \mathrm{cap}_{f,\mathrm{tab}}, \mathrm{cap}_{f,\mathrm{ana}})`$, $`e_g := \max(0, \mathrm{cap}_g)`$, and let $`e_{\mathrm{high}} \ge 0`$ be the maximum over all finite states $`(D, c)`$ ($`D < d_{3m}`$) of the positive higher-motif caps allocated to a host core vertex ($`H_3 \to m`$; $`\frac{1}{2} P_4 \to b, c`$; $`S_3 \to o`$; $`P_5 \to c_{\mathrm{mid}}`$; $`B_5 \to o`$; $`S_4 \to o`$). Since $`T_K`$ has $`C(T)`$ vertices and $`C(T) - 1`$ edges, setting $`E_{\mathrm{tot}} := \max(0, e_f + e_g + e_{\mathrm{high}})`$ gives
```math
\sigma S \big((\mathcal L \Phi)(T) - \lambda \Phi(T)\big) \;\le\; E_{\mathrm{tot}}\, C(T) \qquad \text{for every finite tree } T \text{ with } C(T) \ge 1.
```

For comparability, choose any spanning rooting of $`T_K`$ (purely as an accounting bijection, leaving $`\Phi`$ unchanged): allocate $`f(u)`$ to $`u`$, $`g(p_u, u)`$ to child $`u \ne r`$ ($`0`$ to $`r`$), and each higher motif to its finite-degree host vertex as above. On any core tree $`T_K`$, $`\sum_{u \in K(T)} (c_u - 2) = -2`$, so for any $`\alpha \ge 0`$ (`alpha = 0` for UPPER; any validated $`\alpha \in [0, 1]`$ for LOWER, with $`\alpha = 1/4`$ in `gold1`):
```math
S\Phi(T) \;=\; 2\alpha S + \sum_{u \in K(T)} \big(W(u) + \alpha S (c_u - 2)\big) \;\ge\; m_{\min}\, C(T), \qquad S\Phi(T) \;\le\; m_{\max}\, C(T),
```
where $`m_{\min} > 0`$ and $`m_{\max} > 0`$ are verified over all $`D < D_F`$ (`verify_quintet.cpp`) and analytically over $`D \ge D_F`$ (`verify_tail_bernstein.py`, where no higher motif is hosted; any incident higher motif having a $`D \ge D_F`$ vertex as an outer endpoint is hosted by its finite-degree internal vertex and already accounted for in that host's finite comparability check).

Consequently, for every finite tree $`T`$ with $`C(T) \ge 1`$, setting $`\lambda_{\mathrm{up}}^* := \lambda_{\mathrm{up}} + E_{\mathrm{tot}}/m_{\min}`$ and $`\lambda_{\mathrm{lo}}^* := \lambda_{\mathrm{lo}} - E_{\mathrm{tot}}/m_{\min}`$:
```math
(\mathcal L \Phi_{\mathrm{up}})(T) \;\le\; \lambda_{\mathrm{up}}^*\,\Phi_{\mathrm{up}}(T), \qquad (\mathcal L \Phi_{\mathrm{lo}})(T) \;\ge\; \lambda_{\mathrm{lo}}^*\,\Phi_{\mathrm{lo}}(T).
```

**Non-explosion, stopped Dynkin, and Grönwall deduction:** Starting from the 3-vertex tree $`T_3`$ ($`C(T_3) = 1`$), the total vertex count $`N(t) = |V(T_t)|`$ increases by $`1`$ at total rate $`N(t)`$, so $`N(t)`$ is a non-explosive rate-$`1`$ Yule process with $`\mathbb E[N(t)^2] < \infty`$ for every $`t \ge 0`$, and $`C(T_t) \le N(t)`$. Let $`\tau_M := \inf\{t \ge 0 : N(t) \ge M\}`$. For either potential $`\Phi \in \{\Phi_{\mathrm{lo}}, \Phi_{\mathrm{up}}\}`$, comparability gives $`0 < (m_{\min}/S) C(T) \le \Phi(T) \le (m_{\max}/S) N(T)`$, and every single jump increases $`N`$ by $`1`$ (total jump rate $`N(T)`$), so
```math
|(\mathcal L \Phi)(T)| \;\le\; N(T)\big(\tfrac{m_{\max}}{S}(N(T) + 1) + \tfrac{m_{\max}}{S} N(T)\big) \;=\; \tfrac{m_{\max}}{S} N(T)(2N(T) + 1).
```
On $`[0, t \wedge \tau_M]`$, the chain takes values in the finite set $`\{N \le M\}`$, so Dynkin's formula gives
```math
\mathbb E[\Phi(T_{t \wedge \tau_M})] \;=\; \Phi(T_3) + \mathbb E\left[\int_0^{t \wedge \tau_M} (\mathcal L \Phi)(T_s)\,ds\right].
```
Since $`t \wedge \tau_M \uparrow t`$ a.s., $`0 \le \Phi(T_{t \wedge \tau_M}) \le (m_{\max}/S) N(t) \in L^1`$, and $`\int_0^{t \wedge \tau_M} |(\mathcal L\Phi)(T_s)|\,ds \le (m_{\max}/S)\,t\,N(t)(2N(t)+1) \in L^1`$, dominated convergence as $`M \to \infty`$ yields $`\mathbb E[\Phi(T_t)] = \Phi(T_3) + \int_0^t \mathbb E[(\mathcal L\Phi)(T_s)]\,ds`$ for all $`t \ge 0`$. Hence $`t \mapsto \mathbb E[\Phi(T_t)]`$ is absolutely continuous and satisfies the differential Grönwall inequalities, giving for all $`t \ge 0`$:
```math
\Phi_{\mathrm{lo}}(T_3)\,e^{\lambda_{\mathrm{lo}}^* t} \;\le\; \mathbb E[\Phi_{\mathrm{lo}}(T_t)] \;\le\; \frac{m_{\max,\mathrm{lo}}}{S}\,\mathbb E[C(T_t)], \qquad \frac{m_{\min,\mathrm{up}}}{S}\,\mathbb E[C(T_t)] \;\le\; \mathbb E[\Phi_{\mathrm{up}}(T_t)] \;\le\; \Phi_{\mathrm{up}}(T_3)\,e^{\lambda_{\mathrm{up}}^* t}.
```
Taking $`t^{-1}\log(\cdot)`$ unconditionally proves, from the 3-vertex seed $`T_3`$,
```math
\lambda_{\mathrm{lo}}^* \;\le\; \liminf_{t\to\infty} \frac{1}{t}\log \mathbb E_{T_3}[C(T_t)] \;\le\; \limsup_{t\to\infty} \frac{1}{t}\log \mathbb E_{T_3}[C(T_t)] \;\le\; \lambda_{\mathrm{up}}^*.
```
Passing to the original single-edge seed $`T_2`$ follows directly by conditioning on the first jump time $`\tau_1 \sim \mathrm{Exp}(2)`$ (which transitions deterministically from $`T_2`$ to $`T_3`$, with $`C(T_s) = 0`$ for $`s < \tau_1`$):
```math
\mathbb E_{T_2}[C(T_t)] \;=\; \int_0^t 2 e^{-2s}\,\mathbb E_{T_3}[C(T_{t - s})]\,ds,
```
so the upper bound multiplies by $`2/(2 + \lambda_{\mathrm{up}}^*)`$ and restricting the lower integral to $`s \in [0, 1]`$ (for $`t \ge 1`$) gives a positive constant multiple of $`e^{\lambda_{\mathrm{lo}}^* t}`$, preserving the identical $`\liminf/\limsup`$ bounds from $`T_2`$. (In discrete size $`n \ge 3`$, $`(\mathcal L\Phi)(T_n) = n\,\mathbb E[\Phi(T_{n+1}) - \Phi(T_n) \mid T_n]`$ directly gives $`\Phi_{\mathrm{lo}}(T_3)\prod_{j=3}^{n-1}(1 + \lambda_{\mathrm{lo}}^*/j) \le \mathbb E[\Phi_{\mathrm{lo}}(T_n)]`$ when $`\lambda_{\mathrm{lo}}^* > 0`$ — and when $`\lambda_{\mathrm{lo}}^* \le 0`$, $`C(T_n) \ge 1`$ trivially yields $`\liminf_{n\to\infty}(\log n)^{-1}\log\mathbb E[C(T_n)] \ge 0 \ge \lambda_{\mathrm{lo}}^*`$ — together with $`\mathbb E[\Phi_{\mathrm{up}}(T_n)] \le \Phi_{\mathrm{up}}(T_3)\prod_{j=3}^{n-1}(1 + \lambda_{\mathrm{up}}^*/j)`$, yielding the identical unconditional $`\liminf/\limsup`$ bounds for $`(\log n)^{-1}\log\mathbb E[C(T_n)]`$ as well.) Combining these unconditional drift bounds with the real-variable existence of the mean growth exponent $`\mu = \lim_{t\to\infty} t^{-1}\log \mathbb E[C(T_t)] = \lim_{n\to\infty} (\log n)^{-1}\log \mathbb E[C(T_n)]`$ (`proofs/existence.md` Theorem 1, Astra) yields the certified enclosure $`\lambda_{\mathrm{lo}}^* \le \mu \le \lambda_{\mathrm{up}}^*`$.
