# Rigorous Certified Enclosure of $`\mu`$ via Unified Bernstein-3 Full-Tail Level-4 Motif Potentials (`B3-FullH + L4`)

**Author:** Gemini 4 Neon (Neon lead). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Claim pins:** Neon 1 `12202c21a0eb134f754f5f56b77ce1a036921fc7`; Neon 3 `51f4cfadac675fdf1cef84d44935eebb5456aabb`.

**Review:** GPT-6.2 Sol (Astra and separate cold review) and Claude (Fable) reproduced scoped checks and replays; neither review provides a second full local-row verifier. For exact scope and dependencies, including the existence theorem and the archived GMP witness evidence, see [Neon Claims 1 and 3](../README.md). Source paths are relative to the archived lane root (`round2/proofs/neon`).

## 1. Statement of the Theorem

Let $`(T_n)_{n \ge 2}`$ be the random friend tree (Krapivsky–Redner undirected isotropic redirection): $`T_2`$ is a single edge, $`T_3`$ is the path on $`3`$ vertices rooted at its unique degree-$`2`$ core vertex $`o`$, and at each step $`n \to n+1`$, a vertex $`V_n`$ is chosen uniformly from $`T_n`$, a neighbour $`W_n`$ is chosen uniformly from the neighbours of $`V_n`$, and a new leaf $`n+1`$ is attached to $`W_n`$. Let $`C_n = \#\{v \in T_n : \deg(v) \ge 2\}`$ be the number of non-leaf (core) vertices, and let $`\mu = \lim_{n \to \infty} \frac{\log \mathbb{E}[C_n]}{\log n}`$ (whose existence and almost-sure equality to $`\lim_{n\to\infty} \frac{\log C_n}{\log n}`$ are established in `proofs/existence.md` and `proofs/almost-sure-exponent.md`).

For any finite tree $`T`$ with at least one core vertex ($`C(T) \ge 1`$) whose core subtree $`K(T)`$ is rooted at a core vertex $`o \in K(T)`$, let $`\mathcal{L}`$ be the unnormalized attachment generator:
```math
(\mathcal{L} F)(T) := \sum_{w \in V(T)} r_T(w)\,\big(F(T +_w \text{leaf}) - F(T)\big), \qquad r_T(w) := \sum_{v \sim_T w} \frac{1}{\deg_T(v)},
```
so that $`\mathbb{E}[F(T_{n+1}) - F(T_n) \mid T_n] = \frac{1}{n}(\mathcal{L} F)(T_n)`$ identically for every $`n \ge 3`$.

**Theorem 1 (B3-FullH Drift Certificate).**
Suppose a certificate $`\mathcal{C} = (\text{mode}, P, Q, S, S_r, \text{params}, \text{tables})`$ satisfies the machine-range header and entry limits enforced by `b3f_validate_header` and `read_i64_array` (guaranteeing overflow-free pre-GMP integer preparation) and passes all exact arbitrary-precision integer/rational checks of `round2/cert/cert_b3f_verify_gmp.c`. Then the associated potential $`\Phi_{\mathcal{C}}(T)`$ satisfies:
1. *(Drift inequality)* For every finite tree $`T`$ with $`C(T) \ge 1`$ and core rooted at a core vertex:
   - If $`\text{mode} = 1`$ (UPPER): $`(\mathcal{L}\Phi_{\mathcal{C}})(T) \le \frac{P}{Q}\,\Phi_{\mathcal{C}}(T)`$.
   - If $`\text{mode} = 2`$ (LOWER): $`(\mathcal{L}\Phi_{\mathcal{C}})(T) \ge \frac{P}{Q}\,\Phi_{\mathcal{C}}(T)`$.
2. *(Core comparison)* There exist exact positive rationals $`c_0, C_0, \Phi_{\mathcal{C}}(T_3)`$ computed in GMP such that for every finite tree $`T`$ with $`C(T) \ge 1`$ and core rooted at a core vertex:
   - For UPPER ($`\text{mode} = 1`$): $`c_0 > 0`$, $`\Phi_{\mathcal{C}}(T_3) > 0`$, and $`c_0\,C(T) \le \Phi_{\mathcal{C}}(T) \le C_0\,C(T)`$, hence for all $`n \ge 3`$,
     $`\displaystyle \mathbb{E}[C_n] \le \frac{\Phi_{\mathcal{C}}(T_3)}{c_0} \prod_{j=3}^{n-1} \left(1 + \frac{P/Q}{j}\right), \qquad \text{and therefore } \mu \le \frac{P}{Q}.`$
   - For LOWER ($`\text{mode} = 2`$): $`C_0 > 0`$, $`\Phi_{\mathcal{C}}(T_3) > 0`$, and $`\Phi_{\mathcal{C}}(T) \le C_0\,C(T)`$, hence for all $`n \ge 3`$,
     $`\displaystyle \mathbb{E}[C_n] \ge \frac{\Phi_{\mathcal{C}}(T_3)}{C_0} \prod_{j=3}^{n-1} \left(1 + \frac{P/Q}{j}\right), \qquad \text{and therefore } \mu \ge \frac{P}{Q}.`$

**Corollary (Verified Production Enclosure `B3-FullH-18`).**
For the production parameter set $`(D_0, K_T, L_F, K_F) = (18, 17, 14000, 180)`$ and Level-4 inner cutoffs $`P_4=(10,7)`$, $`Y_{\mathrm{top}}=(9,6)`$, $`Y_{\mathrm{bot}}=(10,7)`$, $`S_3=(8,6)`$ ($`436{,}870{,}169`$ stored potential-table entries, $`436{,}886{,}484`$ total `int64_t` entries including rate/transfer arrays, $`3{,}495{,}091{,}968`$ bytes per certificate file, $`S = 10^{11}`$, $`S_r = 10^9`$), the pure-GMP verifier `round2/cert/cert_b3f_verify_gmp.c` passes with `0` failures on both certificates:
- **Lower certificate** (`b3f18_lo.bin`, SHA-256 `9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec`):
  $`\displaystyle \lambda_{\mathrm{lo}} = \frac{565673}{1000000} = 0.565673, \qquad C_0 = \frac{184124847315}{100000000000}, \qquad \Phi_{\mathrm{lo}}(T_3) = \frac{146042017510}{100000000000},`$
  giving for every $`n \ge 3`$:
  $`\displaystyle \mathbb{E}[C_n] \ge \frac{146042017510}{184124847315} \prod_{j=3}^{n-1} \left(1 + \frac{565673/1000000}{j}\right) \qquad (> 0.793168 \prod_{j=3}^{n-1} (1 + 0.565673/j)).`$
- **Upper certificate** (`b3f18_hi.bin`, SHA-256 `1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41`):
  $`\displaystyle \lambda_{\mathrm{hi}} = \frac{565775}{1000000} = \frac{22631}{40000} = 0.565775, \qquad c_0 = \frac{115462837236}{400000000000}, \qquad \Phi_{\mathrm{hi}}(T_3) = \frac{146044806770}{100000000000},`$
  giving for every $`n \ge 3`$:
  $`\displaystyle \mathbb{E}[C_n] \le \frac{584179227080}{115462837236} \prod_{j=3}^{n-1} \left(1 + \frac{22631/40000}{j}\right) \qquad (< 5.059457 \prod_{j=3}^{n-1} (1 + 0.565775/j)).`$
Consequently, without using any external prior bracket,
```math
\frac{565673}{1000000} \;\le\; \mu \;\le\; \frac{22631}{40000}, \qquad \text{i.e., } \mu \in [0.565673,\; 0.565775] \text{ (width } 0.000102\text{)}.
```

**Corollary 2 (Verified Production Enclosure `B3-FullH-20-L4`, Sub-$`10^{-4}`$ Width).**
For the scaled production parameter set $`(D_0, K_T, L_F, K_F) = (20, 19, 16000, 200)`$ and Level-4 inner cutoffs $`P_4=(12,9)`$, $`Y_{\mathrm{top}}=(11,8)`$, $`Y_{\mathrm{bot}}=(12,9)`$, $`S_3=(10,7)`$ ($`1{,}398{,}838{,}433`$ stored potential-table entries across the eight motif arrays, $`1{,}398{,}860{,}120`$ total `int64_t` entries including the $`21{,}687`$ rate/transfer entries, $`11{,}190{,}881{,}056`$ bytes per certificate file, $`S = 10^{11}`$, $`S_r = 10^9`$), the same hardened pure-GMP verifier `round2/cert/cert_b3f_verify_gmp.c` passes with `0` failures on both certificates (recorded in `round2/outputs/cert_b3f_20_{lo,hi}.log`):
- **Lower certificate** (`b3f20_lo.bin`, SHA-256 `eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469`):
  $`\displaystyle \lambda_{\mathrm{lo}} = \frac{565680}{1000000} = \frac{7071}{12500} = 0.565680, \qquad C_0 = \frac{184219674916}{100000000000}, \qquad \Phi_{\mathrm{lo}}(T_3) = \frac{146041684821}{100000000000},`$
  giving for every $`n \ge 3`$:
  $`\displaystyle \mathbb{E}[C_n] \ge \frac{146041684821}{184219674916} \prod_{j=3}^{n-1} \left(1 + \frac{7071/12500}{j}\right) \qquad (> 0.792758 \prod_{j=3}^{n-1} (1 + 0.565680/j)).`$
- **Upper certificate** (`b3f20_hi.bin`, SHA-256 `55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2`):
  $`\displaystyle \lambda_{\mathrm{hi}} = \frac{565768}{1000000} = \frac{70721}{125000} = 0.565768, \qquad c_0 = \frac{115454594604}{400000000000}, \qquad \Phi_{\mathrm{hi}}(T_3) = \frac{146044513514}{100000000000},`$
  giving for every $`n \ge 3`$:
  $`\displaystyle \mathbb{E}[C_n] \le \frac{584178054056}{115454594604} \prod_{j=3}^{n-1} \left(1 + \frac{70721/125000}{j}\right) \qquad (< 5.059808 \prod_{j=3}^{n-1} (1 + 0.565768/j)).`$
Consequently, without using any external prior bracket,
```math
\frac{7071}{12500} \;\le\; \mu \;\le\; \frac{70721}{125000}, \qquad \text{i.e., } \mu \in [0.565680,\; 0.565768] \text{ (width } 0.000088 < 10^{-4}\text{)}.
```

---

## 2. State Space, Bernstein-3 Tail Embedding, and Potential $`\Phi(T)`$

Root the core tree $`K(T) = \{u \in V(T) : \deg_T(u) \ge 2\}`$ at the initial core vertex $`o`$. Each $`u \in K(T)`$ has:
- $`h_u \in \{0, 1\}`$ ($`0`$ if $`u = o`$, $`1`$ otherwise),
- $`k_u \ge 0`$ core children,
- $`L_u \ge 0`$ leaf children,
- total degree $`D_u = L_u + k_u + h_u \ge 2`$, and core degree $`c_u = k_u + h_u \ge 0`$.

Fix integer parameters $`D_0 \ge 4`$, $`K_T \ge D_0 - 1`$, $`L_F \ge D_0`$, $`K_F \ge K_T + 1`$, and per-motif inner cutoffs $`(D_{P4}, K_{P4})`$, $`(D_{YT}, K_{YT})`$, $`(D_{YB}, K_{YB})`$, $`(D_{S3}, K_{S3})`$ all $`\le (D_0, K_T)`$.

### 2.1 One-Body State Map $`\pi_f`$ and Bernstein-3 Tail Embedding $`\pi_g`$
- **One-body map $`\pi_f(u)`$**: maps $`(L_u, k_u, h_u)`$ to $`\mathrm{VSTAR}_{h_u}`$ if $`k_u > K_F`$, and to $`(\min(L_u, L_F), k_u, h_u)`$ if $`k_u \le K_F`$.
- **Bernstein-3 embedding $`\pi_g(u)`$**:
  - If $`D_u < D_0`$ (*inner*): $`\pi_g(u)`$ is the Dirac mass $`1`$ on the exact state $`(L_u, k_u, h_u)`$.
  - If $`D_u \ge D_0`$ (*tail*): set $`t_u = D_0 / D_u \in (0, 1]`$, $`s_u = 1 - t_u \in [0, 1)`$, and $`\kappa_u = \min(k_u, K_T + 1)`$. Then $`\pi_g(u)`$ is the convex probability combination of 4 tail states:
    $`\displaystyle \pi_g(u) = s_u^3\,\mathrm{INF} + 3 t_u s_u^2\,T_1(\kappa_u, h_u) + 3 t_u^2 s_u\,T_2(\kappa_u, h_u) + t_u^3\,T_3(\kappa_u, h_u).`$
  Assign each basis state $`x \in \mathcal{G}`$ the linear inverse-degree functional
  $`\displaystyle \mathrm{invD}(x) = \begin{cases} 1/D_x & x \text{ inner}, \\ 0 & x = \mathrm{INF}, \\ 1/(3 D_0) & x = T_1(\kappa, h), \\ 2/(3 D_0) & x = T_2(\kappa, h), \\ 1/D_0 & x = T_3(\kappa, h). \end{cases}`$

**Lemma 1 (Exact linear reproduction of $`1/D_u`$).**
For every core vertex $`u \in K(T)`$, $`\sum_{x} \pi_g(u)_x\,\mathrm{invD}(x) = 1/D_u`$ identically.
*Proof.* Immediate for $`D_u < D_0`$. For $`D_u \ge D_0`$,
```math
s_u^3 \cdot 0 + 3 t_u s_u^2 \cdot \frac{1}{3 D_0} + 3 t_u^2 s_u \cdot \frac{2}{3 D_0} + t_u^3 \cdot \frac{1}{D_0} = \frac{t_u}{D_0}\big((1-t_u)^2 + 2 t_u(1-t_u) + t_u^2\big) = \frac{t_u}{D_0} = \frac{1}{D_u}. \qquad \square
```

### 2.2 The Eight Motif Tables Defining $`\Phi(T)`$
Every table value is an exact rational $`q / S`$ ($`q \in \mathbb{Z}`$, $`S > 0`$). Let $`\mathcal{H}_u`$ be the set of all inner AND tail (`T1, T2, T3`) states with $`h = 1, 1 \le k \le K_T`$; let $`\mathcal{C}_u`$ be all inner AND tail (`T1, T2, T3`) states with $`2 \le k \le K_T`$. Note that $`\mathrm{INF} \notin \mathcal{H}_u \cup \mathcal{C}_u`$ (implicitly zero), and since $`K_T \ge D_0 - 1`$, **every** inner state with $`h=1, k \ge 1`$ belongs to $`\mathcal{H}_u`$ and **every** inner state with $`k \ge 2`$ belongs to $`\mathcal{C}_u`$. We require $`g(\mathrm{INF}, \mathrm{INF}) = 0`$.

Define $`\Phi(T) = \Phi_f + \Phi_g + \Phi_{hV} + \Phi_{hC} + \Phi_{P4} + \Phi_{Yt} + \Phi_{Yb} + \Phi_{S3}`$ by multilinearly extending each table via $`\pi_g`$ on every slot:
1. $`\Phi_f(T) = \sum_{u \in K(T)} f(\pi_f(u))`$.
2. $`\Phi_g(T) = \sum_{(p \to w) \in E(K)} g(\pi_g(p), \pi_g(w))`$.
3. $`\Phi_{hV}(T) = \sum_{\substack{p \to u \to w \\ 1 \le k_u \le K_T}} h_V(\pi_g(p), \pi_g(u)|_{\mathcal{H}_u}, \pi_g(w))`$ (full-tail middle $`u`$).
4. $`\Phi_{hC}(T) = \sum_{\substack{u \in K,\, 2 \le k_u \le K_T \\ \{w_1, w_2\} \subset \mathrm{ch}(u)}} h_C(\pi_g(u)|_{\mathcal{C}_u}; \pi_g(w_1), \pi_g(w_2))`$ (full-tail apex $`u`$, symmetric in $`\{w_1, w_2\}`$).
5. $`\Phi_{P4}(T) = \sum_{a \to b \to c \to d} q_{P4}(\pi_g(a), b, c, \pi_g(d))`$ for inner $`b, c \in \mathcal{D}_{P4}`$ ($`h=1, 1 \le k \le K_{P4}, D < D_{P4}`$).
6. $`\Phi_{Yt}(T) = \sum_{p \to u;\, \{w_1,w_2\}\subset\mathrm{ch}(u)} q_{Yt}(\pi_g(p), u; \pi_g(w_1), \pi_g(w_2))`$ for inner $`u \in \mathcal{D}_{YT}`$ ($`h=1, 2 \le k \le K_{YT}, D < D_{YT}`$).
7. $`\Phi_{Yb}(T) = \sum_{u \to w_1;\, u \to w_2 \to c} q_{Yb}(u, w_2; \pi_g(w_1), \pi_g(c))`$ for inner $`u \in \mathcal{D}_{YBU}, w_2 \in \mathcal{D}_{YBW}`$ ($`w_1 \ne w_2`$).
8. $`\Phi_{S3}(T) = \sum_{u;\, \{w_1,w_2,w_3\}\subset\mathrm{ch}(u)} q_{S3}(u; \pi_g(w_1), \pi_g(w_2), \pi_g(w_3))`$ for inner $`u \in \mathcal{D}_{S3}`$ ($`3 \le k \le K_{S3}, D < D_{S3}`$).

---

## 3. Exact Generator Decomposition and Zero-Sum Telescoping Transfers

Under $`\mathcal{L}`$, only two event types change $`K(T)`$:
- **Leaf gain at $`u \in K(T)`$**: rate $`L_u + \rho_u`$, where $`\rho_u = \sum_{v \sim_K u} 1/D_v`$; replaces $`(L_u, k_u, h_u) \mapsto (L_u + 1, k_u, h_u)`$ ($`D_u \mapsto D_u + 1`$).
- **Leaf promotion (budding) at $`u \in K(T)`$** ($`L_u \ge 1`$): rate $`b_u = L_u / D_u`$; replaces $`(L_u, k_u, h_u) \mapsto (L_u - 1, k_u + 1, h_u)`$ ($`D_u`$ unchanged!) and creates a newborn core child $`v_{\mathrm{nb}}`$ of $`u`$ in state $`(1, 0, 1)`$.

Notice that because $`D_u`$ is unchanged on promotion, $`t_u = D_0 / D_u`$ is invariant on promotion: if $`D_u \ge D_0`$ and $`k_u < K_T`$ (so the post-bud core-child count $`k_u + 1 \le K_T`$ remains in $`\mathcal{H}_u, \mathcal{C}_u`$), the newborn `hC(u; w, v_nb)` (for each existing core child $`w`$ of $`u`$) and newborn `hV(p, u, v_nb)` (when $`h_u = 1`$) carry the **exact same** Bernstein-3 weights $`(3 t_u s_u^2, 3 t_u^2 s_u, t_u^3)`$ on $`(T_1, T_2, T_3)`$ as the parent/child edge $`g(u, w)`$ or $`g(p, u)`$. At the boundary $`k_u = K_T`$, the post-bud state $`k_u + 1 = K_T + 1`$ exits $`\mathcal{H}_u \cup \mathcal{C}_u`$ (no surviving newborn restricted triple is created, while the edge $`g`$ transitions to $`\kappa = K_T + 1`$).

### 3.1 Exact Zero-Sum Telescoping Transfers
Fix any rational centering constants $`\bar r_{\mathrm{par}}(u), \bar r_{\mathrm{ch}}(u)`$ (for inner $`u`$) and $`\bar r^{\mathrm{tail}}_{\mathrm{par}}(\kappa, h), \bar r^{\mathrm{tail}}_{\mathrm{ch}}(\kappa, h)`$ (for tail $`0 \le \kappa \le K_T`$), plus any rational arrays $`c_2(\kappa, h), c_3(\kappa, h)`$, $`Y_{p,\mathrm{ref}}(\kappa, h, w)`$, $`Y_{w,\mathrm{ref}}(p, \kappa)`$.

1. **Level 1 $`\to`$ Level 2 (`f -> g`)**:
   - For inner $`u`$ ($`D_u < D_0`$): write $`\rho_u = (h_u \bar r_{\mathrm{par}}(u) + k_u \bar r_{\mathrm{ch}}(u)) + \sum_{v \sim_K u} (\frac{1}{D_v} - \bar r_{u,v})`$. Keep the centered base rate in $`f(u)`$ and transfer $`(\frac{1}{D_v} - \bar r_{u,v}) \Delta_L f(u)`$ to edge $`g`$ incident to $`\{u, v\}`$. By Lemma 1, $`\frac{1}{D_v} - \bar r_{u,v} = \sum_y \pi_g(v)_y (\mathrm{invD}(y) - \bar r_{u,v})`$.
   - For tail $`u`$ ($`D_u \ge D_0, k_u \le K_T`$): set $`\psi_u = \frac{c_2(k_u,h_u)}{D_u^2} + \frac{c_3(k_u,h_u)}{D_u^3}`$. Since
     $`\displaystyle 3 t_u^2 s_u \cdot \frac{c_2}{3 D_0^2} + t_u^3 \left(\frac{c_2}{D_0^2} + \frac{c_3}{D_0^3}\right) = \frac{c_2}{D_u^2}(1 - t_u + t_u) + \frac{c_3}{D_u^3} = \psi_u,`$
     adding $`(\bar R_u - \rho_u)\psi_u`$ to $`f(u)`$ and $`(\mathrm{invD}(y) - \bar r_{u,v})\frac{c_2}{3 D_0^2}`$ to `g(T2)` and $`(\mathrm{invD}(y) - \bar r_{u,v})(\frac{c_2}{D_0^2} + \frac{c_3}{D_0^3})`$ to `g(T3)` sums identically to $`0`$ on every tree!
2. **Level 2 $`\to`$ Level 3 (`g -> hV, hC`, BOTH inner AND tail $`u`$ with $`k_u \le K_T`$)**:
   - For inner $`u`$ ($`D_u < D_0`$): since $`K_T \ge D_0 - 1`$, every external neighbour $`z \sim_K u`$ ($`z \ne w`$) forms an `hV` or `hC` motif with middle/apex $`u`$, so $`(\frac{1}{D_z} - \bar r_{u,z}) \Delta_L g`$ transfers completely to `hV`/`hC`.
   - For tail $`u`$ ($`D_u \ge D_0, k_u \le K_T`$) on edge $`g(u, w)`$ (or $`g(p, u)`$):
     - On `T2` (weight $`3 t_u^2 s_u`$): subtracting $`3 t_u^2 s_u (\rho_{u,\mathrm{ext}} - \bar R_{u,\mathrm{ext}})\frac{F_0 - F_1}{D_0}`$ from `g` and adding $`(\mathrm{invD}(z) - \bar r_{u,z})\frac{F_0 - F_1}{D_0}`$ to `hV(*, T2(u), *)` and `hC(T2(u); *, *)` cancels identically!
     - On `T3` (weight $`t_u^3`$): subtracting $`t_u^3 (\rho_{u,\mathrm{ext}} - \bar R_{u,\mathrm{ext}}) Y_{\mathrm{ref}}`$ from `g` and adding $`(\mathrm{invD}(z) - \bar r_{u,z}) Y_{\mathrm{ref}}`$ to `hV(*, T3(u), *)` and `hC(T3(u); *, *)` cancels identically for EVERY choice of $`Y_{\mathrm{ref}}`$!
3. **Level 3 $`\to`$ Level 4 (`hV, hC -> P4, Ytop, Ybot, S3`)**:
   - Whenever an inner motif vertex in `hV` or `hC` lies in the corresponding inner Level-4 cutoff domain, its external neighbour credit $`(\frac{1}{D_z} - \bar r)`$ is transferred identically to `P4`, `Ytop`, `Ybot`, or `S3`.

Because every transfer pair adds $`+X`$ to one motif and $`-X`$ to its neighbour motif, **the sum of all transferred terms over any tree $`T`$ is identically zero**.

---

## 4. Exact Bernstein-3/4 Tail Enclosures

### 4.1 Tail Vertex with $`k \le K_T`$
For a tail vertex $`x`$ with degree $`D \ge D_0`$, $`c = k + h`$, $`L = D - c \ge 0`$, $`t = D_0 / D \in (0, t_{\max}]`$ where $`t_{\max} = \min(1, D_0 / c)`$, $`s = 1 - t`$, $`t' = D_0 / (D + 1) = \frac{D_0 t}{D_0 + t}`$, leaf-gain rate $`L + \rho`$ ($`\rho = \rho_{\mathrm{in}} + \rho_{\mathrm{ext}}`$), bud rate $`L/D = 1 - c/D`$, and $`n_{\mathrm{in}} \in \{0, 1, 2\}`$ inner motif neighbours receiving credit $`1/D`$ from $`x`$ with values $`(d_0, d_1, d_2, F_{3,\mathrm{up}} - n_{\mathrm{in}} F_3)`$, direct rational algebra (independently verified in exact Python `fractions.Fraction`) gives the **exact identity**:
```math
(L + \rho)\big(B_3(F; t') - B_3(F; t)\big) + \frac{L}{D}\big(B_3(F_b; t) - B_3(F; t)\big) + \frac{1}{D} B_3(d; t) - 3 t^2 s (\rho_{\mathrm{ext}} - \bar R_{\mathrm{ext}})\frac{F_0 - F_1}{D_0} - t^3 (\rho_{\mathrm{ext}} - \bar R_{\mathrm{ext}}) Y_{\mathrm{ref}}
```
```math
= 3 t s^2 C_{T1} + 3 t^2 s C_{T2}(\rho_{T2}) + t^3 C_{T3}(t, \rho_{\mathrm{in}}, \rho_{\mathrm{ext}}),
```
where $`F_{0b} = F_0`$ (the $`\mathrm{INF}`$ value is unchanged on promotion), $`\rho_{T2} = \rho_{\mathrm{in}} + \bar R_{\mathrm{ext}}`$ when the external-credit transfer is active (`g` with $`k \le K_T`$) and $`\rho_{T2} = \rho_{\mathrm{in}} + \rho_{\mathrm{ext}}`$ otherwise (when no transfer is subtracted and $`Y_{\mathrm{ref}} = 0`$), and:
- $`C_{T1} = (F_0 + F_{1b} - 2 F_1) + \frac{d_0}{3 D_0}`$;
- $`C_{T2}(\rho_{T2}) = 2 F_1 + F_{2b} - 3 F_2 - \frac{c}{D_0}(F_{1b} - F_1) - \frac{c + 1 - \rho_{T2}}{D_0}(F_0 - F_1) + \frac{3 d_1 - d_0}{3 D_0}`$ (affine in $`\rho_{T2}`$, hence extremal at its two endpoints);
- $`C_{T3}(t, \rho_{\mathrm{in}}, \rho_{\mathrm{ext}}) = \frac{P(t)}{(D_0 + t)^3} - (\rho_{\mathrm{ext}} - \bar R_{\mathrm{ext}}) Y_{\mathrm{ref}}`$, where $`B(t) = (D_0 + t)^3 > 0`$ and $`P(t)`$ is the explicit degree-$`\le 4`$ polynomial in $`t`$ (affine in $`(\rho_{\mathrm{in}}, \rho_{\mathrm{ext}})`$): setting $`a = F_0 - F_1`$, $`b = F_1 - F_2`$, $`p = F_2 - F_3`$, $`\delta_j = F_{jb} - F_j`$ ($`\delta_0 = 0`$), and $`x = c + 1 - (\rho_{\mathrm{in}} + \rho_{\mathrm{ext}})`$,
  $`\displaystyle P_2(t) = \big(3 D_0^2 p + 3 D_0(a - b)\big) + \big(3 D_0 p + (3 - 3 D_0)a + (6 D_0 - 3)b\big)t + (p - 2a + 4b)t^2,`$
  $`\displaystyle P(t) = (D_0 + t - x t) P_2(t) + s x\left(\frac{3 a(D_0 + 1)}{D_0} - 6 b\right)(D_0 + t)^2 + \left[\delta_3\left(1 - \frac{c t}{D_0}\right) + \frac{3 c s(\delta_1 - \delta_2)}{D_0} + s\frac{d_0 - 3 d_1 + 3 d_2}{D_0} + t\frac{F_{3,\mathrm{up}} - n_{\mathrm{in}} F_3}{D_0}\right](D_0 + t)^3.`$

Write $`u = t / t_{\max} \in (0, 1]`$ and expand both $`P(t) = \sum_{j=0}^4 \binom{4}{j} u^j (1-u)^{4-j} B_j(P)`$ and $`B(t) = \sum_{j=0}^4 \binom{4}{j} u^j (1-u)^{4-j} B_j(B)`$ in the degree-4 Bernstein basis on $`[0, t_{\max}]`$. Since $`B_j(B) > 0`$ for all $`j = 0, \dots, 4`$ and $`\binom{4}{j} u^j (1-u)^{4-j} \ge 0`$, for every $`t \in (0, t_{\max}]`$:
```math
\min_{0 \le j \le 4} \frac{B_j(P)}{B_j(B)} \le \frac{P(t)}{B(t)} \le \max_{0 \le j \le 4} \frac{B_j(P)}{B_j(B)}.
```
Taking extrema over the corners of $`(\rho_{\mathrm{in}}, \rho_{\mathrm{ext}})`$ (affine!) and $`j \in \{0, 1, 2, 3, 4\}`$ yields an exact finite rational bound valid uniformly for all $`D \ge \max(D_0, c)`$.

### 4.2 Tail Overflow ($`k > K_T`$)
When $`k > K_T`$, promotion leaves the uncut overflow slot label $`\kappa = K_T + 1`$ unchanged ($`\Delta_{\mathrm{bud}} = 0`$ on existing higher motifs; newborn `f` and `g` remain accounted for in the one-body equation for `f`). Because every core neighbour has degree $`\ge 2`$, $`\rho \in [0, c/2]`$, hence
```math
0 \le L + \rho = D - c + \rho \le D - \frac{c}{2} < D + 1.
```
Writing $`L + \rho = \alpha(D + 1)`$ with $`\alpha \in [0, 1)`$ and using the exact polynomial identity
```math
(D + 1)\big(B_3(F; t') - B_3(F; t)\big) = 3 t s^2 (F_0 - F_1) + 3 t^2 s \cdot 2(F_1 - F_2) + t^3 \frac{P_2(t)}{(D_0 + t)^2},
```
with the same $`P_2(t)`$ of degree $`\le 2`$ as above, the three terms are bounded uniformly for all $`D \ge D_0`$ and $`\alpha \in [0, 1]`$ by $`\mathrm{ext}(0, F_0 - F_1)`$ on `T1`, $`\mathrm{ext}(0, 2(F_1 - F_2))`$ on `T2`, and $`\mathrm{ext}(\{0\} \cup \{B_j(P_2)/B_j((D_0+t)^2) : j=0,1,2\})`$ on `T3`. Likewise, the inner-neighbour credit satisfies
```math
\frac{1}{D} B_3(d; t) = 3 t s^2 \frac{d_0}{3 D_0} + 3 t^2 s \frac{3 d_1 - d_0}{3 D_0} + t^3 \left[(1-t)\frac{d_0 - 3 d_1 + 3 d_2}{D_0} + t \frac{F_{3,\mathrm{up}} - n_{\mathrm{in}} F_3}{D_0}\right],
```
bounded on `T3` by the extrema of the two affine endpoints at $`t=0`$ and $`t=1`$.

### 4.3 One-Body Capped States ($`L \ge L_F`$ and $`k > K_F`$)
- At $`L \ge L_F, k \le K_F`$: leaf gain has $`\Delta_L f = 0`$; the $`f`$-tail transfer $`(\bar R_u - \rho_u)(c_2 x^2 + c_3 x^3)`$ for $`x = 1/D \in (0, x_0]`$, $`x_0 = 1/(L_F + k + h)`$, is bounded by the degree-3 Bernstein extrema of $`c_2 x^2 + c_3 x^3`$ on $`[0, x_0]`$, namely $`\{0, c_2 x_0^2 / 3, c_2 x_0^2 + c_3 x_0^3\}`$; budding separates $`L = L_F`$ (exact bud transition to $`L_F - 1`$, paired with the bounded transfer above) from $`L \ge L_F + 1`$ (where $`b \in [(L_F+1)/(L_F+1+k+h), 1]`$ and $`g(\pi_g(u'), v_{\mathrm{nb}})`$ is bounded by its 4 de Casteljau Bernstein-3 control points on $`[0, D_0/(L_F+1+k+h)]`$).
- At $`k > K_F`$ (`VSTAR`): $`\Delta_L f = 0`$, $`b \in [0, 1]`$, and $`g(\pi_g(u'), v_{\mathrm{nb}})`$ is bounded by its 4 de Casteljau control points on $`[0, D_0/(K_F+1+h)]`$.

---

## 5. Core Comparison Bounds and Code Map

Every motif in $`\Phi(T)`$ is charged to a unique core vertex $`v \in K(T)`$:
- `f(v)` to $`v`$ ($`1`$ term);
- `g(p, v)`, `hV(p, u, v)`, `qP4(a, b, c, v)` to bottom vertex $`v`$ (at most $`1`$ each, only if $`h_v = 1`$);
- `qYb(u; w1, w2 -> v)` to bottom vertex $`v`$ (at most $`1`$ group of $`k_u - 1`$ siblings sharing $`(u, w_2)`$, only if $`h_v = 1`$);
- `hC(v; w1, w2)` to apex $`v`$ ($`\binom{k_v}{2}`$ pairs when $`2 \le k_v \le K_T`$; $`0`$ otherwise);
- `qYt(p, v; w1, w2)` to apex $`v`$ ($`\binom{k_v}{2}`$ pairs when $`v \in \mathcal{D}_{YT}`$);
- `qS3(v; w1, w2, w3)` to apex $`v`$ ($`\binom{k_v}{3}`$ triples when $`v \in \mathcal{D}_{S3}`$).

Taking the exact integer maximum $`S C_0`$ of the per-vertex charge over all finitely many $`(h, k \le K_F, L \le L_F)`$ plus $`\mathrm{VSTAR}_0, \mathrm{VSTAR}_1`$ gives $`\Phi(T) \le C_0 C(T)`$ for every finite tree $`T`$ with $`C(T) \ge 1`$. For the lower comparison $`c_0 C(T) \le \Phi(T)`$, we apply the exact zero-sum tree-edge re-charging identity
```math
\sum_{u \in K(T)} \frac{1}{4}(k_u - h_u) = \frac{|E(K(T))| - |E(K(T))|}{4} = 0,
```
valid identically on every rooted tree $`T`$ (each core edge contributes $`+1`$ to $`k_p`$ and $`+1`$ to $`h_w`$). Since $`\frac{1}{4}(k_u - h_u) \ge \frac{1}{4}(K_F + 1 - h_u)`$ for every overflow vertex $`k_u > K_F`$, taking the exact integer minimum $`4 S c_0`$ of $`4 S \cdot \mathrm{charge}_{\min}(u) + (k_u - h_u) S`$ (with $`k_u`$ replaced by $`K_F + 1`$ on $`\mathrm{VSTAR}_h`$) yields $`c_0 C(T) \le \Phi(T)`$ for every finite tree $`T`$ with $`C(T) \ge 1`$. Since $`\mathbb{E}[\Phi(T_{n+1}) \mid T_n] = \Phi(T_n) + \frac{1}{n}(\mathcal{L}\Phi)(T_n)`$, induction from $`T_3`$ yields the product bounds of Theorem 1.

### Code Map (`round2/cert/`)
- `cert_b3f_format.h`: header struct, validation (`b3f_validate_header`), layout indices (`b3f_layout_init`), exact file size (`b3f_expected_bytes`).
- `cert_b3f_verify_gmp.c`:
  - `add_tail_gmp`: Section 4.1–4.2 exact pure-`mpz_t` Bernstein-3/4 tail verification.
  - `main` §0–§5: exact verification of `g[INF,INF]==0`, `f`, `g`, `hV`, `hC`, `P4`, `Ytop`, `Ybot`, `S3`.
  - `verify_comparison_bounds`: Section 5 exact `mpz_t` computation of $`S c_0`$, $`S C_0`$, $`S \Phi(T_3)`$.
  - `eval_Sphi_q` / `check_concrete_tree_q`: built-in exact `mpq_t` brute-force generator verification on boundary and random trees.
- `tree_identity_check.py`: independent Python `fractions.Fraction` reader and brute-force tree-drift verifier.
- `negative_controls.py`: automated rejection tests (impossible `--override-pq`, header/truncation corruption, non-zero `g[INF,INF]`, single-entry perturbation in all 12 tables).
