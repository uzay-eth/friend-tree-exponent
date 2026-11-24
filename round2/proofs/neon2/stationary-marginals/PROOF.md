# Mathematical Proof for Certified Stationary-Marginal Enclosures (`FT2_M5V1`)

**Pinned claim commit:** `6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d`.  
**Review and dependencies:** GPT-family full audits and historical replay; Claude supplied independent simulation consistency checks and checked the degree-two sum. See [Claim 10 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-10.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-10.md).

**Author:** Gemini 4 Neon (Neon-2).

This note proves that any pair of binary snapshots (`lower` with $`\sigma = -1`$, `upper` with $`\sigma = +1`$) accepted by the standalone zero-floating-point 128-bit verifier `verify_marginal.cpp` and the exact rational 2D tensor-Bernstein tail verifier `verify_marginal_tail.py` (combined by `cert_marginal.py` / `verify_all_marginals.py`) certifies an exact rational two-sided enclosure on the stationary marginal $`\mathcal M_\mu(\psi) = \pi(\psi)`$ of the continuous-time friend tree, conditional only on a certified rational bracket $`\mu \in [\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$ (supplied by Neon-2 Claim 6 `gold2`, commit `b435b4c`).

---

## 1. Unoriented Five-Vertex Coboundary Functional & Observables

Let $`\mathcal T`$ be the set of isomorphism classes of finite unrooted trees $`T`$ with core vertex set $`K(T) = \{u \in V(T) : D_T(u) \ge 2\}`$ of size $`C(T) = |K(T)| \ge 1`$, evolving under the physical continuous-time friend-tree generator $`\mathcal L`$ (`quintet-cert/PROOF.md` §1). For each $`u \in K(T)`$, write $`D_u = D_T(u) \ge 2`$, core degree $`c_u = |\{v \in K(T) : v \sim u\}|`$, and leaf count $`L_u = D_u - c_u \ge 0`$.

Fix one of the eight bounded observables $`\psi_{\mathrm{id}}(D, c) \in \{0, 1\}`$ (`marginal_format.hpp`):
1. `obs_id = 1`: $`\mathbf 1\{D = 2\}`$ ($`p_2`$),
2. `obs_id = 2`: $`\mathbf 1\{D = 3\}`$ ($`p_3`$),
3. `obs_id = 3`: $`\mathbf 1\{D = 4\}`$ ($`p_4`$),
4. `obs_id = 4`: $`\mathbf 1\{D = 5\}`$ ($`p_5`$),
5. `obs_id = 5`: $`\mathbf 1\{D = 6\}`$ ($`p_6`$),
6. `obs_id = 6`: $`\mathbf 1\{D = 2, c = 1\}`$ ($`\pi(J)`$),
7. `obs_id = 7`: $`\mathbf 1\{D = 2, c = 2\}`$ ($`\pi(2, 2)`$),
8. `obs_id = 8`: $`\mathbf 1\{D = c,\, D < D_F\}`$ ($`\pi(L_o = 0, D_o < D_F)`$; the far tail $`D \ge D_F`$ is bounded separately in §4).

Note that every $`\psi_{\mathrm{id}}`$ vanishes identically on $`D \ge D_F`$ (since $`D_F \ge 20 > 6`$). Write $`A_\psi(T) := \sum_{u \in K(T)} \psi(D_u, c_u)`$.

An `FT2_M5V1` binary snapshot has a 640-byte `MarginalHeader` wrapping a 512-byte `QuintetHeader` (with `base.mode == 1`, `alpha_lo_num == 0`, tick scale $`S = 2^{40}`$, center scale $`S_r = 2^{30}`$), specifying `obs_id` $`\in \{1,\dots,8\}`$, sign $`\sigma = \mathrm{obs\_sign} \in \{-1, +1\}`$, nominal target $`\theta = \theta_P / \theta_Q \in [0, 1]`$, certified $`\mu`$ endpoints $`0 < \mu_{\mathrm{lo}} \le \mu_{\mathrm{hi}} < 1`$, and four signed analytic tail parameters $`(f_{p0}, f_{p1}, f_{p2}, f_{p3}) \in [-16 S, 16 S]^4`$. It defines a **signed** unoriented 5-vertex (diameter-4) local potential $`\Phi : \mathcal T \to \mathbb R`$ by the exact same 8-shape motif expansion as `quintet-cert/PROOF.md` §2:
```math
\Phi(T) \;=\; \sum_{u \in K(T)} f(D_u, c_u) + \sum_{\{u, v\} \in E_K(T)} g(u, v) + \sum_{H_3} h_3 + \sum_{P_4} p_4 + \sum_{S_3} s_3 + \sum_{P_5} p_5 + \sum_{B_5} b_5 + \sum_{S_4} s_4,
```
where for $`D \ge D_F`$, writing $`q = L/D = (D - c)/D \in [0, 1]`$ and $`x = 1/D \in (0, 1/D_F]`$,
```math
S \cdot f(D, c) \;=\; w(q, x) \;:=\; f_{p0} + f_{p1} q + f_{p2} q^2 + f_{p3} q x.
```
Because `base.mode == 1`, all signed upper-bound transfers, 1D/2D/3D envelopes, Fenchel–Young degree-transfer knots, overflow/tail endpoint bounds, and per-core multiplicity aggregations are identical to the upper-bound mode of `quintet-cert/PROOF.md` §§2–4 (which uses no sign assumption on the weights $`f, g, h_3, p_4, s_3, p_5, b_5, s_4`$ themselves; only the direction $`\le`$ of the generator inequality matters).

---

## 2. Deterministic Coboundary Inequality Across $`[\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$ & Linear Boundedness

### Theorem 1 (Exact Deterministic Coboundary & Linear Boundedness)
Suppose a snapshot passes `verify_marginal.cpp` and `verify_marginal_tail.py` with verified endpoint error caps $`E(\mu_{\mathrm{lo}}), E(\mu_{\mathrm{hi}}) \ge 0`$ and linear bound $`M_{\mathrm{ticks}} \ge 0`$. Set $`\varepsilon_* := \max(E(\mu_{\mathrm{lo}}), E(\mu_{\mathrm{hi}})) / S`$ and $`M := M_{\mathrm{ticks}} / S`$. Then for **every** finite tree $`T \in \mathcal T`$ ($`C(T) \ge 1`$) and **every** $`\lambda \in [\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$:
```math
(\mathcal L \Phi - \lambda \Phi)(T) + \sigma\big(A_\psi(T) - \theta C(T)\big) \;\le\; \varepsilon_* C(T),
\qquad\text{(1)}
```
```math
|\Phi(T)| \;\le\; M\, C(T).
\qquad\text{(2)}
```

**Proof.**
1. **Finite rows and analytic tail at each endpoint $`\lambda \in \{\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}\}`$:**
   Decomposing $`(\mathcal L \Phi - \lambda \Phi)(T) + \sigma(A_\psi(T) - \theta C(T))`$ over the 8 motif families exactly as in `quintet-cert/PROOF.md` §3 places the 1-vertex forcing $`\sigma(\psi(D_u, c_u) - \theta)`$ exclusively on the singleton row $`R_f(u)`$ of each core vertex $`u \in K(T)`$.
   - For $`2 \le D < D_F`$ and $`0 \le c \le D`$, `verify_marginal.cpp` checks every finite singleton row in exact directed 128-bit sub-tick arithmetic ($`K_{\mathrm{SUB}} = 1024`$), including the exact 4-term seam $`w((D_F - c)/D_F, 1/D_F)`$ at $`D = D_F - 1`$ and the exact rational forcing $`\sigma(\psi(D, c)\theta_Q - \theta_P) S / \theta_Q`$.
   - For $`D \ge D_F`$, $`\psi(D, c) = 0`$, the gain jump satisfies $`(1+x)^2(w((q+x)/(1+x), x/(1+x)) - w(q,x))/x = \mathrm{gain}(q,x)`$, the promotion jump satisfies $`w(q-x,x) - w(q,x) = \mathrm{bud}(q,x)`$, the newborn edge satisfies $`S\cdot(f(J) + g((D, c+1), J)) \le C_{\mathrm{nb}} + d_g x`$ with $`C_{\mathrm{nb}} = S f(J) + S g(\infty, J)`$ and $`d_g = d_2(\max_{1 \le \kappa \le K_2+1} S g(T_2(\kappa), J) - S g(\infty, J))`$, and the gain rate $`r_u`$ lies in $`[L_u, L_u + c_u/2] = [q/x, (1+q)/(2x)]`$. By affinity in $`r_u`$ and positivity of $`(1+x)^2 \in [1, (1+1/D_F)^2]`$, `verify_marginal_tail.py` bounds the exact numerator polynomial on $`[0, 1] \times [0, 1/D_F]`$ across 16 rational Bernstein panels (plus an independent exact concrete-state oracle check at $`D \in \{D_F, D_F+1, 3D_F+17\}`$), yielding $`\mathrm{cap}_{f,\mathrm{ana}}(\lambda)`$.
   - Every row of $`g, H_3, P_4, S_3, P_5, B_5, S_4`$ and all 7 exact 128-bit cross-level envelopes (`env23` through `env45s4`) are verified verbatim as in `quintet-cert/PROOF.md` §3 at both $`\lambda = \mu_{\mathrm{lo}}`$ and $`\lambda = \mu_{\mathrm{hi}}`$ on the exact same SHA-256-authenticated in-memory buffer.
   - Aggregating per core vertex exactly as in `quintet-cert/PROOF.md` §4 (with $`C(T)`$ singletons, $`C(T)-1 < C(T)`$ core edges, and every positive higher-motif cap charged to its finite internal host(s)), setting $`\mathrm{cap}_f(\lambda) = \max(\mathrm{cap}_{f,\mathrm{in}}(\lambda), \mathrm{cap}_{f,\mathrm{tab}}(\lambda), \mathrm{cap}_{f,\mathrm{ana}}(\lambda))`$ and $`E(\lambda) = \max(0, \mathrm{cap}_f(\lambda) + \max(0, \mathrm{cap}_g(\lambda)) + \mathrm{higher\_per\_core}(\lambda))`$ yields
     $`\displaystyle (\mathcal L\Phi - \lambda\Phi)(T) + \sigma(A_\psi(T) - \theta C(T)) \;\le\; \frac{E(\lambda)}{S}\, C(T) \qquad (\lambda \in \{\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}\}).`$
2. **Affinity in $`\lambda`$:**
   For each fixed $`T \in \mathcal T`$, the left-hand side of (1) is affine in $`\lambda`$. Writing any $`\lambda \in [\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$ as $`(1-\alpha)\mu_{\mathrm{lo}} + \alpha\mu_{\mathrm{hi}}`$ with $`\alpha \in [0, 1]`$ gives (1) immediately with $`\varepsilon_* = \max(E(\mu_{\mathrm{lo}}), E(\mu_{\mathrm{hi}}))/S`$.
3. **Linear boundedness (2):**
   Rooting the core tree arbitrarily at any $`r \in K(T)`$ allocates each edge $`\{p_u, u\}`$ to $`u \ne r`$ and each higher motif to its finite internal host(s) (`quintet-cert/PROOF.md` §4). Vertices with $`D \ge D_F`$ host no higher motif (all higher-motif interior/arm cutoffs $`d_{3m}, d_{4i}, d_{S3}, d_{5i}, d_{5o}, d_{5c}, d_{S4} \le d_2 < D_F`$), so their allocated tick weight lies in $`[w_{\min} + g_{\mathrm{far\_min}}, w_{\max} + g_{\mathrm{far\_max}}]`$ (and in $`[w_{\min}, w_{\max}]`$ when $`C(T) = 1`$). Taking the extrema over all finite $`(D, c)`$ ($`D < D_F`$) and the analytic tail gives $`m_{\mathrm{lo}} C(T) \le S\Phi(T) \le m_{\mathrm{hi}} C(T)`$ for all $`T \in \mathcal T`$, hence $`|\Phi(T)| \le M C(T)`$ with $`M = \max(|m_{\mathrm{lo}}|, |m_{\mathrm{hi}}|)/S`$. $`\square`$

---

## 3. Stopped Dynkin & Exact Stationary Marginal Deduction

### Theorem 2 (Certified Stationary Marginal Enclosures)
Assume $`\mu \in [\mu_{\mathrm{lo}}, \mu_{\mathrm{hi}}]`$ (Neon-2 Claim 6 `gold2`, commit `b435b4c`) and let $`\mathcal M_\mu`$ be the unique stationary fringe probability measure of Neon-2 Claim 3 (`claim14-proofs/README.md`, commit `ea6a093`, which Theorem 1(2) identifies with the downward-fringe marginal $`\pi_{\mathrm{fringe}}`$ of the unique spatial eigenmeasure $`\pi`$ of Claim 4/10). For any snapshot satisfying Theorem 1:
1. **(Discounted Cesàro bound from $`T_3`$, using only positive mean growth $`e^{-\mu t}\mathbb E_{T_3}[C(T_t)] \to c_* \in (0,\infty)`$, without $`R_{\psi,\mu}`$):** Writing $`I(t) := \int_0^t e^{-\mu s}\mathbb E_{T_3}[C(T_s)]\,ds`$,
   $`\displaystyle \sigma \limsup_{t\to\infty}\frac{\int_0^t e^{-\mu s}\mathbb E_{T_3}[A_\psi(T_s)]\,ds}{I(t)} \;\le\; \sigma\theta + \varepsilon_* \quad (\sigma = +1), \qquad \liminf_{t\to\infty}\frac{\int_0^t e^{-\mu s}\mathbb E_{T_3}[A_\psi(T_s)]\,ds}{I(t)} \;\ge\; \theta - \varepsilon_* \quad (\sigma = -1).`$
2. **(Exact stationary marginal bound on $`\mathcal M_\mu(\psi) = \pi(\psi)`$):**
   $`\displaystyle \sigma\big(\mathcal M_\mu(\psi) - \theta\big) \;\le\; \varepsilon_*, \qquad\text{(3)}`$
   i.e. $`\mathcal M_\mu(\psi) \le \theta + \varepsilon_*`$ when $`\sigma = +1`$, and $`\mathcal M_\mu(\psi) \ge \theta - \varepsilon_*`$ when $`\sigma = -1`$.

**Proof.**
Start the continuous-time friend tree $`(T_t)_{t \ge 0}`$ from the 3-vertex tree $`T_3`$ ($`C(T_3) = 1`$). Its total vertex count $`N(t) = |T_t|`$ is a rate-1 Yule process with $`\mathbb E[N(t)] = 3 e^t < \infty`$ and $`\mathbb E[N(t)^2] < \infty`$.

By Neon-2 Claim 3 Corollary 2.1(1) (`claim14-proofs/README.md` eq. (10), contour-free, audited `VERDICT: CORRECT`), since $`\psi(D, c)`$ is a bounded degree score, the explicit coboundary $`R_{\psi, \mu}(T)`$ satisfies **pointwise on every $`T \in \mathcal T`$**:
```math
(\mathcal L R_{\psi, \mu} - \mu R_{\psi, \mu})(T) \;=\; \mathcal M_\mu(\psi)\, C(T) - A_\psi(T), \qquad |R_{\psi, \mu}(T)| \;\le\; K_\psi C(T).
```
Set $`G(T) := \Phi(T) + \sigma R_{\psi, \mu}(T)`$. Then $`|G(T)| \le (M + K_\psi) C(T)`$, and adding $`\sigma`$ times the coboundary identity to (1) at $`\lambda = \mu`$ cancels $`A_\psi(T)`$ **identically on every finite tree $`T \in \mathcal T`$**:
```math
(\mathcal L G - \mu G)(T) + \sigma\big(\mathcal M_\mu(\psi) - \theta\big) C(T) \;\le\; \varepsilon_* C(T).
\qquad\text{(4)}
```
For each integer $`K \ge 4`$, let $`\tau_K := \inf\{t \ge 0 : N(t) \ge K\} \uparrow \infty`$ a.s. On $`[0, t \wedge \tau_K]`$ the chain visits finitely many isomorphism classes, so Dynkin's formula applied to $`e^{-\mu s} G(T_s)`$ and (4) give
```math
\mathbb E_{T_3}\big[e^{-\mu(t \wedge \tau_K)} G(T_{t \wedge \tau_K})\big] - G(T_3) + \sigma\big(\mathcal M_\mu(\psi) - \theta\big)\mathbb E_{T_3}\!\left[\int_0^{t \wedge \tau_K} e^{-\mu s} C(T_s)\,ds\right] \;\le\; \varepsilon_* \mathbb E_{T_3}\!\left[\int_0^{t \wedge \tau_K} e^{-\mu s} C(T_s)\,ds\right].
```
Since $`|e^{-\mu(t \wedge \tau_K)} G(T_{t \wedge \tau_K})| \le (M + K_\psi) N(t) \in L^1`$ and $`0 \le e^{-\mu s} C(T_s) \le N(t) \in L^1`$ for $`s \in [0, t]`$, letting $`K \to \infty`$ by dominated convergence yields
```math
\sigma\big(\mathcal M_\mu(\psi) - \theta\big) I(t) \;\le\; \varepsilon_* I(t) + G(T_3) + (M + K_\psi) e^{-\mu t}\mathbb E_{T_3}[C(T_t)].
```
By Neon-2 Claim 3 Proposition 0 + Theorem 2 (`ea6a093`, contour-free), $`e^{-\mu t}\mathbb E_{T_3}[C(T_t)] \to c_* := \Psi^0(J)\Phi^{\mathrm{eig}}(T_3) \in (0, \infty)`$ as $`t \to \infty`$, so $`I(t) \sim c_* t \to \infty`$ and $`(G(T_3) + (M + K_\psi) e^{-\mu t}\mathbb E_{T_3}[C(T_t)]) / I(t) = O(1/t) \to 0`$. Dividing by $`I(t)`$ and letting $`t \to \infty`$ proves (3). Applying the identical stopped-Dynkin step to $`\Phi`$ alone (without $`R_{\psi, \mu}`$) proves part (1). $`\square`$

---

## 4. Exact Analytic Tail Bound for Bare Cores (`obs_id = 8`, $`L_o = 0, D_o \ge D_F`$)

For `obs_id = 8`, the verified finite-table observable is $`\psi_8(D, c) = \mathbf 1\{D = c,\, D < D_F\}`$, which gives $`\mathcal M_\mu(L_o = 0, D_o < D_F) \le \mathcal M_\mu(L_o = 0)`$ (so the certified lower bound on $`\psi_8`$ is automatically a lower bound on all bare cores $`L_o = 0`$). For the upper bound, we bound the excluded tail $`\mathcal M_\mu(L_o = 0, D_o \ge D_F)`$ analytically:

### Theorem 3 (Four-Line Bare-Core Tail Bound)
For every integer $`D_F \ge 6`$ and any $`\mu \ge \mu_{\mathrm{lo}} > 0`$,
```math
0 \;\le\; \mathcal M_\mu(L_o = 0,\, D_o \ge D_F) \;\le\; \frac{40}{\mu_{\mathrm{lo}}}\, 2^{-D_F}.
\qquad\text{(5)}
```
(For our `high` preset $`D_F = 1600`$, this tail is $`< 10^{-479}`$.)

**Proof.**
Define $`H(T) := \sum_{u \in K(T)} h(D_u, c_u)`$ with
```math
h(D, c) \;:=\; \begin{cases}
0, & 2 \le D \le 5, \\ 2^{-(D_F - D)} 4^{-(D - c)}, & 6 \le D \le D_F, \\ 4^{-(D - c)}, & D \ge D_F.
\end{cases}
```
Clearly $`0 \le h(D, c) \le 1`$, so $`0 \le H(T) \le C(T)`$ on all $`T \in \mathcal T`$, and $`h(2, 1) = 0`$ (newborn core vertices contribute $`0`$). Note that the two formulas agree at $`D = D_F`$. Check the per-vertex drift $`(\mathcal L h - \mu h)(D, c)`$ with $`L = D - c \ge 0`$:
1. **$`D \ge D_F, L = 0`$:** Gain moves $`(D, D) \to (D+1, D)`$ with jump $`4^{-1} - 1 = -3/4 < 0`$ (at rate $`r \ge 0`$), promotion rate is $`L/D = 0`$, so $`(\mathcal L h - \mu h)(D, D) \le -\mu \cdot 1 = -\mu`$.
2. **$`D \ge D_F, L \ge 1`$:** Gain jump is $`-\frac{3}{4} 4^{-L} < 0`$ at rate $`r \ge L`$; promotion jump is $`+3 \cdot 4^{-L}`$ at rate $`L/D \le L/6`$. Hence $`(\mathcal L h - \mu h)(D, c) \le (-\frac{3}{4} + \frac{3}{6}) L 4^{-L} - \mu 4^{-L} < 0`$.
3. **$`6 \le D \le D_F - 1, L = 0`$:** Gain moves $`(D, D) \to (D+1, D)`$, multiplying $`h`$ by $`2 \cdot 4^{-1} = 1/2`$, so gain jump is $`-\frac{1}{2} h(D, D) < 0`$ and promotion rate is $`0`$, giving $`(\mathcal L h - \mu h)(D, D) \le -\mu h(D, D) < 0`$.
4. **$`6 \le D \le D_F - 1, L \ge 1`$:** Gain jump is $`(\frac{2}{4} - 1) h(D, c) = -\frac{1}{2} h(D, c) < 0`$ at rate $`r \ge L`$; promotion jump is $`(4 - 1) h(D, c) = 3 h(D, c) > 0`$ at rate $`L/D \le L/6`$. Hence $`(\mathcal L h - \mu h)(D, c) \le (-\frac{1}{2} + \frac{3}{6}) L h(D, c) - \mu h(D, c) = -\mu h(D, c) \le 0`$.
5. **$`D = 5, 0 \le c \le 5`$:** $`h(5, c) = h(5, c+1) = 0`$. Gain moves $`(5, c) \to (6, c)`$ at rate $`r \le L + c/2 = 5 - c/2`$ with $`h(6, c) = 2^{-(D_F - 6)} 4^{-(6 - c)}`$. Checking all $`c \in \{0,\dots,5\}`$ gives $`(5 - c/2) 4^{-(6 - c)} \in \{\frac{5}{4096}, \frac{9}{2048}, \frac{1}{64}, \frac{7}{128}, \frac{3}{16}, \frac{5}{8}\}`$, maximized at $`c = 5`$ by $`\frac{5}{8}`$. Thus $`(\mathcal L h - \mu h)(5, c) \le \frac{5}{8} 2^{-(D_F - 6)} = 40 \cdot 2^{-D_F}`$.
6. **$`2 \le D \le 4`$:** $`h = 0`$ on $`(D, c)`$, $`(D+1, c)`$, and $`(D, c+1)`$, so drift is $`0`$.

Summing over $`u \in K(T)`$ gives the exact pointwise inequality on every $`T \in \mathcal T`$:
```math
(\mathcal L H - \mu H)(T) + \mu \sum_{u \in K(T)} \mathbf 1\{D_u \ge D_F, L_u = 0\} \;\le\; 40 \cdot 2^{-D_F}\, C(T).
```
Applying the exact same stopped-Dynkin + Cesàro proof as Theorem 2 (combining with $`R_{\psi_{\mathrm{tail}}, \mu}`$ for $`\psi_{\mathrm{tail}}(D, c) = \mathbf 1\{D \ge D_F, L = 0\}`$) and dividing by $`\mu \ge \mu_{\mathrm{lo}}`$ yields (5). $`\square`$

---

## 5. Exact Stationary Balance Identities for Degree-2 Core Vertices

### Theorem 4 (Degree-2 Decomposition & Conditional Neighbor-Reciprocal Expectations)
Let $`f_J(T) := \sum_{u \in K(T)} \mathbf 1\{D_u = 2, c_u = 1\}`$, $`f_{22}(T) := \sum_{u \in K(T)} \mathbf 1\{D_u = 2, c_u = 2\}`$, $`Z_J(T) := \sum_{u \in K(T)} \mathbf 1\{D_u = 2, c_u = 1\}\frac{1}{D_{p_u}}`$ (where $`p_u`$ is the unique core neighbor of $`u`$), and $`Z_{22}(T) := \sum_{u \in K(T)} \mathbf 1\{D_u = 2, c_u = 2\}\big(\frac{1}{D_{v_1(u)}} + \frac{1}{D_{v_2(u)}}\big)`$ (where $`v_1(u), v_2(u)`$ are the two core neighbors of $`u`$). Then:
1. $`\mathcal M_\mu(2, 0) = 0`$, hence $`p_2 = \mathcal M_\mu(D_o = 2) = \mathcal M_\mu(J) + \mathcal M_\mu(2, 2)`$.
2. Under the discounted Cesàro limit from $`T_3`$ (and hence under the spatial eigenmeasure $`\pi`$ of Claim 4/10, where $`\pi(J) = \mathcal M_\mu(J) > 0`$ and $`\pi(2, 2) = \mathcal M_\mu(2, 2) > 0`$),
   $`\displaystyle \lim_{t\to\infty}\frac{\int_0^t e^{-\mu s}\mathbb E_{T_3}[Z_J(T_s)]\,ds}{\int_0^t e^{-\mu s}\mathbb E_{T_3}[f_J(T_s)]\,ds} \;=\; \frac{\mu}{\mathcal M_\mu(J)} - \mu - \frac{3}{2}, \qquad\text{(6)}`$
   $`\displaystyle \lim_{t\to\infty}\frac{\int_0^t e^{-\mu s}\mathbb E_{T_3}[Z_{22}(T_s)]\,ds}{\int_0^t e^{-\mu s}\mathbb E_{T_3}[f_{22}(T_s)]\,ds} \;=\; \frac{\mathcal M_\mu(J)}{2\,\mathcal M_\mu(2, 2)} - \mu. \qquad\text{(7)}`$

**Proof.**
1. On $`\mathcal S`$, every planted tree retains its inert parent slot counted in $`c_o`$, so $`c_o \ge 1`$ everywhere on $`\mathcal S`$ and hence $`\mathcal M_\mu(2, 0) = 0`$ and $`p_2 = \mathcal M_\mu(D_o = 2) = \mathcal M_\mu(J) + \mathcal M_\mu(2, 2)`$ directly. Moreover, on the free tree $`(T_s)_{s \ge 0}`$ started from $`T_3`$, a core vertex has $`c_u = 0`$ if and only if $`C(T_s) = 1`$, which holds only before the first leaf promotion of the single initial core; while $`C(T_s) = 1`$, $`c_o = 0`$ so $`L_o = D_o`$ and the promotion rate is identically $`L_o / D_o = 1`$, giving $`\mathbb P_{T_3}(C(T_s) = 1) = e^{-s}`$ and $`\int_0^\infty e^{-\mu s}\mathbb P_{T_3}(C(T_s) = 1)\,ds = \frac{1}{1+\mu} < \infty`$.
2. For any $`T \in \mathcal T`$ with $`C(T) \ge 2`$, every promotion creates one new $`(2, 1)`$ core (total birth rate $`(\mathcal L C)(T)`$), and no existing core can enter $`(2, 1)`$ (since $`D`$ never decreases and $`c \ge 1`$ for all cores when $`C \ge 2`$). Each existing $`u`$ in state $`(2, 1)`$ leaves $`(2, 1)`$ via own-leaf gain (rate $`1`$), neighbor-induced gain (rate $`1/D_{p_u}`$), or own-leaf promotion to $`(2, 2)`$ (rate $`1/2`$). Hence
   $`\displaystyle (\mathcal L f_J)(T) \;=\; (\mathcal L C)(T) - \frac{3}{2} f_J(T) - Z_J(T) + \mathbf 1\{C(T) = 1, D_o = 2\},`$
   where at $`C(T) = 1, D_o = 2`$, promoting one leaf turns the root into $`(2, 1)`$ and creates a second $`(2, 1)`$ vertex (change $`+2 = (\mathcal L C)(T) + 1`$). Similarly, on every $`T \in \mathcal T`$ ($`C(T) \ge 1`$), a vertex enters $`(2, 2)`$ solely when a $`(2, 1)`$ vertex promotes its leaf (rate $`1/2`$) and leaves $`(2, 2)`$ solely via neighbor-induced gain (rate $`1/D_{v_1(u)} + 1/D_{v_2(u)}`$), giving
   $`\displaystyle (\mathcal L f_{22})(T) \;=\; \frac{1}{2} f_J(T) - Z_{22}(T).`$
   Applying Yule-dominated stopped Dynkin to $`e^{-\mu s} f_J(T_s)`$, $`e^{-\mu s} f_{22}(T_s)`$, and $`e^{-\mu s} C(T_s)`$ on $`[0, t]`$ (each bounded by $`C(T_s)`$), subtracting $`\mu`$ times each integral, using $`\int_0^\infty e^{-\mu s}\mathbb P_{T_3}(C(T_s)=1)\,ds < \infty`$, and dividing by $`I(t) \sim c_* t \to \infty`$ (with Theorem 2(1) and Claim 3 Corollary 2.1 giving $`\int_0^t e^{-\mu s}\mathbb E[f_J(T_s)]\,ds / I(t) \to \mathcal M_\mu(J)`$ and $`\int_0^t e^{-\mu s}\mathbb E[f_{22}(T_s)]\,ds / I(t) \to \mathcal M_\mu(2, 2)`$) yields (6) and (7). Under Claim 4/10, the left-hand sides of (6) and (7) equal $`\mathbb E_\pi[1/D_{\mathrm{parent}} \mid J]`$ and $`\mathbb E_\pi[1/D_{v_1} + 1/D_{v_2} \mid (2, 2)]`$, respectively. $`\square`$
