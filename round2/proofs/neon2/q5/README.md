# Track 5B: Simultaneous Cutoff Limits and Total-Variation Regimes for Addario-Berry et al. Open Question 5

**Pinned claim commit:** `6950241589978a7064ade0f6a2ee97cbbdcc6a34`.  
**Review and dependencies:** GPT-family proof audits; Claude independently derived the negative total-variation implication. See [Claim 5 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-5.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-5.md).

**Author:** Neon-2 lead (Gemini 4 Neon).
**Cross-family cold audit:** GPT-6.2 Sol Ultrafast (`verify_q5_exact.py`, `sim_q5.cpp`); this revision incorporates every repair from the first audit (removing the unproved fixed-shape intensity at scale $`K`$; fixing the dirty-core upward boundary charge, Poisson tail constant, $`Y`$-localization, exponential martingale tail from `almost-sure.md` (8), localized stopping time and $`2+c_x`$ inter-mesh jump bound; finite-horizon vanishing-tolerance Yule concentration; promotion-age split and integer-bin Borel–Cantelli in Step 3.2; simultaneous coupling with single stopped send count and no-descendant-hub property; and exact $`(t_-, t_+]`$ upper-CDF accounting).

**Context & explicit Astra dependencies:** Astra's `round2/open/forest-age/README.md` answered Addario-Berry et al. §7(5) in the **sequential** limit ($`n \to \infty`$ at fixed $`N`$, then $`N \to \infty`$), noting explicitly: *"These results retain the sequential order of limits. A threshold growing with the observed finite tree size requires an additional uniform estimate."* Here we resolve the **simultaneous** regime $`n \to \infty, N = N(n) \to \infty`$, proving:
1. **Full simultaneous range $`N(n) \to \infty, N(n)/n \to 0`$:** exact simultaneous component count and simultaneous rescaled log-volume convergence to $`\mathrm{Uniform}[0, 1]`$ almost surely;
2. **Separated simultaneous total-variation regimes:** a.s. total-variation convergence $`\|\nu_{n, N(n)} - \mu_{N(n)}\|_{\mathrm{TV}} \to 0`$ for logarithmic thresholds $`N(n) \le c_0 \log n`$, versus a.s. non-vanishing TV distance (reaching $`1`$ at $`N(n) = n^{1-o(1)}`$) for polynomial thresholds $`N(n) = n^\alpha`$ with $`\alpha > \mu/(1+\mu)`$.

**Explicit Astra inputs cited:**
- `round2/open/forest-age/README.md`: sequential stationary measure $`M_N(f) = \mathbb E_\pi[\phi_{N,f}]`$, $`d_N = M_N(1) \sim A N^{-\mu}\log N`$, and $`\|\nu_N - \mu_N\|_{\mathrm{TV}} \to 0`$ as $`N \to \infty`$ (where components include isolated surviving original cores of size $`1`$).
- `round2/fringe/almost-sure.md` (Claim 6) & `round2/fringe/finite-kernel.md` (Claim 5) & `round2/fringe/renewal.md`: 19-label weight $`W(S)`$, exponential spectral gap $`\sigma < \mu`$ uniform for signed $`|f| \le W`$ (`almost-sure.md` eq. (4)), dirty-core coupling (`almost-sure.md` §3 eq. (6)), variance bound `almost-sure.md` (8)–(9), and $`C_n/n^\mu \to X \in (0,\infty)`$, $`n e^{-\tau_n} \to Y \in (0,\infty)`$ ($`Y \sim \mathrm{Gamma}(3,1)`$ from $`T_3`$) a.s.
- `round2/fringe/hub-spectrum.md` (Claim 11): common-event uniform intermediate degree law $`N_{\ge k}(n) \sim A X (n/k)^\mu`$ for all $`k \to \infty, k/n \to 0`$ (§4); $`\mathbb E[\sum_{x \in T_t} e^{c_x(t)}] \le K e^{e t}`$ and $`\max_{x \in T_t} c_x(t) = O(t)`$ a.s. (§3); tagged-leaf Yule/gamma-minus-exponential degree comparisons (§2).
- `round2/fringe/hub-children.md` §3: Poisson-marked inert-virgin coupling of a core's promoted children (extended below to simultaneous $`\mathcal F_{s_0}`$-measurable hubs).

---

## 1. Definitions and Main Theorems

For $`n \ge 3`$ and integer $`N \ge 2`$, delete all vertices of the friend tree $`T_n`$ whose original degree exceeds $`N`$. Let $`\mathcal C_{n,N}`$ be the multiset of remaining connected components containing at least one original nonleaf (core) — including isolated surviving original cores of size $`1`$ — each rooted at its highest surviving core (closest to the founder, forgetting the removed parent edge). Set
```math
M_{n, N}(f) \;:=\; \sum_{F \in \mathcal C_{n, N}} f(F), \qquad M_{n, N} \;:=\; M_{n, N}(1), \qquad \nu_{n, N}(f) \;:=\; \frac{M_{n, N}(f)}{M_{n, N}} \quad (M_{n, N} > 0),
```
and let $`\mu_N := \frac{1}{\log N} \int_0^{\log N} P_a^0(J, \cdot)\,da`$, where $`P_a^0(J, \cdot)`$ is the law of the ordinary rooted tree of an inert planted virgin $`J`$ at elapsed age $`a`$ (omitting its inert parent slot). Write $`V(F) = |F| \in \mathbb Z_{\ge 1}`$ for ordinary vertex count and $`A = \lim_{k\to\infty} k^\mu \pi(D_o \ge k) > 0`$.

### Theorem 1 (Simultaneous Count and Log-Volume Law for ALL $`N(n) \to \infty, N(n)/n \to 0`$)
Let $`N(n) \in \mathbb Z_{\ge 2}`$ be **any** deterministic sequence satisfying $`N(n) \to \infty`$ and $`N(n)/n \to 0`$ as $`n \to \infty`$. Then almost surely:
1. **(Simultaneous component count):**
   $`\displaystyle \frac{M_{n, N(n)}}{A X \,(n / N(n))^\mu \log N(n)} \;\longrightarrow\; 1. \qquad\text{(1)}`$
2. **(Simultaneous rescaled log-volume law):** For the ordinary vertex count $`V_{n, N}`$ of a component sampled uniformly from $`\mathcal C_{n, N(n)}`$ given $`T_n`$,
   $`\displaystyle \sup_{x \in [0, 1]} \Big|\nu_{n, N(n)}\!\big(V \le N(n)^x\big) - x\Big| \;\longrightarrow\; 0, \qquad\text{(2)}`$
   i.e., $`\mathcal L(\log V_{n, N(n)} / \log N(n) \mid T_n) \Rightarrow \mathrm{Uniform}[0, 1]`$ almost surely.

### Theorem 2 (Separated Simultaneous Total-Variation Regimes)
1. **(Convergence regime — logarithmic thresholds):** There exists an explicit deterministic constant $`c_0 > 0`$ (depending only on $`\mu`$ and the spectral gap $`\sigma < \mu`$ of Claim 6) such that for **every** deterministic sequence $`N(n) \to \infty`$ with $`N(n) \le c_0 \log n`$,
   $`\displaystyle \big\|\nu_{n, N(n)} - \mu_{N(n)}\big\|_{\mathrm{TV}} \;\longrightarrow\; 0 \qquad \text{almost surely.} \qquad\text{(3)}`$
2. **(Breakdown regime — polynomial thresholds):** For **every** deterministic sequence $`N(n) \to \infty`$ with $`N(n)/n \to 0`$ and $`\alpha := \liminf_{n\to\infty} \frac{\log N(n)}{\log n} > \frac{\mu}{1 + \mu}`$ ($`\approx 0.3613`$),
   $`\displaystyle \liminf_{n\to\infty} \big\|\nu_{n, N(n)} - \mu_{N(n)}\big\|_{\mathrm{TV}} \;\ge\; 1 - \mu\!\left(\frac{1}{\alpha} - 1\right) \;>\; 0 \qquad \text{almost surely.} \qquad\text{(4)}`$
   In particular, whenever $`N(n) = n^{1 - o(1)}`$ with $`N(n)/n \to 0`$ (for example $`N(n) = \lfloor n / \log n \rfloor`$), $`\|\nu_{n, N(n)} - \mu_{N(n)}\|_{\mathrm{TV}} \to 1`$ almost surely — even though (1)–(2) hold.

*(Honest scope: Theorem 2 establishes separated sufficient regimes for TV convergence and TV breakdown; it does not locate a sharp single threshold between $`c_0 \log n`$ and $`n^{\mu/(1+\mu)}`$, nor does it prove volume-only TV convergence below $`\mu/(1+\mu)`$ or fixed-shape intensity at the $`O((n/N)^\mu)`$ scale.)*

---

## 2. Engine B: Uniform Concentration for Small-to-Polynomial Cutoffs $`2 \le N \le n^{\gamma_0}`$

Work in continuous time $`(T_t)_{t \ge 0}`$ started from $`T_3`$, with size-$`n`$ epochs $`\tau_n`$ ($`n e^{-\tau_n} \to Y \sim \mathrm{Gamma}(3,1)`$ a.s., $`X = d_2 Z^*_\infty / Y^\mu`$). Write $`M_{t, N}(f)`$ for the component sum on $`T_t`$. Every component in $`\mathcal C(T_t, N)`$ has a unique top core, which is either the founder $`o`$ or a core child $`w`$ of some deleted parent $`v \in K(T_t)`$ with $`D_v(t) > N`$. On the 19-marked planted state space of Claim 5 (`finite-kernel.md`), define for any $`N \ge 2`$ and $`|f| \le 1`$ the downward-fringe score
```math
\phi_{N, f}(S) \;:=\; \mathbf 1_{\{D_o(S) > N\}} \sum_{w \text{ core child of } o \text{ in } S} \mathbf 1_{\{D_w(S) \le N\}} f\big(F_N(w)\big),
```
using actual degrees (including the top parent slot of $`S`$) and erasing birth marks and the exterior slot from the output component $`F_N(w)`$. Separating the founder $`o`$ (both as a possible component top when $`D_o(t) \le N`$ and as a possible deleted parent of up to $`c_o(t)`$ surviving core children when $`D_o(t) > N`$) gives
```math
\Big|M_{t, N}(f) - \sum_{v \in K(T_t) \setminus \{o\}} \phi_{N, f}(S_v(t))\Big| \;\le\; 1 + c_o(t).
\qquad\text{(5)}
```

### Lemma 2.1 (Uniform $`W`$-Domination and Spectral-Gap Mean Bound)
On the 19-label planted space (`finite-kernel.md` §1, `renewal.md` §3), for every $`N \ge 2`$ and every signed $`|f| \le 1`$, $`|\phi_{N, f}(S)| \le \frac{26}{7} W(S)`$. Consequently, there exists $`K_0 < \infty`$ independent of $`N, f, u, h`$ such that for all $`u \ge 0, h \ge 0`$,
```math
\Big|\mathbb E\big[M_{u+h, N}(f) \mid \mathcal F_u\big] - e^{\mu h} \Psi(T_u)\, d_2 M_N(f)\Big| \;\le\; K_0 \big(C(T_u) + c_o(u) + h + 1\big) e^{\sigma h}.
\qquad\text{(6)}
```

**Proof.** In `finite-kernel.md` §1 and `renewal.md` §3, $`W(S) = (a_S + d)^{-1}\sum_{x \in K(S)} w_x h(c_x)`$ with $`d = 3/20`$, $`a_S \in [1/20, 1/2]`$, $`h(c) \ge 7/10`$, and every immediate core child $`y`$ of the root carries capped birth mark $`a_y \ge 1/20`$, giving $`w_y = a_y / (a_y + d) \ge (1/20)/(4/20) = 1/4`$. Since the root of $`S`$ has $`c_o(S) - 1`$ real core children and all terms are nonnegative,
```math
W(S) \;\ge\; \frac{1}{1/2 + 3/20} \cdot \frac{7}{10} \cdot \frac{c_o(S) - 1}{4} \;=\; \frac{7}{26}\big(c_o(S) - 1\big) \;\ge\; \frac{7}{26}|\phi_{N, f}(S)|.
```
Applying `almost-sure.md` eq. (4) (which explicitly covers signed scores $`|g| \le W`$ with uniform constant) to $`\frac{7}{26}\phi_{N, f}`$ gives the mean expansion for $`\sum_{v \ne o} \phi_{N, f}(S_v(u+h))`$ with uniform error $`O(C(T_u) e^{\sigma h})`$ and leading constant $`d(\phi_{N, f})`$. To identify $`d(\phi_{N, f}) = d_2 M_N(f)`$: at fixed $`N`$, Astra's `open/forest-age/README.md` gives $`M_{t, N}(f)/C(T_t) \to M_N(f)`$ a.s. (for signed $`|f| \le 1`$ by sequential TV) with $`|M_{t, N}(f)/C(T_t)| \le 1`$. By $`C(T_t) \le C_3(1 + \Psi(T_t))`$ and the uniform $`L^2`$ bound of `almost-sure.md` (9), $`e^{-\mu t} C(T_t)`$ is uniformly integrable with mean limit $`d_2 \Psi(T_3) > 0`$, so $`e^{-\mu t}\mathbb E[M_{t, N}(f)] \to d_2 \Psi(T_3) M_N(f)`$. Comparing with `almost-sure.md` eq. (4) (since the founder error in (5) has expectation $`\le 1 + c_o(0) + t = O(t)`$) identifies $`d(\phi_{N, f}) = d_2 M_N(f)`$. The conditional founder bound in (5) is $`\mathbb E[1 + c_o(u+h) \mid \mathcal F_u] \le 1 + c_o(u) + h`$, absorbed into (6) since $`\sigma > 0`$. $`\square`$

### Lemma 2.2 (Uniform Dirty-Core One-Leaf Gradient & Unconditional Conditional-Variance Bound)
There exists a constant $`K_1 < \infty`$ independent of $`N \ge 2`$ and $`|f| \le 1`$ such that for all $`0 \le u \le v`$ with $`h = v - u`$,
```math
\mathbb E\!\left[\big(M_{v, N}(f) - \mathbb E[M_{v, N}(f) \mid \mathcal F_u]\big)^2\right] \;\le\; K_1 (v^2 + 1) e^{u + 4h}.
\qquad\text{(7)}
```

**Proof.** Fix any finite tree $`T`$ with at least one core, $`M(T) := \max_{x \in K(T)} c_x(T) \ge 0`$ (so $`M(T) + 1 \ge 1`$), and any one-leaf extension $`T^+`$ (noting $`|T^+| = |T| + 1`$ and $`M(T^+) \le M(T) + 1`$). Couple $`(T_r, T^+_r)_{r \ge 0}`$ by the dirty-core coupling of `almost-sure.md` §3 (eq. (6)), so the dirty-core sets $`D_{\text{dirty}}^{(1)}(r) \subset K(T_r)`$ and $`D_{\text{dirty}}^{(2)}(r) \subset K(T^+_r)`$ satisfy $`\mathbb E[|D_{\text{dirty}}^{(1)}(r)| + |D_{\text{dirty}}^{(2)}(r)|] \le 3 e^{2r}`$. Count paired actual components before erasing vertex identities: any component $`F`$ in one copy that contains **no** dirty core, has **no** dirty parent above its top, and has **no** dirty child directly below a surviving core of $`F`$ agrees identically with its paired component in the other copy (all vertices in $`F`$ are clean with matching leaf stocks and matching core neighborhoods, and all edges exiting $`F`$ go to paired deleted clean vertices). Thus every unpaired/changed component in copy $`i \in \{1, 2\}`$ either:
1. contains a dirty core (at most $`|D_{\text{dirty}}^{(i)}(r)|`$ by disjointness);
2. has a dirty parent directly above its top (at most $`\sum_{x \in D_{\text{dirty}}^{(i)}(r)} c_x^{(i)}(r)`$); or
3. has a dirty core directly below it across a boundary edge — charged upward to that dirty core's unique parent edge (at most $`|D_{\text{dirty}}^{(i)}(r)|`$).
Hence for every $`N \ge 2`$ and $`|f| \le 1`$,
```math
|M_{r, N}(f)(T^+_r) - M_{r, N}(f)(T_r)| \;\le\; 2\sum_{i=1}^2 \sum_{x \in D_{\text{dirty}}^{(i)}(r)} \big(2 + c_x^{(i)}(r)\big).
```
In either copy $`i \in \{1, 2\}`$ (denote the initial tree by $`\tilde T \in \{T, T^+\}`$), every physical core $`x \in K(\tilde T_r)`$ — whether an initial core, an initial leaf later promoted (entering at core degree $`1`$), or a post-$`0`$ birth — satisfies $`c_x(r) \le M(T) + 1 + N_x(r)`$, where $`N_x(r)`$ counts the rings of $`x`$'s own rate-1 clock during $`[0, r]`$ after $`x`$'s birth ($`N_x(r) \sim \mathrm{Pois}(r)`$ for initial vertices, and conditionally $`\mathrm{Pois}(r - \tau_x) \preceq \mathrm{Pois}(r)`$ for vertices born at $`\tau_x \in (0, r]`$). For any real $`K \ge 1`$,
```math
\sum_{x \in D_{\text{dirty}}^{(i)}(r)} c_x^{(i)}(r) \;\le\; \big(M(T) + 1 + K\big)|D_{\text{dirty}}^{(i)}(r)| + \sum_{x \in \tilde T_r} N_x(r) \mathbf 1_{\{N_x(r) > K\}}.
```
Using the strong Markov property at births (separating initial vertices from post-$`0`$ births, total expected living population $`|\tilde T| e^r`$) and the exact Poisson tail bound for $`Z \sim \mathrm{Pois}(r)`$:
```math
\mathbb E[Z \mathbf 1_{\{Z > K\}}] = r\,\mathbb P(Z \ge \lfloor K \rfloor) \le r\,e^{-\lfloor K \rfloor} e^{(e-1)r} \le C_P e^{2r - K}, \qquad C_P := e \sup_{r \ge 0} r e^{-(3 - e)r} < \infty.
```
Setting $`K := \log|T| + 3r + 1`$ gives $`\mathbb E_{\tilde T}[\sum_{x \in \tilde T_r} N_x(r)\mathbf 1_{\{N_x(r) > K\}}] \le 2 C_P`$, yielding the uniform one-leaf gradient
```math
\big|\mathbb E_{T^+}[M_{r, N}(f)] - \mathbb E_T[M_{r, N}(f)]\big| \;\le\; C_2\big(M(T) + \log|T| + r + 1\big) e^{2r}.
```
Applying backward martingale jump isometry to $`\mathbb E[M_{v, N}(f) \mid \mathcal F_s]`$ on $`s \in [u, v]`$ (justified by finite-size stopping and Yule moments) and taking expectations from $`T_3`$:
```math
\mathbb E\!\left[\big(M_{v, N}(f) - \mathbb E[M_{v, N}(f) \mid \mathcal F_u]\big)^2\right]
\;\le\; C_2^2 \int_u^v e^{4(v - s)} \,\mathbb E\!\left[|T_s|\big(M(T_s) + \log|T_s| + h + 1\big)^2\right] ds.
```
By Cauchy–Schwarz, $`\mathbb E[|T_s| M(T_s)^2] \le (\mathbb E[|T_s|^2])^{1/2}(\mathbb E[M(T_s)^4])^{1/2}`$. From the Yule(3) size process, $`\mathbb E[|T_s|^2] = O(e^{2s})`$ and $`\mathbb P(\log|T_s| > s + z) \le 3 e^{-z}`$, so $`\mathbb E[(\log|T_s|)^4] = O(1 + s^4)`$. From `hub-spectrum.md` §3 ($`\mathbb E[\sum_x e^{c_x(s)}] \le K e^{e s}`$), integrating $`\int_0^\infty 4 y^3 \min(1, K e^{es - y})\,dy`$ split at $`e s + \log\max(1, K) + 1`$ yields $`\mathbb E[M(T_s)^4] = O(1 + s^4)`$. Thus the inner expectation is $`O((v^2 + 1) e^s)`$, and $`\int_u^v e^{4(v-s)+s} ds \le \frac{1}{3} e^{u + 4h}`$ proves (7). $`\square`$

### Proposition 2.3 (Uniform Concentration on Continuous-Time Cutoffs $`2 \le N \le B e^{\gamma_0 (j+1)}`$)
Let $`\eta := \mu - \sigma > 0`$, $`a := \frac{2\mu - 1}{10} > 0`$, and $`\gamma_0 := \frac{\min(\eta a,\, 2\mu - 1)}{20} > 0`$. Fix any deterministic integer $`B \ge 1`$ and, for each integer $`j \ge 1`$ and $`2 \le N \le B e^{\gamma_0(j+1)}`$, let $`\mathcal F_{j, N}`$ be any deterministic finite family of test functions $`|f| \le 1`$ with $`|\mathcal F_{j, N}| \le C_F e^{\gamma_0 j}`$. Then almost surely:
```math
\lim_{j \to \infty} \;\sup_{t \in [j, j+1]} \;\max_{2 \le N \le B e^{\gamma_0(j+1)}} \;\max_{f \in \mathcal F_{j, N}} \; N^\mu \left|\frac{M_{t, N}(f)}{C(T_t)} - M_N(f)\right| \;=\; 0.
\qquad\text{(8)}
```
Consequently, intersecting over countable integers $`B \ge 1`$ (since $`n e^{-\tau_n} \to Y \in (0, \infty)`$ implies eventually $`n^{\gamma_0} \le B e^{\gamma_0(\lfloor \tau_n \rfloor + 1)}`$ for some $`B`$), for any fixed finite family of test functions (such as $`f \equiv 1`$, or $`f_q(F) = \mathbf 1_{\{V(F) \le N^q\}}`$ for a fixed rational $`q \in (0, 1)`$),
```math
\sup_{2 \le N \le n^{\gamma_0}} \; N^\mu \left|\frac{M_{n, N}(f)}{C_n} - M_N(f)\right| \;\longrightarrow\; 0 \qquad \text{almost surely.}
\qquad\text{(9)}
```

**Proof.** Set $`b := 3\gamma_0`$. On each block $`[j, j+1]`$, put $`u_j := (1 - a)j`$ and partition $`[j, j+1]`$ into $`m_j := \lceil e^{b j} \rceil`$ equal subintervals of length $`\delta_j := 1/m_j \le e^{-b j}`$ (yielding $`m_j + 1`$ mesh points).
1. **Mesh concentration:** At each mesh point $`v \in [j, j+1]`$, $`h = v - u_j = a j + O(1)`$, so by Lemma 2.2,
   $`\displaystyle \mathbb E\!\left[\Big(e^{-\mu v}\big(M_{v, N}(f) - \mathbb E[M_{v, N}(f) \mid \mathcal F_{u_j}]\big)\Big)^2\right] \;\le\; K_1 (j+2)^2 e^{(4 - 2\mu)v - 3 u_j} \;=\; O\!\left((j+2)^2 e^{-(2\mu - 1) j + 3 a j}\right).`$
   Note $`3 a = \frac{3}{10}(2\mu - 1)`$, so the decay exponent is $`-\frac{7}{10}(2\mu - 1) j`$. Multiplying by $`N^{2\mu} \le B^{2\mu} e^{2\mu\gamma_0(j+1)}`$ and summing over $`m_j + 1 = O(e^{3\gamma_0 j})`$ mesh points, $`O(e^{\gamma_0 j})`$ cutoffs $`N`$, and $`O(e^{\gamma_0 j})`$ tests $`f \in \mathcal F_{j,N}`$ gives a combined positive exponent of at most $`(5 + 2\mu)\gamma_0 j \le 7\gamma_0 j \le \frac{7}{20}(2\mu - 1)j < \frac{7}{10}(2\mu - 1)j`$. Chebyshev and Borel–Cantelli send the weighted centered mesh error to $`0`$ a.s.
2. **Conditional mean & exponential martingale tail:** By Lemma 2.1, since $`e^{-\mu u_j}(C(T_{u_j}) + c_o(u_j) + j)`$ is bounded a.s. and $`N^\mu e^{-\eta(v - u_j)} \le B^\mu e^{(\mu\gamma_0 - \eta a)j + O(1)} \to 0`$ (as $`\mu\gamma_0 < \eta a/20`$), uniformly on the mesh $`N^\mu |e^{-\mu v} M_{v, N}(f) - d_2 Z^*_{u_j} M_N(f)| \to 0`$ a.s. To replace $`Z^*_{u_j}`$ by $`Z^*_\infty`$ despite $`N^\mu |M_N(f)| \le N^\mu d_N = O(\log N) = O(j)`$: apply `almost-sure.md` (8) at integer $`m`$ with lag $`h = \kappa m`$ where $`0 < 4\kappa < 2\mu - 1`$; using uniform $`\sup_m \mathbb E[(Z^*_m)^2] < \infty`$ (`almost-sure.md` (9)) and $`C \le C_3(1 + \Psi)`$, the orthogonal martingale increments obey $`\mathbb E[(Z^*_{m+1} - Z^*_m)^2] = O(e^{-(2\mu - 1 - 4\kappa) m} + e^{-2\eta\kappa m})`$. Summing over $`m \ge \lfloor u_j \rfloor`$ and applying conditional Jensen between integer times gives $`\mathbb E[(Z^*_\infty - Z^*_{u_j})^2] \le C e^{-2\lambda_0 u_j}`$ for an explicit $`\lambda_0 > 0`$, whence $`(j + 1)|Z^*_{u_j} - Z^*_\infty| \to 0`$ a.s. by Borel–Cantelli! Thus uniformly on the mesh, $`N^\mu |e^{-\mu v} M_{v, N}(f) - d_2 Z^*_\infty M_N(f)| \to 0`$ a.s.
3. **Localized inter-mesh variation without monotonicity:** Fix deterministic integers $`M_* \ge 1, K_* \ge 1`$ and consider the path event $`\mathcal E_{M_*, K_*} := \{\forall t \ge 0,\, C(T_t) \le M_* e^{\mu t} \text{ and } \forall t \ge 1,\, \max_{x \in T_t} c_x(t) \le K_* t\}`$ (which exhausts probability $`1`$ as $`M_*, K_* \to \infty`$ by Claim 6 and `hub-spectrum.md` §3). On each mesh subinterval $`[v, w] \subset [j, j+1]`$ of length $`\delta_j \le e^{-3\gamma_0 j}`$, define the stopping time $`\sigma_{v, K_*} := \inf\{t \ge v : C(T_t) > 2 C(T_v) \text{ or } \max_x c_x(t) > K_*(j+2)\}`$. Because $`e^{-\mu t} C(T_t) \to d_2 Z^*_\infty > 0`$ uniformly on compact steps and $`\delta_j \to 0`$, on $`\mathcal E_{M_*, K_*}`$ we have $`\sigma_{v, K_*} > w`$ eventually for all mesh intervals. A leaf attachment to a vertex of degree $`> N`$ leaves $`\mathcal C(T_t, N)`$ **completely invariant** (the vertex remains deleted and the new leaf is an isolated non-core). Thus $`\mathcal C(T_t, N)`$ changes ONLY at: (i) a promotion anywhere (total rate $`\le C(T_{t-})`$, changing at most $`2`$ component indicators), or (ii) a leaf attachment to a core $`x`$ with $`D_x(t-) \le N`$ (total rate $`\le \sum_{D_x \le N} D_x \le N C(T_{t-})`$, changing at most $`2 + c_x(t-)`$ components: if $`D_x`$ stays $`\le N`$ it modifies only $`x`$'s component, cost $`\le 2`$; if $`D_x`$ crosses $`N`$, it removes $`x`$'s component, possibly leaves one ancestral remnant, and exposes at most $`c_x`$ core-child components). Prior to $`\sigma_{v, K_*}`$, the predictable left-limit rate of relevant events is at most $`2(N+1) C(T_v)`$, so by standard thinning the stopped event count $`E_N(v, w \wedge \sigma_{v, K_*})`$ is conditionally dominated given $`\mathcal F_v`$ by $`\mathrm{Pois}(\Lambda)`$ with $`\Lambda := 2(N+1) C(T_v) \delta_j`$. Poisson Chernoff gives $`\mathbb P(E_N(v, w \wedge \sigma_{v, K_*}) > 2\Lambda + j^2 \mid \mathcal F_v) \le e^{-c j^2}`$, summable over all $`m_j`$ intervals and all $`N \le B e^{\gamma_0(j+1)}`$! On $`\mathcal E_{M_*, K_*}`$ and $`\{E_N(v, w \wedge \sigma_{v, K_*}) \le 2\Lambda + j^2\}`$, for all $`t \in [v, w]`$ and ALL $`|f| \le 1`$ simultaneously:
   $`\displaystyle \big|M_{t, N}(f) - M_{v, N}(f)\big| \;\le\; \big(2 + K_*(j+2)\big)\big(4(N+1) M_* e^{\mu(j+1) - 3\gamma_0 j} + j^2\big) \;=\; o\big(N^{-\mu} e^{\mu j}\big),`$
   since $`N^{1+\mu} e^{-3\gamma_0 j} \le B^{1+\mu} e^{(1+\mu-3)\gamma_0 j + O(1)} \le C e^{-\gamma_0 j}`$. The change in the scaling factor satisfies $`|e^{-\mu t} - e^{-\mu v}| M_{v, N}(1) \le \mu \delta_j e^{-\mu v} C(T_v) = O(e^{-3\gamma_0 j}) = o(N^{-\mu})`$.
4. **Normalization by $`e^{-\mu t} C(T_t)`$:** Because $`N^\mu |M_N(f)| \le N^\mu d_N = O(\log N) = O(j)`$, passing from $`N^\mu |e^{-\mu t} M_{t, N}(f) - d_2 Z^*_\infty M_N(f)| \to 0`$ to (8) requires $`(j + 1)\sup_{t \in [j, j+1]} |e^{-\mu t} C(T_t) - d_2 Z^*_\infty| \to 0`$ a.s. This follows on the exact same mesh $`u_j = (1 - a)j`$, $`\delta_j \le e^{-3\gamma_0 j}`$: by `almost-sure.md` (4) and (7) applied to $`C = A_2`$, the variance at each mesh point $`v`$ multiplied by $`(j+1)^2`$ is $`O((j+1)^2 e^{-(2\mu - 1) j + 3 a j})`$, summable over all $`m_j + 1`$ mesh points; the conditional mean error is $`O((j+1) e^{-\eta a j}) \to 0`$; Step 2 already proved $`(j+1)|Z^*_{u_j} - Z^*_\infty| \to 0`$ a.s.; and monotonicity of $`C(T_t)`$ across $`[v, w]`$ introduces relative factor $`e^{\pm \mu \delta_j} = 1 + O(e^{-3\gamma_0 j})`$, negligible after multiplying by $`(j+1)`$. Dividing $`e^{-\mu t} M_{t, N}(f)`$ by $`e^{-\mu t} C(T_t) = d_2 Z^*_\infty + o(j^{-1})`$ now yields (8) cleanly, and passing to $`\tau_n`$ proves (9). $`\square`$

### Deduction of Theorem 1 on the Engine B Subsequence $`\{n : N(n) \le n^{\gamma_0}\}`$
Along any subsequence with $`N(n) \to \infty`$ and $`N(n) \le n^{\gamma_0}`$: (9) for $`f \equiv 1`$ together with Astra's $`d_N = M_N(1) \sim A N^{-\mu}\log N`$ and $`C_n / n^\mu \to X`$ gives $`M_{n, N(n)} / (C_n d_{N(n)}) \to 1`$, proving (1). For each fixed rational $`q \in (0, 1)`$, (9) for $`f_q(F) = \mathbf 1_{\{V(F) \le N^q\}}`$ gives $`\nu_{n, N(n)}(V \le N(n)^q) - \nu_{N(n)}(V \le N(n)^q) \to 0`$ a.s.; combining with Astra's $`\|\nu_N - \mu_N\|_{\mathrm{TV}} \to 0`$ and $`\mu_N(V \le N^q) \to q`$ (`open/forest-age/README.md` eq. (2)), intersecting over rational $`q \in (0, 1)`$, and using CDF monotonicity with endpoints $`0, 1`$ proves (2)!

---

## 3. Engine A: Large Cutoffs $`N(n) > n^{\gamma_0}`$ with $`N(n)/n \to 0`$

Now let $`N = N(n)`$ be any deterministic sequence with $`N(n) > n^{\gamma_0}`$ and $`N(n)/n \to 0`$ (combining with Engine B on $`\{n : N(n) \le n^{\gamma_0}\}`$ covers EVERY deterministic sequence $`N(n) \to \infty, N(n)/n \to 0`$!). Put $`L := \log N(n) \ge \gamma_0 \log n`$ and $`t_n := \log n`$, so $`L \asymp t_n`$ and $`t_n - L = \log(n/N(n)) \to \infty`$.

Fix any small rational $`\varepsilon \in (0, 1/4)`$ and write $`\omega_\varepsilon := \log(1 + \varepsilon) \in (0, \varepsilon)`$. Since $`Y \sim \mathrm{Gamma}(3, 1)`$ is continuous, $`\mathbb P(Y = (1+\varepsilon)^r) = 0`$ for every $`r \in \mathbb Z`$. On the event $`\{Y \in ((1+\varepsilon)^r, (1+\varepsilon)^{r+1})\}`$ (whose countable union over $`r \in \mathbb Z`$ has probability $`1`$), for all sufficiently large $`n`$ the exact epoch $`\tau_n`$ lies strictly inside the deterministic interval $`[t_-, t_+]`$ of exact width $`\omega_\varepsilon`$:
```math
t_- := \log n - (r + 1)\omega_\varepsilon \;\le\; \tau_n \;\le\; \log n - r\omega_\varepsilon =: t_+ = t_- + \omega_\varepsilon.
```
Set $`b := 4 \log t_- = O(\log \log n) = o(L)`$ (so $`e^b = t_-^4`$) and $`s_0 := t_- - L + b`$ (note $`s_0 \to \infty`$ since $`t_- - L \to \infty`$). Let $`H(t) := \{v \in K(T_t) : D_v(t) > N\}`$; monotonicity of degrees gives $`H(t_-) \subseteq H(\tau_n) \subseteq H(t_+)`$.

### Step 3.1 (Finite-Horizon Hub Localization at $`s_0 = t_- - L + b`$)
Set $`\delta_n := t_-^{-1/4}`$ and define the $`\mathcal F_{s_0}`$-measurable hub sets
```math
H_{s_0}^- := \big\{v \in K(T_{s_0}) : D_v(s_0) > (1 + \delta_n) e^b\big\}, \qquad
H_{s_0}^+ := \big\{v \in K(T_{s_0}) : D_v(s_0) > (1 - \delta_n) e^{-\omega_\varepsilon} e^b\big\}.
```
1. **No late-blooming hubs:** For any $`v`$ alive at $`s_0`$ with $`D_v(s_0) = d \le t_-^2`$, over duration $`h_+ := t_+ - s_0 = L - b + \omega_\varepsilon`$, $`D_v(t_+)`$ is dominated by $`\mathrm{NB}(d, p)`$ with $`p = e^{-h_+}`$. Since a $`\mathrm{Geom}(p)`$ variable satisfies $`\mathbb E[e^{(p/2)\mathrm{Geom}(p)}] = \frac{p e^{p/2}}{1 - (1-p)e^{p/2}} \le 2`$ (as $`(1 - p/2)e^{p/2} \le 1`$ for $`p \in (0, 1]`$), Markov gives $`\mathbb P(D_v(t_+) > N \mid \mathcal F_{s_0}) \le 2^d e^{-p N/2} = 2^d \exp(-\frac{1}{2} e^{b - \omega_\varepsilon}) \le \exp(-\frac{1}{4} t_-^4)`$. For vertices born at $`r \in (s_0, t_+]`$, $`D_v(t_+)`$ is dominated by $`\mathrm{Geom}(e^{-(t_+ - r)})`$, giving $`\mathbb P(D_v(t_+) > N \mid \tau_v = r) \le \exp(-e^{b - \omega_\varepsilon})`$. Compensating births over $`(s_0, t_+]`$ ($`\mathbb E[|T_{t_+}|] = O(n)`$) gives total failure probability $`O(n \exp(-c t_-^4))`$, summable in $`n`$.
2. **Finite-horizon vanishing-tolerance concentration:** Fix deterministic integer $`K_* \ge 1`$ and localize on $`\{\max_{x \in T_{s_0}} c_x(s_0) \le K_* t_-\}`$. For any $`v \in K(T_{s_0})`$ with $`D = D_v(s_0) > t_-^2`$ and $`c = c_v(s_0) \le K_* t_- \ll \delta_n D`$, apply the tagged-leaf construction of `hub-spectrum.md` §2 over finite elapsed time $`h \in \{L - b, L - b + \omega_\varepsilon\}`$: conditional on $`\mathcal F_{s_0}`$, $`e^{-h} D_v(s_0 + h)`$ is sandwiched between a normalized $`\mathrm{Yule}(D - c)`$ minus a normalized rate-1-immigrating Yule population $`E_h`$ ($`\mathbb E[e^{E_h/2}] \le 2`$) and a normalized $`\mathrm{Yule}(D)`$. Each normalized $`\mathrm{Yule}(1)`$ individual $`G = e^{-h}\mathrm{Geom}(e^{-h})`$ has $`\mathbb E[G] = 1`$ and uniform MGF $`\mathbb E[e^{G/2}] \le 2`$, so by Taylor's theorem its centered log-MGF obeys $`\log \mathbb E[e^{\theta(G - 1)}] \le C \theta^2`$ uniformly for $`|\theta| \le 1/4`$ and all $`h > 0`$. Setting $`\theta \propto \delta_n`$ (plus the exponential tail $`\mathbb P(E_h > \delta_n D / 16) \le 2 e^{-\delta_n D / 32}`$) yields $`\mathbb P(|e^{-h}D_v(s_0 + h)/D - 1| > \delta_n / 4 \mid \mathcal F_{s_0}) \le K \exp(-c \delta_n^2 D) \le K \exp(-c t_-^{3/2})`$. Summing over $`\mathbb E[|T_{s_0}|] = O(n)`$ gives $`O(n \exp(-c t_-^{3/2}))`$, summable in $`n`$!
3. **Sandwich & asymptotic hub counts:** By Borel–Cantelli (exhausting $`K_* \to \infty`$), almost surely eventually,
   $`\displaystyle H_{s_0}^- \;\subseteq\; H(t_-) \;\subseteq\; H(\tau_n) \;\subseteq\; H(t_+) \;\subseteq\; H_{s_0}^+. \qquad\text{(10)}`$
   Evaluating the common-event intermediate degree law of Claim 11 (`hub-spectrum.md` §4) at the integer tree size $`|T_{s_0}| \sim Y e^{s_0}`$ with integer ceiling/floor rounding on thresholds $`(1+\delta_n)e^b`$ and $`(1-\delta_n)e^{-\omega_\varepsilon}e^b`$ (using $`e^b \to \infty`$ and $`e^b / e^{s_0} = e^{L - t_-} \to 0`$) yields almost surely:
   $`\displaystyle \frac{|H_{s_0}^-|}{A X (n/N)^\mu} \;\to\; \left(\frac{Y}{(1+\varepsilon)^{r+1}}\right)^\mu \in \big[(1+\varepsilon)^{-\mu}, 1\big], \qquad \frac{|H_{s_0}^+|}{A X (n/N)^\mu} \;\to\; \left(\frac{Y}{(1+\varepsilon)^r}\right)^\mu \in \big[1, (1+\varepsilon)^\mu\big]. \qquad\text{(11)}`$

### Step 3.2 (Negligibility of Pre-$`s_0`$ Children of $`H_{s_0}^+`$ at Scale $`K L`$)
We show $`\sum_{v \in H_{s_0}^+} c_v(s_0) = o\big(L |H_{s_0}^+|\big)`$ a.s. (the founder contributes at most $`c_o(s_0) = O(t_-)`$). Fix deterministic integers $`M_*, K_* \ge 1`$ and localize on $`\mathcal E_{M_*, K_*}`$. Set $`M := \lceil 4 b / \mu \rceil = O(\log t_-) = o(L)`$, and split $`v \in H_{s_0}^+`$ by **core promotion time** $`\tau_v^{\text{core}}`$ (counting the founder separately):
- **Old cores ($`\tau_v^{\text{core}} \le s_0 - M`$; empty if $`s_0 \le M`$):** Let $`P_{\text{old}}(s_0)`$ count all promotions up to $`s_0`$ by cores promoted by time $`s_0 - M`$. Fix any deterministic integer $`D \ge 1`$; by (11), the countable union over $`D \ge 1`$ of the eventual events $`\{|H_{s_0}^+| \ge e^{\mu(s_0 - b)} / D\}`$ has probability $`1`$. Prior to the deterministic localization exit $`\sigma_{M_*} := \inf\{u \ge 0 : C(T_u) > M_* e^{\mu u}\}`$, the predictable compensator of $`P_{\text{old}}(s_0 \wedge \sigma_{M_*})`$ over $`[0, s_0]`$ is bounded deterministically by $`\Lambda_{\text{old}} := M_* \int_0^{s_0 - M} e^{\mu u}\,du + M M_* e^{\mu(s_0 - M)} \le C M_* M e^{\mu(s_0 - M)} = O(M e^{-(4 - \mu)b} e^{\mu(s_0 - b)})`$. Set the **deterministic** threshold $`T_n := \varepsilon L e^{\mu(s_0 - b)} / D`$. Since $`s_0 - b = t_- - L \to \infty`$ and $`L \ge \gamma_0 \log n`$, for all sufficiently large $`n`$ we have $`T_n \ge 2 e \Lambda_{\text{old}}`$ AND $`(\ln 2) T_n \ge 3 \log n`$ deterministically! Hence Poisson Chernoff on the stopped count gives $`\mathbb P(P_{\text{old}}(s_0 \wedge \sigma_{M_*}) \ge T_n) \le 2^{-T_n} = e^{-(\ln 2) T_n} \le n^{-3}`$, summable in $`n`$. On $`\{|H_{s_0}^+| \ge e^{\mu(s_0 - b)}/D\} \cap \mathcal E_{M_*, K_*}`$, we have $`T_n \le \varepsilon L |H_{s_0}^+|`$ and $`C(T_{s_0 - M}) \le M_* e^{\mu(s_0 - M)} = o(|H_{s_0}^+|)`$. Sending rational $`\varepsilon \downarrow 0`$ shows the old-core contribution is $`o(L |H_{s_0}^+|)`$ a.s.!
- **Young cores ($`\tau_v^{\text{core}} \in (s_0 - M, s_0]`$):** Group $`t_- \in [j, j+1]`$ ($`j = \lfloor t_- \rfloor \ge 1`$); then $`b = 4 \log j + O(1/j)`$ and $`s_0 \in [m-1, m]`$ for one of at most $`j + 3`$ integer right-endpoints $`m \in \{0, 1, \dots, j+2\}`$. Put $`M_j := \lceil 4(4\log j + 1)/\mu \rceil + 2`$, so $`(s_0 - M, s_0] \subseteq (m - M_j, m]`$ and $`c_v(s_0) \le c_v(m)`$ by monotonicity. Truncating at $`10 M_j + 1`$:
  $`\displaystyle \sum_{\substack{v \in H_{s_0}^+ \\ \tau_v^{\text{core}} \in (s_0 - M, s_0]}} c_v(s_0) \;\le\; (10 M_j + 1)|H_{s_0}^+| + \sum_{\tau_v^{\text{core}} \in (m - M_j, m]} c_v(m) \mathbf 1_{\{c_v(m) > 10 M_j + 1\}}.`$
  The first term is $`O(M_j |H_{s_0}^+|) = o(L |H_{s_0}^+|)`$. For the tail sum at a fixed integer pair $`(j, m)`$, every core promoted at $`\tau_v^{\text{core}} \in (m - M_j, m]`$ starts at core degree $`1`$ and has $`c_v(m) - 1 \preceq \mathrm{Pois}(M_j)`$ by strong Markov; exact Poisson Chernoff at $`10 M_j`$ has exponent $`-(10\log 10 - 9) M_j < -14 M_j`$, absorbing $`(1 + M_j)`$ to yield expectation $`\le C \mathbb E[C(T_m)] e^{-10 M_j} \le C' e^{\mu m - 40 \log j / \mu}`$. By Markov, the probability that this tail sum exceeds $`e^{\mu(m - 4\log j)}`$ is at most $`C' j^{4\mu - 40/\mu} \le C' j^{-30}`$. Summing over all $`j + 3`$ integer bins $`m`$ gives $`O(j^{-29})`$, **summable over all integers $`j \ge 1`$**! Since $`e^{\mu(m - 4\log j)} = O(e^{\mu(s_0 - b)}) = O(|H_{s_0}^+|) = o(L |H_{s_0}^+|)`$, Borel–Cantelli yields $`\sum_{v \in H_{s_0}^+} c_v(s_0) = o(L |H_{s_0}^+|)`$ a.s.!

### Step 3.3 (Simultaneous Independent Poisson–Inert Coupling on $`[s_0, t_+]`$)
Condition on $`\mathcal F_{s_0}`$. Every hub $`v \in H_{s_0}^+`$ is ALREADY a core at time $`s_0`$, so the leaf branches promoted by distinct $`v \in H_{s_0}^+`$ during $`(s_0, t_+]`$ are **pairwise disjoint** and contain **no** hub in $`H_{s_0}^+`$. Extending `hub-children.md` §3 simultaneously to all $`v \in H_{s_0}^+`$: attach to each $`v`$'s independent rate-1 initiating clock on $`[s_0, t_+]`$ a fresh independent sequence of inert virgin trajectories $`S^{(v, j)}(u)`$ ($`u \ge 0`$), using fresh disjoint internal drivers for each branch while routing outgoing child-to-parent arrows on a separate channel that updates the exterior without altering the ideal branch. Conditional on $`\mathcal F_{s_0}`$, the resulting ideal Poisson–inert marks across all $`v \in H_{s_0}^+`$ are **mutually independent**!

Moreover, on the event (10) (every vertex of degree $`> N`$ at $`\tau_n`$ belongs to $`H_{s_0}^+`$ and was already a core at $`s_0`$), **no** post-$`s_0`$ promoted branch can contain any deleted vertex of degree $`> N`$ as a proper descendant! Thus any post-$`s_0`$ promoted branch that remains unseparated from its ideal trajectory AND satisfies $`|S^{(v,j)}(t_+ - \sigma_{v,j})| \le N`$ has **zero** internal degree deletions throughout $`[\sigma_{v,j}, t_+]`$, so its physical component $`F_N(w; \tau_n)`$ coincides with its full ideal tree $`S^{(v,j)}(\tau_n - \sigma_{v,j})`$ and has monotonically increasing ordinary volume $`V(t_-) \le V(\tau_n) \le V(t_+)`$!

Control the two discrepancy counts without conditioning on any future event:
1. **Single stopped count of hub-to-existing-core sends:** Fix deterministic integer $`K_* \ge 1`$ and define the stopping time $`\eta_{K_*} := \inf\{u \ge s_0 : \max_{x \in T_u} c_x(u) > 2 K_* t_-\}`$. On $`\mathcal E_{M_*, K_*}`$, for all $`u \in [s_0, t_+]`$ and large $`n`$ we have $`u \le t_+ = t_- + \omega_\varepsilon < 2 t_-`$, so $`\max_x c_x(u) \le K_* u < 2 K_* t_-`$ and therefore $`\eta_{K_*} > t_+`$ eventually! Each clock ring of $`v \in H_{s_0}^+`$ that chooses an existing core neighbor causes **at most $`2`$ discrepancies** (at most $`1`$ missed promotion AND at most $`1`$ first separation arrow into an already promoted child branch). Before $`\eta_{K_*}`$, since $`D_v(u) \ge D_v(s_0) \ge \frac{1}{2} t_-^4`$, the total predictable rate of hub-to-existing-core sends across all $`v \in H_{s_0}^+`$ is at most $`|H_{s_0}^+| \cdot 2 K_* t_- / (\frac{1}{2} t_-^4) = 4 K_* t_-^{-3} |H_{s_0}^+|`$. Thus the stopped send count $`S_{\text{bad}}`$ over $`[s_0, t_+ \wedge \eta_{K_*}]`$ is conditionally dominated given $`\mathcal F_{s_0}`$ by $`\mathrm{Pois}(4 K_* (L - b + \omega_\varepsilon) t_-^{-3} |H_{s_0}^+|)`$.
2. **Oversized ideal branches at $`t_+`$:** By independent marking given $`\mathcal F_{s_0}`$, the count $`O_{\text{bad}}`$ of ideal points in $`[s_0, t_+]`$ with $`|S^{(v,j)}(t_+ - \sigma_{v,j})| > N`$ is **genuinely Poisson** of mean $`\le |H_{s_0}^+| \int_0^{L - b + \omega_\varepsilon} \mathbb P(|S_u| > N)\,du \le |H_{s_0}^+| (L - b + \omega_\varepsilon) \cdot 2 e^{\omega_\varepsilon} t_-^{-4}`$ (by Markov on $`\mathbb E[|S_u|] \le 2 e^u \le 2 N e^{\omega_\varepsilon} t_-^{-4}`$).

On the $`\mathcal F_{s_0}`$-measurable event $`\{|H_{s_0}^+| \ge K_0\}`$ (with deterministic $`K_0`$ chosen so $`(\ln 2)\varepsilon K_0 \gamma_0 > 3`$), the threshold $`\varepsilon L |H_{s_0}^+|`$ exceeds $`2 e`$ times both Poisson means for all large $`n`$, so Poisson Chernoff gives $`\mathbb P(S_{\text{bad}} > \varepsilon L |H_{s_0}^+| \text{ or } O_{\text{bad}} > \varepsilon L |H_{s_0}^+|,\, |H_{s_0}^+| \ge K_0 \mid \mathcal F_{s_0}) \le 2 \cdot 2^{-\varepsilon L |H_{s_0}^+|} = 2 e^{-(\ln 2)\varepsilon L |H_{s_0}^+|} \le 2 n^{-3}`$, summable in $`n`$! Since $`|H_{s_0}^+| \to \infty`$ a.s., Borel–Cantelli ensures $`2 S_{\text{bad}} + O_{\text{bad}} \le 3\varepsilon L |H_{s_0}^+|`$ eventually a.s.

### Step 3.4 (Compound-Poisson Concentration & Sandwich to $`\tau_n`$)
Write $`K_- := |H_{s_0}^-|`$, $`K_+ := |H_{s_0}^+|`$, $`L_- := L - b`$, and $`L_+ := L - b + \omega_\varepsilon`$ (so for all large $`n`$, $`L_\pm \ge \frac{1}{2} \gamma_0 \log n`$ and $`L_\pm / L \to 1`$ deterministically). For Poisson Chernoff on any $`\mathrm{Pois}(\lambda)`$ with $`\lambda \le K_\pm L_+`$, we have the uniform additive tail $`\mathbb P(|\mathrm{Pois}(\lambda) - \lambda| > \varepsilon K_\pm L) \le 2 \exp(-c \varepsilon^2 K_\pm L)`$ (with an absolute constant $`c > 0`$ independent of $`\lambda \le K_\pm L_+`$ and $`x`$!). Choose deterministic $`K_0`$ so that $`c \varepsilon^2 K_0 \gamma_0 > 3`$; both $`\{K_- < K_0\}`$ and $`\{K_+ < K_0\}`$ occur only finitely often a.s. by (11).
1. **Component count (1):** On (10) and $`\{\eta_{K_*} > t_+\}`$, $`M_{\tau_n, N}`$ is sandwiched between the ideal point count on $`H_{s_0}^- \times [s_0, t_-]`$ ($`\sim \mathrm{Pois}(K_- L_-)`$) minus $`2 S_{\text{bad}} + O_{\text{bad}}`$, and the ideal point count on $`H_{s_0}^+ \times [s_0, t_+]`$ ($`\sim \mathrm{Pois}(K_+ L_+)`$) plus $`1 + c_o(s_0) + \sum_{v \in H_{s_0}^+} c_v(s_0) + 2 S_{\text{bad}}`$. On $`\{K_\pm \ge K_0\}`$, additive Poisson Chernoff gives failure probability $`\le 2\exp(-c\varepsilon^2 K_\pm L) \le 2 n^{-3}`$, summable in $`n`$. Combining with (11) and Steps 3.2–3.3 and letting rational $`\varepsilon \downarrow 0`$ proves (1)!
2. **Rescaled log-volume law (2):** Fix any rational $`x \in (0, 1)`$. On (10) and $`\{\eta_{K_*} > t_+\}`$:
   - *Lower bound on $`M_{\tau_n, N}(\mathbf 1_{\{V \le N^x\}})`$:* take all ideal marks on $`H_{s_0}^- \times [s_0, t_-]`$ with $`V(S^{(v,j)}(t_+ - \sigma_{v,j})) \le N^x`$ (which implies $`V(\tau_n) \le V(t_+) \le N^x < N`$), minus $`2 S_{\text{bad}} + O_{\text{bad}}`$. Conditional on $`\mathcal F_{s_0}`$, this ideal count is Poisson with mean $`\lambda_- \le K_- L_-`$ equal to $`K_- \int_{\omega_\varepsilon}^{L_+} \mathbb P(|S_u| \le N^x)\,du = (x + o(1)) K_- L`$ (by the Yule comparison $`1 + \mathrm{Yule}(1) \le |S_u| \le \mathrm{Yule}(2)`$ in `open/forest-age/README.md`).
   - *Upper bound on $`M_{\tau_n, N}(\mathbf 1_{\{V \le N^x\}})`$:* take all ideal marks on $`H_{s_0}^+ \times [s_0, t_-]`$ with $`V(S^{(v,j)}(t_- - \sigma_{v,j})) \le N^x`$, **plus ALL ideal births on $`H_{s_0}^+ \times (t_-, t_+]`$** (which are Poisson of mean $`K_+ \omega_\varepsilon \le \varepsilon K_+ L`$, controlled at additive tolerance $`\varepsilon K_+ L`$ by Chernoff!), plus $`1 + c_o(s_0) + \sum_{v \in H_{s_0}^+} c_v(s_0) + 2 S_{\text{bad}} + O_{\text{bad}}`$. The main ideal count on $`[s_0, t_-]`$ is Poisson with mean $`\lambda_+ \le K_+ L_-`$ equal to $`K_+ \int_0^{L_-} \mathbb P(|S_u| \le N^x)\,du = (x + o(1)) K_+ L`$.
   Applying additive Poisson Chernoff on $`\{K_\pm \ge K_0\}`$ (failure probability $`\le 2 n^{-3}`$, summable in $`n`$), dividing by $`M_{\tau_n, N}`$ using (1), intersecting over all rational $`x \in (0, 1)`$ and rational $`\varepsilon \downarrow 0`$, and using CDF monotonicity with boundary values $`0, 1`$ proves (2)! $`\square`$

---

## 4. Proof of Theorem 2 (Separated Simultaneous Total-Variation Regimes)

### Proof of Theorem 2(1) (Simultaneous TV Convergence for $`N(n) \le c_0 \log n`$)
Set $`c_0 := \frac{\gamma_0}{2 \log 4} > 0`$. By canonical recursively sorted balanced-parenthesis encodings (injective on rooted unlabelled trees), the set $`\mathcal R_N`$ of all rooted unlabelled trees of ordinary size $`|R| \le N`$ has cardinality $`|\mathcal R_N| \le \sum_{m=1}^N 4^m \le 4^{N+1}`$. On each countable $`Y`$-bracket $`\{Y \in (B^{-1}, B)\}`$, whenever $`N(n) \le c_0 \log n`$ we have $`N(n) \le c_0(j + 1 + \log B)`$ on $`\tau_n \in [j, j+1]`$ and hence $`|\mathcal R_N| \le 4 B^{c_0 \log 4} e^{(\gamma_0/2) j}`$. Repeat the proof of Proposition 2.3 on $`\mathcal F_{j, N} := \{1\} \cup \{\mathbf 1_{\{F = R\}} : R \in \mathcal R_N\}`$ multiplying the Chebyshev weight on single-tree indicators by $`4^{N+1}`$:
- In Step 1, the squared extra weight adds $`\gamma_0 j + O(1)`$ to the positive exponent, giving at most $`(6 + 2\mu)\gamma_0 j \le 8\gamma_0 j \le \frac{8}{20}(2\mu - 1)j < \frac{7}{10}(2\mu - 1)j`$ (still summable!).
- In Step 2, the conditional mean error picks up $`4^{N+1} = O(e^{(\gamma_0/2) j})`$ while $`e^{-\eta a j}`$ decays at rate $`\eta a \ge 20\gamma_0`$; when replacing $`Z^*_{u_j}`$ by $`Z^*_\infty`$ in the $`\ell^1`$ sum over $`R \in \mathcal R_N`$, use $`\sum_{R \in \mathcal R_N} M_N(\mathbf 1_R) \le d_N = O(N^{-\mu}\log N)`$ directly (no cardinality factor!).
- In Step 3, every single relevant event changes at most $`2 + c_x`$ components and hence changes the $`\ell^1`$ sum $`\sum_{R \in \mathcal R_N} |M_{t, N}(\mathbf 1_R) - M_{v, N}(\mathbf 1_R)|`$ by at most $`2(2 + c_x)`$ (again, NO $`|\mathcal R_N|`$ factor!).
Therefore almost surely,
```math
\frac{M_{n, N(n)}}{C_n d_{N(n)}} \to 1, \qquad \sum_{R \in \mathcal R_{N(n)}} \left|\frac{M_{n, N(n)}(\mathbf 1_R)}{C_n d_{N(n)}} - \frac{M_{N(n)}(\mathbf 1_R)}{d_{N(n)}}\right| \to 0.
```
Normalizing by $`M_{n, N(n)} / (C_n d_{N(n)}) \to 1`$ gives $`\sum_{R \in \mathcal R_{N(n)}} |\nu_{n, N(n)}(\{R\}) - \nu_{N(n)}(\{R\})| \to 0`$ a.s.

For the complement $`\{|F| > N\}`$ (for $`N \ge 3`$ so $`\log N \ge 1`$): under $`\mu_N`$, on $`a \in [0, \log N - \sqrt{\log N}]`$ we have $`\mathbb E_J^0[|S_a|] \le 2 e^a`$, so Markov gives $`P_a^0(J, |F| > N) \le 2 e^{-\sqrt{\log N}}`$, while the remaining interval $`[\log N - \sqrt{\log N}, \log N]`$ has weight $`1/\sqrt{\log N}`$. Hence $`\mu_N(|F| > N) \le \frac{1}{\sqrt{\log N}} + 2 e^{-\sqrt{\log N}} \to 0`$. By Astra's sequential limit $`\|\nu_N - \mu_N\|_{\mathrm{TV}} \to 0`$, $`\nu_N(|F| > N) \to 0`$, whence $`\nu_{n, N(n)}(|F| > N(n)) = 1 - \sum_{R \in \mathcal R_{N(n)}} \nu_{n, N(n)}(\{R\}) \to 0`$ a.s. Thus $`\|\nu_{n, N(n)} - \nu_{N(n)}\|_{\mathrm{TV}} \le \frac{1}{2}\sum_{R \in \mathcal R_{N(n)}} |\nu_{n, N(n)}(\{R\}) - \nu_{N(n)}(\{R\})| + \frac{1}{2}\big(\nu_{n, N(n)}(|F|>N(n)) + \nu_{N(n)}(|F|>N(n))\big) \to 0`$ a.s., and the triangle inequality with $`\|\nu_{N(n)} - \mu_{N(n)}\|_{\mathrm{TV}} \to 0`$ proves (3). $`\square`$

### Proof of Theorem 2(2) (Total-Variation Breakdown via Pure-Birth Holding Time)
In the inert planted virgin process $`S_a`$ started from $`J`$ ($`|S_0| = 2`$), every intrinsic jump increases ordinary size $`|S_a|`$ by $`+1`$ at exact rate $`r_0(S_a) = |S_a| - 1/D_o(S_a) \ge |S_a| - 1/2`$ (`finite-kernel.md` §2, since $`D_o(S_a) \ge 2`$ forever). Because no intrinsic transition preserves size, the single residence time at size $`v \ge 2`$ given the entrance state is $`\mathrm{Exp}(v - 1/D_o)`$ with expectation $`\le 1/(v - 1/2)`$ (every size $`v \ge 2`$ is visited once a.s. with no explosion by Yule domination). Mixing over entrance states gives the exact point-mass bound
```math
\mu_N(V = v) \;=\; \frac{1}{\log N} \int_0^{\log N} \mathbb P_J^0(|S_a| = v)\,da \;\le\; \frac{1}{(v - 1/2)\log N} \qquad \forall v \ge 2.
\qquad\text{(12)}
```
Let $`S_{n, N} := \{V(F) : F \in \mathcal C_{n, N}\} \subset \mathbb Z_{\ge 1}`$ be the empirical support of component vertex counts (which may include size $`1`$); clearly $`|S_{n, N}| \le M_{n, N}`$. For $`N \ge 3`$ (so $`\log\log N > 0`$), set $`v_n := \lceil M_{n, N} \cdot \frac{\log\log N}{\log N} \rceil + 1 \ge 2`$ and define the pathwise event $`E_n := \{F : V(F) \in S_{n, N} \text{ and } V(F) > v_n\}`$ (every $`v \in E_n`$ satisfies $`v \ge 3 > 2`$). Then:
1. By (12),
   $`\displaystyle \mu_N(E_n) \;\le\; \sum_{v \in S_{n, N},\, v > v_n} \frac{1}{(v - 1/2)\log N} \;\le\; \frac{M_{n, N}}{(v_n - 1/2)\log N} \;\le\; \frac{1}{\log\log N} \;\longrightarrow\; 0.`$
2. Under $`\nu_{n, N}`$, every component has $`V(F) \in S_{n, N}`$, so $`\nu_{n, N}(E_n) = \nu_{n, N}(V > v_n)`$.
3. By Theorem 1(1), $`M_{n, N} \sim A X (n/N)^\mu \log N`$ a.s. (diverging to $`\infty`$), so
   $`\displaystyle \frac{\log v_n}{\log N} \;=\; \frac{\mu \log(n/N) + \log(A X) + \log\log\log N + o(1)}{\log N} \;=\; \mu\!\left(\frac{\log n}{\log N} - 1\right) + o(1) \quad \text{a.s.}`$
   When $`\alpha = \liminf \frac{\log N}{\log n} > \frac{\mu}{1+\mu}`$, we have $`\limsup \frac{\log v_n}{\log N} = \mu(\alpha^{-1} - 1) \in [0, 1)`$, so by Theorem 1(2) $`\liminf \nu_{n, N}(V > v_n) \ge 1 - \mu(\alpha^{-1} - 1)`$ a.s. Hence
   $`\displaystyle \liminf_{n\to\infty} \big\|\nu_{n, N} - \mu_N\big\|_{\mathrm{TV}} \;\ge\; \liminf_{n\to\infty} \big(\nu_{n, N}(E_n) - \mu_N(E_n)\big) \;\ge\; 1 - \mu\!\left(\frac{1}{\alpha} - 1\right) \;>\; 0 \qquad \text{a.s.}`$
   When $`\alpha = 1`$ ($`N(n) = n^{1 - o(1)}`$), since TV distance is at most $`1`$, $`\|\nu_{n, N} - \mu_N\|_{\mathrm{TV}} \to 1`$ a.s. $`\square`$

---

## 5. Replay & Exact Verification
```sh
python3 -B round2/proofs/neon2/q5/verify_q5_exact.py
g++ -std=c++17 -O3 -Wall -Wextra -pedantic round2/proofs/neon2/q5/sim_q5.cpp -o ./sim_q5 && ./sim_q5
```
