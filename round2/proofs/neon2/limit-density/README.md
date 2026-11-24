# Neon-2 Claim 1: Absolute Continuity and Everywhere-Positive Lower Semicontinuous Density of $`X = \lim_{n\to\infty} C_n / n^\mu`$

**Pinned claim commit:** `47a79911dcc46ee9d54551aebb97a00a2b83b852`.  
**Review and dependencies:** GPT-6.2 Astra and GPT-6.2 Sol audited the submitted arguments; Claude independently checked absolute continuity and a shorter positivity argument with different premises. See [Claim 1 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-1.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-1.md).

**Authors / cross-family verification:**
- Original theorem, edge-split rate-cancellation proof ($`\lambda_m + 1/D_{\mathrm{root}} = m`$), and countable hub-hub cover: Neon-2 lead (Gemini 4 Neon).
- Independent cold audit, Claim-6-only refinement (removing the `forest-age/inert` dependency for absolute continuity), conditional full-support & everywhere-positive lower semicontinuous density version, multi-cut disjointification & regularity obstruction analysis, and exact rational + Monte Carlo verification (`check_density.py`): GPT-6.2 Sol Ultrafast (`REPORT.md`).
- Final proof audit and script replay: Neon-2 lead (Gemini 4 Neon).

**Explicit dependencies:**
- **Part (a) (Absolute continuity of $`X`$):** Uses ONLY Addario-Berry et al. Theorem 3.1 ($`Z_1 + Z_2 > 0`$ a.s.), Astra Claim 6 ($`C_n / n^\mu \to X \in (0, \infty)`$ a.s., `round2/standalone/claim-6.md`), and Astra's fixed-hub children construction (`round2/fringe/hub-children.md` §§1, 3, proving that a hub produces infinitely many hub children a.s.).
- **Part (b) (Positive compact-local lower bounds and everywhere-positive l.s.c. density version):** Additionally uses the inert-slot planted mean coefficient $`\Phi`$ and uniform one-leaf gradient $`L_0`$ (`round2/open/forest-age/inert/README.md` §§1–3, or equivalently our Track 3 grafted re-rooting bound).

**What is NOT proved:** Boundedness, continuity ($`C^0`$), or higher smoothness ($`C^k / C^\infty`$) of the full density $`p_X`$ on $`(0, \infty)`$ remain open (see Section 4 for the exact summability gap across the multi-cut partition).

---

## Theorem
Let $`X = \lim_{n\to\infty} C_n / n^\mu`$ be the almost-sure scaling limit of the core count in the undirected random friend tree started from one edge $`T_2`$ (Astra Claim 6). Then:
1. **(Absolute continuity):** The law of $`X`$ is absolutely continuous with respect to Lebesgue measure on $`(0, \infty)`$. In particular, $`X`$ has **no atoms** and **no singular continuous component**.
2. **(Everywhere-positive lower semicontinuous density):** There exists a lower semicontinuous extended-valued version $`p_X : (0, \infty) \to (0, \infty]`$ of the density of $`X`$ (finite almost everywhere) which is strictly positive everywhere on $`(0, \infty)`$ and bounded below by a strictly positive constant $`c(K) > 0`$ on every compact subinterval $`K \subset (0, \infty)`$.

---

## 1. Countable Almost-Sure Cover by Uncrossed Core–Core Cuts
Embed the discrete chain $`(T_n)_{n \ge 2}`$ in continuous time with independent rate-one initiating clocks at all vertices (total jump rate at size $`m`$ equals $`m`$; holding times $`E_m \sim \mathrm{Exp}(m)`$ are independent of the discrete jump chain, and $`e^{-t}|T_t| \to W_{\mathrm{size}} \in (0, \infty)`$ a.s.).

By Addario-Berry et al. Theorem 3.1, at least one initial vertex $`u \in \{1, 2\}`$ is a hub ($`Z_u > 0`$, hence $`D_u(t) \sim W_u e^t`$ with $`W_u > 0`$). By `round2/fringe/hub-children.md` §3 (coupling $`u`$'s rate-one clock points before conditioning on $`\{Z_u > 0\}`$ to i.i.d. silent planted newborns with hub probability $`q \ge 1/2`$, while the total sending rate $`c_u(t)/D_u(t) = O((1+t)e^{-t})`$ from $`u`$ to existing core neighbours has finite integral on $`\{Z_u > 0\}`$), $`u`$ almost surely produces infinitely many hub core children $`v`$ ($`Z_v > 0`$).

For any pair of adjacent hubs $`e = (u, v)`$, the directional crossing compensator $`\int_0^\infty (1/D_u(t) + 1/D_v(t))\,dt`$ is finite almost surely, so by localization at bounded compensator levels the edge $`e`$ is crossed by arrows only finitely many times almost surely. Choosing an integer size $`n \ge 4`$ after the last crossing and after both $`u`$ and $`v`$ are cores, the countable union
```math
\bigcup_{n \ge 4}\ \bigcup_{T_0\text{ of size }n}\ \bigcup_{e \in E_{\mathrm{core}}(T_0)} \big(\{T_n = T_0\} \cap \mathcal E_{n, e}\big),
\qquad
\mathcal E_{n, e} := \{\text{no arrow crosses } e \text{ after } \tau_n\},
```
has probability $`1`$, where $`E_{\mathrm{core}}(T_0)`$ is the set of core–core edges of $`T_0`$.

---

## 2. Exact Rate Cancellation and Infinite-Horizon Beta Factorization
Fix $`(n, T_0, e)`$ with $`e = (u, v) \in E_{\mathrm{core}}(T_0)`$. Deleting $`e`$ splits $`T_0`$ into two planted subtrees $`A_0`$ (of ordinary vertex size $`m_A \ge 2`$, rooted at $`u`$ with one exterior slot) and $`B_0`$ (of size $`m_B \ge 2`$, rooted at $`v`$ with one exterior slot), with $`m_A + m_B = n`$.

Run independent local evolutions on $`A`$ and $`B`$ after $`\tau_n`$ until the first attempted crossing of $`e`$ (selecting the slot at $`u`$ or at $`v`$ kills that side). On side $`A`$, at any reachable state $`S`$ of ordinary size $`m`$, the sum of genuine internal transition rates is $`\lambda(S) = m - 1/D_u(S)`$, while the kill rate (choosing the slot at $`u`$) is $`\kappa(S) = 1/D_u(S)`$. **Their sum cancels identically:**
```math
\lambda(S) + \kappa(S) \;=\; \Big(m - \frac{1}{D_u(S)}\Big) + \frac{1}{D_u(S)} \;=\; m.
\qquad\text{(1)}
```
Equivalently, construct side $`A`$ by preassigning independent holding times $`H_m^A \sim \mathrm{Exp}(m)`$ ($`m \ge m_A`$) and independent arrow-choice uniforms (at surviving state $`S`$ of size $`m`$, choosing each internal arrow $`r`$ with probability $`r/m`$ and killing with probability $`\kappa(S)/m`$, writing a cemetery symbol after kill). Let $`E_M^A`$ be survival through step $`M`$ and $`\mathcal E_A = \bigcap_{M \ge m_A} E_M^A`$. For any bounded cylinder $`f`$ of the times $`(H_m^A)`$ and any bounded cylinder $`g`$ of the discrete arrow marks, once $`M`$ exceeds all inspected indices, independence of the preassigned sequences gives $`\mathbb E[f g \mathbf 1_{E_M^A}] = \mathbb E[f]\mathbb E[g \mathbf 1_{E_M^A}]`$. Sending $`M \to \infty`$ by dominated convergence and applying the monotone-class theorem yields: **on $`\mathcal E_A`$, the entire holding-time sequence $`(H_m^A)_{m \ge m_A}`$ retains its exact independent $`\prod_{m \ge m_A} \mathrm{Exp}(m)`$ product law and is completely independent of the entire discrete jump chain $`(A_m)_{m \ge m_A}`$.**

**Strict positivity $`\mathbb P(\mathcal E_A) > 0`$ for every $`A_0`$ of size $`m_A \ge 2`$:** If $`u`$ is bare, a real neighbour sends it a leaf with positive probability; finitely many leaf-to-root arrows bring its leaf stock to $`K \ge 3`$ (all with positive probability even with killing active). From stock $`K`$, the stopped reciprocal-stock supermartingale $`1/L_u`$ (`open/forest-age/inert` §2) bounds the total future kill probability by $`14/(3K) < 1`$ for $`K \ge 5`$, so $`\mathbb P(\mathcal E_A) > 0`$ (and likewise $`\mathbb P(\mathcal E_B) > 0`$).

**Existence and positivity of $`X_A, X_B \in (0, \infty)`$ a.s. on $`\mathcal E_{n, e} = \mathcal E_A \cap \mathcal E_B`$ (Claim 6 alone):**
1. Since $`e`$ is a core–core edge, $`C(T_{\tau_n + t}) = C(A_t) + C(B_t)`$ and $`|T_{\tau_n + t}| = |A_t| + |B_t|`$ for all $`t \ge 0`$ on $`\mathcal E_{n, e}`$. By Claim 6 and the Yule clock, $`e^{-\mu t}(C(A_t) + C(B_t))`$ converges a.s. to a positive finite limit on $`\mathcal E_{n, e}`$.
2. Holding all marks and all times except $`H_{m_A}^A`$ fixed, the post-first-jump core count on side $`A`$ is a deterministic function $`F(t - s)`$ of $`t`$ minus $`s = H_{m_A}^A`$. By Fubini, convergence of $`e^{-\mu t}(C(A_t) + C(B_t))`$ holds for a full-measure set of $`s > 0`$. Picking $`s_1 < s_2`$ there with $`d = s_2 - s_1 > 0`$ and subtracting gives $`H(t) - e^{-\mu d} H(t - d) \to \ell`$ for $`H(t) = e^{-\mu t} F(t)`$, with $`H`$ bounded. Setting $`\tilde H(t) = H(t) - \ell/(1 - e^{-\mu d})`$ gives $`\limsup_{t\to\infty}|\tilde H(t)| \le e^{-\mu d}\limsup_{t\to\infty}|\tilde H(t)|`$, hence $`H(t) \to \ell/(1 - e^{-\mu d})`$. Thus $`Z_A = \lim_{t\to\infty} e^{-\mu t} C(A_t)`$ and $`X_A = Z_A / W_A^\mu = \lim_{m\to\infty} C(A_m)/m^\mu`$ exist in $`[0, \infty)`$ a.s. on $`\mathcal E_A`$ (and symmetrically for $`B`$).
3. To show $`X_A > 0`$ a.s. on $`\mathcal E_A`$, glue two identical copies of $`A_0`$ at their roots into a reachable free tree $`T_{\mathrm{dbl}}`$. On its positive-probability no-cross event, the two side limits $`X_A, X_A'`$ are i.i.d.; if $`\mathbb P(X_A = 0 \mid \mathcal E_A) = p > 0`$, then $`X = 0`$ on $`T_{\mathrm{dbl}}`$ with probability $`\ge \mathbb P(\mathcal E_{\mathrm{dbl}}) p^2 > 0`$, contradicting Claim 6 ($`X > 0`$ a.s. from any reachable core seed).

**Beta mixture & absolute continuity:** On $`\{T_n = T_0\} \cap \mathcal E_{n, e}`$, the exact product $`\mathrm{Exp}(m)`$ laws give $`W_A \sim \mathrm{Gamma}(m_A, 1)`$ and $`W_B \sim \mathrm{Gamma}(m_B, 1)`$ independent of each other and of $`(X_A, X_B)`$, so $`U := W_A / (W_A + W_B) \sim \mathrm{Beta}(m_A, m_B)`$ is independent of $`(X_A, X_B)`$, and
```math
X \;=\; \frac{Z_A + Z_B}{(W_A + W_B)^\mu} \;=\; U^\mu X_A + (1 - U)^\mu X_B.
\qquad\text{(2)}
```
For any fixed $`x_A, x_B > 0`$ and $`\mu \in (0, 1)`$, $`g(u) = u^\mu x_A + (1 - u)^\mu x_B`$ is strictly concave on $`(0, 1)`$ with a single interior critical point $`u_* \in (0, 1)`$ and $`g'(u) \ne 0`$ on $`(0, u_*) \cup (u_*, 1)`$. Exhausting $`(0, u_*)`$ and $`(u_*, 1)`$ by compact intervals on which $`g^{-1}`$ is Lipschitz shows $`\mathbb P(g(U) \in N) = 0`$ for every Lebesgue-null Borel set $`N \subset (0, \infty)`$. Integrating over $`(X_A, X_B)`$ gives $`\mathbb P(X \in N, T_n = T_0, \mathcal E_{n, e}) = 0`$. Summing over the countable cover of Section 1 yields $`\mathbb P(X \in N) = 0`$ (and Claim 6 gives $`\mathbb P(X = 0) = 0`$). This proves Part 1. $`\square`$

---

## 3. Conditional Full Support and Everywhere-Positive Lower Semicontinuous Density
Using the inert-slot mean coefficient $`\Phi`$ and uniform one-leaf gradient $`L_0`$ (`open/forest-age/inert` §§1–3, or Track 3's grafted re-rooting bound), $`Y_m = \frac{\Gamma(m)}{\Gamma(m+\mu)}\Phi(S_m)`$ for the unconditioned inert embedded chain satisfies $`\mathbb E[Y_{m+1} - Y_m \mid S_m] = O(m^{-1-\mu})`$ and centered jumps $`O(m^{-\mu})`$, giving $`\mathbb E[(X_i - Y_m)^2 \mid S_m] \le V m^{1-2\mu}`$ uniformly. For any $`x > 0`$ and open neighbourhood $`I \ni x`$: grow root leaves internally to large size $`N`$ (forcing $`Y_N < x`$), then elongate one arm until $`Y_m`$ first crosses $`x`$ (overshoot $`O(N^{-\mu})`$) while keeping root leaf stock $`\ge N - O(1)`$; for $`N`$ large enough, the conditional failure probabilities for $`X_i \notin I`$ and for a future slot-kill are each $`< 1/4`$, so $`\mathbb P(X_i \in I \mid \mathcal E_i) > 0`$. Thus each one-slot surviving limit has full support $`[0, \infty)`$.

On the positive-probability $`P_4`$ event (size $`n = 4`$, unique core–core edge $`e`$), $`X_A, X_B`$ are independent with full support $`[0, \infty)`$. Given any $`x > 0`$, choosing $`a = b \in (2^{\mu - 1} x, x)`$ places $`x`$ strictly inside $`(a, 2^{1-\mu} a)`$, so for $`(x_A, x_B)`$ in a small neighbourhood of $`(a, b)`$ and $`y`$ near $`x`$, the interior inverse branch has density $`\beta_{2,2}(u)/|g'(u)|`$ uniformly bounded below by a positive constant. Integrating over the positive-mass neighbourhood of $`(a, b)`$ yields a positive constant lower bound on $`p_X`$ locally around every $`x > 0`$, hence on every compact $`K \subset (0, \infty)`$.

Ordering all one-cut candidates $`(n, T_0, e)`$ by start size $`n`$ (placing the labelled $`P_4`$ candidates first at $`n = 4`$) and disjointifying via "current survives and every predecessor fails" preserves the exact product holding-time law (every predecessor test depends only on the history through $`n`$ and internal marked arrows of the current components, never on relative component timing; see `REPORT.md` §3). Extending each open-interval inverse density branch by $`0`$ outside its open image gives a nonnegative lower semicontinuous function for each fixed $`(x_A, x_B)`$; integrating and summing over the disjoint countable partition preserves lower semicontinuity by Fatou and produces an everywhere-positive, lower semicontinuous version $`p_X : (0, \infty) \to (0, \infty]`$ finite a.e., bounded below on every compact $`K \subset (0, \infty)`$. $`\square`$

---

## 4. Multi-Cut Dirichlet Regularity and Exact Open Gap
Cutting $`k`$ core–core edges around a common hub yields a disjoint clock-preserving Dirichlet mixture $`X = \sum_{i=1}^{k+1} U_i^\mu X_i`$ whose characteristic function for any localized piece decays as $`|\xi|^{-r/2}`$ (with $`r = \lfloor k/2 \rfloor`$ pairs), giving $`C^s`$ densities on each localized piece for arbitrarily large $`s`$. However, upgrading this to $`C^0 / C^k / C^\infty`$ or local upper boundedness of the full density $`p_X`$ requires summability of the piece-wise density/derivative bounds across the infinite countable partition, which remains open.

---

## 5. Verification and Replay
Replay command (from the repository root):
```sh
python3 -B round2/proofs/neon2/limit-density/check_density.py
```
Recorded output:
```text
Exact: 873 states by size 2–7 [1, 2, 6, 24, 120, 720]; 5912 path Laplace identities
Controls: 4/4 caught (including the leaf-edge false-core convention)
P4: seed 6142026, 6551/50000 accepted, both sides reach 128; mu=0.56575
Beta(2,2) KS D=0.00953, 1% DKW threshold=0.02011
Pearson r(U, X_A[128])=-0.01604, r(U, X_B[128])=+0.01289
```
See `round2/proofs/neon2/limit-density/REPORT.md` for the complete independent GPT-6.2 Sol Ultrafast audit.
