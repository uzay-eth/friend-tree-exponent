# Track 3: Unexploited Structure in Claim 14, Explicit Stationary Fringe Law, and Shorter Real-Variable Proofs of Claims 6 and 7

**Pinned claim commit:** `ea6a093c247d6b9bfd38dc27346bfe3c54bcf012`.  
**Review and dependencies:** GPT-family proof audits and scoped Claude review of the variance repair, gradient and summability. See [Claim 3 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-3.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-3.md).

**Author:** Neon-2 lead (Gemini 4 Neon).
**Cross-family cold audit:** GPT-6.2 Sol Ultrafast (`verify_track3.py`); this revision incorporates all repairs from both audit rounds (exact tightness constant, bounded generator domain, non-circular root-independence via $`s > 1`$, positive-probability reachability for arbitrary-seed growth, grafted re-rooting bound for inert trees, exact nonlocal coboundary $`Q = -\Psi^0(J) R`$, Cesàro-coboundary identification of $`d_k`$ bypassing marked/inert amplitude ratios, discrete $`L^2`$ + Doob maximal convergence and consistent $`\Psi^0(J)Y_n`$ full-support scale, and strict $`n > m`$ Abel equality in Theorem 3).

**Honest provenance & explicit Astra dependencies (contour-free Section 3 revision):** This work inherits Astra's working copy and explicitly lists every imported input:
- **Exponent existence and bracket $`1/2 < \ell < \mu < 1`$ ($`\ell = 565658/10^6`$):** Real-variable existence of $`\mu`$ is `proofs/existence.md` Theorem 1 (Astra); the numerical bracket is supplied directly by Neon-2 Claim 6 (`round2/proofs/neon2/quintet-cert/PROOF.md`, commit `b04cc50`, standalone zero-floating-point unoriented 5-vertex certificate proving $`0.565702177913 \le \mu \le 0.565748407438`$ with zero Astra numerical/asymptotic dependencies and zero contour integrals; alternatively by Astra's physical bracket using `renewal.md` §5's real constant-factor bounds).
- **Astra Claim 14 contraction only** (`round2/influence/zero/pair/README.md`, finite cap-48 row certificate, zero contour integrals): for $`\Re s \ge \ell := 565658/10^6`$, $`\|\Delta V_s C\|_{\mathrm{TV}} \le q := 0.8600511260 < 1`$ on $`\ell^\infty(\mathcal S)`$. The critical-point identities $`s M_s(J)(1 - U_s(J)) = 1`$ ($`s > \mu`$), $`U_\mu(J) = 1`$, and $`-U'_\mu(J) = 1/(\mu \Psi^0(J)) > 0`$ are **derived directly in Proposition 0 without citing `influence/zero/critical.md` or `spatial-limit.md`**.
- **Astra Claim 5 & Claim 12** (`round2/fringe/finite-kernel.md` §§1–2, `round2/finite/gauge/physical-trigger.md` §§1, 4, `round2/standalone/claim-6.md` §§2, 4 line 149 — purely positive semigroup & real-variable bounds, zero contour integrals): the 19-label positive residual weight $`W`$, the critical crossing $`F(\mu) = 1`$ and untriggered semigroup bound $`\hat R_t W \le C_0 e^{\rho' t} W`$ ($`\rho' < \mu`$), the physical first-moment identity with marked initial nonfounders, and the uniform founder-immigration bound $`I_{T, u} W \le K_I C(T) e^{q_0 u}`$ ($`q_0 = 9/25 < \mu`$).
- **Combinatorial eligibility & inert separation** (`proofs/existence.md` Lemma 2.1, `round2/open/forest-age/inert/README.md` §2 lines 68–83 — elementary tree counting and stopped reciprocal-stock supermartingale, zero contour integrals): at least $`(C(T) + 5)/4`$ eligible cores for $`C(T) \ge 2`$, and separation probability $`\le 14/(3k)`$ for a planted virgin.
- **Astra Claim 4 / Claim 10** (`round2/operator.md`, `round2/fringe/spatial-limit.md`): cited **ONLY** in Theorem 1(2) to identify our explicit measure $`\mathcal M_\mu`$ with Astra's named spatial marginal $`\pi_{\mathrm{fringe}}`$. Section 3 (our purely real-variable proof of Claim 6) does **NOT** cite `spatial-limit.md`, `almost-sure.md`, `asymptotics.md`, `open/limit-law/README.md`, `influence/zero/critical.md`, or Round 1's `proofs/almost-sure-exponent.md`: it uses **zero complex contour integrals, zero Fourier inversion, and no spectral gap $`\sigma < \mu`$**.
- **Section 4 (Theorem 3, summability Claim 7):** depends only on the post-promotion supermartingale $`D_v(n)/n`$ (zero Astra dependencies).

---

## 1. Explicit Contractive Representation of the Stationary Fringe Law $`\pi`$ and Exact Coboundaries

Let $`\mathcal S`$ be the countable set of isomorphism classes of finite unmarked planted trees $`S`$ with at least one core (root has an inert exterior parent slot counted in both $`D_o`$ and $`c_o`$), and $`J \in \mathcal S`$ the virgin (one root, one real leaf, one inert parent slot). Write $`L^0`$ for the intrinsic generator, $`V_s = (s - L^0)^{-1}`$ ($`\Re s > 0`$) on $`\ell^\infty(\mathcal S)`$ (operator-norm holomorphic on $`\{\Re s > 0\}`$), $`b(S) = L_o/D_o`$, $`(C f)(S) = D_o^{-1} \sum_{Y \text{ core child}} f(Y)`$, $`(\Delta f)(S) = f(S^+) - f(S)`$, and $`(T f)(S) = \sum_{v \in K(S)} f(S_v)`$. Let $`\mathcal D(L^0) := \{f \in \ell^\infty(\mathcal S) : L^0 f \in \ell^\infty(\mathcal S)\}`$ denote the bounded pointwise generator domain (equivalently, the range of $`V_s`$ on $`\ell^\infty(\mathcal S)`$, by the first-jump equation and stopped Dynkin).

### Theorem 1 (Explicit Stationary Fringe Measure & Holomorphic Continuation $`\mathcal M_s`$)
For every $`\psi \in \ell^\infty(\mathcal S)`$ and every $`s \in \mathbb C`$ with $`\Re s \ge \ell`$, the resolvent equation $`W = V_s(\psi + C \Delta W)`$ (equivalently $`(s - L^0) W = \psi + C \Delta W`$ with $`W \in \mathcal D(L^0)`$) has the unique bounded solution
```math
H_{\psi, s} := (I - \Delta V_s C)^{-1} \Delta V_s \psi, \qquad
W_{\psi, s} := V_s\big(\psi + C H_{\psi, s}\big),
```
satisfying $`\|H_{\psi, s}\|_\infty \le \frac{2}{\Re s (1 - q)}\|\psi\|_\infty`$ and $`\|W_{\psi, s}\|_\infty \le \frac{1}{\Re s}\big(1 + \frac{2}{\Re s (1 - q)}\big)\|\psi\|_\infty`$.
Define the bounded linear functional $`\mathcal M_s : \ell^\infty(\mathcal S) \to \mathbb C`$ on $`\{\Re s \ge \ell\}`$ (holomorphic in the open half-plane $`\{\Re s > \ell\}`$, continuous on the closed half-plane, and analytically extendable to an open neighbourhood of every boundary point on $`\{\Re s = \ell\}`$ by local invertibility of $`I - \Delta V_s C`$) by
```math
\mathcal M_s(\psi) \;:=\; s \, W_{\psi, s}(J) \;=\; s \Big( V_s \big( I + C(I - \Delta V_s C)^{-1} \Delta V_s \big) \psi \Big)(J).
\qquad\text{(1)}
```
Then:
1. **(Discounted occupation identity for $`s > \mu`$):** For every real $`s > \mu`$ and every $`\psi \in \ell^\infty(\mathcal S)`$, with $`A_\psi(S) := (T \psi)(S) = \sum_{v \in K(S)} \psi(S_v)`$ and $`M_s(J) = \int_0^\infty e^{-st} \mathbb E^0_J[C(S_t)]\,dt`$,
   $`\displaystyle \mathcal M_s(\psi) \;=\; \frac{\int_0^\infty e^{-st} \mathbb E^0_J[A_\psi(S_t)]\,dt}{\int_0^\infty e^{-st} \mathbb E^0_J[C(S_t)]\,dt}. \qquad\text{(2)}`$
2. **(Explicit stationary downward-fringe law at $`s = \mu`$):** As $`s \downarrow \mu`$, $`\mathcal M_\mu`$ is a countably additive probability measure on $`\mathcal S`$ with **no unknown Perron eigenvector**, coinciding with the downward-fringe marginal $`\pi_{\mathrm{fringe}}`$ of the unique spatial eigenmeasure $`\pi`$ of Claim 4 / Claim 10 (and on every $`\psi \in \ell^\infty(\mathcal S)`$, $`\mathcal M_\mu(\psi) = \sum_{S \in \mathcal S} \mathcal M_\mu(\mathbf 1_{\{S\}})\,\psi(S)`$).
3. **(Exact evaluations & bounded-domain generator invariance):** For all $`\Re s \ge \ell`$ and every $`f \in \mathcal D(L^0)`$,
   $`\displaystyle \mathcal M_s(1) = 1, \qquad \mathcal M_s(b) = s U_s(J), \qquad \mathcal M_\mu(U_\mu) = -\mu U'_\mu(J) = \frac{1}{\Psi^0(J)}, \qquad\text{(3)}`$
   $`\displaystyle \mathcal M_s\big(L^0 f + C \Delta f + b f(J) - s f\big) \;=\; s\big(U_s(J) - 1\big) f(J). \qquad\text{(4)}`$

**Proof.**
Uniqueness and norm bounds for $`H_{\psi,s}, W_{\psi,s}`$ follow from $`\|\Delta V_s C\|_{\mathrm{TV}} \le q < 1`$, $`\|C\|_{\infty\to\infty} \le 1`$, $`\|V_s\|_{\infty\to\infty} \le 1/\Re s`$, and $`\|\Delta\|_{\infty\to\infty} \le 2`$.

For any $`w \in \ell^\infty(\mathcal S)`$ and any finite $`S \in \mathcal S`$, summing $`L^0 w(S_v)`$ over all cores $`v \in K(S)`$ accounts for every transition of $`L^0`$ on $`S`$ except: (i) parent-to-child arrows $`p \to v`$ (rate $`1/D_p`$), which add $`\Delta w(S_v)`$ and re-group at $`p`$ into $`(C \Delta w)(S_p)`$; and (ii) promotions at $`v`$ (rate $`b(S_v)`$), which create a new virgin summand $`w(J)`$. Arrows to the exterior slot are suppressed. Hence pointwise on every finite $`S \in \mathcal S`$,
```math
L^0(T w) = T\big(L^0 w + C \Delta w + w(J) b\big).
\qquad\text{(5)}
```
Applying (5) to $`w = W_{\psi,s} \in \mathcal D(L^0)`$ and using $`(s - L^0) W_{\psi,s} = \psi + C \Delta W_{\psi,s}`$ gives $`(s - L^0)(T W_{\psi,s}) = A_\psi - W_{\psi,s}(J) A_b`$ on every finite $`S \in \mathcal S`$.
For real $`s > \mu`$, $`|T W_{\psi,s}| \le \|W_{\psi,s}\|_\infty C`$ and $`\mathbb E^0_J[C(S_t)] = O(e^{\mu t})`$ (by Proposition 0(3) below), so stopped Dynkin on $`[0, t]`$ (using finite-horizon Yule polynomial size domination) followed by $`t \to \infty`$ gives
```math
W_{\psi,s}(J) = (T W_{\psi,s})(J) = \int_0^\infty e^{-st} \mathbb E^0_J[A_\psi(S_t)]\,dt - W_{\psi,s}(J) \int_0^\infty e^{-st} \mathbb E^0_J[A_b(S_t)]\,dt.
```
Since $`A_b = L^0 C`$ and $`C(J) = 1`$, integration by parts gives $`\int_0^\infty e^{-st} \mathbb E^0_J[A_b(S_t)]\,dt = s M_s(J) - 1`$. Collecting the $`W_{\psi,s}(J)`$ terms yields $`s M_s(J) W_{\psi,s}(J) = \int_0^\infty e^{-st} \mathbb E^0_J[A_\psi(S_t)]\,dt`$ (with no division by $`W_{\psi,s}(J)`$). Dividing by $`M_s(J) > 0`$ proves (2).

From (2), $`\mathcal M_s`$ is a probability measure on $`\mathcal S`$ for every real $`s > \mu`$. Every core born at time $`\tau`$ starts with downward ordinary vertex size $`2`$, and its total real growth rate at size $`n`$ is at most $`n + 1/2`$ (the single incoming arrow from its core parent has rate $`\le 1/2`$). By comparison with the pure-birth chain of rate $`n + 1/2`$ started at $`2`$ and Jensen's inequality (`other-marginals.md` §2), after elapsed age $`a`$ we have $`\mathbb E[(|S_v(\tau + a)|)^{1/2} \mid \mathcal F_\tau] \le C_0 e^{a/2}`$ with $`C_0 = \sqrt{5/2}`$. Compensating core promotions (predictable intensity $`A_b`$) and using Tonelli gives, for every real $`s > \mu > 1/2`$,
```math
\int_0^\infty e^{-st} \mathbb E^0_J\!\Big[\sum_{v \in K(S_t)} |S_v(t)|^{1/2}\Big]\,dt
\;\le\; \frac{C_0}{s - 1/2}\!\left(1 + \int_0^\infty e^{-su} \mathbb E^0_J[A_b(S_u)]\,du\right)
\;=\; \frac{C_0 s M_s(J)}{s - 1/2}.
```
Testing (2) on the bounded truncations $`|S_o|^{1/2} \wedge M`$ and sending $`M \to \infty`$ by monotone convergence yields the exact uniform moment bound
```math
\mathcal M_s(|S_o|^{1/2}) \;\le\; \frac{C_0 s}{s - 1/2} \;\le\; \frac{C_0 \mu}{\mu - 1/2} \qquad (s > \mu),
```
since $`s/(s - 1/2)`$ is decreasing for $`s > 1/2`$. As $`s \downarrow \mu`$, $`\mathcal M_s(\psi) \to \mathcal M_\mu(\psi)`$ for every $`\psi \in \ell^\infty(\mathcal S)`$ (in particular $`\mathcal M_\mu(\psi) \ge 0`$ for $`\psi \ge 0`$). For any $`M \ge 2`$, the set $`\mathcal S_M := \{S \in \mathcal S : |S| \le M\}`$ is finite and $`\mathcal M_s(\mathbf 1_{\mathcal S \setminus \mathcal S_M}) \le \frac{C_0 \mu}{(\mu - 1/2)\sqrt M}`$ uniformly for $`s > \mu`$, hence also at $`s = \mu`$. Splitting any $`\psi \in \ell^\infty(\mathcal S)`$ into $`\psi \mathbf 1_{\mathcal S_M} + \psi \mathbf 1_{\mathcal S \setminus \mathcal S_M}`$ and letting $`M \to \infty`$ proves that $`\mathcal M_\mu`$ is a countably additive probability measure on $`\mathcal S`$ with $`\mathcal M_\mu(\psi) = \sum_{S \in \mathcal S} \mathcal M_\mu(\mathbf 1_{\{S\}})\psi(S)`$ for all $`\psi \in \ell^\infty(\mathcal S)`$.

For (3): $`W_{1,s} \equiv 1/s \in \mathcal D(L^0)`$ solves $`(s - L^0)(1/s) = 1 + C\Delta(1/s)`$, so $`\mathcal M_s(1) = 1`$; $`W_{b,s} = U_s`$ gives $`\mathcal M_s(b) = s U_s(J)`$. Differentiating the bounded-resolvent equation $`U_s = V_s(b + C \Delta U_s)`$ in $`\ell^\infty(\mathcal S)`$ at $`s = \mu`$ (using $`V'_s = -V_s^2`$) yields $`-U'_\mu = V_\mu(U_\mu + C \Delta(-U'_\mu))`$, so uniqueness gives $`W_{U_\mu, \mu} = -U'_\mu`$ and $`\mathcal M_\mu(U_\mu) = -\mu U'_\mu(J) = 1/\Psi^0(J)`$ by Proposition 0(3) (with no citation to `critical.md`).

For (4): if $`f \in \mathcal D(L^0)`$, then $`\psi := (s - L^0)f - C \Delta f \in \ell^\infty(\mathcal S)`$ and $`f = V_s(\psi + C \Delta f)`$, so uniqueness gives $`W_{\psi, s} = f`$ and hence $`\mathcal M_s((s - L^0)f - C \Delta f) = s f(J)`$. Subtracting $`f(J) \mathcal M_s(b) = s U_s(J) f(J)`$ proves (4).

Finally, $`\mathcal M_\mu = \pi_{\mathrm{fringe}}`$ follows directly by an Abelian ratio from Claim 10's inert-seed mean convergence (`spatial-limit.md` §1 applied to the inert virgin cell with silent exterior): whenever $`e^{-\mu t}\mathbb E^0_J[A_\psi(S_t)] \to \Psi^0(J)\pi_{\mathrm{fringe}}(\psi)`$ holds for a bounded score $`\psi`$, splitting the Laplace integrals in (2) at a fixed large time and using $`M_s(J) \to \infty`$ as $`s \downarrow \mu`$ (Proposition 0(3)) yields $`\mathcal M_\mu(\psi) = \pi_{\mathrm{fringe}}(\psi)`$. (Section 3 below uses ONLY $`\mathcal M_\mu(\psi_k)`$ directly and never invokes `spatial-limit.md`.) $`\square`$

---

### Proposition 0 (Contour-Free Real Key Renewal, Uniform Trigger Amplitude, and Claim 14 Critical Identities)
Using **only** Claim 5 (`finite-kernel.md` §§1–2), Claim 12 (`physical-trigger.md` §§1, 4), Neon-2 Claim 6 ($`\ell < \mu < 1`$), Claim 14's finite-table contraction $`\|\Delta V_s C\|_{\mathrm{TV}} \le q < 1`$ on $`\{\Re s \ge \ell\}`$ (`pair/README.md`), `proofs/existence.md` Lemma 2.1, and `standalone/claim-6.md` §§2, 4 line 149 (zero complex contour integrals, zero Fourier inversion, no spectral gap $`\sigma < \mu`$):
1. **Real 1D Key Renewal Theorem & uniform physical-tree mean limit:**
   In the 19-label positive residual process of Claim 5 with Claim 12's degree-two trigger split $`R + B_0 J_{\mathrm{mark}} = \hat R + h_{\mathrm{trig}} \otimes \delta_j`$ ($`j = J_{1/2}`$, $`h_{\mathrm{trig}}(S) = b_S \mathbf 1_{\{D_o = 2\}}`$, untriggered positive semigroup bound $`\hat R_t W \le C_0 e^{\rho' t} W`$ with $`\rho' < \mu`$), set $`Z_k(t) := e^{-\mu t}(P_t f_k)(j)`$ and $`z_k(t) := e^{-\mu t}(\hat R_t f_k)(j)`$ for $`f_k(S) = \mathbf 1_{\{D_o \ge k\}}`$ ($`k \ge 2`$; $`0 \le f_k \le 1 < W`$).
   - *Direct Riemann integrability:* By the finite-successor first-event integral and local weight bound of `standalone/claim-6.md` §2, $`t \mapsto (\hat R_t f_k)(j)`$ is continuous on $`[0, \infty)`$ ($`z_2(0) = 1`$), and $`0 \le z_k(t) \le C_0 W(j) e^{-(\mu - \rho') t}`$, so $`z_k`$ is directly Riemann integrable.
   - *Feller's real Key Renewal Theorem:* The positive-semigroup first-trigger Duhamel identity along lineages gives $`Z_k(t) = z_k(t) + \int_0^t Z_k(t - u)\,e^{-\mu u} p(u)\,du`$, where $`e^{-\mu u} p(u)\,du`$ is a continuous nonlattice probability density on $`(0, \infty)`$ ($`\int_0^\infty e^{-\mu u} p(u)\,du = F(\mu) = 1`$, $`p(u) \ge \frac{1}{2} e^{-2u} > 0`$) with finite positive mean $`-F'(\mu) \in (0, \infty)`$. By Feller's real Key Renewal Theorem (Vol. II, Ch. XI.1), $`Z_k(t) \to \tilde d_k := (-F'(\mu))^{-1}\int_0^\infty z_k(u)\,du`$ ($`\tilde d_2 > 0`$) and $`\sup_{t \ge 0} Z_k(t) < \infty`$.
   - *Uniform marked-seed and physical-tree convergence:* For any 19-marked planted seed $`S`$, let $`q_S(u) := e^{-\mu u}(\hat R_u h_{\mathrm{trig}})(S) \le C_0 e^{-(\mu - \rho') u} W(S)`$ and $`H_{\mathrm{trig}}(S) := \int_0^\infty q_S(u)\,du \le \frac{C_0}{\mu - \rho'} W(S)`$. Splitting the first-trigger Duhamel integral $`e^{-\mu t}(P_t f_k)(S) = e^{-\mu t}(\hat R_t f_k)(S) + \int_0^t q_S(u) Z_k(t - u)\,du`$ at $`u = t/2`$ yields $`|e^{-\mu t}(P_t f_k)(S) - H_{\mathrm{trig}}(S)\tilde d_k| \le \varepsilon_k(t) W(S)`$ with $`\varepsilon_k(t) \to 0`$ and $`\sup_{t \ge 0}\varepsilon_k(t) < \infty`$. Applying the physical first-moment identity (`standalone/claim-6.md` §2 eq. (1)) from any core founder $`o \in K(T)`$ with compatible initial nonfounder marks (e.g. $`1/2`$), defining
     $`\displaystyle A_{\mathrm{trig}}(T) \;:=\; \sum_{v \in K(T) \setminus \{o\}} H_{\mathrm{trig}}(\tilde S_v) + \int_0^\infty e^{-\mu u} I_{T, u}(H_{\mathrm{trig}})\,du,`$
     and inserting $`\sum_{v \ne o} W(\tilde S_v) \le K_1 C(T)`$ and $`I_{T, u} W \le K_I C(T) e^{q_0 u}`$ ($`q_0 = 9/25 < \mu`$; `standalone/claim-6.md` §4 line 149, splitting the immigration convolution at $`u = h/2`$) gives
     $`\displaystyle \big|e^{-\mu h}\mathbb E_T[N_{\ge k}(T_h)] - \tilde d_k A_{\mathrm{trig}}(T)\big| \;\le\; \tilde\varepsilon_k(h) C(T), \qquad \tilde\varepsilon_k(h) \to 0 \text{ as } h \to \infty, \qquad\text{(12)}`$
     uniformly over all $`T \in \mathcal T`$.
2. **Direct contour-free lower bound on $`A_{\mathrm{trig}}(T)`$:**
   Under $`\hat R`$, degree-two *direct offspring* are removed as triggers (`h_{\mathrm{trig}} = b \mathbf 1_{\{D_o = 2\}}`), while intrinsic promotions and non-degree-two direct births remain and the only death is the fixed root mark $`a \le 1/2`$. For every marked newborn $`J_{1/i}`$ ($`2 \le i \le 20`$), before its first jump (intrinsic change rate $`3/2`$, death rate $`\le 1/2`$, total rate $`\le 2`$) it emits $`h_{\mathrm{trig}}(J_{1/i}) = 1/2`$, so $`H_{\mathrm{trig}}(J_{1/i}) \ge \int_0^\infty e^{-(\mu + 2) u} \frac{1}{2}\,du = \frac{1}{2\mu + 4} > \frac{1}{2\mu + 5} > 0`$. By the two-step local survival script of `standalone/claim-6.md` §4 / `renewal.md` §5 applied directly inside $`\hat R`$:
   - For any eligible nonfounder $`v`$ with $`c_v \le 2, L_v \ge 1`$, require the original particle to survive through time $`1`$ (cost $`\ge e^{-1/2}`$) and the root's first ring in $`[0, 1]`$ to choose a leaf (probability $`\ge (1 - e^{-1})/3`$, since $`c_v \le 2`$ stays fixed and $`L_v`$ cannot decrease before the root's first ring). Compensating that promotion at its pre-promotion degree $`D`$: if $`D = 2`$, the promotion itself pays immediate trigger toll $`1`$ (`h_{\mathrm{trig}} = b`); if $`D > 2`$, it emits a retained direct-birth newborn $`J_{1/m(D)}`$ paying child trigger potential $`H_{\mathrm{trig}}(J_{1/m(D)}) \ge \frac{1}{2\mu + 5}`$. Both tolls are $`\ge \frac{1}{2\mu + 5}`$ and discounting over $`[0, 1]`$ costs at most $`e^{-\mu}`$, giving $`H_{\mathrm{trig}}(\tilde S_v) \ge c_1' > 0`$.
   - For a bare nonfounder $`v`$ ($`c_v = D_v = 2`$) with a bare core child $`w`$, silence $`v`$ and $`w`$'s other real neighbor and forbid death during $`[0, 1]`$, letting $`w`$'s first ring choose $`v`$ (probability $`1/2`$, since $`w`$ remains degree $`2`$ before its first ring and bare states have $`b = 0`$ so no competing trigger occurs), producing $`D_v \ge 3, c_v = 2, L_v \ge 1`$ and applying the first bound on $`[1, 2]`$.
   Charging bare eligible vertices without a bare child to their bare parent (at most two charges each, excluding at most one founder) via `existence.md` Lemma 2.1 supplies at least $`C(T)/12 - 1`$ nonfounders with $`H_{\mathrm{trig}}(\tilde S_v) \ge c_1 > 0`$, yielding
   $`\displaystyle A_{\mathrm{trig}}(T) \;\ge\; c_1\big(C(T)/12 - 1\big)_+ \qquad (\forall T \in \mathcal T), \qquad A_{\mathrm{trig}}(P_3) > 0`$
   (for $`P_3`$, the no-event interval has strictly positive founder direct-immigration rate to $`J_{1/2}`$).
3. **Contour-free derivation of Claim 14's critical identities via inert-founder first-moment identity:**
   Apply the exact first-moment identity (`standalone/claim-6.md` §2 eq. (1)) to the inert virgin $`J`$ by retaining the inert root $`o`$ explicitly as the zero-mark founder (with its slot-inclusive $`p_o, b_o, c_o`$, zero incoming arrow from its silent parent slot, and root score $`1`$ outside the marked sum):
   $`\displaystyle \mathbb E_J^0[C(S_t)] \;=\; 1 + \int_0^t I_u^J(P_{t - u} 1)\,du.`$
   The exact same weight estimate as `standalone/claim-6.md` §2 eq. (1) and §4 line 149 (with $`d = 3/20`$, $`q_0 = 9/25 < \mu`$, zero top-mark killing, and no incoming top arrow; every promoted child enters with weight $`W(J_{1/m}) \le 1/d`$ and every core child copy is weighted by its capped mark, with $`\mathbb E_J^0[Z(S_u)] \le e^{q_0 u} h(1) = e^{q_0 u}`$) yields $`I_u^J W \le 1/d + e^{q_0 u} \le K_I^0 e^{q_0 u}`$ ($`K_I^0 := 1/d + 1`$). Convolving with the uniform marked-seed KRT $`|e^{-\mu v}(P_v 1)(S) - H_{\mathrm{trig}}(S)\tilde d_2| \le \varepsilon_2(v) W(S)`$ of part (1) (splitting at $`u = t/2`$) gives
   $`\displaystyle e^{-\mu t}\mathbb E_J^0[C(S_t)] \;\longrightarrow\; D_J \;:=\; \tilde d_2 \int_0^\infty e^{-\mu u} I_u^J(H_{\mathrm{trig}})\,du \;\in\; (0, \infty),`$
   where $`D_J < \infty`$ since $`q_0 < \mu`$, and $`D_J > 0`$ because the initial inert virgin $`J`$ promotes directly to $`J_{1/2}`$ at rate $`1/2`$ before its first intrinsic jump (with $`H_{\mathrm{trig}}(J_{1/2}) \ge \frac{1}{2\mu + 5} > 0`$). In particular $`\mathbb E_J^0[C(S_t)] = O(e^{\mu t})`$ and the real Abelian limit gives $`\lim_{s \downarrow \mu} (s - \mu) M_s(J) = D_J \in (0, \infty)`$.

   Applying (5) to $`U_s = W_{b, s}`$ (holomorphic on $`\{\Re s \ge \ell\}`$ by `pair/README.md` alone) gives $`(s - L^0)(T U_s) = (1 - U_s(J)) L^0 C`$; stopped Dynkin on $`[0, t]`$ and $`t \to \infty`$ for real $`s > \mu`$ (using $`\mathbb E_J^0[C(S_t)] = O(e^{\mu t})`$ and integration by parts $`\int_0^\infty e^{-st}\mathbb E_J^0[L^0 C(S_t)]\,dt = s M_s(J) - 1`$) yields
   $`\displaystyle s M_s(J)\big(1 - U_s(J)\big) \;=\; 1 \qquad (\forall s > \mu),`$
   hence $`1 - U_s(J) = \frac{1}{s M_s(J)} > 0`$ for all $`s > \mu`$. Letting $`s \downarrow \mu`$ immediately gives
   $`\displaystyle U_\mu(J) \;=\; 1, \qquad -U'_\mu(J) \;=\; \lim_{s \downarrow \mu} \frac{1 - U_s(J)}{s - \mu} \;=\; \lim_{s \downarrow \mu} \frac{1}{s (s - \mu) M_s(J)} \;=\; \frac{1}{\mu D_J} \;>\; 0.`$
   Setting $`\Psi^0(J) := D_J \in (0, \infty)`$ gives $`\mathcal M_\mu(U_\mu) = -\mu U'_\mu(J) = 1/\Psi^0(J)`$. Finally, for every $`T \in \mathcal T`$, let $`\Phi_s(T; o)`$ be defined by (6)–(7) below (which uses only `pair/README.md`). By (12) at $`k = 2`$, $`m_t^{\mathrm{free}}(T) = O_T(e^{\mu t})`$, so stopped Dynkin on (8) holds for all real $`s > \mu`$ and gives $`\Phi_s(T; o) = (1 - U_s(J))(s M_s^{\mathrm{free}}(T) - C(T)) > 0`$. Taking the real Abelian limit as $`s \downarrow \mu`$ yields
   $`\displaystyle \Phi(T) \;=\; \lim_{s \downarrow \mu} \frac{1 - U_s(J)}{s - \mu}\cdot (s - \mu)\big(s M_s^{\mathrm{free}}(T) - C(T)\big) \;=\; \frac{\tilde d_2}{\Psi^0(J)}\,A_{\mathrm{trig}}(T) \;\ge\; \frac{\tilde d_2 c_1}{\Psi^0(J)}\big(C(T)/12 - 1\big)_+,`$
   establishing **both** $`\tilde d_2 A_{\mathrm{trig}}(T) = \Psi^0(J)\Phi(T)`$ **and** affine lower comparability $`C(T) \le K_4(1 + \Phi(T))`$ purely via real-variable renewal theory and real Abelian limits. $`\square`$

---

## 2. Exact Free-Tree Eigenfunction $`\Phi`$, One-Line Re-Rooting Gradient, and Free-Tree Coboundaries

Let $`\mathcal T`$ be the set of finite free (unrooted) trees $`T`$ with at least one core ($`C(T) \ge 1`$), evolving under the physical continuous-time friend-tree generator $`\mathcal L`$. For each $`T \in \mathcal T`$, write $`m_t^{\mathrm{free}}(T) := \mathbb E_T[C(T_t)]`$ and $`M_s^{\mathrm{free}}(T) := \int_0^\infty e^{-st} m_t^{\mathrm{free}}(T)\,dt`$.

### Theorem 2 (Exact Free-Tree Eigenfunction & One-Line Uniform Gradient)
For any $`T \in \mathcal T`$, any chosen core root $`o \in K(T)`$, and any $`s \in \mathbb C`$ with $`\Re s \ge \ell`$, define
```math
f_s(T, o) \;:=\; \int_0^\infty e^{-st} \mathbb E_{T, o}\!\left[ b_o(t) + \frac{1}{D_o(t)} \sum_{Y \text{ core child of } o \text{ at } t} H_s(Y) \right] dt,
\qquad\text{(6)}
```
```math
\Phi_s(T; o) \;:=\; f_s(T, o) + \sum_{v \in K(T) \setminus \{o\}} U_s(S_v),
\qquad\text{(7)}
```
where $`(T_t, o)`$ evolves under $`\mathcal L`$ tracking the core vertex $`o`$ (no inert slot), and each non-root core $`v \ne o`$ carries its downward planted subtree $`S_v \in \mathcal S`$ away from $`o`$ (retaining its parent slot toward $`o`$). Then:
1. $`\|f_s\|_\infty \le \frac{1 + \|H_s\|_\infty}{\Re s} < \infty`$ uniformly over all $`(T, o)`$.
2. Pointwise on every $`T \in \mathcal T`$ and every $`o \in K(T)`$,
   $`\displaystyle (s - \mathcal L) \Phi_s(T; o) \;=\; \big(1 - U_s(J)\big) (\mathcal L C)(T). \qquad\text{(8)}`$
3. **Root-independence & eigenfunction:** For every $`\Re s \ge \ell`$, $`\Phi_s(T; o)`$ is **independent of the choice of root $`o \in K(T)`$**. For every real $`s > \mu`$, $`\Phi_s(T) = (1 - U_s(J))\big(s M_s^{\mathrm{free}}(T) - C(T)\big) > 0`$. By real continuity as $`s \downarrow \mu`$, $`\Phi(T) := \Phi_\mu(T)`$ satisfies $`\Phi(T) \ge 0`$ and $`\mathcal L \Phi(T) = \mu \Phi(T)`$ pointwise on $`\mathcal T`$, with $`\Phi(T) \le b_0 C(T)`$ ($`b_0 = \|f_\mu\|_\infty + \|U_\mu\|_\infty`$) and affine lower comparability $`C(T) \le K_4(1 + \Phi(T))`$ (and hence strict linear comparability $`a_0 C(T) \le \Phi(T) \le b_0 C(T)`$ for $`a_0 > 0`$ on all $`T \in \mathcal T`$).
4. **(One-line re-rooting uniform one-leaf gradient — free AND inert trees):**
   - For every free tree $`T \in \mathcal T`$ and every $`T^+`$ obtained by attaching one new leaf to ANY vertex $`x \in T`$ (core or leaf),
     $`\displaystyle |\Phi_s(T^+) - \Phi_s(T)| \;\le\; 2\|f_s\|_\infty + |U_s(J)|, \qquad\text{(9)}`$
     and at $`s = \mu`$, $`|\Phi(T^+) - \Phi(T)| \le L := 2\|f_\mu\|_\infty + 1 \le 94`$.
   - **Grafted extension to inert planted trees:** For $`\Phi^0(S) := (T U_\mu)(S) = \sum_{v \in K(S)} U_\mu(S_v)`$ on $`\mathcal S`$, we have $`L^0 \Phi^0 = \mu \Phi^0`$, $`\Phi^0 \ge 0`$, and for EVERY one-leaf attachment $`S^+`$ anywhere in $`S`$,
     $`\displaystyle |\Phi^0(S^+) - \Phi^0(S)| \;\le\; L \;\le\; 94. \qquad\text{(9')}`$

**Proof.**
(1) A tracked core remains a core; at every $`t \ge 0`$, $`0 \le b_o(t) \le 1`$ and $`\big|\frac{1}{D_o(t)} \sum_Y H_s(Y)\big| \le \frac{c_o(t)}{D_o(t)}\|H_s\|_\infty \le \|H_s\|_\infty`$, giving the uniform bound.
(2) By the first-jump equation (stopped Dynkin) for the bounded resolvent (6), $`(s - \mathcal L) f_s(T, o) = b_o + D_o^{-1}\sum_{Y \text{ core child of } o} H_s(Y)`$. Applying $`\mathcal L`$ to $`\sum_{v \ne o} U_s(S_v)`$ produces the exact internal $`L^0 U_s(S_v)`$ terms, parent-to-child transfers $`(C H_s)(S_p)`$ for $`p \ne o`$, root-to-child transfers $`D_o^{-1}\sum_Y H_s(Y)`$ (cancelling the root term in $`(s - \mathcal L) f_s`$), and promotion newborn terms $`b_w U_s(J)`$ over all $`w \in K(T)`$. Since $`L^0 U_s + C H_s = s U_s - b`$, this gives (8).
(3) To prove root-independence **without circularity**, start on real $`s > 1`$: Yule domination $`m_t^{\mathrm{free}}(T) \le \mathbb E_T[|T_t|] = |T| e^t`$ justifies stopped Dynkin on (8) and integration by parts, yielding
```math
\Phi_s(T; o) = \big(1 - U_s(J)\big)\big(s M_s^{\mathrm{free}}(T) - C(T)\big) \qquad (s > 1),
```
whose right-hand side is independent of $`o \in K(T)`$ and strictly positive (as $`1 - U_s(J) > 0`$ by Proposition 0(3) and every finite tree with a core has a leaf adjacent to a core, so $`s M_s^{\mathrm{free}}(T) > C(T)`$). Since $`f_s(T, o)`$ and $`U_s(S_v)`$ are holomorphic on $`\{\Re s > \ell\}`$ and continuous up to $`\{\Re s = \ell\}`$, the identity theorem extends root-independence $`\Phi_s(T; o_1) = \Phi_s(T; o_2)`$ to the entire half-plane $`\{\Re s \ge \ell\}`$. Furthermore, Proposition 0(1) eq. (12) gives $`m_t^{\mathrm{free}}(T) = O_T(e^{\mu t})`$ for every fixed $`T \in \mathcal T`$. Hence the real occupation formula holds for all real $`s > \mu`$, so passing to the limit along real $`s \downarrow \mu`$ gives $`\Phi(T) = \Phi_\mu(T) \ge 0`$ and $`\Phi(T) = \frac{1}{\Psi^0(J)} \lim_{s \downarrow \mu} (s - \mu) M_s^{\mathrm{free}}(T)`$.

Upper comparability $`\Phi(T) \le (\|f_\mu\|_\infty + \|U_\mu\|_\infty) C(T)`$ is immediate from (7). Affine lower comparability $`C(T) \le K_4(1 + \Phi(T))`$ ($`K_4 \ge 1`$) follows immediately from Proposition 0(3): $`\Phi(T) = \frac{\tilde d_2}{\Psi^0(J)} A_{\mathrm{trig}}(T) \ge \frac{\tilde d_2 c_1}{\Psi^0(J)}(C(T)/12 - 1)_+`$ (with $`\Phi \ge 0`$ covering $`C(T) \le 12`$), using **zero contour integrals**. To upgrade this to a uniform positive lower bound $`\Phi(T) \ge c_{\min} > 0`$, set $`N := \lceil 2 K_4 \rceil \ge 2`$. From ANY seed $`T \in \mathcal T`$, pick a core $`x`$ with $`c_x \le 1`$ (the sole core if $`C = 1`$, or a leaf of the core tree if $`C \ge 2`$) and one of its real leaves $`\ell`$. On $`[0, 1]`$, silence $`x`$ and require $`\ell`$ to ring at least $`2N`$ times (since $`x`$ is silent, no leaf of $`x`$ can promote, giving $`L_x(1) \ge 2N + 1`$ and $`c_x(1) \le 1`$). On $`[1, 2]`$, require $`x`$ to ring at least $`N`$ times and each of its first $`N`$ rings to choose a leaf: before each of those $`N`$ rings, $`L_x \ge N + 2`$ and $`c_x \le N`$, so each conditional leaf-choice probability is at least $`(N+2)/(2N+2) > 1/2`$. This event has uniform positive probability $`c_* := e^{-1}\mathbb P(\mathrm{Pois}(1) \ge 2N)\,\mathbb P(\mathrm{Pois}(1) \ge N)\,2^{-N} > 0`$ independent of $`T`$, and forces $`C(T_2) \ge N \ge 2 K_4`$, where $`\Phi(T_2) \ge C(T_2)/K_4 - 1 \ge 1`$. By harmonicity $`\Phi(T) = e^{-2\mu}\mathbb E_T[\Phi(T_2)]`$ (justified by stopped Dynkin and finite-horizon Yule domination, or by Step 3.1's martingale (11) which does not use lower comparability), we obtain $`\Phi(T) \ge e^{-2\mu} c_* =: c_{\min} > 0`$ uniformly over ALL $`T \in \mathcal T`$. Combining $`\Phi(T) \ge c_{\min}`$ with $`C(T) \le K_4(1 + \Phi(T))`$ yields strict linear comparability $`a_0 C(T) \le \Phi(T) \le b_0 C(T)`$ with $`a_0 := \frac{c_{\min}}{K_4(1 + c_{\min})} > 0`$ on ALL $`T \in \mathcal T`$.

(4) **Re-rooting proof of (9):**
- If $`x \in K(T)`$ is a core, choose $`o = x`$ in both $`T`$ and $`T^+`$. Every non-root core $`v \ne x`$ has the identical downward subtree $`S_v`$ in $`T`$ and $`T^+`$, so $`\Phi_s(T^+) - \Phi_s(T) = f_s(T^+, x) - f_s(T, x)`$ (modulus $`\le 2\|f_s\|_\infty`$).
- If $`x`$ is a leaf, let $`p \in K(T)`$ be its unique core neighbour and choose $`o = p`$ in both $`T`$ and $`T^+`$. Attaching a leaf to $`x`$ promotes $`x`$ to a new non-root virgin $`J`$ while keeping $`D_p`$ fixed ($`c_p \to c_p + 1, L_p \to L_p - 1`$) and leaving every other non-root subtree $`S_v`$ ($`v \notin \{p, x\}`$) unchanged. Hence $`\Phi_s(T^+) - \Phi_s(T) = f_s(T^+, p) - f_s(T, p) + U_s(J)`$ (modulus $`\le 2\|f_s\|_\infty + |U_s(J)|`$).
Using $`\ell = 565658/10^6`$ and $`q = 0.8600511260`$, $`\|H_\mu\|_\infty < 25.264265`$, $`\|f_\mu\|_\infty < 46.431350`$, and $`L = 2\|f_\mu\|_\infty + 1 < 93.862700 < 94`$.

**Grafted proof of (9') for inert planted trees:** Given any $`S \in \mathcal S`$ and integer $`M \ge 1`$, form the free tree $`T_M(S) \in \mathcal T`$ by replacing the inert parent slot of $`S`$ with a real adjacent auxiliary core $`a`$ bearing $`M`$ additional real leaves, and choose $`o = a`$ as the root in (7). Because $`S`$ is the sole core child subtree of $`a`$, (7) gives the exact identity
```math
\Phi\big(T_M(S)\big) \;=\; f_\mu\big(T_M(S), a\big) + \Phi^0(S).
```
Under $`\mathcal L`$ started from $`T_M(S)`$, the auxiliary vertex $`a`$ satisfies $`D_a(t) \ge M + 1`$ for all $`t \ge 0`$, and $`c_a(t) \le 1 + N_a(t)`$ where $`N_a(t) \sim \mathrm{Pois}(t)`$ counts the rings of $`a`$'s own rate-one clock (only rings of $`a`$ can promote its adjacent leaves). Hence at every $`t \ge 0`$, since $`b_a(t) = 1 - c_a(t)/D_a(t)`$,
```math
\mathbb E_{T_M(S), a}\!\left|\left(b_a(t) + \frac{1}{D_a(t)}\sum_{Y \text{ core child of } a} H_\mu(Y)\right) - 1\right|
\;\le\; \frac{(1 + \|H_\mu\|_\infty)(1 + t)}{M + 1}.
```
Integrating against $`e^{-\mu t}\,dt`$ gives $`|f_\mu(T_M(S), a) - 1/\mu| \le \frac{(1 + \|H_\mu\|_\infty)(1/\mu + 1/\mu^2)}{M + 1} \to 0`$ as $`M \to \infty`$, uniformly in $`S`$. Applying (9) to $`T_M(S^+)`$ and $`T_M(S)`$ at the same $`M`$ and sending $`M \to \infty`$ proves $`|\Phi^0(S^+) - \Phi^0(S)| \le L`$. Finally, (5) with $`U_\mu(J) = 1`$ gives $`L^0 \Phi^0 = \mu \Phi^0`$. For nonnegativity of $`\Phi^0(S)`$, every finite $`S \in \mathcal S`$ is reachable from $`J`$ at time $`1`$ with positive probability $`p_S := \mathbb P_J^0(S_1 = S) > 0`$ (build $`S`$ from $`J`$ by finite positive-rate leaf attachments and no further jumps before time $`1`$), so by the Markov property $`p_S \mathbb E_S^0[C(S_t)] \le \mathbb E_J^0[C(S_{t+1})] = O(e^{\mu t})`$. Hence stopped Dynkin on $`(s - L^0)(T U_s) = (1 - U_s(J)) L^0 C`$ holds for all real $`s > \mu`$ (with vanishing endpoint as $`t \to \infty`$), yielding $`(T U_s)(S) \ge 0`$ for all $`s > \mu`$ and therefore $`\Phi^0(S) = \lim_{s \downarrow \mu} (T U_s)(S) \ge 0`$. $`\square`$

### Corollary 2.1 (Exact Free-Tree Coboundaries with One-Line Re-Rooting Bound)
1. **(Local degree scores):** For any bounded degree score $`\psi(D, c)`$ (such as $`\psi_k(D, c) = \mathbf 1_{\{D \ge k\}}`$; recall that for every non-root core $`v \ne o`$, its parent slot in $`S_v`$ restores both its true free-tree degree $`D_v`$ and true core degree $`c_v`$), write $`A_\psi(T) := \sum_{v \in K(T)} \psi(D_v, c_v)`$. For $`\Re s \ge \ell`$, set $`r_{\psi, s}(S) := W_{\psi, s}(S) - W_{\psi, s}(J)`$ (so $`r_{\psi, s}(J) = 0`$) and
   $`\displaystyle f_{\psi, s}(T, o) := \int_0^\infty e^{-st} \mathbb E_{T, o}\!\left[\psi(D_o(t), c_o(t)) - \mathcal M_s(\psi) + \frac{1}{D_o(t)}\sum_{Y \text{ core child of } o} H_{\psi, s}(Y)\right] dt,`$
   $`\displaystyle R_{\psi, s}(T) := f_{\psi, s}(T, o) + \sum_{v \in K(T) \setminus \{o\}} r_{\psi, s}(S_v).`$
   Then $`R_{\psi, s}(T)`$ is independent of $`o \in K(T)`$ (proved first for real $`s > 1`$ via Yule domination, then by analyticity on $`\Re s \ge \ell`$), satisfies the exact coboundary identity
   $`\displaystyle (s - \mathcal L) R_{\psi, s}(T) \;=\; A_\psi(T) - \mathcal M_s(\psi) C(T) \qquad\text{(10)}`$
   on every $`T \in \mathcal T`$, obeys $`|R_{\psi, s}(T)| \le (\|f_{\psi, s}\|_\infty + \|r_{\psi, s}\|_\infty) C(T)`$, and has the one-line uniform jump bound $`|R_{\psi, s}(T^+) - R_{\psi, s}(T)| \le 2 \|f_{\psi, s}\|_\infty`$ (with no virgin term since $`r_{\psi, s}(J) = 0`$).
2. **(Exact nonlocal coboundary for $`C - \Psi^0(J)\Phi`$):** Fix $`\psi^*(S) := U_\mu(S) - 1 \in \ell^\infty(\mathcal S)`$ and tracked-root score $`\psi^*_{\mathrm{root}}(T, o) := f_\mu(T, o) - 1`$ (kept fixed while varying $`s`$). Define $`R_{\psi^*, s}(T; o)`$ by replacing $`\psi(D_o(t), c_o(t))`$ in $`f_{\psi, s}`$ with $`\psi^*_{\mathrm{root}}(T_t, o)`$. Then by (7), the exact oriented source in $`(s - \mathcal L) R_{\psi^*, s}(T; o)`$ is
   $`\displaystyle \big(f_\mu(T, o) - 1\big) + \sum_{v \ne o} \big(U_\mu(S_v) - 1\big) - \mathcal M_s(\psi^*) C(T) \;=\; \Phi(T) - \big(1 + \mathcal M_s(\psi^*)\big) C(T),`$
   which is **identically independent of $`o`$ for every $`s`$**. By Yule domination for $`s > 1`$ and analyticity on $`\Re s \ge \ell`$, $`R_{\psi^*, s}(T)`$ is root-independent with uniform jump bound $`\le 2\|f_{\psi^*, s}\|_\infty`$. At $`s = \mu`$, Theorem 1(3) gives $`1 + \mathcal M_\mu(U_\mu - 1) = \mathcal M_\mu(U_\mu) = 1/\Psi^0(J)`$, so $`Q(T) := -\Psi^0(J) R_{\psi^*, \mu}(T)`$ satisfies **with zero remainder**:
   $`\displaystyle (\mu - \mathcal L) Q(T) \;=\; C(T) - \Psi^0(J) \Phi(T), \qquad |Q(T)| \le K_Q C(T), \qquad |Q(T^+) - Q(T)| \le L_Q. \qquad\text{(10')}`$
   Integrating $`e^{-\mu t}`$ against (10') and (10): the discounted jump martingales associated to $`Q`$ and $`R_{\psi_k, \mu}`$ have uniform one-leaf jump bounds ($`L_Q`$ and $`2\|f_{\psi_k, \mu}\|_\infty`$) and hence are $`L^2`$-bounded and converge a.s. by the exact same bracket isometry as (11); both $`|Q(T)| \le K_Q C(T)`$ and $`|R_{\psi_k, \mu}(T)| \le K_k C(T)`$ ensure $`e^{-\mu T}|Q(T_T)|/T \to 0`$ and $`e^{-\mu T}|R_{\psi_k, \mu}(T_T)|/T \to 0`$ a.s. (by $`C \le K_4(1 + \Phi)`$ and $`M_t \to M_\infty`$). Applying (10') first gives the a.s. Cesàro limit $`\frac{1}{T}\int_0^T e^{-\mu t} C(T_t)\,dt \to \Psi^0(J) M_\infty`$, and then applying (10) to $`\psi_k`$ gives $`\frac{1}{T}\int_0^T e^{-\mu t} N_{\ge k}(T_t)\,dt \to \mathcal M_\mu(\psi_k)\Psi^0(J) M_\infty`$ with zero renewal theory (which Theorem 1(2) separately identifies with $`\pi(D_o \ge k)\Psi^0(J) M_\infty`$).

---

## 3. Shorter, Purely Real-Variable Proof of Astra Claim 6 (Zero Contour Integrals)

Start from $`T_2`$ (equivalently from its deterministic first unlabelled successor $`T_3`$; start the continuous-time clock at $`T_3`$ so $`n e^{-\tau_n} \to W_{\mathrm{size}} \sim \mathrm{Gamma}(3, 1)`$). Using only the contour-free inputs listed at the top of this note, we prove that there exists a square-integrable random variable $`X \in (0, \infty)`$ a.s. with full support $`[0, \infty)`$ such that simultaneously for all integers $`k \ge 2`$,
```math
\frac{C_n}{n^\mu} \to X \quad (\text{a.s. and in } L^2), \qquad \frac{N_{\ge k}(n)}{C_n} \to \mathcal M_\mu(\psi_k), \qquad \frac{N_{=k}(n)}{C_n} \to \mathcal M_\mu(\mathbf 1_{\{D_o = k\}}) \quad \text{a.s.},
```
with deterministic mean asymptotics $`\mathbb E[N_{\ge k}(n)] / n^\mu \to c_k := \mathcal M_\mu(\psi_k)\mathbb E[X]`$ and $`c_k/c_2 = \mathcal M_\mu(\psi_k)`$ (which Theorem 1(2) separately identifies with Astra's named spatial marginal $`\pi(D_o \ge k)`$).

### Step 3.1 (Three-line $`L^2`$-bounded martingale convergence — continuous AND discrete)
By Theorem 2, $`M_t := e^{-\mu t} \Phi(T_t)`$ is a nonnegative càdlàg martingale from any $`T_0 \in \mathcal T`$. Every jump attaches one leaf (total rate $`|T_t|`$) and changes $`M_t`$ by at most $`L e^{-\mu t}`$ by (9). Stopping at finite sizes, applying jump-martingale isometry with $`\mathbb E_{T_0}[|T_u|] = |T_0| e^u`$, and sending the stop to infinity gives
```math
\mathbb E_{T_0}\big[(M_t - M_0)^2\big] \;\le\; L^2 \int_0^t e^{-2\mu u} \mathbb E_{T_0}[|T_u|]\,du \;=\; L^2 |T_0| \int_0^t e^{(1 - 2\mu) u}\,du \;\le\; \frac{L^2 |T_0|}{2\mu - 1}
\qquad\text{(11)}
```
since $`2\mu > 1`$. Hence $`M_t \to M_\infty`$ a.s. and in $`L^2`$, with $`\mathbb E_{T_0}[M_\infty] = \Phi(T_0) > 0`$.

For the discrete embedded chain $`(T_n)_{n \ge 3}`$, set $`a_n := \Gamma(n + \mu)/\Gamma(n)`$ ($`a_{n+1}/a_n = 1 + \mu/n`$, $`a_n \sim n^\mu`$). Since the jump rate at size $`n`$ is $`n`$, $`\mathcal L \Phi = \mu \Phi`$ gives $`\mathbb E[\Phi(T_{n+1}) \mid T_n] = (1 + \mu/n)\Phi(T_n)`$, so $`Y_n := \Phi(T_n)/a_n`$ is a nonnegative discrete martingale. Because $`\Phi(T_n) \le b_0 C_n \le b_0 n`$ and $`|\Phi(T_{n+1}) - \Phi(T_n)| \le L`$, its increment obeys $`|Y_{n+1} - Y_n| = a_{n+1}^{-1}|\Phi(T_{n+1}) - \Phi(T_n) - (\mu/n)\Phi(T_n)| \le (L + \mu b_0)/a_{n+1}`$. Hence the conditional increment variances sum uniformly over ALL size-$`n`$ seeds $`T_n`$:
```math
\sum_{j=n}^\infty \mathbb E\big[(Y_{j+1} - Y_j)^2 \mid T_n\big] \;\le\; (L + \mu b_0)^2 \sum_{j=n}^\infty a_{j+1}^{-2} \;=\; O(n^{1 - 2\mu}) \;<\; \infty.
\qquad\text{(11')}
```
Therefore $`(Y_n)_{n \ge 3}`$ is $`L^2`$-bounded, converges a.s. and in $`L^2`$ to a limit $`Y_\infty \in L^2`$, and satisfies the uniform conditional tail-variance bound $`\mathop{\mathrm{Var}}\nolimits (Y_\infty \mid T_n) = \mathbb E[(Y_\infty - Y_n)^2 \mid T_n] \le (L + \mu b_0)^2 \sum_{j=n}^\infty a_{j+1}^{-2} = O(n^{1-2\mu})`$.

### Step 3.2 (Self-contained 10-line Lévy 0–1 proof of $`M_\infty > 0`$ a.s. — no `almost-sure-exponent.md` needed)
For the inert virgin $`J`$, by Theorem 2(4) (eq. (9')), $`M_t^0 := e^{-\mu t}\Phi^0(S_t)`$ is a nonnegative martingale with jump size $`\le L e^{-\mu t}`$ and total jump rate $`\le |S_t|`$ (so $`\mathbb E_J[|S_u|] \le 2 e^u`$). By the identical 3-line isometry (11), $`M_t^0 \to Z_J^0`$ a.s. and in $`L^2`$, with $`\mathbb E_J[Z_J^0] = \Phi^0(J) = U_\mu(J) = 1`$, hence $`p := \mathbb P_J(Z_J^0 > 0) > 0`$.

From ANY finite seed $`T \in \mathcal T`$: pick a core $`x`$ with $`c_x \le 1`$ (the sole core if $`C = 1`$, or a leaf of the core tree if $`C \ge 2`$) and one of its real leaves $`\ell`$.
- On $`[0, 1]`$, require $`x`$ not to ring and $`\ell`$ to ring at least $`k`$ times (probability $`e^{-1}\mathbb P(\mathrm{Pois}(1) \ge k)`$); no leaf of $`x`$ can promote while $`x`$ is silent, so $`L_x \ge k + 1`$ and $`c_x \le 1`$.
- On $`[1, 2]`$, require $`x`$ to ring at least once and its first choice to be a leaf (conditional probability $`\ge (1 - e^{-1})/2`$). This promotes a new virgin child $`v`$ at time $`\tau \le 2`$ while $`L_x \ge k`$. Thus this planting event has probability at least $`a_k := e^{-1}\mathbb P(\mathrm{Pois}(1) \ge k)(1 - e^{-1})/2 > 0`$, independent of $`T`$.
- Couple $`v`$'s downward branch from $`\tau`$ to an inert virgin $`J`$ until the first arrow from $`x`$ into $`v`$. By the stopped reciprocal-stock supermartingale $`1/L_x`$ (`forest-age/inert` §2), the total separation probability over $`[\tau, \infty)`$ is at most $`14/(3k)`$.
Fix $`k`$ with $`14/(3k) \le p/2`$. By a union bound, with probability $`\ge p/2`$ conditional on the planting event, there is no separation and $`Z_J^0 > 0`$. On that event, $`C(T_{\tau + t}) \ge C(S_t^{(v)})`$ for all $`t \ge 0`$, so by $`\Phi^0 \le \|U_\mu\|_\infty C`$ and $`a_0 C \le \Phi`$,
```math
M_\infty \;=\; \lim_{t\to\infty} e^{-\mu(\tau + t)} \Phi(T_{\tau + t})
\;\ge\; \frac{a_0}{\|U_\mu\|_\infty} e^{-\mu \tau} \lim_{t\to\infty} e^{-\mu t} \Phi^0(S_t^{(v)})
\;=\; \frac{a_0}{\|U_\mu\|_\infty} e^{-\mu \tau} Z_J^0 \;>\; 0.
```
Thus $`\mathbb P_T(M_\infty > 0) \ge \delta := a_k p/2 > 0`$ uniformly over ALL $`T \in \mathcal T`$. Since $`\{M_\infty > 0\}`$ is invariant under restarting at integer time $`j`$ (factor $`e^{-\mu j} > 0`$), the Markov property gives $`\mathbb P(M_\infty > 0 \mid \mathcal F_j) \ge \delta`$ for all $`j \ge 1`$; Lévy's 0–1 law forces $`M_\infty > 0`$ almost surely.

### Step 3.3 (Purely Real-Variable Key Renewal + Elementary Fixed-Lag Coupling Variance for $`N_{\ge k}(T_t)`$)
1. **Real-variable uniform mean limit (already established in Proposition 0(1) & (3)):**
   By (12) and the Abelian identification $`\tilde d_2 A_{\mathrm{trig}}(T) = \Psi^0(J)\Phi(T)`$ of Proposition 0(3), we have for every $`k \ge 2`$ and every $`T \in \mathcal T`$:
   $`\displaystyle \big|e^{-\mu h}\mathbb E_T[N_{\ge k}(T_h)] - d_k \Phi(T)\big| \;\le\; \tilde\varepsilon_k(h) C(T), \qquad \tilde\varepsilon_k(h) \to 0 \text{ as } h \to \infty, \qquad\text{(12')}`$
   where $`d_k := (\tilde d_k / \tilde d_2)\Psi^0(J)`$ (with $`d_2 = \Psi^0(J)`$), proved purely via Feller's real Key Renewal Theorem and real Abelian limits (zero complex contour integrals, no spectral gap $`\sigma < \mu`$).
2. **Self-contained elementary fixed-lag one-leaf coupling variance, almost-sure limit, and Cesàro identification of $`d_k`$:**
   Because we take $`t \to \infty`$ at a **fixed** lag $`h > 0`$ before sending $`h \to \infty`$, we need NO uniform-in-lag gradient and NO citation to `open/limit-law/README.md`. Given any $`T \in \mathcal T`$ and any one-leaf addition $`T^+`$ (adding leaf $`x^+`$ to parent $`p`$), couple $`(T_u, T_u^+)_{u \ge 0}`$ on the same vertex label space by letting $`A(u) := (V(T_u) \triangle V(T_u^+)) \cup \{v \in V(T_u) \cap V(T_u^+) : N_{T_u}(v) \ne N_{T_u^+}(v)\}`$ be the asymmetric discrepancy set ($`|A(0)| = 2`$, namely $`\{x^+, p\}`$).
   - *Clean rings preserve $`A(u)`$:* Every common vertex $`w \notin A(u)`$ has identical neighbor sets $`N_{T_u}(w) = N_{T_u^+}(w)`$. Couple its rate-1 clock and uniform neighbor choice $`y \in N(w)`$ identically in both trees, attaching the same fresh leaf label $`\ell`$ to $`y`$ in both trees. Then $`N_{T_u}(\ell) = \{y\} = N_{T_u^+}(\ell)`$ (so $`\ell \notin A`$), and adding the same label $`\ell`$ to both $`N_{T_u}(y)`$ and $`N_{T_u^+}(y)`$ preserves whether $`y \in A`$. Thus clean rings never enlarge $`A(u)`$.
   - *Discrepancy rings:* Each vertex $`v \in A(u)`$ rings at rate $`1`$ (one shared clock if $`v \in V(T_u) \cap V(T_u^+)`$, one private clock if $`v \in V(T_u) \triangle V(T_u^+)`$), attaching at most one leaf to one neighbor in $`T_u`$ (adding at most $`2`$ vertices to $`A`$) and at most one leaf to one neighbor in $`T_u^+`$ (adding at most $`2`$ vertices to $`A`$), so $`|A(u)|`$ increases by at most $`+4`$ at rate at most $`|A(u)|`$.
   Stopping at bounded sizes and applying Yule domination and Grönwall gives $`\mathbb E[|A(u)|] \le 2 e^{4u}`$. Since every vertex outside $`A(u)`$ has identical degree in $`T_u`$ and $`T_u^+`$, we have $`|N_{\ge k}(T_u^+) - N_{\ge k}(T_u)| \le |A(u)|`$ pointwise, hence uniformly for ALL $`T \in \mathcal T`$, ALL $`k \ge 2`$, and ALL $`u \ge 0`$:
   $`\displaystyle \big|\mathbb E_{T^+}[N_{\ge k}(T_u)] - \mathbb E_T[N_{\ge k}(T_u)]\big| \;\le\; 2 e^{4u}. \qquad\text{(13a)}`$
   Applying the backward square-integrable martingale jump isometry to $`r \mapsto \mathbb E[N_{\ge k}(T_t) \mid \mathcal F_r]`$ on $`[t - h, t]`$ (for any fixed lag $`h > 0`$ and $`t > h`$, total jump rate $`|T_r|`$ with each one-leaf jump bounded by $`2 e^{4(t - r)}`$) gives
   $`\displaystyle \mathop{\mathrm{Var}}\nolimits \big(e^{-\mu t} N_{\ge k}(T_t) \mid \mathcal F_{t-h}\big) \;\le\; 4 e^{-2\mu t} \int_{t-h}^t e^{8(t-r)} \mathbb E[|T_r| \mid \mathcal F_{t-h}]\,dr \;=\; 4 e^{-2\mu t} |T_{t-h}| \int_0^h e^{8(h-v) + v}\,dv \;\le\; \frac{4}{7} e^{8h} e^{-2\mu t} |T_{t-h}|. \qquad\text{(13)}`$
   Taking expectations ($`\mathbb E[|T_{t-h}|] = 3 e^{t-h}`$) bounds the second moment of $`e^{-\mu t}(N_{\ge k}(T_t) - \mathbb E[N_{\ge k}(T_t) \mid \mathcal F_{t-h}])`$ by $`\frac{12}{7} e^{7h} e^{-(2\mu - 1) t}`$, which at fixed $`h`$ is summable along any mesh $`t_j = j/m > h`$ because $`2\mu > 1`$. Chebyshev and Borel–Cantelli give $`e^{-\mu t_j}(N_{\ge k}(T_{t_j}) - \mathbb E[N_{\ge k}(T_{t_j}) \mid \mathcal F_{t_j - h}]) \to 0`$ a.s. Combining with (12') and $`C \le K_4(1 + \Phi)`$ gives $`\limsup_{j\to\infty}|e^{-\mu t_j} N_{\ge k}(T_{t_j}) - d_k M_\infty| \le \tilde\varepsilon_k(h) K_4 M_\infty`$ a.s.; intersecting over countable $`k, m, h`$ and letting $`h \to \infty`$ kills the error. Monotonicity of $`N_{\ge k}(T_t)`$ bridges between mesh points with factors $`e^{\pm \mu/m}`$; sending $`m \to \infty`$ gives $`e^{-\mu t} N_{\ge k}(T_t) \to d_k M_\infty`$ at ALL real times $`t`$ simultaneously for all $`k \ge 2`$.
   - *Purely real-variable Cesàro identification of $`d_k`$ (zero `spatial-limit.md` needed):* Since $`e^{-\mu t} N_{\ge k}(T_t) \to d_k M_\infty`$ a.s., its time Cesàro average $`\frac{1}{T}\int_0^T e^{-\mu t} N_{\ge k}(T_t)\,dt`$ converges a.s. to $`d_k M_\infty`$. Matching this with the independently proved Cesàro limit $`\mathcal M_\mu(\psi_k)\Psi^0(J) M_\infty`$ of Corollary 2.1 and cancelling $`M_\infty > 0`$ (Step 3.2) immediately identifies $`d_k = \mathcal M_\mu(\psi_k)\Psi^0(J)`$.
   - *Discrete Yule epochs, $`L^2`$ convergence, and fixed-size mean ratios:* Passing to the Yule epochs $`\tau_n`$ yields $`C_n / n^\mu \to X := \Psi^0(J) W_{\mathrm{size}}^{-\mu} M_\infty = \Psi^0(J) Y_\infty \in (0, \infty)`$ a.s., $`N_{\ge k}(n)/C_n \to \mathcal M_\mu(\psi_k)`$ a.s., and $`N_{=k}(n)/C_n \to \mathcal M_\mu(\mathbf 1_{\{D_o = k\}})`$ by subtraction. Finally, linear comparability gives $`0 \le N_{\ge k}(n)/n^\mu \le C_n / n^\mu \le a_0^{-1}(a_n / n^\mu) Y_n \le C' \sup_{m \ge 3} Y_m`$, where $`\sup_{m \ge 3} Y_m \in L^2`$ by Doob's $`L^2`$ maximal inequality applied to the $`L^2`$-bounded martingale $`(Y_n)`$ of (11'); Dominated Convergence therefore proves $`C_n / n^\mu \to X`$ in $`L^2`$ as well as a.s. and simultaneously yields $`\mathbb E[N_{\ge k}(n)] / n^\mu \to c_k := \mathcal M_\mu(\psi_k)\mathbb E[X]`$ (giving $`c_k/c_2 = \mathcal M_\mu(\psi_k)`$ and $`p_k = \mathcal M_\mu(\mathbf 1_{\{D_o = k\}})`$ without citing Claim 3's Tauberian argument or `spatial-limit.md`).

### Step 3.4 (Full support of $`X`$ on $`[0, \infty)`$ via discrete leaf-relocation)
For any open interval $`(a, b) \subset (0, \infty)`$, pick a target $`x \in (a, b)`$ and margin $`\eta > 0`$ with $`(x - 2\eta, x + 2\eta) \subset (a, b)`$. Working consistently in the scale $`\Psi^0(J) Y_n = \mathbb E[X \mid T_n]`$: at fixed large size $`n \ge 4`$, the star $`S_n`$ has $`\Psi^0(J) Y_n(S_n) \le \Psi^0(J) b_0 / a_n = O(n^{-\mu}) < x - \eta`$ and the path $`P_n`$ has $`\Psi^0(J) Y_n(P_n) \ge \Psi^0(J) a_0(n - 2)/a_n = \Omega(n^{1-\mu}) > x + \eta`$. Transform $`P_n`$ into $`S_n`$ step by step at fixed size $`n \ge 4`$ by relocating one leaf at distance $`\ge 2`$ from a fixed center to that center; each relocation deletes a leaf (leaving a tree of size $`n - 1 \ge 3`$ with at least one core, inside $`\mathcal T`$) and reattaches it, so by (9) it changes $`\Phi`$ by at most $`2L`$ and changes $`\Psi^0(J) Y_n`$ by at most $`2\Psi^0(J) L / a_n < \eta`$ for $`n`$ large. Thus some unlabelled tree shape $`T_*`$ of size $`n`$ satisfies $`|\Psi^0(J) Y_n(T_*) - x| < \eta`$. Every unlabelled tree shape of size $`n`$ has positive probability to be reached from $`T_3`$, and since $`\Phi`$ is label- and root-independent with $`\mathbb E[X \mid T_n] = \Psi^0(J) Y_n`$, (11') and conditional Chebyshev give $`\mathbb P(|X - \Psi^0(J) Y_n(T_*)| \ge \eta \mid T_n = T_*) \le (\Psi^0(J))^2 O(n^{1 - 2\mu})/\eta^2 < 1/2`$ for $`n`$ large. Hence $`\mathbb P(X \in (a, b)) > 0`$, proving full support $`[0, \infty)`$. $`\square`$

---

## 4. Strengthening Astra Claim 7: Summability Criterion $`\sum C_n / n^2 < \infty`$ (No Power Bound Needed)

### Theorem 3 (Hub Mass Conservation Under $`\sum C_n / n^2 < \infty`$ a.s.)
Let $`Z_v = \lim_{n\to\infty} D_v(n)/n`$ (existence from the post-promotion supermartingale $`D_v(n)/n`$), and let $`p_n := \mathcal L(W_n \mid T_n)`$ denote the conditional law of the next attachment target $`W_n`$. Suppose only that
```math
\sum_{n=3}^\infty \frac{C_n}{n^2} \;<\; \infty \qquad \text{almost surely}
\qquad\text{(14)}
```
(for example, $`C_n = O(n / (\log n)^{1+\varepsilon})`$ a.s.; no expectation bound and no power bound $`C_n = O(n^\beta)`$, $`\beta < 1`$, is required). Then almost surely $`\sum_{v \ge 1} Z_v = 1`$, $`\|p_n - (Z_v)\|_{\mathrm{TV}} \to 0`$, and $`\mathcal L(D_{W_n}(n)/n \mid T_n) \Rightarrow \sum_v Z_v \delta_{Z_v}`$.

**Proof.**
First, (14) and monotonicity of $`C_n`$ imply $`C_n/n \to 0`$ a.s.: for every $`n \ge 6`$,
```math
\frac{C_{\lfloor n/2 \rfloor}}{2 n} \;\le\; C_{\lfloor n/2 \rfloor}\!\left(\frac{1}{\lfloor n/2 \rfloor} - \frac{1}{n + 1}\right) \;=\; C_{\lfloor n/2 \rfloor}\sum_{r = \lfloor n/2 \rfloor}^n \frac{1}{r(r+1)} \;\le\; \sum_{r = \lfloor n/2 \rfloor}^n \frac{C_r}{r(r+1)} \;\longrightarrow\; 0;
```
evaluating along even $`n = 2m`$ gives $`C_m / m \to 0`$ a.s.

For fixed $`K \ge 1`$, define the stopping time $`\tau_K := \inf\big\{r \ge 3 : \sum_{j=3}^r \frac{C_j}{j(j+1)} > K\big\}`$ ($`\inf\varnothing = \infty`$) and set $`B_r := C_{r \wedge \tau_K}`$. Then $`0 \le B_r / r \le 1`$ deterministically (as $`C_j \le j`$) and $`B_n / n \to 0`$ everywhere as $`n \to \infty`$ (on $`\{\tau_K < \infty\}`$ because $`B_n = C_{\tau_K}`$ is constant for $`n \ge \tau_K`$; on $`\{\tau_K = \infty\}`$ because $`\sum_{j=3}^\infty C_j/(j(j+1)) \le K`$ implies $`C_n/n \to 0`$ deterministically by the first display). Furthermore, for every $`m \ge 3`$,
```math
0 \;\le\; \sum_{r = m+1}^\infty \frac{B_r}{r(r+1)} \;\le\; K + 1
\qquad\text{(15)}
```
everywhere: on $`\{\tau_K = \infty\}`$ the whole sum is at most $`K`$; on $`\{\tau_K < \infty\}`$ the strictly pre-crossing part ($`r < \tau_K`$) is at most $`K`$ and the frozen tail ($`r \ge \tau_K`$) is at most $`\sum_{r = \tau_K}^\infty \frac{C_{\tau_K}}{r(r+1)} = \frac{C_{\tau_K}}{\tau_K} \le 1`$. Moreover $`\sum_{r=m+1}^\infty \frac{B_r}{r(r+1)} \to 0`$ pointwise everywhere as $`m \to \infty`$.

Fix $`m \ge 3`$ and construct $`X_n^{(m)}`$ ($`n \ge m`$) starting at $`X_m^{(m)} = 0`$: at step $`r - 1 \to r`$ ($`r > m`$), on the $`\mathcal F_{r-1}`$-measurable event $`\{r \le \tau_K\} = \{\tau_K > r - 1\}`$, update the normalized degrees $`D_v(r)/r`$ of already tracked cores promoted after $`m`$ and add $`2/r`$ for any vertex promoted at step $`r`$; freeze $`X^{(m)}`$ once $`r > \tau_K`$ (so $`X_n^{(m)} \equiv 0`$ if $`\tau_K \le m`$). Since $`D_v(r)/r`$ is a supermartingale after promotion and promotions add $`2(B_r - B_{r-1})/r`$ (including the crossing step $`r = \tau_K`$), exact summation by parts and monotonicity of $`B_r`$ ($`B_n/n = B_n\sum_{r=n}^\infty \frac{1}{r(r+1)} \le \sum_{r=n}^\infty \frac{B_r}{r(r+1)}`$) give for every **$`n > m \ge 3`$**,
```math
\mathbb E[X_n^{(m)}] \;\le\; 2\,\mathbb E\!\sum_{r=m+1}^n \frac{B_r - B_{r-1}}{r}
\;=\; 2\,\mathbb E\!\left[\frac{B_n}{n} - \frac{B_m}{m+1} + \sum_{r=m+1}^{n-1} \frac{B_r}{r(r+1)}\right]
\;\le\; 2\,\mathbb E\!\left[\sum_{r=m+1}^\infty \frac{B_r}{r(r+1)}\right],
```
while at $`n = m`$, the upper bound $`\mathbb E[X_m^{(m)}] = 0 \le 2\mathbb E\sum_{r=m+1}^\infty \frac{B_r}{r(r+1)}`$ is trivial.

Let $`q_n(v) := L_v(n)/n`$ on cores ($`0`$ otherwise), so $`\sum_v q_n(v) = (n - C_n)/n`$ and $`q_n(v) \to Z_v`$ for every fixed $`v`$ on $`\{C_n/n \to 0\}`$ (on promoted cores $`c_v(n) \le C_n = o(n)`$; on forever-leaves $`D_v(n) \equiv 1`$ so $`q_n(v) = 0 = Z_v`$). Every core with label $`v > m`$ was promoted after $`m`$ and satisfies $`L_v \le D_v`$, so on $`\{\tau_K = \infty\}`$ we have $`X_n^{(m)} \ge \sum_{v > m} q_n(v) \to 1 - \sum_{v \le m} Z_v \ge 0`$ as $`n \to \infty`$. By Fatou's lemma as $`n \to \infty`$,
```math
\mathbb E\!\left[\mathbf 1_{\{\tau_K = \infty\}}\Big(1 - \sum_{v \le m} Z_v\Big)\right] \;\le\; 2\,\mathbb E\!\left[\sum_{r=m+1}^\infty \frac{B_r}{r(r+1)}\right].
```
Sending $`m \to \infty`$ at fixed $`K`$ via bounded convergence using (15) makes the right-hand side tend to $`0`$, proving $`\sum_{v \ge 1} Z_v = 1`$ a.s. on $`\{\tau_K = \infty\}`$. Letting integer $`K \to \infty`$ exhausts probability $`1`$ by (14). Coordinatewise convergence $`q_n(v) \to Z_v`$ with $`\sum_v q_n(v) \to 1 = \sum_v Z_v`$ yields $`\|q_n - Z\|_1 \to 0`$, hence $`\|p_n - Z\|_1 \le \|p_n - q_n\|_1 + \|q_n - Z\|_1 \le 3 C_n / n + \|q_n - Z\|_1 \to 0`$ a.s., and weak convergence $`\mathcal L(D_{W_n}(n)/n \mid T_n) \Rightarrow \sum_v Z_v \delta_{Z_v}`$ follows by testing bounded continuous functions on $`[0, 1]`$, replacing $`p_n`$ by $`Z`$ in the sum, and applying dominated convergence against the probability vector $`(Z_v)`$. $`\square`$

*(Continuous-time expectation remark: under the stronger expectation integrability $`\int_0^\infty e^{-s}\mathbb E[C(s)]\,ds < \infty`$, compensating promotions of rate $`\le C(s)`$ and conditioning at each promotion time $`\tau_v = s`$ to bound $`\mathbb E[D_v(t) \mid \mathcal F_s] \le 2 e^{t-s}`$ gives $`\sup_{t \ge T} \mathbb E\big[\sum_{v : T < \tau_v \le t} e^{-t} D_v(t)\big] \le 2\int_T^\infty e^{-s}\mathbb E[C(s)]\,ds \to 0`$ as $`T \to \infty`$.)*

---

## 5. Replay
```sh
python3 -B round2/proofs/neon2/claim14-proofs/verify_track3.py
```
