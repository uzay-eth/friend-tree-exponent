CORRECT

# Cold review of `round2/degree-law.md`

Reviewer: Gemini 4 Carbon  
Target: `../degree-law.md` (`degree-law.md`) (using `../operator.md` (`operator.md`) and `../../proofs/existence.md` (`proofs/existence.md`))  
Checker: `operator_carbon.py` (`python3 round2/checks/operator_carbon.py`, Check 3C–3D and Planted Error 8)

## Verdict summary

The proposed theorem and proof in `../degree-law.md` (`degree-law.md`) are **CORRECT**. Given a probability eigenmeasure $`\pi`$ of the core first-moment semigroup $`(Q_t)_{t\ge0}`$ satisfying $`\pi Q_t = e^{\mu t}\pi`$ and $`\pi(e^{\theta c_o}) < \infty`$ for some $`\theta > 0`$ (both supplied by Section 2 of `../operator.md` (`operator.md`), with $`\mu \in (0,1]`$), all three components of the proof—the birth-time renewal / stable-age representation (2), the finiteness and strict positivity of $`A = \mathbb E_\beta[W^\mu]`$ yielding $`\pi(D_o \ge k) \sim A k^{-\mu}`$, and the degree-crossing flux identity yielding $`\pi(D_o = k) \sim \mu A k^{-1-\mu}`$—are complete and mathematically sound. The draft also accurately restricts its claim to the tilted eigenmeasure $`\pi`$ rather than fixed-size empirical degree frequencies.

---

## Detailed audit

### 1. Birth-time renewal and stable-age mixture (Lines 11–22, Equation (2))

- **Cell decomposition and predictable compensation.** In the infinite-volume driver construction of `../operator.md` (`operator.md`) §1, every member of $`H_o(t)`$ is either the persistent old root $`o`$ (contributing $`(S_t f)(G, o) = \mathbb E_{G, o}[f(G_t, o)]`$) or a core vertex born at a unique promotion stopping time $`s \in (0, t]`$ from a parent $`u \in H_o(s-)`$, after which its own root trajectory over the remaining $`t - s`$ time units is tracked by $`S_{t-s} f`$ (its own later children are counted at their own birth times in $`H_o`$). For any bounded nonnegative cylinder $`f`$, the spatial cutoff estimate of `../operator.md` (`operator.md`) eq. (3) (which holds with factor $`1`$ in place of Yule $`L^2`$ norms for the single-root operator $`S_r`$) implies joint continuity of $`(r, (G', u')) \mapsto (B S_r f)(G', u')`$ on $`[0, t] \times \mathcal X`$. Because $`s \mapsto G_{s-}`$ is left-continuous in $`\mathcal X`$ and $`|H_o(s-)| \le |H_o(t)|`$ with $`\mathbb E|H_o(t)| \le e^t < \infty`$, the sum $`\sum_{u \in H_o(s-)}(B S_{t-s} f)(G_{s-}, u)`$ is predictable and integrable on $`[0, t]`$. Conditioning at each promotion stopping time via the strong Markov property of the Poisson drivers and compensating the rate-$`b_u`$ promotion clocks gives
  $`\displaystyle Q_t f = S_t f + \int_0^t Q_s B S_{t-s} f\, ds`$
  for bounded nonnegative cylinders, and hence for all bounded nonnegative Borel $`f`$ by monotone class.
- **Passage to the exponential-age mixture (2).** Integrating against $`\pi`$ (Tonelli's theorem applies by nonnegativity), using $`\pi Q_t = e^{\mu t}\pi`$ and $`\pi Q_s = e^{\mu s}\pi`$, multiplying by $`e^{-\mu t}`$, and changing variables to $`a = t - s`$ gives
  $`\displaystyle \pi f = e^{-\mu t}\pi(S_t f) + \int_0^t e^{-\mu a}\pi(B S_a f)\, da.`$
  Note that $`\mu > 0`$ holds automatically (both from the repository's lower bounds on $`\mu`$ and directly from $`\pi Q_1 = e^\mu \pi`$, since $`Q_1 1(G, o) > 1`$ everywhere on $`\mathcal X`$ implies $`e^\mu = \pi(Q_1 1) > 1`$). Since $`0 \le \pi(S_t f) \le \|f\|_\infty`$, the boundary term $`e^{-\mu t}\pi(S_t f)`$ vanishes as $`t \to \infty`$, and monotone convergence gives (2):
  $`\displaystyle \pi f = \int_0^\infty e^{-\mu a}\pi B S_a f\, da.`$
- **Newborn law $`\beta`$.** Evaluating (2) at $`f \equiv 1`$ gives $`1 = (\pi B 1)/\mu`$, so $`\pi B 1 = \mu`$ (matching $`\pi(b_o) = \mu`$ in `../operator.md` (`operator.md`) §4). Thus $`\beta := \pi B / \mu`$ is a Borel probability on $`\mathcal X`$. Every tree $`(G^{b_o}, v)`$ in the image of $`B`$ is rooted at a freshly promoted leaf $`v`$, which has $`(L_v, c_v) = (1, 1)`$ and $`D_v = 2`$. Under (2), $`\pi`$ is the law of $`(G_{\mathbf a}, o)`$ started from $`\beta`$ at an independent $`\mathbf a \sim \mathrm{Exp}(\mu)`$.

### 2. Tail constant $`0 < A < \infty`$ and $`\pi(D_o \ge k) \sim A k^{-\mu}`$ (Lines 23–30)

- **Supermartingale convergence and Yule $`L^2`$ domination.** Started from any state in the support of $`\beta`$, $`D(0) = 2`$. Promotions at the tracked root replace a leaf by a core neighbour, leaving $`D = L + c`$ unchanged, while leaf gains increment $`D`$ by $`+1`$ at predictable rate $`L + \sum_{w \sim o} 1/D_w \le L + c/2 \le D`$. Enlarging the probability space with an independent Poisson source to supply the shortfall $`Y(t-) - r_o(t-) \ge 0`$ couples $`D(t) \le Y(t)`$ pathwise for all $`t \ge 0`$, where $`Y(t)`$ is a rate-$`Y`$ pure-birth Yule process begun at $`Y(0) = 2`$, uniformly over the newborn environment. In particular $`\mathbb E[D(t) \mid \mathcal F_s] \le D(s)e^{t-s}`$, so $`e^{-t}D(t)`$ is a right-continuous nonnegative supermartingale and converges almost surely to a finite limit $`W \ge 0`$.
- **Exact moments and Doob maximal bound.** By Lemma 1.4 of `../../proofs/existence.md` (`proofs/existence.md`) with $`k = 2`$, $`\mathbb E[Y(t)^2] = 6e^{2t} - 2e^t`$, so the nonnegative martingale $`M_t := e^{-t}Y(t)`$ satisfies $`\mathbb E[M_t^2] = 6 - 2e^{-t} \le 6`$. Doob's $`L^2`$ maximal inequality gives $`\mathbb E[(\sup_{t \ge 0} M_t)^2] \le 24`$. Since $`\mu \le 1 < 2`$ (from $`e^{\mu t} = \pi(Q_t 1) \le e^t`$), Hölder's inequality gives
  $`\displaystyle \mathbb E_\beta\! \left[\sup_{t \ge 0}\big(e^{-t}D(t)\big)^\mu\right] \le \mathbb E\! \left[\Big(\sup_{t \ge 0} M_t\Big)^\mu\right] \le 24^{\mu/2} < \infty,`$
  establishing $`A := \mathbb E_\beta[W^\mu] < \infty`$ as well as uniform integrability of $`\{e^{-t}D(t)\}_{t \ge 0}`$.
- **Strict positivity $`A > 0`$.** On $`(0, 1]`$, requiring the root's promotion clock to be silent (probability $`e^{-1}`$) and its initial tagged leaf to ring at least twice (probability $`1 - 2e^{-1} > 0`$) has probability $`e^{-1}(1 - 2e^{-1}) > 0`$ independently of the exterior, and yields $`L(1) \ge 3`$. Condition on $`\mathcal F_1`$ on this event; at time $`1`$ the shifted process starts from a state $`(G_1, o) \in \mathcal X`$ with $`\ell := L(1) \ge 3`$ and finite total degree $`d_1 := D(1) < \infty`$. At every later state, the instantaneous drift of $`L`$ is
  $`\displaystyle r_o - b_o = L + \rho_o - \frac{L}{D} \ge L - 1,`$
  which holds for **all** $`L \ge 0`$ (including $`L = 0`$, where $`\rho_o - 0 \ge 0 > -1`$). Applying finite-ball compensation and the spatial cutoff limit (dominated in $`L^2`$ by a Yule process begun at $`d_1`$) to $`e^{-s}(L(1+s) - 1)`$ gives $`\mathbb E_{G_1, o}[L(1+s) - 1] \ge e^s(\ell - 1)`$ for all $`s \ge 0`$, hence $`\mathbb E_{G_1, o}[e^{-s}D(1+s)] \ge \ell - 1 \ge 2`$. Since $`\{e^{-s}D(1+s)\}_{s \ge 0}`$ is $`L^2`$-bounded by $`d_1(d_1 + 1)`$ for each fixed $`(G_1, o)`$, uniform integrability gives $`\mathbb E[e W \mid \mathcal F_1] \ge \ell - 1 \ge 2 > 0`$ on $`\{L(1) \ge 3\}`$. Thus $`\Pr_\beta(W > 0) > 0`$ and $`A = \mathbb E_\beta[W^\mu] > 0`$.
- **Cumulative tail asymptotic.** Since $`t \mapsto D(t)`$ is right-continuous and nondecreasing, $`\{t \ge 0 : D(t) \ge k\} = [T_k, \infty)`$ when $`T_k < \infty`$ and $`\emptyset`$ when $`T_k = \infty`$. Substituting $`f = \mathbf 1_{D_o \ge k}`$ into (2) and applying Tonelli's theorem gives
  $`\displaystyle \pi(D_o \ge k) = \int_0^\infty \mu e^{-\mu a}\Pr_\beta(D(a) \ge k)\, da = \mathbb E_\beta\! \left[e^{-\mu T_k};\, T_k < \infty\right],`$
  so $`k^\mu \pi(D_o \ge k) = \mathbb E_\beta[(k e^{-T_k})^\mu \mathbf 1_{T_k < \infty}]`$. Because independent Poisson drivers have no simultaneous jumps almost surely, $`D(t)`$ increases by unit steps $`+1`$, so $`D(T_k) = k`$ almost surely whenever $`T_k < \infty`$ ($`k \ge 2`$). On $`\{\lim_{t \to \infty} D(t) = \infty\}`$, nonexplosion ($`D(t) \le Y(t) < \infty`$) implies $`T_k \uparrow \infty`$, hence $`k e^{-T_k} = D(T_k)e^{-T_k} \to W`$. On $`\{\sup_{t \ge 0} D(t) < \infty\}`$, $`W = 0`$ and $`\mathbf 1_{T_k < \infty} = 0`$ for all sufficiently large $`k`$. Since $`(k e^{-T_k})^\mu \mathbf 1_{T_k < \infty} \le \sup_{t \ge 0}(e^{-t}D(t))^\mu \in L^1(\beta)`$ for every $`k \ge 2`$, dominated convergence yields $`\lim_{k \to \infty} k^\mu \pi(D_o \ge k) = \mathbb E_\beta[W^\mu] = A`$.

### 3. Exact-degree flux identity and $`\pi(D_o = k) \sim \mu A k^{-1-\mu}`$ (Lines 31–37)

- **Validity of Dynkin compensation for the unbounded-root indicator $`\mathbf 1_{D_o \ge k+1}`$.** Fix $`k \ge 2`$ and put $`f_k(G, o) = \mathbf 1_{D_o \ge k+1}`$. A newborn has degree $`2 < k + 1`$, promotions preserve $`D`$, and transitions at $`x \ne o`$ do not change $`D_o`$; only a leaf gain at $`o`$ when $`D_o = k`$ crosses the threshold, at rate $`r_o = L_o + \rho_o = k - c_o + \rho_o`$. Hence
  $`\displaystyle (A f_k)(G, o) = (k - c_o + \rho_o)\mathbf 1_{D_o = k}.`$
  Although $`f_k`$ is not supported on bounded root degrees (so it lies outside the bounded-ball test class of `../operator.md` (`operator.md`) §4), both $`f_k`$ (radius $`1`$, $`\|f_k\|_\infty \le 1`$) and $`A f_k`$ (radius $`2`$, since $`D_w`$ for the $`c_o \le k`$ neighbours of $`o`$ depends only on $`B_2(G, o)`$, with $`0 \le A f_k \le k`$ because $`\rho_o \le c_o/2 \le c_o`$) are **bounded continuous cylinders** on $`\mathcal X`$. Stopping the $`R`$-truncated system ($`R \ge 2`$) at finite total size, removing the stopping bound via $`\mathbb E|H_o^{[R]}(1)| \le e`$, passing $`R \to \infty`$ via `../operator.md` (`operator.md`) eq. (3) applied to $`f_k`$ and $`A f_k`$, and integrating over $`[0, 1]`$ against $`\pi`$ rigorously yields
  $`\displaystyle \mu \pi(D_o \ge k + 1) = \pi\! \left[(k - c_o + \rho_o)\mathbf 1_{D_o = k}\right].`$
- **Upper and lower squeezing.** Since $`D_w \ge 2`$ for every core neighbour $`w \sim o`$, $`0 \le \rho_o \le c_o/2`$, so $`0 \le k - c_o \le k - c_o + \rho_o \le k`$ on $`\{D_o = k\}`$.
  - From $`\mu \pi(D_o \ge k+1) \le k \pi(D_o = k)`$ and $`(k+1)^\mu \pi(D_o \ge k+1) \to A`$, we obtain $`\liminf_{k \to \infty} k^{1+\mu}\pi(D_o = k) \ge \mu A`$.
  - For any $`0 < \delta < 1`$, nonnegativity of $`k - c_o + \rho_o = L_o + \rho_o`$ on $`\{D_o = k, \, c_o > \delta k\}`$ gives
    $`\displaystyle \mu \pi(D_o \ge k+1) \ge (1 - \delta)k\,\pi(D_o = k, \, c_o \le \delta k) \ge (1 - \delta)k\big[\pi(D_o = k) - \pi(c_o > \delta k)\big].`$
    By Markov's inequality, $`k^{1+\mu}\pi(c_o > \delta k) \le k^{1+\mu} e^{-\theta \delta k}\pi(e^{\theta c_o}) \to 0`$ as $`k \to \infty`$. Rearranging gives $`\limsup_{k \to \infty} k^{1+\mu}\pi(D_o = k) \le \mu A / (1 - \delta)`$; letting $`\delta \downarrow 0`$ completes $`\pi(D_o = k) \sim \mu A k^{-1-\mu}`$.

---

## Independent small-tree checks

In `operator_carbon.py` (`python3 round2/checks/operator_carbon.py`), Checks 3C, 3D, and Planted Error 8 independently verify in exact `fractions.Fraction` arithmetic on small marked core trees:
1. $`A(\mathbf 1_{D_o \ge k+1})(G, o) = (k - c_o + \rho_o)\mathbf 1_{D_o = k}`$ across $`k = 2, \dots, 6`$ and six test trees (including singleton, asymmetric edge, leafless type-B $`(0, 2)`$ cores, and branched cores).
2. Exact Taylor-coefficient equality up to order $`3`$ of the cell birth-time renewal equation $`Q_t f = S_t f + \int_0^t Q_s B S_{t-s} f\, ds`$ for cylinder test functions and $`f \equiv 1`$.
3. Negative controls: replacing $`S_{t-s}`$ by $`Q_{t-s}`$ in the renewal equation, or omitting $`-c_o`$ in $`(k - c_o + \rho_o)\mathbf 1_{D_o = k}`$, produces immediate exact rational discrepancies caught by the checker.
