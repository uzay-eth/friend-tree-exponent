# Neon Rival Attack on Astra Claim 17 (Lebesgue Equivalence & Exact Right-Tail Rate $`I`$)

**VERDICT: HOLDS**

**Scope of verdict:** Astra Claim 17's two main strengthenings that had not previously received a non-GPT audit:
- Commit `044f411` (`round2/open/limit-law/density/README.md`): equivalence of $`\mathrm{Law}(X)`$ to Lebesgue measure on $`(0, \infty)`$ (density positive a.e.), conditional on Astra Claim 6 and the inert-seed limit (`51a2a3e`, `open/forest-age/inert/README.md`).
- Commit `d65ead7` (`round2/open/limit-law/right/README.md`): existence of the exact right-tail rate $`I = \lim_{x\to\infty} (-\log\mathbb{P}(X > x))/x^{1/(1-\mu)} \in (0, \infty)`$ and its finite-block variational formula, conditional on Astra Claim 6 and the uniform finite-time free one-leaf gradient (`7ce6253`).
(The original support and matching left/right tail orders at `7ce6253` and `f7e025c` were already independently confirmed in Neon Fact F7 / `round2/open_questions/REPORT.md` §2 and in Fable's `fable-on-astra-17.md`; the left-tail prefactor commit `02a5b00` is outside this report, while its equivalent statement in Fable Claim 2(d) passed Gemini audit in `neon-on-fable-2.md`.)

**Attacker:** Neon lead (Gemini 4 Neon, independent of Astra's GPT family and Fable's Claude family).

---

## 1. Line-by-Line Audit of `density/README.md` (`044f411`)

1. **Finite-crossing adjacent hub edge:** From any core tip $`x`$ ($`c_x \le 1, L_x \ge 1`$), silence of $`x`$ for unit time ($`e^{-1}`$) and $`\ge 2`$ rings of one tagged leaf ($`1 - 2e^{-1}`$) followed by promotion ($`\ge 3/4`$) has probability $`\ge q = \frac{3}{4}e^{-1}(1 - 2e^{-1}) > 0.072`$. Coupling leaf families at both core endpoints $`x, y`$ ($`D_x, D_y \ge 2 \implies 2/D \le 1`$) below by independent linear birth–death chains (rates $`1`$ and $`1/2`$, normalized second moment $`4 - 3e^{-t/2} \le 4`$, survival $`\ge 1/4`$ each) gives conditional probability $`\ge q/16`$ that both reciprocal degree integrals are finite. Lévy's 0–1 law gives probability $`1`$; finite crossing compensator $`\int (D_x^{-1} + D_y^{-1})\,dt < \infty`$ forces finitely many crossings a.s., hence a deterministic integer cut time $`j`$ and core edge $`e`$ with zero subsequent crossings. **Verified.**
2. **Absolute continuity via one initial holding time $`T`$:** The no-crossing subprobability density $`w = \exp(-\int_0^\infty (D_1^{-1} + D_2^{-1})\,dt) \in (0, 1]`$ relative to the independent inert product law $`Q`$ is strictly positive on good paths. Under $`Q`$, varying the first half's independent initial holding time $`T \sim \mathrm{Exp}(m - 1/D_o)`$ gives $`r = e^{-T} \in (0, 1)`$ and combined ratio $`f(r) = (B r^\mu + D)/(A r + C)^\mu`$ with exact derivative $`f'(r) = \mu(BC r^{\mu-1} - AD)/(Ar+C)^{\mu+1}`$, which has at most one zero on $`(0, \infty)`$. Smoothness and local invertibility on each regular branch imply zero mass on every Lebesgue-null set; countable union over $`(T_{\mathrm{seed}}, e, j)`$ gives $`\mathrm{Law}(X) \ll \mathrm{Leb}`$. **Verified.**
3. **Every positive-Lebesgue-measure set has positive probability ($`\mathrm{Leb} \ll \mathrm{Law}(X)`$):** Since $`1 \le p^\mu + (1-p)^\mu \le 2^{1-\mu}`$ for $`p \in (0,1)`$, the combined ratio lies in $`[\min(x_1, x_2), 2^{1-\mu}\max(x_1, x_2)]`$, so full support $`(0, \infty)`$ guarantees reachable planted seeds with positive-probability good paths having intrinsic ratios $`< l`$ and $`> h`$. Joining them by a single edge and varying both independent initial waits $`(T_1, T_2)`$ gives $`r = e^{T_2 - T_1}`$ strictly positive density on $`(0, \infty)`$, sweeping $`f(r)`$ continuously between $`x_2`$ and $`x_1`$ across all of $`(l, h)`$ with $`w > 0`$ everywhere on the good continuations. Exhausting $`(0, \infty)`$ proves equivalence to Lebesgue measure. **Verified.**

---

## 2. Line-by-Line Audit of `right/README.md` (`d65ead7`)

1. **Deterministic join bound $`|\Phi(\text{joined}) - \sum \Phi(\text{free})| \le Gk`$ (§1):** At each vertex with internal degree $`d \ge 1`$ and $`c`$ cut edges, both internal-arrow deficit and extra cross-arrow rate equal $`c/(d+c)`$, summing to $`\le \sum 2c_v = 4k`$. Combining with the finite-time free one-leaf gradient $`L_f e^{\mu u}`$ and integrating $`4k L_f \int_0^t e^{\mu(t-s)}\,ds`$, dividing by $`e^{\mu t}`$ ($`2k`$ endpoint status difference vanishes) gives $`G = 4L_f/\mu`$. **Verified.**
2. **Scaffold setup, balanced Pólya urn, and unconditional root-source trimming (§2):** Degree-$`\le 3`$ scaffold of size $`3k`$ costs $`\le 12^{3k-3}`$; Stage 1 (add 1 leaf per anchor via scaffold neighbor of deg $`\le 4`$, size $`\le 4k`$) gives $`\ge j/(16k)`$; Stage 2 (promote via anchor of deg $`\le 4`$, size $`\le 5k`$) gives $`\ge j/(20k)`$. No-scaffold draw over $`L = k(m-2)`$ steps has probability $`\prod_{i=0}^{L-1} \frac{2k+i}{5k+i} \ge e^{-c_1 k(1+\log m)}`$. Conditional on no scaffold, the $`k`$-color Pólya($`2, \dots,2`$) composition weight is proportional to $`\prod_{r=1}^k (\ell_r + 1)`$, maximized by AM–GM at the balanced composition $`\ell_r = m-2`$ (probability $`\ge \binom{L+k-1}{k-1}^{-1} \ge (em)^{-k}`$). Unconditional root-source count $`R_m = \sum_{j=2}^{m-1} \mathrm{Bern}(1/j)`$ satisfies $`\mathbb{P}(R_m > s_m) \le (e \log m / s_m)^{s_m} \ll e^{-(D+1)m}`$ at $`s_m = \lceil m/\sqrt{\log m}\rceil = o(m)`$, paying at most $`2^{-s_m}`$ per block for the $`d/(d+1) \ge 1/2`$ root-source likelihood ratio. **Verified.**
3. **Matching $`I_+ \le \lim_{\delta\downarrow 0} V(\delta) \le I_-`$ & Hoeffding lower bound (§3, §5):** Fixed-block repetition gives $`I_+ \le V(\delta)`$. Along a sequence attaining $`I_-`$, choosing $`m = \lceil M x^p \rceil`$ uses the exact identity $`2 + p(2\mu - 1) = p`$ ($`p = 1/(1-\mu)`$) so conditional concentration eliminates $`\{\Phi(T_m) \le (1-2\varepsilon) x m^\mu\}`$ at exponent $`c\varepsilon^2 M^{2\mu-1} x^p > I_- x^p`$, yielding $`V(\delta) \le I_-/(1-2\varepsilon)^p`$. In §5, Hoeffding's lemma with conditional range $`2L\gamma_{j+1}`$ and $`\sum_{j \ge m}\gamma_{j+1}^2 \sim m^{1-2\mu}/(2\mu-1)`$ ($`2\mu > 1`$) yields $`I \ge \frac{2\mu-1}{2L^2}\max_t t^{2\mu-1}(1-B t^a)^2`$ with exact unique maximizer $`B t^a = 2\mu - 1`$. **Verified.**

---

## 3. Independent Verification Script

Script: `ft2-neon/round2/outputs/attack_astra_17/audit_claim17.py` (log `audit_claim17.log`, stdlib only):
- `[Claim 17 Density] PASS`: verified $`f'(r)`$, unique turning point, limiting endpoints $`x_1, x_2`$, sandwich $`[1, 2^{1-\mu}]`$, preparation bound $`q > 0.07`$, and caught planted defect $`2^{1-\mu} \to 1`$.
- `[Claim 17 Right Rate] PASS`: verified exact `Fraction` identity $`2 + p(2\mu - 1) = p`$, exact Pólya($`2,2,2`$) balanced-composition maximality, exact Hoeffding optimizer $`B t^a = 2\mu - 1`$, root-source Chernoff decay, and caught planted defect at $`\mu = 1/2`$.
- Also independently replayed Astra's own `round2/open/limit-law/right/check.py` (all 3 assertions passed).

**What is NOT proved by Claim 17 (correctly stated in the claim):** Smoothness/formula for the density of $`X`$, numerical value of $`I`$, or existence of an optimizer at fixed density $`b`$.
