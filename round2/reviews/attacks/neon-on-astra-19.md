# Neon Rival Attack on Astra Claim 19 (Low-Degree Forest at Fixed Cutoff)

**VERDICT: HOLDS**

**Scope of verdict:** Astra Claim 19 pinned to `ft2-astra@b2f38c5` (`round2/open/low-forest/README.md`), conditional on Astra Claim 10 (`spatial-limit.md`: almost-sure empirical convergence to the one-ended decorated core law $`\pi`$ with $`\pi Q_t = e^{\mu t}\pi`$) and $`\mu \ge 0.565`$ (already three-family verified). This audit does NOT re-prove Claim 10 or the later large-cutoff limit (see separate `neon-on-astra-20.md`).

**Attacker:** Neon lead (Gemini 4 Neon, independent of Astra's GPT family and Fable's Claude family). In Fable's `attacks/fable-on-astra-19.md`, the almost-sure component law and quantitative radius proof were noted as checked only within GPT; this attack independently audits every step and constant of §§1–3 directly from the Gemini family.

---

## 1. Line-by-Line Audit of the Proof

### §1. Sampling laws, ordinary-vertex isolation, and component counting
- Deleting every vertex of original degree $`> N`$ leaves every non-isolated vertex inside a component containing at least one surviving original core ($`2 \le D \le N`$); each such core has at most $`N`$ leaf neighbours, so at most $`(N+1)C_n`$ ordinary vertices lie in non-isolated components. Since $`C_n/n \to 0`$ a.s. ($`\mu < 1`$) and at least $`n - C_n`$ vertices survive, a uniformly sampled ordinary (or surviving ordinary) vertex is isolated with limiting probability $`1`$.
- Every surviving component containing $`k_N \ge 1`$ surviving cores is counted once by summing $`1/k_N`$ over its surviving cores, or equivalently at its unique highest surviving core (whose parent has degree $`> N`$, up to at most one founder component contributing $`\le 1/C_n \to 0`$). Any finite rooted component type of $`r`$ ordinary vertices is measurable in the core-radius-$`(r+1)`$ neighbourhood (sufficient to verify all internal degrees $`\le N`$ and all boundary core degrees $`> N`$). Once $`\pi(R_N < \infty) = 1`$, the countable probabilities on finite types (plus cemetery $`\{D_o > N\}`$) sum to $`1`$, so Claim 10's cylinder convergence upgrades to almost-sure total-variation convergence by Scheffé's lemma; weighting by $`1/k_N \in (0, 1]`$ preserves TV convergence of finite measures. **Verified.**

### §2. Backward-window exponential radius bound $`\pi(R_N \ge k) \le 5 e^{-0.08 k}`$ ($`3 \mid k`$, $`k \ge 64 N^3`$)
Independently re-derived every step and constant:
1. **Pull-back & cell second moment:** Using $`\pi Q_t = e^{\mu t}\pi`$, each core at time $`t`$ lies in the cell $`\mathrm{cell}_t(u)`$ of a unique old core $`u`$ at time $`0`$. Since every core promotes at rate $`L_v/D_v \le 1`$, cell size $`Y_t = |\mathrm{cell}_t(o)|`$ is stochastically dominated uniformly in $`T_0`$ by a rate-$`1`$ Yule process from $`1`$ ($`\mathrm{Geom}(e^{-t})`$), giving $`\mathbb{E}[Y_t^2 \mid T_0] \le 2 e^{2t} - e^t \le 2 e^{2t}`$.
2. **Simple-path decomposition ($`k = 3r`$):** Any simple path of length $`3r`$ from $`v \in \mathrm{cell}_t(o)`$ in the surviving low-degree graph splits into at most one recent ascent (length $`a_1`$), one old core-tree segment starting at $`o`$ (length $`m`$), and one recent descent (length $`a_2`$), so $`\max(a_1, m, a_2) \ge r`$.
3. **Old-edge term ($`m \ge r`$):** There are at most $`N(N-1)^{r-1} \le N^r`$ simple old low paths of length $`r`$ from $`o`$. While all $`r+1`$ vertices on a fixed path stay $`\le N`$, the $`2r`$ disjoint directed sends along the path fire at total rate $`\ge 2r/N`$, incrementing the degree sum (initially $`\ge 2(r+1)`$) by at most $`(r+1)(N-2)`$ before exit. For $`t = k/(16N) \ge 4N^2`$, the threshold $`(r+1)(N-2) \le 2rN < 4rN \le rt/N`$ is below half the Poisson mean $`u = 2rt/N`$. Since Cramér's constant $`\frac{1}{2} - \frac{1}{2}\ln 2 \approx 0.1534 > \frac{1}{8}`$, $`\mathbb{P}(\mathrm{Pois}(u) \le u/2) \le e^{-u/8} = e^{-rt/(4N)}`$. For all $`N \ge 2`$, $`4N^2 \ge 8N \ln N`$, so $`N^r e^{-rt/(4N)} \le e^{-rt/(8N)}`$. Cauchy–Schwarz against $`\sqrt{\mathbb{E}[Y_t^2 \mid T_0]} \le \sqrt{2} e^t`$ gives $`\le 2 e^{(1-\mu)t - rt/(16N)}`$. At $`k \ge 64N^3`$ ($`r = k/3 \ge \frac{64}{3}N^3`$), the exponent is $`\le (\frac{0.435}{16N} - \frac{N}{12})k \le -0.153 k < -0.08 k`$ for all $`N \ge 2`$.
4. **Recent-depth term ($`a_1 \ge r`$ or $`a_2 \ge r`$):** When $`m < r`$, the deep recent segment lies in $`\mathrm{cell}_t(o)`$ itself (counted separately even if $`D_o(0) > N`$) or in an old anchor at initial low-degree distance $`1 \le m \le r-1`$ from $`o`$; total anchor count is $`\le 1 + \sum_{m=1}^{r-1} N(N-1)^{m-1} \le 2 N^r`$ for all $`N \ge 2, r \ge 1`$. Iterated predictable compensation of promotion rates $`\le 1`$ gives $`\mathbb{E}[Z_r(t) \mid T_0] \le t^r/r!`$, so $`\mathbb{P}(\text{depth} \ge r \mid T_0) \le t^r/r!`$. With $`Nt/r = 3/16`$ and $`r! \ge (r/e)^r`$, $`N^r t^r/r! \le (3e/16)^r \approx (0.5097)^r < (9/16)^r`$. Cauchy–Schwarz gives $`\le 2\sqrt{2}\,e^{0.435 t}(3/4)^r \le 3 e^{(0.435/(16N) + \frac{1}{3}\ln(3/4))k}`$. At worst case $`N=2`$, the rate is $`0.013594 - 0.095894 = -0.082300 < -0.08`$.
5. Summing both terms gives $`\pi(R_N \ge k) \le 5 e^{-0.08 k}`$ for every multiple of $`3`$ with $`k \ge 64N^3`$. **Verified.**

### §3. Early-age 4-vertex young-tree test
Independently enumerated all directed ordinary arrows for two births from a newborn with parent degree $`d`$ ($`u = 1 + 1/d`$): star coefficient $`u(u+1) + \frac{3}{4}`$, path coefficient $`\frac{7}{6}u`$, exact star ratio $`\frac{u(u+1)+3/4}{u(u+1)+3/4+7u/6}`$, matching $`1161/1651`$ at $`d=20`$ and $`33/47`$ as $`d \to \infty`$, strictly exceeding the ordinary 4-vertex RFT star probability $`2/3`$ for all $`u \ge 1`$ since $`u^2 - \frac{4}{3}u + \frac{3}{4} = (u - \frac{2}{3})^2 + \frac{11}{36} > 0`$. **Verified.**

---

## 2. Independent Verification Script & Planted Defects

Script: `ft2-neon/round2/outputs/attack_astra_1920/audit_claims_19_20.py` (log `audit_claims_19_20.log`, stdlib only):
- `[Claim 19 §2] PASS`: exact check of both backward-window terms for all $`N \in \{2, \dots,60\}`$ at boundary multiples of $`3`$ ($`k \ge 64N^3`$): `worst_old_rate = -0.20487`, `worst_recent_rate = -0.10060 < -0.08`. Caught **2/2 planted defects** (tightening `-0.08` to `-0.09` fails the relaxed bound at $`N=2`$; doubling $`t = k/(8N)`$ makes $`3e/8 > 1`$ blow up the recent term).
- `[Claim 19 §3] PASS`: exact `Fraction` directed-arrow enumeration matches the formula for $`d \in \{2,3,4,5,10,20,100,1000, \infty\}`$ (`1161/1651` at $`d=20`$, `33/47 > 2/3` at $`d=\infty`$); caught planted omission of parent arrow $`1/d`$.
- `[Claim 19 §1] PASS`: checked `2,000` random friend-tree/cutoff pairs ($`n \in [5, 120], N \in \{2,3,5,10\}`$), confirming exact identity between $`1/k_N`$-weighted surviving-core sum and top-core count (+ founder exception) and the $`(N+1)C_n`$ bound.

**What is NOT proved by Claim 19 (correctly stated in the claim):** Uniform-in-$`N`$ bounds, simultaneous $`n, N \to \infty`$ limits, and finiteness of mean component volume.
