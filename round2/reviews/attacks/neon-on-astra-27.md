VERDICT: HOLDS

# Independent Gemini-Family Audit of Astra Claim 27 (All-Label Strict Order for Every Positive Moment $`p > 0`$)

• **Target:** Astra Claim 27 (*Every positive moment of limiting degree share decreases with birth label*): for every real $`p > 0`$,

```math
\mathbb{E}[Z_1^p] = \mathbb{E}[Z_2^p] > \mathbb{E}[Z_3^p] > \mathbb{E}[Z_4^p] > \cdots,
```

with uniform pre-placement conditional advantage $`\mathbb{E}[Z_{n+1}^p - Z_{n+2}^p \mid \mathcal{F}_n] \ge \frac{9}{5000}\frac{\Gamma(1+p)\Gamma(n)}{\Gamma(n+p+2)}`$ for every $`n \ge 2`$, founder gap $`\mathbb{E}[Z_2^p - Z_3^p] \ge \frac{9}{5000(p+1)(p+2)}`$, conditioned Beta-function telescoping $`\mathbb{E}[Z_i^p - Z_j^p \mid \mathcal{F}_m] \ge \frac{9}{5000(p+1)}\{\mathrm{B}(i-1, p+1) - \mathrm{B}(j-1, p+1)\}`$, and extension to nonzero finite positive mixtures of powers.
- **Pinned commit audited:** `ft2-astra@37e7bf685a16eb0c789dba3d7097ebcf5aaf44c4` (building on `239d9db` for $`p \ge 1`$ and `aca9937`).
- **Auditor:** Neon lead (`gemini-4-neon`, pure Gemini-family audit; all 9 proof documents read directly and core identities/barriers/trees independently reconstructed from scratch with zero shared code from `ft2-astra`).
- **Scope:** Full Claim 27 across the entire real half-line $`0 < p < \infty`$, including all four ranges $`(0, 473/512]`$, $`[473/512, 1]`$, $`[1, 3]`$, $`[3, 7]`$, and $`[7, \infty)`$, seed/star boundary cases, the founder seam for both $`0 < p \le 1`$ and $`p \ge 1`$, and Beta-function telescoping.

---

## 1. Independent Mathematical Audit of Every Proof Link

I read all nine proof notes directly (`all-powers/positive/README.md`, `all-powers/README.md`, `mellin/README.md`, `mean-ledger/real/README.md`, `mean-ledger/fractional/README.md`, `fractional-middle/README.md`, `tip-moments/README.md`, `small-powers/README.md`, `mean-ledger/README.md`) and rederived every key step:

1. **Common-clock Gamma factorization:** Starting from any fixed completed tree $`T_m`$ of size $`m = n + 2`$, the total-size Yule clock $`N(t)`$ has holding rates $`j`$ depending only on size and is independent of the embedded jump chain, so $`S = \lim e^{-t} N(t) \sim \mathrm{Gamma}(m, 1)`$ is independent of $`Z_i`$. Hence $`W_i = S Z_i`$ satisfies $`\mathbb{E}[W_i^p \mid T_m]/\Gamma(1+p) = \frac{(m)_p}{\Gamma(1+p)}\mathbb{E}[Z_i^p \mid T_m]`$, and multiplying by the exact two-birth identity gives $`\frac{\Gamma(1+p)}{n(n+1)(n+2)_p} = \frac{\Gamma(1+p)\Gamma(n)}{\Gamma(n+p+2)}`$. At later random stopping times, the proof relies strictly on exponential clock discounting $`e^{-p\tau}`$, never on post-stopping Gamma independence.
2. **Exact discounted harmonic stock boundary ($`0 < p \le 1`$) and shifted submartingale ($`p \ge 1`$):**
   - For $`0 < p \le 1`$, $`\varphi_p(j) = e^{-1}\sum_{k=0}^{j-1} P_p(j-k)/k!`$ has generating function $`H(z) = e^{z-1} z (1 - z)^{-p-1}`$ solving $`z(1-z)H' + (z^2 - (1+p)z - 1)H = 0`$, which yields the exact harmonic equation $`j\varphi_p(j+1) + \varphi_p(j-1) - (j+1+p)\varphi_p(j) = 0`$ for all $`j \ge 1`$. Stopping the alternating $`e^{-1}`$ series at odd index $`17`$ gives a strict lower bound on $`e^{-1}`$.
   - For $`p \ge 1`$, $`\phi_p(j) = (j-1)_p / \Gamma(1+p)`$ has exact lower-generator drift $`p\phi_p(j)(\frac{1}{j-1} - \frac{1}{j+p-2}) \ge 0`$ for all $`j \ge 2`$.
   - Uniform integrability for $`t \to \infty`$ is provided by the pure Yule total-size process dominating $`D(t) \ge j(t)`$.
3. **Grid-free alternating-sum bound for $`p \ge 7`$:**
   - Rederived $`(\mathcal{L} - q) P_q(j) = -P_{q-1}(j)`$ identically for all $`j \ge 1`$ and $`q \ge 1`$ (birth term cancels $`-q P_q(j)`$ identically; death term shifts $`q \to q-1`$).
   - Under $`u = 1 - e^{-t}`$, $`r_q(u) = e^{-qt}\mathbb{E}[P_q(J_t)]`$ satisfies $`r_q'(u) = -r_{q-1}(u)`$, $`r_q(0) = 1`$, $`0 \le r_q \le 1`$. Iterating 7 times yields $`M_p \ge \sum_{j=0}^7 (-1)^j/j! = 103/280`$.
   - Confirmed the old-leaf edge-weight bound $`3 X \le A = \sum_u w_u \ge 1`$ (checked analytically and on all 986 free trees $`n = 2\dots 12`$), yielding exact margin $`103/280 - (4 A_7(1/2) + 2 A_7(1/3))/3 = 22890923/1058158080 > 21/1000 > 9/5000`$.
4. **Euler product upper bound, stock-sensitive tips, low-exponent host boost, and high-degree tails:**
   - Verified Weierstrass product identity $`A_p(a) = \prod_{k \ge 0}\bigl(1 + \frac{p(a-1)}{(k+1)(k+a+p)}\bigr)`$, monotonicity in $`p`$ for both $`a \le 1`$ and $`a > 1`$, and tail bound $`\exp(p(a-1)/T) \le 1/(1 - p(a-1)/T)`$ for $`a > 1`$.
   - Verified `tail(k, j) = (1 + 1/k)\min(1/(j+1), 1/(1 + \mathrm{lo}\,k))` via first-send barrier + Laplace hazard domination, optional stopping at $`\sigma = \tau \wedge T_{\mathrm{stop}}`$, Möbius-function monotonicity of the host-boost race $`\min(x/y, (x+uH)/(y+u))`$, and both $`x \le 1/2`$ and $`x > 1/2`$ regimes of the $`d \ge D`$ analytic tail (including $`x(\frac12 - \frac{d}{d+1}x)_+ \le \frac{1}{4(d+1)}`$ on $`x \ge 1/2`$).
5. **Founder seam for all $`p > 0`$ and Beta telescoping:**
   - For $`0 < p \le 1`$: $`\varphi_p(2) = (2+p)e^{-1}`$ with Jensen $`\Gamma(1+p) = \mathbb{E}[e^{p\log E}] \ge e^{-p}`$ gives $`C_p = (2+p)e^{-1-p}`$; Hölder gives $`L_p = a^{1-p} b^p`$ ($`a = 197/305`$, $`b = 387/1000`$); strict concavity of $`\log(C_p/L_p)`$ with $`e < 11/4`$ gives $`C_p/L_p > 51/50`$, hence $`C_p - L_p \ge b/50 = 387/50000`$; dividing by $`2(3)_p`$ and using $`\Gamma(1+p) \le 1`$ yields $`\mathbb{E}[Z_2^p - Z_3^p] \ge \frac{387/50000}{(p+1)(p+2)} > \frac{9}{5000(p+1)(p+2)}`$.
   - For $`p \ge 1`$: centre $`\ge 1`$, leaf $`\le A_p(1/2) \le 1/2`$ gives $`\frac{1}{2(p+1)(p+2)} > \frac{9}{5000(p+1)(p+2)}`$.

---

## 2. Independent Zero-Shared-Code Audit Evidence

Standalone script: `round2/outputs/attack_astra_27/audit_claim27.py` in `exact/neon`.

Command run:
```sh
python3 -IB round2/outputs/attack_astra_27/audit_claim27.py
```

Observed output (exit 0, wall `8.23s`):
- Verified exact harmonic equation for $`\varphi_p(j)`$, submartingale drift for $`\phi_p(j)`$, and $`(\mathcal{L}-q)P_q = -P_{q-1}`$; recomputed exact $`p \ge 7`$ margin `22890923/1058158080` (`0.021632801`), exact tip-moments ceilings `U_1(2,3) = (0.386883135517, 0.264826886411)` and `U_2(2,3) = (0.301580912392, 0.184391218060)`, founder gap $`\ge 9/[5000(p+1)(p+2)]`$ across all $`p > 0`$, Euler product upper bound against `math.lgamma`, and `900` exact Beta telescopes across 15 powers from $`10^{-18}`$ to $`1000`$.
- **Exact free-tree enumeration (`n = 2..12`, all `986` non-isomorphic free trees x `11` powers = `10,846` exact rational checks, `2.57s`):** Verified $`A \ge 1`$, $`3X \le A`$, and $`N_p(T) > 9/5000`$ on every tree across $`p \in \{1/512, 1/4, 1/2, 473/512, 1, 3/2, 2, 3, 5, 7, 10\}`$ (minimum ratios `12.0x` to `131.2x` above `9/5000`).
- **Planted-error controls (6/6 caught):** even truncation of $`e^{-1}`$ at $`j=16`$ (`+2.66e-15`), dropping return-race $`B_d(2)=0`$ at $`p=1`$ (`-3.38e-2`), wrong Euler endpoint `lo` at $`a=3/2`$ (`-8.08e-4`), dropping the $`(4/3, 2/3)`$ leaf-edge split at $`p=3`$ (`-7.58e-2`), even alternating sum at $`j=6`$ (`+1.76e-4`), and dropping host `boost()` at cell `[0, 1/512]` (`-2.13e-2`).

---

## 3. Complete Clean-Archive Replay (`ft2-astra@37e7bf6` + `aca9937`)

Unpacked clean `git archive` snapshots and executed **all 10 official Claim 27 checkers** (all exited 0):
1. `all-powers/positive/check.py` (at `37e7bf685a16eb0c789dba3d7097ebcf5aaf44c4`): `PASS 9 source pins, 512 fractional boxes joined, all seed/star and founder bounds, 880 exact beta telescopes`.
2. `mean-ledger/fractional/star_bound.py` (at `37e7bf6`): `PASS all 512 boxes: normalized seed and star numerator >= 9/5000` (exact global minima matched: `406542520567/2500000000000`, `2085176012479/12000000000000`, `265041348327659/2012097656250000`).
3. `mean-ledger/fractional/check.py` (at `37e7bf6`): `PASS 473 cells [0,473/512]; normalized numerator >= .033 with at least two cores`.
4. `fractional-middle/check.py`: `PASS 256 boxes` ($`>0.55`$), `224 boxes` (`[1/8,1]`), `39 boxes` (`[473/512,1]`).
5. `all-powers/check.py`: `PASS five writer pins, three uniform margins, nine founder checks and 711 exact beta telescopes`.
6. `mellin/check.py`: `PASS 256 lower-power boxes and 2048 signed boxes on [3,7], minimum 1892360505093/1000000000000000`; `PASS p>=7 normalized serial lower 103/280 signed margin at least 22890923/1058158080`.
7. `mean-ledger/real/check.py`: `PASS counts {'sitewise': 87, 'suppressed': 169}; PASS entire closed exponent interval [1,3], uniform normalized numerator 1/40`.
8. `mean-ledger/check.py`, `tip-moments/check.py`, `small-powers/check.py`: all passed cleanly.

No gaps were found in any part of Claim 27. **VERDICT: HOLDS.**
