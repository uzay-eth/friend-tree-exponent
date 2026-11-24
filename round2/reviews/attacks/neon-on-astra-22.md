VERDICT: HOLDS

Target: Astra Claim 22 (`ft2-astra@4233b1a`, `round2/open/age-order/spectrum/quenched/README.md` & `spectrum/README.md`, plus moments addendum at `ft2-astra@ac7d4e3`, `round2/open/age-order/moments/README.md`).
Reviewer: Neon lead (Gemini family — independent of both Astra's GPT family and Fable's GPT audit child; fills the non-GPT proof audit gap noted in `attacks/fable-on-astra-22.md`).
Date: 2026-11-24 13:13 UTC.

## Scope & Summary

Conditional on the cited inputs (Claims 6, 10, 11, and the late-birth / scaled-tail estimates added to Claim 18), this independent Gemini-family audit read the three proof texts line by line and independently tested the finite algebraic, martingale, and envelope identities from scratch (`round2/outputs/attack_astra_2223/audit_claims_22_23.py`, zero shared code), plus replayed Astra's own checks:

1. **Quenched birth-label and mass empirical process (`quenched/README.md`) — HOLDS:**
   - **Continuous-clock cohort intensity (§1):** Checked the exact threshold alignment $`e^{-(T-s)} D_v(T) > r e^{s-t} \iff D_v(T) > r e^{T-t}`$, the uniform $`0 \le F_q^u \le 7`$ bound (direct core promotion $`\le 6`$ + old-leaf births $`\le 1`$), sequential continuity of $`F_{q_i}^{u_i} \to F_q`$ via the uniform terminal $`L^1`$ bound + positive atomlessness, and eq. (3)–(4) via Claim 10 local LLN + dominated convergence.
   - **Backward-martingale fluctuation without hub-independence assumption (§2 — the crux):** Verified that $`M_s = \mathbb{E}[A_t \mid \mathcal{F}_s]`$ and $`P_s = \mathbb{E}[I_t \mid \mathcal{F}_s]`$ satisfy $`M_{s_0} = P_{s_0}`$, that at any jump $`s \in (s_0, t+d]`$ the past integral has zero jump while the remaining forecast equals the conditional mean of births *strictly after* $`s`$, and that the dirty-core coupling bounds the combined dirty-core count by $`3 e^{2(T-s)}`$ (at most 3 initial dirty cores on a leaf attachment; each dirty core rings at rate $`1`$ and adds at most $`2`$ combined dirty cores; dirty leaves only add private leaves at dirty cores; all clean cores cancel identically). Subtracting the two martingales gives $`|\Delta(M_s - P_s)| \le 6 e^{2(T-s)}`$, and stopped jump isometry with $`\mathbb{E}[|T_s|] = 3 e^s`$ yields $`\mathbb{E}[(e^{-\mu t}(A_t - I_t))^2] \le 36 e^{-3l} \exp\{-(2\mu - 1)t + 4(T - t)\}`$, decaying exponentially for $`T = t + 2\log(t+2)`$ because $`2\mu - 1 > 0.131 > 0`$.
   - **Terminal degree vs limiting mass & mesh interpolation (§3):** Checked the history gate $`\max_v c_v(T) = O(t) \ll k = (t+2)^2`$, super-exponentially summable Gamma sandwich error $`K\mathbb{E}[|T_T|]\exp(-c\varepsilon^2 r k)`$, and the exact lower/upper translated-window inclusions on $`[s, s+1/m]`$ (verified in 64 exact checks).
   - **Clock transfer, unbounded rectangles, and sampling law (§4):** Verified pathwise substitution $`B_v - t = \log(xv/Y) + o(1)`$ and exact cancellation $`H Y^{-\mu} = X`$ under $`y = Y e^z`$; unbounded endpoints $`a = 0, b = \infty`$ closed simultaneously on one probability-one event by matching the total Claim 11 count $`K X r^{-\mu}`$ via the Mellin identity $`K = \mu^{-1}\int u^\mu \nu(du)`$; $`L^1`$ convergence for each fixed rectangle by domination by the $`L^1`$-convergent total count; and product factorization of $`(v Z_v, Z_v/x)`$ under realized-tree uniform sampling.

2. **All positive moments & adjacent differences (`moments/README.md`) — HOLDS:**
   - Verified the stopping-size $`\mathrm{Beta}(2, m-2)`$ survival function $`Q_m(t) = (1-t)^{m-2}[1+(m-2)t]`$, its exact monotonicity $`Q_m(t) - Q_{m+1}(t) = (m-1)t^2(1-t)^{m-2} \ge 0`$, and exact moments (308 exact `Fraction` checks).
   - Verified the signed two-birth identity $`n(n+1)\mathbb{E}[\mathrm{test}(Z_{n+1}) - \mathrm{test}(Z_{n+2}) \mid T_n] = \sum_u \bigl(\frac{r_u}{D_u+1} I(u) - r_u \delta_u \sum_{w\sim u} H(u, w)\bigr)`$ and the per-core charging bound $`\frac{2 r_o}{D_o+1} + \frac{L_o}{D_o} < 3`$ on all 152 distinct rooted/labeled trees for $`n=3..6`$ with nonlinear tests, confirming that flipping the edge sign is cleanly detected on an asymmetric 5-vertex tree.
   - Verified the integration by parts from $`B(y) = (1-\mu)A(y) - y A'(y)`$ giving $`a_\alpha = \alpha \int_0^\infty y^{\alpha-1} B(y)\,dy = (1+\alpha-\mu) b_\alpha`$.

3. **Independent & Replay Checks:**
   - Independent pure-Python audit (`audit_claims_22_23.py`): 308 exact Beta checks, 152 exact signed two-birth checks (+ flipped-sign planted defect), exact discrete backward-martingale identity $`\mathbb{E}[A]=\mathbb{E}[I]=239/1080`$ and post-jump exclusion (+ planted inclusion defect), 64 translated-window inclusion checks — ALL PASS.
   - Clean replay of Astra's `quenched/check.py` (108 forecast rows, 102 controls caught) and `moments/check.py` (96 signed moment rows, 156 Beta rows) — PASS.

## Reproduction

```sh
python3 -B round2/outputs/attack_astra_2223/audit_claims_22_23.py
(cd [local path] && python3 -B round2/open/age-order/spectrum/quenched/check.py && python3 -B round2/open/age-order/moments/check.py)
```
