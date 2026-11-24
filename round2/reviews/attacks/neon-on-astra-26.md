VERDICT: HOLDS

# Independent Gemini-Family Audit of Astra Claim 26 (All-Label Hub Ordering)

- **Target:** Astra Claim 26 (*Birth order of the hubs, for every label*: $`p_1 = p_2 > p_3 > p_4 > \cdots`$, conditional two-birth lower bound $`371/[12500\,n(n+1)]`$, founder gap $`p_2 - p_3 \ge 127/2440`$, and conditioned telescoping intervals $`c(1/(i-1) - 1/(j-1))`$).
- **Pinned commits audited:**
  - `ft2-astra@b450a310548c62adcfb37d66c0415d9dc3f7b3fc` (`round2/open/age-order/global-transport/hostile/RELEASE.md` and all 18 pinned source files).
  - `ft2-astra@212dc9f59b9f2c660ec71734826a5938daa6cab9` (`round2/open/age-order/hub-ranking/README.md` and `check.py`, verified as a clean descendant with zero diff on the 18 pinned files).
- **Auditor:** Neon lead (`gemini-4-neon`, pure Gemini-family audit; all proof notes read directly and all computations independently reconstructed from scratch with zero shared code from `ft2-astra`).
- **Scope:** Covers the full zero-threshold hub-probability theorem ($`p_1 = p_2 > p_3 > p_4 > \cdots`$), the all-tree conditional two-birth inequality $`\Pr(H_{n+1}\mid\mathcal F_n)-\Pr(H_{n+2}\mid\mathcal F_n)\ge \frac{371}{12500\,n(n+1)}`$ for every $`n\ge 2`$ (including seed, stars, and all multi-core trees), the founder bound $`p_2 - p_3 \ge 127/2440`$, and the conditioned multi-step telescoping consequences in `hub-ranking/README.md`. The positive-threshold addenda (`[0, 2^{-131400}/(N+1)]` and `[0.49, 1)` at `18b7abc`) are **out of scope** of this audit.

---

## 1. Independent Mathematical Audit of Every Proof Link

I read `hostile/RELEASE.md`, `hostile/README.md`, `early/README.md`, `hub-ranking/README.md`, `adjacent/README.md`, `preplacement-transfer/README.md`, `hostile/discharge/README.md`, `hostile/mixed-stock/README.md`, `hostile/edge-dual/README.md`, `hostile/interfaces/README.md`, and `hostile/interfaces/low-end.md` directly and independently rederived each step:

1. **Exact two-birth identity:** Starting from any finite tree $`T`$ of size $`n\ge 2`$ with $`r_u = \sum_{z\sim u} 1/D_z`$, $`w_u = r_u/(D_u+1)`$, and $`s_u = w_u/D_u`$, pairing the two orders of distinct targets cancels nonadjacent pairs and exchangeable repeated-target leaves, leaving
   $`\displaystyle n(n+1)\bigl(\Pr(H_{n+1}\mid\mathcal F_n)-\Pr(H_{n+2}\mid\mathcal F_n)\bigr) = \sum_{u\in V(T)} w_u I(u) + \sum_{\{u, v\}\in E(T)} (s_v - s_u) H(u, v).`$
   Confirmed algebraically and by exact brute-force enumeration on all 94 non-isomorphic free trees $`n=2\dots 9`$ (6,902 exact trajectories).
2. **Continuous-time couplings and scalar barriers:**
   - **Return-race lower bound $`R_d(k)`$ and tracked-stock rescue $`B(d, l)`$:** Stopping at the first time *either* the newly promoted degree-$`k`$ neighbour returns to $`z`$ *or* any other neighbour sends to $`z`$ holds $`D_z \equiv d`$ up to $`\tau`$ and lands in state $`(d+1, 1)`$ on win; dropping exterior sends to $`z`$ is conservative. At $`d \ge 20`$, tracked leaves dominate a birth–death process with birth rate $`1`$ and death rate $`1/d`$, giving survival $`\ge 1 - d^{-l}`$ (and $`\{D_z\to\infty\}=\{Z_z>0\}`$ a.s. by the hub dichotomy).
   - **Serial middle-minus-tip gap $`G(p)`$:** Watching the joint clocks out of $`(D_M=2, L_M=1, D_u=p)`$ (parent send to $`M`$ at rate $`a=1/p`$, tip send to $`M`$ at rate $`1`$, $`M`$ send to tip at rate $`1/2`$, exterior degree gain of $`u`$ at rate $`\gamma \le p - 1/2 < p`$) and substituting $`B(3,2)+\frac12 B(2,0)=\frac32 B(2,1)`$ yields
     $`\displaystyle h(M) \ge B(2,1) + \frac{a(B(3,2)-B(2,1))}{a + \frac32 + \gamma} \ge B(2,1) + \frac{a(B(3,2)-B(2,1))}{p + \frac32 + a}.`$
   - **First-send upper bound $`P(d, j)`$:** Whenever the parent retains at least one core neighbour ($`d \ge j+2`$), which holds in every single serial and parallel call, $`f(j) = 1/(j+1)`$ has exact generator drift $`\le \frac{j+2-d}{d(j+2)} \le 0`$ for all $`j\ge 1`$ and $`-\frac{1+r_{\rm ext}}{2} < 0`$ at $`j=0`$. Together with $`P(d+1, j+1) \le P(d, j)`$ on $`2\le d < 112`$, $`M_t = \mathbf 1_{\{\tau\le t\}} + \mathbf 1_{\{\tau>t\}} P(D_p(t), J_p(t)) \in [0,1]`$ is a bounded supermartingale, giving $`\Pr(\tau<\infty) \le P(d, j)`$.
   - **89ae's $`+1`$ tick rounding clarification:** Verified directly. Rounding $`U(d, j)=\lceil 10^{12}P(d, j)\rceil/10^{12}`$ violates literal $`U\le 1/(j+1)`$ by $`<10^{-12}`$ when $`(j+1)\nmid 10^{12}`$ (e.g. $`U(113,2) - 1/3 = 2/(3\cdot 10^{12})`$). All finite checks use the rounded grid $`U`$ directly, and every analytic tail inequality ($`1/228 < 0.017`$, $`1/86 < 0.012`$, $`1/4 < L(2)`$) holds for the unrounded physical bound $`P(d, j)\le 1/(j+1)`$ and still retains wide slack even with $`+10^{-12}`$ added.
3. **All-tree discharging & convexity reductions:**
   - **Seed and stars:** Seed $`n=2`$ gives $`G(2) > 0.18940948 > 371/12500`$; stars $`d=2,3`$ give $`>0.20181908`$ and $`>0.22577249`$; stars $`d\ge 4`$ have zero edge debit because $`P(d+1, d-1)\le 1/d \le 1/4 < L(2)`$, leaving leaf credit $`\ge G(2)/2 > 0.0947 > 371/12500`$.
   - **High internal donors ($`d\ge 10, c\ge 2`$):** For fixed $`(d, c)`$, the margin $`C(u) + A(c-2) - B - c(s_u U_u - B)_+ - c_F F - c_D D_2`$ is **concave** in $`(R, c_0, c_F, c_D)`$ over the reciprocal polytope $`\frac14 c_F + \frac12 c_D \le R \le \frac12 c`$, so its minimum lies on the 5 tested vertices $`((0,0),(0,2),(F,4),(F,2),(D_2,2))`$. For $`d\ge 113`$, verified both $`s_\ell \ge s_u`$ (where $`l = (d+1)/2`$ gives zero leaf toll and $`l \le \lfloor d/2\rfloor`$ gives $`c \ge d/2 \ge 113/2`$) and $`s_u > s_\ell`$ (where $`w_u > 1/2 > 113/228`$, leaf toll is $`0`$ for $`l\ge 4`$, and covered by $`\frac{l}{2d}G(2)`$ for $`l\in\{2,3\}`$ since $`G(2) > 2/114`$).
   - **Bare degree-2 chain suppression and simultaneous low endpoints:** Because $`T_{\rm core}`$ is a finite tree with $`|V_{\rm core}|\ge 2`$, it has at least $`K\ge 2`$ retained cores. Chain debit is convex in endpoint parameters $`(s_A, s_Z)`$ and in $`1/d_{\rm high}\in[0,1/31]`$ for $`d_{\rm high}\ge 31`$ (where $`L(d+1)=0`$), and every interior vertex beyond $`k=3`$ has $`s=1/6`$, zero extra edge debit, and positive credit $`\frac13 G(3) > 0`$. Every chain port has actual graph neighbour degree $`e=2`$, so direct high transfers $`F`$ (`e == 0`) are never double-counted or claimed across chain edges. Summing the $`K\ge 2`$ retained vertex inequalities yields $`N(T) > K B - A\sum (c_u - 2) - (K-1)B = B + 2A = 371/12500`$.
4. **Early pair ($`p_2 - p_3 \ge 127/2440`$) & telescoping synthesis:** Recomputed from scratch: $`f(j)=\min(1, 2/(j+1))`$ has exact drift $`-1/3`$ at $`j=1`$ and $`-(j-2)/[(j+1)(j+2)] \le 0`$ for all $`j\ge 2`$; the 4-step survival ratios are identically $`q = (2/3, 3/4, 6/7, 56/61)`$ with $`v_4 = 3/5`$, giving strictly increasing $`v_0 < v_1 < v_2 < v_3 < v_4`$ and $`v_1 = 108/305`$; hence $`H_\ell \le 197/305`$, $`H_c \ge 1 - (1/2)^2 = 3/4`$, and $`p_2 - p_3 = \frac12(H_c - H_\ell) \ge \frac12(\frac34 - \frac{197}{305}) = \frac{127}{2440}`$. Conditioned telescoping $`\sum_{v=i}^{j-1}\frac{c}{(v-1)v} = c(\frac1{i-1}-\frac1{j-1})`$ follows immediately by the tower property.

---

## 2. Independent Verification & Exact Enumeration Evidence

Standalone script (zero imports from `ft2-astra`): `round2/outputs/attack_astra_26/audit_claim26.py` in `exact/neon`.

Command run:
```sh
python3 -IB round2/outputs/attack_astra_26/audit_claim26.py
```

Observed output (exit 0, wall `35.43s`):
- **Exact constants & scalar barriers:**
  - `c = 371/12500 = 0.02968000` (`B = 17/1000`, `A = 634/100000`).
  - `p_2 - p_3 >= 127/2440 = 0.0520491803` (`v_1 = 108/305`).
  - `6,105` exact `P(d, j)` monotonicity rows passed; `B(2,1)=0.777202011`, `U(2,0)=0.610259660`, `L(2)=0.340723508`, `G(2)=0.189409485`.
- **Independent all-tree discharging recomputation:**
  - `30,900` high-donor simplex corners (`10 <= d <= 112`): all `> 0`, minimum `0.029129651179` at `(d=10, c=9, cost=0, e=0)`, plus `d >= 113` analytic inequalities and spot checks at `d = 113, 114, 200`.
  - Bare degree-2 chains (`k = 1, 2, 3`, low–low / low–high / high–high): max residuals `(0.010693898382, 0.016746143875, 0.011500000000)` at `k=1`, `(0.008932801958, 0.016904844577, 0.016928215518)` at `k=2`, `(-0.013583895045, -0.005611852426, 0.002360190192)` at `k=3` — all `< B = 0.017`.
  - `184,628` simultaneous low-endpoint profiles (`d = 2..9`): all `> 0`, matching per-degree minima (`d=2: 0.000005478682`, `d=3: 0.000020000000`, `d=4: 0.009920000000`, `d=5: 0.002020000000`, `d=6: 0.008360000000`, `d=7: 0.014356549140`, `d=8: 0.018727171760`, `d=9: 0.023728599884`).
- **Exact brute-force two-birth algebraic identity (all free trees $`n = 2\dots 9`$):**
  - `94` non-isomorphic free trees, `6,902` exact rational two-birth trajectories matched identically against an asymmetric rooted-isomorphism test functional.
- **Exact conditional two-birth inequality & skeleton ledger enumeration (ALL non-isomorphic free trees $`n = 2\dots 16`$, OEIS A000055 verified):**
  - **`32,507` free trees** (`32,492` multi-core) checked in exact `Fraction` arithmetic in `17.77s`; **0 failures**.
  - Global minimum over all `32,507` trees: `0.1894094851` at $`n = 2`$ (seed $`K_2`$, margin `+0.1597294851` above `371/12500`).
  - Global multi-core minimum over all `32,492` multi-core trees ($`n = 4\dots 16`$): `0.2740427539` at $`n = 4`$ ($`P_4`$, deg `(2,2,1,1)`, margin `+0.2443627539` above `371/12500`; at $`n = 16`$, `mc_min = 0.2789598015` on deg `(14,2,1^{14})`).
- **Adversarial large hostile trees ($`n`$ up to `461`):**
  - All 6 stress configurations (`hub3x113_on_3_3`, `hub4x115_on_4_4`, `hub113_3_2_chain1_d2term`, `mixed_stress_web`, `bare_path_P14`, `spider_200x2`) passed both `N_phys > 371/12500` and the exact contracted-skeleton ledger with strictly positive surplus (`+0.185` to `+1.376`).
- **Planted-error controls (7/7 caught):**
  1. Drop return-race rescue ($`R_d(2)=0`$): `d=2` terminal ledger drops to `-5.150346e-02`.
  2. Drop direct high-to-low transfer ($`F=0`$): `(3,3)` recipient drops to `-1.066000e-02`.
  3. Illegally credit a bare chain port (`e=2`) with direct transfer $`F`$: conservation imbalance `+3.560000e-03`.
  4. Omit serial path term $`w_u I(u)`$ in two-birth identity on $`T_2`$: exact mismatch `5.208803e-01`.
  5. Flip edge-term sign in two-birth identity on $`T_3`$: exact mismatch `8.741588e-03`.
  6. Truncate $`U(3,0)`$ downward instead of upward: unsound deficit `-1.623967e-13`.
  7. Tighten constant to `372/12500` (`A + 4e-5`): `d=2` terminal fails at `-3.452132e-05`.

---

## 3. Clean-Archive Hash Inventory & Replay

- Verified all **18/18** git blob SHA-1 hashes listed in `RELEASE.md` at `b450a310548c62adcfb37d66c0415d9dc3f7b3fc` (exact match, and zero diff on those 18 files at `212dc9f59b9f2c660ec71734826a5938daa6cab9`).
- Unpacked clean `git archive` snapshots of `b450a310548c62adcfb37d66c0415d9dc3f7b3fc` and `212dc9f59b9f2c660ec71734826a5938daa6cab9` into temporary directories and ran:
  - All **12/12** commands from `RELEASE.md` (`hostile/check.py`, `discharge/check.py` in all 3 modes, `mixed-stock/check.py`, `interfaces/check.py` default & `9`, `interfaces/low-end.py`, `edge-dual/check.py`, `preplacement-discharge-sol.py`, `preplacement-all-sol.py`, `preplacement-all-second.py`): **all 12 exited 0**.
  - `hub-ranking/check.py` (`PASS seed, stars, founder gap, 4753 first-send rows and 4960 label intervals`) and `early/check.py --horizon 4` (`p1=p2; p2-p3 >= 127/2440 > 0.052049180`): **both exited 0**.

No gaps were found in the zero-threshold hub-probability theorem or its pinned proofs and arithmetic. **VERDICT: HOLDS.**
