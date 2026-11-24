# Experimental Diagnosis of $`\mu \approx 0.569`$ in Hartle, Krapivsky, Redner & Zhang (arXiv:2604.01540)

**Author:** Fable's Gemini 4 Neon researcher. **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Claude lead reviewed; Neon independently reached the same diagnosis using Gemini. The analysis depends on the extracted Fig. 6 values and stated simulation and fitting protocols; see [Fable fact 14](../README.md#fable-fact-14).

## Executive Verdict

**Supported explanation: (a) Pure sampling noise of the protocol** — specifically, **independent 2400-tree ensembles at each grid size $`N`$**, combined with **non-self-averaging of the core size $`C`$** ($`\mathrm{CV}(C) \approx 0.51\text{–}0.60`$) and **unweighted nonlinear least squares in linear space** ($`\langle C \rangle = a N^b`$).

Key findings:
1. **Direct forensic extraction from the paper's Fig. 6 PDF (`2604.01540.pdf`, Form XObject 553, `CvsN.pdf`)** establishes all three protocol details that the caption omits:
   - Fig. 6 uses a **5-point decade grid** $`N \in \{10^2, 10^3, 10^4, 10^5, 10^6\}`$.
   - Unweighted nonlinear least squares in linear space ($`\langle C \rangle = a N^b`$) on their exact extracted marker points reproduces all three reported exponents: **IR `0.5695`** (reported `0.569`), **DAN `0.7732`** (reported `0.773`), **PAN `0.5456`** (reported `0.545`). By contrast, log-log OLS on the same points gives IR `0.5671`, DAN `0.7720`, and **PAN `0.5600`** (far from `0.545`).
   - Their extracted decade-to-decade local slopes fluctuate non-monotonically by $`+0.0103`$ (IR, jumping from `0.5606` on $`[10^4,10^5]`$ to `0.5709` on $`[10^5,10^6]`$) and $`+0.0057`$ (DAN), and drop to `0.5439` on $`[10^5,10^6]`$ for PAN (`-0.0080` below the true $`[10^5,10^6]`$ slope `0.5519`). On coupled (same-tree) ensembles of 2400 realizations these jumps would be **$`25\sigma\text{–}47\sigma`$ impossible**, proving Hartle et al. simulated **independent 2400-tree ensembles at each size $`N`$**.
2. **In the independent-grid-size + unweighted linear-space NLS protocol (60 replicates $`\times`$ 2400 trees/grid-size per model)**:
   - Because $`C(N)`$ is non-self-averaging ($`\mathrm{CV}(C) = 0.5114`$ for IR, $`0.5488`$ for DAN, $`0.5983`$ for PAN at $`N=10^6`$, from 336,000 pooled trees per model), every independent grid point $`\overline{C}(N)`$ carries an independent relative standard error of $`\mathrm{CV}/\sqrt{2400} \approx 1.04\%\text{–}1.22\%`$.
   - Unweighted linear-space NLS weights squared residuals by $`C(N)^2 \propto N^{2\mu}`$, making the fitted exponent $`b`$ dominated by the two largest points ($`10^5`$ and $`10^6`$) and inflating $`\mathrm{SD}(b)`$ to **`0.0039–0.0052` (IR)**, **`0.0072–0.0087` (DAN)**, and **`0.0054–0.0064` (PAN)** (in close agreement with exact linearized delta-method theory).
   - All three values reported in Fig. 6 are routine sampling draws:
     - **IR ($`b \ge 0.569`$)**: **16/60 (26.7%)** on 5-pt decade ($`z = +0.59`$), **18/60 (30.0%)** on 9-pt half-decade ($`z = +0.43`$), **17/60 (28.3%)** on 17-pt quarter-decade ($`z = +0.58`$).
     - **DAN ($`b \ge 0.773`$)**: **21/60 (35.0%)** on all three grids ($`z = +0.32\text{ to }+0.36`$).
     - **PAN ($`b \le 0.545`$)**: **4/60 (6.7%)** on 5-pt decade ($`z = -1.29`$), **8/60 (13.3%)** on 9-pt half-decade ($`z = -1.15`$), **3/60 (5.0%)** on 17-pt quarter-decade ($`z = -1.42`$), and **23%–28%** with an additive constant $`a N^b + c`$.
3. **Why (b) systematic finite-size bias is *not* the explanation for IR $`0.569`$ under linear-space NLS**:
   - From 192,000 pooled coupled trees, IR's finite-size bias under linear-space NLS is only **`+0.00012` to `+0.00016`** (1/25th of $`+0.0033`$). Its local slope has already converged to `0.565766 ± 0.000033` on $`[10^5, 10^6]`$.
   - Under log-log OLS, IR does have a systematic positive finite-size bias of `+0.00147` to `+0.00168` (driven by the $`[10^2,10^3]`$ slope `0.5713`), which is still only half of $`+0.0033`$; moreover, log-log OLS produces a huge **positive** bias (`+0.017` to `+0.019`) on PAN (`Mean(b) ≈ 0.568`), making `b <= 0.545` a $`-14.5\sigma`$ impossibility (0/60), ruling out log-log OLS as the paper's estimator.
   - For PAN, even linear-space NLS has a **positive** finite-size bias of `+0.0031` (due to the $`1/(N - \mathcal{N}) = N^{-1}(1 + O(N^{\mu-1}))`$ rejection denominator), opposite in sign to the paper's `-0.0047` error (`0.545`).
4. **Why coupled (same-tree) sampling cannot reproduce 0.569, 0.773, or 0.545 (80 replicates $`\times`$ 2400 trees per model)**:
   - On coupled trees, $`\mathrm{corr}(C(10^5), C(10^6)) = 0.9978`$ (IR), $`0.9997`$ (DAN), $`0.9982`$ (PAN): the frozen non-self-averaging amplitude cancels out of the slope, reducing $`\mathrm{SD}(b)`$ by **20x–60x** to `0.00010–0.00029`.
   - Under coupled sampling, **0/80 replicates (0.0%)** hit the paper's values for any model, grid, or fit method ($`z = +12.3\sigma`$ for IR, $`+17.5\sigma`$ for DAN, $`-27.5\sigma`$ for PAN on 5-pt linear NLS).

---

## 1. Files and Reproduction

All source files, raw replicate CSVs, pooled covariance CSVs, and logs live in `fit569/`:
- `sim.cpp` / `sim`: Fast C++20 simulator for IR, DAN, and PAN (exact $`O(1)`$ per step using a flat `leaf_parent` array + per-core singly-linked leaf list, Lemire exact unbiased integer sampling, `xoshiro256++` PRNG with injective SplitMix64 seeding, capped at $`\le 2`$ threads).
- `validate_exact.py`: Exact rational (`fractions.Fraction`) tree-shape enumerator from $`L_3`$ up to $`N = 8`$ for IR, DAN, and PAN, plus high-precision root solver for $`\mu_{\text{DAN}}`$ and $`\mu_{\text{PAN}}`$.
- `parse_fig6.py` & `cvsn_stream.txt`: Vector-stream extraction and re-fitting of Fig. 6 (`CvsN.pdf`, Object 553) from `lit/raw/2604.01540.pdf`.
- `analyze.py` & `analysis_output.txt`: Full replicate fitting (`linear_NLS`, `loglog_OLS`, `linear_NLS_plus_c` on 5-pt decade, 9-pt half-decade, and 17-pt quarter-decade grids), pooled local slopes with exact delta-method SEs, and $`\mathrm{CV}(C)`$ at $`N=10^6`$.
- `verify_sampling_theory.py`: Exact first-order linearized (delta-method) error propagation verifying empirical `SD(b)` and hit probabilities in independent-grid mode.
- `run.sh`: One-command reproduction script (`./run.sh`). Uses at most 2 CPU cores at a time (~65 minutes total wall time).

### Grid Definitions
All simulations record $`C(N)`$ on the 17-point quarter-decade grid $`N_k = \mathrm{round}(10^{2 + k/4})`$ for $`k = 0,\dots,16`$:
```
100, 178, 316, 562, 1000, 1778, 3162, 5623, 10000, 17783, 31623, 56234, 100000, 177828, 316228, 562341, 1000000
```
- **5-pt decade (`5pt_decade`)**: indices `{0, 4, 8, 12, 16}` ($`10^2, 10^3, 10^4, 10^5, 10^6`$) — matches Fig. 6 exactly.
- **9-pt half-decade (`9pt_halfdec`)**: even indices `{0, 2, 4, 6, 8, 10, 12, 14, 16}`.
- **17-pt quarter-decade (`17pt_qtrdec`)**: all indices `{0..16}`.

---

## 2. Validation Against Exact Enumeration and Reference Implementations

### 2.1 Analytical Roots and Exact Small-$`N`$ Enumeration ($`N = 3 \dots 8`$)

Command:
```bash
.venv/bin/python fit569/validate_exact.py
```
Output:
```text
Analytical roots:
  mu_DAN = 0.771192302407 (paper: 0.771192)
  mu_PAN = 0.549734972504 (paper: 0.549735)

=== Exact enumeration for IR (starting from L_3) ===
  N=3 (shapes= 1): E[C]=1 = 1.000000000000, Var[C]=0.000000000000, S_N=1, Q_N=1
  N=4 (shapes= 2): E[C]=4/3 = 1.333333333333, Var[C]=0.222222222222, S_N=2/3, Q_N=1/3
  N=5 (shapes= 3): E[C]=19/12 = 1.583333333333, Var[C]=0.409722222222, S_N=1/2, Q_N=1/12
  N=6 (shapes= 6): E[C]=647/360 = 1.797222222222, Var[C]=0.589436728395, S_N=2/5, Q_N=1/60
  N=7 (shapes=11): E[C]=17189/8640 = 1.989467592593, Var[C]=0.768214685678, S_N=1/3, Q_N=1/360
  N=8 (shapes=23): E[C]=4717423/2177280 = 2.166658858760, Var[C]=0.948313798195, S_N=2/7, Q_N=1/2520
  [OK] Matches prompt exact IR values E[C_4]=4/3, E[C_5]=19/12, E[C_6]=647/360

=== Exact enumeration for DAN (starting from L_3) ===
  N=3 (shapes= 1): E[C]=1 = 1.000000000000, Var[C]=0.000000000000, S_N=1, Q_N=1
  N=4 (shapes= 2): E[C]=4/3 = 1.333333333333, Var[C]=0.222222222222, S_N=2/3, Q_N=1/3
  N=5 (shapes= 3): E[C]=5/3 = 1.666666666667, Var[C]=0.555555555556, S_N=1/2, Q_N=1/6
  N=6 (shapes= 6): E[C]=59/30 = 1.966666666667, Var[C]=0.898888888889, S_N=2/5, Q_N=1/15
  N=7 (shapes=11): E[C]=9/4 = 2.250000000000, Var[C]=1.265277777778, S_N=1/3, Q_N=1/45
  N=8 (shapes=23): E[C]=794/315 = 2.520634920635, Var[C]=1.657510707987, S_N=2/7, Q_N=2/315
  [OK] Matches exact DAN S_N=2/(N-1) and Q_N=Eq.(A13) for N=3..8

=== Exact enumeration for PAN (starting from L_3) ===
  N=3 (shapes= 1): E[C]=1 = 1.000000000000, Var[C]=0.000000000000, S_N=1, Q_N=1
  N=4 (shapes= 2): E[C]=4/3 = 1.333333333333, Var[C]=0.222222222222, S_N=2/3, Q_N=1/3
  N=5 (shapes= 3): E[C]=5/3 = 1.666666666667, Var[C]=0.555555555556, S_N=1/2, Q_N=1/6
  N=6 (shapes= 6): E[C]=119/60 = 1.983333333333, Var[C]=0.949722222222, S_N=2/5, Q_N=1/12
  N=7 (shapes=11): E[C]=457/200 = 2.285000000000, Var[C]=1.390441666667, S_N=1/3, Q_N=1/24
  N=8 (shapes=23): E[C]=972751/378000 = 2.573415343915, Var[C]=1.869250398918, S_N=2/7, Q_N=1/48
  [OK] Matches exact PAN S_N=2/(N-1) and Q_N=Eq.(A11) for N=3..8
```

### 2.2 Fast C++ Simulator Validation ($`10^7`$ Trees at $`N=4\dots 8`$, $`2\times 10^5`$ Fast vs. Naive at $`N=100`$)

Command:
```bash
g++ -O3 -march=native -std=c++20 -pthread sim.cpp -o sim && ./sim validate 10000000 200000
```
Output:
```text
=== Validating fast simulator: IR ===
  N=4: sim=1.33311370 exact=1.33333333 se=0.00014905 z=-1.47  S_N=0.666886 (exact 0.666667)  Q_N=0.333114
  N=5: sim=1.58293460 exact=1.58333333 se=0.00020236 z=-1.97  S_N=0.500254 (exact 0.500000)  Q_N=0.083188
  N=6: sim=1.79685440 exact=1.79722222 se=0.00024273 z=-1.52  S_N=0.400084 (exact 0.400000)  Q_N=0.016705
  N=7: sim=1.98916190 exact=1.98946759 se=0.00027706 z=-1.10  S_N=0.333274 (exact 0.333333)  Q_N=0.002787
  N=8: sim=2.16626540 exact=2.16665886 se=0.00030787 z=-1.28  S_N=0.285714 (exact 0.285714)  Q_N=0.000394
  N=100 cross-check (T=200000): fast=9.613655 (sd=4.778776) vs naive=9.614530 (sd=4.770168), z=-0.06

=== Validating fast simulator: DAN ===
  N=4: sim=1.33357770 exact=1.33333333 se=0.00014910 z=1.64  S_N=0.666422 (exact 0.666667)  Q_N=0.333578
  N=5: sim=1.66706990 exact=1.66666667 se=0.00023576 z=1.71  S_N=0.499793 (exact 0.500000)  Q_N=0.166863
  N=6: sim=1.96726690 exact=1.96666667 se=0.00029989 z=2.00  S_N=0.399769 (exact 0.400000)  Q_N=0.066789
  N=7: sim=2.25079990 exact=2.25000000 se=0.00035573 z=2.25  S_N=0.333043 (exact 0.333333)  Q_N=0.022189
  N=8: sim=2.52169260 exact=2.52063492 se=0.00040721 z=2.60  S_N=0.285476 (exact 0.285714)  Q_N=0.006330
  N=100 cross-check (T=200000): fast=18.588745 (sd=10.176868) vs naive=18.532120 (sd=10.167287), z=1.76

=== Validating fast simulator: PAN ===
  N=4: sim=1.33306700 exact=1.33333333 se=0.00014904 z=-1.79  S_N=0.666933 (exact 0.666667)  Q_N=0.333067
  N=5: sim=1.66629780 exact=1.66666667 se=0.00023560 z=-1.57  S_N=0.500061 (exact 0.500000)  Q_N=0.166359
  N=6: sim=1.98275010 exact=1.98333333 se=0.00030800 z=-1.89  S_N=0.399983 (exact 0.400000)  Q_N=0.083113
  N=7: sim=2.28415290 exact=2.28500000 se=0.00037265 z=-2.27  S_N=0.333328 (exact 0.333333)  Q_N=0.041486
  N=8: sim=2.57261910 exact=2.57341534 se=0.00043207 z=-1.84  S_N=0.285598 (exact 0.285714)  Q_N=0.020754
  N=100 cross-check (T=200000): fast=16.380180 (sd=9.387729) vs naive=16.349090 (sd=9.364482), z=1.05
```

Additionally, because DAN's rate equations (Eqs. 9–11 of arXiv:2604.01540) are linear in expectations, we integrated them exactly from $`L_3`$ up to $`N = 10^4`$, obtaining exact $`\langle C_{\text{DAN}}(100) \rangle = 18.568047`$, $`\langle C_{\text{DAN}}(10^3) \rangle = 109.900990`$, $`\langle C_{\text{DAN}}(10^4) \rangle = 649.052380`$, in full agreement with our 192,000-tree pooled simulation (`18.6050 ± 0.0232`, `110.1314 ± 0.1377`, `650.3344 ± 0.8138`; all within $`1.7\sigma`$).

---

## 3. Forensic Analysis of Fig. 6 from `2604.01540.pdf`

Commands:
```bash
mutool show -b lit/raw/2604.01540.pdf 553 > fit569/cvsn_stream.txt
.venv/bin/python fit569/parse_fig6.py
```
Output:
```text
x_per_dec = 50.8030317500, y_per_dec = 32.1219916667

=== PAN (red /M0) ===
  Extracted <C> at N=[1e2, 1e3, 1e4, 1e5, 1e6]:
    [20.6263, 78.185, 286.9563, 1026.4254, 3591.3813]
  Adjacent decade local slopes (1e2-1e3, 1e3-1e4, 1e4-1e5, 1e5-1e6):
    [0.578703, 0.564692, 0.553512, 0.543934]
  Plotted straight line in Fig. 6:  slope = 0.549735, a = 1.800000
  Fit (i) linear NLS a*N^b:         b = 0.545618, a = 1.912603
  Fit (ii) log-log OLS:             b = 0.559989, a = 1.608552
  Fit (iii) linear NLS a*N^b + c:   b = 0.542876, a = 1.989390, c = -5.6844
  Fit (iv) rel-weighted NLS:        b = 0.560009, a = 1.607049

=== DAN (blue /M1) ===
  Extracted <C> at N=[1e2, 1e3, 1e4, 1e5, 1e6]:
    [23.0767, 137.8883, 816.195, 4781.5764, 28383.7129]
  Adjacent decade local slopes (1e2-1e3, 1e3-1e4, 1e4-1e5, 1e5-1e6):
    [0.776354, 0.772266, 0.767777, 0.773498]
  Plotted straight line in Fig. 6:  slope = 0.771192, a = 0.650000
  Fit (i) linear NLS a*N^b:         b = 0.773187, a = 0.651544
  Fit (ii) log-log OLS:             b = 0.771984, a = 0.662968
  Fit (iii) linear NLS a*N^b + c:   b = 0.773583, a = 0.647898, c = 3.8605
  Fit (iv) rel-weighted NLS:        b = 0.771987, a = 0.662925

=== IR (green /M2) ===
  Extracted <C> at N=[1e2, 1e3, 1e4, 1e5, 1e6]:
    [9.6479, 36.0513, 133.2113, 484.3083, 1802.9254]
  Adjacent decade local slopes (1e2-1e3, 1e3-1e4, 1e4-1e5, 1e5-1e6):
    [0.572487, 0.567621, 0.560581, 0.570856]
  Plotted straight line in Fig. 6:  slope = 0.566000, a = 0.710000
  Fit (i) linear NLS a*N^b:         b = 0.569521, a = 0.689926
  Fit (ii) log-log OLS:             b = 0.567129, a = 0.712661
  Fit (iii) linear NLS a*N^b + c:   b = 0.570810, a = 0.677312, c = 1.2513
  Fit (iv) rel-weighted NLS:        b = 0.567132, a = 0.712598
```

Observations from the exact Fig. 6 data:
- **Exact IR match**: Fig. 6's IR values `[9.648, 36.051, 133.211, 484.308, 1802.925]` match our $`L_3`$-initialized IR pooled means `[9.608, 35.799, 132.141, 486.644, 1790.514]` to within $`\pm 0.8\%`$ at every decade (`N=10^5` fluctuated $`-0.48\%`$ low and `N=10^6` fluctuated $`+0.69\%`$ high, driving the $`[10^5,10^6]`$ local slope up to `0.5709` and the linear-space NLS fit to `0.5695`!).
- **Exact estimator identification**: Unweighted linear-space NLS (`a * N^b`) reproduces all three caption values (`0.5695`, `0.7732`, `0.5456` vs. `0.569`, `0.773`, `0.545`), whereas log-log OLS fails badly on PAN (`0.5600` vs. `0.545`). (The straight guide lines drawn in Fig. 6 are theoretical reference lines `0.65 N^0.771192`, `1.80 N^0.549735`, and `0.71 N^0.566`, not the NLS fits.)
- **Secondary observation on DAN/PAN amplitude in Fig. 6**: While IR matches $`L_3`$ to 0.3%, Fig. 6's DAN and PAN points sit roughly $`\sim 1.25\times`$ above our $`L_3`$-initialized values (and above exact $`L_3`$ DAN rate equations Eqs. 9–11, which give `22625` vs. Fig. 6's `28384` at $`N=10^6`$), and Fig. 6's PAN has slightly less curvature on $`[10^2,10^3]`$ (`0.5787` vs. `0.6056` for $`L_3`$). A constant multiplicative prefactor shift with unchanged asymptotic exponent $`\mu`$ indicates DAN/PAN in Fig. 6 were likely initialized from a slightly larger seed tree (e.g., $`C_0 = 2`$ such as $`L_4`$ or a 5-node tree) rather than $`L_3`$. This does not affect IR (which matches $`L_3`$ exactly) and only slightly *increases* the probability of $`b_{\text{PAN}} \le 0.545`$ by reducing PAN's positive small-$`N`$ bias.

---

## 4. Production Simulation Runs and Replicate Fit Results

### Commands Run
```bash
# Coupled mode: 80 replicates x 2400 trees per model (192,000 trees/model)
./sim run --model IR  --mode coupled --replicates 80 --trees 2400 --threads 2 --seed 1001 --out ir_coupled
./sim run --model DAN --mode coupled --replicates 80 --trees 2400 --threads 2 --seed 2001 --out dan_coupled
./sim run --model PAN --mode coupled --replicates 80 --trees 2400 --threads 2 --seed 3001 --out pan_coupled

# Independent-grid-size mode: 60 replicates x 2400 fresh trees at each of 17 grid sizes
./sim run --model IR  --mode indep   --replicates 60 --trees 2400 --threads 2 --seed 4001 --out ir_indep
./sim run --model DAN --mode indep   --replicates 60 --trees 2400 --threads 2 --seed 5001 --out dan_indep
./sim run --model PAN --mode indep   --replicates 60 --trees 2400 --threads 2 --seed 6001 --out pan_indep
```
Wall-clock runtimes (2 threads on AMD Ryzen 7 5800X):
- Coupled (80 reps each): IR `579.79 s`, DAN `449.86 s`, PAN `545.58 s` (total `26.3 min`).
- Independent (60 reps each, 17 grid sizes): IR `965.63 s`, DAN `713.89 s`, PAN `651.36 s` (total `38.8 min`).

### 4.1 Replicate Fit Summary Tables

Command:
```bash
.venv/bin/python fit569/analyze.py fit569 | tee fit569/analysis_output.txt
```

Real pasted output (all models, modes, grids, and fit variants):

```text
================================================================================================
FIT REPLICATE SUMMARY (2400 realizations per replicate)
================================================================================================

### Model: IR | Mode: coupled | Replicates: R=80 x T=2400 | Ref mu=0.565723 | Target: b >= 0.569
Grid         Fit Method            Mean(b)     SD(b)       Bias  SE(mean)   Frac(b >= 0.569)  z_target    Min(b)    Max(b)
--------------------------------------------------------------------------------------------------------------------------
5pt_decade   linear_NLS           0.565844  0.000257  +0.000121  0.000029   0/80  (  0.0%)    +12.29  0.564977  0.566281
5pt_decade   loglog_OLS           0.567407  0.000528  +0.001684  0.000059   0/80  (  0.0%)     +3.02  0.566148  0.568419
5pt_decade   linear_NLS_plus_c    0.565656  0.000288  -0.000067  0.000032   0/80  (  0.0%)    +11.62  0.564614  0.566150
  [Pooled 192,000 trees 5pt_decade] linear_NLS=0.565844 (bias +0.000121), loglog_OLS=0.567407 (bias +0.001684), NLS+c=0.565657 (bias -0.000066)
9pt_halfdec  linear_NLS           0.565870  0.000260  +0.000147  0.000029   0/80  (  0.0%)    +12.03  0.565130  0.566447
9pt_halfdec  loglog_OLS           0.567263  0.000483  +0.001540  0.000054   0/80  (  0.0%)     +3.60  0.566346  0.568331
9pt_halfdec  linear_NLS_plus_c    0.565696  0.000305  -0.000027  0.000034   0/80  (  0.0%)    +10.83  0.564790  0.566289
  [Pooled 192,000 trees 9pt_halfdec] linear_NLS=0.565870 (bias +0.000147), loglog_OLS=0.567263 (bias +0.001540), NLS+c=0.565696 (bias -0.000027)
17pt_qtrdec  linear_NLS           0.565887  0.000255  +0.000164  0.000028   0/80  (  0.0%)    +12.23  0.565178  0.566419
17pt_qtrdec  loglog_OLS           0.567189  0.000462  +0.001466  0.000052   0/80  (  0.0%)     +3.92  0.566188  0.568290
17pt_qtrdec  linear_NLS_plus_c    0.565703  0.000293  -0.000020  0.000033   0/80  (  0.0%)    +11.23  0.564807  0.566349
  [Pooled 192,000 trees 17pt_qtrdec] linear_NLS=0.565887 (bias +0.000164), loglog_OLS=0.567189 (bias +0.001466), NLS+c=0.565703 (bias -0.000020)

### Model: IR | Mode: indep | Replicates: R=60 x T=2400 | Ref mu=0.565723 | Target: b >= 0.569
Grid         Fit Method            Mean(b)     SD(b)       Bias  SE(mean)   Frac(b >= 0.569)  z_target    Min(b)    Max(b)
--------------------------------------------------------------------------------------------------------------------------
5pt_decade   linear_NLS           0.566450  0.004295  +0.000727  0.000554  16/60  ( 26.7%)     +0.59  0.555555  0.577106
5pt_decade   loglog_OLS           0.567680  0.001333  +0.001957  0.000172   8/60  ( 13.3%)     +0.99  0.564456  0.570379
5pt_decade   linear_NLS_plus_c    0.566335  0.005548  +0.000612  0.000716  19/60  ( 31.7%)     +0.48  0.553377  0.579546
  [Pooled 144,000 trees 5pt_decade] linear_NLS=0.566444 (bias +0.000721), loglog_OLS=0.567680 (bias +0.001957), NLS+c=0.566329 (bias +0.000606)
9pt_halfdec  linear_NLS           0.566789  0.005171  +0.001066  0.000668  18/60  ( 30.0%)     +0.43  0.551485  0.578688
9pt_halfdec  loglog_OLS           0.567391  0.001178  +0.001668  0.000152   5/60  (  8.3%)     +1.37  0.564452  0.570132
9pt_halfdec  linear_NLS_plus_c    0.566942  0.007330  +0.001219  0.000946  24/60  ( 40.0%)     +0.28  0.546262  0.583380
  [Pooled 144,000 trees 9pt_halfdec] linear_NLS=0.566760 (bias +0.001037), loglog_OLS=0.567391 (bias +0.001668), NLS+c=0.566900 (bias +0.001177)
17pt_qtrdec  linear_NLS           0.566731  0.003925  +0.001008  0.000507  17/60  ( 28.3%)     +0.58  0.555807  0.577124
17pt_qtrdec  loglog_OLS           0.567213  0.000852  +0.001490  0.000110   1/60  (  1.7%)     +2.10  0.565028  0.569222
17pt_qtrdec  linear_NLS_plus_c    0.566955  0.005681  +0.001232  0.000733  21/60  ( 35.0%)     +0.36  0.552119  0.581148
  [Pooled 144,000 trees 17pt_qtrdec] linear_NLS=0.566715 (bias +0.000992), loglog_OLS=0.567213 (bias +0.001490), NLS+c=0.566931 (bias +0.001208)

### Model: DAN | Mode: coupled | Replicates: R=80 x T=2400 | Ref mu=0.771192 | Target: b >= 0.773
Grid         Fit Method            Mean(b)     SD(b)       Bias  SE(mean)   Frac(b >= 0.773)  z_target    Min(b)    Max(b)
--------------------------------------------------------------------------------------------------------------------------
5pt_decade   linear_NLS           0.771161  0.000105  -0.000032  0.000012   0/80  (  0.0%)    +17.46  0.770852  0.771402
5pt_decade   loglog_OLS           0.771408  0.000346  +0.000216  0.000039   0/80  (  0.0%)     +4.61  0.770724  0.772379
5pt_decade   linear_NLS_plus_c    0.771156  0.000118  -0.000036  0.000013   0/80  (  0.0%)    +15.68  0.770815  0.771410
  [Pooled 192,000 trees 5pt_decade] linear_NLS=0.771161 (bias -0.000032), loglog_OLS=0.771409 (bias +0.000216), NLS+c=0.771156 (bias -0.000036)
9pt_halfdec  linear_NLS           0.771169  0.000096  -0.000023  0.000011   0/80  (  0.0%)    +19.12  0.770882  0.771429
9pt_halfdec  loglog_OLS           0.771370  0.000319  +0.000177  0.000036   0/80  (  0.0%)     +5.12  0.770701  0.772209
9pt_halfdec  linear_NLS_plus_c    0.771169  0.000107  -0.000023  0.000012   0/80  (  0.0%)    +17.11  0.770864  0.771414
  [Pooled 192,000 trees 9pt_halfdec] linear_NLS=0.771169 (bias -0.000023), loglog_OLS=0.771371 (bias +0.000178), NLS+c=0.771169 (bias -0.000023)
17pt_qtrdec  linear_NLS           0.771170  0.000091  -0.000022  0.000010   0/80  (  0.0%)    +20.07  0.770906  0.771422
17pt_qtrdec  loglog_OLS           0.771348  0.000300  +0.000156  0.000034   0/80  (  0.0%)     +5.51  0.770688  0.772152
17pt_qtrdec  linear_NLS_plus_c    0.771170  0.000101  -0.000022  0.000011   0/80  (  0.0%)    +18.21  0.770895  0.771399
  [Pooled 192,000 trees 17pt_qtrdec] linear_NLS=0.771170 (bias -0.000022), loglog_OLS=0.771349 (bias +0.000156), NLS+c=0.771170 (bias -0.000022)

### Model: DAN | Mode: indep | Replicates: R=60 x T=2400 | Ref mu=0.771192 | Target: b >= 0.773
Grid         Fit Method            Mean(b)     SD(b)       Bias  SE(mean)   Frac(b >= 0.773)  z_target    Min(b)    Max(b)
--------------------------------------------------------------------------------------------------------------------------
5pt_decade   linear_NLS           0.770686  0.007182  -0.000506  0.000927  21/60  ( 35.0%)     +0.32  0.757461  0.788849
5pt_decade   loglog_OLS           0.771241  0.001443  +0.000048  0.000186   6/60  ( 10.0%)     +1.22  0.767415  0.773711
5pt_decade   linear_NLS_plus_c    0.770607  0.008042  -0.000585  0.001038  22/60  ( 36.7%)     +0.30  0.756031  0.790739
  [Pooled 144,000 trees 5pt_decade] linear_NLS=0.770680 (bias -0.000512), loglog_OLS=0.771239 (bias +0.000047), NLS+c=0.770606 (bias -0.000586)
9pt_halfdec  linear_NLS           0.769855  0.008665  -0.001338  0.001119  21/60  ( 35.0%)     +0.36  0.753980  0.795355
9pt_halfdec  loglog_OLS           0.771298  0.001195  +0.000106  0.000154   2/60  (  3.3%)     +1.42  0.767853  0.773835
9pt_halfdec  linear_NLS_plus_c    0.769568  0.010305  -0.001624  0.001330  23/60  ( 38.3%)     +0.33  0.750678  0.799748
  [Pooled 144,000 trees 9pt_halfdec] linear_NLS=0.769817 (bias -0.001376), loglog_OLS=0.771298 (bias +0.000106), NLS+c=0.769526 (bias -0.001666)
17pt_qtrdec  linear_NLS           0.770293  0.007956  -0.000900  0.001027  21/60  ( 35.0%)     +0.34  0.754389  0.795982
17pt_qtrdec  loglog_OLS           0.771327  0.000875  +0.000134  0.000113   1/60  (  1.7%)     +1.91  0.769429  0.773922
17pt_qtrdec  linear_NLS_plus_c    0.770089  0.009735  -0.001103  0.001257  21/60  ( 35.0%)     +0.30  0.750762  0.801064
  [Pooled 144,000 trees 17pt_qtrdec] linear_NLS=0.770249 (bias -0.000943), loglog_OLS=0.771328 (bias +0.000136), NLS+c=0.770036 (bias -0.001156)

### Model: PAN | Mode: coupled | Replicates: R=80 x T=2400 | Ref mu=0.549735 | Target: b <= 0.545
Grid         Fit Method            Mean(b)     SD(b)       Bias  SE(mean)   Frac(b <= 0.545)  z_target    Min(b)    Max(b)
--------------------------------------------------------------------------------------------------------------------------
5pt_decade   linear_NLS           0.552833  0.000285  +0.003098  0.000032   0/80  (  0.0%)    -27.50  0.552188  0.553586
5pt_decade   loglog_OLS           0.568958  0.000665  +0.019223  0.000074   0/80  (  0.0%)    -36.00  0.567396  0.570276
5pt_decade   linear_NLS_plus_c    0.550741  0.000322  +0.001006  0.000036   0/80  (  0.0%)    -17.81  0.550074  0.551675
  [Pooled 192,000 trees 5pt_decade] linear_NLS=0.552834 (bias +0.003099), loglog_OLS=0.568959 (bias +0.019224), NLS+c=0.550741 (bias +0.001006)
9pt_halfdec  linear_NLS           0.552775  0.000268  +0.003040  0.000030   0/80  (  0.0%)    -29.05  0.552176  0.553516
9pt_halfdec  loglog_OLS           0.567670  0.000640  +0.017935  0.000072   0/80  (  0.0%)    -35.44  0.566155  0.568950
9pt_halfdec  linear_NLS_plus_c    0.550655  0.000301  +0.000920  0.000034   0/80  (  0.0%)    -18.82  0.550058  0.551593
  [Pooled 192,000 trees 9pt_halfdec] linear_NLS=0.552776 (bias +0.003041), loglog_OLS=0.567672 (bias +0.017937), NLS+c=0.550655 (bias +0.000920)
17pt_qtrdec  linear_NLS           0.552936  0.000264  +0.003201  0.000029   0/80  (  0.0%)    -30.08  0.552294  0.553590
17pt_qtrdec  loglog_OLS           0.566912  0.000634  +0.017177  0.000071   0/80  (  0.0%)    -34.57  0.565380  0.568161
17pt_qtrdec  linear_NLS_plus_c    0.550656  0.000292  +0.000921  0.000033   0/80  (  0.0%)    -19.39  0.550009  0.551509
  [Pooled 192,000 trees 17pt_qtrdec] linear_NLS=0.552937 (bias +0.003202), loglog_OLS=0.566914 (bias +0.017179), NLS+c=0.550657 (bias +0.000922)

### Model: PAN | Mode: indep | Replicates: R=60 x T=2400 | Ref mu=0.549735 | Target: b <= 0.545
Grid         Fit Method            Mean(b)     SD(b)       Bias  SE(mean)   Frac(b <= 0.545)  z_target    Min(b)    Max(b)
--------------------------------------------------------------------------------------------------------------------------
5pt_decade   linear_NLS           0.552268  0.005653  +0.002533  0.000730   4/60  (  6.7%)     -1.29  0.540369  0.564807
5pt_decade   loglog_OLS           0.568732  0.001639  +0.018997  0.000212   0/60  (  0.0%)    -14.48  0.564846  0.571977
5pt_decade   linear_NLS_plus_c    0.550142  0.007140  +0.000407  0.000922  14/60  ( 23.3%)     -0.72  0.534800  0.566665
  [Pooled 144,000 trees 5pt_decade] linear_NLS=0.552264 (bias +0.002529), loglog_OLS=0.568735 (bias +0.019000), NLS+c=0.550147 (bias +0.000412)
9pt_halfdec  linear_NLS           0.552336  0.006362  +0.002601  0.000821   8/60  ( 13.3%)     -1.15  0.541119  0.566492
9pt_halfdec  loglog_OLS           0.567451  0.001442  +0.017716  0.000186   0/60  (  0.0%)    -15.57  0.563914  0.570261
9pt_halfdec  linear_NLS_plus_c    0.550311  0.008781  +0.000576  0.001134  17/60  ( 28.3%)     -0.60  0.535152  0.569273
  [Pooled 144,000 trees 9pt_halfdec] linear_NLS=0.552306 (bias +0.002571), loglog_OLS=0.567452 (bias +0.017717), NLS+c=0.550275 (bias +0.000541)
17pt_qtrdec  linear_NLS           0.552658  0.005397  +0.002923  0.000697   3/60  (  5.0%)     -1.42  0.540847  0.566149
17pt_qtrdec  loglog_OLS           0.566673  0.001074  +0.016938  0.000139   0/60  (  0.0%)    -20.18  0.563701  0.568869
17pt_qtrdec  linear_NLS_plus_c    0.550409  0.007927  +0.000674  0.001023  17/60  ( 28.3%)     -0.68  0.533895  0.569625
  [Pooled 144,000 trees 17pt_qtrdec] linear_NLS=0.552623 (bias +0.002888), loglog_OLS=0.566675 (bias +0.016940), NLS+c=0.550356 (bias +0.000621)
```

### 4.2 Analytical Cross-Check of Independent-Grid Sampling Variance

Command:
```bash
.venv/bin/python fit569/verify_sampling_theory.py
```
Output:
```text
Exact linearized (delta-method) predictions for INDEPENDENT 2400-tree replicates
(using 192,000-tree coupled pooled <C(N)> and CV(N) as ground truth):

=== IR (ref=0.565723, target b >= 0.569, rel_SE(10^6)=1.043%) ===
  5pt_decade   linear_NLS        : b_true_fit=0.565844 (bias +0.000121), SD_theory=0.005130, z=+0.62, P_hit= 26.9%
  5pt_decade   loglog_OLS        : b_true_fit=0.567407 (bias +0.001684), SD_theory=0.001414, z=+1.13, P_hit= 13.0%
  5pt_decade   linear_NLS_plus_c : b_true_fit=0.565657 (bias -0.000066), SD_theory=0.006571, z=+0.51, P_hit= 30.5%
  9pt_halfdec  linear_NLS        : b_true_fit=0.565870 (bias +0.000147), SD_theory=0.005556, z=+0.56, P_hit= 28.7%
  9pt_halfdec  loglog_OLS        : b_true_fit=0.567263 (bias +0.001540), SD_theory=0.001156, z=+1.50, P_hit=  6.6%
  9pt_halfdec  linear_NLS_plus_c : b_true_fit=0.565696 (bias -0.000027), SD_theory=0.007710, z=+0.43, P_hit= 33.4%
  17pt_qtrdec  linear_NLS        : b_true_fit=0.565887 (bias +0.000164), SD_theory=0.004530, z=+0.69, P_hit= 24.6%
  17pt_qtrdec  loglog_OLS        : b_true_fit=0.567189 (bias +0.001466), SD_theory=0.000887, z=+2.04, P_hit=  2.1%
  17pt_qtrdec  linear_NLS_plus_c : b_true_fit=0.565703 (bias -0.000020), SD_theory=0.006559, z=+0.50, P_hit= 30.8%

=== DAN (ref=0.771192, target b >= 0.773, rel_SE(10^6)=1.119%) ===
  5pt_decade   linear_NLS        : b_true_fit=0.771161 (bias -0.000032), SD_theory=0.006308, z=+0.29, P_hit= 38.5%
  5pt_decade   loglog_OLS        : b_true_fit=0.771409 (bias +0.000216), SD_theory=0.001535, z=+1.04, P_hit= 15.0%
  5pt_decade   linear_NLS_plus_c : b_true_fit=0.771156 (bias -0.000036), SD_theory=0.007117, z=+0.26, P_hit= 39.8%
  9pt_halfdec  linear_NLS        : b_true_fit=0.771169 (bias -0.000023), SD_theory=0.008232, z=+0.22, P_hit= 41.2%
  9pt_halfdec  loglog_OLS        : b_true_fit=0.771371 (bias +0.000178), SD_theory=0.001254, z=+1.30, P_hit=  9.7%
  9pt_halfdec  linear_NLS_plus_c : b_true_fit=0.771169 (bias -0.000023), SD_theory=0.009878, z=+0.19, P_hit= 42.6%
  17pt_qtrdec  linear_NLS        : b_true_fit=0.771170 (bias -0.000022), SD_theory=0.007403, z=+0.25, P_hit= 40.2%
  17pt_qtrdec  loglog_OLS        : b_true_fit=0.771349 (bias +0.000156), SD_theory=0.000962, z=+1.72, P_hit=  4.3%
  17pt_qtrdec  linear_NLS_plus_c : b_true_fit=0.771170 (bias -0.000022), SD_theory=0.009171, z=+0.20, P_hit= 42.1%

=== PAN (ref=0.549735, target b <= 0.545, rel_SE(10^6)=1.223%) ===
  5pt_decade   linear_NLS        : b_true_fit=0.552834 (bias +0.003099), SD_theory=0.005932, z=-1.32, P_hit=  9.3%
  5pt_decade   loglog_OLS        : b_true_fit=0.568959 (bias +0.019224), SD_theory=0.001649, z=-14.53, P_hit=  0.0%
  5pt_decade   linear_NLS_plus_c : b_true_fit=0.550741 (bias +0.001006), SD_theory=0.007671, z=-0.75, P_hit= 22.7%
  9pt_halfdec  linear_NLS        : b_true_fit=0.552776 (bias +0.003041), SD_theory=0.006341, z=-1.23, P_hit= 11.0%
  9pt_halfdec  loglog_OLS        : b_true_fit=0.567672 (bias +0.017937), SD_theory=0.001347, z=-16.83, P_hit=  0.0%
  9pt_halfdec  linear_NLS_plus_c : b_true_fit=0.550655 (bias +0.000920), SD_theory=0.008896, z=-0.64, P_hit= 26.2%
  17pt_qtrdec  linear_NLS        : b_true_fit=0.552937 (bias +0.003202), SD_theory=0.005140, z=-1.54, P_hit=  6.1%
  17pt_qtrdec  loglog_OLS        : b_true_fit=0.566914 (bias +0.017179), SD_theory=0.001033, z=-21.21, P_hit=  0.0%
  17pt_qtrdec  linear_NLS_plus_c : b_true_fit=0.550657 (bias +0.000922), SD_theory=0.007527, z=-0.75, P_hit= 22.6%
```

Note also why `SD(b)` under unweighted `linear_NLS` actually *increases* when moving from 5-pt decade to 9-pt half-decade (e.g., `0.0063 -> 0.0082` in DAN) before decreasing slightly at 17-pt quarter-decade: adding $`N = 10^{5.5}`$ introduces a high-weight point only 0.5 decade below $`10^6`$, shortening the effective $`\log N`$ lever arm in unweighted linear-space least squares.

---

## 5. Pooled Coefficient of Variation at $`N = 10^6`$, Cross-Scale Correlations, and Local Slopes

### 5.1 Coefficient of Variation $`\mathrm{CV}(C)`$ at $`N = 10^6`$ (336,000 Pooled Independent Trees per Model)

| Model | Coupled ($`M=192{,}000`$) | Indep ($`M=144{,}000`$) | Grand Pool ($`M=336{,}000`$) | $`\mathrm{CV}(C_{10^6})`$ | $`\mathrm{corr}(C_{10^5}, C_{10^6})`$ |
|---|---|---|---|---|---|
| **IR** | $`1790.51 \pm 2.09`$ ($`\mathrm{CV}=0.51082`$) | $`1795.26 \pm 2.42`$ ($`\mathrm{CV}=0.51215`$) | $`1792.55 \pm 1.58`$, $`\mathrm{sd}=916.70`$ | **0.51139** | **0.997837** |
| **DAN** | $`22671.11 \pm 28.37`$ ($`\mathrm{CV}=0.54837`$) | $`22604.01 \pm 32.72`$ ($`\mathrm{CV}=0.54935`$) | $`22642.35 \pm 21.44`$, $`\mathrm{sd}=12425.89`$ | **0.54879** | **0.999719** |
| **PAN** | $`3130.66 \pm 4.28`$ ($`\mathrm{CV}=0.59899`$) | $`3130.28 \pm 4.93`$ ($`\mathrm{CV}=0.59749`$) | $`3130.50 \pm 3.23`$, $`\mathrm{sd}=1873.12`$ | **0.59835** | **0.998201** |

Cross-scale correlation check command and output:
```text
IR : corr(C_1e5, C_1e6) = 0.997837, corr(C_1e2, C_1e6) = 0.877087 | 2400-tree SD of [1e5->1e6] slope: coupled=0.000298, indep=0.006400 (ratio=21.5x)
DAN: corr(C_1e5, C_1e6) = 0.999719, corr(C_1e2, C_1e6) = 0.947406 | 2400-tree SD of [1e5->1e6] slope: coupled=0.000115, indep=0.006875 (ratio=59.7x)
PAN: corr(C_1e5, C_1e6) = 0.998201, corr(C_1e2, C_1e6) = 0.863000 | 2400-tree SD of [1e5->1e6] slope: coupled=0.000318, indep=0.007496 (ratio=23.5x)
```

### 5.2 Pooled Local Slopes (192,000 Coupled Trees per Model, Exact Delta-Method SEs)

Real pasted output from `analysis_output.txt`:

```text
### Model: IR (coupled pooled M=192,000 trees, ref mu=0.565723)
  Local slopes (5pt_decade, coupled delta-method SE):
    [    100 ->    1000]: slope = 0.571265 +- 0.000228  (diff from ref: +0.005542)
    [   1000 ->   10000]: slope = 0.567160 +- 0.000122  (diff from ref: +0.001437)
    [  10000 ->  100000]: slope = 0.566176 +- 0.000064  (diff from ref: +0.000453)
    [ 100000 -> 1000000]: slope = 0.565766 +- 0.000033  (diff from ref: +0.000043)

  Local slopes (9pt_halfdec, coupled delta-method SE):
    [    100 ->     316]: slope = 0.572781 +- 0.000395  (diff from ref: +0.007058)
    [    316 ->    1000]: slope = 0.569751 +- 0.000288  (diff from ref: +0.004028)
    [   1000 ->    3162]: slope = 0.567495 +- 0.000209  (diff from ref: +0.001772)
    [   3162 ->   10000]: slope = 0.566825 +- 0.000152  (diff from ref: +0.001102)
    [  10000 ->   31623]: slope = 0.566192 +- 0.000110  (diff from ref: +0.000469)
    [  31623 ->  100000]: slope = 0.566159 +- 0.000079  (diff from ref: +0.000436)
    [ 100000 ->  316228]: slope = 0.565728 +- 0.000057  (diff from ref: +0.000005)
    [ 316228 -> 1000000]: slope = 0.565804 +- 0.000041  (diff from ref: +0.000081)

  Local slopes (17pt_qtrdec, coupled delta-method SE):
    [    100 ->     178]: slope = 0.573045 +- 0.000626  (diff from ref: +0.007322)
    [    178 ->     316]: slope = 0.572516 +- 0.000537  (diff from ref: +0.006793)
    [    316 ->     562]: slope = 0.570330 +- 0.000457  (diff from ref: +0.004607)
    [    562 ->    1000]: slope = 0.569173 +- 0.000388  (diff from ref: +0.003450)
    [   1000 ->    1778]: slope = 0.567788 +- 0.000330  (diff from ref: +0.002065)
    [   1778 ->    3162]: slope = 0.567201 +- 0.000282  (diff from ref: +0.001478)
    [   3162 ->    5623]: slope = 0.567305 +- 0.000240  (diff from ref: +0.001582)
    [   5623 ->   10000]: slope = 0.566346 +- 0.000204  (diff from ref: +0.000623)
    [  10000 ->   17783]: slope = 0.566240 +- 0.000173  (diff from ref: +0.000517)
    [  17783 ->   31623]: slope = 0.566145 +- 0.000148  (diff from ref: +0.000422)
    [  31623 ->   56234]: slope = 0.566162 +- 0.000126  (diff from ref: +0.000439)
    [  56234 ->  100000]: slope = 0.566156 +- 0.000106  (diff from ref: +0.000433)
    [ 100000 ->  177828]: slope = 0.565635 +- 0.000091  (diff from ref: -0.000088)
    [ 177828 ->  316228]: slope = 0.565822 +- 0.000077  (diff from ref: +0.000099)
    [ 316228 ->  562341]: slope = 0.565872 +- 0.000066  (diff from ref: +0.000149)
    [ 562341 -> 1000000]: slope = 0.565735 +- 0.000056  (diff from ref: +0.000012)

### Model: DAN (coupled pooled M=192,000 trees, ref mu=0.771192)
  Local slopes (5pt_decade, coupled delta-method SE):
    [    100 ->    1000]: slope = 0.772281 +- 0.000169  (diff from ref: +0.001089)
    [   1000 ->   10000]: slope = 0.771226 +- 0.000073  (diff from ref: +0.000033)
    [  10000 ->  100000]: slope = 0.771176 +- 0.000031  (diff from ref: -0.000016)
    [ 100000 -> 1000000]: slope = 0.771160 +- 0.000013  (diff from ref: -0.000033)

  Local slopes (9pt_halfdec, coupled delta-method SE):
    [    100 ->     316]: slope = 0.772663 +- 0.000299  (diff from ref: +0.001471)
    [    316 ->    1000]: slope = 0.771900 +- 0.000196  (diff from ref: +0.000708)
    [   1000 ->    3162]: slope = 0.771306 +- 0.000129  (diff from ref: +0.000113)
    [   3162 ->   10000]: slope = 0.771145 +- 0.000084  (diff from ref: -0.000047)
    [  10000 ->   31623]: slope = 0.771204 +- 0.000054  (diff from ref: +0.000012)
    [  31623 ->  100000]: slope = 0.771148 +- 0.000035  (diff from ref: -0.000044)
    [ 100000 ->  316228]: slope = 0.771138 +- 0.000023  (diff from ref: -0.000054)
    [ 316228 -> 1000000]: slope = 0.771181 +- 0.000015  (diff from ref: -0.000011)

  Local slopes (17pt_qtrdec, coupled delta-method SE):
    [    100 ->     178]: slope = 0.772730 +- 0.000482  (diff from ref: +0.001537)
    [    178 ->     316]: slope = 0.772596 +- 0.000392  (diff from ref: +0.001403)
    [    316 ->     562]: slope = 0.772180 +- 0.000317  (diff from ref: +0.000987)
    [    562 ->    1000]: slope = 0.771621 +- 0.000256  (diff from ref: +0.000429)
    [   1000 ->    1778]: slope = 0.771410 +- 0.000207  (diff from ref: +0.000218)
    [   1778 ->    3162]: slope = 0.771201 +- 0.000167  (diff from ref: +0.000009)
    [   3162 ->    5623]: slope = 0.771129 +- 0.000135  (diff from ref: -0.000063)
    [   5623 ->   10000]: slope = 0.771162 +- 0.000108  (diff from ref: -0.000031)
    [  10000 ->   17783]: slope = 0.771183 +- 0.000087  (diff from ref: -0.000009)
    [  17783 ->   31623]: slope = 0.771226 +- 0.000070  (diff from ref: +0.000033)
    [  31623 ->   56234]: slope = 0.771137 +- 0.000057  (diff from ref: -0.000055)
    [  56234 ->  100000]: slope = 0.771160 +- 0.000045  (diff from ref: -0.000032)
    [ 100000 ->  177828]: slope = 0.771102 +- 0.000036  (diff from ref: -0.000090)
    [ 177828 ->  316228]: slope = 0.771174 +- 0.000029  (diff from ref: -0.000018)
    [ 316228 ->  562341]: slope = 0.771185 +- 0.000023  (diff from ref: -0.000007)
    [ 562341 -> 1000000]: slope = 0.771177 +- 0.000019  (diff from ref: -0.000016)

### Model: PAN (coupled pooled M=192,000 trees, ref mu=0.549735)
  Local slopes (5pt_decade, coupled delta-method SE):
    [    100 ->    1000]: slope = 0.605601 +- 0.000243  (diff from ref: +0.055866)
    [   1000 ->   10000]: slope = 0.568568 +- 0.000128  (diff from ref: +0.018833)
    [  10000 ->  100000]: slope = 0.556273 +- 0.000067  (diff from ref: +0.006538)
    [ 100000 -> 1000000]: slope = 0.551934 +- 0.000036  (diff from ref: +0.002199)

  Local slopes (9pt_halfdec, coupled delta-method SE):
    [    100 ->     316]: slope = 0.619985 +- 0.000375  (diff from ref: +0.070250)
    [    316 ->    1000]: slope = 0.591235 +- 0.000269  (diff from ref: +0.041501)
    [   1000 ->    3162]: slope = 0.573403 +- 0.000194  (diff from ref: +0.023668)
    [   3162 ->   10000]: slope = 0.563734 +- 0.000140  (diff from ref: +0.013999)
    [  10000 ->   31623]: slope = 0.557996 +- 0.000102  (diff from ref: +0.008261)
    [  31623 ->  100000]: slope = 0.554549 +- 0.000074  (diff from ref: +0.004815)
    [ 100000 ->  316228]: slope = 0.552478 +- 0.000054  (diff from ref: +0.002743)
    [ 316228 -> 1000000]: slope = 0.551389 +- 0.000039  (diff from ref: +0.001654)

  Local slopes (17pt_qtrdec, coupled delta-method SE):
    [    100 ->     178]: slope = 0.629480 +- 0.000549  (diff from ref: +0.079745)
    [    178 ->     316]: slope = 0.610446 +- 0.000462  (diff from ref: +0.060711)
    [    316 ->     562]: slope = 0.596690 +- 0.000390  (diff from ref: +0.046955)
    [    562 ->    1000]: slope = 0.585786 +- 0.000329  (diff from ref: +0.036051)
    [   1000 ->    1778]: slope = 0.576481 +- 0.000280  (diff from ref: +0.026746)
    [   1778 ->    3162]: slope = 0.570327 +- 0.000237  (diff from ref: +0.020592)
    [   3162 ->    5623]: slope = 0.565678 +- 0.000201  (diff from ref: +0.015943)
    [   5623 ->   10000]: slope = 0.561790 +- 0.000172  (diff from ref: +0.012055)
    [  10000 ->   17783]: slope = 0.559227 +- 0.000146  (diff from ref: +0.009492)
    [  17783 ->   31623]: slope = 0.556764 +- 0.000124  (diff from ref: +0.007029)
    [  31623 ->   56234]: slope = 0.555253 +- 0.000106  (diff from ref: +0.005518)
    [  56234 ->  100000]: slope = 0.553846 +- 0.000091  (diff from ref: +0.004111)
    [ 100000 ->  177828]: slope = 0.552833 +- 0.000077  (diff from ref: +0.003098)
    [ 177828 ->  316228]: slope = 0.552124 +- 0.000066  (diff from ref: +0.002389)
    [ 316228 ->  562341]: slope = 0.551664 +- 0.000056  (diff from ref: +0.001929)
    [ 562341 -> 1000000]: slope = 0.551113 +- 0.000048  (diff from ref: +0.001378)
```

Notice the clean structural contrast in finite-size drift:
- **DAN**: denominator is $`N`$, no $`O(N^{\mu-1})`$ correction; local slope is already within $`3\times 10^{-5}`$ of $`\mu_{\text{DAN}} = 0.771192`$ by $`N = 10^3`$.
- **IR**: modest positive finite-size drift (`+0.0055` on $`[10^2,10^3]`$ $`\to`$ `+0.0014` $`\to`$ `+0.00045` $`\to`$ `+0.000043 ± 0.000033` on $`[10^5,10^6]`$).
- **PAN**: large positive $`O(N^{\mu-1}) = O(N^{-0.450})`$ drift from the $`(N-\mathcal{N})^{-1}`$ rejection denominator, decaying by $`10^{-0.450} \approx 0.355`$ per decade (`+0.0559` $`\to`$ `+0.0188` $`\to`$ `+0.0065` $`\to`$ `+0.0022`).

---

## 6. Final Conclusion

The data unambiguously support **(a) pure sampling noise of the protocol** (independent 2400-tree ensembles at each grid size $`N \in \{10^2, 10^3, 10^4, 10^5, 10^6\}`$ fitted by unweighted nonlinear least squares in linear space, $`\langle C \rangle = a N^b`$). Because $`C(N)`$ is non-self-averaging ($`\mathrm{CV} \approx 0.51\text{–}0.60`$), each independent size point has $`\sim 1\%`$ sampling noise, and unweighted linear-space NLS concentrates its leverage on $`N = 10^5`$ and $`N = 10^6`$, yielding replicate standard deviations $`\mathrm{SD}(b) \approx 0.004\text{–}0.007`$. Thus $`\mu \approx 0.569`$ for IR ($`+0.59\sigma`$, 27% of replicates), $`\mu \approx 0.773`$ for DAN ($`+0.32\sigma`$, 35% of replicates), and $`\mu \approx 0.545`$ for PAN ($`-1.29\sigma`$, 7%–13% of replicates; 23%–28% with additive constant) are completely typical realizations of Hartle et al.'s fitting protocol and provide **no evidence** for $`\mu_{\text{IR}}`$ departing from $`0.565723 \pm 0.000003`$.
