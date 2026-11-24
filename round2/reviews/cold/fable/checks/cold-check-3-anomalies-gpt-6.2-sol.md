ANOMALY-A: CHANCE-PLAUSIBLE
ANOMALY-B: UNEXPLAINED

A reproduces in the historical sample. Its 3.35-SE contrast was selected among many correlated looks, and the paths added after the earlier 3.80-SE look do not reinforce it. The discrepancy is evidence against using (C) for both full and ANPQ at once. It cannot locate an absolute bias in either one. B reproduces under exact checkpoint pairing and the declared run blocks, including on phase 2 alone. I found no chunk-construction, seed-collision, or reset mechanism. The few disjoint post-60 blocks also have low variance, with little discriminatory power.

This audit reads the production files and writes only under `checks/cold-estimator`. It does not read the late-cohort output. Every compute command was pinned to one core, with a timeout of at most 240 seconds on the long commands. The only builds are scratch test programs under `task3`.

## A. Independent paired arithmetic

I read only complete phase-2 files with the expected per-file row count, exactly 14 columns, the stated layout, matching checkpoint order and each expected run ID. The 120 chronological main chunks end at 06:08:05.897595 UTC, with counts 41, 41 and 38 for workers 0, 1 and 2. The rungs use the stated 52 and 51 complete phase-2 chunks respectively. I formed the two ratios on the same paths and computed the exact leave-one-run jackknife of their *difference*. The centered influence estimates used by the lead agree up to negligible finite-sample differences.

| Full minus ANPQ | Runs | Difference | Paired SE | z |
| --- | ---: | ---: | ---: | ---: |
| Main 20→24 | 12,000 | +6.5663e-8 | 1.9594e-8 | +3.351 |
| Main 16→20 | 12,000 | −1.1060e-7 | 6.1203e-8 | −1.807 |
| Independent rung 16→20 | 25,500 | +9.2459e-8 | 4.2247e-8 | +2.189 |
| Independent rung 12→16 | 260,000 | −7.3080e-8 | 4.0748e-8 | −1.793 |
| Main 20→24, first observed 10,200 | 10,200 | +8.0620e-8 | 2.1189e-8 | +3.805 |
| Main 20→24, next 1,800 only | 1,800 | −1.7452e-8 | 5.1379e-8 | −0.340 |

Both 6,000-path chronological halves and all three worker aggregates are positive. Those are internal consistency checks; they include the signal that led to the selected contrast. The next 1,800 paths have weak power but provide an actual disjoint continuation. Paired residual kurtosis was 4.26, the largest absolute residual was 6.16 SDs, and dropping the 100 largest absolute residuals gives `+7.2946e-8 ± 1.8833e-8`. Outliers do not explain the sign.

The different 16→20 estimates refer to the *same population window and observables*. Their independent difference (rung minus main) is `+2.0306e-7 ± 7.4369e-8` (2.73 SE). A time-varying population-mode explanation alone cannot create that between-cohort discrepancy. The main simulation has an extra t=12 checkpoint, but in the ideal lazy-Yule simulator exponential memorylessness and the Yule semigroup preserve the same marginal joint law at 16 and 20. Added checkpoints can change the realized RNG path under the same seed, which is different from changing the law.

An illustrative Gaussian max-|z| reference preserving empirical within-run covariance gives 1.58% for at least the observed late ANPQ statistic across 20 named contrasts. Including an earlier main look increases this to 2.10% across 32 contrasts. The pointwise two-sided normal tail is 0.0805%. These are references under an all-zero contrast null, not a reconstruction of the actual exploratory stopping and reporting rule. They put a chance excursion in context without making a multiplicity-adjusted claim about the experiment.

The independent five-weight early rungs also reproduce the large early transients. Relative to the working reference 0.56572366, full and ANPQ are `+1.0764e-4` and `+1.4382e-4` at 4→8, then `+2.0205e-6` and `+2.1218e-6` at 8→12. The latter difference is `−1.013e-7 ± 1.175e-7`. Early shrink does not distinguish a small later component from noise.

```text
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/raw.py
main 120 phase2 complete chunks, workers: {0: (41, 500, 540), 1: (41, 500, 540), 2: (38, 500, 537)} last snapshot mtime 2026-11-24 06:08:05.897595+00:00
main(12, 16): n=12000 full=0.5657237489087 ±8.6759e-08; ANPQ=0.5657237177354 ±2.0858e-07; full-ANPQ=+3.1173337e-08 ±1.878559e-07 z=+0.1659
main(16, 20): n=12000 full=0.5657237129174 ±2.7865e-08; ANPQ=0.5657238235174 ±6.7269e-08; full-ANPQ=-1.1059998e-07 ±6.120342e-08 z=-1.8071
main(20, 24): n=12000 full=0.5657236417946 ±8.8292e-09; sane=0.5657236531615 ±1.0726e-08; full-sane=-1.1366848e-08 ±6.00051e-09 z=-1.8943
main(20, 24): n=12000 full=0.5657236417946 ±8.8292e-09; AN=0.5657234991645 ±8.5025e-08; full-AN=+1.4263018e-07 ±8.450769e-08 z=+1.6878
main(20, 24): n=12000 full=0.5657236417946 ±8.8292e-09; ANPQ=0.5657235761316 ±2.1618e-08; full-ANPQ=+6.5663066e-08 ±1.95939e-08 z=+3.3512
main(20, 24): n=12000 full=0.5657236417946 ±8.8292e-09; AR=0.5657235852321 ±6.3257e-08; full-AR=+5.6562517e-08 ±6.253833e-08 z=+0.9044
r1216(12, 16): n=260000 full=0.5657237030344 ±1.8364e-08; ANPQ=0.5657237761147 ±4.4787e-08; full-ANPQ=-7.3080352e-08 ±4.074754e-08 z=-1.7935
r1620(16, 20): n=25500 full=0.5657236423605 ±1.8814e-08; ANPQ=0.5657235499011 ±4.6235e-08; full-ANPQ=+9.2459419e-08 ±4.224745e-08 z=+2.1885
chronological first half: n=6000 full=0.5657236521499 ±1.2425e-08; ANPQ=0.5657235801739 ±3.0217e-08; full-ANPQ=+7.197592e-08 ±2.73927e-08 z=+2.6276
chronological second half: n=6000 full=0.5657236313962 ±1.2549e-08; ANPQ=0.5657235720723 ±3.0933e-08; full-ANPQ=+5.9323882e-08 ±2.802846e-08 z=+2.1166
chronological initial 10,200: n=10200 full=0.5657236454597 ±9.5922e-09; ANPQ=0.5657235648399 ±2.3424e-08; full-ANPQ=+8.0619785e-08 ±2.118858e-08 z=+3.8049
new 10,200 to 12,000: n=1800 full=0.5657236214275 ±2.2598e-08; ANPQ=0.5657236388800 ±5.6118e-08; full-ANPQ=-1.7452464e-08 ±5.137915e-08 z=-0.3397
worker 0: n=4100 full=0.5657236372238 ±1.5122e-08; ANPQ=0.5657235695568 ±3.6951e-08; full-ANPQ=+6.7666986e-08 ±3.371736e-08 z=+2.0069
worker 1: n=4100 full=0.5657236640622 ±1.4993e-08; ANPQ=0.5657236007625 ±3.6452e-08; full-ANPQ=+6.3299711e-08 ±3.295841e-08 z=+1.9206
worker 2: n=3800 full=0.5657236225278 ±1.5793e-08; ANPQ=0.5657235564922 ±3.9056e-08; full-ANPQ=+6.6035637e-08 ±3.523879e-08 z=+1.8739
main comparison: n=12000 full=0.5657237129174 ±2.7865e-08; ANPQ=0.5657238235174 ±6.7269e-08; full-ANPQ=-1.1059998e-07 ±6.120342e-08 z=-1.8071
r1620 comparison: n=25500 full=0.5657236423605 ±1.8814e-08; ANPQ=0.5657235499011 ±4.6235e-08; full-ANPQ=+9.2459419e-08 ±4.224745e-08 z=+2.1885
independent at same window (r1620 minus main): +2.030594e-07 ±7.436872e-08 z=+2.7304
ANPQ paired residual kurtosis=4.262 max/core paired residual SD=6.157
drop largest 10: n=11990 full=0.5657236419885 ±8.8358e-09; ANPQ=0.5657235776226 ±2.1481e-08; full-ANPQ=+6.43659e-08 ±1.942706e-08 z=+3.3132
drop largest 100: n=11900 full=0.5657236402330 ±8.903e-09; ANPQ=0.5657235672865 ±2.0964e-08; full-ANPQ=+7.2946421e-08 ±1.883259e-08 z=+3.8734
Gaussian max-|z| on 20 specified contrasts (all-zero null), individual two-sided p=0.000804622, family chance >= |ANPQ late z| 3.3512: 0.015812; observed max|z|=3.3614
Gaussian max-|z| on 32 specified contrasts (all-zero null), individual two-sided p=0.000804622, family chance >= |ANPQ late z| 3.3512: 0.020992; observed max|z|=3.8049
```

### Consequence for contraction across observables

Let Lf and La be the full and ANPQ population ratios at 20→24, and Ef and Ea the respective ratios at 12→16. If both observables obey both intervening contractions with factor rho, put `c=rho^2`. From `|Lf−mu|<=c|Ef−mu|` we get `|Lf−mu|<=c/(1−c)|Ef−Lf|`. Applying the corresponding ANPQ inequality and substituting gives

```text
|Lf−La| <= c*(1+c)/(1−c)*|Ef−Lf| + c*|Ea−Lf|.
```

This removes the unknown true mu. Using the largest rho allowed by the proved exponent bracket for this check, and conservative approximate simultaneous bounds for these three contrasts, I obtain a data lower bound for Lf−La of `2.397e-8` and a joint-(C) upper bound of `3.877e-9`. This is a selected diagnostic rather than an adjusted confirmatory test. It gives a specific reason to restrict any conditional assumption to the primary observable. The data do not decompose the discrepancy into two absolute biases.

```text
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/secondary.py
dual-C check: rho_upper=0.104350485 rho_squared=0.0108890237
approximate simultaneous bounds for late full−ANPQ: data lower=2.39663605e-08, under joint (C) upper=3.87698647e-09 (Bonferroni over these three contrasts only)
  early full−late full upper abs=1.10020114e-07; early ANPQ−late full upper abs=2.43602856e-07; z values 2.12805,2.39398
```

### P/Q generator, output and weights

I traced all P/Q branches in the shipped source. Own-state gates use L<=2 with c<=3. Pair degrees are integer-capped before indexing and sorted once. For paths, the middle degree must be <7 and the far degree is capped before indexing. A gain at the middle can change its type or remove the term at the threshold. A conversion at the middle can change the positive-leaf flag and adds the new far degree-2 vertex. A far conversion leaves its degree unchanged. A far gain only changes the entry if it crosses a degree bucket. Creation of a center or a neighbor follows the same cases. The newborn's new edge alone adds no path through itself, while the newborn as center sees paths through its parent if the middle passes the cutoff. The exact test below exercises those cases without using the production's incremental P/Q emitters to construct the reference drifts.

The C++ source keeps L and D in uint64 before the caps. In the audited late sample the largest core count is 5,198,452, well below the uint32 range of the core indexing. Large L at a remote hub is either clipped in the feature or absent from P/Q. If a gain changes a P bucket from the zero-leaf class at a large-degree neighbor, its gain rate is the neighbor credit, rather than its large leaf count. A path with a large far degree has a stable bucket. Nothing in these paths switches by time, observable name or core count. The output writes all five dot products in a fixed order on each trajectory. The independent shape/ID check found matching column counts and trajectories, can reject truncated chunks, and agrees with the previous frozen-binary row replays.

For a targeted test I built a scratch copy of the tree evaluator from the current source and evaluated 381 synthetic rooted core trees, using exact Fraction rates from the oriented-pair multiplicities and whole-tree P/Q recomputation for every possible attachment class. Maximum hub L was 2^53. The tests included both sides of 32-bit limits, zero-leaf hubs, and capped path and pair classes. Absolute P and Q drift discrepancies were at most 6.13e-14 and 5.24e-14. This is especially useful alongside the lead's full-size replay because the reference event enumeration and feature definitions are separate. The tests and source inspection find no candidate P/Q defect at the scale of A.

The 100,000 training rows show no exceptional ANPQ coefficient on a rare column. With both plain weights normalized to 1, the largest other absolute weights are 0.5413 (ANPQ, column occurs in 35.95% of rows) and 1.0723 (full, R, 100% of rows). Among columns activated in less than 2.5% of training rows, maxima are 0.1107 (ANPQ) and 0.1382 (full). Here activation means a nonzero differenced generator or value component in a training row, not late-vertex prevalence. The largest ANPQ P coefficient is 0.1180 on a column activated in 10.99% of runs, and its largest Q coefficient is 0.1018 on a column activated in 80.40%. The hub own and constant large-neighbor coefficients are comparable to full. The saturated degrees include moderate degrees and cannot single out the largest hubs. Late denominator scales are also healthy. The mean 20→24 denominator increments are about 793,421 for ANPQ and 789,680 for full. At t=24, den/C is 0.86947 and 0.86537 respectively.

ANPQ is a refit excluding R, rather than full with only P/Q terms perturbed. The maximum coefficient change in A is 0.385, compared with 0.085 in P and 0.053 in Q. A real differential finite-time component can project on any of these changes. The gap between two fitted ratios does not identify a P/Q defect or assign the bias to ANPQ.

```text
$ taskset -c 15 env TMPDIR="$PWD/task3/tmp" timeout 240 g++ -O2 -std=c++17 ../../src/feat_tree.cpp -o task3/feat_tree
$ taskset -c 15 timeout 240 python3 -B task3/pq_hubs.py
381 aggregate explicit-event cases, max tested hub L=9007199254740992
P exercised slots 272 max abs drift error 6.128709089882341e-14 max scaled drift error 3.991758773716465e-15
Q exercised slots 463 max abs drift error 5.238095238095238e-14 max scaled drift error 4.322455191571782e-15
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/weights.py
../../data/theta_full.f64 7d825f142759c83b8890cb61585d5d59489bffa26c33e211188f192169eb688f
../../data/theta_ANPQ.f64 a7cfe2c93a9c07d3cfe2426e512a5eb24ec69680cbc7c6351b730ba2936961c2
transition/value nonzero fractions over 100000 training paths
full
  A     nonzero   83 max|w|=0.4611644 freq of max=1; min freq among weighted=0.35952; index 2 A
  H     nonzero   21 max|w|=0.08311445 freq of max=0.99995; min freq among weighted=0.99995; index 99 H
  N     nonzero  952 max|w|=0.1856315 freq of max=0.24791; min freq among weighted=0.01124; index 224 N
  B     nonzero  186 max|w|=0.151492 freq of max=1; min freq among weighted=0.99508; index 1755 B
  N2    nonzero  222 max|w|=0.04025822 freq of max=0.9958; min freq among weighted=0.01381; index 2153 N2
  B2    nonzero   42 max|w|=0.03782051 freq of max=0.99998; min freq among weighted=0.99985; index 2391 B2
  P     nonzero  383 max|w|=0.09213997 freq of max=0.99999; min freq among weighted=0.01036; index 2573 P(L=0,c=2,ct1=11,ct2=11)
  Q     nonzero  504 max|w|=0.0501656 freq of max=0.80401; min freq among weighted=0.01964; index 3480 Q(L=0,c=3,ct=2,Dfar=6)
  R     nonzero  813 max|w|=1.072279 freq of max=1; min freq among weighted=0.03383; index 4210 R
plain 1.0
  top: 4210 R w=-1.072279 freq=1
  top: 4211 R w=-0.8289877 freq=1
  top: 4216 R w=-0.7084665 freq=1
  top: 4360 R w=-0.5835454 freq=0.99999
  top: 4091 R w=-0.5214825 freq=1
  top: 4107 R w=-0.4941479 freq=0.99979
  top: 2 A w=-0.4611644 freq=1
  top: 4105 R w=-0.4424498 freq=0.99986
  all features with frequency<2.5%: 76 max|w|=0.138245
  hub-related count block 96:120:3 nonzero=7 max|w|=0.0831145
  hub-related count block 1740:2010:3 nonzero=62 max|w|=0.151492
  hub-related count block 2370:2430:3 nonzero=14 max|w|=0.0378205
  P exposure to exactly capped remote/high-neighbor degree: 90 nonzero max|w|=0.09214
  Q exposure to exactly capped remote/high-neighbor degree: 72 nonzero max|w|=0.0265078
ANPQ
  A     nonzero   83 max|w|=0.5412872 freq of max=0.35952; min freq among weighted=0.35952; index 7 A
  H     nonzero   21 max|w|=0.08524824 freq of max=0.99995; min freq among weighted=0.99995; index 99 H
  N     nonzero  952 max|w|=0.1515105 freq of max=0.24791; min freq among weighted=0.01124; index 224 N
  B     nonzero  186 max|w|=0.1374705 freq of max=1; min freq among weighted=0.99508; index 1755 B
  N2    nonzero  222 max|w|=0.08681506 freq of max=0.12091; min freq among weighted=0.01381; index 2113 N2
  B2    nonzero   42 max|w|=0.03769383 freq of max=1; min freq among weighted=0.99985; index 2394 B2
  P     nonzero  383 max|w|=0.1180176 freq of max=0.1099; min freq among weighted=0.01036; index 2681 P(L=0,c=3,ct1=8,ct2=11)
  Q     nonzero  504 max|w|=0.1018461 freq of max=0.80401; min freq among weighted=0.01964; index 3480 Q(L=0,c=3,ct=2,Dfar=6)
R zero
plain 1.0
  top: 7 A w=-0.5412872 freq=0.35952
  top: 2 A w=-0.5033276 freq=1
  top: 4 A w=-0.4313854 freq=0.99964
  top: 6 A w=-0.3503051 freq=0.8341
  top: 5 A w=-0.3162085 freq=0.98993
  top: 15 A w=-0.296058 freq=0.47164
  top: 3 A w=-0.2910925 freq=0.99999
  top: 21 A w=-0.2121997 freq=0.99364
  all features with frequency<2.5%: 76 max|w|=0.110725
  hub-related count block 96:120:3 nonzero=7 max|w|=0.0852482
  hub-related count block 1740:2010:3 nonzero=62 max|w|=0.13747
  hub-related count block 2370:2430:3 nonzero=14 max|w|=0.0376938
  P exposure to exactly capped remote/high-neighbor degree: 90 nonzero max|w|=0.118018
  Q exposure to exactly capped remote/high-neighbor degree: 72 nonzero max|w|=0.0256105
diff A: RMS=0.0584625 max|diff|=0.384773
diff H: RMS=0.00104595 max|diff|=0.00221448
diff N: RMS=0.0169955 max|diff|=0.183765
diff B: RMS=0.0089639 max|diff|=0.0757577
diff N2: RMS=0.00908828 max|diff|=0.0883571
diff B2: RMS=0.00590273 max|diff|=0.0156856
diff P: RMS=0.0140001 max|diff|=0.0851938
diff Q: RMS=0.00661204 max|diff|=0.0530405
diff R: RMS=0.0855704 max|diff|=1.07228
diff plain: RMS=0 max|diff|=0
```

### Cheapest discriminating follow-up

Use the already running, independent stream 9500. I would register the full-minus-ANPQ contrast at **20→24** there before opening its results, alongside the existing **24→28** rule. The same-window check can replicate A without confounding replication with a changed time window. A nonzero contrast can decay before 24→28 even when the original 20→24 contrast was real. A significant same-window repeat plus a paired change at 24→28 would support differential finite-time effects. Componentwise residual diagnostics would be useful only if a specific P/Q code mechanism remains suspect. Neither contrast identifies which observable is closer to mu. Also, `|L3|+2 SE` estimates an upper allowance for a difference of finite-window targets. An absolute bias bound would additionally need a premise about the accuracy of the late target.

### IR-0 complex root

I implemented a separate floating-point degree recursion in C++ using the higher-order boundary coefficients from the independent real calculation. The upper c boundary is zero here, and I varied the cutoff and tail order. At degree 4096 with tail order 6 I obtain `-0.5513165560428824 ± 0.5015275330498485 i`, with a residual of `9.49e-17` in that finite recursion and stable leading digits across orders and horizons. The leading real root remains a useful positive-real check.

The direct Laplace transform is guaranteed by the bounded birth intensity only in Re(s)>0. Negative real parts here use continuation of the resolvent, outside the certified real-s calculation. On the rectangle in the brief I also computed numerical contour winding. The 800-step grid misses the winding, while 1,600 and 3,200 steps give winding 2 for the rational truncated recursions. Their one pole at zero would make that consistent with three zeros. Numerical winding on finite recursions cannot certify an exclusion for the infinite continuation. It does show sensitivity to contour resolution.

The complex pair's real-part gap gives an envelope reduction of about 62.45 over four time units in IR-0. An oscillatory mode can still have a very small projection on one window and a larger absolute projection on the next. Neither the pair nor its envelope proves the per-window magnitude contraction in (C), even in IR-0. A fast negative and a fast positive mode can also cancel on one window, so a late sign change need not imply a slow mode. These facts do not transfer a spectral gap to the friend tree.

```text
$ taskset -c 15 env TMPDIR="$PWD/task3/tmp" timeout 240 g++ -O2 -std=c++17 task3/ir0_complex.cpp -o task3/ir0_complex
$ taskset -c 15 timeout 240 task3/ir0_complex
D=256 ccap=28 order=1 root=(-0.5510605042325062,0.5029999174290268) residual=1.753345416651335e-18 real reference residual=5.617754509607432e-06
D=256 ccap=28 order=3 root=(-0.5513156481464873,0.5015283115142305) residual=7.047288728641609e-18 real reference residual=1.589290595719987e-09
D=1024 ccap=44 order=4 root=(-0.5513165558626979,0.5015275330296293) residual=2.244581760823769e-17 real reference residual=4.633923453289945e-14
D=4096 ccap=72 order=6 root=(-0.5513165560428824,0.5015275330498485) residual=9.485905681187956e-17 real reference residual=1.084202172485504e-19
finite-recursion winding D=256 ccap=32 order=4 vertical steps=800 winding=7.698885176413028e-22 min contour |F|=0.206680605412113 max sampled arg step=2.757501842744406
finite-recursion winding D=256 ccap=32 order=4 vertical steps=1600 winding=2 min contour |F|=0.206680605412113 max sampled arg step=2.358102348664464
finite-recursion winding D=512 ccap=40 order=5 vertical steps=1600 winding=2 min contour |F|=0.2066806054255635 max sampled arg step=2.358042589509524
finite-recursion winding D=512 ccap=40 order=5 vertical steps=3200 winding=2 min contour |F|=0.2066806054255635 max sampled arg step=1.312999959494075
```

## B. Block underdispersion

I used the declared 5,000-path run-index blocks. For a given ratio I computed the path residual `e=Delta num−Rhat*Delta den` once, then applied both the actual block partition and 1,999 uniform run permutations to that same centered vector. The statistic T is `sum of squared block sums / sum of squared individual residuals`. A common permutation across columns preserves correlations. The finite-sample null mean is close to 1−1/m for m blocks, which is also borne out by the permutations.

| Rung 12→16 | Blocks | T plain | T sane | T full | Lower p sane / full |
| --- | ---: | ---: | ---: | ---: | --- |
| Chronological first 60 | 60 | 1.046 | 0.534 | 0.462 | 0.0035 / 0.0005 |
| Phase 2 only | 52 | 1.269 | 0.540 | 0.467 | 0.0035 / 0.0015 |
| Seven disjoint blocks after the first 60 | 7 | 1.843 | 0.439 | 0.461 | 0.2275 / 0.2350 |

Sane and full path residuals have correlation 0.828. The permutation lower-tail probability for the minimum T of the two scores is 0.0015 for either large comparison. This accounts for this two-score look only. The 15 old phase-1 blocks give T=0.543 and 0.422; phase 2 alone gives much the same effect, so mixing two- and five-output files does not produce it. AN, ANPQ and AR in phase 2 have T=1.093, 1.379 and 0.848. The corresponding 16→20 production-process blocks give T near one for both sane and full. On the new early five-weight 8→12 cohort their T values are 0.960 and 0.955. These checks help localize the observation to a shared fitted residual on this rung. They do not identify its cause.

Staggering the block origin raises phase-2 full T from 0.467 to 0.757 at offset 2500, but it is still 0.506 at offset 3750. Short-lag sample correlations were of order a few 1e-3 or less. If the effect were a uniform within-block anticorrelation, approximately −1.1e-4 per pair would suffice. These short-lag summaries cannot exclude that magnitude or a nonuniform long-range pattern.

The continuity of the output run IDs, the persistent effect in phase 2, and per-path re-seeding exclude a simple overlap or checkpoint mismatch within a chunk. The source calls `seed2` on every run, resets the tree, builds a fresh snapshot credit vector, and resets weighted accumulators per checkpoint. A file boundary affects the executable's process bookkeeping; the checked replay evidence is consistent with no process-boundary dependence for a fixed run ID. Each input path still shares an initialization family, so absence of literal collisions does not prove independence. The anomaly is post hoc, with many windows, controls and possible block sizes considered, and the seven unused blocks lack power to resolve it. I would retain the larger run-level SE while treating the pattern as unresolved.

The clean confirmatory test for B is a fresh, prespecified independent-seed cohort at 12→16, inspecting just the 5,000-path sane/full common score and its aligned block statistic. If attribution to seeding matters, compare a separate seed family or an independently salted seeding method. A few blocks can be timed cheaply, but approximately 30–60 would make a replication of T near 0.5 informative. Replaying existing seeds or shuffling existing values can diagnose bookkeeping or quantify the observation, respectively. They cannot provide a new replication. I ran no extra simulated paths for this audit.

The supplied new seed provenance also passes the exact initialized-state check. Across 29 selected input and validation ranges, I computed 7,121,500 pre-seeds with no collisions. For the 335,000 investigated stream-9200 paths I checked for partial overlaps of the four SplitMix seed words, which would require a shift of one, two or three multiples of the SplitMix increment. There were no matches. No starting state matched a state advanced by up to 16 xoshiro steps within these paths. The first few uniform-output block tests give T between 0.77 and 1.21. The tests address short-range overlap and basic seed initialization, with no conclusion about correlations deep into the trajectories.

The standardized sum of sane and full residuals carries the low block variance. In phase 2 it has T=0.432 and an exploratory permutation lower tail of 0.0017. The standardized difference has T=1.263 and lower tail 0.907. This is consistent with one common score direction driving the two reported low T values.

```text
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/rng_blocks.py
r1216 seed 9200, offset matches y_j = y_i + offset*gamma
  SplitMix state seed word window shift 1 matches 0
  SplitMix state seed word window shift 2 matches 0
  SplitMix state seed word window shift 3 matches 0
  xoshiro output 1: exact state matches=0, uniform block T=0.768729
  xoshiro output 2: exact state matches=0, uniform block T=1.12016
  xoshiro output 4: exact state matches=0, uniform block T=1.20903
  xoshiro output 8: exact state matches=0, uniform block T=0.8431
  xoshiro output 16: exact state matches=0, uniform block T=1.00223
phase2 r1216 z_sane + z_full, z_sane - z_full:
  T [0.43187986 1.262964  ] ordinary lower p [0.00166667 0.90666667]
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/provenance.py
r04_08.txt seed 9410, IDs validated from 0 to 1999999
r08_12.txt seed 9420, IDs validated from 0 to 299999
r04_08_five.txt seed 9440, IDs validated from 0 to 999999
  near-reference offsets: [('plain', 7.354646188328928e-05), ('sane', -7.459876174842339e-05), ('full', 0.00010763902303134731), ('AN', 0.00014728914838191276), ('ANPQ', 0.00014382273673618773), ('AR', -4.781254620589426e-05)]
r04_08_five.txt: n=1000000 full=0.5658312990230 ±7.8411e-07; ANPQ=0.5658674827367 ±8.2813e-07; full-ANPQ=-3.6183714e-05 ±2.404343e-07 z=-150.4932
r08_12_five.txt seed 9430, IDs validated from 0 to 299999
  near-reference offsets: [('plain', -4.221501942791583e-06), ('sane', -4.938549029720463e-07), ('full', 2.0204840837090288e-06), ('AN', 2.1264587849456618e-06), ('ANPQ', 2.1217780408200326e-06), ('AR', -9.356195596055272e-07)]
r08_12_five.txt: n=300000 full=0.5657256804841 ±6.3284e-08; ANPQ=0.5657257817780 ±1.339e-07; full-ANPQ=-1.0129396e-07 ±1.174642e-07 z=-0.8623
actual selected inputs and known controls: 29 ranges, 7121500 pre-seeds, collisions=0
sane bc71186625e0bb8ff19d97db83e088285e7fdc9f26118541dce2444b971cbb3c
full 7d825f142759c83b8890cb61585d5d59489bffa26c33e211188f192169eb688f
AN ccf909d9658ecbeb0cc584ffdfbcc896b22de59c495bfa4a341a86e20c20145e
ANPQ a7cfe2c93a9c07d3cfe2426e512a5eb24ec69680cbc7c6351b730ba2936961c2
AR f781b2af42706b24064a475dfe399e940772b7c5c7c82ee6af9c12662be0d179
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/secondary.py
phase2 r1216 (all five and plain): 52 x 5000 paths
  plain T=1.2694 lower p=0.9175
  sane T=0.54019 lower p=0.005
  full T=0.46721 lower p=0.0025
  AN T=1.0927 lower p=0.7175
  ANPQ T=1.3785 lower p=0.9725
  AR T=0.84848 lower p=0.295
phase2 r1620: 51 x 500 paths
  plain T=0.9642 lower p=0.477
  sane T=1.0024 lower p=0.585
  full T=0.95907 lower p=0.491
  AN T=0.9678 lower p=0.495
  ANPQ T=1.0828 lower p=0.712
  AR T=1.0853 lower p=0.706
old r1216: 15 x 5000 paths
  plain T=0.6302 lower p=0.18
  sane T=0.54283 lower p=0.129
  full T=0.42232 lower p=0.041
round2/checks/cold-estimator/task3/raw.py:29: UserWarning: loadtxt: input contained no data: "../../outputs/prod/main_w1_c00545.txt"
  x = np.loadtxt(path,ndmin=2)
main 12:16, 10 consecutive processes per worker block: 11 x 1000 paths
  plain T=1.2165 lower p=0.824
  sane T=0.41829 lower p=0.078
  full T=0.77196 lower p=0.419
  AN T=0.7342 lower p=0.367
  ANPQ T=0.23734 lower p=0.01
  AR T=0.54165 lower p=0.166
early 8:12 five: 60 x 5000 paths
  plain T=0.91944 lower p=0.3425
  sane T=0.95964 lower p=0.45
  full T=0.95469 lower p=0.455
  AN T=1.2745 lower p=0.9375
  ANPQ T=1.0196 lower p=0.66
  AR T=0.63595 lower p=0.0225
```

## New complex renewal scan and weaker premise (C′)

I independently decoded the completed T=17 binary, reproduced the reported eigenvalues, and varied the tail and high-degree bin closures. At the real candidate s=0.5657237 the supplied leading value is 0.997691 with batch jackknife SE 0.00506. At s=0.5657237+0.05i the nearest eigenvalue is 0.971307−0.157111i, with complex jackknife RMS 0.00507. Across the 7×40 sampled nodes, excluding the real leading candidate, the closest eigenvalue distance is 0.05123 and the smallest singular value of I−K is 0.01829. The no-tail, midpoint-bin and alternate-moment schemes also give no closer root on these nodes. In particular all sampled geometric tail ratios have magnitude below one.

### What the scan can support

The scan is a numerical motivation for trying a slower scale in a sensitivity analysis. The absolute bound under (C′) is valid **conditional on the stated finite-window contraction** and positive denominator increments, using the same exact bound algebra as in REPORT2. The factor `exp(-4*(0.565-0.30)) = 0.34645581033` uses the established lower bound for mu so that it is conservative relative to a hypothesized line at 0.30. At the earlier approximate w=0.0942451, the combined-window two-step coefficient `f*rho'/(1-f*rho')` is 0.16464735, and the one-step coefficient `f/(1-f)` is 0.68932938, where `f=w+(1−w)*rho'`. For a single late window with k assumed contractions, use `rho'^k/(1−rho'^k)` against its contrast with the early anchor. Use the final population-weight uncertainty and the nominated contrast when quoting numerical confidence bounds.

Applying the same exploratory joint-contraction diagnostic from A with the weaker factor gives a simultaneous upper value of about `4.60e-8`, above the data lower endpoint `2.40e-8`. That particular diagnostic does not reject the joint weaker assumption.

### What the scan does not establish

- The established uniform renewal bound cited in Astra's proof is for Re(s)>0.51. The extension to about 0.20–0.21 is a numerical tail-rate estimate. At Re(s)=0.30 the measured real last-interval ratio is 0.9074, giving an effective exponent 0.2029. Completing the tail changes the leading real eigenvalue from 3.9406 to 5.2965. The two sampled intervals do not bound the omitted tail, and a single geometric fit shared by all entries is an approximation.
- The complex evaluator lumps large marks. It divides the **truncated** row-a=0 moment by the **tail-corrected** total mass and uses that same conditional mark for every row. The desired interpolated linear operator would integrate the mark weight in each row, with a matched tail if included. The closures I tried do not materially change the observed grid minima. None supplies an operator error bound. Astra's exact capped-label representation would remove the high-degree column interpolation if used.
- There are only seven real parts and forty frequencies, with frequency steps as large as ten. A grid-point separation cannot exclude zeros between nodes or above |Im s|=30. Non-normality matters because closeness of eigenvalues alone does not bound the perturbation of invertibility. A quantitative exclusion would require operator/norm error control, a mesh or winding argument with between-node bounds, and high-frequency control. The batch jackknife at the real leading value measures sampling variability there. It does not control these other errors.
- Even with a certified pole-free region, one must control the initial source and the transformed scores of the particular observables, as well as the resolvent and its inverse-transform remainder. A spectral gap would then motivate an asymptotic **envelope** for the relative error. An envelope places no bound on the ratio of two successive magnitudes when fast modes cancel or oscillate. This issue remains at factor 0.3465, just as it does at 0.104. The steps used by the finite-window claim are the assumptions themselves.

I did not read either newly registered A-test. The earlier verdicts concern the reconstructed pre-registration evidence only.

```text
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 ../../.venv/bin/python -u -B task3/spectrum.py
T=17.0, complete 10 x 2000, rows=49 bins=44, 7 sigma by 40 omega, max omega=30
sigma grid [0.3       0.35      0.4       0.45      0.5       0.55      0.5657237]
real-axis results, output ordering lead eigenvalue and closest distance from 1
  sigma=0.3000000 q=0.90743 alpha_eff=0.20286 second |eig|=0.51251 sv_min(I-K)=0.307175 supplied:5.29647 truncated:3.94064 consistent:5.28128 midpoint:5.27285 lo:5.23721 hi:5.30851
  sigma=0.3500000 q=0.86316 alpha_eff=0.20284 second |eig|=0.40399 sv_min(I-K)=0.366292 supplied:3.11358 truncated:2.74859 consistent:3.10987 midpoint:3.10679 lo:3.0899 hi:3.12372
  sigma=0.4000000 q=0.82104 alpha_eff=0.20282 second |eig|=0.32521 sv_min(I-K)=0.276931 supplied:2.13183 truncated:2.02223 consistent:2.13081 midpoint:2.12954 lo:2.12054 hi:2.13858
  sigma=0.4500000 q=0.78098 alpha_eff=0.20279 second |eig|=0.26555 sv_min(I-K)=0.175587 supplied:1.59631 truncated:1.56154 consistent:1.59601 midpoint:1.59545 lo:1.59033 hi:1.60059
  sigma=0.5000000 q=0.74287 alpha_eff=0.20277 second |eig|=0.21925 sv_min(I-K)=0.0885418 supplied:1.26846 truncated:1.25706 consistent:1.26837 midpoint:1.26812 lo:1.26508 hi:1.27117
  sigma=0.5500000 q=0.70663 alpha_eff=0.20275 second |eig|=0.18272 sv_min(I-K)=0.0182944 supplied:1.05123 truncated:1.04739 consistent:1.0512 midpoint:1.05108 lo:1.04921 hi:1.05295
  sigma=0.5657237 q=0.69560 alpha_eff=0.20274 second |eig|=0.17285 sv_min(I-K)=0.000842102 supplied:0.997691 truncated:0.994958 consistent:0.997671 midpoint:0.997576 lo:0.99597 hi:0.999187
number of real s/bin weights w outside [0,1] for supplied closure: 0 examples: []
complex-grid comparison (excluding the leading real candidate at max sigma):
  supplied: min eigenvalue distance=(np.float64(0.05122551018054455), (0.55, 0.0, np.complex128(1.0512255101805446+7.619900673367022e-76j))), min singular value=(np.float64(0.018294400629176217), (0.55, 0.0))
  supplied closure: |q|>=1 at 0 grid points; max=(np.float64(0.9126612090379601), (0.3, 15.0)); examples=[]
  truncated: min eigenvalue distance=(np.float64(0.04739169795213094), (0.55, 0.0, np.complex128(1.047391697952131+4.347800097999703e-78j))), min singular value=(np.float64(0.016947448339087807), (0.55, 0.0))
  consistent: min eigenvalue distance=(np.float64(0.05119766525494107), (0.55, 0.0, np.complex128(1.051197665254941-3.568127346733964e-79j))), min singular value=(np.float64(0.018286293660139955), (0.55, 0.0))
  midpoint: min eigenvalue distance=(np.float64(0.051076770531873095), (0.55, 0.0, np.complex128(1.051076770531873-3.4003842800783544e-71j))), min singular value=(np.float64(0.018252100289370403), (0.55, 0.0))
at sigma=0.5657237, omega=0 nearest eigenvalue=(0.9976908682106207+0j) complex jackknife RMS=0.0050559
at sigma=0.5657237, omega=0.05 nearest eigenvalue=(0.9713068199349097-0.15711094217482258j) complex jackknife RMS=0.00506592
weaker premise rho=0.34645581033; f=0.40804912604, f*rho=0.14137099062, normalized coefficients 0.16464734952 and 0.68932937512
```
