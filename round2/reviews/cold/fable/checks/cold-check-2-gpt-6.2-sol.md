VERDICT: GAP

The deterministic inequalities in (d) are valid with positive population denominator increments and the stated contractions for the full observable. The reviewed draft still needs a confirmatory data split for nominal statistical coverage, outward factors for strict upper-bound wording, and clearer labels for calibration results. The same-code benchmarks and observed friend-tree contrasts can support a conditional inference. They do not establish the contraction premise.

Reviewed copy: `followup/CLAIM.reviewed.md`, SHA-256 `606ef97a44deab67e0ea388ce0ea72ea750716dd2eeaa23a9eb65dfeba4f8a3b`. I used the archived 05:10 UTC reduced snapshot for the numbers in that copy. The independent raw-file audit reads only the snapshot's completed files. All new artifacts are under this check directory. I did not rebuild the production binary.

## (b) Statistical statement and provenance

The supplied reducer/final estimator structure preserves checkpoint pairing and restricts cross-observable comparisons to common chunks. I replayed the published evaluator on a copy of the archived reduction. My separate script checked checkpoint ordering, shape, expected seed/run IDs, all stored sums, finite values and duplicate initialized RNG states, then recalculated estimates from the raw records. All 158 archived chunks passed. There were 222,100 distinct production initial states and a largest archive-vs-recalculation discrepancy of `5.5e-16` after scaling by `1+abs(sum)`.

The 81-chunk main result is `0.5657236639625 ± 8.8783e-9` (delete-one-chunk). The asymptotic independent-run standard error is `9.3149e-9` (delete-one-run). The estimated ordinary ratio bias is about `9.0e-13`. The overall estimate matches (b). The statistical SE pertains to the population *finite-window* ratio. At fixed windows with finite second moments and a separated mean denominator, either run-level or equal-size independent-chunk jackknife is asymptotically justified. A 95% confidence statement additionally uses normal/large-sample calibration.

The draft discloses that the estimand was chosen after 600 paths were seen. Reusing those paths gives a descriptive all-run SE with no adjustment for selection. Phase 2 (chunks numbered at least 500) gives a straightforward post-selection holdout for the primary estimate and both independent rungs. On the interim archive the holdout main result is `0.5657236585544 ± 1.1336e-8` (51 chunks, 5,100 runs) or `±1.1760e-8` using independent-run errors. The all-run estimate can appear as a sensitivity check. The separation costs no further simulation. A check on the final raw data would use each shard's expected run IDs and confirm that the holdout begins after selection.

The early independent rung's chunk SE, `1.2610e-8` over 39 blocks, is appreciably smaller than its run SE, `2.1263e-8` over 195,000 paths. The ratio of those values is 0.593. This audit found distinct initial states and correct rows. I would use the direct run-level SE for final independent-run data or show both errors. Equal-size clustering by arbitrary file boundaries discards within-file variance information and can make an interim SE noisy.

```text
$ ../../.venv/bin/python ../../src/prod_final.py followup/interim.tsv
--- main: 81 chunks, 8100 runs, 8.271e+09 core vertices at t=24
       12 -> 16  | plain: 0.5657452700 +- 3.3e-05 | sane: 0.5657237057 +- 1.4e-07 | full: 0.5657237100 +- 1.2e-07 | AN: 0.5657230279 +- 1.5e-06 [51] | ANPQ: 0.5657235623 +- 4.0e-07 [51] | AR: 0.5657252960 +- 8.4e-07 [51]
       16 -> 20  | plain: 0.5657258232 +- 1.0e-05 | sane: 0.5657236593 +- 3.7e-08 | full: 0.5657236724 +- 3.0e-08 | AN: 0.5657245676 +- 3.8e-07 [51] | ANPQ: 0.5657237787 +- 9.6e-08 [51] | AR: 0.5657237342 +- 3.1e-07 [51]
       16 -> 24  | plain: 0.5657258617 +- 3.0e-06 | sane: 0.5657236611 +- 1.1e-08 | full: 0.5657236640 +- 8.9e-09 | AN: 0.5657235463 +- 1.1e-07 [51] | ANPQ: 0.5657236138 +- 2.9e-08 [51] | AR: 0.5657235218 +- 9.1e-08 [51]
--- r1216: 39 chunks, 195000 runs, 2.145e+09 core vertices at t=16
       12 -> 16  | plain: 0.5657259925 +- 7.1e-06 | sane: 0.5657236440 +- 1.7e-08 | full: 0.5657236860 +- 1.3e-08 | AN: 0.5657245126 +- 2.1e-07 [24] | ANPQ: 0.5657238145 +- 8.9e-08 [24] | AR: 0.5657236145 +- 1.7e-07 [24]
--- r1620: 38 chunks, 19000 runs, 2.007e+09 core vertices at t=20
       16 -> 20  | plain: 0.5657275917 +- 6.4e-06 | sane: 0.5657236526 +- 2.5e-08 | full: 0.5657236605 +- 2.1e-08 | AN: 0.5657234197 +- 2.9e-07 [23] | ANPQ: 0.5657234819 +- 7.9e-08 [23] | AR: 0.5657236554 +- 2.3e-07 [23]

$ ../../.venv/bin/python -B followup/recompute.py
158 archived chunks validated; 222100 distinct production initial states; max scaled recalc discrepancy=5.5e-16
validated counts by kind and phase2: {('main', False): 30, ('main', True): 51, ('r1216', False): 15, ('r1216', True): 24, ('r1620', False): 15, ('r1620', True): 23}
phase2-only=False
  main full 16:24 estimate=0.5657236639625 chunkSE=8.8782554e-09 runSE=9.3148695e-09 ratioBias=+9.009e-13 chunks=81
  main full 16:20 estimate=0.5657236723985 chunkSE=2.966398e-08 runSE=3.3612728e-08 ratioBias=+3.519e-12 chunks=81
  main full 20:24 estimate=0.5657236630847 chunkSE=1.0886607e-08 runSE=1.0685023e-08 ratioBias=+6.285e-13 chunks=81
  r1216 full 12:16 estimate=0.5657236860494 chunkSE=1.2609781e-08 runSE=2.1263063e-08 ratioBias=-9.21e-14 chunks=39
  r1216 sane 12:16 estimate=0.5657236440458 chunkSE=1.6572323e-08 runSE=2.5546549e-08 ratioBias=-1.54e-13 chunks=39
  r1620 full 16:20 estimate=0.5657236605194 chunkSE=2.1182835e-08 runSE=2.1968648e-08 ratioBias=+4.847e-13 chunks=38
  w=0.09424506665 SE=2.280096e-06 rho(at main plug-in)=0.10404886277
  all-run mean den/core by time [np.longdouble('0.86532846271932066577'), np.longdouble('0.8653745826252739753'), np.longdouble('0.86536840550650898714'), np.longdouble('0.8653678029640806462')]
  f=0.1884878374 f*rho=0.0196119451 k2=0.0200042677 k1=0.2322674214
  two-step: rung-main=+2.2086876e-08 SE=1.5421738e-08; plug-in upper95=1.0464818e-09; Bonferroni upper95=1.1383622e-09 (w upper=0.0942495356, rho upper=0.1043504848)
  one-step: rung-main=-3.4430326e-09 SE=2.296815e-08; plug-in upper95=1.1255628e-08; Bonferroni upper95=1.2780161e-08 (w upper=0.0942495356, rho upper=0.1043504848)
  paired r1216 full-sane=+4.20035384e-08 SE=1.4680532e-08
  full-sane: +2.8206525e-09 +- 6.441973e-09 (81 common chunks); sane/full SE ratio=1.2878, unnormalized run residual SD ratio=1.1043, normalized=1.2206
  full-AN: +1.1225808e-07 +- 1.112824e-07 (51 common chunks); AN/full SE ratio=9.9042, unnormalized run residual SD ratio=10.032, normalized=9.8925
  full-ANPQ: +4.4794375e-08 +- 2.389761e-08 (51 common chunks); ANPQ/full SE ratio=2.5277, unnormalized run residual SD ratio=2.445, normalized=2.4335
  full-AR: +1.3678381e-07 +- 8.884058e-08 (51 common chunks); AR/full SE ratio=8.0385, unnormalized run residual SD ratio=6.8359, normalized=7.469
phase2-only=True
  main full 16:24 estimate=0.5657236585544 chunkSE=1.1336298e-08 runSE=1.1759735e-08 ratioBias=+2.561e-12 chunks=51
  main full 16:20 estimate=0.5657237283599 chunkSE=3.6706005e-08 runSE=4.2792092e-08 ratioBias=+4.249e-12 chunks=51
  main full 20:24 estimate=0.5657236512912 chunkSE=1.3959528e-08 runSE=1.3554941e-08 ratioBias=+2.388e-12 chunks=51
  r1216 full 12:16 estimate=0.5657236978923 chunkSE=1.4411105e-08 runSE=2.7075135e-08 ratioBias=-1.739e-13 chunks=24
  r1216 sane 12:16 estimate=0.5657236462518 chunkSE=1.9009217e-08 runSE=3.2583953e-08 ratioBias=-3.549e-13 chunks=24
  r1620 full 16:20 estimate=0.5657236504029 chunkSE=2.7261876e-08 runSE=2.8049177e-08 ratioBias=-2.534e-13 chunks=23
  w=0.09424342376 SE=2.894493e-06 rho(at main plug-in)=0.10404886502
  all-run mean den/core by time [np.longdouble('0.8653056230178056715'), np.longdouble('0.86540096318759626926'), np.longdouble('0.86538009674572423044'), np.longdouble('0.8653645923870324236')]
  f=0.1884863675 f*rho=0.0196117926 k2=0.0200041090 k1=0.2322651893
  two-step: rung-main=+3.9337846e-08 SE=1.8335528e-08; plug-in upper95=1.5058057e-09; Bonferroni upper95=1.6162218e-09 (w upper=0.0942490969, rho upper=0.1043504848)
  one-step: rung-main=-8.1515105e-09 SE=2.9524931e-08; plug-in upper95=1.5333988e-08; Bonferroni upper95=1.7295405e-08 (w upper=0.0942490969, rho upper=0.1043504848)
  paired r1216 full-sane=+5.16404657e-08 SE=2.0062441e-08
  full-sane: +5.1767162e-09 +- 8.036056e-09 (51 common chunks); sane/full SE ratio=1.2454, unnormalized run residual SD ratio=1.1158, normalized=1.2333
  full-AN: +1.1225808e-07 +- 1.112824e-07 (51 common chunks); AN/full SE ratio=9.9042, unnormalized run residual SD ratio=10.032, normalized=9.8925
  full-ANPQ: +4.4794375e-08 +- 2.389761e-08 (51 common chunks); ANPQ/full SE ratio=2.5277, unnormalized run residual SD ratio=2.445, normalized=2.4335
  full-AR: +1.3678381e-07 +- 8.884058e-08 (51 common chunks); AR/full SE ratio=8.0385, unnormalized run residual SD ratio=6.8359, normalized=7.469
```

## (c) Meaning of the contraction assumption

Use the inverse growth factor `rho = exp(-4*mu)`. The leading growth factor of an exponential denominator is `exp(4*mu) ≈ 9.61`. Write `D_t = E[den(t+4)] - E[den(t)]` and `A_t = beta(t+4)-beta(t)`. Assumption (C) means `|A_(t+4)| <= rho * |D_(t+4)/D_t| * |A_t|`. Non-growing absolute increments of beta have the same interpretation only when the adjacent D ratio equals `exp(4*mu)`. The observed normalization is a ratio of sample sums estimating `E[den]/E[C]`. An asymptotic limit for that ratio is unnecessary for the finite-window algebra. Positive population increments at 16–20 and 20–24 suffice for the convex weights.

Only the full observable needs contractions at starts 12 and 16 for the two-step conclusion. The one-step conclusion only needs the contraction at 16. Hypotheses on all starts and all five observables make a stronger claim than the bound uses.

The archived normalization and the observed four-time-unit ratios of mean denominator increments are below. The measured growth is close to the hypothesized scale. It still leaves the behavior of the residual increments open.

```text
t=12 sane den/C=0.7830386342 full den/C=0.8653284627
t=16 sane den/C=0.7829341568 full den/C=0.8653745826
t=20 sane den/C=0.7829133501 full den/C=0.8653684055
t=24 sane den/C=0.7829112387 full den/C=0.8653678030
sane successive window-denominator growth factors 9.6098770384, 9.6106352076
full successive window-denominator growth factors 9.6098771431, 9.6106349702
exp(4*mu_ref)=9.6108690992
```

## (d) Exact coefficients and inference

For full, set `R0 = R(12,16), R1 = R(16,20), R2 = R(20,24)`, `R = R(16,24)`, `w = D16/(D16+D20)`, and `f = w+(1-w)*rho`. Summing the two numerator and denominator increments gives `R = w*R1 + (1-w)*R2`. Under the contraction at 16, `|R-mu| <= f*|R1-mu|`. With the additional contraction at 12 this is at most `f*rho*|R0-mu|`. Using `Rj-mu=(Rj-R)+(R-mu)` and rearranging gives the two deterministic formulas.

```text
two-step coefficient = (f*rho)/(1-f*rho)
one-step coefficient = f/(1-f)
```

The archive gives `w_hat=0.09424506665 ± 0.00000228010` and the exponent plug-in gives `rho=0.10404886277`. I obtain `f=0.1884878374`, `f*rho=0.0196119451`, two-step coefficient `0.0200042677`, and one-step coefficient `0.2322674214`. Thus 0.1885 is a good outward illustration at these plug-in values. The displayed 0.0196, 0.0200 and 0.232 round downward. Use the exact expressions or round upward when writing an inequality.

The independent two-step rung-minus-main contrast is `+2.20869e-8 ± 1.54217e-8` using chunk SEs. With a fixed plug-in coefficient this gives `1.0465e-9` as an approximate conditional 95% upper value. The one-step contrast is `-3.44303e-9 ± 2.29682e-8` and gives `1.1256e-8` on the same basis. Strictly, `w` is a population quantity. One way to cover its estimation is a 97.5% one-sided upper limit for w and a 97.5% two-sided interval for the independent rung contrast, with Bonferroni coverage even though the main data are shared. The proved lower bracket 0.5650 gives a conservative upper rho of 0.104350485. With those choices I obtain conditional upper values of `1.1384e-9` and `1.2780e-8` from the all-run archive or `1.6162e-9` and `1.7295e-8` from phase 2 alone, using chunk SEs for this replay. The final claim should label an approximate confidence upper bound conditional on (C) and use one nominated contrast, or specify joint coverage if using multiple options for selection.

Using the direct run jackknife for the independent contrast and the same Bonferroni treatment of w gives a more conservative interim comparison, especially for the earlier rung.

```text
all two-step: using independent-run contrast SE=2.3213889e-08, conditional Bonferroni upper95=1.489303e-09
all one-step: using independent-run contrast SE=2.3861858e-08, conditional Bonferroni upper95=1.3246289e-08
holdout two-step: using independent-run contrast SE=2.9518711e-08, conditional Bonferroni upper95=2.1198863e-09
holdout one-step: using independent-run contrast SE=3.0414597e-08, conditional Bonferroni upper95=1.7759405e-08
```

## (e) Evidence checked against the data

The same-weight friend-tree early rungs reproduce the reported estimates. Using the exact archived main estimate as the working reference, the full offsets are `+1.06066e-4` and `+1.95599e-6` (plug-in shrink factor 54.2). Sane gives `-7.40804e-5` and `-6.26842e-7` (118.2). Small differences from the draft's factors 55 and 116 come from the reference value and rounding. Applying the single-geometric-mode extrapolation `b(12)=b(8)^2/b(4)` gives `+3.6071e-8` for full and `-5.3041e-9` for sane. The predicted difference `+4.1375e-8` agrees with the observed paired `+4.2004e-8 ± 1.4681e-8`. The measured contrast has a common-mode ambiguity and the extrapolation is a model for the trend. The first observed contraction starts at t=4, outside the draft's t>=8 domain.

IR-0 reproduces the late agreement with the independently verified root. The full-precision offsets relative to that root are `+1.29082356e-4`, `+1.36300336e-6`, `+7.47255e-9 ± 4.74301e-9` and `-6.256e-11 ± 1.44326e-9`. The last matched contrast between the 20,000-path late windows is `-7.5351e-9 ± 5.3077e-9`. The last two windows cannot establish a shrink factor greater than five. The penultimate offset's approximate 95% interval includes zero, and the latest shift is compatible with zero. State the offsets and SEs instead of dividing noisy estimates near zero. The first two clearly resolve a large shrink.

The existing cold PAN and DAN estimates reproduce their reported offsets. Running the lead's archived PAN vector data through the fitter reproduces `0.549734977 ± 6.73e-9` as a five-fold training cross-fit. The archive contains those vector records and a weight file. Unless a distinct lead evaluation dataset is supplied, label that number a cross-fit and give the cold PAN result as the held-out validation. The fit's printed jackknife holds the fold weights fixed, so it lacks the fixed-weight independent-production justification.

All four main paired observable differences reproduce. The quoted residual multipliers need a definition and provenance. On the identical 51-chunk phase-2 sample, AN, ANPQ and AR have chunk jackknife SE ratios to full of 9.90, 2.53 and 8.04. Their run-level normalized residual SD ratios are 9.89, 2.43 and 7.47. These measure estimation noise, not the amplitudes of systematic modes. If choosing a raw, unnormalized run residual metric, the ratios are 10.03, 2.45 and 6.84.

Astra's held-out t=24 midpoint is `0.56572369741 ± 2.62e-8` (2,048 paths). The current Astra t=28 midpoint report shows `0.56572364724 ± 2.62e-8` (208 paths). The draft's slightly different t=28 value is in the Fable board summary, and I did not locate its original output. Use the current Astra number or give the timestamp and original file for the older value. Astra's own report also flags its inherited high-intensity sampler approximation and describes the rates as finite-time corroboration.

```text
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/rungs_early/r04_08.txt
runs=2000000 checkpoints=[np.float64(4.0), np.float64(8.0)]  mean C at checkpoints: [np.float64(12.2), np.float64(119.1)]
window   4.0 -> 8     | plain: 0.565788202 +- 2.01e-05 | th1: 0.565649584 +- 2.89e-07 | th2: 0.565829730 +- 5.52e-07
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/rungs_early/r08_12.txt
runs=300000 checkpoints=[np.float64(8.0), np.float64(12.0)]  mean C at checkpoints: [np.float64(119.1), np.float64(1146.3)]
window   8.0 -> 12    | plain: 0.565721811 +- 1.67e-05 | th1: 0.565723037 +- 6.44e-08 | th2: 0.565725620 +- 6.32e-08
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/valid/ir0_early.txt
runs=400000 checkpoints=[np.float64(6.0), np.float64(10.0), np.float64(14.0)]  mean C at checkpoints: [np.float64(27.2), np.float64(189.3), np.float64(1305.7)]
window   6.0 -> 10    | plain: 0.482333116 +- 4.18e-05 | th1: 0.482396001 +- 8.34e-07
window   6.0 -> 14    | plain: 0.482272984 +- 1.48e-05 | th1: 0.482284481 +- 1.06e-07
window  10.0 -> 14    | plain: 0.482264249 +- 1.59e-05 | th1: 0.482268281 +- 2.56e-08
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/valid/ir0_ladder.txt
runs=20000 checkpoints=[np.float64(10.0), np.float64(14.0), np.float64(18.0), np.float64(22.0)]  mean C at checkpoints: [np.float64(189.1), np.float64(1304.6), np.float64(8982.5), np.float64(61829.3)]
window  10.0 -> 14    | plain: 0.482194684 +- 7.08e-05 | th1: 0.482268246 +- 1.18e-07
window  10.0 -> 18    | plain: 0.482261137 +- 2.53e-05 | th1: 0.482267093 +- 1.51e-08
window  10.0 -> 22    | plain: 0.482264525 +- 9.54e-06 | th1: 0.482266943 +- 2.40e-09
window  14.0 -> 18    | plain: 0.482270792 +- 2.73e-05 | th1: 0.482266926 +- 4.74e-09
window  14.0 -> 22    | plain: 0.482265812 +- 9.65e-06 | th1: 0.482266919 +- 1.25e-09
window  18.0 -> 22    | plain: 0.482265089 +- 1.04e-05 | th1: 0.482266918 +- 1.44e-09

$ ../../.venv/bin/python -B followup/ladders.py
sane 4:8 n=2000000 ratio=0.5656495835686 off=-7.40803939e-05 SE=2.8902026e-07
full 4:8 n=2000000 ratio=0.5658297300125 off=+0.00010606605 SE=5.5217664e-07
sane 8:12 n=300000 ratio=0.5657230371202 off=-6.26842328e-07 SE=6.4388242e-08
full 8:12 n=300000 ratio=0.5657256199486 off=+1.95598608e-06 SE=6.3231911e-08
sane ref=0.5657236639625 plug-in abs contraction factor=118.1803 model b(12)=b(8)^2/b(4)=-5.3041201e-09
full ref=0.5657236639625 plug-in abs contraction factor=54.22638 model b(12)=b(8)^2/b(4)=+3.6070746e-08
IR-0 6:10 n=400000 ratio=0.4823960006795 off=+0.000129082356 SE=8.3435245e-07
IR-0 10:14 n=400000 ratio=0.4822682813267 off=+1.36300336e-06 SE=2.5611075e-08
IR-0 14:18 n=20000 ratio=0.4822669257959 off=+7.47254864e-09 SE=4.7430082e-09
IR-0 18:22 n=20000 ratio=0.4822669182608 off=-6.25580776e-11 SE=1.4432649e-09
IR-0 paired last contrast late minus early=-7.5351067e-09 SE=5.3076633e-09

+ ../../.venv/bin/python ../../src/fit_theta.py ../../outputs/valid/pan followup/pan-refit.theta --mu0 0.5497349725
training runs=30000 layout p=3189 equations kept q=1971 (minfrac=0.01)
plain (differenced)      : 0.549799924 +- 3.73e-05
cross-fitted GMM (5 folds): 0.549734977 +- 6.73e-09   se ratio vs plain = 5535.61
theta (all runs) written to followup/pan-refit.theta
+ sha256sum followup/pan-refit.theta ../../outputs/valid/pan.theta.bin
c8b2c7f9c489a8a4eea01822a75ffc1d6d5d6fa99c5320ceab726709ed47dcaf  followup/pan-refit.theta
15872c01bd70691ec29c9c3fee3ad81f13809742117ba8dc4691491b1100fed3  ../../outputs/valid/pan.theta.bin

$ rg "t=28 n=208 midpoint" ASTRA_REPO/round2/outputs/numeric-report.log
t=28 n=208 midpoint   0.56572364724 +- 2.62e-08 raw_SE/SE=243.4
$ rg "t=24 independent" ASTRA_REPO/round2/outputs/numeric-confirmation.log
t=24 independent runs=2048 midpoint=0.56572369741 +- 2.62e-08
```

## (f) Scope and additional evidence

A slow shared relative correction can make the friend-tree contrasts small while keeping the primary finite-window ratio away from the asymptotic exponent. It must vary slowly enough across checkpoints, too. Changes across time and independent observables can reveal it when their sensitivities or precision suffice. A second failure mode is cancellation. Even if two individual modes contract faster than rho, their sum need not contract at that rate. For example, `b_j = A*(r^j-s^j)` with `0<s<r<=rho` has zero at j=0 and nonzero magnitude at j=1. This is why (C) has to constrain the total bias, and why the observed late agreement in the calibration models gives limited evidence about a different process.

For the conditional inference itself I would require the final frozen, post-selection scalar estimate and independent rung checks, appropriate variance estimates and a named confidence procedure. Existing production outputs are sufficient to do those checks without another simulation. The final manifest should include model, executable and weight hashes, launch seeds, run ranges and the actual checkpoint layout for every included shard. This also permits replay of the five-observable controls. Sampling SE and conditional systematic upper value should remain separate in the write-up unless a joint interval is explicitly constructed.

The most informative additional Monte Carlo sensitivity check within the remaining budget would be a preregistered late-only independent cohort through t=28, with matched 20, 24 and 28 checkpoints. Evaluate the same frozen full and sane observables and, if available, an independently constructed late observable. Report the paired `R(24,28)-R(20,24)` and the independent comparison with the production main window. Based on the observed t=24 throughput and about 9.6 times more core work at t=28, nine hours on one or two cores would plausibly give order 2,000–4,500 late paths, subject to a timing pilot. This would probe discrepancies at a new timescale. A finite ladder still cannot exclude an arbitrarily slow common correction. Such a check would strengthen the numerical evidence for (C), rather than supply an unconditional proof.

## Independent deterministic IR-0 root

Write `d=L+c` and `b(d,c)=(d-c)/d`. Model 3 discards core-to-core attempts, so the two effective transitions of one core are `(d,c)->(d+1,c)` at rate `d-c` and `(d,c)->(d,c+1)` at rate `b(d,c)`. The second gives a newborn at (d,c)=(2,1). This follows directly from the effective state transitions in the simulator and generator and agrees with the chain supplied in the brief. The initial root's state does not change the newborn Malthusian equation.

I computed the discounted expected birth count `H_mu(d,c)=E integral exp(-mu*t) b(D(t),c(t)) dt`, with `0<=H_mu<=1/mu`. The first-step recurrence for c<d is

```text
H_mu(d,c) = [d*H_mu(d+1,c) + 1 + H_mu(d,c+1)] / [d+1 + mu*d/(d-c)].
```

For c=d the state is absorbing and H is zero. The newborn function `H_mu(2,1)` is continuous and strictly decreasing for mu>0. The initial birth intensity is positive and the Laplace transform weights positive times strictly less as mu increases. The two certified signs below thus bracket a unique root of `H_mu(2,1)=1`.

I used a degree-boundary trial `F(d,c)=1/mu + sum_{j=1..K} q_j(c)/d^j`. Each q is solved by finite polynomial differences to cancel the coefficients through degree K in `x=1/d` of `(mu-G)F-b`. At order j the linear step is `mu+j-Delta_c`, with a finite Neumann inverse on polynomials. For each tested rational mu I computed the cleared residual polynomial exactly using Python Fraction. The numerator after multiplying by `(1+x)^K` has the form `sum r_ij c^i x^j` with `j>=K+1` and `i<=K+1`. With `p=K+1, a=3p` and degree boundary D, its magnitude in the region d>=D is bounded by `M W(d,c)`, where

```text
W(d,c) = ((c+a)/d)^p
M = sum |r_ij| * D^(-(j-p)) * a^(i-p).
```

The function W has nonpositive generator. At c<d, divide the generator by `(d-c)*W` to get the sum of `(d/(d+1))^p-1` and `[(1+1/(c+a))^p-1]/d`. The first term is at most `-1/(d+1)` and the second is at most `(exp(1/3)-1)/d`. The total is nonpositive for d>=1, and the absorbing case has zero generator. Thus `E W(X(t))<=W(X(0))`. The chain is nonexplosive by domination of degree growth by a Yule chain. The trial F is bounded on d>=D, c<=d. Discounted Dynkin integration bounds the error in this region by `(M/mu)*W`.

I put that verified interval on the degree boundary. At the separate c cut I used the unconditional range [0,1/mu], so all potential visits beyond the numerical state cap are included. I propagated lower and upper values with 75-place directed Decimal operations in the positive recurrence. This is a deterministic root enclosure with controlled truncation and numerical rounding. The leading digits agree with the board reference. No Neon solver or code was used.

The degree-16,384 check alone brackets the root in an interval of width `1e-48`. The boundary width in H from cutoff and rounding is about `2.01e-52` at both endpoints.

```text
+ python3 followup/ir0_cert.py --degree 16384 --cap 96 --order 20 --certify 0.482266918323368718268897622913549185185240982589 0.482266918323368718268897622913549185185240982590
D=16384 c_cap=96 expansion_order=20
mu=0.482266918323368718268897622913549185185240982589: H-1 in [+7.3010035804483475097377107000000000000e-49, +7.3030088858605338027197897000000000000e-49] width=2.0053054e-52; residual scale M<=0.221046
mu=0.482266918323368718268897622913549185185240982590: H-1 in [-1.4870186245865906563498165730000000000e-48, -1.4868180940453720270516086630000000000e-48] width=2.0053054e-52; residual scale M<=0.221046
CERTIFIED: 0.482266918323368718268897622913549185185240982589 < mu_IR0 < 0.482266918323368718268897622913549185185240982590
```

At the deeper boundary I enclosed the root in an interval of width `1e-60`. Both endpoint signs were certified by the degree-32,768 resolvent. I also verified the upper endpoint independently using the lower-endpoint value and `-H'(mu) = integral t exp(-mu*t) E[b(X(t))] dt` and restricted the integral to time spent at degrees 2, 3 and 4 before any conversion. On the interval in question the remaining integrand is nonnegative and the restricted derivative magnitude is bounded below by its exact rational value at mu=1/2. In the restricted calculation each visited degree contributes its discounted entrance factor times its discounted holding intensity. Differentiating that product gives the sum of inverse waiting denominators implemented in `ir0_one_sided.py`. The strict upper sign follows from the exact slope bound and the lower-endpoint residual enclosure.

```text
+ python3 followup/ir0_cert.py --degree 32768 --cap 112 --order 24 --certify 0.482266918323368718268897622913549185185240982589329346611091 0.482266918323368718268897622913549185185240982589329346611092
D=32768 c_cap=112 expansion_order=24
mu=0.482266918323368718268897622913549185185240982589329346611091: H-1 in [+2.5841825163855000000000000000000000000e-61, +2.5841834046368000000000000000000000000e-61] width=8.882513e-68; residual scale M<=0.203206
mu=0.482266918323368718268897622913549185185240982589329346611092: H-1 in [-1.9587007309928880000000000000000000000e-60, -1.9587006421677300000000000000000000000e-60] width=8.8825158e-68; residual scale M<=0.203206
CERTIFIED: 0.482266918323368718268897622913549185185240982589329346611091 < mu_IR0 < 0.482266918323368718268897622913549185185240982589329346611092
$ python3 followup/ir0_one_sided.py
partial exact slope at mu=1/2 = 224889/834632 = 0.26944689396045
H(lower)-1 lower > 2.584182516385500e-61
H(upper)-1 upper <= -1.102855349677203e-62
CERTIFIED: 0.482266918323368718268897622913549185185240982589329346611091 < mu_IR0 < 0.482266918323368718268897622913549185185240982589329346611092
Confirmed 1e-60 bracket using one resolvent evaluation plus the exact slope lower bound.
```

## Frozen-binary provenance spot check

I replayed four paths from the first five-observable main chunk and four from each of the first five-observable rung chunks. I used the specified frozen binary with the archived weights and expected seeds and run starts. All 32 checkpoint rows matched bit for bit, including all five numerator/denominator pairs. The main replay runs through t=24. It checks the launch fields and weight order on those raw chunks without recompiling anything.

```text
+ ../../bin/ftprod.frozen-0446 9100 50000 4 followup/replay-main ../../data/theta_sane.f64,../../data/theta_full.f64,../../data/theta_AN.f64,../../data/theta_ANPQ.f64,../../data/theta_AR.f64 12 16 20 24 Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8
+ ../../bin/ftprod.frozen-0446 9200 2500000 4 followup/replay-r1216 ../../data/theta_sane.f64,../../data/theta_full.f64,../../data/theta_AN.f64,../../data/theta_ANPQ.f64,../../data/theta_AR.f64 12 16 Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8
+ ../../bin/ftprod.frozen-0446 9300 250000 4 followup/replay-r1620 ../../data/theta_sane.f64,../../data/theta_full.f64,../../data/theta_AN.f64,../../data/theta_ANPQ.f64,../../data/theta_AR.f64 16 20 Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8
+ ../../.venv/bin/python -B -
main 16 checkpoint rows bitwise equal: True max scaled difference: 0.0
r1216 8 checkpoint rows bitwise equal: True max scaled difference: 0.0
r1620 8 checkpoint rows bitwise equal: True max scaled difference: 0.0
```

The early friend and IR-0 output text gives ordered run IDs and observable values. It contains no seed or binary hashes. The launch metadata would be needed to replay their claimed independence with the same precision. Their reported estimates all reproduce from the supplied files.
