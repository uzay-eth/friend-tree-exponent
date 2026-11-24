# Claim F1 — numerical and conditional; not a theorem

**Author:** Fable lead (Claude Fable 5.2; final production figures and corrections by Claude Opus 5.5). **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** GPT-family cold check and final Astra audit; see [F1](README.md#f1). The finite-window bias conclusions require unproved contraction premises.

The cited generator and ratio proofs are reproduced in [the derivation](estimator-derivation.md). The frozen-weight SHA-256 prefixes refer to the original binary64 bytes, not the public text transcriptions.

*Version 3 (form fixed 06:45 UTC on 24 Nov 2026). The numbers are final: the production runs ended at 15:30 UTC on
24 Nov and `bash src/final.sh` was run at 15:30; the INTERIM values of 06:35 UTC were replaced and nothing else changed.*

## 1. Objects

**Process.** Continuous-time friend tree started from the path on three vertices at t = 0: every vertex rings at
rate 1 and gives a new leaf to a uniformly chosen neighbour. mu is the exponent of the repository's theorem.

**Observables.** F is the feature map of `src/feat2.hpp` with layout keys `Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8`
(5058 features and the plain entry). For a weight vector theta: den(t) = sum over core vertices u of theta.F(u),
num(t) = sum_u theta.a_u, where a_u is the exact generator action. Population quantities:

    D(t1,t2) = E den(t2) - E den(t1),     R(t1,t2) = [E num(t2) - E num(t1)] / D(t1,t2),     b(t) = R(t, t+4) - mu.

**Frozen weight vectors** (sha256 of the files in `data/`): full `7d825f14…`, sane `bc711866…`, AN `ccf909d9…`,
ANPQ `a7cfe2c9…`, AR `f781b2af…`. All were fitted on training streams 811–814 (window 12 -> 16) and never refitted.
Everything below is for theta = full unless another name is given.

**Samples** (independent random streams; one frozen executable, sha256 `8806c87f…`):

| name | stream | checkpoints | what it is |
|---|---|---|---|
| main | 9100, chunk numbers >= 500 | 12, 16, 20, 24 | hold-out for the primary estimate: started 04:30 UTC, after the five weight vectors were frozen and after the primary quantity R_full(16,24) had been fixed (04:13 UTC, on the basis of 600 earlier runs, which are not in this sample) |
| r1216 | 9200, chunk numbers >= 500 | 12, 16 | independent rung |
| r1620 | 9300, chunk numbers >= 500 | 16, 20 | independent rung |
| late | 9500 | 20, 24, 28 | late cohort, quantities registered 06:00 UTC before its first run |
| early | 9410, 9420, 9430, 9440 | 4, 8 and 8, 12 | early rungs |

## 2. Exact identities (proved in RESULT.md, sections 2.1–2.2)

(I1) E num(t) = d/dt E den(t) for every fixed weight vector (Dynkin's formula for the tree chain).

(I2) If D(t1,t2) > 0 then R(t1,t2) - mu = [beta(t2) - beta(t1)] / D(t1,t2), beta(t) = E num(t) - mu E den(t).

(I3) If D(16,20) > 0 and D(20,24) > 0 then R(16,24) = w R(16,20) + (1 - w) R(20,24), w = D(16,20) / D(16,24).

## 3. Statistical statements

They concern the finite-window population ratios R, not mu. They assume independent runs, correct code and
large-sample normality of ratio estimators; "+-" is one standard error from the linearised delete-one-run
jackknife (sqrt(n/(n-1) sum psi_r^2), psi_r = (e_r - delta d_r)/sum d; see `src/final_analysis.py`).

(S1) R_full(16,24) = 0.5657236606 +- 0.0000000062 (main, 19000 runs).

(S2) R_full(12,16) = 0.5657237030 +- 0.000000018 (r1216, 260000 runs);
     R_full(16,20) = 0.5657236424 +- 0.000000019 (r1620, 25500 runs).

(S3) R_full(24,28) = 0.5657236688 +- 0.0000000060 (late, 2840 runs).

(S4) w = 0.0942443 +- 0.0000016 (main).

## 4. Premises — not established, stated so that the conclusions can be read off

With rho = e^{-4 mu} and rho' = e^{-4 (mu - 3/10)}; by the proved bound mu >= 113/200 they satisfy
rho <= 0.1043505 and rho' <= 0.3464559.

(C)   |b(16)| <= rho |b(12)|  and  |b(20)| <= rho |b(16)|.
(C')  the same two inequalities with rho' in place of rho.
(C+)  in addition to (C): |b(24)| <= rho |b(20)|.        (C'+): the same with rho'.

In words: these are finite-window premises on |b|. (C) uses the reciprocal of the four-time-unit exponential growth
factor at rate mu; (C') uses the larger factor rho'. Why 3/10 was chosen is explained in section 6, item 4; that
explanation is a motivation, not part of the premise. Only the observable full is concerned.

## 5. Conclusions

Exact algebra from (I3) and the premises. Standing conditions: D(12,16), D(16,20), D(20,24) > 0 (for K3, K4 also
D(24,28) > 0), so that all ratios are defined and w lies in (0,1); f = w + (1 - w) r for r = rho or rho'.

(K1) Under (C):   |R(16,24) - mu| <= k |R(12,16) - R(16,24)|,   k = f rho / (1 - f rho).
(K2) Under (C'):  the same with k' = f' rho' / (1 - f' rho').
(K3) Under (C+):  |R(24,28) - mu| <= rho^3/(1 - rho^3) |R(12,16) - R(24,28)|;  the coefficient is <= 0.001138.
(K4) Under (C'+): the same with rho'^3/(1 - rho'^3) <= 0.043391.

Numerically k <= 0.020094 and k' <= 0.164649 whenever w <= 0.094249, which is the upper one-sided 97.5% limit of
(S4) (0.094249 is a conservative cap; the final calculated limit is 0.0942474, rounded up); k and k' increase with w and with rho. All coefficients are rounded upward.

(Proof of K1: b = R(16,24) - mu = w b(16) + (1-w) b(20), so |b| <= f |b(16)| <= f rho |b(12)|, and
b(12) = [R(12,16) - R(16,24)] + b. K3: |b(24)| <= rho^3 |b(12)| and b(12) = [R(12,16) - R(24,28)] + b(24).)

With the measured contrasts (independent samples; all bounds rounded upward; each line is an approximate 95% statement on its own, not simultaneously with the others: a one-sided 97.5% limit for w and a
two-sided 97.5% interval for the contrast, i.e. the bound is coefficient x (|estimated contrast| + 2.2414 standard
errors of the contrast)):

    R(12,16) - R(16,24) = +4.2e-8 +- 1.9e-8
    under (C):   |R_full(16,24) - mu| <= 1.73e-9
    under (C'):  |R_full(16,24) - mu| <= 1.42e-8
    R(12,16) - R(24,28) = +3.4e-8 +- 1.9e-8
    under (C+):  |R_full(24,28) - mu| <= 8.83e-11
    under (C'+): |R_full(24,28) - mu| <= 3.37e-9

The nominated statement is K1/K2 for the primary quantity. K3/K4 are given for the late cohort, which was added
later; its quantities were registered before its first run.

## 6. Evidence about the premises (not part of the claim)

For:

1. Early rungs of the same frozen observable (independent runs): R - 0.56572365 = +1.06e-4 at 4 -> 8 and +1.97e-6
   at 8 -> 12 (ratio 54); for sane -7.4e-5 and -6.1e-7 (ratio 120). The first contraction starts at t = 4, outside
   the range of the premises. On unit windows (600,000 runs on a fresh stream, declared before its first run) the
   correction of full is positive at every window from 4 -> 5 to 12 -> 13 and the ratio b(t+4 -> t+5)/b(t -> t+1) is
   0.029(1), 0.018(1), 0.013(1), 0.008(3), 0.014(6) at t = 4, ..., 8; for sane it is 0.0077(6) at t = 4 and the same
   within errors later (`data/fine_ladder.txt`). Premise (C) allows 0.104. The premises concern t = 12, 16, 20,
   where the correction is at or below the noise.
2. Known-truth models, same code path, weights frozen before the runs. IR-0 (exponent certified to 60 digits by the
   cold check): offsets +1.2908e-4 +- 8.3e-7, +1.363e-6 +- 2.6e-8, +7.5e-9 +- 4.7e-9, -0.6e-10 +- 1.4e-9 at 6 -> 10,
   10 -> 14, 14 -> 18, 18 -> 22. DAN (cold check's own runs): +6.4e-8, +0.75e-8, -0.03e-8 (+- 0.3e-8) at 9 -> 12,
   12 -> 14, 14 -> 16. PAN, held-out run of the cold check: +1.75e-8 +- 0.96e-8 at 12 -> 15.
3. The solvable cousin IR-0: its characteristic equation has exactly one zero with real part >= -0.551315 (the
   exponent) and exactly three with real part >= -9/16 (the exponent and a pair near -0.551317 +- 0.501528 i), by an
   interval-arithmetic argument-principle program of a child (GPT family; rerun by the lead, not re-derived by a
   second family); a floating-point scan had found the same pair. Gap 1.03 per unit time. This says nothing about
   the size of the components in a given observable, and nothing about the friend tree itself.
4. The renewal kernel of the friend tree itself: a Monte Carlo evaluation of the next-generation kernel K_s of the
   Astra lead's renewal characterisation on the grid Re s = 0.30, ..., 0.5657, 0 <= Im s <= 30 resolves no further
   eigenvalue near 1 (at the known real root the measured leading eigenvalue is 0.9977 +- 0.0051). This is the
   motivation for the number 3/10 in (C'). The scan cannot exclude crossings between its 280 grid points, at higher
   frequencies, or hidden by tail and closure errors (the geometric tail extrapolation is large at Re s = 0.30); the
   renewal representation is proved by its author only for Re s > 0.51; and a pole-free strip would give an
   asymptotic envelope, not the finite-window inequalities (C').
5. Agreement of the late cohort with the primary: R(24,28) - R(16,24) = +8.2e-9 +- 8.6e-9. This estimates (with a
   standard error; it is not a deterministic bound) the change of the finite-window value between the two windows,
   not the distance to mu.

Against, or unresolved:

6. Control observable ANPQ. On the first 12000 hold-out runs R_full(20,24) - R_ANPQ(20,24) = +6.6e-8 +- 2.0e-8
   (3.4 standard errors), with the opposite sign at 12 -> 16 (-7.3e-8 +- 4.1e-8). If (C) held for both observables
   the difference at 20 -> 24 would be bounded by about 4e-9 (cold check's inequality, itself an approximate
   simultaneous calculation from noisy contrasts); so either this is a chance excursion or (C) fails for at least
   one of the two. As an illustration only, not an exact false-positive probability for the exploratory search: for
   a correlated family of about 30 Gaussian contrasts under an all-zero null, including a second look, the largest
   one exceeds the observed value in about 2% of cases. The weaker joint premise (C') is not contradicted. Two tests of this contrast on runs that did
   not exist when it was found were registered at 06:25 UTC; both came out "not reproduced" (-0.2 and +1.2 standard
   errors; RESULT.md section 3.3; at the registration 120 runs of the late cohort already existed, and without them
   the second test gives +1.3, also not reproduced). Treating the original excursion as chance is a judgement.
7. In the sample r1216 the sums of the residuals over blocks of 5000 consecutive runs vary less than under random
   reordering (variance ratio 0.43 for the standardised sum of the sane and full residuals, permutation p = 0.001,
   found post hoc). No mechanism was found by the lead or by the cold check. A replication on a fresh stream
   (72 blocks of 5000 runs, statistic and rule registered at 06:31 UTC before its first run) gave a ratio of 1.017,
   p = 0.60: not reproduced. The deviation is treated as a chance excursion of that sample; the larger, run-level
   standard errors are used throughout in any case.

## 7. What the claim does not cover

* A component of the mean dynamics that decays more slowly than the premises allow, or several components that
  cancel at one window and not at the next, would violate a premise and remove the corresponding guarantee (the
  bound could still happen to hold) while leaving every contrast above small. No finite ladder of windows can
  exclude this.
* The statistical statements are large-sample statements about finite-window ratios; they say nothing about mu
  without section 4.
* The claim is about the observable full. The other four observables are controls.
