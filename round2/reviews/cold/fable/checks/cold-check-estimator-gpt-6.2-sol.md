VERDICT: GAP

The fixed-weight pipeline passed the independent generator and small-tree simulator checks and estimates a finite-window population ratio. Its jackknife describes sampling uncertainty about that ratio. The proposed friend-tree systematic bound below the sampling SE needs a separate, stated contraction assumption for the frozen full observable. The PAN result is within 1.83 SE of its known root. The DAN benchmark at (9,12) differs from its root by +6.07e-8, or 12.62 sampling SE, and shows why a good generator and a small SE do not establish accuracy for the asymptotic exponent.

## Scope

I used the lead binaries or captured copies for the unmodified code and built the fault-injected copy under this directory. The two validation layouts use the defaults except for `Rpow=3 Rsum=3 Dp=7 Dk=8` where specified. Source inspection and exact arithmetic checks distinguish floating-point roundoff from differences in the modeled transitions.

## Exact model check

I seeded Python's random generator with 5927831 and generated 296 random trees (111 uniform Prüfer trees, 111 caterpillars and 74 friend-chain trees), alongside 37 stars, 37 paths and two two-hub boundary cases. Sizes ran from 4 through 40. The path and caterpillar cases include internal degree-2 vertices with zero leaves. The implementation evaluated the named feature coordinates from scratch on full adjacency before and after each attachment, enumerated all ordered pairs `(V,W)` with `Fraction(1,n*deg(V))`, and pooled identical destinations before computing the exact weighted sum. Core membership, parent relation and newborn contributions were determined by the updated full tree. I compared all active entries of both the exact and printed vectors, including spurious nonzero entries, and checked the printed `n,C,S`.

Commands and output are below. `fail_a` and `fail_b` count discrepancies greater than `1e-12*(1+abs(exact))`. The largest absolute discrepancy was `3.333e-13` (Q, both layouts). The largest scaled discrepancy was `4.537e-15` (R, extended layout). This setup exercised transitions as well as values, including the degree caps and all enabled groups.

```text
$ taskset -c 14 python3 exact_tree_check.py
seed=5927831 shapes={'star': 37, 'path': 37, 'prufer': 111, 'caterpillar': 111, 'friend': 74, 'two-hubs': 1, 'two-mid': 1}
layout p=2262 trees=372 model=0 hub_decay=False literal_N2_B2=False
A  exercised=  63 max_abs_a=4.615e-14 max_abs_b=0.000e+00 max_scaled_a=3.471e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
H  exercised=  21 max_abs_a=3.333e-14 max_abs_b=3.077e-15 max_scaled_a=2.703e-15 max_scaled_b=1.053e-15 fail_a=0 fail_b=0
N  exercised= 430 max_abs_a=4.286e-14 max_abs_b=0.000e+00 max_scaled_a=3.472e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
B  exercised= 102 max_abs_a=4.615e-14 max_abs_b=4.545e-15 max_scaled_a=4.068e-15 max_scaled_b=1.852e-15 fail_a=0 fail_b=0
N2 exercised= 156 max_abs_a=3.333e-14 max_abs_b=4.444e-15 max_scaled_a=2.569e-15 max_scaled_b=2.132e-15 fail_a=0 fail_b=0
B2 exercised=  48 max_abs_a=5.238e-15 max_abs_b=4.444e-15 max_scaled_a=2.210e-15 max_scaled_b=1.818e-15 fail_a=0 fail_b=0
P  exercised= 186 max_abs_a=3.333e-14 max_abs_b=0.000e+00 max_scaled_a=3.008e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
Q  exercised= 156 max_abs_a=3.333e-13 max_abs_b=0.000e+00 max_scaled_a=3.846e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
R  exercised=   0 max_abs_a=0.000e+00 max_abs_b=0.000e+00 max_scaled_a=0.000e+00 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
layout p=4218 trees=372 model=0 hub_decay=False literal_N2_B2=False
A  exercised=  63 max_abs_a=4.615e-14 max_abs_b=0.000e+00 max_scaled_a=3.471e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
H  exercised=  21 max_abs_a=3.333e-14 max_abs_b=3.077e-15 max_scaled_a=2.703e-15 max_scaled_b=1.053e-15 fail_a=0 fail_b=0
N  exercised= 430 max_abs_a=4.286e-14 max_abs_b=0.000e+00 max_scaled_a=3.472e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
B  exercised= 102 max_abs_a=4.615e-14 max_abs_b=4.545e-15 max_scaled_a=4.068e-15 max_scaled_b=1.852e-15 fail_a=0 fail_b=0
N2 exercised= 156 max_abs_a=3.333e-14 max_abs_b=4.444e-15 max_scaled_a=2.569e-15 max_scaled_b=2.132e-15 fail_a=0 fail_b=0
B2 exercised=  48 max_abs_a=5.238e-15 max_abs_b=4.444e-15 max_scaled_a=2.210e-15 max_scaled_b=1.818e-15 fail_a=0 fail_b=0
P  exercised= 272 max_abs_a=3.333e-14 max_abs_b=0.000e+00 max_scaled_a=3.008e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
Q  exercised= 370 max_abs_a=3.333e-13 max_abs_b=0.000e+00 max_scaled_a=4.129e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
R  exercised= 445 max_abs_a=5.185e-14 max_abs_b=4.815e-15 max_scaled_a=4.537e-15 max_scaled_b=2.203e-15 fail_a=0 fail_b=0
layout p=2262 trees=2 model=0 hub_decay=False literal_N2_B2=True
A  exercised=   7 max_abs_a=3.333e-14 max_abs_b=0.000e+00 max_scaled_a=2.941e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
H  exercised=   0 max_abs_a=0.000e+00 max_abs_b=0.000e+00 max_scaled_a=0.000e+00 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
N  exercised=   8 max_abs_a=3.333e-15 max_abs_b=0.000e+00 max_scaled_a=1.250e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
B  exercised=   9 max_abs_a=1.893e-15 max_abs_b=0.000e+00 max_scaled_a=5.672e-16 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
N2 exercised=   7 max_abs_a=4.850e+00 max_abs_b=6.400e-01 max_scaled_a=4.850e+00 max_scaled_b=6.400e-01 fail_a=11 fail_b=3
B2 exercised=   6 max_abs_a=2.498e+00 max_abs_b=6.400e-01 max_scaled_a=2.498e+00 max_scaled_b=6.400e-01 fail_a=6 fail_b=3
P  exercised=   0 max_abs_a=0.000e+00 max_abs_b=0.000e+00 max_scaled_a=0.000e+00 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
Q  exercised=   0 max_abs_a=0.000e+00 max_abs_b=0.000e+00 max_scaled_a=0.000e+00 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
R  exercised=   0 max_abs_a=0.000e+00 max_abs_b=0.000e+00 max_scaled_a=0.000e+00 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
witness N2 a (case, kind, feature, exact, printed): (2, 'two-mid', ('N2', 2, 1, (6, 2, 0)), '0', '-485037037037037/100000000000000')
witness N2 b (case, kind, feature, exact, printed): (2, 'two-mid', ('N2', 2, 1, (6, 2, 0)), '0', '16/25')
witness B2 a (case, kind, feature, exact, printed): (1, 'two-hubs', ('B2', 2, 1, 0, 2), '0', '-3903420781893/1562500000000')
witness B2 b (case, kind, feature, exact, printed): (1, 'two-hubs', ('B2', 2, 1, 0, 0), '0', '16/25')
```

The top specification comment underdescribes N2/B2. It names the extra factor `x = L1cap/L` only. The detailed emitter comment defines both `x` and `x^2`, and the layout reserves both blocks. The first two runs above apply the documented two-block definition. The final run evaluates the opening comment literally and reports the extra N2/B2 values and drifts as discrepancies. For example, two degree-6 cores with five leaves each produce an N2 second-block value of `16/25`. The additional coordinates are internally valid observables with the correct drift. This is a documentation error.

## Independent simulator

The discrete and continuous-time implementations share the explicit-adjacency attachment rule and use the standard library MT19937-64. The CT simulation uses total rate `n`, independent exponential inter-event times, and retains overshoots between checkpoints. It sums the inverse degree of each leaf's sole neighbour for S. The first test takes exact discrete means from column 3 of the supplied file. At n=3 both outputs are identically 1 and there is zero sampling variance. The CT comparison combines the two independent standard errors in quadrature. To recover the beta-check SE for C from its format, I used the same reference runs with `mu0=0,+1,-1` and the identity `Var(C) = (Var(S-C)+Var(S+C))/2 - Var(S)`. The beta output rounds the SEs to six decimals. Checkpoints within a group share trajectories.

```text
$ g++ -O2 -std=c++17 -Wall -Wextra -o naive_sim naive_sim.cpp
$ ./naive_sim discrete 99628376051 1000000
# mode=discrete seed=99628376051 runs=1000000
3 1.000000000 0.000000000
4 1.333353000 0.000471412
5 1.583462000 0.000640146
6 1.797191000 0.000767286
7 1.989570000 0.000875961
8 2.166846000 0.000973354
9 2.333011000 0.001063959
10 2.489787000 0.001148098
11 2.639143000 0.001227608
12 2.782273000 0.001302246
$ ./naive_sim ct 99628376052 500000
# mode=ct seed=99628376052 runs=500000
1 2.117882000 0.001658281 1.335515554 0.000566718
2 3.864486000 0.003273414 2.261925130 0.001475963
3 6.904680000 0.006009851 3.980770547 0.003029811
4 12.260610000 0.010799020 7.021679378 0.005754408
$ ./baseline/beta_check 99628376053 500000 4 0
# seed=99628376053 runs=500000 mu0=0.000000000
# t meanC meanS mean(S-mu0*C) se
1 2.116770 1.335248 1.335248 0.000566
2 3.863052 2.261892 2.261892 0.001476
3 6.901874 3.979909 3.979909 0.003030
4 12.252562 7.018508 7.018508 0.005749
$ ./baseline/beta_check 99628376053 500000 4 1
# seed=99628376053 runs=500000 mu0=1.000000000
# t meanC meanS mean(S-mu0*C) se
1 2.116770 1.335248 -0.781522 0.001219
2 3.863052 2.261892 -1.601160 0.001986
3 6.901874 3.979909 -2.921965 0.003218
4 12.252562 7.018508 -5.234054 0.005323
$ ./baseline/beta_check 99628376053 500000 4 -1
# seed=99628376053 runs=500000 mu0=-1.000000000
# t meanC meanS mean(S-mu0*C) se
1 2.116770 1.335248 3.452018 0.002152
2 3.863052 2.261892 6.124944 0.004666
3 6.901874 3.979909 10.881783 0.008954
4 12.252562 7.018508 19.271070 0.016456
$ python3 simulator_compare.py
Discrete: n, explicit mean, exact mean, explicit SE, z
 3 1.000000000 1.000000000 0.000000000 +0.000
 4 1.333353000 1.333333333 0.000471412 +0.042
 5 1.583462000 1.583333333 0.000640146 +0.201
 6 1.797191000 1.797222222 0.000767286 -0.041
 7 1.989570000 1.989467593 0.000875961 +0.117
 8 2.166846000 2.166658859 0.000973354 +0.192
 9 2.333011000 2.332470537 0.001063959 +0.508
10 2.489787000 2.489218217 0.001148098 +0.495
11 2.639143000 2.638484346 0.001227608 +0.537
12 2.782273000 2.781414022 0.001302246 +0.660
CT: t, metric, explicit mean, beta_check mean, explicit SE, beta_check SE, z
 1 C 2.117882000 2.116770000 0.001658281 0.001654744 +0.475
 1 S 1.335515554 1.335248000 0.000566718 0.000566000 +0.334
 2 C 3.864486000 3.863052000 0.003273414 0.003267920 +0.310
 2 S 2.261925130 2.261892000 0.001475963 0.001476000 +0.016
 3 C 6.904680000 6.901874000 0.006009851 0.006006989 +0.330
 3 S 3.980770547 3.979909000 0.003029811 0.003030000 +0.201
 4 C 12.260610000 12.252562000 0.010799020 0.010794264 +0.527
 4 S 7.021679378 7.018508000 0.005754408 0.005749000 +0.390
```

## Fixed-weight statistics and stream independence

For one fixed training outcome let `X_r = num_r(t2)-num_r(t1)` and `Y_r = den_r(t2)-den_r(t1)`. The production point estimator is `sum(X)/sum(Y)`. Its population target for the chosen window is `R=E[X]/E[Y]`. The estimating equation is unbiased at R, while its ratio solution generally has first-order bias

```text
[R Var(Y) - Cov(X,Y)] / [n E[Y]^2].
```

Independence of training and production does not remove this ratio bias or the finite-time difference `R-mu`. The delete-one-run jackknife in `prod_eval.py` is appropriate for the *conditional* sampling SE if runs are independent, second moments are finite and the mean denominator is well separated from zero. It preserves checkpoint pairing. The reported point estimate has no jackknife bias correction. For comparison, `window_check.py` calculates an equivalent residual-form jackknife SE and estimates first-order ratio bias from the leave-one-out replicates. It also checks denominator leverage and unique row IDs within each shard.

`fit_theta.py` reads only its supplied training prefixes when it writes the all-training theta. The production tool reads the theta file and does no fitting. The fit uses the raw second moment of `da-mu0*db` with a diagonal ridge, which is equivalent to using the population covariance when the residual mean vanishes. The training cross-fit uses a common feature screen selected across all training folds. The leave-one-out diagnostic there does not refit weights, so it does not inherit the fixed-weight production jackknife justification. Show the independent production SE. Selecting weights or windows using a pilot and then running fresh production preserves conditional inference for that chosen finite-window target. Selection among estimates on the final production draws would require separate confirmation or accounting for selection.

The binaries key runs by `(seed,run)`, independently of the process ID. Repeated pairs repeat the initialized stream. Within one seed, the odd multiplier makes `run -> preseed` injective over the uint64 range. Across seeds there are exact collisions, and one can construct them by inverting SplitMix64. The xoshiro streams use no explicit jump-separated substreams. Distinct initial states are not a proof of disjoint consumed segments or statistical independence. I found zero initial-state collisions in the lead's stated training and ladder manifest plus the four independent validation shards. The deliberate matched control/fault simulations reuse a stream.

```text
$ taskset -c 15 python3 seed_check.py 601:0:20000 602:0:20000 7000:0:2000000 7001:0:300000 7002:0:30000 7003:0:3000 62863482315701:0:30000 62863482315702:0:30000 62863482315703:0:30000 62863482315704:0:30000
shards=10 runs=2493000 distinct pre-seeds=2493000 collisions=0
constructed collision: (601, 0) and (1306363514317730328, 1) map to 0x248b6a394cfb5000
identical state: 03923b727a8d68b6 ac546219739f3e26 b125d179cab8754a d1c586cd8e9136e9
```

`prod_eval.py` checks checkpoint row alignment but sees run numbers without seed provenance. It does not detect submitting a shard twice, overlapping `(seed,run)` ranges, or choosing the wrong weights with a matching byte length. The seed and run manifest is the independence check. Retain it for final production with theta and binary hashes.

## Fault injection

I changed only the PAN conversion rate in the copied feature generator from 1 to 0.99 for core vertices with leaves. The copied simulator was left unchanged. Both control and fault used the PAN training seed and first 3,000 runs, the same `12,15` checkpoint schedule, and the extended layout. The genuine plain-moment z-score was +0.89 and the faulty one was -40.37. A friend-only perturbation of the `L/D` branch would leave both PAN and DAN untouched because the model rates use disjoint branches. The exact friend-tree enumeration supplies the direct check of that branch.

```diff
--- a/src/feat2.hpp
+++ b/src/feat2.hpp
@@ -152,7 +152,7 @@
 static inline void rates(const Snap &T, uint32_t v, int model, double &gain, double &conv) {
     double L = (double)T.L[v];
     if (model == 0) { gain = L + T.cred[v]; conv = L / (L + (double)T.c[v]); }
-    else { gain = L + ((model == 2 && T.L[v] == 0) ? 1.0 : 0.0); conv = T.L[v] > 0 ? 1.0 : 0.0; }
+    else { gain = L + ((model == 2 && T.L[v] == 0) ? 1.0 : 0.0); conv = T.L[v] > 0 ? (model == 1 ? 0.99 : 1.0) : 0.0; }
 }
 // ------------------------------------------------------------------ brute force (the specification)
 template <class S> static void feat_bf(const Layout &Y, const Snap &T, uint32_t u, double w, S &s) {
```

```text
+ ./fault/build.sh
built: beta_check feat_tree ftprod ftvec layout_info test_feat2 test_yule validate_small 
+ ./fault/bin/ftvec 62863482315701 3000 pan_fault 12 15 model=1 Rpow=3 Rsum=3 Dp=7 Dk=8 diffout=1
+ ./baseline/ftvec 62863482315701 3000 pan_control 12 15 model=1 Rpow=3 Rsum=3 Dp=7 Dk=8 diffout=1
+ ../../.venv/bin/python ../../src/eq_ztest.py pan_control 0.549734972504018815
runs=3000 equations tested=2751 at mu=0.549734972504
z-scores: mean=+0.014 sd=1.010 max|z|=3.32  fraction |z|>2: 0.050 (normal: 0.046)  |z|>3: 0.0007 (normal: 0.0027)
plain equation z = +0.89
split-sample Hansen J = 140115746.5 for 2751 equations (expect about 2751 +- 74 if all equations hold at this mu)
+ ../../.venv/bin/python ../../src/eq_ztest.py pan_fault 0.549734972504018815
runs=3000 equations tested=2751 at mu=0.549734972504
z-scores: mean=-0.209 sd=1.513 max|z|=40.37  fraction |z|>2: 0.079 (normal: 0.046)  |z|>3: 0.0185 (normal: 0.0027)
plain equation z = -40.37
split-sample Hansen J = 239562047.0 for 2751 equations (expect about 2751 +- 74 if all equations hold at this mu)
$ paired memmap check (same-seed difference records)
dC bitwise equal: True
db bitwise equal: True
max |dS_fault-0.99*dS_control|: 1.1386873666197062e-08
```

The largest marginal z-score and the plain equation detect the fault. The split-sample J line itself is not a calibrated discriminator at this sample size. The clean run reports `J=140115746.5` against its printed nominal `2751 ± 74`. The covariance fit has 1,500 rows for 2,751 columns, exact dependencies such as the appended plain count versus the sum of own-state counts, and rare coordinates. A diagonal ridge makes the matrix invertible but does not restore the declared chi-square null. This is a separate diagnostic defect.

## Finite-time population target

Let `m_theta(t)=E[den_theta(t)]`, `a_theta(t)=E[num_theta(t)]=m'_theta(t)`, and `beta_theta(t)=a_theta(t)-mu*m_theta(t)`. The generator identity holds for any fixed finite weight vector on the ideal embedded chain. Each feature total and its pointwise drift have polynomial bounds in the current population, whose finite-time Yule law has all polynomial moments. The resulting population identity holds at any finite window with nonzero denominator.

```text
R_theta(t1,t2) - mu
  = [beta_theta(t2) - beta_theta(t1)] / [m_theta(t2) - m_theta(t1)].
```

Differencing removes a constant part of `beta_theta`. The plain `E[S]-mu E[C]` is the residual for a different observable. Its value near 0.09 and its uncertainty give no quantitative bound on the fitted residual drift. A residual computed at a reference `mu0` also shifts the displayed window ratio by `mu-mu0`.

The generator identity alone does not yield consistency for arbitrary theta. There is an exact counterexample in this feature map. Give unit weight to every A coordinate with c=0 and to the constant H term with c=0. Give all other coordinates zero weight. The core induces a connected tree, so the resulting observable is precisely `1[C=1]`. The sole core's rate-1 ring always converts a leaf, while other events keep `C=1`. Once `C>1` it cannot go back, so `m_theta(t)=exp(-t)` and both ratios equal -1 at every time/window. The same observable can be represented with a plain weight of 1 by subtracting the sum of all own-state count features and then adding the c=0 indicator. I shared the example with the lead, who agreed to restrict the statement to the selected observable.

The repository proves logarithmic growth `E[C(t)]=exp(mu*t+o(t))`. A sufficient asymptotic condition for the fitted observable at fixed positive window width is `m_theta(t)/exp(mu*t) -> a != 0` and `a_theta(t)/exp(mu*t) -> mu*a`. A stable, nonzero `E[den_theta]/E[C]` is evidence for the leading projection, subject to the growth and drift of `E[C]`. It does not on its own control the derivative or the absolute error in the fitted ratio.

For a numerical claim below `1e-8` I would require the fixed final theta evaluated on independent matched draws across several later windows, with ratio estimates, SEs, their covariance when paired, `E[den]/E[C]`, and the projected residual difference. Vary starting times and widths, evaluate the plain and a separately constructed late observable as controls, and reserve a late confirmatory window from all window and feature selection. The later-window comparison needs precision well below the allowance for finite-time bias. To turn plateau evidence into an absolute bound also requires a justified constraint on remaining modes or on the projected residual. A slow correction can evade close-window contrasts. For example, if a bias is `1e-6` and decays at rate `0.001`, moving both endpoints four time units changes it by only about `4e-9`. A defensible numerical extrapolation should propagate uncertainty in the decay law, test additional or slower modes against withheld data, and state its assumptions. Without that constraint, the ladder measures window-to-window stability and a claimed absolute bound is conditional on an extrapolation model.

## Additional layout and arithmetic checks

The same independent enumeration also passed for the extended layout with `topx=1 Bconst=0`. That check suppresses the plain N/B top-bin entries and the constant B/B2 high-degree-neighbour terms when producing the exact reference. Its maximum discrepancies were `3.333e-13` absolute and `4.537e-15` scaled.

```text
$ taskset -c 14 python3 exact_tree_check.py --hub-decay
seed=5927831 shapes={'star': 37, 'path': 37, 'prufer': 111, 'caterpillar': 111, 'friend': 74, 'two-hubs': 1, 'two-mid': 1}
layout p=4218 trees=372 model=0 hub_decay=True literal_N2_B2=False
A  exercised=  63 max_abs_a=4.615e-14 max_abs_b=0.000e+00 max_scaled_a=3.471e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
H  exercised=  21 max_abs_a=3.333e-14 max_abs_b=3.077e-15 max_scaled_a=2.703e-15 max_scaled_b=1.053e-15 fail_a=0 fail_b=0
N  exercised= 352 max_abs_a=4.286e-14 max_abs_b=0.000e+00 max_scaled_a=3.472e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
B  exercised=  52 max_abs_a=3.333e-14 max_abs_b=4.201e-15 max_scaled_a=2.703e-15 max_scaled_b=1.811e-15 fail_a=0 fail_b=0
N2 exercised= 156 max_abs_a=3.333e-14 max_abs_b=4.444e-15 max_scaled_a=2.569e-15 max_scaled_b=2.132e-15 fail_a=0 fail_b=0
B2 exercised=  32 max_abs_a=5.238e-15 max_abs_b=3.882e-15 max_scaled_a=1.947e-15 max_scaled_b=1.813e-15 fail_a=0 fail_b=0
P  exercised= 272 max_abs_a=3.333e-14 max_abs_b=0.000e+00 max_scaled_a=3.008e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
Q  exercised= 370 max_abs_a=3.333e-13 max_abs_b=0.000e+00 max_scaled_a=4.129e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
R  exercised= 445 max_abs_a=5.185e-14 max_abs_b=4.815e-15 max_scaled_a=4.537e-15 max_scaled_b=2.203e-15 fail_a=0 fail_b=0
```

I copied the lead's new `RFast` production source and added counters for the existing large-lambda Poisson approximation. Using identical trajectories and theta files, I compared the accelerated cubic-polynomial evaluation with the generic feature-emission path. The aggregate ratios differed by at most `2.467e-13` on this sample. Both runs reached the simulator's `lambda > 1e9` normal approximation four times. The approximation changes the ideal Yule sampling law on those trajectories. These runs establish that the branch is live at t=20 and leave its impact on mu unquantified. They do not measure the impact for a final production schedule.

```text
+ g++ -O2 -std=c++17 -o opt-check/ftprod opt-check/ftprod.cpp
+ ./opt-check/ftprod 7003 0 192 opt-check/late-fast ../../outputs/train/thetaX.bin,../../outputs/train/thetaY.bin 16 20 Rsum=3 Rpow=3 Dp=7 Dk=8 rfast=1
Poisson normal calls=4 max lambda=1390249014.02
+ ./opt-check/ftprod 7003 0 192 opt-check/late-generic ../../outputs/train/thetaX.bin,../../outputs/train/thetaY.bin 16 20 Rsum=3 Rpow=3 Dp=7 Dk=8 rfast=0
Poisson normal calls=4 max lambda=1390249014.02
$ ../../.venv/bin/python compare_rfast.py
same metadata, C, S for 384 run-checkpoint rows
theta1 largest abs num delta=1.0012e-08 den delta=1e-07 ratio fast-generic=+1.1102e-14 max abs ratio-scale per-run delta=6.8716e-13
theta2 largest abs num delta=1e-07 den delta=1e-07 ratio fast-generic=-2.4669e-13 max abs ratio-scale per-run delta=2.3243e-12
```

## Completed pilot ladder

The lead supplied fixed weights thetaX and thetaY for the initial ladder. The table reproduces the differenced production evaluator rows. `ladder_final.log` records all checkpoints and denominator/residual diagnostics. `ladder_04_08.log` has the first-rung diagnostics.

```text
$ for file in ../../outputs/ladder1/r*.txt; do ../../.venv/bin/python ../../src/prod_eval.py "$file"; done
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/ladder1/r04_08.txt
runs=2000000 checkpoints=[np.float64(4.0), np.float64(8.0)]  mean C at checkpoints: [np.float64(12.3), np.float64(119.2)]
window   4.0 -> 8     | plain: 0.565753866 +- 2.01e-05 | th1: 0.565878774 +- 6.29e-07 | th2: 0.565650031 +- 2.85e-07
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/ladder1/r08_12.txt
runs=300000 checkpoints=[np.float64(8.0), np.float64(12.0)]  mean C at checkpoints: [np.float64(119.2), np.float64(1147.3)]
window   8.0 -> 12    | plain: 0.565705279 +- 1.68e-05 | th1: 0.565725742 +- 6.74e-08 | th2: 0.565723033 +- 6.82e-08
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/ladder1/r12_16.txt
runs=30000 checkpoints=[np.float64(12.0), np.float64(16.0)]  mean C at checkpoints: [np.float64(1147.2), np.float64(11025.9)]
window  12.0 -> 16    | plain: 0.565696275 +- 1.70e-05 | th1: 0.565723679 +- 5.57e-08 | th2: 0.565723656 +- 6.88e-08
$ ../../.venv/bin/python ../../src/prod_eval.py ../../outputs/ladder1/r16_20.txt
runs=3000 checkpoints=[np.float64(16.0), np.float64(20.0)]  mean C at checkpoints: [np.float64(10741.8), np.float64(103235.8)]
window  16.0 -> 20    | plain: 0.565719944 +- 1.82e-05 | th1: 0.565723789 +- 5.67e-08 | th2: 0.565723836 +- 7.18e-08
```

For thetaY, moving from (12,16) to (16,20) gives a difference of `+1.803e-7 ± 9.945e-8` using the independent-rung SE. Its denominator/core ratio is about 0.7836 at the late checkpoints and its jackknife sampling bias estimates are order `1e-12`. Those diagnostics support a well-conditioned ratio. At the reference `0.56572374`, the (12,16) projected residual difference is `-0.00065023 ± 0.00053284` with mean denominator difference 7740.65. A `1e-8` systematic allowance there would correspond to a residual difference below about `0.0000774` in absolute value. The (16,20) window gives `+0.00697781 ± 0.00520221` with mean denominator difference 72477.63. Its corresponding allowance would be about `0.000725`. The pilot ladder does not resolve the requested absolute bias.

## Sensitivity to the fitting reference

I fitted three versions of the original extended layout on the seed-601 friend-tree shard with `mu0` displaced by `±1e-5` and evaluated the fitted weights on the held-out seed-602 shard. This isolates the training-reference choice with paired evaluation. The observed displacement in the resulting ratio was about 2% of the displacement in the fitting reference, with uncertainty comparable to the displacement in the ratio. This check concerns the older layout and window (12,15). The fit penalizes `E[(da-mu0*db)^2]`, which contains a squared population mean residual, so train/production separation by itself does not establish immunity to the choice of `mu0` when the finite-time moments differ.

```text
+ for m in 0.56571374 0.56572374 0.56573374
+ ../../.venv/bin/python ../../src/fit_theta.py ../../outputs/train/X_a target-check/weights_0.56571374.bin --mu0 0.56571374
training runs=20000 layout p=4218 equations kept q=2677 (minfrac=0.01)
plain (differenced)      : 0.565737132 +- 3.07e-05
cross-fitted GMM (5 folds): 0.565723940 +- 6.90e-07   se ratio vs plain = 44.48
theta (all runs) written to target-check/weights_0.56571374.bin
+ for m in 0.56571374 0.56572374 0.56573374
+ ../../.venv/bin/python ../../src/fit_theta.py ../../outputs/train/X_a target-check/weights_0.56572374.bin --mu0 0.56572374
training runs=20000 layout p=4218 equations kept q=2677 (minfrac=0.01)
plain (differenced)      : 0.565737132 +- 3.07e-05
cross-fitted GMM (5 folds): 0.565723922 +- 1.15e-07   se ratio vs plain = 265.82
theta (all runs) written to target-check/weights_0.56572374.bin
+ for m in 0.56571374 0.56572374 0.56573374
+ ../../.venv/bin/python ../../src/fit_theta.py ../../outputs/train/X_a target-check/weights_0.56573374.bin --mu0 0.56573374
training runs=20000 layout p=4218 equations kept q=2677 (minfrac=0.01)
plain (differenced)      : 0.565737132 +- 3.07e-05
cross-fitted GMM (5 folds): 0.565723926 +- 6.78e-07   se ratio vs plain = 45.24
theta (all runs) written to target-check/weights_0.56573374.bin
$ ../../.venv/bin/python -B target_check.py
target-check/weights_0.56571374.bin: heldout est=0.5657232927919 se=2.28049e-07 meanD=4695.0887
target-check/weights_0.56572374.bin: heldout est=0.5657234970316 se=1.09919e-07 meanD=4693.2347
target-check/weights_0.56573374.bin: heldout est=0.5657237012697 se=2.02972e-07 meanD=4691.3836
paired estimate 0 minus middle: -2.0423967e-07 SE=1.85794e-07
paired estimate 2 minus middle: +2.0423815e-07 SE=1.85792e-07
```

## Independent exact-model pipeline

The supplied local text for Hartle–Krapivsky–Redner–Zhang, arXiv:2604.01540, equations (17) and (22), gives the DAN and PAN roots. The code's two model rates match the described rules: leaf clocks increase the corresponding leaf count, a core with leaves converts at rate one and a core without leaves rejects its ringing attempt (PAN) or gains a leaf (DAN). These validation rates lack friend-tree core-to-core redirection.

The validation jobs use separate seeds for training and production, 30,000 runs for each. The layout is the extended layout above. Fits supply the known exponent as the pilot parameter, as requested.

I evaluated the two roots independently by integrating the power series term by term, using the coefficients `exp(-1)/k!` and positive denominators `1+mu+k` in Decimal arithmetic.

```text
$ python3 exact_model_roots.py
PAN 0.549734972504018815207626667733775464614903108536879799776695
DAN 0.771192302407000256316924280951926934487082438054978762769935
```

### PAN result

The independent fixed-weight result is consistent with the exact root at about two sampling standard errors. The scientific result uses the distinct production seed, as opposed to the training cross-fit line displayed below.

```text
+ ../../bin/ftvec 62863482315701 30000 pan_train 12 15 model=1 Rpow=3 Rsum=3 Dp=7 Dk=8 diffout=1
+ ../../.venv/bin/python ../../src/fit_theta.py pan_train pan.theta.bin --mu0 0.5497349725
training runs=30000 layout p=4218 equations kept q=2751 (minfrac=0.01)
plain (differenced)      : 0.549692616 +- 3.76e-05
cross-fitted GMM (5 folds): 0.549734971 +- 9.76e-09   se ratio vs plain = 3851.63
theta (all runs) written to pan.theta.bin
+ ../../bin/ftprod 62863482315702 0 30000 pan_prod pan.theta.bin 12 15 model=1 Rpow=3 Rsum=3 Dp=7 Dk=8
+ ../../.venv/bin/python ../../src/prod_eval.py pan_prod.txt
runs=30000 checkpoints=[np.float64(12.0), np.float64(15.0)]  mean C at checkpoints: [np.float64(1469.0), np.float64(7645.9)]
window  None -> 12    | plain: 0.550096854 +- 7.68e-05 | th1: 0.549978283 +- 1.08e-06
window  None -> 15    | plain: 0.549804234 +- 3.37e-05 | th1: 0.549781742 +- 2.07e-07
window  12.0 -> 15    | plain: 0.549734642 +- 3.73e-05 | th1: 0.549734990 +- 9.58e-09
+ ../../.venv/bin/python ../../src/eq_ztest.py pan_train 0.549734972504018815
runs=30000 equations tested=2751 at mu=0.549734972504
z-scores: mean=-0.015 sd=0.976 max|z|=3.23  fraction |z|>2: 0.035 (normal: 0.046)  |z|>3: 0.0015 (normal: 0.0027)
plain equation z = -1.13
split-sample Hansen J = 2970.0 for 2751 equations (expect about 2751 +- 74 if all equations hold at this mu)
$ ../../.venv/bin/python window_check.py 0.549734972504018815 pan_prod.txt:30000
pan_prod.txt runs=30000 window=12:15 meanC=(1469.0091,7645.9038) mu_ref=0.549734972504
  plain est=0.5497346420175 se=3.7308e-05 delta=-3.3048651e-07 z=-0.0089 JKratioBias=+8.5103e-10 meanD=6176.8948 nonposD=0 maxAbsDshare=0.000192675
        den/C=(1,1) beta1=+0.53160779 se1=0.112848 beta2=+0.52956641 se2=0.257428 deltaBeta=-0.0020413804 seDelta=0.230444
  th1   est=0.5497349900367 se=9.57878e-09 delta=+1.7532729e-08 z=+1.8304 JKratioBias=+4.6265e-13 meanD=4050.1668 nonposD=0 maxAbsDshare=0.000193257
        den/C=(0.655832856,0.655722503) beta1=+0.23441089 se1=0.000478934 beta2=+0.23448191 se2=0.000479822 deltaBeta=+7.1010478e-05 seDelta=3.87927e-05
```

### DAN result

The independent production estimate is **12.62 sampling SE above the exact DAN root** at the specified (9,12) window. The fitted projected residual grows by `0.00064798136 ± 0.0000513441`. Dividing by the mean denominator difference `10668.067` gives the observed offset `+6.074028e-8`. The estimated ratio bias is `6.75e-13`, too small to explain the discrepancy. This is an empirical counterexample to treating a precise differenced ratio as an asymptotic exponent without a finite-time check.

```text
+ ../../bin/ftvec 62863482315703 30000 dan_train 9 12 model=2 Rpow=3 Rsum=3 Dp=7 Dk=8 diffout=1
+ ../../.venv/bin/python ../../src/fit_theta.py dan_train dan.theta.bin --mu0 0.771192302407000256
training runs=30000 layout p=4218 equations kept q=2792 (minfrac=0.01)
plain (differenced)      : 0.771157558 +- 2.39e-05
cross-fitted GMM (5 folds): 0.771192351 +- 5.12e-09   se ratio vs plain = 4660.94
theta (all runs) written to dan.theta.bin
+ ../../bin/ftprod 62863482315704 0 30000 dan_prod dan.theta.bin 9 12 model=2 Rpow=3 Rsum=3 Dp=7 Dk=8
+ ../../.venv/bin/python ../../src/prod_eval.py dan_prod.txt
runs=30000 checkpoints=[np.float64(9.0), np.float64(12.0)]  mean C at checkpoints: [np.float64(1247.7), np.float64(12613.9)]
window  None -> 9     | plain: 0.771189938 +- 6.65e-05 | th1: 0.771308069 +- 5.01e-07
window  None -> 12    | plain: 0.771172773 +- 2.07e-05 | th1: 0.771203806 +- 4.97e-08
window   9.0 -> 12    | plain: 0.771170889 +- 2.39e-05 | th1: 0.771192363 +- 4.81e-09
+ ../../.venv/bin/python ../../src/eq_ztest.py dan_train 0.771192302407000256
runs=30000 equations tested=2792 at mu=0.771192302407
z-scores: mean=-0.002 sd=1.016 max|z|=4.17  fraction |z|>2: 0.047 (normal: 0.046)  |z|>3: 0.0043 (normal: 0.0027)
plain equation z = -1.46
split-sample Hansen J = 4941.0 for 2792 equations (expect about 2792 +- 75 if all equations hold at this mu)
$ ../../.venv/bin/python window_check.py 0.771192302407000256 dan_prod.txt:30000
dan_prod.txt runs=30000 window=9:12 meanC=(1247.659,12613.947) mu_ref=0.771192302407
  plain est=0.7711708886910 se=2.38948e-05 delta=-2.1413716e-05 z=-0.8962 JKratioBias=+9.7481e-11 meanD=11366.288 nonposD=0 maxAbsDshare=0.000218776
        den/C=(1,1) beta1=-0.0029501621 se1=0.0829134 beta2=-0.24634462 se2=0.261448 deltaBeta=-0.24339446 seDelta=0.271593
  th1   est=0.7711923631473 se=4.81438e-09 delta=+6.074028e-08 z=+12.6164 JKratioBias=+6.7456e-13 meanD=10668.067 nonposD=0 maxAbsDshare=0.000219204
        den/C=(0.938448078,0.938558721) beta1=+0.13554736 se1=4.48309e-05 beta2=+0.13619534 se2=4.28582e-05 deltaBeta=+0.00064798136 seDelta=5.13441e-05
```

I also enumerated PAN and DAN one-step transitions independently from the paper's attachment rules on the same 372 trees, including PAN no-ops. The largest absolute drift discrepancy in each model is `4.815e-14`. There were zero out-of-tolerance entries, as recorded group by group in `exact_tree_pan.log` and `exact_tree_dan.log`. These tests concern the generator applied to the explicit small-tree state. They do not explain the time profile by themselves.

```text
$ python3 exact_tree_check.py --model 1 > exact_tree_pan.log
seed=5927831 shapes={'star': 37, 'path': 37, 'prufer': 111, 'caterpillar': 111, 'friend': 74, 'two-hubs': 1, 'two-mid': 1}
layout p=4218 trees=372 model=1 hub_decay=False literal_N2_B2=False
R  exercised= 445 max_abs_a=4.815e-14 max_abs_b=4.815e-15 max_scaled_a=4.368e-15 max_scaled_b=2.203e-15 fail_a=0 fail_b=0
$ python3 exact_tree_check.py --model 2 > exact_tree_dan.log
seed=5927831 shapes={'star': 37, 'path': 37, 'prufer': 111, 'caterpillar': 111, 'friend': 74, 'two-hubs': 1, 'two-mid': 1}
layout p=4218 trees=372 model=2 hub_decay=False literal_N2_B2=False
R  exercised= 445 max_abs_a=4.815e-14 max_abs_b=4.815e-15 max_scaled_a=4.150e-15 max_scaled_b=2.203e-15 fail_a=0 fail_b=0
```

The validation training processes used the old binary before the large-mean sampler and R-group optimizer update. Both production process executables had the same newer binary hash, copied to `dan-followup/ftprod` for the later windows. The earlier Poisson approximation is replaced in the new code by splitting the mean into parts at most `1e8` and summing independent draws. The mathematical Poisson sum law is exact. The implementation still uses floating-point rejection tests.

```text
$ sha256sum /proc/<active-validation-pid>/exe
seed=62863482315702 process executable c3396f4cfb9fc0f7abbca7ad91bdb04fc862918684c19bd9a8f57feb6e9ba3fb  /proc/559407/exe
seed=62863482315704 process executable c3396f4cfb9fc0f7abbca7ad91bdb04fc862918684c19bd9a8f57feb6e9ba3fb  /proc/566301/exe
```

### DAN time ladder with fixed weights

The first follow-up independently generated 6,000 trajectories at checkpoints 9, 12 and 14 using the existing DAN theta. Both estimates are relative to the externally known root. The paired data show a drop in the bias across these windows, despite the narrower second window. The sample errors are small enough to expose the early bias. The trajectory count, time shifts and widths matter when describing any observed decay.

```text
+ ./dan-followup/ftprod 62863482315705 0 6000 dan-followup/run dan.theta.bin 9 12 14 model=2 Rpow=3 Rsum=3 Dp=7 Dk=8 rfast=1
+ ../../.venv/bin/python ../../src/prod_eval.py dan-followup/run.txt
runs=6000 checkpoints=[np.float64(9.0), np.float64(12.0), np.float64(14.0)]  mean C at checkpoints: [np.float64(1252.1), np.float64(12660.4), np.float64(59185.8)]
window  None -> 9     | plain: 0.771232392 +- 1.46e-04 | th1: 0.771307665 +- 1.12e-06
window  None -> 12    | plain: 0.771191752 +- 4.62e-05 | th1: 0.771203750 +- 1.11e-07
window  None -> 14    | plain: 0.771211943 +- 2.13e-05 | th1: 0.771194758 +- 2.39e-08
window   9.0 -> 12    | plain: 0.771187291 +- 5.27e-05 | th1: 0.771192346 +- 1.09e-08
window   9.0 -> 14    | plain: 0.771211501 +- 2.19e-05 | th1: 0.771192318 +- 3.46e-09
window  12.0 -> 14    | plain: 0.771217438 +- 2.84e-05 | th1: 0.771192311 +- 4.24e-09
$ ../../.venv/bin/python -B dan-followup/compare.py
theta window 9:12 est=0.7711923462535 se=1.0895e-08 delta=+4.3846491e-08 z=+4.0245 JKratioBias=+1.3655e-12 meanD=10707.341 deltaBeta=+0.00046947932 deltaBetaSE=0.00011668
theta window 12:14 est=0.7711923106145 se=4.2404e-09 delta=+8.207499e-09 z=+1.9356 JKratioBias=-3.2016e-13 meanD=43668.248 deltaBeta=+0.0003584071 deltaBetaSE=0.000185212
theta window 9:14 est=0.7711923176323 se=3.45992e-09 delta=+1.5225333e-08 z=+4.4005 JKratioBias=+6.9357e-15 meanD=54375.589 deltaBeta=+0.00082788642 deltaBetaSE=0.000188289
paired late minus early = -3.5638992e-08 se=1.27829e-08 z=-2.7880
```

The extension used the first 2,400 trajectories from that follow-up and added checkpoint 16. The first three checkpoints match bit for bit, so the three-window comparisons share trajectories. The (14,16) ratio is only 0.10 SE below the exact root. The data resolve a drop from (9,12) to (12,14), with too little precision to determine the second drop or fit a general decay law. Both the initial displacement and the late agreement use an independent production seed and the same frozen DAN theta.

```text
+ ./dan-followup/ftprod 62863482315705 0 2400 dan-later/run dan.theta.bin 9 12 14 16 model=2 Rpow=3 Rsum=3 Dp=7 Dk=8 rfast=1
+ ../../.venv/bin/python ../../src/prod_eval.py dan-later/run.txt
runs=2400 checkpoints=[np.float64(9.0), np.float64(12.0), np.float64(14.0), np.float64(16.0)]  mean C at checkpoints: [np.float64(1236.2), np.float64(12507.1), np.float64(58477.3), np.float64(273407.9)]
window  None -> 9     | plain: 0.771291851 +- 2.32e-04 | th1: 0.771309033 +- 1.76e-06
window  None -> 12    | plain: 0.771230180 +- 7.39e-05 | th1: 0.771203896 +- 1.74e-07
window  None -> 14    | plain: 0.771237978 +- 3.31e-05 | th1: 0.771194788 +- 3.76e-08
window  None -> 16    | plain: 0.771203601 +- 1.54e-05 | th1: 0.771192834 +- 8.34e-09
window   9.0 -> 12    | plain: 0.771223416 +- 8.45e-05 | th1: 0.771192366 +- 1.74e-08
window   9.0 -> 14    | plain: 0.771236815 +- 3.40e-05 | th1: 0.771192321 +- 5.47e-09
window   9.0 -> 16    | plain: 0.771203200 +- 1.56e-05 | th1: 0.771192306 +- 2.21e-09
window  12.0 -> 14    | plain: 0.771240100 +- 4.45e-05 | th1: 0.771192310 +- 6.70e-09
window  12.0 -> 16    | plain: 0.771202327 +- 1.64e-05 | th1: 0.771192303 +- 2.32e-09
window  14.0 -> 16    | plain: 0.771194248 +- 2.06e-05 | th1: 0.771192302 +- 2.99e-09
$ ../../.venv/bin/python -B dan-later/compare.py
first three checkpoints bitwise equal to the earlier run for all 2400 trajectories
theta window 9:12 est=0.7711923662958 se=1.74358e-08 delta=+6.3888814e-08 z=+3.6642 JKratioBias=-2.3438e-12 meanD=10578.687 deltaBeta=+0.00067585977 deltaBetaSE=0.000184793
theta window 12:14 est=0.7711923099197 se=6.70078e-09 delta=+7.5127459e-09 z=+1.1212 JKratioBias=+3.391e-12 meanD=43147.592 deltaBeta=+0.00032415689 deltaBetaSE=0.000288937
theta window 14:16 est=0.7711923021126 se=2.99106e-09 delta=-2.9438648e-10 z=-0.0984 JKratioBias=+1.4812e-12 meanD=201731.64 deltaBeta=-5.9387068e-05 deltaBetaSE=0.000603275
adjacent paired late minus early = -5.6376068e-08 se=2.04195e-08 z=-2.7609
adjacent paired late minus early = -7.8071323e-09 se=8.26004e-09 z=-0.9452
```

## Frozen friend-tree production configuration

The lead froze production at 04:04:51 UTC with the `L1cap=8` extended layout, four independent training shards at (12,16), and two fixed weight files. The main estimate is the full weight at (20,24). The sane weight and early rungs are controls. I checked the metadata, byte lengths, actual weight lengths, all dropped slices, normalized plain entry, and file hashes. The given slices drop the top-class unweighted N/B terms and the constant high-degree-neighbour B/B2 terms. The frozen feature generator also passes the 372-tree independent replay, including the largest observed scaled difference `4.537e-15` and largest absolute difference `3.333e-13`.

```text
$ ../../bin/layout_info Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8
A0=0 H0=96 N0=120 B0=1740 N20=2010 B20=2370 P0=2430 Q0=3294 R0=4050 total=5058 NS1=45 NT=36 NC=12 ND=7 NR=16 L1cap=8 c1cap=4
$ training metadata, weight-mask and hash check via Python
W_1: meta=5058 1 25000 25000 1 / 12.000000 16.000000 bytes=2023600000
W_2: meta=5058 1 25000 25000 1 / 12.000000 16.000000 bytes=2023600000
W_3: meta=5058 1 25000 25000 1 / 12.000000 16.000000 bytes=2023600000
W_4: meta=5058 1 25000 25000 1 / 12.000000 16.000000 bytes=2023600000
theta_sane.f64: sha256=bc71186625e0bb8ff19d97db83e088285e7fdc9f26118541dce2444b971cbb3c entries=5059 nonzero=3006 plain=1.0 maxabs=1 dropped-block nonzeros=[0, 0, 0, 0]
theta_full.f64: sha256=7d825f142759c83b8890cb61585d5d59489bffa26c33e211188f192169eb688f entries=5059 nonzero=3207 plain=1.0 maxabs=1.0722795 dropped-block nonzeros=[111, 21, 55, 14]
$ python3 exact_tree_check.py --l1cap 8
seed=5927831 shapes={'star': 37, 'path': 37, 'prufer': 111, 'caterpillar': 111, 'friend': 74, 'two-hubs': 1, 'two-mid': 1}
layout p=5058 trees=372 model=0 hub_decay=False literal_N2_B2=False
A  exercised=  63 max_abs_a=4.615e-14 max_abs_b=0.000e+00 max_scaled_a=3.471e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
H  exercised=  21 max_abs_a=3.333e-14 max_abs_b=3.077e-15 max_scaled_a=2.703e-15 max_scaled_b=1.053e-15 fail_a=0 fail_b=0
N  exercised= 604 max_abs_a=4.762e-14 max_abs_b=0.000e+00 max_scaled_a=3.708e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
B  exercised= 183 max_abs_a=4.615e-14 max_abs_b=4.201e-15 max_scaled_a=4.068e-15 max_scaled_b=1.811e-15 fail_a=0 fail_b=0
N2 exercised= 122 max_abs_a=5.455e-15 max_abs_b=7.692e-16 max_scaled_a=2.204e-15 max_scaled_b=3.460e-16 fail_a=0 fail_b=0
B2 exercised=  48 max_abs_a=3.884e-14 max_abs_b=5.562e-16 max_scaled_a=2.557e-15 max_scaled_b=4.477e-16 fail_a=0 fail_b=0
P  exercised= 272 max_abs_a=3.333e-14 max_abs_b=0.000e+00 max_scaled_a=3.008e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
Q  exercised= 370 max_abs_a=3.333e-13 max_abs_b=0.000e+00 max_scaled_a=4.129e-15 max_scaled_b=0.000e+00 fail_a=0 fail_b=0
R  exercised= 445 max_abs_a=5.185e-14 max_abs_b=4.815e-15 max_scaled_a=4.537e-15 max_scaled_b=2.203e-15 fail_a=0 fail_b=0
```

The launch driver uses distinct declared run-index intervals within seed 9100 across workers, and seeds 9200 and 9300 for the monitor rungs. These are separate from training seeds 811–814. To audit a running directory, `production-check/audit.py` reads only full-sized shards with finite data and exact expected run indices and checkpoint order. It checks actual completed streams against the training and pilot manifests and uses paired jackknife differences when a set of trajectories supplies both windows. The final manifest must be checked after production ends.

The originally supplied `prod_summary.py` can treat a parseable partial chunk as complete. A missing last checkpoint within a non-last file shifts the concatenated checkpoint arrays and can cause an alignment assertion even after global truncation. A whole-run prefix also passes while being labeled a finished chunk. I reproduced both using local synthetic fragments. This affects interim monitoring and any final summary that includes a failed or incomplete chunk. Full-size/ID validation or atomic publication of verified files resolves it.

```text
$ ../../.venv/bin/python production-check/lead_prod_summary.py production-check/partial --names sane,full
    table(rows, f"{label}, {nf} chunks")
  File "checks/cold-estimator/production-check/lead_prod_summary.py", line 28, in table
    for t in ts: assert np.all(by[t][:, 0] == by[ts[0]][:, 0])
                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError
$ ../../.venv/bin/python production-check/lead_prod_summary.py production-check/short --names sane,full
--- MAIN (checkpoints 12 16 20 24), 2 chunks: runs=2, checkpoints=[12.0, 16.0, 20.0, 24.0], mean C=[24.0, 32.0, 40.0, 48.0], den/C at last checkpoint: sane=1.0000, full=1.0000
```

### Assumption needed for the systematic bound

Write `R0, R1, R2` for the full-theta population ratios on (12,16), (16,20), and (20,24) and `bj=Rj-mu`. If the magnitude of the bias contracts by at most `rho = exp(-4*mu) ≈ 0.10405` on the final step, the triangle inequality gives

```text
|b2| <= rho/(1-rho) * |R2-R1|.
```

If the same contraction holds on both steps, the earlier contrast gives

```text
|b2| <= rho^2/(1-rho) * |R1-R0|.
```

These inequalities do not require a common sign. Their premise is a bound on the *total* bias for the final fixed theta. It requires evidence or explicit assumption. A bound on individual decay modes fails to imply contraction of their sum when components cancel. Also `exp(-4*mu)` describes inverse leading-scale growth across these windows, while the numerator of the bias is the difference of two residuals. The residual differences can grow or change sign. The fast early pilot decay was observed for different weights and cannot be transferred as a bound to the frozen full weight. The DAN time ladder shows how window changes can be checked against an external exact root. It does not identify the friend-tree contraction factor.

For an empirical systematic bound I would report one prespecified contrast, its paired or independent SE as appropriate, the confidence level, and the asserted contraction hypothesis. A simple conditional 95% upper value uses `factor*(abs(observed contrast) + 1.96*contrast_SE)`. The alternative contrasts are sensitivity checks unless simultaneous coverage is accounted for. A final standard error is the sampling uncertainty about the frozen finite-window target.

### Frozen production replay and summary patch

The lead revised the monitor to exclude files without the required row count and repeated run ID within each checkpoint block. The revised tool ran successfully on the still-writing production directory. The separate cold audit checks exact IDs, checkpoint order, finite data and initialized-state duplicates as well.

I replayed the first 24 frozen main trajectories at every production checkpoint with the same stored weights and the generic generator-evaluation path. Identical IDs, times, C and S compare the optimized and generic dot products on the same tree states. The largest measured change in the final-window ratio was `1.585e-12`. This check covered t=24 under the updated sampler.

```text
$ ../../.venv/bin/python ../../src/prod_summary.py ../../outputs/prod --names sane,full
--- MAIN (checkpoints 12 16 20 24), 3 chunks: runs=300, checkpoints=[12.0, 16.0, 20.0, 24.0], mean C=[1147.3, 11023.4, 105956.2, 1018240.2], den/C at last checkpoint: sane=0.7829, full=0.8654
--- RUNG 12->16, 3 chunks: runs=15000, checkpoints=[12.0, 16.0], mean C=[1150.6, 11058.0], den/C at last checkpoint: sane=0.7829, full=0.8654
--- RUNG 16->20, 2 chunks: runs=1000, checkpoints=[16.0, 20.0], mean C=[11314.4, 108746.2], den/C at last checkpoint: sane=0.7830, full=0.8654
+ ./dan-followup/ftprod 9100 0 24 production-check/frozen-generic ../../outputs/train/theta_sane.f64,../../outputs/train/theta_full.f64 12 16 20 24 Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8 rfast=0
$ ../../.venv/bin/python -B production-check/compare_rfast.py
identical run IDs, checkpoints, C and S on 96 rows
theta1 largest abs endpoint num delta=1.3e-06 den delta=2e-06 fast-generic ratio=-5.1002e-13
theta2 largest abs endpoint num delta=7.9998e-06 den delta=4e-06 fast-generic ratio=+1.58475e-12
```

I also checked a larger seed envelope before the deadline. It allows 500,000 runs per main worker starting at indices 0, 10,000,000 and 20,000,000, 3,000,000 first-rung runs and 500,000 second-rung runs. This finds zero shared initial states for any eventual manifest within those ranges together with the stated training, pilot and cold-validation inputs.

```text
$ python3 seed_check.py 601:0:20000 602:0:20000 7000:0:2000000 7001:0:300000 7002:0:30000 7003:0:3000 62863482315701:0:30000 62863482315702:0:30000 62863482315703:0:30000 62863482315704:0:30000 62863482315705:0:6000 811:0:25000 812:0:25000 813:0:25000 814:0:25000 9100:0:500000 9100:10000000:500000 9100:20000000:500000 9200:0:3000000 9300:0:500000
shards=20 runs=7599000 distinct pre-seeds=7599000 collisions=0
constructed collision: (601, 0) and (1306363514317730328, 1) map to 0x248b6a394cfb5000
identical state: 03923b727a8d68b6 ac546219739f3e26 b125d179cab8754a d1c586cd8e9136e9
```

### Frozen production at the review cutoff

Production was still running at this review's cutoff around 04:19 UTC. The cold auditor accepted 31 complete shards and verified 45,500 actual production streams against the declared training and pilot ranges. The primary (20,24) ratio was 0.5657236812 with sampling SE `2.4559e-8` on 1,500 runs. The full denominator/core normalization was about 0.8654 at all four checkpoints. The paired sane/full difference on this window was `-5.92e-9 ± 1.64e-8`. The first-order sampling bias estimate for the full ratio was `-2.20e-12`. These numbers are an interim check of the estimator and files.

The independent (16,20) control produced a conditional 95% upper bias estimate of `1.34e-8` at the assumed factor 0.10405. The single-step paired comparison on the main draws gave `2.79e-8`. The earlier rung contrast gave `1.67e-9` under two contraction steps. Those are different statistical contrasts conditioned on stated assumptions. They provide no unconditional error ceiling. The complete audit output is in `production-check/review-cutoff.log`. Final sampling digits and the final conditional upper bounds depend on the completed production shards.

```text
$ ../../.venv/bin/python -B production-check/audit.py ../../outputs/prod
initial-state audit: 45500 production runs across 31 complete shards; 2644500 total unique pre-seeds incl. training and controls
main: 1500 complete runs, mean C=[   1166.72466667   11214.67866667  107764.99       1035706.36666667]
      20->24: 0.5657306319653 SE=7.69135e-06 avg den=927941.38 projected residual at mu_ref=+6.3953398 JKratioBias=-2.691e-10
  sane den/C=[0.78315684 0.78287808 0.78291417 0.7829154 ]
      20->24: 0.5657236752785 SE=2.88165e-08 avg den=726499.73 projected residual at mu_ref=-0.047020129 JKratioBias=-6.505e-12
  full den/C=[0.86546093 0.86531263 0.86536932 0.8653724 ]
      20->24: 0.5657236811996 SE=2.4559e-08 avg den=803015.19 projected residual at mu_ref=-0.04721763 JKratioBias=-2.202e-12
r1216: 40000 complete runs, mean C=[ 1150.576825 11059.0963  ]
  sane den/C=[0.78305675 0.78292324]
  full den/C=[0.86534898 0.86536254]
r1620: 4000 complete runs, mean C=[ 11151.90825 107181.59225]
  sane den/C=[0.7828938  0.78291731]
  full den/C=[0.86533005 0.86537281]
paired main (20,24) sane minus full = -5.9210482e-09 SE=1.643413e-08
sane paired main (16,20) vs (20,24): late minus early=+3.9254871e-08 SE=1.021201e-07 conditional upper95(main bias)=2.780325e-08 if |bias_next| <= 0.1040488*|bias_prev| (factor=0.1161323)
sane independent rung (16,20) vs main (20,24): late minus early=-4.3143641e-08 SE=6.519047e-08 conditional upper95(main bias)=1.984897e-08 if |bias_next| <= 0.1040488*|bias_prev| (factor=0.1161323)
full paired main (16,20) vs (20,24): late minus early=+7.7051893e-08 SE=8.308036e-08 conditional upper95(main bias)=2.78589e-08 if |bias_next| <= 0.1040488*|bias_prev| (factor=0.1161323)
full independent rung (16,20) vs main (20,24): late minus early=-1.0046782e-08 SE=5.388135e-08 conditional upper95(main bias)=1.343119e-08 if |bias_next| <= 0.1040488*|bias_prev| (factor=0.1161323)
sane independent rung (12,16) vs (16,20): late minus early=+8.3016745e-08 SE=8.12367e-08 conditional upper95(main bias)=2.927097e-09 if |bias_next| <= 0.1040488*|bias_prev| (factor=0.01208343)
full independent rung (12,16) vs (16,20): late minus early=+6.4447581e-09 SE=6.716267e-08 conditional upper95(main bias)=1.668523e-09 if |bias_next| <= 0.1040488*|bias_prev| (factor=0.01208343)
```
