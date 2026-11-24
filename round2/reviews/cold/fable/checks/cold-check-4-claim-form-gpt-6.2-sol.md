VERDICT: CORRECT

Sections 1–5 of the revised version 3 state the exact identities, finite-window statistical estimates, and conditional systematic bounds correctly. The interim coefficients and numerical bounds also pass the outward check. This verdict is for the copied wording in `task4/CLAIM.verified.md`, SHA-256 `dbcdb01950103efb541654d56e13c2e2e5451022b7e6e51e455404b124063625`, and the interim moments archive identified below. The final data substitution remains for the later review.

I loaded only definitions from the two Python programs and evaluated the full-observable S1–S4 estimates and nominated contrasts, plus L3. The synthetic paired check is separate. I did not execute `final.sh`, the full `final_analysis.py` program, or any of the registered A-test analyses. All artifact writes were inside this directory.

## Exact statements and algebra

For the ideal chain and fixed finite weights, the generator identity (I1) is justified at finite time by polynomial bounds in the Yule population. Substitution gives (I2), and additivity of the numerator and denominator increments gives (I3). The positive D standing conditions now appear in section 5. They ensure that the convex weights in (I3) apply to K1 and K2.

For either contraction factor r, set `f=w+(1−w)*r` and let B denote `R(16,24)−mu`. Then `|B| <= f*|b(16)| <= f*r*|b(12)|`. With `b(12) = R(12,16)−R(16,24)+B`, rearrangement gives coefficients `f*r/(1−f*r)` for K1 and K2. Because `0<w<1` and `0<r<1`, the coefficient is increasing in both, and the denominator in the bound is positive. For K3 and K4, three contractions give `|b(24)|<=r^3*|b(12)|`. The same rearrangement against `R(12,16)−R(24,28)` gives `r^3/(1−r^3)`.

At the displayed upper values `w=0.094249`, `rho=0.1043505` and `rho'=0.3464559`, I get:

| Bound | Coefficient at the displayed limits | Printed upper coefficient |
| --- | ---: | ---: |
| K1 | 0.020093471958268 | 0.020094 |
| K2 | 0.164648645342862 | 0.164649 |
| K3 | 0.001137567986852 | 0.001138 |
| K4 | 0.043390094309554 | 0.043391 |

The exponent-derived values before rounding rho and rho' are smaller. Section 4 now states the contractions directly on finite-window absolute biases of full. It neither assumes an exact denominator growth relation nor derives a per-window contraction from a spectral envelope.

## Statistical calculation and scope

The reducer gives each run a stacked vector U of the consecutive-window pairs `(e,d)` by observable, where `e=Δnum−MU_REF*d`. Its stored `S=ΣU` and `Q=ΣUUᵀ` preserve the cross-window and cross-observable moments. For the selected window, the analysis uses selectors `ce,cd`, `D=cd·S`, `delta=(ce·S)/D`, and `v=(ce−delta*cd)/D`. The ratio is `MU_REF+delta`. Changing the reference by h replaces e by `e−h*d` and delta by `delta−h`; `e−delta*d` and the ratio are invariant.

For a delete-one perturbation of the sample mean, its first-order shift is `−n/(n−1)*(v·U_r−v·S/n)`. The linearised jackknife variance is therefore

```text
n/(n−1) * [vᵀ Q v − (v·S)^2/n].
```

The exact unlinearised leave-one ratio has a run-specific denominator. The distinction matters, and section 3 now labels the quoted SE accurately. The code handles paired contrasts by using the difference of the two influence vectors in the same Q, and independent contrasts by quadrature. The w calculation is the same first-order ratio calculation on two denominator sums. I tested the production functions on a synthetic fixture with scrambled checkpoint rows, three chunks, correlated observables and differing path sizes. The individually reconstructed sums and cross-products, ratios, paired errors and split-weight error matched. Two centring constants gave the same answers. The fixture also verifies that the code's run-level SE is first-order, since the exact leave-one calculation differs numerically.

`final.sh` specifies phase-2 main and independent rung samples with chunks at least 500. That excludes the 600 earlier trajectories that prompted the main estimand choice. The stats estimate the finite-window ratios on that held-out sample under the advertised run independence and large-sample approximation. The late ratio itself was registered before its first run. The comparison with the early rung is a conditional calculation using that fixed estimate; the wording need not imply that every later comparison was part of the original late-cohort registration. The reducer tests repeated run IDs within each archive sample. Cross-sample independence still comes from the launch manifest and correct simulation, as stated.

The approximate 95% primary bounds use two events: the upper one-sided 97.5% limit for w and a two-sided 97.5% confidence interval for the nominated contrast. The Bonferroni step does not require those two statistics to be independent. K1 and K2 use the same events. K3 and K4 need only the separate late contrast interval, so the same multiplier is conservative for a nominal 95% late statement under the normal approximation. The wording supports these separately scoped statements. A simultaneous 95% coverage claim for all primary and late bounds would need a familywise calculation.

## Interim replay and rounding

All S1–S4 estimates, errors and counts match the displayed precision. The raw upper limit for w from the code is 0.094248162543; 0.094249 is its outward-rounded advertised limit. The table below uses the outward coefficients from the text and the unrounded estimates and run errors. Using the exact normal quantile in this check still puts every bound below the stated upper number.

| Bound | Recomputed pre-display upper | Upper number in claim |
| --- | ---: | ---: |
| K1 | 1.81816074078e-9 | 1.82e-9 |
| K2 | 1.48978972732e-8 | 1.49e-8 |
| K3 | 1.12398984940e-10 | 1.13e-10 |
| K4 | 4.28568045300e-9 | 4.29e-9 |

The revised display helpers give the same four displayed numbers on these inputs. One boundary robustness note for the later substitution: the deliberate subtraction of 1e-9 in the normalized argument to the ceiling lets a value a few ulps above a decimal boundary round down to it. I reproduced that on synthetic near-grid probes. The actual interim numbers are outward and unaffected. Either use exact decimal ceiling (for example Decimal.from_float and ROUND_CEILING) or check the final displayed inequalities against the unrounded values, as I did here. In the secondary pooled estimator, all parts must share a centring reference unless explicitly recentered; the actual main1/main configuration does share one.

```text
$ taskset -c 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 120 ../../.venv/bin/python -B task4/audit.py
synthetic mu_ref=0.5657236: max scaled moment delta=7.26376e-10; ratio and paired SE pass; split-weight SE passes
  exact vs linear run jackknife: ratio 3.94247131e-05 vs 3.79893384e-05; paired 4.38284853e-05 vs 4.03554221e-05
synthetic mu_ref=0.566: max scaled moment delta=3.34398e-11; ratio and paired SE pass; split-weight SE passes
  exact vs linear run jackknife: ratio 3.94247131e-05 vs 3.79893384e-05; paired 4.38284853e-05 vs 4.03554221e-05
Changing mu_ref leaves these ratio and linearised variance calculations unchanged.
coefficients C exp(bound): rho=0.104350484754765 K1/K2=0.020093467464298 K3/K4=0.001137567487701
coefficients C displayed cap: rho=0.104350500000000 K1/K2=0.020093471958268 K3/K4=0.001137567986852
coefficients C_prime exp(bound): rho=0.346455810330057 K1/K2=0.164648557544510 K3/K4=0.043390059156962
coefficients C_prime displayed cap: rho=0.346455900000000 K1/K2=0.164648645342862 K3/K4=0.043390094309554
current archive SHA-256 e95b635c268b55affa4cc8f4c75d45b012fc585f3cbc1edf1cd9d37ee9452247
main R(16,24)=0.56572365682350 run SE=7.273846982e-09 n=13500 chunks=135, Dhat=873720.797
r1216 R(12,16)=0.56572370303439 run SE=1.836374867e-08 n=260000 chunks=52, Dhat=8552.36914
r1620 R(16,20)=0.56572364236053 run SE=1.881414132e-08 n=25500 chunks=51, Dhat=81940.5896
late R(24,28)=0.56572367025551 run SE=2.30123403e-08 n=180 chunks=9, Dhat=7784039.39
S4: w=0.094244523213 SE=1.856800734e-06 with code upper 0.094248162543; z_0.9875=2.241402727605
main: r1216 minus target=+4.621088712e-08 SE=1.975186358e-08 upper |contrast|=9.04827680291e-08
  C: using printed outward coefficient 0.020094, direct upper bias=1.81816074078e-09; updated code internal coeff 0.0200933860173, rounded printed bound 1.82e-09
  C_prime: using printed outward coefficient 0.164649, direct upper bias=1.48978972732e-08; updated code internal coeff 0.164648300342, rounded printed bound 1.49e-08
late: r1216 minus target=+3.277888194e-08 SE=2.944138365e-08 upper |contrast|=9.87688795603e-08
  C: using printed outward coefficient 0.001138, direct upper bias=1.1239898494e-10; updated code internal coeff 0.0011375674877, rounded printed bound 1.13e-10
  C_prime: using printed outward coefficient 0.043391, direct upper bias=4.285680453e-09; updated code internal coeff 0.043390059157, rounded printed bound 4.29e-09
allowed independent L3 (late full minus primary full): +1.343200517e-08 ±2.413455316e-08
rounding helper boundary probe: up(np.float64(0.020094000000000004)) = 0.020094, strict upward = False
rounding helper boundary probe: upsig(np.float64(1.8200000000000001e-09)) = 1.82e-09, strict upward = False
```

## Evidence wording outside the verdict

- Item 3 should say the next complex pair *found* by the numerical IR-0 calculation. The numerical root replay confirms the pair and suggested gap. A complete root exclusion and the absence of slow modes for every observable are stronger claims.
- Item 4 should say no additional near-one candidate was resolved on the sampled grid. The known real-root neighborhood has a measured eigenvalue about 0.9977 with SE about 0.0051. The scan itself cannot exclude off-grid crossings, unsampled high frequencies, or tail/closure errors. Its stated role as motivation for (C') is appropriate.
- Item 5 reproduces as late minus primary `+1.3432e-8 ± 2.4135e-8` and is a finite-window contrast with the right interpretation.
- Item 6's roughly 4e-9 bound was an approximate simultaneous conditional calculation from noisy data. The “about 2%” value came from an illustrative correlated family of Gaussian contrasts under an all-zero null, including a second look. It is not an exact false-positive probability for the exploratory search. Label both quantities in those terms.
- Item 7 accurately describes the block pattern as exploratory and unresolved. Practical absence of a found mechanism is compatible with the statistical assumptions being explicit.
- In section 7, a premise violation would remove the corresponding guarantee. The bounds could still happen to hold in that case.
