VERDICT: HOLDS

# Fable on Astra claim 5 (nineteen newborn labels; mu = largest real zero of det(I_19 - K(s)))

**Scope of the verdict:** The Claude-family lead checked the capped weight calculation and read the determinant argument using Claim 2’s renewal construction, then estimated the capped kernel in finite-horizon Monte Carlo. The lead did not rebuild the rational size-cut matrices or rederive their convergence and inverse-slope bounds. A saved Gemini calculation rebuilt small cutoffs before stopping; a GPT-family child later rebuilt and checked larger finite cutoffs with its own code.

Attacker: Claude Fable 5.2 (a different family from the GPT claimant). The Gemini child’s partial work checked the capped drift and cancellation and reported size cuts 2–4 before its run stopped; see Fable folder `round2/checks/astra-claim5-recomputation-gemini-partial.txt`. The later same-family reconstruction is in Fable folder `round2/checks/claim5-independent-gpt-6.2-sol-ultrafast.md`. Astra folder `round2/fringe/finite-kernel.md` (`2ccd472`) also claims a positive series (2), rational lower/upper matrices (3), explicit convergence (7) and uniform inverse slope; the cap is on the permanent mark, not physical send rates.

## What was attacked and how

1. **Does capping the mark keep every inequality?** The construction needs a_v >= 1/D_parent(t) for all later
   t; 1/min(D_parent(birth), 20) >= 1/D_parent(birth) >= 1/D_parent(t). Yes. The weight estimate needs
   (1 - c/D)[h(c+1) - h(c) + 1/(1 + d min(D,20))] <= q h(c): for D < 20 this is the uncapped inequality
   (checked exactly, see fable-on-astra-2.md); for D >= 20 the newborn factor is 1/4 and the left side is at
   most 0.10 (c = 1), 0.18 (c = 2), 0.17 (c = 3), 0.25 (c >= 4) against q h(c) = 0.36, 0.306, 0.281, 0.252, and
   negative for c = 0. So R W <= (51/100) W survives (the case c >= 4 has slack 0.002 only).
2. **"Largest real zero of det" versus "spectral radius 1".** K(s) is a non-negative 19 x 19 matrix with
   positive entries, decreasing in s; for s > mu its spectral radius is < 1, so det(I - K(s)) != 0; at s = mu
   the Perron eigenvalue is 1. Hence mu is the largest real zero on (51/100, infinity). Correct as stated.
3. **Monte Carlo estimate of the 19 x 19 matrix.** Same own simulator as for claim 2 with the
   capped marks and labels (`RES_VARIANT=4`), horizon 14, 80,000 families, geometric tail (the capped residual
   discounted reward has a fitted last-bin ratio of 0.767 at the central discount):

   | s | 0.5557 | 0.5657237 | 0.5757 |
   |---|---|---|---|
   | r(K(s)) | 1.0344 +- 0.0026 | 0.9995 +- 0.0025 | 0.9672 +- 0.0024 |

   The central estimate is consistent with radius one near the reference estimate 0.56572366.
   A local interpolation gives s ≈ 0.5656 ± 0.0008 (sampling error only). The truncated estimate at
   s = 0.5657237 is 0.98350(238); it estimates the radius of a time-truncated kernel, not a certified lower bound.
   Dividing by the eigenvector coordinate of the auxiliary mark-zero row in the extended Monte Carlo
   matrix gives 0.354, 0.464, 0.543, 0.649, 0.761, 0.867 at j = 2, 3, 4, 6, 10, 20, respectively.
   The mark-zero row is not one of the nineteen labels.

## What was not recomputed by the lead

The lead did not rebuild the finite rational matrices (3) or rederive (7). The stated size-cut error decays as N^(-1/24). In this Monte Carlo sample, the fitted time tail contributes about 1.7% of the sum of the extended kernel entries and about 1.6% of the tail-corrected radius at the central discount. These time-tail percentages do not measure the size-cut error. Neither this simulation nor the finite matrices tested by the children yield an improved enclosure for μ.

## Commands (from Fable folder `round2/`; use its NumPy environment, and create `bin/` and `outputs/residual/` if needed)

    g++ -O2 -std=c++17 -o bin/residual_sim src/residual_sim.cpp
    RES_VARIANT=4 bin/residual_sim 61 16 5000 14 outputs/residual/cap20_T14.txt ; .venv/bin/python src/residual_eval.py outputs/residual/cap20_T14.txt

    horizon T=14, batches=16, families per batch=5000
    s=0.5557000: r(K_s) truncated at T = 1.01508 +- 0.00247;  with geometric tail (ratio per unit time 0.7750) = 1.03444 +- 0.00257
    s=0.5657237: r(K_s) truncated at T = 0.98350 +- 0.00238;  with geometric tail (ratio per unit time 0.7673) = 0.99950 +- 0.00246
    s=0.5757000: r(K_s) truncated at T = 0.95395 +- 0.00230;  with geometric tail (ratio per unit time 0.7596) = 0.96719 +- 0.00236
