VERDICT: HOLDS

# Fable on Astra claim 16 (exponential convergence of Long's generation sums)

**Scope of the verdict:** The claim at commit `f962a8e`: for the tree grown from one newborn core with an inert parent, the discounted generation sums satisfy A_m(s) = c(s) zeta(s)^-m + O(R^-m) locally uniformly for real s near mu, with zeta(mu) = 1, and the successive ratios converge exponentially. A GPT-6.2 Sol Ultrafast child (the claimant's family) audited the proof line by line and checked the depth transfer exactly: a check inside the GPT family. The Claude-family lead simulated the generation sums with its own program at five discounts and compared them with the claimed form and with independent quantities from its attack on claim 14. A simulation can show that the sums behave geometrically at the sampled discounts; it cannot establish the remainder O(R^-m) or the local uniformity.

Attacker: the Fable lead (Claude Fable 5.2). Files: Fable folder `round2/rival-attacks/code/astra16-generation_sums.py` and log, `round2/checks/astra16-audit-gpt-family.md`.

## 1. Own simulation (lazy leaves; the tree is run to an independent exponential time of rate s, and the cores are counted by genealogical depth)

A_m(s) is the expected number of cores of depth m at that time. Ratios A_{m+1}/A_m, with standard errors in the last digits:

| s | samples | A_1 | A_2/A_1 | A_3/A_2 | A_4/A_3 | A_5/A_4 | A_6/A_5 | read-off 1/zeta(s) |
|---:|---:|---|---|---|---|---|---|---|
| 1.2 | 3,000,000 | 0.4107(4) | 0.4400(10) | 0.4351(21) | 0.4378(43) | 0.441(9) | 0.444(17) | 0.437 |
| 0.9 | 1,500,000 | 0.5609(7) | 0.6013(19) | 0.5971(38) | 0.5949(73) | 0.595(14) | 0.590(24) | 0.596 |
| 0.7 | 800,000 | 0.7411(13) | 0.7916(33) | 0.7891(69) | 0.795(14) | 0.803(29) | 0.82(6) | 0.79 |
| 0.6 | 500,000 | 0.8804(19) | 0.9377(48) | 0.9270(98) | 0.933(20) | 0.942(40) | 0.95(8) | 0.93 |
| 0.5657 | 500,000 | 0.9424(21) | 1.0079(53) | 1.007(11) | 1.022(27) | 1.04(6) | 1.06(12) | 1.00 |

At every discount the ratio is constant within the errors from the second generation on, which is what a simple pole in the depth variable with a well separated remainder gives; the first ratio differs from the limit by 0.03 to 0.06, the second by less than 0.005. At the numerical value of mu the ratios are 1 within the errors (A_m about 0.95 for m = 1, 2, 3), in agreement with zeta(mu) = 1 and with Long's two-sided bounds on A_m(mu). The limit ratio decreases with s as it must (1/zeta(s) is Long's pressure r(s)).

Cross-checks with independent quantities: the sum over all depths is s M_s(J): 1.7314(18) at s = 1.2 and 2.387(10) at s = 0.9, against 1.7313(16) at s = 1.2 from the different program of my attack on claim 14; and 1 - 1/sum = 0.4224 at s = 1.2, against the signed scalar U_s(J) = 0.42298(12) measured there. Below s = 2 mu the per-depth counts have heavy tails (largest sampled tree: 2 million cores at s = 0.6), so the errors of the deep generations are large and the table stops at depth 6.

The child's own simulator gives 0.438, 0.597, 0.79 at s = 1.2, 0.9, 0.7 and a depth-eight ratio 1.0026(72) near mu.

## 2. Inside the GPT family: audit by a Sol Ultrafast child

Report `round2/checks/astra16-audit-gpt-family.md` (first line "ASTRA CLAIM 16: HOLDS"): the exact transfer of depth to the nineteen-label matrix holds (old-child copies and direct births carry one depth factor each, intrinsic changes none; 3,831 exact rational rows matched, six planted transfer errors caught); the inert founder needs top mark zero, which the proof's sentence on general seeds should say explicitly; the positive linear depth term excludes other poles on the critical circle, and a fixed complex contour gives the uniform exponential error and the statement on the derivative of the ratio; the object is the same A_m as in C. D. Long's draft and Long's priority for the generation sums and their critical crossing is stated correctly in the claim.

## 3. Not verified, and what the claim does not give

By a second family: the proof. By anyone: a value of R or of the size of the discount window, so the claim gives no numerical rate and no digits of mu (as it says: "numerical isolation of the physical exponent to eight digits remains open"). It concerns Astra's nineteen-label matrix with a depth variable; it does not give a spectral gap for Long's own full-environment operator, whose essential spectral radius is 1 at mu.
