VERDICT: HOLDS

# Fable on Astra claim 12 (positive scalar return equations)

**Scope of the verdict:** The Claude-family lead recomputed the nineteen-label matrix return and radii from finite-horizon Monte Carlo at three discounts. This is a consistency check near the reference, not a test of exact crossing or an entire half-line. Only a GPT-family child read the two scalar proofs and the dormant-response transcendence argument, checked exact finite cuts and replayed Astra’s binary all-depth data. Claims 2 and 5 enter the proof; Claim 9 enters the finer contraction route only.

Attacker: the Fable lead (Claude Fable 5.2; family different from the claimant's GPT) for the numerical test of
section 1 with its own Monte Carlo of the nineteen-label kernel (Fable folder `round2/src/residual_sim.cpp`,
`round2/attacks/astra12/scalar_return.py`). Section 2 is a proof audit by a GPT-6.2 Sol Ultrafast child (same family
as the claimant) and is labelled as a check inside the GPT family.

## 1. Monte Carlo of the kernel

My simulation of the residual process with marks capped at 1/min(D, 20) (80000 residual families, horizon 14
plus a geometric tail, delete-one-batch jackknife) gives the 19 x 19 matrix K_s; with k_s its column of label 2
and C_s the same-size matrix with that column zeroed:

    s            F(s) = e_2^T (I - C_s)^-1 k_s     r(C_s)         r(K_s)
    0.5557       1.398 +- 0.043                    0.921 +- 0.003  1.0344 +- 0.0026
    0.5657237    0.996 +- 0.020                    0.885 +- 0.003  0.9995 +- 0.0025
    0.5757       0.798 +- 0.011                    0.851 +- 0.002  0.9672 +- 0.0024

The three estimates of F are consistent with a crossing near the reference estimate 0.56572366. The estimates of r(C_s) are below one at all three sampled discounts, two of which lie on the claimed half-line s≥0.56564. The adjacent secants for F are about −40 and −20; three points do not determine the local derivative or convexity. Separately, the determinant identity det(I − K_s) = det(I − C_s)(1 − F(s)) holds when I − C_s is invertible.

From Fable folder `round2/`, run `.venv/bin/python rival-attacks/code/astra12-scalar_return.py outputs/residual/cap20_T14.txt`. The raw input is a working output; the simulator command is in the attack on Claim 5.

## 2. Inside the GPT family: proof audit by a Sol Ultrafast child

Report: Fable folder `round2/children/astra12/REPORT.md` and tracked `round2/checks/astra12-audit-gpt-family.md`. Audit of Astra folder `round2/finite/gauge/physical-trigger.md`, `round2/finite/gauge/dormant-return.md` and `round2/finite/gauge/isolated-transfer.md` (`67a1968`):

* The determinant identity and the first-return (excursion) decomposition of F hold; positivity of the path series
  gives strict monotonicity of F and hence a unique crossing, at mu.
* The finer contraction route uses the upper endpoint of Claim 9. A separate rational certificate uses the wide bracket alone and already gives contraction at 0.5641, extending to larger discounts by positivity and monotonicity; the stated half-line does not require Claim 9.
* The binary active/dormant split preserves the physical first moments; the response formula
  k_A + k_beta (1 + d_A)/(1 - d_beta) follows. The child checked the residual-constant arithmetic after replaying the binary all-depth inputs with Astra’s code:
  11562.89 < 11600 at exponent 0.56470 and 4799.18 < 4800 at exponent 0.5630.
* The transcendence argument concerns the dormant response function over C(s); it gives no information on the
  arithmetic nature of mu, as the claim itself says.
* Exact finite checks on size-cut matrices (own code) of the determinant identity and of the excursion sum passed;
  two planted errors were detected.

## 3. Not verified

The numerical test covers three finite-horizon estimates with Monte Carlo errors of about 1–4% in F; it does not certify either scalar equation. The binary model and its constants were audited only inside the GPT family. As the claim says, the responses
still sum over unbounded descendant trees: nothing here gives digits of mu.
