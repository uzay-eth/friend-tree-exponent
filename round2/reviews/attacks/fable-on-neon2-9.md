VERDICT: HOLDS

# Fable on Neon-2 claim 9 (E Z_1^p = E Z_2^p > E Z_3^p > ... for every real p ≥ 1; per-tree gap at least Gamma(n+2)Gamma(p+1) c*(p)/(n(n+1)Gamma(n+2+p)), c*(p) > 1/4000)

**Scope of the verdict.**
- Attacker: the Fable lead, Claude Opus 5.5. Pin: `fa615a1` of the claimant's branch. The `all-moments-order` directory is unchanged from that pin to the claimant's present head (`git diff --stat` is empty).
- Re-derived by me: the analytic skeleton (section 1).
- Checked only by a GPT-6.2 Sol Ultrafast child of this lead (report `round2/checks/neon2-claims8-9-audit-gpt-family.md`, "Claim 9 (`fa615a1`): HOLDS"):
  - the 19 rational p-intervals tiling [1, 12] and their 115,964 rows;
  - the analytic tail p ≥ 12;
  - the first-ring and hazard bounds.
- The child checked the two repairs recorded in the claimant's `AUDIT.md`, where an earlier pass had said GAP. Astra's verdict is in `astra-on-neon2-9.md`.

## 1. Steps re-derived by the lead

**Normalization.** With S ~ Gamma(N, 1) independent of the jump chain:
- E_T[W_v^p] = (Gamma(N+p)/Gamma(N)) E_T[Z_v^p];
- M_p(T, v) := E_T[W_v^p]/Gamma(p+1);
- e^{-pt} M_p(T_t, v) is a closed martingale. Uniform integrability follows from the integer Yule moments of order ceil(p) + 1.

**Lower function** f_p(L) = 1{L ≥ 2} Gamma(L−1+p)/(Gamma(L−1)Gamma(p+1)).
- A gain multiplies f by 1 + p/(L−1), at rate at least L. A loss at L ≥ 3 lowers f by f·p/(L−2+p), at rate at most 1.
- So for L ≥ 3 the drift is at least p f [L/(L−1) − 1/(L+p−2)], which is at least p f exactly when p ≥ 1.
- At L = 2 (where f = 1) the drift is at least 2p − 1 ≥ p.
- So f_p ≤ M_p for p ≥ 1.

**Upper function** W_p(D, L) = Gamma(x+p)/(Gamma(x)Gamma(p+1)), with x = (D+L)/2.
- The gain rate is r_v ≤ L + (D−L)/2 = x, since core neighbours have degree at least 2. A gain raises x by one and multiplies W by 1 + p/x. Promotions lower x by 1/2.
- So the drift is at most p W, and M_p ≤ W_p.

**Founder gap.** With label 2 the centre of the three-vertex tree half of the time:
- E Z_2^p − E Z_3^p = (M_p(centre) − M_p(leaf)) Gamma(p+1)/(2 (3)_p);
- the claimed (1/((p+1)(p+2)))(1 − 2 sqrt(p+1)/(sqrt(3)(2p+1))) follows from M_p(centre) ≥ f_p(2) = 1 and the stated leaf hazard bound. That bound is checked by the child, not by me.

**Consistency.**
- At p = 1 claim 9 reduces to claim 8's scale (see `fable-on-neon2-8.md`).
- Against Astra claim 27 (constant 9/5000 on the same Gamma scale): both are lower bounds. Claim 9's constant is larger on [1, 12]. Its analytic floor is weaker just above 12 and larger from 13 on. No contradiction.

## 2. Checked by the GPT-family child only

- The tiling of [1, 12] by six `full`, four `vx_only` and nine `1v` closed intervals.
- Affine interpolation in (rates, p) inside each box.
- The p = 12 seam inequality 15974400 < 15976009.
- The replay: `ALL_REAL_MOMENTS_EXACT_CHECKS_PASSED`, with 14 out of 14 bundled controls rejected.
- A planted one-tick change at eta1[15][2]: rejected.
- Exhaustive screens of all 47 free trees with n ≤ 8 at ten values of p: no violation.
- The child also found one wrong code comment: verifier line 530 says W_p, but the loop checks f_p.

Not covered: 0 < p < 1 (Astra claim 27 claims it); stochastic dominance.
