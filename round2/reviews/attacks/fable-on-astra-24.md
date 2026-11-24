VERDICT: HOLDS

# Fable on Astra claim 24 (strict stochastic order for all sufficiently late adjacent labels, at every threshold)

**Scope of the verdict:** The proof was audited inside the GPT family only: the claimant is GPT and the line-by-line audit is by a GPT-6.2 Sol Ultrafast child of this lead (report: Fable folder `round2/checks/astra24-audit-gpt-family.md`, first line "ASTRA CLAIM 24: HOLDS"). The Claude-family lead (Claude Opus 5.5) checked the consistency of the asymptotic statements with its own claim 5 (section 2). Conditional on the inputs the claim inherits (claims 6, 10, 11, 18, 22 and the "previously certified spectral constants" for the effective threshold).

## 1. Statement

There is v_0 with P(Z_v > x) > P(Z_(v+1) > x) for all v >= v_0 and all 0 <= x < 1. For v x -> infinity, p_v(x) ~ kappa v^(mu-3/2) x^(-1/2) (1-x)^(v-1/2) and p_v - p_(v+1) ~ x p_v; at bounded v x the signed profile is positive. Effective threshold 2^(2^(2^123)); a pre-placement comparison on every old tree of size n when (n+1) x >= 232.

## 2. Checked by the lead

* Against claim 5(e) of this lead (Z_v is stochastically below Beta(1/2, v - 1/2)): for fixed x and large v, P(Beta(1/2, v - 1/2) > x) ~ x^(-1/2) (1-x)^(v-1/2) v^(-1/2)/sqrt(pi) (Laplace at t = x, and B(1/2, v - 1/2) ~ sqrt(pi) v^(-1/2)). The claimed p_v(x) is smaller by the factor kappa sqrt(pi) v^(mu-1), the order of the hub probability: consistent, and it shows where the Beta bound loses.
* The gap: p_(v+1)/p_v = (1 - x)(1 + O(1/v)), so p_v - p_(v+1) = p_v (x + O(1/v)), which is ~ x p_v exactly when v x -> infinity, as claimed; at bounded v x the two terms are of the same order and the sign has to come from the profile B, which the claim treats separately.
* Combined with claim 5(c) of this lead (domination for u <= 25), domination is now known for u <= 25 and for u >= v_0; the labels in between are open (dead ends D4 to D6 of this lead describe why the local methods do not reach them).

## 3. Inside the GPT family: the child's audit

No gap: the relative equivalents hold through both endpoints (near x = 1 the exact factor (1-x)^(v-1/2) is essential); at small x with v x -> infinity the two-birth identity controls the signed gap on its own scale; one cut-off works for all thresholds by a subsequence argument; the tower follows from the inherited spectral bounds with the corrected rounding of the cylinder estimate; the pre-placement comparison holds for every old tree (independent rational margin check). The child also observes that combining the pre-placement comparison with the compact-cylinder estimate lowers the sufficient threshold to 2^(2^(2^23)) under the same inputs (not checked by the lead). A planted factor-of-two error is rejected; a 24,000-path simulation is only illustrative.

## 4. Not given

Order for the labels between 26 and the threshold; a usable threshold; the constant kappa (the child's heuristic estimate is 0.10 to 0.11).
