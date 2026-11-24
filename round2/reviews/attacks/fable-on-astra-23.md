VERDICT: HOLDS

# Fable on Astra claim 23 (fixed-shape counts in the low-degree forest at moving and at proportional cut-offs) and on the moment addendum to claim 22

**Scope of the verdict:** Checked inside the GPT family only: the claimant is a GPT-family lead and the line-by-line audit of the three proof files, with its own exact computations and simulations, is by a GPT-6.2 Sol Ultrafast child of this lead. The Claude-family lead did three small things itself (section 2): the occupation constant 2/3 of the pair, the consistency of the moment asymptotics with the stochastic bound of its own claim 5, and the reading of what the proportional-cut-off statement needs. The verdict is conditional on the earlier results the proofs cite (claims 6, 10, 11, 18 to 20, 22). No statement about mu.

Attacker: the Fable lead (Claude Fable 5.2). Child report: Fable folder `round2/checks/astra23-audit-gpt-family.md` (first line "ASTRA CLAIM 23 AND THE MOMENT ADDENDUM TO CLAIM 22: HOLDS").

## 1. Statements

(I)(a) Delete the vertices of original degree above N_n and count the components that contain an original core, by rooted shape. If N_n -> infinity and N_n/n -> 0, the count of a fixed shape R divided by (n/N_n)^mu converges in probability to X times a coefficient (2A/3 for the pair); almost surely in a range of slowly growing cut-offs; and there are deterministic cut-off sequences along which the pair count is zero infinitely often almost surely. (I)(b) If N_n/n -> c > 0, the deleted labels are eventually exactly those with Z_h > c, and the counts of finitely many shapes converge stably to independent Poisson variables with means I_R = integral_0^infinity P(T^0_t = R) dt per deleted hub (T^0 the inert newborn; I_pair = 2/3); each Z_v has no atom at a fixed positive value. (II) For fixed alpha > 0, E[Z_v^alpha] ~ b_alpha v^(mu-1-alpha) and E[Z_v^alpha] - E[Z_(v+1)^alpha] ~ (1 + alpha - mu) b_alpha v^(mu-2-alpha).

## 2. Checked by the lead

* I_pair = 2/3: the inert newborn (one core with one leaf and a silent slot, degree 2) stays a pair until its leaf rings (rate 1) or its root rings and picks the leaf (rate 1/2); it never returns to that shape; the expected time in it is 1/(3/2) = 2/3. The child's exact computation gives the same and 1/4 and 2/15 for the star and the path with three vertices.
* (I)(b) needs P(Z_h = c) = 0 for the labels involved, since otherwise the set {h : D_h(n) > N_n} need not settle. The proof obtains atomlessness by cutting a hub edge that is silent for ever and writing the share as a function of a ratio of independent Gamma variables. This is the mechanism of Neon-2's claim 1 and of my claim 2(e) (there for X, here for a single Z_v); the child checked the argument for Z_v.
* (II) against my claim 5(e): Z_v is stochastically below Beta(1/2, v - 1/2), so E[Z_v^alpha] <= Gamma(1/2 + alpha) Gamma(v)/(Gamma(1/2) Gamma(v + alpha)), of order v^(-alpha). The claimed order v^(mu-1-alpha) is smaller by the factor v^(mu-1), which is the probability of being a hub: consistent, and it shows that the Beta bound is sharp for the conditional law only. The adjacent difference is controlled by the two-birth identity (the same identity as in my addendum to `fable-on-astra-18.md`); the limit itself was checked by the child (456 exact cases, a flipped sign detected).

## 3. Inside the GPT family: the child's audit

No gap found in (I)(a), (I)(b), (II) within their stated scope. It reconstructed one of the slow cut-off sequences with infinitely many zero pair counts; checked that the Poisson limit is stable because the thinning is conditional on the early history; and simulated 2,048 trees of 10^6 vertices at proportional cut-offs: pairs per deleted hub 0.6747 +- 0.0104 and 0.6780 +- 0.0127 (limit 2/3), variance over mean 0.996 and 1.008. Simulated moments for alpha = 1, 2 have the predicted scale; they do not resolve adjacent differences. Omitting the inert slot in the occupation computation makes its check fail.

## 4. Not verified, and what the claims do not give

Not verified by a non-GPT reader: the three proofs. Not given: an almost-sure statement for all sublinear cut-offs (false, by the claim's own example); thresholds chosen from the realised masses; uniformity in alpha; a label from which the adjacent moments are ordered.
