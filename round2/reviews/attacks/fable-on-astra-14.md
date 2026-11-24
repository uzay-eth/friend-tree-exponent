VERDICT: HOLDS

# Fable on Astra claim 14 (a signed scalar equation on unmarked shapes)

**Scope of the verdict:** The Claude-family lead re-derived the transfer identity and the scalar identity s M_s(J) (1 - U_s(J)) = 1 by hand, and evaluated the signed scalar U_s(J) by its own Monte Carlo, written from the claim's statement alone: (a) at three discounts above mu both sides of the scalar identity were simulated independently and agree to 0.1%; (b) near the critical discount the simulated U_s(J) crosses 1 at s = 0.5659 +- 0.0002, against mu = 0.5657. Both are numerical tests with stated errors; they do not prove the equation, and the randomised recursion used for U has no proved variance bound. The contraction bound 0.8600511260 (which gives the existence of the bounded solution down to the discount 0.565658), the analytic continuation behind the sign for s < mu, and the constants 471522 and 1/845 were audited and replayed only by a GPT-6.2 Sol Ultrafast child, the same family as the claimant: that part is a check inside the GPT family. The claim gives no digits of mu, as it says itself.

Attacker: the Fable lead (Claude Fable 5.2; family different from the claimant's GPT). Programs and logs: Fable folder `round2/rival-attacks/code/astra14-signed_identity.py`, `astra14-signed_critical.py` and their logs.

## 1. The claim in my words, and the two identities

Shapes are finite planted trees without marks; the top has a silent slot for its parent that counts in its degree. L0 is the physical dynamics inside the shape, b = L_o/D_o the promotion rate at the top, (C f)(S) = (1/D_o) sum of f over the core children of the top, (Delta f)(S) = f(S plus one leaf at the top) - f(S), V_s = (s - L0)^-1. U_s is the bounded solution of (s - L0) U = b + C Delta U. Claim: U_s(J) - 1 has the sign of mu - s for all s >= 0.565658, where J is the bare newborn.

I checked the two identities of section 1 of `influence/zero/critical.md` myself. Write (T U)(S) for the sum of U over the fringes of all core vertices of S. Under the dynamics of S, the fringe of a vertex v sees every event except a leaf sent to v by its parent; that event has rate 1/D_parent and changes the fringe by one leaf at its top, which is the term C Delta U at the parent; a promotion at w adds the new fringe J, the term b U(J) at w. Hence L0 (T U) = T (L0 U + C Delta U + b U(J)), and with the equation for U, (s - L0) T U = (1 - U(J)) T b. Since T b is the total promotion rate, i.e. the rate of increase of the number of cores, integrating for s > mu gives T U_s = (1 - U_s(J)) (s M_s - C_0) with M_s(S) the Laplace transform of the mean core count from the isolated seed S; at S = J this is s M_s(J) (1 - U_s(J)) = 1. Both steps are correct as stated. So above mu the signed scalar is 1 - 1/(s M_s(J)), below 1, and it tends to 1 as s decreases to mu; the content of the claim beyond this is that the bounded solution exists and stays analytic down to 0.565658 (the contraction), and the sign below mu.

## 2. Monte Carlo of both sides of the scalar identity (own program)

s M_s(J) = E[number of cores at an independent exponential time of rate s] by plain simulation of the isolated seed J. U_s(J) by the randomised recursion U(J) = (1/s) E g(S_tau), g(S) = b(S) + 1{a uniformly chosen neighbour of the top is a core child Y} H(Y), H(Y) = (1/s) [g(evolution of Y plus a leaf) - g(evolution of Y)], with independent evolutions.

| s | samples | U_s(J) | s M_s(J) | s M_s(J) (1 - U_s(J)), claimed 1 |
|---:|---:|---|---|---|
| 1.5 | 3,000,000 | 0.33271 +- 0.00009 | 1.49998 +- 0.00088 | 1.00091 +- 0.00060 |
| 1.2 | 3,000,000 | 0.42298 +- 0.00012 | 1.73129 +- 0.00162 | 0.99900 +- 0.00096 |
| 1.0 | 300,000 | 0.51620 +- 0.00049 | 2.06704 +- 0.01444 | 1.00003 +- 0.00706 |

Below s = 2 mu = 1.13 the plain estimator of s M_s(J) has infinite variance; the row s = 1.0 is therefore indicative only (its error estimate is not reliable), and the product was not formed at lower discounts.

## 3. Monte Carlo of U_s(J) at and around the critical discount (own program)

Second program: leaves handled lazily (exact negative-binomial updates between explicit events), the two evolutions of a pair share the exponential time, and the choice of the core child is a Russian roulette (continue with probability min(1, 2c/D), weight (c/D)/q). At s = 1.5 it gives 0.33355 +- 0.00046, in agreement with section 2.

| s | samples | U_s(J) | claim |
|---:|---:|---|---|
| 0.5457 | 500,000 | 1.04440 +- 0.00084 | > 1 (the claim covers s >= 0.565658 only; this point is below that range and is given for the slope) |
| 0.5657 | 2,000,000 | 1.00043 +- 0.00040 | 1 at s = mu = 0.5657237; here s is 2.4e-5 below mu, so U should exceed 1 by about 5e-5 |
| 0.5857 | 500,000 | 0.96006 +- 0.00075 | < 1 |

The three values give a slope of -2.11 per unit of s and a zero of U_s(J) - 1 at s = 0.56590 +- 0.00019 (statistical error only; the error of the middle point divided by the slope). The numerical value of mu is 0.56572. So a program written from the statement of the claim, which knows nothing of mu except the discount at which it is run, puts the crossing within one standard error of mu, at a resolution of 2e-4. The largest single sample is 6.3 in absolute value and the batch standard errors agree with the plain ones, so the roulette estimator is well behaved at these sample sizes; I have no proof that its variance is finite. The expected work per sample is infinite for s <= mu (the tree at an exponential time has infinite mean size there); the runs terminate because the divergence is logarithmic.

## 4. Inside the GPT family: audit and replay by a Sol Ultrafast child

Report: Fable folder `round2/checks/astra14-audit-gpt-family.md` (first line "CLAIM 14: HOLDS"). The child

* audited `influence/zero/critical.md` and the files it rests on line by line (table in the report), supplying the arguments it needed; no substantive gap in the two-sided sign statement for the unlimited signed response;
* replayed the claimant's cap-48 pair certificate: 22,956,300 rows, total-variation bound 0.8600511260 for the operator Delta V_s C at real s >= 0.565658, including all end-point combinations;
* verified with its own exact code, on all shapes of at most 9 vertices, the generator transfers, the weight inequalities and the resolvent identities; two planted faults were rejected;
* recomputed in rational arithmetic the window factor (< 471522) and the bound -U'_mu(J) > 1/845.

The contraction bound is a replay of the claimant's program, not an independent implementation.

## 5. Not verified

By me: the contraction bound and therefore the existence of the bounded solution down to 0.565658; the analytic identity (2) of the proof and the sign for s < mu as a theorem (my section 3 is a numerical observation of that sign at one discount below mu); the constants 471522 and 1/845. My Monte Carlo resolves the crossing to about 2e-4 in s; the claim's own interval evaluation at s = 0.56572365, [0.846, 1.145], still contains 1. Nothing here gives digits of mu beyond what simulation already gives.
