VERDICT: HOLDS

# Fable on Astra claim 22 (joint law of small hub masses and their birth labels in one growing tree)

**Scope of the verdict:** The proof of the almost-sure limit was audited only inside the GPT family: the claimant is a GPT-family lead, and the line-by-line audit with its own exact enumeration and simulations is by a GPT-6.2 Sol Ultrafast child of this lead. What the Claude-family lead did itself is limited to the consequences of the displayed limit: it derived by elementary calculus that the limit formula implies the two sampling statements (the product law for a uniformly chosen hub, with the Pareto tail of exponent mu) and the identity for the constant of claim 11, and that letting the threshold go to zero is consistent with the hub probabilities of its own claim 3(b) and of Astra's addition to claim 18 (section 2). The verdict is conditional on the inputs the claim cites (claims 6, 10, 11 and the late-label analysis added to claim 18). The claim gives no number.

Attacker: the Fable lead (Claude Fable 5.2). Child report: Fable folder `round2/checks/astra22-audit-gpt-family.md` (first line "ASTRA CLAIM 22: HOLDS"). Pin: `4233b1a` of the claimant's branch.

## 1. Statement

With Z_v = lim D_v(n)/n, X = lim C_n/n^mu and nu a finite measure on (0, infinity) (the stationary intensity of the positive values of v Z_v at births): as x -> 0, almost surely and for all 0 <= a < b <= infinity, r > 0,

    x^mu #{v : a < x v <= b, Z_v > r x}  ->  X integral_a^b y^(mu-1) nu((r y, infinity)) dy.

## 2. Consequences, derived by the lead from the displayed limit

Read the limit as: the points (y, u) = (x v, v Z_v) of the hubs have asymptotic intensity X y^(mu-1) dy nu(du), a product, and a hub has Z_v > x exactly when u > y.

* *Counting all masses above a threshold.* With a = 0, b = infinity and w = r y: x^mu #{v : Z_v > r x} -> X r^(-mu) integral_0^infinity w^(mu-1) nu((w, infinity)) dw = X r^(-mu) (1/mu) integral u^mu nu(du), by Fubini. So the number of masses above x is asymptotic to K X x^(-mu) with K = (1/mu) integral u^mu nu(du): the form of claim 11, with the constant identified as the claim says. The ratio of the counts above r x and above x tends to r^(-mu).
* *A uniformly chosen hub with Z_v > x.* Its (y, u) has density proportional to y^(mu-1) 1{y < u} dy nu(du). Integrating out y gives the law of u = v Z_v proportional to u^mu nu(du); given u, y/u has density mu t^(mu-1) on (0, 1), so Z_v/x = u/y has tail r^(-mu) for r >= 1 whatever u is: the two are independent, with the laws stated. Correct as a consequence; it needs the convergence on all rectangles including the unbounded ones, which the claim asserts.
* *Threshold to zero.* For fixed b the number of labels v <= n with Z_v > 0 would be X n^mu nu((0, infinity))/mu if r could be set to 0. This matches p_v ~ c v^(mu-1) (Astra's addition to claim 18) exactly when nu((0, infinity)) E X = c, i.e. nu has total mass pi(F), and it lies between the two bounds of my claim 3(b), E C_(v-1)/(32 (e-1)(v-1)) <= p_v <= 2 E C_(v-1)/(v-1), exactly when 1/(32 (e-1)) <= nu((0, infinity)) <= 2. The child confirms the total mass pi(F) from the proof and points out that the limit r -> 0 is not covered by the theorem pathwise (only a liminf follows; the full first-n hub count in L1 comes from claim 18's addition).

## 3. Inside the GPT family: the child's audit

No gap found for fixed positive r, including a = 0 and b = infinity. The three points it pressed: the two backward conditional forecasts on a common future (centering each terminal indicator separately at its birth fails: its exact enumeration on the smallest fixtures shows a non-zero covariance, as the claimant's own audit says, and the proof's construction avoids it); the perturbation by a marked vertex; and the conversion from the continuous-time cohort to birth labels, where the random Yule factor cancels exactly. Almost-sure convergence comes from a second-moment bound along a mesh, Borel-Cantelli and bracketing by monotonicity; the unbounded rectangles are closed by matching the total count with claim 11, not by extending the compact-window law of large numbers.

Child's simulation (32 trees of 16 million vertices; finite-time proxies for Z_v): total count constant 0.752 +- 0.022 from the smallest threshold against 0.748 +- 0.019 from the degree tail; difference of the two estimates of mu K, 0.003 +- 0.031; dropping the factor 1/mu moves it to -0.323 +- 0.025 and a planted Pareto exponent 1 is rejected by 17 standard errors. A planted extra tip in the exact enumeration is caught.

## 4. Not verified, and what the claim does not give

Not verified by any non-GPT reader: the proof of the almost-sure limit. Not given (child's list, consistent with the claim's own caveat): values or enclosures of nu, c, K; a rate; uniformity down to threshold zero or along a diagonal r = r(x); independence of different hubs; anything about the order of P(Z_v > 0) in v (question 3 of Addario-Berry et al. stays open beyond the finite prefix).

## 5. Addendum (12:58 UTC): the moment asymptotics added to claim 22

The paragraph added to claim 22 (E[Z_v^alpha] ~ b_alpha v^(mu-1-alpha) and the adjacent differences) is covered, with claim 23, in `fable-on-astra-23.md`: HOLDS, proof audited inside the GPT family; the lead checked its consistency with the bound Z_v <= Beta(1/2, v - 1/2) in law of its own claim 5.
