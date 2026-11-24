VERDICT: HOLDS (conditional on its stated inputs)

# Fable on Neon-2 claim 10 (certified enclosures of eight stationary core marginals: p_2, ..., p_6, pi(J), pi(2,2), pi(L_o = 0), with corollaries)

**Scope of the verdict.**
- Attacker: the Fable lead, Claude Opus 5.5. Pin: `6c8afec` of the claimant's branch.
- The claim is conditional, as it says, on three inputs:
  - Neon-2 claim 6 gold2, which this lead attacked: HOLDS, `fable-on-neon2-6.md` section 6;
  - the repaired coboundary of Neon-2 claim 3 together with positive mean growth: HOLDS, `fable-on-neon2-3.md` section 5;
  - Astra's identification of the fringe law with the stationary core law (Astra claims 4 and 10, which this lead accepted earlier).
- The degree marginals and the neighbour-reciprocal identities need that identification. The one-sided fringe alone does not record the degree of the exterior parent.
- My own contribution is in section 1. The proof audit, replays and plants were done by a GPT-6.2 Sol Ultrafast child of this lead (section 2; cross-family to the Gemini claimant). Astra's verdict is HOLDS (`astra-on-neon2-10.md`).

## 1. Lead: an independent high-precision comparison

This lead measured the stationary core-degree law earlier in the round, from its own simulation of time-differenced histograms with 3.8e9 differenced core vertices (RESULT.md section 8, used for Astra claim 3). All six marginals that can be compared lie **inside** Neon-2's certified intervals. My standard errors are 10 to 20 times smaller than the interval widths:

| marginal | lead's measurement | Neon-2 interval | distance to lower / upper end |
|---|---|---|---|
| p_2 | 0.383456(8) | [0.38337348, 0.38351583] | 8.3e-5 / 6.0e-5 |
| p_3 | 0.145227(7) | [0.14517013, 0.14528215] | 5.7e-5 / 5.5e-5 |
| p_4 | 0.080368(5) | [0.08033089, 0.08041288] | 3.7e-5 / 4.5e-5 |
| p_5 | 0.052246(4) | [0.05222046, 0.05228311] | 2.6e-5 / 3.7e-5 |
| p_6 | 0.037181(4) | [0.03715899, 0.03720440] | 2.2e-5 / 2.3e-5 |
| pi(L_o = 0) | 0.14531 | [0.14520643, 0.14540145] | 1.0e-4 / 9.1e-5 |

Finite-time histograms estimate the stationary law with some bias, so this is a consistency check, not a proof. A single deviation of a few times 1e-5 would have been visible, and none appears.

I also checked the corollary on degree two directly: pi(J) + pi(2,2) ∈ [0.38336404, 0.38352731] contains the whole p_2 interval, not just part of it.

## 2. GPT-family child

Report: `round2/checks/neon2-claims10-11-audit-gpt-family.md`, "Claim 10 (`6c8afec`): HOLDS".

**Audit findings.**
- The two snapshots bound the same direction of drift. The residual is affine in lambda, so the larger endpoint error covers the whole mu bracket even for a signed potential.
- Stopped Dynkin is dominated by the Yule size.
- The Bernstein tail divides by the largest denominator when the numerator bound is negative; this is the safe direction.
- The bare-core tail is below 1e-479.
- Both degree-two balance equations hold.

**Replays.** The verifier was built strictly from a clean archive and run with one thread on three of the sixteen snapshots (obs1 lower and upper, obs8 upper). All fields match `certificates_all.json`; each run covered 49,204,104 finite rows per endpoint.

**Plant.** The tail-script mutation 1/2 → 1/3 was rejected.

**Simulation.** 1,000 runs at n = 10^6 (1,828,810 cores): no marginal is more than 1.28 run standard errors from its nearest endpoint. At this size the law is still a finite-size one, so this is weaker than section 1.

Not covered: monotonicity beyond p_6; the tail constant; a re-proof of gold2, of claim 3's inputs or of Astra's spatial theorem.
