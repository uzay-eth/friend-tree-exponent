VERDICT: HOLDS

# Fable on Neon-2 claim 11 (almost surely 5000/5067 ≤ liminf H(T_n)/log n ≤ limsup ≤ 1.462248 for the height, twice these for the diameter, and an explicit finite-n tail bound for the core height)

**Scope of the verdict.**
- Attacker: the Fable lead, Claude Opus 5.5. Pin: `93484ba` of the claimant's branch.
- My own checks are in section 1. The line-by-line audit, the replay, the plants and the simulation are by a GPT-6.2 Sol Ultrafast child of this lead (report `round2/checks/neon2-claims10-11-audit-gpt-family.md`, "Claim 11 (`93484ba`): HOLDS"); this is cross-family to the Gemini claimant.

## 1. Lead

**Upper constant.**
- The potential Φ_z = Σ_{cores} z^{depth} f(D, L) grows at rate at most λ, with z = 31/10 and λ = 165439/100000.
- So P(H_core ≥ c log n) ≤ K n^{−(c log z − λ)}. This decays exactly when c > λ/log z = 1.46224758..., which the claim rounds up to 1.462248.
- The relations H_core − 1 ≤ H ≤ H_core + 2 and H_core ≤ diam ≤ 2H_core + 2 then give the height and diameter upper constants, 1.462248 and 2.924496.

**Lower constant.**
- Depth k is first reached near time β k, where β is the speed of the embedded branching random walk of newborn cores.
- At time log n the depth is about log n/β, so β* < 5067/5000 gives the height lower constant 5000/5067.
- Two independent surviving branches double it for the diameter, giving 10000/5067.
- The time-to-size conversion is τ_n/log n → 1, because n e^{−τ_n} → W > 0 for the Yule size.

**Scope check.** The claim does not assert convergence of H/log n to a limit, and nothing here gives one.

## 2. GPT-family child

**Audit.**
- Both rate extremes on all 5,148 finite rows suffice, because the residual is affine in the gain rate.
- Beyond the cutoff only the newborn term z f(2,1) ≤ λ is positive.
- The CMJ embedding: only the parent can promote a leaf; the driving sets are disjoint and fresh marked clocks are used at stopping times. This gives the history-independent rates.
- The hypotheses of Biggins' theorem hold: i.i.d. reproduction point processes, bounded offspring, mean > 1, finite tilts.
- The repaired two-branch binomial argument holds.
- The child also wrote a Biggins-free block-process argument for the lower endpoints (`tilted-block.md` in the child folder; the child's own supplement, not part of the claim).

**Replay.** `HEIGHT_DIAMETER_EXACT_CHECKS_PASSED`: 5,148 upper rows and 7,500 CMJ rows, the exact prefactor 861558416425/171115606632, and 11 out of 11 built-in negative controls rejected.

**Plants.** Raising z by 10^-5 on the frozen table was rejected (the far-tail condition fails). A 10^-9 increase was accepted, which shows the certificate has slack at that size, not that an error went through.

**Simulation.** 1,000 trees at n = 10^5 and 10^6, compared with the asymptotic window:

| quantity | n | below | inside | above |
|---|---|---:|---:|---:|
| H/log n | 10^6 | 138 | 853 | 9 |
| diam/log n | 10^6 | 362 | 637 | 1 |

These are finite-size values. The window is asymptotic, so trees on either side of it do not contradict the claim.
