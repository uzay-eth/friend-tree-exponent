VERDICT: HOLDS

# Fable on Astra claim 11 (ranked limiting hubs; degrees at intermediate scales)

**Scope of the verdict:** The Claude-family lead recomputed the pooled and per-trajectory intermediate-degree statistics from finite simulations at one time. The ranked-mass, intermediate-degree and uniform moderate-degree proofs were read only by a GPT-family child, taking Claims 3, 6 and 7 as inputs. The finite physical-identity checks shared with Claim 10 do not recompute these limits.

Attacker: the Fable lead (Claude Fable 5.2; family different from the claimant's GPT) for the numerical test of
section 1, with own simulator and counting code (Fable folder `round2/attacks/astra10-11/`). Section 2 is a proof
audit by a GPT-6.2 Sol Ultrafast child (same family as the claimant) and is labelled as a check inside the GPT family.

The intermediate-degree prediction is k^μ N_{≥k}(n)/C_n → A when k→∞ and k/n→0, after cancelling X using C_n/n^μ→X. The ranked claim instead concerns x^μ #{v: Z_v≥x} → AX as x↓0 for the limiting masses. Taking k=xn at fixed x does not meet k/n→0; the proof compares the regimes separately. The uniform moderate-degree estimate (2) is an additional statement in the proof file and was not tested numerically.

## 1. Degrees at least k on single trajectories

3000 trajectories of my continuous-time simulator at t = 20 (ordinary size of order 10^9 in expectation; mean core count 108,388):

    k                          8       32      128     1024    4096    2^14    2^16    2^18    2^20    2^24
    k^mu sum N_{>=k} / sum C   0.8869  0.7803  0.7481  0.7357  0.7337  0.7335  0.7335  0.7328  0.7337  0.7323
    per trajectory, mean       0.8868  0.7802  0.7481  0.7359  0.7349  0.7354  0.7380  0.7418  0.7515  0.8134
    per trajectory, sd         0.0047  0.0068  0.0109  0.0187  0.0274  0.0417  0.0611  0.1023  0.1591  0.4554

The pooled statistic stays between 0.7328 and 0.7337 at the displayed levels from 2^12 to 2^20 (pointwise SE 0.0005 at 2^14); 0.735 was an earlier rough extrapolation from the degree histogram. The unweighted means are separate and rise at the upper thresholds. At k=2^16, the mean count is 149.8; a simple independent-count benchmark would give relative noise 1/√149.8 ≈ 0.082, compared with the observed relative spread 0.0611/0.7380 ≈ 0.083. This comparison does not identify the source of the scatter. The spread is larger at 2^24.

See Claim 10 for the prepared input and evaluation command; the tracked readout is Fable folder `round2/rival-attacks/code/astra10-11-eval_16_20.log`.

## 2. Inside the GPT family: proof audit by a Sol Ultrafast child

Report: Fable folder `round2/children/astra10-11/REPORT.md` and tracked `round2/checks/astra10-11-audit-gpt-family.md`. Line-by-line audit of Astra folder `round2/fringe/hub-spectrum.md` (`cedd91d`), taking
the cited earlier results as inputs: every step "follows" or "follows with a supplied argument"; no gap found.

* The exponent budget of the Borel-Cantelli step was recomputed: with delta = 2 mu - 1, a = delta/8,
  g = eta delta/100, b = 2g the variance exponent is at most 0.5343 delta < delta and the mean condition reduces to
  (1 + mu)/100 < 1/8; both hold with room and use no numerical value of mu (the audit notes that the variance bound
  is even better than the 4a quoted in the proof).
* The passage from fixed thresholds to ranked masses uses generalized inverses with ties, the upper and lower
  Yule/Gamma comparisons for the future normalized degree of a vertex, and a bound on later-born vertices; all
  checked. The intermediate-degree statement holds on one probability-one event for every admissible sequence k.
* The shared enumeration through 12 vertices checked the finite physical identities described in Claim 10, with two detected plants; it does not test the ranked or uniform limits.

## 3. Not verified

The earlier theorems used as inputs (claims 3, 6, 7) were attacked separately and not re-derived here. The
uniform statement (2) with its exponent gamma was audited by the child only, not tested numerically.
