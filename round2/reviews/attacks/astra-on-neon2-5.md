VERDICT: HOLDS

# Astra cold attack on Neon-2 Claim 5

**Scope: the posted claim, given the cited Astra inputs.** This reviews Neon-2 Claim 5 (`neon2.md#claim-5--simultaneous-cutoff-limits--separated-total-variation-regimes-for-addario-berry-et-al-open-question-5-track-5b`), at `ft2-neon2@6950241589978a7064ade0f6a2ee97cbbdcc6a34`, `round2/neon2/q5/README.md`. The `GAP` entries in the pinned `AUDIT.md` concern earlier drafts; the third pass on the pinned proof is `CORRECT`. Astra independently checked the two engines and the TV deductions. Components include isolated surviving original cores and exclude components made only of original leaves.

| Posted conclusion | Verdict |
| --- | --- |
| Theorem 1(1): almost-sure component count asymptotic | **HOLDS** |
| Theorem 1(2): almost-sure uniform log-volume CDF convergence | **HOLDS** |
| Theorem 2(1): almost-sure rooted-shape TV convergence for `N ≤ c₀ log n` | **HOLDS**, with `c₀ = γ₀/(2 log 4)` as in the pin |
| Theorem 2(2): TV lower bound above `μ/(1+μ)`; TV tends to 1 at `α = 1` | **HOLDS**, by the count and CDF results in Theorem 1 |

## Count and log-volume: the two engines

**Engine B (`N ≤ n^γ₀`).** For signed component tests `|f| ≤ 1`, the downward-fringe score is bounded by `(26/7)W`; the score itself need not be bounded. The dirty-core comparison covers both directions across the boundary: a dirty child may change the clean component above it even when the dirty child is deleted. Truncating the own-clock ring counts before summing over dirty vertices avoids any independence assumption about dirtiness. The variance estimate gives a centered mesh error with second-moment decay `poly(j) exp(−0.7(2μ−1)j)` before the union and cutoff weights. A stopped count of promotions and attachments to surviving cores controls variation between mesh points. Attachments to already deleted vertices leave the forest unchanged. The core count on the same mesh gives the `o(1/j)` denominator error needed after weighting by `N^μ`. The fixed-cutoff limit identifies the conditional mean coefficient; summability on the mesh is what proves the almost-sure result for growing cutoffs.

**Engine A (`N > n^γ₀`).** Select the two hub sets at deterministic time `s₀` before exposing the future. The finite-horizon degree comparison sandwiches the hubs at the random size epoch by these sets. With `L = log N` and `K` the size of either selected set, the degree law gives `K → ∞` even when `n/N` grows as slowly as `logloglog n`. The mark counts have means at most `O(KL)`, so an additive tolerance `ρKL` gives conditional failure probability at most `2 exp(−cρ²KL)`. On the starting-history event `K ≥ K₀`, choose deterministic `K₀` with `cρ²γ₀K₀ > 3`; the bound is then summable over all `n`. Remove the gate eventually by `K → ∞`, and apply the bracket containing the random clock epoch only pathwise. This needs no rate of divergence and no regularity in the cutoff schedule.

For children already present at `s₀`, the independent large-cutoff reconstruction bounds the old-core contribution by `C(T_(s₀−M)) c_max(s₀) = o(KL)` for `M = ceil(16 log(t₋)/μ)` when `s₀ ≥ M`, with the founder charged separately. The young-core tail is summed over integer-time bins, not over all observation sizes. Fresh promoted branches are disjoint, and neither their roots nor their descendants can be deleted on the hub sandwich. A hub send can spoil two marks: an unused mark at that ring and an earlier mark whose branch first separates. The upper CDF bound includes all births in `(t₋,t₊]`; the lower bound uses only births by `t₋` tested at their older age at `t₊`. The inert size comparisons give limiting CDF `x` at interior rational `x`, and monotonicity gives uniform convergence including the endpoints.

Together the engines cover every fixed deterministic cutoff schedule with `N → ∞` and `N/n → 0`, including schedules crossing the boundary repeatedly. The verdict has the posted quantifier: for each such schedule, the assertions hold almost surely.

## Total variation

- **Logarithmic cutoffs.** For rooted unlabelled shapes of ordinary size at most `N`, the alphabet has at most `4^(N+1)` types, including the singleton. Localizing the size clock makes the test family deterministic on each time block. At `c₀ = γ₀/(2 log 4)`, squaring the alphabet weight for the mesh bound adds at most `γ₀j + O(1)` in the exponent, within the summability margin. The common martingale and denominator errors are weighted by total mass, and interpolation uses the vector jump bound. Components may have more than `N` vertices; the inert age mixture gives a vanishing tail above `N`, which transfers via sequential TV and then via the empirical mass on the finite alphabet.
- **Polynomial cutoffs.** The inert process leaves size `v` at rate `v − 1/D_root ≥ v − 1/2`, giving `μ_N(V=v) ≤ 1/((v−1/2) log N)` for `v ≥ 2`. Write `M` for the empirical number of components and `v_n = ceil(M loglog N/log N) + 1`. The observed sizes greater than `v_n` have inert mass at most `1/loglog N`, since there are at most `M` of them. Their empirical mass is `ν(V>v_n)`. The count result gives `log v_n/log N = μ(log n/log N−1)+o(1)`, and the uniform CDF result applies at that random threshold. This yields the stated lower bound `1−μ(α⁻¹−1)` and forces TV to 1 when `α=1`.

The Astra inputs used here are the sequential forest-age limits, core and size limits, weighted spectral gap and dirty-core coupling, common-event intermediate degree law, and tagged-leaf and one-hub comparisons. The new uniform concentration and simultaneous mark arguments were checked above. The TV regimes are separated sufficient regimes; the log-volume result is CDF convergence. The withdrawn fixed-shape intensity is outside this verdict.

## Replays and independent review

Evidence is under `round2/attacks/neon2-forest`. All four `primary/` files matched the blob IDs from the pinned commit.

- Original `python3 -B verify_q5_exact.py`: 3,481,218 marked cases, two planted negative controls, holding bounds at sizes 2–8, and 205 hub-attachment invariances plus crossing and singleton controls passed (`exact.log`).
- Original `sim_q5.cpp` built with `g++ -std=c++17 -O3 -Wall -Wextra -pedantic` and ran eight trees through size 524,288 (`simulation.log`). The simulation is diagnostic; its plug-in TV estimates have sampling bias.
- Engine A reviewer `[session ID removed]`: `large/README.md`; `python3 -I -B round2/attacks/neon2-forest/large/check.py` passed 41,113 exact rows, the late-birth and double-mismatch controls, and the retrospective selection control (`1/2 → 9/17`).
- Engine B and TV reviewer `[session ID removed]`: `small/README.md`; `python3 -I -B round2/attacks/neon2-forest/small/check.py` passed 4,362 crossings, 1,232 deleted-target invariances, 16,038 promotions, 11,484 surviving gains, and holding, shape, singleton, boundary and unbounded-score controls.
