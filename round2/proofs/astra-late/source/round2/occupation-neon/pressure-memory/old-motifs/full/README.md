> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# All-degree positive return: direct resets and virgin copies

**Two complete chronological classes of the genuine positive trigger return are bounded without any arm-size cut.** Persistent information doubles the lower bound on the single virgin-copy class. At least `0.6777666584` of the critical return still lies outside even the class allowing arbitrarily many direct resets before that copy. Neither exponent endpoint moves.

## Which return is counted

Use the original [nineteen-mark positive trigger](../../../../../../../astra/source/round2/finite/gauge/physical-trigger.md), starting at `J₂` (the subscript is the mark index, so the killing rate is `1/2`). Its positive chronological paths have intrinsic moves, old-child copies, and nontrigger direct newborn productions; the terminal reward is a degree-two promotion. `F(μ)=1`, and `F` is decreasing in the trial discount `s`. In the class `C₀` no old-child copy occurs on the selected lineage, with any number of direct newborn resets. In `Cᵥ`, a lineage stays in its original particle until one old-child copy while that child is still virgin, and then earns its trigger before the first intrinsic change of the copied particle. It uses no direct newborn reset on that lineage. The shifted copy is impossible in `Cᵥ`: `Y⁺` is not virgin and intrinsic growth cannot make it virgin again. The extended class `C*` allows any number of nontrigger direct resets followed by either an immediate trigger or the `Cᵥ` episode. All these statements refer to the selected lineage; other copies or direct newborns may be emitted, and are not killing clocks on the original particle. Intrinsic promotion, in particular, is always retained.

## Root and child phases, arbitrary arms

For one original particle with exterior mark index `i`, use discount `λ=s+1/i` and root census `(D,k)`, with `k` core children and `L=D−k−1` leaves. Root gain rate is `L+Σ 1/D_y`; promotion rate is `b=L/D`. At each actual promotion start a separate tagged-child episode, at the **post-promotion** root census `(D,k+1)`. The child birth mark is `a=1/min(D,20)`. While that child is still virgin it sends `1/2`, and its own first change has hazard `3/2+1/D` (one leaf gain, a half-rate promotion and the actual parent send). It permanently exits the episode at this first change. Thus, conditional on the full past, its expected future contribution is the discounted integral of `a−1/D` until the exit, times `h_j=1/[2(s+1/j+3/2)]`, with `j=min(D_birth,20)`. Keep the same root gain in the holding denominator and successor. The admissible gain interval in this phase is `[L+1/2,L+k/2]`. The matched no-memory version uses `[L,L+k/2]` with the same exit and source. Other arms are arbitrary: their internal transitions are integrated out into the adapted interval, not killed. A subsequent root promotion keeps the old tagged episode and starts a new one. Branch additivity of expectation permits the sum over actual births without any independence of episodes or their environments.

For `k≥K` crown samples, put the indicator on the **copy-producing** root at the time of the copy, not on its birth. The same child phase then measures a physical discounted occupation under the positive trigger process, weighted by the old-child copy and trigger factors; it is not an unweighted physical Palm crown probability. The identity

```
a − 1/D(t) = Σ_{m>j} 1{D(t)≥m}/[(m−1)m]
```

expresses the rate in terms of actual monotone degree thresholds, including all thresholds above the mark saturation 20. In the tagged-child phase root degree occupation couples to its exact send and disappearance clock. A threshold occupation can equivalently be split at its first crossing while the child is still alive. It cannot be inferred from the frequency with which one samples a large root. The examples below measure the whole indicated class; no sampling frequencies enter the bounds.

## High-degree entrances and every exit

No `D=H` crossing is discarded. The discounted child continuation at every `D≥H` is between `(a−1/H)/(λ+3/2+1/H)` and `a/(λ+3/2)`, directly by monotone degree and the hazard bounds. For a `k≥K` objective the lower expression is used only once `k≥K`, otherwise zero; the upper expression pays the whole tail even if not yet eligible. Thus it includes every future eligibility entrance after the crossing.

In the parent table the tail fee per promotion is bounded by the corresponding child bound with `j=20`, times `h₂₀`. For the lower crown fee use it only when the parent has `k≥K−1`. The function `b=L/D` satisfies `A b≥−b/6` at `k=0` and `A b≥0` once `k≥1`, and these lower bounds remain valid under the true incoming arrows. Indeed a gain increases `b` by `(k+1)/(D(D+1))`, a promotion decreases it by `1/D`, and gains occur at rate at least `L`. Conversely `A(1−b)≥−(1−b)` since gain rate is at most `D` and promotions increase `1−b`. Hence from an actual high-degree state,

```
b/(λ+η) ≤ E ∫ exp(−λt) b(t) dt ≤ 1/λ − (1−b)/(λ+1),
η = 1/6 if k=0, else 0.
```

These bounds account for **all** remaining promotions of the original particle after the first `H` crossing, including births from a bare root that later gains leaves. They apply also to direct-newborn fees after saturation. Below the boundary both tables solve backward (degree or arity increases in every retained action), checking the two gain extrema of each linear fractional quotient and rounding outward. Child exits to nonvirgin carry zero **by definition of this class**. Their possible contribution to other classes of `F` is not estimated to be zero.

## Positive direct-reset hierarchy and remaining error

Let `A` be the true first direct-reset matrix without any old-child-copy traversals; its rows are the expected discounted promotions of the original root restricted to `D>2`, grouped by newborn index. Then `C₀=(I−A)⁻¹h` and `C*=(I−A)⁻¹(h+q)`, where `q` is the per-mark virgin-copy episode. The parent recursion applies `A` to any nonnegative vector without constructing whole-tree rows. A certified upper on `‖A‖∞` comes from mark 20 and the unit payoff: the intrinsic path is mark-independent, while original-root death decreases with the index. At the lower trial `α≤.658088622`. For a partial upper vector `w`, pay `ε=max_i(b_i^+ + (A w)_i^+−w_i)_+/(1−α)`, where `b^+` is the corresponding source bound. Then `w+ε1` supersolves the full hierarchy. This includes all uncounted direct-reset generations rather than dropping them at iteration 12. It does not pay two-copy paths or a first copy on a nonvirgin arm. For arbitrary continuation scores, the same first-production decomposition uses one conditional value per initial side arm plus birth forcing for the new arms, but its child rate/mark value must be bounded for **both** `Y` and `Y⁺`; extending by `q` alone does not price that joint residual.

## Bounds and scale

Outward display at `s=.565676636772`, with `H=48` for `Cᵥ` and `H=24` for the direct-reset hierarchy:

| Class | Lower | Upper |
|---|---:|---:|
| `Cᵥ`, matched no-memory envelope | `.0004948467` | `.0011621039` |
| `Cᵥ`, persistent virgin donor | `.0010051735` | `.0011621039` |
| `Cᵥ`, root `k≥3` at copy | `.0001125926` | `.0002255844` |
| `Cᵥ`, root `k≥5` at copy | `.0000129203` | `.0001041967` |
| `Cᵥ`, root `k≥8` at copy | `.0000002613` | `.0000927175` |
| `C₀`, all direct generations | `.2848637114` | `.3196038258` |
| `C*`, all direct generations then optional virgin copy | `.2864821270` | `.3222333416` |

The upper generation remainder for `C*` is below `.000175510`; the larger `.03576` enclosure width is predominantly physical receive relaxation, not the last generation. By monotonicity, since the lower trial is below `μ`, at least `1−.3222333416=.6777666584` of the critical return comes from the complementary paths. This is an all-degree quantitative coverage barrier for these classes, not an error upper bound for the complement or a bound on higher motif charges. Both endpoint trials give positive subreturns well below one; without independent bounds on the rest of `F` they do not resolve a trial sign. No stationary motif credit or endpoint gain is booked.

## Replay

From the repository root:

```sh
taskset -c 8 nice -n 17 timeout 50s env OPENBLAS_NUM_THREADS=1 ./round2/.venv/bin/python -Bu round2/occupation-neon/pressure-memory/old-motifs/full/first_copy.py --hierarchy 12 --extend-copy --broad 3 5 8
```

Passed eight independent whole-fringe rate fixtures using the original physical transition enumerator, seven rational threshold identities, nested-degree controls, and inward matched comparison for the virgin donor. The recurrence, tail prices and generation residual calculations use exact fractions and directed rounding. No D18 row scan is needed: this route evaluates a chronological subreturn of the original trigger, rather than a subset of the eight signed charge families.
