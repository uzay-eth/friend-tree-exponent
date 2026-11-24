> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Uniform truncation for the Gold2 FY reward and the potential comparison

The unbounded small-payer donor class does not require a degree moment. Root
the paid reward at the host. It has at most 17 neighbours and hence at most
`17·16` ordered payer–donor terms. For each of the finitely many paid
host/payer labels write `F(x)=σ(xz−A_h(x))` in ticks. With the original
piecewise-linear spline on `[0,1/2]`, let `B` be the maximum absolute slope
of all these finitely many `F` segments. The reward is `max_k F(ξ_k)−F(1/D)`.
Replacing the donor reciprocal by zero only when `D>M`, with `M≥12`, changes
this reward by at most `B/M` per term. This does not delete the absorbing
class; it keeps its incidence, its eligibility, and its entire lower
payment. For the alternative reward that already pays the constant
absorbing-class minimum, this replacement is unnecessary: the predicate
`D≥12` is itself a finite-radius test. Neither approximation is used to
change a transition rate in the Bellman proof.

In particular the capped-reciprocal host score `f_M` has fixed read radius
and `‖f−f_M‖∞ ≤ 272 B/M`, uniformly in the physical tree and host. To read
the degree of a neighbouring core, read its core incidences one level
farther out and its leaf mark. Registration weights are only needed for
proving the prices; the total credited rewards read current oriented
occurrences. The other two rewards have bounded role labels and multiplicity.
All are bounded local functions even before truncation, because the small
payer uses the reciprocal, not a factor growing with donor degree.

For the finite process from the three-vertex seed, write
`K(t)=E C(T_t)` and `k=lim e^(−μt)K(t)>0`. The oriented local mean
convergence for each fixed `f_M` follows from the [physical spatial
law](../../../../../astra/source/round2/fringe/spatial-limit.md), including its ancestor construction.
Equivalently, the empirical averages converge almost surely, and one can
multiply them by `e^(−μt)C(T_t)` and take expectations using the uniform
`L²` core bound from the [core eigenmartingale](../../../../../astra/source/round2/fringe/almost-sure.md).
The truncation error in this expected scaled sum is at most
`272 B sup_t e^(−μt)K(t)/M`; the error in `k π(f)` is at most `272 B k/M`.
Letting first `t` and then `M` tend to infinity gives the required mean
limit, with an error uniform in the observation time. This argument does
not bound the number of hosts adjacent to a high-degree donor by the
number of such donors; that counting shortcut need not hold.

For the potential, stopped Dynkin followed by dominated convergence on
each finite interval suffices: the Yule vertex count has a finite second
moment, the total jump rate is the vertex count, and the potential before
or after a jump is bounded by `m_max` times that count plus one. Thus
`Q(t)=e^(−μt)E P(T_t)` is locally absolutely continuous. Integrating
`Q′+(μ−λ)Q ≥ e^(−μt)(E R−E_error K)` gives the stated Cesàro inequality.
The boundary term vanishes because `Q` is bounded by comparability.
Since the right-hand limiting lower bound is strictly positive and `Q`
is nonnegative, first conclude `μ>λ`; only then use
`limsup_T T⁻¹∫₀ᵀ Q ≤ m_max k`. No pointwise limit for the total potential
or its generator is needed. The argument concerns the stationary limit
and passage to expectations; it uses the already certified reward prices
and signed finite-tree comparison as inputs.
