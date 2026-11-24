> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Fixed-motif lineage and normalized signed forcing

**Result.** There is an exact birth-lineage representation of all eight Neon motif measures and their joint reciprocal pressures. A single normalized renewal quotient carries the signed forcing, with aggregate rather than per-exit error. The outer renewal can be truncated to `1e−8` **if the intra-cycle rows are known**. We do not obtain that accuracy for the full signed forcing: the projected motif census and its current pressure are not Markov, and the generic core-count bound for evaluating the missing rows is prohibitively slow. No exponent endpoint is claimed here.

## 1. Motifs have unique birth events

Use the [oriented uniform-core law](../../../../../astra/source/round2/fringe/spatial-limit.md) `π`, with one parent toward the end. For a rooted connected core motif `M`, let `π_M` count each physical embedding at its highest vertex, with siblings unordered unless they have distinct roles. Count all geometrical embeddings for now, allowing the test to enforce its real domain. At each promotion of a leaf to a new core `v` by `u`, every new connected motif consists of `v` and a unique old connected set containing `u`; `v` is a terminal tip at that instant. Its newest label stays unique thereafter even when other children grow. For `|M|=1` the old set is empty. Define the joint post-promotion source measure

```
I_M(F) = π[ b_u Σ_(M embeddings born at the promotion by u) F(post-tree, embedding) ].
```

Here `b_u=L_u/D_u`; the chosen core leaf is exchangeable, and the new core has one leaf. No domain mask is applied at birth: an embedding might enter or exit the coefficient support later. `S_t^M` follows the **same tagged vertices** in the evolving physical environment without creating any replacement tags. Then

```
π_M(F) = ∫_0^∞ exp(−μt) I_M(S_t^M F) dt.                      (1)
```

For example, a sibling pair is born when the second selected sibling is promoted, with one choice per existing selected sibling. A three-leaf star is born with one choice per *unordered pair* of existing siblings. A bent four-vertex motif can arrive through either the terminal tip of its long arm or its short arm. Both types of births enter (1), with their actual host environment. Promoting one of the old motif vertices while following it is an ordinary physical transition, not a renewal of that embedding.

Proof: on a finite tree, sum over the unique birth times of embeddings, applying the Markov property at actual promotions and their predictable compensation. Equivalently, apply `π Q_t=e^(μt)π` to the owner score. The finite-age identity with multiplier `e^(−μt)` has the source integral up to `t`, plus initial embeddings or embeddings touching both old and new cells. There are only a bounded-radius set of old anchor vertices, and new vertices that can appear in a mixed embedding are at bounded distance of the old cell root. Their expected multiplicities grow at most polynomially in `t`, using independent rate-one promotion proposals to dominate each finite-depth tree of births. The old finite-radius moments under `π` are finite: move each configuration to its highest ancestor by mass transport, then use the physical stable-age law and the same finite-depth proposal domination for a newborn. At fixed depth, moments of any fixed order of the proposal count are polynomial in elapsed time (induct on depth using the Poisson moment formula). These bounds justify localization, make the mixed and initial terms vanish after discount, and permit the use of (1) for bounded tests and tests of polynomial local core counts. For ordered labels use the common automorphism convention in the birth source and in the occupation. No probability normalization for a multi-embedding count has been inserted in (1).

### Directed pressure and all-degree moments

For an embedding `E` and its role `i`, write `r_i(E)=Σ_(x adjacent to i, x outside E) 1/D_x`. With every motif census `x_i=(D_i,k_i)` retain the joint pressure measures

```
π_M[ H(x_E) ∏_(i∈E) r_i(E)^a_i ] .                         (2)
```

Expand each power into ordered donor slots, allowing repetitions. Group slots naming the same donor at a receiver; a block of size `j` contributes `D_donor^(−j)` and uses one new vertex. A donor cannot be adjacent to two distinct vertices of a connected subtree in a tree. Each term is consequently a connected augmented embedding, with at most `|M|+Σ a_i` distinct vertices, whose weight is an exact product of all-degree reciprocal moments. The groupings use labelled slots (or their exact multiplicities). If an external parent becomes the highest vertex, reverse the owner by mass transport, retaining the distinguished original motif and its roles. This is also a mass-transport identity for each measurable subset of the original census, not just for global means. For example the one-vertex first moment gives `π(r)=π((k+1)/D)=1−μ`; the restricted law is not a product of these means.

The physical generator of a test of the motif census alone has in gain row `i` the actual rate `L_i+Σ_(j inside adjacent)1/D_j+r_i`; promotion uses `b_i` and appends no tag. Thus true signed Neon forcing (including every `g`/higher transfer and every new-motif term) requires only the **zeroth and first** conditional pressure moments. One can instead evaluate it from uncharged physical motif drift: assign a newly created motif to its old connected set from (1). The charged and uncharged sums differ by exact opposite transfers; do not plan to bound their magnitudes separately. Evaluating the pressure via (2) does need up to **five** core vertices for a four-core row unless the one-donor distribution has already been integrated. Higher moments in (2) retain the correlations if a joint law is wanted: at fixed census the pressure is supported in `Π_i[0,(k_i+1)/2]`, so the polynomial moments determine its finite measure. This is an exact representation, not a closed finite moment system.

## 2. One normalized measure for the entire signed source

Turn a motif sum into a downward fringe test `F` with the same owner for all terms. If its cost requires the exterior parent of that owner, lift its score to that parent's fringe, with the original owner distinguished. This uses the same mass transport as (2). For the pinned drift, the first-body leaf index is capped; on other motifs the weight on each tail coordinate is cubic in `18/D`. A leaf gain changes that coordinate by `O(D⁻²)`, cancelling the possible factor `D` in its rate. Promotion has rate at most one. For the masks, a seam at a finite degree incurs only a finite rate. Exterior pressures at large `k` are bounded by `D/2` and multiply the same gain differences or transferred controls. Vertical motifs have one chain per bottom; fork/bend multiplicities are bounded on their supported branching coordinates. Thus the normalized forcing is a bounded-depth test with bounded per-incidence physical coefficients and can be allocated with bounded multiplicity per vertex, including after lifting one parent. In particular `|F|≤M W` for some finite `M` using the nineteen-mark weight: every edge on this bounded-depth assignment has a factor bounded away from zero. Alternatively retain the reward as a sum of the two nonnegative physical charged sides, used with their original signs. The constant `M` here is a norm bound, not the small stationary mean of the source.

Use the physical identity with the [positive first-moment residual process](../../../../../astra/source/round2/fringe/renewal.md); its root killing and copies remove the exterior parent rate by **exact** signed cancellation summed over owner vertices. Put `G_s=(s−R)^−1` on downward fringe scores. With newborn labels `2,…,20`, let `J` select the newborns, `B` be the direct newborn forcing column operator, and `K_s=JG_sB`. The direct label-two trigger keeps the old particle alive. Remove its newborn column from `K_s`, writing the rest as `C_s`. Define

```
ω_s(F) = e₂ᵀ (I−C_s)^−1 J G_s F,
Z_s    = ω_s(1),
ν(F)   = ω_μ(F)/Z_μ.                                       (3)
```

The normalization is exact: `e₂ᵀ(I−C_μ)^−1` is a positive left vector of `K_μ` by the scalar return value one; the [rank-one limit](../renewal/README.md) identifies the quotient with the fringe law. Formula (3) combines paths through consecutive marked spines while preserving the input to each copy. The label `min(D,20)` caps a permanent mark only; arrows and degree moments still use the actual `D`.

For any candidate number `c`, evaluate `ω_s(F−c)`, not unrelated intervals for numerator and denominator. More generally, for a full trial `U` with the appropriate weight bound, and `A=R+B_nontrigger`, a scalar exact defect identity is

```
ω_s(F−c) = U(J₂) + ω_s(F−c−(s−A)U).                      (4)
```

This follows from the subcritical resolvent, or discounted Dynkin with terminal term vanishing in the common weight. For finite retained states the corresponding stopped identity has a single signed exit measure: push forward and **sum coefficients at identical complete exit states before** pairing with the unknown common continuation. The absolute remainder is bounded by the occupancy of the absolute defect, not by separate continuation intervals on all paths. For a quotient trial certify `|ω_s(F−c)|≤ε Z_s`; then `|ω_s(F)/Z_s−c|≤ε`. This keeps shared normalization and the true sign of the forcing. When using an approximate trial for a five-vertex flux or donor closure, its physical defect must still be included in (4).

The whole uncharged drift satisfies the useful cross-check `ν(ℒΦ−λΦ)=(μ−λ)ν(Φ)` by the mean stationarity identity (integrate the finite-tree mean drift over a unit interval before sending time to infinity). Using this alone to evaluate the drift to `1e−8` is circular if the intended result is to infer `μ` at that scale.

### Discount increment and renewal-generation error

For `s<t` in the convergence domain, all cycles in (3) have a common positive measure on elapsed times and exit states. Write its kernel as `A_u` (first-moment semigroup). In particular

```
ω_s(F)−ω_t(F) = ∫ (exp(−su)−exp(−tu)) (A_u F)(J₂) du.
```

One can apply this to a previously centered source or score. If `|F|≤M W` and `r<s` is in the convergence domain, the absolute difference is at most `(t−s) M ω_r(W)/(s−r)` since `u exp(−(s−r)u)≤1/(s−r)`. For bounded `F`, the normalized quotient difference is at most `osc(F)(Z_s−Z_t)/Z_s` by the positive mixture of the common part and the extra mass. These are **increment** bounds; they give no absolute anchor for a missing part common to both discounts. The same discount-increment argument applies to a positive lineage class after applying one common support cut.

The [physical-trigger comparison](../../../../../astra/source/round2/finite/gauge/physical-trigger.md) supplies `‖C_s‖_g≤θ=.976434`, `min g≥m₀=.12510370`, `max g=1`, throughout the incoming interval. Write `ω_s^{(m)}` for the first `m` terms (`0,…,m−1`) of (3). Since `1≤W`, `(s−R)W≥7/10` for `s≥.51`, `JG_s1≤50/7` and `Z_s≥1/(s+2)` from the first virgin residence. Positivity gives

```
τ_m(1) := (ω_s−ω_s^(m))(1) ≤ (50/7) θ^m/[m₀(1−θ)],
|ω_s(F)/Z_s − ω_s^(m)(F)/ω_s^(m)(1)|
      ≤ osc(F) τ_m(1)/Z_s.                                 (5)
```

For weighted unbounded `F`, replace the last numerator by `τ_m(|F|)+|c_m|τ_m(1)`, `c_m=ω_s^(m)(F)/ω_s^(m)(1)`, with `JG_s W≤5/(s−.51)` using the original weight. In the bounded case, the exact rational thresholds for normalized error `≤1e−8` from (5) are `m=1139` for unit oscillation, `m=753` if an independently verified oscillation is `1e−4`. A small expected signed forcing does not show such a small oscillation. These cuts count only **untriggered direct generations**. They do not compute one residual Green row and cannot certify its internal truncation.

## 3. Why the intra-cycle row is still the obstruction

**A projection counterexample.** Take a virgin tracked root with a bare degree-two parent and give the grandparent degree either two or four by adjusting only the grandparent's leaves. The tracked state, its current pressure `1/2`, and the parent's entire census agree. A gain at the parent changes the pressure by `1/3−1/2`, and has intensity `1/2+1/D_grandparent`. The tracked promotion contributes `+1/4` to instantaneous pressure drift in both cases. The two conditional generators on pressure therefore differ by **`1/24`**. Both configurations satisfy the physical orientation. A projected deterministic transition on current pressure cannot represent both; allowing all conditional mixtures only makes a relaxation. This example does not imply a lower bound of `1/24` on any stationary error. Similarly, attaching a further donor census updates the generator via its own donors; reversal alone does not supply their conditional distribution.

**A count cut that ignores leaf volume, with an explicit cost.** Let `C` count cores in a residual fringe, and use the degree-dependent weight (`round2/finite/residual-weight.md`, historical reference) with `RW≤.448W−h`, not a bound in ordinary tree size. Copies do not increase `C` even after a leaf shift; only physical promotions increase it by one. In the same degree table, take the maximum of `b[Δ_c h+z(mark)]/h` alone. The replay checks every inner census and the capped constant tail: the maximum is `14283051730/35661533271 < .401`. Thus the sum of promotion-successor weights times their rates is at most `(C+.401)W`. Other gains are not included in this estimate. For `V=C^δW`, concavity gives

```
RV/C^δ ≤ .448 W − h + δ(1+.401/C)W.
```

Use `W/C≤h_max/(1/20+.18)` and `.401δ h_max/(1/20+.18)≤h_min` with `δ=.117`. Thus `RV≤.565 V≤LO·V`. Stop each residual lineage at its first `C=N+1` and apply the usual positive stopped account at any `s≥LO`. The sum of discounted exit `W` is at most `(N+1)^−δ W(start)` from a newborn. For a source with `|F|≤MW` the absolute continuation lost at such exits is bounded at a newborn by `M W(start)(N+1)^−δ/(s−.448)`. This bound concerns genuine whole continuations, applies without assuming independence of exits, and is unaffected by arbitrarily many leaves. Insert the analogous positive error for `B g` when truncating renewal rows; it cannot be omitted when using (5). With unit prefactor the core-count envelope would already require a cutoff above `10^68` to pay `1e−8`. This is a limitation of the proved envelope, not a lower bound on the true error or necessary cost of a sharper method.

For an explicit **normalized** error including both row defects, let `C_N≤C_s` and `f_N=JG_{s,N}(F−c)` be the core-killed rows, where `c` is the *single* candidate quotient. Let `f=JG_s(F−c)`, `R_C=(I−C_s)^−1`, and `R_N=(I−C_N)^−1` in label space. Then

```
ω_s(F−c) − e₂ᵀ R_N f_N
  = e₂ᵀ R_C [ f−f_N + (C_s−C_N) R_N f_N ].
```

Write `‖·‖_g` for the weighted max norm on columns (or its induced norm); a nonnegative matrix has norm `‖Mg‖_g`. If `e_f` bounds `‖f−f_N‖_g` and `e_K` bounds `‖(C_s−C_N)g‖_g`, an outward radius around the candidate `c` is

```
(s+2) [ |e₂ᵀ R_N f_N| + g₂/(1−θ) · (e_f + e_K ‖R_N f_N‖_g) ].
```

This uses only `Z_s≥1/(s+2)`; a better positive lower on `Z_s` can replace it. For the generic cut one can take `e_f` from the whole-exit weight bound for `F−c`, and `e_K` from that same bound for `B_nontrigger g`. Alternatively certify a tighter bound on the combined signed column in brackets. This formula makes explicit why high precision in the retained row solve alone, or just iterating more spine cycles, does not suffice.

**High-degree mass cannot just be discarded as small forcing.** Even the true first-body term for a childless nonfounder at saturated leaf stock does not decay to zero. The paired sparse replay evaluates its two actual tail limits and an `A/D+B/D²` enclosure around each limit, retaining the live parent reciprocal. Limits are `7.191102599034e−5` (lower) and `1.6411401538e−5` (upper) in signed units. At the first all-stock point `D=14002`, the signed rows are about `1.7248e−5` and `7.1094e−5`. The proved uniform errors around the limits at `D≥10^6` are below `1.530e−6` per occurrence. To refine this, retain the signed **limit term times its actual tail mass**, and the next reciprocal coefficients times their joint moments, rather than replacing each exit with zero or paying separate row intervals. In this example the exact first-body expression is a polynomial in `1/D` with one coefficient affine in the live parent input, so **including** the corresponding mixed moments removes its stock truncation error identically. Moment transport places those mixed terms on directed edges. This example is not a stationary positive lower bound on the omitted `k=0` population; it tests the row error, not its occupation.

## Practical criterion

The useful next certificate is a two-sided value of the **one** projected signed defect in (4), with tail rows expressed in joint reciprocal moments and an explicit bound on the residual flux from hidden donors. An absolute anchor for this value is required; adjacent-discount differences and exact normalization cannot recover it alone. A rough unconditional cap on the two nonnegative residual sides is not a substitute for that estimate. This route replaces absolute tree size by local ancestry and coupled pressure, but neither the projection nor the present weighted core-count estimate supplies an efficient certified `1e−8` evaluator.

## Replay

From this directory, using the preverified Neon payloads:

```sh
timeout 20s taskset -c 7 nice -n 10 env OPENBLAS_NUM_THREADS=1 \
  ./round2/.venv/bin/python -Bu check.py \
  b3f18_lo.bin b3f18_hi.bin
```

The rational checks verify unique newest-vertex counting against connected subsets, the renewal-generation bounds, the degree-weight table with the promotion-only maximum, the pressure nonclosure example and the sparse actual signed `f` tail. No raw files are copied, rehashed or exhaustively scanned. The newborn-label contraction is inherited.
