> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Original Gold2 upper: joint 3999→4000 entry and edge reserves

**Global result:** `μ ≤ 0.565735923284044`, using the original upper witness and the existing uniform D24+H3 composition, with new Gold2-only edge reserves. This improves its best previous upper endpoint `0.565736253459` by `3.3017495e−7`. In particular `μ ≠ 2087/3689`; `241/426` remains inside the enclosure. The Gold2 solver, witness and existing receipts are unchanged.

**Seam result:** physical stationarity eliminates the whole selected old-edge generator, including exterior gain increments. Pricing the surviving signed entry jointly cuts its upper bound by **53.8%** versus the global far-edge cap times a joint degree/arity crossing moment. It does **not** produce an additive upper-witness credit; the endpoint improvement above comes from the edge reserves below.

## Exact stationary cancellation

Let `U` sum the original upper `g` over physical unoriented edges with an endpoint `D≥4000`. No arity restriction. Let `E_s` be the sum of each selected old edge's full true drift less `s` times its old weight. It includes both sends on the edge, all exterior sends and all old-mark promotion changes. Let `B` be `g(post,J)` on promotions at `D≥4000`, and `I` the **new** selected value of each old edge entering at a first gain `3999→4000`. Then

```
π E_s = (μ−s)π U − π B − π I.
```

Pathwise `L U = E_s + s U + B + I`; the sparse read checks that the whole double-G2-tail face is zero, and checks literal ownership on ordinary-tree arrows. The identity holds with **no** exterior-independence assumption. To justify the stationary step, take finite-tree expectations, integrate with weight `exp(−μt)` and divide by elapsed time. The boundary term vanishes since `|U|` is at most the far-edge cap times core count. At a gain above the seam the changing high-coordinate weight has an exact `1/[D(D+1)]` factor; with at most `c≤D` affected incidences and total rate `≤D` this gives a uniformly bounded per-core high-end gain contribution. High promotions change an edge by `O(1/D)` via the clamped class, and their total rate is bounded by one per core. The other nonzero actions at the low edge endpoint occur only at `D<252` (including entry to 252), so those, `B`, and `I` are bounded per-core local tests after allocation at their active core. The spatial mean limit therefore applies. A fixed-trial comparison must retain `(μ−s)π U`; it is not zero unless `s=μ`.

The entry is not `g(4000,neighbor)−g(3999,neighbor)`. The raw change of the previously unselected edge is still in the low sector. Also, the upper seam face has both signs: its minimum is `−175539908336.917` ticks, maximum `21976.55` ticks. Thus a proof that assumes every entry is negative would be invalid. The bounded dominant magnitude satisfies the exact outward comparisons

```
155105227 ≤ −π I ≤ 3288652241 ticks/core
global-cap comparator ≥ 7124998973 ticks/core.
```

The comparator uses the original absolute far-edge cap, not a fitted average. Both estimates use the same first crossing and full uncapped support. The smaller bound is on a **stationary flux**, not a reduced production cap on every tree. In the raw upper split the term is `−π I>0`, not a rebate; replacing the old-edge subtotal would additionally require independently pricing the accompanying low-sector complement, including the old low-low seam jump. In particular neither the flux magnitude nor the difference between the two flux bounds can be subtracted from the original 463-tick edge cap.

## Live-neighbor calculation

At a nonfounder host's physical birth start at `(2,1)`. Let `τ` be its first crossing of 4000. Discount to `τ` and pay all its actual neighbors' new edge weights. Multiplying the resulting expected payment by `μ` gives the stationary flux: `π B1=μ`, and the rate of the crossing is the true `3999−c+Σ 1/D_neighbor`, including the marked edge. This is compensation against the physical stable-age mixture; it is not a product of degree marginals.

For `(-g)_+` the replay prices each distinct child from its actual registration, retaining both live degrees and the child's arity through degree 32. The two tracked gains include reciprocals in **both** physical directions. Exterior sums and host arity are relaxed at checked rate corners. On an upper child exit above 32, carry the far-arm maximum only through the remaining host gains (at rate at most its degree); on a lower child exit pay zero. The lower sweep may kill the registered child contribution on a later host promotion, with hazard at most one. Host promotions are counted by the outer solve at their actual degree; the host arity states through 44 retain post-promotion continuation. The final registered child is paid at old arity 44, and later registrations are bounded separately. At a common independent exponential age `A` with parameter `μ`, `P(c(A)>k)≤(1+LO)^(-k)` by the rate-one promotion proposal bound. As arity is monotone, the same age bounds the discounted arity at `τ` above the cutoff. The unique mother is allowed any mark and costs at most the seam maximum times the crossing mass. The positive part of `g` is separately charged by the maximum positive word on the seam times the crossing arity moment. All rounding is directed and the corner inequalities are checked using integers.

## Credits that do compose globally

The charged upper `g` rows give four disjoint payments. They use different allocations from the raw seam calculation.

| Payment | Original row slot | Floor (ticks/core) |
|---|---|---:|
| D18 | Entire edge cap minus charged row, both endpoints `D≤18` | `238.0808742657006` |
| Small payer FY | Exact-G2 non-M3 host `19..251,c≤17`; payer `D≤5`; distinct donor, absorbed at `12` | `440547.471096166` |
| Large payer FY | Exact-G2 non-M3 host `D≤150,c≤10`; payer `6..34,c≤7`; distinct donor `2..5`; paying-edge max `>18` | `8313.1113818097` |
| Double tail | Entire 463-tick cap on the subset of edges whose parent was already `D≥252` at birth, once the child also reaches 252 | `0.3606267611` |

D18 retains the true centered singleton transfer and the original H3 newborn/transfer allocation. The two FY sweeps use the unchanged 20-knot debit formula with **upper** outward sign, the original upper payload and the original Gold2 exponent bracket. The knot gap is nonnegative for both signs, and the edge cap still covers the resulting augmented raw row by the original convex knot bound. The underlying registration machinery is unchanged; only in-process parameters and the selected witness differ. The original FY masks and argument are in the [cold audit](../../../../../../reviews/cold/astra/checks/gold2-h3-sol.md).

The double-tail `g` row is zero: it has no own value or successor, no exact endpoint singleton transfer/FY, no M3 transfer and no newborn H3. Degree is monotone. The replay lower-bounds the flux of promotions at `D≥252` by `π(D≥252)−π(c;D≥252)/252`. It then applies a uniform lower first-hit transform from each real child birth to `D=252`, and divides by `HI`, since the edge persists and its expected discounted residence after that hit is `1/μ`. This counts only the retained registration subset, never adds an age-independent product assumption, and even this tiny credit is a separate valid edge-cap reduction.

The existing D24 upper profile used zero edge credit. Its **uniform** trial (not its restricted contradiction margin) implies `(u_uniform−μ)X≥0`, with `0<X≤X_D24` for the original upper potential in ticks. More explicitly, the same trial ledger has its original constant repaired margin plus all four fresh edge credits on the right; the credits preserve the unchanged local H3 and singleton cap allocations. Thus the new upper is `u_uniform−sum(credits)/X_D24`, rounded **upward** to 15 decimals. There remain about `11,185,763` ticks/core to reach `241/426` with this same denominator conversion. For the rational exclusion, `1945/3438 < new upper < 2087/3689 < old upper < 142/251`. The outer two fractions have determinant one and denominator sum 3689, so the new cut excludes exactly one rational through order 3689. No improvement in the minimum possible denominator or in an algebraic degree/height bound is claimed.

## Replay

```sh
nice -n 19 python3 -IBu round2/occupation-quintet-gold2/true-drift/check.py
```

The replay verified the original source HEAD, receipt hashes and full upper payload SHA-256, then read 132,756 sparse words. It checked 60,130 literal arrows, including 12,015 nonzero exterior changes, 14,535 D18 edge classes, and FY solves of 3,585,120 and 3,946,320 states. Last lines:

```
new upper (upward 15 decimals): 0.565735923284044 gain vs best prior: 3.301749568757967e-07 gap to 241/426: 8.223753526692278e-06
more ticks needed from new upper on same denominator: 11185763.463340959
new rational exclusion: 2087/3689 through order 3689 with 241/426 still inside
PASS exact upper-only seam corners and edge composition; sparse words: 132756
```
