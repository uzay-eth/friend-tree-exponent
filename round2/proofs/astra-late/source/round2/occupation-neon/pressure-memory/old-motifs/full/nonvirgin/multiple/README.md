> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../../README.md#claim-21) for exact scope and citations.

# An exponential copy-count tail for the genuine positive trigger

**Uniformly at and above the lower trial, the generating series in the number of old-child copies is analytic for `|z|<1.00132`. The tail from copy number `m` onward is at most `21982.190751·1.00132^(−m)`.** This proves an all-copy error bound, but the constants make it unusable at low copy counts. The exact two-copy pilot lifts the certified trial lower subtotal above `.3395` at both endpoints. Neither endpoint moves.

## Tilt and domination

Use the original nineteen-mark positive process; let `T` be its old-child copy operator, `L` the intrinsic generator, `a_o` the constant exterior death mark, and `B₊` the nontrigger direct newborn production. Define `R_z=L−a_o+zT` and `F_s(z)=[(s−R_z−B₊)⁻¹h](J₂)` with `h=1/2` on the virgin. For complex arguments the definition is the positive-coefficient chronological series `Σ f_k(s) z^k`; each term counts exactly `k` copy traversals, and direct resets and terminal trigger do not get a marker. Treat real `z≥1` first, and keep the same two physical copy destinations.

Use the [original residual weight](../../../../../../../../../astra/source/round2/fringe/renewal.md) with `d=.15`, `q=.36`, `ρ=.51`. Its exact identities are `TW=(a_o+d)W−h₀(c_o)`, `W(Y⁺)=W(Y)`, `L W≤qW`, and `(ρ−R₁)W≥h₀≥.7`. Fix `ζ≥1`, `Δ=(ζ−1)(.5+d)`. To dominate **all resolvent sources**, a pointwise inequality on `W` alone would not justify the operator assertion. Instead conjugate by `W`: `W⁻¹ R_z W` has nonnegative off-diagonal entries. The additional tilted off-diagonal rates sum to at most `Δ`, but this alone would not order their state dependence or the semigroups. The needed comparison is the following **chronological kernel interpolation** alternative.

For `R_z` the weight drift gives `(ρ+Δ−R_z)W≥.7`, and stopped Dynkin yields a bound `W/.7` on the residual integral of a unit rate at that discount. In particular all rows of the nontrigger newborn matrix `C_(ρ+Δ,z)` sum to at most `50/7`. Use the [inherited un-tilted comparison](../../../../../../../../../astra/source/round2/finite/gauge/physical-trigger.md) `C_(.565,1)v≤.988137v` with `max v=1`, `min v≥.126364920098` (its proof only requires the earlier prior `[.565,.5665]`). Let `p=1/n`. A chronological path for a direct nontrigger newborn has measure weighted by `exp(−s time) z^copies`; Hölder on that **same joint measure** including its terminal label gives

```
C_( (1−p)·.565+p(ρ+Δ₀), z₀^p ) v^(1−p)
    ≤ .988137^(1−p)(50/7)^p v^(1−p).
```

The replay checks rational trial interpolation, `z^n<z₀`, and the strict comparison `(.988137)^(n−1)(50/7)<θ^n`. Both monotonicities are in the positive direction. At the actual trial and checked `z`, `(s−R_z)W≥.7`, so the terminal vector `k_z` is bounded componentwise by `W(J_i)/(2·.7)≤25/7`. One may use this bound in the full positive reset telescope and get a slightly larger prefactor than the table below. To get the displayed `5/2` terminal prefactor, use the exact direct trigger rate `h≤ h₀/2` (only in the virgin, where `h₀=1`), hence `(s−R_z)W/2≥h` and `k_z(i)≤W(J_i)/2≤5/2`. This yields `F_s(z)≤ (5/2)/[.126364920098(1−θ)]`. Positivity bounds `Σ_{k≥m}f_k(s)≤F_s(z) z^(−m)` and gives absolute convergence on the closed disk and analyticity inside. This uses no equality at either trial or any subtraction from critical mass.

| `z` | `θ` | Upper prefactor | Per-copy decay | Copy cutoff for `10⁻⁵` | Cutoff for `10⁻⁸` |
|---:|---:|---:|---:|---:|---:|
| `1.00101` | `.994` | `3297.328613` | `.9989910191` | 19430 | 26273 |
| `1.00118` | `.997` | `6594.657226` | `.9988213908` | 17220 | 23077 |
| `1.00132` | `.9991` | `21982.190751` | `.9986817402` | 16307 | 21544 |

These are sufficient cutoffs for the stated bounds, not lower bounds on the true copy depth needed. The first cone gives the smallest total prefactor among the three, while the last gives the strongest long-tail decay.

## Exact propagation through two copies

In the finite sector `D≤5`, at most three cores, uniformize the original generator with 2656 reachable marked types. Factor the positive discounted one-step operator as `P₀+zP₁`. `P₀` contains intrinsic moves, nontrigger resets and uniformization waiting. `P₁` contains the two copies using the same source arm and parent degree; the shifted destination increments root leaf stock but retains the same entire child tuple and every old mark. A subsequent copy sees that retained geometry. Keep integer numerator arrays for coefficient indices 0, 1, 2; an update is `x_k ← P₀ x_k + P₁ x_(k−1) + 1_{k=0} h/(Λ+s)`. The nested fixture checks the exact two-copy coefficient against direct summation through both old-arm destinations.

At each geometry exit also propagate its entire `W` with the copy tally of the crossing, in separate first-exit arrays. For a uniformization-step cutoff, bound the remaining inside `W` by the unmarked retained operator's power. These two accounts are disjoint by first exit. To bound continuation of exactly `k` copies at an arbitrary crossing, use `W·cM^k`. Here `c=(13/14)/(LO+.05)` because the exact zero-copy newborn value is at most `(.25)/(1−.659)<1` and the intrinsic discount is at least `LO+.05`. For `σ=LO+.05−.36`, the intrinsic inverse of `TW` is at most `(.65)/σ·W`. Its newborn source is at most five times that coefficient, and summing all future direct-reset generations costs at most `1/(1−.659)`. Resolving this newborn source back into the old starting geometry gives `M=(.65)/σ·[1+5c/(1−.659)]`. This overpays every finite-sector exit; it never resets an exiting old subtree. The resulting two-copy upper was above `3639` on both trials. The all-copy tilt, after subtracting positive zero/one contributions from its weighted budget, is a sharper but still very loose two-copy upper.

| Trial | Finite `f₀` lower | Finite `f₁` lower | Finite `f₂` lower | Combined all-degree `f₀+f₁+f₂` interval |
|---|---:|---:|---:|---:|
| `.565676636772` | `.2337400866` | `.0078309533` | `.0008341635` | `[.3396547061,3290.8020]` |
| `.565771871520` | `.2337274076` | `.0078294053` | `.0008339351` | `[.3395883495,3290.8019]` |

The lower subtotal uses the disjoint previous all-degree zero/one lower plus the exact new two-copy subreturn. None of these trial constraints use `F(μ)=1`. In the lower-trial run the discounted first-exit `W` accounts at indices 0, 1, 2 are `.69520`, `.150997`, `.0154503`; the remaining inside step-tail `W` is `.00123631`. Thus the geometry seam, rather than integer accumulation, is the first bottleneck. An implementable next target is a local positive two-copy continuation at the exact exiting child geometry, and tighter `C_(s,z)` row comparisons on an optimized weight; longer finite powers alone do not price that seam. No trial sign follows from these enclosures.

## Replay

```sh
taskset -c 8 nice -n 17 timeout 35s python3 -Bu round2/occupation-neon/pressure-memory/old-motifs/full/nonvirgin/multiple/replay.py --degree 5 --cores 3 --steps 96
```

The final replay passed both trials, the rational cone comparisons and tail cutoff checks, and the exact shared-subtree and two-copy fixtures.
