> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../../../README.md#claim-21) for exact scope and citations.

# Remaining-copy transfer on the actual first exit

**The complete `f₂` interval narrows by 8.22×** at both trials relative to the [source-only budget](../exit/README.md). The new upper is `0.112124` at both trials, with no critical-equality subtraction or common source budget. This falls short of 10× and does not move an exponent endpoint.

## Coefficient filter and common discount

In the positive marked walk count old-child traversals only. Fix the earlier sector (at most three cores, all degrees at most five), and stop on the first geometry exit or after 96 uniformized steps, whichever comes first. For exit tally `j`, retain the actual exit geometry and marks, and let `ν_j` be its discounted exit measure. The exact remainder of `f₂` is `Σ_{j=0}² ν_j(f_(2−j))`; the same identity holds with two replaced by zero or one. An exit caused by a shifted copy has already spent one copy. All entrances through the retained walk are included, and all future entrances from the continuation are paid there without re-killing the lineage.

This is a coefficient-filtered log-`z` moment. For a single exit at tally `j`, the second log derivative of `z^j U(z,y)` is `z^j(j² U+2j ∂logz U+∂logz² U)`; its coefficient of `z²` is `4 f_(2−j)(y)`. Using all three remaining-copy coefficients directly avoids paying counts above two or applying an unfiltered second-moment bound. The replay checks this product rule and exact coefficient transfer on an asymmetric marked seam. Use a common lower discount `s*=.565676636772` for continuation envelopes; positivity makes them valid at both trials, while the exact pre-exit transfer uses its own trial discount.

## Positive envelope on the actual exited state

Let `R_n(i,D,k)` be the all-degree table for a detached old root of exterior mark index `i`, with its **existing** child accounts omitted but all future births included. Write `A_n(i,j,D,k)` for the corresponding attached-child continuation, at child mark `j` and original parent mark `i`, with no accrued child account at entry. `A_0` earns a first copy with destination score `R_0`; `A_1` earns a first copy with destination score `R_1` and an extra account for grandchildren born before emission. For both, use emission `a R_n(Y)+p[R_n(Y⁺)−R_n(Y)]` with `a=1/j` and the same actual parent port `p` in the attached gain control. Existing children at an exit tighten that control to `p≤min(1/j,1/D_parent_exit)`. The root tables group only the original parent death discount, using an upper representative; each child's permanent mark remains exact. The original root's virgin promotion launches its branch account at the exact half rate, separately from nontrigger direct resets. Old root degrees never decrease, so the virgin correction does not leak into other root states.

For any copied detached root of mark `m`, the all-degree zero-copy value is at most `V_0(20)/(s*+1/m)`: the code checks the maximum against all inner table rows and the tail. Thus an existing grandchild of mark `m` attached to a root with exterior mark `j` has a one-copy lifetime fee at most `G(j,m)=(1/m)·V_0(20)/[(s*+1/m)(s*+1/j)]`. It covers both destinations and all future descendant geometry. For `A_1`, an already-present grandchild contributes `a/(s*+1/i)·G(j,m)`; a future grandchild adds this fee at its actual promotion, not at the initial time. This compensator follows by integrating the entire future constant emission rate `a` after the promotion, as a positive upper bound. The true detached continuation shares the grandchild on both shifted/unshifted ports. Therefore at the exact exited tree `t` the upper coefficient is the root table plus a sum over its **actual** current children, with their mark, degree and arity, and their **actual** already-born grandchildren with their marks. No current subtree is reset to a virgin. A grandchild's deeper shape is relaxed only inside `G`.

For every high-degree crossing use saturated birth fees, `E∫e^(−λt)b≤1/λ−(1−b)/(λ+1)`, and the same bound inside the attached lifetime integral. A shifted detached root has `1−b⁺=(k+1)/(D+1)`; hence the integrated detached term uses `1/(γλ)−(k+1)/(D+1)/(γ+1)/(λ+1)`. The promotion-account part uses the unshifted parent census. These bounds include all degree exits, including promotion at a tail root, and retain old children separately. The newborn reset vector in each layer has a positive a posteriori remainder correction `max positive residual/(1−α)` using the direct-only row-sum bound. For 12 rounds these corrections were `1.76102e−4`, `9.32651e−5`, `9.44477e−5` for the three layers, applied before their responses are used upstream. This is not a cutoff at 12 direct generations.

## Exact replay and scale

Rounded outward:

| Trial | Complete exactly-two-copy interval | Geometry-exit contribution upper | Inside step remainder upper |
|---|---:|---:|---:|
| `.565676636772` | `[.0008341635188,.1121232998]` | `.1111724014` | `.0001167349` |
| `.565771871520` | `[.0008339351325,.1121058555]` | `.1111552587` | `.0001166617` |

The same replay encloses `f₀` above by `.3177526076` / `.3177284405` and `f₁` by `.1161309475` / `.1161130030`. It integrates the boundary products with integer numerators on 2656 retained states, after directed rounding of the analytic envelopes. Geometry exits dominate; merely increasing the exact step count does not close the remaining gap. A plausible next refinement is to replace the `G(j,m)` bound with a joint old-grandchild age/degree comparison. A separate positive bound on ≥3 copies is still needed for an absolute `F` sign. The prior source-only numbers are used solely to display the width ratio. No D18 or `hV` charge is used in the return estimate.

```sh
taskset -c 8 nice -n 17 timeout 100s python3 -Bu round2/occupation-neon/pressure-memory/old-motifs/full/nonvirgin/multiple/filtered/replay.py --rounds 12 --both
```

Both trials passed; width ratios were `8.22117238` and `8.22262221`. The exact marked-seam and log-derivative fixtures passed.
