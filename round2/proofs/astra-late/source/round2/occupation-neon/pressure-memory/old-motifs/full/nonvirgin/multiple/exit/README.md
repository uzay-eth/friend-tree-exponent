> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../../../README.md#claim-21) for exact scope and citations.

# Source-normalized old-copy exits

**The complete exactly-two-copy coefficient is below `.916` at both existing trials**, compared with the previous upper of about `3290.802`. This is a joint bound on all exits using a source-normalized positive renewal comparison. It does not give a new exponent endpoint or a tight estimate of the second coefficient.

## Harmonic and Cauchy accounts at the unchanged boundary

Fix a trial `s`. Let `U=(s−R̂)⁻¹h` on the *full* physical domain; `U(J₂)=F(s)>0`. In a first-event or uniformized walk, retain the original degree/core and step exit thresholds. Include the copy tally on the exit transition. For the resulting nonnegative discounted path kernel `P`, the harmonic transform uses `P^U(x,y)=P(x,y)U(y)/U(x)` and terminal trigger probability `h̃(x)/U(x)`, where `h̃` is the trigger's one-step discounted reward. Their sum is one by the renewal equation. This is a lineage measure, not a transform of an individual core's physical clock. Multiplication telescopes along the marked spine, preserving the old threshold and giving the exact exit probability contribution `U(y)/U(J₂)` times the original discounted path weight. At cutoffs, use the actual crossed subtree and mark, not the last interior state or a reset virgin.

Let `ν_j` be the original discounted first-exit measure at copy tally `j`, at the geometry seam or at the retained step boundary. Write `f_r(y)` for the full continuation coefficient of copy count `r` starting at `y`. The omitted two-copy fraction under the transform is

```
Σ_{j=0}² ∫ [U(y)/F(s)] [f_(2−j)(y)/U(y)] ν_j(dy).
```

Thus for each `j`, Cauchy bounds this integral above by

```
sqrt( ∫ (U(y)/F(s))² dν_j · ∫ (f_(2−j)(y)/U(y))² dν_j ).
```

Or, dropping the filter, it is at most `sqrt(ν_j(1) ∫ (U/F(s))² dν_j)`. These are inequalities under the *actual marked first-exit measure*, with no independent redistribution of arms or emission ports. Using only `f_r≤U` and a uniform weight bound for `U` in this Cauchy estimate inherits the bad geometry factor. In particular, a global first moment of the harmonic likelihood ratio does not bound its second moment: mass `ε` carrying ratio `1/ε` keeps the first moment at one with second moment `1/ε`. This is an obstruction to that inference, not a claim that the physical exit moment diverges.

The useful comparison here is instead to use **one global first-moment budget** with the source exactly matched. Under the transformed law, complete zero/one-copy paths and retained two-copy paths are disjoint. Their measured positive fluxes give a joint maximum for what remains, including every geometry exit, every step exit, and every copy count above two. This is why no per-exit square or generic `W` multiple is needed for the numerical bound. The exact threshold and two-copy lower subtotal are unchanged. The result does not charge the entire budget to each individual exit.

## Explicit frozen Perron source budget

Use `a=.56564`, `R=1.006`, `θ=.976433602362`. The [nineteen-mark trigger comparison](../../../../../../../../../../astra/source/round2/finite/gauge/physical-trigger.md) gives on the **same** positive Perron vector `g` of `K_a`:

```
(C_a + k_a e_2^T)g = K_a g ≤ R g,
C_a g ≤ θ g.
```

Multiply `k_a g₂≤(I−C_a)g+(R−1)g` by the positive untriggered inverse, then extract the `J₂` row and divide by `g₂`. This gives

```
F(a) ≤ 1 + (R−1)/(1−θ).
```

There is no small-`g₂` penalty. The `g₂` at terminal degree-two renewal and the initial `J₂` cancel exactly. The inherited bounds use a *prior* exponent bracket to certify the Perron comparison; no equality `F(s)=1` is used at either tested `s`. The replay checks the archived rational lower `g₂`, the rank-one contraction arithmetic, and the rational log-convex bound on `R`. It does not recompute the entire older Perron-vector proof.

Before any production from `J₂`, the first holding denominator is `s+2`; the remaining positive numerator is decreasing in `s`. Therefore at both trials

```
F(s) ≤ M_s := [1 + (R−1)/(1−θ)] (a+2)/(s+2).
f₂(s) ≤ M_s − L₀₁(s).
```

Here `L₀₁` is the previous all-degree lower of the **complete** zero/one-copy class on the same chronological generator. Subtracting also the retained exact two-copy lower `L₂` gives the budget for the unresolved two-copy mass, or for the strictly disjoint union of unresolved two-copy and ≥3-copy paths. Normalizing at the actual source yields the upper `1−(L₀₁+L₂)/M_s` on the corresponding harmonic fraction.

| Trial | `F(s)` upper | Complete `f₂(s)` interval, outward | Unresolved two-copy upper | Normalized unresolved upper |
|---|---:|---:|---:|---:|
| `.565676636772` | `1.2545818796` | `[.0008341635188,.9157613370]` | `.9149271735` | `.7292686021` |
| `.565771871520` | `1.2545353128` | `[.0008339351325,.9157808984]` | `.9149469632` | `.7293114462` |

The scalar comparison alone cannot reduce the upper `M_s−L₀₁` appreciably further: nonnegative coefficient accounts can allocate almost all the slack either to `f₂` or to ≥3 copies with identical source total and zero/one bound. The exact fixture constructs both accounts; they are not candidate physical laws. A sharper geometry-dependent result would have to distinguish the two using the post-exit *remaining copy count* or certified moments under `ν_j`. Neither this bound nor the current lower subtotals determine a trial sign.

## Replay

```sh
taskset -c 8 nice -n 17 timeout 8s python3 -Bu round2/occupation-neon/pressure-memory/old-motifs/full/nonvirgin/multiple/exit/replay.py
```

Passed rational harmonic conservation, strict Cauchy and scalar-complement fixtures, and both trial budget checks. The input lower bounds are the outward results from the adjoining exact two-copy and all-degree one-copy replays, not new runs of those state propagators.
