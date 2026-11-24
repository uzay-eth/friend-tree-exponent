# Independent check: current marks give a positive physical identity, not a lumping of the old residual

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 13: Claude for transfer bookkeeping; GPT for the uncapped bound and root transfer; inputs include physical jump geometry and the Round 1 exponent and prior. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The current-mark coefficients in `eval/gauge.py` agree with the physical generator. A useful exact simplification is **one forcing function of unmarked shape, followed by several intrinsic Green applications**. Killing labels still enter the Green operators. The original permanent-mark residual is not strongly lumpable onto unmarked subtrees, even if the root killing label is retained; its newborn kernel also differs from the current-mark one.

## Coefficients and the simultaneous update

Let `q_D=1/min(D,20)`, `p_x=1/D_x`, and let `F(S,a)` be a test on unmarked planted shapes with an exterior root label. On a physical tree sum `F(S_x,q_{D_parent})` over nonfounder cores. For an old core `x`, the intrinsic shape generator holds its exterior label fixed and omits its parent's incoming arrow. Put `r_x⁰=L_x+Σ_y p_y` over old core children `y`. In addition to its intrinsic shape row, the positive active row at `x` has pure death at its exterior label `a` and emits:

| Successor or score | Rate |
|---|---:|
| `J` active at `q_{D_x}` | `L_x/D_x` (direct) |
| `S_y⁺` active at `q_{D_x}` | `p_x` per old child `y` |
| `S_y` active at `q_{D_x}` | `q_{D_x}−p_x` per old child `y` |
| excess `F(S_y,q_{D_x+1})−F(S_y,q_{D_x})` | `r_x⁰` per old child, if the mark changes |
| excess `F(S_z,q_{D_y+1})−F(S_z,q_{D_y})` | `p_x` per old grandchild `z` through `y`, if the mark changes |

The founder uses this table with its **true slot-free** degree and no active killing term; its emissions supply immigration. The child death and the two child-copy terms combine to give the omitted signed incoming shape increment. In one `x→y` arrow the child jumps to `S_y⁺` at unchanged exterior label `q_{D_x}`, its old children change labels at unchanged shapes, and `x`'s own degree is unchanged. All three effects belong to the same arrow. When a child `y` sends back to `x`, `x` gains degree, `y`'s shape does not change, and `y`'s exterior label drops. Charging `p_parent` once in the grandparent row makes `r_x⁰` independent of the exterior, so the current particle needs no variable exterior killing.

## Dormancy, divided differences, and shape-only forcing

For `x<y` a dormant at shape `S` runs `L⁰`, dies at rate `x`, and activates at rate `y−x` into active label `x` with the then-current shape. Thus its holding diagonal uses **`y`**, not `x`. For any common nonnegative forcing `G(S)` and intrinsic Green `V_a=(s+a−L⁰)⁻¹`,

```math
 V_xG-V_yG=(y-x)V_y V_xG.
```

Both resolvents act on the same shape domain. Intrinsic boundary conditions must be identical for both labels (e.g. a common killed domain; otherwise a boundary difference has to be included). Consequently the numerical active responses satisfy `U_a=V_a G`; `G` includes the direct reward, the two active-copy actions and the two dormant-excess actions from the table. Dividing the dormant value by `y−x` and multiplying its emission coefficient by `y−x` gives the same model with no subtractive arithmetic. For the 19-label cap only adjacent mark pairs occur. One can solve on a **shape-only forcing** variable, applying `V_y V_x` for each such pair; no explicit permanent descendant labels or dormant shape storage are logically required. This is an exact Green-equation reduction of the **current-mark** construction, not an assertion that all Green queries collapse to one Markov type.

### Positive full-mean semigroup identity

The scaled dormant model also gives a short uniform **finite-time** bound. Realize its active emission as a birth of an active particle at rate one, retaining the dormant source, and kill the dormant at rate `y`. This implements `L⁰−y` plus active forcing of coefficient one; a rate-one *switch* plus death `y` would be wrong. Write `n` for physical shape size. In an active row the two ordinary copy terms add at most `n/2+1`; scaled child-drop terms add at most `n/3` because `r_x⁰≤D_x` and `D_x(q_D−q_{D+1})≤1/3`; scaled grandchild drops add at most `n/12` from `p_x≤1/2`, `q_{D_y}−q_{D_y+1}≤1/6`, and disjoint grandchild subtrees. Intrinsic grows `n` at rate at most `n`, the direct newborns add size at rate at most two, and deaths reduce it. Thus the generator on total particle size is at most `4` times that size. Dormant rows obey the same bound. Total event rate is linear in total size, so localization proves nonexplosion and finite-time expected size at most `e^(4t)n`.

To obtain the physical mean identity for an active observable `|f|≤n` independent of the exterior label, take the test to be zero on dormant types and use the full positive semigroup `P_t`. Its active rows at two exterior labels have the same shape law and forcing apart from constant killing. Subtracting their backward equations, with zero initial difference, shows that the scaled dormant row of `P_t f` is exactly `[P_t f(S,x)−P_t f(S,y)]/(y−x)`. This also follows from convolution with the constant kill survival factors. On a bounded horizon, `P_t f=O(n)` and its pointwise backward derivative is `O(n²)` by the first-event equation and the rate/size bounds. Apply the time-dependent physical jump formula to the sum of active scores `P_{t-u}f` over nonfounder cores. The coefficient identity cancels every interior term; only the founder immigration table survives. Yule moments on bounded intervals dominate the localization error (the sum over all cores of the derivative and jump-rate bounds is polynomial in physical size). This proves the full physical mean identity without a spectral gap or finite-domain approximation. It does not prove a strict residual growth bound near the target discount.

For a label-dependent observable there is an initial divided difference; initialize each dormant observable to that difference (which need not be nonnegative) to obtain the analogous identity. In particular the core count uses the simpler zero dormant observable.

### A pointwise condition behind the candidate weight

For the current path weight with edge factor `z_D=q_D/(q_D+d)` and root prefactor `1/(a+d)`, the descendant loss when a vertex's degree jumps from `D` to `D+1` is `d k_D` times the unprefactored descendant sum, with `k_D=(z_D−z_{D+1})/d`. A tempting scaled dormant value is the divided difference of these root prefactors times the same unprefactored shape score. Its descendant compensation restores exactly that intrinsic loss when both are charged at the same event. But its own drift is then `(q_H+d)` times the dormant value when the unprefactored shape score grows at rate `q_H`: the factor `+d` follows by writing `1/(x+d)` in units of `1/[(x+d)(y+d)]` and subtracting the dormant diagonal `y`. An active-row certificate has to check this dormant row too, and the edge compensation for external sends must be applied at their actual event rates. See the [cap-36 drift certificate](../weight/README.md) for the full inequalities.

## Why this is not a permanent-mark lumping

Two reachable states have one root leaf and one virgin child, and identical root mark. In the first history the root promotes at degree 2, then gains a leaf; in the second it gains a leaf first and promotes at degree 3. The child permanently has mark `1/2` in the first and `1/3` in the second. Now `p_root=1/3` in both. On the projected unmarked test `1_{shape=J}` the old residual has coefficient `1/6` in the first state and zero in the second; there is no other transition into that class. This is a generator-level obstruction to strong lumping even with an exterior label for the root. A positive permanent excess split cannot discard this datum: the same fixture has different emission intensity of the excess itself. The 19 pure-root survival curves `exp(−a t)` are linearly independent (evaluate at `0,t₀,…,18t₀` to get a Vandermonde). Thus a root-only common forcing does not imply a single scalar time-homogeneous survival phase representing all 19 arbitrary root starts.

There is also a newborn-level distinction: the separate [`physical/check.py`](../physical/check.py) derivative replay reports that the current-minus-permanent **unit-reward** residual offspring intensity derivatives at the mark-`1/2` virgin differ at order 4 by `1/96`, with orders 0–3 equal. Finite-time derivatives here are generator iterates on finite reachable rows; their difference rules out identical Laplace responses as functions of discount. A common physical full mean identity or equality of the critical exponent is not equality of the two residual kernels. The obstruction does not rule out a non-Markovian re-expression with history, a gauge-change identity for selected observables, or an alternative positive scalar renewal reduction preserving just the critical point.

## Replay and scope

Run `timeout 14s taskset -c 8 python3 -B round2/current-mark/check/check.py` from the repository root. The independent oracle enumerates all directed physical arrows on flat trees, materializes their successor, and matches the signed coefficient of **each** `(shape,label)` against both the independently coded formula and `gauge.py`. It checked 122 physical rows (ordinary small trees plus cap seams and simultaneous shape/mark fixtures), and 648 exact finite-domain dormant Green identities. It rejects substituting the changed shape into the same-time mark-drop term, charging the exterior twice, omitting the grandparent transfer, dropping below the cap, and using the lower-rate dormant diagonal. The hidden-child fixture gives exactly `1/6`; the Vandermonde fixture has rank 19. No kernel endpoint, subcritical drift bound for the enlarged residual, or new exponent estimate is certified by this replay.
