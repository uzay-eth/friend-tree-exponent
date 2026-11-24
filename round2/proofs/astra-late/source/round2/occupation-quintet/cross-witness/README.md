> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Original U5: paired old-edge seam renewal

**Result.** On a five-core pendant fringe, the normalized sum of the two original raw actions changes sign at the degree-2000 join for all six tested entrance stocks (path, bent fork and four-arm star, each at two stocks). The conditional physical occupation calculation pins the signed reward in degrees 2000–2003 to narrow intervals. Degrees 2004–2011 alone have a **larger absolute integrated reward than the paid window in all six cases**. The fixed-stock raw score turns positive again far out, so one cannot sign the whole unpaid tail from the near-seam result. This is a conditional renewal calculation and complement obstruction, with no new endpoint.

## Reward and physical slice

Use the original detached `b04cc50c1e391e0420eb7b006ac781d484c68543` source, `gold.json` and the two accepted payload copies under `round2/.work/replay-inputs/`. The checker uses sparse reads through the prior independent wire reader, pins the commit, payload header and golden receipt, and checks file identities before/after. Authentication of the payload bytes is inherited from the full cold audit; this run does not hash the payloads again.

Cut at the birth edge of a nonfounder core `0`. At entry, the younger side has exactly five cores and is one of the labelled trees in `SHAPES`; the only exterior core edge is from `0` to its parent. Core degrees count that edge. Initially `D_0=1999`; the other four degrees are either their physical minima or those minima plus one. Bare vectors, including `D_0`, are `(1999,2,2,2,2)`, `(1999,2,3,2,2)` and `(1999,4,2,2,2)` for P5, B5 and S4. Hold the topology until the first subsequent promotion, or until one of the four downstream degrees exceeds its physical minimum plus one. Changing or aging the exterior environment does not stop the process. The two windows also stop when the focal degree exceeds their upper limit.

For witness `i`, `W_i` sums every potential occurrence wholly within the fringe, with actual physical marks at its boundary. Normalize by `f_i(2,1)=2^40` on both sides. `A` and `B` are the real-valued actions `(L−s)W_i` at the common trial `s=.565724`; no verifier orientation is inserted before addition. The generator in the *reward* includes true post-event weights at every successor, including the new core and all internal newborn motifs on a promotion. Even a promotion that kills the occupation path contributes its instantaneous rate-times-jump to the raw action at the pre-event state. The decomposition `quintets` includes all P5/B5/S4 contributions to the action, including newborns. It is a component of `A+B`, not additional credit. The parent edge and crossing motifs are outside `W_i` by definition.

## Conditional occupation and correlation

For fixed fringe degrees `d`, let `x=1/D_parent ∈ [0,1/2]`. Gain and promotion hazards are `g_i=D_i−c_i+Σ_(v internal neighbor i)1/D_v+1_(i=0)x` and `b_i=(D_i−c_i)/D_i`. The physical law uses the same `x` in A, B, the holding rate and the next-state term. If `V` is a candidate continuation, verify both directions of the inequality against

```
  reward(d,x) + Σ_i g_i(d,x) V(d+e_i)
  -----------------------------------
       δ + Σ_i (g_i(d,x)+b_i(d))
```

`V=0` on exit and `δ` ranges over the authenticated bracket. There is no frozen-parent or Markov assumption on the exterior process. After multiplication by the positive denominator, residuals are affine in `x,δ`; for absolute rewards only the zeros of A, B or A+B add breakpoints. The recursion is triangular; each exact rational corner is checked after outward dyadic rounding. This argument applies to the physical compensator even with time-varying `x`. The literal directed-arrow control independently enumerates ordinary-tree event rates, including an exterior parent of degree 2000, and agrees on 144 family/witness/fixture comparisons.

| Bare entrance | Paid `A+B`, degrees 2000–2003 | Omitted absolute reward, degrees 2004–2011 | Paid five-vertex component |
|---|---:|---:|---:|
| P5 | `[-2.096377, −2.090514]e−7` | `[4.128552, 4.140103]e−7` | `[-1.541940, −1.541555]e−5` |
| B5 | `[-2.102224, −2.096669]e−7` | `[4.134149, 4.145076]e−7` | `[-3.905925, −3.904949]e−6` |
| S4 | `[-2.023478, −2.016809]e−7` | `[3.970546, 3.983646]e−7` | `[-5.421762, −5.420412]e−6` |

The table uses outward-rounded endpoints. Exact fractions, results for all six stocks and separate A/B intervals are in `replay.log`. For bare P5 the paired absolute reward is about `2.093e−7`, versus `2.383e−7` for `|A|+|B|`, a real **12%** cancellation. Pairing barely narrows the signed Bellman interval compared with separately optimizing A and B, so that is not a material gain. For the higher P5 stock the quintet component changes from negative to positive, while the total stays negative; the lower levels cancel most of the five-vertex action at the bare stock. Thus neither dropping the quintets nor substituting a single quintet value independent of pre-seam stock is justified.

## Noncentred tail and unpaid complement

For either fixed downstream stock the whole tail on `D=D_0≥2000` is rational in `y=1/D`. Each affected motif has the focal vertex as a tip or exterior donor, except the singleton. Its weight is affine in `y` in the first case and independent of `y` in the second. The singleton is quadratic. On a focal gain, substitute `y/(1+y)`; the difference has a factor `y²` before multiplying by the focal hazard `y⁻¹−c+ρ`. On other gains the hazard is affine in `y`; promotion hazards are affine too, and post-promotion weights retain the same degree dependence (or become identically zero outside support). Therefore after multiplication by `(1+y)^2` the degree is at most five. Six exact values determine it and five independent degrees check it at each exterior-receive corner. The dependence on `x` is affine. Exact Bernstein hulls of `(P(y)−P(0)(1+y)^2)/y`, divided outward by `(1+y)^2`, give the noncentred pointwise linear enclosure `R∞+m_-/D ≤ R(D,x) ≤ R∞+m_+/D`.

For the old **internal** edge `0–1` alone, the bare P5 limit is about `+0.2834249090` and the slope hull lies in `[-1.280439, -1.278397]` (outward). For the whole bare P5 reward the limit is only `+1.091625e−5`, but the slope lies in `[-0.232648, -0.231877]`. For all six fixed stocks, the full numerator is negative on `[2004,2011]` and positive on `[100000,∞)`. The linear enclosure survives conditional expectation using the actual reciprocal moment, including a state- or holding-time-weighted law. It has a nonzero constant, and using only an unweighted mean degree or Jensen with the constant removed would not represent this reward.

The omitted-window bound is an obstruction in **absolute** reward for the same stopped cohort, not an assertion about the sign of the entire omitted expectation. To quantify actual entry mass, separately force the four promotions and the necessary minimal downstream gains, then force only focal gains to 1999. At each step the desired event's discounted incidence is at least its hazard at `x=0` divided by the total tracked rate at `x=1/2` plus `HI`. The product gives a uniform entry floor from a physical newborn, regardless of the exterior. Multiplying by `LO` and the omitted-window absolute value floor gives per-core floors under the physical stable-age birth identity: **`>8.18e−30` P5, `>2.86e−33` B5, `>7.85e−38` S4**. These very small guaranteed floors establish positive mass and show why this particular entrance proof cannot certify a useful global margin. There is no bound here on the full crossing-motif ledger, all other entrance stocks, or the signed tail after the stopped window.

## Replay

```sh
taskset -c 1 nice -n 19 timeout 90s /usr/bin/time \
  -f 'elapsed=%e user=%U sys=%S rss=%M exit=%x' \
  python3 -IBu round2/occupation-quintet/cross-witness/probe.py
```

`replay.log`: **3,456** exact receive/discount state-corner checks across both windows and three topologies, **144** literal controls, **1,210** sparse words per witness; exit **0**, wall **9.95 s**, peak RSS **38,856 KiB**. All new files are in this directory.
