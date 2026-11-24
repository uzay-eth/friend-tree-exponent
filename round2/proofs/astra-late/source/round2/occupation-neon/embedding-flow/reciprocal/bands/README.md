> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Reciprocal slab balances at the physical D17→18 seam

`BandModel` strengthens the unchanged `ChainModel` with exact physical balance tests distinguished by the slab at each role. The original global constraints, objective columns, chain and edge masks, projection and pressure transport all stay in place. No paid-cost benchmark or endpoint is included here.

## Physical rows

For `y=H/D`, slab 0 is `[0,H/18]` (the physical degrees `D≥18`), and slab 1 is `(H/18,1]` (`H≤D≤17`). Bernstein control nets use the *closed* interval at both ends, but the additional support constraints below exclude a spurious slab-1 atom at `D=18`. Every role of a root, edge, or three-chain can cross this seam. If the slab-1 nonfounder arity is at least 17, the whole cell is empty and all its mass and pressure coordinates are forced to zero. If arity is 16, every slab-1 mass atom has `D=17` and the remainder is forced to zero.

Fix a slab-1 receiver `i` and an index tuple `j` over the other coordinates. Set `K_j=∏_{h≠i}B_{N,j_h}(y_h)/W_h(y_h)`, using the correct physical slab for each factor. With `η=π/μ`, introduce nonnegative variables

```
alpha_j = ∫ 1{D_i=17} K_j dη_occ
Gamma_j = ∫ 1{D_i=17} K_j g_i dη_occ
g_i = (16−k_i) + ∑ internal sender 1/D_sender + exterior incoming field.
```

`Gamma_j` is a gain intensity, not a probability or a separate copy of the exterior field. The stored seam slice in the receiver moment `r,j` is `alpha_j × B_{N,r}(H/17)/W_i(H/17)`, with the Bernstein argument mapped to slab 1. This slice is dominated by the existing receiver-basis mass, summed over its discount lanes. Box `Gamma` on that slice by `(16−k_i)` plus exact internal reciprocal for finite roles and physical lower/upper reciprocals for each high slab; add the actual exterior-port counts times `1/2` for the upper rate. This is a relaxation on the slice, not an exact recovery of conditional exterior fields from whole-cell moments. Global dual-repair budgets use at most `1/μ` occurrence mass per topology (reverse to the terminal, which has one ancestor tuple) for **each receiver role**, and at most `20/μ` for its `D=17` gain flux. Optional cross-topology seam-flux sharing is not imposed.

For the remaining slice `R=M−alpha`, use the inequalities requested for the exact gap to `D≤16`. With `t` the normalized slab-1 coordinate and `z` the normalized location of `H/16`, impose `∫ (t−z)B_{N−1,r}(t) K_j/W_i dR≥0`. In degree `N`, the two coefficients are `−z(N−r)/N` at index `r` and `(1−z)(r+1)/N` at `r+1`. This both preserves physical support and rejects a `D=17` atom hidden in the remainder (without requiring arbitrary continuous LP measures to be discrete everywhere).

The signed old-gain difference is unchanged for within-tail gains. For a gain crossing from `D=17` to `18`, that difference would subtract the continuation in the departing slab. Add `(H/18)^e_i Γ_j` to the departing row and subtract it in the entered row, multiplying by the other roles' exact representations of `W_h y_h^e_h`. No `W_i` is inserted: the receiver coordinate was integrated out of `Gamma`. The finite `H−1→H` seam enters slab 1 at `y=1`. Promotions and virgin insertion leave the old role's degree and slab unchanged; they enter the exact post-promotion arity cell, and the newborn's exact degree-two cell. A projection exit is still an exit from the retained family. All rights-hand sides are zero except the singleton virgin source of one.

## Exact checks and matched comparison

`check.py` uses point atoms and literal sender reciprocals to compute the outgoing and incoming coefficients of a physical gain, existing promotion, and newborn insertion independently of the slab-row construction. It checks root/edge/chain coefficients before and after `17→18`, a two-high edge, `D=16`, the `H−1→H` seam, and the `k=16/17` boundary. Setting alpha to zero on a true `D=17` atom passes simple mass domination but fails the new discrete-gap constraint. All 27 global balance rows of the small model are recovered coefficient-for-coefficient by summing their local rows. Final check: `literal physical atoms 13 nonzero local coefficients 82 global rows recovered 27 planted missing seam rejected H4 atoms 3 H4 coefficients 17 k16 atom/k17 impossible checked`; CPU `2.78+0.08` s.

In `compare.py` the program builds **one** original model, prices it, appends the slab rows, and reprices the **same cost vector on the same support**. The objective is the incidence of eligible chains with a high apex `D≥18`, or the same mass weighted by the true apex reciprocal. Both endpoints use rational dual repair. All reported lower bounds are zero after intersecting with the original enclosure. `arity=1`, chain arities `(1,1,0)`, only apex wide, and `inner=H−1` throughout.

| Matched support | Objective | Original upper per core | Slab upper per core | Reduction |
|---|---|---:|---:|---:|
| H=3, order 1, 3 chain cells | Mass | `0.00579604308403` | `0.00469387674231` | 19.0% |
| Same | Reciprocal | `0.000322002393617` | `0.000205417992923` | 36.2% |
| H=3, order 3, 3 chain cells | Mass | `0.00383976259796` | `0.00273782123245` | 28.7% |
| Same | Reciprocal | `0.000213320144385` | `0.000123883935456` | 41.9% |
| H=4, order 1, 16 chain cells | Mass | `0.00897456323297` | `0.00724226006064` | 19.3% |

At H=3 order 3, the exact repaired slab upper bounds for mass and reciprocal were `1558163959415540234370380170399632520109577544183409/569125529800056016023425940000000000000000000000000000` and `682267723042302310442631959978489463647818822139772109/5507313926769182055855488136192000000000000000000000000000`. The model grew from 2,718 columns / 65 equalities / 1,348 inequalities to 2,876 / 254 / 2,849. Comparison CPU: `1.42+0.30` s (H=3 order 1), `17.32+0.21` s (H=3 order 3), `1.29+0.10` s (H=4 order 1). These are bounds on the stated chain support, not hV charge or exponent amounts. The optional hV CLI calls the unchanged parent cost reader, but the signed paid-cost replay belongs to the other lane and was not duplicated.

## Replay

From the repo root, set `R=round2/occupation-neon/embedding-flow/reciprocal/bands` and `PY=./round2/.venv/bin/python`. Runs used `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 11 nice -n 19 /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M'` and timeouts between 15 and 40 seconds. Use an available niced core for a fresh replay.

```sh
"$PY" -Bu "$R/check.py"
"$PY" -Bu "$R/compare.py"
"$PY" -Bu "$R/compare.py" --order 3
"$PY" -Bu "$R/compare.py" --threshold 4 --inner 3 --mode mass
```
