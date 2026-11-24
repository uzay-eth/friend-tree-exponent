> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Gold2 D24 stationary interval with matched original H3

**Both fixed trials pass:**

```math
\boxed{0.565714675151 < \mu < 0.565736253459}.
```

This is the Gold2 D24 stationary composition, using corrected conditioned singleton/edge inputs and both signs of the same original H3 mask. It clears `99/175` by more than `3.89436e−7` at the lower endpoint. It is separate from the standalone no-LP scalar theorem.

| Profile | Lower (approximately) | Upper (approximately) |
|---|---:|---:|
| Same credits, original scalar denominator, no host savings | 0.565714418717201 | 0.565736412684620 |
| D24: separately priced denominator and host savings | 0.565714675049909 | 0.565736253540235 |
| D24: uniform-prior joint objective | 0.565714675149971 | 0.565736253460852 |
| D24: restricted-prior fixed-trial comparison (exact) | **0.565714675151** | **0.565736253459** |

The scalar control in the first row uses the **same credits** as this D24 profile (conditioned singleton on both witnesses, three lower edge components, and matched H3 on both), so it isolates the effect of the LP. Its exact lower margin at `99/175` is positive (about `240549.4973` ticks/core). It is not the separate smaller-input scalar release coordinate.

## Original middle cap retained

At a physical H3 row with original outward signed residual `q`, write its original nonnegative production cap as `K` (positive-part of the production cap if necessary). The paid amount is `r=[−q]_+` within the frozen nonfounder mother→middle→child mask, zero elsewhere. On that occurrence `q+r=max(q,0)≤K`. Thus **the unchanged local H3 cap at its middle continues to cover the row after rebate**. Do not subtract the local cap again or reinterpret the payment as `K−q`. In the D24 cells, the H3 part of the higher envelope retains `binom(c,2)` times the original family cap (`408` lower, `233` upper) at the eligible middle. The other families keep their original allocations, including the half P4 cap at each interior end. The published aggregate cap and the computed family-cap sum each upper-bound the host sum *after* the rebate, so their pointwise minimum still does. Consequently the H3 floor and `H−h` host savings can be added regardless of overlap of their middle-degree support.

The mask is mother `19≤D≤79,c≤8`, middle `D≤5,c≤4`, child `D≤18,c≤3`; all are actual live censuses. The raw H3 row includes its own original transfers. Singleton `f`, D18 `g`, payer-≤5 FY, and payer->5 FY spend distinct production row slots from H3. No FY is paid at this H3 middle (M3); the transferred `g→H3` gain is charged exactly once by the original decomposition. The composition uses one parent→middle→child H3 price per witness; no separate H3 search floor, reserve, or overlap correction enters. The exact H3 floors in ticks/core are `1084423030952180046412413909/1550425090937500000000` and `3304204122287698944634328103/6201700363750000000000`. The two source logs are checked for matched physical mask, original side and no reserve or ceiling; their normalized fractions are multiplied by `2^40` exactly once.

## Stationary algebra and exact fixed trials

Let `d=+1/−1` for lower/upper; `X=πΦ` is the mean potential in ticks/core, `λ` the original nominal, `H` the original global higher cap, `h` the D24 cell's retained higher cap, and `E` the original total per-core error. Write `B=F+G+R_H3−E` for the signed net of the certified stationary floors. Original drift and the pointwise cap argument above give

```
d (μ−λ) X ≥ B + π[H−h].
```

`conditioned24.json` supplies exact `F`, `G`, host cells, denominator bounds and `π[H−h]` bound. All D24 degrees above the exact cutoff use the original overflow envelope; no H3 host row is cut out or assigned another potential. Since `B` is positive on both signs, the scalar and separate D24 tests divide by the applicable **upper** potential bound. For inward shift `a=d(s−λ)≥0`, subtract `aX` and use the original pointwise unrecharged upper `u` to obtain

```
[d (μ−λ)−a] X ≥ B−a M_max + π[H−h + a(M_max−u)].
```

The sided D24 LP bounds the nonnegative cost on the right in occupation units `π/μ`. Candidate duals are rationally rounded, every original column is repriced, and the certified repair is subtracted using the lower prior bound on `μ`. Multiplication by that same positive prior bound gives ticks/core. With the whole original prior this yields margins about `136.02469`/`107.97563` ticks at the separate D24 trial points; a positive uniform margin divided by the certified D24 upper potential bound gives a further explicit inward step.

For the stronger fixed-trial test, assume the conclusion fails, **including equality**. Then on the lower test `μ∈[LO,s]`; on the upper test `μ∈[s,HI]`. Restrict only the LP discount range and the nonnegative objective normalization to this interval. Hold all original H3, edge and conditioned singleton floors fixed over the full original bracket. The exact repaired fixed-trial margins are positive (about **41.742864** lower and **9.071194** upper ticks/core). Since `X>0`, the left side would be nonpositive under the contrary assumption. This proves strict inequalities. The restricted margin is used only for contradiction at its trial; it is not divided into an additional endpoint step.

## Reproduce

The [exact receipt](receipt.json) pins corrected `conditioned24.json`, the published Gold2 receipt, both H3 logs and the physical LP implementation. This composition reads those receipts and rebuilds the D24 LP; it does not read or re-scan either payload. The input profile's existing files were left intact.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 5 nice -n 19 timeout 45s \
  /usr/bin/time -f 'wall=%e rss=%M exit=%x' round2/.work/venv/bin/python -Bu \
  round2/occupation-quintet-gold2/host/h3-joint/check.py
```

Output: `outward 0.565714675151000 0.565736253459000`; `PASS exact both-sign Gold2 H3+D24 composition`; `wall=6.36 rss=105824 exit=0`. It reproduced `receipt.json` byte-for-byte after the generating `--write` run (`wall=15.51 rss=106832 exit=0`). SHA-256: `ae6f5cde2ef710b4196c3c445a3786580e5ccc9dfd9206b2962f60628866c6ed`.
