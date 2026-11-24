> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Streamed complementary P4 strata

The release sources and receipts in the parent directory are unchanged. These are experiments on the same positive, complementary D18 support; any replacement must displace the earlier payment. This experiment has not booked a replacement or changed an exponent endpoint.

## Architecture and exactness

The projected skeleton is an oriented path census (prefix length one through four) and a bit for whether the specified distinct donor is registered. Histories giving the same census and actual registered donor/mother states share a row, irrespective of which physical promotion order created the two tags. Registration still retains the previous cohort and emits each eligible new alternative on the *same* promotion clock. A donor is never simultaneously registered as the path child. We do not identify unequal rates or different payoff faces on the basis of approximate values.

The reverse sweep stores a face only if at least one of the two authenticated witnesses has nonzero outward reward/slope ticks for the current component. On a row with no face and no path-graph successor reaching one, discount and the triangular successor recursion give exactly zero. Transitions into such zero rows still contribute their full rate to the Bellman residual. Donor and mother internal changes do not alter the skeleton; testing the union of their possible terminal reward faces is sufficient for safe pruning. In the tested complementary masks **all** skeleton rows can reach a payable face, so this mechanism saved no Bellman rows. For example the medium minimum component has 91,794 nonzero faces out of 98,280 terminal candidates, but all 103,070 skeleton rows survive. Small, medium, and the wider joint envelope all have this saturation. Selecting near-zero faces by magnitude would change the objective and is not used.

Values for donor arity `k` at fixed degree `d` are one array over path and mother states. The degree-gain successor reads the **old value at the same slot** before overwriting it; donor promotion reads the already computed `k+1` array. Read order within a slab is decreasing skeleton rank and then decreasing mother rank. Overflow is precomputed and retained until the top degree is complete. Only the virgin donor slab is retained for the unregistered solve. The peak donor slabs are `side_cap+1`, versus `side_cap(side_cap+1)/2` finite-or-overflow slabs in the dense layout, before including the separately streamed unregistered block. The support pass and allocations never multiply the skeleton array by all donor degrees and arities. Memory is finite with exact descendants at each resolved degree; overflow keeps the frozen rate interval `[0,1/(side_cap+1)]` throughout. This covers arbitrary **actual** exterior degrees, but does not give an exact evaluation of an unbounded degree-dependent tail. No degree `>=18` apex reward or 2417 ownership is used.

Witness charge extraction, `10^18` grid, rate denominators, outward integer rounding, and independent signed 128-bit Bellman residual checks are inherited from the frozen solver. True incoming rate appears on the same arrow as the associated affine reward slope. The `plain` engine keeps separate arrows, `merged` adds intervals only if they have identical successor and zero additional slope, and `strata` adds exact support pruning to `merged`. The sum of intervals gives precisely the attainable sum of rates, so merging does not change even the rounded ratio. A checksum hashes the skeleton, donor/mother state, and all four exact values in traversal order. The streamed traversal has a different checksum convention from the frozen recursive traversal.

## Measured matching costs

Every conditional birth-mother tick matched, not just the sum over components. CPU below sums the printed per-component process times (authentication excluded); slots are the maximum number of simultaneous four-value records, excluding small topology/face arrays and the mmap inputs.

| Cut / engine | States | Arrows | Residuals | Peak value slots | Solver CPU s |
|---|---:|---:|---:|---:|---:|
| Small: frozen | 216,294 | 1,915,450 | 865,176 | 57,324 | 0.178 |
| Small: plain streamed | 216,294 | 2,497,486 | 865,176 | 30,996 | 0.153 |
| Small: merged streamed | 216,294 | 1,915,450 | 865,176 | 30,996 | 0.138 |
| Small: strata | 216,294 | 1,915,450 | 865,176 | 30,996 | 0.140 |
| Medium: frozen | 6,017,604 | 58,102,093 | 24,070,416 | 1,375,170 | 5.567 |
| Medium: merged streamed | 6,017,604 | 58,102,093 | 24,070,416 | 771,120 | 4.035 |
| Medium: strata | 6,017,604 | 58,102,093 | 24,070,416 | 771,120 | 3.982 |
| Medium role 0, donor cap 6: frozen | 4,240,641 | 43,990,010 | 16,962,564 | 4,240,641 | 3.621 |
| Same cap 6: strata | 4,240,641 | 43,990,010 | 16,962,564 | 1,339,926 | 2.592 |

Medium strata-only follow-up against the frozen JSON: `11.44` wall, `7.26+1.12` CPU seconds (authentication `4.415` CPU), peak RSS `142780` KiB. Medium combined ablation: `29.86` wall, `16.74+2.70` CPU, peak RSS `231784` KiB. Cap-six pair: `21.85` wall, `9.67+2.74` CPU, peak RSS `257616` KiB. This establishes a memory route to larger donor caps, not a faster asymptotic row count. The cap-six single role increased its uniform floors from `1.33509e-7 / 6.29833e-8` to `1.57951e-7 / 7.98359e-8`; it was a matched ablation on the medium domain, not a wide replacement.

## Joint support bound and birth law

`--parts envelope` reconstructs each original affine row as `minimum + sum of nonnegative donor excesses`, puts the offsets in the source, and puts each signed slope on the matching incoming gain arrow. An opaque donor group at role `i` uses the **sum of the real off-path reciprocals** as a single control in `[0,count_i/2]`; all terms in that group have the same original slope. The real mother's pressure uses its exact state until overflow and carries the same slope into the same apex gain. This relaxation needs no fifth vertex. All component charges appear once; the CLI disallows combining this result with the separate parts. The upper endpoint bounds the entire positive charge on the stated cut (including the real physical processes), but does not bound support outside that cut.

| Joint full-affine cut `(16,9,9,16)` | Lower-witness full interval per core | Upper-witness full interval per core |
|---|---:|---:|
| Arities `(5,5,5,4)`, stateless | `[1.24443e-11,4.35274e-6]` | `[1.23616e-11,4.29166e-6]` |
| Arities `(4,4,4,3)`, mother cap 3, uniform birth law | `[2.44203e-8,3.70443e-6]` | `[1.23821e-7,4.03523e-6]` |
| Same mother cap 3, physical sided birth law | `[4.96153e-8,3.62423e-6]` | `[2.00096e-7,3.90337e-6]` |

These finite-support ceilings do **not** rule out `2e-6/core`, and none of the joint floors improve the release floors. The widening is not a matched estimate of maternal-memory impact. The meaningful matched comparison in the last two rows changes only birth pooling: `pool_envelope.py` applies the existing exact repaired sided LP to the same physical pre-promotion mother classes via `β=x(b·)` and multiplies by `μ` once. JSON retains the four conditional entry ticks per mother class for further repricing by `embedding-flow`. Nothing assumes an unconditioned donor mean is a live receive constraint.

The stateless envelope ran `11.88+2.20` CPU seconds and 1,302,128 KiB RSS (4,105,079 states, 16,420,316 residuals). Active mother envelope: `10.49+1.66` CPU and 910,736 KiB RSS (6,369,201 states, 25,476,804 residuals). Exact birth reprice: `5.08+0.10` CPU, 95,192 KiB RSS. This quantifies why sign-face pruning alone is an unpromising route to a full larger replacement on shared cores; the remaining lever is stronger physical simultaneous incoming control or a separate efficient bound on degree/arity tails, not more zero-row filtering.

## Reproduce

From the project root, set `R=round2/occupation-neon/pressure-memory/path-fast/strata`, `SO=<scratch-path>`, `PY=./round2/.venv/bin/python`. Build with `g++ -std=c++20 -O3 -Wall -Wextra -Werror -shared -fPIC "$R/solver.cc" -o "$SO"`; prepend `taskset -c <available-cpu> nice -n 19` to builds and replays. Use one job at a time. Each solve authenticated the original blobs (not a saved table surrogate); replays used `OPENBLAS_NUM_THREADS=1`.

```sh
replay() { "$PY" -Bu "$R/replay.py" --library "$SO" --lower b3f18_lo.bin --upper b3f18_hi.bin "$@"; }
replay --engines frozen plain merged strata
replay --degrees 9 7 6 10 --arities 3 3 3 2 --engines frozen merged strata --reference "$R/../medium-entries.json"
replay --degrees 9 7 6 10 --arities 3 3 3 2 --side-cap 6 --parts 0 --engines frozen strata
replay --degrees 16 9 9 16 --arities 5 5 5 4 --mother-cap 0 --side-cap 0 --parts envelope
replay --degrees 16 9 9 16 --arities 4 4 4 3 --mother-cap 3 --parts envelope
"$PY" -Bu "$R/pool_envelope.py" "$R/envelope-wide-mother3.json"
```
