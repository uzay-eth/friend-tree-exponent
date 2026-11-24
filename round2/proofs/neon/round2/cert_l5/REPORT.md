# B3-FullH + coarse endpoint L5: exact v2 replay

**Authors:** Gemini 4 Neon (lead) and GPT-6.2 Sol (implementation). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Review and dependencies:** Gemini/GPT development checks of the D10/D12 benchmarks; the separate GPT cold review covers the production profile. The benchmark harness uses the [version-2 drift theorem and GMP checker](proof-b3f-l5.md); see [Neon Claim 4](../../README.md) for scope.

The builder uses binary64 to fit and polish candidate coefficients in place. The verifier reads the integers back and proves their drift and comparison inequalities without floating point. A candidate file alone is not a certificate until it passes the verifier.

## Implementation notes

- The first 96 bytes retain the v1 layout; the 24-byte extension is serialized field by field. File length and the independently recomputed shape count are checked before allocating payload arrays. Two buckets require `Dc < D0`; all four L4 degree caps must be at least `Dc`, and all core-child caps must be at least `Dc−1`. In particular, `KS3=Dc−2` leaves the exact S3 apex `(L,k,h)=(0,Dc−1,0)` without a row to absorb a spectator-side birth. The corrected D0=10 test uses `KS3=5`.
- L5 row order is rooted shape, apex `h`, then dense rank. Identical child shapes are combined by colex multiset rank. Apex BIG uses one state per bucket, although its L4 spectator fiber includes both physical apex `h` values. Fibers enumerate the exact paired L4 control differences; per-thread caches store a lossless integer only after the direction-selected GMP extremum is known. The whole ordered tuple and marked slot are in the cache key.
- L4 exact-promoter births add the L5 value independent of receiver mask. Only eligible exact external-credit intervals change to `c_ext/Dc`; the Bernstein tail and overflow budgets retain their original bounds. L5 uses a GMP numerator over `Q·lcm(1,…,D0)`, including the zero-or-crossing endpoint corner and reciprocal credits for the correct neighbor bucket. A coarse endpoint promotion does not change its bucket. A coarse interior is unsupported.
- Every supported motif is charged to its topmost interior. For an exact charge degree `D`, `D·binom(D−1+4(Dc−2),3)/4` bounds the number of such charges. The global envelope was too loose for the D0=10 upper endpoint despite zero row failures. The verifier separates the extrema by exact charge state and rooted shape, including zero. For each group of `m` identical child shapes, count `binom(k_remaining,m)` choices, then multiply the downward capacities for those children using at most `Dc−2` forward branches. An endpoint apex necessarily occupies the charge vertex's unique parent and costs no choices. The sum of shape capacities is checked against the original `B(D)`. The resulting weighted extrema are added in `mpz_t` before the v1 zero-sum recharge. The Python replay independently enumerates physical subsets and checks every observed per-shape count, grows every concrete tree at every possible recipient, and checks both comparison inequalities, newborn ownership and crossing rates.

## Reproduction

```sh
bash test_cert_l5.sh 10 12
```

The harness pins one CPU, disables OpenMP worker expansion, reports elapsed seconds and KiB RSS, and removes its private scratch directory. It compiles with `-Wall -Wextra -Werror` and uses ASan+UBSan for the D0=10 verifier and controls. Ordinary runs have a virtual-memory cap; ASan reserves a large virtual shadow mapping, so its actual RSS is reported instead. Mutable files and executables live under the temporary directory.

The D0=12 historical no-L5 float width was `0.00013360681999630142`; the comparable split-L5 float width was `0.00011427489906101762073831196175888180732727050781250`. These are diagnostics, not exact intervals. Production-profile timings likewise must be labeled projections because the production layouts exceed the local RAM limit.

## D0=10 results

Corrected `Dc=6`, two buckets, receiver mask `3`, `KS3=5`, `Q=10,000,000`, cushion `5e-7`, 22 fitting sweeps plus five fixed-target sweeps. Both endpoints store 2,127,088 L5 rows and use 145,990,696 bytes per file. Certified interval **[5656316, 5658254] / 10,000,000**, width **0.0001938**.

| Endpoint | `4*S*c0` | `S*C0` | `S*Phi(T3)` | Builder s / KiB | GMP s / KiB | ASan+UBSan s / KiB |
|---|---:|---:|---:|---:|---:|---:|
| Lower | 92796844940 | 183528483027 | 146046735703 | 164.78 / 161636 | 37.80 / 153856 | 81.77 / 787892 |
| Upper | 102046552236 | 194911759971 | 146046935242 | 158.14 / 161632 | 35.93 / 154112 | 93.21 / 787364 |

`S=100,000,000,000`. Both ordinary and ASan+UBSan GMP replays reported `total_fails = 0`, with 34/34 exact concrete trees each. Python independently reduced the comparison bounds and passed 38/38 concrete trees each. Inward overrides by `20/Q` failed 2,042,023 lower and 2,897,705 upper checks. An in-range L5 integer corruption caused five failures. Structural rejects covered `reserved`, truncated and extra bytes, wrong row count, invalid `Dc`/`buckets`/mask and `KS3=Dc−2`; the builder rejected matching invalid CLI cases. Wrong TAIL fiber, inverted source, and omitted L4 newborn caused 137,020, 915,170 and 1,505,119 failures, respectively. The independent Python term tests rejected omitted crossing and newborn. Omitting the adverse crossing term alone relaxes the direction-selected row inequality, so the crossing control checks the actual physical event against the planted interval instead of expecting a fitted row failure. No sanitizer diagnostics were emitted.

## D0=12 results

`Dc=7`, two buckets, mask `3`, 18 fitting sweeps plus five fixed-target sweeps. Both endpoints store 14,215,993 L5 rows and use 773,937,928 bytes per file.

**Certified interval: [5656682, 5657836] / 10,000,000. Exact width: 0.0001154.** It is 13.6272% narrower than the historical no-L5 D0=12 float width and 0.9846% wider than the split-L5 float width. Each endpoint passed all nine GMP tables and 34/34 concrete trees. The independent Python checker passed 38/38 concrete trees per endpoint. Overrides `5656702/10000000` (lower) and `5657816/10000000` (upper) rejected with 4,847,273 and 7,239,621 failures.

| Endpoint | `4*S*c0` | `S*C0` | `S*Phi(T3)` | Builder s / KiB | GMP s / KiB | Python s / KiB |
|---|---:|---:|---:|---:|---:|---:|
| Lower | 76800849604 | 183926565482 | 146044567320 | 795.18 / 781768 | 168.81 / 767616 | 7.25 / 772352 |
| Upper | 91168234836 | 196279029033 | 146045587918 | 778.57 / 781768 | 158.54 / 767232 | 7.44 / 772352 |

The exact constants are the displayed integers divided by `4*S`, `S`, and `S` respectively, with `S=100,000,000,000`.

## D0=20 and D0=24 production projections

Count-only commands are at the end of the harness. These use `Dc=7`, two buckets, `KS3=8`, and the production L4 caps shown there. The byte counts are exact for these configurations. Sweep bytes include the live coefficient arrays, layout, scratch and the larger of the staging or three-worker caches; verifier bytes are the live payload, layout, `dLf` and three-worker L5 caches. Neither is an RSS measurement; runtime, allocator, GMP temporaries, stacks, page cache, empirical initialization and serialization overhead are excluded from the stage-specific structural totals.

| D0 | File bytes per endpoint | Builder sweep bytes | Verifier L5 phase bytes | Projected build s | Projected GMP s |
|---:|---:|---:|---:|---:|---:|
| 20 | 11,370,156,264 | 11,448,307,212 | 11,396,840,972 | 2,901 (~48 min) | 802 (~13 min) |
| 24 | 26,942,980,144 | 27,033,465,520 | 26,969,804,016 | 7,209 (~2.0 h) | 1,899 (~32 min) |

Timing estimates assume ideal three-way scaling. Builder extrapolation uses the mean D0=12 L5 sweep subtotal (373.964 s per endpoint) as a fixed-row term; the rest scales by v1 payload size, with 27 total sweeps at D0=20 and 29 at D0=24. GMP time uses mean D0=12 wall time scaled by total rows. Larger fibers, memory bandwidth, imbalanced shape jobs, serial work and shared-host contention can materially change both estimates. No production-sized tables were allocated.
