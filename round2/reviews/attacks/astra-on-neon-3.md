# Astra cold audit of original Neon Claim 3

**VERDICT: HOLDS.** At `ft2-neon@51f4cfadac675fdf1cef84d44935eebb5456aabb`, both hash-pinned B3-FullH-20-L4 certificates replayed in full. Independent comparison, physical-generator checks and the drift-to-product argument confirmed the bracket **`[0.565680, 0.565768]`** and the stated finite-size bounds for every integer `n ≥ 3`, without an external prior bracket. This is the original Neon claim on `neon.md`.

## Production replay

Fresh `git archive` of the exact pin, then `gcc -O3 -march=native -fopenmp -Wall -Wextra -Werror round2/cert/cert_b3f_verify_gmp.c -lgmp -o verify-gmp`. The two 11,190,881,056-byte payloads matched these full SHA-256 hashes and ran serially with `nice -n 15 taskset -c 2,6 ionice -c 3 stdbuf -oL ./verify-gmp BLOB --threads 2 --tree-check 30`:

| Witness | SHA-256 | Verifier result | Wall (s) | Peak RSS (KiB) |
|---|---|---|---:|---:|
| `b3f20_lo.bin` | `eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469` | `PASS`, `total_fails = 0`, `34/34` trees, exit 0 | 2253.17 | 10921344 |
| `b3f20_hi.bin` | `55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2` | `PASS`, `total_fails = 0`, `34/34` trees, exit 0 | 1750.52 | 10931840 |

Both logged `6,432,404 / 121,693 / 27,746,004 / 13,880,592 / 354,856,788 / 540,925,385 / 331,822,656 / 122,901,075` checks in `f/g/hV/hC/P4/Yt/Yb/S3` order with zero failures in each section. The `g` sweep checks admissible pairs; storage is a dense `523²` matrix. An independent array reducer confirmed **1,398,838,433 potential entries + 21,687 auxiliary entries = 1,398,860,120 total**.

## Comparison, drift and products

The independent reducer recovered the exact verifier integers:

| Witness | `4Sc0` | `SC0` | `SΦ(T3)` |
|---|---:|---:|---:|
| Lower | 113472636536 | 184219674916 | 146041684821 |
| Upper | 115454594604 | 200854054323 | 146044513514 |

The lower charge includes the zero-sum correction `Σ(k−h)/4=0`; using the clipped `KF+1` on overflow preserves a lower bound. Omitting the correction gives a negative raw lower comparison minimum. For the growth law, directed edges from donor `u` each have rate `1/D_u`, so total rate is `n`. Grouping physical events gives gain rate `L+Σ_core-neighbor 1/D` and bud rate `L/D`. A bud preserves the old core vertex's degree and creates a core child `(1,0,1)`. Each centered transfer removes and restores the same donor contribution on the same eligible pre-event domain. The nonnegative Bernstein weights and positive denominator controls extend the row bounds to physical tail degrees; in arity overflow the gain factor `(D−c+ρ)/(D+1)` lies in `[0,1]` regardless of clipped representation.

An independent exact-rational evaluator passed 84 physical fixtures per sign through degree 25,020. These exercised all eight motif families and seams at `D=19/20`, `k=19/20`, `k=200/201` and `L=16000/16001`. On 57 fixtures per sign, direct directed-arrow and grouped gain/bud sums agreed component by component. The planted wrong recipient rate and wrong degree after promotion changed the result. A separate UBSan tail probe matched 168 exact physical interpolation controls; the planted transfer-sign reversal disagreed on all 84 of its controls.

Since the conditional increment is `ℒΦ/n`, the signed bounds iterate from `T3` by the factor `1+λ/j`. Divide the lower potential bound by lower `C0` and the upper by upper `c0`; the resulting factors are exactly `146041684821/184219674916` at `λ=7071/12500` and `584178054056/115454594604` at `λ=70721/125000`. The product is empty at `n=3`. Every exact mean from direct labelled-history enumeration through `n=10` lay between the products. The logarithms of the two products are `λ log n + O(1)`, giving the claim's exponent endpoints under the existing mean-exponent/a.s. equality theorems.

## Repair scope and receipts

The pinned proof and board use the correct literal counts. The pinned source replaces the old fixed pair staging with a full Cartesian-product allocation: the formerly problematic generic Ybot `ksum=24` bucket needs 30,899 slots, beyond the old 16,384. Checked allocations fail before use; the ABI guard enforces the file layout, and original and override headers are validated in that order. Nine small invalid-input probes passed with expected nonzero exits and no UBSan failure. These checks cover those repair paths and the pinned production claim; there was no exhaustive generic-profile build or tighter endpoint claim.

Replay scripts and detailed evidence: `round2/attacks/neon3/README.md`. Production receipt: `production.log` (SHA-256 `1b8f37f67875a75bd4a13dbd40f272d2763280204a2fb6ca70b7cc5d576af514`). The completed independent binary/physical rerun took 143.07 s at 10,973,176 KiB peak RSS, exit 0, on one worker after production.
