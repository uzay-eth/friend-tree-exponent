# Wider low-apex complementary P4 replacement

**Outward enclosure: `[0.565679672983674, 0.565768735586042]`.** The exact endpoints are:

```text
6853566920645426224860150084809360393/12115632305640977999025000000000000000
10848445658578037268151132921697974717/19174699795563760432875000000000000000
```

Start from the [live hV freeze](../../hv-lifetime/release/README.md) at `f7cc2e1c77192ecea8f231519f4afd2bf50aab9b`. Replace only the six-part complementary medium P4 floor, read exactly from the [old ledger](../../combined/donor-lift.log). Per witness the new net numerator is the old net numerator minus the booked medium floor plus the pooled wide floor. The original signed trials and physical denominator ceilings are unchanged, as are root P4, hV, S3, Yb and the older EDGE and memory payments.

| Witness | Subtract medium P4 | Add pooled wide P4 |
|---|---:|---:|
| Lower | `359297001126132137/10^24` | `23575773396079784114802186917/25500000000000000000000000000000000` |
| Upper | `6963149521191559/25000000000000000000000` | `43746697941322995419/56100000000000000000000000` |

## Physical support and certification

The [native producer](../README.md) uses current degrees at most `(16,9,9,16)` and arities at most `(4,4,4,3)`. It prices only the positive six-way split of the original P4 row, outside the frozen higher-reserve support and outside the entire root P4 band. For the spine, let `n=Σk_i−3` and `S_min=Σ(d_i−k_i)+2n`. The gate requires `S_min>20` or `S_min>12` and `n>0`, so it excludes the old size-at-most-12 and moderate three-/four-core supports. It also excludes all rows with degrees within `(12,8,6,16)` and arities `(1,k,1,0)` for `k≤3`. The live reciprocal of a marked donor affects the same gain event and its own reward. The two centered hV-to-P4 transfers are part of the original P4 accounting; they do not rebook the distinct hV exterior-arrow excess.

All six components are pooled against a **single** physical birth-mother law. The refinement replaces only the mother component by degree classes 2 through 12 and the degree-13+ overflow. Other components use their original degree-3 overflow throughout the refined sub-bins. The birth-law LP with cutoff 18 lifts stationarity, sided edge flows and the physical Palm lower kernel; exact rational reduced-cost repair certifies the bound after floating dual proposals. No embedding LP, high-apex calculation, conditional contrast, or other candidate enters this release.

The [separate cold review](../../../../checks/path-fast-sol.md) checked the native coefficient identities and arithmetic bound, full medium UBSan row parity, sampled wide UBSan/Fraction parity, the full native Bellman replays, and exact LP/ledger recombination. The separate isolated producer replay ran both Python parity controls, the medium comparison, the wide head and tail and the refined mother shard. All four conditional JSON files were byte-identical to the pinned originals, and the pooled endpoints agreed exactly with the cold check. An [optional exhaustive Python charge extraction](optional-row-parity-aborted.log) on the wide head was stopped before it completed a component; it supplied no parity receipt. The original exhaustive D18 GMP checks and physical stationary/Palm identities are inherited.

## Reproduce from an archive

Use GCC 13.3.0, CPython 3.12.3 with NumPy 2.3.5 and SciPy 1.16.3, and keep assertions enabled. Set `P` to that Python, `L`/`U` to the original mode-2/mode-1 D18 payload paths, and `OUT` to a fresh empty directory. From the archive root:

```sh
R=round2/occupation-neon/pressure-memory/path-fast/low-release
sha256sum --strict -c "$R/release.sha256"
"$P" -EBu "$R/verify.py" --payloads "$L" "$U"
bash "$R/reproduce.sh" "$P" "$L" "$U" "$OUT"
```

`verify.py` authenticates the inherited manifests and exact receipts, checks all 24 medium endpoints and the cold and producer pooled fractions, and recomputes the signed residuals and final outward rounding. The full replay compiles the native solver, authenticates the payloads, runs the shards sequentially and checks the exact state counters, checksums, prices and conditional files. The wide tail shard needs about 3 GiB resident memory. The compiled library and D18 payloads are external to the release.
