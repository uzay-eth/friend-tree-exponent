# Exact multi-family replay

**Outward enclosure: `[0.565678493279215, 0.565770152011715]`.**

| Sign | Exact endpoint |
|---|---|
| Lower | `6853552627779976665165032629043939443/12115632305640977999025000000000000000` |
| Upper | `328741600548942531604102687021150749/581051508956477588875000000000000000` |

The ledger starts from the authenticated `ead947b` clean-replay receipt. It reconstructs each EDGE net numerator from the signed pointwise margin and frozen reserve, adds the exact existing joint memory floor once, and verifies that this reproduces both exact `ead947b` endpoints. It then uses the same physical potential denominator on each side.

The additional terms are the unpaid root P4 band, the complementary six-part P4 lift, the corrected cap-five S3 floor, and the unused Yb floor. The lift excludes the entire root band on present-time descriptors. Its minimum and five donor contributions partition the instantaneous charge; they are added once. S3 requires `k≥4` and a lower bound of 13 on the actual focal fringe size. Yb uses the independent signed replay's exact disjoint positive floors. The three families are distinct allocated charges. No signed contrast, overlapping diagnostic, or nested mask is a payment.

The [S3 cold audit](../../../checks/s3-lifetime-sol.md) found a missing half-rate corner in the original tagged fourth-child diagnostic. The [repaired producer](../../s3-lifetime/README.md) checks the defect on the physical upper crown and reproduces the independently audited broad floors. Only the `broader tip cap 5` rows enter this ledger. The broad calculation uses the shared host receive and has no tagged-baseline term.

## Verify and reproduce

From the repository root, with Python 3.12 and the original read-only D18 files at absolute paths `L` (mode 2) and `U` (mode 1):

```sh
R=round2/occupation-neon/pressure-memory/combined
python3 -EBu "$R/verify.py" --payloads "$L" "$U"
```

This hashes the complete inputs, checks the inherited manifests, reads exact rational fields from the pinned receipts, and rounds outward at 15 decimal places. The receipt-only variant omits `--payloads`. The original exhaustive D18 GMP checks are inherited.

For an independent archive rerun of the new prices, point `OUT` to an empty scratch directory and set `P` to a Python 3.12 interpreter with NumPy. Use the following from the archive root. Keep assertions enabled and run one at a time if memory is limited.

```sh
R=round2/occupation-neon/pressure-memory/combined
M=round2/occupation-neon/pressure-memory
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$P" -EBu "$R/verify.py" --inputs-only
/usr/bin/time -f 'wall=%e exit=%x' "$P" -EBu "$M/path_band.py" \
  --lower "$L" --upper "$U" --degrees 12 8 6 16 --middle-arity 3 --donor 2 --mother 4 \
  > "$OUT/p4.log" 2>&1
/usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M exit=%x' "$P" -EBu "$M/donor_lift.py" \
  --lower "$L" --upper "$U" --degrees 9 7 6 10 --arities 3 3 3 2 \
  --side-cap 3 --mother-cap 3 --budget 3000000 > "$OUT/donor-lift.log" 2>&1
/usr/bin/time -f 'wall=%e exit=%x' "$P" -EBu round2/occupation-neon/s3-lifetime/replay.py \
  "$L" "$U" --tip-cap 5 --large-only 4 > "$OUT/s3.log" 2>&1
/usr/bin/time -f 'wall=%e exit=%x' "$P" -EBu round2/checks/yb-lifetime-sol.py \
  "$L" "$U" > "$OUT/yb.log" 2>&1
"$P" -EBu "$R/verify.py" --compare-replays "$OUT"
```

Replay verification compares the exact floor pairs and all six lift component intervals, not decimal displays or wall times.

The isolated replay used an `ead947b` archive with only the manifest's new sources overlaid at their original paths. [Root P4](p4.log) finished in 39.74 s; [all six medium lift parts](donor-lift.log) in 437.78 s; [repaired S3](s3.log) in 93.30 s; and the [independent Yb solve](yb.log) in 44.08 s. All exited zero and their exact prices matched the pinned ledger. The [lower](s3-cold-lower.log) and [upper](s3-cold-upper.log) S3 companions independently reproduced the same broad floors on the prerepair source; the broad solver was not changed by the tagged repair.

## Independent checks

- [Root P4](../../../checks/p4-band-sol.md): independent physical row reconstruction, wide census and sign-corner checks, and a matched single-row solve. The wide Bellman floors come from the authenticated producer replay.
- [Complementary P4](../../../checks/p4-lift-sol.md): independent medium mask, labelled arrow and raw score checks, and exact recombination of the original six-part receipt. The isolated producer replay repeats all six medium parts. The smaller matched lift is a control only.
- [S3](../../../checks/s3-lifetime-sol.md): independent broad Bellman comparison on both signs, plus physical size and multiplicity controls. The repaired fourth-child diagnostic is outside the paid ledger.
- [Yb](../../../checks/yb-lifetime-sol.md): independent all-age signed solve on the stated mask. Only the positive floors on unused support enter the ledger; the signed contrast still straddles zero.
