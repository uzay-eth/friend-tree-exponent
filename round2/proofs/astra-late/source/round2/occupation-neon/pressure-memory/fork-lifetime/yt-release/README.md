> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Unpaid Yt exterior-donor installment

This addendum starts from [the wider low-apex P4 release](../../path-fast/low-release/README.md), `f4f7ea3e4afb951dea1ae56727ef2c634043cb2d`, and adds only the [original Yt donor excess](../old/README.md). The **outward enclosure is `[0.565679893536223, 0.565768402553188]`**. Its exact endpoints are:

```text
14157233661424604900081722725993405911089489/25026934531690315946179968825000000000000000
22409328097906915347513239293870120527589541/39608659650801579704658199875000000000000000
```

Add `3092299274719269569703/10328365000000000000000000000` and `9405336866046559940037/20656730000000000000000000000` per core to the corresponding signed net numerators. The respective physical denominators stay `1357493815758092773/1000000000000000000` and `273436004214812983/200000000000000000`. Both original EDGE trials and every earlier payment are kept. Round outward after exact composition.

## Physical charge support

The selected crown is apex, trunk and an unordered pair of distinct tips. Tag the tips in chronological order, and tag each exterior donor at its actual birth, except the mother who is already present. Pay only the excess `c p − min(0,c/2)` of its reciprocal `p` with its own signed port coefficient `c`, once per incidence. The same reciprocal enters the receiver's gain and holding rate. The all-rate box floor is not paid. The four roles in this release are mother, apex child, trunk child and the two tips' child directions pooled by the same recurrence.

The complete planted apex fringe excludes its mother. With `m=Σk_i−3` the lower counts are `S_min=Σ(d_i−k_i)+2m` and `C_min=4+m`; a registered side donor adds `d+k−2` and `k` respectively, retaining only attained counts at overflow. Credit requires the lower counts to prove **`S>20` or (`S>12` and `C≥5`)** on the complete fringe. These monotone bounds exclude the frozen size-12 and size-13–20 three-/four-core reserves. The independently audited true-fringe example has unmarked counts `(9,5)` and becomes eligible at `(13,6)` with a marked `(5,1)` side donor. The proof deliberately leaves ambiguous hidden fringes unpaid.

Incoming centered transfers from `hV` and `hC` are already part of the original `Yt` row. The `hV` release books excess only on its free ports; the centered transfer to this larger motif does not duplicate it. Root and wide complementary `P4`, `S3`, unused `Yb`, and the EDGE and memory rows pay different allocated charges. Neither an exploratory Yt mask nor a diagnostic ceiling or signed contrast is included. The original exhaustive D18 GMP checks and physical stationary/Palm identities are inherited.

## Reproduce

The [independent cold audit](../../../../../../../../reviews/cold/astra/checks/yt-donor-sol.md) accepted this exact scope, comparing 78,650 signed row/port instances and independently solving all four components. The isolated producer replay below matched all sixteen exact interval entries, 1,141,312 states, and 4,565,248 directed inequalities.

From the archive root, set `P` to CPython 3.12 with NumPy and `OUT` to an empty absolute scratch directory. The pinned producer expects the original mode-2 and mode-1 D18 blobs at `b3f18_lo.bin` and `b3f18_hi.bin` and authenticates them on each invocation. Keep assertions enabled.

```sh
M=round2/occupation-neon/pressure-memory/fork-lifetime
R=$M/yt-release
L=b3f18_lo.bin U=b3f18_hi.bin
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
sha256sum --strict -c "$R/release.sha256"
"$P" -EBu "$R/verify.py" --payloads "$L" "$U"
run() {
  tag=$1; shift
  /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M exit=%x' "$P" -EBu "$M/old/replay.py" \
    --caps 16 5 5 5 --arity 5 3 2 2 --parts "$@" > "$OUT/$tag.log" 2>&1
}
run mother-apex mother apex
run trunk-tips trunk tips
"$P" -EBu "$R/verify.py" --payloads "$L" "$U" --compare-replays "$OUT"
```

Run `"$P" -EBu round2/checks/yt-donor-sol.py "$L" "$U"` for the separate cold solve. The release verifier authenticates the inherited ledger and the archived exact receipts and checks both the net-numerator and endpoint-difference formulas against the cold arithmetic.
