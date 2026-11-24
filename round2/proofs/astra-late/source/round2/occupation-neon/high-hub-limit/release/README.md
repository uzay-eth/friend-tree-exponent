> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Claim 21: high101 all-arm addendum

**Outward enclosure: `[0.565680279007146, 0.565768284181773]`.** Baseline: [`aa511af202d0d64cd1a08f900e0516965d9db952`](../../pressure-memory/fork-lifetime/yt-release/README.md). Book only the `all` floors from the two original D18 witnesses on nonfounder apices with **current** degree `D≥101`, `k≤12` actual core children, and a complete first arm of at most seven ordinary vertices. Each old arm carries its continuation once across `100→101`; later promotions register one newborn arm each. The split `old`/`young` rows are diagnostic. Masked ceilings only bound this mask. The D32 total, overlapping younger cut and high P4 `18≤D≤80` are excluded.

| Original witness | Exact `all` floor/core | Unchanged physical denominator |
|---|---:|---:|
| Lower (mode 2) | `52327439425129643019077301/100000000000000000000000000000000` | `1357493815758092773/1000000000000000000` |
| Upper (mode 1) | `8091751687185104316481949937/50000000000000000000000000000000000` | `273436004214812983/200000000000000000` |

Add the lower floor divided by its denominator to the baseline lower endpoint, and subtract the oriented upper floor divided by its denominator from the baseline upper endpoint. The exact endpoints are:

```text
14157243308580160981686611486065459691229739/25026934531690315946179968825000000000000000
44818646818747647275382808904477397643909921/79217319301603159409316399750000000000000000
```

[835c's independent acceptance](../../../../../../../reviews/cold/astra/checks/high-hub-sol.md) checked the nonfounder mask, owner disjointness, seam registration, arbitrary-degree bound and signs. The original bulk lower and upper schedules passed in 98.33 and 59.43 s from an isolated baseline archive with the pinned producer source. Separate exact schedules passed in 31.89 and 34.25 s; full-host physical replay passed in 3.46 s, with 1,092 matches per witness. Fresh companion solves passed in 35.87 and 41.82 s, matching all 36 exact endpoints, 420 extra degree checks per sign and 356 labelled controls. All five exact/physical/cold outputs matched the companion's authenticated receipts byte-for-byte.

The producer's [original manifest](../all-arms/frozen.sha256) is preserved. Its combined text has one extra zero in the denominator of the lower-witness `ceiling old` row total; `lo.log` gives the regenerated sum `16094889408445993677/20000000000000000000000000`. The three family values, both `all` floors and all other rows agree. The released `lo.log`/`hi.log` are the authoritative exact receipts.

## Verify from a clean archive

Use Python 3.12+ with NumPy. From the repository root, with both original payloads available:

```sh
H=round2/occupation-neon/high-hub-limit
R=$H/release
P=python3
L=b3f18_lo.bin U=b3f18_hi.bin
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
sha256sum --strict -c "$R/release.sha256"
"$P" -EBu "$R/verify.py" "$L" "$U"
# Original degree-101 producer schedules (Linux affinity as run):
taskset -c 3 nice -n 15 /usr/bin/time -p timeout 150s "$P" -Bu "$H/bulk.py" "$L"
taskset -c 8 nice -n 15 /usr/bin/time -p timeout 150s "$P" -Bu "$H/bulk.py" "$U"
# Regenerate the exact fractions and replay the cold recurrence:
"$P" -EBu "$H/all-arms/exact_receipt.py" "$L" | diff -u "$R/lo.log" -
"$P" -EBu "$H/all-arms/exact_receipt.py" "$U" | diff -u "$R/hi.log" -
"$P" -EBu "$H/physical_check.py" "$L" "$U" | diff -u "$R/phys.log" -
"$P" -EBu "$R/verify.py" "$L" "$U" --solve lo
"$P" -EBu "$R/verify.py" "$L" "$U" --solve hi
```

The companion's external transport source was authenticated at freeze time for the overlap review. It is not imported or used to price high101; the archive verifier skips only that external file pin. Both original D18 payload hashes and all computational pins remain checked.
