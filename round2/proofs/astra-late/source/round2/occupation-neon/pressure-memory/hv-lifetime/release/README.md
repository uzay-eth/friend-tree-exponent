> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Live hV donor excess after the multi-family freeze

This addendum composes the [e911a1c ledger](../../combined/README.md) with the [six hV donor floors](../old/EAD947B.md). All charges are priced on the original two signed D18 witnesses at their original trials. The physical mean and Palm/age identities and the original exhaustive D18 GMP checks are inherited.

**Outward enclosure: `[0.565679256596412, 0.565769102233118]`.** The two exact endpoints are:

```text
6853561875850472647665032629043939443/12115632305640977999025000000000000000
328740990573504279104102687021150749/581051508956477588875000000000000000
```

The [independent cold audit](../../../../../../../../reviews/cold/astra/checks/hv-atom-sol.md) found the installment correct in scope. It separately checked physical row ownership, tail classification, the parent Palm weight, the six priced role rows and exact recombination. The five production lifetime runs are archived here.

## Support ledger

| Existing or new payment | Allocated charge and support |
|---|---|
| Original higher reserve | Top-owned higher rows on complete focal planted fringes of size at most 12; size 13–20 suffix on three or four cores, further restricted to the interior. |
| Other EDGE and memory payments | First-body and inner-edge reserves, high-parent mixed `g` and virgin `hC`, corrected `f/g` and their residual memory. In particular the `g` promotion that creates an `hV` weight stays with `g`. |
| e911a1c additions | Original `P4` band and complementary six-part `P4` on its excluded band, broad `S3`, and unused positive `Yb`. These are distinct charged embeddings from `hV`. |
| This addendum | Only live exterior donor excess `c p − min(0,c/2)` in the original `hV` row. Pay on a *proved* complete apex fringe size greater than 20, or size greater than 12 with at least five cores. Each qualifying exterior incidence is booked once, with parent, apex-child, middle-child and tip-child distinct. Current apex-degree bands 2–20 and 21–48 are disjoint within the roles that use both. |

The [earlier partition](../../../high-parent/README.md#physical-charge-and-disjointness) is preserved. For spine slots `(d_i,k_i)`, put `n=Σk_i−2`: the lower counts are `Ĉ=3+n` and `Ŝ=Σ(d_i−k_i)+2n`. An exact tracked side donor `(d,k)` adds `k` cores and `d+k−2` vertices to these counts. The gate is increasing, so lower counts that satisfy it certify actual complete-fringe support. The parent is excluded from those counts; its price uses the degree-at-least-three maternal birth class and the conservative Palm fraction `782836785161/1032836785161`. The all-rate `hV` baseline, conditional contrasts, and exploratory `hV` prices are outside this ledger.

## Exact residual and reproduction

From the repository root, with Python 3.12 and the original read-only D18 payloads `L` (mode 2), `U` (mode 1), both specified by absolute paths:

```sh
R=round2/occupation-neon/pressure-memory/hv-lifetime/release
sha256sum --strict -c "$R/release.sha256"
python3 -EBu "$R/verify.py" --payloads "$L" "$U"
```

The verifier authenticates the previous release and reads its negative EDGE trial margins, previously booked reserves, exact joint memory floor and denominator. It checks the sum including all four `e911a1c` additions against the previous exact endpoints, adds the six donor floors to each corresponding net numerator, and checks the equivalent inward move from the previous exact endpoints. Only the final enclosure is rounded outward.

For fresh exact integer donor receipts, set `P` to a Python 3.12 interpreter with NumPy and `OUT` to an existing empty scratch directory. From the archive root, with `R`, `L` and `U` as above:

```sh
S=round2/occupation-neon/pressure-memory/hv-lifetime/old/replay.py
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
run() {
  tag=$1; shift
  /usr/bin/time -f 'wall=%e rss=%M exit=%x' "$P" -EBu "$S" "$L" "$U" --frozen "$@" > "$OUT/$tag.log" 2>&1
}
run apex20 --caps 20 5 5 --arity 6 2 2 --atoms apex --atom-cap 4
run apex48 --caps 48 5 5 --arity 6 2 2 --atoms apex --atom-cap 4 --host-min 21
run parent20 --caps 20 5 5 --arity 6 2 2 --atoms parent
run parent48 --caps 48 5 5 --arity 6 2 2 --atoms parent --host-min 21
run sides --caps 11 10 10 --arity 3 2 2 --atoms middle tip
python3 -EBu "$R/verify.py" --compare-replays "$OUT"
```

To repeat the independent checks on the packaged receipts:

```sh
OPENBLAS_NUM_THREADS=1 "$P" -EBu round2/checks/hv-atom-sol.py "$L" "$U" --archive "$R"
```
