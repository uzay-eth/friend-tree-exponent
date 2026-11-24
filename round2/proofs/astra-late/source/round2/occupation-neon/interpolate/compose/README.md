> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Claim 21 addendum: repaired pointwise witness and all-seven composition

**Outward interval `[0.565676930574, 0.565771634792]`.** The fixed
pointwise trials are `.565675148075` and `.565773278615`. Compared with
the previous high-parent ledger, the lower endpoint moves inward by
`0.000000296083` and the upper endpoint by `0.000000239672`.

| Endpoint | Exact fraction before outward rounding |
|---|---|
| Lower | `15350577873033193312515081376688304908311/27136651758904633839660987238280000000000` |
| Upper | `15470132261475507150976552475495052519651/27343421462238657786289217467400000000000` |

## Accounting

The first-body corrector uses only charged `f` at the original trial
`λ`. All seven other charged rows at that **same** `λ` are still present:
`ε(ℒ−s)Φnew ≥ m(C−1)+Σ non-f`, with `m=s z(J)−δM≥0` and
`Φnew=Φ+ε(Σ_nonfounder(z−z(J))−δM/s)`. The actual physical process,
original transfer terms, and ownership of charges do not change. No
stationary `f` floor is added. Integrating the signed inequality, with
the same normalized Laplace limit as the prior physical release, yields
`ε(μ−s)πΦnew ≥ m + retained`. The constant and one founder vanish after
normalization; `πΦnew` is positive by the pointwise comparison for the
new potential.

| Paid once | Lower | Upper |
|---|---:|---:|
| Six top-owned higher families, size ≤12 | 1.2769028496806e−6 | 1.5010294586750e−6 |
| Size 13–20, three/four core suffix **plus** true inner-edge floor | 5.89618197103865e−7 | 4.66358431549090e−7 |
| High-parent `g/hC`, joint | 5.52030676155105e−7 | 2.79999563548312e−7 |

The suffix's actual `g` inner edges may pay before size 13; its six higher
families may not. The high-parent `g` has a tail parent; its `hC` is
beyond the suffix degree seam and cannot fit the size-12 fringe. Thus
the old [support partition](../../high-parent/README.md#physical-charge-and-disjointness)
applies without edits. The suffix uses the larger physical stable-age
replay, not a sum with the separate 19-label residual price. At the
matched size-18 cut the prior residual prices were `2.48e−7` / `1.89e−7`
versus physical `4.90e−7` / `3.67e−7`. No independent new pointwise
all-eight correction cost is booked: all seven are priced by occupation
after the pointwise `f` correction. In particular, an edge-corrected
witness cannot just keep the entire `g` reserve from this ledger.

For the denominator use the old all-degree stationary ceiling plus an
upper bound on `π[ε(z−z(J))]`. First take a uniform upper `H` over the
finite census and the full arity tail. Then pay the nonnegative gap
`H−ε(z−z(J))` through the physical stable-age census, zeroing its reward
beyond degree 64. This gives signed offsets at most `−4.33365768308e−7`
and `+1.07790193289e−6`, respectively. No constant is subtracted from the
limit; its normalized contribution is exactly zero. Uniform-only gives
`[0.565676930572,0.565771634797]`.

## Upper failure and matched repair

The source of the dramatic withdrawn gain was an upper `k=1` floor
printed as `861071198001/10^15` instead of `861071198001/10^18`.
At the withdrawn trial the planted replay produces a seemingly positive
Bellman margin but exactly this actual whole-tree signed drift:

```
-902875978288208312083656815181293852064686619898589/1097886886904005745377459200000000000000000000000000000
```

Parents are `(-1,0,1,2,1,0,5)` and leaf stocks `(62,63,0,1,4,65,1)`.
At the repaired trial the same tree passes. The replay uses actual
labelled donor arrows and full eight-motif potential differences. All
five tail floors and all 2,079 degree-64 `f` rows **per side** are
re-extracted from the SHA-authenticated binary and compared with the
hash-pinned row/floor source. Matched margins at the repaired upper trial:

| Correction variant | Signed per-core margin |
|---|---:|
| Interior first-body only, zero tail | −1.164091e−6 |
| Interior includes charge deficit, zero tail | −6.035758e−8 |
| Add stock continuation only | −4.907544e−9 |
| Add charge deficit continuation only | −5.550301e−8 |
| Both tails | +8.970118e−13 |

This separates the genuine small tail contribution from the invalid
thousandfold inflation. The cold reviewer separately rebuilt both
Bellman witnesses, checked 4,521 corners and 17 physical all-eight cases
per side, reproduced the exact negative regression, and agreed on both
composed outward endpoints. See the [independent audit](../../../../../../../reviews/cold/astra/checks/interpolate-sol.md)
and its [fixture](../../../checks/interpolate_sol.py).
The producer evidence is from `2f74` and the separate GPT-family audit is
from `7b08`; these are review identifiers, not Git pins.

## One-command archive replay

From a clean archive of this addendum, set `PYTHON`, `D18_LOWER` and
`D18_UPPER` as in the [original release](../../release/README.md#replay).
Use 64-bit Python 3.11+ with NumPy and assertions enabled.

```sh
OPENBLAS_NUM_THREADS=1 nice -n 15 taskset -c 8 "$PYTHON" -Bu round2/occupation-neon/interpolate/compose/replay.py --lower "$D18_LOWER" --upper "$D18_UPPER"
```

The command verifies the [addendum manifest](inputs.txt) and the inherited
manifests, authenticates both originals, runs the independent audit,
re-extracts the signed `f` rows and stock floors, and compares the finite
backward values with a forward discounted-flow enclosure. It replays the
exact negative and repaired positive tree regressions, tail ablations,
all-degree ceilings and all three retained reward blocks, and checks the
exact endpoints above. `--side` is a diagnostic subset of this command.

The earlier complete per-sign transcripts are [lower](lower.log) and
[upper](upper.log). The cold all-seven replay uses the archived base GMP
receipts (`round2/attacks/neon1/outputs/neon1-production-rehearsal.tar.gz`, historical reference)
and inherits the stationary and stable-age results. It does not repeat the
exhaustive GMP row scan. The [immutable Neon proof pin](../../../../../../../reviews/cold/astra/occupation-neon/check/neon-proof-pin.md)
identifies commit `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, path
`round2/proof-b3f-motif4.md`, blob `3bbda54c3cbdbd235a6cf4af8ecf80ca75d16ae7`.
