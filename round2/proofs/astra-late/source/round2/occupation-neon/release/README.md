> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Claim 21: Neon D18 occupation release

**Certified outward interval:** `[0.565675793084, 0.565772420378]`. These endpoints combine the first-body and size-at-most-12 higher-order rewards with the improved all-degree stationary denominator for each original Neon witness.

| Endpoint | Exact fraction before outward rounding |
|---|---|
| Lower | `27411699836178314229824868195233827/48458322189659161283895969899000000` |
| Upper | `322294615524091104807263734480727/569654164670823324688527410960000` |

## Replay

Keep the repository-relative layout in [inputs.txt](inputs.txt). Use 64-bit Python 3.11+ with NumPy and assertions enabled; the release replay uses CPython 3.12.3 and NumPy 2.3.5. Supply the two original, unchanged D18 files, 3,495,091,968 bytes each. From the repository root, with `PYTHON` pointing to that interpreter and `D18_LOWER` and `D18_UPPER` to the payloads:

```sh
"$PYTHON" -Bu round2/occupation-neon/release/verify.py --lower "$D18_LOWER" --upper "$D18_UPPER"
```

Expected final line:

```text
PASS: combined outward interval [0.565675793084, 0.565772420378]
```

The command verifies the source manifest and both full payload hashes, then re-extracts 2,079 first-body rows, five all-stock floors, and seven young-vertex ceilings per sign. It reproduces the independent crown **root** comparison, checks 114 whole-tree drift identities and 20,298 higher Bellman contexts per sign, and recomputes the all-degree potential ceiling with six charge fixtures per sign. It prints the two exact reward fractions, the seven-bin and improved denominator fractions, and the exact combined endpoint for each sign. No frozen extracted table is required.

## Composition and scope

For each sign use the **original** trial `λ` and compute `λ ± (first_body + higher)/ceiling` in exact rational arithmetic, then round outward to 12 places. Both signs have the same reward support: the `f` charge on nonfounders (degree through 64 and all-stock continuation in the five supported bands), plus the top-owned `hV,hC,P4,Yt,Yb,S3` charges on complete descendant fringes of size at most 12 in the stated interior domain. The first-body affine-floor calculation is an alternative bound on the same reward, not another term. The seven-bin and all-degree ceilings bound the same potential per core; denominator deficit credits are not added to the signed rewards. The inherited occupation prior is `[0.565673570322, 0.565775]`.

The archived base GMP receipts (`round2/attacks/neon1/outputs/neon1-production-rehearsal.tar.gz`, historical reference) at `de931f5` record zero failures in all eight sections and 34 exact tree checks for each original payload. This replay checks those receipts and authenticates the supplied files; it does not repeat the exhaustive GMP verification or rebuild the inherited Astra prior. The analytic transfer and ownership arguments are in the [crown](../../../../../../reviews/cold/astra/occupation-neon/check/README.md), [renewal](../renewal/README.md), [higher-order](../higher/README.md), and [denominator](../denom/README.md) notes. Finite fixtures exercise the implementations, not the limiting transfer theorem.

The [first-body review](../../../../../../reviews/cold/astra/checks/neon-occupation-sol.md) and [higher-order and recharge review](../../../../../../reviews/cold/astra/checks/neon-higher-sol.md) record separate checks of the components. The certified interval books no standalone `g` reward and includes no larger-fringe suffix.
