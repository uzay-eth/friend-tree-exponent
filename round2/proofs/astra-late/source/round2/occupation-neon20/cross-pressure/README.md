> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# D20 versus D18 on one physical occupation

On the common all-age `f+g+hV` mask with maximum degrees **16/5/4**, D20's charged pressure is **9.42–9.57% lower** on the lower witnesses and **9.41–9.83% higher** on the upper witnesses (rounded outward). The sum of the two charged pressures drops by **0.93–0.98%**. These conclusions concern the retained rows. The omitted complement leaves the sign of the whole-pressure change unresolved for either pair.

## Matched comparison

Differences below are D20 minus D18 per physical core, in units of `1e−7`; decimal intervals are rounded outward. Both models price the identical mask, with 135 singleton, 140 edge and 324 rooted `hV` pattern classes. “Independent” means subtracting the two marginal enclosures (or adding the two side contrasts in the last row). “Box” still uses the same occupation variables for the two objectives; it boxes the exterior reciprocal rates without shared donor allocation. “Shared” allocates those rates through order three, using order-four extensions capped at degree four.

| Contrast | Box, direct | Shared, independent subtraction | Shared, direct |
|---|---:|---:|---:|
| Lower pair | `[-2.36940, -2.10990]` | `[-3.52414, -0.95136]` | `[-2.30239, -2.17309]` |
| Upper pair | `[1.72203, 1.93391]` | `[0.68534, 2.96756]` | `[1.77368, 1.87921]` |
| Sum of lower and upper changes | `[-0.43550, -0.38786]` | `[-0.52870, -0.29389]` | `[-0.42319, -0.39941]` |

Shared allocation cuts the width of each direct interval by about half relative to the matched box. The shared ratio enclosures from [final.log](final.log) are `[0.904335, 0.905714]` (lower), `[1.094101, 1.098291]` (upper), and `[0.990263, 0.990684]` (sum of D20 pressures divided by sum of D18 pressures). The box-model independent subtraction straddles zero for both pairs, while direct pricing certifies both signs before adding shared constraints. The [uniform degree-four pilot](pilot.log) prices all three families separately: its `g` and `hV` contrasts are of order `1e−12` or less, while the `f` contrasts are of order `1e−7`. The [graded 7/4/3 repeat](wider.log) also finds a negative lower contrast, a positive upper contrast and a negative four-way sum.

## Omitted pressure and proof scope

| Pair | Omitted D20−D18 contribution | Whole D20−D18 (using common `μ`) | Whole D20/D18 ratio |
|---|---:|---:|---:|
| Lower | `[-1.24391e−4, 1.15104e−4]` | `[-1.02170e−4, 9.06807e−5]` | `[0.04117, 4.42434]` |
| Upper | `[-1.25701e−4, 1.15778e−4]` | `[-1.02630e−4, 9.08619e−5]` | `[0.04018, 4.37440]` |

The replay authenticates the four original binaries and checks 2,079 exported D20 singleton rows per sign. It uses the original charged `f`, `g`, and rooted `hV` rows on the same physical occurrences: `f` includes its assigned newborn-edge term, and `g` uses its separate charged row. It checks exact nonnegativity on the retained box and 16 whole-charge identities against direct labelled-arrow generator values, including boundary-stock fixtures. The nonnegative omitted complement relies on the original all-state signed-row certificates (D18 [GMP archive](../../../../../../reviews/cold/astra/attacks/neon1/REHEARSAL.md), D20 [occupation audit](../../../../../../reviews/cold/astra/checks/neon20-occupation-sol.md)). The D20 JSON exports and extraction source are SHA-pinned; their recorded source revision is metadata, not a local checkout.

For a sign `s`, whole pressure is `s(μ−λ)πΦ` with the common discount in the [D20 renewal bracket](../renewal/README.md). The D18 stationary lower comparisons and [upper ceilings](../../occupation-neon/denom/README.md) are inherited; D20 uses the authenticated binary `c0` and the exported recharge ceilings. Each omitted side is bounded by its whole ceiling minus its certified masked floor. The sharper whole contrast evaluates the common-discount extrema with both the potential and masked lower bounds, including the kink where their maximum switches. The common positive conversion from `η=π/μ` cancels in masked ratio tests of `C20−t C18`. Each reported LP endpoint uses exact residual repair over all columns and certified group caps; solver optimality is not required. The mask is part of the original charged sum and provides no additive release credit.

## Replay

From the Astra root, with the four authenticated payloads at `b3f{18,20}_{lo,hi}.bin` (or supply the four `--d18-lo`, `--d18-hi`, `--d20-lo`, `--d20-hi` paths):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 5 nice -n 19 timeout 165s \
  round2/.work/venv/bin/python -Bu round2/occupation-neon20/cross-pressure/replay.py \
  --caps 16:5:4:4 --share 0 3 --seconds 8 --ratios --joint-only --four-way
```

The pinned replay printed `PASS CPU/wall 103.208 112.654`; [the log](final.log) records the exact rational bounds. Per-objective time limits can affect the multipliers found and hence the tightness of the repaired bounds.
