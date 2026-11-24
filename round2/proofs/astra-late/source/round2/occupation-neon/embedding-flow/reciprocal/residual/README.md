> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Reciprocal moments for the frozen post-EDGE residual

On the largest matched cut tested (`H=9`, `R_f` on `D≤8`, `R_g` on
both endpoints `D≤6`), reciprocal moments certify **`1.45443e−7` lower /
`1.08338e−7` upper per core**. These are prices of the corrected
post-EDGE residual. The masks differ from the frozen `ead947b` payment,
so the comparison does not produce a release endpoint.

## Matched results

All entries are joint `R_f+R_g` stationary floors per core; each cell is
`lower witness / upper witness`. Read exact fractions from `matched.log`.
All runs here use the original discount interval `[.565673570322,.565775]`,
tail arity cap four, and one Bernstein slab.

| Census and paid cut `H:F:G` | Degree-only boxes | Degree-only, shared donors | Shared + reciprocal moments (order 1) |
|---|---:|---:|---:|
| `5:4:4` | `8.953e−19 / 1.228e−18` | `6.51750e−8 / 3.74527e−8` | `8.27432e−8 / 4.80707e−8` |
| `7:6:4` | `1.041e−18 / 1.512e−18` | `9.58596e−8 / 6.55017e−8` | `1.08005e−7 / 7.40062e−8` |
| `9:8:6` | `5.08979e−10 / 1.922e−18` | `1.33430e−7 / 9.90912e−8` | `1.45443e−7 / 1.08338e−7` |

On `9:8:6`, moments improve over degree-only **shared** donor rows by
9.0% / 9.3%; most of the difference from the plain boxes already comes
from sharing. A separate smaller-cut diagnostic (`small.log`, explicitly
using the tighter `ead947b` discount interval) also ran order-one with
sharing off: its joint floors were `9.73e−19 / 1.26e−18` at `5:4:4` and
`1.85e−9 / 1.55e−18` at `7:6:4`. Thus the gain on these cuts depends
strongly on matching the donor projections to the moments. `small.log`
also records separately priced `R_f` and `R_g` on that interval.

The existing `ead947b` rewards are `1.03095824e−7 / 5.99909485e−8`, on
`F=18,G=4`. The `7:6:4` floors are already 104.8% / 123.4% of these;
the `9:8:6` floors are 141.1% / 180.6%. The latter masks **are not
nested** with `F=18,G=4`: they omit first-body degrees 9–18 and include
edge degrees 5–6. The original frozen denominator ceilings are
`1.357493815758092773 / 1.367180021074064915`; using the exact ceilings,
the gross `9:8:6` floors divided by the ceilings are `1.07141e−7 /
7.92422e−8`. The positive differences from the already-booked rewards
are comparisons only. Do not add those floors on top of `ead947b`.

## What is priced and checked

* Verify the frozen EDGE, composition and pressure-memory manifests and
  authenticate both full D18 payloads. Recompute `V`, `z`, the census
  margin, retained charges, and the exact frozen denominator. Check the
  exact EDGE base and the outward side of the frozen subsequent endpoint.
* Call the frozen `freward` and `greward`. The root intercept retains the
  census correction, deficit, and newborn continuation. The old-edge
  intercept subtracts the **full** booked inner-edge floor and includes
  only surviving-edge continuation. Evaluate `R_f` on total root receive;
  evaluate `R_g` on the two exterior receives (mother plus sibling pressure
  at the parent, child pressure at the child). The internal reciprocal
  sends are already included in the old-edge intercept.
* Pay only named inner root and existing-edge cells. All exterior degrees
  still enter the constraints via the exact small census and the positive
  Bernstein tail representation in the base reciprocal model. Nothing
  pays a mixed `g` row, high-degree type, or higher-body family. No tail
  payoff extrapolation is used.
* The exact coefficient maps agree by column name across all three
  variants. The first model also evaluates the frozen objective at every
  exterior box corner of every explicit root/edge class, checking affine
  coefficient recovery, nonnegativity, and zero rewards outside the mask:
  300, 1,820 and 6,300 checks **per witness**, respectively.
* HiGHS proposes prices only. The bounded path uses the separate pricing
  module's rational reduced-cost repair on **every** original column with
  the model's certified group caps. The largest moments LP has 13,333
  columns, 1,665 equalities and 7,573 inequalities. Its final unnormalized
  repairs are `1.295e−16 / 1.634e−17`. Stable age converts the nonnegative
  discount-integral floor with the original lower discount multiplier.

The next release comparison needs exact `F=18,G=4` support, or a separately
reviewed full-mask construction that treats the remaining inner first-body
labels and any validated piecewise tail charges. The present results are a
matched constraint study; no floor is added to or substituted in the release.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 10 nice -n 19 timeout 260s \
  round2/.work/venv/bin/python -Bu \
  round2/occupation-neon/embedding-flow/reciprocal/residual/replay.py \
  --lower b3f18_lo.bin --upper b3f18_hi.bin \
  --cuts 5:4:4 7:6:4 9:8:6 --variants box share moments --seconds 28 --interval original
```

`matched.log`: both witnesses, all nine matched variants, exit 0;
181.334 CPU / 235.46 wall seconds including frozen extraction. The largest
moment pair took 10.13 / 14.62 CPU seconds for pricing, respectively.
The earlier `inner.log` outer timeout and `inner6.log` 15-second upper IPM
timeout are superseded by this complete replay. `--audit-only` checks
the objectives and frozen extraction without running a numerical solver.
