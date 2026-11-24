> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Original Gold2 upper: disjoint signed H3 frontier

**Composition correction after 835c's receipt:** The original `D_m≤5` older-mother mask below overlaps 835c's `D_m=c_m=2`, `D_mother=80…480` payment. Treat the original two-mask result as a standalone comparison only. Repricing the older-mother mask with `3≤D_m≤5` certifies `9111370309193113579909683/387606272734375000000` ticks/core (23,506.7669182), with matched ceiling `1973454550357085636285062899/22624212005000000000000` (87,227.5485184). The high-stock `D_mother=19…79` mask is unchanged and disjoint. **Only the corrected union is additive with 835c and the frozen H3**: exact floor `13266167036667732636019911/364805903750000000000` ticks/core (36,365.0009506) and ceiling `4064388205245621840498385821/22624212005000000000000` (179,647.7244974). The [floor](middle3-floor.log) and [ceiling](middle3-ceiling.log) are successful replays with `--receiver 120 --rstock 10 --donor 7 --middle-min 3` and the matched `--ceiling`; each checked 176,562 path states and 50,184 mask tests. The remainder of this document records the earlier standalone comparison.

Two new upper-witness masks give an additive floor of **110,095.35685534128 ticks/core** (`1.0013114374973273e−7` potential units/core). A matched rigorous ceiling for their union is **352,029.51530683535 ticks/core** (decimal display). This is a scope barrier: even earning the entire ceiling would leave more than **11.28 million ticks/core** unaccounted for at `241/426` in the existing D24 *separate* denominator/host-savings calculation. It is not a ceiling on H3 or five-core rows outside these masks.

## Exact masks and results

Both price nonfounder mother → middle → child paths, with actual middle `D≤5,c≤4` and actual child `D≤7,c≤3`. Bounds on `D,c` apply at payment. Extra births update live arity and retain the marked path up to its caps. The mother ranges distinguish each new mask from the frozen upper threshold H3 (mother `D=19…79,c≤8`) and from one another. 8790's `D=80…100` seam work is lower-witness only and cannot duplicate either upper payment.

| Upper Gold2 mother mask | Floor (ticks/core) | Ceiling (ticks/core) |
|---|---:|---:|
| `D=80…120,c≤10` | `120607099996138877416354587/1240340072750000000000` | `146836418285767474325320737/565605300125000000000` |
| `D=19…79,c=9…12` | `79742914676261637533783559/6201700363750000000000` | `1045466827444268102106661461/11312106002500000000000` |
| Disjoint sum | `341389207328478012307778247/3100850181875000000000` | `3982195193159617588613076201/11312106002500000000000` |

## Physical ownership and matched ceiling

Use the authenticated private original Gold2 upper payload, SHA `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`, and original outward upper sign `+1` at the nominal `28286843/50000000`. The replay reuses the pinned original [post-FY physical row and prefix registration](../README.md); it opens no lower payload. Each paid row uses `R=[−H3_out]_+` from zero, with original incoming `g→H3`, original outgoing `H3→P4/S3` subtracted, and newborn weight in the original transmitting row. This preserves the existing nonnegative host cap without taking it as a reward. No FY reward is paid. The literal checks cover 37 H3 rows, 40 edge rows, 64 H3 outgoing and 74 g outgoing incidences, four FY-gap identities, and eight distinct actual-path registrations.

The mother and child are live marks; rates use current `1/D`. Remaining outside incidences are unknown predictable reciprocal sums in `[0,n_ext/2]`, entering both the charge and the matching gain hazard. The affine comparison checks exact integer Bellman inequalities after directed rounding of the reward and port coefficients. If the reward changes sign, the maximum of the zero and uncut-affine comparisons is a lower subsolution; for the ceiling, the supremum over ports commutes with the maximum of those two branches. The original exact `LO,HI` bracket gives `LO V_floor(2,1) ≤ π[R] ≤ HI V_ceiling(2,1)`, via one new value per actual nonfounder birth. Degree and arity **caps** are irreversible exits; the lower **gates** (`D=80` or `c=9`) disable reward but leave continuation intact before entry. Thus the zero cuts also apply to the ceiling on these exact scopes, without asserting anything about paid rows beyond them.

## Inward comparison and next direction

The original scalar ledger needs `176655823289357881/10650000000` ticks/core in total to reach `241/426`. In the separately proved D24 denominator and host-savings profile (already includes the frozen H3 and singleton inputs), the additional target credit is `11634970.46294531` ticks/core. Adding the new floors lowers that separate D24 upper bound to **`0.565736172599` outward**, still about `8.47307e−6` above `241/426`. This is a comparison, not a refreshed D24 fixed-trial LP or a board update.

835c's original upper census ranks the H3 `D=180…1998, max c=9…15` band at approximately `1.025/0.936` ppm pooled/held out, and all H3 at `5.499/5.419` ppm. Those are measurements, not bounds. A more promising H3 direction is an all-age mother-tail continuation across the actual `D=88` H3 endpoint knot, with larger mother stocks and a tracked outside donor. However, even the measured total H3 plus quintet negative-row pools are only about 9 potential ppm on those trees, before subtracting already-paid H3. The current D24 separate target asks for about 10.58 additional potential ppm; this favors working on the denominator or additional families alongside any wider H3 effort rather than only widening these compact masks.

## Reproduction

From `.`, run the sibling [`threshold.sha256`](../threshold.sha256) check, then `python3 -Bu round2/occupation-quintet-gold2/post-fy/upper-frontier/frontier.py` with the flags below. All runs were niced to 19 and CPU-pinned, one at a time, as reflected in the logs.

| Flags | Terminal path states / paid mask tests | Wall, exit |
|---|---:|---|
| `--fixtures --receiver 120 --rstock 10 --donor 7` | `176562 / 56457` | `60.20s`, 0; [stdout transcription](floor.log) |
| `--receiver 120 --rstock 10 --donor 7 --ceiling` | `176562 / 56457` | `33.97s`, 0; [captured](ceiling.log) |
| `--old 19 --receiver 79 --rstock 12 --mother-min-stock 9 --donor 7` | `134793 / 37332` | `15.48s`, 0; [captured](stock-floor.log) |
| Same with `--ceiling` | `134793 / 37332` | `46.86s`, 0; [captured](stock-ceiling.log) |
