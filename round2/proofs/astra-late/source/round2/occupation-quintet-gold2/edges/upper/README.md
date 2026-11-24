> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# Gold2 upper: original edge and small-payer FY budget

The two disjoint physical floors total **`4.07124156e−7/core`** (rounded down) on the original Gold2 upper witness. The tested cuts did not exceed `1e−6/core`. The all-arity ceiling for the same small-payer FY mask is `1.95323307e−6/core` (rounded up), so a better certificate could still reach the target. This is a stationary budget result for the upper witness, separate from the frozen lower `99/175` calculation.

| Payment or bound | Exact normalized value | Approx. per core |
|---|---|---:|
| Full original edge slack, both endpoints degree ≤18 | `1190404372635664108984647/5497558138880000000000000000000000` | `2.16533294e−10` |
| FY floor, all donor degrees, absorb at 15 | `4473996630151725250641/10995116277760000000000000000` | `4.06907623e−7` |
| FY ceiling on the same physical mask, absorb at 8 | `107380123330368644331237/54975581388800000000000000000` | `1.95323306e−6` |

The edge floor recovers **238.0808745 ticks/core**, about 51.4% of the original cap of 463 ticks. The FY floor covers `18<D_host≤251`, all exact `G2` host arities through 22, excluding `M3`, paying neighbor degree ≤5, and any distinct physical donor. Since the host is above 18, the paying edge is outside the full-edge mask. The complementary `payer/upper` calculation uses paying-neighbor degree >5. No higher rebate or Jensen term is included.

## Physical checks and threshold

The outward sign is `+1` for this witness. Every FY reward is the maximum over the **original 20 signed knots** minus the signed score at the physical donor reciprocal. At the absorbing donor threshold `k`, physical degree remains at least `k`; the bound takes the exact extrema on `[0,1/k]` using the knots and the band endpoints. Gains and promotions before entry are tracked from actual physical births. Within the band, internal transitions vanish from the test function; live receive rates use the band's reciprocal interval. The floor omits promotion credit from a prospective mother already in the band. The ceiling retains that mother's possible source, using promotion hazard at most one. All other transitions and exact affine checks use the original Gold2 discount interval `[.565713486333,.565736986679]`. The `--ceiling` mode bounds the upper witness's reward, rather than switching witness signs.

The [floor at donor cut 8](absorb8.log) is `3.83940781e−7/core` in 28.116 CPU seconds with host arity ≤17. [Cut 12](absorb12.log) gives `4.00675592e−7/core` in 53.896 seconds. [Cut 15 with full arity](absorb15.log) gives `4.06907623e−7/core` in 87.297 seconds. These are alternative bounds on overlapping support. The tested widening yields only about 6% over cut 8, so it supplies no evidence that further refinement of the donor cutoff alone will reach `1e−6`. The [ceiling and edge receipt](ceiling8.log) also ends in `PASS`. The replay checks 54 physical charged edge-plus-FY cases, including both cutoff seams, and every computed Bellman inequality.

## Replay

From the project root, restore the authenticated Gold2 private inputs under `round2/.work/replay-inputs/` as described in [the parent edge README](../README.md). The replay binds the published receipt and original upper payload `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876`. Add `--hash` to reauthenticate the binary. No production verifier rescan is needed.

```sh
R=round2/occupation-quintet-gold2/edges/upper/replay.py
taskset -c 5 nice -n 19 python3 -Bu "$R" --donor 15 --stock 22 --fy-only
taskset -c 5 nice -n 19 python3 -Bu "$R" --donor 8 --stock 22 --ceiling
```

The second command also prints the degree-18 floor; use `--edge-only` to replay it by itself. Only the completed receipts linked above are used.
