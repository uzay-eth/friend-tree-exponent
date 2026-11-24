> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../README.md#claim-28) for exact scope and citations.

# Gold2 lower: small-paying-edge FY

The strongest completed floor is **`3.594882105691855e−7/core`** (approximate), exact `39526146756920663734311/109951162777600000000000000000`. The mask is exact-G2, non-M3 host, `18<D_host≤251`, host arity `≤17`, paying neighbor degree `≤5`, distinct donor of any degree. This is a lower-side Gold2-only result, conditional on acceptance of the original Gold2 verifier and its discount bracket. Gold1 supplies only unchanged solver and row-reading code, no numerical floors or discount inputs. The wider Gold2 M3 cutoff (`80/10`) is used. This small-payer lane alone does not meet the threshold needed to exclude `μ=99/175`.

| Variant | Normalized floor, approximate | CPU |
|---|---:|---:|
| Host `120/8`, payer/donor degree ≤5 | `5.68745438e−8` | 2.202 s |
| Host `251/13`, payer/donor degree ≤5 | `9.24038351e−8` | 10.711 s |
| Host `251/12`, payer ≤5, donor degree/arity ≤20/5 | `1.97405701e−7` | 36.405 s |
| Host `251/14`, donor absorbs at degree 8 | `3.25262141e−7` | 17.015 s |
| Host `251/17`, donor absorbs at degree 12 | `3.59488211e−7` | 45.643 s |

Use only the best of these overlap-containing payments. All price the physical ordered FY knot gap recomputed with the original **20** Gold2 knots. A donor evolves under actual promotion and gain from registration. In the absorbing state its degree and stock can keep growing; the tracked band never exits and the internal transitions cancel. Its true reciprocal lies in `[0,1/12]` in the strongest run. The reward floor uses all knot endpoints in this band, and tracked gain hazards use the matching receive interval. If the untracked prospective mother already occupies the absorbing band at host birth, its additional birth source is omitted in the lower solution. Both mother/child and the two sibling orders are represented. All Bellman corner inequalities use exact integer rounding and checks.

The separate full Gold2 edge residual on `max(D_u,D_v)≤18` is **233.95946818135405 ticks/core** (approximate), exact `2339594681813540639638101/10000000000000000000000`, or `2.127848967405535e−10/core` after dividing by `2^40`. It is disjoint from this FY floor, and the complementary `payer>5` stream can retain `max(D_host,D_payer)>18` without overlap. No higher residual or `M3` transfer is included.

## Replay and pins

Restore the private Gold2 files under `round2/.work/replay-inputs`, including the detached source with `certificate_gold2.json`. The source commit is `b435b4c824e2ba5835f58b85e67a27de7d6a0d49`; lower binary SHA-256 `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413`. The replay checks header, certificate binding and file stat stability; `--hash` additionally rehashes the private copy. No production verifier scan is performed.

```sh
R=round2/occupation-quintet-gold2/edges/fy.py
taskset -c 6 nice -n 19 python3 -Bu "$R" --host 251 --stock 17 --donor 12 --absorb
taskset -c 6 nice -n 19 python3 -Bu "$R" --edge18
```

Completed receipts: [absorbing cut 12](absorb12.log), [absorbing cut 8](absorb8.log), [finite mixed cut](mixed20.log), [small-wide](small-wide.log), [small pilot](small.log), and [degree-18 edge](edge18.log). Each ends with `PASS`.
