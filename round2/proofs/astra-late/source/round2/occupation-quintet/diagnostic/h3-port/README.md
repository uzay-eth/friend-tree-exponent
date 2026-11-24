> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Stationary H3 recovery from the high mother's additional children

On the original `b04cc50` witnesses, explicit sibling registration certifies **3.607803827e−7 lower / 3.833407783e−7 upper** normalized potential units per core. This is the incremental negative-row reward over a sibling-box replacement baseline; the upstream donor's actual reciprocal is unchanged in both terms.

| Same target, per core | Lower | Upper |
|---|---:|---:|
| Stationary floor | 3.607803827e−7 | 3.833407783e−7 |
| Pooled age-20 diagnostic | 9.932912500e−7 | 1.051285496e−6 |
| Held-out age-20 diagnostic | 9.894069172e−7 | 1.015580106e−6 |
| Pooled diagnostic, mother degree 19–37 | 8.198502493e−7 | 8.622612568e−7 |
| Pooled diagnostic, mother degree 38–65 | 1.734410007e−7 | 1.890242387e−7 |

The stationary floor is about 36% of the pooled measurement on each witness. The comparison is descriptive, not an estimate of the remaining stationary gap.

## Row and disjoint accounting

For the nonfounder path `h → m → y`, use mother `D=19..65,c=3..7`, middle `D≤5,c=2..4`, young tip `D≤18,c≤3`. Arity counts core neighbors, including the parent. There are `n=c_h−2` other children of `h`. `Reader.h3_physical` returns the signed charged affine row in ticks. Put `b=b_h` and `P_box=P` with the sibling part of the exterior rate replaced by its adverse box corner (each summand `1/D_s` replaced by `1/2` if `b>0`, or zero otherwise). Leave upstream and every other rate unchanged. For each actual sibling define `z_s=max(0,b/2)−b/D_s≥0`. Thus `P=P_box−Σz_s` and the paid increment is

```
Δ = [−P]_+ − [−P_box]_+ = F(P_box, Σz_s),
F(B,z) = [z−B]_+ − [−B]_+.
```

`F` is nonincreasing in `B`, nondecreasing and convex in `z≥0`, with `F(B,0)=0`. Therefore `Σ F(B,z_s) ≤ Δ ≤ Σ F(B,n z_s)/n`. For the floor use the exact maximum of `P_box` over the still-untracked upstream, middle and young-tip reciprocal boxes. This yields a separately priceable nonnegative term for each distinguished sibling, requiring no independence between siblings and no conditioning on empirical means. `price.py` also implements the Jensen majorant with favorable `B` for an optional ceiling; the reported result is the floor.

The four-vertex donor registration has a fork at `h`, but always pays the three-vertex chain H3's high tip. An H3 whose center is the genealogical mother is a different crown row; it is not included. A covered port has zero physical slope and no increment. The paid term overlaps the original full physical negative-row reward and can only be added to a sibling-box baseline. The `8d8f` baseline on 38–65, including the disjoint `c_h=2` direct slice, is aligned: `1_submask [−P_box]_+ + Δ ≤ [−P]_+` everywhere on the union. The resulting joint floor is **5.88485e−7 lower / 6.56968e−7 upper**, with exact joint fractions in `../../higher/full-h3/README.md`. This does not update an exponent endpoint.

## Stationary verification

- Root at a nonfounder birth and extend prefixes in each feasible chronological order `hmys`, `hmsy`, `hsmy`. At the designated parent's promotion, the shorter prefix continues at its actual new arity and one new longer prefix starts at that arity plus the newborn `(2,1)`. Other promotions change only the existing arity. The actual parent pointers determine every reported `(h,m,y,s)` once. No rank-by-degree or inferred birth order is used.
- The distinguished sibling is exact at degree 2; at degree ≥3 it keeps its identity and uses reciprocal `[0,1/3]`. Degree and promotion of this overflow donor leave its abstract state unchanged. Other participants exit only on monotone degree/arity upper-cap violations. The mother degree lower threshold gates reward, not occupation. Additional children outside the registered roles are retained in the rate boxes; in particular, earlier siblings are not killed by another birth.
- Each represented gain hazard is `D−c +` the sum of visible reciprocals plus residual rates in `[0,(c−visible_neighbors)/2]`. The fixed promotion hazard is `(D−c)/D`. External events do not change represented coordinates. The overflow reciprocal is adversarially boxed for both reward and hazard. Relaxing its correlation only weakens the lower bound. At each state the integer verifier checks the exact minimum residual of `reward + L V − HI V` over all hazard boxes. Rewards and values are rounded down on the `2^20` subgrid of witness ticks (`2^40` per potential unit). The sum of nonnegative bounded prefix potentials gives, by normalized stable-age balance, `πΔ ≥ LO·Σ_order V_birth`, using the original exact Gold1 endpoints for the birth/discount interval as in `../../higher/README.md`. No snapshot contributes to this inequality.

## Reproduction and evidence

Run from Astra root on one niced core:

```sh
R=round2/occupation-quintet/diagnostic/h3-port
taskset -c 1 nice -n 19 timeout 110s python3 -IBu "$R/price.py" --side lower --cut 2
taskset -c 1 nice -n 19 timeout 110s python3 -IBu "$R/price.py" --side upper --cut 2
bash "$R/measure.sh"
```

Both floor runs exited 0: lower 70.82 s wall / 46.05 s user / 301392 KiB RSS; upper 79.00 s / 46.61 s / 299116 KiB. The ancestry oracle found 26 distinct registrations versus brute force. The shared terminal table had 411,750 states per side, 114,766 / 202,986 with a positive certified per-donor instantaneous reward. Exact normalized floors:

```
6491591975220984267131640711/17993195545502789105918935040000000
2207208187097745408119311435947/5757822574560892513894059212800000000
```

The diagnostic replay exited 0 (lower/upper 1.35/1.06 s wall); it found 40,542 pooled and 10,823 held-out H3 rows across 896,887 and 235,709 cores. Its two clips reconcile with the physical negative row, and the one-sibling/Jensen bounds bracket the exact empirical increment. Files `sample-lower.tsv` and `sample-upper.tsv` retain per-tree totals for the two mother-degree bands. The preserved blobs are read-only and inherit their original witness authentication.
