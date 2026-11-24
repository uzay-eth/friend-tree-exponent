> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# D20 older crowns: two real arms and the same-edge excess

**On the tested cuts, the multisource `hV` side-edge certificate exceeds the finite-geometry widening certificate on both witnesses.** The size-disjoint older reward is at least **`1.10972285509e−7 / 3.25225638090e−7` per core**. Widening the finite solve only to size 18 instead gives `4.622430884e−8 / 5.082444891e−8`. These are overlapping alternative extensions, not additive. No credit is written into the audited renewal or receipts.

## Physical birth sources and exact ownership

The tags are an apex `A`, a spine child `B`, its child `C`, and a distinct apex side child `U`. They represent two actual arms around the same hub; the side arm is not a prescribed rate or a reset to a virgin at the time the spine completes. Start at the physical Palm birth of `A`. Use all three birth orders `(B,C,U)`, `(B,U,C)`, `(U,B,C)`. A tagged parent's actual promotion continues the old registration and, where it makes the next prescribed adjacent tag, emits exactly one extended registration with the post-promotion census and newborn `(1,0,1)`. In particular an already-grown side sibling survives the introduction of the spine tags, and vice versa. The branching of registrations is linear in expectation and does not duplicate the final `(A,B,C,U)` incidence. Untracked branches remain physically present throughout.

All known sends across these edges are exact in **both directions** at the same instantaneous pair of degrees. An untracked port retains its normal receive box; the terminal paid quantity does not treat it as a known donor. Each intrinsic gain advances one degree, and each promotion advances one arity, so the finite tagged system is acyclic. At an extension birth, `B/F` is the compensator of the actual promotion followed by that next registration, not a mass estimated from selected snapshots. Stop the relevant registration on first overshoot of a tag cap. The resulting value is a lower bound for every physical history, using the upper discount `.565768` and the lower birth-rate normalization `.565680`. No D18 occupation numbers occur in that proof.

Write the original charged `hV` row as `a+Σ_i r_i c_i`, with the actual multiplicity of each free port. Its independent GMP corner floor is `a+Σ_i n_i min(0,c_i/2)≥0`. For the one chosen apex side edge, pay just `c_A/D_U−min(0,c_A/2)`. Summing over all admissible `U` counts at most the full sum of distinct side-edge increments in that row. The row floor, other donor increments, and other motif families are not included. In particular the floor diagnostic in the script is **per tagged incidence** and cannot be paid once per donor as row credit. The signed D20 corner allocation and actual shared degree are used on the same edge; no pressure is reassigned between endpoints. The elbow case `g` has zero such free-port excess on the examined small cut, because the side receive is routed into the higher charge.

The paid support requires `4+Σ_tag L≥17`. Those leaves are distinct and are present in the true planted fringe; surplus branches only increase the true size. Thus every paid `hV` incidence is **strictly outside** the audited size-≤16 `g/h` support. The audited first-body and all-stock payments are different original charged families. We do not infer disjointness from the fact that a donor was separately tagged. The independently widened size-17–18, ≤4-core experiment overlaps this test, compares the **whole** charged `g/h` sum on that finite shell, and is not added to it. The two methods have different larger supports; the comparison does not prove which would win with unlimited enumeration.

## Results and source diagnostics

On the main cut degrees `D_A≤16`, `D_B,D_C,D_U≤6`, arities `k_A≤3` and other `k≤2`, the separate birth-order prices (rounded down) are:

| Source order | Lower witness | Upper witness |
|---|---:|---:|
| `B,C,U` | `4.03428774415e−8` | `1.75177464202e−7` |
| `B,U,C` | `2.76609076108e−8` | `5.90407013113e−8` |
| `U,B,C` | `4.29685004568e−8` | `9.10074725771e−8` |

The summed lower certificate is about **2.40× / 6.40×** the size-18 widening lower certificate. On the smaller cut `D_A≤13`, other `D≤4`, the excess still pays `1.61598344553e−8 / 5.08973748643e−8`; the discounted normalized occupation of **positive-excess tagged incidences** is at least `.0004864795880006 / .0005899645892252`. On that cut the per-incidence repeated-floor diagnostic pays only `2.460813780736e−12 / 2.459179171753e−12`. These independently minimized objectives are not a ratio guarantee under the unknown physical law.

Pointwise controls, using `(L,k,h)` tuples in order `(A,B,C,U)`, show why measuring sources matters:

| Witness and tagged census | D18 frozen corner floor | D20 frozen corner floor | D20 single-edge excess |
|---|---:|---:|---:|
| Lower: `(8,3,1),(0,1,1),(3,0,1),(3,0,1)` | `2.92318e−9` | `3.05397e−9` | `2.0639987e−4` |
| Upper: `(9,2,1),(0,1,1),(3,0,1),(1,0,1)` | `2.99095e−9` | `3.02747e−9` | `4.3454967e−4` |
| Lower: `(11,3,1),(0,1,1),(5,0,1),(5,0,1)` | `3.02251e−9` | `2.90797e−9` | `2.2021954e−4` |
| Upper: `(11,3,1),(0,1,1),(5,0,1),(1,0,1)` | `2.90699e−9` | `3.04969e−9` | `3.6860609e−4` |

The third census has positive old-donor loss but replacing the last donor by a virgin makes the lower excess zero. The fourth has the opposite send preference. D18 here is only a pointwise diagnostic at its own signed trial. Large local multipliers are not assumed to be stationary rewards; the birth-order calculations explicitly pay for reaching the support. With the audited D20 denominator ceilings the new reward scale is roughly `8.18e−8 / 2.38e−7` in endpoint units, conditional on ledger review. The reviewed headline is left intact.

## Replay

```sh
P=./round2/.venv/bin/python
D=round2/occupation-neon20/renewal/older
taskset -c 8 nice -n 17 timeout 70s env OPENBLAS_NUM_THREADS=1 "$P" -Bu "$D/replay.py" b3f20_lo.bin lower --preauthenticated --family hV
taskset -c 8 nice -n 17 timeout 70s env OPENBLAS_NUM_THREADS=1 "$P" -Bu "$D/replay.py" b3f20_hi.bin upper --preauthenticated --family hV
taskset -c 8 nice -n 17 timeout 48s env OPENBLAS_NUM_THREADS=1 "$P" -Bu "$D/compare.py" lower
taskset -c 8 nice -n 17 timeout 48s env OPENBLAS_NUM_THREADS=1 "$P" -Bu "$D/compare.py" upper
```

All four runs completed. The final main replay checks the cross-multiplied Bellman inequality against all independent unknown-port corners by summing their exact minima. The smaller source/floor probes ran with `--cap 13 --tip 4 --objective excess --objective positive-source --objective floor`. Only preauthenticated local payloads were used for these quick diagnostics; `replay.py` can reauthenticate by omitting the flag. No audited file was edited.
