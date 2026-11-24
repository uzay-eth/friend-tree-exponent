> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Exact FY slack from actual neighbor births and ages

The completed `175/5/8` run certifies the following stationary lower bounds, in normalized units per core. The upper witness exceeds the `1e−6` target; the lower witness does not.

| Paying edge / other donor | Lower | Upper |
|---|---:|---:|
| Mother pays / child donates | `2.31931793e−8` | `1.00789711e−7` |
| Child pays / mother donates | `2.14279671e−8` | `6.52863970e−8` |
| Two children, both ordered directions | `4.39574356e−7` | `1.53773836e−6` |
| **Total** | **`4.84195503e−7`** | **`1.70381447e−6`** |

Rounded down. The exact totals in ticks per core are `266189292981106865228871/500000000000000000` and `1873363823474418533418231/1000000000000000000`; divide by `2^40`. The `90/5/5`, `175/5/5` and `175/5/8` runs use nested masks, so only the strongest matched receipt is used. The large lower-only attempts did not finish.

## What is charged

For an exact `G2` host `x` outside `M3` with `c≥2`, an ordered receiver `y` and different donor `w`, put `z=g(x^+,y)-g(x,y)` with `x^+` the host after a degree gain, `t=1/D_w`, and `F(u)=σ(u z−A_x(u))`. The original outward signed exterior contribution at the host is `σ(c−1)A_x(1/D_y)+Σ_w F(t_w)`. The reward is `max_k F(ξ_k)−F(t)`; interpolation within each physical donor knot segment is exact. For a tail receiver the slope is its physical interpolated slope; the maximum of affine functions in that slope is convex, so the original verifier's corner maxima dominate the physical maximum. No additional Jensen reward is taken. Replacing physical FY terms by their maxima in the edge proof leaves `G+Σ FY≤cap_g`; all `M3` and higher transfers are unchanged. The FY reward covers ordered pairs, while the paying edge's underlying row is still counted once. Both endpoint-host contributions of an edge can be included because the original verifier bounds the two endpoint contributions separately.

The earlier [degree-18 payment](../README.md) prices the **full** `cap_g−G` on edges whose endpoints both have degree at most 18. We explicitly found positive overlap (28 lower and 26 upper fixture triples), e.g. the exact arity-seam host `(16,16)` with paying neighbor `(2,1)`. Every FY payment here instead requires `max(D_x,D_y)>18`, irrespective of the other donor’s degree. Thus FY adds to degree-18 only on its complement, as well as to the centered singleton ledger. It uses no `H3` or higher cap residual, including no `M3` `env23` or higher `env34p` envelope slack.

## Physical occupation

Each path is registered at the actual promotion that creates its second edge. For a mother–child path, the mother edge is registered at the host's birth; each subsequent host promotion creates one pair with that old mother edge. For two children, the first edge is registered when its child is promoted, then each later host promotion creates one unordered pair with that older child; its two ordered FY rewards are evaluated separately and added. Every pair of physical children has a unique older member even if their census labels are equal. In both cases the already-tracked neighbor evolves from its actual first promotion state, and the newly added neighbor starts at `J=(2,1)`.

The three stopped Bellman stages are singleton, old edge, and completed path. The degree caps are 175 at the host and 5 on both other vertices, and all three have arity cap 8. Each vertex uses the reciprocals of currently tracked neighbors exactly; its remaining receive hazard is bounded by `[0,(c−tracked)/2]`. This box is only for transition hazards: every nonzero reward uses the actual path state and actual reciprocal, not a worst-case donor value in the reward. Paths stop on leaving the monotone degree/arity cuts. Discount at the proven high endpoint `HI`; compute and round down each continuation and verify the exact worst affine corner at that rounded value. When summing over the physical populations, the edge/path birth terms are equal and opposite to the preceding stage’s output terms. Consequently the combined subsolution has only the source at core birth. Stable age supplies `μ W(J)≥LO W(J)` per core. This books the actual promotion times and all subsequent ages inside the stopped support without assuming independent ages or neighbor degrees.

## Conditional gap versus box width

Examples in **million ticks per ordered pair**, with actual donor degrees substituted into the physical reciprocal. The box width is the difference of the maximum and minimum over all FY knots at that same host and receiver, not a stationary expectation. `(D,c)` denotes the census state.

| Host; paying neighbor | Side | Donor `D=2` | `D=4` | `D=6` | `D=18` | Box width |
|---|---|---:|---:|---:|---:|---:|
| `(38,2); (2,1)` | Lower | 0 | 0.588 | 2.817 | 13.298 | 54.958 |
| `(38,2); (2,1)` | Upper | 53.067 | 18.348 | 7.815 | 1.343 | 53.067 |
| `(38,2); (4,1)` | Lower | 40.743 | 11.042 | 3.175 | 0.195 | 40.743 |
| `(38,2); (4,1)` | Upper | 4.494 | 0.353 | 0.013 | 7.131 | 45.137 |

See [stock8.log](stock8.log) for more conditional examples, tail receivers, exact totals and all three component floors. The independent cold review reported exact agreement on the full degree-18 edge solve, the three small-cut FY roots and 96 selected strong-cut path/edge values per sign, and verified unique birth ownership on increasing trees. The strong-cut totals are bound to the producer's receipt.

## Full replay from a clean archive

From the unpacked project root, keep `round2/occupation-quintet/census/census.py`, `edges/replay.py` and `edges/fy/` in their original layout. Restore the authenticated `quintet_lower.bin`, `quintet_upper.bin` and `gold.json` under `round2/.work/replay-inputs/`. The large binaries are inputs, not generated by this replay. Then run:

```sh
CPU=5 sh round2/occupation-quintet/edges/fy/reproduce.sh --hash
```

This runs the degree-18 edge baseline and then both sides of the strongest matched FY cut, serially on one niced core. Override `CPU` if core 5 is unavailable. `--hash` checks each file against its pinned digest and can be omitted when replaying on the same authenticated private copies. The only full-file read with `--hash` is authentication; there is no production row rescan or solver dependency.
