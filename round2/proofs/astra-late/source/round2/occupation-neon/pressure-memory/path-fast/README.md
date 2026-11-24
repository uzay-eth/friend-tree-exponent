> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Native stationary P4 donor lift

This directory prices the same positive complementary P4 decomposition as [`donor_lift.py`](../donor_lift.py) from both authenticated original D18 witnesses. All numbers are stationary per core. Each broader result **replaces** the previous complementary P4 payment throughout; the previously excluded wide band is still excluded. These are bounds on positive original charge components, with no whole-motif payment from a signed row.

| Complementary P4 calculation | Lower-witness floor | Upper-witness floor |
|---|---:|---:|
| Published `(9,7,6,10)`, arities `(3,3,3,2)`, uniform component sum | `3.59297001126e-7` | `2.78525980848e-7` |
| `(16,9,9,16)`, arities `(4,4,4,3)`, uniform component sum | `8.14986586919e-7` | `5.84251449828e-7` |
| Same wide cut, shared physical birth law | `8.84750046e-7` | `6.67298549e-7` |
| Same wide cut, mother component resolved through degree 12, shared physical birth law | **`9.24540133e-7`** | **`7.79798537e-7`** |

Final incremental floors over the published cut are at least `5.65243e-7` and `5.01272e-7`, respectively. Neither final floor nor its increment meets the `1e-6/core` target. The raw wide upper ceilings are about `3.24048e-6` and `2.80956e-6`; this gap is uncertainty, not available certified material.

## Exact computation

- `replay.py` uses `original()` to SHA-256-check both complete blobs before mapping any input. `charges.hh` reads the original mixed `P4` and `hV` potential tables and directional rate arrays. It forms the P4 baseline from its own potential, gain, promotion, and exactly two incoming centered `hV` transfers (the connected tip deletions). It derives the minimum and five exterior directions on first use of each terminal census, memoized once across components. `--row-parity` compares every requested census with the rational extractor, including rejection, all five slopes, their offsets, and the two minima.
- Path prefixes retain old registrations. At a promotion the sum contains the old registration and each eligible fresh extension, as in Python; if both new candidates originate at the same promotion, they are alternatives, not the same physical child filling both roles. Calendar-time seam exit is irreversible. The passive slots and same-clock reward/successor coupling are unchanged. Identical successor-plus-reward arrows are aggregated by addition of their intervals; the sum of independent interval controls fills exactly the resulting interval. Sloped edges stay separate.
- `solver.cc` ranks by path length, donor presence, then mixed radix on `(degree,arity)` digits at each role, donor digits, and mother degrees. Each prefix block filters the available digits by required existing children and arity reserved for unregistered roles. Degree increases, promotions, prefix/donor registration and upgrades to overflow strictly increase rank on every nonterminal arrow; code checks that fact at runtime. Values are allocated in a flat state-indexed slab per component, with one extra shared slab indexed by complete census for charges. For the documented runs native allocated slots equal the reachable counts.
- Rate comparisons retain `lcm(1e6,2,...,Dmax)` exactly; charge arithmetic uses `lcm(1e9,2,...,17)`. The value grid is still **`10^18`**. The signed 128-bit ratio uses explicit floor/ceiling, even for negative numerators. Every stored value, summed successor, arrow endpoint and received coefficient is checked against range limits. Rate unit `≤2^48`, at most 20 arrows at `≤32*unit`, and `|value|<2^60` bound ratio numerators and the direct drift audit well inside signed 128-bit (even using doubled values in differences). In row extraction each original tick has `|tick|<2^50`; gains/differences are bounded by twice that, and even after scaling by `grid/S=10^7` the two transfers and four gain/promotion terms fit signed 128-bit. All four residuals and enclosure order are checked at every live state. FNV-style word checksums include every reached rank and its four values in deterministic completion order.

## Conditional birth-mother interface

`--entries` writes the exact four virgin entry ticks `[lower_lo,lower_hi,upper_lo,upper_hi]` for each component and physical birth-mother degree, in units of `1/grid` **before** multiplying by `μ`. `overflow` is the combined class for degrees above the cap, never a degree-zero mother. Birth does not alter mother degree. `medium-entries.json` and `wide-{head,tail}-entries.json` cover every component; `wide-mother12-entries.json` replaces only the mother reward, on the exact same terminal support, with a finer class partition. The remaining components' degree-three overflow remains valid on every refined sub-bin. The receipts are ready for projective mother constraints from an independent LP; no `embedding-flow` constraint has been assumed here.

`pool.py` sums all six components **at each same birth-mother class** before applying the existing sided birth-law LP. If `V_j(d)` is a certified discounted endpoint, its class constraint prices `x(b V_j(D))`, where `x=π/μ` and `β=x(b·)` on the pre-promotion census. The two endpoint signs are rounded with the physical `μ` bracket and intersected with the shared-class uniform bound. A successor or birth-heavy joint constraint does not replace a live conditional receive-rate bound. This is also a stateless pooling interface: no tagged trajectories or rewards are recomputed when the birth-mother law changes.

## Replays

From the repository root (pick an available CPU; these runs used CPU11 at nice +15):

```sh
R=round2/occupation-neon/pressure-memory/path-fast
SO=<scratch-path>
PY=./round2/.venv/bin/python
L=b3f18_lo.bin U=b3f18_hi.bin
taskset -c 11 nice -n 15 g++ -std=c++20 -O3 -Wall -Wextra -Werror -shared -fPIC "$R/solver.cc" -o "$SO"
replay() {
  OPENBLAS_NUM_THREADS=1 taskset -c 11 nice -n 15 /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M' \
    "$PY" -Bu "$R/replay.py" --library "$SO" --lower "$L" --upper "$U" "$@"
}
replay --python-parity
replay --python-parity --include-old-band --side-cap 0 --mother-cap 0
replay --degrees 9 7 6 10 --arities 3 3 3 2
replay --degrees 16 9 9 16 --arities 4 4 4 3 --budget 40000000 --parts minimum mother 0 1
replay --degrees 16 9 9 16 --arities 4 4 4 3 --budget 40000000 --parts 2 3
replay --degrees 16 9 9 16 --arities 4 4 4 3 --budget 40000000 --parts mother --mother-cap 12
OPENBLAS_NUM_THREADS=1 taskset -c 11 nice -n 15 "$PY" -Bu "$R/pool.py" \
  "$R/wide-head-entries.json" "$R/wide-tail-entries.json" --refine-mother "$R/wide-mother12-entries.json"
```

| Run | Verification result | User + system CPU s | Peak RSS (KiB) |
|---|---|---:|---:|
| Small `--python-parity` | All endpoints, reachable/census counts and row ticks match Python; 216,294 states | `18.12+0.92` | 148,500 |
| Small old-including zero-memory `--python-parity` | All six components match Python, including immediate donor overflow | `7.07+0.79` | 88,520 |
| Medium | Every exact fraction and count agrees with [archived Python](../p4-donor-medium.log); 6,017,604 states and 24,070,416 drift checks | `8.76+0.89` | 231,996 |
| Wide first four components | 76,893,297 states; 307,573,188 drift checks | `82.23+1.97` | 2,240,404 |
| Wide last two components | 64,145,862 states; 256,583,448 drift checks | `69.06+4.14` | 3,051,488 |
| Wide refined mother | 25,476,804 states; 101,907,216 drift checks | `27.23+3.48` | 1,838,852 |
| Wide pooled, degree 3 | Exact dual repair; 462 rows, 3,820 columns | `5.73+0.14` | 95,748 |
| Wide pooled, mother degree 12 | Exact dual repair; 462 rows, 3,820 columns | `5.06+0.08` | 95,996 |

The original medium run used `366.32+12.34` CPU seconds and 1,757,476 KiB RSS. The matched native run is ~39× faster end to end, including ~4.23 CPU seconds to reauthenticate the originals. Compact arrow aggregation reduces the medium rate scan from 77,856,628 arrows to 58,102,093, while preserving every audited endpoint. The initial single wide invocation timed out at 115 seconds; the shard invocations completed. The broad rows have not all been replayed with Python rational extraction.
