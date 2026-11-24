# Cold audit: all-age signed `Yb` and disjoint positive floors

**Verdict: correct on the declared mask in the stationary core mean.** Both unused-support floors compose with `ead947b`, giving the outward enclosure **`[0.565678169553103, 0.565770393870602]`**. The paired `θ=1` contrast still straddles zero. All-age coverage applies only to the stated row mask; the audit inherits the physical mean and stable-age identities and the original D18 GMP certificates.

| Stationary charge per core | Whole mask (7,056 labels) | Unused mask (5,865 labels) |
|---|---:|---:|
| Lower-witness positive row | `[4.14056375e−7, 3.615561188e−6]` | `[1.53289608e−7, 1.271173510e−6]` |
| Upper-witness positive row | `[2.27167969e−7, 3.534395283e−6]` | `[1.03263130e−7, 1.334529218e−6]` |
| Shared-trajectory lower minus upper | `[−2.866290553e−6, 3.218086824e−6]` | `[−1.043427602e−6, 1.097198185e−6]` |

## Physical ownership and rate controls

- The accounting unit is a complete ordered embedding `(A,S,M,D)`, at its actual apex, where `S,M` are distinct children of `A` and `D` is a child of `M`. Reversing its physical promotion order removes the most recently filled role at each step, giving exactly one path back to the apex birth. This covers `S,M,D`, `M,S,D` and `M,D,S`. An old tag survives a promotion and emits one separate copy per eligible role. If both branches can be middle, the two assignments are distinct embeddings. The companion checks tag multisets on labelled trees (including the extra side birth in all 12 orders).
- On each of the 7,056 labels the companion reconstructs the row from the raw witness tables. It evaluates the motif's own drift and the two incoming centered terms separately: removing `S` leaves the `hV` derivative at `A`, and removing `D` leaves the `hC` derivative at the actual `M`. Each transfer uses the corresponding removed donor's true reciprocal minus the frozen center. The two signed witness rows are each nonnegative at every exterior-rate corner. For `θ=1`, the affine rows are subtracted **before** bounding with the shared rates; 7,019 row labels admit negative box corners and 7,019 admit positive box corners. These counts are controls on the arithmetic, not claims about attainable stationary signs.
- A tracked exterior child registers only when promoted as a physical child of the governor. It is an extra donor on the continuing tag; on a copy that assigns it to a motif role it instead contributes by that internal edge. Its `1/d` enters the governor's gain and affine reward jointly, and its own gain has the fixed `L+1/D_governor` part plus at most `k/2` from its exterior children. Beyond degree 4 the same child remains an exterior donor with reciprocal at most `1/5`; degree monotonicity forbids re-entry. The labelled controls include 359 live and 183 overflow incidences; 3,024 reward labels have nonzero apex-donor coefficients and at least one apex exterior donor.
- **Founder boundary:** a tracked nonfounder mother of degree `m` has gain at least `1/D_A` and at most `1/D_A+m−3/2`. Indeed, if she has `a` other core children, she has `m−a−2` leaves, and her own parent and the `a` children each send at most `1/2`. The upper bound is not valid for a degree-two founder with one leaf and the apex child. Remove founder and founder-mother births before imposing the bound. Their expected counts are at most `1+t`, because the founder promotes at rate at most one. The paid degree/arity mask uniformly bounds the absolute reward per such apex, so their normalized contribution vanishes even for the signed contrast. The stationary birth normalization remains `μ`.
- Degree and arity never decrease. In particular, once `k` plus unfilled child roles exceeds the cap, none of that tag's extensions can qualify. A killed tag therefore has zero remaining reward on this mask. The companion independently enumerates 126,813 states, checks both rounded Bellman inequalities for each of six objectives, and takes extrema of `μ` times both value endpoints over `μ∈[.565673,.565775]`. The signed endpoint multiplication is necessary. The inequalities hold conditional on the physical history because the rate boxes bound each current intensity and each coefficient uses its own gain rate. Summing over surviving tags and integrating with the discount bounds the complete age integral: the retained number of tags per apex is uniformly bounded, so the discounted terminal value vanishes.

## Disjoint support and exact ledger

For a complete tag let `n=Σk_i−3` and `m=Σ(D_i−k_i)+2n`. The planted apex fringe contains the four cores and their leaf stocks and also `n` disjoint extra child subtrees, each with at least two vertices and one core. Thus its size is at least `m` and its core count at least `4+n`. If `m>20`, both old size reserves are excluded. If `12<m≤20`, the size-12 reserve is excluded and the size-20 suffix is excluded either by `n>0` (at least five cores) or by an included degree greater than 16. The latter can pay with **only four cores**: the bare census `((17,2),(2,0),(2,1),(2,0))` has `m=20,n=0` and fails the suffix degree test. Under the audited degree caps `(6,7,5,7)` this degree exception never fires. The companion exercises the degree seam and the size-12 seam separately.

The frozen higher-family reservations (`occupation-neon/high-parent/README.md#physical-charge-and-disjointness`) therefore cannot use a paid `Yb` row. The high-parent `g/hC` reserve, the post-EDGE joint `f+g` residual (`occupation-neon/interpolate/compose/residual/edge-projection/combined/FREEZE.md`), and the separate `hV/P4/Yt` objectives use different charged rows. The same family argument applies to a price restricted to original `P4/S3` rows: sharing a birth, incidence, or derivative across a centered transfer does not imply sharing the allocated row. This audit adds only the two disjoint `Yb` floors to `ead947b`.

The following are exact. `N_EDGE` is the frozen **net** numerator, reconstructed from each exact EDGE endpoint, its original trial and the full denominator; this retains both negative local margins. The companion checks that adding just the frozen joint price gives the exact `ead947b` fractions. The new numerator adds `Y` once on each side, with `+` on the lower endpoint and `−` on the upper endpoint. Neither denominator changes.

| | Lower | Upper |
|---|---|---|
| EDGE trial | `.565677930574` | `.565770634792` |
| Full denominator | `1357493815758092773/1000000000000000000` | `273436004214812983/200000000000000000` |
| New `Y` floor | `1532896084966910913/10000000000000000000000000` | `1032631300867872471/10000000000000000000000000` |
| Exact endpoint | `137070974112668703714891209168897/242312646112819559980500000000000` | `6574834821628274854124818955857/11621030179129551777500000000000` |

The producer's `[0.565678093598318, 0.565770437762087]` also replays: it uses the EDGE-only base, the displayed rounded floors, and the coarser denominator ceilings. It is not the post-memory composition above.

## Replay

From `[local path]`, with the original read-only payloads:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
P=round2/.work/venv/bin/python
nice -n 15 taskset -c 1 timeout 20s "$P" -Bu pressure-memory/yb-lifetime/fixture.py
nice -n 15 taskset -c 10 timeout 150s "$P" -Bu pressure-memory/yb-lifetime/replay.py \
  --lower b3f18_lo.bin --upper b3f18_hi.bin --degrees 6 7 5 7 \
  --arities 3 1 2 1 --donor 4 --mother 3 --verify-sha
nice -n 15 taskset -c 10 timeout 150s "$P" -Bu round2/checks/yb-lifetime-sol.py \
  b3f18_lo.bin b3f18_hi.bin
```

The producer fixture reported `PASS all 12 chronological birth orders, no duplicate tags, 7008 exact labelled arrow rows`. The authenticated producer reported `rows 7056 disjoint 5865 states 126813 live fifth 100791 exact inequalities 2536260`, and printed the six displayed bounds and EDGE-only candidate. The companion checked the producer source pins, the combined release manifest and both payload digests, and reported `independent Bellman states 126813 row labels {'whole': 7056, 'clean': 5865} nonzero fifth-donor reward rows 3024 signed corner counts {'negative': 7019, 'positive': 7019}` and the exact composed endpoints above. It uses the same pinned binary reader and cross-checks against the production row extractor; the raw reconstruction, rate-box solve and exact composition are separate. No GMP verifier rerun.
