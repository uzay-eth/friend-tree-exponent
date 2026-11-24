# All-age physical `hV` occupation

For every top-owned nonfounder `hV` path with all three degrees at most 7 and arities at most 2, the two signed witness rewards and their **joint** contrast are strictly positive in stationary expectation. Keeping the same degree-two side children from birth raises the certified lower-minus-upper margin from **`6.07e−14` to `1.198e−13` per core**. This is a lifetime calculation over all matching descriptors, including a virgin last slot. There is no frozen-row or isolated-entry approximation.

All entries below are intervals in original signed numerator units per core, scaled by `10⁻¹⁶`, with outward endpoints. The mother is untracked (`1/D_mother∈[0,1/2]`) in these comparisons.

| Mask and side registry | States | Lower witness | Upper witness | Lower − upper |
|---|---:|---:|---:|---:|
| `D≤7,k≤2`, free rates | 5,448 | `[1918851,2907795]` | `[1916498,2903164]` | `[607,6376]` |
| `D≤7,k≤2`, degree two | 76,426 | `[1986144,2766554]` | `[1983767,2762071]` | `[1198,5641]` |
| Disjoint subset, free rates | 5,448 | `[81783,217495]` | `[81740,217443]` | `[-318,462]` |
| Disjoint subset, degree two | 76,426 | `[96942,189996]` | `[96893,189936]` | `[-135,278]` |
| `D≤4,k≤2`, free rates | 618 | `[854893,1439955]` | `[852382,1435280]` | `[1959,5330]` |
| `D≤4,k≤2`, degree two | 6,568 | `[880412,1351299]` | `[877831,1346719]` | `[2174,5025]` |
| `D≤4,k≤2`, through degree four | 280,023 | `[891295,1320661]` | `[888679,1316121]` | `[2252,4915]` |

Four finite trees of size 16,000 (seeds `122419..122422`, 570 cores, 50 matching full-mask paths) gave actual reward/core `2.6330620e−10 / 2.6331555e−10` and contrast `−9.350e−15`. The disjoint mask had two matching paths and sample reward/core `1.0673651e−11 / 1.0569171e−11`, contrast `+1.045e−13`. These finite snapshots check exact payoff reconstruction but do not estimate the all-age sign or certify it. In particular the full-mask finite-sample contrast has the opposite sign from the stationary result.

## Physical lifetime and branching registration

From each Palm newborn root, start one length-one prefix. For each leaf promotion of the last tracked vertex of a length-one or length-two prefix, add the continuations of **both** branches: the original prefix, with the newborn registered as its side child, and the extended prefix, with the same newborn registered as its selected child. A promotion away from the terminal role just registers a new side child. The existing registrations, which may share a physical vertex, do not split its clock or create an independent copy of the physical tree. Linearity of the conditional expected sum justifies adding their continuation bounds. Every ordered `hV` has a unique first chosen edge and second chosen edge, hence exactly one completed registration. Only completed registrations pay, and they pay at every eligible time. Cutoffs on degree and arity are irreversible. A truncated prefix with terminal arity equal to the cap cannot generate another live descendant registration; it need not be tracked further.

All path-adjacent receives are exact. On each side the sorted registry is a multiset, with transitions multiplied by multiplicity. The no-memory model uses a free rate in `[0,1/2]` per side child from birth. At degree cap 2, `D=2,k=0` is virgin (gain `1+1/D_path`, promotion `1/2`); `D=2,k=1` is bare (gain in `[1/D_path,1/D_path+1/2]`); `D≥3` is saturated and sends in `[0,1/3]`. At larger caps retain the full known `(D,k)` until the next gain crosses the cap. A known side child sends exactly `1/D`, receives from its path parent at `1/D_path`, and has only its own child rates unresolved. Saturation never restarts the clock. The optional mother cap retains exact reciprocal through that degree and allows every Palm initial class above/below cap. With no mother arity observed, the allowed gain range is `[p,p+D_mother−1]` with `p=1/D_root`, which covers even a mother without a core parent.

## Payoffs and one physical discount

Both blob hashes are checked before extracting any rows. `Local.charge` supplies the signed, transferred all-flow row. Its exact zero-pressure intercept and directional slopes are reconstructed in the replay, checked against the analytic floor, and checked on every matched finite sample against the true charge at the actual adjacent degrees. The transition carrying a rate carries its signed affine credit as well. In the contrast, subtract lower/upper witness row coefficients *before* constructing the Bellman system, on the same path, registry and physical gain rate. Separate witness intervals alone would not certify its sign.

Discount and stable-age prefactor are the same physical `μ∈[LO,HI]` in all three objectives. The replay encloses every fixed `μ` in that interval (a wider relaxation permits the discount to vary by state); the final multiplication by `μ` stays inside the same interval and takes the signed extrema. `--mu-cells` can split it into narrower shared intervals. Each Bellman endpoint optimizes a ratio affine in the independent rate boxes, maintaining the common holding-rate denominator. Integer rates/credits are exact; only the final continuation ticks are outward rounded. A signed corner-enumeration regression checks both slopes, the holding denominator and the discount endpoint. Discounted Dynkin or conditional first-step induction bounds every adapted physical completion, without a Markov assumption on the omitted environment. Apply the stable-age identity to the birth-rooted sum to obtain the stationary results above.

## Disjoint support and booking

For an already registered path let `n=Σ k_i−2` be the number of immediate side children. The apex fringe has at least `3+n` cores and at least `Σ(D_i−k_i)+2n` ordinary vertices: the latter counts the three path cores with their *current* leaf stock plus the independent two-vertex minimum of each side subtree from its birth. No subtree shrinks. Thus `--disjoint` pays only if `3+n≥5` and `Σ(D_i−k_i)+2n>12`. This implements a reliable monotone lower counter without assuming an untracked tail subtree has any particular size. It is disjoint from the [size-at-most-12 release](../../release/README.md) and the size-above-12 suffix accounts that require at most four cores. Its apex degree cap also excludes the older, old-flow and high-apex `hV` masks. The reward is the live `hV` charge and does not spend the `g` promotion that originally created it. The same `hV` objective is not additive with other calculations on the identical mask.

The two positive disjoint lower bounds (`9.6942e−12` lower and `9.6893e−12` upper) can be added at the numerator level to the stated release. Equivalently, apply those numerators with the independently validated conservative [all-degree ceilings](../../denom/README.md) `1.356833021311 / 1.367169995210` to the release's exact fractions, then round outward. This gives the optional combined interval **`[0.565675793091, 0.565772420371]`**. The disjoint contrast remains unresolved. The full-mask calculation overlaps that release and must not be added to it.

## Replay

In the Astra checkout, set `PYTHON` to the existing NumPy environment, `SCRIPT` to this `replay.py`, and `LO_BLOB`/`HI_BLOB` to the original pinned production blobs. The result table used:

```sh
taskset -c 3 nice -n 17 timeout 140s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  "$PYTHON" -Bu "$SCRIPT" "$LO_BLOB" "$HI_BLOB" --donors 0 2
taskset -c 3 nice -n 17 timeout 80s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  "$PYTHON" -Bu "$SCRIPT" "$LO_BLOB" "$HI_BLOB" --donors 0 2 --disjoint
taskset -c 3 nice -n 17 timeout 55s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  "$PYTHON" -Bu "$SCRIPT" "$LO_BLOB" "$HI_BLOB" --degree 4 --donors 0 2 4 --sample 0
```

Respectively: `wall=31.11 rss=131424 exit=0`, `wall=62.55 rss=119868 exit=0`, `wall=50.51 rss=146152 exit=0`. Degree-four side tracking with a wider path cap or mother registry is substantially more expensive than the default; no `D≤7` degree-four side result is needed for the comparisons above.
