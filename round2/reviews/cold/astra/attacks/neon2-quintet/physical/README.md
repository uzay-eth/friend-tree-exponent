# Neon-2 Claim 6: sparse exact physical-tree audit

**PASS on both gold1 endpoint payloads.** This is an independent labelled-graph test of the physical potential, one-birth changes, signed global drift, and comparability. It does not certify the exhaustive local corner inequalities.

The only inputs were the read-only clean detached `round2/.work/neon2-claim6-7M6Ys3rJ/source` at `b04cc50c1e391e0420eb7b006ac781d484c68543`, and the `quintet_lower.bin` and `quintet_upper.bin` copies in its parent. The source blobs used for the wire format, proof and claimed caps are `647a892cfbcbea11efd8c2bc5dc9898f8d564e06`, `708f4a6b95f02a47e74e2a89b192640924ebd1c5` and `59b00308b778ee0142087d2477d9fa7312df3ca2`, respectively. The script checks the source pin, binary mode, nominal rational, file size, layout length, word bounds on reads, and the exact `g(INF,INF)=0` entry. No full-payload digest or table scan is done here.

## Independent physical check

`check.py` implements a sparse little-endian reader (`pread`); it imports no claimant code. It reconstructs wire ranks for each signed-int64 potential table, applies endpoint interpolation in `Fraction` per labelled endpoint, and uses the analytic singleton formulas at `D>=DF`. It skips centering, FY and envelope payload sections instead of using any transferred row to compute physical drift. It checks 1,125 distinct payload words per endpoint.

The first path expands connected vertex **label sets** of size at most five and classifies the induced shapes. A vertex through which an already connected set of size at least two expands becomes internal, so it is safe to prune expansions through vertices outside `M3`; all allowed internal supports are nested in `M3` in these payloads. The actual ordinary adjacency lists enumerate every directed donor slot at rate `1/D_source`. The checker groups slots only when they lead to the *same recipient*, appends a real new leaf, reconstructs degrees and induced cores from scratch, and re-evaluates each labelled motif. It checks that total slot rate is ordinary tree size.

The alternate path independently generates the eight shapes from centers, ordered routes modulo reversal, or combinations of distinct labelled arms. It applies `D↦D+1` at a core gain, or `c↦c+1` at the promoter plus the actual old leaf label as a `(2,1)` core, and reconstructs all motif scores for that update. It checks the exact raw increment for each labelled motif, including old contributions that drop to zero and all newborn contributions, before summing by family. It checks that every newly supported nonsingleton includes the fresh core and has exactly one supported parent on deletion of its label. The exact rate identity is `D−c+Σ_{w core neighbor}1/D_w` at a core recipient; for each individually labelled promotion its rate is `1/D_promoter`, yielding `(D−c)/D` upon summing. The two enumeration paths share the wire lookup, so their equality checks the geometry, weights attached to those embeddings, and physical changes, not independent binary decoders.

The 22 trees include `C=1,2`; paths with equal and unequal internal states, mixed tail ends and a scrambled labelling; a degree-six hub with five isomorphic arms; three isomorphic length-two arms producing brooms and 5-paths; support drops at the H3, S3, P4, P5, B5 and S4 boundaries; exact-degree overflow and gains across the `d2` knot; adjacent physical tails at degrees 199 and 237; and stars and chains crossing `DF−1→DF`, including a chain whose internal vertex crosses the seam. The tests tally both creation and a nonzero raw increment for every family. The equal-arm case counts the five size-four arm subsets separately, even when every lookup tuple is the same.

## Signed hypotheses and comparability

In ticks `W=S·Φ`, the nominal hypothesis is `σ(ℒW−λW)≤E C`, with `σ=-1` for the lower payload and `+1` for the upper. The inspected aggregation uses one singleton per core, `C−1` physical edges, and caps charged to the finite internal hosts for higher motifs (one-half of P4 at each internal). The check uses both the edge-count-sensitive bound and the reported `E=max(0, cap_f+max(0,cap_g)+higher_per_core)`. Overflow/tail interpolation uses true degree, and two joined tails keep both directed reciprocal donor rates. Analytic singleton gains use the true successor states, including the seam, and an analytic vertex may still be an **endpoint** of a higher motif hosted elsewhere. Beyond `DF` the analytic *singleton row* has no higher newborn allocation with that vertex as an internal host; this does not justify deleting those existing higher motifs from the potential.

For the global lower comparison, the lower payload's nonnegative shift is `α=1/4`; the upper has `α=0`. For `C≥2`, charge each edge to its rooted child (zero edge for the root), each higher motif to the specified host, and add `α S(c−2)` to each local lower-accounting row. The identity `Σ(c−2)=-2` gives back **positive** `2α S` in the global lower bound. For `C=1` the separate singleton bound is unshifted. In the analytic region the lower comparison treats `c=0,1,2,≥3` separately, including a possible allocated far-tail edge at `c≥1`. The upper accounting is unshifted. The checker verifies these host inequalities for two root choices on the fixtures, and checks `0<m_min C≤W≤m_max C` on every fixture and every physical successor. This directly tests the claimed global bounds on the sample; the universal step also requires the reported per-host extrema and signed local residual hypotheses.

Since `W≥m_min C`, the signed cap implies `σ(ℒW−(λ+σ E/m_min)W)≤0`; the script recovers that endpoint rational independently and checks its physical sign. The upper and lower endpoint slacks are each strictly positive on all fixtures. No product with a rounded decimal endpoint is used.

| Endpoint | Nominal cap / core | Worst exact signed nominal residual / core | Least exact adjusted sign slack / core | Sample potential / core range (ticks) |
|---|---:|---:|---:|---:|
| Lower | `19186` | `−10072217322379/4900000000` (arity seam) | `49898323235430010681698187/244711634769700000000` | `[6847915552977/7, 1808362214356]` |
| Upper | `19783` | `−228192582809/50000000` (two cores) | `254861973760945713652183/1552367714350000000` | `[6848855011458/7, 1808494615915]` |

Those slacks compare against the exact adjusted lower `2825181729593748861/4994114995300000000` and upper `439124781075521909/776183857175000000`. In particular the sampled worst residuals are already negative at the nominal targets; this says nothing about the worst unsampled tree.

## Exact replay and negative controls

From `[local path]`:

```text
python3 -IB round2/attacks/neon2-quintet/physical/check.py
lower: 22 trees, 20994 arrows, 10519 recipients, 1125 distinct payload words
upper: 22 trees, 20994 arrows, 10519 recipients, 1125 distinct payload words
```

Both modes gave these counts (a nonzero comparison is per recipient and labelled motif, not a count of distinct payload addresses):

| | f | g | H3 | P4 | S3 | P5 | B5 | S4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Initially supported | 196 | 174 | 64 | 22 | 34 | 11 | 5 | 20 |
| Nonzero raw comparisons | 20842 | 24870 | 4219 | 2835 | 461 | 501 | 303 | 217 |
| New supported | 10323 | 10323 | 498 | 226 | 132 | 55 | 191 | 17 |
| Dropped support | 0 | 0 | 1 | 2 | 1 | 3 | 2 | 15 |

On a hub leaf gain, collapsing S4 arm choices by their equal state tuples erases the labelled contribution `(0,1,2,3,5)`: the exact increment is `−1807707` lower and `−6859561` upper, versus mutant `0`. On a gain at the hub with three long arms, the broom `(0,1,2,3,4)` increments by `+10079660` lower and `+3271430` upper, versus the respective negatives in the sign mutant. The script requires those disagreements so a silent mutation fails the run. No mutation touches either input binary.
