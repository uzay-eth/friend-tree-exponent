> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Additive `hV` donor installment on top of `ead947b`

Certified stationary numerator floors are **`10361983749/10^16` lower** and **`14352363253/10^16` upper** per core. Both exceed `10⁻⁶`. With the unchanged denominators and exact endpoints from the [post-EDGE freeze](../../../interpolate/compose/residual/edge-projection/combined/FREEZE.md), the standalone outward enclosure is **`[0.565678819949264, 0.565769419622021]`**. This is before adding the separate `P4/Yb/S3` package; use the original exact endpoints or the exact combined numerators when merging packages.

## Exact reserve exclusion and charge ownership

Let `S,C` be the size in ordinary vertices and number of cores of the **complete planted fringe at the apex** (exclude its parent). At an `hV` occurrence `C≥3` automatically. The original frozen higher reserve uses size ≤12, and the size-20 higher suffix permits 3 or 4 cores. The exact complement of this size-and-core predicate is

```
S > 20  or  (S > 12 and C >= 5).
```

The old geometric interior restrictions only reduce the reserved support; this replay deliberately discards their extra complement. The post-EDGE residual applies only to `f` and existing `g`; it does not consume the original `hV`. The frozen mixed `g` and virgin `hC`, and the upcoming `P4/Yb/S3` installment, are different charge owners. In particular, the birth of the `hV` weight already credited to a promoting `g` is not paid here. None of the tentative lifetime/old-flow `hV` ledgers is an input; any future proposal to combine those requires its own occurrence/port comparison.

Only proven instances of the exact predicate are paid. For spine slots `(d_i,k_i)`, set `n=Σk_i−2`. A lower count of the actual fringe is `Ĉ=3+n`, `Ŝ=Σ(d_i−k_i)+2n`. All immediate side-child subtrees are disjoint and have size at least two. For a tracked side donor with exact `(d,k)`, replace its two-vertex minimum by `d+k`, and add `k` to the core lower count. In saturation retain the unrefined bounds. Evaluate the displayed predicate on `(Ŝ,Ĉ)`; it is increasing in both arguments, so a positive test certifies the exact physical complement. In particular, every paid occurrence at apex degree ≥21 automatically has actual focal size >20. Ambiguous smaller fringes are unpriced, not classified as outside by a proxy. Upward-parent atoms are exterior to the focal fringe and do not change these counts.

The [atom recurrence](README.md#physical-side-atoms-and-all-time-results) splits each live affine row into its nonnegative all-rate baseline and nonnegative per-arrow excess `c p − min(0,c/2)` for actual reciprocal `p`. Only the listed arrow excesses are paid. Side-donor registrations enumerate all legal orders of the spine births and the individual donor birth. Different roles are distinct arrows; the two apex bands are disjoint in **current** apex degree. The cap on degree and arity is a calculational cutoff per tagged role, not a cut on untracked siblings.

For the apex's **parent arrow** only, pay on Palm births whose mother has degree at least 3 immediately before promotion. Her degree stays nondecreasing, so this arrow has `p≤1/3` throughout the lifespan, including every prefix and its future branches. A degree-two mother can give birth only from virgin stock, at rate `1/2`. Virgin status is irreversible and its exit hazard is at least `3/2`, thus `π(virgin)≤μ/(μ+3/2)` by the same stable-age identity. Consequently the Palm weight of births with mother degree at least 3 is at least `1−(1/2)/(LO+3/2)`. Founder intensity vanishes. The replay applies this factor to the nonnegative conditional lower floors and rounds down again; it does not use the conditional contrast as an unconditional one.

## Exact ledger and captured replays

Each entry is a lower bound in `10⁻¹⁶` units/core, with both original blob hashes verified. All runs used `--frozen`, default `--old 8` (`max D≥8`), and:

```sh
taskset -c 3 nice -n 17 timeout 145s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  ./round2/.venv/bin/python -Bu \
  round2/occupation-neon/pressure-memory/hv-lifetime/old/replay.py \
  b3f18_lo.bin b3f18_hi.bin --frozen <flags>
```

| Paid role | Flags | Lower ticks | Upper ticks | States | Wall / RSS KiB / exit |
|---|---|---:|---:|---:|---|
| Apex child | `--caps 20 5 5 --arity 6 2 2 --atoms apex --atom-cap 4` | 5700604254 | 10578591182 | 170010 | `60.87 / 271700 / 0` |
| Apex child, 21–48 | `--caps 48 5 5 --arity 6 2 2 --atoms apex --atom-cap 4 --host-min 21` | 1816107172 | 2889332514 | 463058 | `134.81 / 307348 / 0` |
| Parent, apex ≤20 | `--caps 20 5 5 --arity 6 2 2 --atoms parent` | 1629245400 | 449615529 | 18470 | `48.25 / 249788 / 0` |
| Parent, apex 21–48 | `--caps 48 5 5 --arity 6 2 2 --atoms parent --host-min 21` | 590196762 | 126911006 | 49438 | `78.67 / 179780 / 0` |
| Middle child | `--caps 11 10 10 --arity 3 2 2 --atoms middle tip` | 174005659 | 98372532 | 76382 | `49.81 / 261624 / 0` (joint invocation) |
| Tip child | Same invocation | 451824502 | 209540490 | 85103 | same |
| **Total** | | **10361983749** | **14352363253** | | |

The exact inward gains in the exponent are approximately `7.633171973762067e−7 / 1.0497785976805525e−6`, using the unchanged denominators `1357493815758092773/10^18` and `273436004214812983/(2·10^17)`. Exact standalone endpoints before outward rounding:

- Lower: `137071131711883506705531412118897/242312646112819559980500000000000`.
- Upper: `6574823499856115541816419305857/11621030179129551777500000000000`.

The larger complete-fringe complement is not exhausted. In particular, the price leaves all ambiguous size/core states, larger role caps, the other Palm parent class, and the separate all-rate `hV` baseline unspent. Only the displayed installment is certified here.
