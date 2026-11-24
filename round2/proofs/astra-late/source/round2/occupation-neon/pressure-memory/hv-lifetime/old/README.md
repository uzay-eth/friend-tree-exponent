> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Older `hV`: individual exterior atoms

**One actual donor at a time works.** On the matched path cut `D≤(11,10,10)`, `k≤(3,2,2)`, `max D≥8`, the exterior-child excess has a zero lower bound for both witnesses with free-input controls. Three separately tracked atom families raise the simultaneous lower bounds to **`2.606190617e−7` lower and `4.159567571e−7` upper per core**. No snapshots or assumed entry frequencies are used. The signed lower-minus-upper interval still contains zero. This matched cut does not establish `10⁻⁶` for both witnesses.

## Exact routing and signs, including the unbounded tail

At an extant nonfounder path `(A,B,C)`, put `Δ_i=H(x+gain_i)−H(x)` for the original `hV` weight. The normalized signed external slope is zero if routed, otherwise `sign·Δ_i/S`. Absent ports have zero multiplicity.

| Exterior port | Multiplicity | Routed exactly when |
|---|---:|---|
| Parent of apex | 1 | `A,B` both in `P4` domain |
| Apex side children | `k_A−1` | `A` in `Yb` apex domain and `B` in `Yb` middle domain |
| Middle side children | `k_B−1` | `B` in `Yt` fork domain |
| Tip children | `k_C` | `B,C` both in `P4` domain |

The two apex signs agree whenever both are free. There is no `hV` row above middle arity 17. Tail gains use the exact *difference* of those `hV` weights, not the special asymptotic `f/g` transferred derivative. `--sign-only D_A k_A D_B k_B D_C k_C` classifies all four ports on any individual legal census, with no degree limit. For a tail slot the integer homogeneous Bernstein weights are `((D−18)^3,54(D−18)^2,972(D−18),5832)` in denominator `D³`. Replace the gained slot's weights by the cross-multiplied signed difference at `D+1` and `D`, retaining the union of inner/tail keys at the seam. Contract with the other two slots' homogeneous weights and take the integer sign. The replay checks against the original rational slopes. For example at `(40,5),(2,1),(13,2)` both apex ports and the tip port are negative in the lower witness and positive in the upper; the middle port is absent. This gives a test for every individual census, not a constant sign assigned to an entire degree band.

## A certified small cut for substantial free-rate demand

The **box demand** is `W=Σ n|c|/2`, the span of the independent exterior-rate box. It is a nonnegative census-only reward; its stationary expectation quantifies how much uncertainty the free box admits, even though it is not itself attainable excess. On the nested cut `D≤(d,7,7), k≤(3,2,2), max D≥8`:

| Apex ceiling `d` | Stationary `W`, lower witness | Stationary `W`, upper witness |
|---:|---:|---:|
| 8 or 9 | 0 | 0 |
| 10 | `[0.7117434873,1.9765766608]e−6` | `[0.8080353114,2.2463888044]e−6` |
| 11 | `[1.2600518260,3.3977273431]e−6` | `[1.4218658263,3.8401285238]e−6` |
| 12 | `[1.6892372166,4.4506066159]e−6` | `[1.8971593733,5.0065165629]e−6` |

Thus `d=10` is the first nonzero shell, and `d=11` is the first **certified** above `1e−6` for both in this ladder; the certificate does not exclude the threshold at 10. On the stricter size/core-disjoint subset at 11, apex-side demand alone is at least `4.013804584e−7 / 4.588519809e−7`; upward-parent demand is at least `2.791301021e−7 / 3.115173149e−7`; both other ports route or are absent there. With the wider cap `(11,10,10)`, old middle and tip ports have real demand and positive prices, as below. No global minimality over other arity cuts is claimed.

## Physical side atoms and all-time results

For each exterior donor with reciprocal `p`, pay `c p−min(0,c/2)`. This equals `c⁺p+(−c)⁺(1/2−p)`. The independent GMP baseline is `a+Σ n min(0,c/2)≥0`. It is counted zero times in the atom solve. Affine reconstruction, including a mixed-tail row, is checked exactly. The atom's true reciprocal enters the carrier gain *and* its reward through the same edge variable. The other donors have free gain inputs but no price in this solve.

Start at Palm apex birth. A promotion can continue the current registration, extend the main path, or introduce the distinguished side atom. Its actual rate multiplies the sum of these continuations. One birth never simultaneously fills both labels. A side atom before the selected spine child survives in the extended prefix; a spine child chosen earlier can also acquire a later side atom. Later distinct side births are spawned by the unmarked prefix, not by already marked copies. Thus every legal `(path, exterior child)` is counted once, including atoms that precede the completed path. No incumbent sibling needs a registry. Gains at the atom are live, with its degree-2 virgin/bare sends exact and saturated sends bounded by `1/3`; there is no reset. The registry cap can be raised independently. The lower, upper, and contrast objectives use the same physical age law and control each affine incoming rate jointly with the holding denominator. The exact directed integer Bellman argument is the same as [the young lifetime proof](../README.md).

Matched cut `(11,10,10)` / `(3,2,2)`, all figures in `10⁻¹⁶` units/core:

| Price | Lower witness | Upper witness |
|---|---:|---:|
| Entire charged row with free boxes | `[946766,45381134038]` | `[947274,45719099486]` |
| All directional excess with free boxes | `[0,45379802437]` | `[0,45717772571]` |
| Apex child atoms | `[1637662009,12509520050]` | `[3705515171,17105955411]` |
| Middle child atoms | `[300548104,1649054344]` | `[173974193,1342428230]` |
| Tip child atoms | `[667980504,3659139521]` | `[280078207,2158282042]` |
| **Sum of three child objectives** | **`[2606190617,17817713915]`** | **`[4159567571,20606665683]`** |

The contrast bounds from summing three joint signed enclosures are `[-14376157315,10271201087]×10⁻¹⁶`. On wider `(18,10,10)` / `(3,2,2)`, the middle and tip alone give lower bounds `1.023420818e−7 / 6.58754649e−8`. A different wider apex-only cut `(48,5,5)` / `(6,2,2)` gave `5.484708970e−7 / 1.0919936359e−6` (upper crosses the target; that run emitted the final row but lacked a captured exit receipt). These masks overlap the smaller cut and should not be added to it. On the size/core-disjoint criterion inherited from the young solve, cut `(16,5,5)` / `(5,2,2)` gives positive apex-atom floors `2.808185782e−7 / 5.423674831e−7`. That size/core test guarantees disjointness from the size-12 release and the three/four-core suffix, but the older `(13…17,k≥3)` releases may overlap. No headline endpoint is booked. Do not add separately tracked atoms to a whole-row Bellman price on the same support.

## Replay evidence

Use `./round2/.venv/bin/python -Bu` with `round2/occupation-neon/pressure-memory/hv-lifetime/old/replay.py b3f18_lo.bin b3f18_hi.bin`. Both blob hashes are checked each invocation. The following flags were run under `taskset -c 3 nice -n 17` and bounded timeouts:

| Flags | Wall / exit |
|---|---:|
| `--caps 12 7 7 --ladder 8 9 10 11 12 --track 0 --objective width` | `20.97 s / 0` |
| `--caps 11 7 7 --track 0 --objective apex-demand tip-demand middle-demand parent-demand --disjoint` | `33.83 s / 0` |
| `--caps 11 10 10 --track 0 --objective reward excess` | `47.24 s / 0` |
| `--caps 11 10 10 --atoms apex middle tip` | `53.14 s / 0` |
| `--caps 18 10 10 --atoms middle tip` | `80.16 s / 0` |
| `--caps 16 5 5 --arity 5 2 2 --atoms apex --disjoint` | `19.67 s / 0` |

The next target is a larger shared cut with a sharper lower-witness apex loss bound (longer atom memory or parent-stem loss), adjacent to the measured middle/tip child contribution. A demand certificate is not a substitute for proving that loss or send is actually present.
