> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Older Neon D18 charges: exact census and conditional occupation

The matched replay proves a **small positive stationary price** for old, nonvirgin `hV` paths on moderately crowded hosts: `2.6657003958914333e−9` lower and `2.3601199127282544e−9` upper in the original signed reward units. With the [all-degree denominators](../../denom/README.md) these give inward moves of approximately `1.965e−9` and `1.726e−9`. This is about 0.02–0.03% of a `10⁻⁵` reward benchmark per sign; the ambitious target is not attained. The most promising held-out `P4` and `Yt` descriptors have a quantified missing-send obstruction below rather than a certified occupation.

## Current-time masks and composition

The matched `--older-path` result pays only `hV`, with a nonfounder apex, no virgin in the occurrence, at least five cores and total size above 12 in the apex's **complete** descendant fringe. Its apex mask is `(13≤D≤17, 3≤k≤6)` or `(18≤D≤25, 3≤k≤4)` at replay settings; middle and tip degrees are at most 5 and 4. The upper apex cap can be increased through 64 under the same split. Marking a birth narrows this reward further and never enlarges its support.

| Other calculation | Relationship to the matched price |
|---|---|
| First-body `f`, mixed `g`, and `high-parent` virgin-pair `hC` | Different charge owner family; no overlap. |
| `higher/` size ≤12; `arms/` and `renewal/next/` higher reward on size >12, cores ≤4 | Disjoint by complete-fringe size/core count. |
| 2864 old-flow main `D=18…64, k≥5` apex | Disjoint by the explicit degree/arity exclusion. |
| 2417 `D≥65, k≥5` apex | Disjoint by apex degree cap. |
| Any shared-age or lower-degree certificate for the same `hV` occurrences | Compare by max on overlap by default; add only after proving an exact disjoint support split. |
| `--tagged` preliminary replay in this directory | Overlaps `--older-path` and old-flow, can include virgin tips, and is not in the additive result. |

No combined endpoint is booked here. The numerator and denominator always refer to the same original pinned sign, never to an interpolated witness.

## Exact scanner and sample-only ranking

`Local` reuses the mmap-backed `Cert` and exact table evaluator, keeping only bounded caches of sparse evaluations. It enumerates supported occurrences once. For each occurrence and each incident receiver, it computes the gain and promotion increments once, using actual external and internal reciprocal sends. A new motif is paid in the old connected embedding obtained by deleting its newborn; the newborn `f` goes to the promoting `f`. Separately reconstruct each nonnegative charged row by applying its centered outgoing transfers and incoming tip-deletion corrections. For tail `f` the transferred derivative is `c2/D²+c3/D³`; for tail edge `g` it is the exact T2/T3 combination, using the same Bernstein embedding of the other slot. Above `k=17` these transfers are disabled as in the pinned proof. The physical incident sum and the charge sum agree on every sampled tree. On 28 fixture trees per blob they also agree with a separate calculation that recomputes `Cert.potential` after each complete-tree gain and promotion. The fixtures cover all eight motifs, mixed tails, and the 17/18 arity and degree boundaries. See [the higher-order provenance](../README.md#production-obligations-and-provenance) for the pinned source and full GMP obligations.

The rank uses real finite friend-tree histories: choose an ordinary vertex uniformly, then choose its neighbor uniformly and attach to that recipient. `size > 12` is applied to the complete apex fringe and the founder is excluded from paid rows. The descriptor is `(core-count bin, family, max-degree slot, degree-band start, host k bin, number of virgin motif slots)`; the first maximal slot breaks ties, `5` means `≥5`, degree bands start at 8, 18, 32, 64 with `0` below 8, and `8` for arity means `≥8`. The denominator is the total number of sampled cores. Four matched size-16,000 trees at seed `122419` give 546 cores. All numbers below are **finite-sample ranks**, not stationary estimates usable in a proof.

| `cores≥5, size>12` sector | Lower sample reward/core | Upper sample reward/core |
|---|---:|---:|
| All slots inner | `1.298e−5` | `1.023e−5` |
| Apex tail | `1.063e−5` | `1.558e−5` |
| Apex inner, another slot tail | `2.083e−6` | `2.242e−6` |
| Total, all 6 families | `2.569e−5` | `2.805e−5` |

The sample includes large nonvirgin `hV` rows, as does 34fa's held-out census. The biggest sample `P4/Yt` rows are often driven by the children **of a middle or terminal motif slot**, not by the other children of the old apex. For any fixed slot descriptor, all unresolved receive sums enter affinely, and each is bounded by half its unobserved count. The exact rational extrema give these examples (state triples below are `(L,k,h)`, in motif order):

| Side and fixed descriptor | Sample row | Uniform row floor / ceiling | Dominant contribution to admissible width |
|---|---:|---:|---|
| Upper `Yt`: `(4811,8,1),(0,2,1),(0,1,1),(1,1,1)` | `.001902` | `6.006e−10 / .002097` | First tip's child sends `.001902`; second tip `.000195`. |
| Upper `P4`: `(78,2,1),(0,2,1),(0,1,1),(7,0,1)` | `.001600` | `1.176e−9 / .001600` | First middle's child sends `.001597`. |
| Lower `P4`: `(1260,9,1),(0,1,1),(0,2,1),(1,1,1)` | `.000708` | `6.871e−10 / .001027` | Terminal child sends `.000708`; next middle `.000320`. |
| Lower `hV`: `(34,5,1),(0,1,1),(10,2,1)` | `.000511` | `1.930e−6 / .000737` | Terminal child sends `.000466`; host child sends `.000216`. |

Using these fixed descriptors alone to move an endpoint by `10⁻⁵` at the supplied ceiling, even pricing every occurrence at its ceiling would require incidences/core of at least `.00652`, `.00854`, `.01321`, and `.01841` respectively. Using only the uniform floors would require about `22763`, `11628`, `19747`, and `7.03` incidences/core to suffice. These are required-mass/width diagnostics for those exact descriptors only, not a bound on their stationary prevalence or on the entire omitted sector. Bounding the old host's omitted child sum below by `count/4` leaves the upper `Yt` and `P4` floors unchanged. A useful next projection must retain low-degree senders at the **specific receiving slot** while allowing the host and unrelated branches to grow.

## Conditional physical price

Upon every host promotion, start one marked middle-core valuation `A`, and upon **every** promotion by that middle core start a tip valuation `V`. These clocks are simultaneously live, with no prescribed chronology and no independence. At fixed physical time an oriented `hV` occurrence identifies a unique middle and tip birth, so there is no multiplicity from multiple starts. A host's other children and a middle's other children remain live. The `V` state is `(D,k)` for the host, `(d,j)` for the middle, and `(e,m,v)` for the tip (degree, core-child count, and virgin core-child count). The tip's `v` virgin children each send `1/2`. Their losses have exact aggregate rate `v(3/2+1/e)` and include child's own gain/promotion and the tip's send to it. Loss of virgin status is irreversible. The tip's promotion increments both `m,v` and keeps `e` fixed. Host and middle promotions have their true rates and retain `V`; host/middle gains include the exact known reciprocal arrows and bound remaining receives independently. These are pointwise adapted rate intervals, not conditional independence assumptions. Unknown sends used in the reward can be minimized separately from the transition hazards, which is conservative and is one source of loss.

For each finite projected state, evaluate the discounted recurrence at all vertices of the receive-rate rectangle. The numerator contains the reward floor plus discounted continuation rates, the denominator is `HI` plus the total projected exit rate, and every rounded value satisfies all cross-multiplied inequalities. Coordinatewise monotonicity of degrees and child counts, with virgin-loss transitions decreasing `v` at fixed other coordinates, makes the recursion acyclic. In the middle-only state `A`, a promotion pays the immediate continuation value `V(post-promotion, virgin tip)` **and** continues `A`. A host-only census prices new middle starts by the same compensation. Apply discounted Dynkin at each stage under the full physical filtration and Tonelli over nonnegative births. All bound values are finite and bounded on the cuts. The host starts as a physical newborn under any complete birth environment. Multiply by `LO` under the stable-age identity, with `HI≥μ` in all stages. This prices the actual sum of path incidences, with marks used only as a subreward; it does not assign a stationary probability to a simulated descriptor.

## Replay

Use Python with NumPy, the original payload paths, and the linked Astra helper layout. For each `BLOB` (both sides), run on one assigned CPU:

```sh
taskset -c "$CPU" nice -n 17 timeout 35s "$PYTHON" -Bu "$SCRIPT" "$BLOB" --n 16000 --reps 4
taskset -c "$CPU" nice -n 17 timeout 35s "$PYTHON" -Bu "$SCRIPT" "$BLOB" \
  --reps 0 --older-path --host-cap 25 --kcap 6 --child-cap 5 --tip-cap 4 --verify-hash
```

`SCRIPT` is this directory's `replay.py`. The matched hashed occupation runs exited 0 (17.98 s lower, 15.40 s upper, under 115,000 KiB RSS) and evaluated `51,774` tip states and `2,787` middle-only states per side. Exact outputs are `13328501979457167255346839/5000000000000000000000000000000000` lower and `11800599563641271253668973/5000000000000000000000000000000000` upper. A larger `(32,6,6,5)` exploratory replay finished only on the upper side (`5.88182277606657e−9`, 33.95 s; lower incomplete under the same limit). It is not an extra installment on the matched support.
