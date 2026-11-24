# Independent audit: complementary P4 by distinguished donor

**Correct for the matched release.** With tag degree caps `(6,5,4,6)`, arity caps `(2,3,2,1)`, side/mother memory `3/3`, and the entire previous `(12,8,6,16)` band excluded, the fresh replay reproduces these stationary lower bounds on the two original signed charges per core:

| Witness | Exact floor | Approximate |
|---|---:|---:|
| Lower | `8878929713164459/125000000000000000000000` | `7.1031437705e−8` |
| Upper | `54446048547132641/10^24` | `5.4446048547e−8` |

The inherited inputs are the two original D18 proofs and the physical birth-Palm/stable-age normalization. This reviews the matched run and the method; it does not certify an unpublished wider run.

## One minimum plus one term per actual donor

For each labelled physical path and one signed witness, write the raw charge `a+Σ_y w_y p_y`, where `p_y=1/D_y`. The sum includes exactly one mother incidence at `A` and one incidence for **every** core child outside the path at its attachment vertex. The slope is common within a slot, but the reciprocals need not be equal. The exact decomposition is

```text
m    = a + Σ_y min(w_y,0)/2
h_y  = w_y p_y − min(w_y,0)/2
raw  = m + Σ_y h_y.
```

Here `m≥0` on every priced row (checked by the producer, on both payloads) and `h_y≥0` for `p_y∈[0,1/2]`. The minimum is booked once per path; a donor component books only its own increment. Neither overlapping paths that use the same physical `y` nor two choices of `y` for the same path imply independence. By linearity the separate floors add, and the same donor on two different paths is counted for two different charged embeddings, as required. The formula does not book all of `m` again on each incidence, price just one representative per slot, or multiply a single reciprocal by the donor count. An `h_y` can vanish; zero-clamping a lower bound uses nonnegativity, not strict positivity.

The exact port `p_y` or its overflow interval appears **once** in the Bellman comparison: its receive edge contributes `p_y[V(gain)−V+w_y]`; `−min(w_y,0)/2` belongs to the running intercept. This is an algebraic decomposition of the instantaneous toll, not a second compensated payment at a jump. The mother component uses the same rule, with the mother existing at the apex's actual birth. Slope rounding is outward because `p_y≥0`. Sending from `y` to its attachment increases the recipient's degree while `y` keeps its own degree; arrivals at `y` and its own promotions are separate physical arrows. Grouping identical unseen-rate axes gives a safe box; correlations between different rates, objects, witnesses and discount are relaxed, while the single rate multiplying a donor's combined coefficient is preserved. Independent corner enumeration checked `extrema` on signed rows, including an adversarial cancellation where splitting the shared rate manufactures a negative lower bound.

## Birth order, arity, and exit

Fix a potential physical final pair `(A→B→C→D, y)`, with `y` a child of one of the four path vertices and not a path tag. The extended registrations begin at the promotion of `A`. The path tags must be born in order `ABCD`; `y` must be born after its attachment. There are **10** orders across the four slots:

| Attachment | Legal promotion orders |
|---|---|
| `A` | `AYBCD`, `ABYCD`, `ABCYD`, `ABCDY` |
| `B` | `ABYCD`, `ABCYD`, `ABCDY` |
| `C` | `ABCYD`, `ABCDY` |
| `D` | `ABCDY` |

At a tip promotion the new child can be a path extension **or** the newly distinguished `y` for the old prefix. These are separate summands at one `L/D` rate. The continuation of the old prefix without a distinguished donor remains a third possible summand. There is no branch in which this child fills both roles. When its eventual path child is older than `y`, the already-registered extension gets `y` as a later donor birth. When `y` is older, the already-distinguished shorter prefix carries it into subsequent extensions. Fixed labelled paths and donors identify their unique leaves and promote at unique times, so induction over the arrows gives one registration of each final incidence. The two coincident-branch cases are explicitly exercised in the arrow referee.

For every prefix, `k` is the total physical core-child count at that vertex, including unregistered, older and off-tag children. Known path and donor children each subtract one from the unseen receive count; the short tip does not subtract a child for its still-future extension. If the same tip also owes an unattached donor, it needs two distinct future children, so `k+2≤arity` there. At a promotion that fulfills one obligation, `k` rises by one and the fulfilled demand falls by one. Thus failing a degree or `k+demand` cap for that registration is irreversible, even when a physical child later fulfills one of its obligations. No old child can be adopted retroactively into a registration that was killed.

**Reward eligibility is different from the state cutoff.** A live four-tuple can have zero current reward because it is still in the excluded old band or is still too small to clear the frozen reserves. It can later acquire a price. The referee takes `(3,1),(5,3),(3,1),(5,0)` from old-band exclusion to eligibility by promoting an ordinary leaf of `C` to the new distinguished child; it also feeds a frozen-size descriptor to a newly eligible one. A separate over-degree fixture remains killed. A paid descriptor is not resurrected from an earlier degree-cap exit. Donor and mother overflow are also irreversible, but are retained as live reciprocal intervals rather than killed.

## Disjoint support

Let `j=Σ_i k_i−3`. Each of the `j` exterior core-child branches is disjoint from the path, is below `A`, and has at least two vertices. The true `A` fringe therefore has at least `Σ(D_i−k_i)+2j` vertices excluding its mother, and at least `4+j` cores. Keeping only lower bounds on size greater than 20, or greater than 12 with `j>0`, excludes the size-12 and three/four-core size-13–20 higher-family payments. The separate `g/hC` and residual `f/g` payments do not use raw P4. The previous band exclusion checks all four degree caps and the complete old arity pattern, not only the subset actually priced by the previous unpaid objective. This removes the whole prior band on present-time physical descriptors. Changing the method or the distinguished-donor memory cap alone does not make two prices for the same raw support additive.

## Evidence

From `[local path]`, both commands passed with assertions enabled and original payload hashes checked:

```sh
P=round2/.venv/bin/python
OPENBLAS_NUM_THREADS=1 nice -n 15 taskset -c 3 timeout 85s "$P" -Bu \
  round2/checks/p4-lift-sol.py b3f18_lo.bin b3f18_hi.bin
OPENBLAS_NUM_THREADS=1 nice -n 15 taskset -c 8 timeout 100s "$P" -Bu \
  round2/occupation-neon/pressure-memory/donor_lift.py \
  --lower b3f18_lo.bin --upper b3f18_hi.bin
```

- The referee (`p4-lift-sol.py`) checked all ten orders, 256 independently enumerated Bellman corner boxes, 4,050 mask descriptors, and 531 projected row probes with 2,715 grouped successor/credit keys. It probes the pinned producer row with symbolic successors and compares to rates obtained independently by enumerating directed arrows on the labelled tree. This checks bounded hidden hazards and shared physical credit on the same successor, including short prefixes and both exact and overflow memory.
- On 32 labelled trees (64 signed raw score comparisons), the independent whole-tree potential evaluator and two centered incoming transfers reconstructed the score from 266 signed donor terms. Six slot instances had multiple nonvirgin donors, four had side-port sums at least `5/6`, and omitting a second donor changed 23 signed comparisons. No independence assumption is used for the samples or the proof.
- The matched six-part producer replay checked 216,294 live states across components and 865,176 signed drift inequalities and reproduced both exact floors above. The reviewed `donor_lift.py` SHA-256 is `b9423487db066349cdb8bbd7a64380041c65fe2fae14d05907a4932775bd195b`. The referee pins the earlier band's companion and its inputs, and does not modify them.

## Medium addendum: `(9,7,6,10)` / `(3,3,3,2)`

**Correct in scope: medium physical accounting and exact composition.** The medium run uses the same pinned source and mask predicates, side and mother memory `3/3`, budget `3,000,000` per part, and all six parts. It replaces the smaller complementary price. The producer receipt (`occupation-neon/pressure-memory/p4-donor-medium.log`), SHA-256 `5302e1414a632177f1254001854adc6c171aeb2a0df559cd4924d67b01730994`, reports the completed Bellman run. Those same fractions appeared in the live receipt `astra-p4-lift-97610-3312.log`.

The independent referee reran the labelled and signed comparisons under the medium caps: **98,280** candidate masks, **91,794** eligible supported descriptors, **669** projected physical rows, **4,009** grouped successor/credit keys, and **64** signed raw score reconstructions with **358** signed donor terms. It exercised **32** slots with multiple nonvirgin donors and **8** with combined side pressure at least `5/6`; omitting a second donor changed **47** signed comparisons. The maximum-cap trees exercise up to two off-path core children at each tagged vertex. The same ten birth orders and 256 independent corner boxes also passed.

The exact recombiner read all six component blocks in order, checked the mode/caps, counts, and nonnegative interval orientation, summed all 24 exact price endpoints by witness and direction, and matched the printed totals. It counted **6,017,604** live states and **24,070,416** drift inequalities in that receipt. The audit invocation with `--rerun minimum mother` finished the independent fixture and arithmetic but did not report completion of either selected Bellman part; no medium Bellman component is claimed as freshly rerun by this reviewer.

Separately, efd8 ran all six components to completion in an external `ead947b` archive with explicit overlays for `path_band.py`, `path.py`, and `donor_lift.py`. That replay exited zero after 437.78 s wall. The archive receipt `donor-lift.log` has SHA-256 `a1561153bc9ec3f36905295d74aea094fc333d4e09b13294bb895f67fb5a06f9`. A separate `python3 -B` receipt comparison checked all six per-part counts and all 24 price endpoints exactly against the producer receipt, then repeated the rational composition.

| Witness | Exact stationary floor | Exact upper bound on this support |
|---|---:|---:|
| Lower | `359297001126132137/10^24` | `8788537130454641/5000000000000000000000` |
| Upper | `6963149521191559/25000000000000000000000` | `32373738731142851/20000000000000000000000` |

To compose, start at the exact endpoints in `ead947b`; divide by the full corrected EDGE denominators `1357493815758092773/10^18` and `273436004214812983/(2·10^17)`. The medium-only outward pair is `[0.565678321308780, 0.565770265677644]`. Adding the disjoint prior wide-band floor as well gives **`[0.565678353835574, 0.565770240831497]`**. This is the P4 subtotal for efd8's freeze; no other motif credit enters these fractions:

```text
137071018766641883155149709118897/242312646112819559980500000000000
6574833043156219378139076805857/11621030179129551777500000000000
```

Run the independent accounting and receipt recombination together from `[local path]`:

```sh
P=round2/.venv/bin/python
OPENBLAS_NUM_THREADS=1 nice -n 15 taskset -c 3 timeout 150s "$P" -Bu \
  round2/checks/p4-lift-sol.py b3f18_lo.bin b3f18_hi.bin --medium \
  --receipt round2/occupation-neon/pressure-memory/p4-donor-medium.log
```
