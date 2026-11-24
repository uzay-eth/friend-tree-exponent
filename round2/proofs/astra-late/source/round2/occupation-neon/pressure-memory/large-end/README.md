> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Large physical endpoint sends in the D18 P4 and Yb rows

The two kinds of endpoint terms have different truncation costs. For the **internal** endpoint-to-neighbour transfer in the paired D18 row, the smallest rectangle in the tested envelope for a `1e−6` omitted-price budget is `(degree, arity) = (142611, 40)`. For the **exterior** child-to-endpoint gain term it is `(4097, 40)`. These are certified bounds on omitted terms; the occupations inside either rectangle have not been computed. Four explicit unspent physical cylinders have positive full-row prices for both witnesses, but these particular lower bounds are too small to justify an endpoint move.

## Row and ledger

Use the original `b3f18_lo.bin` and `b3f18_hi.bin` from the [D18 verification](../../higher/README.md), opened read-only. `δ_lo` and `δ_hi` are the separately nonnegative signed charges (`+drift` for the lower witness, `−drift` for the upper). The paired row here is `δ_lo−δ_hi`, at multiplier one. Slot orders are `A→B→C→T` for `P4` and `(A,S,B,T)` for `Yb`, with `A→S`, `A→B→T`. In both blobs, `P4` slots `B,C` and `Yb` slots `A,B` must be in their finite table supports (`D≤9,k≤7`). Thus these slots cannot be large ends. `D=17` is already outside the ≤16 inner-gain calculation, but it has its own gain seam, so this replay starts at `D=18`.

| Large end | Internal send, with connected-deletion owner | Maximum compatible motifs per end |
|---|---|---:|
| `P4:A` | `A→B`, from `hV(B,C,T)` | `49 k_A` |
| `P4:T` | `T→C`, from `hV(A,B,C)` | `1` |
| `Yb:S` | `S→A`, from `hV(A,B,T)` | `42` |
| `Yb:T` | `T→B`, from `hC(A,S,B)` | `6` |

These counts follow from the two arity-seven slots and the unique actual ancestry: for the `Yb` short tip there are at most six other children of `A` and seven choices beneath each; for the deep tip there are only six alternative short children. For an exterior child send to the large end multiply the corresponding count by at most `k_end`. All counts refer to physical incidences, even when ports have equal census labels.

The internal part is `(1/D_end−center) ΔW_deleted` in signed units, where the deleted row reads the **receiver's** gain derivative. Only the `1/D_end` term is the sender intensity budgeted below. The end's other sends already enter the genuine within-motif rate. The exterior child part is `ΔW_full Σ_child,ext 1/D_child`, since four-body rows do not route exterior sends further. The paired potential differences in these two coefficients use the **sum** of lower and upper payload weights before the drift, because the upper charge already carries a minus sign. The discount term must still use each payload's own trial. Neither endpoint channel can by itself be booked as an additional nonnegative full-row charge.

For a concrete unspent subset, insist on a nonfounder apex and `Σ_tag k≥4`, giving at least five cores in the complete fringe, regardless of unobserved descendants. `D_end≥18` forces the complete planted fringe above size 12. This excludes size-12, the ≤16 size-13–20 suffix, and `arms` (≤4 cores). For `P4` additionally require apex `k<5` and some non-apex `D>5`. The arity condition rejects the old-flow, persistent old-motif and high-apex bands; the degree condition excludes the staged `hV/P4/Yt` domains and the four-bare-degree-two `P4`. `Yb` does not belong to those tagged `P4/Yt` additions. The broad bounds below drop these exclusions only when taking an upper bound on an absolute omitted term. In the corrected EDGE ledger only `R_f,R_g` are confirmed unspent for its own objective, so this higher-row result does not amend EDGE.

## Exact degree interpolation and a positive support control

Fix the other labels and the end's actual `k,h`. At `t=18/D`, each payload weight is cubic in the variable end, using its Bernstein controls, including gain and promotion successors. The gain uses `t'=18t/(18+t)`. After multiplying the affine row by `(18+t)^3`, the gain numerator has degree at most six and vanishes at zero, cancelling the possible `1/t` in leaf stock. Other degree-dependent multipliers are affine in `t`. Hence the cleared row has degree at most seven. The connected-deletion transfer's coefficient is constant in the end's degree. For a global exterior slope bound, the derivative of a cubic has magnitude at most three times the maximum absolute adjacent control difference. Thus `|ΔW_full|≤C/[D(D+1)]` with `C=54 max|W_{j+1}−W_j|`, including the full overflow arity control. The other physical tail slot, if present, is a convex Bernstein mixture of its own controls. The replay scans all compatible deleted-row receiver gains for the internal coefficient ceiling, and all tail-control slices for `C`.

The probes keep the same physical exterior donor assignments as degree varies: the apex parent has degree eight, and there is one child of degree three outside the motif at the end. At degree 18 the other labels are `P4:A` = `B(4,1,1), C(1,1,1), T(1,0,1)`; `P4:T` = `A(0,1,1), B(1,1,1), C(1,1,1)`; `Yb` = the non-variable entries of `((1,2,1),(1,0,1),(1,1,1),(1,0,1))`. The endpoint receives one more child than required by the motif topology. Values displayed here are approximate; the computations are rational.

| End | Paired row at 18 | Exterior child slope at 18 | Internal transfer coefficient | Forward curvature `q(20)−2q(19)+q(18)` |
|---|---:|---:|---:|---:|
| `P4:A` | `+1.34634e−8` | `−2.08651e−7` | `+2.89431e−4` | `+3.28644e−10` |
| `P4:T` | `−2.89869e−4` | `−2.72040e−5` | `+1.53046e−2` | `+1.63143e−7` |
| `Yb:S` | `+8.42787e−6` | `−8.78240e−7` | `+7.82641e−4` | `−5.25552e−9` |
| `Yb:T` | `+6.65599e−6` | `−2.57066e−6` | `+1.67589e−3` | `−7.73476e−10` |

The support is positive under the true joint law: from the newborn apex, prescribe one legal sequence of its intrinsic gains and promotions and those of its actual successive descendants, with no other intrinsic fringe jumps along that sequence. Multiply each desired physical rate by the inverse of the total intrinsic jump rate plus `HI+1/2`, then pay a final holding interval with the same denominator. This lower-bounds the discounted physical probability/occupation uniformly even for an adaptive real parent send, and multiplication by `S≤μ` gives a stationary lower bound. The parent port is allowed its entire real range in the reward. For example the `P4:T` cylinder has mass at least `2.1436e−13` per core and full-row lower bounds `5.1980e−19` and `9.0977e−19` in the lower and upper witnesses. This proves attainable unspent support rather than inferring mass from a Cartesian product of stationary census marginals. No price is added to an exponent ledger here.

## Common-age Yule/Poisson bound and cap

Move the sum to each actual end by incidence reversal *before* taking a stationary expectation. Under the [physical birth-Palm identity](../README.md#an-effective-joint-degree-and-child-count-tail), the end's age is exponential with rate `μ`. At fixed elapsed age use independent **proposal** processes `Y` and `N` to dominate its degree and arity on the same trajectory; keep the common age when integrating. Specifically `Y` starts at 2 and gains at rate `d−1/2`; `N` is the end's own unit-rate send clock. Actual edge sends and actual promotions are thinnings conditional on the full past. This does not make physical motif vertices, degree, or arity independent. Couple the observation age below `Exp(S)` with `S=.565673`.

For omission by degree (`D>J`), `π[(k)_m 1{D>J}]≤P_J E[(T')^m]` where `P_J=∏_{d=2}^J(d−1/2)/(d−1/2+S)` and `T'` is an independent `Exp(S)` plus independent exponentials of rates `d−1/2+S` for `2≤d≤J`. Indeed, integrate `N`'s factorial moment above the hitting time for `Y` and use exponential tilting plus memorylessness. For omission by arity (`k>K`), the proposal `N(Exp(S))` is geometric and gives exact polynomial-weighted tails via its memoryless shift by `K+1`. Use `D≥max(18,K+2)` on the arity overflow, `D≥J+1` on degree overflow, and add the two upper bounds. This is an omitted-region union bound, even if the two conditions overlap. It covers the weighted coefficients and all remaining physical correlations, but does not infer the occupation of retained rows.

The objective here is the minimum rectangle label product `(J−17)(K+1)` meeting a specified remainder in this envelope; it is not a claim about an optimal physical cap or actual storage after support pruning. The code searches `18≤J≤262144`, `0≤K≤96`, binary searching `J` only after checking monotonicity of the rounded degree bounds. Table remainders are displayed rounded up.

| Omitted term | Budget | J | K | Rectangle labels | Upper remainder |
|---|---:|---:|---:|---:|---:|
| Raw internal send mass | `1e−3` | 5425 | 28 | 156832 | `9.99995e−4` |
| Paired internal transfer | `1e−5` | 30281 | 35 | 1089504 | `9.99968e−6` |
| Paired internal transfer | `1e−6` | 142611 | 40 | 5846354 | `9.99996e−7` |
| Paired internal transfer | `1e−7` | — | — | — | Not met on grid |
| Paired exterior child term | `1e−6` | 4097 | 40 | 167280 | `9.99678e−7` |
| Paired exterior child term | `1e−7` | 10530 | 46 | 494111 | `9.99966e−8` |
| Paired exterior child term | `1e−8` | 27524 | 51 | 1430364 | `9.99976e−9` |

The internal channel decays only like a reciprocal degree, whereas the exterior slope has a gain difference of order `D⁻²`. The separate omission caps do not budget the entire signed row: center, discount, promotion, and all other role gains need matched treatment for that objective. The next useful calculation is a joint occupation of the **retained** weighted incidences, using the actual endpoints and their birth sources; a tighter internal remainder can also bound coefficient and arity jointly rather than use their separate maxima.

## Replay

```sh
taskset -c 8 nice -n 15 env OPENBLAS_NUM_THREADS=1 round2/.work/venv/bin/python -Bu round2/occupation-neon/pressure-memory/large-end/replay.py
```

Checks include matching the independent `Local.charge` evaluator for both blobs at eight interpolation and three held-out degrees per role, zero constant and linear terms in the cleared exterior slope, 870 chronological core trees for the incidence counts, and three exact small-cap common-age references. The path defaults match the existing verified payloads; this replay does not stream new hashes for the multi-gigabyte inputs.
