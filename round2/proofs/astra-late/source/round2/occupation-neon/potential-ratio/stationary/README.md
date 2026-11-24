> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Physical stationary ratio of the paired D18 potentials

**Exact enclosure:** `0.9995 ≤ πΦ_hi/πΦ_lo ≤ 1.085`. Both ends improve the [global pointwise comparison](../README.md) `[.999,1.12]`; the width falls from `.121` to `.0855` (29.3%). The `.997` side of the requested bracket is covered; the proposed `1.02` upper is not certified. Finite shape histograms are not used as stationary masses. The global comparison stays unchanged (its physical wide-crown fixture already exceeds `1.064`).

## Signed comparison on a shared physical law

For trial `q` use `X=sign·(Φ_hi−q Φ_lo)`, with `sign=+1` at an upper trial and `−1` at a lower trial. Normalize the potentials only once by their common raw scale. Compare the ratio of the two stationary means, not the mean of the pointwise ratio. The lower witness has a strictly positive stationary mean by its original core-count lower comparison.

The existing exact block reducer computes integer controls for this *one* signed source. The layout and both SHA-256 pins are checked before using the payloads. `f` reads its exact joint signed rational value at the same leaf count and arity. `g` starts from a signed cap conditioning on its actual bottom type, and the remaining six families use signed control bounds with the original support masks, shared Bernstein coordinates, and correct physical multiplicities. Thus nearby weights cancel before extrema. We never compare the two individual denominator ceilings to get the ratio.

Let `B(d,k)` be the resulting nonfounder census cap. Edge, vertical motifs and bent motifs are charged to the bottom; sibling pair and triple go to their apex, and Y-tail to its branching middle. The founder exception in this inequality and in transporting gaps disappears in the normalized physical mean: only bounded depth from the founder is involved, and its expected count there grows at most polynomially (each core creates core children at rate at most one). Bounded-row score allocations and their mean limits are as in the [denominator proof](../../denom/README.md).

## All-degree upper and coupled physical gap credit

1. Add a rational `t(k−1)` to the bound at each nonfounder. This is mean-zero by physical mass transport (one parent and mean one child), or by exact telescoping on finite trees with the founder correction. For `d≤48`, contract the **actual** leaf/child census and the allocated endpoint's Bernstein weights. Above that degree, use an all-degree cap, maximizing over legal stored leaf counts and all control groups, with the `f` leaf cap and the exact overflow law. For `k>180` the underlying cap is constant and the additional slope `t≤0` makes the first overflow arity sufficient. No tail mass is dropped or assigned the finite sample's ratio.
2. Bound the discounted integral of the recharged census cap at a physical virgin by a backward supersolution. The promotion rate is exactly `(d−1−k)/d` and the total gain rate is in `[L,L+(k+1)/2]`, allowing arbitrary adaptive exterior. At the degree exit pay the constant `max(M/LOW,M/HIGH)` for the all-degree cap `M`. Test both gain extremes and both discount endpoints on the same rounded value, including signed continuation at exit. This certifies a lifetime upper `U`, hence per-core upper `max(LOW·U,HIGH·U)` by the positive [birth-Palm identity](../../../../../../astra/source/round2/degree-law.md). Float search chooses only the recharge; every supersolution row is checked with rational arithmetic.
3. The difference between each signed family cap and its *actual joint occurrence sum* is nonnegative. Transport these differences along the actual ancestor paths to an owner that sees the whole group. In particular the edge deficit is moved to the actual parent and the Y-tail deficit at the branching middle is moved to its parent. A bend is moved to the grandparent of its bottom so all eligible distinct sibling choices stay together. Vertical P3/P4 go back two/three levels. The sibling pair and triple stay where charged. This transport preserves the sum, with no extra charge for a donor or for the same motif at another owner. First-body already reads exactly and contributes no gap. In a planted fringe of size at most 11, contract the physical values of **all** the groups simultaneously from both original potentials, and assert each gap nonnegative. Rearrange no rates in the potential score itself. In the physical newborn occupation of the gaps the full internal jump sequence does retain actual parent/child feedback: every internal source uses its true reciprocal degree, all promotions change the actual tree, and the same outside parent rate in `[0,1/2]` applies to the root gain. The shared rate enters that one first-event row, including the actual gain successor. Uniformly lower-bound this nonnegative discounted gap at `HIGH`, kill at the size exit with value zero, and multiply by `LOW`. Every omitted group/window stays covered by the initial cap. No residual-renewal copy valuation is added to this physical occupation, and no drift slack or donor reward is added anywhere.

| Trial | Census upper before transport credit (approx.) | Certified credit (approx.) | Final signed upper (approx.) |
|---:|---:|---:|---:|
| lower `0.9995` | `3.90348e−5` | `6.11609e−5` | `−2.21261e−5` |
| upper `1.085` | `0.00968551` | `0.01027324` | `−0.000587732` |
| upper `1.02` | `0.06731683` | `0.00988985` | `+0.05742699` |

A positive upper is **unresolved**, not a physical counterexample. The coarse census envelope has a large `1.02` remainder even after the exactly observed gaps; crediting just the `g` edges resolves `1.09` but not `1.08`, and does not resolve `.9995`. Further narrowing needs conditional restrictions on the *remaining* signed joint potential terms at older parent/child census, particularly mixed/tail groups; the marginal maximum over other motif roles is still retained there. One option is a live parent-degree-shell table for joint `g` and vertical terms, with signed Bernstein support kept on the same reciprocal and exact stationary edge-flow projections. Nothing here claims the true ratio exceeds `1.02` or changes an exponent endpoint.

## Authenticated replay

From this directory:

```sh
timeout 110s taskset -c 7 nice -n 10 env OPENBLAS_NUM_THREADS=1 \
  ./round2/.venv/bin/python -Bu check.py \
  b3f18_lo.bin b3f18_hi.bin .9995 1.085 1.02 \
  --cap 48 --edge-cap 11 --full-credit
```

Both hashes passed, both enclosing signs passed, and `1.02` was unresolved as above. The stopped solve cached 7,812 geometries including its exit layer. Decimal output is display only; the sign checks use rational numbers. The all-degree scan is of **potential coefficients**, not of `LΦ−λΦ` charges.
