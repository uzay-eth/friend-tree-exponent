> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Physical birth-parent law: joint degree, arity and sided flux

**Matched lower `hV`:** the same lifetime price vector now gives **`≥5.229999e−9` per core**, versus the booked `3.06204334292295e−9` (about **71% more**). The new LP reduces the possible mixture range for that vector to **`[5.229999e−9,5.550233e−9]`**. Its width is `3.20234e−10`, vs `7.85868e−9` with an arbitrary mixture of the same conditional prices and bracket discount normalization (about 96% narrower). This is a bound on the expectation of the certified **lower-price vector**, not a two-sided enclosure of the true family reward. No numerator term was added to an existing ledger.

## Law and exact rows

Work on the physical oriented stationary tree. For the ordinary root law `π`, set `x=π/μ`; for the joint parent–child edge incidence law per core `E`, set `e=E/μ`. The birth-mother law in **pre-promotion** `(D,k)` is `β(f)=x(bf)` with `b=(D−k−1)/D`. It is not the uniform-root law or the existing-edge `k` tilt. For a post-promotion mark use `f(D,k+1)` in the same formula. One birth creates one real core with census `(2,0)` and leaves the mother's degree fixed. The rescaling keeps the denominator exact without choosing a separate `μ` for each objective; in particular `x(b)=1`. A true law satisfies

```
x(k−1) = 0,   x(b−μ) = 0,   x(ρ↑+ρ↓+μ−1) = 0,
x[ μ(V−V(2,0)) + (L+ρ↑+ρ↓)(V−V_gain) + b(V−V_promote) ] = 0.
```

Here `ρ↑=Σ_children 1/D_child` and `ρ↓=1/D_mother`. The source is physical stable age, and the sided constraints refine the existing [verified joint-edge flow](../../flow/README.md): the **same** `e` has outgoing count `x k`, incoming count `x`, outgoing parent reciprocal `x k/D`, and incoming child reciprocal `x/D`. Conditional on endpoint bins, its child-reciprocal row marginal is `x ρ↑` and its parent-reciprocal column marginal is `x ρ↓`. Lower and upper port bounds use those opposite endpoint reciprocal coordinates, not the same vertex's own `1/D`. The two unweighted normalized fluxes are distinct: `x(ρ↑)=x(1/D)` and `x(ρ↓)=x(k/D)`. They sum to `(1−μ)/μ` under `x`.

A birth-tilted **receive intensity** on a set of mothers is `x(b ρ^side 1_cell)`. On retained census cells `b` is a constant, so it is a linear objective compatible with the sided flow rows. This is an integrated observation about actual old children or parent at the time of a real birth. It is not a bound on conditional pressure in a later stage of a P4 or fork path. In particular we do not assume that the edge table itself is birth tilted. The code lifts `μ` to bracket endpoints in each column (an outward relaxation of common `μ`), and interpolates each physical port rate separately over its exact box `[0,k/2]` or `[0,1/2]`. The true lift preserves the required affine moments; `b ρ` objectives here are only evaluated on finite degree/arity classes with exact `b`. For a future genuine tail birth-flux objective one would additionally need to lift that product or apply a proved tail cap.

Overflow for degree keeps arities `0,…,6` and `≥7`, not a cemetery. Incoming high-arity child flux and excess core stubs go into nonnegative columns; the weighted outgoing stub identity is still `kp=1−b−p`. The port table uses exact degree ≤12, then reciprocal intervals within `13–18`, `19–cut`, and `>cut`, and arity buckets `0,1,2,≥3`. Each true edge lifts by the two endpoint reciprocal weights in the very same table. All column-class masses are at most `1/LO`: occupation and edge each have mass `1/μ`, excess stubs at most `x k`, tail child receive at most `xρ`, and combined Palm slack at most one. This also validates the exact dual repair: round all proposed multipliers to `10^-12`, inspect every rational reduced-cost inequality, subtract each class's largest nonnegative defect times `1/LO`. Neither a floating primal solution nor a float objective is used as evidence. HiGHS dual simplex is a fallback for an interior-point `Unknown` status.

The lifetime birth constraints `β≥β A` from [the independent P4 audit](../../../../../../../reviews/cold/astra/checks/stationary-p4-sol.md) are imposed alongside the static census and sides. They are directly linear in the birth occupations. Each row of the physical direct-production matrix `A` includes all degrees through 34, with nonnegative tails omitted only in this lower matrix; it uses `HI` and the appropriate monotone cap on birth-parent pressure. This is the physical class kernel, not residual killing or Perron weights.

## Birth-law diagnostics

| Quantity | Scalar, cut 18 | Shared sided ports, cut 18 | Shared sided ports, cut 24 |
|---|---:|---:|---:|
| `β(D=2)` | `[.217545126367,.242051863970]` | `[.217776309272,.242051863970]` | `[.217785391885,.242051863970]` |
| `β(D≥13)` | `[.263815974008,.355189364367]` | `[.273111811000,.338356047416]` | `[.276529282910,.335406889385]` |
| `β(k≥2)` | `[.234631412730,.284765667286]` | `[.241625129359,.282844362987]` | `[.244246720404,.280794141298]` |
| `β(D≥13,k≥2)` | `[.170337426541,.242552675079]` | `[.179636173372,.237558189257]` | `[.185109153887,.234876313767]` |
| Child receive under birth tilt, `13≤D≤cut,k≥2` | `[0,.040818356556]` | `[0,.038790791121]` | `[0,.067282584084]` |
| Mother receive under birth tilt, same mask | `[0,.016973221409]` | `[0,.015809879189]` | `[0,.026294111229]` |

The receive comparison across cuts changes the objective's mask; only the first two columns are matched for those two rows. The prior aggregate `πρ` flow diagnostics cannot be applied as pointwise rates. Example: at cut 18 actual separated child/mother flux sharing shrinks the high-degree high-arity birth-mass width by about 20% relative to scalar flow, and its finite child/mother flux caps by about 5%/7%.

## Lifetime pricing: same support, no whole-row stacking

* **Lower `hV`, `18≤D≤32,5≤k≤7`, both non-apex tags nonvirgin and degree at most three:** at cut 18 the scalar LP gives `[4.994005e−9,5.372885e−9]` and the sided LP `[5.018668e−9,5.301858e−9]`; at cut 24 use `[5.229999e−9,5.550233e−9]`. The old five fixed class floors gave `3.06204334292295e−9`. The shared full-vector price is a replacement, not additional `hV`. The replay uses `LO=.565673` (thus reports the old floor in this run as `3.062040e−9`; the old source uses `.565673570322`). `--fresh` reproduces all five exact conditional prices from the pinned witness payload, without a new payload hash or full GMP replay.
* **Upper `hV`, `18≤D≤64,5≤k≤8`, same tip mask:** at cut 24 `[6.573014e−9,6.766145e−9]`. This raises the booked unstratified `4.00077770896578e−9` floor about 64%; reusing only the five old class floors on the newly completed vector would give `5.233545e−9`. The unconditional upper family still has no full-reward ceiling here.
* **Lower fork `Yt`, degree `18…36`, arity `5,6`, nonvirgin tags of degree at most three** (the source's conservative half-payment): cut 18 gives `[1.1567e−11,1.2304e−11]` for the given vector, compared to the source's `1.49696e−12` universal floor or `6.527e−12` with just old class floors. Absolute improvement is much smaller than `hV`. This is the available old-flow fork forecast; the other writer's broader complete P4 seed masks were left untouched.
* **Specified high-spread P4 cell**, caps `6/12`, same conditional score intervals as the stationary P4 review. At multiplier `1.4`, the old Palm-only contrast `[1.58154e−10,2.813838e−9]` tightens to `[1.76491e−10,2.774745e−9]` at LP cut 24. At `1.5`, the interval is `[-1.91693e−10,2.552074e−9]`. The contrast is affine in the multiplier for the physical expectation, so exact interpolation at `1.44` gives `[2.9218e−11,2.685677e−9]` and proves **stationary odds >1.44** for this cell. A straddling `1.5` interval does not prove an upper odds bound.

## Denominator uncertainty

The inherited witnesses have `πΦ₋ ∈ [.28366341151,1.356833021311]` and `πΦ₊ ∈ [.28865709309,1.367169995210]` (GMP lower comparisons and certified stationary uppers; see [denom](../../denom/README.md) and [donor-flux controls](../../higher/donor-flux/memory/README.md)). Thus even the tightly priced lower `hV` vector normalized by `πΦ₋` lies in the wider `[3.854564e−9,1.9566261e−8]`; upper `hV` likewise lies in `[4.807751e−9,2.3440078e−8]`. For a **fixed replacement difference of certified floors**, dividing by the stationary ceiling gives only about `1.598e−9` (lower) and `1.881e−9` (upper, versus booked universal) in endpoint units. These are scales, not committed endpoint moves. Using the respective floor rather than ceiling changes each normalized number by a factor up to about `4.78` or `4.74`; the physical `μ` normalization varies by only about `0.018%` over the bracket. No equality, joint bound on the two denominators, or separator for the total signed drift follows from the two narrow mixture prices. No expenditure was added to any existing whole-row spend.

## Replay

Let `R=round2/occupation-neon/pressure-memory/birth-law/replay.py` and `PY=./round2/.venv/bin/python`. Commands below use `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 13 nice -n 10 /usr/bin/time -f 'cpu=%U+%S wall=%e'`, with bounded timeouts. Before the final diagnostic print-only edit the matched runs were:

| Arguments after `$PY -Bu $R` | Actual result | CPU user + system |
|---|---|---:|
| `hV --limit 18 --fresh` | Exact five conditional inputs replayed; sided lower floor `5.018668e−9` | `9.43+1.62` |
| `hV probe --limit 18 --scalar` | 194 rows, 1,426 columns; lower `hV` floor `4.994005e−9` | `3.35+0.05` |
| `hV p4 probe --limit 24` | 613 rows, 5,554 columns; both P4 base contrasts, lower `hV` and all birth-law diagnostics passed | `22.69+0.27` |
| `p4 --limit 24` | Both contrasts and exact interpolation at `1.44` passed | `11.84+0.24` |
| `hV --side upper --degree 64 --arity 8 --limit 24` | All five wide upper-witness conditionals passed; floor `6.573014e−9` | `16.04+3.24` |
| `hV Yt --limit 18 --degree 36 --arity 6` | Both matched vectors passed, Yt floor `1.1567e−11` | `22.21+4.09` |

An earlier parallel attempt on CPU12 completed the five lower tagged passes but hit its wall limit before LP output; the CPU13 replays above supply the certificates. The raw D18 payloads are prior authenticated inputs; no new denominator or GMP scan was run.
