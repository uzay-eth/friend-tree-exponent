> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# All-degree stationary denominators for Neon D18

Here `πΦ` is the limit of the normalized mean additive potential per
physical core, with a separate potential for each sign. The stationary
upper bounds are rational, on the full degree domain:

| Witness | Global `C0` | Certified upper on `πΦ` | Reduction |
|---|---:|---:|---:|
| Lower | 1.84124847315 | 1.356833021311 | 26.30% |
| Upper | 1.95820675182 | 1.367169995210 | 30.18% |

Holding the earlier [`check`](../../../../../../reviews/cold/astra/occupation-neon/check/README.md) numerator certificates
fixed gives the outward interval **`[0.565674816769, 0.565773546690]`**.
The inward changes from those endpoints are `0.000000477975` and
`0.000000438646`. This replay does not use new numerator work or iterate
the exponent bracket.

## All-degree charge

Use the accepted physical eight-motif decomposition and the independent
Neon audit binary layout. In raw units let `q(L,k,h)` bound the total
charge allocated to a core. Allocate the edge, P3, P4, and bend terms to
their bottom vertices. For a bend, that vertex belongs to at most `k_u−1`
terms for its fixed grandparent `u`; maximize that multiplicity times the
control value over the active grandparent descriptors, including zero.
Allocate the child-pair, Y-tail, and triple terms to their branching core,
using their exact pair/triple multiplicities. Their controls are maximized
over the other slots. Only the allocated core's Bernstein weights remain;
these are used exactly for its census. For the global cap, the bound is
the maximum same-census sum on all stored leaf rows and an outward convex
tail cap. The lookup clamps the first-body leaf index at `LF`, and the
rest of the charge uses the proper tail or vanishes past its arity cap.
The arity overflow class uses the first overflow arity for recharging;
larger arities only decrease the bound when the recharge `t≤0`.

The exact tree identity `sum(k−h)=0` allows the signed recharge
`t(k−h)` at each core. The script sweeps a rational grid of valid `t`, takes
the corresponding all-degree cap `M`, and uses the *deficit*
`M−q(L,k,h)−t(k−h)≥0` as a reward. In particular, a bin with a loose
charge is never treated as if its potential were known exactly. The
resulting finite-tree inequality is `Φ ≤ M C − sum(deficit)`, in normalized
units. We omit the founder and all deficits outside the degree-64 cut;
neither omission can make the bound optimistic. Each allocated uncharged
potential is a bounded local score because motif multiplicities are capped
as above. The mean (and Laplace mean) per-core limit thus follows from
the oriented renewal limit; equivalently take the Laplace limit directly
in the finite-tree inequality using the bounded truncated deficit.

## Two independent occupation payments

1. For each physical newborn, solve the root census backward Bellman
   inequality to degree 64 with discount `HIGH`. The real gain intensity is
   between `L` and `L+(k+1)/2`, with promotion rate `L/D`. The fractional
   linear minimum need only check the two endpoints and is valid for an
   adaptive physical exterior. The discount is above `μ`; multiply the
   downward-rounded lifetime reward by `LOW`, below the stationary birth
   rate. The resulting deficit credits are at least `.459897103759`
   (lower witness) and `.462522415600` (upper witness).
2. Independently expand the genuine capped nineteen-label residual
   occupation with size 8 and depths 0–10 at `HIGH`. All omitted terms
   are nonnegative. For any normalization with `ell·1=1`, the critical
   left Perron coefficients satisfy `ell K=ell`. The stationary normalized
   root occupation is `ell G f / ell G 1` and `ν(b)=μ` by physical birth
   compensation. Consequently `ell G 1=1/μ`; multiplying the worst of
   the nineteen certified rows by `LOW` is valid without fixing the
   unknown Perron weights. This gives independent credits of
   `.328840332095` (lower) and `.328684023120` (upper). These prefix
   credits are weaker; the two credits **are not added**.

The signed drift argument uses the unchanged earlier guarantee
`σ(λ−μ)πΦ≥reserve`. Substituting a smaller stationary upper bound on `πΦ`
is legitimate because the original sign of `σ(λ−μ)` is already
certified. The displayed exponent endpoints reuse only the older rounded
gain multiplied by its original `C0`, a conservative numerator floor.

## Replay

Inputs are the same accepted payloads as in [`rows`](../rows/README.md)
and [`check`](../../../../../../reviews/cold/astra/occupation-neon/check/README.md); this does not rerun their full GMP
inequality replay. Six labelled test trees per witness exercise all eight
physical motif contributions and check their sum against the local charge
sum, including the leaf and arity overflows. The binary caps and integer
reduction overflow are guarded by the independent audit and this script.

```sh
OPENBLAS_NUM_THREADS=1 taskset -c 8 nice -n 12 /usr/bin/time -f 'wall=%e cpu=%U+%S' "$PYTHON" -Bu round2/occupation-neon/denom/check.py "$D18_UPPER" "$D18_LOWER"
```

The run took `57.52+0.66` CPU seconds. The selected normalized recharges
were `−1/500` (lower) and `−9/1000` (upper), with all-degree caps
`1.81673012507` and `1.82969241081` before occupation credit.
