# Cold review: stationary quintet joint certificate

**Holds on both signs at cut 24.** Independent rational repair reproduces both exact endpoint fractions in the raw joint receipt (`occupation-quintet/stationary/joint24.log`). The outward interval is **`[0.565703060655196, 0.565746333629979]`**. A fresh solve at each outward endpoint also certifies a strictly positive inward trial margin.

| Exact rank checked before rounding | Lower, outward | Upper, outward |
|---|---:|---:|
| Original production receipt | 0.565702177913 | 0.565748407438 |
| Historical host-only, cut 24 | 0.565702368490575 | 0.565748281318382 |
| Matched cut 24, release both edge credits | 0.565702661777197 | 0.565747732733322 |
| Full cut-24 joint certificate | **0.565703060655196** | **0.565746333629979** |

The rank is strictly increasing lower and strictly decreasing upper, in exact arithmetic. At the stated outward trials the fresh repaired signed margins are at least **327.021697154 / 364.273166213 ticks** (lower/upper). Even just transporting the pinned trial certificates gives strictly positive margins `0.001237743365 / 0.000509888245` ticks, rounded down. The fresh trial checks are evidence for the stated release, not an additional endpoint payment.

## Shared coordinates and masks

The script checks the requested source and receipt SHA-256s, the unchanged historical `weights.py`, `check.py`, `host24.log`, the exact export and edge receipts, and the original gold receipt. The payload authentication inside the pinned host reader checks each full digest; the independent [host audit](quintet-stationary-sol.md) rebuilds all 307 cells and matches the signed host arrays. Both readers check physical snapshot identity through the replay. The numerical recomputation calls neither the joint estimator nor the producer's `lp.lower`.

The ordinary occupation coordinates are `d,k`, with `c=k+1`, plus child receive `u` and mother receive `v`. The exact affine singleton residual is `R=A+B*(u+v)` on finite columns through 24. The same `u+v` appears in the gain coefficient, with its sides separately conserved by shared edge marginals. The checker reconstructs all 1,484 exported rows per side from the original physical `f` charge and the value increments, tests both receive endpoints, and checks both rate coordinates in every finite occupation column. There is no birth-rate multiplier on `R`. All 120 degree-25 occurrence columns get **zero `R`** even though their apparent degree has a positive exported intercept. They are a sentinel for every true degree above 24. The saturated band splits out all child receive through `tail-up`; that auxiliary column has zero objective, while finite receivers continue to see the true marginals even with a tail donor. Host and denominator envelopes still apply to every tail occupation.

The 12-by-12 physical birth kernel was independently rebuilt and matched exactly. Its source mark is the **mother's degree at the tagged newborn's physical birth**. Monotone maternal degree justifies the permanent cap `1/m`, including `m=13` for the last class. Conditional on this mark, rate bounds for the tagged child are `L` and `L+k/2+1/m`; its reward is actual promotion in the selected degree class, counted through 34. A larger discount `0.565775 >= HI_quintet >= mu` lowers all these nonnegative birth values. The physical class relation is `beta >= beta*K`; the column coefficients are `b*(delta_ij-K_ij)` and the nonnegative slack columns enter with `−1`. This is a lower matrix on physical birth fractions, with no residual Perron or normalization of retained mass.

## No shared reserve; signed release comparison

Write `x=sign*(mu−lambda)` with `sign=+1` lower and `−1` upper, and let `M` be the mean original witness potential in ticks. If `H` is the higher-family host cap and `H0` its original global cap, then

```
x*M >= LO*V(J) + g18 + FY_outside - err_per_core + pi[H0-H+R_mask].
```

- **Singleton:** The [separate dual review](quintet-dual-sol.md) proves `pi(F)=mu*V(J)+(HI−mu)*pi(V)+pi(R_full)`. The joint price uses `LO*V(J)` and the masked remainder of that one identity. It does not add the standalone `R` price or a second lifetime floor. The centered `f` row retains the original post-promotion newborn edge; the `V` generator evolves only the old vertex.
- **Edge:** The [separate edge/FY review](quintet-fy-sol.md) prices the full `g` slack only on paying edges with both endpoint degrees `<=18`. Outside FY `(175,5,8)` uses only ordered FY gaps at a paying edge with `max(D_host,D_receiver)>18`, irrespective of donor degree. Thus the two credits are disjoint at the time of payment, including on histories that cross 18. They use neither singleton slack nor `H` or projection-envelope slack. The D18 figure here is the exact floor for the **quintet** witness's degree-18 mask. The strongest FY receipt is pinned; its three exact components were summed again, but its full root optimization was not rerun in this integration audit.
- **Higher and denominator:** `H0−H` replaces the single original global higher cap through its host-local savings. The denominator bounds compare the original potential, including the lower-only recharge with zero stationary mean. Neither denominator bound is a credit on a drift row. The original aggregate error, including the negative upper singleton cap, is deducted once.

The independently repaired lower bounds on `pi[H0−H+R_mask]` are **111136.572367432 / 189366.776697011 ticks**, rounded down. At a positive trial `t`, the objective adds `t*(global_hi−W_hi)` to that same score on each occurrence column, and the anchor subtracts `t*global_hi`. For a negative trial the corresponding term and anchor are `(-t)*(W_lo−global_lo)` and `−t*global_lo`. Both full trials are positive, so their initial net floor is divided by the stationary **upper** denominator, as is each strictly positive joint margin. The signed margins at those intermediate trials are at least `17460.691172550 / 44913.143367107` ticks. Dropping both edge floors retains the same occupation model, denominator and `H/R` score, then recomputes the trial penalty and its dual. It does not transport the full endpoint by simply subtracting the removed credits.

For each objective HiGHS supplies multipliers only. The checker rounds to `10^-12`, evaluates all rational column inequalities, and subtracts the maximum nonnegative deficit for each of the five column classes on the per-core scale. Occurrences and edges have mass one on that scale; excess stubs and tail receive are bounded by mean arity and total receive, and total Palm slack is at most the normalized physical birth mass. The positive birth multiplier uses `LO` once. This gives six independent exact repairs per sign: two denominators, the joint savings, the full and edge-released trial margins, and the direct endpoint trial.

## Evidence

From `[local path]` on one niced core:

```sh
taskset -c 8 nice -n 17 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 round2/.venv/bin/python -Bu round2/checks/quintet-joint-sol.py
```

Exit 0. The shared matrix has 613 rows and 5,554 columns (2,420 occupation, 3,120 edge, 12 Palm slack, one stub and one tail-up). Per sign the replay checked 307 host cells, 1,484 affine rows, 120 unpaid tail occurrence columns, and independently repaired all six objectives. Output includes `rank signs PASS; exact repairs 6` on each sign, both exact receipt comparisons, and the signed release margins tabulated above. Only these two audit files were added.
