> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../../README.md#claim-28) for exact scope and citations.

# Gold2 singleton frontier toward `241/426`

**HOLD: the H3-dependent endpoints and budgets below are withdrawn.** The frozen H3 replay mis-scaled its physical reward slopes. The raw singleton roots, singleton bounds, and birth-tag comparisons are unaffected. See [`repair/`](repair/) for corrected lower-side accounting.

The already cold-audited paired baseline is **`[0.565714675151000, 0.565736253459000]`** (`host/h3-joint/receipt.json`). It is stronger than the three-credit Farey profile. With the original Gold2 witness and its existing H3/FY/edge floors retained, the new marked singleton replacement and uniform coupled price give **`[0.565714767317686, 0.565736152862047]`**. This moves the ends inward by at least `9.2166686e-8` and `1.00596953e-7` in exponent units. `241/426` remains inside.

## Two trials, with units

The improvement goal is `2e-6` **from each cold-audited endpoint**, giving trial points `0.565716675151` and `0.565734253459`. The second pair of trials is at `241/426` itself. All entries below are **ticks per core**, with `S=2^40` ticks per potential unit. Each gap range describes a **stationary disjoint extra-credit ledger with the booked non-singleton floors held fixed**. The left number is a necessary floor on the extra credit for this ledger (even with the entire true singleton and all higher-cap savings); the right number is a sufficient extra credit from the solved lower comparison. Use a strict excess to exclude a trial rational. This is an obstruction to closing the target with this allocation of credits, not a lower bound on the true unmeasured residual or an impossibility claim for a different proof.

| Quantity | Lower witness | Upper witness |
|---|---:|---:|
| Already booked singleton floor (replaced) | 207990.215716345 | 296272.606830282 |
| New actual marked singleton floor | 332177.794045450 | 431889.647083866 |
| All-stock upper ceiling on true singleton | 1518437.847417055 | 1613587.014761099 |
| Fixed other floors: H3 plus eligible `g`/FY | 1147273.073216289 | 532790.029908787 |
| Extra credit for `2e-6`: necessary to sufficient | **554828.667314120 – 2593227.512804899** | **719503.102565319 – 2582949.207487477** |
| Extra credit for `241/426`: necessary to sufficient | **12367553.941862948 – 17579455.781670218** | **7899387.993715880 – 11497268.529874944** |

The actual `2e-6` gap lower bounds are about `0.505e-6/0.654e-6` **potential units/core**, and the currently sufficient credits are about `2.358e-6/2.349e-6` potential units/core. These units are not exponent changes. At the two `2e-6` trials the coupled obstruction is stronger than multiplying a separate stationary denominator lower bound by the target displacement. Its negative margin already includes the full original `H` as a possible extra saving, i.e. it is more generous than the existing cap-savings ledger. The same obstruction limits the inward move from the paired baseline to less than `1.592e-6/1.472e-6`, even with perfect recovery of the original singleton. In the sufficient direction the candidate trial margins are `1241.380878383/1323.153480460` ticks/core; their exact division by the stationary `M_hi` yields the endpoints above.

## Which slots can supply the shortfall?

- **Singleton `f`:** no other singleton payment can be added. The new solve replaces the booked one, adding at least `124187.578329105/135617.040253584` ticks/core to its respective floor. Its upper ceiling covers all degrees and all arities, so no untracked maternal mark, extra donor, or high-stock `f` band can invalidate the obstruction. Finer tags can improve a *floor* within the ceiling.
- **Lower `g`/FY:** freeze the booked D18 cap slot and both booked FY masks. Outside those masks use exact ordered `(host, payer, distinct donor)` atoms at the time of payment. Candidate complementary sectors are high-degree mother donors (distinct from the payer), young payers fed by older side children beyond the small mask, and payer `>34` or donor `>5` outside the older-payer mask and its other gates. A larger old-mask solve is a replacement or a pointwise remainder, not another credit. Ordinary non-FY reciprocal gaps and post-FY edge cap gaps must be separated before adding them to any FY contribution.
- **Upper `g`/FY:** the paired ledger books no edge reward. The same real donor-age/degree tags can route a new upper-side solve without colliding with the booked H3; g-to-H3 transfers at M3 are part of H3, not FY. The upper witness also needs roughly `0.72` million extra ticks/core even at the most generous certified singleton ceiling.
- **Non-H3 higher rows:** physical negative-part payments on original `b5`, `s4`, then `p4/s3/p5` are the distinct remaining higher families. To compose with the booked host savings, keep the original post-transfer cap by proving `P+A<=max(P,0)<=max(K,0)` on the actual occurrence. No H3 solve, cap reconstruction, or release update is attempted here. The existing finite-census diagnostic suggests `b5` and `s4` child-donor corrections and lower FY mother-donor/younger-child-payer sectors as places to look, but its empirical row and box gaps are neither certified credits nor ceilings.

## Denser marks and exact coupling

The focal solve uses the pressure comparison at `D<=700,c<=12`, mother live degree/stock through five, and two real newest donors through degree three. It preserves the original signed singleton reward; physical receive drives both reward and gain. The existing analytic suffix and separate arity-overflow bound cover the omitted continuation. The four standard runs each check `2,800,083` physical states and `17,262` explicit corner cases.

For the new tagged birth constraints let `beta_d` be the law of *actual* maternal degree at a birth. A focal virgin with mark `i` has maternal send in `[0,1/i]` thereafter by monotone degree; for `i>12` use `[0,1/13]`. Bound its discounted outgoing birth flux at each actual degree `d=2..12` with directed triangular lifetime comparisons, using `b=(D-c)/D`. Stationary birth balance `beta K_mu=beta` then implies `sum_i beta_i K^-_{id} <= beta_d <= sum_i beta_i K^+_{id}`. In the sided LP these are physical promotion-flux constraints, not empirical frequencies. Each affine inequality has its own nonnegative slack kind, bounded by `1/LO` for exact dual repair because outgoing discounted birth count is at most `1/LO` and the birth mass is one. The augmented LP has 623 rows and 5,564 columns; floating-point optimization only proposes the dual, and every column is checked rationally with directed repair.

With fixed two-donor windows the tags improve the lower floor by **279.744642552 ticks/core** and have no material effect on the upper floor or ceilings (microtick-scale repair). In a separate matched one-donor test at maternal cutoff eight, retaining degrees 6–8 instead of pooling them improved the lower floor by **933.687408916 ticks/core**. That solve's overall `250276.448880140` ticks/core was weaker than the two-donor cutoff-five floor. These are measured conditional comparisons, not a theorem that larger tags can never help. In this lift a longer birth table is a small lever compared with the old-degree receive and higher/edge residual sectors.

For trial `s`, let `delta = sign*(s-lambda)>0`, `B` be the fixed H3/edge floors, `E` the original error, `M` the original comparison potential, and `H-h` the existing unused cap allowance. Price `B-E+pi F+pi(H-h)-delta*pi M` on the same occupation columns. A sufficient bound uses the marked lower value and each upper host comparison; the obstruction uses the marked **upper** value, the lower host comparison, and all `H`. Omitted high-arity singleton reward is charged once on the upper side of this inequality. For the candidate endpoint, uniform (not trial-conditional) pricing lets a positive margin be divided by `M_hi` once. For the `2e-6` ceiling on a possible supported move, the same algebra bounds the supportable displacement because `M<=M_hi`. No old host trial margin or alternate reduced singleton calculation is added.

## Replay

From `.` use the authenticated inputs and existing SciPy environment. Run each witness and each bound:

```sh
for side in lower upper; do
  for bound in floor ceiling; do
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 220s \
      round2/.work/venv/bin/python -Bu round2/occupation-quintet-gold2/lifetime/pressure/frontier/replay.py \
      --side "$side" --bound "$bound" --mother 5
  done
done
```

The one-donor density control uses `--side lower --bound floor --mother 8 --donors 1`. The original receipt, conditioned host input and stronger paired receipt are SHA-256 pinned. The checker reads the audit's already-booked inputs and cells rather than recomputing either host or H3 work.
