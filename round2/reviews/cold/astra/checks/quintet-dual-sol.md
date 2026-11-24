# Cold review: quintet singleton reduced-reward export

**Holds for both signs.** With physical total receive `ρ=u+v`, `πF ≥ LO·V(J)+πR_mask` is an exact inequality. The identity is `πF=μV(J)+(HI−μ)πV+πR_full`, not equality with the masked reward or with LO in the first term. The degree-24 standalone LP gives total singleton floors **391079.994162541737385 / 731984.671964200118015 ticks** (lower/upper, rounded down). This reviews the new singleton export and its sided consumer, not any joint price or new exponent endpoint.

## Which certificate is inherited

The original witness is `b04cc50c1e391e0420eb7b006ac781d484c68543` and the pinned receipt SHA-256 is `d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb`. I inherited the accepted signed full singleton caps and the [original lifetime Bellman proof](quintet-lifetime-sol.md), including the all-degree suffix and arity exit; I did not rerun the production verifier. The new `export.py --replay` actually authenticates both production blobs, reconstructs the **full** frozen `D≤12000,c≤30` sweep with the same analytic continuation, checks its integer inequalities, and reproduces both JSON files byte for byte. The independent audit also authenticates the blobs and reads the raw census with the earlier independent binary oracle. Its new exact residual and grid-quotient checks concern the exported slice, not a new proof based only on low rows that the global residual is nonnegative.

The witness sign is `+1` lower and `−1` upper, applied to the original drift at that witness's own nominal `λ`; add `cap_f=0/−1053` **after** signing. `F` is the resulting nonnegative singleton slack, not the unclipped original signed term or a new comparison denominator. The production singleton promotion term already includes its assigned original newborn singleton and **postpromotion** newborn edge. By contrast `A_phys V` evolves only the old focal core; the virgin corrector is paid once by the stable-age source. The value grid is `10^12` subticks per tick, and production ticks use `2^40` per unit of potential. The LP's `10^7` is temporary conditioning; it does not change either unit.

## Physical rows, seams and overflow

For a nonfounder `c=k+1`, including the mother, `L=D−c`, `J=(2,1)`. If `F=a+B_fρ` including its cap, the exported reduced row is exactly

```
R_mask(D,c,ρ) = a + L·(V_gain−V) + (L/D)·(V_bud−V) − HI·V
                 + ρ·[B_f+(V_gain−V)]
```

on listed cells. Physical gains occur at `L+u+v`; mother and all actual core children send at their own inverse degrees. The production centering switch affects `F`, never the gain clock. Promotions change `c` at fixed `D`, use the actual postpromotion census, and have zero rate at bare states. At degree 64 the row uses the real degree-65 successor. At `c=30` a promotion pays its full rate but exits to zero `V`; above arity 30 the Bellman value vanishes irreversibly, so the full residual there is `F≥0`, not an assertion of zero physical pressure. Above the export's degree cutoff the global `V` and residual can be nonzero. Beyond the original tabular cutoff a downward-rounded seam enters the separately proved `h_c/D` continuation. The normal generator argument is legal: bounded `V`, finite tabular degrees, bounded gain drift for `h_c/D` under rate `≤D`, and promotion rate `≤1`. Discounted Dynkin under the physical stable-age law yields `π A_phys V=μ π V−μ V(J)`; founder omission does not drop the physical send into a retained child. Therefore the displayed bound uses `R_full≥0`, `0≤R_mask≤R_full`, `V≥0`, and `LO≤μ≤HI`.

At LP limit `ℓ`, a finite occurrence column through `D≤ℓ` retains **all** child pressure `u`, including that of a child whose own degree or arity falls in an overflow band; it also retains maternal pressure `v` regardless of the mother's band. Both are the same marginal coordinates used in the gain equation and score. Exact finite states with `c>30` have zero *cost* from the arity mask but retain those rates. The aggregated degree-overflow cell at `ℓ+1` is unpaid by `sided.py`, even if its index numerically equals a listed export degree. Only for saturated tail arity (`k≥7` in the collapsed infinite-degree band) is **all** child receive placed in `tail-up` rather than the occurrence coordinate `u`. It still enters the total receive and shared child-to-parent reciprocal marginal. Its degree gain and saturated promotion stay in the same class, and placing zero score here is a further nonnegative mask. Scoring the tail surrogate as if it were an exact degree row would be wrong, even adding the aggregate tail flux afterward.

## Normalization, physical lift and the kernel

The LP works in `η=π/μ`, with birth forcing exactly one and total occurrence/edge mass `1/μ`. In physical units `π b=μ`, `πρ=1−μ`, `π k=1`. Use the same discount mixture at every finite and grouped tail state. In the saturated tail write `p=1/D`, `b=1−(k+1)p`; the permissible hull has `0≤b≤1−8p`, and `1−b−p=k p` uses the true product as a single observed marginal. Excess stubs carry `k−7`, while `tail-up` carries actual total child receive. Retained exact edge incidences, including edges to a tail donor, lift into reciprocal-interval corners preserving both directional means; no independence of their pressures is needed. Multiplying a **nonnegative** bound on `ηR_mask` by LO yields the reported stationary increment, once. There is no promotion multiplier in this occurrence objective; `b` is used only where the physical balance or birth-Palm model requires it.

The `class_kernel` is for the **birth-Palm distribution of mother degree at the physical promotion**, which is also her degree just afterward. Its classes are `D=2,…,12` and `D≥13`; the last incoming class uses a permanent maternal reciprocal cap `1/13`. Maternal degree is monotone, so for input class with minimum `m`, the actual later gain rate at the tagged child lies between `L` and `L+k/2+1/m` in the full history. Each entry measures a lower bound on discounted **real** births in its target degree class, with zero continuation beyond degree 34, not a probability at an independently sampled mother state. The rows are substochastic; use lower inequalities `w≥wK`, never equality. To see the stationary relation, count the discounted births from each physical virgin classified by her mother's birth-time degree; each later birth has just one mother, and exponential age/transport gives exactly the corresponding promotion-class law at discount `μ`. The inherited kernel uses `HI_kernel=.565775≥HI_quintet`, so nonnegative discounted birth rewards only increase when the discount is reduced to `μ`. Kernel validity depends on the physical virgin, the unchanged directed clocks, and these class/discount conditions, not the witness sign or its nominal trial. Each of the 144 entries was independently rebuilt by checking the conditional one-jump quotients.

The exact dual uses one conservative mass bound `1/LO` for each column group: occurrences and edges by `η(1)`, stubs by `ηk`, tail receive by `ηρ`, and the **sum** of Palm slacks by the total birth mass `η b=1` (kernel entries nonnegative). The checker captures HiGHS's proposal, reconstructs the objective from the JSON and raw `u+v`, and independently repairs the exact sparse reduced costs per group. It also checks that finite gain columns consume the same two pressures, and that the unpaid saturated tail routes its child receive through the separate marginal. The exact extra stationary rewards are:

```
lower 1276783006042809110549818991696409322453507199781/13644757069600662706896800000000000000000000
upper 212108873684137139278951921799899102453819851509/1240432460872787518808800000000000000000000
```

A separate joint LP may include this same masked objective on its physical occurrence lift; the standalone price is not an additional independent payment on top of a joint evaluation of the same rows.

## Bounded replay

From the project root, these commands all exited zero on one niced CPU:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 135 round2/.work/venv/bin/python -Bu round2/occupation-quintet/lifetime/dual/export.py --replay
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 75 round2/.work/venv/bin/python -Bu round2/occupation-quintet/lifetime/dual/sided.py --limit 24
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 17 timeout 100 round2/.work/venv/bin/python -Bu round2/checks/quintet-dual-sol.py
```

Both exports reproduce 1,514 values, 1,484 rows, and 2,968 exact corner checks per side. The independent checker compares the 1,484 Fraction grid quotients and 4,452 residual values (two corners and one interior), and 22 directed-arrow fixtures per side, including nine unpriced seams/overflow states. It detects dropping child receive in 11 fixtures and mother receive in 13 per sign, and verifies two real degree-65 successors. Fixtures include original centering/arity changes 26/27, tabular/tail seams 179/180 and 1999/2000, and 12000/12001 for the original charge. The LP has 613 rows and 5,554 columns; both native and independently repaired replays give extra floors `93573.157772619572338 / 170995.906971746003733` and the full singleton floors above. No degree-64 LP price was run.
