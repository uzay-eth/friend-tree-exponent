# Cold audit: additive hV donor atoms

**Correct in scope.** The six hV payments add **`10361983749/10^16` lower** and **`14352363253/10^16` upper** to the stationary numerators per core. Their strict margins over `10^−6` are `361983749/10^16` and `4352363253/10^16`. This threshold is on **numerator credit**, not on the movement of each exponent endpoint. The outward intervals are:

| Exact starting ledger | With this hV installment |
|---|---|
| `ead947b` post-EDGE memory | `[0.565678819949264, 0.565769419622021]` |
| `e911a1c` post-EDGE memory plus P4/S3/Yb | **`[0.565679256596412, 0.565769102233118]`** |

## D18 per-census corner classifier

For a live oriented hV at `(A,B,C)`, the exterior port multiplicities are `1,k_A−1,k_B−1,k_C` (parent, apex child, middle child, tip child). An absent port contributes nothing; a routed port gives the entire exterior deviation to its applicable larger motif (respectively P4, Yb, Yt, P4). The remaining coefficient is the signed **gain difference** of the hV weight at the receiving tag. This is not the special tail transfer coefficient of `f/g`. If both exterior `A` ports are free, their signs agree, but they are different physical arrows. Middle arity above 17 kills hV; other tail arities use the stored clamp, preserving whether the tag has a parent. The `inf` middle control has zero hV weight.

The integer classifier takes the difference of two rational Bernstein evaluations by cross-multiplying their positive denominators. Below 18 the control is the single inner row. At 18 the cubic is concentrated at `T3`; at a gain from 17 to 18 its key is distinct from the inner row. Above 18 the four numerators are `((D−18)^3,54(D−18)^2,972(D−18),5832)`, each divided by `D³`. Contraction with the two other tags' controls changes no sign after clearing their denominators. This proves the algorithm for **each legal integer census**, without a degree cap; it does not claim one sign works on an unbounded band. The official price still uses finite tag caps.

The referee independently contracted the original binary tables using rational Bernstein weights, including inner routes, the 17/18 seam, mixed and all-tail axes, the middle arity-17 seam, clamping, and degrees above `10^60`. All **252** censuses per witness matched both the classifier and the original rational slope. Exact zero, when it occurs on a free port, is classified by the equality test on the integer numerator.

## Physical ownership, nonvirgin adversary, and founder

Let `j=k_A+k_B+k_C−2`; the physical complete fringe rooted at `A`, excluding its mother, has at least `Σ(D_i−k_i)+2j` vertices and `3+j` cores. Every counted side branch lies below a distinct immediate exterior core child. For the tracked branch with exact `(d,k)`, the donor, its own ordinary leaves, and its `k` children with at least two vertices each contribute at least `d+k` vertices and `1+k` cores. This upgrades only that donor's previous two-vertex/one-core minimum. At overflow keep the unrefined bounds. A mother is never inside this fringe. The exact complement of the frozen size/core predicate is `S>20 or (S>12 and C≥5)`, and it is monotone in both counts; using lower counts can omit payable rows but cannot introduce overlapping ones. For `D_A≥21` even the apex and its direct children already contain `D_A` fringe vertices. An exactly-three-core fringe of size 21 passes, while one of size 20 does not. The low and high apex bands split **current** degree at 20/21.

For a free physical sender `y`, pay `c p_y−min(0,c/2)`, once for its actual incidence. The all-rate baseline of the signed hV row is nonnegative and unbooked here. The parent, apex side, middle side and tip side atoms partition disjoint arrows of that row even if two slots share the same coefficient. A physical donor may participate in different hV occurrences; summing over occurrences is the original charge convention. The labelled referee reconstructs the raw hV from directed arrow gains, promotions with the correct newborn four-core weight, and the two centered incoming `g→hV` transfers; it compares with baseline plus individual terms. In particular the hV *weight* born on an earlier promotion and assigned to an old `g` is not booked again. The tests use two unequal **nonvirgin** donors at each downward slot, with live core children of their own, and include both large-apex and mixed-tail paths. They passed 36 signed raw rows, 144 fringe-mask comparisons, and 84 controls where omitting an actual donor changed the score.

The promotion event can continue the old registration, attach a side donor, or extend the path. The additional registrations use the same physical rate and are separate successors, never one newborn in both labels. The six legal orders for a fixed side donor (3/2/1 by slot) produce exactly one completed marked path. A known donor appears on its own receive edge with the same reciprocal in the continuation and score. The unknown side remainder has its own box but no extra score. Exact corner enumeration agreed with the integer extremizer on 128 signed rows. The degree and arity caps bound the number of retained registrations, and discarded tags are monotone; finite-horizon discounted comparisons pass to all ages by the inherited stable-age formula.

Ignore founder apices and apices with founder mothers. Founder send has rate one and so produces at most linearly many such direct births in mean, negligible against the exponential mean core normalization; paid paths and incidences per apex and their rewards are bounded on each cap. For a nonfounder mother of degree two, `k=0` is the only census with positive promotion rate, `1/2`; the other census has no ordinary child. A virgin exits irreversibly at hazard at least `1+1/2=3/2` from the leaf's send and its own promotion, even before including parent sends. Thus the fraction of physical births from degree-two mothers is at most `(1/2)/(μ+3/2)` by the physical stable-age formula and `πb=μ`. The parent part takes only the complementary birth class; its mother's degree stays `≥3`, so `p≤1/3` even on short prefixes. The common lower weight is exactly `782836785161/1032836785161`, using `LO=.565673570322`. Multiplying the conditional lower ticks by this weight and rounding down gives both archived unconditional parent rows. The parent class weighting is not applied to child atoms.

## Composition and receipts

The raw hV occurrence is untouched by the post-EDGE `R_f/R_g`, its mixed `g/hC` reserve, and the separate `e911a1c` P4/S3/Yb installment. The inherited size-12 and three/four-core size-20 hV payments are excluded as above. Use the two full corrected EDGE denominator ceilings `1357493815758092773/10^18` and `273436004214812983/(2·10^17)`, and add the corresponding nonnegative hV floor to each signed stationary numerator exactly once. The potential and its denominator ceilings do not change. From the exact `e911a1c` endpoints, this yields before directed rounding:

```text
6853561875850472647665032629043939443/12115632305640977999025000000000000000
328740990573504279104102687021150749/581051508956477588875000000000000000
```

The exponent gains are about `7.63317197e−7` and `1.04977860e−6`. No overlapping older hV price or contrast is added. Both the original D18 GMP proof and physical stationary/Palm identities are inherited.

From `[local path]`, the following independent command passed with both complete original payloads and the five code pins authenticated:

```sh
OPENBLAS_NUM_THREADS=1 nice -n 15 taskset -c 3 timeout 70s \
  round2/.venv/bin/python -Bu round2/checks/hv-atom-sol.py \
  b3f18_lo.bin b3f18_hi.bin --archive [local path]
```

efd8 ran the **five** external archive producer batches, covering all **six** role prices. The referee authenticated the five receipts, checked exit zero, cap/role/state fields and every lower tick against the producer ledger (`occupation-neon/pressure-memory/hv-lifetime/old/EAD947B.md`), and recomputed both Palm-weighted floors. Archive wall times for apex20, apex48, parent20, parent48 and sides were `30.39`, `130.11`, `17.17`, `33.43` and `25.55` seconds. The companion (`hv-atom-sol.py`) independently checks the local rates, tail classifier and exact recombination; it does not rerun the six production lifetime prices.
