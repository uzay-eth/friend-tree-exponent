# Directed edge and singleton reserves

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 9: Claude and HY4 Flash on specified components; GPT for full replay; inputs include the Round 1 physical exponent and prior, signed raw witnesses and disjoint occupation credits. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


This header reads the two frozen nominal snapshots. It checks the digest of the opened file and returns integer **floors** in units of $`2^{-40}`$. Include it in a C++20 program, link GMP with `-lgmpxx -lgmp`, and have GNU `sha256sum` on the path. Construct one `motif_reserve::Scores` per single-threaded solver process. Its pair and singleton caches each contain 32,768 fixed slots. A collision recomputes an exact result. [Callable header](reserve.hpp)

```cpp
#include "reserve.hpp"
using namespace motif_reserve;

Scores score(path_to_snapshot);
Core parent{24, 2}, child{2, 1};
auto vertex_credit = score.singleton_floor(parent);
auto edge_credit = score.edge_unused_floor(parent, child);
auto known = score.edge_donor_floor(parent, child, Slot::ParentOfParent, Donor::exact(7));
auto pooled = score.edge_donor_floor(parent, child, Slot::ParentOfParent, Donor::pool(7));
```

`Core{D,c}` means ordinary degree and **total core degree including the parent**. Hence $`k=c-1`$ and $`L=D-c`$. The pair order is chronological parent then child. A parent in this pair must have $`c\ge2`$. Both pair endpoints have known $`2\le D\le72`$. `singleton_floor` supports every valid signed-32-bit input degree, including the analytic tail. `Donor::pool(m)` supplies an interval $`[0,1/m]`$ with $`m\ge2`$. A nonexistent or already centered slot contributes zero. The external slot is one of the following.

| Slot | Contact receiving a leaf | External donor | Available copies |
|---|---|---|---:|
| `ParentOfParent` | Parent $`p`$ | Parent of $`p`$ | 1 |
| `OtherChildOfParent` | Parent $`p`$ | Child of $`p`$ distinct from $`w`$ | $`c_p-2`$ |
| `ChildOfChild` | Child $`w`$ | Child of $`w`$ | $`c_w-1`$ |

On a chronological triple $`a\to b\to d`$, book the donor $`a`$ on edge $`(b,d)`$ with `ParentOfParent` and the donor $`d`$ on edge $`(a,b)`$ with `ChildOfChild`. On a cherry with parent $`a`$ and children $`b,d`$, book $`d`$ on $`(a,b)`$ and $`b`$ on $`(a,d)`$ with `OtherChildOfParent`. Each deletion of an undirected leaf gives one incidence. A pooled leaf may act as donor if both remaining edge members are known. Book the intrinsic edge reward once on the pair and the singleton reward once on its vertex.

## Exact pair row

For this discussion every table and jump is signed by $`\sigma`$ and measured in ticks. Fix the physical states of $`p,w`$, put $`G=g(p,w)`$ and let $`\delta_x`$ be the gain jump of this physical edge at $`x`$. If $`D_x<72`$ its successor is an exact table corner, including degree 72 itself. If $`D_x=72`$, its old mass is the tail knot and its gain jump is $`(G_\infty-G_T)/73`$ with the other endpoint fixed. Promotion keeps the true degree and changes $`k`$ by one. These identities also hold when the old physical $`k`$ exceeds the clamp. Births are evaluated at the exact post-promotion states, projecting the surviving outer endpoint to level three for the corresponding cherry or vertical triple.

Start the row with $`-\lambda G`$. For each endpoint $`x`$ include $`(L_x+i_{\rm other})\delta_x`$ and the true $`L_x/D_x`$ times its old-edge promotion difference plus its allocated triple newborn. An inner parent supplies $`(i_w-r_c(p))\Delta f_p`$ and an inner child supplies $`(i_p-r_p(w))\Delta f_w`$. The incoming FY values are $`A_{cp}(p,i_w)`$ when the parent lacks vertical support, $`(k_p-1)A_{cc}(p,i_w)`$ when it lacks cherry support, and $`k_w A_{pc}(w,i_p)`$ when the child lacks vertical support. Only applicable inner states and positive counts contribute.

For each remaining external slot at a contact $`x`$, its contribution for a donor of reciprocal $`i`$ is a function $`\psi(i)`$.

| Slot treatment | $`\psi(i)`$ |
|---|---|
| Active triple absorbs the centered difference | $`z\delta_x`$, using the original $`x_p,x_{cc},x_c`$ center |
| Inner contact with FY remainder | $`i\delta_x-A(i)`$, using respectively $`A_{pc},A_{cc},A_{cp}`$ in the three slot roles |
| Degree-72 contact | $`i\delta_x`$ |

For each slot set $`U=\sup\psi`$. The **physical** optimization uses $`i\in\{0,1/2,1/3,\ldots\}`$. In the frozen grid the only infeasible knots are $`5/12,11/24`$, both strictly between the consecutive possible reciprocals $`1/3,1/2`$. In every other affine interval the two feasible ends bound the function throughout. The constructor checks this property on the actual grid and the first segment extends to zero. Thus a maximum over feasible knots computes the physical supremum exactly. The **continuous** option uses all knots on $`[0,1/2]`$. Uncentered and centered slots have the same supremum in both options. Summing each $`U`$ with its exact multiplicity gives $`R_{\max}`$ after adding the terms above.

For both endpoints below 72, the continuous row agrees with the frozen exact corner expression before its outward integer rounding. At degree 72 the physical formula retains the finite promotion rate, the actual $`L,k`$ even in overflow, and the known neighbour reciprocal. The corner gain and promotion extrema enclose these choices. If the other endpoint is inner, its gain times the degree-72 neighbour reciprocal uses the actual fixed knot jump, which is at most the archived maximum of the knot and infinity jumps. If both are degree 72, each own-tail gain rate is enclosed separately. Consequently $`R_{\rm physical}\le R_{\rm continuous}\le R_{\rm archived}`$, and every row here is covered by the original directed edge cap.

## Using the returned floors

The edge-row credit is $`\lfloor\max(0,-R_{\max})\rfloor`$. For a known exterior donor the credit is $`\lfloor U-\psi(1/D_y)\rfloor`$. For a pooled donor let $`H`$ be the maximum of $`\psi`$ at zero, $`1/m`$, and every FY knot inside that interval. These points suffice because $`\psi`$ is piecewise affine. In continuous mode the result is $`\lfloor U-H\rfloor`$. In physical mode the result is $`\max(0,\lfloor U-H\rfloor)`$, since every actual donor also satisfies $`\psi(i)\le U`$. In particular a pooled donor need not have the reciprocal $`1/m`$. The hot donor calculation evaluates these differences with signed 128-bit numerators and rounds down. Its largest products use at most $`16S`$, $`D_y\le2^{31}-1`$, a knot denominator of 10,080 and jump denominator at most 73; they stay below $`2^{111}`$. Cold singleton and pair-row arithmetic uses GMP rationals.

For an actual edge, its signed physical row plus any chosen donor credits is at most the same $`R_{\max}`$. Adding its row credit gives at most $`\max(0,R_{\max})`$, below the frozen positive edge cap on either witness. Use the **same exterior mode** for row and donors when claiming these together. For a donor-only calculation using the stationary LP's archived worst-donor row, pass `Exterior::Continuous` as the final argument to `edge_donor_floor`. That stronger donor-only credit is already covered by the archived edge cap.

The continuous mode also proves the instantaneous inequality for every hypothetical real reciprocal in a pooled interval. The physical mode proves it for every physical degree compatible with that interval. An occupation row can use this physical source lower bound and independently minimize its Bellman gain over the whole interval because both inequalities hold at the actual tree state.

For a singleton, $`R_{f,\max}`$ uses the fixed centered rates below degree 72, and the signed maximizing neighbor-rate endpoint at all larger degrees. It uses the actual promoted state and analytic singleton whenever required. The returned score is $`\lfloor\mathrm{cap}_f-R_{f,\max}\rfloor`$, where the global singleton cap is 154,457 ticks for upper and 95,343 for lower. The high-degree singleton donor credit already used by the baseline is the further difference between this maximizing term and the physical rate. Thus both singleton credits may be retained together.

**Combining with the stationary LP.** Its existing vertex and directed-edge marginals already claim unused singleton and pair rows. Count a particular unused row only once, or use separate complete alternatives. The original LP leaves edge-donor credits available. Its experimental FY extension can cover vertical and sibling incidences at inner receiving degrees 17–71. Combine those two donor calculations only on disjoint triples/slots or by taking the stronger of complete lower bounds. The original singleton-donor and quartet-donor returns concern different gain rows.

## Reproduction and measured scales

From the repository root, the following command compiles the callable header and validation driver with undefined-behavior checks, runs the exact fixtures on both witnesses, reproduces the held-out totals, and measures cached and new queries on one core. The defaults read `mc-upper.bin` and `mc-lower.bin` from `round2/outputs`. Use `--snapshots DIR` for different input paths and `--work DIR` for a new or empty output directory. Replay

```sh
python3 -B round2/motif-cert/reserve/replay.py
```

Per witness the independent Fraction checker compared **76,216 integer floors**, comprising 238 singleton states and 567 oriented pairs with exact and pooled donors in both optimization modes. It checked 21,310 physical donor incidences and 198 identities combining the singleton row with its full physical donor credit. The fixtures include degree 72, both sides of the analytic seam, clamped core counts and degree $`2^{31}-1`$. The checks detect using a pool endpoint as a fixed degree and paying a singleton or edge row twice. Fraction checker · Verified output

Pass `--all-pairs` to include the exhaustive finite census. For each snapshot it checks all 979,299 table singleton ceilings and all 6,714,612 ordered physical endpoint pairs, including 365,437 pairs containing degree 72. It compares the exact physical and continuous ceilings with the separately evaluated frozen corner at every pair. The maximum physical ceiling is **520 ticks upper** and **10,646 lower**, within the archived caps of 524 and 10,651. The physical unused-row floor is positive on 3,753,890 upper pairs and 5,869,714 lower pairs. The full two-witness replay with this option ran in 2 min 19 s and peaked at 1,207,936 KiB RSS. `edge_residual_ceiling` exposes the signed ceiling used by this audit. Census · Both exhaustive results

Twelve independently seeded friend trees of 32,768 vertices supplied 3,336 nonfounder cores, 2,289 eligible edges and 6,056 eligible donor incidences. The following values are observed sums of the integer floors divided by $`2^{40}`$ and by the number of nonfounder cores. No lifetime or discount is applied in these measurements. Pilot implementation · Upper totals · Lower totals

| Source per core | Lower witness | Upper witness |
|---|---:|---:|
| Singleton unused row | $`6.04433\times10^{-7}`$ | $`1.11566\times10^{-6}`$ |
| Edge unused row, physical exterior | $`3.26825\times10^{-8}`$ | $`1.79961\times10^{-8}`$ |
| Edge donors, physical exterior | $`1.37450\times10^{-5}`$ | $`1.44287\times10^{-5}`$ |
| Edge donors, archived continuous exterior | $`1.37467\times10^{-5}`$ | $`1.44412\times10^{-5}`$ |
| Physical edge donors, degrees $`\ge7`$ pooled | $`9.41203\times10^{-6}`$ | $`1.02261\times10^{-5}`$ |
| Same pool, receiver $`D\le48,c\le12`$, partner $`D\le8,c\le9`$ | $`5.54002\times10^{-6}`$ | $`7.95411\times10^{-6}`$ |

At the upper snapshot, one core averaged 21 ns per cached FY donor, 28 ns per cached pooled FY donor, 4.4 ns per cached pair row, and 4.1 ns per cached singleton. Across 16,384 distinct pairs, new donor queries averaged 0.47 µs and new pair-row queries averaged 7.1 µs. The pair row's GMP calculation is deferred until its first request. The cached donor path uses the precomputed maximizing knot on each prefix and includes the actual pool endpoint. Benchmark · Timing and checksums · Line query driver
