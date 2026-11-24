# Additional physical occupation sources

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 9: Claude and HY4 Flash on specified components; GPT for full replay; inputs include the Round 1 physical exponent and prior, signed raw witnesses and disjoint occupation credits. See the [Astra ledger](../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The read-only C++ source supplies two extra credits in integer units of $`S=2^{40}`$. The first belongs to a real triple or quartet and uses its exact signed row after maximizing the remaining exterior reciprocals. The second belongs to an active triple and one specified outside donor whose quartet extension is outside support. The second credit is booked on their connected four-set. Both can be added to the frozen higher-row comparison, alongside the previously retained four-core and high-degree singleton donor credits. [API](credit.hpp) · [Frozen transfer and donor proof](../higher-rows.md)

## Exact row

Use the signed convention $`w_M=\sigma S W_M`$, with $`\sigma=+1`$ for the upper witness and $`-1`$ for the lower. A real member has degree $`D`$, total core degree $`c`$, $`k=c-1`$ core children and $`L=D-c`$ leaves. The source evaluates old, gained and promoted motif weights at their physical $`(D,c)`$ values. At a tail endpoint it takes the tensor average using its actual $`d_j/D`$. A member outside the exact internal support gives zero weight. A gain that crosses a support boundary therefore has the correct exit jump.

Fix a supported embedding $`M`$ of size three or four. Write $`d_x=w_M(x^+)-w_M`$, $`b_x=L_x/D_x`$ and $`p_y=1/D_y`$. For every undirected leaf $`y`$ of $`M`$, delete it to obtain $`J`$, and let $`x`$ be its unique contact. The incoming term is exactly
```math
 (p_y-z_{J,x,y})(w_J(x^+)-w_J).
```
The calculation reads the underlying edge or triple table, applies the original transfer center and computes the gain on the smaller physical embedding. The center is at the exact pre-gain state. A donor above the contact uses the parent center. Child-center and sibling-center choices for edge-to-triple transfers follow the frozen cut table.

For each member $`x`$, partition its outside slots into the centered slots $`B_x`$ and uncovered slots $`U_x`$. A slot of a triple is centered exactly when adjoining a parent or child on that side gives an active quartet on the **pre-event** states. The new endpoint's state plays no role in these internal support gates. At a quartet every exterior slot is uncovered. Define
```math
 \begin{aligned}
 R^*(M)={}&-\lambda w_M+
 \sum_{\text{deleted leaves }y}(p_y-z_{J,x,y})\Delta_x w_J\\
 &+\sum_{x\in M}\Big[
    (L_x+\sum_{y\in M:\,y\sim x} p_y+\sum_{s\in B_x} z_s)d_x
    +\tfrac12|U_x|\max(d_x,0)\\
 &\hspace{34mm}+ b_x\big(w_M(x^{\rm bud})-w_M+w_{M+xN}\big)\Big] .
 \end{aligned}
```
When $`L_x=0`$ the final term is zero. For a triple, $`w_{M+xN}`$ is its four-core extension after the promoter changes state. At a quartet it is zero. The six preimage/birth cases are listed in the frozen proof. This expression maximizes the signed row over an independent input in $`[0,1/2]`$ at each uncovered slot. It works for physically specified tail members too, since the actual gain difference is an exact rational even there.

With $`J(d,p)=\max(d,0)/2-pd`$, the signed, physically transferred residual in table units is exactly
```math
 R_{\rm actual}(M)=R^*(M)-\sum_{x\in M}\sum_{y\in U_x}J(d_x,p_y).
```
Here $`y`$ ranges over actual donors. A covered slot already transfers $`(p_y-z_s)d_x`$ into a quartet and retains $`z_s d_x`$ on the triple.

## Compatibility with the frozen cap

Fix the physical member states and numerically replace every uncovered exterior reciprocal by an arbitrary element of $`[0,1/2]`$. In the frozen corner proof the only information required about these exterior values is precisely this interval. In particular the total gain coefficient of a tail member still satisfies the two tail-rate bounds and its known internal-neighbor lower bound. All internal reciprocals, gain differences, birth terms and transfer centers remain fixed. Hence the frozen product interpolation majorizes the row for every such choice, including a choice attaining $`R^*(M)`$. If $`c_M`$ is the nonnegative charge from the frozen positive corner parts allocated to this occurrence, then $`\max(0,R^*(M))\le c_M`$. The per-core counting argument already pays these charges.

The new intrinsic credit is
```math
 I_M=\lfloor\max(0,-R^*(M))\rfloor
     =\max(0,-\lceil R^*(M)\rceil).
```
For any selection of distinct uncovered donors one obtains
```math
 R_{\rm actual}(M)+I_M+
  \sum_{(x,y)\ \rm selected}\lfloor J(d_x,p_y)\rfloor
 \le \max(0,R^*(M))\le c_M .
```
Thus the negative part of the tightened row and the losses from known donors can both be retained. The established quartet donor credit uses its own outside slots on a five-set. The additional triple donor uses an uncovered slot of a triple on a four-set. The high-degree singleton credit concerns the scalar vertex potential. These summands have distinct physical origins.

For a donor pooled at a proven degree lower bound $`m\ge2`$, use $`p\in[0,1/m]`$. The floored lower credit is $`\lfloor d(1/2-1/m)\rfloor`$ for $`d\ge0`$, and zero for $`d<0`$. The underlying triple must still consist of real members. For a selected set with a pooled member the intrinsic credit is zero. Deletion of one pooled leaf of a four-set can still give a real triple and a positive donor credit.

## Chronological interface

The header reads the frozen snapshot and takes an explicit rational target. Each connected set has an apex at input position zero. Every other entry's parent index precedes it. The header matches the supplied order to the frozen rooted shape, including the two distinct roles of the fork. The following is a complete call inside an existing chronological loop. [Implementation](credit.hpp)

```cpp
#include "round2/motif-high/occupation/credit.hpp"

mh_occ::Source source("round2/outputs/mc-upper.bin",56581165,100000000);
// Construct the shape once for this chronological pattern.
int parents[4] = {-1,0,0,2};
auto shape = mh_occ::Shape::from_parents(4,parents);
mh_occ::Members members = {{{4,3},{7,1},{3,2},{4,1}}}; // (degree, total core degree)
int64_t units = source.additional(shape,members.data());
```

Compile with C++17, use the project root as an include directory, and link with `-lgmpxx -lgmp`. A pool is supplied as `{proven_minimum_degree,0}`. Ordinary vertices must have $`2\le D<\mathrm{INT\_MAX}`$, $`1+\#\text{selected children}\le c\le D`$. Each call uses a read-only snapshot on a single thread.

- Call `additional` once at the row for each chronological triple or quartet. A triple receives its intrinsic credit. A quartet receives its intrinsic credit plus the donor losses of each eligible triple obtained by deleting one undirected leaf. Include geometric four-sets outside the quartet table support, since their triples can have uncovered donors.
- `row` and `triple_donors_on_four` expose those summands separately. Each triple/donor incidence uniquely specifies its four-set and the donor to delete.
- `prepare` gives the signed `residual_ceiling`, the payable `intrinsic` and reusable gains. Its `triple_donor(contact,Direction::Parent|Child,degree,pooled)` pays for **one actual exterior donor** in an uncovered slot. The direct `Source::triple_donor` computes only the needed gain. A covered or absent slot, or a quartet queried through this triple-only interface, pays zero. The caller identifies each physical donor once.
- `exact_row` returns the full GMP rational $`R^*`$ for audits. All source-returning calls use integer snapshot units. For the lower snapshot the nominal target is `5656453/10000000`.

## Arithmetic

When all old members are strictly inside their own interaction cutoff, self weights and smaller-table source gains use a single table entry even if a successor reaches the knot. A newborn quartet may interpolate at up to two old endpoints. The fast evaluator forms this newborn numerator exactly and uses integer terms for all other transfers. The frozen 4096-subunit outward filter determines the exact integer row ceiling. An ambiguous ceiling is recomputed with GMP. The fast branch has $`D<128`$, $`|w|\le2^{52}`$, $`|z_{\rm num}|,S_r\le2^{40}`$. A centered numerator is bounded by $`127\cdot 2^{40}\cdot 2^{53}<2^{100}`$. A newborn numerator before its rate is bounded by $`127^2 2^{52}`$. Multiplication by the rate numerator and the 4096 subscale fits signed 128 bits, as does the sum of fewer than 64 such terms. All fast denominators fit signed 64 bits. The direct donor floor uses signed 128-bit products.

As soon as a physical member reaches its own level's cutoff the evaluator uses GMP for the entire row. Tensor products of up to three actual degrees, and the change to an endpoint's denominator after a gain, are computed with arbitrary-precision integers. Exact fractions are also retained for the donor gain. The API checks bounds at each final signed 64-bit conversion. [Source and arithmetic](credit.hpp) · [Frozen integer filter](../rows.hpp)

## Numerical evidence

The finite sweep visits all exact states compatible with a physical rooted embedding, at triple degrees through 16 and quartet degrees through 9, including equal-state siblings. The two snapshots each have 6,235,196 active rows in these windows. The largest intrinsic credit over both shapes of triple and all four shapes of quartet is 133 units for the upper witness and 3,923 units for the lower. After assigning chains and forks to their deepest vertices and sharing cherry, Y and trident credits across the apex's children, the **entire capped intrinsic credit** is at most 403/S per core upper and 7,184/S per core lower. These are pointwise ceilings on the size of the available intrinsic source in these windows. The calculation uses the per-apex maxima of the intrinsic credits in the same counting argument as the frozen residual.

In the same triple window a single uncovered donor of degree 2 through 16 can pay 8,400,900,612 units upper or 6,394,805,380 units lower. Each family log reports its counts, largest credits, and a physical degree/core-degree witness. Upper sweep · Lower sweep · All-GMP upper · All-GMP lower

A separate Python program begins with actual physical rates on decorated trees, subtracts each centered transfer by walking actual outside neighbors, adds incoming raw gains and books the newborns. It then adds the actual losses over the remaining unknown slots to recover $`R^*`$. The same physical descriptor must give the same $`R^*`$ across different exteriors. It compares the exact rational and both integer credits with the C++ result, compares $`R^*`$ with the frozen product corner bound, and checks both implementations of donor-only queries. It also checks the chronological sum on active and inactive four-sets, including a pooled donor and all choices of deleted member. The test includes degree 2,000,000,003, seam degrees and concurrent tail endpoints. Doubling one donor credit produces a planted failure. Independent audit · Upper results · Lower results

Eight new physical random friend trees were grown to 200,000 vertices each, using seeds 2026112401 through 2026112408. A uniform vertex and then a uniform neighbor determine each attachment. Each reported value is the exact sum of the integer credits over the eligible noninitial supports, divided by $`S`$ and by the combined 7,582 cores. The founder and donor slots at the founder were excluded from the rewards. The Python audit independently reproduced every family subtotal from a held-out 1,600-vertex tree. These occupation measurements describe the sampled reward sizes. The chronological reserve supplies the discounted lower bound needed for an exponent conclusion.

| Member degree window for triples and quartets | Witness | Intrinsic per core | Uncovered triple donor per core |
|---|---|---:|---:|
| 16, 9 | Upper | $`1.5666\,10^{-12}`$ | $`9.77991\,10^{-6}`$ |
| 16, 9 | Lower | $`4.0803\,10^{-11}`$ | $`1.96929\,10^{-5}`$ |
| 40, 20 | Upper | $`2.95913\,10^{-6}`$ | $`1.48743\,10^{-5}`$ |
| 40, 20 | Lower | $`3.73920\,10^{-6}`$ | $`2.61367\,10^{-5}`$ |

The raw integer numerators and the range across seeds are retained with the summaries. Held-out measurement · Occupancy totals

From the repository root, run `sh round2/motif-high/occupation/replay.sh` to compile, audit and reproduce the measurements. Adding `--all-gmp` also evaluates the full rational physical formula independently at every active row of both windows. Replay
