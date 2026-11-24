# Physical occupation of the four-vertex outside-donor reward

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 15: Claude on frozen composition and selected inequalities, HY4 Flash on components, GPT on full production; inputs include Claim 9, the five-pass replacement, and the unsigned K7/144 dual. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


For an ordinary connected quartet with frozen score $`H`$, choose a gain site $`x`$ and an ordinary external core neighbour $`z`$. Put $`a=\sigma[H(x^+)-H(x)]`$ in ticks of scale $`S=2^{40}`$, where $`\sigma=1`$ for the upper table and $`-1`$ for the lower table. The physical reward for this incidence is
```math
  J(a,p_z)=\tfrac12\max(a,0)-p_z a,\qquad p_z=1/D_z .
```
All the calculations here apply to the two frozen `mc` tables. Each signed newborn total replaces the archived Q4 amount once. The independent five-pass freeze gives 36,290,118 lower and 32,703,071 upper ticks. [Five-pass audit](split-audit.md) · [Original quartet calculation](../occupation/inside.md) · [Shape definitions](../motif-high/snapshot.hpp)

## Physical measurements

The first sample consists of 12 independently seeded, unconditioned trees grown from the three-vertex tree to 4,194,304 vertices each. Every step uses a uniform initiator and one of its neighbours chosen uniformly. All rates below divide the aggregate ordinary-donor reward by the total core count, include the founder in that count and use parts per million of $`S`$. Shapes or donors touching the founder have a separate boundary convention. Each quartet here consists of four ordinary cores. An incidence with a founder donor is recorded in a separate column. [Fixed seeds and enumeration](census.cpp) · Per-sign, per-family and per-seed output

| Observed ordinary-donor contribution per core, ppm | Lower score | Upper score |
|---|---:|---:|
| All quartet and donor degrees | 28.579025 | 29.227503 |
| All four quartet degrees at most 6 and core degrees at most 5 | 6.152385 | 6.674012 |
| All five sites satisfy those cutoffs | 4.061003 | 4.279431 |
| Some quartet member has degree at least 14 | 17.212486 | 17.984389 |
| Oldest quartet member has degree at least 14 | 10.924810 | 11.959656 |
| Outside parent donor | 7.854242 | 8.042983 |
| Outside child donor | 20.724783 | 21.184520 |
| Additional founder donor | 0.103506 | 0.182800 |

There were 49,213 cores, 13,903 active ordinary quartets and 64,689 ordinary donor incidences. The checkpoint at 524,288 vertices gave 29.780519 and 29.213237 ppm. On a separate continuous-time tree with 235,709 cores, all ordinary quartet donors contributed 30.704867 and 30.621993 ppm. An independently implemented analytic derivative produced the same eight signed family totals on that tree to at least 15 relative digits. These finite-tree quantities describe the measured reward. The occupation bounds below use the integer Bellman inequalities. Independent cohorts · Shared-tree comparison

On the continuous-time tree, take the largest amount any individual candidate can book for each incidence. The original wide mesh covered at most 19.354 lower and 16.836 upper ppm. The per-deletion meshes in the frozen five passes raised those source ceilings to 25.445 and 23.612 ppm. The remaining gaps were 5.260 and 7.010 ppm. Five-sets whose oldest vertex had degree at least 64 accounted for 3.119 and 4.927 ppm of those gaps. At degree at least 512 the values were still 1.622 and 3.057. Internal states missing from every candidate accounted for only 0.480 and 0.400 ppm. These overlapping diagnostics explain which portion of the measured reward a finer old-endpoint description can reach. Incidence-level maxima by quartet family

On the same tree, the degree-64 high window could supply at most a further 1.232217 lower and 1.518053 upper ppm before birth occupation. At degree 180 the corresponding ceilings were 1.151076 and 1.635794 ppm. Taking the largest available amount per incidence across all three marks gave 1.562043 and 2.117281 ppm. These ceilings use the original spent amount at the true state. A cell that contains more than one possible previous amount has to subtract an upper bound when funding a deterministic new score. Marked-source ceiling on the same physical tree


## Physical generator and integer rows

Order each connected set by promotion time. Adding its youngest member at a promotion identifies a unique older prefix. A member with $`k`$ selected children has physical core degree $`c\ge1+k`$. Its total degree is $`D=L+c`$, its promotion rate is $`b=L/D`$, and its gain rate is $`L+\sum 1/D_y`$ over core neighbours. A new child arrives in state $`(D,c)=(2,1)`$, simultaneously with the promoter's change to $`(D,c+1)`$. Partial states for the old prefix use a lower bound on every compatible value of that extension. [Integer comparisons](flow.hpp)

Each chronological position has its own exact $`(D,c)`$ windows. Above its degree window a site can retain its exact core degree up to a threshold, then a lower bound on its core degree. With exact $`c`$ and lower bound $`D\ge m`$, the promotion rate lies in $`[\max(0,1-c/m),1]`$. At a clamped core count a promotion leaves the label fixed and omitted extensions have nonnegative weight. Each internal neighbour on a pooled degree interval supplies a reciprocal in $`[0,1/m]`$. For an exterior slot the comparison uses $`[0,1/2]`$. This permits arbitrary predictable dependence on the rest of the tree.

For an old set $`P`$, write $`A_P`$ for the generator on its own members and $`E_i`$ for the new extension at a promotion at $`i`$. At discount $`u=1133/2000`$, the integer rows give
```math
   (u-A_P)U_P\le F_P+\sum_{i\in P} b_i E_i .
```
A floating point proposal for the nonnegative weight is checked against every rate interval by evaluating its worst endpoint on the proposed signed increment. Each row multiplies by $`2000\cdot 2\prod D_j`$ and checks the residual in signed 128-bit arithmetic. Rows outside a candidate's support have weight zero, with nonnegative continuation and source. Extending, summing over ordinary sets and compensating births gives $`\sum F_P\ge(u-\mathcal L)T+vB`$, where $`v`$ is the singleton newborn value and $`T`$ is the summed auxiliary weight. The expected number of connected sets of any fixed size is $`O(e^{\mu t})`$ by the bounded-order core-degree moment estimate. Every score here is bounded per set. Normalized Laplace integration at $`s\downarrow\mu`$ gives $`[\sum F_P]\ge\mu v`$ since $`u\ge\mu`$. Contributions within fixed core distance of the founder vanish under this normalization. [Localization and growth](../occupation/physical-test.md)

## Physical tail gradients and five passes

All four quartet families have finite internal support and two or three endpoint coordinates. At $`D\ge14`$ a physical endpoint has tail weight $`14/D`$ on its clamped child-count corner and the remaining weight on its infinity corner. Multiple endpoints interpolate by the product of their weights. For a pooled, unchanged endpoint, its corner is shared by the quartet before and after the specified gain. We enumerate all compatible inner states, all possible clamped tail classes, and the largest reciprocal and zero in each tail class. The gain change is affine in each such endpoint separately. Taking the signed value closest to zero over these corners, then rounding towards zero, gives a conservative signed slope. A pooled contact or interior can be stopped. The production split enumerates multiple unresolved endpoints when there are at most 4,096 corner combinations. A negative slope with a pooled donor can pay zero and a positive slope $`a`$ pays at least $`\lfloor a(1/2-1/m)\rfloor`$. [Gradient enclosure](quartet.hpp) · [Exact dedicated endpoint tests](../outputs/five-next-split-slopes.log)

Every undirected leaf deletion of a five-set gives one original quartet and its outside neighbour. The complete decomposition includes deleting the oldest vertex. The primary score of a five-set is
```math
 U_5(x)=\max\{U_{\mathrm{common}}(x),U_{\mathrm{age}}(x),\ \sum_e U_e(x)\},
```
where each $`U_e`$ uses only deletion $`e`$. The two whole-source candidates and the 66 specialized candidates across 24 chronological shapes have independent integer Bellman rows. For any rates in the admissible intervals, a maximizing candidate has the same present value and each of its successors is bounded above by the successor of the maximum. The sum of the deletion candidates receives the sum of their disjoint sources. These facts establish the same robust comparison for the physical maximum.

For each later pass use only the outside gain slack of the immediately preceding auxiliary weights with one through four sites. If $`d`$ is a change in the preceding physical weight, the physical gain exceeds its worst exterior contribution by $`J(-d,1/D_z)`$. Each supported enlarged set has size two through five. With separate auxiliary functions for successive passes, adding the inequalities cancels the earlier auxiliary residuals exactly once. [All mesh cutoffs and pass construction](reserve.cpp) · [Frozen integers](split-reference.tsv)

## Sixth site in the primary auxiliary gain

The previous five passes use sizes one through four of the preceding weights. The primary size-five weight has its own outside-gain residual
```math
   J_6=\sum_{|P|=5}\sum_{x\in P}\sum_{z\sim x,\,z\notin P,\,z\ne o}
       J\bigl(U_5(P)-U_5(P^{x+}),1/D_z\bigr).
```
A six-set leaf deletion recovers the old five-set and contact. The program projects the exact *combined* primary weight, including its sum and maximum, onto all 17,776,584 required exact five-member configurations through $`D=7,\ c=4`$. These cover the pre-gain and post-gain states of every admitted contact. For the sixth-site correction a pooled member of the surviving five-set drops that incidence. On the other rows the resulting signed change is exact. The six-set degree windows in chronological order are 6, 5, 5, 5, 4, 4. A separate sequence of nonnegative scores transports this reward down to its newborn value. It is an additive improvement to the frozen replacement. [Projection](projection.hpp) · Six-site rows · Rational six-site and reversed-sign controls

The six-site auxiliary itself has exterior gain slacks at sizes one through four. Three further passes apply the same signed-gain identity to the immediately preceding auxiliary, summing their separate newborn values with the first six-site value. These passes are funded inside the primary six-site residual. A direct rational sum on 3,573 trees includes a tree made by identifying the founders of two measured fixtures. In that tree the sixth-site gain and the unused original high-degree source are both positive. The check verifies each intermediate source inequality and the resulting physical drift after adding the components. Six-site additional passes · Combined physical check · Combined output

## Remaining original reward at a high-degree oldest member

For each deletion $`e`$ and true five-set state $`x`$, denote by $`C_e(x)`$ the maximum *booked source* across the two whole-source windows and the specialized window for that deletion. A window omitted at the physical state has source zero. Every original primary candidate is funded by $`\sum_e C_e(x)`$, leaving at least $`\sum_e(J_e(x)-C_e(x))`$ as an independent original source. Each auxiliary residual used by the previous calculations is already supported by the first quantity.

Fix an oldest-member threshold $`m`$. For sizes three through five the new descriptors require $`D_{\rm oldest}\ge m`$, with its core degree identified up to ten. The two-member prefix can follow an exact oldest degree below that mark while a younger child evolves. A new five-member cell gives a lower bound $`j_e`$ on the physical $`J_e`$. Project the cell onto every matching previous window to bound its booked source from above by $`c_e`$. If the previous slope is undetermined, this deletion may contribute zero. If the new slope is known negative and the previous donor is pooled, the old payment is identically zero. Otherwise the old donor's known interval and signed slope give an upper bound. The remaining whole-state source takes $`\sum_e\max(0,j_e-c_e)`$. The source lies below the unspent original amount for every compatible realization. Source subtraction and high-degree windows

Several choices of $`m`$ can compete when the third member of a set is born. Each threshold is funded by the same remaining source. At sizes three through five use the pointwise maximum of their scores. At size two use the largest compatible entry in the two available chronological third-member shapes. This preserves the robust Bellman inequality by the pointwise-maximum argument and gives one new newborn value. The final high-degree contribution is the largest of the combined value and the separately funded alternatives. On 3,572 exact trees the verifier measures the original signed physical reward, reconstructs every old matching charge with rational arithmetic and checks both the new source and the full auxiliary physical generator. High-degree verification and double-payment controls
