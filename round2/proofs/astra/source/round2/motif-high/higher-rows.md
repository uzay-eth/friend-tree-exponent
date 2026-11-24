# Four-core motif potentials on the noninitial core tree

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 9: Claude and HY4 Flash on specified components; GPT for full replay; inputs include the Round 1 physical exponent and prior, signed raw witnesses and disjoint occupation credits. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


At $`\lambda_U=56581165/100000000`$ the six higher families contribute at most $`627/2^{40}`$ per core to the signed drift violation. At $`\lambda_L=5656453/10000000`$ they contribute at most $`9081/2^{40}`$. These are bounds on the contribution after the transfers described below. They allow separate positive credits for the observed external donors of four-core rows and for the incoming rate of high-degree vertex rows. The lower-order rows, potential comparison, founder terms and discounted occupation are needed to convert these bounds into an exponent enclosure. Upper scan · Lower scan · [Occupation inequality](../occupation/physical-test.md)

## Physical weights and their support

A noninitial core has one parent, $`k`$ core children, $`L`$ leaves and $`D=L+k+1\ge2`$. Its leaf-gain rate is $`r=L+\sum_{y\sim x}1/D_y`$. Promotion has rate $`b=L/D`$, changes its state to $`(L-1,k+1)`$ and attaches a new core child in state $`N=(1,0)`$. A leaf gain adds one to $`L`$. The bulk calculation uses the parent slot at every noninitial core. The founder boundary is handled in the occupation proof. [Core-coordinate jump law](../../proofs/bracket-motifs.md) · [Founder terms](../occupation/physical-test.md)

A motif occurs once per connected set of core vertices with the indicated rooted shape. In the table, arrows give parentage. Listed internal limits are strict in degree and inclusive in child count. The program reads these limits from each snapshot. The two frozen snapshots have the following limits.

| Family | Physical vertices in checker order | Internal support |
|---|---|---|
| $`V`$ | $`a\to u\to b`$ | $`D_u<24,\ 1\le k_u\le12`$ |
| $`C`$ | $`u\to b,c`$ | $`D_u<22,\ 2\le k_u\le10`$ |
| $`P`$ | $`a\to s\to u\to b`$ | $`D_s,D_u<11,\ D_s+D_u<16,\ 1\le k_s,k_u\le7`$ |
| $`F`$ | $`(s,a,u,b)`$ with $`s\to a,u`$ and $`u\to b`$ | $`D_s,D_u<10,\ D_s+D_u<15,\ 2\le k_s\le6,\ 1\le k_u\le6`$ |
| $`Y`$ | $`a\to u`$, $`u\to b,c`$ | $`D_u<10,\ 2\le k_u\le6`$ |
| $`T`$ | $`u\to b,c,d`$ | $`D_u<8,\ 3\le k_u\le5`$ |

Outside internal support the weight is zero, including at a successor outside support. An existing embedding already has the children required by the lower child-count constraints. Both degree and child count are nondecreasing. Crossing an upper limit excludes that old embedding permanently. Deleting any undirected leaf of an active shape leaves an active lower shape, by the nested degree and child-count bounds. At a promotion, the old embedding obtained by deleting the newborn also passes its lower support before promotion. Every retained internal vertex already had the selected old children, and its earlier degree and child count obey the same upper limits. These facts exclude missing entries into support.

Each degree-one endpoint of a shape has a table label at level $`j`$, where $`j`$ is the number of vertices. At the three interaction levels $`(d_2,d_3,d_4)=(72,28,14)`$ and the respective child clamps are $`(K_2,K_3,K_4)=(24,12,8)`$. If $`D<d_j`$ its label is the exact state. Otherwise its law is
```math
 \mu_j=(1-t)\delta_\infty+t\delta_{T_\kappa},\qquad t=d_j/D,\quad \kappa=\min(k,K_j+1).
```
The reciprocal takes corner values $`0,1/d_j`$, with expectation $`1/D`$. Endpoints above a selected child must have $`k\ge1`$, so their parent-label list omits the zero-child inner and knot labels. Child labels include all states. Different physical vertices use independent interpolation coordinates, including children with the same numerical state. Children of a cherry or trident are unordered. The fork's side child and deepest child have different roles.

The fine-to-coarse projection sends infinity to infinity, a fine exact state below the coarse cutoff to itself, and every other state of represented degree $`d`$ and clamped child count $`k`$ to the coarse tail law with mass $`d_{\rm coarse}/d`$ at its clamped knot. Composition sends the actual fine law to the actual coarse law. Above the fine cutoff the two non-infinite probabilities multiply to $`d_{\rm coarse}/D`$, and successive child-count clamps compose. Thus the physical value of every motif is precisely its interpolated integer weight divided by $`S=2^{40}`$. [Frozen layout](../motif/cert4v_format.h) · [Input view and support](snapshot.hpp)

## Transferring the physical generator

For an old embedding $`H`$, let $`W_H`$ be its physical weight. Let $`\Delta_x W_H`$ and $`\nabla_x W_H`$ denote the change under gain and promotion of member $`x`$, respectively, on the old vertex set. Start with $`\sum_{x\in H}(r_x\Delta_x W_H+b_x\nabla_x W_H)`$. For each promotion, every newly created embedding contains the newborn as a childless endpoint. Deleting that endpoint determines its unique old embedding, which receives the birth term $`b_x W_{H'}`$ with the promoter in its **post-promotion** state. For triples the six cases are as follows.

| Old embedding | Promoter | New four-core term |
|---|---|---|
| $`V(a,u,b)`$ | $`a`$ | $`F(a,N,u,b)`$ |
| $`V(a,u,b)`$ | $`u`$ | $`Y(a,u,b,N)`$ |
| $`V(a,u,b)`$ | $`b`$ | $`P(a,u,b,N)`$ |
| $`C(u,b,c)`$ | $`u`$ | $`T(u,b,c,N)`$ |
| $`C(u,b,c)`$ | $`b`$ | $`F(u,c,b,N)`$ |
| $`C(u,b,c)`$ | $`c`$ | $`F(u,b,c,N)`$ |

Each contribution has multiplicity one because the new vertex identifies the old support. At an edge, a promotion at the top gives a cherry with the old child and a promotion at the bottom gives a vertical triple with the old parent. A newborn edge is booked on the promoting vertex. The birth decomposition stops at four cores.

For a smaller embedding $`J`$, an external donor $`y\sim x\in J`$ contributes $`D_y^{-1}\Delta_xW_J`$. Whenever $`H=J\cup\{y\}`$ belongs to the next supported family, split this into $`z\Delta_x W_J`$ kept on $`J`$ and $`(D_y^{-1}-z)\Delta_x W_J`$ assigned to $`H`$. The contact is unique by the tree property. For an edge extended to $`V`$, use $`x_{\rm par}`$ if the donor is the contact's parent and $`x_{\rm ch}`$ if the donor is its child. For an edge extended to $`C`$ use $`x_{cc}`$. For a triple extended to a quartet use $`r_{\rm par}`$ or $`r_{\rm ch}`$ according to the same donor direction. Each center is the integer at the exact pre-gain contact state divided by $`S_r=2^{30}`$. The smaller edge rows account for their remaining Fenchel–Young terms. [Edge transfer identity](../../proofs/bracket-motifs.md) · [Low-row checker](../motif-cert/low.cpp)

The incoming terms of a larger embedding are therefore indexed by **all its undirected leaves**, including its top when this top has one selected child. Delete the donor and apply the smaller weight's gain difference at its contact. For each remaining outside parent or child slot of a triple, test support of the corresponding quartet before the gain. This test depends only on old vertices that the quartet uses internally. It is independent of the new outside donor's state. At an active extension the triple keeps its centered rate for that slot. At an inactive extension it keeps the physical donor input in $`[0,1/2]`$. Every four-core outside slot keeps its physical donor input. The checker constructs all four-core extensions of each triple by adding a parent above the top or a child at any member, then matching parent relations. Subtracting a centered amount at its source and adding it once at its destination preserves the sum of physical drifts. [Generic extension and cut enumeration](rows.hpp)

The frozen affine envelopes bound each gain difference in weight units. They are verified on the lower level at all fine labels or pairs. With mean $`m`$ and radius $`e`$, this says $`\Pi(m-e)\le S\Delta W_J\le\Pi(m+e)`$, using product projection for a pair. An incoming source at a coarse corner uses the donor factor $`z'=\mathop{\mathrm{inv}}\nolimits (y)-z`$ and contributes $`z'\sigma m+|z'|e`$ to the upper bound on the **signed** drift, with $`\sigma=1`$ for upper and $`-1`$ for lower. The donor coordinate is distinct from every coordinate of $`J`$. Expand its reciprocal at its own corners first, then apply the enclosure to the other coordinates for either sign of $`z'`$. For the symmetric smaller cherry, averaging each two-dimensional envelope with its transpose is a valid enclosure because the two fine gain differences agree after swapping arguments. This makes the corner computation invariant under permutations of the retained children even for an asymmetric supplied envelope. [Exact projection checks](../motif-cert/envelopes.cpp)

For reproducibility, the envelope bindings can be read directly from these cuts. The symbols name the payload arrays in the frozen specification.

| Larger family | Deleted leaf | Enclosed gain difference and arguments |
|---|---|---|
| $`V(a,u,b)`$ | $`a`$, $`b`$ | $`\Delta_u g(u,b)`$ via dguw at $`b`$, $`\Delta_u g(a,u)`$ via dgpu at $`a`$ |
| $`C(u,b,c)`$ | either child | dguw at the other child |
| $`P(a,s,u,b)`$ | $`a`$, $`b`$ | p4 first envelope at $`b`$, p4 second at $`a`$ |
| $`F(s,a,u,b)`$ | $`a`$, $`b`$ | f4 first envelope at $`b`$, f4 second at $`a`$ |
| $`Y(a,u,b,c)`$ | $`a`$, a child | y4 hc at the two children, y4 hv at $`a`$ and the other child |
| $`T(u,b,c,d)`$ | any child | y4 hc for the same apex at the other two children |

## Finite corner inequality

Work in integer weight units and let $`w=\sigma S W_H`$ at a corner. Start the row with its incoming envelope contributions minus $`\lambda w`$. For each exact member add the following terms, writing $`d_x=w(x^+)-w`$.

1. Add $`L_x d_x`$ and $`d_x/D_y`$ for each exact adjacent member $`y`$. An infinite adjacent corner contributes zero. A knot adjacent corner uses the quadratic bound below.
2. For each centered external slot add $`z d_x`$. If $`n`$ external slots remain, add $`n\max(d_x,0)/2`$. The outside parent exists only at the top of the embedding. There are $`k_x`$ minus the number of selected children outside child slots.
3. If $`L_x>0`$, add $`(L_x/D_x)\big(w(x^{\rm bud})-w + w_{\rm new}\big)`$. For a triple, evaluate the new quartet by projection of all preserved endpoints. Its active internal vertices project to single exact states. Its newborn is the exact degree-two state. At most two old interpolated endpoints remain, and every product weight is rational. For a quartet $`w_{\rm new}=0`$. At the exact-to-tail boundary a gain lands at the appropriate knot with probability one.

For a physical tail endpoint $`x`$, let $`v`$ be its unique selected neighbor, necessarily exact. At its infinity corner both own transition contributions vanish. At its knot corner the own gain is bounded above by
```math
 \max_{\gamma\in\{g,1\}}\gamma\,[w(x=\infty)-w(x=T_\kappa)],\qquad
 g=\begin{cases}
 \max(0,d_j-\kappa-1+1/D_v)/(d_j+1), &\kappa\le K_j,\\
 0,&\kappa=K_j+1 .
 \end{cases}
```
Indeed $`t-t^+=t/(D_x+1)`$, so the physical expression has the knot mass $`t`$ times this jump with coefficient $`r_x/(D_x+1)`$. This coefficient is at most one. For fixed $`k_x=\kappa`$ the lower bound $`(D_x-k_x-1+1/D_v)/(D_x+1)`$ increases in $`D_x`$. For $`\kappa\le K_j`$, promotion changes the knot index to $`\kappa+1`$ with coefficient $`b_x`$ inside the same mass $`t`$. Its two coefficient extremes are $`\max(0,d_j-\kappa-1)/d_j`$ and one. Overflow leaves the knot index unchanged. The promoter in an active newborn quartet must be internal and have degree below the fourth-level cutoff.

A tail donor $`y`$ also affects the physical gain rate of its exact neighbor $`x`$. Here the reciprocal and $`d_x`$ involve the *same* interpolation parameter. Let $`d_0,d_1`$ be that gain difference with $`y`$ fixed respectively at infinity and its knot. Since
```math
 (t/d_j)\big((1-t)d_0+t d_1\big)\le(t/d_j)\max(d_0,d_1),
```
put zero at the infinity corner and $`\max(d_0,d_1)/d_j`$ at the knot. At each other coordinate expand in its two nonnegative weights. A maximum of linear functions of those coordinates is bounded by the interpolation of the corner maxima. Thus this reduction applies to each incident tail donor even if two or three differently interpolated donors meet at one vertex. [One-knot and quadratic lemmas](../../proofs/bracket-motifs.md) · Nonzero counterexample to same-corner multiplication

These steps give a finite rational $`R_H(c)`$ at every active corner. For the physical, transferred row $`J_H`$, with product weights $`\omega_c`$,
```math
 \sigma S(J_H-\lambda W_H)\ \le\ \sum_c\omega_c R_H(c).
```
Interior support decisions agree throughout each interpolation cell. If a fine endpoint becomes an internal vertex of a possible extension, its degree must be below the new cutoff. A physically interpolated fine endpoint fails this internal condition at both coarse corners. Exact states above the new cutoff also fail it at both corners. The source enclosure, exact births and the tail bounds above cover the remaining terms.

## Retaining two physical credits

For $`a`$ equal to a signed physical gain difference in integer weight units and $`p=1/D_y`$, define
```math
 q(a,p)=\tfrac12\max(a,0)-pa
        =\max(a,0)(1/2-p)+\max(-a,0)p\ge0 .
```
For a four-core row with four exact states, its external gain contribution at a specified contact is maximized as $`n\max(a,0)/2`$. Subtracting its physical contribution leaves exactly $`\sum_{y\ { m external}}q(a,1/D_y)`$. The occupation calculation takes rows whose needed old and successor weights can be read exactly. For those rows **the same corner upper bound still holds with the retained sum of $`q`$'s added on its left**. The source envelope and other extremal terms retain their original inequality directions. Each selected incidence $`(H,x,y)`$ corresponds uniquely to the connected five-set $`H\cup\{y\}`$ together with its deleted undirected leaf $`y`$. Different leaf deletions of one five-set represent different donor/row incidences. A pooled donor of degree at least $`m`$ has $`p\in[0,1/m]`$. Its exact minimum contribution is $`a(1/2-1/m)`$ for $`a\ge0`$ and zero for $`a<0`$. Downward integer rounding preserves the lower bound.

In a vertex table row with $`72\le D<D_F`$, the low-row checker maximizes the signed gain over total uncentered neighbor input in $`[0,(k+1)/2]`$. If $`a=\sigma S(f(L+1,k)-f(L,k))`$, the difference between this extreme and the physical gain is
```math
 (k+1)\max(a,0)/2-a\sum_{y\sim x}1/D_y=\sum_{y\sim x}q(a,1/D_y).
```
The selected vertex rewards use $`D+1<D_F`$, so both values come from the finite table. Parent and child donor incidences have disjoint ordered slots. A physical edge may contribute in both directions, with a different recipient at each end. These vertex credits come from the drift of $`f`$, while the five-core credits come from the drift of $`h_P,h_F,h_Y,h_T`$. Even if both use the same event, they are distinct additive summands of the linear generator.

Let $`Q_4,Q_1`$ be any sums of these eligible nonnegative contributions, with each incidence used at most once in its own family. Let $`e_{\rm low},e_{\rm high}`$ be the respective signed per-core residual bounds in weight units. Before the founder correction one therefore has
```math
 \sigma S(\mathcal L\Phi-\lambda\Phi) + Q_4+Q_1
 \ \le\ (e_{\rm low}+e_{\rm high})\, C_{\rm core}.
```
Thus positive violations of nominal rows are charged once through their signed caps. The same outward rounding margins hold with the physical credits included. [Four-core and high-degree vertex credits](../occupation/five_reserve.cpp) · [Independent physical credit check](physical_check.py)

## From rows to a charge per core

Write $`u_H=\lceil R_H\rceil`$ in weight units and take a positive part whenever an embedding may be absent. A deepest vertex determines at most one vertical triple and one length-three chain. For the fork it determines its parent and grandparent. The grandparent's other child has at most $`k_s-1`$ choices. Thus a uniform per-deepest-vertex charge for these three families is
```math
 \max(0,\max u_V)+\max(0,\max u_P)
              +\max_{F\ \mathrm{rows}}(k_s-1)\max(0,u_F).
```
For each exact apex state $`s=(L,k)`$, let $`M_C(s),M_Y(s),M_T(s)`$ be the separate maxima of zero and the row ceilings at that apex, with zero for absent support. Its total charge is at most
```math
 A(s)=\binom{k}{2}\big(M_C(s)+M_Y(s)\big)+\binom{k}{3}M_T(s).
```
An exact apex with $`k\ge2`$ can pass $`A(s)/k`$ to each of its children. The number of parent-child edges is at most the total number of cores. Hence $`\max_{s:k\ge2}\lceil A(s)/k\rceil`$ is a valid additional per-core charge. Interpolation at exterior vertices preserves each corresponding supremum.

| Family | Active corners | Upper maximum $`\lceil R\rceil`$ | Lower maximum $`\lceil R\rceil`$ |
|---|---:|---:|---:|
| $`V`$ | 28,329,840 | 167 | 5,647 |
| $`C`$ | 10,398,780 | 125 | 904 |
| $`P`$ | 9,643,480 | 35 | 705 |
| $`F`$ | 5,467,736 | 31 | 985 |
| $`Y`$ | 11,332,200 | 64 | 1,723 |
| $`T`$ | 1,591,659 | 87 | 50 |
| Total row count | 66,763,695 | | |
| Chains and fork, per deepest vertex | | 337 | 7,440 |
| Cherry, Y and trident, per child | | 290 | 1,641 |
| **Higher residual per core in units of $`2^{-40}`$** | | **627** | **9,081** |

The per-anchor profiles give the maxima used in this calculation. A row ceiling of $`1`$ records any strictly positive violation up to one weight unit. Upper profile · Lower profile

## Quantization, arithmetic and independent checks

The two binary files define rational weights exactly. Their generation with a floating proposal incurs no separate approximation allowance, since all inequalities are evaluated again on their stored integers. The loader checks the frozen header, exact file size, little-endian LP64 layout, radii, and arithmetic ranges. For the row formula it permits $`|w|,|m|,e\le2^{52}`$, $`|z_{\rm num}|,S_r\le2^{40}`$, and all finite corner degrees at most 127. A signed difference is at most $`2^{53}`$. A symmetrized source uses a doubled mean or radius, so their combined absolute numerator is at most $`2^{54}`$. Its donor numerator is at most $`128\,2^{40}=2^{47}`$. A projected newborn has at most two old mixture denominators, each at most 127, and its promotion factor adds one such denominator. Every term denominator fits signed 64 bits. Before division the largest scaled numerator is at most $`2^{101}\,4096=2^{113}`$. Each row adds fewer than 64 such terms in signed 128-bit arithmetic. Integer ceiling outputs and fork-charge conversions have range checks. The geometry preflight also bounds index and memory sizes. The emitted sources require the independent affine-envelope pass for the same frozen file. [Arithmetic implementation](rows.hpp) · [File checks](snapshot.hpp) · [Shared envelope verifier](../motif-cert/envelopes.cpp)

For each rational term $`n/d`$ the fast checker adds $`\lceil4096n/d\rceil`$ to an integer $`U`$, and counts the inexact divisions in $`b`$. Thus the exact row lies in $`[(U-b)/4096,U/4096]`$. If the ceilings of these two ends agree, that value is the exact integer ceiling of the row. If they differ, the checker rebuilds the complete row with GMP rationals and takes its exact ceiling. All 66,763,695 active entries are visited once. Symmetric loops include diagonal label choices and the checker compares each count with the frozen payload size. [Scanner and aggregation](check.cpp)

An independent Python reader reconstructs the on-disk offsets and all rooted embeddings. For each decorated tree it differences the physical potential $`\sum(g+h_3+h_4)`$ before and after each possible gain and promotion, explicitly attaching the newborn. Its booked calculation separately walks physical outside neighbors and uses hand-listed incoming deletions. Exact equality is required for each transferred incidence, each newborn support, and the sum of drifts. Each realized higher row, **including its retained four-core credits**, is checked against the product interpolation of independent exact GMP row queries. The test covers changing support, multiple identical child labels, up to three simultaneous tail mixtures, clamped high child counts and degree 1009. It also checks the high-degree singleton gap identity using the stored vertex table. Its inert apex exercises the complete local bulk algebra. The treatment of the actual founder is in the occupation proof. [Independent physical implementation](oracle.py) · [Test](physical_check.py) · Upper results · Lower results

Mutation tests change one entry in each of the six families in a private copy for each nominal exponent and require the corresponding frozen family budget to fail. The test also plants a negative radius, an out-of-range integer, a missing fork transfer, a missing newborn, and duplicate claiming of four-core donor credit. An isolated nonzero fixture shows exactly how direct multiplication at the same tail corner can understate physical drift. The original snapshot digests are checked again after both mutation runs. [Mutation implementation](negative.py) · Outcomes

## Reproduction

From the project root, run the following command for the nominal upper and lower witnesses. It checks the input digests, builds the checker, scans every family, replays physical transfers and tests private mutations. GMP, a C++17 compiler and Python 3 are required. The integration bundle supplies the frozen fixture and the two binary inputs. Every stage runs on one thread and writes its logs under the assigned output prefix. Replay script · Input digests

```sh
sh round2/motif-high/replay.sh
```

For a second pass evaluating the full rational formula with GMP at every corner, use `sh round2/motif-high/replay.sh --all-gmp`. To inspect a chosen row, `mh-check FILE --query --pq P Q` reads the family and its ordered corner labels from standard input. It returns the exact signed residual in weight units for each line.
