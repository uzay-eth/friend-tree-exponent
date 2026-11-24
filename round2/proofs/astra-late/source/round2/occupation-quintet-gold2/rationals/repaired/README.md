> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical-floor and frontier audits and scoped Claude conditional reduction and arithmetic; original Neon-2 Gold2 witnesses and the four reward floors. See the [release index](../../../../../README.md#claim-28) for exact scope and citations.

# A lower bound above 99/175 using four physical rewards and a finite potential bound

For the exponent of the original edge-seeded random friend tree,
```math
 \frac{10043824822832044427}{17754235186087500000}
 \ \le\ \mu\ \le \frac{127993535256830043}{226242120050000000}\, .
```
In particular $`99/175<\mu<142/251`$. Every reduced rational in the interval has denominator at least 426, and $`241/426`$ is inside. Irrationality and transcendence remain open. The proof uses the original Gold2 certificate, the physical stable-age law and the following four stationary reward floors. A calculation on just 20 states already gives the strict comparison at $`99/175`$. [Original Gold2 verification](../../../../../../../reviews/cold/astra/attacks/neon2-quintet/gold2/README.md) · [Physical age and spatial law](../../../../../../astra/source/round2/fringe/spatial-limit.md) · [Exact checker](check.py)

## Four lower-witness rewards

All entries are payload ticks per physical core, with $`2^{40}`$ ticks per unit of potential. The proof rounds the prices down to whole ticks.

| Original row and condition at payment | Floor used |
|---|---:|
| Singleton at the actual receive. The D430 calculation conditions on the live mother and latest child. | 207,990 |
| Ordered edge, paying neighbour of degree at most five and host of degree greater than 18, on the exact two-core rows outside the three-core transfer mask. | 395,261 |
| Ordered edge, paying neighbour of degree 6–34 and distinct donor of degree at most five, with a paying-edge endpoint above 18 on the same eligible row family. | 52,341 |
| Negative part of the signed mother–middle–child row, using an actual additional child of the mother. Mother degree 19–79, middle degree at most five and young child degree at most five, under the arity caps in the receipt. | 251,635 |

The singleton and the paths are different terms of the original signed drift. The two edge terms use disjoint payer degrees. For the path reward, each additional sibling is registered in its actual order of birth. Let $`P_0`$ be the signed row after placing exterior neighbours at their adverse bounds, and let $`z_i\ge0`$ be the improvements from distinct actual siblings. For either sign of $`P_0`$, $`F(z)=[z-P_0]_+-[-P_0]_+`$ is convex, nonnegative on $`z\ge0`$ and satisfies $`F(0)=0`$. Hence $`\sum_i F(z_i)\le F(\sum_i z_i)`$, which fits the negative part of the original row. Thus simultaneous siblings spend the negative part of the physical row once. Its nonnegative local production cap still pays the positive part. The only path price used here is the sibling component. Its independently reproduced value is $`16363762149897743918379249111/65029541606195200000000`$ ticks per core. [Singleton](../../lifetime/README.md) · [Small-payer FY](../../edges/README.md) · [Larger-payer FY](../../payer/README.md) · [Sibling proof](../../diagnostic/lower-sibling/README.md) · [Independent sibling replay](../../../../../../../reviews/cold/astra/checks/gold2-lower-sibling-sol.md)

The four floors sum to 907,227. The original signed error is 19,217 and the net is $`r=888010>0`$. The original nominal trial is $`\lambda=7071421/12500000`$. The stationary drift argument first yields $`\mu>\lambda`$. In detail, set $`Q(t)=e^{-\mu t}\mathbb E P(T_t)`$ for the lower potential measured in ticks. Original comparability bounds $`Q`$ above and below on late times. Integrate $`Q'+(\mu-\lambda)Q`$ and divide by time. The boundary term vanishes. If $`c_E=\lim e^{-\mu t}\mathbb E C(T_t)>0`$, the four stationary rewards give a lower limit of at least $`r c_E`$ on the right. This forces the sign of $`\mu-\lambda`$ before replacing the potential on the left by a ceiling. For the unbounded donor class the FY reward has uniformly bounded reciprocal-degree variation per host. [Cesàro and locality details](../uniformity.md)

## A finite calculation for the potential ceiling

Let $`M_*=1808603193693`$ be the original global upper potential per core. For each actual core of total degree $`D`$ and $`c`$ core neighbours, the original upper allocation has local bound $`u(D,c)`$. The upper allocation assigns each old edge to its younger endpoint, each higher motif to its listed internal vertices, and uses the original singleton. The two internal vertices each get half for a four-vertex path. The upper bound permits an absent edge. Summing bounds dominates the tree potential. Founders contribute a bounded endpoint error if their incoming edge is treated separately. These local upper bounds are obtained from the original weight extrema. The local column named by $`(D,k)`$ has $`k=c-1`$. [Allocation and reconstruction](../../host/README.md)

Fix $`d_0\in\{6,18,24\}`$ and define $`B(D,c)=M_*-u(D,c)`$ for $`2\le D\le d_0`$, with zero reward elsewhere. All retained $`B`$ are nonnegative. The pointwise tree bound gives
```math
 \limsup_{t\to\infty} Q(t)/c_E\le M_*-\pi B ,
```
since $`B`$ is a bounded function of the focal degree and its core-neighbour count. This only uses physical local convergence for a bounded finite-state statistic.

Write $`L=D-c`$, $`p=L/D`$, and let $`s`$ be the exact upper endpoint of the original Gold2 enclosure. Use $`V=0`$ after the first degree above $`d_0`$. The following recurrence defines a nonnegative integer in descending degree and then descending $`c`$:
```math
 V(D,c)=\left\lfloor \min_{g\in\{L,L+c/2\}}
  \frac{B(D,c)+gV(D+1,c)+pV(D,c+1)}{s+g+p}\right\rfloor .
```
The physical gain rate is $`L+\sum_{z}1/D_z`$ over the $`c`$ core neighbours, each of degree at least two. The promotion rate $`p`$ changes $`c`$ while leaving $`D`$ fixed. For each retained state the defining generator inequality is affine in the actual gain rate at the computed value of $`V`$. Its two ends therefore verify every exterior configuration. At exit the continuation is zero, and $`D`$ cannot return to the retained set. All states and contributions here belong to the potential bound.

Discounted integration along a typical newly born core gives a lifetime lower bound $`V(2,1)`$ conditional on its physical birth tree. Each promotion creates one new core in $`(2,1)`$. The physical stationary per-core age identity has birth normalizer $`\mu`$ and discount $`\mu`$. Replacing them by the original Gold2 lower normalizer $`\ell`$ and upper discount $`s`$ is conservative for $`B\ge0`$. Equivalently, stationary balance for the sum of the bounded $`V`$ over physical cores adds one source $`\mu V(2,1)`$ and gives $`\pi B\ge\mu V(2,1)+(s-\mu)\pi V`$. Consequently $`M=\lceil M_*-\ell V(2,1)\rceil`$ is an upper bound in the integrated inequality and $`\mu\ge\lambda+r/M`$.

| Largest retained degree | States | $`V(2,1)`$ | Upper $`M`$ in ticks/core | Margin $`r-(99/175-\lambda)M`$, ticks/core |
|---:|---:|---:|---:|---:|
| 6 | 20 | 619269501408 | 1458274085072 | $`3681682957/781250`$ |
| 18 | 170 | 680528936614 | 1423618796411 | $`321296970031/12500000`$ |
| 24 | 299 | 686326891945 | 1420338814887 | $`2422917810989/87500000`$ |

The upper endpoint in the theorem is Neon-2’s original one. Since $`175\cdot142-99\cdot251=1`$, a reduced $`a/b`$ strictly between the neighbours satisfies $`b=251(175a-99b)+175(142b-251a)\ge426`$. The exact replay checks that the mediant is strictly inside.

## Earlier receipt and replay

The three-credit claim at `6b36cfd` used 699,435.94 ticks/core for the compact H3 row. Its Bellman solver divided a physical reward slope by a rate denominator a second time. The associated 32,325-tick margin and the H3-dependent two-sided receipt at `15c2426` are withdrawn. The present computation uses the actual-sibling price under unshifted gain continuations, a certified singleton price and a separate bound on the potential. It uses the 299 finite local upper entries of the authenticated host export as a table of weights. No linear program, dual multiplier, older H3 receipt, upstream reward, seam reward or additional cap rebate enters the calculation. [Rate counterexample and audit](../../../checks/gold2-lower-h3-alt-sol.md)

From the Astra root, the first command uses only pinned files and the Python standard library. The second also reconstructs all 299 upper entries from the SHA-256-authenticated original lower binary.

```sh
python3 -IB round2/occupation-quintet-gold2/rationals/repaired/check.py
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 round2/.work/venv/bin/python -IB round2/occupation-quintet-gold2/rationals/repaired/check.py --rebuild-host
```
