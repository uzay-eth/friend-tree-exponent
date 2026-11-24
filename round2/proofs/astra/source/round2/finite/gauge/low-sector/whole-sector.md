# Four low cores, first-exit costs and a joint Bellman comparison

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 12: HY4 Flash on the binary and earlier nineteen-label results, Claude on identities/fixtures, GPT on later strengthening; inputs include Claims 2, 5 and 9; wider-prior arguments as specified in the proof. See the [Astra ledger](../../../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## 1. Results

Use the binary positive fringe kernel $`K_s`$ with birth marks $`A=1/2`$ for parent degree below 16 and $`\beta=1/16`$ otherwise. In this note a **low** vertex has total degree $`D\le6`$. A **one-high** state has exactly one vertex of degree at least seven. All entry widths below measure the loss in a specified deterministic enclosure. The true kernel has a single value at each test.

At $`s=0.565809`$, retaining the whole four-core low sector gives $`K_s(\beta,\beta)\in[0.307898452027,3.190580848741]`$. It removes $`0.738815362931`$ of previously charged four-core width, of which $`0.337746999`$ reappears at five low cores and $`0.372195357`$ at four cores with one high vertex. The net width reduction is $`0.018710522`$ to nine decimal places. The complete calculation has 4300 labelled multicore rows. An independent backward solve took 3.30 seconds and 76,296 KiB. Its endpoints agree with the forward calculation within $`2\cdot10^{-12}`$. Backward solve · Forward solve

The most effective three-core extension per added row was a high **root** with $`7\le D\le15`$. It removed $`0.0463764`$ of diagonal width per thousand added labelled rows. The corresponding rates were $`0.0372850`$ for completing three low cores, $`0.0163524`$ for one high vertex below the root, and $`0.0006432`$ for completing five low cores. A joint Bellman calculation on the existing four-core boundary through high degree 24 removed at least $`0.057044995064`$ of width, giving $`[0.314726046818,3.140363448468]`$. The true radius signs near $`\mu`$ remain unresolved by these enclosures. Cost and Bellman replay

## 2. Complete enumeration

A planted state is an unordered rooted core tree with an exact leaf count $`L_v`$ at every core. Its total degree is $`1+L_v+\#\{\text{children}\}`$, including the parent slot. A vertex's degree never decreases. Thus every internal child in an all-low tree was born when its parent had degree at most six and carries permanent mark $`A`$. The exterior root can have either mark. The same argument covers the retained multicore extensions with maximum degree 15. In particular every residual copy from an interior multicore row starts with root mark $`A`$; the $`\beta`$ starting mark only occurs before the first copy. Star rays retain their true newborn marks, including births at degrees at least 16.

Let $`N_k`$ be the number of root-mark-erased all-low geometries on $`k`$ cores. With $`m`$ children there are $`7-\max(2,1+m)`$ choices of root degree. The multiset of the children is counted by
```math
 N_k=\sum_{m=0}^{\min(5,k-1)}(7-\max(2,1+m))
 [x^{k-1}y^m]\prod_{i=1}^{k-1}(1-yx^i)^{-N_i}.
```
This gives $`5,25,185,1530,13915,133751`$ through six cores. Each shape is physically reachable. Grow the root to its desired degree while it still has a leaf, promote its required children, and apply the same construction to each new child. All prescribed finite steps have positive probability. The construction stays within the final size, degree and core-count cutoffs. The extra one-high domains with the high vertex restricted to the root or below it have the same property: every intermediate is all low or has its sole high vertex in the final location. [Independent count](results.py)

All runs keep the existing size-at-most-nine domain, both degrees below 16 in the two-core domain, and the indicated additional sectors. These sets overlap. Completing four low cores adds $`2(1530-136)=2788`$ labelled rows. The forward solver stores each interior geometry once and applies the two root denominators to it. It includes both labelled star rays at degrees 2 through 31. At degree 32 the previous analytic first-promotion tail supplies a certified lower and upper terminal value. This also accounts for old child marks when a high-degree star promotes. [Canonical enumeration and solver](sector.py) · Star and boundary formulas

## 3. A certified first-exit account

For a fixed nonnegative newborn test $`g`$, let $`a`$ be the root mark, $`n`$ the full physical size, $`p=1/D_o`$ and $`b=L_o/D_o`$. Every intrinsic rate is divided by $`s+a+n-p`$. The same denominator is used for the two residual copies of a child $`Y`$, with rates $`a_Y-p`$ for $`Y`$ and $`p`$ for $`Y^+`$. The immediate payment is $`b g(\ell(D_o))/(s+a+n-p)`$, where $`\ell(D)=A`$ below 16 and $`\ell(D)=\beta`$ at or above 16. For an excluded successor $`y`$ with normalized rate $`P_{iy}`$ and continuation interval $`[l_y,h_y]`$, enter $`P_{iy}l_y`$ in the base account and $`P_{iy}(h_y-l_y)`$ in its exit class. The upper and lower endpoints of the finite enclosure solve positive linear equations with the same retained kernel $`P`$.

The forward computation starts unit mass at the specified newborn. Within each copy generation intrinsic transitions strictly increase $`n`$, so an increasing-size pass sums every intrinsic path. An interior copy seeds the next generation. At each transfer the lower and upper masses are respectively rounded down and up in units of $`10^{-23}`$. Multiplication by each state's nonnegative payments then bounds every individual account. This uses the same binary degree tables as the earlier comparison.

For the remaining generations the code checks a finite rational supersolution. If $`F_i^+`$ is the total upper payment in row $`i`$, it verifies $`b_i=W_i-\sum_jP_{ij}W_j>0`$ and rounds $`C=\max_i F_i^+/b_i`$ upwards. Consequently the entire upper solution, and each nonnegative width account, is at most $`CW`$. After 27 generations the residual upper charge is $`C\sum_i m_i^+W_i`$ for the unprocessed masses. It is at most $`10^{-16}`$ in the reported runs. This proof uses the actual finite rows and their exterior payments. It uses no prior interval for $`\mu`$. [Forward algorithm](sector.py) · Raw class contributions

## 4. Cost per unit width

All rows concern $`K_{0.565809}(\beta,\beta)`$. Each timed forward replay computes the $`e_\beta`$ column at both newborns, on one core, with 27 generations. Rows exclude the common 60 labelled stars. The two rightmost columns give the width saved per 1000 **additional** rows and the full replay's user CPU seconds per 0.01 of width saved. Wall times include competition for the shared processor. All numerical savings in this table are rounded to the displayed precision.

| Extension | Compared with | Rows | Wall s | User CPU s | Peak KiB | Width | Saved / 1000 | CPU s / 0.01 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Two low | — | 1290 | 1.29 | 0.62 | 25344 | 2.909670 | — | — |
| Three low | Two low | 1512 | 1.02 | 1.00 | 27520 | 2.901393 | .037285 | 1.21 |
| Four low | Three low | 4300 | 2.91 | 2.86 | 55768 | 2.882682 | .006711 | 1.53 |
| Four low, high root at up to three cores, $`D\le15`$ | Four low | 5016 | 4.92 | 2.13 | 60668 | 2.849477 | .046376 | 0.64 |
| Four low, high below root at up to three cores, $`D\le15`$ | Four low | 5556 | 4.86 | 2.39 | 66976 | 2.862144 | .016352 | 1.16 |
| Four low, either high position at up to three cores, $`D\le15`$ | Four low | 6272 | 3.76 | 3.54 | 71464 | 2.828303 | .027576 | 0.65 |
| Five low | Four low | 31846 | 37.89 | 21.53 | 382104 | 2.864964 | .000643 | 12.15 |
| Five low plus one high at up to four cores, $`D\le10`$ | Four low | 42638 | 29.30 | 23.03 | 484572 | 2.798622 | .002193 | 2.74 |

The last case encloses the diagonal in $`[0.319887257123,3.118508771015]`$. Its one-high five-core exits alone contribute at least $`0.417858801635+0.388896019964`$. Added to its lower endpoint they exceed $`1.12664207`$. Thus the larger solve also moves consequential uncertainty to its next boundary. Matched log summary · Largest run

## 5. The destination of the four-core width

At $`s=0.565809`$, the old and new diagonal width contributions are shown to nine decimal places.

| First exit | Three-low model | Four-low model | Change |
|---|---:|---:|---:|
| Four low | .738815363 | 0 | −.738815363 |
| Five low | .249894504 | .587641503 | +.337746999 |
| Four cores, high root | .068399792 | .213164755 | +.144764963 |
| Four cores, high below root | .051628796 | .279059190 | +.227430394 |
| All remaining classes together | | | +.010162485 |

Of the removed four-core charge, about 97.47% returns elsewhere. Completing five low cores next changes the six-low charge from $`0.046719819413`$ to $`0.344572073703`$ and creates a one-high five-core charge of $`0.261880312885`$. Completing the one-high three-core sector through degree 15 changes its own charge from $`1.399606719817`$ to $`0.256430050828`$. The one-high four-core charge increases from $`0.492223944286`$ to $`1.392416721404`$, and the several-high three-core charge reaches $`0.230008636963`$. Three-low exits · Four-low exits · Five-low exits · One-high three-core exits

Three concrete obstructions follow for the four-low model at the upper Sabba endpoint. Normalise any positive test by $`g_\beta=1`$. The old upper boundary payments are monotone in $`g`$, as are the true continuation values if an exact oracle supplies some exits. A subset's unchanged upper payment for $`g`$ therefore dominates its $`e_\beta`$ upper payment. The complete modified upper value dominates the old $`e_\beta`$ lower value plus the unchanged $`e_\beta`$ widths: the geometry is fixed, the retained equations are positive and linear in these payments, and an exact payment at any other exit is at least its old lower payment.

* The existing one-high three-core widths alone exceed $`1.39960671`$, even with exact values at all other exits.
* With exact values at **every** one-high exit and on the star tail, the other low and several-high widths plus the baseline lower endpoint still exceed $`1.01701699`$.
* With exact values at all five-low and all one-high three-core exits, the remaining widths plus the baseline lower still exceed $`1.20333262`$.

These are lower bounds on the proposed upper certificate for every positive test. They identify which classes must share an improved continuation to reach an upper radius sign. [Exact rational comparisons](results.py)

## 6. Joint Bellman bounds with the actual incoming rates

The four-low model has 518 labelled first exits with exactly three cores and one high degree, 236 with high root and 282 with high vertex below it. Their discounted entrance masses are computed with outward rounding. For a three-core state let $`u(y)`$ be a continuation at a successor, and let $`t^{+w}`$ mean an added leaf at $`w`$. Start the unnormalised Bellman numerator with the root reward $`b_o g`$, the true leaf-gain terms $`\sum_v L_v u(t^{+v})`$, the true promotion terms $`\sum_v b_v u(t^{\mathrm{birth}(v)})`$, and the unshifted copies $`\sum_{y\text{ child of }o} a_y u(y)`$. Start the denominator at $`s+a_o+\sum_v(L_v+b_v)`$.

Each directed internal edge $`v\to w`$ adds $`x_{vw} u(t^{+w})`$ to the numerator and $`x_{vw}`$ to the denominator. Each copy at an old root child $`y`$ further adds $`x_{oy}^{\rm copy}[u(y^+)-u(y)]`$ to the numerator. All true variables sent by $`v`$ equal $`p_v=1/D_v`$. For the rectangular bound give each slot its own value in $`[0,p_v]`$. For the tied bound choose a single value in $`[0,p_v]`$ for every old slot with sender $`v`$, including its copy splits. The third calculation uses the actual $`p_v`$. In every case $`a_y\ge p_o`$ keeps all copy rates nonnegative. At the actual rates the denominator equals $`s+a_o+n-1/D_o`$. Thus all three calculations bound the true continuation, and their admissible controls are nested.

After multiplying each slot coefficient by $`p_v`$, its switch ranges from zero to one. To minimise $`(N+\sum_i a_i z_i)/(D+\sum_i b_i z_i)`$, first include negative numerator-only switches, then sort $`a_i/b_i`$ increasingly for the remaining $`b_i>0`$. Include switches while their ratios are less than the current ratio. The minimum occurs at the resulting prefix by the sign of $`a_i-xb_i`$ at the optimal ratio $`x`$. Reverse both inequalities for the maximum. Tying switches adds their numerator and denominator coefficients before this operation. Use lower successor values in the minimum and upper values in the maximum. Intersect with the existing boundary at each state. The solver rounds its lower and upper values outwards in units of $`10^{-12}`$. Sixty-four whole-tree source/denominator checks and 128 exhaustive small-box extremes independently verify the construction. [Bellman implementation](joint_bellman.py) · [Whole-tree check](check.py) · Check output

One can include every gain while the state still has three cores, one high, and its high degree is at most 24. All retained intrinsic transitions increase the degree sum. Promotions go to four-core states; copies have at most two cores. This gives an acyclic recursion with the old continuation at exit. Incoming old child marks stay fixed. A promotion at parent degree at least 16 receives its true $`\beta`$ mark on crossing this boundary. There are 40 possible rooted low configurations per fixed high degree when the root is high, and 70 when it is below. The cap allows 3960 labelled chamber states before taking repeated copies of the control calculation. Four independent runs by root position and exterior mark took 5.30 seconds wall time in total, 5.12 user CPU seconds, and at most 30,764 KiB each.

Multiplying each pointwise improvement by the lower entrance mass and discarding one unit for the old rounding gives the following guaranteed changes in the $`\beta,\beta`$ enclosure. Unvisited generations have nonnegative improvement. The last two columns compare consecutive schemes directly at each state.

| Full three-core chamber | Lower endpoint increase | Upper endpoint decrease | Extra lower increase | Extra upper decrease |
|---|---:|---:|---:|---:|
| Independent slots | .002442863643 | .048561985590 | | |
| Shared rate per sender | .002453487688 | .049118597455 | .000010624042 | .000556611862 |
| Actual reciprocal degrees | .006827594791 | .050217400273 | .004374107103 | .001098802816 |

Sharing improves the width by at least $`0.000567235904`$ beyond the independent controls. Fixing the actual degrees then improves it by at least another $`0.005472909919`$. This separates the gain due to jointly treating parent rates from the larger gain due to continuing beyond the original frontier. Two Bellman steps already earn $`0.024414352867`$ using the actual rates. [Entrance masses](../../../outputs/finite-low-entrances.tsv) · Two-step comparison · Full chamber sums

## 7. Both no-prior radius tests

The complete four-low domain yields the following true radius enclosures. The lower endpoint uses the Perron root of the entrywise lower $`2\times2`$ matrix. The upper endpoint uses a separate directly evaluated test $`g=(1/5,1)`$ and $`\max_i (K_sg)_i/g_i`$. All inputs come from the two stated spectral parameters.

| $`s`$ | Lower | Upper |
|---|---:|---:|
| 0.565615 | 0.69163337 | 3.37971168 |
| 0.565809 | 0.69135153 | 3.37376332 |

The lower and upper comparison methods each have an obstruction independent of the choice of a positive test. Section 5 treats the upper sign. For the lower sign at $`s=0.565615`$, remove the downward rounding from the **lower** controlled value tables and the star boundary, keeping all column-independent direct barriers and analytic moment coefficients fixed. Call the two newborn outputs $`\Phi(g)`$. The operations defining these lower payments preserve order and positive homogeneity. For the degree-64 tail the coefficient multiplying the lower toll is $`[(1+\lambda)T_{\rm lo}-2M_{\rm hi}]/\lambda\ge0`$; the checker verifies the sign for both root labels. Hence the numerical downward lower certificate is bounded by $`\Phi(g)`$ for every $`g\ge0`$.

Replaying the lower tables with **upward** rounding instead gives a rigorous upper bound on $`\Phi`$ at $`h=(2/3,1)`$. It yields
```math
 \Phi(h)\le(0.429069393997,\ 0.882015181463)\le0.882016h .
```
For any positive $`g`$ let $`r=\max_i g_i/h_i`$ and take an index $`i`$ attaining that maximum. Then $`\Phi(g)_i\le r\Phi(h)_i<g_i`$. Thus the present lower continuation cannot certify $`r(K_{0.565615})>1`$ with any positive test. This bound concerns the lower continuation alone. The actual kernel may have larger action. The lower side needs improved lower boundary payments, while the upper side needs width removed from more than one exit class. [Lower-comparison proof replay](lower_model.py) · Rational bounds on its action · Left upper test · Right upper test

From the repository root, `python3 -B round2/finite/gauge/low-sector/results.py` checks the recorded counts, times, outward intervals, comparison obstructions and Bellman gains. Each raw log begins with its exact command, and the timing line reports peak resident memory. Combined verification
