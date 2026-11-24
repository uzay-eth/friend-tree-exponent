# Physical stationary occupation of overlapping core motifs

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 15: Claude on frozen composition and selected inequalities, HY4 Flash on components, GPT on full production; inputs include Claim 9, the five-pass replacement, and the unsigned K7/144 dual. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


The primary application is the nominal four-vertex witness in `round2/motif-cert/`. The two-sided calculation uses physical vertex and oriented-edge occupation, signed residuals, a per-state potential ceiling and an exact rational LP dual. The earlier three-vertex witness supplies a separate test of the method on a degree-two tip and one hub. All paths below are relative to the repository root.

## Occupation and the founder

Write $`C`$ for core count, $`B=GC`$ for the promotion rate and $`G`$ for the physical generator. For $`z>\mu`$, integrate expected additive observations against $`e^{-zt}\,dt`$ and divide by the corresponding integral of $`C`$. The latter diverges as $`z\downarrow\mu`$ by the renewal result in `round2/fringe/renewal.md`. Along any convergent subsequence denote a normalized integral by brackets. Then $`[1]=1`$ in per-core notation, $`[B]=\mu`$, and $`[GH]=\mu[H]`$ for each bounded-per-core count tested here. The generator identity follows by integration by parts, with vanishing normalized initial value. All coordinates have finite bounds, including the sum of all outgoing slots per core, which tends to one.

The initial centre has at most one promotion per unit time in expectation. Its expected number of core children is $`O(1+t)`$. Each omitted occurrence at the founder involves a bounded number of selected positions and bounded branching at active internal vertices. Counts of incident founder edges or active motifs therefore disappear in the normalized limit. Degree tests are constant or have bounded generator at arbitrarily high degree. The nominal four-vertex physical comparison also uses the founder correction in `round2/occupation/physical-test.md`.

## Shared physical marginals

An ordinary vertex type is $`(d,c)`$, where $`c`$ is total core degree, $`k=c-1`$ is its number of core children, and $`L=d-c`$. Along the physical process, gain takes $`(d,c)`$ to $`(d+1,c)`$ at rate $`L+\sum_{y\sim x}1/d_y`$. Promotion takes it to $`(d,c+1)`$ at rate $`L/d`$ and creates a tip $`(2,1)`$.

`round2/moments/occupation.py` builds nonnegative vertex counts $`v_s`$ and oriented parent-child counts $`e_{ab}`$. Incoming edges sum to $`v_b`$, and outgoing edges of an exact type sum to $`k_a v_a`$. Above the joint cutoff $`K`$, give each vertex an auxiliary colour with reciprocal value zero or $`1/(K+1)`$, choosing the latter with probability $`(K+1)/d`$. Use independent colours on distinct vertices and the same colour wherever a vertex recurs. Their conditional mean is its physical reciprocal. All dynamic tests of a high endpoint aggregate both colours. This preserves arbitrarily high physical degrees, including their contribution to the exact identity $`B/C=1-\sum_{a\to b}(1/d_a+1/d_b)/C`$.

For each finite middle state $`s`$, real triples $`T_{asb}`$ satisfy $`\sum_aT_{asb}=e_{sb}`$ and $`\sum_bT_{asb}=k_s e_{as}`$. The symmetric cherry incidence of child type $`b`$ is $`(k_s-1)e_{sb}`$, with equal labels counted twice. These constraints supply the actual neighbour contributions to both ends of every projected edge. The LP uses the physical vertex and edge generator with the published interval $`0.565\le\mu\le0.5665`$ applied to each nonnegative count.

A bare connector is $`(2,2)`$. Four-site counts $`Q_{asBw}`$ with bare $`B`$ have left marginal $`T_{asB}`$ and right marginal $`T_{sBw}`$. A path with $`n=1,2`$ consecutive bare connectors and arbitrary recorded endpoints loses an internal connector at total rate $`1/d_a+n-1+1/d_w`$. New paths ending in a tip arise when the old terminal tip promotes at rate $`1/2`$. The external endpoint fluxes in these two path populations also project to the same edge flux. This gives stationary balances even when the top is an arbitrarily large hub. Additional one-sided birth balances cover a branching middle. These are equations of one time-independent joint occupation measure.

`round2/moments/marginal_lift.py` keeps every $`(d,c)`$ and its incoming parent and child reciprocal rates up to a second degree cutoff. Their sums, reciprocal-weighted sums, and outgoing-slot counts match the same high-end edge marginals. The remaining degree is unbounded. `round2/moments/full_channel.py` additionally maps every active high middle of the original vertical table to its actual high–bare–donor triple, with constraints from its vertex and parent-edge bounds.

## Nominal four-vertex comparison

The directed integer inequality in `round2/motif-cert/proof.md` has sign $`\sigma=1`$ for upper and $`-1`$ for lower. Let $`S=2^{40}`$, $`m=[\Phi]`$, and $`r_f,r_g`$ denote signed singleton and edge residual caps in unscaled units. Let $`h`$ be the verified higher-family cap per core. The physical credits $`J_p,J_c,J_4`$ use distinct parent-to-receiver and child-to-receiver singleton incidences and quartet-to-fifth-site incidences. In physical occupation,

```math
 \sigma(\lambda-\mu)m
 \ \ge\ [J_p+J_c+J_4]-[r_f+r_g]-h.
```

`round2/moments/v5_profile.cpp` calls the same nominal integer `Low::vertex` and `Low::edge` functions as the upstream checker. Each of the 979,299 finite singleton rows is evaluated with its actual degree and child count. Each of the 6,483,402 paired fine edge rows is visited before aggregating by endpoint class. The nominal exponent is an explicit argument. The profile uses the signed extrema, including negative caps. For a high degree the edge kernel interpolates the class corner and infinity. The paired cap takes the upper bound over every corner supported by that endpoint group, including all eligible clamped and infinite corners.

Singleton drift uses the neighbour-extremized rate inside $`r_f`$. In a table-tail state its physical donor credit is

```math
 J(a,i)=\tfrac12\max(a,0)-ai,\qquad
 a=\sigma\{f(L+1,k)-f(L,k)\}.
```

The LP bounds each directional sum using the same parent and child input marginals. A pooled state uses a lower bound over the signed jump interval. For each direction the physical credit also dominates the existing newborn lifetime reserve times $`\mu`$. The occupation variable names the donor. The archived label names the receiver, so “parent vertex slack” supplies the child-donor variable and “child vertex slack” supplies the parent-donor variable. The LP imposes both lower bounds on the same credit variable. The quartet credit uses the separate fifth-site lifetime reserve. An independently proved edge-donor reserve can be added if the incidences were left uncharged by the nominal edge cap. The signed singleton and paired edge row improvements in this LP already use the corresponding F-ROW and E-ROW sources.

| Witness | Global singleton cap | Global paired edge cap | Higher cap | Parent donor reserve | Child donor reserve | Quartet reserve |
|---|---:|---:|---:|---:|---:|---:|
| Lower | 95,343 | 10,651 | 9,081 | 301,604 | 11,253,928 | 7,402,100 |
| Upper | 154,457 | 524 | 627 | 3,226,994 | 8,148,713 | 7,709,975 |

The caps are per core and the reserves are per newborn, each in ticks of $`2^{-40}`$.

The table `mc-mbar-*.i64` bounds the unsigned potential assigned to each real vertex. Beyond degree 1399 the singleton is bounded by maximizing its analytic expression on the leaf-fraction bin. The high-degree bottom profile is affine in reciprocal degree. Four integer units cover the four individual upper roundings. For each oriented edge the LP subtracts a certified lower bound on the bottom edge maximum minus its true unsigned edge weight. This improves the denominator without using a drift reward.

Tail cells share their actual vertex mass, child-slot count, and incoming reciprocals with the fine vertex marginals. Cells group total degree, core children and physical leaf fraction $`L/d`$. In the strengthened model their integrated birth rates also match. For $`T\ge0`$ and child threshold $`j`$, stationary tests of $`(1-T/d)_+\mathbf1_{k\ge j}`$ use monotonicity to enclose the generator on each cell. A gain at $`d\ge T`$ changes the first factor by $`T/[d(d+1)]`$. A promotion crosses the child threshold only at $`k=j-1`$. The remaining analytic singleton cap uses an exact Bernstein bound separately on each of 16 leaf-fraction panels. The underlying input hashes and component reserves are embedded in the dual.

If $`D`$ is the occupied upper ceiling for $`m`$ and $`N`$ the net lower credit above, a rational proof of $`N-qD\ge0`$ with $`q>0`$ gives $`\sigma(\lambda-\mu)\ge q`$. The approved physical comparison supplies positivity of $`m`$. The LP optimizes the numerator and denominator on one measure.

## Additional V5 overlap credits

For an inner receiving state outside the appropriate vertical or cherry support, the existing edge checker uses the conjugate of its signed FY table $`A`$. For a retained edge jump $`d`$ and an outside donor of reciprocal $`i`$, the available amount is

```math
A^*(d)-id+A(i), \qquad A^*(d)=\max_j\{\xi_jd-A(\xi_j)\}.
```

`round2/moments/v5_donor.cpp` evaluates a lower bound for each such incidence with receiving degree 17–71. In each group it encloses the signed jump by its two extremes over compatible opposite corners. For each donor reciprocal, a chosen knot yields the lower bound $`\min_{d}(\xi_j-i)d+A(i)-A(\xi_j)`$ on that interval. Take the best knot and zero, then minimize over all possible donor labels and middle states. A covered slot contributes zero. At a tail donor $`A`$ is affine on the first knot interval, and the two endpoint tests cover the physical reciprocal. Using separate corner breakpoints for the other endpoint also covers its interpolated jump. The integer divisions round down.

`round2/moments/v5_overlap.py` places these prices on actual vertical triples and unordered sibling pairs in their common middle group. The vertical price covers the parent donor to the descending edge and the child donor to the ascending edge. Each sibling price includes both orientations. The variables project to the same directed edge counts. For groups with a variable child count the marginals are bracketed by its minimum and maximum multiplicities. Where the fine vertex lift covers the group, its total pair count and reciprocal-weighted incidences supply additional constraints. The FY extension overlaps an independently computed E-DONOR reserve on the same incidences.

`round2/moments/v5_score.cpp` subtracts further comparison defects. For a vertical triple with known inner states it bounds its true vertical term and the contributions of forks and chains assigned to its bottom. For a four-site chain with a bare second middle it can further replace the maximum over a free top by the observed top. At a finite cherry apex it bounds cherries, stemmed cherries and three-child stars conditional on two child classes. Each star has three constituent unordered pairs, so a pair receives one third of its bounded defect. Corner minima and maxima respect all three interpolation scales. These deductions concern the unsigned comparison score. The signed edge and donor budgets remain independent of these comparison deductions.

The obstruction uses a lower bound for the actual unsigned potential. `round2/moments/v5_floor.cpp` takes the lower end of the archived bottom, cherry, stem and star comparisons and retains the actual singleton at every finite degree. On the analytic tail the quadratic singleton takes its minimum at a panel endpoint, with the reciprocal correction enclosed separately. The lower score can bound how large a *fixed* numerator divided by the true potential could possibly be. The same rational dual check handles both directions.

### Standalone nominal potential per core and per newborn

Let $`M`$ denote the Abel occupation of the unsigned nominal potential per core and $`B`$ the corresponding birth rate. The mode `--score-only` has two objectives, $`M`$ and $`M/[B]`$, with $`[B]=\mu`$. Its inequalities are the stationary occupation balances, overlap marginals, reciprocal relations and tail bounds on the original published interval. Its inputs are the nominal snapshot, upper comparison profile for that same snapshot, and unsigned score reductions. The shared profile also stores signed drift columns. Those columns have zero coefficients in the standalone denominator objectives and constraints. The solver-free check of this model can run independently of every lifetime-reserve file.

`round2/moments/v5_edge_score.cpp` tightens the unsigned comparison on three families assigned to a bottom edge and three assigned to its parent. Fix an edge with bottom $`w`$ and middle $`m`$ of degree above the joint cutoff $`K`$. For the vertical triple, four-vertex chain and fork, maximize the actual coefficients over the remaining compatible sites with $`m`$ and the class of $`w`$ fixed. A fork has $`k_{\rm top}-1`$ possible side children. Each bottom edge owns its full vertical, chain and fork contribution. For the cherry and stemmed cherry at $`m`$, assign half of each term to each constituent child. For the three-child star assign one third to each child. The resulting conditional upper bounds have multiplicities $`(k_m-1)/2`$, $`(k_m-1)/2`$ and $`(k_m-1)(k_m-2)/6`$ on this outgoing edge. Subtract these upper bounds from the archived upper allowances. The difference is a lower bound on the excess already included in the unsigned comparison.

For $`D_m\le K`$, the original triple and cherry deductions still apply. When the parent of $`m`$ has degree above $`K`$, the original triple deduction used only the vertical component. The extra path row in this case assigns the chain and fork excess conditional on a high parent. At a tail endpoint the physical weight is a convex combination of the clamped corner and infinity. A minimum difference over all applicable corners and all middle states in the observed class is therefore valid. The finite group enumeration reaches degree 73. At that point all active middles have disappeared, all degree bands are represented, and every clamped and infinite outer corner has entered its final band. Integer division rounds each positive deduction down.

In `round2/moments/check_score.py` the independent Python oracle evaluates the physical weights at each actual vertex on deterministic test trees. It checks the edge or path charge against exactly its allocated families, including interpolation on degrees up to 2,000,003. The unsigned price sheets can accompany a previously checked positive numerator because their objectives price only the potential.

The optional `--tail-score` uses the actual birth intensity in the score comparison for a tail cell. Write $`q=L/D`$ and let $`s(D,k)`$ be its unsigned upper comparison in ticks. For each of a small number of integer slopes $`h\ge0`$, `round2/moments/score_tail.py` computes $`a_h=\max\lceil s(D,k)-hq\rceil`$ over every finite state of that cell. For $`D\ge1400`$ it maximizes the upper affine or lower concave-quadratic singleton after subtracting $`hq`$, bounds the affine bottom term at both reciprocal endpoints, and allows four upward roundings. If $`s_0`$ is the constant cap for the cell, $`v`$ its mass and $`b`$ its physical promotion intensity, its excess allowance is at least $`((s_0-a_h)v-hb)/S`$ for every slope. The LP subtracts the maximum of those lower bounds and zero. For the physical occupation this quantity is at most the excess already charged to the cell, and at most $`v\max(0,\max_h(s_0-a_h)/S)`$, the declared variable bound.

`round2/moments/check_tables.py` independently decodes the two snapshots with the Python certificate layout. Its integer checks cover every finite per-state maximum, edge corner and retained unsigned overlap deduction. Each birth-conditioned line is tested by cross-multiplication on every finite state. In each analytic cell the checker tests both reciprocal endpoints and the exact quadratic discriminant when an internal minimum is possible. NumPy is used only for batches of bounded integer operations.

## Edges ending in a degree-two tip

For a state $`s=(D,c)`$ write $`k=c-1`$, $`b_s=(D-c)/D`$ and $`e_s`$ for the count of descending edges ending in a tip $`(2,1)`$. Let $`p_s`$ sum the reciprocal of the hub's parent over these edges. Let $`w_s`$ sum all sibling reciprocals over the same edges. Parent gain carries the edge forward in total degree with flux $`g_s=(D-c+1/2)e_s+p_s+w_s`$. A tip loses its type at rate $`3/2+1/D`$. The exact generator of this edge count is

```math
G e_{D,c}
=g_{D-1,c}-g_{D,c}
 +b_{D,c-1}(v_{D,c-1}+e_{D,c-1})
 -(b_{D,c}+3/2+1/D)e_{D,c}.
```

Impossible preceding states supply zero. Incoming gains at the joint cutoff are the actual path and sibling flux from the smaller edge system. When a hub has $`n`$ tips, its outgoing tip edges receive at least $`n(n-1)/2`$ in total sibling reciprocal from the other tips alone. The LP includes the linear supporting bounds of this quadratic. The sums of $`e_s,p_s,w_s`$ in each degree and child-count band agree with the corresponding edge marginals and are bracketed by the projected triple and sibling reciprocals.

For the descending edge the outer tip is fixed. `round2/moments/v5_tip.cpp` evaluates its paired residual and unsigned score exactly up to the checker rounding, and computes affine lower supports $`a+bi`$ for the available FY donor surplus at its parent. The conjugate uses all archived knots. The lower convex hull spans all possible physical reciprocals, including zero as a limit. The open interval $`(1/3,1/2)`$ contains no reciprocal of a core degree. Both coefficients are rounded down in ticks. Each support gives $`a e_s + b p_s`$ for the donor above the hub or $`a(k-1)e_s + b w_s`$ for sibling donors. A group takes the stronger of this bound and the earlier budget on its tip-touching incidences. Other projected donors retain their earlier weights.

The star refinement places one joint distribution on a hub, its parent and the unordered degree classes of all its children. It covers degrees 17 through the last complete degree band up to 71, and up to seven children. Children of hubs with at most three children are classed as tips, other degrees through four, degrees five through sixteen, and degrees at least seventeen. At larger arity the last two classes are pooled. Parent classes distinguish degree two, degrees three and four, degrees five through sixteen, and degrees at least seventeen. These nonnegative configuration masses project to the existing vertex, tip-edge, oriented-edge, path and unordered sibling measures. Parent, child and tip-conditional reciprocal sums are enclosed in the same configurations. `round2/moments/v5_star.cpp` floors both directed vertical contributions and both directed sibling contributions using the real middle state and the extrema permitted by each outer class. The sum over a configuration uses each incidence once. For each middle group the star budget and the entire previous FY budget are alternative lower bounds on one source.

For the optional `--singleton-gap` test let $`F_{\rm cut}`$ be the actual source function used by the archived two-site lifetime proof. That proof gives $`[F_{\rm cut}]\ge\mu R`$ for each receiving orientation. Any pointwise lower bound for $`J-F_{\rm cut}`$ on the same directed pairs can be added to $`\mu R`$. The replay checks the archived degree and core-degree cutoffs. A receiver beyond $`D=1024`$ or $`c=24`$ contributes zero to that source. An unresolved donor beyond $`D=16`$ or $`c=8`$ uses the pool lower-degree label nine. The excess is at least $`a(1/9-i)`$ for positive signed jump $`a`$ and at least $`-ai`$ for negative $`a`$. An excluded receiver retains the entire physical $`J`$. The LP groups receiver states with a common sign and a rounded-down jump magnitude. Pair counts and sums of actual donor reciprocals project both to the receiver's physical input and to the oriented-edge law. Only table-tail receivers $`72\le D<1399`$ receive this additional credit. Each direction takes the stronger of the full stationary credit and the lifetime reserve plus its own uncharged excess.

Optional `--use-prior` first replays the two degree-32 anchor duals with the original growth bracket. Their exact rational endpoints tighten the coefficient of each stationary balance in a successor LP. The checker verifies the acyclic predecessor chain and its digests on every such replay.

The independent objective `original-M-per-birth` reads the original three-vertex potential in `verify/motif-bracket/reference/`. It assigns the finite singleton and its clamped tail to actual vertex cells and the old edge weight to both projected endpoint classes. At a high middle it charges the possible vertical triple to its bottom edge, and assigns half of each cherry to each of its outgoing edges. These amounts use the actual group of the middle and the old support in that group. The smaller joint states retain their direct path and cherry scores. This supplies a potential-per-birth ceiling to a separate calculation that already has a lower numerator for the original witness. The objective contains no donor reward.

## Commands

`round2/moments/PINS.md` gives the standalone source and data inventory, the separate signed coefficient proof and the SHA-256 files for both.

From the repository root, `bash round2/moments/replay.sh check` verifies the final nominal duals and runs the finite physical tests. `bash round2/moments/replay.sh profiles` recompiles the extractors and compares the raw integer rows, FY and score profiles byte for byte with the retained inputs. `bash round2/moments/replay.sh unsigned-inputs` runs the independent integer array and analytic checks. Its interpreter defaults to the project NumPy environment and can be set with `MOMENTS_PYTHON`. `bash round2/moments/replay.sh donor-api` checks every positive star price on real endpoint pairs against the separate GMP reserve API. `bash round2/moments/replay.sh first-witness` reproduces the original degree-two-tip obstruction and the earlier normalization results. The numerical comparison is in `round2/moments/RESULTS.md`.

Use the following command for an isolated, solver-free check of the degree-5 upper denominator. The data access guard permits exactly the price sheet and its six nominal potential inputs under `round2/outputs/`. Use `lower` in the filename for the other witness.

```sh
python3 -B round2/moments/check_unsigned.py round2/outputs/moments-v5-unsigned-upper-dual.json.gz
```

To reproduce both substitutions into the fixed upstream positive numerators, run `bash round2/moments/replay.sh denominator unsigned`. Each isolated denominator check finishes before the separate numerator is read. To substitute a different independently checked complete numerator with the same nominal residual caps, use `python3 -B round2/moments/denominator.py unsigned --newborn-total LOWER UPPER`, replacing the last two arguments by the total integer ticks per newborn in units of $`2^{-40}`$. The script uses fractions through the final outward rounding.

The higher-resolution unsigned proposal can be regenerated from an existing copy of the project NumPy/SciPy environment with one CPU.

```sh
round2/.work/venv/bin/python -B round2/moments/v5.py solve round2/outputs/moments-reproposed-denominator.json.gz \
  --side upper --cutoff 7 --vertex-cutoff 144 --degree-tests --score-only \
  --score-overlap --edge-score --tail-score --rows M-per-core M-per-birth
```

To propose a fresh upper nominal dual with SciPy, run the following command. Both nominal snapshots, the companion comparison profiles and the named prior reserve logs must be present in `round2/outputs/`.

```sh
round2/.work/venv/bin/python -B round2/moments/v5.py solve round2/outputs/moments-v5-reproposal.json.gz \
  --side upper --cutoff 5 --vertex-cutoff 80 --degree-tests --edge-donor --score-overlap \
  --star-children 7 --use-prior --rows shift edge-donor-min
```

## Rational verification

`round2/moments/dual.py` uses SciPy HiGHS to propose prices. For $`Ax=b, Gx\ge h, 0\le x\le U`$, any prices $`p,y`$ with $`y\ge0`$ prove the lower bound

```math
 b^Tp+h^Ty+\sum_j U_j\min(0,c_j-(A^Tp+G^Ty)_j)
```

on $`c^Tx`$. Every coefficient and correction is evaluated exactly. `round2/moments/proofcheck.py` reconstructs the same inequality independently in a common integer unit. Checking uses Python's standard library. All proposals pin one CPU and cap address space at 6 GiB. No Monte Carlo inputs enter any LP.
