# The four-vertex integer witness

> **Author/model:** Astra (GPT-6.2). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21` ([original blob](../../../SOURCES.tsv)).
> **Review and dependencies:** Claim 9: Claude and HY4 Flash on specified components; GPT for full replay; inputs include the Round 1 physical exponent and prior, signed raw witnesses and disjoint occupation credits. See the [Astra ledger](../../../../../LEDGER.md#astra) for the exact scopes and dependencies.


## Potential and generator

Orient the core away from the founder. An ordinary core has state $`s=(L,k)`$, where $`k`$ counts core children, degree $`D=L+k+1`$, reciprocal $`i=1/D`$, leaf-gain rate $`L+\sum_{y\sim u}i_y`$ and promotion rate $`b=L/D`$. Gain changes $`s`$ to $`s^+=(L+1,k)`$. Promotion changes it to $`s^b=(L-1,k+1)`$ and creates a child with state $`J=(1,0)`$. In the ideal model the founder has one inert parent with zero reciprocal and the same state convention as an ordinary core. Apply the physical theorem away from the founder, with the separate founder lemma for the boundary.

At level $`j=2,3,4`$, use the point mass on the exact state for $`D<d_j`$, otherwise use masses $`t=d_j/D`$ at $`T_j(\min(k,K_j+1))`$ and $`1-t`$ at infinity. An inert parent is a point mass at infinity. A motif value is the expectation of the corner table over the **product** of these endpoint masses. They are interpolation weights of a deterministic function on a fixed tree. Each internal vertex must satisfy the support restrictions in the header. A value outside support is zero. The raw layout enumerates all degrees, internal states, unordered pairs and triples as in the frozen header. Quantized values define the potential exactly.

Sum the singleton $`f`$, edge $`g`$, vertical three-path $`V`$, cherry $`C`$, vertical four-path $`P`$, fork $`F`$ (one branch of length two), stem with a cherry $`Y`$, and three-child star $`T`$. Each embedding is an unordered vertex set with its inherited orientation. The dummy can supply a missing ancestral endpoint. Its score and its own transition rates are zero as a singleton. Write $`\Phi`$ for the potential in physical units and $`S=2^{40}`$. The verifier measures every weight and gain in ticks of $`1/S`$, so its total potential is $`S\Phi`$. Orient the arithmetic by $`\sigma=+1`$ for upper and $`-1`$ for lower. Every local test now bounds an upper residual, generator minus $`\lambda`$ times the signed value.

## Accounting for the existing motifs and births

For an existing embedding the old-vertex drift is the sum of its rate times its change at each physical position. For each newly born embedding remove its new leaf. The remaining vertices are the unique old embedding. A singleton creator receives the newborn singleton and the edge from the **post-promotion** creator to the newborn. At an edge $`(p,w)`$, a promotion of $`p`$ creates a cherry using $`w,J`$, and a promotion of $`w`$ creates a vertical triple using $`p,J`$. At any triple, append the newborn child at the promoting position. The ordered roles are fixed by the old edges, and symmetric siblings are unordered. At four vertices there is no higher newborn. A previously existing embedding cannot enter support through a change of an old vertex. It already has the number of selected children required by the shape, all upper limits on $`D,k`$ and any sum of internal degrees are downward closed, and gain and promotion are nondecreasing in these variables. An old cut embedding of an active newborn embedding lies in the required smaller support, because the high internal supports nest. These observations include supports with a lower arity condition reached for the first time through the newborn.

The edge transfers use signed weights and signed Fenchel–Young (FY) knot values. For $`D<d_2`$ replace the singleton gain rate by $`L+r_p+k r_c`$. Transfer $`(i_y-r_p)\Delta f`$ for the parent and $`(i_y-r_c)\Delta f`$ for each child to their incident edge. Each endpoint is then counted exactly once. For $`d_2\le D`$ retain the entire physical gain rate in the singleton row until its extremization.

For an edge $`p\to w`$, the counterpart's reciprocal in the gain rate is known. If $`p`$ is in the vertical-middle support, retain $`x_p(p)`$ from its parent's rate. Outside that support, use the FY transfer below. If $`p`$ is in the cherry support, retain $`x_{cc}(p)`$ for each of its $`k_p-1`$ other children. If $`w`$ is in vertical-middle support, retain $`x_c(w)`$ for each child. In each retained case transfer $`(i_y-x)\Delta g`$ to the triple consisting of the edge and donor. Only an exact inner vertex can have such support. For the remaining donors at an inner contact $`u`$, move a signed piecewise-linear value $`A(i_y)`$ from the edge receiving the gain to the adjacent edge containing the donor. The three FY arrays correspond respectively to a parent donating to a child edge, a child donating to the parent edge, and one child donating to another child edge. Thus the recipient keeps $`i_y d-A(i_y)`$ with $`d=\Delta g`$. The donor edge receives $`A(i_y)`$ per target edge. The counts at $`u`$ are respectively $`k,1,k-1`$ for these incoming terms. In the dual remainders the respective counts are $`1,k,k-1`$.

On each knot segment $`x d-A(x)`$ interpolates its two endpoint values. It is at most $`A^*(d)=\max_j(\xi_j d-A(\xi_j))`$. This function is convex in $`d`$, so the same bound holds after interpolation over the opposite endpoint. The first knot interval is exactly $`[0,1/d_2]`$, making incoming $`A(i_y)`$ its exact corner interpolation even for an unbounded degree. This argument accepts any signed knot values.

An existing triple retains the appropriate $`r_p,r_c`$ contributions from outside donors precisely when adjoining such a donor gives a supported four-vertex motif. It transfers each remaining centered difference to that motif. The source terms of a triple or a quartet can equivalently be found by deleting each undirected endpoint of its shape. The deleted vertex is the donor, the incident vertex the contact, and the change is the lower-family gain at this contact. Uncovered donors of a triple and all outside donors of a quartet at an inner contact have reciprocals in $`[0,1/2]`$. They can be bounded separately. For a tail contact its entire own-vertex rate is handled by the tail identity below. All transfers are paired between adjacent motif sizes or between two incident edges. Summing them over a tree preserves the physical generator.

## Interpolation and corner inequalities

For a projection from fine cutoff $`d_f`$ to coarse cutoff $`d_c`$, infinity maps to infinity. An exact state of fine degree $`D<d_c`$ maps to itself. Every other finite labeled corner of fine degree $`D`$ maps with weights $`d_c/D,1-d_c/D`$ to its coarse clamped class and coarse infinity. The fine tail corner has degree $`d_f`$ in this formula. Since clamping composes and $`(d_f/D)(d_c/d_f)=d_c/D`$, projection carries the physical fine mass exactly to the coarse mass. The identity holds separately in every coordinate.

For each lower-table gain used as a higher source, the checker enumerates its fine corner values $`d(s)`$. It verifies $`|d(s)-(\Pi m)(s)|\le(\Pi r)(s)`$ with $`r\ge0`$. Use $`\Pi\otimes\Pi`$ in the two-coordinate case. For $`2\to3`$ these gains change either end of $`g`$ at each supported middle. For $`3\to4`$, the two path gains change the respective outer endpoint of $`V`$. A fork uses $`V`$ on its descending branch and $`C`$ on its side branch. The $`Y,T`$ sources use the change of the retained middle of $`V`$ or apex of $`C`$. In a cherry the interchange of its two children gives two valid enclosures. Their average is also valid. The integer builder fits midpoint/radius envelopes and the standalone checker substitutes all fine states or fine ordered pairs by clearing positive denominators. The construction of a candidate envelope has no role in acceptance.

Taking expectations proves the same center/radius enclosure at each physical tuple. The multiplier $`i_y-z`$ is on the **deleted donor**, disjoint from the arguments of the lower-table gain. At physical arguments,

```math
(i_y-z)d\ \le\ (i_y-z)\,E m+|i_y-z|\,E r
\ \le\ E[(i_{j}(y)-z)m+|i_{j}(y)-z|r].
```

Here $`r\ge0`$, $`i_y=E i_j(y)`$, and the last inequality is convexity of the absolute value. The table center changes sign in a lower certificate, while the radius stays nonnegative.

For a tail endpoint let $`H_1,H_0`$ denote the class corner and infinity values, with other labels fixed. The gain rate times its change is $`t\gamma(H_0-H_1)`$ with $`\gamma=(L+\sum i_y)/(D+1)\le1`$. If the class is exact $`(k\le K_j)`$ and an adjacent reciprocal has lower bound $`a`$, then

```math
\gamma\ \ge\ \max(0,(d_j-k-1+a)/(d_j+1)).
```

Indeed the lower expression before substituting $`d_j`$ equals $`1-(k+2-a)/(D+1)`$, increasing in $`D`$. At level 2 take $`a=1/D_y`$ only for an exact inner opposite state and zero otherwise. A higher tail endpoint has a known inner neighbour. A clamped overflow class may use zero. Promotion contributes $`t b(H_{1,\mathrm{next}}-H_1)`$ with $`b\in[\max(0,d_j-k-1)/d_j,1]`$ for exact $`k`$. In the overflow class its difference vanishes. Extremize these rates with the signed jump. At infinity these contributions vanish. A dummy has no physical events and its single infinity label also contributes zero.

If an inner contact sees a tail neighbour, the same neighbour appears in its reciprocal and its own interpolated gain jump. With signed jumps $`d_0,d_1`$, the exact product is $`t((1-t)d_0+t d_1)/d_j\le t\max(d_0,d_1)/d_j`$. The latter is an affine corner charge at the class label only. Apply the argument conditional on all other labels. For other unknown or extremized nonnegative rates, a convex combination of corner maxima bounds the corresponding physical product. Together these identities cover mixed inner/tail pairs, several simultaneous tails, exact cutoff degree, and arbitrarily clamped arities.

## Singleton tail independent of the proposal

For $`d_2\le D<D_F`$ check both gain-rate endpoints $`L`$ and $`L+(k+1)/2`$, the exact promoted state at the **old** degree, and the interpolated newborn edge at that state. At $`D=D_F-1`$ use the analytic gain successor. For $`D\ge D_F`$ put $`x=1/D,q=L/D`$ and write the **unsigned** analytic table in integer units as $`w=a+bq+cq^2+eqx`$. The upper coefficients are $`(A,B,0,0)`$, the lower ones $`(0,B_1,-B_2,-C_D)`$. For a gain and a promotion the new coordinates are respectively $`((q+x)/(1+x),x/(1+x))`$ and $`(q-x,x)`$. The signed residual is bounded using either endpoint $`r=\rho/x`$, $`\rho\in\{q,(1+q)/2\}`$. The newborn $`f(J)+g(s^b,J)`$ is $`C+d_kx`$, where $`C=f(J)+g(\infty,J)`$ and $`d_k=d_2[g(T(\min(k+1,K_2+1)),J)-g(\infty,J)]`$. Use the largest $`d_k`$ for upper and the smallest for lower. This is valid even for $`L=0`$, whose promotion term vanishes.

Multiplying the signed residual by $`(1+x)^2`$ yields the polynomial

```math
\sigma\{\rho N+(1+x)^2[q(U+C+dx)-\lambda w]\},
```

where

```math
N=b(1-q)(1+x)+c(1-q)[2q+(1+q)x]+e[(1-2q)x-qx^2],
\qquad U=-bx+c(-2qx+x^2)-ex^2.
```

These are $`(1+x)^2\Delta^+w/x`$ and $`\Delta^bw`$. On the larger full rectangle $`q\in[0,1],x\in[0,1/D_F]`$, the tensor Bernstein coefficients bound the polynomial by their maximum (the basis is nonnegative and sums to one). Use 16 equal $`q`$-intervals and exact rational power-to-Bernstein conversion. Divide a negative maximum by $`(1+1/D_F)^2`$, a nonnegative one by 1. This gives a finite universal directed bound even when the search's lower-tail Bernstein heuristic or its fallback did not establish a sign. It also proves the founder condition $`D|\Delta^+f|=O(1)`$ by the formula for $`N`$. All table and tail values and promotion jumps are bounded. At a tail endpoint of a motif bounded in absolute value by $`B`$, the leaf-gain jump is at most $`2B d_j/[D(D+1)]`$, by the difference of the old and new interpolation masses. At an exact internal position or exact endpoint the degree is bounded by the relevant cutoff. Thus every affected motif has bounded $`L|\Delta^+h|`$. The number of embeddings touching a fixed founder neighborhood is bounded by a polynomial in the core degrees one further fixed radius out. These are the local growth conditions used by the founder lemma.

## Multiplicities, comparison and the two physical credits

There is one singleton and one incoming edge per core in the ideal model. The highest signed cap for each contributes with its **sign**. Each bottom vertex receives at most one $`V`$, one $`P`$, and at most $`k_a-1`$ forks for its unique grandparent $`a`$. Including zero accounts for optional embeddings. At each real apex with $`k`$ children, there are $`\binom{k}{2}`$ cherries and stemmed cherries and $`\binom{k}{3}`$ three-child stars whenever supported. The high checker groups the positive branch deficits by apex state, divides their total by $`k`$ rounding up, then charges to the children. The sum of child counts is $`C(T)-1`$. Its bound plus the signed singleton and edge caps gives $`eS`$ in the replay table.

The comparison bounds assign the values $`f,C,Y,T`$ to their apex and $`g,V,P,F`$ to their bottom. For each possible bottom label take extrema over the other roles and the allowed multiplicities. Interpolate these bounds in the actual bottom state. For supported apex motifs bound each tuple by the extrema of that specific apex's table times its actual binomial count. The finite-degree comparison scans all $`D<D_F`$. In the analytic tail every bottom bound is affine in $`1/D`$ at fixed class, so the two endpoints $`0,1/D_F`$ bound it. For lower $`f\ge q(B_1-B_2-C_D/D_F)\ge0`$. For fixed $`k`$ use $`q\ge\max(0,1-(k+1)/D_F)`$. Beyond all clamped classes use $`q\ge0`$. The shift $`+\alpha(k-1)`$, $`\alpha=1/4`$ in the lower witness, gives a positive lower comparison in the zero-leaf tail, since $`\sum(k-1)=-1`$. The unshifted upper comparison and the lower bound $`m_0 C\le\Phi`$ follow. The separate upper row profile records the integer upper bound of this same allocation.

For any signed gain jump $`d`$ in ticks and any physical donor reciprocal $`i\in[0,1/2]`$, set

```math
J(d,i)=\tfrac12\max(d,0)-i d\ge0 .
```

**Singleton credit.** In a table-tail singleton the extremized contribution minus the true gain equals $`\sum_{y\sim u}J(\sigma\Delta f_u,i_y)`$. Each incidence is directed by its receiving vertex. A physical calculation can retain both orientations of an edge when the relevant receiving vertices qualify. It may restrict to $`d_2\le D<D_F-1`$ to use two table entries.

**Quartet credit.** At an inner contact $`x`$ of a supported quartet all external donors are uncentered. For $`m`$ such donors, the corner bound contributes $`m E[\max(0,d_{\rm corner})]/2`$, which is at least $`m\max(0,E d_{\rm corner})/2`$. Its excess over the actual external gain is at least $`\sum_y J(\sigma\Delta_x h_4,i_y)`$. When all four positions are inner, the conditional corner interpolation is a point mass. Selecting a fifth vertex uniquely indexes a quartet, its external donor and the contact to which the donor is adjacent. Different deletions of the same five-vertex set give different such incidences. All other relaxations in that row have nonnegative slack. The singleton and quartet credits are paid by separate true-gain terms in different rows. Thus **both may be added before applying the signed additive cap** even when individual rows have positive residual.

```math
\sigma S(\mathcal L\Phi-\lambda\Phi)+\sum J_f+\sum J_4\le (Se)C+E_{\rm founder,ticks}.
```

Selecting any subset of these nonnegative credits or rounding it down preserves the inequality. Here $`Se`$ is the reported integer per-core cap, and all $`J`$ values are in ticks. The physical occupation estimate controls the founder term and converts the discounted credits to exponent bounds.
