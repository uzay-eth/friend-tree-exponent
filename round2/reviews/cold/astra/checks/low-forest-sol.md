# Independent check: the fixed-threshold low forest

**Verdict.** The fixed-$`N`$ radius bound and component limits follow from Claim 10 and its physical spatial semigroup. I found no counterexample. The four-vertex calculation also applies to a physical cut component under the stated young-root and high-parent conditioning. [Addario-Berry et al., §7(5)](https://arxiv.org/html/2403.20185v1#S7) asks about the induced forest and resemblance to young friend trees, without specifying either these sampling laws or an age convention.

## Marked local limit and sample weights

Read original core degree as core valence plus leaf mark. A component with $`k`$ cores is detected inside radius $`k+1`$: walk at most $`k-1`$ edges to a retained core, one more to a boundary core and another to read that boundary core's full star. Adjacent original leaves of deleted parents become singletons; a retained core with no retained neighbours may also be a singleton. Deleted cores never count. Finiteness under the limit and the cemetery state turn pointwise cylinder convergence into total variation on the countable alphabet.

Include the core count in the type until after weighting. A component with $`k`$ cores contributes $`k`$ visits to the uniform-core numerator; applying $`1/k`$ counts it once per $`C_n`$. Transporting this weight to the unique top gives the same intensity, with the retained founder counted in finite trees. The short virgin event with a parent already above threshold supplies positive intensity even on nonisolated core-containing components. The ordinary-vertex and all-component laws concentrate on singleton original leaves since $`C_n=o(n)`$ and at most $`(N+1)C_n`$ vertices lie in core-containing components. Including these singleton leaves in the per-$`C_n`$ count makes the singleton intensity diverge.

## All-finite-threshold radius proof

In the backward window under the spatial eigenmeasure, sum over those cores in the focal old cell, including the old focal core itself. Conditional on the entire old environment, the focal cell core count is dominated by a rate-one Yule process with second moment $`2e^{2t}-e^t`$, since each core's promotion rate is at most one. The same rate bound integrated by depth gives at most $`t^j/j!`$ expected individuals at depth $`j`$ in any given cell. Cauchy–Schwarz against the focal cell count handles the potentially correlated selection of the new root. The factor $`e^{-\mu t}`$ is applied once to the full charge.

For any fixed old path of length $`r`$, take the marks at most $`1/N`$ on each directed rate-one edge driver. While all path vertices stay low these marks are accepted, and each adds a degree unit at an endpoint. On survival the number of base proposals by time $`t`$ cannot exceed $`(r+1)(N-2)`$. This compares an unconditional Poisson proposal count to a deterministic budget; it does not condition the drivers on survival. For $`t\ge4N^2`$ the budget is below half the mean $`2rt/N`$, giving $`e^{-rt/(4N)}`$.

A simple path stays in the focal cell or crosses the old tree from the focal anchor to the terminal anchor. An already-high anchor blocks entering the old tree but need not block a path turning below it within the focal cell. A length-$`3r`$ path has at most two recent segments and one old segment, so some segment has length at least $`r`$. If the old segment is shorter, every candidate terminal anchor other than the focal one belongs to the initially low ball of radius $`r-1`$. There are at most $`2N^r`$ candidate cells including the focal one, whether its old anchor is low or high. The long old charge uses at most $`N^r`$ old prefixes. No independence between cells or charges is required.

For $`k=3r\ge64N^3`$ and $`t=k/(16N)`$, the recent bound is at most $`3\exp(.435t)(3/4)^r`$; the log rate per $`k`$ is at most $`.435/32+\log(3/4)/3<-.08`$. The old term is at most $`2\exp(.435t-rt/(16N))`$ since $`t\ge8N\log N`$; also $`r/(16N)\ge4N^2/3\ge.435+1.28N`$ for every $`N\ge2`$. This gives the advertised exponential bound on the limiting radius. Rounding down to multiples of three costs a constant; bounded degrees give the stated polynomial volume tail. This does not bound the largest finite-tree component. In the successive threshold limit the coupled root component exhausts the locally finite decorated tree; that argument does not give a size law for uniform components after renormalizing their vanishing per-core intensity.

## Early-age Palm test

Sample the physical newborn law with parent degree $`d>N`$ **at birth** and observe at age $`a`$ tending to zero. The parent is still cut at observation by monotonicity. For $`N\ge4`$, exactly two inside births leave all four branch vertices eligible: the greatest possible original degree is four. Fewer births cannot make four, while more have probability $`O(a^3)`$, uniformly in the exterior by the birth-rate bound $`m+1/2`$ at size $`m`$. Hence pruning does not alter the leading coefficients. For each fixed environment the probability of a parent degree change tends to zero, including after a prescribed first inside birth. Bounded inside hazards permit dominated convergence in the two ordered hazard integrals and in the newborn environment mixture. Conditioning stable age to be at most $`\epsilon`$ gives the same limiting ratio because both coefficients carry the same age factor.

Collect target hazards as sums of reciprocal physical neighbour degrees. With $`u=1+1/d`$, the first birth lands at the top or its leaf at rates $`u,1/2`$. After the first target is the top, the second-birth star/path rates are $`u+1,2/3`$; after it is the leaf they are $`3/2, u`$. Thus the ordered weights are $`u(u+1)+3/4`$ and $`7u/6`$. The star weight minus twice the path weight equals $`(u-2/3)^2+11/36>0`$, so the ratio beats the ordinary RFT's $`2/3`$ star probability even after mixing over high-parent birth environments. As $`d\to\infty`$, the exterior-to-root hazard vanishes, but the parent slot continues to enter the focal degree, yielding $`33/47`$. This comparison does not identify the unconditioned size-four law under core or component sampling, or rule out a physically planted and age-mixed comparison law.

## Independent replay

`timeout 8s taskset -c 13 nice -n 10 python3 -B round2/checks/low_forest_sol.py`

```
cut: 5 components; inverse-core weighting and high-anchor turn pass
Palm: {5: Fraction(339, 479), 8: Fraction(67, 95), 20: Fraction(1161, 1651), 64: Fraction(34371, 48931)}; ordinary size-four star=2/3
```

The script uses ordinary adjacency lists and exact target rates on a complete physical fixture, enumerating interior and exterior targets. The tail and small-age arguments above are analytic.
