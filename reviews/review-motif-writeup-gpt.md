# Cold review of the motif product bounds and their checker

VERDICT: CORRECT. The two specified certificates prove the finite-size product bounds, and the revised proof and code map agree at the interpolation boundary.

## From the tree law to finite inequalities

1. **Tree law.** Root at the centre of the three-vertex tree. Summing the inverse source degrees over the neighbours of each potential receiver gives rate $`L+\rho`$ at a core receiver and total rate $`L/D`$ for its leaf children, after multiplication by $`n`$. Summing these rates counts each source once. A bud creates one new core with state $`(1,0,1)`$ and leaves its old parent's degree fixed.
2. **Drift identity.** A leaf gain changes only terms containing the receiver. A bud also creates the newborn vertex term, its incoming edge, at most one incoming vertical path and one cherry with each older core sibling. The two domain-entry cases have respectively zero earlier children or one earlier child, so entry activates no old path or old cherry. Splitting the old terms reproduces the three received knot contributions with multiplicities $`k_w`$, $`1`$ and $`k_p-1`$. A vertical path receives two centred residuals, and an unordered cherry receives one from each of its two incident edges. These assignments recover the event-wise identity upon summation. See [the motif proof, §§2–4](../proofs/bracket-motifs.md).
3. **Tail interpolation.** Write $`t=D_0/D`$. A leaf gain changes $`t`$ by $`-t/(D+1)`$, while a bud changes only the clamped child count. This gives the two separate bounded-rate terms at the finite tail endpoint, with zero contribution at infinity. A credit times a jump depending on the same endpoint is quadratic. Bounding the jump by its endpoint extremum and multiplying by the nonnegative credit gives the checked one-sided inequality. With two tail endpoints, first hold the other endpoint fixed and then interpolate its values. The knot extrema follow from linearity on each knot interval. Convexity for the maximum and concavity for the minimum give the required directions after interpolation.
4. **Finite families.** The vertex rows distinguish exact leaf saturation from leaf overflow, and treat child overflow for both root flags. In each overflow rectangle, taking all four corners covers the actual pair of rates. Tail–tail edges have zero table values and remain tail–tail. The middle of a path has $`k-1`$ unexposed core neighbours and a cherry apex has $`h+k-2`$, giving the two different half-degree allowances. Domain exits zero the affected motif. Equal child *states* at distinct vertices still contribute both endpoint changes. The triangular cherry loop includes its diagonal. See [the state layout](../verify/motif-bracket/cert_format.h) and [the exact checker](../verify/motif-bracket/cert_verify_gmp.c).
5. **Comparison and iteration.** Charging each edge and vertical path to its bottom and each cherry to its apex gives one vertex bound per core. Zero covers absent terms. The checker takes both candidate endpoints in each tail-state comparison, including the zero-weight infinity endpoint at $`D=D_0`$. This bounds each interpolated term in the appropriate direction. The comparison constants and the positive values at the three-vertex start give the two product prefactors below. Iteration uses positive factors $`1+\lambda/n`$ for both directions. The dyadic upper-tail probabilities at any strictly larger exponent are summable, and monotonicity fills in between dyadic sizes. No existence assumption enters either product bound.

## Code map

I traced the four row families and the comparison step to the loops identified in the [code map](../proofs/bracket-motifs-code-map.md). Enumerating the domains independently gives 572 inner states, 555 possible parent states, 290 child states, 84 path middles and 75 cherry apexes. The admissible edge pairs exclude 27 tail parent candidates times 15 tail child candidates. Thus the edge count is $`555\cdot290-27\cdot15=160545`$. The path and cherry counts are $`555\cdot84\cdot290=13519800`$ and $`75\cdot290\cdot291/2=3164625`$. The two root flags, the leaf-saturation cases and the two overflow states give 484138 vertex checks. All four counts match the executable output.

The comparisons with the exponent use rational or integer GMP operations after bounded native-integer arithmetic. For the two specified headers, the middle and apex domains lie inside the interpolation cutoff and the child-overflow projection lies entirely in the tail. The fixed-width factors in the direct tree evaluator use actual degrees. Its largest constructed starting tree has 803 vertices, so one additional birth gives the bound $`804^2`$ for a product of endpoint degrees. The finite inequalities use the separate common denominators in the path and cherry loops.

## Exact replay

I compiled `cert_build.c` and `cert_verify_gmp.c` with `cc -O3 -fopenmp`, linking the builder to `-lm` and the checker to `-lgmp`. I ran `cert_build 1 1133 2000 upper_5665.bin` and `cert_build 2 113 200 lower_5650.bin`, checked each SHA-256 digest, and ran `cert_verify_gmp` on each file with the tree checks enabled.

| File | SHA-256 | Exact prefactor |
|---|---|---|
| `upper_5665.bin` | `3a4aef888a5b567a8ada05f9546f98c979d361da9f9b33ba6e794eb91fb30abe` | $`14606753031/2989767809`$ |
| `lower_5650.bin` | `9e966fba9aee6eae3da659b693d3bb9c481c1cea66635b5218d0887f182a60d9` | $`7304352308/8763936965`$ |

Each run returned the following totals and exited successfully.

```text
  [6/6] Exact brute-force trees: 37 tested, fails: 0 (PASS)
=== SUMMARY: 17329108 exact GMP inequality checks, 0 failures => VERIFIED OK ===
```

## Independent recheck of the released tree test

I ran `python3 verify/tree-drift/motif_identity.py` with the upper and lower files identified by the hashes above. For Identity (I), `brute_drift` sums the inverse source degrees at each receiver and recomputes the full potential after attaching a leaf there. `decomposition` separately sums the four local expressions, including the newborn edge, path and cherries. For (5.0), each expression is compared with its own weighted term in the direction of the certificate. The term values also sum to the original potential. [Released test](../verify/tree-drift/motif_identity.py)

The command exited 0. Seeds 2–5 each ran 43 tree cases per certificate, including the same three fixed trees in each seed. Both certificates together gave 4,460 local checks on those cases and 1,838 on the five boundary trees. Every group reported `identity_fail=0 local_fail=0 compare_fail=0`. An independent count from the parent arrays gave 266 sample edges with only the parent in the tail, 64 with only the child in the tail, and 46 with both endpoints in the tail, matching the revised inventory. [Proof, §8](../proofs/bracket-motifs.md)
