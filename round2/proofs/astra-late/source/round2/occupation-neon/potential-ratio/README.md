> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Pointwise comparison of the two original Neon D18 potentials

For every finite physical tree with a core, with the witnesses rooted at the same founder,

```text
999/1000 ≤ Φ_hi(T)/Φ_lo(T) ≤ 28/25.
```

`lo` is the `565673/1000000` witness (mode 2) and `hi` is the `22631/40000` witness (mode 1). Their original positive comparison certificates ensure that the quotient is defined. This proof does not divide separate absolute comparison constants. The endpoints need not be optimal. A physical crown with 18 short core arms already gives a ratio of about `1.064316611055262`; an upper endpoint within one percent of 1 is impossible globally for these payloads.

## Exact comparison before any bounding

At a rational trial `q=N/B`, form `B*hi-N*lo` entry by entry for every potential table on the common layout. Both scales are the same exact integer `10^11`; there is no rounding of either coefficient and no magnitude division in the checker. The normalized score of the resulting tables is exactly `B Φ_hi-N Φ_lo`: physical Bernstein embeddings act linearly with the same nonnegative weights on both payloads. Take an extremum of the **difference** controls, or of compatible combinations thereof. On a tail use the fixed core-child count and parent flag, all four Bernstein controls (with zero only for masked controls), and the full first-body range over legal leaf counts. At and above the stored leaf cap the first-body entry is constant. The child-count overflow bucket starts at 181, with true child count potentially unbounded.

## Lower certificate

Charge edge, P3, P4 and bent motifs to their bottom, and the sibling pair, Y-tail and triple to the branching vertex. A bottom has at most one vertical path of each length. For a bent motif the bottom has unique parent and grandparent; the apex can designate at most `k_apex-1` other children, so maximize or minimize only after multiplying each apex slice by its exact multiplicity. Optional bottom motifs allow zero; the nonroot incoming edge is mandatory. For the branching motifs the exact multiplicities are binomial coefficients of the actual number of children in the active range. The two tail/inner interpretations of sibling controls use the same physical degree and legal masks as the witness evaluator. The loop scans all legal leaf counts through their cap for each `(k,h)` and then all legal overflow types. For the lower comparison, a zero-sum term of `50000000000*(k-h)` is added to each unscaled integer row. Since `Σ k=Σ h`, it adds nothing to the tree score. Its slope is positive on either overflow type, so testing at 181 suffices there. The smallest recharged row is strictly positive (`22944552835` raw units).

## Upper certificate: recharge along actual edges

The same grouping and scalar recharge gave `5/4`, but fails for `28/25`. To retain some parent/child compatibility, remove the edge and P3 from the bottom bounds. Charge each P3 instead to its middle vertex `x`, grouped on its incoming edge from `p`. For each physical type pair `(p,x)`, bound the edge by the largest compatible difference control and P3 by `k_x` times the largest control compatible with both `p,x`, after maximizing out the bottom. If P3 is disabled or masked at a compared middle control, include zero for that choice. Call this combined bound `G(p,x)` and the rest of the local vertex cap `A(x)`. In the overflow bucket P3 is disabled, so `G` and `A` require no growth factor in the true `k` there. Inner types retain exact `(L,k,h)`; tail types retain `(k,h)` and a common all-leaves cap.

Send a charge `m(p)` to `p` for each child `x` and subtract the same charge from `x`. This zero-sum recharge gives the local caps

```text
founder x:                  A(x) + k_x*m(x)
nonfounder x with parent p: A(x) + k_x*m(x) + G(p,x) - m(p).
```

Generate a candidate by integer max-plus updates `m(p) ← max(m(p), max_x [A(x)+k_x*m(x)+G(p,x)])`. Two updates suffice from the given start. The proof is the final exact verification of all edge inequalities, all founder rows, and the nonpositive overflow slopes, regardless of how the candidate was generated. The two overflow slopes are both `-241991156468`, so the test at `k=181` bounds every larger `k`. The smallest founder reserve is `438056321030`; some edge constraints are tight. In both comparisons the checker range-checks signed products and sums before performing native integer arithmetic; rational interval arithmetic for the scalar recharge uses Python integers. A failed trial is not a counterexample to the ratio.

## Replay

Run from the project root with Python and NumPy. The default run verifies both full original SHA-256 pins, checks identity of the role/embedding layouts, and verifies both endpoints. It also enumerates 2,956 chronological trees through size eight and explicit mixed tails, wide crowns, paths, and capped/overflow fixtures, using the separate exact whole-tree motif-product evaluator from `attacks/neon1/audit.py`. The reported empirical extrema are diagnostic only; assertions use exact rational scores. No witness files are copied.

```sh
nice -n 15 taskset -c 3 timeout 40s python3 -B round2/occupation-neon/potential-ratio/check.py b3f18_lo.bin b3f18_hi.bin
```

The replay passed on CPU 3: 2,984 whole trees, all eight motif families, observed ratios about `0.9999401226701936` to `1.064316611055262`. At the upper empirical extreme the founder state is `(0,18,0)` with 18 newborn core children. These are physical trees with a core, the domain of the original score certificates; the bare size-two edge is not assigned a quotient.
