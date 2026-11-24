# D20 older-arm `hV` excess — independent cold audit

**Conditional HOLDS.** With the authenticated D20 signed-row certificate and the audited D20 occupation baseline as inputs, the older-arm excess may be added to the baseline numerator on each witness. The resulting outward interval is **[0.565681798679387, 0.565766116121802]**. This combines only the baseline and the size-disjoint older-arm payment. The alternative size-17/18 widening was not audited here and is not included. Composition with 8790's corrector is open: it would need a proof of separate charge ownership or a joint residual inequality, with a valid ceiling for the resulting potential.

## Physical incidence and sign before pricing

The paid row is `hV(A,B,C)` with actual oriented edges `A→B→C`. The extra tag `U` is another actual child of `A`. Read the **child receive port at slot A** of this three-role row. This port has multiplicity `k_A−1`, not `k_A`, and has no connection to the parent port or to the receive at B. Let `c` be its already signed slope in physical units (`sign/S` times the potential tick derivative after outgoing transfers). On this cut the transfer to `Yb` makes that slope zero whenever its two restricted slots `A,B` qualify; otherwise it is the signed exact gain difference of the `hV` potential at A. This identity was checked on the charged states against the authenticated D20 tables. In particular, it is not the unsigned difference, the `Yb` row, or a measure of the outgoing arrow from A.

The original corner residual is `a+Σ_i n_i min(0,c_i/2)≥0`; subtracting the corner contribution of the **particular donor** gives `c/D_U−min(0,c/2)`. This is `c/D_U` for `c≥0` and `|c|(1/2−1/D_U)` for `c<0`. Thus a virgin donor gives no negative-slope excess, but can give positive-slope excess. Replacing its degree by the recipient's, or resetting to a virgin, is not valid. On a given row we can sum these increments over any subset of distinct outside children, since each appears once in the original affine sum and all other increments and the corner residual are nonnegative. No row floor is paid again per tag. The exact replay rejects recipient degree, virgin reset, flipped sign, and including positive rows at the size-16 seam on **both** witnesses.

## Birth orders, true size and founder

Fix any labelled physical quadruple. `B` must precede `C`; otherwise only the order of `U` is free. `BCU`, `BUC`, and `UBC` exhaust the possibilities and partition completed incidences. At each matching promotion the old registration follows the actual post-promotion census and additionally emits an extended registration with exactly that new child; the two value terms use **one** physical promotion rate in the generator and **one** contribution to the holding rate. In `BCU`, the spine is already live when the side appears; in `BUC`, the side is already live when B first or later emits the selected C; in `UBC`, the existing side survives both later births. At each match distinct previously existing children retain their identities and current states. Linear summation of registrations does not assert independent histories. The checker also enumerates admissible chronologies of a labelled small tree with several eligible siblings and verifies multiplicity one for every completed incidence.

Known tagged edge sends use actual instantaneous degrees in **both** directions. Other neighbours remain present: a tracked vertex has `k+1−#tagged neighbours` unknown receiving incidences, each at most `1/2`. This covers the actual exterior even conditional on previous births, cap survival and tag history. Degree and arity cuts are irreversible for each registered tag. An exiting gain or promotion still pays its holding-rate cost with zero continuation. For any fixed history the discounted comparison uses the physical generator of the sum of prefixes and the lower bound at every allowed unknown-rate corner. One may stop at finite times and event cutoffs. Expectations of tagged registration counts at these finite depths are bounded by polynomials in elapsed time by rate-one promotion proposals per tracked parent. Discount by positive `HI` controls the remaining boundary term; optional dropping of untracked branches is only in the payoff, not a claim that they disappeared from the physical tree.

Once all four tags exist, `4+Σ L_tag` counts distinct ordinary vertices in the true downward fringe (none of its current leaves is itself a core tag). Requiring it to be **greater than 16** implies the true entire planted fringe exceeds 16 even with other branches present. Smaller registered states still propagate, because stock can later grow and the mask is on present payment, not a birth restriction. The baseline finite `g` and higher price stops at size 16 and four cores; for this comparison the size seam alone is sufficient. The baseline singleton, including its assigned newborn-edge potential term and its all-stock tail, is a different originally charged family. The actual founder has `h=0` and is **never** initialized as a virgin with a parent. All new prices omit it; its actual signed rows are nonnegative by the original certificate. No vanishing-founder asymptotic for this reward is needed. The founder remains in the whole-tree potential comparison and in the baseline accounting.

## Exact reprice and composition

At each state the independent checker solves the cross-multiplied scalar Bellman residual in integer ticks of `10^-24`, using breakpoints of the gain continuations rather than the production greedy quotient. The exact residual is nonnegative at the selected tick and negative at the next. The same inequality is explicitly checked at every independent unknown-rate corner on the size-17 states. The discount is the physical upper trial `HI=.565768`; multiply the discounted newborn value by `LO=.565680` once to get a lower stationary charge per core. The birth-compensation argument from the baseline applies to the apex-born aggregate: histories are counted from actual nonfounder births, and the number of simultaneously paid incidences per apex is bounded by the arity mask. Rewards on the admitted censuses are bounded. Discounted comparison covers every live mother environment. To see the same argument via finite-tree Laplace normalization, the normalized number of births tends to `μ`. Taking a discount down to `μ` in the age integrals of these bounded rewards gives the same uniform lower payment; `μ≤HI` and the reward is nonnegative. There is no Perron mark mixture or empirical birth-frequency normalization.

| Ordered source | Lower newborn ticks | Upper newborn ticks |
|---|---:|---:|
| `BCU` | 71317489466711121 | 309675901926955038 |
| `BUC` | 48898507302423287 | 104371201582794230 |
| `UBC` | 75959023576674456 | 160881545356374804 |

Exact older per-core rewards are `43348549027037952417/390625000000000000000000000` and `508165059516545414139/1562500000000000000000000000`, respectively. These imply the reported downward-decimal rewards `1.10972285509e−7 / 3.25225638090e−7`. The exact baseline inputs are the sums of the first-body/all-stock and joint `g`/higher ticks from the [audited D20 price](neon20-occupation-sol.md); no rounded baseline endpoint is used in the composition:

| | Lower | Upper |
|---|---|---|
| Baseline numerator | `2911747728504795030087/1250000000000000000000000000` | `2808168533038621918443/1250000000000000000000000000` |
| Numerator with older excess | `15252315426956582389107/6250000000000000000000000000` | `16073502903259291248771/6250000000000000000000000000` |
| Original D20 stationary ceiling | `16959459852924085304227707/12500000000000000000000000` | `17064269776474933055949939/12500000000000000000000000` |

The numerator is still the signed drift of the **same** D20 potential at each original trial, lower-bounded by two disjoint nonnegative payments. Its positive potential and its D20 stationary ceiling are unchanged. Add the lower quotient to `.565680`, subtract the upper quotient from `.565768`, and round outward to 15 decimals.

## Replay and freeze

The checker pins the older README and price code, the baseline replay, both exports and the original charged-row extractor by SHA-256. It uses the same two unchanged local D20 payloads independently hashed in the baseline audit; no additional full binary or GMP sweep. The already-audited row evaluator is shared; the signed slot identity, chronological model, exact rounding and ledger are checked independently. Run from `[local path]`:

```sh
taskset -c 8 nice -n 17 round2/.venv/bin/python -Bu round2/checks/neon20-older-sol.py lower
taskset -c 8 nice -n 17 round2/.venv/bin/python -Bu round2/checks/neon20-older-sol.py upper
```

Both exited zero. Per sign: 153,121 / 153,303 / 153,303 cached states including exits in the three orders; 15,744 extension transitions; 10,554 explicit-corner states; 10,542 size-16 seam checks. Positive charged states cover both slope signs; the wrong-recipient test differed on 15,196 lower and 16,634 upper states, and the false size-16 credit on 1,674 lower and 1,872 upper visits. Final outward integer ticks were `565681798679387` and `565766116121802`. Only the two cold-audit files were added.
