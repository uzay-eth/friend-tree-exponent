# Astra cold attack on Neon-2 Claim 1

**VERDICT: HOLDS.** Scope: `ft2-neon2@47a79911dcc46ee9d54551aebb97a00a2b83b852`, `round2/neon2/limit-density/{README,REPORT}.md`: absolute continuity, exact conditional Beta split, and an everywhere-positive extended lower semicontinuous density with compact-local positive lower bounds. The multi-cut smoothing is valid in its stated localized sense; full-density continuity, upper bounds, and smoothness are not established. Dependencies are taken with their stated scope; no new certification of Claim 6 or the planted mean/gradient is implied.

## Attacks on the clock and selection

- **Two different physical directions.** A send from `u` to `v` attaches at `v`, but the kill must be charged to the **sending** side, at `1/D_u`. The reverse kill rate is `1/D_v`, which need not agree. With the slot present in the degree, one side's genuine rate is `m−1/D_u`; adding its own kill gives `m`. The independent local `Exp(m)` waits are those of the *killed construction*, conditional on survival in the local arrow marks. They are not the `Exp(m−1/D_u)` waits of the unmasked inert chain. This also reconciles Neon-2's factorization with Astra's path-dependent no-cross weight on the inert product law. Both are valid descriptions of the same surviving law. There is no conditioning on the edge being a hub–hub edge in the fixed-cut calculation: hub–hub existence only supplies a countable cover.
- **Infinite survival and split.** At a fixed labelled cut of fixed start size, survive-all is an intersection of cylinder tests of local arrows only. Product waits remain independent of all marks conditional on this intersection, by dominated convergence and a monotone class. Conditioning on both sides surviving factors across the two sides at this point. Their size masses therefore are independent Gamma variables, independent of the two embedded shape limits, and the Beta formula follows. For the Claim-6-only route, fixing marks and shifting one first wait gives two convergent discounted totals for almost every choice of the other data. Their difference is the stated contractive difference equation for the eventually bounded discounted side count. Doubling a surviving seed then rules out zero side limit by positive probability and independence. These are conditional statements on the surviving law.
- **Disjoint partition.** At the current start size all predecessor start states and cut edges are already known from the past marked history. A still-eligible predecessor fails later iff one of its edges is crossed; an edge of the current cut cannot cross, and any other such edge is inside a single component. Its old endpoint labels are fixed, so detecting its crossing from local arrow marks does not require interleaving, even though new global birth labels would. This verifies that the partition preserves **joint clock/mark independence**, not mutual independence of the coefficient limits after selection. Ordering by start size is material. A selection by the first global post-cut birth is a counterexample: the mean of the Beta split changes from `1/3` to `3/7` for initial sizes `(2,4)` if the first birth is required to be on side A. That selector is absent from the proof.

## Is the lower regularity genuinely stronger?

**Yes.** Astra's earlier a.e. positivity gives equivalence to Lebesgue measure but does not give a uniform lower bound for the mass per unit length of shrinking intervals near an arbitrary interior point. Here the surviving `P₄` cut alone dominates `c_x dy` on a whole neighborhood of every `x>0`. The support input is obtained under survival by a finite inert script followed by a union bound on the unconditional inert amplitude failure and the future kill failure, not by assuming that full support persists under arbitrary conditioning. With the two surviving coefficients near `a=b` strictly between `2^(μ−1)x` and `x`, a branch stays strictly interior and regular; its positive Beta density divided by the bounded absolute derivative gives the uniform lower bound. Finite subcovers give a compact-local bound.

For the l.s.c. version, each inverse branch is positive and continuous on an open image interval. Extension by zero is l.s.c. including at the omitted fold and output endpoints; taking zero at a fold is permitted since the critical preimage is a singleton. Integrating and summing by Fatou is valid for the countable disjoint sublaws; overlap would not justify a sum of densities. The full surviving `P₄` contribution occurs as a partition piece for each eligible label state, since the single core edge is unique there and distinct start states are mutually exclusive. Its interior lower bound remains pointwise in this version. The density may be infinite on a null set. Localized multi-cut Fourier estimates do not bound the sum of the pieces' derivative norms.

## Independent evidence

`ft2-astra/round2/attacks/neon2-density/check.py` imports no Neon code. From the Astra root:

```text
$ taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/attacks/neon2-density/check.py
PASS: asymmetric hazards 1/2 vs 1/3; 2304 path pairs; mark-only predecessor selection preserves both clock Laplace factors
COUNTEREXAMPLE: first-global-birth selection changes Beta mean 1/3 to 3/7
PASS: regular interior Beta(2,2) inverse branch at output 7/5
```

The independent rational fixtures test finite path likelihood, including an older cut edge becoming internal to the current right component, and do not stand in for the infinite-horizon argument. No Neon-2 files were edited.
